//  ██████   ██████ ██████████  █████████  █████   █████ █████    ███████
// ░░██████ ██████ ░░███░░░░░█ ███░░░░░███░░███   ░░███ ░░███   ███░░░░░███      ███         ███
//  ░███░█████░███  ░███  █ ░ ░███    ░░░  ░███    ░███  ░███  ███     ░░███    ░███        ░███
//  ░███░░███ ░███  ░██████   ░░█████████  ░███████████  ░███ ░███      ░███ ███████████ ███████████
//  ░███ ░░░  ░███  ░███░░█    ░░░░░░░░███ ░███░░░░░███  ░███ ░███      ░███░░░░░███░░░ ░░░░░███░░░
//  ░███      ░███  ░███ ░   █ ███    ░███ ░███    ░███  ░███ ░░███     ███     ░███        ░███
//  █████     █████ ██████████░░█████████  █████   █████ █████ ░░░███████░      ░░░         ░░░
// ░░░░░     ░░░░░ ░░░░░░░░░░  ░░░░░░░░░  ░░░░░   ░░░░░ ░░░░░    ░░░░░░░
//
//
//  License:         MIT License
//                   SPDX-License-Identifier: MIT
//                   meshio++ default license: LICENSE
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//
// iga_evaluate, iga_refine, iga_sample. See operations/iga.hpp for the contract.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/iga.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/cad_json.hpp"
#include "meshioplusplus/formats/mfem.hpp"
#include "meshioplusplus/parallel.hpp"

// Project includes (private, not installed)
#include "../detail/brep_trim.hpp"
#include "../detail/iga_kernels.hpp"
#include "../formats/lagrange_common.hpp"

namespace meshioplusplus {
namespace {

constexpr const char* kIgaPrefix = "meshio++: iga: ";
constexpr int kIgaMaxSubdivisions = 64;

[[noreturn]] void iga_fail(const std::string& rWhat) {
    throw std::invalid_argument(kIgaPrefix + rWhat);
}

// What a patch's arrays say it is.
struct IgaShape {
    std::size_t mDim = 0;
    std::size_t mSdim = 0;
    std::array<std::size_t, 3> mCounts{1, 1, 1};
    std::size_t mNumPoints = 0;
    std::vector<std::size_t> mFieldComps;
};

IgaShape iga_check_patch(const NurbsPatch& rP, const std::string& rWhat) {
    IgaShape sh;
    sh.mDim = rP.mDegrees.size();
    if (sh.mDim < 1 || sh.mDim > 3 || rP.mKnots.size() != sh.mDim)
        iga_fail(rWhat + " needs 1 to 3 degrees and as many knot vectors");
    sh.mNumPoints = 1;
    for (std::size_t d = 0; d < sh.mDim; ++d) {
        const int p = rP.mDegrees[d];
        if (p < 1 || p > detail::kBrepMaxDegree)
            iga_fail(rWhat + " degree " + std::to_string(p) + " is outside 1.." +
                     std::to_string(detail::kBrepMaxDegree));
        const std::vector<double>& k = rP.mKnots[d];
        if (k.size() < 2 * static_cast<std::size_t>(p) + 2)
            iga_fail(rWhat + " knot vector is too short");
        for (std::size_t i = 0; i < k.size(); ++i)
            if (!std::isfinite(k[i]) || (i && k[i] < k[i - 1]))
                iga_fail(rWhat + " knot vector is not finite and non-decreasing");
        sh.mCounts[d] = k.size() - static_cast<std::size_t>(p) - 1;
        if (k[static_cast<std::size_t>(p)] >= k[sh.mCounts[d]])
            iga_fail(rWhat + " knot vector in direction " + std::to_string(d) +
                     " has zero active span length");
        sh.mNumPoints *= sh.mCounts[d];
    }
    const NDArray& cps = rP.mControlPoints;
    if (cps.Dtype() != DType::Float64 || cps.Shape().size() != 2 ||
        cps.Shape()[0] != sh.mNumPoints || cps.Shape()[1] < 1 || cps.Shape()[1] > 3)
        iga_fail(rWhat + " control points are not Float64 (" + std::to_string(sh.mNumPoints) +
                 ", 1..3) for its knot vectors");
    sh.mSdim = cps.Shape()[1];
    if (rP.mWeights.Dtype() != DType::Float64 || rP.mWeights.Size() != sh.mNumPoints)
        iga_fail(rWhat + " has " + std::to_string(rP.mWeights.Size()) + " weights for " +
                 std::to_string(sh.mNumPoints) + " control points");
    for (const NDArray& f : rP.mFields) {
        if (f.Dtype() != DType::Float64 || f.Shape().size() != 2 || f.Shape()[0] != sh.mNumPoints ||
            f.Shape()[1] < 1)
            iga_fail(rWhat + " has a field that is not Float64 (" + std::to_string(sh.mNumPoints) +
                     ", components)");
        sh.mFieldComps.push_back(f.Shape()[1]);
    }
    return sh;
}

bool iga_is_rational(const NurbsPatch& rP) {
    const double* w = rP.mWeights.As<double>();
    for (std::size_t i = 0; i < rP.mWeights.Size(); ++i)
        if (w[i] != 1.0)
            return true;
    return false;
}

void iga_check_weights(const NurbsPatch& rP, const std::string& rWhat) {
    const double* w = rP.mWeights.As<double>();
    for (std::size_t i = 0; i < rP.mWeights.Size(); ++i)
        if (!std::isfinite(w[i]) || w[i] == 0.0)
            iga_fail(rWhat + " has a zero or non-finite weight");
}

// The homogeneous net: (w x, w fields..., w), or (x, fields...) when polynomial.
detail::IgaNet iga_make_net(const NurbsPatch& rP, const IgaShape& rSh, bool Rational) {
    detail::IgaNet net;
    net.mDim = rSh.mDim;
    for (std::size_t d = 0; d < rSh.mDim; ++d) {
        net.mDegrees[d] = rP.mDegrees[d];
        net.mKnots[d] = rP.mKnots[d];
        net.mCounts[d] = rSh.mCounts[d];
    }
    std::size_t comps = rSh.mSdim;
    for (const std::size_t c : rSh.mFieldComps)
        comps += c;
    net.mColumns = comps + (Rational ? 1 : 0);
    net.mValues.resize(rSh.mNumPoints * net.mColumns);
    const double* w = rP.mWeights.As<double>();
    for (std::size_t i = 0; i < rSh.mNumPoints; ++i) {
        double* row = net.mValues.data() + i * net.mColumns;
        std::size_t c = 0;
        const double wi = Rational ? w[i] : 1.0;
        const double* x = rP.mControlPoints.As<double>() + i * rSh.mSdim;
        for (std::size_t d = 0; d < rSh.mSdim; ++d)
            row[c++] = Rational ? wi * x[d] : x[d];
        for (std::size_t f = 0; f < rSh.mFieldComps.size(); ++f) {
            const double* v = rP.mFields[f].As<double>() + i * rSh.mFieldComps[f];
            for (std::size_t d = 0; d < rSh.mFieldComps[f]; ++d)
                row[c++] = Rational ? wi * v[d] : v[d];
        }
        if (Rational)
            row[c] = wi;
    }
    return net;
}

NurbsPatch iga_patch_from_net(const detail::IgaNet& rNet, const IgaShape& rSh, bool Rational) {
    NurbsPatch out;
    const std::size_t n = rNet.NumPoints();
    for (std::size_t d = 0; d < rSh.mDim; ++d) {
        out.mDegrees.push_back(rNet.mDegrees[d]);
        out.mKnots.push_back(rNet.mKnots[d]);
    }
    out.mControlPoints = NDArray(DType::Float64, {n, rSh.mSdim});
    out.mWeights = NDArray(DType::Float64, {n});
    for (const std::size_t c : rSh.mFieldComps)
        out.mFields.emplace_back(DType::Float64, std::vector<std::size_t>{n, c});
    double* px = out.mControlPoints.As<double>();
    double* pw = out.mWeights.As<double>();
    for (std::size_t i = 0; i < n; ++i) {
        const double* row = rNet.mValues.data() + i * rNet.mColumns;
        const double w = Rational ? row[rNet.mColumns - 1] : 1.0;
        std::size_t c = 0;
        for (std::size_t d = 0; d < rSh.mSdim; ++d, ++c)
            px[i * rSh.mSdim + d] = Rational ? row[c] / w : row[c];
        for (std::size_t f = 0; f < rSh.mFieldComps.size(); ++f) {
            double* pv = out.mFields[f].As<double>() + i * rSh.mFieldComps[f];
            for (std::size_t d = 0; d < rSh.mFieldComps[f]; ++d, ++c)
                pv[d] = Rational ? row[c] / w : row[c];
        }
        pw[i] = w;
    }
    return out;
}

void iga_check_direction(const IgaShape& rSh, std::int32_t Direction) {
    if (Direction < 0 || static_cast<std::size_t>(Direction) >= rSh.mDim)
        iga_fail("direction " + std::to_string(Direction) + " is outside the patch's " +
                 std::to_string(rSh.mDim) + " parametric directions");
}

NurbsPatch iga_refine_patch(const NurbsPatch& rP, bool Rational, const IgaRefineOptions& rOpt,
                            const std::string& rWhat) {
    const IgaShape sh = iga_check_patch(rP, rWhat);
    if (Rational)
        iga_check_weights(rP, rWhat);
    // Nothing to do: the patch itself, not its homogeneous round trip (which
    // can move a rational control point by an ulp).
    bool identity = true;
    for (std::size_t d = 0; d < sh.mDim; ++d)
        identity = identity && rOpt.mElevate[d] == 0 && rOpt.mSubdivisions[d] == 1;
    if (identity && (Rational || !iga_is_rational(rP)))
        return rP;
    detail::IgaNet net = iga_make_net(rP, sh, Rational);
    for (std::size_t d = 0; d < sh.mDim; ++d) {
        if (rOpt.mElevate[d] < 0)
            iga_fail("cannot elevate a degree by " + std::to_string(rOpt.mElevate[d]));
        if (net.mDegrees[d] + rOpt.mElevate[d] > detail::kBrepMaxDegree)
            iga_fail(rWhat + ": the degree would pass " + std::to_string(detail::kBrepMaxDegree));
        detail::iga_net_elevate(net, d, rOpt.mElevate[d]);
    }
    for (std::size_t d = 0; d < sh.mDim; ++d) {
        if (rOpt.mSubdivisions[d] < 1 || rOpt.mSubdivisions[d] > kIgaMaxSubdivisions)
            iga_fail("subdivisions must be 1 to " + std::to_string(kIgaMaxSubdivisions));
        detail::iga_net_insert(
            net, d, detail::iga_split_knots(net.mDegrees[d], net.mKnots[d], rOpt.mSubdivisions[d]));
    }
    return iga_patch_from_net(net, sh, Rational);
}

// The values at lattice nodes: the geometry, then the fields, in one table.
struct IgaTable {
    detail::IgaEvalNet mNet;
    std::vector<double> mValues;
    IgaShape mShape;
    bool mRational = false;
};

IgaTable iga_make_table(const NurbsPatch& rP, const std::string& rWhat, bool Rational) {
    IgaTable t;
    t.mShape = iga_check_patch(rP, rWhat);
    t.mRational = Rational;
    if (Rational)
        iga_check_weights(rP, rWhat);
    const IgaShape& sh = t.mShape;
    std::size_t nv = sh.mSdim;
    for (const std::size_t c : sh.mFieldComps)
        nv += c;
    t.mValues.resize(sh.mNumPoints * nv);
    for (std::size_t i = 0; i < sh.mNumPoints; ++i) {
        double* row = t.mValues.data() + i * nv;
        std::size_t c = 0;
        const double* x = rP.mControlPoints.As<double>() + i * sh.mSdim;
        for (std::size_t d = 0; d < sh.mSdim; ++d)
            row[c++] = x[d];
        for (std::size_t f = 0; f < sh.mFieldComps.size(); ++f) {
            const double* v = rP.mFields[f].As<double>() + i * sh.mFieldComps[f];
            for (std::size_t d = 0; d < sh.mFieldComps[f]; ++d)
                row[c++] = v[d];
        }
    }
    t.mNet.mDim = sh.mDim;
    for (std::size_t d = 0; d < sh.mDim; ++d) {
        t.mNet.mDegrees[d] = rP.mDegrees[d];
        t.mNet.mKnots[d] = &rP.mKnots[d];
        t.mNet.mCounts[d] = sh.mCounts[d];
    }
    t.mNet.mValues = t.mValues.data();
    t.mNet.mColumns = nv;
    t.mNet.mWeights = Rational ? rP.mWeights.As<double>() : nullptr;
    return t;
}

IgaEvaluation iga_evaluate_patch(const NurbsPatch& rP, bool Rational, const NDArray& rParams,
                                 std::int32_t Derivatives, const std::string& rWhat) {
    if (Derivatives < 0 || Derivatives > 2)
        iga_fail("derivatives must be 0, 1 or 2");
    const IgaTable t = iga_make_table(rP, rWhat, Rational);
    const IgaShape& sh = t.mShape;
    if (rParams.Dtype() != DType::Float64 ||
        !(rParams.Shape().size() == 2 && rParams.Shape()[1] == sh.mDim) &&
            !(sh.mDim == 1 && rParams.Shape().size() == 1))
        iga_fail("parameters are not Float64 (N, " + std::to_string(sh.mDim) + ")");
    const std::size_t n = rParams.Shape()[0];
    const std::size_t dim = sh.mDim, sdim = sh.mSdim, nv = t.mNet.mColumns;
    std::array<double, 3> lo{0, 0, 0}, hi{0, 0, 0}, tol{0, 0, 0};
    for (std::size_t d = 0; d < dim; ++d) {
        lo[d] = rP.mKnots[d][static_cast<std::size_t>(rP.mDegrees[d])];
        hi[d] = rP.mKnots[d][sh.mCounts[d]];
        tol[d] = 1e-12 * (hi[d] - lo[d]);
    }
    const double* params = rParams.As<double>();
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t d = 0; d < dim; ++d) {
            const double u = params[i * dim + d];
            if (!(u >= lo[d] - tol[d] && u <= hi[d] + tol[d]))
                iga_fail("parameter " + std::to_string(i) + ", direction " + std::to_string(d) +
                         " is outside the knot range [" + std::to_string(lo[d]) + ", " +
                         std::to_string(hi[d]) + "]");
        }
    IgaEvaluation out;
    out.mPoints = NDArray(DType::Float64, {n, sdim});
    std::vector<double> values(n * nv), ders(Derivatives ? n * dim * nv : 0),
        ders2(Derivatives == 2 ? n * dim * dim * nv : 0);
    parallel_for(n, [&](std::size_t i) {
        std::array<std::size_t, 3> span{0, 0, 0};
        std::array<double, 3> u{0, 0, 0};
        for (std::size_t d = 0; d < dim; ++d) {
            u[d] = std::min(std::max(params[i * dim + d], lo[d]), hi[d]);
            span[d] = detail::nurbs_find_span(rP.mDegrees[d], rP.mKnots[d], sh.mCounts[d], u[d]);
        }
        if (Derivatives == 2)
            detail::iga_eval_point_d2(t.mNet, span.data(), u.data(), values.data() + i * nv,
                                      ders.data() + i * dim * nv,
                                      ders2.data() + i * dim * dim * nv);
        else
            detail::iga_eval_point(t.mNet, span.data(), u.data(), Derivatives,
                                   values.data() + i * nv,
                                   Derivatives ? ders.data() + i * dim * nv : nullptr);
    });
    double* px = out.mPoints.As<double>();
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t d = 0; d < sdim; ++d)
            px[i * sdim + d] = values[i * nv + d];
    std::size_t col = sdim;
    for (const std::size_t comps : sh.mFieldComps) {
        NDArray f(DType::Float64, {n, comps});
        double* pf = f.As<double>();
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t d = 0; d < comps; ++d)
                pf[i * comps + d] = values[i * nv + col + d];
        out.mFields.push_back(std::move(f));
        col += comps;
    }
    if (Derivatives) {
        out.mDerivatives = NDArray(DType::Float64, {n, dim, sdim});
        double* pd = out.mDerivatives.As<double>();
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t k = 0; k < dim; ++k)
                for (std::size_t d = 0; d < sdim; ++d)
                    pd[(i * dim + k) * sdim + d] = ders[(i * dim + k) * nv + d];
        if (Derivatives == 2) {
            out.mSecondDerivatives = NDArray(DType::Float64, {n, dim, dim, sdim});
            double* ps = out.mSecondDerivatives.As<double>();
            for (std::size_t i = 0; i < n; ++i)
                for (std::size_t k = 0; k < dim; ++k)
                    for (std::size_t l = 0; l < dim; ++l)
                        for (std::size_t d = 0; d < sdim; ++d)
                            ps[((i * dim + k) * dim + l) * sdim + d] =
                                ders2[((i * dim + k) * dim + l) * nv + d];
        }
        if (dim == 2 && sdim == 3) {
            out.mNormals = NDArray(DType::Float64, {n, std::size_t{3}});
            double* pn = out.mNormals.As<double>();
            for (std::size_t i = 0; i < n; ++i) {
                const double* a = pd + (i * 2) * 3;
                const double* b = pd + (i * 2 + 1) * 3;
                const double cx = a[1] * b[2] - a[2] * b[1];
                const double cy = a[2] * b[0] - a[0] * b[2];
                const double cz = a[0] * b[1] - a[1] * b[0];
                const double len = std::sqrt(cx * cx + cy * cy + cz * cz);
                if (len > 0.0) {
                    pn[i * 3] = cx / len;
                    pn[i * 3 + 1] = cy / len;
                    pn[i * 3 + 2] = cz / len;
                } else {
                    pn[i * 3] = pn[i * 3 + 1] = pn[i * 3 + 2] = 0.0;
                }
            }
        }
    }
    return out;
}

std::string iga_quad_family(std::size_t Dim, int Order) {
    if (Dim == 1)
        return Order <= 1 ? "line" : Order == 2 ? "line3" : "VTK_LAGRANGE_CURVE";
    if (Dim == 2)
        return Order <= 1 ? "quad" : Order == 2 ? "quad9" : "VTK_LAGRANGE_QUADRILATERAL";
    return Order <= 1 ? "hexahedron" : Order == 2 ? "hexahedron27" : "VTK_LAGRANGE_HEXAHEDRON";
}

}  // namespace

NurbsPatch iga_insert_knots(const NurbsPatch& rPatch, std::int32_t Direction,
                            const std::vector<double>& rKnots, bool Rational) {
    const IgaShape sh = iga_check_patch(rPatch, "patch");
    iga_check_direction(sh, Direction);
    if (Rational)
        iga_check_weights(rPatch, "patch");
    detail::IgaNet net = iga_make_net(rPatch, sh, Rational);
    detail::iga_net_insert(net, static_cast<std::size_t>(Direction), rKnots);
    return iga_patch_from_net(net, sh, Rational);
}

NurbsPatch iga_elevate_degree(const NurbsPatch& rPatch, std::int32_t Direction, std::int32_t Times,
                              bool Rational) {
    const IgaShape sh = iga_check_patch(rPatch, "patch");
    iga_check_direction(sh, Direction);
    if (Times < 0)
        iga_fail("cannot elevate a degree by " + std::to_string(Times));
    if (rPatch.mDegrees[static_cast<std::size_t>(Direction)] + Times > detail::kBrepMaxDegree)
        iga_fail("the degree would pass " + std::to_string(detail::kBrepMaxDegree));
    if (Rational)
        iga_check_weights(rPatch, "patch");
    detail::IgaNet net = iga_make_net(rPatch, sh, Rational);
    detail::iga_net_elevate(net, static_cast<std::size_t>(Direction), Times);
    return iga_patch_from_net(net, sh, Rational);
}

NurbsInfo iga_refine(const NurbsInfo& rInfo, const IgaRefineOptions& rOptions) {
    if (rInfo.mPatches.empty())
        iga_fail("the model has no patches");
    if (rInfo.mTopology.mDimension != 0) {
        // The patches' directions follow their own edges, so a per-direction
        // refinement would not agree on the knot vector two patches share.
        const auto dim = static_cast<std::size_t>(rInfo.mTopology.mDimension);
        for (std::size_t d = 1; d < std::min<std::size_t>(dim, 3); ++d)
            if (rOptions.mSubdivisions[d] != rOptions.mSubdivisions[0] ||
                rOptions.mElevate[d] != rOptions.mElevate[0])
                iga_fail(
                    "a model with a patch topology refines every direction alike (its patches "
                    "share knot vectors along differently oriented directions)");
    }
    NurbsInfo out;
    out.mFieldNames = rInfo.mFieldNames;
    out.mTopology = rInfo.mTopology;
    for (std::size_t i = 0; i < rInfo.mPatches.size(); ++i) {
        const NurbsPatch& p = rInfo.mPatches[i];
        const std::string what = "patch " + std::to_string(i);
        iga_check_patch(p, what);
        out.mPatches.push_back(iga_refine_patch(p, iga_is_rational(p), rOptions, what));
    }
    out.mCellPatch = NDArray(DType::Int64, {0});
    out.mCellSpan = NDArray(DType::Int64, {0, 3});
    out.mBoundaryCellPatch = NDArray(DType::Int64, {0});
    return out;
}

BrepInfo iga_refine(const BrepInfo& rInfo, const IgaRefineOptions& rOptions) {
    BrepInfo out = rInfo;
    std::vector<int> flag(rInfo.mPatches.size(), -1);  // -1 unreferenced, 0/1 face flag
    for (std::size_t f = 0; f < rInfo.mFaces.size(); ++f) {
        const BrepFace& face = rInfo.mFaces[f];
        if (face.mPatch < 0 || static_cast<std::size_t>(face.mPatch) >= rInfo.mPatches.size())
            iga_fail("face " + std::to_string(f) + " names no patch");
        int& fl = flag[static_cast<std::size_t>(face.mPatch)];
        const int want = face.mRational ? 1 : 0;
        if (fl >= 0 && fl != want)
            iga_fail("patch " + std::to_string(face.mPatch) +
                     " is shared by faces that disagree on rational");
        fl = want;
    }
    for (std::size_t i = 0; i < rInfo.mPatches.size(); ++i) {
        const NurbsPatch& p = rInfo.mPatches[i];
        const std::string what = "patch " + std::to_string(i);
        const IgaShape sh = iga_check_patch(p, what);
        if (sh.mDim != 2)
            iga_fail(what + " is not a surface patch");
        const bool rational = flag[i] < 0 ? iga_is_rational(p) : flag[i] == 1;
        out.mPatches[i] = iga_refine_patch(p, rational, rOptions, what);
    }
    for (BrepFace& face : out.mFaces)
        face.mControlPointIds = NDArray(DType::Int64, {0});
    return out;
}

IgaEvaluation iga_evaluate(const NurbsInfo& rInfo, std::size_t Patch, const NDArray& rParams,
                           std::int32_t Derivatives) {
    if (Patch >= rInfo.mPatches.size())
        iga_fail("patch " + std::to_string(Patch) + " is outside the model's " +
                 std::to_string(rInfo.mPatches.size()));
    const NurbsPatch& p = rInfo.mPatches[Patch];
    return iga_evaluate_patch(p, iga_is_rational(p), rParams, Derivatives,
                              "patch " + std::to_string(Patch));
}

IgaEvaluation iga_evaluate(const BrepInfo& rInfo, std::size_t Face, const NDArray& rParams,
                           std::int32_t Derivatives, std::int32_t TrimSamples) {
    if (Face >= rInfo.mFaces.size())
        iga_fail("face " + std::to_string(Face) + " is outside the model's " +
                 std::to_string(rInfo.mFaces.size()));
    if (TrimSamples < 2 || TrimSamples > 4096)
        iga_fail("trim_samples must be 2 to 4096");
    const BrepFace& face = rInfo.mFaces[Face];
    const std::string what = "face " + std::to_string(Face);
    if (face.mPatch < 0 || static_cast<std::size_t>(face.mPatch) >= rInfo.mPatches.size())
        iga_fail(what + " names no patch");
    const NurbsPatch& patch = rInfo.mPatches[static_cast<std::size_t>(face.mPatch)];
    if (patch.mDegrees.size() != 2)
        iga_fail(what + " is not a surface patch");
    for (const BrepLoop& loop : face.mLoops)
        for (const BrepTrim& t : loop.mTrims)
            detail::brep_check_curve(t.mCurve, what + " trimming curve", kIgaPrefix);
    IgaEvaluation out = iga_evaluate_patch(patch, face.mRational, rParams, Derivatives, what);
    const std::size_t n = rParams.Shape()[0];
    out.mInside = NDArray(DType::UInt8, {n});
    std::uint8_t* inside = out.mInside.As<std::uint8_t>();
    const bool clip = face.mTrimmed && !face.mLoops.empty();
    std::vector<std::vector<std::array<double, 2>>> loops;
    if (clip)
        for (const BrepLoop& loop : face.mLoops) {
            auto poly = detail::brep_loop_polygon(loop, TrimSamples);
            if (poly.size() >= 3)
                loops.push_back(std::move(poly));
        }
    const double* params = rParams.As<double>();
    for (std::size_t i = 0; i < n; ++i)
        inside[i] = !clip || detail::brep_inside(loops, params[i * 2], params[i * 2 + 1]) ? 1 : 0;
    return out;
}

Mesh iga_sample(const BrepInfo& rInfo, const IgaSampleOptions& rOptions) {
    CadJsonReadOptions opt;
    opt.mSubdivisions = rOptions.mSubdivisions;
    opt.mTrimSamples = rOptions.mTrimSamples;
    opt.mTrim = rOptions.mTrim;
    opt.mEdges = rOptions.mEdges;
    return sample_brep(rInfo, opt);
}

Mesh iga_sample(const NurbsInfo& rInfo, const IgaSampleOptions& rOptions, NurbsInfo* pMapped) {
    if (rOptions.mSubdivisions < 0 || rOptions.mSubdivisions > kIgaMaxSubdivisions)
        iga_fail("subdivisions must be 0 (automatic) to " + std::to_string(kIgaMaxSubdivisions));
    if (rInfo.mPatches.empty())
        iga_fail("the model has no patches");
    if (rOptions.mSubdivisions == 0 && rInfo.mTopology.mDimension != 0) {
        // A model with its MFEM patch topology samples as `read_mfem` does.
        try {
            // The model stays as it is (the sampler regathers its control
            // points from MFEM's numbering, which can move a shared one by an
            // ulp); only the cell maps are the sampling's.
            NurbsInfo derived;
            Mesh mesh = sample_mfem_nurbs(rInfo, pMapped ? &derived : nullptr);
            if (pMapped) {
                *pMapped = rInfo;
                pMapped->mCellPatch = std::move(derived.mCellPatch);
                pMapped->mCellSpan = std::move(derived.mCellSpan);
                pMapped->mBoundaryCellPatch = std::move(derived.mBoundaryCellPatch);
            }
            return mesh;
        } catch (const WriteError& rE) {
            iga_fail(std::string("the patch topology is inconsistent: ") + rE.what());
        }
    }
    const int S = std::max(rOptions.mSubdivisions, 1);
    std::vector<IgaTable> tables;
    int k = 1;
    std::size_t dim = 0, sdim = 0;
    for (std::size_t i = 0; i < rInfo.mPatches.size(); ++i) {
        const NurbsPatch& p = rInfo.mPatches[i];
        const std::string what = "patch " + std::to_string(i);
        tables.push_back(iga_make_table(p, what, iga_is_rational(p)));
        const IgaShape& sh = tables.back().mShape;
        if (i == 0) {
            dim = sh.mDim;
            sdim = sh.mSdim;
        }
        if (sh.mDim != dim || sh.mSdim != sdim)
            iga_fail(what + " differs in parametric or space dimension from patch 0");
        if (p.mFields.size() != rInfo.mFieldNames.size())
            iga_fail(what + " has " + std::to_string(p.mFields.size()) + " fields for " +
                     std::to_string(rInfo.mFieldNames.size()) + " names");
        if (sh.mFieldComps != tables.front().mShape.mFieldComps)
            iga_fail(what + " has fields of different components than patch 0");
        for (std::size_t d = 0; d < dim; ++d)
            k = std::max(k, p.mDegrees[d]);
    }
    const lagrange::Shape shape = dim == 1   ? lagrange::Shape::Line
                                  : dim == 2 ? lagrange::Shape::Quad
                                             : lagrange::Shape::Hexahedron;
    const auto lattice = lagrange::vtk_lattice(shape, k);
    const std::vector<std::size_t>& comps = tables.front().mShape.mFieldComps;
    const std::size_t nv = tables.front().mNet.mColumns;

    struct Pending {
        std::array<std::vector<double>, 3> mG;
        std::array<std::vector<std::size_t>, 3> mSpan;
        std::array<std::size_t, 3> mWidth{1, 1, 1};
        std::vector<std::size_t> mPoints;
    };
    std::vector<Pending> pending(rInfo.mPatches.size());
    std::vector<std::int64_t> conn, cell_patch, cell_span;
    std::size_t num_points = 0, num_cells = 0;
    for (std::size_t pi = 0; pi < rInfo.mPatches.size(); ++pi) {
        const NurbsPatch& p = rInfo.mPatches[pi];
        const IgaShape& sh = tables[pi].mShape;
        Pending& pd = pending[pi];
        std::array<detail::BrepAxis, 3> axis;
        std::array<std::size_t, 3> ncell{1, 1, 1};
        for (std::size_t d = 0; d < dim; ++d) {
            axis[d] = detail::brep_axis(p.mDegrees[d], p.mKnots[d], sh.mCounts[d], S, nullptr);
            ncell[d] = axis[d].mLo.size();
            pd.mG[d] = detail::brep_lattice(axis[d], k);
            pd.mWidth[d] = pd.mG[d].size();
            const auto& sp = axis[d].mSpan;
            for (std::size_t i = 0; i < pd.mG[d].size(); ++i)
                pd.mSpan[d].push_back(
                    sp.empty() ? 0 : sp[std::min(i / static_cast<std::size_t>(k), sp.size() - 1)]);
        }
        std::vector<std::int64_t> local(pd.mWidth[0] * pd.mWidth[1] * pd.mWidth[2], -1);
        for (std::size_t iw = 0; iw < ncell[2]; ++iw)
            for (std::size_t jv = 0; jv < ncell[1]; ++jv)
                for (std::size_t iu = 0; iu < ncell[0]; ++iu) {
                    for (const auto& ijk : lattice) {
                        const std::size_t gi =
                            iu * static_cast<std::size_t>(k) + static_cast<std::size_t>(ijk[0]);
                        const std::size_t gj = dim > 1 ? jv * static_cast<std::size_t>(k) +
                                                             static_cast<std::size_t>(ijk[1])
                                                       : 0;
                        const std::size_t gl = dim > 2 ? iw * static_cast<std::size_t>(k) +
                                                             static_cast<std::size_t>(ijk[2])
                                                       : 0;
                        const std::size_t g = gi + pd.mWidth[0] * (gj + pd.mWidth[1] * gl);
                        if (local[g] < 0) {
                            local[g] = static_cast<std::int64_t>(num_points + pd.mPoints.size());
                            pd.mPoints.push_back(g);
                        }
                        conn.push_back(local[g]);
                    }
                    cell_patch.push_back(static_cast<std::int64_t>(pi));
                    cell_span.push_back(axis[0].mSpanNo[iu]);
                    cell_span.push_back(dim > 1 ? axis[1].mSpanNo[jv] : 0);
                    cell_span.push_back(dim > 2 ? axis[2].mSpanNo[iw] : 0);
                    ++num_cells;
                }
        num_points += pd.mPoints.size();
    }

    const std::size_t cols = sdim == 2 ? 2 : 3;
    NDArray points(DType::Float64, {num_points, cols});
    std::vector<NDArray> fields;
    for (const std::size_t c : comps)
        fields.emplace_back(DType::Float64, c == 1 ? std::vector<std::size_t>{num_points}
                                                   : std::vector<std::size_t>{num_points, c});
    std::size_t offset = 0;
    for (std::size_t pi = 0; pi < rInfo.mPatches.size(); ++pi) {
        const Pending& pd = pending[pi];
        const IgaTable& t = tables[pi];
        double* px = points.As<double>();
        parallel_for(pd.mPoints.size(), [&](std::size_t i) {
            const std::size_t g = pd.mPoints[i];
            const std::size_t gi = g % pd.mWidth[0], gj = (g / pd.mWidth[0]) % pd.mWidth[1],
                              gl = g / (pd.mWidth[0] * pd.mWidth[1]);
            const std::array<std::size_t, 3> at{gi, gj, gl};
            std::array<std::size_t, 3> span{0, 0, 0};
            std::array<double, 3> u{0, 0, 0};
            for (std::size_t d = 0; d < dim; ++d) {
                span[d] = pd.mSpan[d][at[d]];
                u[d] = pd.mG[d][at[d]];
            }
            std::vector<double> row(nv);
            detail::iga_eval_point(t.mNet, span.data(), u.data(), 0, row.data(), nullptr);
            double* out = px + (offset + i) * cols;
            for (std::size_t d = 0; d < cols; ++d)
                out[d] = d < sdim ? row[d] : 0.0;
            std::size_t c0 = sdim;
            for (std::size_t f = 0; f < comps.size(); ++f) {
                double* pf = fields[f].As<double>() + (offset + i) * comps[f];
                for (std::size_t d = 0; d < comps[f]; ++d)
                    pf[d] = row[c0 + d];
                c0 += comps[f];
            }
        });
        offset += pd.mPoints.size();
    }

    Mesh mesh;
    mesh.AssignPoints(std::move(points));
    if (num_cells) {
        NDArray cc(DType::Int64, {num_cells, lattice.size()});
        std::copy(conn.begin(), conn.end(), cc.As<std::int64_t>());
        mesh.AddCellBlock(iga_quad_family(dim, k), std::move(cc));
        NDArray cp(DType::Int64, {num_cells}), cs(DType::Int64, {num_cells, std::size_t{3}});
        std::copy(cell_patch.begin(), cell_patch.end(), cp.As<std::int64_t>());
        std::copy(cell_span.begin(), cell_span.end(), cs.As<std::int64_t>());
        if (pMapped) {
            *pMapped = rInfo;
            pMapped->mCellPatch = cp;
            pMapped->mCellSpan = cs;
            pMapped->mBoundaryCellPatch = NDArray(DType::Int64, {0});
        }
        std::vector<NDArray> a, b;
        a.push_back(std::move(cp));
        b.push_back(std::move(cs));
        mesh.AddCellData(kIgaPatchName, std::move(a));
        mesh.AddCellData(kIgaSpanName, std::move(b));
    } else if (pMapped) {
        *pMapped = rInfo;
        pMapped->mCellPatch = NDArray(DType::Int64, {0});
        pMapped->mCellSpan = NDArray(DType::Int64, {0, 3});
        pMapped->mBoundaryCellPatch = NDArray(DType::Int64, {0});
    }
    for (std::size_t f = 0; f < fields.size(); ++f)
        mesh.AddPointData(rInfo.mFieldNames[f], std::move(fields[f]));
    return mesh;
}

NurbsInfo iga_remove_knots(const NurbsInfo& rInfo, std::size_t Patch, std::int32_t Direction,
                           const std::vector<double>& rKnots, double Tolerance,
                           IgaRemoveKnotsReport* pReport) {
    if (Patch >= rInfo.mPatches.size())
        iga_fail("patch " + std::to_string(Patch) + " is out of range (the model has " +
                 std::to_string(rInfo.mPatches.size()) + ")");
    if (!std::isfinite(Tolerance) || Tolerance < 0.0)
        iga_fail("the tolerance must be finite and not negative");
    const NurbsPatch& patch = rInfo.mPatches[Patch];
    const std::string what = "patch " + std::to_string(Patch);
    const IgaShape sh = iga_check_patch(patch, what);
    iga_check_direction(sh, Direction);
    const bool rational = iga_is_rational(patch);
    if (rational)
        iga_check_weights(patch, what);
    const auto dir = static_cast<std::size_t>(Direction);
    detail::IgaNet net = iga_make_net(patch, sh, rational);

    // The knots to remove: distinct values ascending, each with its count.
    const std::vector<double>& k = net.mKnots[dir];
    const int p = net.mDegrees[dir];
    const double lo = k[static_cast<std::size_t>(p)], hi = k[sh.mCounts[dir]];
    std::vector<double> values;
    std::vector<std::int32_t> asked;
    if (rKnots.empty()) {
        for (std::size_t i = static_cast<std::size_t>(p) + 1; i < sh.mCounts[dir]; ++i) {
            if (!(k[i] > lo && k[i] < hi))
                continue;
            if (!values.empty() && values.back() == k[i])
                ++asked.back();
            else {
                values.push_back(k[i]);
                asked.push_back(1);
            }
        }
    } else {
        std::vector<double> sorted = rKnots;
        std::sort(sorted.begin(), sorted.end());
        for (const double u : sorted) {
            if (!std::isfinite(u))
                iga_fail("a knot to remove is not finite");
            if (!values.empty() && values.back() == u)
                ++asked.back();
            else {
                values.push_back(u);
                asked.push_back(1);
            }
        }
    }
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (!(values[i] > lo && values[i] < hi))
            iga_fail("knot " + std::to_string(values[i]) +
                     " is not an interior knot of direction " + std::to_string(Direction));
        const auto mult = static_cast<std::int32_t>(std::count(k.begin(), k.end(), values[i]));
        if (mult == 0)
            iga_fail("knot " + std::to_string(values[i]) +
                     " is not in the knot vector of direction " + std::to_string(Direction));
        if (asked[i] > mult)
            iga_fail("knot " + std::to_string(values[i]) + " occurs " + std::to_string(mult) +
                     " time(s) in direction " + std::to_string(Direction) + ", not " +
                     std::to_string(asked[i]));
    }

    auto compute_tol = [](const NurbsPatch& rP, const IgaShape& rSh, bool Rat,
                          const detail::IgaNet& rNet, double Tol) {
        double wmin = 1.0, pmax = 0.0, big = 0.0;
        const double* w = rP.mWeights.As<double>();
        const double* x = rP.mControlPoints.As<double>();
        for (std::size_t i = 0; i < rSh.mNumPoints; ++i) {
            if (Rat)
                wmin = i == 0 ? w[0] : std::min(wmin, w[i]);
            double norm = 0.0;
            for (std::size_t d = 0; d < rSh.mSdim; ++d)
                norm += x[i * rSh.mSdim + d] * x[i * rSh.mSdim + d];
            pmax = std::max(pmax, std::sqrt(norm));
        }
        for (const double v : rNet.mValues)
            big = std::max(big, std::abs(v));
        return std::max(Tol * wmin / (1.0 + pmax), 1e-12 * (1.0 + big));
    };

    NurbsInfo out;
    out.mFieldNames = rInfo.mFieldNames;
    out.mTopology = rInfo.mTopology;
    out.mPatches = rInfo.mPatches;

    IgaRemoveKnotsReport report;
    if (rInfo.mTopology.mDimension != 0 && rInfo.mPatches.size() > 1) {
        // Joined patches (roadmap 9.8.7): collect every (patch, direction) sharing
        // this knot vector and remove a knot only to the multiplicity every
        // sharing patch accepts within its tolerance.
        struct SharedDir {
            std::size_t mPatch = 0;
            std::size_t mDir = 0;
            IgaShape mShape;
            bool mRational = false;
            detail::IgaNet mNet;
            double mTol = 0.0;
        };
        std::vector<SharedDir> shared;
        for (std::size_t pi = 0; pi < rInfo.mPatches.size(); ++pi) {
            const NurbsPatch& rp = rInfo.mPatches[pi];
            const std::string pwhat = "patch " + std::to_string(pi);
            const IgaShape rsh = iga_check_patch(rp, pwhat);
            const bool rrat = iga_is_rational(rp);
            if (rrat)
                iga_check_weights(rp, pwhat);
            for (std::size_t d = 0; d < rsh.mDim; ++d) {
                if (rp.mDegrees[d] == p && rp.mKnots[d] == k) {
                    SharedDir sd;
                    sd.mPatch = pi;
                    sd.mDir = d;
                    sd.mShape = rsh;
                    sd.mRational = rrat;
                    sd.mNet = iga_make_net(rp, rsh, rrat);
                    sd.mTol = compute_tol(rp, rsh, rrat, sd.mNet, Tolerance);
                    shared.push_back(std::move(sd));
                }
            }
        }
        for (std::size_t i = 0; i < values.size(); ++i) {
            int agreed = asked[i];
            for (const SharedDir& sd : shared) {
                detail::IgaNet probe = sd.mNet;
                const int done = detail::iga_net_remove_knot(probe, sd.mDir, values[i], agreed,
                                                             sd.mTol);
                agreed = std::min(agreed, done);
                if (agreed == 0)
                    break;
            }
            if (agreed > 0) {
                for (SharedDir& sd : shared)
                    (void)detail::iga_net_remove_knot(sd.mNet, sd.mDir, values[i], agreed, sd.mTol);
            }
            report.mKnots.push_back(values[i]);
            report.mRequested.push_back(asked[i]);
            report.mRemoved.push_back(agreed);
            report.mNumRemoved += agreed;
        }
        for (const SharedDir& sd : shared)
            out.mPatches[sd.mPatch] = iga_patch_from_net(sd.mNet, sd.mShape, sd.mRational);
    } else {
        const double tol = compute_tol(patch, sh, rational, net, Tolerance);
        for (std::size_t i = 0; i < values.size(); ++i) {
            const int done = detail::iga_net_remove_knot(net, dir, values[i], asked[i], tol);
            report.mKnots.push_back(values[i]);
            report.mRequested.push_back(asked[i]);
            report.mRemoved.push_back(done);
            report.mNumRemoved += done;
        }
        out.mPatches[Patch] = iga_patch_from_net(net, sh, rational);
    }

    out.mCellPatch = NDArray(DType::Int64, {0});
    out.mCellSpan = NDArray(DType::Int64, {0, 3});
    out.mBoundaryCellPatch = NDArray(DType::Int64, {0});
    if (pReport)
        *pReport = std::move(report);
    return out;
}

BrepInfo iga_remove_knots(const BrepInfo& rInfo, std::size_t Face, std::int32_t Direction,
                          const std::vector<double>& rKnots, double Tolerance,
                          IgaRemoveKnotsReport* pReport) {
    if (Face >= rInfo.mFaces.size())
        iga_fail("face " + std::to_string(Face) + " is out of range (the model has " +
                 std::to_string(rInfo.mFaces.size()) + ")");
    const std::int64_t pid = rInfo.mFaces[Face].mPatch;
    if (pid < 0 || static_cast<std::size_t>(pid) >= rInfo.mPatches.size())
        iga_fail("face " + std::to_string(Face) + " names no patch");
    NurbsInfo single;
    single.mPatches.push_back(rInfo.mPatches[static_cast<std::size_t>(pid)]);
    NurbsInfo reduced = iga_remove_knots(single, 0, Direction, rKnots, Tolerance, pReport);
    BrepInfo out = rInfo;
    out.mPatches[static_cast<std::size_t>(pid)] = std::move(reduced.mPatches[0]);
    return out;
}

IgaClosestPoint iga_closest_point(const NurbsInfo& rInfo, std::size_t Patch,
                                  const NDArray& rPoints) {
    if (Patch >= rInfo.mPatches.size())
        iga_fail("patch " + std::to_string(Patch) + " is out of range (the model has " +
                 std::to_string(rInfo.mPatches.size()) + ")");
    NurbsPatch patch = rInfo.mPatches[Patch];
    const std::string what = "patch " + std::to_string(Patch);
    patch.mFields.clear();  // only the geometry takes part
    const bool rational = iga_is_rational(patch);
    const IgaTable t = iga_make_table(patch, what, rational);
    const IgaShape& sh = t.mShape;
    const std::size_t dim = sh.mDim, sdim = sh.mSdim;
    if (dim < 1 || dim > 3)
        iga_fail("closest_point expects 1, 2 or 3 parametric directions");
    if (rPoints.Dtype() != DType::Float64 || rPoints.Shape().size() != 2 ||
        rPoints.Shape()[1] != sdim)
        iga_fail("points are not Float64 (N, " + std::to_string(sdim) + ")");
    const std::size_t n = rPoints.Shape()[0];
    const double* q = rPoints.As<double>();
    for (std::size_t i = 0; i < n * sdim; ++i)
        if (!std::isfinite(q[i]))
            iga_fail("a query point has a non-finite coordinate");

    std::array<double, 3> lo{0, 0, 0}, hi{0, 0, 0};
    for (std::size_t d = 0; d < dim; ++d) {
        lo[d] = patch.mKnots[d][static_cast<std::size_t>(patch.mDegrees[d])];
        hi[d] = patch.mKnots[d][sh.mCounts[d]];
    }
    const auto eval = [&](const double* pU, double* pX, double* pD1, double* pD2) {
        std::array<std::size_t, 3> span{0, 0, 0};
        std::array<double, 3> u{0, 0, 0};
        for (std::size_t d = 0; d < dim; ++d) {
            u[d] = std::min(std::max(pU[d], lo[d]), hi[d]);
            span[d] =
                detail::nurbs_find_span(patch.mDegrees[d], patch.mKnots[d], sh.mCounts[d], u[d]);
        }
        if (pD2)
            detail::iga_eval_point_d2(t.mNet, span.data(), u.data(), pX, pD1, pD2);
        else
            detail::iga_eval_point(t.mNet, span.data(), u.data(), 0, pX, nullptr);
    };

    // The seed lattice: degree + 1 samples per nonzero knot span.
    std::array<std::vector<double>, 3> axis;
    for (std::size_t d = 0; d < dim; ++d) {
        const std::vector<double>& k = patch.mKnots[d];
        std::size_t spans = 0;
        for (std::size_t s = static_cast<std::size_t>(patch.mDegrees[d]); s < sh.mCounts[d]; ++s)
            spans += k[s] < k[s + 1] ? 1 : 0;
        if (spans == 0)
            iga_fail("closest_point knot vector in direction " + std::to_string(d) +
                     " has zero active span length");
        const std::size_t cap = dim == 1 ? 4096 : (dim == 2 ? 64 : 16);
        std::size_t per = static_cast<std::size_t>(patch.mDegrees[d]) + 1;
        per = std::max<std::size_t>(1, std::min(per, cap > spans ? (cap - 1) / spans : 1));
        axis[d].push_back(lo[d]);
        for (std::size_t s = static_cast<std::size_t>(patch.mDegrees[d]); s < sh.mCounts[d]; ++s)
            if (k[s] < k[s + 1])
                for (std::size_t j = 1; j <= per; ++j)
                    axis[d].push_back(j == per ? k[s + 1]
                                               : k[s] + (k[s + 1] - k[s]) * static_cast<double>(j) /
                                                            static_cast<double>(per));
    }
    const std::size_t n0 = axis[0].size();
    const std::size_t n1 = dim >= 2 ? axis[1].size() : 1;
    const std::size_t n2 = dim == 3 ? axis[2].size() : 1;
    const std::size_t nseeds = n0 * n1 * n2;
    std::vector<double> seed_x(nseeds * sdim);
    for (std::size_t c = 0; c < n2; ++c)
        for (std::size_t b = 0; b < n1; ++b)
            for (std::size_t a = 0; a < n0; ++a) {
                const double u[3] = {axis[0][a], dim >= 2 ? axis[1][b] : 0.0,
                                     dim == 3 ? axis[2][c] : 0.0};
                const std::size_t idx = a + n0 * (b + n1 * c);
                eval(u, seed_x.data() + idx * sdim, nullptr, nullptr);
            }

    IgaClosestPoint out;
    out.mParams = NDArray(DType::Float64, {n, dim});
    out.mPoints = NDArray(DType::Float64, {n, sdim});
    out.mDistances = NDArray(DType::Float64, {n});
    double* op = out.mParams.As<double>();
    double* ox = out.mPoints.As<double>();
    double* od = out.mDistances.As<double>();
    parallel_for(n, [&](std::size_t i) {
        const double* p = q + i * sdim;
        // The three nearest samples, ascending distance then index.
        std::array<std::size_t, 3> best{0, 0, 0};
        std::array<double, 3> bestd{1e300, 1e300, 1e300};
        for (std::size_t s = 0; s < nseeds; ++s) {
            double d2 = 0.0;
            for (std::size_t d = 0; d < sdim; ++d)
                d2 += (seed_x[s * sdim + d] - p[d]) * (seed_x[s * sdim + d] - p[d]);
            for (std::size_t r = 0; r < 3; ++r)
                if (d2 < bestd[r]) {
                    for (std::size_t m = 2; m > r; --m) {
                        best[m] = best[m - 1];
                        bestd[m] = bestd[m - 1];
                    }
                    best[r] = s;
                    bestd[r] = d2;
                    break;
                }
        }
        double win_f = 1e300;
        std::array<double, 3> win_u{lo[0], dim >= 2 ? lo[1] : 0.0, dim == 3 ? lo[2] : 0.0};
        std::array<double, 3> x{0, 0, 0};
        std::array<double, 9> dd1{};
        std::array<double, 27> dd2{};
        for (std::size_t r = 0; r < 3; ++r) {
            if (bestd[r] >= 1e300)
                break;
            const std::size_t s = best[r];
            const std::size_t a = s % n0;
            const std::size_t b = (s / n0) % n1;
            const std::size_t c = s / (n0 * n1);
            std::array<double, 3> u{axis[0][a], dim >= 2 ? axis[1][b] : 0.0,
                                    dim == 3 ? axis[2][c] : 0.0};
            const auto residual = [&](const std::array<double, 3>& rU, double* pX, double* pRes) {
                eval(rU.data(), pX, nullptr, nullptr);
                double f = 0.0;
                for (std::size_t d = 0; d < sdim; ++d) {
                    pRes[d] = pX[d] - p[d];
                    f += 0.5 * pRes[d] * pRes[d];
                }
                return f;
            };
            std::array<double, 3> res{0, 0, 0};
            double f = residual(u, x.data(), res.data());
            for (int iter = 0; iter < 80; ++iter) {
                eval(u.data(), x.data(), dd1.data(), dd2.data());
                double g[3] = {0, 0, 0};
                double H[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
                double G[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
                for (std::size_t d = 0; d < sdim; ++d)
                    res[d] = x[d] - p[d];
                for (std::size_t k = 0; k < dim; ++k) {
                    for (std::size_t d = 0; d < sdim; ++d)
                        g[k] += res[d] * dd1[k * sdim + d];
                    for (std::size_t l = 0; l < dim; ++l)
                        for (std::size_t d = 0; d < sdim; ++d) {
                            G[k][l] += dd1[k * sdim + d] * dd1[l * sdim + d];
                            H[k][l] += dd1[k * sdim + d] * dd1[l * sdim + d] +
                                       res[d] * dd2[(k * dim + l) * sdim + d];
                        }
                }
                std::size_t free_idx[3];
                std::size_t nfree = 0;
                for (std::size_t k = 0; k < dim; ++k) {
                    if (!((u[k] <= lo[k] && g[k] > 0.0) || (u[k] >= hi[k] && g[k] < 0.0)))
                        free_idx[nfree++] = k;
                }
                if (nfree == 0)
                    break;
                double step[3] = {0, 0, 0};
                auto solve_spd = [&](double A[3][3]) -> bool {
                    double L[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
                    for (std::size_t ri = 0; ri < nfree; ++ri) {
                        for (std::size_t ci = 0; ci <= ri; ++ci) {
                            double sum = A[free_idx[ri]][free_idx[ci]];
                            for (std::size_t mi = 0; mi < ci; ++mi)
                                sum -= L[ri][mi] * L[ci][mi];
                            if (ri == ci) {
                                if (!(sum > 0.0))
                                    return false;
                                L[ri][ri] = std::sqrt(sum);
                            } else {
                                L[ri][ci] = sum / L[ci][ci];
                            }
                        }
                    }
                    double y[3] = {0, 0, 0};
                    for (std::size_t ri = 0; ri < nfree; ++ri) {
                        double sum = -g[free_idx[ri]];
                        for (std::size_t mi = 0; mi < ri; ++mi)
                            sum -= L[ri][mi] * y[mi];
                        y[ri] = sum / L[ri][ri];
                    }
                    double z[3] = {0, 0, 0};
                    for (std::size_t ri = nfree; ri-- > 0;) {
                        double sum = y[ri];
                        for (std::size_t mi = ri + 1; mi < nfree; ++mi)
                            sum -= L[mi][ri] * z[mi];
                        z[ri] = sum / L[ri][ri];
                    }
                    for (std::size_t ri = 0; ri < nfree; ++ri)
                        step[free_idx[ri]] = z[ri];
                    return true;
                };
                if (!solve_spd(H) && !solve_spd(G)) {
                    double scale = 1e-300;
                    for (std::size_t ri = 0; ri < nfree; ++ri)
                        scale += G[free_idx[ri]][free_idx[ri]];
                    for (std::size_t ri = 0; ri < nfree; ++ri)
                        step[free_idx[ri]] = -g[free_idx[ri]] / scale;
                }
                double alpha = 1.0;
                bool moved = false;
                std::array<double, 3> trial = u;
                std::array<double, 3> xt{0, 0, 0}, rt{0, 0, 0};
                double ft = f;
                for (int half = 0; half < 40 && !moved; ++half, alpha *= 0.5) {
                    for (std::size_t k = 0; k < dim; ++k)
                        trial[k] = std::min(std::max(u[k] + alpha * step[k], lo[k]), hi[k]);
                    ft = residual(trial, xt.data(), rt.data());
                    moved = ft <= f + 1e-14 * f;
                }
                if (!moved)
                    break;
                double change = 0.0;
                for (std::size_t k = 0; k < dim; ++k)
                    change = std::max(change, std::abs(trial[k] - u[k]) / (hi[k] - lo[k]));
                u = trial;
                x = xt;
                f = ft;
                if (change < 1e-15)
                    break;
            }
            if (f < win_f) {
                win_f = f;
                win_u = u;
            }
        }
        eval(win_u.data(), x.data(), nullptr, nullptr);
        double dist2 = 0.0;
        for (std::size_t d = 0; d < sdim; ++d) {
            ox[i * sdim + d] = x[d];
            dist2 += (x[d] - p[d]) * (x[d] - p[d]);
        }
        for (std::size_t k = 0; k < dim; ++k)
            op[i * dim + k] = win_u[k];
        od[i] = std::sqrt(dist2);
    });
    return out;
}

IgaClosestPoint iga_closest_point(const BrepInfo& rInfo, std::size_t Face, const NDArray& rPoints) {
    if (Face >= rInfo.mFaces.size())
        iga_fail("face " + std::to_string(Face) + " is out of range (the model has " +
                 std::to_string(rInfo.mFaces.size()) + ")");
    const std::int64_t pid = rInfo.mFaces[Face].mPatch;
    if (pid < 0 || static_cast<std::size_t>(pid) >= rInfo.mPatches.size())
        iga_fail("face " + std::to_string(Face) + " names no patch");
    NurbsInfo single;
    single.mPatches.push_back(rInfo.mPatches[static_cast<std::size_t>(pid)]);
    return iga_closest_point(single, 0, rPoints);
}

Mesh iga_bezier_extract(const NurbsInfo& rInfo, std::size_t Patch) {
    if (Patch >= rInfo.mPatches.size())
        iga_fail("patch " + std::to_string(Patch) + " is out of range (the model has " +
                 std::to_string(rInfo.mPatches.size()) + ")");
    const NurbsPatch& patch = rInfo.mPatches[Patch];
    const std::string what = "patch " + std::to_string(Patch);
    const IgaShape sh = iga_check_patch(patch, what);
    if (patch.mFields.size() != rInfo.mFieldNames.size())
        iga_fail(what + " has " + std::to_string(patch.mFields.size()) + " fields for " +
                 std::to_string(rInfo.mFieldNames.size()) + " names");
    const int p = patch.mDegrees[0];
    for (std::size_t d = 1; d < sh.mDim; ++d)
        if (patch.mDegrees[d] != p)
            iga_fail(what +
                     " has unequal degrees; Bezier cells take one order, elevate with "
                     "iga_elevate_degree first");
    const bool rational = iga_is_rational(patch);
    if (rational)
        iga_check_weights(patch, what);
    detail::IgaNet net = iga_make_net(patch, sh, rational);

    // Raise every interior knot to multiplicity p.
    for (std::size_t d = 0; d < sh.mDim; ++d) {
        const std::vector<double>& k = net.mKnots[d];
        std::vector<double> ins;
        std::size_t i = static_cast<std::size_t>(p) + 1;
        while (i < sh.mCounts[d]) {
            std::size_t j = i;
            while (j < k.size() && k[j] == k[i])
                ++j;
            const int mult = static_cast<int>(j - i);
            for (int m = mult; m < p; ++m)
                ins.push_back(k[i]);
            i = j;
        }
        if (!ins.empty())
            detail::iga_net_insert(net, d, ins);
    }

    // Elements of the refined net: spans of nonzero length.
    std::array<std::vector<std::size_t>, 3> elem;
    std::array<std::vector<std::int64_t>, 3> elem_no;
    for (std::size_t d = 0; d < 3; ++d) {
        if (d >= sh.mDim) {
            elem[d].push_back(static_cast<std::size_t>(p));
            elem_no[d].push_back(0);
            continue;
        }
        const std::vector<double>& k = net.mKnots[d];
        std::int64_t no = 0;
        for (std::size_t s = static_cast<std::size_t>(p); s < net.mCounts[d]; ++s)
            if (k[s] < k[s + 1]) {
                elem[d].push_back(s);
                elem_no[d].push_back(no++);
            }
    }
    const std::size_t dim = sh.mDim, sdim = sh.mSdim;
    const lagrange::Shape shape = dim == 1   ? lagrange::Shape::Line
                                  : dim == 2 ? lagrange::Shape::Quad
                                             : lagrange::Shape::Hexahedron;
    const auto lattice = lagrange::vtk_lattice(shape, p);
    const std::size_t nn = lattice.size();
    const std::size_t ncell = elem[0].size() * elem[1].size() * elem[2].size();
    const std::size_t npts = ncell * nn;
    const std::size_t C = net.mColumns;
    const std::size_t cols = sdim == 2 ? 2 : 3;

    NDArray points(DType::Float64, {npts, cols});
    NDArray weights(DType::Float64, {npts});
    std::vector<NDArray> fields;
    for (const std::size_t c : sh.mFieldComps)
        fields.emplace_back(DType::Float64, c == 1 ? std::vector<std::size_t>{npts}
                                                   : std::vector<std::size_t>{npts, c});
    NDArray conn(DType::Int64, {ncell, nn});
    NDArray degrees(DType::Float64, {ncell, std::size_t{3}});
    NDArray cpatch(DType::Int64, {ncell});
    NDArray cspan(DType::Int64, {ncell, std::size_t{3}});
    std::int64_t* pc = conn.As<std::int64_t>();
    double* pdeg = degrees.As<double>();
    double* px = points.As<double>();
    double* pw = weights.As<double>();
    std::size_t cell = 0;
    for (std::size_t iw = 0; iw < elem[2].size(); ++iw)
        for (std::size_t jv = 0; jv < elem[1].size(); ++jv)
            for (std::size_t iu = 0; iu < elem[0].size(); ++iu, ++cell) {
                for (std::size_t n = 0; n < nn; ++n) {
                    const std::size_t g0 = elem[0][iu] - static_cast<std::size_t>(p) +
                                           static_cast<std::size_t>(lattice[n][0]);
                    const std::size_t g1 = dim > 1 ? elem[1][jv] - static_cast<std::size_t>(p) +
                                                         static_cast<std::size_t>(lattice[n][1])
                                                   : 0;
                    const std::size_t g2 = dim > 2 ? elem[2][iw] - static_cast<std::size_t>(p) +
                                                         static_cast<std::size_t>(lattice[n][2])
                                                   : 0;
                    const double* row =
                        net.mValues.data() + (g0 + net.mCounts[0] * (g1 + net.mCounts[1] * g2)) * C;
                    const std::size_t at = cell * nn + n;
                    pc[at] = static_cast<std::int64_t>(at);
                    const double w = rational ? row[C - 1] : 1.0;
                    pw[at] = w;
                    for (std::size_t d = 0; d < cols; ++d)
                        px[at * cols + d] = d < sdim ? (rational ? row[d] / w : row[d]) : 0.0;
                    std::size_t c = sdim;
                    for (std::size_t f = 0; f < fields.size(); ++f) {
                        double* pf = fields[f].As<double>() + at * sh.mFieldComps[f];
                        for (std::size_t d = 0; d < sh.mFieldComps[f]; ++d, ++c)
                            pf[d] = rational ? row[c] / w : row[c];
                    }
                }
                for (std::size_t d = 0; d < 3; ++d)
                    pdeg[cell * 3 + d] = d < dim ? static_cast<double>(p) : 0.0;
                cpatch.As<std::int64_t>()[cell] = static_cast<std::int64_t>(Patch);
                cspan.As<std::int64_t>()[cell * 3 + 0] = elem_no[0][iu];
                cspan.As<std::int64_t>()[cell * 3 + 1] = dim > 1 ? elem_no[1][jv] : 0;
                cspan.As<std::int64_t>()[cell * 3 + 2] = dim > 2 ? elem_no[2][iw] : 0;
            }

    Mesh mesh;
    mesh.AssignPoints(std::move(points));
    const char* type = dim == 1   ? "VTK_BEZIER_CURVE"
                       : dim == 2 ? "VTK_BEZIER_QUADRILATERAL"
                                  : "VTK_BEZIER_HEXAHEDRON";
    mesh.AddCellBlock(type, std::move(conn));
    std::vector<NDArray> a, b, c;
    a.push_back(std::move(degrees));
    b.push_back(std::move(cpatch));
    c.push_back(std::move(cspan));
    mesh.AddCellData("HigherOrderDegrees", std::move(a));
    mesh.AddCellData(kIgaPatchName, std::move(b));
    mesh.AddCellData(kIgaSpanName, std::move(c));
    if (rational)
        mesh.AddPointData("RationalWeights", std::move(weights));
    for (std::size_t f = 0; f < fields.size(); ++f)
        mesh.AddPointData(rInfo.mFieldNames[f], std::move(fields[f]));
    return mesh;
}

Mesh iga_bezier_extract(const BrepInfo& rInfo, std::size_t Face) {
    if (Face >= rInfo.mFaces.size())
        iga_fail("face " + std::to_string(Face) + " is out of range (the model has " +
                 std::to_string(rInfo.mFaces.size()) + ")");
    const std::int64_t pid = rInfo.mFaces[Face].mPatch;
    if (pid < 0 || static_cast<std::size_t>(pid) >= rInfo.mPatches.size())
        iga_fail("face " + std::to_string(Face) + " names no patch");
    NurbsInfo single;
    single.mPatches.push_back(rInfo.mPatches[static_cast<std::size_t>(pid)]);
    return iga_bezier_extract(single, 0);
}

namespace {

struct IgaExtraction1D {
    std::vector<std::vector<double>> mC;  // each (p+1)x(p+1) row-major: N^e = C^e * B
    std::vector<std::size_t> mSpanStart;  // first B-spline basis index for each element
};

IgaExtraction1D iga_bezier_extraction_1d(const std::vector<double>& rU, int P) {
    const auto p = static_cast<std::size_t>(P);
    const std::size_t m = rU.size();
    const std::size_t np = p + 1;
    auto eye = [&]() {
        std::vector<double> mat(np * np, 0.0);
        for (std::size_t i = 0; i < np; ++i)
            mat[i * np + i] = 1.0;
        return mat;
    };

    IgaExtraction1D res;
    std::size_t a = p;
    std::size_t b = a + 1;
    std::vector<double> cur = eye();
    std::vector<double> alphas(np, 0.0);
    while (b < m) {
        std::vector<double> next = eye();
        const std::size_t i = b;
        while (b + 1 < m && rU[b + 1] == rU[b])
            ++b;
        const std::size_t mult = b - i + 1;
        res.mSpanStart.push_back(a - p);
        if (mult < p) {
            const double numer = rU[b] - rU[a];
            for (std::size_t j = p; j > mult; --j)
                alphas[j - mult - 1] = numer / (rU[a + j] - rU[a]);
            const std::size_t r = p - mult;
            for (std::size_t j = 1; j <= r; ++j) {
                const std::size_t save = r - j;
                const std::size_t s = mult + j;
                for (std::size_t k = p; k >= s; --k) {
                    const double alpha = alphas[k - s];
                    for (std::size_t row = 0; row < np; ++row)
                        cur[row * np + k] =
                            alpha * cur[row * np + k] + (1.0 - alpha) * cur[row * np + (k - 1)];
                }
                if (b + 1 < m) {
                    for (std::size_t l = 0; l <= j; ++l)
                        next[(save + l) * np + save] = cur[(p - j + l) * np + p];
                }
            }
        }
        res.mC.push_back(std::move(cur));
        if (b + 1 < m) {
            a = b;
            ++b;
            cur = std::move(next);
        } else {
            break;
        }
    }
    return res;
}

}  // namespace

IgaBezierExtraction iga_bezier_extract_operators(const NurbsInfo& rInfo, std::size_t Patch) {
    IgaBezierExtraction out;
    out.mMesh = iga_bezier_extract(rInfo, Patch);
    const NurbsPatch& patch = rInfo.mPatches[Patch];
    const IgaShape sh = iga_check_patch(patch, "patch " + std::to_string(Patch));
    const std::size_t dim = sh.mDim;
    const int p = patch.mDegrees[0];
    const auto np = static_cast<std::size_t>(p + 1);

    std::array<IgaExtraction1D, 3> ex;
    for (std::size_t d = 0; d < 3; ++d) {
        if (d < dim) {
            ex[d] = iga_bezier_extraction_1d(patch.mKnots[d], p);
        } else {
            ex[d].mC.push_back({1.0});
            ex[d].mSpanStart.push_back(0);
        }
    }
    const std::size_t nu = dim >= 1 ? np : 1;
    const std::size_t nv = dim >= 2 ? np : 1;
    const std::size_t nw = dim >= 3 ? np : 1;
    const std::size_t nn = nu * nv * nw;

    for (std::size_t ew = 0; ew < ex[2].mC.size(); ++ew) {
        for (std::size_t ev = 0; ev < ex[1].mC.size(); ++ev) {
            for (std::size_t eu = 0; eu < ex[0].mC.size(); ++eu) {
                const auto& cu = ex[0].mC[eu];
                const auto& cv = ex[1].mC[ev];
                const auto& cw = ex[2].mC[ew];
                NDArray op = NDArray::Uninit(DType::Float64, {nn, nn});
                double* po = op.As<double>();
                for (std::size_t rw = 0; rw < nw; ++rw)
                    for (std::size_t rv = 0; rv < nv; ++rv)
                        for (std::size_t ru = 0; ru < nu; ++ru) {
                            const std::size_t row = ru + nu * (rv + nv * rw);
                            for (std::size_t kw = 0; kw < nw; ++kw)
                                for (std::size_t kv = 0; kv < nv; ++kv)
                                    for (std::size_t ku = 0; ku < nu; ++ku) {
                                        const std::size_t col = ku + nu * (kv + nv * kw);
                                        po[row * nn + col] = cu[ru * nu + ku] * cv[rv * nv + kv] *
                                                             cw[rw * nw + kw];
                                    }
                        }
                NDArray cp = NDArray::Uninit(DType::Int64, {nn});
                std::int64_t* pcp = cp.As<std::int64_t>();
                const std::size_t su = ex[0].mSpanStart[eu];
                const std::size_t sv = ex[1].mSpanStart[ev];
                const std::size_t sw = ex[2].mSpanStart[ew];
                for (std::size_t rw = 0; rw < nw; ++rw)
                    for (std::size_t rv = 0; rv < nv; ++rv)
                        for (std::size_t ru = 0; ru < nu; ++ru) {
                            const std::size_t row = ru + nu * (rv + nv * rw);
                            const std::size_t g =
                                (su + ru) + sh.mCounts[0] * ((sv + rv) + sh.mCounts[1] * (sw + rw));
                            pcp[row] = static_cast<std::int64_t>(g);
                        }
                out.mOperators.push_back(std::move(op));
                out.mElementControlPoints.push_back(std::move(cp));
            }
        }
    }
    return out;
}

IgaBezierExtraction iga_bezier_extract_operators(const BrepInfo& rInfo, std::size_t Face) {
    if (Face >= rInfo.mFaces.size())
        iga_fail("face " + std::to_string(Face) + " is out of range (the model has " +
                 std::to_string(rInfo.mFaces.size()) + ")");
    const std::int64_t pid = rInfo.mFaces[Face].mPatch;
    if (pid < 0 || static_cast<std::size_t>(pid) >= rInfo.mPatches.size())
        iga_fail("face " + std::to_string(Face) + " names no patch");
    NurbsInfo single;
    single.mPatches.push_back(rInfo.mPatches[static_cast<std::size_t>(pid)]);
    return iga_bezier_extract_operators(single, 0);
}

}  // namespace meshioplusplus

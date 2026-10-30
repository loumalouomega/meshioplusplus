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
//                   meshio++ default license: LICENSE
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//
// fit_cad. See operations/fit_cad.hpp for the contract.

// System includes
#include <string>

// System includes
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

// Project includes
#include "meshioplusplus/operations/fit_cad.hpp"
#include "meshioplusplus/providers.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/step.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/operations/refine_patches.hpp"
#include "meshioplusplus/operations/segment.hpp"

// Project includes (private, not installed)
#include "../detail/closed_facets.hpp"
#include "../detail/fit_cad_native.hpp"
#include "../detail/cut_graph.hpp"
#include "../detail/freeform_adaptive.hpp"
#include "../detail/patch_flatten.hpp"
#include "../detail/patch_network.hpp"
#include "../detail/primitive_constraints.hpp"
#include "../detail/quad_layout.hpp"
#include "../detail/step_brep.hpp"
#include "../detail/surface_rings.hpp"

// System includes
#include <array>
#include <functional>
#include <map>
#include <numeric>
#include <queue>
#include <set>
#include <tuple>
#include <unordered_map>
#include <unordered_set>


namespace meshioplusplus {

bool fit_cad_has_occt() noexcept {
    const auto p_kernel = cad_kernel();
    return p_kernel && p_kernel->mFitCad;
}

const char* fit_cad_occt_version() noexcept {
    // The registry retains every provider it ever held, so this stays valid.
    const auto p_kernel = cad_kernel();
    return p_kernel ? p_kernel->mVersion.c_str() : "";
}

namespace {

constexpr const char* kFcnPrefix = "fit_cad: ";

// The share of a triangulated surface's area on facet faces (`fit_cad:type`
// 0): `rTri` three point ids per triangle, `rType` per triangle. Summed in
// triangle order, so the value does not depend on the thread count.
double fcn_facet_area_fraction(const std::vector<double>& rXyz,
                               const std::vector<std::int64_t>& rTri,
                               const std::vector<std::int64_t>& rType) {
    double total = 0.0, facet = 0.0;
    for (std::size_t t = 0; t < rType.size(); ++t) {
        const double* a = &rXyz[static_cast<std::size_t>(rTri[3 * t]) * 3];
        const double* b = &rXyz[static_cast<std::size_t>(rTri[3 * t + 1]) * 3];
        const double* c = &rXyz[static_cast<std::size_t>(rTri[3 * t + 2]) * 3];
        const double u[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
        const double v[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
        const double n[3] = {u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2],
                             u[0] * v[1] - u[1] * v[0]};
        const double area = 0.5 * std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
        total += area;
        if (rType[t] == 0)
            facet += area;
    }
    return total > 0.0 ? facet / total : 0.0;
}

// `fit_cad:face_table`: one row per record, the `FitCadFaceRecord` columns in
// declaration order as float64.
NDArray fcn_face_table(const std::vector<FitCadFaceRecord>& rFaces) {
    constexpr std::size_t kCols = 13;
    NDArray a = NDArray::Uninit(DType::Float64, {rFaces.size(), kCols});
    double* p = a.As<double>();
    for (std::size_t i = 0; i < rFaces.size(); ++i) {
        const FitCadFaceRecord& f = rFaces[i];
        const double row[kCols] = {static_cast<double>(f.mPatch),
                                   static_cast<double>(f.mFace),
                                   static_cast<double>(f.mKind),
                                   static_cast<double>(f.mParamMethod),
                                   static_cast<double>(f.mDemotionReason),
                                   static_cast<double>(f.mNumControl),
                                   static_cast<double>(f.mIterations),
                                   static_cast<double>(f.mRefinements),
                                   static_cast<double>(f.mEscalationSteps),
                                   f.mMaxDist,
                                   f.mRmsDist,
                                   f.mReverseMaxDist,
                                   f.mMaxG1AngleDeg};
        std::copy(row, row + kCols, p + i * kCols);
    }
    return a;
}

// The annotated skin: segmented by the surfaces its facets lie on and fitted
// when `segment:patch` is missing, fitted when the `fit:patch_*` table is.
// Returns a pointer to the mesh to work on, owned by rOwned or the input.
const Mesh* fcn_annotate(const Mesh& rMesh, Mesh& rOwned) {
    if (rMesh.HasCellData(kSegmentPatchName)) {
        if (rMesh.HasFieldData("fit:patch_id"))
            return &rMesh;
        log::info("fit_cad: no fit:patch_* table; running fit_primitives");
        rOwned = fit_primitives(rMesh).mMesh;
        return &rOwned;
    }
    log::info(
        "fit_cad: no segment:patch labels; running segment (method 'primitive') and "
        "fit_primitives");
    SegmentOptions seg;
    seg.mMethod = SegmentMethod::Primitive;
    Mesh segmented = segment(rMesh, seg).mMesh;
    // Fit at the tolerance the patches were grown at: an absolute one, the
    // segmentation's share of the diagonal.
    double lo[3], hi[3];
    for (int k = 0; k < 3; ++k) {
        lo[k] = std::numeric_limits<double>::infinity();
        hi[k] = -lo[k];
    }
    const std::vector<double> xyz = detail::surface_xyz(rMesh);
    for (std::size_t i = 0; i + 2 < xyz.size(); i += 3)
        for (std::size_t k = 0; k < 3; ++k) {
            lo[k] = std::min(lo[k], xyz[i + k]);
            hi[k] = std::max(hi[k], xyz[i + k]);
        }
    double diag = 0.0;
    for (int k = 0; k < 3; ++k)
        diag += hi[k] > lo[k] ? (hi[k] - lo[k]) * (hi[k] - lo[k]) : 0.0;
    diag = std::sqrt(diag);
    FitPrimitivesOptions fit;
    if (diag > 0.0) {
        fit.mTolerance = seg.mTolerance * diag;
        fit.mRelative = false;
    }
    rOwned = fit_primitives(segmented, fit).mMesh;
    return &rOwned;
}

// What the refinement and the primitive constraints did, for the report.
using detail::FcnRefined;

// The tolerance a face may keep from its patch before a writer demotes it,
// from the table alone: the requested one, or three times the largest patch
// tolerance (`fcn_tolerance` without its bounding-box fallback).
double fcn_face_tolerance(const Mesh& rMesh, double Requested) {
    if (Requested > 0.0)
        return Requested;
    double largest = 0.0;
    if (rMesh.HasFieldData("fit:patch_tolerance")) {
        const NDArray& tols = rMesh.FieldData("fit:patch_tolerance");
        for (std::size_t r = 0; r < tols.Size(); ++r) {
            const double t = detail::read_double(tols, r);
            if (std::isfinite(t))
                largest = std::max(largest, t);
        }
    }
    return 3.0 * largest;
}

// The annotated skin, refined by `refine_patches` and its primitives
// constrained when asked.
const Mesh* fcn_prepare(const Mesh& rMesh, Mesh& rOwned, const FitCadOptions& rOptions,
                        FcnRefined& rRefined) {
    const Mesh* work = fcn_annotate(rMesh, rOwned);
    if (rOptions.mRefineSegmentation) {
        rRefined.mResult = refine_patches(*work, rOptions.mRefine);
        rRefined.mRan = true;
        log::info("fit_cad: refined the patches: {} -> {} of the area explained by primitives",
                  rRefined.mResult.mInitialExplainedAreaFraction,
                  rRefined.mResult.mExplainedAreaFraction);
        rOwned = std::move(rRefined.mResult.mMesh);
        work = &rOwned;
    }
    if (rOptions.mConstrainPrimitives) {
        const double face_tol = fcn_face_tolerance(*work, rOptions.mTolerance);
        if (face_tol > 0.0) {
            if (work != &rOwned) {
                rOwned = detail::clone_mesh(
                    *work, [](DataLocation, const std::string&, std::string&) { return true; });
                work = &rOwned;
            }
            rRefined.mConstraints = detail::constrain_primitives(rOwned, face_tol);
            log::info("fit_cad: constrained {} primitive(s) in {} group(s)",
                      rRefined.mConstraints.mNumConstrained, rRefined.mConstraints.mNumGroups);
        }
    }
    return work;
}

// The refinement's report and its `refine:*` field data on the result.
void fcn_report_refined(const FcnRefined& rRefined, const Mesh& rWork, FitCadResult& rResult) {
    rResult.mNumConstrainedPrimitives = rRefined.mConstraints.mNumConstrained;
    rResult.mNumConstraintGroups = rRefined.mConstraints.mNumGroups;
    if (!rRefined.mRan)
        return;
    const RefinePatchesResult& r = rRefined.mResult;
    rResult.mRefined = true;
    rResult.mRefinementKept = true;
    rResult.mRefineIterations = r.mNumIterations;
    rResult.mNumRefineSplits = r.mNumSplits;
    rResult.mNumRefineMerges = r.mNumMerges;
    rResult.mInitialExplainedAreaFraction = r.mInitialExplainedAreaFraction;
    rResult.mExplainedAreaFraction = r.mExplainedAreaFraction;
    for (const char* name :
         {kRefinePatchIdName, kRefinePatchProxyName, kRefinePatchIterationsName,
          kRefinePatchRmsName, kRefinePatchMaxName, kRefinePatchCellsName,
          kRefinePatchDemotionName, kRefineExplainedAreaName, kRefineHausdorffName})
        if (rWork.HasFieldData(name))
            rResult.mMesh.AddFieldData(name, rWork.FieldData(name));
}

// The kernel-free backend: the fitted patches as analytic faces and the rest
// as coplanar polygons, written by the B-rep `step` writer.
// What the native backend's adaptive freeform pass (v17.12.0) built: the
// faces for the writer, the mesh to write (its unexplained patches split into
// discs where that was needed), and a record per patch it tried.
struct FcnFreeform {
    std::vector<detail::StepFreeformFace> mFaces;
    std::map<std::int64_t, FitCadFaceRecord> mRecords;  // by label
    Mesh mMesh;  // with the pieces' labels and fit-table rows, when split
    bool mSplit = false;
    std::int64_t mNumSplit = 0;
    std::int64_t mNumResegmented = 0;  // regions the re-segmentation rung cut
    // The patch network's (network and auto modes): its patches' labels, its
    // shared edges and extraordinary vertices, and its largest G1 angle.
    std::set<std::int64_t> mNetworkLabels;
    std::vector<detail::StepEdgeCurve> mEdgeCurves;
    std::int64_t mNumNetworkEdges = 0;
    std::int64_t mNumExtraordinary = 0;
    double mMaxG1AngleDeg = 0.0;
};

// The tolerance the native writer will apply: as `step_brep_write` takes it.
double fcn_tolerance(const Mesh& rMesh, const detail::ClosedFacets& rFc, double Requested) {
    if (Requested > 0.0)
        return Requested;
    double largest = 0.0;
    if (rMesh.HasFieldData("fit:patch_type") && rMesh.HasFieldData("fit:patch_tolerance") &&
        rMesh.HasFieldData("fit:patch_params")) {
        const NDArray& types = rMesh.FieldData("fit:patch_type");
        const NDArray& tols = rMesh.FieldData("fit:patch_tolerance");
        const NDArray& params = rMesh.FieldData("fit:patch_params");
        const std::size_t np = types.Size() > 0 ? params.Size() / types.Size() : 0;
        for (std::size_t r = 0; r < types.Size() && r < tols.Size(); ++r) {
            const std::int64_t t = detail::read_int(types, r);
            bool finite = true;
            for (std::size_t k = 0; k < np; ++k)
                finite = finite && std::isfinite(detail::read_double(params, r * np + k));
            const double pt = detail::read_double(tols, r);
            if (t >= 1 && t <= 5 && finite && std::isfinite(pt))
                largest = std::max(largest, pt);
        }
    }
    if (largest > 0.0)
        return 3.0 * largest;
    double lo[3] = {1e300, 1e300, 1e300}, hi[3] = {-1e300, -1e300, -1e300};
    for (std::int64_t v : rFc.mNodes)
        for (int j = 0; j < 3; ++j) {
            lo[j] = std::min(lo[j], rFc.mXyz[static_cast<std::size_t>(v) * 3 + j]);
            hi[j] = std::max(hi[j], rFc.mXyz[static_cast<std::size_t>(v) * 3 + j]);
        }
    return 1e-3 * std::sqrt((hi[0] - lo[0]) * (hi[0] - lo[0]) + (hi[1] - lo[1]) * (hi[1] - lo[1]) +
                            (hi[2] - lo[2]) * (hi[2] - lo[2]));
}

// Facets as local triangles (polygons fanned), points renumbered in the order
// they are first used; with `pGlobal`, each local point's skin point.
void fcn_local_triangles(const detail::ClosedFacets& rFc, const std::vector<std::size_t>& rFacets,
                         std::vector<detail::PatchVec>& rPts, std::vector<detail::PatchTri>& rTris,
                         std::vector<std::size_t>& rTriFacet, std::vector<std::int64_t>* pGlobal) {
    std::unordered_map<std::int64_t, std::size_t> local;
    const auto vid = [&](std::int64_t V) {
        auto ins = local.emplace(V, rPts.size());
        if (ins.second) {
            const std::size_t v = static_cast<std::size_t>(V);
            rPts.push_back({rFc.mXyz[v * 3], rFc.mXyz[v * 3 + 1], rFc.mXyz[v * 3 + 2]});
            if (pGlobal)
                pGlobal->push_back(V);
        }
        return ins.first->second;
    };
    for (std::size_t f : rFacets) {
        const std::int64_t b = rFc.mStart[f], n = rFc.mStart[f + 1] - b;
        const std::size_t a0 = vid(rFc.mNodes[static_cast<std::size_t>(b)]);
        for (std::int64_t k = 1; k + 1 < n; ++k) {
            rTris.push_back({a0, vid(rFc.mNodes[static_cast<std::size_t>(b + k)]),
                             vid(rFc.mNodes[static_cast<std::size_t>(b + k + 1)])});
            rTriFacet.push_back(f);
        }
    }
}

// A face's sense: its surface normal against the facets', area-weighted;
// `rU`, `rV` the local points' parameters on it.
bool fcn_sense(const detail::ClosedFacets& rFc, const detail::BsplineSurface& rS,
               const std::vector<detail::PatchTri>& rTris,
               const std::vector<std::size_t>& rTriFacet, const std::vector<double>& rU,
               const std::vector<double>& rV) {
    double agree = 0.0;
    for (std::size_t t = 0; t < rTris.size(); ++t) {
        const std::size_t f = rTriFacet[t];
        double u = 0.0, v = 0.0;
        for (std::size_t c : rTris[t]) {
            u += rU[c] / 3.0;
            v += rV[c] / 3.0;
        }
        detail::BsplineEval e;
        detail::bspline_eval(rS, u, v, 1, e);
        const double n[3] = {e.mSu[1] * e.mSv[2] - e.mSu[2] * e.mSv[1],
                             e.mSu[2] * e.mSv[0] - e.mSu[0] * e.mSv[2],
                             e.mSu[0] * e.mSv[1] - e.mSu[1] * e.mSv[0]};
        agree += n[0] * rFc.mNormals[f * 3] + n[1] * rFc.mNormals[f * 3 + 1] +
                 n[2] * rFc.mNormals[f * 3 + 2];
    }
    return agree >= 0.0;
}

// A freeform region's patch network (network and auto modes): the items it
// joins, its facets as local triangles, its quad layout (`mLayout.mOk` false:
// none, and no network), its fitted network, and per layout patch its facets
// (each facet with its first triangle's patch).
using detail::FcnRegionNetwork;

// The items (patches, as facet lists) grouped into regions -- items sharing a
// facet edge -- and each laid out and fitted as one network, in the order of
// each region's first item. Free of OpenCASCADE: both backends build on it.
std::vector<FcnRegionNetwork> fcn_networks(const detail::ClosedFacets& rFc,
                                           const std::vector<std::vector<std::size_t>>& rItems,
                                           const FitCadOptions& rOptions, double Tol) {
    double diag = 0.0;
    {
        double lo[3] = {1e300, 1e300, 1e300}, hi[3] = {-1e300, -1e300, -1e300};
        for (std::int64_t v : rFc.mNodes)
            for (int j = 0; j < 3; ++j) {
                lo[j] = std::min(lo[j], rFc.mXyz[static_cast<std::size_t>(v) * 3 + j]);
                hi[j] = std::max(hi[j], rFc.mXyz[static_cast<std::size_t>(v) * 3 + j]);
            }
        diag = std::sqrt((hi[0] - lo[0]) * (hi[0] - lo[0]) + (hi[1] - lo[1]) * (hi[1] - lo[1]) +
                         (hi[2] - lo[2]) * (hi[2] - lo[2]));
    }
    std::map<std::pair<std::int64_t, std::int64_t>, std::size_t> edge_item;
    std::vector<std::size_t> parent(rItems.size());
    std::iota(parent.begin(), parent.end(), std::size_t{0});
    const auto find = [&](std::size_t A) {
        while (parent[A] != A)
            A = parent[A] = parent[parent[A]];
        return A;
    };
    for (std::size_t i = 0; i < rItems.size(); ++i)
        for (std::size_t f : rItems[i]) {
            const std::int64_t b = rFc.mStart[f], n = rFc.mStart[f + 1] - b;
            for (std::int64_t k = 0; k < n; ++k) {
                const std::int64_t a = rFc.mNodes[static_cast<std::size_t>(b + k)];
                const std::int64_t c = rFc.mNodes[static_cast<std::size_t>(b + (k + 1) % n)];
                const auto key = std::make_pair(std::min(a, c), std::max(a, c));
                const auto ins = edge_item.emplace(key, i);
                if (!ins.second) {
                    const std::size_t x = find(ins.first->second), y = find(i);
                    if (x != y)
                        parent[std::max(x, y)] = std::min(x, y);
                }
            }
        }
    std::map<std::size_t, std::vector<std::size_t>> regions;
    for (std::size_t i = 0; i < rItems.size(); ++i)
        regions[find(i)].push_back(i);
    std::vector<FcnRegionNetwork> out;
    for (const auto& kv : regions) {
        FcnRegionNetwork reg;
        reg.mItems = kv.second;
        std::vector<std::size_t> facets;
        for (std::size_t i : kv.second)
            facets.insert(facets.end(), rItems[i].begin(), rItems[i].end());
        std::sort(facets.begin(), facets.end());
        fcn_local_triangles(rFc, facets, reg.mPts, reg.mTris, reg.mTriFacet, &reg.mGlobal);
        detail::QuadLayoutOptions lo;
        lo.mTargetSize = rOptions.mLayoutTargetPatchSize * diag;
        lo.mMaxValenceDeviation = static_cast<int>(rOptions.mLayoutMaxValenceDeviation);
        lo.mSnapToFeatures = rOptions.mLayoutSnapToFeatures;
        reg.mLayout = detail::quad_layout(reg.mPts, reg.mTris, lo);
        if (!reg.mLayout.mOk) {
            log::debug("fit_cad: no quad layout of a {}-triangle freeform region (failure {})",
                       reg.mTris.size(), static_cast<int>(reg.mLayout.mFailure));
            out.push_back(std::move(reg));
            continue;
        }
        detail::PatchNetworkOptions no;
        no.mTolerance = Tol;
        no.mFreeformTolerance = rOptions.mFreeformTolerance;
        no.mMaxPolesPerDirection = static_cast<std::size_t>(rOptions.mMaxFreeformPoles);
        no.mMaxTotalPoles = static_cast<std::size_t>(rOptions.mFreeformMaxControlPoints);
        no.mG1ToleranceDeg = rOptions.mFreeformG1ToleranceDeg;
        no.mExtraordinary = rOptions.mExtraordinaryVertex;
        no.mValidate = rOptions.mValidateBySurfaceSampling;
        reg.mNet = detail::patch_network_fit(reg.mPts, reg.mTris, reg.mLayout, no);
        // The tangent-plane treatment of the extraordinary vertices kept only
        // where it does better there than leaving them alone (`g0`) and fits
        // as much: next to a vertex where the patches' corners meet at very
        // unequal angles, its constraint can bend the edges between them.
        if (no.mExtraordinary != FitCadExtraordinaryVertex::G0 && reg.mNet.mNumExtraordinary > 0) {
            detail::PatchNetworkOptions g0 = no;
            g0.mExtraordinary = FitCadExtraordinaryVertex::G0;
            detail::PatchNetwork plain =
                detail::patch_network_fit(reg.mPts, reg.mTris, reg.mLayout, g0);
            const auto fitted = [&](const detail::PatchNetwork& rN) {
                std::size_t n = 0;
                for (std::size_t i = 0; i < reg.mLayout.mPatches.size(); ++i)
                    n += rN.mPatches[i].mOk ? reg.mLayout.mPatches[i].mTris.size() : 0;
                return n;
            };
            if (plain.mMaxG1AngleEvDeg < reg.mNet.mMaxG1AngleEvDeg &&
                fitted(plain) >= fitted(reg.mNet)) {
                log::debug(
                    "fit_cad: the extraordinary vertices of a {}-triangle region left G0: "
                    "{} degrees next to them instead of {}",
                    reg.mTris.size(), plain.mMaxG1AngleEvDeg, reg.mNet.mMaxG1AngleEvDeg);
                reg.mNet = std::move(plain);
            }
        }
        log::debug(
            "fit_cad: patch network of a {}-triangle freeform region: {} patch(es), {} "
            "edge(s), {} extraordinary vertex(es), G1 {} degrees",
            reg.mTris.size(), reg.mLayout.mPatches.size(), reg.mNet.mNumEdges,
            reg.mNet.mNumExtraordinary, reg.mNet.mMaxG1AngleDeg);
        reg.mPatchFacets.assign(reg.mLayout.mPatches.size(), {});
        for (std::size_t t = 0; t < reg.mTris.size(); ++t)
            if (t == 0 || reg.mTriFacet[t] != reg.mTriFacet[t - 1])
                reg.mPatchFacets[static_cast<std::size_t>(reg.mLayout.mPatchOfTri[t])].push_back(
                    reg.mTriFacet[t]);
        out.push_back(std::move(reg));
    }
    return out;
}

// The record of an accepted network patch.
FitCadFaceRecord fcn_network_record(std::int64_t Label, const detail::NetworkPatch& rP) {
    FitCadFaceRecord rec;
    rec.mPatch = Label;
    rec.mParamMethod = FitCadParamMethod::Square;
    rec.mNumControl = static_cast<std::int64_t>(rP.mSurface.NumPoles());
    rec.mIterations = rP.mIterations;
    rec.mRefinements = rP.mRefinements;
    rec.mMaxDist = rP.mMaxDist;
    rec.mRmsDist = rP.mRmsDist;
    rec.mReverseMaxDist = rP.mReverseMaxDist;
    rec.mMaxG1AngleDeg = rP.mMaxG1AngleDeg;
    return rec;
}

// The re-segmentation rung: a region the adaptive faces and the network left
// as facets cut into pieces. Where edges sharper than 30 degrees already cut
// it, the pieces between them; otherwise `K` pieces by variational shape
// approximation on the L2,1 normal metric (Cohen-Steiner, Alliez and Desbrun
// 2004, the energy refine_patches' Lloyd step uses): proxies seeded farthest
// apart in normal, grown together through one priority queue that crosses an
// edge (in steps of ten degrees) only once every gentler one is taken, so the
// pieces meet at the sharpest ridge between their seeds, re-seeded at the
// facet nearest their proxy, for at most ten Lloyd steps, then the cuts
// smoothed. A cell fanned into several facets stays whole. Empty when no cut
// leaves every piece at least four facets.
// `OnlyAcrossCreases`: none for a region with fewer than three sharp edges
// on curves that close or end on its boundary.
std::vector<std::vector<std::size_t>> fcn_resegment(const detail::ClosedFacets& rFc,
                                                    const std::vector<std::size_t>& rFacets,
                                                    std::size_t K, bool OnlyAcrossCreases = false) {
    const std::size_t n = rFacets.size();
    if (K < 2 || n < 4 * K)
        return {};
    const auto normal = [&](std::size_t I) { return &rFc.mNormals[rFacets[I] * 3]; };
    std::vector<double> area(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t f = rFacets[i];
        const std::int64_t b = rFc.mStart[f], m = rFc.mStart[f + 1] - b;
        const double* p0 =
            &rFc.mXyz[static_cast<std::size_t>(rFc.mNodes[static_cast<std::size_t>(b)]) * 3];
        for (std::int64_t k = 1; k + 1 < m; ++k) {
            const double* p1 =
                &rFc.mXyz[static_cast<std::size_t>(rFc.mNodes[static_cast<std::size_t>(b + k)]) *
                          3];
            const double* p2 = &rFc.mXyz[static_cast<std::size_t>(
                                             rFc.mNodes[static_cast<std::size_t>(b + k + 1)]) *
                                         3];
            const double u[3] = {p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};
            const double v[3] = {p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]};
            const double c[3] = {u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2],
                                 u[0] * v[1] - u[1] * v[0]};
            area[i] += 0.5 * std::sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]);
        }
    }
    // Neighbours across shared edges, with the edge's dihedral angle in
    // tens of degrees (from 3, 30 degrees, it is a crease).
    static constexpr double kCosTens[9] = {
        0.98480775301220806, 0.93969262078590838, 0.86602540378443865,
        0.76604444311897804, 0.64278760968653933, 0.5,
        0.34202014332566873, 0.17364817766693035, 0.0};
    std::vector<std::vector<std::pair<std::size_t, std::uint8_t>>> adj(n);
    std::size_t num_creases = 0;
    {
        std::unordered_map<std::uint64_t, std::size_t> first;
        std::vector<std::pair<std::uint64_t, std::uint64_t>> crease_edges;
        for (std::size_t i = 0; i < n; ++i) {
            const std::size_t f = rFacets[i];
            const std::int64_t b = rFc.mStart[f], m = rFc.mStart[f + 1] - b;
            for (std::int64_t k = 0; k < m; ++k) {
                const std::uint64_t a =
                    static_cast<std::uint64_t>(rFc.mNodes[static_cast<std::size_t>(b + k)]);
                const std::uint64_t c = static_cast<std::uint64_t>(
                    rFc.mNodes[static_cast<std::size_t>(b + (k + 1) % m)]);
                const std::uint64_t key = (std::min(a, c) << 32) | std::max(a, c);
                const auto ins = first.emplace(key, i);
                if (ins.second || ins.first->second == i)
                    continue;
                const std::size_t j = ins.first->second;
                const double* ni = normal(i);
                const double* nj = normal(j);
                const double dot = ni[0] * nj[0] + ni[1] * nj[1] + ni[2] * nj[2];
                std::uint8_t tens = 0;
                while (tens < 9 && dot < kCosTens[tens])
                    ++tens;
                adj[i].emplace_back(j, tens);
                adj[j].emplace_back(i, tens);
                ins.first->second = n;  // matched: an inner edge
                if (tens >= 3)
                    crease_edges.emplace_back(a, c);
            }
        }
        // A crease must run somewhere: a group of connected sharp edges that
        // neither reaches the region's boundary nor closes (a coarse or noisy
        // skin's scattered ones) is not counted; one that starts at the
        // boundary and fades out inside (a crease partly under the angle) is.
        std::unordered_set<std::uint64_t> on_boundary;
        for (const auto& kv : first)
            if (kv.second != n) {
                on_boundary.insert(kv.first >> 32);
                on_boundary.insert(kv.first & 0xffffffffu);
            }
        std::unordered_map<std::uint64_t, std::uint64_t> parent;
        const auto find = [&](std::uint64_t V) {
            auto it = parent.emplace(V, V).first;
            while (it->second != it->first) {
                auto up = parent.find(it->second);
                it->second = parent.find(up->second)->first;
                it = up;
            }
            return it->first;
        };
        for (const auto& e : crease_edges) {
            const std::uint64_t ra = find(e.first), rb = find(e.second);
            if (ra != rb)
                parent[std::max(ra, rb)] = std::min(ra, rb);
        }
        std::unordered_map<std::uint64_t, std::size_t> edges_of, verts_of;
        std::unordered_set<std::uint64_t> touches;
        for (const auto& e : crease_edges)
            ++edges_of[find(e.first)];
        for (const auto& kv : parent) {
            const std::uint64_t r = find(kv.first);
            ++verts_of[r];
            if (on_boundary.count(kv.first))
                touches.insert(r);
        }
        for (const auto& kv : edges_of)
            if (touches.count(kv.first) || kv.second >= verts_of[kv.first])  // or a cycle
                num_creases += kv.second;
    }
    // A fanned cell takes its first facet's piece.
    const auto pieces = [&](std::vector<std::size_t> Label, std::size_t NumPieces) {
        std::unordered_map<std::int64_t, std::size_t> of_cell;
        for (std::size_t i = 0; i < n; ++i)
            Label[i] = of_cell.emplace(rFc.mCell[rFacets[i]], Label[i]).first->second;
        std::vector<std::vector<std::size_t>> out(NumPieces);
        for (std::size_t i = 0; i < n; ++i)
            out[Label[i]].push_back(rFacets[i]);
        for (const auto& piece : out)
            if (piece.size() < 4)
                return std::vector<std::vector<std::size_t>>{};
        return out;
    };
    // A region whose creases already cut it: its pieces between them.
    {
        std::vector<std::size_t> comp(n, n);
        std::size_t num = 0;
        for (std::size_t i = 0; i < n; ++i) {
            if (comp[i] != n)
                continue;
            std::vector<std::size_t> stack{i};
            comp[i] = num;
            while (!stack.empty()) {
                const std::size_t a = stack.back();
                stack.pop_back();
                for (const auto& [j, tens] : adj[a])
                    if (tens < 3 && comp[j] == n) {
                        comp[j] = num;
                        stack.push_back(j);
                    }
            }
            ++num;
        }
        // A piece under four facets (a corner between two creases) joins the
        // piece it touches first, then the pieces are numbered again.
        std::vector<std::size_t> size(num, 0);
        for (std::size_t c : comp)
            ++size[c];
        for (bool changed = true; changed;) {
            changed = false;
            for (std::size_t i = 0; i < n; ++i)
                if (size[comp[i]] < 4)
                    for (const auto& nb : adj[i])
                        if (comp[nb.first] != comp[i] && size[comp[nb.first]] >= 4) {
                            const std::size_t from = comp[i], to = comp[nb.first];
                            for (std::size_t& c : comp)
                                if (c == from)
                                    c = to;
                            size[to] += size[from];
                            size[from] = 0;
                            changed = true;
                            break;
                        }
        }
        std::vector<std::size_t> renum(num, num);
        std::size_t kept = 0;
        for (std::size_t& c : comp) {
            if (renum[c] == num)
                renum[c] = kept++;
            c = renum[c];
        }
        num = kept;
        if (num >= 2) {
            auto cut = pieces(std::move(comp), num);
            if (!cut.empty())
                return cut;
        }
    }
    if (OnlyAcrossCreases && num_creases < 3)
        return {};
    const auto error = [&](std::size_t I, const std::array<double, 3>& rN) {
        const double* q = normal(I);
        const double d[3] = {q[0] - rN[0], q[1] - rN[1], q[2] - rN[2]};
        return area[I] * (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
    };
    // Seeds: farthest from the mean normal, then farthest from those taken.
    std::array<double, 3> mean = {0.0, 0.0, 0.0};
    for (std::size_t i = 0; i < n; ++i)
        for (int j = 0; j < 3; ++j)
            mean[j] += area[i] * normal(i)[j];
    const double len = std::sqrt(mean[0] * mean[0] + mean[1] * mean[1] + mean[2] * mean[2]);
    if (!(len > 0.0))
        return {};
    for (double& c : mean)
        c /= len;
    std::vector<std::array<double, 3>> proxy;
    std::vector<std::size_t> seed;
    std::vector<double> nearest(n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
        nearest[i] = error(i, mean) / std::max(area[i], 1e-300);
    while (seed.size() < K) {
        const std::size_t s = static_cast<std::size_t>(
            std::max_element(nearest.begin(), nearest.end()) - nearest.begin());
        if (!(nearest[s] > 1e-12))
            return {};
        seed.push_back(s);
        const double* q = normal(s);
        proxy.push_back({q[0], q[1], q[2]});
        for (std::size_t i = 0; i < n; ++i)
            nearest[i] = std::min(nearest[i], error(i, proxy.back()) / std::max(area[i], 1e-300));
    }
    std::vector<std::size_t> label(n);
    for (int step = 0; step < 10; ++step) {
        // Grow: (sharpest edge crossed, error, facet, proxy), smallest first,
        // so the pieces meet at the sharpest ridge between their seeds.
        using Entry = std::tuple<std::uint8_t, double, std::size_t, std::size_t>;
        std::priority_queue<Entry, std::vector<Entry>, std::greater<Entry>> queue;
        std::vector<std::size_t> grown(n, K);
        for (std::size_t p = 0; p < K; ++p)
            queue.emplace(std::uint8_t{0}, 0.0, seed[p], p);
        while (!queue.empty()) {
            const auto [sharpest, e, i, p] = queue.top();
            queue.pop();
            if (grown[i] != K)
                continue;
            grown[i] = p;
            for (const auto& [j, tens] : adj[i])
                if (grown[j] == K)
                    queue.emplace(std::max(sharpest, tens), error(j, proxy[p]), j, p);
        }
        // A facet no seed reaches (a region in pieces) joins the first proxy.
        for (std::size_t& g : grown)
            if (g == K)
                g = 0;
        const bool same = step > 0 && grown == label;
        label = std::move(grown);
        if (same)
            break;
        // Refit the proxies, and re-seed each at its facet nearest to it.
        for (std::size_t p = 0; p < K; ++p) {
            std::array<double, 3> s = {0.0, 0.0, 0.0};
            for (std::size_t i = 0; i < n; ++i)
                if (label[i] == p)
                    for (int j = 0; j < 3; ++j)
                        s[j] += area[i] * normal(i)[j];
            const double l = std::sqrt(s[0] * s[0] + s[1] * s[1] + s[2] * s[2]);
            if (!(l > 0.0))
                return {};
            for (double& c : s)
                c /= l;
            proxy[p] = s;
            double best = std::numeric_limits<double>::infinity();
            for (std::size_t i = 0; i < n; ++i)
                if (label[i] == p && error(i, s) / std::max(area[i], 1e-300) < best) {
                    best = error(i, s) / std::max(area[i], 1e-300);
                    seed[p] = i;
                }
        }
    }
    // Smooth the cuts: a facet with more sides on one other piece than on its
    // own moves there (a few passes), so no piece keeps a jagged spur.
    for (int pass = 0; pass < 4; ++pass) {
        bool changed = false;
        for (std::size_t i = 0; i < n; ++i) {
            std::vector<std::size_t> votes(K, 0);
            for (const auto& nb : adj[i])
                ++votes[label[nb.first]];
            for (std::size_t p = 0; p < K; ++p)
                if (p != label[i] && votes[p] > votes[label[i]] && votes[p] >= 2) {
                    label[i] = p;
                    changed = true;
                    break;
                }
        }
        if (!changed)
            break;
    }
    return pieces(std::move(label), K);
}

// The adaptive freeform pass of the native backend: every patch the fit table
// leaves unexplained (type 0, or a refinement's freeform or facets verdict)
// fitted by detail/freeform_adaptive; a patch that is not a disc cut into
// discs first (detail/cut_graph), one that is but misses halved, up to
// `mMaxFreeformSplits` times. The pieces get new labels after the largest.
FcnFreeform fcn_freeform(const Mesh& rSource, const detail::ClosedFacets& rFc,
                         const std::vector<std::int64_t>& rCandidates,
                         const FitCadOptions& rOptions, double Tol) {
    FcnFreeform out;
    if (!rSource.HasCellData("segment:patch") || !rSource.HasFieldData("fit:patch_id") ||
        !rSource.HasFieldData("fit:patch_type"))
        return out;
    const NDArray& ids = rSource.FieldData("fit:patch_id");
    const NDArray& types = rSource.FieldData("fit:patch_type");
    const std::vector<std::int64_t> bases = detail::block_bases(rSource);
    std::vector<std::int64_t> cell_label(static_cast<std::size_t>(detail::total_cells(bases)), -1);
    for (std::size_t b = 0; b < rSource.NumCellBlocks(); ++b) {
        const NDArray& a = rSource.CellData("segment:patch", b);
        for (std::size_t c = 0; c < rSource.Cells(b).NumCells() && c < a.Size(); ++c)
            cell_label[static_cast<std::size_t>(bases[b]) + c] = detail::read_int(a, c);
    }
    std::map<std::int64_t, std::vector<std::size_t>> facets_of;  // by label
    for (std::size_t f = 0; f < rFc.NumFacets(); ++f)
        facets_of[cell_label[static_cast<std::size_t>(rFc.mCell[f])]].push_back(f);
    std::int64_t next_label = 0;
    for (std::size_t r = 0; r < ids.Size(); ++r)
        next_label = std::max(next_label, detail::read_int(ids, r) + 1);

    // A work list of (label, facets, halvings so far); pieces join its end.
    struct Item {
        std::int64_t mLabel;
        std::vector<std::size_t> mFacets;
        std::int64_t mDepth;
        std::int64_t mEscalations;
        std::int64_t mCuts = 0;  // early re-segmentation cuts in its line
    };
    std::vector<std::pair<std::int64_t, std::vector<std::size_t>>> pieces;  // new labels
    std::map<std::int64_t, double> piece_tol;

    // The patch network (network and auto modes): the given patches grouped
    // into regions, each laid out and fitted as one network; a patch of it
    // within the tolerance becomes a face under a new label, and the rest
    // come back (as new labels and their facets).
    const auto network =
        [&](const std::vector<std::pair<std::int64_t, std::vector<std::size_t>>>& rItems) {
            std::vector<std::pair<std::int64_t, std::vector<std::size_t>>> left;
            std::vector<std::vector<std::size_t>> item_facets;
            for (const auto& it : rItems)
                item_facets.push_back(it.second);
            for (const FcnRegionNetwork& reg : fcn_networks(rFc, item_facets, rOptions, Tol)) {
                if (!reg.mLayout.mOk) {
                    for (std::size_t i : reg.mItems)
                        left.push_back(rItems[i]);
                    continue;
                }
                const detail::QuadLayout& layout = reg.mLayout;
                const detail::PatchNetwork& net = reg.mNet;
                bool any = false;
                for (std::size_t i = 0; i < layout.mPatches.size(); ++i) {
                    const std::vector<std::size_t>& pf = reg.mPatchFacets[i];
                    if (pf.empty())
                        continue;
                    const detail::NetworkPatch& np = net.mPatches[i];
                    const std::int64_t label = next_label++;
                    pieces.emplace_back(label, pf);
                    if (!np.mOk) {
                        left.emplace_back(label, pf);
                        continue;
                    }
                    any = true;
                    std::vector<detail::PatchTri> ptris;
                    std::vector<std::size_t> pfac;
                    for (std::size_t t : layout.mPatches[i].mTris) {
                        ptris.push_back(reg.mTris[t]);
                        pfac.push_back(reg.mTriFacet[t]);
                    }
                    detail::StepFreeformFace face;
                    face.mLabel = label;
                    face.mSurface = np.mSurface;
                    face.mSense = fcn_sense(rFc, np.mSurface, ptris, pfac, np.mU, np.mV);
                    face.mNetwork = true;
                    out.mFaces.push_back(std::move(face));
                    out.mRecords[label] = fcn_network_record(label, np);
                    out.mNetworkLabels.insert(label);
                }
                if (any) {
                    // Its shared edge curves, for the writer.
                    for (std::size_t c = 0; c < layout.mChains.size(); ++c) {
                        if (net.mChains[c].mPoles.size() < 2)
                            continue;
                        detail::StepEdgeCurve ec;
                        for (std::size_t v : layout.mChains[c].mVerts)
                            ec.mVerts.push_back(reg.mGlobal[v]);
                        ec.mParams = net.mChainParams[c];
                        ec.mCurve = net.mChains[c];
                        out.mEdgeCurves.push_back(std::move(ec));
                    }
                    out.mNumNetworkEdges += static_cast<std::int64_t>(net.mNumEdges);
                    out.mNumExtraordinary += static_cast<std::int64_t>(net.mNumExtraordinary);
                    out.mMaxG1AngleDeg = std::max(out.mMaxG1AngleDeg, net.mMaxG1AngleDeg);
                }
            }
            return left;
        };

    std::vector<Item> work;
    if (rOptions.mFreeformMode == FitCadFreeformMode::Network) {
        std::vector<std::pair<std::int64_t, std::vector<std::size_t>>> items;
        for (std::int64_t label : rCandidates) {
            const auto it = facets_of.find(label);
            if (it != facets_of.end())
                items.emplace_back(label, it->second);
        }
        for (auto& it : network(items))
            work.push_back({it.first, std::move(it.second), 0, 1, rOptions.mResegmentRounds});
    } else {
        for (std::int64_t label : rCandidates) {
            const auto it = facets_of.find(label);
            if (it != facets_of.end())
                work.push_back({label, it->second, 0, 0});
        }
    }
    detail::FreeformAdaptiveOptions o;
    o.mTolerance = Tol;
    o.mFreeformTolerance = rOptions.mFreeformTolerance;
    o.mSmoothing = rOptions.mFreeformSmoothing;
    o.mMaxPoles = static_cast<std::size_t>(rOptions.mMaxFreeformPoles);
    o.mMaxIterations = static_cast<int>(rOptions.mMaxFreeformIterations);
    o.mMaxTotalPoles = static_cast<std::size_t>(rOptions.mFreeformMaxControlPoints);
    o.mParameterization = static_cast<detail::ParamLadder>(rOptions.mFreeformParameterization);
    o.mValidate = rOptions.mValidateBySurfaceSampling;
    o.mThb = rOptions.mLocalRefine == FitCadLocalRefine::Thb;
    // The work list from `From` on: each item a face, or cut into pieces that
    // join the list's end, or facets with its reason.
    const auto run = [&](std::size_t From) {
        for (std::size_t w = From; w < work.size(); ++w) {
            const Item item = work[w];
            // The patch as local triangles (polygons fanned), points renumbered.
            std::vector<detail::PatchVec> pts;
            std::vector<detail::PatchTri> tris;
            std::vector<std::size_t> tri_facet;
            fcn_local_triangles(rFc, item.mFacets, pts, tris, tri_facet, nullptr);
            const detail::FreeformAdaptiveResult fit = detail::freeform_adaptive_fit(pts, tris, o);
            log::debug(
                "fit_cad: adaptive freeform fit of patch {} ({} triangles): {}, {} pole(s), vertex "
                "{} / centroid {} / samples {}, allowance {} (noise {}), {} escalation(s)",
                item.mLabel, tris.size(), fit.mOk ? "accepted" : "rejected",
                fit.mSurface.NumPoles(), fit.mMaxDist, fit.mMaxCentroidDist, fit.mReverseMaxDist,
                fit.mAllowance, fit.mSigma, fit.mEscalations);
            FitCadFaceRecord rec;
            rec.mPatch = item.mLabel;
            rec.mParamMethod = fit.mParam;
            rec.mNumControl = static_cast<std::int64_t>(fit.mNumDofs > 0 ? fit.mNumDofs
                                                                         : fit.mSurface.NumPoles());
            rec.mIterations = fit.mIterations;
            rec.mRefinements = fit.mRefinements;
            rec.mEscalationSteps = item.mEscalations + fit.mEscalations;
            rec.mMaxDist = fit.mMaxDist;
            rec.mRmsDist = fit.mRmsDist;
            rec.mReverseMaxDist = fit.mReverseMaxDist;
            if (fit.mOk) {
                detail::StepFreeformFace face;
                face.mLabel = item.mLabel;
                face.mSurface = fit.mSurface;
                face.mSense = fcn_sense(rFc, fit.mSurface, tris, tri_facet, fit.mU, fit.mV);
                out.mFaces.push_back(std::move(face));
                out.mRecords[item.mLabel] = rec;
                continue;
            }
            // Not a face: cut into discs, else halve, else facets.
            std::size_t num = 0;
            std::vector<std::int64_t> piece;
            bool cut_here = false;
            if (item.mCuts < rOptions.mResegmentRounds) {
                // A patch across sharp edges is cut along them first (at most
                // `mResegmentRounds` times in its line): the re-segmentation
                // rung, where the segmentation left a crease. Its pieces keep
                // the halvings left.
                const auto cut = fcn_resegment(rFc, item.mFacets, 2, true);
                if (!cut.empty()) {
                    std::unordered_map<std::size_t, std::int64_t> of_facet;
                    for (std::size_t k = 0; k < cut.size(); ++k)
                        for (std::size_t f : cut[k])
                            of_facet[f] = static_cast<std::int64_t>(k);
                    piece.resize(tris.size());
                    for (std::size_t t = 0; t < tris.size(); ++t)
                        piece[t] = of_facet[tri_facet[t]];
                    num = cut.size();
                    cut_here = true;
                    ++out.mNumResegmented;
                }
            }
            if (num < 2 && item.mDepth < rOptions.mMaxFreeformSplits && tris.size() >= 32) {
                piece = detail::patch_split_to_discs(pts, tris, 24, num);
                if (num < 2)
                    piece = detail::patch_bisect(pts, tris, num);
            }
            if (num < 2) {
                rec.mDemotionReason = fit.mReason == FitCadDemotionReason::None
                                          ? FitCadDemotionReason::Other
                                          : fit.mReason;
                out.mRecords[item.mLabel] = rec;
                continue;
            }
            std::vector<std::vector<std::size_t>> pf(num);
            for (std::size_t t = 0; t < tris.size(); ++t) {
                if (t > 0 && tri_facet[t] == tri_facet[t - 1])
                    continue;
                pf[static_cast<std::size_t>(piece[t])].push_back(tri_facet[t]);
            }
            for (std::vector<std::size_t>& facets : pf) {
                const std::int64_t label = next_label++;
                work.push_back({label, facets, cut_here ? item.mDepth : item.mDepth + 1,
                                rec.mEscalationSteps + 1, item.mCuts + (cut_here ? 1 : 0)});
                pieces.emplace_back(label, std::move(facets));
            }
            out.mNumSplit += static_cast<std::int64_t>(num);
        }
    };
    run(0);
    // Auto: what the adaptive faces leave goes through the network.
    if (rOptions.mFreeformMode == FitCadFreeformMode::Auto) {
        std::vector<std::pair<std::int64_t, std::vector<std::size_t>>> items;
        for (const Item& item : work) {
            const auto r = out.mRecords.find(item.mLabel);
            if (r != out.mRecords.end() && r->second.mDemotionReason != FitCadDemotionReason::None)
                items.emplace_back(item.mLabel, item.mFacets);
        }
        if (!items.empty()) {
            const std::size_t before = out.mNetworkLabels.size();
            const auto left = network(items);
            if (out.mNetworkLabels.size() > before) {
                // The regions the network took are its pieces now; what it
                // could not fit stays facets, with the adaptive face's reason.
                std::set<std::int64_t> kept;
                for (const auto& l : left)
                    kept.insert(l.first);
                for (const auto& it : items)
                    if (!kept.count(it.first))
                        out.mRecords.erase(it.first);
                for (const auto& l : left)
                    if (!out.mRecords.count(l.first)) {
                        FitCadFaceRecord rec;
                        rec.mPatch = l.first;
                        rec.mEscalationSteps = 1;
                        rec.mDemotionReason = FitCadDemotionReason::HausdorffExceeded;
                        out.mRecords[l.first] = rec;
                    }
            }
        }
    }
    // The re-segmentation rung: up to `mResegmentRounds` times, every region
    // still left as facets is cut by its normals (fcn_resegment, two pieces)
    // and the pieces climb the adaptive rungs again; one that fails is cut
    // again in the next round, never halved.
    for (std::int64_t round = 0; round < rOptions.mResegmentRounds; ++round) {
        std::map<std::int64_t, const std::vector<std::size_t>*> facets_by_label;
        for (const auto& it : facets_of)
            facets_by_label[it.first] = &it.second;
        for (const auto& p : pieces)
            facets_by_label[p.first] = &p.second;
        std::vector<std::tuple<std::int64_t, std::int64_t, std::vector<std::vector<std::size_t>>>>
            cuts;  // (label, escalation steps, pieces)
        for (const auto& [label, rec] : out.mRecords) {
            if (rec.mDemotionReason == FitCadDemotionReason::None)
                continue;
            const auto f = facets_by_label.find(label);
            if (f == facets_by_label.end())
                continue;
            auto cut = fcn_resegment(rFc, *f->second, 2);
            if (!cut.empty())
                cuts.emplace_back(label, rec.mEscalationSteps, std::move(cut));
        }
        if (cuts.empty())
            break;
        const std::size_t from = work.size();
        for (auto& [label, steps, cut] : cuts) {
            log::debug("fit_cad: re-segmentation round {} cuts patch {} into {} piece(s)",
                       round + 1, label, cut.size());
            out.mRecords.erase(label);
            ++out.mNumResegmented;
            for (std::vector<std::size_t>& facets : cut) {
                const std::int64_t piece = next_label++;
                work.push_back({piece, facets, rOptions.mMaxFreeformSplits, steps + 1,
                                rOptions.mResegmentRounds});
                pieces.emplace_back(piece, std::move(facets));
            }
        }
        run(from);
    }
    if (pieces.empty())
        return out;

    // The mesh to write: the pieces under their labels, each a fit-table row.
    out.mSplit = true;
    out.mMesh = detail::clone_mesh(
        rSource, [](DataLocation, const std::string&, std::string&) { return true; });
    std::vector<std::int64_t> new_label = cell_label;
    for (const auto& p : pieces)
        for (std::size_t f : p.second)
            new_label[static_cast<std::size_t>(rFc.mCell[f])] = p.first;
    std::vector<NDArray> blocks;
    for (std::size_t b = 0; b < rSource.NumCellBlocks(); ++b) {
        const std::size_t nc = rSource.Cells(b).NumCells();
        NDArray a = NDArray::Uninit(DType::Int64, {nc});
        for (std::size_t c = 0; c < nc; ++c)
            a.As<std::int64_t>()[c] = new_label[static_cast<std::size_t>(bases[b]) + c];
        blocks.push_back(std::move(a));
    }
    out.mMesh.AddCellData("segment:patch", std::move(blocks));
    const NDArray& params = rSource.FieldData("fit:patch_params");
    const NDArray& tols = rSource.FieldData("fit:patch_tolerance");
    const std::size_t m = ids.Size(), mp = m + pieces.size();
    const std::size_t np = m > 0 ? params.Size() / m : 10;
    NDArray nid = NDArray::Uninit(DType::Int64, {mp});
    NDArray ntype = NDArray::Uninit(DType::Int64, {mp});
    NDArray nparam = NDArray::Uninit(DType::Float64, {mp, np});
    NDArray ntol = NDArray::Uninit(DType::Float64, {mp});
    for (std::size_t r = 0; r < m; ++r) {
        nid.As<std::int64_t>()[r] = detail::read_int(ids, r);
        ntype.As<std::int64_t>()[r] = detail::read_int(types, r);
        for (std::size_t k = 0; k < np; ++k)
            nparam.As<double>()[r * np + k] = detail::read_double(params, r * np + k);
        ntol.As<double>()[r] = detail::read_double(tols, r);
    }
    for (std::size_t i = 0; i < pieces.size(); ++i) {
        nid.As<std::int64_t>()[m + i] = pieces[i].first;
        ntype.As<std::int64_t>()[m + i] = 0;
        for (std::size_t k = 0; k < np; ++k)
            nparam.As<double>()[(m + i) * np + k] = std::numeric_limits<double>::quiet_NaN();
        ntol.As<double>()[m + i] = Tol;
    }
    out.mMesh.AddFieldData("fit:patch_id", std::move(nid));
    out.mMesh.AddFieldData("fit:patch_type", std::move(ntype));
    out.mMesh.AddFieldData("fit:patch_params", std::move(nparam));
    out.mMesh.AddFieldData("fit:patch_tolerance", std::move(ntol));
    return out;
}

FitCadResult fcn_native_fit_cad(const Mesh& rMesh, const std::string& rStepPath,
                                const FitCadOptions& rOptions) {
    // A write entry point: the STEP header's provenance block holds the notes
    // this call raises (its own segmentation and fit included), never ones an
    // earlier, unrelated operation on this thread left behind.
    detail::provenance_begin_write();
    Mesh owned;
    FcnRefined refined;
    const Mesh& work = *fcn_prepare(rMesh, owned, rOptions, refined);
    const detail::ClosedFacets fc = detail::closed_facets(work, kFcnPrefix);
    if (fc.NumFacets() == 0)
        throw std::invalid_argument(std::string(kFcnPrefix) + "no surface facets");
    if (fc.mNumSkipped > 0)
        log::warn("fit_cad: {} cell(s) that are not surface facets (lines, vertices) are ignored",
                  fc.mNumSkipped);

    // The patches whose type is left out, or that a refinement found no
    // primitive for, are written as coplanar groups.
    Mesh typed;
    const Mesh* pSource = &work;
    const bool verdicts = refined.mRan && work.HasFieldData(kRefinePatchProxyName) &&
                          work.HasFieldData(kRefinePatchIdName) &&
                          work.HasFieldData("fit:patch_id");
    if (((rOptions.mTypes & kFitAllTypes) != kFitAllTypes || verdicts) &&
        work.HasFieldData("fit:patch_type")) {
        typed = detail::clone_mesh(
            work, [](DataLocation, const std::string&, std::string&) { return true; });
        const NDArray& types = work.FieldData("fit:patch_type");
        std::unordered_map<std::int64_t, std::int64_t> proxy_of;
        if (verdicts) {
            const NDArray& kinds = work.FieldData(kRefinePatchProxyName);
            const NDArray& labels = work.FieldData(kRefinePatchIdName);
            for (std::size_t i = 0; i < std::min(kinds.Size(), labels.Size()); ++i)
                proxy_of.emplace(detail::read_int(labels, i), detail::read_int(kinds, i));
        }
        NDArray kept = NDArray::Uninit(DType::Int64, {types.Size()});
        for (std::size_t r = 0; r < types.Size(); ++r) {
            const std::int64_t t = detail::read_int(types, r);
            bool drop = t >= 1 && t <= 5 && !(rOptions.mTypes & (1u << static_cast<unsigned>(t)));
            if (verdicts && r < work.FieldData("fit:patch_id").Size()) {
                const auto k = proxy_of.find(detail::read_int(work.FieldData("fit:patch_id"), r));
                drop = drop || (k != proxy_of.end() &&
                                (k->second == static_cast<std::int64_t>(RefineProxy::Facets) ||
                                 k->second == static_cast<std::int64_t>(RefineProxy::Freeform)));
            }
            kept.As<std::int64_t>()[r] = drop ? 0 : t;
        }
        typed.AddFieldData("fit:patch_type", std::move(kept));
        pSource = &typed;
    }
    // Outside the legacy mode, the patches the fit table leaves unexplained (its
    // type 0, or a refinement's freeform or facets verdict; not a type left
    // out of `mTypes`) are tried as adaptive freeform faces (v17.12.0).
    FcnFreeform freeform;
    if (rOptions.mFreeformMode != FitCadFreeformMode::Legacy && work.HasFieldData("fit:patch_id") &&
        work.HasFieldData("fit:patch_type")) {
        std::unordered_map<std::int64_t, std::int64_t> proxy_of;
        if (work.HasFieldData(kRefinePatchProxyName) && work.HasFieldData(kRefinePatchIdName)) {
            const NDArray& kinds = work.FieldData(kRefinePatchProxyName);
            const NDArray& labels = work.FieldData(kRefinePatchIdName);
            for (std::size_t i = 0; i < std::min(kinds.Size(), labels.Size()); ++i)
                proxy_of.emplace(detail::read_int(labels, i), detail::read_int(kinds, i));
        }
        const NDArray& ids = work.FieldData("fit:patch_id");
        const NDArray& types = work.FieldData("fit:patch_type");
        std::vector<std::int64_t> candidates;
        for (std::size_t r = 0; r < std::min(ids.Size(), types.Size()); ++r) {
            const std::int64_t label = detail::read_int(ids, r);
            const auto k = proxy_of.find(label);
            if (detail::read_int(types, r) == 0 ||
                (k != proxy_of.end() &&
                 (k->second == static_cast<std::int64_t>(RefineProxy::Facets) ||
                  k->second == static_cast<std::int64_t>(RefineProxy::Freeform))))
                candidates.push_back(label);
        }
        if (rOptions.mFreeform && !candidates.empty()) {
            const double wtol = fcn_tolerance(*pSource, fc, rOptions.mTolerance);
            freeform = fcn_freeform(*pSource, fc, candidates, rOptions, wtol);
            if (freeform.mSplit)
                pSource = &freeform.mMesh;
        }
    }
    StepWriteOptions opts;
    opts.mFaces = StepFaces::Fitted;
    opts.mTolerance = rOptions.mTolerance;
    detail::StepBrepResult brep;
    detail::write_step_report(rStepPath, *pSource, opts, &brep,
                              freeform.mFaces.empty() ? nullptr : &freeform.mFaces,
                              freeform.mEdgeCurves.empty() ? nullptr : &freeform.mEdgeCurves,
                              rOptions.mFreeformMode != FitCadFreeformMode::Legacy);

    FitCadResult result;
    result.mBackend = FitCadBackend::Native;
    // Outside the legacy mode, the solid as NURBS too (v17.12.0), for
    // `write_cad_json`; empty, with a warning, where a face cannot be one.
    if (rOptions.mFreeformMode != FitCadFreeformMode::Legacy) {
        result.mBrep = std::move(brep.mBrep);
        if (!brep.mBrepWhy.empty())
            log::warn("fit_cad: the native solid has no NURBS form ({}); the B-rep is empty",
                      brep.mBrepWhy);
    }
    result.mValid = brep.mValid;
    result.mTolerance = brep.mTolerance;
    result.mSkinVolume = brep.mVolume;
    result.mVolume = std::numeric_limits<double>::quiet_NaN();
    result.mVolumeRatio = std::numeric_limits<double>::quiet_NaN();
    result.mMaxDeviation = brep.mMaxDeviation;
    result.mNumFaces = brep.mNumFaces;
    result.mNumPlane = brep.mNumPlane;
    result.mNumSphere = brep.mNumSphere;
    result.mNumCylinder = brep.mNumCylinder;
    result.mNumCone = brep.mNumCone;
    result.mNumTorus = brep.mNumTorus;
    result.mNumBspline = brep.mNumBspline;
    result.mNumFreeformAdaptive = brep.mNumBspline;
    result.mNumNetworkEdges = freeform.mNumNetworkEdges;
    result.mNumExtraordinaryVertices = freeform.mNumExtraordinary;
    result.mMaxG1AngleDeg = freeform.mMaxG1AngleDeg;
    result.mNumSplitPatches = freeform.mNumSplit;
    result.mNumResegmented = freeform.mNumResegmented;
    result.mNumFacetedFaces = brep.mNumFacetedFaces;
    result.mNumPatches = brep.mNumPatches;
    result.mNumDemoted = brep.mNumDemoted;
    result.mNumExactEdges = brep.mNumLineEdges + brep.mNumCircleEdges;
    result.mNumApproxEdges = brep.mNumSplineEdges;
    result.mNumSolids = brep.mNumSolids;
    result.mClosed = brep.mValid;

    // The skin itself, each cell with the face it was written into.
    result.mMesh = detail::clone_mesh(
        work, [](DataLocation, const std::string&, std::string&) { return true; });
    const std::vector<std::int64_t> bases = detail::block_bases(work);
    const std::size_t ncells = static_cast<std::size_t>(detail::total_cells(bases));
    std::vector<std::int64_t> face(ncells, -1), type(ncells, 0), patch(ncells, -1);
    for (std::size_t f = fc.NumFacets(); f-- > 0;) {
        const std::size_t c = static_cast<std::size_t>(fc.mCell[f]);
        face[c] = brep.mFaceOfFacet[f];
        type[c] = brep.mTypeOfFacet[f];
        patch[c] = brep.mPatchOfFacet[f];
    }
    std::vector<NDArray> fb, tb, pb;
    for (std::size_t b = 0; b < work.NumCellBlocks(); ++b) {
        const std::size_t base = static_cast<std::size_t>(bases[b]);
        const std::size_t nc = work.Cells(b).NumCells();
        NDArray a = NDArray::Uninit(DType::Int64, {nc});
        NDArray t = NDArray::Uninit(DType::Int64, {nc});
        NDArray q = NDArray::Uninit(DType::Int64, {nc});
        for (std::size_t c = 0; c < nc; ++c) {
            a.As<std::int64_t>()[c] = face[base + c];
            t.As<std::int64_t>()[c] = type[base + c];
            q.As<std::int64_t>()[c] = patch[base + c];
        }
        fb.push_back(std::move(a));
        tb.push_back(std::move(t));
        pb.push_back(std::move(q));
    }
    result.mMesh.AddCellData(kFitCadFaceName, std::move(fb));
    result.mMesh.AddCellData(kFitCadTypeName, std::move(tb));
    result.mMesh.AddCellData(kFitCadPatchName, std::move(pb));

    // The freeform-mode report (v17.12.0): the facet area share over the skin
    // fanned into triangles, and a record per patch tried as a surface.
    result.mFreeformMode = rOptions.mFreeformMode;
    {
        std::vector<std::int64_t> tri, tri_type;
        for (std::size_t f = 0; f < fc.NumFacets(); ++f)
            for (std::int64_t k = fc.mStart[f] + 1; k + 1 < fc.mStart[f + 1]; ++k) {
                tri.insert(tri.end(), {fc.mNodes[static_cast<std::size_t>(fc.mStart[f])],
                                       fc.mNodes[static_cast<std::size_t>(k)],
                                       fc.mNodes[static_cast<std::size_t>(k + 1)]});
                tri_type.push_back(brep.mTypeOfFacet[f]);
            }
        result.mFacetAreaFraction = fcn_facet_area_fraction(fc.mXyz, tri, tri_type);
    }
    {
        std::unordered_map<std::int64_t, std::pair<std::int64_t, std::int64_t>> first_face;
        for (std::size_t f = 0; f < fc.NumFacets(); ++f)
            if (brep.mTypeOfFacet[f] != 0)
                first_face.emplace(brep.mPatchOfFacet[f],
                                   std::make_pair(brep.mFaceOfFacet[f], brep.mTypeOfFacet[f]));
        std::unordered_map<std::int64_t, std::size_t> num_facets;
        for (std::size_t f = 0; f < fc.NumFacets(); ++f)
            ++num_facets[brep.mPatchOfFacet[f]];
        const NDArray* pTypes = pSource->HasFieldData("fit:patch_type")
                                    ? &pSource->FieldData("fit:patch_type")
                                    : nullptr;
        double freeform_sq = 0.0;
        std::size_t freeform_n = 0;
        for (std::size_t r = 0; r < brep.mRowLabel.size(); ++r) {
            const std::int64_t lbl = brep.mRowLabel[r];
            const std::int64_t t = pTypes && r < pTypes->Size() ? detail::read_int(*pTypes, r) : 0;
            if (num_facets.find(lbl) == num_facets.end())
                continue;
            const auto ff = freeform.mRecords.find(lbl);
            if (ff != freeform.mRecords.end()) {
                // A patch the adaptive pass tried: its fit's record, and the
                // writer's verdict on it.
                FitCadFaceRecord rec = ff->second;
                const auto it = first_face.find(lbl);
                if (rec.mDemotionReason == FitCadDemotionReason::None && it != first_face.end() &&
                    brep.mRowReason[r] == 0) {
                    rec.mFace = it->second.first;
                    rec.mKind = it->second.second;
                } else if (rec.mDemotionReason == FitCadDemotionReason::None) {
                    rec.mDemotionReason =
                        brep.mRowReason[r] != 0
                            ? static_cast<FitCadDemotionReason>(brep.mRowReason[r])
                            : FitCadDemotionReason::Other;
                }
                result.mNumEscalations += rec.mEscalationSteps;
                if (rec.mKind == kFitCadBsplineType && freeform.mNetworkLabels.count(lbl)) {
                    ++result.mNumNetworkPatches;
                    --result.mNumFreeformAdaptive;
                }
                if (rec.mKind == kFitCadBsplineType) {
                    result.mMaxFreeformDistance =
                        std::max(result.mMaxFreeformDistance, rec.mMaxDist);
                    freeform_sq += rec.mRmsDist * rec.mRmsDist;
                    ++freeform_n;
                }
                result.mFaces.push_back(rec);
                continue;
            }
            if (t < 1 || t > 5)
                continue;
            FitCadFaceRecord rec;
            rec.mPatch = lbl;
            const auto it = first_face.find(lbl);
            if (it != first_face.end() && brep.mRowReason[r] == 0) {
                rec.mFace = it->second.first;
                rec.mKind = it->second.second;
            } else {
                rec.mDemotionReason = static_cast<FitCadDemotionReason>(
                    brep.mRowReason[r] != 0 ? brep.mRowReason[r] : 10);
            }
            result.mFaces.push_back(rec);
        }
        if (freeform_n > 0)
            result.mRmsFreeformDistance = std::sqrt(freeform_sq / static_cast<double>(freeform_n));
    }
    if (rOptions.mFreeformMode != FitCadFreeformMode::Legacy)
        result.mMesh.AddFieldData(kFitCadFaceTableName, fcn_face_table(result.mFaces));
    if (result.mNumDemoted > 0)
        log::info("fit_cad: {} fitted patch(es) written as coplanar facets", result.mNumDemoted);
    fcn_report_refined(refined, work, result);
    return result;
}

}  // namespace

namespace detail {

// The native backend's helpers the OpenCASCADE backend builds on, exported for
// providers/occt (see detail/fit_cad_native.hpp).

double native_facet_area_fraction(const std::vector<double>& rXyz,
                                  const std::vector<std::int64_t>& rTri,
                                  const std::vector<std::int64_t>& rType) {
    return fcn_facet_area_fraction(rXyz, rTri, rType);
}

NDArray native_face_table(const std::vector<FitCadFaceRecord>& rFaces) {
    return fcn_face_table(rFaces);
}

const Mesh* native_prepare(const Mesh& rMesh, Mesh& rOwned, const FitCadOptions& rOptions,
                           FcnRefined& rRefined) {
    return fcn_prepare(rMesh, rOwned, rOptions, rRefined);
}

void native_report_refined(const FcnRefined& rRefined, const Mesh& rWork,
                           FitCadResult& rResult) {
    fcn_report_refined(rRefined, rWork, rResult);
}

std::vector<FcnRegionNetwork> native_networks(const ClosedFacets& rFc,
                                              const std::vector<std::vector<std::size_t>>& rItems,
                                              const FitCadOptions& rOptions, double Tol) {
    return fcn_networks(rFc, rItems, rOptions, Tol);
}

std::vector<std::vector<std::size_t>> native_resegment(const ClosedFacets& rFc,
                                                       const std::vector<std::size_t>& rFacets,
                                                       std::size_t K, bool OnlyAcrossCreases) {
    return fcn_resegment(rFc, rFacets, K, OnlyAcrossCreases);
}

}  // namespace detail


// One build of the solid, as `fit_cad` makes it: the options checked, the
// backend chosen, the freeform modes' general fit run first. `fit_cad` adds the
// constrained-primitives choice on top of it.
static FitCadResult fit_cad_build(const Mesh& rMesh, const std::string& rStepPath,
                                  const FitCadOptions& rOptions) {
    if (std::isnan(rOptions.mTolerance) || std::isinf(rOptions.mTolerance))
        throw std::invalid_argument(std::string(kFcnPrefix) + "tolerance must be finite");
    if ((rOptions.mTypes & ~kFitAllTypes) != 0)
        throw std::invalid_argument(std::string(kFcnPrefix) +
                                    "types may name only plane, sphere, cylinder, cone and torus");
    if (rOptions.mBackend != FitCadBackend::Auto && rOptions.mBackend != FitCadBackend::Occt &&
        rOptions.mBackend != FitCadBackend::Native)
        throw std::invalid_argument(std::string(kFcnPrefix) +
                                    "backend must be 'auto', 'occt' or 'native'");
    // The general freeform options (v17.11.0), checked whether or not they are used.
    if (!(rOptions.mFreeformTolerance >= 0.0) || std::isinf(rOptions.mFreeformTolerance))
        throw std::invalid_argument(std::string(kFcnPrefix) +
                                    "freeform_tolerance must be finite and >= 0");
    if (!(rOptions.mFreeformSmoothing >= 0.0) || std::isinf(rOptions.mFreeformSmoothing))
        throw std::invalid_argument(std::string(kFcnPrefix) +
                                    "freeform_smoothing must be finite and >= 0");
    if (rOptions.mMaxFreeformPoles < 4 || rOptions.mMaxFreeformPoles > 1024)
        throw std::invalid_argument(std::string(kFcnPrefix) +
                                    "max_freeform_poles must be 4 to 1024");
    if (rOptions.mMaxFreeformIterations < 1)
        throw std::invalid_argument(std::string(kFcnPrefix) +
                                    "max_freeform_iterations must be >= 1");
    if (rOptions.mMaxFreeformSplits < 0)
        throw std::invalid_argument(std::string(kFcnPrefix) + "max_freeform_splits must be >= 0");
    if (!(rOptions.mSewTolerance >= 0.0) || std::isinf(rOptions.mSewTolerance))
        throw std::invalid_argument(std::string(kFcnPrefix) +
                                    "sew_tolerance must be finite and >= 0");
    // The freeform modes (v17.12.0), checked whether or not they are used.
    if (static_cast<std::uint32_t>(rOptions.mFreeformMode) > 3u)
        throw std::invalid_argument(
            std::string(kFcnPrefix) +
            "freeform_mode must be 'legacy', 'adaptive', 'network' or 'auto'");
    if (static_cast<std::uint32_t>(rOptions.mFreeformParameterization) > 3u)
        throw std::invalid_argument(
            std::string(kFcnPrefix) +
            "freeform_parameterization must be 'auto', 'projection', 'lscm' or 'harmonic'");
    if (static_cast<std::uint32_t>(rOptions.mExtraordinaryVertex) > 3u)
        throw std::invalid_argument(std::string(kFcnPrefix) +
                                    "extraordinary_vertex must be 'g0', 'g1', 'g1x' or 'g1xx'");
    if (static_cast<std::uint32_t>(rOptions.mLocalRefine) > 1u)
        throw std::invalid_argument(std::string(kFcnPrefix) +
                                    "local_refine must be 'knots' or 'thb'");
    if (!(rOptions.mFreeformG1ToleranceDeg > 0.0) || std::isinf(rOptions.mFreeformG1ToleranceDeg))
        throw std::invalid_argument(std::string(kFcnPrefix) +
                                    "freeform_g1_tolerance_deg must be finite and > 0");
    if (rOptions.mFreeformMaxControlPoints < 16 || rOptions.mFreeformMaxControlPoints > 1048576)
        throw std::invalid_argument(std::string(kFcnPrefix) +
                                    "freeform_max_control_points must be 16 to 1048576");
    if (!(rOptions.mLayoutTargetPatchSize > 0.0) || std::isinf(rOptions.mLayoutTargetPatchSize))
        throw std::invalid_argument(std::string(kFcnPrefix) +
                                    "layout_target_patch_size must be finite and > 0");
    if (rOptions.mLayoutMaxValenceDeviation < 1 || rOptions.mLayoutMaxValenceDeviation > 3)
        throw std::invalid_argument(std::string(kFcnPrefix) +
                                    "layout_max_valence_deviation must be 1 to 3");
    if (rOptions.mResegmentRounds < 0)
        throw std::invalid_argument(std::string(kFcnPrefix) + "resegment_rounds must be >= 0");
    // The freeform modes build on the general fit, and their faces run on
    // past their patches so an edge with an analytic neighbour is still an
    // intersection.
    if (rOptions.mFreeformMode != FitCadFreeformMode::Legacy &&
        (!rOptions.mFreeformGeneral || !rOptions.mExtendSurfaces)) {
        FitCadOptions general = rOptions;
        general.mFreeformGeneral = true;
        general.mExtendSurfaces = true;
        return fit_cad_build(rMesh, rStepPath, general);
    }
    const auto p_kernel = cad_kernel();
    if (rOptions.mBackend != FitCadBackend::Native && p_kernel && p_kernel->mFitCad)
        return p_kernel->mFitCad(rMesh, rStepPath, rOptions);
    if (rOptions.mBackend == FitCadBackend::Occt)
        throw Unsupported(
            "fit_cad: no OpenCASCADE backend is available (the release binaries carry it and "
            "the Python wheels provide it through the meshioplusplus-occt package, pip install "
            "meshioplusplus[occt]); a source build needs -DMESHIOPLUSPLUS_WITH_OCCT=ON "
            "(OpenCASCADE 7.6 or later, build/build-occt.sh builds one), or use the native "
            "backend");
    // OpenCASCADE is what the result is meant to be built with: say once that
    // this build falls back, and what that costs. An atomic flag rather than
    // std::call_once, whose libstdc++ thread-locals the xeus-cpp JIT running
    // the C++ notebooks on the single header cannot resolve.
    if (rOptions.mBackend == FitCadBackend::Auto) {
        static std::atomic<bool> s_fcn_fallback_warned{false};
        if (!s_fcn_fallback_warned.exchange(true))
            log::warn(
                "fit_cad: no OpenCASCADE is available, so the kernel-free native backend "
                "is used: no freeform, fillet or thread faces, no B-rep check and no solid "
                "volume. The release binaries carry OpenCASCADE and the Python wheels get it "
                "from the meshioplusplus-occt package; a source build finds one with "
                "-DMESHIOPLUSPLUS_WITH_OCCT=ON (the default)");
    }
    return fcn_native_fit_cad(rMesh, rStepPath, rOptions);
}

FitCadResult fit_cad(const Mesh& rMesh, const std::string& rStepPath,
                     const FitCadOptions& rOptions) {
    if (!rOptions.mConstrainPrimitives)
        return fit_cad_build(rMesh, rStepPath, rOptions);
    // The constrained surfaces are fitted within tolerance, but the solid built
    // from them can lose a junction the unconstrained surfaces kept (a fillet
    // tangent to its wall, say), so the constrained solid is kept only where it
    // demotes no more patches than the unconstrained one. The constrained solid
    // is written to `rStepPath`; the unconstrained one is built unwritten, and
    // written only when it is the one kept.
    FitCadResult on = fit_cad_build(rMesh, rStepPath, rOptions);
    FitCadOptions plain = rOptions;
    plain.mConstrainPrimitives = false;
    const FitCadResult off = fit_cad_build(rMesh, "", plain);
    const bool keep_on = !(off.mValid && (!on.mValid || off.mNumDemoted < on.mNumDemoted));
    if (keep_on)
        return on;
    if (!on.mValid)
        log::warn("fit_cad: the constrained primitives give a solid that is not valid; the "
                  "unconstrained solid is returned");
    else
        log::warn("fit_cad: the constrained primitives demote {} patch(es) against {} for the "
                  "unconstrained solid; the unconstrained solid is returned",
                  on.mNumDemoted, off.mNumDemoted);
    FitCadResult out = fit_cad_build(rMesh, rStepPath, plain);
    out.mNumConstrainedPrimitives = 0;
    out.mNumConstraintGroups = 0;
    return out;
}

}  // namespace meshioplusplus

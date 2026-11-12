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

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

// Project includes
#include "meshioplusplus/detail/surface_distance.hpp"
#include "meshioplusplus/operations/sdf.hpp"
#include "brep_deviation.hpp"
#include "cdt2d.hpp"
#include "geometry_utils.hpp"
#include "hausdorff_samples.hpp"
#include "iga_kernels.hpp"
#include "math_utils.hpp"
#include "nurbs_eval.hpp"

namespace meshioplusplus::detail {

namespace {

constexpr double kBdPi = 3.14159265358979323846;

// A face prepared for evaluation: its surface as an evaluation net, its knot
// domain and its trims as parameter-space polygons.
struct BdFace {
    bool mOk = false;
    IgaEvalNet mNet;
    std::size_t mNu = 0, mNv = 0;
    int mP = 0, mQ = 0;
    double mU0 = 0.0, mU1 = 0.0, mV0 = 0.0, mV1 = 0.0;
    std::vector<std::vector<double>> mLu, mLv;
    /// The trims on the surface: per loop, closed polylines in 3D (each
    /// parameter-space segment subdivided) and their parameters.
    std::vector<std::vector<Vec3>> mBoundary;
    std::vector<std::vector<std::array<double, 2>>> mBoundaryUv;
    /// Its grid nodes inside the trims and their parameters: the starts of
    /// an exhaustive projection.
    std::vector<Vec3> mNodes;
    std::vector<std::array<double, 2>> mNodeUv;
};

// A trimming curve as polyline points in its loop's direction, appended to
// (rU, rV): the control points of a polyline, else samples per knot span.
void bd_trim_polyline(const BrepTrim& rTrim, std::size_t PerSpan, std::vector<double>& rU,
                      std::vector<double>& rV) {
    const BrepCurve& c = rTrim.mCurve;
    if (c.mControlPoints.Shape().size() != 2 || c.mControlPoints.Shape()[1] < 2)
        return;
    const std::size_t n = c.mControlPoints.Shape()[0];
    const std::size_t dim = c.mControlPoints.Shape()[1];
    const double* cp = c.mControlPoints.As<double>();
    std::vector<double> pu, pv;
    const bool weighted = c.mRational && c.mWeights.Size() == n;
    if (c.mDegree <= 1 && !weighted) {
        for (std::size_t i = 0; i < n; ++i) {
            pu.push_back(cp[i * dim]);
            pv.push_back(cp[i * dim + 1]);
        }
    } else if (n > static_cast<std::size_t>(c.mDegree) &&
               c.mKnots.size() == n + static_cast<std::size_t>(c.mDegree) + 1) {
        const double* w = weighted ? c.mWeights.As<double>() : nullptr;
        double t0 = c.mKnots[static_cast<std::size_t>(c.mDegree)], t1 = c.mKnots[n];
        if (c.mActiveRange.size() == 2) {
            t0 = std::max(t0, c.mActiveRange[0]);
            t1 = std::min(t1, c.mActiveRange[1]);
        }
        const std::size_t spans = std::max<std::size_t>(1, n - static_cast<std::size_t>(c.mDegree));
        const std::size_t m = std::max<std::size_t>(2, spans * PerSpan + 1);
        std::vector<double> out(dim);
        for (std::size_t k = 0; k < m; ++k) {
            const double t = t0 + (t1 - t0) * static_cast<double>(k) / static_cast<double>(m - 1);
            const std::size_t s = nurbs_find_span(c.mDegree, c.mKnots, n, t);
            nurbs_curve_point_in_span(c.mDegree, c.mKnots, cp, w, dim, s, t, out.data());
            pu.push_back(out[0]);
            pv.push_back(out[1]);
        }
    }
    if (!rTrim.mDirection) {
        std::reverse(pu.begin(), pu.end());
        std::reverse(pv.begin(), pv.end());
    }
    rU.insert(rU.end(), pu.begin(), pu.end());
    rV.insert(rV.end(), pv.begin(), pv.end());
}

BdFace bd_face(const BrepInfo& rBrep, const BrepFace& rFace, const BrepDeviationOptions& rOptions) {
    BdFace f;
    if (rFace.mPatch < 0 || static_cast<std::size_t>(rFace.mPatch) >= rBrep.mPatches.size())
        return f;
    const NurbsPatch& patch = rBrep.mPatches[static_cast<std::size_t>(rFace.mPatch)];
    if (patch.mDegrees.size() != 2 || patch.mKnots.size() != 2 ||
        patch.mControlPoints.Shape().size() != 2 || patch.mControlPoints.Shape()[1] != 3)
        return f;
    f.mP = patch.mDegrees[0];
    f.mQ = patch.mDegrees[1];
    if (f.mP < 1 || f.mQ < 1 || patch.mKnots[0].size() <= 2 * static_cast<std::size_t>(f.mP) + 1 ||
        patch.mKnots[1].size() <= 2 * static_cast<std::size_t>(f.mQ) + 1)
        return f;
    f.mNu = patch.mKnots[0].size() - static_cast<std::size_t>(f.mP) - 1;
    f.mNv = patch.mKnots[1].size() - static_cast<std::size_t>(f.mQ) - 1;
    if (patch.mControlPoints.Shape()[0] != f.mNu * f.mNv)
        return f;
    f.mNet.mDim = 2;
    f.mNet.mDegrees = {f.mP, f.mQ, 0};
    f.mNet.mKnots = {&patch.mKnots[0], &patch.mKnots[1], nullptr};
    f.mNet.mCounts = {f.mNu, f.mNv, 1};
    f.mNet.mValues = patch.mControlPoints.As<double>();
    f.mNet.mColumns = 3;
    f.mNet.mWeights =
        patch.mWeights.Size() == f.mNu * f.mNv ? patch.mWeights.As<double>() : nullptr;
    f.mU0 = patch.mKnots[0][static_cast<std::size_t>(f.mP)];
    f.mU1 = patch.mKnots[0][f.mNu];
    f.mV0 = patch.mKnots[1][static_cast<std::size_t>(f.mQ)];
    f.mV1 = patch.mKnots[1][f.mNv];
    if (!(f.mU1 > f.mU0) || !(f.mV1 > f.mV0))
        return f;
    if (rFace.mTrimmed)
        for (const BrepLoop& loop : rFace.mLoops) {
            std::vector<double> lu, lv;
            for (const BrepTrim& t : loop.mTrims)
                bd_trim_polyline(t, rOptions.mTrimSamplesPerSpan, lu, lv);
            if (lu.size() >= 3) {
                f.mLu.push_back(std::move(lu));
                f.mLv.push_back(std::move(lv));
            }
        }
    f.mOk = true;
    return f;
}

bool bd_inside(const BdFace& rF, double U, double V) {
    if (rF.mLu.empty())
        return true;
    bool in = false;
    for (std::size_t i = 0; i < rF.mLu.size(); ++i)
        in = in != GeometryUtils::PointInPolygon2D(rF.mLu[i], rF.mLv[i], U, V);
    return in;
}

// The surface point and its first derivatives at (U, V), clamped to the domain.
void bd_eval(const BdFace& rF, double U, double V, Vec3& rS, Vec3& rSu, Vec3& rSv) {
    U = std::clamp(U, rF.mU0, rF.mU1);
    V = std::clamp(V, rF.mV0, rF.mV1);
    const std::size_t span[2] = {nurbs_find_span(rF.mP, *rF.mNet.mKnots[0], rF.mNu, U),
                                 nurbs_find_span(rF.mQ, *rF.mNet.mKnots[1], rF.mNv, V)};
    const double uv[2] = {U, V};
    double out[3], der[6];
    iga_eval_point(rF.mNet, span, uv, 1, out, der);
    rS = {out[0], out[1], out[2]};
    rSu = {der[0], der[1], der[2]};
    rSv = {der[3], der[4], der[5]};
}

// The trims of a face as closed 3D polylines on its surface: every
// parameter-space segment cut in two while its midpoint on the surface lies
// further than `Eps` from its chord (up to `kBdBoundaryDepth` times), so a
// trim that is straight in (u, v) follows the surface between its points.
constexpr int kBdBoundaryDepth = 10;

void bd_boundary_segment(const BdFace& rF, const std::array<double, 2>& rA, const Vec3& rPa,
                         const std::array<double, 2>& rB, const Vec3& rPb, double Eps, int Depth,
                         std::vector<Vec3>& rPts, std::vector<std::array<double, 2>>& rUv) {
    const std::array<double, 2> m{0.5 * (rA[0] + rB[0]), 0.5 * (rA[1] + rB[1])};
    Vec3 pm, su, sv;
    bd_eval(rF, m[0], m[1], pm, su, sv);
    const Vec3 chord{0.5 * (rPa[0] + rPb[0]), 0.5 * (rPa[1] + rPb[1]), 0.5 * (rPa[2] + rPb[2])};
    if (Depth < kBdBoundaryDepth && MathUtils::Norm(vec3_sub(pm, chord)) > Eps) {
        bd_boundary_segment(rF, rA, rPa, m, pm, Eps, Depth + 1, rPts, rUv);
        bd_boundary_segment(rF, m, pm, rB, rPb, Eps, Depth + 1, rPts, rUv);
        return;
    }
    // The segment's start; its end starts the next one.
    rPts.push_back(rPa);
    rUv.push_back(rA);
}

void bd_boundary(BdFace& rF, double Eps) {
    for (std::size_t l = 0; l < rF.mLu.size(); ++l) {
        std::vector<Vec3> pts;
        std::vector<std::array<double, 2>> uvs;
        const std::size_t n = rF.mLu[l].size();
        std::vector<Vec3> at(n);
        for (std::size_t i = 0; i < n; ++i) {
            Vec3 su, sv;
            bd_eval(rF, rF.mLu[l][i], rF.mLv[l][i], at[i], su, sv);
        }
        for (std::size_t i = 0; i < n; ++i) {
            const std::size_t j = (i + 1) % n;
            bd_boundary_segment(rF, {rF.mLu[l][i], rF.mLv[l][i]}, at[i],
                                {rF.mLu[l][j], rF.mLv[l][j]}, at[j], Eps, 0, pts, uvs);
        }
        rF.mBoundary.push_back(std::move(pts));
        rF.mBoundaryUv.push_back(std::move(uvs));
    }
}

// The distance from rP to a face's boundary polylines, with the parameters of
// the nearest polyline point (interpolated along its segment).
double bd_boundary_distance(const BdFace& rF, const Vec3& rP, double& rU, double& rV) {
    double best = std::numeric_limits<double>::infinity();
    for (std::size_t l = 0; l < rF.mBoundary.size(); ++l) {
        const std::vector<Vec3>& pts = rF.mBoundary[l];
        const std::size_t n = pts.size();
        for (std::size_t i = 0; i < n; ++i) {
            const Vec3& a = pts[i];
            const Vec3& b = pts[(i + 1) % n];
            const Vec3 ab = vec3_sub(b, a);
            const double len2 = MathUtils::Dot(ab, ab);
            double t = len2 > 0.0 ? MathUtils::Dot(vec3_sub(rP, a), ab) / len2 : 0.0;
            t = std::clamp(t, 0.0, 1.0);
            const Vec3 q{a[0] + t * ab[0], a[1] + t * ab[1], a[2] + t * ab[2]};
            const double d = MathUtils::Norm(vec3_sub(rP, q));
            if (d < best) {
                best = d;
                const auto& ua = rF.mBoundaryUv[l][i];
                const auto& ub = rF.mBoundaryUv[l][(i + 1) % n];
                rU = ua[0] + t * (ub[0] - ua[0]);
                rV = ua[1] + t * (ub[1] - ua[1]);
            }
        }
    }
    return best;
}

// The trimmed domain of a face triangulated in its parameters, appended to
// the soup with each triangle's corner parameters and face; false (nothing
// appended) for an untrimmed face or trims the triangulation refuses.
bool bd_triangulate(const BdFace& rF, std::size_t Face, const std::vector<double>& rNu,
                    const std::vector<double>& rNv, const std::vector<char>& rInside,
                    TriangleSoup& rSoup, std::vector<std::array<double, 6>>& rUv,
                    std::vector<std::size_t>& rFace) {
    if (rF.mLu.empty())
        return false;
    double lo[2] = {rF.mLu[0][0], rF.mLv[0][0]}, hi[2] = {lo[0], lo[1]};
    for (std::size_t l = 0; l < rF.mLu.size(); ++l)
        for (std::size_t k = 0; k < rF.mLu[l].size(); ++k) {
            if (!cdt2d_in_range(rF.mLu[l][k]) || !cdt2d_in_range(rF.mLv[l][k]))
                return false;
            lo[0] = std::min(lo[0], rF.mLu[l][k]);
            hi[0] = std::max(hi[0], rF.mLu[l][k]);
            lo[1] = std::min(lo[1], rF.mLv[l][k]);
            hi[1] = std::max(hi[1], rF.mLv[l][k]);
        }
    if (!(hi[0] > lo[0]) || !(hi[1] > lo[1]))
        return false;
    Cdt2d cdt(lo[0], lo[1], hi[0], hi[1]);
    for (std::size_t l = 0; l < rF.mLu.size(); ++l) {
        std::vector<Cdt2d::Vec2> pts;
        for (std::size_t k = 0; k < rF.mLu[l].size(); ++k)
            pts.push_back({rF.mLu[l][k], rF.mLv[l][k]});
        const std::vector<std::int64_t> ids = cdt.InsertBatch(pts);
        for (std::size_t k = 0; k < ids.size(); ++k) {
            const std::int64_t a = ids[k], b = ids[(k + 1) % ids.size()];
            if (a != b && cdt.InsertSegment(a, b) != Cdt2dStatus::Ok)
                return false;
        }
    }
    if (cdt.Classify() != Cdt2dStatus::Ok)
        return false;
    for (std::size_t k = 0; k < rInside.size(); ++k)
        if (rInside[k] && cdt2d_in_range(rNu[k]) && cdt2d_in_range(rNv[k]))
            cdt.InsertSteiner({rNu[k], rNv[k]});
    const std::vector<Cdt2d::Tri> tris = cdt.InsideTriangles();
    if (tris.empty())
        return false;
    std::vector<Vec3> at(static_cast<std::size_t>(cdt.NumVertices()));
    std::vector<char> done(at.size(), 0);
    const auto point = [&](std::int64_t V) {
        const std::size_t i = static_cast<std::size_t>(V);
        if (!done[i]) {
            Vec3 su, sv;
            bd_eval(rF, cdt.Point(V)[0], cdt.Point(V)[1], at[i], su, sv);
            done[i] = 1;
        }
        return at[i];
    };
    for (const Cdt2d::Tri& t : tris) {
        const std::int64_t base = static_cast<std::int64_t>(rSoup.mPoints.size());
        for (std::int64_t v : t) {
            const Vec3 p = point(v);
            rSoup.mPoints.push_back(p);
            rSoup.mCorners.push_back(p);
        }
        rSoup.mVertices.push_back({base, base + 1, base + 2});
        rSoup.mSourceCell.push_back(static_cast<std::int64_t>(rFace.size()));
        rUv.push_back({cdt.Point(t[0])[0], cdt.Point(t[0])[1], cdt.Point(t[1])[0],
                       cdt.Point(t[1])[1], cdt.Point(t[2])[0], cdt.Point(t[2])[1]});
        rFace.push_back(Face);
    }
    return true;
}

// The unsigned angle, in degrees within [0, 90], between two directions; 0
// when either vanishes.
double bd_angle_deg(const Vec3& rA, const Vec3& rB) {
    const double la = MathUtils::Norm(rA), lb = MathUtils::Norm(rB);
    if (!(la > 0.0) || !(lb > 0.0))
        return 0.0;
    const double c = std::min(1.0, std::abs(MathUtils::Dot(rA, rB)) / (la * lb));
    return std::acos(c) * 180.0 / kBdPi;
}

// Gauss-Newton foot point of rP on the face from (rU, rV); returns the
// distance and leaves the parameters and the surface normal at the foot.
double bd_project(const BdFace& rF, const Vec3& rP, std::size_t MaxIterations, double& rU,
                  double& rV, Vec3& rNormal) {
    Vec3 s, su, sv;
    bd_eval(rF, rU, rV, s, su, sv);
    double best = MathUtils::Norm(vec3_sub(rP, s));
    for (std::size_t it = 0; it < MaxIterations; ++it) {
        const Vec3 r = vec3_sub(rP, s);
        const double a = MathUtils::Dot(su, su), b = MathUtils::Dot(su, sv),
                     c = MathUtils::Dot(sv, sv);
        const double g0 = MathUtils::Dot(su, r), g1 = MathUtils::Dot(sv, r);
        const double det = a * c - b * b;
        if (!(det > 1e-300))
            break;
        const double du = (c * g0 - b * g1) / det, dv = (a * g1 - b * g0) / det;
        const double nu = std::clamp(rU + du, rF.mU0, rF.mU1);
        const double nv = std::clamp(rV + dv, rF.mV0, rF.mV1);
        Vec3 s2, su2, sv2;
        bd_eval(rF, nu, nv, s2, su2, sv2);
        const double d2 = MathUtils::Norm(vec3_sub(rP, s2));
        if (!(d2 < best))
            break;
        const bool small = std::abs(nu - rU) <= 1e-14 * (rF.mU1 - rF.mU0) &&
                           std::abs(nv - rV) <= 1e-14 * (rF.mV1 - rF.mV0);
        rU = nu;
        rV = nv;
        s = s2;
        su = su2;
        sv = sv2;
        best = d2;
        if (small)
            break;
    }
    rNormal = vec3_cross(su, sv);
    return best;
}

}  // namespace

BrepDeviation brep_deviation(const BrepInfo& rBrep, const std::vector<double>& rXyz,
                             const std::vector<std::int64_t>& rTri,
                             const BrepDeviationOptions& rOptions,
                             const std::vector<std::int64_t>* pTriFace) {
    BrepDeviation out;
    out.mFaces.assign(rBrep.mFaces.size(), BrepFaceDeviation{});
    if (rBrep.mFaces.empty()) {
        out.mWhy = "the B-rep has no face";
        return out;
    }
    if (rTri.size() < 3) {
        out.mWhy = "the skin has no triangle";
        return out;
    }
    const std::size_t N = std::max<std::size_t>(2, rOptions.mSamplesPerDirection);

    // The skin, as a soup to measure the faces' samples against.
    TriangleSoup skin;
    skin.mPoints.reserve(rXyz.size() / 3);
    for (std::size_t i = 0; i + 2 < rXyz.size(); i += 3)
        skin.mPoints.push_back({rXyz[i], rXyz[i + 1], rXyz[i + 2]});
    const auto pt = [&](std::int64_t V) { return skin.mPoints[static_cast<std::size_t>(V)]; };
    for (std::size_t t = 0; t + 2 < rTri.size(); t += 3) {
        skin.mVertices.push_back({rTri[t], rTri[t + 1], rTri[t + 2]});
        skin.mCorners.push_back(pt(rTri[t]));
        skin.mCorners.push_back(pt(rTri[t + 1]));
        skin.mCorners.push_back(pt(rTri[t + 2]));
        skin.mSourceCell.push_back(static_cast<std::int64_t>(t / 3));
    }
    SurfaceDistanceOptions dopt;
    dopt.mSign = SdfSign::Unsigned;
    dopt.mWatertightCheck = SdfWatertightCheck::Off;
    const DistanceQuery skin_query = build_distance_query(skin, dopt);

    // Every face on its grid: the nodes, which are inside the trims, and the
    // cells (two triangles each) a corner of which is.
    TriangleSoup cad;
    std::vector<std::array<double, 6>> cad_uv;  // per CAD triangle: its corners' (u, v)
    std::vector<std::size_t> cad_face;
    std::vector<Vec3> samples;
    std::vector<std::size_t> sample_face;
    std::vector<BdFace> faces;
    faces.reserve(rBrep.mFaces.size());
    // The trim curves follow their surfaces to a millionth of the skin's size.
    double boundary_eps = 0.0;
    {
        double lo[3] = {1e300, 1e300, 1e300}, hi[3] = {-1e300, -1e300, -1e300};
        for (const Vec3& p : skin.mPoints)
            for (int d = 0; d < 3; ++d) {
                lo[d] = std::min(lo[d], p[static_cast<std::size_t>(d)]);
                hi[d] = std::max(hi[d], p[static_cast<std::size_t>(d)]);
            }
        boundary_eps =
            1e-6 * std::sqrt((hi[0] - lo[0]) * (hi[0] - lo[0]) + (hi[1] - lo[1]) * (hi[1] - lo[1]) +
                             (hi[2] - lo[2]) * (hi[2] - lo[2]));
    }
    for (std::size_t fi = 0; fi < rBrep.mFaces.size(); ++fi) {
        faces.push_back(bd_face(rBrep, rBrep.mFaces[fi], rOptions));
        if (faces.back().mOk)
            bd_boundary(faces.back(), boundary_eps);
        const BdFace& f = faces.back();
        if (!f.mOk)
            continue;
        double lo_u = f.mU0, hi_u = f.mU1, lo_v = f.mV0, hi_v = f.mV1;
        if (!f.mLu.empty()) {
            lo_u = hi_u = f.mLu[0][0];
            lo_v = hi_v = f.mLv[0][0];
            for (std::size_t l = 0; l < f.mLu.size(); ++l)
                for (std::size_t k = 0; k < f.mLu[l].size(); ++k) {
                    lo_u = std::min(lo_u, f.mLu[l][k]);
                    hi_u = std::max(hi_u, f.mLu[l][k]);
                    lo_v = std::min(lo_v, f.mLv[l][k]);
                    hi_v = std::max(hi_v, f.mLv[l][k]);
                }
            lo_u = std::max(lo_u, f.mU0);
            hi_u = std::min(hi_u, f.mU1);
            lo_v = std::max(lo_v, f.mV0);
            hi_v = std::min(hi_v, f.mV1);
            if (!(hi_u > lo_u) || !(hi_v > lo_v))
                continue;
        }
        // A face of more knot spans than `N` (a thread's many turns) gets a
        // node per span, so no node is a turn away from its neighbours.
        const std::size_t Gu =
            std::min<std::size_t>(1024, std::max(N, f.mNu > static_cast<std::size_t>(f.mP)
                                                        ? f.mNu - static_cast<std::size_t>(f.mP)
                                                        : std::size_t{0}));
        const std::size_t Gv =
            std::min<std::size_t>(1024, std::max(N, f.mNv > static_cast<std::size_t>(f.mQ)
                                                        ? f.mNv - static_cast<std::size_t>(f.mQ)
                                                        : std::size_t{0}));
        std::vector<Vec3> node((Gu + 1) * (Gv + 1));
        std::vector<char> inside((Gu + 1) * (Gv + 1), 0);
        std::vector<double> nu((Gu + 1) * (Gv + 1)), nv((Gu + 1) * (Gv + 1));
        for (std::size_t j = 0; j <= Gv; ++j)
            for (std::size_t i = 0; i <= Gu; ++i) {
                const std::size_t k = j * (Gu + 1) + i;
                nu[k] = lo_u + (hi_u - lo_u) * static_cast<double>(i) / static_cast<double>(Gu);
                nv[k] = lo_v + (hi_v - lo_v) * static_cast<double>(j) / static_cast<double>(Gv);
                Vec3 su, sv;
                bd_eval(f, nu[k], nv[k], node[k], su, sv);
                inside[k] = bd_inside(f, nu[k], nv[k]) ? 1 : 0;
                if (inside[k]) {
                    samples.push_back(node[k]);
                    sample_face.push_back(fi);
                    faces.back().mNodes.push_back(node[k]);
                    faces.back().mNodeUv.push_back({nu[k], nv[k]});
                }
            }
        // The face's triangles for the mesh-to-CAD side: its trimmed domain
        // triangulated in its parameters (the trims as constraints, the
        // inside grid nodes as interior points), so they cover the face and
        // nothing past its trims; failing that (trims that cross), the grid
        // cells a node of which is inside or a trim runs through.
        if (bd_triangulate(f, fi, nu, nv, inside, cad, cad_uv, cad_face))
            continue;
        std::vector<char> crossed(Gu * Gv, 0);
        for (std::size_t l = 0; l < f.mBoundaryUv.size(); ++l)
            for (const auto& uv : f.mBoundaryUv[l]) {
                const double fu = (uv[0] - lo_u) / (hi_u - lo_u) * static_cast<double>(Gu);
                const double fv = (uv[1] - lo_v) / (hi_v - lo_v) * static_cast<double>(Gv);
                if (!(fu >= 0.0) || !(fv >= 0.0) || fu > static_cast<double>(Gu) ||
                    fv > static_cast<double>(Gv))
                    continue;
                const std::size_t ci = std::min(Gu - 1, static_cast<std::size_t>(fu));
                const std::size_t cj = std::min(Gv - 1, static_cast<std::size_t>(fv));
                crossed[cj * Gu + ci] = 1;
            }
        for (std::size_t j = 0; j < Gv; ++j)
            for (std::size_t i = 0; i < Gu; ++i) {
                const std::size_t a = j * (Gu + 1) + i, b = a + 1, c = a + Gu + 1, d = c + 1;
                if (!(inside[a] || inside[b] || inside[c] || inside[d] || crossed[j * Gu + i]))
                    continue;
                for (const std::array<std::size_t, 3>& tri :
                     {std::array<std::size_t, 3>{a, b, d}, std::array<std::size_t, 3>{a, d, c}}) {
                    const std::int64_t base = static_cast<std::int64_t>(cad.mPoints.size());
                    for (std::size_t q : tri)
                        cad.mPoints.push_back(node[q]);
                    cad.mVertices.push_back({base, base + 1, base + 2});
                    for (std::size_t q : tri)
                        cad.mCorners.push_back(node[q]);
                    cad.mSourceCell.push_back(static_cast<std::int64_t>(cad_face.size()));
                    cad_uv.push_back(
                        {nu[tri[0]], nv[tri[0]], nu[tri[1]], nv[tri[1]], nu[tri[2]], nv[tri[2]]});
                    cad_face.push_back(fi);
                }
            }
    }
    if (samples.empty() || cad_face.empty()) {
        out.mWhy = "no face could be sampled";
        return out;
    }

    // CAD to mesh.
    std::vector<double> c2m(samples.size());
    {
        const std::vector<ClosestPointHit> hits = query_closest_points(skin_query, samples);
        for (std::size_t s = 0; s < samples.size(); ++s) {
            c2m[s] = hits[s].mFound ? hits[s].mDistance : 0.0;
            BrepFaceDeviation& fd = out.mFaces[sample_face[s]];
            fd.mCadToMeshMax = std::max(fd.mCadToMeshMax, c2m[s]);
            ++fd.mNumSamples;
        }
        std::vector<std::vector<double>> per_face(rBrep.mFaces.size());
        for (std::size_t s = 0; s < samples.size(); ++s)
            per_face[sample_face[s]].push_back(c2m[s]);
        for (std::size_t fi = 0; fi < per_face.size(); ++fi)
            if (!per_face[fi].empty())
                out.mFaces[fi].mCadToMeshP95 = coverage_quantile(per_face[fi], 0.95);
    }

    // Mesh to CAD: the skin's vertices, then its facet centres (whose facet
    // normal gives the normal deviation).
    std::vector<Vec3> queries = skin.mPoints;
    const std::size_t nv = queries.size();
    for (std::size_t t = 0; t < skin.NumTriangles(); ++t) {
        const Vec3& a = skin.mCorners[3 * t];
        const Vec3& b = skin.mCorners[3 * t + 1];
        const Vec3& c = skin.mCorners[3 * t + 2];
        queries.push_back(
            {(a[0] + b[0] + c[0]) / 3.0, (a[1] + b[1] + c[1]) / 3.0, (a[2] + b[2] + c[2]) / 3.0});
    }
    // A skin point no triangle uses (an orphan) is not part of the skin.
    std::vector<char> used(nv, 0);
    for (const std::array<std::int64_t, 3>& v : skin.mVertices)
        for (std::int64_t q : v)
            used[static_cast<std::size_t>(q)] = 1;
    SurfaceDistanceOptions copt;
    copt.mSign = SdfSign::Unsigned;
    copt.mWatertightCheck = SdfWatertightCheck::Off;
    const DistanceQuery cad_query = build_distance_query(cad, copt);
    // Each query resolved against the face of the nearest grid triangle: its
    // foot inside the trims, else its nearest trim point. A foot outside the
    // trims (a point over a face's hole, where the grid's cells overhang it)
    // tries the face of the nearest triangle of another face, up to
    // `kBdTries` faces; the nearest result stands.
    struct BdHit {
        double mDist = std::numeric_limits<double>::infinity();
        Vec3 mNormal{0.0, 0.0, 0.0};
        std::size_t mFace = 0;
        bool mInside = false;
    };
    const auto resolve = [&](const Vec3& rP, const SurfaceProjection& rH) {
        BdHit r;
        const std::size_t t = static_cast<std::size_t>(rH.mTriangle);
        r.mFace = cad_face[t];
        const BdFace& f = faces[r.mFace];
        // The foot's parameters, interpolated over the grid triangle.
        const Vec3& a = cad.mCorners[3 * t];
        const Vec3& b = cad.mCorners[3 * t + 1];
        const Vec3& c = cad.mCorners[3 * t + 2];
        const Vec3 e0 = vec3_sub(b, a), e1 = vec3_sub(c, a), ep = vec3_sub(rH.mPoint, a);
        const double d00 = MathUtils::Dot(e0, e0), d01 = MathUtils::Dot(e0, e1),
                     d11 = MathUtils::Dot(e1, e1);
        const double d20 = MathUtils::Dot(ep, e0), d21 = MathUtils::Dot(ep, e1);
        const double den = d00 * d11 - d01 * d01;
        double wb = 0.0, wc = 0.0;
        if (den > 0.0) {
            wb = (d11 * d20 - d01 * d21) / den;
            wc = (d00 * d21 - d01 * d20) / den;
        }
        const std::array<double, 6>& uv = cad_uv[t];
        double u = (1.0 - wb - wc) * uv[0] + wb * uv[2] + wc * uv[4];
        double v = (1.0 - wb - wc) * uv[1] + wb * uv[3] + wc * uv[5];
        Vec3 n_foot{0.0, 0.0, 0.0};
        const double dn = bd_project(f, rP, rOptions.mMaxNewtonIterations, u, v, n_foot);
        if (bd_inside(f, u, v)) {
            r.mDist = dn;
            r.mNormal = n_foot;
            r.mInside = true;
        }
        if (!f.mBoundary.empty()) {
            double bu = 0.0, bv = 0.0;
            const double db = bd_boundary_distance(f, rP, bu, bv);
            if (db < r.mDist) {
                r.mDist = db;
                Vec3 sp, su, sv;
                bd_eval(f, bu, bv, sp, su, sv);
                r.mNormal = vec3_cross(su, sv);
            }
        }
        if (!std::isfinite(r.mDist)) {
            r.mDist = rH.mDistance;
            r.mNormal = rH.mNormal;
        }
        return r;
    };
    constexpr int kBdTries = 3;
    std::vector<BdHit> best(queries.size());
    std::vector<std::vector<std::size_t>> tried(queries.size());
    std::vector<std::size_t> pending;
    // A query measured to one face: its foot from the face's nearest sample
    // node, kept inside the trims, else the nearest trim point.
    const auto on_face = [&](const Vec3& rP, std::size_t Fi) {
        BdHit r;
        r.mFace = Fi;
        const BdFace& f = faces[Fi];
        if (!f.mOk)
            return r;
        if (!f.mNodes.empty()) {
            // Gauss-Newton from the few nearest nodes: one start can stall in
            // a valley of a curved face.
            constexpr std::size_t kStarts = 6;
            std::vector<std::pair<double, std::size_t>> by;
            by.reserve(f.mNodes.size());
            for (std::size_t k = 0; k < f.mNodes.size(); ++k)
                by.emplace_back(MathUtils::Norm(vec3_sub(rP, f.mNodes[k])), k);
            const std::size_t ns = std::min(kStarts, by.size());
            std::partial_sort(by.begin(), by.begin() + static_cast<std::ptrdiff_t>(ns), by.end());
            for (std::size_t i = 0; i < ns; ++i) {
                double u = f.mNodeUv[by[i].second][0], v = f.mNodeUv[by[i].second][1];
                Vec3 n_foot{0.0, 0.0, 0.0};
                const double dn = bd_project(f, rP, rOptions.mMaxNewtonIterations, u, v, n_foot);
                if (dn < r.mDist && bd_inside(f, u, v)) {
                    r.mDist = dn;
                    r.mNormal = n_foot;
                    r.mInside = true;
                }
            }
        }
        if (!f.mBoundary.empty()) {
            double bu = 0.0, bv = 0.0;
            const double db = bd_boundary_distance(f, rP, bu, bv);
            if (db < r.mDist) {
                r.mDist = db;
                Vec3 sp, su, sv;
                bd_eval(f, bu, bv, sp, su, sv);
                r.mNormal = vec3_cross(su, sv);
            }
        }
        return r;
    };
    const bool hinted = pTriFace && pTriFace->size() == skin.NumTriangles();
    if (hinted) {
        // Each vertex's faces, from its facets.
        std::vector<std::vector<std::size_t>> vface(nv);
        for (std::size_t t = 0; t < skin.NumTriangles(); ++t) {
            const std::int64_t fi = (*pTriFace)[t];
            if (fi < 0 || static_cast<std::size_t>(fi) >= faces.size())
                continue;
            for (std::int64_t v : skin.mVertices[t]) {
                auto& list = vface[static_cast<std::size_t>(v)];
                if (std::find(list.begin(), list.end(), static_cast<std::size_t>(fi)) == list.end())
                    list.push_back(static_cast<std::size_t>(fi));
            }
        }
        for (std::size_t q = 0; q < queries.size(); ++q) {
            std::vector<std::size_t> cand;
            if (q < nv) {
                cand = vface[q];
            } else {
                const std::int64_t fi = (*pTriFace)[q - nv];
                if (fi >= 0 && static_cast<std::size_t>(fi) < faces.size())
                    cand.push_back(static_cast<std::size_t>(fi));
            }
            for (std::size_t fi : cand) {
                const BdHit h = on_face(queries[q], fi);
                tried[q].push_back(fi);
                if (h.mDist < best[q].mDist)
                    best[q] = h;
            }
            if (!tried[q].empty() && !std::isfinite(best[q].mDist))
                tried[q].clear();  // its face could not be measured
        }
    } else {
        const std::vector<SurfaceProjection> proj = query_surface_projections(cad_query, queries);
        for (std::size_t q = 0; q < queries.size(); ++q) {
            if ((q < nv && !used[q]) || !proj[q].mFound || proj[q].mTriangle < 0)
                continue;
            best[q] = resolve(queries[q], proj[q]);
            tried[q].push_back(best[q].mFace);
            if (!best[q].mInside)
                pending.push_back(q);
        }
    }
    for (int pass = 1; pass < kBdTries && !pending.empty(); ++pass) {
        std::vector<Vec3> pts;
        for (std::size_t q : pending)
            pts.push_back(queries[q]);
        const std::vector<SurfaceProjection> proj = query_surface_projections_filtered(
            cad_query, pts,
            [&](std::size_t K, std::int64_t T) {
                const std::vector<std::size_t>& no = tried[pending[K]];
                return std::find(no.begin(), no.end(), cad_face[static_cast<std::size_t>(T)]) ==
                       no.end();
            },
            std::numeric_limits<double>::infinity());
        std::vector<std::size_t> next;
        for (std::size_t k = 0; k < pending.size(); ++k) {
            const std::size_t q = pending[k];
            // The grid's triangles sag from a curved face, so their distance
            // bounds nothing: the next face is resolved whatever it is.
            if (!proj[k].mFound || proj[k].mTriangle < 0)
                continue;
            const BdHit h = resolve(queries[q], proj[k]);
            tried[q].push_back(h.mFace);
            if (h.mDist < best[q].mDist)
                best[q] = h;
            if (!h.mInside)
                next.push_back(q);
        }
        pending = std::move(next);
    }
    // What is still not over any face it tried: every face, each from its
    // nearest grid node.
    for (std::size_t q : pending) {
        for (std::size_t fi = 0; fi < faces.size(); ++fi) {
            const BdFace& f = faces[fi];
            if (!f.mOk || f.mNodes.empty() ||
                std::find(tried[q].begin(), tried[q].end(), fi) != tried[q].end())
                continue;
            std::size_t near = 0;
            double nd = std::numeric_limits<double>::infinity();
            for (std::size_t k = 0; k < f.mNodes.size(); ++k) {
                const double d = MathUtils::Norm(vec3_sub(queries[q], f.mNodes[k]));
                if (d < nd) {
                    nd = d;
                    near = k;
                }
            }
            double u = f.mNodeUv[near][0], v = f.mNodeUv[near][1];
            Vec3 n_foot{0.0, 0.0, 0.0};
            const double dn =
                bd_project(f, queries[q], rOptions.mMaxNewtonIterations, u, v, n_foot);
            if (dn < best[q].mDist && bd_inside(f, u, v)) {
                best[q].mDist = dn;
                best[q].mNormal = n_foot;
                best[q].mFace = fi;
                best[q].mInside = true;
            }
        }
    }
    std::vector<double> m2c, normal_deg;
    m2c.reserve(queries.size());
    for (std::size_t q = 0; q < queries.size(); ++q) {
        if (tried[q].empty())
            continue;
        const BdHit& h = best[q];
        m2c.push_back(h.mDist);
        BrepFaceDeviation& fd = out.mFaces[h.mFace];
        fd.mMeshToCadMax = std::max(fd.mMeshToCadMax, h.mDist);
        if (q >= nv) {
            const std::size_t st = q - nv;
            const double ang = bd_angle_deg(
                GeometryUtils::TriangleNormal(skin.mCorners[3 * st], skin.mCorners[3 * st + 1],
                                              skin.mCorners[3 * st + 2]),
                h.mNormal);
            normal_deg.push_back(ang);
            fd.mNormalMaxDeg = std::max(fd.mNormalMaxDeg, ang);
        }
    }

    out.mCadToMeshMax = *std::max_element(c2m.begin(), c2m.end());
    out.mCadToMeshP95 = coverage_quantile(c2m, 0.95);
    out.mCadToMeshP99 = coverage_quantile(c2m, 0.99);
    if (!m2c.empty()) {
        out.mMeshToCadMax = *std::max_element(m2c.begin(), m2c.end());
        out.mMeshToCadP95 = coverage_quantile(m2c, 0.95);
        out.mMeshToCadP99 = coverage_quantile(m2c, 0.99);
    }
    if (!normal_deg.empty()) {
        out.mNormalMaxDeg = *std::max_element(normal_deg.begin(), normal_deg.end());
        out.mNormalP95Deg = coverage_quantile(normal_deg, 0.95);
    }
    out.mOk = true;
    return out;
}

}  // namespace meshioplusplus::detail

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
// parameter-space segment cut into `Sub` pieces, so a straight trim in (u, v)
// follows the surface between its points.
void bd_boundary(BdFace& rF, std::size_t Sub) {
    for (std::size_t l = 0; l < rF.mLu.size(); ++l) {
        std::vector<Vec3> pts;
        std::vector<std::array<double, 2>> uvs;
        const std::size_t n = rF.mLu[l].size();
        for (std::size_t i = 0; i < n; ++i) {
            const std::size_t j = (i + 1) % n;
            for (std::size_t k = 0; k < Sub; ++k) {
                const double t = static_cast<double>(k) / static_cast<double>(Sub);
                const double u = rF.mLu[l][i] + t * (rF.mLu[l][j] - rF.mLu[l][i]);
                const double v = rF.mLv[l][i] + t * (rF.mLv[l][j] - rF.mLv[l][i]);
                Vec3 s, su, sv;
                bd_eval(rF, u, v, s, su, sv);
                pts.push_back(s);
                uvs.push_back({u, v});
            }
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
                             const BrepDeviationOptions& rOptions) {
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
    for (std::size_t fi = 0; fi < rBrep.mFaces.size(); ++fi) {
        faces.push_back(bd_face(rBrep, rBrep.mFaces[fi], rOptions));
        if (faces.back().mOk)
            bd_boundary(faces.back(), 4);
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
        std::vector<Vec3> node((N + 1) * (N + 1));
        std::vector<char> inside((N + 1) * (N + 1), 0);
        std::vector<double> nu((N + 1) * (N + 1)), nv((N + 1) * (N + 1));
        for (std::size_t j = 0; j <= N; ++j)
            for (std::size_t i = 0; i <= N; ++i) {
                const std::size_t k = j * (N + 1) + i;
                nu[k] = lo_u + (hi_u - lo_u) * static_cast<double>(i) / static_cast<double>(N);
                nv[k] = lo_v + (hi_v - lo_v) * static_cast<double>(j) / static_cast<double>(N);
                Vec3 su, sv;
                bd_eval(f, nu[k], nv[k], node[k], su, sv);
                inside[k] = bd_inside(f, nu[k], nv[k]) ? 1 : 0;
                if (inside[k]) {
                    samples.push_back(node[k]);
                    sample_face.push_back(fi);
                }
            }
        // A cell a trim runs through joins too, so the face's triangles reach
        // its boundary where no node of the cell falls inside.
        std::vector<char> crossed(N * N, 0);
        for (std::size_t l = 0; l < f.mBoundaryUv.size(); ++l)
            for (const auto& uv : f.mBoundaryUv[l]) {
                const double fu = (uv[0] - lo_u) / (hi_u - lo_u) * static_cast<double>(N);
                const double fv = (uv[1] - lo_v) / (hi_v - lo_v) * static_cast<double>(N);
                if (!(fu >= 0.0) || !(fv >= 0.0) || fu > static_cast<double>(N) ||
                    fv > static_cast<double>(N))
                    continue;
                const std::size_t ci = std::min(N - 1, static_cast<std::size_t>(fu));
                const std::size_t cj = std::min(N - 1, static_cast<std::size_t>(fv));
                crossed[cj * N + ci] = 1;
            }
        for (std::size_t j = 0; j < N; ++j)
            for (std::size_t i = 0; i < N; ++i) {
                const std::size_t a = j * (N + 1) + i, b = a + 1, c = a + N + 1, d = c + 1;
                if (!(inside[a] || inside[b] || inside[c] || inside[d] || crossed[j * N + i]))
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
    const std::vector<SurfaceProjection> proj = query_surface_projections(cad_query, queries);
    std::vector<double> m2c, normal_deg;
    m2c.reserve(queries.size());
    for (std::size_t q = 0; q < queries.size(); ++q) {
        if (q < nv && !used[q])
            continue;
        const SurfaceProjection& h = proj[q];
        if (!h.mFound || h.mTriangle < 0)
            continue;
        const std::size_t t = static_cast<std::size_t>(h.mTriangle);
        const std::size_t fi = cad_face[t];
        const BdFace& f = faces[fi];
        // The foot's parameters, interpolated over the CAD triangle.
        const Vec3& a = cad.mCorners[3 * t];
        const Vec3& b = cad.mCorners[3 * t + 1];
        const Vec3& c = cad.mCorners[3 * t + 2];
        const Vec3 e0 = vec3_sub(b, a), e1 = vec3_sub(c, a), ep = vec3_sub(h.mPoint, a);
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
        // The nearest point of a trimmed face is its foot inside the trims or
        // a point of its boundary: the triangles only choose the face and the
        // start. An untrimmed face's foot always counts.
        Vec3 n_surface{0.0, 0.0, 0.0};
        double d = std::numeric_limits<double>::infinity();
        Vec3 n_foot{0.0, 0.0, 0.0};
        const double dn = bd_project(f, queries[q], rOptions.mMaxNewtonIterations, u, v, n_foot);
        if (bd_inside(f, u, v)) {
            d = dn;
            n_surface = n_foot;
        }
        if (!f.mBoundary.empty()) {
            double bu = 0.0, bv = 0.0;
            const double db = bd_boundary_distance(f, queries[q], bu, bv);
            if (db < d) {
                d = db;
                Vec3 s, su, sv;
                bd_eval(f, bu, bv, s, su, sv);
                n_surface = vec3_cross(su, sv);
            }
        }
        if (!std::isfinite(d)) {
            d = h.mDistance;
            n_surface = h.mNormal;
        }
        m2c.push_back(d);
        BrepFaceDeviation& fd = out.mFaces[fi];
        fd.mMeshToCadMax = std::max(fd.mMeshToCadMax, d);
        if (q >= nv) {
            const std::size_t st = q - nv;
            const double ang = bd_angle_deg(
                GeometryUtils::TriangleNormal(skin.mCorners[3 * st], skin.mCorners[3 * st + 1],
                                              skin.mCorners[3 * st + 2]),
                n_surface);
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

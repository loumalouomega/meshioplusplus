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
#include <cmath>

// Project includes
#include "meshioplusplus/detail/surface_distance.hpp"
#include "meshioplusplus/operations/sdf.hpp"
#include "lsq_fit.hpp"
#include "nurbs_eval.hpp"
#include "param_disc.hpp"
#include "surface_validate.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

// The median of a non-empty vector (the upper one of an even count).
double sv_median(std::vector<double> rV) {
    const std::size_t mid = rV.size() / 2;
    std::nth_element(rV.begin(), rV.begin() + static_cast<std::ptrdiff_t>(mid), rV.end());
    return rV[mid];
}

// Whether (U, V) lies inside the closed polygon (even-odd rule).
bool sv_inside(const std::vector<double>& rPu, const std::vector<double>& rPv, double U, double V) {
    bool in = false;
    const std::size_t n = rPu.size();
    for (std::size_t i = 0, j = n - 1; i < n; j = i++) {
        if ((rPv[i] > V) != (rPv[j] > V) &&
            U < (rPu[j] - rPu[i]) * (V - rPv[i]) / (rPv[j] - rPv[i]) + rPu[i])
            in = !in;
    }
    return in;
}

}  // namespace

double patch_noise_sigma(const std::vector<PatchVec>& rPts, const std::vector<PatchTri>& rTris) {
    std::vector<std::vector<std::size_t>> ring(rPts.size());
    for (const PatchTri& t : rTris)
        for (int k = 0; k < 3; ++k) {
            ring[t[static_cast<std::size_t>(k)]].push_back(
                t[static_cast<std::size_t>((k + 1) % 3)]);
            ring[t[static_cast<std::size_t>((k + 1) % 3)]].push_back(
                t[static_cast<std::size_t>(k)]);
        }
    for (auto& r : ring) {
        std::sort(r.begin(), r.end());
        r.erase(std::unique(r.begin(), r.end()), r.end());
    }
    bool ok = true;
    std::vector<std::uint8_t> boundary(rPts.size(), 0);
    for (const auto& loop : param_boundary_loops(rTris, ok))
        for (std::size_t v : loop)
            boundary[v] = 1;
    std::vector<double> dev, sizes, plane, ring_sizes;
    std::vector<std::size_t> near;
    for (std::size_t v = 0; v < rPts.size(); ++v) {
        // Interior vertices whose one ring is interior too: a boundary ring
        // is one-sided.
        if (ring[v].size() < 4 || boundary[v])
            continue;
        // The two-ring, without the vertex itself.
        near = ring[v];
        for (std::size_t w : ring[v])
            near.insert(near.end(), ring[w].begin(), ring[w].end());
        std::sort(near.begin(), near.end());
        near.erase(std::unique(near.begin(), near.end()), near.end());
        near.erase(std::remove(near.begin(), near.end(), v), near.end());
        if (near.size() < 9)
            continue;
        // A frame: the normal of the one ring's best-fit plane.
        double c[3] = {0.0, 0.0, 0.0};
        for (std::size_t w : ring[v])
            for (int j = 0; j < 3; ++j)
                c[j] += rPts[w][static_cast<std::size_t>(j)] / static_cast<double>(ring[v].size());
        double cov[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        for (std::size_t w : ring[v]) {
            const double d[3] = {rPts[w][0] - c[0], rPts[w][1] - c[1], rPts[w][2] - c[2]};
            for (int a = 0; a < 3; ++a)
                for (int b = 0; b < 3; ++b)
                    cov[a * 3 + b] += d[a] * d[b];
        }
        double evals[3], evecs[9];
        lsq_sym3_eigen(cov, evals, evecs);  // ascending: the first is the normal
        const double* n = evecs;
        plane.push_back(std::fabs((rPts[v][0] - c[0]) * n[0] + (rPts[v][1] - c[1]) * n[1] +
                                  (rPts[v][2] - c[2]) * n[2]));
        ring_sizes.push_back(static_cast<double>(ring[v].size()));
        const double* e1 = evecs + 6;
        const double e2[3] = {n[1] * e1[2] - n[2] * e1[1], n[2] * e1[0] - n[0] * e1[2],
                              n[0] * e1[1] - n[1] * e1[0]};
        // The height quadric z = a x^2 + b x y + c y^2 + d x + e y + f through
        // the two-ring, by least squares; the vertex's residual from it.
        std::vector<double> ata(36, 0.0), atb(6, 0.0);
        for (std::size_t w : near) {
            const double q[3] = {rPts[w][0] - rPts[v][0], rPts[w][1] - rPts[v][1],
                                 rPts[w][2] - rPts[v][2]};
            const double x = q[0] * e1[0] + q[1] * e1[1] + q[2] * e1[2];
            const double y = q[0] * e2[0] + q[1] * e2[1] + q[2] * e2[2];
            const double z = q[0] * n[0] + q[1] * n[1] + q[2] * n[2];
            const double row[6] = {x * x, x * y, y * y, x, y, 1.0};
            for (int a = 0; a < 6; ++a) {
                atb[static_cast<std::size_t>(a)] += row[a] * z;
                for (int b = 0; b < 6; ++b)
                    ata[static_cast<std::size_t>(a * 6 + b)] += row[a] * row[b];
            }
        }
        if (!lsq_solve(ata, atb, 6))
            continue;
        // The vertex sits at the frame's origin: its residual is the constant.
        dev.push_back(std::fabs(atb[5]));
        sizes.push_back(static_cast<double>(near.size()));
    }
    if (dev.size() < 8)
        return 0.0;
    // The residual's deviation is sigma sqrt(1 + 6 / k) for k neighbours
    // fitting six coefficients; the median absolute value over 0.6745.
    const double k = sv_median(sizes);
    const double sigma = sv_median(dev) / 0.6745 / std::sqrt(1.0 + 6.0 / k);
    // Noise reads the same at both scales: the one ring's plane deviation is
    // sigma sqrt(1 + 1 / k1) too. Curvature the tessellation under-resolves
    // reads far larger against the plane than against the quadric, and is
    // not noise: such a patch is taken as clean (its faces held to the
    // tolerance itself, the conservative side).
    const double k1 = sv_median(ring_sizes);
    const double sigma_plane = sv_median(plane) / 0.6745 / std::sqrt(1.0 + 1.0 / k1);
    return sigma_plane > 1.5 * sigma ? 0.0 : sigma;
}

SurfaceValidation surface_validate(const BsplineSurface& rS, const std::vector<PatchVec>& rPts,
                                   const std::vector<PatchTri>& rTris,
                                   const std::vector<double>& rU, const std::vector<double>& rV,
                                   int SamplesPerSpan, const std::vector<PatchTri>* pDistanceTris) {
    SurfaceValidation out;
    if (rTris.empty() || rS.NumPoles() == 0)
        return out;

    // The patch (or the triangles given) as a soup to measure distances to.
    const std::vector<PatchTri>& dist_tris = pDistanceTris ? *pDistanceTris : rTris;
    TriangleSoup soup;
    soup.mPoints.reserve(rPts.size());
    for (const PatchVec& p : rPts)
        soup.mPoints.push_back({p[0], p[1], p[2]});
    for (std::size_t t = 0; t < dist_tris.size(); ++t) {
        soup.mVertices.push_back({static_cast<std::int64_t>(dist_tris[t][0]),
                                  static_cast<std::int64_t>(dist_tris[t][1]),
                                  static_cast<std::int64_t>(dist_tris[t][2])});
        for (std::size_t c : dist_tris[t])
            soup.mCorners.push_back({rPts[c][0], rPts[c][1], rPts[c][2]});
        soup.mSourceCell.push_back(static_cast<std::int64_t>(t));
    }
    SurfaceDistanceOptions dopt;
    dopt.mSign = SdfSign::Unsigned;
    const DistanceQuery query = build_distance_query(soup, dopt);

    // The trim: the patch's boundary loops in parameter space (the longest
    // loop the outer one; a sample inside an inner loop is outside the trim).
    bool ok = true;
    const std::vector<std::vector<std::size_t>> loops = param_boundary_loops(rTris, ok);
    std::vector<std::vector<double>> lu, lv;
    for (const auto& loop : loops) {
        std::vector<double> pu, pv;
        for (std::size_t v : loop) {
            pu.push_back(rU[v]);
            pv.push_back(rV[v]);
        }
        lu.push_back(std::move(pu));
        lv.push_back(std::move(pv));
    }
    const auto in_trim = [&](double U, double V) {
        if (lu.empty())
            return true;  // a closed patch: everywhere
        bool in = false;
        for (std::size_t i = 0; i < lu.size(); ++i)
            in = in != sv_inside(lu[i], lv[i], U, V);
        return in;
    };

    // Forward: Greville points and a K x K grid per span, inside the trim.
    std::vector<Vec3> samples;
    const auto add = [&](double U, double V) {
        if (!in_trim(U, V))
            return;
        BsplineEval e;
        bspline_eval(rS, U, V, 0, e);
        samples.push_back({e.mS[0], e.mS[1], e.mS[2]});
    };
    const auto greville = [](const std::vector<double>& rK, int P, std::size_t N) {
        std::vector<double> g(N, 0.0);
        for (std::size_t i = 0; i < N; ++i) {
            for (int k = 1; k <= P; ++k)
                g[i] += rK[i + static_cast<std::size_t>(k)];
            g[i] /= static_cast<double>(P);
        }
        return g;
    };
    const std::vector<double> gu = greville(rS.mKnotsU, rS.mDegreeU, rS.mNu);
    const std::vector<double> gv = greville(rS.mKnotsV, rS.mDegreeV, rS.mNv);
    for (double v : gv)
        for (double u : gu)
            add(u, v);
    const int k = std::max(1, SamplesPerSpan);
    const std::vector<std::size_t> su = nurbs_spans(rS.mDegreeU, rS.mKnotsU, rS.mNu);
    const std::vector<std::size_t> sv = nurbs_spans(rS.mDegreeV, rS.mKnotsV, rS.mNv);
    for (std::size_t jv : sv)
        for (int b = 0; b < k; ++b) {
            const double v = rS.mKnotsV[jv] + (b + 0.5) / k * (rS.mKnotsV[jv + 1] - rS.mKnotsV[jv]);
            for (std::size_t ju : su)
                for (int a = 0; a < k; ++a)
                    add(rS.mKnotsU[ju] + (a + 0.5) / k * (rS.mKnotsU[ju + 1] - rS.mKnotsU[ju]), v);
        }
    if (!samples.empty()) {
        const std::vector<DistanceHit> hits = query_distances(query, samples, dopt);
        double sum = 0.0;
        for (const DistanceHit& h : hits) {
            const double d = std::fabs(h.mSignedDistance);
            out.mMaxForward = std::max(out.mMaxForward, d);
            sum += d * d;
        }
        out.mNumSamples = samples.size();
        out.mRmsForward = std::sqrt(sum / static_cast<double>(samples.size()));
    }

    // Reverse: every vertex the triangles use, projected onto the surface.
    std::vector<std::uint8_t> used(rPts.size(), 0);
    for (const PatchTri& t : rTris)
        for (std::size_t v : t)
            used[v] = 1;
    double sum = 0.0;
    std::size_t n = 0;
    for (std::size_t v = 0; v < rPts.size(); ++v) {
        if (!used[v])
            continue;
        double u = rU[v], w = rV[v];
        const double d = bspline_project(rS, rPts[v].data(), u, w);
        out.mMaxReverse = std::max(out.mMaxReverse, d);
        sum += d * d;
        ++n;
    }
    if (n > 0)
        out.mRmsReverse = std::sqrt(sum / static_cast<double>(n));
    return out;
}

}  // namespace detail
}  // namespace meshioplusplus

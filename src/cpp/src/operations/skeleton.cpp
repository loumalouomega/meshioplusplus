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

#include "meshioplusplus/operations/skeleton.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <queue>
#include <set>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/surface_distance.hpp"
#include "meshioplusplus/mesh_api.hpp"

#include "../detail/sparse_solve.hpp"
#include "../detail/surface_fem.hpp"

namespace meshioplusplus {

namespace {

using detail::Vec3;

double sk_signed_volume(const std::vector<Vec3>& rPts,
                        const std::vector<std::array<std::int64_t, 3>>& rTris) {
    double vol6 = 0.0;
    for (const auto& t : rTris) {
        const Vec3& a = rPts[static_cast<std::size_t>(t[0])];
        const Vec3& b = rPts[static_cast<std::size_t>(t[1])];
        const Vec3& c = rPts[static_cast<std::size_t>(t[2])];
        vol6 += detail::triple_product(a, b, c);
    }
    return vol6 / 6.0;
}

void sk_assemble_cot_laplacian(const std::vector<Vec3>& rPts,
                               const std::vector<std::array<std::int64_t, 3>>& rTris,
                               std::vector<double>& rArea,
                               std::vector<detail::SparseTriplet>& rLap) {
    const std::size_t n = rPts.size();
    rArea.assign(n, 0.0);
    rLap.clear();
    rLap.reserve(rTris.size() * 9);
    for (const auto& t : rTris) {
        const auto i0 = static_cast<std::size_t>(t[0]);
        const auto i1 = static_cast<std::size_t>(t[1]);
        const auto i2 = static_cast<std::size_t>(t[2]);
        const Vec3 p[3] = {rPts[i0], rPts[i1], rPts[i2]};
        const Vec3 e0 = detail::vec3_sub(p[2], p[1]);
        const Vec3 e1 = detail::vec3_sub(p[0], p[2]);
        const Vec3 e2 = detail::vec3_sub(p[1], p[0]);
        const double dbl_area = detail::vec3_norm(detail::vec3_cross(e2, detail::vec3_scale(e1, -1.0)));
        if (!(dbl_area > 1e-18))
            continue;
        const double a_third = dbl_area / 6.0;
        rArea[i0] += a_third;
        rArea[i1] += a_third;
        rArea[i2] += a_third;

        // Clamp cotangent weights to [0, 1e4] for unconditional SPD stability during contraction.
        auto cot_w = [&](const Vec3& rU, const Vec3& rV) {
            const double c = detail::vec3_dot(rU, rV) / dbl_area;
            return std::clamp(0.5 * c, 0.0, 1e4);
        };
        const double w0 = cot_w(detail::vec3_sub(p[1], p[0]), detail::vec3_sub(p[2], p[0]));  // opp (1,2)
        const double w1 = cot_w(detail::vec3_sub(p[2], p[1]), detail::vec3_sub(p[0], p[1]));  // opp (2,0)
        const double w2 = cot_w(detail::vec3_sub(p[0], p[2]), detail::vec3_sub(p[1], p[2]));  // opp (0,1)

        const std::size_t idx[3] = {i0, i1, i2};
        const double w[3] = {w0, w1, w2};
        for (int k = 0; k < 3; ++k) {
            const std::size_t a = idx[(k + 1) % 3];
            const std::size_t b = idx[(k + 2) % 3];
            const double wk = w[k];
            if (wk > 0.0) {
                rLap.push_back({a, a, wk});
                rLap.push_back({b, b, wk});
                rLap.push_back({a, b, -wk});
                rLap.push_back({b, a, -wk});
            }
        }
    }
}

struct SkUf {
    std::vector<std::int64_t> mParent;
    explicit SkUf(std::size_t N) : mParent(N) { std::iota(mParent.begin(), mParent.end(), 0); }
    std::int64_t Find(std::int64_t X) {
        while (mParent[static_cast<std::size_t>(X)] != X) {
            mParent[static_cast<std::size_t>(X)] =
                mParent[static_cast<std::size_t>(mParent[static_cast<std::size_t>(X)])];
            X = mParent[static_cast<std::size_t>(X)];
        }
        return X;
    }
    void UniteInto(std::int64_t Keep, std::int64_t Kill) {
        Keep = Find(Keep);
        Kill = Find(Kill);
        if (Keep != Kill)
            mParent[static_cast<std::size_t>(Kill)] = Keep;
    }
};

}  // namespace

SkeletonResult extract_skeleton(const Mesh& rMesh, const SkeletonOptions& rOptions) {
    if (rOptions.mMaxIterations < 1)
        throw std::invalid_argument("extract_skeleton: max_iterations must be >= 1");
    if (!(rOptions.mContractionWeight > 0.0) || !std::isfinite(rOptions.mContractionWeight))
        throw std::invalid_argument("extract_skeleton: contraction_weight must be > 0");
    if (!(rOptions.mAttractionWeight > 0.0) || !std::isfinite(rOptions.mAttractionWeight))
        throw std::invalid_argument("extract_skeleton: attraction_weight must be > 0");
    if (!(rOptions.mVolumeTolerance > 0.0) || !std::isfinite(rOptions.mVolumeTolerance))
        throw std::invalid_argument("extract_skeleton: volume_tolerance must be > 0");

    const detail::SurfaceFem fem = detail::build_surface_fem(rMesh, "extract_skeleton");
    const std::size_t n = fem.mPoints.size();
    std::vector<std::array<std::int64_t, 3>> tris;
    tris.reserve(fem.mTriangles.size());
    for (const auto& t : fem.mTriangles)
        tris.push_back(t.mNodes);
    if (tris.empty())
        throw std::invalid_argument("extract_skeleton: mesh has no valid triangles");

    const std::vector<Vec3> orig_xyz = fem.mPoints;
    std::vector<Vec3> xyz = orig_xyz;

    const double vol0 = std::abs(sk_signed_volume(orig_xyz, tris));
    std::vector<double> area0, area_cur;
    std::vector<detail::SparseTriplet> lap_trip;
    sk_assemble_cot_laplacian(orig_xyz, tris, area0, lap_trip);

    // 1. Implicit mean-curvature flow contraction: (diag(w_H * M) + w_L * L) X' = diag(w_H * M) X.
    double w_l = rOptions.mContractionWeight * fem.mMeanEdge * fem.mMeanEdge;
    double vol_ratio = 1.0;
    std::int64_t iters_run = 0;

    for (std::int64_t iter = 0; iter < rOptions.mMaxIterations; ++iter) {
        sk_assemble_cot_laplacian(xyz, tris, area_cur, lap_trip);
        std::vector<detail::SparseTriplet> sys;
        sys.reserve(lap_trip.size() + n);
        for (const auto& t : lap_trip)
            sys.push_back({t.mRow, t.mCol, w_l * t.mValue});

        std::array<std::vector<double>, 3> rhs;
        std::array<std::vector<double>, 3> sol;
        for (int d = 0; d < 3; ++d) {
            rhs[static_cast<std::size_t>(d)].resize(n);
            sol[static_cast<std::size_t>(d)].resize(n);
        }

        for (std::size_t i = 0; i < n; ++i) {
            const double a0 = std::max(area0[i], 1e-18);
            const double ac = std::max(area_cur[i], 1e-18);
            // Medial-axis attraction increases as local area collapses.
            const double wh = rOptions.mAttractionWeight * std::min(std::sqrt(a0 / ac), 1e3);
            const double diag_m = (area0[i] > 0.0) ? (wh * ac) : 1.0;
            sys.push_back({i, i, diag_m});
            for (int d = 0; d < 3; ++d) {
                rhs[static_cast<std::size_t>(d)][i] = diag_m * xyz[i][static_cast<std::size_t>(d)];
                sol[static_cast<std::size_t>(d)][i] = xyz[i][static_cast<std::size_t>(d)];
            }
        }

        const detail::SparseCsr a_csr = detail::csr_from_triplets(n, n, std::move(sys));
        detail::CsrLdlt ldlt;
        const double* rhs_ptrs[3] = {rhs[0].data(), rhs[1].data(), rhs[2].data()};
        double* sol_ptrs[3] = {sol[0].data(), sol[1].data(), sol[2].data()};
        if (ldlt.Factor(a_csr)) {
            ldlt.SolveBlock(3, rhs_ptrs, sol_ptrs);
        } else {
            detail::CsrIc0 precond;
            precond.Build(a_csr);
            detail::CsrPcgWorkspace ws;
            detail::csr_ic0_pcg_block(a_csr, precond, 3, rhs_ptrs, sol_ptrs, ws, 1e-10, n + 200);
        }

        for (std::size_t i = 0; i < n; ++i)
            xyz[i] = {sol[0][i], sol[1][i], sol[2][i]};

        ++iters_run;
        const double vol_cur = std::abs(sk_signed_volume(xyz, tris));
        vol_ratio = vol0 > 1e-30 ? (vol_cur / vol0) : 0.0;
        if (vol_ratio <= rOptions.mVolumeTolerance)
            break;
        w_l *= 2.0;
    }

    // 2. Connectivity thinning: collapse shortest triangle edges until 0 live triangles remain.
    SkUf uf(n);
    std::vector<std::uint8_t> tri_alive(tris.size(), 1);
    std::vector<std::vector<std::size_t>> v_tris(n);
    for (std::size_t f = 0; f < tris.size(); ++f)
        for (int k = 0; k < 3; ++k)
            v_tris[static_cast<std::size_t>(tris[f][static_cast<std::size_t>(k)])].push_back(f);

    // Deterministic min-heap: (edge_length, u, v, version_u, version_v).
    using HeapItem = std::tuple<double, std::int64_t, std::int64_t, std::uint32_t, std::uint32_t>;
    std::priority_queue<HeapItem, std::vector<HeapItem>, std::greater<HeapItem>> pq;
    std::vector<std::uint32_t> ver(n, 0);

    auto push_edge = [&](std::int64_t a, std::int64_t b) {
        a = uf.Find(a);
        b = uf.Find(b);
        if (a == b)
            return;
        if (a > b)
            std::swap(a, b);
        const double len = detail::vec3_norm(
            detail::vec3_sub(xyz[static_cast<std::size_t>(a)], xyz[static_cast<std::size_t>(b)]));
        pq.push({len, a, b, ver[static_cast<std::size_t>(a)], ver[static_cast<std::size_t>(b)]});
    };

    for (std::size_t f = 0; f < tris.size(); ++f) {
        push_edge(tris[f][0], tris[f][1]);
        push_edge(tris[f][1], tris[f][2]);
        push_edge(tris[f][2], tris[f][0]);
    }

    // Track 1-D graph edges created when triangles collapse.
    std::set<std::pair<std::int64_t, std::int64_t>> graph_edges;
    for (const auto& t : tris) {
        for (int k = 0; k < 3; ++k) {
            std::int64_t a = t[static_cast<std::size_t>(k)];
            std::int64_t b = t[static_cast<std::size_t>((k + 1) % 3)];
            if (a > b)
                std::swap(a, b);
            graph_edges.insert({a, b});
        }
    }

    std::size_t live_tris = tris.size();
    while (!pq.empty() && live_tris > 0) {
        const auto [len, u0, v0, vu, vv] = pq.top();
        pq.pop();
        (void)len;
        const std::int64_t u = uf.Find(u0);
        const std::int64_t v = uf.Find(v0);
        if (u == v || u != u0 || v != v0)
            continue;
        if (ver[static_cast<std::size_t>(u)] != vu || ver[static_cast<std::size_t>(v)] != vv)
            continue;

        // Only collapse if edge (u, v) belongs to at least one live triangle.
        bool in_live_tri = false;
        for (std::size_t f : v_tris[static_cast<std::size_t>(u)]) {
            if (!tri_alive[f])
                continue;
            const std::int64_t a = uf.Find(tris[f][0]);
            const std::int64_t b = uf.Find(tris[f][1]);
            const std::int64_t c = uf.Find(tris[f][2]);
            if ((a == u || b == u || c == u) && (a == v || b == v || c == v)) {
                in_live_tri = true;
                break;
            }
        }
        if (!in_live_tri)
            continue;

        // Collapse v -> u, midpoint placement.
        xyz[static_cast<std::size_t>(u)] = detail::vec3_scale(
            detail::vec3_add(xyz[static_cast<std::size_t>(u)], xyz[static_cast<std::size_t>(v)]),
            0.5);
        uf.UniteInto(u, v);
        ++ver[static_cast<std::size_t>(u)];

        // Migrate v's triangles to u, killing degenerate ones.
        auto& ut = v_tris[static_cast<std::size_t>(u)];
        for (std::size_t f : v_tris[static_cast<std::size_t>(v)]) {
            if (!tri_alive[f])
                continue;
            const std::int64_t a = uf.Find(tris[f][0]);
            const std::int64_t b = uf.Find(tris[f][1]);
            const std::int64_t c = uf.Find(tris[f][2]);
            if (a == b || b == c || c == a) {
                tri_alive[f] = 0;
                --live_tris;
            } else {
                ut.push_back(f);
            }
        }
        v_tris[static_cast<std::size_t>(v)].clear();

        // Re-check u's triangles and push updated edges.
        for (std::size_t f : ut) {
            if (!tri_alive[f])
                continue;
            const std::int64_t a = uf.Find(tris[f][0]);
            const std::int64_t b = uf.Find(tris[f][1]);
            const std::int64_t c = uf.Find(tris[f][2]);
            if (a == b || b == c || c == a) {
                tri_alive[f] = 0;
                --live_tris;
            } else {
                push_edge(a, b);
                push_edge(b, c);
                push_edge(c, a);
            }
        }
    }

    // 3. Optional short-edge merge on the 1-D graph.
    const double merge_tol =
        rOptions.mMergeTolerance > 0.0 ? rOptions.mMergeTolerance : 0.25 * fem.mMeanEdge;
    if (merge_tol > 0.0) {
        bool merged = true;
        while (merged) {
            merged = false;
            std::set<std::pair<std::int64_t, std::int64_t>> remapped;
            for (const auto& [a0, b0] : graph_edges) {
                std::int64_t a = uf.Find(a0);
                std::int64_t b = uf.Find(b0);
                if (a == b)
                    continue;
                if (a > b)
                    std::swap(a, b);
                remapped.insert({a, b});
            }
            graph_edges.swap(remapped);
            for (const auto& [a, b] : graph_edges) {
                const double len = detail::vec3_norm(detail::vec3_sub(
                    xyz[static_cast<std::size_t>(a)], xyz[static_cast<std::size_t>(b)]));
                if (len < merge_tol) {
                    xyz[static_cast<std::size_t>(a)] = detail::vec3_scale(
                        detail::vec3_add(xyz[static_cast<std::size_t>(a)],
                                         xyz[static_cast<std::size_t>(b)]),
                        0.5);
                    uf.UniteInto(a, b);
                    merged = true;
                    break;
                }
            }
        }
    }

    // 4. Compact surviving skeleton vertices & edges.
    std::vector<std::int64_t> old_to_new(n, -1);
    std::vector<Vec3> sk_pts;
    std::vector<std::pair<std::int64_t, std::int64_t>> sk_edges;

    // Collect deduplicated edges over canonical roots.
    std::set<std::pair<std::int64_t, std::int64_t>> final_edges;
    for (const auto& [a0, b0] : graph_edges) {
        std::int64_t a = uf.Find(a0);
        std::int64_t b = uf.Find(b0);
        if (a == b)
            continue;
        if (a > b)
            std::swap(a, b);
        final_edges.insert({a, b});
    }

    auto ensure_node = [&](std::int64_t root) -> std::int64_t {
        const auto ru = static_cast<std::size_t>(root);
        if (old_to_new[ru] < 0) {
            old_to_new[ru] = static_cast<std::int64_t>(sk_pts.size());
            sk_pts.push_back(xyz[ru]);
        }
        return old_to_new[ru];
    };

    for (const auto& [a, b] : final_edges) {
        const std::int64_t na = ensure_node(a);
        const std::int64_t nb = ensure_node(b);
        sk_edges.emplace_back(na, nb);
    }
    // Ensure every input vertex maps to a skeleton node (even if isolated / collapsed to 1 point).
    for (std::size_t i = 0; i < n; ++i)
        ensure_node(uf.Find(static_cast<std::int64_t>(i)));

    const std::size_t n_sk = sk_pts.size();
    std::vector<std::int64_t> surf_map(n);
    for (std::size_t i = 0; i < n; ++i)
        surf_map[i] = old_to_new[static_cast<std::size_t>(uf.Find(static_cast<std::int64_t>(i)))];

    // 5. Optional recentering & inscribed radius via DistanceQuery on original surface.
    const detail::TriangleSoup soup = detail::build_triangle_soup(rMesh, "");
    const detail::DistanceQuery dq =
        detail::build_distance_query(soup, SurfaceDistanceOptions{});
    if (rOptions.mRecenter && n_sk > 0) {
        std::vector<Vec3> sum_p(n_sk, {0.0, 0.0, 0.0});
        std::vector<std::size_t> cnt_p(n_sk, 0);
        for (std::size_t i = 0; i < n; ++i) {
            const auto sid = static_cast<std::size_t>(surf_map[i]);
            sum_p[sid] = detail::vec3_add(sum_p[sid], orig_xyz[i]);
            ++cnt_p[sid];
        }
        for (std::size_t s = 0; s < n_sk; ++s) {
            if (cnt_p[s] > 0) {
                const Vec3 avg = detail::vec3_scale(sum_p[s], 1.0 / static_cast<double>(cnt_p[s]));
                // Recenter to the centroid of the mapped original surface ring (Au et al. 2008).
                sk_pts[s] = avg;
            }
        }
    }

    const std::vector<detail::ClosestPointHit> hits = detail::query_closest_points(dq, sk_pts);

    // 6. Compute degrees and build output 1-D line Mesh.
    std::vector<std::int64_t> deg(n_sk, 0);
    for (const auto& [a, b] : sk_edges) {
        ++deg[static_cast<std::size_t>(a)];
        ++deg[static_cast<std::size_t>(b)];
    }

    SkeletonResult result;
    result.mContractionIterations = iters_run;
    result.mVolumeRatio = vol_ratio;
    result.mNumNodes = static_cast<std::int64_t>(n_sk);
    result.mNumSegments = static_cast<std::int64_t>(sk_edges.size());
    for (std::int64_t d : deg) {
        if (d == 1)
            ++result.mNumEndpoints;
        else if (d >= 3)
            ++result.mNumBranches;
    }

    Mesh out;
    NDArray pts_arr(DType::Float64, {n_sk, std::size_t{3}});
    double* pp = pts_arr.As<double>();
    for (std::size_t i = 0; i < n_sk; ++i) {
        pp[i * 3 + 0] = sk_pts[i][0];
        pp[i * 3 + 1] = sk_pts[i][1];
        pp[i * 3 + 2] = sk_pts[i][2];
    }
    out.AssignPoints(std::move(pts_arr));

    if (!sk_edges.empty()) {
        NDArray conn(DType::Int64, {sk_edges.size(), std::size_t{2}});
        std::int64_t* pc = conn.As<std::int64_t>();
        for (std::size_t i = 0; i < sk_edges.size(); ++i) {
            pc[i * 2 + 0] = sk_edges[i].first;
            pc[i * 2 + 1] = sk_edges[i].second;
        }
        out.AddCellBlock(cell_type_name(CellType::Line), std::move(conn));
    }

    NDArray rad_arr(DType::Float64, {n_sk});
    double* pr = rad_arr.As<double>();
    for (std::size_t i = 0; i < n_sk; ++i)
        pr[i] = hits[i].mFound ? hits[i].mDistance : 0.0;
    out.AddPointData(kSkeletonRadiusName, std::move(rad_arr));

    NDArray deg_arr(DType::Int64, {n_sk});
    std::copy(deg.begin(), deg.end(), deg_arr.As<std::int64_t>());
    out.AddPointData(kSkeletonDegreeName, std::move(deg_arr));

    NDArray map_arr(DType::Int64, {n});
    std::copy(surf_map.begin(), surf_map.end(), map_arr.As<std::int64_t>());
    out.AddFieldData(kSkeletonSurfaceMapName, std::move(map_arr));

    result.mMesh = std::move(out);
    return result;
}

}  // namespace meshioplusplus

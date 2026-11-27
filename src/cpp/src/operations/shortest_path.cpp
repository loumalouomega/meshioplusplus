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

#include "meshioplusplus/operations/shortest_path.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <queue>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/node_adjacency.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/mesh_api.hpp"
#include "meshioplusplus/operations/geodesic.hpp"

#include "../detail/crease_edges.hpp"
#include "../detail/point_knn.hpp"
#include "../detail/surface_rings.hpp"

namespace meshioplusplus {

namespace {

using detail::Vec3;

std::string sp_lower(const std::string& rS) {
    std::string out = rS;
    for (char& c : out)
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return out;
}

std::vector<Vec3> sp_read_xyz(const Mesh& rMesh) {
    const std::size_t n = rMesh.NumPoints();
    const std::size_t pdim = rMesh.PointDim();
    const NDArray& pts = rMesh.Points();
    std::vector<Vec3> xyz(n);
    for (std::size_t i = 0; i < n; ++i) {
        xyz[i] = {detail::read_double(pts, i * pdim + 0),
                  pdim > 1 ? detail::read_double(pts, i * pdim + 1) : 0.0,
                  pdim > 2 ? detail::read_double(pts, i * pdim + 2) : 0.0};
    }
    return xyz;
}

std::int64_t sp_nearest_vertex(const std::vector<Vec3>& rXyz, const Vec3& rP) {
    double best = std::numeric_limits<double>::infinity();
    std::int64_t best_id = -1;
    for (std::size_t i = 0; i < rXyz.size(); ++i) {
        const double d2 = detail::vec3_norm_sq(detail::vec3_sub(rXyz[i], rP));
        if (d2 < best) {
            best = d2;
            best_id = static_cast<std::int64_t>(i);
        }
    }
    return best_id;
}

Mesh sp_build_line_mesh(const std::vector<Vec3>& rPts, const std::vector<double>& rDists,
                        const std::vector<std::int64_t>& rSrcCells) {
    Mesh m;
    const std::size_t n = rPts.size();
    NDArray pts(DType::Float64, {n, std::size_t{3}});
    double* pp = pts.As<double>();
    for (std::size_t i = 0; i < n; ++i) {
        pp[i * 3 + 0] = rPts[i][0];
        pp[i * 3 + 1] = rPts[i][1];
        pp[i * 3 + 2] = rPts[i][2];
    }
    m.AssignPoints(std::move(pts));

    const std::size_t nseg = n >= 2 ? n - 1 : 0;
    if (nseg > 0) {
        NDArray conn(DType::Int64, {nseg, std::size_t{2}});
        std::int64_t* pc = conn.As<std::int64_t>();
        for (std::size_t i = 0; i < nseg; ++i) {
            pc[i * 2 + 0] = static_cast<std::int64_t>(i);
            pc[i * 2 + 1] = static_cast<std::int64_t>(i + 1);
        }
        m.AddCellBlock(cell_type_name(CellType::Line), std::move(conn));

        NDArray sc(DType::Int64, {nseg});
        std::int64_t* psc = sc.As<std::int64_t>();
        for (std::size_t i = 0; i < nseg; ++i)
            psc[i] = i < rSrcCells.size() ? rSrcCells[i] : -1;
        m.AddCellData(kPathSourceCellName, {std::move(sc)});
    }

    NDArray da(DType::Float64, {n});
    std::copy(rDists.begin(), rDists.end(), da.As<double>());
    m.AddPointData(kPathDistanceName, std::move(da));
    return m;
}

// Lookup an incident global cell id for edge (u, v) using a sorted edge->cell table.
std::int64_t sp_find_edge_cell(const Mesh& rMesh, std::int64_t U, std::int64_t V,
                               const std::vector<std::pair<std::uint64_t, std::int64_t>>& rMap) {
    const std::int64_t lo = std::min(U, V);
    const std::int64_t hi = std::max(U, V);
    const std::uint64_t key = (static_cast<std::uint64_t>(static_cast<std::uint32_t>(lo)) << 32) |
                              static_cast<std::uint64_t>(static_cast<std::uint32_t>(hi));
    auto it = std::lower_bound(rMap.begin(), rMap.end(), key,
                               [](const auto& rP, std::uint64_t K) { return rP.first < K; });
    if (it != rMap.end() && it->first == key)
        return it->second;
    (void)rMesh;
    return -1;
}

std::vector<std::pair<std::uint64_t, std::int64_t>> sp_build_edge_cell_map(const Mesh& rMesh) {
    std::vector<std::pair<std::uint64_t, std::int64_t>> out;
    std::int64_t gcell = 0;
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto cb = rMesh.Cells(b);
        if (cb.IsRagged()) {
            gcell += static_cast<std::int64_t>(cb.NumCells());
            continue;
        }
        const std::size_t k = cb.NodesPerCell();
        const NDArray& conn = cb.Conn();
        for (std::size_t c = 0; c < cb.NumCells(); ++c, ++gcell) {
            for (std::size_t i = 0; i < k; ++i) {
                const std::int64_t u = detail::read_int(conn, c * k + i);
                const std::int64_t v = detail::read_int(conn, c * k + ((i + 1) % k));
                if (u == v)
                    continue;
                const std::int64_t lo = std::min(u, v);
                const std::int64_t hi = std::max(u, v);
                const std::uint64_t key =
                    (static_cast<std::uint64_t>(static_cast<std::uint32_t>(lo)) << 32) |
                    static_cast<std::uint64_t>(static_cast<std::uint32_t>(hi));
                out.emplace_back(key, gcell);
            }
        }
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end(),
                          [](const auto& a, const auto& b) { return a.first == b.first; }),
              out.end());
    return out;
}

}  // namespace

ShortestPathMethod shortest_path_method_from_name(const std::string& rName) {
    const std::string low = sp_lower(rName);
    if (low == "dijkstra")
        return ShortestPathMethod::Dijkstra;
    if (low == "heat_trace" || low == "heat-trace" || low == "heattrace")
        return ShortestPathMethod::HeatTrace;
    throw std::invalid_argument("shortest_path: unknown method '" + rName +
                                "' (expected 'dijkstra' or 'heat_trace')");
}

ShortestPathResult shortest_path(const Mesh& rMesh, const ShortestPathOptions& rOptions) {
    const std::size_t n = rMesh.NumPoints();
    if (n == 0)
        throw std::invalid_argument("shortest_path: mesh has no points");
    const std::vector<Vec3> xyz = sp_read_xyz(rMesh);

    const std::int64_t src =
        rOptions.mHasSourcePoint ? sp_nearest_vertex(xyz, rOptions.mSourcePoint) : rOptions.mSource;
    const std::int64_t dst =
        rOptions.mHasTargetPoint ? sp_nearest_vertex(xyz, rOptions.mTargetPoint) : rOptions.mTarget;
    if (src < 0 || static_cast<std::size_t>(src) >= n)
        throw std::invalid_argument("shortest_path: source vertex out of range");
    if (dst < 0 || static_cast<std::size_t>(dst) >= n)
        throw std::invalid_argument("shortest_path: target vertex out of range");

    ShortestPathResult result;
    result.mSourceVertex = src;
    result.mTargetVertex = dst;

    if (src == dst) {
        result.mReached = true;
        result.mTotalLength = 0.0;
        result.mNumSegments = 0;
        result.mMesh = sp_build_line_mesh({xyz[static_cast<std::size_t>(src)]}, {0.0}, {});
        return result;
    }

    if (rOptions.mMethod == ShortestPathMethod::HeatTrace) {
        // Run heat geodesic_distance from src, then back-trace -grad(d) across triangles.
        GeodesicOptions gopt;
        gopt.mSources = {src};
        gopt.mTimeFactor = rOptions.mTimeFactor;
        const GeodesicResult gres = geodesic_distance(rMesh, gopt);
        const NDArray& d_arr = gres.mMesh.PointData(kGeodesicDistanceName);
        const double* pd = d_arr.As<double>();
        if (!std::isfinite(pd[static_cast<std::size_t>(dst)])) {
            result.mReached = false;
            result.mTotalLength = std::numeric_limits<double>::quiet_NaN();
            return result;
        }

        // Build triangle rings and face-edge adjacency.
        const auto bases = detail::block_bases(rMesh);
        const std::size_t total_cells = bases.empty() ? 0 : static_cast<std::size_t>(bases.back());
        const detail::SurfaceRings rings = detail::surface_rings(
            rMesh, std::vector<std::uint8_t>(total_cells, 1), "shortest_path");
        const std::size_t nf = rings.NumFaces();
        std::vector<std::vector<std::size_t>> v_faces(n);
        for (std::size_t f = 0; f < nf; ++f) {
            for (std::int64_t p = rings.mStart[f]; p < rings.mStart[f + 1]; ++p)
                v_faces[static_cast<std::size_t>(rings.mNodes[static_cast<std::size_t>(p)])]
                    .push_back(f);
        }
        // Map undirected edge -> opposite face for triangles.
        std::vector<std::array<std::int64_t, 3>> opp_face(nf, {-1, -1, -1});
        {
            const auto uses = detail::crease_edge_uses(rings.mStart, rings.mNodes);
            std::size_t i = 0;
            while (i < uses.size()) {
                std::size_t j = i + 1;
                while (j < uses.size() && uses[j].mLo == uses[i].mLo && uses[j].mHi == uses[i].mHi)
                    ++j;
                if (j - i == 2) {
                    const auto f0 = static_cast<std::size_t>(uses[i].mFace);
                    const auto f1 = static_cast<std::size_t>(uses[i + 1].mFace);
                    auto set_opp = [&](std::size_t f, std::int64_t lo, std::int64_t hi,
                                       std::size_t other) {
                        const std::int64_t s = rings.mStart[f];
                        if (rings.mStart[f + 1] - s != 3)
                            return;
                        for (int e = 0; e < 3; ++e) {
                            const std::int64_t a = rings.mNodes[static_cast<std::size_t>(s + e)];
                            const std::int64_t b =
                                rings.mNodes[static_cast<std::size_t>(s + ((e + 1) % 3))];
                            if (std::min(a, b) == lo && std::max(a, b) == hi)
                                opp_face[f][static_cast<std::size_t>(e)] =
                                    static_cast<std::int64_t>(other);
                        }
                    };
                    set_opp(f0, uses[i].mLo, uses[i].mHi, f1);
                    set_opp(f1, uses[i].mLo, uses[i].mHi, f0);
                }
                i = j;
            }
        }

        // Back-trace from dst -> src along -grad(d).
        std::vector<Vec3> rev_pts;
        std::vector<std::int64_t> rev_cells;
        Vec3 cur_pos = xyz[static_cast<std::size_t>(dst)];
        double cur_d = pd[static_cast<std::size_t>(dst)];
        std::int64_t cur_v = dst;
        std::int64_t cur_f = -1;
        rev_pts.push_back(cur_pos);

        const std::size_t max_steps = std::max<std::size_t>(nf * 4, 1000);
        for (std::size_t step = 0; step < max_steps; ++step) {
            if (cur_v == src)
                break;

            // If at a vertex, pick the incident triangle with steepest descent in -grad(d).
            if (cur_v >= 0) {
                const auto& vf = v_faces[static_cast<std::size_t>(cur_v)];
                std::int64_t best_f = -1;
                double best_slope = 0.0;
                Vec3 best_exit = cur_pos;
                std::int64_t best_next_v = -1;
                std::int64_t best_next_f = -1;
                double best_next_d = cur_d;

                for (std::size_t f : vf) {
                    const std::int64_t s = rings.mStart[f];
                    if (rings.mStart[f + 1] - s != 3)
                        continue;
                    std::int64_t v[3] = {rings.mNodes[static_cast<std::size_t>(s)],
                                         rings.mNodes[static_cast<std::size_t>(s + 1)],
                                         rings.mNodes[static_cast<std::size_t>(s + 2)]};
                    int k = (v[0] == cur_v) ? 0 : (v[1] == cur_v ? 1 : 2);
                    const std::int64_t va = v[k];
                    const std::int64_t vb = v[(k + 1) % 3];
                    const std::int64_t vc = v[(k + 2) % 3];
                    if (vb == src || vc == src) {
                        best_f = static_cast<std::int64_t>(f);
                        best_exit = xyz[static_cast<std::size_t>(src)];
                        best_next_v = src;
                        best_next_f = -1;
                        best_next_d = 0.0;
                        best_slope = std::numeric_limits<double>::infinity();
                        break;
                    }
                    // Compute in-plane gradient of d on triangle (va, vb, vc).
                    const Vec3 p0 = xyz[static_cast<std::size_t>(va)];
                    const Vec3 p1 = xyz[static_cast<std::size_t>(vb)];
                    const Vec3 p2 = xyz[static_cast<std::size_t>(vc)];
                    const Vec3 e1 = detail::vec3_sub(p1, p0);
                    const Vec3 e2 = detail::vec3_sub(p2, p0);
                    const Vec3 nrm = detail::vec3_cross(e1, e2);
                    const double a2 = detail::vec3_norm_sq(nrm);
                    if (a2 > 1e-30) {
                        const double inv_a2 = 1.0 / a2;
                        const Vec3 g1 = detail::vec3_scale(detail::vec3_cross(nrm, e2), -inv_a2);
                        const Vec3 g2 = detail::vec3_scale(detail::vec3_cross(nrm, e1), inv_a2);
                        const double d1 =
                            pd[static_cast<std::size_t>(vb)] - pd[static_cast<std::size_t>(va)];
                        const double d2 =
                            pd[static_cast<std::size_t>(vc)] - pd[static_cast<std::size_t>(va)];
                        const Vec3 dir =
                            detail::vec3_scale(detail::vec3_add(detail::vec3_scale(g1, d1),
                                                                detail::vec3_scale(g2, d2)),
                                               -1.0);
                        // Express dir in (e1, e2) coordinates: dir = alpha e1 + beta e2.
                        const double alpha = detail::vec3_dot(dir, g1);
                        const double beta = detail::vec3_dot(dir, g2);
                        if (alpha > 1e-12 && beta > 1e-12) {
                            const double t = 1.0 / (alpha + beta);
                            const double u = beta * t;  // fraction along edge (vb, vc)
                            const Vec3 exit_p = detail::vec3_add(
                                p1, detail::vec3_scale(detail::vec3_sub(p2, p1), u));
                            const double exit_d = (1.0 - u) * pd[static_cast<std::size_t>(vb)] +
                                                  u * pd[static_cast<std::size_t>(vc)];
                            const double seg_l = detail::vec3_norm(detail::vec3_sub(exit_p, p0));
                            const double slope = seg_l > 1e-15 ? (cur_d - exit_d) / seg_l : 0.0;
                            if (slope > best_slope) {
                                best_slope = slope;
                                best_f = static_cast<std::int64_t>(f);
                                best_exit = exit_p;
                                best_next_v = -1;
                                best_next_f = opp_face[f][static_cast<std::size_t>((k + 1) % 3)];
                                best_next_d = exit_d;
                            }
                        }
                    }
                    // Also consider descending along edges (va->vb) and (va->vc).
                    for (std::int64_t nb : {vb, vc}) {
                        const double dn = pd[static_cast<std::size_t>(nb)];
                        const double elen = detail::vec3_norm(
                            detail::vec3_sub(xyz[static_cast<std::size_t>(nb)], p0));
                        const double slope = elen > 1e-15 ? (cur_d - dn) / elen : 0.0;
                        if (slope > best_slope) {
                            best_slope = slope;
                            best_f = static_cast<std::int64_t>(f);
                            best_exit = xyz[static_cast<std::size_t>(nb)];
                            best_next_v = nb;
                            best_next_f = -1;
                            best_next_d = dn;
                        }
                    }
                }
                if (best_f < 0)
                    break;
                rev_pts.push_back(best_exit);
                rev_cells.push_back(rings.mCell[static_cast<std::size_t>(best_f)]);
                cur_pos = best_exit;
                cur_d = best_next_d;
                cur_v = best_next_v;
                cur_f = best_next_f;
                if (cur_v < 0 && cur_f < 0) {
                    // Hit boundary edge: snap to lower-distance vertex of the triangle.
                    const std::int64_t s = rings.mStart[static_cast<std::size_t>(best_f)];
                    std::int64_t best_bv = rings.mNodes[static_cast<std::size_t>(s)];
                    for (int e = 1; e < 3; ++e) {
                        const std::int64_t cand = rings.mNodes[static_cast<std::size_t>(s + e)];
                        if (pd[static_cast<std::size_t>(cand)] <
                            pd[static_cast<std::size_t>(best_bv)])
                            best_bv = cand;
                    }
                    cur_v = best_bv;
                    cur_pos = xyz[static_cast<std::size_t>(cur_v)];
                    cur_d = pd[static_cast<std::size_t>(cur_v)];
                    rev_pts.push_back(cur_pos);
                    rev_cells.push_back(rings.mCell[static_cast<std::size_t>(best_f)]);
                }
            } else {
                // Inside triangle cur_f at an edge point cur_pos; trace -grad(d) to exit edge.
                const std::size_t f = static_cast<std::size_t>(cur_f);
                const std::int64_t s = rings.mStart[f];
                if (rings.mStart[f + 1] - s != 3)
                    break;
                std::int64_t v[3] = {rings.mNodes[static_cast<std::size_t>(s)],
                                     rings.mNodes[static_cast<std::size_t>(s + 1)],
                                     rings.mNodes[static_cast<std::size_t>(s + 2)]};
                // If src is a vertex of this triangle, step directly to src.
                bool hit_src = false;
                for (int e = 0; e < 3; ++e) {
                    if (v[e] == src) {
                        rev_pts.push_back(xyz[static_cast<std::size_t>(src)]);
                        rev_cells.push_back(rings.mCell[f]);
                        cur_v = src;
                        hit_src = true;
                        break;
                    }
                }
                if (hit_src)
                    break;

                const Vec3 p[3] = {xyz[static_cast<std::size_t>(v[0])],
                                   xyz[static_cast<std::size_t>(v[1])],
                                   xyz[static_cast<std::size_t>(v[2])]};
                const Vec3 e1 = detail::vec3_sub(p[1], p[0]);
                const Vec3 e2 = detail::vec3_sub(p[2], p[0]);
                const Vec3 nrm = detail::vec3_cross(e1, e2);
                const double a2 = detail::vec3_norm_sq(nrm);
                bool advanced = false;
                if (a2 > 1e-30) {
                    const double inv_a2 = 1.0 / a2;
                    const Vec3 g1 = detail::vec3_scale(detail::vec3_cross(nrm, e2), -inv_a2);
                    const Vec3 g2 = detail::vec3_scale(detail::vec3_cross(nrm, e1), inv_a2);
                    const Vec3 rel = detail::vec3_sub(cur_pos, p[0]);
                    const double b1 = detail::vec3_dot(rel, g1);
                    const double b2 = detail::vec3_dot(rel, g2);
                    const double b0 = 1.0 - b1 - b2;
                    const double d1 =
                        pd[static_cast<std::size_t>(v[1])] - pd[static_cast<std::size_t>(v[0])];
                    const double d2 =
                        pd[static_cast<std::size_t>(v[2])] - pd[static_cast<std::size_t>(v[0])];
                    const Vec3 dir = detail::vec3_scale(
                        detail::vec3_add(detail::vec3_scale(g1, d1), detail::vec3_scale(g2, d2)),
                        -1.0);
                    const double db1 = detail::vec3_dot(dir, g1);
                    const double db2 = detail::vec3_dot(dir, g2);
                    const double db0 = -(db1 + db2);
                    const double b[3] = {b0, b1, b2};
                    const double db[3] = {db0, db1, db2};
                    double best_t = std::numeric_limits<double>::infinity();
                    int exit_opp = -1;
                    for (int k = 0; k < 3; ++k) {
                        if (db[k] < -1e-12) {
                            const double t = -b[k] / db[k];
                            if (t > 1e-9 && t < best_t) {
                                best_t = t;
                                exit_opp = k;
                            }
                        }
                    }
                    if (exit_opp >= 0 && std::isfinite(best_t)) {
                        const double nb[3] = {b[0] + best_t * db[0], b[1] + best_t * db[1],
                                              b[2] + best_t * db[2]};
                        const Vec3 next_p = {
                            nb[0] * p[0][0] + nb[1] * p[1][0] + nb[2] * p[2][0],
                            nb[0] * p[0][1] + nb[1] * p[1][1] + nb[2] * p[2][1],
                            nb[0] * p[0][2] + nb[1] * p[1][2] + nb[2] * p[2][2],
                        };
                        const double next_d = nb[0] * pd[static_cast<std::size_t>(v[0])] +
                                              nb[1] * pd[static_cast<std::size_t>(v[1])] +
                                              nb[2] * pd[static_cast<std::size_t>(v[2])];
                        if (next_d < cur_d - 1e-14) {
                            rev_pts.push_back(next_p);
                            rev_cells.push_back(rings.mCell[f]);
                            cur_pos = next_p;
                            cur_d = next_d;
                            // Edge opposite vertex exit_opp is edge ((exit_opp + 1) % 3).
                            const int edge_idx = (exit_opp + 1) % 3;
                            const std::int64_t nf_next =
                                opp_face[f][static_cast<std::size_t>(edge_idx)];
                            if (nf_next >= 0) {
                                cur_f = nf_next;
                                cur_v = -1;
                            } else {
                                // Boundary edge: snap to lower-distance endpoint.
                                const std::int64_t ea = v[edge_idx];
                                const std::int64_t eb = v[(edge_idx + 1) % 3];
                                cur_v = pd[static_cast<std::size_t>(ea)] <=
                                                pd[static_cast<std::size_t>(eb)]
                                            ? ea
                                            : eb;
                                cur_f = -1;
                            }
                            advanced = true;
                        }
                    }
                }
                if (!advanced) {
                    // Fallback: step to lowest-distance vertex of cur_f.
                    std::int64_t best_bv = v[0];
                    for (int e = 1; e < 3; ++e)
                        if (pd[static_cast<std::size_t>(v[e])] <
                            pd[static_cast<std::size_t>(best_bv)])
                            best_bv = v[e];
                    if (pd[static_cast<std::size_t>(best_bv)] >= cur_d - 1e-14)
                        break;
                    cur_v = best_bv;
                    cur_pos = xyz[static_cast<std::size_t>(cur_v)];
                    cur_d = pd[static_cast<std::size_t>(cur_v)];
                    cur_f = -1;
                    rev_pts.push_back(cur_pos);
                    rev_cells.push_back(rings.mCell[f]);
                }
            }
        }

        if (cur_v == src) {
            std::reverse(rev_pts.begin(), rev_pts.end());
            std::reverse(rev_cells.begin(), rev_cells.end());
            std::vector<double> cum(rev_pts.size(), 0.0);
            for (std::size_t i = 1; i < rev_pts.size(); ++i)
                cum[i] =
                    cum[i - 1] + detail::vec3_norm(detail::vec3_sub(rev_pts[i], rev_pts[i - 1]));
            result.mReached = true;
            result.mTotalLength = cum.back();
            result.mNumSegments = static_cast<std::int64_t>(rev_pts.size() - 1);
            result.mMesh = sp_build_line_mesh(rev_pts, cum, rev_cells);
            return result;
        }
        // If heat_trace stalled before reaching src, fall through to deterministic Dijkstra.
    }

    // Deterministic edge-graph Dijkstra.
    const detail::NodeAdjacency adj =
        detail::build_node_adjacency(rMesh, n, detail::NodeAdjacencyKind::Edge);
    constexpr double kInf = std::numeric_limits<double>::infinity();
    std::vector<double> dist(n, kInf);
    std::vector<std::int64_t> parent(n, -1);

    using State = std::pair<double, std::int64_t>;
    std::priority_queue<State, std::vector<State>, std::greater<State>> pq;
    dist[static_cast<std::size_t>(src)] = 0.0;
    pq.push({0.0, src});

    while (!pq.empty()) {
        const auto [du, u] = pq.top();
        pq.pop();
        if (du > dist[static_cast<std::size_t>(u)])
            continue;
        if (u == dst)
            break;
        const std::int64_t beg = adj.mXadj[static_cast<std::size_t>(u)];
        const std::int64_t end = adj.mXadj[static_cast<std::size_t>(u) + 1];
        for (std::int64_t e = beg; e < end; ++e) {
            const std::int64_t v = adj.mAdj[static_cast<std::size_t>(e)];
            const double w = detail::vec3_norm(detail::vec3_sub(xyz[static_cast<std::size_t>(v)],
                                                                xyz[static_cast<std::size_t>(u)]));
            const double nd = du + w;
            const std::size_t vu = static_cast<std::size_t>(v);
            if (nd < dist[vu] || (nd == dist[vu] && (parent[vu] < 0 || u < parent[vu]))) {
                dist[vu] = nd;
                parent[vu] = u;
                pq.push({nd, v});
            }
        }
    }

    if (!std::isfinite(dist[static_cast<std::size_t>(dst)])) {
        result.mReached = false;
        result.mTotalLength = std::numeric_limits<double>::quiet_NaN();
        return result;
    }

    std::vector<std::int64_t> path_ids;
    for (std::int64_t v = dst; v >= 0; v = parent[static_cast<std::size_t>(v)]) {
        path_ids.push_back(v);
        if (v == src)
            break;
    }
    std::reverse(path_ids.begin(), path_ids.end());

    const auto edge_map = sp_build_edge_cell_map(rMesh);
    std::vector<Vec3> path_pts;
    std::vector<double> path_d;
    std::vector<std::int64_t> src_cells;
    path_pts.reserve(path_ids.size());
    path_d.reserve(path_ids.size());
    for (std::size_t i = 0; i < path_ids.size(); ++i) {
        const auto vid = static_cast<std::size_t>(path_ids[i]);
        path_pts.push_back(xyz[vid]);
        path_d.push_back(dist[vid]);
        if (i + 1 < path_ids.size())
            src_cells.push_back(sp_find_edge_cell(rMesh, path_ids[i], path_ids[i + 1], edge_map));
    }

    result.mReached = true;
    result.mTotalLength = path_d.back();
    result.mNumSegments = static_cast<std::int64_t>(path_ids.size() - 1);
    result.mMesh = sp_build_line_mesh(path_pts, path_d, src_cells);
    return result;
}

GeodesicVoronoiResult geodesic_voronoi(const Mesh& rMesh, const GeodesicVoronoiOptions& rOptions) {
    const std::size_t n = rMesh.NumPoints();
    if (n == 0)
        throw std::invalid_argument("geodesic_voronoi: mesh has no points");
    if (rOptions.mSeedPoints.size() % 3 != 0)
        throw std::invalid_argument("geodesic_voronoi: seed_points length must be a multiple of 3");

    const std::vector<Vec3> xyz = sp_read_xyz(rMesh);
    std::vector<std::int64_t> seeds;
    seeds.reserve(rOptions.mSources.size() + rOptions.mSeedPoints.size() / 3);
    for (std::int64_t s : rOptions.mSources) {
        if (s < 0 || static_cast<std::size_t>(s) >= n)
            throw std::invalid_argument("geodesic_voronoi: source vertex out of range");
        seeds.push_back(s);
    }
    for (std::size_t i = 0; i < rOptions.mSeedPoints.size(); i += 3) {
        const Vec3 p{rOptions.mSeedPoints[i], rOptions.mSeedPoints[i + 1],
                     rOptions.mSeedPoints[i + 2]};
        seeds.push_back(sp_nearest_vertex(xyz, p));
    }
    if (seeds.empty())
        throw std::invalid_argument(
            "geodesic_voronoi: at least one source or seed_point is required");

    const detail::NodeAdjacency adj =
        detail::build_node_adjacency(rMesh, n, detail::NodeAdjacencyKind::Edge);
    constexpr double kInf = std::numeric_limits<double>::infinity();
    std::vector<double> dist(n, kInf);
    std::vector<std::int64_t> label(n, -1);

    // Priority tuple: (distance, seed_label, vertex_id) for strict determinism.
    using State = std::tuple<double, std::int64_t, std::int64_t>;
    std::priority_queue<State, std::vector<State>, std::greater<State>> pq;
    for (std::size_t sidx = 0; sidx < seeds.size(); ++sidx) {
        const auto v = static_cast<std::size_t>(seeds[sidx]);
        if (dist[v] > 0.0) {
            dist[v] = 0.0;
            label[v] = static_cast<std::int64_t>(sidx);
            pq.push({0.0, static_cast<std::int64_t>(sidx), seeds[sidx]});
        }
    }

    while (!pq.empty()) {
        const auto [du, su, u] = pq.top();
        pq.pop();
        const auto uu = static_cast<std::size_t>(u);
        if (du > dist[uu] || (du == dist[uu] && su > label[uu]))
            continue;
        const std::int64_t beg = adj.mXadj[uu];
        const std::int64_t end = adj.mXadj[uu + 1];
        for (std::int64_t e = beg; e < end; ++e) {
            const std::int64_t v = adj.mAdj[static_cast<std::size_t>(e)];
            const auto vu = static_cast<std::size_t>(v);
            const double w = detail::vec3_norm(detail::vec3_sub(xyz[vu], xyz[uu]));
            const double nd = du + w;
            if (nd < dist[vu] || (nd == dist[vu] && (label[vu] < 0 || su < label[vu]))) {
                dist[vu] = nd;
                label[vu] = su;
                pq.push({nd, su, v});
            }
        }
    }

    GeodesicVoronoiResult result;
    result.mNumSources = static_cast<std::int64_t>(seeds.size());
    result.mMesh = detail::clone_mesh(rMesh);

    NDArray d_arr(DType::Float64, {n});
    NDArray l_arr(DType::Int64, {n});
    double* pd = d_arr.As<double>();
    std::int64_t* pl = l_arr.As<std::int64_t>();
    for (std::size_t i = 0; i < n; ++i) {
        if (std::isfinite(dist[i])) {
            pd[i] = dist[i];
            pl[i] = label[i];
            result.mMaxDistance = std::max(result.mMaxDistance, dist[i]);
        } else {
            pd[i] = std::numeric_limits<double>::quiet_NaN();
            pl[i] = -1;
            ++result.mNumUnreached;
        }
    }
    result.mMesh.AddPointData(kGeodesicDistanceName, std::move(d_arr));
    result.mMesh.AddPointData(kGeodesicNearestSourceName, std::move(l_arr));

    if (rOptions.mAttachCellLabels && rMesh.NumCellBlocks() > 0) {
        std::vector<NDArray> cell_labels;
        cell_labels.reserve(rMesh.NumCellBlocks());
        for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
            const auto cb = rMesh.Cells(b);
            NDArray cl(DType::Int64, {cb.NumCells()});
            std::int64_t* pcl = cl.As<std::int64_t>();
            for (std::size_t c = 0; c < cb.NumCells(); ++c) {
                std::int64_t best_lbl = -1;
                double best_d = kInf;
                if (cb.IsRagged()) {
                    const std::int64_t* row = cb.Row(c);
                    const std::size_t sz = cb.RowSize(c);
                    for (std::size_t i = 0; i < sz; ++i) {
                        const auto vid = static_cast<std::size_t>(row[i]);
                        if (label[vid] >= 0 &&
                            (dist[vid] < best_d ||
                             (dist[vid] == best_d && (best_lbl < 0 || label[vid] < best_lbl)))) {
                            best_d = dist[vid];
                            best_lbl = label[vid];
                        }
                    }
                } else {
                    const std::size_t k = cb.NodesPerCell();
                    const NDArray& conn = cb.Conn();
                    for (std::size_t i = 0; i < k; ++i) {
                        const auto vid =
                            static_cast<std::size_t>(detail::read_int(conn, c * k + i));
                        if (label[vid] >= 0 &&
                            (dist[vid] < best_d ||
                             (dist[vid] == best_d && (best_lbl < 0 || label[vid] < best_lbl)))) {
                            best_d = dist[vid];
                            best_lbl = label[vid];
                        }
                    }
                }
                pcl[c] = best_lbl;
            }
            cell_labels.push_back(std::move(cl));
        }
        result.mMesh.AddCellData(kGeodesicNearestSourceName, std::move(cell_labels));
    }

    return result;
}

}  // namespace meshioplusplus

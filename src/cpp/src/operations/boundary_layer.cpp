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
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "../detail/layer_sweep.hpp"
#include "../detail/predicates.hpp"
#include "../detail/tri_bvh.hpp"
#include "../detail/tri_tri.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/surface_distance.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/operations/boundary_layer.hpp"
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/region.hpp"

namespace meshioplusplus {

namespace {

using BlV3 = std::array<double, 3>;

constexpr const char* kBoundaryLayerPrefix = "meshio++: boundary_layer: ";
// The least cosine at which a point still sees all its faces.
constexpr double kBlVisible = 0.05;
// Offsets are divided by the least cosine, at most this many times.
constexpr double kBlMaxStretch = 4.0;
// A cell is valid when its least scaled corner Jacobian passes this.
constexpr double kBlMinJacobian = 1e-6;

[[noreturn]] void bl_fail(const std::string& rWhy) {
    throw std::invalid_argument(std::string(kBoundaryLayerPrefix) + rWhy);
}

BlV3 bl_sub(const BlV3& rA, const BlV3& rB) {
    return {rA[0] - rB[0], rA[1] - rB[1], rA[2] - rB[2]};
}
double bl_dot(const BlV3& rA, const BlV3& rB) {
    return rA[0] * rB[0] + rA[1] * rB[1] + rA[2] * rB[2];
}
BlV3 bl_cross(const BlV3& rA, const BlV3& rB) {
    return {rA[1] * rB[2] - rA[2] * rB[1], rA[2] * rB[0] - rA[0] * rB[2],
            rA[0] * rB[1] - rA[1] * rB[0]};
}
double bl_norm(const BlV3& rA) {
    return std::sqrt(bl_dot(rA, rA));
}
BlV3 bl_unit(const BlV3& rA) {
    const double l = bl_norm(rA);
    return l > 0.0 ? BlV3{rA[0] / l, rA[1] / l, rA[2] / l} : BlV3{0.0, 0.0, 0.0};
}

// Row I of a C-contiguous array, as bytes.
std::size_t bl_row_bytes(const NDArray& rA) {
    std::size_t r = dtype_size(rA.Dtype());
    for (std::size_t d = 1; d < rA.Ndim(); ++d)
        r *= rA.Shape()[d];
    return r;
}

// The rows `rFrom` of `rA`, as a new array.
NDArray bl_gather(const NDArray& rA, const std::vector<std::int64_t>& rFrom) {
    std::vector<std::size_t> shape = rA.Shape();
    if (shape.empty())
        shape.push_back(0);
    shape[0] = rFrom.size();
    NDArray out = NDArray::Uninit(rA.Dtype(), shape);
    const std::size_t row = bl_row_bytes(rA);
    for (std::size_t i = 0; i < rFrom.size(); ++i)
        std::memcpy(out.Data() + i * row, rA.Data() + static_cast<std::size_t>(rFrom[i]) * row,
                    row);
    return out;
}

// The front test (6.9.1, after Pirzadeh 1996's advancing-layer test): every
// layer cell against the surface's faces and the last layer's (quads by
// their first diagonal), over a bounding-volume hierarchy, by the exact
// triangle-triangle kernel: a face that crosses one of the cell's faces, or
// lies inside the cell, meets it. A face whose column shares a surface point
// with the cell is its neighbour and is left out. Marks the surface points
// of each such cell, and of a last-layer face it meets, in `rBad`; returns
// the pairs found. Cells are tested in parallel, each pair once, so the
// count does not depend on threads.
std::int64_t bl_front_crossings(const detail::LayerSweepResult& rSweep, std::vector<char>& rBad) {
    detail::TriangleSoup soup;
    const std::size_t npts = rSweep.mOrigin.size();
    soup.mPoints.resize(npts);
    for (std::size_t i = 0; i < npts; ++i)
        soup.mPoints[i] = {rSweep.mPoints[3 * i], rSweep.mPoints[3 * i + 1],
                           rSweep.mPoints[3 * i + 2]};
    std::vector<char> is_cap;
    const auto add = [&](std::int64_t A, std::int64_t B, std::int64_t C, bool Cap) {
        soup.mVertices.push_back({A, B, C});
        for (const std::int64_t v : {A, B, C})
            soup.mCorners.push_back(soup.mPoints[static_cast<std::size_t>(v)]);
        soup.mSourceCell.push_back(static_cast<std::int64_t>(is_cap.size()));
        is_cap.push_back(Cap ? 1 : 0);
    };
    for (const bool cap : {true, false}) {
        for (const auto& r_c : cap ? rSweep.mCapTriangles : rSweep.mWallTriangles)
            add(r_c.mNodes[0], r_c.mNodes[1], r_c.mNodes[2], cap);
        for (const auto& r_c : cap ? rSweep.mCapQuads : rSweep.mWallQuads) {
            add(r_c.mNodes[0], r_c.mNodes[1], r_c.mNodes[2], cap);
            add(r_c.mNodes[0], r_c.mNodes[2], r_c.mNodes[3], cap);
        }
    }
    const std::size_t nt = is_cap.size();
    if (nt == 0)
        return 0;
    const detail::TriBvh tree = detail::build_tri_bvh(soup);
    // The cells, each with its kind (its node count), and the triangles of
    // its faces by local node.
    struct BlCell {
        const detail::LayerSweepCell* pCell;
        int mKind;
    };
    std::vector<BlCell> cells;
    for (const auto& r_c : rSweep.mHexahedra)
        cells.push_back({&r_c, 8});
    for (const auto& r_c : rSweep.mWedges)
        cells.push_back({&r_c, 6});
    for (const auto& r_c : rSweep.mPyramids)
        cells.push_back({&r_c, 5});
    for (const auto& r_c : rSweep.mTetrahedra)
        cells.push_back({&r_c, 4});
    using Tri = std::array<int, 3>;
    static const std::vector<Tri> k_hex = {{0, 1, 2}, {0, 2, 3}, {4, 5, 6}, {4, 6, 7},
                                           {0, 1, 5}, {0, 5, 4}, {1, 2, 6}, {1, 6, 5},
                                           {2, 3, 7}, {2, 7, 6}, {3, 0, 4}, {3, 4, 7}};
    static const std::vector<Tri> k_wedge = {{0, 1, 2}, {3, 4, 5}, {0, 1, 4}, {0, 4, 3},
                                             {1, 2, 5}, {1, 5, 4}, {2, 0, 3}, {2, 3, 5}};
    static const std::vector<Tri> k_pyramid = {{0, 1, 2}, {0, 2, 3}, {0, 1, 4},
                                               {1, 2, 4}, {2, 3, 4}, {3, 0, 4}};
    static const std::vector<Tri> k_tet = {{0, 1, 2}, {0, 1, 3}, {1, 2, 3}, {2, 0, 3}};
    const auto faces_of = [&](int Kind) -> const std::vector<Tri>& {
        return Kind == 8 ? k_hex : Kind == 6 ? k_wedge : Kind == 5 ? k_pyramid : k_tet;
    };
    const auto pt = [&](std::int64_t I) { return rSweep.mPoints.data() + 3 * I; };
    std::vector<std::vector<std::size_t>> hits(cells.size());
    parallel_for(
        cells.size(),
        [&](std::size_t C) {
            const detail::LayerSweepCell& r_c = *cells[C].pCell;
            const int kind = cells[C].mKind;
            const auto nn = static_cast<std::size_t>(kind);
            BlV3 lo{std::numeric_limits<double>::infinity(),
                    std::numeric_limits<double>::infinity(),
                    std::numeric_limits<double>::infinity()};
            BlV3 hi{-lo[0], -lo[1], -lo[2]};
            BlV3 mid{0.0, 0.0, 0.0};
            for (std::size_t k = 0; k < nn; ++k) {
                const double* x = pt(r_c.mNodes[k]);
                for (std::size_t d = 0; d < 3; ++d) {
                    lo[d] = std::min(lo[d], x[d]);
                    hi[d] = std::max(hi[d], x[d]);
                    mid[d] += x[d] / double(nn);
                }
            }
            const auto shares = [&](std::size_t U) {
                for (std::size_t a = 0; a < nn; ++a)
                    for (std::size_t b = 0; b < 3; ++b)
                        if (rSweep.mOrigin[static_cast<std::size_t>(r_c.mNodes[a])] ==
                            rSweep.mOrigin[static_cast<std::size_t>(soup.mVertices[U][b])])
                            return true;
                return false;
            };
            // Inside the cell: on no face's far side from its centroid (the
            // cells are valid, so nearly convex).
            const auto inside = [&](const double* pP) {
                for (const Tri& r_f : faces_of(kind)) {
                    const double* a = pt(r_c.mNodes[static_cast<std::size_t>(r_f[0])]);
                    const double* b = pt(r_c.mNodes[static_cast<std::size_t>(r_f[1])]);
                    const double* c = pt(r_c.mNodes[static_cast<std::size_t>(r_f[2])]);
                    const double s_mid = detail::orient3d(a, b, c, mid.data());
                    const double s_p = detail::orient3d(a, b, c, pP);
                    if ((s_mid > 0.0 && s_p < 0.0) || (s_mid < 0.0 && s_p > 0.0))
                        return false;
                }
                return true;
            };
            std::vector<std::int32_t> stack{0};
            while (!stack.empty()) {
                const auto n = static_cast<std::size_t>(stack.back());
                stack.pop_back();
                const detail::TriBvhNode& r_n = tree.mNodes[n];
                bool apart = false;
                for (std::size_t d = 0; d < 3 && !apart; ++d)
                    apart = r_n.mHi[d] < lo[d] || r_n.mLo[d] > hi[d];
                if (apart)
                    continue;
                if (r_n.mRight >= 0) {
                    stack.push_back(r_n.mRight);
                    stack.push_back(static_cast<std::int32_t>(n + 1));
                    continue;
                }
                for (std::int32_t i = 0; i < r_n.mCount; ++i) {
                    const auto u = static_cast<std::size_t>(
                        tree.mOrder[static_cast<std::size_t>(r_n.mBegin + i)]);
                    if (shares(u))
                        continue;
                    const double* d0 = soup.mCorners[3 * u].data();
                    const double* d1 = soup.mCorners[3 * u + 1].data();
                    const double* d2 = soup.mCorners[3 * u + 2].data();
                    bool meets = inside(d0);
                    for (const Tri& r_f : faces_of(kind)) {
                        if (meets)
                            break;
                        meets = detail::tri_tri_intersect(
                            pt(r_c.mNodes[static_cast<std::size_t>(r_f[0])]),
                            pt(r_c.mNodes[static_cast<std::size_t>(r_f[1])]),
                            pt(r_c.mNodes[static_cast<std::size_t>(r_f[2])]), d0, d1, d2);
                    }
                    if (meets)
                        hits[C].push_back(u);
                }
            }
        },
        16);
    std::int64_t count = 0;
    for (std::size_t c = 0; c < cells.size(); ++c)
        for (const std::size_t u : hits[c]) {
            ++count;
            const detail::LayerSweepCell& r_c = *cells[c].pCell;
            for (int k = 0; k < cells[c].mKind; ++k)
                rBad[static_cast<std::size_t>(rSweep.mOrigin[static_cast<std::size_t>(
                    r_c.mNodes[static_cast<std::size_t>(k)])])] = 1;
            if (is_cap[u])
                for (std::size_t k = 0; k < 3; ++k)
                    rBad[static_cast<std::size_t>(
                        rSweep.mOrigin[static_cast<std::size_t>(soup.mVertices[u][k])])] = 1;
        }
    return count;
}

}  // namespace

BoundaryLayerResult boundary_layer(const Mesh& rSurface, const BoundaryLayerOptions& rOptions) {
    const BoundaryLayerOptions& o = rOptions;
    if (!(std::isfinite(o.mFirstHeight) && o.mFirstHeight > 0.0))
        bl_fail("first_height must be positive and finite");
    if (!(std::isfinite(o.mGrowth) && o.mGrowth >= 1.0 && o.mGrowth <= 10.0))
        bl_fail("growth must be between 1 and 10");
    if (o.mNumLayers < 1 || o.mNumLayers > 1000)
        bl_fail("num_layers must be between 1 and 1000");
    if (!(o.mCollisionFraction > 0.0 && o.mCollisionFraction <= 1.0))
        bl_fail("collision_fraction must be in (0, 1]");
    if (o.mSmoothIterations < 0 || o.mSmoothIterations > 1000)
        bl_fail("smooth_iterations must be between 0 and 1000");
    if (o.mMaxShrink < 0 || o.mMaxShrink > 60)
        bl_fail("max_shrink must be between 0 and 60");
    if (o.mCollapse != BoundaryLayerCollapse::Points &&
        o.mCollapse != BoundaryLayerCollapse::Layers)
        bl_fail("collapse must be points or layers");
    if (o.mNormalSmoothIterations < 0 || o.mNormalSmoothIterations > 1000)
        bl_fail("normal_smooth_iterations must be between 0 and 1000");
    const std::size_t np = rSurface.NumPoints();
    const NDArray& r_pts = rSurface.Points();
    if (np > 0 && (r_pts.Ndim() != 2 || r_pts.Shape()[1] != 3))
        bl_fail("the surface's points must be three-dimensional");
    std::vector<double> x(3 * np);
    for (std::size_t i = 0; i < 3 * np; ++i)
        x[i] = detail::read_double(r_pts, i);
    const auto point = [&](std::int64_t P) {
        const auto i = static_cast<std::size_t>(P);
        return BlV3{x[3 * i], x[3 * i + 1], x[3 * i + 2]};
    };

    // The faces, and where each comes from.
    std::vector<std::array<std::int64_t, 4>> faces;
    std::vector<std::array<std::size_t, 2>> source;  // (block, cell)
    std::vector<char> face_block(rSurface.NumCellBlocks(), 0);
    for (std::size_t b = 0; b < rSurface.NumCellBlocks(); ++b) {
        const auto cells = rSurface.Cells(b);
        const std::string type(cells.Type());
        if (type.rfind("line", 0) == 0 || type == "vertex")
            continue;
        if ((type != "triangle" && type != "quad") || cells.IsRagged())
            bl_fail("a '" + type + "' block: layers grow from triangles and quads (linearize or " +
                    "extract_surface first)");
        face_block[b] = 1;
        const std::size_t k = type == "triangle" ? 3 : 4;
        for (std::size_t c = 0; c < cells.NumCells(); ++c) {
            std::array<std::int64_t, 4> f{-1, -1, -1, -1};
            for (std::size_t i = 0; i < k; ++i) {
                f[i] = detail::read_int(cells.Conn(), k * c + i);
                if (f[i] < 0 || static_cast<std::size_t>(f[i]) >= np)
                    bl_fail("a face names a point out of range");
            }
            faces.push_back(f);
            source.push_back({b, c});
        }
    }
    if (faces.empty())
        bl_fail("the surface has no triangle or quad");
    const std::size_t nf = faces.size();
    const auto corners = [&](std::size_t F) { return faces[F][3] >= 0 ? 4 : 3; };
    // Manifold and consistently oriented.
    {
        std::map<std::pair<std::int64_t, std::int64_t>, std::array<int, 2>> uses;
        for (std::size_t f = 0; f < nf; ++f)
            for (int i = 0; i < corners(f); ++i) {
                const std::int64_t a = faces[f][static_cast<std::size_t>(i)];
                const std::int64_t b = faces[f][static_cast<std::size_t>((i + 1) % corners(f))];
                auto& r_u = uses[{std::min(a, b), std::max(a, b)}];
                ++r_u[0];
                r_u[1] += a < b ? 1 : -1;
            }
        for (const auto& [e, r_u] : uses) {
            if (r_u[0] > 2)
                bl_fail("edge (" + std::to_string(e.first) + ", " + std::to_string(e.second) +
                        ") is used by more than two faces");
            if (r_u[0] == 2 && r_u[1] != 0)
                bl_fail("the faces at edge (" + std::to_string(e.first) + ", " +
                        std::to_string(e.second) + ") disagree on its direction (orient them)");
        }
    }
    // Unit face normals along the growth, and each point's faces with its
    // corner angle there.
    const double sign = o.mReverse ? -1.0 : 1.0;
    std::vector<BlV3> fn(nf);
    std::vector<std::vector<std::pair<std::size_t, double>>> around(np);
    std::vector<std::vector<std::int64_t>> nbrs(np);
    for (std::size_t f = 0; f < nf; ++f) {
        const int k = corners(f);
        BlV3 n{0.0, 0.0, 0.0};
        for (int i = 0; i < k; ++i) {
            const BlV3 a = point(faces[f][static_cast<std::size_t>(i)]);
            const BlV3 b = point(faces[f][static_cast<std::size_t>((i + 1) % k)]);
            n = {n[0] + (a[1] - b[1]) * (a[2] + b[2]), n[1] + (a[2] - b[2]) * (a[0] + b[0]),
                 n[2] + (a[0] - b[0]) * (a[1] + b[1])};
        }
        n = bl_unit(n);
        fn[f] = {sign * n[0], sign * n[1], sign * n[2]};
        for (int i = 0; i < k; ++i) {
            const std::int64_t p = faces[f][static_cast<std::size_t>(i)];
            const std::int64_t q = faces[f][static_cast<std::size_t>((i + 1) % k)];
            const std::int64_t r = faces[f][static_cast<std::size_t>((i + k - 1) % k)];
            const BlV3 u = bl_unit(bl_sub(point(q), point(p)));
            const BlV3 v = bl_unit(bl_sub(point(r), point(p)));
            const double angle = std::acos(std::clamp(bl_dot(u, v), -1.0, 1.0));
            around[static_cast<std::size_t>(p)].push_back({f, angle});
            nbrs[static_cast<std::size_t>(p)].push_back(q);
            nbrs[static_cast<std::size_t>(p)].push_back(r);
        }
    }
    for (auto& r_n : nbrs) {
        std::sort(r_n.begin(), r_n.end());
        r_n.erase(std::unique(r_n.begin(), r_n.end()), r_n.end());
    }
    // Each point's most normal normal (Aubry and Löhner 2009): from the
    // angle-weighted mean, steps toward the face that sees it least.
    std::vector<BlV3> dir(np, BlV3{0.0, 0.0, 0.0});
    std::vector<double> stretch(np, 0.0), sigma(np, 0.0), seen(np, 0.0);
    std::vector<char> used(np, 0), not_visible(np, 0);
    parallel_for(
        np,
        [&](std::size_t P) {
            const auto& r_a = around[P];
            if (r_a.empty())
                return;
            used[P] = 1;
            BlV3 m{0.0, 0.0, 0.0};
            for (const auto& [f, w] : r_a)
                m = {m[0] + w * fn[f][0], m[1] + w * fn[f][1], m[2] + w * fn[f][2]};
            BlV3 n = bl_unit(m);
            const auto least = [&](const BlV3& rN, std::size_t& rWho) {
                double lo = std::numeric_limits<double>::infinity();
                for (const auto& [f, w] : r_a) {
                    (void)w;
                    const double d = bl_dot(rN, fn[f]);
                    if (d < lo) {
                        lo = d;
                        rWho = f;
                    }
                }
                return lo;
            };
            std::size_t who = 0;
            BlV3 best = n;
            double best_lo = least(n, who);
            for (int it = 0; it < 100 && r_a.size() > 1; ++it) {
                const double eta = 1.0 / double(it + 2);
                n = bl_unit(
                    {n[0] + eta * fn[who][0], n[1] + eta * fn[who][1], n[2] + eta * fn[who][2]});
                const double lo = least(n, who);
                if (lo > best_lo) {
                    best_lo = lo;
                    best = n;
                }
            }
            dir[P] = best;
            if (!(best_lo > kBlVisible)) {
                not_visible[P] = 1;
                return;
            }
            stretch[P] = 1.0 / std::max(best_lo, 1.0 / kBlMaxStretch);
            sigma[P] = 1.0;
            seen[P] = best_lo;
        },
        64);
    // Smoothed directions (6.9.1): each moved toward its neighbours' mean,
    // kept while every face around it still sees it at more than half the
    // least cosine its most normal normal had (Jacobi passes, so the result
    // does not depend on the order).
    const auto least_cos = [&](std::size_t P, const BlV3& rN) {
        double lo = std::numeric_limits<double>::infinity();
        for (const auto& [f, w] : around[P]) {
            (void)w;
            lo = std::min(lo, bl_dot(rN, fn[f]));
        }
        return lo;
    };
    for (int it = 0; it < o.mNormalSmoothIterations; ++it) {
        std::vector<BlV3> next = dir;
        parallel_for(
            np,
            [&](std::size_t P) {
                if (!(sigma[P] > 0.0) || nbrs[P].empty())
                    return;
                BlV3 m = dir[P];
                for (const std::int64_t q : nbrs[P]) {
                    const auto& r_q = dir[static_cast<std::size_t>(q)];
                    if (sigma[static_cast<std::size_t>(q)] > 0.0)
                        m = {m[0] + r_q[0], m[1] + r_q[1], m[2] + r_q[2]};
                }
                const BlV3 n = bl_unit(m);
                const double lo = least_cos(P, n);
                if (lo > 0.5 * seen[P] && lo > kBlVisible)
                    next[P] = n;
            },
            64);
        dir = std::move(next);
    }
    if (o.mNormalSmoothIterations > 0)
        for (std::size_t p = 0; p < np; ++p)
            if (sigma[p] > 0.0)
                stretch[p] = 1.0 / std::max(least_cos(p, dir[p]), 1.0 / kBlMaxStretch);
    // The levels' offsets, before each point's own fraction and stretch.
    const auto nl = static_cast<std::size_t>(o.mNumLayers);
    std::vector<double> level(nl);
    {
        double h = o.mFirstHeight, s = 0.0;
        for (std::size_t k = 0; k < nl; ++k) {
            s += h;
            level[k] = s;
            h *= o.mGrowth;
        }
    }
    const double total = level.back();
    // Collisions: a ray along each point's direction against the surface.
    detail::TriangleSoup soup;
    soup.mPoints.resize(np);
    for (std::size_t i = 0; i < np; ++i)
        soup.mPoints[i] = point(static_cast<std::int64_t>(i));
    for (std::size_t f = 0; f < nf; ++f) {
        const auto& r_f = faces[f];
        const auto tri = [&](std::int64_t A, std::int64_t B, std::int64_t C) {
            soup.mVertices.push_back({A, B, C});
            soup.mCorners.push_back(point(A));
            soup.mCorners.push_back(point(B));
            soup.mCorners.push_back(point(C));
            soup.mSourceCell.push_back(static_cast<std::int64_t>(f));
        };
        tri(r_f[0], r_f[1], r_f[2]);
        if (r_f[3] >= 0)
            tri(r_f[0], r_f[2], r_f[3]);
    }
    const detail::TriBvh tree = detail::build_tri_bvh(soup);
    parallel_for(
        np,
        [&](std::size_t P) {
            if (!(sigma[P] > 0.0))
                return;
            double shortest = std::numeric_limits<double>::infinity();
            for (const std::int64_t q : nbrs[P])
                shortest = std::min(shortest,
                                    bl_norm(bl_sub(point(q), point(static_cast<std::int64_t>(P)))));
            const double want = total * stretch[P];
            const detail::TriRayHit hit =
                detail::tri_bvh_ray(tree, point(static_cast<std::int64_t>(P)), dir[P],
                                    1e-6 * shortest, want / o.mCollisionFraction);
            if (hit.mTriangle >= 0)
                sigma[P] = std::min(1.0, o.mCollisionFraction * hit.mT / want);
        },
        64);
    // Tapered: each fraction at most the mean of itself and its neighbours'.
    for (int it = 0; it < o.mSmoothIterations; ++it) {
        std::vector<double> next = sigma;
        for (std::size_t p = 0; p < np; ++p) {
            if (!used[p] || nbrs[p].empty())
                continue;
            double mean = 0.0;
            for (const std::int64_t q : nbrs[p])
                mean += sigma[static_cast<std::size_t>(q)];
            mean /= double(nbrs[p].size());
            next[p] = std::min(sigma[p], 0.5 * (sigma[p] + mean));
        }
        sigma = std::move(next);
    }
    // Sweep, judge, shrink the points of the cells that fold, and again.
    detail::LayerSweepInput in;
    in.pPoints = x.data();
    in.mNumPoints = static_cast<std::int64_t>(np);
    // The sweep orients a cell by its face's normal: against it, each face
    // runs the other way.
    in.mFaces = faces;
    if (o.mReverse)
        for (auto& r_f : in.mFaces) {
            if (r_f[3] >= 0)
                std::swap(r_f[1], r_f[3]);
            else
                std::swap(r_f[1], r_f[2]);
        }
    in.mLayers = o.mNumLayers;
    in.mDirections.resize(3 * np);
    for (std::size_t p = 0; p < np; ++p)
        for (std::size_t c = 0; c < 3; ++c)
            in.mDirections[3 * p + c] = dir[p][c];
    std::vector<int> halvings(np, 0);
    // With Layers, the layers each point grows (6.9.1).
    std::vector<std::size_t> grown(np, nl);
    const bool by_layers = o.mCollapse == BoundaryLayerCollapse::Layers;
    detail::LayerSweepResult sweep;
    double min_sj = std::numeric_limits<double>::infinity();
    std::int64_t front_crossings = 0;
    const int max_rounds = std::max(64 * (o.mMaxShrink + 2), by_layers ? 4 * o.mNumLayers + 64 : 0);
    for (int round = 0;; ++round) {
        in.mOffsets.assign(np * nl, 0.0);
        for (std::size_t p = 0; p < np; ++p)
            if (sigma[p] > 0.0 && grown[p] > 0)
                for (std::size_t k = 0; k < nl; ++k)
                    in.mOffsets[p * nl + k] =
                        sigma[p] * stretch[p] * level[std::min(k, grown[p] - 1)];
        sweep = detail::layer_sweep(in);
        std::vector<char> bad_point(np, 0);
        bool bad = false;
        min_sj = std::numeric_limits<double>::infinity();
        const auto judge = [&](const std::vector<detail::LayerSweepCell>& rCells, int Kind) {
            for (const auto& r_c : rCells) {
                const double sj =
                    detail::layer_cell_scaled_jacobian(sweep.mPoints.data(), r_c, Kind);
                min_sj = std::min(min_sj, sj);
                if (sj > kBlMinJacobian)
                    continue;
                bad = true;
                for (int i = 0; i < Kind; ++i)
                    bad_point[static_cast<std::size_t>(sweep.mOrigin[static_cast<std::size_t>(
                        r_c.mNodes[static_cast<std::size_t>(i)])])] = 1;
            }
        };
        judge(sweep.mHexahedra, 8);
        judge(sweep.mWedges, 6);
        judge(sweep.mPyramids, 5);
        judge(sweep.mTetrahedra, 4);
        if (!bad && o.mFrontCheck) {
            const std::int64_t crossings = bl_front_crossings(sweep, bad_point);
            front_crossings += crossings;
            bad = crossings > 0;
        }
        if (!bad)
            break;
        if (round > max_rounds)
            throw std::logic_error(std::string(kBoundaryLayerPrefix) +
                                   "internal error: the layers do not settle");
        for (std::size_t p = 0; p < np; ++p) {
            if (!bad_point[p] || !(sigma[p] > 0.0) || grown[p] == 0)
                continue;
            if (by_layers) {
                if (--grown[p] == 0)
                    sigma[p] = 0.0;
            } else if (halvings[p] < o.mMaxShrink) {
                sigma[p] *= 0.5;
                ++halvings[p];
            } else {
                sigma[p] = 0.0;
            }
        }
    }

    // The mesh: the cells, then the wall and the cap.
    BoundaryLayerResult result;
    const auto warn = [&](const std::string& rMsg) {
        log::warn("boundary_layer: {}", rMsg);
        result.mWarnings.push_back(rMsg);
    };
    Mesh& r_m = result.mMesh;
    const std::size_t nout = sweep.mOrigin.size();
    NDArray pts = NDArray::Uninit(DType::Float64, {nout, std::size_t{3}});
    std::copy(sweep.mPoints.begin(), sweep.mPoints.end(), pts.As<double>());
    r_m.AssignPoints(std::move(pts));
    struct Block {
        const char* pType;
        const std::vector<detail::LayerSweepCell>* pCells;
        std::size_t mNodes;
        int mRole;  // 3 cells, 0 wall, 1 cap
    };
    const std::vector<Block> blocks = {
        {"wedge", &sweep.mWedges, 6, 3},           {"hexahedron", &sweep.mHexahedra, 8, 3},
        {"pyramid", &sweep.mPyramids, 5, 3},       {"tetra", &sweep.mTetrahedra, 4, 3},
        {"triangle", &sweep.mWallTriangles, 3, 0}, {"quad", &sweep.mWallQuads, 4, 0},
        {"triangle", &sweep.mCapTriangles, 3, 1},  {"quad", &sweep.mCapQuads, 4, 1}};
    std::vector<const Block*> written;
    std::vector<std::int64_t> wall_ids, cap_ids, cell_ids;
    std::int64_t global = 0;
    for (const Block& r_b : blocks) {
        if (r_b.pCells->empty())
            continue;
        NDArray conn = NDArray::Uninit(DType::Int64, {r_b.pCells->size(), r_b.mNodes});
        std::int64_t* pc = conn.As<std::int64_t>();
        for (std::size_t c = 0; c < r_b.pCells->size(); ++c) {
            for (std::size_t k = 0; k < r_b.mNodes; ++k)
                pc[r_b.mNodes * c + k] = (*r_b.pCells)[c].mNodes[k];
            (r_b.mRole == 3 ? cell_ids : r_b.mRole == 0 ? wall_ids : cap_ids).push_back(global++);
        }
        r_m.AddCellBlock(r_b.pType, std::move(conn));
        written.push_back(&r_b);
    }
    // The layer of every cell and face.
    {
        std::vector<NDArray> arrays;
        for (const Block* p_b : written) {
            NDArray a = NDArray::Uninit(DType::Int32, {p_b->pCells->size()});
            for (std::size_t c = 0; c < p_b->pCells->size(); ++c)
                a.As<std::int32_t>()[c] = (*p_b->pCells)[c].mLayer;
            arrays.push_back(std::move(a));
        }
        r_m.AddCellData(kBoundaryLayerName, std::move(arrays));
    }
    // The surface's own cell data, carried from each face to what grew from it.
    for (const std::string& name : rSurface.CellDataNames()) {
        if (name == kBoundaryLayerName)
            continue;
        const NDArray* p_first = nullptr;
        bool ok = true;
        for (std::size_t b = 0; b < rSurface.NumCellBlocks() && ok; ++b) {
            if (!face_block[b])
                continue;
            const NDArray& r_a = rSurface.CellData(name, b);
            if (!p_first)
                p_first = &r_a;
            else if (r_a.Dtype() != p_first->Dtype() ||
                     bl_row_bytes(r_a) != bl_row_bytes(*p_first) || r_a.Ndim() != p_first->Ndim())
                ok = false;
        }
        if (!p_first || !ok) {
            warn("cell data '" + name + "' differs between the surface's blocks and is dropped");
            continue;
        }
        // One array of every face's rows, then gathered per written block.
        std::vector<std::size_t> shape = p_first->Shape();
        shape[0] = nf;
        NDArray all = NDArray::Uninit(p_first->Dtype(), shape);
        const std::size_t row = bl_row_bytes(*p_first);
        for (std::size_t f = 0; f < nf; ++f)
            std::memcpy(all.Data() + f * row,
                        rSurface.CellData(name, source[f][0]).Data() + source[f][1] * row, row);
        std::vector<NDArray> arrays;
        for (const Block* p_b : written) {
            std::vector<std::int64_t> from(p_b->pCells->size());
            for (std::size_t c = 0; c < from.size(); ++c)
                from[c] = (*p_b->pCells)[c].mFace;
            arrays.push_back(bl_gather(all, from));
        }
        r_m.AddCellData(name, std::move(arrays));
    }
    for (const std::string& name : rSurface.PointDataNames())
        r_m.AddPointData(name, bl_gather(rSurface.PointData(name), sweep.mOrigin));
    for (const std::string& name : rSurface.FieldDataNames())
        r_m.AddFieldData(name, detail::data_owned_copy(rSurface.FieldData(name)));
    const auto column = [](const std::vector<std::int64_t>& rV) {
        NDArray a = NDArray::Uninit(DType::Int64, {rV.size()});
        std::copy(rV.begin(), rV.end(), a.As<std::int64_t>());
        return a;
    };
    // The surface's regions: a Cell region holds what grew from its faces
    // (cells, wall and cap), a Point region every copy of its points.
    {
        std::vector<std::int64_t> block_start(rSurface.NumCellBlocks() + 1, 0);
        for (std::size_t b = 0; b < rSurface.NumCellBlocks(); ++b)
            block_start[b + 1] =
                block_start[b] + static_cast<std::int64_t>(rSurface.Cells(b).NumCells());
        std::vector<std::int64_t> face_of_global(static_cast<std::size_t>(block_start.back()), -1);
        for (std::size_t f = 0; f < nf; ++f)
            face_of_global[static_cast<std::size_t>(block_start[source[f][0]]) + source[f][1]] =
                static_cast<std::int64_t>(f);
        std::size_t dropped = 0;
        for (std::size_t i = 0; i < rSurface.NumRegions(); ++i) {
            const Region& r_r = rSurface.Region(i);
            const NDArray& r_e = r_r.mEntries;
            std::vector<std::int64_t> entries;
            if (r_r.mKind == RegionKind::Cell) {
                std::vector<char> in_region(nf, 0);
                for (std::size_t k = 0; k < r_r.NumEntries(); ++k) {
                    const std::int64_t g = detail::read_int(r_e, k);
                    if (g >= 0 && static_cast<std::size_t>(g) < face_of_global.size() &&
                        face_of_global[static_cast<std::size_t>(g)] >= 0)
                        in_region[static_cast<std::size_t>(
                            face_of_global[static_cast<std::size_t>(g)])] = 1;
                }
                std::int64_t id = 0;
                for (const Block* p_b : written)
                    for (const auto& r_c : *p_b->pCells) {
                        if (in_region[static_cast<std::size_t>(r_c.mFace)])
                            entries.push_back(id);
                        ++id;
                    }
            } else if (r_r.mKind == RegionKind::Point) {
                std::vector<char> in_region(np, 0);
                for (std::size_t k = 0; k < r_r.NumEntries(); ++k) {
                    const std::int64_t q = detail::read_int(r_e, k);
                    if (q >= 0 && static_cast<std::size_t>(q) < np)
                        in_region[static_cast<std::size_t>(q)] = 1;
                }
                for (std::size_t q = 0; q < nout; ++q)
                    if (in_region[static_cast<std::size_t>(sweep.mOrigin[q])])
                        entries.push_back(static_cast<std::int64_t>(q));
            } else {
                ++dropped;
                continue;
            }
            r_m.AddRegion(Region(r_r.mName, r_r.mKind,
                                 r_r.mKind == RegionKind::Cell ? -1 : r_r.mDim, r_r.mTag,
                                 column(entries)));
        }
        if (dropped > 0)
            warn(std::to_string(dropped) + " side region(s) of the surface are dropped");
    }
    r_m.AddRegion(Region(kBoundaryLayerCells, RegionKind::Cell, 3, -1, column(cell_ids)));
    r_m.AddRegion(Region(kBoundaryLayerWall, RegionKind::Cell, 2, -1, column(wall_ids)));
    r_m.AddRegion(Region(kBoundaryLayerCap, RegionKind::Cell, 2, -1, column(cap_ids)));

    // The report.
    result.mNumWedges = static_cast<std::int64_t>(sweep.mWedges.size());
    result.mNumHexahedra = static_cast<std::int64_t>(sweep.mHexahedra.size());
    result.mNumPyramids = static_cast<std::int64_t>(sweep.mPyramids.size());
    result.mNumTetrahedra = static_cast<std::int64_t>(sweep.mTetrahedra.size());
    result.mNumCells =
        result.mNumWedges + result.mNumHexahedra + result.mNumPyramids + result.mNumTetrahedra;
    result.mMinScaledJacobian = result.mNumCells > 0 ? min_sj : 0.0;
    result.mMinThickness = std::numeric_limits<double>::infinity();
    result.mMaxThickness = 0.0;
    for (std::size_t p = 0; p < np; ++p) {
        if (!used[p])
            continue;
        const bool fewer = grown[p] < nl && grown[p] > 0 && sigma[p] > 0.0;
        result.mNumShrunk += sigma[p] < 1.0 || grown[p] < nl ? 1 : 0;
        result.mNumCollapsed += sigma[p] > 0.0 && grown[p] > 0 ? 0 : 1;
        result.mNumNotVisible += not_visible[p] ? 1 : 0;
        result.mNumFewerLayers += fewer ? 1 : 0;
        const double thick = grown[p] > 0 ? sigma[p] * level[grown[p] - 1] : 0.0;
        result.mMinThickness = std::min(result.mMinThickness, thick);
        result.mMaxThickness = std::max(result.mMaxThickness, thick);
    }
    result.mNumFrontCrossings = front_crossings;
    if (!std::isfinite(result.mMinThickness))
        result.mMinThickness = 0.0;
    {
        std::map<std::pair<std::int64_t, std::int64_t>, int> uses;
        for (const auto* p_cap : {&sweep.mCapTriangles, &sweep.mCapQuads})
            for (const auto& r_c : *p_cap) {
                const int k = r_c.mNodes[3] >= 0 ? 4 : 3;
                for (int i = 0; i < k; ++i) {
                    const std::int64_t a = r_c.mNodes[static_cast<std::size_t>(i)];
                    const std::int64_t b = r_c.mNodes[static_cast<std::size_t>((i + 1) % k)];
                    ++uses[{std::min(a, b), std::max(a, b)}];
                }
            }
        result.mCapClosed = !uses.empty();
        for (const auto& [e, n] : uses) {
            (void)e;
            if (n != 2)
                result.mCapClosed = false;
        }
    }
    if (result.mNumFrontCrossings > 0)
        log::debug("boundary_layer: {} crossings of the last layer's faces given way to",
                   result.mNumFrontCrossings);
    if (result.mNumFewerLayers > 0)
        warn(std::to_string(result.mNumFewerLayers) +
             " points grow fewer layers than asked (their top layers would fold or cross)");
    if (result.mNumCollapsed > 0)
        warn(std::to_string(result.mNumCollapsed) +
             " points do not grow (their layers would fold or no direction sees all their "
             "faces); their cells are pyramids and tetrahedra");
    return result;
}

}  // namespace meshioplusplus

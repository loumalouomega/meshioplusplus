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

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/subdivide_surface.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/region_remap.hpp"
#include "meshioplusplus/mesh_api.hpp"

// Core-private includes
#include "../detail/crease_edges.hpp"
#include "../detail/surface_rings.hpp"

namespace meshioplusplus {

namespace {

constexpr const char* kSubdPrefix = "meshio++: subdivide_surface: ";
constexpr double kSubdPi = 3.141592653589793238462643383279502884;

std::uint64_t subd_edge_key(std::int64_t a, std::int64_t b) {
    const std::uint64_t lo = static_cast<std::uint64_t>(std::min(a, b));
    const std::uint64_t hi = static_cast<std::uint64_t>(std::max(a, b));
    return (lo << 32) | hi;
}

std::string subd_lower(std::string s) {
    for (char& c : s)
        if (c >= 'A' && c <= 'Z')
            c = static_cast<char>(c - 'A' + 'a');
    return s;
}

struct SubdEdgeInfo {
    std::int64_t mLo = 0;
    std::int64_t mHi = 0;
    std::vector<std::int64_t> mFaces;
    std::vector<std::int64_t> mOppVerts;  // for Loop: opposite vertex per incident triangle
    bool mIsCrease = false;
    bool mIsBoundary = false;
};

struct SubdFloatField {
    std::string mName;
    DType mOrigDtype = DType::Float64;
    std::vector<std::size_t> mTailShape;
    std::size_t mWidth = 1;
    std::vector<double> mData;  // n_points * mWidth
};

struct SubdOtherField {
    std::string mName;
    DType mDtype = DType::Int64;
    std::vector<std::size_t> mTailShape;
    std::size_t mRowBytes = 0;
    std::vector<std::uint8_t> mData;  // n_points * mRowBytes
};

std::vector<SubdEdgeInfo> subd_build_edges(const std::vector<std::vector<std::int64_t>>& rFaces,
                                           const std::unordered_set<std::uint64_t>& rCreases,
                                           bool PreserveBoundary,
                                           std::unordered_map<std::uint64_t, std::size_t>& rEdgeMap,
                                           std::int64_t& rNumCrease, std::int64_t& rNumBoundary) {
    std::vector<SubdEdgeInfo> edges;
    rEdgeMap.clear();
    rNumCrease = 0;
    rNumBoundary = 0;

    for (std::size_t f = 0; f < rFaces.size(); ++f) {
        const auto& ring = rFaces[f];
        const std::size_t m = ring.size();
        for (std::size_t k = 0; k < m; ++k) {
            const std::int64_t a = ring[k];
            const std::int64_t b = ring[(k + 1) % m];
            if (a == b)
                continue;
            const std::uint64_t key = subd_edge_key(a, b);
            auto it = rEdgeMap.find(key);
            std::size_t eid = 0;
            if (it == rEdgeMap.end()) {
                eid = edges.size();
                rEdgeMap.emplace(key, eid);
                SubdEdgeInfo e;
                e.mLo = std::min(a, b);
                e.mHi = std::max(a, b);
                edges.push_back(std::move(e));
            } else {
                eid = it->second;
            }
            edges[eid].mFaces.push_back(static_cast<std::int64_t>(f));
            if (m == 3)
                edges[eid].mOppVerts.push_back(ring[(k + 2) % 3]);
        }
    }

    // Sort edges deterministically by (mLo, mHi) so new edge-point ids are independent of hash
    // order.
    std::sort(edges.begin(), edges.end(), [](const SubdEdgeInfo& x, const SubdEdgeInfo& y) {
        return x.mLo != y.mLo ? x.mLo < y.mLo : x.mHi < y.mHi;
    });
    rEdgeMap.clear();
    for (std::size_t i = 0; i < edges.size(); ++i) {
        const std::uint64_t key = subd_edge_key(edges[i].mLo, edges[i].mHi);
        rEdgeMap.emplace(key, i);
        const bool is_bnd = edges[i].mFaces.size() == 1;
        const bool is_nm = edges[i].mFaces.size() >= 3;
        edges[i].mIsBoundary = is_bnd;
        edges[i].mIsCrease = rCreases.count(key) > 0 || is_nm || (PreserveBoundary && is_bnd);
        if (is_bnd)
            ++rNumBoundary;
        if (edges[i].mIsCrease)
            ++rNumCrease;
    }
    return edges;
}

void subd_step_loop(std::vector<double>& rPts, std::vector<std::vector<std::int64_t>>& rFaces,
                    std::vector<std::int64_t>& rParents,
                    std::unordered_set<std::uint64_t>& rCreases,
                    std::vector<SubdFloatField>& rFloatFields,
                    std::vector<SubdOtherField>& rOtherFields, bool PreserveBoundary,
                    std::int64_t& rNumCrease, std::int64_t& rNumBoundary) {
    const std::size_t n_old = rPts.size() / 3;
    std::unordered_map<std::uint64_t, std::size_t> edge_map;
    std::vector<SubdEdgeInfo> edges =
        subd_build_edges(rFaces, rCreases, PreserveBoundary, edge_map, rNumCrease, rNumBoundary);
    const std::size_t n_edges = edges.size();
    const std::size_t n_new = n_old + n_edges;

    // Incident neighbors & incident crease neighbors per old vertex.
    std::vector<std::vector<std::int64_t>> adj(n_old);
    std::vector<std::vector<std::int64_t>> crease_adj(n_old);
    std::vector<bool> non_manifold_vert(n_old, false);
    for (const auto& e : edges) {
        adj[static_cast<std::size_t>(e.mLo)].push_back(e.mHi);
        adj[static_cast<std::size_t>(e.mHi)].push_back(e.mLo);
        if (e.mFaces.size() >= 3) {
            non_manifold_vert[static_cast<std::size_t>(e.mLo)] = true;
            non_manifold_vert[static_cast<std::size_t>(e.mHi)] = true;
        }
        if (e.mIsCrease) {
            crease_adj[static_cast<std::size_t>(e.mLo)].push_back(e.mHi);
            crease_adj[static_cast<std::size_t>(e.mHi)].push_back(e.mLo);
        }
    }

    auto apply_stencil = [&](const std::vector<double>& rSrc, std::size_t width) {
        std::vector<double> dst(n_new * width, 0.0);
        // 1. Old vertices
        for (std::size_t v = 0; v < n_old; ++v) {
            const auto& cadj = crease_adj[v];
            if (non_manifold_vert[v] || cadj.size() == 1 || cadj.size() >= 3) {
                for (std::size_t c = 0; c < width; ++c)
                    dst[v * width + c] = rSrc[v * width + c];
            } else if (cadj.size() == 2) {
                const std::size_t a = static_cast<std::size_t>(cadj[0]);
                const std::size_t b = static_cast<std::size_t>(cadj[1]);
                for (std::size_t c = 0; c < width; ++c)
                    dst[v * width + c] = 0.75 * rSrc[v * width + c] +
                                         0.125 * (rSrc[a * width + c] + rSrc[b * width + c]);
            } else {
                const auto& nb = adj[v];
                const std::size_t k = nb.size();
                if (k == 0) {
                    for (std::size_t c = 0; c < width; ++c)
                        dst[v * width + c] = rSrc[v * width + c];
                } else {
                    const double cn =
                        0.375 + 0.25 * std::cos(2.0 * kSubdPi / static_cast<double>(k));
                    const double beta = (0.625 - cn * cn) / static_cast<double>(k);
                    const double self_w = 1.0 - static_cast<double>(k) * beta;
                    for (std::size_t c = 0; c < width; ++c) {
                        double sum = 0.0;
                        for (std::int64_t u : nb)
                            sum += rSrc[static_cast<std::size_t>(u) * width + c];
                        dst[v * width + c] = self_w * rSrc[v * width + c] + beta * sum;
                    }
                }
            }
        }
        // 2. Edge points
        for (std::size_t ei = 0; ei < n_edges; ++ei) {
            const auto& e = edges[ei];
            const std::size_t out_v = n_old + ei;
            const std::size_t u0 = static_cast<std::size_t>(e.mLo);
            const std::size_t u1 = static_cast<std::size_t>(e.mHi);
            if (e.mIsCrease || e.mOppVerts.size() != 2) {
                for (std::size_t c = 0; c < width; ++c)
                    dst[out_v * width + c] = 0.5 * (rSrc[u0 * width + c] + rSrc[u1 * width + c]);
            } else {
                const std::size_t w0 = static_cast<std::size_t>(e.mOppVerts[0]);
                const std::size_t w1 = static_cast<std::size_t>(e.mOppVerts[1]);
                for (std::size_t c = 0; c < width; ++c)
                    dst[out_v * width + c] = 0.375 * (rSrc[u0 * width + c] + rSrc[u1 * width + c]) +
                                             0.125 * (rSrc[w0 * width + c] + rSrc[w1 * width + c]);
            }
        }
        return dst;
    };

    rPts = apply_stencil(rPts, 3);
    for (auto& ff : rFloatFields)
        ff.mData = apply_stencil(ff.mData, ff.mWidth);
    for (auto& of : rOtherFields) {
        std::vector<std::uint8_t> dst(n_new * of.mRowBytes, 0);
        if (n_old > 0 && of.mRowBytes > 0)
            std::memcpy(dst.data(), of.mData.data(), n_old * of.mRowBytes);
        for (std::size_t ei = 0; ei < n_edges; ++ei) {
            if (of.mRowBytes > 0)
                std::memcpy(
                    dst.data() + (n_old + ei) * of.mRowBytes,
                    of.mData.data() + static_cast<std::size_t>(edges[ei].mLo) * of.mRowBytes,
                    of.mRowBytes);
        }
        of.mData = std::move(dst);
    }

    // Propagate creases to the two half-edges of every crease edge.
    std::unordered_set<std::uint64_t> next_creases;
    for (std::size_t ei = 0; ei < n_edges; ++ei) {
        if (edges[ei].mIsCrease) {
            const std::int64_t mid = static_cast<std::int64_t>(n_old + ei);
            next_creases.insert(subd_edge_key(edges[ei].mLo, mid));
            next_creases.insert(subd_edge_key(mid, edges[ei].mHi));
        }
    }
    rCreases = std::move(next_creases);

    std::vector<std::vector<std::int64_t>> next_faces;
    std::vector<std::int64_t> next_parents;
    next_faces.reserve(rFaces.size() * 4);
    next_parents.reserve(rFaces.size() * 4);
    for (std::size_t f = 0; f < rFaces.size(); ++f) {
        const std::int64_t a = rFaces[f][0];
        const std::int64_t b = rFaces[f][1];
        const std::int64_t c = rFaces[f][2];
        const std::int64_t eab =
            static_cast<std::int64_t>(n_old + edge_map.at(subd_edge_key(a, b)));
        const std::int64_t ebc =
            static_cast<std::int64_t>(n_old + edge_map.at(subd_edge_key(b, c)));
        const std::int64_t eca =
            static_cast<std::int64_t>(n_old + edge_map.at(subd_edge_key(c, a)));
        const std::int64_t par = rParents[f];
        next_faces.push_back({a, eab, eca});
        next_faces.push_back({b, ebc, eab});
        next_faces.push_back({c, eca, ebc});
        next_faces.push_back({eab, ebc, eca});
        for (int k = 0; k < 4; ++k)
            next_parents.push_back(par);
    }
    rFaces = std::move(next_faces);
    rParents = std::move(next_parents);
}

void subd_step_catmull_clark(std::vector<double>& rPts,
                             std::vector<std::vector<std::int64_t>>& rFaces,
                             std::vector<std::int64_t>& rParents,
                             std::unordered_set<std::uint64_t>& rCreases,
                             std::vector<SubdFloatField>& rFloatFields,
                             std::vector<SubdOtherField>& rOtherFields, bool PreserveBoundary,
                             std::int64_t& rNumCrease, std::int64_t& rNumBoundary) {
    const std::size_t n_old = rPts.size() / 3;
    const std::size_t n_faces = rFaces.size();
    std::unordered_map<std::uint64_t, std::size_t> edge_map;
    std::vector<SubdEdgeInfo> edges =
        subd_build_edges(rFaces, rCreases, PreserveBoundary, edge_map, rNumCrease, rNumBoundary);
    const std::size_t n_edges = edges.size();
    const std::size_t face_base = n_old;
    const std::size_t edge_base = n_old + n_faces;
    const std::size_t n_new = n_old + n_faces + n_edges;

    std::vector<std::vector<std::size_t>> v_faces(n_old);
    for (std::size_t f = 0; f < n_faces; ++f) {
        for (std::int64_t v : rFaces[f])
            v_faces[static_cast<std::size_t>(v)].push_back(f);
    }
    std::vector<std::vector<std::int64_t>> adj(n_old);
    std::vector<std::vector<std::int64_t>> crease_adj(n_old);
    std::vector<bool> non_manifold_vert(n_old, false);
    for (const auto& e : edges) {
        adj[static_cast<std::size_t>(e.mLo)].push_back(e.mHi);
        adj[static_cast<std::size_t>(e.mHi)].push_back(e.mLo);
        if (e.mFaces.size() >= 3) {
            non_manifold_vert[static_cast<std::size_t>(e.mLo)] = true;
            non_manifold_vert[static_cast<std::size_t>(e.mHi)] = true;
        }
        if (e.mIsCrease) {
            crease_adj[static_cast<std::size_t>(e.mLo)].push_back(e.mHi);
            crease_adj[static_cast<std::size_t>(e.mHi)].push_back(e.mLo);
        }
    }

    auto apply_stencil = [&](const std::vector<double>& rSrc, std::size_t width) {
        std::vector<double> dst(n_new * width, 0.0);
        // 1. Face points
        for (std::size_t f = 0; f < n_faces; ++f) {
            const auto& ring = rFaces[f];
            const double inv = 1.0 / static_cast<double>(ring.size());
            for (std::int64_t v : ring) {
                const std::size_t sv = static_cast<std::size_t>(v);
                for (std::size_t c = 0; c < width; ++c)
                    dst[(face_base + f) * width + c] += rSrc[sv * width + c] * inv;
            }
        }
        // 2. Edge points
        for (std::size_t ei = 0; ei < n_edges; ++ei) {
            const auto& e = edges[ei];
            const std::size_t u0 = static_cast<std::size_t>(e.mLo);
            const std::size_t u1 = static_cast<std::size_t>(e.mHi);
            const std::size_t out_v = edge_base + ei;
            if (e.mIsCrease || e.mFaces.size() != 2) {
                for (std::size_t c = 0; c < width; ++c)
                    dst[out_v * width + c] = 0.5 * (rSrc[u0 * width + c] + rSrc[u1 * width + c]);
            } else {
                const std::size_t f0 = face_base + static_cast<std::size_t>(e.mFaces[0]);
                const std::size_t f1 = face_base + static_cast<std::size_t>(e.mFaces[1]);
                for (std::size_t c = 0; c < width; ++c)
                    dst[out_v * width + c] = 0.25 * (rSrc[u0 * width + c] + rSrc[u1 * width + c] +
                                                     dst[f0 * width + c] + dst[f1 * width + c]);
            }
        }
        // 3. Old vertices
        for (std::size_t v = 0; v < n_old; ++v) {
            const auto& cadj = crease_adj[v];
            if (non_manifold_vert[v] || cadj.size() == 1 || cadj.size() >= 3) {
                for (std::size_t c = 0; c < width; ++c)
                    dst[v * width + c] = rSrc[v * width + c];
            } else if (cadj.size() == 2) {
                const std::size_t a = static_cast<std::size_t>(cadj[0]);
                const std::size_t b = static_cast<std::size_t>(cadj[1]);
                for (std::size_t c = 0; c < width; ++c)
                    dst[v * width + c] = 0.75 * rSrc[v * width + c] +
                                         0.125 * (rSrc[a * width + c] + rSrc[b * width + c]);
            } else {
                const std::size_t nf = v_faces[v].size();
                const std::size_t ne = adj[v].size();
                if (nf == 0 || ne == 0 || nf != ne) {
                    for (std::size_t c = 0; c < width; ++c)
                        dst[v * width + c] = rSrc[v * width + c];
                } else {
                    const double n = static_cast<double>(nf);
                    for (std::size_t c = 0; c < width; ++c) {
                        double f_avg = 0.0;
                        for (std::size_t f : v_faces[v])
                            f_avg += dst[(face_base + f) * width + c];
                        f_avg /= n;
                        double r_avg = 0.0;
                        for (std::int64_t u : adj[v])
                            r_avg += 0.5 * (rSrc[v * width + c] +
                                            rSrc[static_cast<std::size_t>(u) * width + c]);
                        r_avg /= n;
                        dst[v * width + c] =
                            (f_avg + 2.0 * r_avg + (n - 3.0) * rSrc[v * width + c]) / n;
                    }
                }
            }
        }
        return dst;
    };

    rPts = apply_stencil(rPts, 3);
    for (auto& ff : rFloatFields)
        ff.mData = apply_stencil(ff.mData, ff.mWidth);
    for (auto& of : rOtherFields) {
        std::vector<std::uint8_t> dst(n_new * of.mRowBytes, 0);
        if (n_old > 0 && of.mRowBytes > 0)
            std::memcpy(dst.data(), of.mData.data(), n_old * of.mRowBytes);
        for (std::size_t f = 0; f < n_faces; ++f) {
            if (of.mRowBytes > 0 && !rFaces[f].empty())
                std::memcpy(dst.data() + (face_base + f) * of.mRowBytes,
                            of.mData.data() + static_cast<std::size_t>(rFaces[f][0]) * of.mRowBytes,
                            of.mRowBytes);
        }
        for (std::size_t ei = 0; ei < n_edges; ++ei) {
            if (of.mRowBytes > 0)
                std::memcpy(
                    dst.data() + (edge_base + ei) * of.mRowBytes,
                    of.mData.data() + static_cast<std::size_t>(edges[ei].mLo) * of.mRowBytes,
                    of.mRowBytes);
        }
        of.mData = std::move(dst);
    }

    std::unordered_set<std::uint64_t> next_creases;
    for (std::size_t ei = 0; ei < n_edges; ++ei) {
        if (edges[ei].mIsCrease) {
            const std::int64_t mid = static_cast<std::int64_t>(edge_base + ei);
            next_creases.insert(subd_edge_key(edges[ei].mLo, mid));
            next_creases.insert(subd_edge_key(mid, edges[ei].mHi));
        }
    }
    rCreases = std::move(next_creases);

    std::vector<std::vector<std::int64_t>> next_faces;
    std::vector<std::int64_t> next_parents;
    for (std::size_t f = 0; f < n_faces; ++f) {
        const auto& ring = rFaces[f];
        const std::size_t m = ring.size();
        const std::int64_t fp = static_cast<std::int64_t>(face_base + f);
        const std::int64_t par = rParents[f];
        for (std::size_t k = 0; k < m; ++k) {
            const std::int64_t vk = ring[k];
            const std::int64_t v_next = ring[(k + 1) % m];
            const std::int64_t v_prev = ring[(k + m - 1) % m];
            const std::int64_t e_next =
                static_cast<std::int64_t>(edge_base + edge_map.at(subd_edge_key(vk, v_next)));
            const std::int64_t e_prev =
                static_cast<std::int64_t>(edge_base + edge_map.at(subd_edge_key(v_prev, vk)));
            next_faces.push_back({vk, e_next, fp, e_prev});
            next_parents.push_back(par);
        }
    }
    rFaces = std::move(next_faces);
    rParents = std::move(next_parents);
}

}  // namespace

SubdivisionScheme subdivision_scheme_from_name(const std::string& rName) {
    const std::string low = subd_lower(rName);
    if (low == "loop")
        return SubdivisionScheme::Loop;
    if (low == "catmull_clark" || low == "catmull-clark" || low == "catmullclark" || low == "cc")
        return SubdivisionScheme::CatmullClark;
    throw std::invalid_argument(std::string(kSubdPrefix) + "unknown scheme '" + rName +
                                "' (expected 'loop' or 'catmull_clark')");
}

SubdivideSurfaceResult subdivide_surface(const Mesh& rMesh,
                                         const SubdivideSurfaceOptions& rOptions) {
    if (rOptions.mLevels < 0)
        throw std::invalid_argument(std::string(kSubdPrefix) + "levels must be >= 0");
    if (!std::isfinite(rOptions.mCreaseAngle))
        throw std::invalid_argument(std::string(kSubdPrefix) + "crease_angle must be finite");

    SubdivideSurfaceResult res;
    if (rOptions.mLevels == 0) {
        res.mMesh = detail::clone_mesh(rMesh);
        return res;
    }

    const std::size_t n_pts_in = rMesh.NumPoints();
    std::size_t total_cells = 0;
    for (const auto cb : rMesh.CellRange())
        total_cells += cb.NumCells();

    std::vector<std::uint8_t> selected(total_cells, 1);
    detail::SurfaceRings rings = detail::surface_rings(rMesh, selected, kSubdPrefix);
    if (rings.NumFaces() == 0)
        throw std::invalid_argument(std::string(kSubdPrefix) +
                                    "mesh contains no 2-D surface cells");

    std::vector<std::vector<std::int64_t>> faces(rings.NumFaces());
    std::vector<std::int64_t> parents(rings.NumFaces());
    for (std::size_t f = 0; f < rings.NumFaces(); ++f) {
        const std::int64_t s0 = rings.mStart[f];
        const std::int64_t s1 = rings.mStart[f + 1];
        const std::size_t m = static_cast<std::size_t>(s1 - s0);
        if (rOptions.mScheme == SubdivisionScheme::Loop && m != 3)
            throw std::invalid_argument(std::string(kSubdPrefix) +
                                        "Loop subdivision requires triangles only (use "
                                        "'catmull_clark' or 'simplexify' first)");
        faces[f].assign(rings.mNodes.begin() + s0, rings.mNodes.begin() + s1);
        parents[f] = rings.mCell[f];
    }

    std::vector<double> pts = detail::surface_xyz(rMesh);

    // Initial creases via shared detail::crease_edges (takes per-face unit normals).
    std::unordered_set<std::uint64_t> creases;
    if (rOptions.mCreaseAngle > 0.0 && rOptions.mCreaseAngle < 180.0) {
        std::vector<double> face_normals(faces.size() * 3, 0.0);
        for (std::size_t f = 0; f < faces.size(); ++f) {
            const auto& ring = faces[f];
            const std::size_t m = ring.size();
            double nx = 0.0, ny = 0.0, nz = 0.0;
            for (std::size_t k = 0; k < m; ++k) {
                const std::size_t i0 = static_cast<std::size_t>(ring[k]);
                const std::size_t i1 = static_cast<std::size_t>(ring[(k + 1) % m]);
                const double x0 = pts[3 * i0 + 0], y0 = pts[3 * i0 + 1], z0 = pts[3 * i0 + 2];
                const double x1 = pts[3 * i1 + 0], y1 = pts[3 * i1 + 1], z1 = pts[3 * i1 + 2];
                nx += (y0 - y1) * (z0 + z1);
                ny += (z0 - z1) * (x0 + x1);
                nz += (x0 - x1) * (y0 + y1);
            }
            const double len = std::sqrt(nx * nx + ny * ny + nz * nz);
            if (len > 0.0) {
                face_normals[3 * f + 0] = nx / len;
                face_normals[3 * f + 1] = ny / len;
                face_normals[3 * f + 2] = nz / len;
            }
        }
        const std::vector<detail::CreaseEdge> ce =
            detail::crease_edges(rings.mStart, rings.mNodes, face_normals, rOptions.mCreaseAngle);
        for (const auto& e : ce) {
            if (e.IsCrease())
                creases.insert(subd_edge_key(e.mLo, e.mHi));
        }
    }

    // Extract point_data for interpolation / inheritance.
    std::vector<SubdFloatField> float_fields;
    std::vector<SubdOtherField> other_fields;
    for (const std::string& name : rMesh.PointDataNames()) {
        const NDArray& arr = rMesh.PointData(name);
        if (arr.Ndim() == 0 || arr.Shape()[0] != n_pts_in)
            continue;
        std::vector<std::size_t> tail(arr.Shape().begin() + 1, arr.Shape().end());
        std::size_t width = 1;
        for (std::size_t d : tail)
            width *= d;
        if (rOptions.mInterpolatePointData &&
            (arr.Dtype() == DType::Float64 || arr.Dtype() == DType::Float32)) {
            SubdFloatField ff;
            ff.mName = name;
            ff.mOrigDtype = arr.Dtype();
            ff.mTailShape = std::move(tail);
            ff.mWidth = width;
            ff.mData.resize(n_pts_in * width);
            if (arr.Dtype() == DType::Float64) {
                std::memcpy(ff.mData.data(), arr.As<double>(), n_pts_in * width * sizeof(double));
            } else {
                const float* sp = arr.As<float>();
                for (std::size_t i = 0; i < n_pts_in * width; ++i)
                    ff.mData[i] = static_cast<double>(sp[i]);
            }
            float_fields.push_back(std::move(ff));
        } else {
            SubdOtherField of;
            of.mName = name;
            of.mDtype = arr.Dtype();
            of.mTailShape = std::move(tail);
            of.mRowBytes = width * dtype_size(arr.Dtype());
            of.mData.resize(n_pts_in * of.mRowBytes);
            if (n_pts_in * of.mRowBytes > 0)
                std::memcpy(of.mData.data(), arr.Data(), n_pts_in * of.mRowBytes);
            other_fields.push_back(std::move(of));
        }
    }

    std::int64_t n_crease = 0;
    std::int64_t n_bnd = 0;
    for (int lvl = 0; lvl < rOptions.mLevels; ++lvl) {
        if (rOptions.mScheme == SubdivisionScheme::Loop) {
            subd_step_loop(pts, faces, parents, creases, float_fields, other_fields,
                           rOptions.mPreserveBoundary, n_crease, n_bnd);
        } else {
            subd_step_catmull_clark(pts, faces, parents, creases, float_fields, other_fields,
                                    rOptions.mPreserveBoundary, n_crease, n_bnd);
        }
    }
    res.mNumCreaseEdges = n_crease;
    res.mNumBoundaryEdges = n_bnd;

    // Assemble output Mesh.
    const std::size_t n_out_pts = pts.size() / 3;
    const std::size_t out_dim = rMesh.PointDim() == 2 ? 2 : 3;
    NDArray out_pts = NDArray::Uninit(DType::Float64, {n_out_pts, out_dim});
    auto* op = out_pts.As<double>();
    for (std::size_t i = 0; i < n_out_pts; ++i) {
        op[i * out_dim + 0] = pts[3 * i + 0];
        op[i * out_dim + 1] = pts[3 * i + 1];
        if (out_dim == 3)
            op[i * out_dim + 2] = pts[3 * i + 2];
    }

    Mesh out;
    out.AssignPoints(std::move(out_pts));

    const std::size_t n_out_faces = faces.size();
    const std::size_t k_nodes = rOptions.mScheme == SubdivisionScheme::Loop ? 3 : 4;
    const char* cell_type = rOptions.mScheme == SubdivisionScheme::Loop ? "triangle" : "quad";
    NDArray conn = NDArray::Uninit(DType::Int64, {n_out_faces, k_nodes});
    auto* cp = conn.As<std::int64_t>();
    for (std::size_t f = 0; f < n_out_faces; ++f) {
        for (std::size_t k = 0; k < k_nodes; ++k)
            cp[f * k_nodes + k] = faces[f][k];
    }
    out.AddCellBlock(cell_type, std::move(conn));

    // Attach point_data.
    for (const auto& ff : float_fields) {
        std::vector<std::size_t> shape = {n_out_pts};
        shape.insert(shape.end(), ff.mTailShape.begin(), ff.mTailShape.end());
        NDArray arr = NDArray::Uninit(ff.mOrigDtype, shape);
        if (ff.mOrigDtype == DType::Float64) {
            std::memcpy(arr.As<double>(), ff.mData.data(), n_out_pts * ff.mWidth * sizeof(double));
        } else {
            float* dp = arr.As<float>();
            for (std::size_t i = 0; i < n_out_pts * ff.mWidth; ++i)
                dp[i] = static_cast<float>(ff.mData[i]);
        }
        out.AddPointData(ff.mName, std::move(arr));
    }
    for (const auto& of : other_fields) {
        std::vector<std::size_t> shape = {n_out_pts};
        shape.insert(shape.end(), of.mTailShape.begin(), of.mTailShape.end());
        NDArray arr = NDArray::Uninit(of.mDtype, shape);
        if (n_out_pts * of.mRowBytes > 0)
            std::memcpy(arr.Data(), of.mData.data(), n_out_pts * of.mRowBytes);
        out.AddPointData(of.mName, std::move(arr));
    }

    // Replicate cell_data from parent cells.
    if (rMesh.NumCellBlocks() > 0) {
        std::vector<std::size_t> blk_base(rMesh.NumCellBlocks() + 1, 0);
        std::size_t bi = 0;
        for (const auto cb : rMesh.CellRange()) {
            blk_base[bi + 1] = blk_base[bi] + cb.NumCells();
            ++bi;
        }
        for (const std::string& name : rMesh.CellDataNames()) {
            if (rMesh.CellDataNumBlocks(name) != rMesh.NumCellBlocks())
                continue;
            const NDArray& first = rMesh.CellData(name, 0);
            std::vector<std::size_t> shape = {n_out_faces};
            if (first.Ndim() > 1)
                shape.insert(shape.end(), first.Shape().begin() + 1, first.Shape().end());
            std::size_t width = 1;
            for (std::size_t d = 1; d < shape.size(); ++d)
                width *= shape[d];
            const std::size_t row_bytes = width * dtype_size(first.Dtype());
            NDArray dst = NDArray::Uninit(first.Dtype(), shape);
            std::byte* dp = dst.Data();
            bool ok = true;
            for (std::size_t f = 0; f < n_out_faces; ++f) {
                const std::size_t g = static_cast<std::size_t>(parents[f]);
                const auto it = std::upper_bound(blk_base.begin(), blk_base.end(), g);
                const std::size_t b = static_cast<std::size_t>((it - blk_base.begin()) - 1);
                const std::size_t loc = g - blk_base[b];
                const NDArray& src = rMesh.CellData(name, b);
                if (src.Dtype() != first.Dtype()) {
                    ok = false;
                    break;
                }
                if (row_bytes > 0)
                    std::memcpy(dp + f * row_bytes, src.Data() + loc * row_bytes, row_bytes);
            }
            if (ok)
                out.AddCellData(name, {std::move(dst)});
        }
    }

    if (rOptions.mRecordParentIds) {
        NDArray par = NDArray::Uninit(DType::Int64, {n_out_faces});
        std::copy(parents.begin(), parents.end(), par.As<std::int64_t>());
        out.AddCellData(kSubdivideParentCellName, {std::move(par)});
    }

    for (const std::string& name : rMesh.FieldDataNames())
        out.AddFieldData(name, detail::data_owned_copy(rMesh.FieldData(name)));

    // Remap Point and Cell regions.
    std::vector<std::vector<std::int64_t>> parent_to_children(total_cells);
    for (std::size_t f = 0; f < n_out_faces; ++f)
        parent_to_children[static_cast<std::size_t>(parents[f])].push_back(
            static_cast<std::int64_t>(f));
    for (std::size_t r = 0; r < rMesh.NumRegions(); ++r) {
        const Region& reg = rMesh.Region(r);
        if (reg.mKind == RegionKind::Point) {
            out.AddRegion(Region(reg.mName, reg.mKind, reg.mDim, reg.mTag,
                                 detail::data_owned_copy(reg.mEntries)));
        } else if (reg.mKind == RegionKind::Cell) {
            std::vector<std::int64_t> mapped;
            const std::int64_t* ents = reg.Entries();
            for (std::size_t e = 0; e < reg.NumEntries(); ++e) {
                const std::int64_t c = ents[e];
                if (c >= 0 && static_cast<std::size_t>(c) < total_cells) {
                    const auto& ch = parent_to_children[static_cast<std::size_t>(c)];
                    mapped.insert(mapped.end(), ch.begin(), ch.end());
                }
            }
            NDArray arr = NDArray::Uninit(DType::Int64, {mapped.size()});
            std::copy(mapped.begin(), mapped.end(), arr.As<std::int64_t>());
            out.AddRegion(Region(reg.mName, RegionKind::Cell, reg.mDim, reg.mTag, std::move(arr)));
        }
    }

    res.mMesh = std::move(out);
    return res;
}

}  // namespace meshioplusplus

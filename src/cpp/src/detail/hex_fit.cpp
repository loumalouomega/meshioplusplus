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

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <utility>
#include <vector>

// Project includes
#include "hex_fit.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

// The faces of a hexahedron in the `cell_faces` order (outward normals).
constexpr int kHexfitFace[6][4] = {{0, 4, 7, 3}, {1, 2, 6, 5}, {0, 1, 5, 4},
                                   {3, 7, 6, 2}, {0, 3, 2, 1}, {4, 5, 6, 7}};

// Each corner's three edge neighbours, right-handed.
constexpr int kHexfitCorner[8][3] = {{1, 3, 4}, {2, 0, 5}, {3, 1, 6}, {0, 2, 7},
                                     {7, 5, 0}, {4, 6, 1}, {5, 7, 2}, {6, 4, 3}};

// The face-connected components of the cells of a 2x2x2 block in @p Mask.
int hexfit_components(unsigned Mask) {
    unsigned left = Mask;
    int n = 0;
    while (left != 0) {
        ++n;
        unsigned front = left & (~left + 1u);
        unsigned seen = front;
        while (front != 0) {
            unsigned next = 0;
            for (int b = 0; b < 8; ++b)
                if (front & (1u << b))
                    for (int axis = 0; axis < 3; ++axis)
                        next |= 1u << (b ^ (1 << axis));
            next &= left & ~seen;
            seen |= next;
            front = next;
        }
        left &= ~seen;
    }
    return n;
}

// Each VTK corner's position in the unit cube, as bits (x, y, z) = (1, 2, 4).
constexpr int kHexfitBits[8] = {0, 1, 3, 2, 4, 5, 7, 6};

// The VTK corner at each bit pattern (the inverse of kHexfitBits).
constexpr int kHexfitAtBits[8] = {0, 1, 3, 2, 4, 5, 7, 6};

// Per mask of kept cells in a 2x2x2 block, whether the boundary is manifold there.
std::array<std::uint8_t, 256> hexfit_valid_table() {
    std::array<std::uint8_t, 256> t{};
    for (unsigned m = 0; m < 256; ++m)
        t[m] =
            m == 0 || m == 255 || (hexfit_components(m) == 1 && hexfit_components(~m & 255u) == 1);
    return t;
}

double hexfit_corner(const std::vector<Vec3>& rP, const std::array<std::int64_t, 8>& rH, int C) {
    const Vec3& o = rP[static_cast<std::size_t>(rH[static_cast<std::size_t>(C)])];
    const Vec3 a = vec3_sub(rP[static_cast<std::size_t>(rH[kHexfitCorner[C][0]])], o);
    const Vec3 b = vec3_sub(rP[static_cast<std::size_t>(rH[kHexfitCorner[C][1]])], o);
    const Vec3 c = vec3_sub(rP[static_cast<std::size_t>(rH[kHexfitCorner[C][2]])], o);
    const double l = vec3_norm(a) * vec3_norm(b) * vec3_norm(c);
    if (!(l > 0.0))
        return -1.0;
    return triple_product(a, b, c) / l;
}

}  // namespace

std::int64_t hex_grid_make_manifold(const std::array<std::int64_t, 3>& rN,
                                    std::vector<std::uint8_t>& rKept,
                                    const std::vector<double>& rScore) {
    static const std::array<std::uint8_t, 256> valid = hexfit_valid_table();
    const std::int64_t nx = rN[0], ny = rN[1], nz = rN[2];
    auto cell = [&](std::int64_t I, std::int64_t J, std::int64_t K) -> std::int64_t {
        if (I < 0 || J < 0 || K < 0 || I >= nx || J >= ny || K >= nz)
            return -1;
        return (K * ny + J) * nx + I;
    };
    std::int64_t added = 0;
    for (bool changed = true; changed;) {
        changed = false;
        for (std::int64_t k = 0; k <= nz; ++k)
            for (std::int64_t j = 0; j <= ny; ++j)
                for (std::int64_t i = 0; i <= nx; ++i) {
                    std::array<std::int64_t, 8> ids{};
                    unsigned mask = 0;
                    for (int b = 0; b < 8; ++b) {
                        ids[static_cast<std::size_t>(b)] =
                            cell(i - 1 + (b & 1), j - 1 + ((b >> 1) & 1), k - 1 + ((b >> 2) & 1));
                        const std::int64_t c = ids[static_cast<std::size_t>(b)];
                        if (c >= 0 && rKept[static_cast<std::size_t>(c)])
                            mask |= 1u << b;
                    }
                    while (!valid[mask]) {
                        int best = -1;
                        for (int b = 0; b < 8; ++b) {
                            const std::int64_t c = ids[static_cast<std::size_t>(b)];
                            if (c < 0 || (mask & (1u << b)))
                                continue;
                            if (best < 0 || rScore[static_cast<std::size_t>(c)] >
                                                rScore[static_cast<std::size_t>(
                                                    ids[static_cast<std::size_t>(best)])])
                                best = b;
                        }
                        if (best < 0)
                            break;  // only cells outside the grid left to add
                        rKept[static_cast<std::size_t>(ids[static_cast<std::size_t>(best)])] = 1;
                        mask |= 1u << best;
                        ++added;
                        changed = true;
                    }
                }
    }
    return added;
}

namespace {

// `hex_make_manifold` and `hex_shrink_manifold`: around every node, add the
// empty hexahedron with the largest score (among @p pAllowed), or, with
// @p Shrink, remove the kept one with the smallest, until the node passes.
std::int64_t hexfit_star_repair(const HexMesh& rMesh, std::vector<std::uint8_t>& rKept,
                                const std::vector<double>& rScore,
                                const std::vector<std::uint8_t>* pAllowed, bool Shrink) {
    const std::size_t nh = rMesh.mHexes.size(), np = rMesh.mPoints.size();
    // The hexahedron across each facet (-1: none).
    std::vector<std::array<std::int64_t, 6>> across(nh);
    {
        std::vector<std::pair<std::array<std::int64_t, 4>, std::int64_t>> faces;
        faces.reserve(nh * 6);
        for (std::size_t h = 0; h < nh; ++h) {
            across[h].fill(-1);
            for (int f = 0; f < 6; ++f) {
                std::array<std::int64_t, 4> key{};
                for (int k = 0; k < 4; ++k)
                    key[static_cast<std::size_t>(k)] =
                        rMesh.mHexes[h][static_cast<std::size_t>(kHexfitFace[f][k])];
                std::sort(key.begin(), key.end());
                faces.emplace_back(key,
                                   static_cast<std::int64_t>(h * 6 + static_cast<std::size_t>(f)));
            }
        }
        std::sort(faces.begin(), faces.end());
        for (std::size_t i = 0; i + 1 < faces.size(); ++i)
            if (faces[i].first == faces[i + 1].first) {
                const std::int64_t a = faces[i].second, b = faces[i + 1].second;
                across[static_cast<std::size_t>(a / 6)][static_cast<std::size_t>(a % 6)] = b / 6;
                across[static_cast<std::size_t>(b / 6)][static_cast<std::size_t>(b % 6)] = a / 6;
            }
    }
    std::vector<std::int64_t> start(np + 1, 0), at;
    for (const auto& h : rMesh.mHexes)
        for (std::int64_t v : h)
            ++start[static_cast<std::size_t>(v) + 1];
    for (std::size_t v = 0; v < np; ++v)
        start[v + 1] += start[v];
    at.resize(static_cast<std::size_t>(start[np]));
    {
        std::vector<std::int64_t> fill(start.begin(), start.end() - 1);
        for (std::size_t h = 0; h < nh; ++h)
            for (std::int64_t v : rMesh.mHexes[h])
                at[static_cast<std::size_t>(fill[static_cast<std::size_t>(v)]++)] =
                    static_cast<std::int64_t>(h);
    }
    std::int64_t added = 0;
    std::vector<std::vector<std::size_t>> link;
    std::vector<std::uint8_t> seen;
    // The components of the star's cells with `Kept` (the outside, index n,
    // counting as empty).
    auto components = [&](const std::vector<std::int64_t>& rStar, bool Kept, bool Outside) {
        const std::size_t n = rStar.size();
        auto member = [&](std::size_t I) {
            if (I == n)
                return !Kept && Outside;
            return (rKept[static_cast<std::size_t>(rStar[I])] != 0) == Kept;
        };
        seen.assign(n + 1, 0);
        int count = 0;
        for (std::size_t i = 0; i <= n; ++i) {
            if (seen[i] || !member(i))
                continue;
            ++count;
            std::vector<std::size_t> front{i};
            seen[i] = 1;
            while (!front.empty()) {
                const std::size_t u = front.back();
                front.pop_back();
                for (std::size_t w : link[u])
                    if (!seen[w] && member(w)) {
                        seen[w] = 1;
                        front.push_back(w);
                    }
            }
        }
        return count;
    };
    for (bool changed = true; changed;) {
        changed = false;
        for (std::size_t p = 0; p < np; ++p) {
            const std::vector<std::int64_t> star(at.begin() + start[p], at.begin() + start[p + 1]);
            const std::size_t n = star.size();
            if (n == 0)
                continue;
            link.assign(n + 1, {});
            bool outside = false;
            for (std::size_t i = 0; i < n; ++i) {
                const auto& hex = rMesh.mHexes[static_cast<std::size_t>(star[i])];
                for (int f = 0; f < 6; ++f) {
                    bool has = false;
                    for (int k = 0; k < 4 && !has; ++k)
                        has = hex[static_cast<std::size_t>(kHexfitFace[f][k])] ==
                              static_cast<std::int64_t>(p);
                    if (!has)
                        continue;
                    const std::int64_t g =
                        across[static_cast<std::size_t>(star[i])][static_cast<std::size_t>(f)];
                    const std::size_t j =
                        g < 0 ? n
                              : static_cast<std::size_t>(std::find(star.begin(), star.end(), g) -
                                                         star.begin());
                    outside = outside || g < 0;
                    link[i].push_back(j);
                    link[j].push_back(i);
                }
            }
            for (;;) {
                if (components(star, true, outside) <= 1 && components(star, false, outside) <= 1)
                    break;
                std::int64_t best = -1;
                for (std::int64_t h : star) {
                    const auto u = static_cast<std::size_t>(h);
                    if (Shrink) {
                        if (rKept[u] &&
                            (best < 0 || rScore[u] < rScore[static_cast<std::size_t>(best)]))
                            best = h;
                    } else if (!rKept[u] && (pAllowed == nullptr || (*pAllowed)[u]) &&
                               (best < 0 || rScore[u] > rScore[static_cast<std::size_t>(best)])) {
                        best = h;
                    }
                }
                if (best < 0)
                    break;  // only the outside (or nothing allowed) left to add
                rKept[static_cast<std::size_t>(best)] = Shrink ? 0 : 1;
                ++added;
                changed = true;
            }
        }
    }
    return added;
}

}  // namespace

std::int64_t hex_make_manifold(const HexMesh& rMesh, std::vector<std::uint8_t>& rKept,
                               const std::vector<double>& rScore,
                               const std::vector<std::uint8_t>* pAllowed) {
    return hexfit_star_repair(rMesh, rKept, rScore, pAllowed, false);
}

std::int64_t hex_shrink_manifold(const HexMesh& rMesh, std::vector<std::uint8_t>& rKept,
                                 const std::vector<double>& rScore) {
    return hexfit_star_repair(rMesh, rKept, rScore, nullptr, true);
}

HexInterface hex_pillow_interface(HexMesh& rMesh, const std::vector<std::uint8_t>& rInner) {
    const std::size_t nh = rMesh.mHexes.size(), np = rMesh.mPoints.size();
    HexInterface out;
    out.mShell.assign(np, -1);
    out.mInner.assign(np, -1);
    out.mFirstHex = nh;
    std::vector<std::pair<std::array<std::int64_t, 4>, std::int64_t>> faces;
    faces.reserve(nh * 6);
    for (std::size_t h = 0; h < nh; ++h)
        for (std::int64_t f = 0; f < 6; ++f) {
            std::array<std::int64_t, 4> key{};
            for (std::size_t c = 0; c < 4; ++c)
                key[c] = rMesh.mHexes[h][static_cast<std::size_t>(kHexfitFace[f][c])];
            std::sort(key.begin(), key.end());
            faces.emplace_back(key, static_cast<std::int64_t>(h) * 6 + f);
        }
    parallel_sort(faces.begin(), faces.end());
    // The faces between an inner hexahedron and an outer one, from the inner.
    for (std::size_t i = 0; i + 1 < faces.size(); ++i) {
        if (faces[i].first != faces[i + 1].first)
            continue;
        std::int64_t a = faces[i].second, b = faces[i + 1].second;
        const bool ia = rInner[static_cast<std::size_t>(a / 6)] != 0;
        const bool ib = rInner[static_cast<std::size_t>(b / 6)] != 0;
        if (ia == ib)
            continue;
        if (!ia)
            std::swap(a, b);
        const auto h = static_cast<std::size_t>(a / 6);
        const auto f = static_cast<std::size_t>(a % 6);
        std::array<std::int64_t, 4> q{};
        for (std::size_t c = 0; c < 4; ++c)
            q[c] = rMesh.mHexes[h][static_cast<std::size_t>(kHexfitFace[f][c])];
        out.mQuads.push_back(q);
        out.mInnerHex.push_back(static_cast<std::int64_t>(h));
    }
    for (const auto& q : out.mQuads)
        for (std::int64_t v : q)
            out.mShell[static_cast<std::size_t>(v)] = 0;
    for (std::size_t v = 0; v < np; ++v)
        if (out.mShell[v] == 0) {
            out.mShell[v] = static_cast<std::int64_t>(rMesh.mPoints.size());
            rMesh.mPoints.push_back(rMesh.mPoints[v]);
            out.mInner[v] = static_cast<std::int64_t>(rMesh.mPoints.size());
            rMesh.mPoints.push_back(rMesh.mPoints[v]);
        }
    // The inner hexahedra move onto the inner copies (all of them at a node:
    // the inner set is manifold, so they meet the outer ones only there).
    for (std::size_t h = 0; h < nh; ++h)
        if (rInner[h])
            for (std::int64_t& v : rMesh.mHexes[h])
                if (out.mInner[static_cast<std::size_t>(v)] >= 0)
                    v = out.mInner[static_cast<std::size_t>(v)];
    // A face points out of its inner hexahedron, towards the outer side: the
    // inner layer runs from the inner copies to the shell nodes, the outer one
    // from the shell nodes to the outer nodes, both on the face's positive side.
    for (const auto& q : out.mQuads) {
        std::array<std::int64_t, 8> inner{}, outer{};
        for (std::size_t c = 0; c < 4; ++c) {
            const auto v = static_cast<std::size_t>(q[c]);
            inner[c] = out.mInner[v];
            inner[c + 4] = out.mShell[v];
            outer[c] = out.mShell[v];
            outer[c + 4] = q[c];
        }
        rMesh.mHexes.push_back(inner);
        rMesh.mHexes.push_back(outer);
    }
    out.mShell.resize(rMesh.mPoints.size(), -1);
    out.mInner.resize(rMesh.mPoints.size(), -1);
    return out;
}

HexBoundary hex_boundary(const HexMesh& rMesh) {
    std::vector<std::pair<std::array<std::int64_t, 4>, std::int64_t>> faces;
    faces.reserve(rMesh.mHexes.size() * 6);
    for (std::size_t h = 0; h < rMesh.mHexes.size(); ++h)
        for (std::int64_t f = 0; f < 6; ++f) {
            std::array<std::int64_t, 4> key{};
            for (std::size_t c = 0; c < 4; ++c)
                key[c] = rMesh.mHexes[h][static_cast<std::size_t>(kHexfitFace[f][c])];
            std::sort(key.begin(), key.end());
            faces.emplace_back(key, static_cast<std::int64_t>(h) * 6 + f);
        }
    parallel_sort(faces.begin(), faces.end());
    HexBoundary out;
    for (std::size_t i = 0; i < faces.size();) {
        std::size_t j = i;
        while (j < faces.size() && faces[j].first == faces[i].first)
            ++j;
        if (j - i == 1) {
            const auto h = static_cast<std::size_t>(faces[i].second / 6);
            const auto f = static_cast<std::size_t>(faces[i].second % 6);
            std::array<std::int64_t, 4> q{};
            for (std::size_t c = 0; c < 4; ++c)
                q[c] = rMesh.mHexes[h][static_cast<std::size_t>(kHexfitFace[f][c])];
            out.mQuads.push_back(q);
            out.mSides.push_back({static_cast<std::int64_t>(h), static_cast<std::int64_t>(f)});
        }
        i = j;
    }
    std::vector<std::array<std::int64_t, 2>> edges;
    edges.reserve(out.mQuads.size() * 4);
    for (const auto& q : out.mQuads)
        for (std::size_t e = 0; e < 4; ++e)
            edges.push_back({std::min(q[e], q[(e + 1) % 4]), std::max(q[e], q[(e + 1) % 4])});
    std::sort(edges.begin(), edges.end());
    for (std::size_t i = 0; i < edges.size();) {
        std::size_t j = i;
        while (j < edges.size() && edges[j] == edges[i])
            ++j;
        if (j - i == 1)
            ++out.mNumOpenEdges;
        else if (j - i > 2)
            ++out.mNumNonManifoldEdges;
        i = j;
    }
    return out;
}

std::vector<std::int64_t> hex_pillow(HexMesh& rMesh, const HexBoundary& rBoundary) {
    const std::size_t np = rMesh.mPoints.size();
    std::vector<std::int64_t> copy(np, -1);
    for (const auto& q : rBoundary.mQuads)
        for (std::int64_t v : q)
            copy[static_cast<std::size_t>(v)] = 0;
    for (std::size_t p = 0; p < np; ++p)
        if (copy[p] == 0) {
            copy[p] = static_cast<std::int64_t>(rMesh.mPoints.size());
            rMesh.mPoints.push_back(rMesh.mPoints[p]);
        }
    // A boundary face's normal points out of its hexahedron, so its copy lies
    // on the positive side of it: (a, b, c, d) is the new hexahedron's bottom.
    for (const auto& q : rBoundary.mQuads) {
        std::array<std::int64_t, 8> h{};
        for (std::size_t c = 0; c < 4; ++c) {
            h[c] = q[c];
            h[c + 4] = copy[static_cast<std::size_t>(q[c])];
        }
        rMesh.mHexes.push_back(h);
    }
    copy.resize(rMesh.mPoints.size(), -1);
    return copy;
}

double hex_min_scaled_jacobian(const std::vector<Vec3>& rPoints,
                               const std::array<std::int64_t, 8>& rHex) {
    double lo = 1.0;
    for (int c = 0; c < 8; ++c)
        lo = std::min(lo, hexfit_corner(rPoints, rHex, c));
    return lo;
}

std::array<double, 8> hex_corner_scaled_jacobians(const std::vector<Vec3>& rPoints,
                                                  const std::array<std::int64_t, 8>& rHex) {
    std::array<double, 8> out{};
    for (int c = 0; c < 8; ++c)
        out[static_cast<std::size_t>(c)] = hexfit_corner(rPoints, rHex, c);
    return out;
}

std::vector<Vec3> hex_boundary_normals(const std::vector<Vec3>& rPoints,
                                       const std::vector<std::array<std::int64_t, 4>>& rQuads,
                                       int Smooth) {
    const std::size_t np = rPoints.size();
    std::vector<Vec3> n(np, Vec3{0.0, 0.0, 0.0});
    std::vector<std::vector<std::int64_t>> nbr(np);
    for (const auto& q : rQuads) {
        const Vec3 d1 = vec3_sub(rPoints[static_cast<std::size_t>(q[2])],
                                 rPoints[static_cast<std::size_t>(q[0])]);
        const Vec3 d2 = vec3_sub(rPoints[static_cast<std::size_t>(q[3])],
                                 rPoints[static_cast<std::size_t>(q[1])]);
        const Vec3 area = vec3_scale(vec3_cross(d1, d2), 0.5);
        for (std::size_t e = 0; e < 4; ++e) {
            const auto v = static_cast<std::size_t>(q[e]);
            n[v] = vec3_add(n[v], area);
            nbr[v].push_back(q[(e + 1) % 4]);
            nbr[v].push_back(q[(e + 3) % 4]);
        }
    }
    for (auto& a : nbr) {
        std::sort(a.begin(), a.end());
        a.erase(std::unique(a.begin(), a.end()), a.end());
    }
    for (std::size_t v = 0; v < np; ++v)
        if (vec3_norm(n[v]) > 0.0)
            n[v] = vec3_normalize(n[v]);
    for (int s = 0; s < Smooth; ++s) {
        std::vector<Vec3> next(np, Vec3{0.0, 0.0, 0.0});
        for (std::size_t v = 0; v < np; ++v) {
            if (nbr[v].empty())
                continue;
            Vec3 sum = n[v];
            for (std::int64_t u : nbr[v])
                sum = vec3_add(sum, n[static_cast<std::size_t>(u)]);
            next[v] = vec3_norm(sum) > 0.0 ? vec3_normalize(sum) : Vec3{0.0, 0.0, 0.0};
        }
        n.swap(next);
    }
    return n;
}

HexUntangleReport hex_untangle(std::vector<Vec3>& rPoints,
                               const std::vector<std::array<std::int64_t, 8>>& rHexes,
                               const std::vector<HexNodeMode>& rMode, const HexProjector& rProject,
                               const HexUntangleOptions& rOptions) {
    const std::size_t np = rPoints.size();
    // Incident hexahedra and edge neighbours of every node.
    std::vector<std::vector<std::int64_t>> inc(np), nbr(np);
    for (std::size_t h = 0; h < rHexes.size(); ++h)
        for (int c = 0; c < 8; ++c) {
            const auto v = static_cast<std::size_t>(rHexes[h][static_cast<std::size_t>(c)]);
            inc[v].push_back(static_cast<std::int64_t>(h));
            for (int e = 0; e < 3; ++e)
                nbr[v].push_back(rHexes[h][static_cast<std::size_t>(kHexfitCorner[c][e])]);
        }
    for (auto& n : nbr) {
        std::sort(n.begin(), n.end());
        n.erase(std::unique(n.begin(), n.end()), n.end());
    }
    auto worst = [&](std::size_t V) {
        double w = 1.0;
        for (std::int64_t h : inc[V])
            w = std::min(w, hex_min_scaled_jacobian(rPoints, rHexes[static_cast<std::size_t>(h)]));
        return w;
    };
    // The robust tie-break: the corners around node @p V, each capped at the
    // target, summed. A hexahedron with two poor corners at different nodes
    // keeps its minimum when either node moves alone, so a strict rise of
    // the minimum never comes; a move that keeps the minimum and raises this
    // sum still makes progress, and the pair only grows, so it cannot cycle.
    auto capped = [&](std::size_t V) {
        double sum = 0.0;
        for (std::int64_t h : inc[V])
            for (int c = 0; c < 8; ++c)
                sum += std::min(rOptions.mTarget,
                                hexfit_corner(rPoints, rHexes[static_cast<std::size_t>(h)], c));
        return sum;
    };
    // The mean of the trilinear reflections of node @p V: per incident
    // hexahedron, the position that cancels its twist (the alternating sum
    // of its corners), which a parallelepiped's corner satisfies.
    auto reflection = [&](std::size_t V) {
        Vec3 sum{0.0, 0.0, 0.0};
        for (std::int64_t h : inc[V]) {
            const auto& hex = rHexes[static_cast<std::size_t>(h)];
            int at = 0;
            while (static_cast<std::size_t>(hex[static_cast<std::size_t>(at)]) != V)
                ++at;
            const int b0 = kHexfitBits[at];
            Vec3 p{0.0, 0.0, 0.0};
            for (int b = 0; b < 8; ++b) {
                if (b == b0)
                    continue;
                int flips = 0;
                for (int k = 0; k < 3; ++k)
                    flips += ((b ^ b0) >> k) & 1;
                const Vec3& x = rPoints[static_cast<std::size_t>(
                    hex[static_cast<std::size_t>(kHexfitAtBits[b])])];
                p = flips % 2 ? vec3_add(p, x) : vec3_sub(p, x);
            }
            sum = vec3_add(sum, p);
        }
        return vec3_scale(sum, 1.0 / static_cast<double>(inc[V].size()));
    };
    HexUntangleReport rep;
    // Only nodes whose neighbourhood changed (or that are still below the
    // target) can move: the others would be skipped again, so the result is
    // that of visiting every node.
    std::vector<std::uint8_t> active(np, 1), next(np, 0);
    std::int64_t step_base = 0;  // the sweep the step size counts from
    bool restarted = false;
    std::vector<Vec3> cand;
    for (std::int64_t sweep = 0; sweep < rOptions.mIterations; ++sweep) {
        ++rep.mSweeps;
        const double scale = 0.2 * std::pow(0.5, static_cast<double>((sweep - step_base) / 10));
        std::int64_t moves = 0;
        double low = 1.0;  // the worst corner of the nodes visited
        std::fill(next.begin(), next.end(), 0);
        for (std::size_t v = 0; v < np; ++v) {
            if (!active[v] || rMode[v] == HexNodeMode::Fixed || inc[v].empty())
                continue;
            const double now = worst(v);
            if (now >= rOptions.mTarget)
                continue;
            next[v] = 1;
            const Vec3 keep = rPoints[v];
            double edge = 0.0;
            Vec3 mean{0.0, 0.0, 0.0};
            for (std::int64_t q : nbr[v]) {
                const Vec3& p = rPoints[static_cast<std::size_t>(q)];
                edge += vec3_norm(vec3_sub(p, keep));
                mean = vec3_add(mean, p);
            }
            edge /= static_cast<double>(nbr[v].size());
            mean = vec3_scale(mean, 1.0 / static_cast<double>(nbr[v].size()));
            cand.clear();
            cand.push_back(mean);
            const bool wide = rOptions.mRobust && now < 0.2;
            if (rOptions.mRobust)
                cand.push_back(reflection(v));
            if (wide) {
                for (int dz = -1; dz <= 1; ++dz)
                    for (int dy = -1; dy <= 1; ++dy)
                        for (int dx = -1; dx <= 1; ++dx) {
                            if (dx == 0 && dy == 0 && dz == 0)
                                continue;
                            const Vec3 d = vec3_normalize(Vec3{static_cast<double>(dx),
                                                               static_cast<double>(dy),
                                                               static_cast<double>(dz)});
                            cand.push_back(vec3_add(keep, vec3_scale(d, scale * edge)));
                        }
            } else {
                for (std::size_t k = 0; k < 6; ++k) {
                    Vec3 c = keep;
                    c[k / 2] += (k % 2 ? -1.0 : 1.0) * scale * edge;
                    cand.push_back(c);
                }
            }
            double best = now;
            double best_sum = rOptions.mRobust ? capped(v) : 0.0;
            Vec3 best_at = keep;
            auto attempt = [&](Vec3 C) {
                if (rMode[v] == HexNodeMode::Slide)
                    C = rProject(v, C);
                rPoints[v] = C;
                const double w = worst(v);
                if (w > best + 1e-12) {
                    best = w;
                    best_at = C;
                    if (rOptions.mRobust)
                        best_sum = capped(v);
                    return true;
                }
                if (rOptions.mRobust && w >= best) {
                    const double sum = capped(v);
                    if (sum > best_sum + 1e-12) {
                        best = w;
                        best_sum = sum;
                        best_at = C;
                        return true;
                    }
                }
                return false;
            };
            for (const Vec3& c : cand)
                attempt(c);
            if (wide && best_at != keep) {
                // A line search along the best step: further, then shorter.
                const Vec3 d = vec3_sub(best_at, keep);
                if (!attempt(vec3_add(keep, vec3_scale(d, 2.0))))
                    attempt(vec3_add(keep, vec3_scale(d, 0.5)));
            }
            rPoints[v] = best_at;
            low = std::min(low, best);
            if (best_at != keep) {
                ++moves;
                for (std::int64_t h : inc[v])
                    for (std::int64_t q : rHexes[static_cast<std::size_t>(h)]) {
                        active[static_cast<std::size_t>(q)] = 1;  // later in this sweep
                        next[static_cast<std::size_t>(q)] = 1;
                    }
            }
        }
        active.swap(next);
        rep.mMoves += moves;
        if (moves == 0) {
            if (!rOptions.mRobust || restarted || low >= 0.2)
                break;
            // Stalled with a poor corner left: start the steps again, once.
            restarted = true;
            step_base = sweep + 1;
        }
    }
    return rep;
}

}  // namespace detail
}  // namespace meshioplusplus

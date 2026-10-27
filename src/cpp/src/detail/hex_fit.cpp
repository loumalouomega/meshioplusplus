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

std::int64_t hex_make_manifold(const HexMesh& rMesh, std::vector<std::uint8_t>& rKept,
                               const std::vector<double>& rScore) {
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
                for (std::int64_t h : star)
                    if (!rKept[static_cast<std::size_t>(h)] &&
                        (best < 0 || rScore[static_cast<std::size_t>(h)] >
                                         rScore[static_cast<std::size_t>(best)]))
                        best = h;
                if (best < 0)
                    break;  // only the outside left to add
                rKept[static_cast<std::size_t>(best)] = 1;
                ++added;
                changed = true;
            }
        }
    }
    return added;
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

HexUntangleReport hex_untangle(std::vector<Vec3>& rPoints,
                               const std::vector<std::array<std::int64_t, 8>>& rHexes,
                               const std::vector<HexNodeMode>& rMode,
                               const std::function<Vec3(const Vec3&)>& rProject,
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
    HexUntangleReport rep;
    // Only nodes whose neighbourhood changed (or that are still below the
    // target) can move: the others would be skipped again, so the result is
    // that of visiting every node.
    std::vector<std::uint8_t> active(np, 1), next(np, 0);
    for (std::int64_t sweep = 0; sweep < rOptions.mIterations; ++sweep) {
        ++rep.mSweeps;
        const double scale = 0.2 * std::pow(0.5, static_cast<double>(sweep / 10));
        std::int64_t moves = 0;
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
            std::array<Vec3, 7> cand{};
            cand[0] = mean;
            for (std::size_t k = 0; k < 6; ++k) {
                cand[k + 1] = keep;
                cand[k + 1][k / 2] += (k % 2 ? -1.0 : 1.0) * scale * edge;
            }
            double best = now;
            Vec3 best_at = keep;
            for (Vec3 c : cand) {
                if (rMode[v] == HexNodeMode::Slide)
                    c = rProject(c);
                rPoints[v] = c;
                const double w = worst(v);
                if (w > best + 1e-12) {
                    best = w;
                    best_at = c;
                }
            }
            rPoints[v] = best_at;
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
        if (moves == 0)
            break;
    }
    return rep;
}

}  // namespace detail
}  // namespace meshioplusplus

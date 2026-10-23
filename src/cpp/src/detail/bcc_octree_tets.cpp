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
// The graded BCC tetrahedralisation of detail/bcc_octree_tets.hpp.
// Anonymous-namespace helpers are prefixed `bcco_` (the amalgamation
// concatenates every translation unit into one).

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "bcc_octree_tets.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

using BccoPoint = std::array<std::int64_t, 3>;

constexpr double kBccoRadToDeg = 57.29577951308232087680;

std::uint64_t bcco_key(const BccoPoint& rP) {
    return linoct_key(rP[0], rP[1], rP[2]);
}

// One tetrahedron before numbering: its vertices' doubled coordinates and the
// types of its non-corner vertices.
struct BccoTet {
    BccoPoint mV[4];
    BccVertexType mType[4];
    std::uint8_t mShape = 0;
};

// The exact orientation of four integer points: the sign of
// `(b - a) · ((c - a) × (d - a))`, computed on the differences divided by
// their common power of two so that nothing overflows.
int bcco_orientation(const BccoPoint* pV) {
    std::int64_t d[3][3];
    std::int64_t g = 0;
    for (std::size_t r = 0; r < 3; ++r)
        for (std::size_t k = 0; k < 3; ++k) {
            d[r][k] = pV[r + 1][k] - pV[0][k];
            g = std::gcd(g, d[r][k] < 0 ? -d[r][k] : d[r][k]);
        }
    if (g > 1)
        for (auto& row : d)
            for (std::int64_t& v : row)
                v /= g;
    const std::int64_t det = d[0][0] * (d[1][1] * d[2][2] - d[1][2] * d[2][1]) -
                             d[0][1] * (d[1][0] * d[2][2] - d[1][2] * d[2][0]) +
                             d[0][2] * (d[1][0] * d[2][1] - d[1][1] * d[2][0]);
    return (det > 0) - (det < 0);
}

// Emit `(apex0, apex1, a, b)` positively oriented.
void bcco_emit(std::vector<BccoTet>& rOut, const BccoPoint& rP0, BccVertexType T0,
               const BccoPoint& rP1, BccVertexType T1, const BccoPoint& rA, const BccoPoint& rB,
               std::uint8_t Shape) {
    BccoTet t;
    t.mV[0] = rP0;
    t.mV[1] = rP1;
    t.mV[2] = rA;
    t.mV[3] = rB;
    t.mType[0] = T0;
    t.mType[1] = T1;
    t.mType[2] = BccVertexType::Corner;
    t.mType[3] = BccVertexType::Corner;
    t.mShape = Shape;
    const int o = bcco_orientation(t.mV);
    if (o == 0)
        throw std::logic_error("meshio++: bcc_octree_lattice: a degenerate tetrahedron");
    if (o < 0) {
        std::swap(t.mV[2], t.mV[3]);
        std::swap(t.mType[2], t.mType[3]);
    }
    rOut.push_back(t);
}

}  // namespace

std::array<double, 6> tet_dihedral_angles(const Vec3& rA, const Vec3& rB, const Vec3& rC,
                                          const Vec3& rD) {
    const Vec3 p[4] = {rA, rB, rC, rD};
    // The edge (i, j) has the opposite edge (k, l); its dihedral angle is the
    // angle between the normals of faces (i, j, k) and (i, j, l).
    constexpr int kEdge[6][4] = {{0, 1, 2, 3}, {0, 2, 1, 3}, {0, 3, 1, 2},
                                 {1, 2, 0, 3}, {1, 3, 0, 2}, {2, 3, 0, 1}};
    std::array<double, 6> out{};
    for (std::size_t e = 0; e < 6; ++e) {
        const Vec3& a = p[kEdge[e][0]];
        const Vec3& b = p[kEdge[e][1]];
        const Vec3 ab = vec3_sub(b, a);
        const Vec3 n1 = vec3_cross(ab, vec3_sub(p[kEdge[e][2]], a));
        const Vec3 n2 = vec3_cross(ab, vec3_sub(p[kEdge[e][3]], a));
        const double c = vec3_dot(n1, n2) / std::sqrt(vec3_norm_sq(n1) * vec3_norm_sq(n2));
        out[e] = std::acos(std::clamp(c, -1.0, 1.0)) * kBccoRadToDeg;
    }
    std::sort(out.begin(), out.end());
    return out;
}

BccLattice bcc_octree_lattice(const LinearOctree& rTree, std::int64_t MaxTets) {
    for (std::size_t k = 0; k < 3; ++k)
        if (rTree.Extent(k) > (std::int64_t{1} << 20))
            throw std::invalid_argument(
                "meshio++: bcc_octree_lattice: the lattice is too fine for doubled keys");
    const std::size_t n = rTree.Size();
    // Every leaf corner, on the doubled lattice.
    std::vector<std::uint64_t> corner(n * 8);
    parallel_for(n, [&](std::size_t i) {
        const auto lo = rTree.Lo(i);
        const std::int64_t s2 = 2 * rTree.Span(i);
        for (std::size_t c = 0; c < 8; ++c)
            corner[i * 8 + c] =
                bcco_key({2 * lo[0] + ((c & 1) ? s2 : 0), 2 * lo[1] + ((c & 2) ? s2 : 0),
                          2 * lo[2] + ((c & 4) ? s2 : 0)});
    });
    parallel_sort(corner.begin(), corner.end());
    corner.erase(std::unique(corner.begin(), corner.end()), corner.end());
    auto is_corner = [&](const BccoPoint& rP) {
        return std::binary_search(corner.begin(), corner.end(), bcco_key(rP));
    };
    const auto faces = linear_octree_directions(OctreeAdjacency::Face);

    // Every leaf's tetrahedra, in a fixed order.
    std::vector<std::vector<BccoTet>> per_leaf(n);
    parallel_for(
        n,
        [&](std::size_t i) {
            std::vector<BccoTet>& out = per_leaf[i];
            const auto lo = rTree.Lo(i);
            const std::int64_t s = rTree.Span(i);  // finest units
            const BccoPoint centre{2 * lo[0] + s, 2 * lo[1] + s, 2 * lo[2] + s};
            for (const auto& d : faces) {
                const std::size_t k = d[0] ? 0 : (d[1] ? 1 : 2);
                const int sg = d[0] + d[1] + d[2];
                const std::size_t u = (k + 1) % 3, v = (k + 2) % 3;
                const std::int64_t plane = 2 * lo[k] + (sg > 0 ? 2 * s : 0);
                const auto nb = linear_octree_neighbours(rTree, i, d[0], d[1], d[2]);
                auto at = [&](std::int64_t A, std::int64_t B) {
                    BccoPoint p{};
                    p[k] = plane;
                    p[u] = A;
                    p[v] = B;
                    return p;
                };
                // The boundary of a square [u0, u0 + w] x [v0, v0 + w] on the
                // face, as segments, each cut at its midpoint when a corner
                // lies there; `Cut` receives whether it was.
                auto loop = [&](std::int64_t U0, std::int64_t V0, std::int64_t W,
                                std::vector<std::pair<BccoPoint, BccoPoint>>& rSeg,
                                std::vector<std::uint8_t>& rCut) {
                    const BccoPoint q[4] = {at(U0, V0), at(U0 + W, V0), at(U0 + W, V0 + W),
                                            at(U0, V0 + W)};
                    for (std::size_t e = 0; e < 4; ++e) {
                        const BccoPoint& a = q[e];
                        const BccoPoint& b = q[(e + 1) % 4];
                        BccoPoint m{};
                        for (std::size_t c = 0; c < 3; ++c)
                            m[c] = (a[c] + b[c]) / 2;
                        if (W % 2 == 0 && is_corner(m)) {
                            // A second cut would need leaves two levels apart.
                            BccoPoint m1{}, m2{};
                            for (std::size_t c = 0; c < 3; ++c) {
                                m1[c] = (a[c] + m[c]) / 2;
                                m2[c] = (m[c] + b[c]) / 2;
                            }
                            if ((W / 2) % 2 == 0 && (is_corner(m1) || is_corner(m2)))
                                throw std::logic_error(
                                    "meshio++: bcc_octree_lattice: the tree is not 2:1-balanced "
                                    "across edges");
                            rSeg.emplace_back(a, m);
                            rSeg.emplace_back(m, b);
                            rCut.push_back(1);
                            rCut.push_back(1);
                        } else {
                            rSeg.emplace_back(a, b);
                            rCut.push_back(0);
                        }
                    }
                };
                const std::int64_t u0 = 2 * lo[u], v0 = 2 * lo[v], w = 2 * s;
                if (nb.size() == 1 &&
                    rTree.mLevel[static_cast<std::size_t>(nb[0])] == rTree.mLevel[i]) {
                    // The same size across: the BCC tetrahedra around the axis
                    // between the two centres, emitted by the lower leaf.
                    if (sg < 0)
                        continue;
                    BccoPoint other = centre;
                    other[k] += 2 * s;
                    std::vector<std::pair<BccoPoint, BccoPoint>> seg;
                    std::vector<std::uint8_t> cut;
                    loop(u0, v0, w, seg, cut);
                    for (std::size_t e = 0; e < seg.size(); ++e)
                        bcco_emit(out, centre, BccVertexType::Centre, other, BccVertexType::Centre,
                                  seg[e].first, seg[e].second, cut[e] ? 2 : 1);
                    continue;
                }
                if (nb.size() <= 1) {
                    // The domain boundary, or a coarser leaf across: one portion.
                    BccoPoint q = at(u0 + s, v0 + s);
                    std::vector<std::pair<BccoPoint, BccoPoint>> seg;
                    std::vector<std::uint8_t> cut;
                    loop(u0, v0, w, seg, cut);
                    for (std::size_t e = 0; e < seg.size(); ++e)
                        bcco_emit(out, centre, BccVertexType::Centre, q, BccVertexType::FaceCentre,
                                  seg[e].first, seg[e].second, cut[e] ? 4 : 3);
                    continue;
                }
                // Four finer leaves across: four portions, their faces.
                if (nb.size() != 4)
                    throw std::logic_error(
                        "meshio++: bcc_octree_lattice: the tree is not 2:1-balanced across faces");
                const std::int64_t h = s;  // the quarter square's edge, doubled
                for (std::int64_t qv = 0; qv < 2; ++qv)
                    for (std::int64_t qu = 0; qu < 2; ++qu) {
                        const std::int64_t a0 = u0 + qu * h, b0 = v0 + qv * h;
                        const BccoPoint q = at(a0 + h / 2, b0 + h / 2);
                        std::vector<std::pair<BccoPoint, BccoPoint>> seg;
                        std::vector<std::uint8_t> cut;
                        loop(a0, b0, h, seg, cut);
                        for (std::size_t e = 0; e < seg.size(); ++e) {
                            if (cut[e])
                                throw std::logic_error(
                                    "meshio++: bcc_octree_lattice: the tree is not 2:1-balanced "
                                    "across edges");
                            // Outer when the segment lies on the leaf's own edge.
                            const BccoPoint& a = seg[e].first;
                            const BccoPoint& b = seg[e].second;
                            const bool outer = (a[u] == b[u] && (a[u] == u0 || a[u] == u0 + w)) ||
                                               (a[v] == b[v] && (a[v] == v0 || a[v] == v0 + w));
                            bcco_emit(out, centre, BccVertexType::Centre, q,
                                      BccVertexType::FaceCentre, a, b, outer ? 5 : 6);
                        }
                    }
            }
        },
        16);

    // Count, check the budget, then number the vertices by key.
    std::vector<std::int64_t> count(n), offset(n);
    for (std::size_t i = 0; i < n; ++i)
        count[i] = static_cast<std::int64_t>(per_leaf[i].size());
    const std::int64_t total =
        parallel_exclusive_scan(count.data(), n, offset.data(), std::int64_t{0});
    if (MaxTets > 0 && total > MaxTets)
        throw std::invalid_argument("meshio++: octree_mesh: the lattice would have " +
                                    std::to_string(total) + " tetrahedra, more than max_tets (" +
                                    std::to_string(MaxTets) + ")");
    std::vector<std::pair<std::uint64_t, std::uint8_t>> keys(static_cast<std::size_t>(total) * 4);
    parallel_for(n, [&](std::size_t i) {
        const auto at = static_cast<std::size_t>(offset[i]) * 4;
        for (std::size_t t = 0; t < per_leaf[i].size(); ++t)
            for (std::size_t c = 0; c < 4; ++c)
                keys[at + t * 4 + c] = {bcco_key(per_leaf[i][t].mV[c]),
                                        static_cast<std::uint8_t>(per_leaf[i][t].mType[c])};
    });
    std::vector<std::pair<std::uint64_t, std::uint8_t>> unique_keys = keys;
    parallel_sort(unique_keys.begin(), unique_keys.end());
    unique_keys.erase(
        std::unique(unique_keys.begin(), unique_keys.end(),
                    [](const auto& rA, const auto& rB) { return rA.first == rB.first; }),
        unique_keys.end());
    BccLattice out;
    const std::size_t nv = unique_keys.size();
    out.mKey.resize(nv);
    out.mType.resize(nv);
    out.mIjk2.resize(nv);
    for (std::size_t p = 0; p < nv; ++p) {
        out.mKey[p] = unique_keys[p].first;
        out.mType[p] = static_cast<BccVertexType>(unique_keys[p].second);
        out.mIjk2[p] = linoct_decode(unique_keys[p].first);
    }
    out.mTets.resize(static_cast<std::size_t>(total));
    out.mShape.resize(static_cast<std::size_t>(total));
    out.mLeaf.resize(static_cast<std::size_t>(total));
    parallel_for(n, [&](std::size_t i) {
        const auto at = static_cast<std::size_t>(offset[i]);
        for (std::size_t t = 0; t < per_leaf[i].size(); ++t) {
            for (std::size_t c = 0; c < 4; ++c)
                out.mTets[at + t][c] =
                    static_cast<std::int64_t>(std::lower_bound(out.mKey.begin(), out.mKey.end(),
                                                               keys[(at + t) * 4 + c].first) -
                                              out.mKey.begin());
            out.mShape[at + t] = per_leaf[i][t].mShape;
            out.mLeaf[at + t] = static_cast<std::int64_t>(i);
        }
    });
    return out;
}

}  // namespace detail
}  // namespace meshioplusplus

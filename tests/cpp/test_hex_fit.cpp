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
// detail/hex_fit: manifold repair leaves every grid vertex's 2x2x2 block with
// face-connected kept and empty cells (checked by brute force) and only adds
// cells; the repaired cells' boundary is closed and edge-manifold and
// encloses their volume; a pillow layer keeps every hexahedron positive with
// one face on the new boundary; the untangler recovers an inverted node; and
// its robust mode leaves a plateau the plain pattern search stalls on, where
// a hexahedron has flat corners at two nodes; and the boundary normals of a
// box point out of its faces, edges and corners.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <random>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "../../src/cpp/src/detail/hex_fit.hpp"

namespace {

namespace det = meshioplusplus::detail;
using det::Vec3;

// Face-connected components of the cells of a 2x2x2 block, by flood fill.
int hf_components(const std::array<bool, 8>& rIn) {
    std::array<int, 8> label{};
    label.fill(-1);
    int n = 0;
    for (int s = 0; s < 8; ++s) {
        if (!rIn[s] || label[s] >= 0)
            continue;
        std::vector<int> stack{s};
        label[s] = n;
        while (!stack.empty()) {
            const int b = stack.back();
            stack.pop_back();
            for (int axis = 0; axis < 3; ++axis) {
                const int o = b ^ (1 << axis);
                if (rIn[o] && label[o] < 0) {
                    label[o] = n;
                    stack.push_back(o);
                }
            }
        }
        ++n;
    }
    return n;
}

det::HexMesh hf_grid(const std::array<std::int64_t, 3>& rN,
                     const std::vector<std::uint8_t>& rKept) {
    det::HexMesh m;
    const std::int64_t vx = rN[0] + 1, vy = rN[1] + 1;
    for (std::int64_t k = 0; k <= rN[2]; ++k)
        for (std::int64_t j = 0; j <= rN[1]; ++j)
            for (std::int64_t i = 0; i <= rN[0]; ++i)
                m.mPoints.push_back(
                    {static_cast<double>(i), static_cast<double>(j), static_cast<double>(k)});
    auto v = [&](std::int64_t I, std::int64_t J, std::int64_t K) { return (K * vy + J) * vx + I; };
    for (std::int64_t k = 0; k < rN[2]; ++k)
        for (std::int64_t j = 0; j < rN[1]; ++j)
            for (std::int64_t i = 0; i < rN[0]; ++i)
                if (rKept[static_cast<std::size_t>((k * rN[1] + j) * rN[0] + i)])
                    m.mHexes.push_back({v(i, j, k), v(i + 1, j, k), v(i + 1, j + 1, k),
                                        v(i, j + 1, k), v(i, j, k + 1), v(i + 1, j, k + 1),
                                        v(i + 1, j + 1, k + 1), v(i, j + 1, k + 1)});
    return m;
}

// The volume a closed, outward quad boundary encloses (each quad as two triangles).
double hf_enclosed(const det::HexMesh& rM, const det::HexBoundary& rB) {
    double vol = 0.0;
    for (const auto& q : rB.mQuads) {
        const Vec3& a = rM.mPoints[static_cast<std::size_t>(q[0])];
        const Vec3& b = rM.mPoints[static_cast<std::size_t>(q[1])];
        const Vec3& c = rM.mPoints[static_cast<std::size_t>(q[2])];
        const Vec3& d = rM.mPoints[static_cast<std::size_t>(q[3])];
        vol += (det::triple_product(a, b, c) + det::triple_product(a, c, d)) / 6.0;
    }
    return vol;
}

}  // namespace

TEST(HexFit, ManifoldRepairLeavesEveryBlockManifoldAndOnlyAdds) {
    std::mt19937 rng(11u);
    const std::array<std::int64_t, 3> n{{7, 6, 8}};
    for (int trial = 0; trial < 40; ++trial) {
        const std::size_t num = static_cast<std::size_t>(n[0] * n[1] * n[2]);
        std::vector<std::uint8_t> kept(num, 0);
        std::vector<double> score(num);
        for (std::size_t c = 0; c < num; ++c) {
            kept[c] = rng() % 100 < static_cast<unsigned>(20 + trial);
            score[c] = static_cast<double>(rng() % 1000);
        }
        const auto before = kept;
        const std::int64_t added = det::hex_grid_make_manifold(n, kept, score);
        std::int64_t diff = 0;
        for (std::size_t c = 0; c < num; ++c) {
            EXPECT_GE(kept[c], before[c]);
            diff += kept[c] - before[c];
        }
        EXPECT_EQ(diff, added);
        // Every interior vertex's block (blocks at the grid's edge may keep a
        // non-manifold pattern against the outside, which cannot be filled).
        for (std::int64_t k = 1; k < n[2]; ++k)
            for (std::int64_t j = 1; j < n[1]; ++j)
                for (std::int64_t i = 1; i < n[0]; ++i) {
                    std::array<bool, 8> in{}, out{};
                    int count = 0;
                    for (int b = 0; b < 8; ++b) {
                        const std::int64_t c =
                            ((k - 1 + (b >> 2)) * n[1] + (j - 1 + ((b >> 1) & 1))) * n[0] +
                            (i - 1 + (b & 1));
                        in[b] = kept[static_cast<std::size_t>(c)] != 0;
                        out[b] = !in[b];
                        count += in[b];
                    }
                    if (count == 0 || count == 8)
                        continue;
                    EXPECT_EQ(hf_components(in), 1) << "trial " << trial;
                    EXPECT_EQ(hf_components(out), 1) << "trial " << trial;
                }
        // With an empty border, the repaired cells' boundary is closed and manifold.
        std::vector<std::uint8_t> inner(num, 0);
        for (std::int64_t k = 1; k + 1 < n[2]; ++k)
            for (std::int64_t j = 1; j + 1 < n[1]; ++j)
                for (std::int64_t i = 1; i + 1 < n[0]; ++i) {
                    const std::size_t c = static_cast<std::size_t>((k * n[1] + j) * n[0] + i);
                    inner[c] = before[c];
                }
        for (std::size_t c = 0; c < num; ++c)
            score[c] = 0.0;
        det::hex_grid_make_manifold(n, inner, score);
        const det::HexMesh m = hf_grid(n, inner);
        if (m.mHexes.empty())
            continue;
        const det::HexBoundary b = det::hex_boundary(m);
        EXPECT_EQ(b.mNumOpenEdges, 0);
        EXPECT_EQ(b.mNumNonManifoldEdges, 0) << "trial " << trial;
        EXPECT_NEAR(hf_enclosed(m, b), static_cast<double>(m.mHexes.size()), 1e-9);
    }
}

TEST(HexFit, PillowKeepsEveryHexahedronPositiveWithOneBoundaryFace) {
    const std::array<std::int64_t, 3> n{{5, 5, 5}};
    std::vector<std::uint8_t> kept(125, 0);
    for (std::int64_t k = 1; k < 4; ++k)
        for (std::int64_t j = 1; j < 4; ++j)
            for (std::int64_t i = 1; i < 4; ++i)
                kept[static_cast<std::size_t>((k * 5 + j) * 5 + i)] = 1;
    det::HexMesh m = hf_grid(n, kept);
    const det::HexBoundary b0 = det::hex_boundary(m);
    const std::size_t grid_points = m.mPoints.size();
    const std::size_t grid_hexes = m.mHexes.size();
    const std::vector<std::int64_t> copy = det::hex_pillow(m, b0);
    EXPECT_EQ(m.mHexes.size(), grid_hexes + b0.mQuads.size());
    // Push every copy out of the box's centre by a third of a cell.
    for (std::size_t p = 0; p < grid_points; ++p) {
        if (copy[p] < 0)
            continue;
        Vec3& c = m.mPoints[static_cast<std::size_t>(copy[p])];
        for (std::size_t k = 0; k < 3; ++k)
            c[k] += c[k] < 2.5 ? -1.0 / 3.0 : 1.0 / 3.0;
    }
    for (const auto& h : m.mHexes)
        EXPECT_GT(det::hex_min_scaled_jacobian(m.mPoints, h), 0.0);
    const det::HexBoundary b1 = det::hex_boundary(m);
    EXPECT_EQ(b1.mNumOpenEdges, 0);
    EXPECT_EQ(b1.mNumNonManifoldEdges, 0);
    EXPECT_EQ(b1.mQuads.size(), b0.mQuads.size());
    std::vector<int> faces_on(m.mHexes.size(), 0);
    for (const auto& s : b1.mSides) {
        EXPECT_GE(static_cast<std::size_t>(s[0]), grid_hexes);
        EXPECT_EQ(s[1], 5);  // the copies' face
        ++faces_on[static_cast<std::size_t>(s[0])];
    }
    for (std::size_t h = grid_hexes; h < m.mHexes.size(); ++h)
        EXPECT_EQ(faces_on[h], 1);
    // The outer box spans [1 - 1/3, 4 + 1/3] on every axis.
    const double edge = 11.0 / 3.0;
    EXPECT_NEAR(hf_enclosed(m, b1), edge * edge * edge, 1e-9);
}

TEST(HexFit, UntangleRecoversAnInvertedNode) {
    const std::array<std::int64_t, 3> n{{3, 3, 3}};
    std::vector<std::uint8_t> kept(27, 1);
    det::HexMesh m = hf_grid(n, kept);
    // The interior node (1, 1, 1) pushed through its neighbours.
    const std::size_t centre = (1 * 4 + 1) * 4 + 1;
    m.mPoints[centre] = {2.6, 2.4, 2.7};
    double before = 1.0;
    for (const auto& h : m.mHexes)
        before = std::min(before, det::hex_min_scaled_jacobian(m.mPoints, h));
    EXPECT_LT(before, 0.0);
    std::vector<det::HexNodeMode> mode(m.mPoints.size(), det::HexNodeMode::Fixed);
    mode[centre] = det::HexNodeMode::Free;
    const auto rep = det::hex_untangle(
        m.mPoints, m.mHexes, mode, [](std::size_t, const Vec3& rP) { return rP; },
        det::HexUntangleOptions{});
    double after = 1.0;
    for (const auto& h : m.mHexes)
        after = std::min(after, det::hex_min_scaled_jacobian(m.mPoints, h));
    EXPECT_GT(after, 0.5);
    EXPECT_GT(rep.mMoves, 0);
    // A fixed node never moves.
    EXPECT_EQ(m.mPoints[0], (Vec3{0.0, 0.0, 0.0}));
}

// Two nodes of one hexahedron each put on the line through two of their edge
// neighbours: both corners are flat, so moving either node alone cannot raise
// the hexahedron's smallest corner, and the plain search, which takes a move
// only when that rises, stalls. The robust search also takes a move that keeps
// the smallest corner and raises the corners around the node.
TEST(HexFit, RobustUntangleLeavesATwoNodePlateau) {
    const std::array<std::int64_t, 3> n{{3, 3, 3}};
    std::vector<std::uint8_t> kept(27, 1);
    const det::HexMesh grid = hf_grid(n, kept);
    auto node = [](int I, int J, int K) { return static_cast<std::size_t>((K * 4 + J) * 4 + I); };
    const std::size_t a = node(1, 1, 1), b = node(1, 1, 2);
    auto run = [&](bool Robust) {
        det::HexMesh m = grid;
        m.mPoints[a] = {1.5, 1.5, 1.0};
        m.mPoints[b] = {1.5, 1.5, 2.0};
        std::vector<det::HexNodeMode> mode(m.mPoints.size(), det::HexNodeMode::Fixed);
        mode[a] = mode[b] = det::HexNodeMode::Free;
        det::HexUntangleOptions o;
        o.mRobust = Robust;
        det::hex_untangle(
            m.mPoints, m.mHexes, mode, [](std::size_t, const Vec3& rP) { return rP; }, o);
        double worst = 1.0;
        for (const auto& h : m.mHexes)
            worst = std::min(worst, det::hex_min_scaled_jacobian(m.mPoints, h));
        return worst;
    };
    const double plain = run(false), robust = run(true);
    EXPECT_LE(plain, 0.0);
    EXPECT_GT(robust, 0.5);
}

// The normals of a box's boundary: a face's interior nodes point along its
// axis, an edge's along the sum of its two faces', a corner's along the sum
// of its three; smoothing keeps them unit and outward.
TEST(HexFit, BoundaryNormalsPointOutOfABox) {
    const std::array<std::int64_t, 3> n{{3, 3, 3}};
    std::vector<std::uint8_t> kept(27, 1);
    const det::HexMesh m = hf_grid(n, kept);
    const det::HexBoundary b = det::hex_boundary(m);
    auto node = [](int I, int J, int K) { return static_cast<std::size_t>((K * 4 + J) * 4 + I); };
    const std::vector<Vec3> raw = det::hex_boundary_normals(m.mPoints, b.mQuads, 0);
    const double r3 = 1.0 / std::sqrt(3.0), r2 = 1.0 / std::sqrt(2.0);
    auto near = [](const Vec3& rA, const Vec3& rB) {
        for (std::size_t k = 0; k < 3; ++k)
            if (std::abs(rA[k] - rB[k]) > 1e-12)
                return false;
        return true;
    };
    EXPECT_TRUE(near(raw[node(1, 1, 0)], Vec3{0.0, 0.0, -1.0}));
    EXPECT_TRUE(near(raw[node(3, 1, 0)], Vec3{r2, 0.0, -r2}));
    EXPECT_TRUE(near(raw[node(3, 3, 3)], Vec3{r3, r3, r3}));
    EXPECT_TRUE(near(raw[node(1, 1, 1)], Vec3{0.0, 0.0, 0.0}));  // interior
    const std::vector<Vec3> smooth = det::hex_boundary_normals(m.mPoints, b.mQuads, 2);
    for (const auto& q : b.mQuads)
        for (std::int64_t v : q) {
            const Vec3& x = smooth[static_cast<std::size_t>(v)];
            EXPECT_NEAR(det::vec3_norm(x), 1.0, 1e-12);
            const Vec3 out = det::vec3_sub(m.mPoints[static_cast<std::size_t>(v)], {1.5, 1.5, 1.5});
            EXPECT_GT(det::vec3_dot(x, out), 0.0);
        }
}

// The interface pillow: on a 4x4x4 block whose inner 2x2x2 is the inner set,
// every hexahedron stays positive, the inner ones move onto the copies, and
// each new hexahedron has exactly one face on the shell nodes' surface.
TEST(HexFit, InterfacePillowKeepsEveryHexahedronPositive) {
    const std::array<std::int64_t, 3> n{{4, 4, 4}};
    std::vector<std::uint8_t> kept(64, 1);
    det::HexMesh m = hf_grid(n, kept);
    std::vector<std::uint8_t> inner(m.mHexes.size(), 0);
    for (std::int64_t k = 1; k < 3; ++k)
        for (std::int64_t j = 1; j < 3; ++j)
            for (std::int64_t i = 1; i < 3; ++i)
                inner[static_cast<std::size_t>((k * 4 + j) * 4 + i)] = 1;
    const std::size_t before = m.mHexes.size();
    const det::HexInterface itf = det::hex_pillow_interface(m, inner);
    EXPECT_EQ(itf.mQuads.size(), 24u);  // the inner cube's 6 faces of 4 quads
    EXPECT_EQ(m.mHexes.size(), before + 2 * 24);
    // Move the layers apart as a fit would: shell nodes stay, inner copies go
    // a quarter towards the centre, outer nodes a quarter away.
    const Vec3 centre{2.0, 2.0, 2.0};
    for (std::size_t v = 0; v < itf.mShell.size(); ++v) {
        if (itf.mShell[v] < 0)
            continue;
        const Vec3 x = m.mPoints[v];
        const Vec3 d = det::vec3_sub(x, centre);
        m.mPoints[static_cast<std::size_t>(itf.mInner[v])] =
            det::vec3_sub(x, det::vec3_scale(d, 0.25 / det::vec3_norm(d)));
        m.mPoints[v] = det::vec3_add(x, det::vec3_scale(d, 0.25 / det::vec3_norm(d)));
    }
    for (const auto& h : m.mHexes)
        EXPECT_GT(det::hex_min_scaled_jacobian(m.mPoints, h), 0.0);
    // The inner hexahedra use no original interface node any more.
    for (std::size_t h = 0; h < before; ++h)
        if (inner[h])
            for (std::int64_t v : m.mHexes[h])
                EXPECT_LT(itf.mShell[static_cast<std::size_t>(v)], 0);
    EXPECT_EQ(det::hex_boundary(m).mNumOpenEdges, 0);
}

// The masked repair adds only allowed hexahedra; without a mask it is the
// plain one; and the shrinking repair always ends manifold.
TEST(HexFit, MaskedAndShrinkingRepairs) {
    const std::array<std::int64_t, 3> n{{4, 4, 4}};
    std::vector<std::uint8_t> all(64, 1);
    const det::HexMesh m = hf_grid(n, all);
    std::mt19937 rng(7);
    for (int trial = 0; trial < 20; ++trial) {
        std::vector<std::uint8_t> kept(64), allowed(64);
        std::vector<double> score(64);
        for (std::size_t c = 0; c < 64; ++c) {
            kept[c] = rng() % 3 == 0;
            allowed[c] = rng() % 2 == 0;
            score[c] = static_cast<double>(rng() % 1000);
        }
        std::vector<std::uint8_t> a = kept, b = kept, s = kept;
        det::hex_make_manifold(m, a, score);
        det::hex_make_manifold(m, b, score, nullptr);
        EXPECT_EQ(a, b);
        std::vector<std::uint8_t> masked = kept;
        det::hex_make_manifold(m, masked, score, &allowed);
        for (std::size_t c = 0; c < 64; ++c)
            if (masked[c] && !kept[c])
                EXPECT_TRUE(allowed[c]);
        det::hex_shrink_manifold(m, s, score);
        for (std::size_t c = 0; c < 64; ++c)
            EXPECT_LE(s[c], kept[c]);
        // Shrunk: a second repair finds nothing to do.
        std::vector<std::uint8_t> again = s;
        EXPECT_EQ(det::hex_make_manifold(m, again, score, &s), 0);
        EXPECT_EQ(again, s);
    }
}

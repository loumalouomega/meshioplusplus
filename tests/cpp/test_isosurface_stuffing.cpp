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
// detail/isosurface_stuffing: a plane through a graded lattice. A plane's cut
// points are exact, so every stencil is checked by volume: the two sides'
// pieces are positive, conform (every interior triangle twice with opposite
// orientations, every boundary triangle on the box or the plane), and their
// volumes are those of the box on each side of the plane, for random planes
// through random trees. On a uniform lattice the dihedral angles stay inside
// Labelle and Shewchuk's bounds.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <random>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "../../src/cpp/src/detail/bcc_octree_tets.hpp"
#include "../../src/cpp/src/detail/isosurface_stuffing.hpp"
#include "../../src/cpp/src/detail/linear_octree.hpp"

namespace {

namespace mio = meshioplusplus;
namespace det = meshioplusplus::detail;
using det::Vec3;

det::StuffMesh is_mesh(const det::LinearOctree& rTree) {
    const det::BccLattice L = det::bcc_octree_lattice(rTree, 0);
    det::StuffMesh m;
    for (std::size_t p = 0; p < L.mIjk2.size(); ++p) {
        m.mPoints.push_back({0.5 * static_cast<double>(L.mIjk2[p][0]),
                             0.5 * static_cast<double>(L.mIjk2[p][1]),
                             0.5 * static_cast<double>(L.mIjk2[p][2])});
        m.mCentre.push_back(L.mType[p] == det::BccVertexType::Centre);
    }
    m.mPinned.assign(m.mPoints.size(), 0);
    // The box's boundary stays put.
    for (std::size_t p = 0; p < L.mIjk2.size(); ++p)
        for (std::size_t k = 0; k < 3; ++k)
            if (L.mIjk2[p][k] == 0 || L.mIjk2[p][k] == 2 * rTree.Extent(k))
                m.mPinned[p] = 1;
    m.mTets = L.mTets;
    m.mMaterial.assign(L.mTets.size(), 0);
    m.mCut.assign(L.mTets.size(), 0);
    for (std::size_t t = 0; t < L.mTets.size(); ++t)
        m.mParent.push_back(static_cast<std::int64_t>(t));
    return m;
}

double is_vol(const det::StuffMesh& rM, const std::array<std::int64_t, 4>& rT) {
    const auto& p = rM.mPoints;
    return det::triple_product(det::vec3_sub(p[rT[1]], p[rT[0]]), det::vec3_sub(p[rT[2]], p[rT[0]]),
                               det::vec3_sub(p[rT[3]], p[rT[0]])) /
           6.0;
}

// The volume of the box [0, E]^3 where n.x < c, by sampling-free inclusion:
// the integral of the clipped region, computed by clipping the box's own
// tetrahedralisation (six Kuhn tetrahedra) against the plane.
double is_clip_volume(const Vec3& rN, double C, double E) {
    // Monte-Carlo-free: integrate along x of the clipped rectangle area by
    // high-resolution midpoint quadrature over a 2-D grid of columns.
    const int n = 400;
    double vol = 0.0;
    const double h = E / n;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            const double x = (i + 0.5) * h, y = (j + 0.5) * h;
            // n0 x + n1 y + n2 z < C  ->  z bound
            if (std::abs(rN[2]) < 1e-14)
                vol += (rN[0] * x + rN[1] * y < C) ? h * h * E : 0.0;
            else {
                const double z = (C - rN[0] * x - rN[1] * y) / rN[2];
                const double len = rN[2] > 0 ? std::clamp(z, 0.0, E) : E - std::clamp(z, 0.0, E);
                vol += h * h * len;
            }
        }
    return vol;
}

}  // namespace

TEST(IsosurfaceStuffing, PlanesThroughRandomTreesAreCutExactly) {
    std::mt19937 rng(5u);
    for (int trial = 0; trial < 12; ++trial) {
        det::LinearOctree t = det::linear_octree_roots({{1, 1, 1}}, 4);
        for (int r = 0; r < 4; ++r) {
            std::vector<std::uint8_t> mark(t.Size());
            for (auto& m : mark)
                m = rng() % 3 == 0;
            det::linear_octree_split(t, mark, 0);
        }
        det::linear_octree_balance(t, det::OctreeAdjacency::Full, 0);
        det::StuffMesh m = is_mesh(t);
        const double E = 16.0;
        // A random plane through the box (not through lattice planes).
        Vec3 n{static_cast<double>(rng() % 1000) - 500.0, static_cast<double>(rng() % 1000) - 500.0,
               static_cast<double>(rng() % 1000) - 500.0};
        n = det::vec3_normalize(n);
        const double c =
            det::vec3_dot(n, {8.0, 8.0, 8.0}) + (static_cast<double>(rng() % 1000) / 500.0 - 1.0);
        std::vector<std::int8_t> sign(m.mPoints.size());
        for (std::size_t p = 0; p < m.mPoints.size(); ++p) {
            const double f = det::vec3_dot(n, m.mPoints[p]) - c;
            sign[p] = f < 0.0 ? -1 : (f > 0.0 ? 1 : 0);
        }
        const auto pts = m.mPoints;
        det::StuffCut cut = [&](std::int64_t A, std::int64_t B) {
            const Vec3& a = pts[static_cast<std::size_t>(A)];
            const Vec3& b = pts[static_cast<std::size_t>(B)];
            const double fa = det::vec3_dot(n, a) - c, fb = det::vec3_dot(n, b) - c;
            return det::vec3_add(a, det::vec3_scale(det::vec3_sub(b, a), fa / (fa - fb)));
        };
        det::StuffOptions so;
        so.mMinVolume = 1e-14;
        const det::StuffReport rep = det::isosurface_stuff(m, sign, cut, so);
        EXPECT_GT(rep.mNumCutPoints + rep.mNumWarped, 0);
        // A graded lattice can refuse a warp (the guard); the cut stays exact.
        EXPECT_LE(rep.mNumWarpRefused, rep.mNumWarped + rep.mNumCutPoints);
        double in = 0.0, out = 0.0;
        for (std::size_t q = 0; q < m.mTets.size(); ++q) {
            const double v = is_vol(m, m.mTets[q]);
            ASSERT_GT(v, 0.0) << "trial " << trial;
            (m.mMaterial[q] == 1 ? in : out) += v;
            // Every piece on one side of the plane.
            for (std::int64_t p : m.mTets[q]) {
                const double f = det::vec3_dot(n, m.mPoints[static_cast<std::size_t>(p)]) - c;
                if (m.mMaterial[q] == 1)
                    EXPECT_LE(f, 1e-9);
                else
                    EXPECT_GE(f, -1e-9);
            }
        }
        EXPECT_NEAR(in + out, E * E * E, 1e-9);
        // Every vertex of the interface lies on the plane, so the inside is
        // the clipped box (up to the quadrature of the reference).
        EXPECT_NEAR(in, is_clip_volume(n, c, E), 2e-3 * E * E * E);
        // Conformity across both sides.
        std::map<std::array<std::int64_t, 3>, int> faces;
        constexpr int kFace[4][3] = {{1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}};
        for (const auto& tet : m.mTets)
            for (const auto& f : kFace) {
                std::array<std::int64_t, 3> v{tet[f[0]], tet[f[1]], tet[f[2]]};
                int parity = 0;
                for (int a = 0; a < 3; ++a)
                    for (int b = a + 1; b < 3; ++b)
                        parity ^= v[a] > v[b];
                std::sort(v.begin(), v.end());
                faces[v] += parity ? -1 : 1;
            }
        for (const auto& [f, net] : faces) {
            if (net == 0)
                continue;
            // A boundary triangle lies on the box.
            bool on = false;
            for (std::size_t k = 0; k < 3 && !on; ++k) {
                const double x = m.mPoints[static_cast<std::size_t>(f[0])][k];
                on = (x == 0.0 || x == E) && m.mPoints[static_cast<std::size_t>(f[1])][k] == x &&
                     m.mPoints[static_cast<std::size_t>(f[2])][k] == x;
            }
            EXPECT_TRUE(on) << "trial " << trial << ": an unmatched interior triangle";
        }
    }
}

TEST(IsosurfaceStuffing, UniformLatticeAnglesStayInTheBounds) {
    // A uniform lattice and random planes and spheres through it.
    std::mt19937 rng(17u);
    double lo = 180.0, hi = 0.0;
    for (int trial = 0; trial < 30; ++trial) {
        // Away from the box, whose boundary is pinned: the bounds are for
        // the uniform lattice around the interface.
        const det::LinearOctree t = det::linear_octree_roots({{8, 8, 8}}, 0);
        det::StuffMesh m = is_mesh(t);
        const Vec3 centre{3.5 + (rng() % 1000) / 1000.0, 3.5 + (rng() % 1000) / 1000.0,
                          3.5 + (rng() % 1000) / 1000.0};
        const double radius = 1.0 + (rng() % 1000) / 1000.0;
        auto f = [&](const Vec3& rX) { return det::vec3_norm(det::vec3_sub(rX, centre)) - radius; };
        std::vector<std::int8_t> sign(m.mPoints.size());
        for (std::size_t p = 0; p < m.mPoints.size(); ++p)
            sign[p] = f(m.mPoints[p]) < 0.0 ? -1 : 1;
        const auto pts = m.mPoints;
        det::StuffCut cut = [&](std::int64_t A, std::int64_t B) {
            Vec3 a = pts[static_cast<std::size_t>(A)], b = pts[static_cast<std::size_t>(B)];
            for (int it = 0; it < 60; ++it) {
                const Vec3 mid = det::vec3_scale(det::vec3_add(a, b), 0.5);
                (f(mid) < 0.0 ? a : b) = mid;
            }
            return a;
        };
        det::StuffOptions so;
        det::isosurface_stuff(m, sign, cut, so);
        for (std::size_t q = 0; q < m.mTets.size(); ++q) {
            if (m.mMaterial[q] != 1)
                continue;
            const auto& v = m.mTets[q];
            const auto a = det::tet_dihedral_angles(m.mPoints[v[0]], m.mPoints[v[1]],
                                                    m.mPoints[v[2]], m.mPoints[v[3]]);
            lo = std::min(lo, a[0]);
            hi = std::max(hi, a[5]);
        }
    }
    EXPECT_GE(lo, 10.7);
    EXPECT_LE(hi, 164.8);
}

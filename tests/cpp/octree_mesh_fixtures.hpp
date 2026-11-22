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
// Checks and fixtures shared by the octree_mesh tests (`test_octree_mesh.cpp`, the full
// ones, and `test_quick_octree_mesh.cpp`, their small versions that run in the pull-request
// coverage job): the cell boxes and regions of the castellated output, the validity checks
// of the tetrahedral and hexahedral outputs, the option presets and the nested-sphere skins.

#pragma once

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "meshioplusplus/detail/surface_distance.hpp"
#include "meshioplusplus/operations/convert_cells.hpp"
#include "meshioplusplus/operations/merge.hpp"
#include "meshioplusplus/operations/octree_mesh.hpp"
#include "meshioplusplus/operations/primitives.hpp"
#include "meshioplusplus/operations/surface.hpp"
#include "../../src/cpp/src/detail/bcc_octree_tets.hpp"
#include "../../src/cpp/src/detail/hex_fit.hpp"
#include "../../src/cpp/src/detail/tri_bvh.hpp"

namespace omfix {

namespace mio = meshioplusplus;
using mio::detail::Vec3;

struct OmBox {
    Vec3 mLo, mHi;
};

inline std::vector<OmBox> om_boxes(const mio::Mesh& rM) {
    const double* p = rM.Points().As<double>();
    const auto block = rM.Cells(0);
    const std::int64_t* c = block.Conn().As<std::int64_t>();
    std::vector<OmBox> out(block.NumCells());
    for (std::size_t n = 0; n < out.size(); ++n) {
        out[n].mLo = {1e300, 1e300, 1e300};
        out[n].mHi = {-1e300, -1e300, -1e300};
        for (std::size_t v = 0; v < 8; ++v)
            for (std::size_t k = 0; k < 3; ++k) {
                const double x = p[static_cast<std::size_t>(c[n * 8 + v]) * 3 + k];
                out[n].mLo[k] = std::min(out[n].mLo[k], x);
                out[n].mHi[k] = std::max(out[n].mHi[k], x);
            }
    }
    return out;
}

inline double om_volume(const OmBox& rB) {
    return (rB.mHi[0] - rB.mLo[0]) * (rB.mHi[1] - rB.mLo[1]) * (rB.mHi[2] - rB.mLo[2]);
}

inline std::vector<std::int64_t> om_region(const mio::Mesh& rM, const char* pName,
                                           mio::RegionKind Kind) {
    const std::size_t r = rM.FindRegion(pName, Kind);
    EXPECT_NE(r, mio::Mesh::npos) << pName;
    const auto& e = rM.Region(r).mEntries;
    return std::vector<std::int64_t>(e.As<std::int64_t>(), e.As<std::int64_t>() + e.Size());
}

inline mio::OctreeMeshOptions om_castellated(double MaxSize, double MinSize) {
    mio::OctreeMeshOptions o;
    o.mOutput = mio::OctreeMeshOutput::Castellated;
    o.mMaxSize = MaxSize;
    o.mMinSize = MinSize;
    o.mPaddingRelative = 0.05;
    return o;
}

// The checks every tetrahedral output must pass: positive volumes, every
// interior triangle shared by two tetrahedra with opposite orientations, and
// a boundary that is a closed 2-manifold (every boundary edge on exactly two
// boundary triangles).
struct OmTetCheck {
    double mVolume = 0.0;
    double mMinDihedral = 180.0, mMaxDihedral = 0.0;
    std::size_t mBoundaryFaces = 0;
};

inline OmTetCheck om_check_tets(const mio::Mesh& rM) {
    OmTetCheck out;
    const double* p = rM.Points().As<double>();
    const auto block = rM.Cells(0);
    EXPECT_EQ(block.Type(), "tetra");
    const std::int64_t* c = block.Conn().As<std::int64_t>();
    auto at = [&](std::int64_t v) { return Vec3{p[3 * v], p[3 * v + 1], p[3 * v + 2]}; };
    std::map<std::array<std::int64_t, 3>, int> faces;
    constexpr int kFace[4][3] = {{1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}};
    for (std::size_t t = 0; t < block.NumCells(); ++t) {
        const std::int64_t* v = c + 4 * t;
        const double vol = mio::detail::triple_product(mio::detail::vec3_sub(at(v[1]), at(v[0])),
                                                       mio::detail::vec3_sub(at(v[2]), at(v[0])),
                                                       mio::detail::vec3_sub(at(v[3]), at(v[0]))) /
                           6.0;
        EXPECT_GT(vol, 0.0) << "tet " << t;
        out.mVolume += vol;
        const auto a = mio::detail::tet_dihedral_angles(at(v[0]), at(v[1]), at(v[2]), at(v[3]));
        out.mMinDihedral = std::min(out.mMinDihedral, a[0]);
        out.mMaxDihedral = std::max(out.mMaxDihedral, a[5]);
        for (const auto& f : kFace) {
            std::array<std::int64_t, 3> k{v[f[0]], v[f[1]], v[f[2]]};
            int parity = 0;
            for (int x = 0; x < 3; ++x)
                for (int y = x + 1; y < 3; ++y)
                    parity ^= k[x] > k[y];
            std::sort(k.begin(), k.end());
            const int s = parity ? -1 : 1;
            auto [it, fresh] = faces.emplace(k, s);
            if (!fresh) {
                EXPECT_EQ(it->second, -s) << "a triangle seen twice the same way";
                it->second = 2;  // matched
            }
        }
    }
    std::map<std::array<std::int64_t, 2>, int> edges;
    for (const auto& [k, s] : faces) {
        if (s == 2)
            continue;
        ++out.mBoundaryFaces;
        for (int e = 0; e < 3; ++e)
            ++edges[{std::min(k[e], k[(e + 1) % 3]), std::max(k[e], k[(e + 1) % 3])}];
    }
    for (const auto& [e, n] : edges)
        EXPECT_EQ(n, 2) << "boundary edge " << e[0] << "-" << e[1] << " on " << n << " triangles";
    return out;
}

inline mio::OctreeMeshOptions om_tet(double MaxSize, double MinSize, bool Uniform) {
    mio::OctreeMeshOptions o;
    o.mOutput = mio::OctreeMeshOutput::Tetra;
    o.mMaxSize = MaxSize;
    o.mMinSize = MinSize;
    o.mPaddingRelative = 0.05;
    if (Uniform) {
        o.mSkinGrading = mio::OctreeMeshSkinGrading::Uniform;
        o.mSurfaceSize = 2.0 * MinSize;
    }
    return o;
}

// The volume a closed, outward surface encloses.
inline double om_enclosed(const mio::Mesh& rSkin) {
    const auto soup = mio::detail::build_triangle_soup(rSkin, "");
    double v = 0.0;
    for (std::size_t t = 0; t < soup.NumTriangles(); ++t)
        v += mio::detail::vec3_dot(
                 soup.mCorners[3 * t],
                 mio::detail::vec3_cross(soup.mCorners[3 * t + 1], soup.mCorners[3 * t + 2])) /
             6.0;
    return v;
}

struct OmHexCheck {
    double mVolume = 0.0;
    double mMinScaledJacobian = 1.0;
};

// Every cell a positive hexahedron; the boundary closed and edge-manifold.
inline OmHexCheck om_check_hexes(const mio::Mesh& rM) {
    OmHexCheck out;
    EXPECT_EQ(rM.NumCellBlocks(), 1u);
    const auto block = rM.Cells(0);
    EXPECT_EQ(block.Type(), "hexahedron");
    mio::detail::HexMesh hm;
    const double* p = rM.Points().As<double>();
    for (std::size_t v = 0; v < rM.NumPoints(); ++v)
        hm.mPoints.push_back({p[3 * v], p[3 * v + 1], p[3 * v + 2]});
    const std::int64_t* c = block.Conn().As<std::int64_t>();
    for (std::size_t h = 0; h < block.NumCells(); ++h) {
        std::array<std::int64_t, 8> hex{};
        for (std::size_t k = 0; k < 8; ++k)
            hex[k] = c[8 * h + k];
        hm.mHexes.push_back(hex);
        out.mMinScaledJacobian =
            std::min(out.mMinScaledJacobian, mio::detail::hex_min_scaled_jacobian(hm.mPoints, hex));
    }
    const mio::detail::HexBoundary b = mio::detail::hex_boundary(hm);
    EXPECT_EQ(b.mNumOpenEdges, 0);
    EXPECT_EQ(b.mNumNonManifoldEdges, 0);
    for (const auto& q : b.mQuads) {
        const Vec3& a = hm.mPoints[static_cast<std::size_t>(q[0])];
        out.mVolume += (mio::detail::triple_product(a, hm.mPoints[static_cast<std::size_t>(q[1])],
                                                    hm.mPoints[static_cast<std::size_t>(q[2])]) +
                        mio::detail::triple_product(a, hm.mPoints[static_cast<std::size_t>(q[2])],
                                                    hm.mPoints[static_cast<std::size_t>(q[3])])) /
                       6.0;
    }
    return out;
}

inline mio::OctreeMeshOptions om_hex(double SurfaceSize) {
    mio::OctreeMeshOptions o;
    o.mOutput = mio::OctreeMeshOutput::Hexa;
    o.mSurfaceSize = SurfaceSize;
    return o;
}

// The hexahedra's corners with a non-positive scaled Jacobian.
inline std::int64_t om_flat_corners(const mio::Mesh& rM) {
    std::vector<Vec3> pts;
    const double* p = rM.Points().As<double>();
    for (std::size_t v = 0; v < rM.NumPoints(); ++v)
        pts.push_back({p[3 * v], p[3 * v + 1], p[3 * v + 2]});
    const auto block = rM.Cells(0);
    const std::int64_t* c = block.Conn().As<std::int64_t>();
    std::int64_t flat = 0;
    for (std::size_t h = 0; h < block.NumCells(); ++h) {
        std::array<std::int64_t, 8> hex{};
        for (std::size_t k = 0; k < 8; ++k)
            hex[k] = c[8 * h + k];
        for (double j : mio::detail::hex_corner_scaled_jacobians(pts, hex))
            flat += j <= 0.0;
    }
    return flat;
}

// A unit box on the lattice's planes, its faces split into triangles.
inline mio::Mesh om_unit_box() {
    return mio::convert_cells(mio::extract_surface(mio::box({{4, 4, 4}}, {{1.0, 1.0, 1.0}})),
                              {mio::ConvertCellsMode::Simplexify, false})
        .mMesh;
}

// The volume of each material of a hexahedral output (each hexahedron as six
// tetrahedra around its 0-6 diagonal), indexed by material.
inline std::vector<double> om_material_volumes(const mio::Mesh& rM) {
    const double* p = rM.Points().As<double>();
    const std::int64_t* c = rM.Cells(0).Conn().As<std::int64_t>();
    const std::int64_t* mat = rM.CellData(mio::kOctreeMeshMaterialName, 0).As<std::int64_t>();
    std::vector<double> out;
    static constexpr int kTets[6][4] = {{0, 1, 2, 6}, {0, 2, 3, 6}, {0, 3, 7, 6},
                                        {0, 7, 4, 6}, {0, 4, 5, 6}, {0, 5, 1, 6}};
    for (std::size_t h = 0; h < rM.Cells(0).NumCells(); ++h) {
        const auto m = static_cast<std::size_t>(mat[h]);
        if (out.size() <= m)
            out.resize(m + 1, 0.0);
        auto at = [&](int K) {
            const std::int64_t v = c[8 * h + static_cast<std::size_t>(K)];
            return Vec3{p[3 * v], p[3 * v + 1], p[3 * v + 2]};
        };
        for (const auto& t : kTets) {
            const Vec3 o = at(t[0]);
            out[m] += mio::detail::triple_product(mio::detail::vec3_sub(at(t[1]), o),
                                                  mio::detail::vec3_sub(at(t[2]), o),
                                                  mio::detail::vec3_sub(at(t[3]), o)) /
                      6.0;
        }
    }
    return out;
}

inline mio::Mesh om_nested_spheres(const std::vector<double>& rRadii, int Level = 4) {
    std::vector<mio::Mesh> parts;
    for (double r : rRadii)
        parts.push_back(mio::icosa_sphere(r, Level));
    std::vector<const mio::Mesh*> ptrs;
    for (const auto& m : parts)
        ptrs.push_back(&m);
    return mio::merge(ptrs, mio::MergeOptions{}).mMesh;
}

// A box of `N` squares a side whose face x = 0.49 has a patch of points pushed
// out by up to `Amplitude`, by an integer hash (the same on every platform): a
// scan's noisy creases on a part with real edges.
inline mio::Mesh om_jittered_bump_box(int N, double Amplitude) {
    mio::Mesh box =
        mio::convert_cells(
            mio::extract_surface(mio::box({{N, N, N}}, {{1.0, 1.0, 1.0}}, {{-0.51, -0.53, -0.55}})),
            {mio::ConvertCellsMode::Simplexify, false})
            .mMesh;
    const std::size_t np = box.NumPoints();
    mio::NDArray pts = mio::NDArray::Uninit(mio::DType::Float64, {np, std::size_t{3}});
    const double* src = box.Points().As<double>();
    double* p = pts.As<double>();
    std::copy(src, src + 3 * np, p);
    std::uint64_t h = 0x9e3779b97f4a7c15ull;
    for (std::size_t v = 0; v < np; ++v) {
        if (!(std::abs(p[3 * v] - 0.49) < 1e-9 && std::abs(p[3 * v + 1] + 0.03) < 0.4 &&
              std::abs(p[3 * v + 2] + 0.05) < 0.4))
            continue;
        h ^= h >> 33;
        h *= 0xff51afd7ed558ccdull;
        h ^= h >> 33;
        p[3 * v] += Amplitude * static_cast<double>(h % 1000003u) / 1000003.0;
    }
    box.AssignPoints(std::move(pts));
    return box;
}

}  // namespace omfix

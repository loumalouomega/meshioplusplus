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
#include <array>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/operations/slice.hpp"

namespace {

using meshioplusplus::Mesh;
using meshioplusplus::NDArray;
using meshioplusplus::slice;
using meshioplusplus::SliceOptions;

using P3 = std::array<double, 3>;

// Newell normal of a polygon corner ring (twin of test_skin.cpp's helper).
P3 slice_test_newell(const std::vector<P3>& rRing) {
    P3 n = {0, 0, 0};
    const std::size_t k = rRing.size();
    for (std::size_t i = 0; i < k; ++i) {
        const P3& a = rRing[i];
        const P3& b = rRing[(i + 1) % k];
        n[0] += (a[1] - b[1]) * (a[2] + b[2]);
        n[1] += (a[2] - b[2]) * (a[0] + b[0]);
        n[2] += (a[0] - b[0]) * (a[1] + b[1]);
    }
    return n;
}

// A single unit cube [0,1]^3 as one hexahedron (VTK node order).
Mesh unit_hex() {
    return mt::make_mesh(
        {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}},
        "hexahedron", {{0, 1, 2, 3, 4, 5, 6, 7}});
}

// Total area of a triangle/quad surface mesh (fan area of each face).
double section_area(const Mesh& rMesh) {
    const double* p = rMesh.Points().As<double>();
    const std::size_t dim = rMesh.PointDim();
    double area = 0.0;
    for (const auto cb : rMesh.CellRange()) {
        const std::int64_t* conn = cb.Conn().As<std::int64_t>();
        const std::size_t npc = cb.NodesPerCell();
        for (std::size_t c = 0; c < cb.NumCells(); ++c) {
            std::vector<P3> ring;
            for (std::size_t k = 0; k < npc; ++k) {
                const std::size_t n = static_cast<std::size_t>(conn[c * npc + k]);
                P3 v = {0, 0, 0};
                for (std::size_t d = 0; d < dim; ++d)
                    v[d] = p[n * dim + d];
                ring.push_back(v);
            }
            const P3 nrm = slice_test_newell(ring);
            area += 0.5 * std::sqrt(nrm[0] * nrm[0] + nrm[1] * nrm[1] + nrm[2] * nrm[2]);
        }
    }
    return area;
}

}  // namespace

TEST(Slice, HexMidPlaneAreaAndWinding) {
    Mesh m = unit_hex();
    SliceOptions opt;
    opt.mOrigin = {0.0, 0.0, 0.5};
    opt.mNormal = {0.0, 0.0, 1.0};
    Mesh sec = slice(m, opt);

    // The mid-plane section of the unit cube is a unit square.
    EXPECT_NEAR(section_area(sec), 1.0, 1e-12);

    // Marching-tetrahedra winding invariant: every section face's Newell normal
    // points toward the +normal side (dot with (0,0,1) > 0).
    const double* p = sec.Points().As<double>();
    const std::size_t dim = sec.PointDim();
    ASSERT_GT(sec.NumCellBlocks(), 0u);
    for (const auto cb : sec.CellRange()) {
        const std::int64_t* conn = cb.Conn().As<std::int64_t>();
        const std::size_t npc = cb.NodesPerCell();
        for (std::size_t c = 0; c < cb.NumCells(); ++c) {
            std::vector<P3> ring;
            for (std::size_t k = 0; k < npc; ++k) {
                const std::size_t n = static_cast<std::size_t>(conn[c * npc + k]);
                P3 v = {0, 0, 0};
                for (std::size_t d = 0; d < dim; ++d)
                    v[d] = p[n * dim + d];
                ring.push_back(v);
            }
            const P3 nrm = slice_test_newell(ring);
            EXPECT_GT(nrm[2], 1e-12) << "section face is not wound toward +normal";
            // The section lies exactly in the plane z = 0.5.
            for (const P3& v : ring)
                EXPECT_NEAR(v[2], 0.5, 1e-12);
        }
    }
}

TEST(Slice, PointDataInterpolatedExactly) {
    Mesh m = unit_hex();
    // f = z at the eight corners.
    NDArray f = NDArray::Uninit(meshioplusplus::DType::Float64, {8, 1});
    double* fd = f.As<double>();
    const double* p = m.Points().As<double>();
    for (std::size_t i = 0; i < 8; ++i)
        fd[i] = p[i * 3 + 2];
    m.AddPointData("f", std::move(f));

    SliceOptions opt;
    opt.mOrigin = {0.0, 0.0, 0.5};
    opt.mNormal = {0.0, 0.0, 1.0};
    Mesh sec = slice(m, opt);

    ASSERT_TRUE(sec.HasPointData("f"));
    const NDArray& out = sec.PointData("f");
    const double* od = out.As<double>();
    for (std::size_t i = 0; i < out.Size(); ++i)
        EXPECT_NEAR(od[i], 0.5, 1e-12);
}

TEST(Slice, MissingPlaneIsEmpty) {
    Mesh m = unit_hex();
    SliceOptions opt;
    opt.mOrigin = {0.0, 0.0, 5.0};
    opt.mNormal = {0.0, 0.0, 1.0};
    Mesh sec = slice(m, opt);
    EXPECT_EQ(sec.NumCellBlocks(), 0u);
    EXPECT_EQ(sec.NumPoints(), 0u);
}

TEST(Slice, TwoDMeshYieldsLines) {
    // Two unit quads in the xy-plane, cut by y = 0.5.
    Mesh m = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {2, 0, 0}, {2, 1, 0}},
                           "quad", {{0, 1, 2, 3}, {1, 4, 5, 2}});
    SliceOptions opt;
    opt.mOrigin = {0.0, 0.5, 0.0};
    opt.mNormal = {0.0, 1.0, 0.0};
    Mesh sec = slice(m, opt);

    ASSERT_EQ(sec.NumCellBlocks(), 1u);
    EXPECT_EQ(sec.Cells(0).Type(), "line");
    const double* p = sec.Points().As<double>();
    const std::size_t dim = sec.PointDim();
    const std::int64_t* conn = sec.Cells(0).Conn().As<std::int64_t>();
    double total = 0.0;
    for (std::size_t c = 0; c < sec.Cells(0).NumCells(); ++c) {
        const std::size_t a = static_cast<std::size_t>(conn[c * 2 + 0]);
        const std::size_t b = static_cast<std::size_t>(conn[c * 2 + 1]);
        double d2 = 0.0;
        for (std::size_t d = 0; d < dim; ++d) {
            const double dd = p[a * dim + d] - p[b * dim + d];
            d2 += dd * dd;
        }
        total += std::sqrt(d2);
    }
    EXPECT_NEAR(total, 2.0, 1e-12);
}

TEST(Slice, GrazingSharedFaceNotDoubleEmitted) {
    // Two hexes stacked in z sharing the internal face at z = 0.5; the plane
    // lies exactly on it. Emitted once (area 1.0), not twice (2.0).
    Mesh m = mt::make_mesh({{0, 0, 0},
                            {1, 0, 0},
                            {1, 1, 0},
                            {0, 1, 0},
                            {0, 0, 0.5},
                            {1, 0, 0.5},
                            {1, 1, 0.5},
                            {0, 1, 0.5},
                            {0, 0, 1},
                            {1, 0, 1},
                            {1, 1, 1},
                            {0, 1, 1}},
                           "hexahedron", {{0, 1, 2, 3, 4, 5, 6, 7}, {4, 5, 6, 7, 8, 9, 10, 11}});
    SliceOptions opt;
    opt.mOrigin = {0.0, 0.0, 0.5};
    opt.mNormal = {0.0, 0.0, 1.0};
    Mesh sec = slice(m, opt);
    EXPECT_NEAR(section_area(sec), 1.0, 1e-12);
}

TEST(Slice, ZeroNormalThrows) {
    Mesh m = unit_hex();
    SliceOptions opt;
    opt.mNormal = {0.0, 0.0, 0.0};
    EXPECT_THROW(slice(m, opt), std::invalid_argument);
}

TEST(Slice, ParentCellProvenance) {
    Mesh m = unit_hex();
    SliceOptions opt;
    opt.mOrigin = {0.0, 0.0, 0.5};
    opt.mNormal = {0.0, 0.0, 1.0};
    opt.mRecordParentIds = true;
    Mesh sec = slice(m, opt);
    ASSERT_TRUE(sec.HasCellData("slice:parent_cell"));
    // The single input hex is the parent of every section face.
    for (std::size_t b = 0; b < sec.CellDataNumBlocks("slice:parent_cell"); ++b) {
        const NDArray& a = sec.CellData("slice:parent_cell", b);
        const std::int64_t* d = a.As<std::int64_t>();
        for (std::size_t i = 0; i < a.Size(); ++i)
            EXPECT_EQ(d[i], 0);
    }
}

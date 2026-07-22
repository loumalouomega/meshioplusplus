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
#include <set>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/operations/isosurface.hpp"

namespace {

using meshioplusplus::DType;
using meshioplusplus::isosurface;
using meshioplusplus::IsosurfaceOptions;
using meshioplusplus::Mesh;
using meshioplusplus::NDArray;

using P3 = std::array<double, 3>;

// Newell normal of a polygon corner ring (twin of test_slice.cpp's helper).
P3 iso_test_newell(const std::vector<P3>& rRing) {
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

// The corner ring of one cell of one block, padded to three components.
std::vector<P3> iso_test_ring(const Mesh& rMesh, const std::int64_t* pConn, std::size_t Npc,
                              std::size_t Cell) {
    const double* p = rMesh.Points().As<double>();
    const std::size_t dim = rMesh.PointDim();
    std::vector<P3> ring;
    for (std::size_t k = 0; k < Npc; ++k) {
        const std::size_t n = static_cast<std::size_t>(pConn[Cell * Npc + k]);
        P3 v = {0, 0, 0};
        for (std::size_t d = 0; d < dim; ++d)
            v[d] = p[n * dim + d];
        ring.push_back(v);
    }
    return ring;
}

// Total area of a triangle/quad contour mesh (fan area of each face).
double contour_area(const Mesh& rMesh) {
    double area = 0.0;
    for (const auto cb : rMesh.CellRange()) {
        const std::int64_t* conn = cb.Conn().As<std::int64_t>();
        for (std::size_t c = 0; c < cb.NumCells(); ++c) {
            const P3 nrm = iso_test_newell(iso_test_ring(rMesh, conn, cb.NodesPerCell(), c));
            area += 0.5 * std::sqrt(nrm[0] * nrm[0] + nrm[1] * nrm[1] + nrm[2] * nrm[2]);
        }
    }
    return area;
}

// A single unit cube [0,1]^3 as one hexahedron (VTK node order) carrying the
// linear field f = x on its corners.
Mesh unit_hex_with_fx() {
    Mesh m = mt::make_mesh(
        {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}},
        "hexahedron", {{0, 1, 2, 3, 4, 5, 6, 7}});
    m.AddPointData("fx", mt::data_array({0, 1, 1, 0, 0, 1, 1, 0}));
    return m;
}

// Two hexes stacked in z sharing the internal face at z = 0.5, with a field
// equal to the isovalue exactly on that shared face.
Mesh stacked_hexes_with_fz() {
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
    m.AddPointData("fz", mt::data_array({0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0}));
    return m;
}

IsosurfaceOptions opts(const std::string& rArray, std::vector<double> Values) {
    IsosurfaceOptions o;
    o.mArrayName = rArray;
    o.mIsovalues = std::move(Values);
    return o;
}

}  // namespace

TEST(Isosurface, LinearFieldLevelSetIsPlanar) {
    Mesh iso = isosurface(unit_hex_with_fx(), opts("fx", {0.5}));

    // f = x is linear, so its 0.5 level set is exactly the plane x = 0.5, whose
    // cross-section of the unit cube is a unit square.
    EXPECT_NEAR(contour_area(iso), 1.0, 1e-12);
    ASSERT_GT(iso.NumCellBlocks(), 0u);
    const double* p = iso.Points().As<double>();
    for (std::size_t i = 0; i < iso.NumPoints(); ++i)
        EXPECT_NEAR(p[i * iso.PointDim()], 0.5, 1e-12);
}

TEST(Isosurface, ContouredFieldIsExactlyTheIsovalue) {
    Mesh iso = isosurface(unit_hex_with_fx(), opts("fx", {0.5}));
    ASSERT_TRUE(iso.HasPointData("fx"));
    const NDArray& a = iso.PointData("fx");
    ASSERT_GT(a.Size(), 0u);
    const double* v = a.As<double>();
    for (std::size_t i = 0; i < a.Size(); ++i)
        EXPECT_EQ(v[i], 0.5) << "the contoured field must be exact, not merely close";
}

TEST(Isosurface, FacesAreWoundTowardIncreasingField) {
    Mesh iso = isosurface(unit_hex_with_fx(), opts("fx", {0.5}));
    ASSERT_GT(iso.NumCellBlocks(), 0u);
    for (const auto cb : iso.CellRange()) {
        const std::int64_t* conn = cb.Conn().As<std::int64_t>();
        for (std::size_t c = 0; c < cb.NumCells(); ++c) {
            const P3 nrm = iso_test_newell(iso_test_ring(iso, conn, cb.NodesPerCell(), c));
            // f = x increases along +x, so every contour normal points that way.
            EXPECT_GT(nrm[0], 1e-12) << "contour face is not wound toward increasing field";
        }
    }
}

TEST(Isosurface, MultipleIsovaluesAreTaggedAscending) {
    Mesh iso = isosurface(unit_hex_with_fx(), opts("fx", {0.75, 0.25}));
    ASSERT_TRUE(iso.HasCellData("iso:value"));
    ASSERT_TRUE(iso.HasCellData("iso:index"));
    ASSERT_EQ(iso.CellDataNumBlocks("iso:value"), iso.NumCellBlocks());

    std::set<std::pair<std::int64_t, double>> tags;
    for (std::size_t b = 0; b < iso.NumCellBlocks(); ++b) {
        const NDArray& v = iso.CellData("iso:value", b);
        const NDArray& k = iso.CellData("iso:index", b);
        ASSERT_EQ(v.Dtype(), DType::Float64);
        ASSERT_EQ(k.Dtype(), DType::Int64);
        for (std::size_t c = 0; c < v.Size(); ++c)
            tags.emplace(k.As<std::int64_t>()[c], v.As<double>()[c]);
    }
    // Sorted ascending on the way in, so index 0 is the smaller value.
    const std::set<std::pair<std::int64_t, double>> want = {{0, 0.25}, {1, 0.75}};
    EXPECT_EQ(tags, want);
    // Two parallel unit squares.
    EXPECT_NEAR(contour_area(iso), 2.0, 1e-12);
}

TEST(Isosurface, DuplicateIsovaluesAreCutOnce) {
    EXPECT_NEAR(contour_area(isosurface(unit_hex_with_fx(), opts("fx", {0.5, 0.5}))), 1.0, 1e-12);
}

TEST(Isosurface, PlateauAtTheIsovalueIsEmittedOnce) {
    // `d >= 0` is the positive side, so the sign mask is total and the shared
    // face lying exactly at the isovalue is emitted once (area 1), not twice.
    Mesh iso = isosurface(stacked_hexes_with_fz(), opts("fz", {0.5}));
    EXPECT_NEAR(contour_area(iso), 1.0, 1e-12);
}

TEST(Isosurface, OutOfRangeIsovalueIsEmptyNotAnError) {
    Mesh iso = isosurface(unit_hex_with_fx(), opts("fx", {9.0}));
    EXPECT_EQ(iso.NumCellBlocks(), 0u);
    EXPECT_EQ(iso.NumPoints(), 0u);
}

TEST(Isosurface, ParentCellProvenance) {
    IsosurfaceOptions o = opts("fz", {0.25});
    o.mRecordParentIds = true;
    Mesh iso = isosurface(stacked_hexes_with_fz(), o);
    ASSERT_TRUE(iso.HasCellData("iso:parent_cell"));
    // fz = 0.25 straddles the lower hex only.
    for (std::size_t b = 0; b < iso.NumCellBlocks(); ++b) {
        const NDArray& a = iso.CellData("iso:parent_cell", b);
        ASSERT_EQ(a.Dtype(), DType::Int64);
        for (std::size_t c = 0; c < a.Size(); ++c)
            EXPECT_EQ(a.As<std::int64_t>()[c], 0);
    }
}

TEST(Isosurface, CellDataFieldThrowsByName) {
    Mesh m = unit_hex_with_fx();
    std::vector<NDArray> blocks;
    blocks.push_back(mt::data_array({7.0}));
    m.AddCellData("mat", std::move(blocks));
    try {
        isosurface(m, opts("mat", {1.0}));
        FAIL() << "contouring cell_data must throw";
    } catch (const std::invalid_argument& e) {
        const std::string msg = e.what();
        EXPECT_NE(msg.find("mat"), std::string::npos);
        // The message must name the conversion that makes the field contourable.
        EXPECT_NE(msg.find("to-point"), std::string::npos);
    }
}

TEST(Isosurface, BadArgumentsThrow) {
    Mesh m = unit_hex_with_fx();
    EXPECT_THROW(isosurface(m, opts("nope", {0.5})), std::invalid_argument);
    EXPECT_THROW(isosurface(m, opts("fx", {})), std::invalid_argument);
    EXPECT_THROW(isosurface(m, opts("", {0.5})), std::invalid_argument);
    EXPECT_THROW(isosurface(m, opts("fx", {std::nan("")})), std::invalid_argument);
    IsosurfaceOptions o = opts("fx", {0.5});
    o.mComponent = 7;
    EXPECT_THROW(isosurface(m, o), std::invalid_argument);
}

TEST(Isosurface, TwoDMeshYieldsLines) {
    // A single unit quad in the xy-plane with h = x: the 0.5 level set is one
    // segment of unit length.
    Mesh m = mt::make_mesh({{0, 0}, {1, 0}, {1, 1}, {0, 1}}, "quad", {{0, 1, 2, 3}});
    m.AddPointData("h", mt::data_array({0, 1, 1, 0}));
    Mesh iso = isosurface(m, opts("h", {0.5}));
    ASSERT_EQ(iso.NumCellBlocks(), 1u);
    EXPECT_EQ(std::string(iso.Cells(0).Type()), "line");

    double total = 0.0;
    const double* p = iso.Points().As<double>();
    const std::size_t dim = iso.PointDim();
    const auto cb = iso.Cells(0);
    const std::int64_t* conn = cb.Conn().As<std::int64_t>();
    for (std::size_t c = 0; c < cb.NumCells(); ++c) {
        const std::size_t a = static_cast<std::size_t>(conn[c * 2]);
        const std::size_t b = static_cast<std::size_t>(conn[c * 2 + 1]);
        double s = 0.0;
        for (std::size_t k = 0; k < dim; ++k) {
            const double d = p[a * dim + k] - p[b * dim + k];
            s += d * d;
        }
        total += std::sqrt(s);
    }
    EXPECT_NEAR(total, 1.0, 1e-12);
}

TEST(Isosurface, DeterministicAcrossRuns) {
    Mesh a = isosurface(unit_hex_with_fx(), opts("fx", {0.25, 0.6}));
    Mesh b = isosurface(unit_hex_with_fx(), opts("fx", {0.25, 0.6}));
    ASSERT_EQ(a.NumPoints(), b.NumPoints());
    ASSERT_EQ(a.NumCellBlocks(), b.NumCellBlocks());
    const double* pa = a.Points().As<double>();
    const double* pb = b.Points().As<double>();
    for (std::size_t i = 0; i < a.NumPoints() * a.PointDim(); ++i)
        EXPECT_EQ(pa[i], pb[i]);
    for (std::size_t k = 0; k < a.NumCellBlocks(); ++k) {
        const auto ca = a.Cells(k);
        const auto cbk = b.Cells(k);
        ASSERT_EQ(ca.NumCells(), cbk.NumCells());
        const std::int64_t* xa = ca.Conn().As<std::int64_t>();
        const std::int64_t* xb = cbk.Conn().As<std::int64_t>();
        for (std::size_t i = 0; i < ca.NumCells() * ca.NodesPerCell(); ++i)
            EXPECT_EQ(xa[i], xb[i]);
    }
}

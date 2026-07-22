// SPDX-License-Identifier: MIT
/// @file test_view_payload.cpp
/// @brief Invariants of the mesh -> Polyscope mapping.
///
/// These deliberately mirror the assertions in `tests/python/test_viewer.py`, which is
/// what makes the C++ and Python payload builders a *parity* pair rather than
/// two unrelated implementations that happen to coexist. When one changes,
/// this file and that one should change together.
///
/// Nothing here needs OpenGL, a window, or Polyscope itself: `view_payload.cpp`
/// is deliberately free of Polyscope headers, so these run in the default
/// build.

#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "mesh_fixtures.hpp"
#include "../../src/cpp/cli/view_payload.hpp"

using meshioplusplus::Mesh;
using meshioplusplus::NDArray;
using meshioplusplus::cli::build_view_payload;
using meshioplusplus::cli::kPadIndex;
using meshioplusplus::cli::QuantityKind;
using meshioplusplus::cli::QuantityOn;
using meshioplusplus::cli::ViewKind;

namespace {

/// The one quantity with this name, or nullptr.
const meshioplusplus::cli::Quantity* find_quantity(
    const meshioplusplus::cli::ViewPayload& rPayload, const std::string& rName) {
    for (const auto& q : rPayload.mQuantities)
        if (q.mName == rName)
            return &q;
    return nullptr;
}

/// A cell-data array over one block, as Float64.
NDArray scalar_block(const std::vector<double>& rValues) {
    NDArray a(meshioplusplus::DType::Float64, {rValues.size()});
    auto* p = reinterpret_cast<double*>(a.Data());
    for (std::size_t i = 0; i < rValues.size(); ++i)
        p[i] = rValues[i];
    return a;
}

}  // namespace

// --- kind routing ---------------------------------------------------------

TEST(ViewPayload, SurfaceMeshRoutesToSurface) {
    const auto payload = build_view_payload(mt::tri_mesh());
    EXPECT_EQ(payload.mKind, ViewKind::Surface);
    EXPECT_EQ(payload.mFaces.size(), 2u);
    EXPECT_TRUE(payload.mTets.empty());
    EXPECT_TRUE(payload.mMixedCells.empty());
}

TEST(ViewPayload, TetMeshUsesTheTetFastPath) {
    const auto payload = build_view_payload(mt::tet_mesh());
    EXPECT_EQ(payload.mKind, ViewKind::Volume);
    EXPECT_FALSE(payload.mTets.empty());
    EXPECT_TRUE(payload.mMixedCells.empty());
}

TEST(ViewPayload, HexMeshKeepsItsHexesAndSaysNothing) {
    // Decomposing a hexahedral mesh would destroy exactly the element structure
    // a user opens a viewer to look at, so it must not happen silently -- or
    // at all.
    const auto payload = build_view_payload(mt::hex_mesh());
    EXPECT_EQ(payload.mKind, ViewKind::Volume);
    EXPECT_EQ(payload.mHexes.size(), 1u);
    EXPECT_TRUE(payload.mTets.empty());
    EXPECT_TRUE(payload.mNotes.empty());
}

TEST(ViewPayload, VolumeKindOnASurfaceMeshThrows) {
    EXPECT_THROW(build_view_payload(mt::tri_mesh(), ViewKind::Volume),
                 std::invalid_argument);
}

// --- geometry -------------------------------------------------------------

TEST(ViewPayload, WedgeIsSimplexifiedWithANote) {
    const auto payload = build_view_payload(mt::wedge_mesh());
    EXPECT_EQ(payload.mKind, ViewKind::Volume);
    // A wedge splits into three tetrahedra.
    EXPECT_EQ(payload.mTets.size(), 3u);
    ASSERT_FALSE(payload.mNotes.empty());
    EXPECT_NE(payload.mNotes[0].find("split into tetrahedra"), std::string::npos);
}

TEST(ViewPayload, MixedTetAndHexArePaddedAndKeepBlockOrder) {
    // The row order must be *ours*: Polyscope's own registerTetHexMesh
    // concatenates tets-then-hexes regardless of input order, which would
    // silently re-map every per-cell quantity.
    Mesh mesh = mt::hex_mesh();
    mesh.AddCellBlock("tetra", mt::conn_from({{4, 5, 6, 0}}));

    const auto payload = build_view_payload(mesh);
    ASSERT_EQ(payload.mMixedCells.size(), 2u);
    EXPECT_TRUE(payload.mTets.empty());
    EXPECT_TRUE(payload.mHexes.empty());

    // Block order: the hexahedron first, then the tetrahedron.
    for (std::size_t k = 0; k < 8; ++k)
        EXPECT_NE(payload.mMixedCells[0][k], kPadIndex) << "hexahedron slot " << k;
    for (std::size_t k = 4; k < 8; ++k)
        EXPECT_EQ(payload.mMixedCells[1][k], kPadIndex) << "tetra pad slot " << k;

    // Both types are natively renderable, so nothing was converted.
    EXPECT_TRUE(payload.mNotes.empty());
}

TEST(ViewPayload, SurfaceKindExtractsTheBoundaryOfASolid) {
    const auto payload = build_view_payload(mt::hex_mesh(), ViewKind::Surface);
    EXPECT_EQ(payload.mKind, ViewKind::Surface);
    ASSERT_EQ(payload.mFaces.size(), 6u);
    for (const auto& face : payload.mFaces)
        EXPECT_EQ(face.size(), 4u);
}

// --- cell data alignment (the reason the global cell index exists) --------

TEST(ViewPayload, CellDataFollowsBlockOrder) {
    Mesh mesh = mt::tri_quad_mesh();
    std::vector<NDArray> blocks;
    // One value per block, distinct so a mis-ordering is visible.
    double next = 10.0;
    for (const auto& block : mesh.CellRange()) {
        std::vector<double> values(block.NumCells());
        for (auto& v : values)
            v = next++;
        blocks.push_back(scalar_block(values));
    }
    mesh.AddCellData("tag", std::move(blocks));

    const auto payload = build_view_payload(mesh);
    const auto* tag = find_quantity(payload, "tag");
    ASSERT_NE(tag, nullptr);
    ASSERT_EQ(tag->mValues.size(), payload.mFaces.size());
    EXPECT_EQ(tag->mOn, QuantityOn::Faces);
    // Values must come out in block order, ascending, exactly as assigned.
    for (std::size_t i = 1; i < tag->mValues.size(); ++i)
        EXPECT_LT(tag->mValues[i - 1], tag->mValues[i]);
}

TEST(ViewPayload, ASkippedLowerDimensionalBlockDoesNotOffsetCellData) {
    // A gmsh-style boundary-marker `line` block sits before the triangles. A
    // naive per-dimension concatenation would read the line block's values as
    // the triangles' and colour the whole mesh wrong -- while still rendering.
    Mesh combined = mt::tri_mesh();
    combined.AddCellBlock("line", mt::conn_from({{0, 1}, {1, 2}}));

    std::vector<NDArray> blocks;
    blocks.push_back(scalar_block({7.0, 8.0}));    // the triangles
    blocks.push_back(scalar_block({-1.0, -2.0}));  // the lines
    combined.AddCellData("tag", std::move(blocks));

    const auto payload = build_view_payload(combined);
    ASSERT_EQ(payload.mFaces.size(), 2u);
    const auto* tag = find_quantity(payload, "tag");
    ASSERT_NE(tag, nullptr);
    ASSERT_EQ(tag->mValues.size(), 2u);
    EXPECT_DOUBLE_EQ(tag->mValues[0], 7.0);
    EXPECT_DOUBLE_EQ(tag->mValues[1], 8.0);
}

TEST(ViewPayload, SurfaceCellDataIsSampledFromTheOwningVolumeCell) {
    Mesh mesh = mt::hex_mesh();
    std::vector<NDArray> blocks;
    blocks.push_back(scalar_block({42.0}));
    mesh.AddCellData("material", std::move(blocks));

    const auto payload = build_view_payload(mesh, ViewKind::Surface);
    const auto* material = find_quantity(payload, "material");
    ASSERT_NE(material, nullptr);
    ASSERT_EQ(material->mValues.size(), 6u);
    for (double v : material->mValues)
        EXPECT_DOUBLE_EQ(v, 42.0);

    // Provenance is plumbing; it must never reach a colour-by menu.
    EXPECT_EQ(find_quantity(payload, "surface:parent_cell"), nullptr);
}

// --- quantity mapping rules ----------------------------------------------

TEST(ViewPayload, ScalarPointDataBecomesOneScalarQuantity) {
    Mesh mesh = mt::tri_mesh();
    NDArray t(meshioplusplus::DType::Float64, {mesh.NumPoints()});
    auto* p = reinterpret_cast<double*>(t.Data());
    for (std::size_t i = 0; i < mesh.NumPoints(); ++i)
        p[i] = static_cast<double>(i);
    mesh.AddPointData("t", std::move(t));

    const auto payload = build_view_payload(mesh);
    const auto* q = find_quantity(payload, "t");
    ASSERT_NE(q, nullptr);
    EXPECT_EQ(q->mKind, QuantityKind::Scalar);
    EXPECT_EQ(q->mOn, QuantityOn::Vertices);
    EXPECT_EQ(q->mValues.size(), mesh.NumPoints());
}

TEST(ViewPayload, ThreeWideBecomesAVectorPlusAMagnitude) {
    Mesh mesh = mt::tri_mesh();
    const std::size_t n = mesh.NumPoints();
    NDArray v(meshioplusplus::DType::Float64, {n, 3});
    auto* p = reinterpret_cast<double*>(v.Data());
    for (std::size_t i = 0; i < n; ++i) {
        p[i * 3] = 3.0;
        p[i * 3 + 1] = 4.0;
        p[i * 3 + 2] = 0.0;
    }
    mesh.AddPointData("disp", std::move(v));

    const auto payload = build_view_payload(mesh);
    const auto* vec = find_quantity(payload, "disp");
    const auto* mag = find_quantity(payload, "disp:magnitude");
    ASSERT_NE(vec, nullptr);
    ASSERT_NE(mag, nullptr);
    EXPECT_EQ(vec->mKind, QuantityKind::Vector);
    EXPECT_EQ(mag->mKind, QuantityKind::Scalar);
    EXPECT_DOUBLE_EQ(mag->mValues[0], 5.0);
}

TEST(ViewPayload, ColourIsDetectedByNameNotByRange) {
    // A normalized field in [0,1] must not silently render as RGB.
    auto with_named_array = [](const std::string& rName) {
        Mesh mesh = mt::tri_mesh();
        const std::size_t n = mesh.NumPoints();
        NDArray v(meshioplusplus::DType::Float64, {n, 3});
        auto* p = reinterpret_cast<double*>(v.Data());
        for (std::size_t i = 0; i < n * 3; ++i)
            p[i] = 0.25;
        mesh.AddPointData(rName, std::move(v));
        return build_view_payload(mesh);
    };

    const auto asVector = with_named_array("normalized_displacement");
    ASSERT_NE(find_quantity(asVector, "normalized_displacement"), nullptr);
    EXPECT_EQ(find_quantity(asVector, "normalized_displacement")->mKind,
              QuantityKind::Vector);

    const auto asColor = with_named_array("color");
    ASSERT_NE(find_quantity(asColor, "color"), nullptr);
    EXPECT_EQ(find_quantity(asColor, "color")->mKind, QuantityKind::Color);
}

TEST(ViewPayload, NothingLossyMeansNoNotes) {
    EXPECT_TRUE(build_view_payload(mt::tri_mesh()).mNotes.empty());
}

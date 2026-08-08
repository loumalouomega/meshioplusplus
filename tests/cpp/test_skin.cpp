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
#include <cstdint>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/cell_faces.hpp"
#include "meshioplusplus/skin.hpp"

namespace {

using meshioplusplus::cell_type_from_name;
using meshioplusplus::cell_type_num_nodes;
using meshioplusplus::CellType;
using meshioplusplus::detail::cell_faces;
using meshioplusplus::detail::CellFaceDef;

using P3 = std::array<double, 3>;

// ---- reference elements ----
//
// Corner coordinates of each supported volume type's reference element,
// plus the type's edge list (mid-edge nodes are generated as exact edge
// midpoints) and face-center definitions, so the invariant tests below can
// verify both the outward winding of every face row *and* the mid-node /
// face-center index assignment of the higher-order tables.

struct SkinTestRefElement {
    std::vector<P3> mPoints;
};

P3 skin_test_mid(const P3& rA, const P3& rB) {
    return {0.5 * (rA[0] + rB[0]), 0.5 * (rA[1] + rB[1]), 0.5 * (rA[2] + rB[2])};
}

P3 skin_test_avg(const std::vector<P3>& rPts, const std::vector<int>& rIds) {
    P3 c = {0, 0, 0};
    for (int i : rIds)
        for (int k = 0; k < 3; ++k)
            c[k] += rPts[i][k];
    for (int k = 0; k < 3; ++k)
        c[k] /= static_cast<double>(rIds.size());
    return c;
}

// Build the reference element of `type`: corners, then mid-edge nodes for
// each (a, b) in `edges`, then one node per id-list in `face_centers`.
SkinTestRefElement skin_test_ref(std::vector<P3> corners,
                                 const std::vector<std::array<int, 2>>& rEdges,
                                 const std::vector<std::vector<int>>& rFaceCenters = {}) {
    SkinTestRefElement e;
    e.mPoints = std::move(corners);
    const std::vector<P3> base = e.mPoints;
    for (const auto& ed : rEdges)
        e.mPoints.push_back(skin_test_mid(base[ed[0]], base[ed[1]]));
    for (const auto& fc : rFaceCenters)
        e.mPoints.push_back(skin_test_avg(base, fc));
    return e;
}

// The reference elements, keyed by meshio type name. Corner layouts follow
// this repo's conventions (see cell_faces.hpp): tetra/pyramid base normal
// toward the apex, hexahedron/wedge base normal toward the top face.
std::map<std::string, SkinTestRefElement> skin_test_reference_elements() {
    std::map<std::string, SkinTestRefElement> out;

    const std::vector<P3> tet = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    const std::vector<std::array<int, 2>> tet_edges = {{0, 1}, {1, 2}, {2, 0},
                                                       {0, 3}, {1, 3}, {2, 3}};
    out["tetra"] = skin_test_ref(tet, {});
    out["tetra10"] = skin_test_ref(tet, tet_edges);

    const std::vector<P3> hex = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                                 {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
    const std::vector<std::array<int, 2>> hex_edges = {{0, 1}, {1, 2}, {2, 3}, {3, 0},
                                                       {4, 5}, {5, 6}, {6, 7}, {7, 4},
                                                       {0, 4}, {1, 5}, {2, 6}, {3, 7}};
    // hexahedron27 nodes 20-25, in vtkTriQuadraticHexahedron::Faces order:
    // x-min, x-max, y-min, y-max, bottom, top. Corrected in v9.9.0 alongside
    // cell_faces.cpp's own table -- the geometric assertion in
    // SkinFaceTables.OutwardWindingAndNodeCounts below (mNodes[8] must sit at
    // the face's corner average) is only a real oracle if this agrees with the
    // real VTK convention rather than with whatever the table happens to say.
    const std::vector<std::vector<int>> hex_face_centers = {
        {0, 4, 7, 3}, {1, 2, 6, 5}, {0, 1, 5, 4}, {3, 7, 6, 2}, {0, 1, 2, 3}, {4, 5, 6, 7}};
    out["hexahedron"] = skin_test_ref(hex, {});
    out["hexahedron20"] = skin_test_ref(hex, hex_edges);
    // hexahedron27 additionally has the body center as node 26.
    out["hexahedron27"] = skin_test_ref(hex, hex_edges, hex_face_centers);
    out["hexahedron27"].mPoints.push_back(skin_test_avg(hex, {0, 1, 2, 3, 4, 5, 6, 7}));

    const std::vector<P3> wedge = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0},
                                   {0, 0, 1}, {1, 0, 1}, {0, 1, 1}};
    const std::vector<std::array<int, 2>> wedge_edges = {{0, 1}, {1, 2}, {2, 0}, {3, 4}, {4, 5},
                                                         {5, 3}, {0, 3}, {1, 4}, {2, 5}};
    // wedge18 additionally has a center node on each of the three quad faces
    // (nodes 15-17); the two triangle faces gain none.
    const std::vector<std::vector<int>> wedge_face_centers = {
        {0, 1, 4, 3}, {1, 2, 5, 4}, {2, 0, 3, 5}};
    out["wedge"] = skin_test_ref(wedge, {});
    out["wedge15"] = skin_test_ref(wedge, wedge_edges);
    out["wedge18"] = skin_test_ref(wedge, wedge_edges, wedge_face_centers);

    const std::vector<P3> pyr = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0.5, 0.5, 1}};
    const std::vector<std::array<int, 2>> pyr_edges = {{0, 1}, {1, 2}, {2, 3}, {3, 0},
                                                       {0, 4}, {1, 4}, {2, 4}, {3, 4}};
    out["pyramid"] = skin_test_ref(pyr, {});
    out["pyramid13"] = skin_test_ref(pyr, pyr_edges);
    out["pyramid14"] = skin_test_ref(pyr, pyr_edges, {{0, 1, 2, 3}});

    return out;
}

P3 skin_test_sub(const P3& rA, const P3& rB) {
    return {rA[0] - rB[0], rA[1] - rB[1], rA[2] - rB[2]};
}

// Newell normal of a polygon's corner ring.
P3 skin_test_newell(const std::vector<P3>& rRing) {
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

double skin_test_dot(const P3& rA, const P3& rB) {
    return rA[0] * rB[0] + rA[1] * rB[1] + rA[2] * rB[2];
}

// Count cells per type name across all blocks of a mesh.
std::map<std::string, std::size_t> skin_test_cell_counts(const mt::Mesh& rMesh) {
    std::map<std::string, std::size_t> counts;
    for (const auto cb : rMesh.CellRange())
        counts[cb.Type()] += cb.NumCells();
    return counts;
}

// A unit cube split into 6 pyramids whose shared apex (node 8) is the cube
// center: the skin is the 6 cube faces, and compaction must drop the center.
mt::Mesh skin_test_pyramid_cube() {
    std::vector<std::vector<double>> pts = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0},
                                            {0, 1, 0}, {0, 0, 1}, {1, 0, 1},
                                            {1, 1, 1}, {0, 1, 1}, {0.5, 0.5, 0.5}};
    // Each base is a cube face wound so its normal points toward the apex.
    std::vector<std::vector<std::int64_t>> cells = {
        {0, 1, 2, 3, 8},  // bottom (z=0), normal +z
        {7, 6, 5, 4, 8},  // top (z=1), normal -z
        {0, 4, 5, 1, 8},  // front (y=0), normal +y
        {2, 6, 7, 3, 8},  // back (y=1), normal -y
        {0, 3, 7, 4, 8},  // left (x=0), normal +x
        {1, 5, 6, 2, 8},  // right (x=1), normal -x
    };
    return mt::make_mesh(pts, "pyramid", cells);
}

mt::Mesh skin_test_pyramid_mesh() {
    return mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0.5, 0.5, 1}}, "pyramid",
                         {{0, 1, 2, 3, 4}});
}

}  // namespace

// ---- face-table invariants ----

TEST(SkinFaceTables, OutwardWindingAndNodeCounts) {
    const auto refs = skin_test_reference_elements();
    for (const auto& [name, ref] : refs) {
        const CellType ct = cell_type_from_name(name);
        ASSERT_NE(ct, CellType::Custom) << name;
        const std::vector<CellFaceDef>& faces = cell_faces(ct);
        ASSERT_FALSE(faces.empty()) << name;
        ASSERT_EQ(static_cast<int>(ref.mPoints.size()), cell_type_num_nodes(ct)) << name;

        // Cell centroid from the corner nodes, which are the leading nodes
        // of the reference element: 4 (tetra), 8 (hexa), 6 (wedge), 5 (pyramid).
        std::size_t nc = 0;
        switch (ct) {
            case CellType::Tetra:
            case CellType::Tetra10:
                nc = 4;
                break;
            case CellType::Hexahedron:
            case CellType::Hexahedron20:
            case CellType::Hexahedron27:
                nc = 8;
                break;
            case CellType::Wedge:
            case CellType::Wedge15:
            case CellType::Wedge18:
                nc = 6;
                break;
            default:
                nc = 5;  // pyramids
                break;
        }
        P3 centroid = {0, 0, 0};
        for (std::size_t i = 0; i < nc; ++i)
            for (int k = 0; k < 3; ++k)
                centroid[k] += ref.mPoints[i][k];
        for (int k = 0; k < 3; ++k)
            centroid[k] /= static_cast<double>(nc);

        for (const CellFaceDef& fd : faces) {
            // Node-count consistency with the face type itself.
            EXPECT_EQ(static_cast<int>(fd.mNumNodes), cell_type_num_nodes(fd.mFaceType)) << name;
            EXPECT_TRUE(fd.mNumCorners == 3 || fd.mNumCorners == 4) << name;
            for (std::uint8_t k = 0; k < fd.mNumNodes; ++k)
                ASSERT_LT(fd.mNodes[k], ref.mPoints.size()) << name;

            // Outward winding: Newell normal of the corner ring must point
            // away from the cell centroid.
            std::vector<P3> ring;
            for (std::uint8_t k = 0; k < fd.mNumCorners; ++k)
                ring.push_back(ref.mPoints[fd.mNodes[k]]);
            const P3 normal = skin_test_newell(ring);
            P3 face_centroid = {0, 0, 0};
            for (const P3& p : ring)
                for (int k = 0; k < 3; ++k)
                    face_centroid[k] += p[k];
            for (int k = 0; k < 3; ++k)
                face_centroid[k] /= static_cast<double>(ring.size());
            EXPECT_GT(skin_test_dot(normal, skin_test_sub(face_centroid, centroid)), 0.0)
                << name << " face is not wound outward";

            // Mid-edge nodes: node mNumCorners+k must sit exactly at the
            // midpoint of corners k and k+1 on the reference element.
            const std::uint8_t num_mid =
                fd.mNumNodes > fd.mNumCorners ? fd.mNumCorners : std::uint8_t{0};
            for (std::uint8_t k = 0; k < num_mid; ++k) {
                const P3& a = ref.mPoints[fd.mNodes[k]];
                const P3& b = ref.mPoints[fd.mNodes[(k + 1) % fd.mNumCorners]];
                const P3 expect = skin_test_mid(a, b);
                const P3& got = ref.mPoints[fd.mNodes[fd.mNumCorners + k]];
                for (int c = 0; c < 3; ++c)
                    EXPECT_DOUBLE_EQ(got[c], expect[c])
                        << name << " mid-edge node mismatch on face";
            }

            // Face-center node (quad9): must be the corner average.
            if (fd.mNumNodes == 9) {
                std::vector<P3> corner_ring;
                for (std::uint8_t k = 0; k < 4; ++k)
                    corner_ring.push_back(ref.mPoints[fd.mNodes[k]]);
                P3 expect = {0, 0, 0};
                for (const P3& p : corner_ring)
                    for (int c = 0; c < 3; ++c)
                        expect[c] += 0.25 * p[c];
                const P3& got = ref.mPoints[fd.mNodes[8]];
                for (int c = 0; c < 3; ++c)
                    EXPECT_DOUBLE_EQ(got[c], expect[c]) << name << " face-center node mismatch";
            }
        }
    }
}

// ---- extraction ----

TEST(Skin, TetMesh) {
    // 2 tets sharing one interior face: 8 faces total, 6 on the boundary.
    const mt::Mesh skin = meshioplusplus::extract_skin(mt::tet_mesh());
    const auto counts = skin_test_cell_counts(skin);
    ASSERT_EQ(counts.size(), 1u);
    EXPECT_EQ(counts.at("triangle"), 6u);
    EXPECT_EQ(skin.NumPoints(), 5u);  // every node is on the boundary
}

TEST(Skin, HexMesh) {
    const mt::Mesh skin = meshioplusplus::extract_skin(mt::hex_mesh());
    const auto counts = skin_test_cell_counts(skin);
    ASSERT_EQ(counts.size(), 1u);
    EXPECT_EQ(counts.at("quad"), 6u);
    EXPECT_EQ(skin.NumPoints(), 8u);
}

TEST(Skin, WedgeMesh) {
    const mt::Mesh skin = meshioplusplus::extract_skin(mt::wedge_mesh());
    const auto counts = skin_test_cell_counts(skin);
    EXPECT_EQ(counts.at("triangle"), 2u);
    EXPECT_EQ(counts.at("quad"), 3u);
}

// wedge18 became skinnable when its row was added to cell_faces.cpp (v7.5.0);
// before that it warned and was skipped.
TEST(Skin, Wedge18Mesh) {
    const std::vector<std::vector<double>> pts = {
        {0, 0, 0},   {1, 0, 0},     {0, 1, 0},   {0, 0, 1},     {1, 0, 1},       {0, 1, 1},
        {0.5, 0, 0}, {0.5, 0.5, 0}, {0, 0.5, 0}, {0.5, 0, 1},   {0.5, 0.5, 1},   {0, 0.5, 1},
        {0, 0, 0.5}, {1, 0, 0.5},   {0, 1, 0.5}, {0.5, 0, 0.5}, {0.5, 0.5, 0.5}, {0, 0.5, 0.5}};
    const mt::Mesh skin = meshioplusplus::extract_skin(mt::make_mesh(
        pts, "wedge18", {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17}}));
    const auto counts = skin_test_cell_counts(skin);
    EXPECT_EQ(counts.at("triangle6"), 2u);
    EXPECT_EQ(counts.at("quad9"), 3u);
}

TEST(Skin, PyramidMesh) {
    const mt::Mesh skin = meshioplusplus::extract_skin(skin_test_pyramid_mesh());
    const auto counts = skin_test_cell_counts(skin);
    EXPECT_EQ(counts.at("triangle"), 4u);
    EXPECT_EQ(counts.at("quad"), 1u);
}

TEST(Skin, Tet10Mesh) {
    const mt::Mesh skin = meshioplusplus::extract_skin(mt::tet10_mesh());
    const auto counts = skin_test_cell_counts(skin);
    ASSERT_EQ(counts.size(), 1u);
    EXPECT_EQ(counts.at("triangle6"), 4u);
    EXPECT_EQ(skin.NumPoints(), 10u);
}

TEST(Skin, Hex20Mesh) {
    const mt::Mesh skin = meshioplusplus::extract_skin(mt::hex20_mesh());
    const auto counts = skin_test_cell_counts(skin);
    ASSERT_EQ(counts.size(), 1u);
    EXPECT_EQ(counts.at("quad8"), 6u);
}

TEST(Skin, PyramidCubeCompaction) {
    // 6 pyramids filling a cube: the skin is the 6 cube faces, and the
    // shared apex (the cube center, node 8) must be compacted away.
    const mt::Mesh skin = meshioplusplus::extract_skin(skin_test_pyramid_cube());
    const auto counts = skin_test_cell_counts(skin);
    ASSERT_EQ(counts.size(), 1u);
    EXPECT_EQ(counts.at("quad"), 6u);
    EXPECT_EQ(skin.NumPoints(), 8u);

    // Compaction keeps ascending original order: node i of the skin is
    // cube corner i.
    const meshioplusplus::NDArray& pts = skin.Points();
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(pts, 0 * 3 + 0), 0.0);
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(pts, 6 * 3 + 0), 1.0);
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(pts, 6 * 3 + 2), 1.0);
}

TEST(Skin, PointDataSubset) {
    mt::Mesh volume = skin_test_pyramid_cube();
    meshioplusplus::NDArray pd(meshioplusplus::DType::Float64, {9});
    for (int i = 0; i < 9; ++i)
        pd.As<double>()[i] = 10.0 * i;
    volume.AddPointData("tag", std::move(pd));

    const mt::Mesh skin = meshioplusplus::extract_skin(volume);
    ASSERT_TRUE(skin.HasPointData("tag"));
    const meshioplusplus::NDArray& out = skin.PointData("tag");
    ASSERT_EQ(out.Shape()[0], 8u);
    for (int i = 0; i < 8; ++i)
        EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(out, i), 10.0 * i);
}

TEST(Skin, Linearize) {
    // tetra10 -> corner-only triangles; the 6 mid-edge nodes compact away.
    const mt::Mesh skin = meshioplusplus::extract_skin(mt::tet10_mesh(), /*linearize=*/true);
    const auto counts = skin_test_cell_counts(skin);
    ASSERT_EQ(counts.size(), 1u);
    EXPECT_EQ(counts.at("triangle"), 4u);
    EXPECT_EQ(skin.NumPoints(), 4u);
}

TEST(Skin, SurfaceOnlyThrows) {
    EXPECT_THROW(meshioplusplus::extract_skin(mt::tri_mesh()), std::invalid_argument);
    EXPECT_FALSE(meshioplusplus::has_skinnable_cells(mt::tri_mesh()));
    EXPECT_TRUE(meshioplusplus::has_skinnable_cells(mt::tet_mesh()));
}

TEST(Skin, HasSkinnableCellsRecognizesAPolyhedronOnlyMesh) {
    // A polyhedron block carries its own faces and IS skinnable (extract_skin
    // supports it since v9.16.0) -- has_skinnable_cells must agree, or every
    // consumer gated on it (convertSurfaceOps/convert_surface, and the
    // STL/PLY/SVG/TikZ writers' skin=true default) silently takes the
    // unskinned fallback for a mesh with only a polyhedron block.
    mt::Mesh m;
    m.AssignPoints(mt::points_from({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}}));
    m.AddPolyhedronBlock("polyhedron4", {{{0, 2, 1}, {0, 1, 3}, {1, 2, 3}, {2, 0, 3}}});
    EXPECT_TRUE(meshioplusplus::has_skinnable_cells(m));
    EXPECT_NO_THROW(meshioplusplus::extract_skin(m));
}

TEST(Skin, SurfaceBlocksIgnored) {
    // A mesh with both a tetra block and a pre-existing triangle block: the
    // skin comes from the volume cells only.
    mt::Mesh m = mt::tet_mesh();
    m.AddCellBlock("triangle", mt::conn_from({{0, 1, 2}}));
    const mt::Mesh skin = meshioplusplus::extract_skin(m);
    const auto counts = skin_test_cell_counts(skin);
    ASSERT_EQ(counts.size(), 1u);
    EXPECT_EQ(counts.at("triangle"), 6u);
}

TEST(Skin, TwoHexStack) {
    // Two stacked unit cubes: 12 faces total, 2 interior -> 10 boundary quads.
    std::vector<std::vector<double>> pts;
    for (int layer = 0; layer < 3; ++layer)
        for (auto& xy : std::vector<std::array<double, 2>>{{0, 0}, {1, 0}, {1, 1}, {0, 1}})
            pts.push_back({xy[0], xy[1], static_cast<double>(layer)});
    const mt::Mesh volume =
        mt::make_mesh(pts, "hexahedron", {{0, 1, 2, 3, 4, 5, 6, 7}, {4, 5, 6, 7, 8, 9, 10, 11}});
    const mt::Mesh skin = meshioplusplus::extract_skin(volume);
    const auto counts = skin_test_cell_counts(skin);
    ASSERT_EQ(counts.size(), 1u);
    EXPECT_EQ(counts.at("quad"), 10u);
    EXPECT_EQ(skin.NumPoints(), 12u);
}

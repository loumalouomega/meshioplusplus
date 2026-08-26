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

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/cell_subdivision.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/gid.hpp"

#ifdef MESHIOPLUSPLUS_HAS_GIDPOST

#include <algorithm>
#include <array>
#include <fstream>
#include <map>
#include <sstream>
#include <utility>
#include <string>
#include <vector>

using meshioplusplus::CellType;
using meshioplusplus::DType;
using meshioplusplus::Mesh;
using meshioplusplus::NDArray;

namespace {

// A minimal, independent parser of the ASCII `.post.msh` gidpost writes --
// deliberately NOT going through any meshio++ reading code (there is none:
// gidpost has no read functions). This is what makes it a bytes oracle
// rather than a round trip: it checks WHERE THE FILE SAYS a node is,
// geometrically, against an independently-derived edge/midpoint table.
struct ParsedMeshBlock {
    std::string mName;
    bool mHadCoordinatesRows = false;
    // Each row is [id, n1, n2, ...] (material column, if any, included --
    // no test here writes a material id, so rows are always plain).
    std::vector<std::vector<long long>> mElementRows;
};

struct ParsedGidFile {
    std::map<long long, std::array<double, 3>> mCoords;  // 1-based id -> xyz
    std::vector<ParsedMeshBlock> mMeshes;
};

ParsedGidFile gid_parse_ascii_msh(const std::string& rPath) {
    std::ifstream f(rPath);
    ParsedGidFile out;
    std::string line;
    bool in_coords = false;
    bool in_elems = false;
    while (std::getline(f, line)) {
        if (line.rfind("MESH ", 0) == 0) {
            out.mMeshes.emplace_back();
            const auto q1 = line.find('"');
            const auto q2 = line.find('"', q1 + 1);
            out.mMeshes.back().mName = line.substr(q1 + 1, q2 - q1 - 1);
            continue;
        }
        if (line == "Coordinates") {
            in_coords = true;
            continue;
        }
        if (line.rfind("End Coordinates", 0) == 0) {
            in_coords = false;
            continue;
        }
        if (line == "Elements") {
            in_elems = true;
            continue;
        }
        if (line.rfind("End Elements", 0) == 0) {
            in_elems = false;
            continue;
        }
        if (in_coords) {
            std::istringstream iss(line);
            long long id = 0;
            double x = 0, y = 0, z = 0;
            if (iss >> id >> x >> y >> z) {
                out.mCoords[id] = {x, y, z};
                if (!out.mMeshes.empty())
                    out.mMeshes.back().mHadCoordinatesRows = true;
            }
        } else if (in_elems && !out.mMeshes.empty()) {
            std::istringstream iss(line);
            std::vector<long long> row;
            long long v = 0;
            while (iss >> v)
                row.push_back(v);
            if (!row.empty())
                out.mMeshes.back().mElementRows.push_back(row);
        }
    }
    return out;
}

// Builds `nCorners` corners plus one point per edge of `Type` (per
// detail::cell_refine_edges -- meshio++'s OWN edge table, read here, never
// retyped), each at the true midpoint of its edge's two corners. Returns the
// full point list; corner k is at `corners[k]`.
std::vector<std::vector<double>> corners_plus_edge_midpoints(
    const std::vector<std::vector<double>>& rCorners, CellType type) {
    std::vector<std::vector<double>> pts = rCorners;
    for (const auto& e : meshioplusplus::detail::cell_refine_edges(type)) {
        const auto& a = rCorners[e[0]];
        const auto& b = rCorners[e[1]];
        std::vector<double> mid(a.size());
        for (std::size_t c = 0; c < a.size(); ++c)
            mid[c] = 0.5 * (a[c] + b[c]);
        pts.push_back(mid);
    }
    return pts;
}

void expect_point_near(const std::array<double, 3>& rGot, const std::vector<double>& rExpected) {
    for (std::size_t c = 0; c < rExpected.size(); ++c)
        EXPECT_NEAR(rGot[c], rExpected[c], 1e-9);
}

/**
 * @brief `corners_plus_edge_midpoints`'s hexahedron27-only extension: also
 * appends one point per `cell_refine_quad_faces(Hexahedron)` row (the mean of
 * that face's 4 corners) and one body-centre point (the mean of all 8
 * corners). The test builds and controls every one of the 27 point VALUES
 * directly -- this is not "how meshio++ would compute node 26 elsewhere",
 * merely a convenient, geometrically-distinguishable fixture.
 */
std::vector<std::vector<double>> hex27_fixture_points(
    const std::vector<std::vector<double>>& rCorners) {
    std::vector<std::vector<double>> pts =
        corners_plus_edge_midpoints(rCorners, CellType::Hexahedron);  // 0-19
    for (const auto& face : meshioplusplus::detail::cell_refine_quad_faces(CellType::Hexahedron)) {
        std::vector<double> mid(3, 0.0);
        for (std::uint8_t corner : face)
            for (std::size_t c = 0; c < 3; ++c)
                mid[c] += rCorners[corner][c] / 4.0;
        pts.push_back(mid);  // 20-25
    }
    std::vector<double> body(3, 0.0);
    for (const auto& corner : rCorners)
        for (std::size_t c = 0; c < 3; ++c)
            body[c] += corner[c] / 8.0;
    pts.push_back(body);  // 26
    return pts;
}

}  // namespace

// ---------------------------------------------------------------------------
// GidOrdering: the written-bytes *permutation* oracle. Each test builds a cell
// whose points sit at geometrically-known positions with IDENTITY
// connectivity, writes it, then re-parses the RAW file and checks where each
// referenced point actually is -- with no round trip involved, since gidpost
// has no reader.
//
// WHAT THIS DOES AND DOES NOT PIN. It pins that the writer emits slots in the
// order it was handed them: any permutation gid.cpp applied (or later starts
// applying) fails here. It does NOT pin that meshio++'s order IS GiD's order,
// and cannot: `pts` is built from meshio++'s own edge table, so for an
// identity mapping the assertion reduces to pts[slot] == pts[slot]. Only an
// external oracle -- a file GiD itself wrote, or GiD rendering our output --
// can settle the convention, and neither is available here. See gid.cpp's
// cell-type table for the evidence the mapping actually rests on, and
// doc/formats/gid.md for the one case where the available sources conflict.
// ---------------------------------------------------------------------------

TEST(GidOrdering, Hexahedron20MatchesGidGeometry) {
    const std::vector<std::vector<double>> corners = {
        {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1},
    };
    const auto pts = corners_plus_edge_midpoints(corners, CellType::Hexahedron);
    ASSERT_EQ(pts.size(), 20u);
    std::vector<std::int64_t> row(20);
    for (int i = 0; i < 20; ++i)
        row[i] = i;
    const Mesh m = mt::make_mesh(pts, "hexahedron20", {row});

    const std::string path = mt::temp_path(".post.msh");
    meshioplusplus::write_gid(path, m, meshioplusplus::GidMode::Ascii);

    const auto parsed = gid_parse_ascii_msh(path);
    ASSERT_EQ(parsed.mMeshes.size(), 1u);
    ASSERT_EQ(parsed.mMeshes[0].mElementRows.size(), 1u);
    const auto& row_out = parsed.mMeshes[0].mElementRows[0];
    ASSERT_EQ(row_out.size(), 21u);  // id + 20 nodes

    for (std::size_t slot = 0; slot < 20; ++slot) {
        const long long node_id = row_out[slot + 1];
        expect_point_near(parsed.mCoords.at(node_id), pts[slot]);
    }
}

TEST(GidOrdering, Tetra10MatchesGidGeometry) {
    const std::vector<std::vector<double>> corners = {
        {0, 0, 0},
        {1, 0, 0},
        {0, 1, 0},
        {0, 0, 1},
    };
    const auto pts = corners_plus_edge_midpoints(corners, CellType::Tetra);
    ASSERT_EQ(pts.size(), 10u);
    std::vector<std::int64_t> row(10);
    for (int i = 0; i < 10; ++i)
        row[i] = i;
    const Mesh m = mt::make_mesh(pts, "tetra10", {row});

    const std::string path = mt::temp_path(".post.msh");
    meshioplusplus::write_gid(path, m, meshioplusplus::GidMode::Ascii);

    const auto parsed = gid_parse_ascii_msh(path);
    ASSERT_EQ(parsed.mMeshes.size(), 1u);
    const auto& row_out = parsed.mMeshes[0].mElementRows.at(0);
    ASSERT_EQ(row_out.size(), 11u);
    for (std::size_t slot = 0; slot < 10; ++slot)
        expect_point_near(parsed.mCoords.at(row_out[slot + 1]), pts[slot]);
}

TEST(GidOrdering, Triangle6MatchesGidGeometry) {
    const std::vector<std::vector<double>> corners = {
        {0, 0, 0},
        {1, 0, 0},
        {0, 1, 0},
    };
    const auto pts = corners_plus_edge_midpoints(corners, CellType::Triangle);
    ASSERT_EQ(pts.size(), 6u);
    std::vector<std::int64_t> row = {0, 1, 2, 3, 4, 5};
    const Mesh m = mt::make_mesh(pts, "triangle6", {row});

    const std::string path = mt::temp_path(".post.msh");
    meshioplusplus::write_gid(path, m, meshioplusplus::GidMode::Ascii);

    const auto parsed = gid_parse_ascii_msh(path);
    const auto& row_out = parsed.mMeshes.at(0).mElementRows.at(0);
    ASSERT_EQ(row_out.size(), 7u);
    for (std::size_t slot = 0; slot < 6; ++slot)
        expect_point_near(parsed.mCoords.at(row_out[slot + 1]), pts[slot]);
}

TEST(GidOrdering, Quad8MatchesGidGeometry) {
    const std::vector<std::vector<double>> corners = {
        {0, 0, 0},
        {1, 0, 0},
        {1, 1, 0},
        {0, 1, 0},
    };
    const auto pts = corners_plus_edge_midpoints(corners, CellType::Quad);
    ASSERT_EQ(pts.size(), 8u);
    std::vector<std::int64_t> row = {0, 1, 2, 3, 4, 5, 6, 7};
    const Mesh m = mt::make_mesh(pts, "quad8", {row});

    const std::string path = mt::temp_path(".post.msh");
    meshioplusplus::write_gid(path, m, meshioplusplus::GidMode::Ascii);

    const auto parsed = gid_parse_ascii_msh(path);
    const auto& row_out = parsed.mMeshes.at(0).mElementRows.at(0);
    ASSERT_EQ(row_out.size(), 9u);
    for (std::size_t slot = 0; slot < 8; ++slot)
        expect_point_near(parsed.mCoords.at(row_out[slot + 1]), pts[slot]);
}

// The three orderings that were "not independently verified" until Kratos's
// own geometry classes (kratos/geometries/hexahedra_3d_27.h, prism_3d_15.h,
// pyramid_3d_13.h) were read directly and cross-checked against a second,
// Element-agnostic Kratos source (kratos/input_output/vtk_output.cpp's
// Kratos-to-VTK conversion, mirrored in ensight_output.cpp). See
// gid_common.hpp's `gid_cell_perm_table()` for the full derivation.
//
// The expected permutation is written out LITERALLY here rather than by
// calling gid_detail::gid_cell_perm() -- a bug in that table must not be able
// to validate itself. `perm[slot]` names which point in `pts` (built purely
// from meshio++'s own edge/face tables, in meshio++'s own node order) the
// GiD file's slot `slot` must hold.

TEST(GidOrdering, Hexahedron27MatchesGidGeometry) {
    const std::vector<std::vector<double>> corners = {
        {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1},
    };
    const auto pts = hex27_fixture_points(corners);
    ASSERT_EQ(pts.size(), 27u);
    std::vector<std::int64_t> row(27);
    for (int i = 0; i < 27; ++i)
        row[i] = i;
    const Mesh m = mt::make_mesh(pts, "hexahedron27", {row});

    const std::string path = mt::temp_path(".post.msh");
    meshioplusplus::write_gid(path, m, meshioplusplus::GidMode::Ascii);

    const auto parsed = gid_parse_ascii_msh(path);
    ASSERT_EQ(parsed.mMeshes.size(), 1u);
    const auto& row_out = parsed.mMeshes[0].mElementRows.at(0);
    ASSERT_EQ(row_out.size(), 28u);  // id + 27 nodes

    // Kratos's internal order (== what its GiD writer emits, since no
    // Conditions-only swap exists for this type): corners and bottom-ring
    // edges match meshio++'s own table; slots 12-15 <-> 16-19 (top ring vs.
    // verticals) and the face-centre block 20-25 are permuted; body centre
    // 26 is fixed. Independently confirmed against Kratos's own
    // vtk_output.cpp array (which expresses the identical mapping the other
    // way round).
    const std::array<std::size_t, 27> perm = {
        0,  1,  2,  3,           // corners 0-3
        4,  5,  6,  7,           // corners 4-7
        8,  9,  10, 11,          // bottom-ring edges (identity)
        16, 17, 18, 19,          // slot 12-15 (top ring)   <- meshio++'s verticals
        12, 13, 14, 15,          // slot 16-19 (verticals)  <- meshio++'s top ring
        24, 22, 21, 23, 20, 25,  // face centres 20-25
        26,                      // body centre
    };
    for (std::size_t slot = 0; slot < 27; ++slot)
        expect_point_near(parsed.mCoords.at(row_out[slot + 1]), pts[perm[slot]]);
}

TEST(GidOrdering, Wedge15MatchesGidGeometry) {
    const std::vector<std::vector<double>> corners = {
        {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {0, 1, 1},
    };
    const auto pts = corners_plus_edge_midpoints(corners, CellType::Wedge);
    ASSERT_EQ(pts.size(), 15u);
    std::vector<std::int64_t> row(15);
    for (int i = 0; i < 15; ++i)
        row[i] = i;
    const Mesh m = mt::make_mesh(pts, "wedge15", {row});

    const std::string path = mt::temp_path(".post.msh");
    meshioplusplus::write_gid(path, m, meshioplusplus::GidMode::Ascii);

    const auto parsed = gid_parse_ascii_msh(path);
    ASSERT_EQ(parsed.mMeshes.size(), 1u);
    const auto& row_out = parsed.mMeshes[0].mElementRows.at(0);
    ASSERT_EQ(row_out.size(), 16u);  // id + 15 nodes

    // Corners and bottom-triangle edges match meshio++'s own table; slots
    // 9-11 <-> 12-14 (top triangle vs. verticals) are swapped, the same
    // reverse-split hexahedron27 has.
    const std::array<std::size_t, 15> perm = {
        0, 1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 9, 10, 11,
    };
    for (std::size_t slot = 0; slot < 15; ++slot)
        expect_point_near(parsed.mCoords.at(row_out[slot + 1]), pts[perm[slot]]);
}

TEST(GidOrdering, Pyramid13MatchesGidGeometry) {
    const std::vector<std::vector<double>> corners = {
        {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0.5, 0.5, 1},
    };
    const auto pts = corners_plus_edge_midpoints(corners, CellType::Pyramid);
    ASSERT_EQ(pts.size(), 13u);
    std::vector<std::int64_t> row(13);
    for (int i = 0; i < 13; ++i)
        row[i] = i;
    const Mesh m = mt::make_mesh(pts, "pyramid13", {row});

    const std::string path = mt::temp_path(".post.msh");
    meshioplusplus::write_gid(path, m, meshioplusplus::GidMode::Ascii);

    const auto parsed = gid_parse_ascii_msh(path);
    ASSERT_EQ(parsed.mMeshes.size(), 1u);
    const auto& row_out = parsed.mMeshes[0].mElementRows.at(0);
    ASSERT_EQ(row_out.size(), 14u);  // id + 13 nodes

    // Kratos's Pyramid3D13 order is IDENTICAL to meshio++'s own (confirmed
    // via vtk_output.cpp explicitly skipping any conversion for this type) --
    // this is a regression lock, not a live risk the way the other two are.
    for (std::size_t slot = 0; slot < 13; ++slot)
        expect_point_near(parsed.mCoords.at(row_out[slot + 1]), pts[slot]);
}

// ---------------------------------------------------------------------------
// Structural tests: multi-block id uniqueness, the empty-coordinates-block
// requirement for meshes 2..n, Gauss-point-set referencing, and the
// unsupported-type/compiled-out error paths.
// ---------------------------------------------------------------------------

TEST(GidWrite, MultiBlockElementIdsAreGloballyUnique) {
    Mesh m;
    m.AssignPoints(mt::points_from({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0.5, 0.5, 1}}));
    m.AddCellBlock("triangle", mt::conn_from({{0, 1, 2}, {0, 2, 3}}));
    m.AddCellBlock("quad", mt::conn_from({{0, 1, 2, 3}}));

    const std::string path = mt::temp_path(".post.msh");
    meshioplusplus::write_gid(path, m, meshioplusplus::GidMode::Ascii);

    const auto parsed = gid_parse_ascii_msh(path);
    ASSERT_EQ(parsed.mMeshes.size(), 2u);
    std::vector<long long> ids;
    for (const auto& block : parsed.mMeshes)
        for (const auto& row : block.mElementRows)
            ids.push_back(row[0]);
    ASSERT_EQ(ids.size(), 3u);
    // Globally unique AND contiguous 1..n (the whole point of a shared
    // element-id space across blocks).
    std::sort(ids.begin(), ids.end());
    for (std::size_t i = 0; i < ids.size(); ++i)
        EXPECT_EQ(ids[i], static_cast<long long>(i + 1));
}

TEST(GidWrite, SecondMeshHasEmptyCoordinatesBlock) {
    Mesh m;
    m.AssignPoints(mt::points_from({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}}));
    m.AddCellBlock("triangle", mt::conn_from({{0, 1, 2}}));
    m.AddCellBlock("triangle", mt::conn_from({{0, 2, 3}}));

    const std::string path = mt::temp_path(".post.msh");
    meshioplusplus::write_gid(path, m, meshioplusplus::GidMode::Ascii);

    const auto parsed = gid_parse_ascii_msh(path);
    ASSERT_EQ(parsed.mMeshes.size(), 2u);
    EXPECT_TRUE(parsed.mMeshes[0].mHadCoordinatesRows);
    EXPECT_FALSE(parsed.mMeshes[1].mHadCoordinatesRows);
}

TEST(GidWrite, UnsupportedCellTypeThrowsByName) {
    Mesh m;
    m.AssignPoints(mt::points_from({{0, 0}, {1, 0}, {1, 1}, {0, 1}, {0.5, 0.5}}));
    m.AddPolygonBlock("polygon", {{0, 1, 2, 3, 4}});

    const std::string path = mt::temp_path(".post.msh");
    try {
        meshioplusplus::write_gid(path, m, meshioplusplus::GidMode::Ascii);
        FAIL() << "expected a WriteError";
    } catch (const meshioplusplus::WriteError& e) {
        EXPECT_NE(std::string(e.what()).find("polygon"), std::string::npos) << e.what();
    }
}

TEST(GidWrite, ResultsReferenceGaussPointSets) {
    Mesh m;
    m.AssignPoints(mt::points_from({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}}));
    m.AddCellBlock("quad", mt::conn_from({{0, 1, 2, 3}}));
    std::vector<NDArray> pressure;
    pressure.push_back(mt::data_array({42.0}));
    m.AddCellData("pressure", std::move(pressure));

    const std::string path = mt::temp_path(".post.msh");
    meshioplusplus::write_gid(path, m, meshioplusplus::GidMode::Ascii);

    const std::string res_path = path.substr(0, path.size() - 3) + "res";
    std::ifstream rf(res_path);
    std::string content((std::istreambuf_iterator<char>(rf)), std::istreambuf_iterator<char>());
    EXPECT_NE(content.find("GaussPoints \"gp_quad_0\""), std::string::npos);
    EXPECT_NE(content.find("\"pressure\""), std::string::npos);
    EXPECT_NE(content.find("OnGaussPoints \"gp_quad_0\""), std::string::npos);
}

// ---------------------------------------------------------------------------
// Result types (Matrix, Complex, ...). The declaration rides a field_data key
// (kGidResultTypePrefix) because the registry's (path, mesh) writers cannot
// carry a side-channel struct -- see gid.hpp. Component counts are gidpost's
// own; component ORDER is quoted verbatim from CIMNE's Customization Manual,
// and is pinned literally below because the two complex families use
// opposite conventions (ComplexVector interleaved, ComplexMatrix blocked).
// ---------------------------------------------------------------------------

namespace {

/// A 2-triangle mesh with one declared point_data array whose values encode
/// their own (row, component) position, so a permutation anywhere in the
/// write path is visible rather than merely plausible.
Mesh gid_typed_mesh(meshioplusplus::GidResultType type, std::size_t k) {
    Mesh m;
    m.AssignPoints(mt::points_from({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}}));
    m.AddCellBlock("triangle", mt::conn_from({{0, 1, 2}, {1, 3, 2}}));
    std::vector<double> vals(4 * k);
    for (std::size_t i = 0; i < vals.size(); ++i)
        vals[i] = static_cast<double>(i);
    m.AddPointData("f", mt::data_array(vals, k));
    NDArray decl(DType::Int64, {1});
    decl.As<std::int64_t>()[0] = static_cast<std::int64_t>(type);
    m.AddFieldData(std::string(meshioplusplus::kGidResultTypePrefix) + "f", std::move(decl));
    return m;
}

/// The numeric rows of the single Result block in a written .post.res.
std::vector<std::vector<double>> gid_read_values_block(const std::string& rResPath) {
    std::ifstream f(rResPath);
    std::string line;
    bool inside = false;
    std::vector<std::vector<double>> rows;
    while (std::getline(f, line)) {
        std::istringstream iss(line);
        std::string first;
        iss >> first;
        std::string lower = first;
        std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
        if (lower == "values") {
            inside = true;
            continue;
        }
        if (lower == "end" && inside)
            break;
        if (!inside)
            continue;
        std::vector<double> row;
        std::istringstream row_stream(line);
        double v;
        while (row_stream >> v)
            row.push_back(v);
        if (!row.empty())
            rows.push_back(row);
    }
    return rows;
}

}  // namespace

TEST(GidWrite, DeclaredResultTypeReachesTheFileHeader) {
    const std::string path = mt::temp_path(".post.msh");
    meshioplusplus::write_gid(path, gid_typed_mesh(meshioplusplus::GidResultType::Matrix, 6),
                              meshioplusplus::GidMode::Ascii);
    std::ifstream rf(path.substr(0, path.size() - 3) + "res");
    std::string content((std::istreambuf_iterator<char>(rf)), std::istreambuf_iterator<char>());
    EXPECT_NE(content.find("\"f\" \"meshio++\" 1 Matrix OnNodes"), std::string::npos) << content;
}

TEST(GidWrite, ComplexVectorIsInterleavedAndComplexMatrixIsBlocked) {
    // ComplexVector:6 = x_re x_im y_re y_im z_re z_im (interleaved);
    // ComplexMatrix:12 = 6 reals then 6 imaginaries (blocked). Values are
    // stored verbatim in GiD's order -- this pins that "verbatim" is a
    // checkable claim, not a hope.
    for (const auto& [type, k] : std::vector<std::pair<meshioplusplus::GidResultType, std::size_t>>{
             {meshioplusplus::GidResultType::ComplexVector, 6},
             {meshioplusplus::GidResultType::ComplexMatrix, 12}}) {
        const std::string path = mt::temp_path(".post.msh");
        meshioplusplus::write_gid(path, gid_typed_mesh(type, k), meshioplusplus::GidMode::Ascii);
        const auto rows = gid_read_values_block(path.substr(0, path.size() - 3) + "res");
        ASSERT_EQ(rows.size(), 4u);
        for (std::size_t r = 0; r < rows.size(); ++r) {
            ASSERT_EQ(rows[r].size(), k + 1) << "row " << r;
            EXPECT_DOUBLE_EQ(rows[r][0], static_cast<double>(r + 1));  // node id
            for (std::size_t c = 0; c < k; ++c)
                EXPECT_DOUBLE_EQ(rows[r][c + 1], static_cast<double>(r * k + c))
                    << "type=" << static_cast<int>(type) << " row=" << r << " comp=" << c;
        }
    }
}

TEST(GidWrite, MatrixNeedsNoPermutation) {
    // GiD's Matrix:6 is Sxx Syy Szz Sxy Syz Sxz -- already meshio/VTK's
    // symmetric-tensor order -- so values must reach the file untouched.
    const std::string path = mt::temp_path(".post.msh");
    meshioplusplus::write_gid(path, gid_typed_mesh(meshioplusplus::GidResultType::Matrix, 6),
                              meshioplusplus::GidMode::Ascii);
    const auto rows = gid_read_values_block(path.substr(0, path.size() - 3) + "res");
    ASSERT_FALSE(rows.empty());
    ASSERT_EQ(rows[0].size(), 7u);
    for (std::size_t c = 0; c < 6; ++c)
        EXPECT_DOUBLE_EQ(rows[0][c + 1], static_cast<double>(c));
}

TEST(GidWrite, IllegalComponentCountErrorsByName) {
    const std::string path = mt::temp_path(".post.msh");
    try {
        meshioplusplus::write_gid(path, gid_typed_mesh(meshioplusplus::GidResultType::Matrix, 5),
                                  meshioplusplus::GidMode::Ascii);
        FAIL() << "expected a WriteError";
    } catch (const meshioplusplus::WriteError& e) {
        const std::string msg = e.what();
        EXPECT_NE(msg.find("'f'"), std::string::npos) << msg;
        EXPECT_NE(msg.find("accepts 3, 6"), std::string::npos) << msg;
    }
}

TEST(GidWrite, OutOfRangeResultTypeErrorsByName) {
    Mesh m = gid_typed_mesh(meshioplusplus::GidResultType::Matrix, 6);
    NDArray bad(DType::Int64, {1});
    bad.As<std::int64_t>()[0] = 42;
    m.AddFieldData(std::string(meshioplusplus::kGidResultTypePrefix) + "f", std::move(bad));
    const std::string path = mt::temp_path(".post.msh");
    try {
        meshioplusplus::write_gid(path, m, meshioplusplus::GidMode::Ascii);
        FAIL() << "expected a WriteError";
    } catch (const meshioplusplus::WriteError& e) {
        const std::string msg = e.what();
        EXPECT_NE(msg.find("42"), std::string::npos) << msg;
        EXPECT_NE(msg.find("0..8"), std::string::npos) << msg;
    }
}

TEST(GidWrite, UndeclaredSixComponentArrayStillSplits) {
    // No declaration -> the historical inference is unaffected: (n,6) is
    // genuinely ambiguous (Matrix:6? ComplexMatrix:3? ComplexVector:6?), so it
    // must NOT silently become a Matrix just because the count matches.
    Mesh m;
    m.AssignPoints(mt::points_from({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}}));
    m.AddCellBlock("triangle", mt::conn_from({{0, 1, 2}, {1, 3, 2}}));
    m.AddPointData("s6", mt::data_array({0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11,
                                         12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23},
                                        6));
    const std::string path = mt::temp_path(".post.msh");
    meshioplusplus::write_gid(path, m, meshioplusplus::GidMode::Ascii);
    std::ifstream rf(path.substr(0, path.size() - 3) + "res");
    std::string content((std::istreambuf_iterator<char>(rf)), std::istreambuf_iterator<char>());
    EXPECT_EQ(content.find("Matrix"), std::string::npos) << content;
    for (int c = 1; c <= 6; ++c)
        EXPECT_NE(content.find("\"s6_" + std::to_string(c) + "\""), std::string::npos);
}

// ---------------------------------------------------------------------------
// Reader. Note the GidOrdering suite above is deliberately KEPT as a raw-bytes
// oracle now that a reader exists: a reader+writer round trip is a weak oracle
// that a consistently-wrong permutation survives, so it cannot replace it.
// ---------------------------------------------------------------------------

TEST(GidRead, AsciiRoundTrip) {
    const Mesh in = mt::tri_mesh();
    const std::string path = mt::temp_path(".post.msh");
    meshioplusplus::write_gid(path, in, meshioplusplus::GidMode::Ascii);

    const Mesh out = meshioplusplus::read_gid(path);
    ASSERT_EQ(out.NumPoints(), in.NumPoints());
    ASSERT_EQ(out.NumCellBlocks(), in.NumCellBlocks());
    EXPECT_EQ(out.Cells(0).Type(), in.Cells(0).Type());
    EXPECT_EQ(out.Cells(0).NumCells(), in.Cells(0).NumCells());
    for (std::size_t i = 0; i < in.NumPoints() * in.PointDim(); ++i)
        EXPECT_NEAR(meshioplusplus::detail::read_double(out.Points(), i),
                    meshioplusplus::detail::read_double(in.Points(), i), 1e-8);
}

TEST(GidRead, RepeatedNodeTablesAreDeduplicated) {
    // Real files (Kratos's own GiD output) repeat the FULL node table in every
    // MESH block; meshio++'s writer emits it once and writes empty Coordinates
    // pairs thereafter. This shape is unreachable through our own writer, so it
    // is spelled out by hand.
    const std::string path = mt::temp_path(".post.msh");
    {
        std::ofstream f(path);
        f << "MESH \"a\" dimension 3 ElemType Triangle Nnode 3\n"
             "Coordinates\n1 0 0 0\n2 1 0 0\n3 0 1 0\nEnd Coordinates\n"
             "Elements\n1 1 2 3\nEnd Elements\n"
             "MESH \"b\" dimension 3 ElemType Triangle Nnode 3\n"
             "Coordinates\n1 0 0 0\n2 1 0 0\n3 0 1 0\nEnd Coordinates\n"
             "Elements\n1 1 2 3\nEnd Elements\n";
    }
    const Mesh out = meshioplusplus::read_gid(path);
    EXPECT_EQ(out.NumPoints(), 3u);  // 3, not 6
    EXPECT_EQ(out.NumCellBlocks(), 2u);
}

TEST(GidRead, ElementIdsMayRestartPerBlock) {
    // Element ids are NOT globally unique in real files -- gidpost's own block
    // writer numbers 1..n per mesh. A reader keying results off one global map
    // would mis-associate them.
    const std::string path = mt::temp_path(".post.msh");
    {
        std::ofstream f(path);
        f << "MESH \"a\" dimension 3 ElemType Triangle Nnode 3\n"
             "Coordinates\n1 0 0 0\n2 1 0 0\n3 0 1 0\n4 1 1 0\nEnd Coordinates\n"
             "Elements\n1 1 2 3\nEnd Elements\n"
             "MESH \"b\" dimension 3 ElemType Triangle Nnode 3\n"
             "Coordinates\nEnd Coordinates\n"
             "Elements\n1 2 4 3\n2 1 2 4\nEnd Elements\n";
    }
    const Mesh out = meshioplusplus::read_gid(path);
    ASSERT_EQ(out.NumCellBlocks(), 2u);
    EXPECT_EQ(out.Cells(0).NumCells(), 1u);
    EXPECT_EQ(out.Cells(1).NumCells(), 2u);
}

TEST(GidRead, MaterialColumnIsDisambiguatedByNnode) {
    // There is NO separator between connectivity and the optional trailing
    // material id, so the row width against Nnode is the only disambiguator.
    const std::string with_mat = mt::temp_path(".post.msh");
    {
        std::ofstream f(with_mat);
        f << "MESH \"a\" dimension 3 ElemType Triangle Nnode 3\n"
             "Coordinates\n1 0 0 0\n2 1 0 0\n3 0 1 0\nEnd Coordinates\n"
             "Elements\n1 1 2 3 42\nEnd Elements\n";
    }
    const Mesh a = meshioplusplus::read_gid(with_mat);
    ASSERT_TRUE(a.HasCellData("gmsh:physical"));
    EXPECT_EQ(meshioplusplus::detail::read_int(a.CellData("gmsh:physical", 0), 0), 42);

    const std::string without = mt::temp_path(".post.msh");
    {
        std::ofstream f(without);
        f << "MESH \"a\" dimension 3 ElemType Triangle Nnode 3\n"
             "Coordinates\n1 0 0 0\n2 1 0 0\n3 0 1 0\nEnd Coordinates\n"
             "Elements\n1 1 2 3\nEnd Elements\n";
    }
    // No material information at all -- inventing an all-zero array would be
    // data the file never carried.
    EXPECT_FALSE(meshioplusplus::read_gid(without).HasCellData("gmsh:physical"));
}

TEST(GidRead, SuppressedRepeatedIdIsAContinuationRow) {
    // gidpost omits a Values row's id when it repeats the previous row's, so a
    // G>1 Gauss-point result writes the id once and the following rows begin
    // with whitespace. Parsing must survive that (the result is then dropped,
    // below, because meshio++ cannot represent per-point values).
    const std::string mesh_path = mt::temp_path(".post.msh");
    const std::string res_path = mesh_path.substr(0, mesh_path.size() - 3) + "res";
    {
        std::ofstream f(mesh_path);
        f << "MESH \"s\" dimension 3 ElemType Triangle Nnode 3\n"
             "Coordinates\n1 0 0 0\n2 1 0 0\n3 0 1 0\nEnd Coordinates\n"
             "Elements\n1 1 2 3\nEnd Elements\n";
    }
    {
        std::ofstream f(res_path);
        f << "GiD Post Results File 1.2\n"
             "GaussPoints \"g2\" ElemType Triangle \"s\"\n"
             "Number Of Gauss Points: 2\n"
             "Natural Coordinates: Internal\nEnd GaussPoints\n"
             "Result \"two\" \"a\" 1 Scalar OnGaussPoints \"g2\"\n"
             "Values\n1 5\n 6\nEnd Values\n";
    }
    const Mesh out = meshioplusplus::read_gid(mesh_path);
    EXPECT_EQ(out.NumPoints(), 3u);
    // Parsed without error, then dropped rather than averaged or truncated.
    EXPECT_FALSE(out.HasCellData("two"));
}

TEST(GidRead, ResultsSiblingIsOptional) {
    const std::string path = mt::temp_path(".post.msh");
    {
        std::ofstream f(path);
        f << "MESH \"s\" dimension 3 ElemType Triangle Nnode 3\n"
             "Coordinates\n1 0 0 0\n2 1 0 0\n3 0 1 0\nEnd Coordinates\n"
             "Elements\n1 1 2 3\nEnd Elements\n";
    }
    const Mesh out = meshioplusplus::read_gid(path);  // no .post.res exists
    EXPECT_EQ(out.NumPoints(), 3u);
    EXPECT_EQ(out.NumPointData(), 0u);
}

TEST(GidRead, ElementTypesWithNoMeshioCounterpartAreRefusedByName) {
    // hexahedron27/wedge15/pyramid13 used to be exactly this test's fixture --
    // now Kratos-derived orderings, they read successfully (see the
    // GidOrdering suite above), so this is repurposed to a type that is
    // genuinely and permanently unsupported: Sphere rows carry a radius, not
    // a node list, an entirely different row shape meshio++ has never mapped,
    // not merely an unverified ordering.
    const std::string path = mt::temp_path(".post.msh");
    {
        std::ofstream f(path);
        f << "MESH \"s\" dimension 3 ElemType Sphere Nnode 1\n"
             "Coordinates\n1 0 0 0\nEnd Coordinates\nElements\nEnd Elements\n";
    }
    try {
        meshioplusplus::read_gid(path);
        FAIL() << "expected a ReadError";
    } catch (const meshioplusplus::ReadError& e) {
        EXPECT_NE(std::string(e.what()).find("Sphere"), std::string::npos) << e.what();
    }
}

TEST(GidRead, MissingGeometryFileIsAnError) {
    EXPECT_THROW(meshioplusplus::read_gid(mt::temp_path(".post.msh")), meshioplusplus::ReadError);
}

// --- Grammar conformance against CIMNE's published specification. ----------
// None of these is reachable through our own writer (gidpost emits one fixed
// casing, always writes a mesh name, and always spells the 1-D type
// "Linear"), so a round trip cannot exercise them.

TEST(GidRead, KeywordsAreCaseInsensitive) {
    // GiD states this explicitly. The manual's own worked example opens a block
    // with `Coordinates` and closes it with `end coordinates`, so a
    // case-sensitive reader rejects the specification's own example file.
    const std::string path = mt::temp_path(".post.msh");
    {
        std::ofstream f(path);
        f << "mesh \"lower\" DIMENSION 3 elemtype Triangle nnode 3\n"
             "coordinates\n1 0 0 0\n2 1 0 0\n3 0 1 0\nEND COORDINATES\n"
             "Elements\n1 1 2 3\nend elements\n";
    }
    const Mesh out = meshioplusplus::read_gid(path);
    EXPECT_EQ(out.NumPoints(), 3u);
    ASSERT_EQ(out.NumCellBlocks(), 1u);
    EXPECT_EQ(out.Cells(0).Type(), "triangle");
}

TEST(GidRead, ANamelessMeshIsNotNamedDimension) {
    // The mesh name is optional. Parsing a nameless header alone cannot tell
    // the two readings apart -- the ElemType/Nnode scan finds its keywords
    // either way -- so this pins the one place the name is USED: binding a
    // GaussPoints set to a block. Taking the token after MESH unconditionally
    // names the mesh "dimension", and a set declared OnMesh "dimension" then
    // binds to it and attaches a result that does not belong to it.
    const std::string mesh_path = mt::temp_path(".post.msh");
    const std::string res_path = mesh_path.substr(0, mesh_path.size() - 9) + ".post.res";
    {
        std::ofstream f(mesh_path);
        f << "MESH dimension 3 ElemType Triangle Nnode 3\n"
             "Coordinates\n1 0 0 0\n2 1 0 0\n3 0 1 0\nEnd Coordinates\n"
             "Elements\n1 1 2 3\nEnd Elements\n";
    }
    {
        std::ofstream f(res_path);
        f << "GiD Post Results File 1.2\n"
             "GaussPoints \"gp\" ElemType Triangle \"dimension\"\n"
             "Number Of Gauss Points: 1\n"
             "Natural Coordinates: Internal\n"
             "End GaussPoints\n"
             "Result \"q\" \"a\" 1 Scalar OnGaussPoints \"gp\"\n"
             "Values\n1 5\nEnd Values\n";
    }
    const Mesh out = meshioplusplus::read_gid(mesh_path);
    EXPECT_FALSE(out.HasCellData("q"));
}

TEST(GidRead, BothSpellingsOfTheOneDimensionalTypeAreRead) {
    // gidpost emits "Linear"; CIMNE's current published grammar names the type
    // "Line". Both are CIMNE's, so files in the wild carry either.
    for (const char* spelling : {"Linear", "Line"}) {
        const std::string path = mt::temp_path(".post.msh");
        {
            std::ofstream f(path);
            f << "MESH \"l\" dimension 3 ElemType " << spelling
              << " Nnode 2\n"
                 "Coordinates\n1 0 0 0\n2 1 0 0\nEnd Coordinates\n"
                 "Elements\n1 1 2\nEnd Elements\n";
        }
        const Mesh out = meshioplusplus::read_gid(path);
        ASSERT_EQ(out.NumCellBlocks(), 1u) << spelling;
        EXPECT_EQ(out.Cells(0).Type(), "line") << spelling;
    }
}

TEST(GidRead, IsReadableRegardlessOfTheWriter) {
    // The reader needs no gidpost at all, so Ascii is readable in every build.
    EXPECT_TRUE(meshioplusplus::gid_readable(meshioplusplus::GidMode::Ascii));
}

#else  // !MESHIOPLUSPLUS_HAS_GIDPOST

TEST(GidRead, AsciiStaysReadableWithoutGidpost) {
    // The whole point of keeping the reader outside gid.cpp's gidpost guard:
    // a build that cannot WRITE GiD at all still reads the ascii flavour.
    EXPECT_TRUE(meshioplusplus::gid_readable(meshioplusplus::GidMode::Ascii));
}

TEST(GidWrite, CompiledOutThrowsNamingTheFlag) {
    meshioplusplus::Mesh m;
    try {
        meshioplusplus::write_gid("out.post.msh", m);
        FAIL() << "expected a WriteError";
    } catch (const meshioplusplus::WriteError& e) {
        const std::string msg = e.what();
        EXPECT_NE(msg.find("MESHIOPLUSPLUS_WITH_GIDPOST"), std::string::npos) << msg;
        EXPECT_NE(msg.find("MESHIOPLUSPLUS_WITH_ZLIB"), std::string::npos) << msg;
    }
}

#endif  // MESHIOPLUSPLUS_HAS_GIDPOST

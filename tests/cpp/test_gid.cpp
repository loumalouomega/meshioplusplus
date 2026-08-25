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
#include <string>
#include <vector>

using meshioplusplus::CellType;
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

}  // namespace

// ---------------------------------------------------------------------------
// GidOrdering: the written-bytes ordering oracle. Each test builds a cell
// whose points sit at geometrically-known positions with IDENTITY
// connectivity, writes it, then re-parses the RAW file and checks where each
// referenced point actually is -- independent of gid.cpp's own mapping
// table, and (since gidpost has no reader) with no round trip possible.
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
        {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1},
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
                    meshioplusplus::detail::read_double(in.Points(), i),
                    1e-8);
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

TEST(GidRead, UnverifiedOrderingIsRefusedByName) {
    // hexahedron27/wedge15/pyramid13 are exactly the orderings the WRITER
    // refuses; refusing them on read keeps the position consistent instead of
    // guessing a permutation in one direction we decline to guess in the other.
    const std::string path = mt::temp_path(".post.msh");
    {
        std::ofstream f(path);
        f << "MESH \"h\" dimension 3 ElemType Hexahedra Nnode 27\n"
             "Coordinates\n1 0 0 0\nEnd Coordinates\nElements\nEnd Elements\n";
    }
    try {
        meshioplusplus::read_gid(path);
        FAIL() << "expected a ReadError";
    } catch (const meshioplusplus::ReadError& e) {
        EXPECT_NE(std::string(e.what()).find("Hexahedra"), std::string::npos) << e.what();
    }
}

TEST(GidRead, MissingGeometryFileIsAnError) {
    EXPECT_THROW(meshioplusplus::read_gid(mt::temp_path(".post.msh")),
                 meshioplusplus::ReadError);
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

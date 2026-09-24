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
// Tests for the Elmer mesh directory reader and writer (formats/elmer).

// External includes
#include <gtest/gtest.h>

// System includes
#include <algorithm>
#include <cmath>
#include <cstring>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/elmer.hpp"
#include "meshioplusplus/operations/sniff.hpp"
#include "meshioplusplus/read_options.hpp"
#include "meshioplusplus/region.hpp"

using namespace meshioplusplus;
using namespace mt;
namespace fs = std::filesystem;

namespace {

fs::path fresh_dir() {
    const fs::path dir = temp_path("_elmer");
    fs::remove_all(dir);
    fs::create_directories(dir);
    return dir;
}

void put(const fs::path& rFile, const std::string& rText) {
    fs::create_directories(rFile.parent_path());
    std::ofstream(rFile, std::ios::binary) << rText;
}

std::string slurp(const fs::path& rFile) {
    std::ifstream in(rFile, std::ios::binary);
    std::stringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

std::vector<std::string> lines(const fs::path& rFile) {
    std::vector<std::string> out;
    std::istringstream in(slurp(rFile));
    std::string line;
    while (std::getline(in, line))
        out.push_back(line);
    return out;
}

const meshioplusplus::Region* find_region(const Mesh& rMesh, const std::string& rName) {
    for (std::size_t r = 0; r < rMesh.NumRegions(); ++r)
        if (rMesh.Region(r).mName == rName)
            return &rMesh.Region(r);
    return nullptr;
}

// Two unit tets glued on the face (1, 2, 3), in two bodies, with the outer
// face (0, 1, 3) and the interface as boundary triangles.
Mesh two_tets() {
    Mesh m;
    m.AssignPoints(points_from({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 1, 1}}));
    m.AddCellBlock("tetra", conn_from({{0, 1, 2, 3}, {1, 3, 2, 4}}));
    m.AddCellBlock("triangle", conn_from({{0, 1, 3}, {1, 2, 3}}));
    auto ids = [](std::vector<std::int64_t> v) {
        NDArray a(DType::Int64, {v.size()});
        std::copy(v.begin(), v.end(), a.As<std::int64_t>());
        return a;
    };
    m.AddRegion(meshioplusplus::Region("left", RegionKind::Cell, 3, 1, ids({0})));
    m.AddRegion(meshioplusplus::Region("right", RegionKind::Cell, 3, 2, ids({1})));
    m.AddRegion(meshioplusplus::Region("outer", RegionKind::Cell, 2, 7, ids({2})));
    m.AddRegion(meshioplusplus::Region("interface", RegionKind::Cell, 2, 8, ids({3})));
    return m;
}

}  // namespace

TEST(Elmer, WritesTheFilesElmerGridWrites) {
    const fs::path dir = fresh_dir() / "mesh";
    write_elmer(dir.string(), two_tets());
    EXPECT_EQ(lines(dir / "mesh.header"),
              (std::vector<std::string>{"5      2      2     ", "2     ", "303    2     ",
                                        "504    2     "}));
    EXPECT_EQ(lines(dir / "mesh.nodes")[4], "5 -1 1 1 1");
    EXPECT_EQ(lines(dir / "mesh.elements"),
              (std::vector<std::string>{"1 1 504 1 2 3 4", "2 2 504 2 4 3 5"}));
    // The outer face has one parent, the interface two.
    EXPECT_EQ(lines(dir / "mesh.boundary"),
              (std::vector<std::string>{"1 7 1 0 303 1 2 4", "2 8 1 2 303 2 3 4"}));
    const std::string names = slurp(dir / "mesh.names");
    EXPECT_NE(names.find("$ left = 1\n$ right = 2\n"), std::string::npos) << names;
    EXPECT_NE(names.find("$ outer = 7\n$ interface = 8\n"), std::string::npos) << names;
    EXPECT_EQ(sniff_format(dir.string()), "elmer");
    fs::remove_all(dir.parent_path());
}

TEST(Elmer, ProvenanceLivesInMeshNames) {
    const fs::path dir = fresh_dir();
    write_elmer(dir.string(), two_tets());
    // ElmerSolver parses only lines holding both `$` and `=`: the block holds neither `$`.
    for (const std::string& rLine : lines(dir / "mesh.names"))
        if (rLine.rfind("! ", 0) == 0)
            EXPECT_EQ(rLine.find('$'), std::string::npos) << rLine;
    for (const fs::path& rPath : {dir, dir / "mesh.header"}) {
        const auto found = detail::read_provenance_lines(rPath.string());
        EXPECT_TRUE(found.mRecognised) << rPath;
    }
    fs::remove_all(dir);
}

TEST(Elmer, RoundTripsBodiesAndBoundariesAsRegions) {
    const fs::path dir = fresh_dir();
    write_elmer(dir.string(), two_tets());
    const Mesh back = read_elmer((dir / "mesh.header").string());
    expect_points_close(two_tets(), back, 0.0);
    EXPECT_EQ(cell_rows(back), cell_rows(two_tets()));
    for (const auto& [name, tag, dim] : std::vector<std::tuple<std::string, int, int>>{
             {"left", 1, 3}, {"right", 2, 3}, {"outer", 7, 2}, {"interface", 8, 2}}) {
        const auto* reg = find_region(back, name);
        ASSERT_NE(reg, nullptr) << name;
        EXPECT_EQ(reg->mTag, tag);
        EXPECT_EQ(reg->mDim, dim);
        EXPECT_EQ(reg->mKind, RegionKind::Cell);
        EXPECT_EQ(reg->NumEntries(), 1u);
    }
    EXPECT_FALSE(back.HasCellData("partition:part"));
    fs::remove_all(dir);
}

TEST(Elmer, SideRegionsBecomeBoundaryElements) {
    Mesh m;
    m.AssignPoints(points_from({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}}));
    m.AddCellBlock("quad", conn_from({{0, 1, 2, 3}}));
    NDArray side(DType::Int64, {1, 2});
    side.As<std::int64_t>()[0] = 0;
    side.As<std::int64_t>()[1] = 1;  // quad edge 1 = (1, 2)
    m.AddRegion(meshioplusplus::Region("right", RegionKind::Side, 1, 4, std::move(side)));
    const fs::path dir = fresh_dir();
    write_elmer(dir.string(), m);
    EXPECT_EQ(lines(dir / "mesh.boundary"), (std::vector<std::string>{"1 4 1 0 202 2 3"}));
    const Mesh back = read_elmer(dir.string());
    const auto* reg = find_region(back, "right");
    ASSERT_NE(reg, nullptr);
    EXPECT_EQ(reg->mKind, RegionKind::Cell);
    EXPECT_EQ(reg->mDim, 1);
    fs::remove_all(dir);
}

TEST(Elmer, ReadsHexahedron20InElmerOrder) {
    // One 820 brick on [0,2]^3, its nodes listed in Elmer's order: corners,
    // bottom mid-edges, vertical mid-edges, top mid-edges.
    const double c[8][3] = {{0, 0, 0}, {2, 0, 0}, {2, 2, 0}, {0, 2, 0},
                            {0, 0, 2}, {2, 0, 2}, {2, 2, 2}, {0, 2, 2}};
    const int bottom[4][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
    std::ostringstream nodes;
    int id = 0;
    auto node = [&](double x, double y, double z) {
        nodes << ++id << " -1 " << x << ' ' << y << ' ' << z << '\n';
    };
    for (const auto& p : c)
        node(p[0], p[1], p[2]);
    for (const auto& e : bottom)
        node((c[e[0]][0] + c[e[1]][0]) / 2, (c[e[0]][1] + c[e[1]][1]) / 2, 0);
    for (int k = 0; k < 4; ++k)
        node(c[k][0], c[k][1], 1);
    for (const auto& e : bottom)
        node((c[e[0]][0] + c[e[1]][0]) / 2, (c[e[0]][1] + c[e[1]][1]) / 2, 2);
    const fs::path dir = fresh_dir();
    put(dir / "mesh.header", "20 1 0\n1\n820 1\n");
    put(dir / "mesh.nodes", nodes.str());
    std::string element = "1 1 820";
    for (int k = 1; k <= 20; ++k)
        element += " " + std::to_string(k);
    put(dir / "mesh.elements", element + "\n");
    put(dir / "mesh.boundary", "");
    const Mesh mesh = read_elmer(dir.string());
    ASSERT_EQ(mesh.NumCellBlocks(), 1u);
    EXPECT_EQ(mesh.Cells(0).Type(), "hexahedron20");
    // meshio++ (VTK) hexahedron20 edges: 8-11 bottom, 12-15 top, 16-19 vertical.
    const int edges[12][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6},
                              {6, 7}, {7, 4}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};
    const NDArray& conn = mesh.Cells(0).Conn();
    const NDArray& pts = mesh.Points();
    for (int e = 0; e < 12; ++e)
        for (int d = 0; d < 3; ++d) {
            const auto a = detail::read_int(conn, static_cast<std::size_t>(edges[e][0]));
            const auto b = detail::read_int(conn, static_cast<std::size_t>(edges[e][1]));
            const auto m = detail::read_int(conn, static_cast<std::size_t>(8 + e));
            EXPECT_DOUBLE_EQ(detail::read_double(pts, static_cast<std::size_t>(m * 3 + d)),
                             (detail::read_double(pts, static_cast<std::size_t>(a * 3 + d)) +
                              detail::read_double(pts, static_cast<std::size_t>(b * 3 + d))) /
                                 2)
                << "edge " << e;
        }
    // And back out in Elmer's order.
    const fs::path out = fresh_dir();
    write_elmer(out.string(), mesh);
    EXPECT_EQ(lines(out / "mesh.elements")[0], element);
    fs::remove_all(dir);
    fs::remove_all(out);
}

TEST(Elmer, NamesAndDefaults) {
    const fs::path dir = fresh_dir();
    put(dir / "mesh.header", "3 1 1\n2\n202 1\n303 1\n");
    put(dir / "mesh.nodes", "1 -1 0 0 0\n2 -1 1 0 0\n3 -1 0 1 0\n");
    put(dir / "mesh.elements", "1 3 303 1 2 3\n");
    put(dir / "mesh.boundary", "1 5 1 0 202 1 2\n");
    put(dir / "mesh.names",
        "! ----- names for bodies -----\n$ plate = 3\n! ----- names for boundaries -----\n"
        "$ plate = 5\n");
    Mesh mesh = read_elmer(dir.string());
    EXPECT_NE(find_region(mesh, "body:plate"), nullptr);
    EXPECT_NE(find_region(mesh, "boundary:plate"), nullptr);
    fs::remove(dir / "mesh.names");
    mesh = read_elmer(dir.string());
    EXPECT_EQ(find_region(mesh, "body_3")->mTag, 3);
    EXPECT_EQ(find_region(mesh, "boundary_5")->mTag, 5);
    fs::remove_all(dir);
}

TEST(Elmer, UnknownTypeCodesFailOrAreSkipped) {
    const fs::path dir = fresh_dir();
    put(dir / "mesh.header", "3 1 1\n2\n102 1\n303 1\n");
    put(dir / "mesh.nodes", "1 -1 0 0 0\n2 -1 1 0 0\n3 -1 0 1 0\n");
    put(dir / "mesh.elements", "1 1 303 1 2 3\n");
    put(dir / "mesh.boundary", "1 1 0 0 102 1 2\n");
    EXPECT_THROW(read_elmer(dir.string()), ReadError);
    ReadOptions lenient;
    lenient.mLenient = true;
    const Mesh mesh = read_elmer(dir.string(), lenient);
    EXPECT_EQ(mesh.NumCellBlocks(), 1u);
    EXPECT_THROW(read_elmer((dir / "nope").string()), ReadError);
    fs::remove_all(dir);
}

TEST(Elmer, MergesPartitionsAndLabelsCells) {
    // Two triangles split over two parts, sharing nodes 2 and 3; part 2 also
    // holds a halo copy of element 1, owned by part 1.
    const fs::path root = fresh_dir();
    const fs::path parts = root / "partitioning.2";
    put(parts / "part.1.header", "3 1 1\n2\n202 1\n303 1\n2 0\n");
    put(parts / "part.1.nodes", "1 -1 0 0 0\n2 -1 1 0 0\n3 -1 0 1 0\n");
    put(parts / "part.1.elements", "1 1 303 1 2 3\n");
    put(parts / "part.1.boundary", "1 4 1 0 202 1 2\n");
    put(parts / "part.2.header", "4 2 0\n1\n303 2\n2 0\n");
    put(parts / "part.2.nodes", "2 -1 1 0 0\n3 -1 0 1 0\n4 -1 1 1 0\n1 -1 0 0 0\n");
    put(parts / "part.2.elements", "2 1 303 2 4 3\n1/1 1 303 1 2 3\n");
    put(parts / "part.2.boundary", "");
    for (const fs::path& rPath : {parts, root}) {
        const Mesh mesh = read_elmer(rPath.string());
        EXPECT_EQ(mesh.NumPoints(), 4u);
        ASSERT_EQ(mesh.NumCellBlocks(), 2u);
        EXPECT_EQ(mesh.Cells(0).NumCells(), 2u);
        const NDArray& labels = mesh.CellData("partition:part", 0);
        EXPECT_EQ(detail::read_int(labels, 0), 0);
        EXPECT_EQ(detail::read_int(labels, 1), 1);
        EXPECT_EQ(detail::read_int(mesh.CellData("partition:part", 1), 0), 0);
    }
    ReadOptions one;
    one.mPieceSet = true;
    one.mPiece = -1;
    const Mesh second = read_elmer(parts.string(), one);
    EXPECT_EQ(second.Cells(0).NumCells(), 2u);  // its own element and the halo copy
    one.mPiece = 2;
    EXPECT_THROW(read_elmer(parts.string(), one), ReadError);
    EXPECT_EQ(sniff_format(parts.string()), "elmer");

    // A serial mesh next to its partitioning is read serially, labelled.
    put(root / "mesh.header", "4 2 1\n2\n202 1\n303 2\n");
    put(root / "mesh.nodes", "1 -1 0 0 0\n2 -1 1 0 0\n3 -1 0 1 0\n4 -1 1 1 0\n");
    put(root / "mesh.elements", "1 1 303 1 2 3\n2 1 303 2 4 3\n");
    put(root / "mesh.boundary", "1 4 2 0 202 2 4\n");
    const Mesh serial = read_elmer(root.string());
    EXPECT_EQ(detail::read_int(serial.CellData("partition:part", 0), 1), 1);
    EXPECT_EQ(detail::read_int(serial.CellData("partition:part", 1), 0), 1);
    fs::remove_all(root);
}

TEST(Elmer, WriterRejectsWhatElmerCannotHold) {
    Mesh points_only;
    points_only.AssignPoints(points_from({{0, 0, 0}}));
    points_only.AddCellBlock("vertex", conn_from({{0}}));
    const fs::path dir = fresh_dir();
    EXPECT_THROW(write_elmer(dir.string(), points_only), WriteError);
    Mesh polygon;
    polygon.AssignPoints(points_from({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 2, 0}}));
    polygon.AddCellBlock("polygon5", conn_from({{0, 1, 2, 3, 4}}));
    EXPECT_THROW(write_elmer(dir.string(), polygon), WriteError);
    fs::remove_all(dir);
}

TEST(Elmer, ReadsBinaryMeshesInEitherByteOrder) {
    // Write the text mesh, then its binary twin the way ElmerGrid -bin does
    // (fem/src/MeshIO.F90's stream layout), in both byte orders.
    const fs::path dir = fresh_dir();
    write_elmer((dir / "text").string(), two_tets());
    const Mesh text = read_elmer((dir / "text").string());
    for (bool big : {false, true}) {
        const fs::path bin = dir / (big ? "big" : "little");
        fs::create_directories(bin);
        fs::copy_file(dir / "text" / "mesh.header", bin / "mesh.header");
        fs::copy_file(dir / "text" / "mesh.names", bin / "mesh.names");
        auto put_bytes = [&](std::string& rOut, const void* pV, std::size_t N) {
            char b[8];
            std::memcpy(b, pV, N);
            if (big)
                std::reverse(b, b + N);
            rOut.append(b, N);
        };
        auto int32 = [&](std::string& rOut, std::int64_t V) {
            const std::int32_t v = static_cast<std::int32_t>(V);
            put_bytes(rOut, &v, 4);
        };
        std::string nodes, elements, boundary;
        for (const std::string& line : lines(dir / "text" / "mesh.nodes")) {
            std::istringstream in(line);
            std::int64_t id, part;
            double xyz[3];
            in >> id >> part >> xyz[0] >> xyz[1] >> xyz[2];
            int32(nodes, id);
            for (double v : xyz)
                put_bytes(nodes, &v, 8);
        }
        for (const auto& [file, out] : {std::make_pair("mesh.elements", &elements),
                                        std::make_pair("mesh.boundary", &boundary)}) {
            for (const std::string& line : lines(dir / "text" / file)) {
                std::istringstream in(line);
                std::vector<std::int64_t> v;
                for (std::int64_t x; in >> x;)
                    v.push_back(x);
                int32(*out, v[0]);
                int32(*out, -1);  // the owning part, -1 in a serial mesh
                for (std::size_t k = 1; k < v.size(); ++k)
                    int32(*out, v[k]);
            }
        }
        put(bin / "mesh.nodes.bin", nodes);
        put(bin / "mesh.elements.bin", elements);
        put(bin / "mesh.boundary.bin", boundary);
        const Mesh m = read_elmer(bin.string());
        ASSERT_EQ(m.NumCellBlocks(), text.NumCellBlocks());
        for (std::size_t k = 0; k < 15; ++k)
            EXPECT_EQ(detail::read_double(m.Points(), k), detail::read_double(text.Points(), k));
        EXPECT_EQ(m.NumRegions(), text.NumRegions());
    }
    fs::remove_all(dir);
}

TEST(Elmer, WritesPartitionsFromPartitionPart) {
    Mesh m = two_tets();
    auto labels = [](std::vector<std::int64_t> v) {
        NDArray a(DType::Int64, {v.size()});
        std::copy(v.begin(), v.end(), a.As<std::int64_t>());
        return a;
    };
    m.AddCellData("partition:part", {labels({0, 1}), labels({0, 0})});
    const fs::path dir = fresh_dir() / "mesh";
    write_elmer(dir.string(), m);
    const fs::path pdir = dir / "partitioning.2";
    EXPECT_EQ(lines(pdir / "part.1.elements"), (std::vector<std::string>{"1 1 504 1 2 3 4"}));
    EXPECT_EQ(lines(pdir / "part.2.elements"), (std::vector<std::string>{"2 2 504 2 4 3 5"}));
    // Nodes 2, 3, 4 are shared; part 1 owns them.
    EXPECT_EQ(lines(pdir / "part.2.shared"),
              (std::vector<std::string>{"2 2 1 2", "3 2 1 2", "4 2 1 2"}));
    // The interface keeps the parent of each side only.
    EXPECT_EQ(lines(pdir / "part.1.boundary"),
              (std::vector<std::string>{"1 7 1 0 303 1 2 4", "2 8 1 0 303 2 3 4"}));
    EXPECT_EQ(lines(pdir / "part.2.boundary"), (std::vector<std::string>{"2 8 0 2 303 2 3 4"}));
    EXPECT_EQ(lines(pdir / "part.2.header"),
              (std::vector<std::string>{"4      1      1     ", "2     ", "504    1     ",
                                        "303    1     ", "3      0     "}));
    // Read back merged: the same cells and the labels.
    const Mesh back = read_elmer(pdir.string());
    ASSERT_TRUE(back.HasCellData("partition:part"));
    EXPECT_EQ(detail::read_int(back.CellData("partition:part", 0), 1), 1);
    fs::remove_all(dir.parent_path());
}

TEST(Elmer, WritesTheHaloLayer) {
    // The two tetrahedra share a face, so each is the other part's halo (ElmerGrid
    // -halo): copied as `id/owner`, its nodes added to the part and shared with it,
    // the halo part listed after the owner.
    Mesh m = two_tets();
    NDArray tet_parts(DType::Int64, {2});
    tet_parts.As<std::int64_t>()[0] = 0;
    tet_parts.As<std::int64_t>()[1] = 1;
    NDArray tri_parts(DType::Int64, {2});
    tri_parts.As<std::int64_t>()[0] = 0;
    tri_parts.As<std::int64_t>()[1] = 0;
    m.AddCellData("partition:part", {std::move(tet_parts), std::move(tri_parts)});
    const fs::path dir = fresh_dir() / "mesh";
    write_elmer(dir.string(), m, true);
    const fs::path pdir = dir / "partitioning.2";
    EXPECT_EQ(lines(pdir / "part.1.elements"),
              (std::vector<std::string>{"1 1 504 1 2 3 4", "2/2 2 504 2 4 3 5"}));
    EXPECT_EQ(lines(pdir / "part.2.elements"),
              (std::vector<std::string>{"1/1 1 504 1 2 3 4", "2 2 504 2 4 3 5"}));
    EXPECT_EQ(lines(pdir / "part.1.shared"),
              (std::vector<std::string>{"1 2 1 2", "2 2 1 2", "3 2 1 2", "4 2 1 2", "5 2 2 1"}));
    EXPECT_EQ(lines(pdir / "part.1.header"),
              (std::vector<std::string>{"5      2      2     ", "2     ", "504    2     ",
                                        "303    2     ", "5      0     "}));
    // The merged read keeps each halo element once, with its owner's label.
    const Mesh back = read_elmer(pdir.string());
    ASSERT_EQ(back.Cells(0).NumCells(), 2u);
    EXPECT_EQ(detail::read_int(back.CellData("partition:part", 0), 1), 1);
    // Without the flag the same mesh writes no halo.
    const fs::path plain = dir.parent_path() / "plain";
    write_elmer(plain.string(), m, false);
    EXPECT_EQ(lines(plain / "partitioning.2" / "part.1.elements"),
              (std::vector<std::string>{"1 1 504 1 2 3 4"}));
    fs::remove_all(dir.parent_path());
}

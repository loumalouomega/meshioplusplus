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
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/elmer.hpp"
#include "meshioplusplus/read_options.hpp"
#include "meshioplusplus/region.hpp"

// Pins how the text Elmer reader reads hand-edited `mesh.nodes`,
// `mesh.elements`, `mesh.boundary` and `mesh.names`, so a change of its line
// walker (roadmap 3.1.1.1) can be checked against the reader it replaces. A
// mesh written by `write_elmer` is the reference; each formatting variant of
// its files must read to the same mesh, and the odd cases record what the
// reader did when this file was written, not what ElmerGrid asks for.

namespace fs = std::filesystem;

namespace {

using meshioplusplus::detail::read_double;
using meshioplusplus::detail::read_int;

fs::path elm_dir() {
    const fs::path dir = mt::temp_path("_elmtext");
    fs::remove_all(dir);
    fs::create_directories(dir);
    return dir;
}

std::string slurp(const fs::path& rPath) {
    std::ifstream f(rPath, std::ios::binary);
    return std::string((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
}

void spit(const fs::path& rPath, const std::string& rText) {
    fs::create_directories(rPath.parent_path());
    std::ofstream f(rPath, std::ios::binary);
    f << rText;
}

// Two tetrahedra in two bodies glued on a face, with two boundary triangles.
mt::Mesh two_tets() {
    mt::Mesh m;
    m.AssignPoints(mt::points_from({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 1, 1}}));
    m.AddCellBlock("tetra", mt::conn_from({{0, 1, 2, 3}, {1, 3, 2, 4}}));
    m.AddCellBlock("triangle", mt::conn_from({{0, 1, 3}, {1, 2, 3}}));
    auto ids = [](std::vector<std::int64_t> v) {
        meshioplusplus::NDArray a(meshioplusplus::DType::Int64, {v.size()});
        std::copy(v.begin(), v.end(), a.As<std::int64_t>());
        return a;
    };
    using meshioplusplus::RegionKind;
    m.AddRegion(meshioplusplus::Region("left", RegionKind::Cell, 3, 1, ids({0})));
    m.AddRegion(meshioplusplus::Region("right", RegionKind::Cell, 3, 2, ids({1})));
    m.AddRegion(meshioplusplus::Region("outer", RegionKind::Cell, 2, 7, ids({2})));
    m.AddRegion(meshioplusplus::Region("interface", RegionKind::Cell, 2, 8, ids({3})));
    return m;
}

// A reference mesh directory on disk.
class ElmCase {
public:
    ElmCase() : mBase(elm_dir()), mDir(mBase / "mesh") {
        meshioplusplus::write_elmer(mDir.string(), two_tets());
    }
    ~ElmCase() {
        std::error_code ec;
        fs::remove_all(mBase, ec);
    }
    const fs::path& Dir() const { return mDir; }
    std::string Get(const char* pName) const { return slurp(mDir / pName); }
    void Put(const char* pName, const std::string& rText) const { spit(mDir / pName, rText); }

private:
    fs::path mBase, mDir;
};

std::string describe(const mt::Mesh& rMesh) {
    std::string d = "P" + std::to_string(rMesh.NumPoints());
    for (std::size_t i = 0; i < rMesh.NumPoints() * rMesh.PointDim(); ++i)
        d += " " + std::to_string(read_double(rMesh.Points(), i));
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto cb = rMesh.Cells(b);
        d += " [" + std::string(cb.Type()) + " " + std::to_string(cb.NumCells()) + ":";
        for (std::size_t i = 0; i < cb.Conn().Size(); ++i)
            d += " " + std::to_string(read_int(cb.Conn(), i));
        d += "]";
    }
    for (std::size_t r = 0; r < rMesh.NumRegions(); ++r)
        d += " {" + rMesh.Region(r).mName + "}";
    return d;
}

std::string scrub(std::string Text, const fs::path& rDir) {
    const std::string dir = rDir.string();
    for (std::size_t at = 0; (at = Text.find(dir, at)) != std::string::npos; at += 5)
        Text.replace(at, dir.size(), "<dir>");
    return Text;
}

// What reading the directory gives: a description of the mesh or the error text.
std::string outcome(const fs::path& rDir, bool Lenient = false) {
    meshioplusplus::ReadOptions opts;
    opts.mLenient = Lenient;
    try {
        return describe(meshioplusplus::read_elmer(rDir.string(), opts));
    } catch (const meshioplusplus::ReadError& rExc) {
        return scrub(std::string("ReadError: ") + rExc.what(), rDir);
    } catch (const std::exception& rExc) {
        return scrub(std::string("exception: ") + rExc.what(), rDir);
    }
}

std::string restyle(const std::string& rText, const std::string& rLead, const std::string& rTrail,
                    const std::string& rEol, bool BlankAfter) {
    std::string out;
    std::size_t pos = 0;
    while (pos < rText.size()) {
        std::size_t eol = rText.find('\n', pos);
        if (eol == std::string::npos)
            eol = rText.size();
        out += rLead + rText.substr(pos, eol - pos) + rTrail + rEol;
        if (BlankAfter)
            out += rEol;
        pos = eol + 1;
    }
    return out;
}

void restyle_case(const ElmCase& rCase, const std::string& rLead, const std::string& rTrail,
                  const std::string& rEol, bool BlankAfter) {
    for (const char* name : {"mesh.nodes", "mesh.elements", "mesh.boundary", "mesh.names"})
        rCase.Put(name, restyle(rCase.Get(name), rLead, rTrail, rEol, BlankAfter));
}

void edit(const ElmCase& rCase, const char* pName, const std::string& rFrom,
          const std::string& rTo) {
    std::string t = rCase.Get(pName);
    const std::size_t at = t.find(rFrom);
    ASSERT_NE(at, std::string::npos) << pName << ": no '" << rFrom << "' in\n" << t;
    t.replace(at, rFrom.size(), rTo);
    rCase.Put(pName, t);
}

}  // namespace

TEST(ElmerText, TheReferenceMeshReadsWithItsRegions) {
    const ElmCase c;
    const std::string ref = outcome(c.Dir());
    EXPECT_EQ(ref.rfind("P5 ", 0), 0u) << ref;
    EXPECT_NE(ref.find("[tetra 2:"), std::string::npos) << ref;
    EXPECT_NE(ref.find("{left}"), std::string::npos) << ref;
    EXPECT_NE(ref.find("{interface}"), std::string::npos) << ref;
}

TEST(ElmerText, CrlfLineEndings) {
    const ElmCase c;
    const std::string ref = outcome(c.Dir());
    restyle_case(c, "", "", "\r\n", false);
    EXPECT_EQ(outcome(c.Dir()), ref);
}

TEST(ElmerText, TabsAndLeadingBlanks) {
    const ElmCase c;
    const std::string ref = outcome(c.Dir());
    restyle_case(c, " \t ", "\t ", "\n", false);
    EXPECT_EQ(outcome(c.Dir()), ref);
}

TEST(ElmerText, BlankLinesBetweenRecords) {
    const ElmCase c;
    const std::string ref = outcome(c.Dir());
    restyle_case(c, "", "", "\n", true);
    EXPECT_EQ(outcome(c.Dir()), ref);
}

TEST(ElmerText, AFinalLineWithNoNewline) {
    const ElmCase c;
    const std::string ref = outcome(c.Dir());
    for (const char* name : {"mesh.nodes", "mesh.elements", "mesh.boundary"}) {
        std::string t = c.Get(name);
        while (!t.empty() && t.back() == '\n')
            t.pop_back();
        c.Put(name, t);
    }
    EXPECT_EQ(outcome(c.Dir()), ref);
}

TEST(ElmerText, FormFeedsAndVerticalTabsAreNotSeparators) {
    // Only a space, a tab and a carriage return separate tokens here, so a form
    // feed glues two numbers into one bad token.
    const ElmCase c;
    edit(c, "mesh.nodes", "2 -1 1 0 0", "2 -1 1\f0 0");
    EXPECT_EQ(outcome(c.Dir()), "ReadError: Elmer mesh: <dir>/mesh.nodes:2: a node line needs `id part x y z`");
}

TEST(ElmerText, CoordinatesWithExponentsAndSigns) {
    const ElmCase c;
    const std::string ref = outcome(c.Dir());
    edit(c, "mesh.nodes", "2 -1 1 0 0", "2 -1 +1.0e0 -0 .0");
    const std::string got = outcome(c.Dir());
    EXPECT_EQ(got.rfind("P5 ", 0), 0u) << got;
}

TEST(ElmerText, ABadCoordinateIsNamed) {
    const ElmCase c;
    edit(c, "mesh.nodes", "2 -1 1 0 0", "2 -1 1x 0 0");
    EXPECT_EQ(outcome(c.Dir()), "ReadError: Elmer mesh: <dir>/mesh.nodes:2: bad coordinate '1x'");
}

TEST(ElmerText, NonFiniteAndHexadecimalCoordinates) {
    const ElmCase c;
    edit(c, "mesh.nodes", "2 -1 1 0 0", "2 -1 inf nan 0x1p3");
    EXPECT_EQ(outcome(c.Dir()), "P5 0.000000 0.000000 0.000000 inf nan 8.000000 0.000000 1.000000 0.000000 0.000000 0.000000 1.000000 1.000000 1.000000 1.000000 [tetra 2: 0 1 2 3 1 3 2 4] [triangle 2: 0 1 3 1 2 3] {interface} {left} {outer} {right}");
}

TEST(ElmerText, AnOverflowingCoordinate) {
    const ElmCase c;
    edit(c, "mesh.nodes", "2 -1 1 0 0", "2 -1 1e999 -1e999 1e-999");
    EXPECT_EQ(outcome(c.Dir()), "P5 0.000000 0.000000 0.000000 inf -inf 0.000000 0.000000 1.000000 0.000000 0.000000 0.000000 1.000000 1.000000 1.000000 1.000000 [tetra 2: 0 1 2 3 1 3 2 4] [triangle 2: 0 1 3 1 2 3] {interface} {left} {outer} {right}");
}

TEST(ElmerText, ANodeLineNeedsFiveTokens) {
    const ElmCase c;
    edit(c, "mesh.nodes", "2 -1 1 0 0", "2 -1 1 0");
    EXPECT_EQ(outcome(c.Dir()), "ReadError: Elmer mesh: <dir>/mesh.nodes:2: a node line needs `id part x y z`");
}

TEST(ElmerText, ANodeLineWithExtraTokensKeepsTheFirstFive) {
    const ElmCase c;
    const std::string ref = outcome(c.Dir());
    edit(c, "mesh.nodes", "2 -1 1 0 0", "2 -1 1 0 0 9 9");
    EXPECT_EQ(outcome(c.Dir()), ref);
}

TEST(ElmerText, ABadNodeId) {
    const ElmCase c;
    edit(c, "mesh.nodes", "2 -1 1 0 0", "two -1 1 0 0");
    EXPECT_EQ(outcome(c.Dir()), "ReadError: Elmer mesh: <dir>/mesh.nodes:2: bad node id 'two'");
}

TEST(ElmerText, ANodeIdWithTrailingText) {
    const ElmCase c;
    edit(c, "mesh.nodes", "2 -1 1 0 0", "2x -1 1 0 0");
    EXPECT_EQ(outcome(c.Dir()), "ReadError: Elmer mesh: <dir>/mesh.nodes:2: bad node id '2x'");
}

TEST(ElmerText, AnElementLineNeedsItsNodes) {
    const ElmCase c;
    edit(c, "mesh.elements", "1 1 504 1 2 3 4", "1 1 504 1 2 3");
    EXPECT_EQ(outcome(c.Dir()), "ReadError: Elmer mesh: <dir>/mesh.elements:1: type 504 needs 4 nodes");
}

TEST(ElmerText, AnElementLineWithTooManyNodes) {
    const ElmCase c;
    edit(c, "mesh.elements", "1 1 504 1 2 3 4", "1 1 504 1 2 3 4 5");
    EXPECT_EQ(outcome(c.Dir()), "ReadError: Elmer mesh: <dir>/mesh.elements:1: type 504 needs 4 nodes");
}

TEST(ElmerText, AnElementLineTooShort) {
    const ElmCase c;
    edit(c, "mesh.elements", "1 1 504 1 2 3 4", "1 1");
    EXPECT_EQ(outcome(c.Dir()), "ReadError: Elmer mesh: <dir>/mesh.elements:1: an element line needs `id body type nodes`");
}

TEST(ElmerText, AMalformedElementLine) {
    const ElmCase c;
    edit(c, "mesh.elements", "1 1 504 1 2 3 4", "1 x 504 1 2 3 4");
    EXPECT_EQ(outcome(c.Dir()), "ReadError: Elmer mesh: <dir>/mesh.elements:1: malformed line");
}

TEST(ElmerText, ABadNodeIdInAnElement) {
    const ElmCase c;
    edit(c, "mesh.elements", "1 1 504 1 2 3 4", "1 1 504 1 2 three 4");
    EXPECT_EQ(outcome(c.Dir()), "ReadError: Elmer mesh: <dir>/mesh.elements:1: bad node id 'three'");
}

TEST(ElmerText, AnElementIdWithAnOwningPart) {
    const ElmCase c;
    edit(c, "mesh.elements", "1 1 504", "1/1 1 504");
    const std::string got = outcome(c.Dir());
    EXPECT_EQ(got, "P5 0.000000 0.000000 0.000000 1.000000 0.000000 0.000000 0.000000 1.000000 0.000000 0.000000 0.000000 1.000000 1.000000 1.000000 1.000000 [tetra 2: 0 1 2 3 1 3 2 4] [triangle 2: 0 1 3 1 2 3] {interface} {left} {outer} {right}");
}

TEST(ElmerText, AnUnknownElementTypeFailsOrIsSkipped) {
    // Type 999 has no cell type; its line carries its 99 nodes.
    const ElmCase c;
    std::string line = "2 2 999";
    for (int k = 0; k < 99; ++k)
        line += " 1";
    edit(c, "mesh.elements", "2 2 504 2 4 3 5", line);
    EXPECT_EQ(outcome(c.Dir()),
              "ReadError: Elmer mesh: <dir>/mesh.elements:2: element type 999 has no meshio++ cell type");
    EXPECT_EQ(outcome(c.Dir(), true),
              "P5 0.000000 0.000000 0.000000 1.000000 0.000000 0.000000 0.000000 1.000000 0.000000 0.000000 0.000000 1.000000 1.000000 1.000000 1.000000 [tetra 1: 0 1 2 3] [triangle 2: 0 1 3 1 2 3] {interface} {left} {outer}");
}

TEST(ElmerText, ABoundaryLineNeedsItsFields) {
    const ElmCase c;
    edit(c, "mesh.boundary", "1 7 1 0 303 1 2 4", "1 7 1 0");
    EXPECT_EQ(outcome(c.Dir()), "ReadError: Elmer mesh: <dir>/mesh.boundary:1: a boundary line needs `id boundary parent1 parent2 type nodes`");
}

TEST(ElmerText, ABoundaryLineWithTheWrongNodeCount) {
    const ElmCase c;
    edit(c, "mesh.boundary", "1 7 1 0 303 1 2 4", "1 7 1 0 303 1 2");
    EXPECT_EQ(outcome(c.Dir()), "ReadError: Elmer mesh: <dir>/mesh.boundary:1: type 303 needs 3 nodes");
}

TEST(ElmerText, ABadParentIsZero) {
    const ElmCase c;
    const std::string ref = outcome(c.Dir());
    edit(c, "mesh.boundary", "1 7 1 0 303", "1 7 x 0 303");
    EXPECT_EQ(outcome(c.Dir()), ref);
}

TEST(ElmerText, AMissingNodesFile) {
    const ElmCase c;
    fs::remove(c.Dir() / "mesh.nodes");
    EXPECT_EQ(outcome(c.Dir()), "ReadError: Elmer mesh: missing <dir>/mesh.nodes");
}

TEST(ElmerText, AMissingElementsFile) {
    const ElmCase c;
    fs::remove(c.Dir() / "mesh.elements");
    EXPECT_EQ(outcome(c.Dir()), "ReadError: Elmer mesh: missing <dir>/mesh.elements");
}

TEST(ElmerText, AMissingBoundaryFileIsAMeshWithoutBoundaries) {
    const ElmCase c;
    fs::remove(c.Dir() / "mesh.boundary");
    EXPECT_EQ(outcome(c.Dir()), "ReadError: Elmer mesh: missing <dir>/mesh.boundary");
}

TEST(ElmerText, ACommentLineBeforeTheNamesIsIgnored) {
    const ElmCase c;
    const std::string ref = outcome(c.Dir());
    c.Put("mesh.names", "! a comment line\n" + c.Get("mesh.names"));
    EXPECT_EQ(outcome(c.Dir()), ref);
}

TEST(ElmerText, NamesWithTabsAndUppercaseHeaders) {
    const ElmCase c;
    c.Put("mesh.names",
          "! Names for bodies:\n$\tleft\t=\t1\n$right=2\n"
          "! NAMES FOR BOUNDARIES:\n$ outer = 7 \r\n$   interface   =   8\n");
    EXPECT_EQ(outcome(c.Dir()), "P5 0.000000 0.000000 0.000000 1.000000 0.000000 0.000000 0.000000 1.000000 0.000000 0.000000 0.000000 1.000000 1.000000 1.000000 1.000000 [tetra 2: 0 1 2 3 1 3 2 4] [triangle 2: 0 1 3 1 2 3] {interface} {left} {outer} {right}");
}

TEST(ElmerText, ANamesLineWithoutAnIdOrNameIsSkipped) {
    const ElmCase c;
    c.Put("mesh.names",
          "! names for bodies\n$ left = x\n$ = 2\n$ right = 2 trailing\n$ nodollar = 3\n"
          "! names for boundaries\n$ outer = 7\n");
    EXPECT_EQ(outcome(c.Dir()), "P5 0.000000 0.000000 0.000000 1.000000 0.000000 0.000000 0.000000 1.000000 0.000000 0.000000 0.000000 1.000000 1.000000 1.000000 1.000000 [tetra 2: 0 1 2 3 1 3 2 4] [triangle 2: 0 1 3 1 2 3] {body_1} {boundary_8} {outer} {right}");
}

TEST(ElmerText, ManyRecordsKeepTheirOrder) {
    const fs::path base = elm_dir();
    const fs::path dir = base / "big";
    constexpr std::size_t n = 3000;
    std::string nodes, elements;
    for (std::size_t i = 1; i <= n; ++i)
        nodes += std::to_string(i) + " -1 " + std::to_string(i) + " " + std::to_string(2 * i) +
                 " 0.5\n";
    for (std::size_t i = 1; i + 2 <= n; ++i)
        elements += std::to_string(i) + " 1 303 " + std::to_string(i) + " " + std::to_string(i + 1) +
                    " " + std::to_string(i + 2) + "\n";
    spit(dir / "mesh.header",
         std::to_string(n) + " " + std::to_string(n - 2) + " 0\n1\n303 " + std::to_string(n - 2) +
             "\n");
    spit(dir / "mesh.nodes", nodes);
    spit(dir / "mesh.elements", elements);
    spit(dir / "mesh.boundary", "");
    const std::string out = outcome(dir);
    EXPECT_EQ(out.rfind("P3000 ", 0), 0u) << out.substr(0, 80);
    EXPECT_NE(out.find("[triangle 2998:"), std::string::npos);
    fs::remove_all(base);
}

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
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/openfoam.hpp"

// Pins how the ASCII OpenFOAM reader reads hand-edited text, so a change of its
// tokenizer (roadmap 3.1.1.1) can be checked against the reader it replaces.
// A case written by `write_openfoam` is the reference; each formatting variant
// of its files (comments, CRLF, tabs, blank lines, header layouts) must read to
// the same mesh. The odd cases record what the reader did when this file was
// written -- lines it skips, counts it ignores -- not what the format asks.

namespace fs = std::filesystem;

namespace {

using meshioplusplus::detail::read_double;
using meshioplusplus::detail::read_int;

fs::path temp_case() {
    static std::atomic<unsigned> counter{0};
    return fs::temp_directory_path() / ("meshio_ofpin_" + std::to_string(counter++));
}

// Two hexahedra side by side, so `neighbour` holds one internal face.
mt::Mesh two_hexes() {
    mt::Mesh m;
    std::vector<std::vector<double>> pts;
    for (int k = 0; k < 2; ++k)
        for (int j = 0; j < 2; ++j)
            for (int i = 0; i < 3; ++i)
                pts.push_back({double(i), double(j), double(k)});
    auto id = [](int i, int j, int k) { return std::int64_t(i + 3 * j + 6 * k); };
    std::vector<std::vector<std::int64_t>> rows;
    for (int i = 0; i < 2; ++i)
        rows.push_back({id(i, 0, 0), id(i + 1, 0, 0), id(i + 1, 1, 0), id(i, 1, 0), id(i, 0, 1),
                        id(i + 1, 0, 1), id(i + 1, 1, 1), id(i, 1, 1)});
    m.AssignPoints(mt::points_from(pts));
    m.AddCellBlock("hexahedron", mt::conn_from(rows));
    return m;
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

// A reference case on disk; `Poly()` is its polyMesh directory.
class FoamCase {
public:
    FoamCase() : mBase(temp_case()) {
        meshioplusplus::write_openfoam((mBase / "case.foam").string(), two_hexes(),
                                       meshioplusplus::OpenFoamInfo{});
    }
    ~FoamCase() {
        std::error_code ec;
        fs::remove_all(mBase, ec);
    }
    fs::path Poly() const { return mBase / "constant" / "polyMesh"; }
    fs::path File() const { return mBase / "case.foam"; }
    const fs::path& Base() const { return mBase; }
    std::string Get(const char* pName) const { return slurp(Poly() / pName); }
    void Put(const char* pName, const std::string& rText) const { spit(Poly() / pName, rText); }

private:
    fs::path mBase;
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
    return d;
}

// The case directory's temporary name made stable, for messages that carry a path.
std::string scrub(std::string Text, const fs::path& rBase) {
    const std::string base = rBase.string();
    for (std::size_t at = 0; (at = Text.find(base, at)) != std::string::npos; at += 6)
        Text.replace(at, base.size(), "<case>");
    return Text;
}

// What reading the case gives: a description of the mesh or the error text.
std::string outcome(const fs::path& rFile) {
    meshioplusplus::OpenFoamInfo info;
    const fs::path base = rFile.parent_path();
    try {
        return describe(meshioplusplus::read_openfoam(rFile.string(), info));
    } catch (const meshioplusplus::ReadError& rExc) {
        return scrub(std::string("ReadError: ") + rExc.what(), base);
    } catch (const std::exception& rExc) {
        return scrub(std::string("exception: ") + rExc.what(), base);
    }
}

// Every line of @p rText re-laid-out: a prefix, a suffix, an end of line and a
// blank line after it.
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

void restyle_case(const FoamCase& rCase, const std::string& rLead, const std::string& rTrail,
                  const std::string& rEol, bool BlankAfter) {
    for (const char* name : {"points", "faces", "owner", "neighbour", "boundary"})
        rCase.Put(name, restyle(rCase.Get(name), rLead, rTrail, rEol, BlankAfter));
}

// Replace the first occurrence of @p rFrom in a case file.
void edit(const FoamCase& rCase, const char* pName, const std::string& rFrom,
          const std::string& rTo) {
    std::string t = rCase.Get(pName);
    const std::size_t at = t.find(rFrom);
    ASSERT_NE(at, std::string::npos) << pName << ": no '" << rFrom << "' in\n" << t;
    t.replace(at, rFrom.size(), rTo);
    rCase.Put(pName, t);
}

}  // namespace

TEST(OpenFoamText, TheReferenceCaseReadsAsTwoHexahedra) {
    const FoamCase c;
    const std::string ref = outcome(c.File());
    EXPECT_EQ(ref.rfind("P12 ", 0), 0u) << ref;
    EXPECT_NE(ref.find("[hexahedron 2:"), std::string::npos) << ref;
}

TEST(OpenFoamText, CrlfLineEndings) {
    const FoamCase c;
    const std::string ref = outcome(c.File());
    restyle_case(c, "", "", "\r\n", false);
    EXPECT_EQ(outcome(c.File()), ref);
}

TEST(OpenFoamText, TabsAndLeadingBlanks) {
    const FoamCase c;
    const std::string ref = outcome(c.File());
    restyle_case(c, "\t  ", "  \t", "\n", false);
    EXPECT_EQ(outcome(c.File()), ref);
}

TEST(OpenFoamText, BlankLinesBetweenEntries) {
    const FoamCase c;
    const std::string ref = outcome(c.File());
    restyle_case(c, "", "", "\n", true);
    EXPECT_EQ(outcome(c.File()), ref);
}

TEST(OpenFoamText, TrailingLineComments) {
    const FoamCase c;
    const std::string ref = outcome(c.File());
    restyle_case(c, "", " // a note ( with ) parentheses 12", "\n", false);
    EXPECT_EQ(outcome(c.File()), ref);
}

TEST(OpenFoamText, TrailingBlockComments) {
    const FoamCase c;
    const std::string ref = outcome(c.File());
    restyle_case(c, "", " /* a note ( with ) { braces } 7 */", "\n", false);
    EXPECT_EQ(outcome(c.File()), ref);
}

TEST(OpenFoamText, ABlockCommentSpanningLines) {
    const FoamCase c;
    const std::string ref = outcome(c.File());
    c.Put("points", "/* leading\n   banner\n*/\n" + c.Get("points") + "\n/* trailing\nlines */\n");
    c.Put("owner", "// only a line comment\n" + c.Get("owner"));
    EXPECT_EQ(outcome(c.File()), ref);
}

TEST(OpenFoamText, ACommentInsideAListBetweenItsEntries) {
    const FoamCase c;
    const std::string ref = outcome(c.File());
    edit(c, "owner", "(\n", "(\n// first\n/* second */\n");
    EXPECT_EQ(outcome(c.File()), ref);
}

TEST(OpenFoamText, TheFoamFileBlockWithItsBraceOnTheNextLine) {
    const FoamCase c;
    const std::string ref = outcome(c.File());
    for (const char* name : {"points", "faces", "owner", "neighbour", "boundary"}) {
        std::string t = c.Get(name);
        const std::size_t at = t.find("FoamFile");
        ASSERT_NE(at, std::string::npos) << name;
        c.Put(name, t.substr(0, at) + "FoamFile\n" + t.substr(at + 8));
    }
    EXPECT_EQ(outcome(c.File()), ref);
}

TEST(OpenFoamText, NoFoamFileHeaderAtAll) {
    const FoamCase c;
    const std::string ref = outcome(c.File());
    for (const char* name : {"points", "faces", "owner", "neighbour", "boundary"}) {
        std::string t = c.Get(name);
        const std::size_t open = t.find('{');
        int depth = 0;
        std::size_t close = std::string::npos;
        for (std::size_t p = open; p < t.size(); ++p) {
            if (t[p] == '{')
                ++depth;
            else if (t[p] == '}' && --depth == 0) {
                close = p;
                break;
            }
        }
        ASSERT_NE(close, std::string::npos) << name;
        c.Put(name, t.substr(close + 1));
    }
    EXPECT_EQ(outcome(c.File()), ref);
}

TEST(OpenFoamText, AWrongCountOnlySizesTheReservation) {
    const FoamCase c;
    const std::string ref = outcome(c.File());
    edit(c, "points", "\n12\n", "\n99999999\n");
    EXPECT_EQ(outcome(c.File()), ref);
    edit(c, "owner", "\n11\n", "\n3\n");
    EXPECT_EQ(outcome(c.File()), ref);
}

TEST(OpenFoamText, OwnerEntriesMaySharePhysicalLines) {
    const FoamCase c;
    const std::string ref = outcome(c.File());
    std::string t = c.Get("owner");
    const std::size_t open = t.find("\n(\n");
    ASSERT_NE(open, std::string::npos);
    std::string body = t.substr(open + 3);
    for (char& ch : body)
        if (ch == '\n')
            ch = ' ';
    // The closing parenthesis keeps a line of its own.
    const std::size_t close = body.rfind(')');
    ASSERT_NE(close, std::string::npos);
    c.Put("owner", t.substr(0, open + 3) + body.substr(0, close) + "\n)\n");
    EXPECT_EQ(outcome(c.File()), ref);
}

TEST(OpenFoamText, APointsListOnOneLineReadsNoPoints) {
    const FoamCase c;
    std::string t = c.Get("points");
    const std::size_t open = t.find("\n12\n");
    ASSERT_NE(open, std::string::npos);
    std::string flat = t.substr(open + 1);
    for (char& ch : flat)
        if (ch == '\n')
            ch = ' ';
    c.Put("points", t.substr(0, open + 1) + flat);
    EXPECT_EQ(outcome(c.File()), "ReadError: OpenFOAM: face 0 names point 1, but the mesh has 0 points");
}

TEST(OpenFoamText, AMissingPointsFileIsAnOpenError) {
    const FoamCase c;
    fs::remove(c.Poly() / "points");
    EXPECT_EQ(outcome(c.File()), "ReadError: Could not open OpenFOAM file: <case>/constant/polyMesh/points");
}

TEST(OpenFoamText, ABlockCommentOpenedBeforeTheListEndsAtTheBanner) {
    // Every file `write_openfoam` writes starts with a banner comment, so a
    // comment opened ahead of it closes at the banner's own `*/`.
    const FoamCase c;
    const std::string ref = outcome(c.File());
    c.Put("owner", "/* opens before the banner\n" + c.Get("owner"));
    EXPECT_EQ(outcome(c.File()), ref);
}

TEST(OpenFoamText, AnUnterminatedBlockCommentSwallowsTheRest) {
    // The closing parenthesis is swallowed with it; the list still reads.
    const FoamCase c;
    const std::string ref = outcome(c.File());
    edit(c, "owner", "\n)\n", "\n/* never closed\n)\n");
    EXPECT_EQ(outcome(c.File()), ref);
}

TEST(OpenFoamText, ABigEndianArchIsRefusedByName) {
    const FoamCase c;
    edit(c, "points", "format", "arch \"BSB;label=32;scalar=64\";\n format");
    EXPECT_EQ(outcome(c.File()), "ReadError: OpenFOAM: big-endian ('BSB') binary files are not supported, only little-endian ('LSB')");
}

TEST(OpenFoamText, AFacesLineWithoutParenthesesIsSkipped) {
    const FoamCase c;
    edit(c, "faces", "4(", "4 ");
    EXPECT_EQ(outcome(c.File()), "ReadError: OpenFOAM: owner lists 11 faces, but the mesh has 10");
}

TEST(OpenFoamText, APointsLineWithTwoNumbersIsSkipped) {
    const FoamCase c;
    edit(c, "points", "(0 0 0)", "(0 0)");
    EXPECT_EQ(outcome(c.File()), "ReadError: OpenFOAM: face 6 names point 11, but the mesh has 11 points");
}

TEST(OpenFoamText, APointsLineWithAWordIsSkipped) {
    const FoamCase c;
    edit(c, "points", "(0 0 0)", "(0 zero 0)");
    EXPECT_EQ(outcome(c.File()), "ReadError: OpenFOAM: face 6 names point 11, but the mesh has 11 points");
}

TEST(OpenFoamText, APointsLineWithFourNumbersKeepsTheFirstThree) {
    const FoamCase c;
    const std::string ref = outcome(c.File());
    edit(c, "points", "(0 0 0)", "(0 0 0 9)");
    EXPECT_EQ(outcome(c.File()), ref);
}

TEST(OpenFoamText, PointCoordinatesWithExponentsAndSigns) {
    const FoamCase c;
    const std::string ref = outcome(c.File());
    edit(c, "points", "(1 0 0)", "(+1.0e0 0 0.)");
    EXPECT_EQ(outcome(c.File()), ref);
}

TEST(OpenFoamText, ANonFiniteCoordinateSkipsThePoint) {
    const FoamCase c;
    edit(c, "points", "(2 1 1)", "(inf 1 1)");
    EXPECT_EQ(outcome(c.File()), "ReadError: OpenFOAM: face 6 names point 11, but the mesh has 11 points");
}

TEST(OpenFoamText, AFaceCountMayBeSeparatedFromItsParenthesis) {
    const FoamCase c;
    edit(c, "faces", "4(", "4  (");
    EXPECT_EQ(outcome(c.File()), outcome(FoamCase().File()));
}

TEST(OpenFoamText, TheBoundaryWithCommentsAndNestedBraces) {
    const FoamCase c;
    const std::string ref = outcome(c.File());
    std::string t = c.Get("boundary");
    const std::size_t at = t.find("type");
    ASSERT_NE(at, std::string::npos);
    t.insert(at, "// the first patch\n        inGroups List<word> 1(wall);\n        sub { a { b 1; } }\n        ");
    c.Put("boundary", t);
    EXPECT_EQ(outcome(c.File()), ref);
}

// ---- field files --------------------------------------------------------

namespace {

// A `0/<name>` field file; the case is read with its time step 0.
std::string field_outcome(const FoamCase& rCase, const std::string& rName,
                          const std::string& rBody, const std::string& rClass) {
    spit(rCase.Base() / "0" / rName,
         "FoamFile\n{\n format ascii;\n class " + rClass + ";\n object " + rName + ";\n}\n" +
             rBody);
    meshioplusplus::OpenFoamInfo info;
    try {
        const mt::Mesh m =
            meshioplusplus::read_openfoam(rCase.File().string(), meshioplusplus::ReadOptions{}, info);
        std::string d;
        for (std::size_t b = 0; b < m.NumCellBlocks(); ++b) {
            if (!m.HasCellData(rName)) {
                if (!m.HasPointData(rName))
                    return "none";
                break;
            }
            const auto& a = m.CellData(rName, b);
            d += "[" + std::string(m.Cells(b).Type()) + ":";
            for (std::size_t i = 0; i < a.Size(); ++i) {
                const double v = read_double(a, i);
                d += " " + (std::isnan(v) ? std::string("nan") : std::to_string(v));
            }
            d += "]";
        }
        return d;
    } catch (const meshioplusplus::ReadError& rExc) {
        return scrub(std::string("ReadError: ") + rExc.what(), rCase.Base());
    } catch (const std::exception& rExc) {
        return scrub(std::string("exception: ") + rExc.what(), rCase.Base());
    }
}

}  // namespace

TEST(OpenFoamTextFields, UniformScalarAndVector) {
    const FoamCase c;
    EXPECT_EQ(field_outcome(c, "p", "internalField   uniform 42;\n", "volScalarField"), "[hexahedron: 42.000000 42.000000][quad: nan nan nan nan nan nan nan nan nan nan]");
    EXPECT_EQ(field_outcome(c, "U", "internalField   uniform (1 2 3);\n", "volVectorField"), "[hexahedron: 1.000000 2.000000 3.000000 1.000000 2.000000 3.000000][quad: nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan]");
}

TEST(OpenFoamTextFields, UniformScalarWithTrailingTokens) {
    const FoamCase c;
    EXPECT_EQ(field_outcome(c, "p", "internalField   uniform 4.5e1 junk;\n", "volScalarField"),
              "[hexahedron: 45.000000 45.000000][quad: nan nan nan nan nan nan nan nan nan nan]");
    EXPECT_EQ(field_outcome(c, "q", "internalField   uniform junk;\n", "volScalarField"), "[hexahedron: 0.000000 0.000000][quad: nan nan nan nan nan nan nan nan nan nan]");
}

TEST(OpenFoamTextFields, NonuniformScalarWithBlankLinesAndComments) {
    const FoamCase c;
    EXPECT_EQ(field_outcome(c, "p",
                            "internalField   nonuniform List<scalar>\n\n2\n\n(\n1.5\n\n 2.5 \n)\n;\n",
                            "volScalarField"),
              "[hexahedron: 1.500000 2.500000][quad: nan nan nan nan nan nan nan nan nan nan]");
}

TEST(OpenFoamTextFields, NonuniformVectorOnePerLine) {
    const FoamCase c;
    EXPECT_EQ(field_outcome(c, "U",
                            "internalField   nonuniform List<vector>\n2\n(\n(1 2 3)\n  (4 5 6)  \n)\n;\n",
                            "volVectorField"),
              "[hexahedron: 1.000000 2.000000 3.000000 4.000000 5.000000 6.000000][quad: nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan nan]");
}

TEST(OpenFoamTextFields, NonuniformScalarLinesWithTrailingText) {
    const FoamCase c;
    EXPECT_EQ(field_outcome(c, "p",
                            "internalField   nonuniform List<scalar>\n2\n(\n1.5 extra\n2x\n)\n;\n",
                            "volScalarField"),
              "[hexahedron: 1.500000 2.000000][quad: nan nan nan nan nan nan nan nan nan nan]");
}

TEST(OpenFoamTextFields, NonuniformWithAShortList) {
    const FoamCase c;
    EXPECT_EQ(field_outcome(c, "p",
                            "internalField   nonuniform List<scalar>\n2\n(\n1.5\n)\n;\n",
                            "volScalarField"),
              "[hexahedron: 1.500000 0.000000][quad: nan nan nan nan nan nan nan nan nan nan]");
}

TEST(OpenFoamTextFields, NonuniformCrlf) {
    const FoamCase c;
    EXPECT_EQ(field_outcome(c, "p",
                            "internalField   nonuniform List<scalar>\r\n2\r\n(\r\n1.5\r\n2.5\r\n)\r\n;\r\n",
                            "volScalarField"),
              "[hexahedron: 1.500000 2.500000][quad: nan nan nan nan nan nan nan nan nan nan]");
}

TEST(OpenFoamTextFields, NeitherUniformNorNonuniform) {
    const FoamCase c;
    EXPECT_EQ(field_outcome(c, "p", "internalField   calculated;\n", "volScalarField"), "ReadError: OpenFOAM: internalField is neither uniform nor nonuniform: <case>/0/p");
}

TEST(OpenFoamTextFields, NoInternalField) {
    const FoamCase c;
    EXPECT_EQ(field_outcome(c, "p", "boundaryField {}\n", "volScalarField"), "ReadError: OpenFOAM: field file has no internalField: <case>/0/p");
}

// ---- a polyMesh whose lists disagree is refused, not read past -----------

TEST(OpenFoamCheck, AFaceNamingAMissingPoint) {
    const FoamCase c;
    edit(c, "faces", "4(0 6 9 3)", "4(0 6 9 99)");
    EXPECT_EQ(outcome(c.File()), "ReadError: OpenFOAM: face 1 names point 99, but the mesh has 12 points");
}

TEST(OpenFoamCheck, AFaceNamingANegativePoint) {
    const FoamCase c;
    edit(c, "faces", "4(0 6 9 3)", "4(0 6 9 -3)");
    EXPECT_EQ(outcome(c.File()), "ReadError: OpenFOAM: face 1 names point -3, but the mesh has 12 points");
}

TEST(OpenFoamCheck, AnOwnerListLongerThanTheFaces) {
    const FoamCase c;
    edit(c, "owner", "\n)\n", "\n0\n)\n");
    EXPECT_EQ(outcome(c.File()), "ReadError: OpenFOAM: owner lists 12 faces, but the mesh has 11");
}

TEST(OpenFoamCheck, ANegativeOwner) {
    const FoamCase c;
    edit(c, "owner", "(\n0", "(\n-4");
    EXPECT_EQ(outcome(c.File()), "ReadError: OpenFOAM: face 0 has owner -4");
}

TEST(OpenFoamCheck, ACellIdBeyondWhatTheListsCanAccountFor) {
    const FoamCase c;
    edit(c, "owner", "(\n0", "(\n4000000000");
    EXPECT_EQ(outcome(c.File()), "ReadError: OpenFOAM: face 0 has owner 4000000000");
}

TEST(OpenFoamCheck, ANeighbourBelowMinusOne) {
    const FoamCase c;
    edit(c, "neighbour", "(\n1", "(\n-2");
    EXPECT_EQ(outcome(c.File()), "ReadError: OpenFOAM: face 0 has neighbour -2");
}

TEST(OpenFoamCheck, ANeighbourListLongerThanTheFaces) {
    const FoamCase c;
    std::string body = "\n";
    for (int i = 0; i < 12; ++i)
        body += "1\n";
    edit(c, "neighbour", "\n1\n)\n", body + ")\n");
    EXPECT_EQ(outcome(c.File()), "ReadError: OpenFOAM: neighbour lists 12 faces, but the mesh has 11");
}

TEST(OpenFoamCheck, APatchReachingPastTheLastFace) {
    const FoamCase c;
    std::string t = c.Get("boundary");
    const std::size_t at = t.find("nFaces");
    ASSERT_NE(at, std::string::npos);
    const std::size_t semi = t.find(';', at);
    t.replace(at, semi - at, "nFaces 500");
    c.Put("boundary", t);
    const std::string out = outcome(c.File());
    EXPECT_EQ(out.rfind("ReadError: OpenFOAM: patch '", 0), 0u) << out;
    EXPECT_NE(out.find("but the mesh has 11 faces"), std::string::npos) << out;
}

TEST(OpenFoamCheck, AValidCaseStillReads) {
    const FoamCase c;
    EXPECT_EQ(outcome(c.File()).rfind("P12 ", 0), 0u);
}

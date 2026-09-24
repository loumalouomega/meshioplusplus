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
// MSC Marc input decks and formatted post files, read by the C++ core directly
// (the Python suite, tests/python/test_marc.py, runs both engines on the
// generated fixtures).

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/marc.hpp"
#include "meshioplusplus/region.hpp"
#include "meshioplusplus/registry.hpp"

namespace {

std::string write_temp(const std::string& rContents, const char* pSuffix) {
    const std::string path = mt::temp_path(pSuffix);
    std::ofstream(path, std::ios::binary) << rContents;
    return path;
}

// One 20-node brick (type 21) in EXTENDED fields with Mentat's reals (no
// exponent letter, touching fields), a brick with collapsed nodes (a wedge),
// an element of a type meshio++ does not read, and sets.
const char* const kDeck = R"(title               one hex20
$....a comment
extended
sizing                       0         3        22         0
end
connectivity
         3         0         1
        10        21         1         2         3         4         5         6         7         8         9        10        11        12        13        14
        15        16        17        18        19        20
        20         7         2        21        21         3         6        22        22         7
        30       116         1         2         3         4
coordinates
         3        22         0         1
         1 0.000000000000000+0 0.000000000000000+0 0.000000000000000+0
         2 1.000000000000000+0 0.000000000000000+0 0.000000000000000+0
         3 1.000000000000000+0 1.000000000000000+0 0.000000000000000+0
         4 0.000000000000000+0 1.000000000000000+0 0.000000000000000+0
         5 0.000000000000000+0 0.000000000000000+0 1.000000000000000+0
         6 1.000000000000000+0 0.000000000000000+0 1.000000000000000+0
         7 1.000000000000000+0 1.000000000000000+0 1.000000000000000+0
         8 0.000000000000000+0 1.000000000000000+0 1.000000000000000+0
         9 5.000000000000000-1 0.000000000000000+0 0.000000000000000+0
        10 1.000000000000000+0 5.000000000000000-1 0.000000000000000+0
        11 5.000000000000000-1 1.000000000000000+0 0.000000000000000+0
        12 0.000000000000000+0 5.000000000000000-1 0.000000000000000+0
        13 5.000000000000000-1 0.000000000000000+0 1.000000000000000+0
        14 1.000000000000000+0 5.000000000000000-1 1.000000000000000+0
        15 5.000000000000000-1 1.000000000000000+0 1.000000000000000+0
        16 0.000000000000000+0 5.000000000000000-1 1.000000000000000+0
        17 0.000000000000000+0 0.000000000000000+0 5.000000000000000-1
        18 1.000000000000000+0 0.000000000000000+0 5.000000000000000-1
        19 1.000000000000000+0 1.000000000000000+0 5.000000000000000-1
        20 0.000000000000000+0 1.000000000000000+0 5.000000000000000-1
        21 2.000000000000000+0 0.000000000000000+0 0.000000000000000+0
        22 2.000000000000000+0 0.000000000000000+0 1.000000000000000+0
         1-1.000000000000000-1 0.000000000000000+0 0.000000000000000+0
define              node                set                 bottom
         1 to         4   c
         9 to        12
define              node                set                 corners
         1 to         8 except bottom
define              element             set                 bricks
        10        20        30
isotropic
         1
bricks
end option
)";

// A post file of one 8-node brick: one increment, the equivalent stress and the
// stress tensor at one integration point, displacements and a node set.
const char* const kPost = R"(=beg=50100 (Analysis Title)
          job
=end=
=beg=50200 (Analysis Verification Data)
            7            8            1            3            1            9
            1            0            3            8            2            0
            0           12            0            1            0            0
            0            0            0            0            0            0
            0            0            0            0            0            0
=end=
=beg=50600 (Element Variable Postcodes)
           17
          311
          312
          313
          314
          315
          316
=end=
=beg=50700 (Element Connectivities)
            5            7            8            1            2            3
            4            5            6            7            8
=end=
=beg=50800 (Nodal Coordinates)
            1 0.000000E+00 0.000000E+00 0.000000E+00
            2 0.100000E+01 0.000000E+00 0.000000E+00
            3 0.100000E+01 0.100000E+01 0.000000E+00
            4 0.000000E+00 0.100000E+01 0.000000E+00
            5 0.000000E+00 0.000000E+00 0.100000E+01
            6 0.100000E+01 0.000000E+00 0.100000E+01
            7 0.100000E+01 0.100000E+01 0.100000E+01
            8 0.000000E+00 0.100000E+01 0.100000E+01
=end=
=beg=51000 (Nodal Codes and Transformation ID)
            3            3            3            3            3            3
            3            3
=end=
=beg=51301 (Set Definitions)
            1
base
            2            1
            1            2
=end=
****
=beg=51701 (Integer Increment Verification Data)
            0            3            0          102            1            0
            0            0            0            0            0            0
=end=
=beg=51801 (Real Increment Verification Data)
           24
 0.250000E+01 0.000000E+00 0.000000E+00 0.000000E+00 0.000000E+00 0.000000E+00
 0.000000E+00 0.000000E+00 0.000000E+00 0.000000E+00 0.000000E+00 0.000000E+00
 0.000000E+00 0.000000E+00 0.000000E+00 0.000000E+00 0.000000E+00 0.000000E+00
 0.000000E+00 0.000000E+00 0.000000E+00 0.000000E+00 0.000000E+00 0.000000E+00
=end=
=beg=52300 (Element Integration Point Values)
 0.288818E+03 0.480487E+02 0.276662E+02-0.630963E+02 0.327900E+02-0.152387E+03
 0.267586E+01
=end=
=beg=52401 (Nodal Results)
            1            3
Displacement
            1            0            0            3            0            0
           -1            0            0            0            0            0
 0.100000E+01 0.000000E+00 0.000000E+00 0.200000E+01 0.000000E+00 0.000000E+00
 0.300000E+01 0.000000E+00 0.000000E+00 0.400000E+01 0.000000E+00 0.000000E+00
 0.500000E+01 0.000000E+00 0.000000E+00 0.600000E+01 0.000000E+00 0.000000E+00
 0.700000E+01 0.000000E+00 0.000000E+00 0.800000E+01-0.100000E+01 0.000000E+00
=end=
----
++++
)";

}  // namespace

TEST(Marc, ContentCheck) {
    EXPECT_TRUE(meshioplusplus::is_marc_deck(kDeck));
    EXPECT_TRUE(meshioplusplus::is_marc_deck("$ comment\nsizing, 0, 1,\nconnectivity\n"));
    EXPECT_FALSE(meshioplusplus::is_marc_deck("TITLE = \"x\"\nVARIABLES = x\nend\n"));
    EXPECT_FALSE(meshioplusplus::is_marc_deck("title  x\nsizing 0\n"));  // no END yet
    EXPECT_FALSE(meshioplusplus::is_marc_deck("ZONE T=\"a\"\n"));
}

TEST(Marc, ExtendedDeck) {
    using meshioplusplus::detail::read_double;
    using meshioplusplus::detail::read_int;
    const std::string path = write_temp(kDeck, ".dat");
    const meshioplusplus::Mesh mesh = meshioplusplus::read_marc(path);
    ASSERT_EQ(mesh.NumCellBlocks(), 2u);
    EXPECT_EQ(mesh.Cells(0).Type(), "hexahedron20");
    EXPECT_EQ(mesh.Cells(1).Type(), "wedge");  // 2 21 21 3 / 6 22 22 7
    EXPECT_EQ(mesh.NumPoints(), 22u);
    for (std::size_t k = 0; k < 20; ++k)
        EXPECT_EQ(read_int(mesh.Cells(0).Conn(), k), static_cast<std::int64_t>(k));
    EXPECT_EQ(read_int(mesh.CellData("marc:element", 0), 0), 10);
    EXPECT_EQ(read_int(mesh.CellData("marc:type", 1), 0), 7);
    // Node 1's second definition wins; 5.000000000000000-1 is 0.5.
    EXPECT_DOUBLE_EQ(read_double(mesh.Points(), 0), -0.1);
    EXPECT_DOUBLE_EQ(read_double(mesh.Points(), 8 * 3), 0.5);
    const auto region = [&](const char* pName, meshioplusplus::RegionKind Kind) {
        const std::size_t k = mesh.FindRegion(pName, Kind);
        EXPECT_NE(k, meshioplusplus::Mesh::npos) << pName;
        std::vector<std::int64_t> out;
        const auto& entries = mesh.Region(k).mEntries;
        for (std::size_t i = 0; i < entries.Size(); ++i)
            out.push_back(read_int(entries, i));
        return out;
    };
    using meshioplusplus::RegionKind;
    EXPECT_EQ(region("bottom", RegionKind::Point),
              (std::vector<std::int64_t>{0, 1, 2, 3, 8, 9, 10, 11}));
    EXPECT_EQ(region("corners", RegionKind::Point), (std::vector<std::int64_t>{4, 5, 6, 7}));
    // Element 30 has no cell: the set keeps the other two.
    EXPECT_EQ(region("bricks", RegionKind::Cell), (std::vector<std::int64_t>{0, 1}));
    EXPECT_EQ(mesh.Region(mesh.FindRegion("bricks", RegionKind::Cell)).mDim, 3);
    // resolve_format: an existing Marc deck named .dat is Marc's.
    EXPECT_EQ(meshioplusplus::resolve_format(path, ""), "marc");
    std::remove(path.c_str());
}

TEST(Marc, DeckRefusals) {
    const std::string short_element =
        write_temp("title x\nend\nconnectivity\n    0\n    1    7    1    2\n", ".dat");
    EXPECT_THROW(meshioplusplus::read_marc(short_element), meshioplusplus::ReadError);
    std::remove(short_element.c_str());
    const std::string undefined = write_temp(
        "title x\nend\nconnectivity\n    0\n    1  134    1    2    3    9\ncoordinates\n"
        "    3    3\n    1        0.        0.        0.\n    2        1.        0.        0.\n"
        "    3        0.        1.        0.\n",
        ".dat");
    EXPECT_THROW(meshioplusplus::read_marc(undefined), meshioplusplus::ReadError);
    std::remove(undefined.c_str());
    const std::string tecplot = write_temp("TITLE = \"x\"\nVARIABLES = x y\n", ".dat");
    EXPECT_THROW(meshioplusplus::read_marc(tecplot), meshioplusplus::ReadError);
    EXPECT_EQ(meshioplusplus::resolve_format(tecplot, ""), "tecplot");
    std::remove(tecplot.c_str());
}

TEST(Marc, PostFile) {
    using meshioplusplus::detail::read_double;
    const std::string path = write_temp(kPost, ".t19");
    const auto meta = meshioplusplus::read_marc_t19_metadata(path);
    EXPECT_EQ(meta.mTimeValues, (std::vector<double>{2.5}));
    const meshioplusplus::Mesh mesh = meshioplusplus::read_marc_t19(path);
    ASSERT_EQ(mesh.NumCellBlocks(), 1u);
    EXPECT_EQ(mesh.Cells(0).Type(), "hexahedron");
    EXPECT_DOUBLE_EQ(read_double(mesh.FieldData("meshio:time"), 0), 2.5);
    EXPECT_EQ(meshioplusplus::detail::read_int(mesh.FieldData("marc:increment"), 0), 3);
    // One integration point: the point axis is dropped. The tensor codes
    // 311..316 are xx yy zz xy yz zx: their von Mises stress is code 17's.
    const auto& stress = mesh.CellData("Stress", 0);
    ASSERT_EQ(stress.Shape(), (std::vector<std::size_t>{1, 6}));
    double t[6];
    for (std::size_t c = 0; c < 6; ++c)
        t[c] = read_double(stress, c);
    const double mises =
        std::sqrt(0.5 * ((t[0] - t[1]) * (t[0] - t[1]) + (t[1] - t[2]) * (t[1] - t[2]) +
                         (t[2] - t[0]) * (t[2] - t[0])) +
                  3.0 * (t[3] * t[3] + t[4] * t[4] + t[5] * t[5]));
    EXPECT_NEAR(mises, read_double(mesh.CellData("Equivalent Von Mises Stress", 0), 0), 1e-3);
    const auto& u = mesh.PointData("Displacement");
    ASSERT_EQ(u.Shape(), (std::vector<std::size_t>{8, 3}));
    EXPECT_DOUBLE_EQ(read_double(u, 7 * 3), 8.0);
    EXPECT_DOUBLE_EQ(read_double(u, 7 * 3 + 1), -1.0);
    EXPECT_NE(mesh.FindRegion("base", meshioplusplus::RegionKind::Point),
              meshioplusplus::Mesh::npos);

    meshioplusplus::ReadOptions options;
    options.mDataArrays = std::vector<std::string>{"Displacement"};
    const meshioplusplus::Mesh only = meshioplusplus::read_marc_t19(path, options);
    EXPECT_FALSE(only.HasCellData("Stress"));
    EXPECT_TRUE(only.HasPointData("Displacement"));
    options = {};
    options.mTimeStep = 1;
    EXPECT_THROW(meshioplusplus::read_marc_t19(path, options), meshioplusplus::ReadError);

    // A remeshing increment is refused.
    std::string text = kPost;
    const std::string head = "(Integer Increment Verification Data)\n";
    const std::size_t at = text.find(head) + head.size();
    text.replace(at, 13, "            1");
    const std::string remesh = write_temp(text, ".t19");
    EXPECT_THROW(meshioplusplus::read_marc_t19(remesh), meshioplusplus::ReadError);
    std::remove(remesh.c_str());
    std::remove(path.c_str());
}

// License: MIT License, meshio++ default license: LICENSE

#include <filesystem>
#include <cmath>
#include <string>

#include <gtest/gtest.h>

#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/abaqus.hpp"
#include "meshioplusplus/formats/dex.hpp"
#include "meshioplusplus/formats/nastran.hpp"
#include "meshioplusplus/formats/flac3d.hpp"
#include "meshioplusplus/formats/medit.hpp"
#include "meshioplusplus/formats/femap.hpp"
#include "meshioplusplus/formats/gid.hpp"
#include "meshioplusplus/formats/gmsh.hpp"
#include "meshioplusplus/formats/lsdyna.hpp"
#include "meshioplusplus/formats/obj_off.hpp"
#include "meshioplusplus/formats/pcd.hpp"
#include "meshioplusplus/formats/radioss.hpp"
#include "meshioplusplus/formats/unv.hpp"
#include "meshioplusplus/formats/xyz.hpp"

namespace {

template <class Reader>
mt::Mesh tiv_read(const std::string& rText, const std::string& rSuffix, Reader&& read) {
    const std::string path = mt::temp_path(rSuffix);
    {
        auto out = meshioplusplus::detail::make_classic_ofstream(path, std::ios::binary);
        out << rText;
    }
    auto mesh = read(path);
    std::filesystem::remove(path);
    {
        auto out = meshioplusplus::detail::make_classic_ofstream(path, std::ios::binary);
        out << std::string(rText.size(), '?');
    }
    std::filesystem::remove(path);
    return mesh;
}

}  // namespace

TEST(TextIoViews, OffStreamExtractionAndOwnership) {
    for (const std::size_t padding : {0u, 16u << 20}) {
        const auto mesh = tiv_read("OFF\r\n#" + std::string(padding, ' ') +
                                       "\r\n\t\r\n3 1 0 ignored\r\n"
                                       "+0 -0 0\v1e0 0 0\f0 1.0 0\r\n+3 0 1 2 trailing",
                                   ".off", meshioplusplus::read_off);
        EXPECT_EQ(mesh.NumPoints(), 3u);
        ASSERT_EQ(mesh.NumCellBlocks(), 1u);
        EXPECT_EQ(mesh.Cells(0).Type(), "triangle");
        EXPECT_DOUBLE_EQ(mesh.Points().As<double>()[3], 1.0);
        EXPECT_EQ(mesh.Cells(0).Conn().As<std::int64_t>()[2], 2);
    }
}

TEST(TextIoViews, MeditPrefixesCommentsAndOwnership) {
    for (const std::size_t padding : {0u, 16u << 20}) {
        const auto mesh = tiv_read("#" + std::string(padding, ' ') +
                                       "\nMeshVersionFormatted 1\r\nDimension 3\n"
                                       "Vertices\n1\n+1suffix 2e0suffix nope +7suffix\nEnd",
                                   ".mesh", meshioplusplus::read_medit_ascii);
        ASSERT_EQ(mesh.NumPoints(), 1u);
        EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(mesh.Points(), 0), 1.0);
        EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(mesh.Points(), 1), 2.0);
        EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(mesh.Points(), 2), 0.0);
        EXPECT_EQ(mesh.PointData("medit:ref").As<std::int64_t>()[0], 7);
    }
}

TEST(TextIoViews, XyzDelimiterAndOwnership) {
    for (const std::size_t padding : {0u, 16u << 20}) {
        const auto mesh =
            tiv_read("#" + std::string(padding, ' ') + "\r\n1 :: +2 :: 3e0 ::", ".xyz",
                     [](const std::string& rPath) {
                         meshioplusplus::XyzReadOptions options;
                         options.mDelimiter = "::";
                         return meshioplusplus::read_xyz(rPath, options);
                     });
        ASSERT_EQ(mesh.NumPoints(), 1u);
        EXPECT_DOUBLE_EQ(mesh.Points().As<double>()[0], 1.0);
        EXPECT_DOUBLE_EQ(mesh.Points().As<double>()[1], 2.0);
        EXPECT_DOUBLE_EQ(mesh.Points().As<double>()[2], 3.0);
    }
}

TEST(TextIoViews, UnvShortAndLongFortranRealPrefixes) {
    // Free-format nodes exercise both the stack field and owned long-token path.
    for (const std::size_t zeros : {0u, 80u}) {
        const auto mesh = tiv_read(
            "    -1\n 2411\n1 1 1 11\n1." + std::string(zeros, '0') + "D+00suffix 2d0 +3.0\n    -1",
            ".unv", [](const std::string& rPath) { return meshioplusplus::read_unv(rPath); });
        ASSERT_EQ(mesh.NumPoints(), 1u);
        EXPECT_DOUBLE_EQ(mesh.Points().As<double>()[0], 1.0);
        EXPECT_DOUBLE_EQ(mesh.Points().As<double>()[1], 2.0);
        EXPECT_DOUBLE_EQ(mesh.Points().As<double>()[2], 3.0);
    }
}

TEST(TextIoViews, PcdBoundedLongNumbersAndOwnership) {
    const auto mesh = tiv_read(
        "VERSION .7\nFIELDS x y z\nSIZE 8 8 8\nTYPE F F F\n"
        "COUNT 1 1 1\nWIDTH 1\nHEIGHT 1\nPOINTS 1\nDATA ascii\n1." +
            std::string(80, '0') + "e0 +2 3",
        ".pcd", [](const std::string& rPath) { return meshioplusplus::read_pcd(rPath); });
    ASSERT_EQ(mesh.NumPoints(), 1u);
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(mesh.Points(), 0), 1.0);
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(mesh.Points(), 1), 2.0);
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(mesh.Points(), 2), 3.0);
}

TEST(TextIoViews, FemapFieldsAndLongNumbersOwnTheirResults) {
    for (const std::size_t padding : {0u, 16u << 20}) {
        const auto mesh = tiv_read(
            std::string(padding, ' ') +
                "\n -1\n100\n<NULL>\n9.3,\n-1\n -1\n403\n"
                "1,0,0,1,46,0,0,0,0,0,0,1." +
                std::string(80, '0') + "e0,2,3,0,1,\n-1",
            ".neu", [](const std::string& rPath) { return meshioplusplus::read_femap(rPath); });
        ASSERT_EQ(mesh.NumPoints(), 1u);
        EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(mesh.Points(), 0), 1.0);
        EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(mesh.Points(), 1), 2.0);
        EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(mesh.Points(), 2), 3.0);
    }
}

TEST(TextIoViews, RadiossEngineBoundedNumberTokens) {
    const auto mesh = tiv_read(
        "#RADIOSS ENGINE\n/ANIM/DT\n1." + std::string(80, '0') + "e0 +2e0 junk 3suffix\n/END",
        "_0001.rad", meshioplusplus::read_radioss);
    ASSERT_TRUE(mesh.HasFieldData("radioss:engine:ANIM/DT"));
    const auto& values = mesh.FieldData("radioss:engine:ANIM/DT");
    ASSERT_EQ(values.Size(), 2u);
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(values, 0), 1.0);
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(values, 1), 2.0);
}

TEST(TextIoViews, GidQuotedNamesAndNumericPrefixes) {
    const auto mesh = tiv_read(
        "MESH \"a long quoted mesh name\" dimension 3 ElemType Point Nnode 1\r\n"
        "Coordinates\r\n+1suffix +1suffix 2e0suffix 3\r\nEnd Coordinates\r\n"
        "Elements\r\n+1suffix 1suffix +7suffix\r\nEnd Elements",
        ".post.msh", [](const std::string& rPath) { return meshioplusplus::read_gid(rPath); });
    ASSERT_EQ(mesh.NumPoints(), 1u);
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(mesh.Points(), 0), 1.0);
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(mesh.Points(), 1), 2.0);
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(mesh.Points(), 2), 3.0);
}

TEST(TextIoViews, GmshFloatSpelledIdsAndHexCoordinates) {
    const auto mesh = tiv_read(
        "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n2\n"
        "2.0 1e0 -0 0x1p1\n1e3 3 4 5\n$EndNodes\n"
        "$Elements\n0\n$EndElements",
        ".msh", [](const std::string& rPath) { return meshioplusplus::read_gmsh(rPath); });
    ASSERT_EQ(mesh.NumPoints(), 2u);
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(mesh.Points(), 0), 1.0);
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(mesh.Points(), 2), 2.0);
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(mesh.Points(), 3), 3.0);
}

TEST(TextIoViews, AbaqusTokensRetainNativePrefixSemantics) {
    const std::string body = "*NODE\n+10suffix,1." + std::string(80, '0') + "suffix,2" +
                             std::string("\0ignored,junk,", 14);
    const auto mesh = tiv_read(body, ".inp", meshioplusplus::read_abaqus);
    ASSERT_EQ(mesh.NumPoints(), 1u);
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(mesh.Points(), 0), 1.0);
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(mesh.Points(), 1), 2.0);
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(mesh.Points(), 2), 0.0);
}

TEST(TextIoViews, DexNormalizesShortAndLongFortranExponents) {
    const std::string body =
        "# NB_POINTS = 1\n# NB_COMP = 1\n1D0 +2d0 -0D0 3." + std::string(100, '0') + "D0suffix";
    const auto mesh = tiv_read(body, ".dex", meshioplusplus::read_dex);
    ASSERT_EQ(mesh.NumPoints(), 1u);
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(mesh.Points(), 0), 1.0);
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(mesh.Points(), 1), 2.0);
    EXPECT_TRUE(std::signbit(meshioplusplus::detail::read_double(mesh.Points(), 2)));
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(mesh.PointData("dex:field"), 0), 3.0);
}

TEST(TextIoViews, NastranMixedContinuationOwnsNoncontiguousJoins) {
    std::string body = "BEGIN BULK\n";
    for (int i = 1; i <= 300; ++i) {
        const std::string id = std::to_string(i);
        body += "GRID*   " + id + std::string(16 - id.size(), ' ') +
                "0               1.25            2               \n*CONT,3,.\n";
    }
    body += "ENDDATA";
    const auto mesh = tiv_read(
        body, ".bdf", [](const std::string& rPath) { return meshioplusplus::read_nastran(rPath); });
    ASSERT_EQ(mesh.NumPoints(), 300u);
    for (std::size_t i = 0; i < 300; ++i) {
        EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(mesh.Points(), 3 * i), 1.25);
        EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(mesh.Points(), 3 * i + 1), 2.0);
        EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(mesh.Points(), 3 * i + 2), 3.0);
    }
}

TEST(TextIoViews, Flac3dReusedTokensKeepGroupNamesAndPrefixes) {
    const auto mesh = tiv_read(
        "*FLAC3D\r\nG 1suffix 0 0 0\r\nG 2suffix 1 0 0\r\nG 3suffix 0 1 0\r\n"
        "F T3 +1suffix 1suffix 2suffix 3suffix\r\n"
        "FGROUP 'a long quoted group name' SLOT 'a long owned slot'\r\n+1suffix",
        ".f3grid", meshioplusplus::read_flac3d);
    ASSERT_EQ(mesh.NumPoints(), 3u);
    ASSERT_EQ(mesh.NumRegions(), 1u);
    EXPECT_EQ(mesh.Region(0).mName, "face:a long quoted group name:a long owned slot");
    ASSERT_EQ(mesh.Region(0).NumEntries(), 1u);
    EXPECT_EQ(mesh.Region(0).Entries()[0], 0);
}

TEST(TextIoViews, LsdynaLineViewsKeepFinalEmptyRecord) {
    for (const auto ending : {"", "\n", "\r\n"}) {
        const auto mesh = tiv_read("*KEYWORD\r\n*NODE\r\n1,1,2,3\r\n*END" + std::string(ending),
                                   ".k", meshioplusplus::read_lsdyna);
        ASSERT_EQ(mesh.NumPoints(), 1u);
        EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(mesh.Points(), 0), 1.0);
        EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(mesh.Points(), 1), 2.0);
        EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(mesh.Points(), 2), 3.0);
    }
}

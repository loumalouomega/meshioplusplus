// SPDX-License-Identifier: MIT
#include <cstdio>
#include <string>
#include <gtest/gtest.h>
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/registry.hpp"

namespace {
std::string vtkpieces_temp(const std::string& rText, const std::string& rFormat) {
    const auto path = mt::temp_path("." + rFormat);
    auto out = meshioplusplus::detail::make_classic_ofstream(path, std::ios::binary);
    out << rText;
    return path;
}
}  // namespace

TEST(VtkPieces, NativeRegistryReadsAllSerialXmlPieces) {
    using namespace meshioplusplus;
    for (const std::string format : {"vtp", "vts", "vtr", "vti"}) {
        const bool poly = format == "vtp";
        const std::string type = poly              ? "PolyData"
                                 : format == "vts" ? "StructuredGrid"
                                 : format == "vtr" ? "RectilinearGrid"
                                                   : "ImageData";
        auto array = [](const std::string& rName, const std::string& rValues, int Nc = 1) {
            return "<DataArray type='Int64' Name='" + rName + "' NumberOfComponents='" +
                   std::to_string(Nc) + "' format='ascii'>" + rValues + "</DataArray>";
        };
        std::string xml =
            "<VTKFile type='" + type + "'><" + type + (poly ? ">" : " WholeExtent='0 2 0 1 0 1'>");
        for (int part = 0; part < 2; ++part) {
            const std::string extent =
                std::to_string(part) + " " + std::to_string(part + 1) + " 0 1 0 1";
            xml += "<Piece " +
                   (poly ? std::string("NumberOfPoints='4'") : "Extent='" + extent + "'") + ">";
            if (poly || format == "vts")
                xml += "<Points>" +
                       array("Points",
                             poly ? "0 0 0 1 0 0 0 1 0 1 1 0"
                                  : "0 0 0 1 0 0 0 1 0 1 1 0 0 0 1 1 0 1 0 1 1 1 1 1",
                             3) +
                       "</Points>";
            if (poly)
                xml += "<Polys>" + array("connectivity", "0 1 3 2") + array("offsets", "4") +
                       "</Polys>";
            if (format == "vtr")
                xml += "<Coordinates>" +
                       array("X", std::to_string(part) + " " + std::to_string(part + 1)) +
                       array("Y", "0 1") + array("Z", "0 1") + "</Coordinates>";
            xml += "<CellData>" + array("material", std::to_string(10 + part)) + "</CellData>";
            xml += "<FieldData>" + array("region:cell:body", "0") + "</FieldData></Piece>";
        }
        xml += "</" + type + "></VTKFile>";
        const auto path = vtkpieces_temp(xml, format);
        const auto mesh = registry_read(path, format, {});
        EXPECT_EQ(mesh.NumPoints(), poly ? 8 : 16);
        ASSERT_EQ(mesh.NumCellBlocks(), 1);
        EXPECT_EQ(mesh.Cells(0).NumCells(), 2);
        EXPECT_EQ(detail::read_int(mesh.CellData("material", 0), 1), 11);
        EXPECT_GE(detail::read_int(mesh.Cells(0).Conn(), mesh.Cells(0).NodesPerCell()),
                  poly ? 4 : 8);
        ASSERT_EQ(mesh.NumRegions(), 1);
        EXPECT_EQ(mesh.Region(0).mEntries.Size(), 2);
        ReadOptions opts;
        opts.mPointsOnly = true;
        const auto geometry = registry_read(path, format, opts);
        EXPECT_EQ(geometry.NumCellData(), 0);
        EXPECT_EQ(geometry.NumRegions(), 1);
        const auto meta = registry_read_metadata(path, format, {});
        EXPECT_EQ(meta.mNumPoints, mesh.NumPoints());
        ASSERT_EQ(meta.mCellBlocks.size(), 1);
        EXPECT_EQ(meta.mCellBlocks[0].mNumCells, 2);
        std::remove(path.c_str());
    }
}

TEST(VtkPieces, LegacyStructuredAttributesAndAxisPermutations) {
    using namespace meshioplusplus;
    for (const std::string dims : {"3 2 2", "3 2 1", "1 3 2", "1 1 4"}) {
        const auto path = vtkpieces_temp(
            "# vtk DataFile Version 5.1\nindependent\nASCII\n"
            "DATASET STRUCTURED_POINTS\nDIMENSIONS " +
                dims + "\nORIGIN 0 0 0\nSPACING 1 2 3\n",
            "vtk");
        const auto mesh = registry_read(path, "vtk", {});
        ASSERT_EQ(mesh.NumCellBlocks(), 1);
        const auto block = mesh.Cells(0);
        EXPECT_EQ(block.Type(), dims == "3 2 2" ? "hexahedron" : dims == "1 1 4" ? "line" : "quad");
        EXPECT_EQ(detail::read_int(block.Conn(), 0), 0);
        std::remove(path.c_str());
    }
}

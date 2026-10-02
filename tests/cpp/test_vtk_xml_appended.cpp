// SPDX-License-Identifier: MIT
// Independent byte framing, exercised through native registry paths (no fallback).
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/detail/vtu_binary.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/registry.hpp"

namespace {

std::string vtk_appended_fixture(const std::string& rFormat, bool Base64, bool UInt64,
                                 bool BigEndian) {
    const std::string type = rFormat == "vtp"   ? "PolyData"
                             : rFormat == "vts" ? "StructuredGrid"
                             : rFormat == "vtr" ? "RectilinearGrid"
                                                : "ImageData";
    std::vector<unsigned char> payload;
    auto pack = [&](std::uint64_t value, std::size_t size, std::vector<unsigned char>& rOut) {
        for (std::size_t i = 0; i < size; ++i)
            rOut.push_back(
                static_cast<unsigned char>(value >> (8 * (BigEndian ? size - 1 - i : i))));
    };
    auto array = [&](const std::string& rName, const std::vector<std::int64_t>& rValues,
                     int Nc = 0) {
        const auto offset = payload.size();
        std::vector<unsigned char> header, body;
        pack(rValues.size() * 8, UInt64 ? 8 : 4, header);
        for (const auto value : rValues)
            pack(static_cast<std::uint64_t>(value), 8, body);
        if (Base64) {
            const auto text = meshioplusplus::detail::b64encode(header.data(), header.size()) +
                              meshioplusplus::detail::b64encode(body.data(), body.size());
            payload.insert(payload.end(), text.begin(), text.end());
        } else {
            payload.insert(payload.end(), header.begin(), header.end());
            payload.insert(payload.end(), body.begin(), body.end());
        }
        return "<DataArray type='Int64' Name='" + rName + "' format='appended' offset='" +
               std::to_string(offset) + "'" +
               (Nc ? " NumberOfComponents='" + std::to_string(Nc) + "'" : "") + "/>";
    };
    std::string xml = "<VTKFile type='" + type + "' byte_order='" +
                      (BigEndian ? "BigEndian" : "LittleEndian") + "' header_type='" +
                      (UInt64 ? "UInt64" : "UInt32") + "'><" + type;
    const std::string extent = "0 1 0 1 0 1";
    if (rFormat != "vtp")
        xml += " WholeExtent='" + extent + "'";
    xml += "><Piece " +
           (rFormat == "vtp" ? std::string("NumberOfPoints='4' NumberOfPolys='1'")
                             : "Extent='" + extent + "'") +
           ">";
    const std::size_t n = rFormat == "vtp" ? 4 : 8;
    std::vector<std::int64_t> values(n), points(n * 3);
    for (std::size_t i = 0; i < n; ++i) {
        values[i] = static_cast<std::int64_t>(i) - 2;
        points[i * 3] = i % 2;
        points[i * 3 + 1] = (i / 2) % 2;
        points[i * 3 + 2] = i / 4;
    }
    xml += "<PointData>" + array("tag", values) + "</PointData><CellData>" + array("mat", {19}) +
           "</CellData>";
    if (rFormat == "vtp" || rFormat == "vts")
        xml += "<Points>" + array("Points", points, 3) + "</Points>";
    if (rFormat == "vtp")
        xml += "<Polys>" + array("connectivity", {0, 1, 3, 2}) + array("offsets", {4}) + "</Polys>";
    if (rFormat == "vtr")
        xml += "<Coordinates>" + array("X", {0, 1}) + array("Y", {0, 1}) + array("Z", {0, 1}) +
               "</Coordinates>";
    xml += "</Piece></" + type + "><AppendedData encoding='" + (Base64 ? "base64" : "raw") + "'>_";
    xml.append(reinterpret_cast<const char*>(payload.data()), payload.size());
    return xml + "</AppendedData></VTKFile>";
}

std::string vtk_appended_temp(const std::string& rText, const std::string& rFormat) {
    const auto path = mt::temp_path("." + rFormat);
    auto out = meshioplusplus::detail::make_classic_ofstream(path, std::ios::binary);
    out.write(rText.data(), static_cast<std::streamsize>(rText.size()));
    return path;
}

}  // namespace

TEST(VtkXmlAppended, NativeRegistryReadsAndSummaries) {
    namespace mio = meshioplusplus;
    for (const std::string format : {"vtp", "vts", "vtr", "vti"})
        for (const bool base64 : {false, true})
            for (const bool uint64 : {false, true})
                for (const bool big_endian : {false, true}) {
                    SCOPED_TRACE(format + (base64 ? " base64" : " raw") +
                                 (uint64 ? " UInt64" : " UInt32") + (big_endian ? " BE" : " LE"));
                    const auto path = vtk_appended_temp(
                        vtk_appended_fixture(format, base64, uint64, big_endian), format);
                    const auto mesh = mio::registry_read(path, format, {});
                    const auto meta = mio::registry_read_metadata(path, format, {});
                    std::remove(path.c_str());
                    EXPECT_EQ(mesh.NumPoints(), format == "vtp" ? 4u : 8u);
                    ASSERT_EQ(mesh.NumCellBlocks(), 1u);
                    EXPECT_EQ(mesh.Cells(0).NumCells(), 1u);
                    EXPECT_EQ(meta.mNumPoints, mesh.NumPoints());
                    EXPECT_EQ(meta.NumCells(), 1u);
                    EXPECT_EQ(meta.mPointDataNames, (std::vector<std::string>{"tag"}));
                    EXPECT_EQ(mio::detail::read_int(mesh.PointData("tag"), 0), -2);
                    EXPECT_EQ(mio::detail::read_int(mesh.CellData("mat", 0), 0), 19);
                    EXPECT_EQ(mio::detail::read_double(mesh.Points(), 3), 1.0);
                }
}

TEST(VtkXmlAppended, InvalidOffsetsAndTruncationRaiseReadError) {
    namespace mio = meshioplusplus;
    for (const std::string format : {"vtp", "vts", "vtr", "vti"}) {
        for (const std::string offset : {"-1", "bogus", "18446744073709551616", "999999"}) {
            auto text = vtk_appended_fixture(format, false, false, false);
            const auto at = text.find("offset='0'");
            text.replace(at, 10, "offset='" + offset + "'");
            const auto path = vtk_appended_temp(text, format);
            EXPECT_THROW(mio::registry_read(path, format, {}), mio::ReadError);
            std::remove(path.c_str());
        }
        auto text = vtk_appended_fixture(format, false, false, false);
        text.erase(text.find("'>_", text.find("<AppendedData")) + 3);
        text.append(2, '\0');
        text += "</AppendedData></VTKFile>";
        const auto path = vtk_appended_temp(text, format);
        EXPECT_THROW(mio::registry_read(path, format, {}), mio::ReadError);
        std::remove(path.c_str());
    }
}

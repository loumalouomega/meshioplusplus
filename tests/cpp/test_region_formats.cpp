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
// Named regions in the formats that gained them in v16.26.0: XDMF <Set>s and
// the VTU/VTP <FieldData> convention (detail/region_field_data.hpp).

// System includes
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>
#include <tuple>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/formats/vtp.hpp"
#include "meshioplusplus/formats/vtu.hpp"
#include "meshioplusplus/formats/xdmf.hpp"
#include "meshioplusplus/region.hpp"
#include "mesh_fixtures.hpp"

using meshioplusplus::Mesh;
using meshioplusplus::NDArray;
using meshioplusplus::Region;
using meshioplusplus::RegionKind;

namespace {

NDArray rf_ids(const std::vector<std::int64_t>& rV, std::size_t Stride = 1) {
    NDArray a = Stride == 1 ? NDArray(meshioplusplus::DType::Int64, {rV.size()})
                            : NDArray(meshioplusplus::DType::Int64, {rV.size() / 2, 2});
    std::copy(rV.begin(), rV.end(), a.As<std::int64_t>());
    return a;
}

NDArray rf_conn(const std::vector<std::int64_t>& rV, std::size_t Cols) {
    NDArray a(meshioplusplus::DType::Int64, {rV.size() / Cols, Cols});
    std::copy(rV.begin(), rV.end(), a.As<std::int64_t>());
    return a;
}

using RegionRow = std::tuple<std::string, int, int, std::int64_t, std::vector<std::int64_t>>;

std::vector<RegionRow> rf_rows(const Mesh& rMesh) {
    std::vector<RegionRow> out;
    for (std::size_t i = 0; i < rMesh.NumRegions(); ++i) {
        const Region& r = rMesh.Region(i);
        std::vector<std::int64_t> e(r.Entries(), r.Entries() + r.NumEntries() * r.Stride());
        out.emplace_back(r.mName, static_cast<int>(r.mKind), r.mDim, r.mTag, e);
    }
    return out;
}

/// Two tetrahedra and a triangle, carrying every region kind.
Mesh rf_mesh() {
    Mesh m = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0.5, 0.5, 1}}, "tetra",
                           {{0, 1, 2, 4}, {0, 2, 3, 4}});
    m.AddCellBlock("triangle", rf_conn({0, 1, 2}, 3));
    m.AddRegion(Region("inlet", RegionKind::Point, 0, 7, rf_ids({0, 4})));
    m.AddRegion(Region("vol", RegionKind::Cell, rf_ids({1, 2})));
    m.AddRegion(Region("wall", RegionKind::Side, rf_ids({0, 3, 1, 3, 2, 0}, 2)));
    m.AddRegion(Region("empty", RegionKind::Cell, rf_ids({})));
    return m;
}

}  // namespace

TEST(RegionFormats, VtuRoundTripsEveryKind) {
    const Mesh m = rf_mesh();
    for (const bool binary : {false, true}) {
        const std::string path = mt::temp_path(".vtu");
        meshioplusplus::write_vtu(path, m, binary, binary);
        const Mesh back = meshioplusplus::read_vtu(path);
        EXPECT_EQ(rf_rows(back), rf_rows(m)) << "binary=" << binary;
        EXPECT_EQ(back.NumFieldData(), 0u) << "region arrays are not field data";
        std::error_code ec;
        std::filesystem::remove(path, ec);
    }
}

TEST(RegionFormats, VtpNumbersCellsInFileOrder) {
    // Blocks triangle, line, vertex are written Verts, Lines, Polys: the
    // region arrays name file cells, and read back onto the reordered blocks.
    Mesh m = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}}, "triangle",
                           {{0, 1, 2}, {0, 2, 3}});
    m.AddCellBlock("line", rf_conn({0, 1}, 2));
    m.AddCellBlock("vertex", rf_conn({3}, 1));
    m.AddRegion(Region("tri", RegionKind::Cell, rf_ids({1})));
    m.AddRegion(Region("dot", RegionKind::Cell, 0, 9, rf_ids({3})));
    m.AddRegion(Region("edge", RegionKind::Side, rf_ids({0, 0, 1, 2}, 2)));
    const std::string path = mt::temp_path(".vtp");
    meshioplusplus::write_vtp(path, m, true, true);
    const Mesh back = meshioplusplus::read_vtp(path);
    ASSERT_EQ(back.NumCellBlocks(), 3u);
    EXPECT_EQ(std::string(back.Cells(0).Type()), "vertex");
    // vertex 0, line 1, triangles 2 and 3.
    std::vector<RegionRow> want{
        {"dot", 1, 0, 9, {0}}, {"tri", 1, -1, -1, {3}}, {"edge", 2, -1, -1, {2, 0, 3, 2}}};
    EXPECT_EQ(rf_rows(back), want);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(RegionFormats, AMalformedRegionArrayStaysFieldData) {
    Mesh m = mt::tet_mesh();
    // A cell index past the file's cells is not a region; the writer never
    // emits one, so the file is edited by hand.
    const std::string path = mt::temp_path(".vtu");
    meshioplusplus::write_vtu(path, m, false, false);
    std::ifstream in(path);
    std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    in.close();
    const std::string grid = "<UnstructuredGrid>\n";
    text.insert(text.find(grid) + grid.size(),
                "<FieldData>\n<DataArray type=\"Int64\" Name=\"region:cell:bad\" "
                "NumberOfTuples=\"1\" format=\"ascii\">\n99\n</DataArray>\n</FieldData>\n");
    std::ofstream(path) << text;
    const Mesh back = meshioplusplus::read_vtu(path);
    EXPECT_EQ(back.NumRegions(), 0u);
    EXPECT_TRUE(back.HasFieldData("region:cell:bad"));
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(RegionFormats, XdmfSetsRoundTripEveryKind) {
    const Mesh m = rf_mesh();
    const std::string path = mt::temp_path(".xdmf");
    meshioplusplus::write_xdmf(path, m, "XML", -1);
    const Mesh back = meshioplusplus::read_xdmf(path);
    EXPECT_EQ(rf_rows(back), rf_rows(m));
    // The Side region spans a tetra face and a triangle edge: one Face set
    // and one Edge set, merged back into one region.
    std::ifstream in(path);
    const std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    EXPECT_NE(text.find("SetType=\"Face\""), std::string::npos);
    EXPECT_NE(text.find("SetType=\"Edge\""), std::string::npos);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

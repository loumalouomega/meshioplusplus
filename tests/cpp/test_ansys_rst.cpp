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
#include <cstring>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/ansys_model.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/ansys_rst.hpp"
#include "meshioplusplus/region.hpp"

namespace {

// Records as MAPDL lays them out: length in words, flags, data, trailer. The
// twin of the Python suite's synthetic file (tests/python/test_ansys_rst.py).
struct RstWriter {
    std::vector<std::int32_t> mWords;

    std::size_t Raw(const std::vector<std::int32_t>& rData, std::uint32_t Flags) {
        const std::size_t ptr = mWords.size();
        mWords.push_back(static_cast<std::int32_t>(rData.size()));
        mWords.push_back(static_cast<std::int32_t>(Flags << 24));
        mWords.insert(mWords.end(), rData.begin(), rData.end());
        mWords.push_back(0);
        return ptr;
    }
    std::size_t Ints(std::vector<std::int32_t> Values, std::size_t Size = 0) {
        if (Values.size() < Size)
            Values.resize(Size, 0);
        return Raw(Values, 0x80);
    }
    static std::vector<std::int32_t> AsWords(const std::vector<double>& rValues) {
        std::vector<std::int32_t> out(rValues.size() * 2);
        std::memcpy(out.data(), rValues.data(), rValues.size() * 8);
        return out;
    }
    std::size_t Doubles(const std::vector<double>& rValues) { return Raw(AsWords(rValues), 0x00); }
    std::size_t Shorts(std::vector<std::int16_t> Values) {
        if (Values.size() % 2)
            Values.push_back(0);
        std::vector<std::int32_t> data(Values.size() / 2);
        std::memcpy(data.data(), Values.data(), Values.size() * 2);
        return Raw(data, 0xC0);
    }
    std::size_t BsparseDoubles(const std::vector<double>& rValues) {
        std::uint32_t bits = 0;
        std::vector<double> packed;
        for (std::size_t k = 0; k < rValues.size(); ++k)
            if (rValues[k] != 0.0) {
                bits |= 1u << k;
                packed.push_back(rValues[k]);
            }
        std::vector<std::int32_t> data{static_cast<std::int32_t>(rValues.size()),
                                       static_cast<std::int32_t>(bits)};
        const auto words = AsWords(packed);
        data.insert(data.end(), words.begin(), words.end());
        return Raw(data, 0x08);
    }
    void Patch(std::size_t Ptr, std::size_t K, std::int64_t Value) {
        mWords[Ptr + 2 + K] = static_cast<std::int32_t>(static_cast<std::uint32_t>(Value));
    }
    void Save(const std::string& rPath) const {
        std::ofstream(rPath, std::ios::binary)
            .write(reinterpret_cast<const char*>(mWords.data()),
                   static_cast<std::streamsize>(mWords.size() * 4));
    }
};

constexpr double kAngles[3] = {30.0, 45.0, 60.0};  // node 3: THXY, THYZ, THZX
const std::vector<std::int32_t> kNeqv = {8, 7, 6, 5, 4, 3, 2, 1};
const double kHex[8][3] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                           {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};

std::int32_t neqv_position(std::int32_t Node) {
    for (std::size_t k = 0; k < kNeqv.size(); ++k)
        if (kNeqv[k] == Node)
            return static_cast<std::int32_t>(k) + 1;
    return 0;
}

std::string write_synthetic(bool Zlib = false, std::int32_t GlobalNodes = 0) {
    RstWriter w;
    w.Ints({12}, 100);  // the standard header: file 12
    const std::size_t header = w.Ints({12, 8, 8, 10, 4, 1, 1, 0, 2}, 80);
    w.Patch(header, 20, 1);
    w.Patch(header, 48, GlobalNodes);
    w.Patch(header, 14, static_cast<std::int64_t>(w.Ints(kNeqv)));

    const std::size_t geometry = w.Ints({0, 1, 0, 8, 1}, 80);
    w.Patch(header, 15, static_cast<std::int64_t>(geometry));
    const std::size_t ety = w.Ints({0});
    w.Patch(geometry, 20, static_cast<std::int64_t>(ety));
    const std::size_t solid185 = w.Ints({1, 185, 0}, 100);
    w.Patch(solid185, 60, 8);  // nodelm
    w.Patch(ety, 0, static_cast<std::int64_t>(solid185 - ety));
    std::size_t nodes = 0;
    for (int n = 1; n <= 8; ++n) {
        std::vector<double> values = {
            static_cast<double>(n), kHex[n - 1][0], kHex[n - 1][1], kHex[n - 1][2], 0.0, 0.0, 0.0};
        if (n == 3)
            std::copy(kAngles, kAngles + 3, values.begin() + 4);
        const std::size_t ptr = n == 1 ? w.BsparseDoubles(values) : w.Doubles(values);
        if (n == 1)
            nodes = ptr;
    }
    w.Patch(geometry, 26, static_cast<std::int64_t>(nodes));
    const std::size_t eid = w.Ints({0, 0});
    w.Patch(geometry, 28, static_cast<std::int64_t>(eid));
    const std::size_t element = w.Shorts({2, 1, 3, 4, 0, 1, 0, 0, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8});
    w.Patch(eid, 0, static_cast<std::int64_t>(element - eid));
    std::vector<std::int32_t> comp = {1, ('F' << 24) | ('A' << 16) | ('C' << 8) | 'E'};
    comp.resize(9, 0x20202020);
    comp.push_back(1);
    comp.push_back(-4);
    w.Patch(geometry, 50, static_cast<std::int64_t>(w.Ints(comp)));
    w.Patch(geometry, 48, 1);

    std::vector<std::int32_t> sets;
    for (int s = 0; s < 2; ++s) {
        const std::size_t base = w.Ints({0, 1, 8}, 150);
        w.Patch(base, 19, 4);
        const std::int32_t codes[4] = {1, 2, 3, 20};
        for (std::size_t k = 0; k < 4; ++k)
            w.Patch(base, 20 + k, codes[k]);
        std::vector<double> rows;
        std::size_t nsl = 0;
        if (s == 0) {  // every node: U = (n, 0, 0), TEMP = 1000 + n
            for (std::int32_t n : kNeqv)
                rows.insert(rows.end(), {double(n), 0.0, 0.0, 1000.0 + n});
            nsl = w.Doubles(rows);
        } else {  // nodes 5 and 2 only, node 5's TEMP undefined
            rows = {0.0, 5.0, 0.0, std::ldexp(1.0, 100), 0.0, 2.0, 0.0, 1002.0};
            nsl = w.Doubles(rows);
            w.Ints({neqv_position(5), neqv_position(2)});
        }
        w.Patch(base, 104, static_cast<std::int64_t>(nsl - base));
        sets.push_back(static_cast<std::int32_t>(base));
    }
    const std::size_t dsi = w.Ints(sets, 20);
    const std::size_t tim = w.Doubles({0.5, 1.0, 0, 0, 0, 0, 0, 0, 0, 0});
    const std::size_t lsp = w.Ints({1, 1, 1, 1, 2, 2}, 30);
    w.Patch(header, 10, static_cast<std::int64_t>(dsi));
    w.Patch(header, 11, static_cast<std::int64_t>(tim));
    w.Patch(header, 12, static_cast<std::int64_t>(lsp));
    if (Zlib)
        w.mWords[tim + 1] = static_cast<std::int32_t>(0x20u << 24);
    const std::string path = mt::temp_path(".rst");
    w.Save(path);
    return path;
}

void remove_file(const std::string& rPath) {
    std::error_code ec;
    std::filesystem::remove(rPath, ec);
}

double read_point(const meshioplusplus::Mesh& rMesh, const std::string& rName, std::size_t K) {
    return meshioplusplus::detail::read_double(rMesh.PointData(rName), K);
}

}  // namespace

TEST(AnsysRst, SyntheticModelAndFirstSet) {
    const std::string path = write_synthetic();
    const meshioplusplus::Mesh mesh = meshioplusplus::read_ansys_rst(path);
    ASSERT_EQ(mesh.NumCellBlocks(), 1u);
    EXPECT_EQ(mesh.Cells(0).Type(), "hexahedron");
    EXPECT_EQ(mesh.NumPoints(), 8u);
    for (std::size_t k = 0; k < 8; ++k)
        EXPECT_EQ(meshioplusplus::detail::read_int(mesh.Cells(0).Conn(), k),
                  static_cast<std::int64_t>(k));
    EXPECT_EQ(meshioplusplus::detail::read_int(mesh.CellData("ansys:mat", 0), 0), 2);
    EXPECT_EQ(meshioplusplus::detail::read_int(mesh.CellData("ansys:secnum", 0), 0), 4);
    const std::size_t face = mesh.FindRegion("FACE", meshioplusplus::RegionKind::Point);
    ASSERT_NE(face, meshioplusplus::Mesh::npos);
    EXPECT_EQ(mesh.Region(face).NumEntries(), 4u);
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(mesh.FieldData("meshio:time"), 0), 0.5);
    for (std::size_t p = 0; p < 8; ++p)
        EXPECT_DOUBLE_EQ(read_point(mesh, "TEMP", p), 1001.0 + static_cast<double>(p));

    // Node 3 (point 2): U = (3, 0, 0) in axes turned about Z by 30, then the new
    // X by 45, then the newest Y by 60 degrees: v = Rz Rx Ry (3, 0, 0).
    const double d = std::acos(-1.0) / 180.0;
    const double cx = std::cos(kAngles[1] * d), sx = std::sin(kAngles[1] * d);
    const double cy = std::cos(kAngles[2] * d), sy = std::sin(kAngles[2] * d);
    const double cz = std::cos(kAngles[0] * d), sz = std::sin(kAngles[0] * d);
    double v[3] = {3.0 * cy, 0.0, -3.0 * sy};                                         // Ry
    double t[3] = {v[0], cx * v[1] - sx * v[2], sx * v[1] + cx * v[2]};               // Rx
    const double expected[3] = {cz * t[0] - sz * t[1], sz * t[0] + cz * t[1], t[2]};  // Rz
    for (std::size_t c = 0; c < 3; ++c)
        EXPECT_NEAR(read_point(mesh, "U", 2 * 3 + c), expected[c], 1e-14);
    EXPECT_DOUBLE_EQ(read_point(mesh, "U", 4 * 3), 5.0);  // node 5 is not rotated
    remove_file(path);
}

TEST(AnsysRst, PartialSetAndUndefinedValues) {
    const std::string path = write_synthetic();
    meshioplusplus::ReadOptions options;
    options.mTimeStep = 1;
    const meshioplusplus::Mesh mesh = meshioplusplus::read_ansys_rst(path, options);
    EXPECT_EQ(meshioplusplus::detail::read_int(mesh.FieldData("ansys:substep"), 0), 2);
    EXPECT_DOUBLE_EQ(read_point(mesh, "TEMP", 1), 1002.0);
    EXPECT_TRUE(std::isnan(read_point(mesh, "TEMP", 4)));  // node 5: undefined
    EXPECT_TRUE(std::isnan(read_point(mesh, "TEMP", 0)));  // node 1: no solution
    EXPECT_DOUBLE_EQ(read_point(mesh, "U", 1 * 3 + 1), 2.0);
    EXPECT_DOUBLE_EQ(read_point(mesh, "U", 4 * 3 + 1), 5.0);
    EXPECT_TRUE(std::isnan(read_point(mesh, "U", 7 * 3)));

    const auto meta = meshioplusplus::read_ansys_rst_metadata(path);
    EXPECT_EQ(meta.mTimeValues, (std::vector<double>{0.5, 1.0}));
    options.mTimeStep = 2;
    EXPECT_THROW(meshioplusplus::read_ansys_rst(path, options), meshioplusplus::ReadError);
    remove_file(path);
}

TEST(AnsysRst, Refusals) {
    for (const auto& path : {write_synthetic(true), write_synthetic(false, 16)}) {
        EXPECT_THROW(meshioplusplus::read_ansys_rst(path), meshioplusplus::ReadError);
        remove_file(path);
    }
    const std::string junk = mt::temp_path(".rst");
    std::ofstream(junk, std::ios::binary) << std::string(64, '\0');
    EXPECT_THROW(meshioplusplus::read_ansys_rst(junk), meshioplusplus::ReadError);
    remove_file(junk);
}

TEST(AnsysRst, ElementCategories) {
    using meshioplusplus::detail::AnsysCategory;
    EXPECT_EQ(meshioplusplus::detail::ansys_category(186), AnsysCategory::Brick);
    EXPECT_EQ(meshioplusplus::detail::ansys_category(187), AnsysCategory::Tet);
    EXPECT_EQ(meshioplusplus::detail::ansys_category(181), AnsysCategory::Shell);
    EXPECT_EQ(meshioplusplus::detail::ansys_category(188), AnsysCategory::LinearLine);
    EXPECT_EQ(meshioplusplus::detail::ansys_category(999), AnsysCategory::Skip);
    EXPECT_EQ(meshioplusplus::detail::ansys_mesh200_category(6), AnsysCategory::Shell);
}

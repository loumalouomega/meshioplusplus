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
//                   SPDX-License-Identifier: MIT
//                   meshio++ default license: LICENSE
//
//  Main authors:    Vicente Mataix Ferrandiz
//
// External includes
#include <gtest/gtest.h>

// System includes
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <sstream>
#include <string>

// Project includes
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/registry.hpp"
#include "mesh_fixtures.hpp"

using namespace meshioplusplus;

namespace {

class IpCerrCapture {
public:
    IpCerrCapture() : mOld(std::cerr.rdbuf(mBuf.rdbuf())) {}
    ~IpCerrCapture() { std::cerr.rdbuf(mOld); }
    std::string str() const { return mBuf.str(); }

private:
    std::ostringstream mBuf;
    std::streambuf* mOld;
};

Mesh ip_mesh(std::int64_t Value, DType Dtype = DType::Int64) {
    Mesh m = mt::tri_mesh();
    NDArray a(Dtype, {m.NumPoints()});
    for (std::size_t i = 0; i < m.NumPoints(); ++i) {
        if (Dtype == DType::Float64)
            a.As<double>()[i] = static_cast<double>(Value);
        else
            a.As<std::int64_t>()[i] = Value;
    }
    m.AddPointData("tag", std::move(a));
    return m;
}

std::string ip_write_warnings(const std::string& rFormat, const Mesh& rMesh) {
    const std::filesystem::path dir = std::filesystem::temp_directory_path() / "mio_ip_test";
    std::filesystem::create_directories(dir);
    IpCerrCapture cap;
    try {
        registry_writers().at(rFormat)((dir / ("m." + rFormat)).string(), rMesh);
    } catch (const std::exception&) {
        // Only the warning is under test; a writer may reject the fixture.
    }
    return cap.str();
}

}  // namespace

TEST(IntegerPrecision, WarnsOnceForAnIntegerArrayAFloatFormatRounds) {
    const std::string out = ip_write_warnings("ensight", ip_mesh((1LL << 24) + 1));
    EXPECT_NE(out.find("the integer array 'tag' is stored as floating point"), std::string::npos)
        << out;
    const std::size_t first = out.find("'tag'");
    ASSERT_NE(first, std::string::npos);
    EXPECT_EQ(out.find("'tag'", first + 1), std::string::npos) << "once per array\n" << out;
}

TEST(IntegerPrecision, IsSilentForSmallIntegersFloatsAndIntegerFormats) {
    EXPECT_EQ(ip_write_warnings("ensight", ip_mesh(1LL << 24)).find("stored as floating point"),
              std::string::npos);
    EXPECT_EQ(ip_write_warnings("ensight", ip_mesh(1LL << 40, DType::Float64))
                  .find("stored as floating point"),
              std::string::npos);
    EXPECT_EQ(ip_write_warnings("vtu", ip_mesh(1LL << 40)).find("stored as floating point"),
              std::string::npos);
}

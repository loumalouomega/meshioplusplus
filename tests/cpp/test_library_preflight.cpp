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

// Library container preflight: corrupt magic and superblock versions refuse as
// ReadError before HDF5/netCDF sees them, so a library abort becomes a
// diagnosable refusal.

// System includes
#include <filesystem>
#include <fstream>
#include <string>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/exodus.hpp"
#include "meshioplusplus/formats/vtkhdf.hpp"

namespace {

std::string preflight_tmp(const std::string& rName) {
    return (std::filesystem::temp_directory_path() / rName).string();
}

void preflight_write(const std::string& rPath, const std::string& rBytes) {
    std::ofstream f(rPath, std::ios::binary);
    f.write(rBytes.data(), static_cast<std::streamsize>(rBytes.size()));
}

}  // namespace

TEST(LibraryPreflight, WrongMagicRefused) {
    const std::string path = preflight_tmp("mio-preflight-magic.vtkhdf");
    preflight_write(path, "not an hdf5 container at all........");
#ifdef MESHIOPLUSPLUS_HAS_HDF5
    EXPECT_THROW(meshioplusplus::read_vtkhdf(path), meshioplusplus::ReadError);
#else
    EXPECT_THROW(meshioplusplus::read_vtkhdf(path), meshioplusplus::ReadError);
#endif
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(LibraryPreflight, BadSuperblockRefused) {
    const std::string path = preflight_tmp("mio-preflight-superblock.vtkhdf");
    std::string bytes("\x89HDF\r\n\x1a\n", 8);
    bytes.push_back(static_cast<char>(9));  // no HDF5 release writes version 9
    bytes += std::string(32, '\0');
    preflight_write(path, bytes);
    EXPECT_THROW(meshioplusplus::read_vtkhdf(path), meshioplusplus::ReadError);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

#ifdef MESHIOPLUSPLUS_HAS_NETCDF
TEST(LibraryPreflight, ExodusWrongMagicRefused) {
    const std::string path = preflight_tmp("mio-preflight-exodus.e");
    preflight_write(path, "not netcdf nor hdf5...................");
    EXPECT_THROW(meshioplusplus::read_exodus(path), meshioplusplus::ReadError);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}
#endif

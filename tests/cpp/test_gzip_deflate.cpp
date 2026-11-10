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
#include <string>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "meshioplusplus/detail/compressed_input.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "../../src/cpp/src/detail/gzip_deflate.hpp"

namespace md = meshioplusplus::detail;

TEST(GzipDeflate, RoundTripsThroughTheReadingSide) {
    if (!md::gzip_compress_available()) {
        EXPECT_THROW(md::gzip_compress("x", "test"), meshioplusplus::WriteError);
        GTEST_SKIP() << "this build has no zlib";
    }
    std::string text;
    for (int i = 0; i < 2000; ++i)
        text += "(" + std::to_string(i) + " 0.5 1e-3)\n";
    const std::string gz = md::gzip_compress(text, "test");
    EXPECT_LT(gz.size(), text.size());
    EXPECT_EQ(md::decompress_bytes(md::Compression::Gzip, gz), text);
    EXPECT_EQ(md::decompress_bytes(md::Compression::Gzip, md::gzip_compress("", "test")), "");
}

TEST(GzipDeflate, TheHeaderIsCPythonsWithMtimeZero) {
    if (!md::gzip_compress_available())
        GTEST_SKIP() << "this build has no zlib";
    // gzip.compress(data, mtime=0): magic, deflate, no flags, mtime 0, XFL 2
    // (level 9), OS 255 -- so a native writer and its Python twin agree.
    const std::string gz = md::gzip_compress("abc", "test");
    ASSERT_GE(gz.size(), 18u);
    EXPECT_EQ(gz.substr(0, 10), std::string("\x1f\x8b\x08\x00\x00\x00\x00\x00\x02\xff", 10));
    // The trailer is CRC-32 then the input size, little-endian.
    EXPECT_EQ(gz.substr(gz.size() - 4), std::string("\x03\x00\x00\x00", 4));
}

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

// External includes
#include <gtest/gtest.h>

// System includes
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <system_error>

// Project includes
#include "meshioplusplus/detail/vtu_binary.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/vtp.hpp"
#include "meshioplusplus/formats/vtu.hpp"
#include "mesh_fixtures.hpp"

using namespace meshioplusplus;
using meshioplusplus::detail::VtkCodec;

// The optional zstd/lz4 codecs. Everything here is guarded so the suite is
// meaningful in a pure build too: the availability and error-message tests run
// unconditionally, and only the round-trips are #ifdef'd.

namespace {

void codec_remove(const std::string& rPath) {
    std::error_code ec;
    std::filesystem::remove(rPath, ec);
}

std::string codec_slurp(const std::string& rPath) {
    std::ifstream in(rPath, std::ios::binary);
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

}  // namespace

TEST(Codecs, ZlibIsAlwaysTheDefaultAndIsAvailable) {
    // The non-negotiable baseline: a build may lack zstd/lz4, but the default
    // write path must keep working.
    EXPECT_TRUE(detail::vtk_codec_available(VtkCodec::None));
#ifdef MESHIOPLUSPLUS_HAS_ZLIB
    EXPECT_TRUE(detail::vtk_codec_available(VtkCodec::Zlib));
#endif
}

TEST(Codecs, CompressorAttributeNames) {
    // LZ4's name is a real VTK class; ZSTD's is our documented extension.
    EXPECT_STREQ(detail::vtk_codec_compressor(VtkCodec::Zlib), "vtkZLibDataCompressor");
    EXPECT_STREQ(detail::vtk_codec_compressor(VtkCodec::LZ4), "vtkLZ4DataCompressor");
    EXPECT_STREQ(detail::vtk_codec_compressor(VtkCodec::ZSTD), "vtkZSTDDataCompressor");
    EXPECT_STREQ(detail::vtk_codec_compressor(VtkCodec::None), "");
}

TEST(Codecs, MissingCodecErrorNamesTheBuildOption) {
    // The whole point of the message: tell the user what to turn on.
    const std::string zstd_msg = detail::vtk_codec_missing_message(VtkCodec::ZSTD, false);
    EXPECT_NE(zstd_msg.find("MESHIOPLUSPLUS_WITH_ZSTD=ON"), std::string::npos) << zstd_msg;
    const std::string lz4_msg = detail::vtk_codec_missing_message(VtkCodec::LZ4, true);
    EXPECT_NE(lz4_msg.find("MESHIOPLUSPLUS_WITH_LZ4=ON"), std::string::npos) << lz4_msg;
    // lzma is recognized but never implemented, so it must not advertise a flag.
    const std::string lzma_msg = detail::vtk_codec_missing_message(VtkCodec::LZMA, false);
    EXPECT_NE(lzma_msg.find("not implemented"), std::string::npos) << lzma_msg;
}

TEST(Codecs, UnavailableCodecThrowsRatherThanFailingToLink) {
    // A codec compiled out must surface as a catchable error -- that is what
    // routes callers to the Python fallback instead of breaking the build.
    for (VtkCodec codec : {VtkCodec::Zlib, VtkCodec::LZ4, VtkCodec::ZSTD}) {
        if (detail::vtk_codec_available(codec))
            continue;
        EXPECT_THROW(detail::vtk_codec_require_read(codec), ReadError);
        EXPECT_THROW(detail::vtk_codec_require_write(codec), WriteError);
    }
}

// A pure build must keep reading what it always wrote. This is the guarantee
// that makes the codecs safe to add at all.
TEST(Codecs, ZlibFilesRoundTripRegardlessOfOptionalCodecs) {
#ifdef MESHIOPLUSPLUS_HAS_ZLIB
    const Mesh source = mt::data_mesh();
    const std::string path = mt::temp_path(".vtu");
    write_vtu(path, source, /*binary=*/true, /*zlib=*/true);

    const std::string text = codec_slurp(path);
    EXPECT_NE(text.find("vtkZLibDataCompressor"), std::string::npos);

    mt::expect_same_geometry(read_vtu(path), source);
    codec_remove(path);
#else
    GTEST_SKIP() << "build has no zlib";
#endif
}

namespace {

/** @brief Write with @p codec, read back, and compare against the source. */
void codec_expect_roundtrip(VtkCodec codec, const char* pCompressor) {
    const Mesh source = mt::data_mesh();
    const std::string path = mt::temp_path(".vtu");
    write_vtu_codec(path, source, /*binary=*/true, codec);

    const std::string text = codec_slurp(path);
    EXPECT_NE(text.find(pCompressor), std::string::npos) << "compressor attribute not recorded";

    const Mesh back = read_vtu(path);
    mt::expect_same_geometry(back, source);
    EXPECT_EQ(back.PointDataNames(), source.PointDataNames());
    EXPECT_EQ(back.CellDataNames(), source.CellDataNames());

    codec_remove(path);
}

}  // namespace

TEST(Codecs, ZstdRoundTrip) {
#ifdef MESHIOPLUSPLUS_HAS_ZSTD
    codec_expect_roundtrip(VtkCodec::ZSTD, "vtkZSTDDataCompressor");
#else
    GTEST_SKIP() << "build has no zstd";
#endif
}

TEST(Codecs, Lz4RoundTrip) {
#ifdef MESHIOPLUSPLUS_HAS_LZ4
    codec_expect_roundtrip(VtkCodec::LZ4, "vtkLZ4DataCompressor");
#else
    GTEST_SKIP() << "build has no lz4";
#endif
}

TEST(Codecs, VtpRoundTripsWithOptionalCodecs) {
    const Mesh source = mt::tri_mesh();
#ifdef MESHIOPLUSPLUS_HAS_LZ4
    {
        const std::string path = mt::temp_path(".vtp");
        write_vtp_codec(path, source, /*binary=*/true, VtkCodec::LZ4);
        mt::expect_same_geometry(read_vtp(path), source);
        codec_remove(path);
    }
#endif
#ifdef MESHIOPLUSPLUS_HAS_ZSTD
    {
        const std::string path = mt::temp_path(".vtp");
        write_vtp_codec(path, source, /*binary=*/true, VtkCodec::ZSTD);
        mt::expect_same_geometry(read_vtp(path), source);
        codec_remove(path);
    }
#endif
    SUCCEED();
}

// Compression must change only the framing. Every codec's payload has to match
// the uncompressed variant's exactly, which is what proves the block framing is
// shared rather than reimplemented per codec.
TEST(Codecs, EveryCodecDecodesToTheSameBytes) {
    const std::string data(100000, '\0');
    std::string payload = data;
    for (std::size_t i = 0; i < payload.size(); ++i)
        payload[i] = static_cast<char>(i * 31 + (i >> 3));
    const unsigned char* bytes = reinterpret_cast<const unsigned char*>(payload.data());

    const std::string raw = detail::vtu_encode_binary(bytes, payload.size(), VtkCodec::None);
    const std::vector<unsigned char> expected =
        detail::vtu_decode_uncompressed(raw.c_str(), raw.size(), 4);
    ASSERT_EQ(expected.size(), payload.size());

    for (VtkCodec codec : {VtkCodec::Zlib, VtkCodec::LZ4, VtkCodec::ZSTD}) {
        if (!detail::vtk_codec_available(codec))
            continue;
        const std::string encoded = detail::vtu_encode_binary(bytes, payload.size(), codec);
        const std::vector<unsigned char> decoded =
            detail::vtu_decode_blocks(encoded.c_str(), encoded.size(), 4, codec);
        EXPECT_EQ(decoded, expected)
            << "codec " << detail::vtk_codec_name(codec) << " changed the payload";
        // Multi-block: 100000 bytes exceeds the fixed 32 KiB block size.
        EXPECT_LT(encoded.size(), raw.size())
            << "codec " << detail::vtk_codec_name(codec) << " did not compress";
    }
}

TEST(Codecs, ReadingAFileNeedingAnAbsentCodecReportsTheOption) {
    // Hand-built VTU declaring a codec; if this build lacks it, the reader must
    // say which option to enable rather than failing obscurely.
    for (VtkCodec codec : {VtkCodec::LZ4, VtkCodec::ZSTD}) {
        if (detail::vtk_codec_available(codec))
            continue;
        const std::string path = mt::temp_path(".vtu");
        {
            std::ofstream os(path);
            os << "<VTKFile type=\"UnstructuredGrid\" compressor=\""
               << detail::vtk_codec_compressor(codec) << "\"><UnstructuredGrid>"
               << "<Piece NumberOfPoints=\"0\"/></UnstructuredGrid></VTKFile>";
        }
        try {
            read_vtu(path);
            ADD_FAILURE() << "expected a ReadError for " << detail::vtk_codec_name(codec);
        } catch (const ReadError& e) {
            EXPECT_NE(std::string(e.what()).find("MESHIOPLUSPLUS_WITH_"), std::string::npos)
                << e.what();
        }
        codec_remove(path);
    }
    SUCCEED();
}

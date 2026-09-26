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
/**
 * @file test_vtu_decode.cpp
 * @brief The private in-place VTU binary decode (`src/cpp/src/detail/
 *        vtu_decode.hpp`) against the buffered one it replaces.
 */

#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "../../src/cpp/src/detail/vtu_decode.hpp"
#include "meshioplusplus/detail/vtu_binary.hpp"
#include "meshioplusplus/exceptions.hpp"

namespace {

namespace det = meshioplusplus::detail;
using meshioplusplus::DType;
using meshioplusplus::NDArray;

std::vector<unsigned char> vd_bytes(std::size_t n) {
    std::vector<unsigned char> b(n);
    std::uint32_t s = 7;
    for (auto& c : b) {
        s = s * 1664525u + 1013904223u;
        c = static_cast<unsigned char>(s >> 24);
    }
    return b;
}

// The pre-v16.19 result: decode into a buffer, then copy into a zeroed array.
NDArray vd_reference(const std::string& rText, std::size_t hsz, det::VtkCodec codec, DType dt) {
    const std::vector<unsigned char> bytes =
        codec == det::VtkCodec::None
            ? det::vtu_decode_uncompressed(rText.data(), rText.size(), hsz)
            : det::vtu_decode_blocks(rText.data(), rText.size(), hsz, codec);
    const std::size_t isz = meshioplusplus::dtype_size(dt);
    NDArray a(dt, {bytes.size() / isz});
    if (a.Size())
        std::memcpy(a.Data(), bytes.data(), a.Nbytes());
    return a;
}

void vd_expect_same(const NDArray& rA, const NDArray& rB) {
    ASSERT_EQ(rA.Shape(), rB.Shape());
    ASSERT_EQ(rA.Dtype(), rB.Dtype());
    if (rA.Nbytes() > 0)  // an empty array's Data() may be null
        EXPECT_EQ(0, std::memcmp(rA.Data(), rB.Data(), rA.Nbytes()));
}

}  // namespace

TEST(VtuDecode, Base64WindowsMatchTheWholeDecode) {
    const std::vector<unsigned char> raw = vd_bytes(3 * 1000 + 2);
    std::string text = det::b64encode(raw.data(), raw.size());
    for (std::size_t i = 76; i < text.size(); i += 77)  // line-wrapped, as some writers do
        text.insert(i, "\n");
    const std::vector<unsigned char> whole = det::b64decode(text.data(), text.size());
    ASSERT_EQ(whole, raw);
    const det::VtubB64 stream(text.data(), text.size());
    ASSERT_EQ(stream.Size(), raw.size());
    for (std::size_t lo : {std::size_t{0}, std::size_t{1}, std::size_t{4}, std::size_t{8},
                           std::size_t{1499}, std::size_t{3000}})
        for (std::size_t hi : {lo, lo + 1, lo + 5, raw.size(), raw.size() + 10}) {
            std::vector<unsigned char> out(raw.size() + 16, 0xAB);
            stream.Decode(lo, hi, out.data());
            const std::size_t end = std::min(hi, raw.size());
            for (std::size_t k = lo; k < end; ++k)
                ASSERT_EQ(out[k - lo], raw[k]) << lo << " " << hi << " " << k;
            if (end >= lo)
                EXPECT_EQ(out[end - lo], 0xAB);  // nothing written past the window
        }
}

TEST(VtuDecode, InPlaceMatchesTheBufferedDecode) {
    std::vector<det::VtkCodec> codecs = {det::VtkCodec::None};
    for (det::VtkCodec c : {det::VtkCodec::Zlib, det::VtkCodec::LZ4, det::VtkCodec::ZSTD})
        if (det::vtk_codec_available(c))
            codecs.push_back(c);
    for (det::VtkCodec codec : codecs)
        for (std::size_t hsz : {std::size_t{4}, std::size_t{8}})
            for (std::size_t nbytes : {std::size_t{0}, std::size_t{8}, std::size_t{8 * 50001},
                                       std::size_t{8 * 50001 + 3}}) {
                const std::vector<unsigned char> raw = vd_bytes(nbytes);
                const std::string text = det::vtu_encode_binary(raw.data(), raw.size(), codec, hsz);
                vd_expect_same(
                    det::vtu_decode_ndarray(text.data(), text.size(), hsz, codec, DType::Float64),
                    vd_reference(text, hsz, codec, DType::Float64));
            }
}

TEST(VtuDecode, InPlaceRefusesWhatTheBufferedDecodeRefuses) {
    const std::vector<unsigned char> raw = vd_bytes(800);
    std::string text = det::vtu_encode_binary(raw.data(), raw.size(), det::VtkCodec::None, 4);
    const std::string truncated = text.substr(0, text.size() / 2);
    EXPECT_THROW(det::vtu_decode_ndarray(truncated.data(), truncated.size(), 4, det::VtkCodec::None,
                                         DType::Float64),
                 meshioplusplus::ReadError);
    EXPECT_THROW(det::vtu_decode_ndarray("QQ==", 4, 4, det::VtkCodec::None, DType::Float64),
                 meshioplusplus::ReadError);
}

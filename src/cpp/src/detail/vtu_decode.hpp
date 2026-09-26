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

#pragma once

/**
 * @file detail/vtu_decode.hpp
 * @brief A VTK XML "binary" DataArray decoded straight into its array.
 *
 * A **core-private** header (the `slot_runs.hpp` precedent): no installed
 * header names it, and it adds nothing to the API or the ABI; the exported
 * `vtu_decode_uncompressed`/`vtu_decode_blocks`/`vtu_parse_binary` keep their
 * signatures and results.
 *
 * The readers stripped the element text into a copy, decoded its base64 into
 * one buffer (and, uncompressed, copied the payload out of it), or
 * decompressed every block into a buffer of its own and copied it on, and
 * then zero-filled an array and copied the bytes in once more. Here the text
 * is a view of the document's own characters, the base64 stream decodes any
 * byte window straight to a destination (`VtubB64`), and every compressed
 * block decompresses into its place in an uninitialised array. Roadmap §3,
 * "The VTU binary read copies each payload five or six times".
 */

// System includes
#include <cctype>
#include <cstddef>
#include <string_view>
#include <vector>

// Project includes
#include "meshioplusplus/detail/vtu_binary.hpp"
#include "meshioplusplus/ndarray.hpp"

namespace meshioplusplus {
namespace detail {

/**
 * @brief A base64 text, its valid characters counted once per fixed chunk:
 * `Decode` then writes any byte window of the decoded stream straight to a
 * destination, the chunks that overlap it in parallel. Characters outside the
 * alphabet (padding, whitespace, line breaks) are skipped, as `b64decode`
 * skips them.
 */
class VtubB64 {
public:
    VtubB64(const char* pS, std::size_t len);

    /// Decoded bytes in the whole stream.
    std::size_t Size() const { return mSize; }

    /// Writes the decoded bytes [Lo, Hi) (clamped to `Size()`) to pOut[0 ..).
    void Decode(std::size_t Lo, std::size_t Hi, unsigned char* pOut) const;

private:
    const char* mpS;
    std::size_t mLen;
    std::vector<std::size_t> mValid;  ///< valid characters before each chunk
    std::size_t mSize = 0;
};

/**
 * @brief `vtu_parse_binary`'s result -- the DataArray's payload as `dt`
 * elements, a trailing partial element dropped -- decoded in place, with the
 * checks and messages of `vtu_decode_uncompressed`/`vtu_decode_blocks`.
 */
NDArray vtu_decode_ndarray(const char* pText, std::size_t len, std::size_t hsz, VtkCodec codec,
                           DType dt);

/**
 * @brief Appends to @p rOut what `vtu_encode_binary` base64-encodes, as raw
 * bytes: the `hsz`-byte size header and the payload, or with @p codec the
 * block header and the compressed blocks -- one array of a raw
 * `<AppendedData>` section.
 */
void vtu_encode_raw(const unsigned char* pData, std::size_t nbytes, VtkCodec codec, std::size_t hsz,
                    std::vector<unsigned char>& rOut);

/// `vtu_parse_binary(vtu_strip(...), ...)` over a view of the element's text.
inline NDArray vtu_decode_bin_view(std::string_view Text, DType dt, VtkCodec codec,
                                   std::size_t hsz) {
    return vtu_decode_ndarray(Text.data(), Text.size(), hsz, codec, dt);
}

/// `vtu_strip` as a view of @p pS itself, with no copy.
inline std::string_view vtu_strip_view(const char* pS) {
    std::string_view t = pS ? std::string_view(pS) : std::string_view();
    std::size_t b = 0, e = t.size();
    while (b < e && std::isspace(static_cast<unsigned char>(t[b])))
        ++b;
    while (e > b && std::isspace(static_cast<unsigned char>(t[e - 1])))
        --e;
    return t.substr(b, e - b);
}

}  // namespace detail
}  // namespace meshioplusplus

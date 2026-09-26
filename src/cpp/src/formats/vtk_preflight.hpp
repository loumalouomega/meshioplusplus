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
 * @file vtk_preflight.hpp
 * @brief Refuse a VTK XML file the C++ readers cannot decode from its start
 *        tag, before its DOM is loaded.
 *
 * A **core-private** header beside the formats (the `face_cells_common.hpp`
 * precedent): it adds nothing to the installed headers or the ABI.
 *
 * The five VTK XML readers (`.vtu`, `.vtp`, `.vti`, `.vtr`, `.vts`) refused an
 * lzma compressor, or an lz4 or zstd one the build lacks, only after pugixml
 * had loaded the whole document, base64 bodies included -- and the caller then
 * read the file again through the Python reader. The compressor is an
 * attribute of the `<VTKFile ...>` start tag, so the first few kilobytes
 * decide it. The pre-flight is loose on purpose: whenever it cannot find the
 * tag, or the tag names another file type, it decides nothing and the full
 * parse gives its usual answer. Roadmap §3, "A declined C++ read".
 */

// System includes
#include <cstddef>
#include <string>
#include <string_view>

// Project includes
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/vtu_binary.hpp"
#include "meshioplusplus/exceptions.hpp"

namespace meshioplusplus {
namespace detail {

/// The value of attribute @p Name in the start tag @p Tag, or "" without one.
inline std::string vtk_preflight_attribute(std::string_view Tag, std::string_view Name) {
    std::size_t at = 0;
    while ((at = Tag.find(Name, at)) != std::string_view::npos) {
        const bool starts = at > 0 && (Tag[at - 1] == ' ' || Tag[at - 1] == '\t' ||
                                       Tag[at - 1] == '\n' || Tag[at - 1] == '\r');
        std::size_t p = at + Name.size();
        at = p;
        if (!starts)
            continue;
        while (p < Tag.size() &&
               (Tag[p] == ' ' || Tag[p] == '\t' || Tag[p] == '\n' || Tag[p] == '\r'))
            ++p;
        if (p >= Tag.size() || Tag[p] != '=')
            continue;
        ++p;
        while (p < Tag.size() &&
               (Tag[p] == ' ' || Tag[p] == '\t' || Tag[p] == '\n' || Tag[p] == '\r'))
            ++p;
        if (p >= Tag.size() || (Tag[p] != '"' && Tag[p] != '\''))
            continue;
        const char quote = Tag[p++];
        const std::size_t close = Tag.find(quote, p);
        if (close == std::string_view::npos)
            return {};
        return std::string(Tag.substr(p, close - p));
    }
    return {};
}

/**
 * @brief Throws, before the document is parsed, the error the reader's header
 * check would give for the file's compressor: @p pLzmaMessage for lzma, and
 * `vtk_codec_require_read`'s for an lz4 or zstd codec this build lacks.
 * Decides nothing when the `<VTKFile>` start tag is not within the first
 * 4 KiB or its `type` is not @p pType.
 */
inline void vtk_preflight(const std::string& rPath, const char* pType, const char* pLzmaMessage) {
    auto in = make_classic_ifstream(rPath, std::ios::binary);
    if (!in)
        return;  // the full parse reports an unreadable file
    std::string head(4096, '\0');
    in.read(head.data(), static_cast<std::streamsize>(head.size()));
    head.resize(static_cast<std::size_t>(in.gcount()));
    const std::size_t tag = head.find("<VTKFile");
    if (tag == std::string::npos)
        return;
    const std::size_t end = head.find('>', tag);
    if (end == std::string::npos)
        return;
    const std::string_view start_tag(head.data() + tag, end - tag);
    if (vtk_preflight_attribute(start_tag, "type") != pType)
        return;
    const std::string compressor = vtk_preflight_attribute(start_tag, "compressor");
    if (compressor == vtk_codec_compressor(VtkCodec::LZMA))
        throw ReadError(pLzmaMessage);
    if (compressor == vtk_codec_compressor(VtkCodec::LZ4))
        vtk_codec_require_read(VtkCodec::LZ4);
    else if (compressor == vtk_codec_compressor(VtkCodec::ZSTD))
        vtk_codec_require_read(VtkCodec::ZSTD);
}

}  // namespace detail
}  // namespace meshioplusplus

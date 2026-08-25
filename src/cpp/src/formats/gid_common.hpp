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
 * @file formats/gid_common.hpp
 * @brief Path/flavour helpers shared by the GiD writer and reader.
 *
 * A **format-private** header (the `formats/xdmf_doc.hpp` precedent, included
 * by both `xdmf.cpp` and `xdmf_time_series.cpp`): it lives beside the `.cpp`
 * files rather than under `src/cpp/include/`, because nothing outside this
 * format needs it and no installed header may name it. The amalgamator
 * resolves a quoted include against the including file's own directory first,
 * so it is inlined exactly once with no generator change.
 *
 * These three helpers used to sit in `gid.cpp`'s anonymous namespace **inside**
 * its `#ifdef MESHIOPLUSPLUS_HAS_GIDPOST` guard. They moved here for two
 * reasons, and both matter:
 *
 *  1. The **reader needs them but must not be behind that guard.** gidpost is
 *     a write-only library; reading a GiD file is ordinary text/record parsing
 *     whose real dependencies are per-flavour (ASCII: none, binary: zlib,
 *     HDF5: HDF5), not "was gidpost compiled in". Leaving them behind the
 *     guard would have forced the reader behind it too, and a zlib-less build
 *     (the statically-linked release CLI binaries, the Windows wheels) would
 *     then carry no GiD reader for no reason at all.
 *  2. They are in a **named** namespace, not an anonymous one, because the
 *     amalgamation concatenates every `src/cpp/src/**.cpp` into a single
 *     translation unit -- two files each carrying an anonymous-namespace
 *     `gid_has_suffix` would be a redefinition there even though they compile
 *     fine separately.
 */

// System includes
#include <string>
#include <utility>

// Project includes
#include "meshioplusplus/formats/gid.hpp"

namespace meshioplusplus {
namespace gid_detail {

/// True when @p rS ends with @p rSuffix.
inline bool gid_has_suffix(const std::string& rS, const std::string& rSuffix) {
    return rS.size() >= rSuffix.size() &&
           rS.compare(rS.size() - rSuffix.size(), rSuffix.size(), rSuffix) == 0;
}

/**
 * @brief Resolves `GidMode::Auto` against a path's extension.
 *
 * `.post.bin` -> Binary, `.post.h5` -> Hdf5, anything else (including
 * `.post.msh`/`.post.res`) -> Ascii. A non-`Auto` mode passes through.
 *
 * @note On the READ side this is only the first guess: `read_gid` additionally
 * sniffs the leading bytes, because a `.post.msh` may legitimately be gzipped
 * (gidpost's `GiD_PostAsciiZipped` writes the same ASCII text through
 * `gzprintf`), and the extension alone cannot say so.
 */
inline GidMode gid_resolve_mode(const std::string& rPath, GidMode mode) {
    if (mode != GidMode::Auto)
        return mode;
    if (gid_has_suffix(rPath, ".post.bin"))
        return GidMode::Binary;
    if (gid_has_suffix(rPath, ".post.h5"))
        return GidMode::Hdf5;
    return GidMode::Ascii;
}

/**
 * @brief `"<stem>.post.msh"` <-> `"<stem>.post.res"`.
 *
 * An arbitrary path is treated as the stem itself (the
 * `ensight_case_geo_paths` precedent, generalized to a third case since gid's
 * two known suffixes share the same length).
 *
 * @return `{mesh path, result path}`.
 */
inline std::pair<std::string, std::string> gid_ascii_paths(const std::string& rPath) {
    if (gid_has_suffix(rPath, ".post.msh"))
        return {rPath, rPath.substr(0, rPath.size() - 3) + "res"};
    if (gid_has_suffix(rPath, ".post.res"))
        return {rPath.substr(0, rPath.size() - 3) + "msh", rPath};
    return {rPath + ".post.msh", rPath + ".post.res"};
}

}  // namespace gid_detail
}  // namespace meshioplusplus

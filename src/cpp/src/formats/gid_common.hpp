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
#include <array>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

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

/**
 * @brief The one `{GidResultType, GiD spelling, legal component counts}`
 * table, shared by the writer and the reader.
 *
 * Counts are gidpost's own `_ResultTypeInfo` (`gidpostInt.c`); spellings are
 * what `GetResultTypeName` emits and what a `.post.res` `Result` header
 * therefore carries. A single table because the reader validates counts it
 * parsed and the writer validates counts it is about to emit -- two
 * transcriptions of an irregular nine-row table would drift, and a wrong row
 * produces a plausible-looking file rather than an error.
 *
 * A zero terminates the count list (gidpost's own convention).
 */
struct GidResultTypeEntry {
    GidResultType mType;
    const char* mName;
    std::array<std::size_t, 3> mDims;  // zero-terminated
};

inline const std::vector<GidResultTypeEntry>& gid_result_type_table() {
    static const std::vector<GidResultTypeEntry> table = {
        {GidResultType::Scalar, "Scalar", {1, 0, 0}},
        {GidResultType::Vector, "Vector", {2, 3, 4}},
        {GidResultType::Matrix, "Matrix", {3, 6, 0}},
        {GidResultType::PlainDeformationMatrix, "PlainDeformationMatrix", {4, 0, 0}},
        {GidResultType::MainMatrix, "MainMatrix", {12, 0, 0}},
        {GidResultType::LocalAxes, "LocalAxes", {3, 0, 0}},
        {GidResultType::ComplexScalar, "ComplexScalar", {2, 0, 0}},
        {GidResultType::ComplexVector, "ComplexVector", {4, 6, 0}},
        {GidResultType::ComplexMatrix, "ComplexMatrix", {6, 12, 0}},
    };
    return table;
}

/// The entry for @p type, or nullptr when @p type is out of range (which a
/// caller-supplied `field_data` value can be).
inline const GidResultTypeEntry* gid_find_result_type(GidResultType type) {
    for (const GidResultTypeEntry& e : gid_result_type_table())
        if (e.mType == type)
            return &e;
    return nullptr;
}

/// Human-readable legal counts for @p type, for an error message.
inline std::string gid_legal_dims_text(const GidResultTypeEntry& rEntry) {
    std::string out;
    for (std::size_t d : rEntry.mDims) {
        if (d == 0)
            break;
        if (!out.empty())
            out += ", ";
        out += std::to_string(d);
    }
    return out;
}

/**
 * @brief The type the writer picks for a @p k -component array with no
 * declaration, or nullptr when it would split the array into k scalars.
 *
 * Shared so the reader can decide whether a declaration carries information:
 * it records one only when the file's declared type differs from what this
 * would have inferred, which is what keeps an ordinary scalar/3-vector round
 * trip free of `gid:result_type:*` noise (the all-zero-material-column rule).
 */
inline const GidResultTypeEntry* gid_inferred_result_type(std::size_t k) {
    if (k == 1)
        return gid_find_result_type(GidResultType::Scalar);
    if (k == 2 || k == 3)
        return gid_find_result_type(GidResultType::Vector);
    return nullptr;
}

}  // namespace gid_detail
}  // namespace meshioplusplus

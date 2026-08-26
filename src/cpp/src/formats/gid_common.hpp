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

/**
 * @brief Node-slot permutations for the three quadratic cell types whose GiD
 * ordering was, until now, "not independently verified" -- `hexahedron27`,
 * `wedge15`, `pyramid13`.
 *
 * Shared between the writer and the reader for the same reason the result
 * type table above is: a permutation transcribed twice risks two silently
 * different tables, and a wrong entry here does not fail loudly -- it writes
 * or reads a plausible file with quietly transposed quadratic nodes.
 *
 * ## Derivation
 *
 * Every entry meshio++ already supported was cross-checked against Kratos
 * Multiphysics's production GiD writer; these three were refused because no
 * such cross-check existed. Read from Kratos's own geometry classes
 * (`kratos/geometries/hexahedra_3d_27.h`, `prism_3d_15.h`, `pyramid_3d_13.h`)
 * and independently confirmed against Kratos's Kratos-to-VTK conversion
 * utility (`kratos/input_output/vtk_output.cpp`, mirrored in
 * `ensight_output.cpp`) -- an Element-agnostic source, unlike the
 * Hexahedra3D20 reorder in `gid_mesh_container.h`, which the research
 * confirming this table found lives only in that file's *Conditions*-writing
 * path (elements are written with no reorder at all). That does not change
 * `hexahedron20`'s own identity mapping below -- `vtk_output.cpp`'s
 * conversion, which every Hexahedra3D20 element or condition goes through for
 * VTK/EnSight output, independently reproduces the identical swap -- but it
 * is the more precise citation and the one these three lean on, since two of
 * them (`hexahedron27`, `wedge15`) have no GiD-specific Kratos precedent at
 * all beyond "no reorder is applied", and the third (`pyramid13`) has no
 * Kratos GiD precedent whatsoever -- Kratos never registers a GiD mesh
 * container for any pyramid, of any order, so its ordering rests on Kratos's
 * internal geometry convention alone (confirmed to already equal VTK's own
 * convention by `vtk_output.cpp` explicitly skipping it: `Pyramid3D13` falls
 * through that function's `else` branch, needing no conversion).
 *
 * `hexahedron27`: meshio++'s own table (`cell_subdivision.cpp`) orders edges
 * 8-11 bottom ring, 12-15 TOP ring, 16-19 verticals; Kratos's internal order
 * is 8-11 bottom ring, 12-15 VERTICALS, 16-19 top ring -- the reverse split.
 * Face centres 20-25 also disagree: meshio++ orders them x-min/x-max/y-min/
 * y-max/bottom/top, Kratos orders them bottom/y-min/x-max/y-max/x-min/top.
 * Body centre 26 agrees. The permutation below is `dst[c] = src[p[c]]` (the
 * `med_node_perm()`/`flatten_f`/`unflatten_f` convention already used
 * elsewhere in this repo, reused rather than reinvented) mapping a
 * GiD/Kratos slot to the meshio++ index holding the same edge/face/body; it
 * is self-inverse (an involution), like every other permutation this
 * convention has produced so far, and matches `vtk_output.cpp`'s own
 * Kratos-index array verbatim.
 *
 * `wedge15`: meshio++ orders edges 6-8 bottom triangle, 9-11 TOP triangle,
 * 12-14 verticals; Kratos orders 6-8 bottom triangle, 9-11 VERTICALS, 12-14
 * top triangle -- the same reverse-split pattern, and again self-inverse.
 *
 * `pyramid13`: identical in both conventions (base-ring edges 5-8, apex
 * edges 9-12) -- no permutation, `mPerm == nullptr`.
 */
struct GidCellPermEntry {
    const char* mMeshioName;
    std::size_t mNumNodes;
    const int* mPerm;  // nullptr = identity
};

/// `dst[c] = src[p[c]]`, self-inverse. See `hexahedron27`'s derivation above;
/// independently confirmed against Kratos's `vtk_output.cpp` array verbatim.
/// Indexed and grouped deliberately, not a flat literal, so a slipped digit
/// here is easy to catch by eye against the derivation comment above:
///   corners 0-7, bottom-ring edges 8-11: identity
///   top-ring edges 12-15 <- meshio++'s verticals 16-19
///   verticals 16-19      <- meshio++'s top-ring edges 12-15
///   face centres: 20<-24(bottom), 21<-22(y-min), 22<-21(x-max),
///                 23<-23(y-max, fixed), 24<-20(x-min), 25<-25(top, fixed)
///   body centre 26: fixed
inline constexpr int kGidHexahedron27Perm[27] = {
    0,  1,  2,  3,  4,  5,  6, 7,  // corners
    8,  9,  10, 11,                // bottom-ring edges
    16, 17, 18, 19,                // slot 12-15 (top ring)   <- verticals
    12, 13, 14, 15,                // slot 16-19 (verticals)  <- top ring
    24, 22, 21, 23, 20, 25,        // face centres 20-25
    26,                            // body centre
};

/// `dst[c] = src[p[c]]`, self-inverse. See `wedge15`'s derivation above:
///   corners 0-5, bottom-triangle edges 6-8: identity
///   top-triangle edges 9-11 <- meshio++'s verticals 12-14
///   verticals 12-14         <- meshio++'s top-triangle edges 9-11
inline constexpr int kGidWedge15Perm[15] = {
    0,  1,  2,  3, 4, 5,  // corners
    6,  7,  8,            // bottom-triangle edges
    12, 13, 14,           // slot 9-11 (top triangle) <- verticals
    9,  10, 11,           // slot 12-14 (verticals)   <- top triangle
};

inline const std::vector<GidCellPermEntry>& gid_cell_perm_table() {
    static const std::vector<GidCellPermEntry> table = {
        {"hexahedron27", 27, kGidHexahedron27Perm},
        {"wedge15", 15, kGidWedge15Perm},
        {"pyramid13", 13, nullptr},
    };
    return table;
}

/// The permutation for @p rMeshioName / @p nnode, or nullptr for identity
/// (either because the type needs none, like `pyramid13`, or because it is
/// not in this table at all -- every OTHER GiD-supported type is identity).
inline const int* gid_cell_perm(const std::string& rMeshioName, std::size_t nnode) {
    for (const GidCellPermEntry& e : gid_cell_perm_table())
        if (e.mNumNodes == nnode && rMeshioName == e.mMeshioName)
            return e.mPerm;
    return nullptr;
}

}  // namespace gid_detail
}  // namespace meshioplusplus

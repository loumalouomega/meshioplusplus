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
 * @file detail/mesh_carve.hpp
 * @brief The shared "split a mesh into named cell-index parts" helper used by every
 * composite/multi-part writer (VTKHDF's composite types, EnSight Gold's multi-`part`
 * geometry).
 *
 * A part is a name plus a sorted list of global (block-major) cell indices. When a
 * mesh's Cell regions exactly partition every cell -- no overlap, no gap -- each
 * region becomes one part, in `(tag, name)` order, which is also the mesh's own
 * region storage order. Otherwise the mesh falls back to one part per cell block,
 * named `block_<i>`, with a `log::warn` naming the reason.
 */

#include <cstddef>
#include <string>
#include <vector>

#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {
namespace detail {

/// One named part: a sorted list of global (block-major) cell indices.
struct MeshPart {
    std::string mName;
    std::vector<std::size_t> mCells;
};

/**
 * @brief Splits @p rMesh's cells into named parts, preferring its Cell regions.
 *
 * @param rMesh The mesh to carve. Must have at least one cell (`WriteError` otherwise,
 *        naming @p rCaller).
 * @param rCaller Short name of the calling format (`"vtkhdf"`, `"ensight"`, ...), used
 *        in every exception/warning message.
 * @param StrictNames When `true`, a region-name collision (or, with @p RejectSlash, a
 *        `/` in a name) is a `WriteError` naming @p rCaller (VTKHDF's own composite
 *        blocks need this: a soft fallback would silently rename the caller's Assembly
 *        links). When `false` (the default), it is a `log::warn` and a fallback to one
 *        part per cell block, the same fallback used when the regions do not exactly
 *        partition the mesh either way.
 * @param RejectSlash When `true`, a part name containing `/` is a `WriteError` (only
 *        meaningful together with @p StrictNames; VTKHDF's Assembly links are a path
 *        hierarchy where `/` would be misread as a separator).
 * @return One part per Cell region when they exactly partition every cell (an empty
 *         region becomes an empty part), else one part per cell block.
 */
MESHIOPLUSPLUS_API std::vector<MeshPart> carve_by_region(const Mesh& rMesh,
                                                         const std::string& rCaller,
                                                         bool StrictNames = false,
                                                         bool RejectSlash = false);

}  // namespace detail
}  // namespace meshioplusplus

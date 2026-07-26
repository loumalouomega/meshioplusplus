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
 * @file detail/region_remap.hpp
 * @brief The one shared "carry a mesh's regions through an operation" helper:
 * drop removed entries, remap survivors, expand parents to children, and decide
 * whether a side-set facet still exists.
 *
 * Every operation that renumbers points or cells already returns index maps
 * (`mPointMap`, `mCellMaps`) so its caller can carry per-entity information
 * across. Before regions existed, nine operation shims each re-implemented that
 * carry in Python, subtly differently. This header is the single owner of it,
 * in the repo's shared-`detail/` idiom (`subset.hpp`, `node_adjacency.hpp`,
 * `spatial_hash.hpp`, `marching.hpp`, `cell_subdivision.hpp`,
 * `space_filling.hpp`).
 *
 * ## The three cell-map shapes
 *
 * Operations do not all describe their cell maps the same way, and pretending
 * otherwise would silently corrupt regions, so the shape is explicit:
 *
 *  - `CellMapKind::Direct` — per input block, `map[c]` is *the* output cell's
 *    index within the corresponding output block, or -1 if it was dropped.
 *    (crop, split, clean, partition, reorder.)
 *  - `CellMapKind::FirstChild` — per input block, `map[c]` is the index of the
 *    **first** of a contiguous run of children. The run ends at the next
 *    non-negative map entry, or at the end of the output block. (convert_cells
 *    under `Simplexify`, refine, decimate.) This shape also covers 1:1 maps
 *    correctly, but `Direct` must not be replaced by it: a permutation's map is
 *    not monotone, so the "next entry" rule would invent nonsense ranges.
 *  - `CellMapKind::Global` — one flat array indexed by input **global** cell
 *    index, yielding an output global cell index. (merge, whose maps are per
 *    input *mesh* rather than per block.)
 *
 * ## Side regions
 *
 * A side entry `(global cell, local facet)` survives only if its cell survives
 * **and** the facet still exists: the output cell must have the same type as
 * the input cell, and the facet index must still be in range for that type.
 * Under `FirstChild` a parent's children are new cells of a different or
 * subdivided topology, so there is no facet correspondence to preserve at all
 * and side regions are dropped by name. `mDropSideRegions` forces that for any
 * operation that knows its facets do not survive.
 *
 * Free functions in `meshioplusplus::detail`, called once per operation rather
 * than per element, so the bodies live in `src/cpp/src/detail/region_remap.cpp`.
 * Built on the uniform mesh API only, so it works under every backend.
 */

// System includes
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/region.hpp"

namespace meshioplusplus {
namespace detail {

/// How an operation's cell map describes the input -> output correspondence.
enum class CellMapKind {
    Direct,      ///< per block: map[c] = the output cell, -1 = dropped
    FirstChild,  ///< per block: map[c] = first of a contiguous child run
    Global,      ///< flat: map[global in cell] = global out cell, -1 = dropped
};

/// The index maps an operation hands `remap_regions`.
struct RegionRemap {
    /// Int64 `(num_points_in,)`, input point -> output point (-1 = dropped).
    /// Null means points were not renumbered (identity).
    const NDArray* pPointMap = nullptr;

    /// Which shape `pCellMaps` / `pGlobalCellMap` is in.
    CellMapKind mCellMapKind = CellMapKind::Direct;

    /// Per input block, Int64 `(num_cells_in_block,)`. Used by `Direct` and
    /// `FirstChild`. Null means cells were not renumbered (identity).
    const std::vector<NDArray>* pCellMaps = nullptr;

    /// Flat Int64 map used by `CellMapKind::Global`.
    const NDArray* pGlobalCellMap = nullptr;

    /// Input block -> output block. Empty means the identity, which is the case
    /// for every operation that keeps its block structure 1:1; `split` (which
    /// drops empty blocks) is the one that must fill it in. An entry of
    /// `kBlockDropped` means the block has no output counterpart.
    std::vector<std::size_t> mBlockMap;

    /// Force side regions to be dropped even under a `Direct` map.
    bool mDropSideRegions = false;

    /// Operation name, used in the "regions dropped" warning.
    std::string mOpName;
};

/// `RegionRemap::mBlockMap` entry meaning "this input block has no output block".
inline constexpr std::size_t kBlockDropped = static_cast<std::size_t>(-1);

/**
 * @brief Number of facets a cell type has: faces for a 3-D type
 * (`cell_faces.hpp`), edges for a 2-D one (`cell_edges.hpp`).
 * @param rType The meshio cell-type name.
 * @return The facet count, or 0 for a type with no facet table (which makes
 *         every side entry on it invalid, and so dropped).
 */
MESHIOPLUSPLUS_API std::size_t region_num_facets(const std::string& rType);

/**
 * @brief Carry one region across, without adding it to the output mesh.
 *
 * The building block `remap_regions` loops over. It is public because `merge`
 * needs to rename a region *between* remapping it and storing it: two inputs
 * may carry the same region name, and adding both under that name would make
 * the second replace the first (they share a `(kind, name, dim, tag)` key).
 *
 * @param rIn The operation's input mesh.
 * @param rOut The operation's output mesh (read for block bases and cell types).
 * @param rRegion The input region to carry.
 * @param rMaps The index maps.
 * @param rResult Receives the carried region on success.
 * @return `false` when nothing survived, or when the kind cannot be carried at
 *         all (in which case a `log::warn` has already been emitted).
 */
MESHIOPLUSPLUS_API bool remap_region(const Mesh& rIn, const Mesh& rOut, const Region& rRegion,
                  const RegionRemap& rMaps, Region& rResult);

/**
 * @brief Carry every region of @p rIn onto @p rOut through @p rMaps.
 *
 * Entries whose entity did not survive are dropped; a region left with no
 * entries at all is not added to the output. Regions whose kind cannot be
 * carried (see the file docs) are dropped with a single `log::warn` naming the
 * region and the operation.
 *
 * @param rIn The operation's input mesh.
 * @param rOut The operation's output mesh; regions are added to it.
 * @param rMaps The index maps.
 */
MESHIOPLUSPLUS_API void remap_regions(const Mesh& rIn, Mesh& rOut, const RegionRemap& rMaps);

/**
 * @brief Drop every region with one warning, for operations whose output has no
 * entity correspondence with their input at all (slice, isosurface, surface
 * extraction).
 *
 * A deliberate no-op when the input carries no regions, so the common case
 * stays silent.
 *
 * @param rIn The operation's input mesh.
 * @param rOpName The operation name, for the warning.
 */
MESHIOPLUSPLUS_API void warn_regions_dropped(const Mesh& rIn, const std::string& rOpName);

}  // namespace detail
}  // namespace meshioplusplus

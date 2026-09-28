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
 * A side entry `(global cell, local facet)` whose cell keeps its identity --
 * one output cell of the same type, reached by no other input cell, whose
 * facet at that number still has the same nodes -- keeps its number. Any
 * other entry (since v16.26.0) is found again by what the facet is made of:
 *
 *  - **refined** (several children): every child facet lying within it -- each
 *    of the child facet's nodes the image of one of the facet's, or a point
 *    the operation created that lies on it;
 *  - **merged or retyped** (one child that others share, or of another type):
 *    the output facet containing it -- each surviving node one of that
 *    facet's, and a node the operation removed lying on it (a second pass
 *    also accepts a surviving node that only lies on it, for a merged face
 *    that dropped an interior vertex).
 *
 * When neither finds it under a `FirstChild` map, a facet whose corners all
 * survive is looked up across the whole output by those corners (`FacetIndex`):
 * a decimation collapse hands a boundary edge to a neighbouring cell. (Not
 * under `Direct`: there the facet's own cell is gone, and handing an interface
 * facet to the neighbour would flip its orientation.) Node identity decides wherever it
 * can; geometry (a relative 1e-9 plane or segment test) only places the
 * points an operation created or removed, which are averages of the facet's
 * own corners. An entry with no counterpart is
 * dropped, with one warning per region counting them. Polygon edges and
 * polyhedron faces count as facets too (edge k of a polygon runs from node k
 * to node k + 1; face k of a polyhedron is its k-th face). `mDropSideRegions`
 * still drops every side region, by name, for an operation that asks.
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

    /// Drop every side region, by name, instead of carrying it.
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
 * @brief Carry regions onto a mesh of facets extracted from @p rIn
 * (`extract_surface`, `extract_skin`).
 *
 * A Side region becomes a Cell region of the same name, dim and tag naming the
 * output cells its facets became; a facet that was not extracted (an interior
 * one) is lost, with one warning per region counting them. A Point region
 * follows @p rPointMap. A Cell region names input cells, which are not in the
 * output, so it is dropped with a warning.
 *
 * @param rIn The input mesh.
 * @param rOut The facet mesh; regions are added to it.
 * @param rOutParent Per output cell (global, block-major), the input global
 *        cell it is a facet of.
 * @param rOutFacet Per output cell, which facet of that cell it is (the Side
 *        numbering: `cell_faces`/`cell_edges`, a polyhedron's face index).
 * @param rPointMap Input point -> output point, -1 for a point not kept.
 * @param rOpName The operation name, for the warnings.
 */
MESHIOPLUSPLUS_API void carry_regions_to_facet_mesh(const Mesh& rIn, Mesh& rOut,
                                                    const std::vector<std::int64_t>& rOutParent,
                                                    const std::vector<std::int64_t>& rOutFacet,
                                                    const std::vector<std::int64_t>& rPointMap,
                                                    const std::string& rOpName);

/**
 * @brief Drop every region with one warning, for operations whose output has no
 * entity correspondence with their input at all (slice, isosurface).
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

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
 * @file operations/region_ops.hpp
 * @brief Set algebra and bookkeeping on named regions: union, intersection,
 * difference, rename, retag and delete.
 *
 * Every consumer that builds boundary conditions from regions ends up
 * combining them -- "the inlet is these two gmsh groups", "the wall is the
 * skin minus the inlet and outlet", "call it `Inlet` and give it tag 10" --
 * and before v16.23.0 each one hand-rolled it. `edit_regions` applies a list of
 * such edits in order to a copy of the mesh. It is a data-only operation:
 * points, cells, data and property sets are untouched.
 *
 * ### Selecting a region
 *
 * A region's identity is `(kind, name, dim, tag)` (`region.hpp`), and a file
 * may carry the same name twice -- gmsh allows a physical name per dimension.
 * A `RegionSelector` names the region and optionally pins the kind, dimension
 * and tag; it must match **exactly one** region, and an ambiguous or missing
 * match throws with the candidates listed, rather than picking the first.
 *
 * ### The edits
 *
 *  - **Union / Intersection / Difference** (two or more inputs, all of one
 *    kind; difference is the first minus the rest) add a region named
 *    `mOutputName`. Side regions combine their (cell, facet) pairs. The result
 *    is kept even when empty -- a named group is information. Its dimension is
 *    the inputs' shared one (-1 if they differ) and its tag -1, unless
 *    `mOutputDim` / `mOutputTag` set them. With `mKeepInputs` off the inputs are
 *    removed.
 *  - **Rename** (one input) moves the region to `mOutputName`, keeping its
 *    kind, entries, dimension and tag unless `mOutputDim` / `mOutputTag` are set.
 *  - **Retag** (one input) sets its tag and/or dimension; the name is kept
 *    unless `mOutputName` is given. A gmsh write takes a cell's physical tag
 *    from the `gmsh:physical` cell data when present, so retagging a Cell
 *    region of such a mesh warns that the array still carries the old tag.
 *  - **Delete** (one or more inputs) removes them.
 *
 * An edit whose output key collides with an existing region other than one of
 * its own inputs throws: an edit never silently overwrites an unrelated region.
 */

// System includes
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// The edits `edit_regions` knows.
enum class RegionOp : std::int32_t {
    Union = 0,
    Intersection = 1,
    Difference = 2,
    Rename = 3,
    Retag = 4,
    Delete = 5,
};

/// The name of @p Op as the CLIs, pipelines and bindings spell it:
/// `union`, `intersection`, `difference`, `rename`, `retag`, `delete`.
MESHIOPLUSPLUS_API const char* region_op_name(RegionOp Op);

/// The `RegionOp` named @p rName (also accepting `intersect`).
/// @throws std::invalid_argument on an unknown name.
MESHIOPLUSPLUS_API RegionOp region_op_from_name(const std::string& rName);

/// Sentinel for "any" in a `RegionSelector`, and "inherit" in a `RegionEdit`.
inline constexpr std::int64_t kRegionAny = -2;

/// Names one region of a mesh.
struct RegionSelector {
    std::string mName;
    /// A `RegionKind` value (0 point, 1 cell, 2 side), or -1 for any kind.
    std::int32_t mKind = -1;
    /// The dimension to match, or `kRegionAny`.
    std::int64_t mDim = kRegionAny;
    /// The tag to match, or `kRegionAny`.
    std::int64_t mTag = kRegionAny;
};

/// One edit; see the file comment for what each operation does.
struct RegionEdit {
    RegionOp mOp = RegionOp::Union;
    std::vector<RegionSelector> mInputs;
    /// The result's name (union/intersection/difference/rename; optional for
    /// retag).
    std::string mOutputName;
    /// The result's dimension, or `kRegionAny` to inherit.
    std::int64_t mOutputDim = kRegionAny;
    /// The result's tag, or `kRegionAny` to inherit (-1 for a set operation).
    std::int64_t mOutputTag = kRegionAny;
    /// Set operations only: keep the input regions.
    bool mKeepInputs = true;
};

/**
 * @brief The index (`Region(i)` numbering) of the one region @p rSelector matches.
 * @throws std::invalid_argument when none or several match.
 */
MESHIOPLUSPLUS_API std::size_t find_region(const Mesh& rMesh, const RegionSelector& rSelector);

/**
 * @brief Apply @p rEdits, in order, to a copy of @p rMesh's regions.
 * @throws std::invalid_argument on a bad edit (see the file comment); the input
 *         is never modified.
 */
MESHIOPLUSPLUS_API Mesh edit_regions(const Mesh& rMesh, const std::vector<RegionEdit>& rEdits);

}  // namespace meshioplusplus

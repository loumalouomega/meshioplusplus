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
 * @file operations/voxelize.hpp
 * @brief Regular hexahedron grids: `grid` builds one from nothing, `voxelize`
 * builds one around a mesh and optionally marks which cells the mesh occupies.
 *
 * ### The output is an ordinary `Mesh`, and that is the whole design
 *
 * A voxel grid is one `hexahedron` cell block over a shared corner lattice --
 * not a bespoke in-memory object with its own accessors. Once it *is* a `Mesh`,
 * every writer already handles it, `view`/`screenshot` already render it,
 * `crop`/`split`/`data calc` already operate on it, `--color-by` already colours
 * it, and `isosurface` already contours whatever field it carries. None of that
 * needed a line of new code, and none of it would have been available from a
 * dedicated grid type.
 *
 * `custom` was considered for the adaptive case and rejected: `CellType::Custom`
 * reports -1 nodes and -1 dimension, which makes the block invisible to `stats`,
 * `quality`, `surface`, `gradient` and `refine` and unwritable by most formats --
 * forfeiting exactly the property that motivates the choice. A new `voxel` cell
 * type was rejected too: VTK's type 11 is deliberately unmapped in this codebase,
 * and adding it would mean a row in all 42 format tables to buy an implicit node
 * ordering nothing needs.
 *
 * ### `grid` is a primitive constructor, not a byproduct
 *
 * The lattice `voxelize` needs is also the mesh-generation primitive the library
 * has never had -- everything else transforms a mesh you already have. It is
 * exposed as `grid` for that reason, over the same
 * `detail/grid_lattice.hpp` that owns the numbering.
 *
 * @see doc/voxelize.md, detail/grid_lattice.hpp
 */

// System includes
#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/operations/sdf.hpp"

namespace meshioplusplus {

/// The Int64 0/1 `cell_data` array marking cells the input mesh occupies.
inline constexpr const char* kVoxelOccupancyName = "voxel:occupancy";

/// Which cells of the lattice to keep.
enum class VoxelFill {
    /// Every cell of the bounding box. The input mesh contributes only its
    /// bounding box, so this is a background grid rather than a shape.
    All = 0,
    /// Only cells a surface triangle passes through, by exact triangle/box
    /// overlap. Needs no distance and no closed surface -- it works on an open
    /// sheet, which `Inside` cannot.
    Surface = 1,
    /// Only cells whose centre is inside the surface. Needs a sign, and
    /// therefore a surface that is closed enough for the chosen `SdfSign`.
    Inside = 2,
};

/// Parse a `VoxelFill` name (`all`, `surface`, `inside`).
MESHIOPLUSPLUS_API VoxelFill voxel_fill_from_name(const std::string& rName);

/// Options for `voxelize`.
struct VoxelOptions {
    /// Cell counts per axis. Exactly one of this and `mCellSize` must be set.
    std::optional<std::array<std::int64_t, 3>> mResolution;
    /// Cell size (cubic). Exactly one of this and `mResolution` must be set.
    std::optional<double> mCellSize;
    /// Explicit bounds `{lo[3], hi[3]}`; unset uses the input's bounding box.
    std::optional<std::array<double, 6>> mBounds;
    /// Padding added to every side, in world units.
    double mPadding = 0.0;
    /// Padding added to every side, as a fraction of the bounding-box diagonal.
    double mPaddingRelative = 0.0;
    /// Which cells to keep.
    VoxelFill mFill = VoxelFill::All;
    /// Attach `voxel:occupancy` instead of dropping the cells it would mark.
    /// With `VoxelFill::All` this is the only way occupancy is reported at all.
    bool mAttachOccupancy = false;
    /// Refuse to generate more cells than this, naming the option. The default
    /// is a little above 256^3, which is ~1.5 GB of points and connectivity;
    /// 512^3 would be ~11.8 GB, which is why there is a limit at all.
    std::int64_t mMaxCells = 20000000;
    /// How `VoxelFill::Inside` decides what is inside.
    SurfaceDistanceOptions mDistance;
};

/// The result of `voxelize`: the grid and the lattice it came from.
struct VoxelResult {
    Mesh mMesh;                                       ///< The grid.
    std::array<std::int64_t, 3> mDims{{0, 0, 0}};     ///< Cell counts per axis.
    std::array<double, 3> mOrigin{{0.0, 0.0, 0.0}};   ///< Lattice lo corner.
    std::array<double, 3> mSpacing{{0.0, 0.0, 0.0}};  ///< Cell size per axis.
    /// How many cells the fill rule kept (equal to the total for `All`).
    std::int64_t mNumOccupied = 0;
};

/**
 * @brief Build a regular hexahedron lattice from nothing.
 * @param rDims cell counts per axis; a zero on any axis yields an empty mesh.
 * @param rOrigin the lo corner.
 * @param rSpacing cell size per axis.
 * @return a mesh with one `hexahedron` block, points x-fastest.
 * @throws std::invalid_argument on a negative count or a non-positive spacing,
 *         and when the cell count exceeds @p MaxCells.
 */
MESHIOPLUSPLUS_API Mesh grid(const std::array<std::int64_t, 3>& rDims,
                             const std::array<double, 3>& rOrigin = {{0.0, 0.0, 0.0}},
                             const std::array<double, 3>& rSpacing = {{1.0, 1.0, 1.0}},
                             std::int64_t MaxCells = 20000000);

/**
 * @brief Build a regular grid over @p rMesh and mark or keep the cells it fills.
 * @param rMesh the mesh to voxelize; used for its bounding box, and for its
 *        surface when `mFill` is not `All`.
 * @param rOptions resolution/cell size, bounds, padding and the fill rule.
 * @return the grid plus the lattice geometry and the occupied-cell count.
 * @throws std::invalid_argument when neither or both of `mResolution` and
 *         `mCellSize` are set, on a non-positive resolution or cell size, on
 *         inverted bounds, and when the lattice would exceed `mMaxCells`.
 */
MESHIOPLUSPLUS_API VoxelResult voxelize(const Mesh& rMesh, const VoxelOptions& rOptions = {});

}  // namespace meshioplusplus

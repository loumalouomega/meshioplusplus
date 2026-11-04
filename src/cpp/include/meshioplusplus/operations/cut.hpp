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
 * @file operations/cut.hpp
 * @brief Cut a mesh by level sets: the pieces of its cells on one side of one
 * or several level sets, or on their zero sets, and quadrature rules on those
 * pieces (cut-cell quadrature for unfitted finite elements).
 *
 * The cutting is CutCells' (S. Claus, MIT), an optional dependency
 * (`MESHIOPLUSPLUS_WITH_CUTCELLS`; `cut_has_cutcells()`). Without it both
 * operations throw `Unsupported` naming the option.
 *
 *  - **Level sets.** Each is a scalar `point_data` array (`mLevelSets`), one
 *    value per node, or a plane (`mPlanes`: `(x - origin) . normal`, named
 *    `plane0`, `plane1`, ...). On `triangle6`/`tetra10` cells a level set is
 *    quadratic (all ten or six nodes); on every other cell it is linear in the
 *    cell's corners.
 *  - **Selection.** `mSelection` is CutCells' expression over the level sets'
 *    names: `phi < 0`, `phi > 0`, `phi = 0` (the interface), joined with `and`
 *    and `or`, such as `"wall < 0 and hole > 0"`. Empty means `<first> < 0`.
 *  - **Cells.** The cells of the highest topological dimension are cut
 *    (`line`, `triangle`, `quad`, `tetra`, `hexahedron`, `wedge`, `pyramid`,
 *    and their higher-order forms through their corners); lower-dimensional
 *    blocks are left out with a warning. Polygons and polyhedra are refused.
 *    Geometry is straight-sided: a curved higher-order cell is cut through
 *    its corners.
 *  - **Backends.** `"lut"` cuts by lookup tables into straight pieces of the
 *    cell's own family (a hexahedron cut by a plane through two of its faces
 *    stays one hexahedron); `"quadrays"` integrates curved interfaces to high
 *    order and gives Lagrange pieces of degree `mDegree`. `"auto"` takes
 *    quadrays for a quadratic level set, the lookup tables otherwise.
 */

// System includes
#include <array>
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// A plane level set, `(x - mOrigin) . mNormal` (the normal is normalised).
struct CutPlane {
    std::array<double, 3> mOrigin = {0.0, 0.0, 0.0};
    std::array<double, 3> mNormal = {1.0, 0.0, 0.0};
};

struct CutOptions {
    /// Scalar `point_data` arrays that are level sets; their names are the
    /// names `mSelection` uses.
    std::vector<std::string> mLevelSets;
    /// Plane level sets, named `plane0`, `plane1`, ... after `mLevelSets`.
    std::vector<CutPlane> mPlanes;
    /// The part to keep; empty means `<first level set> < 0`.
    std::string mSelection;
    /// Keep the whole cells that lie wholly in the part, besides the pieces of
    /// the cut ones.
    bool mIncludeUncut = true;
    /// `auto`, `lut` or `quadrays` (see the file comment).
    std::string mBackend = "auto";
    /// quadrays only: the Lagrange degree of the pieces `cut_level_set` gives.
    std::int32_t mDegree = 2;
    /// lut only: split the pieces, and the whole cells kept with them, into
    /// simplices (tetrahedra or triangles).
    bool mTriangulate = false;
    /// `cut_quadrature` only: the polynomial degree the rules integrate exactly.
    std::int32_t mOrder = 2;
    /// Write `cut:parent_cell` (`cut_level_set`) and `quadrature:parent_cell`.
    bool mRecordParentIds = true;
};

struct CutLevelSetResult {
    /// The pieces, and the whole cells with `mIncludeUncut`: point data
    /// interpolated in each piece's parent cell, cell data and cell regions
    /// carried from it, cell data `cut:parent_cell` (the input cell, global
    /// block-major) and `cut:cut` (1 on a piece of a cut cell).
    Mesh mMesh;
    std::int64_t mNumCells = 0;       ///< input cells of the cut dimension
    std::int64_t mNumCutCells = 0;    ///< of them, those a level set crosses
    std::int64_t mNumUncutCells = 0;  ///< of them, those wholly in the part
    std::int64_t mNumPieces = 0;      ///< output cells
    std::int32_t mDimension = 0;      ///< of the part: the cells' or one less
    double mMeasure = 0.0;            ///< volume, area or length of the part
    std::string mSelection;           ///< as applied
    std::string mBackend;             ///< as applied: `lut` or `quadrays`
    std::vector<std::string> mWarnings;
};

struct CutQuadratureResult {
    /// One `vertex` per quadrature point, at its physical position, with point
    /// data `quadrature:weight` (physical), `quadrature:reference` (in the
    /// parent cell's reference element), `quadrature:rule` and, with
    /// `mRecordParentIds`, `quadrature:parent_cell`.
    Mesh mMesh;
    std::int64_t mNumRules = 0;  ///< one per cell piece the rules cover
    std::int64_t mNumPoints = 0;
    double mSumWeights = 0.0;  ///< the measure of the part, to the rules' accuracy
    std::int32_t mDimension = 0;
    std::string mSelection;
    std::string mBackend;
    std::vector<std::string> mWarnings;
};

/// Whether this build carries CutCells.
MESHIOPLUSPLUS_API bool cut_has_cutcells();

/**
 * @brief The part of `rMesh` that `rOptions.mSelection` picks (see the file
 * comment).
 * @throws std::invalid_argument for no level set, a missing or non-scalar
 *         array, a bad selection, an unknown backend, or no cell to cut.
 * @throws Unsupported without CutCells, or for polygons and polyhedra.
 */
MESHIOPLUSPLUS_API CutLevelSetResult cut_level_set(const Mesh& rMesh,
                                                   const CutOptions& rOptions = {});

/**
 * @brief Quadrature rules of degree `rOptions.mOrder` on that part.
 * @throws as `cut_level_set`.
 */
MESHIOPLUSPLUS_API CutQuadratureResult cut_quadrature(const Mesh& rMesh,
                                                      const CutOptions& rOptions = {});

}  // namespace meshioplusplus

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
 * @file operations/select.hpp
 * @brief Write a named `Region` (of points, cells or cell sides) chosen by
 * geometry, by a data or quality threshold, by connected component, or by
 * growing a smooth patch from a seed point -- so boundary-condition and contact
 * patches are reproducible from a script rather than painted by hand.
 *
 * ### Contract
 *
 * `select_region` returns a copy of the input with **one region added**:
 * points, connectivity, block order, data and every other region are
 * bit-identical. A region of the same name and kind is replaced (whatever its
 * dimension or tag), so running the same selection twice leaves one region.
 * The region's dimension is `-1` for `Point` and `Cell` regions and the facet
 * dimension (2 for faces, 1 for edges) for `Side` regions; its tag is
 * `SelectOptions::mTag`.
 *
 * ### The candidates
 *
 *  - `Point`: every point.
 *  - `Cell`: every cell (global, block-major index).
 *  - `Side`: the facets of the cells of the mesh's highest dimension -- faces
 *    of 3-D cells, or edges of 2-D cells when there are none --, named
 *    `(cell, local facet)` in the numbering `doc/regions.md` gives. With
 *    `mBoundaryOnly` (the default) only facets owned by one such cell are
 *    candidates. Cells whose type has no facet table (the VTK Lagrange types)
 *    contribute no sides.
 *
 * ### The criteria
 *
 * At least one is required; a candidate is selected when it passes **every**
 * criterion that is set.
 *
 *  - **Box, plane, sphere** (`mBox`, `mPlane`, `mSphere`): together they
 *    form one region, the intersection of those set. A point is inside when
 *    it lies in the box (`lo <= p <= hi`), on the non-negative side of the
 *    plane (`(p - point) . normal >= 0`) and in the closed ball. A cell or
 *    side passes by `mMode`: all of its nodes inside (`All`), any of them
 *    (`Any`), or its centroid, the mean of its distinct nodes (`Centroid`).
 *  - **Where** (`mWhereArray`, `mWhereOp`, `mWhereValue`): a comparison on a
 *    scalar `point_data` array for points, or a scalar `cell_data` array for
 *    cells and sides (a side takes its owner cell's value). The comparison is
 *    `refine_compare_value`'s, so a NaN never matches.
 *  - **Quality** (`mQualityMetric`, ...): the same comparison on one of
 *    `compute_quality`'s metrics (`scaled_jacobian`, `min_dihedral`, ...),
 *    computed on the fly (with or without the `quality:` prefix `attach_quality`
 *    gives it); cells and sides only. A metric that does not apply
 *    to a cell is NaN there, so the cell never matches.
 *  - **Component** (`mComponent >= 0`): the cells of the k-th node-connected
 *    component, numbered in the order of each component's lowest cell -- the
 *    numbering `split(by="component")` gives its pieces. Points pass when a
 *    cell of the component uses them, sides when their owner is in it.
 *  - **Grow** (`mSeed`): a patch grown from the facet nearest `mSeed` (by
 *    centroid; ties: the lowest index) across the edges `feature_edges` would
 *    not draw -- an edge is crossed when exactly two facets use it and their
 *    dihedral is at most `mFeatureAngleDeg` --, and, when
 *    `mNormalToleranceDeg > 0`, only onto facets whose normal is within that
 *    angle of the seed's. The facets grown over are the 2-D cells for `Cell`
 *    regions; the boundary faces of the 3-D cells (or every candidate face
 *    with `mBoundaryOnly` off) for `Side` regions; and, for `Point` regions,
 *    the 2-D cells or, when there are none, the boundary faces, whose nodes
 *    are selected. Growing needs faces: a `Side` selection on a mesh whose
 *    highest dimension is 2 (whose sides are edges) refuses it.
 *
 * Every criterion is evaluated over all candidates, independently, and the
 * results intersected -- a patch is grown over the whole surface even when a
 * box also restricts the selection.
 *
 * ### Determinism
 *
 * Per-candidate tests run in parallel into disjoint slots; the growth is a
 * serial breadth-first walk with neighbours in ascending order. The region is
 * the same on every backend and thread count.
 *
 * ### Feeding the interface operations
 *
 * `Cell` regions are what `find_interface`, `region_adjacency` and
 * `contact_pairs` select by; `Side` regions are what `split_interface` takes
 * and what solver input decks call surfaces or sidesets.
 *
 * Everything is standard C++ over the uniform mesh API, so it compiles under
 * every mesh backend.
 */

// System includes
#include <array>
#include <cstdint>
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/operations/refine.hpp"
#include "meshioplusplus/region.hpp"

namespace meshioplusplus {

/// How a cell or side reduces its nodes to one verdict for a box, plane or sphere.
enum class SelectMode : std::int32_t {
    All = 0,       ///< every node inside
    Any = 1,       ///< at least one node inside
    Centroid = 2,  ///< the mean of its distinct nodes inside
};

/**
 * @brief Parse a `SelectMode` from `"all"`, `"any"` or `"centroid"`.
 * @throws std::invalid_argument naming the accepted values.
 */
MESHIOPLUSPLUS_API SelectMode select_mode_from_name(const std::string& rName);

/// The lower-case spelling of a `SelectMode`.
MESHIOPLUSPLUS_API const char* select_mode_name(SelectMode Mode);

/// What `select_region` selects and how.
struct SelectOptions {
    /// The region's name; must not be empty.
    std::string mName;
    /// The kind of region written.
    RegionKind mKind = RegionKind::Cell;
    /// The region's tag (`Region::mTag`); `-1` for none.
    std::int64_t mTag = -1;
    /// How a cell or side reduces its nodes for the box, plane and sphere.
    SelectMode mMode = SelectMode::All;
    /// `Side` regions: only facets owned by one cell are candidates.
    bool mBoundaryOnly = true;

    /// Box criterion: `mBoxLo <= p <= mBoxHi` component-wise.
    bool mBox = false;
    std::array<double, 3> mBoxLo{0.0, 0.0, 0.0};
    std::array<double, 3> mBoxHi{0.0, 0.0, 0.0};

    /// Half-space criterion: `(p - mPlanePoint) . mPlaneNormal >= 0`.
    bool mPlane = false;
    std::array<double, 3> mPlanePoint{0.0, 0.0, 0.0};
    std::array<double, 3> mPlaneNormal{0.0, 0.0, 1.0};

    /// Ball criterion: `|p - mSphereCenter| <= mSphereRadius`.
    bool mSphere = false;
    std::array<double, 3> mSphereCenter{0.0, 0.0, 0.0};
    double mSphereRadius = 0.0;

    /// Data criterion: a scalar array compared with `mWhereValue`; empty disables it.
    std::string mWhereArray;
    RefineCompare mWhereOp = RefineCompare::Less;
    double mWhereValue = 0.0;

    /// Quality criterion: a `compute_quality` metric; empty disables it.
    std::string mQualityMetric;
    RefineCompare mQualityOp = RefineCompare::Less;
    double mQualityValue = 0.0;

    /// Component criterion: the component number; negative disables it.
    std::int64_t mComponent = -1;

    /// Grow criterion: a patch grown from the facet nearest `mSeedPoint`.
    bool mSeed = false;
    std::array<double, 3> mSeedPoint{0.0, 0.0, 0.0};
    /// The largest dihedral angle, in degrees, the growth crosses; in `[0, 180]`.
    double mFeatureAngleDeg = 30.0;
    /// The largest angle, in degrees, between a facet's normal and the seed's;
    /// `0` disables the test. In `[0, 180]`.
    double mNormalToleranceDeg = 0.0;
};

/// The mesh with its new region, and what the selection found.
struct SelectResult {
    /// The input with the region added (or replaced).
    Mesh mMesh;
    /// The region written, canonical (sorted, no duplicates).
    meshioplusplus::Region mRegion;
    /// Candidates considered (points, cells or sides).
    std::int64_t mNumCandidates = 0;
    /// Candidates selected: the region's size.
    std::int64_t mNumSelected = 0;
    /// Facets the grown patch covers, or -1 when no seed was given.
    std::int64_t mNumGrown = -1;
    /// Node-connected components, or -1 when no component was asked for.
    std::int64_t mNumComponents = -1;
};

/**
 * @brief Select points, cells or sides and write them as a named region.
 *
 * @param rMesh the input mesh.
 * @param rOptions see `SelectOptions`.
 * @return the mesh with the region, and counts.
 * @throws std::invalid_argument on an empty name, no criterion, an option out
 *         of range, an unknown or non-scalar array, an unknown quality metric,
 *         a quality criterion on a `Point` region, a seed with no facet to grow
 *         over, or a `Side` region on a mesh with no 2-D or 3-D cells.
 */
MESHIOPLUSPLUS_API SelectResult select_region(const Mesh& rMesh, const SelectOptions& rOptions);

}  // namespace meshioplusplus

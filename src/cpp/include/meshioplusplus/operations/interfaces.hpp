//  ██████   ██████ ██████████  █████████  █████   █████ █████    ███████
// ░░██████ ██████ ░░███░░░░░█ ███░░░░░███░░███   ░░███ ░░███   ███░░░░░███
//  ░███░█████ ░███  ░███  █ ░ ░███    ░░░  ░███    ░███  ░███  ███     ░░███
//  ░███░░███ ░███  ░██████   ░░█████████  ░███████████  ░███ ░███      ░███
//  ░███ ░░░  ░███  ░███░░█    ░░░░░░░░███ ░███░░░░░███ ░███ ░███      ░███
//  ░███      ░███  ░███ ░   █ ███    ░███ ░███    ░███ ░███ ░░███     ███
//  █████     █████ ██████████░░█████████  █████   █████ █████ ░░░███████░
// ░░░░░     ░░░░░ ░░░░░░░░░░  ░░░░░░░░░  ░░░░░   ░░░░░ ░░░░░    ░░░░░░░
//
//  License:         MIT License
//                   meshio++ default license: LICENSE
//
//  Main authors:    Vicente Mataix Ferrandiz
//
#pragma once

/**
 * @file operations/interfaces.hpp
 * @brief Adjacency, interface search, contact projection and interface splitting.
 *
 * `region_adjacency` emits conforming shared facets. `find_interface` locates
 * matching boundary facets either exactly or by geometric proximity.
 *
 * The returned mesh retains the input points and point/field data. Each output
 * facet carries the selected-region indices, representative source cell/facet
 * pairs, the number of participating source cells, and its measure. Indices
 * refer to the selector order (or the sorted Cell-region order when selectors
 * are omitted). One Cell region named `adjacency:<a>:<b>` groups each pair's
 * output facets. Point regions are retained; input Cell and Side regions are
 * not copied because their entity indices refer to the source mesh.
 */

// System includes
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/operations/region_ops.hpp"

namespace meshioplusplus {

/**
 * @brief Build the shared-facet mesh for pairs of selected Cell regions.
 *
 * @param rMesh Source mesh. Rectangular linear and higher-order cells,
 *        polygons and polyhedra contribute their corner facets.
 * @param rRegions Cell-region selectors. Empty selects every Cell region in
 *        the mesh's canonical region order. If fewer than two Cell regions
 *        exist, multiple cell blocks become groups. Selectors can pin kind,
 *        dimension and tag when names repeat.
 * @return A mesh containing one output facet for each conforming shared facet
 *         and selected region pair. An empty result is valid.
 * @throws std::invalid_argument when fewer than two regions are selected, a
 *         selector is missing/ambiguous, or a selected region is not a Cell
 *         region.
 */
MESHIOPLUSPLUS_API Mesh region_adjacency(const Mesh& rMesh,
                                         const std::vector<RegionSelector>& rRegions = {});

/// How `find_interface` matches two selected parts.
enum class InterfaceMode { Conforming = 0, Proximity = 1 };

/// Which side supplies the output facets of `find_interface`.
enum class InterfaceMaster { A = 0, B = 1 };

/// Options for `find_interface`.
struct FindInterfaceOptions {
    InterfaceMode mMode = InterfaceMode::Conforming;
    InterfaceMaster mMaster = InterfaceMaster::A;
    /// Maximum closest-point gap. Zero derives 1% of the mean boundary edge.
    double mGapTolerance = 0.0;
    /// Maximum deviation from opposing normals, in degrees.
    double mAngleTolerance = 30.0;
    /// Additional distance admitted for small overlaps. Non-negative.
    double mOverlapTolerance = 0.0;
};

/// Summary returned by `find_interface`.
struct InterfaceReport {
    std::int64_t mNumPairs = 0;
    double mArea = 0.0;
    double mMaxGap = 0.0;
    std::int64_t mUnmatchedA = 0;
    std::int64_t mUnmatchedB = 0;
};

/// Matched master facets plus Side regions addressing both input meshes.
struct FindInterfaceResult {
    Mesh mMesh;
    /// `Side` entries `(global_cell, local_facet)` on the A input mesh.
    Region mSideA;
    /// `Side` entries `(global_cell, local_facet)` on the B input mesh.
    Region mSideB;
    InterfaceReport mReport;
};

/**
 * @brief Find interfaces between two Cell regions of one mesh.
 * @throws std::invalid_argument for invalid regions, modes or tolerances.
 */
MESHIOPLUSPLUS_API FindInterfaceResult find_interface(const Mesh& rMesh,
                                                      const RegionSelector& rRegionA,
                                                      const RegionSelector& rRegionB,
                                                      const FindInterfaceOptions& rOptions = {});

/**
 * @brief Find interfaces between Cell regions of two meshes.
 *
 * In conforming mode the regions must already share point ids (typically after
 * welding/merging). Proximity mode supports independently numbered surfaces.
 */
MESHIOPLUSPLUS_API FindInterfaceResult find_interface(const Mesh& rMeshA,
                                                      const RegionSelector& rRegionA,
                                                      const Mesh& rMeshB,
                                                      const RegionSelector& rRegionB,
                                                      const FindInterfaceOptions& rOptions = {});

/// Options for the node-to-facet projection performed by `contact_pairs`.
struct ContactPairsOptions {
    /// Maximum accepted distance; zero derives 1% of the master mean edge.
    double mTolerance = 0.0;
    /// If true, a slave point with no master facet in range is an error.
    bool mRequireComplete = false;
};

/// The nearest master-facet projection for each selected slave point.
struct ContactPairsResult {
    NDArray mSlavePoint;        ///< Int64 `(n,)`, source slave point ids.
    NDArray mMasterCell;        ///< Int64 `(n,)`, global master cell ids (`-1` unmatched).
    NDArray mMasterFacet;       ///< Int64 `(n,)`, local master facet ids (`-1` unmatched).
    NDArray mMasterSubfacet;    ///< Int64 `(n,)`, fan-triangle ordinal (`0` for edges/triangles).
    NDArray mLocalCoordinates;  ///< Float64 `(n, 3)`, barycentrics/edge weights in subfacet.
    NDArray mClosestPoint;      ///< Float64 `(n, 3)`, xyz-padded projection.
    NDArray mGap;               ///< Float64 `(n,)`, signed along the master normal.
    NDArray mNormal;            ///< Float64 `(n, 3)`, unit master normal.
    NDArray mUnmatched;         ///< Int64 `(k,)`, unmatched slave point ids.
};

/**
 * @brief Project points in a Point region to the closest facets of a Cell region.
 * The two meshes may be the same object. A volume Cell region contributes its
 * boundary; a 3-D surface Cell region contributes its own cells.
 */
MESHIOPLUSPLUS_API ContactPairsResult contact_pairs(const Mesh& rSlave,
                                                    const RegionSelector& rSlavePoints,
                                                    const Mesh& rMaster,
                                                    const RegionSelector& rMasterCells,
                                                    const ContactPairsOptions& rOptions = {});

/// Cohesive-element insertion options.
struct SplitInterfaceOptions {
    bool mAddCohesive = false;
};

/// The split mesh and deterministic duplication counts.
struct SplitInterfaceResult {
    Mesh mMesh;
    std::int64_t mNumDuplicatedPoints = 0;
    std::int64_t mNumCohesiveCells = 0;
};

/**
 * @brief Duplicate point fans along the facets named by a Side region.
 *
 * Polyhedra are rejected. With cohesive insertion, triangle/quad interfaces
 * become wedge/hexahedron cells; a 2-D cohesive `line` stores the opposite
 * duplicated trace as `cohesive:trace_b` (Int64, two components per cell).
 */
MESHIOPLUSPLUS_API SplitInterfaceResult split_interface(const Mesh& rMesh,
                                                        const RegionSelector& rSide,
                                                        const SplitInterfaceOptions& rOptions = {});

}  // namespace meshioplusplus

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
 * @file operations/check_conformity.hpp
 * @brief Whether a mesh is conforming: hanging nodes, T-junctions and facets
 * that only partly overlap, reported with a flag per point.
 *
 * A mesh is conforming when every two cells meet in a whole shared facet, edge
 * or corner, with the same point ids on both sides. This check reports the ways
 * that fails on the mesh's top-dimensional cells (volume cells when there are
 * any, else 2-D cells; lines and vertices are ignored):
 *
 *  - a **hanging node** (3-D): a point some cell uses that lies on a face of a
 *    volume cell -- in its interior or on one of its edges -- without being one
 *    of that cell's nodes;
 *  - a **T-junction** (2-D): a point some cell uses that lies on the interior of
 *    an edge of a 2-D cell without being one of its nodes;
 *  - a **partial overlap**: two unmatched facets (faces in 3-D, edges in 2-D)
 *    that lie in one plane (on one line) and overlap over a positive area
 *    (length) without being the same facet;
 *  - an **unwelded** pair: two unmatched facets whose corners coincide but whose
 *    point ids differ.
 *
 * Only unmatched facets -- a face of one volume cell, an edge of one 2-D cell
 * -- are examined: a node hanging on a face is always a corner of an unmatched
 * facet of the finer side, and a face matched by its neighbour cannot hold one.
 * A cell's own nodes are never flagged against it, so a quadratic cell's
 * mid-edge nodes are not hanging. Two points closer than the tolerance are one
 * location: a point on a facet's corner is an unwelded duplicate, not a hanging
 * node.
 *
 * The verdict tells a caller which `find_interface` mode suits the mesh:
 * `"conforming"` when nothing is reported, `"weld"` when the only defect is
 * unwelded facets (merge the coincident points, then `"conforming"`), and
 * `"proximity"` when facets partly overlap. Points, connectivity, block order,
 * data and regions are returned bit-identical; only the flag array is added.
 *
 * Deterministic: the flags are computed per point and the counts summed in
 * facet order, so the result is identical across parallel backends and thread
 * counts.
 *
 * @see doc/check_conformity.md
 */

// System includes
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// Point data (Int64): 0, or `ConformityFlag::HangingNode`/`TJunction`.
inline constexpr const char* kConformityHangingName = "conformity:hanging";

/// The values of `conformity:hanging`.
enum class ConformityFlag : std::int32_t {
    Conforming = 0,
    HangingNode = 1,  ///< on a face of a volume cell that does not use it
    TJunction = 2,    ///< on an edge of a 2-D cell that does not use it
};

/// Options for `check_conformity`.
struct ConformityOptions {
    /// Distances at or below this fraction of the bounding-box diagonal of the
    /// examined cells' points count as zero; must be finite and non-negative.
    double mTolerance = 1e-9;
    /// Attach `conformity:hanging` to the returned mesh.
    bool mAttachFlags = true;
};

/// The verdict and the counts behind it.
struct ConformityResult {
    /// The input, bit-identical, with `conformity:hanging` when asked for.
    Mesh mMesh;
    /// 3 for volume cells, 2 for 2-D cells, 0 when the mesh has neither.
    std::int64_t mDimension = 0;
    std::int64_t mNumHangingNodes = 0;       ///< points flagged `HangingNode`
    std::int64_t mNumTJunctions = 0;         ///< points flagged `TJunction`
    std::int64_t mNumPartialOverlaps = 0;    ///< pairs of partly overlapping facets
    std::int64_t mNumUnweldedFacets = 0;     ///< pairs of coincident, unwelded facets
    std::int64_t mNumBoundaryFacets = 0;     ///< facets used by one cell
    std::int64_t mNumNonManifoldFacets = 0;  ///< facets used by three or more cells
    /// No hanging node, T-junction, partial overlap or unwelded pair.
    bool mConforming = true;
    /// `"conforming"`, `"weld"` or `"proximity"`: see the file comment.
    std::string mInterfaceMode = "conforming";
    /// What the check could not examine, one entry per line.
    std::vector<std::string> mWarnings;
};

/**
 * @brief Report the hanging nodes, T-junctions, partial overlaps and unwelded
 * facets of @p rMesh.
 * @param rMesh any mesh; its volume cells (else its 2-D cells) are examined.
 * @param rOptions the tolerance and whether to attach the flags.
 * @throws std::invalid_argument on a negative or non-finite tolerance.
 */
MESHIOPLUSPLUS_API ConformityResult check_conformity(const Mesh& rMesh,
                                                     const ConformityOptions& rOptions = {});

}  // namespace meshioplusplus

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
 * @file detail/facet_index.hpp
 * @brief Find the cell facet a file names by its nodes.
 *
 * Several formats name a boundary facet by its node list rather than by
 * (cell, local facet): LS-DYNA segments, FEBio surfaces, Elmer boundary
 * elements. `FacetIndex` keys every facet of a mesh by its sorted corner nodes
 * (at most four: every cell type's faces and edges are triangles, quads or
 * lines) and answers, for a node list, which cells own that facet and under which
 * local facet number — the numbering a `Side` region uses (`cell_faces` for a
 * 3-D cell, `cell_edges` for a 2-D one; see `doc/regions.md`). An interior
 * facet has two owners; the first two in block-major order are kept, with the
 * total count.
 *
 * `facet_nodes` is the other direction: the node list of one (cell, facet).
 */

// System includes
#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/detail/polyhedron.hpp"

namespace meshioplusplus {
namespace detail {

/// Which facets a `FacetIndex` holds.
struct FacetIndexOptions {
    /// The faces of every 3-D cell (`cell_faces`).
    bool mSolidFaces = true;
    /// The edges of every 2-D cell (`cell_edges`).
    bool mSurfaceEdges = true;
    /// A `triangle`/`quad` cell's own face, as facet 0. Not a `Side` region
    /// facet (a 2-D cell's facets are its edges): LS-DYNA's shell segments
    /// and FEBio's shell surfaces use it to recognise a shell by its nodes.
    bool mSurfaceSelf = false;
};

/// One owner of a facet: a global (block-major) cell index and a local facet.
struct FacetOwner {
    std::int64_t mCell = -1;
    std::int64_t mFacet = -1;
};

/// The owners of one facet: the first two in block-major order, and how many.
struct FacetHit {
    FacetOwner mFirst;
    FacetOwner mSecond;
    std::size_t mCount = 0;
};

/// Sorted-corner lookup of every facet of a mesh.
class MESHIOPLUSPLUS_API FacetIndex {
public:
    /**
     * @brief Index the facets of @p rMesh.
     * @param rMesh The mesh; ragged and polyhedral blocks are skipped.
     * @param rOptions Which facets to hold.
     */
    explicit FacetIndex(const Mesh& rMesh, const FacetIndexOptions& rOptions = {});

    /**
     * @brief The owners of the facet whose corners are @p pCorners.
     * @param pCorners Corner node ids, in any order (mid-side nodes excluded).
     * @param N Number of corners.
     * @return The hit, or `nullptr` when no indexed facet has these corners.
     */
    const FacetHit* Find(const std::int64_t* pCorners, std::size_t N) const;

    /// Number of distinct facets held.
    std::size_t Size() const { return mHits.size(); }

private:
    // Sorted facet keys -- the corner count, then the sorted corners,
    // zero-padded -- and each one's hit, found by binary search. A sorted
    // table, not a hash map, since v16.16.0 (ABI 18): it is built by a
    // parallel sort with the same hits a serial insert gave.
    std::vector<std::array<std::int64_t, 5>> mKeys;
    std::vector<FacetHit> mHits;
};

/**
 * @brief The nodes of one `Side` facet: its cell type and node ids, corners
 * first, then mid-side and centre nodes as `cell_faces`/`cell_edges` list them.
 * @param rMesh The mesh.
 * @param Cell Global (block-major) cell index.
 * @param Facet Local facet: a face of a 3-D cell, an edge of a 2-D one.
 * @param rType Receives the facet's cell type.
 * @param rNodes Receives the facet's node ids.
 * @return `false` when the cell or facet does not exist (nothing is written).
 */
MESHIOPLUSPLUS_API bool facet_nodes(const Mesh& rMesh, std::int64_t Cell, std::int64_t Facet,
                                    CellType& rType, std::vector<std::int64_t>& rNodes);

}  // namespace detail
}  // namespace meshioplusplus

#pragma once

/**
 * @file detail/hex_dual_schemes.hpp
 * @brief The all-hexahedral dual of a paired, strongly balanced octree, after
 * the transition schemes of Livesu, Pitzalis and Cherchi ("Optimal dual
 * schemes for adaptive grid based hexmeshing", ACM TOG 41(2), 2022), ported
 * from CinoLib (MIT).
 *
 * The leaves are the primal mesh: hexahedra with hanging nodes where a leaf
 * meets leaves one level finer. A **transition vertex** is the centre of a
 * face of a parent whose eight children are leaves and whose neighbour across
 * that face is split two levels down on that side; every coarse leaf around
 * the transitions is replaced by one of eight polyhedral schemes, placed by
 * the transition vertices it holds. The result is a conforming polyhedral mesh
 * whose every interior vertex has eight cells around it, so its dual is all
 * hexahedra: one per interior primal vertex, with a corner at the centroid of
 * each cell around it.
 *
 * The port keeps CinoLib's detection, marking order and orientation rules,
 * with every coordinate on an integer lattice of 1/20 of the finest cell (the
 * schemes use 0, 0.3, 0.5 and 1 of a cell), so vertices weld by exact keys
 * and every rotation and reflection is a signed permutation. Every assertion
 * of the original is an `Unsupported` error naming the configuration, and the
 * schemes for weakly balanced grids are not ported: the tree must be balanced
 * across faces, edges and corners and paired (`linear_octree_pair`).
 *
 * A **core-private** header (the `cut_graph.hpp` precedent): nothing installed
 * names it.
 */

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "linear_octree.hpp"
#include "meshioplusplus/detail/geometry.hpp"

namespace meshioplusplus {
namespace detail {

/// The schemes for strongly balanced grids, in CinoLib's order.
enum class HexTransition {
    Flat,
    FlatConvex,
    Convex1,
    Convex2,
    Convex3,
    Edge,
    VertCenter,
    VertSide,
};

constexpr std::size_t kNumHexTransitions = 8;

/// One scheme: a unit cube cut into polyhedra, coordinates times 20.
struct HexDualScheme {
    const std::int8_t* mVerts;  ///< 3 per vertex, each in [0, 20]
    std::size_t mNumVerts;
    const std::uint16_t* mFaces;     ///< vertex ids of every face, concatenated
    const std::uint16_t* mFaceEnds;  ///< one past each face's last entry
    std::size_t mNumFaces;
    const std::uint16_t* mPolys;     ///< face ids of every cell, concatenated
    const std::uint16_t* mPolyEnds;  ///< one past each cell's last entry
    std::size_t mNumPolys;
    const std::uint8_t* mWinding;  ///< per entry of `mPolys`
};

const HexDualScheme& hex_dual_scheme(HexTransition Type);

/// The dual mesh, in lattice units (1/20 of the finest cell) from the forest's
/// lo corner.
struct HexDualMesh {
    std::vector<Vec3> mPoints;
    std::vector<std::array<std::int64_t, 8>> mHexes;  ///< VTK order, positive volume
    /// The primal vertex each hexahedron is the dual of.
    std::vector<std::array<std::int64_t, 3>> mCentre;
    /// The finest level among the leaves around that vertex.
    std::vector<std::int32_t> mLevel;
    std::int64_t mNumTransitionVertices = 0;
    std::array<std::int64_t, kNumHexTransitions> mNumSchemes{};
    std::int64_t mNumPrimalCells = 0;
};

/// The lattice units per finest cell.
constexpr std::int64_t kHexDualUnit = 20;

/**
 * @brief The all-hexahedral dual of @p rTree, which must be paired and
 * balanced across faces, edges and corners. Vertices on the forest's outer
 * boundary have no hexahedron.
 * @throws Unsupported naming the configuration when a scheme cannot be
 * placed, the schemes do not conform, or an interior vertex's dual is not a
 * hexahedron; never returns a partial mesh.
 */
HexDualMesh hex_dual_mesh(const LinearOctree& rTree);

}  // namespace detail
}  // namespace meshioplusplus

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

#pragma once

/**
 * @file detail/tet_remesh.hpp
 * @brief The editable tetrahedral store behind `adapt_volume`: tetrahedra with
 * tombstones, each vertex's incident tetrahedra, the tagged boundary faces
 * (outer boundary and material interfaces), the ridge edges, and the
 * classification of every vertex as a corner, a ridge vertex (on a chain), a
 * surface vertex (in a patch) or an interior vertex.
 *
 * The store is plain data with small, exact queries; the operation composes
 * its edits from them (precondition, a simulated edit whose postconditions are
 * checked, then the commit). Orientation decisions use the exact `orient3d`
 * predicate. Hash maps hold the boundary and ridge lookups only: nothing that
 * leaves the store is produced by iterating them, so the output is
 * deterministic.
 *
 * A **core-private** header: nothing installed names it.
 */

// System includes
#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

// Project includes
#include "meshioplusplus/detail/geometry.hpp"

namespace meshioplusplus {
namespace detail {

/// What a vertex may do: lower kinds are more constrained.
enum class TetVertexKind : std::uint8_t {
    Corner = 0,    ///< never moves, never removed
    Ridge = 1,     ///< moves along its chain only
    Surface = 2,   ///< moves within its patch only
    Interior = 3,  ///< free
};

/// The four corners of a tetrahedron, positively oriented in the store.
using TetCorners = std::array<std::int64_t, 4>;
/// A triangle as its three sorted vertex ids.
using TriKey = std::array<std::int64_t, 3>;

/// The sorted key of the triangle (a, b, c).
TriKey tri_key(std::int64_t A, std::int64_t B, std::int64_t C);

/// A 64-bit key for the undirected edge (a, b); both ids must be below 2^32.
inline std::uint64_t edge_key64(std::int64_t A, std::int64_t B) {
    const std::uint64_t lo = static_cast<std::uint64_t>(A < B ? A : B);
    const std::uint64_t hi = static_cast<std::uint64_t>(A < B ? B : A);
    return (lo << 32) | hi;
}

/// A hash of a `TriKey` (a 64-bit mix of its three ids).
struct TriKeyHash {
    std::size_t operator()(const TriKey& rKey) const noexcept;
};

/// A boundary face: on the outer boundary or between two material classes.
struct TetBoundaryFace {
    std::int32_t mPatch = -1;   ///< the surface patch it belongs to
    std::int32_t mTag = 0;      ///< caller's tag (side-region membership), 0 for none
    std::int32_t mClassA = -1;  ///< the class it faces (the lower one on an interface)
    std::int32_t mClassB = -1;  ///< the other class, -1 on the outer boundary
};

/// A ridge chain as an arc-length parametrised polyline of the input boundary.
struct TetRidgeChain {
    std::int64_t mFront = -1;  ///< the corner it starts at
    std::int64_t mBack = -1;   ///< the corner it ends at
    std::vector<Vec3> mPoints;
    std::vector<double> mArc;  ///< cumulative length at each point

    double Length() const { return mArc.empty() ? 0.0 : mArc.back(); }
    /// The point at arc length @p S, clamped to the chain.
    Vec3 Eval(double S) const;
    /// The arc length of the point of the chain nearest @p rP.
    double Project(const Vec3& rP) const;
};

/// The classification counts.
struct TetClassifyReport {
    std::int64_t mNumCorners = 0;
    std::int64_t mNumRidgeVertices = 0;
    std::int64_t mNumChains = 0;
    std::int64_t mNumPatches = 0;
};

/**
 * @brief The editable tetrahedral mesh. Members are public plain arrays
 * (struct of arrays) so the operation reads them without accessor noise; the
 * methods keep the incidence lists, the boundary map and the ridge map
 * consistent.
 */
class TetRemeshStore {
public:
    // --- vertices ---------------------------------------------------------
    std::vector<double> mXyz;                ///< 3 per vertex
    std::vector<double> mSize;               ///< target edge length
    std::vector<std::uint8_t> mKind;         ///< TetVertexKind
    std::vector<std::int32_t> mFeature;      ///< patch (Surface) or chain (Ridge), else -1
    std::vector<double> mChainParam;         ///< arc length on the chain (Ridge)
    std::vector<std::uint8_t> mFrozen;       ///< never moves, never removed
    std::vector<std::uint8_t> mVertexAlive;  ///< 0 once collapsed away
    std::vector<std::int64_t> mSource;       ///< input point whose integer data it keeps
    std::vector<std::vector<std::int64_t>> mVertexTets;  ///< alive incident tets, ascending

    // --- tetrahedra -------------------------------------------------------
    std::vector<TetCorners> mTets;
    std::vector<std::int32_t> mClass;     ///< material class
    std::vector<std::uint8_t> mTetAlive;  ///< tombstones

    // --- boundary ---------------------------------------------------------
    std::unordered_map<TriKey, TetBoundaryFace, TriKeyHash> mBoundary;
    std::unordered_map<std::uint64_t, std::int32_t> mRidgeEdges;  ///< edge -> chain
    std::vector<TetRidgeChain> mChains;
    std::int32_t mNumPatches = 0;

    /**
     * @brief Load the tetrahedra @p rTets over the points @p rXyz (3 per
     * point). Each tetrahedron is reoriented to a non-negative exact
     * orientation; the faces with one tetrahedron and those between two of
     * different @p rClass become boundary faces, tagged from @p rFaceTags
     * (sorted triangle -> tag) where present.
     * @throws std::invalid_argument on a face shared by more than two tetrahedra.
     */
    void Build(std::vector<double> Xyz, std::vector<TetCorners> Tets,
               std::vector<std::int32_t> Class,
               const std::unordered_map<TriKey, std::int32_t, TriKeyHash>& rFaceTags);

    /**
     * @brief Classify every vertex, number the patches and trace the chains.
     * A boundary edge is a ridge when it does not have exactly two boundary
     * faces, when its two faces differ in tag or classes, or (with
     * @p PreserveFeatures) when their normals turn by more than
     * @p FeatureAngleDeg. Corners are where the ridge graph branches, ends or
     * turns by more than @p CornerAngleDeg, and where chains of different
     * patches meet.
     */
    TetClassifyReport Classify(double FeatureAngleDeg, double CornerAngleDeg,
                               bool PreserveFeatures);

    std::size_t NumVertices() const { return mVertexAlive.size(); }
    std::size_t NumTets() const { return mTets.size(); }
    std::int64_t NumAliveTets() const;

    Vec3 Point(std::int64_t V) const {
        const std::size_t i = static_cast<std::size_t>(V) * 3;
        return {mXyz[i], mXyz[i + 1], mXyz[i + 2]};
    }
    void SetPoint(std::int64_t V, const Vec3& rP) {
        const std::size_t i = static_cast<std::size_t>(V) * 3;
        mXyz[i] = rP[0];
        mXyz[i + 1] = rP[1];
        mXyz[i + 2] = rP[2];
    }
    TetVertexKind Kind(std::int64_t V) const {
        return static_cast<TetVertexKind>(mKind[static_cast<std::size_t>(V)]);
    }

    /// True when the four points are strictly positively oriented (exact).
    static bool Positive(const Vec3& rA, const Vec3& rB, const Vec3& rC, const Vec3& rD);
    /// The exact orientation sign of the four points: 1, 0 or -1.
    static int OrientSign(const Vec3& rA, const Vec3& rB, const Vec3& rC, const Vec3& rD);

    /// The alive tetrahedra containing both @p A and @p B, ascending.
    void SharedTets(std::int64_t A, std::int64_t B, std::vector<std::int64_t>& rOut) const;
    /// The sorted vertices sharing a tetrahedron with @p V, without @p V.
    void VertexRing(std::int64_t V, std::vector<std::int64_t>& rOut) const;
    /// The boundary faces containing @p V, sorted.
    void BoundaryFacesOf(std::int64_t V, std::vector<TriKey>& rOut) const;
    /// True when @p V lies on a boundary face.
    bool OnBoundary(std::int64_t V) const;
    /// The boundary face (a, b, c), or null.
    const TetBoundaryFace* FindFace(std::int64_t A, std::int64_t B, std::int64_t C) const;
    /// The chain of the ridge edge (a, b), or -1.
    std::int32_t RidgeChain(std::int64_t A, std::int64_t B) const;
    /// The arc length of @p V on chain @p Chain (a corner reads as the end
    /// nearer @p Near).
    double ChainParam(std::int64_t V, std::int32_t Chain, double Near) const;
    /// The boundary face (a, b, c) wound outward from its class-A tetrahedron.
    std::array<std::int64_t, 3> OrientedFace(const TriKey& rKey) const;

    /// Append a vertex and return its id.
    std::int64_t AddVertex(const Vec3& rP, double Size, TetVertexKind Kind, std::int32_t Feature,
                           double ChainParam, std::int64_t Source);
    /// Append a tetrahedron and return its id.
    std::int64_t AddTet(const TetCorners& rCorners, std::int32_t Class);
    /// Tombstone a tetrahedron.
    void KillTet(std::int64_t T);

    /**
     * @brief Check the invariants: every alive tetrahedron strictly positive
     * (exact) or flagged in @p rAllowFlat, no duplicate tetrahedra, consistent
     * incidence lists, no face with three tetrahedra, every boundary or
     * interface face in the boundary map and every boundary map entry a
     * boundary face, and an outer boundary with no open edge.
     * @param rWhy set to the first violation found.
     */
    bool CheckInvariants(std::string& rWhy, const std::vector<std::uint8_t>& rAllowFlat = {}) const;

private:
    void BuildVertexTets();
};

}  // namespace detail
}  // namespace meshioplusplus

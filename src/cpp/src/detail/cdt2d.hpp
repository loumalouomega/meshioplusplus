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
 * @file detail/cdt2d.hpp
 * @brief A constrained Delaunay triangulation of a planar straight-line graph,
 * the substrate of `mesh_cad_surface` (one per face, in its parameters) and
 * of `triangulate_polygon_2d`.
 *
 *  - **Insertion.** Points go into a triangle enclosing their box (the
 *    "super triangle", vertices 0, 1, 2), each located by a walk from the
 *    last triangle made, split in three (or in four on an edge), and the
 *    Delaunay property restored by Lawson's flips (C. L. Lawson, "Software for
 *    C1 surface interpolation", 1977). A batch is inserted in Hilbert order
 *    (`detail/space_filling.hpp`), so walks are short.
 *  - **Constraints.** A segment that is not an edge is recovered as in M. V.
 *    Anglada, "An improved incremental algorithm for constructing restricted
 *    Delaunay triangulations", Computers & Graphics 21(2), 1997: the
 *    triangles it crosses are removed and the two pseudo-polygons on either
 *    side re-triangulated, each by the vertex no other vertex of it lies in
 *    the circle of. A segment through a vertex is split there; one that
 *    crosses a constrained edge is refused as `SelfIntersecting`.
 *  - **Domain.** `Classify` marks the triangles reached from the super
 *    triangle across an odd number of constrained edges as inside (the
 *    even-odd rule of `detail/brep_trim.hpp`'s inside test), so holes need no
 *    orientation; a dangling segment, which makes the parity inconsistent,
 *    is refused as `OpenLoop`.
 *  - **Steiner points.** After `Classify`, `InsertSteiner` adds a point
 *    strictly inside the domain and restores the constrained Delaunay
 *    property, never flipping a constrained edge; a point outside, on a
 *    constrained edge or on a vertex is refused.
 *  - **Segments inside the domain.** After `Classify`, `Insert` and
 *    `InsertSegment` add a curve embedded in the domain (`mesh_cad_surface`'s
 *    embedded edges): it is recovered like a boundary segment, but the
 *    classification has already been made from the loops alone, so it never
 *    changes which triangles are inside. `Where` says, without inserting,
 *    whether a point is inside, outside, on a constrained edge or on a
 *    vertex, so a caller can check a curve before committing to it.
 *  - **Crossings.** `cdt2d_find_crossings` lists the pairs of segments of a
 *    straight-line graph that meet anywhere but at a shared endpoint (a
 *    proper crossing, a touch, an overlap, a fold-back), with exact
 *    orientation tests over a bucket grid, so a caller can split them before
 *    the triangulation would refuse them as `SelfIntersecting`.
 *
 * Every orientation and in-circle decision is one of the exact predicates of
 * `detail/predicates.hpp`, and no coordinate is ever compared to another with
 * a tolerance. The predicates are exact only while no product inside them
 * overflows or underflows (Shewchuk 1997), and `incircle` raises coordinate
 * differences to the fourth power, so the coordinates must lie in the range
 * `cdt2d_in_range` accepts (zero, or a magnitude in [2^-128, 2^128]); callers
 * refuse anything outside it by name: a point is a duplicate only when it coincides exactly
 * (merging near-duplicates is the caller's choice, as `cdt2d_triangulate_polygon` does). Triangles
 * are counter-clockwise. Everything is sequential and the output depends only on the input, so it
 * is deterministic.
 *
 * A **core-private** header: nothing installed names it.
 */

// System includes
#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace meshioplusplus {
namespace detail {

enum class Cdt2dStatus : int {
    Ok = 0,
    SelfIntersecting = 1,  ///< two constrained segments cross
    OpenLoop = 2,          ///< the constraints do not bound a domain consistently
    Degenerate = 3,        ///< too few distinct points, or a zero-area domain
};

const char* cdt2d_status_name(Cdt2dStatus Status);

/// Whether a coordinate is one the triangulation stays exact for: zero, or
/// finite with a magnitude in [2^-128, 2^128] (about 2.9e-39 to 3.4e38).
/// Beyond it the super triangle's fourth powers inside `incircle` overflow;
/// below it, the tails of its products underflow.
bool cdt2d_in_range(double X);

/// Where a point lies in a classified triangulation (`Cdt2d::Where`).
enum class Cdt2dPlace : int {
    Outside = 0,       ///< outside the domain (or in a hole)
    Inside = 1,        ///< strictly inside an inside triangle, or on an unconstrained edge
    OnConstraint = 2,  ///< on a constrained edge, away from its ends
    OnVertex = 3,      ///< exactly on a vertex
};

struct Cdt2dWhere {
    Cdt2dPlace mPlace = Cdt2dPlace::Outside;
    std::int64_t mVertex = -1;  ///< the vertex, for `OnVertex`
};

class Cdt2d {
public:
    using Vec2 = std::array<double, 2>;
    using Tri = std::array<std::int64_t, 3>;

    /// The first vertex id of an inserted point (0, 1 and 2 are the super
    /// triangle's).
    static constexpr std::int64_t kFirstVertex = 3;

    /// A triangulation whose super triangle encloses the box, generously.
    Cdt2d(double MinU, double MinV, double MaxU, double MaxV);

    /// Insert one point; returns its vertex id, or the id of the vertex it
    /// coincides with exactly. The point must lie in the box given at
    /// construction (or close to it).
    std::int64_t Insert(const Vec2& rP);

    /// Insert many, in Hilbert order; the ids come back in input order.
    std::vector<std::int64_t> InsertBatch(const std::vector<Vec2>& rPoints);

    /// Make the segment between two vertices an edge and constrain it.
    Cdt2dStatus InsertSegment(std::int64_t A, std::int64_t B);

    /// Mark the inside triangles by the even-odd rule over the constraints.
    Cdt2dStatus Classify();

    /// After `Classify`: insert a point strictly inside the domain, or
    /// return -1 (outside, on a constrained edge or on a vertex, or outside
    /// `cdt2d_in_range`).
    std::int64_t InsertSteiner(const Vec2& rP);

    /// After `Classify`: where a point lies, decided exactly; inserts nothing.
    Cdt2dWhere Where(const Vec2& rP);

    /// The inside triangles (after `Classify`), counter-clockwise, in
    /// triangle-slot order.
    std::vector<Tri> InsideTriangles() const;

    /// Every live triangle, the super triangle's vertices included.
    std::vector<Tri> AllTriangles() const;

    /// Whether the edge between two vertices exists and is constrained.
    bool IsConstrainedEdge(std::int64_t A, std::int64_t B) const;

    const Vec2& Point(std::int64_t V) const { return mPoints[static_cast<std::size_t>(V)]; }
    std::int64_t NumVertices() const { return static_cast<std::int64_t>(mPoints.size()); }

private:
    struct Triangle {
        std::array<std::int64_t, 3> mV{};  ///< counter-clockwise
        std::array<std::int64_t, 3> mN{};  ///< the neighbour across the edge opposite mV[i]
        std::array<bool, 3> mC{};          ///< whether that edge is constrained
        bool mInside = false;
        bool mAlive = true;
    };
    struct Location {
        std::int64_t mTri = -1;
        int mKind = 0;  ///< 0 strictly inside, 1 on edge mIndex, 2 on vertex mVertex
        int mIndex = -1;
        std::int64_t mVertex = -1;
    };
    struct Boundary {
        std::int64_t mFrom = -1;
        std::int64_t mTo = -1;
        std::int64_t mOut = -1;
        bool mConstrained = false;
    };

    std::vector<Vec2> mPoints;
    std::vector<Triangle> mTris;
    std::vector<std::int64_t> mVertTri;
    std::vector<std::int64_t> mFree;
    std::int64_t mLast = 0;
    std::uint64_t mWalkState = 0x9E3779B97F4A7C15ull;
    bool mClassified = false;

    Location Locate(const Vec2& rP);
    Location LocateBruteForce(const Vec2& rP) const;
    std::int64_t InsertAt(const Vec2& rP, const Location& rLoc, bool Inside);
    std::vector<std::int64_t> ReplaceCavity(const std::vector<std::int64_t>& rRemoved,
                                            const std::vector<Tri>& rNew, bool Inside);
    void Legalize(std::vector<std::pair<std::int64_t, int>>& rStack);
    void Flip(std::int64_t T, int E);
    int IndexOf(std::int64_t T, std::int64_t V) const;
    int NeighbourIndex(std::int64_t T, std::int64_t Of) const;
    void SetConstrained(std::int64_t T, int E);
    void TriangulatePseudoPolygon(std::int64_t A, std::int64_t B,
                                  const std::vector<std::int64_t>& rChain,
                                  std::vector<Tri>& rOut) const;
    double Orient(std::int64_t A, std::int64_t B, std::int64_t C) const;
    double Orient(std::int64_t A, std::int64_t B, const Vec2& rC) const;
    double InCircle(std::int64_t A, std::int64_t B, std::int64_t C, std::int64_t D) const;
};

/**
 * @brief The pairs of segments (indices into `rSegments`, the first smaller,
 * ascending) that meet anywhere but at an endpoint they share: a proper
 * crossing, an endpoint on the other's interior, a collinear overlap, or two
 * segments sharing an endpoint and folding back along each other. Two
 * distinct points at the same place touch. A segment of zero length meets
 * nothing. Exact (`orient2d` and coordinate comparisons only); the
 * candidates come from a uniform bucket grid over the segments' box.
 */
std::vector<std::array<std::int64_t, 2>> cdt2d_find_crossings(
    const std::vector<Cdt2d::Vec2>& rPoints,
    const std::vector<std::array<std::int64_t, 2>>& rSegments);

/// A polygon with holes, triangulated: `rLoops[0]` the outer boundary, the
/// rest holes (any orientation; the even-odd rule decides the domain). Each
/// loop is a closed polyline whose last point is not repeated.
struct Cdt2dPolygon {
    Cdt2dStatus mStatus = Cdt2dStatus::Ok;
    std::string mWhy;
    std::vector<Cdt2d::Vec2> mPoints;    ///< boundary points first, then interior ones
    std::vector<Cdt2d::Tri> mTriangles;  ///< counter-clockwise, into mPoints
    std::vector<std::array<std::int64_t, 2>> mSegments;  ///< the constrained boundary edges
    std::int64_t mNumBoundary = 0;  ///< how many of mPoints are on the boundary
};

/**
 * @brief Triangulate a polygon with holes.
 * @param TargetSize 0 keeps the boundary as given and adds no interior point;
 *        `h > 0` splits every boundary segment longer than `h` into equal
 *        parts and fills the inside with an equilateral lattice of spacing
 *        `h`, keeping only points farther than `h / 2` from the boundary
 *        points and from every point already kept, and than `0.4 h` from
 *        the boundary segments.
 * @param MergeTolerance points closer than this are merged into the first
 *        (negative: `1e-12` times the box diagonal).
 */
Cdt2dPolygon cdt2d_triangulate_polygon(const std::vector<std::vector<Cdt2d::Vec2>>& rLoops,
                                       double TargetSize, double MergeTolerance = -1.0);

}  // namespace detail
}  // namespace meshioplusplus

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
 *  - **Refinement.** `Refine` is Ruppert's algorithm (J. Ruppert, "A
 *    Delaunay refinement algorithm for quality 2-dimensional mesh
 *    generation", J. Algorithms 18(3), 1995) with Chew's circumcentre
 *    insertion (L. P. Chew, "Guaranteed-quality triangular meshes", 1989;
 *    "Guaranteed-quality mesh generation for curved surfaces", SoCG 1993):
 *    a triangle whose circumradius-to-shortest-edge ratio exceeds
 *    `B = 1 / (2 sin(min_angle))` gets its circumcentre, unless that point
 *    encroaches a subsegment (lies strictly inside its diametral circle,
 *    the exact `diametral` predicate) or lies beyond one, which is split
 *    instead; a subsegment encroached by a vertex is split first. With
 *    input angles of at least 60 degrees this bounds every angle by about
 *    20.7 degrees at `B = sqrt(2)` (Ruppert 1995; J. R. Shewchuk, "Delaunay
 *    refinement algorithms for triangular mesh generation", Comput. Geom.
 *    22, 2002); past that `min_angle` is a target with a reported shortfall
 *    (it usually succeeds to about 33 degrees). Bad triangles are taken worst first
 *    (Shewchuk 2002), from a priority queue keyed by the ratio and broken by
 *    the sorted vertex ids, never by a slot or a counter, so the order is
 *    the input's alone. Input angles below 60 degrees (Ruppert's algorithm
 *    need not terminate there) are handled as in Shewchuk (2002) and G. L.
 *    Miller, S. E. Pav and N. J. Walkington, "When and why Ruppert's
 *    algorithm works" (IMR 2003): a subsegment at the apex of such an angle
 *    splits on concentric shells (at a power-of-two distance from the apex,
 *    so the subsegments of a cluster stay equal), a skinny triangle whose
 *    shortest edge spans the angle (or that holds it) is left as the angle
 *    forces it, and a subsegment of the cluster is not split when its
 *    halves would be shorter than the skinny triangle's shortest edge. The
 *    smallest angle those forced triangles keep is measured and reported
 *    per corner (`Cdt2dSmallAngle::mMinForcedDeg`). In place of shells,
 *    `mLopCorners` is corner lopping (S. E. Pav and N. J. Walkington,
 *    "Delaunay refinement by corner lopping", 14th IMR, 2005, Algorithm 2):
 *    an input point at an acute corner gets a ball of radius
 *    d = (sqrt 2 - 1) lfs0 (lfs0 its distance to the nearest other input
 *    point), its segments split where they leave it; in priority, a vertex
 *    inside a ball shrinks it by beta in [1/2, 1) and splits them again, an
 *    encroached subsegment is split (one at a corner by shrinking its ball),
 *    and a bad triangle is refined unless its point lies in a ball. It
 *    terminates for angles under `cdt2d_lop_bound_deg` (Corollary 3), and
 *    every triangle left under the bound has its vertices within 2 d of a
 *    corner (Lemma 7), checked and reported per corner. Every
 *    insertion counts against `max_points`. A caller that shares its
 *    boundary with others (`mesh_cad_surface`'s edges) asks for the splits
 *    instead of making them, and may judge triangles by its own measure
 *    (the surface's); the plane stays where the points are found.
 *  - **A metric.** With `SetMetric`, the flips after an insertion decide by
 *    the in-circumellipse test of H. Borouchaki, P. L. George, F. Hecht, P.
 *    Laug and E. Saltel, "Delaunay mesh generation governed by metric
 *    specifications. Part I. Algorithms", Finite Elem. Anal. Des. 25, 1997:
 *    the opposite vertex lies in the triangle's circumcircle measured in the
 *    metric at the shared edge's midpoint (the four points mapped by its
 *    Cholesky factor, then the exact `incircle` on the mapped coordinates).
 *    A flip still needs a strictly convex quadrilateral (exact `orient2d`),
 *    and, as without a metric, only the edges opposite the new vertex are
 *    tested, each flip adding an edge to its star that is never flipped
 *    away, so the flips end. Refinement's encroachment is measured in the
 *    same metric (the diametral circle at the subsegment's midpoint). The
 *    triangulation is no longer Delaunay in the plane; insertion and
 *    refinement need only a valid one.
 *  - **Moving a vertex.** `MoveVertex` moves an interior vertex inside its
 *    star, refusing a move that would invert a triangle (exact `orient2d`),
 *    then restores the Delaunay property by flips (smoothing).
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
#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>
#include <utility>
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

/// What `Cdt2d::Refine` is asked for.
struct Cdt2dRefineOptions {
    /// The angle sought, in degrees (0, 60): a triangle is bad when its
    /// circumradius-to-shortest-edge ratio exceeds 1 / (2 sin(min_angle)).
    double mMinAngleDeg = 20.7;
    /// The most vertices the pass may insert (splits and circumcentres).
    std::int64_t mMaxPoints = 1000000;
    /// Split an encroached subsegment here (a planar domain), or only ask
    /// for it (`Cdt2dRefineReport::mRequests`), when the caller owns the
    /// boundary (a curved edge shared with other faces).
    bool mSplitSegments = true;
    /// The ratio to judge a triangle (its vertex ids) by; null: in the
    /// plane. A negative value exempts the triangle (counted as forced).
    std::function<double(const std::array<std::int64_t, 3>&)> mRatio;
    /// The point to insert for a bad triangle; null: its circumcentre in
    /// the plane. (A caller judging triangles by another measure gives the
    /// plane image of that measure's circumcentre.)
    std::function<std::array<double, 2>(const std::array<std::int64_t, 3>&)> mCentre;
    /// Whether to insert the point found for a bad triangle, given the
    /// vertices of the cavity it would open; null: always. A declined
    /// triangle is left (counted unresolved).
    std::function<bool(const std::array<std::int64_t, 3>&, const std::array<double, 2>&,
                       const std::vector<std::int64_t>&)>
        mAccept;
    /// Whether to split (or ask for) the subsegment `A`-`B` a bad triangle's
    /// point encroaches; null: always. A declined triangle is left
    /// (counted unresolved).
    std::function<bool(std::int64_t, std::int64_t, const std::array<std::int64_t, 3>&)>
        mAcceptSplit;
    /// The smallest angle (degrees) to report a forced triangle by; null:
    /// in the plane.
    std::function<double(const std::array<std::int64_t, 3>&)> mMinAngle;
    /// Corner lopping (Pav and Walkington 2005, Algorithm 2) in place of
    /// concentric shells: every input point at an acute corner (two input
    /// segments meeting at under 90 degrees inside the domain) is protected by
    /// a ball of radius (sqrt 2 - 1) lfs0, its segments split where they leave
    /// it; a vertex inside a ball shrinks it by `mLopBeta` and splits them
    /// again; a bad triangle whose circumcentre lies in a ball is left.
    bool mLopCorners = false;
    /// The factor a ball shrinks by, in [1/2, 1).
    double mLopBeta = 0.5;
    /// The input points (vertex ids) beside every end of a constrained edge
    /// at an acute corner at the start (a kink inside a curve is one); empty:
    /// every end of a constrained edge at the start.
    std::vector<std::int64_t> mCorners;
    /// Vertices never input points (a pole's, which has its own handling).
    std::vector<std::int64_t> mLopSkip;
    /// A ball's radius carried from an earlier pass (vertex id, radius), used
    /// in place of (sqrt 2 - 1) lfs0, so balls only shrink.
    std::vector<std::pair<std::int64_t, double>> mLopRadius;
    /// The distance from input point `P` to vertex `V`, or, with `V` < 0, to
    /// the plane point `rX`, by which balls and lfs0 are measured; null: in
    /// the plane.
    std::function<double(std::int64_t P, std::int64_t V, const std::array<double, 2>& rX)>
        mDistance;
};

/// An acute input corner lopped (`Cdt2dRefineOptions::mLopCorners`).
struct Cdt2dLoppedCorner {
    std::int64_t mVertex = -1;
    double mAngleDeg = 0.0;  ///< the smallest acute wedge at it, in the plane
    double mLfs0 = 0.0;      ///< the distance to the nearest other input point
    double mRadius = 0.0;    ///< its ball's radius d at the end
    /// The bad triangles left because their circumcentre lies in the ball,
    /// and the smallest angle among them (180 for none).
    std::int64_t mNumForced = 0;
    double mMinForcedDeg = 180.0;
    /// Whether every one of them lies within 2 d of the corner (Lemma 7).
    bool mWithinTwoRadii = true;
};

/// A corner where two segments meet at less than 60 degrees inside the
/// domain (in the plane): Ruppert's small input angle.
struct Cdt2dSmallAngle {
    std::int64_t mVertex = -1;
    double mAngleDeg = 0.0;  ///< the smallest wedge at it
    /// The smallest angle of the bad triangles it forces (holding it, or
    /// whose shortest edge spans it) after refinement; 180 for none.
    double mMinForcedDeg = 180.0;
};

/// A subsegment the caller is asked to split (`mSplitSegments` false), its
/// vertices in ascending order; `mApex` the end at a small input angle that
/// its split should keep a power-of-two distance from, or -1.
struct Cdt2dSplitRequest {
    std::int64_t mA = -1;
    std::int64_t mB = -1;
    std::int64_t mApex = -1;
    /// Corner lopping: split where the distance from `mApex` is this (the
    /// ball's radius); 0: the subsegment's midpoint, or its shell.
    double mDistance = 0.0;
    /// By subsegment, then apex: corner lopping may ask a subsegment
    /// between two lopped corners to split at each ball.
    bool operator<(const Cdt2dSplitRequest& rO) const {
        if (mA != rO.mA)
            return mA < rO.mA;
        return mB != rO.mB ? mB < rO.mB : mApex < rO.mApex;
    }
    bool operator==(const Cdt2dSplitRequest& rO) const {
        return mA == rO.mA && mB == rO.mB && mApex == rO.mApex;
    }
};

struct Cdt2dRefineReport {
    std::int64_t mNumCircumcentres = 0;         ///< circumcentres inserted
    std::int64_t mNumSplits = 0;                ///< subsegments split here
    std::vector<Cdt2dSplitRequest> mRequests;   ///< sorted, distinct
    std::vector<Cdt2dSmallAngle> mSmallAngles;  ///< by vertex
    /// Bad triangles left because a small input angle forces them (or the
    /// caller exempted them), and those left otherwise (the budget, a
    /// split asked for, a circumcentre on a vertex).
    std::int64_t mNumForced = 0;
    std::int64_t mNumUnresolved = 0;
    bool mBudgetReached = false;
    /// The bad triangles a small input angle forces (vertex ids ascending),
    /// with that angle's apex, sorted.
    std::vector<std::pair<std::array<std::int64_t, 3>, std::int64_t>> mForced;
    /// Corner lopping: the acute corners, by vertex.
    std::vector<Cdt2dLoppedCorner> mLopped;
};

/**
 * @brief The largest angle bound (degrees) Pav and Walkington (2005,
 * Corollary 3) prove corner lopping terminates for, given `DeltaT`, the most
 * any (sub)curve turns: kappa < arcsin(1 / (2 sqrt 2 mu)) and
 * arcsin(1 / (2 zeta eta)), with eta = (1 + tan(dt/2)) / (1 - tan(dt/2)),
 * mu = 1 / (1 - tan(dt/2)), zeta = sqrt 2 sin(dt/2) / (sqrt(1 + sin dt) - 1).
 * arcsin(1 / 2 sqrt 2), 20.7 degrees, for straight segments; 0 past pi/3,
 * where the paper assumes no curve goes.
 */
double cdt2d_lop_bound_deg(double DeltaT);

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

    /// After `Classify`: Delaunay refinement of the inside (see the file
    /// comment). New vertices take ids from `NumVertices()` on.
    Cdt2dRefineReport Refine(const Cdt2dRefineOptions& rOptions);

    /// After `Classify`: move vertex `V` (not on a constrained edge) to
    /// `rP`; refused (false, nothing changed) when `rP` is outside the star
    /// of `V`'s triangles or would invert one of them. Lawson's flips
    /// follow, unless `Legalize` is false (the caller flips by its own
    /// measure, `FlipEdge`).
    bool MoveVertex(std::int64_t V, const Vec2& rP, bool Legalize = true);

    /// Every constrained edge, its vertices ascending, sorted.
    std::vector<std::array<std::int64_t, 2>> ConstrainedEdges() const;

    /// Flip the unconstrained edge between `A` and `B` when the two
    /// triangles beside it form a strictly convex quadrilateral (exact
    /// orient2d); false (nothing changed) otherwise. No flips follow: the
    /// caller decides by its own measure (a surface's angles).
    bool FlipEdge(std::int64_t A, std::int64_t B);

    /// Whether `V` is an endpoint of a constrained edge.
    bool OnConstraint(std::int64_t V) const;

    /// The vertices joined to `V` by an edge, counter-clockwise.
    std::vector<std::int64_t> Neighbours(std::int64_t V) const;

    /// The triangles around `V`, counter-clockwise.
    std::vector<Tri> Star(std::int64_t V) const;

    /// The inside triangles (after `Classify`) around `V`.
    std::vector<Tri> InsideStar(std::int64_t V) const;

    /// The inside triangles (after `Classify`), counter-clockwise, in
    /// triangle-slot order.
    std::vector<Tri> InsideTriangles() const;

    /// Every live triangle, the super triangle's vertices included.
    std::vector<Tri> AllTriangles() const;

    /// Whether the edge between two vertices exists and is constrained.
    bool IsConstrainedEdge(std::int64_t A, std::int64_t B) const;

    /// A metric over the plane, (m11, m12, m22) at a point, symmetric
    /// positive definite; one that is not (or not finite) where it is asked
    /// falls back to the plane's in-circle test there.
    using Metric = std::function<std::array<double, 3>(const Vec2&)>;

    /// From now on, decide the flips after each insertion in this metric
    /// (see the file comment); an empty one restores the plane's test.
    void SetMetric(Metric M) { mMetric = std::move(M); }

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
    Metric mMetric;

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
    /// Whether Lawson's flip of the edge (B, C) of triangle (A, B, C), D
    /// across it, is due: the plane's in-circle test, or the metric's.
    bool FlipDue(std::int64_t A, std::int64_t B, std::int64_t C, std::int64_t D) const;
    /// The first `N` points mapped by the metric's Cholesky factor at
    /// `rAt`, relative to it; false (unchanged) without a metric, where it
    /// is not positive definite, or where a mapped coordinate leaves
    /// `cdt2d_in_range`.
    bool MetricMapped(const Vec2& rAt, std::array<Vec2, 4>& rPoints, std::size_t N = 4) const;
    std::vector<std::int64_t> TrianglesAround(std::int64_t V) const;
    bool FindEdge(std::int64_t A, std::int64_t B, std::int64_t& rT, int& rE) const;
    std::int64_t SplitConstrained(std::int64_t T, int E, const Vec2& rP);

    friend class Cdt2dRefiner;
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
    /// The constrained boundary edges (after refinement, its subsegments).
    std::vector<std::array<std::int64_t, 2>> mSegments;
    std::int64_t mNumBoundary = 0;  ///< how many of mPoints are on the given boundary
    Cdt2dRefineReport mRefine;      ///< what refinement did, when asked for
};

/// Quality asked of `cdt2d_triangulate_polygon`: 0 keeps today's output.
struct Cdt2dPolygonQuality {
    double mMinAngleDeg = 0.0;
    std::int64_t mMaxPoints = 1000000;  ///< the most points refinement adds
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

/// The same, then refined (`Cdt2d::Refine`, splitting the boundary here)
/// when `rQuality.mMinAngleDeg` > 0; the points it adds follow the lattice's.
Cdt2dPolygon cdt2d_triangulate_polygon(const std::vector<std::vector<Cdt2d::Vec2>>& rLoops,
                                       double TargetSize, double MergeTolerance,
                                       const Cdt2dPolygonQuality& rQuality);

}  // namespace detail
}  // namespace meshioplusplus

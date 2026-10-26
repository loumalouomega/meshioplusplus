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
// mesh_cad_surface and triangulate_polygon_2d: see operations/mesh_cad.hpp
// for the contract and doc/mesh_cad.md for the algorithm. Anonymous-namespace
// helpers are prefixed `mcad_` (the amalgamation concatenates every
// translation unit into one).

// System includes
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "../detail/cdt2d.hpp"
#include "../detail/curved_triangle.hpp"
#include "../detail/edge_sampler.hpp"
#include "../detail/json_value.hpp"
#include "../detail/mesh_cad_reparam.hpp"
#include "../detail/mesh_cad_size.hpp"
#include "mesh_cad_testing.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/operations/mesh_cad.hpp"
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/region.hpp"
#include "meshioplusplus/detail/value_io.hpp"

namespace meshioplusplus {

namespace {

constexpr const char* kMcadPrefix = "meshio++: mesh_cad_surface: ";
constexpr double kMcadPi = 3.14159265358979323846;

std::atomic<std::uint64_t> g_mcad_face_order_seed{0};

/// The order the faces' parallel loops visit them in: their index order,
/// or (the test hook) a permutation of it.
std::vector<std::size_t> mcad_face_order(std::size_t N) {
    std::vector<std::size_t> order(N);
    std::iota(order.begin(), order.end(), std::size_t{0});
    std::uint64_t x = g_mcad_face_order_seed.load();
    if (x == 0)
        return order;
    for (std::size_t i = N; i > 1; --i) {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        std::swap(order[i - 1], order[static_cast<std::size_t>(x % i)]);
    }
    return order;
}

using Vec2 = detail::EdgeVec2;
using Vec3 = detail::EdgeVec3;

[[noreturn]] void mcad_fail(const std::string& rWhat) {
    throw std::invalid_argument(std::string(kMcadPrefix) + rWhat);
}

/// An exception's text without the operation's prefix.
std::string mcad_why(const std::exception& rE) {
    std::string s = rE.what();
    const std::string p = kMcadPrefix;
    if (s.rfind(p, 0) == 0)
        s.erase(0, p.size());
    return s;
}

double mcad_dist(const Vec3& rA, const Vec3& rB) {
    return std::sqrt((rA[0] - rB[0]) * (rA[0] - rB[0]) + (rA[1] - rB[1]) * (rA[1] - rB[1]) +
                     (rA[2] - rB[2]) * (rA[2] - rB[2]));
}

double mcad_segment_distance(const Vec3& rP, const Vec3& rA, const Vec3& rB) {
    const Vec3 d{rB[0] - rA[0], rB[1] - rA[1], rB[2] - rA[2]};
    const double len2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
    double t = 0.0;
    if (len2 > 0.0)
        t = std::clamp(
            ((rP[0] - rA[0]) * d[0] + (rP[1] - rA[1]) * d[1] + (rP[2] - rA[2]) * d[2]) / len2, 0.0,
            1.0);
    return mcad_dist(rP, {rA[0] + t * d[0], rA[1] + t * d[1], rA[2] + t * d[2]});
}

/// A straight trim in the parameters, for the boundary of an untrimmed face.
BrepCurve mcad_line(double U0, double V0, double U1, double V1) {
    BrepCurve c;
    c.mDegree = 1;
    c.mKnots = {0.0, 0.0, 1.0, 1.0};
    c.mControlPoints = NDArray(DType::Float64, {2, 3});
    double* p = c.mControlPoints.As<double>();
    p[0] = U0;
    p[1] = V0;
    p[2] = 0.0;
    p[3] = U1;
    p[4] = V1;
    p[5] = 0.0;
    c.mWeights = NDArray(DType::Float64, {2});
    c.mWeights.As<double>()[0] = 1.0;
    c.mWeights.As<double>()[1] = 1.0;
    return c;
}

/// The order to walk a loop's trims in: head to tail by their ends in the
/// face's parameters. The listed order is kept wherever it connects (to a
/// billionth of the loop's extent); elsewhere the next trim is the one whose
/// start is nearest the current end, so a loop listed out of order (a hole
/// whose trims run backwards but are listed forwards) still closes. A trim
/// that cannot be evaluated keeps the listed order; the trace fails later.
std::vector<std::size_t> mcad_chain(const BrepLoop& rLoop) {
    const std::size_t n = rLoop.mTrims.size();
    std::vector<std::size_t> order(n);
    std::iota(order.begin(), order.end(), std::size_t{0});
    if (n < 3)
        return order;
    std::vector<Vec2> start(n);
    std::vector<Vec2> end(n);
    double lo[2] = {std::numeric_limits<double>::infinity(),
                    std::numeric_limits<double>::infinity()};
    double hi[2] = {-lo[0], -lo[1]};
    try {
        for (std::size_t k = 0; k < n; ++k) {
            const detail::BrepCurveEval c(rLoop.mTrims[k].mCurve);
            double a[3];
            double b[3];
            c.Eval(c.T0(), a, nullptr);
            c.Eval(c.T1(), b, nullptr);
            start[k] = {a[0], a[1]};
            end[k] = {b[0], b[1]};
            if (!rLoop.mTrims[k].mDirection)
                std::swap(start[k], end[k]);
            for (const Vec2* p : {&start[k], &end[k]})
                for (int d = 0; d < 2; ++d) {
                    lo[d] = std::min(lo[d], (*p)[static_cast<std::size_t>(d)]);
                    hi[d] = std::max(hi[d], (*p)[static_cast<std::size_t>(d)]);
                }
        }
    } catch (const std::exception&) {
        return order;
    }
    const double tol = 1e-9 * std::hypot(hi[0] - lo[0], hi[1] - lo[1]);
    const auto gap = [&](std::size_t A, std::size_t B) {
        return std::hypot(end[A][0] - start[B][0], end[A][1] - start[B][1]);
    };
    std::vector<bool> used(n, false);
    order.clear();
    order.push_back(0);
    used[0] = true;
    for (std::size_t step = 1; step < n; ++step) {
        const std::size_t cur = order.back();
        std::size_t next = cur + 1 < n && !used[cur + 1] && gap(cur, cur + 1) <= tol ? cur + 1 : n;
        if (next == n) {
            double best = std::numeric_limits<double>::infinity();
            for (std::size_t k = 0; k < n; ++k)
                if (!used[k] && gap(cur, k) < best) {
                    best = gap(cur, k);
                    next = k;
                }
        }
        used[next] = true;
        order.push_back(next);
    }
    return order;
}

class McadUnionFind {
public:
    explicit McadUnionFind(std::size_t N) : mParent(N) {
        std::iota(mParent.begin(), mParent.end(), 0);
    }
    std::size_t Find(std::size_t X) {
        while (mParent[X] != X) {
            mParent[X] = mParent[mParent[X]];
            X = mParent[X];
        }
        return X;
    }
    // The smaller root wins, so the classes do not depend on the union order.
    void Join(std::size_t A, std::size_t B) {
        A = Find(A);
        B = Find(B);
        if (A != B)
            mParent[std::max(A, B)] = std::min(A, B);
    }

private:
    std::vector<std::size_t> mParent;
};

/// A uniform hash grid of 3-D points and segments, for the spacing test.
class McadGrid {
public:
    explicit McadGrid(double Cell) : mCell(Cell) {}
    using Key = std::array<std::int64_t, 3>;
    Key KeyOf(const Vec3& rP) const {
        return {static_cast<std::int64_t>(std::floor(rP[0] / mCell)),
                static_cast<std::int64_t>(std::floor(rP[1] / mCell)),
                static_cast<std::int64_t>(std::floor(rP[2] / mCell))};
    }
    void AddPoint(const Vec3& rP) { mPoints[KeyOf(rP)].push_back(rP); }
    void AddSegment(const Vec3& rA, const Vec3& rB) {
        const Vec3 m{0.5 * (rA[0] + rB[0]), 0.5 * (rA[1] + rB[1]), 0.5 * (rA[2] + rB[2])};
        mSegments[KeyOf(m)].push_back({rA, rB});
        mLongest = std::max(mLongest, mcad_dist(rA, rB));
    }
    /// Whether a point lies within `Rp` of a point or `Rs` of a segment.
    bool Near(const Vec3& rP, double Rp, double Rs) const {
        const Key k = KeyOf(rP);
        const auto reach = [&](double R) {
            return static_cast<std::int64_t>(std::ceil(R / mCell));
        };
        const std::int64_t rp = reach(Rp);
        for (std::int64_t i = -rp; i <= rp; ++i)
            for (std::int64_t j = -rp; j <= rp; ++j)
                for (std::int64_t l = -rp; l <= rp; ++l) {
                    const auto it = mPoints.find({k[0] + i, k[1] + j, k[2] + l});
                    if (it == mPoints.end())
                        continue;
                    for (const Vec3& r_q : it->second)
                        if (mcad_dist(r_q, rP) < Rp)
                            return true;
                }
        // A segment within Rs has its midpoint within Rs + half its length.
        const std::int64_t rs = reach(Rs + 0.5 * mLongest);
        for (std::int64_t i = -rs; i <= rs; ++i)
            for (std::int64_t j = -rs; j <= rs; ++j)
                for (std::int64_t l = -rs; l <= rs; ++l) {
                    const auto it = mSegments.find({k[0] + i, k[1] + j, k[2] + l});
                    if (it == mSegments.end())
                        continue;
                    for (const auto& r_s : it->second)
                        if (mcad_segment_distance(rP, r_s[0], r_s[1]) < Rs)
                            return true;
                }
        return false;
    }

private:
    struct KeyHash {
        std::size_t operator()(const Key& rK) const {
            std::uint64_t h = 1469598103934665603ull;
            for (const std::int64_t v : rK) {
                h ^= static_cast<std::uint64_t>(v);
                h *= 1099511628211ull;
            }
            return static_cast<std::size_t>(h);
        }
    };
    double mCell;
    double mLongest = 0.0;
    std::unordered_map<Key, std::vector<Vec3>, KeyHash> mPoints;
    std::unordered_map<Key, std::vector<std::array<Vec3, 2>>, KeyHash> mSegments;
};

/// One trim of a face's loop: a use of an edge.
struct McadUse {
    std::size_t mFace = 0;
    std::int64_t mTrimIndex = -1;
    const BrepCurve* pCurve = nullptr;
    bool mForward = true;
    std::size_t mLoopSize = 1;  ///< trims in its loop
    bool mEmbedded = false;     ///< a curve inside the face, not on a loop
    std::int64_t mEdge = -1;
    bool mSame = true;  ///< runs the way its edge's master does
    /// A closed use rotated onto its master's corner: the arc length along
    /// its trace where it now starts and ends; -1 for none.
    double mStartLength = -1.0;
    std::vector<double> mParams;     ///< along this use, ends included
    std::vector<std::int64_t> mIds;  ///< mesh points along this use, ends included
    MeshCadFaceStatus mFail = MeshCadFaceStatus::Ok;
    std::string mWhy;
};

struct McadEdge {
    std::int64_t mBrepEdge = -1;
    std::vector<std::size_t> mUses;  ///< ascending
    /// The use sampled, whose points the others share: the first, unless a
    /// closed edge's uses start apart and a later one cannot be rotated.
    std::size_t mMaster = 0;
    detail::EdgeSamples mSamples;
    std::int64_t mFirstInterior = 0;
    MeshCadFaceStatus mFail = MeshCadFaceStatus::Ok;
    std::string mWhy;
};

/// A boundary piece corner lopping asks to split (6.2.3) where its distance
/// from the lopped corner at one end (1 the piece's start, 2 its end) is the
/// corner's ball radius.
struct McadLopAsk {
    std::size_t mEdge = 0;
    std::size_t mPiece = 0;
    std::size_t mEnd = 0;
    double mDistance = 0.0;
};

struct McadFace {
    bool mUsable = false;
    std::unique_ptr<detail::BrepSurfaceEval> pSurface;
    std::vector<BrepCurve> mSynthetic;
    std::vector<std::vector<std::size_t>> mLoops;  ///< use indices per loop
    std::vector<std::size_t> mEmbedded;            ///< use indices of its embedded edges
    /// Embedded edges left out, with why, warned in face order.
    std::vector<std::pair<std::size_t, std::string>> mSkipped;
    std::size_t mSkippedBefore = std::string::npos;  ///< mSkipped's size before meshing
    MeshCadFaceRecord mRecord;
    std::vector<Vec3> mInterior;
    std::vector<std::array<std::int64_t, 3>> mTriangles;  ///< < 0: -(1 + interior index)
    /// Boundary pieces refinement asks to split: (edge, master piece, the
    /// end at a sharp corner: 0 none, 1 the piece's start, 2 its end).
    std::vector<std::array<std::size_t, 3>> mAsks;
    /// Corner lopping's asks (6.2.3), and the corners the last meshing
    /// lopped, by position, with what is reported of each.
    std::vector<McadLopAsk> mLopAsks;
    std::vector<std::pair<Vec3, MeshCadLoppedCorner>> mLopped;
    std::array<std::int64_t, kMeshCadAngleBins> mHistogram{};
    /// Triangles under min_angle that refinement leaves where the face's
    /// metric departs too far from its plane's, warned in face order.
    std::int64_t mNumStretched = 0;
    /// Refinement spent the face's share of max_points; the face is then
    /// meshed unrefined (`mUnrefined`) and warned.
    bool mRefineBudgetReached = false;
    bool mUnrefined = false;
    /// Points refinement inserted in the face over every round so far: each
    /// round meshes the face again, so the bound is on their sum.
    std::int64_t mRefinedSoFar = 0;
    /// The face asked for boundary splits, so it is meshed again and its
    /// smoothing was left for that meshing; and `mRefinedSoFar` before it,
    /// so the same meshing can be made again, smoothed.
    bool mSmoothingDeferred = false;
    std::int64_t mRefinedBefore = 0;
    /// Curved elements (order 2 and 3, roadmap 6.6), per triangle of
    /// `mTriangles`: each edge k (corner k to corner k + 1) runs along a
    /// boundary piece (its nodes the edge's, made at assembly) or holds
    /// nodes of the face's own (`mHighOrder`, from corner k); the centre
    /// (order 3); and the face's normal there, the sign of the element's
    /// Jacobian.
    struct CurvedEdge {
        std::int64_t mEdge = -1;                     ///< boundary: the McadEdge
        std::int64_t mPiece = -1;                    ///< its master piece
        bool mReversed = false;                      ///< runs against the master
        std::array<std::int64_t, 2> mNodes{-1, -1};  ///< interior: into mHighOrder
    };
    std::vector<std::array<CurvedEdge, 3>> mCurvedEdges;
    std::vector<std::int64_t> mCurvedCentre;
    std::vector<Vec3> mCurvedNormal;
    std::vector<Vec3> mHighOrder;
    /// Pole rings (6.4.2) were placed in the face's latest meshing.
    bool mPoleRingsPlaced = false;
    /// In some meshing its unrefined triangulation folded on the surface
    /// (6.4.3): it asks for boundary splits even when no angle or size is
    /// sought, so the rounds run.
    bool mNeedsRounds = false;
    /// Every triangle of its latest meshing collapsed (6.5.7).
    bool mCollapsed = false;
    /// By boundary point: a corner virtual edges suppress (6.5.2), left out
    /// of the face's loops; null when there are none.
    const std::vector<char>* pSuppressed = nullptr;
    /// Mesh it with a small interior grid whatever its size: set after a
    /// meshing collapsed, kept for every later one.
    bool mForceInterior = false;
    /// Its latest meshing left a fold (6.4.3), and how many meshings in a
    /// row did.
    bool mFolded = false;
    int mFoldStreak = 0;
    /// Mesh it in its mean metric whatever the options: set after a
    /// meshing in the varying metric folded (6.4.1), kept for every later
    /// one.
    bool mMeanMetric = false;
    /// The first meshing of a sized run (6.3): each point of the face with
    /// the surface's normal and largest principal curvature there.
    std::vector<detail::MeshCadSizeSample> mSizeSamples;
    /// The discrete reparametrisation (6.4.4): measured on the face's first
    /// meshing, then kept; a face that takes it and later cannot (a
    /// flattening that folds once its boundary is split) is declined for
    /// every later meshing.
    enum class Route : std::uint8_t { Unmeasured, Below, Taken, Declined };
    Route mRoute = Route::Unmeasured;
    double mRouteCond = 0.0;
    double mRouteCondDisc = 0.0;
    std::int32_t mRouteMethod = 0;
    std::string mRouteWhy;
    /// Per loop and use: walked against the use (a merged sliver's, 6.5.2);
    /// empty for every use forward.
    std::vector<std::vector<char>> mLoopReversed;
    /// Virtual faces (6.5.2): on a host, the sliver meshed with it and their
    /// merged loop (the host's with the shared run replaced by the sliver's
    /// other uses); on a sliver, its host. -1 for none.
    std::int64_t mMergedWith = -1;
    std::int64_t mMergedInto = -1;
    std::vector<std::vector<std::size_t>> mMergedLoops;
    std::vector<std::vector<char>> mMergedReversed;
    /// The `BrepEdge`s a host and its sliver share, left out of the mesh.
    std::vector<std::size_t> mSharedEdges;
    /// Why a merge found was undone (on the host it was found with).
    std::string mMergeWhy;
    /// A merged meshing's triangles: 0 the host's, 1 the sliver's (by the
    /// coarse triangle holding each centroid).
    std::vector<std::int32_t> mTriOwner;
};

/// What refinement and smoothing ask of every face.
struct McadQuality {
    double mMinAngleDeg = 0.0;
    std::int64_t mSmoothIterations = 0;
    MeshCadSmoothing mSmoothMethod = MeshCadSmoothing::Laplacian;
    std::vector<char> mPole;     ///< by mesh point: a pole's corner
    std::int64_t mNumFaces = 1;  ///< faces meshed, sharing max_points
    /// A face that asks for boundary splits is not smoothed: the next round
    /// meshes it again (or, when the rounds stop, `mDeferSmoothing` off).
    bool mDeferSmoothing = false;
    /// The size field of a sized run (6.3); null: `target_size` everywhere.
    const detail::MeshCadSizeField* pSize = nullptr;
    /// Record each face's size samples (the first meshing of a sized run).
    bool mCollectSizes = false;
    /// Do not refine for `pSize` (a face meshed unrefined after spending its
    /// points); the lattice still follows it.
    bool mNoSizeRefine = false;
    /// Corner lopping (6.2.3) in place of concentric shells, and the balls'
    /// radii carried between meshings by corner position: the smallest any
    /// face's meshing or any split left.
    bool mLopCorners = false;
    const std::map<Vec3, double>* pLopRadius = nullptr;
};

/// The angles of a triangle on the surface, in degrees.
std::array<double, 3> mcad_angles(const Vec3& rA, const Vec3& rB, const Vec3& rC) {
    const Vec3* p[3] = {&rA, &rB, &rC};
    std::array<double, 3> out{};
    for (int k = 0; k < 3; ++k) {
        const Vec3& o = *p[k];
        const Vec3& q = *p[(k + 1) % 3];
        const Vec3& r = *p[(k + 2) % 3];
        const Vec3 a{q[0] - o[0], q[1] - o[1], q[2] - o[2]};
        const Vec3 b{r[0] - o[0], r[1] - o[1], r[2] - o[2]};
        const Vec3 c{a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
                     a[0] * b[1] - a[1] * b[0]};
        out[static_cast<std::size_t>(k)] =
            std::atan2(std::sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]),
                       a[0] * b[0] + a[1] * b[1] + a[2] * b[2]) *
            (180.0 / 3.14159265358979323846);
    }
    return out;
}

double mcad_min_angle(const Vec3& rA, const Vec3& rB, const Vec3& rC) {
    const auto a = mcad_angles(rA, rB, rC);
    return std::min({a[0], a[1], a[2]});
}

/// Circumradius over shortest edge of a triangle on the surface: the
/// product of its two longer edges over twice its area's double.
double mcad_ratio(const Vec3& rA, const Vec3& rB, const Vec3& rC) {
    const double la = mcad_dist(rB, rC);
    const double lb = mcad_dist(rC, rA);
    const double lc = mcad_dist(rA, rB);
    const Vec3 u{rB[0] - rA[0], rB[1] - rA[1], rB[2] - rA[2]};
    const Vec3 v{rC[0] - rA[0], rC[1] - rA[1], rC[2] - rA[2]};
    const Vec3 c{u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]};
    const double twice_area = std::sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]);
    if (!(twice_area > 0.0))
        return std::numeric_limits<double>::infinity();
    return la * lb * lc / (2.0 * twice_area * std::min({la, lb, lc}));
}

/// How far a trace turns between two parameters (radians): the angles
/// between its tangents at nine equal steps of the parameter, summed (Pav and
/// Walkington's total variation of a curve, 6.2.3).
double mcad_piece_turning(const detail::TrimTrace& rTrace, double T0, double T1) {
    double sum = 0.0;
    Vec3 last{};
    for (int j = 0; j <= 8; ++j) {
        const Vec3 v = rTrace.Velocity(T0 + (T1 - T0) * (j / 8.0));
        if (j > 0) {
            const Vec3 c{last[1] * v[2] - last[2] * v[1], last[2] * v[0] - last[0] * v[2],
                         last[0] * v[1] - last[1] * v[0]};
            sum += std::atan2(std::sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]),
                              last[0] * v[0] + last[1] * v[1] + last[2] * v[2]);
        }
        last = v;
    }
    return sum;
}

/// Curved elements whose Jacobian is not positive everywhere (6.6, 6.6.1)
/// made straight-sided: their edge nodes moved to the chord points of their
/// corners and the centre to the corners' centroid. The nodes are shared, so
/// the elements around a moved node are judged again, until none changes
/// (nodes only ever move to straight). `rElems` holds the node ids, `Order`
/// gives their count; returns the elements straightened, and in
/// `rStillInvalid` those invalid even straight.
std::int64_t mcad_straighten(int Order, const std::vector<std::array<std::int64_t, 10>>& rElems,
                             const std::vector<Vec3>& rNormal, double* pX,
                             std::int64_t& rStillInvalid) {
    const int nn = detail::curved_triangle_num_nodes(Order);
    // Elements around each high-order node.
    std::unordered_map<std::int64_t, std::vector<std::size_t>> around;
    for (std::size_t e = 0; e < rElems.size(); ++e)
        for (int i = 3; i < nn; ++i)
            around[rElems[e][static_cast<std::size_t>(i)]].push_back(e);
    std::vector<char> straight(rElems.size(), 0);
    std::set<std::size_t> work;
    for (std::size_t e = 0; e < rElems.size(); ++e)
        work.insert(e);
    std::int64_t count = 0;
    rStillInvalid = 0;
    std::array<detail::CurvedVec3, 10> nodes{};
    while (!work.empty()) {
        const std::size_t e = *work.begin();
        work.erase(work.begin());
        const auto& r_el = rElems[e];
        for (int i = 0; i < nn; ++i)
            for (std::size_t k = 0; k < 3; ++k)
                nodes[static_cast<std::size_t>(i)][k] =
                    pX[3 * static_cast<std::size_t>(r_el[static_cast<std::size_t>(i)]) + k];
        if (detail::curved_triangle_check(Order, nodes.data(), rNormal[e]).mValid)
            continue;
        if (straight[e] != 0) {
            ++rStillInvalid;
            continue;
        }
        straight[e] = 1;
        ++count;
        // Edge k from corner k to corner k + 1, its nodes from corner k.
        int slot = 3;
        for (std::size_t k = 0; k < 3; ++k) {
            const auto& a = nodes[k];
            const auto& b = nodes[(k + 1) % 3];
            for (int j = 1; j < Order; ++j, ++slot) {
                const double f = double(j) / Order;
                const auto id = static_cast<std::size_t>(r_el[static_cast<std::size_t>(slot)]);
                for (std::size_t c = 0; c < 3; ++c)
                    pX[3 * id + c] = a[c] + f * (b[c] - a[c]);
            }
        }
        if (Order == 3) {
            const auto id = static_cast<std::size_t>(r_el[9]);
            for (std::size_t c = 0; c < 3; ++c)
                pX[3 * id + c] = (nodes[0][c] + nodes[1][c] + nodes[2][c]) / 3.0;
        }
        for (int i = 3; i < nn; ++i)
            for (const std::size_t o : around[r_el[static_cast<std::size_t>(i)]])
                if (o != e)
                    work.insert(o);
    }
    return count;
}

void mcad_check_options(const MeshCadSurfaceOptions& rO) {
    if (!(std::isfinite(rO.mTargetSize) && rO.mTargetSize > 0.0))
        mcad_fail("target_size must be positive and finite");
    if (!(std::isfinite(rO.mMaxChordalDeviation) && rO.mMaxChordalDeviation >= 0.0))
        mcad_fail("max_chordal_deviation must be finite and not negative");
    if (rO.mMinPointsPerEdge < 0)
        mcad_fail("min_points_per_edge must not be negative");
    if (!(std::isfinite(rO.mGeometricTolerance) && rO.mGeometricTolerance >= 0.0))
        mcad_fail("tolerance must be finite and not negative");
    if (rO.mMaxPoints < 1)
        mcad_fail("max_points must be positive");
    if (!(std::isfinite(rO.mMinAngleDeg) && rO.mMinAngleDeg >= 0.0 && rO.mMinAngleDeg < 60.0))
        mcad_fail("min_angle must lie in [0, 60) degrees");
    if (rO.mSmoothIterations < 0)
        mcad_fail("smooth_iterations must not be negative");
    if (rO.mSmoothMethod != MeshCadSmoothing::Laplacian &&
        rO.mSmoothMethod != MeshCadSmoothing::Odt)
        mcad_fail("smooth_method must be laplacian or odt");
    // In the layout, not yet implemented: refused rather than ignored.
    const auto not_yet = [](const char* pOption, const char* pItem) {
        mcad_fail(std::string(pOption) + " is not implemented yet (roadmap " + pItem +
                  "); leave it at its default");
    };
    if (rO.mInsertion != MeshCadInsertion::Circumcentre)
        not_yet("insertion=\"offcentre\"", "6.2.9");
    if (!(std::isfinite(rO.mGrading) && (rO.mGrading == 0.0 || rO.mGrading > 1.0)))
        mcad_fail("grading must be 0 (the default, 1.5) or greater than 1");
    if (!(std::isfinite(rO.mCurvatureDeviation) && rO.mCurvatureDeviation >= 0.0))
        mcad_fail("curvature_deviation must be finite and not negative");
    if (!(std::isfinite(rO.mCellsAcrossGaps) && rO.mCellsAcrossGaps >= 0.0))
        mcad_fail("cells_across_gaps must be finite and not negative");
    if (rO.mMetric != MeshCadMetric::Mean && rO.mMetric != MeshCadMetric::Varying)
        mcad_fail("metric must be \"mean\" or \"varying\"");
    if (rO.mPoleRings < 0 || rO.mPoleRings > 8)
        mcad_fail("pole_rings must lie in [0, 8]");
    if (!(std::isfinite(rO.mReparamConditionNumber) &&
          (rO.mReparamConditionNumber == 0.0 || rO.mReparamConditionNumber > 1.0)))
        mcad_fail("reparam_condition_number must be 0 (off) or greater than 1");
    if (rO.mOrder < 1 || rO.mOrder > 3)
        mcad_fail("order must be 1, 2 or 3");
    if (rO.mReparamConditionNumber > 0.0 && rO.mOrder > 1)
        mcad_fail("reparam_condition_number is order 1 only (roadmap 6.4.4)");
}

/// One loop sample of a face: its mesh point, its parameters, and the use
/// and sample index it comes from.
struct McadEntry {
    std::int64_t mId;
    Vec2 mUv;
    std::size_t mUse;
    std::size_t mK;
    /// Walked against its use (a merged sliver's loop, 6.5.2): the piece to
    /// the next sample is `mK - 1`.
    bool mReversed = false;
};

/// A face's loops as closed polylines in its parameters, and the plane they
/// are triangulated in: x = L (u, v), L the Cholesky factor of the mean first
/// fundamental form over the loops' box, so a triangulation Delaunay in x is
/// Delaunay in (an average of) the surface's own metric. In raw parameters,
/// whose scales can differ by orders of magnitude, a Delaunay triangulation
/// links points far apart on the surface.
struct McadPlane {
    std::vector<std::vector<McadEntry>> mLoops;
    std::size_t mTotal = 0;
    std::vector<std::array<Vec3, 2>> mJac;  ///< (S_u, S_v) on a 9 x 9 grid of the box
    std::vector<Vec3> mJacAt;               ///< S at each of them
    double mL11 = 1.0;
    double mL12 = 0.0;
    double mL22 = 1.0;
    std::vector<Vec2> mXs;  ///< the loops' samples in x, loop after loop
    double mXlo[2] = {0.0, 0.0};
    double mXhi[2] = {0.0, 0.0};
    bool mInRange = true;  ///< every mapped sample within the exact predicates' range
    /// Meshed in a flattening of its coarse mesh (6.4.4): `ToUv` gives the
    /// flattening's coordinates, which the disc carries to the parameters.
    const detail::McadDisc* pDisc = nullptr;

    Vec2 ToX(const Vec2& rUv) const { return {mL11 * rUv[0] + mL12 * rUv[1], mL22 * rUv[1]}; }
    Vec2 ToUv(const Vec2& rX) const {
        const double v = rX[1] / mL22;
        return {(rX[0] - mL12 * v) / mL11, v};
    }
};

/// Whether the surface's metric at a point, (S_u, S_v) there, departs more
/// than twofold from the face's plane: the condition number of J L^-1 above
/// 2, where one affine map cannot follow the surface (a varying metric,
/// roadmap 6.4.1, is what can).
bool mcad_stretched(const Vec3& rSu, const Vec3& rSv, const McadPlane& rPl) {
    Vec3 c1{};
    Vec3 c2{};
    for (std::size_t i = 0; i < 3; ++i) {
        c1[i] = rSu[i] / rPl.mL11;
        c2[i] = rSv[i] / rPl.mL22 - rSu[i] * rPl.mL12 / (rPl.mL11 * rPl.mL22);
    }
    const double a = c1[0] * c1[0] + c1[1] * c1[1] + c1[2] * c1[2];
    const double b = c1[0] * c2[0] + c1[1] * c2[1] + c1[2] * c2[2];
    const double d = c2[0] * c2[0] + c2[1] * c2[1] + c2[2] * c2[2];
    const double root = std::sqrt((a - d) * (a - d) + 4.0 * b * b);
    const double lmax = 0.5 * (a + d + root);
    const double lmin = 0.5 * (a + d - root);
    return !(lmin > 0.0 && lmax <= 4.0 * lmin);
}

/// Build a face's plane from its uses' current samples (`mParams`, `mIds`).
McadPlane mcad_face_plane(const McadFace& rF, const std::vector<McadUse>& rUses,
                          const std::vector<std::optional<detail::TrimTrace>>& rTraces) {
    McadPlane pl;
    double lo[2] = {std::numeric_limits<double>::infinity(),
                    std::numeric_limits<double>::infinity()};
    double hi[2] = {-lo[0], -lo[1]};
    for (std::size_t l = 0; l < rF.mLoops.size(); ++l) {
        const auto& r_loop = rF.mLoops[l];
        std::vector<McadEntry> poly;
        for (std::size_t i = 0; i < r_loop.size(); ++i) {
            const std::size_t u = r_loop[i];
            const McadUse& r_u = rUses[u];
            const detail::TrimTrace& r_t = *rTraces[u];
            const bool back = !rF.mLoopReversed.empty() && rF.mLoopReversed[l][i] != 0;
            for (std::size_t s = 0; s + 1 < r_u.mIds.size(); ++s) {
                const std::size_t j = back ? r_u.mIds.size() - 1 - s : s;
                const std::int64_t id = r_u.mIds[j];
                if (rF.pSuppressed && id >= 0 &&
                    static_cast<std::size_t>(id) < rF.pSuppressed->size() &&
                    (*rF.pSuppressed)[static_cast<std::size_t>(id)])
                    continue;
                const Vec2 uv = r_t.Uv(r_u.mParams[j]);
                poly.push_back({id, uv, u, j, back});
                for (int k = 0; k < 2; ++k) {
                    lo[k] = std::min(lo[k], uv[k]);
                    hi[k] = std::max(hi[k], uv[k]);
                }
            }
        }
        if (!poly.empty()) {
            pl.mTotal += poly.size();
            pl.mLoops.push_back(std::move(poly));
        }
    }
    if (pl.mTotal < 3)
        return pl;

    const detail::BrepSurfaceEval& r_s = *rF.pSurface;
    constexpr int kSamples = 9;
    double me = 0.0;
    double mf = 0.0;
    double mg = 0.0;
    for (int j = 0; j < kSamples; ++j)
        for (int i = 0; i < kSamples; ++i) {
            const double u = lo[0] + (hi[0] - lo[0]) * double(i) / double(kSamples - 1);
            const double v = lo[1] + (hi[1] - lo[1]) * double(j) / double(kSamples - 1);
            double x[3];
            Vec3 su{};
            Vec3 sv{};
            r_s.Eval(u, v, x, su.data(), sv.data());
            pl.mJac.push_back({su, sv});
            pl.mJacAt.push_back({x[0], x[1], x[2]});
            me += su[0] * su[0] + su[1] * su[1] + su[2] * su[2];
            mf += su[0] * sv[0] + su[1] * sv[1] + su[2] * sv[2];
            mg += sv[0] * sv[0] + sv[1] * sv[1] + sv[2] * sv[2];
        }
    const double ns = double(pl.mJac.size());
    me /= ns;
    mf /= ns;
    mg /= ns;
    if (me > 0.0 && mg > 0.0 && me - mf * mf / mg > 1e-12 * me) {
        pl.mL11 = std::sqrt(me);
        pl.mL12 = mf / pl.mL11;
        pl.mL22 = std::sqrt(mg - pl.mL12 * pl.mL12);
    }
    pl.mXs.reserve(pl.mTotal);
    pl.mXlo[0] = pl.mXlo[1] = std::numeric_limits<double>::infinity();
    pl.mXhi[0] = pl.mXhi[1] = -std::numeric_limits<double>::infinity();
    for (const auto& r_poly : pl.mLoops)
        for (const McadEntry& r_e : r_poly) {
            pl.mXs.push_back(pl.ToX(r_e.mUv));
            pl.mInRange = pl.mInRange && detail::cdt2d_in_range(pl.mXs.back()[0]) &&
                          detail::cdt2d_in_range(pl.mXs.back()[1]);
            for (int k = 0; k < 2; ++k) {
                pl.mXlo[k] = std::min(pl.mXlo[k], pl.mXs.back()[static_cast<std::size_t>(k)]);
                pl.mXhi[k] = std::max(pl.mXhi[k], pl.mXs.back()[static_cast<std::size_t>(k)]);
            }
        }
    if (!pl.mInRange) {
        // Real CAD writes parameters a hair off zero (7e-44 where 0 was
        // meant), under the 2^-128 the predicates hold for (roadmap 6.5.7):
        // a coordinate 1e-30 of the plane's size from zero is zero. Only a
        // face that would fail otherwise takes this, so every other face is
        // triangulated on the coordinates it always was.
        double size = 0.0;
        for (const Vec2& r_x : pl.mXs)
            size = std::max({size, std::fabs(r_x[0]), std::fabs(r_x[1])});
        const double flush = 1e-30 * size;
        pl.mInRange = true;
        for (Vec2& r_x : pl.mXs)
            for (double& r_c : r_x) {
                if (std::fabs(r_c) < flush)
                    r_c = 0.0;
                pl.mInRange = pl.mInRange && detail::cdt2d_in_range(r_c);
            }
        for (int k = 0; k < 2; ++k) {
            pl.mXlo[k] = std::min(pl.mXlo[k], 0.0);
            pl.mXhi[k] = std::max(pl.mXhi[k], 0.0);
        }
    }
    return pl;
}

/// A face's plane in the flattening of its coarse mesh (6.4.4): the loops
/// keep their parameters (the trim-range check reads them), and their samples
/// are their coarse vertices' disc points mapped by the disc's own plane.
/// False when a sample is not in the coarse mesh.
bool mcad_disc_plane(McadPlane& rPl, const detail::McadDisc& rDisc) {
    const detail::McadCoarseMesh& r_c = *rDisc.pCoarse;
    rPl.pDisc = &rDisc;
    rPl.mL11 = rDisc.mL11;
    rPl.mL12 = rDisc.mL12;
    rPl.mL22 = rDisc.mL22;
    rPl.mJac = rDisc.mJac;
    rPl.mJacAt = rDisc.mJacAt;
    rPl.mXs.clear();
    rPl.mXlo[0] = rPl.mXlo[1] = std::numeric_limits<double>::infinity();
    rPl.mXhi[0] = rPl.mXhi[1] = -std::numeric_limits<double>::infinity();
    rPl.mInRange = true;
    for (const auto& r_poly : rPl.mLoops)
        for (const McadEntry& r_e : r_poly) {
            const auto it = r_c.mOfSample.find({r_e.mUse, r_e.mK});
            if (it == r_c.mOfSample.end())
                return false;
            rPl.mXs.push_back(rPl.ToX(rDisc.mP[it->second]));
            rPl.mInRange = rPl.mInRange && detail::cdt2d_in_range(rPl.mXs.back()[0]) &&
                           detail::cdt2d_in_range(rPl.mXs.back()[1]);
            for (int k = 0; k < 2; ++k) {
                rPl.mXlo[k] = std::min(rPl.mXlo[k], rPl.mXs.back()[static_cast<std::size_t>(k)]);
                rPl.mXhi[k] = std::max(rPl.mXhi[k], rPl.mXs.back()[static_cast<std::size_t>(k)]);
            }
        }
    return true;
}

/// A host's and a sliver's coarse meshes glued along the edges they share
/// (6.5.2): the sliver's vertices at the shared samples (`rShared`, mesh ids)
/// become the host's, its triangles owner 1, turned round when they run the
/// same way as the host's along a shared edge so the two are oriented alike.
/// Empty `rWhy` on success.
detail::McadCoarseMesh mcad_glue(const detail::McadCoarseMesh& rHost,
                                 const detail::McadCoarseMesh& rSliver,
                                 const std::set<std::int64_t>& rShared, std::string& rWhy) {
    detail::McadCoarseMesh out = rHost;
    std::map<std::int64_t, std::size_t> host_at;
    for (std::size_t v = 0; v < rHost.mId.size(); ++v)
        if (rShared.count(rHost.mId[v]) && !host_at.emplace(rHost.mId[v], v).second) {
            rWhy = "a shared point appears twice in the host's coarse mesh";
            return out;
        }
    std::vector<std::size_t> to(rSliver.mX.size());
    for (std::size_t v = 0; v < rSliver.mX.size(); ++v) {
        const auto it = host_at.find(rSliver.mId[v]);
        if (rShared.count(rSliver.mId[v]) && it != host_at.end()) {
            to[v] = it->second;
        } else {
            to[v] = out.mX.size();
            out.mX.push_back(rSliver.mX[v]);
            out.mId.push_back(rSliver.mId[v]);
        }
    }
    // The orientation along the shared edges: a directed edge both meshes
    // run the same way means the sliver's triangles turn the other way.
    std::set<std::pair<std::size_t, std::size_t>> host_dir;
    for (const auto& r_t : rHost.mTris)
        for (std::size_t k = 0; k < 3; ++k)
            host_dir.emplace(r_t[k], r_t[(k + 1) % 3]);
    std::size_t same = 0, opposite = 0;
    for (const auto& r_t : rSliver.mTris)
        for (std::size_t k = 0; k < 3; ++k) {
            const std::size_t a = to[r_t[k]], b = to[r_t[(k + 1) % 3]];
            if (host_dir.count({a, b}))
                ++same;
            if (host_dir.count({b, a}))
                ++opposite;
        }
    if (same + opposite == 0) {
        rWhy = "the two coarse meshes share no edge";
        return out;
    }
    if (same > 0 && opposite > 0) {
        rWhy = "the two faces are not oriented alike along their shared edges";
        return out;
    }
    const bool turn = same > 0;
    for (std::size_t t = 0; t < rSliver.mTris.size(); ++t) {
        std::array<std::size_t, 3> tri{to[rSliver.mTris[t][0]], to[rSliver.mTris[t][1]],
                                       to[rSliver.mTris[t][2]]};
        auto uv = rSliver.mCornerUv[t];
        std::int8_t pole = rSliver.mPoleCorner[t];
        if (tri[0] == tri[1] || tri[1] == tri[2] || tri[0] == tri[2])
            continue;
        if (turn) {
            std::swap(tri[1], tri[2]);
            std::swap(uv[1], uv[2]);
            if (pole == 1 || pole == 2)
                pole = static_cast<std::int8_t>(3 - pole);
        }
        out.mTris.push_back(tri);
        out.mCornerUv.push_back(uv);
        out.mOwner.push_back(1);
        out.mPoleCorner.push_back(pole);
    }
    out.mReversed = {rHost.mReversed.empty() ? char(0) : rHost.mReversed[0], char(turn ? 1 : 0)};
    for (const auto& [key, seg] : rSliver.mPoleSegment)
        out.mPoleSegment[{1, to[key.second]}] = seg;
    for (const auto& [key, v] : rSliver.mOfSample)
        out.mOfSample[key] = to[v];
    if (!rSliver.mWhy.empty())
        rWhy = rSliver.mWhy;
    return out;
}

/// What the crossing search asks of one face.
struct McadCrossings {
    std::vector<std::array<std::size_t, 2>> mAsks;  ///< (edge, master piece) to halve
    std::vector<std::size_t> mDrop;                 ///< embedded uses that cannot be separated
    bool mStuck = false;  ///< its loops meet where the curves themselves meet
};

/**
 * @brief Find where a face's sampled loops and embedded edges meet in its
 * plane, so they can be split before the triangulation refuses them.
 *
 * Each segment of a crossing pair of loops asks for its edge's master piece
 * to be halved. When a piece involved is already within the tolerance, the
 * curves themselves meet there and the face is stuck (it stays
 * self_intersecting). A pair with an embedded edge is not split: that edge
 * (the later of two embedded edges) is dropped, since refining the loops
 * around a curve that really crosses them or another embedded curve (two
 * stiffeners crossing, a curve leaving its face) would only crowd points on
 * boundaries the neighbouring faces share. A face whose loop uses failed is
 * skipped, as is an embedded use that failed. Two samples are one vertex
 * when they are one mesh point at the same place, as the triangulation
 * would merge them.
 */
McadCrossings mcad_face_crossings(const McadFace& rF, const std::vector<McadUse>& rUses,
                                  const std::vector<std::optional<detail::TrimTrace>>& rTraces,
                                  const std::vector<McadEdge>& rEdges, double Tol) {
    McadCrossings out;
    for (const auto& r_loop : rF.mLoops)
        for (const std::size_t u : r_loop)
            if (rUses[u].mFail != MeshCadFaceStatus::Ok)
                return out;
    const McadPlane pl = mcad_face_plane(rF, rUses, rTraces);
    if (pl.mTotal < 3 || !pl.mInRange)
        return out;  // the triangulation refuses the face
    std::vector<Vec2> pts;
    std::vector<std::array<std::int64_t, 2>> segs;
    std::vector<std::array<std::size_t, 2>> origin;  // (use, sample) a segment starts at
    std::map<std::pair<std::int64_t, Vec2>, std::int64_t> vertex;
    const auto vertex_of = [&](std::int64_t Id, const Vec2& rX) {
        const auto [it, fresh] =
            vertex.emplace(std::make_pair(Id, rX), static_cast<std::int64_t>(pts.size()));
        if (fresh)
            pts.push_back(rX);
        return it->second;
    };
    std::size_t at = 0;
    for (const auto& r_poly : pl.mLoops) {
        std::vector<std::int64_t> idx(r_poly.size());
        for (std::size_t i = 0; i < r_poly.size(); ++i)
            idx[i] = vertex_of(r_poly[i].mId, pl.mXs[at + i]);
        for (std::size_t i = 0; i < r_poly.size(); ++i) {
            const std::int64_t a = idx[i];
            const std::int64_t b = idx[(i + 1) % r_poly.size()];
            if (a == b)
                continue;
            segs.push_back({a, b});
            origin.push_back({r_poly[i].mUse, r_poly[i].mK});
        }
        at += r_poly.size();
    }
    for (const std::size_t u : rF.mEmbedded) {
        const McadUse& r_u = rUses[u];
        if (!rTraces[u] || r_u.mFail != MeshCadFaceStatus::Ok || r_u.mIds.empty())
            continue;
        std::int64_t prev = -1;
        for (std::size_t j = 0; j < r_u.mIds.size(); ++j) {
            const Vec2 x = pl.ToX(rTraces[u]->Uv(r_u.mParams[j]));
            if (!detail::cdt2d_in_range(x[0]) || !detail::cdt2d_in_range(x[1]))
                break;  // left out when the face is meshed
            const std::int64_t v = vertex_of(r_u.mIds[j], x);
            if (j > 0 && v != prev) {
                segs.push_back({prev, v});
                origin.push_back({u, j - 1});
            }
            prev = v;
        }
    }
    for (const auto& r_c : detail::cdt2d_find_crossings(pts, segs)) {
        std::array<std::size_t, 2> ask[2];
        bool tiny[2];
        bool embedded[2];
        std::size_t use[2];
        for (int i = 0; i < 2; ++i) {
            const auto& r_o = origin[static_cast<std::size_t>(r_c[static_cast<std::size_t>(i)])];
            const McadUse& r_u = rUses[r_o[0]];
            const auto edge = static_cast<std::size_t>(r_u.mEdge);
            const std::size_t n = r_u.mParams.size();
            const std::size_t piece = r_u.mSame ? r_o[1] : n - 2 - r_o[1];
            const auto& r_pts = rEdges[edge].mSamples.mPoints;
            ask[i] = {edge, piece};
            tiny[i] = !(mcad_dist(r_pts[piece], r_pts[piece + 1]) > Tol);
            embedded[i] = r_u.mEmbedded;
            use[i] = r_o[0];
        }
        if (embedded[0] || embedded[1]) {
            out.mDrop.push_back(embedded[1] && (!embedded[0] || use[1] > use[0]) ? use[1] : use[0]);
            continue;
        }
        if (tiny[0] || tiny[1]) {
            out.mStuck = true;
            continue;
        }
        out.mAsks.push_back(ask[0]);
        out.mAsks.push_back(ask[1]);
    }
    if (out.mStuck || !out.mDrop.empty())
        out.mAsks.clear();
    return out;
}

/// Mesh one face in its parameters (see the file comment of mesh_cad.hpp).
/// With `pSnapshot`, stop at the unrefined lattice and record it there as
/// the face's coarse mesh (6.4.4); with `pDisc`, mesh in the flattening of
/// that coarse mesh instead of the parameters.
void mcad_mesh_face(McadFace& rF, const BrepFace& rBrepFace, const std::vector<McadUse>& rUses,
                    const std::vector<std::optional<detail::TrimTrace>>& rTraces,
                    const std::vector<Vec3>& rBoundary, const MeshCadSurfaceOptions& rO,
                    const McadQuality& rQ, detail::McadCoarseMesh* pSnapshot = nullptr,
                    const detail::McadDisc* pDisc = nullptr) {
    MeshCadFaceRecord& r_rec = rF.mRecord;
    const auto fail = [&](MeshCadFaceStatus S, const std::string& rWhy) {
        r_rec.mStatus = S;
        r_rec.mReason = rWhy;
        rF.mTriangles.clear();
        rF.mInterior.clear();
        rF.mAsks.clear();
        rF.mLopAsks.clear();
        rF.mLopped.clear();
        rF.mCurvedEdges.clear();
        rF.mCurvedCentre.clear();
        rF.mCurvedNormal.clear();
        rF.mHighOrder.clear();
        rF.mTriOwner.clear();
    };
    for (const auto& r_loop : rF.mLoops)
        for (const std::size_t u : r_loop)
            if (rUses[u].mFail != MeshCadFaceStatus::Ok)
                return fail(rUses[u].mFail, rUses[u].mWhy);
    const double h = rO.mTargetSize;
    rF.mPoleRingsPlaced = false;
    rF.mCollapsed = false;
    rF.mFolded = false;

    McadPlane pl_built = mcad_face_plane(rF, rUses, rTraces);
    if (pDisc && pl_built.mTotal >= 3 && !mcad_disc_plane(pl_built, *pDisc))
        return fail(MeshCadFaceStatus::Other, "a loop sample is not in its coarse mesh");
    const McadPlane& pl = pl_built;
    if (pl.mTotal < 3)
        return fail(MeshCadFaceStatus::Empty, "its loops have fewer than three points");
    if (!pl.mInRange)
        return fail(MeshCadFaceStatus::Unsupported,
                    "its parameters lie outside the range the exact predicates hold in");
    {
        // The surface exists on its knot range only: past it the patch is
        // evaluated clamped, so a trim that leaves the range by more than a
        // thousandth of it would put points off the face.
        const auto& r_box = rF.pSurface->Box();
        const double du = 1e-3 * (r_box[1] - r_box[0]);
        const double dv = 1e-3 * (r_box[3] - r_box[2]);
        for (const auto& r_poly : pl.mLoops)
            for (const McadEntry& r_e : r_poly)
                if (static_cast<std::int64_t>(rUses[r_e.mUse].mFace) == r_rec.mFace &&
                    !(r_e.mUv[0] >= r_box[0] - du && r_e.mUv[0] <= r_box[1] + du &&
                      r_e.mUv[1] >= r_box[2] - dv && r_e.mUv[1] <= r_box[3] + dv))
                    return fail(MeshCadFaceStatus::Unsupported,
                                "its trims leave its surface's parameter range");
    }
    const auto& loops = pl.mLoops;
    const auto& jac = pl.mJac;
    const double l11 = pl.mL11;
    const double l12 = pl.mL12;
    const double l22 = pl.mL22;
    const double* xlo = pl.mXlo;
    const double* xhi = pl.mXhi;
    const detail::BrepSurfaceEval& r_s = *rF.pSurface;
    // The surface at a point of the plane, and its derivatives along the
    // plane's coordinates before L: in the parameters, or, in a flattening
    // (6.4.4), through the coarse triangle holding the point to its owner's
    // parameters, the derivatives by the chain rule. False when a flattening
    // has no triangle there (the nearest one's clamped values are used).
    const auto surface_x = [&](const detail::Cdt2d::Vec2& rX, double* pX, double* pSu,
                               double* pSv) {
        const Vec2 q = pl.ToUv(rX);
        if (!pl.pDisc) {
            r_s.Eval(q[0], q[1], pX, pSu, pSv);
            return true;
        }
        const detail::McadDisc::Hit hit = pl.pDisc->Locate(q);
        Vec3 su{}, sv{};
        pl.pDisc->mSurfaces[static_cast<std::size_t>(hit.mOwner)]->Eval(
            hit.mUv[0], hit.mUv[1], pX, pSu ? su.data() : nullptr, pSv ? sv.data() : nullptr);
        if (pSu && pSv)
            for (std::size_t k = 0; k < 3; ++k) {
                pSu[k] = su[k] * hit.mD[0] + sv[k] * hit.mD[2];
                pSv[k] = su[k] * hit.mD[1] + sv[k] * hit.mD[3];
            }
        return hit.mInside;
    };
    // The surface and parameters at a point of the plane (its owner's, in a
    // flattening).
    const auto params_x =
        [&](const detail::Cdt2d::Vec2& rX) -> std::pair<const detail::BrepSurfaceEval*, Vec2> {
        const Vec2 q = pl.ToUv(rX);
        if (!pl.pDisc)
            return {&r_s, q};
        const detail::McadDisc::Hit hit = pl.pDisc->Locate(q);
        return {pl.pDisc->mSurfaces[static_cast<std::size_t>(hit.mOwner)], hit.mUv};
    };
    detail::Cdt2d cdt(xlo[0], xlo[1], xhi[0], xhi[1]);
    // A varying metric (6.4.1): Lawson's flips decide in the surface's own
    // metric at each edge, the first fundamental form in the plane's
    // coordinates, (J L^-1)^T (J L^-1); the size is a scale the in-circle
    // test does not see. A flattening (6.4.4) is meshed in its mean metric:
    // it is close to conformal, which an in-circle test does not see either.
    const bool varying = rO.mMetric == MeshCadMetric::Varying && !rF.mMeanMetric && !pl.pDisc;
    if (varying)
        cdt.SetMetric([&](const detail::Cdt2d::Vec2& rX) {
            const Vec2 uv = pl.ToUv(rX);
            Vec3 x{}, su{}, sv{};
            r_s.Eval(uv[0], uv[1], x.data(), su.data(), sv.data());
            Vec3 c1{}, c2{};
            for (std::size_t k = 0; k < 3; ++k) {
                c1[k] = su[k] / l11;
                c2[k] = sv[k] / l22 - su[k] * l12 / (l11 * l22);
            }
            return std::array<double, 3>{c1[0] * c1[0] + c1[1] * c1[1] + c1[2] * c1[2],
                                         c1[0] * c2[0] + c1[1] * c2[1] + c1[2] * c2[2],
                                         c2[0] * c2[0] + c2[1] * c2[1] + c2[2] * c2[2]};
        });
    const std::vector<std::int64_t> vid = cdt.InsertBatch(pl.mXs);
    // Each triangulation vertex is one mesh point: a boundary vertex holds
    // its id, an interior one -(1 + its index in mInterior).
    std::unordered_map<std::int64_t, std::int64_t> mesh_of;
    {
        std::size_t at = 0;
        for (const auto& r_poly : loops)
            for (const McadEntry& r_e : r_poly) {
                const auto [it, fresh] = mesh_of.emplace(vid[at++], r_e.mId);
                if (!fresh && it->second != r_e.mId)
                    return fail(MeshCadFaceStatus::Degenerate,
                                "two distinct boundary points have the same parameters");
            }
    }
    // The use and sample each boundary segment starts at, by its vertices:
    // where a split refinement asks for goes.
    struct Origin {
        std::size_t mUse;
        std::size_t mK;
        std::int64_t mStart;  ///< the vertex at sample k
    };
    std::map<std::pair<std::int64_t, std::int64_t>, Origin> origin;
    const auto seg_key = [](std::int64_t A, std::int64_t B) {
        return std::make_pair(std::min(A, B), std::max(A, B));
    };
    {
        std::size_t at = 0;
        for (const auto& r_poly : loops) {
            for (std::size_t i = 0; i < r_poly.size(); ++i) {
                const std::int64_t a = vid[at + i];
                const std::int64_t b = vid[at + (i + 1) % r_poly.size()];
                if (a == b)
                    continue;
                // A use walked backwards runs from its sample k to k - 1.
                if (r_poly[i].mReversed)
                    origin.emplace(seg_key(a, b), Origin{r_poly[i].mUse, r_poly[i].mK - 1, b});
                else
                    origin.emplace(seg_key(a, b), Origin{r_poly[i].mUse, r_poly[i].mK, a});
                const detail::Cdt2dStatus st = cdt.InsertSegment(a, b);
                if (st == detail::Cdt2dStatus::SelfIntersecting)
                    return fail(MeshCadFaceStatus::SelfIntersecting,
                                "its boundary crosses itself in the face's parameters at this "
                                "target size");
                if (st != detail::Cdt2dStatus::Ok)
                    return fail(MeshCadFaceStatus::Other, detail::cdt2d_status_name(st));
            }
            at += r_poly.size();
        }
    }
    if (cdt.Classify() != detail::Cdt2dStatus::Ok)
        return fail(MeshCadFaceStatus::OpenLoop, "its loops do not bound a domain consistently");
    if (cdt.InsideTriangles().empty())
        return fail(MeshCadFaceStatus::Empty, "its loops enclose no area");

    // Embedded edges, after the classification so the even-odd rule counted
    // the loops alone: every sample strictly inside the domain or on the
    // loop vertex that is its own mesh point, then recovered as constrained
    // segments like the loops (Anglada 1997). One that is not is left out,
    // and the face is meshed without it.
    std::vector<std::array<Vec3, 2>> embedded_segments;
    // (use, sample) of each embedded sample: its vertex and parameters, for
    // the coarse mesh.
    std::map<std::pair<std::size_t, std::size_t>, std::pair<std::int64_t, Vec2>> embedded_vertex;
    for (const std::size_t u : rF.mEmbedded) {
        const McadUse& r_u = rUses[u];
        if (!rTraces[u])
            continue;  // already in mSkipped
        if (r_u.mFail != MeshCadFaceStatus::Ok) {
            rF.mSkipped.emplace_back(u, r_u.mWhy);
            continue;
        }
        std::vector<Vec2> xs;
        bool inside = true;
        for (std::size_t j = 0; j < r_u.mIds.size() && inside; ++j) {
            if (pl.pDisc) {
                // Its samples are vertices of the coarse mesh, or it was left
                // out there.
                const auto it = pl.pDisc->pCoarse->mOfSample.find({u, j});
                if (it == pl.pDisc->pCoarse->mOfSample.end()) {
                    inside = false;
                    break;
                }
                xs.push_back(pl.ToX(pl.pDisc->mP[it->second]));
            } else {
                xs.push_back(pl.ToX(rTraces[u]->Uv(r_u.mParams[j])));
            }
            if (!detail::cdt2d_in_range(xs.back()[0]) || !detail::cdt2d_in_range(xs.back()[1])) {
                inside = false;
                break;
            }
            const detail::Cdt2dWhere w = cdt.Where(xs.back());
            if (w.mPlace == detail::Cdt2dPlace::OnVertex) {
                const auto it = mesh_of.find(w.mVertex);
                inside = it != mesh_of.end() && it->second == r_u.mIds[j];
            } else {
                inside = w.mPlace == detail::Cdt2dPlace::Inside;
            }
        }
        if (!inside) {
            rF.mSkipped.emplace_back(u,
                                     "it leaves the face's domain, or meets its boundary "
                                     "away from a corner");
            continue;
        }
        std::int64_t prev = -1;
        for (std::size_t j = 0; j < xs.size(); ++j) {
            const std::int64_t v = cdt.Insert(xs[j]);
            if (pSnapshot)
                embedded_vertex.emplace(std::make_pair(u, j),
                                        std::make_pair(v, rTraces[u]->Uv(r_u.mParams[j])));
            const auto [it, fresh] = mesh_of.emplace(v, r_u.mIds[j]);
            if (!fresh && it->second != r_u.mIds[j])
                return fail(MeshCadFaceStatus::Other,
                            "an embedded edge has two points at one place in the parameters");
            if (prev >= 0 && v != prev && cdt.InsertSegment(prev, v) != detail::Cdt2dStatus::Ok)
                return fail(MeshCadFaceStatus::Other, "an embedded edge could not be recovered");
            if (prev >= 0 && v != prev)
                origin.emplace(seg_key(prev, v), Origin{u, j - 1, prev});
            if (j > 0)
                embedded_segments.push_back({rBoundary[static_cast<std::size_t>(r_u.mIds[j - 1])],
                                             rBoundary[static_cast<std::size_t>(r_u.mIds[j])]});
            prev = v;
        }
    }

    // Interior points: an equilateral lattice in the mapped plane, its
    // spacing h over the largest stretch of the surface against the map (the
    // largest singular value of J L^-1 on the samples), so no lattice step is
    // longer than h on the surface.
    double kappa = 0.0;
    for (const auto& r_j : jac) {
        Vec3 c1{};
        Vec3 c2{};
        for (std::size_t k = 0; k < 3; ++k) {
            c1[k] = r_j[0][k] / l11;
            c2[k] = r_j[1][k] / l22 - r_j[0][k] * l12 / (l11 * l22);
        }
        const double a = c1[0] * c1[0] + c1[1] * c1[1] + c1[2] * c1[2];
        const double b = c1[0] * c2[0] + c1[1] * c2[1] + c1[2] * c2[2];
        const double d = c2[0] * c2[0] + c2[1] * c2[1] + c2[2] * c2[2];
        const double lmax = 0.5 * (a + d + std::sqrt((a - d) * (a - d) + 4.0 * b * b));
        kappa = std::max(kappa, std::sqrt(lmax));
    }
    if (kappa > 0.0 && xhi[0] > xlo[0] && xhi[1] > xlo[1]) {
        const double dx = h / kappa;
        const double dy = dx * std::sqrt(3.0) / 2.0;
        const double rows = std::floor((xhi[1] - xlo[1]) / dy);
        const double cols = std::floor((xhi[0] - xlo[0]) / dx) + 1.0;
        if (rows * cols > double(rO.mMaxPoints))
            return fail(MeshCadFaceStatus::BudgetExceeded,
                        "its interior lattice would exceed max_points");
        McadGrid grid(h);
        for (const auto& r_poly : loops)
            for (std::size_t i = 0; i < r_poly.size(); ++i) {
                const Vec3& a = rBoundary[static_cast<std::size_t>(r_poly[i].mId)];
                const Vec3& b =
                    rBoundary[static_cast<std::size_t>(r_poly[(i + 1) % r_poly.size()].mId)];
                grid.AddPoint(a);
                grid.AddSegment(a, b);
            }
        for (const auto& r_seg : embedded_segments) {
            grid.AddPoint(r_seg[0]);
            grid.AddPoint(r_seg[1]);
            grid.AddSegment(r_seg[0], r_seg[1]);
        }
        // Graded poles (6.4.2): rings of points at 3-D distances h, 2h, ...
        // from each pole, along the face's parameter lines out of the pole's
        // degenerate segment, each as many points as its length holds at h,
        // so the pole is the centre of a star of near-equilateral triangles
        // instead of a fan's apex. The lattice keeps out of the rings' band.
        std::vector<std::pair<Vec3, double>> ring_band;
        // Not in a coarse mesh (it would depend on them) nor in a
        // flattening, which has no pole.
        if (rO.mPoleRings > 0 && !pSnapshot && !pl.pDisc) {
            std::map<std::int64_t, std::vector<Vec2>> pole_uv;
            Vec2 lo{std::numeric_limits<double>::infinity(),
                    std::numeric_limits<double>::infinity()};
            Vec2 hi{-lo[0], -lo[1]};
            Vec2 centre{0.0, 0.0};
            double count = 0.0;
            for (const auto& r_poly : loops)
                for (const McadEntry& r_e : r_poly) {
                    for (int d = 0; d < 2; ++d) {
                        lo[d] = std::min(lo[d], r_e.mUv[d]);
                        hi[d] = std::max(hi[d], r_e.mUv[d]);
                        centre[d] += r_e.mUv[d];
                    }
                    count += 1.0;
                    if (r_e.mId >= 0 && static_cast<std::size_t>(r_e.mId) < rQ.mPole.size() &&
                        rQ.mPole[static_cast<std::size_t>(r_e.mId)] != 0)
                        pole_uv[r_e.mId].push_back(r_e.mUv);
                }
            centre = {centre[0] / count, centre[1] / count};
            for (const auto& [id, r_uv] : pole_uv) {
                // The degenerate segment, and the way into the face across it.
                Vec2 a = r_uv.front(), b = r_uv.front();
                double far = 0.0;
                for (std::size_t i = 0; i < r_uv.size(); ++i)
                    for (std::size_t j = i + 1; j < r_uv.size(); ++j) {
                        const double d =
                            std::hypot(r_uv[i][0] - r_uv[j][0], r_uv[i][1] - r_uv[j][1]);
                        if (d > far) {
                            far = d;
                            a = r_uv[i];
                            b = r_uv[j];
                        }
                    }
                if (!(far > 0.0))
                    continue;  // a pole at one parameter point has no lines out of it
                Vec2 in{-(b[1] - a[1]) / far, (b[0] - a[0]) / far};
                const Vec2 mid{0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1])};
                if (in[0] * (centre[0] - mid[0]) + in[1] * (centre[1] - mid[1]) < 0.0)
                    in = {-in[0], -in[1]};
                double t_max = std::numeric_limits<double>::infinity();
                for (int d = 0; d < 2; ++d) {
                    if (in[d] > 0.0)
                        t_max = std::min(t_max, (hi[d] - std::min(a[d], b[d])) / in[d]);
                    else if (in[d] < 0.0)
                        t_max = std::min(t_max, (lo[d] - std::max(a[d], b[d])) / in[d]);
                }
                if (!(t_max > 0.0) || !std::isfinite(t_max))
                    continue;
                const Vec3& r_p = rBoundary[static_cast<std::size_t>(id)];
                const double hp = rQ.pSize ? rQ.pSize->At(r_p) : h;
                const auto at = [&](double S, double T) {
                    return Vec2{a[0] + S * (b[0] - a[0]) + T * in[0],
                                a[1] + S * (b[1] - a[1]) + T * in[1]};
                };
                const auto surface = [&](const Vec2& rUv) {
                    Vec3 x{};
                    r_s.Eval(rUv[0], rUv[1], x.data(), nullptr, nullptr);
                    return x;
                };
                // The parameter offset along the line at S where the surface
                // is R from the pole: the first sign change on a scan, then
                // bisection; NaN where the line never gets that far.
                const auto offset = [&](double S, double R) {
                    constexpr int kScan = 64;
                    double t0 = 0.0;
                    for (int j = 1; j <= kScan; ++j) {
                        double t1 = t_max * double(j) / kScan;
                        if (mcad_dist(surface(at(S, t1)), r_p) >= R) {
                            for (int it = 0; it < 60; ++it) {
                                const double tm = 0.5 * (t0 + t1);
                                if (mcad_dist(surface(at(S, tm)), r_p) >= R)
                                    t1 = tm;
                                else
                                    t0 = tm;
                            }
                            return 0.5 * (t0 + t1);
                        }
                        t0 = t1;
                    }
                    return std::numeric_limits<double>::quiet_NaN();
                };
                int placed = 0;
                for (std::int32_t k = 1; k <= rO.mPoleRings; ++k) {
                    const double r = double(k) * hp;
                    constexpr int kTrace = 64;
                    std::vector<double> cum{0.0};
                    Vec3 prev{};
                    bool ok = true;
                    for (int i = 0; i <= kTrace; ++i) {
                        const double t = offset(double(i) / kTrace, r);
                        if (!std::isfinite(t)) {
                            ok = false;
                            break;
                        }
                        const Vec3 x = surface(at(double(i) / kTrace, t));
                        if (i > 0)
                            cum.push_back(cum.back() + mcad_dist(prev, x));
                        prev = x;
                    }
                    if (!ok)
                        break;
                    const auto n = std::max<std::int64_t>(3, std::llround(cum.back() / hp));
                    for (std::int64_t j = 1; j < n; ++j) {
                        const double target = cum.back() * double(j) / double(n);
                        std::size_t i = 0;
                        while (i + 1 < cum.size() && cum[i + 1] < target)
                            ++i;
                        const double span = cum[i + 1] - cum[i];
                        const double sv =
                            (double(i) + (span > 0.0 ? (target - cum[i]) / span : 0.0)) / kTrace;
                        const double t = offset(sv, r);
                        if (!std::isfinite(t))
                            continue;
                        const Vec2 uv = at(sv, t);
                        const Vec3 x = surface(uv);
                        if (grid.Near(x, 0.4 * hp, 0.3 * hp))
                            continue;
                        const std::int64_t v_id = cdt.InsertSteiner(pl.ToX(uv));
                        if (v_id < 0)
                            continue;
                        mesh_of[v_id] = -(1 + static_cast<std::int64_t>(rF.mInterior.size()));
                        rF.mInterior.push_back(x);
                        grid.AddPoint(x);
                    }
                    placed = k;
                }
                if (placed > 0)
                    ring_band.emplace_back(r_p, (double(placed) + 0.5) * hp);
            }
        }
        rF.mPoleRingsPlaced = !ring_band.empty();
        const auto nrows = static_cast<std::int64_t>(rows);
        const auto ncols = static_cast<std::int64_t>(cols);
        for (std::int64_t j = 1; j <= nrows; ++j) {
            const double y = xlo[1] + double(j) * dy;
            for (std::int64_t i = 0; i <= ncols; ++i) {
                const Vec2 xp{xlo[0] + (double(i) + 0.5 * double(j % 2)) * dx, y};
                if (!(xp[0] > xlo[0] && xp[0] < xhi[0]))
                    continue;
                Vec3 x{};
                if (!surface_x(xp, x.data(), nullptr, nullptr))
                    continue;  // outside a flattening's coarse mesh
                const double hx = rQ.pSize ? rQ.pSize->At(x) : h;
                if (grid.Near(x, 0.5 * hx, 0.4 * hx))
                    continue;
                if (std::any_of(ring_band.begin(), ring_band.end(),
                                [&](const auto& rB) { return mcad_dist(x, rB.first) < rB.second; }))
                    continue;
                const std::int64_t v_id = cdt.InsertSteiner(xp);
                if (v_id < 0)
                    continue;
                mesh_of[v_id] = -(1 + static_cast<std::int64_t>(rF.mInterior.size()));
                rF.mInterior.push_back(x);
                grid.AddPoint(x);
            }
        }
        // After a meshing in which every triangle collapsed: three by three
        // points at the quarters of the plane's box, whatever h, so the face
        // has an inside (those outside its domain are not inserted).
        if (rF.mForceInterior)
            for (int j = 1; j <= 3; ++j)
                for (int i = 1; i <= 3; ++i) {
                    const Vec2 xp{xlo[0] + 0.25 * double(i) * (xhi[0] - xlo[0]),
                                  xlo[1] + 0.25 * double(j) * (xhi[1] - xlo[1])};
                    const std::int64_t v_id = cdt.InsertSteiner(xp);
                    if (v_id < 0 || mesh_of.count(v_id) > 0)
                        continue;
                    const Vec3 x = [&]() {
                        Vec3 y{};
                        surface_x(xp, y.data(), nullptr, nullptr);
                        return y;
                    }();
                    mesh_of[v_id] = -(1 + static_cast<std::int64_t>(rF.mInterior.size()));
                    rF.mInterior.push_back(x);
                    grid.AddPoint(x);
                }
    }

    const auto pos = [&](std::int64_t Id) -> const Vec3& {
        return Id >= 0 ? rBoundary[static_cast<std::size_t>(Id)]
                       : rF.mInterior[static_cast<std::size_t>(-Id - 1)];
    };
    if (pSnapshot) {
        // The coarse mesh (6.4.4): this unrefined triangulation, each
        // triangle with its corners' own parameters. The consecutive loop
        // samples of one mesh point (a pole: a degenerate trim's two ends)
        // are one vertex, the triangle they collapse is dropped, and the
        // pole's degenerate segment is kept for the collapsed coordinates;
        // the seam's two copies of a point stay apart, as in the parameters.
        detail::McadCoarseMesh& r_c = *pSnapshot;
        r_c = detail::McadCoarseMesh{};
        r_c.mL11 = l11;
        r_c.mL12 = l12;
        r_c.mL22 = l22;
        const std::size_t nv = static_cast<std::size_t>(cdt.NumVertices());
        std::vector<std::int64_t> rep(nv, -1);  // vertex -> its run's first vertex
        std::vector<char> exact(nv, 0);
        std::vector<Vec2> uv_of(nv);
        std::map<std::int64_t, std::vector<Vec2>> run_uv;  // run vertex -> its parameters
        {
            std::size_t at = 0;
            for (const auto& r_poly : loops) {
                const std::size_t n = r_poly.size();
                for (std::size_t i = 0; i < n; ++i) {
                    const auto v = static_cast<std::size_t>(vid[at + i]);
                    uv_of[v] = r_poly[i].mUv;
                    exact[v] = 1;
                }
                // Runs of one mesh id, cyclically: start after a change.
                std::size_t start = 0;
                while (start < n && r_poly[start].mId == r_poly[(start + n - 1) % n].mId)
                    ++start;
                if (start == n)
                    return;  // a loop of one point: nothing to flatten
                for (std::size_t i = 0; i < n;) {
                    const std::size_t a = (start + i) % n;
                    std::size_t len = 1;
                    while (len < n && r_poly[(a + len) % n].mId == r_poly[a].mId)
                        ++len;
                    const std::int64_t head = vid[at + a];
                    for (std::size_t j = 0; j < len; ++j) {
                        const std::int64_t w = vid[at + (a + j) % n];
                        rep[static_cast<std::size_t>(w)] = head;
                        if (len > 1)
                            run_uv[head].push_back(r_poly[(a + j) % n].mUv);
                    }
                    i += len;
                }
                at += n;
            }
        }
        for (const auto& [key, r_v] : embedded_vertex) {
            const auto w = static_cast<std::size_t>(r_v.first);
            uv_of[w] = r_v.second;
            exact[w] = 1;
        }
        std::vector<std::int64_t> coarse_of(nv, -1);
        const auto vertex = [&](std::int64_t V) -> std::size_t {
            std::int64_t r = rep[static_cast<std::size_t>(V)];
            if (r < 0)
                r = V;
            auto& r_slot = coarse_of[static_cast<std::size_t>(r)];
            if (r_slot < 0) {
                r_slot = static_cast<std::int64_t>(r_c.mX.size());
                const std::int64_t id = mesh_of.at(r);
                r_c.mX.push_back(pos(id));
                r_c.mId.push_back(id);
            }
            return static_cast<std::size_t>(r_slot);
        };
        for (const auto& r_t : cdt.InsideTriangles()) {
            std::array<std::size_t, 3> tri{};
            std::array<Vec2, 3> uv{};
            std::int8_t pole_corner = -1;
            int poles = 0;
            for (std::size_t k = 0; k < 3; ++k) {
                tri[k] = vertex(r_t[k]);
                const auto w = static_cast<std::size_t>(r_t[k]);
                uv[k] = exact[w] ? uv_of[w] : pl.ToUv(cdt.Point(r_t[k]));
                const std::int64_t r = rep[w] < 0 ? r_t[k] : rep[w];
                if (run_uv.count(r)) {
                    pole_corner = static_cast<std::int8_t>(k);
                    ++poles;
                }
            }
            if (tri[0] == tri[1] || tri[1] == tri[2] || tri[0] == tri[2])
                continue;  // a pole collapses it
            if (poles > 1) {
                r_c.mWhy = "a triangle of its coarse mesh touches two poles";
                pole_corner = -1;
            }
            r_c.mTris.push_back(tri);
            r_c.mCornerUv.push_back(uv);
            r_c.mOwner.push_back(0);
            r_c.mPoleCorner.push_back(pole_corner);
        }
        r_c.mReversed = {0};
        for (const auto& [head, r_uv] : run_uv) {
            // The segment: its two samples farthest apart.
            std::array<Vec2, 2> seg{r_uv.front(), r_uv.front()};
            double far = -1.0;
            for (std::size_t i = 0; i < r_uv.size(); ++i)
                for (std::size_t j = i + 1; j < r_uv.size(); ++j) {
                    const double d = std::hypot(r_uv[i][0] - r_uv[j][0], r_uv[i][1] - r_uv[j][1]);
                    if (d > far) {
                        far = d;
                        seg = {r_uv[i], r_uv[j]};
                    }
                }
            const std::int64_t c = coarse_of[static_cast<std::size_t>(head)];
            if (c >= 0)
                r_c.mPoleSegment[{0, static_cast<std::size_t>(c)}] = seg;
        }
        {
            std::size_t at = 0;
            for (const auto& r_poly : loops) {
                const auto coarse_at = [&](std::size_t I) {
                    const std::int64_t w = vid[at + I];
                    const std::int64_t r =
                        rep[static_cast<std::size_t>(w)] < 0 ? w : rep[static_cast<std::size_t>(w)];
                    return coarse_of[static_cast<std::size_t>(r)];
                };
                for (std::size_t i = 0; i < r_poly.size(); ++i) {
                    const std::int64_t c = coarse_at(i);
                    if (c >= 0)
                        r_c.mOfSample[{r_poly[i].mUse, r_poly[i].mK}] = static_cast<std::size_t>(c);
                    // A use's last sample is the next use's first (a loop
                    // walked against its uses reads it, 6.5.2).
                    const std::size_t next = (i + 1) % r_poly.size();
                    const std::int64_t cn = coarse_at(next);
                    if (r_poly[next].mUse != r_poly[i].mUse && cn >= 0)
                        r_c.mOfSample[{r_poly[i].mUse, rUses[r_poly[i].mUse].mIds.size() - 1}] =
                            static_cast<std::size_t>(cn);
                }
                at += r_poly.size();
            }
        }
        for (const auto& [key, r_v] : embedded_vertex) {
            const std::int64_t c = coarse_of[static_cast<std::size_t>(r_v.first)];
            if (c >= 0)
                r_c.mOfSample[key] = static_cast<std::size_t>(c);
        }
        // A pole corner whose segment is a point needs no collapsed
        // coordinates.
        for (std::size_t t = 0; t < r_c.mTris.size(); ++t) {
            const int c = r_c.mPoleCorner[t];
            if (c >= 0 && !r_c.mPoleSegment.count({0, r_c.mTris[t][static_cast<std::size_t>(c)]}))
                r_c.mPoleCorner[t] = -1;
        }
        r_c.mFilled = true;
        return;
    }
    const auto lift = [&](const detail::Cdt2d::Vec2& rX) {
        Vec3 x{};
        surface_x(rX, x.data(), nullptr, nullptr);
        return x;
    };
    // A flattening (6.4.4) has no pole: the pole is one vertex of it.
    const auto pole = [&](std::int64_t Id) {
        return !pl.pDisc && Id >= 0 && static_cast<std::size_t>(Id) < rQ.mPole.size() &&
               rQ.mPole[static_cast<std::size_t>(Id)] != 0;
    };
    // A pole is a fan's apex, where the face's parameters collapse: the
    // triangles within h of it (on the surface) are left to graded poles
    // (roadmap 6.4.2) and reported with the pole. Refinement there only
    // crowds the pole with points, since no point in the plane makes the
    // fan's thin triangles round.
    std::vector<std::int64_t> face_poles;
    for (const auto& r_poly : loops)
        for (const McadEntry& r_e : r_poly)
            if (pole(r_e.mId) &&
                std::find(face_poles.begin(), face_poles.end(), r_e.mId) == face_poles.end())
                face_poles.push_back(r_e.mId);
    std::sort(face_poles.begin(), face_poles.end());
    // The reach of a pole: h, or farther where the surface's metric departs
    // more than twofold from the plane's (the condition number of J L^-1
    // over 2 on the Jacobian samples), the region one affine map cannot
    // follow (a varying metric, roadmap 6.4.1, is what can).
    double pole_reach = h;
    if (rF.mPoleRingsPlaced) {
        // Its rings make a pole a star's centre: only the triangles at the
        // pole itself, whose third corners the rings place, are left alone.
        pole_reach = 0.5 * h;
    } else if (!face_poles.empty())
        for (std::size_t k = 0; k < pl.mJac.size(); ++k) {
            if (!mcad_stretched(pl.mJac[k][0], pl.mJac[k][1], pl))
                continue;
            double nearest = std::numeric_limits<double>::infinity();
            for (const std::int64_t p : face_poles)
                nearest = std::min(nearest, mcad_dist(pl.mJacAt[k], pos(p)));
            pole_reach = std::max(pole_reach, nearest);
        }
    const auto pole_near = [&](std::int64_t Id) -> std::int64_t {
        std::int64_t best = -1;
        double best_d = std::numeric_limits<double>::infinity();
        for (const std::int64_t p : face_poles) {
            const double d = Id == p ? 0.0 : mcad_dist(pos(Id), pos(p));
            if (d <= pole_reach && d < best_d) {
                best = p;
                best_d = d;
            }
        }
        return best;
    };
    // Whether a triangle lies where the surface stretches its plane more
    // than twofold (at its centroid). There a point placed for the surface
    // can leave triangles as bad on the surface as the one it replaced,
    // each smaller, and refinement need not stop (a planar bilinear patch
    // whose corner stretches 2.4-fold ran away to max_points). Such a
    // triangle is refined only while its shortest edge on the surface is
    // h / 4 or more, and no point it adds there lands closer than h / 4 to
    // its cavity: finitely many such points fit, and what stays under the
    // bound there is reported, until a varying metric (6.4.1).
    const auto stretched_at = [&](const detail::Cdt2d::Vec2& rX) {
        Vec3 x{};
        Vec3 su{};
        Vec3 sv{};
        surface_x(rX, x.data(), su.data(), sv.data());
        return mcad_stretched(su, sv, pl);
    };
    std::map<std::array<std::int64_t, 3>, bool> stretched_cache;
    const auto stretched = [&](const detail::Cdt2d::Tri& rT) {
        std::array<std::int64_t, 3> key = rT;
        std::sort(key.begin(), key.end());
        const auto it = stretched_cache.find(key);
        if (it != stretched_cache.end())
            return it->second;
        const auto& xa = cdt.Point(rT[0]);
        const auto& xb = cdt.Point(rT[1]);
        const auto& xc = cdt.Point(rT[2]);
        const bool out =
            stretched_at({(xa[0] + xb[0] + xc[0]) / 3.0, (xa[1] + xb[1] + xc[1]) / 3.0});
        stretched_cache.emplace(key, out);
        return out;
    };
    const auto left_stretched = [&](const detail::Cdt2d::Tri& rT,
                                    const std::array<std::int64_t, 3>& rIds) {
        if (!(rQ.mMinAngleDeg > 0.0) || varying)
            return false;
        const double shortest =
            std::min({mcad_dist(pos(rIds[0]), pos(rIds[1])), mcad_dist(pos(rIds[1]), pos(rIds[2])),
                      mcad_dist(pos(rIds[2]), pos(rIds[0]))});
        return shortest < 0.25 * h && stretched(rT);
    };
    // A triangle judged on the surface; one a pole collapses, within h of a
    // pole, or under h / 4 where the surface stretches the plane, is exempt
    // (negative).
    const auto exempt = [&](const detail::Cdt2d::Tri& rT, std::array<std::int64_t, 3>& rIds) {
        for (int k = 0; k < 3; ++k)
            rIds[static_cast<std::size_t>(k)] = mesh_of.at(rT[static_cast<std::size_t>(k)]);
        if (rIds[0] == rIds[1] || rIds[1] == rIds[2] || rIds[0] == rIds[2] ||
            pole_near(rIds[0]) >= 0 || pole_near(rIds[1]) >= 0 || pole_near(rIds[2]) >= 0)
            return true;
        return left_stretched(rT, rIds);
    };
    // Delaunay in the plane is not Delaunay on the surface where the
    // face's metric departs from its mean: an edge whose flip raises the
    // smaller of its two triangles' smallest angles on the surface is
    // flipped (when the plane quadrilateral is convex, exact orient2d),
    // and refinement runs again on the new triangles.
    // Every inside edge (a < b): the apex left of a->b and the apex right,
    // and the last quadrilateral judged and not flipped across it, which is
    // judged the same until one of its corners moves.
    // Kept from one call to the next: refinement only inserts points, and
    // every triangle it makes lies around a point it inserted (insertion
    // and Lawson's flips make only those), so the stars of the points
    // inserted since the last call hold every edge whose quadrilateral
    // changed, and only those are written again. An edge refinement
    // removed keeps a stale entry nothing reads.
    struct Side {
        std::array<std::int64_t, 2> mApex{-1, -1};
        std::array<std::int64_t, 4> mKept{-1, -1, -1, -1};
    };
    struct SideHash {
        std::size_t operator()(const std::pair<std::int64_t, std::int64_t>& rE) const {
            return std::hash<std::int64_t>()(rE.first) * 1099511628211ull ^
                   std::hash<std::int64_t>()(rE.second);
        }
    };
    std::unordered_map<std::pair<std::int64_t, std::int64_t>, Side, SideHash> sides;
    // The edges a sweep judges, in key order: those written since they were
    // last judged, then those written again or put off; every other edge's
    // quadrilateral is kept, so a sweep over every edge judges the same
    // edges in the same order.
    std::set<std::pair<std::int64_t, std::int64_t>> pending;
    // The points the sides describe; -1: build them from every triangle (at
    // first, after smoothing moves points, and after a call that stopped
    // with edges still pending, which refinement may since have removed).
    std::int64_t sides_upto = -1;
    // With pole rings a pole is a star's centre whose triangles flips may
    // improve on the surface: only those it collapses are left alone.
    const auto flip_exempt = [&](const detail::Cdt2d::Tri& rT, std::array<std::int64_t, 3>& rIds) {
        if (!rF.mPoleRingsPlaced)
            return exempt(rT, rIds);
        for (int k = 0; k < 3; ++k)
            rIds[static_cast<std::size_t>(k)] = mesh_of.at(rT[static_cast<std::size_t>(k)]);
        if (rIds[0] == rIds[1] || rIds[1] == rIds[2] || rIds[0] == rIds[2])
            return true;
        return left_stretched(rT, rIds);
    };
    const auto surface_flips = [&]() {
        std::int64_t flips = 0;
        const auto write = [&](const detail::Cdt2d::Tri& rT) {
            for (int k = 0; k < 3; ++k) {
                const std::int64_t x = rT[static_cast<std::size_t>((k + 1) % 3)];
                const std::int64_t y = rT[static_cast<std::size_t>((k + 2) % 3)];
                sides[seg_key(x, y)].mApex[x < y ? 0 : 1] = rT[static_cast<std::size_t>(k)];
                pending.insert(seg_key(x, y));
            }
        };
        if (sides_upto < 0) {
            const std::vector<detail::Cdt2d::Tri> inside = cdt.InsideTriangles();
            sides.clear();
            sides.reserve(2 * inside.size());
            pending.clear();
            for (const auto& r_t : inside)
                write(r_t);
        } else {
            for (std::int64_t v = sides_upto; v < cdt.NumVertices(); ++v)
                for (const auto& r_t : cdt.InsideStar(v))
                    write(r_t);
        }
        sides_upto = cdt.NumVertices();
        const auto fixed = cdt.ConstrainedEdges();
        for (int sweep = 0; sweep < 32; ++sweep) {
            // A flip changes the two triangles beside its edge, so the sides
            // recorded for the other edges of its quadrilateral are stale
            // until the next sweep.
            std::vector<char> touched(static_cast<std::size_t>(cdt.NumVertices()), 0);
            const auto was_touched = [&](std::int64_t V) {
                return V >= 0 && static_cast<std::size_t>(V) < touched.size() &&
                       touched[static_cast<std::size_t>(V)] != 0;
            };
            std::vector<std::array<std::int64_t, 4>> flipped;
            const std::vector<std::pair<std::int64_t, std::int64_t>> judge(pending.begin(),
                                                                           pending.end());
            pending.clear();
            for (const auto& r_e : judge) {
                Side& r_side = sides.at(r_e);
                const auto& r_s = r_side.mApex;
                if (was_touched(r_e.first) && was_touched(r_e.second)) {
                    pending.insert(r_e);
                    continue;
                }
                const std::array<std::int64_t, 4> quad{r_e.first, r_e.second, r_s[0], r_s[1]};
                if (r_side.mKept == quad)
                    continue;
                if (r_s[0] < 0 || r_s[1] < 0 ||
                    std::binary_search(fixed.begin(), fixed.end(),
                                       std::array<std::int64_t, 2>{r_e.first, r_e.second})) {
                    r_side.mKept = quad;
                    continue;
                }
                const detail::Cdt2d::Tri t1{r_s[0], r_e.first, r_e.second};
                const detail::Cdt2d::Tri t2{r_s[1], r_e.second, r_e.first};
                const detail::Cdt2d::Tri n1{r_s[0], r_e.first, r_s[1]};
                const detail::Cdt2d::Tri n2{r_s[0], r_s[1], r_e.second};
                std::array<std::int64_t, 3> i1{}, i2{}, j1{}, j2{};
                if (flip_exempt(t1, i1) || flip_exempt(t2, i2) || flip_exempt(n1, j1) ||
                    flip_exempt(n2, j2)) {
                    r_side.mKept = quad;
                    continue;
                }
                const double now = std::min(mcad_min_angle(pos(i1[0]), pos(i1[1]), pos(i1[2])),
                                            mcad_min_angle(pos(i2[0]), pos(i2[1]), pos(i2[2])));
                const double then = std::min(mcad_min_angle(pos(j1[0]), pos(j1[1]), pos(j1[2])),
                                             mcad_min_angle(pos(j2[0]), pos(j2[1]), pos(j2[2])));
                if (then > now + 1e-9 && cdt.FlipEdge(r_e.first, r_e.second)) {
                    flipped.push_back(quad);
                    for (const std::int64_t v : quad)
                        if (v >= 0 && static_cast<std::size_t>(v) < touched.size())
                            touched[static_cast<std::size_t>(v)] = 1;
                } else {
                    r_side.mKept = quad;
                }
            }
            flips += static_cast<std::int64_t>(flipped.size());
            if (flipped.empty())
                break;
            // The next sweep's sides: a flipped edge is gone, and its two new
            // triangles (the new diagonal and the four sides of its
            // quadrilateral) are written again. Those are the only edges a
            // flip changes, and no two flips of a sweep share a triangle (an
            // edge both of whose ends a flip touched is put off).
            for (const auto& r_q : flipped) {
                sides.erase(seg_key(r_q[0], r_q[1]));
                pending.erase(seg_key(r_q[0], r_q[1]));
            }
            for (const auto& r_q : flipped) {
                write({r_q[2], r_q[0], r_q[3]});
                write({r_q[2], r_q[3], r_q[1]});
            }
        }
        if (!pending.empty())
            sides_upto = -1;
        return flips;
    };
    // Folds (6.4.3): on a small, strongly curved face a triangle valid in
    // the plane can lift with its normal against S_u x S_v (both taken with
    // the plane's counter-clockwise winding, which the plane's map keeps).
    // One that turns more than 120 degrees from the surface's normal at its
    // centroid is folded; the margin keeps a sliver across a narrow strip,
    // whose normal is near perpendicular to the surface's, from counting on
    // rounding. A pole's collapsed triangles are not folds.
    // The fold verdict of a triangle depends on its vertices alone, and the
    // surface normal at a vertex on that vertex alone, so both are kept for
    // this meshing: refinement judges a triangle several times (its first
    // pass, its last, the fold check) and a corner is shared by its star.
    std::map<std::array<std::int64_t, 3>, bool> folded_cache;
    std::vector<std::optional<Vec3>> corner_normal;
    const auto folded = [&](const detail::Cdt2d::Tri& rT) {
        std::array<std::int64_t, 3> ids{};
        for (int k = 0; k < 3; ++k)
            ids[static_cast<std::size_t>(k)] = mesh_of.at(rT[static_cast<std::size_t>(k)]);
        if (ids[0] == ids[1] || ids[1] == ids[2] || ids[0] == ids[2])
            return false;
        std::array<std::int64_t, 3> key = rT;
        std::sort(key.begin(), key.end());
        const auto hit = folded_cache.find(key);
        if (hit != folded_cache.end())
            return hit->second;
        const Vec3& pa = pos(ids[0]);
        const Vec3& pb = pos(ids[1]);
        const Vec3& pc = pos(ids[2]);
        const Vec3 e1{pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2]};
        const Vec3 e2{pc[0] - pa[0], pc[1] - pa[1], pc[2] - pa[2]};
        const Vec3 nt{e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2],
                      e1[0] * e2[1] - e1[1] * e2[0]};
        const double nt2 = nt[0] * nt[0] + nt[1] * nt[1] + nt[2] * nt[2];
        const auto against = [&](const Vec3& ns) {
            const double dot = nt[0] * ns[0] + nt[1] * ns[1] + nt[2] * ns[2];
            return dot < -0.5 * std::sqrt(nt2 * (ns[0] * ns[0] + ns[1] * ns[1] + ns[2] * ns[2]));
        };
        // The surface's normal at a vertex's point, once per vertex.
        const auto normal_at = [&](std::int64_t v) -> Vec3 {
            if (corner_normal.size() < static_cast<std::size_t>(cdt.NumVertices()))
                corner_normal.resize(static_cast<std::size_t>(cdt.NumVertices()));
            auto& r_slot = corner_normal[static_cast<std::size_t>(v)];
            if (!r_slot) {
                Vec3 x{}, su{}, sv{};
                surface_x(cdt.Point(v), x.data(), su.data(), sv.data());
                r_slot = Vec3{su[1] * sv[2] - su[2] * sv[1], su[2] * sv[0] - su[0] * sv[2],
                              su[0] * sv[1] - su[1] * sv[0]};
            }
            return *r_slot;
        };
        const auto& xa = cdt.Point(rT[0]);
        const auto& xb = cdt.Point(rT[1]);
        const auto& xc = cdt.Point(rT[2]);
        // Against the surface at its centroid in the parameters, or at every
        // corner: a triangle spanning most of a closed curve can lift onto
        // the short way round, backwards there, while its parameter
        // centroid lies on the far side.
        const Vec2 centre{(xa[0] + xb[0] + xc[0]) / 3.0, (xa[1] + xb[1] + xc[1]) / 3.0};
        Vec3 cx{}, csu{}, csv{};
        surface_x(centre, cx.data(), csu.data(), csv.data());
        const Vec3 cns{csu[1] * csv[2] - csu[2] * csv[1], csu[2] * csv[0] - csu[0] * csv[2],
                       csu[0] * csv[1] - csu[1] * csv[0]};
        const bool out = against(cns) || (against(normal_at(rT[0])) && against(normal_at(rT[1])) &&
                                          against(normal_at(rT[2])));
        folded_cache.emplace(key, out);
        return out;
    };
    // A face that folds is refined even with no angle or size sought: a
    // folded triangle is bad, split at its centroid in the plane, until none
    // is left (the real STEP models of doc/benchmarks.md#cad-corpus fold on
    // small fillets and holes; none of the fixtures does).
    bool fold_refine = false;
    for (const auto& r_t : cdt.InsideTriangles())
        if (folded(r_t)) {
            fold_refine = true;
            break;
        }
    rF.mNeedsRounds = rF.mNeedsRounds || fold_refine;
    detail::Cdt2dRefineReport refined;
    std::set<std::array<std::int64_t, 3>> forced;  // mesh ids, ascending
    std::map<std::int64_t, double> corner_angle;   // mesh point -> plane angle
    std::map<std::int64_t, detail::Cdt2d::Vec2> corner_x;
    // A size field refines too (6.3): with no angle sought, a bound so
    // small no triangle misses it, so only size makes a triangle bad.
    const double refine_deg = rQ.mMinAngleDeg > 0.0 ? rQ.mMinAngleDeg : 1e-3;
    // A face meshed unrefined after spending its points is not refined for
    // its size either (only a fold still is).
    const detail::MeshCadSizeField* const p_refine_size = rQ.mNoSizeRefine ? nullptr : rQ.pSize;
    if (rQ.mMinAngleDeg > 0.0 || p_refine_size || fold_refine) {
        // Ruppert's refinement in the plane, every triangle judged on the
        // surface (Chew 1993; Boissonnat and Oudot 2005): the plane is only
        // where the point is found. The boundary is the edges', shared with
        // other faces, so its splits are asked for, not made.
        detail::Cdt2dRefineOptions ro;
        ro.mMinAngleDeg = refine_deg;
        ro.mSplitSegments = false;
        // An equal share of what max_points leaves after the boundary, so
        // the faces together stay within it whatever each one needs; and at
        // most twice the face's unrefined points (1024 at least) summed over
        // every round, since each round meshes the face again. On a surface
        // termination is not the plane's theorem, and a face that runs away
        // (the sphere at 30 degrees, past the 20.7 Ruppert guarantees) would
        // otherwise fill max_points, 50 million by default, in every round.
        // A face that converges adds at most 1.46 times its points over all
        // rounds where it has 256 or more (1.07 on the fixtures), and at most
        // 521 points on any face (3.07 times on the coarse cylinder's 136 at
        // 30 degrees). A face that reaches either is meshed unrefined and
        // warned.
        const auto unrefined_points = static_cast<std::int64_t>(rF.mInterior.size() + pl.mTotal);
        // A size field asks for about its own count, the area over an
        // equilateral triangle of the local size, per two triangles a point.
        std::int64_t sized_points = 0;
        if (p_refine_size) {
            double count = 0.0;
            for (const auto& r_t : cdt.InsideTriangles()) {
                std::array<std::int64_t, 3> ids{};
                for (int k = 0; k < 3; ++k)
                    ids[static_cast<std::size_t>(k)] = mesh_of.at(r_t[static_cast<std::size_t>(k)]);
                const Vec3& pa = pos(ids[0]);
                const Vec3& pb = pos(ids[1]);
                const Vec3& pc = pos(ids[2]);
                const Vec3 e1{pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2]};
                const Vec3 e2{pc[0] - pa[0], pc[1] - pa[1], pc[2] - pa[2]};
                const Vec3 n{e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2],
                             e1[0] * e2[1] - e1[1] * e2[0]};
                const double area = 0.5 * std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
                // The field integrated over the triangle (Dunavant's 7-point
                // rule), so a layer narrower than the triangle, which its
                // centroid can miss (a size from a field's Hessian, 6.11),
                // still counts.
                static constexpr std::array<std::array<double, 3>, 7> kRule = {
                    {{1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0},
                     {0.059715871789769820, 0.470142064105115090, 0.470142064105115090},
                     {0.470142064105115090, 0.059715871789769820, 0.470142064105115090},
                     {0.470142064105115090, 0.470142064105115090, 0.059715871789769820},
                     {0.797426985353087322, 0.101286507323456339, 0.101286507323456339},
                     {0.101286507323456339, 0.797426985353087322, 0.101286507323456339},
                     {0.101286507323456339, 0.101286507323456339, 0.797426985353087322}}};
                static constexpr std::array<double, 7> kWeight = {0.225,
                                                                  0.132394152788506181,
                                                                  0.132394152788506181,
                                                                  0.132394152788506181,
                                                                  0.125939180544827153,
                                                                  0.125939180544827153,
                                                                  0.125939180544827153};
                double inverse_area = 0.0;
                for (std::size_t q = 0; q < 7; ++q) {
                    Vec3 xq{};
                    for (std::size_t k = 0; k < 3; ++k)
                        xq[k] = kRule[q][0] * pa[k] + kRule[q][1] * pb[k] + kRule[q][2] * pc[k];
                    const double hx = p_refine_size->At(xq);
                    inverse_area += kWeight[q] / (hx * hx);
                }
                count += area * inverse_area / 0.4330127018922193;
            }
            // `count` is the triangles the field asks for: half as many points.
            sized_points = static_cast<std::int64_t>(std::ceil(0.5 * count));
        }
        // A size field asks for a bounded count of points, but each round
        // meshes the face again and re-spends them: a meshing may spend four
        // times the count (twice its triangles), and the rounds together
        // sixteen times it, where a face sized down to a narrow layer needs
        // ten (the 6.11 probe; nine meshes it unrefined).
        const std::int64_t cap_total =
            std::max<std::int64_t>({1024, 2 * unrefined_points, 16 * sized_points}) -
            rF.mRefinedSoFar;
        const std::int64_t cap_left =
            sized_points > 0
                ? std::min(cap_total,
                           std::max<std::int64_t>({1024, 2 * unrefined_points, 4 * sized_points}))
                : cap_total;
        if (cap_left <= 0)
            rF.mRefineBudgetReached = true;
        ro.mMaxPoints = std::min<std::int64_t>(
            std::max<std::int64_t>(0,
                                   (rO.mMaxPoints - static_cast<std::int64_t>(rBoundary.size())) /
                                           std::max<std::int64_t>(1, rQ.mNumFaces) -
                                       static_cast<std::int64_t>(rF.mInterior.size())),
            std::max<std::int64_t>(0, cap_left));
        const double ratio_bound = 1.0 / (2.0 * std::sin(refine_deg * (kMcadPi / 180.0)));
        // A triangle whose longest edge on the surface exceeds the local
        // size (at its centroid) by the slack is bad by its size, reported
        // as worse than the angle bound by that factor so the worst go
        // first (6.3).
        const auto oversize = [&](const std::array<std::int64_t, 3>& rIds) {
            if (!p_refine_size)
                return 0.0;
            const Vec3& pa = pos(rIds[0]);
            const Vec3& pb = pos(rIds[1]);
            const Vec3& pc = pos(rIds[2]);
            const double longest =
                std::max({mcad_dist(pb, pc), mcad_dist(pc, pa), mcad_dist(pa, pb)});
            const Vec3 centre{(pa[0] + pb[0] + pc[0]) / 3.0, (pa[1] + pb[1] + pc[1]) / 3.0,
                              (pa[2] + pb[2] + pc[2]) / 3.0};
            const double s = longest / (detail::kMeshCadSizeSlack * p_refine_size->At(centre));
            return s > 1.0 ? s : 0.0;
        };
        ro.mRatio = [&](const detail::Cdt2d::Tri& rT) {
            std::array<std::int64_t, 3> ids{};
            // A fold first of all, even where the triangles a pole or a
            // stretched region keeps are exempt.
            if (folded(rT)) {
                for (int k = 0; k < 3; ++k)
                    ids[static_cast<std::size_t>(k)] = mesh_of.at(rT[static_cast<std::size_t>(k)]);
                return std::max(mcad_ratio(pos(ids[0]), pos(ids[1]), pos(ids[2])),
                                4.0 * ratio_bound);
            }
            if (exempt(rT, ids))
                return -1.0;
            const double q = mcad_ratio(pos(ids[0]), pos(ids[1]), pos(ids[2]));
            const double s = oversize(ids);
            return s > 0.0 ? std::max(q, ratio_bound * (1.0 + s)) : q;
        };
        // The point: the lifted triangle's own circumcentre, carried to the
        // plane by the triangle's affine map (its barycentric coordinates
        // on the surface triangle, applied to the plane triangle): the
        // circumcentre in the triangle's own metric. The plane's
        // circumcentre is not that where the surface stretches the plane
        // unevenly (near a pole), and refinement on it need not stop there.
        ro.mCentre = [&](const detail::Cdt2d::Tri& rT) -> detail::Cdt2d::Vec2 {
            // A fold's lifted circumcentre means nothing: its centroid in
            // the plane, always inside it, splits it.
            if (folded(rT)) {
                const auto& xa = cdt.Point(rT[0]);
                const auto& xb = cdt.Point(rT[1]);
                const auto& xc = cdt.Point(rT[2]);
                return {(xa[0] + xb[0] + xc[0]) / 3.0, (xa[1] + xb[1] + xc[1]) / 3.0};
            }
            const Vec3& pa = pos(mesh_of.at(rT[0]));
            const Vec3& pb = pos(mesh_of.at(rT[1]));
            const Vec3& pc = pos(mesh_of.at(rT[2]));
            const Vec3 e1{pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2]};
            const Vec3 e2{pc[0] - pa[0], pc[1] - pa[1], pc[2] - pa[2]};
            const double g11 = e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2];
            const double g12 = e1[0] * e2[0] + e1[1] * e2[1] + e1[2] * e2[2];
            const double g22 = e2[0] * e2[0] + e2[1] * e2[1] + e2[2] * e2[2];
            const double det = g11 * g22 - g12 * g12;
            if (!(det > 0.0))
                return {std::numeric_limits<double>::quiet_NaN(), 0.0};
            // (O - A) . e_i = |e_i|^2 / 2 for O = A + al e1 + be e2.
            const double al = 0.5 * g22 * (g11 - g12) / det;
            const double be = 0.5 * g11 * (g22 - g12) / det;
            const auto& xa = cdt.Point(rT[0]);
            const auto& xb = cdt.Point(rT[1]);
            const auto& xc = cdt.Point(rT[2]);
            return {xa[0] + al * (xb[0] - xa[0]) + be * (xc[0] - xa[0]),
                    xa[1] + al * (xb[1] - xa[1]) + be * (xc[1] - xa[1])};
        };
        // Ruppert's insertion radii, which bound how close refinement puts
        // points (Ruppert 1995; Shewchuk 2002): with B the ratio bound and l
        // a skinny triangle's shortest edge, its circumcentre lies at least
        // its circumradius (at least B l) from every vertex of its cavity,
        // and the midpoint of a subsegment it encroaches at least B l /
        // sqrt(2) from that subsegment's ends. In a plane that is the
        // surface these always hold, and on a curved one nearly; where the
        // plane's Delaunay cavity is not the surface's (one affine map
        // cannot follow a metric that varies too much, as near a pole) a
        // point can land much closer: below half the radius it is declined,
        // so points do not pile up there, and the triangle is reported.
        const double slack = 0.5;
        // Where the surface stretches the plane, no point closer than h / 4
        // (see `stretched`), wherever its triangle is: refinement there does
        // not make the triangles it would then have to leave.
        const auto stretched_floor = [&](const detail::Cdt2d::Vec2& rX) {
            return !varying && stretched_at(rX) ? 0.25 * h : 0.0;
        };
        // A triangle bad by its size alone keeps its new point a third of
        // the local size from every vertex of its cavity, which its
        // circumcentre, at least half its longest edge from its corners,
        // meets.
        const auto size_only = [&](const detail::Cdt2d::Tri& rT) {
            if (!p_refine_size)
                return false;
            std::array<std::int64_t, 3> ids{};
            for (int k = 0; k < 3; ++k)
                ids[static_cast<std::size_t>(k)] = mesh_of.at(rT[static_cast<std::size_t>(k)]);
            return oversize(ids) > 0.0 &&
                   !(mcad_ratio(pos(ids[0]), pos(ids[1]), pos(ids[2])) > ratio_bound);
        };
        ro.mAccept = [&](const detail::Cdt2d::Tri& rT, const detail::Cdt2d::Vec2& rC,
                         const std::vector<std::int64_t>& rAround) {
            const Vec3& pa = pos(mesh_of.at(rT[0]));
            const Vec3& pb = pos(mesh_of.at(rT[1]));
            const Vec3& pc = pos(mesh_of.at(rT[2]));
            const double shortest =
                std::min({mcad_dist(pa, pb), mcad_dist(pb, pc), mcad_dist(pc, pa)});
            const Vec3 x = lift(rC);
            // A fold's centroid keeps a quarter of its shortest edge from
            // the cavity's vertices.
            const double radius =
                folded(rT)      ? std::max(0.25 * shortest, stretched_floor(rC))
                : size_only(rT) ? std::max(p_refine_size->At(x) / 3.0, stretched_floor(rC))
                                : std::max(slack * ratio_bound * shortest, stretched_floor(rC));
            for (const std::int64_t w : rAround) {
                const auto it = mesh_of.find(w);
                if (it != mesh_of.end() && mcad_dist(pos(it->second), x) < radius)
                    return false;
            }
            return true;
        };
        ro.mAcceptSplit = [&](std::int64_t A, std::int64_t B, const detail::Cdt2d::Tri& rT) {
            // A fold whose centroid encroaches a boundary piece splits it,
            // down to a quarter of the fold's shortest edge.
            if (folded(rT)) {
                const Vec3& pa = pos(mesh_of.at(rT[0]));
                const Vec3& pb = pos(mesh_of.at(rT[1]));
                const Vec3& pc = pos(mesh_of.at(rT[2]));
                const double shortest =
                    std::min({mcad_dist(pa, pb), mcad_dist(pb, pc), mcad_dist(pc, pa)});
                const auto& xa = cdt.Point(A);
                const auto& xb = cdt.Point(B);
                const detail::Cdt2d::Vec2 mid{0.5 * (xa[0] + xb[0]), 0.5 * (xa[1] + xb[1])};
                return 0.5 * mcad_dist(pos(mesh_of.at(A)), pos(mesh_of.at(B))) >=
                       std::max(0.25 * shortest, stretched_floor(mid));
            }
            // A boundary piece is sized by the field already: a triangle bad
            // by its size alone asks to split one only when the piece is
            // itself half again too long.
            if (size_only(rT)) {
                const Vec3& r_a = pos(mesh_of.at(A));
                const Vec3& r_b = pos(mesh_of.at(B));
                const Vec3 mid{0.5 * (r_a[0] + r_b[0]), 0.5 * (r_a[1] + r_b[1]),
                               0.5 * (r_a[2] + r_b[2])};
                return mcad_dist(r_a, r_b) > detail::kMeshCadSizeSlack * p_refine_size->At(mid);
            }
            const Vec3& pa = pos(mesh_of.at(rT[0]));
            const Vec3& pb = pos(mesh_of.at(rT[1]));
            const Vec3& pc = pos(mesh_of.at(rT[2]));
            const double shortest =
                std::min({mcad_dist(pa, pb), mcad_dist(pb, pc), mcad_dist(pc, pa)});
            const auto& xa = cdt.Point(A);
            const auto& xb = cdt.Point(B);
            const detail::Cdt2d::Vec2 mid{0.5 * (xa[0] + xb[0]), 0.5 * (xa[1] + xb[1])};
            return 0.5 * mcad_dist(pos(mesh_of.at(A)), pos(mesh_of.at(B))) >=
                   std::max(slack * ratio_bound * shortest / std::sqrt(2.0), stretched_floor(mid));
        };
        ro.mMinAngle = [&](const detail::Cdt2d::Tri& rT) {
            std::array<std::int64_t, 3> ids{};
            for (int k = 0; k < 3; ++k)
                ids[static_cast<std::size_t>(k)] = mesh_of.at(rT[static_cast<std::size_t>(k)]);
            return mcad_min_angle(pos(ids[0]), pos(ids[1]), pos(ids[2]));
        };
        // Points enter the map as they are inserted, so the measures can
        // lift them: refinement inserts, then asks the ratio of its star.
        std::int64_t next_lift = cdt.NumVertices();
        const auto catch_up = [&]() {
            for (; next_lift < cdt.NumVertices(); ++next_lift) {
                mesh_of[next_lift] = -(1 + static_cast<std::int64_t>(rF.mInterior.size()));
                rF.mInterior.push_back(lift(cdt.Point(next_lift)));
            }
        };
        const auto ratio = ro.mRatio;
        ro.mRatio = [&, ratio](const detail::Cdt2d::Tri& rT) {
            catch_up();
            return ratio(rT);
        };
        const auto min_angle = ro.mMinAngle;
        ro.mMinAngle = [&, min_angle](const detail::Cdt2d::Tri& rT) {
            catch_up();
            return min_angle(rT);
        };
        // Corner lopping (6.2.3): the input points are the ends of the
        // face's edges and any acute kink inside one (a pole's aside, which
        // has its own handling), their balls measured on the
        // surface, in 3-D, and carried from the earlier meshings, so a ball
        // only shrinks; each flip cycle carries them to the next.
        std::map<std::int64_t, double> lop_radius;  // cdt vertex -> radius
        if (rQ.mLopCorners) {
            std::set<std::int64_t> ends;
            const auto add_ends = [&](std::size_t U) {
                const auto& r_ids = rUses[U].mIds;
                if (!r_ids.empty()) {
                    ends.insert(r_ids.front());
                    ends.insert(r_ids.back());
                }
            };
            for (const auto& r_loop : rF.mLoops)
                for (const std::size_t u : r_loop)
                    add_ends(u);
            for (const std::size_t u : rF.mEmbedded)
                add_ends(u);
            for (const auto& [v, id] : mesh_of) {
                if (id >= 0 && pole_near(id) >= 0)
                    ro.mLopSkip.push_back(v);
                else if (id >= 0 && ends.count(id) != 0)
                    ro.mCorners.push_back(v);
            }
            std::sort(ro.mCorners.begin(), ro.mCorners.end());
            std::sort(ro.mLopSkip.begin(), ro.mLopSkip.end());
            // No corner but poles: nothing to lop (an empty list would make
            // every boundary sample an input point).
            ro.mLopCorners = !ro.mCorners.empty();
            if (rQ.pLopRadius)
                for (const std::int64_t v : ro.mCorners) {
                    const auto it = rQ.pLopRadius->find(pos(mesh_of.at(v)));
                    if (it != rQ.pLopRadius->end())
                        lop_radius[v] = it->second;
                }
            ro.mDistance = [&](std::int64_t P, std::int64_t V, const detail::Cdt2d::Vec2& rX) {
                catch_up();
                return mcad_dist(pos(mesh_of.at(P)), V >= 0 ? pos(mesh_of.at(V)) : lift(rX));
            };
        }
        std::set<detail::Cdt2dSplitRequest> asked;
        for (int cycle = 0; cycle < 8; ++cycle) {
            ro.mLopRadius.assign(lop_radius.begin(), lop_radius.end());
            refined = cdt.Refine(ro);
            catch_up();
            for (const auto& r_l : refined.mLopped)
                lop_radius[r_l.mVertex] = r_l.mRadius;
            asked.insert(refined.mRequests.begin(), refined.mRequests.end());
            r_rec.mNumRefinementPoints += refined.mNumCircumcentres;
            rF.mRefinedSoFar += refined.mNumCircumcentres;
            ro.mMaxPoints -= refined.mNumCircumcentres;
            if (refined.mBudgetReached && refined.mNumCircumcentres > 0)
                rF.mRefineBudgetReached = true;
            if (ro.mMaxPoints <= 0 || surface_flips() == 0)
                break;
        }
        refined.mRequests.assign(asked.begin(), asked.end());
        for (const auto& r_q : refined.mRequests) {
            const auto it = origin.find(seg_key(r_q.mA, r_q.mB));
            if (it == origin.end() || pole_near(mesh_of.at(r_q.mA)) >= 0 ||
                pole_near(mesh_of.at(r_q.mB)) >= 0)
                continue;
            const Origin& r_o = it->second;
            const McadUse& r_u = rUses[r_o.mUse];
            const std::size_t n = r_u.mParams.size();
            const std::size_t piece = r_u.mSame ? r_o.mK : n - 2 - r_o.mK;
            std::size_t end = 0;
            if (r_q.mApex >= 0) {
                // The apex at the use's sample k is its piece's start when
                // the use runs with the master.
                const bool at_k = r_o.mStart == r_q.mApex;
                end = (at_k == r_u.mSame) ? 1 : 2;
            }
            if (r_q.mDistance > 0.0) {
                // At the ball's last radius: balls only shrink, so a piece
                // that reached past it at the ask still does.
                const auto it_r = lop_radius.find(r_q.mApex);
                rF.mLopAsks.push_back({static_cast<std::size_t>(r_u.mEdge), piece, end,
                                       it_r != lop_radius.end() ? it_r->second : r_q.mDistance});
                continue;
            }
            rF.mAsks.push_back({static_cast<std::size_t>(r_u.mEdge), piece, end});
        }
        for (const auto& [r_t, apex] : refined.mForced) {
            std::array<std::int64_t, 3> ids{};
            for (int k = 0; k < 3; ++k)
                ids[static_cast<std::size_t>(k)] = mesh_of.at(r_t[static_cast<std::size_t>(k)]);
            std::sort(ids.begin(), ids.end());
            forced.insert(ids);
            (void)apex;
        }
        for (const auto& r_a : refined.mSmallAngles) {
            const std::int64_t id = mesh_of.at(r_a.mVertex);
            corner_angle[id] = r_a.mAngleDeg;
            corner_x[id] = cdt.Point(r_a.mVertex);
        }
    }

    // Guarded smoothing (6.2.5): each interior point toward its star's
    // centroid (Laplacian) or area-weighted circumcentre (ODT, Chen and Xu
    // 2004), in the plane, kept only when no triangle of its star inverts
    // (MoveVertex's exact orient2d) and the star's smallest angle on the
    // surface does not drop (Field 1988). Lawson's flips follow, judged on
    // the surface (the larger smallest angle), so a pass never undoes what
    // the guard kept. Boundary points stay, so shared edges stay shared; a
    // star within h of a pole is left.
    const std::int64_t smooth_passes =
        rQ.mDeferSmoothing && !(rF.mAsks.empty() && rF.mLopAsks.empty()) ? 0 : rQ.mSmoothIterations;
    rF.mSmoothingDeferred = smooth_passes < rQ.mSmoothIterations;
    for (std::int64_t pass = 0; pass < smooth_passes; ++pass) {
        for (std::int64_t v = detail::Cdt2d::kFirstVertex; v < cdt.NumVertices(); ++v) {
            const auto it = mesh_of.find(v);
            if (it == mesh_of.end() || it->second >= 0 || cdt.OnConstraint(v))
                continue;
            const std::vector<detail::Cdt2d::Tri> star = cdt.Star(v);
            bool skip = false;
            double before = 180.0;
            for (const auto& r_t : star) {
                std::array<std::int64_t, 3> ids{};
                if (exempt(r_t, ids)) {
                    skip = true;
                    break;
                }
                before = std::min(before, mcad_min_angle(pos(ids[0]), pos(ids[1]), pos(ids[2])));
            }
            if (skip || star.empty())
                continue;
            detail::Cdt2d::Vec2 target{0.0, 0.0};
            if (rQ.mSmoothMethod == MeshCadSmoothing::Laplacian) {
                const std::vector<std::int64_t> nbrs = cdt.Neighbours(v);
                for (const std::int64_t w : nbrs) {
                    target[0] += cdt.Point(w)[0];
                    target[1] += cdt.Point(w)[1];
                }
                target[0] /= double(nbrs.size());
                target[1] /= double(nbrs.size());
            } else {
                double total = 0.0;
                for (const auto& r_t : star) {
                    const auto& a = cdt.Point(r_t[0]);
                    const auto& b = cdt.Point(r_t[1]);
                    const auto& c = cdt.Point(r_t[2]);
                    const double bx = b[0] - a[0], by = b[1] - a[1];
                    const double cx = c[0] - a[0], cy = c[1] - a[1];
                    const double d = 2.0 * (bx * cy - by * cx);
                    if (!(d > 0.0))
                        continue;
                    const double b2 = bx * bx + by * by;
                    const double c2 = cx * cx + cy * cy;
                    const double area = 0.25 * d;
                    target[0] += area * (a[0] + (cy * b2 - by * c2) / d);
                    target[1] += area * (a[1] + (bx * c2 - cx * b2) / d);
                    total += area;
                }
                if (!(total > 0.0))
                    continue;
                target[0] /= total;
                target[1] /= total;
            }
            const Vec3 moved = lift(target);
            const std::int64_t self = it->second;
            double after = 180.0;
            for (const auto& r_t : star) {
                Vec3 p[3];
                for (int k = 0; k < 3; ++k)
                    p[k] = r_t[static_cast<std::size_t>(k)] == v
                               ? moved
                               : pos(mesh_of.at(r_t[static_cast<std::size_t>(k)]));
                after = std::min(after, mcad_min_angle(p[0], p[1], p[2]));
            }
            if (!(after >= before) || !cdt.MoveVertex(v, target, false))
                continue;
            rF.mInterior[static_cast<std::size_t>(-self - 1)] = moved;
            ++r_rec.mNumSmoothed;
        }
        sides_upto = -1;  // points moved
        surface_flips();
    }

    // A fold refinement could not remove fails the face as degenerate
    // rather than being written folded. The boundary splits it asked for
    // are kept: the next round meshes it again on them, and only a fold
    // left after the last round (one past the face's points, or that no
    // split reaches) stays a failure.
    {
        std::int64_t folds = 0;
        for (const auto& r_t : cdt.InsideTriangles())
            if (folded(r_t))
                ++folds;
        if (folds > 0) {
            std::vector<std::array<std::size_t, 3>> asks = std::move(rF.mAsks);
            std::vector<McadLopAsk> lop_asks = std::move(rF.mLopAsks);
            fail(MeshCadFaceStatus::Degenerate,
                 std::to_string(folds) +
                     " triangles fold on the surface (their normals turn more than 120 "
                     "degrees from the surface's)");
            rF.mAsks = std::move(asks);
            rF.mLopAsks = std::move(lop_asks);
            rF.mFolded = true;
            return;
        }
    }

    // Lift, drop what a pole collapses, orient along the face's normal.
    double min_angle = 180.0;
    double max_angle = 0.0;
    double min_edge = std::numeric_limits<double>::infinity();
    double max_edge = 0.0;
    // The corners that keep triangles under the bound: refinement's small
    // input angles and the poles; the smallest angle each one forces, and
    // that of the rest (the bound the face reached).
    std::map<std::int64_t, double> pole_min;
    double bound = 180.0;
    // Curved elements (6.6): an interior edge's nodes are the surface at the
    // parameter points between its ends (1/2, or 1/3 and 2/3), kept once per
    // edge (by its plane vertices, in key order) for both its triangles; the
    // plane is a linear image of the parameters, so the fractions are the
    // parameters' too.
    // At a pole the parameters collapse: its plane vertices lie along a
    // degenerate segment, each with a parameter that means nothing on the
    // surface. An edge or element there takes the pole's parameters as the
    // point of that segment nearest the other end's (on a sphere, that end's
    // meridian), or its parameter midpoint would twist the element.
    const int order = rO.mOrder;
    std::map<std::int64_t, std::array<Vec2, 2>> pole_segment;
    if (order > 1) {
        std::map<std::int64_t, std::vector<Vec2>> at_pole;
        for (std::int64_t v = detail::Cdt2d::kFirstVertex; v < cdt.NumVertices(); ++v) {
            const auto it = mesh_of.find(v);
            if (it != mesh_of.end() && pole(it->second))
                at_pole[it->second].push_back(pl.ToUv(cdt.Point(v)));
        }
        for (const auto& [id, r_uv] : at_pole) {
            std::array<Vec2, 2> seg{r_uv.front(), r_uv.front()};
            double far = -1.0;
            for (std::size_t i = 0; i < r_uv.size(); ++i)
                for (std::size_t j = i + 1; j < r_uv.size(); ++j) {
                    const double d = std::hypot(r_uv[i][0] - r_uv[j][0], r_uv[i][1] - r_uv[j][1]);
                    if (d > far) {
                        far = d;
                        seg = {r_uv[i], r_uv[j]};
                    }
                }
            pole_segment.emplace(id, seg);
        }
    }
    const auto param_of = [&](std::int64_t V) { return pl.ToUv(cdt.Point(V)); };
    // The parameters of plane vertex `V`, seen from parameters `rToward`.
    const auto param_toward = [&](std::int64_t V, const Vec2& rToward) {
        const auto it = pole_segment.find(mesh_of.at(V));
        if (it == pole_segment.end())
            return param_of(V);
        const Vec2& a = it->second[0];
        const Vec2& b = it->second[1];
        const double dx = b[0] - a[0], dy = b[1] - a[1];
        const double len2 = dx * dx + dy * dy;
        if (!(len2 > 0.0))
            return a;
        const double t =
            std::clamp(((rToward[0] - a[0]) * dx + (rToward[1] - a[1]) * dy) / len2, 0.0, 1.0);
        return Vec2{a[0] + t * dx, a[1] + t * dy};
    };
    const auto surface_at = [&](const Vec2& rUv) {
        Vec3 x{};
        r_s.Eval(rUv[0], rUv[1], x.data(), nullptr, nullptr);
        return x;
    };
    std::map<std::pair<std::int64_t, std::int64_t>, std::array<std::int64_t, 2>> edge_nodes;
    const auto interior_nodes = [&](std::int64_t A, std::int64_t B) {
        const auto key = seg_key(A, B);
        auto it = edge_nodes.find(key);
        if (it == edge_nodes.end()) {
            std::array<std::int64_t, 2> ids{-1, -1};
            const Vec2 ua = param_toward(key.first, param_of(key.second));
            const Vec2 ub = param_toward(key.second, param_of(key.first));
            for (int j = 1; j < order; ++j) {
                const double f = double(j) / order;
                ids[static_cast<std::size_t>(j - 1)] =
                    static_cast<std::int64_t>(rF.mHighOrder.size());
                rF.mHighOrder.push_back(
                    surface_at({ua[0] + f * (ub[0] - ua[0]), ua[1] + f * (ub[1] - ua[1])}));
            }
            it = edge_nodes.emplace(key, ids).first;
        }
        std::array<std::int64_t, 2> ids = it->second;
        if (A != key.first && order == 3)
            std::swap(ids[0], ids[1]);
        return ids;
    };
    for (const auto& r_t : cdt.InsideTriangles()) {
        std::array<std::int64_t, 3> tri{};
        for (int k = 0; k < 3; ++k)
            tri[k] = mesh_of.at(r_t[k]);
        if (tri[0] == tri[1] || tri[1] == tri[2] || tri[0] == tri[2]) {
            ++r_rec.mNumDegenerateDropped;
            continue;
        }
        {
            const auto ang = mcad_angles(pos(tri[0]), pos(tri[1]), pos(tri[2]));
            for (const double a : ang)
                ++rF.mHistogram[std::min<std::size_t>(kMeshCadAngleBins - 1,
                                                      static_cast<std::size_t>(a / 10.0))];
            const double m = std::min({ang[0], ang[1], ang[2]});
            std::array<std::int64_t, 3> key = tri;
            std::sort(key.begin(), key.end());
            bool at_corner = forced.count(key) > 0;
            for (const std::int64_t id : tri) {
                const std::int64_t p = pole_near(id);
                if (p >= 0) {
                    auto [it, fresh] = pole_min.emplace(p, m);
                    if (!fresh)
                        it->second = std::min(it->second, m);
                    at_corner = true;
                }
                at_corner = at_corner || corner_angle.count(id) > 0;
            }
            if (!at_corner && !varying && rQ.mMinAngleDeg > 0.0 && m < rQ.mMinAngleDeg &&
                stretched(r_t)) {
                at_corner = true;
                ++rF.mNumStretched;
            }
            if (!at_corner)
                bound = std::min(bound, m);
            if (rQ.mMinAngleDeg > 0.0 && m < rQ.mMinAngleDeg)
                r_rec.mBelowBound = true;
        }
        std::array<std::int64_t, 3> cv = r_t;
        if (rBrepFace.mSwappedNormal) {
            std::swap(tri[1], tri[2]);
            std::swap(cv[1], cv[2]);
        }
        if (order > 1) {
            std::array<McadFace::CurvedEdge, 3> ce{};
            for (std::size_t k = 0; k < 3; ++k) {
                const std::int64_t a = cv[k];
                const std::int64_t b = cv[(k + 1) % 3];
                const auto it = origin.find(seg_key(a, b));
                if (it != origin.end()) {
                    // A boundary piece: the triangle runs from the use's
                    // sample k to k + 1 when it starts at mStart, and the
                    // use runs with the master when mSame.
                    const McadUse& r_u = rUses[it->second.mUse];
                    const std::size_t n = r_u.mParams.size();
                    ce[k].mEdge = static_cast<std::int64_t>(r_u.mEdge);
                    ce[k].mPiece = static_cast<std::int64_t>(r_u.mSame ? it->second.mK
                                                                       : n - 2 - it->second.mK);
                    ce[k].mReversed = (a == it->second.mStart) != r_u.mSame;
                } else {
                    ce[k].mNodes = interior_nodes(a, b);
                }
            }
            // The centroid's parameters, a pole corner taken toward the
            // other two.
            Vec2 uv{0.0, 0.0};
            for (std::size_t k = 0; k < 3; ++k) {
                const Vec2 ob = param_of(cv[(k + 1) % 3]);
                const Vec2 oc = param_of(cv[(k + 2) % 3]);
                const Vec2 c = param_toward(cv[k], {0.5 * (ob[0] + oc[0]), 0.5 * (ob[1] + oc[1])});
                uv[0] += c[0] / 3.0;
                uv[1] += c[1] / 3.0;
            }
            Vec3 xs{}, su{}, sv{};
            r_s.Eval(uv[0], uv[1], xs.data(), su.data(), sv.data());
            Vec3 nrm{su[1] * sv[2] - su[2] * sv[1], su[2] * sv[0] - su[0] * sv[2],
                     su[0] * sv[1] - su[1] * sv[0]};
            if (rBrepFace.mSwappedNormal)
                nrm = {-nrm[0], -nrm[1], -nrm[2]};
            rF.mCurvedEdges.push_back(ce);
            rF.mCurvedNormal.push_back(nrm);
            if (order == 3) {
                rF.mCurvedCentre.push_back(static_cast<std::int64_t>(rF.mHighOrder.size()));
                rF.mHighOrder.push_back(xs);
            }
        }
        for (int k = 0; k < 3; ++k) {
            const Vec3& o = pos(tri[k]);
            const Vec3& p = pos(tri[(k + 1) % 3]);
            const Vec3& q = pos(tri[(k + 2) % 3]);
            const Vec3 a{p[0] - o[0], p[1] - o[1], p[2] - o[2]};
            const Vec3 b{q[0] - o[0], q[1] - o[1], q[2] - o[2]};
            const Vec3 c{a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
                         a[0] * b[1] - a[1] * b[0]};
            const double ang = std::atan2(std::sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]),
                                          a[0] * b[0] + a[1] * b[1] + a[2] * b[2]) *
                               (180.0 / 3.14159265358979323846);
            min_angle = std::min(min_angle, ang);
            max_angle = std::max(max_angle, ang);
            const double e = mcad_dist(o, p);
            min_edge = std::min(min_edge, e);
            max_edge = std::max(max_edge, e);
        }
        rF.mTriangles.push_back(tri);
        // Merged faces (6.5.2): which face's coarse triangle holds it.
        if (pl.pDisc && pl.pDisc->mSurfaces.size() > 1) {
            const auto& xa = cdt.Point(r_t[0]);
            const auto& xb = cdt.Point(r_t[1]);
            const auto& xc = cdt.Point(r_t[2]);
            rF.mTriOwner.push_back(pl.pDisc
                                       ->Locate(pl.ToUv({(xa[0] + xb[0] + xc[0]) / 3.0,
                                                         (xa[1] + xb[1] + xc[1]) / 3.0}))
                                       .mOwner);
        }
    }
    if (rF.mTriangles.empty()) {
        // A face smaller than h whose loops are poles and a seam (a small
        // surface of revolution, roadmap 6.5.7) has no lattice point, and
        // every triangle of its boundary alone touches a pole or the seam
        // twice: the driver meshes it once more with an interior grid.
        rF.mCollapsed = true;
        return fail(MeshCadFaceStatus::Empty, "every triangle collapses on the surface");
    }
    r_rec.mNumTriangles = static_cast<std::int64_t>(rF.mTriangles.size());
    r_rec.mNumInteriorPoints = static_cast<std::int64_t>(rF.mInterior.size());
    r_rec.mMinAngleDeg = min_angle;
    r_rec.mMaxAngleDeg = max_angle;
    r_rec.mMinEdge = min_edge;
    r_rec.mMaxEdge = max_edge;
    r_rec.mBoundDeg = bound;
    if (rQ.mCollectSizes) {
        // Every point of the face once, in vertex order, with the surface's
        // outward normal and largest principal curvature where it sits.
        std::vector<char> seen(static_cast<std::size_t>(cdt.NumVertices()), 0);
        for (const auto& r_t : cdt.InsideTriangles())
            for (const std::int64_t v : r_t) {
                if (seen[static_cast<std::size_t>(v)])
                    continue;
                seen[static_cast<std::size_t>(v)] = 1;
                const std::int64_t id = mesh_of.at(v);
                Vec3 x{}, su{}, sv{};
                surface_x(cdt.Point(v), x.data(), su.data(), sv.data());
                detail::MeshCadSizeSample smp;
                smp.mX = pos(id);
                smp.mFace = r_rec.mFace;
                Vec3 n{su[1] * sv[2] - su[2] * sv[1], su[2] * sv[0] - su[0] * sv[2],
                       su[0] * sv[1] - su[1] * sv[0]};
                const double len = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
                if (len > 0.0)
                    for (double& r_c : n)
                        r_c *= (rBrepFace.mSwappedNormal ? -1.0 : 1.0) / len;
                smp.mNormal = len > 0.0 ? n : Vec3{0.0, 0.0, 0.0};
                const auto [p_surface, uv] = params_x(cdt.Point(v));
                smp.mCurvature = p_surface->MaxCurvature(uv[0], uv[1]);
                rF.mSizeSamples.push_back(smp);
            }
    }
    if (!refined.mLopped.empty()) {
        // Corner lopping's corners (6.2.3), with what Pav and Walkington
        // (2005) state of them: Corollary 3's bound for the most any
        // boundary piece of the face turns (sampled tangents, in 3-D), on a
        // planar face, and whether the run is under it, finished, and kept
        // every triangle under min_angle within twice a ball's radius of its
        // corner (Lemma 7).
        double turning = 0.0;
        const auto add_turning = [&](std::size_t U) {
            const auto& r_tr = rTraces[U];
            const auto& r_p = rUses[U].mParams;
            if (!r_tr)
                return;
            for (std::size_t k = 0; k + 1 < r_p.size(); ++k)
                turning = std::max(turning, mcad_piece_turning(*r_tr, r_p[k], r_p[k + 1]));
        };
        for (const auto& r_loop : rF.mLoops)
            for (const std::size_t u : r_loop)
                add_turning(u);
        for (const std::size_t u : rF.mEmbedded)
            add_turning(u);
        // Planar: no principal curvature over a 7 by 7 grid of the patch
        // past one in a billion target sizes.
        bool planar = true;
        if (rF.pSurface) {
            const auto& r_box = rF.pSurface->Box();
            for (int i = 0; i <= 6 && planar; ++i)
                for (int j = 0; j <= 6 && planar; ++j)
                    planar =
                        !(rF.pSurface->MaxCurvature(r_box[0] + (r_box[1] - r_box[0]) * i / 6.0,
                                                    r_box[2] + (r_box[3] - r_box[2]) * j / 6.0) *
                              h >
                          1e-9);
        }
        const double bound = planar ? detail::cdt2d_lop_bound_deg(turning) : 0.0;
        for (const auto& r_l : refined.mLopped) {
            MeshCadLoppedCorner c;
            c.mFace = r_rec.mFace;
            c.mPoint = mesh_of.at(r_l.mVertex);
            c.mAngleDeg = r_l.mAngleDeg;
            c.mLfs0 = r_l.mLfs0;
            c.mLopRadius = r_l.mRadius;
            c.mTurningDeg = turning * (180.0 / kMcadPi);
            c.mStatedBoundDeg = bound;
            c.mNumForced = r_l.mNumForced;
            c.mMinForcedDeg = r_l.mMinForcedDeg;
            c.mWithinTwoRadii = r_l.mWithinTwoRadii;
            if (!planar)
                c.mReason = "the face is curved; Corollary 3 is proved for a planar domain";
            else if (!(bound > 0.0))
                c.mReason = "a boundary piece turns more than 60 degrees";
            else if (!(rQ.mMinAngleDeg < bound))
                c.mReason = "min_angle is not under the bound for the face's turning";
            else if (rF.mRefineBudgetReached || refined.mBudgetReached ||
                     refined.mNumUnresolved > 0)
                c.mReason = "refinement did not finish";
            else if (!r_l.mWithinTwoRadii)
                c.mReason = "a triangle under min_angle reaches past twice the ball's radius";
            c.mGuaranteed = c.mReason.empty();
            rF.mLopped.emplace_back(pos(c.mPoint), c);
            if (rQ.mMinAngleDeg > 0.0 && r_l.mAngleDeg < 60.0) {
                // Still a small input angle, as concentric shells report one.
                corner_angle[c.mPoint] = r_l.mAngleDeg;
                corner_x[c.mPoint] = cdt.Point(r_l.mVertex);
            }
        }
    }
    if (rQ.mMinAngleDeg > 0.0) {
        // The corners, in mesh-point order: small input angles with the
        // smallest angle of the triangles they force on the surface, and
        // the poles the face's triangles touch.
        std::map<std::int64_t, double> forced_min;
        for (const auto& [r_t, apex] : refined.mForced) {
            std::array<std::int64_t, 3> ids{};
            for (int k = 0; k < 3; ++k)
                ids[static_cast<std::size_t>(k)] = mesh_of.at(r_t[static_cast<std::size_t>(k)]);
            const double m = mcad_min_angle(pos(ids[0]), pos(ids[1]), pos(ids[2]));
            auto [it, fresh] = forced_min.emplace(mesh_of.at(apex), m);
            if (!fresh)
                it->second = std::min(it->second, m);
        }
        for (const auto& [id, deg] : corner_angle) {
            MeshCadCorner c;
            c.mFace = r_rec.mFace;
            c.mPoint = id;
            const Vec2 uv = params_x(corner_x.at(id)).second;
            c.mU = uv[0];
            c.mV = uv[1];
            c.mX = pos(id);
            c.mAngleDeg = deg;
            const auto it = forced_min.find(id);
            if (it != forced_min.end())
                c.mMinForcedDeg = it->second;
            r_rec.mMinInputAngleDeg = std::min(r_rec.mMinInputAngleDeg, deg);
            r_rec.mSmallAngleCorners.push_back(c);
        }
        for (const auto& [id, m] : pole_min) {
            MeshCadCorner c;
            c.mFace = r_rec.mFace;
            c.mPoint = id;
            for (const auto& r_poly : loops)
                for (const McadEntry& r_e : r_poly)
                    if (r_e.mId == id) {
                        c.mU = r_e.mUv[0];
                        c.mV = r_e.mUv[1];
                    }
            c.mX = pos(id);
            c.mMinForcedDeg = m;
            c.mReach = pole_reach;
            c.mPole = true;
            r_rec.mSmallAngleCorners.push_back(c);
        }
        std::sort(
            r_rec.mSmallAngleCorners.begin(), r_rec.mSmallAngleCorners.end(),
            [](const MeshCadCorner& rA, const MeshCadCorner& rB) { return rA.mPoint < rB.mPoint; });
    }
}

}  // namespace

const char* mesh_cad_face_status_name(MeshCadFaceStatus Status) {
    switch (Status) {
        case MeshCadFaceStatus::Ok:
            return "ok";
        case MeshCadFaceStatus::Empty:
            return "empty";
        case MeshCadFaceStatus::OpenLoop:
            return "open_loop";
        case MeshCadFaceStatus::SelfIntersecting:
            return "self_intersecting";
        case MeshCadFaceStatus::Degenerate:
            return "degenerate";
        case MeshCadFaceStatus::ProjectionFailed:
            return "projection_failed";
        case MeshCadFaceStatus::BudgetExceeded:
            return "budget_exceeded";
        case MeshCadFaceStatus::Unsupported:
            return "unsupported";
        case MeshCadFaceStatus::Other:
            return "other";
    }
    return "other";
}

void detail::mesh_cad_set_face_order_seed(std::uint64_t Seed) {
    g_mcad_face_order_seed.store(Seed);
}

const char* mesh_cad_smoothing_name(MeshCadSmoothing Method) {
    return Method == MeshCadSmoothing::Odt ? "odt" : "laplacian";
}

namespace {

/// One meshing of the model: with @p pSize, sized by that field; with
/// @p pSamples, the size samples of every face appended in face order.
/// The given sizes' refusals (6.3.1).
void mcad_check_sources(const BrepInfo& rBrep, const MeshCadSurfaceOptions& rOptions,
                        const MeshCadSizeSources& rSources) {
    if (rSources.mVirtualEdges && rOptions.mOrder > 1)
        mcad_fail(
            "virtual_edges: curved elements (order 2 or 3) follow each edge's own "
            "curve, so a chord across a suppressed corner cannot be curved; use order 1");
    const auto size_ok = [](double S) { return std::isfinite(S) && S > 0.0; };
    for (const auto& [f, s] : rSources.mFaceSizes) {
        if (f < 0 || static_cast<std::size_t>(f) >= rBrep.mFaces.size())
            mcad_fail("face_sizes: face " + std::to_string(f) + " is not a face of the B-rep");
        if (!size_ok(s))
            mcad_fail("face_sizes: the size of face " + std::to_string(f) +
                      " must be finite and positive");
    }
    for (const auto& [e, s] : rSources.mEdgeSizes) {
        if (e < 0 || static_cast<std::size_t>(e) >= rBrep.mEdges.size())
            mcad_fail("edge_sizes: edge " + std::to_string(e) + " is not an edge of the B-rep");
        if (!size_ok(s))
            mcad_fail("edge_sizes: the size of edge " + std::to_string(e) +
                      " must be finite and positive");
    }
    if (rSources.pBackground) {
        const Mesh& r_bg = *rSources.pBackground;
        if (!r_bg.HasPointData(rSources.mBackgroundArray))
            mcad_fail("background: it has no point_data '" + rSources.mBackgroundArray + "'");
        const NDArray& r_a = r_bg.PointData(rSources.mBackgroundArray);
        if (r_a.Size() != r_bg.NumPoints())
            mcad_fail("background: its point_data '" + rSources.mBackgroundArray +
                      "' must have one component");
    }
}

MeshCadSurfaceResult mcad_surface_run(const BrepInfo& rBrep, const MeshCadSurfaceOptions& rOptions,
                                      const detail::MeshCadSizeField* pSize,
                                      std::vector<detail::MeshCadSizeSample>* pSamples,
                                      bool VirtualEdges = false,
                                      std::vector<std::vector<std::int64_t>>* pChains = nullptr,
                                      const MeshCadSurfaceExtras* pExtras = nullptr,
                                      MeshCadSurfaceExtrasReport* pExtrasReport = nullptr) {
    if (rBrep.mFaces.empty())
        mcad_fail("the B-rep has no faces");
    const MeshCadSurfaceExtras extras = pExtras ? *pExtras : MeshCadSurfaceExtras{};
    const double h = rOptions.mTargetSize;
    MeshCadSurfaceResult result;
    const auto warn = [&](const std::string& rMsg) {
        log::warn("mesh_cad_surface: {}", rMsg);
        result.mWarnings.push_back(rMsg);
    };

    // The tolerance: relative to the control nets' box unless given.
    double lo[3] = {std::numeric_limits<double>::infinity(),
                    std::numeric_limits<double>::infinity(),
                    std::numeric_limits<double>::infinity()};
    double hi[3] = {-lo[0], -lo[1], -lo[2]};
    for (const NurbsPatch& r_p : rBrep.mPatches) {
        if (r_p.mControlPoints.Dtype() != DType::Float64 ||
            r_p.mControlPoints.Shape().size() != 2 || r_p.mControlPoints.Shape()[1] != 3)
            continue;
        const double* c = r_p.mControlPoints.As<double>();
        for (std::size_t i = 0; i < r_p.mControlPoints.Shape()[0]; ++i)
            for (int k = 0; k < 3; ++k)
                if (std::isfinite(c[3 * i + k])) {
                    lo[k] = std::min(lo[k], c[3 * i + k]);
                    hi[k] = std::max(hi[k], c[3 * i + k]);
                }
    }
    double diag = 0.0;
    if (lo[0] <= hi[0])
        diag = std::sqrt((hi[0] - lo[0]) * (hi[0] - lo[0]) + (hi[1] - lo[1]) * (hi[1] - lo[1]) +
                         (hi[2] - lo[2]) * (hi[2] - lo[2]));
    const double tol = rOptions.mGeometricTolerance > 0.0 ? rOptions.mGeometricTolerance
                                                          : (diag > 0.0 ? 1e-7 * diag : 1e-12);
    result.mTolerance = tol;
    // Two uses of an edge or corner farther apart than this mean the B-rep's
    // topology does not describe its geometry: the faces using them fail.
    const double gap_fail = std::max(0.25 * h, 10.0 * tol);

    // Faces, their surfaces and their trims (one use each).
    const std::size_t nf = rBrep.mFaces.size();
    std::vector<McadFace> faces(nf);
    std::vector<McadUse> uses;
    for (std::size_t f = 0; f < nf; ++f) {
        const BrepFace& r_bf = rBrep.mFaces[f];
        McadFace& r_f = faces[f];
        r_f.mRecord.mFace = static_cast<std::int64_t>(f);
        r_f.mRecord.mId = r_bf.mId;
        const auto unsupported = [&](const std::string& rWhy) {
            r_f.mRecord.mStatus = MeshCadFaceStatus::Unsupported;
            r_f.mRecord.mReason = rWhy;
        };
        if (r_bf.mPatch < 0 || static_cast<std::size_t>(r_bf.mPatch) >= rBrep.mPatches.size()) {
            unsupported("its patch index is out of range");
            continue;
        }
        try {
            r_f.pSurface = std::make_unique<detail::BrepSurfaceEval>(
                rBrep.mPatches[static_cast<std::size_t>(r_bf.mPatch)], r_bf.mRational);
        } catch (const std::invalid_argument& r_e) {
            unsupported(mcad_why(r_e));
            continue;
        }
        const auto add_use = [&](const BrepCurve* pC, bool Forward, std::int64_t TrimIndex) {
            McadUse u;
            u.mFace = f;
            u.pCurve = pC;
            u.mForward = Forward;
            u.mTrimIndex = TrimIndex;
            uses.push_back(std::move(u));
            return uses.size() - 1;
        };
        bool any_trim = false;
        for (const BrepLoop& r_l : r_bf.mLoops)
            any_trim = any_trim || !r_l.mTrims.empty();
        if (!any_trim) {
            // An untrimmed face: its parameter rectangle, counter-clockwise.
            const auto& b = r_f.pSurface->Box();
            r_f.mSynthetic = {mcad_line(b[0], b[2], b[1], b[2]), mcad_line(b[1], b[2], b[1], b[3]),
                              mcad_line(b[1], b[3], b[0], b[3]), mcad_line(b[0], b[3], b[0], b[2])};
            std::vector<std::size_t> loop;
            for (const BrepCurve& r_c : r_f.mSynthetic)
                loop.push_back(add_use(&r_c, true, -1));
            for (const std::size_t u : loop)
                uses[u].mLoopSize = loop.size();
            r_f.mLoops.push_back(std::move(loop));
        } else {
            for (const BrepLoop& r_l : r_bf.mLoops) {
                std::vector<std::size_t> loop;
                for (const std::size_t k : mcad_chain(r_l))
                    loop.push_back(add_use(&r_l.mTrims[k].mCurve, r_l.mTrims[k].mDirection,
                                           r_l.mTrims[k].mTrimIndex));
                for (const std::size_t u : loop)
                    uses[u].mLoopSize = loop.size();
                if (!loop.empty())
                    r_f.mLoops.push_back(std::move(loop));
            }
        }
        // Embedded edges (a coupling or load curve): uses with no loop,
        // inserted as constrained segments once the loops are classified.
        for (const BrepTrim& r_t : r_bf.mEmbeddedEdges) {
            const std::size_t u = add_use(&r_t.mCurve, r_t.mDirection, r_t.mTrimIndex);
            uses[u].mEmbedded = true;
            r_f.mEmbedded.push_back(u);
        }
        r_f.mUsable = true;
    }
    std::vector<std::optional<detail::TrimTrace>> traces(uses.size());
    std::vector<std::string> trace_why(uses.size());
    parallel_for(
        uses.size(),
        [&](std::size_t u) {
            try {
                traces[u].emplace(faces[uses[u].mFace].pSurface.get(),
                                  detail::BrepCurveEval(*uses[u].pCurve), uses[u].mForward);
            } catch (const std::exception& r_e) {
                trace_why[u] = mcad_why(r_e);
            }
        },
        1);
    for (std::size_t u = 0; u < uses.size(); ++u) {
        McadFace& r_f = faces[uses[u].mFace];
        if (trace_why[u].empty())
            continue;
        if (uses[u].mEmbedded) {
            // An invalid embedded curve is left out; the face is meshed.
            r_f.mSkipped.emplace_back(u, trace_why[u]);
        } else if (r_f.mUsable) {
            r_f.mUsable = false;
            r_f.mRecord.mStatus = MeshCadFaceStatus::Unsupported;
            r_f.mRecord.mReason = trace_why[u];
        }
    }
    const auto live = [&](std::size_t U) {
        return faces[uses[U].mFace].mUsable && traces[U].has_value();
    };

    // Edges: the B-rep's own, then the trims none names (paired when their
    // images coincide), each with its uses in ascending order.
    std::unordered_map<std::int64_t, std::size_t> face_by_id;
    std::unordered_map<std::string, std::size_t> face_by_name;
    for (std::size_t f = 0; f < nf; ++f) {
        if (rBrep.mFaces[f].mId >= 0)
            face_by_id.emplace(rBrep.mFaces[f].mId, f);
        if (!rBrep.mFaces[f].mName.empty())
            face_by_name.emplace(rBrep.mFaces[f].mName, f);
    }
    std::map<std::pair<std::size_t, std::int64_t>, std::size_t> use_by_trim;
    for (std::size_t u = 0; u < uses.size(); ++u)
        if (uses[u].mTrimIndex >= 0)
            use_by_trim.emplace(std::make_pair(uses[u].mFace, uses[u].mTrimIndex), u);
    std::vector<McadEdge> edges;
    for (std::size_t e = 0; e < rBrep.mEdges.size(); ++e) {
        McadEdge ed;
        ed.mBrepEdge = static_cast<std::int64_t>(e);
        for (const BrepRef& r_r : rBrep.mEdges[e].mTopology) {
            std::optional<std::size_t> f;
            if (r_r.mId >= 0) {
                const auto it = face_by_id.find(r_r.mId);
                if (it != face_by_id.end())
                    f = it->second;
            } else if (!r_r.mName.empty()) {
                const auto it = face_by_name.find(r_r.mName);
                if (it != face_by_name.end())
                    f = it->second;
            }
            if (!f) {
                warn("edge " + std::to_string(e) + " refers to a face that does not exist");
                continue;
            }
            // A face already refused has no trims; its refusal says why.
            if (faces[*f].mRecord.mStatus != MeshCadFaceStatus::Ok)
                continue;
            const auto it = use_by_trim.find({*f, r_r.mTrimIndex});
            if (it == use_by_trim.end()) {
                warn("edge " + std::to_string(e) + " refers to trim " +
                     std::to_string(r_r.mTrimIndex) + " of face " + std::to_string(*f) +
                     ", which does not exist");
                continue;
            }
            McadUse& r_u = uses[it->second];
            if (!live(it->second))
                continue;
            if (r_u.mEdge >= 0) {
                warn("trim " + std::to_string(r_r.mTrimIndex) + " of face " + std::to_string(*f) +
                     " is used by two edges; the first is kept");
                continue;
            }
            r_u.mEdge = static_cast<std::int64_t>(edges.size());
            ed.mUses.push_back(it->second);
        }
        if (ed.mUses.empty())
            continue;
        std::sort(ed.mUses.begin(), ed.mUses.end());
        ed.mMaster = ed.mUses.front();
        edges.push_back(std::move(ed));
    }
    {
        std::vector<std::size_t> free_uses;
        for (std::size_t u = 0; u < uses.size(); ++u)
            if (live(u) && uses[u].mEdge < 0)
                free_uses.push_back(u);
        struct Probe {
            Vec3 mStart, mMid, mEnd;
            bool mDegenerate;
        };
        std::vector<Probe> probe(free_uses.size());
        for (std::size_t i = 0; i < free_uses.size(); ++i) {
            const detail::TrimTrace& r_t = *traces[free_uses[i]];
            probe[i] = {r_t.Point(r_t.Start()), r_t.Point(r_t.ParamAtLength(0.5 * r_t.Length())),
                        r_t.Point(r_t.End()), !(r_t.Length() > tol)};
        }
        std::vector<bool> taken(free_uses.size(), false);
        for (std::size_t i = 0; i < free_uses.size(); ++i) {
            if (taken[i])
                continue;
            taken[i] = true;
            McadEdge ed;
            ed.mUses.push_back(free_uses[i]);
            if (rOptions.mMatchUnreferencedTrims && !probe[i].mDegenerate) {
                for (std::size_t j = i + 1; j < free_uses.size(); ++j) {
                    if (taken[j] || probe[j].mDegenerate)
                        continue;
                    const Probe& a = probe[i];
                    const Probe& b = probe[j];
                    const bool ends =
                        (mcad_dist(a.mStart, b.mStart) <= tol &&
                         mcad_dist(a.mEnd, b.mEnd) <= tol) ||
                        (mcad_dist(a.mStart, b.mEnd) <= tol && mcad_dist(a.mEnd, b.mStart) <= tol);
                    if (ends && mcad_dist(a.mMid, b.mMid) <= tol) {
                        taken[j] = true;
                        ed.mUses.push_back(free_uses[j]);
                    }
                }
            }
            for (const std::size_t u : ed.mUses)
                uses[u].mEdge = static_cast<std::int64_t>(edges.size());
            ed.mMaster = ed.mUses.front();
            edges.push_back(std::move(ed));
        }
    }

    // A model's own tolerances (6.5.1): an edge's "tolerance" and its ends'
    // "end_tolerances" in its extra JSON (a STEP or IGES read writes
    // OpenCASCADE's), each taken no smaller than the global one, decide for
    // that edge and its corners where its uses start at one point and how
    // far apart they may lie. Degeneracy and sampling keep the global one:
    // a declared tolerance never makes an edge a pole.
    std::vector<double> edge_tol(edges.size(), tol);
    std::vector<double> corner_tol(edges.size(), tol);
    // (edge, where): the largest gap a declared tolerance accepted there.
    std::map<std::pair<std::size_t, std::string>, double> tol_accepted;
    std::vector<double> projection_accepted(edges.size(), 0.0);
    if (extras.mModelTolerances)
        for (std::size_t e = 0; e < edges.size(); ++e) {
            const std::int64_t be = edges[e].mBrepEdge;
            if (be < 0 || rBrep.mEdges[static_cast<std::size_t>(be)].mExtra.empty())
                continue;
            const std::string& r_extra = rBrep.mEdges[static_cast<std::size_t>(be)].mExtra;
            const auto bad = [&]() {
                warn("edge " + std::to_string(be) +
                     ": its declared tolerance is not a finite, non-negative number; the "
                     "global one is used");
            };
            detail::JsonValue doc;
            try {
                doc = detail::parse_json(r_extra, "mesh_cad_surface");
            } catch (const std::exception&) {
                continue;  // not JSON: nothing declared
            }
            if (!doc.IsObject())
                continue;
            const auto number = [](const detail::JsonValue* pV, double& rOut) {
                if (!pV || !pV->IsNumber() || !std::isfinite(pV->mNumber) || pV->mNumber < 0.0)
                    return false;
                rOut = pV->mNumber;
                return true;
            };
            if (const detail::JsonValue* p_t = doc.Find("tolerance")) {
                double t = 0.0;
                if (number(p_t, t))
                    edge_tol[e] = std::max(tol, t);
                else
                    bad();
            }
            corner_tol[e] = edge_tol[e];
            if (const detail::JsonValue* p_ends = doc.Find("end_tolerances")) {
                if (!p_ends->IsArray()) {
                    bad();
                } else {
                    for (const detail::JsonValue& r_v : p_ends->mItems) {
                        double t = 0.0;
                        if (number(&r_v, t))
                            corner_tol[e] = std::max(corner_tol[e], t);
                        else
                            bad();
                    }
                }
            }
        }
    const auto accepted = [&](std::int64_t E, double Gap, const char* pWhere) {
        if (E < 0)
            return;
        double& r_g = tol_accepted[{static_cast<std::size_t>(E), pWhere}];
        r_g = std::max(r_g, Gap);
    };

    // Closed edges whose uses start at different points (a cap's circle
    // listed from another angle than the cylinder's rim): a use alone in its
    // loop and closed in its face's parameters can start anywhere, so it is
    // rotated onto the corner where the edge's anchor starts. The anchor is
    // the first use that cannot be rotated (its loop has other trims), or the
    // first use. Models whose uses start together are not touched.
    const auto closed_3d = [&](std::size_t U) {
        const detail::TrimTrace& r_t = *traces[U];
        return r_t.Length() > tol && mcad_dist(r_t.Point(r_t.Start()), r_t.Point(r_t.End())) <= tol;
    };
    const auto rotatable = [&](std::size_t U) {
        if (uses[U].mLoopSize != 1 || !closed_3d(U))
            return false;
        const detail::TrimTrace& r_t = *traces[U];
        const Vec2 a = r_t.Uv(r_t.Start());
        const Vec2 b = r_t.Uv(r_t.End());
        const auto& r_box = faces[uses[U].mFace].pSurface->Box();
        const double extent = std::hypot(r_box[1] - r_box[0], r_box[3] - r_box[2]);
        return std::hypot(a[0] - b[0], a[1] - b[1]) <= 1e-9 * extent;
    };
    for (std::size_t e = 0; e < edges.size(); ++e) {
        McadEdge& r_e = edges[e];
        if (r_e.mUses.size() < 2 || !closed_3d(r_e.mUses.front()))
            continue;
        const auto start_of = [&](std::size_t U) { return traces[U]->Point(traces[U]->Start()); };
        const double te = edge_tol[e];
        bool apart = false;
        for (const std::size_t u : r_e.mUses)
            apart = apart || mcad_dist(start_of(u), start_of(r_e.mUses.front())) > tol;
        if (!apart)
            continue;
        std::size_t anchor = r_e.mUses.front();
        for (const std::size_t u : r_e.mUses)
            if (!rotatable(u)) {
                anchor = u;
                break;
            }
        r_e.mMaster = anchor;
        const Vec3 p = start_of(anchor);
        for (const std::size_t u : r_e.mUses) {
            const double d = mcad_dist(start_of(u), p);
            if (u == anchor || d <= tol)
                continue;
            if (rotatable(u)) {
                uses[u].mStartLength = detail::edge_closest_length(*traces[u], p);
                continue;
            }
            // A use that cannot be rotated starts where the anchor does
            // within the edge's declared tolerance (6.5.1): accepted, and
            // reported, where the global tolerance fails its face.
            if (d <= te) {
                accepted(static_cast<std::int64_t>(e), d, "closed_start");
                continue;
            }
            // Its whole face fails with this reason, rather than with the
            // corner gaps it causes along its loop.
            for (const auto& r_loop : faces[uses[u].mFace].mLoops)
                for (const std::size_t l : r_loop)
                    if (uses[l].mFail == MeshCadFaceStatus::Ok) {
                        uses[l].mFail = MeshCadFaceStatus::ProjectionFailed;
                        uses[l].mWhy =
                            "the uses of a closed edge start at different points, and "
                            "this one cannot be started elsewhere (its loop has other "
                            "trims)";
                    }
        }
    }
    // Where a use's ends are: a rotated use starts and ends where it was
    // rotated to.
    const auto end_point = [&](std::size_t U, int End) {
        const detail::TrimTrace& r_t = *traces[U];
        if (uses[U].mStartLength >= 0.0)
            return r_t.Point(r_t.ParamAtLength(uses[U].mStartLength));
        return r_t.Point(End == 0 ? r_t.Start() : r_t.End());
    };

    // Corners from the topology: a trim's end is the next trim's start; an
    // edge's ends are the same at every use; a degenerate edge is one corner.
    const auto slot = [](std::size_t U, int End) { return 2 * U + static_cast<std::size_t>(End); };
    McadUnionFind corners(2 * uses.size());
    for (const McadFace& r_f : faces) {
        if (!r_f.mUsable)
            continue;
        for (const auto& r_loop : r_f.mLoops)
            for (std::size_t k = 0; k < r_loop.size(); ++k)
                corners.Join(slot(r_loop[k], 1), slot(r_loop[(k + 1) % r_loop.size()], 0));
    }
    for (McadEdge& r_e : edges) {
        const std::size_t m = r_e.mMaster;
        const detail::TrimTrace& r_m = *traces[m];
        const bool degenerate = !(r_m.Length() > tol);
        if (degenerate)
            corners.Join(slot(m, 0), slot(m, 1));
        for (const std::size_t o : r_e.mUses) {
            if (o == m)
                continue;
            uses[o].mSame =
                degenerate || detail::edge_same_direction(r_m, *traces[o], uses[o].mStartLength);
            corners.Join(slot(m, 0), slot(o, uses[o].mSame ? 0 : 1));
            corners.Join(slot(m, 1), slot(o, uses[o].mSame ? 1 : 0));
            if (degenerate)
                corners.Join(slot(o, 0), slot(o, 1));
        }
    }
    // An embedded edge has no loop to close it: a closed one is one corner,
    // and an end within the tolerance of a corner of its face's loops is that
    // corner (a curve from boundary to boundary).
    for (const McadFace& r_f : faces) {
        if (!r_f.mUsable)
            continue;
        for (const std::size_t u : r_f.mEmbedded) {
            if (!live(u))
                continue;
            if (closed_3d(u))
                corners.Join(slot(u, 0), slot(u, 1));
            for (int end = 0; end < 2; ++end) {
                const Vec3 p = end_point(u, end);
                bool joined = false;
                for (const auto& r_loop : r_f.mLoops)
                    for (const std::size_t l : r_loop)
                        for (int k = 0; k < 2 && !joined; ++k)
                            if (mcad_dist(end_point(l, k), p) <= tol) {
                                corners.Join(slot(u, end), slot(l, k));
                                joined = true;
                            }
            }
        }
    }
    std::vector<std::int64_t> corner_id(2 * uses.size(), -1);
    std::vector<Vec3> boundary;
    for (std::size_t s = 0; s < 2 * uses.size(); ++s) {
        if (!live(s / 2))
            continue;
        const std::size_t root = corners.Find(s);
        if (corner_id[root] < 0) {
            corner_id[root] = static_cast<std::int64_t>(boundary.size());
            boundary.push_back(end_point(s / 2, static_cast<int>(s % 2)));
        }
        corner_id[s] = corner_id[root];
    }
    double max_gap = 0.0;
    for (std::size_t s = 0; s < 2 * uses.size(); ++s) {
        if (!live(s / 2))
            continue;
        const double gap = mcad_dist(end_point(s / 2, static_cast<int>(s % 2)),
                                     boundary[static_cast<std::size_t>(corner_id[s])]);
        max_gap = std::max(max_gap, gap);
        const std::int64_t se = uses[s / 2].mEdge;
        const double fail_at =
            se >= 0 ? std::max(gap_fail, 10.0 * corner_tol[static_cast<std::size_t>(se)])
                    : gap_fail;
        if (gap > gap_fail && gap <= fail_at)
            accepted(se, gap, "corner");
        if (gap > fail_at && uses[s / 2].mFail == MeshCadFaceStatus::Ok) {
            uses[s / 2].mFail = MeshCadFaceStatus::OpenLoop;
            uses[s / 2].mWhy = "a corner of its loop is " + std::to_string(gap) +
                               " from the corner the topology joins it to";
        }
    }

    // Sample each edge once on its master, then carry the samples to the
    // other uses.
    detail::EdgeSampleOptions so;
    so.mTargetSize = h;
    if (pSize)
        so.mSize = [pSize](const detail::EdgeVec3& rX) { return pSize->At(rX); };
    so.mMaxChordalDeviation = rOptions.mMaxChordalDeviation;
    so.mMinSegments = rOptions.mMinPointsPerEdge + 1;
    so.mTolerance = tol;
    so.mMaxSegments = rOptions.mMaxPoints;
    std::vector<double> edge_gap(edges.size(), 0.0);
    const auto carry = [&](std::size_t E) {
        McadEdge& r_e = edges[E];
        const std::size_t m = r_e.mMaster;
        uses[m].mParams = r_e.mSamples.mParams;
        for (const std::size_t o : r_e.mUses) {
            if (o == m)
                continue;
            const detail::EdgeProjection p =
                detail::edge_project(*traces[o], r_e.mSamples, uses[o].mSame, uses[o].mStartLength);
            uses[o].mParams = p.mParams;
            edge_gap[E] = std::max(edge_gap[E], p.mMaxGap);
            const double fail_at = std::max(gap_fail, 10.0 * edge_tol[E]);
            if (p.mMaxGap > gap_fail && p.mMaxGap <= fail_at)
                projection_accepted[E] = std::max(projection_accepted[E], p.mMaxGap);
            if (p.mMaxGap > fail_at && uses[o].mFail == MeshCadFaceStatus::Ok) {
                uses[o].mFail = MeshCadFaceStatus::ProjectionFailed;
                uses[o].mWhy = "the trims of one edge are " + std::to_string(p.mMaxGap) + " apart";
            }
        }
    };
    parallel_for(
        edges.size(),
        [&](std::size_t e) {
            McadEdge& r_e = edges[e];
            try {
                // Enough pieces that every loop through the edge has three
                // points: a closed edge alone in its loop needs three.
                detail::EdgeSampleOptions eo = so;
                // An open embedded curve needs one piece.
                for (const std::size_t u : r_e.mUses) {
                    const std::size_t loop =
                        uses[u].mEmbedded && !closed_3d(u) ? 3 : uses[u].mLoopSize;
                    eo.mMinSegments = std::max<std::int64_t>(
                        eo.mMinSegments, static_cast<std::int64_t>((3 + loop - 1) / loop));
                }
                r_e.mSamples = detail::edge_sample(*traces[r_e.mMaster], eo);
                carry(e);
            } catch (const std::length_error& r_x) {
                r_e.mFail = MeshCadFaceStatus::BudgetExceeded;
                r_e.mWhy = mcad_why(r_x);
            } catch (const std::exception& r_x) {
                r_e.mFail = MeshCadFaceStatus::Other;
                r_e.mWhy = mcad_why(r_x);
            }
        },
        1);

    // Point ids: corners, then each edge's interior samples.
    const std::size_t num_corners = boundary.size();
    const auto assign_ids = [&]() {
        boundary.resize(num_corners);
        for (McadEdge& r_e : edges) {
            r_e.mFirstInterior = static_cast<std::int64_t>(boundary.size());
            if (r_e.mFail != MeshCadFaceStatus::Ok) {
                for (const std::size_t u : r_e.mUses)
                    if (uses[u].mFail == MeshCadFaceStatus::Ok) {
                        uses[u].mFail = r_e.mFail;
                        uses[u].mWhy = r_e.mWhy;
                    }
                continue;
            }
            const auto& r_pts = r_e.mSamples.mPoints;
            if (static_cast<std::int64_t>(boundary.size() + r_pts.size()) > rOptions.mMaxPoints)
                throw std::length_error(std::string(kMcadPrefix) +
                                        "the edges alone need more than " +
                                        std::to_string(rOptions.mMaxPoints) + " points");
            for (std::size_t k = 1; k + 1 < r_pts.size(); ++k)
                boundary.push_back(r_pts[k]);
            const std::size_t n = r_pts.size();
            for (const std::size_t u : r_e.mUses) {
                McadUse& r_u = uses[u];
                r_u.mIds.assign(n, 0);
                r_u.mIds.front() = corner_id[slot(u, 0)];
                r_u.mIds.back() = corner_id[slot(u, 1)];
                for (std::size_t k = 1; k + 1 < n; ++k) {
                    const std::size_t j = r_u.mSame ? k : n - 1 - k;
                    r_u.mIds[k] = r_e.mFirstInterior + static_cast<std::int64_t>(j - 1);
                }
            }
        }
    };
    assign_ids();

    // Virtual edges (6.5.2): a corner where exactly two edges meet, used by
    // the same faces, turning by under 30 degrees, at least one of them
    // shorter than h / 2, is suppressed: no face's loop and no line cell
    // goes through it, so a sliver of an edge forces no tiny triangle. The
    // chord across it stays within a tenth of its length of the corner and
    // within 1.5 h; a corner next to one already suppressed across a single
    // piece is kept, so no chord skips more than one corner. The model is not
    // changed; the chains are reported.
    std::vector<char> suppressed;
    if (VirtualEdges) {
        suppressed.assign(num_corners, 0);
        struct Touch {
            std::size_t mEdge;
            std::size_t mFace;
        };
        std::vector<std::vector<Touch>> touching(num_corners);
        for (std::size_t u = 0; u < uses.size(); ++u) {
            if (!live(u) || uses[u].mEmbedded || uses[u].mEdge < 0)
                continue;
            for (int end = 0; end < 2; ++end)
                touching[static_cast<std::size_t>(corner_id[slot(u, end)])].push_back(
                    {static_cast<std::size_t>(uses[u].mEdge), uses[u].mFace});
        }
        // The sample next to corner C along edge E, and E's other corner.
        const auto beside = [&](std::size_t E, std::int64_t C, Vec3& rNext, std::int64_t& rOther,
                                bool& rOnePiece) {
            const McadEdge& r_e = edges[E];
            const auto& r_pts = r_e.mSamples.mPoints;
            const std::size_t m = r_e.mMaster;
            const bool at_start = corner_id[slot(m, 0)] == C;
            rNext = at_start ? r_pts[1] : r_pts[r_pts.size() - 2];
            rOther = corner_id[slot(m, at_start ? 1 : 0)];
            rOnePiece = r_pts.size() == 2;
        };
        for (std::size_t c = 0; c < num_corners; ++c) {
            const auto& r_t = touching[c];
            std::vector<std::size_t> touched;
            for (const Touch& r_x : r_t)
                touched.push_back(r_x.mEdge);
            std::sort(touched.begin(), touched.end());
            touched.erase(std::unique(touched.begin(), touched.end()), touched.end());
            if (touched.size() != 2)
                continue;
            const std::size_t e1 = touched[0], e2 = touched[1];
            const auto usable = [&](std::size_t E) {
                const McadEdge& r_e = edges[E];
                return r_e.mFail == MeshCadFaceStatus::Ok && !r_e.mSamples.mDegenerate &&
                       r_e.mSamples.mPoints.size() >= 2 &&
                       corner_id[slot(r_e.mMaster, 0)] != corner_id[slot(r_e.mMaster, 1)];
            };
            if (!usable(e1) || !usable(e2))
                continue;
            // The same faces, each using each edge once.
            std::vector<std::size_t> f1, f2;
            for (const std::size_t u : edges[e1].mUses)
                f1.push_back(uses[u].mFace);
            for (const std::size_t u : edges[e2].mUses)
                f2.push_back(uses[u].mFace);
            std::sort(f1.begin(), f1.end());
            std::sort(f2.begin(), f2.end());
            if (f1 != f2 || std::adjacent_find(f1.begin(), f1.end()) != f1.end())
                continue;
            if (!(std::min(edges[e1].mSamples.mLength, edges[e2].mSamples.mLength) < 0.5 * h))
                continue;
            Vec3 p{}, q{};
            std::int64_t o1 = -1, o2 = -1;
            bool one1 = false, one2 = false;
            beside(e1, static_cast<std::int64_t>(c), p, o1, one1);
            beside(e2, static_cast<std::int64_t>(c), q, o2, one2);
            if ((one1 && o1 >= 0 && suppressed[static_cast<std::size_t>(o1)]) ||
                (one2 && o2 >= 0 && suppressed[static_cast<std::size_t>(o2)]))
                continue;
            const Vec3& x = boundary[c];
            const Vec3 a{x[0] - p[0], x[1] - p[1], x[2] - p[2]};
            const Vec3 b{q[0] - x[0], q[1] - x[1], q[2] - x[2]};
            const double la = mcad_dist(p, x), lb = mcad_dist(x, q), chord = mcad_dist(p, q);
            const double cos_turn = (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]) / (la * lb);
            if (!(la > 0.0 && lb > 0.0) || !(cos_turn > std::cos(30.0 * kMcadPi / 180.0)))
                continue;
            const Vec3 mid{0.5 * (p[0] + q[0]), 0.5 * (p[1] + q[1]), 0.5 * (p[2] + q[2])};
            const double local = pSize ? pSize->At(mid) : h;
            // The corner's distance from the chord p-q.
            const Vec3 pq{q[0] - p[0], q[1] - p[1], q[2] - p[2]};
            const Vec3 px{x[0] - p[0], x[1] - p[1], x[2] - p[2]};
            const Vec3 cr{pq[1] * px[2] - pq[2] * px[1], pq[2] * px[0] - pq[0] * px[2],
                          pq[0] * px[1] - pq[1] * px[0]};
            const double off = std::sqrt(cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]) / chord;
            if (!(chord <= 1.5 * local) || !(off <= 0.1 * chord))
                continue;
            suppressed[c] = 1;
        }
        bool any = false;
        for (const char r_c : suppressed)
            any = any || r_c;
        if (any) {
            for (McadFace& r_f : faces)
                r_f.pSuppressed = &suppressed;
            // The chains: edges joined through suppressed corners, walked
            // from an end, in edge-index order of their first edge.
            if (pChains) {
                std::map<std::size_t, std::vector<std::size_t>> link;
                for (std::size_t c = 0; c < num_corners; ++c) {
                    if (!suppressed[c])
                        continue;
                    std::vector<std::size_t> touched;
                    for (const Touch& r_x : touching[c])
                        touched.push_back(r_x.mEdge);
                    std::sort(touched.begin(), touched.end());
                    touched.erase(std::unique(touched.begin(), touched.end()), touched.end());
                    link[touched[0]].push_back(touched[1]);
                    link[touched[1]].push_back(touched[0]);
                }
                std::set<std::size_t> done;
                for (const auto& [e, r_n] : link) {
                    if (done.count(e) || r_n.size() != 1)
                        continue;
                    std::vector<std::int64_t> chain;
                    std::size_t at = e, from = edges.size();
                    while (true) {
                        done.insert(at);
                        chain.push_back(edges[at].mBrepEdge);
                        std::size_t next = edges.size();
                        for (const std::size_t n : link[at])
                            if (n != from && !done.count(n))
                                next = n;
                        if (next == edges.size())
                            break;
                        from = at;
                        at = next;
                    }
                    pChains->push_back(std::move(chain));
                }
            }
        }
    }

    // Boundaries that cross at this size: the sampled loops of a face
    // can cross in its parameters where the exact curves do not (two close
    // loops, a thin slot). Every crossing found exactly is separated by
    // halving, by arc length, the master pieces of the edges involved, in
    // every face that uses them, until no face has one. The requests of all
    // faces are merged in (edge, piece) order, so the result does not depend
    // on the thread count. A crossing whose pieces are already within the
    // tolerance is where the curves themselves meet: that face stays
    // self_intersecting, as does one still crossing after the last round.
    const auto fail_face = [&](std::size_t F, MeshCadFaceStatus S, const std::string& rWhy) {
        for (const auto& r_loop : faces[F].mLoops)
            for (const std::size_t u : r_loop)
                if (uses[u].mFail == MeshCadFaceStatus::Ok) {
                    uses[u].mFail = S;
                    uses[u].mWhy = rWhy;
                }
    };
    const auto separate = [&](std::vector<char> check) {
        constexpr int kMaxRounds = 64;
        for (int round = 0;; ++round) {
            std::vector<McadCrossings> found(nf);
            const std::vector<std::size_t> order = mcad_face_order(nf);
            parallel_for(
                nf,
                [&](std::size_t i) {
                    const std::size_t f = order[i];
                    if (check[f] && faces[f].mUsable)
                        found[f] = mcad_face_crossings(faces[f], uses, traces, edges, tol);
                },
                1);
            std::map<std::size_t, std::vector<std::size_t>> split;
            std::vector<std::size_t> asking;
            std::fill(check.begin(), check.end(), 0);
            for (std::size_t f = 0; f < nf; ++f) {
                const McadCrossings& r_c = found[f];
                if (r_c.mStuck) {
                    fail_face(f, MeshCadFaceStatus::SelfIntersecting,
                              "its boundary meets itself in the face's parameters where the "
                              "curves themselves meet (no split separates them)");
                    continue;
                }
                if (!r_c.mDrop.empty()) {
                    // Left out, then the face is searched again without it.
                    for (const std::size_t u : r_c.mDrop)
                        if (uses[u].mFail == MeshCadFaceStatus::Ok) {
                            uses[u].mFail = MeshCadFaceStatus::SelfIntersecting;
                            uses[u].mWhy =
                                "it crosses or touches the face's boundary or "
                                "another embedded edge at this size";
                        }
                    check[f] = 1;
                    continue;
                }
                if (r_c.mAsks.empty())
                    continue;
                asking.push_back(f);
                for (const auto& r_a : r_c.mAsks)
                    split[r_a[0]].push_back(r_a[1]);
            }
            if (split.empty()) {
                if (std::find(check.begin(), check.end(), 1) == check.end())
                    break;
                continue;
            }
            std::size_t added = 0;
            for (auto& [e, r_pieces] : split) {
                std::sort(r_pieces.begin(), r_pieces.end());
                r_pieces.erase(std::unique(r_pieces.begin(), r_pieces.end()), r_pieces.end());
                added += r_pieces.size();
            }
            if (round >= kMaxRounds ||
                static_cast<std::int64_t>(boundary.size() + added) > rOptions.mMaxPoints) {
                for (const std::size_t f : asking)
                    if (round >= kMaxRounds)
                        fail_face(f, MeshCadFaceStatus::SelfIntersecting,
                                  "its boundary still crosses itself in the face's parameters "
                                  "after " +
                                      std::to_string(kMaxRounds) + " rounds of splits");
                    else
                        fail_face(f, MeshCadFaceStatus::BudgetExceeded,
                                  "separating its crossing boundary would exceed max_points");
                break;
            }
            for (const auto& [e, r_pieces] : split) {
                McadEdge& r_e = edges[e];
                detail::edge_refine(*traces[r_e.mMaster], r_e.mSamples, r_pieces);
                carry(e);
                for (const std::size_t u : r_e.mUses)
                    check[uses[u].mFace] = 1;
            }
            assign_ids();
        }
    };
    separate(std::vector<char>(nf, 1));
    for (const double g : edge_gap)
        max_gap = std::max(max_gap, g);
    result.mMaxEdgeGap = max_gap;
    if (max_gap > 0.01 * h)
        warn("the trims of some edges or corners disagree by up to " + std::to_string(max_gap) +
             " (the mesh shares their points regardless)");

    // Faces, each in its own parameters.
    McadQuality quality;
    quality.mMinAngleDeg = rOptions.mMinAngleDeg;
    quality.mSmoothIterations = rOptions.mSmoothIterations;
    quality.mSmoothMethod = rOptions.mSmoothMethod;
    quality.pSize = pSize;
    quality.mCollectSizes = pSamples != nullptr;
    // Corner lopping (6.2.3): each corner's ball, the smallest radius any
    // face's meshing or any split left, by the corner's position (its id
    // stays, but the position names it whatever is renumbered).
    std::map<Vec3, double> lop_radius;
    quality.mLopCorners = extras.mLopCorners;
    quality.pLopRadius = &lop_radius;
    // A pole is the corner of a degenerate edge; corner ids come first, so
    // they do not move when edges are split.
    quality.mNumFaces = static_cast<std::int64_t>(
        std::count_if(faces.begin(), faces.end(), [](const McadFace& rF) { return rF.mUsable; }));
    quality.mPole.assign(num_corners, 0);
    for (const McadEdge& r_e : edges)
        if (r_e.mFail == MeshCadFaceStatus::Ok && r_e.mSamples.mDegenerate)
            for (const std::size_t u : r_e.mUses)
                if (!uses[u].mIds.empty() && uses[u].mIds.front() >= 0 &&
                    static_cast<std::size_t>(uses[u].mIds.front()) < num_corners)
                    quality.mPole[static_cast<std::size_t>(uses[u].mIds.front())] = 1;
    // A face's state before each meshing.
    const auto reset_face = [&](std::size_t F) {
        McadFace& r_g = faces[F];
        r_g.mSmoothingDeferred = false;
        r_g.mRecord = MeshCadFaceRecord{};
        r_g.mRecord.mFace = static_cast<std::int64_t>(F);
        r_g.mRecord.mId = rBrep.mFaces[F].mId;
        r_g.mInterior.clear();
        r_g.mTriangles.clear();
        r_g.mAsks.clear();
        r_g.mLopAsks.clear();
        r_g.mLopped.clear();
        r_g.mCurvedEdges.clear();
        r_g.mCurvedCentre.clear();
        r_g.mCurvedNormal.clear();
        r_g.mHighOrder.clear();
        r_g.mTriOwner.clear();
        r_g.mSizeSamples.clear();
        r_g.mHistogram.fill(0);
        r_g.mNumStretched = 0;
        r_g.mRefineBudgetReached = false;
        // What meshing adds to the skipped embedded edges is added
        // again by the next meshing.
        if (r_g.mSkippedBefore == std::string::npos)
            r_g.mSkippedBefore = r_g.mSkipped.size();
        r_g.mSkipped.resize(r_g.mSkippedBefore);
    };
    // Every usable face, or with `Again` those whose last meshing left its
    // smoothing for the next one, meshed once more from where that one
    // started (the same boundary and refinement budget give the same mesh).
    const auto mesh_faces = [&](bool Defer, bool Again, const std::vector<char>* pOnly = nullptr) {
        const std::vector<std::size_t> order = mcad_face_order(nf);
        McadQuality q = quality;
        q.mDeferSmoothing = Defer && (q.mMinAngleDeg > 0.0 || q.pSize) && q.mSmoothIterations > 0;
        parallel_for(
            nf,
            [&](std::size_t i) {
                const std::size_t f = order[i];
                McadFace& r_f = faces[f];
                // A merged sliver is meshed by its host's task (6.5.2).
                if (!r_f.mUsable || r_f.mMergedInto >= 0 || (Again && !r_f.mSmoothingDeferred) ||
                    (pOnly && !(*pOnly)[f]))
                    return;
                if (Again)
                    r_f.mRefinedSoFar = r_f.mRefinedBefore;
                r_f.mRefinedBefore = r_f.mRefinedSoFar;
                const auto reset = [&]() { reset_face(f); };
                // A flattening's coarse mesh is the uniform lattice, even in a
                // sized run: it only carries the surface's map, while a lattice
                // at the local size would let the conformal solve, whose
                // iterations grow with the mesh, outweigh the meshing it serves
                // (found by the fuzz target on a roof sized to a quarter of h).
                const auto coarse_quality = [](const McadQuality& rQq) {
                    McadQuality c = rQq;
                    c.pSize = nullptr;
                    return c;
                };
                // The discrete reparametrisation (6.4.4): the face's unrefined
                // lattice is its coarse mesh, measured on its first meshing;
                // a face over the threshold is meshed in the flattening of
                // that mesh, rebuilt each meshing since splits move its
                // boundary, and in its parameters when it cannot be.
                // Virtual faces (6.5.2): a host and its sliver meshed as one,
                // in the flattening of their glued coarse meshes, the merged
                // loop leaving out the edges they share. False (the merge
                // undone, both faces meshed on their own) when they cannot be.
                const auto mesh_merged = [&](const McadQuality& rQq) {
                    const auto s = static_cast<std::size_t>(r_f.mMergedWith);
                    McadFace& r_sl = faces[s];
                    std::string why;
                    detail::McadCoarseMesh ch, cs;
                    mcad_mesh_face(r_f, rBrep.mFaces[f], uses, traces, boundary, rOptions,
                                   coarse_quality(rQq), &ch);
                    reset();
                    mcad_mesh_face(r_sl, rBrep.mFaces[s], uses, traces, boundary, rOptions,
                                   coarse_quality(rQq), &cs);
                    reset_face(s);
                    if (!ch.mFilled || !cs.mFilled)
                        why = "one of the two has no unrefined triangulation to flatten";
                    else if (!ch.mWhy.empty())
                        why = ch.mWhy;
                    detail::McadCoarseMesh glued;
                    if (why.empty()) {
                        std::set<std::int64_t> shared;
                        for (const std::size_t e : r_f.mSharedEdges)
                            for (const std::size_t u : edges[e].mUses)
                                shared.insert(uses[u].mIds.begin(), uses[u].mIds.end());
                        glued = mcad_glue(ch, cs, shared, why);
                    }
                    detail::McadDisc disc;
                    if (why.empty())
                        why = detail::mcad_build_disc(
                            glued, {r_f.pSurface.get(), r_sl.pSurface.get()}, disc);
                    if (why.empty()) {
                        std::swap(r_f.mLoops, r_f.mMergedLoops);
                        std::swap(r_f.mLoopReversed, r_f.mMergedReversed);
                        mcad_mesh_face(r_f, rBrep.mFaces[f], uses, traces, boundary, rOptions, rQq,
                                       nullptr, &disc);
                        std::swap(r_f.mLoops, r_f.mMergedLoops);
                        std::swap(r_f.mLoopReversed, r_f.mMergedReversed);
                        const MeshCadFaceStatus st = r_f.mRecord.mStatus;
                        if (st == MeshCadFaceStatus::Ok || r_f.mFolded ||
                            st == MeshCadFaceStatus::BudgetExceeded)
                            return true;
                        why = "meshing their flattening failed: " + r_f.mRecord.mReason;
                        reset();
                    }
                    r_f.mMergedWith = -1;
                    r_sl.mMergedInto = -1;
                    r_f.mMergeWhy = why;
                    mcad_mesh_face(r_sl, rBrep.mFaces[s], uses, traces, boundary, rOptions, rQq);
                    return false;
                };
                const auto mesh_one = [&](const McadQuality& rQq) {
                    if (r_f.mMergedWith >= 0 && mesh_merged(rQq))
                        return;
                    const double threshold = rOptions.mReparamConditionNumber;
                    using Route = McadFace::Route;
                    if (!(threshold > 0.0) || r_f.mRoute == Route::Below ||
                        r_f.mRoute == Route::Declined) {
                        mcad_mesh_face(r_f, rBrep.mFaces[f], uses, traces, boundary, rOptions, rQq);
                        return;
                    }
                    detail::McadCoarseMesh coarse;
                    mcad_mesh_face(r_f, rBrep.mFaces[f], uses, traces, boundary, rOptions,
                                   coarse_quality(rQq), &coarse);
                    reset();
                    const auto decline = [&](const std::string& rWhy) {
                        r_f.mRoute = Route::Declined;
                        r_f.mRouteWhy = rWhy;
                        mcad_mesh_face(r_f, rBrep.mFaces[f], uses, traces, boundary, rOptions, rQq);
                    };
                    if (!coarse.mFilled) {
                        // It fails, or has no lattice to flatten, either way.
                        if (r_f.mRoute == Route::Unmeasured)
                            r_f.mRoute = Route::Below;
                        mcad_mesh_face(r_f, rBrep.mFaces[f], uses, traces, boundary, rOptions, rQq);
                        return;
                    }
                    const std::vector<const detail::BrepSurfaceEval*> surfaces{r_f.pSurface.get()};
                    if (r_f.mRoute == Route::Unmeasured) {
                        r_f.mRouteCond = detail::mcad_coarse_condition(
                            coarse, surfaces, coarse.mL11, coarse.mL12, coarse.mL22);
                        if (!(r_f.mRouteCond > threshold)) {
                            r_f.mRoute = Route::Below;
                            mcad_mesh_face(r_f, rBrep.mFaces[f], uses, traces, boundary, rOptions,
                                           rQq);
                            return;
                        }
                        r_f.mRoute = Route::Taken;
                    }
                    if (!coarse.mWhy.empty())
                        return decline(coarse.mWhy);
                    detail::McadDisc disc;
                    const std::string why = detail::mcad_build_disc(coarse, surfaces, disc);
                    if (!why.empty())
                        return decline(why);
                    const double cond_disc = detail::mcad_disc_condition(disc);
                    if (r_f.mRouteMethod == 0) {
                        // The first flattening is judged against the face's
                        // own parameters; later ones only follow its
                        // boundary splits.
                        r_f.mRouteCondDisc = cond_disc;
                        if (!(cond_disc < r_f.mRouteCond))
                            return decline("its flattening is no better than its parameters");
                    }
                    r_f.mRouteMethod = disc.mMethod;
                    mcad_mesh_face(r_f, rBrep.mFaces[f], uses, traces, boundary, rOptions, rQq,
                                   nullptr, &disc);
                    const MeshCadFaceStatus st = r_f.mRecord.mStatus;
                    const bool structural = st != MeshCadFaceStatus::Ok && !r_f.mFolded &&
                                            st != MeshCadFaceStatus::BudgetExceeded;
                    if (structural) {
                        const std::string what = r_f.mRecord.mReason;
                        reset();
                        decline("meshing its flattening failed: " + what);
                    }
                };
                const auto once = [&]() {
                    reset();
                    try {
                        if (r_f.mUnrefined) {
                            McadQuality unrefined = quality;
                            unrefined.mMinAngleDeg = 0.0;
                            unrefined.mNoSizeRefine = true;
                            mesh_one(unrefined);
                        } else {
                            mesh_one(q);
                        }
                    } catch (const std::exception& r_x) {
                        r_f.mRecord.mStatus = MeshCadFaceStatus::Other;
                        r_f.mRecord.mReason = mcad_why(r_x);
                        r_f.mTriangles.clear();
                        r_f.mInterior.clear();
                        r_f.mAsks.clear();
                        r_f.mLopAsks.clear();
                        r_f.mLopped.clear();
                        r_f.mCurvedEdges.clear();
                        r_f.mCurvedCentre.clear();
                        r_f.mCurvedNormal.clear();
                        r_f.mHighOrder.clear();
                    }
                };
                once();
                r_f.mFoldStreak = r_f.mFolded ? r_f.mFoldStreak + 1 : 0;
                // Every triangle collapsed: once more with an interior grid
                // (6.5.7), kept for every later meshing of the face.
                if (r_f.mCollapsed && !r_f.mForceInterior) {
                    r_f.mForceInterior = true;
                    r_f.mRefinedSoFar = r_f.mRefinedBefore;
                    once();
                }
            },
            1);
        // The balls the faces left, after them all (the smallest, so the
        // order does not matter), for the next meshing.
        for (const McadFace& r_f : faces)
            for (const auto& [x, r_c] : r_f.mLopped) {
                auto [it, fresh] = lop_radius.emplace(x, r_c.mLopRadius);
                if (!fresh)
                    it->second = std::min(it->second, r_c.mLopRadius);
            }
    };
    mesh_faces(true, false);

    // Virtual faces (6.5.2): a face narrower than h (twice its area over its
    // perimeter, on its first meshing), which then has no interior point and
    // forces triangles as thin as itself, is meshed with the neighbour it
    // shares the most edge length with, along manifold edges only, when the
    // two meet tangentially there (their normals within 30 degrees at every
    // shared sample, the virtual edges' turn) and the shared edges are one
    // run of each loop. Their glued coarse meshes are flattened as one
    // (6.4.4), so the host's triangulation spans both; each triangle keeps
    // the face its centroid lies in. The model is not changed (Sheffer et
    // al. 2000 build the same virtual topology).
    std::vector<MeshCadMergedFace> merged_found;
    if (extras.mVirtualFaces) {
        const auto point_of = [&](const McadFace& rF, std::int64_t Id) -> const Vec3& {
            return Id >= 0 ? boundary[static_cast<std::size_t>(Id)]
                           : rF.mInterior[static_cast<std::size_t>(-Id - 1)];
        };
        std::vector<double> width(nf, std::numeric_limits<double>::infinity());
        std::vector<char> sliver(nf, 0);
        for (std::size_t f = 0; f < nf; ++f) {
            const McadFace& r_f = faces[f];
            if (!r_f.mUsable || r_f.mRecord.mStatus != MeshCadFaceStatus::Ok ||
                r_f.mTriangles.empty())
                continue;
            double area = 0.0;
            Vec3 centre{0.0, 0.0, 0.0};
            for (const auto& r_t : r_f.mTriangles) {
                const Vec3& a = point_of(r_f, r_t[0]);
                const Vec3& b = point_of(r_f, r_t[1]);
                const Vec3& c = point_of(r_f, r_t[2]);
                const Vec3 e1{b[0] - a[0], b[1] - a[1], b[2] - a[2]};
                const Vec3 e2{c[0] - a[0], c[1] - a[1], c[2] - a[2]};
                const Vec3 n{e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2],
                             e1[0] * e2[1] - e1[1] * e2[0]};
                area += 0.5 * std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
                for (std::size_t k = 0; k < 3; ++k)
                    centre[k] += (a[k] + b[k] + c[k]) / (3.0 * double(r_f.mTriangles.size()));
            }
            double perimeter = 0.0;
            for (const auto& r_loop : r_f.mLoops)
                for (const std::size_t u : r_loop)
                    for (std::size_t k = 0; k + 1 < uses[u].mIds.size(); ++k)
                        perimeter +=
                            mcad_dist(boundary[static_cast<std::size_t>(uses[u].mIds[k])],
                                      boundary[static_cast<std::size_t>(uses[u].mIds[k + 1])]);
            if (!(perimeter > 0.0))
                continue;
            width[f] = 2.0 * area / perimeter;
            const double hf = pSize ? pSize->At(centre) : h;
            sliver[f] = width[f] < hf ? 1 : 0;
        }
        // The outward normal of face F at the parameters of use U's sample K.
        const auto normal_at = [&](std::size_t F, std::size_t U, std::size_t K) {
            const Vec2 uv = traces[U]->Uv(uses[U].mParams[K]);
            Vec3 x{}, su{}, sv{};
            faces[F].pSurface->Eval(uv[0], uv[1], x.data(), su.data(), sv.data());
            Vec3 n{su[1] * sv[2] - su[2] * sv[1], su[2] * sv[0] - su[0] * sv[2],
                   su[0] * sv[1] - su[1] * sv[0]};
            if (rBrep.mFaces[F].mSwappedNormal)
                n = {-n[0], -n[1], -n[2]};
            return n;
        };
        std::vector<char> hosting(nf, 0);
        for (std::size_t s = 0; s < nf; ++s) {
            if (!sliver[s])
                continue;
            MeshCadMergedFace rec;
            rec.mSliver = static_cast<std::int64_t>(s);
            rec.mWidth = width[s];
            const McadFace& r_s = faces[s];
            const auto decline = [&](const std::string& rWhy) {
                rec.mReason = rWhy;
                merged_found.push_back(rec);
            };
            if (r_s.mLoops.size() != 1 || !r_s.mEmbedded.empty()) {
                decline("it has holes or embedded edges");
                continue;
            }
            // The neighbour across the most edge length.
            std::map<std::size_t, double> along;
            for (const std::size_t u : r_s.mLoops[0]) {
                const std::int64_t e = uses[u].mEdge;
                if (e < 0)
                    continue;
                const McadEdge& r_e = edges[static_cast<std::size_t>(e)];
                if (r_e.mUses.size() != 2 || r_e.mFail != MeshCadFaceStatus::Ok ||
                    r_e.mSamples.mDegenerate)
                    continue;
                const std::size_t o = r_e.mUses[0] == u ? r_e.mUses[1] : r_e.mUses[0];
                if (uses[o].mFace != s && !uses[o].mEmbedded)
                    along[uses[o].mFace] += r_e.mSamples.mLength;
            }
            // The neighbours by the edge length they share, most first (ties
            // to the lower index), each tried until one meets the sliver
            // tangentially along one run of both loops.
            std::vector<std::pair<double, std::size_t>> order;
            for (const auto& [g, len] : along)
                if (!sliver[g] && !hosting[g] && faces[g].mUsable &&
                    faces[g].mRecord.mStatus == MeshCadFaceStatus::Ok)
                    order.emplace_back(-len, g);
            std::sort(order.begin(), order.end());
            if (order.empty()) {
                decline("no neighbour that is not itself a sliver or a host");
                continue;
            }
            // The shared run in each loop: the uses whose edge's other use is
            // the other face's.
            const auto shares = [&](std::size_t U, std::size_t Other) {
                const std::int64_t e = uses[U].mEdge;
                if (e < 0 || uses[U].mEmbedded)
                    return false;
                const McadEdge& r_e = edges[static_cast<std::size_t>(e)];
                if (r_e.mUses.size() != 2 || r_e.mFail != MeshCadFaceStatus::Ok ||
                    r_e.mSamples.mDegenerate)
                    return false;
                const std::size_t o = r_e.mUses[0] == U ? r_e.mUses[1] : r_e.mUses[0];
                return uses[o].mFace == Other;
            };
            // The first position of a cyclic run of the uses shared with
            // `Other` in `rLoop`, and its length; false when they are not one
            // run (or are the whole loop).
            const auto run_of = [&](const std::vector<std::size_t>& rLoop, std::size_t Other,
                                    std::size_t& rFirst, std::size_t& rLen) {
                const std::size_t n = rLoop.size();
                std::size_t count = 0, starts = 0;
                for (std::size_t i = 0; i < n; ++i) {
                    const bool in = shares(rLoop[i], Other);
                    count += in ? 1 : 0;
                    if (in && !shares(rLoop[(i + n - 1) % n], Other)) {
                        ++starts;
                        rFirst = i;
                    }
                }
                rLen = count;
                return starts == 1 && count > 0 && count < n;
            };
            // Merge with `G` when it can be; otherwise why not.
            const auto try_host = [&](std::size_t G) -> std::string {
                std::size_t ps = 0, ns = 0;
                if (!run_of(r_s.mLoops[0], G, ps, ns))
                    return "the edges it shares with its neighbour are not one run of its loop";
                std::size_t lh = faces[G].mLoops.size(), ph = 0, nh = 0;
                for (std::size_t l = 0; l < faces[G].mLoops.size(); ++l) {
                    std::size_t first = 0, len = 0;
                    if (run_of(faces[G].mLoops[l], s, first, len)) {
                        lh = l;
                        ph = first;
                        nh = len;
                    }
                }
                if (lh == faces[G].mLoops.size() || nh != ns)
                    return "the edges it shares are not one run of its neighbour's loop";
                // Tangency at every shared sample.
                std::vector<std::size_t> shared;
                for (std::size_t i = 0; i < ns; ++i) {
                    const std::size_t u = r_s.mLoops[0][(ps + i) % r_s.mLoops[0].size()];
                    const McadEdge& r_e = edges[static_cast<std::size_t>(uses[u].mEdge)];
                    shared.push_back(static_cast<std::size_t>(uses[u].mEdge));
                    const std::size_t o = r_e.mUses[0] == u ? r_e.mUses[1] : r_e.mUses[0];
                    const std::size_t n = uses[u].mIds.size();
                    for (std::size_t k = 0; k < n; ++k) {
                        const std::size_t ko = uses[o].mIds[k] == uses[u].mIds[k] ? k : n - 1 - k;
                        const Vec3 a = normal_at(s, u, k);
                        const Vec3 b = normal_at(G, o, ko);
                        const double la = std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
                        const double lb = std::sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]);
                        if (!(la > 0.0) || !(lb > 0.0))
                            continue;
                        if (!((a[0] * b[0] + a[1] * b[1] + a[2] * b[2]) / (la * lb) >
                              std::cos(30.0 * kMcadPi / 180.0)))
                            return "it does not meet its neighbour tangentially (the normals "
                                   "turn 30 degrees or more)";
                    }
                }
                // The merged loop: the host's from after the run round to
                // its start (corner B to corner A), then the sliver's other
                // uses from A to B, forwards or backwards.
                const auto& r_hl = faces[G].mLoops[lh];
                const auto& r_sl = r_s.mLoops[0];
                const std::int64_t corner_a = uses[r_hl[ph]].mIds.front();
                std::vector<std::size_t> loop;
                std::vector<char> back;
                for (std::size_t i = nh; i < r_hl.size(); ++i) {
                    loop.push_back(r_hl[(ph + i) % r_hl.size()]);
                    back.push_back(0);
                }
                const std::size_t rest = r_sl.size() - ns;
                const std::size_t after = (ps + ns) % r_sl.size();
                const std::size_t before = (ps + r_sl.size() - 1) % r_sl.size();
                if (uses[r_sl[after]].mIds.front() == corner_a) {
                    for (std::size_t i = 0; i < rest; ++i) {
                        loop.push_back(r_sl[(after + i) % r_sl.size()]);
                        back.push_back(0);
                    }
                } else if (uses[r_sl[before]].mIds.back() == corner_a) {
                    for (std::size_t i = 0; i < rest; ++i) {
                        loop.push_back(r_sl[(before + r_sl.size() - i) % r_sl.size()]);
                        back.push_back(1);
                    }
                } else {
                    return "the edges it shares do not end at its neighbour's corners";
                }
                McadFace& r_h = faces[G];
                r_h.mMergedLoops = r_h.mLoops;
                r_h.mMergedLoops[lh] = std::move(loop);
                r_h.mMergedReversed.assign(r_h.mLoops.size(), {});
                for (std::size_t l = 0; l < r_h.mLoops.size(); ++l)
                    r_h.mMergedReversed[l].assign(r_h.mMergedLoops[l].size(), 0);
                r_h.mMergedReversed[lh] = std::move(back);
                r_h.mMergedWith = static_cast<std::int64_t>(s);
                r_h.mSharedEdges = shared;
                faces[s].mMergedInto = static_cast<std::int64_t>(G);
                hosting[G] = 1;
                rec.mHost = static_cast<std::int64_t>(G);
                for (const std::size_t e : shared)
                    rec.mEdges.push_back(edges[e].mBrepEdge);
                return {};
            };
            std::string why;
            for (const auto& r_o : order) {
                why = try_host(r_o.second);
                if (why.empty())
                    break;
            }
            if (!why.empty()) {
                decline(why);
                continue;
            }
            rec.mMerged = true;
            merged_found.push_back(rec);
        }
        std::vector<char> again(nf, 0);
        bool any = false;
        for (std::size_t f = 0; f < nf; ++f)
            if (faces[f].mMergedWith >= 0) {
                again[f] = 1;
                any = true;
            }
        if (any)
            mesh_faces(true, false, &again);
    }

    // Boundary splits refinement asks for (6.2.2): a boundary piece is an
    // edge's, shared with its neighbours, so it is halved by arc length in
    // every use (at a power-of-two distance from a sharp corner: concentric
    // shells), the crossings are searched again, and every face is meshed
    // again, until none asks. The requests of all faces are merged in
    // (edge, piece) order, so the result does not depend on the thread
    // count. A piece within the tolerance, or a degenerate edge, is not
    // split: its face keeps the triangles it has and reports them.
    std::int64_t boundary_refinement = 0;
    // A face that folds asks for splits as refinement does (6.4.3).
    const bool any_fold =
        std::any_of(faces.begin(), faces.end(), [](const McadFace& rF) { return rF.mNeedsRounds; });
    if (rOptions.mMinAngleDeg > 0.0 || pSize || any_fold) {
        constexpr int kMaxRefineRounds = 16;
        // Folded in the varying metric (6.4.1) for `Streak` meshings in a
        // row, or when no split is left to make: the metric at an edge's
        // midpoint is the surface's to first order only, and on a small,
        // strongly curved face it can keep a diagonal that lifts backwards
        // whatever the boundary. Such a face is meshed again in its mean
        // metric, kept for every later meshing of it, and the rounds go on.
        const auto to_mean = [&](int Streak) {
            if (rOptions.mMetric != MeshCadMetric::Varying)
                return false;
            std::vector<char> again(nf, 0);
            bool any = false;
            for (std::size_t f = 0; f < nf; ++f)
                if (faces[f].mFolded && !faces[f].mMeanMetric && faces[f].mFoldStreak >= Streak) {
                    faces[f].mMeanMetric = true;
                    again[f] = 1;
                    any = true;
                }
            if (any)
                mesh_faces(true, false, &again);
            return any;
        };
        // A host and its sliver (6.5.2) folded in their flattening for
        // `Streak` meshings in a row, or when no split is left: the merge is
        // undone, both faces are meshed on their own from then on, and the
        // rounds go on, so a pair is not failed that apart would mesh.
        const auto unmerge = [&](int Streak) {
            std::vector<char> again(nf, 0);
            bool any = false;
            for (std::size_t f = 0; f < nf; ++f) {
                McadFace& r_h = faces[f];
                if (r_h.mMergedWith < 0 || !r_h.mFolded || r_h.mFoldStreak < Streak)
                    continue;
                const auto s = static_cast<std::size_t>(r_h.mMergedWith);
                faces[s].mMergedInto = -1;
                r_h.mMergedWith = -1;
                r_h.mMergeWhy = "their flattening kept a triangle folded";
                again[f] = again[s] = 1;
                any = true;
            }
            if (any)
                mesh_faces(true, false, &again);
            return any;
        };
        for (int round = 0;; ++round) {
            // A face whose refinement spent its points is not refined
            // further, since every round would spend them again: it is
            // meshed unrefined on its current boundary, which stays shared,
            // and warned; the rounds go on for the others, whose asks are
            // still to be split (a fold's among them).
            bool spent = false;
            for (std::size_t f = 0; f < nf; ++f)
                if (faces[f].mRefineBudgetReached && !faces[f].mUnrefined) {
                    faces[f].mUnrefined = true;
                    spent = true;
                    warn("face " + std::to_string(f) +
                         ": refinement did not finish within its points (its share of "
                         "max_points, at most twice its unrefined ones); the face is "
                         "meshed unrefined");
                }
            if (spent) {
                mesh_faces(true, false);
                continue;
            }
            if (to_mean(3) || unmerge(3))
                continue;
            std::map<std::size_t, std::map<std::size_t, std::size_t>> split;
            for (std::size_t f = 0; f < nf; ++f)
                for (const auto& r_a : faces[f].mAsks) {
                    const McadEdge& r_e = edges[r_a[0]];
                    if (r_e.mFail != MeshCadFaceStatus::Ok || r_e.mSamples.mDegenerate ||
                        !(mcad_dist(r_e.mSamples.mPoints[r_a[1]],
                                    r_e.mSamples.mPoints[r_a[1] + 1]) > tol))
                        continue;
                    split[r_a[0]][r_a[1]] |= r_a[2];
                }
            // Corner lopping (6.2.3): a piece is split where it leaves a
            // corner's ball (the smallest radius asked at the start's end,
            // else at the end's), and a piece turning more than pi / 3, past
            // what Pav and Walkington assume of a curve, at its midpoint.
            std::map<std::size_t, std::map<std::size_t, McadLopAsk>> lop_split;
            if (extras.mLopCorners) {
                for (std::size_t f = 0; f < nf; ++f)
                    for (const McadLopAsk& r_a : faces[f].mLopAsks) {
                        const McadEdge& r_e = edges[r_a.mEdge];
                        if (r_e.mFail != MeshCadFaceStatus::Ok || r_e.mSamples.mDegenerate ||
                            r_a.mPiece + 1 >= r_e.mSamples.mPoints.size() ||
                            (r_a.mEnd != 1 && r_a.mEnd != 2))
                            continue;
                        split[r_a.mEdge][r_a.mPiece] |= 0;
                        auto [it, fresh] = lop_split[r_a.mEdge].emplace(r_a.mPiece, r_a);
                        if (!fresh &&
                            (r_a.mEnd < it->second.mEnd ||
                             (r_a.mEnd == it->second.mEnd && r_a.mDistance < it->second.mDistance)))
                            it->second = r_a;
                    }
                for (std::size_t e = 0; e < edges.size(); ++e) {
                    const McadEdge& r_e = edges[e];
                    if (r_e.mFail != MeshCadFaceStatus::Ok || r_e.mSamples.mDegenerate ||
                        !traces[r_e.mMaster])
                        continue;
                    const auto& r_p = r_e.mSamples.mParams;
                    for (std::size_t k = 0; k + 1 < r_p.size(); ++k)
                        if (mcad_piece_turning(*traces[r_e.mMaster], r_p[k], r_p[k + 1]) >
                                kMcadPi / 3.0 &&
                            mcad_dist(r_e.mSamples.mPoints[k], r_e.mSamples.mPoints[k + 1]) > tol)
                            split[e][k] |= 0;
                }
            }
            if (split.empty()) {
                if (!to_mean(1) && !unmerge(1))
                    break;
                continue;
            }
            std::size_t added = 0;
            for (const auto& [e, r_pieces] : split)
                added += r_pieces.size();
            if (round >= kMaxRefineRounds) {
                warn("refinement stopped after " + std::to_string(kMaxRefineRounds) +
                     " rounds of boundary splits");
                unmerge(1);
                break;
            }
            if (static_cast<std::int64_t>(boundary.size() + added) > rOptions.mMaxPoints) {
                warn("refinement stopped: its boundary splits would exceed max_points");
                unmerge(1);
                break;
            }
            for (const auto& [e, r_pieces] : split) {
                McadEdge& r_e = edges[e];
                std::vector<std::size_t> pieces;
                std::vector<double> at;
                const auto lop_it = lop_split.find(e);
                for (const auto& [piece, end] : r_pieces) {
                    double a = 0.5;
                    const McadLopAsk* p_lop = nullptr;
                    if (lop_it != lop_split.end()) {
                        const auto it = lop_it->second.find(piece);
                        if (it != lop_it->second.end())
                            p_lop = &it->second;
                    }
                    if (p_lop) {
                        // Where the 3-D distance from the corner is the
                        // ball's radius: bisected along the trace, then the
                        // corner's radius is that split's own distance.
                        const detail::TrimTrace& r_tr = *traces[r_e.mMaster];
                        const auto& r_fr = r_e.mSamples.mFractions;
                        const std::int64_t corner_id =
                            uses[r_e.mMaster].mIds[p_lop->mEnd == 1 ? piece : piece + 1];
                        const Vec3 corner = boundary[static_cast<std::size_t>(corner_id)];
                        const auto at_fraction = [&](double A) {
                            const double f = (1.0 - A) * r_fr[piece] + A * r_fr[piece + 1];
                            return r_tr.Point(r_tr.ParamAtLength(f * r_e.mSamples.mLength));
                        };
                        const auto from_corner = [&](double A) {
                            return mcad_dist(corner, at_fraction(p_lop->mEnd == 1 ? A : 1.0 - A));
                        };
                        // A piece that ends inside the ball is the shrink
                        // rule's, not split here.
                        if (!(from_corner(1.0) > p_lop->mDistance))
                            continue;
                        double lo = 0.0;
                        double hi = 1.0;
                        for (int k = 0; k < 64 && hi - lo > 1e-15; ++k) {
                            const double mid = 0.5 * (lo + hi);
                            (from_corner(mid) < p_lop->mDistance ? lo : hi) = mid;
                        }
                        const double b = 0.5 * (lo + hi);
                        a = p_lop->mEnd == 1 ? b : 1.0 - b;
                        const double d = from_corner(b);
                        auto [it, fresh] = lop_radius.emplace(corner, d);
                        if (!fresh)
                            it->second = std::min(it->second, d);
                        pieces.push_back(piece);
                        at.push_back(a);
                        continue;
                    }
                    if (end == 1 || end == 2) {
                        // A power of two in (L / 3, 2 L / 3] from the corner.
                        const double len =
                            (r_e.mSamples.mFractions[piece + 1] - r_e.mSamples.mFractions[piece]) *
                            r_e.mSamples.mLength;
                        int ex = 0;
                        std::frexp(2.0 * len / 3.0, &ex);
                        const double r = std::ldexp(1.0, ex - 1);
                        if (len > 0.0 && r > len / 3.0 && r <= 2.0 * len / 3.0)
                            a = end == 1 ? r / len : 1.0 - r / len;
                    }
                    pieces.push_back(piece);
                    at.push_back(a);
                }
                detail::edge_refine(*traces[r_e.mMaster], r_e.mSamples, pieces, at);
                carry(e);
            }
            boundary_refinement += static_cast<std::int64_t>(added);
            assign_ids();
            std::vector<char> check(nf, 0);
            for (const auto& [e, r_pieces] : split)
                for (const std::size_t u : edges[e].mUses)
                    check[uses[u].mFace] = 1;
            separate(check);
            mesh_faces(true, false);
        }
        // The rounds stopped: a face that asked for splits none made (or
        // after the last round) is meshed again, smoothed.
        mesh_faces(false, true);
    }

    // Assemble in face order.
    std::size_t np = boundary.size();
    std::vector<std::int64_t> interior_base(nf, 0);
    for (std::size_t f = 0; f < nf; ++f) {
        interior_base[f] = static_cast<std::int64_t>(np);
        if (faces[f].mRecord.mStatus == MeshCadFaceStatus::Ok)
            np += faces[f].mInterior.size();
    }
    // Curved elements (6.6): after the points above, each edge piece's own
    // nodes (on the edge, by arc length, shared by every face and line that
    // runs along it), then each face's (the surface at parameter points).
    const int order = rOptions.mOrder;
    const std::size_t np_linear = np;
    std::vector<std::vector<std::int64_t>> piece_base(edges.size());
    std::vector<Vec3> piece_points;
    std::vector<std::int64_t> high_base(nf, 0);
    if (order > 1) {
        for (std::size_t e = 0; e < edges.size(); ++e) {
            const McadEdge& r_e = edges[e];
            if (r_e.mFail != MeshCadFaceStatus::Ok || r_e.mSamples.mDegenerate)
                continue;
            const std::vector<std::int64_t>& r_ids = uses[r_e.mMaster].mIds;
            const auto pts = detail::edge_piece_points(*traces[r_e.mMaster], r_e.mSamples, order);
            piece_base[e].assign(r_ids.size() - 1, -1);
            for (std::size_t k = 0; k + 1 < r_ids.size(); ++k) {
                if (r_ids[k] == r_ids[k + 1])
                    continue;
                piece_base[e][k] = static_cast<std::int64_t>(np + piece_points.size());
                for (int j = 0; j + 1 < order; ++j)
                    piece_points.push_back(
                        pts[k * static_cast<std::size_t>(order - 1) + static_cast<std::size_t>(j)]);
            }
        }
        np += piece_points.size();
        for (std::size_t f = 0; f < nf; ++f) {
            high_base[f] = static_cast<std::int64_t>(np);
            if (faces[f].mRecord.mStatus == MeshCadFaceStatus::Ok)
                np += faces[f].mHighOrder.size();
        }
    }
    if (static_cast<std::int64_t>(np) > rOptions.mMaxPoints)
        throw std::length_error(std::string(kMcadPrefix) + "the mesh would have more than " +
                                std::to_string(rOptions.mMaxPoints) + " points");
    NDArray points = NDArray::Uninit(DType::Float64, {np, std::size_t{3}});
    double* px = points.As<double>();
    for (std::size_t i = 0; i < boundary.size(); ++i)
        for (int k = 0; k < 3; ++k)
            px[3 * i + static_cast<std::size_t>(k)] = boundary[i][static_cast<std::size_t>(k)];
    for (std::size_t i = 0; i < piece_points.size(); ++i)
        for (std::size_t k = 0; k < 3; ++k)
            px[3 * (np_linear + i) + k] = piece_points[i][k];
    // A curved element's nodes, in triangle6 / triangle10 order, and the
    // normal its Jacobian is signed against.
    std::vector<std::array<std::int64_t, 10>> curved;
    std::vector<Vec3> curved_normal;
    std::vector<std::array<std::int64_t, 3>> tris;
    std::vector<std::int64_t> tri_face;
    std::vector<std::vector<std::int64_t>> face_cells(nf);
    result.mMinAngleDeg = 180.0;
    result.mMinEdge = std::numeric_limits<double>::infinity();
    // A merged pair's records (6.5.2): the triangles, angles and edges of
    // each face's own share of the host's meshing.
    std::vector<char> dropped_edge(edges.size(), 0);
    for (std::size_t f = 0; f < nf; ++f) {
        McadFace& r_h = faces[f];
        if (r_h.mMergedWith < 0)
            continue;
        if (r_h.mRecord.mStatus != MeshCadFaceStatus::Ok) {
            // The sliver shares its host's outcome.
            MeshCadFaceRecord& r_sr = faces[static_cast<std::size_t>(r_h.mMergedWith)].mRecord;
            r_sr.mStatus = r_h.mRecord.mStatus;
            r_sr.mReason = "meshed with face " + std::to_string(f) + ": " + r_h.mRecord.mReason;
            continue;
        }
        for (const std::size_t e : r_h.mSharedEdges)
            dropped_edge[e] = 1;
        MeshCadFaceRecord& r_sr = faces[static_cast<std::size_t>(r_h.mMergedWith)].mRecord;
        std::array<MeshCadFaceRecord*, 2> rec{&r_h.mRecord, &r_sr};
        std::array<std::int64_t, 2> count{0, 0};
        std::array<double, 2> amin{180.0, 180.0}, amax{0.0, 0.0};
        std::array<double, 2> emin{std::numeric_limits<double>::infinity(),
                                   std::numeric_limits<double>::infinity()};
        std::array<double, 2> emax{0.0, 0.0};
        for (std::size_t t = 0; t < r_h.mTriangles.size(); ++t) {
            const std::size_t o =
                t < r_h.mTriOwner.size() ? static_cast<std::size_t>(r_h.mTriOwner[t]) : 0;
            const auto& r_t = r_h.mTriangles[t];
            const auto at = [&](std::int64_t Id) -> const Vec3& {
                return Id >= 0 ? boundary[static_cast<std::size_t>(Id)]
                               : r_h.mInterior[static_cast<std::size_t>(-Id - 1)];
            };
            const auto ang = mcad_angles(at(r_t[0]), at(r_t[1]), at(r_t[2]));
            ++count[o];
            for (std::size_t k = 0; k < 3; ++k) {
                amin[o] = std::min(amin[o], ang[k]);
                amax[o] = std::max(amax[o], ang[k]);
                const double e = mcad_dist(at(r_t[k]), at(r_t[(k + 1) % 3]));
                emin[o] = std::min(emin[o], e);
                emax[o] = std::max(emax[o], e);
            }
        }
        const MeshCadFaceRecord whole = r_h.mRecord;
        for (std::size_t o = 0; o < 2; ++o) {
            if (count[o] == 0) {
                // Its share is empty: the pair's own extremes.
                rec[o]->mMinAngleDeg = whole.mMinAngleDeg;
                rec[o]->mMaxAngleDeg = whole.mMaxAngleDeg;
                rec[o]->mMinEdge = whole.mMinEdge;
                rec[o]->mMaxEdge = whole.mMaxEdge;
                rec[o]->mBoundDeg = whole.mBoundDeg;
                rec[o]->mNumTriangles = 0;
                continue;
            }
            rec[o]->mNumTriangles = count[o];
            rec[o]->mMinAngleDeg = amin[o];
            rec[o]->mMaxAngleDeg = amax[o];
            rec[o]->mMinEdge = emin[o];
            rec[o]->mMaxEdge = emax[o];
        }
        r_sr.mBoundDeg = count[1] > 0 ? amin[1] : whole.mBoundDeg;
        r_sr.mBelowBound =
            rOptions.mMinAngleDeg > 0.0 && count[1] > 0 && amin[1] < rOptions.mMinAngleDeg;
    }
    for (std::size_t f = 0; f < nf; ++f) {
        McadFace& r_f = faces[f];
        for (const auto& [u, r_why] : r_f.mSkipped)
            warn("face " + std::to_string(f) + ": an embedded edge (trim " +
                 std::to_string(uses[u].mTrimIndex) + ") is not meshed: " + r_why);
        if (r_f.mMeanMetric)
            warn("face " + std::to_string(f) +
                 ": folded in the varying metric, so meshed in its mean metric "
                 "(roadmap 6.4.1)");
        if (r_f.mRoute == McadFace::Route::Declined)
            warn("face " + std::to_string(f) + ": not reparametrised (condition number " +
                 std::to_string(r_f.mRouteCond) + "): " + r_f.mRouteWhy + " (roadmap 6.4.4)");
        if (pExtrasReport &&
            (r_f.mRoute == McadFace::Route::Taken || r_f.mRoute == McadFace::Route::Declined)) {
            MeshCadReparamFace rec;
            rec.mFace = static_cast<std::int64_t>(f);
            rec.mConditionNumber = r_f.mRouteCond;
            rec.mConditionNumberDisc = r_f.mRouteCondDisc;
            rec.mTaken = r_f.mRoute == McadFace::Route::Taken;
            rec.mMethod = rec.mTaken ? r_f.mRouteMethod : 0;
            rec.mReason = rec.mTaken ? std::string() : r_f.mRouteWhy;
            pExtrasReport->mReparametrised.push_back(std::move(rec));
        }
        if (r_f.mNumStretched > 0)
            warn("face " + std::to_string(f) + ": " + std::to_string(r_f.mNumStretched) +
                 " triangles stay under min_angle where its parameters stretch the surface "
                 "more than twofold against the face's plane; metric = \"varying\" meshes them "
                 "in the surface's own metric (roadmap 6.4.1)");
        if (r_f.mRecord.mStatus != MeshCadFaceStatus::Ok) {
            ++result.mNumFailed;
            if (r_f.mRecord.mStatus == MeshCadFaceStatus::Empty && r_f.mRecord.mReason.empty())
                r_f.mRecord.mReason = "nothing to mesh";
            warn("face " + std::to_string(f) + " is not meshed (" +
                 mesh_cad_face_status_name(r_f.mRecord.mStatus) + "): " + r_f.mRecord.mReason);
            continue;
        }
        const std::int64_t base = interior_base[f];
        for (std::size_t i = 0; i < r_f.mInterior.size(); ++i)
            for (int k = 0; k < 3; ++k)
                px[3 * (static_cast<std::size_t>(base) + i) + static_cast<std::size_t>(k)] =
                    r_f.mInterior[i][static_cast<std::size_t>(k)];
        if (order > 1)
            for (std::size_t i = 0; i < r_f.mHighOrder.size(); ++i)
                for (std::size_t k = 0; k < 3; ++k)
                    px[3 * (static_cast<std::size_t>(high_base[f]) + i) + k] = r_f.mHighOrder[i][k];
        for (std::size_t t = 0; t < r_f.mTriangles.size(); ++t) {
            auto tri = r_f.mTriangles[t];
            for (std::int64_t& r_v : tri)
                if (r_v < 0)
                    r_v = base + (-r_v - 1);
            // A merged pair's triangle (6.5.2) keeps the face it lies in.
            const std::size_t owner =
                r_f.mMergedWith >= 0 && t < r_f.mTriOwner.size() && r_f.mTriOwner[t] == 1
                    ? static_cast<std::size_t>(r_f.mMergedWith)
                    : f;
            face_cells[owner].push_back(static_cast<std::int64_t>(tris.size()));
            tris.push_back(tri);
            tri_face.push_back(static_cast<std::int64_t>(owner));
            if (order == 1)
                continue;
            std::array<std::int64_t, 10> el{};
            el.fill(-1);
            std::size_t slot = 0;
            for (; slot < 3; ++slot)
                el[slot] = tri[slot];
            for (std::size_t k = 0; k < 3; ++k) {
                const McadFace::CurvedEdge& r_c = r_f.mCurvedEdges[t][k];
                std::array<std::int64_t, 2> ids{-1, -1};
                if (r_c.mEdge >= 0) {
                    const std::int64_t first = piece_base[static_cast<std::size_t>(r_c.mEdge)]
                                                         [static_cast<std::size_t>(r_c.mPiece)];
                    if (first < 0)
                        throw std::logic_error(std::string(kMcadPrefix) +
                                               "internal error: a curved element runs along "
                                               "a collapsed piece");
                    for (int j = 0; j + 1 < order; ++j)
                        ids[static_cast<std::size_t>(j)] = first + j;
                    if (r_c.mReversed && order == 3)
                        std::swap(ids[0], ids[1]);
                } else {
                    for (int j = 0; j + 1 < order; ++j)
                        ids[static_cast<std::size_t>(j)] =
                            high_base[f] + r_c.mNodes[static_cast<std::size_t>(j)];
                }
                for (int j = 0; j + 1 < order; ++j)
                    el[slot++] = ids[static_cast<std::size_t>(j)];
            }
            if (order == 3)
                el[slot] = high_base[f] + r_f.mCurvedCentre[t];
            curved.push_back(el);
            curved_normal.push_back(r_f.mCurvedNormal[t]);
        }
        result.mNumDegenerateDropped += r_f.mRecord.mNumDegenerateDropped;
        result.mNumRefinementPoints += r_f.mRecord.mNumRefinementPoints;
        result.mNumSmoothed += r_f.mRecord.mNumSmoothed;
        result.mNumFacesBelowBound += r_f.mRecord.mBelowBound ? 1 : 0;
        for (std::size_t b = 0; b < kMeshCadAngleBins; ++b)
            result.mAngleHistogram[b] += r_f.mHistogram[b];
        result.mMinAngleDeg = std::min(result.mMinAngleDeg, r_f.mRecord.mMinAngleDeg);
        result.mMaxAngleDeg = std::max(result.mMaxAngleDeg, r_f.mRecord.mMaxAngleDeg);
        result.mMinEdge = std::min(result.mMinEdge, r_f.mRecord.mMinEdge);
        result.mMaxEdge = std::max(result.mMaxEdge, r_f.mRecord.mMaxEdge);
    }
    if (tris.empty()) {
        result.mMinAngleDeg = 0.0;
        result.mMinEdge = 0.0;
    }

    if (order > 1) {
        std::int64_t still = 0;
        result.mNumCurvedStraightened = mcad_straighten(order, curved, curved_normal, px, still);
        if (result.mNumCurvedStraightened > 0)
            warn(std::to_string(result.mNumCurvedStraightened) +
                 " curved elements had a Jacobian that is not positive everywhere and are kept "
                 "straight-sided");
        if (still > 0)
            warn(std::to_string(still) +
                 " elements stay invalid even straight-sided against their face's normal");
    }
    Mesh& r_m = result.mMesh;
    r_m.AssignPoints(std::move(points));
    if (order == 1) {
        NDArray conn = NDArray::Uninit(DType::Int64, {tris.size(), std::size_t{3}});
        std::int64_t* pc = conn.As<std::int64_t>();
        for (std::size_t t = 0; t < tris.size(); ++t)
            for (int k = 0; k < 3; ++k)
                pc[3 * t + static_cast<std::size_t>(k)] = tris[t][static_cast<std::size_t>(k)];
        r_m.AddCellBlock("triangle", std::move(conn));
    } else {
        const auto nn = static_cast<std::size_t>(detail::curved_triangle_num_nodes(order));
        NDArray conn = NDArray::Uninit(DType::Int64, {curved.size(), nn});
        std::int64_t* pc = conn.As<std::int64_t>();
        for (std::size_t t = 0; t < curved.size(); ++t)
            for (std::size_t k = 0; k < nn; ++k)
                pc[nn * t + k] = curved[t][k];
        r_m.AddCellBlock(order == 2 ? "triangle6" : "triangle10", std::move(conn));
    }
    std::vector<std::array<std::int64_t, 4>> segs;
    std::vector<std::int64_t> seg_edge;
    if (rOptions.mEdges) {
        for (std::size_t e = 0; e < edges.size(); ++e) {
            const McadEdge& r_e = edges[e];
            // An edge a merged pair shares (6.5.2) is inside their mesh.
            if (r_e.mFail != MeshCadFaceStatus::Ok || dropped_edge[e])
                continue;
            const std::vector<std::int64_t>& r_ids = uses[r_e.mMaster].mIds;
            for (std::size_t k = 0; k + 1 < r_ids.size(); ++k) {
                if (r_ids[k] == r_ids[k + 1])
                    continue;
                std::array<std::int64_t, 4> seg{r_ids[k], r_ids[k + 1], -1, -1};
                for (int j = 0; j + 1 < order; ++j)
                    seg[static_cast<std::size_t>(2 + j)] = piece_base[e][k] + j;
                segs.push_back(seg);
                seg_edge.push_back(r_e.mBrepEdge);
            }
        }
        // A suppressed corner (6.5.2): its two line cells become one, the
        // first edge's, straight across it.
        if (!suppressed.empty()) {
            std::map<std::int64_t, std::vector<std::size_t>> at_corner;
            for (std::size_t k = 0; k < segs.size(); ++k)
                for (int j = 0; j < 2; ++j)
                    if (segs[k][j] >= 0 &&
                        static_cast<std::size_t>(segs[k][j]) < suppressed.size() &&
                        suppressed[static_cast<std::size_t>(segs[k][j])])
                        at_corner[segs[k][j]].push_back(k);
            std::vector<char> dropped(segs.size(), 0);
            for (const auto& [c, r_k] : at_corner) {
                if (r_k.size() != 2)
                    continue;
                auto& r_a = segs[r_k[0]];
                const auto& r_b = segs[r_k[1]];
                const std::int64_t far = r_b[0] == c ? r_b[1] : r_b[0];
                (r_a[0] == c ? r_a[0] : r_a[1]) = far;
                dropped[r_k[1]] = 1;
            }
            std::vector<std::array<std::int64_t, 4>> kept;
            std::vector<std::int64_t> kept_edge;
            for (std::size_t k = 0; k < segs.size(); ++k)
                if (!dropped[k]) {
                    kept.push_back(segs[k]);
                    kept_edge.push_back(seg_edge[k]);
                }
            segs = std::move(kept);
            seg_edge = std::move(kept_edge);
        }
        const auto nn = static_cast<std::size_t>(order + 1);
        NDArray lines = NDArray::Uninit(DType::Int64, {segs.size(), nn});
        std::int64_t* pl = lines.As<std::int64_t>();
        for (std::size_t s = 0; s < segs.size(); ++s)
            for (std::size_t k = 0; k < nn; ++k)
                pl[nn * s + k] = segs[s][k];
        r_m.AddCellBlock(order == 1 ? "line" : order == 2 ? "line3" : "line4", std::move(lines));
    }
    const auto column = [](const std::vector<std::int64_t>& rV) {
        NDArray a = NDArray::Uninit(DType::Int64, {rV.size()});
        std::copy(rV.begin(), rV.end(), a.As<std::int64_t>());
        return a;
    };
    if (rOptions.mRecordFaceIds) {
        std::vector<NDArray> blocks;
        blocks.push_back(column(tri_face));
        if (rOptions.mEdges)
            blocks.push_back(column(std::vector<std::int64_t>(segs.size(), -1)));
        r_m.AddCellData(kMeshCadFaceName, std::move(blocks));
    }
    if (rOptions.mEdges) {
        std::vector<NDArray> blocks;
        blocks.push_back(column(std::vector<std::int64_t>(tris.size(), -1)));
        blocks.push_back(column(seg_edge));
        r_m.AddCellData(kMeshCadEdgeName, std::move(blocks));
    }
    if (rOptions.mRecordFaceIds)
        for (std::size_t f = 0; f < nf; ++f)
            r_m.AddRegion(Region(kMeshCadFaceRegionPrefix + std::to_string(f), RegionKind::Cell, 2,
                                 -1, column(face_cells[f])));

    // Conformity: every triangle edge's uses.
    std::unordered_map<std::uint64_t, std::array<std::int64_t, 2>> uses_of;  // count, net direction
    uses_of.reserve(3 * tris.size());
    for (const auto& r_t : tris)
        for (int k = 0; k < 3; ++k) {
            const std::int64_t a = r_t[static_cast<std::size_t>(k)];
            const std::int64_t b = r_t[static_cast<std::size_t>((k + 1) % 3)];
            const auto key = (static_cast<std::uint64_t>(std::min(a, b)) << 32) ^
                             static_cast<std::uint64_t>(std::max(a, b));
            auto& r_c = uses_of[key];
            r_c[0] += 1;
            r_c[1] += a < b ? 1 : -1;
        }
    for (const auto& [key, r_c] : uses_of) {
        (void)key;
        if (r_c[0] == 1)
            ++result.mNumOpenEdges;
        else if (r_c[0] > 2)
            ++result.mNumNonManifoldEdges;
        else if (r_c[1] != 0)
            ++result.mNumMisorientedEdges;
    }
    result.mNumRefinementPoints += boundary_refinement;
    result.mNumPoints = static_cast<std::int64_t>(np);
    result.mNumTriangles = static_cast<std::int64_t>(tris.size());
    result.mNumEdgeSegments = static_cast<std::int64_t>(segs.size());
    result.mNumFaces = static_cast<std::int64_t>(nf);
    result.mNumEdges = static_cast<std::int64_t>(edges.size());
    for (McadFace& r_f : faces) {
        result.mFaces.push_back(std::move(r_f.mRecord));
        if (pSamples)
            pSamples->insert(pSamples->end(), r_f.mSizeSamples.begin(), r_f.mSizeSamples.end());
    }
    if (pExtrasReport) {
        // Where a declared tolerance accepted what the global one would not
        // (6.5.1), in (edge, where) order.
        for (std::size_t e = 0; e < edges.size(); ++e)
            if (projection_accepted[e] > 0.0)
                accepted(static_cast<std::int64_t>(e), projection_accepted[e], "projection");
        for (const auto& [key, gap] : tol_accepted) {
            MeshCadToleranceUse use;
            use.mEdge = edges[key.first].mBrepEdge;
            use.mDeclared = key.second == "corner" ? corner_tol[key.first] : edge_tol[key.first];
            use.mGap = gap;
            use.mWhere = key.second;
            pExtrasReport->mToleranceUses.push_back(std::move(use));
        }
        // The slivers found (6.5.2): a merge undone later says why.
        for (MeshCadMergedFace rec : merged_found) {
            if (rec.mMerged) {
                const McadFace& r_h = faces[static_cast<std::size_t>(rec.mHost)];
                if (r_h.mMergedWith != rec.mSliver) {
                    rec.mMerged = false;
                    rec.mReason = r_h.mMergeWhy;
                }
            }
            pExtrasReport->mMergedFaces.push_back(std::move(rec));
        }
        // The corners lopped (6.2.3), by face then point, as each face's
        // last meshing left them.
        for (const McadFace& r_f : faces)
            if (r_f.mRecord.mStatus == MeshCadFaceStatus::Ok)
                for (const auto& r_l : r_f.mLopped)
                    pExtrasReport->mLoppedCorners.push_back(r_l.second);
        std::sort(pExtrasReport->mLoppedCorners.begin(), pExtrasReport->mLoppedCorners.end(),
                  [](const MeshCadLoppedCorner& rA, const MeshCadLoppedCorner& rB) {
                      return rA.mFace != rB.mFace ? rA.mFace < rB.mFace : rA.mPoint < rB.mPoint;
                  });
    }
    if (rOptions.mStrict && result.mNumFailed > 0)
        for (const MeshCadFaceRecord& r_r : result.mFaces)
            if (r_r.mStatus != MeshCadFaceStatus::Ok)
                throw std::runtime_error(std::string(kMcadPrefix) + "face " +
                                         std::to_string(r_r.mFace) + " is not meshed (" +
                                         mesh_cad_face_status_name(r_r.mStatus) +
                                         "): " + r_r.mReason);
    return result;
}

}  // namespace

MeshCadSurfaceResult mesh_cad_surface(const BrepInfo& rBrep,
                                      const MeshCadSurfaceOptions& rOptions) {
    return mesh_cad_surface(rBrep, rOptions, MeshCadSizeSources{}, nullptr);
}

MeshCadSurfaceResult mesh_cad_surface(const BrepInfo& rBrep, const MeshCadSurfaceOptions& rOptions,
                                      const MeshCadSizeSources& rSources,
                                      MeshCadSizeReport* pReport) {
    return mesh_cad_surface(rBrep, rOptions, rSources, MeshCadSurfaceExtras{}, pReport, nullptr);
}

MeshCadSurfaceResult mesh_cad_surface(const BrepInfo& rBrep, const MeshCadSurfaceOptions& rOptions,
                                      const MeshCadSizeSources& rSources,
                                      const MeshCadSurfaceExtras& rExtras,
                                      MeshCadSizeReport* pReport,
                                      MeshCadSurfaceExtrasReport* pExtrasReport) {
    mcad_check_options(rOptions);
    mcad_check_sources(rBrep, rOptions, rSources);
    if (rExtras.mVirtualFaces && rOptions.mOrder > 1)
        mcad_fail("virtual_faces is order 1 only (roadmap 6.5.2)");
    const std::size_t nf = rBrep.mFaces.size();
    if (pReport) {
        *pReport = MeshCadSizeReport{};
        pReport->mLimitedByEdges.resize(nf);
    }
    if (pExtrasReport)
        *pExtrasReport = MeshCadSurfaceExtrasReport{};
    const bool given =
        !rSources.mFaceSizes.empty() || !rSources.mEdgeSizes.empty() || rSources.pBackground;
    std::vector<std::vector<std::int64_t>>* p_chains = pReport ? &pReport->mVirtualEdges : nullptr;
    if (!(rOptions.mCurvatureDeviation > 0.0 || rOptions.mCellsAcrossGaps > 0.0 || given))
        return mcad_surface_run(rBrep, rOptions, nullptr, nullptr, rSources.mVirtualEdges, p_chains,
                                &rExtras, pExtrasReport);
    // A sized run (6.3): a first, uniform meshing at target_size gives the
    // samples (every point on every face it lies on, with the surface's
    // normal and curvature there, and the edges' points); the field they
    // make sizes the second, which is the result.
    MeshCadSurfaceOptions first = rOptions;
    first.mCurvatureDeviation = 0.0;
    first.mCellsAcrossGaps = 0.0;
    first.mGrading = 0.0;
    first.mOrder = 1;
    first.mSmoothIterations = 0;
    first.mRecordFaceIds = true;
    first.mEdges = true;
    first.mStrict = false;
    std::vector<detail::MeshCadSizeSample> samples;
    const MeshCadSurfaceResult uniform = mcad_surface_run(
        rBrep, first, nullptr, &samples, rSources.mVirtualEdges, nullptr, &rExtras, nullptr);
    const Mesh& r_first = uniform.mMesh;
    // The faces at each point of the first mesh, from its triangles; faces
    // are neighbours where both touch a point: a gap is measured to neither
    // a point's own face nor those.
    std::vector<std::vector<std::int64_t>> at(r_first.NumPoints());
    for (std::size_t b = 0; b < r_first.NumCellBlocks(); ++b) {
        const auto cells = r_first.Cells(b);
        if (cells.Type() != "triangle" || !r_first.HasCellData(kMeshCadFaceName))
            continue;
        const NDArray& face = r_first.CellData(kMeshCadFaceName, b);
        const NDArray& conn = cells.Conn();
        for (std::size_t c = 0; c < cells.NumCells(); ++c) {
            const std::int64_t f = detail::read_int(face, c);
            for (std::size_t k = 0; k < 3; ++k)
                at[static_cast<std::size_t>(detail::read_int(conn, 3 * c + k))].push_back(f);
        }
    }
    std::vector<std::vector<std::int64_t>> adjacent(nf);
    for (auto& r_faces : at) {
        std::sort(r_faces.begin(), r_faces.end());
        r_faces.erase(std::unique(r_faces.begin(), r_faces.end()), r_faces.end());
        for (const std::int64_t a : r_faces)
            for (const std::int64_t b : r_faces)
                if (a != b && a >= 0 && static_cast<std::size_t>(a) < nf)
                    adjacent[static_cast<std::size_t>(a)].push_back(b);
    }
    for (auto& r_adj : adjacent) {
        std::sort(r_adj.begin(), r_adj.end());
        r_adj.erase(std::unique(r_adj.begin(), r_adj.end()), r_adj.end());
    }
    // Each edge's points on the first mesh, in the order its line cells run,
    // by BrepEdge index.
    std::map<std::int64_t, std::vector<std::array<std::int64_t, 2>>> edge_cells;
    for (std::size_t b = 0; b < r_first.NumCellBlocks(); ++b) {
        const auto cells = r_first.Cells(b);
        if (cells.Type() != "line" || !r_first.HasCellData(kMeshCadEdgeName))
            continue;
        const NDArray& edge = r_first.CellData(kMeshCadEdgeName, b);
        const NDArray& conn = cells.Conn();
        for (std::size_t c = 0; c < cells.NumCells(); ++c) {
            const std::int64_t e = detail::read_int(edge, c);
            if (e >= 0)
                edge_cells[e].push_back(
                    {detail::read_int(conn, 2 * c), detail::read_int(conn, 2 * c + 1)});
        }
    }
    const double h = rOptions.mTargetSize;
    // At most six halvings below target_size.
    const double floor = h / 64.0;
    std::vector<double> sizes =
        detail::mesh_cad_sample_sizes(samples, r_first, adjacent, h, floor,
                                      rOptions.mCurvatureDeviation, rOptions.mCellsAcrossGaps);
    std::vector<Vec3> xs;
    xs.reserve(samples.size());
    for (const detail::MeshCadSizeSample& r_s : samples)
        xs.push_back(r_s.mX);
    // Given sizes (6.3.1): a face's at every point of the face; an edge's,
    // and its own curvature's, at every point of the edge (added samples).
    std::vector<double> face_size(nf, std::numeric_limits<double>::infinity());
    for (const auto& [f, s] : rSources.mFaceSizes)
        face_size[static_cast<std::size_t>(f)] =
            std::min(face_size[static_cast<std::size_t>(f)], s);
    for (std::size_t i = 0; i < samples.size(); ++i) {
        const std::int64_t f = samples[i].mFace;
        if (f >= 0 && static_cast<std::size_t>(f) < nf)
            sizes[i] = std::min(sizes[i], face_size[static_cast<std::size_t>(f)]);
    }
    std::map<std::int64_t, double> edge_size;
    for (const auto& [e, s] : rSources.mEdgeSizes) {
        auto [it, fresh] = edge_size.emplace(e, s);
        if (!fresh)
            it->second = std::min(it->second, s);
    }
    const NDArray& r_points = r_first.Points();
    const std::size_t dim = r_first.PointDim();
    const auto point_of = [&](std::int64_t P) {
        Vec3 x{0.0, 0.0, 0.0};
        for (std::size_t k = 0; k < std::min<std::size_t>(dim, 3); ++k)
            x[k] = detail::read_double(r_points, static_cast<std::size_t>(P) * dim + k);
        return x;
    };
    for (const auto& [e, r_cells] : edge_cells) {
        const auto given_it = edge_size.find(e);
        if (given_it != edge_size.end())
            for (const auto& r_c : r_cells)
                for (const std::int64_t p : r_c) {
                    xs.push_back(point_of(p));
                    sizes.push_back(given_it->second);
                }
        if (!(rOptions.mCurvatureDeviation > 0.0))
            continue;
        // The edge's own bend: at a point between two of its cells, the
        // circle through it and its two neighbours, whose radius is the
        // curve's where they lie on a circle (6.3.1).
        std::map<std::int64_t, std::vector<std::int64_t>> next;
        for (const auto& r_c : r_cells) {
            next[r_c[0]].push_back(r_c[1]);
            next[r_c[1]].push_back(r_c[0]);
        }
        const double d = rOptions.mCurvatureDeviation;
        for (const auto& [p, r_n] : next) {
            if (r_n.size() != 2)
                continue;
            const Vec3 a = point_of(r_n[0]);
            const Vec3 b = point_of(p);
            const Vec3 c = point_of(r_n[1]);
            const double ab = mcad_dist(a, b), bc = mcad_dist(b, c), ca = mcad_dist(c, a);
            const Vec3 u{b[0] - a[0], b[1] - a[1], b[2] - a[2]};
            const Vec3 v{c[0] - a[0], c[1] - a[1], c[2] - a[2]};
            const Vec3 w{u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2],
                         u[0] * v[1] - u[1] * v[0]};
            const double twice_area = std::sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);
            if (!(twice_area > 0.0))
                continue;
            const double k = 2.0 * twice_area / (ab * bc * ca);  // 1 / circumradius
            const double q = 2.0 * d / k - d * d;
            if (q > 0.0) {
                xs.push_back(b);
                sizes.push_back(2.0 * std::sqrt(q) / detail::kMeshCadSizeSlack);
            }
        }
    }
    for (double& r_s : sizes)
        r_s = std::max(floor, std::min(h, r_s));
    const double slope = rOptions.mGrading > 0.0 ? rOptions.mGrading - 1.0 : 0.5;
    // The first mesh's points are about h apart: each sample's size holds
    // over the patch of surface around it, three quarters of h out.
    detail::MeshCadSizeField field(h, floor, slope, 0.75 * h, xs, sizes);
    // A background's size is read at every point the field is (6.3.1): it
    // is a size function already, so it is taken as given, not sampled.
    std::unique_ptr<detail::MeshCadBackgroundSize> background;
    if (rSources.pBackground) {
        background = std::make_unique<detail::MeshCadBackgroundSize>(*rSources.pBackground,
                                                                     rSources.mBackgroundArray);
        field.SetBackground(background.get());
    }
    // Who decided (6.3.8): an edge of a face given a size, along which the
    // field is finer than that size allows, limits the face.
    if (pReport) {
        // The faces of each edge, by its topology (ids to indices).
        std::unordered_map<std::int64_t, std::int64_t> face_of_id;
        for (std::size_t f = 0; f < nf; ++f)
            face_of_id.emplace(rBrep.mFaces[f].mId, static_cast<std::int64_t>(f));
        for (const auto& [e, r_cells] : edge_cells) {
            std::set<std::int64_t> faces;
            for (const BrepRef& r_ref : rBrep.mEdges[static_cast<std::size_t>(e)].mTopology) {
                const auto it = face_of_id.find(r_ref.mId);
                if (it != face_of_id.end())
                    faces.insert(it->second);
            }
            // The field's median along the edge: a corner a finer face
            // reaches does not make the whole edge fine.
            std::vector<double> along;
            for (const auto& r_c : r_cells)
                for (const std::int64_t p : r_c)
                    along.push_back(field.At(point_of(p)));
            std::nth_element(along.begin(), along.begin() + along.size() / 2, along.end());
            const double finest = along[along.size() / 2];
            for (const std::int64_t f : faces)
                if (f >= 0 && static_cast<std::size_t>(f) < nf &&
                    std::isfinite(face_size[static_cast<std::size_t>(f)]) &&
                    finest * detail::kMeshCadSizeSlack < face_size[static_cast<std::size_t>(f)])
                    pReport->mLimitedByEdges[static_cast<std::size_t>(f)].push_back(e);
        }
    }
    // No sample asks for less than h (a flat model with no thin part): the
    // uniform mesh, byte for byte.
    MeshCadSurfaceResult r = mcad_surface_run(
        rBrep, rOptions, field.NumSamples() > 0 || field.HasBackground() ? &field : nullptr,
        nullptr, rSources.mVirtualEdges, p_chains, &rExtras, pExtrasReport);
    if (pReport)
        for (std::size_t f = 0; f < nf; ++f)
            for (const std::int64_t e : pReport->mLimitedByEdges[f])
                r.mWarnings.push_back("face " + std::to_string(f) +
                                      ": its size is limited by edge " + std::to_string(e));
    return r;
}

Mesh triangulate_polygon_2d(const std::vector<std::vector<std::array<double, 2>>>& rLoops,
                            double TargetSize) {
    TriangulatePolygon2dOptions o;
    o.mTargetSize = TargetSize;
    return triangulate_polygon_2d(rLoops, o);
}

Mesh triangulate_polygon_2d(const std::vector<std::vector<std::array<double, 2>>>& rLoops,
                            const TriangulatePolygon2dOptions& rOptions) {
    const double TargetSize = rOptions.mTargetSize;
    if (!(std::isfinite(TargetSize) && TargetSize >= 0.0))
        throw std::invalid_argument(
            "meshio++: triangulate_polygon_2d: target_size must be finite and not negative");
    if (!(std::isfinite(rOptions.mMinAngleDeg) && rOptions.mMinAngleDeg >= 0.0 &&
          rOptions.mMinAngleDeg < 60.0))
        throw std::invalid_argument(
            "meshio++: triangulate_polygon_2d: min_angle must lie in [0, 60) degrees");
    if (rOptions.mMaxPoints < 0)
        throw std::invalid_argument(
            "meshio++: triangulate_polygon_2d: max_points must not be negative");
    detail::Cdt2dPolygonQuality quality;
    quality.mMinAngleDeg = rOptions.mMinAngleDeg;
    quality.mMaxPoints = rOptions.mMaxPoints;
    const detail::Cdt2dPolygon t =
        detail::cdt2d_triangulate_polygon(rLoops, TargetSize, -1.0, quality);
    if (t.mStatus != detail::Cdt2dStatus::Ok)
        throw std::invalid_argument("meshio++: triangulate_polygon_2d: " + t.mWhy);
    Mesh m;
    NDArray points = NDArray::Uninit(DType::Float64, {t.mPoints.size(), std::size_t{3}});
    double* px = points.As<double>();
    for (std::size_t i = 0; i < t.mPoints.size(); ++i) {
        px[3 * i] = t.mPoints[i][0];
        px[3 * i + 1] = t.mPoints[i][1];
        px[3 * i + 2] = 0.0;
    }
    m.AssignPoints(std::move(points));
    NDArray conn = NDArray::Uninit(DType::Int64, {t.mTriangles.size(), std::size_t{3}});
    std::int64_t* pc = conn.As<std::int64_t>();
    for (std::size_t i = 0; i < t.mTriangles.size(); ++i)
        for (int k = 0; k < 3; ++k)
            pc[3 * i + static_cast<std::size_t>(k)] = t.mTriangles[i][static_cast<std::size_t>(k)];
    m.AddCellBlock("triangle", std::move(conn));
    return m;
}

}  // namespace meshioplusplus

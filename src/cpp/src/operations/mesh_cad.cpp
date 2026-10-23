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
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "../detail/cdt2d.hpp"
#include "../detail/edge_sampler.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/operations/mesh_cad.hpp"
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/region.hpp"

namespace meshioplusplus {

namespace {

constexpr const char* kMcadPrefix = "meshio++: mesh_cad_surface: ";

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
    std::int64_t mEdge = -1;
    bool mSame = true;               ///< runs the way its edge's master does
    std::vector<double> mParams;     ///< along this use, ends included
    std::vector<std::int64_t> mIds;  ///< mesh points along this use, ends included
    MeshCadFaceStatus mFail = MeshCadFaceStatus::Ok;
    std::string mWhy;
};

struct McadEdge {
    std::int64_t mBrepEdge = -1;
    std::vector<std::size_t> mUses;  ///< ascending; the first is the master
    detail::EdgeSamples mSamples;
    std::int64_t mFirstInterior = 0;
    MeshCadFaceStatus mFail = MeshCadFaceStatus::Ok;
    std::string mWhy;
};

struct McadFace {
    bool mUsable = false;
    std::unique_ptr<detail::BrepSurfaceEval> pSurface;
    std::vector<BrepCurve> mSynthetic;
    std::vector<std::vector<std::size_t>> mLoops;  ///< use indices per loop
    MeshCadFaceRecord mRecord;
    std::vector<Vec3> mInterior;
    std::vector<std::array<std::int64_t, 3>> mTriangles;  ///< < 0: -(1 + interior index)
};

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
}

/// Mesh one face in its parameters (see the file comment of mesh_cad.hpp).
void mcad_mesh_face(McadFace& rF, const BrepFace& rBrepFace, const std::vector<McadUse>& rUses,
                    const std::vector<std::optional<detail::TrimTrace>>& rTraces,
                    const std::vector<Vec3>& rBoundary, const MeshCadSurfaceOptions& rO) {
    MeshCadFaceRecord& r_rec = rF.mRecord;
    const auto fail = [&](MeshCadFaceStatus S, const std::string& rWhy) {
        r_rec.mStatus = S;
        r_rec.mReason = rWhy;
        rF.mTriangles.clear();
        rF.mInterior.clear();
    };
    for (const auto& r_loop : rF.mLoops)
        for (const std::size_t u : r_loop)
            if (rUses[u].mFail != MeshCadFaceStatus::Ok)
                return fail(rUses[u].mFail, rUses[u].mWhy);
    const double h = rO.mTargetSize;

    // The loops as closed polylines of (mesh point, parameters).
    struct Entry {
        std::int64_t mId;
        Vec2 mUv;
    };
    std::vector<std::vector<Entry>> loops;
    double lo[2] = {std::numeric_limits<double>::infinity(),
                    std::numeric_limits<double>::infinity()};
    double hi[2] = {-lo[0], -lo[1]};
    std::size_t total = 0;
    for (const auto& r_loop : rF.mLoops) {
        std::vector<Entry> poly;
        for (const std::size_t u : r_loop) {
            const McadUse& r_u = rUses[u];
            const detail::TrimTrace& r_t = *rTraces[u];
            for (std::size_t j = 0; j + 1 < r_u.mIds.size(); ++j) {
                const Vec2 uv = r_t.Uv(r_u.mParams[j]);
                poly.push_back({r_u.mIds[j], uv});
                for (int k = 0; k < 2; ++k) {
                    lo[k] = std::min(lo[k], uv[k]);
                    hi[k] = std::max(hi[k], uv[k]);
                }
            }
        }
        if (!poly.empty()) {
            total += poly.size();
            loops.push_back(std::move(poly));
        }
    }
    if (total < 3)
        return fail(MeshCadFaceStatus::Empty, "its loops have fewer than three points");

    // The parameters through one affine map x = L (u, v), L the Cholesky
    // factor of the mean first fundamental form over the face's box, so the
    // triangulation is Delaunay in (an average of) the surface's own metric:
    // in raw parameters, whose scales can differ by orders of magnitude, a
    // Delaunay triangulation links points far apart on the surface.
    const detail::BrepSurfaceEval& r_s = *rF.pSurface;
    constexpr int kSamples = 9;
    std::vector<std::array<Vec3, 2>> jac;
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
            jac.push_back({su, sv});
            me += su[0] * su[0] + su[1] * su[1] + su[2] * su[2];
            mf += su[0] * sv[0] + su[1] * sv[1] + su[2] * sv[2];
            mg += sv[0] * sv[0] + sv[1] * sv[1] + sv[2] * sv[2];
        }
    const double ns = double(jac.size());
    me /= ns;
    mf /= ns;
    mg /= ns;
    double l11 = 1.0;
    double l12 = 0.0;
    double l22 = 1.0;
    if (me > 0.0 && mg > 0.0 && me - mf * mf / mg > 1e-12 * me) {
        l11 = std::sqrt(me);
        l12 = mf / l11;
        l22 = std::sqrt(mg - l12 * l12);
    }
    const auto to_x = [&](const Vec2& rUv) -> Vec2 {
        return {l11 * rUv[0] + l12 * rUv[1], l22 * rUv[1]};
    };
    const auto to_uv = [&](const Vec2& rX) -> Vec2 {
        const double v = rX[1] / l22;
        return {(rX[0] - l12 * v) / l11, v};
    };
    std::vector<Vec2> xs;
    xs.reserve(total);
    double xlo[2] = {std::numeric_limits<double>::infinity(),
                     std::numeric_limits<double>::infinity()};
    double xhi[2] = {-xlo[0], -xlo[1]};
    for (const auto& r_poly : loops)
        for (const Entry& r_e : r_poly) {
            xs.push_back(to_x(r_e.mUv));
            for (int k = 0; k < 2; ++k) {
                xlo[k] = std::min(xlo[k], xs.back()[static_cast<std::size_t>(k)]);
                xhi[k] = std::max(xhi[k], xs.back()[static_cast<std::size_t>(k)]);
            }
        }
    detail::Cdt2d cdt(xlo[0], xlo[1], xhi[0], xhi[1]);
    const std::vector<std::int64_t> vid = cdt.InsertBatch(xs);
    // Each triangulation vertex is one mesh point: a boundary vertex holds
    // its id, an interior one -(1 + its index in mInterior).
    std::unordered_map<std::int64_t, std::int64_t> mesh_of;
    {
        std::size_t at = 0;
        for (const auto& r_poly : loops)
            for (const Entry& r_e : r_poly) {
                const auto [it, fresh] = mesh_of.emplace(vid[at++], r_e.mId);
                if (!fresh && it->second != r_e.mId)
                    return fail(MeshCadFaceStatus::Degenerate,
                                "two distinct boundary points have the same parameters");
            }
    }
    {
        std::size_t at = 0;
        for (const auto& r_poly : loops) {
            for (std::size_t i = 0; i < r_poly.size(); ++i) {
                const std::int64_t a = vid[at + i];
                const std::int64_t b = vid[at + (i + 1) % r_poly.size()];
                if (a == b)
                    continue;
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
        const auto nrows = static_cast<std::int64_t>(rows);
        const auto ncols = static_cast<std::int64_t>(cols);
        for (std::int64_t j = 1; j <= nrows; ++j) {
            const double y = xlo[1] + double(j) * dy;
            for (std::int64_t i = 0; i <= ncols; ++i) {
                const Vec2 xp{xlo[0] + (double(i) + 0.5 * double(j % 2)) * dx, y};
                if (!(xp[0] > xlo[0] && xp[0] < xhi[0]))
                    continue;
                const Vec2 uv = to_uv(xp);
                Vec3 x{};
                r_s.Eval(uv[0], uv[1], x.data(), nullptr, nullptr);
                if (grid.Near(x, 0.5 * h, 0.4 * h))
                    continue;
                const std::int64_t v_id = cdt.InsertSteiner(xp);
                if (v_id < 0)
                    continue;
                mesh_of[v_id] = -(1 + static_cast<std::int64_t>(rF.mInterior.size()));
                rF.mInterior.push_back(x);
                grid.AddPoint(x);
            }
        }
    }

    // Lift, drop what a pole collapses, orient along the face's normal.
    const auto pos = [&](std::int64_t Id) -> const Vec3& {
        return Id >= 0 ? rBoundary[static_cast<std::size_t>(Id)]
                       : rF.mInterior[static_cast<std::size_t>(-Id - 1)];
    };
    double min_angle = 180.0;
    double max_angle = 0.0;
    double min_edge = std::numeric_limits<double>::infinity();
    double max_edge = 0.0;
    for (const auto& r_t : cdt.InsideTriangles()) {
        std::array<std::int64_t, 3> tri{};
        for (int k = 0; k < 3; ++k)
            tri[k] = mesh_of.at(r_t[k]);
        if (tri[0] == tri[1] || tri[1] == tri[2] || tri[0] == tri[2]) {
            ++r_rec.mNumDegenerateDropped;
            continue;
        }
        if (rBrepFace.mSwappedNormal)
            std::swap(tri[1], tri[2]);
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
    }
    if (rF.mTriangles.empty())
        return fail(MeshCadFaceStatus::Empty, "every triangle collapses on the surface");
    r_rec.mNumTriangles = static_cast<std::int64_t>(rF.mTriangles.size());
    r_rec.mNumInteriorPoints = static_cast<std::int64_t>(rF.mInterior.size());
    r_rec.mMinAngleDeg = min_angle;
    r_rec.mMaxAngleDeg = max_angle;
    r_rec.mMinEdge = min_edge;
    r_rec.mMaxEdge = max_edge;
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

MeshCadSurfaceResult mesh_cad_surface(const BrepInfo& rBrep,
                                      const MeshCadSurfaceOptions& rOptions) {
    mcad_check_options(rOptions);
    if (rBrep.mFaces.empty())
        mcad_fail("the B-rep has no faces");
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
        if (!r_bf.mEmbeddedEdges.empty())
            warn("face " + std::to_string(f) + ": its embedded edges are not meshed");
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
        if (!trace_why[u].empty() && r_f.mUsable) {
            r_f.mUsable = false;
            r_f.mRecord.mStatus = MeshCadFaceStatus::Unsupported;
            r_f.mRecord.mReason = trace_why[u];
        }
    }
    const auto live = [&](std::size_t U) { return faces[uses[U].mFace].mUsable; };

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
            edges.push_back(std::move(ed));
        }
    }

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
        const std::size_t m = r_e.mUses.front();
        const detail::TrimTrace& r_m = *traces[m];
        const bool degenerate = !(r_m.Length() > tol);
        if (degenerate)
            corners.Join(slot(m, 0), slot(m, 1));
        for (std::size_t i = 1; i < r_e.mUses.size(); ++i) {
            const std::size_t o = r_e.mUses[i];
            uses[o].mSame = degenerate || detail::edge_same_direction(r_m, *traces[o]);
            corners.Join(slot(m, 0), slot(o, uses[o].mSame ? 0 : 1));
            corners.Join(slot(m, 1), slot(o, uses[o].mSame ? 1 : 0));
            if (degenerate)
                corners.Join(slot(o, 0), slot(o, 1));
        }
    }
    std::vector<std::int64_t> corner_id(2 * uses.size(), -1);
    std::vector<Vec3> boundary;
    for (std::size_t s = 0; s < 2 * uses.size(); ++s) {
        if (!live(s / 2))
            continue;
        const std::size_t root = corners.Find(s);
        if (corner_id[root] < 0) {
            const detail::TrimTrace& r_t = *traces[s / 2];
            corner_id[root] = static_cast<std::int64_t>(boundary.size());
            boundary.push_back(r_t.Point(s % 2 == 0 ? r_t.Start() : r_t.End()));
        }
        corner_id[s] = corner_id[root];
    }
    double max_gap = 0.0;
    for (std::size_t s = 0; s < 2 * uses.size(); ++s) {
        if (!live(s / 2))
            continue;
        const detail::TrimTrace& r_t = *traces[s / 2];
        const double gap = mcad_dist(r_t.Point(s % 2 == 0 ? r_t.Start() : r_t.End()),
                                     boundary[static_cast<std::size_t>(corner_id[s])]);
        max_gap = std::max(max_gap, gap);
        if (gap > gap_fail && uses[s / 2].mFail == MeshCadFaceStatus::Ok) {
            uses[s / 2].mFail = MeshCadFaceStatus::OpenLoop;
            uses[s / 2].mWhy = "a corner of its loop is " + std::to_string(gap) +
                               " from the corner the topology joins it to";
        }
    }

    // Sample each edge once on its master, then carry the samples to the
    // other uses.
    detail::EdgeSampleOptions so;
    so.mTargetSize = h;
    so.mMaxChordalDeviation = rOptions.mMaxChordalDeviation;
    so.mMinSegments = rOptions.mMinPointsPerEdge + 1;
    so.mTolerance = tol;
    so.mMaxSegments = rOptions.mMaxPoints;
    std::vector<double> edge_gap(edges.size(), 0.0);
    parallel_for(
        edges.size(),
        [&](std::size_t e) {
            McadEdge& r_e = edges[e];
            try {
                const std::size_t m = r_e.mUses.front();
                // Enough pieces that every loop through the edge has three
                // points: a closed edge alone in its loop needs three.
                detail::EdgeSampleOptions eo = so;
                for (const std::size_t u : r_e.mUses)
                    eo.mMinSegments = std::max<std::int64_t>(
                        eo.mMinSegments,
                        static_cast<std::int64_t>((3 + uses[u].mLoopSize - 1) / uses[u].mLoopSize));
                r_e.mSamples = detail::edge_sample(*traces[m], eo);
                uses[m].mParams = r_e.mSamples.mParams;
                for (std::size_t i = 1; i < r_e.mUses.size(); ++i) {
                    const std::size_t o = r_e.mUses[i];
                    const detail::EdgeProjection p =
                        detail::edge_project(*traces[o], r_e.mSamples, uses[o].mSame);
                    uses[o].mParams = p.mParams;
                    edge_gap[e] = std::max(edge_gap[e], p.mMaxGap);
                    if (p.mMaxGap > gap_fail) {
                        uses[o].mFail = MeshCadFaceStatus::ProjectionFailed;
                        uses[o].mWhy = "the trims of one edge are " + std::to_string(p.mMaxGap) +
                                       " apart (are the uses of a closed edge started at "
                                       "different points?)";
                    }
                }
            } catch (const std::length_error& r_x) {
                r_e.mFail = MeshCadFaceStatus::BudgetExceeded;
                r_e.mWhy = mcad_why(r_x);
            } catch (const std::exception& r_x) {
                r_e.mFail = MeshCadFaceStatus::Other;
                r_e.mWhy = mcad_why(r_x);
            }
        },
        1);
    for (const double g : edge_gap)
        max_gap = std::max(max_gap, g);
    result.mMaxEdgeGap = max_gap;
    if (max_gap > 0.01 * h)
        warn("the trims of some edges or corners disagree by up to " + std::to_string(max_gap) +
             " (the mesh shares their points regardless)");

    // Point ids: corners, then each edge's interior samples.
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
            throw std::length_error(std::string(kMcadPrefix) + "the edges alone need more than " +
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

    // Faces, each in its own parameters.
    parallel_for(
        nf,
        [&](std::size_t f) {
            McadFace& r_f = faces[f];
            if (!r_f.mUsable)
                return;
            try {
                mcad_mesh_face(r_f, rBrep.mFaces[f], uses, traces, boundary, rOptions);
            } catch (const std::exception& r_x) {
                r_f.mRecord.mStatus = MeshCadFaceStatus::Other;
                r_f.mRecord.mReason = mcad_why(r_x);
                r_f.mTriangles.clear();
                r_f.mInterior.clear();
            }
        },
        1);

    // Assemble in face order.
    std::size_t np = boundary.size();
    std::vector<std::int64_t> interior_base(nf, 0);
    for (std::size_t f = 0; f < nf; ++f) {
        interior_base[f] = static_cast<std::int64_t>(np);
        if (faces[f].mRecord.mStatus == MeshCadFaceStatus::Ok)
            np += faces[f].mInterior.size();
    }
    if (static_cast<std::int64_t>(np) > rOptions.mMaxPoints)
        throw std::length_error(std::string(kMcadPrefix) + "the mesh would have more than " +
                                std::to_string(rOptions.mMaxPoints) + " points");
    NDArray points = NDArray::Uninit(DType::Float64, {np, std::size_t{3}});
    double* px = points.As<double>();
    for (std::size_t i = 0; i < boundary.size(); ++i)
        for (int k = 0; k < 3; ++k)
            px[3 * i + static_cast<std::size_t>(k)] = boundary[i][static_cast<std::size_t>(k)];
    std::vector<std::array<std::int64_t, 3>> tris;
    std::vector<std::int64_t> tri_face;
    std::vector<std::vector<std::int64_t>> face_cells(nf);
    result.mMinAngleDeg = 180.0;
    result.mMinEdge = std::numeric_limits<double>::infinity();
    for (std::size_t f = 0; f < nf; ++f) {
        McadFace& r_f = faces[f];
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
        for (auto tri : r_f.mTriangles) {
            for (std::int64_t& r_v : tri)
                if (r_v < 0)
                    r_v = base + (-r_v - 1);
            face_cells[f].push_back(static_cast<std::int64_t>(tris.size()));
            tris.push_back(tri);
            tri_face.push_back(static_cast<std::int64_t>(f));
        }
        result.mNumDegenerateDropped += r_f.mRecord.mNumDegenerateDropped;
        result.mMinAngleDeg = std::min(result.mMinAngleDeg, r_f.mRecord.mMinAngleDeg);
        result.mMaxAngleDeg = std::max(result.mMaxAngleDeg, r_f.mRecord.mMaxAngleDeg);
        result.mMinEdge = std::min(result.mMinEdge, r_f.mRecord.mMinEdge);
        result.mMaxEdge = std::max(result.mMaxEdge, r_f.mRecord.mMaxEdge);
    }
    if (tris.empty()) {
        result.mMinAngleDeg = 0.0;
        result.mMinEdge = 0.0;
    }

    Mesh& r_m = result.mMesh;
    r_m.AssignPoints(std::move(points));
    NDArray conn = NDArray::Uninit(DType::Int64, {tris.size(), std::size_t{3}});
    std::int64_t* pc = conn.As<std::int64_t>();
    for (std::size_t t = 0; t < tris.size(); ++t)
        for (int k = 0; k < 3; ++k)
            pc[3 * t + static_cast<std::size_t>(k)] = tris[t][static_cast<std::size_t>(k)];
    r_m.AddCellBlock("triangle", std::move(conn));
    std::vector<std::array<std::int64_t, 2>> segs;
    std::vector<std::int64_t> seg_edge;
    if (rOptions.mEdges) {
        for (const McadEdge& r_e : edges) {
            if (r_e.mFail != MeshCadFaceStatus::Ok)
                continue;
            const std::vector<std::int64_t>& r_ids = uses[r_e.mUses.front()].mIds;
            for (std::size_t k = 0; k + 1 < r_ids.size(); ++k) {
                if (r_ids[k] == r_ids[k + 1])
                    continue;
                segs.push_back({r_ids[k], r_ids[k + 1]});
                seg_edge.push_back(r_e.mBrepEdge);
            }
        }
        NDArray lines = NDArray::Uninit(DType::Int64, {segs.size(), std::size_t{2}});
        std::int64_t* pl = lines.As<std::int64_t>();
        for (std::size_t s = 0; s < segs.size(); ++s) {
            pl[2 * s] = segs[s][0];
            pl[2 * s + 1] = segs[s][1];
        }
        r_m.AddCellBlock("line", std::move(lines));
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
    result.mNumPoints = static_cast<std::int64_t>(np);
    result.mNumTriangles = static_cast<std::int64_t>(tris.size());
    result.mNumEdgeSegments = static_cast<std::int64_t>(segs.size());
    result.mNumFaces = static_cast<std::int64_t>(nf);
    result.mNumEdges = static_cast<std::int64_t>(edges.size());
    for (McadFace& r_f : faces)
        result.mFaces.push_back(std::move(r_f.mRecord));
    if (rOptions.mStrict && result.mNumFailed > 0)
        for (const MeshCadFaceRecord& r_r : result.mFaces)
            if (r_r.mStatus != MeshCadFaceStatus::Ok)
                throw std::runtime_error(std::string(kMcadPrefix) + "face " +
                                         std::to_string(r_r.mFace) + " is not meshed (" +
                                         mesh_cad_face_status_name(r_r.mStatus) +
                                         "): " + r_r.mReason);
    return result;
}

Mesh triangulate_polygon_2d(const std::vector<std::vector<std::array<double, 2>>>& rLoops,
                            double TargetSize) {
    if (!(std::isfinite(TargetSize) && TargetSize >= 0.0))
        throw std::invalid_argument(
            "meshio++: triangulate_polygon_2d: target_size must be finite and not negative");
    const detail::Cdt2dPolygon t = detail::cdt2d_triangulate_polygon(rLoops, TargetSize);
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

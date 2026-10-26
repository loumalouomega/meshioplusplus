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
// The constrained Delaunay triangulation of detail/cdt2d.hpp.
// Anonymous-namespace helpers are prefixed `cdt2d_` (the amalgamation
// concatenates every translation unit into one).

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "cdt2d.hpp"
#include "predicates.hpp"
#include "meshioplusplus/detail/space_filling.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

constexpr int cdt2d_next(int I) {
    return I == 2 ? 0 : I + 1;
}
constexpr int cdt2d_prev(int I) {
    return I == 0 ? 2 : I - 1;
}

[[noreturn]] void cdt2d_internal(const char* pWhat) {
    throw std::logic_error(std::string("meshio++: cdt2d: internal error: ") + pWhat);
}

/// A uniform hash grid of 2-D points for the merge and spacing tests.
class Cdt2dGrid {
public:
    Cdt2dGrid(double MinU, double MinV, double Cell) : mU0(MinU), mV0(MinV), mCell(Cell) {}
    std::int64_t Key(double U, double V) const {
        const auto i = static_cast<std::int64_t>(std::floor((U - mU0) / mCell));
        const auto j = static_cast<std::int64_t>(std::floor((V - mV0) / mCell));
        return (i << 32) ^ (j & 0xffffffffll);
    }
    std::pair<std::int64_t, std::int64_t> Cell(double U, double V) const {
        return {static_cast<std::int64_t>(std::floor((U - mU0) / mCell)),
                static_cast<std::int64_t>(std::floor((V - mV0) / mCell))};
    }
    static std::int64_t Join(std::int64_t I, std::int64_t J) {
        return (I << 32) ^ (J & 0xffffffffll);
    }
    void Add(double U, double V, std::int64_t Id) { mCells[Key(U, V)].push_back(Id); }
    const std::vector<std::int64_t>* At(std::int64_t K) const {
        const auto it = mCells.find(K);
        return it == mCells.end() ? nullptr : &it->second;
    }

private:
    double mU0;
    double mV0;
    double mCell;
    std::unordered_map<std::int64_t, std::vector<std::int64_t>> mCells;
};

double cdt2d_segment_distance2(const Cdt2d::Vec2& rP, const Cdt2d::Vec2& rA,
                               const Cdt2d::Vec2& rB) {
    const double dx = rB[0] - rA[0];
    const double dy = rB[1] - rA[1];
    const double len2 = dx * dx + dy * dy;
    double t = 0.0;
    if (len2 > 0.0)
        t = std::clamp(((rP[0] - rA[0]) * dx + (rP[1] - rA[1]) * dy) / len2, 0.0, 1.0);
    const double ex = rA[0] + t * dx - rP[0];
    const double ey = rA[1] + t * dy - rP[1];
    return ex * ex + ey * ey;
}

/// Whether `rP`, collinear with the segment (A, B), lies on its closed box.
bool cdt2d_in_box(const Cdt2d::Vec2& rA, const Cdt2d::Vec2& rB, const Cdt2d::Vec2& rP) {
    return std::min(rA[0], rB[0]) <= rP[0] && rP[0] <= std::max(rA[0], rB[0]) &&
           std::min(rA[1], rB[1]) <= rP[1] && rP[1] <= std::max(rA[1], rB[1]);
}

/// Whether D, collinear with A and B and distinct from A, lies on the ray
/// from A through B: exact comparisons of coordinates, no subtraction.
bool cdt2d_same_ray(const Cdt2d::Vec2& rA, const Cdt2d::Vec2& rB, const Cdt2d::Vec2& rD) {
    const int k = rA[0] != rB[0] ? 0 : 1;
    if (rD[k] == rA[k])
        return false;
    return (rB[k] > rA[k]) == (rD[k] > rA[k]);
}

/// Whether two segments (indices into rPoints) meet anywhere but at an
/// endpoint they share (see cdt2d_find_crossings).
bool cdt2d_segments_meet(const std::vector<Cdt2d::Vec2>& rPoints,
                         const std::array<std::int64_t, 2>& rS,
                         const std::array<std::int64_t, 2>& rT) {
    const auto at = [&](std::int64_t I) -> const Cdt2d::Vec2& {
        return rPoints[static_cast<std::size_t>(I)];
    };
    const Cdt2d::Vec2& a = at(rS[0]);
    const Cdt2d::Vec2& b = at(rS[1]);
    const Cdt2d::Vec2& c = at(rT[0]);
    const Cdt2d::Vec2& d = at(rT[1]);
    if (a == b || c == d)
        return false;
    const int shared =
        int(rS[0] == rT[0]) + int(rS[0] == rT[1]) + int(rS[1] == rT[0]) + int(rS[1] == rT[1]);
    if (shared >= 2)
        return true;
    if (shared == 1) {
        // Only a fold-back: the other ends on one ray from the shared one.
        const std::int64_t p = (rS[0] == rT[0] || rS[0] == rT[1]) ? rS[0] : rS[1];
        const Cdt2d::Vec2& q = at(rS[0] == p ? rS[1] : rS[0]);
        const Cdt2d::Vec2& r = at(rT[0] == p ? rT[1] : rT[0]);
        return orient2d(at(p).data(), q.data(), r.data()) == 0.0 && cdt2d_same_ray(at(p), q, r);
    }
    const double o1 = orient2d(a.data(), b.data(), c.data());
    const double o2 = orient2d(a.data(), b.data(), d.data());
    const double o3 = orient2d(c.data(), d.data(), a.data());
    const double o4 = orient2d(c.data(), d.data(), b.data());
    if (((o1 > 0.0 && o2 < 0.0) || (o1 < 0.0 && o2 > 0.0)) &&
        ((o3 > 0.0 && o4 < 0.0) || (o3 < 0.0 && o4 > 0.0)))
        return true;
    return (o1 == 0.0 && cdt2d_in_box(a, b, c)) || (o2 == 0.0 && cdt2d_in_box(a, b, d)) ||
           (o3 == 0.0 && cdt2d_in_box(c, d, a)) || (o4 == 0.0 && cdt2d_in_box(c, d, b));
}

}  // namespace

bool cdt2d_in_range(double X) {
    const double a = std::fabs(X);
    return a == 0.0 || (a >= 0x1p-128 && a <= 0x1p128);
}

const char* cdt2d_status_name(Cdt2dStatus Status) {
    switch (Status) {
        case Cdt2dStatus::Ok:
            return "ok";
        case Cdt2dStatus::SelfIntersecting:
            return "self-intersecting boundary";
        case Cdt2dStatus::OpenLoop:
            return "open boundary loop";
        case Cdt2dStatus::Degenerate:
            return "degenerate domain";
    }
    return "unknown";
}

Cdt2d::Cdt2d(double MinU, double MinV, double MaxU, double MaxV) {
    const double cu = 0.5 * (MinU + MaxU);
    const double cv = 0.5 * (MinV + MaxV);
    double m = std::max(MaxU - MinU, MaxV - MinV);
    if (!(m > 0.0) || !std::isfinite(m))
        m = std::max({1.0, std::fabs(cu), std::fabs(cv)});
    const double r = 32.0 * m;
    mPoints = {{cu - r, cv - r}, {cu + r, cv - r}, {cu, cv + r}};
    Triangle t;
    t.mV = {0, 1, 2};
    t.mN = {-1, -1, -1};
    t.mC = {false, false, false};
    mTris.push_back(t);
    mVertTri = {0, 0, 0};
}

double Cdt2d::Orient(std::int64_t A, std::int64_t B, std::int64_t C) const {
    return orient2d(Point(A).data(), Point(B).data(), Point(C).data());
}

double Cdt2d::Orient(std::int64_t A, std::int64_t B, const Vec2& rC) const {
    return orient2d(Point(A).data(), Point(B).data(), rC.data());
}

double Cdt2d::InCircle(std::int64_t A, std::int64_t B, std::int64_t C, std::int64_t D) const {
    return incircle(Point(A).data(), Point(B).data(), Point(C).data(), Point(D).data());
}

bool Cdt2d::FlipDue(std::int64_t A, std::int64_t B, std::int64_t C, std::int64_t D) const {
    // Without a metric, or with the super triangle's vertices (outside the
    // domain, far away), the plane's test.
    if (!mMetric || A < kFirstVertex || B < kFirstVertex || C < kFirstVertex || D < kFirstVertex)
        return InCircle(A, B, C, D) > 0.0;
    const Vec2& pb = Point(B);
    const Vec2& pc = Point(C);
    std::array<Vec2, 4> m{Point(A), pb, pc, Point(D)};
    if (!MetricMapped({0.5 * (pb[0] + pc[0]), 0.5 * (pb[1] + pc[1])}, m))
        return InCircle(A, B, C, D) > 0.0;
    // Lawson's flip of (B, C) to (A, D) needs the quadrilateral A, B, D, C
    // strictly convex: the new edge crosses the old one.
    if (!(Orient(A, B, D) > 0.0) || !(Orient(A, D, C) > 0.0))
        return false;
    return incircle(m[0].data(), m[1].data(), m[2].data(), m[3].data()) > 0.0;
}

bool Cdt2d::MetricMapped(const Vec2& rAt, std::array<Vec2, 4>& rPoints, std::size_t N) const {
    if (!mMetric)
        return false;
    const std::array<double, 3> m = mMetric(rAt);
    // M = R^T R, R = [[r11, r12], [0, r22]] upper triangular: det R > 0, so
    // the mapping keeps every orientation.
    const double r11 = m[0] > 0.0 ? std::sqrt(m[0]) : 0.0;
    const double r12 = r11 > 0.0 ? m[1] / r11 : 0.0;
    const double q = m[2] - r12 * r12;
    const double r22 = q > 0.0 ? std::sqrt(q) : 0.0;
    if (!(r11 > 0.0 && r22 > 0.0 && std::isfinite(r11) && std::isfinite(r12) && std::isfinite(r22)))
        return false;
    std::array<Vec2, 4> out{};
    for (std::size_t k = 0; k < N; ++k) {
        const double x = rPoints[k][0] - rAt[0], y = rPoints[k][1] - rAt[1];
        out[k] = {r11 * x + r12 * y, r22 * y};
        if (!cdt2d_in_range(out[k][0]) || !cdt2d_in_range(out[k][1]))
            return false;
    }
    for (std::size_t k = 0; k < N; ++k)
        rPoints[k] = out[k];
    return true;
}

int Cdt2d::IndexOf(std::int64_t T, std::int64_t V) const {
    const auto& r_t = mTris[static_cast<std::size_t>(T)];
    for (int i = 0; i < 3; ++i)
        if (r_t.mV[i] == V)
            return i;
    return -1;
}

int Cdt2d::NeighbourIndex(std::int64_t T, std::int64_t Of) const {
    const auto& r_t = mTris[static_cast<std::size_t>(T)];
    for (int i = 0; i < 3; ++i)
        if (r_t.mN[i] == Of)
            return i;
    return -1;
}

void Cdt2d::SetConstrained(std::int64_t T, int E) {
    auto& r_t = mTris[static_cast<std::size_t>(T)];
    r_t.mC[E] = true;
    const std::int64_t u = r_t.mN[E];
    if (u >= 0) {
        const int f = NeighbourIndex(u, T);
        if (f < 0)
            cdt2d_internal("broken adjacency");
        mTris[static_cast<std::size_t>(u)].mC[f] = true;
    }
}

Cdt2d::Location Cdt2d::LocateBruteForce(const Vec2& rP) const {
    for (std::size_t t = 0; t < mTris.size(); ++t) {
        const auto& r_t = mTris[t];
        if (!r_t.mAlive)
            continue;
        double o[3];
        bool in = true;
        for (int e = 0; e < 3 && in; ++e) {
            o[e] = Orient(r_t.mV[cdt2d_next(e)], r_t.mV[cdt2d_prev(e)], rP);
            in = o[e] >= 0.0;
        }
        if (!in)
            continue;
        Location loc;
        loc.mTri = static_cast<std::int64_t>(t);
        int zeros = 0;
        for (int e = 0; e < 3; ++e)
            zeros += o[e] == 0.0;
        if (zeros == 0) {
            loc.mKind = 0;
        } else if (zeros == 1) {
            loc.mKind = 1;
            for (int e = 0; e < 3; ++e)
                if (o[e] == 0.0)
                    loc.mIndex = e;
        } else {
            loc.mKind = 2;
            for (int e = 0; e < 3; ++e)
                if (o[e] != 0.0)
                    loc.mVertex = r_t.mV[e];
        }
        return loc;
    }
    cdt2d_internal("point outside the super triangle");
}

Cdt2d::Location Cdt2d::Locate(const Vec2& rP) {
    std::int64_t t = mLast;
    if (t < 0 || t >= static_cast<std::int64_t>(mTris.size()) ||
        !mTris[static_cast<std::size_t>(t)].mAlive) {
        t = 0;
        while (!mTris[static_cast<std::size_t>(t)].mAlive)
            ++t;
    }
    const std::size_t budget = 4 * mTris.size() + 64;
    for (std::size_t step = 0; step < budget; ++step) {
        const auto& r_t = mTris[static_cast<std::size_t>(t)];
        // A deterministic pseudo-random first edge (a stochastic walk), so the
        // walk cannot cycle in a triangulation that is not Delaunay.
        mWalkState ^= mWalkState << 13;
        mWalkState ^= mWalkState >> 7;
        mWalkState ^= mWalkState << 17;
        const int first = static_cast<int>(mWalkState % 3);
        bool moved = false;
        double o[3];
        for (int k = 0; k < 3; ++k) {
            const int e = (first + k) % 3;
            o[e] = Orient(r_t.mV[cdt2d_next(e)], r_t.mV[cdt2d_prev(e)], rP);
            if (o[e] < 0.0) {
                if (r_t.mN[e] < 0)
                    cdt2d_internal("point outside the super triangle");
                t = r_t.mN[e];
                moved = true;
                break;
            }
        }
        if (moved)
            continue;
        Location loc;
        loc.mTri = t;
        int zeros = 0;
        for (int e = 0; e < 3; ++e)
            zeros += o[e] == 0.0;
        if (zeros == 0) {
            loc.mKind = 0;
        } else if (zeros == 1) {
            loc.mKind = 1;
            for (int e = 0; e < 3; ++e)
                if (o[e] == 0.0)
                    loc.mIndex = e;
        } else {
            loc.mKind = 2;
            for (int e = 0; e < 3; ++e)
                if (o[e] != 0.0)
                    loc.mVertex = r_t.mV[e];
        }
        return loc;
    }
    return LocateBruteForce(rP);
}

std::vector<std::int64_t> Cdt2d::ReplaceCavity(const std::vector<std::int64_t>& rRemoved,
                                               const std::vector<Tri>& rNew, bool Inside) {
    // The cavity's boundary: the removed triangles' edges whose neighbour
    // was not removed, with that neighbour and the edge's constraint.
    std::vector<Boundary> boundary;
    boundary.reserve(3 * rRemoved.size());
    for (const std::int64_t r : rRemoved) {
        const auto& r_t = mTris[static_cast<std::size_t>(r)];
        for (int e = 0; e < 3; ++e) {
            const std::int64_t nb = r_t.mN[e];
            if (nb >= 0 && std::find(rRemoved.begin(), rRemoved.end(), nb) != rRemoved.end())
                continue;
            boundary.push_back({r_t.mV[cdt2d_next(e)], r_t.mV[cdt2d_prev(e)], nb, r_t.mC[e]});
        }
    }
    for (const std::int64_t r : rRemoved)
        mTris[static_cast<std::size_t>(r)].mAlive = false;
    std::vector<std::int64_t> ids;
    ids.reserve(rNew.size());
    std::size_t reuse = 0;
    for (const Tri& r_v : rNew) {
        std::int64_t id = -1;
        if (reuse < rRemoved.size()) {
            id = rRemoved[reuse++];
        } else if (!mFree.empty()) {
            id = mFree.back();
            mFree.pop_back();
        } else {
            id = static_cast<std::int64_t>(mTris.size());
            mTris.emplace_back();
        }
        auto& r_t = mTris[static_cast<std::size_t>(id)];
        r_t.mV = r_v;
        r_t.mN = {-1, -1, -1};
        r_t.mC = {false, false, false};
        r_t.mInside = Inside;
        r_t.mAlive = true;
        ids.push_back(id);
    }
    for (; reuse < rRemoved.size(); ++reuse)
        mFree.push_back(rRemoved[reuse]);
    // Link the new triangles to each other (twin edges) and to the outside.
    for (std::size_t a = 0; a < ids.size(); ++a) {
        auto& r_t = mTris[static_cast<std::size_t>(ids[a])];
        for (int e = 0; e < 3; ++e) {
            if (r_t.mN[e] >= 0)
                continue;
            const std::int64_t x = r_t.mV[cdt2d_next(e)];
            const std::int64_t y = r_t.mV[cdt2d_prev(e)];
            bool linked = false;
            for (std::size_t b = 0; b < ids.size() && !linked; ++b) {
                if (b == a)
                    continue;
                auto& r_u = mTris[static_cast<std::size_t>(ids[b])];
                for (int f = 0; f < 3; ++f) {
                    if (r_u.mV[cdt2d_next(f)] == y && r_u.mV[cdt2d_prev(f)] == x) {
                        r_t.mN[e] = ids[b];
                        r_u.mN[f] = ids[a];
                        linked = true;
                        break;
                    }
                }
            }
            if (linked)
                continue;
            for (const Boundary& r_b : boundary) {
                if (r_b.mFrom != x || r_b.mTo != y)
                    continue;
                r_t.mN[e] = r_b.mOut;
                r_t.mC[e] = r_b.mConstrained;
                if (r_b.mOut >= 0) {
                    auto& r_o = mTris[static_cast<std::size_t>(r_b.mOut)];
                    int g = -1;
                    for (int k = 0; k < 3; ++k)
                        if (r_o.mV[cdt2d_next(k)] == y && r_o.mV[cdt2d_prev(k)] == x)
                            g = k;
                    if (g < 0)
                        cdt2d_internal("cavity boundary does not match");
                    r_o.mN[g] = ids[a];
                }
                linked = true;
                break;
            }
            if (!linked)
                cdt2d_internal("cavity edge without a match");
        }
        for (int i = 0; i < 3; ++i)
            mVertTri[static_cast<std::size_t>(r_t.mV[i])] = ids[a];
    }
    if (!ids.empty())
        mLast = ids.back();
    return ids;
}

void Cdt2d::Flip(std::int64_t T, int E) {
    const std::int64_t u = mTris[static_cast<std::size_t>(T)].mN[E];
    const int f = NeighbourIndex(u, T);
    if (f < 0)
        cdt2d_internal("broken adjacency in a flip");
    Triangle& r_t = mTris[static_cast<std::size_t>(T)];
    Triangle& r_u = mTris[static_cast<std::size_t>(u)];
    // t = (a, b, c) with the flipped edge (b, c) opposite a; u = (d, c, b).
    const std::int64_t a = r_t.mV[E];
    const std::int64_t b = r_t.mV[cdt2d_next(E)];
    const std::int64_t c = r_t.mV[cdt2d_prev(E)];
    const std::int64_t d = r_u.mV[f];
    const std::int64_t n_ca = r_t.mN[cdt2d_next(E)];
    const bool c_ca = r_t.mC[cdt2d_next(E)];
    const std::int64_t n_ab = r_t.mN[cdt2d_prev(E)];
    const bool c_ab = r_t.mC[cdt2d_prev(E)];
    const std::int64_t n_bd = r_u.mN[cdt2d_next(f)];
    const bool c_bd = r_u.mC[cdt2d_next(f)];
    const std::int64_t n_dc = r_u.mN[cdt2d_prev(f)];
    const bool c_dc = r_u.mC[cdt2d_prev(f)];
    // t = (a, b, d), u = (a, d, c).
    r_t.mV = {a, b, d};
    r_t.mN = {n_bd, u, n_ab};
    r_t.mC = {c_bd, false, c_ab};
    r_u.mV = {a, d, c};
    r_u.mN = {n_dc, n_ca, T};
    r_u.mC = {c_dc, c_ca, false};
    if (n_bd >= 0) {
        const int g = NeighbourIndex(n_bd, u);
        if (g < 0)
            cdt2d_internal("broken adjacency in a flip");
        mTris[static_cast<std::size_t>(n_bd)].mN[g] = T;
    }
    if (n_ca >= 0) {
        const int g = NeighbourIndex(n_ca, T);
        if (g < 0)
            cdt2d_internal("broken adjacency in a flip");
        mTris[static_cast<std::size_t>(n_ca)].mN[g] = u;
    }
    mVertTri[static_cast<std::size_t>(a)] = T;
    mVertTri[static_cast<std::size_t>(b)] = T;
    mVertTri[static_cast<std::size_t>(d)] = T;
    mVertTri[static_cast<std::size_t>(c)] = u;
}

void Cdt2d::Legalize(std::vector<std::pair<std::int64_t, int>>& rStack) {
    while (!rStack.empty()) {
        const auto [t, e] = rStack.back();
        rStack.pop_back();
        const auto& r_t = mTris[static_cast<std::size_t>(t)];
        if (!r_t.mAlive || r_t.mC[e])
            continue;
        const std::int64_t u = r_t.mN[e];
        if (u < 0)
            continue;
        const int f = NeighbourIndex(u, t);
        if (f < 0)
            cdt2d_internal("broken adjacency while legalising");
        const std::int64_t d = mTris[static_cast<std::size_t>(u)].mV[f];
        if (FlipDue(r_t.mV[e], r_t.mV[cdt2d_next(e)], r_t.mV[cdt2d_prev(e)], d)) {
            Flip(t, e);
            rStack.emplace_back(t, 0);
            rStack.emplace_back(u, 0);
        }
    }
}

std::int64_t Cdt2d::InsertAt(const Vec2& rP, const Location& rLoc, bool Inside) {
    const auto v = static_cast<std::int64_t>(mPoints.size());
    mPoints.push_back(rP);
    mVertTri.push_back(rLoc.mTri);
    // Reserved for the common case: allocations dominated insertion.
    std::vector<std::int64_t> removed;
    removed.reserve(2);
    removed.push_back(rLoc.mTri);
    std::vector<Tri> fresh;
    fresh.reserve(4);
    std::int64_t split_a = -1;
    std::int64_t split_b = -1;
    bool split_constrained = false;
    const auto& r_t = mTris[static_cast<std::size_t>(rLoc.mTri)];
    // A constrained edge split after the classification divides inside
    // from outside (or a hole): each side keeps its own flag.
    bool side_inside[2] = {Inside, Inside};
    if (rLoc.mKind == 0) {
        for (int e = 0; e < 3; ++e)
            fresh.push_back({r_t.mV[cdt2d_next(e)], r_t.mV[cdt2d_prev(e)], v});
    } else {
        const int e = rLoc.mIndex;
        split_a = r_t.mV[cdt2d_next(e)];
        split_b = r_t.mV[cdt2d_prev(e)];
        split_constrained = r_t.mC[e];
        const std::int64_t u = r_t.mN[e];
        if (mClassified && split_constrained) {
            side_inside[0] = r_t.mInside;
            side_inside[1] = u >= 0 && mTris[static_cast<std::size_t>(u)].mInside;
        }
        // The two edges of t that stay, then u's.
        fresh.push_back({r_t.mV[e], split_a, v});
        fresh.push_back({split_b, r_t.mV[e], v});
        if (u >= 0) {
            removed.push_back(u);
            const int f = NeighbourIndex(u, rLoc.mTri);
            const auto& r_u = mTris[static_cast<std::size_t>(u)];
            fresh.push_back({r_u.mV[f], split_b, v});
            fresh.push_back({split_a, r_u.mV[f], v});
        } else {
            fresh.push_back({split_a, split_b, v});  // never: points stay inside the super triangle
        }
    }
    const std::vector<std::int64_t> ids = ReplaceCavity(removed, fresh, Inside);
    if (mClassified && split_constrained)
        for (std::size_t k = 0; k < ids.size(); ++k)
            mTris[static_cast<std::size_t>(ids[k])].mInside = side_inside[k < 2 ? 0 : 1];
    if (split_constrained) {
        for (const std::int64_t id : ids) {
            for (int e = 0; e < 3; ++e) {
                const auto& r_n = mTris[static_cast<std::size_t>(id)];
                const std::int64_t x = r_n.mV[cdt2d_next(e)];
                const std::int64_t y = r_n.mV[cdt2d_prev(e)];
                if ((x == v && (y == split_a || y == split_b)) ||
                    (y == v && (x == split_a || x == split_b)))
                    SetConstrained(id, e);
            }
        }
    }
    std::vector<std::pair<std::int64_t, int>> stack;
    stack.reserve(16);
    for (const std::int64_t id : ids)
        stack.emplace_back(id, 2);
    Legalize(stack);
    mLast = mVertTri[static_cast<std::size_t>(v)];
    return v;
}

std::int64_t Cdt2d::Insert(const Vec2& rP) {
    const Location loc = Locate(rP);
    if (loc.mKind == 2)
        return loc.mVertex;
    return InsertAt(rP, loc,
                    mClassified ? mTris[static_cast<std::size_t>(loc.mTri)].mInside : false);
}

std::vector<std::int64_t> Cdt2d::InsertBatch(const std::vector<Vec2>& rPoints) {
    std::vector<std::int64_t> ids(rPoints.size(), -1);
    if (rPoints.empty())
        return ids;
    double lo[2] = {rPoints[0][0], rPoints[0][1]};
    double hi[2] = {lo[0], lo[1]};
    for (const Vec2& r_p : rPoints)
        for (int k = 0; k < 2; ++k) {
            lo[k] = std::min(lo[k], r_p[k]);
            hi[k] = std::max(hi[k], r_p[k]);
        }
    const double span = std::max({hi[0] - lo[0], hi[1] - lo[1], 0.0});
    const double scale = span > 0.0 ? double((1u << sfc_bits) - 1) / span : 0.0;
    std::vector<std::uint64_t> keys(rPoints.size());
    for (std::size_t i = 0; i < rPoints.size(); ++i) {
        const std::uint32_t q[3] = {static_cast<std::uint32_t>((rPoints[i][0] - lo[0]) * scale),
                                    static_cast<std::uint32_t>((rPoints[i][1] - lo[1]) * scale),
                                    0u};
        keys[i] = sfc_hilbert_key(q, sfc_bits);
    }
    for (const std::int64_t i : sfc_stable_argsort(keys))
        ids[static_cast<std::size_t>(i)] = Insert(rPoints[static_cast<std::size_t>(i)]);
    return ids;
}

void Cdt2d::TriangulatePseudoPolygon(std::int64_t A, std::int64_t B,
                                     const std::vector<std::int64_t>& rChain,
                                     std::vector<Tri>& rOut) const {
    // The chain lies left of A -> B, in order from A to B.
    if (rChain.empty())
        return;
    std::size_t c = 0;
    for (std::size_t i = 1; i < rChain.size(); ++i)
        if (InCircle(A, B, rChain[c], rChain[i]) > 0.0)
            c = i;
    const std::vector<std::int64_t> before(rChain.begin(), rChain.begin() + std::ptrdiff_t(c));
    const std::vector<std::int64_t> after(rChain.begin() + std::ptrdiff_t(c) + 1, rChain.end());
    TriangulatePseudoPolygon(A, rChain[c], before, rOut);
    TriangulatePseudoPolygon(rChain[c], B, after, rOut);
    rOut.push_back({A, B, rChain[c]});
}

bool Cdt2d::IsConstrainedEdge(std::int64_t A, std::int64_t B) const {
    for (const Triangle& r_t : mTris) {
        if (!r_t.mAlive)
            continue;
        for (int e = 0; e < 3; ++e)
            if (r_t.mV[cdt2d_next(e)] == A && r_t.mV[cdt2d_prev(e)] == B)
                return r_t.mC[e];
    }
    return false;
}

Cdt2dStatus Cdt2d::InsertSegment(std::int64_t A, std::int64_t B) {
    while (A != B) {
        // Find the triangle around A whose wedge holds the direction to B.
        const std::int64_t t0 = mVertTri[static_cast<std::size_t>(A)];
        std::int64_t t = t0;
        std::int64_t next_a = -1;
        bool crossing = false;
        int ia = -1;
        for (std::size_t guard = 0; guard < mTris.size() + 3; ++guard) {
            ia = IndexOf(t, A);
            if (ia < 0)
                cdt2d_internal("vertex-triangle map is stale");
            const auto& r_t = mTris[static_cast<std::size_t>(t)];
            const std::int64_t p = r_t.mV[cdt2d_next(ia)];
            const std::int64_t q = r_t.mV[cdt2d_prev(ia)];
            if (p == B || q == B) {
                SetConstrained(t, p == B ? cdt2d_prev(ia) : cdt2d_next(ia));
                return Cdt2dStatus::Ok;
            }
            const double o1 = Orient(A, p, B);
            const double o2 = Orient(A, q, B);
            const Vec2& r_a = Point(A);
            const Vec2& r_b = Point(B);
            const auto ahead = [&](std::int64_t W) {
                const Vec2& r_w = Point(W);
                return (r_w[0] - r_a[0]) * (r_b[0] - r_a[0]) +
                           (r_w[1] - r_a[1]) * (r_b[1] - r_a[1]) >
                       0.0;
            };
            if (o1 == 0.0 && ahead(p)) {
                SetConstrained(t, cdt2d_prev(ia));
                next_a = p;
                break;
            }
            if (o2 == 0.0 && ahead(q)) {
                SetConstrained(t, cdt2d_next(ia));
                next_a = q;
                break;
            }
            if (o1 > 0.0 && o2 < 0.0) {
                crossing = true;
                break;
            }
            // Rotate counter-clockwise around A: across the edge (q, A).
            t = r_t.mN[cdt2d_next(ia)];
            if (t < 0 || t == t0)
                cdt2d_internal("no triangle around a vertex faces the segment");
        }
        if (next_a >= 0) {
            A = next_a;
            continue;
        }
        if (!crossing)
            cdt2d_internal("segment walk did not start");
        // Walk the triangles the segment crosses, collecting both chains.
        std::vector<std::int64_t> crossed{t};
        std::vector<std::int64_t> left;
        std::vector<std::int64_t> right;
        const auto& r_t0 = mTris[static_cast<std::size_t>(t)];
        std::int64_t p = r_t0.mV[cdt2d_next(ia)];  // right of A -> B
        std::int64_t q = r_t0.mV[cdt2d_prev(ia)];  // left of A -> B
        if (r_t0.mC[ia])
            return Cdt2dStatus::SelfIntersecting;
        right.push_back(p);
        left.push_back(q);
        std::int64_t cur = t;
        int edge = ia;
        std::int64_t end = -1;
        for (std::size_t guard = 0; guard < mTris.size() + 3; ++guard) {
            const std::int64_t u = mTris[static_cast<std::size_t>(cur)].mN[edge];
            if (u < 0)
                cdt2d_internal("segment left the triangulation");
            const int f = NeighbourIndex(u, cur);
            const std::int64_t w = mTris[static_cast<std::size_t>(u)].mV[f];
            crossed.push_back(u);
            if (w == B) {
                end = B;
                break;
            }
            const double ow = Orient(A, B, w);
            if (ow == 0.0) {
                end = w;
                break;
            }
            // u = (w, q, p) counter-clockwise: the segment leaves through
            // (p, w) when w is left of A -> B, else through (w, q).
            const auto& r_u = mTris[static_cast<std::size_t>(u)];
            int next_edge = -1;
            if (ow > 0.0) {
                left.push_back(w);
                // The edge (p, w) is opposite q.
                next_edge = IndexOf(u, q);
                q = w;
            } else {
                right.push_back(w);
                next_edge = IndexOf(u, p);
                p = w;
            }
            if (next_edge < 0)
                cdt2d_internal("segment walk lost its edge");
            if (r_u.mC[next_edge])
                return Cdt2dStatus::SelfIntersecting;
            cur = u;
            edge = next_edge;
        }
        if (end < 0)
            cdt2d_internal("segment walk did not end");
        // Re-triangulate both sides of (A, end).
        std::vector<Tri> fresh;
        TriangulatePseudoPolygon(A, end, left, fresh);
        std::vector<std::int64_t> right_rev(right.rbegin(), right.rend());
        TriangulatePseudoPolygon(end, A, right_rev, fresh);
        const bool inside = mClassified && mTris[static_cast<std::size_t>(t)].mInside;
        const std::vector<std::int64_t> ids = ReplaceCavity(crossed, fresh, inside);
        bool marked = false;
        for (const std::int64_t id : ids) {
            for (int e = 0; e < 3 && !marked; ++e) {
                const auto& r_n = mTris[static_cast<std::size_t>(id)];
                if (r_n.mV[cdt2d_next(e)] == A && r_n.mV[cdt2d_prev(e)] == end) {
                    SetConstrained(id, e);
                    marked = true;
                }
            }
        }
        if (!marked)
            cdt2d_internal("recovered segment is not an edge");
        A = end;
    }
    return Cdt2dStatus::Ok;
}

Cdt2dStatus Cdt2d::Classify() {
    std::vector<int> parity(mTris.size(), -1);
    std::vector<std::int64_t> queue;
    for (std::size_t t = 0; t < mTris.size(); ++t) {
        const auto& r_t = mTris[t];
        if (!r_t.mAlive)
            continue;
        if (r_t.mV[0] < kFirstVertex || r_t.mV[1] < kFirstVertex || r_t.mV[2] < kFirstVertex) {
            parity[t] = 0;
            queue.push_back(static_cast<std::int64_t>(t));
        }
    }
    Cdt2dStatus status = Cdt2dStatus::Ok;
    for (std::size_t head = 0; head < queue.size(); ++head) {
        const std::int64_t t = queue[head];
        const auto& r_t = mTris[static_cast<std::size_t>(t)];
        for (int e = 0; e < 3; ++e) {
            const std::int64_t u = r_t.mN[e];
            if (u < 0)
                continue;
            const int want = parity[static_cast<std::size_t>(t)] ^ (r_t.mC[e] ? 1 : 0);
            int& r_p = parity[static_cast<std::size_t>(u)];
            if (r_p < 0) {
                r_p = want;
                queue.push_back(u);
            } else if (r_p != want) {
                status = Cdt2dStatus::OpenLoop;
            }
        }
    }
    for (std::size_t t = 0; t < mTris.size(); ++t)
        mTris[t].mInside = mTris[t].mAlive && parity[t] == 1;
    mClassified = true;
    return status;
}

Cdt2dWhere Cdt2d::Where(const Vec2& rP) {
    if (!mClassified)
        cdt2d_internal("Where before Classify");
    const Location loc = Locate(rP);
    Cdt2dWhere w;
    if (loc.mKind == 2) {
        w.mPlace = Cdt2dPlace::OnVertex;
        w.mVertex = loc.mVertex;
        return w;
    }
    const auto& r_t = mTris[static_cast<std::size_t>(loc.mTri)];
    if (loc.mKind == 1 && r_t.mC[loc.mIndex])
        w.mPlace = Cdt2dPlace::OnConstraint;
    else
        w.mPlace = r_t.mInside ? Cdt2dPlace::Inside : Cdt2dPlace::Outside;
    return w;
}

std::int64_t Cdt2d::InsertSteiner(const Vec2& rP) {
    if (!mClassified)
        cdt2d_internal("InsertSteiner before Classify");
    if (!cdt2d_in_range(rP[0]) || !cdt2d_in_range(rP[1]))
        return -1;
    const Location loc = Locate(rP);
    if (loc.mKind == 2)
        return -1;
    const auto& r_t = mTris[static_cast<std::size_t>(loc.mTri)];
    if (!r_t.mInside)
        return -1;
    if (loc.mKind == 1 && r_t.mC[loc.mIndex])
        return -1;
    return InsertAt(rP, loc, true);
}

std::vector<std::int64_t> Cdt2d::TrianglesAround(std::int64_t V) const {
    std::vector<std::int64_t> out;
    out.reserve(8);
    const std::int64_t t0 = mVertTri[static_cast<std::size_t>(V)];
    std::int64_t t = t0;
    for (std::size_t guard = 0; guard < mTris.size() + 3; ++guard) {
        out.push_back(t);
        const int i = IndexOf(t, V);
        if (i < 0)
            cdt2d_internal("vertex-triangle map is stale");
        // Counter-clockwise: across the edge from V to the previous vertex.
        t = mTris[static_cast<std::size_t>(t)].mN[cdt2d_next(i)];
        if (t < 0)
            cdt2d_internal("a vertex's star is not closed");
        if (t == t0)
            return out;
    }
    cdt2d_internal("a vertex's star does not close");
}

std::vector<std::int64_t> Cdt2d::Neighbours(std::int64_t V) const {
    std::vector<std::int64_t> out;
    out.reserve(8);
    for (const std::int64_t t : TrianglesAround(V))
        out.push_back(mTris[static_cast<std::size_t>(t)].mV[cdt2d_next(IndexOf(t, V))]);
    return out;
}

std::vector<Cdt2d::Tri> Cdt2d::Star(std::int64_t V) const {
    std::vector<Tri> out;
    out.reserve(8);
    for (const std::int64_t t : TrianglesAround(V))
        out.push_back(mTris[static_cast<std::size_t>(t)].mV);
    return out;
}

std::vector<Cdt2d::Tri> Cdt2d::InsideStar(std::int64_t V) const {
    std::vector<Tri> out;
    out.reserve(8);
    for (const std::int64_t t : TrianglesAround(V))
        if (mTris[static_cast<std::size_t>(t)].mInside)
            out.push_back(mTris[static_cast<std::size_t>(t)].mV);
    return out;
}

bool Cdt2d::FindEdge(std::int64_t A, std::int64_t B, std::int64_t& rT, int& rE) const {
    for (const std::int64_t t : TrianglesAround(A)) {
        const int i = IndexOf(t, A);
        if (mTris[static_cast<std::size_t>(t)].mV[cdt2d_next(i)] == B) {
            rT = t;
            rE = cdt2d_prev(i);
            return true;
        }
    }
    return false;
}

bool Cdt2d::FlipEdge(std::int64_t A, std::int64_t B) {
    std::int64_t t = -1;
    int e = -1;
    if (!FindEdge(A, B, t, e))
        return false;
    const auto& r_t = mTris[static_cast<std::size_t>(t)];
    const std::int64_t u = r_t.mN[e];
    if (r_t.mC[e] || u < 0)
        return false;
    const std::int64_t a = r_t.mV[e];
    const std::int64_t d = mTris[static_cast<std::size_t>(u)].mV[NeighbourIndex(u, t)];
    // t = (a, A, B), u = (d, B, A): the new edge (a, d) must cross (A, B).
    if (!(Orient(a, A, d) > 0.0) || !(Orient(a, d, B) > 0.0))
        return false;
    Flip(t, e);
    return true;
}

bool Cdt2d::OnConstraint(std::int64_t V) const {
    for (const std::int64_t t : TrianglesAround(V)) {
        const int i = IndexOf(t, V);
        const auto& r_t = mTris[static_cast<std::size_t>(t)];
        if (r_t.mC[cdt2d_next(i)] || r_t.mC[cdt2d_prev(i)])
            return true;
    }
    return false;
}

std::vector<std::array<std::int64_t, 2>> Cdt2d::ConstrainedEdges() const {
    std::vector<std::array<std::int64_t, 2>> out;
    for (const Triangle& r_t : mTris) {
        if (!r_t.mAlive)
            continue;
        for (int e = 0; e < 3; ++e) {
            const std::int64_t a = r_t.mV[cdt2d_next(e)];
            const std::int64_t b = r_t.mV[cdt2d_prev(e)];
            if (r_t.mC[e] && a < b)
                out.push_back({a, b});
        }
    }
    std::sort(out.begin(), out.end());
    return out;
}

std::int64_t Cdt2d::SplitConstrained(std::int64_t T, int E, const Vec2& rP) {
    // The point is the subsegment's split, computed in floating point, so it
    // may lie an ulp off the line: it splits the edge regardless, provided
    // the four new triangles turn counter-clockwise (exact orient2d).
    if (!cdt2d_in_range(rP[0]) || !cdt2d_in_range(rP[1]))
        return -1;
    const auto& r_t = mTris[static_cast<std::size_t>(T)];
    const std::int64_t o = r_t.mV[E];
    const std::int64_t a = r_t.mV[cdt2d_next(E)];
    const std::int64_t b = r_t.mV[cdt2d_prev(E)];
    const std::int64_t u = r_t.mN[E];
    if (u < 0 || !(Orient(o, a, rP) > 0.0) || !(Orient(b, o, rP) > 0.0))
        return -1;
    const std::int64_t w = mTris[static_cast<std::size_t>(u)].mV[NeighbourIndex(u, T)];
    if (!(Orient(w, b, rP) > 0.0) || !(Orient(a, w, rP) > 0.0))
        return -1;
    Location loc;
    loc.mTri = T;
    loc.mKind = 1;
    loc.mIndex = E;
    return InsertAt(rP, loc, r_t.mInside);
}

bool Cdt2d::MoveVertex(std::int64_t V, const Vec2& rP, bool Legalize) {
    if (!mClassified)
        cdt2d_internal("MoveVertex before Classify");
    if (V < kFirstVertex || V >= NumVertices() || OnConstraint(V))
        return false;
    if (!cdt2d_in_range(rP[0]) || !cdt2d_in_range(rP[1]))
        return false;
    const std::vector<std::int64_t> star = TrianglesAround(V);
    for (const std::int64_t t : star) {
        const auto& r_t = mTris[static_cast<std::size_t>(t)];
        const int i = IndexOf(t, V);
        // The new place must see every edge of the link strictly inside.
        if (!(Orient(r_t.mV[cdt2d_next(i)], r_t.mV[cdt2d_prev(i)], rP) > 0.0))
            return false;
    }
    mPoints[static_cast<std::size_t>(V)] = rP;
    if (!Legalize)
        return true;
    // Lawson's flips over every edge the move touched, each flip's two
    // triangles checked again on all their edges.
    std::vector<std::pair<std::int64_t, int>> stack;
    for (const std::int64_t t : star)
        for (int e = 0; e < 3; ++e)
            stack.emplace_back(t, e);
    std::size_t budget = 64 * mTris.size() + 64;
    while (!stack.empty()) {
        if (budget-- == 0)
            cdt2d_internal("flips after a move did not end");
        const auto [t, e] = stack.back();
        stack.pop_back();
        const auto& r_t = mTris[static_cast<std::size_t>(t)];
        if (!r_t.mAlive || r_t.mC[e])
            continue;
        const std::int64_t u = r_t.mN[e];
        if (u < 0)
            continue;
        const int f = NeighbourIndex(u, t);
        if (f < 0)
            cdt2d_internal("broken adjacency while legalising");
        const std::int64_t d = mTris[static_cast<std::size_t>(u)].mV[f];
        if (InCircle(r_t.mV[0], r_t.mV[1], r_t.mV[2], d) > 0.0) {
            Flip(t, e);
            for (int k = 0; k < 3; ++k) {
                stack.emplace_back(t, k);
                stack.emplace_back(u, k);
            }
        }
    }
    mLast = mVertTri[static_cast<std::size_t>(V)];
    return true;
}

namespace {

/// A bad triangle in the refinement queue.
struct Cdt2dBad {
    double mRatio = 0.0;
    Cdt2d::Tri mKey{};  ///< its vertex ids, ascending
    std::int64_t mSlot = -1;
};

/// The worst ratio first; among equal ratios the smallest vertex ids, so
/// the order depends on the input alone.
struct Cdt2dBadOrder {
    bool operator()(const Cdt2dBad& rA, const Cdt2dBad& rB) const {
        if (rA.mRatio != rB.mRatio)
            return rA.mRatio < rB.mRatio;
        return rA.mKey > rB.mKey;
    }
};

std::uint64_t cdt2d_edge_key(std::int64_t A, std::int64_t B) {
    return (static_cast<std::uint64_t>(std::min(A, B)) << 32) ^
           static_cast<std::uint64_t>(std::max(A, B));
}

Cdt2d::Tri cdt2d_sorted(Cdt2d::Tri T) {
    std::sort(T.begin(), T.end());
    return T;
}

double cdt2d_length(const Cdt2d::Vec2& rA, const Cdt2d::Vec2& rB) {
    return std::hypot(rB[0] - rA[0], rB[1] - rA[1]);
}

}  // namespace

/// Ruppert's refinement over a classified Cdt2d (see Cdt2d::Refine).
class Cdt2dRefiner {
public:
    Cdt2dRefiner(Cdt2d& rC, const Cdt2dRefineOptions& rO) : mC(rC), mO(rO) {}

    Cdt2dRefineReport Run();

private:
    using Vec2 = Cdt2d::Vec2;
    using Tri = Cdt2d::Tri;

    Cdt2d& mC;
    const Cdt2dRefineOptions& mO;
    double mBound = 0.0;
    std::int64_t mBudget = 0;
    Cdt2dRefineReport mRep;
    /// The input segments (the constrained edges at the start), the input
    /// segment of every subsegment, and of every vertex a split made.
    std::vector<std::array<std::int64_t, 2>> mInput;
    std::unordered_map<std::uint64_t, std::int64_t> mSegOf;
    std::unordered_map<std::int64_t, std::int64_t> mVertexSeg;
    /// Small input angles: the apex's smallest wedge, and per input segment
    /// whether each end is the apex of a small wedge it bounds.
    std::map<std::int64_t, double> mApex;
    std::vector<std::array<bool, 2>> mCluster;
    std::priority_queue<Cdt2dBad, std::vector<Cdt2dBad>, Cdt2dBadOrder> mBad;
    std::set<std::array<std::int64_t, 2>> mEncroached;
    std::set<Cdt2dSplitRequest> mRequests;
    std::map<Tri, std::int64_t> mDeclined;  ///< left by the terminator, with the apex
    /// Corner lopping (Pav and Walkington 2005): the balls' radii by corner,
    /// the vertices last found inside each, the balls to shrink, and what
    /// is reported of each.
    std::map<std::int64_t, double> mBall;
    std::map<std::int64_t, std::vector<std::int64_t>> mInsideBall;
    std::set<std::int64_t> mShrink;
    std::map<std::int64_t, Cdt2dLoppedCorner> mLopInfo;
    /// The bad triangles left because their point lay in a ball, queued
    /// again when it shrinks; the balls too small to shrink further; the
    /// balls shrunk for an asked split this pass (once is enough, the
    /// caller makes it); the corners by plane abscissa, with the largest
    /// radius, to find the balls near a point (plane distances only).
    std::map<std::int64_t, std::vector<Cdt2dBad>> mLopDeclined;
    std::set<std::int64_t> mFloored;
    std::set<std::int64_t> mShrunkByAsk;
    std::vector<std::pair<double, std::int64_t>> mBallX;
    double mMaxBall = 0.0;

    const Vec2& P(std::int64_t V) const { return mC.Point(V); }
    bool Real(const Tri& rT) const {
        return rT[0] >= Cdt2d::kFirstVertex && rT[1] >= Cdt2d::kFirstVertex &&
               rT[2] >= Cdt2d::kFirstVertex;
    }
    double PlaneRatio(const Tri& rT) const;
    double Ratio(const Tri& rT) const { return mO.mRatio ? mO.mRatio(rT) : PlaneRatio(rT); }
    double Shortest(const Tri& rT) const;
    /// Ruppert's encroachment: strictly inside the diametral circle. (A
    /// diametral lens admits more circumcentres, but on the suite's
    /// domains it let refinement past 30 degrees run away where the circle
    /// does not.)
    bool Encroaches(std::int64_t A, std::int64_t B, const Vec2& rP) const {
        // With a metric, the diametral circle measured in it at the
        // subsegment's midpoint: a subsegment long in the plane but short
        // in the metric is not encroached by the points its own size puts
        // beside it.
        std::array<Vec2, 4> m{P(A), P(B), rP, Vec2{}};
        if (mC.MetricMapped({0.5 * (P(A)[0] + P(B)[0]), 0.5 * (P(A)[1] + P(B)[1])}, m, 3))
            return diametral(m[0].data(), m[1].data(), m[2].data()) < 0.0;
        return diametral(P(A).data(), P(B).data(), rP.data()) < 0.0;
    }
    std::int64_t SegmentOf(std::int64_t A, std::int64_t B) const {
        const auto it = mSegOf.find(cdt2d_edge_key(A, B));
        return it == mSegOf.end() ? -1 : it->second;
    }
    /// The end of a subsegment at a small input angle of its segment, or -1.
    std::int64_t ApexOf(std::int64_t A, std::int64_t B) const;
    /// The small-angle apexes a vertex lies on a cluster segment of.
    std::vector<std::int64_t> ClustersOf(std::int64_t V) const;
    /// Corner lopping. The distance from input point `C` to vertex `V` (or,
    /// `V` < 0, to `rX`), as the caller measures it.
    double Dist(std::int64_t C, std::int64_t V, const Vec2& rX) const {
        if (mO.mDistance)
            return mO.mDistance(C, V, rX);
        return cdt2d_length(P(C), V >= 0 ? P(V) : rX);
    }
    /// The corners whose ball may hold `rX`, ascending.
    std::vector<std::int64_t> BallsNear(const Vec2& rX) const;
    /// The lowest corner whose ball holds `rX` strictly, or -1.
    std::int64_t BallOf(const Vec2& rX) const;
    /// Whether vertex `V` lies inside corner `C`'s ball (not on it, where
    /// the ball's own splits are).
    bool InsideBall(std::int64_t C, std::int64_t V) const;
    /// The constrained neighbours of `V`, counter-clockwise, and the inside
    /// wedges between consecutive ones.
    std::vector<std::int64_t> ConstrainedNeighbours(std::int64_t V) const;
    std::vector<std::array<std::int64_t, 2>> InsideWedges(std::int64_t V) const;
    void SetUpBalls();
    /// SPLITBALL(C, 1): split (or ask to split) every subsegment at `C`
    /// where it leaves the ball; returns whether one was split or asked.
    bool SplitIncident(std::int64_t C);
    /// SPLITBALL(C, beta): shrink the ball, split where it now ends, and
    /// queue it again while a vertex is still inside.
    bool ShrinkBall(std::int64_t C);
    /// Queue the balls a new vertex lies in.
    void CheckBalls(std::int64_t V);
    /// Split constrained edge `A`-`B` at `rM` (book-keeping included).
    bool SplitAt(std::int64_t A, std::int64_t B, const Vec2& rM);
    /// The point on subsegment `C`-`W` at distance `D` from `C`.
    Vec2 PointAtDistance(std::int64_t C, std::int64_t W, double D) const;
    /// Ask for the split of `A`-`B` at distance `D` from corner `C`.
    void AskAtBall(std::int64_t A, std::int64_t B, std::int64_t C, double D);
    /// The small-angle apex that forces a triangle (it holds the angle, or
    /// its shortest edge spans it), or -1.
    std::int64_t ForcedBy(const Tri& rT) const;
    bool Forced(const Tri& rT) const { return ForcedBy(rT) >= 0; }
    double MinAngleDeg(const Tri& rT) const;
    void FindSmallAngles();
    void PushTriangle(std::int64_t Slot);
    void CheckSubsegment(std::int64_t A, std::int64_t B);
    void AfterInsert(std::int64_t V);
    bool SplitOrAsk(std::int64_t A, std::int64_t B);
    void HandleEncroached(const std::vector<std::array<std::int64_t, 2>>& rSegments,
                          const Tri& rBad, std::int64_t Slot);
    /// The point to insert for a bad triangle (its circumcentre, or the
    /// caller's).
    Vec2 CentreOf(const Tri& rT) const;
    void Attempt(const Cdt2dBad& rBad);
};

double Cdt2dRefiner::PlaneRatio(const Tri& rT) const {
    const double la = cdt2d_length(P(rT[1]), P(rT[2]));
    const double lb = cdt2d_length(P(rT[2]), P(rT[0]));
    const double lc = cdt2d_length(P(rT[0]), P(rT[1]));
    const double twice_area =
        std::fabs(orient2d(P(rT[0]).data(), P(rT[1]).data(), P(rT[2]).data()));
    if (!(twice_area > 0.0))
        return std::numeric_limits<double>::infinity();
    // R = la lb lc / (4 area): the ratio R / shortest is the product of the
    // two longer edges over twice the doubled area.
    return la * lb * lc / (2.0 * twice_area * std::min({la, lb, lc}));
}

double Cdt2dRefiner::Shortest(const Tri& rT) const {
    return std::min({cdt2d_length(P(rT[1]), P(rT[2])), cdt2d_length(P(rT[2]), P(rT[0])),
                     cdt2d_length(P(rT[0]), P(rT[1]))});
}

std::int64_t Cdt2dRefiner::ApexOf(std::int64_t A, std::int64_t B) const {
    const std::int64_t s = SegmentOf(A, B);
    if (s < 0)
        return -1;
    const auto& r_in = mInput[static_cast<std::size_t>(s)];
    const bool at0 = mCluster[static_cast<std::size_t>(s)][0] && (A == r_in[0] || B == r_in[0]);
    const bool at1 = mCluster[static_cast<std::size_t>(s)][1] && (A == r_in[1] || B == r_in[1]);
    if (at0 == at1)
        return -1;  // neither, or the whole segment between two apexes: its midpoint
    return at0 ? r_in[0] : r_in[1];
}

std::vector<std::int64_t> Cdt2dRefiner::ClustersOf(std::int64_t V) const {
    std::vector<std::int64_t> out;
    const auto it = mVertexSeg.find(V);
    if (it == mVertexSeg.end())
        return out;
    const auto s = static_cast<std::size_t>(it->second);
    for (int k = 0; k < 2; ++k)
        if (mCluster[s][static_cast<std::size_t>(k)])
            out.push_back(mInput[s][static_cast<std::size_t>(k)]);
    return out;
}

std::vector<std::int64_t> Cdt2dRefiner::ConstrainedNeighbours(std::int64_t V) const {
    std::vector<std::int64_t> out;
    for (const std::int64_t t : mC.TrianglesAround(V)) {
        const auto& r_t = mC.mTris[static_cast<std::size_t>(t)];
        const int i = mC.IndexOf(t, V);
        // The edge (V, p) opens this triangle, p its next vertex.
        if (r_t.mC[cdt2d_prev(i)])
            out.push_back(r_t.mV[cdt2d_next(i)]);
    }
    return out;
}

std::vector<std::array<std::int64_t, 2>> Cdt2dRefiner::InsideWedges(std::int64_t V) const {
    // The constrained neighbours of V, counter-clockwise, and whether the
    // wedge after each (to the next) is inside the domain.
    std::vector<std::int64_t> nbrs;
    std::vector<bool> wedge_inside;
    for (const std::int64_t t : mC.TrianglesAround(V)) {
        const auto& r_t = mC.mTris[static_cast<std::size_t>(t)];
        const int i = mC.IndexOf(t, V);
        if (r_t.mC[cdt2d_prev(i)]) {
            nbrs.push_back(r_t.mV[cdt2d_next(i)]);
            wedge_inside.push_back(r_t.mInside);
        }
    }
    std::vector<std::array<std::int64_t, 2>> out;
    const std::size_t n = nbrs.size();
    for (std::size_t k = 0; k < n && n >= 2; ++k)
        if (wedge_inside[k])
            out.push_back({nbrs[k], nbrs[(k + 1) % n]});
    return out;
}

// A vertex closer to a corner than this fraction of its ball's radius under
// it is inside; the ball's own splits, at the radius, are not.
constexpr double kCdt2dLopTol = 1e-9;
// A ball is not shrunk below this fraction of its corner's lfs0 (a vertex
// that close is one rounding put there).
constexpr double kCdt2dLopFloor = 1e-12;

std::vector<std::int64_t> Cdt2dRefiner::BallsNear(const Vec2& rX) const {
    std::vector<std::int64_t> out;
    if (mO.mDistance) {
        for (const auto& r_b : mBall)
            out.push_back(r_b.first);
        return out;
    }
    // Balls only shrink, so the largest radius at the start bounds them all.
    auto it = std::lower_bound(mBallX.begin(), mBallX.end(),
                               std::make_pair(rX[0] - mMaxBall, std::int64_t{-1}));
    for (; it != mBallX.end() && it->first <= rX[0] + mMaxBall; ++it)
        out.push_back(it->second);
    std::sort(out.begin(), out.end());
    return out;
}

std::int64_t Cdt2dRefiner::BallOf(const Vec2& rX) const {
    for (const std::int64_t c : BallsNear(rX))
        if (Dist(c, -1, rX) < mBall.at(c))
            return c;
    return -1;
}

bool Cdt2dRefiner::InsideBall(std::int64_t C, std::int64_t V) const {
    return V != C && Dist(C, V, P(V)) < mBall.at(C) * (1.0 - kCdt2dLopTol);
}

void Cdt2dRefiner::CheckBalls(std::int64_t V) {
    if (mBall.empty())
        return;
    for (const std::int64_t c : BallsNear(P(V)))
        if (mFloored.count(c) == 0 && InsideBall(c, V)) {
            mInsideBall[c].push_back(V);
            mShrink.insert(c);
        }
}

Cdt2dRefiner::Vec2 Cdt2dRefiner::PointAtDistance(std::int64_t C, std::int64_t W, double D) const {
    const Vec2& a = P(C);
    const Vec2& b = P(W);
    const auto at = [&](double F) {
        return Vec2{a[0] + F * (b[0] - a[0]), a[1] + F * (b[1] - a[1])};
    };
    if (!mO.mDistance)
        return at(D / cdt2d_length(a, b));
    // The caller's distance grows along the subsegment: bisect for it.
    double lo = 0.0;
    double hi = 1.0;
    for (int k = 0; k < 64 && hi - lo > 1e-15; ++k) {
        const double mid = 0.5 * (lo + hi);
        (Dist(C, -1, at(mid)) < D ? lo : hi) = mid;
    }
    return at(0.5 * (lo + hi));
}

void Cdt2dRefiner::AskAtBall(std::int64_t A, std::int64_t B, std::int64_t C, double D) {
    Cdt2dSplitRequest r{std::min(A, B), std::max(A, B), C, D};
    mRequests.erase(r);  // the ball's last radius replaces an earlier one
    mRequests.insert(r);
}

bool Cdt2dRefiner::SplitIncident(std::int64_t C) {
    const double d = mBall.at(C);
    bool acted = false;
    for (const std::int64_t w : ConstrainedNeighbours(C)) {
        // A neighbour inside (or on) the ball is the shrink rule's.
        if (!(Dist(C, w, P(w)) > d * (1.0 + kCdt2dLopTol)))
            continue;
        if (!mO.mSplitSegments) {
            AskAtBall(C, w, C, d);
            acted = true;
            continue;
        }
        if (mBudget <= 0)
            break;
        acted = SplitAt(C, w, PointAtDistance(C, w, d)) || acted;
    }
    return acted;
}

bool Cdt2dRefiner::ShrinkBall(std::int64_t C) {
    auto& r_d = mBall.at(C);
    auto& r_info = mLopInfo.at(C);
    if (!(mO.mLopBeta * r_d > kCdt2dLopFloor * r_info.mLfs0)) {
        mFloored.insert(C);
        mInsideBall.erase(C);
        return false;
    }
    r_d *= mO.mLopBeta;
    r_info.mRadius = r_d;
    // The triangles left for this ball may now be refined.
    const auto it = mLopDeclined.find(C);
    if (it != mLopDeclined.end()) {
        for (const Cdt2dBad& r_b : it->second) {
            const auto& r_t = mC.mTris[static_cast<std::size_t>(r_b.mSlot)];
            if (r_t.mAlive && r_t.mInside && cdt2d_sorted(r_t.mV) == r_b.mKey) {
                mDeclined.erase(r_b.mKey);
                mBad.push(r_b);
            }
        }
        mLopDeclined.erase(it);
    }
    const bool acted = SplitIncident(C);
    auto& r_in = mInsideBall[C];
    r_in.erase(
        std::remove_if(r_in.begin(), r_in.end(), [&](std::int64_t V) { return !InsideBall(C, V); }),
        r_in.end());
    if (!r_in.empty())
        mShrink.insert(C);
    return acted;
}

void Cdt2dRefiner::SetUpBalls() {
    if (!(mO.mLopBeta >= 0.5 && mO.mLopBeta < 1.0))
        throw std::invalid_argument("meshio++: cdt2d: lop_beta must lie in [1/2, 1)");
    // The acute inside wedge at a vertex (degrees; 180 for none): turning
    // left at under 90 degrees.
    const auto acute_deg = [&](std::int64_t V) {
        double deg = 180.0;
        for (const auto& r_w : InsideWedges(V)) {
            const double* pv = P(V).data();
            if (!(orient2d(pv, P(r_w[0]).data(), P(r_w[1]).data()) > 0.0) ||
                !(diametral(P(r_w[0]).data(), P(r_w[1]).data(), pv) > 0.0))
                continue;
            const double ax = P(r_w[0])[0] - P(V)[0], ay = P(r_w[0])[1] - P(V)[1];
            const double bx = P(r_w[1])[0] - P(V)[0], by = P(r_w[1])[1] - P(V)[1];
            deg = std::min(deg, std::atan2(ax * by - ay * bx, ax * bx + ay * by) *
                                    (180.0 / 3.14159265358979323846));
        }
        return deg;
    };
    // The input points P: the caller's, and every end of an input segment at
    // an acute corner (all of them when the caller gives none).
    std::vector<std::int64_t> pts = mO.mCorners;
    for (const std::int64_t p : pts)
        if (p < Cdt2d::kFirstVertex || p >= mC.NumVertices())
            throw std::invalid_argument("meshio++: cdt2d: a lopped corner is not a vertex");
    const bool all_ends = pts.empty();
    std::vector<std::int64_t> ends;
    for (const auto& r_s : mInput)
        ends.insert(ends.end(), r_s.begin(), r_s.end());
    std::sort(ends.begin(), ends.end());
    ends.erase(std::unique(ends.begin(), ends.end()), ends.end());
    for (const std::int64_t v : ends)
        if (all_ends || acute_deg(v) < 90.0)
            pts.push_back(v);
    std::sort(pts.begin(), pts.end());
    pts.erase(std::unique(pts.begin(), pts.end()), pts.end());
    if (!mO.mLopSkip.empty()) {
        const std::set<std::int64_t> skip(mO.mLopSkip.begin(), mO.mLopSkip.end());
        pts.erase(std::remove_if(pts.begin(), pts.end(),
                                 [&](std::int64_t V) { return skip.count(V) != 0; }),
                  pts.end());
    }
    std::map<std::int64_t, double> carried;
    for (const auto& [v, d] : mO.mLopRadius)
        if (d > 0.0 && std::isfinite(d)) {
            auto [it, fresh] = carried.emplace(v, d);
            if (!fresh)
                it->second = std::min(it->second, d);
        }
    // The plane abscissae, to find the nearest input point (lfs0) with the
    // plane's distance; the caller's is searched in full.
    std::vector<std::pair<double, std::int64_t>> by_x;
    for (const std::int64_t p : pts)
        by_x.emplace_back(P(p)[0], p);
    std::sort(by_x.begin(), by_x.end());
    const auto lfs0 = [&](std::int64_t Q) {
        double best = std::numeric_limits<double>::infinity();
        if (mO.mDistance) {
            for (const std::int64_t q : pts)
                if (q != Q)
                    best = std::min(best, Dist(Q, q, P(q)));
            return best;
        }
        const auto at = std::lower_bound(by_x.begin(), by_x.end(), std::make_pair(P(Q)[0], Q));
        for (auto it = at; it != by_x.end() && it->first - P(Q)[0] < best; ++it)
            if (it->second != Q)
                best = std::min(best, Dist(Q, it->second, P(it->second)));
        for (auto it = at; it != by_x.begin();) {
            --it;
            if (P(Q)[0] - it->first >= best)
                break;
            best = std::min(best, Dist(Q, it->second, P(it->second)));
        }
        return best;
    };
    for (const std::int64_t p : pts) {
        const double deg = acute_deg(p);
        if (!(deg < 90.0))
            continue;
        const double l = lfs0(p);
        if (!(l > 0.0) || !std::isfinite(l))
            continue;
        double d = (std::sqrt(2.0) - 1.0) * l;
        const auto it = carried.find(p);
        if (it != carried.end())
            d = std::min(d, it->second);
        mBall[p] = d;
        mLopInfo[p] = {p, deg, l, d, 0, 180.0, true};
        mBallX.emplace_back(P(p)[0], p);
        mMaxBall = std::max(mMaxBall, d);
    }
    std::sort(mBallX.begin(), mBallX.end());
    if (mBall.empty())
        return;
    // The vertices already inside a ball (those of the inside triangles),
    // then the balls' splits.
    std::vector<bool> seen(static_cast<std::size_t>(mC.NumVertices()), false);
    for (const auto& r_t : mC.mTris) {
        if (!r_t.mAlive || !r_t.mInside || !Real(r_t.mV))
            continue;
        for (const std::int64_t v : r_t.mV)
            if (!seen[static_cast<std::size_t>(v)]) {
                seen[static_cast<std::size_t>(v)] = true;
                CheckBalls(v);
            }
    }
    for (const auto& r_b : mBall)
        SplitIncident(r_b.first);
}

double Cdt2dRefiner::MinAngleDeg(const Tri& rT) const {
    if (mO.mMinAngle)
        return mO.mMinAngle(rT);
    double m = 180.0;
    for (int k = 0; k < 3; ++k) {
        const Vec2& o = P(rT[static_cast<std::size_t>(k)]);
        const Vec2& p = P(rT[static_cast<std::size_t>((k + 1) % 3)]);
        const Vec2& q = P(rT[static_cast<std::size_t>((k + 2) % 3)]);
        const double ux = p[0] - o[0], uy = p[1] - o[1];
        const double vx = q[0] - o[0], vy = q[1] - o[1];
        m = std::min(m, std::atan2(std::fabs(ux * vy - uy * vx), ux * vx + uy * vy) *
                            (180.0 / 3.14159265358979323846));
    }
    return m;
}

std::int64_t Cdt2dRefiner::ForcedBy(const Tri& rT) const {
    if (mApex.empty())
        return -1;
    // The triangle holding a small input angle: its two edges at the apex
    // are subsegments of the cluster.
    for (int k = 0; k < 3; ++k) {
        const std::int64_t v = rT[static_cast<std::size_t>(k)];
        if (mApex.count(v) == 0)
            continue;
        const std::int64_t a = rT[static_cast<std::size_t>((k + 1) % 3)];
        const std::int64_t b = rT[static_cast<std::size_t>((k + 2) % 3)];
        if (ApexOf(v, a) == v && ApexOf(v, b) == v)
            return v;
    }
    // A shortest edge spanning a small input angle: its ends on two
    // segments of one cluster.
    int sk = 0;
    double best = std::numeric_limits<double>::infinity();
    for (int k = 0; k < 3; ++k) {
        const double l = cdt2d_length(P(rT[static_cast<std::size_t>((k + 1) % 3)]),
                                      P(rT[static_cast<std::size_t>((k + 2) % 3)]));
        if (l < best) {
            best = l;
            sk = k;
        }
    }
    const std::int64_t p = rT[static_cast<std::size_t>((sk + 1) % 3)];
    const std::int64_t q = rT[static_cast<std::size_t>((sk + 2) % 3)];
    const auto it_p = mVertexSeg.find(p);
    const auto it_q = mVertexSeg.find(q);
    if (it_p == mVertexSeg.end() || it_q == mVertexSeg.end() || it_p->second == it_q->second)
        return -1;
    for (const std::int64_t a : ClustersOf(p))
        for (const std::int64_t b : ClustersOf(q))
            if (a == b)
                return a;
    return -1;
}

void Cdt2dRefiner::FindSmallAngles() {
    std::map<std::int64_t, std::vector<std::int64_t>> incident;
    for (const auto& r_e : mC.ConstrainedEdges()) {
        const auto s = static_cast<std::int64_t>(mInput.size());
        mInput.push_back(r_e);
        mSegOf[cdt2d_edge_key(r_e[0], r_e[1])] = s;
        incident[r_e[0]].push_back(s);
        incident[r_e[1]].push_back(s);
    }
    mCluster.assign(mInput.size(), {false, false});
    if (mO.mLopCorners)
        return;  // corner lopping protects small angles in place of shells
    for (const auto& [v, r_segs] : incident) {
        if (r_segs.size() < 2)
            continue;
        for (const auto& r_w : InsideWedges(v)) {
            const std::int64_t a = r_w[0];
            const std::int64_t b = r_w[1];
            const double* pv = P(v).data();
            // Below 60 degrees: turning left (under 180), acute, and in the
            // lens_angle region (below 60 or above 120).
            if (!(orient2d(pv, P(a).data(), P(b).data()) > 0.0) ||
                !(diametral(P(a).data(), P(b).data(), pv) > 0.0) ||
                !(lens_angle(P(a).data(), P(b).data(), pv) > 0.0))
                continue;
            const double ax = P(a)[0] - P(v)[0], ay = P(a)[1] - P(v)[1];
            const double bx = P(b)[0] - P(v)[0], by = P(b)[1] - P(v)[1];
            const double deg =
                std::atan2(ax * by - ay * bx, ax * bx + ay * by) * (180.0 / 3.14159265358979323846);
            auto [it, fresh] = mApex.emplace(v, deg);
            if (!fresh)
                it->second = std::min(it->second, deg);
            for (const std::int64_t w : {a, b}) {
                const std::int64_t s = SegmentOf(v, w);
                if (s >= 0)
                    mCluster[static_cast<std::size_t>(s)]
                            [mInput[static_cast<std::size_t>(s)][0] == v ? 0 : 1] = true;
            }
        }
    }
    for (const auto& [v, deg] : mApex)
        mRep.mSmallAngles.push_back({v, deg});
}

void Cdt2dRefiner::PushTriangle(std::int64_t Slot) {
    const auto& r_t = mC.mTris[static_cast<std::size_t>(Slot)];
    if (!r_t.mAlive || !r_t.mInside || !Real(r_t.mV))
        return;
    const double r = Ratio(r_t.mV);
    if (r > mBound)
        mBad.push({r, cdt2d_sorted(r_t.mV), Slot});
}

void Cdt2dRefiner::CheckSubsegment(std::int64_t A, std::int64_t B) {
    std::int64_t t = -1;
    int e = -1;
    if (!mC.FindEdge(A, B, t, e))
        return;
    const auto& r_t = mC.mTris[static_cast<std::size_t>(t)];
    if (!r_t.mC[e])
        return;
    for (const std::int64_t side : {t, r_t.mN[e]}) {
        if (side < 0)
            continue;
        const auto& r_s = mC.mTris[static_cast<std::size_t>(side)];
        if (!r_s.mInside || !Real(r_s.mV))
            continue;
        for (const std::int64_t apex : r_s.mV)
            if (apex != A && apex != B && Encroaches(A, B, P(apex))) {
                mEncroached.insert({std::min(A, B), std::max(A, B)});
                return;
            }
    }
}

void Cdt2dRefiner::AfterInsert(std::int64_t V) {
    for (const std::int64_t t : mC.TrianglesAround(V)) {
        PushTriangle(t);
        const auto& r_t = mC.mTris[static_cast<std::size_t>(t)];
        const int i = mC.IndexOf(t, V);
        const std::int64_t a = r_t.mV[cdt2d_next(i)];
        const std::int64_t b = r_t.mV[cdt2d_prev(i)];
        if (r_t.mC[i])
            CheckSubsegment(a, b);  // the link edge, now facing V
        if (r_t.mC[cdt2d_prev(i)])
            CheckSubsegment(V, a);  // a subsegment at V
    }
    CheckBalls(V);
}

bool Cdt2dRefiner::SplitAt(std::int64_t A, std::int64_t B, const Vec2& rM) {
    std::int64_t t = -1;
    int e = -1;
    if (!mC.FindEdge(A, B, t, e) || !mC.mTris[static_cast<std::size_t>(t)].mC[e])
        return false;
    const std::int64_t s = SegmentOf(A, B);
    const std::int64_t v = mC.SplitConstrained(t, e, rM);
    if (v < 0)
        return false;
    --mBudget;
    ++mRep.mNumSplits;
    if (s >= 0) {
        mSegOf.erase(cdt2d_edge_key(A, B));
        mSegOf[cdt2d_edge_key(A, v)] = s;
        mSegOf[cdt2d_edge_key(v, B)] = s;
        mVertexSeg[v] = s;
    }
    AfterInsert(v);
    return true;
}

bool Cdt2dRefiner::SplitOrAsk(std::int64_t A, std::int64_t B) {
    if (!mBall.empty()) {
        // A subsegment at a lopped corner: one that reaches past the ball is
        // split where it leaves it; the ball's own is lopped again,
        // SPLITBALL(c, beta) (at beta = 1/2, its midpoint).
        bool at_ball = false;
        bool acted = false;
        for (const std::int64_t c : {std::min(A, B), std::max(A, B)}) {
            const auto it = mBall.find(c);
            if (it == mBall.end() || mFloored.count(c) != 0)
                continue;
            at_ball = true;
            const std::int64_t w = c == A ? B : A;
            if (Dist(c, w, P(w)) > it->second * (1.0 + kCdt2dLopTol)) {
                if (!mO.mSplitSegments) {
                    AskAtBall(c, w, c, it->second);
                    acted = true;
                } else {
                    acted = SplitAt(c, w, PointAtDistance(c, w, it->second)) || acted;
                }
            } else if (mO.mSplitSegments || mShrunkByAsk.insert(c).second) {
                acted = ShrinkBall(c) || acted;
            }
        }
        if (at_ball)
            return acted && mO.mSplitSegments;
    }
    const std::int64_t apex = ApexOf(A, B);
    if (!mO.mSplitSegments) {
        mRequests.insert({std::min(A, B), std::max(A, B), apex});
        return false;
    }
    Vec2 m{0.5 * (P(A)[0] + P(B)[0]), 0.5 * (P(A)[1] + P(B)[1])};
    if (apex >= 0) {
        // Concentric shells: a power-of-two distance from the apex in
        // (L / 3, 2 L / 3], so the subsegments of a cluster stay equal.
        const std::int64_t far = apex == A ? B : A;
        const double len = cdt2d_length(P(apex), P(far));
        int ex = 0;
        std::frexp(2.0 * len / 3.0, &ex);
        const double r = std::ldexp(1.0, ex - 1);
        if (r > len / 3.0 && r <= 2.0 * len / 3.0) {
            const double f = r / len;
            m = {P(apex)[0] + f * (P(far)[0] - P(apex)[0]),
                 P(apex)[1] + f * (P(far)[1] - P(apex)[1])};
        }
    }
    return SplitAt(A, B, m);
}

void Cdt2dRefiner::HandleEncroached(const std::vector<std::array<std::int64_t, 2>>& rSegments,
                                    const Tri& rBad, std::int64_t Slot) {
    const double shortest = Shortest(rBad);
    bool acted = false;
    std::int64_t declined_apex = -1;
    for (const auto& r_s : rSegments) {
        // The terminator: a subsegment of a small-angle cluster is not split
        // when its pieces would be shorter than the skinny triangle's
        // shortest edge; the triangle is what the angle forces.
        const std::int64_t s = SegmentOf(r_s[0], r_s[1]);
        const bool cluster = s >= 0 && (mCluster[static_cast<std::size_t>(s)][0] ||
                                        mCluster[static_cast<std::size_t>(s)][1]);
        if (mO.mAcceptSplit && !mO.mAcceptSplit(r_s[0], r_s[1], rBad))
            continue;
        if (cluster && 0.5 * cdt2d_length(P(r_s[0]), P(r_s[1])) < shortest) {
            const auto& r_in = mInput[static_cast<std::size_t>(s)];
            declined_apex = mCluster[static_cast<std::size_t>(s)][0] ? r_in[0] : r_in[1];
            continue;
        }
        // A split the triangulation refuses (its midpoint rounds off the
        // segment) leaves the triangle declined, not queued again: the same
        // attempt would refuse it forever. Asking always acts.
        if (SplitOrAsk(r_s[0], r_s[1]) || !mO.mSplitSegments)
            acted = true;
        if (mBudget <= 0)
            return;
    }
    if (!acted) {
        mDeclined.emplace(cdt2d_sorted(rBad), declined_apex);
        return;
    }
    if (mO.mSplitSegments) {
        const auto& r_t = mC.mTris[static_cast<std::size_t>(Slot)];
        if (r_t.mAlive && cdt2d_sorted(r_t.mV) == cdt2d_sorted(rBad))
            PushTriangle(Slot);
    }
}

Cdt2dRefiner::Vec2 Cdt2dRefiner::CentreOf(const Tri& rT) const {
    if (mO.mCentre)
        return mO.mCentre(rT);
    // The circumcentre, relative to the first vertex.
    const Vec2& a = P(rT[0]);
    const double bx = P(rT[1])[0] - a[0], by = P(rT[1])[1] - a[1];
    const double cx = P(rT[2])[0] - a[0], cy = P(rT[2])[1] - a[1];
    const double d = 2.0 * orient2d(a.data(), P(rT[1]).data(), P(rT[2]).data());
    const double b2 = bx * bx + by * by;
    const double c2 = cx * cx + cy * cy;
    return {a[0] + (cy * b2 - by * c2) / d, a[1] + (bx * c2 - cx * b2) / d};
}

void Cdt2dRefiner::Attempt(const Cdt2dBad& rBad) {
    const auto& r_t = mC.mTris[static_cast<std::size_t>(rBad.mSlot)];
    if (!r_t.mAlive || !r_t.mInside || cdt2d_sorted(r_t.mV) != rBad.mKey)
        return;  // gone since it was queued
    const Tri tri = r_t.mV;
    if (Forced(tri))
        return;
    const Vec2& a = P(tri[0]);
    const Vec2 c = CentreOf(tri);
    if (!std::isfinite(c[0]) || !std::isfinite(c[1]) || !cdt2d_in_range(c[0]) ||
        !cdt2d_in_range(c[1]))
        return;
    // Corner lopping: a point in a ball is left (until the ball shrinks).
    if (!mBall.empty()) {
        const std::int64_t q = BallOf(c);
        if (q >= 0) {
            mDeclined[rBad.mKey] = q;
            mLopDeclined[q].push_back(rBad);
            return;
        }
    }
    // Walk from the centroid toward it; a constrained edge on the way hides
    // it, and is split instead (Chew 1993; Shewchuk 2002).
    const Vec2 g{(a[0] + P(tri[1])[0] + P(tri[2])[0]) / 3.0,
                 (a[1] + P(tri[1])[1] + P(tri[2])[1]) / 3.0};
    std::int64_t cur = rBad.mSlot;
    std::int64_t prev = -1;
    std::int64_t reached = -1;
    for (std::size_t guard = 0; guard < mC.mTris.size() + 3 && reached < 0; ++guard) {
        const auto& r_c = mC.mTris[static_cast<std::size_t>(cur)];
        // Leave through an edge the point lies beyond, preferring the one
        // the line from the centroid crosses, then one not walked back over.
        int exit = -1;
        int forward = -1;
        int any = -1;
        for (int e = 0; e < 3; ++e) {
            const std::int64_t x = r_c.mV[cdt2d_next(e)];
            const std::int64_t y = r_c.mV[cdt2d_prev(e)];
            if (!(mC.Orient(x, y, c) < 0.0))
                continue;
            if (any < 0)
                any = e;
            if (r_c.mN[e] == prev)
                continue;
            if (forward < 0)
                forward = e;
            const double ox = orient2d(g.data(), c.data(), P(x).data());
            const double oy = orient2d(g.data(), c.data(), P(y).data());
            if (exit < 0 && ((ox <= 0.0 && oy >= 0.0) || (ox >= 0.0 && oy <= 0.0)))
                exit = e;
        }
        if (any < 0) {
            reached = cur;
            break;
        }
        if (exit < 0)
            exit = forward >= 0 ? forward : any;
        if (r_c.mC[exit]) {
            HandleEncroached({{r_c.mV[cdt2d_next(exit)], r_c.mV[cdt2d_prev(exit)]}}, tri,
                             rBad.mSlot);
            return;
        }
        prev = cur;
        cur = r_c.mN[exit];
        if (cur < 0)
            return;
    }
    if (reached < 0)
        return;
    const auto& r_r = mC.mTris[static_cast<std::size_t>(reached)];
    if (!r_r.mInside)
        return;
    int zeros = 0;
    int on = -1;
    for (int e = 0; e < 3; ++e)
        if (mC.Orient(r_r.mV[cdt2d_next(e)], r_r.mV[cdt2d_prev(e)], c) == 0.0) {
            ++zeros;
            on = e;
        }
    if (zeros >= 2)
        return;  // on a vertex
    if (zeros == 1 && r_r.mC[on]) {
        HandleEncroached({{r_r.mV[cdt2d_next(on)], r_r.mV[cdt2d_prev(on)]}}, tri, rBad.mSlot);
        return;
    }
    // The cavity the point would open (triangles whose circumcircle holds
    // it, not across a constraint): the subsegments on its boundary are the
    // ones it can encroach.
    std::vector<std::int64_t> cavity{reached};
    std::vector<std::array<std::int64_t, 2>> encroached;
    for (std::size_t k = 0; k < cavity.size(); ++k) {
        const auto& r_k = mC.mTris[static_cast<std::size_t>(cavity[k])];
        for (int e = 0; e < 3; ++e) {
            const std::int64_t x = r_k.mV[cdt2d_next(e)];
            const std::int64_t y = r_k.mV[cdt2d_prev(e)];
            if (r_k.mC[e]) {
                if (Encroaches(x, y, c))
                    encroached.push_back({std::min(x, y), std::max(x, y)});
                continue;
            }
            const std::int64_t n = r_k.mN[e];
            if (n < 0 || std::find(cavity.begin(), cavity.end(), n) != cavity.end())
                continue;
            const auto& r_n = mC.mTris[static_cast<std::size_t>(n)];
            if (incircle(P(r_n.mV[0]).data(), P(r_n.mV[1]).data(), P(r_n.mV[2]).data(), c.data()) >
                0.0)
                cavity.push_back(n);
        }
    }
    if (!encroached.empty()) {
        std::sort(encroached.begin(), encroached.end());
        encroached.erase(std::unique(encroached.begin(), encroached.end()), encroached.end());
        HandleEncroached(encroached, tri, rBad.mSlot);
        return;
    }
    if (mO.mAccept) {
        std::vector<std::int64_t> around;
        for (const std::int64_t t : cavity)
            for (const std::int64_t w : mC.mTris[static_cast<std::size_t>(t)].mV)
                around.push_back(w);
        std::sort(around.begin(), around.end());
        around.erase(std::unique(around.begin(), around.end()), around.end());
        if (!mO.mAccept(tri, c, around))
            return;
    }
    const std::int64_t v = mC.InsertSteiner(c);
    if (v < 0)
        return;
    --mBudget;
    ++mRep.mNumCircumcentres;
    AfterInsert(v);
}

Cdt2dRefineReport Cdt2dRefiner::Run() {
    if (!mC.mClassified)
        cdt2d_internal("Refine before Classify");
    if (!(mO.mMinAngleDeg > 0.0 && mO.mMinAngleDeg < 60.0))
        throw std::invalid_argument("meshio++: cdt2d: min_angle must lie in (0, 60) degrees");
    mBound = 1.0 / (2.0 * std::sin(mO.mMinAngleDeg * (3.14159265358979323846 / 180.0)));
    mBudget = mO.mMaxPoints;
    FindSmallAngles();
    if (mO.mLopCorners)
        SetUpBalls();
    for (const auto& r_s : mInput)
        CheckSubsegment(r_s[0], r_s[1]);
    for (std::size_t t = 0; t < mC.mTris.size(); ++t)
        PushTriangle(static_cast<std::int64_t>(t));
    while (true) {
        if (mBudget <= 0) {
            mRep.mBudgetReached = true;
            break;
        }
        // Corner lopping's rules in their priority: a vertex inside a ball,
        // an encroached subsegment, a bad triangle.
        if (!mShrink.empty()) {
            const std::int64_t c = *mShrink.begin();
            mShrink.erase(mShrink.begin());
            auto& r_in = mInsideBall[c];
            r_in.erase(std::remove_if(r_in.begin(), r_in.end(),
                                      [&](std::int64_t V) { return !InsideBall(c, V); }),
                       r_in.end());
            if (!r_in.empty() && mFloored.count(c) == 0)
                ShrinkBall(c);
            continue;
        }
        if (!mEncroached.empty()) {
            const auto s = *mEncroached.begin();
            mEncroached.erase(mEncroached.begin());
            // Still a subsegment, and still encroached by an apex?
            const std::size_t before = mEncroached.size();
            CheckSubsegment(s[0], s[1]);
            if (mEncroached.size() > before) {
                mEncroached.erase(s);
                SplitOrAsk(s[0], s[1]);
            }
            continue;
        }
        if (mBad.empty())
            break;
        const Cdt2dBad bad = mBad.top();
        mBad.pop();
        Attempt(bad);
    }
    std::map<std::int64_t, double> forced_min;
    for (const auto& r_t : mC.mTris) {
        if (!r_t.mAlive || !r_t.mInside || !Real(r_t.mV))
            continue;
        const double r = Ratio(r_t.mV);
        if (r < 0.0) {
            ++mRep.mNumForced;
            continue;
        }
        if (!(r > mBound))
            continue;
        std::int64_t apex = ForcedBy(r_t.mV);
        if (apex < 0 && !mBall.empty()) {
            // Corner lopping forces a bad triangle whose point lies in a
            // ball at the end; Lemma 7 has its vertices within twice the
            // radius of the corner.
            const Vec2 c = CentreOf(r_t.mV);
            if (std::isfinite(c[0]) && std::isfinite(c[1]))
                apex = BallOf(c);
            if (apex >= 0) {
                auto& r_info = mLopInfo.at(apex);
                ++r_info.mNumForced;
                r_info.mMinForcedDeg = std::min(r_info.mMinForcedDeg, MinAngleDeg(r_t.mV));
                for (const std::int64_t v : r_t.mV)
                    if (Dist(apex, v, P(v)) > 2.0 * mBall.at(apex) * (1.0 + kCdt2dLopTol))
                        r_info.mWithinTwoRadii = false;
            }
        } else if (apex < 0) {
            const auto it = mDeclined.find(cdt2d_sorted(r_t.mV));
            if (it != mDeclined.end())
                apex = it->second;
        }
        if (apex < 0) {
            ++mRep.mNumUnresolved;
            continue;
        }
        ++mRep.mNumForced;
        mRep.mForced.emplace_back(cdt2d_sorted(r_t.mV), apex);
        const double m = MinAngleDeg(r_t.mV);
        auto [it, fresh] = forced_min.emplace(apex, m);
        if (!fresh)
            it->second = std::min(it->second, m);
    }
    std::sort(mRep.mForced.begin(), mRep.mForced.end());
    for (Cdt2dSmallAngle& r_a : mRep.mSmallAngles) {
        const auto it = forced_min.find(r_a.mVertex);
        if (it != forced_min.end())
            r_a.mMinForcedDeg = it->second;
    }
    for (const auto& r_l : mLopInfo)
        mRep.mLopped.push_back(r_l.second);
    mRep.mRequests.assign(mRequests.begin(), mRequests.end());
    return mRep;
}

double cdt2d_lop_bound_deg(double DeltaT) {
    constexpr double k_pi = 3.14159265358979323846;
    if (!(DeltaT >= 0.0) || DeltaT > k_pi / 3.0)
        return 0.0;
    const double t = std::tan(0.5 * DeltaT);
    const double eta = (1.0 + t) / (1.0 - t);
    const double mu = 1.0 / (1.0 - t);
    // zeta tends to sqrt 2 as the turning vanishes (0 / 0 at 0).
    const double zeta = DeltaT < 1e-6 ? std::sqrt(2.0)
                                      : std::sqrt(2.0) * std::sin(0.5 * DeltaT) /
                                            (std::sqrt(1.0 + std::sin(DeltaT)) - 1.0);
    const double k = std::min(std::asin(1.0 / (2.0 * std::sqrt(2.0) * mu)),
                              std::asin(std::min(1.0, 1.0 / (2.0 * zeta * eta))));
    return k * (180.0 / k_pi);
}

Cdt2dRefineReport Cdt2d::Refine(const Cdt2dRefineOptions& rOptions) {
    Cdt2dRefiner refiner(*this, rOptions);
    return refiner.Run();
}

std::vector<Cdt2d::Tri> Cdt2d::InsideTriangles() const {
    std::vector<Tri> out;
    for (const Triangle& r_t : mTris)
        if (r_t.mAlive && r_t.mInside && r_t.mV[0] >= kFirstVertex && r_t.mV[1] >= kFirstVertex &&
            r_t.mV[2] >= kFirstVertex)
            out.push_back(r_t.mV);
    return out;
}

std::vector<Cdt2d::Tri> Cdt2d::AllTriangles() const {
    std::vector<Tri> out;
    for (const Triangle& r_t : mTris)
        if (r_t.mAlive)
            out.push_back(r_t.mV);
    return out;
}

std::vector<std::array<std::int64_t, 2>> cdt2d_find_crossings(
    const std::vector<Cdt2d::Vec2>& rPoints,
    const std::vector<std::array<std::int64_t, 2>>& rSegments) {
    using Vec2 = Cdt2d::Vec2;
    std::vector<std::array<std::int64_t, 2>> out;
    const std::size_t ns = rSegments.size();
    if (ns < 2)
        return out;
    double lo[2] = {std::numeric_limits<double>::infinity(),
                    std::numeric_limits<double>::infinity()};
    double hi[2] = {-lo[0], -lo[1]};
    for (const auto& r_s : rSegments)
        for (const std::int64_t v : r_s)
            for (int k = 0; k < 2; ++k) {
                lo[k] = std::min(lo[k], rPoints[static_cast<std::size_t>(v)][k]);
                hi[k] = std::max(hi[k], rPoints[static_cast<std::size_t>(v)][k]);
            }
    // Square cells, about sqrt(n) of them along the longer side; a segment
    // is in every cell its box touches, so two that meet share a cell.
    const double per = std::ceil(std::sqrt(double(ns)));
    const double side = std::max(hi[0] - lo[0], hi[1] - lo[1]);
    const double cell = side > 0.0 && std::isfinite(side) ? side / per : 1.0;
    const auto index = [&](double X, int K) {
        const double c = std::floor((X - lo[K]) / cell);
        return static_cast<std::int64_t>(std::clamp(c, 0.0, per));
    };
    const auto nc = static_cast<std::int64_t>(per) + 1;
    std::vector<std::vector<std::uint32_t>> cells(static_cast<std::size_t>(nc * nc));
    for (std::size_t s = 0; s < ns; ++s) {
        const Vec2& a = rPoints[static_cast<std::size_t>(rSegments[s][0])];
        const Vec2& b = rPoints[static_cast<std::size_t>(rSegments[s][1])];
        const std::int64_t i0 = index(std::min(a[0], b[0]), 0);
        const std::int64_t i1 = index(std::max(a[0], b[0]), 0);
        const std::int64_t j0 = index(std::min(a[1], b[1]), 1);
        const std::int64_t j1 = index(std::max(a[1], b[1]), 1);
        for (std::int64_t j = j0; j <= j1; ++j)
            for (std::int64_t i = i0; i <= i1; ++i)
                cells[static_cast<std::size_t>(j * nc + i)].push_back(
                    static_cast<std::uint32_t>(s));
    }
    std::vector<std::uint64_t> pairs;
    for (const auto& r_c : cells)
        for (std::size_t x = 0; x < r_c.size(); ++x)
            for (std::size_t y = x + 1; y < r_c.size(); ++y)
                pairs.push_back((static_cast<std::uint64_t>(r_c[x]) << 32) | r_c[y]);
    std::sort(pairs.begin(), pairs.end());
    pairs.erase(std::unique(pairs.begin(), pairs.end()), pairs.end());
    for (const std::uint64_t key : pairs) {
        const auto s = static_cast<std::size_t>(key >> 32);
        const auto t = static_cast<std::size_t>(key & 0xffffffffull);
        if (cdt2d_segments_meet(rPoints, rSegments[s], rSegments[t]))
            out.push_back({static_cast<std::int64_t>(s), static_cast<std::int64_t>(t)});
    }
    return out;
}

Cdt2dPolygon cdt2d_triangulate_polygon(const std::vector<std::vector<Cdt2d::Vec2>>& rLoops,
                                       double TargetSize, double MergeTolerance) {
    return cdt2d_triangulate_polygon(rLoops, TargetSize, MergeTolerance, Cdt2dPolygonQuality{});
}

Cdt2dPolygon cdt2d_triangulate_polygon(const std::vector<std::vector<Cdt2d::Vec2>>& rLoops,
                                       double TargetSize, double MergeTolerance,
                                       const Cdt2dPolygonQuality& rQuality) {
    using Vec2 = Cdt2d::Vec2;
    Cdt2dPolygon out;
    double lo[2] = {std::numeric_limits<double>::infinity(),
                    std::numeric_limits<double>::infinity()};
    double hi[2] = {-lo[0], -lo[1]};
    std::size_t total = 0;
    for (const auto& r_loop : rLoops) {
        for (const Vec2& r_p : r_loop) {
            if (!std::isfinite(r_p[0]) || !std::isfinite(r_p[1]))
                throw std::invalid_argument("meshio++: triangulate_polygon_2d: a non-finite point");
            if (!cdt2d_in_range(r_p[0]) || !cdt2d_in_range(r_p[1]))
                throw std::invalid_argument(
                    "meshio++: triangulate_polygon_2d: a coordinate outside the range the exact "
                    "predicates hold in (zero, or a magnitude in [2^-128, 2^128])");
            for (int k = 0; k < 2; ++k) {
                lo[k] = std::min(lo[k], r_p[k]);
                hi[k] = std::max(hi[k], r_p[k]);
            }
        }
        total += r_loop.size();
    }
    if (total < 3) {
        out.mStatus = Cdt2dStatus::Degenerate;
        out.mWhy = "fewer than three boundary points";
        return out;
    }
    const double diag = std::hypot(hi[0] - lo[0], hi[1] - lo[1]);
    const double tol = MergeTolerance >= 0.0 ? MergeTolerance : 1e-12 * diag;

    // Merge near-duplicates (first point wins), then split long segments.
    std::vector<Vec2>& r_pts = out.mPoints;
    Cdt2dGrid merge(lo[0], lo[1], tol > 0.0 ? tol : std::max(diag, 1.0));
    const auto add_point = [&](const Vec2& rP) -> std::int64_t {
        if (tol > 0.0) {
            const auto [ci, cj] = merge.Cell(rP[0], rP[1]);
            for (std::int64_t di = -1; di <= 1; ++di)
                for (std::int64_t dj = -1; dj <= 1; ++dj)
                    if (const auto* p_ids = merge.At(Cdt2dGrid::Join(ci + di, cj + dj)))
                        for (const std::int64_t id : *p_ids) {
                            const Vec2& r_q = r_pts[static_cast<std::size_t>(id)];
                            if (std::hypot(r_q[0] - rP[0], r_q[1] - rP[1]) <= tol)
                                return id;
                        }
        }
        const auto id = static_cast<std::int64_t>(r_pts.size());
        r_pts.push_back(rP);
        if (tol > 0.0)
            merge.Add(rP[0], rP[1], id);
        return id;
    };
    std::vector<std::array<std::int64_t, 2>> segments;
    for (const auto& r_loop : rLoops) {
        std::vector<std::int64_t> ids;
        for (const Vec2& r_p : r_loop) {
            const std::int64_t id = add_point(r_p);
            if (ids.empty() || ids.back() != id)
                ids.push_back(id);
        }
        while (ids.size() > 1 && ids.front() == ids.back())
            ids.pop_back();
        if (ids.size() < 3)
            continue;
        for (std::size_t i = 0; i < ids.size(); ++i) {
            const std::int64_t a = ids[i];
            const std::int64_t b = ids[(i + 1) % ids.size()];
            const Vec2 pa = r_pts[static_cast<std::size_t>(a)];
            const Vec2 pb = r_pts[static_cast<std::size_t>(b)];
            const double len = std::hypot(pb[0] - pa[0], pb[1] - pa[1]);
            std::int64_t n = 1;
            if (TargetSize > 0.0)
                n = std::max<std::int64_t>(1,
                                           static_cast<std::int64_t>(std::ceil(len / TargetSize)));
            std::int64_t prev = a;
            for (std::int64_t k = 1; k <= n; ++k) {
                std::int64_t id = b;
                if (k < n) {
                    const double s = double(k) / double(n);
                    id = add_point({pa[0] + s * (pb[0] - pa[0]), pa[1] + s * (pb[1] - pa[1])});
                }
                if (id != prev)
                    segments.push_back({prev, id});
                prev = id;
            }
        }
    }
    out.mNumBoundary = static_cast<std::int64_t>(r_pts.size());
    if (segments.size() < 3) {
        out.mStatus = Cdt2dStatus::Degenerate;
        out.mWhy = "fewer than three distinct boundary segments";
        return out;
    }
    bool flat = true;
    for (std::size_t i = 2; i < r_pts.size() && flat; ++i)
        flat = orient2d(r_pts[0].data(), r_pts[1].data(), r_pts[i].data()) == 0.0;
    if (flat) {
        out.mStatus = Cdt2dStatus::Degenerate;
        out.mWhy = "the boundary points are collinear";
        return out;
    }

    Cdt2d cdt(lo[0], lo[1], hi[0], hi[1]);
    std::vector<std::int64_t> vid = cdt.InsertBatch(r_pts);
    std::vector<std::int64_t> point_of(static_cast<std::size_t>(cdt.NumVertices()), -1);
    for (std::size_t i = 0; i < vid.size(); ++i) {
        // Exact duplicates left by the merge (tolerance 0) share a vertex.
        if (point_of[static_cast<std::size_t>(vid[i])] < 0)
            point_of[static_cast<std::size_t>(vid[i])] = static_cast<std::int64_t>(i);
    }
    for (const auto& r_s : segments) {
        const std::int64_t a = vid[static_cast<std::size_t>(r_s[0])];
        const std::int64_t b = vid[static_cast<std::size_t>(r_s[1])];
        if (a == b)
            continue;
        const Cdt2dStatus st = cdt.InsertSegment(a, b);
        if (st != Cdt2dStatus::Ok) {
            out.mStatus = st;
            out.mWhy = cdt2d_status_name(st);
            out.mTriangles.clear();
            return out;
        }
    }
    const Cdt2dStatus st = cdt.Classify();
    if (st != Cdt2dStatus::Ok) {
        out.mStatus = st;
        out.mWhy = cdt2d_status_name(st);
        return out;
    }

    if (TargetSize > 0.0) {
        // An equilateral lattice, kept farther than h / 2 from the boundary
        // points and every point already kept, and than 0.4 h from the
        // boundary segments (closer would leave edges near h / 2 at the
        // middle of a segment; farther, gaps longer than 1.5 h by it).
        const double h = TargetSize;
        const double r2 = 0.25 * h * h;
        const double rs2 = 0.16 * h * h;
        Cdt2dGrid near(lo[0], lo[1], h);
        for (std::size_t i = 0; i < r_pts.size(); ++i)
            near.Add(r_pts[i][0], r_pts[i][1], static_cast<std::int64_t>(i));
        Cdt2dGrid segs(lo[0], lo[1], h);
        for (std::size_t s = 0; s < segments.size(); ++s) {
            const Vec2& r_a = r_pts[static_cast<std::size_t>(segments[s][0])];
            const Vec2& r_b = r_pts[static_cast<std::size_t>(segments[s][1])];
            segs.Add(0.5 * (r_a[0] + r_b[0]), 0.5 * (r_a[1] + r_b[1]),
                     static_cast<std::int64_t>(s));
        }
        const double dy = h * std::sqrt(3.0) / 2.0;
        const auto rows = static_cast<std::int64_t>(std::floor((hi[1] - lo[1]) / dy));
        const auto cols = static_cast<std::int64_t>(std::floor((hi[0] - lo[0]) / h)) + 1;
        for (std::int64_t j = 1; j <= rows; ++j) {
            const double v = lo[1] + double(j) * dy;
            for (std::int64_t i = 0; i <= cols; ++i) {
                const Vec2 p{lo[0] + (double(i) + 0.5 * double(j % 2)) * h, v};
                if (p[0] <= lo[0] || p[0] >= hi[0])
                    continue;
                const auto [ci, cj] = near.Cell(p[0], p[1]);
                bool close = false;
                for (std::int64_t di = -1; di <= 1 && !close; ++di)
                    for (std::int64_t dj = -1; dj <= 1 && !close; ++dj) {
                        if (const auto* p_ids = near.At(Cdt2dGrid::Join(ci + di, cj + dj)))
                            for (const std::int64_t id : *p_ids) {
                                const Vec2& r_q = r_pts[static_cast<std::size_t>(id)];
                                const double ex = r_q[0] - p[0];
                                const double ey = r_q[1] - p[1];
                                if (ex * ex + ey * ey < r2) {
                                    close = true;
                                    break;
                                }
                            }
                        if (close)
                            break;
                        if (const auto* p_ids = segs.At(Cdt2dGrid::Join(ci + di, cj + dj)))
                            for (const std::int64_t s : *p_ids) {
                                const auto& r_s = segments[static_cast<std::size_t>(s)];
                                if (cdt2d_segment_distance2(
                                        p, r_pts[static_cast<std::size_t>(r_s[0])],
                                        r_pts[static_cast<std::size_t>(r_s[1])]) < rs2) {
                                    close = true;
                                    break;
                                }
                            }
                    }
                if (close)
                    continue;
                const std::int64_t v_id = cdt.InsertSteiner(p);
                if (v_id < 0)
                    continue;
                const auto id = static_cast<std::int64_t>(r_pts.size());
                r_pts.push_back(p);
                near.Add(p[0], p[1], id);
                if (point_of.size() <= static_cast<std::size_t>(v_id))
                    point_of.resize(static_cast<std::size_t>(v_id) + 1, -1);
                point_of[static_cast<std::size_t>(v_id)] = id;
            }
        }
    }

    if (rQuality.mMinAngleDeg > 0.0) {
        Cdt2dRefineOptions ro;
        ro.mMinAngleDeg = rQuality.mMinAngleDeg;
        ro.mMaxPoints = rQuality.mMaxPoints;
        const std::int64_t first = cdt.NumVertices();
        out.mRefine = cdt.Refine(ro);
        if (point_of.size() < static_cast<std::size_t>(cdt.NumVertices()))
            point_of.resize(static_cast<std::size_t>(cdt.NumVertices()), -1);
        for (std::int64_t v = first; v < cdt.NumVertices(); ++v) {
            point_of[static_cast<std::size_t>(v)] = static_cast<std::int64_t>(r_pts.size());
            r_pts.push_back(cdt.Point(v));
        }
        // The report in point ids, like the rest of the result.
        for (Cdt2dSmallAngle& r_a : out.mRefine.mSmallAngles)
            r_a.mVertex = point_of[static_cast<std::size_t>(r_a.mVertex)];
        segments.clear();
        for (const auto& r_e : cdt.ConstrainedEdges())
            segments.push_back({point_of[static_cast<std::size_t>(r_e[0])],
                                point_of[static_cast<std::size_t>(r_e[1])]});
    }

    for (const auto& r_t : cdt.InsideTriangles()) {
        Cdt2d::Tri tri{};
        for (int k = 0; k < 3; ++k)
            tri[k] = point_of[static_cast<std::size_t>(r_t[k])];
        out.mTriangles.push_back(tri);
    }
    for (const auto& r_s : segments)
        if (r_s[0] != r_s[1])
            out.mSegments.push_back(r_s);
    if (out.mTriangles.empty()) {
        out.mStatus = Cdt2dStatus::Degenerate;
        out.mWhy = "the boundary encloses no area";
    }
    return out;
}

}  // namespace detail
}  // namespace meshioplusplus

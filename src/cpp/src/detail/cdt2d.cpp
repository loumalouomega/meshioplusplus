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

}  // namespace

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
        if (InCircle(r_t.mV[0], r_t.mV[1], r_t.mV[2], d) > 0.0) {
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
    std::vector<std::int64_t> removed{rLoc.mTri};
    std::vector<Tri> fresh;
    std::int64_t split_a = -1;
    std::int64_t split_b = -1;
    bool split_constrained = false;
    const auto& r_t = mTris[static_cast<std::size_t>(rLoc.mTri)];
    if (rLoc.mKind == 0) {
        for (int e = 0; e < 3; ++e)
            fresh.push_back({r_t.mV[cdt2d_next(e)], r_t.mV[cdt2d_prev(e)], v});
    } else {
        const int e = rLoc.mIndex;
        split_a = r_t.mV[cdt2d_next(e)];
        split_b = r_t.mV[cdt2d_prev(e)];
        split_constrained = r_t.mC[e];
        const std::int64_t u = r_t.mN[e];
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

std::int64_t Cdt2d::InsertSteiner(const Vec2& rP) {
    if (!mClassified)
        cdt2d_internal("InsertSteiner before Classify");
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

Cdt2dPolygon cdt2d_triangulate_polygon(const std::vector<std::vector<Cdt2d::Vec2>>& rLoops,
                                       double TargetSize, double MergeTolerance) {
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

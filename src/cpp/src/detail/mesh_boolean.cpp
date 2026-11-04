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
// The surface-boolean spike (roadmap 7.8). See detail/mesh_boolean.hpp.
// Anonymous helpers are prefixed `mb_` (the amalgamation concatenates every
// translation unit into one).

// System includes
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <map>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/detail/surface_distance.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/operations/sdf.hpp"

// Project includes (private, not installed)
#include "expansion.hpp"
#include "mesh_boolean.hpp"
#include "predicates.hpp"
#include "tri_pairs.hpp"
#include "tri_tri.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

using MbVec = std::array<double, 3>;
using MbTri = std::array<std::int64_t, 3>;

/// A point as exact homogeneous coordinates; `mApprox` its rounded position.
struct MbPoint {
    std::array<Expansion, 3> mX;
    Expansion mW;
    MbVec mApprox{};
    bool mExplicit = false;
};

/// The failure the spike reports instead of a result.
struct MbFailure {
    std::string mWhy;
};

struct MbSurface {
    std::vector<MbVec> mP;
    std::vector<MbTri> mT;
};

double mb_seconds_since(const std::chrono::steady_clock::time_point& rT0) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - rT0).count();
}

MbSurface mb_surface(const Mesh& rMesh, const char* pName) {
    MbSurface s;
    const std::size_t np = rMesh.NumPoints();
    const NDArray& r_p = rMesh.Points();
    if (np > 0 && (r_p.Ndim() != 2 || r_p.Shape()[1] != 3))
        throw MbFailure{std::string(pName) + ": points must be three-dimensional"};
    s.mP.resize(np);
    for (std::size_t i = 0; i < np; ++i)
        for (std::size_t k = 0; k < 3; ++k) {
            s.mP[i][k] = read_double(r_p, 3 * i + k);
            if (!std::isfinite(s.mP[i][k]))
                throw MbFailure{std::string(pName) + ": a non-finite point"};
        }
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto cells = rMesh.Cells(b);
        if (cells.Type() != "triangle" || cells.IsRagged())
            throw MbFailure{std::string(pName) + ": not a triangle surface ('" +
                            std::string(cells.Type()) + "' block)"};
        for (std::size_t c = 0; c < cells.NumCells(); ++c)
            s.mT.push_back({read_int(cells.Conn(), 3 * c), read_int(cells.Conn(), 3 * c + 1),
                            read_int(cells.Conn(), 3 * c + 2)});
    }
    if (s.mT.empty())
        throw MbFailure{std::string(pName) + ": no triangle"};
    const SurfaceQuality q = surface_watertight_check(rMesh);
    if (!q.mWatertight)
        throw MbFailure{std::string(pName) + ": not closed, manifold and oriented"};
    return s;
}

MbPoint mb_explicit(const MbVec& rP) {
    MbPoint q;
    for (std::size_t k = 0; k < 3; ++k)
        q.mX[k] = exp_of(rP[k]);
    q.mW = exp_of(1.0);
    q.mApprox = rP;
    q.mExplicit = true;
    return q;
}

/// The exact normal `(b - a) x (c - a)`.
std::array<Expansion, 3> mb_normal(const MbVec& rA, const MbVec& rB, const MbVec& rC) {
    std::array<Expansion, 3> u, v, n;
    for (std::size_t k = 0; k < 3; ++k) {
        u[k] = exp_diff(rB[k], rA[k]);
        v[k] = exp_diff(rC[k], rA[k]);
    }
    for (std::size_t k = 0; k < 3; ++k) {
        const std::size_t i = (k + 1) % 3, j = (k + 2) % 3;
        n[k] = exp_sub(exp_mul(u[i], v[j]), exp_mul(u[j], v[i]));
    }
    return n;
}

/// The edge `pq` meets the plane of `abc` (the caller has checked that `p`
/// and `q` lie strictly on opposite sides).
MbPoint mb_lpi(const MbVec& rP, const MbVec& rQ, const MbVec& rA, const MbVec& rB,
               const MbVec& rC) {
    const std::array<Expansion, 3> n = mb_normal(rA, rB, rC);
    Expansion dp, dq;
    for (std::size_t k = 0; k < 3; ++k) {
        dp = exp_sum(dp, exp_mul(n[k], exp_diff(rP[k], rA[k])));
        dq = exp_sum(dq, exp_mul(n[k], exp_diff(rQ[k], rA[k])));
    }
    MbPoint out;
    for (std::size_t k = 0; k < 3; ++k)
        out.mX[k] = exp_sub(exp_scale(dp, rQ[k]), exp_scale(dq, rP[k]));
    out.mW = exp_sub(dp, dq);
    const double w = exp_estimate(out.mW);
    for (std::size_t k = 0; k < 3; ++k)
        out.mApprox[k] = exp_estimate(out.mX[k]) / w;
    return out;
}

bool mb_equal(const MbPoint& rP, const MbPoint& rQ) {
    if (rP.mExplicit && rQ.mExplicit)
        return rP.mApprox == rQ.mApprox;
    for (std::size_t k = 0; k < 3; ++k)
        if (exp_sign(exp_sub(exp_mul(rP.mX[k], rQ.mW), exp_mul(rQ.mX[k], rP.mW))) != 0)
            return false;
    return true;
}

/// The sign of `P_k - Q_k`, exactly.
int mb_cmp_coord(const MbPoint& rP, const MbPoint& rQ, std::size_t K) {
    if (rP.mExplicit && rQ.mExplicit)
        return rP.mApprox[K] < rQ.mApprox[K] ? -1 : (rP.mApprox[K] > rQ.mApprox[K] ? 1 : 0);
    const int s = exp_sign(exp_sub(exp_mul(rP.mX[K], rQ.mW), exp_mul(rQ.mX[K], rP.mW)));
    return s * exp_sign(rP.mW) * exp_sign(rQ.mW);
}

/// The sign of `orient2d` of three points with axis `Ax` dropped.
int mb_orient2d(const MbPoint& r1, const MbPoint& r2, const MbPoint& r3, std::size_t Ax) {
    const std::size_t u = (Ax + 1) % 3, v = (Ax + 2) % 3;
    if (r1.mExplicit && r2.mExplicit && r3.mExplicit) {
        const double a[2] = {r1.mApprox[u], r1.mApprox[v]};
        const double b[2] = {r2.mApprox[u], r2.mApprox[v]};
        const double c[2] = {r3.mApprox[u], r3.mApprox[v]};
        const double o = orient2d(a, b, c);
        return o > 0.0 ? 1 : (o < 0.0 ? -1 : 0);
    }
    const auto minor = [&](const MbPoint& rI, const MbPoint& rJ) {
        return exp_sub(exp_mul(rI.mX[v], rJ.mW), exp_mul(rJ.mX[v], rI.mW));
    };
    Expansion det = exp_mul(r1.mX[u], minor(r2, r3));
    det = exp_sub(det, exp_mul(r2.mX[u], minor(r1, r3)));
    det = exp_sum(det, exp_mul(r3.mX[u], minor(r1, r2)));
    return exp_sign(det) * exp_sign(r1.mW) * exp_sign(r2.mW) * exp_sign(r3.mW);
}

int mb_sign(double V) {
    return V > 0.0 ? 1 : (V < 0.0 ? -1 : 0);
}

/// The points, deduplicated exactly: a bucket grid on the rounded positions,
/// then an exact comparison with the candidates in the 27 cells around.
class MbStore {
public:
    explicit MbStore(double Cell) : mCell(Cell) {}

    std::int64_t Add(MbPoint&& rP) {
        const auto key = Key(rP.mApprox);
        for (std::int64_t dx = -1; dx <= 1; ++dx)
            for (std::int64_t dy = -1; dy <= 1; ++dy)
                for (std::int64_t dz = -1; dz <= 1; ++dz) {
                    const auto it = mGrid.find(Hash({key[0] + dx, key[1] + dy, key[2] + dz}));
                    if (it == mGrid.end())
                        continue;
                    for (const std::int64_t id : it->second)
                        if (mb_equal(mPoints[static_cast<std::size_t>(id)], rP))
                            return id;
                }
        const auto id = static_cast<std::int64_t>(mPoints.size());
        mGrid[Hash(key)].push_back(id);
        mPoints.push_back(std::move(rP));
        return id;
    }

    const MbPoint& operator[](std::int64_t Id) const {
        return mPoints[static_cast<std::size_t>(Id)];
    }
    std::size_t Size() const { return mPoints.size(); }

private:
    std::array<std::int64_t, 3> Key(const MbVec& rP) const {
        return {static_cast<std::int64_t>(std::floor(rP[0] / mCell)),
                static_cast<std::int64_t>(std::floor(rP[1] / mCell)),
                static_cast<std::int64_t>(std::floor(rP[2] / mCell))};
    }
    static std::uint64_t Hash(const std::array<std::int64_t, 3>& rK) {
        return static_cast<std::uint64_t>(rK[0]) * 73856093ULL ^
               static_cast<std::uint64_t>(rK[1]) * 19349663ULL ^
               static_cast<std::uint64_t>(rK[2]) * 83492791ULL;
    }

    double mCell;
    std::vector<MbPoint> mPoints;
    std::unordered_map<std::uint64_t, std::vector<std::int64_t>> mGrid;
};

/// One triangle of a surface as the arrangement cuts it.
struct MbCut {
    std::vector<std::int64_t> mPoints;                             ///< global ids
    std::vector<std::pair<std::int64_t, std::int64_t>> mSegments;  ///< global ids
};

std::uint64_t mb_edge_key(std::int64_t A, std::int64_t B) {
    if (A > B)
        std::swap(A, B);
    return (static_cast<std::uint64_t>(A) << 32) ^ static_cast<std::uint64_t>(B);
}

/// The dropped axis of a triangle (its normal's largest component) and the
/// sign that makes its own corners turn positively in that projection.
std::pair<std::size_t, int> mb_axis(const MbVec& rA, const MbVec& rB, const MbVec& rC) {
    const MbVec u{rB[0] - rA[0], rB[1] - rA[1], rB[2] - rA[2]};
    const MbVec v{rC[0] - rA[0], rC[1] - rA[1], rC[2] - rA[2]};
    const MbVec n{u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]};
    std::size_t ax = 0;
    for (std::size_t k = 1; k < 3; ++k)
        if (std::fabs(n[k]) > std::fabs(n[ax]))
            ax = k;
    const std::size_t i = (ax + 1) % 3, j = (ax + 2) % 3;
    const double a[2] = {rA[i], rA[j]}, b[2] = {rB[i], rB[j]}, c[2] = {rC[i], rC[j]};
    return {ax, mb_sign(orient2d(a, b, c))};
}

/// A small constrained triangulation of one triangle in its plane: points
/// located exactly, segments recovered by flips. `mT` holds local ids, every
/// triangle positive (`Orient` > 0).
class MbTriangulation {
public:
    MbTriangulation(const MbStore& rStore, std::size_t Ax, int Sign,
                    const std::array<std::int64_t, 3>& rCorners)
        : mStore(rStore), mAx(Ax), mSign(Sign) {
        mV = {rCorners[0], rCorners[1], rCorners[2]};
        mT.push_back({0, 1, 2});
    }

    int Orient(int A, int B, int C) const {
        return mSign * mb_orient2d(mStore[mV[static_cast<std::size_t>(A)]],
                                   mStore[mV[static_cast<std::size_t>(B)]],
                                   mStore[mV[static_cast<std::size_t>(C)]], mAx);
    }

    /// The local id of global point `Id`, inserting it.
    int Insert(std::int64_t Id) {
        for (std::size_t i = 0; i < mV.size(); ++i)
            if (mV[i] == Id)
                return static_cast<int>(i);
        const int p = static_cast<int>(mV.size());
        mV.push_back(Id);
        for (std::size_t t = 0; t < mT.size(); ++t) {
            const auto tri = mT[t];
            int o[3];
            for (int k = 0; k < 3; ++k)
                o[k] = Orient(tri[k], tri[(k + 1) % 3], p);
            if (o[0] < 0 || o[1] < 0 || o[2] < 0)
                continue;
            const int zeros = (o[0] == 0) + (o[1] == 0) + (o[2] == 0);
            if (zeros == 0) {
                mT[t] = {tri[0], tri[1], p};
                mT.push_back({tri[1], tri[2], p});
                mT.push_back({tri[2], tri[0], p});
                return p;
            }
            if (zeros == 1) {
                const int k = o[0] == 0 ? 0 : (o[1] == 0 ? 1 : 2);
                SplitEdge(tri[k], tri[(k + 1) % 3], p);
                return p;
            }
            throw MbFailure{"a point coincides with a vertex it was not identified with"};
        }
        throw MbFailure{"a point lies outside the triangle it was found in"};
    }

    /// Recover the segment between local `A` and `B` (both inserted), through
    /// the points lying on it.
    void Constrain(int A, int B, int Depth = 0) {
        if (A == B)
            return;
        if (Depth > 64)
            throw MbFailure{"segment splitting did not terminate"};
        // Points strictly inside the segment split it.
        // The projected coordinate the segment moves most along.
        const std::size_t u = (mAx + 1) % 3, v = (mAx + 2) % 3;
        const std::size_t axis = std::fabs(Pt(A).mApprox[u] - Pt(B).mApprox[u]) >=
                                         std::fabs(Pt(A).mApprox[v] - Pt(B).mApprox[v])
                                     ? u
                                     : v;
        const int dir = mb_cmp_coord(Pt(B), Pt(A), axis);
        if (dir == 0)
            throw MbFailure{"a segment's ends coincide in its triangle's plane"};
        std::vector<int> on;
        for (int q = 0; q < static_cast<int>(mV.size()); ++q) {
            if (q == A || q == B || Orient(A, B, q) != 0)
                continue;
            if (mb_cmp_coord(Pt(q), Pt(A), axis) == dir && mb_cmp_coord(Pt(B), Pt(q), axis) == dir)
                on.push_back(q);
        }
        if (!on.empty()) {
            std::sort(on.begin(), on.end(),
                      [&](int P, int Q) { return mb_cmp_coord(Pt(P), Pt(Q), axis) == dir; });
            int prev = A;
            for (const int q : on) {
                Constrain(prev, q, Depth + 1);
                prev = q;
            }
            Constrain(prev, B, Depth + 1);
            return;
        }
        for (int guard = 0; guard < 100000; ++guard) {
            if (HasEdge(A, B)) {
                mConstrained.insert(
                    mb_edge_key(mV[static_cast<std::size_t>(A)], mV[static_cast<std::size_t>(B)]));
                return;
            }
            bool flipped = false;
            for (std::size_t t = 0; t < mT.size() && !flipped; ++t)
                for (int e = 0; e < 3 && !flipped; ++e) {
                    const int c = mT[t][e], d = mT[t][(e + 1) % 3], x = mT[t][(e + 2) % 3];
                    if (c == A || c == B || d == A || d == B)
                        continue;
                    if (Orient(A, B, c) * Orient(A, B, d) >= 0 ||
                        Orient(c, d, A) * Orient(c, d, B) >= 0)
                        continue;
                    // The triangle across (c, d): (d, c, y).
                    for (std::size_t u = 0; u < mT.size(); ++u) {
                        int y = -1;
                        for (int f = 0; f < 3; ++f)
                            if (mT[u][f] == d && mT[u][(f + 1) % 3] == c)
                                y = mT[u][(f + 2) % 3];
                        if (y < 0)
                            continue;
                        if (Orient(x, y, c) * Orient(x, y, d) >= 0)
                            break;  // not convex: another crossing edge first
                        mT[t] = {c, y, x};
                        mT[u] = {y, d, x};
                        flipped = true;
                        break;
                    }
                }
            if (!flipped)
                throw MbFailure{"a segment could not be recovered by flips"};
        }
        throw MbFailure{"segment recovery did not finish"};
    }

    /// Every triangle positive, or the step named fails.
    void Check(const char* pStep) const {
        for (const auto& r_t : mT)
            if (Orient(r_t[0], r_t[1], r_t[2]) <= 0)
                throw MbFailure{std::string("a non-positive piece after ") + pStep};
    }

    std::vector<MbTri> Triangles() const {
        std::vector<MbTri> out;
        out.reserve(mT.size());
        for (const auto& r_t : mT)
            out.push_back({mV[static_cast<std::size_t>(r_t[0])],
                           mV[static_cast<std::size_t>(r_t[1])],
                           mV[static_cast<std::size_t>(r_t[2])]});
        return out;
    }

    const std::unordered_set<std::uint64_t>& Constrained() const { return mConstrained; }

private:
    const MbPoint& Pt(int L) const { return mStore[mV[static_cast<std::size_t>(L)]]; }

    bool HasEdge(int A, int B) const {
        for (const auto& r_t : mT)
            for (int k = 0; k < 3; ++k)
                if ((r_t[k] == A && r_t[(k + 1) % 3] == B) ||
                    (r_t[k] == B && r_t[(k + 1) % 3] == A))
                    return true;
        return false;
    }

    void SplitEdge(int A, int B, int P) {
        // Every triangle with the edge (A, B) or (B, A) is split at P.
        const std::size_t n = mT.size();
        for (std::size_t t = 0; t < n; ++t)
            for (int k = 0; k < 3; ++k) {
                const int u = mT[t][k], v = mT[t][(k + 1) % 3], w = mT[t][(k + 2) % 3];
                if ((u == A && v == B) || (u == B && v == A)) {
                    mT[t] = {u, P, w};
                    mT.push_back({P, v, w});
                    break;
                }
            }
    }

    const MbStore& mStore;
    std::size_t mAx;
    int mSign;
    std::vector<std::int64_t> mV;
    std::vector<std::array<int, 3>> mT;
    std::unordered_set<std::uint64_t> mConstrained;
};

/// Whether the segment pq crosses triangle abc properly (exact); `rDegenerate`
/// is set when it touches an edge, a vertex or the plane at an end.
bool mb_ray_crosses(const MbVec& rP, const MbVec& rQ, const MbVec& rA, const MbVec& rB,
                    const MbVec& rC, bool& rDegenerate) {
    const int s1 = mb_sign(orient3d(rA.data(), rB.data(), rC.data(), rP.data()));
    const int s2 = mb_sign(orient3d(rA.data(), rB.data(), rC.data(), rQ.data()));
    if (s1 == s2 && s1 != 0)
        return false;
    if (s1 == 0 || s2 == 0) {
        rDegenerate = true;
        return false;
    }
    const int t1 = mb_sign(orient3d(rP.data(), rQ.data(), rA.data(), rB.data()));
    const int t2 = mb_sign(orient3d(rP.data(), rQ.data(), rB.data(), rC.data()));
    const int t3 = mb_sign(orient3d(rP.data(), rQ.data(), rC.data(), rA.data()));
    const bool pos = t1 >= 0 && t2 >= 0 && t3 >= 0, neg = t1 <= 0 && t2 <= 0 && t3 <= 0;
    if (!pos && !neg)
        return false;
    if (t1 == 0 || t2 == 0 || t3 == 0) {
        rDegenerate = true;
        return false;
    }
    return true;
}

/// Whether `rP` lies inside the closed surface `rS` (exact ray parity; a ray
/// that touches an edge, a vertex or a plane is cast again elsewhere).
bool mb_inside(const MbVec& rP, const MbSurface& rS, const MbVec& rLo, const MbVec& rHi,
               std::int64_t& rRetries) {
    static const double kDirs[8][3] = {
        {0.5773502691896258, 0.6123724356957945, 0.5400617248673217},
        {-0.4472135954999579, 0.7745966692414834, 0.4472135954999579},
        {0.8017837257372732, -0.2672612419124244, 0.5345224838248488},
        {0.1543033499620919, 0.6172133998483676, -0.7715167498104595},
        {-0.6666666666666666, -0.3333333333333333, 0.6666666666666666},
        {0.2357022603955158, -0.9428090415820634, 0.2357022603955158},
        {0.7071067811865476, 0.5, -0.5},
        {-0.3015113445777636, 0.3015113445777636, 0.9045340337332909}};
    double diag = 0.0;
    for (std::size_t k = 0; k < 3; ++k)
        diag += (rHi[k] - rLo[k]) * (rHi[k] - rLo[k]);
    const double len = 4.0 * std::sqrt(diag) + 1.0;
    for (const auto& r_d : kDirs) {
        const MbVec q{rP[0] + len * r_d[0], rP[1] + len * r_d[1], rP[2] + len * r_d[2]};
        MbVec lo{}, hi{};
        for (std::size_t k = 0; k < 3; ++k) {
            lo[k] = std::min(rP[k], q[k]);
            hi[k] = std::max(rP[k], q[k]);
        }
        bool degenerate = false;
        std::int64_t crossings = 0;
        for (const auto& r_t : rS.mT) {
            const MbVec& a = rS.mP[static_cast<std::size_t>(r_t[0])];
            const MbVec& b = rS.mP[static_cast<std::size_t>(r_t[1])];
            const MbVec& c = rS.mP[static_cast<std::size_t>(r_t[2])];
            bool outside = false;
            for (std::size_t k = 0; k < 3 && !outside; ++k)
                outside =
                    std::max({a[k], b[k], c[k]}) < lo[k] || std::min({a[k], b[k], c[k]}) > hi[k];
            if (outside)
                continue;
            if (mb_ray_crosses(rP, q, a, b, c, degenerate))
                ++crossings;
            if (degenerate)
                break;
        }
        if (!degenerate)
            return crossings % 2 == 1;
        ++rRetries;
    }
    throw MbFailure{"every ray from a patch's sample touched the other surface"};
}

/// The pieces of one surface after the arrangement, and their patches.
struct MbPieces {
    std::vector<MbTri> mT;             ///< global point ids
    std::vector<MbVec> mNormal;        ///< the parent triangle's normal, rounded
    std::vector<std::int64_t> mPatch;  ///< patch of each piece
    std::int64_t mNumPatches = 0;
};

void mb_patches(MbPieces& rPieces, const std::unordered_set<std::uint64_t>& rCurves) {
    std::unordered_map<std::uint64_t, std::vector<std::size_t>> edges;
    for (std::size_t t = 0; t < rPieces.mT.size(); ++t)
        for (int k = 0; k < 3; ++k)
            edges[mb_edge_key(rPieces.mT[t][k], rPieces.mT[t][(k + 1) % 3])].push_back(t);
    rPieces.mPatch.assign(rPieces.mT.size(), -1);
    std::vector<std::size_t> stack;
    for (std::size_t s = 0; s < rPieces.mT.size(); ++s) {
        if (rPieces.mPatch[s] >= 0)
            continue;
        const std::int64_t id = rPieces.mNumPatches++;
        rPieces.mPatch[s] = id;
        stack.push_back(s);
        while (!stack.empty()) {
            const std::size_t t = stack.back();
            stack.pop_back();
            for (int k = 0; k < 3; ++k) {
                const std::uint64_t key = mb_edge_key(rPieces.mT[t][k], rPieces.mT[t][(k + 1) % 3]);
                if (rCurves.count(key))
                    continue;
                for (const std::size_t u : edges[key])
                    if (rPieces.mPatch[u] < 0) {
                        rPieces.mPatch[u] = id;
                        stack.push_back(u);
                    }
            }
        }
    }
}

}  // namespace

double closed_surface_volume(const Mesh& rSurface) {
    const std::size_t np = rSurface.NumPoints();
    const NDArray& r_p = rSurface.Points();
    MbVec c{};
    for (std::size_t i = 0; i < np; ++i)
        for (std::size_t k = 0; k < 3; ++k)
            c[k] += read_double(r_p, 3 * i + k) / static_cast<double>(np);
    double vol = 0.0;
    for (std::size_t b = 0; b < rSurface.NumCellBlocks(); ++b) {
        const auto cells = rSurface.Cells(b);
        if (cells.Type() != "triangle")
            continue;
        for (std::size_t t = 0; t < cells.NumCells(); ++t) {
            MbVec q[3];
            for (std::size_t j = 0; j < 3; ++j) {
                const auto v = static_cast<std::size_t>(read_int(cells.Conn(), 3 * t + j));
                for (std::size_t k = 0; k < 3; ++k)
                    q[j][k] = read_double(r_p, 3 * v + k) - c[k];
            }
            vol += (q[0][0] * (q[1][1] * q[2][2] - q[1][2] * q[2][1]) -
                    q[0][1] * (q[1][0] * q[2][2] - q[1][2] * q[2][0]) +
                    q[0][2] * (q[1][0] * q[2][1] - q[1][1] * q[2][0])) /
                   6.0;
        }
    }
    return vol;
}

BooleanResult mesh_boolean_spike(const Mesh& rA, const Mesh& rB, BooleanOp Op) {
    BooleanResult result;
    BooleanReport& r_rep = result.mReport;
    try {
        const MbSurface sa = mb_surface(rA, "A");
        const MbSurface sb = mb_surface(rB, "B");

        // The bounding box of both, for the point grid and the rays.
        MbVec lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
        for (const MbSurface* p_s : {&sa, &sb})
            for (const MbVec& r_p : p_s->mP)
                for (std::size_t k = 0; k < 3; ++k) {
                    lo[k] = std::min(lo[k], r_p[k]);
                    hi[k] = std::max(hi[k], r_p[k]);
                }
        const double diag =
            std::sqrt((hi[0] - lo[0]) * (hi[0] - lo[0]) + (hi[1] - lo[1]) * (hi[1] - lo[1]) +
                      (hi[2] - lo[2]) * (hi[2] - lo[2]));
        MbStore store(std::max(diag * 1e-9, 1e-300));
        std::vector<std::int64_t> ida(sa.mP.size()), idb(sb.mP.size());
        for (std::size_t i = 0; i < sa.mP.size(); ++i)
            ida[i] = store.Add(mb_explicit(sa.mP[i]));
        for (std::size_t i = 0; i < sb.mP.size(); ++i)
            idb[i] = store.Add(mb_explicit(sb.mP[i]));

        // The candidate pairs.
        auto t0 = std::chrono::steady_clock::now();
        const TriangleSoup soup_a = build_triangle_soup(rA, "");
        const TriangleSoup soup_b = build_triangle_soup(rB, "");
        std::vector<std::vector<std::int64_t>> cand(sa.mT.size());
        soup_box_pairs(soup_a, &soup_b,
                       [&](std::size_t I, const std::vector<std::int64_t>& rJs) { cand[I] = rJs; });
        std::vector<std::pair<std::size_t, std::size_t>> pairs;
        for (std::size_t i = 0; i < sa.mT.size(); ++i) {
            r_rep.mNumCandidatePairs += static_cast<std::int64_t>(cand[i].size());
            const auto& ta = sa.mT[i];
            for (const std::int64_t jj : cand[i]) {
                const auto j = static_cast<std::size_t>(jj);
                const auto& tb = sb.mT[j];
                if (tri_tri_intersect(sa.mP[ta[0]].data(), sa.mP[ta[1]].data(), sa.mP[ta[2]].data(),
                                      sb.mP[tb[0]].data(), sb.mP[tb[1]].data(),
                                      sb.mP[tb[2]].data()))
                    pairs.emplace_back(i, j);
            }
        }
        r_rep.mNumIntersectingPairs = static_cast<std::int64_t>(pairs.size());
        r_rep.mSecondsPairs = mb_seconds_since(t0);

        // The arrangement: each pair's exact segment.
        t0 = std::chrono::steady_clock::now();
        std::vector<MbCut> cut_a(sa.mT.size()), cut_b(sb.mT.size());
        std::map<std::tuple<int, std::int64_t, std::int64_t, std::int64_t>, std::int64_t> lpi_ids;
        const auto lpi = [&](int Mesh_, const MbSurface& rS, std::int64_t P, std::int64_t Q,
                             const MbSurface& rO, std::size_t Tri) {
            const auto key = std::make_tuple(Mesh_, std::min(P, Q), std::max(P, Q),
                                             static_cast<std::int64_t>(Tri));
            const auto it = lpi_ids.find(key);
            if (it != lpi_ids.end())
                return it->second;
            const auto& t = rO.mT[Tri];
            const std::int64_t id = store.Add(
                mb_lpi(rS.mP[static_cast<std::size_t>(P)], rS.mP[static_cast<std::size_t>(Q)],
                       rO.mP[static_cast<std::size_t>(t[0])], rO.mP[static_cast<std::size_t>(t[1])],
                       rO.mP[static_cast<std::size_t>(t[2])]));
            ++r_rep.mNumImplicitPoints;
            lpi_ids.emplace(key, id);
            return id;
        };
        // Whether global point Id (on the plane of `rS`'s triangle Tri) lies in it.
        const auto in_tri = [&](std::int64_t Id, const MbSurface& rS,
                                const std::vector<std::int64_t>& rIds, std::size_t Tri) {
            const auto& t = rS.mT[Tri];
            const auto [ax, sign] = mb_axis(rS.mP[static_cast<std::size_t>(t[0])],
                                            rS.mP[static_cast<std::size_t>(t[1])],
                                            rS.mP[static_cast<std::size_t>(t[2])]);
            for (int k = 0; k < 3; ++k) {
                const int o =
                    sign * mb_orient2d(store[rIds[static_cast<std::size_t>(t[k])]],
                                       store[rIds[static_cast<std::size_t>(t[(k + 1) % 3])]],
                                       store[Id], ax);
                if (o < 0)
                    return false;
            }
            return true;
        };
        for (const auto& [i, j] : pairs) {
            const auto& ta = sa.mT[i];
            const auto& tb = sb.mT[j];
            int s_a[3], s_b[3];
            for (int k = 0; k < 3; ++k) {
                s_a[k] = mb_sign(orient3d(sb.mP[tb[0]].data(), sb.mP[tb[1]].data(),
                                          sb.mP[tb[2]].data(), sa.mP[ta[k]].data()));
                s_b[k] = mb_sign(orient3d(sa.mP[ta[0]].data(), sa.mP[ta[1]].data(),
                                          sa.mP[ta[2]].data(), sb.mP[tb[k]].data()));
            }
            if (s_a[0] == 0 && s_a[1] == 0 && s_a[2] == 0) {
                ++r_rep.mNumCoplanarPairs;
                continue;
            }
            std::vector<std::int64_t> pts;
            const auto collect = [&](const int* pS, const MbSurface& rS,
                                     const std::vector<std::int64_t>& rIds, const MbTri& rT,
                                     int Mesh_, const MbSurface& rO,
                                     const std::vector<std::int64_t>& rOIds, std::size_t OTri) {
                for (int k = 0; k < 3; ++k) {
                    const int l = (k + 1) % 3;
                    std::int64_t id = -1;
                    if (pS[k] == 0)
                        id = rIds[static_cast<std::size_t>(rT[k])];
                    else if (pS[k] * pS[l] < 0)
                        id = lpi(Mesh_, rS, rT[k], rT[l], rO, OTri);
                    if (id >= 0 && in_tri(id, rO, rOIds, OTri))
                        pts.push_back(id);
                }
            };
            collect(s_a, sa, ida, ta, 0, sb, idb, j);
            collect(s_b, sb, idb, tb, 1, sa, ida, i);
            std::sort(pts.begin(), pts.end());
            pts.erase(std::unique(pts.begin(), pts.end()), pts.end());
            if (pts.empty())
                continue;
            if (pts.size() > 2) {
                // Collinear along the planes' line: ordered by the exact
                // coordinate the line moves most along.
                const MbPoint& p0 = store[pts[0]];
                std::size_t axis = 0;
                double far = -1.0;
                for (std::size_t k = 0; k < 3; ++k)
                    for (const std::int64_t q : pts)
                        if (std::fabs(store[q].mApprox[k] - p0.mApprox[k]) > far) {
                            far = std::fabs(store[q].mApprox[k] - p0.mApprox[k]);
                            axis = k;
                        }
                std::sort(pts.begin(), pts.end(), [&](std::int64_t P, std::int64_t Q) {
                    return mb_cmp_coord(store[P], store[Q], axis) < 0;
                });
            }
            for (const std::int64_t q : pts) {
                cut_a[i].mPoints.push_back(q);
                cut_b[j].mPoints.push_back(q);
            }
            for (std::size_t k = 0; k + 1 < pts.size(); ++k) {
                cut_a[i].mSegments.emplace_back(pts[k], pts[k + 1]);
                cut_b[j].mSegments.emplace_back(pts[k], pts[k + 1]);
                ++r_rep.mNumSegments;
            }
        }
        if (r_rep.mNumCoplanarPairs > 0)
            throw MbFailure{"coplanar contact"};

        // Points on an edge go to both triangles beside it.
        const auto share_edges = [&](const MbSurface& rS, const std::vector<std::int64_t>& rIds,
                                     std::vector<MbCut>& rCut) {
            std::unordered_map<std::uint64_t, std::vector<std::int64_t>> on_edge;
            for (std::size_t t = 0; t < rS.mT.size(); ++t) {
                if (rCut[t].mPoints.empty())
                    continue;
                const auto& tri = rS.mT[t];
                const auto [ax, sign] = mb_axis(rS.mP[static_cast<std::size_t>(tri[0])],
                                                rS.mP[static_cast<std::size_t>(tri[1])],
                                                rS.mP[static_cast<std::size_t>(tri[2])]);
                (void)sign;
                for (const std::int64_t q : rCut[t].mPoints)
                    for (int k = 0; k < 3; ++k) {
                        const std::int64_t u = rIds[static_cast<std::size_t>(tri[k])];
                        const std::int64_t v = rIds[static_cast<std::size_t>(tri[(k + 1) % 3])];
                        if (q != u && q != v && mb_orient2d(store[u], store[v], store[q], ax) == 0)
                            on_edge[mb_edge_key(u, v)].push_back(q);
                    }
            }
            for (std::size_t t = 0; t < rS.mT.size(); ++t) {
                const auto& tri = rS.mT[t];
                for (int k = 0; k < 3; ++k) {
                    const auto it =
                        on_edge.find(mb_edge_key(rIds[static_cast<std::size_t>(tri[k])],
                                                 rIds[static_cast<std::size_t>(tri[(k + 1) % 3])]));
                    if (it != on_edge.end())
                        rCut[t].mPoints.insert(rCut[t].mPoints.end(), it->second.begin(),
                                               it->second.end());
                }
            }
        };
        share_edges(sa, ida, cut_a);
        share_edges(sb, idb, cut_b);
        r_rep.mSecondsArrangement = mb_seconds_since(t0);

        // Retriangulation.
        t0 = std::chrono::steady_clock::now();
        std::unordered_set<std::uint64_t> curves;
        const auto pieces = [&](const MbSurface& rS, const std::vector<std::int64_t>& rIds,
                                std::vector<MbCut>& rCut) {
            MbPieces out;
            for (std::size_t t = 0; t < rS.mT.size(); ++t) {
                const auto& tri = rS.mT[t];
                const std::array<std::int64_t, 3> corners{rIds[static_cast<std::size_t>(tri[0])],
                                                          rIds[static_cast<std::size_t>(tri[1])],
                                                          rIds[static_cast<std::size_t>(tri[2])]};
                const MbVec& pa_ = rS.mP[static_cast<std::size_t>(tri[0])];
                const MbVec& pb_ = rS.mP[static_cast<std::size_t>(tri[1])];
                const MbVec& pc_ = rS.mP[static_cast<std::size_t>(tri[2])];
                const MbVec e1{pb_[0] - pa_[0], pb_[1] - pa_[1], pb_[2] - pa_[2]};
                const MbVec e2{pc_[0] - pa_[0], pc_[1] - pa_[1], pc_[2] - pa_[2]};
                const MbVec nrm{e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2],
                                e1[0] * e2[1] - e1[1] * e2[0]};
                if (rCut[t].mPoints.empty()) {
                    out.mT.push_back(corners);
                    out.mNormal.push_back(nrm);
                    continue;
                }
                ++r_rep.mNumCutTriangles;
                const auto [ax, sign] = mb_axis(rS.mP[static_cast<std::size_t>(tri[0])],
                                                rS.mP[static_cast<std::size_t>(tri[1])],
                                                rS.mP[static_cast<std::size_t>(tri[2])]);
                if (sign == 0)
                    throw MbFailure{"a degenerate input triangle is cut"};
                MbTriangulation tr(store, ax, sign, corners);
                std::vector<std::int64_t>& r_pts = rCut[t].mPoints;
                std::sort(r_pts.begin(), r_pts.end());
                r_pts.erase(std::unique(r_pts.begin(), r_pts.end()), r_pts.end());
                for (const std::int64_t q : r_pts)
                    tr.Insert(q);
                for (const auto& [u, v] : rCut[t].mSegments)
                    tr.Constrain(tr.Insert(u), tr.Insert(v));
                tr.Check("retriangulation");
                curves.insert(tr.Constrained().begin(), tr.Constrained().end());
                for (const MbTri& r_sub : tr.Triangles()) {
                    out.mT.push_back(r_sub);
                    out.mNormal.push_back(nrm);
                }
            }
            r_rep.mNumSubTriangles += static_cast<std::int64_t>(out.mT.size());
            return out;
        };
        MbPieces pa = pieces(sa, ida, cut_a);
        MbPieces pb = pieces(sb, idb, cut_b);
        r_rep.mSecondsTriangulate = mb_seconds_since(t0);

        // Classification.
        t0 = std::chrono::steady_clock::now();
        mb_patches(pa, curves);
        mb_patches(pb, curves);
        r_rep.mNumPatchesA = pa.mNumPatches;
        r_rep.mNumPatchesB = pb.mNumPatches;
        const auto classify = [&](const MbPieces& rP, const MbSurface& rOther) {
            // One sample per patch: the centroid of its largest piece.
            std::vector<double> best(static_cast<std::size_t>(rP.mNumPatches), -1.0);
            std::vector<MbVec> sample(static_cast<std::size_t>(rP.mNumPatches));
            for (std::size_t t = 0; t < rP.mT.size(); ++t) {
                const MbVec& a = store[rP.mT[t][0]].mApprox;
                const MbVec& b = store[rP.mT[t][1]].mApprox;
                const MbVec& c = store[rP.mT[t][2]].mApprox;
                const MbVec u{b[0] - a[0], b[1] - a[1], b[2] - a[2]};
                const MbVec v{c[0] - a[0], c[1] - a[1], c[2] - a[2]};
                const double area2 = std::pow(u[1] * v[2] - u[2] * v[1], 2) +
                                     std::pow(u[2] * v[0] - u[0] * v[2], 2) +
                                     std::pow(u[0] * v[1] - u[1] * v[0], 2);
                const auto p = static_cast<std::size_t>(rP.mPatch[t]);
                if (area2 > best[p]) {
                    best[p] = area2;
                    sample[p] = {(a[0] + b[0] + c[0]) / 3.0, (a[1] + b[1] + c[1]) / 3.0,
                                 (a[2] + b[2] + c[2]) / 3.0};
                }
            }
            std::vector<char> inside(sample.size());
            for (std::size_t p = 0; p < sample.size(); ++p)
                inside[p] = mb_inside(sample[p], rOther, lo, hi, r_rep.mNumRayRetries) ? 1 : 0;
            return inside;
        };
        const std::vector<char> in_a = classify(pa, sb);
        const std::vector<char> in_b = classify(pb, sa);
        // Neighbours across a curve should differ (a touching contact need not).
        const auto conflicts = [&](const MbPieces& rP, const std::vector<char>& rIn) {
            std::unordered_map<std::uint64_t, std::vector<std::size_t>> edges;
            for (std::size_t t = 0; t < rP.mT.size(); ++t)
                for (int k = 0; k < 3; ++k) {
                    const std::uint64_t key = mb_edge_key(rP.mT[t][k], rP.mT[t][(k + 1) % 3]);
                    if (curves.count(key))
                        edges[key].push_back(t);
                }
            for (const auto& r_e : edges)
                if (r_e.second.size() == 2 &&
                    rIn[static_cast<std::size_t>(rP.mPatch[r_e.second[0]])] ==
                        rIn[static_cast<std::size_t>(rP.mPatch[r_e.second[1]])])
                    ++r_rep.mNumLabelConflicts;
        };
        conflicts(pa, in_a);
        conflicts(pb, in_b);
        r_rep.mSecondsClassify = mb_seconds_since(t0);

        // Selection.
        std::vector<MbTri> out;
        std::vector<MbVec> parent_normal;
        std::vector<std::int64_t> source;
        const bool keep_a_inside = Op == BooleanOp::Intersection;
        const bool keep_b_inside = Op != BooleanOp::Union;
        for (std::size_t t = 0; t < pa.mT.size(); ++t)
            if ((in_a[static_cast<std::size_t>(pa.mPatch[t])] != 0) == keep_a_inside) {
                out.push_back(pa.mT[t]);
                parent_normal.push_back(pa.mNormal[t]);
                source.push_back(0);
            }
        for (std::size_t t = 0; t < pb.mT.size(); ++t)
            if ((in_b[static_cast<std::size_t>(pb.mPatch[t])] != 0) == keep_b_inside) {
                MbTri tri = pb.mT[t];
                MbVec nrm = pb.mNormal[t];
                if (Op == BooleanOp::Difference) {
                    std::swap(tri[1], tri[2]);
                    nrm = {-nrm[0], -nrm[1], -nrm[2]};
                }
                out.push_back(tri);
                parent_normal.push_back(nrm);
                source.push_back(1);
            }

        // The mesh: the points used, rounded.
        std::unordered_map<std::int64_t, std::int64_t> local;
        std::vector<std::int64_t> order;
        for (const MbTri& r_t : out)
            for (const std::int64_t v : r_t)
                if (local.emplace(v, static_cast<std::int64_t>(order.size())).second)
                    order.push_back(v);
        NDArray p = NDArray::Uninit(DType::Float64, {order.size(), std::size_t{3}});
        for (std::size_t i = 0; i < order.size(); ++i)
            for (std::size_t k = 0; k < 3; ++k)
                p.As<double>()[3 * i + k] = store[order[i]].mApprox[k];
        result.mMesh.AssignPoints(std::move(p));
        NDArray c = NDArray::Uninit(DType::Int64, {out.size(), std::size_t{3}});
        for (std::size_t t = 0; t < out.size(); ++t)
            for (std::size_t k = 0; k < 3; ++k)
                c.As<std::int64_t>()[3 * t + k] = local.at(out[t][k]);
        result.mMesh.AddCellBlock("triangle", std::move(c));
        NDArray src = NDArray::Uninit(DType::Int64, {source.size()});
        std::copy(source.begin(), source.end(), src.As<std::int64_t>());
        std::vector<NDArray> blocks;
        blocks.push_back(std::move(src));
        result.mMesh.AddCellData("boolean:source", std::move(blocks));
        r_rep.mNumOutputTriangles = static_cast<std::int64_t>(out.size());
        // What rounding the points to doubles cost: pieces whose rounded
        // corners are collinear or turn against the exact piece (its parent's
        // normal), though every piece is positive exactly.
        for (std::size_t t = 0; t < out.size(); ++t) {
            const MbVec& a = store[out[t][0]].mApprox;
            const MbVec& b = store[out[t][1]].mApprox;
            const MbVec& c = store[out[t][2]].mApprox;
            const MbVec u{b[0] - a[0], b[1] - a[1], b[2] - a[2]};
            const MbVec v{c[0] - a[0], c[1] - a[1], c[2] - a[2]};
            const MbVec n{u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2],
                          u[0] * v[1] - u[1] * v[0]};
            const MbVec& rn = parent_normal[t];
            if (n[0] * rn[0] + n[1] * rn[1] + n[2] * rn[2] <= 0.0)
                ++r_rep.mNumRoundingInverted;
        }
        result.mOk = true;
    } catch (const MbFailure& rF) {
        result.mOk = false;
        result.mFailure = rF.mWhy;
    }
    return result;
}

}  // namespace detail
}  // namespace meshioplusplus

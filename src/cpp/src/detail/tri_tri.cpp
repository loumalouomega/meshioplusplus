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
// System includes
#include <array>
#include <cmath>
#include <cstddef>

// Project includes
#include "predicates.hpp"
#include "tri_tri.hpp"

namespace meshioplusplus::detail {

namespace {

int ttr_sign(double V) {
    return (V > 0.0) - (V < 0.0);
}

// Two signs strictly opposite.
bool ttr_opposite(int A, int B) {
    return A * B < 0;
}

// In the plane: the two coordinates the axis `Drop` leaves.
std::array<double, 2> ttr_project(const double* pP, int Drop) {
    return Drop == 0   ? std::array<double, 2>{pP[1], pP[2]}
           : Drop == 1 ? std::array<double, 2>{pP[2], pP[0]}
                       : std::array<double, 2>{pP[0], pP[1]};
}

// Whether 2-D segments ab and cd meet (touching included).
bool ttr_segments_2d(const std::array<double, 2>& rA, const std::array<double, 2>& rB,
                     const std::array<double, 2>& rC, const std::array<double, 2>& rD) {
    const int o1 = ttr_sign(orient2d(rA.data(), rB.data(), rC.data()));
    const int o2 = ttr_sign(orient2d(rA.data(), rB.data(), rD.data()));
    const int o3 = ttr_sign(orient2d(rC.data(), rD.data(), rA.data()));
    const int o4 = ttr_sign(orient2d(rC.data(), rD.data(), rB.data()));
    if (o1 == 0 && o2 == 0 && o3 == 0 && o4 == 0) {
        // On one line: their extents along it overlap.
        const int axis = std::fabs(rA[0] - rB[0]) + std::fabs(rC[0] - rD[0]) >=
                                 std::fabs(rA[1] - rB[1]) + std::fabs(rC[1] - rD[1])
                             ? 0
                             : 1;
        const double lo1 = std::fmin(rA[axis], rB[axis]), hi1 = std::fmax(rA[axis], rB[axis]);
        const double lo2 = std::fmin(rC[axis], rD[axis]), hi2 = std::fmax(rC[axis], rD[axis]);
        return lo1 <= hi2 && lo2 <= hi1;
    }
    return o1 * o2 <= 0 && o3 * o4 <= 0;
}

// Whether 2-D point p lies in or on triangle abc (either orientation). A
// degenerate triangle has no inside: its segments are the edges' tests.
bool ttr_inside_2d(const std::array<double, 2>& rP, const std::array<double, 2>& rA,
                   const std::array<double, 2>& rB, const std::array<double, 2>& rC) {
    if (orient2d(rA.data(), rB.data(), rC.data()) == 0.0)
        return false;
    const int s1 = ttr_sign(orient2d(rA.data(), rB.data(), rP.data()));
    const int s2 = ttr_sign(orient2d(rB.data(), rC.data(), rP.data()));
    const int s3 = ttr_sign(orient2d(rC.data(), rA.data(), rP.data()));
    return !(ttr_opposite(s1, s2) || ttr_opposite(s2, s3) || ttr_opposite(s1, s3));
}

bool ttr_coplanar(const double* const* pT1, const double* const* pT2) {
    // Drop the axis the plane's normal is largest along.
    double n[3] = {0.0, 0.0, 0.0};
    const auto accumulate = [&](const double* const* pT) {
        const double u[3] = {pT[1][0] - pT[0][0], pT[1][1] - pT[0][1], pT[1][2] - pT[0][2]};
        const double v[3] = {pT[2][0] - pT[0][0], pT[2][1] - pT[0][1], pT[2][2] - pT[0][2]};
        const double c[3] = {u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2],
                             u[0] * v[1] - u[1] * v[0]};
        for (int k = 0; k < 3; ++k)
            n[k] += std::fabs(c[k]);
    };
    accumulate(pT1);
    accumulate(pT2);
    const int drop = n[0] >= n[1] && n[0] >= n[2] ? 0 : (n[1] >= n[2] ? 1 : 2);
    std::array<std::array<double, 2>, 3> a{}, b{};
    for (int k = 0; k < 3; ++k) {
        a[static_cast<std::size_t>(k)] = ttr_project(pT1[k], drop);
        b[static_cast<std::size_t>(k)] = ttr_project(pT2[k], drop);
    }
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            if (ttr_segments_2d(a[i], a[(i + 1) % 3], b[j], b[(j + 1) % 3]))
                return true;
    return ttr_inside_2d(a[0], b[0], b[1], b[2]) || ttr_inside_2d(b[0], a[0], a[1], a[2]);
}

bool ttr_same(const double* pA, const double* pB) {
    return pA[0] == pB[0] && pA[1] == pB[1] && pA[2] == pB[2];
}

// The axis the normal of triangle pqr is largest along.
int ttr_drop_axis(const double* pP, const double* pQ, const double* pR) {
    const double u[3] = {pQ[0] - pP[0], pQ[1] - pP[1], pQ[2] - pP[2]};
    const double v[3] = {pR[0] - pP[0], pR[1] - pP[1], pR[2] - pP[2]};
    const double n[3] = {std::fabs(u[1] * v[2] - u[2] * v[1]), std::fabs(u[2] * v[0] - u[0] * v[2]),
                         std::fabs(u[0] * v[1] - u[1] * v[0])};
    return n[0] >= n[1] && n[0] >= n[2] ? 0 : (n[1] >= n[2] ? 1 : 2);
}

// Whether the edge from corner v towards a enters triangle (v, c, d) beyond v:
// a lies in the triangle's plane and the direction v -> a in its corner at v.
bool ttr_enters_corner(const double* pV, const double* pA, const double* pC, const double* pD) {
    if (orient3d(pV, pC, pD, pA) != 0.0)
        return false;
    const int drop = ttr_drop_axis(pV, pC, pD);
    const auto v = ttr_project(pV, drop), a = ttr_project(pA, drop);
    const auto c = ttr_project(pC, drop), d = ttr_project(pD, drop);
    const int o = ttr_sign(orient2d(v.data(), c.data(), d.data()));
    if (o == 0)
        return false;
    return ttr_sign(orient2d(v.data(), c.data(), a.data())) * o >= 0 &&
           ttr_sign(orient2d(v.data(), a.data(), d.data())) * o >= 0;
}

// The points where the edges of triangle t cross the plane whose orient3d
// values at t's corners are rD (the corners in it included).
int ttr_plane_points(const double* const* pT, const double* pD, double (*pOut)[3]) {
    int n = 0;
    for (int i = 0; i < 3 && n < 3; ++i) {
        if (pD[i] == 0.0) {
            for (int k = 0; k < 3; ++k)
                pOut[n][k] = pT[i][k];
            ++n;
        }
        const int j = (i + 1) % 3;
        if ((pD[i] > 0.0 && pD[j] < 0.0) || (pD[i] < 0.0 && pD[j] > 0.0)) {
            const double w = pD[i] / (pD[i] - pD[j]);
            for (int k = 0; k < 3; ++k)
                pOut[n][k] = pT[i][k] + w * (pT[j][k] - pT[i][k]);
            ++n;
        }
    }
    return n;
}

}  // namespace

bool segment_tri_intersect(const double* pA, const double* pB, const double* pP, const double* pQ,
                           const double* pR) {
    const int sa = ttr_sign(orient3d(pP, pQ, pR, pA));
    const int sb = ttr_sign(orient3d(pP, pQ, pR, pB));
    if (sa != 0 && sa == sb)
        return false;
    if (sa == 0 && sb == 0) {
        // In the triangle's plane: as a degenerate triangle.
        const double* t1[3] = {pA, pB, pB};
        const double* t2[3] = {pP, pQ, pR};
        return ttr_coplanar(t1, t2);
    }
    const int s1 = ttr_sign(orient3d(pA, pB, pP, pQ));
    const int s2 = ttr_sign(orient3d(pA, pB, pQ, pR));
    const int s3 = ttr_sign(orient3d(pA, pB, pR, pP));
    return !(ttr_opposite(s1, s2) || ttr_opposite(s2, s3) || ttr_opposite(s1, s3));
}

bool tri_tri_intersect(const double* pA, const double* pB, const double* pC, const double* pD,
                       const double* pE, const double* pF) {
    const int d1 = ttr_sign(orient3d(pD, pE, pF, pA));
    const int d2 = ttr_sign(orient3d(pD, pE, pF, pB));
    const int d3 = ttr_sign(orient3d(pD, pE, pF, pC));
    if (d1 != 0 && d1 == d2 && d2 == d3)
        return false;
    const int e1 = ttr_sign(orient3d(pA, pB, pC, pD));
    const int e2 = ttr_sign(orient3d(pA, pB, pC, pE));
    const int e3 = ttr_sign(orient3d(pA, pB, pC, pF));
    if (e1 != 0 && e1 == e2 && e2 == e3)
        return false;
    const double* t1[3] = {pA, pB, pC};
    const double* t2[3] = {pD, pE, pF};
    if (d1 == 0 && d2 == 0 && d3 == 0)
        return ttr_coplanar(t1, t2);
    for (int i = 0; i < 3; ++i)
        if (segment_tri_intersect(t1[i], t1[(i + 1) % 3], pD, pE, pF) ||
            segment_tri_intersect(t2[i], t2[(i + 1) % 3], pA, pB, pC))
            return true;
    return false;
}

bool tri_tri_intersect_unshared(const double* pA, const double* pB, const double* pC,
                                const double* pD, const double* pE, const double* pF) {
    const double* t1[3] = {pA, pB, pC};
    const double* t2[3] = {pD, pE, pF};
    // Which corner of t2 each corner of t1 is, or -1.
    int match[3] = {-1, -1, -1};
    int shared = 0;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            if (match[i] < 0 && ttr_same(t1[i], t2[j])) {
                match[i] = j;
                ++shared;
                break;
            }
    if (shared == 0)
        return tri_tri_intersect(pA, pB, pC, pD, pE, pF);
    if (shared >= 3)
        return true;
    if (shared == 2) {
        // On one edge pq: they overlap beyond it only folded flat onto each
        // other, their third corners on one side of it in their common plane.
        int r = 0;
        while (match[r] >= 0)
            ++r;
        const double* p = t1[(r + 1) % 3];
        const double* q = t1[(r + 2) % 3];
        int s = 0;
        while (s < 3 && (ttr_same(t2[s], p) || ttr_same(t2[s], q)))
            ++s;
        if (s == 3)
            return false;
        if (orient3d(p, q, t1[r], t2[s]) != 0.0)
            return false;
        const int drop = ttr_drop_axis(p, q, t1[r]);
        const auto p2 = ttr_project(p, drop), q2 = ttr_project(q, drop);
        const auto r2 = ttr_project(t1[r], drop), s2 = ttr_project(t2[s], drop);
        const int sr = ttr_sign(orient2d(p2.data(), q2.data(), r2.data()));
        const int ss = ttr_sign(orient2d(p2.data(), q2.data(), s2.data()));
        return sr != 0 && sr == ss;
    }
    // On one corner v: beyond it, an opposite edge crosses the other triangle,
    // or an edge from v runs into the other's corner at v.
    int i = 0;
    while (match[i] < 0)
        ++i;
    const int j = match[i];
    const double* v = t1[i];
    const double* a = t1[(i + 1) % 3];
    const double* b = t1[(i + 2) % 3];
    const double* c = t2[(j + 1) % 3];
    const double* d = t2[(j + 2) % 3];
    return segment_tri_intersect(a, b, pD, pE, pF) || segment_tri_intersect(c, d, pA, pB, pC) ||
           ttr_enters_corner(v, a, c, d) || ttr_enters_corner(v, b, c, d) ||
           ttr_enters_corner(v, c, a, b) || ttr_enters_corner(v, d, a, b);
}

TriTriSegment tri_tri_segment(const double* pA, const double* pB, const double* pC,
                              const double* pD, const double* pE, const double* pF) {
    TriTriSegment out;
    if (!tri_tri_intersect(pA, pB, pC, pD, pE, pF))
        return out;
    const double* t1[3] = {pA, pB, pC};
    const double* t2[3] = {pD, pE, pF};
    const double d1[3] = {orient3d(pD, pE, pF, pA), orient3d(pD, pE, pF, pB),
                          orient3d(pD, pE, pF, pC)};
    if (d1[0] == 0.0 && d1[1] == 0.0 && d1[2] == 0.0) {
        out.mKind = TriTriContact::Coplanar;
        return out;
    }
    const double d2[3] = {orient3d(pA, pB, pC, pD), orient3d(pA, pB, pC, pE),
                          orient3d(pA, pB, pC, pF)};
    double p1[3][3], p2[3][3];
    const int n1 = ttr_plane_points(t1, d1, p1);
    const int n2 = ttr_plane_points(t2, d2, p2);
    // Both sets lie on the line the planes share: overlap their extents along it.
    const auto normal = [](const double* const* pT, double* pN) {
        const double u[3] = {pT[1][0] - pT[0][0], pT[1][1] - pT[0][1], pT[1][2] - pT[0][2]};
        const double w[3] = {pT[2][0] - pT[0][0], pT[2][1] - pT[0][1], pT[2][2] - pT[0][2]};
        pN[0] = u[1] * w[2] - u[2] * w[1];
        pN[1] = u[2] * w[0] - u[0] * w[2];
        pN[2] = u[0] * w[1] - u[1] * w[0];
    };
    double na[3], nb[3];
    normal(t1, na);
    normal(t2, nb);
    double dir[3] = {na[1] * nb[2] - na[2] * nb[1], na[2] * nb[0] - na[0] * nb[2],
                     na[0] * nb[1] - na[1] * nb[0]};
    if (dir[0] == 0.0 && dir[1] == 0.0 && dir[2] == 0.0) {
        // Planes parallel in double although not exactly: take the spread of
        // the points found.
        const double(*src)[3] = n1 >= 2 ? p1 : p2;
        for (int k = 0; k < 3; ++k)
            dir[k] = src[(n1 >= 2 ? n1 : n2) - 1][k] - src[0][k];
    }
    const auto along = [&](const double* pX) {
        return pX[0] * dir[0] + pX[1] * dir[1] + pX[2] * dir[2];
    };
    const auto extent = [&](const double(*pPts)[3], int N, int& rLo, int& rHi) {
        rLo = rHi = 0;
        for (int k = 1; k < N; ++k) {
            if (along(pPts[k]) < along(pPts[rLo]))
                rLo = k;
            if (along(pPts[k]) > along(pPts[rHi]))
                rHi = k;
        }
    };
    if (n1 == 0 || n2 == 0) {
        // Rounding lost the crossing the exact test saw: report a point.
        const double* x = n1 > 0 ? p1[0] : (n2 > 0 ? p2[0] : pA);
        out.mKind = TriTriContact::Point;
        for (int k = 0; k < 3; ++k)
            out.mP[k] = out.mQ[k] = x[k];
        return out;
    }
    int lo1, hi1, lo2, hi2;
    extent(p1, n1, lo1, hi1);
    extent(p2, n2, lo2, hi2);
    const double* lo = along(p1[lo1]) >= along(p2[lo2]) ? p1[lo1] : p2[lo2];
    const double* hi = along(p1[hi1]) <= along(p2[hi2]) ? p1[hi1] : p2[hi2];
    if (along(lo) > along(hi)) {
        // Rounding inverted a touching contact: report its middle.
        for (int k = 0; k < 3; ++k)
            out.mP[k] = out.mQ[k] = 0.5 * (lo[k] + hi[k]);
        out.mKind = TriTriContact::Point;
        return out;
    }
    for (int k = 0; k < 3; ++k) {
        out.mP[k] = lo[k];
        out.mQ[k] = hi[k];
    }
    out.mKind = ttr_same(out.mP, out.mQ) ? TriTriContact::Point : TriTriContact::Segment;
    return out;
}

}  // namespace meshioplusplus::detail

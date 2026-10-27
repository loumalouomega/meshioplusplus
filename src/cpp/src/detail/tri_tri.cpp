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

}  // namespace meshioplusplus::detail

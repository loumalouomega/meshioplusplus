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
//                   SPDX-License-Identifier: MIT
//                   meshio++ default license: LICENSE
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//
// Moved out of formats/frd.cpp (v15.3.0) so operations/tensor_invariants.cpp
// can share the same eigensolver instead of duplicating it.

// System includes
#include <algorithm>
#include <cmath>

// Project includes
#include "meshioplusplus/detail/sym3_eigen.hpp"

namespace meshioplusplus {
namespace detail {

double sym3_mises(const double* pT) {
    const double xx = pT[0], yy = pT[1], zz = pT[2], xy = pT[3], yz = pT[4], xz = pT[5];
    return std::sqrt(0.5 * ((xx - yy) * (xx - yy) + (yy - zz) * (yy - zz) + (zz - xx) * (zz - xx) +
                            6.0 * (xy * xy + yz * yz + xz * xz)));
}

namespace {

/// Cyclic Jacobi on the symmetric matrix of @p pT: leaves the (diagonalised)
/// matrix in @p a and, when @p pV is non-null, accumulates the rotations into
/// it (the columns of @p pV are the eigenvectors, unsorted). The matrix update
/// does not depend on @p pV, so `sym3_principal`'s values are the same with or
/// without it.
void sym3_jacobi(const double* pT, double a[3][3], double (*pV)[3]) {
    a[0][0] = pT[0];
    a[0][1] = a[1][0] = pT[3];
    a[0][2] = a[2][0] = pT[5];
    a[1][1] = pT[1];
    a[1][2] = a[2][1] = pT[4];
    a[2][2] = pT[2];
    if (pV != nullptr)
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                pV[i][j] = i == j ? 1.0 : 0.0;
    for (int sweep = 0; sweep < 60; ++sweep) {
        const double off = std::fabs(a[0][1]) + std::fabs(a[0][2]) + std::fabs(a[1][2]);
        const double diag = std::fabs(a[0][0]) + std::fabs(a[1][1]) + std::fabs(a[2][2]);
        if (off <= 1e-17 * diag || off == 0.0)
            break;
        for (int p = 0; p < 2; ++p) {
            for (int q = p + 1; q < 3; ++q) {
                if (a[p][q] == 0.0)
                    continue;
                const double theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                const double t = (theta >= 0.0 ? 1.0 : -1.0) /
                                 (std::fabs(theta) + std::sqrt(theta * theta + 1.0));
                const double c = 1.0 / std::sqrt(t * t + 1.0);
                const double s = t * c;
                for (int k = 0; k < 3; ++k) {
                    const double akp = a[k][p];
                    const double akq = a[k][q];
                    a[k][p] = c * akp - s * akq;
                    a[k][q] = s * akp + c * akq;
                }
                for (int k = 0; k < 3; ++k) {
                    const double apk = a[p][k];
                    const double aqk = a[q][k];
                    a[p][k] = c * apk - s * aqk;
                    a[q][k] = s * apk + c * aqk;
                }
                if (pV != nullptr) {
                    for (int k = 0; k < 3; ++k) {
                        const double vkp = pV[k][p];
                        const double vkq = pV[k][q];
                        pV[k][p] = c * vkp - s * vkq;
                        pV[k][q] = s * vkp + c * vkq;
                    }
                }
            }
        }
    }
}

}  // namespace

void sym3_principal(const double* pT, double* pOut) {
    double a[3][3];
    sym3_jacobi(pT, a, nullptr);
    pOut[0] = a[0][0];
    pOut[1] = a[1][1];
    pOut[2] = a[2][2];
    std::sort(pOut, pOut + 3);
}

void sym3_eigen(const double* pT, double* pValues, double* pVectors) {
    double a[3][3];
    double v[3][3];
    sym3_jacobi(pT, a, v);
    int order[3] = {0, 1, 2};
    // Ascending by value; the index breaks a tie, so equal eigenvalues keep a
    // fixed order.
    std::sort(order, order + 3,
              [&](int x, int y) { return a[x][x] < a[y][y] || (a[x][x] == a[y][y] && x < y); });
    for (int k = 0; k < 3; ++k) {
        pValues[k] = a[order[k]][order[k]];
        double col[3] = {v[0][order[k]], v[1][order[k]], v[2][order[k]]};
        // The sign of an eigenvector is arbitrary: make its first component of
        // magnitude above 1e-8 positive, so the result does not depend on how
        // the rotations happened to land.
        for (int i = 0; i < 3; ++i) {
            if (std::fabs(col[i]) > 1e-8) {
                if (col[i] < 0.0)
                    for (int j = 0; j < 3; ++j)
                        col[j] = -col[j];
                break;
            }
        }
        for (int i = 0; i < 3; ++i)
            pVectors[i * 3 + k] = col[i];
    }
}

}  // namespace detail
}  // namespace meshioplusplus

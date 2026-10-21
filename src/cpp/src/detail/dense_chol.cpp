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
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

// Project includes (private, not installed)
#include "dense_chol.hpp"

namespace meshioplusplus {
namespace detail {

bool chol_factor(std::vector<double>& rA, std::size_t N) {
    for (std::size_t j = 0; j < N; ++j) {
        double d = rA[j * N + j];
        for (std::size_t k = 0; k < j; ++k)
            d -= rA[j * N + k] * rA[j * N + k];
        if (!(d > 0.0) || !std::isfinite(d))
            return false;
        const double ljj = std::sqrt(d);
        rA[j * N + j] = ljj;
        for (std::size_t i = j + 1; i < N; ++i) {
            double s = rA[i * N + j];
            for (std::size_t k = 0; k < j; ++k)
                s -= rA[i * N + k] * rA[j * N + k];
            rA[i * N + j] = s / ljj;
        }
    }
    return true;
}

void chol_solve(const std::vector<double>& rL, std::size_t N, double* pB) {
    for (std::size_t i = 0; i < N; ++i) {
        double s = pB[i];
        for (std::size_t k = 0; k < i; ++k)
            s -= rL[i * N + k] * pB[k];
        pB[i] = s / rL[i * N + i];
    }
    for (std::size_t ii = N; ii-- > 0;) {
        double s = pB[ii];
        for (std::size_t k = ii + 1; k < N; ++k)
            s -= rL[k * N + ii] * pB[k];
        pB[ii] = s / rL[ii * N + ii];
    }
}

bool ldlt_factor(std::vector<double>& rA, std::size_t N, double Shift) {
    std::vector<double> g(N);
    for (std::size_t i = 0; i < N; ++i) {
        // g[j] = L_ij D_j for j < i, built row by row.
        for (std::size_t j = 0; j < i; ++j) {
            double s = rA[i * N + j];
            for (std::size_t k = 0; k < j; ++k)
                s -= g[k] * rA[j * N + k];
            g[j] = s;
        }
        double d = rA[i * N + i] + Shift;
        for (std::size_t j = 0; j < i; ++j) {
            const double lij = g[j] / rA[j * N + j];
            d -= g[j] * lij;
            rA[i * N + j] = lij;
        }
        if (d == 0.0 || !std::isfinite(d))
            return false;
        rA[i * N + i] = d;
    }
    return true;
}

void ldlt_solve(const std::vector<double>& rLd, std::size_t N, double* pB) {
    for (std::size_t i = 0; i < N; ++i) {
        double s = pB[i];
        for (std::size_t k = 0; k < i; ++k)
            s -= rLd[i * N + k] * pB[k];
        pB[i] = s;
    }
    for (std::size_t i = 0; i < N; ++i)
        pB[i] /= rLd[i * N + i];
    for (std::size_t ii = N; ii-- > 0;) {
        double s = pB[ii];
        for (std::size_t k = ii + 1; k < N; ++k)
            s -= rLd[k * N + ii] * pB[k];
        pB[ii] = s;
    }
}

SkylineMatrix skyline_zero(const std::vector<std::size_t>& rFirst) {
    SkylineMatrix a;
    a.mN = rFirst.size();
    a.mFirst.resize(a.mN);
    a.mStart.assign(a.mN + 1, 0);
    for (std::size_t i = 0; i < a.mN; ++i) {
        a.mFirst[i] = std::min(rFirst[i], i);
        a.mStart[i + 1] = a.mStart[i] + (i - a.mFirst[i] + 1);
    }
    a.mValues.assign(a.mStart[a.mN], 0.0);
    return a;
}

bool skyline_ldlt(SkylineMatrix& rA, double Shift) {
    const std::size_t n = rA.mN;
    std::vector<double> g(n);
    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t fi = rA.mFirst[i];
        double* row_i = rA.mValues.data() + rA.mStart[i];
        // g[j] = L_ij D_j for fi <= j < i; L_ik = 0 for k < fi.
        for (std::size_t j = fi; j < i; ++j) {
            const std::size_t fj = rA.mFirst[j];
            const double* row_j = rA.mValues.data() + rA.mStart[j];
            double s = row_i[j - fi];
            for (std::size_t k = std::max(fi, fj); k < j; ++k)
                s -= g[k] * row_j[k - fj];
            g[j] = s;
        }
        double d = row_i[i - fi] + Shift;
        for (std::size_t j = fi; j < i; ++j) {
            const double dj = rA.mValues[rA.mStart[j + 1] - 1];
            const double lij = g[j] / dj;
            d -= g[j] * lij;
            row_i[j - fi] = lij;
        }
        if (!(d > 0.0) || !std::isfinite(d))
            return false;
        row_i[i - fi] = d;
    }
    return true;
}

void skyline_solve(const SkylineMatrix& rLd, double* pB) {
    const std::size_t n = rLd.mN;
    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t fi = rLd.mFirst[i];
        const double* row = rLd.mValues.data() + rLd.mStart[i];
        double s = pB[i];
        for (std::size_t k = fi; k < i; ++k)
            s -= row[k - fi] * pB[k];
        pB[i] = s;
    }
    for (std::size_t i = 0; i < n; ++i)
        pB[i] /= rLd.mValues[rLd.mStart[i + 1] - 1];
    // Back substitution by columns of Lᵀ: row i of L scatters into earlier unknowns.
    for (std::size_t ii = n; ii-- > 0;) {
        const std::size_t fi = rLd.mFirst[ii];
        const double* row = rLd.mValues.data() + rLd.mStart[ii];
        const double xi = pB[ii];
        for (std::size_t k = fi; k < ii; ++k)
            pB[k] -= row[k - fi] * xi;
    }
}

}  // namespace detail
}  // namespace meshioplusplus

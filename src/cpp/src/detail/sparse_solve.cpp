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
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/parallel.hpp"

// Project includes (private, not installed)
#include "sparse_solve.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

// Vector loops this short run serially: a parallel region costs more than
// they do, and far more on a loaded machine.
constexpr std::size_t kSpsGrain = std::size_t{1} << 15;

// y = a x + y, element-wise in parallel (independent entries).
void sps_axpy(double A, const double* pX, double* pY, std::size_t N) {
    parallel_for(N, [&](std::size_t i) { pY[i] += A * pX[i]; }, kSpsGrain);
}

// The relative residual ||b - A x|| / ||b|| (0 for a zero b), into rR.
double sps_residual(const SparseCsr& rA, const double* pB, const double* pX,
                    std::vector<double>& rR, double BNorm) {
    const std::size_t n = rA.mRows;
    csr_multiply(rA, pX, rR.data());
    parallel_for(n, [&](std::size_t i) { rR[i] = pB[i] - rR[i]; }, kSpsGrain);
    const double rn = std::sqrt(sparse_dot(rR.data(), rR.data(), n));
    return BNorm > 0.0 ? rn / BNorm : rn;
}

// The inverse diagonal, 1 where the diagonal is zero or not finite.
std::vector<double> sps_jacobi(const SparseCsr& rA) {
    std::vector<double> d = csr_diagonal(rA);
    for (double& v : d)
        v = (v != 0.0 && std::isfinite(v)) ? 1.0 / v : 1.0;
    return d;
}

void sps_check_square(const SparseCsr& rA, const char* pWho) {
    if (rA.mRows != rA.mCols)
        throw std::invalid_argument(std::string("meshio++: ") + pWho +
                                    ": the matrix is not square");
}

}  // namespace

SparseCsr csr_from_triplets(std::size_t Rows, std::size_t Cols,
                            std::vector<SparseTriplet> rTriplets) {
    for (const SparseTriplet& t : rTriplets)
        if (t.mRow >= Rows || t.mCol >= Cols)
            throw std::invalid_argument(
                "meshio++: csr_from_triplets: an entry lies outside the matrix");
    std::stable_sort(rTriplets.begin(), rTriplets.end(),
                     [](const SparseTriplet& rL, const SparseTriplet& rR) {
                         return rL.mRow != rR.mRow ? rL.mRow < rR.mRow : rL.mCol < rR.mCol;
                     });
    SparseCsr a;
    a.mRows = Rows;
    a.mCols = Cols;
    a.mStart.assign(Rows + 1, 0);
    a.mCol.reserve(rTriplets.size());
    a.mValue.reserve(rTriplets.size());
    for (std::size_t k = 0; k < rTriplets.size(); ++k) {
        const SparseTriplet& t = rTriplets[k];
        if (k && t.mRow == rTriplets[k - 1].mRow && t.mCol == rTriplets[k - 1].mCol) {
            a.mValue.back() += t.mValue;
            continue;
        }
        a.mCol.push_back(t.mCol);
        a.mValue.push_back(t.mValue);
        ++a.mStart[t.mRow + 1];
    }
    for (std::size_t r = 0; r < Rows; ++r)
        a.mStart[r + 1] += a.mStart[r];
    return a;
}

void csr_multiply(const SparseCsr& rA, const double* pX, double* pY) {
    parallel_for(
        rA.mRows,
        [&](std::size_t r) {
            double s = 0.0;
            for (std::size_t k = rA.mStart[r]; k < rA.mStart[r + 1]; ++k)
                s += rA.mValue[k] * pX[rA.mCol[k]];
            pY[r] = s;
        },
        kSpsGrain);
}

std::vector<double> csr_diagonal(const SparseCsr& rA) {
    std::vector<double> d(std::min(rA.mRows, rA.mCols), 0.0);
    for (std::size_t r = 0; r < d.size(); ++r)
        for (std::size_t k = rA.mStart[r]; k < rA.mStart[r + 1]; ++k)
            if (rA.mCol[k] == r)
                d[r] = rA.mValue[k];
    return d;
}

double sparse_dot(const double* pA, const double* pB, std::size_t N) {
    const auto chunk = [&](std::size_t Begin, std::size_t End) {
        double s = 0.0;
        for (std::size_t i = Begin; i < End; ++i)
            s += pA[i] * pB[i];
        return s;
    };
    if (N <= kSpsGrain) {
        // The same chunks, summed the same way, without a parallel region.
        double acc = 0.0;
        for (std::size_t b = 0; b < N; b += kSparseChunk)
            acc = acc + chunk(b, std::min(N, b + kSparseChunk));
        return acc;
    }
    return parallel_reduce(N, kSparseChunk, 0.0, chunk, [](double L, double R) { return L + R; });
}

SparseSolveReport csr_pcg(const SparseCsr& rA, const double* pB, double* pX,
                          double RelativeTolerance, std::size_t MaxIterations) {
    sps_check_square(rA, "csr_pcg");
    const std::size_t n = rA.mRows;
    SparseSolveReport rep;
    if (n == 0) {
        rep.mConverged = true;
        return rep;
    }
    if (MaxIterations == 0)
        MaxIterations = 2 * n + 10;
    const std::vector<double> minv = sps_jacobi(rA);
    const double bnorm = std::sqrt(sparse_dot(pB, pB, n));
    std::vector<double> r(n), z(n), p(n), ap(n);
    rep.mRelativeResidual = sps_residual(rA, pB, pX, r, bnorm);
    if (bnorm == 0.0) {
        // b = 0: the solution of an SPD system is x = 0.
        std::fill(pX, pX + n, 0.0);
        rep.mConverged = true;
        rep.mRelativeResidual = 0.0;
        return rep;
    }
    if (rep.mRelativeResidual <= RelativeTolerance) {
        rep.mConverged = true;
        return rep;
    }
    parallel_for(n, [&](std::size_t i) { z[i] = minv[i] * r[i]; }, kSpsGrain);
    p = z;
    double rz = sparse_dot(r.data(), z.data(), n);
    for (std::size_t it = 1; it <= MaxIterations; ++it) {
        csr_multiply(rA, p.data(), ap.data());
        const double pap = sparse_dot(p.data(), ap.data(), n);
        if (!(pap > 0.0) || !std::isfinite(pap))
            break;
        const double alpha = rz / pap;
        sps_axpy(alpha, p.data(), pX, n);
        sps_axpy(-alpha, ap.data(), r.data(), n);
        rep.mIterations = it;
        rep.mRelativeResidual = std::sqrt(sparse_dot(r.data(), r.data(), n)) / bnorm;
        if (rep.mRelativeResidual <= RelativeTolerance) {
            rep.mConverged = true;
            break;
        }
        parallel_for(n, [&](std::size_t i) { z[i] = minv[i] * r[i]; }, kSpsGrain);
        const double rz_new = sparse_dot(r.data(), z.data(), n);
        const double beta = rz_new / rz;
        rz = rz_new;
        parallel_for(n, [&](std::size_t i) { p[i] = z[i] + beta * p[i]; }, kSpsGrain);
    }
    // The recurrence's residual drifts from the true one; report the true one.
    rep.mRelativeResidual = sps_residual(rA, pB, pX, r, bnorm);
    rep.mConverged = rep.mRelativeResidual <= 10.0 * RelativeTolerance;
    return rep;
}

SparseSolveReport csr_bicgstab(const SparseCsr& rA, const double* pB, double* pX,
                               double RelativeTolerance, std::size_t MaxIterations) {
    sps_check_square(rA, "csr_bicgstab");
    const std::size_t n = rA.mRows;
    SparseSolveReport rep;
    if (n == 0) {
        rep.mConverged = true;
        return rep;
    }
    if (MaxIterations == 0)
        MaxIterations = 2 * n + 10;
    const std::vector<double> minv = sps_jacobi(rA);
    const double bnorm = std::sqrt(sparse_dot(pB, pB, n));
    std::vector<double> r(n), r0(n), p(n, 0.0), v(n, 0.0), s(n), t(n), ph(n), sh(n);
    rep.mRelativeResidual = sps_residual(rA, pB, pX, r, bnorm);
    if (rep.mRelativeResidual <= RelativeTolerance) {
        rep.mConverged = true;
        return rep;
    }
    r0 = r;
    double rho = 1.0, alpha = 1.0, omega = 1.0;
    std::vector<double> best(pX, pX + n);
    double best_res = rep.mRelativeResidual;
    for (std::size_t it = 1; it <= MaxIterations; ++it) {
        const double rho_new = sparse_dot(r0.data(), r.data(), n);
        if (rho_new == 0.0 || !std::isfinite(rho_new) || omega == 0.0)
            break;
        const double beta = (rho_new / rho) * (alpha / omega);
        rho = rho_new;
        parallel_for(
            n, [&](std::size_t i) { p[i] = r[i] + beta * (p[i] - omega * v[i]); }, kSpsGrain);
        parallel_for(n, [&](std::size_t i) { ph[i] = minv[i] * p[i]; }, kSpsGrain);
        csr_multiply(rA, ph.data(), v.data());
        const double r0v = sparse_dot(r0.data(), v.data(), n);
        if (r0v == 0.0 || !std::isfinite(r0v))
            break;
        alpha = rho / r0v;
        parallel_for(n, [&](std::size_t i) { s[i] = r[i] - alpha * v[i]; }, kSpsGrain);
        rep.mIterations = it;
        const double snorm = std::sqrt(sparse_dot(s.data(), s.data(), n));
        if (snorm / bnorm <= RelativeTolerance) {
            sps_axpy(alpha, ph.data(), pX, n);
            break;
        }
        parallel_for(n, [&](std::size_t i) { sh[i] = minv[i] * s[i]; }, kSpsGrain);
        csr_multiply(rA, sh.data(), t.data());
        const double tt = sparse_dot(t.data(), t.data(), n);
        if (!(tt > 0.0) || !std::isfinite(tt))
            break;
        omega = sparse_dot(t.data(), s.data(), n) / tt;
        parallel_for(n, [&](std::size_t i) { pX[i] += alpha * ph[i] + omega * sh[i]; }, kSpsGrain);
        parallel_for(n, [&](std::size_t i) { r[i] = s[i] - omega * t[i]; }, kSpsGrain);
        const double res = std::sqrt(sparse_dot(r.data(), r.data(), n)) / bnorm;
        if (res < best_res) {
            best_res = res;
            std::copy(pX, pX + n, best.begin());
        }
        if (res <= RelativeTolerance)
            break;
    }
    std::vector<double> tmp(n);
    double res = sps_residual(rA, pB, pX, tmp, bnorm);
    if (best_res < res) {
        std::copy(best.begin(), best.end(), pX);
        res = sps_residual(rA, pB, pX, tmp, bnorm);
    }
    rep.mRelativeResidual = res;
    rep.mConverged = res <= 10.0 * RelativeTolerance;
    return rep;
}

}  // namespace detail
}  // namespace meshioplusplus

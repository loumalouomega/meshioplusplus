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
#include <cstdint>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/eigenmodes.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/parallel.hpp"

// Project includes (private, not installed)
#include "../detail/jacobi_eigen.hpp"
#include "../detail/sparse_solve.hpp"
#include "../detail/surface_fem.hpp"

namespace meshioplusplus {

namespace {

using LapColumn = std::vector<double>;

/// A deterministic value in [-1, 1) for (row, column): splitmix64.
double lap_hash(std::uint64_t Row, std::uint64_t Col) {
    std::uint64_t z = Row * 0x9E3779B97F4A7C15ull + Col * 0xD1B54A32D192ED03ull + 1;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    z ^= z >> 31;
    return static_cast<double>(z >> 11) * (2.0 / 9007199254740992.0) - 1.0;
}

/// `sum a_i M_i b_i`, the M-inner product (M diagonal), in fixed chunks.
double lap_mdot(const LapColumn& rA, const LapColumn& rB, const std::vector<double>& rMass) {
    const std::size_t n = rA.size();
    const std::size_t chunks = (n + detail::kSparseChunk - 1) / detail::kSparseChunk;
    std::vector<double> part(chunks, 0.0);
    parallel_for(chunks, [&](std::size_t c) {
        const std::size_t lo = c * detail::kSparseChunk;
        const std::size_t hi = std::min(n, lo + detail::kSparseChunk);
        double s = 0.0;
        for (std::size_t i = lo; i < hi; ++i)
            s += rA[i] * rMass[i] * rB[i];
        part[c] = s;
    });
    double s = 0.0;
    for (const double v : part)
        s += v;
    return s;
}

/// Eigen-decomposition of a symmetric P x P matrix (row-major) by cyclic
/// Jacobi: eigenvalues in `rW`, eigenvectors as the columns of `rV`.
void lap_jacobi(std::vector<double> A, std::size_t P, std::vector<double>& rW,
                std::vector<double>& rV) {
    rV.resize(P * P);
    rW.resize(P);
    detail::jacobi_eigen_sym(A.data(), P, rW.data(), rV.data(), 1e-32);
}

/// M-orthonormalise the block in place (modified Gram-Schmidt, twice); a
/// column that collapses is replaced by a fresh deterministic one.
void lap_orthonormalise(std::vector<LapColumn>& rY, const std::vector<double>& rMass,
                        std::uint64_t& rSeed) {
    const std::size_t n = rMass.size();
    for (std::size_t j = 0; j < rY.size(); ++j) {
        for (int attempt = 0;; ++attempt) {
            const double before = std::sqrt(lap_mdot(rY[j], rY[j], rMass));
            for (int pass = 0; pass < 2; ++pass)
                for (std::size_t i = 0; i < j; ++i) {
                    const double c = lap_mdot(rY[i], rY[j], rMass);
                    for (std::size_t r = 0; r < n; ++r)
                        rY[j][r] -= c * rY[i][r];
                }
            const double norm = std::sqrt(lap_mdot(rY[j], rY[j], rMass));
            if (norm > 1e-10 * before && norm > 0.0) {
                for (double& v : rY[j])
                    v /= norm;
                break;
            }
            if (attempt > 8)
                throw std::runtime_error(
                    "laplacian_eigenmodes: internal error: the block lost its rank");
            ++rSeed;
            for (std::size_t r = 0; r < n; ++r)
                rY[j][r] = lap_hash(r, rSeed + 1000003ull);
        }
    }
}

}  // namespace

EigenmodeResult laplacian_eigenmodes(const Mesh& rMesh, const EigenmodeOptions& rOptions) {
    if (!(rOptions.mTolerance > 0.0))
        throw std::invalid_argument("laplacian_eigenmodes: the tolerance must be positive");
    if (rOptions.mNumModes < 1)
        throw std::invalid_argument("laplacian_eigenmodes: num_modes must be at least 1");
    const detail::SurfaceFem fem = detail::build_surface_fem(rMesh, "laplacian_eigenmodes");
    const std::size_t n_all = rMesh.NumPoints();
    std::vector<double> mass_all;
    std::vector<detail::SparseTriplet> stiff_all;
    detail::surface_fem_assemble(fem, mass_all, stiff_all);

    // The problem lives on the surface's points alone.
    std::vector<std::int64_t> local(n_all, -1);
    std::vector<std::size_t> global;
    for (std::size_t i = 0; i < n_all; ++i)
        if (mass_all[i] > 0.0) {
            local[i] = static_cast<std::int64_t>(global.size());
            global.push_back(i);
        }
    const std::size_t n = global.size();
    const std::size_t k = static_cast<std::size_t>(rOptions.mNumModes);
    if (k > n)
        throw std::invalid_argument("laplacian_eigenmodes: num_modes (" + std::to_string(k) +
                                    ") exceeds the surface's " + std::to_string(n) + " points");
    std::vector<double> mass(n);
    double area = 0.0;
    for (std::size_t r = 0; r < n; ++r) {
        mass[r] = mass_all[global[r]];
        area += mass[r];
    }
    const double shift = 1.0 / area;
    std::vector<detail::SparseTriplet> kt, at;
    kt.reserve(stiff_all.size());
    at.reserve(stiff_all.size() + n);
    for (const detail::SparseTriplet& e : stiff_all) {
        const detail::SparseTriplet t{static_cast<std::size_t>(local[e.mRow]),
                                      static_cast<std::size_t>(local[e.mCol]), e.mValue};
        kt.push_back(t);
        at.push_back(t);
    }
    for (std::size_t r = 0; r < n; ++r)
        at.push_back({r, r, shift * mass[r]});
    const detail::SparseCsr stiffness = detail::csr_from_triplets(n, n, std::move(kt));
    const detail::SparseCsr shifted = detail::csr_from_triplets(n, n, std::move(at));

    const std::size_t p = std::min(n, std::max(2 * k, k + 8));
    std::vector<LapColumn> x(p, LapColumn(n));
    for (std::size_t j = 0; j < p; ++j)
        for (std::size_t r = 0; r < n; ++r)
            x[j][r] = lap_hash(global[r], j);
    std::uint64_t seed = p;
    lap_orthonormalise(x, mass, seed);
    std::vector<double> theta(p, std::numeric_limits<double>::quiet_NaN());

    EigenmodeResult result;
    result.mNumDegenerate = fem.mNumDegenerate;
    const std::int64_t cap = rOptions.mMaxIterations > 0 ? rOptions.mMaxIterations : 300;
    std::vector<LapColumn> y(p, LapColumn(n)), ky(p, LapColumn(n));
    std::vector<double> rhs(n);
    // Inexact inverse iteration: each solve only as accurate as the block
    // already is (Rayleigh-Ritz then judges the pairs exactly), tightening
    // as the residual falls.
    double solve_tol = 1e-4;
    for (std::int64_t it = 1; it <= cap; ++it) {
        // y_j = (K + s M)^-1 M x_j, warm-started from x_j / (theta_j + s).
        for (std::size_t j = 0; j < p; ++j) {
            for (std::size_t r = 0; r < n; ++r) {
                rhs[r] = mass[r] * x[j][r];
                y[j][r] = std::isfinite(theta[j]) ? x[j][r] / (theta[j] + shift) : 0.0;
            }
            detail::csr_pcg(shifted, rhs.data(), y[j].data(), solve_tol, 20 * n + 100);
        }
        lap_orthonormalise(y, mass, seed);
        // Rayleigh-Ritz: the projected stiffness in the M-orthonormal basis.
        for (std::size_t j = 0; j < p; ++j)
            detail::csr_multiply(stiffness, y[j].data(), ky[j].data());
        std::vector<double> kp(p * p);
        for (std::size_t i = 0; i < p; ++i)
            for (std::size_t j = i; j < p; ++j)
                kp[i * p + j] = kp[j * p + i] = detail::sparse_dot(y[i].data(), ky[j].data(), n);
        std::vector<double> w, v;
        lap_jacobi(kp, p, w, v);
        std::vector<std::size_t> order(p);
        std::iota(order.begin(), order.end(), std::size_t{0});
        std::stable_sort(order.begin(), order.end(),
                         [&](std::size_t A, std::size_t B) { return w[A] < w[B]; });
        for (std::size_t j = 0; j < p; ++j) {
            const std::size_t c = order[j];
            theta[j] = w[c];
            std::fill(x[j].begin(), x[j].end(), 0.0);
            for (std::size_t i = 0; i < p; ++i) {
                const double coef = v[i * p + c];
                for (std::size_t r = 0; r < n; ++r)
                    x[j][r] += coef * y[i][r];
            }
        }
        // The wanted pairs' residuals |K x - theta M x|_{M^-1} / (theta + s).
        double worst = 0.0;
        std::vector<double> kx(n);
        for (std::size_t j = 0; j < k; ++j) {
            detail::csr_multiply(stiffness, x[j].data(), kx.data());
            double s = 0.0;
            for (std::size_t r = 0; r < n; ++r) {
                const double d = kx[r] - theta[j] * mass[r] * x[j][r];
                s += d * d / mass[r];
            }
            worst = std::max(worst, std::sqrt(s) / (std::abs(theta[j]) + shift));
        }
        result.mIterations = it;
        result.mMaxResidual = worst;
        if (worst <= rOptions.mTolerance) {
            result.mConverged = true;
            break;
        }
        solve_tol = std::max(1e-14, std::min(solve_tol, 1e-2 * worst));
    }

    NDArray modes(DType::Float64, {n_all, k});
    double* pm = modes.As<double>();
    std::fill(pm, pm + n_all * k, std::numeric_limits<double>::quiet_NaN());
    NDArray values(DType::Float64, {k});
    for (std::size_t j = 0; j < k; ++j) {
        std::size_t big = 0;
        for (std::size_t r = 1; r < n; ++r)
            if (std::abs(x[j][r]) > std::abs(x[j][big]))
                big = r;
        const double sign = x[j][big] < 0.0 ? -1.0 : 1.0;
        for (std::size_t r = 0; r < n; ++r)
            pm[global[r] * k + j] = sign * x[j][r];
        // A zero eigenvalue comes out as rounding either side of it.
        const double lambda = std::abs(theta[j]) <= 1e-12 * shift ? 0.0 : theta[j];
        values.As<double>()[j] = lambda;
        result.mEigenvalues.push_back(lambda);
    }
    result.mMesh = detail::clone_mesh(rMesh);
    result.mMesh.AddPointData(kLaplacianModesName, std::move(modes));
    result.mMesh.AddFieldData(kLaplacianEigenvaluesName, std::move(values));
    return result;
}

}  // namespace meshioplusplus

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

// Dense helpers for fit_primitives. See lsq_fit.hpp.

// System includes
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

// Project includes (private, not installed)
#include "lsq_fit.hpp"

namespace meshioplusplus {
namespace detail {

void lsq_sym3_eigen(const double* pA, double* pW, double* pV) {
    double a[3][3];
    double v[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            a[i][j] = pA[i * 3 + j];
    static const int kPq[3][2] = {{0, 1}, {0, 2}, {1, 2}};
    for (int sweep = 0; sweep < 50; ++sweep) {
        const double off = a[0][1] * a[0][1] + a[0][2] * a[0][2] + a[1][2] * a[1][2];
        const double diag = a[0][0] * a[0][0] + a[1][1] * a[1][1] + a[2][2] * a[2][2];
        if (off == 0.0 || off <= 1e-32 * diag)
            break;
        for (const auto& pq : kPq) {
            const int p = pq[0], q = pq[1];
            if (a[p][q] == 0.0)
                continue;
            const double theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
            const double t =
                (theta >= 0.0 ? 1.0 : -1.0) / (std::fabs(theta) + std::sqrt(theta * theta + 1.0));
            const double c = 1.0 / std::sqrt(t * t + 1.0);
            const double s = t * c;
            for (int k = 0; k < 3; ++k) {
                const double akp = a[k][p], akq = a[k][q];
                a[k][p] = c * akp - s * akq;
                a[k][q] = s * akp + c * akq;
            }
            for (int k = 0; k < 3; ++k) {
                const double apk = a[p][k], aqk = a[q][k];
                a[p][k] = c * apk - s * aqk;
                a[q][k] = s * apk + c * aqk;
            }
            a[p][q] = 0.0;
            a[q][p] = 0.0;
            for (int k = 0; k < 3; ++k) {
                const double vkp = v[k][p], vkq = v[k][q];
                v[k][p] = c * vkp - s * vkq;
                v[k][q] = s * vkp + c * vkq;
            }
        }
    }
    int order[3] = {0, 1, 2};
    for (int i = 1; i < 3; ++i)
        for (int j = i; j > 0 && a[order[j]][order[j]] < a[order[j - 1]][order[j - 1]]; --j)
            std::swap(order[j], order[j - 1]);
    for (int r = 0; r < 3; ++r) {
        const int col = order[r];
        pW[r] = a[col][col];
        int big = 0;
        for (int k = 1; k < 3; ++k)
            if (std::fabs(v[k][col]) > std::fabs(v[big][col]))
                big = k;
        const double sign = v[big][col] < 0.0 ? -1.0 : 1.0;
        for (int k = 0; k < 3; ++k)
            pV[r * 3 + k] = sign * v[k][col];
    }
}

bool lsq_solve(std::vector<double>& rA, std::vector<double>& rB, std::size_t K) {
    for (std::size_t col = 0; col < K; ++col) {
        std::size_t piv = col;
        for (std::size_t r = col + 1; r < K; ++r)
            if (std::fabs(rA[r * K + col]) > std::fabs(rA[piv * K + col]))
                piv = r;
        const double pv = rA[piv * K + col];
        if (pv == 0.0 || !std::isfinite(pv))
            return false;
        if (piv != col) {
            for (std::size_t k = 0; k < K; ++k)
                std::swap(rA[piv * K + k], rA[col * K + k]);
            std::swap(rB[piv], rB[col]);
        }
        for (std::size_t r = col + 1; r < K; ++r) {
            const double f = rA[r * K + col] / pv;
            if (f == 0.0)
                continue;
            for (std::size_t k = col; k < K; ++k)
                rA[r * K + k] -= f * rA[col * K + k];
            rB[r] -= f * rB[col];
        }
    }
    for (std::size_t i = K; i-- > 0;) {
        double s = rB[i];
        for (std::size_t k = i + 1; k < K; ++k)
            s -= rA[i * K + k] * rB[k];
        rB[i] = s / rA[i * K + i];
        if (!std::isfinite(rB[i]))
            return false;
    }
    return true;
}

void lsq_frame(const double* pE3, double* pE1, double* pE2) {
    int axis = 0;
    for (int k = 1; k < 3; ++k)
        if (std::fabs(pE3[k]) < std::fabs(pE3[axis]))
            axis = k;
    double h[3] = {0.0, 0.0, 0.0};
    h[axis] = 1.0;
    const double d = h[0] * pE3[0] + h[1] * pE3[1] + h[2] * pE3[2];
    double e1[3] = {h[0] - d * pE3[0], h[1] - d * pE3[1], h[2] - d * pE3[2]};
    const double len = std::sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]);
    for (int k = 0; k < 3; ++k)
        pE1[k] = e1[k] / len;
    pE2[0] = pE3[1] * pE1[2] - pE3[2] * pE1[1];
    pE2[1] = pE3[2] * pE1[0] - pE3[0] * pE1[2];
    pE2[2] = pE3[0] * pE1[1] - pE3[1] * pE1[0];
}

namespace {

double lsq_cost(const std::vector<double>& rR) {
    double s = 0.0;
    for (double x : rR)
        s += x * x;
    return s;
}

}  // namespace

void lsq_levenberg_marquardt(std::vector<double>& rTheta, const std::vector<double>& rTypical,
                             std::size_t N, const LsqResidual& rResidual, int MaxIterations) {
    const std::size_t k = rTheta.size();
    if (N < k || k == 0)
        return;
    std::vector<double> r(N), rn(N), rj(N), jac(N * k);
    rResidual(rTheta, r.data());
    double cost = lsq_cost(r);
    if (!std::isfinite(cost))
        return;
    double lambda = 1e-3;
    std::vector<double> h(k * k), g(k), m, delta, trial(k), step(k);
    for (int it = 0; it < MaxIterations; ++it) {
        for (std::size_t j = 0; j < k; ++j) {
            std::vector<double> t = rTheta;
            const double mag =
                std::fabs(rTheta[j]) > rTypical[j] ? std::fabs(rTheta[j]) : rTypical[j];
            step[j] = 1.4901161193847656e-08 * mag;
            t[j] += step[j];
            rResidual(t, rj.data());
            for (std::size_t i = 0; i < N; ++i)
                jac[i * k + j] = (rj[i] - r[i]) / step[j];
        }
        for (std::size_t a = 0; a < k; ++a) {
            double ga = 0.0;
            for (std::size_t i = 0; i < N; ++i)
                ga += jac[i * k + a] * r[i];
            g[a] = ga;
            for (std::size_t b = 0; b < k; ++b) {
                double hab = 0.0;
                for (std::size_t i = 0; i < N; ++i)
                    hab += jac[i * k + a] * jac[i * k + b];
                h[a * k + b] = hab;
            }
        }
        bool accepted = false;
        double cost_new = cost;
        for (int attempt = 0; attempt < 12; ++attempt) {
            m = h;
            delta.assign(k, 0.0);
            for (std::size_t a = 0; a < k; ++a) {
                m[a * k + a] += lambda * (h[a * k + a] > 1e-300 ? h[a * k + a] : 1e-300);
                delta[a] = -g[a];
            }
            if (!lsq_solve(m, delta, k)) {
                lambda *= 10.0;
                continue;
            }
            for (std::size_t a = 0; a < k; ++a)
                trial[a] = rTheta[a] + delta[a];
            rResidual(trial, rn.data());
            cost_new = lsq_cost(rn);
            if (std::isfinite(cost_new) && cost_new < cost) {
                accepted = true;
                break;
            }
            lambda *= 10.0;
        }
        if (!accepted)
            break;
        const double before = cost;
        rTheta = trial;
        r.swap(rn);
        cost = cost_new;
        lambda = lambda * 0.1 > 1e-15 ? lambda * 0.1 : 1e-15;
        if (before - cost <= 1e-14 * before)
            break;
    }
}

}  // namespace detail
}  // namespace meshioplusplus

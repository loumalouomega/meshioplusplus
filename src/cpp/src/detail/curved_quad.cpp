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
#include <array>
#include <cmath>
#include <cstddef>

// Project includes
#include "curved_quad.hpp"

namespace meshioplusplus::detail {

namespace {

// The quadratic Lagrange basis at 0, 1/2 and 1, and its derivative.
std::array<double, 3> cq_basis(double T) {
    return {2.0 * T * T - 3.0 * T + 1.0, 4.0 * T - 4.0 * T * T, 2.0 * T * T - T};
}

std::array<double, 3> cq_basis_d(double T) {
    return {4.0 * T - 3.0, 4.0 - 8.0 * T, 4.0 * T - 1.0};
}

// Node k of quad9 at lattice position (i, j), i along xi, j along eta.
constexpr std::array<std::array<int, 2>, 9> kCqLattice{
    {{0, 0}, {2, 0}, {2, 2}, {0, 2}, {1, 0}, {2, 1}, {1, 2}, {0, 1}, {1, 1}}};

// The inverse of the cubic Bernstein basis at 0, 1/3, 2/3 and 1: the
// coefficients from the values, row by row.
const std::array<std::array<double, 4>, 4>& cq_bernstein_inverse() {
    static const std::array<std::array<double, 4>, 4> inv = [] {
        std::array<std::array<double, 8>, 4> m{};
        const double binom[4] = {1.0, 3.0, 3.0, 1.0};
        for (int i = 0; i < 4; ++i) {
            const double t = double(i) / 3.0;
            for (int k = 0; k < 4; ++k)
                m[i][k] = binom[k] * std::pow(t, k) * std::pow(1.0 - t, 3 - k);
            m[i][4 + i] = 1.0;
        }
        // Gauss-Jordan with partial pivoting on [M | I].
        for (int c = 0; c < 4; ++c) {
            int p = c;
            for (int r = c + 1; r < 4; ++r)
                if (std::fabs(m[r][c]) > std::fabs(m[p][c]))
                    p = r;
            std::swap(m[c], m[p]);
            const double d = m[c][c];
            for (double& r_v : m[c])
                r_v /= d;
            for (int r = 0; r < 4; ++r)
                if (r != c) {
                    const double f = m[r][c];
                    for (int k = 0; k < 8; ++k)
                        m[r][k] -= f * m[c][k];
                }
        }
        std::array<std::array<double, 4>, 4> out{};
        for (int i = 0; i < 4; ++i)
            for (int k = 0; k < 4; ++k)
                out[i][k] = m[i][4 + k];
        return out;
    }();
    return inv;
}

// Judge the part over [U0, U1] x [V0, V1]: 1 valid, 0 invalid, -1
// undecided at the depth.
int cq_judge(const CurvedVec3* pNodes, const CurvedVec3& rN, double U0, double U1, double V0,
             double V1, int Depth, int MaxDepth) {
    std::array<std::array<double, 4>, 4> vals{};
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) {
            vals[i][j] = curved_quad_jacobian(pNodes, rN, U0 + (U1 - U0) * double(i) / 3.0,
                                              V0 + (V1 - V0) * double(j) / 3.0);
            if (!(vals[i][j] > 0.0))
                return 0;
        }
    // C = Minv V Minv^T.
    const auto& r_inv = cq_bernstein_inverse();
    std::array<std::array<double, 4>, 4> tmp{}, coef{};
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            for (int k = 0; k < 4; ++k)
                tmp[i][j] += r_inv[i][k] * vals[k][j];
    double least = tmp[0][0];
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) {
            for (int k = 0; k < 4; ++k)
                coef[i][j] += tmp[i][k] * r_inv[j][k];
            least = std::min(least, coef[i][j]);
        }
    if (least > 0.0)
        return 1;
    if (Depth >= MaxDepth)
        return -1;
    const double um = 0.5 * (U0 + U1), vm = 0.5 * (V0 + V1);
    int verdict = 1;
    for (const auto& r_b :
         {std::array<double, 4>{U0, um, V0, vm}, std::array<double, 4>{um, U1, V0, vm},
          std::array<double, 4>{U0, um, vm, V1}, std::array<double, 4>{um, U1, vm, V1}}) {
        const int v = cq_judge(pNodes, rN, r_b[0], r_b[1], r_b[2], r_b[3], Depth + 1, MaxDepth);
        if (v == 0)
            return 0;
        if (v < 0)
            verdict = -1;
    }
    return verdict;
}

}  // namespace

CurvedVec3 curved_quad_point(const CurvedVec3* pNodes, double Xi, double Eta) {
    const auto a = cq_basis(Xi), b = cq_basis(Eta);
    CurvedVec3 x{0.0, 0.0, 0.0};
    for (std::size_t n = 0; n < 9; ++n) {
        const double w = a[static_cast<std::size_t>(kCqLattice[n][0])] *
                         b[static_cast<std::size_t>(kCqLattice[n][1])];
        for (std::size_t k = 0; k < 3; ++k)
            x[k] += w * pNodes[n][k];
    }
    return x;
}

double curved_quad_jacobian(const CurvedVec3* pNodes, const CurvedVec3& rNormal, double Xi,
                            double Eta) {
    const auto a = cq_basis(Xi), b = cq_basis(Eta);
    const auto da = cq_basis_d(Xi), db = cq_basis_d(Eta);
    CurvedVec3 xu{0.0, 0.0, 0.0}, xv{0.0, 0.0, 0.0};
    for (std::size_t n = 0; n < 9; ++n) {
        const auto i = static_cast<std::size_t>(kCqLattice[n][0]);
        const auto j = static_cast<std::size_t>(kCqLattice[n][1]);
        for (std::size_t k = 0; k < 3; ++k) {
            xu[k] += da[i] * b[j] * pNodes[n][k];
            xv[k] += a[i] * db[j] * pNodes[n][k];
        }
    }
    const CurvedVec3 c{xu[1] * xv[2] - xu[2] * xv[1], xu[2] * xv[0] - xu[0] * xv[2],
                       xu[0] * xv[1] - xu[1] * xv[0]};
    return c[0] * rNormal[0] + c[1] * rNormal[1] + c[2] * rNormal[2];
}

CurvedTriangleCheck curved_quad_check(const CurvedVec3* pNodes, const CurvedVec3& rNormal,
                                      int MaxDepth) {
    CurvedTriangleCheck out;
    out.mMinNodal = curved_quad_jacobian(pNodes, rNormal, 0.0, 0.0);
    for (std::size_t n = 1; n < 9; ++n)
        out.mMinNodal = std::min(
            out.mMinNodal,
            curved_quad_jacobian(pNodes, rNormal, 0.5 * kCqLattice[n][0], 0.5 * kCqLattice[n][1]));
    const int v = cq_judge(pNodes, rNormal, 0.0, 1.0, 0.0, 1.0, 0, MaxDepth);
    out.mValid = v == 1;
    out.mDecided = v >= 0;
    return out;
}

}  // namespace meshioplusplus::detail

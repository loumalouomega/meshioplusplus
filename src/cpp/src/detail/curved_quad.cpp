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
#include <limits>

// Project includes
#include "curved_quad.hpp"

namespace meshioplusplus::detail {

namespace {

// The Lagrange basis of degree Order (2 or 3) at the equispaced points
// 0, 1/Order, ..., 1, and its derivative.
std::array<double, 4> cq_basis(int Order, double T) {
    if (Order == 2)
        return {2.0 * T * T - 3.0 * T + 1.0, 4.0 * T - 4.0 * T * T, 2.0 * T * T - T, 0.0};
    const double a = 3.0 * T, b = 3.0 * T - 1.0, c = 3.0 * T - 2.0, d = 3.0 * T - 3.0;
    return {-b * c * d / 6.0, a * c * d / 2.0, -a * b * d / 2.0, a * b * c / 6.0};
}

std::array<double, 4> cq_basis_d(int Order, double T) {
    if (Order == 2)
        return {4.0 * T - 3.0, 4.0 - 8.0 * T, 4.0 * T - 1.0, 0.0};
    const double a = 3.0 * T, b = 3.0 * T - 1.0, c = 3.0 * T - 2.0, d = 3.0 * T - 3.0;
    // d/dT of each product of three factors, each factor's derivative 3.
    return {-0.5 * (c * d + b * d + b * c), 1.5 * (c * d + a * d + a * c),
            -1.5 * (b * d + a * d + a * b), 0.5 * (b * c + a * c + a * b)};
}

// Node k of quad9 and quad16 at lattice position (i, j), i along xi, j along
// eta, in steps of 1/Order: Gmsh's layout (corners, then each edge's nodes
// from its first corner, then the inside, itself a quad's corners).
constexpr std::array<std::array<int, 2>, 9> kCqLattice2{
    {{0, 0}, {2, 0}, {2, 2}, {0, 2}, {1, 0}, {2, 1}, {1, 2}, {0, 1}, {1, 1}}};
constexpr std::array<std::array<int, 2>, 16> kCqLattice3{{{0, 0},
                                                          {3, 0},
                                                          {3, 3},
                                                          {0, 3},
                                                          {1, 0},
                                                          {2, 0},
                                                          {3, 1},
                                                          {3, 2},
                                                          {2, 3},
                                                          {1, 3},
                                                          {0, 2},
                                                          {0, 1},
                                                          {1, 1},
                                                          {2, 1},
                                                          {2, 2},
                                                          {1, 2}}};

std::array<int, 2> cq_lattice(int Order, std::size_t N) {
    return Order == 2 ? kCqLattice2[N] : kCqLattice3[N];
}

// The Jacobian's degree in each direction: 2 Order - 1.
constexpr int kCqMaxDegree = 5;
using CqMatrix = std::array<std::array<double, kCqMaxDegree + 1>, kCqMaxDegree + 1>;

// The inverse of the Bernstein basis of degree Degree at the equispaced
// points 0, 1/Degree, ..., 1: the coefficients from the values, row by row.
CqMatrix cq_make_bernstein_inverse(int Degree) {
    const int n = Degree + 1;
    std::array<std::array<double, 2 * (kCqMaxDegree + 1)>, kCqMaxDegree + 1> m{};
    for (int i = 0; i < n; ++i) {
        const double t = double(i) / double(Degree);
        double binom = 1.0;
        for (int k = 0; k < n; ++k) {
            m[i][k] = binom * std::pow(t, k) * std::pow(1.0 - t, Degree - k);
            binom = binom * double(Degree - k) / double(k + 1);
        }
        m[i][n + i] = 1.0;
    }
    // Gauss-Jordan with partial pivoting on [M | I].
    for (int c = 0; c < n; ++c) {
        int p = c;
        for (int r = c + 1; r < n; ++r)
            if (std::fabs(m[r][c]) > std::fabs(m[p][c]))
                p = r;
        std::swap(m[c], m[p]);
        const double d = m[c][c];
        for (int k = 0; k < 2 * n; ++k)
            m[c][k] /= d;
        for (int r = 0; r < n; ++r)
            if (r != c) {
                const double f = m[r][c];
                for (int k = 0; k < 2 * n; ++k)
                    m[r][k] -= f * m[c][k];
            }
    }
    CqMatrix out{};
    for (int i = 0; i < n; ++i)
        for (int k = 0; k < n; ++k)
            out[i][k] = m[i][n + k];
    return out;
}

const CqMatrix& cq_bernstein_inverse(int Degree) {
    static const CqMatrix inv3 = cq_make_bernstein_inverse(3);
    static const CqMatrix inv5 = cq_make_bernstein_inverse(5);
    return Degree == 3 ? inv3 : inv5;
}

// Judge the part over [U0, U1] x [V0, V1]: 1 valid, 0 invalid, -1
// undecided at the depth.
int cq_judge(int Order, const CurvedVec3* pNodes, const CurvedVec3& rN, double U0, double U1,
             double V0, double V1, int Depth, int MaxDepth) {
    const int degree = 2 * Order - 1;
    const int n = degree + 1;
    CqMatrix vals{};
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            vals[i][j] =
                curved_quad_jacobian(Order, pNodes, rN, U0 + (U1 - U0) * double(i) / degree,
                                     V0 + (V1 - V0) * double(j) / degree);
            if (!(vals[i][j] > 0.0))
                return 0;
        }
    // C = Minv V Minv^T.
    const CqMatrix& r_inv = cq_bernstein_inverse(degree);
    CqMatrix tmp{}, coef{};
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; ++k)
                tmp[i][j] += r_inv[i][k] * vals[k][j];
    double least = tmp[0][0];
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k)
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
        const int v =
            cq_judge(Order, pNodes, rN, r_b[0], r_b[1], r_b[2], r_b[3], Depth + 1, MaxDepth);
        if (v == 0)
            return 0;
        if (v < 0)
            verdict = -1;
    }
    return verdict;
}

}  // namespace

std::size_t curved_quad_num_nodes(int Order) {
    return static_cast<std::size_t>((Order + 1) * (Order + 1));
}

std::array<double, 2> curved_quad_node_position(int Order, std::size_t Node) {
    const auto l = cq_lattice(Order, Node);
    return {double(l[0]) / Order, double(l[1]) / Order};
}

CurvedVec3 curved_quad_point(int Order, const CurvedVec3* pNodes, double Xi, double Eta) {
    const auto a = cq_basis(Order, Xi), b = cq_basis(Order, Eta);
    CurvedVec3 x{0.0, 0.0, 0.0};
    for (std::size_t n = 0; n < curved_quad_num_nodes(Order); ++n) {
        const auto l = cq_lattice(Order, n);
        const double w = a[static_cast<std::size_t>(l[0])] * b[static_cast<std::size_t>(l[1])];
        for (std::size_t k = 0; k < 3; ++k)
            x[k] += w * pNodes[n][k];
    }
    return x;
}

double curved_quad_jacobian(int Order, const CurvedVec3* pNodes, const CurvedVec3& rNormal,
                            double Xi, double Eta) {
    const auto a = cq_basis(Order, Xi), b = cq_basis(Order, Eta);
    const auto da = cq_basis_d(Order, Xi), db = cq_basis_d(Order, Eta);
    CurvedVec3 xu{0.0, 0.0, 0.0}, xv{0.0, 0.0, 0.0};
    for (std::size_t n = 0; n < curved_quad_num_nodes(Order); ++n) {
        const auto l = cq_lattice(Order, n);
        const auto i = static_cast<std::size_t>(l[0]);
        const auto j = static_cast<std::size_t>(l[1]);
        for (std::size_t k = 0; k < 3; ++k) {
            xu[k] += da[i] * b[j] * pNodes[n][k];
            xv[k] += a[i] * db[j] * pNodes[n][k];
        }
    }
    const CurvedVec3 c{xu[1] * xv[2] - xu[2] * xv[1], xu[2] * xv[0] - xu[0] * xv[2],
                       xu[0] * xv[1] - xu[1] * xv[0]};
    return c[0] * rNormal[0] + c[1] * rNormal[1] + c[2] * rNormal[2];
}

CurvedTriangleCheck curved_quad_check(int Order, const CurvedVec3* pNodes,
                                      const CurvedVec3& rNormal, int MaxDepth) {
    CurvedTriangleCheck out;
    out.mMinNodal = curved_quad_jacobian(Order, pNodes, rNormal, 0.0, 0.0);
    for (std::size_t n = 1; n < curved_quad_num_nodes(Order); ++n) {
        const auto p = curved_quad_node_position(Order, n);
        out.mMinNodal =
            std::min(out.mMinNodal, curved_quad_jacobian(Order, pNodes, rNormal, p[0], p[1]));
    }
    const int v = cq_judge(Order, pNodes, rNormal, 0.0, 1.0, 0.0, 1.0, 0, MaxDepth);
    out.mValid = v == 1;
    out.mDecided = v >= 0;
    return out;
}

double curved_quad_scaled_bound(int Order, const CurvedVec3* pNodes, const CurvedVec3& rNormal) {
    // As curved_triangle_scaled_bound: the least Bernstein coefficient of J
    // over the largest in magnitude, on the tensor-product basis.
    const int degree = 2 * Order - 1;
    const int n = degree + 1;
    CqMatrix vals{};
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            vals[i][j] = curved_quad_jacobian(Order, pNodes, rNormal, double(i) / degree,
                                              double(j) / degree);
    const CqMatrix& r_inv = cq_bernstein_inverse(degree);
    CqMatrix tmp{};
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; ++k)
                tmp[i][j] += r_inv[i][k] * vals[k][j];
    double least = std::numeric_limits<double>::infinity(), largest = 0.0;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            double c = 0.0;
            for (int k = 0; k < n; ++k)
                c += tmp[i][k] * r_inv[j][k];
            least = std::min(least, c);
            largest = std::max(largest, std::fabs(c));
        }
    return largest > 0.0 ? least / largest : -1.0;
}

CurvedVec3 curved_quad_point(const CurvedVec3* pNodes, double Xi, double Eta) {
    return curved_quad_point(2, pNodes, Xi, Eta);
}

double curved_quad_jacobian(const CurvedVec3* pNodes, const CurvedVec3& rNormal, double Xi,
                            double Eta) {
    return curved_quad_jacobian(2, pNodes, rNormal, Xi, Eta);
}

CurvedTriangleCheck curved_quad_check(const CurvedVec3* pNodes, const CurvedVec3& rNormal,
                                      int MaxDepth) {
    return curved_quad_check(2, pNodes, rNormal, MaxDepth);
}

}  // namespace meshioplusplus::detail

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
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "curved_triangle.hpp"

namespace meshioplusplus::detail {

namespace {

constexpr int kCtMaxNodes = 10;

/// The shape functions of the Lagrange triangle of `Order` at `(Xi, Eta)`
/// and their derivatives, through the barycentric coordinates L1 = 1 - xi -
/// eta, L2 = xi, L3 = eta (corners 0, 1, 2).
void ct_shape(int Order, double Xi, double Eta, double* pN, double* pDXi, double* pDEta) {
    const double l[3] = {1.0 - Xi - Eta, Xi, Eta};
    const double dl[3][2] = {{-1.0, -1.0}, {1.0, 0.0}, {0.0, 1.0}};
    // N, dN/dL_a and dN/dL_b for a function of two barycentrics, and the
    // chain rule into (xi, eta).
    int n = 0;
    const auto put = [&](double N, std::array<double, 3> DL) {
        if (pN)
            pN[n] = N;
        if (pDXi)
            pDXi[n] = DL[0] * dl[0][0] + DL[1] * dl[1][0] + DL[2] * dl[2][0];
        if (pDEta)
            pDEta[n] = DL[0] * dl[0][1] + DL[1] * dl[1][1] + DL[2] * dl[2][1];
        ++n;
    };
    const auto unit = [](int A, double V) {
        std::array<double, 3> d{0.0, 0.0, 0.0};
        d[static_cast<std::size_t>(A)] = V;
        return d;
    };
    const auto pair = [](int A, double Va, int B, double Vb) {
        std::array<double, 3> d{0.0, 0.0, 0.0};
        d[static_cast<std::size_t>(A)] += Va;
        d[static_cast<std::size_t>(B)] += Vb;
        return d;
    };
    if (Order == 1) {
        for (int a = 0; a < 3; ++a)
            put(l[a], unit(a, 1.0));
    } else if (Order == 2) {
        for (int a = 0; a < 3; ++a)
            put(l[a] * (2.0 * l[a] - 1.0), unit(a, 4.0 * l[a] - 1.0));
        // Mid-edge nodes of edges (0, 1), (1, 2), (2, 0).
        for (int a = 0; a < 3; ++a) {
            const int b = (a + 1) % 3;
            put(4.0 * l[a] * l[b], pair(a, 4.0 * l[b], b, 4.0 * l[a]));
        }
    } else if (Order == 3) {
        for (int a = 0; a < 3; ++a)
            put(0.5 * l[a] * (3.0 * l[a] - 1.0) * (3.0 * l[a] - 2.0),
                unit(a, 0.5 * (27.0 * l[a] * l[a] - 18.0 * l[a] + 2.0)));
        // Each edge (a, b): the node at L_a = 2/3 then the one at L_b = 2/3;
        // the node nearer `c` is (9/2) L_c L_o (3 L_c - 1).
        for (int a = 0; a < 3; ++a) {
            const int b = (a + 1) % 3;
            for (const auto& [c, o] : {std::pair<int, int>{a, b}, std::pair<int, int>{b, a}})
                put(4.5 * l[c] * l[o] * (3.0 * l[c] - 1.0),
                    pair(c, 4.5 * l[o] * (6.0 * l[c] - 1.0), o, 4.5 * l[c] * (3.0 * l[c] - 1.0)));
        }
        std::array<double, 3> d{27.0 * l[1] * l[2], 27.0 * l[0] * l[2], 27.0 * l[0] * l[1]};
        put(27.0 * l[0] * l[1] * l[2], d);
    } else {
        throw std::invalid_argument("meshio++: curved_triangle: order must be 1, 2 or 3, not " +
                                    std::to_string(Order));
    }
}

double ct_multinomial(int D, int I, int J, int K) {
    double f = 1.0;
    for (int k = 2; k <= D; ++k)
        f *= k;
    for (const int m : {I, J, K})
        for (int k = 2; k <= m; ++k)
            f /= k;
    return f;
}

/// The lattice of degree `D`: `(i, j)`, `i` outer, `i + j <= D`.
std::vector<std::array<int, 2>> ct_lattice(int D) {
    std::vector<std::array<int, 2>> out;
    for (int i = 0; i <= D; ++i)
        for (int j = 0; i + j <= D; ++j)
            out.push_back({i, j});
    return out;
}

/// The inverse of the Bernstein collocation matrix of degree `D` (values at
/// the lattice to coefficients), row-major, computed once.
const std::vector<double>& ct_bernstein_inverse(int D) {
    const auto build = [](int Deg) {
        const auto lat = ct_lattice(Deg);
        const std::size_t n = lat.size();
        std::vector<double> m(n * n), inv(n * n, 0.0);
        for (std::size_t a = 0; a < n; ++a) {
            const double x = double(lat[a][0]) / Deg, y = double(lat[a][1]) / Deg;
            const double z = 1.0 - x - y;
            for (std::size_t b = 0; b < n; ++b) {
                const int i = lat[b][0], j = lat[b][1], k = Deg - i - j;
                m[a * n + b] =
                    ct_multinomial(Deg, i, j, k) * std::pow(x, i) * std::pow(y, j) * std::pow(z, k);
            }
            inv[a * n + a] = 1.0;
        }
        // Gauss-Jordan with partial pivoting.
        for (std::size_t c = 0; c < n; ++c) {
            std::size_t p = c;
            for (std::size_t r = c + 1; r < n; ++r)
                if (std::fabs(m[r * n + c]) > std::fabs(m[p * n + c]))
                    p = r;
            for (std::size_t k = 0; k < n; ++k) {
                std::swap(m[c * n + k], m[p * n + k]);
                std::swap(inv[c * n + k], inv[p * n + k]);
            }
            const double piv = m[c * n + c];
            for (std::size_t k = 0; k < n; ++k) {
                m[c * n + k] /= piv;
                inv[c * n + k] /= piv;
            }
            for (std::size_t r = 0; r < n; ++r) {
                if (r == c)
                    continue;
                const double f = m[r * n + c];
                if (f == 0.0)
                    continue;
                for (std::size_t k = 0; k < n; ++k) {
                    m[r * n + k] -= f * m[c * n + k];
                    inv[r * n + k] -= f * inv[c * n + k];
                }
            }
        }
        return inv;
    };
    static const std::vector<double> kInv2 = build(2);
    static const std::vector<double> kInv4 = build(4);
    if (D == 2)
        return kInv2;
    if (D == 4)
        return kInv4;
    throw std::invalid_argument("meshio++: curved_triangle: Bernstein degree must be 2 or 4");
}

using CtPoint = std::array<double, 2>;

/// Judge the part of the element over the reference sub-triangle `A`, `B`,
/// `C`: 1 valid, 0 invalid, -1 undecided at the depth.
int ct_judge(int Order, const CurvedVec3* pNodes, const CurvedVec3& rN, const CtPoint& rA,
             const CtPoint& rB, const CtPoint& rC, int Depth, int MaxDepth) {
    const int d = 2 * (Order - 1);
    const auto lat = ct_lattice(d);
    std::array<double, 15> vals{};
    for (std::size_t a = 0; a < lat.size(); ++a) {
        const double s = double(lat[a][0]) / d, t = double(lat[a][1]) / d;
        const double xi = rA[0] + s * (rB[0] - rA[0]) + t * (rC[0] - rA[0]);
        const double eta = rA[1] + s * (rB[1] - rA[1]) + t * (rC[1] - rA[1]);
        vals[a] = curved_triangle_jacobian(Order, pNodes, rN, xi, eta);
        if (!(vals[a] > 0.0))
            return 0;
    }
    std::array<double, 15> coef{};
    curved_triangle_bernstein(d, vals.data(), coef.data());
    if (*std::min_element(coef.begin(), coef.begin() + static_cast<std::ptrdiff_t>(lat.size())) >
        0.0)
        return 1;
    if (Depth >= MaxDepth)
        return -1;
    const CtPoint ab{0.5 * (rA[0] + rB[0]), 0.5 * (rA[1] + rB[1])};
    const CtPoint bc{0.5 * (rB[0] + rC[0]), 0.5 * (rB[1] + rC[1])};
    const CtPoint ca{0.5 * (rC[0] + rA[0]), 0.5 * (rC[1] + rA[1])};
    int verdict = 1;
    for (const auto& r_t :
         {std::array<CtPoint, 3>{rA, ab, ca}, std::array<CtPoint, 3>{ab, rB, bc},
          std::array<CtPoint, 3>{ca, bc, rC}, std::array<CtPoint, 3>{ab, bc, ca}}) {
        const int v = ct_judge(Order, pNodes, rN, r_t[0], r_t[1], r_t[2], Depth + 1, MaxDepth);
        if (v == 0)
            return 0;
        if (v < 0)
            verdict = -1;
    }
    return verdict;
}

}  // namespace

int curved_triangle_num_nodes(int Order) {
    if (Order < 1 || Order > 3)
        throw std::invalid_argument("meshio++: curved_triangle: order must be 1, 2 or 3, not " +
                                    std::to_string(Order));
    return (Order + 1) * (Order + 2) / 2;
}

CurvedVec3 curved_triangle_point(int Order, const CurvedVec3* pNodes, double Xi, double Eta) {
    double shape[kCtMaxNodes];
    ct_shape(Order, Xi, Eta, shape, nullptr, nullptr);
    const int n = curved_triangle_num_nodes(Order);
    CurvedVec3 x{0.0, 0.0, 0.0};
    for (int i = 0; i < n; ++i)
        for (std::size_t k = 0; k < 3; ++k)
            x[k] += shape[i] * pNodes[i][k];
    return x;
}

double curved_triangle_jacobian(int Order, const CurvedVec3* pNodes, const CurvedVec3& rNormal,
                                double Xi, double Eta) {
    double dxi[kCtMaxNodes], deta[kCtMaxNodes];
    ct_shape(Order, Xi, Eta, nullptr, dxi, deta);
    const int n = curved_triangle_num_nodes(Order);
    CurvedVec3 a{0.0, 0.0, 0.0}, b{0.0, 0.0, 0.0};
    for (int i = 0; i < n; ++i)
        for (std::size_t k = 0; k < 3; ++k) {
            a[k] += dxi[i] * pNodes[i][k];
            b[k] += deta[i] * pNodes[i][k];
        }
    const CurvedVec3 c{a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
                       a[0] * b[1] - a[1] * b[0]};
    return c[0] * rNormal[0] + c[1] * rNormal[1] + c[2] * rNormal[2];
}

void curved_triangle_bernstein(int Degree, const double* pValues, double* pCoefficients) {
    const std::vector<double>& inv = ct_bernstein_inverse(Degree);
    const std::size_t n = static_cast<std::size_t>((Degree + 1) * (Degree + 2) / 2);
    for (std::size_t r = 0; r < n; ++r) {
        double s = 0.0;
        for (std::size_t c = 0; c < n; ++c)
            s += inv[r * n + c] * pValues[c];
        pCoefficients[r] = s;
    }
}

CurvedTriangleCheck curved_triangle_check(int Order, const CurvedVec3* pNodes,
                                          const CurvedVec3& rNormal, int MaxDepth) {
    CurvedTriangleCheck out;
    // The nodes' reference positions: the lattice of the order, in the
    // element's node layout.
    static const std::array<CtPoint, 10> kP2{
        {{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}, {0.5, 0.0}, {0.5, 0.5}, {0.0, 0.5}}};
    static const std::array<CtPoint, 10> kP3{{{0.0, 0.0},
                                              {1.0, 0.0},
                                              {0.0, 1.0},
                                              {1.0 / 3.0, 0.0},
                                              {2.0 / 3.0, 0.0},
                                              {2.0 / 3.0, 1.0 / 3.0},
                                              {1.0 / 3.0, 2.0 / 3.0},
                                              {0.0, 2.0 / 3.0},
                                              {0.0, 1.0 / 3.0},
                                              {1.0 / 3.0, 1.0 / 3.0}}};
    const int n = curved_triangle_num_nodes(Order);
    const auto& r_ref = Order == 3 ? kP3 : kP2;
    out.mMinNodal = curved_triangle_jacobian(Order, pNodes, rNormal, 0.0, 0.0);
    for (int i = 0; i < n; ++i)
        out.mMinNodal = std::min(out.mMinNodal, curved_triangle_jacobian(Order, pNodes, rNormal,
                                                                         r_ref[i][0], r_ref[i][1]));
    if (!(out.mMinNodal > 0.0)) {
        out.mValid = false;
        return out;
    }
    if (Order == 1) {
        out.mValid = true;  // J is constant
        return out;
    }
    const int v = ct_judge(Order, pNodes, rNormal, {0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}, 0, MaxDepth);
    out.mValid = v == 1;
    out.mDecided = v >= 0;
    return out;
}

}  // namespace meshioplusplus::detail

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

#pragma once

/**
 * @file formats/lagrange_common.hpp
 * @brief Arbitrary-order Lagrange cells: the `VTK_LAGRANGE_*` node layout and
 *        nodal interpolation between point sets of one polynomial space.
 *
 * - `vtk_lattice` lists the nodes of a VTK Lagrange cell of order `p` as
 *   integer lattice coordinates (`x = i/p` ...) on the reference cell whose
 *   corners are meshio++'s (and MFEM's): the unit segment, the triangle
 *   (0,0) (1,0) (0,1), the unit square, the tetrahedron at the origin and the
 *   unit axes, the unit cube, and the prism of that triangle over [0, 1]. The
 *   order follows VTK: corners, edges, faces, interior (the hexahedron's
 *   vertical edges in VTK 9.1's order). Ported from MFEM's `mesh/vtk.cpp`
 *   (BSD-3-Clause, Lawrence Livermore National Security), which writes the
 *   files ParaView reads.
 * - `interpolation_matrix` maps values at one set of nodes of the cell's
 *   polynomial space (P_p on simplices, Q_p on the square and cube, P_p x P_p
 *   on the prism) to values at other points, through an orthogonal basis
 *   (Legendre products, Dubiner's collapsed bases) so the Vandermonde solve is
 *   well conditioned.
 * - `lattice_weights` names a lattice node by the corners it is a combination
 *   of, with integer weights on a common denominator `p^3`, so the same node
 *   of two cells sharing an edge or face gets the same key.
 *
 * A **format-private** header, beside the `.cpp` files that use it (MFEM
 * today), in a named namespace for the amalgamation. Twin of
 * `src/python/meshioplusplus/_lagrange.py`.
 */

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

namespace meshioplusplus {
namespace lagrange {

enum class Shape { Line, Triangle, Quad, Tetra, Hexahedron, Wedge };

inline int shape_dim(Shape S) {
    switch (S) {
        case Shape::Line:
            return 1;
        case Shape::Triangle:
        case Shape::Quad:
            return 2;
        default:
            return 3;
    }
}

/// Number of nodes of an order-`P` cell.
inline std::size_t num_nodes(Shape S, int P) {
    const std::size_t p = static_cast<std::size_t>(P);
    switch (S) {
        case Shape::Line:
            return p + 1;
        case Shape::Triangle:
            return (p + 1) * (p + 2) / 2;
        case Shape::Quad:
            return (p + 1) * (p + 1);
        case Shape::Tetra:
            return (p + 1) * (p + 2) * (p + 3) / 6;
        case Shape::Hexahedron:
            return (p + 1) * (p + 1) * (p + 1);
        case Shape::Wedge:
            return (p + 1) * (p + 1) * (p + 2) / 2;
    }
    return 0;
}

/// The Gauss-Lobatto points of order `P` on [0, 1], ascending.
inline std::vector<double> gll_points(int P) {
    std::vector<double> x(static_cast<std::size_t>(P) + 1);
    x.front() = 0.0;
    x.back() = 1.0;
    // Interior points: the roots of P'_P on [-1, 1], by Newton from the
    // Chebyshev-Gauss-Lobatto guesses; computed for the lower half and
    // mirrored, so the set is exactly symmetric.
    for (int k = 1; 2 * k <= P; ++k) {
        long double t = -std::cos(3.14159265358979323846264338327950288L * k / P);
        for (int it = 0; it < 100; ++it) {
            // Legendre P_P(t), P'_P(t) by recurrence.
            long double p0 = 1.0L, p1 = t;
            for (int n = 1; n < P; ++n) {
                const long double p2 = ((2 * n + 1) * t * p1 - n * p0) / (n + 1);
                p0 = p1;
                p1 = p2;
            }
            const long double d1 = P * (t * p1 - p0) / (t * t - 1.0L);
            const long double d2 = (2.0L * t * d1 - P * (P + 1) * p1) / (1.0L - t * t);
            const long double step = d1 / d2;
            t -= step;
            if (std::fabs(static_cast<double>(step)) < 1e-19)
                break;
        }
        const double v = static_cast<double>((1.0L + t) / 2.0L);
        x[static_cast<std::size_t>(k)] = v;
        x[static_cast<std::size_t>(P - k)] = 1.0 - v;
    }
    if (P % 2 == 0 && P > 0)
        x[static_cast<std::size_t>(P / 2)] = 0.5;
    return x;
}

/// Equispaced points `i/P` on [0, 1].
inline std::vector<double> uniform_points(int P) {
    std::vector<double> x(static_cast<std::size_t>(P) + 1);
    for (int i = 0; i <= P; ++i)
        x[static_cast<std::size_t>(i)] = static_cast<double>(i) / P;
    return x;
}

// --- VTK node order (after MFEM's mesh/vtk.cpp) ---------------------------------------

inline int vtk_triangle_index(int* pB, int Ref) {
    int max = Ref, min = 0;
    const int bmin = std::min(std::min(pB[0], pB[1]), pB[2]);
    int idx = 0;
    while (bmin > min) {
        idx += 3 * Ref;
        max -= 2;
        ++min;
        Ref -= 3;
    }
    for (int d = 0; d < 3; ++d) {
        if (pB[(d + 2) % 3] == max)
            return idx;
        ++idx;
    }
    for (int d = 0; d < 3; ++d) {
        if (pB[(d + 1) % 3] == min)
            return idx + pB[d] - (min + 1);
        idx += max - (min + 1);
    }
    return idx;
}

inline int vtk_tetra_index(int* pB, int Ref) {
    int idx = 0, max = Ref, min = 0;
    const int bmin = std::min(std::min(std::min(pB[0], pB[1]), pB[2]), pB[3]);
    while (bmin > min) {
        idx += 2 * (Ref * Ref + 1);
        max -= 3;
        ++min;
        Ref -= 4;
    }
    static const int kVertexMax[4] = {3, 0, 1, 2};
    static const int kEdgeMin[6][2] = {{1, 2}, {2, 3}, {0, 2}, {0, 1}, {1, 3}, {0, 3}};
    static const int kEdgeCount[6] = {0, 1, 3, 2, 2, 2};
    static const int kFaceMin[4] = {1, 3, 0, 2};
    static const int kFaceB[4][3] = {{0, 2, 3}, {2, 0, 1}, {2, 1, 3}, {1, 0, 3}};
    for (int v = 0; v < 4; ++v) {
        if (pB[kVertexMax[v]] == max)
            return idx;
        ++idx;
    }
    for (int e = 0; e < 6; ++e) {
        if (pB[kEdgeMin[e][0]] == min && pB[kEdgeMin[e][1]] == min)
            return idx + pB[kEdgeCount[e]] - (min + 1);
        idx += max - (min + 1);
    }
    for (int f = 0; f < 4; ++f) {
        if (pB[kFaceMin[f]] == min) {
            int projected[3];
            for (int i = 0; i < 3; ++i)
                projected[i] = pB[kFaceB[f][i]] - min;
            return idx + vtk_triangle_index(projected, Ref) - 3 * Ref;
        }
        idx += (Ref + 1) * (Ref + 2) / 2 - 3 * Ref;
    }
    return idx;
}

inline int vtk_triangle_offset(int Ref, int I, int J) {
    return I + Ref * (J - 1) - (J * (J + 1)) / 2;
}

inline int vtk_wedge_index(int I, int J, int K, int Ref) {
    const int om1 = Ref - 1;
    const int ibdr = I == 0, jbdr = J == 0, ijbdr = I + J == Ref, kbdr = K == 0 || K == Ref;
    const int nbdr = ibdr + jbdr + ijbdr + kbdr;
    if (nbdr == 3)
        return (ibdr && jbdr ? 0 : (jbdr && ijbdr ? 1 : 2)) + (K ? 3 : 0);
    int offset = 6;
    if (nbdr == 2) {
        if (!kbdr) {
            offset += om1 * 6;
            return offset + (K - 1) + ((ibdr && jbdr) ? 0 : (jbdr && ijbdr ? 1 : 2)) * om1;
        }
        offset += K == Ref ? 3 * om1 : 0;
        if (jbdr)
            return offset + I - 1;
        offset += om1;
        if (ijbdr)
            return offset + J - 1;
        offset += om1;
        return offset + (Ref - J - 1);
    }
    offset += 9 * om1;
    const int ntf = (om1 - 1) * om1 / 2;
    const int nqf = om1 * om1;
    if (nbdr == 1) {
        if (kbdr) {
            if (K > 0)
                offset += ntf;
            return offset + vtk_triangle_offset(Ref, I, J);
        }
        offset += 2 * ntf;
        if (jbdr)
            return offset + (I - 1) + om1 * (K - 1);
        offset += nqf;
        if (ijbdr)
            return offset + (J - 1) + om1 * (K - 1);
        offset += nqf;
        return offset + (Ref - J - 1) + om1 * (K - 1);
    }
    offset += 2 * ntf + 3 * nqf;
    return offset + vtk_triangle_offset(Ref, I, J) + ntf * (K - 1);
}

inline int vtk_tensor_index(int Idx, int Ref, int Dim) {
    const int n = Ref + 1;
    if (Dim == 1) {
        if (Idx == 0 || Idx == Ref)
            return Idx ? 1 : 0;
        return Idx + 1;
    }
    if (Dim == 2) {
        const int i = Idx % n, j = Idx / n;
        const bool ibdr = i == 0 || i == Ref, jbdr = j == 0 || j == Ref;
        if (ibdr && jbdr)
            return i ? (j ? 2 : 1) : (j ? 3 : 0);
        int offset = 4;
        if (jbdr)
            return (i - 1) + (j ? Ref - 1 + Ref - 1 : 0) + offset;
        if (ibdr)
            return (j - 1) + (i ? Ref - 1 : 2 * (Ref - 1) + Ref - 1) + offset;
        offset += 2 * (Ref - 1 + Ref - 1);
        return offset + (i - 1) + (Ref - 1) * (j - 1);
    }
    const int i = Idx % n, j = (Idx / n) % n, k = Idx / (n * n);
    const bool ibdr = i == 0 || i == Ref, jbdr = j == 0 || j == Ref, kbdr = k == 0 || k == Ref;
    const int nbdr = (ibdr ? 1 : 0) + (jbdr ? 1 : 0) + (kbdr ? 1 : 0);
    if (nbdr == 3)
        return (i ? (j ? 2 : 1) : (j ? 3 : 0)) + (k ? 4 : 0);
    int offset = 8;
    if (nbdr == 2) {
        if (!ibdr)
            return (i - 1) + (j ? Ref - 1 + Ref - 1 : 0) + (k ? 2 * (Ref - 1 + Ref - 1) : 0) +
                   offset;
        if (!jbdr)
            return (j - 1) + (i ? Ref - 1 : 2 * (Ref - 1) + Ref - 1) +
                   (k ? 2 * (Ref - 1 + Ref - 1) : 0) + offset;
        offset += 4 * (Ref - 1) + 4 * (Ref - 1);
        return (k - 1) + (Ref - 1) * (i ? (j ? 2 : 1) : (j ? 3 : 0)) + offset;
    }
    offset += 4 * (Ref - 1 + Ref - 1 + Ref - 1);
    if (nbdr == 1) {
        if (ibdr)
            return (j - 1) + (Ref - 1) * (k - 1) + (i ? (Ref - 1) * (Ref - 1) : 0) + offset;
        offset += 2 * (Ref - 1) * (Ref - 1);
        if (jbdr)
            return (i - 1) + (Ref - 1) * (k - 1) + (j ? (Ref - 1) * (Ref - 1) : 0) + offset;
        offset += 2 * (Ref - 1) * (Ref - 1);
        return (i - 1) + (Ref - 1) * (j - 1) + (k ? (Ref - 1) * (Ref - 1) : 0) + offset;
    }
    offset += 2 * 3 * (Ref - 1) * (Ref - 1);
    return offset + (i - 1) + (Ref - 1) * ((j - 1) + (Ref - 1) * (k - 1));
}

/// The nodes of an order-`P` VTK Lagrange cell, in VTK order, as integer
/// lattice coordinates (`x = i / P`).
inline std::vector<std::array<int, 3>> vtk_lattice(Shape S, int P) {
    std::vector<std::array<int, 3>> out(num_nodes(S, P));
    switch (S) {
        case Shape::Triangle: {
            int b[3];
            for (b[1] = 0; b[1] <= P; ++b[1])
                for (b[0] = 0; b[0] <= P - b[1]; ++b[0]) {
                    b[2] = P - b[0] - b[1];
                    out[static_cast<std::size_t>(vtk_triangle_index(b, P))] = {b[0], b[1], 0};
                }
            break;
        }
        case Shape::Tetra: {
            int b[4];
            for (b[2] = 0; b[2] <= P; ++b[2])
                for (b[1] = 0; b[1] <= P - b[2]; ++b[1])
                    for (b[0] = 0; b[0] <= P - b[1] - b[2]; ++b[0]) {
                        b[3] = P - b[0] - b[1] - b[2];
                        out[static_cast<std::size_t>(vtk_tetra_index(b, P))] = {b[0], b[1], b[2]};
                    }
            break;
        }
        case Shape::Wedge:
            for (int k = 0; k <= P; ++k)
                for (int j = 0; j <= P; ++j)
                    for (int i = 0; i <= P - j; ++i)
                        out[static_cast<std::size_t>(vtk_wedge_index(i, j, k, P))] = {i, j, k};
            break;
        default: {
            const int dim = shape_dim(S);
            const int n = P + 1;
            const int total = dim == 1 ? n : (dim == 2 ? n * n : n * n * n);
            for (int idx = 0; idx < total; ++idx) {
                const std::array<int, 3> ijk = {idx % n, dim > 1 ? (idx / n) % n : 0,
                                                dim > 2 ? idx / (n * n) : 0};
                out[static_cast<std::size_t>(vtk_tensor_index(idx, P, dim))] = ijk;
            }
            break;
        }
    }
    return out;
}

/// Corner weights of lattice node `rIjk` of an order-`P` cell, over the common
/// denominator `P^3`: `(corner, weight)` for every corner with a nonzero one,
/// by corner. Barycentric on simplices, multilinear on the square and cube,
/// their product on the prism -- so a node shared by two cells has the same
/// weights over the corners they share.
inline std::vector<std::pair<int, std::int64_t>> lattice_weights(Shape S, int P,
                                                                 const std::array<int, 3>& rIjk) {
    const std::int64_t p = P, i = rIjk[0], j = rIjk[1], k = rIjk[2];
    std::vector<std::int64_t> w;
    switch (S) {
        case Shape::Line:
            w = {(p - i) * p * p, i * p * p};
            break;
        case Shape::Triangle:
            w = {(p - i - j) * p * p, i * p * p, j * p * p};
            break;
        case Shape::Quad:
            w = {(p - i) * (p - j) * p, i * (p - j) * p, i * j * p, (p - i) * j * p};
            break;
        case Shape::Tetra:
            w = {(p - i - j - k) * p * p, i * p * p, j * p * p, k * p * p};
            break;
        case Shape::Hexahedron:
            w = {(p - i) * (p - j) * (p - k),
                 i * (p - j) * (p - k),
                 i * j * (p - k),
                 (p - i) * j * (p - k),
                 (p - i) * (p - j) * k,
                 i * (p - j) * k,
                 i * j * k,
                 (p - i) * j * k};
            break;
        case Shape::Wedge:
            w = {(p - i - j) * (p - k) * p, i * (p - k) * p, j * (p - k) * p,
                 (p - i - j) * k * p,       i * k * p,       j * k * p};
            break;
    }
    std::vector<std::pair<int, std::int64_t>> out;
    for (std::size_t c = 0; c < w.size(); ++c)
        if (w[c] != 0)
            out.emplace_back(static_cast<int>(c), w[c]);
    return out;
}

// --- orthogonal bases ---------------------------------------------------------------

/// Jacobi P_n^{(A,0)}(x) for n = 0..N.
inline std::vector<double> jacobi(int N, double A, double X) {
    std::vector<double> v(static_cast<std::size_t>(N) + 1);
    v[0] = 1.0;
    if (N >= 1)
        v[1] = ((A + 2.0) * X + A) / 2.0;
    for (int n = 1; n < N; ++n) {
        const double a1 = 2.0 * (n + 1) * (n + A + 1) * (2 * n + A);
        const double a2 = (2 * n + A + 1) * A * A;
        const double a3 = (2 * n + A) * (2 * n + A + 1) * (2 * n + A + 2);
        const double a4 = 2.0 * (n + A) * n * (2 * n + A + 2);
        v[static_cast<std::size_t>(n) + 1] = ((a2 + a3 * X) * v[static_cast<std::size_t>(n)] -
                                              a4 * v[static_cast<std::size_t>(n) - 1]) /
                                             a1;
    }
    return v;
}

/// The orthogonal basis of the cell's order-`P` space at reference point `rX`.
inline std::vector<double> basis(Shape S, int P, const std::array<double, 3>& rX) {
    std::vector<double> out;
    const double r = 2.0 * rX[0] - 1.0, s = 2.0 * rX[1] - 1.0, t = 2.0 * rX[2] - 1.0;
    switch (S) {
        case Shape::Line: {
            out = jacobi(P, 0.0, r);
            break;
        }
        case Shape::Quad: {
            const auto a = jacobi(P, 0.0, r), b = jacobi(P, 0.0, s);
            for (int j = 0; j <= P; ++j)
                for (int i = 0; i <= P; ++i)
                    out.push_back(a[static_cast<std::size_t>(i)] * b[static_cast<std::size_t>(j)]);
            break;
        }
        case Shape::Hexahedron: {
            const auto a = jacobi(P, 0.0, r), b = jacobi(P, 0.0, s), c = jacobi(P, 0.0, t);
            for (int k = 0; k <= P; ++k)
                for (int j = 0; j <= P; ++j)
                    for (int i = 0; i <= P; ++i)
                        out.push_back(a[static_cast<std::size_t>(i)] *
                                      b[static_cast<std::size_t>(j)] *
                                      c[static_cast<std::size_t>(k)]);
            break;
        }
        case Shape::Triangle:
        case Shape::Wedge: {
            const double a = std::fabs(1.0 - s) > 1e-14 ? 2.0 * (1.0 + r) / (1.0 - s) - 1.0 : -1.0;
            const auto pa = jacobi(P, 0.0, a);
            std::vector<double> tri;
            for (int i = 0; i <= P; ++i) {
                const auto pb = jacobi(P - i, 2.0 * i + 1.0, s);
                const double f = std::pow((1.0 - s) / 2.0, i);
                for (int j = 0; i + j <= P; ++j)
                    tri.push_back(pa[static_cast<std::size_t>(i)] * f *
                                  pb[static_cast<std::size_t>(j)]);
            }
            if (S == Shape::Triangle) {
                out = std::move(tri);
            } else {
                const auto c = jacobi(P, 0.0, t);
                for (int k = 0; k <= P; ++k)
                    for (double v : tri)
                        out.push_back(v * c[static_cast<std::size_t>(k)]);
            }
            break;
        }
        case Shape::Tetra: {
            const double a = std::fabs(s + t) > 1e-14 ? 2.0 * (1.0 + r) / (-s - t) - 1.0 : -1.0;
            const double b = std::fabs(1.0 - t) > 1e-14 ? 2.0 * (1.0 + s) / (1.0 - t) - 1.0 : -1.0;
            const auto pa = jacobi(P, 0.0, a);
            for (int i = 0; i <= P; ++i) {
                const auto pb = jacobi(P - i, 2.0 * i + 1.0, b);
                const double fb = std::pow((1.0 - b) / 2.0, i);
                for (int j = 0; i + j <= P; ++j) {
                    const auto pc = jacobi(P - i - j, 2.0 * (i + j) + 2.0, t);
                    const double fc = std::pow((1.0 - t) / 2.0, i + j);
                    for (int k = 0; i + j + k <= P; ++k)
                        out.push_back(pa[static_cast<std::size_t>(i)] * fb *
                                      pb[static_cast<std::size_t>(j)] * fc *
                                      pc[static_cast<std::size_t>(k)]);
                }
            }
            break;
        }
    }
    return out;
}

/// The matrix (row-major, targets x nodes) mapping values at `rNodes` -- a
/// unisolvent node set of the order-`P` space -- to values at `rTargets`.
inline std::vector<double> interpolation_matrix(
    Shape S, int P, const std::vector<std::array<double, 3>>& rNodes,
    const std::vector<std::array<double, 3>>& rTargets) {
    const std::size_t n = rNodes.size();
    if (n != num_nodes(S, P))
        throw std::runtime_error("lagrange: node count does not match the space");
    // A = V^T (n x n): A[b][m] = basis_b(node_m); solve A X = Phi^T for X
    // (n x targets); the matrix is X^T.
    std::vector<double> a(n * n);
    for (std::size_t m = 0; m < n; ++m) {
        const auto v = basis(S, P, rNodes[m]);
        for (std::size_t b = 0; b < n; ++b)
            a[b * n + m] = v[b];
    }
    const std::size_t nt = rTargets.size();
    std::vector<double> x(n * nt);
    for (std::size_t t = 0; t < nt; ++t) {
        const auto v = basis(S, P, rTargets[t]);
        for (std::size_t b = 0; b < n; ++b)
            x[b * nt + t] = v[b];
    }
    // LU with partial pivoting, applied to every right-hand side at once.
    for (std::size_t c = 0; c < n; ++c) {
        std::size_t piv = c;
        for (std::size_t r = c + 1; r < n; ++r)
            if (std::fabs(a[r * n + c]) > std::fabs(a[piv * n + c]))
                piv = r;
        if (std::fabs(a[piv * n + c]) < 1e-300)
            throw std::runtime_error("lagrange: the node set is not unisolvent");
        if (piv != c) {
            for (std::size_t k = 0; k < n; ++k)
                std::swap(a[c * n + k], a[piv * n + k]);
            for (std::size_t k = 0; k < nt; ++k)
                std::swap(x[c * nt + k], x[piv * nt + k]);
        }
        for (std::size_t r = c + 1; r < n; ++r) {
            const double f = a[r * n + c] / a[c * n + c];
            if (f == 0.0)
                continue;
            for (std::size_t k = c; k < n; ++k)
                a[r * n + k] -= f * a[c * n + k];
            for (std::size_t k = 0; k < nt; ++k)
                x[r * nt + k] -= f * x[c * nt + k];
        }
    }
    for (std::size_t c = n; c-- > 0;) {
        for (std::size_t k = 0; k < nt; ++k) {
            double v = x[c * nt + k];
            for (std::size_t q = c + 1; q < n; ++q)
                v -= a[c * n + q] * x[q * nt + k];
            x[c * nt + k] = v / a[c * n + c];
        }
    }
    std::vector<double> out(nt * n);
    for (std::size_t t = 0; t < nt; ++t)
        for (std::size_t m = 0; m < n; ++m)
            out[t * n + m] = x[m * nt + t];
    // Snap exact reproductions (a target that is a node) to 0/1.
    for (double& v : out)
        if (std::fabs(v) < 1e-14)
            v = 0.0;
        else if (std::fabs(v - 1.0) < 1e-14)
            v = 1.0;
    return out;
}

}  // namespace lagrange
}  // namespace meshioplusplus

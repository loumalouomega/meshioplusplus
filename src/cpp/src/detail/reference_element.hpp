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

#pragma once

/**
 * @file detail/reference_element.hpp
 * @brief Shared reference-element shape functions, derivatives, Jacobians,
 *        arbitrary-order Lagrange lattices, and Bernstein-Lagrange change of
 *        basis (roadmap 3.4.3 / 15.1 / 15.3).
 *
 * Moved from `formats/lagrange_common.hpp`; twin of
 * `src/python/meshioplusplus/_lagrange.py`.
 */

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/types.hpp"

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

// --- Gmsh's node order (twin of _lagrange.gmsh_lattice) -------------------------------------
//
// The fixed-order types beyond VTK's own (triangle10, quad16, tetra20, wedge40,
// hexahedron64, ...) are held in Gmsh's order: the corners, each edge's
// interior nodes in the direction listed, each face's interior nodes as an
// element of the same family (recursively), then the volume's likewise -- a
// prism's being the inner triangle's nodes, each with its column of nodes
// along the third direction in Gmsh's line order. Checked against Gmsh 4.15
// for every order it has (tests/python/test_lagrange_gmsh_order.py).

using LatticeNode = std::array<int, 3>;

inline LatticeNode lattice_add(LatticeNode A, const LatticeNode& rB, int Scale = 1) {
    for (int k = 0; k < 3; ++k)
        A[static_cast<std::size_t>(k)] += Scale * rB[static_cast<std::size_t>(k)];
    return A;
}

/// One lattice step from `rFrom` towards `rTo` on an order-`P` cell.
inline LatticeNode lattice_step(const LatticeNode& rFrom, const LatticeNode& rTo, int P) {
    LatticeNode d;
    for (std::size_t k = 0; k < 3; ++k)
        d[k] = (rTo[k] - rFrom[k]) / P;  // exact: corners are P steps apart
    return d;
}

inline const std::vector<std::pair<int, int>>& gmsh_edges(Shape S) {
    static const std::vector<std::pair<int, int>> line = {{0, 1}};
    static const std::vector<std::pair<int, int>> tri = {{0, 1}, {1, 2}, {2, 0}};
    static const std::vector<std::pair<int, int>> quad = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
    static const std::vector<std::pair<int, int>> tet = {{0, 1}, {1, 2}, {2, 0},
                                                         {3, 0}, {3, 2}, {3, 1}};
    static const std::vector<std::pair<int, int>> hex = {{0, 1}, {0, 3}, {0, 4}, {1, 2},
                                                         {1, 5}, {2, 3}, {2, 6}, {3, 7},
                                                         {4, 5}, {4, 7}, {5, 6}, {6, 7}};
    static const std::vector<std::pair<int, int>> wedge = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 4},
                                                           {2, 5}, {3, 4}, {3, 5}, {4, 5}};
    switch (S) {
        case Shape::Line:
            return line;
        case Shape::Triangle:
            return tri;
        case Shape::Quad:
            return quad;
        case Shape::Tetra:
            return tet;
        case Shape::Hexahedron:
            return hex;
        case Shape::Wedge:
            return wedge;
    }
    return line;
}

inline const std::vector<std::vector<int>>& gmsh_faces(Shape S) {
    static const std::vector<std::vector<int>> none;
    static const std::vector<std::vector<int>> tet = {{0, 2, 1}, {0, 1, 3}, {0, 3, 2}, {3, 1, 2}};
    static const std::vector<std::vector<int>> hex = {{0, 3, 2, 1}, {0, 1, 5, 4}, {0, 4, 7, 3},
                                                      {1, 2, 6, 5}, {2, 3, 7, 6}, {4, 5, 6, 7}};
    static const std::vector<std::vector<int>> wedge = {
        {0, 2, 1}, {3, 4, 5}, {0, 1, 4, 3}, {0, 3, 5, 2}, {1, 2, 5, 4}};
    switch (S) {
        case Shape::Tetra:
            return tet;
        case Shape::Hexahedron:
            return hex;
        case Shape::Wedge:
            return wedge;
        default:
            return none;
    }
}

/// Each corner moved one lattice step towards each of its `rAround` corners.
inline std::vector<LatticeNode> lattice_inset(const std::vector<LatticeNode>& rCorners, int P,
                                              const std::vector<std::vector<int>>& rAround) {
    std::vector<LatticeNode> out;
    for (std::size_t k = 0; k < rCorners.size(); ++k) {
        LatticeNode c = rCorners[k];
        for (const int m : rAround[k])
            c = lattice_add(c, lattice_step(rCorners[k], rCorners[static_cast<std::size_t>(m)], P));
        out.push_back(c);
    }
    return out;
}

/// The two neighbours of each corner of an `N`-gon.
inline std::vector<std::vector<int>> lattice_ring(int N) {
    std::vector<std::vector<int>> out;
    for (int k = 0; k < N; ++k)
        out.push_back({(k + 1) % N, (k + N - 1) % N});
    return out;
}

inline void gmsh_nodes(Shape S, const std::vector<LatticeNode>& rC, int P,
                       std::vector<LatticeNode>& rOut) {
    if (P == 0) {
        rOut.push_back(rC[0]);
        return;
    }
    rOut.insert(rOut.end(), rC.begin(), rC.end());
    for (const auto& [a, b] : gmsh_edges(S)) {
        const auto& ca = rC[static_cast<std::size_t>(a)];
        const LatticeNode d = lattice_step(ca, rC[static_cast<std::size_t>(b)], P);
        for (int t = 1; t < P; ++t)
            rOut.push_back(lattice_add(ca, d, t));
    }
    for (const auto& face : gmsh_faces(S)) {
        std::vector<LatticeNode> v;
        for (const int i : face)
            v.push_back(rC[static_cast<std::size_t>(i)]);
        const int n = static_cast<int>(face.size());
        const int order = n == 3 ? P - 3 : P - 2;
        if (order >= 0)
            gmsh_nodes(n == 3 ? Shape::Triangle : Shape::Quad, lattice_inset(v, P, lattice_ring(n)),
                       order, rOut);
    }
    switch (S) {
        case Shape::Triangle:
            if (P >= 3)
                gmsh_nodes(S, lattice_inset(rC, P, lattice_ring(3)), P - 3, rOut);
            break;
        case Shape::Quad:
            if (P >= 2)
                gmsh_nodes(S, lattice_inset(rC, P, lattice_ring(4)), P - 2, rOut);
            break;
        case Shape::Tetra:
            if (P >= 4)
                gmsh_nodes(S, lattice_inset(rC, P, {{1, 2, 3}, {0, 2, 3}, {0, 1, 3}, {0, 1, 2}}),
                           P - 4, rOut);
            break;
        case Shape::Hexahedron:
            if (P >= 2)
                gmsh_nodes(S,
                           lattice_inset(rC, P,
                                         {{1, 3, 4},
                                          {0, 2, 5},
                                          {1, 3, 6},
                                          {0, 2, 7},
                                          {0, 5, 7},
                                          {1, 4, 6},
                                          {2, 5, 7},
                                          {3, 4, 6}}),
                           P - 2, rOut);
            break;
        case Shape::Wedge:
            if (P >= 3) {
                std::vector<LatticeNode> flat = {
                    {rC[0][0], rC[0][1], 0}, {rC[1][0], rC[1][1], 0}, {rC[2][0], rC[2][1], 0}};
                const LatticeNode dz = lattice_step(rC[0], rC[3], P);
                std::vector<LatticeNode> column, tri;
                gmsh_nodes(Shape::Line, {lattice_add(rC[0], dz), lattice_add(rC[3], dz, -1)}, P - 2,
                           column);
                gmsh_nodes(Shape::Triangle, lattice_inset(flat, P, lattice_ring(3)), P - 3, tri);
                for (const auto& t : tri)
                    for (const auto& z : column)
                        rOut.push_back({t[0], t[1], z[2]});
            }
            break;
        default:
            break;
    }
}

/// The nodes of an order-`P` cell in Gmsh's order, in `vtk_lattice`'s frame.
inline std::vector<LatticeNode> gmsh_lattice(Shape S, int P) {
    std::vector<LatticeNode> c;
    switch (S) {
        case Shape::Line:
            c = {{0, 0, 0}, {P, 0, 0}};
            break;
        case Shape::Triangle:
            c = {{0, 0, 0}, {P, 0, 0}, {0, P, 0}};
            break;
        case Shape::Quad:
            c = {{0, 0, 0}, {P, 0, 0}, {P, P, 0}, {0, P, 0}};
            break;
        case Shape::Tetra:
            c = {{0, 0, 0}, {P, 0, 0}, {0, P, 0}, {0, 0, P}};
            break;
        case Shape::Hexahedron:
            c = {{0, 0, 0}, {P, 0, 0}, {P, P, 0}, {0, P, 0},
                 {0, 0, P}, {P, 0, P}, {P, P, P}, {0, P, P}};
            break;
        case Shape::Wedge:
            c = {{0, 0, 0}, {P, 0, 0}, {0, P, 0}, {0, 0, P}, {P, 0, P}, {0, P, P}};
            break;
    }
    std::vector<LatticeNode> out;
    gmsh_nodes(S, c, P, out);
    return out;
}

/// `perm` with VTK Lagrange node `j` = Gmsh node `perm[j]` (twin of
/// `_lagrange.gmsh_to_vtk_order`).
inline std::vector<int> gmsh_to_vtk_order(Shape S, int P) {
    const auto gmsh = gmsh_lattice(S, P);
    const auto vtk = vtk_lattice(S, P);
    std::vector<int> perm(vtk.size());
    for (std::size_t j = 0; j < vtk.size(); ++j) {
        const auto it = std::find(gmsh.begin(), gmsh.end(), vtk[j]);
        perm[j] = static_cast<int>(it - gmsh.begin());
    }
    return perm;
}

/// A fixed-order type VTK has no id of its own for, written as a VTK
/// Lagrange cell (twin of `_vtk_common.FIXED_LAGRANGE`).
struct FixedLagrange {
    int mVtkType = 0;
    Shape mShape = Shape::Line;
    int mOrder = 0;
};

/// `rType` (`triangle10`, `hexahedron64`, ...) as a VTK Lagrange cell, for
/// orders 3 and up whose Gmsh order was checked against Gmsh: to 10, the
/// hexahedron and prism to 9.
inline bool fixed_lagrange(const std::string& rType, FixedLagrange& rOut) {
    static const std::pair<const char*, std::pair<Shape, int>> shapes[] = {
        {"line", {Shape::Line, 68}},
        {"triangle", {Shape::Triangle, 69}},
        {"quad", {Shape::Quad, 70}},
        {"tetra", {Shape::Tetra, 71}},
        {"hexahedron", {Shape::Hexahedron, 72}},
        {"wedge", {Shape::Wedge, 73}}};
    for (const auto& [name, info] : shapes) {
        const std::string prefix(name);
        if (rType.size() <= prefix.size() || rType.compare(0, prefix.size(), prefix) != 0)
            continue;
        const std::string digits = rType.substr(prefix.size());
        if (digits.find_first_not_of("0123456789") != std::string::npos || digits.size() > 4)
            return false;
        const std::size_t count = static_cast<std::size_t>(std::stoi(digits));
        const int max_order =
            info.first == Shape::Hexahedron || info.first == Shape::Wedge ? 9 : 10;
        for (int p = 3; p <= max_order; ++p)
            if (num_nodes(info.first, p) == count) {
                rOut = {info.second, info.first, p};
                return true;
            }
        return false;
    }
    return false;
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

// --- Bernstein basis & Lagrange <-> Bezier change of basis (roadmap 15.3) ---

/// Binomial coefficient C(n, k) as double.
inline double binomial(int N, int K) {
    if (K < 0 || K > N)
        return 0.0;
    if (K > N - K)
        K = N - K;
    double res = 1.0;
    for (int i = 1; i <= K; ++i)
        res = res * static_cast<double>(N - K + i) / static_cast<double>(i);
    return res;
}

/// 1-D Bernstein polynomials B_{i,P}(u) and their derivatives dB/du for i=0..P on [0, 1].
inline void bernstein1d_and_der(int P, double U, std::vector<double>& rVal,
                                std::vector<double>& rDer) {
    const std::size_t n = static_cast<std::size_t>(P) + 1;
    rVal.assign(n, 0.0);
    rDer.assign(n, 0.0);
    if (P == 0) {
        rVal[0] = 1.0;
        return;
    }
    const double om = 1.0 - U;
    std::vector<double> pu(n, 1.0), pom(n, 1.0);
    for (int i = 1; i <= P; ++i) {
        pu[static_cast<std::size_t>(i)] = pu[static_cast<std::size_t>(i - 1)] * U;
        pom[static_cast<std::size_t>(i)] = pom[static_cast<std::size_t>(i - 1)] * om;
    }
    for (int i = 0; i <= P; ++i) {
        const double c = binomial(P, i);
        rVal[static_cast<std::size_t>(i)] =
            c * pu[static_cast<std::size_t>(i)] * pom[static_cast<std::size_t>(P - i)];
        double d = 0.0;
        if (i > 0)
            d += static_cast<double>(i) * pu[static_cast<std::size_t>(i - 1)] *
                 pom[static_cast<std::size_t>(P - i)];
        if (i < P)
            d -= static_cast<double>(P - i) * pu[static_cast<std::size_t>(i)] *
                 pom[static_cast<std::size_t>(P - i - 1)];
        rDer[static_cast<std::size_t>(i)] = c * d;
    }
}

/// Evaluate Bernstein basis functions and exact reference derivatives in VTK lattice order.
inline void bernstein_basis_and_ders(Shape S, int P, const std::array<double, 3>& rX,
                                     std::vector<double>& rVal,
                                     std::vector<std::array<double, 3>>& rGrad) {
    const auto lattice = vtk_lattice(S, P);
    const std::size_t n = lattice.size();
    rVal.assign(n, 0.0);
    rGrad.assign(n, {0.0, 0.0, 0.0});
    if (P == 0) {
        if (n > 0)
            rVal[0] = 1.0;
        return;
    }
    const double u = rX[0], v = rX[1], w = rX[2];
    switch (S) {
        case Shape::Line: {
            std::vector<double> bu, dbu;
            bernstein1d_and_der(P, u, bu, dbu);
            for (std::size_t m = 0; m < n; ++m) {
                const auto i = static_cast<std::size_t>(lattice[m][0]);
                rVal[m] = bu[i];
                rGrad[m] = {dbu[i], 0.0, 0.0};
            }
            break;
        }
        case Shape::Quad: {
            std::vector<double> bu, dbu, bv, dbv;
            bernstein1d_and_der(P, u, bu, dbu);
            bernstein1d_and_der(P, v, bv, dbv);
            for (std::size_t m = 0; m < n; ++m) {
                const auto i = static_cast<std::size_t>(lattice[m][0]);
                const auto j = static_cast<std::size_t>(lattice[m][1]);
                rVal[m] = bu[i] * bv[j];
                rGrad[m] = {dbu[i] * bv[j], bu[i] * dbv[j], 0.0};
            }
            break;
        }
        case Shape::Hexahedron: {
            std::vector<double> bu, dbu, bv, dbv, bw, dbw;
            bernstein1d_and_der(P, u, bu, dbu);
            bernstein1d_and_der(P, v, bv, dbv);
            bernstein1d_and_der(P, w, bw, dbw);
            for (std::size_t m = 0; m < n; ++m) {
                const auto i = static_cast<std::size_t>(lattice[m][0]);
                const auto j = static_cast<std::size_t>(lattice[m][1]);
                const auto k = static_cast<std::size_t>(lattice[m][2]);
                rVal[m] = bu[i] * bv[j] * bw[k];
                rGrad[m] = {dbu[i] * bv[j] * bw[k], bu[i] * dbv[j] * bw[k],
                            bu[i] * bv[j] * dbw[k]};
            }
            break;
        }
        case Shape::Triangle: {
            const double l0 = 1.0 - u - v;
            std::vector<double> fact(static_cast<std::size_t>(P) + 1, 1.0);
            for (int k = 1; k <= P; ++k)
                fact[static_cast<std::size_t>(k)] =
                    fact[static_cast<std::size_t>(k - 1)] * static_cast<double>(k);
            const auto pow_int = [](double Base, int Exp) {
                double r = 1.0;
                for (int e = 0; e < Exp; ++e)
                    r *= Base;
                return r;
            };
            for (std::size_t m = 0; m < n; ++m) {
                const int i = lattice[m][0];
                const int j = lattice[m][1];
                const int k = P - i - j;
                const double c = fact[static_cast<std::size_t>(P)] /
                                 (fact[static_cast<std::size_t>(i)] *
                                  fact[static_cast<std::size_t>(j)] *
                                  fact[static_cast<std::size_t>(k)]);
                const double pu = pow_int(u, i);
                const double pv = pow_int(v, j);
                const double pl = pow_int(l0, k);
                rVal[m] = c * pu * pv * pl;
                double du = 0.0, dv = 0.0;
                if (i > 0)
                    du += static_cast<double>(i) * pow_int(u, i - 1) * pv * pl;
                if (k > 0)
                    du -= static_cast<double>(k) * pu * pv * pow_int(l0, k - 1);
                if (j > 0)
                    dv += static_cast<double>(j) * pu * pow_int(v, j - 1) * pl;
                if (k > 0)
                    dv -= static_cast<double>(k) * pu * pv * pow_int(l0, k - 1);
                rGrad[m] = {c * du, c * dv, 0.0};
            }
            break;
        }
        case Shape::Tetra: {
            const double l0 = 1.0 - u - v - w;
            std::vector<double> fact(static_cast<std::size_t>(P) + 1, 1.0);
            for (int k = 1; k <= P; ++k)
                fact[static_cast<std::size_t>(k)] =
                    fact[static_cast<std::size_t>(k - 1)] * static_cast<double>(k);
            const auto pow_int = [](double Base, int Exp) {
                double r = 1.0;
                for (int e = 0; e < Exp; ++e)
                    r *= Base;
                return r;
            };
            for (std::size_t m = 0; m < n; ++m) {
                const int i = lattice[m][0];
                const int j = lattice[m][1];
                const int k = lattice[m][2];
                const int l = P - i - j - k;
                const double c = fact[static_cast<std::size_t>(P)] /
                                 (fact[static_cast<std::size_t>(i)] *
                                  fact[static_cast<std::size_t>(j)] *
                                  fact[static_cast<std::size_t>(k)] *
                                  fact[static_cast<std::size_t>(l)]);
                const double pu = pow_int(u, i);
                const double pv = pow_int(v, j);
                const double pw = pow_int(w, k);
                const double pl = pow_int(l0, l);
                rVal[m] = c * pu * pv * pw * pl;
                double du = 0.0, dv = 0.0, dw = 0.0;
                if (i > 0)
                    du += static_cast<double>(i) * pow_int(u, i - 1) * pv * pw * pl;
                if (l > 0)
                    du -= static_cast<double>(l) * pu * pv * pw * pow_int(l0, l - 1);
                if (j > 0)
                    dv += static_cast<double>(j) * pu * pow_int(v, j - 1) * pw * pl;
                if (l > 0)
                    dv -= static_cast<double>(l) * pu * pv * pw * pow_int(l0, l - 1);
                if (k > 0)
                    dw += static_cast<double>(k) * pu * pv * pow_int(w, k - 1) * pl;
                if (l > 0)
                    dw -= static_cast<double>(l) * pu * pv * pw * pow_int(l0, l - 1);
                rGrad[m] = {c * du, c * dv, c * dw};
            }
            break;
        }
        case Shape::Wedge: {
            const double l0 = 1.0 - u - v;
            std::vector<double> fact(static_cast<std::size_t>(P) + 1, 1.0);
            for (int k = 1; k <= P; ++k)
                fact[static_cast<std::size_t>(k)] =
                    fact[static_cast<std::size_t>(k - 1)] * static_cast<double>(k);
            const auto pow_int = [](double Base, int Exp) {
                double r = 1.0;
                for (int e = 0; e < Exp; ++e)
                    r *= Base;
                return r;
            };
            std::vector<double> bw, dbw;
            bernstein1d_and_der(P, w, bw, dbw);
            for (std::size_t m = 0; m < n; ++m) {
                const int i = lattice[m][0];
                const int j = lattice[m][1];
                const int k = lattice[m][2];
                const int l = P - i - j;
                const double c = fact[static_cast<std::size_t>(P)] /
                                 (fact[static_cast<std::size_t>(i)] *
                                  fact[static_cast<std::size_t>(j)] *
                                  fact[static_cast<std::size_t>(l)]);
                const double pu = pow_int(u, i);
                const double pv = pow_int(v, j);
                const double pl = pow_int(l0, l);
                const double tri_v = c * pu * pv * pl;
                double du = 0.0, dv = 0.0;
                if (i > 0)
                    du += static_cast<double>(i) * pow_int(u, i - 1) * pv * pl;
                if (l > 0)
                    du -= static_cast<double>(l) * pu * pv * pow_int(l0, l - 1);
                if (j > 0)
                    dv += static_cast<double>(j) * pu * pow_int(v, j - 1) * pl;
                if (l > 0)
                    dv -= static_cast<double>(l) * pu * pv * pow_int(l0, l - 1);
                const auto kw = static_cast<std::size_t>(k);
                rVal[m] = tri_v * bw[kw];
                rGrad[m] = {c * du * bw[kw], c * dv * bw[kw], tri_v * dbw[kw]};
            }
            break;
        }
    }
}

/// Invert a small dense row-major N x N matrix via Gaussian elimination with partial pivoting.
inline std::vector<double> invert_small_matrix(std::vector<double> A, std::size_t N) {
    std::vector<double> inv(N * N, 0.0);
    for (std::size_t i = 0; i < N; ++i)
        inv[i * N + i] = 1.0;
    for (std::size_t c = 0; c < N; ++c) {
        std::size_t piv = c;
        for (std::size_t r = c + 1; r < N; ++r)
            if (std::fabs(A[r * N + c]) > std::fabs(A[piv * N + c]))
                piv = r;
        if (std::fabs(A[piv * N + c]) < 1e-300)
            throw std::runtime_error("reference_element: singular change-of-basis matrix");
        if (piv != c) {
            for (std::size_t k = 0; k < N; ++k) {
                std::swap(A[c * N + k], A[piv * N + k]);
                std::swap(inv[c * N + k], inv[piv * N + k]);
            }
        }
        const double diag = A[c * N + c];
        for (std::size_t k = 0; k < N; ++k) {
            A[c * N + k] /= diag;
            inv[c * N + k] /= diag;
        }
        for (std::size_t r = 0; r < N; ++r) {
            if (r == c)
                continue;
            const double f = A[r * N + c];
            if (f == 0.0)
                continue;
            for (std::size_t k = 0; k < N; ++k) {
                A[r * N + k] -= f * A[c * N + k];
                inv[r * N + k] -= f * inv[c * N + k];
            }
        }
    }
    return inv;
}

/// Row-major N x N matrix mapping Bezier control values -> VTK Lagrange nodal values.
inline std::vector<double> bezier_to_lagrange_matrix(Shape S, int P) {
    const auto lattice = vtk_lattice(S, P);
    const std::size_t n = lattice.size();
    std::vector<double> m(n * n, 0.0);
    std::vector<double> bv;
    std::vector<std::array<double, 3>> bg;
    const double inv_p = P > 0 ? 1.0 / static_cast<double>(P) : 0.0;
    for (std::size_t row = 0; row < n; ++row) {
        const std::array<double, 3> xi{lattice[row][0] * inv_p, lattice[row][1] * inv_p,
                                       lattice[row][2] * inv_p};
        bernstein_basis_and_ders(S, P, xi, bv, bg);
        for (std::size_t col = 0; col < n; ++col)
            m[row * n + col] = bv[col];
    }
    return m;
}

/// Row-major N x N matrix mapping VTK Lagrange nodal values -> Bezier control values.
inline std::vector<double> lagrange_to_bezier_matrix(Shape S, int P) {
    const std::size_t n = num_nodes(S, P);
    return invert_small_matrix(bezier_to_lagrange_matrix(S, P), n);
}

/// Evaluate VTK Lagrange shape functions and exact reference derivatives at `rX`.
inline void lagrange_shape_and_ders(Shape S, int P, const std::array<double, 3>& rX,
                                    std::vector<double>& rVal,
                                    std::vector<std::array<double, 3>>& rGrad) {
    const std::size_t n = num_nodes(S, P);
    std::vector<double> bv;
    std::vector<std::array<double, 3>> bg;
    bernstein_basis_and_ders(S, P, rX, bv, bg);
    const auto m_inv = lagrange_to_bezier_matrix(S, P);
    // f(x) = sum_b c_b B_b(x), c = M^{-1} f_nodal => N_j(x) = sum_b (M^{-1})_{b,j} B_b(x).
    rVal.assign(n, 0.0);
    rGrad.assign(n, {0.0, 0.0, 0.0});
    for (std::size_t b = 0; b < n; ++b) {
        const double bval = bv[b];
        const auto& bgrd = bg[b];
        for (std::size_t j = 0; j < n; ++j) {
            const double w = m_inv[b * n + j];
            rVal[j] += w * bval;
            rGrad[j][0] += w * bgrd[0];
            rGrad[j][1] += w * bgrd[1];
            rGrad[j][2] += w * bgrd[2];
        }
    }
}

}  // namespace lagrange

namespace detail {

/// Metadata describing a cell type's reference element for shape-function evaluation.
struct ReferenceElementInfo {
    lagrange::Shape mShape = lagrange::Shape::Line;
    int mOrder = 1;
    int mDim = 1;
    bool mIsSerendipity = false;
    bool mIsPyramid = false;
    bool mIsBezier = false;
    std::size_t mNumNodes = 0;
};

/// Infer uniform polynomial order P for a VTK Lagrange / Bezier cell from its node count.
inline int infer_uniform_order(lagrange::Shape S, std::size_t NumNodes) {
    for (int p = 1; p <= 10; ++p) {
        if (lagrange::num_nodes(S, p) == NumNodes)
            return p;
    }
    return 0;
}

/// Resolve a CellType + node count to its ReferenceElementInfo; returns false if unsupported.
inline bool reference_element_info(CellType Type, std::size_t NumNodes,
                                   ReferenceElementInfo& rOut) {
    rOut = ReferenceElementInfo{};
    rOut.mNumNodes = NumNodes;
    switch (Type) {
        case CellType::Line:
            rOut.mShape = lagrange::Shape::Line;
            rOut.mOrder = 1;
            rOut.mDim = 1;
            return NumNodes == 2;
        case CellType::Line3:
            rOut.mShape = lagrange::Shape::Line;
            rOut.mOrder = 2;
            rOut.mDim = 1;
            return NumNodes == 3;
        case CellType::Triangle:
            rOut.mShape = lagrange::Shape::Triangle;
            rOut.mOrder = 1;
            rOut.mDim = 2;
            return NumNodes == 3;
        case CellType::Triangle6:
            rOut.mShape = lagrange::Shape::Triangle;
            rOut.mOrder = 2;
            rOut.mDim = 2;
            return NumNodes == 6;
        case CellType::Quad:
            rOut.mShape = lagrange::Shape::Quad;
            rOut.mOrder = 1;
            rOut.mDim = 2;
            return NumNodes == 4;
        case CellType::Quad8:
            rOut.mShape = lagrange::Shape::Quad;
            rOut.mOrder = 2;
            rOut.mDim = 2;
            rOut.mIsSerendipity = true;
            return NumNodes == 8;
        case CellType::Quad9:
            rOut.mShape = lagrange::Shape::Quad;
            rOut.mOrder = 2;
            rOut.mDim = 2;
            return NumNodes == 9;
        case CellType::Tetra:
            rOut.mShape = lagrange::Shape::Tetra;
            rOut.mOrder = 1;
            rOut.mDim = 3;
            return NumNodes == 4;
        case CellType::Tetra10:
            rOut.mShape = lagrange::Shape::Tetra;
            rOut.mOrder = 2;
            rOut.mDim = 3;
            return NumNodes == 10;
        case CellType::Hexahedron:
            rOut.mShape = lagrange::Shape::Hexahedron;
            rOut.mOrder = 1;
            rOut.mDim = 3;
            return NumNodes == 8;
        case CellType::Hexahedron20:
            rOut.mShape = lagrange::Shape::Hexahedron;
            rOut.mOrder = 2;
            rOut.mDim = 3;
            rOut.mIsSerendipity = true;
            return NumNodes == 20;
        case CellType::Hexahedron27:
            rOut.mShape = lagrange::Shape::Hexahedron;
            rOut.mOrder = 2;
            rOut.mDim = 3;
            return NumNodes == 27;
        case CellType::Wedge:
            rOut.mShape = lagrange::Shape::Wedge;
            rOut.mOrder = 1;
            rOut.mDim = 3;
            return NumNodes == 6;
        case CellType::Pyramid:
            rOut.mShape = lagrange::Shape::Hexahedron;
            rOut.mOrder = 1;
            rOut.mDim = 3;
            rOut.mIsPyramid = true;
            return NumNodes == 5;
        case CellType::VtkLagrangeCurve:
        case CellType::VtkBezierCurve:
            rOut.mShape = lagrange::Shape::Line;
            rOut.mDim = 1;
            rOut.mIsBezier = (Type == CellType::VtkBezierCurve);
            rOut.mOrder = infer_uniform_order(rOut.mShape, NumNodes);
            return rOut.mOrder > 0;
        case CellType::VtkLagrangeTriangle:
        case CellType::VtkBezierTriangle:
            rOut.mShape = lagrange::Shape::Triangle;
            rOut.mDim = 2;
            rOut.mIsBezier = (Type == CellType::VtkBezierTriangle);
            rOut.mOrder = infer_uniform_order(rOut.mShape, NumNodes);
            return rOut.mOrder > 0;
        case CellType::VtkLagrangeQuadrilateral:
        case CellType::VtkBezierQuadrilateral:
            rOut.mShape = lagrange::Shape::Quad;
            rOut.mDim = 2;
            rOut.mIsBezier = (Type == CellType::VtkBezierQuadrilateral);
            rOut.mOrder = infer_uniform_order(rOut.mShape, NumNodes);
            return rOut.mOrder > 0;
        case CellType::VtkLagrangeTetrahedron:
        case CellType::VtkBezierTetrahedron:
            rOut.mShape = lagrange::Shape::Tetra;
            rOut.mDim = 3;
            rOut.mIsBezier = (Type == CellType::VtkBezierTetrahedron);
            rOut.mOrder = infer_uniform_order(rOut.mShape, NumNodes);
            return rOut.mOrder > 0;
        case CellType::VtkLagrangeHexahedron:
        case CellType::VtkBezierHexahedron:
            rOut.mShape = lagrange::Shape::Hexahedron;
            rOut.mDim = 3;
            rOut.mIsBezier = (Type == CellType::VtkBezierHexahedron);
            rOut.mOrder = infer_uniform_order(rOut.mShape, NumNodes);
            return rOut.mOrder > 0;
        case CellType::VtkLagrangeWedge:
        case CellType::VtkBezierWedge:
            rOut.mShape = lagrange::Shape::Wedge;
            rOut.mDim = 3;
            rOut.mIsBezier = (Type == CellType::VtkBezierWedge);
            rOut.mOrder = infer_uniform_order(rOut.mShape, NumNodes);
            return rOut.mOrder > 0;
        default: {
            lagrange::FixedLagrange fl{};
            if (lagrange::fixed_lagrange(cell_type_name(Type), fl) &&
                lagrange::num_nodes(fl.mShape, fl.mOrder) == NumNodes) {
                rOut.mShape = fl.mShape;
                rOut.mOrder = fl.mOrder;
                rOut.mDim = lagrange::shape_dim(fl.mShape);
                return true;
            }
            return false;
        }
    }
}

/**
 * @brief Evaluate shape functions `rVal` and reference derivatives `rGrad` at
 * `[0, 1]` reference point `rXi` for any supported cell type.
 */
inline bool eval_reference_shape(CellType Type, std::size_t NumNodes,
                                 const std::array<double, 3>& rXi, std::vector<double>& rVal,
                                 std::vector<std::array<double, 3>>& rGrad) {
    ReferenceElementInfo info;
    if (!reference_element_info(Type, NumNodes, info))
        return false;

    const double u = rXi[0], v = rXi[1], w = rXi[2];

    if (info.mIsPyramid) {
        // Linear 5-node pyramid on 0 <= w <= 1, 0 <= u,v <= 1-w.
        const double omw = 1.0 - w;
        const double denom = std::fabs(omw) > 1e-14 ? omw : 1e-14;
        const double a = u / denom;
        const double b = v / denom;
        rVal = {(1.0 - a) * (1.0 - b) * omw, a * (1.0 - b) * omw, a * b * omw,
                (1.0 - a) * b * omw, w};
        // Derivatives via quotient rule: N0 = (1-w-u)(1-w-v)/(1-w), etc.
        const double du0 = -(1.0 - b), dv0 = -(1.0 - a);
        const double du1 = (1.0 - b), dv1 = -a;
        const double du2 = b, dv2 = a;
        const double du3 = -b, dv3 = (1.0 - a);
        // dN/dw = d/dw [ f(u/omw, v/omw) * omw ] = -f + (a*df/da + b*df/db)
        const double dw0 = -(1.0 - a) * (1.0 - b) + a * (-(1.0 - b)) + b * (-(1.0 - a));
        const double dw1 = -a * (1.0 - b) + a * (1.0 - b) + b * (-a);
        const double dw2 = -a * b + a * b + b * a;
        const double dw3 = -(1.0 - a) * b + a * (-b) + b * (1.0 - a);
        rGrad = {{du0, dv0, dw0}, {du1, dv1, dw1}, {du2, dv2, dw2}, {du3, dv3, dw3}, {0, 0, 1}};
        return true;
    }

    if (Type == CellType::Quad8) {
        // Serendipity 8-node quad on [0, 1]^2 (mapped to r = 2u-1, s = 2v-1 in [-1, 1]^2).
        const double r = 2.0 * u - 1.0, s = 2.0 * v - 1.0;
        const double mr = 1.0 - r, pr = 1.0 + r, ms = 1.0 - s, ps = 1.0 + s;
        rVal = {0.25 * mr * ms * (-r - s - 1.0), 0.25 * pr * ms * (r - s - 1.0),
                0.25 * pr * ps * (r + s - 1.0),  0.25 * mr * ps * (-r + s - 1.0),
                0.5 * (1.0 - r * r) * ms,        0.5 * pr * (1.0 - s * s),
                0.5 * (1.0 - r * r) * ps,        0.5 * mr * (1.0 - s * s)};
        // dr/du = 2, ds/dv = 2
        rGrad = {{0.5 * ms * (2.0 * r + s), 0.5 * mr * (r + 2.0 * s), 0.0},
                 {0.5 * ms * (2.0 * r - s), 0.5 * pr * (-r + 2.0 * s), 0.0},
                 {0.5 * ps * (2.0 * r + s), 0.5 * pr * (r + 2.0 * s), 0.0},
                 {0.5 * ps * (2.0 * r - s), 0.5 * mr * (-r + 2.0 * s), 0.0},
                 {-2.0 * r * ms, -(1.0 - r * r), 0.0},
                 {(1.0 - s * s), -2.0 * pr * s, 0.0},
                 {-2.0 * r * ps, (1.0 - r * r), 0.0},
                 {-(1.0 - s * s), -2.0 * mr * s, 0.0}};
        return true;
    }

    if (Type == CellType::Hexahedron20) {
        // Serendipity 20-node hex on [0, 1]^3 (mapped to r=2u-1, s=2v-1, t=2w-1).
        const double r = 2.0 * u - 1.0, s = 2.0 * v - 1.0, t = 2.0 * w - 1.0;
        constexpr double kC[8][3] = {{-1, -1, -1}, {1, -1, -1}, {1, 1, -1}, {-1, 1, -1},
                                     {-1, -1, 1},  {1, -1, 1},  {1, 1, 1},  {-1, 1, 1}};
        // Edge midpoints in meshio/VTK order:
        // 8..11 bottom, 12..15 top, 16..19 vertical (0-4, 1-5, 2-6, 3-7).
        constexpr double kE[12][3] = {{0, -1, -1}, {1, 0, -1}, {0, 1, -1}, {-1, 0, -1},
                                      {0, -1, 1},  {1, 0, 1},  {0, 1, 1},  {-1, 0, 1},
                                      {-1, -1, 0}, {1, -1, 0}, {1, 1, 0},  {-1, 1, 0}};
        rVal.resize(20);
        rGrad.resize(20);
        for (int i = 0; i < 8; ++i) {
            const double ri = kC[i][0], si = kC[i][1], ti = kC[i][2];
            const double fr = 1.0 + ri * r, fs = 1.0 + si * s, ft = 1.0 + ti * t;
            const double g = ri * r + si * s + ti * t - 2.0;
            rVal[static_cast<std::size_t>(i)] = 0.125 * fr * fs * ft * g;
            // d/du = 2 d/dr
            rGrad[static_cast<std::size_t>(i)] = {
                0.25 * ri * fs * ft * (g + fr),
                0.25 * si * fr * ft * (g + fs),
                0.25 * ti * fr * fs * (g + ft),
            };
        }
        for (int i = 0; i < 12; ++i) {
            const double ri = kE[i][0], si = kE[i][1], ti = kE[i][2];
            double val = 0.0, dr = 0.0, ds = 0.0, dt = 0.0;
            if (ri == 0.0) {
                const double fs = 1.0 + si * s, ft = 1.0 + ti * t;
                val = 0.25 * (1.0 - r * r) * fs * ft;
                dr = -0.5 * r * fs * ft;
                ds = 0.25 * (1.0 - r * r) * si * ft;
                dt = 0.25 * (1.0 - r * r) * fs * ti;
            } else if (si == 0.0) {
                const double fr = 1.0 + ri * r, ft = 1.0 + ti * t;
                val = 0.25 * fr * (1.0 - s * s) * ft;
                dr = 0.25 * ri * (1.0 - s * s) * ft;
                ds = -0.5 * s * fr * ft;
                dt = 0.25 * fr * (1.0 - s * s) * ti;
            } else {
                const double fr = 1.0 + ri * r, fs = 1.0 + si * s;
                val = 0.25 * fr * fs * (1.0 - t * t);
                dr = 0.25 * ri * fs * (1.0 - t * t);
                ds = 0.25 * fr * si * (1.0 - t * t);
                dt = -0.5 * t * fr * fs;
            }
            rVal[static_cast<std::size_t>(8 + i)] = val;
            rGrad[static_cast<std::size_t>(8 + i)] = {2.0 * dr, 2.0 * ds, 2.0 * dt};
        }
        return true;
    }

    if (Type == CellType::Hexahedron27) {
        // Meshio/VTK Hexahedron27 has vertical edges 16..19 in cyclic order (0-4, 1-5, 2-6, 3-7),
        // whereas VTK_LAGRANGE_HEXAHEDRON (p=2) has nodes 18 and 19 swapped!
        lagrange::lagrange_shape_and_ders(lagrange::Shape::Hexahedron, 2, rXi, rVal, rGrad);
        std::swap(rVal[18], rVal[19]);
        std::swap(rGrad[18], rGrad[19]);
        return true;
    }

    if (info.mIsBezier) {
        lagrange::bernstein_basis_and_ders(info.mShape, info.mOrder, rXi, rVal, rGrad);
        return true;
    }

    lagrange::lagrange_shape_and_ders(info.mShape, info.mOrder, rXi, rVal, rGrad);
    return true;
}

/// Physical Jacobian measure `|det J|` (or surface area / curve length factor) and
/// physical shape gradients `nabla_x N_i` at reference point `rXi`.
inline bool shape_eval_and_jacobian(CellType Type, const std::vector<Vec3>& rCellPoints,
                                    const std::array<double, 3>& rXi, std::vector<double>& rVal,
                                    std::vector<Vec3>& rPhysGrad, double& rMeasure) {
    ReferenceElementInfo info;
    if (!reference_element_info(Type, rCellPoints.size(), info))
        return false;
    std::vector<std::array<double, 3>> dndxi;
    if (!eval_reference_shape(Type, rCellPoints.size(), rXi, rVal, dndxi))
        return false;

    const std::size_t n = rCellPoints.size();
    // Jacobian columns J[:, d] = sum_i x_i * dN_i / dxi_d
    Vec3 J[3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    for (std::size_t i = 0; i < n; ++i) {
        for (int d = 0; d < info.mDim; ++d) {
            const double g = dndxi[i][static_cast<std::size_t>(d)];
            J[d][0] += rCellPoints[i][0] * g;
            J[d][1] += rCellPoints[i][1] * g;
            J[d][2] += rCellPoints[i][2] * g;
        }
    }

    rPhysGrad.assign(n, {0.0, 0.0, 0.0});
    if (info.mDim == 1) {
        const double len2 = vec3_norm_sq(J[0]);
        rMeasure = std::sqrt(len2);
        if (!(len2 > 1e-30))
            return false;
        const Vec3 t = vec3_scale(J[0], 1.0 / len2);
        for (std::size_t i = 0; i < n; ++i)
            rPhysGrad[i] = vec3_scale(t, dndxi[i][0]);
        return true;
    }
    if (info.mDim == 2) {
        const Vec3 normal = vec3_cross(J[0], J[1]);
        const double area2 = vec3_norm_sq(normal);
        rMeasure = std::sqrt(area2);
        if (!(area2 > 1e-30))
            return false;
        // Metric tensor G = J^T J (2x2), det(G) = |J0 x J1|^2 = area2
        const double g00 = vec3_dot(J[0], J[0]);
        const double g01 = vec3_dot(J[0], J[1]);
        const double g11 = vec3_dot(J[1], J[1]);
        const double inv_det = 1.0 / area2;
        // Dual basis vectors e^0 = (g11 J0 - g01 J1) / det, e^1 = (-g01 J0 + g00 J1) / det
        const Vec3 e0 = vec3_scale(vec3_sub(vec3_scale(J[0], g11), vec3_scale(J[1], g01)), inv_det);
        const Vec3 e1 = vec3_scale(vec3_sub(vec3_scale(J[1], g00), vec3_scale(J[0], g01)), inv_det);
        for (std::size_t i = 0; i < n; ++i) {
            rPhysGrad[i] = vec3_add(vec3_scale(e0, dndxi[i][0]), vec3_scale(e1, dndxi[i][1]));
        }
        return true;
    }
    // 3-D volume cell: J = [J[0], J[1], J[2]], det(J) = triple_product(J[0], J[1], J[2]).
    const double det = triple_product(J[0], J[1], J[2]);
    rMeasure = std::abs(det);
    if (!(rMeasure > 1e-30))
        return false;
    const double inv_det = 1.0 / det;
    // Columns of J^{-T} are (J1 x J2)/det, (J2 x J0)/det, (J0 x J1)/det
    const Vec3 r0 = vec3_scale(vec3_cross(J[1], J[2]), inv_det);
    const Vec3 r1 = vec3_scale(vec3_cross(J[2], J[0]), inv_det);
    const Vec3 r2 = vec3_scale(vec3_cross(J[0], J[1]), inv_det);
    for (std::size_t i = 0; i < n; ++i) {
        rPhysGrad[i] = {r0[0] * dndxi[i][0] + r1[0] * dndxi[i][1] + r2[0] * dndxi[i][2],
                        r0[1] * dndxi[i][0] + r1[1] * dndxi[i][1] + r2[1] * dndxi[i][2],
                        r0[2] * dndxi[i][0] + r1[2] * dndxi[i][1] + r2[2] * dndxi[i][2]};
    }
    return true;
}

}  // namespace detail
}  // namespace meshioplusplus

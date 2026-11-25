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
// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/untangle.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/cell_faces.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/parallel.hpp"

// Project includes (private, not installed)
#include "../detail/typed_view.hpp"

namespace meshioplusplus {

namespace {

using detail::Vec3;

/// A corner simplex: the corner node, its three edge ends (right-handed for a
/// valid cell), and whether it is a tetrahedron (ideal: regular) or a
/// hexahedron corner (ideal: a cube corner).
struct UtgCorner {
    std::array<std::int64_t, 4> mNodes;
    bool mTet;
};

/// One corner incident to a vertex, with the 4-bit mask of local slots (0..3)
/// where that vertex appears in `UtgCorner::mNodes`.
struct UtgNodeCorner {
    std::uint32_t mCorner;
    std::uint8_t mMask;
};

// Knupp's hexahedron corners: corner c, then its three edge neighbours.
constexpr int kUtgHexCorner[8][4] = {{0, 1, 3, 4}, {1, 2, 0, 5}, {2, 3, 1, 6}, {3, 0, 2, 7},
                                     {4, 7, 5, 0}, {5, 4, 6, 1}, {6, 5, 7, 2}, {7, 6, 4, 3}};

// The inverse of the regular tetrahedron's edge matrix (columns: its three
// edges from corner 0, unit length), row-major.
constexpr double kUtgTetWinv[3][3] = {{1.0, -0.57735026918962573, -0.40824829046386302},
                                      {0.0, 1.1547005383792515, -0.40824829046386302},
                                      {0.0, 0.0, 1.2247448713915890}};

// Derivative weights alpha_m[c] such that dS_{r,c}/d(x_{n_m})_i = delta_{r,i} * inv_unit *
// alpha_m[c]: rows 1..3 are rows 0..2 of W^-1, row 0 is their negated sum.
constexpr double kUtgTetAlpha[4][3] = {
    {-(kUtgTetWinv[0][0] + kUtgTetWinv[1][0] + kUtgTetWinv[2][0]),
     -(kUtgTetWinv[0][1] + kUtgTetWinv[1][1] + kUtgTetWinv[2][1]),
     -(kUtgTetWinv[0][2] + kUtgTetWinv[1][2] + kUtgTetWinv[2][2])},
    {kUtgTetWinv[0][0], kUtgTetWinv[0][1], kUtgTetWinv[0][2]},
    {kUtgTetWinv[1][0], kUtgTetWinv[1][1], kUtgTetWinv[1][2]},
    {kUtgTetWinv[2][0], kUtgTetWinv[2][1], kUtgTetWinv[2][2]},
};

constexpr double kUtgHexAlpha[4][3] = {
    {-1.0, -1.0, -1.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};

struct UtgMeasure {
    double mDet;    ///< det S
    double mFrob2;  ///< |S|_F^2
};

inline void utg_shape_matrix(const std::vector<Vec3>& rX, const UtgCorner& rC, double InvScale,
                             double rS[3][3]) {
    const Vec3& p0 = rX[static_cast<std::size_t>(rC.mNodes[0])];
    const Vec3& p1 = rX[static_cast<std::size_t>(rC.mNodes[1])];
    const Vec3& p2 = rX[static_cast<std::size_t>(rC.mNodes[2])];
    const Vec3& p3 = rX[static_cast<std::size_t>(rC.mNodes[3])];
    const double a00 = (p1[0] - p0[0]) * InvScale;
    const double a10 = (p1[1] - p0[1]) * InvScale;
    const double a20 = (p1[2] - p0[2]) * InvScale;
    const double a01 = (p2[0] - p0[0]) * InvScale;
    const double a11 = (p2[1] - p0[1]) * InvScale;
    const double a21 = (p2[2] - p0[2]) * InvScale;
    const double a02 = (p3[0] - p0[0]) * InvScale;
    const double a12 = (p3[1] - p0[1]) * InvScale;
    const double a22 = (p3[2] - p0[2]) * InvScale;
    if (rC.mTet) {
        // W^-1 is upper triangular (w10 = w20 = w21 = 0).
        rS[0][0] = a00;
        rS[1][0] = a10;
        rS[2][0] = a20;
        rS[0][1] = a00 * kUtgTetWinv[0][1] + a01 * kUtgTetWinv[1][1];
        rS[1][1] = a10 * kUtgTetWinv[0][1] + a11 * kUtgTetWinv[1][1];
        rS[2][1] = a20 * kUtgTetWinv[0][1] + a21 * kUtgTetWinv[1][1];
        rS[0][2] = a00 * kUtgTetWinv[0][2] + a01 * kUtgTetWinv[1][2] + a02 * kUtgTetWinv[2][2];
        rS[1][2] = a10 * kUtgTetWinv[0][2] + a11 * kUtgTetWinv[1][2] + a12 * kUtgTetWinv[2][2];
        rS[2][2] = a20 * kUtgTetWinv[0][2] + a21 * kUtgTetWinv[1][2] + a22 * kUtgTetWinv[2][2];
    } else {
        rS[0][0] = a00;
        rS[0][1] = a01;
        rS[0][2] = a02;
        rS[1][0] = a10;
        rS[1][1] = a11;
        rS[1][2] = a12;
        rS[2][0] = a20;
        rS[2][1] = a21;
        rS[2][2] = a22;
    }
}

inline UtgMeasure utg_measure(const std::vector<Vec3>& rX, const UtgCorner& rC, double InvScale) {
    double s[3][3];
    utg_shape_matrix(rX, rC, InvScale, s);
    UtgMeasure m{};
    m.mDet = s[0][0] * (s[1][1] * s[2][2] - s[1][2] * s[2][1]) -
             s[0][1] * (s[1][0] * s[2][2] - s[1][2] * s[2][0]) +
             s[0][2] * (s[1][0] * s[2][1] - s[1][1] * s[2][0]);
    m.mFrob2 = s[0][0] * s[0][0] + s[0][1] * s[0][1] + s[0][2] * s[0][2] + s[1][0] * s[1][0] +
               s[1][1] * s[1][1] + s[1][2] * s[1][2] + s[2][0] * s[2][0] + s[2][1] * s[2][1] +
               s[2][2] * s[2][2];
    return m;
}

// det(kUtgTetWinv) = 1.0 * (2/sqrt(3)) * sqrt(3/2) = sqrt(2).
constexpr double kUtgTetDetWinv = kUtgTetWinv[0][0] * kUtgTetWinv[1][1] * kUtgTetWinv[2][2];

inline double utg_det(const std::vector<Vec3>& rX, const UtgCorner& rC, double InvScale3) {
    const Vec3& p0 = rX[static_cast<std::size_t>(rC.mNodes[0])];
    const Vec3& p1 = rX[static_cast<std::size_t>(rC.mNodes[1])];
    const Vec3& p2 = rX[static_cast<std::size_t>(rC.mNodes[2])];
    const Vec3& p3 = rX[static_cast<std::size_t>(rC.mNodes[3])];
    const double e00 = p1[0] - p0[0], e10 = p1[1] - p0[1], e20 = p1[2] - p0[2];
    const double e01 = p2[0] - p0[0], e11 = p2[1] - p0[1], e21 = p2[2] - p0[2];
    const double e02 = p3[0] - p0[0], e12 = p3[1] - p0[1], e22 = p3[2] - p0[2];
    const double det_a = e00 * (e11 * e22 - e12 * e21) - e01 * (e10 * e22 - e12 * e20) +
                         e02 * (e10 * e21 - e11 * e20);
    return det_a * (rC.mTet ? InvScale3 * kUtgTetDetWinv : InvScale3);
}

inline bool utg_corner_positive(const std::vector<Vec3>& rX, const UtgCorner& rC) {
    const Vec3& p0 = rX[static_cast<std::size_t>(rC.mNodes[0])];
    const Vec3& p1 = rX[static_cast<std::size_t>(rC.mNodes[1])];
    const Vec3& p2 = rX[static_cast<std::size_t>(rC.mNodes[2])];
    const Vec3& p3 = rX[static_cast<std::size_t>(rC.mNodes[3])];
    const double e00 = p1[0] - p0[0], e10 = p1[1] - p0[1], e20 = p1[2] - p0[2];
    const double e01 = p2[0] - p0[0], e11 = p2[1] - p0[1], e21 = p2[2] - p0[2];
    const double e02 = p3[0] - p0[0], e12 = p3[1] - p0[1], e22 = p3[2] - p0[2];
    const double det_a = e00 * (e11 * e22 - e12 * e21) - e01 * (e10 * e22 - e12 * e20) +
                         e02 * (e10 * e21 - e11 * e20);
    return det_a > 0.0;
}

/// The corner's scaled Jacobian: det of the edges over their lengths' product.
double utg_scaled_jacobian(const std::vector<Vec3>& rX, const UtgCorner& rC) {
    const Vec3& p0 = rX[static_cast<std::size_t>(rC.mNodes[0])];
    Vec3 e[3];
    double len = 1.0;
    for (int k = 0; k < 3; ++k) {
        const Vec3& p = rX[static_cast<std::size_t>(rC.mNodes[k + 1])];
        e[k] = {p[0] - p0[0], p[1] - p0[1], p[2] - p0[2]};
        len *= std::sqrt(e[k][0] * e[k][0] + e[k][1] * e[k][1] + e[k][2] * e[k][2]);
    }
    if (!(len > 0.0))
        return -1.0;
    const double det = e[0][0] * (e[1][1] * e[2][2] - e[1][2] * e[2][1]) -
                       e[1][0] * (e[0][1] * e[2][2] - e[0][2] * e[2][1]) +
                       e[2][0] * (e[0][1] * e[1][2] - e[0][2] * e[1][1]);
    return det / len;
}

// Fast full-precision cube root for strictly positive normal IEEE-754 doubles:
// bit-level initial guess (<= 3.2% error) + 2 cubic-convergent Halley steps
// (< 1 ULP, 53-bit exact) without glibc frexp/ldexp/scalbn overhead.
inline double utg_cbrt_pos(double X) {
    if (X < 1e-300 || !std::isfinite(X))
        return std::cbrt(X);
    std::uint64_t bits = 0;
    std::memcpy(&bits, &X, sizeof bits);
    bits = bits / 3 + 0x2a9f7893e6b00000ull;
    double y = 0.0;
    std::memcpy(&y, &bits, sizeof y);
    double y3 = y * y * y;
    y *= (y3 + 2.0 * X) / (2.0 * y3 + X);
    y3 = y * y * y;
    y *= (y3 + 2.0 * X) / (2.0 * y3 + X);
    return y;
}

inline double utg_term(const UtgMeasure& rM, double Delta) {
    const double h = 0.5 * (rM.mDet + std::sqrt(rM.mDet * rM.mDet + 4.0 * Delta * Delta));
    if (!(h > 0.0))
        return std::numeric_limits<double>::infinity();
    const double h13 = utg_cbrt_pos(h);
    return rM.mFrob2 / (3.0 * h13 * h13);
}

/// Closed-form Escobar-Knupp term and analytic gradient with respect to the
/// corner slots selected by `SlotMask`: `f = |S|_F^2 / (3 h^{2/3})`,
/// `d f / d S = (2 / (3 h^{2/3})) S - (2 f / (3 r)) cof(S)` with
/// `r = sqrt(det^2 + 4 delta^2)` and `h = (det + r) / 2`.
bool utg_term_and_grad(const std::vector<Vec3>& rX, const UtgCorner& rC, double InvScale,
                       double Delta, std::uint8_t SlotMask, double& rF, Vec3& rGrad) {
    double s[3][3];
    utg_shape_matrix(rX, rC, InvScale, s);
    const double c00 = s[1][1] * s[2][2] - s[1][2] * s[2][1];
    const double c01 = s[1][2] * s[2][0] - s[1][0] * s[2][2];
    const double c02 = s[1][0] * s[2][1] - s[1][1] * s[2][0];
    const double det = s[0][0] * c00 + s[0][1] * c01 + s[0][2] * c02;
    const double r = std::sqrt(det * det + 4.0 * Delta * Delta);
    const double h = 0.5 * (det + r);
    if (!(h > 0.0) || !(r > 0.0))
        return false;
    const double frob2 = s[0][0] * s[0][0] + s[0][1] * s[0][1] + s[0][2] * s[0][2] +
                         s[1][0] * s[1][0] + s[1][1] * s[1][1] + s[1][2] * s[1][2] +
                         s[2][0] * s[2][0] + s[2][1] * s[2][1] + s[2][2] * s[2][2];
    const double h13 = utg_cbrt_pos(h);
    const double inv_3h23 = 1.0 / (3.0 * h13 * h13);
    const double f = frob2 * inv_3h23;
    if (!std::isfinite(f))
        return false;
    rF += f;

    const double cof[3][3] = {
        {c00, c01, c02},
        {s[0][2] * s[2][1] - s[0][1] * s[2][2], s[0][0] * s[2][2] - s[0][2] * s[2][0],
         s[0][1] * s[2][0] - s[0][0] * s[2][1]},
        {s[0][1] * s[1][2] - s[0][2] * s[1][1], s[0][2] * s[1][0] - s[0][0] * s[1][2],
         s[0][0] * s[1][1] - s[0][1] * s[1][0]},
    };
    const double cs = 2.0 * inv_3h23 * InvScale;
    const double cd = (2.0 * f / (3.0 * r)) * InvScale;
    const auto& alpha = rC.mTet ? kUtgTetAlpha : kUtgHexAlpha;
    for (int m = 0; m < 4; ++m) {
        if ((SlotMask & (1u << m)) == 0)
            continue;
        const double* w = alpha[m];
        for (int row = 0; row < 3; ++row) {
            const double sw = s[row][0] * w[0] + s[row][1] * w[1] + s[row][2] * w[2];
            const double cw = cof[row][0] * w[0] + cof[row][1] * w[1] + cof[row][2] * w[2];
            rGrad[static_cast<std::size_t>(row)] += cs * sw - cd * cw;
        }
    }
    return true;
}

}  // namespace

UntangleResult untangle(const Mesh& rMesh, const UntangleOptions& rOptions) {
    if (rOptions.mMaxSweeps < 1)
        throw std::invalid_argument("untangle: max_sweeps must be at least 1");
    if (rOptions.mRings < 0)
        throw std::invalid_argument("untangle: rings must not be negative");
    const std::size_t n = rMesh.NumPoints();
    const std::size_t dim = rMesh.PointDim();
    if (!rOptions.mFrozen.empty() && rOptions.mFrozen.size() != n)
        throw std::invalid_argument("untangle: the frozen mask must have one entry per point");

    // Cells, their corners, and the nodes no measure covers.
    std::vector<UtgCorner> corners;
    std::vector<std::int64_t> cell_first;  // per measured cell, its first corner
    std::vector<std::size_t> cell_node_deg(n, 0);
    std::vector<std::pair<std::int64_t, std::int64_t>> node_cell_pairs;
    std::vector<char> pinned(n, 0);
    std::vector<std::array<std::int64_t, 4>> face_keys;

    for (const auto cb : rMesh.CellRange()) {
        const std::string& type = cb.Type();
        const bool tet = type == "tetra", hex = type == "hexahedron";
        const CellType ct = cb.IsRagged() ? CellType::Custom : cell_type_from_name(type);
        const int d = cb.IsRagged() ? (cb.IsPolyhedron() ? 3 : 2) : cell_type_dimension(ct);
        if (d < 3)
            continue;
        const std::size_t nc = cb.NumCells();
        if (!tet && !hex) {
            // Not measured: its nodes stay where they are.
            if (cb.IsRagged()) {
                for (std::size_t c = 0; c < nc; ++c)
                    for (std::size_t f = 0; f < cb.NumFaces(c); ++f) {
                        const auto face = cb.Face(c, f);
                        for (std::size_t k = 0; k < face.second; ++k)
                            pinned[static_cast<std::size_t>(face.first[k])] = 1;
                    }
            } else {
                const detail::Int64View conn_v(cb.Conn());
                const std::size_t npc = cb.NodesPerCell();
                for (std::size_t c = 0; c < nc; ++c)
                    for (std::size_t k = 0; k < npc; ++k)
                        pinned[static_cast<std::size_t>(conn_v[c * npc + k])] = 1;
            }
            continue;
        }
        const detail::Int64View conn_v(cb.Conn());
        const std::int64_t* pconn = conn_v.Data();
        const std::size_t npc = cb.NodesPerCell();
        const std::vector<detail::CellFaceDef>& faces = detail::cell_faces(ct);
        corners.reserve(corners.size() + nc * (tet ? 1 : 8));
        cell_first.reserve(cell_first.size() + nc + 1);
        node_cell_pairs.reserve(node_cell_pairs.size() + nc * npc);
        if (rOptions.mFixBoundary)
            face_keys.reserve(face_keys.size() + nc * faces.size());

        for (std::size_t c = 0; c < nc; ++c) {
            const std::int64_t* v = pconn + c * npc;
            const auto id = static_cast<std::int64_t>(cell_first.size());
            cell_first.push_back(static_cast<std::int64_t>(corners.size()));
            if (tet) {
                corners.push_back({{v[0], v[1], v[2], v[3]}, true});
            } else {
                for (const auto& row : kUtgHexCorner)
                    corners.push_back(
                        {{v[static_cast<std::size_t>(row[0])], v[static_cast<std::size_t>(row[1])],
                          v[static_cast<std::size_t>(row[2])], v[static_cast<std::size_t>(row[3])]},
                         false});
            }
            for (std::size_t k = 0; k < npc; ++k) {
                ++cell_node_deg[static_cast<std::size_t>(v[k])];
                node_cell_pairs.push_back({v[k], id});
            }
            if (rOptions.mFixBoundary) {
                for (const auto& f : faces) {
                    std::array<std::int64_t, 4> key{-1, -1, -1, -1};
                    for (std::size_t k = 0; k < f.mNumCorners; ++k)
                        key[k] = v[f.mNodes[k]];
                    std::sort(key.begin(), key.begin() + f.mNumCorners);
                    face_keys.push_back(key);
                }
            }
        }
    }
    cell_first.push_back(static_cast<std::int64_t>(corners.size()));
    const std::size_t num_cells = cell_first.size() - 1;
    if (num_cells == 0)
        throw std::invalid_argument("untangle: the mesh has no tetrahedron or hexahedron");

    // Flat CSR: node -> unique incident measured cells.
    std::vector<std::size_t> node_cell_start(n + 1, 0);
    for (std::size_t i = 0; i < n; ++i)
        node_cell_start[i + 1] = node_cell_start[i] + cell_node_deg[i];
    std::vector<std::int64_t> node_cell_ids(node_cell_start[n]);
    {
        std::vector<std::size_t> cursor = node_cell_start;
        for (const auto& pair : node_cell_pairs)
            node_cell_ids[cursor[static_cast<std::size_t>(pair.first)]++] = pair.second;
        std::vector<std::pair<std::int64_t, std::int64_t>>().swap(node_cell_pairs);
        std::vector<std::size_t> new_start(n + 1, 0);
        std::size_t out_pos = 0;
        for (std::size_t i = 0; i < n; ++i) {
            const std::size_t lo = node_cell_start[i], hi = node_cell_start[i + 1];
            new_start[i] = out_pos;
            for (std::size_t j = lo; j < hi; ++j) {
                if (j == lo || node_cell_ids[j] != node_cell_ids[j - 1])
                    node_cell_ids[out_pos++] = node_cell_ids[j];
            }
        }
        new_start[n] = out_pos;
        node_cell_ids.resize(out_pos);
        node_cell_start.swap(new_start);
    }

    // Flat CSR: node -> only the corners that actually contain that node
    // (for a hexahedron, 4 of its 8 corners instead of all 8).
    std::vector<std::size_t> node_corner_start(n + 1, 0);
    for (const UtgCorner& c : corners) {
        for (int m = 0; m < 4; ++m) {
            bool dup = false;
            for (int prev = 0; prev < m; ++prev)
                if (c.mNodes[static_cast<std::size_t>(prev)] ==
                    c.mNodes[static_cast<std::size_t>(m)])
                    dup = true;
            if (!dup)
                ++node_corner_start[static_cast<std::size_t>(
                                        c.mNodes[static_cast<std::size_t>(m)]) +
                                    1];
        }
    }
    for (std::size_t i = 0; i < n; ++i)
        node_corner_start[i + 1] += node_corner_start[i];
    std::vector<UtgNodeCorner> node_corners(node_corner_start[n]);
    {
        std::vector<std::size_t> cursor = node_corner_start;
        for (std::size_t ci = 0; ci < corners.size(); ++ci) {
            const UtgCorner& c = corners[ci];
            for (int m = 0; m < 4; ++m) {
                const std::int64_t nid = c.mNodes[static_cast<std::size_t>(m)];
                bool dup = false;
                for (int prev = 0; prev < m; ++prev)
                    if (c.mNodes[static_cast<std::size_t>(prev)] == nid)
                        dup = true;
                if (dup)
                    continue;
                std::uint8_t mask = 0;
                for (int k = m; k < 4; ++k)
                    if (c.mNodes[static_cast<std::size_t>(k)] == nid)
                        mask = static_cast<std::uint8_t>(mask | (1u << k));
                node_corners[cursor[static_cast<std::size_t>(nid)]++] = {
                    static_cast<std::uint32_t>(ci), mask};
            }
        }
    }

    if (rOptions.mFixBoundary) {
        std::sort(face_keys.begin(), face_keys.end());
        for (std::size_t i = 0; i < face_keys.size();) {
            std::size_t j = i;
            while (j < face_keys.size() && face_keys[j] == face_keys[i])
                ++j;
            if (j - i == 1)
                for (const std::int64_t v : face_keys[i])
                    if (v >= 0)
                        pinned[static_cast<std::size_t>(v)] = 1;
            i = j;
        }
    }
    for (std::size_t i = 0; i < n; ++i)
        if ((!rOptions.mFrozen.empty() && rOptions.mFrozen[i]) ||
            node_cell_start[i] == node_cell_start[i + 1])
            pinned[i] = 1;

    std::vector<Vec3> x(n, Vec3{0.0, 0.0, 0.0});
    {
        const detail::DoubleView pts_v(rMesh.Points());
        const double* ppts = pts_v.Data();
        const std::size_t pcols = dim < 3 ? dim : 3;
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t k = 0; k < pcols; ++k)
                x[i][k] = ppts[i * dim + k];
    }
    // The length unit: the mean edge of the corners.
    double edge_sum = 0.0;
    std::size_t edge_count = 0;
    for (const UtgCorner& c : corners)
        for (int k = 1; k < 4; ++k) {
            const Vec3& a = x[static_cast<std::size_t>(c.mNodes[0])];
            const Vec3& b = x[static_cast<std::size_t>(c.mNodes[static_cast<std::size_t>(k)])];
            edge_sum += std::sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) +
                                  (a[2] - b[2]) * (a[2] - b[2]));
            ++edge_count;
        }
    const double unit = edge_sum / static_cast<double>(std::max<std::size_t>(edge_count, 1));
    if (!(unit > 0.0))
        throw std::invalid_argument("untangle: every edge has zero length");
    const double inv_unit = 1.0 / unit;
    const double inv_unit3 = inv_unit * inv_unit * inv_unit;

    const auto cell_inverted = [&](std::size_t Cell) {
        for (auto c = cell_first[Cell]; c < cell_first[Cell + 1]; ++c)
            if (!utg_corner_positive(x, corners[static_cast<std::size_t>(c)]))
                return true;
        return false;
    };
    const auto count_and_quality = [&](std::int64_t& rInverted, double& rMinQ) {
        struct UtgCq {
            std::int64_t mInv = 0;
            double mMinQ = std::numeric_limits<double>::infinity();
        };
        const UtgCq res = parallel_reduce(
            num_cells, std::size_t{1024}, UtgCq{},
            [&](std::size_t lo, std::size_t hi) {
                UtgCq part;
                for (std::size_t cell = lo; cell < hi; ++cell) {
                    if (cell_inverted(cell))
                        ++part.mInv;
                    for (auto c = cell_first[cell]; c < cell_first[cell + 1]; ++c)
                        part.mMinQ =
                            std::min(part.mMinQ,
                                     utg_scaled_jacobian(x, corners[static_cast<std::size_t>(c)]));
                }
                return part;
            },
            [](const UtgCq& rA, const UtgCq& rB) {
                return UtgCq{rA.mInv + rB.mInv, std::min(rA.mMinQ, rB.mMinQ)};
            });
        rInverted = res.mInv;
        rMinQ = res.mMinQ;
    };
    UntangleResult result;
    count_and_quality(result.mNumInvertedBefore, result.mMinQualityBefore);

    // Every utg_term is non-negative, so backtracking against `Limit` (`f0`)
    // can exit as soon as the partial sum reaches or exceeds `Limit`.
    const auto node_objective_below = [&](std::size_t V, double Delta, double Limit) {
        double f = 0.0;
        for (std::size_t i = node_corner_start[V]; i < node_corner_start[V + 1]; ++i) {
            f += utg_term(utg_measure(x, corners[node_corners[i].mCorner], inv_unit), Delta);
            if (!(f < Limit))
                return false;
        }
        return f < Limit;
    };
    // Simultaneous objective and closed-form analytic gradient at vertex V.
    const auto node_obj_and_grad = [&](std::size_t V, double Delta, double& rF, Vec3& rG) {
        rF = 0.0;
        rG = {0.0, 0.0, 0.0};
        for (std::size_t i = node_corner_start[V]; i < node_corner_start[V + 1]; ++i) {
            const UtgNodeCorner nc = node_corners[i];
            if (!utg_term_and_grad(x, corners[nc.mCorner], inv_unit, Delta, nc.mMask, rF, rG))
                return false;
        }
        return std::isfinite(rF) && std::isfinite(rG[0]) && std::isfinite(rG[1]) &&
               std::isfinite(rG[2]);
    };
    // Up to 4 steepest-descent + backtracking steps for one vertex, warm-
    // starting the trial step length from the previous accepted step.
    const auto move_vertex = [&](std::size_t V, double Delta) -> std::uint8_t {
        std::uint8_t moves = 0;
        double max_step = 0.5 * unit;
        for (int step = 0; step < 4; ++step) {
            double f0 = 0.0;
            Vec3 g{0.0, 0.0, 0.0};
            if (!node_obj_and_grad(V, Delta, f0, g))
                break;
            const double gn = std::sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]);
            if (!(gn > 0.0))
                break;
            const Vec3 keep = x[V];
            double step_len = max_step;
            bool moved = false;
            for (int tries = 0; tries < 30; ++tries, step_len *= 0.5) {
                const double alpha = step_len / gn;
                for (std::size_t k = 0; k < 3; ++k)
                    x[V][k] = keep[k] - alpha * g[k];
                if (node_objective_below(V, Delta, f0)) {
                    moved = true;
                    max_step = std::min(0.5 * unit, 2.0 * step_len);
                    break;
                }
            }
            if (!moved) {
                x[V] = keep;
                break;
            }
            ++moves;
        }
        return moves;
    };

    // Deterministic greedy graph colouring of unpinned vertices (built lazily
    // on the first active sweep). Two vertices conflict iff they share a corner,
    // so every vertex within one colour batch updates `x[v]` independently.
    bool colored = false;
    std::vector<std::size_t> color_start;
    std::vector<std::size_t> color_nodes;
    const auto ensure_colored = [&]() {
        if (colored)
            return;
        colored = true;
        std::vector<std::int32_t> color(n, -1);
        std::vector<std::size_t> stamp;
        std::size_t num_colors = 0;
        for (std::size_t v = 0; v < n; ++v) {
            if (pinned[v])
                continue;
            const std::size_t mark = v + 1;
            for (std::size_t i = node_corner_start[v]; i < node_corner_start[v + 1]; ++i) {
                const UtgCorner& c = corners[node_corners[i].mCorner];
                for (const std::int64_t nid : c.mNodes) {
                    const auto u = static_cast<std::size_t>(nid);
                    if (u < v && color[u] >= 0) {
                        const auto col = static_cast<std::size_t>(color[u]);
                        if (col >= stamp.size())
                            stamp.resize(col + 1, 0);
                        stamp[col] = mark;
                    }
                }
            }
            std::size_t chosen = 0;
            while (chosen < stamp.size() && stamp[chosen] == mark)
                ++chosen;
            color[v] = static_cast<std::int32_t>(chosen);
            if (chosen + 1 > num_colors)
                num_colors = chosen + 1;
        }
        color_start.assign(num_colors + 1, 0);
        for (std::size_t v = 0; v < n; ++v)
            if (color[v] >= 0)
                ++color_start[static_cast<std::size_t>(color[v]) + 1];
        for (std::size_t c = 0; c < num_colors; ++c)
            color_start[c + 1] += color_start[c];
        color_nodes.resize(color_start[num_colors]);
        std::vector<std::size_t> cursor = color_start;
        for (std::size_t v = 0; v < n; ++v)
            if (color[v] >= 0)
                color_nodes[cursor[static_cast<std::size_t>(color[v])]++] = v;
    };

    double previous_total = std::numeric_limits<double>::infinity();
    std::vector<std::size_t> batch;
    std::vector<std::uint8_t> batch_moves;
    for (std::int64_t sweep = 1; sweep <= rOptions.mMaxSweeps; ++sweep) {
        // delta from the most negative corner (Escobar et al.'s rule).
        const double min_det = parallel_reduce(
            corners.size(), std::size_t{4096}, std::numeric_limits<double>::infinity(),
            [&](std::size_t lo, std::size_t hi) {
                double best = std::numeric_limits<double>::infinity();
                for (std::size_t i = lo; i < hi; ++i)
                    best = std::min(best, utg_det(x, corners[i], inv_unit3));
                return best;
            },
            [](double a, double b) { return std::min(a, b); });
        const double eps = 1e-3;
        const double delta = min_det < eps ? std::sqrt(eps * (eps - min_det)) : 0.0;
        // Which nodes move this sweep.
        std::vector<char> active(n, 0);
        if (rOptions.mSmooth) {
            std::fill(active.begin(), active.end(), 1);
        } else {
            for (std::size_t cell = 0; cell < num_cells; ++cell)
                if (cell_inverted(cell))
                    for (auto c = cell_first[cell]; c < cell_first[cell + 1]; ++c)
                        for (const std::int64_t v : corners[static_cast<std::size_t>(c)].mNodes)
                            active[static_cast<std::size_t>(v)] = 1;
            for (std::int64_t ring = 0; ring < rOptions.mRings; ++ring) {
                std::vector<char> grown = active;
                for (std::size_t v = 0; v < n; ++v)
                    if (active[v])
                        for (std::size_t ci = node_cell_start[v]; ci < node_cell_start[v + 1];
                             ++ci) {
                            const auto cell = static_cast<std::size_t>(node_cell_ids[ci]);
                            for (auto c = cell_first[cell]; c < cell_first[cell + 1]; ++c)
                                for (const std::int64_t q :
                                     corners[static_cast<std::size_t>(c)].mNodes)
                                    grown[static_cast<std::size_t>(q)] = 1;
                        }
                active.swap(grown);
            }
        }
        bool any = false;
        for (std::size_t v = 0; v < n; ++v)
            if (active[v] && !pinned[v]) {
                any = true;
                break;
            }
        if (!any)
            break;
        result.mSweeps = sweep;
        ensure_colored();
        const std::size_t num_colors = color_start.empty() ? 0 : color_start.size() - 1;
        for (std::size_t col = 0; col < num_colors; ++col) {
            batch.clear();
            for (std::size_t idx = color_start[col]; idx < color_start[col + 1]; ++idx) {
                const std::size_t v = color_nodes[idx];
                if (active[v])
                    batch.push_back(v);
            }
            if (batch.empty())
                continue;
            batch_moves.assign(batch.size(), 0);
            parallel_for(
                batch.size(),
                [&](std::size_t bi) { batch_moves[bi] = move_vertex(batch[bi], delta); },
                std::size_t{16});
            for (std::uint8_t mv : batch_moves)
                result.mMoves += mv;
        }
        bool any_inverted = false;
        for (std::size_t cell = 0; cell < num_cells; ++cell)
            if (cell_inverted(cell)) {
                any_inverted = true;
                break;
            }
        if (!any_inverted) {
            if (!rOptions.mSmooth)
                break;
            const double total = parallel_reduce(
                corners.size(), std::size_t{4096}, 0.0,
                [&](std::size_t lo, std::size_t hi) {
                    double s = 0.0;
                    for (std::size_t i = lo; i < hi; ++i)
                        s += utg_term(utg_measure(x, corners[i], inv_unit), 0.0);
                    return s;
                },
                [](double a, double b) { return a + b; });
            if (std::isfinite(previous_total) &&
                previous_total - total <= 1e-6 * std::abs(previous_total))
                break;
            previous_total = total;
        }
    }
    count_and_quality(result.mNumInvertedAfter, result.mMinQualityAfter);

    NDArray pts(DType::Float64, {n, dim});
    double* pp = pts.As<double>();
    const detail::DoubleView orig_pts(rMesh.Points());
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t k = 0; k < dim; ++k)
            pp[i * dim + k] = k < 3 ? x[i][k] : orig_pts[i * dim + k];
    result.mMesh = detail::clone_mesh(rMesh);
    result.mMesh.AssignPoints(std::move(pts));
    return result;
}

}  // namespace meshioplusplus

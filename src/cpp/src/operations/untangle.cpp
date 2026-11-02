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
#include <cstdint>
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
#include "meshioplusplus/detail/value_io.hpp"

namespace meshioplusplus {

namespace {

using UtgVec = std::array<double, 3>;

/// A corner simplex: the corner node, its three edge ends (right-handed for a
/// valid cell), and whether it is a tetrahedron (ideal: regular) or a
/// hexahedron corner (ideal: a cube corner).
struct UtgCorner {
    std::array<std::int64_t, 4> mNodes;
    bool mTet;
};

// Knupp's hexahedron corners: corner c, then its three edge neighbours.
constexpr int kUtgHexCorner[8][4] = {{0, 1, 3, 4}, {1, 2, 0, 5}, {2, 3, 1, 6}, {3, 0, 2, 7},
                                     {4, 7, 5, 0}, {5, 4, 6, 1}, {6, 5, 7, 2}, {7, 6, 4, 3}};

// The inverse of the regular tetrahedron's edge matrix (columns: its three
// edges from corner 0, unit length), row-major.
const double kUtgTetWinv[3][3] = {{1.0, -0.57735026918962573, -0.40824829046386302},
                                  {0.0, 1.1547005383792515, -0.40824829046386302},
                                  {0.0, 0.0, 1.2247448713915890}};

struct UtgMeasure {
    double mDet;    ///< det S
    double mFrob2;  ///< |S|_F^2
};

UtgMeasure utg_measure(const std::vector<UtgVec>& rX, const UtgCorner& rC, double InvScale) {
    const UtgVec& p0 = rX[static_cast<std::size_t>(rC.mNodes[0])];
    double a[3][3];  // columns: the three edges, scaled
    for (int k = 0; k < 3; ++k) {
        const UtgVec& p = rX[static_cast<std::size_t>(rC.mNodes[k + 1])];
        for (int r = 0; r < 3; ++r)
            a[r][k] = (p[r] - p0[r]) * InvScale;
    }
    double s[3][3];
    if (rC.mTet) {
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c) {
                double v = 0.0;
                for (int k = 0; k < 3; ++k)
                    v += a[r][k] * kUtgTetWinv[k][c];
                s[r][c] = v;
            }
    } else {
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c)
                s[r][c] = a[r][c];
    }
    UtgMeasure m{};
    m.mDet = s[0][0] * (s[1][1] * s[2][2] - s[1][2] * s[2][1]) -
             s[0][1] * (s[1][0] * s[2][2] - s[1][2] * s[2][0]) +
             s[0][2] * (s[1][0] * s[2][1] - s[1][1] * s[2][0]);
    m.mFrob2 = 0.0;
    for (int r = 0; r < 3; ++r)
        for (int c = 0; c < 3; ++c)
            m.mFrob2 += s[r][c] * s[r][c];
    return m;
}

/// The corner's scaled Jacobian: det of the edges over their lengths' product.
double utg_scaled_jacobian(const std::vector<UtgVec>& rX, const UtgCorner& rC) {
    const UtgVec& p0 = rX[static_cast<std::size_t>(rC.mNodes[0])];
    UtgVec e[3];
    double len = 1.0;
    for (int k = 0; k < 3; ++k) {
        const UtgVec& p = rX[static_cast<std::size_t>(rC.mNodes[k + 1])];
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

double utg_term(const UtgMeasure& rM, double Delta) {
    const double h = 0.5 * (rM.mDet + std::sqrt(rM.mDet * rM.mDet + 4.0 * Delta * Delta));
    if (!(h > 0.0))
        return std::numeric_limits<double>::infinity();
    return rM.mFrob2 / (3.0 * std::cbrt(h * h));
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
    std::vector<std::vector<std::int64_t>> node_cells(n);
    std::vector<char> pinned(n, 0);
    std::vector<std::array<std::int64_t, 4>> face_keys;
    for (const auto cb : rMesh.CellRange()) {
        const std::string& type = cb.Type();
        const bool tet = type == "tetra", hex = type == "hexahedron";
        const int d = cb.IsRagged() ? (cb.IsPolyhedron() ? 3 : 2)
                                    : cell_type_dimension(cell_type_from_name(type));
        if (d < 3)
            continue;
        const std::size_t npc = cb.IsRagged() ? 0 : cb.NodesPerCell();
        for (std::size_t c = 0; c < cb.NumCells(); ++c) {
            if (!tet && !hex) {
                // Not measured: its nodes stay where they are.
                if (cb.IsRagged()) {
                    for (std::size_t f = 0; f < cb.NumFaces(c); ++f) {
                        const auto face = cb.Face(c, f);
                        for (std::size_t k = 0; k < face.second; ++k)
                            pinned[static_cast<std::size_t>(face.first[k])] = 1;
                    }
                } else {
                    for (std::size_t k = 0; k < npc; ++k)
                        pinned[static_cast<std::size_t>(detail::read_int(cb.Conn(), c * npc + k))] =
                            1;
                }
                continue;
            }
            std::array<std::int64_t, 8> v{};
            for (std::size_t k = 0; k < npc; ++k)
                v[k] = detail::read_int(cb.Conn(), c * npc + k);
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
            for (std::size_t k = 0; k < npc; ++k)
                node_cells[static_cast<std::size_t>(v[k])].push_back(id);
            for (const auto& f : detail::cell_faces(cell_type_from_name(type))) {
                std::array<std::int64_t, 4> key{-1, -1, -1, -1};
                for (std::size_t k = 0; k < f.mNumCorners; ++k)
                    key[k] = v[f.mNodes[k]];
                std::sort(key.begin(), key.begin() + f.mNumCorners);
                face_keys.push_back(key);
            }
        }
    }
    cell_first.push_back(static_cast<std::int64_t>(corners.size()));
    const std::size_t num_cells = cell_first.size() - 1;
    if (num_cells == 0)
        throw std::invalid_argument("untangle: the mesh has no tetrahedron or hexahedron");
    for (auto& cells : node_cells) {
        std::sort(cells.begin(), cells.end());
        cells.erase(std::unique(cells.begin(), cells.end()), cells.end());
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
        if ((!rOptions.mFrozen.empty() && rOptions.mFrozen[i]) || node_cells[i].empty())
            pinned[i] = 1;

    std::vector<UtgVec> x(n, UtgVec{0.0, 0.0, 0.0});
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t k = 0; k < dim && k < 3; ++k)
            x[i][k] = detail::read_double(rMesh.Points(), i * dim + k);
    // The length unit: the mean edge of the corners.
    double edge_sum = 0.0;
    std::size_t edge_count = 0;
    for (const UtgCorner& c : corners)
        for (int k = 1; k < 4; ++k) {
            const UtgVec& a = x[static_cast<std::size_t>(c.mNodes[0])];
            const UtgVec& b = x[static_cast<std::size_t>(c.mNodes[k])];
            edge_sum += std::sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) +
                                  (a[2] - b[2]) * (a[2] - b[2]));
            ++edge_count;
        }
    const double unit = edge_sum / static_cast<double>(std::max<std::size_t>(edge_count, 1));
    if (!(unit > 0.0))
        throw std::invalid_argument("untangle: every edge has zero length");
    const double inv_unit = 1.0 / unit;

    const auto cell_inverted = [&](std::size_t Cell) {
        for (auto c = cell_first[Cell]; c < cell_first[Cell + 1]; ++c)
            if (!(utg_measure(x, corners[static_cast<std::size_t>(c)], inv_unit).mDet > 0.0))
                return true;
        return false;
    };
    const auto count_and_quality = [&](std::int64_t& rInverted, double& rMinQ) {
        rInverted = 0;
        rMinQ = std::numeric_limits<double>::infinity();
        for (std::size_t cell = 0; cell < num_cells; ++cell) {
            if (cell_inverted(cell))
                ++rInverted;
            for (auto c = cell_first[cell]; c < cell_first[cell + 1]; ++c)
                rMinQ =
                    std::min(rMinQ, utg_scaled_jacobian(x, corners[static_cast<std::size_t>(c)]));
        }
    };
    UntangleResult result;
    count_and_quality(result.mNumInvertedBefore, result.mMinQualityBefore);

    // Each node's objective: the terms of every corner of its cells.
    const auto node_objective = [&](std::size_t V, double Delta) {
        double f = 0.0;
        for (const std::int64_t cell : node_cells[V])
            for (auto c = cell_first[static_cast<std::size_t>(cell)];
                 c < cell_first[static_cast<std::size_t>(cell) + 1]; ++c)
                f +=
                    utg_term(utg_measure(x, corners[static_cast<std::size_t>(c)], inv_unit), Delta);
        return f;
    };
    // The nodes of each cell, for the rings.
    const auto cell_nodes = [&](std::size_t Cell, std::vector<std::int64_t>& rOut) {
        rOut.clear();
        for (auto c = cell_first[Cell]; c < cell_first[Cell + 1]; ++c)
            for (const std::int64_t v : corners[static_cast<std::size_t>(c)].mNodes)
                rOut.push_back(v);
    };

    double previous_total = std::numeric_limits<double>::infinity();
    std::vector<std::int64_t> nodes;
    for (std::int64_t sweep = 1; sweep <= rOptions.mMaxSweeps; ++sweep) {
        // delta from the most negative corner (Escobar et al.'s rule).
        double min_det = std::numeric_limits<double>::infinity();
        for (const UtgCorner& c : corners)
            min_det = std::min(min_det, utg_measure(x, c, inv_unit).mDet);
        const double eps = 1e-3;
        const double delta = min_det < eps ? std::sqrt(eps * (eps - min_det)) : 0.0;
        // Which nodes move this sweep.
        std::vector<char> active(n, 0);
        if (rOptions.mSmooth) {
            std::fill(active.begin(), active.end(), 1);
        } else {
            for (std::size_t cell = 0; cell < num_cells; ++cell)
                if (cell_inverted(cell)) {
                    cell_nodes(cell, nodes);
                    for (const std::int64_t v : nodes)
                        active[static_cast<std::size_t>(v)] = 1;
                }
            for (std::int64_t ring = 0; ring < rOptions.mRings; ++ring) {
                std::vector<char> grown = active;
                for (std::size_t v = 0; v < n; ++v)
                    if (active[v])
                        for (const std::int64_t cell : node_cells[v]) {
                            cell_nodes(static_cast<std::size_t>(cell), nodes);
                            for (const std::int64_t q : nodes)
                                grown[static_cast<std::size_t>(q)] = 1;
                        }
                active.swap(grown);
            }
        }
        bool any = false;
        for (std::size_t v = 0; v < n; ++v)
            any = any || (active[v] && !pinned[v]);
        if (!any)
            break;
        result.mSweeps = sweep;
        for (std::size_t v = 0; v < n; ++v) {
            if (!active[v] || pinned[v])
                continue;
            const double h = 1e-6 * unit;
            for (int step = 0; step < 4; ++step) {
                const double f0 = node_objective(v, delta);
                UtgVec g{0.0, 0.0, 0.0};
                bool finite = std::isfinite(f0);
                for (std::size_t k = 0; k < 3 && finite; ++k) {
                    const double keep = x[v][k];
                    x[v][k] = keep + h;
                    const double fp = node_objective(v, delta);
                    x[v][k] = keep - h;
                    const double fm = node_objective(v, delta);
                    x[v][k] = keep;
                    g[k] = (fp - fm) / (2.0 * h);
                    finite = std::isfinite(g[k]);
                }
                const double gn = std::sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]);
                if (!finite || !(gn > 0.0))
                    break;
                // Backtracking from a step of half the length unit.
                const UtgVec keep = x[v];
                double alpha = 0.5 * unit / gn;
                bool moved = false;
                for (int tries = 0; tries < 30; ++tries, alpha *= 0.5) {
                    for (std::size_t k = 0; k < 3; ++k)
                        x[v][k] = keep[k] - alpha * g[k];
                    if (node_objective(v, delta) < f0) {
                        moved = true;
                        break;
                    }
                }
                if (!moved) {
                    x[v] = keep;
                    break;
                }
                ++result.mMoves;
            }
        }
        std::int64_t inverted = 0;
        for (std::size_t cell = 0; cell < num_cells && inverted == 0; ++cell)
            if (cell_inverted(cell))
                ++inverted;
        if (inverted == 0) {
            if (!rOptions.mSmooth)
                break;
            double total = 0.0;
            for (const UtgCorner& c : corners)
                total += utg_term(utg_measure(x, c, inv_unit), 0.0);
            if (std::isfinite(previous_total) &&
                previous_total - total <= 1e-6 * std::abs(previous_total))
                break;
            previous_total = total;
        }
    }
    count_and_quality(result.mNumInvertedAfter, result.mMinQualityAfter);

    NDArray pts(DType::Float64, {n, dim});
    double* pp = pts.As<double>();
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t k = 0; k < dim; ++k)
            pp[i * dim + k] = k < 3 ? x[i][k] : detail::read_double(rMesh.Points(), i * dim + k);
    result.mMesh = detail::clone_mesh(rMesh);
    result.mMesh.AssignPoints(std::move(pts));
    return result;
}

}  // namespace meshioplusplus

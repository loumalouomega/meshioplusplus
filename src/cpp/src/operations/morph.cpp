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
#include "meshioplusplus/operations/morph.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/cell_faces.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/polyhedron.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/operations/convert_cells.hpp"
#include "meshioplusplus/detail/geometry.hpp"

// Project includes (private, not installed)
#include "../detail/sparse_solve.hpp"

namespace meshioplusplus {

namespace {

using MphVec = std::array<double, 3>;

/// The cell types `morph` moves: linear volume and surface cells, and the
/// lower-dimensional ones it carries along.
bool mph_supported(const std::string& rType) {
    static const char* const kTypes[] = {"tetra", "hexahedron", "wedge", "pyramid", "triangle",
                                         "quad",  "polygon",    "line",  "vertex"};
    return std::any_of(std::begin(kTypes), std::end(kTypes),
                       [&](const char* pT) { return rType == pT; });
}

int mph_dimension(const Mesh::CellView& rCb) {
    if (rCb.IsRagged())
        return 2;
    return cell_type_dimension(cell_type_from_name(rCb.Type()));
}

/// The corner ids of row `C` of a 2-D block, in ring order.
void mph_ring(const Mesh::CellView& rCb, std::size_t C, std::vector<std::int64_t>& rOut) {
    rOut.clear();
    if (rCb.IsRagged()) {
        const std::int64_t* row = rCb.Row(C);
        rOut.assign(row, row + rCb.RowSize(C));
        return;
    }
    const std::size_t npc = rCb.NodesPerCell();
    for (std::size_t k = 0; k < npc; ++k)
        rOut.push_back(detail::read_int(rCb.Conn(), C * npc + k));
}

/// The points of facets used by exactly one top-dimensional cell.
std::vector<char> mph_boundary(const Mesh& rMesh, int D, std::size_t N) {
    std::vector<std::array<std::int64_t, 4>> keys;
    std::vector<std::int64_t> ring;
    for (const auto cb : rMesh.CellRange()) {
        if (mph_dimension(cb) != D)
            continue;
        if (D == 3) {
            const CellType type = cell_type_from_name(cb.Type());
            const auto& faces = detail::cell_faces(type);
            const std::size_t npc = cb.NodesPerCell();
            for (std::size_t c = 0; c < cb.NumCells(); ++c)
                for (const auto& f : faces) {
                    std::array<std::int64_t, 4> key{-1, -1, -1, -1};
                    for (std::size_t k = 0; k < f.mNumCorners; ++k)
                        key[k] = detail::read_int(cb.Conn(), c * npc + f.mNodes[k]);
                    std::sort(key.begin(), key.begin() + f.mNumCorners);
                    keys.push_back(key);
                }
        } else {
            for (std::size_t c = 0; c < cb.NumCells(); ++c) {
                mph_ring(cb, c, ring);
                for (std::size_t k = 0; k < ring.size(); ++k) {
                    std::int64_t a = ring[k], b = ring[(k + 1) % ring.size()];
                    if (a > b)
                        std::swap(a, b);
                    keys.push_back({a, b, -1, -1});
                }
            }
        }
    }
    std::sort(keys.begin(), keys.end());
    std::vector<char> out(N, 0);
    for (std::size_t i = 0; i < keys.size();) {
        std::size_t j = i;
        while (j < keys.size() && keys[j] == keys[i])
            ++j;
        if (j - i == 1)
            for (const std::int64_t v : keys[i])
                if (v >= 0 && static_cast<std::size_t>(v) < N)
                    out[static_cast<std::size_t>(v)] = 1;
        i = j;
    }
    return out;
}

/// Per cell, its orientation measure: a volume cell's signed volume, a
/// surface cell's Newell normal (NaN when not measured).
struct MphOrientation {
    std::vector<double> mVolume;
    std::vector<MphVec> mNormal;
};

MphOrientation mph_orientation(const Mesh& rMesh, const std::vector<MphVec>& rXyz) {
    MphOrientation out;
    std::vector<std::int64_t> ring;
    std::vector<detail::Vec3> corners;
    for (const auto cb : rMesh.CellRange()) {
        const int d = mph_dimension(cb);
        for (std::size_t c = 0; c < cb.NumCells(); ++c) {
            double vol = std::numeric_limits<double>::quiet_NaN();
            MphVec nrm{0.0, 0.0, 0.0};
            if (d == 3) {
                const std::size_t npc = cb.NodesPerCell();
                corners.resize(npc);
                for (std::size_t k = 0; k < npc; ++k) {
                    const MphVec& p =
                        rXyz[static_cast<std::size_t>(detail::read_int(cb.Conn(), c * npc + k))];
                    corners[k] = {p[0], p[1], p[2]};
                }
                vol = detail::cell_volume_from_corners(corners.data(),
                                                       cell_type_from_name(cb.Type()));
            } else if (d == 2) {
                mph_ring(cb, c, ring);
                for (std::size_t k = 0; k < ring.size(); ++k) {
                    const MphVec& a = rXyz[static_cast<std::size_t>(ring[k])];
                    const MphVec& b = rXyz[static_cast<std::size_t>(ring[(k + 1) % ring.size()])];
                    nrm[0] += (a[1] - b[1]) * (a[2] + b[2]);
                    nrm[1] += (a[2] - b[2]) * (a[0] + b[0]);
                    nrm[2] += (a[0] - b[0]) * (a[1] + b[1]);
                }
            }
            out.mVolume.push_back(vol);
            out.mNormal.push_back(nrm);
        }
    }
    return out;
}

/// P1 stiffness of a simplex (D + 1 corners), times `Weight`, into triplets;
/// returns its measure (0 when degenerate, nothing added).
double mph_simplex_measure(const std::vector<MphVec>& rXyz, const std::int64_t* pNodes, int D) {
    const MphVec& p0 = rXyz[static_cast<std::size_t>(pNodes[0])];
    MphVec e[3];
    for (int k = 0; k < D; ++k)
        e[k] = detail::vec3_sub(rXyz[static_cast<std::size_t>(pNodes[k + 1])], p0);
    if (D == 2) {
        const double g00 = detail::vec3_dot(e[0], e[0]), g01 = detail::vec3_dot(e[0], e[1]),
                     g11 = detail::vec3_dot(e[1], e[1]);
        const double det = g00 * g11 - g01 * g01;
        return det > 0.0 ? 0.5 * std::sqrt(det) : 0.0;
    }
    const MphVec x{e[1][1] * e[2][2] - e[1][2] * e[2][1], e[1][2] * e[2][0] - e[1][0] * e[2][2],
                   e[1][0] * e[2][1] - e[1][1] * e[2][0]};
    return std::abs(detail::vec3_dot(e[0], x)) / 6.0;
}

void mph_simplex_stiffness(const std::vector<MphVec>& rXyz, const std::int64_t* pNodes, int D,
                           double Measure, double Weight,
                           std::vector<detail::SparseTriplet>& rOut) {
    if (D != 2 && D != 3)
        return;  // triangles and tetrahedra only, which is all the caller passes
    const MphVec& p0 = rXyz[static_cast<std::size_t>(pNodes[0])];
    MphVec e[3];
    for (int k = 0; k < D; ++k)
        e[k] = detail::vec3_sub(rXyz[static_cast<std::size_t>(pNodes[k + 1])], p0);
    // G = E^T E, the gradients' Gram matrix is its inverse: K = |c| B G^-1 B^T.
    double g[3][3], inv[3][3];
    for (int a = 0; a < D; ++a)
        for (int b = 0; b < D; ++b)
            g[a][b] = detail::vec3_dot(e[a], e[b]);
    if (D == 2) {
        const double det = g[0][0] * g[1][1] - g[0][1] * g[1][0];
        inv[0][0] = g[1][1] / det;
        inv[1][1] = g[0][0] / det;
        inv[0][1] = inv[1][0] = -g[0][1] / det;
    } else {
        const double det = g[0][0] * (g[1][1] * g[2][2] - g[1][2] * g[2][1]) -
                           g[0][1] * (g[1][0] * g[2][2] - g[1][2] * g[2][0]) +
                           g[0][2] * (g[1][0] * g[2][1] - g[1][1] * g[2][0]);
        inv[0][0] = (g[1][1] * g[2][2] - g[1][2] * g[2][1]) / det;
        inv[0][1] = (g[0][2] * g[2][1] - g[0][1] * g[2][2]) / det;
        inv[0][2] = (g[0][1] * g[1][2] - g[0][2] * g[1][1]) / det;
        inv[1][0] = (g[1][2] * g[2][0] - g[1][0] * g[2][2]) / det;
        inv[1][1] = (g[0][0] * g[2][2] - g[0][2] * g[2][0]) / det;
        inv[1][2] = (g[0][2] * g[1][0] - g[0][0] * g[1][2]) / det;
        inv[2][0] = (g[1][0] * g[2][1] - g[1][1] * g[2][0]) / det;
        inv[2][1] = (g[0][1] * g[2][0] - g[0][0] * g[2][1]) / det;
        inv[2][2] = (g[0][0] * g[1][1] - g[0][1] * g[1][0]) / det;
    }
    // Local index 0 is the base corner: its row is minus the sum of the others.
    double k[4][4];
    for (int a = 0; a < D; ++a)
        for (int b = 0; b < D; ++b)
            k[a + 1][b + 1] = inv[a][b];
    for (int a = 1; a <= D; ++a) {
        double s = 0.0;
        for (int b = 1; b <= D; ++b)
            s += k[a][b];
        k[a][0] = k[0][a] = -s;
    }
    double s00 = 0.0;
    for (int a = 1; a <= D; ++a)
        s00 -= k[a][0];
    k[0][0] = s00;
    const double scale = Measure * Weight;
    for (int a = 0; a <= D; ++a)
        for (int b = 0; b <= D; ++b)
            rOut.push_back({static_cast<std::size_t>(pNodes[a]),
                            static_cast<std::size_t>(pNodes[b]), scale * k[a][b]});
}

}  // namespace

MorphResult morph(const Mesh& rMesh, const MorphOptions& rOptions) {
    if (rOptions.mArray.empty())
        throw std::invalid_argument("morph: name the displacement point data ('array')");
    if (!rMesh.HasPointData(rOptions.mArray))
        throw std::invalid_argument("morph: no point data named '" + rOptions.mArray + "'");
    if (!(rOptions.mStiffening >= 0.0) || !std::isfinite(rOptions.mStiffening))
        throw std::invalid_argument("morph: stiffening must not be negative");
    if (rOptions.mSteps < 1)
        throw std::invalid_argument("morph: steps must be at least 1");
    const std::size_t n = rMesh.NumPoints();
    const std::size_t dim = rMesh.PointDim();
    const NDArray& disp = rMesh.PointData(rOptions.mArray);
    const std::size_t comps = n > 0 ? disp.Size() / n : 0;
    if (n == 0 || disp.Size() != n * comps || (comps != dim && comps != 3))
        throw std::invalid_argument("morph: '" + rOptions.mArray +
                                    "' must be (n, dim) or (n, 3), one row per point");
    if (!rOptions.mConstrained.empty() && rOptions.mConstrained.size() != n)
        throw std::invalid_argument("morph: the constrained mask must have one entry per point");

    int top = 0;
    for (const auto cb : rMesh.CellRange()) {
        if (cb.IsPolyhedron())
            throw std::invalid_argument(
                "morph: polyhedron cells have no simplex template; decompose them with "
                "poly_to_tet first");
        if (!cb.IsRagged() && !mph_supported(cb.Type()))
            throw std::invalid_argument("morph: '" + cb.Type() +
                                        "' cells are not linear; linearize the mesh first");
        top = std::max(top, mph_dimension(cb));
    }
    if (top < 2)
        throw std::invalid_argument("morph: the mesh has no surface or volume cell to morph");

    // The simplices of the top dimension (points unchanged by simplexify).
    std::vector<std::int64_t> simplices;
    {
        ConvertCellsOptions co;
        co.mMode = ConvertCellsMode::Simplexify;
        const ConvertCellsResult split = convert_cells(detail::clone_geometry(rMesh), co);
        const std::string want = top == 3 ? "tetra" : "triangle";
        for (const auto cb : split.mMesh.CellRange()) {
            if (cb.Type() != want || cb.IsRagged())
                continue;
            const std::size_t npc = cb.NodesPerCell();
            for (std::size_t i = 0; i < cb.NumCells() * npc; ++i)
                simplices.push_back(detail::read_int(cb.Conn(), i));
        }
    }
    const std::size_t corners = static_cast<std::size_t>(top) + 1;
    const std::size_t num_simplices = simplices.size() / corners;

    std::vector<MphVec> xyz(n, MphVec{0.0, 0.0, 0.0});
    std::vector<MphVec> target(n, MphVec{0.0, 0.0, 0.0});
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t k = 0; k < dim && k < 3; ++k)
            xyz[i][k] = detail::read_double(rMesh.Points(), i * dim + k);
        for (std::size_t k = 0; k < comps && k < 3; ++k) {
            const double v = detail::read_double(disp, i * comps + k);
            target[i][k] = k < dim ? v : 0.0;
        }
    }
    const std::vector<MphVec> start = xyz;

    // Constrained: boundary, the masks, and every point no simplex uses.
    std::vector<char> fixed(n, 0);
    if (rOptions.mFixBoundary)
        fixed = mph_boundary(rMesh, top, n);
    for (std::size_t i = 0; i < n; ++i)
        if (!rOptions.mConstrained.empty() && rOptions.mConstrained[i])
            fixed[i] = 1;
    if (!rOptions.mConstrainedArray.empty()) {
        if (!rMesh.HasPointData(rOptions.mConstrainedArray))
            throw std::invalid_argument("morph: no point data named '" +
                                        rOptions.mConstrainedArray + "'");
        const NDArray& mask = rMesh.PointData(rOptions.mConstrainedArray);
        if (mask.Size() != n)
            throw std::invalid_argument("morph: '" + rOptions.mConstrainedArray +
                                        "' must have one entry per point");
        for (std::size_t i = 0; i < n; ++i)
            if (detail::read_double(mask, i) != 0.0)
                fixed[i] = 1;
    }
    std::vector<char> used(n, 0);
    for (const std::int64_t v : simplices)
        used[static_cast<std::size_t>(v)] = 1;
    std::vector<std::int64_t> slot(n, -1);
    std::size_t num_free = 0;
    MorphResult result;
    for (std::size_t i = 0; i < n; ++i) {
        if (!used[i])
            fixed[i] = 1;
        if (fixed[i])
            ++result.mNumConstrained;
        else
            slot[i] = static_cast<std::int64_t>(num_free++);
    }
    result.mNumFree = static_cast<std::int64_t>(num_free);
    const std::size_t cap = rOptions.mMaxIterations > 0
                                ? static_cast<std::size_t>(rOptions.mMaxIterations)
                                : std::max<std::size_t>(2 * num_free + 10, 1000);

    const double frac = 1.0 / static_cast<double>(rOptions.mSteps);
    std::vector<double> measure(num_simplices);
    for (std::int64_t step = 0; step < rOptions.mSteps; ++step) {
        // The stiffness of the current geometry.
        double mean = 0.0;
        std::size_t live = 0;
        for (std::size_t s = 0; s < num_simplices; ++s) {
            measure[s] = mph_simplex_measure(xyz, &simplices[s * corners], top);
            if (measure[s] > 0.0) {
                mean += measure[s];
                ++live;
            }
        }
        if (live == 0)
            throw std::invalid_argument("morph: every cell has zero measure");
        mean /= static_cast<double>(live);
        std::vector<detail::SparseTriplet> local;
        std::vector<detail::SparseTriplet> kff;
        std::vector<std::array<double, 3>> rhs(num_free, std::array<double, 3>{0.0, 0.0, 0.0});
        for (std::size_t s = 0; s < num_simplices; ++s) {
            // A cell flattened to rounding (relative to the mean) adds nothing.
            if (!(measure[s] > 1e-14 * mean))
                continue;
            const double weight = rOptions.mStiffening > 0.0
                                      ? std::pow(mean / measure[s], rOptions.mStiffening)
                                      : 1.0;
            local.clear();
            mph_simplex_stiffness(xyz, &simplices[s * corners], top, measure[s], weight, local);
            for (const detail::SparseTriplet& t : local) {
                if (fixed[t.mRow])
                    continue;
                const auto r = static_cast<std::size_t>(slot[t.mRow]);
                if (fixed[t.mCol]) {
                    for (std::size_t a = 0; a < 3; ++a)
                        rhs[r][a] -= t.mValue * frac * target[t.mCol][a];
                } else {
                    kff.push_back({r, static_cast<std::size_t>(slot[t.mCol]), t.mValue});
                }
            }
        }
        if (num_free > 0) {
            const detail::SparseCsr a =
                detail::csr_from_triplets(num_free, num_free, std::move(kff));
            std::vector<double> b(num_free), u(num_free);
            for (std::size_t c = 0; c < dim && c < 3; ++c) {
                for (std::size_t r = 0; r < num_free; ++r)
                    b[r] = rhs[r][c];
                std::fill(u.begin(), u.end(), 0.0);
                const detail::SparseSolveReport rep =
                    detail::csr_pcg(a, b.data(), u.data(), rOptions.mTolerance, cap);
                result.mIterations += static_cast<std::int64_t>(rep.mIterations);
                result.mConverged = result.mConverged && rep.mConverged;
                for (std::size_t i = 0; i < n; ++i)
                    if (!fixed[i])
                        xyz[i][c] += u[static_cast<std::size_t>(slot[i])];
            }
        }
        for (std::size_t i = 0; i < n; ++i)
            if (fixed[i])
                for (std::size_t c = 0; c < dim && c < 3; ++c)
                    xyz[i][c] += frac * target[i][c];
    }
    // A constrained point lands exactly on its target, whatever the steps.
    for (std::size_t i = 0; i < n; ++i)
        if (fixed[i])
            for (std::size_t c = 0; c < dim && c < 3; ++c)
                xyz[i][c] = start[i][c] + target[i][c];

    // The inversion guard, against the input's orientations.
    const MphOrientation before = mph_orientation(rMesh, start);
    const MphOrientation after = mph_orientation(rMesh, xyz);
    for (std::size_t c = 0; c < before.mVolume.size(); ++c) {
        const double v0 = before.mVolume[c], v1 = after.mVolume[c];
        if (std::isfinite(v0) && std::isfinite(v1) && v0 != 0.0 && !(v0 * v1 > 0.0))
            ++result.mNumInverted;
        const MphVec& n0 = before.mNormal[c];
        const MphVec& n1 = after.mNormal[c];
        if (detail::vec3_dot(n0, n0) > 0.0 && !(detail::vec3_dot(n0, n1) > 0.0))
            ++result.mNumInverted;
    }
    if (result.mNumInverted > 0 && !rOptions.mAllowInverted)
        throw std::invalid_argument(
            "morph: " + std::to_string(result.mNumInverted) +
            " cell(s) turned inside out; apply the displacement in more steps, raise the "
            "stiffening, or pass allow_inverted to keep the result");

    NDArray pts(DType::Float64, {n, dim});
    double* pp = pts.As<double>();
    for (std::size_t i = 0; i < n; ++i) {
        double d2 = 0.0;
        for (std::size_t k = 0; k < dim; ++k) {
            pp[i * dim + k] = k < 3 ? xyz[i][k] : detail::read_double(rMesh.Points(), i * dim + k);
            const double dk = k < 3 ? xyz[i][k] - start[i][k] : 0.0;
            d2 += dk * dk;
        }
        result.mMaxDisplacement = std::max(result.mMaxDisplacement, std::sqrt(d2));
    }
    result.mMesh = detail::clone_mesh(rMesh);
    result.mMesh.AssignPoints(std::move(pts));
    return result;
}

}  // namespace meshioplusplus

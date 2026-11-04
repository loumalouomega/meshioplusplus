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
// Level-set cutting and cut-cell quadrature through CutCells (doc/cut.md). The
// mesh is handed over as a CutCells MeshView of its top-dimensional cells
// (corners in VTK order), every level set as a Pk function on those cells, and
// the part CutCells selects comes back as pieces (straight lookup-table pieces
// or curved quadrays leaves) or as quadrature rules. Point data is
// interpolated in each piece's parent cell; cell data and cell regions follow
// the parent. Without CutCells (MESHIOPLUSPLUS_HAS_CUTCELLS unset) both
// operations throw Unsupported naming the CMake option.

// System includes
#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/detail/subset.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/operations/clean.hpp"
#include "meshioplusplus/operations/convert_cells.hpp"
#include "meshioplusplus/operations/cut.hpp"
#include "meshioplusplus/region.hpp"
#include "meshioplusplus/vtk_common.hpp"

#ifdef MESHIOPLUSPLUS_HAS_CUTCELLS
// External includes
#include <cutcells/cell_types.h>
#include <cutcells/level_set.h>
#include <cutcells/lut/cut_mesh.h>
#include <cutcells/mesh_view.h>
#include <cutcells/part/cut_result.h>
#include <cutcells/part/mesh_part.h>
#include <cutcells/part/output.h>
#include <cutcells/quadrature.h>
#include <cutcells/quadrays/leaves.h>
#include <cutcells/reference_cell.h>
#endif

namespace meshioplusplus {

namespace {

constexpr const char* kCutPrefix = "meshio++: cut: ";

#ifdef MESHIOPLUSPLUS_HAS_CUTCELLS

[[noreturn]] void cut_fail(const std::string& rWhy) {
    throw std::invalid_argument(std::string(kCutPrefix) + rWhy);
}

/// A cuttable cell type: its CutCells family (cutcells::cell::type's value),
/// topological dimension, corner count, node count when the level set may be
/// quadratic on it (0: linear only), its corners' VTK type, and whether its
/// meshio corner order differs from VTK's (the wedge).
struct CutKind {
    int mFamily = 0;
    int mTdim = 0;
    int mCorners = 0;
    int mQuadraticNodes = 0;
    int mVtk = 0;
    bool mWedge = false;
};

std::optional<CutKind> cut_kind(const std::string& rType) {
    // cutcells::cell::type: interval 1, triangle 2, tetrahedron 3,
    // quadrilateral 4, hexahedron 5, prism 6, pyramid 7.
    if (rType == "line" || rType == "line3")
        return CutKind{1, 1, 2, 0, 3, false};
    if (rType == "triangle" || rType == "triangle7")
        return CutKind{2, 2, 3, 0, 5, false};
    if (rType == "triangle6")
        return CutKind{2, 2, 3, 6, 5, false};
    if (rType == "quad" || rType == "quad8" || rType == "quad9")
        return CutKind{4, 2, 4, 0, 9, false};
    if (rType == "tetra")
        return CutKind{3, 3, 4, 0, 10, false};
    if (rType == "tetra10")
        return CutKind{3, 3, 4, 10, 10, false};
    if (rType == "hexahedron" || rType == "hexahedron20" || rType == "hexahedron27")
        return CutKind{5, 3, 8, 0, 12, false};
    if (rType == "wedge" || rType == "wedge15" || rType == "wedge18")
        return CutKind{6, 3, 6, 0, 13, true};
    if (rType == "pyramid" || rType == "pyramid13" || rType == "pyramid14")
        return CutKind{7, 3, 5, 0, 14, false};
    return std::nullopt;
}

/// meshio wedge corners -> VTK order (meshio_to_vtk_order("wedge")); its own inverse.
constexpr int kCutWedgeOrder[6] = {0, 2, 1, 3, 5, 4};

/// Per CutCells family, VTK corner j is Basix corner kCutVtkToBasix[f][j] and
/// back (cutcells::cell::vtk_to_basix_vertex_permutation; each is its own
/// inverse). CutCells' straight pieces and whole cells come in Basix order.
constexpr int kCutVtkToBasix[8][8] = {{0},
                                      {0, 1},
                                      {0, 1, 2},
                                      {0, 1, 2, 3},
                                      {0, 1, 3, 2},
                                      {0, 1, 3, 2, 4, 5, 7, 6},
                                      {0, 1, 2, 3, 4, 5},
                                      {0, 1, 3, 2, 4}};

// --- reference maps of the parent cells, for interpolating point data ------
//
// Each family is parametrised over its own reference cell (corners in VTK
// order): the map x(p) = sum N_i(p) X_i with the cell's linear shape
// functions. A piece's vertex is located in its parent by Gauss-Newton, and
// point data is evaluated with the same functions (quadratic ones on the
// triangle6/tetra10 parents of a quadratic level set).

void cut_linear_shape(int Family, const double* pP, double* pN) {
    const double u = pP[0], v = pP[1], w = pP[2];
    switch (Family) {
        case 1:
            pN[0] = 1.0 - u;
            pN[1] = u;
            return;
        case 2:
            pN[0] = 1.0 - u - v;
            pN[1] = u;
            pN[2] = v;
            return;
        case 3:
            pN[0] = 1.0 - u - v - w;
            pN[1] = u;
            pN[2] = v;
            pN[3] = w;
            return;
        case 4:
            pN[0] = (1.0 - u) * (1.0 - v);
            pN[1] = u * (1.0 - v);
            pN[2] = u * v;
            pN[3] = (1.0 - u) * v;
            return;
        case 5: {
            const double a[8] = {(1 - u) * (1 - v), u * (1 - v), u * v, (1 - u) * v,
                                 (1 - u) * (1 - v), u * (1 - v), u * v, (1 - u) * v};
            for (int i = 0; i < 8; ++i)
                pN[i] = a[i] * (i < 4 ? 1.0 - w : w);
            return;
        }
        case 6: {
            const double l[3] = {1.0 - u - v, u, v};
            for (int i = 0; i < 3; ++i) {
                pN[i] = l[i] * (1.0 - w);
                pN[i + 3] = l[i] * w;
            }
            return;
        }
        case 7: {
            // The base bilinear, shrunk towards the apex: a collapsed hexahedron.
            const double b[4] = {(1 - u) * (1 - v), u * (1 - v), u * v, (1 - u) * v};
            for (int i = 0; i < 4; ++i)
                pN[i] = b[i] * (1.0 - w);
            pN[4] = w;
            return;
        }
        default:
            return;
    }
}

int cut_family_tdim(int Family) {
    return Family == 1 ? 1 : (Family == 2 || Family == 4) ? 2 : 3;
}

int cut_family_corners(int Family) {
    static constexpr int kCorners[8] = {1, 2, 3, 4, 4, 8, 6, 5};
    return kCorners[Family];
}

/// The reference coordinates of `pX` in the cell with corners `rX` (gdim
/// each), by Gauss-Newton from the centroid with a difference Jacobian.
std::array<double, 3> cut_locate(int Family, const std::vector<double>& rX, std::size_t Gdim,
                                 const double* pX) {
    const int tdim = cut_family_tdim(Family);
    const int nc = cut_family_corners(Family);
    std::array<double, 3> p = {0.25, 0.25, 0.25};
    if (Family == 4 || Family == 5)
        p = {0.5, 0.5, 0.5};
    if (Family == 1)
        p = {0.5, 0.0, 0.0};
    if (Family == 6 || Family == 7)
        p = {0.25, 0.25, 0.5};
    double n[8];
    const auto map = [&](const std::array<double, 3>& rP, double* pOut) {
        cut_linear_shape(Family, rP.data(), n);
        for (std::size_t k = 0; k < Gdim; ++k) {
            pOut[k] = 0.0;
            for (int i = 0; i < nc; ++i)
                pOut[k] += n[i] * rX[static_cast<std::size_t>(i) * Gdim + k];
        }
    };
    double scale = 0.0;
    for (int i = 1; i < nc; ++i)
        for (std::size_t k = 0; k < Gdim; ++k)
            scale = std::max(scale, std::fabs(rX[static_cast<std::size_t>(i) * Gdim + k] - rX[k]));
    if (scale == 0.0)
        return p;
    double x[3], xp[3], xm[3];
    for (int it = 0; it < 40; ++it) {
        map(p, x);
        double r[3] = {0, 0, 0};
        double rnorm = 0.0;
        for (std::size_t k = 0; k < Gdim; ++k) {
            r[k] = pX[k] - x[k];
            rnorm = std::max(rnorm, std::fabs(r[k]));
        }
        if (rnorm <= 1e-15 * scale)
            break;
        double jac[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};  // jac[k][j] = dx_k / dp_j
        const double h = 1e-6;
        for (int j = 0; j < tdim; ++j) {
            std::array<double, 3> pp = p, pm = p;
            pp[static_cast<std::size_t>(j)] += h;
            pm[static_cast<std::size_t>(j)] -= h;
            map(pp, xp);
            map(pm, xm);
            for (std::size_t k = 0; k < Gdim; ++k)
                jac[k][j] = (xp[k] - xm[k]) / (2.0 * h);
        }
        // Normal equations (J^T J) d = J^T r, tdim <= 3, by Gaussian elimination.
        double a[3][4] = {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};
        for (int i = 0; i < tdim; ++i) {
            for (int j = 0; j < tdim; ++j)
                for (std::size_t k = 0; k < Gdim; ++k)
                    a[i][j] += jac[k][i] * jac[k][j];
            for (std::size_t k = 0; k < Gdim; ++k)
                a[i][3] += jac[k][i] * r[k];
        }
        bool singular = false;
        for (int c = 0; c < tdim; ++c) {
            int piv = c;
            for (int i = c + 1; i < tdim; ++i)
                if (std::fabs(a[i][c]) > std::fabs(a[piv][c]))
                    piv = i;
            if (std::fabs(a[piv][c]) < 1e-300) {
                singular = true;
                break;
            }
            for (int j = 0; j < 4; ++j)
                std::swap(a[c][j], a[piv][j]);
            for (int i = 0; i < tdim; ++i) {
                if (i == c)
                    continue;
                const double f = a[i][c] / a[c][c];
                for (int j = c; j < 4; ++j)
                    a[i][j] -= f * a[c][j];
            }
        }
        if (singular)
            break;
        double step = 0.0;
        for (int i = 0; i < tdim; ++i) {
            const double d = a[i][3] / a[i][i];
            p[static_cast<std::size_t>(i)] += d;
            step = std::max(step, std::fabs(d));
        }
        if (step < 1e-15)
            break;
    }
    return p;
}

/// The weights of the parent's nodes at reference point `rP`: its quadratic
/// Lagrange functions when `Quadratic` (triangle6, tetra10, VTK edge order),
/// else its linear ones on the corners.
void cut_node_weights(int Family, bool Quadratic, const std::array<double, 3>& rP,
                      std::vector<double>& rW) {
    if (!Quadratic) {
        rW.assign(static_cast<std::size_t>(cut_family_corners(Family)), 0.0);
        cut_linear_shape(Family, rP.data(), rW.data());
        return;
    }
    if (Family == 2) {
        const double l[3] = {1.0 - rP[0] - rP[1], rP[0], rP[1]};
        static constexpr int kE[3][2] = {{0, 1}, {1, 2}, {0, 2}};
        rW.assign(6, 0.0);
        for (int i = 0; i < 3; ++i)
            rW[static_cast<std::size_t>(i)] = l[i] * (2.0 * l[i] - 1.0);
        for (int e = 0; e < 3; ++e)
            rW[static_cast<std::size_t>(3 + e)] = 4.0 * l[kE[e][0]] * l[kE[e][1]];
        return;
    }
    const double l[4] = {1.0 - rP[0] - rP[1] - rP[2], rP[0], rP[1], rP[2]};
    static constexpr int kE[6][2] = {{0, 1}, {1, 2}, {0, 2}, {0, 3}, {1, 3}, {2, 3}};
    rW.assign(10, 0.0);
    for (int i = 0; i < 4; ++i)
        rW[static_cast<std::size_t>(i)] = l[i] * (2.0 * l[i] - 1.0);
    for (int e = 0; e < 6; ++e)
        rW[static_cast<std::size_t>(4 + e)] = 4.0 * l[kE[e][0]] * l[kE[e][1]];
}

/// The cells being cut, flattened: the top-dimensional blocks' cells.
struct CutInput {
    std::size_t mGdim = 3;
    int mTdim = 0;
    bool mQuadratic = false;       // every cut cell is triangle6/tetra10
    std::vector<double> mPoints;   // n * mGdim
    std::vector<double> mPoints3;  // n * 3
    std::vector<int> mCorners;     // CSR, VTK corner order
    std::vector<int> mOffsets;     // num cells + 1
    std::vector<int> mFamily;      // per cell
    std::vector<int> mDofs;        // CSR, Basix order, the level set's dofs
    std::vector<int> mDofOffsets;
    std::vector<std::int64_t> mNodes;  // CSR: every node of the cell, meshio order
    std::vector<std::int64_t> mNodeOffsets;
    std::vector<std::int64_t> mGlobal;         // per cell, its global input cell (block-major)
    std::vector<std::size_t> mBlock;           // per cell, its input block
    std::vector<std::size_t> mLocal;           // per cell, its index in that block
    std::vector<std::size_t> mTopBlocks;       // the cut blocks
    std::vector<std::string> mNames;           // level-set names
    std::vector<std::vector<double>> mValues;  // per level set, one value per point
    std::string mSelection;
    std::vector<std::string> mWarnings;
};

std::string cut_backend(const CutOptions& rOptions, bool Quadratic) {
    const std::string& b = rOptions.mBackend;
    if (b == "auto" || b.empty())
        return Quadratic ? "quadrays" : "lut";
    if (b == "lut" || b == "quadrays")
        return b;
    cut_fail("unknown backend '" + b + "' (auto, lut or quadrays)");
}

CutInput cut_prepare(const Mesh& rMesh, const CutOptions& rOptions) {
    CutInput in;
    if (rOptions.mLevelSets.empty() && rOptions.mPlanes.empty())
        cut_fail("give a level set: a point_data array (level_sets) or a plane (planes)");
    if (rOptions.mOrder < 0 || rOptions.mOrder > 20)
        cut_fail("order must be 0 to 20, not " + std::to_string(rOptions.mOrder));
    if (rOptions.mDegree < 1 || rOptions.mDegree > 8)
        cut_fail("degree must be 1 to 8, not " + std::to_string(rOptions.mDegree));
    const std::size_t np = rMesh.NumPoints();
    const std::size_t gdim = rMesh.PointDim();
    if (np == 0)
        cut_fail("the mesh has no points");
    if (gdim < 1 || gdim > 3)
        cut_fail("points must have 1 to 3 coordinates, not " + std::to_string(gdim));
    in.mGdim = gdim;
    const NDArray& r_p = rMesh.Points();
    in.mPoints.resize(np * gdim);
    in.mPoints3.assign(np * 3, 0.0);
    for (std::size_t i = 0; i < np; ++i)
        for (std::size_t k = 0; k < gdim; ++k) {
            const double x = detail::read_double(r_p, i * gdim + k);
            in.mPoints[i * gdim + k] = x;
            in.mPoints3[i * 3 + k] = x;
        }

    // The blocks of the highest dimension a cuttable type reaches.
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto cells = rMesh.Cells(b);
        if (const auto kind = cut_kind(cells.Type()))
            in.mTdim = std::max(in.mTdim, kind->mTdim);
    }
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto cells = rMesh.Cells(b);
        const std::string type = cells.Type();
        const bool polyhedron = cells.IsPolyhedron() || type.rfind("polyhedron", 0) == 0;
        const bool polygon = !polyhedron && type.rfind("polygon", 0) == 0;
        if (cells.NumCells() > 0 && (polyhedron || (polygon && in.mTdim <= 2)))
            throw Unsupported(std::string(kCutPrefix) +
                              "polygon and polyhedron cells are not cut ('" + type +
                              "'); convert_cells(mode='simplexify') first");
    }
    if (in.mTdim == 0)
        cut_fail(
            "no cell to cut: lines, triangles, quads, tetrahedra, hexahedra, wedges or pyramids");
    bool all_quadratic = true;
    std::int64_t global = 0;
    std::set<std::string> skipped;
    in.mOffsets.push_back(0);
    in.mNodeOffsets.push_back(0);
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto cells = rMesh.Cells(b);
        const std::size_t nc = cells.NumCells();
        const auto kind = cut_kind(cells.Type());
        if (!kind || kind->mTdim != in.mTdim || cells.IsRagged()) {
            if (nc > 0 && cells.Type() != "vertex")
                skipped.insert(cells.Type());
            global += static_cast<std::int64_t>(nc);
            continue;
        }
        in.mTopBlocks.push_back(b);
        if (kind->mQuadraticNodes == 0)
            all_quadratic = false;
        const NDArray& r_c = cells.Conn();
        const std::size_t npc = cells.NodesPerCell();
        if (npc < static_cast<std::size_t>(kind->mCorners))
            cut_fail("a '" + cells.Type() + "' block has " + std::to_string(npc) +
                     " nodes per cell");
        for (std::size_t c = 0; c < nc; ++c) {
            int corner[8];
            for (int v = 0; v < kind->mCorners; ++v) {
                const int src = kind->mWedge ? kCutWedgeOrder[v] : v;
                const std::int64_t node =
                    detail::read_int(r_c, c * npc + static_cast<std::size_t>(src));
                if (node < 0 || static_cast<std::size_t>(node) >= np)
                    cut_fail("a '" + cells.Type() + "' cell names node " + std::to_string(node) +
                             " of " + std::to_string(np));
                corner[v] = static_cast<int>(node);
                in.mCorners.push_back(corner[v]);
            }
            in.mOffsets.push_back(static_cast<int>(in.mCorners.size()));
            for (std::size_t k = 0; k < npc; ++k)
                in.mNodes.push_back(detail::read_int(r_c, c * npc + k));
            in.mNodeOffsets.push_back(static_cast<std::int64_t>(in.mNodes.size()));
            in.mFamily.push_back(kind->mFamily);
            in.mGlobal.push_back(global + static_cast<std::int64_t>(c));
            in.mBlock.push_back(b);
            in.mLocal.push_back(c);
        }
        global += static_cast<std::int64_t>(nc);
    }
    if (in.mFamily.empty())
        cut_fail("no cell to cut in the blocks of dimension " + std::to_string(in.mTdim));
    for (const std::string& r_t : skipped)
        in.mWarnings.push_back(std::string(kCutPrefix) + "the '" + r_t +
                               "' cells are not of dimension " + std::to_string(in.mTdim) +
                               " and are left out");
    in.mQuadratic = all_quadratic;

    // The level sets' names and values, at every point.
    for (const std::string& r_name : rOptions.mLevelSets) {
        if (!rMesh.HasPointData(r_name))
            cut_fail("no point_data array '" + r_name + "' to cut by");
        const NDArray& r_a = rMesh.PointData(r_name);
        if (detail::rows(r_a) != np || detail::cols(r_a) != 1)
            cut_fail("the level set '" + r_name + "' must be one value per point");
        std::vector<double> v(np);
        for (std::size_t i = 0; i < np; ++i) {
            v[i] = detail::read_double(r_a, i);
            if (!std::isfinite(v[i]))
                cut_fail("the level set '" + r_name + "' is not finite at point " +
                         std::to_string(i));
        }
        in.mNames.push_back(r_name);
        in.mValues.push_back(std::move(v));
    }
    for (std::size_t q = 0; q < rOptions.mPlanes.size(); ++q) {
        const CutPlane& r_pl = rOptions.mPlanes[q];
        const double len =
            std::sqrt(r_pl.mNormal[0] * r_pl.mNormal[0] + r_pl.mNormal[1] * r_pl.mNormal[1] +
                      r_pl.mNormal[2] * r_pl.mNormal[2]);
        if (!(len > 0.0) || !std::isfinite(len))
            cut_fail("plane " + std::to_string(q) + " has no normal");
        std::vector<double> v(np);
        for (std::size_t i = 0; i < np; ++i) {
            double s = 0.0;
            for (std::size_t k = 0; k < 3; ++k)
                s += (in.mPoints3[i * 3 + k] - r_pl.mOrigin[k]) * r_pl.mNormal[k];
            v[i] = s / len;
        }
        in.mNames.push_back("plane" + std::to_string(q));
        in.mValues.push_back(std::move(v));
    }
    for (std::size_t i = 0; i < in.mNames.size(); ++i) {
        const std::string& r_n = in.mNames[i];
        if (r_n.empty() || !(std::isalpha(static_cast<unsigned char>(r_n[0])) || r_n[0] == '_'))
            cut_fail("the level set '" + r_n +
                     "' needs a name that starts with a letter, to be selected");
        for (const char ch : r_n)
            if (!(std::isalnum(static_cast<unsigned char>(ch)) || ch == '_'))
                cut_fail("the level set '" + r_n +
                         "' needs a name of letters, digits and '_', to be selected");
        for (std::size_t j = 0; j < i; ++j)
            if (in.mNames[j] == r_n)
                cut_fail("two level sets are named '" + r_n + "'");
    }
    in.mSelection = rOptions.mSelection.empty() ? in.mNames.front() + " < 0" : rOptions.mSelection;

    // The level set's dofs on each cell, in Basix order: the corners
    // (permuted from VTK's), or for a quadratic one all ten or six nodes
    // (Basix numbers a P2 tetrahedron's edges 23, 13, 12, 03, 02, 01 and a
    // P2 triangle's 12, 02, 01; VTK 01, 12, 02, 03, 13, 23 and 01, 12, 02).
    static constexpr int kTet10[10] = {0, 1, 2, 3, 9, 8, 5, 7, 6, 4};
    static constexpr int kTri6[6] = {0, 1, 2, 4, 5, 3};
    in.mDofOffsets.push_back(0);
    for (std::size_t c = 0; c < in.mFamily.size(); ++c) {
        const int fam = in.mFamily[c];
        if (in.mQuadratic) {
            const std::int64_t* p_n = in.mNodes.data() + in.mNodeOffsets[c];
            const int n = fam == 2 ? 6 : 10;
            for (int i = 0; i < n; ++i)
                in.mDofs.push_back(static_cast<int>(p_n[fam == 2 ? kTri6[i] : kTet10[i]]));
        } else {
            const int* p_c = in.mCorners.data() + in.mOffsets[c];
            for (int i = 0; i < cut_family_corners(fam); ++i)
                in.mDofs.push_back(p_c[kCutVtkToBasix[fam][i]]);
        }
        in.mDofOffsets.push_back(static_cast<int>(in.mDofs.size()));
    }
    return in;
}

/// meshio name and VTK -> meshio node order of a CutCells piece family.
std::string cut_family_name(int Family) {
    static const char* kNames[8] = {"vertex", "line",       "triangle", "tetra",
                                    "quad",   "hexahedron", "wedge",    "pyramid"};
    return kNames[Family];
}

/// Raw pieces: points (3 per node), cells (CSR, meshio order) with their
/// meshio type and parent (an index into CutInput's cells).
struct CutPieces {
    std::vector<double> mPoints;
    std::vector<std::int64_t> mConn;
    std::vector<std::int64_t> mOffsets = {0};
    std::vector<std::string> mTypes;
    std::vector<std::int64_t> mParents;
};

/// The output mesh of the pieces: welded points, blocks by type, point data
/// interpolated in the parents, cell data, cell regions and field data
/// carried, the cut:* arrays.
Mesh cut_build_mesh(const Mesh& rMesh, const CutInput& rIn, const CutPieces& rPieces,
                    const std::vector<char>& rIsCut, bool RecordParentIds,
                    std::vector<std::string>& rWarnings) {
    const std::size_t gdim = rIn.mGdim;
    const std::size_t npts = rPieces.mPoints.size() / 3;
    const std::size_t ncells = rPieces.mTypes.size();

    // Blocks: one per (type, node count), in order of first appearance.
    std::vector<std::pair<std::string, std::size_t>> keys;
    std::vector<std::vector<std::size_t>> members;
    for (std::size_t c = 0; c < ncells; ++c) {
        const std::size_t nn =
            static_cast<std::size_t>(rPieces.mOffsets[c + 1] - rPieces.mOffsets[c]);
        std::size_t k = 0;
        while (k < keys.size() && !(keys[k].first == rPieces.mTypes[c] && keys[k].second == nn))
            ++k;
        if (k == keys.size()) {
            keys.emplace_back(rPieces.mTypes[c], nn);
            members.emplace_back();
        }
        members[k].push_back(c);
    }

    Mesh out;
    NDArray pts = NDArray::Uninit(DType::Float64, {npts, gdim});
    for (std::size_t i = 0; i < npts; ++i)
        for (std::size_t k = 0; k < gdim; ++k)
            pts.As<double>()[i * gdim + k] = rPieces.mPoints[i * 3 + k];
    out.AssignPoints(std::move(pts));
    std::vector<std::int64_t> order;  // output cells in block order -> piece
    for (std::size_t k = 0; k < keys.size(); ++k) {
        const std::size_t nn = keys[k].second;
        NDArray conn = NDArray::Uninit(DType::Int64, {members[k].size(), nn});
        for (std::size_t r = 0; r < members[k].size(); ++r) {
            const std::size_t c = members[k][r];
            for (std::size_t j = 0; j < nn; ++j)
                conn.As<std::int64_t>()[r * nn + j] =
                    rPieces.mConn[static_cast<std::size_t>(rPieces.mOffsets[c]) + j];
            order.push_back(static_cast<std::int64_t>(c));
        }
        out.AddCellBlock(keys[k].first, std::move(conn));
    }

    // Point data, interpolated in each point's first parent.
    std::vector<std::int64_t> owner(npts, -1);
    for (std::size_t c = 0; c < ncells; ++c)
        for (std::int64_t j = rPieces.mOffsets[c]; j < rPieces.mOffsets[c + 1]; ++j) {
            const std::size_t node =
                static_cast<std::size_t>(rPieces.mConn[static_cast<std::size_t>(j)]);
            if (owner[node] < 0)
                owner[node] = rPieces.mParents[c];
        }
    const std::vector<std::string> pnames = rMesh.PointDataNames();
    if (!pnames.empty()) {
        std::vector<std::vector<double>> weights(npts);
        std::vector<double> corners;
        for (std::size_t i = 0; i < npts; ++i) {
            if (owner[i] < 0)
                continue;
            const std::size_t c = static_cast<std::size_t>(owner[i]);
            const int fam = rIn.mFamily[c];
            corners.clear();
            for (int j = rIn.mOffsets[c]; j < rIn.mOffsets[c + 1]; ++j)
                for (std::size_t k = 0; k < gdim; ++k)
                    corners.push_back(rIn.mPoints[static_cast<std::size_t>(
                                                      rIn.mCorners[static_cast<std::size_t>(j)]) *
                                                      gdim +
                                                  k]);
            const auto ref = cut_locate(fam, corners, gdim, &rPieces.mPoints[i * 3]);
            cut_node_weights(fam, rIn.mQuadratic, ref, weights[i]);
        }
        for (const std::string& r_name : pnames) {
            const NDArray& r_a = rMesh.PointData(r_name);
            const std::size_t comps = detail::cols(r_a);
            const bool is_float = detail::is_float_dtype(r_a.Dtype());
            std::vector<std::size_t> shape = r_a.Shape();
            shape[0] = npts;
            NDArray a = NDArray::Uninit(is_float ? DType::Float64 : r_a.Dtype(), shape);
            for (std::size_t i = 0; i < npts; ++i) {
                const std::size_t c = owner[i] < 0 ? 0 : static_cast<std::size_t>(owner[i]);
                // The parent's nodes the weights refer to: all of them for a
                // quadratic parent, else its corners (meshio order -> VTK order).
                std::vector<std::int64_t> nodes;
                if (rIn.mQuadratic) {
                    for (std::int64_t j = rIn.mNodeOffsets[c]; j < rIn.mNodeOffsets[c + 1]; ++j)
                        nodes.push_back(rIn.mNodes[static_cast<std::size_t>(j)]);
                } else {
                    for (int j = rIn.mOffsets[c]; j < rIn.mOffsets[c + 1]; ++j)
                        nodes.push_back(rIn.mCorners[static_cast<std::size_t>(j)]);
                }
                const std::vector<double>& r_w = weights[i];
                for (std::size_t k = 0; k < comps; ++k) {
                    if (is_float) {
                        double s = 0.0;
                        for (std::size_t j = 0; j < nodes.size() && j < r_w.size(); ++j)
                            s += r_w[j] * detail::read_double(
                                              r_a, static_cast<std::size_t>(nodes[j]) * comps + k);
                        a.As<double>()[i * comps + k] = s;
                    } else {
                        std::size_t best = 0;
                        for (std::size_t j = 1; j < nodes.size() && j < r_w.size(); ++j)
                            if (r_w[j] > r_w[best])
                                best = j;
                        detail::write_int(
                            a, i * comps + k,
                            detail::read_int(r_a,
                                             static_cast<std::size_t>(nodes[best]) * comps + k));
                    }
                }
            }
            out.AddPointData(r_name, std::move(a));
        }
    }

    // Cell data, from each piece's parent.
    for (const std::string& r_name : rMesh.CellDataNames()) {
        const NDArray& r_first = rMesh.CellData(r_name, rIn.mTopBlocks.front());
        bool uniform = true;
        for (const std::size_t b : rIn.mTopBlocks) {
            const NDArray& r_b = rMesh.CellData(r_name, b);
            if (r_b.Dtype() != r_first.Dtype() || detail::cols(r_b) != detail::cols(r_first))
                uniform = false;
        }
        if (!uniform) {
            rWarnings.push_back(
                std::string(kCutPrefix) + "cell_data '" + r_name +
                "' differs in type or width between the cut blocks and is left out");
            continue;
        }
        const std::size_t row_bytes = detail::cols(r_first) * dtype_size(r_first.Dtype());
        std::vector<NDArray> blocks;
        std::size_t at = 0;
        for (std::size_t k = 0; k < keys.size(); ++k) {
            std::vector<std::size_t> shape = r_first.Shape();
            shape[0] = members[k].size();
            NDArray a = NDArray::Uninit(r_first.Dtype(), shape);
            for (std::size_t r = 0; r < members[k].size(); ++r, ++at) {
                const std::size_t parent =
                    static_cast<std::size_t>(rPieces.mParents[static_cast<std::size_t>(order[at])]);
                const NDArray& r_src = rMesh.CellData(r_name, rIn.mBlock[parent]);
                std::memcpy(a.As<std::byte>() + r * row_bytes,
                            r_src.As<std::byte>() + rIn.mLocal[parent] * row_bytes, row_bytes);
            }
            blocks.push_back(std::move(a));
        }
        out.AddCellData(r_name, std::move(blocks));
    }
    {
        std::vector<NDArray> parent_blocks, cut_blocks;
        std::size_t at = 0;
        for (std::size_t k = 0; k < keys.size(); ++k) {
            NDArray pa = NDArray::Uninit(DType::Int64, {members[k].size()});
            NDArray ca = NDArray::Uninit(DType::Int8, {members[k].size()});
            for (std::size_t r = 0; r < members[k].size(); ++r, ++at) {
                const std::size_t parent =
                    static_cast<std::size_t>(rPieces.mParents[static_cast<std::size_t>(order[at])]);
                pa.As<std::int64_t>()[r] = rIn.mGlobal[parent];
                ca.As<std::int8_t>()[r] = rIsCut[parent] ? 1 : 0;
            }
            parent_blocks.push_back(std::move(pa));
            cut_blocks.push_back(std::move(ca));
        }
        if (RecordParentIds)
            out.AddCellData("cut:parent_cell", std::move(parent_blocks));
        out.AddCellData("cut:cut", std::move(cut_blocks));
    }

    // Cell regions follow the parents; point and side regions do not survive.
    std::vector<std::int64_t> parent_global(order.size());
    for (std::size_t i = 0; i < order.size(); ++i)
        parent_global[i] = rIn.mGlobal[static_cast<std::size_t>(
            rPieces.mParents[static_cast<std::size_t>(order[i])])];
    std::string dropped;
    for (std::size_t r = 0; r < rMesh.NumRegions(); ++r) {
        const meshioplusplus::Region& r_reg = rMesh.Region(r);
        if (r_reg.mKind != RegionKind::Cell) {
            dropped += (dropped.empty() ? "'" : ", '") + r_reg.mName + "'";
            continue;
        }
        std::set<std::int64_t> in_reg(r_reg.Entries(), r_reg.Entries() + r_reg.NumEntries());
        std::vector<std::int64_t> kept;
        for (std::size_t i = 0; i < parent_global.size(); ++i)
            if (in_reg.count(parent_global[i]))
                kept.push_back(static_cast<std::int64_t>(i));
        NDArray entries = NDArray::Uninit(DType::Int64, {kept.size()});
        std::copy(kept.begin(), kept.end(), entries.As<std::int64_t>());
        out.AddRegion(meshioplusplus::Region(r_reg.mName, RegionKind::Cell, r_reg.mDim, r_reg.mTag,
                                             std::move(entries)));
    }
    if (!dropped.empty())
        rWarnings.push_back(std::string(kCutPrefix) + "the point and side regions " + dropped +
                            " are not carried to the pieces");
    for (const std::string& r_name : rMesh.FieldDataNames())
        out.AddFieldData(r_name, rMesh.FieldData(r_name));

    // Weld the pieces' coincident points (each piece has its own).
    double diag = 0.0;
    {
        std::array<double, 3> lo = {0, 0, 0}, hi = {0, 0, 0};
        for (std::size_t i = 0; i < rIn.mPoints3.size() / 3; ++i)
            for (std::size_t k = 0; k < 3; ++k) {
                const double x = rIn.mPoints3[i * 3 + k];
                lo[k] = i == 0 ? x : std::min(lo[k], x);
                hi[k] = i == 0 ? x : std::max(hi[k], x);
            }
        for (std::size_t k = 0; k < 3; ++k)
            diag += (hi[k] - lo[k]) * (hi[k] - lo[k]);
        diag = std::sqrt(diag);
    }
    CleanOptions co;
    co.weld = true;
    co.atol = 1e-10 * (diag > 0.0 ? diag : 1.0);
    co.remove_orphans = true;
    co.drop_degenerate = false;
    co.drop_duplicate_cells = false;
    return std::move(clean(out, co).mMesh);
}

using CutLs = cutcells::LevelSetFunction<double, int>;

/// CutCells' view of the cut cells and its classification of them, kept
/// together because the part points into both.
struct CutRun {
    CutInput mIn;
    std::vector<cutcells::cell::type> mTypes;
    cutcells::MeshView<double, int> mView;
    std::vector<CutLs> mLevelSets;
    cutcells::part::CutResult<double, int> mResult;
    cutcells::part::MeshPart<double, int> mPart;
    std::string mBackend;
};

std::unique_ptr<CutRun> cut_run(const Mesh& rMesh, const CutOptions& rOptions) {
    auto run = std::make_unique<CutRun>();
    run->mIn = cut_prepare(rMesh, rOptions);
    CutInput& in = run->mIn;
    run->mBackend = cut_backend(rOptions, in.mQuadratic);
    for (const int f : in.mFamily)
        run->mTypes.push_back(static_cast<cutcells::cell::type>(f));
    auto& view = run->mView;
    view.gdim = static_cast<int>(in.mGdim);
    view.tdim = in.mTdim;
    view.coordinates = std::span<const double>(in.mPoints);
    view.coordinate_stride = in.mGdim;
    view.connectivity = std::span<const int>(in.mCorners);
    view.offsets = std::span<const int>(in.mOffsets);
    view.cell_count = static_cast<int>(in.mFamily.size());
    view.cell_types = std::span<const cutcells::cell::type>(run->mTypes);
    view.vtk_vertex_order = true;
    try {
        for (std::size_t i = 0; i < in.mNames.size(); ++i) {
            auto data = cutcells::create_level_set_mesh_data<double, int>(
                static_cast<int>(in.mGdim), in.mTdim, in.mQuadratic ? 2 : 1,
                std::span<const double>(in.mPoints), std::span<const int>(in.mDofs),
                std::span<const int>(in.mDofOffsets),
                std::span<const cutcells::cell::type>(run->mTypes));
            run->mLevelSets.push_back(cutcells::create_level_set_function<double, int>(
                std::move(data), std::span<const double>(in.mValues[i]), in.mNames[i]));
        }
        run->mResult = cutcells::part::cut<double, int>(
            view, std::span<const CutLs>(run->mLevelSets), cutcells::part::ClassifyOptions{});
    } catch (const std::exception& rE) {
        throw std::runtime_error(std::string(kCutPrefix) +
                                 "CutCells could not classify the cells: " + rE.what());
    }
    try {
        run->mPart = cutcells::part::select(run->mResult, std::string_view(in.mSelection));
    } catch (const std::exception& rE) {
        cut_fail("the selection '" + in.mSelection + "': " + rE.what() + " (level sets: " + [&] {
            std::string s;
            for (const std::string& r_n : in.mNames)
                s += (s.empty() ? "" : ", ") + r_n;
            return s;
        }() + ")");
    }
    return run;
}

cutcells::quadrature::QuadratureRules<double> cut_rules(const CutRun& rRun,
                                                        const CutOptions& rOptions, int Order) {
    if (rRun.mBackend == "lut") {
        cutcells::lut::Options o;
        o.triangulate = rOptions.mTriangulate;
        return cutcells::part::quadrature_rules(rRun.mPart, Order, rOptions.mIncludeUncut, o);
    }
    cutcells::quadrays::Options o;
    return cutcells::part::quadrature_rules(rRun.mPart, Order, rOptions.mIncludeUncut, o);
}

#endif

}  // namespace

bool cut_has_cutcells() {
#ifdef MESHIOPLUSPLUS_HAS_CUTCELLS
    return true;
#else
    return false;
#endif
}

#ifndef MESHIOPLUSPLUS_HAS_CUTCELLS
namespace {
[[noreturn]] void cut_unavailable() {
    throw Unsupported(
        std::string(kCutPrefix) +
        "this build has no CutCells: build meshio++ with -DMESHIOPLUSPLUS_WITH_CUTCELLS=ON "
        "(CutCells, build/build-cutcells.sh), see doc/cut.md");
}
}  // namespace

CutLevelSetResult cut_level_set(const Mesh& rMesh, const CutOptions& rOptions) {
    (void)rMesh;
    (void)rOptions;
    cut_unavailable();
}

CutQuadratureResult cut_quadrature(const Mesh& rMesh, const CutOptions& rOptions) {
    (void)rMesh;
    (void)rOptions;
    cut_unavailable();
}

#else

CutLevelSetResult cut_level_set(const Mesh& rMesh, const CutOptions& rOptions) {
    const auto run = cut_run(rMesh, rOptions);
    const CutInput& in = run->mIn;
    CutLevelSetResult result;
    result.mWarnings = in.mWarnings;
    result.mSelection = in.mSelection;
    result.mBackend = run->mBackend;
    result.mNumCells = static_cast<std::int64_t>(in.mFamily.size());
    result.mNumCutCells = static_cast<std::int64_t>(run->mResult.cut_cells.size());
    result.mNumUncutCells = static_cast<std::int64_t>(run->mPart.n_uncut_cells());
    result.mDimension = run->mPart.dim;
    std::vector<char> is_cut(in.mFamily.size(), 0);
    for (const int c : run->mResult.cut_cells)
        is_cut[static_cast<std::size_t>(c)] = 1;

    CutPieces pieces;
    try {
        if (run->mBackend == "lut") {
            cutcells::lut::Options o;
            o.triangulate = rOptions.mTriangulate;
            const auto cm =
                cutcells::part::visualization_mesh(run->mPart, rOptions.mIncludeUncut, o);
            const std::size_t g = static_cast<std::size_t>(cm._gdim);
            for (int i = 0; i < cm._num_vertices; ++i)
                for (std::size_t k = 0; k < 3; ++k)
                    pieces.mPoints.push_back(
                        k < g ? cm._vertex_coords[static_cast<std::size_t>(i) * g + k] : 0.0);
            for (int c = 0; c < cm._num_cells; ++c) {
                const int fam = static_cast<int>(cm._types[static_cast<std::size_t>(c)]);
                const int begin = cm._offset[static_cast<std::size_t>(c)];
                const int count = cm._offset[static_cast<std::size_t>(c) + 1] - begin;
                for (int j = 0; j < count; ++j) {
                    // meshio position j <- VTK position (the wedge's swap) <- Basix position.
                    const int vtk = fam == 6 ? kCutWedgeOrder[j] : j;
                    const int src =
                        count == cut_family_corners(fam) ? kCutVtkToBasix[fam][vtk] : vtk;
                    pieces.mConn.push_back(cm._connectivity[static_cast<std::size_t>(begin + src)]);
                }
                pieces.mOffsets.push_back(static_cast<std::int64_t>(pieces.mConn.size()));
                pieces.mTypes.push_back(cut_family_name(fam));
                pieces.mParents.push_back(cm._parent_map[static_cast<std::size_t>(c)]);
            }
        } else {
            cutcells::quadrays::Options o;
            const auto lm = cutcells::part::visualization_mesh(run->mPart, rOptions.mDegree,
                                                               rOptions.mIncludeUncut, o);
            pieces.mPoints = lm.points;
            const auto& vmap = vtk_to_meshio_type();
            for (int c = 0; c < lm.n_cells(); ++c) {
                const int vtk = static_cast<int>(lm.vtk_types[static_cast<std::size_t>(c)]);
                const auto it = vmap.find(vtk);
                if (it == vmap.end())
                    throw std::runtime_error(
                        std::string(kCutPrefix) + "CutCells gave a VTK cell of type " +
                        std::to_string(vtk) + ", which meshio++ has no name for");
                const std::vector<int> perm = vtk_to_meshio_order(vtk);
                const int begin = lm.offsets[static_cast<std::size_t>(c)];
                const int count = lm.offsets[static_cast<std::size_t>(c) + 1] - begin;
                for (int j = 0; j < count; ++j) {
                    const int src = perm.empty() ? j : perm[static_cast<std::size_t>(j)];
                    pieces.mConn.push_back(lm.connectivity[static_cast<std::size_t>(begin + src)]);
                }
                pieces.mOffsets.push_back(static_cast<std::int64_t>(pieces.mConn.size()));
                pieces.mTypes.push_back(it->second);
                pieces.mParents.push_back(lm.parent[static_cast<std::size_t>(c)]);
            }
        }
        const auto rules = cut_rules(*run, rOptions, 1);
        for (const double w : rules._weights)
            result.mMeasure += w;
    } catch (const std::invalid_argument&) {
        throw;
    } catch (const std::exception& rE) {
        throw std::runtime_error(std::string(kCutPrefix) + "CutCells failed: " + rE.what());
    }
    result.mNumPieces = static_cast<std::int64_t>(pieces.mTypes.size());
    result.mMesh =
        cut_build_mesh(rMesh, in, pieces, is_cut, rOptions.mRecordParentIds, result.mWarnings);
    if (rOptions.mTriangulate && run->mBackend == "lut") {
        // CutCells' tables split only the interface quadrilaterals: the
        // volume pieces of a hexahedron or a prism, and the cells wholly
        // inside, keep their shape. Split them here; the data, the parent ids
        // and the regions follow each piece to its simplices.
        ConvertCellsOptions co;
        co.mMode = ConvertCellsMode::Simplexify;
        result.mMesh = std::move(convert_cells(result.mMesh, co).mMesh);
    }
    for (const std::string& r_w : result.mWarnings)
        log::warn("{}", r_w);
    return result;
}

CutQuadratureResult cut_quadrature(const Mesh& rMesh, const CutOptions& rOptions) {
    const auto run = cut_run(rMesh, rOptions);
    const CutInput& in = run->mIn;
    CutQuadratureResult result;
    result.mWarnings = in.mWarnings;
    result.mSelection = in.mSelection;
    result.mBackend = run->mBackend;
    result.mDimension = run->mPart.dim;
    cutcells::quadrature::QuadratureRules<double> rules;
    try {
        rules = cut_rules(*run, rOptions, rOptions.mOrder);
    } catch (const std::exception& rE) {
        throw std::runtime_error(std::string(kCutPrefix) + "CutCells failed: " + rE.what());
    }
    const std::size_t n = rules._weights.size();
    const std::size_t tdim = static_cast<std::size_t>(rules._tdim);
    const std::size_t gdim = in.mGdim;
    // The points' physical positions: each rule's reference points (in the
    // parent's Basix reference cell) through the parent's shape functions.
    // The reference cells are this file's parametrisations, except the
    // pyramid's, whose collapsed coordinates come first.
    std::vector<double> phys(n * 3, 0.0);
    {
        double shape[8];
        for (std::size_t r = 0; r < rules._parent_map.size(); ++r) {
            const std::size_t c = static_cast<std::size_t>(rules._parent_map[r]);
            const int fam = in.mFamily[c];
            for (std::int32_t i = rules._offset[r]; i < rules._offset[r + 1]; ++i) {
                const std::size_t ii = static_cast<std::size_t>(i);
                std::array<double, 3> ref = {0.0, 0.0, 0.0};
                for (std::size_t k = 0; k < tdim; ++k)
                    ref[k] = rules._points[ii * tdim + k];
                if (fam == 7)
                    ref = ref[2] < 1.0 - 1e-14
                              ? std::array<double, 3>{ref[0] / (1.0 - ref[2]),
                                                      ref[1] / (1.0 - ref[2]), ref[2]}
                              : std::array<double, 3>{0.0, 0.0, 1.0};
                cut_linear_shape(fam, ref.data(), shape);
                for (int j = in.mOffsets[c]; j < in.mOffsets[c + 1]; ++j) {
                    const double wgt = shape[j - in.mOffsets[c]];
                    const std::size_t node =
                        static_cast<std::size_t>(in.mCorners[static_cast<std::size_t>(j)]);
                    for (std::size_t k = 0; k < 3; ++k)
                        phys[ii * 3 + k] += wgt * in.mPoints3[node * 3 + k];
                }
            }
        }
    }
    Mesh& out = result.mMesh;
    NDArray pts = NDArray::Uninit(DType::Float64, {n, gdim});
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t k = 0; k < gdim; ++k)
            pts.As<double>()[i * gdim + k] = phys[i * 3 + k];
    out.AssignPoints(std::move(pts));
    NDArray verts = NDArray::Uninit(DType::Int64, {n, std::size_t{1}});
    for (std::size_t i = 0; i < n; ++i)
        verts.As<std::int64_t>()[i] = static_cast<std::int64_t>(i);
    out.AddCellBlock("vertex", std::move(verts));
    NDArray w = NDArray::Uninit(DType::Float64, {n});
    NDArray ref = NDArray::Uninit(DType::Float64, {n, tdim});
    NDArray rule = NDArray::Uninit(DType::Int64, {n});
    NDArray parent = NDArray::Uninit(DType::Int64, {n});
    const std::size_t nrules = rules._parent_map.size();
    for (std::size_t r = 0; r < nrules; ++r)
        for (std::int32_t i = rules._offset[r]; i < rules._offset[r + 1]; ++i) {
            const std::size_t ii = static_cast<std::size_t>(i);
            w.As<double>()[ii] = rules._weights[ii];
            result.mSumWeights += rules._weights[ii];
            for (std::size_t k = 0; k < tdim; ++k)
                ref.As<double>()[ii * tdim + k] = rules._points[ii * tdim + k];
            rule.As<std::int64_t>()[ii] = static_cast<std::int64_t>(r);
            parent.As<std::int64_t>()[ii] =
                in.mGlobal[static_cast<std::size_t>(rules._parent_map[r])];
        }
    out.AddPointData("quadrature:weight", std::move(w));
    out.AddPointData("quadrature:reference", std::move(ref));
    out.AddPointData("quadrature:rule", std::move(rule));
    if (rOptions.mRecordParentIds)
        out.AddPointData("quadrature:parent_cell", std::move(parent));
    result.mNumRules = static_cast<std::int64_t>(nrules);
    result.mNumPoints = static_cast<std::int64_t>(n);
    for (const std::string& r_w : result.mWarnings)
        log::warn("{}", r_w);
    return result;
}

#endif

}  // namespace meshioplusplus

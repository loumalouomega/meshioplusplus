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
#include <cmath>
#include <cstddef>
#include <stdexcept>

// Project includes
#include "surface_fem.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/value_io.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

using SfVec = std::array<double, 3>;

SfVec sf_sub(const SfVec& rA, const SfVec& rB) {
    return {rA[0] - rB[0], rA[1] - rB[1], rA[2] - rB[2]};
}

double sf_dot(const SfVec& rA, const SfVec& rB) {
    return rA[0] * rB[0] + rA[1] * rB[1] + rA[2] * rB[2];
}

SfVec sf_cross(const SfVec& rA, const SfVec& rB) {
    return {rA[1] * rB[2] - rA[2] * rB[1], rA[2] * rB[0] - rA[0] * rB[2],
            rA[0] * rB[1] - rA[1] * rB[0]};
}

double sf_norm(const SfVec& rA) {
    return std::sqrt(sf_dot(rA, rA));
}

}  // namespace

SurfaceFem build_surface_fem(const Mesh& rMesh, const std::string& rOp) {
    const std::size_t n = rMesh.NumPoints();
    const std::size_t dim = rMesh.PointDim();

    // --- the surface: every triangle block, nothing else of dimension >= 2 ---
    std::vector<SurfaceFemTriangle> tris;
    for (const auto cb : rMesh.CellRange()) {
        const std::string& type = cb.Type();
        if (cb.IsPolyhedron())
            throw std::invalid_argument(
                rOp + ": a volume mesh; take its surface with extract_surface first");
        const int d = cb.IsRagged() ? 2 : cell_type_dimension(cell_type_from_name(type));
        if (type == "triangle" && !cb.IsRagged()) {
            const NDArray& conn = cb.Conn();
            for (std::size_t c = 0; c < cb.NumCells(); ++c) {
                SurfaceFemTriangle t{};
                for (std::size_t k = 0; k < 3; ++k) {
                    t.mNodes[k] = read_int(conn, c * 3 + k);
                    if (t.mNodes[k] < 0 || static_cast<std::size_t>(t.mNodes[k]) >= n)
                        throw std::invalid_argument(rOp +
                                                    ": a triangle names a point outside the mesh");
                }
                tris.push_back(t);
            }
        } else if (d >= 3) {
            throw std::invalid_argument(rOp + ": '" + type +
                                        "' cells make a volume; take its surface with "
                                        "extract_surface first");
        } else if (d == 2 || type.rfind("polygon", 0) == 0) {
            throw std::invalid_argument(
                rOp + ": '" + type +
                "' cells are not linear triangles; triangulate them with "
                "convert_cells(mode=\"simplexify\") (or \"linearize\" a higher-order one) first");
        }
    }
    if (tris.empty())
        throw std::invalid_argument(rOp + ": the mesh has no triangle");

    SurfaceFem fem;
    const NDArray& points = rMesh.Points();
    fem.mPoints.assign(n, SfVec{0.0, 0.0, 0.0});
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t k = 0; k < dim && k < 3; ++k)
            fem.mPoints[i][k] = read_double(points, i * dim + k);

    // --- cotangents, areas, the mean edge length ---------------------------
    fem.mTriangles.reserve(tris.size());
    double edge_sum = 0.0;
    std::size_t edge_count = 0;
    for (SurfaceFemTriangle t : tris) {
        const SfVec& a = fem.mPoints[t.mNodes[0]];
        const SfVec& b = fem.mPoints[t.mNodes[1]];
        const SfVec& c = fem.mPoints[t.mNodes[2]];
        const double ab = sf_norm(sf_sub(b, a)), bc = sf_norm(sf_sub(c, b)),
                     ca = sf_norm(sf_sub(a, c));
        t.mDoubleArea = sf_norm(sf_cross(sf_sub(b, a), sf_sub(c, a)));
        const double longest = std::max({ab, bc, ca});
        if (!(t.mDoubleArea > 1e-14 * longest * longest)) {
            ++fem.mNumDegenerate;
            continue;
        }
        t.mCot[0] = sf_dot(sf_sub(b, a), sf_sub(c, a)) / t.mDoubleArea;
        t.mCot[1] = sf_dot(sf_sub(a, b), sf_sub(c, b)) / t.mDoubleArea;
        t.mCot[2] = sf_dot(sf_sub(a, c), sf_sub(b, c)) / t.mDoubleArea;
        edge_sum += ab + bc + ca;
        edge_count += 3;
        fem.mTriangles.push_back(t);
    }
    if (fem.mTriangles.empty())
        throw std::invalid_argument(rOp + ": every triangle has zero area");
    fem.mMeanEdge = edge_sum / static_cast<double>(edge_count);
    return fem;
}

void surface_fem_assemble(const SurfaceFem& rFem, std::vector<double>& rMass,
                          std::vector<SparseTriplet>& rStiffness) {
    rMass.assign(rFem.mPoints.size(), 0.0);
    rStiffness.reserve(rStiffness.size() + rFem.mTriangles.size() * 12);
    for (const SurfaceFemTriangle& t : rFem.mTriangles) {
        for (std::size_t k = 0; k < 3; ++k) {
            rMass[t.mNodes[k]] += t.mDoubleArea / 6.0;
            // The edge opposite corner k, weighted by half its cotangent.
            const std::size_t i = static_cast<std::size_t>(t.mNodes[(k + 1) % 3]);
            const std::size_t j = static_cast<std::size_t>(t.mNodes[(k + 2) % 3]);
            const double w = 0.5 * t.mCot[k];
            rStiffness.push_back({i, j, -w});
            rStiffness.push_back({j, i, -w});
            rStiffness.push_back({i, i, w});
            rStiffness.push_back({j, j, w});
        }
    }
}

}  // namespace detail
}  // namespace meshioplusplus

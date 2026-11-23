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
#include <array>
#include <stdexcept>
#include <vector>

// Project includes
#include "meshioplusplus/operations/mesh_operator.hpp"
#include "meshioplusplus/detail/geometry.hpp"

// Project includes (private, not installed)
#include "../detail/sparse_solve.hpp"
#include "../detail/surface_fem.hpp"

namespace meshioplusplus {

MeshOperatorResult mesh_operator(const Mesh& rMesh, const MeshOperatorOptions& rOptions) {
    using detail::SparseTriplet;
    using detail::Vec3;
    const detail::SurfaceFem fem = detail::build_surface_fem(rMesh, "mesh_operator");
    const std::size_t n = fem.mPoints.size();

    std::size_t rows = n;
    std::vector<SparseTriplet> triplets;
    switch (rOptions.mKind) {
        case MeshOperatorKind::Laplacian: {
            std::vector<double> mass;
            detail::surface_fem_assemble(fem, mass, triplets);
            break;
        }
        case MeshOperatorKind::MassLumped: {
            std::vector<double> mass;
            std::vector<SparseTriplet> stiffness;
            detail::surface_fem_assemble(fem, mass, stiffness);
            for (std::size_t i = 0; i < n; ++i)
                if (mass[i] != 0.0)
                    triplets.push_back({i, i, mass[i]});
            break;
        }
        case MeshOperatorKind::MassConsistent: {
            for (const detail::SurfaceFemTriangle& t : fem.mTriangles)
                for (std::size_t a = 0; a < 3; ++a)
                    for (std::size_t b = 0; b < 3; ++b)
                        triplets.push_back({static_cast<std::size_t>(t.mNodes[a]),
                                            static_cast<std::size_t>(t.mNodes[b]),
                                            t.mDoubleArea / 24.0 * (a == b ? 2.0 : 1.0)});
            break;
        }
        case MeshOperatorKind::Gradient: {
            rows = 3 * static_cast<std::size_t>(fem.mNumInput);
            for (const detail::SurfaceFemTriangle& t : fem.mTriangles) {
                const std::array<Vec3, 3> p = {fem.mPoints[t.mNodes[0]], fem.mPoints[t.mNodes[1]],
                                               fem.mPoints[t.mNodes[2]]};
                const Vec3 normal =
                    detail::vec3_cross(detail::vec3_sub(p[1], p[0]), detail::vec3_sub(p[2], p[0]));
                // The unit normal; |normal| is the doubled area.
                const Vec3 unit = {normal[0] / t.mDoubleArea, normal[1] / t.mDoubleArea,
                                   normal[2] / t.mDoubleArea};
                for (std::size_t i = 0; i < 3; ++i) {
                    // grad(phi_i) = n x e_i / (2 A), e_i the opposite edge, counter-clockwise.
                    const Vec3 edge = detail::vec3_sub(p[(i + 2) % 3], p[(i + 1) % 3]);
                    const Vec3 g = detail::vec3_cross(unit, edge);
                    for (std::size_t k = 0; k < 3; ++k)
                        triplets.push_back({3 * static_cast<std::size_t>(t.mIndex) + k,
                                            static_cast<std::size_t>(t.mNodes[i]),
                                            g[k] / t.mDoubleArea});
                }
            }
            break;
        }
        default:
            throw std::invalid_argument("mesh_operator: unknown kind");
    }

    const detail::SparseCsr csr = detail::csr_from_triplets(rows, n, std::move(triplets));
    MeshOperatorResult out;
    out.mRows = static_cast<std::int64_t>(csr.mRows);
    out.mCols = static_cast<std::int64_t>(csr.mCols);
    out.mIndptr.assign(csr.mStart.begin(), csr.mStart.end());
    out.mIndices.assign(csr.mCol.begin(), csr.mCol.end());
    out.mValues = csr.mValue;
    out.mNumDegenerate = fem.mNumDegenerate;
    return out;
}

}  // namespace meshioplusplus

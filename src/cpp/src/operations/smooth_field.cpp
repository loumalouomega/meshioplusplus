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
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/operations/smooth_field.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/operations/data_common.hpp"

// Project includes (private, not installed)
#include "../detail/sparse_solve.hpp"
#include "../detail/surface_fem.hpp"

namespace meshioplusplus {

SmoothFieldResult smooth_field(const Mesh& rMesh, const SmoothFieldOptions& rOptions) {
    const auto fail = [](const std::string& rMessage) {
        throw std::invalid_argument("meshio++: smooth_field: " + rMessage);
    };
    if (rOptions.mArray.empty())
        fail("an array name is required");
    if (!rMesh.HasPointData(rOptions.mArray))
        fail(data_unknown_key_message(rMesh, DataLocation::Point, rOptions.mArray));
    if (!(rOptions.mLength >= 0.0) || !std::isfinite(rOptions.mLength))
        fail("the length must be finite and not negative");
    if (!(rOptions.mTolerance > 0.0))
        fail("the tolerance must be positive");
    const NDArray& field = rMesh.PointData(rOptions.mArray);
    const std::size_t n = rMesh.NumPoints();
    const std::size_t comps = data_num_components(field);
    if (detail::rows(field) != n)
        fail("'" + rOptions.mArray + "' has " + std::to_string(detail::rows(field)) +
             " rows on a mesh of " + std::to_string(n) + " points");

    const detail::SurfaceFem fem = detail::build_surface_fem(rMesh, "smooth_field");
    std::vector<double> mass;
    std::vector<detail::SparseTriplet> triplets;
    detail::surface_fem_assemble(fem, mass, triplets);
    const double l2 = rOptions.mLength * rOptions.mLength;
    for (detail::SparseTriplet& t : triplets)
        t.mValue *= l2;
    for (std::size_t i = 0; i < n; ++i)
        // A point in no triangle has no mass and no stiffness: an identity row keeps it.
        triplets.push_back({i, i, mass[i] > 0.0 ? mass[i] : 1.0});
    const detail::SparseCsr a = detail::csr_from_triplets(n, n, std::move(triplets));

    NDArray out(DType::Float64, field.Shape());
    double* po = out.As<double>();
    std::vector<double> rhs(n), x(n);
    SmoothFieldResult result;
    const std::size_t cap =
        rOptions.mMaxIterations > 0 ? static_cast<std::size_t>(rOptions.mMaxIterations) : n + 100;
    for (std::size_t c = 0; c < comps; ++c) {
        for (std::size_t i = 0; i < n; ++i) {
            const double f = detail::read_double(field, i * comps + c);
            rhs[i] = (mass[i] > 0.0 ? mass[i] : 1.0) * f;
            x[i] = f;
        }
        const detail::SparseSolveReport report =
            detail::csr_pcg(a, rhs.data(), x.data(), rOptions.mTolerance, cap);
        result.mIterations += static_cast<std::int64_t>(report.mIterations);
        result.mConverged = result.mConverged && report.mConverged;
        for (std::size_t i = 0; i < n; ++i)
            po[i * comps + c] = x[i];
    }
    result.mNumDegenerate = fem.mNumDegenerate;
    result.mMesh = detail::clone_mesh(rMesh);
    result.mMesh.AddPointData(
        rOptions.mOutput.empty() ? rOptions.mArray + kSmoothFieldSuffix : rOptions.mOutput,
        std::move(out));
    return result;
}

}  // namespace meshioplusplus

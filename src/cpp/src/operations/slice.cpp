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
// Planar cross-section of a mesh: the signed distance to the plane, handed to
// the shared marching-tetrahedra cutter. Everything below the distance is in
// `detail/marching.hpp`, hoisted out of this file in v7.15.0 (output unchanged,
// byte for byte) so `operations/isosurface.cpp` can cut on a field instead. See
// operations/slice.hpp for the contract, in particular the degeneracy/winding/
// determinism rules the pure-numpy twin (_slice.py) replicates expression for
// expression.

// System includes
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

// Project includes
#include "meshioplusplus/operations/slice.hpp"
#include "meshioplusplus/detail/region_remap.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/marching.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {

Mesh slice(const Mesh& rMesh, const SliceOptions& rOptions) {
    using detail::Vec3;

    const Vec3 origin = {rOptions.mOrigin[0], rOptions.mOrigin[1], rOptions.mOrigin[2]};
    const Vec3 normal = {rOptions.mNormal[0], rOptions.mNormal[1], rOptions.mNormal[2]};
    if (detail::vec3_norm_sq(normal) <= 0.0)
        throw std::invalid_argument("meshio++: slice: normal must be non-zero");

    // Simplexify: every 3D cell becomes a tetra, every 2D cell a triangle.
    const detail::MarchingInput input = detail::marching_prepare(rMesh);
    const Mesh& simp = input.Simp();
    const std::size_t dim = simp.PointDim();

    // Signed distance of every simplexified node to the plane (the hot loop).
    const std::size_t nnodes_in = simp.NumPoints();
    std::vector<double> dist(nnodes_in, 0.0);
    {
        const NDArray& points = simp.Points();
        parallel_for(nnodes_in, [&](std::size_t g) {
            dist[g] = detail::vec3_dot(
                detail::vec3_sub(detail::read_point(points, dim, static_cast<std::int64_t>(g)),
                                 origin),
                normal);
        });
    }

    detail::MarchingOptions options;
    options.mProvenanceName = "slice:parent_cell";
    options.mRecordParentIds = rOptions.mRecordParentIds;
    // A plane has a global orientation, so every section face is wound toward
    // the user's +normal side (an isosurface has none and uses FieldGradient).
    options.mOrientation = detail::MarchingOrientation::FixedDirection;
    options.mDirection = rOptions.mNormal;

    detail::warn_regions_dropped(rMesh, "slice");
    return detail::marching_cut(input, dist, options);
}

}  // namespace meshioplusplus

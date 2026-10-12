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
 * @file nurbs.hpp
 * @brief `NurbsInfo`: the NURBS patches a reader sampled into Lagrange cells,
 * kept beside the `Mesh` (the `GmshInfo`/`MdpaInfo` side-channel pattern).
 *
 * A `Mesh` holds linear and Lagrange polytopes only, so a reader of a NURBS
 * file (MFEM's `MFEM NURBS mesh`, since v16.40.0) evaluates each patch at the
 * nodes of its knot-span elements. The overload that takes a `NurbsInfo` also
 * returns what the sampling drops -- every patch's degrees, knot vectors,
 * control net and weights, the control values of the NURBS fields read with
 * it, the file's patch topology, and the map from each sampled cell back to its
 * patch and knot span -- and a writer that knows NURBS (MFEM's) writes the
 * patches back from it. The shared registry passes none: a plain read drops it
 * and a plain write writes the sampled cells. Python, WASM and the C API's
 * read/write-with-info paths (also Fortran, Julia and R) carry it explicitly.
 *
 * Nothing here is evaluated or checked when the struct is built; the writer
 * validates it. Every index is 0-based. See doc/formats/mfem.md.
 */

// System includes
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/ndarray.hpp"

namespace meshioplusplus {

/**
 * @brief One tensor-product NURBS patch of parametric dimension
 * `mDegrees.size()` (1, 2 or 3).
 *
 * The control net has `n0 * n1 * n2` points, `nd = mKnots[d].size() -
 * mDegrees[d] - 1` per direction, listed with direction 0 fastest. The knot
 * vectors run along the patch's own parametric directions (MFEM orients each
 * patch's directions by its first vertex: the edges 0-1, 0-3 and 0-4).
 */
struct MESHIOPLUSPLUS_API NurbsPatch {
    /// Polynomial degree per parametric direction.
    std::vector<std::int32_t> mDegrees;
    /// The knot vector of each direction (non-decreasing, `n + degree + 1` knots).
    std::vector<std::vector<double>> mKnots;
    /// Float64 `(n0 * n1 * n2, space dimension)`: Cartesian control points.
    NDArray mControlPoints{DType::Float64, {0, 3}};
    /// Float64 `(n0 * n1 * n2)`: the weight of each control point.
    NDArray mWeights{DType::Float64, {0}};
    /// Float64 `(n0 * n1 * n2, components)` per `NurbsInfo::mFieldNames` entry:
    /// the field's control values (its NURBS degrees of freedom).
    std::vector<NDArray> mFields;
};

/**
 * @brief The patch topology of an MFEM NURBS mesh, as its file states it
 * (`mDimension` 0 when the source has none).
 *
 * Patch and boundary vertices are MFEM's topological vertices `0 ..
 * mNumVertices - 1`, which are also the first `mNumVertices` points of the
 * sampled mesh. Boundary rows are in the orientation the reader normalises
 * them to, so a rewrite is idempotent.
 */
struct MESHIOPLUSPLUS_API NurbsTopology {
    /// Parametric dimension of the patches (1, 2 or 3); 0 for no topology.
    std::int32_t mDimension = 0;
    /// Number of topological vertices (patch corners).
    std::int64_t mNumVertices = 0;
    /// Int64 `(patches, 2^dimension)`: each patch's corner vertices, MFEM order.
    NDArray mPatchVertices{DType::Int64, {0, 0}};
    /// Int64 `(patches)`: each patch's attribute.
    NDArray mPatchAttributes{DType::Int64, {0}};
    /// Int64 `(boundary patches, 2^(dimension - 1))`: their corner vertices.
    NDArray mBoundaryVertices{DType::Int64, {0, 0}};
    /// Int64 `(boundary patches)`: each boundary patch's attribute.
    NDArray mBoundaryAttributes{DType::Int64, {0}};
    /// Int64 `(edges, 3)`: the `edges` section, `(knot vector, vertex, vertex)`.
    NDArray mEdges{DType::Int64, {0, 3}};
};

/**
 * @brief The NURBS description of a sampled mesh: patches, fields, topology
 * and the cell-to-patch map. Empty (`mPatches` empty) for a mesh that is not
 * NURBS.
 */
struct MESHIOPLUSPLUS_API NurbsInfo {
    /// The patches, in file order.
    std::vector<NurbsPatch> mPatches;
    /// The NURBS fields carried per patch (`NurbsPatch::mFields`), by name.
    std::vector<std::string> mFieldNames;
    /// The patch topology the file states.
    NurbsTopology mTopology;
    /// Int64 `(element cells)`: the patch of each cell of the first block.
    NDArray mCellPatch{DType::Int64, {0}};
    /// Int64 `(element cells, 3)`: each such cell's knot span per direction,
    /// counted over the spans of nonzero length (0 for an unused direction).
    NDArray mCellSpan{DType::Int64, {0, 3}};
    /// Int64 `(boundary cells)`: the boundary patch of each cell after the
    /// first block, in block order.
    NDArray mBoundaryCellPatch{DType::Int64, {0}};
};

}  // namespace meshioplusplus

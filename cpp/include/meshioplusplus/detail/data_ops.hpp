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
 * @file data_ops.hpp
 * @brief Header-only helpers shared by the *data* operations (`data_manage`,
 * `data_average`, `data_calc`, `data_condition`, `data_info`).
 *
 * The centrepiece is `clone_mesh`, which rebuilds a mesh through the uniform
 * API while a caller-supplied filter decides, per data array, whether it is
 * kept and under what name. Every data operation is "clone the geometry
 * verbatim, then rewrite some data arrays", so this is the one place that has
 * to get the geometry copy right — including the ragged/polyhedron cases.
 *
 * A fresh `Mesh` is built rather than mutated because the uniform mesh API is
 * strictly additive (no backend offers a remove-or-rename for data arrays) and
 * the KRATOS backend's `Mesh` is not copy-constructible, so `Mesh out = in;`
 * does not compile. This mirrors every other operation, all of which return a
 * new mesh.
 *
 * These are `detail::` inline helpers, so they are exempt from the
 * unique-prefix rule the anonymous-namespace helpers in `cpp/src/**.cpp` follow
 * for the single-header amalgamation.
 */

// System includes
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/cell_faces.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/operations/data_common.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {
namespace detail {

/// A deep copy of an `NDArray` that always owns its buffer (safe even when the
/// source is a view over foreign memory, as it is on the Python boundary).
inline NDArray data_owned_copy(const NDArray& rArray) {
    NDArray c = rArray;
    c.MakeOwned();
    return c;
}

/**
 * @brief Copies only the geometry of @p rMesh — points and every cell block,
 * including ragged polygon and polyhedron blocks — into a fresh mesh.
 *
 * No data arrays are carried; the caller adds whichever it wants afterwards.
 * The polyhedron branch is tested before the ragged one because a polyhedron
 * block is also ragged.
 * @param rMesh the mesh whose geometry is copied.
 * @return a new mesh with identical geometry and no data.
 */
inline Mesh clone_geometry(const Mesh& rMesh) {
    Mesh out;
    out.AssignPoints(data_owned_copy(rMesh.Points()));
    for (const auto cb : rMesh.CellRange()) {
        if (cb.IsPolyhedron()) {
            std::vector<std::vector<std::vector<std::int64_t>>> cells(cb.NumCells());
            for (std::size_t c = 0; c < cb.NumCells(); ++c) {
                cells[c].resize(cb.NumFaces(c));
                for (std::size_t f = 0; f < cb.NumFaces(c); ++f) {
                    auto face = cb.Face(c, f);
                    cells[c][f].assign(face.first, face.first + face.second);
                }
            }
            out.AddPolyhedronBlock(std::string(cb.Type()), std::move(cells));
        } else if (cb.IsRagged()) {
            std::vector<std::vector<std::int64_t>> rows(cb.NumCells());
            for (std::size_t c = 0; c < cb.NumCells(); ++c)
                rows[c].assign(cb.Row(c), cb.Row(c) + cb.RowSize(c));
            out.AddPolygonBlock(std::string(cb.Type()), std::move(rows));
        } else {
            out.AddCellBlock(std::string(cb.Type()), data_owned_copy(cb.Conn()));
        }
    }
    return out;
}

/**
 * @brief Clones @p rMesh, letting @p rFilter decide the fate of each data array.
 *
 * The filter is invoked once per array as
 * `bool(DataLocation location, const std::string& name, std::string& newName)`.
 * Returning `false` drops the array; returning `true` keeps it under
 * `newName`, which starts out equal to `name` and may be rewritten in place to
 * rename. Geometry is always copied verbatim.
 *
 * Arrays are visited in the sorted order the uniform API guarantees, so the
 * result is byte-identical across backends.
 * @tparam TFilter the filter callable.
 * @param rMesh the mesh to clone.
 * @param rFilter the per-array keep/rename decision.
 * @return the rewritten mesh.
 */
template <class TFilter>
Mesh clone_mesh(const Mesh& rMesh, TFilter&& rFilter) {
    Mesh out = clone_geometry(rMesh);
    for (const std::string& name : rMesh.PointDataNames()) {
        std::string target = name;
        if (rFilter(DataLocation::Point, name, target))
            out.AddPointData(target, data_owned_copy(rMesh.PointData(name)));
    }
    for (const std::string& name : rMesh.CellDataNames()) {
        std::string target = name;
        if (rFilter(DataLocation::Cell, name, target)) {
            std::vector<NDArray> blocks;
            blocks.reserve(rMesh.CellDataNumBlocks(name));
            for (std::size_t b = 0; b < rMesh.CellDataNumBlocks(name); ++b)
                blocks.push_back(data_owned_copy(rMesh.CellData(name, b)));
            out.AddCellData(target, std::move(blocks));
        }
    }
    for (const std::string& name : rMesh.FieldDataNames()) {
        std::string target = name;
        if (rFilter(DataLocation::Field, name, target))
            out.AddFieldData(target, data_owned_copy(rMesh.FieldData(name)));
    }
    return out;
}

/// Clones @p rMesh whole — geometry plus every data array, unchanged.
inline Mesh clone_mesh(const Mesh& rMesh) {
    return clone_mesh(rMesh, [](DataLocation, const std::string&, std::string&) { return true; });
}

/**
 * @brief Running min / max / sum over the *finite* values of a data array.
 *
 * Non-finite values are counted but never contribute to the reduction, which
 * is the policy documented in `operations/data_common.hpp`. Instances combine
 * associatively via `Merge`, so a parallel chunked reduction can fold per-chunk
 * instances serially afterwards and stay deterministic.
 */
struct FiniteStats {
    double mMin = 0.0;            ///< Smallest finite value (valid iff mNumFinite > 0).
    double mMax = 0.0;            ///< Largest finite value (valid iff mNumFinite > 0).
    double mSum = 0.0;            ///< Sum of the finite values.
    double mSumSq = 0.0;          ///< Sum of the squares of the finite values.
    std::int64_t mNumFinite = 0;  ///< Count of finite values seen.
    std::int64_t mNumNan = 0;     ///< Count of NaN values seen.
    std::int64_t mNumInf = 0;     ///< Count of +/-inf values seen.

    /// Folds one value in.
    void Add(double v) {
        if (std::isnan(v)) {
            ++mNumNan;
            return;
        }
        if (std::isinf(v)) {
            ++mNumInf;
            return;
        }
        if (mNumFinite == 0) {
            mMin = v;
            mMax = v;
        } else {
            if (v < mMin)
                mMin = v;
            if (v > mMax)
                mMax = v;
        }
        mSum += v;
        mSumSq += v * v;
        ++mNumFinite;
    }

    /// Folds another (independently accumulated) instance in.
    void Merge(const FiniteStats& rOther) {
        if (rOther.mNumFinite > 0) {
            if (mNumFinite == 0) {
                mMin = rOther.mMin;
                mMax = rOther.mMax;
            } else {
                if (rOther.mMin < mMin)
                    mMin = rOther.mMin;
                if (rOther.mMax > mMax)
                    mMax = rOther.mMax;
            }
            mSum += rOther.mSum;
            mSumSq += rOther.mSumSq;
            mNumFinite += rOther.mNumFinite;
        }
        mNumNan += rOther.mNumNan;
        mNumInf += rOther.mNumInf;
    }

    /// Mean of the finite values, or NaN when there were none.
    double Mean() const {
        return mNumFinite > 0 ? mSum / static_cast<double>(mNumFinite) : std::nan("");
    }

    /// Population standard deviation (1/N) of the finite values, or NaN.
    double StdDev() const {
        if (mNumFinite <= 0)
            return std::nan("");
        const double n = static_cast<double>(mNumFinite);
        const double mean = mSum / n;
        const double var = mSumSq / n - mean * mean;
        return var > 0.0 ? std::sqrt(var) : 0.0;
    }

    /// Smallest finite value, or NaN when there were none.
    double Min() const { return mNumFinite > 0 ? mMin : std::nan(""); }

    /// Largest finite value, or NaN when there were none.
    double Max() const { return mNumFinite > 0 ? mMax : std::nan(""); }
};

/**
 * @brief Reduces @p rArray into per-component `FiniteStats`.
 *
 * Chunked with `parallel_for` and combined serially, mirroring the bounding-box
 * reduction in `operations/stats.cpp`, so the result does not depend on the
 * thread count.
 * @param rArray the array to reduce.
 * @param NumComponents its component count (see `data_num_components`).
 * @param rStats per-component accumulators, resized to @p NumComponents and
 *        *folded into* (not reset), so several arrays can share one reduction.
 */
inline void accumulate_stats(const NDArray& rArray, std::size_t NumComponents,
                             std::vector<FiniteStats>& rStats) {
    if (rStats.size() < NumComponents)
        rStats.resize(NumComponents);
    const std::size_t total = rArray.Size();
    if (total == 0 || NumComponents == 0)
        return;
    const std::size_t nrows = total / NumComponents;

    const std::size_t grain = 4096;
    const std::size_t nchunks = (nrows + grain - 1) / grain;
    std::vector<std::vector<FiniteStats>> partial(nchunks);
    parallel_for(
        nchunks,
        [&](std::size_t ci) {
            std::vector<FiniteStats> local(NumComponents);
            const std::size_t begin = ci * grain;
            const std::size_t end = std::min(begin + grain, nrows);
            for (std::size_t r = begin; r < end; ++r)
                for (std::size_t k = 0; k < NumComponents; ++k)
                    local[k].Add(read_double(rArray, r * NumComponents + k));
            partial[ci] = std::move(local);
        },
        1);
    for (const std::vector<FiniteStats>& chunk : partial)
        for (std::size_t k = 0; k < NumComponents && k < chunk.size(); ++k)
            rStats[k].Merge(chunk[k]);
}

/// Collapses per-component stats into one whole-array accumulator.
inline FiniteStats combine_components(const std::vector<FiniteStats>& rStats) {
    FiniteStats all;
    for (const FiniteStats& s : rStats)
        all.Merge(s);
    return all;
}

/**
 * @brief Unsigned area of a corner polygon (triangle / quad) via the Newell
 * normal. Mirrors `stats_area` in `operations/stats.cpp`.
 * @param rCoords the corner coordinates.
 * @param Corners how many of them form the polygon.
 * @return the unsigned area, or 0 for fewer than 3 corners.
 */
inline double polygon_area(const std::vector<Vec3>& rCoords, int Corners) {
    if (Corners < 3)
        return 0.0;
    Vec3 s = {0, 0, 0};
    for (int i = 0; i < Corners; ++i)
        s = vec3_add(s, vec3_cross(rCoords[i], rCoords[(i + 1) % Corners]));
    return 0.5 * vec3_norm(s);
}

/**
 * @brief Signed volume of a 3D cell via the divergence theorem over its
 * outward-wound boundary faces. Mirrors `stats_signed_volume`.
 * @param rCoords the cell's corner coordinates.
 * @param Type the cell type, which supplies the face table.
 * @return the signed volume, or NaN for a type with no face table.
 */
inline double cell_signed_volume(const std::vector<Vec3>& rCoords, CellType Type) {
    const std::vector<CellFaceDef>& faces = cell_faces(Type);
    if (faces.empty())
        return std::nan("");
    double vol6 = 0.0;
    for (const CellFaceDef& f : faces) {
        const Vec3 a = rCoords[f.mNodes[0]];
        for (int i = 1; i + 1 < f.mNumCorners; ++i)
            vol6 += triple_product(a, rCoords[f.mNodes[i]], rCoords[f.mNodes[i + 1]]);
    }
    return vol6 / 6.0;
}

/**
 * @brief The |measure| of one cell — length for a 1D cell, area for a 2D cell,
 * volume for a 3D one — used to weight the cell-to-point average.
 *
 * Returns NaN when the measure is not computable (a ragged or polyhedron block,
 * or a type with no face table); callers fall back to a unit weight.
 * @param rPoints the mesh points.
 * @param PointDim the point dimension.
 * @param rCell the cell block view.
 * @param Index the cell index within the block.
 * @return the unsigned measure, or NaN.
 */
inline double cell_measure(const NDArray& rPoints, std::size_t PointDim,
                           const Mesh::CellView& rCell, std::size_t Index) {
    if (rCell.IsRagged() || rCell.IsPolyhedron())
        return std::nan("");
    const CellType ct = cell_type_from_name(rCell.Type());
    const int corners = cell_corner_count(ct);
    if (corners <= 0)
        return std::nan("");
    std::vector<Vec3> coords;
    read_corner_coords(rPoints, PointDim, rCell.Conn(), Index * rCell.NodesPerCell(),
                       static_cast<std::size_t>(corners), coords);
    const int dim = cell_type_dimension(ct);
    if (dim == 3) {
        const double v = cell_signed_volume(coords, ct);
        return std::isnan(v) ? v : std::fabs(v);
    }
    if (dim == 2)
        return polygon_area(coords, corners);
    if (dim == 1 && corners >= 2)
        return vec3_norm(vec3_sub(coords[1], coords[0]));
    return std::nan("");
}

}  // namespace detail
}  // namespace meshioplusplus

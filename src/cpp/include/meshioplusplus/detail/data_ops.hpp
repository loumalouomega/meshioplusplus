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
 * `clone_mesh<TFilter>` is a template (the filter is a caller-supplied lambda
 * type) and so must stay defined here; `FiniteStats::Add`/`Merge` are called
 * per element inside `accumulate_stats`'s inner loop, so that struct stays
 * inline too. Every other function is called once per array or once per cell
 * (never once per scalar), so those bodies live in
 * `src/cpp/src/detail/data_ops.cpp` instead.
 *
 * These are `detail::` helpers, so they are exempt from the unique-prefix rule
 * the anonymous-namespace helpers in `src/cpp/src/**.cpp` follow for the
 * single-header amalgamation.
 */

// System includes
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/operations/data_common.hpp"

namespace meshioplusplus {
namespace detail {

/// A deep copy of an `NDArray` that always owns its buffer (safe even when the
/// source is a view over foreign memory, as it is on the Python boundary).
MESHIOPLUSPLUS_API NDArray data_owned_copy(const NDArray& rArray);

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
MESHIOPLUSPLUS_API Mesh clone_geometry(const Mesh& rMesh);

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
 *
 * `Add`/`Merge` are called once per element (respectively once per chunk)
 * inside `accumulate_stats`'s hot loop, so they stay inline here rather than
 * moving to a `.cpp` -- at that call frequency the function-call boundary
 * itself would be measurable.
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
MESHIOPLUSPLUS_API void accumulate_stats(const NDArray& rArray, std::size_t NumComponents,
                      std::vector<FiniteStats>& rStats);

/// Collapses per-component stats into one whole-array accumulator.
MESHIOPLUSPLUS_API FiniteStats combine_components(const std::vector<FiniteStats>& rStats);

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
MESHIOPLUSPLUS_API double cell_measure(const NDArray& rPoints, std::size_t PointDim, const Mesh::CellView& rCell,
                    std::size_t Index);

/**
 * @brief A per-component weighted-sum accumulator: `sum(value * weight)`,
 * `sum(weight)`, and a count of rows excluded because the *value* at this
 * component was non-finite (`Add` is only ever called for a row whose
 * weight is already known finite and positive -- geometric exclusion, e.g.
 * an unmeasurable cell, is the caller's job, exactly as `FiniteStats`
 * assumes nothing about where its own values came from).
 */
struct WeightedSum {
    double mSumVW = 0.0;           ///< sum(value * weight) over finite values.
    double mSumW = 0.0;            ///< sum(weight) over rows with a finite value.
    std::int64_t mNumNan = 0;      ///< Rows skipped because the value was non-finite.

    /// Folds one (value, weight) pair in; `weight` must already be finite and > 0.
    void Add(double Value, double Weight) {
        if (!std::isfinite(Value)) {
            ++mNumNan;
            return;
        }
        mSumVW += Value * Weight;
        mSumW += Weight;
    }

    /// Folds another (independently accumulated) instance in.
    void Merge(const WeightedSum& rOther) {
        mSumVW += rOther.mSumVW;
        mSumW += rOther.mSumW;
        mNumNan += rOther.mNumNan;
    }

    /// The weighted mean, or NaN when nothing finite contributed.
    double Mean() const { return mSumW > 0.0 ? mSumVW / mSumW : std::nan(""); }
};

/**
 * @brief Reduces @p rArray into per-component `WeightedSum`, weighting row
 * `r` by @p rWeights[r].
 *
 * A row whose weight is not finite and positive contributes to no component
 * at all (silently skipped here -- callers that need to count such rows,
 * e.g. as "geometrically unmeasurable", do so once themselves, since that
 * count is the same for every component and every array sharing the same
 * weights).
 *
 * Chunked with `parallel_for` and combined serially, mirroring
 * `accumulate_stats`, so the result does not depend on the thread count.
 * @param rArray the array to reduce; row count must equal `rWeights.size()`.
 * @param NumComponents its component count (see `data_num_components`).
 * @param rWeights one weight per row.
 * @param rStats per-component accumulators, resized to @p NumComponents and
 *        *folded into* (not reset), so several arrays can share one reduction.
 */
MESHIOPLUSPLUS_API void accumulate_weighted(const NDArray& rArray, std::size_t NumComponents,
                        const std::vector<double>& rWeights, std::vector<WeightedSum>& rStats);

}  // namespace detail
}  // namespace meshioplusplus

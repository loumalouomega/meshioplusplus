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
// Value conditioning for data arrays: clamp / normalize / standardize, either
// per component or by row magnitude. Statistics come from the shared
// detail::FiniteStats chunked reducer (deterministic: parallel chunks, serial
// combine) and, for cell_data, are gathered across ALL blocks before a single
// transform is applied to each -- normalizing one block against its own
// extremes would make the pieces mutually incomparable.
//
// Geometry is never modified. See operations/data_condition.hpp for the
// contract.

// System includes
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/data_condition.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/operations/data_common.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {

namespace {

/// The affine map `y = scale * x + offset` a mode/scope resolves to, per
/// component (or a single entry when conditioning by magnitude).
struct DcondTransform {
    double mScale = 1.0;
    double mOffset = 0.0;
    bool mClamp = false;  ///< Clamp mode: saturate into [mLo, mHi] instead.
    double mLo = 0.0;
    double mHi = 1.0;
};

/// Applies the NaN policy to one produced value.
double dcond_apply_nan_policy(double v, const DataConditionOptions& rOpts, const std::string& rName,
                              std::size_t index) {
    if (std::isfinite(v))
        return v;
    switch (rOpts.nan_policy) {
        case NanPolicy::Ignore:
            return v;
        case NanPolicy::Replace:
            return rOpts.nan_replacement;
        case NanPolicy::Fail:
            throw std::invalid_argument("meshio++: data_condition: non-finite value at index " +
                                        std::to_string(index) + " of '" + rName + "'");
    }
    return v;
}

/// Builds the transform for one component from its statistics.
DcondTransform dcond_build(const detail::FiniteStats& rStats, const DataConditionOptions& rOpts,
                           const std::string& rName, bool& rWarnedConstant) {
    DcondTransform t;
    switch (rOpts.mode) {
        case ConditionMode::Clamp:
            t.mClamp = true;
            t.mLo = rOpts.lo;
            t.mHi = rOpts.hi;
            return t;
        case ConditionMode::Normalize: {
            const double lo = rStats.Min();
            const double hi = rStats.Max();
            if (rStats.mNumFinite == 0 || !(hi > lo)) {
                // Constant (or entirely non-finite) array: everything maps to
                // the target lower bound. Warn once per array.
                if (!rWarnedConstant) {
                    rWarnedConstant = true;
                    log::warn(
                        "data_condition: '{}' has no value range to normalize (constant or "
                        "non-finite); filling with the target lower bound",
                        rName);
                }
                t.mScale = 0.0;
                t.mOffset = rOpts.lo;
                return t;
            }
            t.mScale = (rOpts.hi - rOpts.lo) / (hi - lo);
            t.mOffset = rOpts.lo - t.mScale * lo;
            return t;
        }
        case ConditionMode::Standardize: {
            const double sd = rStats.StdDev();
            const double mean = rStats.Mean();
            if (rStats.mNumFinite == 0 || !(sd > 0.0)) {
                if (!rWarnedConstant) {
                    rWarnedConstant = true;
                    log::warn("data_condition: '{}' has zero standard deviation; filling with 0",
                              rName);
                }
                t.mScale = 0.0;
                t.mOffset = 0.0;
                return t;
            }
            t.mScale = 1.0 / sd;
            t.mOffset = -mean / sd;
            return t;
        }
    }
    return t;
}

/// Applies one transform to a single value.
double dcond_apply(const DcondTransform& rT, double v) {
    if (!std::isfinite(v))
        return v;  // policy handled separately; never fold NaN into arithmetic
    if (rT.mClamp)
        return v < rT.mLo ? rT.mLo : (v > rT.mHi ? rT.mHi : v);
    return rT.mScale * v + rT.mOffset;
}

/// Reduces the row magnitudes of an array into one accumulator.
void dcond_accumulate_magnitude(const NDArray& rArray, std::size_t ncomp,
                                detail::FiniteStats& rStats) {
    const std::size_t total = rArray.Size();
    if (total == 0 || ncomp == 0)
        return;
    const std::size_t nrows = total / ncomp;
    const std::size_t grain = 4096;
    const std::size_t nchunks = (nrows + grain - 1) / grain;
    std::vector<detail::FiniteStats> partial(nchunks);
    parallel_for(
        nchunks,
        [&](std::size_t ci) {
            detail::FiniteStats local;
            const std::size_t begin = ci * grain;
            const std::size_t end = std::min(begin + grain, nrows);
            for (std::size_t r = begin; r < end; ++r) {
                double acc = 0.0;
                bool finite = true;
                for (std::size_t k = 0; k < ncomp; ++k) {
                    const double v = detail::read_double(rArray, r * ncomp + k);
                    if (!std::isfinite(v)) {
                        finite = false;
                        break;
                    }
                    acc += v * v;
                }
                local.Add(finite ? std::sqrt(acc) : std::nan(""));
            }
            partial[ci] = local;
        },
        1);
    for (const detail::FiniteStats& p : partial)
        rStats.Merge(p);
}

/// The dtype the output array should carry.
DType dcond_output_dtype(const NDArray& rSource, const DataConditionOptions& rOpts) {
    if (rOpts.mode == ConditionMode::Clamp && rOpts.preserve_dtype)
        return rSource.Dtype();
    return DType::Float64;
}

/// Conditions one array (component scope).
NDArray dcond_transform_component(const NDArray& rSource, std::size_t ncomp,
                                  const std::vector<DcondTransform>& rTransforms,
                                  const DataConditionOptions& rOpts, const std::string& rName) {
    NDArray dst = NDArray::Uninit(dcond_output_dtype(rSource, rOpts), rSource.Shape());
    const std::size_t total = rSource.Size();
    // Elementwise and bandwidth-bound.
    parallel_for_bw(total, [&](std::size_t i) {
        const double v = detail::read_double(rSource, i);
        detail::write_double(dst, i, dcond_apply(rTransforms[i % ncomp], v));
    });
    // The NaN policy runs serially so a `Fail` throw is deterministic.
    if (rOpts.nan_policy != NanPolicy::Ignore)
        for (std::size_t i = 0; i < total; ++i) {
            const double v = detail::read_double(dst, i);
            detail::write_double(dst, i, dcond_apply_nan_policy(v, rOpts, rName, i));
        }
    return dst;
}

/// Conditions one array (magnitude scope): each row is rescaled by the factor
/// that its magnitude would have received, so the direction is preserved.
NDArray dcond_transform_magnitude(const NDArray& rSource, std::size_t ncomp,
                                  const DcondTransform& rT, const DataConditionOptions& rOpts,
                                  const std::string& rName) {
    NDArray dst = NDArray::Uninit(dcond_output_dtype(rSource, rOpts), rSource.Shape());
    const std::size_t total = rSource.Size();
    const std::size_t nrows = ncomp > 0 ? total / ncomp : 0;
    parallel_for_bw(nrows, [&](std::size_t r) {
        double acc = 0.0;
        bool finite = true;
        for (std::size_t k = 0; k < ncomp; ++k) {
            const double v = detail::read_double(rSource, r * ncomp + k);
            if (!std::isfinite(v)) {
                finite = false;
                break;
            }
            acc += v * v;
        }
        const double mag = finite ? std::sqrt(acc) : std::nan("");
        // Scale the whole row so its magnitude becomes the conditioned value.
        // A zero-magnitude row has no direction to preserve and stays zero.
        const double target = dcond_apply(rT, mag);
        const double factor = (finite && mag > 0.0) ? target / mag : (finite ? 0.0 : mag);
        for (std::size_t k = 0; k < ncomp; ++k) {
            const double v = detail::read_double(rSource, r * ncomp + k);
            detail::write_double(dst, r * ncomp + k, finite ? v * factor : v);
        }
    });
    if (rOpts.nan_policy != NanPolicy::Ignore)
        for (std::size_t i = 0; i < total; ++i) {
            const double v = detail::read_double(dst, i);
            detail::write_double(dst, i, dcond_apply_nan_policy(v, rOpts, rName, i));
        }
    return dst;
}

/// Resolves the requested names, validating each.
std::vector<std::string> dcond_resolve_names(const Mesh& rMesh, const DataConditionOptions& rOpts) {
    if (rOpts.names.empty())
        return data_names(rMesh, rOpts.location);
    for (const std::string& n : rOpts.names)
        if (!data_has(rMesh, rOpts.location, n))
            throw std::invalid_argument(data_unknown_key_message(rMesh, rOpts.location, n));
    return rOpts.names;
}

}  // namespace

ConditionMode condition_mode_from_name(const std::string& rName) {
    if (rName == "clamp")
        return ConditionMode::Clamp;
    if (rName == "normalize" || rName == "rescale")
        return ConditionMode::Normalize;
    if (rName == "standardize" || rName == "zscore")
        return ConditionMode::Standardize;
    throw std::invalid_argument("meshio++: unknown conditioning mode '" + rName +
                                "' (expected 'clamp', 'normalize' or 'standardize')");
}

ConditionScope condition_scope_from_name(const std::string& rName) {
    if (rName == "component")
        return ConditionScope::Component;
    if (rName == "magnitude")
        return ConditionScope::Magnitude;
    throw std::invalid_argument("meshio++: unknown conditioning scope '" + rName +
                                "' (expected 'component' or 'magnitude')");
}

Mesh data_condition(const Mesh& rMesh, const DataConditionOptions& rOpts) {
    if (rOpts.mode == ConditionMode::Clamp && rOpts.lo > rOpts.hi)
        throw std::invalid_argument(
            "meshio++: data_condition: clamp bounds are inverted (lo=" + std::to_string(rOpts.lo) +
            " > hi=" + std::to_string(rOpts.hi) + ")");
    if (rOpts.mode == ConditionMode::Normalize && rOpts.lo > rOpts.hi)
        throw std::invalid_argument(
            "meshio++: data_condition: normalize target range is inverted (lo=" +
            std::to_string(rOpts.lo) + " > hi=" + std::to_string(rOpts.hi) + ")");

    const std::vector<std::string> names = dcond_resolve_names(rMesh, rOpts);
    Mesh out = detail::clone_mesh(rMesh);
    const bool magnitude = rOpts.scope == ConditionScope::Magnitude;

    for (const std::string& name : names) {
        const std::string target = name + rOpts.suffix;

        if (rOpts.location == DataLocation::Cell) {
            const std::size_t nblocks = rMesh.NumCellBlocks();
            if (rMesh.CellDataNumBlocks(name) != nblocks)
                throw std::invalid_argument(
                    "meshio++: data_condition: cell_data '" + name + "' has " +
                    std::to_string(rMesh.CellDataNumBlocks(name)) + " block(s) but the mesh has " +
                    std::to_string(nblocks) + " cell block(s)");
            if (nblocks == 0)
                continue;
            const std::size_t ncomp = data_num_components(rMesh.CellData(name, 0));

            // Statistics across ALL blocks first.
            std::vector<detail::FiniteStats> stats;
            detail::FiniteStats magstats;
            for (std::size_t b = 0; b < nblocks; ++b) {
                const NDArray& src = rMesh.CellData(name, b);
                if (data_num_components(src) != ncomp)
                    throw std::invalid_argument(
                        "meshio++: data_condition: cell_data '" + name +
                        "' has inconsistent component counts across blocks");
                if (magnitude)
                    dcond_accumulate_magnitude(src, ncomp, magstats);
                else
                    detail::accumulate_stats(src, ncomp, stats);
            }

            bool warned = false;
            std::vector<DcondTransform> transforms;
            DcondTransform magt;
            if (magnitude) {
                magt = dcond_build(magstats, rOpts, name, warned);
            } else {
                transforms.reserve(ncomp);
                for (std::size_t k = 0; k < ncomp; ++k)
                    transforms.push_back(dcond_build(stats[k], rOpts, name, warned));
            }

            std::vector<NDArray> blocks;
            blocks.reserve(nblocks);
            for (std::size_t b = 0; b < nblocks; ++b) {
                const NDArray& src = rMesh.CellData(name, b);
                blocks.push_back(
                    magnitude ? dcond_transform_magnitude(src, ncomp, magt, rOpts, name)
                              : dcond_transform_component(src, ncomp, transforms, rOpts, name));
            }
            out.AddCellData(target, std::move(blocks));
            continue;
        }

        const NDArray& src =
            rOpts.location == DataLocation::Point ? rMesh.PointData(name) : rMesh.FieldData(name);
        const std::size_t ncomp = data_num_components(src);
        bool warned = false;
        NDArray dst;
        if (magnitude) {
            detail::FiniteStats magstats;
            dcond_accumulate_magnitude(src, ncomp, magstats);
            const DcondTransform t = dcond_build(magstats, rOpts, name, warned);
            dst = dcond_transform_magnitude(src, ncomp, t, rOpts, name);
        } else {
            std::vector<detail::FiniteStats> stats;
            detail::accumulate_stats(src, ncomp, stats);
            std::vector<DcondTransform> transforms;
            transforms.reserve(ncomp);
            for (std::size_t k = 0; k < ncomp; ++k)
                transforms.push_back(dcond_build(stats[k], rOpts, name, warned));
            dst = dcond_transform_component(src, ncomp, transforms, rOpts, name);
        }
        if (rOpts.location == DataLocation::Point)
            out.AddPointData(target, std::move(dst));
        else
            out.AddFieldData(target, std::move(dst));
    }
    return out;
}

}  // namespace meshioplusplus

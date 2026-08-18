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
// The Zienkiewicz-Zhu (ZZ) recovery-based error indicator, plus a marking
// pass. Deliberately a COMPOSITION of shipped operations (gradient +
// cell_data_to_point_data + point_data_to_cell_data), not a new numerical
// kernel; see operations/error.hpp for the full contract.
//
// The three composed meshes (the raw gradient, the recovered points, the
// recovered-back-to-cells result) are private working state under fixed
// internal names -- never returned to the caller -- so a name collision with
// a real array in the input mesh cannot lose data: the actual result mesh is
// built separately from a fresh clone of the ORIGINAL input.

// System includes
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/operations/error.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/operations/data_average.hpp"
#include "meshioplusplus/operations/data_common.hpp"
#include "meshioplusplus/operations/gradient.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {

namespace {

constexpr const char* kErrPrefix = "meshio++: estimate_error: ";
// Private working names: never reach the returned mesh (see the header
// comment above), so any clash with a user array only affects internal
// intermediate copies, not the caller's data.
constexpr const char* kErrRawName = "__estimate_error_raw__";
constexpr const char* kErrRecSuffix = "_recovered";

/// One evaluable cell's indicator plus its global (block-major) index, for the
/// Fraction/Dorfler marking passes.
struct ErrEntry {
    double mValue;
    std::int64_t mGlobal;
};

/// Descending by value, ascending by global index on a tie -- a strict weak
/// order, so the result does not depend on how the per-cell work scheduled.
bool err_greater(const ErrEntry& rA, const ErrEntry& rB) {
    if (rA.mValue != rB.mValue)
        return rA.mValue > rB.mValue;
    return rA.mGlobal < rB.mGlobal;
}

}  // namespace

ErrorMethod error_method_from_name(const std::string& rName) {
    if (rName == "zz" || rName.empty())
        return ErrorMethod::Zz;
    throw std::invalid_argument(std::string(kErrPrefix) + "unknown method '" + rName +
                                "' (expected 'zz')");
}

ErrorMarking error_marking_from_name(const std::string& rName) {
    if (rName == "none" || rName.empty())
        return ErrorMarking::None;
    if (rName == "absolute" || rName == "abs")
        return ErrorMarking::Absolute;
    if (rName == "fraction" || rName == "frac")
        return ErrorMarking::Fraction;
    if (rName == "dorfler" || rName == "bulk")
        return ErrorMarking::Dorfler;
    throw std::invalid_argument(std::string(kErrPrefix) + "unknown marking policy '" + rName +
                                "' (expected 'none', 'absolute', 'fraction', or 'dorfler')");
}

ErrorResult estimate_error(const Mesh& rMesh, const ErrorOptions& rOptions) {
    // --- validation, all before any work -------------------------------------
    if (rOptions.mArrayName.empty())
        throw std::invalid_argument(std::string(kErrPrefix) + "an array name is required");
    if (!rMesh.HasPointData(rOptions.mArrayName)) {
        // Same restriction gradient states, for the same reason.
        if (rMesh.HasCellData(rOptions.mArrayName))
            throw std::invalid_argument(
                std::string(kErrPrefix) + "'" + rOptions.mArrayName +
                "' is cell_data, which is piecewise constant and has no derivative to recover; "
                "convert it first with cell_data_to_point_data (CLI: meshioplusplus data "
                "to-point)");
        throw std::invalid_argument(
            std::string(kErrPrefix) +
            data_unknown_key_message(rMesh, DataLocation::Point, rOptions.mArrayName));
    }
    if (rOptions.mMarking == ErrorMarking::Fraction || rOptions.mMarking == ErrorMarking::Dorfler) {
        if (!(rOptions.mMarkingValue > 0.0) || rOptions.mMarkingValue > 1.0)
            throw std::invalid_argument(
                std::string(kErrPrefix) +
                "mMarkingValue must be in (0, 1] for the 'fraction'/'dorfler' marking policy");
    }

    const std::string zz_name = rOptions.mOutputName.empty() ? kErrorZzName : rOptions.mOutputName;
    const std::string marked_name =
        rOptions.mMarkedName.empty() ? kErrorMarkedName : rOptions.mMarkedName;
    if (!rOptions.mOverwrite) {
        if (rMesh.HasCellData(zz_name))
            throw std::invalid_argument(std::string(kErrPrefix) + "'" + zz_name +
                                        "' already exists in cell_data (pass overwrite=true to "
                                        "replace it)");
        if (rOptions.mMarking != ErrorMarking::None && rMesh.HasCellData(marked_name))
            throw std::invalid_argument(std::string(kErrPrefix) + "'" + marked_name +
                                        "' already exists in cell_data (pass overwrite=true to "
                                        "replace it)");
    }

    // --- recovery: gradient -> cell_data_to_point_data -> point_data_to_cell_data ---
    // The recovery operator IS this existing averaging round trip -- reused,
    // not reimplemented, exactly as gradient(mLocation=Point) itself composes
    // cell_data_to_point_data. Every intermediate mesh below is private working
    // state (see the file comment); only `out`, built from a fresh clone of
    // rMesh further down, is ever returned.
    GradientOptions grad_opts;
    grad_opts.mArrayName = rOptions.mArrayName;
    grad_opts.mLocation = DataLocation::Cell;
    grad_opts.mOutputName = kErrRawName;
    grad_opts.mOverwrite = true;
    const GradientResult raw = gradient(rMesh, grad_opts);

    DataAverageOptions to_point;
    to_point.names = {kErrRawName};
    to_point.weight = CellPointWeight::Measure;
    to_point.overwrite = true;
    const Mesh recovered_pt = cell_data_to_point_data(raw.mMesh, to_point);

    DataAverageOptions to_cell;
    to_cell.names = {kErrRawName};
    to_cell.suffix = kErrRecSuffix;
    to_cell.overwrite = true;
    const Mesh recovered = point_data_to_cell_data(recovered_pt, to_cell);
    const std::string rec_name = std::string(kErrRawName) + kErrRecSuffix;

    // --- per-cell indicator ----------------------------------------------------
    const std::size_t nblocks = rMesh.NumCellBlocks();
    Mesh out = detail::clone_mesh(rMesh);
    if (nblocks == 0)
        return ErrorResult{std::move(out), 0.0, 0, 0};

    const std::vector<std::int64_t> bases = detail::block_bases(rMesh);
    const std::size_t total = static_cast<std::size_t>(detail::total_cells(bases));
    const NDArray& points = rMesh.Points();
    const std::size_t pdim = rMesh.PointDim();

    // A flat, block-major buffer of every cell's indicator (NaN where
    // unevaluable), used both as the per-block cell_data source and as the
    // input to the deterministic FiniteStats reduction below.
    NDArray global_zz(DType::Float64, {total});
    double* p_global = global_zz.As<double>();

    std::vector<NDArray> zz_blocks;
    zz_blocks.reserve(nblocks);
    std::int64_t skipped = 0;

    for (std::size_t b = 0; b < nblocks; ++b) {
        const auto cb = rMesh.Cells(b);
        const std::size_t ncells = cb.NumCells();
        const NDArray& raw_arr = recovered.CellData(kErrRawName, b);
        const NDArray& rec_arr = recovered.CellData(rec_name, b);
        const std::size_t comp = data_num_components(raw_arr);
        const std::size_t base = static_cast<std::size_t>(bases[b]);

        NDArray dst(DType::Float64, {ncells});
        double* pdst = dst.As<double>();
        // Per-cell flags reduced serially afterwards -- data_average's "apply
        // the policy serially" idiom, so mNumSkipped does not depend on
        // scheduling.
        std::vector<std::uint8_t> flag_skipped(ncells, 0);

        parallel_for(ncells, [&](std::size_t c) {
            bool finite_pair = true;
            double sumsq = 0.0;
            for (std::size_t k = 0; k < comp; ++k) {
                const double rv = detail::read_double(raw_arr, c * comp + k);
                const double gv = detail::read_double(rec_arr, c * comp + k);
                if (!std::isfinite(rv) || !std::isfinite(gv)) {
                    finite_pair = false;
                    break;
                }
                const double d = gv - rv;
                sumsq += d * d;
            }
            double zz = std::nan("");
            if (finite_pair) {
                // The same primitive data_average's CellPointWeight::Measure
                // just used for the recovery weighting, so the two agree by
                // construction.
                const double measure = std::fabs(detail::cell_measure(points, pdim, cb, c));
                if (std::isfinite(measure))
                    zz = std::sqrt(measure * sumsq);
            }
            if (!std::isfinite(zz))
                flag_skipped[c] = 1;
            pdst[c] = zz;
            p_global[base + c] = zz;
        });

        for (std::size_t c = 0; c < ncells; ++c)
            skipped += flag_skipped[c];
        zz_blocks.push_back(std::move(dst));
    }

    out.AddCellData(zz_name, std::move(zz_blocks));

    // --- global error, via the existing deterministic reduction ----------------
    std::vector<detail::FiniteStats> stats;
    detail::accumulate_stats(global_zz, 1, stats);
    const detail::FiniteStats combined = detail::combine_components(stats);
    const double global_error = combined.mNumFinite > 0 ? std::sqrt(combined.mSumSq) : 0.0;

    // --- marking -----------------------------------------------------------
    std::int64_t marked_count = 0;
    if (rOptions.mMarking != ErrorMarking::None) {
        std::vector<std::uint8_t> mark(total, 0);

        if (rOptions.mMarking == ErrorMarking::Absolute) {
            parallel_for_bw(total, [&](std::size_t i) {
                const double v = p_global[i];
                mark[i] = (std::isfinite(v) && v > rOptions.mMarkingValue) ? 1 : 0;
            });
        } else {
            // Fraction / Dorfler both need the full ranking, so build it once.
            std::vector<ErrEntry> entries;
            entries.reserve(static_cast<std::size_t>(combined.mNumFinite));
            for (std::size_t i = 0; i < total; ++i) {
                const double v = p_global[i];
                if (std::isfinite(v))
                    entries.push_back(ErrEntry{v, static_cast<std::int64_t>(i)});
            }
            std::sort(entries.begin(), entries.end(), err_greater);

            std::size_t take = 0;
            if (rOptions.mMarking == ErrorMarking::Fraction) {
                // int64(ratio * F + 0.5) truncation -- decimate's own
                // target-count rounding convention, never llround.
                take = static_cast<std::size_t>(
                    rOptions.mMarkingValue * static_cast<double>(entries.size()) + 0.5);
                take = std::min(take, entries.size());
            } else {
                // Dorfler: smallest prefix whose cumulative eta_K^2 (added
                // BEFORE the check below fires) already covers the target.
                const double target = rOptions.mMarkingValue * combined.mSumSq;
                double cum = 0.0;
                for (; take < entries.size(); ++take) {
                    if (cum >= target)
                        break;
                    cum += entries[take].mValue * entries[take].mValue;
                }
            }
            for (std::size_t i = 0; i < take; ++i)
                mark[static_cast<std::size_t>(entries[i].mGlobal)] = 1;
        }

        std::vector<NDArray> marked_blocks;
        marked_blocks.reserve(nblocks);
        for (std::size_t b = 0; b < nblocks; ++b) {
            const auto cb = rMesh.Cells(b);
            const std::size_t ncells = cb.NumCells();
            const std::size_t base = static_cast<std::size_t>(bases[b]);
            NDArray dst(DType::Int64, {ncells});
            std::int64_t* pdst = dst.As<std::int64_t>();
            parallel_for_bw(ncells, [&](std::size_t c) { pdst[c] = mark[base + c]; });
            marked_blocks.push_back(std::move(dst));
        }
        out.AddCellData(marked_name, std::move(marked_blocks));

        for (std::uint8_t m : mark)
            marked_count += m;
    }

    return ErrorResult{std::move(out), global_error, skipped, marked_count};
}

}  // namespace meshioplusplus

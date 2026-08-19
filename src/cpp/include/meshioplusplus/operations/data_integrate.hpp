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
 * @file operations/data_integrate.hpp
 * @brief Cell-measure-weighted reductions of a `cell_data` field: a physical
 * total (`sum(value * measure)`), the measure-weighted mean that total
 * implies, and both computed independently for every named `Cell` region --
 * `gradient`'s natural companion (`gradient` differentiates a field; this
 * integrates one).
 *
 * **Why it lives in the mesh-operations layer rather than the `data_*`
 * bundle.** The five `data_*` operations are *defined* by never touching
 * geometry (`operations/data_common.hpp`, `doc/data_operations.md`); this one
 * reads `detail::cell_measure` (face/volume tables), so it belongs here, with
 * the `mPascalCase` option fields the geometry operations use -- the exact
 * reasoning `gradient.hpp`'s own file doc comment states for itself. It is
 * still *reachable* as `meshioplusplus data integrate` in both CLIs, because
 * that is where a user looks for it.
 *
 * ### Contract
 *
 * - **`cell_data` only.** A `point_data` name throws by name, pointing at
 *   `point_data_to_cell_data` (CLI `data to-cell`) -- the mirror image of
 *   `gradient`'s own `point_data`-only contract, which throws on `cell_data`
 *   pointing at `data to-point`. Proper `point_data` integration needs
 *   shape-function quadrature, a separate and larger job.
 * - **Weighting.** Every sum is weighted by `|detail::cell_measure(...)|`
 *   (length / area / volume by the cell's own topological dimension). A cell
 *   whose measure is not computable (a ragged or unsupported-type block, or a
 *   degenerate one) is excluded from **both** the numerator and the
 *   denominator -- never given a fallback weight of 1, unlike
 *   `cell_data_to_point_data`'s own `Measure` weighting, because a silent
 *   unit-weight substitution would corrupt a physical total in a way it only
 *   softens an average. Excluded cells are counted in
 *   `FieldIntegralRegion::mNumSkipped`.
 * - **Non-finite values** exclude a cell from that *component's* numerator
 *   **and** denominator (not just zeroed into the numerator, or the mean
 *   would be silently biased), counted per component in
 *   `FieldIntegralRegion::mNumNanPerComponent`. This mirrors the geometric
 *   exclusion rule applied one level down, component by component.
 * - **Per component only** -- there is no whole-array (cross-component)
 *   total or mean, unlike `DataArrayInfo`'s collapsed statistics. Summing a
 *   vector field's components into one number has no general physical
 *   meaning, unlike a flattened numeric mean, which needs only "these are
 *   numbers".
 * - **No `NanPolicy`.** Like `gradient` and `data_info`, this is a reduction
 *   with nothing to exclude a value *from*; non-finite values are always
 *   excluded and always counted, with no configurable policy.
 * - **Regions are not filtered, and not a partition.** Every named `Cell`
 *   region present on the mesh gets its own independent entry (mirroring
 *   `split`'s own "split by region" contract) -- a cell belonging to two
 *   regions contributes fully to both; a cell in none contributes to
 *   neither. `Point`/`Side` regions are skipped, exactly as `split` skips
 *   them for the same reason (no cell-measure meaning for a point or facet
 *   group).
 * - The mesh is **never modified**; this returns a read-only report, like
 *   `data_info` and `compute_stats`.
 *
 * Output is byte-identical across the three mesh backends, across thread
 * counts, and across the C++/numpy boundary: per-cell measures are computed
 * once (shared across every array and every region) via `parallel_for`, and
 * every weighted-sum reduction is chunked-`parallel_for`-then-serially-merged
 * (`detail::accumulate_weighted`), the same idiom `data_info`'s
 * `detail::accumulate_stats` already uses.
 *
 * Everything is standard C++ and the uniform mesh API only, so it compiles
 * under every mesh backend. This is an operation, not a file format -- it is
 * not in the format registry.
 */

// System includes
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// Cell-measure-weighted reduction of one array over one set of cells (the
/// whole mesh, or one named `Cell` region) -- one entry per component.
struct FieldIntegralRegion {
    /// The region's name; empty for the whole-mesh entry.
    std::string mName;
    /// Cells with a computable measure (candidates for every component).
    std::int64_t mNumCells = 0;
    /// Cells excluded because their measure was not computable (ragged,
    /// unsupported type, or degenerate) -- the same count for every
    /// component, since it is purely geometric.
    std::int64_t mNumSkipped = 0;
    /// `sum(|measure|)` over cells finite in component k, size `mNumComponents`.
    std::vector<double> mDomainMeasurePerComponent;
    /// `sum(value * |measure|)` over the same cells, size `mNumComponents`.
    std::vector<double> mTotalPerComponent;
    /// `mTotalPerComponent[k] / mDomainMeasurePerComponent[k]`, or NaN when
    /// that denominator is zero.
    std::vector<double> mMeanPerComponent;
    /// Measurable cells excluded from component k because its value was
    /// non-finite there.
    std::vector<std::int64_t> mNumNanPerComponent;
};

/// One array's whole-mesh integral, plus one entry per named `Cell` region.
struct FieldIntegralArray {
    std::string mName;
    std::int64_t mNumComponents = 1;
    /// The whole-mesh reduction (`mName` empty).
    FieldIntegralRegion mDomain;
    /// One independent entry per named `Cell` region present on the mesh,
    /// in `Mesh::RegionNames()`'s sorted order.
    std::vector<FieldIntegralRegion> mRegions;
};

/// The full report (see `data_integrate`).
struct DataIntegrateReport {
    std::vector<FieldIntegralArray> mArrays;
};

/// Options for `data_integrate`.
struct DataIntegrateOptions {
    /// `cell_data` array names to integrate; empty means every `cell_data`
    /// array, in sorted name order.
    std::vector<std::string> mArrayNames;
};

/**
 * @brief Cell-measure-weighted total/mean of one or more `cell_data` arrays,
 * for the whole mesh and independently for every named `Cell` region.
 * @param rMesh the mesh to integrate over (never modified).
 * @param rOptions which arrays to integrate.
 * @return the populated report.
 * @throws std::invalid_argument on an unknown array name (listing what
 *         exists), a `point_data`-only name (naming `point_data_to_cell_data`
 *         as the fix), or a `cell_data` array whose block count disagrees
 *         with the mesh.
 */
MESHIOPLUSPLUS_API DataIntegrateReport data_integrate(const Mesh& rMesh,
                                                      const DataIntegrateOptions& rOptions = {});

}  // namespace meshioplusplus

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
 * @file vtkhdf_time_series.hpp
 * @brief Transient VTKHDF writing: the geometry once, then one step at a time.
 *
 * The C++ twin of Python's `meshioplusplus.vtkhdf.TimeSeriesWriter` and the
 * VTKHDF counterpart of `XdmfTimeSeriesWriter`, deliberately shaped the same way
 * (`WritePointsCells`, `WriteData`, `Flush`, `Finalize`, move-only) so the
 * sequence engine drives either through one code path. It is a stateful
 * multi-call object rather than a `(path, mesh)` format, and therefore
 * deliberately **not** in `registry.cpp`: there is no single call for the
 * registry to name.
 *
 * ### On-disk shape
 *
 * Geometry x 1, fields x N -- exactly what VTKHDF's `Steps` offset tables are for.
 * `NumberOfPoints`/`Points`/`Connectivity`/`Offsets`/`Types` are written once;
 * every `PointData/<name>`, `CellData/<name>` and `FieldData/<name>` is a chunked,
 * extendable dataset that gains one step's rows per `WriteData`, and every
 * `Steps/*` table gains one entry. `PartOffsets`, `PointOffsets`, `CellOffsets`
 * and `ConnectivityIdOffsets` are all zero (every step points at the same
 * geometry); `PointDataOffsets/<name>` counts `k * NumberOfPoints`.
 *
 * ### Differences from `XdmfTimeSeriesWriter`, which are improvements
 *
 * - **Nothing is buffered.** XDMF's light data is written once, by `Finalize()`;
 *   here every step lands in the HDF5 file as it is written and `NSteps` is
 *   rewritten after each, so a run that is killed leaves a file ParaView opens,
 *   covering every completed step.
 * - So `Flush()` is `H5Fflush`, cheap and **not** quadratic in the step count, and
 *   auto-flush defaults to **on** (XDMF's defaults to off for exactly that reason).
 *
 * ### The fixed name set
 *
 * The names of the point/cell/field arrays are fixed at the first `WriteData`. A
 * later step introducing a new name, dropping one, or changing an array's shape
 * is a `WriteError` naming it: VTK indexes `PointDataOffsets/<name>` by step, so a
 * short table would be a silent truncation.
 */

#ifdef MESHIOPLUSPLUS_HAS_HDF5

// System includes
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Whether a series starts fresh or continues an existing one.
 *
 * Appending to a path that does **not** exist yet is not an error -- it is exactly
 * `Truncate`, so a solver can pass `Append` unconditionally.
 */
enum class VtkhdfSeriesMode { Truncate, Append };

/**
 * @brief Writes a transient VTKHDF `UnstructuredGrid`: one static grid, then one
 * step at a time.
 *
 * @code
 * meshioplusplus::VtkhdfTimeSeriesWriter w("out.vtkhdf");
 * w.WritePointsCells(mesh);
 * for (int k = 0; k < nsteps; ++k) {
 *     solve(mesh);
 *     w.WriteData(k * dt, mesh);   // point_data, cell_data and field_data are consumed
 * }
 * w.Finalize();                     // the destructor would do this too
 * @endcode
 *
 * Move-only. A moved-from writer answers its observers as for a finished series
 * and throws `WriteError` from `WritePointsCells` and both `WriteData` overloads,
 * because a caller writing a step must never be left believing it landed.
 */
class MESHIOPLUSPLUS_API VtkhdfTimeSeriesWriter {
public:
    /** @brief One solver array for the `NamedArray` `WriteData` overload; row-major. */
    struct NamedArray {
        std::string mName;
        std::size_t mNumComponents = 1;
        std::vector<double> mValues;
    };

    /**
     * @brief Open a series for writing.
     *
     * @param rPath Path of the `.vtkhdf` file.
     * @param GzipLevel Deflate level 0-9 for the datasets, or negative (the default)
     *        for none. Appended datasets are always chunked; the level only adds the
     *        filter.
     * @param Mode `Truncate` starts a fresh series; `Append` continues the one at
     *        `rPath` if there is one.
     * @throws WriteError if `Append` was asked for and the existing file is not a
     *         transient VTKHDF `UnstructuredGrid` this writer can continue.
     */
    explicit VtkhdfTimeSeriesWriter(const std::string& rPath, int GzipLevel = -1,
                                    VtkhdfSeriesMode Mode = VtkhdfSeriesMode::Truncate);

    /**
     * @brief Finalizes the series if it has not been finalized already.
     *
     * @warning This closes the file. Delete the output only after the writer is gone
     *          (or after `Finalize()`), never while it is alive.
     */
    ~VtkhdfTimeSeriesWriter();

    VtkhdfTimeSeriesWriter(const VtkhdfTimeSeriesWriter&) = delete;
    VtkhdfTimeSeriesWriter& operator=(const VtkhdfTimeSeriesWriter&) = delete;
    VtkhdfTimeSeriesWriter(VtkhdfTimeSeriesWriter&&) noexcept;
    VtkhdfTimeSeriesWriter& operator=(VtkhdfTimeSeriesWriter&&) noexcept;

    /**
     * @brief Write the static grid: the points and cells every step shares.
     *
     * Only geometry and connectivity are consumed. Call exactly once, before the
     * first `WriteData`. When continuing an existing series (`Append`) the grid is
     * already on disk: this call then only checks that `rMesh` has the same point
     * and cell counts, and writes nothing, so a driver can call it unconditionally.
     * @throws WriteError if called twice or after `Finalize()`, if the points exceed
     *         dimension 3, if a cell type has no VTK id, or if the counts of an
     *         appended series differ.
     */
    void WritePointsCells(const Mesh& rMesh);

    /**
     * @brief Write one step's `point_data`, `cell_data` and `field_data`.
     *
     * `cell_data` is concatenated across blocks (so `rMesh` must have the block
     * structure passed to `WritePointsCells`); the reserved `meshio:time` field is
     * skipped, since `Time` already records it. Arrays are written in sorted-name
     * order and keep their dtype.
     * @param Time The step's simulation time, appended to `Steps/Values`.
     * @throws WriteError if `WritePointsCells` has not run, the series is finalized,
     *         or the arrays differ from the first step's (see the file comment).
     */
    void WriteData(double Time, const Mesh& rMesh);

    /**
     * @brief Write one step from raw solver arrays, with no `Mesh` in between.
     *
     * @param rPointData Nodal arrays; each must hold `NumPoints * mNumComponents` values.
     * @param rCellData Cell arrays; each must hold `NumCells * mNumComponents` values,
     *        `NumCells` being the total across every block.
     * @throws WriteError as for the `Mesh` overload, or if an array's length does
     *         not match -- stricter than the `Mesh` overload, which reads the counts
     *         off the mesh.
     */
    void WriteData(double Time, const std::vector<NamedArray>& rPointData,
                   const std::vector<NamedArray>& rCellData = {});

    /** @brief `H5Fflush`: everything written so far is durable. Cheap; a no-op once finalized. */
    void Flush();

    /** @brief Flush automatically after every `WriteData` (default: **on**, see above). */
    void SetAutoFlush(bool Enable);
    /** @brief Whether auto-flush is enabled. */
    bool AutoFlush() const;

    /** @brief Flush and close the file. Idempotent; the destructor calls it. */
    void Finalize();

    /** @brief How many steps have been written so far. */
    std::size_t NumSteps() const;

    /** @brief Whether `Finalize()` has already run. */
    bool Finalized() const;

private:
    struct Impl;
    std::unique_ptr<Impl> mImpl;
};

}  // namespace meshioplusplus

#endif  // MESHIOPLUSPLUS_HAS_HDF5

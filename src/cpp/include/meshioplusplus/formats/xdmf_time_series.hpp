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
 * @file xdmf_time_series.hpp
 * @brief Transient (time-series) XDMF3 writing — one static grid plus one
 *        `<Grid>` per time step inside a temporal collection.
 *
 * The C++ twin of Python's `meshioplusplus.xdmf.TimeSeriesWriter`, and the write
 * half of what `ReadOptions::mTimeStep` / `MeshMetadata::mTimeValues` already
 * expose on the read side. A simulation writes the mesh **once** and then one
 * cheap step per solve, which is what makes this a stateful multi-call object
 * rather than a `(path, mesh)` format — and therefore why it is deliberately
 * **not** in `registry.cpp`: there is no single call for the registry to name.
 *
 * The emitted document is
 * @code{.xml}
 * <Xdmf Version="3.0" xmlns:xi="http://www.w3.org/2003/XInclude">
 *   <Domain>
 *     <Grid Name="TimeSeries_meshio" GridType="Collection" CollectionType="Temporal">
 *       <Grid>
 *         <xi:include xpointer="xpointer(...)"/>
 *         <Time Value="0"/>
 *         <Attribute Name="u" Center="Node"><DataItem .../></Attribute>
 *       </Grid>
 *       ...
 *     </Grid>
 *     <Grid Name="mesh" GridType="Uniform">
 *       <Geometry GeometryType="XYZ"><DataItem .../></Geometry>
 *       <Topology TopologyType="Tetrahedron" ...><DataItem .../></Topology>
 *     </Grid>
 *   </Domain>
 * </Xdmf>
 * @endcode
 * i.e. structurally what the Python writer emits, so the same files are readable
 * by ParaView, by Python's `TimeSeriesReader`, and by this repo's own
 * `read_xdmf` (which resolves the collection structurally rather than running an
 * XInclude/XPointer pass).
 *
 * Heavy data goes through the same `xdmfcommon::DataItemStore` the single-step
 * `write_xdmf` uses, so dataset naming (`data0`, `data1`, ... in one flat
 * namespace, spanning every step), the sibling `.h5`/`.bin` file naming and the
 * inline-XML number formatting cannot drift between the two writers. One
 * consequence worth stating: for `"HDF"` the companion file is a **sibling of
 * the `.xdmf`**, not a file in the process's working directory — the Python
 * writer's `filename.stem + ".h5"` resolves against the CWD, which breaks the
 * moment the output path has a directory component.
 *
 * The XML light data is buffered in memory and written **once**, by `Finalize()`
 * (or the destructor). A series is therefore only readable after the writer is
 * done with it; that is inherent to the format, since the collection element has
 * to enclose every step.
 */

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
 * @brief Writes a transient XDMF3 series: one static grid, then one step at a
 * time.
 *
 * Usage is write-the-mesh-then-write-the-steps:
 * @code
 * meshioplusplus::XdmfTimeSeriesWriter w("out.xdmf");
 * w.WritePointsCells(mesh);
 * for (int k = 0; k < nsteps; ++k) {
 *     solve(mesh);
 *     w.WriteData(k * dt, mesh);   // only point_data/cell_data are consumed
 * }
 * w.Finalize();                     // the destructor would do this too
 * @endcode
 *
 * Move-only (it owns an open heavy-data container and an unwritten XML
 * document); copying it would write the same file twice.
 */
/**
 * @brief Whether a series starts fresh or continues an existing one.
 *
 * `Append` is what a restartable analysis needs: a run that was killed, hit a
 * node failure, or simply finished a stage leaves an `.xdmf` plus its heavy
 * data, and the next run must continue that temporal collection rather than
 * overwrite it.
 *
 * Appending to a path that does **not** exist yet is not an error -- it is
 * exactly `Truncate`. That is deliberate: a solver can then pass `Append`
 * unconditionally instead of having to probe the filesystem and branch.
 */
enum class XdmfSeriesMode { Truncate, Append };

class MESHIOPLUSPLUS_API XdmfTimeSeriesWriter {
public:
    /**
     * @brief One solver array for the `NamedArray` `WriteData` overload.
     *
     * Values are row-major, `mNumComponents` per entity: `{n}` for a scalar
     * field, `{n, 3}` for a vector, and so on.
     */
    struct NamedArray {
        std::string mName;
        std::size_t mNumComponents = 1;
        std::vector<double> mValues;
    };

    /**
     * @brief Open a series for writing.
     *
     * Nothing is written yet: the `.xdmf` appears at `Finalize()` (or at the
     * first `Flush()`), and the `"HDF"` companion at the first stored array.
     *
     * @param rPath Path of the `.xdmf`/`.xmf` light-data file to write.
     * @param rDataFormat `"HDF"` (default; companion `<base>.h5`, needs an
     *        HDF5-enabled build), `"XML"` (numbers inline) or `"Binary"`
     *        (sibling `<base><n>.bin` raw files).
     * @param GzipLevel gzip level for `"HDF"` datasets; negative (the default)
     *        means uncompressed. Ignored for the other formats.
     * @param Mode `Truncate` (the default) starts a fresh series; `Append`
     *        continues the one already at `rPath`, if any (see #XdmfSeriesMode).
     * @throws WriteError if `rDataFormat` is unrecognized, is `"HDF"` on a
     *         build without HDF5 support, or if `Append` was asked for and the
     *         existing file is not a series this writer can continue.
     */
    explicit XdmfTimeSeriesWriter(const std::string& rPath, const std::string& rDataFormat = "HDF",
                                  int GzipLevel = -1,
                                  XdmfSeriesMode Mode = XdmfSeriesMode::Truncate);

    /**
     * @brief Finalizes the series if it has not been finalized already.
     *
     * @warning Because this **writes** `rPath`, deleting the output while the
     *          writer is still alive recreates it:
     *          @code
     *          {
     *              XdmfTimeSeriesWriter w(path, "XML");
     *              w.WritePointsCells(m);
     *              w.WriteData(0.0, {u});
     *              std::filesystem::remove(path);   // "clean up"
     *          }   // <-- the destructor finalizes here, recreating path
     *          @endcode
     *          The second half is what makes this bite: with
     *          `XdmfSeriesMode::Append` the *next* run then continues a series
     *          you believed deleted, so the failure surfaces one run later as a
     *          wrong step count rather than at the delete. Destroy the writer
     *          before removing its output -- end the scope, or call `Finalize()`
     *          explicitly and then delete.
     */
    ~XdmfTimeSeriesWriter();

    XdmfTimeSeriesWriter(const XdmfTimeSeriesWriter&) = delete;
    XdmfTimeSeriesWriter& operator=(const XdmfTimeSeriesWriter&) = delete;
    // Moving leaves the source empty. Its contract, so a moved-from writer is
    // diagnosable rather than undefined: the observers and the idempotent
    // operations (`SetAutoFlush`, `AutoFlush`, `NumSteps`, `Finalized`,
    // `Flush`, `Finalize`) are safe no-ops answering as for a finished series,
    // while the three that cannot silently do nothing -- `WritePointsCells` and
    // both `WriteData` overloads -- throw `WriteError`, because a caller
    // writing a step must never be left believing it landed.
    XdmfTimeSeriesWriter(XdmfTimeSeriesWriter&&) noexcept;
    XdmfTimeSeriesWriter& operator=(XdmfTimeSeriesWriter&&) noexcept;

    /**
     * @brief Write the static grid: the points and cell blocks every step shares.
     *
     * Only geometry and connectivity are consumed; any data the mesh carries is
     * ignored here, because in a transient series data belongs to a step. Call
     * exactly once, before the first `WriteData`.
     *
     * @param rMesh The mesh whose points/cells define the series.
     * @throws WriteError if called twice, if the points exceed dimension 3, or
     *         if a cell type has no XDMF spelling.
     */
    void WritePointsCells(const Mesh& rMesh);

    /**
     * @brief Write one time step's `point_data` and `cell_data`.
     *
     * Arrays are emitted in the uniform API's sorted-name order, exactly as
     * `write_xdmf` does, so a series is deterministic. `cell_data` is
     * concatenated across blocks into XDMF's one-array-per-name raw layout
     * (`raw_from_cell_data`), which means the mesh passed here must have the
     * same block structure as the one passed to `WritePointsCells`.
     *
     * @param Time The step's simulation time, emitted as `<Time Value=...>`.
     * @param rMesh The mesh carrying this step's data; its geometry is ignored.
     * @throws WriteError if `WritePointsCells` has not been called, or if the
     *         series is already finalized.
     */
    void WriteData(double Time, const Mesh& rMesh);

    /**
     * @brief Write one step from raw solver arrays, with no `Mesh` in between.
     *
     * The granularity a solver actually has: after `WritePointsCells` the
     * geometry and connectivity are fixed, and each step produces nodal and
     * elemental values, not a new mesh. Under the KRATOS mesh backend the
     * difference is not cosmetic -- handing over a `Mesh` per step re-stages the
     * whole ModelPart every output step when only the values changed.
     *
     * Arrays are emitted in the order given: unlike the `Mesh` overload, whose
     * sorted-name order comes from the uniform API, here the caller's order *is*
     * the deterministic order.
     *
     * @param Time The step's simulation time.
     * @param rPointData Nodal arrays; each must hold `NumPoints * mNumComponents`
     *        values.
     * @param rCellData Cell arrays; each must hold `NumCells * mNumComponents`
     *        values, `NumCells` being the total across every block.
     * @throws WriteError if `WritePointsCells` has not been called, if the series
     *         is already finalized, or if an array's length does not match --
     *         deliberately stricter than the `Mesh` overload, which cannot get
     *         this wrong because it reads the counts off the mesh.
     */
    void WriteData(double Time, const std::vector<NamedArray>& rPointData,
                   const std::vector<NamedArray>& rCellData = {});

    /**
     * @brief Write the light data as it currently stands, without finalizing.
     *
     * The point of the whole feature: without it the `.xdmf` appears only at
     * `Finalize()`, so a run that is killed, crashes, or is simply still going
     * leaves heavy data on disk and **nothing readable**. After a `Flush()` the
     * file opens in ParaView and covers every step written so far.
     *
     * The document is written to a sibling temporary file and renamed over the
     * target, so a crash *during* a flush cannot truncate the previous one; the
     * heavy data is flushed first, so the `.xdmf` never points at a dataset that
     * is not on disk yet.
     *
     * Safe to call any number of times, and a no-op once finalized.
     *
     * @throws WriteError if the `.xdmf` cannot be written.
     */
    void Flush();

    /**
     * @brief Flush automatically after every `WriteData` (default: off).
     *
     * Off by default because a flush re-serializes the whole document, making a
     * per-step flush quadratic in the step count -- and for `"XML"`, whose heavy
     * data lives *in* that document, quadratic in the data volume too. Turn it
     * on for a long solve where crash-resilience is worth more than the writing
     * cost, or call `Flush()` on your own cadence (every n-th step, say).
     */
    void SetAutoFlush(bool Enable);
    /** @brief Whether auto-flush is enabled. */
    bool AutoFlush() const;

    /**
     * @brief Write the `.xdmf` light data and close the heavy-data container.
     *
     * Idempotent, and called by the destructor — an explicit call exists so a
     * failure surfaces as an exception rather than being swallowed during stack
     * unwinding.
     *
     * @throws WriteError if the `.xdmf` cannot be written.
     */
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

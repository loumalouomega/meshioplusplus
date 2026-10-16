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
 * @file vtx.hpp
 * @brief DOLFINx VTX output (an ADIOS2 `.bp` directory) reader, the `vtx`
 *        format.
 *
 * `dolfinx.io.VTXWriter` writes an ADIOS2 BP4/BP5 directory holding a
 * `vtk.xml` attribute (the schema ParaView's `ADIOS2VTXReader` follows) and,
 * per step, the variables it names: `geometry` (points, n x 3), `connectivity`
 * (cells x (1 + nodes), each row prefixed with its node count, VTK legacy
 * style), `types` (the VTK cell type), the point and cell data arrays, and the
 * `step` time. Every MPI rank writes one block with its own point numbering and
 * ghost points (`vtkGhostType`, `vtkOriginalPointIds`).
 *
 * ### Reading
 *
 * - Every step is one mesh; `ReadOptions::mTimeStep` selects it and its time is
 *   attached as `field_data["meshio:time"]`. A step without mesh variables (a
 *   `VTXMeshPolicy.reuse` file writes the mesh in step 0 only) takes the mesh of
 *   the last step that has one.
 * - The rank blocks are concatenated in block order, as ParaView merges them.
 *   `ReadOptions::mGhosts == GhostPolicy::Drop` instead welds every ghost point
 *   onto its owner through `vtkOriginalPointIds`, which gives the serial mesh,
 *   and removes the two ghost arrays.
 * - `VTK_LAGRANGE_*` cells whose node count is a linear or quadratic VTK cell's
 *   with the same node order (every degree-1 cell, and the 3-node curve,
 *   6-node triangle, 9-node quadrilateral and 10-node tetrahedron) are read as
 *   that cell; other degrees stay `VTK_LAGRANGE_*`.
 * - A point or cell array block whose length matches no rank's point or cell
 *   count is skipped (DOLFINx 0.11 writes a stray one-value block of the ghost
 *   arrays in the first step of a `reuse` file).
 * - A `.bp` without a `vtk.xml` attribute (an `adios4dolfinx` checkpoint, Fides
 *   output) is a `ReadError`.
 *
 * Compiled only with `MESHIOPLUSPLUS_WITH_ADIOS2=ON`. See doc/formats/vtx.md.
 */

#ifdef MESHIOPLUSPLUS_HAS_ADIOS2

// System includes
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/read_options.hpp"

namespace meshioplusplus {

/**
 * @brief Read one step of a DOLFINx VTX `.bp` directory.
 *
 * Honoured: `mTimeStep`, `mPointsOnly`, `mDataArrays` and `mGhosts`.
 * @throws ReadError if the path is not an ADIOS2 file with a `vtk.xml` schema,
 *         the schema is not an `UnstructuredGrid`, no step up to the requested
 *         one carries a mesh, or `mTimeStep` is out of range.
 */
MESHIOPLUSPLUS_API Mesh read_vtx(const std::string& rPath, const ReadOptions& rOpts = {});

/** @brief The `step` time of every step, or the step indices where it is absent. */
MESHIOPLUSPLUS_API std::vector<double> vtx_time_values(const std::string& rPath);

/**
 * @brief Summarize a VTX file from its schema and block sizes, without reading
 *        points or data: counts and cell blocks of step 0, array names and
 *        `mTimeValues`.
 */
MESHIOPLUSPLUS_API MeshMetadata read_vtx_metadata(const std::string& rPath,
                                                  const ReadOptions& rOpts = {});

/**
 * @brief Write a DOLFINx VTX `.bp` directory (BP5) holding one step: `VTXWriter`'s
 *        layout for one rank (a `vtk.xml` schema attribute, `step`, `geometry`,
 *        `connectivity`, `types`, `NumberOfNodes`, `NumberOfEntities` and one
 *        variable per point or cell array). The time is `meshio:time` (0 when
 *        absent). Cells are written as the `VTK_LAGRANGE_*` types DOLFINx uses
 *        where the node order allows. The Python twin writes the same variables.
 * @throws WriteError when the mesh has several cell types (as DOLFINx's files
 *         hold one), no cells, or the file cannot be written.
 */
MESHIOPLUSPLUS_API void write_vtx(const std::string& rPath, const Mesh& rMesh);

/**
 * @brief A VTX file with a step per mesh, the mesh written at every step
 *        (`VTXMeshPolicy.update`). Every step must have the first step's points,
 *        cells and array names. The Python twin is `meshioplusplus.vtx.SeriesWriter`.
 */
class MESHIOPLUSPLUS_API VtxSeriesWriter {
public:
    explicit VtxSeriesWriter(const std::string& rPath);
    ~VtxSeriesWriter();

    VtxSeriesWriter(const VtxSeriesWriter&) = delete;
    VtxSeriesWriter& operator=(const VtxSeriesWriter&) = delete;
    VtxSeriesWriter(VtxSeriesWriter&&) noexcept;
    VtxSeriesWriter& operator=(VtxSeriesWriter&&) noexcept;

    /// Add one step. @throws WriteError on a mismatch with the first step or a failure.
    void Write(double Time, const Mesh& rMesh);

    /// The number of steps written so far.
    std::size_t NumSteps() const noexcept;

    /// Close the file. Idempotent.
    void Finalize();

private:
    struct Impl;
    std::unique_ptr<Impl> mpImpl;
};

}  // namespace meshioplusplus

#endif  // MESHIOPLUSPLUS_HAS_ADIOS2

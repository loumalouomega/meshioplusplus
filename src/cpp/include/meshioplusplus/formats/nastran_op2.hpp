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
 * @file nastran_op2.hpp
 * @brief Nastran OP2 result file reader (MSC and NX/Simcenter, 32- and 64-bit).
 *
 * An OP2 is a Fortran unformatted file (`detail/fortran_records.hpp`; 4-byte
 * record markers, 4- or 8-byte words) of named tables: a name, a header record,
 * then records separated by marker triples `[-k, 1, 0]` and closed by `0`.
 *
 * The mesh comes from the geometry tables, as `nastran_h5` builds it
 * (`detail/nastran_model.hpp`): `GEOM1` GRIDs are the points (moved to basic
 * through the GEOM1 CORD records when `CP != 0`; `nastran:cp`/`nastran:cd`),
 * `GEOM2` element records the cells (`nastran:eid`, `nastran:pid`; CONM2 as
 * `vertex`; records with no cell type are named in a warning) and `EPT` the
 * property card of each region `<PTYPE>_<pid>`. Without GRID records the input
 * deck beside the file (`<stem>.bdf`, `.dat`, `.nas` or `.blk`) gives the
 * mesh, read by the bulk-data reader, with `nastran:eid`.
 *
 * Every (subcase, analysis, mode/time/frequency/load step) of the supported
 * tables is a step: `time_step` selects one (0 = first, negative from the
 * end); `field_data["meshio:time"]` is its frequency or time, the eigenvalue of
 * a mode, 0 for a static subcase, with `nastran:subcase`, `nastran:analysis`
 * and `nastran:mode`. Read:
 *  - `OUG*`/`BOUG*`, `OQG*`, `OQMG*`, `OPG*` real SORT1 tables as point data
 *    `DISPLACEMENT`, `EIGENVECTOR`, `VELOCITY`, `ACCELERATION`, `SPC_FORCE`,
 *    `MPC_FORCE`, `APPLIED_LOAD` (each with a `_ROT` twin) and `TEMPERATURE`,
 *    NaN for points without a value, rotated from each GRID's `CD` to basic
 *    (not `BOUG*`, which is basic already);
 *  - `OES*`/`OSTR*` real SORT1 stress and strain of rods (1, 3, 10), shear
 *    panels (4), bars (34), shells (33, 74, and the centre of 64, 70, 75, 82,
 *    144) and solids (39, 67, 68, 255) as cell data `STRESS:<M>`/`STRAIN:<M>`
 *    named like `nastran_h5`'s members (`X1`, `TXY1`, `X`, `TZX`, `A`...) plus
 *    the derived values (`VON_MISES1`, `MAJOR1`, `PRINCIPAL_A`...), the centre
 *    value NaN on cells without it; corner, ply (95-98, 232, 233) and station
 *    (CBEAM 2, CBAR 100) values as `<name>@corner|@ply|@station` with
 *    `nastran:layout:<name>`, as `nastran_h5`;
 *  - `OGPFB*` grid point forces as `GRID_FORCE:<M>` per element node and
 *    `GRID_FORCE:<label>:<M>` point data.
 * Complex, random and SORT2 tables, other element types and other tables are
 * skipped with a warning naming them.
 *
 * The writer (`write_nastran_op2`, `NastranOp2SeriesWriter`) is the reader's
 * inverse for the core of the format: 32-bit little-endian words, the MSC
 * dialect, SORT1 real tables. See `write_nastran_op2`.
 *
 * See doc/formats/nastran_op2.md.
 */

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
 * @brief Read one step of a Nastran OP2 file.
 * @param rPath filesystem path to read
 * @param rOpts `mTimeStep` selects the step; `mArrays`/`mPointsOnly` narrow
 *        the data
 * @throws ReadError if the file is not an OP2, is truncated, has no geometry
 *         and no sibling deck, or `mTimeStep` is out of range
 */
MESHIOPLUSPLUS_API Mesh read_nastran_op2(const std::string& rPath, const ReadOptions& rOpts = {});

/**
 * @brief The time (frequency, time, eigenvalue, or 0) of every step, in file order.
 */
MESHIOPLUSPLUS_API std::vector<double> nastran_op2_time_values(const std::string& rPath);

/**
 * @brief Metadata (counts, data names, `mTimeValues`) of an OP2 file.
 */
MESHIOPLUSPLUS_API MeshMetadata read_nastran_op2_metadata(const std::string& rPath,
                                                          const ReadOptions& rOpts = {});

/**
 * @brief Write a mesh as a one-step OP2 result file.
 *
 * The model is the points (`GEOM1` GRIDs, numbered from 1), the cells (`GEOM2`
 * element cards: `triangle` CTRIA3, `quad` CQUAD4, `tetra` CTETRA, `pyramid`
 * CPYRAM, `wedge` CPENTA, `hexahedron` CHEXA, their quadratic forms, a `line`
 * as CROD, CTUBE (property `PTUBE_<pid>`) or, with no property, CONROD, and a
 * `vertex` as CONM2), with `nastran:eid` and `nastran:pid` when the mesh has
 * them (else numbered from 1, and one property id per family), one `EPT`
 * property card per property id and a unit `MAT1`. The results are the point
 * data `DISPLACEMENT`, `EIGENVECTOR`, `VELOCITY`, `ACCELERATION`, `SPC_FORCE`,
 * `MPC_FORCE` and `APPLIED_LOAD` (each `(n, 3)` with a `_ROT` twin) and
 * `TEMPERATURE`, and the centre `STRESS:<M>`/`STRAIN:<M>` cell data of the
 * linear shells (33, 74), solids (39, 67, 68, 255; the corners from
 * `<name>@corner`, else the centre value) and rods (1, 3, 10), under the
 * reader's member names. A member the mesh lacks, or a NaN in a written row, is
 * written as 0; a row with every member NaN, or a point with every component
 * NaN, is left out. The step is `nastran:subcase`, `nastran:analysis` and
 * `nastran:mode` of the field data (default 1, 1 or 6 by the time, 1) and
 * `meshio:time` (0 when absent). Other arrays are dropped with a warning. The
 * Python twin writes the same bytes.
 *
 * @throws WriteError on a non-positive or repeated element id, a property id
 *         shared by shells, solids and rods, an array of the wrong size, or a
 *         file that cannot be opened
 */
MESHIOPLUSPLUS_API void write_nastran_op2(const std::string& rPath, const Mesh& rMesh);

/**
 * @brief An OP2 file written step by step: the model from the first step, then
 *        a block per step in each result table.
 *
 * The tables hold every step, so the steps are kept and the file is written by
 * `Finalize`. Every step must have the first step's points count, cell blocks
 * and result arrays (a `WriteError` otherwise), and a subcase, analysis and
 * mode or time of its own. The time of a step is the argument, not its
 * `meshio:time`. One step is `write_nastran_op2`'s bytes. The Python twin,
 * `meshioplusplus.nastran_op2.SeriesWriter`, writes the same bytes.
 */
class MESHIOPLUSPLUS_API NastranOp2SeriesWriter {
public:
    /// @throws WriteError when @p rPath cannot be opened.
    explicit NastranOp2SeriesWriter(const std::string& rPath);
    ~NastranOp2SeriesWriter();

    NastranOp2SeriesWriter(const NastranOp2SeriesWriter&) = delete;
    NastranOp2SeriesWriter& operator=(const NastranOp2SeriesWriter&) = delete;
    NastranOp2SeriesWriter(NastranOp2SeriesWriter&&) noexcept;
    NastranOp2SeriesWriter& operator=(NastranOp2SeriesWriter&&) noexcept;

    /**
     * @brief Add one step.
     * @throws WriteError when the step's points, cells or result arrays differ
     *         from the first step's, its subcase, analysis and mode or time
     *         repeat an earlier step's, or on a moved-from writer.
     */
    void Write(double Time, const Mesh& rMesh);

    /// The number of steps added so far.
    std::size_t NumSteps() const noexcept;

    /// Write the file. Idempotent. @throws WriteError when no step was added.
    void Finalize();

private:
    struct Impl;
    std::unique_ptr<Impl> mpImpl;
};

}  // namespace meshioplusplus

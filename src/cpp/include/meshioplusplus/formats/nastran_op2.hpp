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
 * See doc/formats/nastran_op2.md.
 */

// System includes
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

}  // namespace meshioplusplus

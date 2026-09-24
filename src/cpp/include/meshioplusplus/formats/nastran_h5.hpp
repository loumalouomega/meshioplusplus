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
 * @file nastran_h5.hpp
 * @brief MSC Nastran HDF5 result database (`.h5`) C++ reader.
 *
 * The file MSC Nastran writes with `MDLPRM,HDF5,...` (since 2016): compound-type
 * tables under `/NASTRAN/INPUT` (the model) and `/NASTRAN/RESULT` (the results),
 * with `/INDEX/NASTRAN/RESULT/...` giving each result table's row range per
 * result domain. Read-only.
 *
 *  - `/NASTRAN/INPUT/NODE/GRID` gives the points, in file order, moved to the
 *    basic system through the `COORDINATE_SYSTEM` CORD1R/C/S and CORD2R/C/S
 *    tables when `CP != 0` (`detail/nastran_model.hpp`); nodal vectors and grid
 *    point forces of a GRID with `CD != 0` are rotated to basic. `CP`/`CD`
 *    become the point data `nastran:cp`/`nastran:cd` when any is non-zero.
 *    SPOINT/EPOINT scalar points are not points.
 *  - Every `/NASTRAN/INPUT/ELEMENT/<CARD>` table with a cell type becomes up to two
 *    cell blocks (linear, then quadratic: a quadratic card with no mid-side nodes
 *    is linear, one with only some is read as linear with a warning). Integer
 *    cell data `nastran:eid` and `nastran:pid` (-1 when the card has none); one
 *    cell region per property id, named `<PTYPE>_<pid>` after the property table
 *    that defines it. Scalar-point cards (CELAS, CDAMP, CMASS) are skipped.
 *  - Each result domain (one row of `/NASTRAN/RESULT/DOMAINS`) that an INDEX table
 *    references is one step for the sequence engine: `ReadOptions::mTimeStep`
 *    selects it, and only that domain's rows are read. `meshio:time` is the
 *    domain's `TIME_FREQ_EIGR` (time, frequency or eigenvalue, as written);
 *    `nastran:domain`, `nastran:subcase`, `nastran:step`, `nastran:analysis`,
 *    `nastran:mode` and `nastran:eigi` are field data.
 *  - `/NASTRAN/RESULT/NODAL/<T>` tables with one row per node become point data:
 *    `X,Y,Z` -> `<T>`, `RX,RY,RZ` -> `<T>_ROT`, the `_CPLX` tables' `XR..`/`XI..`
 *    -> `<T>_real`/`<T>_imag`, a lone `VALUE` -> `<T>`, any other float member
 *    `M` -> `<T>:<M>`; NaN at the points a table has no row for.
 *  - `/NASTRAN/RESULT/ELEMENTAL/<G>/<T>` tables with one row per element become
 *    cell data `<G>:<M>` per float member (the first entry of an array member:
 *    the centre of a solid's or a corner-output shell's values), NaN elsewhere.
 *  - The other values of an element are `(cells, columns)` cell data with
 *    `field_data["nastran:layout:<name>"]`: `<G>:<M>@corner` (the corners, at
 *    the GRID's position in the cell), `<G>:<M>@ply` (`_COMP` plies),
 *    `<G>:<M>@station` (beam and bar stations) and `GRID_FORCE:<M>` (per
 *    element node); `GRID_FORCE` rows of no element are point data
 *    `GRID_FORCE:<label>:<M>`. Other repeated rows are skipped with a warning.
 *
 * A `.h5` without `/NASTRAN/INPUT/NODE/GRID`, or one whose `/NASTRAN` `VERSION`
 * attribute does not name MSC, is refused: other vendors' HDF5 schemas differ.
 * See doc/formats/nastran_h5.md.
 */

#ifdef MESHIOPLUSPLUS_HAS_HDF5

// System includes
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/read_options.hpp"

namespace meshioplusplus {

/**
 * @brief Read an MSC Nastran HDF5 result database.
 *
 * @param rPath filesystem path to read
 * @param rOpts `mTimeStep` selects the result domain (`ResolveTimeStep`);
 *        `mPointsOnly` and `mDataArrays` narrow the result tables that are read
 * @return the model, with the selected domain's results as point and cell data
 * @throws ReadError if the file is not an MSC Nastran HDF5 file, an element
 *         references an undefined GRID, ids repeat, a result table names a
 *         domain `/NASTRAN/RESULT/DOMAINS` lacks, or the step is out of range
 */
MESHIOPLUSPLUS_API Mesh read_nastran_h5(const std::string& rPath, const ReadOptions& rOpts = {});

/**
 * @brief Summarize an MSC Nastran HDF5 file: the model, step 0's result names,
 *        and every result domain's `TIME_FREQ_EIGR` as `mTimeValues`.
 */
MESHIOPLUSPLUS_API MeshMetadata read_nastran_h5_metadata(const std::string& rPath,
                                                         const ReadOptions& rOpts = {});

}  // namespace meshioplusplus

#endif  // MESHIOPLUSPLUS_HAS_HDF5

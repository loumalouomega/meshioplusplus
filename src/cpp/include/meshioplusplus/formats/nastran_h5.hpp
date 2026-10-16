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
 * result domain. Written too (`write_nastran_h5`): the model and the nodal and
 * centre element results, one domain per mesh or per series step.
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
#include <memory>
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

/**
 * @brief Write a mesh as an MSC Nastran HDF5 result database: one domain.
 *
 * The reader's inverse for the core of the schema: `/NASTRAN/INPUT` holds the
 * GRIDs (ids 1..n), one table per element card (`nastran:eid`, `nastran:pid`
 * give the ids; a quadratic solid's node list is zero-padded and permuted to
 * Nastran's order), the property and material tables the elements name; the
 * `DOMAINS` table has one row from `meshio:time` and the `nastran:*` field
 * data; point data becomes the nodal tables (`<T>` and `<T>_ROT` of 3
 * components, a scalar as `VALUE`, `<T>:<M>` members, `_real`/`_imag` pairs as
 * `<T>_CPLX`) and `<G>:<M>` cell data the element tables of the centre values,
 * each result table with its `/INDEX` twin. The file has the `SCHEMA` attribute
 * and no `VERSION` one. The Python twin (`h5py`) writes the same tables.
 *
 * @throws WriteError on a non-positive or repeated element id, a property id
 *         shared by shells, solids and rods, an array of the wrong size, or a
 *         file that cannot be created
 */
MESHIOPLUSPLUS_API void write_nastran_h5(const std::string& rPath, const Mesh& rMesh);

/**
 * @brief A results file written step by step: the model from the first step,
 *        a domain per step. The file is written by `Finalize`.
 *
 * Every step must have the first step's points count, cell blocks and result
 * tables (a `WriteError` otherwise). The time of a step is the argument, not
 * its `meshio:time`. One step is `write_nastran_h5`'s tables. The Python twin,
 * `meshioplusplus.nastran_h5.SeriesWriter`, writes the same tables.
 */
class MESHIOPLUSPLUS_API NastranH5SeriesWriter {
public:
    /// @throws WriteError when @p rPath cannot be created.
    explicit NastranH5SeriesWriter(const std::string& rPath);
    ~NastranH5SeriesWriter();

    NastranH5SeriesWriter(const NastranH5SeriesWriter&) = delete;
    NastranH5SeriesWriter& operator=(const NastranH5SeriesWriter&) = delete;
    NastranH5SeriesWriter(NastranH5SeriesWriter&&) noexcept;
    NastranH5SeriesWriter& operator=(NastranH5SeriesWriter&&) noexcept;

    /// Add one step. @throws WriteError when it differs from the first step's model or tables.
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

#endif  // MESHIOPLUSPLUS_HAS_HDF5

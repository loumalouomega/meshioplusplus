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
 * @file abaqus_fil.hpp
 * @brief Abaqus results file (`.fil`) reader, binary and ASCII.
 *
 * The results file Abaqus writes on request (`*NODE FILE`, `*EL FILE`,
 * `*FILE FORMAT`): the one route into Abaqus results that needs neither the ODB
 * API nor an Abaqus install. It is a sequence of records `[length, key,
 * attributes...]` of 8-byte words. Binary files hold them in 512-word blocks,
 * each framed as a Fortran record (`detail/fortran_records.hpp`), in the
 * writing machine's byte order; ASCII files (`*FILE FORMAT, ASCII`) start each
 * record with `*` and tag each item `I` (two-digit width, then the digits), `D`
 * (a 22-character real) or `A` (8 characters), on 80-column lines.
 *
 * The mesh comes from the model records: 1901 nodes (their labels kept as
 * `point_data["abaqus:id"]`), 1900/1990 elements (labels as
 * `cell_data["abaqus:id"]`; the type by `detail/abaqus_types.hpp`, element
 * types with no meshio++ cell -- user elements, connectors -- are skipped with
 * a warning), 1931/1932 node sets and 1933/1934 element sets as `Point` and
 * `Cell` regions. A set label longer than 8 characters is an integer resolved
 * through the 1940 cross-reference records; assembly labels keep Abaqus's
 * `ASSEMBLY_INSTANCE_SET` spelling.
 *
 * Every increment (2000 ... 2001) is a step: `time_step` selects one (0 =
 * first, negative from the end), its total time is `field_data["meshio:time"]`
 * with `abaqus:step`, `abaqus:increment`, `abaqus:step_time` and
 * `abaqus:procedure`. Its records become data named by the Abaqus identifier
 * (`U`, `RF`, `S`, `E`, `SINV`, `PEEQ`, ...; `key_<n>` for keys without one):
 *  - nodal records (after a 1911 nodal output request) -> `point_data`, NaN for
 *    nodes without a value;
 *  - element records follow a record-1 header (element, integration or node
 *    point, section point, location): at integration points or element nodes
 *    -> `cell_data` of shape `(cells, points * components)`, point-major, the
 *    same width in every block (NaN-padded) and its `(points, components)` in
 *    `field_data["abaqus:layout:<name>"]`; at the centroid or for the whole
 *    element -> `(cells, components)`; averaged at the nodes -> `point_data`.
 *    A single column drops its axis. Components keep the file's order. Section
 *    points above 1 (shell and beam layers; continuum elements write 0) get
 *    `@sp<k>` appended to the name. Rebar, contact,
 * modal-generalised, element matrix and substructure records are skipped.
 *
 * See doc/formats/abaqus_fil.md.
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
 * @brief Read an Abaqus results file.
 * @param rPath filesystem path to read
 * @param rOpts `mTimeStep` selects the increment; `mArrays`/`mPointsOnly`
 *        narrow the data
 * @return the mesh with the selected increment's results
 * @throws ReadError if the file can't be read or is neither framing, a record
 *         is truncated, or `mTimeStep` is out of range
 */
MESHIOPLUSPLUS_API Mesh read_abaqus_fil(const std::string& rPath, const ReadOptions& rOpts = {});

/**
 * @brief The total time of every increment, in file order.
 */
MESHIOPLUSPLUS_API std::vector<double> abaqus_fil_time_values(const std::string& rPath);

/**
 * @brief Metadata (counts, data names, `mTimeValues`) of an Abaqus results file.
 */
MESHIOPLUSPLUS_API MeshMetadata read_abaqus_fil_metadata(const std::string& rPath,
                                                         const ReadOptions& rOpts);

/**
 * @brief Write a mesh as a one-increment Abaqus results file.
 *
 * The inverse of `read_abaqus_fil`. The model is the points (`abaqus:id` labels,
 * else 1 to n), the cells (`abaqus:id` labels over the written cells, else 1 to
 * n; a cell type becomes the Abaqus element `T3D2`, `S3`, `S4R`, `STRI65`,
 * `S8R`, `S9R5`, `C3D4`, `C3D10`, `C3D5`, `C3D13`, `C3D6`, `C3D15`, `C3D8R` or
 * `C3D20R`) and the point and cell regions as node and element sets (a name
 * over 8 characters through a 1940 label). The increment holds the point data
 * named like an Abaqus nodal result (`U`, `V`, `A`, `RF`, `key_<n>`...) and the
 * cell data named like an element result (`S`, `E`, `MISES`, `NFORC`,
 * `key_<n>`...): an array with an `abaqus:layout:<name>` of
 * `(points, components)` is per integration point (per element node for
 * `NFORC` and the `@el` names), any other one set of values per element; a
 * point-data array named like an element result is averaged at the nodes. A
 * `@sp<k>` suffix is the section point, `@rebar:<name>` a rebar. The
 * increment's time is `meshio:time` (0 when absent), its step, number, step time
 * and procedure the `abaqus:step`, `abaqus:increment`, `abaqus:step_time` and
 * `abaqus:procedure` field data, and a `abaqus:mode` (with its eigenvalue,
 * generalized mass, damping and participation factors) a 1980 record. NaN rows
 * are left out, trailing NaN components are trimmed and an inner NaN is
 * written as 0. With @p Binary the file is the little-endian binary layout
 * (512-word blocks), else ASCII. Cells with no Abaqus element, side regions and
 * other arrays are dropped with a warning. The Python twin writes the same
 * bytes.
 *
 * @throws WriteError on a label that is not positive and unique, a cell or
 *         region naming a point the mesh does not have, a value outside what a
 *         `D22.15` item holds, a layout that does not match its array, or a
 *         file that cannot be opened
 */
MESHIOPLUSPLUS_API void write_abaqus_fil(const std::string& rPath, const Mesh& rMesh,
                                         bool Binary = false);

/**
 * @brief Options of the `.fil` writers: ASCII or binary, and for the binary layout its byte
 *        order and the size of its record markers (4 or 8 bytes). Both orders and both marker
 *        sizes are read.
 */
struct AbaqusFilWriteOptions {
    bool mBinary = false;
    bool mBigEndian = false;
    int mMarkerSize = 4;
};

/// `write_abaqus_fil` with the binary layout's byte order and marker size.
MESHIOPLUSPLUS_API void write_abaqus_fil(const std::string& rPath, const Mesh& rMesh,
                                         const AbaqusFilWriteOptions& rOptions);

/**
 * @brief A results file written step by step: the model with the first step,
 *        then an increment per step.
 *
 * The file is written by `Finalize`. Every step must have the first step's
 * points and cell blocks (a `WriteError` otherwise). The time of a step is the
 * argument, not its `meshio:time`. One step is `write_abaqus_fil`'s bytes. The
 * Python twin, `meshioplusplus.abaqus_fil.SeriesWriter`, writes the same bytes.
 */
class MESHIOPLUSPLUS_API AbaqusFilSeriesWriter {
public:
    /// @throws WriteError when @p rPath cannot be opened.
    explicit AbaqusFilSeriesWriter(const std::string& rPath, bool Binary = false);

    /// As above, with the binary layout's byte order and marker size.
    AbaqusFilSeriesWriter(const std::string& rPath, const AbaqusFilWriteOptions& rOptions);
    ~AbaqusFilSeriesWriter();

    AbaqusFilSeriesWriter(const AbaqusFilSeriesWriter&) = delete;
    AbaqusFilSeriesWriter& operator=(const AbaqusFilSeriesWriter&) = delete;
    AbaqusFilSeriesWriter(AbaqusFilSeriesWriter&&) noexcept;
    AbaqusFilSeriesWriter& operator=(AbaqusFilSeriesWriter&&) noexcept;

    /// Add one step. @throws WriteError as above, or on a moved-from writer.
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

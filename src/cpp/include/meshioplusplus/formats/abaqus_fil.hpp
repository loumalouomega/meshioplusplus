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

}  // namespace meshioplusplus

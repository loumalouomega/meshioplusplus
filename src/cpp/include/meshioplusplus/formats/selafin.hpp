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
 * @file formats/selafin.hpp
 * @brief TELEMAC Selafin (`.slf`) results (roadmap 2.4.3.4).
 *
 * Big-endian Fortran sequential records (`detail/fortran_records.hpp`): an
 * 80-character title (its last eight `SERAFIN ` or `SERAFIND`), `NBV1 NBV2`, one
 * 32-character name-and-unit record per variable, `IPARAM` (ten integers; a date
 * record follows when `IPARAM[9] == 1`), `NELEM NPOIN NDP 1`, the 1-based `IKLE`
 * connectivity, `IPOBO` and the `X` and `Y` coordinates; then per step a time
 * record and one record of `NPOIN` values per variable. Single or double precision
 * is the size of the `X` record.
 *
 * `NDP` 3 is a triangle block, 4 a quad block (a tetra block when `IPARAM[6]`, the
 * plane count, exceeds 1), 6 a wedge block (TELEMAC-3D prisms). Each variable is
 * point data under its name (trailing blanks trimmed), its unit in
 * `field_data["selafin:unit:<name>"]` (bytes); in a 3-D file the `ELEVATION Z`
 * (or `COTE Z`) variable of the chosen step is also the points' z. `IPOBO` is the
 * Int32 point data `selafin:ipobo`, the title `field_data["selafin:title"]`.
 * `ReadOptions::mTimeStep` picks the step, its time `field_data["meshio:time"]`.
 *
 * The writer emits a 2-D double-precision file (`SERAFIND`) with one step: one
 * triangle or quad block, the floating-point point data as variables (names cut
 * to 16 characters), `selafin:ipobo` when present, the time from `meshio:time`
 * (0 otherwise). The Python twin is `meshioplusplus.selafin`.
 * @see doc/formats/selafin.md
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/read_options.hpp"

namespace meshioplusplus {

/** @brief Read a Selafin file's mesh and one time step. */
MESHIOPLUSPLUS_API Mesh read_selafin(const std::string& rPath, const ReadOptions& rOptions = {});

/** @brief The summary, with every step's time. */
MESHIOPLUSPLUS_API MeshMetadata read_selafin_metadata(const std::string& rPath,
                                                      const ReadOptions& rOptions);

/** @brief Write a 2-D triangle or quad mesh as a one-step double-precision Selafin file. */
MESHIOPLUSPLUS_API void write_selafin(const std::string& rPath, const Mesh& rMesh);

}  // namespace meshioplusplus

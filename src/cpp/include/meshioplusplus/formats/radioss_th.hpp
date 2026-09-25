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
 * @file radioss_th.hpp
 * @brief OpenRadioss / Radioss time-history file (`<run>T01`, `T02`...)
 *        reader.
 *
 * A time-history file is a big-endian Fortran unformatted file: a title, the
 * run's date, the part, material, property, subset and TH-group descriptions
 * (ids, titles and variable codes), then one output per time, each the time,
 * the global variables, the part and subset variables and one record per
 * TH group (entities x variables). The layout follows OpenRadioss's
 * `th_to_csv` converter (MIT); no code is copied.
 *
 * Each output is a step with no points (an empty `vertex` block) whose field
 * data holds `meshio:time`, the global variables
 * (`radioss_th:global:<name>`), the part and subset variables
 * (`radioss_th:part:<id>:<code name>`, `radioss_th:subset:<id>:...`) and,
 * per TH group, `radioss_th:<kind>:<id>` (entities x variables) with its
 * `:ids` and `:variables` (the variable codes).
 * See doc/formats/radioss_th.md.
 */

// System includes
#include <cstddef>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/read_options.hpp"

namespace meshioplusplus {

/**
 * @brief Whether @p rPath's basename is a time-history file's: a stem, then
 *        `T` and two digits (`crashT01`), no extension.
 */
MESHIOPLUSPLUS_API bool is_radioss_th_filename(const std::string& rPath);

/**
 * @brief Whether the first bytes of a file are a time-history file's title
 *        record: an 84-byte big-endian record holding a plausible file
 *        version and 80 characters.
 */
MESHIOPLUSPLUS_API bool is_radioss_th_head(const char* pHead, std::size_t Size);

/**
 * @brief Read one output of a time-history file.
 * @param rPath the file
 * @param rOpts `mTimeStep` selects the output; `mArrays`/`mPointsOnly` narrow
 *        the data
 * @throws ReadError if the file is not a time-history file, is truncated
 *         before its first output, or `mTimeStep` is out of range
 */
MESHIOPLUSPLUS_API Mesh read_radioss_th(const std::string& rPath, const ReadOptions& rOpts = {});

/**
 * @brief The time of every complete output (a run stopped mid-write leaves a
 *        last output that is not read).
 */
MESHIOPLUSPLUS_API std::vector<double> radioss_th_time_values(const std::string& rPath);

/**
 * @brief Metadata (counts, data names, `mTimeValues`) of a time-history file.
 */
MESHIOPLUSPLUS_API MeshMetadata read_radioss_th_metadata(const std::string& rPath,
                                                         const ReadOptions& rOpts = {});

}  // namespace meshioplusplus

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
 *        reader and writer.
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
 *
 * The writer (`write_radioss_th`, `RadiossThSeriesWriter`) is the reader's
 * inverse: it takes those same field-data names back to the variable codes and
 * groups, and ignores points and cells. The framing is the shared Fortran
 * record framer (`detail/fortran_records.hpp`).
 * See doc/formats/radioss_th.md.
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

/**
 * @brief Write a mesh's `radioss_th:*` field data as a one-output time-history
 *        file.
 *
 * Points and cells are not stored (a time-history file holds none). The output's
 * time is `meshio:time` (0 when absent). The variables are the reader's inverse:
 * `radioss_th:global:<name>`, `radioss_th:part:<id>:<name>` and
 * `radioss_th:subset:<id>:<name>` scalars (a name outside the reader's tables is
 * `var<code>`), and per TH group `radioss_th:<kind>:<id>` (entities x variables)
 * with its optional `:ids` and `:variables` (else 1..n). Everything is stored
 * as big-endian `float32`; ids and codes as `int32`. With
 * `radioss_th:unit_factors` (three values) the file is version 4021, else 3040
 * (the version OpenRadioss writes).
 * Parts, subsets and groups are written sorted by id (groups by kind code, then
 * id) and their variables by code, so the bytes depend on the names only.
 * Other field data is dropped with a warning. The Python twin writes the same
 * bytes.
 *
 * @throws WriteError if a group array is not 2-D or disagrees with its `:ids`
 *         or `:variables`, a value is not an integer in `int32` range, or the
 *         file cannot be opened
 */
MESHIOPLUSPLUS_API void write_radioss_th(const std::string& rPath, const Mesh& rMesh);

/**
 * @brief A time-history file written output by output: the descriptions from
 *        the first step, then one output per `Write`.
 *
 * Every step must have the first step's variables (the same names, ids and
 * codes); a step that does not is a `WriteError` (write one file per step with
 * `{step}` instead). The time of a step is the argument, not its `meshio:time`.
 * One step is `write_radioss_th`'s bytes. The Python twin,
 * `meshioplusplus.radioss_th.SeriesWriter`, writes the same bytes.
 */
class MESHIOPLUSPLUS_API RadiossThSeriesWriter {
public:
    /// @throws WriteError when @p rPath cannot be opened.
    explicit RadiossThSeriesWriter(const std::string& rPath);
    ~RadiossThSeriesWriter();

    RadiossThSeriesWriter(const RadiossThSeriesWriter&) = delete;
    RadiossThSeriesWriter& operator=(const RadiossThSeriesWriter&) = delete;
    RadiossThSeriesWriter(RadiossThSeriesWriter&&) noexcept;
    RadiossThSeriesWriter& operator=(RadiossThSeriesWriter&&) noexcept;

    /**
     * @brief Write one output (the first also writes the descriptions).
     * @throws WriteError when the step's variables differ from the first
     *         step's, on a write failure, or on a moved-from writer.
     */
    void Write(double Time, const Mesh& rMesh);

    /// The number of steps written so far.
    std::size_t NumSteps() const noexcept;

    /// Close the file. Idempotent. @throws WriteError when no step was written.
    void Finalize();

private:
    struct Impl;
    std::unique_ptr<Impl> mpImpl;
};

}  // namespace meshioplusplus

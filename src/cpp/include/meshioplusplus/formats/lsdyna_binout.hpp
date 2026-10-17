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
 * @file lsdyna_binout.hpp
 * @brief LS-DYNA binary output database (`binout`) reader.
 *
 * A binout is an LSDA file (Livermore Software Data Archival): a symbol table
 * of directories and typed variables whose values sit in DATA records, with
 * continuation files `binout%001`... Each ASCII database LS-DYNA would
 * otherwise write (`nodout`, `glstat`, `matsum`, `elout`...) is a directory
 * holding `metadata` (ids, title...) and one folder per output, `d000001`,
 * `d000002`..., each with its `time` and variables.
 *
 * The steps are `nodout`'s outputs (else the first database's): each is a
 * point cloud (one `vertex` block) of the `nodout` nodes at their current
 * coordinates, with `lsdyna:nid` their ids and `displacement`, `rotation`,
 * `velocity`, `rotational_velocity`, `acceleration` and
 * `rotational_acceleration` as point data. Every other database's latest
 * output at or before the step's time is field data
 * `binout:<database>:<variable>` (`binout:elout/beam:axial`), with its
 * `binout:<database>:time` and, from its metadata, `binout:<database>:ids`.
 *
 * The writer (`write_lsdyna_binout`, `LsdynaBinoutSeriesWriter`) is the
 * reader's inverse: one little-endian LSDA file with 8-byte lengths and
 * offsets, reals stored as `float32` and integers as `int32`.
 *
 * See doc/formats/lsdyna_binout.md.
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
 * @brief Read one output of a binout.
 * @param rPath the binout (its `%001`... continuations are found beside it)
 * @param rOpts `mTimeStep` selects the output; `mArrays`/`mPointsOnly` narrow
 *        the data
 * @throws ReadError if the file is not an LSDA file, has no database with
 *         outputs, is truncated, or `mTimeStep` is out of range
 */
MESHIOPLUSPLUS_API Mesh read_lsdyna_binout(const std::string& rPath, const ReadOptions& rOpts = {});

/**
 * @brief The time of every step (`nodout`'s outputs, else the first
 *        database's).
 */
MESHIOPLUSPLUS_API std::vector<double> lsdyna_binout_time_values(const std::string& rPath);

/**
 * @brief Metadata (counts, data names, `mTimeValues`) of a binout.
 */
MESHIOPLUSPLUS_API MeshMetadata read_lsdyna_binout_metadata(const std::string& rPath,
                                                            const ReadOptions& rOpts = {});

/**
 * @brief Whether a path names a binout: its file name is `binout`, or an MPP
 *        run's `binout0000`... (any case).
 */
MESHIOPLUSPLUS_API bool is_binout_filename(const std::string& rPath);

/**
 * @brief Whether the first bytes of a file are an LSDA header followed by its
 *        symbol-table offset command.
 */
MESHIOPLUSPLUS_API bool is_binout_head(const char* pHead, std::size_t Size);

/**
 * @brief Write a mesh as a one-output binout.
 *
 * `nodout` is the points (`x_coordinate`, `y_coordinate`, `z_coordinate`),
 * their ids (`lsdyna:nid`, else 1 to n), the vectors `displacement`,
 * `rotation`, `velocity`, `rotational_velocity`, `acceleration` and
 * `rotational_acceleration` (`n x 3`), every other one-dimensional per-point
 * array (for 2 or more points) under its own name, and the
 * `binout:nodout:<variable>` field data. Every other database is
 * `binout:<database>:<variable>` field data (`<database>` may nest with `/`,
 * `binout:elout/beam:axial`), with `binout:<database>:time` its output's time
 * (else the step's) and `binout:<database>:ids` its metadata ids. Reals are
 * stored as `float32` and integers as `int32`; the output's time is
 * `meshio:time` (0 when absent). Cells, regions and other arrays are dropped
 * with a warning. The Python twin writes the same bytes.
 *
 * @throws WriteError on an id outside the 32-bit range, a variable or
 *         database name a binout cannot hold, a variable repeated, a mesh with
 *         no point and no `binout:*` field data, or a file that cannot be
 *         opened
 */
MESHIOPLUSPLUS_API void write_lsdyna_binout(const std::string& rPath, const Mesh& rMesh);

/**
 * @brief Options of the binout writers.
 *
 * `mSplitSize`: the data bytes after which a new continuation file (`<path>%001`,
 * `%002`...) starts (0: one file). `mWordSize`: the bytes of every length and offset in
 * the LSDA records, 8 or 4 (a file under 4 GiB). `mBigEndian`: the byte order of the
 * records and the data. The Python twin takes the same options.
 */
struct LsdynaBinoutWriteOptions {
    std::size_t mSplitSize = 0;
    int mWordSize = 8;
    bool mBigEndian = false;
};

/// `write_lsdyna_binout` with options.
MESHIOPLUSPLUS_API void write_lsdyna_binout(const std::string& rPath, const Mesh& rMesh,
                                            const LsdynaBinoutWriteOptions& rOptions);

/**
 * @brief `write_lsdyna_binout` split into continuation files.
 *
 * When @p SplitSize is not 0, a new file (`<path>%001`, `%002`...) starts once the
 * data written to the current one would pass @p SplitSize bytes; an output folder
 * stays whole in one file (so one folder larger than the split makes a file of its
 * own), and the symbol tables add a little more. `read_lsdyna_binout` merges the
 * files. 0 writes one file.
 */
MESHIOPLUSPLUS_API void write_lsdyna_binout(const std::string& rPath, const Mesh& rMesh,
                                            std::size_t SplitSize);

/**
 * @brief A binout written step by step: `nodout`'s ids once and an output
 *        folder per step, and any other database's output whenever its own
 *        `binout:<database>:time` changes (so a database written less often
 *        than `nodout` keeps its own cadence).
 *
 * The file is written by `Finalize`. Every step must have the first step's
 * points and ids and, for each database, its first output's variables (a
 * database may be absent from a step and appear later); a step that does not
 * is a `WriteError`. The time of a step is the argument, not its
 * `meshio:time`. One step is `write_lsdyna_binout`'s bytes. The Python twin,
 * `meshioplusplus.lsdyna_binout.SeriesWriter`, writes the same bytes.
 */
class MESHIOPLUSPLUS_API LsdynaBinoutSeriesWriter {
public:
    /// @throws WriteError when @p rPath cannot be opened.
    explicit LsdynaBinoutSeriesWriter(const std::string& rPath);

    /// As above, splitting into `%001`... continuation files at @p SplitSize data bytes
    /// (see the split `write_lsdyna_binout`).
    LsdynaBinoutSeriesWriter(const std::string& rPath, std::size_t SplitSize);

    /// As above, with the word size and byte order of `LsdynaBinoutWriteOptions`.
    LsdynaBinoutSeriesWriter(const std::string& rPath, const LsdynaBinoutWriteOptions& rOptions);
    ~LsdynaBinoutSeriesWriter();

    LsdynaBinoutSeriesWriter(const LsdynaBinoutSeriesWriter&) = delete;
    LsdynaBinoutSeriesWriter& operator=(const LsdynaBinoutSeriesWriter&) = delete;
    LsdynaBinoutSeriesWriter(LsdynaBinoutSeriesWriter&&) noexcept;
    LsdynaBinoutSeriesWriter& operator=(LsdynaBinoutSeriesWriter&&) noexcept;

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

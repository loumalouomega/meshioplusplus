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
 * See doc/formats/lsdyna_binout.md.
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

}  // namespace meshioplusplus

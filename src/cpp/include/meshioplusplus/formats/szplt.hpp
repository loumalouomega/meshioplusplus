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
 * @file szplt.hpp
 * @brief Tecplot SZL (`.szplt`) reader through a user-installed TecIO, the
 *        `szplt` format.
 *
 * SZL is Tecplot's undocumented subzone-loadable format; TecIO, Tecplot's own
 * library, is its only reader. meshio++ never redistributes it: the reader is
 * compiled only with `MESHIOPLUSPLUS_WITH_TECIO=ON` and `TECIO_ROOT` pointing at
 * a TecIO the user already has (a Tecplot 360 install, or a build of the TecIO
 * source Tecplot distributes).
 *
 * TecIO's reader fills the same zone model the ASCII and `.plt` readers decode
 * into (see tecplot.hpp), so a `.szplt` reads exactly as the same data saved as
 * `.plt`: one cell block and Cell region per zone of the selected step, point
 * or cell data by each variable's location, sharing and passive variables
 * resolved, `ReadOptions::mTimeStep` over the distinct solution times.
 * FEPOLYGON/FEPOLYHEDRON zones are refused by name.
 * See doc/formats/szplt.md.
 */

#ifdef MESHIOPLUSPLUS_HAS_TECIO

// System includes
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/read_options.hpp"

namespace meshioplusplus {

/**
 * @brief Read one step of a `.szplt` file through TecIO.
 * @throws ReadError if TecIO cannot open the file, a zone is face-based, `X`
 *         or `Y` is missing, or `mTimeStep` is out of range
 */
MESHIOPLUSPLUS_API Mesh read_szplt(const std::string& rPath, const ReadOptions& rOptions = {});

/** @brief The distinct solution times, ascending; empty for a static file. */
MESHIOPLUSPLUS_API std::vector<double> szplt_time_values(const std::string& rPath);

/** @brief Metadata of step 0 (a points-only read) and `mTimeValues`. */
MESHIOPLUSPLUS_API MeshMetadata read_szplt_metadata(const std::string& rPath,
                                                    const ReadOptions& rOptions = {});

}  // namespace meshioplusplus

#endif  // MESHIOPLUSPLUS_HAS_TECIO

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
 * @file formats/adcirc.hpp
 * @brief ADCIRC `fort.14` grids, with `fort.63`/`fort.64` ASCII results as steps
 * (roadmap 2.4.3.3).
 *
 * `fort.14`: a title, `NE NP`, `NP` lines `id x y depth`, `NE` lines
 * `id 3 n1 n2 n3`, then the open boundaries (`NOPE`, `NETA`, per boundary
 * `NVDLL [IBTYPEE]` and its node ids) and the land boundaries (`NBOU`, `NVEL`, per
 * boundary `NVELL IBTYPE` and its nodes; the barrier columns after a node are
 * skipped). Text after the numbers on a line (ADCIRC's `! comment` or `= NE NP`
 * labels) is ignored. The grid becomes points `(x, y, 0)` with `depth` point data,
 * one triangle block, and point regions `open_boundary_<k>` / `land_boundary_<k>`
 * whose tag is the boundary type. The title is `field_data["adcirc:title"]` (bytes).
 *
 * Results: `fort.63` (elevation, `zeta`) and `fort.64` (depth-averaged velocity,
 * `velocity`, two components) beside a `fort.14` are read automatically, or any
 * `fort.NN` through `AdcircReadOptions::mResults` (`fort.73` `pressure`, `fort.74`
 * `wind_velocity`, others `adcirc:fort.NN`); full and sparse (`TIME IT NNONDEF
 * DEFAULT`) record layouts. `mTimeStep` (and `ReadOptions::mTimeStep` through the
 * registry) picks the dataset; its time is `field_data["meshio:time"]`.
 *
 * The writer emits `fort.14` from a triangle mesh: `depth` point data (0 when
 * absent), the `open_boundary_*` and `land_boundary_*` point regions (a land type
 * whose records carry barrier columns is refused by name). Results are not written.
 * The Python twin is `meshioplusplus.adcirc`.
 * @see doc/formats/adcirc.md
 */

// System includes
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/read_options.hpp"

namespace meshioplusplus {

/** @brief What to read beside an ADCIRC grid. */
struct AdcircReadOptions {
    std::vector<std::string> mResults;  ///< result files; empty reads the fort.63/64 siblings
    bool mSiblings = true;              ///< read fort.63/fort.64 beside a file named fort.14
    long long mTimeStep = 0;            ///< 0 the first dataset, negatives from the end
};

/** @brief Read an ADCIRC grid and one dataset of its results. */
MESHIOPLUSPLUS_API Mesh read_adcirc(const std::string& rPath, const AdcircReadOptions& rOptions);

/** @brief Read an ADCIRC grid through the registry's options (the time step). */
MESHIOPLUSPLUS_API Mesh read_adcirc(const std::string& rPath, const ReadOptions& rOptions = {});

/** @brief The summary, with every dataset's time. */
MESHIOPLUSPLUS_API MeshMetadata read_adcirc_metadata(const std::string& rPath,
                                                     const ReadOptions& rOptions);

/** @brief Write a triangle mesh as an ADCIRC `fort.14` grid. */
MESHIOPLUSPLUS_API void write_adcirc(const std::string& rPath, const Mesh& rMesh);

}  // namespace meshioplusplus

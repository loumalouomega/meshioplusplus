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
 * @file formats/las.hpp
 * @brief ASPRS LAS point clouds, versions 1.2 to 1.4, point formats 0 to 10
 * (roadmap 2.4.6.1); what airborne and mobile scanners write.
 *
 * The reader takes the public header (scale, offset, counts, the 1.4 64-bit point
 * count), keeps the variable-length records verbatim (`field_data["las:vlrs"]`,
 * bytes; the CRS GeoKeys or WKT ride there) and decodes every point record into
 * points (`X Y Z` scaled to doubles) in one `vertex` block, and point data under
 * laspy's names: `intensity`, `return_number`, `number_of_returns`,
 * `scan_direction_flag`, `edge_of_flight_line`, `classification`, `synthetic`,
 * `key_point`, `withheld`, `overlap` and `scanner_channel` (formats 6 to 10),
 * `scan_angle_rank` (0 to 5) or `scan_angle` (6 to 10), `user_data`,
 * `point_source_id`, `gps_time`, `red`, `green`, `blue` and `nir`, each in its LAS
 * type. Wave-packet fields and bytes past the standard record are skipped.
 * `field_data` also carries `las:scale`, `las:offset`, `las:point_format` and
 * `las:version`. `LasReadOptions::mEvery` keeps every n-th point; the
 * `ReadOptions` array filter picks the fields.
 *
 * The writer emits LAS 1.4 with point format 6, 7 (when the mesh has `red`,
 * `green` and `blue`) or 8 (and `nir`): every point of the mesh (its cells are not
 * written), the scale and offset from `las:scale`/`las:offset` (0.001 and the
 * floored minimum otherwise), the kept VLRs, and the fields above the mesh carries
 * (`scan_angle_rank` converted to `scan_angle`'s 0.006-degree steps). The Python
 * twin is `meshioplusplus.las`, which writes the same bytes; laspy reads both.
 * @see doc/formats/las.md
 */

// System includes
#include <cstddef>
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/read_options.hpp"

namespace meshioplusplus {

/** @brief What to keep from a LAS file. */
struct LasReadOptions {
    std::size_t mEvery = 1;  ///< keep every n-th point (1 keeps all)
};

/** @brief Read a LAS point cloud. */
MESHIOPLUSPLUS_API Mesh read_las(const std::string& rPath, const ReadOptions& rOptions,
                                 const LasReadOptions& rLas);

/** @brief Read a LAS point cloud with the registry's options. */
MESHIOPLUSPLUS_API Mesh read_las(const std::string& rPath, const ReadOptions& rOptions = {});

/** @brief Write a mesh's points as LAS 1.4 (point format 6, 7 or 8). */
MESHIOPLUSPLUS_API void write_las(const std::string& rPath, const Mesh& rMesh);

}  // namespace meshioplusplus

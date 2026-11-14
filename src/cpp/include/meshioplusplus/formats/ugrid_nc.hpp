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
 * @file formats/ugrid_nc.hpp
 * @brief UGRID-conventions NetCDF meshes (roadmap 2.4.3.2): what Delft3D-FM, D-Flow,
 * SCHISM, FVCOM and xugrid write.
 *
 * The reader takes the first variable with `cf_role = "mesh_topology"`: its
 * `node_coordinates` (x, y and an optional third), its `face_node_connectivity`
 * (`start_index` honoured, `_FillValue` padding of mixed polygons dropped, the face
 * dimension from `face_dimension` when the array is stored transposed) and its
 * `edge_node_connectivity`. Faces become `triangle`, `quad` and one `polygon` block
 * per larger node count, in that order (each keeping the file's face order); edges a
 * trailing `line` block. Data variables naming the mesh (`mesh`) with `location`
 * `node`, `face` or `edge` become point and cell data (a face array is NaN or 0 on
 * the line block and the reverse); a leading time dimension is a step, chosen with
 * `ReadOptions::mTimeStep`, its time in `field_data["meshio:time"]`, and the
 * metadata reader lists every step's time. Volume topologies are refused by name.
 *
 * The writer emits a netCDF-4 file with one 2-D topology `mesh` (1-D when the mesh
 * has only lines): faces from the triangle, quad and polygon blocks (`_FillValue`
 * -1 padding, `start_index` 0), edges from the line blocks, point data as `node`
 * variables and cell data as `face` variables (and `<name>_edge` `edge` variables
 * carrying `meshioplusplus_name` for their line-block part). Other cells are refused
 * by name. The Python twin is `meshioplusplus.ugrid_nc` (netCDF4).
 * @see doc/formats/ugrid_nc.md
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/read_options.hpp"

namespace meshioplusplus {

/** @brief Read a UGRID NetCDF mesh and one time step of its data. */
MESHIOPLUSPLUS_API Mesh read_ugrid_nc(const std::string& rPath, const ReadOptions& rOptions = {});

/** @brief The mesh's summary, with every time step's value. */
MESHIOPLUSPLUS_API MeshMetadata read_ugrid_nc_metadata(const std::string& rPath,
                                                       const ReadOptions& rOptions);

/** @brief Write a 1-D or 2-D mesh as a UGRID NetCDF file. */
MESHIOPLUSPLUS_API void write_ugrid_nc(const std::string& rPath, const Mesh& rMesh);

}  // namespace meshioplusplus

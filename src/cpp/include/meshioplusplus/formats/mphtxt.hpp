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
 * @file mphtxt.hpp
 * @brief COMSOL native mesh files: text (.mphtxt) and binary (.mphbin).
 *
 * Both serialise the same sequence: a version pair (`0 1`), tag and type
 * tables, then objects (`0 0 1` and a class name). Read are `Mesh` objects
 * (versions 4 and older, whose element records also carry parameter values
 * and up/down pairs) and `Selection` objects; the first object of any other
 * class stops the read with a warning. Several Mesh objects are merged, each
 * one's cells a region named by its tag. `.mphbin` stores integers as
 * little-endian int32, doubles as float64, and a string as its length and
 * one int32 code point per character; in `.mphtxt` a string is its length,
 * a blank and the characters, and `#` starts a comment.
 *
 * Every element's geometric entity index becomes `cell_data["mphtxt:geom"]`
 * (domains count from 1, boundaries, edges and points from 0). A Selection
 * becomes a cell Region: the elements of its dimension whose entity it
 * lists. Node order goes through the node-ordering registry (format key
 * `mphtxt`): COMSOL lists the corners in tensor order, then the other nodes
 * of the quadratic lattice lexicographically. Types: `vtx`, `edg`/`edg2`,
 * `tri`/`tri2`, `quad`/`quad2`, `tet`/`tet2`, `pyr`/`pyr2`,
 * `prism`/`prism2`, `hex`/`hex2` (`quad9`, `pyramid14`, `wedge18` and
 * `hexahedron27` for the quadratic ones). See doc/formats/mphtxt.md.
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Write a Mesh to a COMSOL text mesh (.mphtxt) file.
 *
 * Writes a version 4 Mesh object. Entity indices come from
 * `cell_data["mphtxt:geom"]` when present; otherwise, per dimension, each
 * pairwise-disjoint cell region is an entity in turn and the remaining cells
 * one more (domains from 1, lower dimensions from 0). A cell region that is
 * a union of whole entities of one dimension is written as a Selection;
 * other regions are dropped with a warning, as are data arrays.
 *
 * @param rPath filesystem path to the .mphtxt file to create/overwrite
 * @param rMesh the mesh to write
 * @throws WriteError on a cell type with no COMSOL equivalent, or a region
 *         naming a cell that does not exist
 */
MESHIOPLUSPLUS_API void write_mphtxt(const std::string& rPath, const Mesh& rMesh);

/**
 * @brief Read a COMSOL text mesh (.mphtxt) file into a Mesh.
 *
 * @param rPath filesystem path to the .mphtxt file to read
 * @return the Mesh, with `cell_data["mphtxt:geom"]` and one cell Region per
 *         Selection
 * @throws ReadError on a malformed file, an unknown element type, or a file
 *         without a Mesh object
 */
MESHIOPLUSPLUS_API Mesh read_mphtxt(const std::string& rPath);

/**
 * @brief Write a Mesh to a COMSOL binary mesh (.mphbin) file; the same
 *        content as #write_mphtxt. There is no comment slot, so no
 *        provenance is written.
 * @param rPath filesystem path to the .mphbin file to create/overwrite
 * @param rMesh the mesh to write
 * @throws WriteError as #write_mphtxt, or on a value beyond 32 bits
 */
MESHIOPLUSPLUS_API void write_mphbin(const std::string& rPath, const Mesh& rMesh);

/**
 * @brief Read a COMSOL binary mesh (.mphbin) file into a Mesh; the same
 *        content as #read_mphtxt.
 * @param rPath filesystem path to the .mphbin file to read
 * @return the Mesh
 * @throws ReadError as #read_mphtxt
 */
MESHIOPLUSPLUS_API Mesh read_mphbin(const std::string& rPath);

}  // namespace meshioplusplus

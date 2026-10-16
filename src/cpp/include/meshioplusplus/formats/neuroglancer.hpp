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
 * @file neuroglancer.hpp
 * @brief Neuroglancer precomputed single-resolution mesh fragments (no extension).
 *
 * Little-endian binary: a `uint32` vertex count `n`, `n * 3` `float32` coordinates,
 * then the rest of the file as `uint32` triangle indices (three per triangle).
 *
 * **Mesh mapping.** float32 points plus one `triangle` block (int64 indices). The
 * writer keeps the triangle blocks only, warns about the other cell types, and narrows
 * the points to float32; data arrays are not stored.
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Read a Neuroglancer precomputed mesh fragment.
 * @throws ReadError on a truncated vertex table, a triangle table that is not a
 *         whole number of triangles, or an index outside the vertex table
 */
MESHIOPLUSPLUS_API Mesh read_neuroglancer(const std::string& rPath);

/**
 * @brief Write the triangle blocks of a mesh as a Neuroglancer fragment.
 * @throws WriteError if the file cannot be opened, there are more than 2^32 - 1
 *         points, or a point index does not fit in 32 bits
 */
MESHIOPLUSPLUS_API void write_neuroglancer(const std::string& rPath, const Mesh& rMesh);

}  // namespace meshioplusplus

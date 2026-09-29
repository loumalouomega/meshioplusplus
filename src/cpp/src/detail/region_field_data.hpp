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
 * @file detail/region_field_data.hpp
 * @brief Named regions in a VTK XML file's `<FieldData>` (`.vtu`, `.vtp`).
 *
 * VTK has no named-set concept, so meshio++ writes each region as one
 * dataset-level Int64 field array (v16.27.0, doc/regions.md):
 *
 *  - `region:point:<name>` — point indices, one component;
 *  - `region:cell:<name>` — cell indices, one component;
 *  - `region:side:<name>` — `(cell, local facet)` pairs, two components;
 *  - `region-meta:<kind>:<name>` — `[dim, tag]`, written only when either is
 *    not -1.
 *
 * Cell indices are in the **file's** cell order -- what any VTK tool numbers
 * cells by -- so a writer whose file order is not block-major (`.vtp` groups
 * Verts, Lines, Polys) translates through a map, as does a reader whose
 * blocks are not (`.vtu` buckets polyhedra by node count). A reader removes
 * the arrays it understands from `field_data`; a malformed one (wrong dtype,
 * component count or range) is warned about and left as ordinary field data.
 *
 * Private to the core (not installed): the two formats' readers and writers
 * are its only users.
 */

// System includes
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/read_options.hpp"

namespace meshioplusplus {
namespace detail {

/// Whether @p rName is one of the region arrays (`region:` / `region-meta:`).
bool is_region_field_name(const std::string& rName);

/**
 * @brief The field arrays holding @p rMesh's regions, in canonical region order.
 * @param rMesh The mesh.
 * @param pGlobalToFile Global (block-major) cell -> file cell, or null for
 *        the identity.
 * @return `(name, array)` pairs, each array Int64.
 */
std::vector<std::pair<std::string, NDArray>> regions_to_field_arrays(
    const Mesh& rMesh, const std::vector<std::int64_t>* pGlobalToFile);

/**
 * @brief Turn the region arrays among @p rArrays into regions of @p rMesh.
 *
 * Call after the cells are built. Every array that is not a well-formed
 * region array is added to @p rMesh as field data unchanged.
 *
 * @param rMesh The mesh read so far (points and cells in place).
 * @param rArrays The file's field arrays, in file order; consumed.
 * @param pFileToGlobal File cell -> global cell, or null for the identity.
 * @param pFormat The format name, for the warnings.
 */
void regions_from_field_arrays(Mesh& rMesh, std::vector<std::pair<std::string, NDArray>>& rArrays,
                               const std::vector<std::int64_t>* pFileToGlobal, const char* pFormat);

/**
 * @brief The region summaries a metadata read reports, from the field arrays'
 * names and declared tuple counts alone (no payload is decoded).
 * @param rArrays `(name, NumberOfTuples)` of every field array.
 * @return One summary per region array, dim/tag -1 (they live in the payload).
 */
std::vector<RegionSummary> region_summaries_from_field_names(
    const std::vector<std::pair<std::string, std::size_t>>& rArrays);

}  // namespace detail
}  // namespace meshioplusplus

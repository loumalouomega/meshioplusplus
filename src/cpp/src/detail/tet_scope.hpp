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

#pragma once

/**
 * @file detail/tet_scope.hpp
 * @brief The tet-only scope check `optimize_volume` and `adapt_volume` share:
 * every block linear `tetra`, at least one of them.
 *
 * A **core-private** header: nothing installed names it.
 */

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {
namespace detail {

/**
 * @brief Throw `std::invalid_argument`, its message starting with @p pPrefix,
 * unless every block of @p rMesh is a linear `tetra` block and one exists: a
 * polyhedron or ragged block, another 3D block (pointing at
 * `convert_cells(mode='simplexify')`) and a lower-dimensional block (pointing
 * at `split`) are each named.
 */
void require_tet_only_blocks(const Mesh& rMesh, const char* pPrefix);

}  // namespace detail
}  // namespace meshioplusplus

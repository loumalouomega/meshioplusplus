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
 * @file operations/adapt_volume_internal.hpp
 * @brief Test hooks into `adapt_volume`'s guards: each runs the operation's
 * own setup and evaluates one candidate edit without applying it.
 *
 * Operation-private (not installed).
 */

// System includes
#include <cstdint>
#include <string>

// Project includes
#include "meshioplusplus/operations/adapt_volume.hpp"

namespace meshioplusplus {
namespace detail {

/**
 * @brief The guard that refuses the half-edge collapse of point @p A onto
 * point @p B of @p rMesh under @p rOptions: `"none"` when it would be applied,
 * else `"not_edge"`, `"link"`, `"inversion"`, `"quality"`, `"length"`,
 * `"class"`, `"hausdorff"` or `"normal"`.
 */
std::string adapt_volume_probe_collapse(const Mesh& rMesh, const AdaptVolumeOptions& rOptions,
                                        std::int64_t A, std::int64_t B);

}  // namespace detail
}  // namespace meshioplusplus

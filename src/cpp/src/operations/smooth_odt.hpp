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
 * @file operations/smooth_odt.hpp
 * @brief `smooth`'s ODT relocation pass, callable without building a mesh.
 *
 * An **operation-private** header (the `formats/gid_common.hpp` precedent):
 * it lives beside the `.cpp` files, nothing outside `smooth.cpp` and
 * `optimize_volume.cpp` includes it, and no installed header may name it, so
 * it adds nothing to the API or the ABI.
 *
 * `optimize_volume` used to relocate by building a fresh `Mesh` and calling
 * `smooth()` once per sweep, which rebuilt the node adjacency (unused by ODT),
 * the boundary facet hash, the cell table and the incidence, and copied the
 * points and connectivity in and out. The flips never change the boundary, so
 * the pin mask is computed once; each sweep then rebuilds only the tet table
 * and incidence (the flips change the tets) and runs the same Jacobi pass,
 * with the same arithmetic in the same order, straight on the caller's
 * coordinate buffer (roadmap §4).
 */

// System includes
#include <array>
#include <cstdint>
#include <vector>

// Project includes
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/operations/smooth.hpp"

namespace meshioplusplus {
namespace detail {

/**
 * @brief The nodes `smooth(rTetMesh, rOptions)` would hold still: the
 * caller's `mFrozen`, the boundary when `mFixBoundary`, and feature nodes
 * (a subset of the boundary) -- phase 2 of `smooth()`, with its warnings.
 * @throws std::invalid_argument when `mFrozen` has the wrong length.
 */
std::vector<std::uint8_t> smooth_odt_pin_mask(const Mesh& rTetMesh, const SmoothOptions& rOptions);

/**
 * @brief One ODT pass (`smooth` with `SmoothMethod::Odt`, one iteration, the
 * lambda from `rOptions`) over the tets `rTets`, moving `rXyz` (a flat `(n, 3)`
 * buffer) in place, with `rPinned` held still and the inversion guard as
 * `rOptions.mGuardInversion` says. Bit-identical to `smooth()` on the mesh of
 * those points and tets.
 * @return The node moves the inversion guard rejected.
 */
std::int64_t smooth_odt_pass(std::vector<double>& rXyz,
                             const std::vector<std::array<std::int64_t, 4>>& rTets,
                             const std::vector<std::uint8_t>& rPinned,
                             const SmoothOptions& rOptions);

}  // namespace detail
}  // namespace meshioplusplus

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
 * @file detail/closed_facets.hpp
 * @brief The surface cells of a mesh as planar facets of a closed,
 * consistently wound shell -- what the faceted `step` writer and `fit_cad`
 * both build a solid from.
 *
 * A **core-private** header (the `surface_rings.hpp` precedent): no installed
 * header names it, and it adds nothing to the API or the ABI.
 *
 * The rings come from `surface_rings`. A polygon whose corners lie within
 * `1e-9` of its extent of their Newell plane is kept whole; any other is
 * fanned into triangles from its first corner. Every facet must have a
 * direction (non-zero area), and the facets must close: no open, non-manifold
 * or inconsistently wound edge (`crease_edges`).
 */

// System includes
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {
namespace detail {

/// Planar facets in CSR form, with their source cells, points and normals.
struct ClosedFacets {
    std::vector<std::int64_t> mStart{0};  ///< offsets into `mNodes`, one more than facets
    std::vector<std::int64_t> mNodes;     ///< the facets' point ids, concatenated
    /// The global (block-major) cell id of each facet; a fanned cell repeats it.
    std::vector<std::int64_t> mCell;
    std::vector<double> mXyz;      ///< the mesh's points, three per point
    std::vector<double> mNormals;  ///< the unit normal of each facet, three per facet
    /// Cells that are not surface facets (lines, vertices), skipped.
    std::size_t mNumSkipped = 0;

    std::size_t NumFacets() const { return mCell.size(); }
};

/**
 * @brief Collect the surface cells of @p rSurface as the facets of a closed shell.
 * @param rPrefix prepended to every error message.
 * @return the facets; empty (with `mNumSkipped` set) when there is no surface
 *         cell, which the caller reports.
 * @throws std::invalid_argument on a volume block, a non-finite coordinate, a
 *         degenerate facet, or a surface that is not closed and consistently
 *         wound.
 */
ClosedFacets closed_facets(const Mesh& rSurface, const std::string& rPrefix);

}  // namespace detail
}  // namespace meshioplusplus

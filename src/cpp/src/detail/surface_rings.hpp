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
 * @file detail/surface_rings.hpp
 * @brief The 2-D cells of a mesh as face rings tagged with their global cell
 * ids -- the surface `segment` grows patches over and `fit_primitives` fits.
 *
 * A **core-private** header (the `crease_edges.hpp` precedent): no installed
 * header names it, and it adds nothing to the API or the ABI.
 *
 * Triangles, quads and their quadratic variants contribute their corners (the
 * leading node of each `cell_edges` entry), polygons their whole row. Lines,
 * vertices, VTK Lagrange surfaces and cells naming a missing point are
 * skipped; a volume block (or polyhedron) throws `std::invalid_argument`.
 * The numpy twin is `_segment.surface_rings`.
 */

// System includes
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {
namespace detail {

/// Face rings in CSR form, one per kept surface cell.
struct SurfaceRings {
    std::vector<std::int64_t> mStart{0};  ///< offsets into `mNodes`, one more than faces
    std::vector<std::int64_t> mNodes;     ///< the rings' point ids, concatenated
    std::vector<std::int64_t> mCell;      ///< the global (block-major) cell id of each ring

    std::size_t NumFaces() const { return mCell.size(); }
};

/**
 * @brief Collect the selected 2-D cells as face rings.
 * @param rSelected per global cell, non-zero when the cell is taken.
 * @param rErrorPrefix prepended to the volume-block error message.
 * @throws std::invalid_argument on a volume or polyhedron block.
 */
SurfaceRings surface_rings(const Mesh& rMesh, const std::vector<std::uint8_t>& rSelected,
                           const std::string& rErrorPrefix);

/// The mesh's points as a flat `n * 3` array, padded with zeros below 3-D.
std::vector<double> surface_xyz(const Mesh& rMesh);

}  // namespace detail
}  // namespace meshioplusplus

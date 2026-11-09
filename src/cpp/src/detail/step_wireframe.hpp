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
 * @file detail/step_wireframe.hpp
 * @brief Curves to STEP (roadmap 7.6.5): a line mesh written as a
 * `GEOMETRIC_CURVE_SET` in a `GEOMETRICALLY_BOUNDED_WIREFRAME_SHAPE_REPRESENTATION`.
 *
 * A `fit_curves` sketch is written as its pieces: a line as a `TRIMMED_CURVE`
 * of a `LINE` (unit direction, trimmed at 0 and its length), an arc as a
 * `TRIMMED_CURVE` of a `CIRCLE` whose axis is the sketch normal (its reverse
 * for a clockwise arc, so every trim runs counter-clockwise), a circle as a
 * `CIRCLE`. Any other line mesh is chained (`detail/polyline_chain.hpp`) and
 * each chain written as a `POLYLINE` (a loop's first point repeated).
 *
 * A **core-private** header: nothing installed names it.
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {
namespace detail {

/// Whether @p rMesh has line cells and no cell of two or three dimensions.
bool step_is_wireframe(const Mesh& rMesh);

/// Write @p rMesh's curves (see the file comment).
/// @throws WriteError on a mesh with no curve to write, or a file that cannot be written.
void step_write_wireframe(const std::string& rPath, const Mesh& rMesh);

}  // namespace detail
}  // namespace meshioplusplus

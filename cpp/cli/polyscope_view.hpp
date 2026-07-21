// SPDX-License-Identifier: MIT
/**
 * @file polyscope_view.hpp
 * @brief The CLI's native viewer entry points.
 *
 * Both functions are **always declared and always defined**, so the `view` and
 * `screenshot` verbs exist in every build and are listed in `--help`. Without
 * `-DMESHIOPLUSPLUS_WITH_POLYSCOPE=ON` they throw naming that flag rather than
 * not existing — the same contract `partition_kahip_parts` follows: never a
 * link error, never a silent no-op.
 *
 * Polyscope is attached to the CLI target only, never to
 * `meshioplusplus_core_obj`, so nothing else in the project can acquire an
 * OpenGL dependency through it.
 */
#ifndef MESHIOPLUSPLUS_CLI_POLYSCOPE_VIEW_HPP
#define MESHIOPLUSPLUS_CLI_POLYSCOPE_VIEW_HPP

#include <string>

#include "meshioplusplus/mesh.hpp"
#include "view_payload.hpp"

namespace meshioplusplus::cli {

/// Whether this build can actually show a window.
bool has_polyscope() noexcept;

/**
 * @brief Open an interactive window on `rMesh` and block until it is closed.
 *
 * @param rMesh the mesh to show
 * @param Kind which representation to build
 * @param rColorBy name of a quantity to enable on load, or empty
 * @param rName structure name shown in Polyscope's UI
 * @throws std::runtime_error when built without Polyscope
 * @throws std::invalid_argument on an unmappable mesh, or an unknown `rColorBy`
 */
void view_mesh(const Mesh& rMesh, ViewKind Kind, const std::string& rColorBy,
               const std::string& rName);

/**
 * @brief Render `rMesh` to a PNG without opening a window.
 *
 * @throws std::runtime_error when built without Polyscope, or when no headless
 *         rendering backend is available
 */
void screenshot_mesh(const Mesh& rMesh, ViewKind Kind, const std::string& rColorBy,
                     const std::string& rName, const std::string& rOutPath, int Width,
                     int Height, bool Transparent);

}  // namespace meshioplusplus::cli

#endif  // MESHIOPLUSPLUS_CLI_POLYSCOPE_VIEW_HPP

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
 * @file formats/gltf.hpp
 * @brief glTF 2.0 writer (`.glb` binary container, `.gltf` + `.bin` beside it).
 *
 * A web-native surface export: the browser viewer, three.js, Blender and every
 * dashboard that loads glTF can show a mesh with no meshio++ in the loop. The
 * writer is write-only -- reading glTF is out of scope -- and exports the
 * *surface* of whatever it is given:
 *
 *  - the skin of volume blocks (`extract_skin`, linearised), with a 2-D cell
 *    that coincides with a skin facet winning over it, so a boundary patch such
 *    as `wall` keeps its own cells;
 *  - 2-D cells, triangulated with the same fan `convert_cells(simplexify)` uses
 *    and linearised;
 *  - `line` cells as `LINES` primitives and `vertex` cells (or a cell-less mesh,
 *    which is how a `.pcd`/`.xyz` point cloud arrives) as `POINTS`.
 *
 * ### Normals
 *
 * glTF normals are per vertex, so a crease needs its position twice. Triangles
 * are grouped into smooth fans exactly as `compute_normals` does
 * (`detail::vertex_normal_groups`), and each fan becomes its own vertex.
 *
 * ### Conventions
 *
 * `float32`, metres, right-handed, Y-up. Source data is usually Z-up (CAE) and
 * rarely Y-up, so the axis change is a rotation on the root node rather than a
 * rewrite of the coordinates, and so is the recentring offset (the bounding-box
 * centre is subtracted from every position so a model with large coordinates
 * keeps `float32` precision): the transform is exactly reversible, because it
 * is a permutation with signs, and the coordinates stay small.
 *
 * ### Fields
 *
 * glTF has no scalar-field concept. Every `point_data` array of one to four
 * components is exported raw as an underscore-prefixed custom attribute
 * (`temperature` -> `_TEMPERATURE`), cast to `float32`. Choosing a field with
 * `mColorBy` additionally bakes it through a colormap into `COLOR_0` and marks
 * the material `KHR_materials_unlit`, because shading distorts a colour reading.
 * `COLOR_0` is linear, so the colormap's sRGB bytes go through the inverse
 * transfer function first.
 *
 * ### Structure
 *
 * One glTF node (and mesh) per cell region, named; a cell that is in several
 * regions goes to the smallest one, and a cell in none goes to `unassigned`.
 * A mesh with no cell regions is one node named `mesh`.
 *
 * ### Determinism
 *
 * The JSON is emitted by a small hand-written writer with a fixed key order
 * and `%.17g` numbers, and the Python reference (`gltf/_gltf.py`) does the same
 * arithmetic in the same order, so the two engines write identical bytes.
 */

// System includes
#include <cstdint>
#include <optional>
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/operations/sdf.hpp"

namespace meshioplusplus {

/// The file layout. `Auto` follows the suffix: `.glb` is binary, anything else JSON.
enum class GltfContainer : std::uint8_t { Auto = 0, Binary = 1, Json = 2 };

/// Which source axis points up. `Auto` is Y for a flat mesh (2-D points, or all
/// `z` within 1e-14 of zero -- the SVG writer's rule) and Z otherwise.
enum class GltfUpAxis : std::uint8_t { Auto = 0, Z = 1, Y = 2, X = 3 };

/// Parse a container name: `auto`, `glb`/`binary`, `gltf`/`json`.
MESHIOPLUSPLUS_API GltfContainer gltf_container_from_name(const std::string& rName);

/// Parse an up-axis name: `auto`, `x`, `y`, `z`.
MESHIOPLUSPLUS_API GltfUpAxis gltf_up_axis_from_name(const std::string& rName);

/// What `write_gltf` should do beyond the defaults.
struct GltfWriteOptions {
    GltfContainer mContainer = GltfContainer::Auto;
    GltfUpAxis mUpAxis = GltfUpAxis::Auto;
    /// How incident faces are weighted into a vertex normal.
    SdfPseudonormalWeight mNormalWeight = SdfPseudonormalWeight::Angle;
    /// Write `NORMAL` (and split vertices at creases). Off shares vertices.
    bool mNormals = true;
    /// Export every one-to-four component `point_data` array as `_NAME`.
    bool mFields = true;
    /// Subtract the bounding-box centre from the positions (carried in the node).
    bool mRecenter = true;
    /// One node per cell region; off writes a single node.
    bool mByRegion = true;
    /// Mark the material `KHR_materials_unlit` when colouring.
    bool mUnlit = true;
    /// The largest dihedral angle, in degrees, still smooth; `[0, 180]`.
    double mSplitAngle = 30.0;
    /// Source unit -> metres, carried as the node scale.
    double mScale = 1.0;
    /// A point_data or cell_data array to bake into `COLOR_0`; empty disables.
    std::string mColorBy;
    /// The component to colour by; the magnitude when unset.
    std::optional<int> mComponent;
    /// The colormap: `viridis`, `coolwarm` or `turbo`.
    std::string mCmap = "viridis";
    std::optional<double> mVMin;
    std::optional<double> mVMax;
    /// The colour of a non-finite value, `#rrggbb`.
    std::string mNanColor = "#808080";
};

/**
 * @brief Write @p rMesh as glTF 2.0.
 *
 * @param rPath the file; `.glb` writes the binary container, anything else the
 *        JSON one with `<path minus .gltf>.bin` beside it (see `mContainer`).
 * @throws WriteError on a non-finite coordinate, a file that cannot be opened,
 *         or a container beyond glTF's 32-bit sizes.
 * @throws std::invalid_argument on a bad option (split angle, scale, colormap,
 *         range, colour, component) or an unknown `mColorBy` array.
 */
MESHIOPLUSPLUS_API void write_gltf(const std::string& rPath, const Mesh& rMesh,
                                   const GltfWriteOptions& rOptions = {});

}  // namespace meshioplusplus

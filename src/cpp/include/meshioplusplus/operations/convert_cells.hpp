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
 * @file operations/convert_cells.hpp
 * @brief Dependency-free conversion of a mesh's *element representation*: drop
 * higher-order nodes, decompose cells into simplices, or promote linear cells
 * to serendipity quadratic ones.
 *
 *  - **Linearize**: every higher-order cell becomes its linear base
 *    (`tetra10` -> `tetra`, `hexahedron27` -> `hexahedron`, ...). Corner
 *    connectivity is kept verbatim; nodes that no cell references any more are
 *    pruned and the connectivity / `point_data` remapped. Cell count is
 *    unchanged, so `cell_data` passes straight through.
 *  - **Simplexify**: every cell is decomposed into simplices of the *same*
 *    topological dimension (quad -> 2 triangles, polygon(n) -> (n-2)-triangle
 *    fan, hexahedron -> 6 tetra, wedge -> 3 tetra, pyramid -> 2 tetra). Points
 *    are untouched -- the children reuse the parent's own corner nodes -- and
 *    each parent's `cell_data` row is replicated to its children. Higher-order
 *    input is linearized first.
 *  - **Elevate**: linear cells are promoted to their serendipity (edge-only)
 *    quadratic counterpart (`triangle` -> `triangle6`, `hexahedron` ->
 *    `hexahedron20`, ...), creating one new node per unique edge at the edge
 *    midpoint with `point_data` set to the mean of the edge's endpoints. Cell
 *    count is unchanged; point count grows.
 *
 * All three modes are **idempotent on cells they do not apply to**: a linear
 * block passes through Linearize unchanged, an already-simplex block passes
 * through Simplexify unchanged, and an already-quadratic block passes through
 * Elevate unchanged. That is what makes the operation safe on a mixed-order
 * mesh. Only genuinely unsupported constructs throw: the full-Lagrange targets
 * that need face/body centres (`quad9`, `hexahedron27`) are an explicit
 * non-goal of this version.
 *
 * **Simplexify decomposes a polyhedron** (since v9.17.0) into one tetrahedron
 * per (face, edge-of-that-face), from the face's corner average to the cell's
 * -- the same fan `detail/polyhedron.hpp` measures, so the children's total
 * volume equals the parent's exactly. That adds `1 + numFaces` points per cell
 * rather than a single centroid: a one-point fan about each face's first node
 * is diagonal-dependent for non-planar faces and inverts on any cell whose
 * faces are not star-shaped about that node. A cell whose faces are not a
 * closed orientable surface throws, naming the cell.
 *
 * **Block structure is preserved 1:1 in every mode** -- the output has exactly
 * as many cell blocks as the input, in the same order -- which is what keeps
 * the `cell_data` correspondence (one array per block) trivially correct under
 * every backend.
 *
 * Output is deterministic: the decomposition templates are fixed, and the
 * Elevate mid-edge numbering is assigned by a serial pass over a
 * parallel-filled record buffer (`src/cpp/src/operations/surface.cpp`'s phase-split
 * idiom), so results are byte-identical across mesh backends and thread counts.
 *
 * Everything is standard C++ and the uniform mesh API only, so it compiles
 * under every mesh backend. This is an operation, not a file format -- it is
 * not in the format registry.
 */

// System includes
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/ndarray.hpp"

namespace meshioplusplus {

/// Which element-representation conversion `convert_cells` performs.
enum class ConvertCellsMode {
    Linearize,   ///< Higher-order cells -> their linear base (drops mid nodes).
    Simplexify,  ///< Cells -> simplices of the same topological dimension.
    Elevate,     ///< Linear cells -> serendipity quadratic (adds mid-edge nodes).
};

/**
 * @brief Parse a mode name into a `ConvertCellsMode`.
 * @param rName one of `"linearize"`, `"simplexify"`, `"elevate"`.
 * @return the matching mode.
 * @throws std::invalid_argument on an unknown name.
 */
MESHIOPLUSPLUS_API ConvertCellsMode convert_cells_mode_from_name(const std::string& rName);

/// Options for `convert_cells`.
struct ConvertCellsOptions {
    /// The conversion to perform.
    ConvertCellsMode mMode = ConvertCellsMode::Linearize;
    /// Attach an Int64 `convert:parent_cell` `cell_data` array recording, per
    /// output cell, the index of the input cell it came from *within its own
    /// block* (blocks correspond 1:1). Mainly useful under `Simplexify`, where
    /// several children share a parent.
    bool mRecordParentIds = false;
};

/// The result of `convert_cells`: the converted mesh plus the index maps.
struct ConvertCellsResult {
    /// The converted mesh.
    Mesh mMesh;
    /// Int64 shape `(num_points_in,)`, input point index -> output point index
    /// (-1 when the point was pruned, which only happens under `Linearize`).
    NDArray mPointMap;
    /// Per input block, Int64 shape `(num_cells_in_block,)`, input cell -> the
    /// index of its **first** child in the corresponding output block. Under
    /// `Linearize`/`Elevate` there is exactly one child per parent, so this is
    /// the identity; under `Simplexify` a parent's children occupy the
    /// contiguous range starting here.
    std::vector<NDArray> mCellMaps;
};

/**
 * @brief Convert the element representation of a mesh.
 * @param rMesh the mesh to convert.
 * @param rOptions the mode and flags (defaults to `Linearize`, no parent ids).
 * @return the converted mesh plus the point/cell index maps.
 * @throws std::invalid_argument when a block cannot be converted in the
 *         requested mode (a full-Lagrange target such as `quad9`/`hexahedron27`
 *         under `Elevate`, or a polyhedron cell under `Simplexify` whose faces
 *         are not a closed orientable surface).
 */
MESHIOPLUSPLUS_API ConvertCellsResult convert_cells(const Mesh& rMesh,
                                                    const ConvertCellsOptions& rOptions = {});

}  // namespace meshioplusplus

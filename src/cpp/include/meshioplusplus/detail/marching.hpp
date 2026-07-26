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
 * @file marching.hpp
 * @brief The marching-tetrahedra cutter shared by `operations/slice.cpp` (the
 * planar cross-section) and `operations/isosurface.cpp` (the level set of a
 * scalar field) — hoisted verbatim from slice in v7.15.0, the way
 * `spatial_hash.hpp` was hoisted from merge in v7.13.0 and `space_filling.hpp`
 * from reorder in v7.6.0. Slice's output is byte-identical across the hoist.
 *
 * Nothing about marching tetrahedra is specific to a *plane*: the algorithm
 * cuts wherever a **per-node signed scalar** changes sign. The two consumers
 * differ only in how they compute that scalar — slice uses
 * `d_i = dot(x_i - origin, normal)`, isosurface `d_i = f_i - isovalue` — and in
 * the name of the provenance array they ask for.
 *
 * The cutter is deliberately **two-phase**:
 *
 *  - `marching_prepare(mesh)` does everything independent of the cut scalar:
 *    `convert_cells(Simplexify)` (so every 3D cell is a tetrahedron and every
 *    2D cell a triangle — each cell's section is then a convex, unambiguously
 *    ordered primitive), the per-block input cell bases, and the mesh's maximum
 *    topological dimension. The caller resolves its scalar against the
 *    **simplexified** mesh (`Simp()`), which is what lets isosurface contour a
 *    `point_data` array that rode through the linearization.
 *  - `marching_cut(input, values, options)` performs one cut and returns the
 *    section mesh. Isosurface calls it once per isovalue against the *same*
 *    prepared input, so several contours cost one simplexification.
 *
 * Contracts carried over from slice verbatim, and relied upon by both:
 *
 *  - **Crossing** — an edge `(i,j)` is cut iff its endpoint values straddle
 *    zero; the crossing is `x_lo + t*(x_hi - x_lo)` with
 *    `t = d_lo/(d_lo - d_hi)` computed from the **sorted** endpoint pair, so
 *    two simplices sharing the edge agree bit for bit.
 *  - **Watertight** — crossings are keyed by that sorted simplexified-node pair
 *    and deduped, so a shared edge yields exactly one output node.
 *  - **Degeneracy** — a node whose value is exactly zero counts as *positive*
 *    (`d >= 0`), which makes the 4-bit (tetra) / 3-bit (triangle) sign mask
 *    total: a cut grazing a shared face is emitted by exactly one incident
 *    simplex, never two. Primitives whose crossings collapse (area/length below
 *    a bbox-relative tolerance) are dropped.
 *  - **Determinism** — the cut and the node numbering are a single **serial**
 *    first-seen sweep in a fixed (block, cell, ring-edge) order over a fixed
 *    case table; only the coordinate/`point_data` passes are parallel. Output
 *    is byte-identical across the three mesh backends, thread counts and the
 *    C++/numpy boundary (`_marching.py` is the expression-for-expression twin).
 *
 * The one genuinely parameterized rule is the **winding** (see
 * `MarchingOrientation`): slice orients every section face toward a fixed user
 * direction, while an isosurface has no global direction and orients toward
 * increasing field instead.
 *
 * `detail/` header: exempt from the operations layer's anon-namespace prefix
 * rule.
 */

// System includes
#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/operations/convert_cells.hpp"

namespace meshioplusplus {
namespace detail {

/**
 * @brief Everything the cutter needs that does not depend on the cut scalar.
 *
 * Owns the simplexified mesh; move it, never copy it (the KRATOS `Mesh` is not
 * copy-constructible).
 */
struct MarchingInput {
    /// The simplexified mesh plus its index maps (`convert:parent_cell` is
    /// recorded on `mMesh`'s `cell_data`).
    ConvertCellsResult mConvert;

    /// Per simplexified block (which maps 1:1 to an input block), the global
    /// block-major index of that input block's first cell.
    std::vector<std::int64_t> mOrigBlockBase;

    /// The maximum topological dimension present: 3 cuts tetrahedra into
    /// triangles/quads, 2 cuts triangles into line segments, anything else
    /// yields nothing.
    int mMaxDim = 0;

    /// The simplexified mesh the caller must resolve its scalar against.
    const Mesh& Simp() const { return mConvert.mMesh; }
};

/// How a section primitive's winding is resolved.
enum class MarchingOrientation {
    /// Flip each face so its Newell normal points along `MarchingOptions`'
    /// fixed `mDirection` — slice's rule, against the user's plane normal.
    FixedDirection,

    /// Flip each face so its Newell normal points toward **increasing field**:
    /// the per-simplex reference direction is the centroid of the corners with
    /// `d >= 0` minus the centroid of those with `d < 0` (both are non-empty
    /// whenever the simplex is cut at all). No matrix solve, and exactly
    /// reproducible in numpy — isosurface's rule.
    FieldGradient
};

/// Options for one cut.
struct MarchingOptions {
    /// Name of the Int64 `cell_data` array recording, per section cell, the
    /// global (block-major) index of the **input** cell it was cut from. Only
    /// written when `mRecordParentIds` is set.
    std::string mProvenanceName;

    /// Whether to attach `mProvenanceName`.
    bool mRecordParentIds = false;

    /// How to wind the section faces.
    MarchingOrientation mOrientation = MarchingOrientation::FixedDirection;

    /// The reference direction for `MarchingOrientation::FixedDirection`
    /// (ignored otherwise; need not be unit length).
    std::array<double, 3> mDirection{{0.0, 0.0, 1.0}};
};

/**
 * @brief Simplexifies @p rMesh and collects everything independent of the cut.
 * @param rMesh the mesh to cut (never modified).
 * @return the prepared input; `Simp()` is what a caller resolves its per-node
 *         scalar against, and `Simp().NumPoints()` is the length
 *         `marching_cut` expects of that scalar.
 * @throws std::invalid_argument when a block cannot be simplexified (a
 *         polyhedron block).
 */
MESHIOPLUSPLUS_API MarchingInput marching_prepare(const Mesh& rMesh);

/**
 * @brief Cuts the prepared mesh where @p rNodeValues changes sign.
 * @param rInput the prepared (simplexified) input.
 * @param rNodeValues one signed value per point of `rInput.Simp()`; the cut is
 *        its zero level set, with zero itself counting as positive.
 * @param rOptions the winding rule and the provenance array name.
 * @return a new mesh one topological dimension below the cut cells:
 *         `triangle`/`quad` blocks (in that fixed order) for a volume mesh, a
 *         `line` block for a 2D surface mesh, or an empty mesh when nothing is
 *         cut. Every `point_data` array of `rInput.Simp()` is interpolated at
 *         the crossings (promoted to Float64) and every per-cell `cell_data`
 *         array is replicated from the parent cell; `convert:parent_cell` is
 *         consumed rather than forwarded.
 */
MESHIOPLUSPLUS_API Mesh marching_cut(const MarchingInput& rInput, const std::vector<double>& rNodeValues,
                  const MarchingOptions& rOptions);

}  // namespace detail
}  // namespace meshioplusplus

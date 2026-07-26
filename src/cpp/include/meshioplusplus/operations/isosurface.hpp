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
 * @file isosurface.hpp
 * @brief Isosurfaces / contours: the level set of a scalar field, one
 * topological dimension below the cut cells. This is an operation, not a file
 * format — deliberately not in the registry — built entirely through the
 * uniform mesh API so it compiles under every mesh backend.
 *
 * `isosurface(mesh, options)` returns a **new mesh** that is the locus where a
 * `point_data` array equals each requested isovalue:
 *
 * - a **3D volume** mesh yields a **surface** mesh (`triangle`/`quad` faces);
 * - a **2D surface** mesh yields a **line** mesh (contour segments).
 *
 * This is the **data-driven sibling of `slice`**: slice cuts where
 * `dot(x - origin, normal) = 0`, isosurface where `f(x) - isovalue = 0`. They
 * share one cutter — `detail/marching.hpp` — so the watertightness, winding,
 * degeneracy and determinism guarantees below are literally the same code.
 *
 * **The field must be `point_data`.** `cell_data` is piecewise constant, so
 * there is no crossing to locate inside a cell and no level set to draw;
 * naming a `cell_data` array throws, pointing at `cell_data_to_point_data`
 * (CLI: `meshioplusplus data to-point`) as the fix. A multi-component array
 * reduces to `mComponent`, or to the row magnitude when that is unset.
 *
 * **Several isovalues land in one mesh.** They are cut in **ascending** order
 * (sorted, exact duplicates dropped) and concatenated, with the section blocks
 * merged by cell type in the canonical `triangle`, `quad`, `line` order, so a
 * single-isovalue call has exactly slice's block structure. Two `cell_data`
 * tags identify each contour: Float64 **`iso:value`** (the isovalue) and Int64
 * **`iso:index`** (its ordinal in the ascending list). Splitting into one mesh
 * per contour is the caller's job via `split(SplitBy::Tag, "iso:index")` (CLI
 * `split --by region --tag iso:index`) — the tag criterion needs an integer
 * array, which is what `iso:index` is for.
 *
 * **The contoured field is exactly the isovalue on the cut points.** Every
 * other `point_data` array is interpolated at the crossing parameter
 * `t = d_lo/(d_lo - d_hi)` (promoted to Float64, exact for a linear field); the
 * contoured array's own value would only reach the isovalue to within round-off,
 * so it is written exactly instead. The one exception is a multi-component array
 * reduced by **magnitude**: `|lerp(v)| != lerp(|v|)` mathematically, not merely
 * in floating point, so that case stays approximate and is left interpolated.
 *
 * Degeneracy rule (inherited, uniform with slice): a node whose value is exactly
 * the isovalue counts as being on the **positive** side (`d >= 0`), which makes
 * the sign mask total — a plateau lying exactly at the isovalue therefore emits
 * its boundary **once**, not twice from both incident cells. Primitives whose
 * crossings collapse (area/length below a bbox-relative tolerance) are dropped.
 * Section faces are wound so their Newell normal points toward **increasing
 * field**.
 *
 * An isovalue outside the field's range yields an empty contour for that value
 * rather than an error; all values outside yields an empty mesh.
 *
 * Each section cell inherits its parent cell's `cell_data` row (replicated, like
 * simplexify); with `mRecordParentIds` an Int64 `iso:parent_cell` `cell_data`
 * array records the originating **input** cell (global, block-major over the
 * input mesh). Section points are all new, so `point_sets`/`cell_sets`,
 * `mesh.info` and `gmsh_periodic` are **not** carried.
 *
 * Output is byte-identical across the three mesh backends, thread counts and the
 * C++/numpy boundary (`tests/python/test_isosurface.py::test_cpp_matches_python`).
 */

// System includes
#include <optional>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// Options for `isosurface`.
struct IsosurfaceOptions {
    /// Name of the `point_data` array to contour (required).
    std::string mArrayName;

    /// The level values. At least one; sorted ascending and de-duplicated
    /// before cutting, so contours are emitted in ascending order.
    std::vector<double> mIsovalues;

    /// Component of a multi-component array to contour; the row magnitude when
    /// unset (in which case the exactly-the-isovalue guarantee does not apply).
    std::optional<int> mComponent;

    /// Attach an Int64 `iso:parent_cell` `cell_data` array recording, per
    /// section cell, the global (block-major) index of the input cell it was
    /// cut from.
    bool mRecordParentIds = false;
};

/**
 * @brief Computes the isosurfaces (level sets) of a scalar `point_data` field.
 * @param rMesh the mesh to contour (never modified).
 * @param rOptions the array name, the isovalues, the component and the
 *        provenance flag.
 * @return a new mesh one topological dimension below the cut cells: a
 *         `triangle`/`quad` surface for a volume mesh, a `line` mesh for a 2D
 *         surface mesh, or an empty mesh when no isovalue crosses the field.
 *         Every contour cell carries a Float64 `iso:value` and an Int64
 *         `iso:index` `cell_data` entry.
 * @throws std::invalid_argument when the array name is empty, names a
 *         `cell_data` array (contour a point field instead — see
 *         `cell_data_to_point_data`) or no array at all, when `mIsovalues` is
 *         empty or holds a non-finite value, when `mComponent` is out of range,
 *         or when a cell block cannot be simplexified (a polyhedron block).
 */
MESHIOPLUSPLUS_API Mesh isosurface(const Mesh& rMesh, const IsosurfaceOptions& rOptions);

}  // namespace meshioplusplus

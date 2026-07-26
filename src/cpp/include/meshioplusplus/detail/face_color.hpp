/**
 * @file face_color.hpp
 * @brief Resolve a named data array into one color per drawn face.
 *
 * This is the shared middle of the SVG and TikZ writers' data-driven coloring.
 * There are four emission sites for that feature -- `src/cpp/src/formats/svg.cpp`,
 * `src/cpp/src/formats/tikz.cpp` and their pure-Python twins -- and everything
 * except the format's own color vocabulary lives here, so the logic exists
 * exactly twice (once per language) rather than four times. The Python twin is
 * `src/python/meshioplusplus/_facecolor.py`, and the two must stay in step: the two
 * writers are pinned byte-for-byte against each other by `tests/python/test_svg.py`
 * and `tests/python/test_tikz.py`, so every arithmetic step here has an
 * expression-for-expression counterpart there.
 *
 * The resolution rules, in one place:
 *  - The name is looked up in the *source* mesh's `point_data` first, then its
 *    `cell_data`; present in neither is an error naming what is available.
 *  - Point data colors a face by the mean of its corner values; cell data by
 *    the value of the cell that owns it.
 *  - A face's owning cell is found through `"surface:parent_cell"` when the
 *    drawn mesh is a boundary extracted from the source (see
 *    `detail::surface_extract`), and directly otherwise.
 *  - Multi-component arrays reduce to a component (`mComponent`) or to their
 *    magnitude.
 *  - Non-finite values are excluded from the auto range and are drawn in the
 *    writer's `nan_color` rather than mapped through the colormap.
 */
#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "meshioplusplus/export.hpp"
#include "meshioplusplus/detail/colormap.hpp"
#include "meshioplusplus/detail/projection.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {
namespace detail {

/** @brief What the writer was asked to color by. */
struct ColorSpec {
    std::string mColorBy;           ///< Data array name; empty disables coloring.
    std::optional<int> mComponent;  ///< Component to take; magnitude when unset.
    std::string mCmap = "viridis";  ///< Built-in colormap name.
    std::optional<double> mVMin;    ///< Range low end; auto when unset.
    std::optional<double> mVMax;    ///< Range high end; auto when unset.
    bool mColorbar = false;         ///< Whether to draw the colorbar.

    /** @brief Whether coloring was requested at all. */
    bool Active() const { return !mColorBy.empty(); }
};

/**
 * @brief One drawn face, as the color resolver needs to see it.
 *
 * Deliberately the same shape as the drawing side minus the geometry, so the
 * projected path (`ProjectedFace`) and the flat path can share one resolver.
 */
struct ColorFace {
    std::array<std::int64_t, 4> mNodes;  ///< Corner node ids into the drawn mesh.
    std::uint8_t mNumNodes;              ///< Corners actually used (2, 3 or 4).
    std::int64_t mSourceCell;            ///< Cell index in the drawn mesh, block-major.
};

/** @brief Per-face scalars and the resolved range, ready for lookup. */
struct FaceColors {
    bool mActive = false;                   ///< False when no coloring was requested.
    std::vector<double> mValues;            ///< One per face; non-finite => nan_color.
    double mVMin = 0.0;                     ///< Low end of the mapped range.
    double mVMax = 1.0;                     ///< High end of the mapped range.
    const std::uint8_t* mpTable = nullptr;  ///< The resolved colormap table.

    /**
     * @brief The color for one face.
     * @param Face index into `mValues`
     * @return the mapped color, or `std::nullopt` when the face's value is
     *         non-finite and the caller should use its `nan_color`
     */
    std::optional<Rgb> Color(std::size_t Face) const;
};

/**
 * @brief Normalize a value into [0, 1] over the range, clamping.
 *
 * A degenerate range (`vMin == vMax`, which an all-constant array produces)
 * maps everything to the middle of the colormap rather than dividing by zero.
 */
MESHIOPLUSPLUS_API double color_param(double v, double vMin, double vMax);

/** @brief The faces a projected surface draws, in emission order. */
MESHIOPLUSPLUS_API std::vector<ColorFace> color_faces_from_projection(const std::vector<ProjectedFace>& rFaces);

/**
 * @brief The faces the flat 2D path draws, in its emission order.
 *
 * Enumerates `line`/`triangle`/`quad` blocks block-major then cell-major --
 * the same rule, and the same order, both flat writers loop over.
 */
MESHIOPLUSPLUS_API std::vector<ColorFace> color_faces_flat(const Mesh& rMesh);

/**
 * @brief Resolve `rSpec` into one scalar per face plus the mapped range.
 *
 * @param rSpec what to color by; an inactive spec yields an inactive result
 * @param rSource the writer's input mesh -- the one the array name refers to
 * @param rDraw the mesh actually being drawn; equal to `rSource` unless a
 *        boundary was extracted, in which case its `"surface:parent_cell"`
 *        carries the provenance back to `rSource`'s cells
 * @param rFaces the drawn faces, in emission order
 * @return per-face values and range; `mValues` is parallel to `rFaces`
 * @throws std::invalid_argument for an unknown array name, an unknown
 *         colormap, an out-of-range component, or `vmin > vmax`
 */
MESHIOPLUSPLUS_API FaceColors resolve_face_colors(const ColorSpec& rSpec, const Mesh& rSource, const Mesh& rDraw,
                               const std::vector<ColorFace>& rFaces);

}  // namespace detail
}  // namespace meshioplusplus

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
 * @file projection.hpp
 * @brief Orthographic camera projection + painter's-algorithm face ordering
 * for the SVG/TikZ 3D rendering path.
 *
 * The camera is parameterized by azimuth/elevation/roll in degrees: the view
 * direction (from the scene toward the camera) is
 * `w = (cos el * cos az, cos el * sin az, sin el)`, the screen-right axis is
 * `u = normalize(cross((0,0,1), w))` (falling back to the y axis as the up
 * reference when looking straight along z), the screen-up axis is
 * `v = cross(w, u)`, and `roll` rotates `(u, v)` about `w`. The defaults
 * `azimuth = 45`, `elevation = 35.264389682754654` (`atan(1/sqrt(2))`) give
 * `w` proportional to `(1,1,1)` — the classic CAD isometric view. Projection
 * is orthographic (`x' = p.u`, `y' = p.v`), and `depth = p.w` orders faces
 * back-to-front (painter's algorithm; larger depth = closer to the camera).
 *
 * KEEP IN SYNC: `src/python/meshioplusplus/_projection.py` is the Python twin used
 * by the pure-Python fallback writers — the arithmetic (down to expression
 * order, which fixes the floating-point rounding) must stay identical so
 * TikZ output is byte-identical across the two implementations.
 *
 * `camera_basis` is called once per write, and `project_surface` once per
 * mesh (its own per-point/per-face loops are unaffected by where the
 * function itself is compiled), so both bodies live in
 * `src/cpp/src/detail/projection.cpp` rather than inline here.
 */

// System includes
#include <array>
#include <cstdint>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {
namespace detail {

/** @brief Orthonormal camera basis: screen-right `u`, screen-up `v`, and the
 * view direction `w` (scene toward camera). */
struct CameraBasis {
    std::array<double, 3> mU;
    std::array<double, 3> mV;
    std::array<double, 3> mW;
};

/**
 * @brief Build the orthographic camera basis for azimuth/elevation/roll.
 * @param azimuth Camera azimuth in degrees (rotation about z, from +x).
 * @param elevation Camera elevation in degrees above the xy plane.
 * @param roll In-screen rotation in degrees about the view direction.
 * @return The `{u, v, w}` basis described in the file header.
 */
MESHIOPLUSPLUS_API CameraBasis camera_basis(double azimuth, double elevation, double roll);

/** @brief One drawable face of a projected surface: 2-4 corner point ids. */
struct ProjectedFace {
    std::array<std::int64_t, 4> mNodes;
    std::uint8_t mNumNodes;  // 2 (line), 3 (triangle), or 4 (quad)
    bool mIsLine;
    double mDepth;  // view-space depth of the face centroid

    // Index of the cell this face came from, counted block-major over every
    // cell of every block of the projected mesh -- including blocks this
    // projection skips, matching "surface:parent_cell"'s convention.
    //
    // It rides on the face rather than in a parallel array because the faces
    // are stable_sorted below: a side array would have to be permuted in
    // lockstep, whereas carrying the id inside the sorted element makes the
    // correspondence impossible to break.
    std::int64_t mSourceCell;
};

/** @brief A surface mesh projected to screen space, faces sorted
 * back-to-front. */
struct ProjectedSurface {
    std::vector<double> mX;
    std::vector<double> mY;
    std::vector<ProjectedFace> mFaces;
};

/**
 * @brief Project a surface mesh's points with an orthographic camera and
 * gather its drawable faces sorted back-to-front (painter's algorithm).
 *
 * Drawable blocks are `line`, `triangle`, `quad`, and the corner-linearized
 * higher-order surface types (`triangle6`, `quad8`, `quad9` — corner nodes
 * only); anything else is skipped. Faces keep enumeration order within
 * equal depths (`std::stable_sort`, matching numpy's stable argsort).
 *
 * @param rMesh The surface mesh to project (already skin-extracted when it
 *              came from a volume mesh).
 * @param azimuth Camera azimuth in degrees.
 * @param elevation Camera elevation in degrees.
 * @param roll In-screen rotation in degrees.
 * @return Projected screen coordinates per point and the sorted face list.
 */
MESHIOPLUSPLUS_API ProjectedSurface project_surface(const Mesh& rMesh, double azimuth, double elevation, double roll);

}  // namespace detail
}  // namespace meshioplusplus

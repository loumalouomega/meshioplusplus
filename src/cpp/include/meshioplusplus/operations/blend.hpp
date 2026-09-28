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
 * @file operations/blend.hpp
 * @brief Linear interpolation between two steps of one mesh -- the kernel time
 * resampling of a sequence is built on.
 *
 * `blend_steps(a, b, w)` returns `a` with every floating-point data array
 * replaced by `(1 - w) * a + w * b`: point data, cell data and field data.
 * The two steps must share a topology -- the same point count, the same cell
 * blocks (type and count) and the same arrays with the same shapes -- which is
 * what a solver's output series has; a mismatch throws naming it rather than
 * blending unrelated rows. Points, connectivity and regions come from `a`,
 * bit for bit, unless `mBlendPoints` also moves the points (a series on a
 * moving mesh). An integer array (a material id, a flag) is not blended: it is
 * taken from the nearer step (`a` when `w < 0.5`). A Float32 array stays
 * Float32; the blend is computed in double.
 *
 * `w` outside `[0, 1]` extrapolates linearly; the sequence resampler never
 * passes one (see `SequenceResample`).
 */

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// How `blend_steps` blends.
struct BlendOptions {
    /// Also blend the point coordinates (a moving mesh).
    bool mBlendPoints = false;
};

/**
 * @brief `a` with its floating-point data linearly blended toward `b` by `w`.
 * @throws std::invalid_argument when the two steps differ in point count, cell
 *         blocks, or the names or shapes of their data arrays.
 */
MESHIOPLUSPLUS_API Mesh blend_steps(const Mesh& rA, const Mesh& rB, double W,
                                    const BlendOptions& rOptions = {});

}  // namespace meshioplusplus

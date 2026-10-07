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
 * @file detail/render_field.hpp
 * @brief The pure helpers of field rendering: value scales, percentiles,
 * legend ticks, the categorical palette and marching on one triangle.
 *
 * A **core-private** header (the `crease_edges.hpp` precedent): no installed
 * header names it, and it adds nothing to the API or the ABI. Everything here
 * is a function of its arguments, so it is tested without a mesh.
 */

// System includes
#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

// Project includes
#include "meshioplusplus/operations/render.hpp"

namespace meshioplusplus {
namespace detail {

/**
 * @brief A value on a colormap's scale: `v` itself, its base-10 logarithm, or a
 * signed logarithm that is linear within `Threshold` of zero.
 * @return the scaled value, or NaN where the scale is undefined (a log scale at
 *         a value that is not positive, or any scale at NaN)
 */
double render_scale_forward(RenderScale Scale, double Threshold, double v);

/**
 * @brief The `P`-th percentile (in [0, 100]) of values, by linear interpolation
 * between the closest ranks. @p rValues is sorted in place; non-finite values
 * must already be removed. NaN for an empty list.
 */
double render_percentile(std::vector<double>& rValues, double P);

/**
 * @brief Tick values inside `[VMin, VMax]` for a legend: round numbers (1, 2 or
 * 5 times a power of ten) about `Target` of them on a linear scale, whole decades
 * on a log scale, and `0` plus decades either side on a symlog scale.
 */
std::vector<double> render_legend_ticks(RenderScale Scale, double Threshold, double VMin,
                                        double VMax, int Target);

/// The qualitative palette of categorical rendering: the Okabe-Ito colours
/// (Okabe and Ito, 2002, designed to survive colour-vision deficiency), then
/// two greys, cycling past ten categories.
const std::array<std::array<std::uint8_t, 4>, 10>& render_category_palette();

/**
 * @brief One contour segment of the linear interpolant on a triangle.
 * @param pXyz the three corners, 3 doubles each
 * @param pValue the three corner values
 * @param Level the contour level
 * @param pOut the two end points, 6 doubles
 * @return whether the level crosses the triangle. A corner exactly on the level
 *         counts as above it, so two triangles sharing an edge agree and a level
 *         through a vertex never draws a stray point.
 */
bool render_triangle_contour(const double* pXyz, const double* pValue, double Level, double* pOut);

}  // namespace detail
}  // namespace meshioplusplus

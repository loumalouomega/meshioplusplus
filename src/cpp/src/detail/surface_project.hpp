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
 * @file detail/surface_project.hpp
 * @brief The nearest point of a NURBS surface to a point (roadmap 6.10.1):
 * point inversion.
 *
 * Newton on the first-order condition `(S - p) . S_u = (S - p) . S_v = 0`,
 * whose Hessian is `J^T J + (S - p) . S_2` with the surface's second
 * derivatives; where that is not positive definite the step is Gauss-
 * Newton's (the first fundamental form), and it is halved until the
 * distance does not grow; the parameters stay in the knot range. Newton
 * converges only from inside the right basin, so it starts from the nearest
 * few of a grid of the patch's points, `(degree + 1)` per knot span in each
 * direction (8 to 32), the subdivision start of Ma and Hewitt (2003,
 * *Point inversion and projection for NURBS curve and surface: control
 * polygon approach*), and the nearest result is kept. The same as
 * `detail::bspline_project` on a `BsplineSurface`, on a rational patch.
 *
 * A **core-private** header: nothing installed names it.
 */

// System includes
#include <array>
#include <vector>

// Project includes (private, not installed)
#include "edge_sampler.hpp"

namespace meshioplusplus {

struct NurbsPatch;

namespace detail {

class SurfaceProjector {
public:
    /// `rSurface` and `rPatch` (the patch it evaluates) outlive it.
    SurfaceProjector(const BrepSurfaceEval& rSurface, const NurbsPatch& rPatch);
    /// The parameters of the surface point nearest `rP`, and in `rX` that
    /// point. Safe to call from many threads.
    std::array<double, 2> Project(const EdgeVec3& rP, EdgeVec3& rX) const;
    /// Newton from `rUv`: the distance reached, `rUv` updated.
    double Polish(const EdgeVec3& rP, std::array<double, 2>& rUv) const;

private:
    const BrepSurfaceEval* pSurface = nullptr;
    std::vector<std::array<double, 2>> mSeedUv;
    std::vector<EdgeVec3> mSeedX;
};

}  // namespace detail
}  // namespace meshioplusplus

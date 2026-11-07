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

#pragma once

/**
 * @file detail/bounding_box_utils.hpp
 * @brief `BoundingBoxUtils`: axis-aligned boxes of point sets, the loop that
 * operations used to write out by hand (doc/geometry_kernels.md, roadmap §8).
 *
 * A Kratos-style utility class (cf. `Kratos::BoundingBox`): a plain
 * `BoundingBox3` value and static inline members. `FromPoints` starts from the
 * first point and widens with `std::min(lo, q)` / `std::max(hi, q)`, the form
 * the copies it replaced used: it keeps `lo` on a tie or a NaN coordinate, so
 * it matches `if (q < lo) lo = q` too, but not a loop seeded with +-infinity
 * (which differs on NaN and on an empty set) -- those stay as they are.
 *
 * A **core-private** header (not installed): no API or ABI surface.
 */

// System includes
#include <algorithm>
#include <vector>

// Project includes
#include "meshioplusplus/detail/geometry.hpp"

namespace meshioplusplus::detail {

/// An axis-aligned box, `mMin` and `mMax` corners.
struct BoundingBox3 {
    Vec3 mMin{0.0, 0.0, 0.0};
    Vec3 mMax{0.0, 0.0, 0.0};
};

/// Static axis-aligned bounding-box helpers.
struct BoundingBoxUtils {
    /// Widens `rBox` to hold `rQ`: `std::min(lo, q)`, `std::max(hi, q)` per axis.
    static void Expand(BoundingBox3& rBox, const Vec3& rQ) {
        for (std::size_t a = 0; a < 3; ++a) {
            rBox.mMin[a] = std::min(rBox.mMin[a], rQ[a]);
            rBox.mMax[a] = std::max(rBox.mMax[a], rQ[a]);
        }
    }

    /**
     * @brief The box of `rPoints`, seeded with the first point. An empty set
     * gives the degenerate box at the origin.
     */
    static BoundingBox3 FromPoints(const std::vector<Vec3>& rPoints) {
        BoundingBox3 box;
        if (rPoints.empty())
            return box;
        box.mMin = rPoints[0];
        box.mMax = rPoints[0];
        for (const Vec3& r_q : rPoints)
            Expand(box, r_q);
        return box;
    }

    /**
     * @brief The box of `N >= 1` coordinates supplied by `rCoord(i) -> Vec3`,
     * written to the flat triples `pLo` and `pHi`: seeded with coordinate 0,
     * then widened by `if (c < lo) lo = c` and `if (c > hi) hi = c` (the same
     * as `Expand`'s `std::min`/`std::max`). For callers that gather coordinates
     * on the fly (`interpolate`'s point, centroid and per-simplex boxes).
     */
    template <class TCoord>
    static void FromCoords(std::size_t N, TCoord&& rCoord, double* pLo, double* pHi) {
        const Vec3 first = rCoord(0);
        for (int d = 0; d < 3; ++d)
            pLo[d] = pHi[d] = first[d];
        for (std::size_t i = 1; i < N; ++i) {
            const Vec3 c = rCoord(i);
            for (int d = 0; d < 3; ++d) {
                if (c[d] < pLo[d])
                    pLo[d] = c[d];
                if (c[d] > pHi[d])
                    pHi[d] = c[d];
            }
        }
    }

    /// `max - min` per axis.
    static Vec3 Extent(const BoundingBox3& rBox) { return vec3_sub(rBox.mMax, rBox.mMin); }

    /// The largest of the three extents.
    static double MaxExtent(const BoundingBox3& rBox) {
        const Vec3 e = Extent(rBox);
        return std::max({e[0], e[1], e[2]});
    }

    /// Half the box's surface area, `ex ey + ey ez + ez ex`.
    static double HalfSurfaceArea(const BoundingBox3& rBox) {
        const Vec3 e = Extent(rBox);
        return e[0] * e[1] + e[1] * e[2] + e[2] * e[0];
    }

    /// The diagonal's length, `|max - min|`.
    static double Diagonal(const BoundingBox3& rBox) { return vec3_norm(Extent(rBox)); }
};

}  // namespace meshioplusplus::detail

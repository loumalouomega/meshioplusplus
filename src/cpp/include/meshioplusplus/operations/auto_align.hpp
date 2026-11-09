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
 * @file operations/auto_align.hpp
 * @brief A consistent pose from a scan alone (roadmap 7.3.5): the rotation
 * and translation that lay a part on its largest flat face, or stand a
 * cylindrical part on the end of its axis, without a nominal to register to.
 *
 * The facet normals of the input's surface cells, weighted by area, are
 * binned on the Gauss sphere: a fixed subdivided icosahedron
 * (`mSubdivisions` 1-to-4 splits; 642 bins at the default 3), each facet in
 * the bin of its nearest vertex. A bin scores the area of every bin within
 * `mAngleToleranceDeg` of it; the best (ties to the lowest bin) gives the
 * dominant direction, refined to the area-weighted mean of the facet normals
 * within the tolerance, and the best bin at least `90 - tolerance` degrees
 * from it the second direction (Gram--Schmidt against the first).
 *
 * An open scan's opening is told from a face by the vector area of its
 * boundary: the facets' summed vector area is the opening's (by Stokes'
 * theorem, the vector area of the boundary loops), so when it exceeds 0.1%
 * of the total area, a pseudo-facet of that area, facing away from the
 * facets' sum, joins the histogram: a part scanned on a turntable lies on
 * its opening. The winding must be consistent for this, as for the faces.
 *
 * Two poses are scored by the share of the area whose normals agree with
 * them within the tolerance: `face` (along any of the three axes of the
 * frame above) and `axis` (the eigenvector of the least eigenvalue of the
 * area-weighted normal scatter, a cylinder's axis: normals perpendicular to
 * it or along it). `axis` is taken when it beats `face` by more than 0.05.
 * `face` maps the dominant direction to -Z (the face down) and the second to
 * +X; `axis` maps the axis to +Z, the end with the larger area facing along
 * the axis down, and the dominant bin around the axis to +X. The result is
 * then translated so its bounding box is centred on the Z axis with its
 * lowest point at z = 0.
 *
 * Only points and surface cells are read; data is carried, vector and
 * tensor point data rotated with the points when `mRotateVectorData` is set
 * (as `transform`). Deterministic: no sampling, every tie to the lowest
 * index.
 * @see doc/registration.md
 */

// System includes
#include <array>
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// Options for `auto_align`.
struct AutoAlignOptions {
    /// The angle within which a normal agrees with a direction, in degrees,
    /// in (0, 45).
    double mAngleToleranceDeg = 10.0;
    /// The icosphere's 1-to-4 splits (the bins: 12, 42, 162, 642, 2562, ...);
    /// 0 to 5.
    std::int64_t mSubdivisions = 3;
    /// Only the surface cells of this cell region vote; empty uses all.
    std::string mRegion;
    /// Rotate vector and tensor point data with the points.
    bool mRotateVectorData = false;
};

/// The aligned mesh and how it was posed.
struct AutoAlignResult {
    /// The input, moved by `mMatrix`.
    Mesh mMesh;
    /// Row-major 4x4: input coordinates onto the aligned ones.
    std::array<double, 16> mMatrix{1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
    /// `"face"` or `"axis"`.
    std::string mMode;
    /// In the input's frame: the direction sent to -Z (face) or +Z (axis),
    /// and the one sent to +X.
    std::array<double, 3> mPrimary{{0.0, 0.0, 0.0}};
    std::array<double, 3> mSecondary{{0.0, 0.0, 0.0}};
    /// The chosen pose's agreeing share of the area, and both poses'.
    double mShare = 0.0;
    double mFaceShare = 0.0;
    double mAxisShare = 0.0;
    /// The opening's area added to the histogram (0 for a closed surface).
    double mOpeningArea = 0.0;
    std::vector<std::string> mWarnings;
};

/**
 * @brief Pose @p rMesh on its largest flat face, or a cylinder on its axis
 * (see the file comment).
 * @throws std::invalid_argument on a bad option, an unknown region, a
 *         non-finite point, or no surface cell of positive area.
 */
MESHIOPLUSPLUS_API AutoAlignResult auto_align(const Mesh& rMesh,
                                              const AutoAlignOptions& rOptions = {});

}  // namespace meshioplusplus

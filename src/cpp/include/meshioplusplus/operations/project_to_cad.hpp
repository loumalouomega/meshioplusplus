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
 * @file project_to_cad.hpp
 * @brief Move a mesh's points back onto the CAD model it was made from
 * (roadmap 6.10.1).
 *
 * A mesh `mesh_cad_surface` wrote carries, per cell, the CAD face it lies on
 * (`cad:face` on its triangles and quads) and the CAD edge it runs along
 * (`cad:edge` on its lines). An operation that knows no CAD (smoothing, an
 * adaptation, a remesh that keeps those arrays) moves the points off the
 * surfaces; `project_to_cad` puts them back:
 *
 *  - a point on the lines of one edge goes to the nearest point of that
 *    edge: its own 3-D curve (`BrepEdge::mCurve`) where the model gives one
 *    that its first use's trace follows (within a millionth of its length),
 *    else that trace, its trim curve lifted through its face's surface (a
 *    scan of 64 steps, then Newton on `(X - p) . X' = 0`); a curve the
 *    trace leaves is warned;
 *  - a point on the lines of two or more edges is a corner of the model:
 *    it goes to the nearest point of the lowest of them, which is the edge's
 *    end for a point near it, and is counted;
 *  - a point whose triangles and quads all name one face goes to the nearest
 *    point of that face's surface: Newton on the first-order condition
 *    `(S - p) . S_u = (S - p) . S_v = 0` with its second derivatives,
 *    Gauss-Newton where that Hessian is not positive, the step halved until
 *    the distance does not grow, started from the nearest few of a grid of
 *    the patch's points (the subdivision start of Ma and Hewitt 2003,
 *    *Point inversion and projection for NURBS curve and surface*);
 *  - a point on two faces and no edge is ambiguous: it stays, counted.
 *
 * A projection that lands outside the face's trimmed region (by the trim
 * loops sampled `mTrimSamples` times per span, even-odd) is still the
 * nearest surface point and is kept, but reported, never clamped. Points no
 * face or edge cell uses stay where they are. Data, regions and property
 * sets pass through; only coordinates change. Points are projected
 * independently and in parallel, so the result does not depend on threads.
 */

// System includes
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/brep.hpp"
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

struct ProjectToCadOptions {
    /// The cell data naming each cell's face (on triangles and quads) and
    /// edge (on lines), as indices into the model's faces and edges; -1 is
    /// none.
    std::string mFaceArray = "cad:face";
    std::string mEdgeArray = "cad:edge";
    /// A point farther than this from where it would go stays where it is
    /// and is counted (0 is no limit).
    double mMaxDistance = 0.0;
    /// Samples per knot span when testing a face's trimmed region.
    std::int32_t mTrimSamples = 16;
};

struct ProjectToCadResult {
    Mesh mMesh;
    std::int64_t mNumProjected = 0;    ///< points moved onto a face or an edge
    std::int64_t mNumOnEdges = 0;      ///< of those, onto an edge
    std::int64_t mNumCorners = 0;      ///< of those, on two or more edges
    std::int64_t mNumAmbiguous = 0;    ///< on two faces and no edge: left
    std::int64_t mNumTooFar = 0;       ///< beyond `mMaxDistance`: left
    std::int64_t mNumFailed = 0;       ///< whose face or edge is unusable: left
    std::int64_t mNumOutsideTrim = 0;  ///< projected outside their face's trims
    double mMaxMove = 0.0;             ///< the largest distance a point moved
    std::int64_t mNumOnCurves = 0;     ///< of those on edges, onto an edge's 3-D curve
    /// The points projected outside their face's trims, ascending.
    std::vector<std::int64_t> mOutsideTrim;
    std::vector<std::string> mWarnings;
};

/**
 * @brief Move every point of `rMesh` onto the face or edge of `rBrep` its
 * cells name (see the file comment).
 * @throws std::invalid_argument when `rMesh` has no `mFaceArray` cell data,
 * its points are not three-dimensional, a face or edge index is out of
 * range, or an option is invalid.
 */
MESHIOPLUSPLUS_API ProjectToCadResult project_to_cad(const Mesh& rMesh, const BrepInfo& rBrep,
                                                     const ProjectToCadOptions& rOptions = {});

}  // namespace meshioplusplus

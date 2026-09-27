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
 * @file operations/feature_edges.hpp
 * @brief The sharp, open, non-manifold and inconsistently wound edges of a
 * surface, as a mesh of `line` cells.
 *
 * The surface is the mesh's 2-D cells, or -- when it has volume cells -- the
 * boundary skin of those. Each edge shared by two faces is compared by the
 * dihedral angle between their normals; an edge used once is open, three or
 * more times non-manifold, and a pair of faces that walk their shared edge the
 * same way disagree about which side is out. This is the crease test
 * `decimate`, `decimate_volume` and `smooth` pin nodes with (since v16.23.0),
 * exposed for inspection, for picking the edges a boundary condition lives on,
 * and for checking what those operations will hold still.
 *
 * ### Output
 *
 * The input's points, verbatim and not compacted (so a line's node ids are the
 * input's own), one `line` block in ascending `(low, high)` endpoint order, and
 * two cell-data arrays:
 *
 *  - `feature:kind` (Int32): 1 feature, 2 boundary, 3 non-manifold,
 *    4 inconsistent. An edge in several categories is labelled by the first
 *    selected one in the order non-manifold, boundary, inconsistent, feature.
 *  - `feature:angle` (Float64): the dihedral in degrees, measured after
 *    reorienting an inconsistent pair; NaN on edges that are not a pair of two
 *    non-degenerate faces.
 *
 * Point data, field data and Point regions ride through; Cell and Side regions
 * name cells that no longer exist and are dropped with a warning.
 */

// System includes
#include <cstdint>
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// Cell data (Int32): the category of each output edge; see `FeatureEdgeKind`.
inline constexpr const char* kFeatureKindName = "feature:kind";
/// Cell data (Float64): the dihedral angle of each output edge, in degrees.
inline constexpr const char* kFeatureAngleName = "feature:angle";

/// The values of `feature:kind`.
enum class FeatureEdgeKind : std::int32_t {
    Feature = 1,
    Boundary = 2,
    NonManifold = 3,
    Inconsistent = 4,
};

/// Which edges `feature_edges` reports.
struct FeatureEdgeOptions {
    /// The largest dihedral angle, in degrees, still treated as smooth; must
    /// lie in `[0, 180]`.
    double mFeatureAngleDeg = 30.0;
    /// Report edges whose dihedral exceeds `mFeatureAngleDeg`.
    bool mFeature = true;
    /// Report open edges (used by one face).
    bool mBoundary = true;
    /// Report edges used by three or more faces.
    bool mNonManifold = true;
    /// Report pairs of faces that walk their shared edge the same way.
    bool mInconsistent = true;
    /// Restrict to the cells of this named `Cell` region (surface cells, or the
    /// volume cells whose skin is taken); empty takes every cell.
    std::string mRegion;
};

/// The reported edges, and how many edges of each category the surface has.
struct FeatureEdgeResult {
    /// The input's points and one `line` block, with `feature:kind` and
    /// `feature:angle`.
    Mesh mMesh;
    /// Per category, the number of such edges on the surface, whether or not
    /// the category was selected for output.
    std::int64_t mNumFeature = 0;
    std::int64_t mNumBoundary = 0;
    std::int64_t mNumNonManifold = 0;
    std::int64_t mNumInconsistent = 0;
};

/**
 * @brief The feature edges of a surface (or of a volume mesh's skin).
 *
 * @param rMesh a surface mesh (triangles, quads, polygons and their quadratic
 *        variants, which contribute their corners), or a volume mesh; lines and
 *        vertices are ignored. With any volume cell present, only the skin of
 *        the volume cells is examined.
 * @param rOptions what to report; see `FeatureEdgeOptions`.
 * @throws std::invalid_argument on a feature angle outside `[0, 180]` or an
 *         unknown region name.
 */
MESHIOPLUSPLUS_API FeatureEdgeResult feature_edges(const Mesh& rMesh,
                                                   const FeatureEdgeOptions& rOptions = {});

}  // namespace meshioplusplus

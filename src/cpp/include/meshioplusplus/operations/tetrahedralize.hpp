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
 * @file tetrahedralize.hpp
 * @brief Tetrahedra filling a closed triangle surface, its triangles kept
 * (roadmap 6.7), through an optional backend.
 *
 * The input is a closed, manifold, consistently oriented triangle surface,
 * such as a `mesh_cad_surface` result or a boundary layer's cap: an open,
 * non-manifold, misoriented or degenerate one is refused by name, with the
 * counts of `surface_watertight_check`. Its connected shells are oriented
 * outward (a shell whose triangles face in is turned), and each shell bounds
 * a region: the inside of the shell minus the shells directly inside it
 * (found by the winding number of one of their points). A box inside a box
 * gives two regions, the shell between them and the inner box. Each region
 * is meshed on its own, from its shell and its inner shells turned inward,
 * so neighbouring regions share the inner shell's triangles exactly and the
 * tetrahedra conform across it.
 *
 * Two backends. Netgen's nglib (J. Schöberl, "NETGEN: an advancing front
 * 2D/3D-mesh generator based on abstract rules", Comput. Visual. Sci. 1,
 * 1997; LGPL-2.1, built with `MESHIOPLUSPLUS_WITH_NETGEN`, linked
 * dynamically, never vendored) keeps the surface's triangles and points
 * exactly, which is checked, and `mOptimizeSteps` passes of its volume
 * optimisation follow; calls into nglib are serialised (its state is
 * global). The native octree backend, always built, is `octree_mesh`'s
 * isosurface stuffing at the surface's own size (F. Labelle and J. R.
 * Shewchuk, "Isosurface stuffing", ACM TOG 26, 2007): its boundary is fitted
 * to the surface within the reported Hausdorff distance rather than kept,
 * its dihedral angles bounded where the stuffing's guarantee holds, and
 * two regions do not share points across their interface. `Auto` takes
 * Netgen when it is built in and the octree backend otherwise, with a
 * warning; asking for Netgen without it throws `Unsupported` naming the
 * flag. Each region is meshed in one thread, so the result is
 * deterministic.
 *
 * The result holds the points, then a `tetra` block (positively oriented,
 * region by region) and a `triangle` block: with Netgen the input's points
 * first (unchanged, at their indices) and the surface's triangles as given;
 * with the octree backend new points and the tetrahedra's boundary facets,
 * outward. The cell data `tetrahedralize:region` gives the region of each
 * tetrahedron (-1 on the triangles) and, when the surface carries
 * `mFaceArray` (`cad:face`), that array is on the triangles (-1 on the
 * tetrahedra; a boundary facet of the octree backend takes the value of the
 * surface triangle nearest its centroid); Cell regions `region:<k>` and, per
 * value of `mFaceArray`, Side regions `face:<v>` on the tetrahedra facets
 * the triangles are.
 */

// System includes
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

inline constexpr const char* kTetrahedralizeRegion = "tetrahedralize:region";

enum class TetrahedralizeBackend : std::int32_t {
    Auto = 0,    ///< Netgen when built in, else the octree backend
    Netgen = 1,  ///< Netgen's nglib (MESHIOPLUSPLUS_WITH_NETGEN): the surface kept
    Octree = 2,  ///< native isosurface stuffing (always built): the surface fitted
};

struct TetrahedralizeOptions {
    TetrahedralizeBackend mBackend = TetrahedralizeBackend::Auto;
    /// The largest edge asked of the tetrahedra (0: none; the surface's
    /// triangles then size the inside).
    double mMaxSize = 0.0;
    /// The largest volume asked of a tetrahedron (0: none): a regular
    /// tetrahedron's edge for it, `(6 sqrt 2 V)^(1/3)`, joins `mMaxSize`.
    double mMaxVolume = 0.0;
    /// How fast the size grows away from the surface, in [0, 1] (Netgen).
    double mGrading = 0.3;
    /// Passes of the backend's volume optimisation (0 to 100).
    std::int32_t mOptimizeSteps = 3;
    /// Refuse a result of more tetrahedra (by `std::length_error`).
    std::int64_t mMaxCells = 50000000;
    /// The cell data naming each triangle's CAD face, carried to the
    /// triangles and made Side regions when present.
    std::string mFaceArray = "cad:face";
};

struct TetrahedralizeResult {
    Mesh mMesh;
    std::int64_t mNumTetrahedra = 0;
    std::int64_t mNumPoints = 0;     ///< every point of the result
    std::int64_t mNumNewPoints = 0;  ///< of those, inside the surface
    std::int64_t mNumRegions = 0;    ///< enclosed regions, one per shell
    std::int64_t mNumBoundaryTriangles = 0;
    /// Every surface triangle and its points kept exactly (checked).
    bool mBoundaryPreserved = false;
    /// The two-sided distance between the result's boundary and the surface:
    /// zero when the boundary is preserved (Netgen), measured by the octree
    /// backend.
    double mHausdorffToSurface = 0.0;
    double mMinDihedralDeg = 0.0;
    double mMaxDihedralDeg = 0.0;
    double mVolume = 0.0;  ///< the tetrahedra's total
    std::string mBackend;  ///< the backend that ran
    std::vector<std::string> mWarnings;
};

/// Whether the Netgen backend is built in (the octree backend always is).
MESHIOPLUSPLUS_API bool tetrahedralize_has_netgen();

/**
 * @brief Fill the closed triangle surface `rSurface` with tetrahedra (see
 * the file comment).
 * @throws std::invalid_argument for an option out of range, a surface with
 * no triangle, a quad, polygon, higher-order or volume block (lines and
 * vertices are skipped), or a surface that is open, non-manifold,
 * misoriented or has degenerate triangles; `Unsupported` for Netgen when it
 * is not built in;
 * `std::length_error` past `mMaxCells`; `std::runtime_error` when the
 * backend fails (a self-intersecting surface, for one).
 */
MESHIOPLUSPLUS_API TetrahedralizeResult tetrahedralize(const Mesh& rSurface,
                                                       const TetrahedralizeOptions& rOptions = {});

}  // namespace meshioplusplus

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
 * @file operations/mass_properties.hpp
 * @brief Mass, centre of mass and inertia of a body (roadmap 9.4).
 *
 * `compute_stats` reports the point mean as its `centroid` and `measure` is
 * GD&T; neither is the centre of mass of a body. `mass_properties(mesh,
 * options)` integrates a density over the mesh and returns, for the whole
 * selection and for every named `Cell` region, the measure (volume, or area for
 * a shell), the mass, the centre of mass, the inertia tensor, its principal
 * moments and axes, and the radii of gyration.
 *
 * ### What is integrated
 *
 * - **Volume cells** (every tabulated 3D type with a face table, and
 *   polyhedra): each cell is the corner-average fan of `detail/polyhedron.hpp`
 *   (the decomposition `poly_measure` and OpenFOAM use), integrated exactly per
 *   tetrahedron for a constant density. A curved (higher-order) cell is
 *   measured on its corners, as everywhere else in the library; a cell with no
 *   face table is skipped with a warning.
 * - **A closed surface** (the highest-dimensional cells are triangles or
 *   polygons that form a closed, consistently wound shell): the body they
 *   enclose, by the divergence theorem, exact for a polyhedron. An inward-wound
 *   shell is turned outward. An open or inconsistently wound surface throws by
 *   name, pointing at `repair`, unless `mThickness` is set.
 * - **A shell** (`mThickness` above zero on surface cells): a thin plate of that
 *   thickness, so the mass is `density * thickness * area`. The measure is the
 *   area. Open surfaces are allowed.
 *
 * ### Density
 *
 * `mDensity` is the density of every cell; `mDensityArray` names a scalar
 * `cell_data` array that replaces it per cell (a closed surface has no cells of
 * its own to carry a density, so the array is refused there). A non-finite
 * density throws.
 *
 * ### The inertia tensor
 *
 * `I = integral(rho * (|r|^2 * identity - r r^T))` as six components in the
 * `xx yy zz xy yz zx` order of the data model, so `Ixy = -integral(rho x y)`.
 * `mAbout` chooses the reference point: the centre of mass (default) or the
 * origin. The principal moments are ascending, the axes the columns of a
 * row-major 3x3 (`detail::sym3_eigen`'s convention), the radii of gyration
 * `sqrt(I_k / mass)`.
 *
 * The mesh is never modified. Per-cell contributions are computed in
 * parallel and summed serially in cell order, so the report is identical on
 * every backend and thread count.
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

/// The point the inertia tensor is taken about.
enum class MassAbout : std::int32_t {
    Centroid = 0,  ///< the centre of mass
    Origin = 1,    ///< the coordinate origin
};

/// What `mass_properties` integrates.
struct MassPropertiesOptions {
    /// Density of every cell (volume cells, a closed surface's body, a shell).
    double mDensity = 1.0;
    /// A scalar `cell_data` array of per-cell densities; empty uses `mDensity`.
    std::string mDensityArray;
    /// Restrict the whole-mesh entry to this `Cell` region; empty is every cell.
    std::string mRegion;
    /// Above zero, surface cells are a shell of this thickness (open allowed).
    double mThickness = 0.0;
    MassAbout mAbout = MassAbout::Centroid;
};

/// The properties of one set of cells.
struct MassProperties {
    /// The region's name; empty for the whole-mesh entry.
    std::string mName;
    /// Cells that took part.
    std::int64_t mNumCells = 0;
    /// Volume (solid) or area (shell).
    double mMeasure = 0.0;
    double mMass = 0.0;
    /// Centre of mass.
    std::array<double, 3> mCentre{0.0, 0.0, 0.0};
    /// Inertia tensor `xx yy zz xy yz zx` about `mAbout`.
    std::array<double, 6> mInertia{};
    /// Principal moments, ascending.
    std::array<double, 3> mPrincipalMoments{};
    /// Principal axes: a row-major 3x3 whose column k belongs to moment k.
    std::array<double, 9> mPrincipalAxes{};
    /// `sqrt(principal moment / mass)`.
    std::array<double, 3> mGyrationRadii{};
};

/// What `mass_properties` computed.
struct MassPropertiesReport {
    /// `"solid"` (volume cells or a closed surface) or `"shell"`.
    std::string mKind;
    /// The whole selection (`mName` empty).
    MassProperties mTotal;
    /// One entry per named `Cell` region, in `Mesh::RegionNames()` order; empty
    /// when `mRegion` restricts the selection.
    std::vector<MassProperties> mRegions;
    std::vector<std::string> mWarnings;
};

/** @brief `"centroid"` or `"origin"`; anything else throws naming them. */
MESHIOPLUSPLUS_API MassAbout mass_about_from_name(const std::string& rName);

/**
 * @brief Mass, centre of mass and inertia of a mesh (see the file comment).
 * @throws std::invalid_argument on an unknown region or density array, a
 *         non-finite density, no cell that can be integrated, an open or
 *         inconsistently wound surface without `mThickness`, a negative
 *         thickness, or a density array on a closed surface.
 */
MESHIOPLUSPLUS_API MassPropertiesReport mass_properties(const Mesh& rMesh,
                                                        const MassPropertiesOptions& rOptions = {});

}  // namespace meshioplusplus

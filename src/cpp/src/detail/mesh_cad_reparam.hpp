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
 * @file detail/mesh_cad_reparam.hpp
 * @brief The discrete parameterisation of `mesh_cad_surface` (roadmap 6.4.4,
 * and the common one of merged faces, 6.5.2): a face whose parameters no
 * metric rescues is meshed in the plane of a flattening of its own coarse
 * mesh instead (Remacle, Geuzaine, Compère and Marchandise 2010), and every
 * new point is carried back through its coarse triangle to the face's
 * parameters and lifted with the face's own surface, so it lies on the exact
 * geometry.
 *
 * The coarse mesh is the face's unrefined triangulation in its parameters,
 * each triangle keeping its corners' parameters. The consecutive loop samples
 * of one mesh point (a pole, the corner a degenerate trim collapses) are one
 * coarse vertex, and the triangle a pole collapses is dropped; in the
 * triangles that keep a corner at the pole the parameters are interpolated in
 * collapsed coordinates (Duffy 1982), the pole's corner taken at the point of
 * its degenerate segment nearest the opposite side, so the map stays
 * continuous across the edges out of the pole and covers the region the
 * dropped triangle spanned.
 *
 * The flattening is the least-squares conformal map (Lévy, Petitjean, Ray and
 * Maillot 2002), which keeps corner angles and needs no disc, accepted when it
 * flips no triangle and its loops neither cross nor touch (exactly, by
 * `cdt2d_find_crossings`): a locally injective map whose boundary images are
 * simple and disjoint is injective. Otherwise mean-value coordinates with the
 * boundary on a circle (Floater 2003), bijective for any disc, holes filled by
 * fans, but only where every boundary corner is at least 60 degrees on the
 * surface, since the circle flattens corners.
 *
 * A **core-private** header (the `mesh_cad_size.hpp` precedent): nothing
 * installed names it. Deterministic: the solvers are `sparse_solve.hpp`'s and
 * every loop runs in index order.
 */

// System includes
#include <array>
#include <cstddef>
#include <cstdint>
#include <map>
#include <string>
#include <utility>
#include <vector>

// Project includes (private, not installed)
#include "edge_sampler.hpp"

namespace meshioplusplus {
namespace detail {

/// A face's coarse mesh, in its parameters (or, merged, several faces').
struct McadCoarseMesh {
    using V2 = std::array<double, 2>;
    using V3 = std::array<double, 3>;
    std::vector<V3> mX;             ///< per vertex: its point
    std::vector<std::int64_t> mId;  ///< per vertex: its mesh point id (< 0: interior)
    /// Per triangle: its vertices, counter-clockwise in its owner's
    /// parameters (or, `mReversed` owner, clockwise), its corners'
    /// parameters, its owner (an index into the surfaces) and the corner at
    /// a pole (-1 for none).
    std::vector<std::array<std::size_t, 3>> mTris;
    std::vector<std::array<V2, 3>> mCornerUv;
    std::vector<std::int32_t> mOwner;
    std::vector<std::int8_t> mPoleCorner;
    /// Per owner: its triangles run clockwise in its parameters (merged
    /// faces oriented alike on the surface, 6.5.2).
    std::vector<char> mReversed;
    /// (owner, vertex) of a pole: its degenerate segment in the parameters.
    std::map<std::pair<std::int32_t, std::size_t>, std::array<V2, 2>> mPoleSegment;
    /// (use, sample) of every loop and embedded sample: its vertex.
    std::map<std::pair<std::size_t, std::size_t>, std::size_t> mOfSample;
    /// The first owner's plane (the Cholesky factor of its mean first
    /// fundamental form), which its condition number is measured against.
    double mL11 = 1.0;
    double mL12 = 0.0;
    double mL22 = 1.0;
    /// Filled by a meshing that reached its lattice; empty `mWhy` when it can
    /// be flattened at all.
    bool mFilled = false;
    std::string mWhy;
};

/// The coarse mesh flattened, and the way back from its plane.
class McadDisc {
public:
    using V2 = std::array<double, 2>;
    using V3 = std::array<double, 3>;

    /// Where a point of the disc lies: the coarse triangle's owner, the
    /// parameters there and their derivative along the disc,
    /// `{du/dp, du/dq, dv/dp, dv/dq}`; `mInside` false when the point is
    /// outside every triangle (the nearest one's clamped values).
    struct Hit {
        bool mInside = false;
        std::int32_t mOwner = 0;
        std::size_t mTri = 0;
        V2 mUv{0.0, 0.0};
        std::array<double, 4> mD{1.0, 0.0, 0.0, 1.0};
    };

    const McadCoarseMesh* pCoarse = nullptr;
    std::vector<const BrepSurfaceEval*> mSurfaces;  ///< per owner
    std::vector<V2> mP;                             ///< per coarse vertex: its disc point
    /// The Cholesky factor of the area-weighted mean of `J'^T J'` (J' the
    /// surface's derivative along the disc): the face's plane is `L' p`.
    double mL11 = 1.0;
    double mL12 = 0.0;
    double mL22 = 1.0;
    /// `(S_p, S_q)` and `S` at each coarse triangle's centroid, for the
    /// lattice's spacing and the stretched test.
    std::vector<std::array<V3, 2>> mJac;
    std::vector<V3> mJacAt;
    int mMethod = 0;  ///< 1 least-squares conformal, 2 mean-value on a circle

    Hit Locate(const V2& rP) const;
    /// The parameters at barycentrics `rL` of triangle `T`.
    V2 UvAt(std::size_t T, const std::array<double, 3>& rL) const;

    /// Build the bucket grid (after `mP`).
    void Index();

private:
    std::array<double, 3> Bary(std::size_t T, const V2& rP) const;
    V2 mLo{0.0, 0.0};
    double mCell = 1.0;
    std::size_t mNx = 1;
    std::size_t mNy = 1;
    std::vector<std::vector<std::size_t>> mBuckets;
};

/// The face's condition number on its coarse mesh (roadmap 6.4.4): the
/// largest stretch of `J L^-1` over its smallest, `J` the surface's
/// derivative at each triangle's centroid and `L` the face's plane, so a
/// metric that varies in direction or in size counts. Infinite where
/// `J L^-1` is singular.
double mcad_coarse_condition(const McadCoarseMesh& rCoarse,
                             const std::vector<const BrepSurfaceEval*>& rSurfaces, double L11,
                             double L12, double L22);

/// The same measure on a disc, against its own plane `L'`.
double mcad_disc_condition(const McadDisc& rDisc);

/// Flatten the coarse mesh; empty on success, else why not.
std::string mcad_build_disc(const McadCoarseMesh& rCoarse,
                            const std::vector<const BrepSurfaceEval*>& rSurfaces, McadDisc& rDisc);

}  // namespace detail
}  // namespace meshioplusplus

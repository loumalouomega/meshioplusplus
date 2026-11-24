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
 * @file detail/point_locator.hpp
 * @brief Which simplex of a simplexified mesh contains a point, and where in it.
 *
 * A **core-private** header: nothing installed names it. It is the search
 * `interpolate`'s barycentric mode ran inline (roadmap 9.3), moved here so
 * `locate_points` shares it instead of repeating it:
 *
 * - the simplices are the mesh's tetrahedra, or, with none, its triangles
 *   (`convert_cells(Simplexify)` has already reduced every cell to them), in
 *   ascending block-then-cell order -- the simplex index every tie-break
 *   refers to;
 * - each simplex's bounding box is quantised onto a grid sized by
 *   `point_grid_cell_size` and grouped by `build_bucket_table`, every bucket
 *   ascending;
 * - a query takes the first containing simplex of its bucket in ascending
 *   index, by barycentric weights that must all be at least `-tolerance`
 *   (a tetrahedron in three dimensions, a triangle in the xy-plane, `z`
 *   ignored: the documented planar assumption of the triangle path).
 *
 * Header-only and inline on purpose: `interpolate` keeps its instruction
 * count and its digests, which a call across a translation unit would risk.
 * `interpolate`'s numpy twin repeats the same expressions.
 */

// System includes
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/detail/bucket_table.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/parallel.hpp"

// Project includes (private, not installed)
#include "bounding_box_utils.hpp"
#include "bucket_build.hpp"
#include "typed_view.hpp"

namespace meshioplusplus {
namespace detail {

/// The default barycentric containment tolerance: dimensionless (the weights
/// are volume ratios), hence scale-invariant. The numpy twin of `interpolate`
/// has the same literal.
inline constexpr double kPointLocatorTolerance = 1e-12;

/**
 * @brief The grid cell size for @p N items spread over the bounding box
 * `[pLo, pHi]`: the largest axis extent divided by the smallest integer R with
 * `R*R*R >= N`. Deliberately no cbrt or pow: the integer loop is exact in C++
 * and Python, and the single division is IEEE-identical, which keeps the numpy
 * fallback's grid byte-compatible.
 */
inline double point_grid_cell_size(const double* pLo, const double* pHi, std::size_t N) {
    double max_extent = 0.0;
    for (int d = 0; d < 3; ++d) {
        const double e = pHi[d] - pLo[d];
        if (e > max_extent)
            max_extent = e;
    }
    std::int64_t r = 1;
    while (r * r * r < static_cast<std::int64_t>(N))
        ++r;
    return max_extent > 0.0 ? max_extent / static_cast<double>(r) : 1.0;
}

/// Barycentric weights of @p rP in the tetrahedron `pNodes[0..3]`; true iff
/// contained (every weight at least `-Tolerance`). A degenerate tetrahedron
/// (exactly zero volume) contains nothing.
inline bool tet_contains(const NDArray& rPts, std::size_t PointDim, const std::int64_t* pNodes,
                         const Vec3& rP, double* pW, double Tolerance) {
    const Vec3 a = read_point(rPts, PointDim, pNodes[0]);
    const Vec3 b = read_point(rPts, PointDim, pNodes[1]);
    const Vec3 c = read_point(rPts, PointDim, pNodes[2]);
    const Vec3 d = read_point(rPts, PointDim, pNodes[3]);
    const double v = triple_product(vec3_sub(b, a), vec3_sub(c, a), vec3_sub(d, a));
    if (v == 0.0)
        return false;
    const double wa = triple_product(vec3_sub(b, rP), vec3_sub(c, rP), vec3_sub(d, rP)) / v;
    const double wb = triple_product(vec3_sub(rP, a), vec3_sub(c, a), vec3_sub(d, a)) / v;
    const double wc = triple_product(vec3_sub(b, a), vec3_sub(rP, a), vec3_sub(d, a)) / v;
    const double wd = triple_product(vec3_sub(b, a), vec3_sub(c, a), vec3_sub(rP, a)) / v;
    if (wa < -Tolerance || wb < -Tolerance || wc < -Tolerance || wd < -Tolerance)
        return false;
    pW[0] = wa;
    pW[1] = wb;
    pW[2] = wc;
    pW[3] = wd;
    return true;
}

/// Twice the signed area of the triangle `p0 p1 p2` **in the xy-plane** (the
/// documented planar assumption of the triangle path; z is ignored).
inline double xy_area2(const Vec3& rP0, const Vec3& rP1, const Vec3& rP2) {
    return (rP1[0] - rP0[0]) * (rP2[1] - rP0[1]) - (rP1[1] - rP0[1]) * (rP2[0] - rP0[0]);
}

/// Barycentric weights of @p rP in the (xy-plane) triangle `pNodes[0..2]`;
/// true iff contained.
inline bool tri_contains_xy(const NDArray& rPts, std::size_t PointDim, const std::int64_t* pNodes,
                            const Vec3& rP, double* pW, double Tolerance) {
    const Vec3 a = read_point(rPts, PointDim, pNodes[0]);
    const Vec3 b = read_point(rPts, PointDim, pNodes[1]);
    const Vec3 c = read_point(rPts, PointDim, pNodes[2]);
    const double v = xy_area2(a, b, c);
    if (v == 0.0)
        return false;
    const double wa = xy_area2(rP, b, c) / v;
    const double wb = xy_area2(a, rP, c) / v;
    const double wc = xy_area2(a, b, rP) / v;
    if (wa < -Tolerance || wb < -Tolerance || wc < -Tolerance)
        return false;
    pW[0] = wa;
    pW[1] = wb;
    pW[2] = wc;
    return true;
}

/// The search structure over one simplexified mesh (see the file comment).
class PointLocator {
public:
    /**
     * @param rSimplified a mesh already reduced to simplices by
     *        `convert_cells(Simplexify)`; it must outlive the locator.
     * @param Tolerance the barycentric containment tolerance.
     * @param rWho the operation's name for the error message.
     * @throws std::invalid_argument when the mesh has neither tetrahedra nor
     *         triangles.
     */
    PointLocator(const Mesh& rSimplified, double Tolerance, const std::string& rWho)
        : mpPoints(&rSimplified.Points()),
          mPointDim(rSimplified.PointDim()),
          mTolerance(Tolerance) {
        bool has_tet = false;
        bool has_tri = false;
        for (const auto cb : rSimplified.CellRange()) {
            if (cb.NumCells() == 0)
                continue;
            if (cb.Type() == "tetra")
                has_tet = true;
            else if (cb.Type() == "triangle")
                has_tri = true;
        }
        if (!has_tet && !has_tri)
            throw std::invalid_argument(
                "meshio++: " + rWho +
                ": barycentric needs a source with triangle or tetrahedron cells after "
                "simplexification");
        mNodesPer = has_tet ? 4 : 3;
        const std::string want = has_tet ? "tetra" : "triangle";

        // Flat simplex connectivity, ascending block-then-cell order.
        for (const auto cb : rSimplified.CellRange()) {
            if (cb.Type() != want || cb.NumCells() == 0)
                continue;
            const NDArray& conn = cb.Conn();
            const Int64View conn_v(conn);
            const std::size_t n = cb.NumCells() * mNodesPer;
            const std::size_t base = mConn.size();
            mConn.resize(base + n);
            for (std::size_t j = 0; j < n; ++j)
                mConn[base + j] = conn_v[j];
        }
        const std::size_t nsimp = mConn.size() / mNodesPer;

        const std::size_t ns = rSimplified.NumPoints();
        const NDArray& spts = *mpPoints;
        const std::size_t sdim = mPointDim;
        auto scoord = [&](std::int64_t i) { return read_point(spts, sdim, i); };
        BoundingBoxUtils::FromCoords(
            ns, [&](std::size_t i) { return scoord(static_cast<std::int64_t>(i)); }, mLo, mHi);
        const double cell = point_grid_cell_size(mLo, mHi, nsimp);

        // Per-simplex quantised bbox keys, filled in parallel and grouped into
        // the frozen grid: each bucket ascending, as serial inserts in simplex
        // order left it.
        std::vector<GridKey> box_lo(nsimp);
        std::vector<GridKey> box_hi(nsimp);
        const auto key_of = [cell](const double* pC) {
            return GridKey{grid_quantize(pC[0], cell), grid_quantize(pC[1], cell),
                           grid_quantize(pC[2], cell)};
        };
        const std::size_t nv = mNodesPer;
        parallel_for(nsimp, [&](std::size_t s) {
            double blo[3], bhi[3];
            BoundingBoxUtils::FromCoords(
                nv, [&](std::size_t k) { return scoord(mConn[s * nv + k]); }, blo, bhi);
            box_lo[s] = key_of(blo);
            box_hi[s] = key_of(bhi);
        });
        mTable = build_bucket_table(cell, box_lo, box_hi);
    }

    /// 3 (triangles) or 4 (tetrahedra).
    std::size_t NodesPerSimplex() const { return mNodesPer; }
    std::size_t NumSimplices() const { return mConn.size() / mNodesPer; }
    /// The flat connectivity: simplex `s` is `[s * NodesPerSimplex(), ...)`.
    const std::vector<std::int64_t>& Connectivity() const { return mConn; }
    /// The bounding box of the simplexified mesh's points.
    const double* Lo() const { return mLo; }
    const double* Hi() const { return mHi; }

    /**
     * @brief The first simplex containing @p rP, in ascending index.
     * @param pW receives its barycentric weights (`NodesPerSimplex()` of them).
     * @return the simplex index, or -1 when no simplex contains the point.
     */
    std::int64_t Locate(const Vec3& rP, double* pW) const {
        BucketView cand(nullptr, 0);
        if (!mTable.Find(mTable.KeyOf(rP.data()), cand))
            return -1;
        const std::size_t nv = mNodesPer;
        for (std::int64_t s : cand) {
            const std::int64_t* nodes = &mConn[static_cast<std::size_t>(s) * nv];
            const bool inside = nv == 4 ? tet_contains(*mpPoints, mPointDim, nodes, rP, pW, mTolerance)
                                        : tri_contains_xy(*mpPoints, mPointDim, nodes, rP, pW, mTolerance);
            if (inside)
                return s;
        }
        return -1;
    }

private:
    const NDArray* mpPoints;
    std::size_t mPointDim;
    double mTolerance;
    std::size_t mNodesPer = 3;
    std::vector<std::int64_t> mConn;
    double mLo[3] = {0.0, 0.0, 0.0};
    double mHi[3] = {0.0, 0.0, 0.0};
    BucketTable mTable;
};

}  // namespace detail
}  // namespace meshioplusplus

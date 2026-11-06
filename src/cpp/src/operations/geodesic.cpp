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
// Geodesic distance on a triangle surface by the heat method. See
// operations/geodesic.hpp for the contract.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/geodesic.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/facet_index.hpp"
#include "meshioplusplus/detail/node_adjacency.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/operations/region_ops.hpp"
#include "meshioplusplus/region.hpp"
#include "meshioplusplus/detail/geometry.hpp"

// Project includes (private, not installed)
#include "../detail/sparse_solve.hpp"
#include "../detail/surface_fem.hpp"

namespace meshioplusplus {

namespace {

using GeoVec = std::array<double, 3>;

// One surface triangle and its cotangents (at corners 0, 1, 2).
using GeoTriangle = detail::SurfaceFemTriangle;

// The points a region names (the periodic operation's rule).
std::vector<std::int64_t> geo_region_points(const Mesh& rMesh, const Region& rR) {
    const std::size_t n = rMesh.NumPoints();
    std::vector<std::int64_t> out, ids;
    const std::size_t count = rR.NumEntries();
    if (rR.mKind == RegionKind::Point) {
        for (std::size_t i = 0; i < count; ++i)
            out.push_back(detail::read_int(rR.mEntries, i));
    } else if (rR.mKind == RegionKind::Cell) {
        const std::vector<std::int64_t> bases = detail::block_bases(rMesh);
        for (std::size_t i = 0; i < count; ++i) {
            const auto [block, row] =
                detail::global_to_block_row(bases, detail::read_int(rR.mEntries, i));
            if (block == static_cast<std::size_t>(-1))
                continue;
            detail::cell_node_ids(rMesh.Cells(block), static_cast<std::size_t>(row), n, ids);
            out.insert(out.end(), ids.begin(), ids.end());
        }
    } else {
        CellType type;
        for (std::size_t i = 0; i < count; ++i)
            if (detail::facet_nodes(rMesh, detail::read_int(rR.mEntries, i * 2),
                                    detail::read_int(rR.mEntries, i * 2 + 1), type, ids))
                out.insert(out.end(), ids.begin(), ids.end());
    }
    return out;
}

std::size_t geo_find(std::vector<std::size_t>& rParent, std::size_t X) {
    while (rParent[X] != X) {
        rParent[X] = rParent[rParent[X]];
        X = rParent[X];
    }
    return X;
}

// Gauss-Seidel sweeps on @p rA x = @p rB from x = 0, in point order (serial, so
// the result does not depend on a schedule). Converges for any symmetric
// positive definite matrix, obtuse triangles' negative weights included. True
// once no value changed by more than `Tolerance` relative to itself in a sweep;
// false on the cap or a value that is not finite.
bool geo_gauss_seidel_heat(const detail::SparseCsr& rA, const std::vector<double>& rB,
                           std::vector<double>& rX, double Tolerance, std::size_t Cap,
                           std::size_t& rSweeps) {
    const std::size_t n = rA.mRows;
    std::vector<double> diag(n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t k = rA.mStart[i]; k < rA.mStart[i + 1]; ++k)
            if (rA.mCol[k] == i)
                diag[i] = rA.mValue[k];
    for (rSweeps = 1; rSweeps <= Cap; ++rSweeps) {
        bool moved = false;
        for (std::size_t i = 0; i < n; ++i) {
            double sum = rB[i];
            for (std::size_t k = rA.mStart[i]; k < rA.mStart[i + 1]; ++k)
                if (rA.mCol[k] != i)
                    sum -= rA.mValue[k] * rX[rA.mCol[k]];
            const double v = sum / diag[i];
            if (!std::isfinite(v))
                return false;
            if (std::fabs(v - rX[i]) > Tolerance * std::fabs(v))
                moved = true;
            rX[i] = v;
        }
        if (!moved)
            return true;
    }
    return false;
}

}  // namespace

GeodesicResult geodesic_distance(const Mesh& rMesh, const GeodesicOptions& rOptions) {
    if (!(rOptions.mTimeFactor > 0.0) || !std::isfinite(rOptions.mTimeFactor))
        throw std::invalid_argument("geodesic_distance: the time factor must be positive");
    if (!(rOptions.mTolerance > 0.0))
        throw std::invalid_argument("geodesic_distance: the tolerance must be positive");
    const std::size_t n = rMesh.NumPoints();
    const detail::SurfaceFem fem = detail::build_surface_fem(rMesh, "geodesic_distance");
    const std::vector<GeoVec>& xyz = fem.mPoints;
    const std::vector<GeoTriangle>& good = fem.mTriangles;
    GeodesicResult result;
    result.mNumDegenerate = fem.mNumDegenerate;
    const double h = fem.mMeanEdge;
    result.mTimeStep = rOptions.mTimeFactor * h * h;

    // --- on the surface, and its connected pieces -------------------------
    std::vector<char> on_surface(n, 0);
    std::vector<std::size_t> parent(n);
    std::iota(parent.begin(), parent.end(), std::size_t{0});
    for (const GeoTriangle& t : good)
        for (std::size_t k = 0; k < 3; ++k) {
            on_surface[t.mNodes[k]] = 1;
            const std::size_t r0 = geo_find(parent, t.mNodes[0]);
            const std::size_t rk = geo_find(parent, t.mNodes[k]);
            if (r0 != rk)
                parent[std::max(r0, rk)] = std::min(r0, rk);
        }

    // --- the sources --------------------------------------------------------
    std::vector<std::int64_t> sources = rOptions.mSources;
    if (!rOptions.mSourceRegion.empty()) {
        RegionSelector selector;
        selector.mName = rOptions.mSourceRegion;
        const std::vector<std::int64_t> extra =
            geo_region_points(rMesh, rMesh.Region(find_region(rMesh, selector)));
        sources.insert(sources.end(), extra.begin(), extra.end());
    }
    std::sort(sources.begin(), sources.end());
    sources.erase(std::unique(sources.begin(), sources.end()), sources.end());
    if (sources.empty())
        throw std::invalid_argument("geodesic_distance: give at least one source point");
    for (const std::int64_t s : sources) {
        if (s < 0 || static_cast<std::size_t>(s) >= n)
            throw std::invalid_argument("geodesic_distance: source point " + std::to_string(s) +
                                        " is outside the " + std::to_string(n) + " points");
        if (!on_surface[static_cast<std::size_t>(s)])
            throw std::invalid_argument("geodesic_distance: source point " + std::to_string(s) +
                                        " is in no triangle");
    }
    result.mNumSources = static_cast<std::int64_t>(sources.size());

    // --- 1. heat: (M + t K) u = delta ---------------------------------------
    std::vector<double> mass;
    std::vector<detail::SparseTriplet> stiffness;
    detail::surface_fem_assemble(fem, mass, stiffness);
    const std::size_t cap = rOptions.mMaxIterations > 0
                                ? static_cast<std::size_t>(rOptions.mMaxIterations)
                                : std::max<std::size_t>(2 * n + 10, 1000);
    std::vector<double> u(n, 0.0);
    {
        std::vector<detail::SparseTriplet> heat;
        heat.reserve(stiffness.size() + n);
        for (const detail::SparseTriplet& e : stiffness)
            heat.push_back({e.mRow, e.mCol, result.mTimeStep * e.mValue});
        std::vector<double> rhs(n, 0.0);
        for (std::size_t i = 0; i < n; ++i)
            // A point in no triangle is an identity row: it stays at 0.
            heat.push_back({i, i, on_surface[i] ? mass[i] : 1.0});
        for (const std::int64_t s : sources)
            rhs[static_cast<std::size_t>(s)] = 1.0;
        const detail::SparseCsr a = detail::csr_from_triplets(n, n, std::move(heat));
        // Gauss-Seidel sweeps, not conjugate gradients: far from the sources
        // the heat falls to 1e-30 and below, which a residual-based stop
        // resolves to noise, while a sweep that stops on every value's own
        // relative change resolves each value however small (a direct
        // factorisation's property, which Crane et al. rely on).
        std::size_t sweeps = 0;
        const bool settled = geo_gauss_seidel_heat(a, rhs, u, rOptions.mTolerance,
                                                   std::max<std::size_t>(cap, 64), sweeps);
        result.mHeatIterations = static_cast<std::int64_t>(sweeps);
        result.mConverged = result.mConverged && settled;
    }

    // --- 2. the normalised direction, and its integrated divergence -------
    std::vector<double> div(n, 0.0);
    for (const GeoTriangle& t : good) {
        const GeoVec& p0 = xyz[t.mNodes[0]];
        const GeoVec& p1 = xyz[t.mNodes[1]];
        const GeoVec& p2 = xyz[t.mNodes[2]];
        const GeoVec normal =
            detail::vec3_cross(detail::vec3_sub(p1, p0), detail::vec3_sub(p2, p0));
        const double area2 = t.mDoubleArea;
        GeoVec unit_n{normal[0] / area2, normal[1] / area2, normal[2] / area2};
        // grad u = sum_i u_i (N x e_i) / (2 A), e_i the edge opposite corner i,
        // oriented counter-clockwise.
        GeoVec grad{0.0, 0.0, 0.0};
        const std::array<const GeoVec*, 3> p = {&p0, &p1, &p2};
        for (std::size_t k = 0; k < 3; ++k) {
            const GeoVec e = detail::vec3_sub(*p[(k + 2) % 3], *p[(k + 1) % 3]);
            const GeoVec ne = detail::vec3_cross(unit_n, e);
            const double uk = u[t.mNodes[k]];
            for (std::size_t c = 0; c < 3; ++c)
                grad[c] += uk * ne[c] / area2;
        }
        const double len = detail::vec3_norm(grad);
        if (!(len > 0.0))
            continue;
        const GeoVec x{-grad[0] / len, -grad[1] / len, -grad[2] / len};
        for (std::size_t k = 0; k < 3; ++k) {
            // Corner k, its two neighbours j (next) and l (after): the edges from
            // k to them, weighted by the cotangents of the angles opposite them.
            const std::size_t j = (k + 1) % 3, l = (k + 2) % 3;
            const GeoVec e1 = detail::vec3_sub(*p[j], *p[k]);
            const GeoVec e2 = detail::vec3_sub(*p[l], *p[k]);
            div[t.mNodes[k]] +=
                0.5 * (t.mCot[l] * detail::vec3_dot(e1, x) + t.mCot[j] * detail::vec3_dot(e2, x));
        }
    }

    // --- 3. Poisson: K phi = -div, one point of each piece held at 0 -------
    std::vector<std::size_t> pin_of_root(n, static_cast<std::size_t>(-1));
    std::vector<char> root_has_source(n, 0);
    for (const std::int64_t s : sources) {
        const std::size_t root = geo_find(parent, static_cast<std::size_t>(s));
        root_has_source[root] = 1;
        if (pin_of_root[root] == static_cast<std::size_t>(-1))
            pin_of_root[root] = static_cast<std::size_t>(s);
    }
    std::vector<char> pinned(n, 0);
    for (std::size_t i = 0; i < n; ++i) {
        if (!on_surface[i]) {
            pinned[i] = 1;
            continue;
        }
        const std::size_t root = geo_find(parent, i);
        if (pin_of_root[root] == static_cast<std::size_t>(-1))
            pin_of_root[root] = i;  // a piece with no source: any point, then NaN
        if (pin_of_root[root] == i)
            pinned[i] = 1;
    }
    std::vector<double> phi(n, 0.0);
    {
        std::vector<detail::SparseTriplet> poisson;
        poisson.reserve(stiffness.size() + n);
        for (const detail::SparseTriplet& e : stiffness)
            if (!pinned[e.mRow] && !pinned[e.mCol])
                poisson.push_back(e);
        std::vector<double> rhs(n, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            if (pinned[i])
                poisson.push_back({i, i, 1.0});
            else
                rhs[i] = -div[i];
        }
        const detail::SparseCsr a = detail::csr_from_triplets(n, n, std::move(poisson));
        const detail::SparseSolveReport r =
            detail::csr_pcg(a, rhs.data(), phi.data(), rOptions.mTolerance, cap);
        result.mPoissonIterations = static_cast<std::int64_t>(r.mIterations);
        result.mConverged = result.mConverged && r.mConverged;
    }
    if (!result.mConverged)
        log::warn("geodesic_distance: a solve hit its iteration cap; the distance is approximate");

    // --- the least value of each piece is 0; sourceless pieces are NaN ----
    std::vector<double> least(n, std::numeric_limits<double>::infinity());
    for (std::size_t i = 0; i < n; ++i)
        if (on_surface[i]) {
            const std::size_t root = geo_find(parent, i);
            least[root] = std::min(least[root], phi[i]);
        }
    NDArray distance(DType::Float64, {n});
    double* pd = distance.As<double>();
    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t root = on_surface[i] ? geo_find(parent, i) : 0;
        if (!on_surface[i] || !root_has_source[root]) {
            pd[i] = std::numeric_limits<double>::quiet_NaN();
            ++result.mNumUnreached;
            continue;
        }
        pd[i] = phi[i] - least[root];
        result.mMaxDistance = std::max(result.mMaxDistance, pd[i]);
    }
    result.mMesh = detail::clone_mesh(rMesh);
    result.mMesh.AddPointData(kGeodesicDistanceName, std::move(distance));
    return result;
}

}  // namespace meshioplusplus

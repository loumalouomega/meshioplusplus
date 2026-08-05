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
// Signed distance to a surface. See operations/sdf.hpp for the contract and
// detail/surface_distance.hpp for the kernel this file drives.
//
// Anonymous-namespace helpers are prefixed `sdfop_`; `sdf_` is left free for the
// detail-layer names so a reader can tell which layer a symbol belongs to.

// System includes
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/sdf.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/surface_distance.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {

namespace {

constexpr const char* kSdfPrefix = "meshio++: sdf: ";

// The query points a mesh contributes, as a flat Float64 (n, 3) array: its own
// points, or its cell centroids.
std::vector<detail::Vec3> sdfop_query_points(const Mesh& rMesh, SdfLocation Location) {
    std::vector<detail::Vec3> out;
    const std::size_t dim = rMesh.PointDim();
    if (Location == SdfLocation::Corner) {
        const std::size_t n = rMesh.NumPoints();
        out.resize(n);
        const NDArray& points = rMesh.Points();
        parallel_for_bw(n, [&](std::size_t p) {
            out[p] = detail::read_point(points, dim, static_cast<std::int64_t>(p));
        });
        return out;
    }

    // Cell centres: the corner average, accumulated in connectivity order so the
    // numpy twin can reproduce it without a summation-order argument.
    const NDArray& points = rMesh.Points();
    for (const auto cb : rMesh.CellRange()) {
        const std::size_t ncells = cb.NumCells();
        const std::size_t base = out.size();
        out.resize(base + ncells);
        if (cb.IsRagged()) {
            for (std::size_t c = 0; c < ncells; ++c) {
                detail::Vec3 sum{0.0, 0.0, 0.0};
                std::size_t count = 0;
                if (cb.IsPolyhedron()) {
                    for (std::size_t f = 0; f < cb.NumFaces(c); ++f) {
                        const auto face = cb.Face(c, f);
                        for (std::size_t i = 0; i < face.second; ++i) {
                            sum = detail::vec3_add(
                                sum, detail::read_point(points, dim, face.first[i]));
                            ++count;
                        }
                    }
                } else {
                    const std::size_t n = cb.RowSize(c);
                    const std::int64_t* row = cb.Row(c);
                    for (std::size_t i = 0; i < n; ++i) {
                        sum = detail::vec3_add(sum, detail::read_point(points, dim, row[i]));
                        ++count;
                    }
                }
                out[base + c] =
                    count == 0 ? sum : detail::vec3_scale(sum, 1.0 / static_cast<double>(count));
            }
            continue;
        }
        const NDArray& conn = cb.Conn();
        const std::size_t npc = cb.NodesPerCell();
        const double inv = npc == 0 ? 0.0 : 1.0 / static_cast<double>(npc);
        parallel_for_bw(ncells, [&](std::size_t c) {
            detail::Vec3 sum{0.0, 0.0, 0.0};
            for (std::size_t i = 0; i < npc; ++i)
                sum = detail::vec3_add(
                    sum, detail::read_point(points, dim, detail::read_int(conn, c * npc + i)));
            out[base + c] = detail::vec3_scale(sum, inv);
        });
    }
    return out;
}

// Report a surface's defects at the requested severity.
void sdfop_report_quality(const SurfaceQuality& rQuality, SdfWatertightCheck Check) {
    if (Check == SdfWatertightCheck::Off || rQuality.mWatertight)
        return;
    const std::string what =
        "the surface is not watertight: " + std::to_string(rQuality.mBoundaryEdges) +
        " boundary edge(s), " + std::to_string(rQuality.mNonManifoldEdges) +
        " non-manifold edge(s), " + std::to_string(rQuality.mInconsistentPairs) +
        " inconsistently wound pair(s), " + std::to_string(rQuality.mDegenerateTriangles) +
        " degenerate triangle(s)";
    if (Check == SdfWatertightCheck::Error)
        throw std::invalid_argument(std::string(kSdfPrefix) + what);
    log::warn("{}{} -- the sign may be wrong near the defects; consider sign='winding-number'.",
              kSdfPrefix, what);
}

}  // namespace

SdfSign sdf_sign_from_name(const std::string& rName) {
    if (rName == "unsigned")
        return SdfSign::Unsigned;
    if (rName == "pseudonormal")
        return SdfSign::Pseudonormal;
    if (rName == "winding-number" || rName == "winding_number")
        return SdfSign::WindingNumber;
    throw std::invalid_argument(std::string(kSdfPrefix) + "unknown sign '" + rName +
                                "' (expected one of: unsigned, pseudonormal, winding-number)");
}

SdfPseudonormalWeight sdf_weight_from_name(const std::string& rName) {
    if (rName == "angle")
        return SdfPseudonormalWeight::Angle;
    if (rName == "area")
        return SdfPseudonormalWeight::Area;
    throw std::invalid_argument(std::string(kSdfPrefix) + "unknown weight '" + rName +
                                "' (expected one of: angle, area)");
}

SdfLocation sdf_location_from_name(const std::string& rName) {
    if (rName == "corner" || rName == "point")
        return SdfLocation::Corner;
    if (rName == "center" || rName == "centre" || rName == "cell")
        return SdfLocation::Center;
    throw std::invalid_argument(std::string(kSdfPrefix) + "unknown location '" + rName +
                                "' (expected one of: corner, center)");
}

SdfStructure sdf_structure_from_name(const std::string& rName) {
    if (rName == "voxel")
        return SdfStructure::Voxel;
    if (rName == "octree")
        return SdfStructure::Octree;
    throw std::invalid_argument(std::string(kSdfPrefix) + "unknown structure '" + rName +
                                "' (expected one of: voxel, octree)");
}

SdfWatertightCheck sdf_watertight_check_from_name(const std::string& rName) {
    if (rName == "off")
        return SdfWatertightCheck::Off;
    if (rName == "warn")
        return SdfWatertightCheck::Warn;
    if (rName == "error")
        return SdfWatertightCheck::Error;
    throw std::invalid_argument(std::string(kSdfPrefix) + "unknown watertight check '" + rName +
                                "' (expected one of: off, warn, error)");
}

SurfaceQuality surface_watertight_check(const Mesh& rSurface) {
    const detail::TriangleSoup soup = detail::build_triangle_soup(rSurface, "");
    return detail::soup_quality(soup);
}

NDArray sample_distance(const Mesh& rSurface, const NDArray& rPoints,
                        const SurfaceDistanceOptions& rOptions) {
    const std::size_t dim = detail::cols(rPoints);
    if (rPoints.Ndim() != 2 || (dim != 2 && dim != 3))
        throw std::invalid_argument(std::string(kSdfPrefix) +
                                    "query points must be a 2-D (n, 2) or (n, 3) array");
    const std::size_t n = detail::rows(rPoints);
    std::vector<detail::Vec3> queries(n);
    parallel_for_bw(n, [&](std::size_t p) {
        queries[p] = detail::read_point(rPoints, dim, static_cast<std::int64_t>(p));
    });

    const detail::TriangleSoup soup =
        detail::build_triangle_soup(rSurface, rOptions.mSurfaceRegion);
    sdfop_report_quality(detail::soup_quality(soup), rOptions.mWatertightCheck);

    detail::DistanceQuery query = detail::build_distance_query(soup, rOptions);
    std::vector<detail::DistanceHit> hits = detail::query_distances(query, queries, rOptions);

    NDArray out = NDArray::Uninit(DType::Float64, {n});
    double* dst = out.As<double>();
    parallel_for_bw(n, [&](std::size_t p) { dst[p] = hits[p].mSignedDistance; });
    return out;
}

SurfaceDistanceResult distance_to_surface(const Mesh& rQuery, const Mesh& rSurface,
                                          const SurfaceDistanceOptions& rOptions) {
    SurfaceDistanceResult result;
    result.mMesh = detail::clone_mesh(rQuery, [](DataLocation, const std::string&, std::string&) {
        return true;  // geometry and every existing array pass through unchanged
    });

    const detail::TriangleSoup soup =
        detail::build_triangle_soup(rSurface, rOptions.mSurfaceRegion);
    result.mQuality = detail::soup_quality(soup);
    sdfop_report_quality(result.mQuality, rOptions.mWatertightCheck);

    const std::vector<detail::Vec3> queries = sdfop_query_points(rQuery, rOptions.mLocation);
    detail::DistanceQuery query = detail::build_distance_query(soup, rOptions);
    const std::vector<detail::DistanceHit> hits =
        detail::query_distances(query, queries, rOptions);
    const std::size_t n = hits.size();

    NDArray dist = NDArray::Uninit(DType::Float64, {n});
    double* pd = dist.As<double>();
    parallel_for_bw(n, [&](std::size_t p) { pd[p] = hits[p].mSignedDistance; });

    const bool banded = rOptions.mBand > 0.0;
    NDArray band;
    if (banded) {
        band = NDArray::Uninit(DType::Int64, {n});
        std::int64_t* pb = band.As<std::int64_t>();
        parallel_for_bw(n, [&](std::size_t p) { pb[p] = hits[p].mInBand ? 1 : 0; });
        for (std::size_t p = 0; p < n; ++p)
            result.mNumBanded += hits[p].mInBand ? 0 : 1;
    }
    NDArray inside;
    if (rOptions.mRecordInside) {
        inside = NDArray::Uninit(DType::Int64, {n});
        std::int64_t* pi = inside.As<std::int64_t>();
        parallel_for_bw(n, [&](std::size_t p) {
            pi[p] = hits[p].mSignedDistance < 0.0 ? 1 : 0;
        });
    }
    NDArray closest;
    if (rOptions.mRecordClosestCell) {
        closest = NDArray::Uninit(DType::Int64, {n});
        std::int64_t* pc = closest.As<std::int64_t>();
        parallel_for_bw(n, [&](std::size_t p) { pc[p] = hits[p].mSourceCell; });
    }

    if (rOptions.mLocation == SdfLocation::Corner) {
        result.mMesh.AddPointData(kSdfDistanceName, std::move(dist));
        if (banded)
            result.mMesh.AddPointData(kSdfBandName, std::move(band));
        if (rOptions.mRecordInside)
            result.mMesh.AddPointData(kSdfInsideName, std::move(inside));
        if (rOptions.mRecordClosestCell)
            result.mMesh.AddPointData(kSdfClosestCellName, std::move(closest));
        return result;
    }

    // Cell data is per block, so split the flat run back up along the block
    // boundaries it was built from.
    const auto split = [&](NDArray& rFlat) {
        std::vector<NDArray> blocks;
        std::size_t base = 0;
        for (const auto cb : result.mMesh.CellRange()) {
            const std::size_t ncells = cb.NumCells();
            NDArray b = NDArray::Uninit(rFlat.Dtype(), {ncells});
            const std::size_t width = rFlat.Nbytes() / (rFlat.Size() == 0 ? 1 : rFlat.Size());
            std::memcpy(b.Data(), rFlat.Data() + base * width, ncells * width);
            blocks.push_back(std::move(b));
            base += ncells;
        }
        return blocks;
    };
    result.mMesh.AddCellData(kSdfDistanceName, split(dist));
    if (banded)
        result.mMesh.AddCellData(kSdfBandName, split(band));
    if (rOptions.mRecordInside)
        result.mMesh.AddCellData(kSdfInsideName, split(inside));
    if (rOptions.mRecordClosestCell)
        result.mMesh.AddCellData(kSdfClosestCellName, split(closest));
    return result;
}

SdfResult compute_sdf(const Mesh& rSurface, const SdfOptions& rOptions) {
    (void)rSurface;
    (void)rOptions;
    // The types are final from v9.24.0 precisely so that adding this body later
    // is a pure .cpp change; see operations/sdf.hpp on why the layout ships
    // ahead of the implementation.
    throw std::invalid_argument(
        std::string(kSdfPrefix) +
        "compute_sdf is not implemented yet (it lands in v9.25.0); compose voxelize() with "
        "distance_to_surface() in the meantime");
}

}  // namespace meshioplusplus

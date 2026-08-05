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
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/sdf.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/grid_lattice.hpp"
#include "meshioplusplus/detail/region_remap.hpp"
#include "meshioplusplus/detail/surface_distance.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/operations/refine.hpp"
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
                            sum = detail::vec3_add(sum,
                                                   detail::read_point(points, dim, face.first[i]));
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

// The per-cell diagonal of every cell of a hexahedral grid, from the cell's own
// corner bounding box rather than from its level. Deriving it from `refine:level`
// would work for an octree built here and be wrong for anything else, and the
// bbox is exact for an axis-aligned box, which is all this grid ever contains.
std::vector<double> sdfop_cell_diagonals(const Mesh& rMesh) {
    std::vector<double> out;
    const NDArray& points = rMesh.Points();
    const std::size_t dim = rMesh.PointDim();
    for (const auto cb : rMesh.CellRange()) {
        const std::size_t ncells = cb.NumCells();
        const std::size_t base = out.size();
        out.resize(base + ncells);
        if (cb.IsRagged()) {
            // An octree pass never produces one (refine refuses ragged blocks),
            // so a zero diagonal here simply means "never in the band".
            for (std::size_t c = 0; c < ncells; ++c)
                out[base + c] = 0.0;
            continue;
        }
        const NDArray& conn = cb.Conn();
        const std::size_t npc = cb.NodesPerCell();
        parallel_for_bw(ncells, [&](std::size_t c) {
            detail::Vec3 lo{{0.0, 0.0, 0.0}};
            detail::Vec3 hi{{0.0, 0.0, 0.0}};
            for (std::size_t i = 0; i < npc; ++i) {
                const detail::Vec3 p =
                    detail::read_point(points, dim, detail::read_int(conn, c * npc + i));
                for (std::size_t d = 0; d < 3; ++d) {
                    if (i == 0 || p[d] < lo[d])
                        lo[d] = p[d];
                    if (i == 0 || p[d] > hi[d])
                        hi[d] = p[d];
                }
            }
            double s = 0.0;
            for (std::size_t d = 0; d < 3; ++d) {
                const double e = hi[d] - lo[d];
                s += e * e;
            }
            out[base + c] = std::sqrt(s);
        });
    }
    return out;
}

// A Float64 field_data entry from a small fixed run of doubles.
NDArray sdfop_field_f64(const double* pValues, std::size_t Count) {
    NDArray a = NDArray::Uninit(DType::Float64, {Count});
    std::memcpy(a.Data(), pValues, Count * sizeof(double));
    return a;
}

// An Int64 field_data entry from a small fixed run.
NDArray sdfop_field_i64(const std::int64_t* pValues, std::size_t Count) {
    NDArray a = NDArray::Uninit(DType::Int64, {Count});
    std::memcpy(a.Data(), pValues, Count * sizeof(std::int64_t));
    return a;
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
    const std::vector<detail::DistanceHit> hits = detail::query_distances(query, queries, rOptions);
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
        parallel_for_bw(n, [&](std::size_t p) { pi[p] = hits[p].mSignedDistance < 0.0 ? 1 : 0; });
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
    // The grid the field will live on. Voxel resolves the caller's resolution or
    // cell size directly; octree resolves a *root* lattice of mRootResolution
    // cells per axis over the same padded box and then refines it.
    detail::LatticeRequest req;
    req.mBounds = rOptions.mBounds;
    req.mPadding = rOptions.mPadding;
    req.mPaddingRelative = rOptions.mPaddingRelative;
    req.mMaxCells = rOptions.mMaxCells;
    if (rOptions.mStructure == SdfStructure::Octree) {
        if (rOptions.mRootResolution <= 0)
            throw std::invalid_argument(std::string(kSdfPrefix) +
                                        "root_resolution must be positive, got " +
                                        std::to_string(rOptions.mRootResolution));
        if (rOptions.mMaxDepth < 0)
            throw std::invalid_argument(std::string(kSdfPrefix) +
                                        "max_depth must not be negative, got " +
                                        std::to_string(rOptions.mMaxDepth));
        if (!(rOptions.mBandCells > 0.0))
            throw std::invalid_argument(std::string(kSdfPrefix) +
                                        "band_cells must be positive, got " +
                                        std::to_string(rOptions.mBandCells));
        if (rOptions.mResolution.has_value() || rOptions.mCellSize.has_value())
            throw std::invalid_argument(
                std::string(kSdfPrefix) +
                "structure 'octree' sizes itself from root_resolution and max_depth; "
                "resolution and cell_size apply to structure 'voxel' only");
        req.mResolution = std::array<std::int64_t, 3>{
            {rOptions.mRootResolution, rOptions.mRootResolution, rOptions.mRootResolution}};
    } else {
        req.mResolution = rOptions.mResolution;
        req.mCellSize = rOptions.mCellSize;
    }
    const detail::LatticeSpec root = detail::lattice_resolve(rSurface, req, kSdfPrefix);

    SdfResult result;
    result.mDims = root.mDims;
    result.mOrigin = root.mOrigin;
    result.mSpacing = root.mSpacing;

    // Every point and cell of the output is new, so nothing can be remapped --
    // the same situation voxelize/slice/isosurface/extract_surface are in.
    detail::warn_regions_dropped(rSurface, "compute_sdf");

    Mesh mesh = detail::lattice_build_mesh(root);

    if (rOptions.mStructure == SdfStructure::Octree) {
        // The soup and its accelerator depend only on the surface, so they are
        // built once and reused by every selection pass.
        const detail::TriangleSoup soup =
            detail::build_triangle_soup(rSurface, rOptions.mDistance.mSurfaceRegion);
        // Selection needs |d| only, so it runs UNSIGNED: the magnitude is the
        // same either way, and an unsigned query neither depends on the surface
        // being closed nor pays for a winding-number sign it would throw away.
        SurfaceDistanceOptions sel_opts = rOptions.mDistance;
        sel_opts.mSign = SdfSign::Unsigned;
        sel_opts.mBand = 0.0;
        sel_opts.mWatertightCheck = SdfWatertightCheck::Off;  // reported once, below
        const detail::DistanceQuery query = detail::build_distance_query(soup, sel_opts);

        for (std::int64_t pass = 0; pass < rOptions.mMaxDepth; ++pass) {
            // Distances at the CURRENT mesh's cell centres. Recomputing per pass
            // is the point: a selection carried over from the previous pass would
            // name cells of a mesh that no longer exists, since these are global
            // block-major indices into whatever `refine` last returned.
            const std::vector<detail::Vec3> centres = sdfop_query_points(mesh, SdfLocation::Center);
            const std::vector<double> diag = sdfop_cell_diagonals(mesh);
            const std::vector<detail::DistanceHit> hits =
                detail::query_distances(query, centres, sel_opts);

            std::vector<std::int64_t> selected;
            for (std::size_t c = 0; c < hits.size(); ++c)
                if (std::fabs(hits[c].mSignedDistance) <= rOptions.mBandCells * diag[c])
                    selected.push_back(static_cast<std::int64_t>(c));
            if (selected.empty())
                break;  // nothing near the surface: the tree has converged

            RefineOptions ro;
            ro.mLevels = 1;
            ro.mCells = std::move(selected);
            // Balanced is the only closure whose cost is bounded by the
            // selection; red-green and propagate both spread the refinement well
            // beyond the band. The output is 1-irregular -- see refine.hpp.
            ro.mClosure = RefineClosure::Balanced;
            ro.mRecordLevels = rOptions.mRecordLevels;
            mesh = refine(mesh, ro).mMesh;
            ++result.mMaxDepth;
            for (std::size_t d = 0; d < 3; ++d)
                result.mSpacing[d] *= 0.5;

            std::size_t total = 0;
            for (const auto cb : mesh.CellRange())
                total += cb.NumCells();
            if (rOptions.mMaxCells > 0 && static_cast<std::int64_t>(total) > rOptions.mMaxCells)
                throw std::invalid_argument(
                    std::string(kSdfPrefix) + "the octree reached " + std::to_string(total) +
                    " cells after " + std::to_string(result.mMaxDepth) +
                    " pass(es), above the limit of " + std::to_string(rOptions.mMaxCells) +
                    " (raise max_cells, lower max_depth or band_cells)");
        }
    }

    // The field is computed ONCE, on the final mesh. Attaching it earlier and
    // letting `refine` carry it would replace every value on a refined cell with
    // an interpolation of its parent's -- smooth, plausible, and wrong.
    SurfaceDistanceResult filled = distance_to_surface(mesh, rSurface, rOptions.mDistance);
    result.mMesh = std::move(filled.mMesh);
    result.mQuality = filled.mQuality;
    result.mNumBanded = filled.mNumBanded;

    // The header, so a caller holding only the mesh still knows what grid it is.
    // No format persists arbitrary field_data, so this survives in memory only --
    // `.vti` round-trips the geometry instead, and `detail::lattice_from_mesh`
    // recovers a dense lattice from its points. See doc/sdf.md.
    const double origin[3] = {result.mOrigin[0], result.mOrigin[1], result.mOrigin[2]};
    const double spacing[3] = {result.mSpacing[0], result.mSpacing[1], result.mSpacing[2]};
    const std::int64_t dims[3] = {result.mDims[0], result.mDims[1], result.mDims[2]};
    const double bounds[6] = {
        root.mOrigin[0],
        root.mOrigin[1],
        root.mOrigin[2],
        root.mOrigin[0] + static_cast<double>(root.mDims[0]) * root.mSpacing[0],
        root.mOrigin[1] + static_cast<double>(root.mDims[1]) * root.mSpacing[1],
        root.mOrigin[2] + static_cast<double>(root.mDims[2]) * root.mSpacing[2]};
    const std::int64_t depth = result.mMaxDepth;
    const std::int64_t structure = rOptions.mStructure == SdfStructure::Octree ? 1 : 0;
    result.mMesh.AddFieldData(kSdfOriginName, sdfop_field_f64(origin, 3));
    result.mMesh.AddFieldData(kSdfSpacingName, sdfop_field_f64(spacing, 3));
    result.mMesh.AddFieldData(kSdfDimsName, sdfop_field_i64(dims, 3));
    result.mMesh.AddFieldData(kSdfBoundsName, sdfop_field_f64(bounds, 6));
    result.mMesh.AddFieldData(kSdfMaxDepthName, sdfop_field_i64(&depth, 1));
    result.mMesh.AddFieldData(kSdfStructureName, sdfop_field_i64(&structure, 1));
    return result;
}

}  // namespace meshioplusplus

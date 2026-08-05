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
// Regular grids: `grid` from nothing, `voxelize` around a mesh. See
// operations/voxelize.hpp for the contract and detail/grid_lattice.hpp for the
// numbering, which this file does not re-derive.
//
// Anonymous-namespace helpers are prefixed `vox_` (the amalgamation concatenates
// every translation unit into one). Note `grid_` was NOT available: it belongs to
// detail/spatial_hash.hpp's grid_quantize, which is a bucket hash rather than a
// mesh, and reusing it would make two unrelated things look related.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/voxelize.hpp"
#include "meshioplusplus/detail/grid_lattice.hpp"
#include "meshioplusplus/detail/surface_distance.hpp"
#include "meshioplusplus/detail/tri_box.hpp"
#include "meshioplusplus/detail/region_remap.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {

namespace {

constexpr const char* kVoxPrefix = "meshio++: voxelize: ";

// The cell budget check for grid(), which has no LatticeRequest to resolve. It
// is a *named* refusal rather than an allocation failure because the failure
// mode it guards is mundane and predictable: 512^3 is 134 million cells and
// ~11.8 GB of points and connectivity, which a user asks for by typo far more
// often than on purpose. `lattice_resolve` applies the identical rule (and the
// identical message) to the request-shaped callers.
void vox_check_budget(std::int64_t Cells, std::int64_t MaxCells, const char* pPrefix) {
    if (MaxCells <= 0 || Cells <= MaxCells)
        return;
    throw std::invalid_argument(std::string(pPrefix) + "the requested grid has " +
                                std::to_string(Cells) + " cells, above the limit of " +
                                std::to_string(MaxCells) +
                                " (raise max_cells, coarsen the resolution, or use a band)");
}

// The lattice a voxelization runs on. Exactly one of resolution and cell size
// must be given: defaulting one of them would silently pick a grid the caller
// did not choose, and for an object whose cost is cubic in that choice that is
// not a kindness. The rule itself lives in detail/grid_lattice.hpp, since
// `compute_sdf` resolves the identical six fields the identical way.
detail::LatticeSpec vox_resolve_lattice(const Mesh& rMesh, const VoxelOptions& rOptions) {
    detail::LatticeRequest req;
    req.mResolution = rOptions.mResolution;
    req.mCellSize = rOptions.mCellSize;
    req.mBounds = rOptions.mBounds;
    req.mPadding = rOptions.mPadding;
    req.mPaddingRelative = rOptions.mPaddingRelative;
    req.mMaxCells = rOptions.mMaxCells;
    return detail::lattice_resolve(rMesh, req, kVoxPrefix);
}

// Which cells the fill rule keeps, as a per-cell flag. An empty result means
// "all of them", which lets the caller skip the subsetting pass entirely.
std::vector<char> vox_occupancy(const Mesh& rMesh, const detail::LatticeSpec& rSpec,
                                const VoxelOptions& rOptions) {
    if (rOptions.mFill == VoxelFill::All)
        return {};

    const std::int64_t ncells = detail::lattice_num_cells(rSpec);
    std::vector<char> occupied(static_cast<std::size_t>(ncells), 0);
    const std::int64_t nx = rSpec.mDims[0];
    const std::int64_t ny = rSpec.mDims[1];

    if (rOptions.mFill == VoxelFill::Surface) {
        // Exact triangle/box overlap, not a bounding-box test: a long diagonal
        // triangle overlaps far more boxes than it enters.
        const detail::TriangleSoup soup =
            detail::build_triangle_soup(rMesh, rOptions.mDistance.mSurfaceRegion);
        const detail::Vec3 half{
            {rSpec.mSpacing[0] * 0.5, rSpec.mSpacing[1] * 0.5, rSpec.mSpacing[2] * 0.5}};
        // Serial over triangles, each touching only the voxels of its own
        // bounding box: a parallel pass would race on the shared flag array for
        // no real gain, since the inner box is small.
        for (std::size_t t = 0; t < soup.NumTriangles(); ++t) {
            const detail::Vec3& a = soup.mCorners[t * 3 + 0];
            const detail::Vec3& b = soup.mCorners[t * 3 + 1];
            const detail::Vec3& c = soup.mCorners[t * 3 + 2];
            std::array<std::int64_t, 3> lo{}, hi{};
            bool skip = false;
            for (std::size_t k = 0; k < 3; ++k) {
                if (!(rSpec.mSpacing[k] > 0.0)) {
                    skip = true;
                    break;
                }
                const double tlo = std::min(a[k], std::min(b[k], c[k]));
                const double thi = std::max(a[k], std::max(b[k], c[k]));
                double flo = std::floor((tlo - rSpec.mOrigin[k]) / rSpec.mSpacing[k]);
                double fhi = std::floor((thi - rSpec.mOrigin[k]) / rSpec.mSpacing[k]);
                flo = flo < 0.0 ? 0.0 : flo;
                const double last = static_cast<double>(rSpec.mDims[k] - 1);
                fhi = fhi > last ? last : fhi;
                if (fhi < flo) {
                    skip = true;
                    break;
                }
                lo[k] = static_cast<std::int64_t>(flo);
                hi[k] = static_cast<std::int64_t>(fhi);
            }
            if (skip)
                continue;
            for (std::int64_t k = lo[2]; k <= hi[2]; ++k)
                for (std::int64_t j = lo[1]; j <= hi[1]; ++j)
                    for (std::int64_t i = lo[0]; i <= hi[0]; ++i) {
                        const std::int64_t cid = (k * ny + j) * nx + i;
                        if (occupied[static_cast<std::size_t>(cid)])
                            continue;
                        const detail::Vec3 centre{
                            {rSpec.mOrigin[0] + (static_cast<double>(i) + 0.5) * rSpec.mSpacing[0],
                             rSpec.mOrigin[1] + (static_cast<double>(j) + 0.5) * rSpec.mSpacing[1],
                             rSpec.mOrigin[2] +
                                 (static_cast<double>(k) + 0.5) * rSpec.mSpacing[2]}};
                        if (detail::tri_box_overlap(centre, half, a, b, c))
                            occupied[static_cast<std::size_t>(cid)] = 1;
                    }
        }
        return occupied;
    }

    // Inside: sign the distance at each cell centre. This is the only fill that
    // needs a closed surface, and it says so through the watertight check rather
    // than by quietly producing a hollow result.
    const detail::TriangleSoup soup =
        detail::build_triangle_soup(rMesh, rOptions.mDistance.mSurfaceRegion);
    std::vector<detail::Vec3> centres(static_cast<std::size_t>(ncells));
    parallel_for_bw(static_cast<std::size_t>(ncells), [&](std::size_t c) {
        const std::int64_t g = static_cast<std::int64_t>(c);
        const std::int64_t i = g % nx;
        const std::int64_t j = (g / nx) % ny;
        const std::int64_t k = g / (nx * ny);
        centres[c] =
            detail::Vec3{{rSpec.mOrigin[0] + (static_cast<double>(i) + 0.5) * rSpec.mSpacing[0],
                          rSpec.mOrigin[1] + (static_cast<double>(j) + 0.5) * rSpec.mSpacing[1],
                          rSpec.mOrigin[2] + (static_cast<double>(k) + 0.5) * rSpec.mSpacing[2]}};
    });

    SurfaceDistanceOptions dopts = rOptions.mDistance;
    if (dopts.mSign == SdfSign::Unsigned)
        throw std::invalid_argument(std::string(kVoxPrefix) +
                                    "fill 'inside' needs a sign, but sign='unsigned' was given");
    dopts.mBand = 0.0;  // a band would clamp the very sign this fill depends on
    const detail::DistanceQuery query = detail::build_distance_query(soup, dopts);
    const std::vector<detail::DistanceHit> hits = detail::query_distances(query, centres, dopts);
    parallel_for_bw(static_cast<std::size_t>(ncells),
                    [&](std::size_t c) { occupied[c] = hits[c].mSignedDistance < 0.0 ? 1 : 0; });
    return occupied;
}

// Emit only the flagged cells, compacting the points they reference in ascending
// order so the output does not depend on traversal.
Mesh vox_build_subset(const detail::LatticeSpec& rSpec, const std::vector<char>& rOccupied,
                      std::int64_t Kept) {
    Mesh out;
    const std::int64_t nx = rSpec.mDims[0];
    const std::int64_t ny = rSpec.mDims[1];
    const std::int64_t px = nx + 1;
    const std::int64_t py = ny + 1;
    const std::int64_t npoints = detail::lattice_num_points(rSpec);
    const std::int64_t ncells = detail::lattice_num_cells(rSpec);

    // Two passes: mark which points the kept cells reference, then number them in
    // ASCENDING original order. That is surface.cpp's used/remap pattern, and the
    // ascending part is load-bearing rather than incidental -- numbering on first
    // encounter instead would make the point ids depend on the hexahedron's node
    // order, which is a traversal detail no caller should be able to observe (and
    // which the numpy twin would then have to replicate rather than simply sort).
    std::vector<char> used(static_cast<std::size_t>(npoints), 0);
    std::vector<std::int64_t> kept_cells;
    kept_cells.reserve(static_cast<std::size_t>(Kept));
    for (std::int64_t c = 0; c < ncells; ++c) {
        if (!rOccupied[static_cast<std::size_t>(c)])
            continue;
        kept_cells.push_back(c);
        const std::int64_t i = c % nx;
        const std::int64_t j = (c / nx) % ny;
        const std::int64_t k = c / (nx * ny);
        const std::int64_t base = (k * py + j) * px + i;
        const std::int64_t top = base + px * py;
        const std::int64_t nodes[8] = {base, base + 1, base + px + 1, base + px,
                                       top,  top + 1,  top + px + 1,  top + px};
        for (std::int64_t nd : nodes)
            used[static_cast<std::size_t>(nd)] = 1;
    }

    std::vector<std::int64_t> remap(static_cast<std::size_t>(npoints), -1);
    std::int64_t next = 0;
    for (std::int64_t p = 0; p < npoints; ++p)
        if (used[static_cast<std::size_t>(p)])
            remap[static_cast<std::size_t>(p)] = next++;

    std::vector<std::int64_t> conn;
    conn.reserve(static_cast<std::size_t>(Kept) * 8);
    for (std::int64_t c : kept_cells) {
        const std::int64_t i = c % nx;
        const std::int64_t j = (c / nx) % ny;
        const std::int64_t k = c / (nx * ny);
        const std::int64_t base = (k * py + j) * px + i;
        const std::int64_t top = base + px * py;
        const std::int64_t nodes[8] = {base, base + 1, base + px + 1, base + px,
                                       top,  top + 1,  top + px + 1,  top + px};
        for (std::int64_t nd : nodes)
            conn.push_back(remap[static_cast<std::size_t>(nd)]);
    }

    NDArray points =
        NDArray::Uninit(DType::Float64, {static_cast<std::size_t>(next), std::size_t{3}});
    double* pdst = points.As<double>();
    parallel_for_bw(static_cast<std::size_t>(npoints), [&](std::size_t p) {
        const std::int64_t slot = remap[p];
        if (slot < 0)
            return;
        const std::int64_t g = static_cast<std::int64_t>(p);
        const std::int64_t i = g % px;
        const std::int64_t j = (g / px) % py;
        const std::int64_t k = g / (px * py);
        pdst[slot * 3 + 0] = rSpec.mOrigin[0] + static_cast<double>(i) * rSpec.mSpacing[0];
        pdst[slot * 3 + 1] = rSpec.mOrigin[1] + static_cast<double>(j) * rSpec.mSpacing[1];
        pdst[slot * 3 + 2] = rSpec.mOrigin[2] + static_cast<double>(k) * rSpec.mSpacing[2];
    });
    out.AssignPoints(std::move(points));

    if (Kept > 0) {
        NDArray block =
            NDArray::Uninit(DType::Int64, {static_cast<std::size_t>(Kept), std::size_t{8}});
        std::memcpy(block.Data(), conn.data(), conn.size() * sizeof(std::int64_t));
        out.AddCellBlock("hexahedron", std::move(block));
    }
    return out;
}

}  // namespace

VoxelFill voxel_fill_from_name(const std::string& rName) {
    if (rName == "all")
        return VoxelFill::All;
    if (rName == "surface")
        return VoxelFill::Surface;
    if (rName == "inside")
        return VoxelFill::Inside;
    throw std::invalid_argument(std::string(kVoxPrefix) + "unknown fill '" + rName +
                                "' (expected one of: all, surface, inside)");
}

Mesh grid(const std::array<std::int64_t, 3>& rDims, const std::array<double, 3>& rOrigin,
          const std::array<double, 3>& rSpacing, std::int64_t MaxCells) {
    constexpr const char* prefix = "meshio++: grid: ";
    for (std::size_t k = 0; k < 3; ++k) {
        if (rDims[k] < 0)
            throw std::invalid_argument(std::string(prefix) +
                                        "cell counts must not be negative, got " +
                                        std::to_string(rDims[k]) + " on axis " + std::to_string(k));
        if (rDims[k] > 0 && !(rSpacing[k] > 0.0))
            throw std::invalid_argument(
                std::string(prefix) + "spacing must be positive on every axis with cells, got " +
                std::to_string(rSpacing[k]) + " on axis " + std::to_string(k));
    }
    detail::LatticeSpec spec;
    spec.mOrigin = rOrigin;
    spec.mSpacing = rSpacing;
    spec.mDims = rDims;
    vox_check_budget(detail::lattice_num_cells(spec), MaxCells, prefix);
    return detail::lattice_build_mesh(spec);
}

VoxelResult voxelize(const Mesh& rMesh, const VoxelOptions& rOptions) {
    const detail::LatticeSpec spec = vox_resolve_lattice(rMesh, rOptions);

    VoxelResult result;
    result.mDims = spec.mDims;
    result.mOrigin = spec.mOrigin;
    result.mSpacing = spec.mSpacing;

    // Every point and cell of the output is new, so nothing can be remapped --
    // the same situation slice/isosurface/extract_surface are in.
    detail::warn_regions_dropped(rMesh, "voxelize");

    const std::int64_t ncells = detail::lattice_num_cells(spec);
    const std::vector<char> occupied = vox_occupancy(rMesh, spec, rOptions);

    if (rOptions.mFill == VoxelFill::All || occupied.empty()) {
        result.mMesh = detail::lattice_build_mesh(spec);
        result.mNumOccupied = ncells;
        if (rOptions.mAttachOccupancy && ncells > 0) {
            // Every cell is kept, so occupancy is all ones -- constant, but
            // emitted anyway so that a caller can switch fill modes without
            // their downstream pipeline discovering the array has vanished.
            NDArray occ = NDArray::Uninit(DType::Int64, {static_cast<std::size_t>(ncells)});
            std::int64_t* dst = occ.As<std::int64_t>();
            parallel_for_bw(static_cast<std::size_t>(ncells), [&](std::size_t c) { dst[c] = 1; });
            std::vector<NDArray> blocks;
            blocks.push_back(std::move(occ));
            result.mMesh.AddCellData(kVoxelOccupancyName, std::move(blocks));
        }
        return result;
    }

    // A selective fill keeps a subset of the cells. Building the full lattice and
    // then subsetting it would cost the memory the budget check exists to avoid,
    // so the kept cells are emitted directly and the points they reference are
    // compacted in ascending order -- surface.cpp's used/remap pattern, which is
    // what makes the output independent of which thread found which cell.
    std::int64_t nkept = 0;
    for (std::int64_t c = 0; c < ncells; ++c)
        nkept += occupied[static_cast<std::size_t>(c)] ? 1 : 0;
    result.mNumOccupied = nkept;
    result.mMesh = vox_build_subset(spec, occupied, nkept);
    if (rOptions.mAttachOccupancy && nkept > 0) {
        NDArray occ = NDArray::Uninit(DType::Int64, {static_cast<std::size_t>(nkept)});
        std::int64_t* dst = occ.As<std::int64_t>();
        parallel_for_bw(static_cast<std::size_t>(nkept), [&](std::size_t c) { dst[c] = 1; });
        std::vector<NDArray> blocks;
        blocks.push_back(std::move(occ));
        result.mMesh.AddCellData(kVoxelOccupancyName, std::move(blocks));
    }
    return result;
}

}  // namespace meshioplusplus

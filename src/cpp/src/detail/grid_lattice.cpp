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
// The regular hexahedron lattice. See detail/grid_lattice.hpp for the numbering
// contract and why it is inherited rather than chosen.
//
// Everything here is index arithmetic, so there is no determinism argument to
// make beyond "do not accumulate": a point's coordinate is origin + index *
// spacing, evaluated independently per point, which is what lets the numpy twin
// be bit-identical without replicating a traversal order.

// System includes
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/detail/grid_lattice.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

// The number of cells an axis needs to cover `extent` with cells of at most
// `cell`. A non-positive extent needs none; a non-positive cell size is the
// caller's error and is reported by lattice_from_cell_size, not here.
std::int64_t lat_axis_count(double extent, double cell) {
    if (!(extent > 0.0) || !(cell > 0.0))
        return 0;
    const double n = std::ceil(extent / cell);
    // A hostile ratio (a denormal cell size against a huge extent) would
    // overflow the cast; clamping here keeps the error a named one from the
    // caller's own budget check rather than undefined behaviour.
    if (!(n < 9.0e18))
        return std::numeric_limits<std::int64_t>::max();
    return static_cast<std::int64_t>(n);
}

}  // namespace

std::int64_t lattice_num_points(const LatticeSpec& rSpec) {
    if (rSpec.mDims[0] < 0 || rSpec.mDims[1] < 0 || rSpec.mDims[2] < 0)
        return 0;
    return (rSpec.mDims[0] + 1) * (rSpec.mDims[1] + 1) * (rSpec.mDims[2] + 1);
}

std::int64_t lattice_num_cells(const LatticeSpec& rSpec) {
    if (rSpec.mDims[0] <= 0 || rSpec.mDims[1] <= 0 || rSpec.mDims[2] <= 0)
        return 0;
    return rSpec.mDims[0] * rSpec.mDims[1] * rSpec.mDims[2];
}

LatticeSpec lattice_from_bounds(const std::array<double, 3>& rLo, const std::array<double, 3>& rHi,
                                const std::array<std::int64_t, 3>& rDims) {
    LatticeSpec spec;
    spec.mOrigin = rLo;
    spec.mDims = rDims;
    for (std::size_t k = 0; k < 3; ++k) {
        const std::int64_t n = rDims[k] < 0 ? 0 : rDims[k];
        spec.mDims[k] = n;
        spec.mSpacing[k] = n > 0 ? (rHi[k] - rLo[k]) / static_cast<double>(n) : 0.0;
    }
    return spec;
}

LatticeSpec lattice_from_cell_size(const std::array<double, 3>& rLo,
                                   const std::array<double, 3>& rHi,
                                   const std::array<double, 3>& rCellSize) {
    LatticeSpec spec;
    spec.mOrigin = rLo;
    for (std::size_t k = 0; k < 3; ++k) {
        spec.mDims[k] = lat_axis_count(rHi[k] - rLo[k], rCellSize[k]);
        // The requested size is honoured exactly and the box grows to fit, so
        // the lattice covers the input rather than clipping it. Recomputing the
        // spacing from the (possibly larger) extent instead would silently hand
        // back cells of a different size than were asked for.
        spec.mSpacing[k] = spec.mDims[k] > 0 ? rCellSize[k] : 0.0;
    }
    return spec;
}

Mesh lattice_build_mesh(const LatticeSpec& rSpec) {
    Mesh out;
    const std::int64_t ncells = lattice_num_cells(rSpec);
    if (ncells <= 0) {
        // No cells means no lattice. Emitting an empty points array and no block
        // at all is what lets a caller test NumCellBlocks() rather than a row
        // count.
        out.AssignPoints(NDArray(DType::Float64, {std::size_t{0}, std::size_t{3}}));
        return out;
    }

    const std::int64_t nx = rSpec.mDims[0];
    const std::int64_t ny = rSpec.mDims[1];
    const std::int64_t nz = rSpec.mDims[2];
    const std::int64_t px = nx + 1;
    const std::int64_t py = ny + 1;
    const std::int64_t npoints = lattice_num_points(rSpec);

    NDArray points =
        NDArray::Uninit(DType::Float64, {static_cast<std::size_t>(npoints), std::size_t{3}});
    {
        double* dst = points.As<double>();
        const double ox = rSpec.mOrigin[0], oy = rSpec.mOrigin[1], oz = rSpec.mOrigin[2];
        const double hx = rSpec.mSpacing[0], hy = rSpec.mSpacing[1], hz = rSpec.mSpacing[2];
        parallel_for_bw(static_cast<std::size_t>(npoints), [&](std::size_t p) {
            const std::int64_t g = static_cast<std::int64_t>(p);
            const std::int64_t i = g % px;
            const std::int64_t j = (g / px) % py;
            const std::int64_t k = g / (px * py);
            // origin + index * spacing, never an accumulation: the coordinate of
            // a point must not depend on how many points precede it.
            dst[p * 3 + 0] = ox + static_cast<double>(i) * hx;
            dst[p * 3 + 1] = oy + static_cast<double>(j) * hy;
            dst[p * 3 + 2] = oz + static_cast<double>(k) * hz;
        });
    }
    out.AssignPoints(std::move(points));

    NDArray conn =
        NDArray::Uninit(DType::Int64, {static_cast<std::size_t>(ncells), std::size_t{8}});
    {
        std::int64_t* dst = conn.As<std::int64_t>();
        parallel_for_bw(static_cast<std::size_t>(ncells), [&](std::size_t c) {
            const std::int64_t g = static_cast<std::int64_t>(c);
            const std::int64_t i = g % nx;
            const std::int64_t j = (g / nx) % ny;
            const std::int64_t k = g / (nx * ny);
            const std::int64_t base = (k * py + j) * px + i;
            const std::int64_t top = base + px * py;
            // The meshio/VTK hexahedron winding: bottom ring counter-clockwise,
            // then the top ring, so the base normal points at the top face.
            dst[c * 8 + 0] = base;
            dst[c * 8 + 1] = base + 1;
            dst[c * 8 + 2] = base + px + 1;
            dst[c * 8 + 3] = base + px;
            dst[c * 8 + 4] = top;
            dst[c * 8 + 5] = top + 1;
            dst[c * 8 + 6] = top + px + 1;
            dst[c * 8 + 7] = top + px;
        });
    }
    out.AddCellBlock("hexahedron", std::move(conn));
    return out;
}

LatticeSpec lattice_resolve(const Mesh& rMesh, const LatticeRequest& rRequest,
                            const char* pPrefix) {
    const std::string prefix(pPrefix);
    if (rRequest.mResolution.has_value() == rRequest.mCellSize.has_value())
        throw std::invalid_argument(prefix + "give exactly one of resolution and cell_size");

    // The box: the caller's, or the mesh's own, in both cases grown by the
    // padding.
    std::array<double, 3> lo{{0.0, 0.0, 0.0}};
    std::array<double, 3> hi{{0.0, 0.0, 0.0}};
    if (rRequest.mBounds.has_value()) {
        const std::array<double, 6>& b = *rRequest.mBounds;
        for (std::size_t k = 0; k < 3; ++k) {
            lo[k] = b[k];
            hi[k] = b[k + 3];
            if (!(hi[k] >= lo[k]))
                throw std::invalid_argument(prefix + "bounds are inverted on axis " +
                                            std::to_string(k) + " (lo " + std::to_string(lo[k]) +
                                            " > hi " + std::to_string(hi[k]) + ")");
        }
    } else if (!point_bbox(rMesh, lo, hi)) {
        throw std::invalid_argument(prefix +
                                    "the mesh has no points, so it has no bounding box to "
                                    "cover (pass explicit bounds)");
    }

    double diag = 0.0;
    for (std::size_t k = 0; k < 3; ++k) {
        const double e = hi[k] - lo[k];
        diag += e * e;
    }
    diag = std::sqrt(diag);
    const double pad = rRequest.mPadding + rRequest.mPaddingRelative * diag;
    if (pad < 0.0)
        throw std::invalid_argument(prefix + "padding is negative");
    for (std::size_t k = 0; k < 3; ++k) {
        lo[k] -= pad;
        hi[k] += pad;
    }

    LatticeSpec spec;
    if (rRequest.mResolution.has_value()) {
        const std::array<std::int64_t, 3>& r = *rRequest.mResolution;
        for (std::size_t k = 0; k < 3; ++k)
            if (r[k] <= 0)
                throw std::invalid_argument(prefix +
                                            "resolution must be positive on every axis, got " +
                                            std::to_string(r[k]) + " on axis " + std::to_string(k));
        spec = lattice_from_bounds(lo, hi, r);
    } else {
        const double cell = *rRequest.mCellSize;
        if (!(cell > 0.0))
            throw std::invalid_argument(prefix + "cell_size must be positive");
        spec = lattice_from_cell_size(lo, hi, {{cell, cell, cell}});
        for (std::size_t k = 0; k < 3; ++k)
            if (spec.mDims[k] <= 0)
                throw std::invalid_argument(
                    prefix + "the bounding box is degenerate on axis " + std::to_string(k) +
                    ", so a cell size cannot fill it (pass an explicit resolution or bounds)");
    }

    const std::int64_t cells = lattice_num_cells(spec);
    if (rRequest.mMaxCells > 0 && cells > rRequest.mMaxCells)
        throw std::invalid_argument(prefix + "the requested grid has " + std::to_string(cells) +
                                    " cells, above the limit of " +
                                    std::to_string(rRequest.mMaxCells) +
                                    " (raise max_cells, coarsen the resolution, or use a band)");
    return spec;
}

bool lattice_from_mesh(const Mesh& rMesh, LatticeSpec& rSpec) {
    // Exactly one hexahedron block, and nothing else.
    if (rMesh.NumCellBlocks() != 1)
        return false;
    const auto cb = rMesh.Cells(0);
    if (cb.Type() != std::string("hexahedron") || cb.IsRagged())
        return false;

    const std::size_t npoints = rMesh.NumPoints();
    const std::size_t dim = rMesh.PointDim();
    if (npoints == 0 || dim < 3)
        return false;
    const NDArray& points = rMesh.Points();

    // The distinct plane positions per axis. An exact sort-and-unique is correct
    // here rather than merely convenient: lattice_build_mesh evaluates
    // `origin + index * spacing` independently per point, so every point on a
    // given plane carries the identical double.
    std::array<std::vector<double>, 3> planes;
    for (std::size_t d = 0; d < 3; ++d) {
        planes[d].resize(npoints);
        for (std::size_t p = 0; p < npoints; ++p)
            planes[d][p] = read_double(points, p * dim + d);
        std::sort(planes[d].begin(), planes[d].end());
        planes[d].erase(std::unique(planes[d].begin(), planes[d].end()), planes[d].end());
        if (planes[d].size() < 2)
            return false;  // a single plane is a sheet, not a lattice
    }

    // A dense lattice has exactly the product of its plane counts as points, and
    // the product of its cell counts as cells. Either mismatch means a subset (a
    // `surface`/`inside` fill, an octree) or something else entirely.
    if (planes[0].size() * planes[1].size() * planes[2].size() != npoints)
        return false;
    const std::int64_t nx = static_cast<std::int64_t>(planes[0].size()) - 1;
    const std::int64_t ny = static_cast<std::int64_t>(planes[1].size()) - 1;
    const std::int64_t nz = static_cast<std::int64_t>(planes[2].size()) - 1;
    if (static_cast<std::size_t>(nx * ny * nz) != cb.NumCells())
        return false;

    LatticeSpec spec;
    spec.mDims = {{nx, ny, nz}};
    for (std::size_t d = 0; d < 3; ++d) {
        const std::vector<double>& v = planes[d];
        const std::int64_t n = static_cast<std::int64_t>(v.size()) - 1;
        spec.mOrigin[d] = v.front();
        // (last - first) / n, not the first gap: it is the least-error estimate
        // and it is what the writer's own `origin + index * spacing` inverts.
        spec.mSpacing[d] = (v.back() - v.front()) / static_cast<double>(n);
        if (!(spec.mSpacing[d] > 0.0))
            return false;
        // Uniformity to a relative 1e-9. `origin + i*h` gaps differ in the last
        // bits, so this cannot be exact, but a genuinely graded mesh is rejected
        // by orders of magnitude rather than by ulps.
        for (std::int64_t i = 0; i < n; ++i) {
            const double gap = v[static_cast<std::size_t>(i) + 1] - v[static_cast<std::size_t>(i)];
            if (std::fabs(gap - spec.mSpacing[d]) > 1.0e-9 * spec.mSpacing[d])
                return false;
        }
    }

    const std::int64_t px = nx + 1;
    const std::int64_t py = ny + 1;

    // The points must be in the lattice's own x-fastest order, not merely occupy
    // its plane positions -- a permuted grid has identical plane sets and is a
    // different mesh. The comparison is EXACT because `planes[d]` holds the very
    // doubles the points carry, so this is an equality test rather than a fit.
    for (std::size_t p = 0; p < npoints; ++p) {
        const std::int64_t g = static_cast<std::int64_t>(p);
        const std::int64_t idx[3] = {g % px, (g / px) % py, g / (px * py)};
        for (std::size_t d = 0; d < 3; ++d)
            if (read_double(points, p * dim + d) != planes[d][static_cast<std::size_t>(idx[d])])
                return false;
    }

    // The connectivity must be the lattice's own, or two different meshes would
    // both claim to be this box. Checking the index formula per cell is O(n) and
    // is the difference between "has lattice-shaped points" and "is a lattice".
    const NDArray& conn = cb.Conn();
    if (cb.NodesPerCell() != 8)
        return false;
    for (std::int64_t c = 0; c < nx * ny * nz; ++c) {
        const std::int64_t i = c % nx;
        const std::int64_t j = (c / nx) % ny;
        const std::int64_t k = c / (nx * ny);
        const std::int64_t base = (k * py + j) * px + i;
        const std::int64_t top = base + px * py;
        const std::int64_t expect[8] = {base, base + 1, base + px + 1, base + px,
                                        top,  top + 1,  top + px + 1,  top + px};
        for (std::size_t m = 0; m < 8; ++m)
            if (read_int(conn, static_cast<std::size_t>(c) * 8 + m) != expect[m])
                return false;
    }

    rSpec = spec;
    return true;
}

bool point_bbox(const Mesh& rMesh, std::array<double, 3>& rLo, std::array<double, 3>& rHi) {
    rLo = {0.0, 0.0, 0.0};
    rHi = {0.0, 0.0, 0.0};
    const std::size_t n = rMesh.NumPoints();
    if (n == 0)
        return false;

    const NDArray& points = rMesh.Points();
    const std::size_t dim = rMesh.PointDim();
    const std::size_t ddim = dim < 3 ? dim : 3;
    if (ddim == 0)
        return false;

    // Chunked parallel reduction, then a serial combine. min/max are associative
    // and exact, so the chunking is not observable -- which is exactly why only
    // the bbox was hoisted here and stats.cpp's centroid sum was not.
    const std::size_t grain = 4096;
    const std::size_t nchunks = (n + grain - 1) / grain;
    std::vector<std::array<double, 3>> pmin(nchunks), pmax(nchunks);
    parallel_for(
        nchunks,
        [&](std::size_t ci) {
            std::array<double, 3> lmin = {std::numeric_limits<double>::infinity(),
                                          std::numeric_limits<double>::infinity(),
                                          std::numeric_limits<double>::infinity()};
            std::array<double, 3> lmax = {-std::numeric_limits<double>::infinity(),
                                          -std::numeric_limits<double>::infinity(),
                                          -std::numeric_limits<double>::infinity()};
            const std::size_t start = ci * grain;
            const std::size_t stop = n < start + grain ? n : start + grain;
            for (std::size_t g = start; g < stop; ++g)
                for (std::size_t d = 0; d < ddim; ++d) {
                    const double v = read_double(points, g * dim + d);
                    lmin[d] = lmin[d] < v ? lmin[d] : v;
                    lmax[d] = lmax[d] > v ? lmax[d] : v;
                }
            pmin[ci] = lmin;
            pmax[ci] = lmax;
        },
        1);

    std::array<double, 3> gmin = {std::numeric_limits<double>::infinity(),
                                  std::numeric_limits<double>::infinity(),
                                  std::numeric_limits<double>::infinity()};
    std::array<double, 3> gmax = {-std::numeric_limits<double>::infinity(),
                                  -std::numeric_limits<double>::infinity(),
                                  -std::numeric_limits<double>::infinity()};
    for (std::size_t ci = 0; ci < nchunks; ++ci)
        for (std::size_t d = 0; d < 3; ++d) {
            gmin[d] = gmin[d] < pmin[ci][d] ? gmin[d] : pmin[ci][d];
            gmax[d] = gmax[d] > pmax[ci][d] ? gmax[d] : pmax[ci][d];
        }
    for (std::size_t d = 0; d < ddim; ++d) {
        rLo[d] = gmin[d];
        rHi[d] = gmax[d];
    }
    return true;
}

}  // namespace detail
}  // namespace meshioplusplus

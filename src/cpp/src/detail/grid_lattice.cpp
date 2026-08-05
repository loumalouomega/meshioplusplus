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
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
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

    NDArray points = NDArray::Uninit(DType::Float64, {static_cast<std::size_t>(npoints),
                                                      std::size_t{3}});
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

    NDArray conn = NDArray::Uninit(DType::Int64, {static_cast<std::size_t>(ncells),
                                                  std::size_t{8}});
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

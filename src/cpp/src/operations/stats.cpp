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
// Geometric statistics of a mesh (read-only): bbox / extent / centroid, per-cell-
// type counts, total area (2D cells + 3D boundary), signed/unsigned volume, and
// inverted-cell count. Reuses detail::geometry and the surface operation for the
// 3D boundary. Built entirely through the uniform mesh API so it compiles under
// every mesh backend. See operations/stats.hpp for the contract.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

// Project includes
#include "meshioplusplus/operations/stats.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/cell_faces.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/polyhedron.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/operations/surface.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {

namespace {

using detail::Vec3;

// Sum the areas of the triangle/quad facets of a (linearized) surface mesh.
double stats_surface_area(const Mesh& rSurf) {
    const NDArray& points = rSurf.Points();
    const std::size_t pdim = rSurf.PointDim();
    double area = 0.0;
    for (const auto cb : rSurf.CellRange()) {
        if (cb.IsPolyhedron())
            continue;  // a volume cell, not a facet of one
        if (cb.IsRagged()) {
            for (std::size_t i = 0; i < cb.NumCells(); ++i) {
                const std::size_t n = cb.RowSize(i);
                std::vector<Vec3> coords(n);
                for (std::size_t k = 0; k < n; ++k)
                    coords[k] = detail::read_point(points, pdim, cb.Row(i)[k]);
                area += detail::polygon_area(coords.data(), n);
            }
            continue;
        }
        const CellType ct = cell_type_from_name(cb.Type());
        const int cc = detail::cell_corner_count(ct);
        if (cell_type_dimension(ct) != 2 || (cc != 3 && cc != 4))
            continue;
        const NDArray& conn = cb.Conn();
        const std::size_t npc = cb.NodesPerCell();
        const std::size_t nc = cb.NumCells();
        std::vector<double> vals(nc);
        parallel_for(nc, [&](std::size_t i) {
            std::vector<Vec3> coords;
            detail::read_corner_coords(points, pdim, conn, i * npc, static_cast<std::size_t>(cc),
                                       coords);
            vals[i] = detail::polygon_area(coords.data(), static_cast<std::size_t>(cc));
        });
        for (double v : vals)
            area += v;
    }
    return area;
}

}  // namespace

StatsReport compute_stats(const Mesh& rMesh) {
    StatsReport rep;
    const NDArray& points = rMesh.Points();
    const std::size_t n = rMesh.NumPoints();
    const std::size_t dim = rMesh.PointDim();
    rep.mNumPoints = static_cast<std::int64_t>(n);

    // --- bounding box / centroid (chunked parallel reduction) ---------------
    if (n > 0 && dim > 0) {
        const std::size_t ddim = std::min<std::size_t>(dim, 3);
        const std::size_t grain = 4096;
        const std::size_t nchunks = (n + grain - 1) / grain;
        std::vector<std::array<double, 3>> pmin(nchunks), pmax(nchunks), psum(nchunks);
        parallel_for(
            nchunks,
            [&](std::size_t ci) {
                std::array<double, 3> lmin = {std::numeric_limits<double>::infinity(),
                                              std::numeric_limits<double>::infinity(),
                                              std::numeric_limits<double>::infinity()};
                std::array<double, 3> lmax = {-std::numeric_limits<double>::infinity(),
                                              -std::numeric_limits<double>::infinity(),
                                              -std::numeric_limits<double>::infinity()};
                std::array<double, 3> lsum = {0, 0, 0};
                const std::size_t start = ci * grain;
                const std::size_t stop = std::min(n, start + grain);
                for (std::size_t g = start; g < stop; ++g)
                    for (std::size_t d = 0; d < ddim; ++d) {
                        const double v = detail::read_double(points, g * dim + d);
                        lmin[d] = std::min(lmin[d], v);
                        lmax[d] = std::max(lmax[d], v);
                        lsum[d] += v;
                    }
                pmin[ci] = lmin;
                pmax[ci] = lmax;
                psum[ci] = lsum;
            },
            1);
        std::array<double, 3> gmin = {std::numeric_limits<double>::infinity(),
                                      std::numeric_limits<double>::infinity(),
                                      std::numeric_limits<double>::infinity()};
        std::array<double, 3> gmax = {-std::numeric_limits<double>::infinity(),
                                      -std::numeric_limits<double>::infinity(),
                                      -std::numeric_limits<double>::infinity()};
        std::array<double, 3> gsum = {0, 0, 0};
        for (std::size_t ci = 0; ci < nchunks; ++ci)
            for (std::size_t d = 0; d < 3; ++d) {
                gmin[d] = std::min(gmin[d], pmin[ci][d]);
                gmax[d] = std::max(gmax[d], pmax[ci][d]);
                gsum[d] += psum[ci][d];
            }
        for (std::size_t d = 0; d < ddim; ++d) {
            rep.mBBoxMin[d] = gmin[d];
            rep.mBBoxMax[d] = gmax[d];
            rep.mExtent[d] = gmax[d] - gmin[d];
            rep.mCentroid[d] = gsum[d] / static_cast<double>(n);
        }
    }

    // --- per-cell-type counts + areas / volumes -----------------------------
    std::unordered_map<std::string, std::size_t> type_index;
    bool any_3d = false;
    for (const auto cb : rMesh.CellRange()) {
        const std::string type(cb.Type());
        const CellType ct = cell_type_from_name(type);
        const std::size_t nc = cb.NumCells();
        rep.mNumCells += static_cast<std::int64_t>(nc);

        auto it = type_index.find(type);
        if (it == type_index.end()) {
            type_index[type] = rep.mCellTypeCounts.size();
            rep.mCellTypeCounts.emplace_back(type, static_cast<std::int64_t>(nc));
        } else {
            rep.mCellTypeCounts[it->second].second += static_cast<std::int64_t>(nc);
        }

        if (cb.IsPolyhedron()) {
            // A polyhedron block's face winding is NOT a contract -- meshio++
            // repairs it rather than requiring it (doc/polyhedra.md) -- so
            // "inverted" is not a question one can ask of such a cell: there is
            // no convention for it to violate. Orient, then take the (now
            // positive) volume. What IS meaningful is an unorientable cell,
            // which is the polyhedral analogue of degenerate and is counted as
            // such rather than silently contributing zero.
            any_3d = true;
            std::vector<double> vals(nc);
            parallel_for(nc, [&](std::size_t i) {
                detail::CellRings rings;
                std::vector<Vec3> coords;
                if (!detail::cell_rings(cb, i, points, dim, rings, coords)) {
                    vals[i] = std::nan("");
                    return;
                }
                if (detail::orient_rings(rings, coords.data()) ==
                    detail::RingOrientation::Unorientable) {
                    vals[i] = std::nan("");
                    return;
                }
                vals[i] = std::abs(detail::poly_measure(rings, coords.data()).mVolume);
            });
            std::size_t n_bad = 0;
            for (double v : vals) {
                if (std::isnan(v)) {
                    ++n_bad;
                    continue;
                }
                rep.mSignedVolume += v;
                rep.mUnsignedVolume += v;
            }
            if (n_bad > 0)
                log::warn(
                    "stats: {} polyhedron cell(s) in block '{}' are not closed orientable "
                    "surfaces; their volume is not defined and is excluded",
                    n_bad, type);
            continue;
        }
        if (cb.IsRagged()) {
            // 1-level ragged (jagged polygons): a surface, so it contributes
            // area. Before v9.16.0 such a block contributed nothing at all.
            std::vector<double> vals(nc);
            parallel_for(nc, [&](std::size_t i) {
                const std::size_t n = cb.RowSize(i);
                std::vector<Vec3> coords(n);
                for (std::size_t k = 0; k < n; ++k)
                    coords[k] = detail::read_point(points, dim, cb.Row(i)[k]);
                vals[i] = detail::polygon_area(coords.data(), n);
            });
            for (double v : vals)
                rep.mTotalArea += v;
            continue;
        }
        const int cc = detail::cell_corner_count(ct);
        const int cdim = cell_type_dimension(ct);
        if (cc <= 0)
            continue;
        const NDArray& conn = cb.Conn();
        const std::size_t npc = cb.NodesPerCell();

        if (cdim == 2 && (cc == 3 || cc == 4)) {
            std::vector<double> vals(nc);
            parallel_for(nc, [&](std::size_t i) {
                std::vector<Vec3> coords;
                detail::read_corner_coords(points, dim, conn, i * npc, static_cast<std::size_t>(cc),
                                           coords);
                vals[i] = detail::polygon_area(coords.data(), static_cast<std::size_t>(cc));
            });
            for (double v : vals)
                rep.mTotalArea += v;
        } else if (cdim == 3) {
            // Deliberately NOT oriented: `cell_faces`' rows are canonically
            // outward-wound, so the sign here genuinely measures whether the
            // CELL is inverted. Repairing the winding first would make every
            // cell positive and silently zero `mNumInverted`.
            any_3d = true;
            std::vector<double> vals(nc);
            parallel_for(nc, [&](std::size_t i) {
                detail::CellRings rings;
                std::vector<Vec3> coords;
                vals[i] = detail::cell_rings(cb, i, points, dim, rings, coords)
                              ? detail::poly_measure(rings, coords.data()).mVolume
                              : std::nan("");
            });
            for (double v : vals) {
                if (std::isnan(v))
                    continue;
                rep.mSignedVolume += v;
                rep.mUnsignedVolume += std::abs(v);
                if (v < 0.0)
                    ++rep.mNumInverted;
            }
        }
    }

    // --- boundary area of 3D cells (via the surface operation) --------------
    if (any_3d) {
        try {
            rep.mTotalArea += stats_surface_area(extract_surface(rMesh));
        } catch (...) {
            // Boundary extraction failed (unsupported cell mix) — leave 3D
            // boundary area out rather than failing the whole stats report.
        }
    }

    return rep;
}

}  // namespace meshioplusplus

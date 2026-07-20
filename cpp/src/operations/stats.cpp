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
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/operations/surface.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {

namespace {

using detail::Vec3;

// Unsigned area of a corner polygon (triangle / quad) via the Newell normal.
double stats_area(const std::vector<Vec3>& rC, int corners) {
    if (corners < 3)
        return 0.0;
    Vec3 s = {0, 0, 0};
    for (int i = 0; i < corners; ++i)
        s = detail::vec3_add(s, detail::vec3_cross(rC[i], rC[(i + 1) % corners]));
    return 0.5 * detail::vec3_norm(s);
}

// Signed volume of a 3D cell via the divergence theorem over its outward-wound
// boundary faces (triangulated fan). NaN for types without a face table.
double stats_signed_volume(const std::vector<Vec3>& rC, CellType ct) {
    const std::vector<detail::CellFaceDef>& faces = detail::cell_faces(ct);
    if (faces.empty())
        return std::nan("");
    double vol6 = 0.0;
    for (const detail::CellFaceDef& f : faces) {
        const Vec3 a = rC[f.mNodes[0]];
        for (int i = 1; i + 1 < f.mNumCorners; ++i)
            vol6 += detail::triple_product(a, rC[f.mNodes[i]], rC[f.mNodes[i + 1]]);
    }
    return vol6 / 6.0;
}

// Sum the areas of the triangle/quad facets of a (linearized) surface mesh.
double stats_surface_area(const Mesh& rSurf) {
    const NDArray& points = rSurf.Points();
    const std::size_t pdim = rSurf.PointDim();
    double area = 0.0;
    for (const auto cb : rSurf.CellRange()) {
        if (cb.IsRagged() || cb.IsPolyhedron())
            continue;
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
            vals[i] = stats_area(coords, cc);
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

        if (cb.IsRagged() || cb.IsPolyhedron())
            continue;
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
                vals[i] = stats_area(coords, cc);
            });
            for (double v : vals)
                rep.mTotalArea += v;
        } else if (cdim == 3) {
            any_3d = true;
            std::vector<double> vals(nc);
            parallel_for(nc, [&](std::size_t i) {
                std::vector<Vec3> coords;
                detail::read_corner_coords(points, dim, conn, i * npc, static_cast<std::size_t>(cc),
                                           coords);
                vals[i] = stats_signed_volume(coords, ct);
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

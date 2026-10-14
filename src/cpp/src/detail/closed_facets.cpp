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
// The surface cells as the facets of a closed shell. See closed_facets.hpp.

// System includes
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

// Project includes (private, not installed)
#include "closed_facets.hpp"
#include "crease_edges.hpp"
#include "surface_rings.hpp"

namespace meshioplusplus {
namespace detail {

ClosedFacets closed_facets(const Mesh& rSurface, const std::string& rPrefix) {
    ClosedFacets out;
    std::size_t total = 0;
    for (const auto cb : rSurface.CellRange())
        total += cb.NumCells();
    const SurfaceRings in = surface_rings(rSurface, std::vector<std::uint8_t>(total, 1), rPrefix);
    out.mNumSkipped = total - in.NumFaces();
    if (in.NumFaces() == 0)
        return out;

    out.mXyz = surface_xyz(rSurface);
    const std::vector<double>& xyz = out.mXyz;
    for (double v : xyz)
        if (!std::isfinite(v))
            throw std::invalid_argument(rPrefix + "a coordinate is not finite");

    // Fan every non-planar facet into triangles from its first corner.
    std::vector<std::int64_t>& start = out.mStart;
    std::vector<std::int64_t>& nodes = out.mNodes;
    for (std::size_t f = 0; f < in.NumFaces(); ++f) {
        const std::int64_t* ring = in.mNodes.data() + in.mStart[f];
        const std::size_t size = static_cast<std::size_t>(in.mStart[f + 1] - in.mStart[f]);
        bool planar = true;
        if (size > 3) {
            double n[3];
            ring_unit_normal(xyz.data(), ring, size, n);
            const double* p0 = &xyz[static_cast<std::size_t>(ring[0]) * 3];
            double dev = 0.0, extent = 0.0;
            for (std::size_t k = 1; k < size; ++k) {
                const double* p = &xyz[static_cast<std::size_t>(ring[k]) * 3];
                const double d[3] = {p[0] - p0[0], p[1] - p0[1], p[2] - p0[2]};
                const double h = std::fabs(d[0] * n[0] + d[1] * n[1] + d[2] * n[2]);
                const double l = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
                dev = h > dev ? h : dev;
                extent = l > extent ? l : extent;
            }
            planar = !(dev > 1e-9 * extent);
        }
        if (planar) {
            nodes.insert(nodes.end(), ring, ring + size);
            start.push_back(static_cast<std::int64_t>(nodes.size()));
            out.mCell.push_back(in.mCell[f]);
        } else {
            for (std::size_t k = 1; k + 1 < size; ++k) {
                nodes.insert(nodes.end(), {ring[0], ring[k], ring[k + 1]});
                start.push_back(static_cast<std::int64_t>(nodes.size()));
                out.mCell.push_back(in.mCell[f]);
            }
        }
    }
    const std::size_t nf = out.NumFacets();
    std::vector<double>& normals = out.mNormals;
    normals.assign(nf * 3, 0.0);
    for (std::size_t f = 0; f < nf; ++f) {
        ring_unit_normal(xyz.data(), nodes.data() + start[f],
                         static_cast<std::size_t>(start[f + 1] - start[f]), &normals[f * 3]);
        if (normals[f * 3] == 0.0 && normals[f * 3 + 1] == 0.0 && normals[f * 3 + 2] == 0.0)
            throw std::invalid_argument(rPrefix + "facet " + std::to_string(f) +
                                        " is degenerate (zero area)");
    }

    // A solid's shell is closed and consistently wound.
    const std::vector<CreaseEdge> bad = crease_edges(start, nodes, normals, 180.0);
    std::size_t open = 0, nonmanifold = 0, inconsistent = 0;
    std::string first;
    std::size_t listed = 0;
    for (const CreaseEdge& e : bad) {
        const bool b = e.IsBoundary(), m = e.IsNonManifold(), i = e.mInconsistent;
        if (!b && !m && !i)
            continue;
        open += b ? 1 : 0;
        nonmanifold += m ? 1 : 0;
        inconsistent += i ? 1 : 0;
        if (listed < 5) {
            first +=
                (listed ? ", (" : "(") + std::to_string(e.mLo) + ", " + std::to_string(e.mHi) + ")";
            ++listed;
        }
    }
    if (open + nonmanifold + inconsistent > 0)
        throw std::invalid_argument(
            rPrefix + "the surface is not a closed, consistently wound solid: " +
            std::to_string(open) + " open, " + std::to_string(nonmanifold) + " non-manifold and " +
            std::to_string(inconsistent) + " inconsistently wound edge(s), first " + first +
            " (point ids); run repair first");
    return out;
}

}  // namespace detail
}  // namespace meshioplusplus

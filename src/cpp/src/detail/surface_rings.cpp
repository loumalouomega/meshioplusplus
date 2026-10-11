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

// The surface cells as face rings. See surface_rings.hpp.

// System includes
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/cell_edges.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/parallel.hpp"

// Project includes (private, not installed)
#include "surface_rings.hpp"

namespace meshioplusplus {
namespace detail {

SurfaceRings surface_rings(const Mesh& rMesh, const std::vector<std::uint8_t>& rSelected,
                           const std::string& rErrorPrefix) {
    SurfaceRings out;
    const std::size_t n = rMesh.NumPoints();
    const std::vector<std::int64_t> bases = block_bases(rMesh);
    std::size_t bi = 0;
    std::vector<std::int64_t> ids;
    for (const auto cb : rMesh.CellRange()) {
        const std::size_t base = static_cast<std::size_t>(bases[bi++]);
        const std::string type(cb.Type());
        const CellType ct = cell_type_from_name(type);
        if (cb.IsPolyhedron() || cell_type_dimension(ct) == 3)
            throw std::invalid_argument(rErrorPrefix + "block '" + type +
                                        "' is a volume block; segment a surface (take the "
                                        "skin first with extract_skin)");
        const bool polygon = type.rfind("polygon", 0) == 0;
        std::vector<std::size_t> corners;
        if (!polygon) {
            if (cell_type_dimension(ct) != 2 || cb.IsRagged())
                continue;
            for (const CellEdgeDef& ed : cell_edges(ct))
                corners.push_back(ed.mNodes[0]);
            if (corners.size() < 3)
                continue;
        }
        const std::size_t npc = cb.IsRagged() ? 0 : cb.NodesPerCell();
        for (std::size_t c = 0; c < cb.NumCells(); ++c) {
            if (!rSelected[base + c])
                continue;
            ids.clear();
            if (polygon && cb.IsRagged()) {
                ids.assign(cb.Row(c), cb.Row(c) + cb.RowSize(c));
            } else if (polygon) {
                // A polygon block can be stored rectangular: one row per cell.
                for (std::size_t k = 0; k < npc; ++k)
                    ids.push_back(read_int(cb.Conn(), c * npc + k));
            } else {
                for (std::size_t k : corners)
                    ids.push_back(read_int(cb.Conn(), c * npc + k));
            }
            if (ids.size() < 3)
                continue;
            bool ok = true;
            for (std::int64_t v : ids)
                ok = ok && v >= 0 && static_cast<std::size_t>(v) < n;
            if (!ok)
                continue;
            out.mNodes.insert(out.mNodes.end(), ids.begin(), ids.end());
            out.mStart.push_back(static_cast<std::int64_t>(out.mNodes.size()));
            out.mCell.push_back(static_cast<std::int64_t>(base + c));
        }
    }
    return out;
}

std::vector<double> surface_xyz(const Mesh& rMesh) {
    const std::size_t n = rMesh.NumPoints();
    const std::size_t dim = rMesh.PointDim();
    std::vector<double> xyz(n * 3, 0.0);
    const NDArray& pts = rMesh.Points();
    parallel_for_bw(n, [&](std::size_t i) {
        for (std::size_t d = 0; d < dim && d < 3; ++d)
            xyz[i * 3 + d] = read_double(pts, i * dim + d);
    });
    return xyz;
}

}  // namespace detail
}  // namespace meshioplusplus

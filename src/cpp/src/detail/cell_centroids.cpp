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

// Per-cell centroids of a whole mesh (detail/cell_centroids.hpp).

// System includes
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

// Project includes
#include "cell_centroids.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/parallel.hpp"
#include "typed_view.hpp"

namespace meshioplusplus::detail {

std::vector<double> cell_centroids(const Mesh& rMesh, std::size_t TotalCells) {
    const NDArray& points = rMesh.Points();
    const DoubleView points_v(points);
    const std::size_t pdim = rMesh.PointDim();
    const std::size_t dim = std::min<std::size_t>(pdim, 3);

    std::vector<double> cent(3 * TotalCells, 0.0);
    std::size_t base = 0;
    for (const auto cb : rMesh.CellRange()) {
        const std::size_t ncells = cb.NumCells();
        if (cb.IsPolyhedron()) {
            parallel_for(ncells, [&, base](std::size_t c) {
                static thread_local std::vector<std::int64_t> nodes;
                nodes.clear();
                for (std::size_t f = 0; f < cb.NumFaces(c); ++f) {
                    const auto face = cb.Face(c, f);
                    nodes.insert(nodes.end(), face.first, face.first + face.second);
                }
                std::sort(nodes.begin(), nodes.end());
                nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());
                double* out = cent.data() + 3 * (base + c);
                for (std::int64_t node : nodes)
                    for (std::size_t d = 0; d < dim; ++d)
                        out[d] += points_v[static_cast<std::size_t>(node) * pdim + d];
                if (!nodes.empty())
                    for (std::size_t d = 0; d < dim; ++d)
                        out[d] /= static_cast<double>(nodes.size());
            });
        } else if (cb.IsRagged()) {
            parallel_for(ncells, [&, base](std::size_t c) {
                const std::int64_t* row = cb.Row(c);
                const std::size_t nn = cb.RowSize(c);
                double* out = cent.data() + 3 * (base + c);
                for (std::size_t k = 0; k < nn; ++k)
                    for (std::size_t d = 0; d < dim; ++d)
                        out[d] += points_v[static_cast<std::size_t>(row[k]) * pdim + d];
                if (nn > 0)
                    for (std::size_t d = 0; d < dim; ++d)
                        out[d] /= static_cast<double>(nn);
            });
        } else {
            const NDArray& conn = cb.Conn();
            const Int64View conn_v(conn);
            const std::size_t npc = cb.NodesPerCell();
            parallel_for(ncells, [&, base](std::size_t c) {
                double* out = cent.data() + 3 * (base + c);
                for (std::size_t k = 0; k < npc; ++k) {
                    const std::size_t node = static_cast<std::size_t>(conn_v[c * npc + k]);
                    for (std::size_t d = 0; d < dim; ++d)
                        out[d] += points_v[node * pdim + d];
                }
                if (npc > 0)
                    for (std::size_t d = 0; d < dim; ++d)
                        out[d] /= static_cast<double>(npc);
            });
        }
        base += ncells;
    }
    return cent;
}

}  // namespace meshioplusplus::detail

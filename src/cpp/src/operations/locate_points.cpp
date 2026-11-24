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
// Point location: see operations/locate_points.hpp for the contract. The
// simplexified mesh and the search are `detail::PointLocator`, the one
// `interpolate`'s barycentric mode uses.

// System includes
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/operations/locate_points.hpp"
#include "meshioplusplus/operations/convert_cells.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/parallel.hpp"

// Project includes (private, not installed)
#include "../detail/point_locator.hpp"

namespace meshioplusplus {

LocatePointsResult locate_points(const Mesh& rMesh, const std::vector<double>& rPoints,
                                 const LocatePointsOptions& rOptions) {
    static const char* const kPrefix = "meshio++: locate_points: ";
    if (rPoints.size() % 3 != 0)
        throw std::invalid_argument(std::string(kPrefix) +
                                    "points must hold three numbers per point");
    if (!std::isfinite(rOptions.mTolerance) || rOptions.mTolerance < 0.0)
        throw std::invalid_argument(std::string(kPrefix) +
                                    "tolerance must be a finite number, zero or above");
    for (std::size_t k = 0; k < rPoints.size(); ++k)
        if (!std::isfinite(rPoints[k]))
            throw std::invalid_argument(std::string(kPrefix) + "point " + std::to_string(k / 3) +
                                        " is not finite");
    const std::size_t n = rPoints.size() / 3;

    // Reduce to simplices, remembering each one's input cell.
    ConvertCellsOptions cco;
    cco.mMode = ConvertCellsMode::Simplexify;
    cco.mRecordParentIds = true;
    const ConvertCellsResult simp = convert_cells(rMesh, cco);
    const Mesh& smesh = simp.mMesh;
    const detail::PointLocator locator(smesh, rOptions.mTolerance, "locate_points");
    const std::size_t nv = locator.NodesPerSimplex();
    const std::vector<std::int64_t>& sconn = locator.Connectivity();

    // Simplex index -> global input cell: the locator takes the blocks of one
    // type in ascending order, and block b of the result is input block b, so
    // the parent id (local to its block) plus the input block's base is it.
    const std::string want = nv == 4 ? "tetra" : "triangle";
    const std::vector<std::int64_t> bases = detail::block_bases(rMesh);
    std::vector<std::int64_t> source_cell;
    source_cell.reserve(locator.NumSimplices());
    std::size_t block = 0;
    for (const auto cb : smesh.CellRange()) {
        if (cb.Type() == want && cb.NumCells() > 0) {
            const NDArray& parents = smesh.CellData("convert:parent_cell", block);
            for (std::size_t c = 0; c < cb.NumCells(); ++c)
                source_cell.push_back(bases[block] + detail::read_int(parents, c));
        }
        ++block;
    }

    LocatePointsResult out;
    out.mCell.assign(n, -1);
    out.mNodes.assign(4 * n, -1);
    out.mBarycentric.assign(4 * n, 0.0);
    parallel_for(n, [&](std::size_t i) {
        const detail::Vec3 p{rPoints[3 * i], rPoints[3 * i + 1], rPoints[3 * i + 2]};
        double w[4] = {0.0, 0.0, 0.0, 0.0};
        const std::int64_t s = locator.Locate(p, w);
        if (s < 0)
            return;
        out.mCell[i] = source_cell[static_cast<std::size_t>(s)];
        for (std::size_t k = 0; k < nv; ++k) {
            out.mNodes[4 * i + k] = sconn[static_cast<std::size_t>(s) * nv + k];
            out.mBarycentric[4 * i + k] = w[k];
        }
    });
    for (std::size_t i = 0; i < n; ++i)
        out.mNumLocated += out.mCell[i] >= 0 ? 1 : 0;
    return out;
}

}  // namespace meshioplusplus

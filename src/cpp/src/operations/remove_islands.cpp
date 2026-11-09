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
// Drop a scan's crumbs (roadmap 7.8.1). See operations/remove_islands.hpp.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/operations/remove_islands.hpp"
#include "../detail/cell_components.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/subset.hpp"
#include "meshioplusplus/detail/value_io.hpp"

namespace meshioplusplus {

namespace {

constexpr const char* kRiPrefix = "meshio++: remove_islands: ";

NDArray ri_int64(const std::vector<std::int64_t>& rV) {
    NDArray a = NDArray::Uninit(DType::Int64, {rV.size()});
    std::copy(rV.begin(), rV.end(), a.As<std::int64_t>());
    return a;
}

}  // namespace

IslandConnectivity island_connectivity_from_name(const std::string& rName) {
    if (rName == "edge")
        return IslandConnectivity::Edge;
    if (rName == "point")
        return IslandConnectivity::Point;
    throw std::invalid_argument(std::string(kRiPrefix) +
                                "connectivity must be 'edge' or 'point', not '" + rName + "'");
}

RemoveIslandsResult remove_islands(const Mesh& rMesh, const RemoveIslandsOptions& rOptions) {
    for (const double v : {rOptions.mMinDiagonal, rOptions.mMinArea})
        if (!std::isfinite(v) || v < 0.0)
            throw std::invalid_argument(std::string(kRiPrefix) +
                                        "min_diagonal and min_area must be finite and >= 0");
    if (rOptions.mMinCells < 0)
        throw std::invalid_argument(std::string(kRiPrefix) + "min_cells must be >= 0");
    if (rOptions.mConnectivity != IslandConnectivity::Edge &&
        rOptions.mConnectivity != IslandConnectivity::Point)
        throw std::invalid_argument(std::string(kRiPrefix) + "unknown connectivity");

    const detail::CellComponents comps = rOptions.mConnectivity == IslandConnectivity::Edge
                                             ? detail::cell_components_by_edge(rMesh)
                                             : detail::cell_components(rMesh);
    const auto ncomp = static_cast<std::size_t>(comps.mCount);
    RemoveIslandsResult result;
    result.mNumComponents = comps.mCount;
    result.mComponentCells.assign(ncomp, 0);
    result.mComponentArea.assign(ncomp, 0.0);
    constexpr double kInf = std::numeric_limits<double>::infinity();
    const std::array<double, 6> empty{kInf, kInf, kInf, -kInf, -kInf, -kInf};
    std::vector<std::array<double, 6>> box(ncomp, empty);
    std::array<double, 6> whole = empty;
    const NDArray& pts = rMesh.Points();
    const std::size_t dim = rMesh.PointDim();
    auto point = [&](std::int64_t V) {
        std::array<double, 3> p{};
        for (std::size_t k = 0; k < dim && k < 3; ++k)
            p[k] = detail::read_double(pts, static_cast<std::size_t>(V) * dim + k);
        return p;
    };
    double total_area = 0.0;
    std::size_t g = 0;
    for (const auto cb : rMesh.CellRange()) {
        const std::string type(cb.Type());
        const bool face =
            !cb.IsPolyhedron() &&
            (type.rfind("polygon", 0) == 0 || cell_type_dimension(cell_type_from_name(type)) == 2);
        for (std::size_t c = 0; c < cb.NumCells(); ++c, ++g) {
            const auto k = static_cast<std::size_t>(comps.mLabel[g]);
            ++result.mComponentCells[k];
            const std::size_t n = cb.IsRagged() ? cb.RowSize(c) : cb.NodesPerCell();
            std::vector<std::array<double, 3>> corner(n);
            for (std::size_t i = 0; i < n; ++i) {
                const std::int64_t v = cb.IsRagged()
                                           ? cb.Row(c)[i]
                                           : detail::read_int(cb.Conn(), c * cb.NodesPerCell() + i);
                corner[i] = point(v);
                for (std::size_t q = 0; q < 3; ++q) {
                    box[k][q] = std::min(box[k][q], corner[i][q]);
                    box[k][3 + q] = std::max(box[k][3 + q], corner[i][q]);
                    whole[q] = std::min(whole[q], corner[i][q]);
                    whole[3 + q] = std::max(whole[3 + q], corner[i][q]);
                }
            }
            if (!face || n < 3)
                continue;
            // The fan's vector area (the polygon's, planar or not).
            std::array<double, 3> s{0.0, 0.0, 0.0};
            for (std::size_t i = 1; i + 1 < n; ++i) {
                const std::array<double, 3> a{corner[i][0] - corner[0][0],
                                              corner[i][1] - corner[0][1],
                                              corner[i][2] - corner[0][2]};
                const std::array<double, 3> b{corner[i + 1][0] - corner[0][0],
                                              corner[i + 1][1] - corner[0][1],
                                              corner[i + 1][2] - corner[0][2]};
                s[0] += a[1] * b[2] - a[2] * b[1];
                s[1] += a[2] * b[0] - a[0] * b[2];
                s[2] += a[0] * b[1] - a[1] * b[0];
            }
            const double area = 0.5 * std::sqrt(s[0] * s[0] + s[1] * s[1] + s[2] * s[2]);
            result.mComponentArea[k] += area;
            total_area += area;
        }
    }
    const double whole_diag = g == 0 ? 0.0
                                     : std::sqrt((whole[3] - whole[0]) * (whole[3] - whole[0]) +
                                                 (whole[4] - whole[1]) * (whole[4] - whole[1]) +
                                                 (whole[5] - whole[2]) * (whole[5] - whole[2]));
    const double min_diag =
        rOptions.mRelative ? rOptions.mMinDiagonal * whole_diag : rOptions.mMinDiagonal;
    const double min_area = rOptions.mRelative ? rOptions.mMinArea * total_area : rOptions.mMinArea;
    result.mComponentDiagonal.assign(ncomp, 0.0);
    result.mComponentRemoved.assign(ncomp, 0);
    for (std::size_t k = 0; k < ncomp; ++k) {
        const auto& b = box[k];
        result.mComponentDiagonal[k] =
            std::sqrt((b[3] - b[0]) * (b[3] - b[0]) + (b[4] - b[1]) * (b[4] - b[1]) +
                      (b[5] - b[2]) * (b[5] - b[2]));
        const bool small = result.mComponentDiagonal[k] < min_diag ||
                           result.mComponentArea[k] < min_area ||
                           result.mComponentCells[k] < rOptions.mMinCells;
        if (small) {
            result.mComponentRemoved[k] = 1;
            ++result.mNumRemoved;
            result.mNumCellsRemoved += result.mComponentCells[k];
        }
    }

    if (rOptions.mDryRun) {
        result.mMesh = detail::clone_mesh(
            rMesh, [](DataLocation, const std::string&, std::string&) { return true; });
        std::vector<NDArray> comp_blocks, gone_blocks;
        std::size_t base = 0;
        for (const auto cb : rMesh.CellRange()) {
            std::vector<std::int64_t> lab(cb.NumCells()), gone(cb.NumCells());
            for (std::size_t c = 0; c < cb.NumCells(); ++c) {
                lab[c] = comps.mLabel[base + c];
                gone[c] = result.mComponentRemoved[static_cast<std::size_t>(lab[c])];
            }
            comp_blocks.push_back(ri_int64(lab));
            gone_blocks.push_back(ri_int64(gone));
            base += cb.NumCells();
        }
        result.mMesh.AddCellData("islands:component", std::move(comp_blocks));
        result.mMesh.AddCellData("islands:removed", std::move(gone_blocks));
        return result;
    }
    std::vector<std::vector<std::int64_t>> kept(rMesh.NumCellBlocks());
    std::size_t base = 0;
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const std::size_t nc = rMesh.Cells(b).NumCells();
        for (std::size_t c = 0; c < nc; ++c)
            if (!result.mComponentRemoved[static_cast<std::size_t>(comps.mLabel[base + c])])
                kept[b].push_back(static_cast<std::int64_t>(c));
        base += nc;
    }
    detail::SubsetResult sub =
        detail::build_cell_subset(rMesh, kept, "", "", true, "remove_islands");
    result.mNumPointsRemoved = static_cast<std::int64_t>(rMesh.NumPoints()) -
                               static_cast<std::int64_t>(sub.mMesh.NumPoints());
    result.mMesh = std::move(sub.mMesh);
    return result;
}

}  // namespace meshioplusplus

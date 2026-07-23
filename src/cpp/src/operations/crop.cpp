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
// Spatial subsetting by bounding box or half-space. A per-point inside mask is
// built in parallel, cells are kept by an all/any policy, then the shared
// detail::build_cell_subset prunes points and remaps connectivity + data. Built
// entirely through the uniform mesh API so it compiles under every mesh backend.
// See operations/crop.hpp for the contract.

// System includes
#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/operations/crop.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/subset.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {

namespace {

using detail::Vec3;

// Build the per-point inside mask for a predicate `inside(point)`.
template <class Predicate>
std::vector<char> crop_point_mask(const Mesh& rMesh, Predicate inside) {
    const std::size_t n = rMesh.NumPoints();
    const std::size_t dim = rMesh.PointDim();
    const NDArray& points = rMesh.Points();
    std::vector<char> mask(n, 0);
    parallel_for(n, [&](std::size_t g) {
        mask[g] = inside(detail::read_point(points, dim, static_cast<std::int64_t>(g))) ? 1 : 0;
    });
    return mask;
}

// Collect the kept cells per block under the all/any policy.
std::vector<std::vector<std::int64_t>> crop_kept_cells(const Mesh& rMesh,
                                                       const std::vector<char>& rMask,
                                                       CropMode mode) {
    std::vector<std::vector<std::int64_t>> kept;
    kept.reserve(rMesh.NumCellBlocks());
    const bool all = mode == CropMode::All;
    for (const auto cb : rMesh.CellRange()) {
        const std::size_t nc = cb.NumCells();
        // Phase 1: per-cell keep flag in parallel (disjoint slots); phase 2:
        // serial compaction, preserving the ascending kept order.
        std::vector<char> keep(nc, 0);
        parallel_for(nc, [&](std::size_t c) {
            bool any_in = false, all_in = true;
            auto visit = [&](std::int64_t v) {
                const bool in = rMask[static_cast<std::size_t>(v)] != 0;
                any_in = any_in || in;
                all_in = all_in && in;
            };
            if (cb.IsPolyhedron()) {
                for (std::size_t f = 0; f < cb.NumFaces(c); ++f) {
                    auto face = cb.Face(c, f);
                    for (std::size_t k = 0; k < face.second; ++k)
                        visit(face.first[k]);
                }
            } else if (cb.IsRagged()) {
                for (std::size_t k = 0; k < cb.RowSize(c); ++k)
                    visit(cb.Row(c)[k]);
            } else {
                const NDArray& conn = cb.Conn();
                const std::size_t npc = cb.NodesPerCell();
                for (std::size_t k = 0; k < npc; ++k)
                    visit(detail::read_int(conn, c * npc + k));
            }
            keep[c] = (all ? all_in : any_in) ? 1 : 0;
        });
        std::vector<std::int64_t> block;
        for (std::size_t c = 0; c < nc; ++c)
            if (keep[c])
                block.push_back(static_cast<std::int64_t>(c));
        kept.push_back(std::move(block));
    }
    return kept;
}

CropResult crop_finish(const Mesh& rMesh, const std::vector<std::vector<std::int64_t>>& rKept,
                       bool record_ids) {
    detail::SubsetResult sub =
        detail::build_cell_subset(rMesh, rKept, record_ids ? "crop:original_point_id" : "",
                                  record_ids ? "crop:original_cell_id" : "",
                                  /*drop_empty_blocks=*/false, "crop");
    CropResult res;
    res.mMesh = std::move(sub.mMesh);
    res.mPointMap = std::move(sub.mPointMap);
    res.mCellMaps = std::move(sub.mCellMaps);
    return res;
}

}  // namespace

CropResult crop_bbox(const Mesh& rMesh, const double* pLo, const double* pHi, CropMode mode,
                     bool record_ids) {
    const std::array<double, 3> lo = {pLo[0], pLo[1], pLo[2]};
    const std::array<double, 3> hi = {pHi[0], pHi[1], pHi[2]};
    auto mask = crop_point_mask(rMesh, [&](const Vec3& p) {
        for (int d = 0; d < 3; ++d)
            if (p[d] < lo[d] || p[d] > hi[d])
                return false;
        return true;
    });
    return crop_finish(rMesh, crop_kept_cells(rMesh, mask, mode), record_ids);
}

CropResult crop_halfspace(const Mesh& rMesh, const double* pPoint, const double* pNormal,
                          CropMode mode, bool record_ids) {
    const Vec3 p0 = {pPoint[0], pPoint[1], pPoint[2]};
    const Vec3 nrm = {pNormal[0], pNormal[1], pNormal[2]};
    auto mask = crop_point_mask(rMesh, [&](const Vec3& p) {
        return detail::vec3_dot(detail::vec3_sub(p, p0), nrm) >= 0.0;
    });
    return crop_finish(rMesh, crop_kept_cells(rMesh, mask, mode), record_ids);
}

}  // namespace meshioplusplus

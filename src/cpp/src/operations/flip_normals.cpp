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
// Reverse the winding of cells (roadmap 7.8.2). See operations/flip_normals.hpp.

// System includes
#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/flip_normals.hpp"
#include "../detail/cell_flip.hpp"
#include "../detail/cell_region_mask.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/region.hpp"

namespace meshioplusplus {

FlipNormalsResult flip_normals(const Mesh& rMesh, const std::string& rRegion) {
    std::vector<std::int64_t> bases{0};
    for (const auto cb : rMesh.CellRange())
        bases.push_back(bases.back() + static_cast<std::int64_t>(cb.NumCells()));
    const auto ncells = static_cast<std::size_t>(bases.back());
    const std::vector<std::uint8_t> mask =
        detail::cell_region_mask(rMesh, rRegion, ncells, "meshio++: flip_normals: ");
    std::vector<detail::BlockFlip> flips;
    for (const auto cb : rMesh.CellRange())
        flips.push_back(detail::block_flip(cb, "flip_normals"));

    FlipNormalsResult result;
    Mesh& out = result.mMesh;
    out.AssignPoints(rMesh.Points());
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto cb = rMesh.Cells(b);
        const detail::BlockFlip& r_flip = flips[b];
        const std::string type(cb.Type());
        const std::size_t nc = cb.NumCells();
        const auto flipped = [&](std::size_t c) {
            return mask[static_cast<std::size_t>(bases[b]) + c] != 0;
        };
        if (r_flip.mPolyhedron) {
            std::vector<std::int64_t> flat, rows{0}, faces{0};
            for (std::size_t c = 0; c < nc; ++c) {
                const bool flip = flipped(c);
                result.mNumFlipped += flip ? 1 : 0;
                for (std::size_t f = 0; f < cb.NumFaces(c); ++f) {
                    const auto face = cb.Face(c, f);
                    for (std::size_t k = 0; k < face.second; ++k)
                        flat.push_back(face.first[flip && k > 0 ? face.second - k : k]);
                    rows.push_back(static_cast<std::int64_t>(flat.size()));
                }
                faces.push_back(static_cast<std::int64_t>(rows.size() - 1));
            }
            out.AddPolyhedronBlock(type, std::move(flat), std::move(rows), std::move(faces));
        } else if (r_flip.mRaggedPolygon) {
            std::vector<std::int64_t> flat, rows{0};
            for (std::size_t c = 0; c < nc; ++c) {
                const bool flip = flipped(c);
                result.mNumFlipped += flip ? 1 : 0;
                const std::int64_t* row = cb.Row(c);
                const std::size_t k = cb.RowSize(c);
                for (std::size_t i = 0; i < k; ++i)
                    flat.push_back(row[flip && i > 0 ? k - i : i]);
                rows.push_back(static_cast<std::int64_t>(flat.size()));
            }
            out.AddPolygonBlock(type, std::move(flat), std::move(rows));
        } else {
            const std::vector<int>& perm = r_flip.Perm();
            const std::size_t npc = cb.NodesPerCell();
            const NDArray& conn = cb.Conn();
            NDArray c_out = NDArray::Uninit(DType::Int64, {nc, npc});
            std::int64_t* p_out = c_out.As<std::int64_t>();
            for (std::size_t c = 0; c < nc; ++c) {
                const bool flip = flipped(c);
                result.mNumFlipped += flip ? 1 : 0;
                for (std::size_t k = 0; k < npc; ++k)
                    p_out[c * npc + k] = detail::read_int(
                        conn, c * npc + (flip ? static_cast<std::size_t>(perm[k]) : k));
            }
            out.AddCellBlock(type, std::move(c_out));
        }
    }
    for (const std::string& r_name : rMesh.PointDataNames())
        out.AddPointData(r_name, rMesh.PointData(r_name));
    for (const std::string& r_name : rMesh.CellDataNames()) {
        std::vector<NDArray> blocks;
        for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b)
            blocks.push_back(rMesh.CellData(r_name, b));
        out.AddCellData(r_name, std::move(blocks));
    }
    for (const std::string& r_name : rMesh.FieldDataNames())
        out.AddFieldData(r_name, rMesh.FieldData(r_name));
    // Regions: a flipped cell's Side entries move to the facet the flip gives.
    for (std::size_t k = 0; k < rMesh.NumRegions(); ++k) {
        Region r = rMesh.Region(k);
        if (r.mKind == RegionKind::Side && r.mEntries.Size() > 0) {
            const std::size_t ne = r.mEntries.Size() / 2;
            NDArray e = NDArray::Uninit(DType::Int64, {ne, std::size_t{2}});
            std::int64_t* p_e = e.As<std::int64_t>();
            for (std::size_t i = 0; i < ne; ++i) {
                const std::int64_t g = detail::read_int(r.mEntries, 2 * i);
                std::int64_t f = detail::read_int(r.mEntries, 2 * i + 1);
                const auto [b, c] = detail::global_to_block_row(bases, g);
                if (b < rMesh.NumCellBlocks() && g >= 0 && static_cast<std::size_t>(g) < ncells &&
                    mask[static_cast<std::size_t>(g)]) {
                    const auto cb = rMesh.Cells(b);
                    const auto uc = static_cast<std::size_t>(c);
                    const std::int64_t nfac =
                        flips[b].mPolyhedron      ? static_cast<std::int64_t>(cb.NumFaces(uc))
                        : flips[b].mRaggedPolygon ? static_cast<std::int64_t>(cb.RowSize(uc))
                                                  : static_cast<std::int64_t>(cb.NodesPerCell());
                    const std::int64_t nf = detail::block_flip_facet(flips[b], f, nfac);
                    if (nf >= 0)
                        f = nf;
                }
                p_e[2 * i] = g;
                p_e[2 * i + 1] = f;
            }
            r.mEntries = std::move(e);
        }
        out.AddRegion(std::move(r));
    }
    return result;
}

}  // namespace meshioplusplus

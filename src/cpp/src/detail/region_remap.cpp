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
// The shared region remapper. See detail/region_remap.hpp for the contract and
// for why the cell-map shape has to be stated explicitly rather than inferred.

// System includes
#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/detail/region_remap.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/cell_edges.hpp"
#include "meshioplusplus/detail/cell_faces.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/detail/provenance.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

/// Read entry `i` of a possibly-absent map array, or the identity.
std::int64_t rremap_lookup(const NDArray* pMap, std::int64_t i) {
    if (!pMap)
        return i;  // null map == identity
    if (i < 0 || static_cast<std::size_t>(i) >= pMap->Size())
        return -1;
    return read_int(*pMap, static_cast<std::size_t>(i));
}

/// Output block for input block @p b under @p rMaps (identity when unset).
std::size_t rremap_out_block(const RegionRemap& rMaps, std::size_t b) {
    if (rMaps.mBlockMap.empty())
        return b;
    return b < rMaps.mBlockMap.size() ? rMaps.mBlockMap[b] : kBlockDropped;
}

/// Build the Int64 entry array a remapped region owns.
NDArray rremap_entries(const std::vector<std::int64_t>& rFlat, std::size_t stride) {
    std::vector<std::size_t> shape;
    if (stride == 1)
        shape = {rFlat.size()};
    else
        shape = {rFlat.size() / stride, stride};
    NDArray out = NDArray::Uninit(DType::Int64, std::move(shape));
    for (std::size_t i = 0; i < rFlat.size(); ++i)
        out.As<std::int64_t>()[i] = rFlat[i];
    return out;
}

/**
 * @brief Every output global cell an input global cell maps to.
 *
 * The `FirstChild` scan walks forward to the next non-negative map entry, which
 * is also correct when intermediate parents were dropped (a dropped parent has
 * no children of its own, so the run still contains exactly this parent's).
 */
void rremap_children(const RegionRemap& rMaps, const std::vector<std::int64_t>& rInBases,
                     const std::vector<std::int64_t>& rOutBases, std::int64_t inGlobal,
                     std::vector<std::int64_t>& rChildren) {
    rChildren.clear();

    if (rMaps.mCellMapKind == CellMapKind::Global) {
        const std::int64_t g = rremap_lookup(rMaps.pGlobalCellMap, inGlobal);
        if (g >= 0 && g < total_cells(rOutBases))
            rChildren.push_back(g);
        return;
    }

    const auto [b, row] = global_to_block_row(rInBases, inGlobal);
    if (b == static_cast<std::size_t>(-1))
        return;
    const std::size_t ob = rremap_out_block(rMaps, b);
    if (ob == kBlockDropped || ob + 1 >= rOutBases.size())
        return;

    const NDArray* p_map =
        (rMaps.pCellMaps && b < rMaps.pCellMaps->size()) ? &(*rMaps.pCellMaps)[b] : nullptr;
    const std::int64_t first = rremap_lookup(p_map, row);
    if (first < 0)
        return;

    const std::int64_t out_block_cells = rOutBases[ob + 1] - rOutBases[ob];
    if (rMaps.mCellMapKind == CellMapKind::Direct) {
        if (first < out_block_cells)
            rChildren.push_back(rOutBases[ob] + first);
        return;
    }

    // FirstChild: the run ends at the next non-negative map entry.
    std::int64_t end = out_block_cells;
    if (p_map) {
        const std::int64_t n_in = static_cast<std::int64_t>(p_map->Size());
        for (std::int64_t c = row + 1; c < n_in; ++c) {
            const std::int64_t nxt = read_int(*p_map, static_cast<std::size_t>(c));
            if (nxt >= 0) {
                end = nxt;
                break;
            }
        }
    } else {
        end = first + 1;
    }
    for (std::int64_t k = first; k < end && k < out_block_cells; ++k)
        rChildren.push_back(rOutBases[ob] + k);
}

/// The meshio type name of an output global cell, or an empty string.
std::string rremap_type_at(const Mesh& rMesh, const std::vector<std::int64_t>& rBases,
                           std::int64_t global) {
    const auto [b, row] = global_to_block_row(rBases, global);
    (void)row;
    if (b == static_cast<std::size_t>(-1))
        return {};
    return std::string(rMesh.Cells(b).Type());
}

}  // namespace

std::size_t region_num_facets(const std::string& rType) {
    const CellType type = cell_type_from_name(rType);
    const int dim = cell_type_dimension(type);
    if (dim == 3)
        return cell_faces(type).size();
    if (dim == 2)
        return cell_edges(type).size();
    return 0;
}

bool remap_region(const Mesh& rIn, const Mesh& rOut, const Region& rRegion,
                  const RegionRemap& rMaps, Region& rResult) {
    const std::vector<std::int64_t> in_bases = block_bases(rIn);
    const std::vector<std::int64_t> out_bases = block_bases(rOut);
    const auto n_out_points = static_cast<std::int64_t>(rOut.NumPoints());

    const std::int64_t* entries = rRegion.Entries();
    const std::size_t n = rRegion.NumEntries();
    std::vector<std::int64_t> children;
    std::vector<std::int64_t> out;

    if (rRegion.mKind == RegionKind::Point) {
        out.reserve(n);
        for (std::size_t k = 0; k < n; ++k) {
            const std::int64_t p = rremap_lookup(rMaps.pPointMap, entries[k]);
            if (p >= 0 && p < n_out_points)
                out.push_back(p);
        }
    } else if (rRegion.mKind == RegionKind::Cell) {
        out.reserve(n);
        for (std::size_t k = 0; k < n; ++k) {
            rremap_children(rMaps, in_bases, out_bases, entries[k], children);
            out.insert(out.end(), children.begin(), children.end());
        }
    } else {
        // Side: the facet must still exist, which needs a 1:1 cell map onto a
        // cell of the same type. Anything else has no facet correspondence.
        if (rMaps.mDropSideRegions || rMaps.mCellMapKind == CellMapKind::FirstChild) {
            log::warn(
                "{}: side region '{}' dropped — this operation does not preserve "
                "facet identity",
                rMaps.mOpName.empty() ? "operation" : rMaps.mOpName, rRegion.mName);
            return false;
        }
        out.reserve(n * 2);
        for (std::size_t k = 0; k < n; ++k) {
            const std::int64_t cell = entries[k * 2];
            const std::int64_t facet = entries[k * 2 + 1];
            rremap_children(rMaps, in_bases, out_bases, cell, children);
            if (children.size() != 1)
                continue;
            const std::string in_type = rremap_type_at(rIn, in_bases, cell);
            const std::string out_type = rremap_type_at(rOut, out_bases, children[0]);
            if (in_type.empty() || in_type != out_type)
                continue;
            if (facet < 0 || static_cast<std::size_t>(facet) >= region_num_facets(out_type))
                continue;
            out.push_back(children[0]);
            out.push_back(facet);
        }
    }

    // A region that lost every entry is still carried, as an empty group. The
    // name is information in its own right — a crop or a partition piece that
    // happens to contain none of "inlet"'s cells has still been told that an
    // "inlet" exists — and it is what the per-operation Python shims did before
    // this helper replaced them.
    rResult = Region(rRegion.mName, rRegion.mKind, rRegion.mDim, rRegion.mTag,
                     rremap_entries(out, rRegion.Stride()));
    return true;
}

void remap_regions(const Mesh& rIn, Mesh& rOut, const RegionRemap& rMaps) {
    const std::size_t nregions = rIn.NumRegions();
    for (std::size_t i = 0; i < nregions; ++i) {
        Region carried;
        if (remap_region(rIn, rOut, rIn.Region(i), rMaps, carried))
            rOut.AddRegion(std::move(carried));
    }
}

void warn_regions_dropped(const Mesh& rIn, const std::string& rOpName) {
    const std::size_t n = rIn.NumRegions();
    if (n == 0)
        return;
    log::warn(
        "{}: {} named region(s) dropped — the output cells and points are newly "
        "created and have no correspondence with the input's",
        rOpName, n);
    provenance_note("regions-dropped",
                    rOpName + ": " + std::to_string(n) +
                        " named region(s) dropped -- the output cells and points are "
                        "newly created and have no correspondence with the input's");
}

}  // namespace detail
}  // namespace meshioplusplus

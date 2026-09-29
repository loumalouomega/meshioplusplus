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
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/detail/region_remap.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/cell_edges.hpp"
#include "meshioplusplus/detail/cell_faces.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/facet_index.hpp"
#include "meshioplusplus/detail/geometry.hpp"
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

// ---------------------------------------------------------------------------
// Side facets by containment (v16.27.0)
//
// A side entry names a facet by (cell, local number). When an operation keeps
// the cell 1:1 the number still names the same facet; otherwise the facet is
// found again by what it is made of: a child facet belongs to a parent facet
// when every one of its nodes is the image of a node of that facet, or a new
// point lying on it (refinement); a source facet belongs to the output facet
// that contains every one of its surviving nodes, a removed node having to lie
// on it (coarsening, and a type change). Node identity decides wherever it
// can; geometry only places the points an operation created or removed.
// ---------------------------------------------------------------------------

/// One facet's nodes, and how many of them are corners (the leading ones).
struct RremapFacet {
    std::vector<std::int64_t> mNodes;
    std::size_t mNumCorners = 0;
};

/// Facet count of one cell: its type's for a rectangular block, its edge
/// count for a polygon, its face count for a polyhedron.
std::size_t rremap_facet_count(const Mesh& rMesh, const std::vector<std::int64_t>& rBases,
                               std::int64_t Global) {
    const auto [b, row] = global_to_block_row(rBases, Global);
    if (b == static_cast<std::size_t>(-1))
        return 0;
    const auto cb = rMesh.Cells(b);
    if (!cb.IsRagged())
        return region_num_facets(std::string(cb.Type()));
    const auto r = static_cast<std::size_t>(row);
    return cb.IsPolyhedron() ? cb.NumFaces(r) : cb.RowSize(r);
}

/// The nodes of facet @p Facet of global cell @p Global; false if it has none.
bool rremap_facet(const Mesh& rMesh, const std::vector<std::int64_t>& rBases, std::int64_t Global,
                  std::int64_t Facet, RremapFacet& rOut) {
    rOut.mNodes.clear();
    rOut.mNumCorners = 0;
    const auto [b, row] = global_to_block_row(rBases, Global);
    if (b == static_cast<std::size_t>(-1) || Facet < 0)
        return false;
    const auto cb = rMesh.Cells(b);
    const auto r = static_cast<std::size_t>(row);
    const auto f = static_cast<std::size_t>(Facet);
    if (cb.IsRagged()) {
        if (cb.IsPolyhedron()) {
            if (f >= cb.NumFaces(r))
                return false;
            const auto [p_face, face_size] = cb.Face(r, f);
            rOut.mNodes.assign(p_face, p_face + face_size);
        } else {
            // A polygon's edge k runs from its node k to node k + 1.
            const std::size_t n = cb.RowSize(r);
            if (f >= n)
                return false;
            const auto nodes = cb.Row(r);
            rOut.mNodes = {nodes[f], nodes[(f + 1) % n]};
        }
        rOut.mNumCorners = rOut.mNodes.size();
        return rOut.mNodes.size() >= 2;
    }
    const CellType type = cell_type_from_name(std::string(cb.Type()));
    const NDArray& conn = cb.Conn();
    const std::size_t k = cb.NodesPerCell();
    const int dim = cell_type_dimension(type);
    if (dim == 3) {
        const auto& faces = cell_faces(type);
        if (f >= faces.size())
            return false;
        for (std::size_t c = 0; c < faces[f].mNumNodes; ++c)
            rOut.mNodes.push_back(read_int(conn, r * k + faces[f].mNodes[c]));
        rOut.mNumCorners = faces[f].mNumCorners;
        return true;
    }
    if (dim == 2) {
        const auto& edges = cell_edges(type);
        if (f >= edges.size())
            return false;
        for (std::size_t c = 0; c < edges[f].mNumNodes; ++c)
            rOut.mNodes.push_back(read_int(conn, r * k + edges[f].mNodes[c]));
        rOut.mNumCorners = edges[f].mNumCorners;
        return true;
    }
    return false;
}

std::array<double, 3> rremap_coords(const Mesh& rMesh, std::int64_t Point) {
    std::array<double, 3> x{0.0, 0.0, 0.0};
    const NDArray& pts = rMesh.Points();
    const std::size_t dim = std::min<std::size_t>(rMesh.PointDim(), 3);
    for (std::size_t d = 0; d < dim; ++d)
        x[d] = read_double(pts, static_cast<std::size_t>(Point) * rMesh.PointDim() + d);
    return x;
}

/**
 * @brief Whether @p rQ lies on the facet whose corners are @p rCorners.
 *
 * Two corners: on the segment. More: on the plane through them (Newell's
 * normal, so a slightly warped quad still has one). The tolerance is relative
 * to the facet's size, since an operation's new points are averages of the
 * corners and sit on it to round-off.
 */
bool rremap_on_facet(const std::array<double, 3>& rQ,
                     const std::vector<std::array<double, 3>>& rCorners) {
    const std::size_t n = rCorners.size();
    if (n < 2)
        return false;
    double size = 0.0;
    for (std::size_t i = 1; i < n; ++i) {
        double d2 = 0.0;
        for (int d = 0; d < 3; ++d)
            d2 += (rCorners[i][d] - rCorners[0][d]) * (rCorners[i][d] - rCorners[0][d]);
        size = std::max(size, std::sqrt(d2));
    }
    if (size == 0.0)
        return false;
    const double tol = 1e-9 * size;
    if (n == 2) {
        std::array<double, 3> e{}, w{};
        double ee = 0.0, we = 0.0;
        for (int d = 0; d < 3; ++d) {
            e[d] = rCorners[1][d] - rCorners[0][d];
            w[d] = rQ[d] - rCorners[0][d];
            ee += e[d] * e[d];
            we += w[d] * e[d];
        }
        const double t = we / ee;
        if (t < -1e-9 || t > 1.0 + 1e-9)
            return false;
        double dist2 = 0.0;
        for (int d = 0; d < 3; ++d)
            dist2 += (w[d] - t * e[d]) * (w[d] - t * e[d]);
        return std::sqrt(dist2) <= tol;
    }
    std::array<double, 3> normal{0.0, 0.0, 0.0}, centre{0.0, 0.0, 0.0};
    for (std::size_t i = 0; i < n; ++i) {
        const auto& a = rCorners[i];
        const auto& b = rCorners[(i + 1) % n];
        normal[0] += (a[1] - b[1]) * (a[2] + b[2]);
        normal[1] += (a[2] - b[2]) * (a[0] + b[0]);
        normal[2] += (a[0] - b[0]) * (a[1] + b[1]);
        for (int d = 0; d < 3; ++d)
            centre[d] += a[d] / static_cast<double>(n);
    }
    const double len =
        std::sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
    if (len == 0.0)
        return false;
    double dist = 0.0;
    for (int d = 0; d < 3; ++d)
        dist += (rQ[d] - centre[d]) * normal[d] / len;
    return std::abs(dist) <= tol;
}

/// The corners of @p rFacet as coordinates of @p rMesh.
std::vector<std::array<double, 3>> rremap_corner_coords(const Mesh& rMesh,
                                                        const RremapFacet& rFacet) {
    std::vector<std::array<double, 3>> out;
    out.reserve(rFacet.mNumCorners);
    for (std::size_t i = 0; i < rFacet.mNumCorners; ++i)
        out.push_back(rremap_coords(rMesh, rFacet.mNodes[i]));
    return out;
}

/**
 * @brief Whether facet @p rSmall of @p rSmallMesh lies within facet @p rBig of
 * @p rBigMesh.
 *
 * @p SmallToBig maps a node of the small facet's mesh to the big one's, -1 for
 * a node with no counterpart there (created or removed by the operation). A
 * node with a counterpart must be a node of @p rBig -- unless @p Geometric,
 * the second pass for a merged face that dropped an interior vertex, where it
 * need only lie on it; a node without one must lie on it.
 */
template <class TMap>
bool rremap_facet_within(const Mesh& rSmallMesh, const RremapFacet& rSmall, const Mesh& rBigMesh,
                         const RremapFacet& rBig, TMap&& SmallToBig, bool Geometric) {
    std::vector<std::array<double, 3>> big_corners;
    for (std::int64_t node : rSmall.mNodes) {
        const std::int64_t there = SmallToBig(node);
        if (there >= 0 &&
            std::find(rBig.mNodes.begin(), rBig.mNodes.end(), there) != rBig.mNodes.end())
            continue;
        if (there >= 0 && !Geometric)
            return false;
        if (big_corners.empty())
            big_corners = rremap_corner_coords(rBigMesh, rBig);
        if (!rremap_on_facet(rremap_coords(rSmallMesh, node), big_corners))
            return false;
    }
    return true;
}

/// Output point -> input point: the inverse of `pPointMap` (identity for the
/// input's points when it is null), -1 for a point the operation created.
std::vector<std::int64_t> rremap_inverse_points(const Mesh& rIn, const Mesh& rOut,
                                                const RegionRemap& rMaps) {
    const std::size_t n_out = rOut.NumPoints();
    std::vector<std::int64_t> inv(n_out, -1);
    const auto n_in = static_cast<std::int64_t>(rIn.NumPoints());
    for (std::int64_t p = 0; p < n_in; ++p) {
        const std::int64_t q = rremap_lookup(rMaps.pPointMap, p);
        if (q >= 0 && static_cast<std::size_t>(q) < n_out && inv[static_cast<std::size_t>(q)] < 0)
            inv[static_cast<std::size_t>(q)] = p;
    }
    return inv;
}

/// How many input cells map onto each output cell (a merge has more than one).
std::vector<std::int32_t> rremap_preimage_counts(const RegionRemap& rMaps,
                                                 const std::vector<std::int64_t>& rInBases,
                                                 const std::vector<std::int64_t>& rOutBases) {
    std::vector<std::int32_t> counts(static_cast<std::size_t>(total_cells(rOutBases)), 0);
    std::vector<std::int64_t> children;
    for (std::int64_t c = 0; c < total_cells(rInBases); ++c) {
        rremap_children(rMaps, rInBases, rOutBases, c, children);
        for (std::int64_t ch : children)
            ++counts[static_cast<std::size_t>(ch)];
    }
    return counts;
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
        if (rMaps.mDropSideRegions) {
            log::warn(
                "{}: side region '{}' dropped — this operation does not preserve "
                "facet identity",
                rMaps.mOpName.empty() ? "operation" : rMaps.mOpName, rRegion.mName);
            return false;
        }
        // Side: see "Side facets by containment" above. The inverse point map
        // and the preimage counts are built only when an entry needs them.
        std::vector<std::int64_t> inv_points;
        std::vector<std::int32_t> preimages;
        bool have_inv = false, have_pre = false;
        auto out_to_in = [&](std::int64_t q) -> std::int64_t {
            return q >= 0 && static_cast<std::size_t>(q) < inv_points.size()
                       ? inv_points[static_cast<std::size_t>(q)]
                       : -1;
        };
        auto in_to_out = [&](std::int64_t p) -> std::int64_t {
            const std::int64_t q = rremap_lookup(rMaps.pPointMap, p);
            return q < n_out_points ? q : -1;
        };
        // Last resort: the facet still exists, by its (surviving) corners,
        // on another output cell -- a decimation collapse hands a boundary
        // edge to a neighbour. Built on first use.
        std::unique_ptr<FacetIndex> out_facets;
        std::vector<std::int64_t> mapped_corners;
        RremapFacet source, target;
        std::size_t lost = 0;
        out.reserve(n * 2);
        for (std::size_t k = 0; k < n; ++k) {
            const std::int64_t cell = entries[k * 2];
            const std::int64_t facet = entries[k * 2 + 1];
            const std::size_t before = out.size();
            rremap_children(rMaps, in_bases, out_bases, cell, children);
            if (!rremap_facet(rIn, in_bases, cell, facet, source)) {
                ++lost;
                continue;
            }
            if (children.empty()) {
                // The cell is gone; only the last resort below can place it.
            } else if (children.size() == 1) {
                const std::int64_t child = children[0];
                if (!have_pre && rMaps.mCellMapKind != CellMapKind::FirstChild) {
                    preimages = rremap_preimage_counts(rMaps, in_bases, out_bases);
                    have_pre = true;
                }
                const bool one_to_one =
                    !have_pre || preimages[static_cast<std::size_t>(child)] == 1;
                const std::string in_type = rremap_type_at(rIn, in_bases, cell);
                // The cell kept its identity: the local number still names
                // the facet (the only rule before v16.27.0) -- provided its
                // nodes are still the facet's, which a flipped or collapsed
                // cell's are not.
                if (one_to_one && in_type == rremap_type_at(rOut, out_bases, child) &&
                    rremap_facet(rOut, out_bases, child, facet, target) &&
                    target.mNodes.size() == source.mNodes.size() &&
                    rremap_facet_within(rIn, source, rOut, target, in_to_out,
                                        /*Geometric=*/false)) {
                    out.push_back(child);
                    out.push_back(facet);
                    continue;
                }
                // Merged or retyped: the output facet that contains it.
                const std::size_t nf = rremap_facet_count(rOut, out_bases, child);
                for (int pass = 0; pass < 2 && out.size() == before; ++pass)
                    for (std::size_t j = 0; j < nf; ++j) {
                        if (!rremap_facet(rOut, out_bases, child, static_cast<std::int64_t>(j),
                                          target))
                            continue;
                        if (rremap_facet_within(rIn, source, rOut, target, in_to_out,
                                                /*Geometric=*/pass == 1)) {
                            out.push_back(child);
                            out.push_back(static_cast<std::int64_t>(j));
                        }
                    }
            } else {
                // Refined: every child facet lying within it.
                if (!have_inv) {
                    inv_points = rremap_inverse_points(rIn, rOut, rMaps);
                    have_inv = true;
                }
                for (std::int64_t child : children) {
                    const std::size_t nf = rremap_facet_count(rOut, out_bases, child);
                    for (std::size_t j = 0; j < nf; ++j) {
                        if (!rremap_facet(rOut, out_bases, child, static_cast<std::int64_t>(j),
                                          target))
                            continue;
                        if (rremap_facet_within(rOut, target, rIn, source, out_to_in,
                                                /*Geometric=*/false)) {
                            out.push_back(child);
                            out.push_back(static_cast<std::int64_t>(j));
                        }
                    }
                }
            }
            if (out.size() == before) {
                mapped_corners.clear();
                for (std::size_t c = 0; c < source.mNumCorners; ++c)
                    mapped_corners.push_back(in_to_out(source.mNodes[c]));
                const bool all_survive = std::none_of(mapped_corners.begin(), mapped_corners.end(),
                                                      [](std::int64_t q) { return q < 0; });
                // Only across a FirstChild map (decimate, repair): under a
                // Direct one (crop, split) the facet's own cell is simply
                // gone, and handing an interface facet to the neighbour would
                // flip its orientation.
                if (rMaps.mCellMapKind == CellMapKind::FirstChild && all_survive &&
                    mapped_corners.size() <= 4) {
                    if (!out_facets)
                        out_facets = std::make_unique<FacetIndex>(rOut);
                    const FacetHit* p_hit =
                        out_facets->Find(mapped_corners.data(), mapped_corners.size());
                    if (p_hit) {
                        out.push_back(p_hit->mFirst.mCell);
                        out.push_back(p_hit->mFirst.mFacet);
                    }
                }
            }
            if (out.size() == before)
                ++lost;
        }
        if (lost > 0)
            log::warn(
                "{}: side region '{}' lost {} of its {} facet(s) — they have no "
                "counterpart in the output",
                rMaps.mOpName.empty() ? "operation" : rMaps.mOpName, rRegion.mName, lost, n);
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

void carry_regions_to_facet_mesh(const Mesh& rIn, Mesh& rOut,
                                 const std::vector<std::int64_t>& rOutParent,
                                 const std::vector<std::int64_t>& rOutFacet,
                                 const std::vector<std::int64_t>& rPointMap,
                                 const std::string& rOpName) {
    const std::size_t nregions = rIn.NumRegions();
    if (nregions == 0)
        return;
    // (input cell, facet) -> output cell, sorted for a binary search: built
    // once, and only when a Side region needs it.
    std::vector<std::array<std::int64_t, 3>> facet_to_out;
    auto lookup = [&](std::int64_t Cell, std::int64_t Facet) -> std::int64_t {
        if (facet_to_out.empty()) {
            facet_to_out.reserve(rOutParent.size());
            for (std::size_t c = 0; c < rOutParent.size() && c < rOutFacet.size(); ++c)
                facet_to_out.push_back({rOutParent[c], rOutFacet[c], static_cast<std::int64_t>(c)});
            std::sort(facet_to_out.begin(), facet_to_out.end());
        }
        const std::array<std::int64_t, 3> key{Cell, Facet, -1};
        const auto it = std::lower_bound(facet_to_out.begin(), facet_to_out.end(), key);
        return it != facet_to_out.end() && (*it)[0] == Cell && (*it)[1] == Facet ? (*it)[2] : -1;
    };
    const std::string op = rOpName.empty() ? "operation" : rOpName;
    for (std::size_t i = 0; i < nregions; ++i) {
        const Region& r_region = rIn.Region(i);
        const std::int64_t* entries = r_region.Entries();
        const std::size_t n = r_region.NumEntries();
        std::vector<std::int64_t> out;
        if (r_region.mKind == RegionKind::Point) {
            for (std::size_t k = 0; k < n; ++k) {
                const std::int64_t p = entries[k];
                if (p >= 0 && static_cast<std::size_t>(p) < rPointMap.size() &&
                    rPointMap[static_cast<std::size_t>(p)] >= 0)
                    out.push_back(rPointMap[static_cast<std::size_t>(p)]);
            }
            rOut.AddRegion(Region(r_region.mName, RegionKind::Point, r_region.mDim, r_region.mTag,
                                  rremap_entries(out, 1)));
        } else if (r_region.mKind == RegionKind::Side) {
            std::size_t lost = 0;
            for (std::size_t k = 0; k < n; ++k) {
                const std::int64_t c = lookup(entries[2 * k], entries[2 * k + 1]);
                if (c >= 0)
                    out.push_back(c);
                else
                    ++lost;
            }
            if (lost > 0)
                log::warn(
                    "{}: side region '{}' lost {} of its {} facet(s) — they are not on the "
                    "extracted boundary",
                    op, r_region.mName, lost, n);
            rOut.AddRegion(Region(r_region.mName, RegionKind::Cell, r_region.mDim, r_region.mTag,
                                  rremap_entries(out, 1)));
        } else {
            log::warn(
                "{}: cell region '{}' dropped — it names input cells, and the output "
                "holds only their facets",
                op, r_region.mName);
            provenance_note("regions-dropped", op + ": cell region '" + r_region.mName +
                                                   "' dropped -- the output holds only facets");
        }
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

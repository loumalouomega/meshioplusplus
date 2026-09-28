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

// feature_edges. See operations/feature_edges.hpp for the contract and
// detail/crease_edges.hpp for the test itself.

// System includes
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/feature_edges.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/cell_edges.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/face_mesh.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/region.hpp"

// Project includes (private, not installed)
#include "../detail/crease_edges.hpp"

namespace meshioplusplus {
namespace {

constexpr const char* kFePrefix = "meshio++: feature_edges: ";

// Per global cell, whether the region selects it (all cells when unnamed).
std::vector<std::uint8_t> fe_selection(const Mesh& rMesh, const std::string& rRegion,
                                       std::size_t NumCells) {
    if (rRegion.empty())
        return std::vector<std::uint8_t>(NumCells, 1);
    const std::size_t r = rMesh.FindRegion(rRegion, RegionKind::Cell);
    if (r == Mesh::npos) {
        std::string names;
        for (const std::string& n : rMesh.RegionNames())
            names += (names.empty() ? "" : ", ") + n;
        throw std::invalid_argument(std::string(kFePrefix) + "no cell region named '" + rRegion +
                                    "' (available: " + (names.empty() ? "none" : names) + ")");
    }
    std::vector<std::uint8_t> selected(NumCells, 0);
    const NDArray& entries = rMesh.Region(r).mEntries;
    for (std::size_t e = 0; e < entries.Size(); ++e) {
        const std::int64_t c = detail::read_int(entries, e);
        if (c >= 0 && static_cast<std::size_t>(c) < NumCells)
            selected[static_cast<std::size_t>(c)] = 1;
    }
    return selected;
}

// The face rings of the surface examined: the skin of the selected volume
// cells when there are any volume cells, else the selected 2-D cells.
void fe_rings(const Mesh& rMesh, const std::vector<std::uint8_t>& rSelected,
              std::vector<std::int64_t>& rStart, std::vector<std::int64_t>& rNodes) {
    const std::size_t n = rMesh.NumPoints();
    rStart.assign(1, 0);
    rNodes.clear();
    auto push = [&](const std::int64_t* pIds, std::size_t Size, bool Reverse) {
        for (std::size_t k = 0; k < Size; ++k)
            if (pIds[k] < 0 || static_cast<std::size_t>(pIds[k]) >= n)
                return;
        for (std::size_t k = 0; k < Size; ++k)
            rNodes.push_back(Reverse ? pIds[Size - 1 - k] : pIds[k]);
        rStart.push_back(static_cast<std::int64_t>(rNodes.size()));
    };

    const detail::GlobalFaces gf = detail::build_global_faces(rMesh);
    if (gf.NumCells() > 0) {
        // A face is on the skin of the selection when exactly one side of it
        // is selected; it is wound out of that side.
        auto in = [&](std::int64_t Compact) {
            return Compact >= 0 && rSelected[static_cast<std::size_t>(
                                       gf.mCellToGlobal[static_cast<std::size_t>(Compact)])] != 0;
        };
        for (std::size_t f = 0; f < gf.NumFaces(); ++f) {
            const bool own = in(gf.mOwner[f]);
            const bool nb = in(gf.mNeighbour[f]);
            if (own == nb)
                continue;
            push(gf.Face(f), gf.FaceSize(f), !own);
        }
        return;
    }

    const std::vector<std::int64_t> bases = detail::block_bases(rMesh);
    std::size_t bi = 0;
    std::vector<std::int64_t> ids;
    for (const auto cb : rMesh.CellRange()) {
        const std::size_t base = static_cast<std::size_t>(bases[bi++]);
        if (cb.IsPolyhedron())
            continue;
        const std::string type(cb.Type());
        const bool polygon = cb.IsRagged() || type.rfind("polygon", 0) == 0;
        std::vector<std::size_t> corners;
        if (!polygon) {
            const CellType ct = cell_type_from_name(type);
            if (cell_type_dimension(ct) != 2)
                continue;
            for (const detail::CellEdgeDef& ed : detail::cell_edges(ct))
                corners.push_back(ed.mNodes[0]);
            if (corners.size() < 3) {
                log::warn("feature_edges: cell type '{}' has no known edge topology; skipped",
                          type);
                continue;
            }
        }
        for (std::size_t c = 0; c < cb.NumCells(); ++c) {
            if (!rSelected[base + c])
                continue;
            ids.clear();
            if (polygon) {
                ids.assign(cb.Row(c), cb.Row(c) + cb.RowSize(c));
            } else {
                const std::size_t npc = cb.NodesPerCell();
                for (std::size_t k : corners)
                    ids.push_back(detail::read_int(cb.Conn(), c * npc + k));
            }
            if (ids.size() >= 3)
                push(ids.data(), ids.size(), false);
        }
    }
}

}  // namespace

FeatureEdgeResult feature_edges(const Mesh& rMesh, const FeatureEdgeOptions& rOptions) {
    if (!(rOptions.mFeatureAngleDeg >= 0.0 && rOptions.mFeatureAngleDeg <= 180.0))
        throw std::invalid_argument(std::string(kFePrefix) +
                                    "the feature angle must lie in [0, 180] degrees");
    const std::vector<std::int64_t> bases = detail::block_bases(rMesh);
    const std::size_t ncells = static_cast<std::size_t>(detail::total_cells(bases));
    const std::vector<std::uint8_t> selected = fe_selection(rMesh, rOptions.mRegion, ncells);

    std::vector<std::int64_t> start;
    std::vector<std::int64_t> nodes;
    fe_rings(rMesh, selected, start, nodes);

    // Flat coordinates, padded to 3-D.
    const std::size_t n = rMesh.NumPoints();
    const std::size_t dim = rMesh.PointDim();
    std::vector<double> xyz(n * 3, 0.0);
    {
        const NDArray& pts = rMesh.Points();
        parallel_for_bw(n, [&](std::size_t i) {
            for (std::size_t d = 0; d < dim && d < 3; ++d)
                xyz[i * 3 + d] = detail::read_double(pts, i * dim + d);
        });
    }
    const std::size_t nf = start.size() - 1;
    std::vector<double> normals(nf * 3);
    parallel_for(nf, [&](std::size_t f) {
        detail::ring_unit_normal(xyz.data(), nodes.data() + start[f],
                                 static_cast<std::size_t>(start[f + 1] - start[f]),
                                 &normals[f * 3]);
    });
    const std::vector<detail::CreaseEdge> edges =
        detail::crease_edges(start, nodes, normals, rOptions.mFeatureAngleDeg);

    FeatureEdgeResult result;
    std::vector<std::int64_t> conn;
    std::vector<std::int64_t> kind;
    std::vector<double> angle;
    for (const detail::CreaseEdge& e : edges) {
        result.mNumNonManifold += e.IsNonManifold() ? 1 : 0;
        result.mNumBoundary += e.IsBoundary() ? 1 : 0;
        result.mNumInconsistent += e.mInconsistent ? 1 : 0;
        result.mNumFeature += e.mSharp ? 1 : 0;
        FeatureEdgeKind k;
        if (e.IsNonManifold() && rOptions.mNonManifold)
            k = FeatureEdgeKind::NonManifold;
        else if (e.IsBoundary() && rOptions.mBoundary)
            k = FeatureEdgeKind::Boundary;
        else if (e.mInconsistent && rOptions.mInconsistent)
            k = FeatureEdgeKind::Inconsistent;
        else if (e.mSharp && rOptions.mFeature)
            k = FeatureEdgeKind::Feature;
        else
            continue;
        conn.push_back(e.mLo);
        conn.push_back(e.mHi);
        kind.push_back(static_cast<std::int64_t>(k));
        angle.push_back(e.mAngleDeg);
    }

    Mesh& out = result.mMesh;
    out.AssignPoints(detail::data_owned_copy(rMesh.Points()));
    const std::size_t ne = kind.size();
    NDArray line = NDArray::Uninit(DType::Int64, {ne, std::size_t{2}});
    std::copy(conn.begin(), conn.end(), line.As<std::int64_t>());
    out.AddCellBlock("line", std::move(line));
    // Int64, not the enum's Int32 backing type -- see the header comment: it
    // keeps the dtype identical whether NativeMesh/KratosMesh canonicalize it or
    // MeshioMesh leaves it alone.
    NDArray kind_a = NDArray::Uninit(DType::Int64, {ne});
    std::copy(kind.begin(), kind.end(), kind_a.As<std::int64_t>());
    NDArray angle_a = NDArray::Uninit(DType::Float64, {ne});
    std::copy(angle.begin(), angle.end(), angle_a.As<double>());
    std::vector<NDArray> kind_blocks;
    kind_blocks.push_back(std::move(kind_a));
    std::vector<NDArray> angle_blocks;
    angle_blocks.push_back(std::move(angle_a));
    out.AddCellData(kFeatureKindName, std::move(kind_blocks));
    out.AddCellData(kFeatureAngleName, std::move(angle_blocks));

    for (const std::string& name : rMesh.PointDataNames())
        out.AddPointData(name, detail::data_owned_copy(rMesh.PointData(name)));
    for (const std::string& name : rMesh.FieldDataNames())
        out.AddFieldData(name, detail::data_owned_copy(rMesh.FieldData(name)));
    bool dropped = false;
    for (std::size_t i = 0; i < rMesh.NumRegions(); ++i) {
        const meshioplusplus::Region& r = rMesh.Region(i);
        if (r.mKind == RegionKind::Point)
            out.AddRegion(r);
        else
            dropped = true;
    }
    if (dropped)
        log::warn(
            "feature_edges: Cell and Side regions name input cells, which the edge mesh does "
            "not have; they are dropped (Point regions are kept)");
    return result;
}

}  // namespace meshioplusplus

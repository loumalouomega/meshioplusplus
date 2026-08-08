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
// Polyhedral refinement: one polyhedral child per face, connected to a new
// interior point. Built entirely through detail::cell_rings/orient_rings, so
// it needs no per-type subdivision table -- see operations/subdivide.hpp for
// the contract.

// System includes
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/operations/subdivide.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/cell_faces.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/polyhedron.hpp"
#include "meshioplusplus/detail/region_remap.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {

namespace {

constexpr const char* kSubPrefix = "meshio++: subdivide: ";
constexpr const char* kSubParentName = "subdivide:parent_cell";

/// A staged output block: either an unchanged copy of an input block (any of
/// its three storage shapes), or a freshly built polyhedron block.
struct SubOutBlock {
    std::string mType;
    std::size_t mNodesPerCell = 0;  // rectangular pass-through only
    std::vector<std::int64_t> mConn;
    std::vector<std::vector<std::int64_t>> mPolygonRows;
    std::vector<std::vector<std::vector<std::int64_t>>> mPolyhedronCells;
    bool mIsRagged = false;
    bool mIsPolyhedron = false;
};

/// Stage an input block unchanged -- the pass-through case, one branch per
/// existing storage shape. Mirrors convert_cells.cpp's own staging helper;
/// not shared with it directly, since that one is file-private there too.
SubOutBlock sub_stage_passthrough(const Mesh::CellView& rBlock) {
    SubOutBlock out;
    out.mType = std::string(rBlock.Type());
    if (rBlock.IsPolyhedron()) {
        out.mIsRagged = true;
        out.mIsPolyhedron = true;
        out.mPolyhedronCells.resize(rBlock.NumCells());
        for (std::size_t c = 0; c < rBlock.NumCells(); ++c) {
            for (std::size_t f = 0; f < rBlock.NumFaces(c); ++f) {
                const auto face = rBlock.Face(c, f);
                out.mPolyhedronCells[c].emplace_back(face.first, face.first + face.second);
            }
        }
    } else if (rBlock.IsRagged()) {
        out.mIsRagged = true;
        out.mPolygonRows.resize(rBlock.NumCells());
        for (std::size_t c = 0; c < rBlock.NumCells(); ++c) {
            const std::int64_t* row = rBlock.Row(c);
            out.mPolygonRows[c].assign(row, row + rBlock.RowSize(c));
        }
    } else {
        const NDArray& conn = rBlock.Conn();
        const std::size_t npc = rBlock.NodesPerCell();
        out.mNodesPerCell = npc;
        out.mConn.resize(rBlock.NumCells() * npc);
        for (std::size_t c = 0; c < rBlock.NumCells(); ++c)
            for (std::size_t k = 0; k < npc; ++k)
                out.mConn[c * npc + k] = detail::read_int(conn, c * npc + k);
    }
    return out;
}

/// Add a staged block to `rOut`.
void sub_emit_block(Mesh& rOut, SubOutBlock& rBlock) {
    if (rBlock.mIsPolyhedron) {
        rOut.AddPolyhedronBlock(rBlock.mType, std::move(rBlock.mPolyhedronCells));
    } else if (rBlock.mIsRagged) {
        rOut.AddPolygonBlock(rBlock.mType, std::move(rBlock.mPolygonRows));
    } else {
        const std::size_t npc = rBlock.mNodesPerCell;
        const std::size_t ncells = npc == 0 ? 0 : rBlock.mConn.size() / npc;
        NDArray conn = NDArray::Uninit(DType::Int64, {ncells, npc});
        std::int64_t* dst = conn.As<std::int64_t>();
        for (std::size_t i = 0; i < rBlock.mConn.size(); ++i)
            dst[i] = rBlock.mConn[i];
        rOut.AddCellBlock(rBlock.mType, std::move(conn));
    }
}

/// Whether a block can be polyhedrally subdivided: a 3D type with a
/// `cell_faces` row (reduces to corners for any quadratic variant), or an
/// existing polyhedron block. 2D/1D blocks and the 3D Lagrange family (no
/// `cell_faces` row) pass through unchanged instead.
bool sub_eligible(const Mesh::CellView& rBlock) {
    if (rBlock.IsPolyhedron())
        return true;
    if (rBlock.IsRagged())
        return false;  // a ragged (polygon) block has no volume to subdivide
    const CellType t = cell_type_from_name(rBlock.Type());
    if (cell_type_dimension(t) != 3)
        return false;
    return !detail::cell_faces(t).empty();
}

}  // namespace

SubdivideResult subdivide(const Mesh& rMesh, const SubdivideOptions& rOptions) {
    const std::size_t nblocks = rMesh.NumCellBlocks();
    const std::size_t num_points = rMesh.NumPoints();
    const NDArray& points = rMesh.Points();
    const std::size_t pdim = rMesh.PointDim();

    // Every new point's coordinates are the plain average of a list of
    // existing node ids -- resolved in one batch pass after the main loop,
    // exactly convert_cells.cpp's own polyhedron-branch convention.
    std::vector<std::vector<std::int64_t>> new_point_src;

    std::vector<SubOutBlock> staged;
    staged.reserve(nblocks);
    std::vector<std::vector<std::int64_t>> first_child(nblocks);
    std::vector<std::vector<std::int64_t>> parent_of_child(nblocks);

    detail::CellRings rings;
    std::vector<detail::Vec3> coords;

    for (std::size_t b = 0; b < nblocks; ++b) {
        const auto cb = rMesh.Cells(b);
        const std::size_t ncells = cb.NumCells();
        std::vector<std::int64_t>& firsts = first_child[b];
        std::vector<std::int64_t>& parents = parent_of_child[b];
        firsts.resize(ncells);

        if (!sub_eligible(cb)) {
            staged.push_back(sub_stage_passthrough(cb));
            for (std::size_t c = 0; c < ncells; ++c) {
                firsts[c] = static_cast<std::int64_t>(c);
                parents.push_back(static_cast<std::int64_t>(c));
            }
            continue;
        }

        SubOutBlock out;
        out.mType = "polyhedron";
        out.mIsRagged = true;
        out.mIsPolyhedron = true;
        out.mPolyhedronCells.reserve(ncells);

        for (std::size_t c = 0; c < ncells; ++c) {
            firsts[c] = static_cast<std::int64_t>(out.mPolyhedronCells.size());
            if (!detail::cell_rings(cb, c, points, pdim, rings, coords))
                continue;  // no face topology at all: contributes nothing
            if (detail::orient_rings(rings, coords.data()) == detail::RingOrientation::Unorientable)
                throw std::invalid_argument(
                    std::string(kSubPrefix) + "cannot subdivide cell " + std::to_string(c) +
                    " of block '" + std::string(cb.Type()) +
                    "': its faces are not a closed orientable surface, so it bounds no volume");

            // One new interior point per cell: the plain average of the
            // cell's own corner nodes -- deliberately not poly_measure()'s
            // volume centroid, a different point. See the file docs.
            const std::int64_t apex = static_cast<std::int64_t>(num_points + new_point_src.size());
            new_point_src.push_back(rings.mNodes);

            for (std::size_t f = 0; f < rings.NumFaces(); ++f) {
                const std::uint32_t* ring = rings.Face(f);
                const std::size_t m = rings.FaceSize(f);

                std::vector<std::vector<std::int64_t>> child;
                child.reserve(1 + m);
                // The original face, unchanged (same global ids, same
                // winding), so a neighbouring cell across it still sees the
                // identical face -- this is what keeps the result conforming.
                std::vector<std::int64_t> face_nodes(m);
                for (std::size_t k = 0; k < m; ++k)
                    face_nodes[k] = rings.mNodes[ring[k]];
                child.push_back(std::move(face_nodes));
                // One new triangle per face edge, back to the apex.
                for (std::size_t k = 0; k < m; ++k) {
                    const std::int64_t a = rings.mNodes[ring[k]];
                    const std::int64_t nb = rings.mNodes[ring[(k + 1) % m]];
                    child.push_back({apex, a, nb});
                }
                out.mPolyhedronCells.push_back(std::move(child));
                parents.push_back(static_cast<std::int64_t>(c));
            }
        }
        staged.push_back(std::move(out));
    }

    Mesh out;
    if (new_point_src.empty()) {
        out.AssignPoints(detail::data_owned_copy(points));
    } else {
        const std::size_t dim = pdim;
        NDArray np = NDArray::Uninit(points.Dtype(), {num_points + new_point_src.size(), dim});
        std::memcpy(np.Data(), points.Data(), points.Nbytes());
        parallel_for_bw(new_point_src.size(), [&](std::size_t i) {
            const std::vector<std::int64_t>& src = new_point_src[i];
            for (std::size_t d = 0; d < dim; ++d) {
                double sum = 0.0;
                for (std::int64_t nid : src)
                    sum += detail::read_double(points, static_cast<std::size_t>(nid) * dim + d);
                detail::write_double(np, (num_points + i) * dim + d,
                                     sum / static_cast<double>(src.size()));
            }
        });
        out.AssignPoints(std::move(np));
    }

    for (SubOutBlock& block : staged)
        sub_emit_block(out, block);

    // point_data: originals copied verbatim; new (apex) points get the mean
    // of their source nodes' values, the same convention as coordinates.
    for (const std::string& name : rMesh.PointDataNames()) {
        const NDArray& a = rMesh.PointData(name);
        if (detail::rows(a) != num_points || new_point_src.empty()) {
            out.AddPointData(name, detail::data_owned_copy(a));
            continue;
        }
        const std::size_t ncomp = num_points == 0 ? 0 : a.Size() / num_points;
        std::vector<std::size_t> shape = a.Shape();
        shape[0] = num_points + new_point_src.size();
        NDArray b = NDArray::Uninit(a.Dtype(), std::move(shape));
        std::memcpy(b.Data(), a.Data(), a.Nbytes());
        parallel_for_bw(new_point_src.size(), [&](std::size_t i) {
            const std::vector<std::int64_t>& src = new_point_src[i];
            for (std::size_t k = 0; k < ncomp; ++k) {
                double sum = 0.0;
                for (std::int64_t nid : src)
                    sum += detail::read_double(a, static_cast<std::size_t>(nid) * ncomp + k);
                detail::write_double(b, (num_points + i) * ncomp + k,
                                     sum / static_cast<double>(src.size()));
            }
        });
        out.AddPointData(name, std::move(b));
    }

    // cell_data: replicate each parent's row to its children, by raw byte copy
    // (dtype-agnostic, no precision loss) -- convert_cells.cpp's own pattern
    // for the identical "each child inherits its parent's row" shape. An
    // array whose block count does not match the mesh is not per-cell data,
    // so it is copied verbatim rather than mangled.
    for (const std::string& name : rMesh.CellDataNames()) {
        const std::size_t ndata = rMesh.CellDataNumBlocks(name);
        std::vector<NDArray> blocks;
        blocks.reserve(ndata);
        for (std::size_t b = 0; b < ndata; ++b) {
            const NDArray& src = rMesh.CellData(name, b);
            const std::size_t in_rows = detail::rows(src);
            if (ndata != nblocks || in_rows == 0 || parent_of_child[b].empty()) {
                blocks.push_back(detail::data_owned_copy(src));
                continue;
            }
            const std::vector<std::int64_t>& parents = parent_of_child[b];
            std::vector<std::size_t> shape = src.Shape();
            shape[0] = parents.size();
            const std::size_t row_bytes = src.Nbytes() / in_rows;
            NDArray dst = NDArray::Uninit(src.Dtype(), std::move(shape));
            const std::byte* p_src = src.Data();
            std::byte* p_dst = dst.Data();
            parallel_for_bw(parents.size(), [&](std::size_t i) {
                std::memcpy(p_dst + i * row_bytes,
                            p_src + static_cast<std::size_t>(parents[i]) * row_bytes, row_bytes);
            });
            blocks.push_back(std::move(dst));
        }
        out.AddCellData(name, std::move(blocks));
    }
    for (const std::string& name : rMesh.FieldDataNames())
        out.AddFieldData(name, detail::data_owned_copy(rMesh.FieldData(name)));

    if (rOptions.mRecordParentIds) {
        std::vector<NDArray> blocks;
        blocks.reserve(nblocks);
        for (const std::vector<std::int64_t>& parents : parent_of_child) {
            NDArray dst(DType::Int64, {parents.size()});
            std::int64_t* p = dst.As<std::int64_t>();
            for (std::size_t i = 0; i < parents.size(); ++i)
                p[i] = parents[i];
            blocks.push_back(std::move(dst));
        }
        out.AddCellData(kSubParentName, std::move(blocks));
    }

    SubdivideResult res;
    res.mMesh = std::move(out);
    res.mCellMaps.reserve(nblocks);
    for (std::vector<std::int64_t>& firsts : first_child) {
        NDArray m(DType::Int64, {firsts.size()});
        std::int64_t* p = m.As<std::int64_t>();
        for (std::size_t i = 0; i < firsts.size(); ++i)
            p[i] = firsts[i];
        res.mCellMaps.push_back(std::move(m));
    }

    detail::RegionRemap rmap;
    rmap.mCellMapKind = detail::CellMapKind::FirstChild;
    rmap.pCellMaps = &res.mCellMaps;
    rmap.mOpName = "subdivide";
    detail::remap_regions(rMesh, res.mMesh, rmap);

    return res;
}

}  // namespace meshioplusplus

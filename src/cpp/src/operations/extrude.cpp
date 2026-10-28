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
// Sweep a surface into a volume or a curve into a surface, along a vector
// (extrude, roadmap 3.1.2) or around an axis (revolve, 3.1.3), on the layer
// sweep `boundary_layer` grows its layers with. See operations/extrude.hpp.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <functional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "../detail/layer_sweep.hpp"
#include "meshioplusplus/detail/cell_edges.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/facet_index.hpp"
#include "meshioplusplus/detail/subset.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/operations/extrude.hpp"
#include "meshioplusplus/operations/transform.hpp"
#include "meshioplusplus/region.hpp"

namespace meshioplusplus {

namespace {

using ExVec = std::array<double, 3>;

ExVec ex_sub(const ExVec& rA, const ExVec& rB) {
    return {rA[0] - rB[0], rA[1] - rB[1], rA[2] - rB[2]};
}
ExVec ex_cross(const ExVec& rA, const ExVec& rB) {
    return {rA[1] * rB[2] - rA[2] * rB[1], rA[2] * rB[0] - rA[0] * rB[2],
            rA[0] * rB[1] - rA[1] * rB[0]};
}
double ex_dot(const ExVec& rA, const ExVec& rB) {
    return rA[0] * rB[0] + rA[1] * rB[1] + rA[2] * rB[2];
}
double ex_norm(const ExVec& rA) {
    return std::sqrt(ex_dot(rA, rA));
}

// The profile: the cells that are swept, and what else the mesh holds.
struct ExProfile {
    bool mCurve = false;                               // lines sweep to quads
    std::vector<std::array<std::int64_t, 4>> mFaces;   // a triangle has -1 last
    std::vector<std::array<std::int64_t, 2>> mLines;
    std::vector<std::size_t> mSourceBlock;             // per element
    std::vector<std::int64_t> mSourceRow;              // per element
    std::vector<std::vector<std::int64_t>> mElementOf; // per input block and row, or -1
    std::vector<char> mIsProfileBlock;                 // per input block
    std::vector<char> mIsLowerBlock;                   // per input block
    std::vector<double> mXyz;                          // three coordinates per point
    std::vector<char> mUsed;                           // per point, used by an element
    std::size_t NumElements() const { return mCurve ? mLines.size() : mFaces.size(); }
    ExVec Point(std::int64_t P) const {
        const auto p = static_cast<std::size_t>(P);
        return {mXyz[3 * p], mXyz[3 * p + 1], mXyz[3 * p + 2]};
    }
};

ExProfile ex_gather(const Mesh& rMesh, const std::string& rOp) {
    ExProfile pr;
    const std::size_t n = rMesh.NumPoints();
    const std::size_t dim = rMesh.PointDim();
    pr.mXyz.assign(3 * n, 0.0);
    for (std::size_t p = 0; p < n; ++p)
        for (std::size_t d = 0; d < dim && d < 3; ++d)
            pr.mXyz[3 * p + d] = detail::read_double(rMesh.Points(), p * dim + d);
    pr.mUsed.assign(n, 0);
    bool has2 = false, has1 = false;
    for (const auto cb : rMesh.CellRange()) {
        const std::string type(cb.Type());
        const int d = cell_type_dimension(cell_type_from_name(type));
        if (d == 3 || cb.IsPolyhedron())
            throw Unsupported("meshio++: " + rOp + ": volume cells ('" + type +
                              "') cannot be swept; pass a surface or a curve");
        if (cb.NumCells() == 0)
            continue;
        has2 = has2 || d == 2 || (cb.IsRagged() && !cb.IsPolyhedron());
        has1 = has1 || d == 1;
    }
    const int profile_dim = has2 ? 2 : (has1 ? 1 : 0);
    if (profile_dim == 0)
        throw std::invalid_argument("meshio++: " + rOp +
                                    ": no triangle, quad or line cells to sweep");
    pr.mCurve = profile_dim == 1;
    const std::size_t nb = rMesh.NumCellBlocks();
    pr.mElementOf.resize(nb);
    pr.mIsProfileBlock.assign(nb, 0);
    pr.mIsLowerBlock.assign(nb, 0);
    for (std::size_t b = 0; b < nb; ++b) {
        const auto cb = rMesh.Cells(b);
        const std::string type(cb.Type());
        const CellType ct = cell_type_from_name(type);
        const int d = cb.IsRagged() ? 2 : cell_type_dimension(ct);
        pr.mElementOf[b].assign(cb.NumCells(), -1);
        if (d < profile_dim) {
            pr.mIsLowerBlock[b] = 1;
            continue;
        }
        const CellType want = pr.mCurve ? CellType::Line : CellType::Triangle;
        if (cb.IsRagged() || (ct != want && !(ct == CellType::Quad && !pr.mCurve)))
            throw Unsupported("meshio++: " + rOp + ": cannot sweep '" + type + "' cells (" +
                              (pr.mCurve ? std::string("line")
                                         : std::string("triangle and quad")) +
                              " only; convert_cells linearizes or triangulates them)");
        pr.mIsProfileBlock[b] = 1;
        const std::size_t npc = cb.NodesPerCell();
        const NDArray& conn = cb.Conn();
        for (std::size_t c = 0; c < cb.NumCells(); ++c) {
            pr.mElementOf[b][c] = static_cast<std::int64_t>(pr.NumElements());
            pr.mSourceBlock.push_back(b);
            pr.mSourceRow.push_back(static_cast<std::int64_t>(c));
            std::array<std::int64_t, 4> v{-1, -1, -1, -1};
            for (std::size_t k = 0; k < npc; ++k) {
                v[k] = detail::read_int(conn, c * npc + k);
                if (v[k] < 0 || static_cast<std::size_t>(v[k]) >= n)
                    throw std::invalid_argument("meshio++: " + rOp + ": cell " +
                                                std::to_string(c) + " of block " +
                                                std::to_string(b) + " names a missing point");
                pr.mUsed[static_cast<std::size_t>(v[k])] = 1;
            }
            if (pr.mCurve)
                pr.mLines.push_back({v[0], v[1]});
            else
                pr.mFaces.push_back(v);
        }
    }
    return pr;
}

// A face's Newell normal (twice its area vector).
ExVec ex_face_normal(const ExProfile& rPr, const std::array<std::int64_t, 4>& rF) {
    const int k = rF[3] >= 0 ? 4 : 3;
    ExVec nrm{0.0, 0.0, 0.0};
    for (int i = 0; i < k; ++i) {
        const ExVec a = rPr.Point(rF[static_cast<std::size_t>(i)]);
        const ExVec b = rPr.Point(rF[static_cast<std::size_t>((i + 1) % k)]);
        nrm[0] += (a[1] - b[1]) * (a[2] + b[2]);
        nrm[1] += (a[2] - b[2]) * (a[0] + b[0]);
        nrm[2] += (a[0] - b[0]) * (a[1] + b[1]);
    }
    return nrm;
}

void ex_flip(std::array<std::int64_t, 4>& rF) {
    if (rF[3] >= 0)
        std::swap(rF[1], rF[3]);
    else
        std::swap(rF[1], rF[2]);
}

// How the sweep is laid out: the offsets the collapse rule reads, where each
// new point goes, and the names.
struct ExPlan {
    std::string mOp;
    int mLayers = 1;
    std::vector<double> mOffsets;  // point-major, n * L
    bool mPeriodic = false;
    bool mCollapsedEdgeWedges = false;
    std::function<void(std::size_t, int, double*)> mPlace;  // (point, level, out xyz)
    std::string mLayerName;
    bool mRecordLayers = false;
    std::string mFirstName, mLastName;
    std::int64_t mMaxCells = 0;
};

// One output cell list, sorted layer by layer.
struct ExOutBlock {
    const char* pType;
    int mNodes;
    std::vector<detail::LayerSweepCell> mCells;
};

void ex_sweep(const Mesh& rMesh, const ExProfile& rPr, const ExPlan& rPlan, SweepResult& rRes) {
    const std::size_t n = rMesh.NumPoints();
    const int L = rPlan.mLayers;
    const std::size_t ne = rPr.NumElements();
    if (rPlan.mMaxCells > 0 &&
        static_cast<double>(ne) * static_cast<double>(L) > static_cast<double>(rPlan.mMaxCells))
        throw std::invalid_argument("meshio++: " + rPlan.mOp + ": " + std::to_string(ne) +
                                    " profile cells times " + std::to_string(L) +
                                    " layers exceed max_cells = " +
                                    std::to_string(rPlan.mMaxCells));
    detail::LayerSweepInput in;
    in.pPoints = rPr.mXyz.data();
    in.mNumPoints = static_cast<std::int64_t>(n);
    in.mDirections.assign(3 * n, 0.0);
    in.mOffsets = rPlan.mOffsets;
    in.mLayers = L;
    in.mPeriodic = rPlan.mPeriodic;
    in.mCollapsedEdgeWedges = rPlan.mCollapsedEdgeWedges;
    if (rPr.mCurve)
        in.mEdges = rPr.mLines;
    else
        in.mFaces = rPr.mFaces;
    detail::LayerSweepResult sw = detail::layer_sweep(in);
    const auto level_id = [&](std::int64_t P, int K) -> std::int64_t {
        return K == 0 ? P
                      : sw.mLevelIds[static_cast<std::size_t>(K - 1) * n +
                                     static_cast<std::size_t>(P)];
    };
    // Put each level's own points where the sweep belongs.
    for (int k = 1; k <= L; ++k)
        for (std::size_t p = 0; p < n; ++p) {
            const std::int64_t id = level_id(static_cast<std::int64_t>(p), k);
            if (id >= static_cast<std::int64_t>(n) &&
                id != level_id(static_cast<std::int64_t>(p), k - 1))
                rPlan.mPlace(p, k, &sw.mPoints[3 * static_cast<std::size_t>(id)]);
        }

    std::vector<ExOutBlock> blocks;
    if (rPr.mCurve) {
        blocks = {{"quad", 4, std::move(sw.mSideQuads)},
                  {"triangle", 3, std::move(sw.mSideTriangles)}};
    } else {
        blocks = {{"wedge", 6, std::move(sw.mWedges)},
                  {"hexahedron", 8, std::move(sw.mHexahedra)},
                  {"pyramid", 5, std::move(sw.mPyramids)},
                  {"tetra", 4, std::move(sw.mTetrahedra)}};
    }
    for (auto& r_b : blocks)
        std::stable_sort(r_b.mCells.begin(), r_b.mCells.end(),
                         [](const detail::LayerSweepCell& rA, const detail::LayerSweepCell& rB) {
                             return rA.mLayer < rB.mLayer;
                         });
    blocks.erase(std::remove_if(blocks.begin(), blocks.end(),
                                [](const ExOutBlock& rB) { return rB.mCells.empty(); }),
                 blocks.end());

    // The mesh.
    Mesh& out = rRes.mMesh;
    const std::size_t npts = sw.mPoints.size() / 3;
    {
        NDArray pts(DType::Float64, {npts, 3});
        std::copy(sw.mPoints.begin(), sw.mPoints.end(), pts.As<double>());
        out.AssignPoints(std::move(pts));
    }
    std::vector<std::int64_t> out_element;  // per output cell (global)
    rRes.mMinScaledJacobian = rPr.mCurve ? 0.0 : 1.0;
    for (const ExOutBlock& r_b : blocks) {
        const auto k = static_cast<std::size_t>(r_b.mNodes);
        NDArray conn(DType::Int64, {r_b.mCells.size(), k});
        std::int64_t* q = conn.As<std::int64_t>();
        for (std::size_t c = 0; c < r_b.mCells.size(); ++c) {
            for (std::size_t i = 0; i < k; ++i)
                q[c * k + i] = r_b.mCells[c].mNodes[i];
            out_element.push_back(r_b.mCells[c].mFace);
            if (!rPr.mCurve)
                rRes.mMinScaledJacobian =
                    std::min(rRes.mMinScaledJacobian,
                             detail::layer_cell_scaled_jacobian(sw.mPoints.data(), r_b.mCells[c],
                                                                r_b.mNodes));
        }
        const std::string type(r_b.pType);
        const auto count = static_cast<std::int64_t>(r_b.mCells.size());
        if (type == "wedge")
            rRes.mNumWedges = count;
        else if (type == "hexahedron")
            rRes.mNumHexahedra = count;
        else if (type == "pyramid")
            rRes.mNumPyramids = count;
        else if (type == "tetra")
            rRes.mNumTetrahedra = count;
        else if (type == "quad")
            rRes.mNumQuads = count;
        else
            rRes.mNumTriangles = count;
        rRes.mNumCells += count;
        out.AddCellBlock(type, std::move(conn));
    }

    // Data. Point data per copy; cell data per layer, gathered across the
    // profile's blocks.
    for (const std::string& name : rMesh.PointDataNames())
        out.AddPointData(name, detail::subset_gather_rows(rMesh.PointData(name), sw.mOrigin));
    bool lower_data = false;
    for (const std::string& name : rMesh.CellDataNames()) {
        const NDArray* p_first = nullptr;
        bool ok = true;
        std::size_t cols = 0;
        for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
            const std::size_t nc = rMesh.Cells(b).NumCells();
            if (rPr.mIsLowerBlock[b]) {
                lower_data = lower_data || (nc > 0 && b < rMesh.CellDataNumBlocks(name));
                continue;
            }
            if (!rPr.mIsProfileBlock[b] || nc == 0)
                continue;
            if (b >= rMesh.CellDataNumBlocks(name)) {
                ok = false;
                break;
            }
            const NDArray& a = rMesh.CellData(name, b);
            const std::size_t c = a.Size() / nc;
            if (p_first == nullptr) {
                p_first = &a;
                cols = c;
            } else if (a.Dtype() != p_first->Dtype() || c != cols) {
                ok = false;
            }
        }
        if (!ok || p_first == nullptr) {
            if (!ok)
                rRes.mWarnings.push_back("cell data '" + name +
                                         "' differs between the profile's blocks and was dropped");
            continue;
        }
        const std::size_t row_bytes = cols * dtype_size(p_first->Dtype());
        std::vector<std::size_t> shape = p_first->Shape();
        std::vector<NDArray> arrays;
        std::size_t g = 0;
        for (const ExOutBlock& r_b : blocks) {
            if (shape.empty())
                shape = {r_b.mCells.size()};
            else
                shape[0] = r_b.mCells.size();
            NDArray a = NDArray::Uninit(p_first->Dtype(), shape);
            for (std::size_t c = 0; c < r_b.mCells.size(); ++c, ++g) {
                const auto e = static_cast<std::size_t>(out_element[g]);
                const NDArray& src = rMesh.CellData(name, rPr.mSourceBlock[e]);
                if (row_bytes > 0)
                    std::memcpy(a.Data() + c * row_bytes,
                                src.Data() + static_cast<std::size_t>(rPr.mSourceRow[e]) * row_bytes,
                                row_bytes);
            }
            arrays.push_back(std::move(a));
        }
        out.AddCellData(name, std::move(arrays));
    }
    if (lower_data)
        rRes.mWarnings.push_back(
            "cell data of the cells that are not swept (lower-dimensional ones) is not carried");
    if (rPlan.mRecordLayers) {
        std::vector<NDArray> arrays;
        for (const ExOutBlock& r_b : blocks) {
            NDArray a(DType::Int32, {r_b.mCells.size()});
            for (std::size_t c = 0; c < r_b.mCells.size(); ++c)
                a.As<std::int32_t>()[c] = r_b.mCells[c].mLayer;
            arrays.push_back(std::move(a));
        }
        out.AddCellData(rPlan.mLayerName, std::move(arrays));
    }
    for (const std::string& name : rMesh.FieldDataNames()) {
        NDArray a = rMesh.FieldData(name);
        a.MakeOwned();
        out.AddFieldData(name, std::move(a));
    }

    // Regions. Every facet is found by its corners among the output's facets;
    // an interior one (two owners) is given to the cell swept from WANT.
    detail::FacetIndexOptions fio;
    fio.mSolidFaces = !rPr.mCurve;
    fio.mSurfaceEdges = rPr.mCurve;
    const detail::FacetIndex index(out, fio);
    const auto find = [&](std::vector<std::int64_t> Corners, std::int64_t Want,
                          std::vector<std::int64_t>& rEntries) {
        std::sort(Corners.begin(), Corners.end());
        Corners.erase(std::unique(Corners.begin(), Corners.end()), Corners.end());
        const detail::FacetHit* hit = index.Find(Corners.data(), Corners.size());
        if (hit == nullptr || hit->mCount == 0)
            return false;
        detail::FacetOwner o = hit->mFirst;
        if (hit->mCount > 1 && Want >= 0 &&
            out_element[static_cast<std::size_t>(o.mCell)] != Want &&
            out_element[static_cast<std::size_t>(hit->mSecond.mCell)] == Want)
            o = hit->mSecond;
        rEntries.push_back(o.mCell);
        rEntries.push_back(o.mFacet);
        return true;
    };
    // What a profile facet (or a lower-dimensional cell) with corners C sweeps
    // to, layer by layer; false when one of its pieces is not found.
    const auto swept = [&](const std::vector<std::int64_t>& rC, std::int64_t Want,
                           std::vector<std::int64_t>& rEntries) {
        bool all = true;
        for (int k = 1; k <= L; ++k) {
            std::vector<std::int64_t> c;
            for (const std::int64_t v : rC) {
                c.push_back(level_id(v, k - 1));
                c.push_back(level_id(v, k));
            }
            std::vector<std::int64_t> u = c;
            std::sort(u.begin(), u.end());
            u.erase(std::unique(u.begin(), u.end()), u.end());
            if (u.size() <= rC.size())
                continue;  // collapsed to nothing in this layer
            all = find(std::move(c), Want, rEntries) && all;
        }
        return all;
    };
    // Per element, its output cells.
    std::vector<std::vector<std::int64_t>> cells_of(ne);
    for (std::size_t g = 0; g < out_element.size(); ++g)
        cells_of[static_cast<std::size_t>(out_element[g])].push_back(static_cast<std::int64_t>(g));
    const std::vector<std::int64_t> in_bases = detail::block_bases(rMesh);
    const auto corners_of_cell = [&](std::size_t B, std::int64_t Row) {
        const auto cb = rMesh.Cells(B);
        std::vector<std::int64_t> c;
        const CellType ct = cell_type_from_name(std::string(cb.Type()));
        if (cb.IsRagged()) {
            c.assign(cb.Row(static_cast<std::size_t>(Row)),
                     cb.Row(static_cast<std::size_t>(Row)) + cb.RowSize(static_cast<std::size_t>(Row)));
            return c;
        }
        // Corners only: a line's two ends, a vertex's node.
        const std::size_t npc = cb.NodesPerCell();
        const std::size_t nc = cell_type_dimension(ct) == 1 ? std::min<std::size_t>(npc, 2)
                                                             : std::min<std::size_t>(npc, 1);
        for (std::size_t k = 0; k < nc; ++k)
            c.push_back(detail::read_int(cb.Conn(), static_cast<std::size_t>(Row) * npc + k));
        return c;
    };
    const int out_dim = rPr.mCurve ? 2 : 3;
    std::vector<std::string> taken;
    for (std::size_t i = 0; i < rMesh.NumRegions(); ++i) {
        const Region& r_in = rMesh.Region(i);
        const std::int64_t* e = r_in.Entries();
        const std::size_t nent = r_in.NumEntries();
        std::vector<std::int64_t> flat, side;
        std::int64_t lost = 0;
        RegionKind kind = r_in.mKind;
        if (r_in.mKind == RegionKind::Point) {
            for (std::size_t k = 0; k < nent; ++k) {
                if (e[k] < 0 || static_cast<std::size_t>(e[k]) >= n) {
                    ++lost;
                    continue;
                }
                for (int l = 0; l <= L; ++l)
                    flat.push_back(level_id(e[k], l));
            }
        } else {
            const bool is_side = r_in.mKind == RegionKind::Side;
            for (std::size_t k = 0; k < nent; ++k) {
                const std::int64_t g = is_side ? e[2 * k] : e[k];
                const auto [b, row] = detail::global_to_block_row(in_bases, g);
                if (b >= rMesh.NumCellBlocks()) {
                    ++lost;
                    continue;
                }
                const std::int64_t el = rPr.mElementOf[b][static_cast<std::size_t>(row)];
                if (!is_side && el >= 0) {
                    const auto& r_cells = cells_of[static_cast<std::size_t>(el)];
                    flat.insert(flat.end(), r_cells.begin(), r_cells.end());
                } else if (!is_side && rPr.mIsLowerBlock[b]) {
                    if (!swept(corners_of_cell(b, row), -1, side))
                        ++lost;
                } else if (is_side && el >= 0 && !rPr.mCurve) {
                    const std::int64_t f = e[2 * k + 1];
                    const auto& r_edges =
                        detail::cell_edges(cell_type_from_name(std::string(rMesh.Cells(b).Type())));
                    if (f < 0 || static_cast<std::size_t>(f) >= r_edges.size()) {
                        ++lost;
                        continue;
                    }
                    const auto cb = rMesh.Cells(b);
                    const std::size_t npc = cb.NodesPerCell();
                    std::vector<std::int64_t> c;
                    for (std::size_t j = 0; j < 2; ++j)
                        c.push_back(detail::read_int(
                            cb.Conn(), static_cast<std::size_t>(row) * npc +
                                           r_edges[static_cast<std::size_t>(f)].mNodes[j]));
                    if (!swept(c, el, side))
                        ++lost;
                } else {
                    ++lost;
                }
            }
            if (!is_side && flat.empty() && !side.empty())
                kind = RegionKind::Side;
            else if (!is_side && !side.empty())
                lost += static_cast<std::int64_t>(side.size() / 2);  // a mixed region keeps its cells
            if (is_side || kind == RegionKind::Side)
                flat = std::move(side);
        }
        if (lost > 0)
            rRes.mWarnings.push_back("region '" + r_in.mName + "': " + std::to_string(lost) +
                                     " entr(ies) could not be carried and were dropped");
        const int dim = r_in.mDim >= 0 && r_in.mKind != RegionKind::Point ? r_in.mDim + 1
                                                                          : r_in.mDim;
        const bool pair = kind == RegionKind::Side;
        NDArray entries = pair ? NDArray(DType::Int64, {flat.size() / 2, 2})
                               : NDArray(DType::Int64, {flat.size()});
        std::copy(flat.begin(), flat.end(), entries.As<std::int64_t>());
        out.AddRegion(Region(r_in.mName, kind, dim, r_in.mTag, std::move(entries)));
        taken.push_back(r_in.mName);
    }

    // The generated side regions.
    const auto add_generated = [&](const std::string& rName, std::vector<std::int64_t> Entries) {
        if (rName.empty())
            return;
        if (std::find(taken.begin(), taken.end(), rName) != taken.end()) {
            rRes.mWarnings.push_back("the input has a region '" + rName +
                                     "', so the generated one of that name is not added");
            return;
        }
        NDArray a(DType::Int64, {Entries.size() / 2, 2});
        std::copy(Entries.begin(), Entries.end(), a.As<std::int64_t>());
        out.AddRegion(Region(rName, RegionKind::Side, out_dim - 1, -1, std::move(a)));
    };
    if (!rPlan.mPeriodic) {
        std::vector<std::int64_t> first, last;
        if (rPr.mCurve) {
            for (std::size_t l = 0; l < rPr.mLines.size(); ++l) {
                const auto& r_l = rPr.mLines[l];
                find({r_l[0], r_l[1]}, static_cast<std::int64_t>(l), first);
                find({level_id(r_l[0], L), level_id(r_l[1], L)}, static_cast<std::int64_t>(l),
                     last);
            }
        } else {
            const auto faces = [&](const std::vector<detail::LayerSweepCell>& rList, int K,
                                   std::vector<std::int64_t>& rInto) {
                for (const auto& r_c : rList) {
                    std::vector<std::int64_t> c(r_c.mNodes.begin(),
                                                r_c.mNodes.begin() + (K == 4 ? 4 : 3));
                    find(std::move(c), r_c.mFace, rInto);
                }
            };
            faces(sw.mWallTriangles, 3, first);
            faces(sw.mWallQuads, 4, first);
            faces(sw.mCapTriangles, 3, last);
            faces(sw.mCapQuads, 4, last);
        }
        add_generated(rPlan.mFirstName, std::move(first));
        add_generated(rPlan.mLastName, std::move(last));
    }
    {
        // The open border: edges of one face (in its winding), or ends of one line.
        std::vector<std::int64_t> sides;
        if (rPr.mCurve) {
            std::vector<std::int64_t> uses(n, 0), line_of(n, -1);
            for (std::size_t l = 0; l < rPr.mLines.size(); ++l)
                for (const std::int64_t v : rPr.mLines[l]) {
                    ++uses[static_cast<std::size_t>(v)];
                    line_of[static_cast<std::size_t>(v)] = static_cast<std::int64_t>(l);
                }
            for (std::size_t p = 0; p < n; ++p)
                if (uses[p] == 1)
                    swept({static_cast<std::int64_t>(p)}, line_of[p], sides);
        } else {
            struct ExEdge {
                std::int64_t mLo, mHi, mFace;
                std::array<std::int64_t, 2> mDirected;
            };
            std::vector<ExEdge> edges;
            for (std::size_t f = 0; f < rPr.mFaces.size(); ++f) {
                const auto& r_f = rPr.mFaces[f];
                const int k = r_f[3] >= 0 ? 4 : 3;
                for (int i = 0; i < k; ++i) {
                    const std::int64_t a = r_f[static_cast<std::size_t>(i)];
                    const std::int64_t b = r_f[static_cast<std::size_t>((i + 1) % k)];
                    edges.push_back({std::min(a, b), std::max(a, b), static_cast<std::int64_t>(f),
                                     {a, b}});
                }
            }
            std::stable_sort(edges.begin(), edges.end(), [](const ExEdge& rA, const ExEdge& rB) {
                return rA.mLo != rB.mLo ? rA.mLo < rB.mLo : rA.mHi < rB.mHi;
            });
            for (std::size_t i = 0; i < edges.size();) {
                std::size_t j = i + 1;
                while (j < edges.size() && edges[j].mLo == edges[i].mLo &&
                       edges[j].mHi == edges[i].mHi)
                    ++j;
                if (j == i + 1)
                    swept({edges[i].mDirected[0], edges[i].mDirected[1]}, edges[i].mFace, sides);
                i = j;
            }
        }
        add_generated("side", std::move(sides));
    }
}

void ex_check_count(std::int64_t Layers, const std::string& rOp) {
    if (Layers < 1)
        throw std::invalid_argument("meshio++: " + rOp + ": num_layers must be at least 1");
}

}  // namespace

SweepResult extrude(const Mesh& rMesh, const ExtrudeOptions& rOptions) {
    const ExtrudeOptions& o = rOptions;
    const std::string op = "extrude";
    const ExVec d = o.mDirection;
    for (const double c : d)
        if (!std::isfinite(c))
            throw std::invalid_argument("meshio++: extrude: the direction must be finite");
    const double dl = ex_norm(d);
    if (!(dl > 0.0))
        throw std::invalid_argument("meshio++: extrude: the direction must not be zero");
    if (!std::isfinite(o.mRatio) || !(o.mRatio > 0.0))
        throw std::invalid_argument("meshio++: extrude: ratio must be positive and finite");
    std::vector<double> t;  // level fractions, 1..L
    if (!o.mOffsets.empty()) {
        if (o.mRatio != 1.0)
            throw std::invalid_argument(
                "meshio++: extrude: give either offsets or a ratio, not both");
        if (o.mNumLayers != 1 && static_cast<std::size_t>(o.mNumLayers) != o.mOffsets.size())
            throw std::invalid_argument("meshio++: extrude: num_layers (" +
                                        std::to_string(o.mNumLayers) + ") differs from the " +
                                        std::to_string(o.mOffsets.size()) + " offsets");
        double prev = 0.0;
        for (const double v : o.mOffsets) {
            if (!std::isfinite(v) || !(v > prev))
                throw std::invalid_argument(
                    "meshio++: extrude: offsets must be positive and increasing");
            prev = v;
        }
        if (prev != 1.0)
            throw std::invalid_argument("meshio++: extrude: the last offset must be 1");
        t = o.mOffsets;
    } else {
        ex_check_count(o.mNumLayers, op);
        const int L = o.mNumLayers;
        for (int k = 1; k <= L; ++k) {
            if (o.mRatio == 1.0)
                t.push_back(static_cast<double>(k) / static_cast<double>(L));
            else
                t.push_back((std::pow(o.mRatio, k) - 1.0) / (std::pow(o.mRatio, L) - 1.0));
        }
        t.back() = 1.0;
    }
    const int L = static_cast<int>(t.size());
    const bool uniform = o.mOffsets.empty() && o.mRatio == 1.0;

    SweepResult res;
    ExProfile pr = ex_gather(rMesh, op);
    if (pr.mCurve) {
        for (std::size_t l = 0; l < pr.mLines.size(); ++l) {
            const ExVec e = ex_sub(pr.Point(pr.mLines[l][1]), pr.Point(pr.mLines[l][0]));
            if (ex_norm(ex_cross(e, d)) <= 1e-9 * ex_norm(e) * dl)
                throw std::invalid_argument("meshio++: extrude: line " + std::to_string(l) +
                                            " is parallel to the direction (or has no length)");
        }
    } else {
        for (std::size_t f = 0; f < pr.mFaces.size(); ++f) {
            const ExVec nrm = ex_face_normal(pr, pr.mFaces[f]);
            const double s = ex_dot(nrm, d);
            if (std::abs(s) <= 1e-9 * ex_norm(nrm) * dl)
                throw std::invalid_argument("meshio++: extrude: face " + std::to_string(f) +
                                            " is parallel to the direction (or has no area)");
            if (s < 0.0) {
                ex_flip(pr.mFaces[f]);
                ++res.mNumFlipped;
            }
        }
    }
    const std::size_t n = rMesh.NumPoints();
    ExPlan plan;
    plan.mOp = op;
    plan.mLayers = L;
    plan.mOffsets.assign(n * static_cast<std::size_t>(L), 0.0);
    for (std::size_t p = 0; p < n; ++p)
        if (pr.mUsed[p])
            for (int k = 0; k < L; ++k)
                plan.mOffsets[p * static_cast<std::size_t>(L) + static_cast<std::size_t>(k)] =
                    t[static_cast<std::size_t>(k)];
    const std::vector<double>& xyz = pr.mXyz;
    // Uniform levels are computed as `x + (d k) / L`, the way `box` computes
    // its coordinates, so extruding `rectangle` gives `box` bit for bit.
    plan.mPlace = [&, uniform, L](std::size_t P, int K, double* pOut) {
        for (std::size_t c = 0; c < 3; ++c)
            pOut[c] = uniform ? xyz[3 * P + c] + d[c] * static_cast<double>(K) / static_cast<double>(L)
                              : xyz[3 * P + c] + d[c] * t[static_cast<std::size_t>(K - 1)];
    };
    plan.mLayerName = kExtrudeLayerName;
    plan.mRecordLayers = o.mRecordLayers;
    plan.mFirstName = "bottom";
    plan.mLastName = "top";
    plan.mMaxCells = o.mMaxCells;
    ex_sweep(rMesh, pr, plan, res);
    return res;
}

SweepResult revolve(const Mesh& rMesh, const RevolveOptions& rOptions) {
    const RevolveOptions& o = rOptions;
    const std::string op = "revolve";
    ExVec axis = o.mAxis;
    const ExVec origin = o.mAxisPoint;
    for (std::size_t c = 0; c < 3; ++c)
        if (!std::isfinite(axis[c]) || !std::isfinite(origin[c]))
            throw std::invalid_argument("meshio++: revolve: the axis must be finite");
    const double al = ex_norm(axis);
    if (!(al > 0.0))
        throw std::invalid_argument("meshio++: revolve: the axis must not be zero");
    if (!std::isfinite(o.mTolerance) || o.mTolerance < 0.0)
        throw std::invalid_argument("meshio++: revolve: tolerance must be finite and non-negative");
    if (!std::isfinite(o.mAngle) || o.mAngle == 0.0 || std::abs(o.mAngle) > 360.0)
        throw std::invalid_argument("meshio++: revolve: angle must be in (0, 360] degrees, either sign");
    ex_check_count(o.mNumLayers, op);
    const int L = o.mNumLayers;
    const bool full = std::abs(o.mAngle) == 360.0;
    if (full && L < 3)
        throw std::invalid_argument("meshio++: revolve: a full turn needs at least 3 steps");
    double angle = o.mAngle;
    if (angle < 0.0) {
        angle = -angle;
        for (double& c : axis)
            c = -c;
    }
    if (angle / L >= 180.0)
        throw std::invalid_argument("meshio++: revolve: each step must turn less than 180 degrees");
    for (double& c : axis)
        c /= al;

    SweepResult res;
    res.mFull = full;
    ExProfile pr = ex_gather(rMesh, op);
    const std::size_t n = rMesh.NumPoints();
    // On the axis: within the tolerance of the used points' bounding-box diagonal.
    ExVec lo{0.0, 0.0, 0.0}, hi{0.0, 0.0, 0.0};
    bool any = false;
    for (std::size_t p = 0; p < n; ++p) {
        if (!pr.mUsed[p])
            continue;
        const ExVec x = pr.Point(static_cast<std::int64_t>(p));
        for (std::size_t c = 0; c < 3; ++c) {
            lo[c] = any ? std::min(lo[c], x[c]) : x[c];
            hi[c] = any ? std::max(hi[c], x[c]) : x[c];
        }
        any = true;
    }
    const double tol = o.mTolerance * ex_norm(ex_sub(hi, lo));
    std::vector<char> on_axis(n, 0);
    const auto radial = [&](std::int64_t P) {
        const ExVec r = ex_sub(pr.Point(P), origin);
        return ex_cross(axis, r);  // the motion at P; its length is the distance to the axis
    };
    for (std::size_t p = 0; p < n; ++p)
        if (pr.mUsed[p] && ex_norm(radial(static_cast<std::int64_t>(p))) <= tol) {
            on_axis[p] = 1;
            ++res.mNumOnAxis;
        }
    if (pr.mCurve) {
        std::int64_t on = 0;
        for (std::size_t l = 0; l < pr.mLines.size(); ++l) {
            const auto& r_l = pr.mLines[l];
            if (on_axis[static_cast<std::size_t>(r_l[0])] && on_axis[static_cast<std::size_t>(r_l[1])]) {
                ++on;
                continue;
            }
            const ExVec a = pr.Point(r_l[0]), b = pr.Point(r_l[1]);
            const ExVec e = ex_sub(b, a);
            const ExVec mid{0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1]), 0.5 * (a[2] + b[2])};
            const ExVec m = ex_cross(axis, ex_sub(mid, origin));
            if (ex_norm(ex_cross(e, m)) <= 1e-9 * ex_norm(e) * ex_norm(m))
                throw std::invalid_argument("meshio++: revolve: line " + std::to_string(l) +
                                            " moves along itself (or has no length)");
        }
        if (on > 0)
            throw std::invalid_argument("meshio++: revolve: " + std::to_string(on) +
                                        " line(s) lie on the axis and would sweep to nothing");
    } else {
        std::int64_t kept = 0;
        for (std::size_t f = 0; f < pr.mFaces.size(); ++f) {
            const ExVec nrm = ex_face_normal(pr, pr.mFaces[f]);
            const int k = pr.mFaces[f][3] >= 0 ? 4 : 3;
            int pos = 0, neg = 0;
            for (int i = 0; i < k; ++i) {
                const std::int64_t v = pr.mFaces[f][static_cast<std::size_t>(i)];
                if (on_axis[static_cast<std::size_t>(v)])
                    continue;
                const ExVec m = radial(v);
                const double s = ex_dot(nrm, m);
                if (std::abs(s) <= 1e-9 * ex_norm(nrm) * ex_norm(m))
                    continue;
                (s > 0.0 ? pos : neg) += 1;
            }
            if (pos > 0 && neg > 0)
                throw std::invalid_argument("meshio++: revolve: face " + std::to_string(f) +
                                            " crosses the axis");
            if (pos == 0 && neg == 0)
                throw std::invalid_argument("meshio++: revolve: face " + std::to_string(f) +
                                            " turns within its own plane (or has no area)");
            if (neg > 0) {
                ex_flip(pr.mFaces[f]);
                ++res.mNumFlipped;
            } else {
                ++kept;
            }
        }
        // Faces turned both ways lie on both sides of the axis: past half a
        // turn, each side sweeps through the other.
        if (kept > 0 && res.mNumFlipped > 0) {
            if (angle > 180.0)
                throw std::invalid_argument(
                    "meshio++: revolve: the profile lies on both sides of the axis, so a turn "
                    "of more than 180 degrees would sweep it through itself");
            res.mWarnings.push_back("the profile lies on both sides of the axis");
        }
    }
    ExPlan plan;
    plan.mOp = op;
    plan.mLayers = L;
    plan.mPeriodic = full;
    plan.mCollapsedEdgeWedges = true;
    plan.mOffsets.assign(n * static_cast<std::size_t>(L), 0.0);
    for (std::size_t p = 0; p < n; ++p)
        if (pr.mUsed[p] && !on_axis[p])
            for (int k = 0; k < L; ++k)
                plan.mOffsets[p * static_cast<std::size_t>(L) + static_cast<std::size_t>(k)] = k + 1;
    const double rad = angle * 3.14159265358979323846 / 180.0;
    std::vector<std::array<double, 16>> rot;
    for (int k = 1; k <= L; ++k)
        rot.push_back(transform_rotation(axis[0], axis[1], axis[2], rad * k / L).mMatrix);
    const std::vector<double>& xyz = pr.mXyz;
    plan.mPlace = [&](std::size_t P, int K, double* pOut) {
        const auto& m = rot[static_cast<std::size_t>(K - 1)];
        const double r[3] = {xyz[3 * P] - origin[0], xyz[3 * P + 1] - origin[1],
                             xyz[3 * P + 2] - origin[2]};
        for (std::size_t c = 0; c < 3; ++c)
            pOut[c] = origin[c] + m[c * 4] * r[0] + m[c * 4 + 1] * r[1] + m[c * 4 + 2] * r[2];
    };
    plan.mLayerName = kRevolveLayerName;
    plan.mRecordLayers = o.mRecordLayers;
    plan.mFirstName = "start";
    plan.mLastName = "end";
    plan.mMaxCells = o.mMaxCells;
    ex_sweep(rMesh, pr, plan, res);
    return res;
}

}  // namespace meshioplusplus

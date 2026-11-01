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
// select_region: write one named region of points, cells or sides chosen by
// geometry, a data or quality threshold, a component, or a patch grown from a
// seed. Each criterion fills a per-candidate pass mask (in parallel, disjoint
// slots); the masks are intersected and the survivors written as the region.
// The growth is a serial breadth-first walk over the facets, crossing only the
// edges the shared crease test (detail/crease_edges) calls smooth. See
// operations/select.hpp for the contract.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/select.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/operations/data_common.hpp"
#include "meshioplusplus/operations/quality.hpp"
#include "meshioplusplus/detail/cell_edges.hpp"
#include "meshioplusplus/detail/cell_faces.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/parallel.hpp"

// Project includes (private, not installed)
#include "../detail/cell_components.hpp"
#include "../detail/crease_edges.hpp"

namespace meshioplusplus {

namespace {

using detail::Vec3;

constexpr double kSelPi = 3.14159265358979323846;

// The facets a Side selection (or a growth) runs over: CSR node lists, corners
// first, with each facet's owner cell, local facet number and corner count.
struct SelFacets {
    std::vector<std::int64_t> mOwner;
    std::vector<std::int64_t> mLocal;
    std::vector<std::int64_t> mStart{0};
    std::vector<std::int64_t> mNodes;
    std::vector<std::int32_t> mCorners;

    std::size_t Size() const { return mOwner.size(); }
    void Add(std::int64_t Owner, std::int64_t Local, const std::int64_t* pNodes, std::size_t N,
             std::size_t Corners) {
        mOwner.push_back(Owner);
        mLocal.push_back(Local);
        mNodes.insert(mNodes.end(), pNodes, pNodes + N);
        mStart.push_back(static_cast<std::int64_t>(mNodes.size()));
        mCorners.push_back(static_cast<std::int32_t>(Corners));
    }
};

// A block's topological dimension; polygons are 2, polyhedra 3.
int sel_block_dim(const Mesh::CellView& rCb) {
    if (rCb.IsPolyhedron())
        return 3;
    return cell_type_dimension(cell_type_from_name(std::string(rCb.Type())));
}

// Every node of one cell, as stored (a polyhedron's nodes once per face).
void sel_cell_nodes(const Mesh::CellView& rCb, const std::int64_t* pConn, std::size_t c,
                    std::vector<std::int64_t>& rOut) {
    rOut.clear();
    if (rCb.IsPolyhedron()) {
        for (std::size_t f = 0; f < rCb.NumFaces(c); ++f) {
            const auto face = rCb.Face(c, f);
            rOut.insert(rOut.end(), face.first, face.first + face.second);
        }
    } else if (rCb.IsRagged()) {
        const std::int64_t* row = rCb.Row(c);
        rOut.assign(row, row + rCb.RowSize(c));
    } else {
        const std::size_t npc = rCb.NodesPerCell();
        rOut.assign(pConn + c * npc, pConn + (c + 1) * npc);
    }
}

// A rectangular block's connectivity as Int64 (a copy only when it is not).
NDArray sel_int64_conn(const Mesh::CellView& rCb) {
    if (rCb.IsPolyhedron() || rCb.IsRagged())
        return NDArray();
    const NDArray& conn = rCb.Conn();
    if (conn.Dtype() == DType::Int64)
        return conn;
    NDArray out = NDArray::Uninit(DType::Int64, {conn.Size()});
    std::int64_t* dst = out.As<std::int64_t>();
    for (std::size_t i = 0; i < conn.Size(); ++i)
        dst[i] = detail::read_int(conn, i);
    return out;
}

// The facets of every cell of dimension `Dim` (faces for 3, edges for 2).
SelFacets sel_collect_facets(const Mesh& rMesh, int Dim) {
    SelFacets out;
    std::int64_t base = 0;
    std::vector<std::int64_t> nodes;
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto cb = rMesh.Cells(b);
        const std::size_t nc = cb.NumCells();
        if (sel_block_dim(cb) != Dim) {
            base += static_cast<std::int64_t>(nc);
            continue;
        }
        const NDArray conn = sel_int64_conn(cb);
        const std::int64_t* pc = conn.Size() ? conn.As<std::int64_t>() : nullptr;
        const CellType type =
            cb.IsPolyhedron() ? CellType::Polyhedron : cell_type_from_name(std::string(cb.Type()));
        for (std::size_t c = 0; c < nc; ++c) {
            const std::int64_t cell = base + static_cast<std::int64_t>(c);
            if (cb.IsPolyhedron()) {
                for (std::size_t f = 0; f < cb.NumFaces(c); ++f) {
                    const auto face = cb.Face(c, f);
                    out.Add(cell, static_cast<std::int64_t>(f), face.first, face.second,
                            face.second);
                }
                continue;
            }
            if (cb.IsRagged() || type == CellType::Polygon) {
                sel_cell_nodes(cb, pc, c, nodes);
                const std::size_t n = nodes.size();
                for (std::size_t k = 0; k < n; ++k) {
                    const std::int64_t e[2] = {nodes[k], nodes[(k + 1) % n]};
                    out.Add(cell, static_cast<std::int64_t>(k), e, 2, 2);
                }
                continue;
            }
            const std::size_t npc = cb.NodesPerCell();
            const std::int64_t* row = pc + c * npc;
            std::int64_t buf[9];
            if (Dim == 3) {
                const auto& defs = detail::cell_faces(type);
                for (std::size_t f = 0; f < defs.size(); ++f) {
                    for (std::size_t k = 0; k < defs[f].mNumNodes; ++k)
                        buf[k] = row[defs[f].mNodes[k]];
                    out.Add(cell, static_cast<std::int64_t>(f), buf, defs[f].mNumNodes,
                            defs[f].mNumCorners);
                }
            } else {
                const auto& defs = detail::cell_edges(type);
                for (std::size_t f = 0; f < defs.size(); ++f) {
                    for (std::size_t k = 0; k < defs[f].mNumNodes; ++k)
                        buf[k] = row[defs[f].mNodes[k]];
                    out.Add(cell, static_cast<std::int64_t>(f), buf, defs[f].mNumNodes,
                            defs[f].mNumCorners);
                }
            }
        }
        base += static_cast<std::int64_t>(nc);
    }
    return out;
}

// How many facets share each facet's corner set (1 on the boundary).
std::vector<std::int64_t> sel_facet_multiplicity(const SelFacets& rF) {
    const std::size_t n = rF.Size();
    std::vector<std::int64_t> keys(rF.mNodes.size());
    parallel_for(n, [&](std::size_t f) {
        const auto s = static_cast<std::size_t>(rF.mStart[f]);
        const auto nc = static_cast<std::size_t>(rF.mCorners[f]);
        std::copy(rF.mNodes.begin() + static_cast<std::ptrdiff_t>(s),
                  rF.mNodes.begin() + static_cast<std::ptrdiff_t>(s + nc),
                  keys.begin() + static_cast<std::ptrdiff_t>(s));
        std::sort(keys.begin() + static_cast<std::ptrdiff_t>(s),
                  keys.begin() + static_cast<std::ptrdiff_t>(s + nc));
    });
    const auto key_less = [&](std::size_t a, std::size_t b) {
        const auto na = static_cast<std::size_t>(rF.mCorners[a]);
        const auto nb = static_cast<std::size_t>(rF.mCorners[b]);
        if (na != nb)
            return na < nb;
        const std::int64_t* ka = keys.data() + rF.mStart[a];
        const std::int64_t* kb = keys.data() + rF.mStart[b];
        for (std::size_t k = 0; k < na; ++k)
            if (ka[k] != kb[k])
                return ka[k] < kb[k];
        return a < b;
    };
    const auto key_equal = [&](std::size_t a, std::size_t b) {
        const auto na = static_cast<std::size_t>(rF.mCorners[a]);
        if (na != static_cast<std::size_t>(rF.mCorners[b]))
            return false;
        return std::equal(keys.data() + rF.mStart[a], keys.data() + rF.mStart[a] + na,
                          keys.data() + rF.mStart[b]);
    };
    std::vector<std::size_t> order(n);
    for (std::size_t i = 0; i < n; ++i)
        order[i] = i;
    parallel_sort(order.begin(), order.end(), key_less);
    std::vector<std::int64_t> mult(n, 0);
    for (std::size_t i = 0; i < n;) {
        std::size_t j = i + 1;
        while (j < n && key_equal(order[i], order[j]))
            ++j;
        for (std::size_t k = i; k < j; ++k)
            mult[order[k]] = static_cast<std::int64_t>(j - i);
        i = j;
    }
    return mult;
}

// Keep only the facets whose flag is set, in order.
SelFacets sel_filter_facets(const SelFacets& rF, const std::vector<char>& rKeep) {
    SelFacets out;
    for (std::size_t f = 0; f < rF.Size(); ++f) {
        if (!rKeep[f])
            continue;
        const auto s = static_cast<std::size_t>(rF.mStart[f]);
        const auto e = static_cast<std::size_t>(rF.mStart[f + 1]);
        out.Add(rF.mOwner[f], rF.mLocal[f], rF.mNodes.data() + s, e - s,
                static_cast<std::size_t>(rF.mCorners[f]));
    }
    return out;
}

// The facets of the 2-D cells, as faces to grow over (rings of corners).
SelFacets sel_surface_cells(const Mesh& rMesh) {
    SelFacets out;
    std::int64_t base = 0;
    std::vector<std::int64_t> nodes;
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto cb = rMesh.Cells(b);
        const std::size_t nc = cb.NumCells();
        if (sel_block_dim(cb) == 2) {
            const NDArray conn = sel_int64_conn(cb);
            const std::int64_t* pc = conn.Size() ? conn.As<std::int64_t>() : nullptr;
            const CellType type = cell_type_from_name(std::string(cb.Type()));
            const bool polygon = cb.IsRagged() || type == CellType::Polygon;
            const std::size_t corners = polygon ? 0 : detail::cell_edges(type).size();
            for (std::size_t c = 0; c < nc; ++c) {
                sel_cell_nodes(cb, pc, c, nodes);
                const std::size_t nk = polygon ? nodes.size() : corners;
                if (nk < 3)
                    continue;  // a 2-D type with no edge table has no ring
                out.Add(base + static_cast<std::int64_t>(c), -1, nodes.data(), nodes.size(), nk);
            }
        }
        base += static_cast<std::int64_t>(nc);
    }
    return out;
}

// The grown patch: a flag per facet of rF.
std::vector<char> sel_grow(const SelFacets& rF, const std::vector<double>& rXyz,
                           const SelectOptions& rOpt, std::int64_t& rNumGrown) {
    const std::size_t n = rF.Size();
    if (n == 0)
        throw std::invalid_argument("select: a seed needs faces to grow over, and there are none");
    // Rings of corners, their normals and centroids.
    std::vector<std::int64_t> start(n + 1, 0), ring;
    for (std::size_t f = 0; f < n; ++f) {
        const auto s = static_cast<std::size_t>(rF.mStart[f]);
        ring.insert(ring.end(), rF.mNodes.begin() + static_cast<std::ptrdiff_t>(s),
                    rF.mNodes.begin() + static_cast<std::ptrdiff_t>(s) + rF.mCorners[f]);
        start[f + 1] = static_cast<std::int64_t>(ring.size());
    }
    std::vector<double> normal(3 * n, 0.0), centroid(3 * n, 0.0);
    parallel_for(n, [&](std::size_t f) {
        const auto s = static_cast<std::size_t>(start[f]);
        const auto k = static_cast<std::size_t>(start[f + 1] - start[f]);
        detail::ring_unit_normal(rXyz.data(), ring.data() + s, k, normal.data() + 3 * f);
        for (std::size_t i = 0; i < k; ++i)
            for (int d = 0; d < 3; ++d)
                centroid[3 * f + d] += rXyz[3 * static_cast<std::size_t>(ring[s + i]) + d];
        for (int d = 0; d < 3; ++d)
            centroid[3 * f + d] /= static_cast<double>(k);
    });

    // The seed: the facet whose centroid is nearest (ties: the lowest index).
    std::size_t seed = 0;
    double best = std::numeric_limits<double>::infinity();
    for (std::size_t f = 0; f < n; ++f) {
        double d2 = 0.0;
        for (int d = 0; d < 3; ++d) {
            const double t = centroid[3 * f + d] - rOpt.mSeedPoint[static_cast<std::size_t>(d)];
            d2 += t * t;
        }
        if (d2 < best) {
            best = d2;
            seed = f;
        }
    }

    // Smooth edges: used by exactly two faces and not sharp.
    const std::vector<detail::CreaseEdge> creases =
        detail::crease_edges(start, ring, normal, rOpt.mFeatureAngleDeg);
    struct Use {
        std::int64_t mLo, mHi, mFace;
    };
    std::vector<Use> uses;
    uses.reserve(ring.size());
    for (std::size_t f = 0; f < n; ++f) {
        const auto s = static_cast<std::size_t>(start[f]);
        const auto k = static_cast<std::size_t>(start[f + 1] - start[f]);
        for (std::size_t i = 0; i < k; ++i) {
            std::int64_t a = ring[s + i], b = ring[s + (i + 1) % k];
            if (a == b)
                continue;
            if (a > b)
                std::swap(a, b);
            uses.push_back({a, b, static_cast<std::int64_t>(f)});
        }
    }
    parallel_sort(uses.begin(), uses.end(), [](const Use& rA, const Use& rB) {
        if (rA.mLo != rB.mLo)
            return rA.mLo < rB.mLo;
        if (rA.mHi != rB.mHi)
            return rA.mHi < rB.mHi;
        return rA.mFace < rB.mFace;
    });
    std::vector<std::vector<std::int64_t>> adj(n);
    for (std::size_t i = 0; i < uses.size();) {
        std::size_t j = i + 1;
        while (j < uses.size() && uses[j].mLo == uses[i].mLo && uses[j].mHi == uses[i].mHi)
            ++j;
        if (j - i == 2 && uses[i].mFace != uses[i + 1].mFace) {
            const auto it =
                std::lower_bound(creases.begin(), creases.end(), uses[i],
                                 [](const detail::CreaseEdge& rE, const Use& rU) {
                                     return rE.mLo != rU.mLo ? rE.mLo < rU.mLo : rE.mHi < rU.mHi;
                                 });
            const bool listed =
                it != creases.end() && it->mLo == uses[i].mLo && it->mHi == uses[i].mHi;
            if (!listed || !it->IsCrease()) {
                adj[static_cast<std::size_t>(uses[i].mFace)].push_back(uses[i + 1].mFace);
                adj[static_cast<std::size_t>(uses[i + 1].mFace)].push_back(uses[i].mFace);
            }
        }
        i = j;
    }
    parallel_for(n, [&](std::size_t f) { std::sort(adj[f].begin(), adj[f].end()); });

    const bool normal_test = rOpt.mNormalToleranceDeg > 0.0;
    const double cos_tol = std::cos(rOpt.mNormalToleranceDeg * kSelPi / 180.0);
    const double* ns = normal.data() + 3 * seed;
    std::vector<char> grown(n, 0);
    std::deque<std::size_t> queue{seed};
    grown[seed] = 1;
    rNumGrown = 1;
    while (!queue.empty()) {
        const std::size_t f = queue.front();
        queue.pop_front();
        for (const std::int64_t g64 : adj[f]) {
            const auto g = static_cast<std::size_t>(g64);
            if (grown[g])
                continue;
            if (normal_test) {
                const double* ng = normal.data() + 3 * g;
                if (ns[0] * ng[0] + ns[1] * ng[1] + ns[2] * ng[2] < cos_tol)
                    continue;
            }
            grown[g] = 1;
            ++rNumGrown;
            queue.push_back(g);
        }
    }
    return grown;
}

// A scalar per global cell from a cell_data array covering every block.
std::vector<double> sel_cell_values(const Mesh& rMesh, const std::string& rArray) {
    if (!rMesh.HasCellData(rArray)) {
        if (rMesh.HasPointData(rArray))
            throw std::invalid_argument(
                "select: '" + rArray +
                "' is point_data; a cell or side selection compares cell_data (convert it "
                "first with `data to-cell`, or select points)");
        throw std::invalid_argument("select: " +
                                    data_unknown_key_message(rMesh, DataLocation::Cell, rArray));
    }
    const std::size_t nblocks = rMesh.NumCellBlocks();
    if (rMesh.CellDataNumBlocks(rArray) != nblocks)
        throw std::invalid_argument("select: cell_data '" + rArray +
                                    "' does not cover every cell block");
    std::vector<double> out;
    for (std::size_t b = 0; b < nblocks; ++b) {
        const NDArray& a = rMesh.CellData(rArray, b);
        const std::size_t nc = rMesh.Cells(b).NumCells();
        if (detail::rows(a) != nc)
            throw std::invalid_argument("select: cell_data '" + rArray + "' has " +
                                        std::to_string(detail::rows(a)) + " rows on a block of " +
                                        std::to_string(nc) + " cells");
        if (nc != 0 && data_num_components(a) != 1)
            throw std::invalid_argument("select: cell_data '" + rArray +
                                        "' must be scalar (one value per cell)");
        const std::size_t off = out.size();
        out.resize(off + nc);
        parallel_for(nc, [&](std::size_t c) { out[off + c] = detail::read_double(a, c); });
    }
    return out;
}

void sel_check_angle(double Value, const char* pName) {
    if (!(Value >= 0.0 && Value <= 180.0))
        throw std::invalid_argument(std::string("select: ") + pName + " must be in [0, 180]");
}

}  // namespace

SelectMode select_mode_from_name(const std::string& rName) {
    if (rName == "all")
        return SelectMode::All;
    if (rName == "any")
        return SelectMode::Any;
    if (rName == "centroid")
        return SelectMode::Centroid;
    throw std::invalid_argument("select: unknown mode '" + rName +
                                "' (expected 'all', 'any' or 'centroid')");
}

const char* select_mode_name(SelectMode Mode) {
    switch (Mode) {
        case SelectMode::Any:
            return "any";
        case SelectMode::Centroid:
            return "centroid";
        case SelectMode::All:
        default:
            return "all";
    }
}

SelectResult select_region(const Mesh& rMesh, const SelectOptions& rOpt) {
    if (rOpt.mName.empty())
        throw std::invalid_argument("select: the region name must not be empty");
    const bool geometric = rOpt.mBox || rOpt.mPlane || rOpt.mSphere;
    const bool where = !rOpt.mWhereArray.empty();
    const bool quality = !rOpt.mQualityMetric.empty();
    const bool component = rOpt.mComponent >= 0;
    if (!geometric && !where && !quality && !component && !rOpt.mSeed)
        throw std::invalid_argument(
            "select: no criterion given (box, plane, sphere, where, quality, component or seed)");
    if (rOpt.mMode != SelectMode::All && rOpt.mMode != SelectMode::Any &&
        rOpt.mMode != SelectMode::Centroid)
        throw std::invalid_argument("select: unknown mode");
    if (rOpt.mKind != RegionKind::Point && rOpt.mKind != RegionKind::Cell &&
        rOpt.mKind != RegionKind::Side)
        throw std::invalid_argument("select: unknown region kind");
    if (rOpt.mSphere && !(rOpt.mSphereRadius >= 0.0))
        throw std::invalid_argument("select: the sphere radius must be non-negative");
    if (rOpt.mPlane && rOpt.mPlaneNormal[0] == 0.0 && rOpt.mPlaneNormal[1] == 0.0 &&
        rOpt.mPlaneNormal[2] == 0.0)
        throw std::invalid_argument("select: the plane normal must not be zero");
    if (quality && rOpt.mKind == RegionKind::Point)
        throw std::invalid_argument(
            "select: a quality criterion scores cells; select cells or sides, not points");
    if (rOpt.mSeed) {
        sel_check_angle(rOpt.mFeatureAngleDeg, "feature_angle");
        sel_check_angle(rOpt.mNormalToleranceDeg, "normal_tolerance");
    }

    const std::size_t npoints = rMesh.NumPoints();
    const std::size_t pdim = rMesh.PointDim();
    std::vector<double> xyz(3 * npoints, 0.0);
    {
        const NDArray& pts = rMesh.Points();
        parallel_for(npoints, [&](std::size_t i) {
            const Vec3 p = detail::read_point(pts, pdim, static_cast<std::int64_t>(i));
            xyz[3 * i] = p[0];
            xyz[3 * i + 1] = p[1];
            xyz[3 * i + 2] = p[2];
        });
    }
    const auto point_at = [&](std::int64_t v) {
        const auto i = static_cast<std::size_t>(v);
        return Vec3{xyz[3 * i], xyz[3 * i + 1], xyz[3 * i + 2]};
    };
    const auto inside = [&](const Vec3& p) {
        if (rOpt.mBox)
            for (std::size_t d = 0; d < 3; ++d)
                if (p[d] < rOpt.mBoxLo[d] || p[d] > rOpt.mBoxHi[d])
                    return false;
        if (rOpt.mPlane) {
            double s = 0.0;
            for (std::size_t d = 0; d < 3; ++d)
                s += (p[d] - rOpt.mPlanePoint[d]) * rOpt.mPlaneNormal[d];
            if (s < 0.0)
                return false;
        }
        if (rOpt.mSphere) {
            double s = 0.0;
            for (std::size_t d = 0; d < 3; ++d) {
                const double t = p[d] - rOpt.mSphereCenter[d];
                s += t * t;
            }
            if (s > rOpt.mSphereRadius * rOpt.mSphereRadius)
                return false;
        }
        return true;
    };
    // A node list's verdict under the mode.
    const auto nodes_inside = [&](const std::int64_t* pNodes, std::size_t N,
                                  std::vector<std::int64_t>& rScratch) {
        if (N == 0)
            return false;
        if (rOpt.mMode == SelectMode::Centroid) {
            rScratch.assign(pNodes, pNodes + N);
            std::sort(rScratch.begin(), rScratch.end());
            rScratch.erase(std::unique(rScratch.begin(), rScratch.end()), rScratch.end());
            Vec3 c{0.0, 0.0, 0.0};
            for (const std::int64_t v : rScratch) {
                const Vec3 p = point_at(v);
                for (std::size_t d = 0; d < 3; ++d)
                    c[d] += p[d];
            }
            for (std::size_t d = 0; d < 3; ++d)
                c[d] /= static_cast<double>(rScratch.size());
            return inside(c);
        }
        const bool all = rOpt.mMode == SelectMode::All;
        for (std::size_t k = 0; k < N; ++k) {
            const bool in = inside(point_at(pNodes[k]));
            if (all && !in)
                return false;
            if (!all && in)
                return true;
        }
        return all;
    };

    // Global cell numbering and dimensions.
    const std::size_t nblocks = rMesh.NumCellBlocks();
    std::vector<std::int64_t> base(nblocks + 1, 0);
    int top_dim = -1;
    bool has_surface = false;
    for (std::size_t b = 0; b < nblocks; ++b) {
        const auto cb = rMesh.Cells(b);
        base[b + 1] = base[b] + static_cast<std::int64_t>(cb.NumCells());
        const int d = sel_block_dim(cb);
        if (cb.NumCells() > 0) {
            top_dim = std::max(top_dim, d);
            has_surface = has_surface || d == 2;
        }
    }
    const auto ncells = static_cast<std::size_t>(base[nblocks]);

    SelectResult result;
    // The candidates, and for sides their facets.
    SelFacets sides;
    std::size_t ncand = 0;
    if (rOpt.mKind == RegionKind::Point) {
        ncand = npoints;
    } else if (rOpt.mKind == RegionKind::Cell) {
        ncand = ncells;
    } else {
        if (top_dim < 2)
            throw std::invalid_argument(
                "select: side regions name faces of 3-D cells or edges of 2-D cells, and this "
                "mesh has neither");
        sides = sel_collect_facets(rMesh, top_dim);
        if (rOpt.mBoundaryOnly) {
            const std::vector<std::int64_t> mult = sel_facet_multiplicity(sides);
            std::vector<char> keep(sides.Size());
            for (std::size_t f = 0; f < sides.Size(); ++f)
                keep[f] = mult[f] == 1 ? 1 : 0;
            sides = sel_filter_facets(sides, keep);
        }
        ncand = sides.Size();
    }
    result.mNumCandidates = static_cast<std::int64_t>(ncand);
    std::vector<char> pass(ncand, 1);

    // The owner cell of a candidate, for the per-cell criteria.
    const auto owner_of = [&](std::size_t i) -> std::size_t {
        return rOpt.mKind == RegionKind::Side ? static_cast<std::size_t>(sides.mOwner[i]) : i;
    };

    // Geometry.
    if (geometric) {
        if (rOpt.mKind == RegionKind::Point) {
            parallel_for(ncand, [&](std::size_t i) {
                if (!inside(point_at(static_cast<std::int64_t>(i))))
                    pass[i] = 0;
            });
        } else if (rOpt.mKind == RegionKind::Side) {
            parallel_for(ncand, [&](std::size_t i) {
                std::vector<std::int64_t> scratch;
                const auto s = static_cast<std::size_t>(sides.mStart[i]);
                const auto e = static_cast<std::size_t>(sides.mStart[i + 1]);
                if (!nodes_inside(sides.mNodes.data() + s, e - s, scratch))
                    pass[i] = 0;
            });
        } else {
            for (std::size_t b = 0; b < nblocks; ++b) {
                const auto cb = rMesh.Cells(b);
                const NDArray conn = sel_int64_conn(cb);
                const std::int64_t* pc = conn.Size() ? conn.As<std::int64_t>() : nullptr;
                const auto off = static_cast<std::size_t>(base[b]);
                parallel_for(cb.NumCells(), [&](std::size_t c) {
                    std::vector<std::int64_t> nodes, scratch;
                    sel_cell_nodes(cb, pc, c, nodes);
                    if (!nodes_inside(nodes.data(), nodes.size(), scratch))
                        pass[off + c] = 0;
                });
            }
        }
    }

    // Data.
    if (where) {
        if (rOpt.mKind == RegionKind::Point) {
            if (!rMesh.HasPointData(rOpt.mWhereArray)) {
                if (rMesh.HasCellData(rOpt.mWhereArray))
                    throw std::invalid_argument(
                        "select: '" + rOpt.mWhereArray +
                        "' is cell_data; a point selection compares point_data (convert it first "
                        "with `data to-point`, or select cells)");
                throw std::invalid_argument(
                    "select: " +
                    data_unknown_key_message(rMesh, DataLocation::Point, rOpt.mWhereArray));
            }
            const NDArray& a = rMesh.PointData(rOpt.mWhereArray);
            if (detail::rows(a) != npoints || (npoints != 0 && data_num_components(a) != 1))
                throw std::invalid_argument("select: point_data '" + rOpt.mWhereArray +
                                            "' must be scalar (one value per point)");
            parallel_for(ncand, [&](std::size_t i) {
                if (!refine_compare_value(detail::read_double(a, i), rOpt.mWhereOp,
                                          rOpt.mWhereValue))
                    pass[i] = 0;
            });
        } else {
            const std::vector<double> v = sel_cell_values(rMesh, rOpt.mWhereArray);
            parallel_for(ncand, [&](std::size_t i) {
                if (!refine_compare_value(v[owner_of(i)], rOpt.mWhereOp, rOpt.mWhereValue))
                    pass[i] = 0;
            });
        }
    }

    // Quality.
    if (quality) {
        const QualityReport report = compute_quality(rMesh);
        // The metric as compute_quality names it, with or without the
        // `quality:` prefix attach_quality's arrays carry.
        static const std::string kPrefix = "quality:";
        const std::string& want = rOpt.mQualityMetric;
        const std::string bare =
            want.compare(0, kPrefix.size(), kPrefix) == 0 ? want.substr(kPrefix.size()) : want;
        const std::vector<NDArray>* arrays = nullptr;
        std::string known;
        for (const auto& entry : report.mCellArrays) {
            const std::string name = entry.first.compare(0, kPrefix.size(), kPrefix) == 0
                                         ? entry.first.substr(kPrefix.size())
                                         : entry.first;
            if (name == bare)
                arrays = &entry.second;
            known += (known.empty() ? "" : ", ") + name;
        }
        if (arrays == nullptr)
            throw std::invalid_argument("select: unknown quality metric '" + rOpt.mQualityMetric +
                                        "' (known: " + known + ")");
        std::vector<double> v(ncells, std::numeric_limits<double>::quiet_NaN());
        for (std::size_t b = 0; b < nblocks && b < arrays->size(); ++b) {
            const NDArray& a = (*arrays)[b];
            const auto off = static_cast<std::size_t>(base[b]);
            const std::size_t nc = static_cast<std::size_t>(base[b + 1] - base[b]);
            if (detail::rows(a) != nc)
                continue;
            parallel_for(nc, [&](std::size_t c) { v[off + c] = detail::read_double(a, c); });
        }
        parallel_for(ncand, [&](std::size_t i) {
            if (!refine_compare_value(v[owner_of(i)], rOpt.mQualityOp, rOpt.mQualityValue))
                pass[i] = 0;
        });
    }

    // Component.
    if (component) {
        const detail::CellComponents comps = detail::cell_components(rMesh);
        result.mNumComponents = comps.mCount;
        if (rOpt.mKind == RegionKind::Point) {
            std::vector<char> used(npoints, 0);
            std::vector<std::int64_t> nodes;
            for (std::size_t b = 0; b < nblocks; ++b) {
                const auto cb = rMesh.Cells(b);
                const NDArray conn = sel_int64_conn(cb);
                const std::int64_t* pc = conn.Size() ? conn.As<std::int64_t>() : nullptr;
                for (std::size_t c = 0; c < cb.NumCells(); ++c) {
                    if (comps.mLabel[static_cast<std::size_t>(base[b]) + c] != rOpt.mComponent)
                        continue;
                    sel_cell_nodes(cb, pc, c, nodes);
                    for (const std::int64_t v : nodes)
                        used[static_cast<std::size_t>(v)] = 1;
                }
            }
            parallel_for(ncand, [&](std::size_t i) {
                if (!used[i])
                    pass[i] = 0;
            });
        } else {
            parallel_for(ncand, [&](std::size_t i) {
                if (comps.mLabel[owner_of(i)] != rOpt.mComponent)
                    pass[i] = 0;
            });
        }
    }

    // Growth.
    if (rOpt.mSeed) {
        if (rOpt.mKind == RegionKind::Cell) {
            const SelFacets faces = sel_surface_cells(rMesh);
            if (faces.Size() == 0)
                throw std::invalid_argument(
                    "select: growing a cell region needs surface (2-D) cells; for the boundary "
                    "of a volume mesh select sides instead");
            const std::vector<char> grown = sel_grow(faces, xyz, rOpt, result.mNumGrown);
            std::vector<char> in_patch(ncells, 0);
            for (std::size_t f = 0; f < faces.Size(); ++f)
                if (grown[f])
                    in_patch[static_cast<std::size_t>(faces.mOwner[f])] = 1;
            parallel_for(ncand, [&](std::size_t i) {
                if (!in_patch[i])
                    pass[i] = 0;
            });
        } else if (rOpt.mKind == RegionKind::Side) {
            if (top_dim != 3)
                throw std::invalid_argument(
                    "select: growing needs faces; the sides of this mesh are edges (select "
                    "cells of the surface instead)");
            const std::vector<char> grown = sel_grow(sides, xyz, rOpt, result.mNumGrown);
            parallel_for(ncand, [&](std::size_t i) {
                if (!grown[i])
                    pass[i] = 0;
            });
        } else {
            SelFacets faces;
            if (has_surface) {
                faces = sel_surface_cells(rMesh);
            } else if (top_dim == 3) {
                SelFacets all = sel_collect_facets(rMesh, 3);
                const std::vector<std::int64_t> mult = sel_facet_multiplicity(all);
                std::vector<char> keep(all.Size());
                for (std::size_t f = 0; f < all.Size(); ++f)
                    keep[f] = mult[f] == 1 ? 1 : 0;
                faces = sel_filter_facets(all, keep);
            }
            if (faces.Size() == 0)
                throw std::invalid_argument(
                    "select: growing needs surface cells or the faces of volume cells, and this "
                    "mesh has neither");
            const std::vector<char> grown = sel_grow(faces, xyz, rOpt, result.mNumGrown);
            std::vector<char> used(npoints, 0);
            for (std::size_t f = 0; f < faces.Size(); ++f) {
                if (!grown[f])
                    continue;
                for (auto k = faces.mStart[f]; k < faces.mStart[f + 1]; ++k)
                    used[static_cast<std::size_t>(faces.mNodes[static_cast<std::size_t>(k)])] = 1;
            }
            parallel_for(ncand, [&](std::size_t i) {
                if (!used[i])
                    pass[i] = 0;
            });
        }
    }

    // The region.
    std::size_t nsel = 0;
    for (std::size_t i = 0; i < ncand; ++i)
        nsel += pass[i] ? 1 : 0;
    result.mNumSelected = static_cast<std::int64_t>(nsel);
    const bool side = rOpt.mKind == RegionKind::Side;
    NDArray entries = side ? NDArray::Uninit(DType::Int64, {nsel, std::size_t{2}})
                           : NDArray::Uninit(DType::Int64, {nsel});
    if (nsel > 0) {
        std::int64_t* dst = entries.As<std::int64_t>();
        std::size_t j = 0;
        for (std::size_t i = 0; i < ncand; ++i) {
            if (!pass[i])
                continue;
            if (side) {
                dst[2 * j] = sides.mOwner[i];
                dst[2 * j + 1] = sides.mLocal[i];
            } else {
                dst[j] = static_cast<std::int64_t>(i);
            }
            ++j;
        }
    }
    meshioplusplus::Region region(rOpt.mName, rOpt.mKind, side ? top_dim - 1 : -1, rOpt.mTag,
                                  std::move(entries));
    region.Canonicalize();

    Mesh& out = result.mMesh;
    out = detail::clone_mesh(rMesh,
                             [](DataLocation, const std::string&, std::string&) { return true; });
    for (std::size_t i = out.NumRegions(); i-- > 0;) {
        const meshioplusplus::Region& r = out.Region(i);
        if (r.mKind == rOpt.mKind && r.mName == rOpt.mName)
            out.RemoveRegion(i);
    }
    out.AddRegion(region);
    result.mRegion = std::move(region);
    return result;
}

}  // namespace meshioplusplus

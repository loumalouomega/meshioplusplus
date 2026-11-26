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
//                   SPDX-License-Identifier: MIT
//                   meshio++ default license: LICENSE
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//
// check_conformity. See operations/check_conformity.hpp for the contract.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/check_conformity.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/cell_edges.hpp"
#include "meshioplusplus/detail/cell_faces.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/face_mesh.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/point_triangle.hpp"
#include "meshioplusplus/detail/bucket_table.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/parallel.hpp"

// Project includes (private, not installed)
#include "../detail/bucket_build.hpp"
#include "../detail/surface_rings.hpp"
#include "../detail/typed_view.hpp"

namespace meshioplusplus {
namespace {

constexpr const char* kCcPrefix = "meshio++: check_conformity: ";

// The unmatched facets of the examined cells, as corner rings (CSR), each with
// the global cell that owns it.
struct CcFacets {
    std::vector<std::int64_t> mStart{0};
    std::vector<std::int64_t> mNodes;
    std::vector<std::int64_t> mOwner;
    std::int64_t mNumNonManifold = 0;

    std::size_t Size() const { return mOwner.size(); }
    std::size_t Corners(std::size_t Facet) const {
        return static_cast<std::size_t>(mStart[Facet + 1] - mStart[Facet]);
    }
    const std::int64_t* Ring(std::size_t Facet) const { return mNodes.data() + mStart[Facet]; }
    void Push(const std::int64_t* pIds, std::size_t Size, std::int64_t Owner) {
        mNodes.insert(mNodes.end(), pIds, pIds + Size);
        mStart.push_back(static_cast<std::int64_t>(mNodes.size()));
        mOwner.push_back(Owner);
    }
};

// The topological dimension of a block: 3, 2, or -1 for anything else.
template <class TView>
int cc_block_dimension(const TView& rBlock) {
    if (rBlock.IsPolyhedron())
        return 3;
    const std::string type(rBlock.Type());
    if (rBlock.IsRagged() || type.rfind("polygon", 0) == 0)
        return 2;
    const int dim = cell_type_dimension(cell_type_from_name(type));
    return dim == 2 || dim == 3 ? dim : -1;
}

// Every node of a cell (all faces' nodes of a polyhedron), unsorted.
template <class TView>
void cc_cell_nodes(const TView& rBlock, std::size_t Row, std::vector<std::int64_t>& rOut) {
    rOut.clear();
    if (rBlock.IsPolyhedron()) {
        for (std::size_t f = 0; f < rBlock.NumFaces(Row); ++f) {
            const auto face = rBlock.Face(Row, f);
            rOut.insert(rOut.end(), face.first, face.first + face.second);
        }
    } else if (rBlock.IsRagged()) {
        rOut.assign(rBlock.Row(Row), rBlock.Row(Row) + rBlock.RowSize(Row));
    } else {
        const std::size_t npc = rBlock.NodesPerCell();
        for (std::size_t k = 0; k < npc; ++k)
            rOut.push_back(detail::read_int(rBlock.Conn(), Row * npc + k));
    }
}

// 3-D: the faces of the volume cells that only one cell uses.
CcFacets cc_volume_facets(const Mesh& rMesh, std::vector<std::string>& rWarnings) {
    CcFacets out;
    bool all_tabulated = true;
    for (const auto cb : rMesh.CellRange()) {
        if (cc_block_dimension(cb) != 3)
            continue;
        if (cb.IsPolyhedron() || cb.IsRagged()) {
            all_tabulated = false;
            break;
        }
        const auto& faces = detail::cell_faces(cell_type_from_name(cb.Type()));
        for (const auto& f : faces)
            if (f.mNumCorners > 4) {
                all_tabulated = false;
                break;
            }
        if (!all_tabulated)
            break;
    }

    if (all_tabulated) {
        // Direct O(F) vertex-bucketed extraction by minimum face vertex id:
        // skips orientation, interior face copies, and global face numbering.
        const std::size_t npts = rMesh.NumPoints();
        const std::vector<std::int64_t> bases = detail::block_bases(rMesh);
        struct CcBlockInfo {
            std::size_t mBlockIdx = 0;
            std::size_t mNumCells = 0;
            std::size_t mNpc = 0;
            std::int64_t mBase = 0;
            const std::vector<detail::CellFaceDef>* mpFaces = nullptr;
            std::unique_ptr<detail::Int64View> mConn;
        };
        std::vector<CcBlockInfo> blocks;
        std::size_t b = 0;
        for (const auto cb : rMesh.CellRange()) {
            const std::size_t block = b++;
            if (cc_block_dimension(cb) != 3)
                continue;
            const auto& faces = detail::cell_faces(cell_type_from_name(cb.Type()));
            if (faces.empty()) {
                if (cb.NumCells() > 0)
                    rWarnings.push_back("block " + std::to_string(block) + " ('" +
                                        std::string(cb.Type()) +
                                        "') has no known face topology and was not examined");
                continue;
            }
            if (cb.NumCells() == 0)
                continue;
            blocks.push_back({block, cb.NumCells(), cb.NodesPerCell(), bases[block], &faces,
                              std::make_unique<detail::Int64View>(cb.Conn())});
        }
        if (blocks.empty())
            return out;

        auto sort3 = [](std::int64_t& a, std::int64_t& b_, std::int64_t& c) {
            if (a > b_)
                std::swap(a, b_);
            if (b_ > c)
                std::swap(b_, c);
            if (a > b_)
                std::swap(a, b_);
        };
        auto sort4 = [](std::int64_t& a, std::int64_t& b_, std::int64_t& c, std::int64_t& d) {
            if (a > b_)
                std::swap(a, b_);
            if (c > d)
                std::swap(c, d);
            if (a > c)
                std::swap(a, c);
            if (b_ > d)
                std::swap(b_, d);
            if (b_ > c)
                std::swap(b_, c);
        };

        // Pass 1: validate indices and count faces per minimum vertex id.
        std::vector<std::size_t> head(npts + 1, 0);
        for (const CcBlockInfo& bi : blocks) {
            const std::int64_t* pconn = bi.mConn->Data();
            const std::size_t npc = bi.mNpc;
            const auto& faces = *bi.mpFaces;
            for (std::size_t c = 0; c < bi.mNumCells; ++c) {
                const std::int64_t* row = pconn + c * npc;
                for (const auto& f : faces) {
                    std::int64_t lo = row[f.mNodes[0]];
                    for (std::size_t k = 0; k < f.mNumCorners; ++k) {
                        const std::int64_t v = row[f.mNodes[k]];
                        if (v < 0 || static_cast<std::size_t>(v) >= npts)
                            throw std::invalid_argument(
                                std::string(kCcPrefix) +
                                "a cell references a point the mesh does not have");
                        if (v < lo)
                            lo = v;
                    }
                    ++head[static_cast<std::size_t>(lo) + 1];
                }
            }
        }
        for (std::size_t i = 0; i < npts; ++i)
            head[i + 1] += head[i];

        struct CcSlot {
            std::int64_t mK1, mK2, mK3;
            std::int64_t mOwner;
            std::uint32_t mBi;
            std::uint32_t mRow;
            std::uint8_t mFace;
            std::uint8_t mNc;
        };
        std::vector<CcSlot> slots(head[npts]);
        std::vector<std::size_t> cursor(head);
        for (std::size_t bi_idx = 0; bi_idx < blocks.size(); ++bi_idx) {
            const CcBlockInfo& bi = blocks[bi_idx];
            const std::int64_t* pconn = bi.mConn->Data();
            const std::size_t npc = bi.mNpc;
            const auto& faces = *bi.mpFaces;
            const std::size_t nf_cell = faces.size();
            for (std::size_t c = 0; c < bi.mNumCells; ++c) {
                const std::int64_t* row = pconn + c * npc;
                const std::int64_t owner = bi.mBase + static_cast<std::int64_t>(c);
                for (std::size_t fi = 0; fi < nf_cell; ++fi) {
                    const auto& f = faces[fi];
                    std::int64_t k0 = row[f.mNodes[0]];
                    std::int64_t k1 = row[f.mNodes[1]];
                    std::int64_t k2 = row[f.mNodes[2]];
                    std::int64_t k3 = -1;
                    if (f.mNumCorners == 3) {
                        sort3(k0, k1, k2);
                    } else {
                        k3 = row[f.mNodes[3]];
                        sort4(k0, k1, k2, k3);
                    }
                    const std::size_t dst = cursor[static_cast<std::size_t>(k0)]++;
                    slots[dst] = {k1,
                                  k2,
                                  k3,
                                  owner,
                                  static_cast<std::uint32_t>(bi_idx),
                                  static_cast<std::uint32_t>(c),
                                  static_cast<std::uint8_t>(fi),
                                  static_cast<std::uint8_t>(f.mNumCorners)};
                }
            }
        }

        // Pass 2: sort each tiny per-vertex bucket and emit single-use facets.
        auto emit_boundary = [&](const CcSlot& rS) {
            const CcBlockInfo& bi = blocks[rS.mBi];
            const std::int64_t* row =
                bi.mConn->Data() + static_cast<std::size_t>(rS.mRow) * bi.mNpc;
            const auto& f = (*bi.mpFaces)[rS.mFace];
            std::int64_t ring[4];
            for (std::size_t k = 0; k < rS.mNc; ++k)
                ring[k] = row[f.mNodes[k]];
            out.Push(ring, rS.mNc, rS.mOwner);
        };
        for (std::size_t v = 0; v < npts; ++v) {
            const std::size_t beg = head[v], end = head[v + 1];
            if (beg == end)
                continue;
            if (end - beg == 1) {
                emit_boundary(slots[beg]);
                continue;
            }
            std::sort(slots.begin() + static_cast<std::ptrdiff_t>(beg),
                      slots.begin() + static_cast<std::ptrdiff_t>(end),
                      [](const CcSlot& rA, const CcSlot& rB) {
                          if (rA.mK1 != rB.mK1)
                              return rA.mK1 < rB.mK1;
                          if (rA.mK2 != rB.mK2)
                              return rA.mK2 < rB.mK2;
                          return rA.mK3 < rB.mK3;
                      });
            for (std::size_t i = beg; i < end;) {
                std::size_t j = i + 1;
                while (j < end && slots[j].mK1 == slots[i].mK1 && slots[j].mK2 == slots[i].mK2 &&
                       slots[j].mK3 == slots[i].mK3)
                    ++j;
                if (j - i == 1)
                    emit_boundary(slots[i]);
                else if (j - i >= 3)
                    ++out.mNumNonManifold;
                i = j;
            }
        }
        return out;
    }

    const detail::GlobalFaces gf = detail::build_global_faces(rMesh);
    out.mNumNonManifold = gf.mNumNonManifold;
    for (std::size_t f = 0; f < gf.NumFaces(); ++f)
        if (gf.mNeighbour[f] < 0)
            out.Push(gf.Face(f), gf.FaceSize(f),
                     gf.mCellToGlobal[static_cast<std::size_t>(gf.mOwner[f])]);
    std::size_t b = 0;
    for (const auto cb : rMesh.CellRange()) {
        if (cc_block_dimension(cb) == 3 &&
            std::find(gf.mNonCellBlocks.begin(), gf.mNonCellBlocks.end(), b) !=
                gf.mNonCellBlocks.end() &&
            cb.NumCells() > 0)
            rWarnings.push_back("block " + std::to_string(b) + " ('" + std::string(cb.Type()) +
                                "') has no known face topology and was not examined");
        ++b;
    }
    if (gf.mNumUnorientable > 0)
        rWarnings.push_back(std::to_string(gf.mNumUnorientable) +
                            " polyhedra are not closed and orientable; their faces were "
                            "examined as given");
    return out;
}

// 2-D: the edges of the 2-D cells that only one cell uses, in ascending
// (low, high) endpoint order.
CcFacets cc_surface_facets(const Mesh& rMesh, std::vector<std::string>& rWarnings) {
    struct Use {
        std::int64_t mLo, mHi, mCell, mA, mB;
    };
    std::vector<Use> uses;
    const std::vector<std::int64_t> bases = detail::block_bases(rMesh);
    const std::size_t npts = rMesh.NumPoints();
    std::vector<std::int64_t> ring;
    std::size_t b = 0;
    for (const auto cb : rMesh.CellRange()) {
        const std::size_t block = b++;
        if (cc_block_dimension(cb) != 2)
            continue;
        const std::string type(cb.Type());
        const bool polygon = cb.IsRagged() || type.rfind("polygon", 0) == 0;
        std::vector<std::size_t> corners;
        if (!polygon) {
            for (const detail::CellEdgeDef& ed : detail::cell_edges(cell_type_from_name(type)))
                corners.push_back(ed.mNodes[0]);
            if (corners.size() < 3) {
                if (cb.NumCells() > 0)
                    rWarnings.push_back("block " + std::to_string(block) + " ('" + type +
                                        "') has no known edge topology and was not examined");
                continue;
            }
        }
        const std::size_t npc = cb.NodesPerCell();
        for (std::size_t c = 0; c < cb.NumCells(); ++c) {
            ring.clear();
            if (cb.IsRagged())
                ring.assign(cb.Row(c), cb.Row(c) + cb.RowSize(c));
            else if (polygon)
                for (std::size_t k = 0; k < npc; ++k)
                    ring.push_back(detail::read_int(cb.Conn(), c * npc + k));
            else
                for (std::size_t k : corners)
                    ring.push_back(detail::read_int(cb.Conn(), c * npc + k));
            const std::int64_t g = bases[block] + static_cast<std::int64_t>(c);
            const std::size_t n = ring.size();
            for (std::size_t k = 0; k < n; ++k) {
                const std::int64_t a = ring[k], bb = ring[(k + 1) % n];
                if (a < 0 || bb < 0 || static_cast<std::size_t>(a) >= npts ||
                    static_cast<std::size_t>(bb) >= npts)
                    throw std::invalid_argument(std::string(kCcPrefix) +
                                                "a cell references a point the mesh does not have");
                if (a == bb)
                    continue;
                uses.push_back(Use{std::min(a, bb), std::max(a, bb), g, a, bb});
            }
        }
    }
    std::sort(uses.begin(), uses.end(), [](const Use& rX, const Use& rY) {
        if (rX.mLo != rY.mLo)
            return rX.mLo < rY.mLo;
        if (rX.mHi != rY.mHi)
            return rX.mHi < rY.mHi;
        return rX.mCell < rY.mCell;
    });
    CcFacets out;
    for (std::size_t i = 0; i < uses.size();) {
        std::size_t j = i + 1;
        while (j < uses.size() && uses[j].mLo == uses[i].mLo && uses[j].mHi == uses[i].mHi)
            ++j;
        if (j - i == 1) {
            const std::int64_t ids[2] = {uses[i].mA, uses[i].mB};
            out.Push(ids, 2, uses[i].mCell);
        } else if (j - i >= 3) {
            ++out.mNumNonManifold;
        }
        i = j;
    }
    return out;
}

detail::Vec3 cc_point(const std::vector<double>& rXyz, std::int64_t Id) {
    const std::size_t i = static_cast<std::size_t>(Id) * 3;
    return {{rXyz[i], rXyz[i + 1], rXyz[i + 2]}};
}

// The squared distance from P to segment AB.
double cc_segment_distance_sq(const detail::Vec3& rP, const detail::Vec3& rA,
                              const detail::Vec3& rB) {
    const detail::Vec3 ab = detail::vec3_sub(rB, rA);
    const double len_sq = detail::vec3_norm_sq(ab);
    double t = len_sq > 0.0 ? detail::vec3_dot(detail::vec3_sub(rP, rA), ab) / len_sq : 0.0;
    t = std::min(1.0, std::max(0.0, t));
    return detail::vec3_norm_sq(
        detail::vec3_sub(rP, detail::vec3_add(rA, detail::vec3_scale(ab, t))));
}

// The squared distance from P to a facet: a segment, or a polygon as the fan
// of triangles from its first corner.
double cc_facet_distance_sq(const detail::Vec3& rP, const std::vector<double>& rXyz,
                            const std::int64_t* pRing, std::size_t Corners) {
    if (Corners == 2)
        return cc_segment_distance_sq(rP, cc_point(rXyz, pRing[0]), cc_point(rXyz, pRing[1]));
    double best = std::numeric_limits<double>::infinity();
    const detail::Vec3 a = cc_point(rXyz, pRing[0]);
    for (std::size_t k = 1; k + 1 < Corners; ++k)
        best = std::min(best, detail::closest_point_on_triangle(rP, a, cc_point(rXyz, pRing[k]),
                                                                cc_point(rXyz, pRing[k + 1]))
                                  .mDistanceSq);
    return best;
}

// A facet's plane: its centroid, unit (Newell) normal and an in-plane frame.
struct CcPlane {
    detail::Vec3 mCentre{{0.0, 0.0, 0.0}};
    detail::Vec3 mNormal{{0.0, 0.0, 0.0}};
    detail::Vec3 mU{{0.0, 0.0, 0.0}};
    detail::Vec3 mV{{0.0, 0.0, 0.0}};
    bool mValid = false;
};

CcPlane cc_plane(const std::vector<double>& rXyz, const std::int64_t* pRing, std::size_t Corners) {
    CcPlane plane;
    for (std::size_t k = 0; k < Corners; ++k)
        plane.mCentre = detail::vec3_add(plane.mCentre, cc_point(rXyz, pRing[k]));
    plane.mCentre = detail::vec3_scale(plane.mCentre, 1.0 / static_cast<double>(Corners));
    detail::Vec3 n{{0.0, 0.0, 0.0}};
    for (std::size_t k = 0; k < Corners; ++k) {
        const detail::Vec3 a = cc_point(rXyz, pRing[k]);
        const detail::Vec3 b = cc_point(rXyz, pRing[(k + 1) % Corners]);
        n[0] += (a[1] - b[1]) * (a[2] + b[2]);
        n[1] += (a[2] - b[2]) * (a[0] + b[0]);
        n[2] += (a[0] - b[0]) * (a[1] + b[1]);
    }
    const double len = detail::vec3_norm(n);
    if (!(len > 0.0))
        return plane;
    plane.mNormal = detail::vec3_scale(n, 1.0 / len);
    const detail::Vec3 e = detail::vec3_sub(cc_point(rXyz, pRing[1]), cc_point(rXyz, pRing[0]));
    const double elen = detail::vec3_norm(e);
    if (!(elen > 0.0))
        return plane;
    plane.mU = detail::vec3_scale(e, 1.0 / elen);
    plane.mV = detail::vec3_cross(plane.mNormal, plane.mU);
    plane.mValid = true;
    return plane;
}

using CcUv = std::array<double, 2>;

double cc_area(const CcUv* pPoly, std::size_t N) {
    double a = 0.0;
    for (std::size_t i = 0; i < N; ++i) {
        const CcUv& p = pPoly[i];
        const CcUv& q = pPoly[i + 1 < N ? i + 1 : 0];
        a += p[0] * q[1] - q[0] * p[1];
    }
    return 0.5 * a;
}

// The area of the intersection of two triangles, both counter-clockwise, by
// clipping the first against the second's three edges (Sutherland and Hodgman).
// Clipping a 3-vertex polygon against 3 half-planes produces at most 6 vertices
// (capacity 8 on the stack, zero heap allocations).
double cc_triangle_overlap(const std::array<CcUv, 3>& rA, const std::array<CcUv, 3>& rB) {
    const double a_umin = std::min({rA[0][0], rA[1][0], rA[2][0]});
    const double a_umax = std::max({rA[0][0], rA[1][0], rA[2][0]});
    const double b_umin = std::min({rB[0][0], rB[1][0], rB[2][0]});
    const double b_umax = std::max({rB[0][0], rB[1][0], rB[2][0]});
    if (a_umax <= b_umin || b_umax <= a_umin)
        return 0.0;
    const double a_vmin = std::min({rA[0][1], rA[1][1], rA[2][1]});
    const double a_vmax = std::max({rA[0][1], rA[1][1], rA[2][1]});
    const double b_vmin = std::min({rB[0][1], rB[1][1], rB[2][1]});
    const double b_vmax = std::max({rB[0][1], rB[1][1], rB[2][1]});
    if (a_vmax <= b_vmin || b_vmax <= a_vmin)
        return 0.0;

    std::array<CcUv, 8> poly{rA[0], rA[1], rA[2]};
    std::array<CcUv, 8> next{};
    std::size_t n_poly = 3;
    for (std::size_t e = 0; e < 3 && n_poly > 0; ++e) {
        const CcUv& p = rB[e];
        const CcUv& q = rB[e + 1 < 3 ? e + 1 : 0];
        const double ex = q[0] - p[0];
        const double ey = q[1] - p[1];
        auto side = [&](const CcUv& rX) { return ex * (rX[1] - p[1]) - ey * (rX[0] - p[0]); };
        std::size_t n_next = 0;
        for (std::size_t i = 0; i < n_poly; ++i) {
            const CcUv& a = poly[i];
            const CcUv& b = poly[i + 1 < n_poly ? i + 1 : 0];
            const double sa = side(a), sb = side(b);
            if (sa >= 0.0)
                next[n_next++] = a;
            if ((sa >= 0.0) != (sb >= 0.0)) {
                const double t = sa / (sa - sb);
                next[n_next++] = {{a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])}};
            }
        }
        poly = next;
        n_poly = n_next;
    }
    return n_poly < 3 ? 0.0 : std::abs(cc_area(poly.data(), n_poly));
}

// Orient one fan triangle counter-clockwise; returns false if degenerate.
bool cc_orient_tri(const CcUv& rP0, const CcUv& rP1, const CcUv& rP2, std::array<CcUv, 3>& rOut) {
    rOut = {{rP0, rP1, rP2}};
    const double a = (rOut[1][0] - rOut[0][0]) * (rOut[2][1] - rOut[0][1]) -
                     (rOut[1][1] - rOut[0][1]) * (rOut[2][0] - rOut[0][0]);
    if (a < 0.0)
        std::swap(rOut[1], rOut[2]);
    return a != 0.0;
}

// A ring as the counter-clockwise fan of triangles from its first corner.
std::vector<std::array<CcUv, 3>> cc_fan(const CcUv* pRing, std::size_t N) {
    std::vector<std::array<CcUv, 3>> fan;
    if (N >= 3)
        fan.reserve(N - 2);
    for (std::size_t k = 1; k + 1 < N; ++k) {
        std::array<CcUv, 3> t;
        if (cc_orient_tri(pRing[0], pRing[k], pRing[k + 1], t))
            fan.push_back(t);
    }
    return fan;
}

// How two unmatched facets relate.
enum class CcPair { Apart, Unwelded, Overlap };

CcPair cc_classify_polygons(const std::vector<double>& rXyz, const std::int64_t* pF, std::size_t Nf,
                            const CcPlane& rPf, const std::int64_t* pG, std::size_t Ng,
                            const CcPlane& rPg, double Tol) {
    if (!rPf.mValid || !rPg.mValid)
        return CcPair::Apart;
    for (std::size_t k = 0; k < Ng; ++k)
        if (std::abs(detail::vec3_dot(rPf.mNormal,
                                      detail::vec3_sub(cc_point(rXyz, pG[k]), rPf.mCentre))) > Tol)
            return CcPair::Apart;
    for (std::size_t k = 0; k < Nf; ++k)
        if (std::abs(detail::vec3_dot(rPg.mNormal,
                                      detail::vec3_sub(cc_point(rXyz, pF[k]), rPg.mCentre))) > Tol)
            return CcPair::Apart;
    if (Nf == Ng) {
        bool same = true;
        for (std::size_t i = 0; i < Nf && same; ++i) {
            bool found = false;
            for (std::size_t j = 0; j < Ng && !found; ++j)
                found = detail::vec3_norm_sq(detail::vec3_sub(cc_point(rXyz, pF[i]),
                                                              cc_point(rXyz, pG[j]))) <= Tol * Tol;
            same = found;
        }
        if (same)
            return CcPair::Unwelded;
    }
    if (Nf == 3 && Ng == 3) {
        std::size_t shared = 0;
        std::int64_t s0 = -1, s1 = -1;
        std::size_t f_rem = 3;
        for (std::size_t i = 0; i < 3; ++i) {
            if (pF[i] == pG[0] || pF[i] == pG[1] || pF[i] == pG[2]) {
                if (shared == 0)
                    s0 = pF[i];
                else if (shared == 1)
                    s1 = pF[i];
                ++shared;
            } else {
                f_rem = i;
            }
        }
        if (shared == 2 && f_rem < 3) {
            std::size_t g_rem = 3;
            for (std::size_t j = 0; j < 3; ++j)
                if (pG[j] != s0 && pG[j] != s1) {
                    g_rem = j;
                    break;
                }
            if (g_rem < 3) {
                const detail::Vec3 p0 = cc_point(rXyz, s0);
                const detail::Vec3 e = detail::vec3_sub(cc_point(rXyz, s1), p0);
                const double sf = detail::vec3_dot(
                    detail::vec3_cross(e, detail::vec3_sub(cc_point(rXyz, pF[f_rem]), p0)),
                    rPf.mNormal);
                const double sg = detail::vec3_dot(
                    detail::vec3_cross(e, detail::vec3_sub(cc_point(rXyz, pG[g_rem]), p0)),
                    rPf.mNormal);
                if (sf * sg <= 0.0)
                    return CcPair::Apart;
            }
        }
    }
    auto uv = [&](std::int64_t Id) {
        const detail::Vec3 d = detail::vec3_sub(cc_point(rXyz, Id), rPf.mCentre);
        return CcUv{{detail::vec3_dot(d, rPf.mU), detail::vec3_dot(d, rPf.mV)}};
    };
    constexpr std::size_t kCcStackRing = 8;
    std::array<CcUv, kCcStackRing> f_stack{}, g_stack{};
    std::vector<CcUv> f_heap, g_heap;
    CcUv* pf = f_stack.data();
    CcUv* pg = g_stack.data();
    if (Nf > kCcStackRing) {
        f_heap.resize(Nf);
        pf = f_heap.data();
    }
    if (Ng > kCcStackRing) {
        g_heap.resize(Ng);
        pg = g_heap.data();
    }
    double f_umin = std::numeric_limits<double>::infinity();
    double f_umax = -f_umin, f_vmin = f_umin, f_vmax = -f_umin;
    for (std::size_t k = 0; k < Nf; ++k) {
        const CcUv p = uv(pF[k]);
        pf[k] = p;
        f_umin = std::min(f_umin, p[0]);
        f_umax = std::max(f_umax, p[0]);
        f_vmin = std::min(f_vmin, p[1]);
        f_vmax = std::max(f_vmax, p[1]);
    }
    double g_umin = std::numeric_limits<double>::infinity();
    double g_umax = -g_umin, g_vmin = g_umin, g_vmax = -g_umin;
    for (std::size_t k = 0; k < Ng; ++k) {
        const CcUv p = uv(pG[k]);
        pg[k] = p;
        g_umin = std::min(g_umin, p[0]);
        g_umax = std::max(g_umax, p[0]);
        g_vmin = std::min(g_vmin, p[1]);
        g_vmax = std::max(g_vmax, p[1]);
    }
    if (f_umax <= g_umin || g_umax <= f_umin || f_vmax <= g_vmin || g_vmax <= f_vmin)
        return CcPair::Apart;

    // The facets overlap when their intersection is a strip wider than the
    // tolerance: its area beyond the tolerance times the smaller facet's
    // diameter. Facets that only share a boundary intersect in a segment.
    double area = 0.0;
    if (Nf <= kCcStackRing && Ng <= kCcStackRing) {
        std::array<std::array<CcUv, 3>, kCcStackRing - 2> ff{}, gf{};
        std::size_t nff = 0, ngf = 0;
        for (std::size_t k = 1; k + 1 < Nf; ++k)
            if (cc_orient_tri(pf[0], pf[k], pf[k + 1], ff[nff]))
                ++nff;
        for (std::size_t k = 1; k + 1 < Ng; ++k)
            if (cc_orient_tri(pg[0], pg[k], pg[k + 1], gf[ngf]))
                ++ngf;
        for (std::size_t a = 0; a < nff; ++a)
            for (std::size_t b = 0; b < ngf; ++b)
                area += cc_triangle_overlap(ff[a], gf[b]);
    } else {
        const std::vector<std::array<CcUv, 3>> ff = cc_fan(pf, Nf), gf = cc_fan(pg, Ng);
        for (const auto& a : ff)
            for (const auto& b : gf)
                area += cc_triangle_overlap(a, b);
    }
    if (!(area > 0.0))
        return CcPair::Apart;
    auto diameter = [](const CcUv* pRing, std::size_t N) {
        double d = 0.0;
        for (std::size_t i = 0; i < N; ++i)
            for (std::size_t j = i + 1; j < N; ++j)
                d = std::max(d, std::hypot(pRing[i][0] - pRing[j][0], pRing[i][1] - pRing[j][1]));
        return d;
    };
    return area > Tol * std::min(diameter(pf, Nf), diameter(pg, Ng)) ? CcPair::Overlap
                                                                     : CcPair::Apart;
}

CcPair cc_classify_segments(const std::vector<double>& rXyz, const std::int64_t* pF,
                            const std::int64_t* pG, double Tol) {
    const detail::Vec3 a = cc_point(rXyz, pF[0]), b = cc_point(rXyz, pF[1]);
    const detail::Vec3 c = cc_point(rXyz, pG[0]), d = cc_point(rXyz, pG[1]);
    const detail::Vec3 ab = detail::vec3_sub(b, a), cd = detail::vec3_sub(d, c);
    const double lf = detail::vec3_norm(ab), lg = detail::vec3_norm(cd);
    if (!(lf > 0.0) || !(lg > 0.0))
        return CcPair::Apart;
    // Distance of a point from the infinite line through P0 along unit U.
    auto off_line = [&](const detail::Vec3& rP, const detail::Vec3& rP0, const detail::Vec3& rU) {
        const detail::Vec3 w = detail::vec3_sub(rP, rP0);
        return detail::vec3_norm(detail::vec3_cross(w, rU));
    };
    const detail::Vec3 uf = detail::vec3_scale(ab, 1.0 / lf), ug = detail::vec3_scale(cd, 1.0 / lg);
    if (off_line(c, a, uf) > Tol || off_line(d, a, uf) > Tol || off_line(a, c, ug) > Tol ||
        off_line(b, c, ug) > Tol)
        return CcPair::Apart;
    const double tol_sq = Tol * Tol;
    auto near = [&](const detail::Vec3& rX, const detail::Vec3& rY) {
        return detail::vec3_norm_sq(detail::vec3_sub(rX, rY)) <= tol_sq;
    };
    if ((near(a, c) && near(b, d)) || (near(a, d) && near(b, c)))
        return CcPair::Unwelded;
    double t0 = detail::vec3_dot(detail::vec3_sub(c, a), uf);
    double t1 = detail::vec3_dot(detail::vec3_sub(d, a), uf);
    if (t0 > t1)
        std::swap(t0, t1);
    return std::min(lf, t1) - std::max(0.0, t0) > Tol ? CcPair::Overlap : CcPair::Apart;
}

}  // namespace

ConformityResult check_conformity(const Mesh& rMesh, const ConformityOptions& rOptions) {
    if (!(rOptions.mTolerance >= 0.0) || !std::isfinite(rOptions.mTolerance))
        throw std::invalid_argument(std::string(kCcPrefix) +
                                    "the tolerance must be finite and non-negative");
    ConformityResult result;

    int dim = -1;
    for (const auto cb : rMesh.CellRange())
        if (cb.NumCells() > 0)
            dim = std::max(dim, cc_block_dimension(cb));
    CcFacets facets;
    if (dim == 3)
        facets = cc_volume_facets(rMesh, result.mWarnings);
    else if (dim == 2)
        facets = cc_surface_facets(rMesh, result.mWarnings);
    else
        result.mWarnings.push_back("the mesh has no 2-D or 3-D cells; nothing was examined");
    result.mDimension = dim < 2 ? 0 : dim;
    result.mNumNonManifoldFacets = facets.mNumNonManifold;
    const std::size_t nf = facets.Size();
    result.mNumBoundaryFacets = static_cast<std::int64_t>(nf);

    const std::size_t npts = rMesh.NumPoints();
    std::vector<std::int64_t> flags(npts, 0);
    if (nf > 0) {
        const std::vector<double> xyz = detail::surface_xyz(rMesh);

        // The facets' boxes and the tolerance, from the box of their corners.
        std::vector<std::array<double, 6>> box(nf);
        const double inf = std::numeric_limits<double>::infinity();
        std::array<double, 3> lo{{inf, inf, inf}};
        std::array<double, 3> hi{{-inf, -inf, -inf}};
        double extent_sum = 0.0;
        for (std::size_t f = 0; f < nf; ++f) {
            std::array<double, 6>& r_box = box[f];
            r_box = {{inf, inf, inf, -inf, -inf, -inf}};
            const std::int64_t* ring = facets.Ring(f);
            for (std::size_t k = 0; k < facets.Corners(f); ++k) {
                const detail::Vec3 p = cc_point(xyz, ring[k]);
                for (std::size_t d = 0; d < 3; ++d) {
                    r_box[d] = std::min(r_box[d], p[d]);
                    r_box[3 + d] = std::max(r_box[3 + d], p[d]);
                }
            }
            for (std::size_t d = 0; d < 3; ++d) {
                lo[d] = std::min(lo[d], r_box[d]);
                hi[d] = std::max(hi[d], r_box[3 + d]);
            }
            extent_sum += std::max({r_box[3] - r_box[0], r_box[4] - r_box[1], r_box[5] - r_box[2]});
        }
        const double diag =
            std::sqrt((hi[0] - lo[0]) * (hi[0] - lo[0]) + (hi[1] - lo[1]) * (hi[1] - lo[1]) +
                      (hi[2] - lo[2]) * (hi[2] - lo[2]));
        const double tol = rOptions.mTolerance * diag;
        for (std::array<double, 6>& r_box : box)
            for (std::size_t d = 0; d < 3; ++d) {
                r_box[d] -= tol;
                r_box[3 + d] += tol;
            }

        // A bucket grid over the facets' boxes, sized by their mean extent and
        // never finer than 1/256 of the diagonal, so a large facet among small
        // ones spans a bounded number of buckets.
        double cell = std::max(extent_sum / static_cast<double>(nf), diag / 256.0);
        if (!(cell > 0.0))
            cell = 1.0;
        std::vector<detail::GridKey> box_lo(nf);
        std::vector<detail::GridKey> box_hi(nf);
        parallel_for(nf, [&](std::size_t f) {
            box_lo[f] = detail::bucket_key(box[f].data(), cell);
            box_hi[f] = detail::bucket_key(box[f].data() + 3, cell);
        });
        const detail::BucketTable grid = detail::build_bucket_table(cell, box_lo, box_hi);

        // One view per block, and the nodes of a facet's owner on demand.
        const std::vector<std::int64_t> bases = detail::block_bases(rMesh);
        auto owner_has = [&](std::size_t Facet, std::int64_t Point,
                             std::vector<std::int64_t>& rScratch) {
            const auto br = detail::global_to_block_row(bases, facets.mOwner[Facet]);
            cc_cell_nodes(rMesh.Cells(br.first), static_cast<std::size_t>(br.second), rScratch);
            return std::find(rScratch.begin(), rScratch.end(), Point) != rScratch.end();
        };

        // Candidates: the nodes of the cells that own an unmatched facet.
        std::vector<std::int64_t> owners(facets.mOwner);
        std::sort(owners.begin(), owners.end());
        owners.erase(std::unique(owners.begin(), owners.end()), owners.end());
        std::vector<std::int64_t> candidates;
        {
            std::vector<std::int64_t> nodes;
            for (std::int64_t g : owners) {
                const auto br = detail::global_to_block_row(bases, g);
                cc_cell_nodes(rMesh.Cells(br.first), static_cast<std::size_t>(br.second), nodes);
                for (std::int64_t p : nodes)
                    if (p >= 0 && static_cast<std::size_t>(p) < npts)
                        candidates.push_back(p);
            }
        }
        std::sort(candidates.begin(), candidates.end());
        candidates.erase(std::unique(candidates.begin(), candidates.end()), candidates.end());

        const std::int64_t flag = static_cast<std::int64_t>(dim == 3 ? ConformityFlag::HangingNode
                                                                     : ConformityFlag::TJunction);
        const double tol_sq = tol * tol;
        parallel_for(candidates.size(), [&](std::size_t i) {
            const std::int64_t p = candidates[i];
            const detail::Vec3 x = cc_point(xyz, p);
            detail::BucketView bucket(nullptr, 0);
            if (!grid.Find(grid.KeyOf(x.data()), bucket))
                return;
            static thread_local std::vector<std::int64_t> scratch;
            for (std::int64_t id : bucket) {
                const std::size_t f = static_cast<std::size_t>(id);
                const std::array<double, 6>& b = box[f];
                if (x[0] < b[0] || x[1] < b[1] || x[2] < b[2] || x[0] > b[3] || x[1] > b[4] ||
                    x[2] > b[5])
                    continue;
                const std::int64_t* ring = facets.Ring(f);
                const std::size_t nc = facets.Corners(f);
                if (std::find(ring, ring + nc, p) != ring + nc)
                    continue;
                if (cc_facet_distance_sq(x, xyz, ring, nc) > tol_sq)
                    continue;
                // On a corner's location: an unwelded duplicate, not a hanging node.
                bool on_corner = false;
                for (std::size_t k = 0; k < nc && !on_corner; ++k)
                    on_corner =
                        detail::vec3_norm_sq(detail::vec3_sub(x, cc_point(xyz, ring[k]))) <= tol_sq;
                if (on_corner || owner_has(f, p, scratch))
                    continue;
                flags[static_cast<std::size_t>(p)] = flag;
                return;
            }
        });

        // Facet pairs, each counted once at its lower facet.
        std::vector<CcPlane> planes(nf);
        if (dim == 3)
            parallel_for(nf, [&](std::size_t f) {
                planes[f] = cc_plane(xyz, facets.Ring(f), facets.Corners(f));
            });
        std::vector<std::int64_t> partial(nf, 0), unwelded(nf, 0);
        parallel_for(nf, [&](std::size_t f) {
            const std::array<double, 6>& bf = box[f];
            auto test_pair = [&](std::size_t g) {
                const std::array<double, 6>& bg = box[g];
                if (bf[0] > bg[3] || bg[0] > bf[3] || bf[1] > bg[4] || bg[1] > bf[4] ||
                    bf[2] > bg[5] || bg[2] > bf[5])
                    return;
                const CcPair kind =
                    dim == 3
                        ? cc_classify_polygons(xyz, facets.Ring(f), facets.Corners(f), planes[f],
                                               facets.Ring(g), facets.Corners(g), planes[g], tol)
                        : cc_classify_segments(xyz, facets.Ring(f), facets.Ring(g), tol);
                if (kind == CcPair::Overlap)
                    ++partial[f];
                else if (kind == CcPair::Unwelded)
                    ++unwelded[f];
            };
            if (box_lo[f] == box_hi[f]) {
                detail::BucketView bucket(nullptr, 0);
                if (grid.Find(box_lo[f], bucket))
                    for (std::int64_t id : bucket)
                        if (static_cast<std::size_t>(id) > f)
                            test_pair(static_cast<std::size_t>(id));
                return;
            }
            static thread_local std::vector<std::int64_t> near;
            near.clear();
            grid.ForEachInBox(box_lo[f], box_hi[f], [&](const detail::BucketView& rIds) {
                for (std::int64_t id : rIds)
                    if (static_cast<std::size_t>(id) > f)
                        near.push_back(id);
            });
            std::sort(near.begin(), near.end());
            near.erase(std::unique(near.begin(), near.end()), near.end());
            for (std::int64_t id : near)
                test_pair(static_cast<std::size_t>(id));
        });
        for (std::size_t f = 0; f < nf; ++f) {
            result.mNumPartialOverlaps += partial[f];
            result.mNumUnweldedFacets += unwelded[f];
        }
    }

    for (std::int64_t v : flags) {
        result.mNumHangingNodes += v == static_cast<std::int64_t>(ConformityFlag::HangingNode);
        result.mNumTJunctions += v == static_cast<std::int64_t>(ConformityFlag::TJunction);
    }
    result.mConforming = result.mNumHangingNodes == 0 && result.mNumTJunctions == 0 &&
                         result.mNumPartialOverlaps == 0 && result.mNumUnweldedFacets == 0;
    if (result.mNumPartialOverlaps > 0 || result.mNumHangingNodes > 0 || result.mNumTJunctions > 0)
        result.mInterfaceMode = "proximity";
    else if (result.mNumUnweldedFacets > 0)
        result.mInterfaceMode = "weld";

    result.mMesh = detail::clone_mesh(
        rMesh, [&](DataLocation Location, const std::string& rName, std::string&) {
            return !(rOptions.mAttachFlags && Location == DataLocation::Point &&
                     rName == kConformityHangingName);
        });
    if (rOptions.mAttachFlags) {
        NDArray out = NDArray::Uninit(DType::Int64, {npts});
        std::copy(flags.begin(), flags.end(), out.As<std::int64_t>());
        result.mMesh.AddPointData(kConformityHangingName, std::move(out));
    }
    for (std::size_t i = 0; i < rMesh.NumRegions(); ++i)
        result.mMesh.AddRegion(rMesh.Region(i));
    return result;
}

}  // namespace meshioplusplus

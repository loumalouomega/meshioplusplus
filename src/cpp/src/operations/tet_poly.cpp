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
// Tetrahedral <-> polyhedral conversion: the median dual (`tet_to_poly`) and
// simplexify with exact-tetrahedron recognition (`poly_to_tet`) -- see
// operations/tet_poly.hpp for the contract.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/tet_poly.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/face_mesh.hpp"
#include "meshioplusplus/detail/region_remap.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/parallel.hpp"

// Project includes (private, not installed)
#include "../detail/simplexify.hpp"
#include "../detail/typed_view.hpp"

namespace meshioplusplus {

namespace {

constexpr const char* kTpPrefix = "meshio++: tet_to_poly: ";
constexpr double kTpPi = 3.14159265358979323846;

using TpEdge = std::pair<std::int64_t, std::int64_t>;

TpEdge tp_edge(std::int64_t a, std::int64_t b) {
    return a < b ? TpEdge{a, b} : TpEdge{b, a};
}

/// Everything the per-vertex construction reads, built once.
struct TpContext {
    const detail::GlobalFaces* mpFaces = nullptr;
    std::size_t mNumPoints = 0;
    std::vector<std::size_t> mBlockOf;                  // compact tet -> input block
    std::vector<std::int64_t> mTetNodes;                // 4 per compact tet
    std::vector<std::int64_t> mTetFaceId;               // 4 global face ids per tet
    std::vector<std::array<std::int64_t, 3>> mTetRing;  // 4 rings per tet, wound outward
    std::vector<std::int64_t> mVertTetStart;            // CSR vertex -> tets (ascending)
    std::vector<std::int64_t> mVertTets;
    std::vector<std::int64_t> mFaceTag;  // face -> input triangle block, or -1
    std::vector<TpEdge> mMidEdges;       // sorted: edges that get a midpoint
    const double* mpXyz = nullptr;       // points as float64, 3 per point
    double mCosFeature = 1.0;

    std::int64_t MidPoint(std::int64_t a, std::int64_t b) const {
        const TpEdge e = tp_edge(a, b);
        const auto it = std::lower_bound(mMidEdges.begin(), mMidEdges.end(), e);
        if (it == mMidEdges.end() || *it != e)
            throw std::logic_error(std::string(kTpPrefix) + "missing boundary edge midpoint");
        return static_cast<std::int64_t>(mNumPoints) + (it - mMidEdges.begin());
    }
    std::int64_t FacePoint(std::int64_t f) const {
        return static_cast<std::int64_t>(mNumPoints + mMidEdges.size()) + f;
    }
    std::int64_t TetPoint(std::int64_t t) const {
        return static_cast<std::int64_t>(mNumPoints + mMidEdges.size() + mpFaces->NumFaces()) + t;
    }
    /// The tet on the other side of face @p f from tet @p t; -1 on the boundary.
    std::int64_t Across(std::int64_t f, std::int64_t t) const {
        const auto uf = static_cast<std::size_t>(f);
        return mpFaces->mOwner[uf] == t ? mpFaces->mNeighbour[uf] : mpFaces->mOwner[uf];
    }
    /// True when both sides of @p f are tets of the same block.
    bool Interior(std::int64_t f) const {
        const auto uf = static_cast<std::size_t>(f);
        const std::int64_t nb = mpFaces->mNeighbour[uf];
        return nb >= 0 && mBlockOf[static_cast<std::size_t>(mpFaces->mOwner[uf])] ==
                              mBlockOf[static_cast<std::size_t>(nb)];
    }
    /// Local face @p k of tet @p t: its global id and its ring, wound outward
    /// in the mesh's one consistent orientation.
    std::int64_t FaceId(std::int64_t t, int k) const {
        return mTetFaceId[4 * static_cast<std::size_t>(t) + static_cast<std::size_t>(k)];
    }
    const std::array<std::int64_t, 3>& Ring(std::int64_t t, int k) const {
        return mTetRing[4 * static_cast<std::size_t>(t) + static_cast<std::size_t>(k)];
    }
};

bool tp_contains(const std::array<std::int64_t, 3>& rTri, std::int64_t v) {
    return rTri[0] == v || rTri[1] == v || rTri[2] == v;
}

/// True when the cyclic triangle @p rTri traverses the directed edge v -> w.
bool tp_has_directed(const std::array<std::int64_t, 3>& rTri, std::int64_t v, std::int64_t w) {
    for (int k = 0; k < 3; ++k)
        if (rTri[k] == v && rTri[(k + 1) % 3] == w)
            return true;
    return false;
}

/// A boundary triangle around the vertex, wound outward from the component's
/// side and rotated so the vertex comes first: (v, a, b).
struct TpBoundaryTri {
    std::int64_t mFace = 0;
    std::int64_t mA = 0;
    std::int64_t mB = 0;
    std::int64_t mKind = -1;  // block on the other side, or -1 = mesh boundary
    std::int64_t mTag = -1;   // input triangle block, or -1
    std::int64_t mNext = -1;  // local index of the next triangle around v
    std::int64_t mPrev = -1;
    std::array<double, 3> mNormal{0.0, 0.0, 0.0};
};

/// One dual cell, as local face CSR.
struct TpCell {
    std::size_t mBlock = 0;
    std::vector<std::int64_t> mNodes;
    std::vector<std::int64_t> mFaceStart{0};
    std::vector<std::int64_t> mTets;  // compact ids, ascending
    std::int64_t mNumBoundaryFaces = 0;

    void AddFace(const std::vector<std::int64_t>& rRing) {
        mNodes.insert(mNodes.end(), rRing.begin(), rRing.end());
        mFaceStart.push_back(static_cast<std::int64_t>(mNodes.size()));
    }
};

/// One dual boundary polygon destined for an output polygon block.
struct TpPolygon {
    std::size_t mTag = 0;
    std::vector<std::int64_t> mRing;
    std::vector<std::int64_t> mFaces;  // the input triangles it covers
};

struct TpVertexOut {
    std::vector<TpCell> mCells;
    std::vector<TpPolygon> mPolygons;
};

std::array<double, 3> tp_unit_normal(const TpContext& rCtx, std::int64_t v, std::int64_t a,
                                     std::int64_t b) {
    const double* p = rCtx.mpXyz + 3 * v;
    const double* q = rCtx.mpXyz + 3 * a;
    const double* r = rCtx.mpXyz + 3 * b;
    const double u[3] = {q[0] - p[0], q[1] - p[1], q[2] - p[2]};
    const double w[3] = {r[0] - p[0], r[1] - p[1], r[2] - p[2]};
    std::array<double, 3> n{u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2],
                            u[0] * w[1] - u[1] * w[0]};
    const double len = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    if (len > 0.0)
        for (double& c : n)
            c /= len;
    return n;
}

/// Signed volume of the tet over the four point ids @p pN (connectivity order).
double tp_signed_volume(const double* pXyz, const std::int64_t* pN) {
    const double* a = pXyz + 3 * pN[0];
    double u[3][3];
    for (int r = 0; r < 3; ++r)
        for (int d = 0; d < 3; ++d)
            u[r][d] = pXyz[3 * pN[r + 1] + d] - a[d];
    return (u[0][0] * (u[1][1] * u[2][2] - u[1][2] * u[2][1]) -
            u[0][1] * (u[1][0] * u[2][2] - u[1][2] * u[2][0]) +
            u[0][2] * (u[1][0] * u[2][1] - u[1][1] * u[2][0])) /
           6.0;
}

/// Find-with-path-halving over a small local forest.
std::size_t tp_find(std::vector<std::size_t>& rParent, std::size_t i) {
    while (rParent[i] != i) {
        rParent[i] = rParent[rParent[i]];
        i = rParent[i];
    }
    return i;
}

/// Build every dual cell (and destined boundary polygon) around vertex @p v.
TpVertexOut tp_build_vertex(const TpContext& rCtx, std::int64_t v) {
    TpVertexOut out;
    const auto uv = static_cast<std::size_t>(v);
    const std::int64_t* tets = rCtx.mVertTets.data() + rCtx.mVertTetStart[uv];
    const std::size_t ntets =
        static_cast<std::size_t>(rCtx.mVertTetStart[uv + 1] - rCtx.mVertTetStart[uv]);
    if (ntets == 0)
        return out;
    auto local_of = [&](std::int64_t t) -> std::size_t {
        const std::int64_t* it = std::lower_bound(tets, tets + ntets, t);
        if (it == tets + ntets || *it != t)
            throw std::logic_error(std::string(kTpPrefix) + "tet adjacency out of step");
        return static_cast<std::size_t>(it - tets);
    };

    // Components: tets around v joined through faces at v interior to a block.
    std::vector<std::size_t> parent(ntets);
    std::iota(parent.begin(), parent.end(), std::size_t{0});
    for (std::size_t i = 0; i < ntets; ++i) {
        const std::int64_t t = tets[i];
        for (int k = 0; k < 4; ++k) {
            const std::int64_t f = rCtx.FaceId(t, k);
            if (!tp_contains(rCtx.Ring(t, k), v) || !rCtx.Interior(f))
                continue;
            const std::size_t j = local_of(rCtx.Across(f, t));
            const std::size_t ri = tp_find(parent, i);
            const std::size_t rj = tp_find(parent, j);
            if (ri != rj)
                parent[std::max(ri, rj)] = std::min(ri, rj);
        }
    }
    // Roots are each component's smallest local index, so walking i upward
    // meets the components in ascending order of their smallest tet.
    std::vector<std::vector<std::size_t>> comps;
    std::vector<std::int64_t> comp_of_root(ntets, -1);
    for (std::size_t i = 0; i < ntets; ++i) {
        const std::size_t r = tp_find(parent, i);
        if (comp_of_root[r] < 0) {
            comp_of_root[r] = static_cast<std::int64_t>(comps.size());
            comps.emplace_back();
        }
        comps[static_cast<std::size_t>(comp_of_root[r])].push_back(i);
    }

    std::vector<std::int64_t> ring;
    for (const std::vector<std::size_t>& members : comps) {
        TpCell cell;
        const std::int64_t t0 = tets[members.front()];
        cell.mBlock = rCtx.mBlockOf[static_cast<std::size_t>(t0)];
        for (const std::size_t i : members)
            cell.mTets.push_back(tets[i]);

        // The boundary triangles of this component at v.
        std::vector<TpBoundaryTri> bts;
        auto bt_of_face = [&](std::int64_t f) -> std::int64_t {
            for (std::size_t k = 0; k < bts.size(); ++k)
                if (bts[k].mFace == f)
                    return static_cast<std::int64_t>(k);
            return -1;
        };
        std::vector<std::int64_t> others;  // the w of every edge (v, w)
        for (const std::int64_t t : cell.mTets) {
            for (int k = 0; k < 4; ++k) {
                const std::int64_t n = rCtx.mTetNodes[4 * static_cast<std::size_t>(t) + k];
                if (n != v)
                    others.push_back(n);
                const std::int64_t f = rCtx.FaceId(t, k);
                const std::array<std::int64_t, 3>& tri = rCtx.Ring(t, k);
                if (!tp_contains(tri, v) || rCtx.Interior(f))
                    continue;
                TpBoundaryTri bt;
                bt.mFace = f;
                const int at = tri[0] == v ? 0 : (tri[1] == v ? 1 : 2);
                bt.mA = tri[(at + 1) % 3];
                bt.mB = tri[(at + 2) % 3];
                const std::int64_t other = rCtx.Across(f, t);
                bt.mKind =
                    other < 0
                        ? -1
                        : static_cast<std::int64_t>(rCtx.mBlockOf[static_cast<std::size_t>(other)]);
                bt.mTag = rCtx.mFaceTag[static_cast<std::size_t>(f)];
                bt.mNormal = tp_unit_normal(rCtx, v, bt.mA, bt.mB);
                bts.push_back(bt);
            }
        }
        std::sort(others.begin(), others.end());
        others.erase(std::unique(others.begin(), others.end()), others.end());

        // --- interior faces: one per fan of tets around each edge (v, w) ---
        for (const std::int64_t w : others) {
            // Each candidate tet's two faces through (v, w): the one its outward
            // winding traverses v -> w (counter-clockwise about v -> w) and the
            // one traversing w -> v.
            struct Cand {
                std::int64_t mTet;
                std::int64_t mCcw;
                std::int64_t mCw;
                bool mSeen;
            };
            std::vector<Cand> cands;
            for (const std::int64_t t : cell.mTets) {
                const std::int64_t* tn = rCtx.mTetNodes.data() + 4 * static_cast<std::size_t>(t);
                if (tn[0] != w && tn[1] != w && tn[2] != w && tn[3] != w)
                    continue;
                Cand c{t, -1, -1, false};
                for (int k = 0; k < 4; ++k) {
                    const std::array<std::int64_t, 3>& tri = rCtx.Ring(t, k);
                    if (tp_has_directed(tri, v, w))
                        c.mCcw = rCtx.FaceId(t, k);
                    else if (tp_has_directed(tri, w, v))
                        c.mCw = rCtx.FaceId(t, k);
                }
                if (c.mCcw < 0 || c.mCw < 0)
                    throw std::logic_error(std::string(kTpPrefix) + "tetrahedron face winding");
                cands.push_back(c);
            }
            auto cand_of = [&](std::int64_t t) -> Cand& {
                for (Cand& c : cands)
                    if (c.mTet == t)
                        return c;
                throw std::logic_error(std::string(kTpPrefix) + "fan left its component");
            };
            // Open fans first: from a tet whose clockwise face is a boundary.
            for (Cand& start : cands) {
                if (start.mSeen || rCtx.Interior(start.mCw))
                    continue;
                ring.clear();
                ring.push_back(rCtx.FacePoint(start.mCw));
                Cand* cur = &start;
                for (std::size_t guard = 0; guard <= cands.size(); ++guard) {
                    cur->mSeen = true;
                    ring.push_back(rCtx.TetPoint(cur->mTet));
                    ring.push_back(rCtx.FacePoint(cur->mCcw));
                    if (!rCtx.Interior(cur->mCcw)) {
                        const std::int64_t s = bt_of_face(start.mCw);
                        const std::int64_t e = bt_of_face(cur->mCcw);
                        if (s >= 0 && e >= 0) {
                            bts[static_cast<std::size_t>(s)].mNext = e;
                            bts[static_cast<std::size_t>(e)].mPrev = s;
                        }
                        break;
                    }
                    cur = &cand_of(rCtx.Across(cur->mCcw, cur->mTet));
                }
                ring.push_back(rCtx.MidPoint(v, w));
                cell.AddFace(ring);
            }
            // Then the closed fans that remain.
            for (Cand& start : cands) {
                if (start.mSeen)
                    continue;
                ring.clear();
                Cand* cur = &start;
                for (std::size_t guard = 0; guard <= cands.size(); ++guard) {
                    cur->mSeen = true;
                    ring.push_back(rCtx.TetPoint(cur->mTet));
                    ring.push_back(rCtx.FacePoint(cur->mCcw));
                    cur = &cand_of(rCtx.Across(cur->mCcw, cur->mTet));
                    if (cur->mSeen)
                        break;
                }
                cell.AddFace(ring);
            }
        }

        // --- boundary faces: runs of boundary quads around v ---
        std::vector<bool> visited(bts.size(), false);
        auto brk = [&](const TpBoundaryTri& rX, const TpBoundaryTri& rY) {
            if (rX.mKind != rY.mKind || rX.mTag != rY.mTag)
                return true;
            const double d = rX.mNormal[0] * rY.mNormal[0] + rX.mNormal[1] * rY.mNormal[1] +
                             rX.mNormal[2] * rY.mNormal[2];
            return d < rCtx.mCosFeature;
        };
        auto emit_run = [&](const std::vector<std::size_t>& rRun, bool Closed) {
            ring.clear();
            if (!Closed)
                ring.push_back(v);
            for (const std::size_t k : rRun) {
                ring.push_back(rCtx.MidPoint(v, bts[k].mA));
                ring.push_back(rCtx.FacePoint(bts[k].mFace));
            }
            if (!Closed)
                ring.push_back(rCtx.MidPoint(v, bts[rRun.back()].mB));
            cell.AddFace(ring);
            ++cell.mNumBoundaryFaces;
            const TpBoundaryTri& first = bts[rRun.front()];
            if (first.mTag >= 0 &&
                (first.mKind < 0 || static_cast<std::int64_t>(cell.mBlock) < first.mKind)) {
                TpPolygon poly;
                poly.mTag = static_cast<std::size_t>(first.mTag);
                poly.mRing = ring;
                for (const std::size_t k : rRun)
                    poly.mFaces.push_back(bts[k].mFace);
                out.mPolygons.push_back(std::move(poly));
            }
        };
        auto emit_sequence = [&](const std::vector<std::size_t>& rSeq, bool Cyclic) {
            const std::size_t n = rSeq.size();
            std::vector<bool> cut(n, false);  // cut[k]: a break after rSeq[k]
            bool any = !Cyclic;
            for (std::size_t k = 0; k + 1 < n || (Cyclic && k < n); ++k) {
                cut[k] = brk(bts[rSeq[k]], bts[rSeq[(k + 1) % n]]);
                any = any || cut[k];
            }
            if (!any) {
                emit_run(rSeq, /*Closed=*/true);
                return;
            }
            // Start just after a break; a chain starts at its own first element.
            std::size_t start = 0;
            if (Cyclic) {
                while (!cut[start])
                    ++start;
                start = (start + 1) % n;
            }
            std::vector<std::size_t> run;
            for (std::size_t s = 0; s < n; ++s) {
                const std::size_t k = (start + s) % n;
                run.push_back(rSeq[k]);
                if (s + 1 == n || cut[k]) {
                    emit_run(run, /*Closed=*/false);
                    run.clear();
                }
            }
        };
        std::vector<std::size_t> seq;
        // Chains (a triangle with no predecessor) first, then the cycles.
        for (int pass = 0; pass < 2; ++pass) {
            for (std::size_t k0 = 0; k0 < bts.size(); ++k0) {
                if (visited[k0] || (pass == 0 && bts[k0].mPrev >= 0))
                    continue;
                seq.clear();
                std::int64_t k = static_cast<std::int64_t>(k0);
                while (k >= 0 && !visited[static_cast<std::size_t>(k)]) {
                    visited[static_cast<std::size_t>(k)] = true;
                    seq.push_back(static_cast<std::size_t>(k));
                    k = bts[static_cast<std::size_t>(k)].mNext;
                }
                const bool cyclic = pass == 1 && k == static_cast<std::int64_t>(k0);
                emit_sequence(seq, cyclic);
            }
        }
        out.mCells.push_back(std::move(cell));
    }
    return out;
}

/// |volume| of compact tet @p t.
double tp_tet_volume(const TpContext& rCtx, std::size_t t) {
    return std::abs(tp_signed_volume(rCtx.mpXyz, rCtx.mTetNodes.data() + 4 * t));
}

/// Area of global face @p f.
double tp_face_area(const TpContext& rCtx, std::size_t f) {
    const std::int64_t* n = rCtx.mpFaces->Face(f);
    const std::array<double, 3> raw = [&] {
        const double* p = rCtx.mpXyz + 3 * n[0];
        const double* q = rCtx.mpXyz + 3 * n[1];
        const double* r = rCtx.mpXyz + 3 * n[2];
        const double u[3] = {q[0] - p[0], q[1] - p[1], q[2] - p[2]};
        const double w[3] = {r[0] - p[0], r[1] - p[1], r[2] - p[2]};
        return std::array<double, 3>{u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2],
                                     u[0] * w[1] - u[1] * w[0]};
    }();
    return 0.5 * std::sqrt(raw[0] * raw[0] + raw[1] * raw[1] + raw[2] * raw[2]);
}

bool tp_is_float(DType Dt) {
    return Dt == DType::Float32 || Dt == DType::Float64;
}

/// Reduce @p rSrc's rows over each output cell's members: a weighted mean for a
/// floating-point array, the heaviest member's row (ties: the first) otherwise.
/// @p rMembers is CSR over output cells of source rows, @p rWeights parallel.
NDArray tp_reduce_rows(const NDArray& rSrc, const std::vector<std::int64_t>& rStart,
                       const std::vector<std::int64_t>& rRows,
                       const std::vector<double>& rWeights) {
    const std::size_t nout = rStart.empty() ? 0 : rStart.size() - 1;
    const std::size_t in_rows = detail::rows(rSrc);
    const std::size_t ncomp = in_rows == 0 ? 0 : rSrc.Size() / in_rows;
    std::vector<std::size_t> shape = rSrc.Shape();
    shape[0] = nout;
    NDArray dst = NDArray::Uninit(rSrc.Dtype(), std::move(shape));
    if (tp_is_float(rSrc.Dtype())) {
        const detail::DoubleView src(rSrc);
        detail::DoubleSink sink(dst, 0, nout * ncomp);
        double* const vals = sink.Data();
        parallel_for(nout, [&](std::size_t c) {
            const auto b = static_cast<std::size_t>(rStart[c]);
            const auto e = static_cast<std::size_t>(rStart[c + 1]);
            double wsum = 0.0;
            for (std::size_t j = b; j < e; ++j)
                wsum += rWeights[j];
            for (std::size_t k = 0; k < ncomp; ++k) {
                double acc = 0.0;
                for (std::size_t j = b; j < e; ++j) {
                    const double wj = wsum > 0.0 ? rWeights[j] : 1.0;
                    acc += wj * src[static_cast<std::size_t>(rRows[j]) * ncomp + k];
                }
                const double denom = wsum > 0.0 ? wsum : static_cast<double>(e - b);
                vals[c * ncomp + k] = denom > 0.0 ? acc / denom : 0.0;
            }
        });
        sink.Commit();
    } else {
        const std::size_t row_bytes = in_rows == 0 ? 0 : rSrc.Nbytes() / in_rows;
        const std::byte* src = rSrc.Data();
        std::byte* out = dst.Data();
        parallel_for(nout, [&](std::size_t c) {
            const auto b = static_cast<std::size_t>(rStart[c]);
            const auto e = static_cast<std::size_t>(rStart[c + 1]);
            std::size_t best = b;
            for (std::size_t j = b + 1; j < e; ++j)
                if (rWeights[j] > rWeights[best])
                    best = j;
            if (b < e)
                std::memcpy(out + c * row_bytes,
                            src + static_cast<std::size_t>(rRows[best]) * row_bytes, row_bytes);
        });
    }
    return dst;
}

/// Gather rows of @p rSrc: output row i is source row @p rRows[i].
NDArray tp_gather_rows(const NDArray& rSrc, const std::vector<std::int64_t>& rRows) {
    const std::size_t in_rows = detail::rows(rSrc);
    std::vector<std::size_t> shape = rSrc.Shape();
    shape[0] = rRows.size();
    NDArray dst = NDArray::Uninit(rSrc.Dtype(), std::move(shape));
    const std::size_t row_bytes = in_rows == 0 ? 0 : rSrc.Nbytes() / in_rows;
    const std::byte* src = rSrc.Data();
    std::byte* out = dst.Data();
    parallel_for(rRows.size(), [&](std::size_t i) {
        std::memcpy(out + i * row_bytes, src + static_cast<std::size_t>(rRows[i]) * row_bytes,
                    row_bytes);
    });
    return dst;
}

/// One staged output block, with what its cell_data reduces over.
struct TpOutBlock {
    std::size_t mInBlock = 0;
    bool mIsPolygon = false;
    std::vector<std::int64_t> mFlat;
    std::vector<std::int64_t> mFaceStart{0};  // polygon: row offsets
    std::vector<std::int64_t> mCellStart{0};  // polyhedron: cell -> face offsets
    std::vector<std::int64_t> mVertex;        // per output cell
    std::vector<std::int64_t> mMemberStart{0};
    std::vector<std::int64_t> mMemberRow;  // source row in mInBlock
    std::vector<double> mMemberWeight;
    std::size_t NumCells() const { return mVertex.size(); }
};

}  // namespace

// --- tet_to_poly -------------------------------------------------------------

TetToPolyResult tet_to_poly(const Mesh& rMesh, const TetToPolyOptions& rOptions) {
    if (!(rOptions.mFeatureAngleDeg >= 0.0 && rOptions.mFeatureAngleDeg < 180.0))
        throw std::invalid_argument(std::string(kTpPrefix) +
                                    "mFeatureAngleDeg must lie in [0, 180)");

    // --- validate: every volume block linear tetra --------------------------
    const std::size_t nblocks_in = rMesh.NumCellBlocks();
    std::vector<bool> is_tet(nblocks_in, false);
    std::vector<bool> is_tri(nblocks_in, false);
    std::size_t num_tet_blocks = 0;
    {
        std::size_t bi = 0;
        for (const auto cb : rMesh.CellRange()) {
            const std::string type(cb.Type());
            if (cb.IsPolyhedron())
                throw std::invalid_argument(
                    std::string(kTpPrefix) + "cell block '" + type +
                    "' is already polyhedral; tet_to_poly needs a tetrahedral mesh (use "
                    "poly_to_tet first)");
            if (!cb.IsRagged()) {
                const CellType ct = cell_type_from_name(type);
                if (ct == CellType::Tetra) {
                    is_tet[bi] = true;
                    ++num_tet_blocks;
                } else if (ct == CellType::Triangle) {
                    is_tri[bi] = true;
                } else if (cell_type_dimension(ct) == 3) {
                    throw std::invalid_argument(
                        std::string(kTpPrefix) + "cell block '" + type + "' is not linear tetra; " +
                        (ct == CellType::Tetra10
                             ? "convert_cells(mode='linearize') it first"
                             : "decompose it with poly_to_tet or convert_cells(mode='simplexify') "
                               "first"));
                }
            }
            ++bi;
        }
    }
    if (num_tet_blocks == 0)
        throw std::invalid_argument(std::string(kTpPrefix) + "the mesh has no tetra block");
    if (rMesh.PointDim() != 3)
        throw std::invalid_argument(std::string(kTpPrefix) + "points must be 3D, got " +
                                    std::to_string(rMesh.PointDim()) + "D");

    std::vector<std::size_t> tet_blocks;
    for (std::size_t b = 0; b < nblocks_in; ++b)
        if (is_tet[b])
            tet_blocks.push_back(b);
    const detail::GlobalFaces gf = detail::build_global_faces(rMesh, tet_blocks);
    if (gf.mNumNonManifold > 0)
        throw std::invalid_argument(std::string(kTpPrefix) + "mesh contains " +
                                    std::to_string(gf.mNumNonManifold) +
                                    " face(s) shared by three or more tetrahedra (non-manifold)");
    if (gf.mNumUnorientable > 0)
        throw std::invalid_argument(std::string(kTpPrefix) + "mesh contains " +
                                    std::to_string(gf.mNumUnorientable) +
                                    " degenerate tetrahedra (repeated nodes)");

    const std::size_t num_points = rMesh.NumPoints();
    const std::size_t nc = gf.NumCells();
    const std::size_t nf = gf.NumFaces();
    const std::vector<std::int64_t> in_bases = detail::block_bases(rMesh);

    TpContext ctx;
    ctx.mpFaces = &gf;
    ctx.mNumPoints = num_points;
    ctx.mCosFeature = std::cos(rOptions.mFeatureAngleDeg * kTpPi / 180.0) - 1e-12;
    const detail::DoubleView xyz(rMesh.Points());
    ctx.mpXyz = xyz.Data();

    // Each tet's nodes in connectivity order and its four faces as rings in
    // that order's orientation: (b, c, d), (a, d, c), (a, b, d), (a, c, b),
    // outward when (a, b, c, d) is positively oriented. build_global_faces
    // repairs each cell's winding by its own geometry, which leaves an
    // inverted tet disagreeing with its neighbours; the dual's fans need one
    // orientation consistent across shared faces instead, chosen below.
    ctx.mBlockOf.resize(nc);
    ctx.mTetNodes.resize(4 * nc);
    ctx.mTetFaceId.assign(4 * nc, -1);
    ctx.mTetRing.resize(4 * nc);
    std::vector<std::int64_t> tet_row(nc);
    std::vector<double> signed_volume(nc);
    parallel_for(nc, [&](std::size_t c) {
        const auto [blk, row] = detail::global_to_block_row(in_bases, gf.mCellToGlobal[c]);
        ctx.mBlockOf[c] = blk;
        tet_row[c] = row;
        if (gf.NumCellFaces(c) != 4)
            throw std::logic_error(std::string(kTpPrefix) + "a tetrahedron without four faces");
        const NDArray& conn = rMesh.Cells(blk).Conn();
        std::int64_t* n = ctx.mTetNodes.data() + 4 * c;
        for (std::size_t k = 0; k < 4; ++k)
            n[k] = detail::read_int(conn, 4 * static_cast<std::size_t>(row) + k);
        if (n[0] == n[1] || n[0] == n[2] || n[0] == n[3] || n[1] == n[2] || n[1] == n[3] ||
            n[2] == n[3])
            throw std::invalid_argument(std::string(kTpPrefix) +
                                        "a tetrahedron repeats a node (degenerate)");
        const std::array<std::array<std::int64_t, 3>, 4> rings{
            {{n[1], n[2], n[3]}, {n[0], n[3], n[2]}, {n[0], n[1], n[3]}, {n[0], n[2], n[1]}}};
        const std::int64_t* row_faces = gf.CellFaces(c);
        for (std::size_t k = 0; k < 4; ++k) {
            ctx.mTetRing[4 * c + k] = rings[k];
            std::array<std::int64_t, 3> key = rings[k];
            std::sort(key.begin(), key.end());
            for (std::size_t j = 0; j < 4; ++j) {
                const std::int64_t sf = row_faces[j];
                const std::size_t f = static_cast<std::size_t>((sf > 0 ? sf : -sf) - 1);
                std::array<std::int64_t, 3> fk{gf.Face(f)[0], gf.Face(f)[1], gf.Face(f)[2]};
                std::sort(fk.begin(), fk.end());
                if (fk == key)
                    ctx.mTetFaceId[4 * c + k] = static_cast<std::int64_t>(f);
            }
            if (ctx.mTetFaceId[4 * c + k] < 0)
                throw std::logic_error(std::string(kTpPrefix) + "tetrahedron faces out of step");
        }
        signed_volume[c] = tp_signed_volume(xyz.Data(), n);
    });

    // One consistent orientation: a face shared by two tets must be traversed
    // in opposite directions by them. Propagate it across shared faces from
    // the lowest tet of each connected component, then flip a component whose
    // total signed volume comes out negative.
    {
        std::vector<std::int8_t> sign(nc, 0);
        std::vector<std::size_t> queue;
        std::vector<std::size_t> members;
        auto local_face = [&](std::size_t t, std::int64_t f) {
            for (int k = 0; k < 4; ++k)
                if (ctx.mTetFaceId[4 * t + static_cast<std::size_t>(k)] == f)
                    return k;
            return -1;
        };
        for (std::size_t seed = 0; seed < nc; ++seed) {
            if (sign[seed] != 0)
                continue;
            sign[seed] = 1;
            queue.assign(1, seed);
            members.clear();
            double volume = 0.0;
            for (std::size_t head = 0; head < queue.size(); ++head) {
                const std::size_t t = queue[head];
                members.push_back(t);
                volume += sign[t] * signed_volume[t];
                for (std::size_t k = 0; k < 4; ++k) {
                    const std::int64_t f = ctx.mTetFaceId[4 * t + k];
                    const std::int64_t other = ctx.Across(f, static_cast<std::int64_t>(t));
                    if (other < 0)
                        continue;
                    const auto u = static_cast<std::size_t>(other);
                    const int ku = local_face(u, f);
                    const auto& rt = ctx.mTetRing[4 * t + k];
                    const auto& ru = ctx.mTetRing[4 * u + static_cast<std::size_t>(ku)];
                    // Same cyclic direction as stored -> the signs must differ.
                    const bool same = tp_has_directed(ru, rt[0], rt[1]);
                    const std::int8_t want = static_cast<std::int8_t>(same ? -sign[t] : sign[t]);
                    if (sign[u] == 0) {
                        sign[u] = want;
                        queue.push_back(u);
                    } else if (sign[u] != want) {
                        throw std::invalid_argument(
                            std::string(kTpPrefix) +
                            "the tetrahedra cannot be oriented consistently (a face is shared "
                            "inconsistently)");
                    }
                }
            }
            if (volume < 0.0)
                for (const std::size_t t : members)
                    sign[t] = static_cast<std::int8_t>(-sign[t]);
        }
        parallel_for(nc, [&](std::size_t t) {
            if (sign[t] < 0)
                for (std::size_t k = 0; k < 4; ++k)
                    std::swap(ctx.mTetRing[4 * t + k][1], ctx.mTetRing[4 * t + k][2]);
        });
    }

    // Vertex -> tets, ascending (a counting pass over ascending tets).
    ctx.mVertTetStart.assign(num_points + 1, 0);
    for (std::size_t i = 0; i < 4 * nc; ++i)
        ++ctx.mVertTetStart[static_cast<std::size_t>(ctx.mTetNodes[i]) + 1];
    for (std::size_t p = 0; p < num_points; ++p)
        ctx.mVertTetStart[p + 1] += ctx.mVertTetStart[p];
    ctx.mVertTets.resize(4 * nc);
    {
        std::vector<std::int64_t> fill(ctx.mVertTetStart.begin(), ctx.mVertTetStart.end() - 1);
        for (std::size_t c = 0; c < nc; ++c)
            for (int k = 0; k < 4; ++k) {
                const auto p = static_cast<std::size_t>(ctx.mTetNodes[4 * c + k]);
                ctx.mVertTets[static_cast<std::size_t>(fill[p]++)] = static_cast<std::int64_t>(c);
            }
    }

    // Input triangle blocks tag the faces they cover.
    ctx.mFaceTag.assign(nf, -1);
    std::vector<std::int64_t> face_tri_row(nf, -1);
    {
        const detail::FaceLookup lookup(gf);
        std::size_t dropped = 0;
        std::size_t bi = 0;
        for (const auto cb : rMesh.CellRange()) {
            if (is_tri[bi]) {
                const NDArray& conn = cb.Conn();
                const detail::Int64View cv(conn);
                for (std::size_t r = 0; r < cb.NumCells(); ++r) {
                    const std::int64_t ids[3] = {cv[3 * r], cv[3 * r + 1], cv[3 * r + 2]};
                    const std::int64_t f = lookup.Find(ids, 3);
                    if (f < 0 || ctx.Interior(f)) {
                        ++dropped;
                        continue;
                    }
                    if (ctx.mFaceTag[static_cast<std::size_t>(f)] < 0) {
                        ctx.mFaceTag[static_cast<std::size_t>(f)] = static_cast<std::int64_t>(bi);
                        face_tri_row[static_cast<std::size_t>(f)] = static_cast<std::int64_t>(r);
                    }
                }
            }
            ++bi;
        }
        if (dropped > 0)
            log::warn("{}{} triangle(s) not on the boundary of a tetra block dropped", kTpPrefix,
                      dropped);
    }

    // Edges of every boundary or interface face get a midpoint.
    for (std::size_t f = 0; f < nf; ++f) {
        if (ctx.Interior(static_cast<std::int64_t>(f)))
            continue;
        const std::int64_t* n = gf.Face(f);
        for (int k = 0; k < 3; ++k)
            ctx.mMidEdges.push_back(tp_edge(n[k], n[(k + 1) % 3]));
    }
    std::sort(ctx.mMidEdges.begin(), ctx.mMidEdges.end());
    ctx.mMidEdges.erase(std::unique(ctx.mMidEdges.begin(), ctx.mMidEdges.end()),
                        ctx.mMidEdges.end());

    // --- the dual cells, vertex by vertex ------------------------------------
    std::vector<TpVertexOut> per_vertex(num_points);
    parallel_for(num_points, [&](std::size_t p) {
        per_vertex[p] = tp_build_vertex(ctx, static_cast<std::int64_t>(p));
    });

    std::vector<double> tet_volume(nc);
    parallel_for(nc, [&](std::size_t c) { tet_volume[c] = tp_tet_volume(ctx, c); });

    // --- stage the output blocks in input order ------------------------------
    std::vector<TpOutBlock> staged;
    std::vector<std::int64_t> out_of_in(nblocks_in, -1);
    std::int64_t num_boundary_faces = 0;
    std::int64_t num_split = 0;
    for (std::size_t p = 0; p < num_points; ++p)
        if (per_vertex[p].mCells.size() > 1)
            ++num_split;
    {
        std::size_t bi = 0;
        for (const auto cb : rMesh.CellRange()) {
            if (!is_tet[bi] && !is_tri[bi]) {
                log::warn(
                    "{}cell block '{}' is neither tetra nor a boundary triangle block; "
                    "dropped",
                    kTpPrefix, cb.Type());
                ++bi;
                continue;
            }
            TpOutBlock ob;
            ob.mInBlock = bi;
            ob.mIsPolygon = is_tri[bi];
            for (std::size_t p = 0; p < num_points; ++p) {
                const TpVertexOut& vo = per_vertex[p];
                if (!ob.mIsPolygon) {
                    for (const TpCell& cell : vo.mCells) {
                        if (cell.mBlock != bi)
                            continue;
                        const auto face_base = static_cast<std::int64_t>(ob.mFlat.size());
                        ob.mFlat.insert(ob.mFlat.end(), cell.mNodes.begin(), cell.mNodes.end());
                        for (std::size_t k = 1; k < cell.mFaceStart.size(); ++k)
                            ob.mFaceStart.push_back(face_base + cell.mFaceStart[k]);
                        ob.mCellStart.push_back(static_cast<std::int64_t>(ob.mFaceStart.size()) -
                                                1);
                        ob.mVertex.push_back(static_cast<std::int64_t>(p));
                        for (const std::int64_t t : cell.mTets) {
                            ob.mMemberRow.push_back(tet_row[static_cast<std::size_t>(t)]);
                            ob.mMemberWeight.push_back(tet_volume[static_cast<std::size_t>(t)]);
                        }
                        ob.mMemberStart.push_back(static_cast<std::int64_t>(ob.mMemberRow.size()));
                        num_boundary_faces += cell.mNumBoundaryFaces;
                    }
                } else {
                    for (const TpPolygon& poly : vo.mPolygons) {
                        if (poly.mTag != bi)
                            continue;
                        ob.mFlat.insert(ob.mFlat.end(), poly.mRing.begin(), poly.mRing.end());
                        ob.mFaceStart.push_back(static_cast<std::int64_t>(ob.mFlat.size()));
                        ob.mVertex.push_back(static_cast<std::int64_t>(p));
                        for (const std::int64_t f : poly.mFaces) {
                            ob.mMemberRow.push_back(face_tri_row[static_cast<std::size_t>(f)]);
                            ob.mMemberWeight.push_back(
                                tp_face_area(ctx, static_cast<std::size_t>(f)));
                        }
                        ob.mMemberStart.push_back(static_cast<std::int64_t>(ob.mMemberRow.size()));
                    }
                }
            }
            out_of_in[bi] = static_cast<std::int64_t>(staged.size());
            staged.push_back(std::move(ob));
            ++bi;
        }
    }
    std::vector<TpVertexOut>().swap(per_vertex);

    // --- points: originals, then edge midpoints, face and tet centroids ------
    const std::size_t n_mid = ctx.mMidEdges.size();
    const std::size_t num_new = n_mid + nf + nc;
    auto sources = [&](std::size_t i, std::int64_t* pOut) -> std::size_t {
        if (i < n_mid) {
            pOut[0] = ctx.mMidEdges[i].first;
            pOut[1] = ctx.mMidEdges[i].second;
            return 2;
        }
        i -= n_mid;
        if (i < nf) {
            const std::int64_t* n = gf.Face(i);
            std::copy(n, n + 3, pOut);
            return 3;
        }
        i -= nf;
        std::copy(ctx.mTetNodes.data() + 4 * i, ctx.mTetNodes.data() + 4 * i + 4, pOut);
        return 4;
    };
    // Append the source-node mean of every new point to a per-point array.
    auto extend = [&](const NDArray& rA) {
        const std::size_t ncomp = num_points == 0 ? 0 : rA.Size() / num_points;
        std::vector<std::size_t> shape = rA.Shape();
        shape[0] = num_points + num_new;
        NDArray b = NDArray::Uninit(rA.Dtype(), std::move(shape));
        std::memcpy(b.Data(), rA.Data(), rA.Nbytes());
        const detail::DoubleView av(rA);
        detail::DoubleSink sink(b, num_points * ncomp, num_new * ncomp);
        double* const vals = sink.Data();
        parallel_for(num_new, [&](std::size_t i) {
            std::int64_t src[4];
            const std::size_t count = sources(i, src);
            for (std::size_t k = 0; k < ncomp; ++k) {
                double sum = 0.0;
                for (std::size_t j = 0; j < count; ++j)
                    sum += av[static_cast<std::size_t>(src[j]) * ncomp + k];
                vals[i * ncomp + k] = sum / static_cast<double>(count);
            }
        });
        sink.Commit();
        return b;
    };

    Mesh out;
    out.AssignPoints(extend(rMesh.Points()));
    for (TpOutBlock& ob : staged) {
        if (ob.mIsPolygon)
            out.AddPolygonBlock("polygon", ob.mFlat, ob.mFaceStart);
        else
            out.AddPolyhedronBlock("polyhedron", ob.mFlat, ob.mFaceStart, ob.mCellStart);
        std::vector<std::int64_t>().swap(ob.mFlat);
        std::vector<std::int64_t>().swap(ob.mFaceStart);
        std::vector<std::int64_t>().swap(ob.mCellStart);
    }

    for (const std::string& name : rMesh.PointDataNames()) {
        const NDArray& a = rMesh.PointData(name);
        if (detail::rows(a) != num_points) {
            out.AddPointData(name, detail::data_owned_copy(a));
            continue;
        }
        out.AddPointData(name, extend(a));
    }
    for (const std::string& name : rMesh.FieldDataNames())
        out.AddFieldData(name, detail::data_owned_copy(rMesh.FieldData(name)));

    // --- cell_data -----------------------------------------------------------
    std::vector<std::string> cell_names = rMesh.CellDataNames();
    for (const std::string& name : cell_names) {
        bool ok = rMesh.CellDataNumBlocks(name) == nblocks_in;
        for (std::size_t b = 0; ok && b < nblocks_in; ++b)
            if (out_of_in[b] >= 0 &&
                detail::rows(rMesh.CellData(name, b)) != rMesh.Cells(b).NumCells())
                ok = false;
        if (!ok) {
            log::warn("{}cell_data '{}' does not have one row per cell of every block; dropped",
                      kTpPrefix, name);
            continue;
        }
        std::vector<NDArray> blocks;
        blocks.reserve(staged.size());
        for (const TpOutBlock& ob : staged)
            blocks.push_back(tp_reduce_rows(rMesh.CellData(name, ob.mInBlock), ob.mMemberStart,
                                            ob.mMemberRow, ob.mMemberWeight));
        if (!blocks.empty())
            out.AddCellData(name, std::move(blocks));
    }
    if (rOptions.mPointDataToCellData) {
        for (const std::string& name : rMesh.PointDataNames()) {
            const NDArray& a = rMesh.PointData(name);
            if (detail::rows(a) != num_points)
                continue;
            if (std::find(cell_names.begin(), cell_names.end(), name) != cell_names.end()) {
                log::warn("{}point_data '{}' not copied to cell_data: the name is taken", kTpPrefix,
                          name);
                continue;
            }
            std::vector<NDArray> blocks;
            for (const TpOutBlock& ob : staged)
                blocks.push_back(tp_gather_rows(a, ob.mVertex));
            if (!blocks.empty())
                out.AddCellData(name, std::move(blocks));
        }
    }
    if (rOptions.mRecordVertexIds) {
        std::vector<NDArray> blocks;
        for (const TpOutBlock& ob : staged) {
            NDArray ids = NDArray::Uninit(DType::Int64, {ob.NumCells()});
            if (ob.NumCells() != 0)
                std::memcpy(ids.Data(), ob.mVertex.data(), ob.NumCells() * sizeof(std::int64_t));
            blocks.push_back(std::move(ids));
        }
        if (!blocks.empty())
            out.AddCellData("dual:vertex", std::move(blocks));
    }

    TetToPolyResult res;
    {
        std::size_t total = 0;
        for (const TpOutBlock& ob : staged)
            total += ob.NumCells();
        res.mVertexMap = NDArray::Uninit(DType::Int64, {total});
        std::int64_t* vm = res.mVertexMap.As<std::int64_t>();
        for (const TpOutBlock& ob : staged)
            for (const std::int64_t v : ob.mVertex)
                *vm++ = v;
    }
    res.mMesh = std::move(out);
    res.mNumNewPoints = static_cast<std::int64_t>(num_new);
    res.mNumBoundaryFaces = num_boundary_faces;
    res.mNumSplitVertices = num_split;

    // Point regions keep their ids; a cell or side region has no counterpart.
    NDArray none =
        NDArray::Uninit(DType::Int64, {static_cast<std::size_t>(detail::total_cells(in_bases))});
    std::fill(none.As<std::int64_t>(), none.As<std::int64_t>() + none.Size(), std::int64_t{-1});
    detail::RegionRemap rmap;
    rmap.mCellMapKind = detail::CellMapKind::Global;
    rmap.pGlobalCellMap = &none;
    rmap.mOpName = "tet_to_poly";
    detail::remap_regions(rMesh, res.mMesh, rmap);
    return res;
}

// --- poly_to_tet -------------------------------------------------------------

PolyToTetResult poly_to_tet(const Mesh& rMesh, const PolyToTetOptions& rOptions) {
    detail::SimplexifyCounts counts;
    ConvertCellsResult cc =
        detail::simplexify_ex(rMesh, rOptions.mRecordParentIds, /*ExactTets=*/true, &counts);
    PolyToTetResult res;
    res.mMesh = std::move(cc.mMesh);
    res.mPointMap = std::move(cc.mPointMap);
    res.mCellMaps = std::move(cc.mCellMaps);
    res.mNumExactTets = counts.mNumExactTets;
    res.mNumFannedCells = counts.mNumFannedCells;
    return res;
}

}  // namespace meshioplusplus

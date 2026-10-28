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
// The all-hexahedral dual of a paired, strongly balanced octree. See
// detail/hex_dual_schemes.hpp.
//
// The scheme detection, marking and orientation are ported from CinoLib's
// include/cinolib/hex_transition_install.cpp and hex_transition_orient.cpp,
// and the dual from include/cinolib/dual_mesh.cpp. CinoLib carries this
// notice:
//
//   Copyright(C) 2016: Marco Livesu
//   Author(s): Luca Pitzalis, Marco Livesu
//
//   Permission is hereby granted, free of charge, to any person obtaining a
//   copy of this software and associated documentation files (the "Software"),
//   to deal in the Software without restriction, including without limitation
//   the rights to use, copy, modify, merge, publish, distribute, sublicense,
//   and/or sell copies of the Software, and to permit persons to whom the
//   Software is furnished to do so, subject to the following conditions:
//
//   The above copyright notice and this permission notice shall be included in
//   all copies or substantial portions of the Software.
//
//   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//   FITNESS FOR A PARTICULAR PURPOSE AND NON INFRINGEMENT. IN NO EVENT SHALL THE
//   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
//   IN THE SOFTWARE.

// System includes
#include <algorithm>
#include <array>
#include <cstdint>
#include <deque>
#include <initializer_list>
#include <iterator>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "hex_dual_schemes.hpp"
#include "meshioplusplus/exceptions.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

constexpr const char* kHdsPrefix = "meshio++: octree_mesh: graded hexahedra: ";

using HdsPoint = std::array<std::int64_t, 3>;

struct HdsPointHash {
    std::size_t operator()(const HdsPoint& rP) const {
        std::uint64_t h = static_cast<std::uint64_t>(rP[0]) * 0x9E3779B97F4A7C15ULL;
        h ^= static_cast<std::uint64_t>(rP[1]) + 0x7F4A7C159E3779B9ULL + (h << 6) + (h >> 2);
        h ^= static_cast<std::uint64_t>(rP[2]) + 0x94D049BB133111EBULL + (h << 6) + (h >> 2);
        return static_cast<std::size_t>(h);
    }
};

std::string hds_str(const HdsPoint& rP) {
    return "(" + std::to_string(rP[0]) + ", " + std::to_string(rP[1]) + ", " +
           std::to_string(rP[2]) + ")";
}

[[noreturn]] void hds_fail(const std::string& rWhy) {
    throw Unsupported(std::string(kHdsPrefix) + rWhy);
}

// Directions, with CinoLib's values (symbols.h), so the ported case tables read
// the same and sort in the same order.
enum : int {
    kHdsPlusX = 10,
    kHdsPlusY = 11,
    kHdsPlusZ = 12,
    kHdsMinusX = 13,
    kHdsMinusY = 14,
    kHdsMinusZ = 15,
};

enum class HdsKind { FlatS, ConcS, CornS, ConvS };

struct HdsInfo {
    HexTransition mType = HexTransition::Flat;
    HdsKind mKind = HdsKind::FlatS;
    std::int64_t mScale = 0;  // the cell's edge, lattice units
    std::vector<HdsPoint> mTVerts;
    std::map<int, bool> mCuts;
    std::vector<int> mOrientations;

    bool Cut(int Dir) const {
        const auto it = mCuts.find(Dir);
        return it != mCuts.end() && it->second;
    }
};

// ---- the primal grid: the leaves as hexahedra with hanging nodes ----------

// Hexahedron corners by bits (x, y, z), its faces and its edges.
constexpr int kHdsFaces[6][4] = {{0, 2, 6, 4}, {1, 3, 7, 5}, {0, 1, 5, 4},
                                 {2, 3, 7, 6}, {0, 1, 3, 2}, {4, 5, 7, 6}};
constexpr int kHdsEdges[12][2] = {{0, 1}, {2, 3}, {4, 5}, {6, 7}, {0, 2}, {1, 3},
                                  {4, 6}, {5, 7}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};

struct HdsGrid {
    std::vector<HdsPoint> mVerts;
    std::unordered_map<HdsPoint, std::int32_t, HdsPointHash> mVertId;
    std::vector<std::array<std::int32_t, 8>> mPolyVerts;
    std::vector<std::array<std::int32_t, 6>> mPolyFaces;
    std::vector<std::array<std::int32_t, 12>> mPolyEdges;
    std::vector<std::int32_t> mLabel;  // the leaf's level
    std::vector<std::array<std::int32_t, 4>> mFaceVerts;
    std::vector<std::vector<std::int32_t>> mFacePolys;
    std::vector<std::array<std::int32_t, 2>> mEdgeVerts;
    std::vector<std::vector<std::int32_t>> mEdgePolys;
    std::unordered_map<std::uint64_t, std::int32_t> mEdgeId;
    std::vector<std::vector<std::int32_t>> mV2P, mV2E, mV2V, mP2P;

    static std::uint64_t PairKey(std::int32_t A, std::int32_t B) {
        if (A > B)
            std::swap(A, B);
        return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(A)) << 32) |
               static_cast<std::uint32_t>(B);
    }
    std::int32_t EdgeId(std::int32_t A, std::int32_t B) const {
        const auto it = mEdgeId.find(PairKey(A, B));
        return it == mEdgeId.end() ? -1 : it->second;
    }
    bool PolyHasVert(std::int32_t P, std::int32_t V) const {
        const auto& r = mPolyVerts[static_cast<std::size_t>(P)];
        return std::find(r.begin(), r.end(), V) != r.end();
    }
    bool PolyHasEdge(std::int32_t P, std::int32_t E) const {
        const auto& r = mPolyEdges[static_cast<std::size_t>(P)];
        return std::find(r.begin(), r.end(), E) != r.end();
    }
    const HdsPoint& Vert(std::int32_t V) const { return mVerts[static_cast<std::size_t>(V)]; }
    std::int32_t Label(std::int32_t P) const { return mLabel[static_cast<std::size_t>(P)]; }
    std::int64_t EdgeLength(std::int32_t E) const {
        const auto& e = mEdgeVerts[static_cast<std::size_t>(E)];
        const HdsPoint &a = Vert(e[0]), &b = Vert(e[1]);
        return std::abs(a[0] - b[0]) + std::abs(a[1] - b[1]) + std::abs(a[2] - b[2]);
    }
    std::int64_t PolyScale(std::int32_t P) const {
        return EdgeLength(mPolyEdges[static_cast<std::size_t>(P)][0]);
    }
    HdsPoint EdgeMid(std::int32_t E) const {
        const auto& e = mEdgeVerts[static_cast<std::size_t>(E)];
        const HdsPoint &a = Vert(e[0]), &b = Vert(e[1]);
        return {(a[0] + b[0]) / 2, (a[1] + b[1]) / 2, (a[2] + b[2]) / 2};
    }
    // The axis direction of edge `E` leaving vertex `V`.
    int Orientation(std::int32_t V, std::int32_t E) const {
        const auto& e = mEdgeVerts[static_cast<std::size_t>(E)];
        const HdsPoint& a = Vert(V);
        const HdsPoint& b = Vert(e[0] == V ? e[1] : e[0]);
        for (int k = 0; k < 3; ++k) {
            const std::int64_t d = b[static_cast<std::size_t>(k)] - a[static_cast<std::size_t>(k)];
            if (d > 0)
                return kHdsPlusX + k;
            if (d < 0)
                return kHdsMinusX + k;
        }
        hds_fail("a degenerate edge at " + hds_str(a));
    }
    int MinRef(std::int32_t V) const {
        int r = 1 << 30;
        for (std::int32_t p : mV2P[static_cast<std::size_t>(V)])
            r = std::min(r, Label(p));
        return r;
    }
    int MaxRef(std::int32_t V) const {
        int r = -1;
        for (std::int32_t p : mV2P[static_cast<std::size_t>(V)])
            r = std::max(r, Label(p));
        return r;
    }
    std::int32_t SharedFace(std::int32_t P, std::int32_t Q) const {
        for (std::int32_t f : mPolyFaces[static_cast<std::size_t>(P)])
            for (std::int32_t g : mPolyFaces[static_cast<std::size_t>(Q)])
                if (f == g)
                    return f;
        return -1;
    }
    // The vertex at `rAt` among the neighbours of `V`, or -1.
    std::int32_t NeighbourAt(std::int32_t V, const HdsPoint& rAt) const {
        for (std::int32_t w : mV2V[static_cast<std::size_t>(V)])
            if (Vert(w) == rAt)
                return w;
        return -1;
    }
};

std::int32_t hds_vertex(HdsGrid& rG, const HdsPoint& rP) {
    const auto it = rG.mVertId.find(rP);
    if (it != rG.mVertId.end())
        return it->second;
    const auto id = static_cast<std::int32_t>(rG.mVerts.size());
    rG.mVerts.push_back(rP);
    rG.mVertId.emplace(rP, id);
    return id;
}

HdsGrid hds_grid(const LinearOctree& rTree) {
    HdsGrid g;
    const std::size_t n = rTree.Size();
    g.mPolyVerts.resize(n);
    g.mPolyFaces.resize(n);
    g.mPolyEdges.resize(n);
    g.mLabel.resize(n);
    std::map<std::array<std::int32_t, 4>, std::int32_t> face_id;
    for (std::size_t i = 0; i < n; ++i) {
        const auto lo = rTree.Lo(i);
        const std::int64_t s = rTree.Span(i) * kHexDualUnit;
        g.mLabel[i] = rTree.mLevel[i];
        for (int c = 0; c < 8; ++c)
            g.mPolyVerts[i][static_cast<std::size_t>(c)] = hds_vertex(
                g, {lo[0] * kHexDualUnit + (c & 1) * s, lo[1] * kHexDualUnit + ((c >> 1) & 1) * s,
                    lo[2] * kHexDualUnit + ((c >> 2) & 1) * s});
        const auto& pv = g.mPolyVerts[i];
        for (int f = 0; f < 6; ++f) {
            std::array<std::int32_t, 4> q{};
            for (int k = 0; k < 4; ++k)
                q[static_cast<std::size_t>(k)] = pv[static_cast<std::size_t>(kHdsFaces[f][k])];
            std::array<std::int32_t, 4> key = q;
            std::sort(key.begin(), key.end());
            auto it = face_id.find(key);
            if (it == face_id.end()) {
                it = face_id.emplace(key, static_cast<std::int32_t>(g.mFaceVerts.size())).first;
                g.mFaceVerts.push_back(q);
                g.mFacePolys.emplace_back();
            }
            g.mPolyFaces[i][static_cast<std::size_t>(f)] = it->second;
            g.mFacePolys[static_cast<std::size_t>(it->second)].push_back(
                static_cast<std::int32_t>(i));
        }
        for (int e = 0; e < 12; ++e) {
            const std::int32_t a = pv[static_cast<std::size_t>(kHdsEdges[e][0])];
            const std::int32_t b = pv[static_cast<std::size_t>(kHdsEdges[e][1])];
            auto it = g.mEdgeId.find(HdsGrid::PairKey(a, b));
            if (it == g.mEdgeId.end()) {
                it = g.mEdgeId
                         .emplace(HdsGrid::PairKey(a, b),
                                  static_cast<std::int32_t>(g.mEdgeVerts.size()))
                         .first;
                g.mEdgeVerts.push_back({a, b});
                g.mEdgePolys.emplace_back();
            }
            g.mPolyEdges[i][static_cast<std::size_t>(e)] = it->second;
            g.mEdgePolys[static_cast<std::size_t>(it->second)].push_back(
                static_cast<std::int32_t>(i));
        }
    }
    const std::size_t nv = g.mVerts.size();
    g.mV2P.resize(nv);
    g.mV2E.resize(nv);
    g.mV2V.resize(nv);
    for (std::size_t i = 0; i < n; ++i)
        for (std::int32_t v : g.mPolyVerts[i])
            g.mV2P[static_cast<std::size_t>(v)].push_back(static_cast<std::int32_t>(i));
    for (std::size_t e = 0; e < g.mEdgeVerts.size(); ++e) {
        const auto& ev = g.mEdgeVerts[e];
        g.mV2E[static_cast<std::size_t>(ev[0])].push_back(static_cast<std::int32_t>(e));
        g.mV2E[static_cast<std::size_t>(ev[1])].push_back(static_cast<std::int32_t>(e));
        g.mV2V[static_cast<std::size_t>(ev[0])].push_back(ev[1]);
        g.mV2V[static_cast<std::size_t>(ev[1])].push_back(ev[0]);
    }
    g.mP2P.resize(n);
    for (const auto& fp : g.mFacePolys)
        if (fp.size() == 2) {
            g.mP2P[static_cast<std::size_t>(fp[0])].push_back(fp[1]);
            g.mP2P[static_cast<std::size_t>(fp[1])].push_back(fp[0]);
        }
    for (auto& r : g.mP2P)
        std::sort(r.begin(), r.end());
    return g;
}

// ---- transition vertices --------------------------------------------------

// The centre of each face of a parent whose eight children are leaves, when
// the leaves across that face at its centre are one level finer than the
// children; its direction points into the parent. Checked against CinoLib's
// rule (the longest edge shared by four cells), which must agree.
void hds_transition_vertices(const LinearOctree& rTree, const HdsGrid& rG,
                             std::vector<std::uint8_t>& rIsT, std::vector<int>& rDir) {
    rIsT.assign(rG.mVerts.size(), 0);
    rDir.assign(rG.mVerts.size(), 0);
    for (std::size_t i = 0; i < rTree.Size(); ++i) {
        const auto lo = rTree.Lo(i);
        const std::int64_t s = rTree.Span(i);
        const int level = rTree.mLevel[i];
        if (level == 0 || lo[0] % (2 * s) != 0 || lo[1] % (2 * s) != 0 || lo[2] % (2 * s) != 0)
            continue;  // not a parent's first child
        for (int axis = 0; axis < 3; ++axis)
            for (int side = 0; side < 2; ++side) {
                // The four finest cells across the face, around its centre.
                const auto ax = static_cast<std::size_t>(axis);
                const std::size_t u = (ax + 1) % 3, w = (ax + 2) % 3;
                std::array<std::int64_t, 3> centre{lo[0] + s, lo[1] + s, lo[2] + s};
                centre[ax] = side == 0 ? lo[ax] : lo[ax] + 2 * s;
                int finer = 0, same = 0;
                for (int q = 0; q < 4; ++q) {
                    std::array<std::int64_t, 3> c = centre;
                    c[ax] = side == 0 ? centre[ax] - 1 : centre[ax];
                    c[u] += (q & 1) ? 0 : -1;
                    c[w] += (q & 2) ? 0 : -1;
                    if (!rTree.Contains(c[0], c[1], c[2]))
                        continue;
                    const std::int64_t j = rTree.Locate(c[0], c[1], c[2]);
                    const int lj = rTree.mLevel[static_cast<std::size_t>(j)];
                    if (lj == level + 1)
                        ++finer;
                    else
                        ++same;
                }
                if (finer == 0)
                    continue;
                const HdsPoint at{centre[0] * kHexDualUnit, centre[1] * kHexDualUnit,
                                  centre[2] * kHexDualUnit};
                if (finer != 4 || same != 0)
                    hds_fail("the leaves across the face centre " + hds_str(at) +
                             " are not all one level finer; pair the tree first");
                const std::int32_t v = rG.mVertId.at(at);
                rIsT[static_cast<std::size_t>(v)] = 1;
                rDir[static_cast<std::size_t>(v)] =
                    side == 0 ? kHdsPlusX + axis : kHdsMinusX + axis;
            }
    }
    for (std::size_t v = 0; v < rIsT.size(); ++v) {
        if (!rIsT[v])
            continue;
        const auto vid = static_cast<std::int32_t>(v);
        if (rG.mV2P[v].size() != 8)
            hds_fail("the transition vertex " + hds_str(rG.Vert(vid)) + " has " +
                     std::to_string(rG.mV2P[v].size()) + " cells around it, not 8");
        std::int64_t length = 0;
        int dir = -1;
        for (std::int32_t e : rG.mV2E[v])
            if (rG.mEdgePolys[static_cast<std::size_t>(e)].size() == 4 &&
                rG.EdgeLength(e) > length) {
                length = rG.EdgeLength(e);
                dir = rG.Orientation(vid, e);
            }
        if (dir != rDir[v])
            hds_fail("the transition vertex " + hds_str(rG.Vert(vid)) +
                     " points into a cell its edges do not");
    }
}

// ---- marking (CinoLib's hex_transition_install) ---------------------------

using HdsSchemes = std::map<std::int32_t, HdsInfo>;

bool hds_is(const HdsSchemes& rMap, std::int32_t P, std::initializer_list<HexTransition> Types) {
    const auto it = rMap.find(P);
    if (it == rMap.end())
        return false;
    return std::find(Types.begin(), Types.end(), it->second.mType) != Types.end();
}

// Whether the vertex at the middle of edge `V`-`W` has cells of two levels
// around it: a weakly balanced corner, which the ported schemes do not cover.
bool hds_middle_is_cut(const HdsGrid& rG, std::int32_t From, std::int32_t V, std::int32_t W) {
    const std::int32_t e = rG.EdgeId(V, W);
    const HdsPoint mid = rG.EdgeMid(e);
    const std::int32_t m = rG.NeighbourAt(From, mid);
    if (m < 0)
        hds_fail("no vertex at the middle " + hds_str(mid) + " of a transition edge");
    std::set<int> refs;
    for (std::int32_t p : rG.mV2P[static_cast<std::size_t>(m)])
        refs.insert(rG.Label(p));
    return refs.size() > 1;
}

void hds_mark_concave_vert(const HdsGrid& rG, const std::vector<std::int32_t>& rT,
                           const std::vector<int>& rDir, HdsSchemes& rMap) {
    const int min_ref = rG.MinRef(rT[0]);
    std::set<std::int32_t> cluster;
    for (std::int32_t v : rT)
        for (std::int32_t p : rG.mV2P[static_cast<std::size_t>(v)])
            if (rG.Label(p) == min_ref)
                cluster.insert(p);
    std::set<std::int32_t> flats = cluster;
    std::int32_t concave = -1, conc_vert = -1;
    std::vector<std::int32_t> adjs;
    for (std::int32_t p : cluster) {
        if (!(rG.PolyHasVert(p, rT[0]) && rG.PolyHasVert(p, rT[1]) && rG.PolyHasVert(p, rT[2])))
            continue;
        concave = p;
        const auto& pv = rG.mPolyVerts[static_cast<std::size_t>(p)];
        std::set<std::int32_t> rest(pv.begin(), pv.end());
        for (std::int32_t t : rT) {
            for (std::int32_t v : rG.mV2V[static_cast<std::size_t>(t)])
                if (rG.PolyHasVert(p, v))
                    rest.erase(v);
            rest.erase(t);
        }
        if (rest.size() != 1)
            hds_fail("a concave corner cell at " + hds_str(rG.Vert(rT[0])) + " has " +
                     std::to_string(rest.size()) + " candidate corners, not 1");
        conc_vert = *rest.begin();
        for (std::int32_t v : rG.mV2V[static_cast<std::size_t>(conc_vert)])
            if (rG.PolyHasVert(concave, v))
                adjs.push_back(v);
        if (adjs.size() != 3)
            hds_fail("a concave corner at " + hds_str(rG.Vert(conc_vert)) + " has " +
                     std::to_string(adjs.size()) + " edges in its cell, not 3");
        break;
    }
    if (concave < 0)
        hds_fail("no cell holds the three transition vertices of the concave corner at " +
                 hds_str(rG.Vert(rT[0])));
    std::map<std::int32_t, bool> v2cut;
    if (!rMap.count(concave) ||
        hds_is(rMap, concave,
               {HexTransition::Flat, HexTransition::FlatConvex, HexTransition::VertSide})) {
        HdsInfo info;
        info.mType = HexTransition::VertCenter;
        info.mScale = rG.PolyScale(concave);
        for (std::int32_t v : adjs) {
            const bool cut = hds_middle_is_cut(rG, v, conc_vert, v);
            if (cut)
                hds_fail("the concave corner at " + hds_str(rG.Vert(conc_vert)) +
                         " is only weakly balanced");
            v2cut[v] = false;
            int ori = rG.Orientation(conc_vert, rG.EdgeId(conc_vert, v));
            ori = ori > 12 ? ori - 3 : ori;
            info.mCuts[ori] = false;
        }
        info.mKind = HdsKind::CornS;
        for (std::int32_t t : rT)
            info.mOrientations.push_back(rDir[static_cast<std::size_t>(t)]);
        info.mTVerts.push_back(rG.Vert(conc_vert));
        rMap[concave] = info;
    }
    flats.erase(concave);

    int lateral_y = kHdsMinusY;
    for (std::int32_t t : rT)
        if (rDir[static_cast<std::size_t>(t)] == kHdsPlusY)
            lateral_y = kHdsPlusY;

    for (std::int32_t p : rG.mP2P[static_cast<std::size_t>(concave)]) {
        const auto it = rMap.find(p);
        if (it == rMap.end() ||
            (it->second.mKind == HdsKind::CornS && it->second.mType != HexTransition::VertCenter)) {
            HdsInfo info;
            bool typed = false;
            for (const auto& kv : v2cut) {
                if (!rG.PolyHasVert(p, kv.first))
                    continue;
                for (std::int32_t v : rG.mV2V[static_cast<std::size_t>(kv.first)])
                    if (rG.PolyHasVert(p, v) && hds_middle_is_cut(rG, kv.first, kv.first, v))
                        hds_fail("the cell beside the concave corner at " +
                                 hds_str(rG.Vert(conc_vert)) + " is only weakly balanced");
                info.mType = HexTransition::VertSide;
                typed = true;
            }
            if (!typed)
                hds_fail("a cell beside the concave corner at " + hds_str(rG.Vert(conc_vert)) +
                         " touches none of its edges");
            info.mKind = HdsKind::CornS;
            info.mScale = rG.PolyScale(concave);
            const std::int32_t f = rG.SharedFace(p, concave);
            for (std::int32_t v : rG.mFaceVerts[static_cast<std::size_t>(f)]) {
                if (v == rT[0] || v == rT[1] || v == rT[2]) {
                    info.mOrientations.push_back(rDir[static_cast<std::size_t>(v)]);
                } else {
                    const int cur = rG.Label(p);
                    std::set<int> refs;
                    for (std::int32_t q : rG.mV2P[static_cast<std::size_t>(v)])
                        if (rG.Label(q) >= cur)
                            refs.insert(rG.Label(q));
                    if (refs.size() >= 2)
                        info.mTVerts.push_back(rG.Vert(v));
                }
            }
            info.mOrientations.push_back(lateral_y);
            rMap[p] = info;
        }
        flats.erase(p);
    }

    for (std::int32_t p : flats) {
        if (rMap.count(p))
            continue;
        HdsInfo info;
        info.mType = HexTransition::Flat;
        info.mKind = HdsKind::CornS;
        info.mScale = rG.PolyScale(concave);
        for (std::size_t k = 0; k < 3; ++k)
            if (rG.PolyHasVert(p, rT[k])) {
                info.mOrientations.push_back(rDir[static_cast<std::size_t>(rT[k])]);
                for (std::size_t j = 0; j < 3; ++j)
                    if (j != k)
                        info.mOrientations.push_back(rDir[static_cast<std::size_t>(rT[j])]);
                info.mTVerts.push_back(rG.Vert(rT[k]));
                break;
            }
        rMap[p] = info;
    }
}

void hds_mark_concave_edge(const HdsGrid& rG, const std::vector<std::int32_t>& rT,
                           const std::vector<int>& rDir, HdsSchemes& rMap) {
    std::vector<std::int32_t> a = rG.mV2V[static_cast<std::size_t>(rT[0])];
    std::vector<std::int32_t> b = rG.mV2V[static_cast<std::size_t>(rT[1])];
    std::sort(a.begin(), a.end());
    std::sort(b.begin(), b.end());
    std::vector<std::int32_t> common;
    std::set_intersection(a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(common));
    if (common.size() != 1 && common.size() != 2)
        hds_fail("the concave edge between " + hds_str(rG.Vert(rT[0])) + " and " +
                 hds_str(rG.Vert(rT[1])) + " has " + std::to_string(common.size()) +
                 " common neighbours");
    int ref = 0;
    for (std::int32_t p : rG.mV2P[static_cast<std::size_t>(rT[0])])
        if (rG.PolyHasVert(p, rT[1])) {
            ref = rG.Label(p);
            break;
        }
    std::int32_t conc = -1;
    for (std::int32_t v : common) {
        std::set<int> refs;
        for (std::int32_t p : rG.mV2P[static_cast<std::size_t>(v)])
            if (ref <= rG.Label(p))
                refs.insert(rG.Label(p));
        if (refs.size() == 2) {
            conc = v;
            break;
        }
        if (refs.size() == 3)
            hds_fail("the concave edge at " + hds_str(rG.Vert(v)) + " is only weakly balanced");
    }
    if (conc < 0)
        hds_fail("no concave edge vertex between " + hds_str(rG.Vert(rT[0])) + " and " +
                 hds_str(rG.Vert(rT[1])));

    std::vector<std::int32_t> concaves;
    const int min_ref = rG.MinRef(conc);
    for (std::int32_t p : rG.mV2P[static_cast<std::size_t>(conc)]) {
        if (!(rG.PolyHasVert(p, rT[0]) && rG.PolyHasVert(p, rT[1])) || rG.Label(p) != min_ref)
            continue;
        if (rMap.count(p) && !hds_is(rMap, p, {HexTransition::Flat, HexTransition::FlatConvex}))
            continue;
        int orient = -1;
        std::int32_t check = -1;
        for (std::int32_t e : rG.mV2E[static_cast<std::size_t>(conc)]) {
            const int o = rG.Orientation(conc, e);
            if (rG.PolyHasEdge(p, e) && o != rDir[static_cast<std::size_t>(rT[0])] &&
                o != rDir[static_cast<std::size_t>(rT[1])]) {
                orient = o;
                check = e;
                break;
            }
        }
        if (check < 0)
            hds_fail("the concave edge cell at " + hds_str(rG.Vert(conc)) +
                     " has no edge along it");
        const auto& ev = rG.mEdgeVerts[static_cast<std::size_t>(check)];
        if (hds_middle_is_cut(rG, conc, ev[0], ev[1]))
            hds_fail("the concave edge at " + hds_str(rG.Vert(conc)) + " is only weakly balanced");
        concaves.push_back(p);
        HdsInfo info;
        info.mType = HexTransition::Edge;
        info.mScale = rG.PolyScale(p);
        info.mTVerts.push_back(rG.Vert(conc));
        info.mOrientations.push_back(orient);
        info.mKind = HdsKind::ConcS;
        rMap[p] = info;
    }

    for (std::size_t i = 0; i < rT.size(); ++i)
        for (std::int32_t p : rG.mV2P[static_cast<std::size_t>(rT[i])]) {
            if (rG.PolyHasVert(p, conc) || rG.Label(p) != min_ref || rMap.count(p))
                continue;
            HdsInfo info;
            info.mType = HexTransition::Flat;
            for (std::int32_t q : rG.mP2P[static_cast<std::size_t>(p)]) {
                if (rG.PolyHasVert(q, rT[i]))
                    continue;
                if (hds_is(
                        rMap, q,
                        {HexTransition::Edge, HexTransition::VertSide, HexTransition::FlatConvex}))
                    info.mType = HexTransition::FlatConvex;
            }
            info.mScale = rG.PolyScale(p);
            info.mTVerts.push_back(rG.Vert(rT[i]));
            info.mOrientations.push_back(rDir[static_cast<std::size_t>(rT[i])]);
            info.mOrientations.push_back(rDir[static_cast<std::size_t>(rT[(i + 1) % rT.size()])]);
            info.mKind = HdsKind::ConcS;
            rMap[p] = info;
        }

    for (std::int32_t c : concaves)
        for (std::int32_t q : rG.mP2P[static_cast<std::size_t>(c)]) {
            const auto it = rMap.find(q);
            if (it == rMap.end() || rG.PolyHasVert(q, rT[0]) || rG.PolyHasVert(q, rT[1]))
                continue;
            if (it->second.mType == HexTransition::Flat)
                it->second.mType = HexTransition::FlatConvex;
        }
}

void hds_mark_convex(const HdsGrid& rG, const std::vector<std::int32_t>& rT,
                     const std::vector<int>& rDir, HdsSchemes& rMap) {
    const std::int32_t vert = rT.back();
    const int min_ref = rG.MinRef(vert);
    std::vector<std::int32_t> a = rG.mV2P[static_cast<std::size_t>(rT[0])];
    std::vector<std::int32_t> b = rG.mV2P[static_cast<std::size_t>(rT[1])];
    std::sort(a.begin(), a.end());
    std::sort(b.begin(), b.end());
    std::vector<std::int32_t> common;
    std::set_intersection(a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(common));
    if (!common.empty())
        return;
    for (std::int32_t p : rG.mV2P[static_cast<std::size_t>(vert)]) {
        if (rG.Label(p) != min_ref || rG.PolyHasVert(p, rT[0]) || rG.PolyHasVert(p, rT[1]))
            continue;
        int orient = -1;
        for (std::int32_t e : rG.mV2E[static_cast<std::size_t>(vert)]) {
            const int o = rG.Orientation(vert, e);
            if (rG.PolyHasEdge(p, e) && o != rDir[static_cast<std::size_t>(rT[0])] &&
                o != rDir[static_cast<std::size_t>(rT[1])]) {
                orient = o;
                break;
            }
        }
        if (orient < 0)
            hds_fail("the convex cell at " + hds_str(rG.Vert(vert)) + " has no edge to cut");
        auto it = rMap.find(p);
        if (it == rMap.end()) {
            HdsInfo info;
            info.mType = HexTransition::Convex1;
            info.mKind = HdsKind::ConvS;
            info.mCuts[orient] = true;
            info.mScale = rG.PolyScale(p);
            info.mTVerts.push_back(rG.Vert(vert));
            info.mOrientations.push_back(orient);
            rMap[p] = info;
        } else if (it->second.mType == HexTransition::Convex1) {
            it->second.mType = HexTransition::Convex2;
            it->second.mCuts[orient] = true;
            it->second.mTVerts.push_back(rG.Vert(vert));
        } else if (it->second.mType == HexTransition::Convex2) {
            it->second.mType = HexTransition::Convex3;
            it->second.mCuts[orient] = true;
            it->second.mTVerts.push_back(rG.Vert(vert));
        }
    }
}

void hds_mark_flat(const HdsGrid& rG, std::int32_t T, const std::vector<int>& rDir,
                   HdsSchemes& rMap) {
    const int min_ref = rG.MinRef(T);
    for (std::int32_t p : rG.mV2P[static_cast<std::size_t>(T)]) {
        if (rG.Label(p) != min_ref)
            continue;
        HdsInfo info;
        info.mScale = rG.PolyScale(p);
        info.mOrientations.push_back(rDir[static_cast<std::size_t>(T)]);
        info.mTVerts.push_back(rG.Vert(T));
        info.mKind = HdsKind::FlatS;
        if (!rMap.count(p)) {
            info.mType = HexTransition::Flat;
            rMap[p] = info;
        } else if (hds_is(rMap, p,
                          {HexTransition::Flat, HexTransition::Convex1, HexTransition::Convex2,
                           HexTransition::Convex3})) {
            info.mType = HexTransition::FlatConvex;
            rMap[p] = info;
            for (std::int32_t q : rG.mP2P[static_cast<std::size_t>(p)]) {
                const auto it = rMap.find(q);
                if (it != rMap.end() && it->second.mType == HexTransition::Flat &&
                    (it->second.mKind == HdsKind::ConcS || it->second.mKind == HdsKind::CornS))
                    it->second.mType = HexTransition::FlatConvex;
            }
        }
    }
}

void hds_cut_flats(const HdsGrid& rG, HdsSchemes& rMap) {
    std::deque<std::int32_t> queue;
    for (const auto& kv : rMap)
        queue.push_back(kv.first);
    while (!queue.empty()) {
        const std::int32_t id = queue.front();
        queue.pop_front();
        const HdsInfo p = rMap.at(id);
        if (p.mType == HexTransition::FlatConvex)
            for (std::int32_t q : rG.mP2P[static_cast<std::size_t>(id)]) {
                const auto it = rMap.find(q);
                if (it != rMap.end() && it->second.mType == HexTransition::Flat &&
                    (it->second.mKind == HdsKind::ConcS || it->second.mKind == HdsKind::CornS) &&
                    p.mTVerts[0] != it->second.mTVerts[0]) {
                    it->second.mType = HexTransition::FlatConvex;
                    queue.push_back(q);
                }
            }
        if (p.mType == HexTransition::Flat && p.mKind == HdsKind::CornS) {
            int sides = 0;
            for (std::int32_t q : rG.mP2P[static_cast<std::size_t>(id)])
                if (hds_is(rMap, q, {HexTransition::VertSide}))
                    ++sides;
            if (sides > 2) {
                rMap.at(id).mType = HexTransition::FlatConvex;
                queue.push_back(id);
            }
        }
    }
}

HdsSchemes hds_mark(const HdsGrid& rG, const std::vector<std::uint8_t>& rIsT,
                    const std::vector<int>& rDir) {
    std::set<std::int32_t> flats;
    std::set<std::vector<std::int32_t>> concaves, corners, convexes;
    std::set<std::int32_t> not_flats;
    for (std::size_t p = 0; p < rG.mPolyVerts.size(); ++p) {
        std::vector<std::int32_t> ids;
        for (std::int32_t v : rG.mPolyVerts[p])
            if (rIsT[static_cast<std::size_t>(v)])
                ids.push_back(v);
        std::sort(ids.begin(), ids.end());
        bool flat = ids.size() == 1;
        if (ids.size() == 2)
            flat = flat || rG.EdgeId(ids[0], ids[1]) != -1;
        if (flat) {
            for (std::int32_t v : ids)
                if (!not_flats.count(v))
                    flats.insert(v);
        } else if (ids.size() == 2) {
            for (std::int32_t f : rG.mPolyFaces[p]) {
                const auto& fv = rG.mFaceVerts[static_cast<std::size_t>(f)];
                if (std::find(fv.begin(), fv.end(), ids[0]) == fv.end() ||
                    std::find(fv.begin(), fv.end(), ids[1]) == fv.end())
                    continue;
                if (rG.MaxRef(ids[0]) == rG.MaxRef(ids[1])) {
                    concaves.insert(ids);
                    for (std::int32_t v : ids) {
                        not_flats.insert(v);
                        flats.erase(v);
                    }
                    break;
                }
            }
        } else if (ids.size() == 3) {
            corners.insert(ids);
            for (std::int32_t v : ids) {
                not_flats.insert(v);
                flats.erase(v);
            }
        } else if (ids.size() > 3) {
            hds_fail("a cell at " + hds_str(rG.Vert(ids[0])) + " holds " +
                     std::to_string(ids.size()) + " transition vertices");
        }
    }
    for (std::size_t v = 0; v < rG.mVerts.size(); ++v) {
        if (rIsT[v])
            continue;
        std::vector<std::int32_t> adj;
        for (std::int32_t w : rG.mV2V[v])
            if (rIsT[static_cast<std::size_t>(w)])
                adj.push_back(w);
        if (adj.size() == 3) {
            int top = 0;
            for (std::int32_t w : adj)
                top = std::max(rG.MaxRef(w), top);
            for (std::size_t i = 0; i < 3; ++i)
                if (rG.MaxRef(adj[i]) != top) {
                    adj.erase(adj.begin() + static_cast<std::ptrdiff_t>(i));
                    break;
                }
        }
        if (adj.size() == 2 &&
            rDir[static_cast<std::size_t>(adj[0])] != rDir[static_cast<std::size_t>(adj[1])]) {
            std::sort(adj.begin(), adj.end());
            if (!concaves.count(adj)) {
                adj.push_back(static_cast<std::int32_t>(v));
                convexes.insert(adj);
            }
        }
    }
    HdsSchemes map;
    for (const auto& c : corners)
        hds_mark_concave_vert(rG, c, rDir, map);
    for (const auto& c : concaves)
        hds_mark_concave_edge(rG, c, rDir, map);
    for (const auto& c : convexes)
        hds_mark_convex(rG, c, rDir, map);
    for (std::int32_t f : flats)
        hds_mark_flat(rG, f, rDir, map);
    hds_cut_flats(rG, map);
    return map;
}

// ---- orientation (CinoLib's hex_transition_orient) ------------------------

using HdsMat = std::array<std::array<std::int64_t, 3>, 3>;

// A rotation by `Quarters` quarter turns about an axis, as CinoLib's
// `mat3d::ROT_3D` builds it (right-handed).
HdsMat hds_rotation(char Axis, int Quarters) {
    const std::int64_t u = Axis == 'x', v = Axis == 'y', w = Axis == 'z';
    const int q = ((Quarters % 4) + 4) % 4;
    const std::int64_t c = q == 0 ? 1 : (q == 2 ? -1 : 0);
    const std::int64_t s = q == 1 ? 1 : (q == 3 ? -1 : 0);
    HdsMat m{};
    m[0][0] = c + u * u * (1 - c);
    m[1][0] = w * s + v * u * (1 - c);
    m[2][0] = -v * s + w * u * (1 - c);
    m[0][1] = -w * s + u * v * (1 - c);
    m[1][1] = c + v * v * (1 - c);
    m[2][1] = u * s + w * v * (1 - c);
    m[0][2] = v * s + u * w * (1 - c);
    m[1][2] = -u * s + v * w * (1 - c);
    m[2][2] = c + w * w * (1 - c);
    return m;
}

void hds_apply(std::vector<HdsPoint>& rV, const HdsMat& rM) {
    for (HdsPoint& p : rV) {
        const HdsPoint q = p;
        for (std::size_t i = 0; i < 3; ++i)
            p[i] = rM[i][0] * q[0] + rM[i][1] * q[1] + rM[i][2] * q[2];
    }
}

void hds_rotate(std::vector<HdsPoint>& rV, char Axis, int Quarters) {
    hds_apply(rV, hds_rotation(Axis, Quarters));
}

// CinoLib's `reflect`: negate every axis the string does not name.
void hds_reflect(std::vector<HdsPoint>& rV, const std::string& rAxes) {
    for (HdsPoint& p : rV)
        for (std::size_t k = 0; k < 3; ++k)
            if (rAxes.find(static_cast<char>('x' + k)) == std::string::npos)
                p[k] = -p[k];
}

struct HdsPlaced {
    HexTransition mType;
    std::vector<HdsPoint> mVerts;  // centred, units of the cell's edge / 20
    HdsPoint mCentroid;
    std::int64_t mSpan;  // the cell's edge / 20, lattice units

    HdsPoint At(std::size_t I) const {
        return {mCentroid[0] + mVerts[I][0] * mSpan, mCentroid[1] + mVerts[I][1] * mSpan,
                mCentroid[2] + mVerts[I][2] * mSpan};
    }
};

HdsPlaced hds_load(HexTransition Type, const HdsPoint& rCentroid, std::int64_t Scale) {
    const HexDualScheme& s = hex_dual_scheme(Type);
    HdsPlaced out{Type, std::vector<HdsPoint>(s.mNumVerts), rCentroid, Scale / kHexDualUnit};
    for (std::size_t v = 0; v < s.mNumVerts; ++v)
        for (std::size_t k = 0; k < 3; ++k)
            out.mVerts[v][k] = s.mVerts[v * 3 + k] - kHexDualUnit / 2;
    return out;
}

void hds_try_reflections(HdsPlaced& rS, const HdsInfo& rInfo) {
    for (const char* axes : {"x", "y", "z", "xy", "xz", "yz"}) {
        HdsPlaced tmp = rS;
        hds_reflect(tmp.mVerts, axes);
        if (tmp.At(0) == rInfo.mTVerts[0]) {
            rS = tmp;
            return;
        }
    }
    hds_fail("no reflection puts a scheme on the transition vertex " + hds_str(rInfo.mTVerts[0]));
}

void hds_try_rotations(HdsPlaced& rS, const HdsInfo& rInfo) {
    HdsPlaced tmp = rS;
    for (int x = 0; x < 4; ++x) {
        for (int y = 0; y < 4; ++y) {
            for (int z = 0; z < 4; ++z) {
                hds_rotate(tmp.mVerts, 'z', 1);
                if ((rInfo.mTVerts[0] == tmp.At(0) && rInfo.mTVerts[1] == tmp.At(4)) ||
                    (rInfo.mTVerts[0] == tmp.At(4) && rInfo.mTVerts[1] == tmp.At(0))) {
                    rS = tmp;
                    return;
                }
            }
            hds_rotate(tmp.mVerts, 'y', 1);
        }
        hds_rotate(tmp.mVerts, 'x', 1);
    }
    hds_fail("no rotation puts a two-turn convex scheme on " + hds_str(rInfo.mTVerts[0]) + " and " +
             hds_str(rInfo.mTVerts[1]));
}

void hds_turn_to(HdsPlaced& rS, int Dir) {
    switch (Dir) {
        case kHdsPlusY:
            break;
        case kHdsPlusX:
            hds_rotate(rS.mVerts, 'z', -1);
            break;
        case kHdsPlusZ:
            hds_rotate(rS.mVerts, 'x', 1);
            break;
        case kHdsMinusX:
            hds_rotate(rS.mVerts, 'z', 1);
            break;
        case kHdsMinusY:
            hds_rotate(rS.mVerts, 'z', 2);
            break;
        case kHdsMinusZ:
            hds_rotate(rS.mVerts, 'x', -1);
            break;
        default:
            hds_fail("a scheme with no direction");
    }
}

HdsPlaced hds_orient_flat(HdsInfo Info, const HdsPoint& rCentroid) {
    HdsPlaced s = hds_load(Info.mType, rCentroid, Info.mScale);
    std::vector<int>& o = Info.mOrientations;
    auto& v = s.mVerts;
    if (Info.mKind == HdsKind::FlatS) {
        hds_turn_to(s, o[0]);
    } else if (Info.mKind == HdsKind::ConcS) {
        switch (o[0]) {
            case kHdsPlusX:
                hds_rotate(v, 'z', -1);
                if (o[1] == kHdsPlusY)
                    break;
                if (o[1] == kHdsPlusZ)
                    hds_rotate(v, 'x', -1);
                if (o[1] == kHdsMinusY)
                    hds_rotate(v, 'x', 2);
                if (o[1] == kHdsMinusZ)
                    hds_rotate(v, 'x', 1);
                break;
            case kHdsPlusY:
                if (o[1] == kHdsPlusX)
                    break;
                if (o[1] == kHdsPlusZ)
                    hds_rotate(v, 'y', 1);
                if (o[1] == kHdsMinusX)
                    hds_rotate(v, 'y', 2);
                if (o[1] == kHdsMinusZ)
                    hds_rotate(v, 'y', -1);
                break;
            case kHdsPlusZ:
                hds_rotate(v, 'x', 1);
                if (o[1] == kHdsPlusX)
                    break;
                if (o[1] == kHdsPlusY)
                    hds_rotate(v, 'z', -1);
                if (o[1] == kHdsMinusX)
                    hds_rotate(v, 'z', 2);
                if (o[1] == kHdsMinusY)
                    hds_rotate(v, 'z', 1);
                break;
            case kHdsMinusX:
                hds_rotate(v, 'z', 1);
                if (o[1] == kHdsPlusY)
                    break;
                if (o[1] == kHdsPlusZ)
                    hds_rotate(v, 'x', 1);
                if (o[1] == kHdsMinusY)
                    hds_rotate(v, 'x', 2);
                if (o[1] == kHdsMinusZ)
                    hds_rotate(v, 'x', -1);
                break;
            case kHdsMinusY:
                hds_rotate(v, 'z', 2);
                if (o[1] == kHdsMinusX)
                    break;
                if (o[1] == kHdsMinusZ)
                    hds_rotate(v, 'y', 1);
                if (o[1] == kHdsPlusX)
                    hds_rotate(v, 'y', 2);
                if (o[1] == kHdsPlusZ)
                    hds_rotate(v, 'y', -1);
                break;
            case kHdsMinusZ:
                hds_rotate(v, 'x', -1);
                if (o[1] == kHdsMinusX)
                    break;
                if (o[1] == kHdsMinusY)
                    hds_rotate(v, 'z', -1);
                if (o[1] == kHdsPlusX)
                    hds_rotate(v, 'z', 2);
                if (o[1] == kHdsPlusY)
                    hds_rotate(v, 'z', 1);
                break;
            default:
                hds_fail("a flat scheme with no direction");
        }
    } else {
        std::sort(o.begin() + 1, o.end());
        auto is = [&](int A, int B) { return o[1] == A && o[2] == B; };
        switch (o[0]) {
            case kHdsPlusX:
                hds_rotate(v, 'z', -1);
                hds_reflect(v, "xz");
                if (is(kHdsPlusY, kHdsMinusZ))
                    hds_reflect(v, "xy");
                if (is(kHdsPlusZ, kHdsMinusY))
                    hds_reflect(v, "xz");
                if (is(kHdsMinusY, kHdsMinusZ)) {
                    hds_reflect(v, "xy");
                    hds_reflect(v, "xz");
                }
                break;
            case kHdsPlusY:
                hds_rotate(v, 'y', 1);
                hds_reflect(v, "xy");
                if (is(kHdsPlusX, kHdsMinusZ))
                    hds_reflect(v, "xy");
                if (is(kHdsMinusX, kHdsMinusZ))
                    hds_reflect(v, "y");
                if (is(kHdsPlusZ, kHdsMinusX))
                    hds_reflect(v, "yz");
                break;
            case kHdsPlusZ:
                hds_rotate(v, 'x', 1);
                hds_reflect(v, "xz");
                if (is(kHdsPlusY, kHdsMinusX))
                    hds_reflect(v, "yz");
                if (is(kHdsPlusX, kHdsMinusY))
                    hds_reflect(v, "xz");
                if (is(kHdsMinusX, kHdsMinusY)) {
                    hds_reflect(v, "yz");
                    hds_reflect(v, "xz");
                }
                break;
            case kHdsMinusX:
                hds_rotate(v, 'z', 1);
                if (is(kHdsPlusY, kHdsMinusZ))
                    hds_reflect(v, "xy");
                if (is(kHdsPlusZ, kHdsMinusY))
                    hds_reflect(v, "xz");
                if (is(kHdsMinusY, kHdsMinusZ)) {
                    hds_reflect(v, "yz");
                    hds_reflect(v, "xz");
                }
                break;
            case kHdsMinusY:
                hds_rotate(v, 'y', 1);
                hds_reflect(v, "xy");
                hds_reflect(v, "xz");
                if (is(kHdsPlusX, kHdsMinusZ))
                    hds_reflect(v, "xy");
                if (is(kHdsMinusX, kHdsMinusZ))
                    hds_reflect(v, "y");
                if (is(kHdsPlusZ, kHdsMinusX))
                    hds_reflect(v, "yz");
                break;
            case kHdsMinusZ:
                hds_rotate(v, 'x', -1);
                if (is(kHdsPlusY, kHdsMinusX))
                    hds_reflect(v, "yz");
                if (is(kHdsPlusX, kHdsMinusY))
                    hds_reflect(v, "xz");
                if (is(kHdsMinusX, kHdsMinusY)) {
                    hds_reflect(v, "yz");
                    hds_reflect(v, "xz");
                }
                break;
            default:
                hds_fail("a corner flat scheme with no direction");
        }
    }
    if (s.At(0) != Info.mTVerts[0])
        hds_try_reflections(s, Info);
    return s;
}

HdsPlaced hds_orient_convex(const HdsInfo& rInfo, const HdsPoint& rCentroid) {
    HdsPlaced s = hds_load(rInfo.mType, rCentroid, rInfo.mScale);
    auto& v = s.mVerts;
    if (rInfo.mType == HexTransition::Convex1) {
        hds_turn_to(s, rInfo.mOrientations[0]);
        if (s.At(0) != rInfo.mTVerts[0])
            hds_try_reflections(s, rInfo);
    } else if (rInfo.mType == HexTransition::Convex2) {
        if (!((rInfo.mTVerts[0] == s.At(0) && rInfo.mTVerts[1] == s.At(4)) ||
              (rInfo.mTVerts[0] == s.At(4) && rInfo.mTVerts[1] == s.At(0))))
            hds_try_rotations(s, rInfo);
    } else {
        auto cuts = [&](int A, int B, int C) {
            return rInfo.Cut(A) && rInfo.Cut(B) && rInfo.Cut(C);
        };
        if (cuts(kHdsPlusX, kHdsPlusY, kHdsPlusZ)) {
            hds_rotate(v, 'y', 2);
        } else if (cuts(kHdsPlusX, kHdsPlusZ, kHdsMinusY)) {
            hds_rotate(v, 'x', -1);
            hds_rotate(v, 'y', 2);
        } else if (cuts(kHdsPlusZ, kHdsMinusX, kHdsMinusY)) {
            hds_rotate(v, 'y', 1);
            hds_rotate(v, 'z', 1);
        } else if (cuts(kHdsPlusX, kHdsMinusY, kHdsMinusZ)) {
            hds_rotate(v, 'x', 2);
            hds_rotate(v, 'y', 2);
        } else if (cuts(kHdsPlusX, kHdsPlusY, kHdsMinusZ)) {
            hds_rotate(v, 'y', -1);
        } else if (cuts(kHdsPlusY, kHdsPlusZ, kHdsMinusX)) {
            hds_rotate(v, 'x', 1);
        } else if (cuts(kHdsMinusX, kHdsMinusY, kHdsMinusZ)) {
            hds_rotate(v, 'z', 1);
        } else if (!cuts(kHdsPlusY, kHdsMinusX, kHdsMinusZ)) {
            hds_fail("a three-turn convex scheme at " + hds_str(rInfo.mTVerts[0]) +
                     " whose cuts are not three axes of one corner");
        }
    }
    return s;
}

HdsPlaced hds_orient_concave_edge(const HdsInfo& rInfo, const HdsPoint& rCentroid) {
    HdsPlaced s = hds_load(rInfo.mType, rCentroid, rInfo.mScale);
    auto& v = s.mVerts;
    hds_reflect(v, "xz");
    switch (rInfo.mOrientations[0]) {
        case kHdsPlusY:
            hds_reflect(v, "xz");
            break;
        case kHdsPlusX:
            hds_rotate(v, 'z', -1);
            hds_reflect(v, "yz");
            break;
        case kHdsPlusZ:
            hds_rotate(v, 'x', 1);
            hds_reflect(v, "xy");
            break;
        case kHdsMinusX:
            hds_rotate(v, 'z', 1);
            hds_reflect(v, "yz");
            break;
        case kHdsMinusY:
            hds_rotate(v, 'z', 2);
            hds_reflect(v, "xz");
            break;
        case kHdsMinusZ:
            hds_rotate(v, 'x', -1);
            hds_reflect(v, "xy");
            break;
        default:
            hds_fail("a concave edge scheme with no direction");
    }
    if (s.At(0) != rInfo.mTVerts[0])
        hds_try_reflections(s, rInfo);
    return s;
}

HdsPlaced hds_orient_concave_vert(HdsInfo Info, const HdsPoint& rCentroid) {
    HdsPlaced s = hds_load(Info.mType, rCentroid, Info.mScale);
    auto& v = s.mVerts;
    std::vector<int>& o = Info.mOrientations;
    std::sort(o.begin(), o.end());
    auto is = [&](int A, int B, int C) { return o[0] == A && o[1] == B && o[2] == C; };
    if (is(kHdsPlusX, kHdsPlusY, kHdsMinusZ)) {
        hds_reflect(v, "xy");
    } else if (is(kHdsPlusY, kHdsMinusX, kHdsMinusZ)) {
        hds_reflect(v, "y");
    } else if (is(kHdsPlusY, kHdsPlusZ, kHdsMinusX)) {
        hds_reflect(v, "yz");
    } else if (is(kHdsPlusX, kHdsPlusZ, kHdsMinusY)) {
        hds_reflect(v, "xz");
    } else if (is(kHdsPlusX, kHdsMinusY, kHdsMinusZ)) {
        hds_reflect(v, "xz");
        hds_reflect(v, "xy");
    } else if (is(kHdsMinusX, kHdsMinusY, kHdsMinusZ)) {
        hds_reflect(v, "xz");
        hds_reflect(v, "y");
    } else if (is(kHdsPlusZ, kHdsMinusX, kHdsMinusY)) {
        hds_reflect(v, "xz");
        hds_reflect(v, "yz");
    }
    return s;
}

HdsPlaced hds_orient_concave_vert_side(HdsInfo Info, const HdsPoint& rCentroid) {
    HdsPlaced s = hds_load(Info.mType, rCentroid, Info.mScale);
    auto& v = s.mVerts;
    const std::size_t tv = 2;
    std::vector<int>& o = Info.mOrientations;
    if (o.size() != 3 || Info.mTVerts.empty())
        hds_fail("a cell beside a concave corner has " + std::to_string(Info.mTVerts.size()) +
                 " transition vertices and " + std::to_string(o.size()) + " directions");
    std::sort(o.begin(), o.begin() + 2);
    auto is = [&](int A, int B) { return o[0] == A && o[1] == B; };
    auto unless_at = [&](const char* rAxes) {
        if (s.At(tv) != Info.mTVerts[0])
            hds_reflect(v, rAxes);
    };
    if (is(kHdsPlusX, kHdsPlusY)) {
        hds_rotate(v, 'y', 1);
        hds_rotate(v, 'x', 1);
        unless_at("xy");
    } else if (is(kHdsPlusX, kHdsPlusZ)) {
        if (o[2] == kHdsMinusY)
            hds_reflect(v, "xz");
    } else if (is(kHdsPlusX, kHdsMinusY)) {
        hds_rotate(v, 'y', 1);
        hds_rotate(v, 'x', 1);
        hds_reflect(v, "xz");
        unless_at("xy");
    } else if (is(kHdsPlusX, kHdsMinusZ)) {
        hds_reflect(v, "xy");
        if (o[2] == kHdsMinusY)
            hds_reflect(v, "xz");
    } else if (is(kHdsPlusY, kHdsPlusZ)) {
        hds_rotate(v, 'y', -1);
        hds_rotate(v, 'z', -1);
        unless_at("yz");
    } else if (is(kHdsPlusY, kHdsMinusX)) {
        hds_rotate(v, 'y', 1);
        hds_rotate(v, 'x', 1);
        hds_reflect(v, "yz");
        unless_at("xy");
    } else if (is(kHdsPlusY, kHdsMinusZ)) {
        hds_rotate(v, 'y', -1);
        hds_rotate(v, 'z', -1);
        hds_reflect(v, "xy");
        unless_at("yz");
    } else if (is(kHdsPlusZ, kHdsMinusX)) {
        hds_reflect(v, "yz");
        if (o[2] == kHdsMinusY)
            hds_reflect(v, "xz");
    } else if (is(kHdsPlusZ, kHdsMinusY)) {
        hds_rotate(v, 'y', -1);
        hds_rotate(v, 'z', -1);
        hds_reflect(v, "xz");
        unless_at("yz");
    } else if (is(kHdsMinusX, kHdsMinusY)) {
        hds_rotate(v, 'y', 1);
        hds_rotate(v, 'x', 1);
        hds_reflect(v, "xz");
        hds_reflect(v, "yz");
        unless_at("xy");
    } else if (is(kHdsMinusX, kHdsMinusZ)) {
        hds_reflect(v, "y");
        if (o[2] == kHdsMinusY)
            hds_reflect(v, "xz");
    } else if (is(kHdsMinusY, kHdsMinusZ)) {
        hds_rotate(v, 'y', -1);
        hds_rotate(v, 'z', -1);
        hds_reflect(v, "xz");
        hds_reflect(v, "xy");
        unless_at("yz");
    }
    return s;
}

HdsPlaced hds_orient(const HdsInfo& rInfo, const HdsPoint& rCentroid) {
    if (rInfo.mScale % kHexDualUnit != 0 || rInfo.mTVerts.empty())
        hds_fail("a scheme with no transition vertex or a cell off the lattice");
    switch (rInfo.mType) {
        case HexTransition::Flat:
        case HexTransition::FlatConvex:
            return hds_orient_flat(rInfo, rCentroid);
        case HexTransition::Convex1:
        case HexTransition::Convex2:
        case HexTransition::Convex3:
            return hds_orient_convex(rInfo, rCentroid);
        case HexTransition::Edge:
            return hds_orient_concave_edge(rInfo, rCentroid);
        case HexTransition::VertCenter:
            return hds_orient_concave_vert(rInfo, rCentroid);
        case HexTransition::VertSide:
            return hds_orient_concave_vert_side(rInfo, rCentroid);
    }
    hds_fail("an unknown scheme");
}

// ---- the conforming primal mesh -------------------------------------------

struct HdsPrimal {
    std::vector<HdsPoint> mVerts;
    std::unordered_map<HdsPoint, std::int32_t, HdsPointHash> mVertId;
    std::vector<std::vector<std::int32_t>> mFaces;
    std::map<std::vector<std::int32_t>, std::int32_t> mFaceId;
    std::vector<std::vector<std::int32_t>> mPolys;  // face ids
    std::vector<std::int32_t> mLabel;
    std::map<std::vector<std::int32_t>, std::int32_t> mPolyId;

    std::int32_t Vertex(const HdsPoint& rP) {
        const auto it = mVertId.find(rP);
        if (it != mVertId.end())
            return it->second;
        const auto id = static_cast<std::int32_t>(mVerts.size());
        mVerts.push_back(rP);
        mVertId.emplace(rP, id);
        return id;
    }
    std::int32_t Face(const std::vector<std::int32_t>& rVerts) {
        std::vector<std::int32_t> key = rVerts;
        std::sort(key.begin(), key.end());
        const auto it = mFaceId.find(key);
        if (it != mFaceId.end())
            return it->second;
        const auto id = static_cast<std::int32_t>(mFaces.size());
        mFaces.push_back(rVerts);
        mFaceId.emplace(std::move(key), id);
        return id;
    }
    void Poly(const std::vector<std::int32_t>& rFaces, std::int32_t Label) {
        std::vector<std::int32_t> key = rFaces;
        std::sort(key.begin(), key.end());
        if (mPolyId.count(key))
            return;
        mPolyId.emplace(std::move(key), static_cast<std::int32_t>(mPolys.size()));
        mPolys.push_back(rFaces);
        mLabel.push_back(Label);
    }
};

HdsPrimal hds_primal(const HdsGrid& rG, const HdsSchemes& rMap,
                     std::array<std::int64_t, kNumHexTransitions>& rCount) {
    HdsPrimal m;
    for (const HdsPoint& p : rG.mVerts)
        m.Vertex(p);
    for (std::size_t p = 0; p < rG.mPolyVerts.size(); ++p) {
        const auto pid = static_cast<std::int32_t>(p);
        const auto it = rMap.find(pid);
        if (it == rMap.end()) {
            std::vector<std::int32_t> faces;
            for (std::int32_t f : rG.mPolyFaces[p]) {
                const auto& fv = rG.mFaceVerts[static_cast<std::size_t>(f)];
                faces.push_back(m.Face({fv[0], fv[1], fv[2], fv[3]}));
            }
            m.Poly(faces, rG.Label(pid));
            continue;
        }
        const auto& pv = rG.mPolyVerts[p];
        const HdsPoint& lo = rG.Vert(pv[0]);
        const HdsPoint& hi = rG.Vert(pv[7]);
        const HdsPoint centroid{(lo[0] + hi[0]) / 2, (lo[1] + hi[1]) / 2, (lo[2] + hi[2]) / 2};
        const HdsPlaced s = hds_orient(it->second, centroid);
        ++rCount[static_cast<std::size_t>(it->second.mType)];
        const HexDualScheme& t = hex_dual_scheme(s.mType);
        std::vector<std::int32_t> vid(t.mNumVerts);
        for (std::size_t v = 0; v < t.mNumVerts; ++v)
            vid[v] = m.Vertex(s.At(v));
        std::vector<std::int32_t> fid(t.mNumFaces);
        for (std::size_t f = 0, at = 0; f < t.mNumFaces; ++f) {
            std::vector<std::int32_t> fv;
            for (; at < t.mFaceEnds[f]; ++at)
                fv.push_back(vid[t.mFaces[at]]);
            fid[f] = m.Face(fv);
        }
        for (std::size_t q = 0, at = 0; q < t.mNumPolys; ++q) {
            std::vector<std::int32_t> faces;
            for (; at < t.mPolyEnds[q]; ++at)
                faces.push_back(fid[t.mPolys[at]]);
            m.Poly(faces, rG.Label(pid));
        }
    }
    return m;
}

// ---- the dual (CinoLib's dual_mesh, without clipped cells) ----------------

double hds_hex_volume(const std::vector<Vec3>& rP, const std::array<std::int64_t, 8>& rH) {
    static constexpr int kTets[6][4] = {{0, 1, 2, 6}, {0, 2, 3, 6}, {0, 3, 7, 6},
                                        {0, 7, 4, 6}, {0, 4, 5, 6}, {0, 5, 1, 6}};
    double vol = 0.0;
    for (const auto& t : kTets) {
        const Vec3& a = rP[static_cast<std::size_t>(rH[static_cast<std::size_t>(t[0])])];
        const Vec3 b =
            vec3_sub(rP[static_cast<std::size_t>(rH[static_cast<std::size_t>(t[1])])], a);
        const Vec3 c =
            vec3_sub(rP[static_cast<std::size_t>(rH[static_cast<std::size_t>(t[2])])], a);
        const Vec3 d =
            vec3_sub(rP[static_cast<std::size_t>(rH[static_cast<std::size_t>(t[3])])], a);
        vol += vec3_dot(b, vec3_cross(c, d)) / 6.0;
    }
    return vol;
}

HexDualMesh hds_dual(const HdsPrimal& rM, const LinearOctree& rTree) {
    HexDualMesh out;
    const std::size_t nv = rM.mVerts.size(), nf = rM.mFaces.size(), np = rM.mPolys.size();
    std::vector<std::vector<std::int32_t>> f2p(nf);
    for (std::size_t p = 0; p < np; ++p)
        for (std::int32_t f : rM.mPolys[p])
            f2p[static_cast<std::size_t>(f)].push_back(static_cast<std::int32_t>(p));
    std::array<std::int64_t, 3> extent{};
    for (std::size_t k = 0; k < 3; ++k)
        extent[k] = rTree.Extent(k) * kHexDualUnit;
    auto on_boundary = [&](const HdsPoint& rP) {
        for (std::size_t k = 0; k < 3; ++k)
            if (rP[k] == 0 || rP[k] == extent[k])
                return true;
        return false;
    };
    // Every face has two cells, except on the forest's boundary.
    for (std::size_t f = 0; f < nf; ++f) {
        const std::size_t c = f2p[f].size();
        if (c == 2 || c == 0)
            continue;
        // A face of one cell must lie on one of the forest's boundary planes.
        bool outer = false;
        for (std::size_t k = 0; k < 3 && c == 1 && !outer; ++k) {
            bool lo = true, hi = true;
            for (std::int32_t v : rM.mFaces[f]) {
                lo = lo && rM.mVerts[static_cast<std::size_t>(v)][k] == 0;
                hi = hi && rM.mVerts[static_cast<std::size_t>(v)][k] == extent[k];
            }
            outer = lo || hi;
        }
        if (!outer)
            hds_fail("the schemes do not conform: a face at " +
                     hds_str(rM.mVerts[static_cast<std::size_t>(rM.mFaces[f][0])]) + " has " +
                     std::to_string(c) + " cells");
    }
    // Each cell's dual point, the average of its vertices.
    std::vector<Vec3> centroid(np);
    std::vector<std::vector<std::int32_t>> v2p(nv);
    for (std::size_t p = 0; p < np; ++p) {
        std::vector<std::int32_t> verts;
        for (std::int32_t f : rM.mPolys[p])
            for (std::int32_t v : rM.mFaces[static_cast<std::size_t>(f)])
                verts.push_back(v);
        std::sort(verts.begin(), verts.end());
        verts.erase(std::unique(verts.begin(), verts.end()), verts.end());
        Vec3 c{0.0, 0.0, 0.0};
        for (std::int32_t v : verts) {
            const HdsPoint& q = rM.mVerts[static_cast<std::size_t>(v)];
            for (std::size_t k = 0; k < 3; ++k)
                c[k] += static_cast<double>(q[k]);
            v2p[static_cast<std::size_t>(v)].push_back(static_cast<std::int32_t>(p));
        }
        for (std::size_t k = 0; k < 3; ++k)
            c[k] /= static_cast<double>(verts.size());
        centroid[p] = c;
    }
    std::vector<std::vector<std::int32_t>> v2f(nv);
    for (std::size_t f = 0; f < nf; ++f)
        for (std::int32_t v : rM.mFaces[f])
            v2f[static_cast<std::size_t>(v)].push_back(static_cast<std::int32_t>(f));

    std::vector<std::int64_t> point_of(np, -1);
    for (std::size_t v = 0; v < nv; ++v) {
        const HdsPoint& at = rM.mVerts[v];
        if (on_boundary(at) || v2p[v].empty())
            continue;
        const std::vector<std::int32_t>& cells = v2p[v];
        if (cells.size() != 8)
            hds_fail("the interior vertex " + hds_str(at) + " has " + std::to_string(cells.size()) +
                     " cells around it, so its dual is not a hexahedron");
        // The cells' adjacency across the faces at the vertex: a cube's graph.
        std::array<std::array<int, 3>, 8> nb{};
        std::array<int, 8> deg{};
        auto local = [&](std::int32_t P) {
            return static_cast<int>(std::find(cells.begin(), cells.end(), P) - cells.begin());
        };
        bool ok = true;
        for (std::int32_t f : v2f[v]) {
            const auto& fp = f2p[static_cast<std::size_t>(f)];
            if (fp.size() != 2) {
                ok = false;
                break;
            }
            const int a = local(fp[0]), b = local(fp[1]);
            if (deg[static_cast<std::size_t>(a)] == 3 || deg[static_cast<std::size_t>(b)] == 3) {
                ok = false;
                break;
            }
            nb[static_cast<std::size_t>(a)]
              [static_cast<std::size_t>(deg[static_cast<std::size_t>(a)]++)] = b;
            nb[static_cast<std::size_t>(b)]
              [static_cast<std::size_t>(deg[static_cast<std::size_t>(b)]++)] = a;
        }
        std::array<int, 8> label;
        label.fill(-1);
        if (ok) {
            label[0] = 0;
            for (int k = 0; k < 3 && ok; ++k) {
                if (deg[0] != 3) {
                    ok = false;
                    break;
                }
                label[static_cast<std::size_t>(nb[0][static_cast<std::size_t>(k)])] = 1 << k;
            }
            // Two passes: the cells beside two labelled ones, then the last.
            for (int pass = 0; pass < 2 && ok; ++pass)
                for (std::size_t c = 0; c < 8; ++c) {
                    if (label[c] >= 0)
                        continue;
                    int bits = 0, count = 0;
                    for (int k = 0; k < deg[c]; ++k) {
                        const int l =
                            label[static_cast<std::size_t>(nb[c][static_cast<std::size_t>(k)])];
                        if (l > 0 && (pass == 0 ? (l & (l - 1)) == 0 : (l & (l - 1)) != 0)) {
                            bits |= l;
                            ++count;
                        }
                    }
                    if (count >= 2)
                        label[c] = bits;
                }
            std::array<int, 8> seen{};
            for (std::size_t c = 0; c < 8 && ok; ++c) {
                if (label[c] < 0 || label[c] > 7 || deg[c] != 3 ||
                    seen[static_cast<std::size_t>(label[c])]++) {
                    ok = false;
                    break;
                }
                for (int k = 0; k < 3; ++k) {
                    const int x =
                        label[c] ^
                        label[static_cast<std::size_t>(nb[c][static_cast<std::size_t>(k)])];
                    if (x != 1 && x != 2 && x != 4)
                        ok = false;
                }
            }
        }
        if (!ok)
            hds_fail("the cells around the interior vertex " + hds_str(at) +
                     " do not meet as a cube's eight corners, so its dual is not a hexahedron");
        std::array<std::int32_t, 8> by_label{};
        for (std::size_t c = 0; c < 8; ++c)
            by_label[static_cast<std::size_t>(label[c])] = cells[c];
        std::array<std::int64_t, 8> hex{};
        static constexpr int kVtk[8] = {0, 1, 3, 2, 4, 5, 7, 6};
        static constexpr int kMirror[8] = {0, 1, 5, 4, 2, 3, 7, 6};
        for (std::size_t k = 0; k < 8; ++k) {
            const std::int32_t p = by_label[static_cast<std::size_t>(kVtk[k])];
            std::int64_t& r = point_of[static_cast<std::size_t>(p)];
            if (r < 0) {
                r = static_cast<std::int64_t>(out.mPoints.size());
                out.mPoints.push_back(centroid[static_cast<std::size_t>(p)]);
            }
            hex[k] = r;
        }
        if (hds_hex_volume(out.mPoints, hex) < 0.0)
            for (std::size_t k = 0; k < 8; ++k)
                hex[k] = point_of[static_cast<std::size_t>(
                    by_label[static_cast<std::size_t>(kMirror[k])])];
        int level = 0;
        for (std::int32_t p : cells)
            level = std::max(level, rM.mLabel[static_cast<std::size_t>(p)]);
        out.mHexes.push_back(hex);
        out.mCentre.push_back(at);
        out.mLevel.push_back(level);
    }
    return out;
}

}  // namespace

HexDualMesh hex_dual_mesh(const LinearOctree& rTree) {
    const OctreeAdjacency full = OctreeAdjacency::Full;
    std::string why;
    if (!linear_octree_check(rTree, &full, &why))
        hds_fail("the tree is not balanced across faces, edges and corners: " + why);
    if (!linear_octree_is_paired(rTree, &why))
        hds_fail("the tree is not paired: " + why);
    const HdsGrid grid = hds_grid(rTree);
    std::vector<std::uint8_t> is_t;
    std::vector<int> dir;
    hds_transition_vertices(rTree, grid, is_t, dir);
    const HdsSchemes map = hds_mark(grid, is_t, dir);
    std::array<std::int64_t, kNumHexTransitions> count{};
    const HdsPrimal primal = hds_primal(grid, map, count);
    HexDualMesh out = hds_dual(primal, rTree);
    out.mNumSchemes = count;
    out.mNumPrimalCells = static_cast<std::int64_t>(primal.mPolys.size());
    for (std::uint8_t t : is_t)
        out.mNumTransitionVertices += t;
    return out;
}

}  // namespace detail
}  // namespace meshioplusplus

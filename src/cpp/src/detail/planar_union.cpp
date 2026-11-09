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
// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

// Project includes
#include "planar_union.hpp"
#include "bucket_build.hpp"
#include "cdt2d.hpp"
#include "predicates.hpp"
#include "meshioplusplus/detail/bucket_table.hpp"
#include "meshioplusplus/detail/spatial_hash.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

using PuVec2 = std::array<double, 2>;
using PuSeg = std::array<std::int64_t, 2>;
using PuTri = std::array<std::int64_t, 3>;

constexpr int kPuMaxRounds = 32;

// The merged point set: a point within the tolerance of an earlier one is that
// one (the lowest id when several are), found through a grid of the
// tolerance's size; with no tolerance only exact coincidence merges.
class PuPoints {
public:
    explicit PuPoints(double Tolerance) : mTol(Tolerance) {}

    std::int64_t Add(const PuVec2& rP) {
        if (!(mTol > 0.0)) {
            const auto it = mExact.find(rP);
            if (it != mExact.end())
                return it->second;
            const auto id = static_cast<std::int64_t>(mXy.size());
            mExact.emplace(rP, id);
            mXy.push_back(rP);
            return id;
        }
        const GridKey k{grid_quantize(rP[0], mTol), grid_quantize(rP[1], mTol), 0};
        std::int64_t best = -1;
        for (std::int64_t dx = -1; dx <= 1; ++dx)
            for (std::int64_t dy = -1; dy <= 1; ++dy) {
                const auto it = mGrid.find(GridKey{k.x + dx, k.y + dy, 0});
                if (it == mGrid.end())
                    continue;
                for (const std::int64_t id : it->second) {
                    const PuVec2& q = mXy[static_cast<std::size_t>(id)];
                    const double ex = q[0] - rP[0], ey = q[1] - rP[1];
                    if (ex * ex + ey * ey <= mTol * mTol && (best < 0 || id < best))
                        best = id;
                }
            }
        if (best >= 0)
            return best;
        const auto id = static_cast<std::int64_t>(mXy.size());
        mGrid[k].push_back(id);
        mXy.push_back(rP);
        return id;
    }

    const PuVec2& operator[](std::int64_t I) const { return mXy[static_cast<std::size_t>(I)]; }
    const std::vector<PuVec2>& All() const { return mXy; }

private:
    double mTol;
    std::vector<PuVec2> mXy;
    std::map<PuVec2, std::int64_t> mExact;
    std::unordered_map<GridKey, std::vector<std::int64_t>, GridKeyHash> mGrid;
};

std::uint64_t pu_key(std::int64_t A, std::int64_t B) {
    return (static_cast<std::uint64_t>(A) << 32) ^ static_cast<std::uint64_t>(B);
}

// A soup's non-degenerate triangles over merged ids, counter-clockwise, and
// its candidate boundary edges (see the file comment), appended to rSegs.
std::vector<PuTri> pu_soup(const PlanarSoup& rSoup, PuPoints& rPts, std::vector<PuSeg>& rSegs) {
    std::vector<std::int64_t> id(rSoup.mPoints.size());
    for (std::size_t i = 0; i < rSoup.mPoints.size(); ++i)
        id[i] = rPts.Add(rSoup.mPoints[i]);
    std::vector<PuTri> tris;
    struct PuEdgeUse {
        int mCount = 0;
        int mLeft = 0;   // triangles on the left of min -> max
        int mRight = 0;  // and on its right
    };
    std::map<PuSeg, PuEdgeUse> uses;
    for (const PuTri& r_t : rSoup.mTriangles) {
        PuTri t{};
        for (std::size_t k = 0; k < 3; ++k) {
            const std::int64_t v = r_t[k];
            if (v < 0 || v >= static_cast<std::int64_t>(id.size()))
                return {};
            t[k] = id[static_cast<std::size_t>(v)];
        }
        if (t[0] == t[1] || t[1] == t[2] || t[0] == t[2])
            continue;
        const double o = orient2d(rPts[t[0]].data(), rPts[t[1]].data(), rPts[t[2]].data());
        if (o == 0.0)
            continue;
        if (o < 0.0)
            std::swap(t[1], t[2]);
        tris.push_back(t);
        for (std::size_t k = 0; k < 3; ++k) {
            const std::int64_t x = t[k], y = t[(k + 1) % 3];
            PuEdgeUse& r_u = uses[PuSeg{std::min(x, y), std::max(x, y)}];
            ++r_u.mCount;
            // The triangle is on the left of x -> y.
            if (x < y)
                ++r_u.mLeft;
            else
                ++r_u.mRight;
        }
    }
    for (const auto& [seg, use] : uses)
        if (!(use.mCount == 2 && use.mLeft == 1 && use.mRight == 1))
            rSegs.push_back(seg);
    return tris;
}

// The parameter of point rQ along the segment rA -> rB, on its dominant axis.
double pu_param(const PuVec2& rA, const PuVec2& rB, const PuVec2& rQ) {
    const std::size_t ax = std::fabs(rB[0] - rA[0]) >= std::fabs(rB[1] - rA[1]) ? 0 : 1;
    return (rQ[ax] - rA[ax]) / (rB[ax] - rA[ax]);
}

bool pu_opposite(double A, double B) {
    return (A > 0.0 && B < 0.0) || (A < 0.0 && B > 0.0);
}

// Splits two meeting segments: at a touching endpoint, at a collinear
// overlap's endpoints, or at a proper crossing's point.
void pu_split_pair(const PuSeg& rP, const PuSeg& rQ, PuPoints& rPts,
                   std::vector<std::int64_t>& rSplitP, std::vector<std::int64_t>& rSplitQ) {
    const PuVec2 p0 = rPts[rP[0]], p1 = rPts[rP[1]], q0 = rPts[rQ[0]], q1 = rPts[rQ[1]];
    const double o1 = orient2d(p0.data(), p1.data(), q0.data());
    const double o2 = orient2d(p0.data(), p1.data(), q1.data());
    const double o3 = orient2d(q0.data(), q1.data(), p0.data());
    const double o4 = orient2d(q0.data(), q1.data(), p1.data());
    auto inside = [](const PuVec2& rA, const PuVec2& rB, const PuVec2& rX) {
        const double t = pu_param(rA, rB, rX);
        return t > 0.0 && t < 1.0;
    };
    if (pu_opposite(o1, o2) && pu_opposite(o3, o4)) {
        const double t = o3 / (o3 - o4);
        const std::int64_t x =
            rPts.Add(PuVec2{p0[0] + t * (p1[0] - p0[0]), p0[1] + t * (p1[1] - p0[1])});
        rSplitP.push_back(x);
        rSplitQ.push_back(x);
        return;
    }
    if (o1 == 0.0 && rQ[0] != rP[0] && rQ[0] != rP[1] && inside(p0, p1, q0))
        rSplitP.push_back(rQ[0]);
    if (o2 == 0.0 && rQ[1] != rP[0] && rQ[1] != rP[1] && inside(p0, p1, q1))
        rSplitP.push_back(rQ[1]);
    if (o3 == 0.0 && rP[0] != rQ[0] && rP[0] != rQ[1] && inside(q0, q1, p0))
        rSplitQ.push_back(rP[0]);
    if (o4 == 0.0 && rP[1] != rQ[0] && rP[1] != rQ[1] && inside(q0, q1, p1))
        rSplitQ.push_back(rP[1]);
}

// A uniform bucket grid over a soup's triangles, for the coverage test.
class PuCover {
public:
    PuCover(const std::vector<PuTri>& rTris, const PuPoints& rPts) : mTris(rTris), mPts(rPts) {
        if (mTris.empty())
            return;
        double extent = 0.0;
        std::vector<std::array<double, 4>> boxes(mTris.size());
        for (std::size_t i = 0; i < mTris.size(); ++i) {
            std::array<double, 4> b{rPts[mTris[i][0]][0], rPts[mTris[i][0]][1],
                                    rPts[mTris[i][0]][0], rPts[mTris[i][0]][1]};
            for (std::size_t k = 1; k < 3; ++k) {
                const PuVec2& r_p = rPts[mTris[i][k]];
                b[0] = std::min(b[0], r_p[0]);
                b[1] = std::min(b[1], r_p[1]);
                b[2] = std::max(b[2], r_p[0]);
                b[3] = std::max(b[3], r_p[1]);
            }
            extent += std::max(b[2] - b[0], b[3] - b[1]);
            boxes[i] = b;
        }
        double cell = extent / static_cast<double>(mTris.size());
        if (!(cell > 0.0) || !std::isfinite(cell))
            cell = 1.0;
        std::vector<GridKey> lo(mTris.size()), hi(mTris.size());
        for (std::size_t i = 0; i < mTris.size(); ++i) {
            lo[i] = GridKey{grid_quantize(boxes[i][0], cell), grid_quantize(boxes[i][1], cell), 0};
            hi[i] = GridKey{grid_quantize(boxes[i][2], cell), grid_quantize(boxes[i][3], cell), 0};
        }
        mTable = build_bucket_table(cell, lo, hi);
    }

    bool Covers(const PuVec2& rP) const {
        if (mTris.empty())
            return false;
        const double xyz[3] = {rP[0], rP[1], 0.0};
        BucketView ids(nullptr, 0);
        if (!mTable.Find(mTable.KeyOf(xyz), ids))
            return false;
        for (const std::int64_t i : ids) {
            const PuTri& r_t = mTris[static_cast<std::size_t>(i)];
            const double* a = mPts[r_t[0]].data();
            const double* b = mPts[r_t[1]].data();
            const double* c = mPts[r_t[2]].data();
            if (orient2d(a, b, rP.data()) >= 0.0 && orient2d(b, c, rP.data()) >= 0.0 &&
                orient2d(c, a, rP.data()) >= 0.0)
                return true;
        }
        return false;
    }

private:
    const std::vector<PuTri>& mTris;
    const PuPoints& mPts;
    BucketTable mTable;
};

PlanarRegion pu_fail(const std::string& rWhy) {
    PlanarRegion r;
    r.mOk = false;
    r.mWhy = rWhy;
    return r;
}

}  // namespace

double planar_loop_twice_area(const std::vector<std::array<double, 2>>& rPoints,
                              const std::vector<std::int64_t>& rLoop) {
    double a = 0.0;
    for (std::size_t i = 0; i < rLoop.size(); ++i) {
        const auto& r_p = rPoints[static_cast<std::size_t>(rLoop[i])];
        const auto& r_q = rPoints[static_cast<std::size_t>(rLoop[(i + 1) % rLoop.size()])];
        a += r_p[0] * r_q[1] - r_q[0] * r_p[1];
    }
    return a;
}

PlanarRegion planar_region(const PlanarSoup& rAdd, const PlanarSoup& rSubtract,
                           double MergeTolerance) {
    for (const PlanarSoup* p_s : {&rAdd, &rSubtract})
        for (const PuVec2& r_p : p_s->mPoints)
            if (!cdt2d_in_range(r_p[0]) || !cdt2d_in_range(r_p[1]))
                return pu_fail(
                    "a projected coordinate is outside the range the exact predicates "
                    "hold for (zero, or a magnitude in [2^-128, 2^128])");
    PuPoints pts(MergeTolerance);
    std::vector<PuSeg> segs;
    const std::vector<PuTri> add = pu_soup(rAdd, pts, segs);
    const std::vector<PuTri> sub = pu_soup(rSubtract, pts, segs);
    std::sort(segs.begin(), segs.end());
    segs.erase(std::unique(segs.begin(), segs.end()), segs.end());
    if (add.empty())
        return {};

    // Split the candidate edges where they meet, until none does.
    for (int round = 0;; ++round) {
        const auto crossings = cdt2d_find_crossings(pts.All(), segs);
        if (crossings.empty())
            break;
        if (round == kPuMaxRounds)
            return pu_fail("the projected edges still cross after " + std::to_string(kPuMaxRounds) +
                           " rounds of splitting");
        std::vector<std::vector<std::int64_t>> splits(segs.size());
        for (const auto& r_c : crossings)
            pu_split_pair(
                segs[static_cast<std::size_t>(r_c[0])], segs[static_cast<std::size_t>(r_c[1])], pts,
                splits[static_cast<std::size_t>(r_c[0])], splits[static_cast<std::size_t>(r_c[1])]);
        std::vector<PuSeg> next;
        next.reserve(segs.size() + 2 * crossings.size());
        for (std::size_t s = 0; s < segs.size(); ++s) {
            const PuSeg& r_s = segs[s];
            std::vector<std::pair<double, std::int64_t>> along{{0.0, r_s[0]}, {1.0, r_s[1]}};
            for (const std::int64_t x : splits[s])
                if (x != r_s[0] && x != r_s[1])
                    along.emplace_back(pu_param(pts[r_s[0]], pts[r_s[1]], pts[x]), x);
            std::sort(along.begin(), along.end());
            for (std::size_t k = 1; k < along.size(); ++k) {
                const std::int64_t a = along[k - 1].second, b = along[k].second;
                if (a != b)
                    next.push_back(PuSeg{std::min(a, b), std::max(a, b)});
            }
        }
        std::sort(next.begin(), next.end());
        next.erase(std::unique(next.begin(), next.end()), next.end());
        segs = std::move(next);
    }

    // The constrained triangulation of the split edges.
    std::vector<std::int64_t> used;
    for (const PuSeg& r_s : segs)
        used.insert(used.end(), r_s.begin(), r_s.end());
    std::sort(used.begin(), used.end());
    used.erase(std::unique(used.begin(), used.end()), used.end());
    if (used.size() < 3)
        return {};
    double lo_u = pts[used[0]][0], lo_v = pts[used[0]][1], hi_u = lo_u, hi_v = lo_v;
    std::vector<Cdt2d::Vec2> batch;
    batch.reserve(used.size());
    for (const std::int64_t i : used) {
        const PuVec2& r_p = pts[i];
        lo_u = std::min(lo_u, r_p[0]);
        lo_v = std::min(lo_v, r_p[1]);
        hi_u = std::max(hi_u, r_p[0]);
        hi_v = std::max(hi_v, r_p[1]);
        batch.push_back(r_p);
    }
    Cdt2d cdt(lo_u, lo_v, hi_u, hi_v);
    const std::vector<std::int64_t> vid = cdt.InsertBatch(batch);
    std::unordered_map<std::int64_t, std::int64_t> vertex_of;
    for (std::size_t i = 0; i < used.size(); ++i)
        vertex_of.emplace(used[i], vid[i]);
    for (const PuSeg& r_s : segs) {
        const Cdt2dStatus st = cdt.InsertSegment(vertex_of.at(r_s[0]), vertex_of.at(r_s[1]));
        if (st != Cdt2dStatus::Ok)
            return pu_fail(std::string("the planar triangulation refused an edge (") +
                           cdt2d_status_name(st) + ")");
    }

    // Components of triangles joined across unconstrained edges.
    const std::vector<Cdt2d::Tri> tris = cdt.AllTriangles();
    std::unordered_map<std::uint64_t, std::int64_t> tri_of;  // directed edge -> triangle
    tri_of.reserve(3 * tris.size());
    for (std::size_t t = 0; t < tris.size(); ++t)
        for (std::size_t k = 0; k < 3; ++k)
            tri_of.emplace(pu_key(tris[t][k], tris[t][(k + 1) % 3]), static_cast<std::int64_t>(t));
    std::unordered_set<std::uint64_t> constrained;
    for (const auto& r_e : cdt.ConstrainedEdges())
        constrained.insert(pu_key(r_e[0], r_e[1]));
    std::vector<std::int64_t> comp(tris.size(), -1);
    std::vector<char> comp_super;
    std::vector<std::int64_t> comp_best;
    std::vector<double> comp_area;
    for (std::size_t seed = 0; seed < tris.size(); ++seed) {
        if (comp[seed] >= 0)
            continue;
        const auto c = static_cast<std::int64_t>(comp_super.size());
        comp_super.push_back(0);
        comp_best.push_back(-1);
        comp_area.push_back(-1.0);
        std::vector<std::size_t> stack{seed};
        comp[seed] = c;
        while (!stack.empty()) {
            const std::size_t t = stack.back();
            stack.pop_back();
            const Cdt2d::Tri& r_t = tris[t];
            if (r_t[0] < Cdt2d::kFirstVertex || r_t[1] < Cdt2d::kFirstVertex ||
                r_t[2] < Cdt2d::kFirstVertex) {
                comp_super[static_cast<std::size_t>(c)] = 1;
            } else {
                const double a = orient2d(cdt.Point(r_t[0]).data(), cdt.Point(r_t[1]).data(),
                                          cdt.Point(r_t[2]).data());
                if (a > comp_area[static_cast<std::size_t>(c)]) {
                    comp_area[static_cast<std::size_t>(c)] = a;
                    comp_best[static_cast<std::size_t>(c)] = static_cast<std::int64_t>(t);
                }
            }
            for (std::size_t k = 0; k < 3; ++k) {
                const std::int64_t x = r_t[k], y = r_t[(k + 1) % 3];
                if (constrained.count(pu_key(std::min(x, y), std::max(x, y))))
                    continue;
                const auto it = tri_of.find(pu_key(y, x));
                if (it == tri_of.end() || comp[static_cast<std::size_t>(it->second)] >= 0)
                    continue;
                comp[static_cast<std::size_t>(it->second)] = c;
                stack.push_back(static_cast<std::size_t>(it->second));
            }
        }
    }
    const PuCover cover_add(add, pts), cover_sub(sub, pts);
    std::vector<char> comp_covered(comp_super.size(), 0);
    for (std::size_t c = 0; c < comp_super.size(); ++c) {
        if (comp_super[c] || comp_best[c] < 0)
            continue;
        const Cdt2d::Tri& r_t = tris[static_cast<std::size_t>(comp_best[c])];
        PuVec2 g{0.0, 0.0};
        for (const std::int64_t v : r_t) {
            g[0] += cdt.Point(v)[0] / 3.0;
            g[1] += cdt.Point(v)[1] / 3.0;
        }
        comp_covered[c] = cover_add.Covers(g) && !cover_sub.Covers(g) ? 1 : 0;
    }
    auto covered = [&](std::int64_t T) {
        return comp_covered[static_cast<std::size_t>(comp[static_cast<std::size_t>(T)])] != 0;
    };

    // The boundary half-edges, walked into loops with the region on the left.
    std::vector<std::pair<PuSeg, std::int64_t>> boundary;  // (x -> y), its triangle
    for (std::size_t t = 0; t < tris.size(); ++t) {
        if (!covered(static_cast<std::int64_t>(t)))
            continue;
        for (std::size_t k = 0; k < 3; ++k) {
            const std::int64_t x = tris[t][k], y = tris[t][(k + 1) % 3];
            const auto it = tri_of.find(pu_key(y, x));
            if (it == tri_of.end() || !covered(it->second))
                boundary.push_back({PuSeg{x, y}, static_cast<std::int64_t>(t)});
        }
    }
    std::sort(boundary.begin(), boundary.end());
    std::unordered_map<std::uint64_t, std::size_t> boundary_index;
    for (std::size_t i = 0; i < boundary.size(); ++i)
        boundary_index.emplace(pu_key(boundary[i].first[0], boundary[i].first[1]), i);
    auto next_of = [&](std::size_t I) -> std::int64_t {
        const std::int64_t y = boundary[I].first[1];
        std::int64_t cur = boundary[I].second;
        for (std::size_t guard = 0; guard <= tris.size(); ++guard) {
            const Cdt2d::Tri& r_t = tris[static_cast<std::size_t>(cur)];
            const std::size_t at = r_t[0] == y ? 0 : (r_t[1] == y ? 1 : 2);
            const std::int64_t w = r_t[(at + 1) % 3];
            const auto it = tri_of.find(pu_key(w, y));
            if (it == tri_of.end() || !covered(it->second)) {
                const auto b = boundary_index.find(pu_key(y, w));
                return b == boundary_index.end() ? -1 : static_cast<std::int64_t>(b->second);
            }
            cur = it->second;
        }
        return -1;
    };
    std::vector<char> walked(boundary.size(), 0);
    std::vector<std::vector<std::int64_t>> loops;
    for (std::size_t s = 0; s < boundary.size(); ++s) {
        if (walked[s])
            continue;
        std::vector<std::int64_t> loop;
        std::int64_t i = static_cast<std::int64_t>(s);
        while (i >= 0 && !walked[static_cast<std::size_t>(i)]) {
            walked[static_cast<std::size_t>(i)] = 1;
            loop.push_back(boundary[static_cast<std::size_t>(i)].first[0]);
            i = next_of(static_cast<std::size_t>(i));
        }
        if (i != static_cast<std::int64_t>(s))
            return pu_fail("a boundary loop of the covered region does not close");
        const auto first = std::min_element(loop.begin(), loop.end());
        std::rotate(loop.begin(), first, loop.end());
        loops.push_back(std::move(loop));
    }
    std::sort(loops.begin(), loops.end());

    // Compact the vertices used, in ascending triangulation id.
    PlanarRegion out;
    std::vector<std::int64_t> keep;
    for (std::size_t t = 0; t < tris.size(); ++t)
        if (covered(static_cast<std::int64_t>(t)))
            keep.insert(keep.end(), tris[t].begin(), tris[t].end());
    std::sort(keep.begin(), keep.end());
    keep.erase(std::unique(keep.begin(), keep.end()), keep.end());
    std::unordered_map<std::int64_t, std::int64_t> new_id;
    for (std::size_t i = 0; i < keep.size(); ++i) {
        new_id.emplace(keep[i], static_cast<std::int64_t>(i));
        out.mPoints.push_back(cdt.Point(keep[i]));
    }
    for (std::size_t t = 0; t < tris.size(); ++t)
        if (covered(static_cast<std::int64_t>(t)))
            out.mTriangles.push_back(
                PuTri{new_id.at(tris[t][0]), new_id.at(tris[t][1]), new_id.at(tris[t][2])});
    for (auto& r_l : loops) {
        for (std::int64_t& r_v : r_l)
            r_v = new_id.at(r_v);
        out.mLoops.push_back(std::move(r_l));
    }
    return out;
}

}  // namespace detail
}  // namespace meshioplusplus

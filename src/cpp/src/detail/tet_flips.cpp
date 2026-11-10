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

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

// Project includes
#include "tet_flips.hpp"

// Project includes (private, not installed)
#include "../operations/smooth_odt.hpp"

namespace meshioplusplus {
namespace detail {

// --- 3-2 sub-pass ------------------------------------------------------------
// Every interior edge shared by exactly three tets forming a ring is a
// candidate: replace the three tets with two capping the ring triangle, iff the
// two caps tile the same trigonal bipyramid (a signed-volume tiling test) and
// the worse cap beats the worst of the three old tets.
std::int64_t tet_flip_pass_32(std::vector<double>& rXyz, std::vector<TetFlipTet>& rTets,
                              double min_improve, double eps_len, double eps_vol6) {
    const std::size_t nt = rTets.size();
    std::vector<std::uint8_t> alive(nt, 1);

    // Every (edge key, tet) sorted (a total order: parallel_sort gives one
    // sequence on every backend): each run is one edge's tets in ascending
    // tet order -- what the former hash map's per-edge vector held -- and the
    // runs come in ascending key order, the order that map's keys were
    // sorted into for a deterministic sweep.
    static const int kEdges[6][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
    struct EdgeRec {
        std::array<std::int64_t, 2> mKey;
        std::size_t mTet;
    };
    std::vector<EdgeRec> recs;
    recs.reserve(nt * 6);
    for (std::size_t t = 0; t < nt; ++t)
        for (auto& ev : kEdges)
            recs.push_back({tet_flip_edge_key(rTets[t][ev[0]], rTets[t][ev[1]]), t});
    parallel_sort(recs.begin(), recs.end(), [](const EdgeRec& a, const EdgeRec& b) {
        return a.mKey != b.mKey ? a.mKey < b.mKey : a.mTet < b.mTet;
    });

    std::vector<TetFlipTet> new_tets;
    std::int64_t applied = 0;
    for (std::size_t r0 = 0, r1 = 0; r0 < recs.size(); r0 = r1) {
        for (r1 = r0 + 1; r1 < recs.size() && recs[r1].mKey == recs[r0].mKey; ++r1) {
        }
        const std::int64_t uu = recs[r0].mKey[0], vv = recs[r0].mKey[1];
        // the alive tets on this edge
        std::array<std::size_t, 3> ring_tets{};
        int nr = 0;
        bool too_many = false;
        for (std::size_t r = r0; r < r1; ++r) {
            const std::size_t t = recs[r].mTet;
            if (!alive[t])
                continue;
            if (nr >= 3) {
                too_many = true;
                break;
            }
            ring_tets[nr++] = t;
        }
        if (too_many || nr != 3)
            continue;
        // ring vertices: nodes != u,v across the three tets; must be exactly
        // three distinct, each appearing in exactly two tets.
        std::array<std::int64_t, 3> others{};
        std::array<int, 3> ocount{0, 0, 0};
        int no = 0;
        bool bad = false;
        for (int r = 0; r < 3 && !bad; ++r)
            for (int k = 0; k < 4; ++k) {
                const std::int64_t nd = rTets[ring_tets[r]][k];
                if (nd == uu || nd == vv)
                    continue;
                int idx = -1;
                for (int o = 0; o < no; ++o)
                    if (others[o] == nd)
                        idx = o;
                if (idx < 0) {
                    if (no >= 3) {
                        bad = true;
                        break;
                    }
                    others[no] = nd;
                    ocount[no] = 1;
                    ++no;
                } else {
                    ++ocount[idx];
                }
            }
        if (bad || no != 3 || ocount[0] != 2 || ocount[1] != 2 || ocount[2] != 2)
            continue;

        std::int64_t a = others[0], b = others[1], c = others[2];
        const Vec3 pu = tet_flip_point(rXyz, uu), pv = tet_flip_point(rXyz, vv);
        Vec3 pa = tet_flip_point(rXyz, a), pb = tet_flip_point(rXyz, b),
             pc = tet_flip_point(rXyz, c);
        if (detail::tet_volume6(pa, pb, pc, pu) < 0.0) {
            std::swap(b, c);
            std::swap(pb, pc);
        }
        // u on the positive side of (a,b,c); a valid removal needs v strictly on
        // the far side, and the two caps must tile the three old tets.
        if (detail::tet_volume6(pa, pb, pc, pv) >= -eps_vol6)
            continue;
        const double vcap1 = detail::tet_volume6(pa, pb, pc, pu);  // > 0
        const double vcap2 =
            detail::tet_volume6(pa, pc, pb, pv);  // stored orientation, should be > 0
        if (vcap1 <= eps_vol6 || vcap2 <= eps_vol6)
            continue;
        double old_sum = 0.0;
        for (int r = 0; r < 3; ++r)
            old_sum += std::fabs(detail::tet_volume6(tet_flip_point(rXyz, rTets[ring_tets[r]][0]),
                                                     tet_flip_point(rXyz, rTets[ring_tets[r]][1]),
                                                     tet_flip_point(rXyz, rTets[ring_tets[r]][2]),
                                                     tet_flip_point(rXyz, rTets[ring_tets[r]][3])));
        const double new_sum = vcap1 + vcap2;
        if (std::fabs(new_sum - old_sum) > 1e-9 * old_sum + eps_vol6)
            continue;  // not a clean tiling -> non-convex / invalid
        const TetFlipTet cap1{a, b, c, uu}, cap2{a, c, b, vv};
        const double q_new = std::min(tet_flip_quality_tet(rXyz, cap1, eps_len),
                                      tet_flip_quality_tet(rXyz, cap2, eps_len));
        double q_old = 1.0;
        for (int r = 0; r < 3; ++r)
            q_old = std::min(q_old, tet_flip_quality_tet(rXyz, rTets[ring_tets[r]], eps_len));
        if (q_new <= q_old + min_improve)
            continue;
        for (int r = 0; r < 3; ++r)
            alive[ring_tets[r]] = 0;
        new_tets.push_back(cap1);
        new_tets.push_back(cap2);
        ++applied;
    }

    std::vector<TetFlipTet> out;
    out.reserve(nt);
    for (std::size_t t = 0; t < nt; ++t)
        if (alive[t])
            out.push_back(rTets[t]);
    for (const TetFlipTet& t : new_tets)
        out.push_back(t);
    rTets.swap(out);
    return applied;
}

// --- the worst-cell objective (OptimizeVolumeObjective::Worst) ---------------
// The objective is the least dihedral angle. Every move below is accepted only
// when it raises the least angle of the tets it touches, so the least angle of
// the mesh never drops: the relocation pass is guarded vertex by vertex, the
// topological moves compare the worst new tet with the worst old one, and the
// worst-cell smoothing accepts a step only when its vertex's star improves.

constexpr double kTetFlipRadToDeg = 57.295779513082320876798;

// The least dihedral angle of a tet, in degrees, whatever its orientation; 0
// for a degenerate tet. The shared kernel (detail/fem_kernels.hpp), moved
// there verbatim so adapt_volume measures exactly what this operation does.
inline double tet_flip_min_dihedral(const Vec3& rA, const Vec3& rB, const Vec3& rC,
                                    const Vec3& rD) {
    return detail::tet_min_dihedral(rA, rB, rC, rD);
}

double tet_flip_min_dihedral_tet(const std::vector<double>& rXyz, const TetFlipTet& t) {
    return tet_flip_min_dihedral(tet_flip_point(rXyz, t[0]), tet_flip_point(rXyz, t[1]),
                                 tet_flip_point(rXyz, t[2]), tet_flip_point(rXyz, t[3]));
}

// The worst-cell quality of a tet as stored: its least dihedral angle when it
// is positively oriented by more than the volume tolerance, else -1, so an
// inverted or flat tet is never accepted and never mistaken for a good one.
double tet_flip_worst_q(const std::vector<double>& rXyz, const TetFlipTet& t, double eps_vol6) {
    const Vec3 a = tet_flip_point(rXyz, t[0]), b = tet_flip_point(rXyz, t[1]),
               c = tet_flip_point(rXyz, t[2]), d = tet_flip_point(rXyz, t[3]);
    if (detail::tet_volume6(a, b, c, d) <= eps_vol6)
        return -1.0;
    return tet_flip_min_dihedral(a, b, c, d);
}

// The tets around each vertex (CSR, ascending tet order).
TetFlipStar tet_flip_star(std::size_t N, const std::vector<TetFlipTet>& rTets) {
    TetFlipStar s;
    s.mStart.assign(N + 1, 0);
    for (const TetFlipTet& t : rTets)
        for (std::int64_t v : t)
            ++s.mStart[static_cast<std::size_t>(v) + 1];
    for (std::size_t i = 0; i < N; ++i)
        s.mStart[i + 1] += s.mStart[i];
    s.mTets.resize(static_cast<std::size_t>(s.mStart[N]));
    std::vector<std::int64_t> fill(s.mStart.begin(), s.mStart.end() - 1);
    for (std::size_t t = 0; t < rTets.size(); ++t)
        for (std::int64_t v : rTets[t])
            s.mTets[static_cast<std::size_t>(fill[static_cast<std::size_t>(v)]++)] =
                static_cast<std::int64_t>(t);
    return s;
}

double tet_flip_star_min(const std::vector<double>& rXyz, const std::vector<TetFlipTet>& rTets,
                         const TetFlipStar& rStar, std::size_t V, double eps_vol6) {
    double m = 180.0;
    for (std::int64_t k = rStar.mStart[V]; k < rStar.mStart[V + 1]; ++k)
        m = std::min(
            m, tet_flip_worst_q(
                   rXyz, rTets[static_cast<std::size_t>(rStar.mTets[static_cast<std::size_t>(k)])],
                   eps_vol6));
    return m;
}

// One ODT pass, then the guard: a vertex whose star's least angle fell below
// what it was before the pass goes back, in ascending order, until no moved
// vertex's star is worse than before. Every tet that changed has a moved
// vertex, so no tet ends below the least angle the pass started from.
std::int64_t tet_flip_guarded_odt(std::vector<double>& rXyz, const std::vector<TetFlipTet>& rTets,
                                  const std::vector<std::uint8_t>& rPinned,
                                  const SmoothOptions& rOdt, double eps_vol6) {
    const std::size_t n = rXyz.size() / 3;
    const std::vector<double> before = rXyz;
    detail::smooth_odt_pass(rXyz, rTets, rPinned, rOdt);
    std::vector<std::uint8_t> moved(n, 0);
    std::vector<std::int64_t> work;
    for (std::size_t i = 0; i < n; ++i)
        if (rXyz[3 * i] != before[3 * i] || rXyz[3 * i + 1] != before[3 * i + 1] ||
            rXyz[3 * i + 2] != before[3 * i + 2]) {
            moved[i] = 1;
            work.push_back(static_cast<std::int64_t>(i));
        }
    if (work.empty())
        return 0;
    const TetFlipStar star = tet_flip_star(n, rTets);
    std::vector<double> floor(n, 0.0);
    parallel_for(work.size(), [&](std::size_t k) {
        const auto v = static_cast<std::size_t>(work[k]);
        floor[v] = tet_flip_star_min(before, rTets, star, v, eps_vol6);
    });
    while (!work.empty()) {
        std::vector<std::int64_t> next;
        for (const std::int64_t v64 : work) {
            const auto v = static_cast<std::size_t>(v64);
            if (!moved[v] || tet_flip_star_min(rXyz, rTets, star, v, eps_vol6) >= floor[v])
                continue;
            for (int d = 0; d < 3; ++d)
                rXyz[3 * v + static_cast<std::size_t>(d)] =
                    before[3 * v + static_cast<std::size_t>(d)];
            moved[v] = 0;
            for (std::int64_t k = star.mStart[v]; k < star.mStart[v + 1]; ++k)
                for (const std::int64_t w :
                     rTets[static_cast<std::size_t>(star.mTets[static_cast<std::size_t>(k)])])
                    if (moved[static_cast<std::size_t>(w)])
                        next.push_back(w);
        }
        std::sort(next.begin(), next.end());
        next.erase(std::unique(next.begin(), next.end()), next.end());
        work.swap(next);
    }
    std::int64_t count = 0;
    for (std::size_t i = 0; i < n; ++i)
        count += moved[i];
    return count;
}

// Edge removal (Klingner and Shewchuk 2007): an interior edge (u, v) whose k
// tets form a closed ring r_0 .. r_{k-1} is replaced by the 2 (k - 2) tets
// over the triangulation of the ring that maximises their least angle, found
// by the dynamic program over the triangulations of a polygon. Each triangle
// (r_i, r_j, r_l) gives the tets (r_i, r_j, r_l, u) and (r_i, r_l, r_j, v),
// both of which must be positively oriented. k = 3 is the 3-2 flip.
std::int64_t tet_flip_pass_edge_removal(std::vector<double>& rXyz, std::vector<TetFlipTet>& rTets,
                                        double min_improve, double eps_vol6, int max_ring,
                                        double SkipAbove, std::int64_t& rNum32) {
    const std::size_t nt = rTets.size();
    std::vector<std::uint8_t> alive(nt, 1);
    static const int kEdges[6][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
    struct EdgeRec {
        std::array<std::int64_t, 2> mKey;
        std::size_t mTet;
    };
    std::vector<EdgeRec> recs;
    recs.reserve(nt * 6);
    for (std::size_t t = 0; t < nt; ++t)
        for (auto& ev : kEdges)
            recs.push_back({tet_flip_edge_key(rTets[t][ev[0]], rTets[t][ev[1]]), t});
    parallel_sort(recs.begin(), recs.end(), [](const EdgeRec& a, const EdgeRec& b) {
        return a.mKey != b.mKey ? a.mKey < b.mKey : a.mTet < b.mTet;
    });

    std::vector<TetFlipTet> new_tets;
    std::int64_t applied = 0;
    std::vector<std::size_t> ring_tets;
    std::vector<std::int64_t> ring;
    for (std::size_t r0 = 0, r1 = 0; r0 < recs.size(); r0 = r1) {
        for (r1 = r0 + 1; r1 < recs.size() && recs[r1].mKey == recs[r0].mKey; ++r1) {
        }
        const std::int64_t uu = recs[r0].mKey[0], vv = recs[r0].mKey[1];
        ring_tets.clear();
        bool dead = false;
        for (std::size_t r = r0; r < r1; ++r) {
            if (!alive[recs[r].mTet])
                dead = true;
            ring_tets.push_back(recs[r].mTet);
        }
        const int k = static_cast<int>(ring_tets.size());
        if (dead || k < 3 || k > max_ring)
            continue;
        double q_old = 180.0;
        for (const std::size_t t : ring_tets)
            q_old = std::min(q_old, tet_flip_worst_q(rXyz, rTets[t], eps_vol6));
        if (q_old >= SkipAbove)
            continue;  // every tet around the edge is already good enough
        // Each tet holds one ring edge (a, b); chain them into a closed cycle.
        std::vector<std::array<std::int64_t, 2>> link(static_cast<std::size_t>(k));
        for (int i = 0; i < k; ++i) {
            int w = 0;
            for (std::int64_t nd : rTets[ring_tets[static_cast<std::size_t>(i)]])
                if (nd != uu && nd != vv && w < 2)
                    link[static_cast<std::size_t>(i)][static_cast<std::size_t>(w++)] = nd;
        }
        ring.assign(1, link[0][0]);
        std::int64_t prev_tet = 0, cur = link[0][1];
        bool closed = false;
        for (int step = 1; step <= k; ++step) {
            if (cur == ring[0]) {
                closed = step == k;
                break;
            }
            ring.push_back(cur);
            int found = -1;
            for (int i = 0; i < k; ++i)
                if (i != prev_tet && (link[static_cast<std::size_t>(i)][0] == cur ||
                                      link[static_cast<std::size_t>(i)][1] == cur)) {
                    if (found >= 0) {
                        found = -2;  // a ring vertex on three tets: not a simple ring
                        break;
                    }
                    found = i;
                }
            if (found < 0)
                break;
            prev_tet = found;
            const auto& l = link[static_cast<std::size_t>(found)];
            cur = l[0] == cur ? l[1] : l[0];
        }
        if (!closed || static_cast<int>(ring.size()) != k)
            continue;  // a boundary edge (open fan) or a non-manifold one
        // Orient the ring so that (r_i, r_j, r_l, u) is positive for i < j < l.
        const Vec3 pu = tet_flip_point(rXyz, uu), pv = tet_flip_point(rXyz, vv);
        if (detail::tet_volume6(pu, pv, tet_flip_point(rXyz, ring[0]),
                                tet_flip_point(rXyz, ring[1])) > 0.0)
            std::reverse(ring.begin() + 1, ring.end());

        double old_sum = 0.0;
        for (const std::size_t t : ring_tets) {
            const TetFlipTet& o = rTets[t];
            old_sum += std::fabs(
                detail::tet_volume6(tet_flip_point(rXyz, o[0]), tet_flip_point(rXyz, o[1]),
                                    tet_flip_point(rXyz, o[2]), tet_flip_point(rXyz, o[3])));
        }
        // best[i][l]: the best least angle over triangulations of r_i .. r_l.
        const auto K = static_cast<std::size_t>(k);
        constexpr double kNone = -std::numeric_limits<double>::infinity();
        std::vector<double> best(K * K, std::numeric_limits<double>::infinity());
        std::vector<int> pick(K * K, -1);
        const auto tri_q = [&](int i, int j, int l) {
            const TetFlipTet a{ring[static_cast<std::size_t>(i)], ring[static_cast<std::size_t>(j)],
                               ring[static_cast<std::size_t>(l)], uu};
            const TetFlipTet b{ring[static_cast<std::size_t>(i)], ring[static_cast<std::size_t>(l)],
                               ring[static_cast<std::size_t>(j)], vv};
            const double qa = tet_flip_worst_q(rXyz, a, eps_vol6);
            const double qb = tet_flip_worst_q(rXyz, b, eps_vol6);
            return qa < 0.0 || qb < 0.0 ? kNone : std::min(qa, qb);
        };
        for (int len = 2; len < k; ++len)
            for (int i = 0; i + len < k; ++i) {
                const int l = i + len;
                double b = kNone;
                int arg = -1;
                for (int j = i + 1; j < l; ++j) {
                    const double q = std::min(
                        {best[static_cast<std::size_t>(i) * K + static_cast<std::size_t>(j)],
                         best[static_cast<std::size_t>(j) * K + static_cast<std::size_t>(l)],
                         tri_q(i, j, l)});
                    if (q > b) {
                        b = q;
                        arg = j;
                    }
                }
                best[static_cast<std::size_t>(i) * K + static_cast<std::size_t>(l)] = b;
                pick[static_cast<std::size_t>(i) * K + static_cast<std::size_t>(l)] = arg;
            }
        const double q_new = best[K - 1];
        if (!(q_new > q_old + min_improve))
            continue;
        std::vector<TetFlipTet> made;
        std::vector<std::array<int, 2>> todo{{0, k - 1}};
        while (!todo.empty()) {
            const auto [i, l] = todo.back();
            todo.pop_back();
            if (l - i < 2)
                continue;
            const int j = pick[static_cast<std::size_t>(i) * K + static_cast<std::size_t>(l)];
            made.push_back({ring[static_cast<std::size_t>(i)], ring[static_cast<std::size_t>(j)],
                            ring[static_cast<std::size_t>(l)], uu});
            made.push_back({ring[static_cast<std::size_t>(i)], ring[static_cast<std::size_t>(l)],
                            ring[static_cast<std::size_t>(j)], vv});
            todo.push_back({i, j});
            todo.push_back({j, l});
        }
        double new_sum = 0.0;
        for (const TetFlipTet& o : made)
            new_sum += detail::tet_volume6(tet_flip_point(rXyz, o[0]), tet_flip_point(rXyz, o[1]),
                                           tet_flip_point(rXyz, o[2]), tet_flip_point(rXyz, o[3]));
        if (std::fabs(new_sum - old_sum) > 1e-9 * old_sum + eps_vol6)
            continue;  // not a clean tiling of the ring
        for (const std::size_t t : ring_tets)
            alive[t] = 0;
        new_tets.insert(new_tets.end(), made.begin(), made.end());
        ++applied;
        if (k == 3)
            ++rNum32;
    }

    std::vector<TetFlipTet> out;
    out.reserve(nt + new_tets.size());
    for (std::size_t t = 0; t < nt; ++t)
        if (alive[t])
            out.push_back(rTets[t]);
    out.insert(out.end(), new_tets.begin(), new_tets.end());
    rTets.swap(out);
    return applied;
}

// Multi-face removal (Klingner and Shewchuk 2007), the inverse of edge
// removal: the faces sandwiched between two vertices a and b -- each face f
// with the tets (f, a) and (f, b) -- that form a disc bounded by a ring of
// k <= max_ring vertices are replaced by the k tets around the new edge (a, b).
// The disc is grown from an interior face across its edges as long as the
// neighbouring face is sandwiched by the same two vertices; one face alone is
// the 2-3 flip, which its own pass handles.
std::int64_t tet_flip_pass_multi_face(std::vector<double>& rXyz, std::vector<TetFlipTet>& rTets,
                                      double min_improve, double eps_vol6, int max_ring,
                                      double SkipAbove) {
    const std::size_t nt = rTets.size();
    std::vector<std::uint8_t> alive(nt, 1);
    struct FaceRec {
        std::array<std::int64_t, 3> mKey;
        std::size_t mTet;
        std::int64_t mApex;
    };
    std::vector<FaceRec> recs(nt * 4);
    for (std::size_t t = 0; t < nt; ++t)
        for (int lf = 0; lf < 4; ++lf) {
            const auto fo = tet_flip_face_ordered(rTets[t], lf);
            recs[t * 4 + static_cast<std::size_t>(lf)] = {tet_flip_face_key(fo[0], fo[1], fo[2]), t,
                                                          rTets[t][lf]};
        }
    parallel_sort(recs.begin(), recs.end(), [](const FaceRec& a, const FaceRec& b) {
        return a.mKey != b.mKey ? a.mKey < b.mKey : a.mTet < b.mTet;
    });
    std::vector<std::array<std::int64_t, 2>> edges;
    edges.reserve(nt * 6);
    for (const TetFlipTet& t : rTets)
        for (int i = 0; i < 4; ++i)
            for (int j = i + 1; j < 4; ++j)
                edges.push_back(tet_flip_edge_key(t[i], t[j]));
    parallel_sort(edges.begin(), edges.end());
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
    const auto has_edge = [&](std::int64_t a, std::int64_t b) {
        return std::binary_search(edges.begin(), edges.end(), tet_flip_edge_key(a, b));
    };
    // The (tet, apex) pairs of an interior face, or nullptr.
    const auto face_pair = [&](const std::array<std::int64_t, 3>& rKey) -> const FaceRec* {
        auto it = std::lower_bound(
            recs.begin(), recs.end(), rKey,
            [](const FaceRec& r, const std::array<std::int64_t, 3>& k) { return r.mKey < k; });
        if (it == recs.end() || it->mKey != rKey || it + 1 == recs.end() || (it + 1)->mKey != rKey)
            return nullptr;
        if (it + 2 != recs.end() && (it + 2)->mKey == rKey)
            return nullptr;  // non-manifold
        return &*it;
    };
    // The apex across the face (x, y, apex) from the tet `from`, or -1.
    const auto across = [&](std::int64_t x, std::int64_t y, std::int64_t apex,
                            std::size_t from) -> std::pair<std::int64_t, std::size_t> {
        const FaceRec* p = face_pair(tet_flip_face_key(x, y, apex));
        if (p == nullptr)
            return {-1, 0};
        const FaceRec& other = p[0].mTet == from ? p[1] : p[0];
        if (other.mTet == from)
            return {-1, 0};
        return {other.mApex, other.mTet};
    };

    std::vector<TetFlipTet> new_tets;
    std::int64_t applied = 0;
    for (std::size_t r = 0; r + 1 < recs.size(); ++r) {
        if (recs[r].mKey != recs[r + 1].mKey || (r > 0 && recs[r - 1].mKey == recs[r].mKey))
            continue;  // the first record of an interior face only
        if (r + 2 < recs.size() && recs[r + 2].mKey == recs[r].mKey)
            continue;
        const std::size_t ta = recs[r].mTet, tb = recs[r + 1].mTet;
        if (!alive[ta] || !alive[tb])
            continue;
        const std::int64_t a = recs[r].mApex, b = recs[r + 1].mApex;
        if (a == b ||
            std::min(tet_flip_worst_q(rXyz, rTets[ta], eps_vol6),
                     tet_flip_worst_q(rXyz, rTets[tb], eps_vol6)) >= SkipAbove ||
            has_edge(a, b))
            continue;
        // Grow the sandwiched disc.
        struct SFace {
            std::array<std::int64_t, 3> mV;
            std::size_t mTa, mTb;
        };
        std::vector<SFace> disc{{recs[r].mKey, ta, tb}};
        std::vector<std::array<std::int64_t, 3>> seen{recs[r].mKey};
        bool too_big = false;
        for (std::size_t q = 0; q < disc.size() && !too_big; ++q) {
            const SFace f = disc[q];
            for (int e = 0; e < 3 && !too_big; ++e) {
                const std::int64_t x = f.mV[static_cast<std::size_t>(e)];
                const std::int64_t y = f.mV[static_cast<std::size_t>((e + 1) % 3)];
                const auto [wa, na] = across(x, y, a, f.mTa);
                const auto [wb, nb] = across(x, y, b, f.mTb);
                if (wa < 0 || wa != wb || wa == a || wa == b || !alive[na] || !alive[nb])
                    continue;
                const auto key = tet_flip_face_key(x, y, wa);
                if (std::find(seen.begin(), seen.end(), key) != seen.end())
                    continue;
                seen.push_back(key);
                disc.push_back({key, na, nb});
                if (static_cast<int>(disc.size()) + 2 > max_ring)
                    too_big = true;
            }
        }
        if (too_big || disc.size() < 2)
            continue;
        // The disc's boundary: the edges one face of it uses, chained into a ring.
        std::vector<std::array<std::int64_t, 2>> bnd;
        for (const SFace& f : disc)
            for (int e = 0; e < 3; ++e)
                bnd.push_back(tet_flip_edge_key(f.mV[static_cast<std::size_t>(e)],
                                                f.mV[static_cast<std::size_t>((e + 1) % 3)]));
        std::sort(bnd.begin(), bnd.end());
        std::vector<std::array<std::int64_t, 2>> once;
        for (std::size_t i = 0; i < bnd.size();) {
            std::size_t j = i + 1;
            while (j < bnd.size() && bnd[j] == bnd[i])
                ++j;
            if (j - i == 1)
                once.push_back(bnd[i]);
            i = j;
        }
        const std::size_t k = once.size();
        if (k != disc.size() + 2)
            continue;  // not a disc
        std::vector<std::int64_t> ring{once[0][0], once[0][1]};
        std::vector<std::uint8_t> used(k, 0);
        used[0] = 1;
        bool ok = true;
        while (ring.size() < k && ok) {
            ok = false;
            for (std::size_t i = 0; i < k; ++i)
                if (!used[i] && (once[i][0] == ring.back() || once[i][1] == ring.back())) {
                    const std::int64_t nxt = once[i][0] == ring.back() ? once[i][1] : once[i][0];
                    if (std::find(ring.begin(), ring.end(), nxt) != ring.end())
                        break;
                    used[i] = 1;
                    ring.push_back(nxt);
                    ok = true;
                    break;
                }
        }
        if (ring.size() != k)
            continue;
        const Vec3 pa = tet_flip_point(rXyz, a), pb = tet_flip_point(rXyz, b);
        const double sign = detail::tet_volume6(pa, pb, tet_flip_point(rXyz, ring[0]),
                                                tet_flip_point(rXyz, ring[1])) > 0.0
                                ? 1.0
                                : -1.0;
        std::vector<TetFlipTet> made;
        double new_sum = 0.0;
        bool valid = true;
        for (std::size_t i = 0; i < k && valid; ++i) {
            TetFlipTet t{a, b, ring[i], ring[(i + 1) % k]};
            const double s6 =
                detail::tet_volume6(pa, pb, tet_flip_point(rXyz, t[2]), tet_flip_point(rXyz, t[3]));
            if (s6 * sign <= eps_vol6)
                valid = false;
            if (sign < 0.0)
                std::swap(t[2], t[3]);
            new_sum += std::fabs(s6);
            made.push_back(t);
        }
        if (!valid)
            continue;
        double q_old = 180.0, old_sum = 0.0;
        for (const SFace& f : disc)
            for (const std::size_t t : {f.mTa, f.mTb}) {
                q_old = std::min(q_old, tet_flip_worst_q(rXyz, rTets[t], eps_vol6));
                const TetFlipTet& o = rTets[t];
                old_sum += std::fabs(
                    detail::tet_volume6(tet_flip_point(rXyz, o[0]), tet_flip_point(rXyz, o[1]),
                                        tet_flip_point(rXyz, o[2]), tet_flip_point(rXyz, o[3])));
            }
        if (std::fabs(new_sum - old_sum) > 1e-9 * old_sum + eps_vol6)
            continue;
        double q_new = 180.0;
        for (const TetFlipTet& t : made)
            q_new = std::min(q_new, tet_flip_worst_q(rXyz, t, eps_vol6));
        if (!(q_new > q_old + min_improve))
            continue;
        for (const SFace& f : disc) {
            alive[f.mTa] = 0;
            alive[f.mTb] = 0;
        }
        new_tets.insert(new_tets.end(), made.begin(), made.end());
        ++applied;
    }

    std::vector<TetFlipTet> out;
    out.reserve(nt + new_tets.size());
    for (std::size_t t = 0; t < nt; ++t)
        if (alive[t])
            out.push_back(rTets[t]);
    out.insert(out.end(), new_tets.begin(), new_tets.end());
    rTets.swap(out);
    return applied;
}

// The six dihedral angles of a tet in a fixed edge order, in degrees.
std::array<double, 6> tet_flip_dihedrals(const Vec3& rA, const Vec3& rB, const Vec3& rC,
                                         const Vec3& rD) {
    const Vec3 p[4] = {rA, rB, rC, rD};
    constexpr int kEdge[6][4] = {{0, 1, 2, 3}, {0, 2, 1, 3}, {0, 3, 1, 2},
                                 {1, 2, 0, 3}, {1, 3, 0, 2}, {2, 3, 0, 1}};
    std::array<double, 6> out{};
    for (std::size_t e = 0; e < 6; ++e) {
        const Vec3& a = p[kEdge[e][0]];
        const Vec3 ab = detail::vec3_sub(p[kEdge[e][1]], a);
        const Vec3 n1 = detail::vec3_cross(ab, detail::vec3_sub(p[kEdge[e][2]], a));
        const Vec3 n2 = detail::vec3_cross(ab, detail::vec3_sub(p[kEdge[e][3]], a));
        const double den = std::sqrt(detail::vec3_dot(n1, n1) * detail::vec3_dot(n2, n2));
        out[e] = den > 0.0 ? std::acos(std::clamp(detail::vec3_dot(n1, n2) / den, -1.0, 1.0)) *
                                 kTetFlipRadToDeg
                           : 0.0;
    }
    return out;
}

// Worst-cell smoothing (Freitag and Ollivier-Gooch 1997): a non-smooth ascent
// on the least dihedral angle of one vertex's star. The search direction is
// the shortest vector in the convex hull of the gradients of the angles within
// a small band of the least one (central differences); a backtracking line
// search accepts a step only when the star's least angle rises. Vertices of
// tets below `Threshold` degrees are visited worst first.
std::int64_t tet_flip_smooth_worst(std::vector<double>& rXyz, const std::vector<TetFlipTet>& rTets,
                                   const std::vector<std::uint8_t>& rFixed, double eps_vol6,
                                   double Threshold) {
    const std::size_t n = rXyz.size() / 3;
    const TetFlipStar star = tet_flip_star(n, rTets);
    std::vector<double> q(rTets.size());
    parallel_for(rTets.size(),
                 [&](std::size_t t) { q[t] = tet_flip_worst_q(rXyz, rTets[t], eps_vol6); });
    std::vector<std::pair<double, std::int64_t>> order;
    for (std::size_t v = 0; v < n; ++v) {
        if (rFixed[v] || star.mStart[v] == star.mStart[v + 1])
            continue;
        double m = 180.0;
        for (std::int64_t k = star.mStart[v]; k < star.mStart[v + 1]; ++k)
            m = std::min(m, q[static_cast<std::size_t>(star.mTets[static_cast<std::size_t>(k)])]);
        if (m < Threshold && m >= 0.0)
            order.push_back({m, static_cast<std::int64_t>(v)});
    }
    std::sort(order.begin(), order.end());

    std::int64_t improved = 0;
    for (const auto& [m0, v64] : order) {
        (void)m0;
        const auto v = static_cast<std::size_t>(v64);
        // The star's tets with v as corner 0, so a move of v moves corner 0.
        std::vector<TetFlipTet> local;
        double shortest = std::numeric_limits<double>::infinity();
        for (std::int64_t k = star.mStart[v]; k < star.mStart[v + 1]; ++k) {
            const TetFlipTet& t =
                rTets[static_cast<std::size_t>(star.mTets[static_cast<std::size_t>(k)])];
            local.push_back(t);
            for (std::int64_t w : t)
                if (w != v64) {
                    const Vec3 d =
                        detail::vec3_sub(tet_flip_point(rXyz, w), tet_flip_point(rXyz, v64));
                    shortest = std::min(shortest, std::sqrt(detail::vec3_dot(d, d)));
                }
        }
        const auto f_at = [&](const Vec3& x) {
            double m = 180.0;
            for (const TetFlipTet& t : local) {
                Vec3 p[4];
                for (int c = 0; c < 4; ++c)
                    p[c] = t[static_cast<std::size_t>(c)] == v64
                               ? x
                               : tet_flip_point(rXyz, t[static_cast<std::size_t>(c)]);
                if (detail::tet_volume6(p[0], p[1], p[2], p[3]) <= eps_vol6)
                    return -1.0;
                m = std::min(m, tet_flip_min_dihedral(p[0], p[1], p[2], p[3]));
            }
            return m;
        };
        // Every angle of the star as a function of the vertex position.
        const auto angles_at = [&](const Vec3& x, std::vector<double>& rOut) {
            rOut.clear();
            for (const TetFlipTet& t : local) {
                Vec3 p[4];
                for (int c = 0; c < 4; ++c)
                    p[c] = t[static_cast<std::size_t>(c)] == v64
                               ? x
                               : tet_flip_point(rXyz, t[static_cast<std::size_t>(c)]);
                const auto a = tet_flip_dihedrals(p[0], p[1], p[2], p[3]);
                rOut.insert(rOut.end(), a.begin(), a.end());
            }
        };
        Vec3 x = tet_flip_point(rXyz, v64);
        const Vec3 x0 = x;
        double fx = f_at(x);
        if (fx < 0.0 || !std::isfinite(shortest) || !(shortest > 0.0))
            continue;
        const double h = 1e-6 * shortest;
        std::vector<double> base;
        for (int iter = 0; iter < 8; ++iter) {
            angles_at(x, base);
            // The active set: angles within a band of the least; only their
            // tets are differentiated (central differences).
            const double band = 1e-3 + 1e-6 * fx;
            std::vector<std::array<double, 3>> act;
            for (std::size_t a = 0; a < base.size(); ++a) {
                if (base[a] > fx + band)
                    continue;
                const TetFlipTet& t = local[a / 6];
                std::array<double, 3> g{};
                for (int d = 0; d < 3; ++d) {
                    Vec3 xp = x, xm = x;
                    xp[static_cast<std::size_t>(d)] += h;
                    xm[static_cast<std::size_t>(d)] -= h;
                    Vec3 pp[4], pm[4];
                    for (int c = 0; c < 4; ++c) {
                        const std::int64_t w = t[static_cast<std::size_t>(c)];
                        pp[c] = w == v64 ? xp : tet_flip_point(rXyz, w);
                        pm[c] = w == v64 ? xm : tet_flip_point(rXyz, w);
                    }
                    const auto ap = tet_flip_dihedrals(pp[0], pp[1], pp[2], pp[3]);
                    const auto am = tet_flip_dihedrals(pm[0], pm[1], pm[2], pm[3]);
                    g[static_cast<std::size_t>(d)] = (ap[a % 6] - am[a % 6]) / (2.0 * h);
                }
                act.push_back(g);
            }
            if (act.empty())
                break;
            // The shortest vector in their convex hull (Frank-Wolfe).
            std::array<double, 3> g = act[0];
            for (int it = 0; it < 64 && act.size() > 1; ++it) {
                std::size_t best_i = 0;
                double best_dot = std::numeric_limits<double>::infinity();
                for (std::size_t i = 0; i < act.size(); ++i) {
                    const double dt = g[0] * act[i][0] + g[1] * act[i][1] + g[2] * act[i][2];
                    if (dt < best_dot) {
                        best_dot = dt;
                        best_i = i;
                    }
                }
                const std::array<double, 3> s{act[best_i][0] - g[0], act[best_i][1] - g[1],
                                              act[best_i][2] - g[2]};
                const double ss = s[0] * s[0] + s[1] * s[1] + s[2] * s[2];
                const double gs = -(g[0] * s[0] + g[1] * s[1] + g[2] * s[2]);
                if (ss <= 0.0 || gs <= 1e-15 * ss)
                    break;
                const double lam = std::min(1.0, gs / ss);
                for (int c = 0; c < 3; ++c)
                    g[static_cast<std::size_t>(c)] += lam * s[static_cast<std::size_t>(c)];
            }
            const double gn = std::sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]);
            if (!(gn > 0.0) || !std::isfinite(gn))
                break;  // stationary: no direction raises every active angle
            bool stepped = false;
            for (int k = 0; k < 14; ++k) {
                const double alpha = 0.2 * shortest * std::ldexp(1.0, -k);
                const Vec3 y{x[0] + alpha * g[0] / gn, x[1] + alpha * g[1] / gn,
                             x[2] + alpha * g[2] / gn};
                const double fy = f_at(y);
                if (fy > fx + 1e-9) {
                    x = y;
                    fx = fy;
                    stepped = true;
                    break;
                }
            }
            if (!stepped)
                break;
        }
        if (x != x0) {
            for (int d = 0; d < 3; ++d)
                rXyz[3 * v + static_cast<std::size_t>(d)] = x[static_cast<std::size_t>(d)];
            ++improved;
        }
    }
    return improved;
}

// The boundary vertices: those on a face only one tet has.
std::vector<std::uint8_t> tet_flip_boundary_mask(std::size_t N,
                                                 const std::vector<TetFlipTet>& rTets) {
    std::vector<std::array<std::int64_t, 3>> faces;
    faces.reserve(rTets.size() * 4);
    for (const TetFlipTet& t : rTets)
        for (int lf = 0; lf < 4; ++lf) {
            const auto fo = tet_flip_face_ordered(t, lf);
            faces.push_back(tet_flip_face_key(fo[0], fo[1], fo[2]));
        }
    parallel_sort(faces.begin(), faces.end());
    std::vector<std::uint8_t> on(N, 0);
    for (std::size_t i = 0; i < faces.size();) {
        std::size_t j = i + 1;
        while (j < faces.size() && faces[j] == faces[i])
            ++j;
        if (j - i == 1)
            for (std::int64_t v : faces[i])
                on[static_cast<std::size_t>(v)] = 1;
        i = j;
    }
    return on;
}

}  // namespace detail
}  // namespace meshioplusplus

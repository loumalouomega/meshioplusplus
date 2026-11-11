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

#pragma once

/**
 * @file detail/tet_flips.hpp
 * @brief The topological moves and guarded smoothing of `optimize_volume`,
 * shared with `adapt_volume`: the 2-3 and 3-2 flips, edge removal and
 * multi-face removal (Shewchuk 2002; Klingner and Shewchuk 2007), the guarded
 * ODT relocation and the worst-cell smoothing (Freitag and Ollivier-Gooch
 * 1997), over a flat xyz buffer and a vector of tetrahedra.
 *
 * Moved verbatim from `optimize_volume.cpp` (only the names changed: the
 * `optvol_` prefix is `tet_flip_`), so `optimize_volume`'s golden digests and
 * the determinism baselines pin it. A move touches interior faces and edges
 * only: a face with one tetrahedron in the vector given is never flipped, so a
 * caller keeps a face fixed by leaving one side of it out of the vector.
 *
 * A **core-private** header: nothing installed names it.
 */

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
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/operations/smooth.hpp"
#include "meshioplusplus/parallel.hpp"

// Project includes (private, not installed)
#include "fem_kernels.hpp"

namespace meshioplusplus {
namespace detail {

// --- tiny geometry over a flat xyz buffer ------------------------------------

inline Vec3 tet_flip_point(const std::vector<double>& rXyz, std::int64_t v) {
    const std::size_t o = static_cast<std::size_t>(v) * 3;
    return {rXyz[o], rXyz[o + 1], rXyz[o + 2]};
}

// Six times the signed volume of a tetra (`detail::tet_volume6`) is the sign
// of its orientation and 6*|volume|; it serves both the convexity/validity
// tests and the tiling checks (no division, no in-sphere test anywhere). The
// per-corner factor of the scaled Jacobian is `detail::det_unit3`, the kernel
// `quality` reports with.

// Orientation-normalised scaled Jacobian in [0, 1]: the tet is measured on its
// positively-oriented reordering, so a negatively-stored but well-shaped tet
// reads as good and a degenerate one as 0. This is exactly the metric
// `compute_quality`/`quality` reports for a tetra (`sqrt(2) * min corner det of
// unit outgoing edges`), so the before/after numbers this operation reports
// line up with what a user sees from `quality`.
inline double tet_flip_quality(const Vec3& q0, const Vec3& q1, const Vec3& q2, const Vec3& q3,
                               double eps_len) {
    Vec3 p0 = q0, p1 = q1, p2 = q2, p3 = q3;
    if (detail::tet_volume6(p0, p1, p2, p3) < 0.0)
        std::swap(p2, p3);  // reorder to positive orientation
    const double j0 = detail::det_unit3(detail::vec3_sub(p1, p0), detail::vec3_sub(p2, p0),
                                        detail::vec3_sub(p3, p0), eps_len);
    const double j1 = detail::det_unit3(detail::vec3_sub(p2, p1), detail::vec3_sub(p0, p1),
                                        detail::vec3_sub(p3, p1), eps_len);
    const double j2 = detail::det_unit3(detail::vec3_sub(p0, p2), detail::vec3_sub(p1, p2),
                                        detail::vec3_sub(p3, p2), eps_len);
    const double j3 = detail::det_unit3(detail::vec3_sub(p0, p3), detail::vec3_sub(p2, p3),
                                        detail::vec3_sub(p1, p3), eps_len);
    if (!(std::isfinite(j0) && std::isfinite(j1) && std::isfinite(j2) && std::isfinite(j3)))
        return 0.0;
    const double sj = std::sqrt(2.0) * std::min({j0, j1, j2, j3});
    return std::isfinite(sj) && sj > 0.0 ? sj : 0.0;
}

/// A tetrahedron as its four point ids (the same layout as `TetCorners`).
using TetFlipTet = std::array<std::int64_t, 4>;

inline double tet_flip_quality_tet(const std::vector<double>& rXyz, const TetFlipTet& t,
                                   double eps_len) {
    return tet_flip_quality(tet_flip_point(rXyz, t[0]), tet_flip_point(rXyz, t[1]),
                            tet_flip_point(rXyz, t[2]), tet_flip_point(rXyz, t[3]), eps_len);
}

// Store a tet in positive orientation (swap the last two corners if needed).
inline TetFlipTet tet_flip_orient_positive(const std::vector<double>& rXyz, TetFlipTet t) {
    if (detail::tet_volume6(tet_flip_point(rXyz, t[0]), tet_flip_point(rXyz, t[1]),
                            tet_flip_point(rXyz, t[2]), tet_flip_point(rXyz, t[3])) < 0.0)
        std::swap(t[2], t[3]);
    return t;
}

inline std::array<std::int64_t, 3> tet_flip_face_key(std::int64_t a, std::int64_t b,
                                                     std::int64_t c) {
    std::array<std::int64_t, 3> f{a, b, c};
    std::sort(f.begin(), f.end());
    return f;
}
inline std::array<std::int64_t, 2> tet_flip_edge_key(std::int64_t a, std::int64_t b) {
    return a < b ? std::array<std::int64_t, 2>{a, b} : std::array<std::int64_t, 2>{b, a};
}

// The three corners of tet `t` other than local corner `lf`, IN t's own order
// (so their winding relative to the apex t[lf] is preserved).
inline std::array<std::int64_t, 3> tet_flip_face_ordered(const TetFlipTet& t, int lf) {
    std::array<std::int64_t, 3> f{};
    int w = 0;
    for (int k = 0; k < 4; ++k)
        if (k != lf)
            f[w++] = t[k];
    return f;
}

// --- 2-3 sub-pass ------------------------------------------------------------
// Every interior triangular face shared by two tets is a candidate: replace the
// two tets with three around the segment joining their apexes, iff the segment
// lies inside the (convex) union -- a pure signed-volume test -- and the worst
// of the three new tets is strictly better than the worst of the two old ones.
template <class TQuality>
std::int64_t tet_flip_pass_23(std::vector<double>& rXyz, std::vector<TetFlipTet>& rTets,
                              double min_improve, double eps_vol6, TQuality&& rQuality,
                              double SkipAbove = std::numeric_limits<double>::infinity()) {
    const std::size_t nt = rTets.size();
    std::vector<std::uint8_t> alive(nt, 1);

    // Every (face key, slot) of every tet, slot = 4 t + local face, sorted (a
    // total order, so parallel_sort gives one sequence on every backend):
    // equal faces form runs in slot order -- the order the former hash map
    // saw them in -- so a run gives each face's count and its first two
    // (tet, apex) occurrences, with no hashing. Counts are kept modulo 256
    // and the two recorded occurrences are the last ones at each residue,
    // exactly as the former `std::uint8_t` counter left them.
    struct FaceRec {
        std::array<std::int64_t, 3> mKey;
        std::size_t mSlot;
    };
    struct FaceRun {
        std::uint8_t mCount = 0;
        std::array<std::pair<std::size_t, std::int64_t>, 2> mEntries{};
    };
    std::vector<FaceRec> recs(nt * 4);
    for (std::size_t t = 0; t < nt; ++t)
        for (int lf = 0; lf < 4; ++lf) {
            const auto fo = tet_flip_face_ordered(rTets[t], lf);
            recs[t * 4 + static_cast<std::size_t>(lf)] = {tet_flip_face_key(fo[0], fo[1], fo[2]),
                                                          t * 4 + static_cast<std::size_t>(lf)};
        }
    parallel_sort(recs.begin(), recs.end(), [](const FaceRec& a, const FaceRec& b) {
        return a.mKey != b.mKey ? a.mKey < b.mKey : a.mSlot < b.mSlot;
    });
    std::vector<FaceRun> runs;
    std::vector<std::uint32_t> run_of(nt * 4);
    for (std::size_t i = 0; i < recs.size();) {
        std::size_t j = i;
        FaceRun run;
        for (; j < recs.size() && recs[j].mKey == recs[i].mKey; ++j) {
            const std::size_t slot = recs[j].mSlot;
            if (run.mCount < 2)
                run.mEntries[run.mCount] = {slot / 4, rTets[slot / 4][slot % 4]};
            ++run.mCount;
            run_of[slot] = static_cast<std::uint32_t>(runs.size());
        }
        runs.push_back(run);
        i = j;
    }

    std::vector<TetFlipTet> new_tets;
    std::int64_t applied = 0;
    for (std::size_t t = 0; t < nt; ++t) {
        if (!alive[t])
            continue;
        for (int lf = 0; lf < 4; ++lf) {
            const auto fo = tet_flip_face_ordered(rTets[t], lf);
            const FaceRun& run = runs[run_of[t * 4 + static_cast<std::size_t>(lf)]];
            if (run.mCount != 2)
                continue;  // boundary or non-manifold face
            const auto& pair = run.mEntries;
            std::size_t u = pair[0].first == t ? pair[1].first : pair[0].first;
            std::int64_t e = pair[0].first == t ? pair[1].second : pair[0].second;
            if (u == t || u < t || !alive[u])
                continue;  // process each interior face once, from the lower-id tet
            const std::int64_t d = rTets[t][lf];

            // Orient the face so (f0,f1,f2,d) is positive; then a valid flip
            // needs e strictly on the far side (svol6 negative). The three new
            // tets are (f0,f1,e,d), (f1,f2,e,d), (f2,f0,e,d) -- positively
            // oriented in a convex union, and at least one non-positive in a
            // non-convex one (the predicate-free convexity test). Winding pinned
            // by OptimizeVolume.TwentyThreeFlipWindingIsPositiveOnAConvexUnion.
            std::int64_t f0 = fo[0], f1 = fo[1], f2 = fo[2];
            const Vec3 pf0 = tet_flip_point(rXyz, f0), pf1 = tet_flip_point(rXyz, f1),
                       pf2 = tet_flip_point(rXyz, f2), pd = tet_flip_point(rXyz, d),
                       pe = tet_flip_point(rXyz, e);
            if (detail::tet_volume6(pf0, pf1, pf2, pd) < 0.0) {
                std::swap(f1, f2);
            }
            const Vec3 g0 = tet_flip_point(rXyz, f0), g1 = tet_flip_point(rXyz, f1),
                       g2 = tet_flip_point(rXyz, f2);
            if (detail::tet_volume6(g0, g1, g2, pe) >= -eps_vol6)
                continue;  // apex e not strictly on the opposite side
            const TetFlipTet A{f0, f1, e, d}, B{f1, f2, e, d}, C{f2, f0, e, d};
            const double vA = detail::tet_volume6(g0, g1, pe, pd);
            const double vB = detail::tet_volume6(g1, g2, pe, pd);
            const double vC = detail::tet_volume6(g2, g0, pe, pd);
            if (vA <= eps_vol6 || vB <= eps_vol6 || vC <= eps_vol6)
                continue;  // non-convex union -> invalid flip
            const double q_old = std::min(rQuality(rTets[t]), rQuality(rTets[u]));
            if (q_old >= SkipAbove)
                continue;  // both tets already good enough: not a worst cell
            const double q_new = std::min({rQuality(A), rQuality(B), rQuality(C)});
            if (q_new <= q_old + min_improve)
                continue;
            alive[t] = 0;
            alive[u] = 0;
            new_tets.push_back(A);
            new_tets.push_back(B);
            new_tets.push_back(C);
            ++applied;
            break;  // t is now dead
        }
    }

    // compact: surviving originals (ascending) then this pass's new tets
    std::vector<TetFlipTet> out;
    out.reserve(nt + new_tets.size());
    for (std::size_t t = 0; t < nt; ++t)
        if (alive[t])
            out.push_back(rTets[t]);
    for (const TetFlipTet& t : new_tets)
        out.push_back(t);
    rTets.swap(out);
    return applied;
}

std::int64_t tet_flip_pass_32(std::vector<double>& rXyz, std::vector<TetFlipTet>& rTets,
                              double min_improve, double eps_len, double eps_vol6);

double tet_flip_min_dihedral_tet(const std::vector<double>& rXyz, const TetFlipTet& t);

double tet_flip_worst_q(const std::vector<double>& rXyz, const TetFlipTet& t, double eps_vol6);

struct TetFlipStar {
    std::vector<std::int64_t> mStart;
    std::vector<std::int64_t> mTets;
};

TetFlipStar tet_flip_star(std::size_t N, const std::vector<TetFlipTet>& rTets);

double tet_flip_star_min(const std::vector<double>& rXyz, const std::vector<TetFlipTet>& rTets,
                         const TetFlipStar& rStar, std::size_t V, double eps_vol6);

std::int64_t tet_flip_guarded_odt(std::vector<double>& rXyz, const std::vector<TetFlipTet>& rTets,
                                  const std::vector<std::uint8_t>& rPinned,
                                  const SmoothOptions& rOdt, double eps_vol6);

std::int64_t tet_flip_pass_edge_removal(std::vector<double>& rXyz, std::vector<TetFlipTet>& rTets,
                                        double min_improve, double eps_vol6, int max_ring,
                                        double SkipAbove, std::int64_t& rNum32);

std::int64_t tet_flip_pass_multi_face(std::vector<double>& rXyz, std::vector<TetFlipTet>& rTets,
                                      double min_improve, double eps_vol6, int max_ring,
                                      double SkipAbove);

std::int64_t tet_flip_smooth_worst(std::vector<double>& rXyz, const std::vector<TetFlipTet>& rTets,
                                   const std::vector<std::uint8_t>& rFixed, double eps_vol6,
                                   double Threshold);

std::array<double, 6> tet_flip_dihedrals(const Vec3& rA, const Vec3& rB, const Vec3& rC,
                                         const Vec3& rD);

std::vector<std::uint8_t> tet_flip_boundary_mask(std::size_t N,
                                                 const std::vector<TetFlipTet>& rTets);

}  // namespace detail
}  // namespace meshioplusplus

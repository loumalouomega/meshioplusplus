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
// Predicate-free ODT remeshing of a tetrahedral mesh: alternate ODT vertex
// relocation (reused verbatim from `smooth`'s SmoothMethod::Odt) with
// quality-improving 2-3 / 3-2 topological flips. See
// operations/optimize_volume.hpp for the contract, and doc/optimize_volume.md
// for the algorithm and attribution (Freitag & Ollivier-Gooch's improvement
// rule; no in-sphere/Delaunay predicate).

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/optimize_volume.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/operations/smooth.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {

namespace {

using detail::Vec3;

constexpr const char* kOvPrefix = "meshio++: optimize_volume: ";

// --- tet-only scope check (the smooth ODT / decimate_volume shape) -----------

void optvol_check_blocks(const Mesh& rMesh) {
    bool has_tet = false;
    for (const auto cb : rMesh.CellRange()) {
        if (cb.IsPolyhedron())
            throw std::invalid_argument(
                std::string(kOvPrefix) +
                "operates on tet-only meshes, but the mesh contains a polyhedron cell block -- "
                "run convert_cells(mode='simplexify') first");
        if (cb.IsRagged())
            throw std::invalid_argument(
                std::string(kOvPrefix) +
                "operates on tet-only meshes, but the mesh contains ragged cell block '" +
                std::string(cb.Type()) + "'");
        const CellType ct = cell_type_from_name(cb.Type());
        if (ct == CellType::Tetra) {
            has_tet = true;
            continue;
        }
        const int dim = cell_type_dimension(ct);
        if (dim == 3)
            throw std::invalid_argument(
                std::string(kOvPrefix) +
                "operates on tet-only meshes, but the mesh contains 3D cell block '" +
                std::string(cb.Type()) +
                "' that is not linear tetra -- run convert_cells(mode='simplexify') first");
        throw std::invalid_argument(
            std::string(kOvPrefix) +
            "operates on tet-only meshes, but the mesh contains non-3D cell block '" +
            std::string(cb.Type()) +
            "' alongside its tets (its nodes would dangle; drop it first, e.g. via split)");
    }
    if (!has_tet)
        throw std::invalid_argument(std::string(kOvPrefix) +
                                    "requires at least one tetra cell block");
}

// --- tiny geometry over a flat xyz buffer ------------------------------------

inline Vec3 optvol_point(const std::vector<double>& rXyz, std::int64_t v) {
    const std::size_t o = static_cast<std::size_t>(v) * 3;
    return {rXyz[o], rXyz[o + 1], rXyz[o + 2]};
}

// Six times the signed volume of tetra (p0,p1,p2,p3): triple product of the
// three edges out of p0. Sign is the orientation, magnitude 6*|volume|; used
// for both the convexity/validity tests and the tiling checks (no division,
// no in-sphere test anywhere).
inline double optvol_svol6(const Vec3& p0, const Vec3& p1, const Vec3& p2, const Vec3& p3) {
    return detail::triple_product(detail::vec3_sub(p1, p0), detail::vec3_sub(p2, p0),
                                  detail::vec3_sub(p3, p0));
}

// det of the three UNIT outgoing edges at a corner -- the per-corner factor of
// the scaled Jacobian (quality.cpp's quality_det_unit3, re-derived here rather
// than reached into its anonymous namespace). NaN on a zero-length edge.
inline double optvol_det_unit3(const Vec3& u, const Vec3& v, const Vec3& w, double eps) {
    const double lu = detail::vec3_norm(u), lv = detail::vec3_norm(v), lw = detail::vec3_norm(w);
    if (lu < eps || lv < eps || lw < eps)
        return std::numeric_limits<double>::quiet_NaN();
    return detail::triple_product(detail::vec3_scale(u, 1.0 / lu), detail::vec3_scale(v, 1.0 / lv),
                                  detail::vec3_scale(w, 1.0 / lw));
}

// Orientation-normalised scaled Jacobian in [0, 1]: the tet is measured on its
// positively-oriented reordering, so a negatively-stored but well-shaped tet
// reads as good and a degenerate one as 0. This is exactly the metric
// `compute_quality`/`quality` reports for a tetra (`sqrt(2) * min corner det of
// unit outgoing edges`), so the before/after numbers this operation reports
// line up with what a user sees from `quality`.
double optvol_quality(const Vec3& q0, const Vec3& q1, const Vec3& q2, const Vec3& q3,
                      double eps_len) {
    Vec3 p0 = q0, p1 = q1, p2 = q2, p3 = q3;
    if (optvol_svol6(p0, p1, p2, p3) < 0.0)
        std::swap(p2, p3);  // reorder to positive orientation
    const double j0 = optvol_det_unit3(detail::vec3_sub(p1, p0), detail::vec3_sub(p2, p0),
                                       detail::vec3_sub(p3, p0), eps_len);
    const double j1 = optvol_det_unit3(detail::vec3_sub(p2, p1), detail::vec3_sub(p0, p1),
                                       detail::vec3_sub(p3, p1), eps_len);
    const double j2 = optvol_det_unit3(detail::vec3_sub(p0, p2), detail::vec3_sub(p1, p2),
                                       detail::vec3_sub(p3, p2), eps_len);
    const double j3 = optvol_det_unit3(detail::vec3_sub(p0, p3), detail::vec3_sub(p2, p3),
                                       detail::vec3_sub(p1, p3), eps_len);
    if (!(std::isfinite(j0) && std::isfinite(j1) && std::isfinite(j2) && std::isfinite(j3)))
        return 0.0;
    const double sj = std::sqrt(2.0) * std::min({j0, j1, j2, j3});
    return std::isfinite(sj) && sj > 0.0 ? sj : 0.0;
}

using Tet = std::array<std::int64_t, 4>;

double optvol_quality_tet(const std::vector<double>& rXyz, const Tet& t, double eps_len) {
    return optvol_quality(optvol_point(rXyz, t[0]), optvol_point(rXyz, t[1]),
                          optvol_point(rXyz, t[2]), optvol_point(rXyz, t[3]), eps_len);
}

// Store a tet in positive orientation (swap the last two corners if needed).
inline Tet optvol_orient_positive(const std::vector<double>& rXyz, Tet t) {
    if (optvol_svol6(optvol_point(rXyz, t[0]), optvol_point(rXyz, t[1]), optvol_point(rXyz, t[2]),
                     optvol_point(rXyz, t[3])) < 0.0)
        std::swap(t[2], t[3]);
    return t;
}

// FNV-1a hash for small int64 arrays (face/edge keys).
template <std::size_t K>
struct OptvolArrHash {
    std::size_t operator()(const std::array<std::int64_t, K>& a) const {
        std::size_t h = 1469598103934665603ull;
        for (std::size_t i = 0; i < K; ++i) {
            h ^= static_cast<std::size_t>(a[i]);
            h *= 1099511628211ull;
        }
        return h;
    }
};

std::array<std::int64_t, 3> optvol_face_key(std::int64_t a, std::int64_t b, std::int64_t c) {
    std::array<std::int64_t, 3> f{a, b, c};
    std::sort(f.begin(), f.end());
    return f;
}
std::array<std::int64_t, 2> optvol_edge_key(std::int64_t a, std::int64_t b) {
    return a < b ? std::array<std::int64_t, 2>{a, b} : std::array<std::int64_t, 2>{b, a};
}

// The three corners of tet `t` other than local corner `lf`, IN t's own order
// (so their winding relative to the apex t[lf] is preserved).
std::array<std::int64_t, 3> optvol_face_ordered(const Tet& t, int lf) {
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
std::int64_t optvol_pass_23(std::vector<double>& rXyz, std::vector<Tet>& rTets,
                            double min_improve, double eps_len, double eps_vol6) {
    const std::size_t nt = rTets.size();
    std::vector<std::uint8_t> alive(nt, 1);

    // face key -> up to two (tet index, apex node id)
    std::unordered_map<std::array<std::int64_t, 3>, std::array<std::pair<std::size_t, std::int64_t>, 2>,
                       OptvolArrHash<3>>
        face_map;
    std::unordered_map<std::array<std::int64_t, 3>, std::uint8_t, OptvolArrHash<3>> face_count;
    face_map.reserve(nt * 4);
    face_count.reserve(nt * 4);
    for (std::size_t t = 0; t < nt; ++t)
        for (int lf = 0; lf < 4; ++lf) {
            const auto fo = optvol_face_ordered(rTets[t], lf);
            const auto key = optvol_face_key(fo[0], fo[1], fo[2]);
            std::uint8_t& cnt = face_count[key];
            if (cnt < 2)
                face_map[key][cnt] = {t, rTets[t][lf]};
            ++cnt;
        }

    std::vector<Tet> new_tets;
    std::int64_t applied = 0;
    for (std::size_t t = 0; t < nt; ++t) {
        if (!alive[t])
            continue;
        for (int lf = 0; lf < 4; ++lf) {
            const auto fo = optvol_face_ordered(rTets[t], lf);
            const auto key = optvol_face_key(fo[0], fo[1], fo[2]);
            if (face_count[key] != 2)
                continue;  // boundary or non-manifold face
            const auto& pair = face_map[key];
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
            const Vec3 pf0 = optvol_point(rXyz, f0), pf1 = optvol_point(rXyz, f1),
                       pf2 = optvol_point(rXyz, f2), pd = optvol_point(rXyz, d),
                       pe = optvol_point(rXyz, e);
            if (optvol_svol6(pf0, pf1, pf2, pd) < 0.0) {
                std::swap(f1, f2);
            }
            const Vec3 g0 = optvol_point(rXyz, f0), g1 = optvol_point(rXyz, f1),
                       g2 = optvol_point(rXyz, f2);
            if (optvol_svol6(g0, g1, g2, pe) >= -eps_vol6)
                continue;  // apex e not strictly on the opposite side
            const Tet A{f0, f1, e, d}, B{f1, f2, e, d}, C{f2, f0, e, d};
            const double vA = optvol_svol6(g0, g1, pe, pd);
            const double vB = optvol_svol6(g1, g2, pe, pd);
            const double vC = optvol_svol6(g2, g0, pe, pd);
            if (vA <= eps_vol6 || vB <= eps_vol6 || vC <= eps_vol6)
                continue;  // non-convex union -> invalid flip
            const double q_new = std::min(
                {optvol_quality_tet(rXyz, A, eps_len), optvol_quality_tet(rXyz, B, eps_len),
                 optvol_quality_tet(rXyz, C, eps_len)});
            const double q_old = std::min(optvol_quality_tet(rXyz, rTets[t], eps_len),
                                          optvol_quality_tet(rXyz, rTets[u], eps_len));
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
    std::vector<Tet> out;
    out.reserve(nt + new_tets.size());
    for (std::size_t t = 0; t < nt; ++t)
        if (alive[t])
            out.push_back(rTets[t]);
    for (const Tet& t : new_tets)
        out.push_back(t);
    rTets.swap(out);
    return applied;
}

// --- 3-2 sub-pass ------------------------------------------------------------
// Every interior edge shared by exactly three tets forming a ring is a
// candidate: replace the three tets with two capping the ring triangle, iff the
// two caps tile the same trigonal bipyramid (a signed-volume tiling test) and
// the worse cap beats the worst of the three old tets.
std::int64_t optvol_pass_32(std::vector<double>& rXyz, std::vector<Tet>& rTets,
                            double min_improve, double eps_len, double eps_vol6) {
    const std::size_t nt = rTets.size();
    std::vector<std::uint8_t> alive(nt, 1);

    std::unordered_map<std::array<std::int64_t, 2>, std::vector<std::size_t>, OptvolArrHash<2>>
        edge_map;
    edge_map.reserve(nt * 6);
    static const int kEdges[6][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
    for (std::size_t t = 0; t < nt; ++t)
        for (auto& ev : kEdges)
            edge_map[optvol_edge_key(rTets[t][ev[0]], rTets[t][ev[1]])].push_back(t);

    // deterministic iteration order over edge keys
    std::vector<std::array<std::int64_t, 2>> keys;
    keys.reserve(edge_map.size());
    for (const auto& kv : edge_map)
        keys.push_back(kv.first);
    std::sort(keys.begin(), keys.end());

    std::vector<Tet> new_tets;
    std::int64_t applied = 0;
    for (const auto& ek : keys) {
        const std::int64_t uu = ek[0], vv = ek[1];
        // the alive tets on this edge
        std::array<std::size_t, 3> ring_tets{};
        int nr = 0;
        bool too_many = false;
        for (std::size_t t : edge_map[ek]) {
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
        const Vec3 pu = optvol_point(rXyz, uu), pv = optvol_point(rXyz, vv);
        Vec3 pa = optvol_point(rXyz, a), pb = optvol_point(rXyz, b), pc = optvol_point(rXyz, c);
        if (optvol_svol6(pa, pb, pc, pu) < 0.0) {
            std::swap(b, c);
            std::swap(pb, pc);
        }
        // u on the positive side of (a,b,c); a valid removal needs v strictly on
        // the far side, and the two caps must tile the three old tets.
        if (optvol_svol6(pa, pb, pc, pv) >= -eps_vol6)
            continue;
        const double vcap1 = optvol_svol6(pa, pb, pc, pu);  // > 0
        const double vcap2 = optvol_svol6(pa, pc, pb, pv);  // stored orientation, should be > 0
        if (vcap1 <= eps_vol6 || vcap2 <= eps_vol6)
            continue;
        double old_sum = 0.0;
        for (int r = 0; r < 3; ++r)
            old_sum += std::fabs(optvol_svol6(
                optvol_point(rXyz, rTets[ring_tets[r]][0]), optvol_point(rXyz, rTets[ring_tets[r]][1]),
                optvol_point(rXyz, rTets[ring_tets[r]][2]),
                optvol_point(rXyz, rTets[ring_tets[r]][3])));
        const double new_sum = vcap1 + vcap2;
        if (std::fabs(new_sum - old_sum) > 1e-9 * old_sum + eps_vol6)
            continue;  // not a clean tiling -> non-convex / invalid
        const Tet cap1{a, b, c, uu}, cap2{a, c, b, vv};
        const double q_new =
            std::min(optvol_quality_tet(rXyz, cap1, eps_len), optvol_quality_tet(rXyz, cap2, eps_len));
        double q_old = 1.0;
        for (int r = 0; r < 3; ++r)
            q_old = std::min(q_old, optvol_quality_tet(rXyz, rTets[ring_tets[r]], eps_len));
        if (q_new <= q_old + min_improve)
            continue;
        for (int r = 0; r < 3; ++r)
            alive[ring_tets[r]] = 0;
        new_tets.push_back(cap1);
        new_tets.push_back(cap2);
        ++applied;
    }

    std::vector<Tet> out;
    out.reserve(nt);
    for (std::size_t t = 0; t < nt; ++t)
        if (alive[t])
            out.push_back(rTets[t]);
    for (const Tet& t : new_tets)
        out.push_back(t);
    rTets.swap(out);
    return applied;
}

// --- ODT relocation half (delegates to smooth's SmoothMethod::Odt) -----------
// Build a single-tetra-block Float64 mesh from the working buffers, run one ODT
// relocation pass through `smooth` (reusing all of its boundary/feature/frozen
// pinning and inversion guard), and read the moved points back into `rXyz`.
// Float64 throughout so a Float32 input mesh does not accumulate one rounding
// per sweep -- the single cast happens once, at final emission.
void optvol_relocate(std::vector<double>& rXyz, const std::vector<Tet>& rTets,
                     const OptimizeVolumeOptions& rOptions) {
    const std::size_t n = rXyz.size() / 3;
    Mesh tmp;
    NDArray pts = NDArray::Uninit(DType::Float64, {n, 3});
    std::memcpy(pts.Data(), rXyz.data(), rXyz.size() * sizeof(double));
    tmp.AssignPoints(std::move(pts));
    NDArray conn = NDArray::Uninit(DType::Int64, {rTets.size(), 4});
    std::int64_t* cd = conn.As<std::int64_t>();
    for (std::size_t t = 0; t < rTets.size(); ++t)
        for (int k = 0; k < 4; ++k)
            cd[t * 4 + k] = rTets[t][k];
    tmp.AddCellBlock(cell_type_name(CellType::Tetra), std::move(conn));

    SmoothOptions so;
    so.mMethod = SmoothMethod::Odt;
    so.mIterations = 1;
    so.mFixBoundary = rOptions.mPreserveBoundary;
    so.mGuardInversion = true;
    so.mFrozen = rOptions.mFrozen;
    SmoothResult sr = smooth(tmp, so);

    const NDArray& mp = sr.mMesh.Points();
    for (std::size_t i = 0; i < n * 3; ++i)
        rXyz[i] = detail::read_double(mp, i);
}

}  // namespace

OptimizeVolumeResult optimize_volume(const Mesh& rMesh, const OptimizeVolumeOptions& rOptions) {
    optvol_check_blocks(rMesh);

    OptimizeVolumeResult result;
    const std::size_t n = static_cast<std::size_t>(rMesh.NumPoints());
    const std::size_t dim = rMesh.PointDim();

    // --- working buffers ------------------------------------------------------
    std::vector<double> xyz(n * 3, 0.0);
    {
        const NDArray& points = rMesh.Points();
        parallel_for_bw(n, [&](std::size_t i) {
            for (std::size_t d = 0; d < dim && d < 3; ++d)
                xyz[i * 3 + d] = detail::read_double(points, i * dim + d);
        });
    }
    const std::vector<double> xyz0 = xyz;  // for the moved-count at the end

    std::vector<Tet> tets;
    {
        const std::vector<std::int64_t> bases = detail::block_bases(rMesh);
        tets.resize(static_cast<std::size_t>(detail::total_cells(bases)));
        std::size_t bi = 0;
        for (const auto cb : rMesh.CellRange()) {
            const NDArray& conn = cb.Conn();
            const std::size_t nc = cb.NumCells();
            const std::size_t base = static_cast<std::size_t>(bases[bi]);
            for (std::size_t c = 0; c < nc; ++c)
                for (std::size_t k = 0; k < 4; ++k)
                    tets[base + c][k] = detail::read_int(conn, c * 4 + k);
            ++bi;
        }
    }

    // Scale-relative epsilons (bbox diagonal), so validity is scale-invariant.
    double lo[3] = {0, 0, 0}, hi[3] = {0, 0, 0};
    if (n > 0) {
        for (int d = 0; d < 3; ++d)
            lo[d] = hi[d] = xyz[d];
        for (std::size_t i = 1; i < n; ++i)
            for (int d = 0; d < 3; ++d) {
                lo[d] = std::min(lo[d], xyz[i * 3 + d]);
                hi[d] = std::max(hi[d], xyz[i * 3 + d]);
            }
    }
    const double L = std::sqrt((hi[0] - lo[0]) * (hi[0] - lo[0]) + (hi[1] - lo[1]) * (hi[1] - lo[1]) +
                               (hi[2] - lo[2]) * (hi[2] - lo[2]));
    const double eps_len = (L > 0.0 ? L : 1.0) * 1e-14;
    const double eps_vol6 = (L > 0.0 ? L * L * L : 1.0) * 1e-14;

    auto min_quality = [&](const std::vector<Tet>& ts) {
        double m = 1.0;
        for (const Tet& t : ts)
            m = std::min(m, optvol_quality_tet(xyz, t, eps_len));
        return ts.empty() ? 0.0 : m;
    };
    result.mMinQualityBefore = min_quality(tets);

    // --- the optimisation loop ------------------------------------------------
    for (int sweep = 0; sweep < rOptions.mMaxIterations; ++sweep) {
        std::int64_t moved_before = 0;
        if (rOptions.mRelocate) {
            const std::vector<double> before = xyz;
            optvol_relocate(xyz, tets, rOptions);
            for (std::size_t i = 0; i < n; ++i) {
                const double dx = xyz[i * 3] - before[i * 3];
                const double dy = xyz[i * 3 + 1] - before[i * 3 + 1];
                const double dz = xyz[i * 3 + 2] - before[i * 3 + 2];
                if (std::sqrt(dx * dx + dy * dy + dz * dz) > 1e-12 * (L > 0.0 ? L : 1.0))
                    ++moved_before;
            }
        }
        std::int64_t flips_this_sweep = 0;
        if (rOptions.mFlip) {
            const std::int64_t n23 =
                optvol_pass_23(xyz, tets, rOptions.mMinImprovement, eps_len, eps_vol6);
            const std::int64_t n32 =
                optvol_pass_32(xyz, tets, rOptions.mMinImprovement, eps_len, eps_vol6);
            result.mNum23Flips += n23;
            result.mNum32Flips += n32;
            flips_this_sweep = n23 + n32;
        }
        result.mNumFlips += flips_this_sweep;
        if (moved_before == 0 && flips_this_sweep == 0)
            break;  // fixed point
    }

    result.mMinQualityAfter = min_quality(tets);
    result.mNumTets = static_cast<std::int64_t>(tets.size());
    for (std::size_t i = 0; i < n; ++i) {
        const double dx = xyz[i * 3] - xyz0[i * 3];
        const double dy = xyz[i * 3 + 1] - xyz0[i * 3 + 1];
        const double dz = xyz[i * 3 + 2] - xyz0[i * 3 + 2];
        if (std::sqrt(dx * dx + dy * dy + dz * dz) > 1e-9 * (L > 0.0 ? L : 1.0))
            ++result.mNumVerticesMoved;
    }

    // --- emit: single tetra block, points invariant --------------------------
    Mesh& out = result.mMesh;
    {
        const NDArray& points = rMesh.Points();
        NDArray np = NDArray::Uninit(points.Dtype(), {n, dim});
        detail::dispatch_dtype(points.Dtype(), [&]<class T>() {
            T* dst = np.As<T>();
            for (std::size_t i = 0; i < n; ++i)
                for (std::size_t d = 0; d < dim; ++d)
                    dst[i * dim + d] = d < 3 ? static_cast<T>(xyz[i * 3 + d])
                                             : static_cast<T>(detail::read_double(points, i * dim + d));
        });
        out.AssignPoints(std::move(np));
    }
    {
        NDArray block = NDArray::Uninit(DType::Int64, {tets.size(), 4});
        std::int64_t* dst = block.As<std::int64_t>();
        for (std::size_t t = 0; t < tets.size(); ++t) {
            const Tet ot = optvol_orient_positive(xyz, tets[t]);
            for (int k = 0; k < 4; ++k)
                dst[t * 4 + k] = ot[k];
        }
        out.AddCellBlock(cell_type_name(CellType::Tetra), std::move(block));
    }

    // point_data + field_data + Point regions carry (point set is invariant);
    // cell_data + Cell/Side regions are dropped -- flips have no cell map.
    for (const std::string& name : rMesh.PointDataNames())
        out.AddPointData(name, detail::data_owned_copy(rMesh.PointData(name)));
    for (const std::string& name : rMesh.FieldDataNames())
        out.AddFieldData(name, detail::data_owned_copy(rMesh.FieldData(name)));

    std::size_t dropped_cell_regions = 0;
    for (std::size_t i = 0; i < static_cast<std::size_t>(rMesh.NumRegions()); ++i) {
        const Region& r = rMesh.Region(i);
        if (r.mKind == RegionKind::Point)
            out.AddRegion(Region(r.mName, r.mKind, r.mDim, r.mTag, detail::data_owned_copy(r.mEntries)));
        else
            ++dropped_cell_regions;
    }
    if (!rMesh.CellDataNames().empty())
        log::warn("{}dropped {} cell_data array(s): a flip has no cell correspondence",
                  kOvPrefix, rMesh.CellDataNames().size());
    if (dropped_cell_regions > 0)
        log::warn("{}dropped {} Cell/Side region(s): a flip has no cell correspondence", kOvPrefix,
                  dropped_cell_regions);

    return result;
}

}  // namespace meshioplusplus

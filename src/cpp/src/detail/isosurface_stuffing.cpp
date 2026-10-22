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
// Isosurface stuffing of one interface (detail/isosurface_stuffing.hpp).
// Anonymous-namespace helpers are prefixed `isostuff_` (the amalgamation
// concatenates every translation unit into one).

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

// Project includes
#include "isosurface_stuffing.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

using IsostuffEdge = std::pair<std::int64_t, std::int64_t>;

double isostuff_volume(const Vec3& rA, const Vec3& rB, const Vec3& rC, const Vec3& rD) {
    return triple_product(vec3_sub(rB, rA), vec3_sub(rC, rA), vec3_sub(rD, rA)) / 6.0;
}

// The pieces of one tetrahedron on the side whose vertices have sign `Target`.
struct IsostuffCutter {
    const std::vector<Vec3>* pPoints = nullptr;
    const std::vector<IsostuffEdge>* pEdges = nullptr;  // the surviving cut edges, sorted
    std::int64_t mFirstCut = 0;                         // the id of the first cut point
    const std::vector<Vec3>* pCutPos = nullptr;

    std::int64_t Cut(std::int64_t A, std::int64_t B) const {
        const IsostuffEdge e{std::min(A, B), std::max(A, B)};
        const auto it = std::lower_bound(pEdges->begin(), pEdges->end(), e);
        return mFirstCut + static_cast<std::int64_t>(it - pEdges->begin());
    }
    const Vec3& Pos(std::int64_t V) const {
        return V < mFirstCut ? (*pPoints)[static_cast<std::size_t>(V)]
                             : (*pCutPos)[static_cast<std::size_t>(V - mFirstCut)];
    }

    void Pieces(const std::array<std::int64_t, 4>& rV, const std::array<std::int8_t, 4>& rS,
                std::int8_t Target, std::vector<std::array<std::int64_t, 4>>& rOut) const {
        std::array<std::int64_t, 4> nr{}, zr{}, fr{};
        std::size_t nn = 0, nz = 0, nf = 0;
        for (std::size_t i = 0; i < 4; ++i) {
            if (rS[i] == Target)
                nr[nn++] = rV[i];
            else if (rS[i] == 0)
                zr[nz++] = rV[i];
            else
                fr[nf++] = rV[i];
        }
        std::sort(nr.begin(), nr.begin() + static_cast<std::ptrdiff_t>(nn));
        std::sort(zr.begin(), zr.begin() + static_cast<std::ptrdiff_t>(nz));
        std::sort(fr.begin(), fr.begin() + static_cast<std::ptrdiff_t>(nf));
        if (nn == 0)
            return;
        if (nf == 0) {
            rOut.push_back(rV);
            return;
        }
        if (nn == 1) {
            const std::int64_t a = nr[0];
            if (nz == 0)
                rOut.push_back({a, Cut(a, fr[0]), Cut(a, fr[1]), Cut(a, fr[2])});
            else if (nz == 1)
                rOut.push_back({a, zr[0], Cut(a, fr[0]), Cut(a, fr[1])});
            else
                rOut.push_back({a, zr[0], zr[1], Cut(a, fr[0])});
            return;
        }
        if (nn == 2 && nf == 2) {
            // A prism: the mesh faces' quadrilaterals take the diagonal through
            // the lower near vertex `x`; the interface's its shorter one.
            const std::int64_t x = nr[0], y = nr[1], c = fr[0], d = fr[1];
            const std::int64_t xc = Cut(x, c), xd = Cut(x, d), yc = Cut(y, c), yd = Cut(y, d);
            rOut.push_back({x, y, yc, yd});
            const double l1 = vec3_norm_sq(vec3_sub(Pos(xc), Pos(yd)));
            const double l2 = vec3_norm_sq(vec3_sub(Pos(yc), Pos(xd)));
            const bool first = l1 < l2 || (l1 == l2 && std::min(xc, yd) < std::min(yc, xd));
            if (first) {
                rOut.push_back({x, xc, yc, yd});
                rOut.push_back({x, xc, yd, xd});
            } else {
                rOut.push_back({x, xc, yc, xd});
                rOut.push_back({x, yc, yd, xd});
            }
            return;
        }
        if (nn == 2) {  // (2, 1, 1): a pyramid on the zero vertex
            const std::int64_t x = nr[0], y = nr[1], z = zr[0], d = fr[0];
            const std::int64_t xd = Cut(x, d), yd = Cut(y, d);
            rOut.push_back({z, x, y, yd});
            rOut.push_back({z, x, yd, xd});
            return;
        }
        // (3, 0, 1): a prism under the far vertex.
        const std::int64_t a = nr[0], b = nr[1], c = nr[2], d = fr[0];
        const std::int64_t ad = Cut(a, d), bd = Cut(b, d), cd = Cut(c, d);
        rOut.push_back({a, b, c, cd});
        rOut.push_back({a, b, cd, bd});
        rOut.push_back({a, ad, bd, cd});
    }
};

}  // namespace

StuffReport isosurface_stuff(StuffMesh& rMesh, std::vector<std::int8_t>& rSign,
                             const StuffCut& rCut, const StuffOptions& rOptions) {
    StuffReport report;
    const std::size_t nt = rMesh.mTets.size();
    const std::size_t np = rMesh.mPoints.size();
    // The crossing edges of the material being split, sorted.
    std::vector<IsostuffEdge> edges;
    for (std::size_t t = 0; t < nt; ++t) {
        if (rMesh.mMaterial[t] != rOptions.mMaterial)
            continue;
        const auto& v = rMesh.mTets[t];
        for (std::size_t a = 0; a < 4; ++a)
            for (std::size_t b = a + 1; b < 4; ++b) {
                const auto sa = rSign[static_cast<std::size_t>(v[a])];
                const auto sb = rSign[static_cast<std::size_t>(v[b])];
                if (sa * sb < 0)
                    edges.emplace_back(std::min(v[a], v[b]), std::max(v[a], v[b]));
            }
    }
    parallel_sort(edges.begin(), edges.end());
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
    // Their cut points, and where along each edge (from its lower end) they lie.
    std::vector<Vec3> cut(edges.size());
    std::vector<double> frac(edges.size());
    parallel_for(
        edges.size(),
        [&](std::size_t e) {
            const auto [lo, hi] = edges[e];
            const bool lo_neg = rSign[static_cast<std::size_t>(lo)] < 0;
            cut[e] = lo_neg ? rCut(lo, hi) : rCut(hi, lo);
            const Vec3& a = rMesh.mPoints[static_cast<std::size_t>(lo)];
            const Vec3 ab = vec3_sub(rMesh.mPoints[static_cast<std::size_t>(hi)], a);
            const double len2 = vec3_norm_sq(ab);
            frac[e] =
                len2 > 0.0 ? std::clamp(vec3_dot(vec3_sub(cut[e], a), ab) / len2, 0.0, 1.0) : 0.5;
        },
        64);
    // Which vertices each cut point violates.
    std::vector<std::pair<std::int64_t, std::int64_t>> violates;  // (vertex, edge)
    for (std::size_t e = 0; e < edges.size(); ++e) {
        const auto [lo, hi] = edges[e];
        const bool short_edge = rMesh.mCentre[static_cast<std::size_t>(lo)] !=
                                rMesh.mCentre[static_cast<std::size_t>(hi)];
        const double alpha = short_edge ? rOptions.mAlphaShort : rOptions.mAlphaLong;
        if (frac[e] < alpha && !rMesh.mPinned[static_cast<std::size_t>(lo)])
            violates.emplace_back(lo, static_cast<std::int64_t>(e));
        if (1.0 - frac[e] < alpha && !rMesh.mPinned[static_cast<std::size_t>(hi)])
            violates.emplace_back(hi, static_cast<std::int64_t>(e));
    }
    std::sort(violates.begin(), violates.end());
    // Every tetrahedron around each vertex, for the inversion guard.
    std::vector<std::int64_t> inc_start(np + 1, 0), inc;
    for (const auto& t : rMesh.mTets)
        for (std::int64_t v : t)
            ++inc_start[static_cast<std::size_t>(v) + 1];
    for (std::size_t p = 0; p < np; ++p)
        inc_start[p + 1] += inc_start[p];
    inc.resize(static_cast<std::size_t>(inc_start[np]));
    {
        std::vector<std::int64_t> fill(inc_start.begin(), inc_start.end() - 1);
        for (std::size_t t = 0; t < nt; ++t)
            for (std::int64_t v : rMesh.mTets[t])
                inc[static_cast<std::size_t>(fill[static_cast<std::size_t>(v)]++)] =
                    static_cast<std::int64_t>(t);
    }
    // Warp, vertex by vertex in ascending order, onto the nearest live cut point.
    for (std::size_t k = 0; k < violates.size();) {
        const std::int64_t v = violates[k].first;
        std::size_t end = k;
        while (end < violates.size() && violates[end].first == v)
            ++end;
        const auto vi = static_cast<std::size_t>(v);
        std::int64_t best = -1;
        double best_d = 0.0;
        for (std::size_t j = k; j < end; ++j) {
            const auto e = static_cast<std::size_t>(violates[j].second);
            const auto [lo, hi] = edges[e];
            if (rSign[static_cast<std::size_t>(lo)] * rSign[static_cast<std::size_t>(hi)] >= 0)
                continue;  // a warp at its other end removed it
            const double d = vec3_norm_sq(vec3_sub(cut[e], rMesh.mPoints[vi]));
            if (best < 0 || d < best_d) {
                best = static_cast<std::int64_t>(e);
                best_d = d;
            }
        }
        k = end;
        if (best < 0 || rSign[vi] == 0)
            continue;
        const Vec3 target = cut[static_cast<std::size_t>(best)];
        bool ok = true;
        for (std::int64_t q = inc_start[vi]; q < inc_start[vi + 1] && ok; ++q) {
            const auto& t = rMesh.mTets[static_cast<std::size_t>(inc[static_cast<std::size_t>(q)])];
            Vec3 p[4];
            for (std::size_t c = 0; c < 4; ++c)
                p[c] = t[c] == v ? target : rMesh.mPoints[static_cast<std::size_t>(t[c])];
            ok = isostuff_volume(p[0], p[1], p[2], p[3]) > rOptions.mMinVolume;
        }
        if (!ok) {
            ++report.mNumWarpRefused;
            continue;
        }
        rMesh.mPoints[vi] = target;
        rSign[vi] = 0;
        rMesh.mPinned[vi] = 1;
        ++report.mNumWarped;
    }
    // The surviving cut points become points.
    std::vector<IsostuffEdge> live;
    std::vector<Vec3> live_pos;
    for (std::size_t e = 0; e < edges.size(); ++e) {
        const auto [lo, hi] = edges[e];
        if (rSign[static_cast<std::size_t>(lo)] * rSign[static_cast<std::size_t>(hi)] < 0) {
            live.push_back(edges[e]);
            live_pos.push_back(cut[e]);
        }
    }
    report.mNumCutPoints = static_cast<std::int64_t>(live.size());
    IsostuffCutter cutter;
    cutter.pPoints = &rMesh.mPoints;
    cutter.pEdges = &live;
    cutter.mFirstCut = static_cast<std::int64_t>(np);
    cutter.pCutPos = &live_pos;
    // The pieces, tetrahedron by tetrahedron.
    std::vector<std::vector<std::array<std::int64_t, 4>>> near(nt), far(nt);
    parallel_for(
        nt,
        [&](std::size_t t) {
            if (rMesh.mMaterial[t] != rOptions.mMaterial)
                return;
            const auto& v = rMesh.mTets[t];
            std::array<std::int8_t, 4> s{};
            std::size_t nonzero = 0;
            for (std::size_t c = 0; c < 4; ++c) {
                s[c] = rSign[static_cast<std::size_t>(v[c])];
                nonzero += s[c] != 0;
            }
            if (nonzero == 0) {
                near[t].push_back(v);  // on the interface: the near side's
                return;
            }
            cutter.Pieces(v, s, -1, near[t]);
            cutter.Pieces(v, s, +1, far[t]);
        },
        256);
    StuffMesh out;
    out.mPoints = rMesh.mPoints;
    out.mPoints.insert(out.mPoints.end(), live_pos.begin(), live_pos.end());
    out.mCentre = rMesh.mCentre;
    out.mCentre.resize(out.mPoints.size(), 0);
    out.mPinned = rMesh.mPinned;
    out.mPinned.resize(out.mPoints.size(), 1);
    for (std::size_t t = 0; t < nt; ++t) {
        auto push = [&](std::array<std::int64_t, 4> rT, std::int32_t Material, bool Whole) {
            Vec3 p[4];
            for (std::size_t c = 0; c < 4; ++c)
                p[c] = out.mPoints[static_cast<std::size_t>(rT[c])];
            const double vol = isostuff_volume(p[0], p[1], p[2], p[3]);
            if (vol < 0.0)
                std::swap(rT[2], rT[3]);
            if (std::abs(vol) <= rOptions.mMinVolume)
                ++report.mNumDegenerate;
            out.mTets.push_back(rT);
            out.mMaterial.push_back(Material);
            out.mParent.push_back(rMesh.mParent[t]);
            out.mCut.push_back(Whole ? rMesh.mCut[t] : std::uint8_t{1});
        };
        if (rMesh.mMaterial[t] != rOptions.mMaterial) {
            push(rMesh.mTets[t], rMesh.mMaterial[t], true);
            continue;
        }
        const bool whole_near =
            near[t].size() == 1 && far[t].empty() && near[t][0] == rMesh.mTets[t];
        const bool whole_far = far[t].size() == 1 && near[t].empty() && far[t][0] == rMesh.mTets[t];
        for (const auto& piece : near[t])
            push(piece, rOptions.mMaterial + 1, whole_near);
        for (const auto& piece : far[t])
            push(piece, rOptions.mMaterial, whole_far);
    }
    rSign.resize(out.mPoints.size(), 0);
    rMesh = std::move(out);
    return report;
}

}  // namespace detail
}  // namespace meshioplusplus

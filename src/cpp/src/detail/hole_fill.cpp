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
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <utility>

// Project includes
#include "hole_fill.hpp"
#include "sparse_solve.hpp"

namespace meshioplusplus::detail {

namespace {

using HfTri = std::array<std::size_t, 3>;

const double kHfPi = std::acos(-1.0);

double hf_area(const Vec3& rA, const Vec3& rB, const Vec3& rC) {
    return 0.5 * vec3_norm(vec3_cross(vec3_sub(rB, rA), vec3_sub(rC, rA)));
}

// The angle at A of the triangle A, B, C.
double hf_angle(const Vec3& rA, const Vec3& rB, const Vec3& rC) {
    const Vec3 u = vec3_sub(rB, rA);
    const Vec3 v = vec3_sub(rC, rA);
    const double nu = vec3_norm(u), nv = vec3_norm(v);
    if (!(nu > 0.0 && nv > 0.0))
        return 0.0;
    return std::acos(std::clamp(vec3_dot(u, v) / (nu * nv), -1.0, 1.0));
}

// The minimum-area triangulation of the polygon 0 .. n-1 (ties: the lowest
// split vertex).
std::vector<HfTri> hf_min_area(const std::vector<Vec3>& rP, std::size_t N) {
    std::vector<double> w(N * N, 0.0);
    std::vector<std::size_t> split(N * N, 0);
    for (std::size_t len = 2; len < N; ++len)
        for (std::size_t i = 0; i + len < N; ++i) {
            const std::size_t j = i + len;
            double best = std::numeric_limits<double>::infinity();
            std::size_t arg = i + 1;
            for (std::size_t m = i + 1; m < j; ++m) {
                const double c = w[i * N + m] + w[m * N + j] + hf_area(rP[i], rP[m], rP[j]);
                if (c < best) {
                    best = c;
                    arg = m;
                }
            }
            w[i * N + j] = best;
            split[i * N + j] = arg;
        }
    std::vector<HfTri> tris;
    std::vector<std::pair<std::size_t, std::size_t>> stack{{0, N - 1}};
    while (!stack.empty()) {
        const auto [i, j] = stack.back();
        stack.pop_back();
        if (j < i + 2)
            continue;
        const std::size_t m = split[i * N + j];
        tris.push_back({i, m, j});
        stack.push_back({m, j});
        stack.push_back({i, m});
    }
    return tris;
}

// Flip the interior edges to the Delaunay criterion until none flips (at
// most 64 passes); an edge whose flip would duplicate an edge stays.
void hf_relax(const std::vector<Vec3>& rP, std::vector<HfTri>& rTris) {
    for (int pass = 0; pass < 64; ++pass) {
        bool flipped = false;
        std::map<std::pair<std::size_t, std::size_t>, std::vector<std::size_t>> uses;
        for (std::size_t t = 0; t < rTris.size(); ++t)
            for (std::size_t e = 0; e < 3; ++e) {
                const std::size_t a = rTris[t][e], b = rTris[t][(e + 1) % 3];
                uses[{std::min(a, b), std::max(a, b)}].push_back(t);
            }
        std::vector<char> touched(rTris.size(), 0);
        for (const auto& [edge, r_ts] : uses) {
            if (r_ts.size() != 2 || touched[r_ts[0]] || touched[r_ts[1]])
                continue;
            const HfTri& t1 = rTris[r_ts[0]];
            const HfTri& t2 = rTris[r_ts[1]];
            // t1 = (a, b, c) traverses a -> b; t2 = (b, a, d).
            std::size_t k = 0;
            while (!(t1[k] != edge.first && t1[k] != edge.second))
                ++k;
            const std::size_t c = t1[k], a = t1[(k + 1) % 3], b = t1[(k + 2) % 3];
            std::size_t d = t2[0];
            for (const std::size_t v : t2)
                if (v != a && v != b)
                    d = v;
            if (uses.count({std::min(c, d), std::max(c, d)}))
                continue;
            if (hf_angle(rP[c], rP[a], rP[b]) + hf_angle(rP[d], rP[a], rP[b]) <= kHfPi + 1e-12)
                continue;
            rTris[r_ts[0]] = {c, a, d};
            rTris[r_ts[1]] = {d, b, c};
            touched[r_ts[0]] = touched[r_ts[1]] = 1;
            flipped = true;
        }
        if (!flipped)
            return;
    }
}

// Split triangles at their centroids to the spacing (Liepa's density rule),
// relaxing after each pass; at most 32 passes.
void hf_refine(std::vector<Vec3>& rP, std::vector<double>& rSigma, std::vector<HfTri>& rTris) {
    const double alpha = std::sqrt(2.0);
    for (int pass = 0; pass < 32; ++pass) {
        bool split = false;
        const std::size_t nt = rTris.size();
        for (std::size_t t = 0; t < nt; ++t) {
            const HfTri tri = rTris[t];
            const Vec3 c =
                vec3_scale(vec3_add(vec3_add(rP[tri[0]], rP[tri[1]]), rP[tri[2]]), 1.0 / 3.0);
            const double sc = (rSigma[tri[0]] + rSigma[tri[1]] + rSigma[tri[2]]) / 3.0;
            bool ok = true;
            for (const std::size_t v : tri) {
                const double d = alpha * vec3_norm(vec3_sub(c, rP[v]));
                ok = ok && d > sc && d > rSigma[v];
            }
            if (!ok)
                continue;
            const std::size_t p = rP.size();
            rP.push_back(c);
            rSigma.push_back(sc);
            rTris[t] = {tri[0], tri[1], p};
            rTris.push_back({tri[1], tri[2], p});
            rTris.push_back({tri[2], tri[0], p});
            split = true;
        }
        if (!split)
            return;
        hf_relax(rP, rTris);
    }
}

}  // namespace

HoleFillPatch hole_fill(const HoleFillInput& rIn, std::int64_t Base) {
    const std::vector<Vec3>& r_all = *rIn.pPoints;
    const std::size_t n = rIn.mLoop.size();
    std::vector<Vec3> p(n);
    std::vector<double> sigma(n);
    for (std::size_t i = 0; i < n; ++i) {
        p[i] = r_all[static_cast<std::size_t>(rIn.mLoop[i])];
        sigma[i] = rIn.mRimSpacing[i];
    }
    std::vector<HfTri> tris = hf_min_area(p, n);
    if (rIn.mRefine) {
        hf_relax(p, tris);
        hf_refine(p, sigma, tris);
    }
    const std::size_t np = p.size();
    // Global ids: the loop's own, then Base + k.
    const auto gid = [&](std::size_t L) {
        return L < n ? rIn.mLoop[L] : Base + static_cast<std::int64_t>(L - n);
    };
    if (rIn.mFairing != HoleFairing::None && np > n) {
        // Neighbours, as local ids; an outside point is np + its index.
        std::vector<std::vector<std::size_t>> nbr(np);
        for (const HfTri& r_t : tris)
            for (std::size_t e = 0; e < 3; ++e) {
                nbr[r_t[e]].push_back(r_t[(e + 1) % 3]);
                nbr[r_t[e]].push_back(r_t[(e + 2) % 3]);
            }
        std::map<std::int64_t, std::size_t> local;
        for (std::size_t i = 0; i < n; ++i)
            local[rIn.mLoop[i]] = i;
        std::vector<Vec3> outside;
        for (std::size_t i = 0; i < n; ++i)
            for (const std::int64_t g : rIn.mRimRing[i]) {
                const auto it = local.find(g);
                if (it != local.end()) {
                    nbr[i].push_back(it->second);
                    continue;
                }
                local[g] = np + outside.size();
                nbr[i].push_back(np + outside.size());
                outside.push_back(r_all[static_cast<std::size_t>(g)]);
            }
        for (auto& r_n : nbr) {
            std::sort(r_n.begin(), r_n.end());
            r_n.erase(std::unique(r_n.begin(), r_n.end()), r_n.end());
        }
        const auto pos = [&](std::size_t L) { return L < np ? p[L] : outside[L - np]; };
        const std::size_t nu = np - n;  // unknowns: local n .. np-1
        // Rows: (column, coefficient) per unknown, and the known part.
        std::vector<SparseTriplet> trip;
        std::vector<std::array<double, 3>> rhs(nu, {0.0, 0.0, 0.0});
        const auto add_row = [&](std::size_t V, bool Squared) {
            // U_V = mean of neighbours - x_V, linear in the unknowns.
            const double inv = 1.0 / static_cast<double>(nbr[V].size());
            std::vector<std::pair<std::size_t, double>> coef;
            Vec3 known{0.0, 0.0, 0.0};
            const auto term = [&](std::size_t L, double C) {
                if (L >= n && L < np)
                    coef.push_back({L - n, C});
                else
                    known = vec3_add(known, vec3_scale(pos(L), C));
            };
            for (const std::size_t q : nbr[V])
                term(q, inv);
            term(V, -1.0);
            if (!Squared) {
                // Membrane: deg(V) U_V = 0 for an unknown V, the graph
                // Laplacian's row (symmetric, positive definite with the rim
                // fixed).
                const double deg = static_cast<double>(nbr[V].size());
                for (const auto& [c, w] : coef)
                    trip.push_back({V - n, c, -w * deg});
                for (std::size_t k = 0; k < 3; ++k)
                    rhs[V - n][k] += known[k] * deg;
                return;
            }
            for (const auto& [ci, wi] : coef) {
                for (const auto& [cj, wj] : coef)
                    trip.push_back({ci, cj, wi * wj});
                for (std::size_t k = 0; k < 3; ++k)
                    rhs[ci][k] -= wi * known[k];
            }
        };
        if (rIn.mFairing == HoleFairing::Membrane) {
            for (std::size_t v = n; v < np; ++v)
                add_row(v, false);
        } else {
            for (std::size_t v = 0; v < np; ++v)
                add_row(v, true);
        }
        const SparseCsr a = csr_from_triplets(nu, nu, std::move(trip));
        for (std::size_t k = 0; k < 3; ++k) {
            std::vector<double> b(nu), x(nu);
            for (std::size_t u = 0; u < nu; ++u) {
                b[u] = rhs[u][k];
                x[u] = p[n + u][k];
            }
            csr_pcg(a, b.data(), x.data(), 1e-13, 0);
            for (std::size_t u = 0; u < nu; ++u)
                p[n + u][k] = x[u];
        }
    }
    HoleFillPatch out;
    out.mPoints.assign(p.begin() + static_cast<std::ptrdiff_t>(n), p.end());
    for (const HfTri& r_t : tris)
        out.mTriangles.push_back({gid(r_t[0]), gid(r_t[1]), gid(r_t[2])});
    return out;
}

}  // namespace meshioplusplus::detail

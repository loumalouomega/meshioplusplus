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
// The triangle BVH of detail/tri_bvh.hpp: fast winding numbers and ray casts.
// Anonymous-namespace helpers are prefixed `tribvh_` (the amalgamation
// concatenates every translation unit into one).

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
#include "tri_bvh.hpp"
#include "meshioplusplus/detail/space_filling.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

constexpr double kTribvhFourPi = 12.566370614359172953850573533118;

std::int32_t tribvh_build_nodes(TriBvh& rTree, std::int32_t Begin, std::int32_t End,
                                std::int32_t LeafSize) {
    const auto at = static_cast<std::int32_t>(rTree.mNodes.size());
    rTree.mNodes.emplace_back();
    rTree.mNodes.back().mBegin = Begin;
    rTree.mNodes.back().mCount = End - Begin;
    if (End - Begin <= LeafSize)
        return at;
    const std::int32_t mid = Begin + (End - Begin) / 2;
    tribvh_build_nodes(rTree, Begin, mid, LeafSize);
    const std::int32_t right = tribvh_build_nodes(rTree, mid, End, LeafSize);
    rTree.mNodes[static_cast<std::size_t>(at)].mRight = right;
    return at;
}

// The far-field solid angle of node `I` from a point at `r = p - q`.
double tribvh_far(const TriBvh& rTree, std::size_t I, const Vec3& rR) {
    const double d2 = vec3_norm_sq(rR);
    const double d = std::sqrt(d2);
    const double inv3 = 1.0 / (d2 * d);
    const double inv5 = inv3 / d2;
    const Vec3& n = rTree.mNormalSum[I];
    double omega = (n[0] * rR[0] + n[1] * rR[1] + n[2] * rR[2]) * inv3;
    const auto& m = rTree.mMoment[I];
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j) {
            const double jij = (i == j ? inv3 : 0.0) - 3.0 * rR[i] * rR[j] * inv5;
            omega += m[i * 3 + j] * jij;
        }
    // ½ H:C with H_ijk = -3 (δ_ij r_k + δ_ik r_j + δ_jk r_i) / |r|^5
    // + 15 r_i r_j r_k / |r|^7, C symmetric in (j, k).
    const auto& c = rTree.mSecond[I];
    const double inv7 = inv5 / d2;
    double diag = 0.0, cubic = 0.0;
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t k = 0; k < 3; ++k) {
            diag += 2.0 * c[9 * i + 3 * i + k] * rR[k] + c[9 * k + 3 * i + i] * rR[k];
            for (std::size_t j = 0; j < 3; ++j)
                cubic += c[9 * i + 3 * j + k] * rR[i] * rR[j] * rR[k];
        }
    omega += 0.5 * (-3.0 * inv5 * diag + 15.0 * inv7 * cubic);
    return omega;
}

// Whether the ray reaches the box no later than `TMax`, and from where.
bool tribvh_slab(const Vec3& rLo, const Vec3& rHi, const Vec3& rOrigin, const Vec3& rDir,
                 double TMin, double TMax) {
    double t0 = TMin, t1 = TMax;
    for (std::size_t k = 0; k < 3; ++k) {
        if (rDir[k] == 0.0) {
            if (rOrigin[k] < rLo[k] || rOrigin[k] > rHi[k])
                return false;
            continue;
        }
        const double inv = 1.0 / rDir[k];
        double a = (rLo[k] - rOrigin[k]) * inv, b = (rHi[k] - rOrigin[k]) * inv;
        if (a > b)
            std::swap(a, b);
        t0 = std::max(t0, a);
        t1 = std::min(t1, b);
        if (t0 > t1)
            return false;
    }
    return true;
}

}  // namespace

double tri_solid_angle(const Vec3& rA, const Vec3& rB, const Vec3& rC, const Vec3& rQ) {
    const Vec3 a = vec3_sub(rA, rQ), b = vec3_sub(rB, rQ), c = vec3_sub(rC, rQ);
    const double la = vec3_norm(a), lb = vec3_norm(b), lc = vec3_norm(c);
    const double num = vec3_dot(a, vec3_cross(b, c));
    const double den =
        la * lb * lc + vec3_dot(a, b) * lc + vec3_dot(a, c) * lb + vec3_dot(b, c) * la;
    return 2.0 * std::atan2(num, den);
}

TriBvh build_tri_bvh(const TriangleSoup& rSoup, std::int32_t LeafSize) {
    TriBvh tree;
    tree.mpSoup = &rSoup;
    const std::size_t n = rSoup.NumTriangles();
    if (n == 0)
        return tree;
    // Morton order of the centroids over their bounding box.
    std::vector<Vec3> centroid(n);
    Vec3 lo{std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(),
            std::numeric_limits<double>::infinity()};
    Vec3 hi{-lo[0], -lo[1], -lo[2]};
    for (std::size_t t = 0; t < n; ++t) {
        const Vec3* c = &rSoup.mCorners[3 * t];
        for (std::size_t k = 0; k < 3; ++k) {
            centroid[t][k] = (c[0][k] + c[1][k] + c[2][k]) / 3.0;
            lo[k] = std::min(lo[k], centroid[t][k]);
            hi[k] = std::max(hi[k], centroid[t][k]);
        }
    }
    const double scale = static_cast<double>((std::uint32_t{1} << sfc_bits) - 1);
    std::vector<std::uint64_t> key(n);
    for (std::size_t t = 0; t < n; ++t) {
        std::uint32_t q[3];
        for (std::size_t k = 0; k < 3; ++k) {
            const double ext = hi[k] - lo[k];
            const double f = ext > 0.0 ? (centroid[t][k] - lo[k]) / ext : 0.0;
            q[k] = static_cast<std::uint32_t>(std::clamp(f, 0.0, 1.0) * scale);
        }
        key[t] = sfc_morton_key(q);
    }
    const std::vector<std::int64_t> order = sfc_stable_argsort(key);
    tree.mOrder.resize(n);
    for (std::size_t i = 0; i < n; ++i)
        tree.mOrder[i] = static_cast<std::int32_t>(order[i]);
    tribvh_build_nodes(tree, 0, static_cast<std::int32_t>(n), std::max<std::int32_t>(LeafSize, 1));

    // Node data, children before parents (they come later in preorder).
    const std::size_t nn = tree.mNodes.size();
    tree.mCentre.assign(nn, {0.0, 0.0, 0.0});
    tree.mNormalSum.assign(nn, {0.0, 0.0, 0.0});
    tree.mMoment.assign(nn, {});
    tree.mSecond.assign(nn, {});
    tree.mRadius.assign(nn, 0.0);
    std::vector<double> area(nn, 0.0);
    for (std::size_t i = nn; i-- > 0;) {
        TriBvhNode& node = tree.mNodes[i];
        if (node.mRight < 0) {
            Vec3 sum{0.0, 0.0, 0.0}, plain{0.0, 0.0, 0.0}, nsum{0.0, 0.0, 0.0};
            Vec3 blo = rSoup.mCorners[3 * static_cast<std::size_t>(tree.mOrder[node.mBegin])];
            Vec3 bhi = blo;
            double a = 0.0;
            for (std::int32_t k = 0; k < node.mCount; ++k) {
                const auto t = static_cast<std::size_t>(tree.mOrder[node.mBegin + k]);
                const Vec3* c = &rSoup.mCorners[3 * t];
                const Vec3 nt =
                    vec3_scale(vec3_cross(vec3_sub(c[1], c[0]), vec3_sub(c[2], c[0])), 0.5);
                const double at = vec3_norm(nt);
                a += at;
                nsum = vec3_add(nsum, nt);
                sum = vec3_add(sum, vec3_scale(centroid[t], at));
                plain = vec3_add(plain, centroid[t]);
                for (std::size_t v = 0; v < 3; ++v)
                    for (std::size_t d = 0; d < 3; ++d) {
                        blo[d] = std::min(blo[d], c[v][d]);
                        bhi[d] = std::max(bhi[d], c[v][d]);
                    }
            }
            const Vec3 p = a > 0.0 ? vec3_scale(sum, 1.0 / a)
                                   : vec3_scale(plain, 1.0 / static_cast<double>(node.mCount));
            std::array<double, 9> m{};
            std::array<double, 27> sec{};
            double r = 0.0;
            for (std::int32_t k = 0; k < node.mCount; ++k) {
                const auto t = static_cast<std::size_t>(tree.mOrder[node.mBegin + k]);
                const Vec3* c = &rSoup.mCorners[3 * t];
                const Vec3 nt =
                    vec3_scale(vec3_cross(vec3_sub(c[1], c[0]), vec3_sub(c[2], c[0])), 0.5);
                const double at = vec3_norm(nt);
                const Vec3 dc = vec3_sub(centroid[t], p);
                // ∫ (x - p)(x - p)^T dA / a = dc dc^T + (1/12) Σ_v (v - g)(v - g)^T.
                std::array<double, 9> s{};
                for (std::size_t y = 0; y < 3; ++y)
                    for (std::size_t z = 0; z < 3; ++z)
                        s[y * 3 + z] = dc[y] * dc[z];
                for (std::size_t v = 0; v < 3; ++v) {
                    const Vec3 e = vec3_sub(c[v], centroid[t]);
                    for (std::size_t y = 0; y < 3; ++y)
                        for (std::size_t z = 0; z < 3; ++z)
                            s[y * 3 + z] += e[y] * e[z] / 12.0;
                }
                const Vec3 un = at > 0.0 ? vec3_scale(nt, 1.0 / at) : Vec3{0.0, 0.0, 0.0};
                for (std::size_t x = 0; x < 3; ++x) {
                    for (std::size_t y = 0; y < 3; ++y) {
                        m[x * 3 + y] += nt[x] * dc[y];
                        for (std::size_t z = 0; z < 3; ++z)
                            sec[9 * x + 3 * y + z] += un[x] * at * s[y * 3 + z];
                    }
                }
                for (std::size_t v = 0; v < 3; ++v)
                    r = std::max(r, vec3_norm(vec3_sub(c[v], p)));
            }
            node.mLo = blo;
            node.mHi = bhi;
            area[i] = a;
            tree.mCentre[i] = p;
            tree.mNormalSum[i] = nsum;
            tree.mMoment[i] = m;
            tree.mSecond[i] = sec;
            tree.mRadius[i] = r;
            continue;
        }
        const std::size_t l = i + 1, rr = static_cast<std::size_t>(node.mRight);
        const double a = area[l] + area[rr];
        const Vec3 p = a > 0.0 ? vec3_scale(vec3_add(vec3_scale(tree.mCentre[l], area[l]),
                                                     vec3_scale(tree.mCentre[rr], area[rr])),
                                            1.0 / a)
                               : vec3_scale(vec3_add(tree.mCentre[l], tree.mCentre[rr]), 0.5);
        std::array<double, 9> m{};
        std::array<double, 27> sec{};
        for (std::size_t child : {l, rr}) {
            // With y = x - p_child and s = p_child - p: ∫ n (y + s)(y + s)^T
            // = C + M s^T + (M s^T)^T + N s s^T, M and N the child's.
            const Vec3 shift = vec3_sub(tree.mCentre[child], p);
            const Vec3& nc = tree.mNormalSum[child];
            const auto& mc = tree.mMoment[child];
            const auto& cc = tree.mSecond[child];
            for (std::size_t x = 0; x < 3; ++x)
                for (std::size_t y = 0; y < 3; ++y) {
                    m[x * 3 + y] += mc[x * 3 + y] + nc[x] * shift[y];
                    for (std::size_t z = 0; z < 3; ++z)
                        sec[9 * x + 3 * y + z] += cc[9 * x + 3 * y + z] + mc[x * 3 + y] * shift[z] +
                                                  mc[x * 3 + z] * shift[y] +
                                                  nc[x] * shift[y] * shift[z];
                }
        }
        for (std::size_t d = 0; d < 3; ++d) {
            node.mLo[d] = std::min(tree.mNodes[l].mLo[d], tree.mNodes[rr].mLo[d]);
            node.mHi[d] = std::max(tree.mNodes[l].mHi[d], tree.mNodes[rr].mHi[d]);
        }
        area[i] = a;
        tree.mCentre[i] = p;
        tree.mNormalSum[i] = vec3_add(tree.mNormalSum[l], tree.mNormalSum[rr]);
        tree.mMoment[i] = m;
        tree.mSecond[i] = sec;
        tree.mRadius[i] = std::max(vec3_norm(vec3_sub(tree.mCentre[l], p)) + tree.mRadius[l],
                                   vec3_norm(vec3_sub(tree.mCentre[rr], p)) + tree.mRadius[rr]);
    }
    return tree;
}

double fast_winding(const TriBvh& rTree, const Vec3& rQ, double Beta) {
    if (rTree.mNodes.empty())
        return 0.0;
    const TriangleSoup& soup = *rTree.mpSoup;
    double omega = 0.0;
    std::int32_t stack[128];
    int top = 0;
    stack[top++] = 0;
    while (top > 0) {
        const auto i = static_cast<std::size_t>(stack[--top]);
        const TriBvhNode& node = rTree.mNodes[i];
        if (Beta > 0.0) {
            const Vec3 r = vec3_sub(rTree.mCentre[i], rQ);
            if (vec3_norm(r) > Beta * rTree.mRadius[i]) {
                omega += tribvh_far(rTree, i, r);
                continue;
            }
        }
        if (node.mRight < 0) {
            for (std::int32_t k = 0; k < node.mCount; ++k) {
                const Vec3* c =
                    &soup.mCorners[3 * static_cast<std::size_t>(rTree.mOrder[node.mBegin + k])];
                omega += tri_solid_angle(c[0], c[1], c[2], rQ);
            }
            continue;
        }
        stack[top++] = node.mRight;
        stack[top++] = static_cast<std::int32_t>(i + 1);
    }
    return omega / kTribvhFourPi;
}

std::vector<double> fast_windings(const TriBvh& rTree, const std::vector<Vec3>& rPoints,
                                  double Beta) {
    std::vector<double> out(rPoints.size(), 0.0);
    parallel_for(
        rPoints.size(), [&](std::size_t i) { out[i] = fast_winding(rTree, rPoints[i], Beta); },
        256);
    return out;
}

TriRayHit tri_bvh_ray(const TriBvh& rTree, const Vec3& rOrigin, const Vec3& rDir, double TMin,
                      double TMax, std::int64_t Skip) {
    TriRayHit best;
    if (rTree.mNodes.empty())
        return best;
    const TriangleSoup& soup = *rTree.mpSoup;
    double best_t = TMax;
    std::int32_t stack[128];
    int top = 0;
    stack[top++] = 0;
    while (top > 0) {
        const auto i = static_cast<std::size_t>(stack[--top]);
        const TriBvhNode& node = rTree.mNodes[i];
        if (!tribvh_slab(node.mLo, node.mHi, rOrigin, rDir, TMin, best_t))
            continue;
        if (node.mRight >= 0) {
            stack[top++] = node.mRight;
            stack[top++] = static_cast<std::int32_t>(i + 1);
            continue;
        }
        for (std::int32_t k = 0; k < node.mCount; ++k) {
            const std::int64_t t = rTree.mOrder[node.mBegin + k];
            if (t == Skip)
                continue;
            const Vec3* c = &soup.mCorners[3 * static_cast<std::size_t>(t)];
            const Vec3 e1 = vec3_sub(c[1], c[0]), e2 = vec3_sub(c[2], c[0]);
            const Vec3 p = vec3_cross(rDir, e2);
            const double det = vec3_dot(e1, p);
            if (std::abs(det) < 1e-300)
                continue;
            const double inv = 1.0 / det;
            const Vec3 s = vec3_sub(rOrigin, c[0]);
            // A ray through a shared edge or vertex must not slip between
            // the triangles there: the barycentric tests allow a sliver.
            constexpr double kEps = 1e-10;
            const double u = vec3_dot(s, p) * inv;
            if (u < -kEps || u > 1.0 + kEps)
                continue;
            const Vec3 q = vec3_cross(s, e1);
            const double v = vec3_dot(rDir, q) * inv;
            if (v < -kEps || u + v > 1.0 + kEps)
                continue;
            const double th = vec3_dot(e2, q) * inv;
            if (th <= TMin || th > best_t)
                continue;
            if (th < best_t || best.mTriangle < 0 || t < best.mTriangle) {
                best_t = th;
                best.mT = th;
                best.mTriangle = t;
                best.mNormal = vec3_cross(e1, e2);
            }
        }
    }
    return best;
}

std::vector<TriSegmentHit> tri_bvh_segment_hits(const TriBvh& rTree, const Vec3& rA,
                                                const Vec3& rB) {
    std::vector<std::pair<double, double>> raw;  // (t, side)
    if (rTree.mNodes.empty())
        return {};
    const TriangleSoup& soup = *rTree.mpSoup;
    const Vec3 dir = vec3_sub(rB, rA);
    std::int32_t stack[128];
    int top = 0;
    stack[top++] = 0;
    constexpr double kEps = 1e-10;
    while (top > 0) {
        const auto i = static_cast<std::size_t>(stack[--top]);
        const TriBvhNode& node = rTree.mNodes[i];
        if (!tribvh_slab(node.mLo, node.mHi, rA, dir, 0.0, 1.0))
            continue;
        if (node.mRight >= 0) {
            stack[top++] = node.mRight;
            stack[top++] = static_cast<std::int32_t>(i + 1);
            continue;
        }
        for (std::int32_t k = 0; k < node.mCount; ++k) {
            const Vec3* c =
                &soup.mCorners[3 * static_cast<std::size_t>(rTree.mOrder[node.mBegin + k])];
            const Vec3 e1 = vec3_sub(c[1], c[0]), e2 = vec3_sub(c[2], c[0]);
            const Vec3 p = vec3_cross(dir, e2);
            const double det = vec3_dot(e1, p);
            if (std::abs(det) < 1e-300)
                continue;
            const double inv = 1.0 / det;
            const Vec3 s = vec3_sub(rA, c[0]);
            const double u = vec3_dot(s, p) * inv;
            if (u < -kEps || u > 1.0 + kEps)
                continue;
            const Vec3 q = vec3_cross(s, e1);
            const double v = vec3_dot(dir, q) * inv;
            if (v < -kEps || u + v > 1.0 + kEps)
                continue;
            const double th = vec3_dot(e2, q) * inv;
            if (th < 0.0 || th > 1.0)
                continue;
            raw.emplace_back(th, vec3_dot(vec3_cross(e1, e2), dir) > 0.0 ? 1.0 : -1.0);
        }
    }
    std::sort(raw.begin(), raw.end());
    std::vector<TriSegmentHit> out;
    for (std::size_t i = 0; i < raw.size();) {
        std::size_t j = i;
        double side = 0.0;
        while (j < raw.size() && raw[j].first - raw[i].first <= 1e-9) {
            side += raw[j].second;
            ++j;
        }
        if (side != 0.0)
            out.push_back({raw[i].first, side > 0.0 ? 1.0 : -1.0});
        i = j;
    }
    return out;
}

}  // namespace detail
}  // namespace meshioplusplus

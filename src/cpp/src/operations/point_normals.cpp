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
#include <limits>
#include <queue>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/point_normals.hpp"
#include "../detail/bounding_box_utils.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/sym3_eigen.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/operations/neighbors.hpp"
#include "meshioplusplus/operations/normals.hpp"
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/detail/geometry.hpp"

namespace meshioplusplus {

namespace {

using PtnVec = std::array<double, 3>;

/// The sign that makes the largest-magnitude component positive (the lowest
/// axis on a tie).
double ptn_canonical_sign(const PtnVec& rN) {
    std::size_t k = 0;
    for (std::size_t j = 1; j < 3; ++j)
        if (std::abs(rN[j]) > std::abs(rN[k]))
            k = j;
    return rN[k] < 0.0 ? -1.0 : 1.0;
}

/// The plane fit of one neighbourhood: false when it spans no plane.
bool ptn_fit(const std::vector<PtnVec>& rP, std::size_t Self, const std::int64_t* pNbrs,
             std::size_t K, PtnVec& rNormal, double& rVariation) {
    PtnVec c{0.0, 0.0, 0.0};
    const auto add = [&](const PtnVec& rQ) {
        c[0] += rQ[0];
        c[1] += rQ[1];
        c[2] += rQ[2];
    };
    add(rP[Self]);
    for (std::size_t j = 0; j < K; ++j)
        add(rP[static_cast<std::size_t>(pNbrs[j])]);
    const double m = static_cast<double>(K + 1);
    c = {c[0] / m, c[1] / m, c[2] / m};
    double t[6] = {0, 0, 0, 0, 0, 0};  // xx yy zz xy yz zx
    const auto acc = [&](const PtnVec& rQ) {
        const double x = rQ[0] - c[0], y = rQ[1] - c[1], z = rQ[2] - c[2];
        t[0] += x * x;
        t[1] += y * y;
        t[2] += z * z;
        t[3] += x * y;
        t[4] += y * z;
        t[5] += z * x;
    };
    acc(rP[Self]);
    for (std::size_t j = 0; j < K; ++j)
        acc(rP[static_cast<std::size_t>(pNbrs[j])]);
    double l[3];
    detail::sym3_principal(t, l);
    const double trace = l[0] + l[1] + l[2];
    rVariation = trace > 0.0 ? std::max(l[0], 0.0) / trace : 0.0;
    if (!(trace > 0.0) || !(l[1] - l[0] > 1e-12 * l[2]))
        return false;
    // The null vector of C - l0 I: the largest cross product of two rows.
    const PtnVec r0{t[0] - l[0], t[3], t[5]};
    const PtnVec r1{t[3], t[1] - l[0], t[4]};
    const PtnVec r2{t[5], t[4], t[2] - l[0]};
    const std::array<PtnVec, 3> cand{detail::vec3_cross(r0, r1), detail::vec3_cross(r0, r2),
                                     detail::vec3_cross(r1, r2)};
    std::size_t best = 0;
    for (std::size_t j = 1; j < 3; ++j)
        if (detail::vec3_dot(cand[j], cand[j]) > detail::vec3_dot(cand[best], cand[best]))
            best = j;
    const double len = std::sqrt(detail::vec3_dot(cand[best], cand[best]));
    if (!(len > 0.0))
        return false;
    rNormal = {cand[best][0] / len, cand[best][1] / len, cand[best][2] / len};
    return true;
}

}  // namespace

NormalOrientation normal_orientation_from_name(const std::string& rName) {
    if (rName == "consistent")
        return NormalOrientation::Consistent;
    if (rName == "viewpoint")
        return NormalOrientation::Viewpoint;
    if (rName == "none")
        return NormalOrientation::None;
    throw std::invalid_argument("meshio++: estimate_normals: unknown orientation '" + rName +
                                "' (expected 'consistent', 'viewpoint' or 'none')");
}

PointNormalsResult estimate_normals(const Mesh& rMesh, const PointNormalsOptions& rOptions) {
    const std::size_t n = rMesh.NumPoints();
    const std::size_t dim = rMesh.PointDim();
    if (n < 3)
        throw std::invalid_argument(
            "meshio++: estimate_normals: needs at least three points, not " + std::to_string(n));
    if (dim < 2)
        throw std::invalid_argument(
            "meshio++: estimate_normals: the points must have 2 or 3 coordinates");
    if (rOptions.mK < 2)
        throw std::invalid_argument("meshio++: estimate_normals: k must be at least 2, not " +
                                    std::to_string(rOptions.mK));
    std::vector<PtnVec> p(n, PtnVec{0.0, 0.0, 0.0});
    NDArray coords(DType::Float64, {n, std::size_t{3}});
    double* pc = coords.As<double>();
    const NDArray& points = rMesh.Points();
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t k = 0; k < dim && k < 3; ++k) {
            const double v = detail::read_double(points, i * dim + k);
            if (!std::isfinite(v))
                throw std::invalid_argument("meshio++: estimate_normals: point " +
                                            std::to_string(i) + " has a non-finite coordinate");
            p[i][k] = v;
            pc[i * 3 + k] = v;
        }

    const std::size_t k = static_cast<std::size_t>(
        std::min<std::int64_t>(rOptions.mK, static_cast<std::int64_t>(n) - 1));
    NeighborOptions search;
    search.mMethod = NeighborMethod::KNearest;
    search.mK = static_cast<std::int64_t>(k);
    {
        // The automatic lattice assumes the points fill their box; a scan
        // samples a surface, so size the cells for about k/2 points each over
        // an area estimated from the box (half its surface: a plane's area, a
        // sphere's within 5%). The cell size changes the search's speed,
        // never its answer.
        const detail::BoundingBox3 box = detail::BoundingBoxUtils::FromPoints(p);
        const double area = detail::BoundingBoxUtils::HalfSurfaceArea(box);
        const double extent = detail::BoundingBoxUtils::MaxExtent(box);
        if (area > 0.0)
            search.mCellSize =
                std::max(std::sqrt(area / static_cast<double>(n) * static_cast<double>(k) * 0.5),
                         extent * 1e-9);
    }
    const NeighborPairs pairs = neighbor_pairs(coords, search);
    const std::int64_t* pTarget = pairs.mTarget.As<std::int64_t>();
    if (pairs.mTarget.Size() != n * k)
        throw std::runtime_error("meshio++: estimate_normals: internal error: neighbour count");

    std::vector<PtnVec> normal(n, PtnVec{0.0, 0.0, 0.0});
    std::vector<double> variation(n, 0.0);
    std::vector<char> ok(n, 0);
    parallel_for(n, [&](std::size_t i) {
        PtnVec v;
        ok[i] = ptn_fit(p, i, pTarget + i * k, k, v, variation[i]) ? 1 : 0;
        if (ok[i]) {
            const double s = ptn_canonical_sign(v);
            normal[i] = {s * v[0], s * v[1], s * v[2]};
        }
    });

    PointNormalsResult result;
    if (rOptions.mOrientation == NormalOrientation::Viewpoint) {
        for (std::size_t i = 0; i < n; ++i) {
            const PtnVec to{rOptions.mViewpoint[0] - p[i][0], rOptions.mViewpoint[1] - p[i][1],
                            rOptions.mViewpoint[2] - p[i][2]};
            if (ok[i] && detail::vec3_dot(normal[i], to) < 0.0)
                normal[i] = {-normal[i][0], -normal[i][1], -normal[i][2]};
        }
    } else if (rOptions.mOrientation == NormalOrientation::Consistent) {
        // The symmetrised neighbour graph over the fitted points, as CSR:
        // each row the point's own neighbours and those that list it,
        // sorted and distinct.
        std::vector<std::int64_t> start(n + 1, 0);
        for (std::size_t i = 0; i < n; ++i) {
            if (!ok[i])
                continue;
            for (std::size_t j = 0; j < k; ++j) {
                const auto t = static_cast<std::size_t>(pTarget[i * k + j]);
                if (ok[t] && t != i) {
                    ++start[i + 1];
                    ++start[t + 1];
                }
            }
        }
        for (std::size_t i = 0; i < n; ++i)
            start[i + 1] += start[i];
        std::vector<std::int64_t> nbr(static_cast<std::size_t>(start[n]));
        {
            std::vector<std::int64_t> fill(start.begin(), start.end() - 1);
            for (std::size_t i = 0; i < n; ++i) {
                if (!ok[i])
                    continue;
                for (std::size_t j = 0; j < k; ++j) {
                    const auto t = static_cast<std::size_t>(pTarget[i * k + j]);
                    if (ok[t] && t != i) {
                        nbr[static_cast<std::size_t>(fill[i]++)] = static_cast<std::int64_t>(t);
                        nbr[static_cast<std::size_t>(fill[t]++)] = static_cast<std::int64_t>(i);
                    }
                }
            }
        }
        std::vector<std::int64_t> degree(n, 0);
        parallel_for(n, [&](std::size_t i) {
            const auto b = nbr.begin() + start[i], e = nbr.begin() + start[i + 1];
            std::sort(b, e);
            degree[i] = std::unique(b, e) - b;
        });
        const auto row = [&](std::size_t I) {
            return std::make_pair(nbr.data() + start[I], nbr.data() + start[I] + degree[I]);
        };
        std::vector<std::int64_t> piece(n, -1);
        std::vector<char> done(n, 0);
        std::vector<double> key(n, std::numeric_limits<double>::infinity());
        std::vector<std::int64_t> key_src(n, std::numeric_limits<std::int64_t>::max());
        for (std::size_t s0 = 0; s0 < n; ++s0) {
            if (!ok[s0] || piece[s0] >= 0)
                continue;
            // The piece, by a breadth-first sweep; its highest point seeds it.
            const auto id = result.mNumComponents++;
            std::vector<std::size_t> members{s0};
            piece[s0] = id;
            for (std::size_t q = 0; q < members.size(); ++q)
                for (auto [it, end] = row(members[q]); it != end; ++it)
                    if (piece[static_cast<std::size_t>(*it)] < 0) {
                        piece[static_cast<std::size_t>(*it)] = id;
                        members.push_back(static_cast<std::size_t>(*it));
                    }
            std::size_t seed = members.front();
            for (const std::size_t m : members)
                if (p[m][2] > p[seed][2] || (p[m][2] == p[seed][2] && m < seed))
                    seed = m;
            // Prim's tree from the seed, cheapest edge first (1 - |n_i . n_j|),
            // ties to the lower target, then the lower source. An edge is
            // queued only when it betters its target's best so far, which
            // pops the same tree as queueing every edge.
            using Edge = std::tuple<double, std::int64_t, std::int64_t>;
            std::priority_queue<Edge, std::vector<Edge>, std::greater<Edge>> heap;
            const auto push_from = [&](std::size_t S) {
                for (auto [it, end] = row(S); it != end; ++it) {
                    const auto t = static_cast<std::size_t>(*it);
                    if (done[t])
                        continue;
                    const double w = 1.0 - std::abs(detail::vec3_dot(normal[S], normal[t]));
                    const auto src = static_cast<std::int64_t>(S);
                    if (w < key[t] || (w == key[t] && src < key_src[t])) {
                        key[t] = w;
                        key_src[t] = src;
                        heap.emplace(w, *it, src);
                    }
                }
            };
            done[seed] = 1;
            push_from(seed);
            while (!heap.empty()) {
                const auto [w, t, s] = heap.top();
                heap.pop();
                const auto ti = static_cast<std::size_t>(t);
                if (done[ti] || w != key[ti] || s != key_src[ti])
                    continue;
                done[ti] = 1;
                if (detail::vec3_dot(normal[ti], normal[static_cast<std::size_t>(s)]) < 0.0)
                    normal[ti] = {-normal[ti][0], -normal[ti][1], -normal[ti][2]};
                push_from(ti);
            }
            // The whole piece outward from its centroid, if that decides;
            // else its seed's normal up, else its canonical sign.
            PtnVec c{0.0, 0.0, 0.0};
            for (const std::size_t m : members)
                for (std::size_t a = 0; a < 3; ++a)
                    c[a] += p[m][a];
            for (std::size_t a = 0; a < 3; ++a)
                c[a] /= static_cast<double>(members.size());
            double outward = 0.0, spread = 0.0;
            for (const std::size_t m : members) {
                const PtnVec d{p[m][0] - c[0], p[m][1] - c[1], p[m][2] - c[2]};
                outward += detail::vec3_dot(normal[m], d);
                spread += std::sqrt(detail::vec3_dot(d, d));
            }
            double flip = 1.0;
            if (std::abs(outward) > 1e-6 * spread)
                flip = outward < 0.0 ? -1.0 : 1.0;
            else if (std::abs(normal[seed][2]) > 1e-9)
                flip = normal[seed][2] < 0.0 ? -1.0 : 1.0;
            else
                flip = ptn_canonical_sign(normal[seed]);
            if (flip < 0.0)
                for (const std::size_t m : members)
                    normal[m] = {-normal[m][0], -normal[m][1], -normal[m][2]};
        }
    }

    NDArray out(DType::Float64, {n, std::size_t{3}});
    double* po = out.As<double>();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    for (std::size_t i = 0; i < n; ++i) {
        if (!ok[i])
            ++result.mNumDegenerate;
        for (std::size_t a = 0; a < 3; ++a)
            po[i * 3 + a] = ok[i] ? normal[i][a] : nan;
    }
    result.mMesh = detail::clone_mesh(rMesh);
    result.mMesh.AddPointData(kNormalsName, std::move(out));
    if (rOptions.mRecordVariation) {
        NDArray var(DType::Float64, {n});
        std::copy(variation.begin(), variation.end(), var.As<double>());
        result.mMesh.AddPointData(kNormalsVariationName, std::move(var));
    }
    return result;
}

}  // namespace meshioplusplus

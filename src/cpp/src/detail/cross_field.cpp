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
#include <complex>
#include <map>

// Project includes
#include "cross_field.hpp"
#include "sparse_solve.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

constexpr double kCfPi = 3.14159265358979323846;

using CfComplex = std::complex<double>;

PatchVec cf_sub(const PatchVec& rA, const PatchVec& rB) {
    return {rA[0] - rB[0], rA[1] - rB[1], rA[2] - rB[2]};
}
double cf_dot(const PatchVec& rA, const PatchVec& rB) {
    return rA[0] * rB[0] + rA[1] * rB[1] + rA[2] * rB[2];
}
PatchVec cf_cross(const PatchVec& rA, const PatchVec& rB) {
    return {rA[1] * rB[2] - rA[2] * rB[1], rA[2] * rB[0] - rA[0] * rB[2],
            rA[0] * rB[1] - rA[1] * rB[0]};
}
PatchVec cf_unit(const PatchVec& rA) {
    const double l = std::sqrt(cf_dot(rA, rA));
    return l > 0.0 ? PatchVec{rA[0] / l, rA[1] / l, rA[2] / l} : PatchVec{0.0, 0.0, 0.0};
}

// The angle of the 3-D vector `rD` in triangle F's frame.
double cf_angle(const CrossField& rF, std::size_t F, const PatchVec& rD) {
    return std::atan2(cf_dot(rD, rF.mE2[F]), cf_dot(rD, rF.mE1[F]));
}

// `A` wrapped into (-pi/4, pi/4].
double cf_wrap_quarter(double A) {
    const double q = 0.5 * kCfPi;
    A = std::fmod(A, q);
    if (A > 0.5 * q)
        A -= q;
    else if (A <= -0.5 * q)
        A += q;
    return A;
}

using CfEdge = std::pair<std::size_t, std::size_t>;

CfEdge cf_edge(std::size_t A, std::size_t B) {
    return {std::min(A, B), std::max(A, B)};
}

}  // namespace

PatchVec cross_branch(const CrossField& rField, std::size_t F, int K) {
    const double t = rField.mTheta[F] + 0.5 * kCfPi * static_cast<double>(((K % 4) + 4) % 4);
    const double c = std::cos(t), s = std::sin(t);
    const PatchVec& e1 = rField.mE1[F];
    const PatchVec& e2 = rField.mE2[F];
    return {c * e1[0] + s * e2[0], c * e1[1] + s * e2[1], c * e1[2] + s * e2[2]};
}

CrossField cross_field(const std::vector<PatchVec>& rPts, const std::vector<PatchTri>& rTris,
                       const CrossFieldOptions& rOptions) {
    CrossField out;
    const std::size_t nt = rTris.size();
    if (nt == 0)
        return out;
    out.mE1.resize(nt);
    out.mE2.resize(nt);
    out.mNormal.resize(nt);
    out.mTheta.assign(nt, 0.0);
    std::vector<double> area(nt, 0.0);
    for (std::size_t f = 0; f < nt; ++f) {
        const PatchVec& a = rPts[rTris[f][0]];
        const PatchVec u = cf_sub(rPts[rTris[f][1]], a), v = cf_sub(rPts[rTris[f][2]], a);
        const PatchVec n = cf_cross(u, v);
        area[f] = 0.5 * std::sqrt(cf_dot(n, n));
        out.mNormal[f] = cf_unit(n);
        out.mE1[f] = cf_unit(u);
        out.mE2[f] = cf_cross(out.mNormal[f], out.mE1[f]);
    }

    // Edges and the triangles beside each.
    std::map<CfEdge, std::vector<std::size_t>> edge_tris;
    for (std::size_t f = 0; f < nt; ++f)
        for (int k = 0; k < 3; ++k)
            edge_tris[cf_edge(rTris[f][static_cast<std::size_t>(k)],
                              rTris[f][static_cast<std::size_t>((k + 1) % 3)])]
                .push_back(f);
    std::vector<CfEdge> features = rOptions.mFeatureEdges;
    for (auto& e : features)
        e = cf_edge(e.first, e.second);
    std::sort(features.begin(), features.end());

    // Constraints: a triangle on a boundary or feature edge holds that edge's
    // direction (the sum of its constrained edges' crosses, when it has more).
    std::vector<CfComplex> fixed(nt, CfComplex(0.0, 0.0));
    std::vector<std::uint8_t> is_fixed(nt, 0);
    for (const auto& kv : edge_tris) {
        const bool boundary = kv.second.size() == 1;
        const bool feature = std::binary_search(features.begin(), features.end(), kv.first);
        if (!(boundary && rOptions.mConstrainBoundary) && !feature)
            continue;
        const PatchVec d = cf_sub(rPts[kv.first.second], rPts[kv.first.first]);
        for (std::size_t f : kv.second) {
            fixed[f] += std::polar(1.0, 4.0 * cf_angle(out, f, d));
            is_fixed[f] = 1;
        }
    }
    for (std::size_t f = 0; f < nt; ++f)
        if (is_fixed[f]) {
            const double m = std::abs(fixed[f]);
            fixed[f] = m > 1e-9 ? fixed[f] / m : CfComplex(1.0, 0.0);
            ++out.mNumConstrained;
        }

    // The soft pull towards the principal curvature directions: each
    // triangle's bending tensor sum_e beta_e |e| t_e t_e^T over its edges
    // (beta the dihedral angle, t the in-plane direction across the edge).
    std::vector<CfComplex> pull(nt, CfComplex(0.0, 0.0));
    std::vector<double> pull_w(nt, 0.0);
    if (rOptions.mCurvatureWeight > 0.0)
        for (std::size_t f = 0; f < nt; ++f) {
            double t00 = 0.0, t01 = 0.0, t11 = 0.0;
            for (int k = 0; k < 3; ++k) {
                const std::size_t a = rTris[f][static_cast<std::size_t>(k)];
                const std::size_t b = rTris[f][static_cast<std::size_t>((k + 1) % 3)];
                const auto& adj = edge_tris[cf_edge(a, b)];
                if (adj.size() != 2)
                    continue;
                const std::size_t g = adj[0] == f ? adj[1] : adj[0];
                const double beta = std::acos(
                    std::max(-1.0, std::min(1.0, cf_dot(out.mNormal[f], out.mNormal[g]))));
                const PatchVec e = cf_sub(rPts[b], rPts[a]);
                const double len = std::sqrt(cf_dot(e, e));
                const PatchVec across = cf_unit(cf_cross(out.mNormal[f], e));
                const double x = cf_dot(across, out.mE1[f]), y = cf_dot(across, out.mE2[f]);
                t00 += beta * len * x * x;
                t01 += beta * len * x * y;
                t11 += beta * len * y * y;
            }
            const double tr = t00 + t11;
            const double diff = std::sqrt((t00 - t11) * (t00 - t11) + 4.0 * t01 * t01);
            if (!(tr > 0.0))
                continue;
            // The largest eigenvector's angle: 2 phi = atan2(2 t01, t00 - t11).
            const double phi = 0.5 * std::atan2(2.0 * t01, t00 - t11);
            pull[f] = std::polar(1.0, 4.0 * phi);
            pull_w[f] = rOptions.mCurvatureWeight * diff / tr;
        }

    // The free triangles, numbered, and the Hermitian system as its real twin.
    std::vector<std::size_t> slot(nt, nt);
    std::size_t nfree = 0;
    for (std::size_t f = 0; f < nt; ++f)
        if (!is_fixed[f])
            slot[f] = nfree++;
    std::vector<SparseTriplet> trip;
    std::vector<double> rhs(2 * nfree, 0.0);
    const auto add = [&](std::size_t R, std::size_t S, CfComplex H) {
        trip.push_back({2 * R, 2 * S, H.real()});
        trip.push_back({2 * R, 2 * S + 1, -H.imag()});
        trip.push_back({2 * R + 1, 2 * S, H.imag()});
        trip.push_back({2 * R + 1, 2 * S + 1, H.real()});
    };
    bool any_anchor = false;
    for (const auto& kv : edge_tris) {
        if (kv.second.size() != 2)
            continue;
        const std::size_t f = kv.second[0], g = kv.second[1];
        const PatchVec d = cf_sub(rPts[kv.first.second], rPts[kv.first.first]);
        const double af = cf_angle(out, f, d), ag = cf_angle(out, g, d);
        // z_f ~ c z_g: the cross carried from g into f's frame.
        const CfComplex c = std::polar(1.0, 4.0 * (af - ag));
        const bool ff = !is_fixed[f], fg = !is_fixed[g];
        if (ff && fg) {
            add(slot[f], slot[f], 1.0);
            add(slot[g], slot[g], 1.0);
            add(slot[f], slot[g], -c);
            add(slot[g], slot[f], -std::conj(c));
        } else if (ff) {
            add(slot[f], slot[f], 1.0);
            const CfComplex b = c * fixed[g];
            rhs[2 * slot[f]] += b.real();
            rhs[2 * slot[f] + 1] += b.imag();
            any_anchor = true;
        } else if (fg) {
            add(slot[g], slot[g], 1.0);
            const CfComplex b = std::conj(c) * fixed[f];
            rhs[2 * slot[g]] += b.real();
            rhs[2 * slot[g] + 1] += b.imag();
            any_anchor = true;
        }
    }
    for (std::size_t f = 0; f < nt; ++f)
        if (!is_fixed[f] && pull_w[f] > 0.0) {
            add(slot[f], slot[f], pull_w[f]);
            rhs[2 * slot[f]] += pull_w[f] * pull[f].real();
            rhs[2 * slot[f] + 1] += pull_w[f] * pull[f].imag();
            any_anchor = true;
        }
    std::vector<CfComplex> z(nt, CfComplex(1.0, 0.0));
    for (std::size_t f = 0; f < nt; ++f)
        if (is_fixed[f])
            z[f] = fixed[f];
    if (nfree > 0) {
        std::vector<double> x(2 * nfree, 0.0);
        if (any_anchor) {
            // A tiny shift keeps a free component the anchors do not reach
            // solvable; it does not move an anchored one measurably.
            for (std::size_t i = 0; i < 2 * nfree; ++i)
                trip.push_back({i, i, 1e-10});
            const SparseCsr a = csr_from_triplets(2 * nfree, 2 * nfree, trip);
            csr_pcg(a, rhs.data(), x.data(), 1e-12, 20 * nfree + 100);
        } else {
            // Nothing pins the field: the smallest eigenvector, by inverse
            // iteration from a fixed start.
            for (std::size_t i = 0; i < 2 * nfree; ++i)
                trip.push_back({i, i, 1e-8});
            const SparseCsr a = csr_from_triplets(2 * nfree, 2 * nfree, trip);
            std::vector<double> v(2 * nfree, 0.0);
            for (std::size_t i = 0; i < nfree; ++i)
                v[2 * i] = 1.0;
            for (int it = 0; it < rOptions.mPowerIterations; ++it) {
                std::vector<double> y(2 * nfree, 0.0);
                csr_pcg(a, v.data(), y.data(), 1e-12, 20 * nfree + 100);
                const double n = std::sqrt(sparse_dot(y.data(), y.data(), y.size()));
                if (!(n > 0.0))
                    break;
                for (std::size_t i = 0; i < y.size(); ++i)
                    v[i] = y[i] / n;
            }
            x = v;
        }
        for (std::size_t f = 0; f < nt; ++f)
            if (!is_fixed[f])
                z[f] = CfComplex(x[2 * slot[f]], x[2 * slot[f] + 1]);
    }
    for (std::size_t f = 0; f < nt; ++f) {
        if (!(std::abs(z[f]) > 1e-12)) {
            // A triangle no constraint reaches keeps branch 0 along its edge.
            out.mTheta[f] = 0.0;
            continue;
        }
        out.mTheta[f] = 0.25 * std::arg(z[f]);
    }

    // Indices: the cross's turning around each interior vertex.
    std::size_t npts = 0;
    for (const PatchTri& t : rTris)
        npts = std::max({npts, t[0] + 1, t[1] + 1, t[2] + 1});
    out.mIndex.assign(std::max(npts, rPts.size()), 0);
    std::vector<std::vector<std::size_t>> vert_tris(out.mIndex.size());
    for (std::size_t f = 0; f < nt; ++f)
        for (std::size_t v : rTris[f])
            vert_tris[v].push_back(f);
    for (std::size_t v = 0; v < vert_tris.size(); ++v) {
        const std::vector<std::size_t>& ring = vert_tris[v];
        if (ring.size() < 3)
            continue;
        // Walk the fan: from a triangle, across its edge (v, the vertex before
        // v), to the triangle beyond.
        std::size_t f = ring[0];
        double turn = 0.0, angles = 0.0;
        bool closed = false;
        for (std::size_t step = 0; step <= ring.size(); ++step) {
            const PatchTri& t = rTris[f];
            const int iv = t[0] == v ? 0 : (t[1] == v ? 1 : 2);
            const std::size_t prev = t[static_cast<std::size_t>((iv + 2) % 3)];
            const std::size_t next = t[static_cast<std::size_t>((iv + 1) % 3)];
            const PatchVec a = cf_unit(cf_sub(rPts[next], rPts[v]));
            const PatchVec b = cf_unit(cf_sub(rPts[prev], rPts[v]));
            angles += std::acos(std::max(-1.0, std::min(1.0, cf_dot(a, b))));
            const auto& adj = edge_tris[cf_edge(v, prev)];
            if (adj.size() != 2)
                break;  // a boundary vertex
            const std::size_t g = adj[0] == f ? adj[1] : adj[0];
            const PatchVec d = cf_sub(rPts[prev], rPts[v]);
            const double carried = out.mTheta[f] + cf_angle(out, g, d) - cf_angle(out, f, d);
            turn += cf_wrap_quarter(out.mTheta[g] - carried);
            f = g;
            if (f == ring[0]) {
                closed = step + 1 == ring.size();
                break;
            }
        }
        if (!closed)
            continue;
        const double index = (turn + (2.0 * kCfPi - angles)) / (2.0 * kCfPi);
        out.mIndex[v] = static_cast<int>(std::lround(4.0 * index));
        if (out.mIndex[v] != 0)
            out.mSingular.push_back(v);
    }
    out.mOk = true;
    return out;
}

}  // namespace detail
}  // namespace meshioplusplus

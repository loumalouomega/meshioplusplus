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
#include <limits>
#include <map>
#include <queue>
#include <utility>

// Project includes
#include "feature_fit.hpp"
#include "lsq_fit.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

constexpr double kFfPi = 3.14159265358979323846;

double ff_dot(const PatchVec& rA, const PatchVec& rB) {
    return rA[0] * rB[0] + rA[1] * rB[1] + rA[2] * rB[2];
}

PatchVec ff_sub(const PatchVec& rA, const PatchVec& rB) {
    return {rA[0] - rB[0], rA[1] - rB[1], rA[2] - rB[2]};
}

PatchVec ff_axpy(const PatchVec& rA, double S, const PatchVec& rB) {
    return {rA[0] + S * rB[0], rA[1] + S * rB[1], rA[2] + S * rB[2]};
}

double ff_norm(const PatchVec& rA) {
    return std::sqrt(ff_dot(rA, rA));
}

// The vertex's height along and distance from the axis through pP[0..2]
// along the unit pP[3..5].
void ff_axial(const double* pP, const PatchVec& rX, double& rH, double& rRho) {
    const PatchVec d{rX[0] - pP[0], rX[1] - pP[1], rX[2] - pP[2]};
    const PatchVec a{pP[3], pP[4], pP[5]};
    rH = ff_dot(d, a);
    rRho = ff_norm(ff_axpy(d, -rH, a));
}

// Outward unit vertex normals (area-weighted) and which vertices no
// boundary edge touches.
void ff_vertex_normals(const std::vector<PatchVec>& rPts, const std::vector<PatchTri>& rTris,
                       std::vector<PatchVec>& rNormals, std::vector<std::uint8_t>& rInterior) {
    rNormals.assign(rPts.size(), PatchVec{0.0, 0.0, 0.0});
    std::map<std::pair<std::size_t, std::size_t>, int> uses;
    for (const PatchTri& t : rTris) {
        const PatchVec e1 = ff_sub(rPts[t[1]], rPts[t[0]]);
        const PatchVec e2 = ff_sub(rPts[t[2]], rPts[t[0]]);
        const PatchVec n{e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2],
                         e1[0] * e2[1] - e1[1] * e2[0]};
        for (std::size_t v : t)
            rNormals[v] = ff_axpy(rNormals[v], 1.0, n);
        for (int k = 0; k < 3; ++k) {
            const std::size_t a = t[static_cast<std::size_t>(k)];
            const std::size_t b = t[static_cast<std::size_t>((k + 1) % 3)];
            ++uses[{std::min(a, b), std::max(a, b)}];
        }
    }
    for (PatchVec& n : rNormals) {
        const double len = ff_norm(n);
        if (len > 0.0)
            n = {n[0] / len, n[1] / len, n[2] / len};
    }
    rInterior.assign(rPts.size(), 0);
    for (const PatchTri& t : rTris)
        for (std::size_t v : t)
            rInterior[v] = 1;
    for (const auto& kv : uses)
        if (kv.second == 1)
            rInterior[kv.first.first] = rInterior[kv.first.second] = 0;
}

}  // namespace

double primitive_distance(int Type, const double* pP, const PatchVec& rX) {
    switch (Type) {
        case 1: {
            const PatchVec d{rX[0] - pP[0], rX[1] - pP[1], rX[2] - pP[2]};
            return std::fabs(d[0] * pP[3] + d[1] * pP[4] + d[2] * pP[5]);
        }
        case 2: {
            const PatchVec d{rX[0] - pP[0], rX[1] - pP[1], rX[2] - pP[2]};
            return std::fabs(ff_norm(d) - pP[6]);
        }
        case 3: {
            double h = 0.0, rho = 0.0;
            ff_axial(pP, rX, h, rho);
            return std::fabs(rho - pP[6]);
        }
        case 4: {
            double h = 0.0, rho = 0.0;
            ff_axial(pP, rX, h, rho);
            const double phi = pP[8] * kFfPi / 180.0;
            if (h * std::cos(phi) + rho * std::sin(phi) < 0.0)
                return std::sqrt(h * h + rho * rho);  // behind the apex
            return std::fabs(rho * std::cos(phi) - h * std::sin(phi));
        }
        case 5: {
            double h = 0.0, rho = 0.0;
            ff_axial(pP, rX, h, rho);
            return std::fabs(std::hypot(rho - pP[6], h) - pP[7]);
        }
        default:
            return std::numeric_limits<double>::infinity();
    }
}

FilletFit fillet_fit(const std::vector<PatchVec>& rPts, const std::vector<PatchTri>& rTris,
                     int Type1, const double* pParams1, int Type2, const double* pParams2,
                     double RMin, double RMax) {
    FilletFit out;
    if (!(RMin > 0.0) || !(RMax > RMin))
        return out;
    std::vector<PatchVec> normals;
    std::vector<std::uint8_t> interior;
    ff_vertex_normals(rPts, rTris, normals, interior);
    std::vector<std::size_t> verts;
    for (std::size_t v = 0; v < rPts.size(); ++v)
        if (interior[v] && ff_norm(normals[v]) > 0.5)
            verts.push_back(v);
    if (verts.size() < 4)
        return out;
    auto misfit = [&](int Side, double R) {
        double s = 0.0;
        for (std::size_t v : verts) {
            const PatchVec c = ff_axpy(rPts[v], -Side * R, normals[v]);
            const double d1 = primitive_distance(Type1, pParams1, c) - R;
            const double d2 = primitive_distance(Type2, pParams2, c) - R;
            s += d1 * d1 + d2 * d2;
        }
        return s / static_cast<double>(2 * verts.size());
    };
    constexpr int kScan = 48;
    const double lr0 = std::log(RMin), lr1 = std::log(RMax);
    double best = std::numeric_limits<double>::infinity();
    for (int side : {1, -1}) {
        int arg = 0;
        double low = std::numeric_limits<double>::infinity();
        for (int i = 0; i < kScan; ++i) {
            const double f = misfit(side, std::exp(lr0 + (lr1 - lr0) * i / (kScan - 1)));
            if (f < low) {
                low = f;
                arg = i;
            }
        }
        // Golden section between the best sample's neighbours.
        double a = lr0 + (lr1 - lr0) * std::max(arg - 1, 0) / (kScan - 1);
        double b = lr0 + (lr1 - lr0) * std::min(arg + 1, kScan - 1) / (kScan - 1);
        const double g = 0.5 * (std::sqrt(5.0) - 1.0);
        double x1 = b - g * (b - a), x2 = a + g * (b - a);
        double f1 = misfit(side, std::exp(x1)), f2 = misfit(side, std::exp(x2));
        for (int it = 0; it < 80; ++it) {
            if (f1 <= f2) {
                b = x2;
                x2 = x1;
                f2 = f1;
                x1 = b - g * (b - a);
                f1 = misfit(side, std::exp(x1));
            } else {
                a = x1;
                x1 = x2;
                f1 = f2;
                x2 = a + g * (b - a);
                f2 = misfit(side, std::exp(x2));
            }
        }
        const double r = std::exp(0.5 * (a + b));
        const double f = misfit(side, r);
        if (f < best) {
            best = f;
            out.mRadius = r;
            out.mSide = side;
        }
    }
    out.mRms = std::sqrt(best);
    out.mOk = std::isfinite(out.mRms);
    for (std::size_t v : verts)
        out.mCentres.push_back(ff_axpy(rPts[v], -out.mSide * out.mRadius, normals[v]));
    return out;
}

PatchVec HelicoidFit::At(double Rho, double Theta) const {
    PatchVec p = ff_axpy(mOrigin, Rho * std::cos(Theta), mE1);
    p = ff_axpy(p, Rho * std::sin(Theta), mE2);
    return ff_axpy(p, mC + mK * Theta + mT * Rho, mAxis);
}

double HelicoidFit::Pitch() const {
    return 2.0 * kFfPi * mK;
}

HelicoidFit helicoid_fit(const std::vector<PatchVec>& rPts, const std::vector<PatchTri>& rTris,
                         const PatchVec& rOrigin, const PatchVec& rAxis) {
    HelicoidFit out;
    const std::size_t n = rPts.size();
    const double alen = ff_norm(rAxis);
    if (n < 4 || rTris.empty() || !(alen > 0.0))
        return out;
    out.mOrigin = rOrigin;
    out.mAxis = {rAxis[0] / alen, rAxis[1] / alen, rAxis[2] / alen};
    lsq_frame(out.mAxis.data(), out.mE1.data(), out.mE2.data());
    std::vector<double> z(n), rho(n), phi(n);
    double rho_max = 0.0;
    for (std::size_t v = 0; v < n; ++v) {
        const PatchVec d = ff_sub(rPts[v], rOrigin);
        z[v] = ff_dot(d, out.mAxis);
        const PatchVec radial = ff_axpy(d, -z[v], out.mAxis);
        rho[v] = ff_norm(radial);
        phi[v] = std::atan2(ff_dot(radial, out.mE2), ff_dot(radial, out.mE1));
        rho_max = std::max(rho_max, rho[v]);
    }
    for (double r : rho)
        if (!(r > 1e-9 * rho_max))
            return out;  // the patch touches the axis
    // Unwrap the angle over the patch's edges.
    std::vector<std::vector<std::size_t>> adj(n);
    for (const PatchTri& t : rTris)
        for (int k = 0; k < 3; ++k) {
            adj[t[static_cast<std::size_t>(k)]].push_back(t[static_cast<std::size_t>((k + 1) % 3)]);
            adj[t[static_cast<std::size_t>((k + 1) % 3)]].push_back(t[static_cast<std::size_t>(k)]);
        }
    for (auto& a : adj)
        std::sort(a.begin(), a.end());
    std::vector<double> theta(n, 0.0);
    std::vector<std::uint8_t> seen(n, 0);
    std::queue<std::size_t> todo;
    std::size_t start = 0;
    while (start < n && adj[start].empty())
        ++start;
    if (start == n)
        return out;
    theta[start] = phi[start];
    seen[start] = 1;
    todo.push(start);
    while (!todo.empty()) {
        const std::size_t u = todo.front();
        todo.pop();
        for (std::size_t v : adj[u]) {
            double step = phi[v] - phi[u];
            step -= 2.0 * kFfPi * std::round(step / (2.0 * kFfPi));
            if (!seen[v]) {
                seen[v] = 1;
                theta[v] = theta[u] + step;
                todo.push(v);
            } else if (std::fabs(theta[v] - theta[u] - step) > 1e-6) {
                return out;  // the angle does not unwrap: the patch winds round the axis
            }
        }
    }
    for (std::size_t v = 0; v < n; ++v)
        if (!seen[v] && !adj[v].empty())
            return out;  // more than one piece
    // z = C + K theta + T rho, by least squares.
    std::vector<double> a(9, 0.0), b(3, 0.0);
    for (std::size_t v = 0; v < n; ++v) {
        if (adj[v].empty())
            continue;
        const double row[3] = {1.0, theta[v], rho[v]};
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j)
                a[static_cast<std::size_t>(i * 3 + j)] += row[i] * row[j];
            b[static_cast<std::size_t>(i)] += row[i] * z[v];
        }
    }
    if (!lsq_solve(a, b, 3))
        return out;
    out.mC = b[0];
    out.mK = b[1];
    out.mT = b[2];
    out.mThetaMin = out.mRhoMin = std::numeric_limits<double>::infinity();
    out.mThetaMax = out.mRhoMax = -std::numeric_limits<double>::infinity();
    double sum = 0.0;
    std::size_t used = 0;
    for (std::size_t v = 0; v < n; ++v) {
        if (adj[v].empty())
            continue;
        const double e = z[v] - (out.mC + out.mK * theta[v] + out.mT * rho[v]);
        const double grad =
            std::sqrt(1.0 + out.mT * out.mT + (out.mK / rho[v]) * (out.mK / rho[v]));
        const double d = std::fabs(e) / grad;
        sum += d * d;
        out.mMax = std::max(out.mMax, d);
        ++used;
        out.mThetaMin = std::min(out.mThetaMin, theta[v]);
        out.mThetaMax = std::max(out.mThetaMax, theta[v]);
        out.mRhoMin = std::min(out.mRhoMin, rho[v]);
        out.mRhoMax = std::max(out.mRhoMax, rho[v]);
    }
    out.mRms = std::sqrt(sum / static_cast<double>(used));
    // Half a turn at least: less of a helicoid is as well a cone or a plane.
    out.mOk = out.mThetaMax - out.mThetaMin >= kFfPi && std::isfinite(out.mRms) &&
              out.mRhoMax > out.mRhoMin;
    return out;
}

}  // namespace detail
}  // namespace meshioplusplus

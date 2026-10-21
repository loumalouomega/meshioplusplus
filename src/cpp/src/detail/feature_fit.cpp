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
#include <functional>
#include <limits>
#include <map>
#include <queue>
#include <utility>

// Project includes
#include "dense_chol.hpp"
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
    double z = mC + mK * Theta;
    if (mQ != 0.0)
        z += mQ * (Theta - mThetaRef) * (Theta - mThetaRef);
    return ff_axpy(p, z + mT * Rho, mAxis);
}

double HelicoidFit::Pitch() const {
    return 2.0 * kFfPi * mK;
}

double HelicoidFit::PitchAt(double Theta) const {
    return 2.0 * kFfPi * (mK + 2.0 * mQ * (Theta - mThetaRef));
}

HelicoidFit helicoid_fit(const std::vector<PatchVec>& rPts, const std::vector<PatchTri>& rTris,
                         const PatchVec& rOrigin, const PatchVec& rAxis, bool VariablePitch) {
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
    // z = C + K theta + T rho (+ Q (theta - ThetaRef)^2), by least squares.
    if (VariablePitch) {
        double sum = 0.0;
        std::size_t count = 0;
        for (std::size_t v = 0; v < n; ++v)
            if (!adj[v].empty()) {
                sum += theta[v];
                ++count;
            }
        out.mThetaRef = sum / static_cast<double>(count);
    }
    const int cols = VariablePitch ? 4 : 3;
    std::vector<double> a(static_cast<std::size_t>(cols * cols), 0.0),
        b(static_cast<std::size_t>(cols), 0.0);
    for (std::size_t v = 0; v < n; ++v) {
        if (adj[v].empty())
            continue;
        const double dt = theta[v] - out.mThetaRef;
        const double row[4] = {1.0, theta[v], rho[v], dt * dt};
        for (int i = 0; i < cols; ++i) {
            for (int j = 0; j < cols; ++j)
                a[static_cast<std::size_t>(i * cols + j)] += row[i] * row[j];
            b[static_cast<std::size_t>(i)] += row[i] * z[v];
        }
    }
    if (!lsq_solve(a, b, static_cast<std::size_t>(cols)))
        return out;
    out.mC = b[0];
    out.mK = b[1];
    out.mT = b[2];
    if (VariablePitch)
        out.mQ = b[3];
    out.mThetaMin = out.mRhoMin = std::numeric_limits<double>::infinity();
    out.mThetaMax = out.mRhoMax = -std::numeric_limits<double>::infinity();
    double sum = 0.0;
    std::size_t used = 0;
    for (std::size_t v = 0; v < n; ++v) {
        if (adj[v].empty())
            continue;
        double e = z[v] - (out.mC + out.mK * theta[v] + out.mT * rho[v]);
        double dz = out.mK;  // the axial term's derivative in theta
        if (VariablePitch) {
            const double dt = theta[v] - out.mThetaRef;
            e -= out.mQ * dt * dt;
            dz += 2.0 * out.mQ * dt;
        }
        const double grad = std::sqrt(1.0 + out.mT * out.mT + (dz / rho[v]) * (dz / rho[v]));
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

// ---------------------------------------------------------------------------
// Laws, rolling-ball blends of varying radius between any two surfaces, and
// helical surfaces of any profile (roadmap 6.3.1).
// ---------------------------------------------------------------------------

namespace {

PatchVec ff_cross(const PatchVec& rA, const PatchVec& rB) {
    return {rA[1] * rB[2] - rA[2] * rB[1], rA[2] * rB[0] - rA[0] * rB[2],
            rA[0] * rB[1] - rA[1] * rB[0]};
}

PatchVec ff_unit(const PatchVec& rA) {
    const double len = ff_norm(rA);
    return len > 0.0 ? PatchVec{rA[0] / len, rA[1] / len, rA[2] / len} : PatchVec{0.0, 0.0, 0.0};
}

std::size_t ff_law_size(bool Periodic, int Degree) {
    return Periodic ? static_cast<std::size_t>(2 * Degree + 1)
                    : static_cast<std::size_t>(Degree + 1);
}

// The law's basis functions at T, and their derivatives in T when pD is set.
void ff_law_basis(bool Periodic, double T0, double T1, std::size_t N, double T, double* pB,
                  double* pD) {
    if (Periodic) {
        const double w = 2.0 * kFfPi / (T1 - T0), x = w * (T - T0);
        pB[0] = 1.0;
        if (pD)
            pD[0] = 0.0;
        for (std::size_t j = 1; 2 * j - 1 < N; ++j) {
            const double jd = static_cast<double>(j);
            const double c = std::cos(jd * x), s = std::sin(jd * x);
            pB[2 * j - 1] = c;
            if (pD)
                pD[2 * j - 1] = -jd * w * s;
            if (2 * j < N) {
                pB[2 * j] = s;
                if (pD)
                    pD[2 * j] = jd * w * c;
            }
        }
        return;
    }
    const double du = 2.0 / (T1 - T0), u = du * (T - T0) - 1.0;
    pB[0] = 1.0;
    if (pD)
        pD[0] = 0.0;
    if (N > 1) {
        pB[1] = u;
        if (pD)
            pD[1] = du;
    }
    for (std::size_t k = 2; k < N; ++k) {
        pB[k] = 2.0 * u * pB[k - 1] - pB[k - 2];
        if (pD)
            pD[k] = 2.0 * du * pB[k - 1] + 2.0 * u * pD[k - 1] - pD[k - 2];
    }
}

// Every edge of a triangle set, counted: an edge used once is a boundary edge.
std::map<std::pair<std::size_t, std::size_t>, int> ff_edge_uses(
    const std::vector<PatchTri>& rTris) {
    std::map<std::pair<std::size_t, std::size_t>, int> uses;
    for (const PatchTri& t : rTris)
        for (int k = 0; k < 3; ++k) {
            const std::size_t a = t[static_cast<std::size_t>(k)];
            const std::size_t b = t[static_cast<std::size_t>((k + 1) % 3)];
            ++uses[{std::min(a, b), std::max(a, b)}];
        }
    return uses;
}

// At most Cap entries of rIn, evenly spread by position.
std::vector<std::size_t> ff_spread(const std::vector<std::size_t>& rIn, std::size_t Cap) {
    if (rIn.size() <= Cap)
        return rIn;
    std::vector<std::size_t> out;
    out.reserve(Cap);
    for (std::size_t i = 0; i < Cap; ++i)
        out.push_back(rIn[i * rIn.size() / Cap]);
    return out;
}

// The radius in [RMin, RMax] minimising F, by a log-spaced scan refined by
// golden section between the best sample's neighbours.
double ff_min_radius(const std::function<double(double)>& rF, double RMin, double RMax, int Scan,
                     int Iterations, double& rBest) {
    const double lr0 = std::log(RMin), lr1 = std::log(RMax);
    int arg = 0;
    double low = std::numeric_limits<double>::infinity();
    for (int i = 0; i < Scan; ++i) {
        const double f = rF(std::exp(lr0 + (lr1 - lr0) * i / (Scan - 1)));
        if (f < low) {
            low = f;
            arg = i;
        }
    }
    double a = lr0 + (lr1 - lr0) * std::max(arg - 1, 0) / (Scan - 1);
    double b = lr0 + (lr1 - lr0) * std::min(arg + 1, Scan - 1) / (Scan - 1);
    const double g = 0.5 * (std::sqrt(5.0) - 1.0);
    double x1 = b - g * (b - a), x2 = a + g * (b - a);
    double f1 = rF(std::exp(x1)), f2 = rF(std::exp(x2));
    for (int it = 0; it < Iterations; ++it) {
        if (f1 <= f2) {
            b = x2;
            x2 = x1;
            f2 = f1;
            x1 = b - g * (b - a);
            f1 = rF(std::exp(x1));
        } else {
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = a + g * (b - a);
            f2 = rF(std::exp(x2));
        }
    }
    const double r = std::exp(0.5 * (a + b));
    rBest = rF(r);
    return r;
}

// The nearest point of the segment [rA, rB] to rP, as the fraction along it.
double ff_segment_fraction(const PatchVec& rP, const PatchVec& rA, const PatchVec& rB) {
    const PatchVec d = ff_sub(rB, rA);
    const double dd = ff_dot(d, d);
    if (!(dd > 0.0))
        return 0.0;
    return std::clamp(ff_dot(ff_sub(rP, rA), d) / dd, 0.0, 1.0);
}

double ff_wrap(double A) {
    return A - 2.0 * kFfPi * std::round(A / (2.0 * kFfPi));
}

}  // namespace

double FeatureLaw::At(double T) const {
    const std::size_t n = mCoef.size();
    if (n == 0)
        return 0.0;
    std::vector<double> b(n);
    ff_law_basis(mPeriodic, mT0, mT1, n, T, b.data(), nullptr);
    double s = 0.0;
    for (std::size_t k = 0; k < n; ++k)
        s += mCoef[k] * b[k];
    return s;
}

double FeatureLaw::Slope(double T) const {
    const std::size_t n = mCoef.size();
    if (n == 0)
        return 0.0;
    std::vector<double> b(n), d(n);
    ff_law_basis(mPeriodic, mT0, mT1, n, T, b.data(), d.data());
    double s = 0.0;
    for (std::size_t k = 0; k < n; ++k)
        s += mCoef[k] * d[k];
    return s;
}

bool feature_law_fit(const std::vector<double>& rT, const std::vector<double>& rY, bool Periodic,
                     double T0, double T1, int Degree, FeatureLaw& rOut) {
    const std::size_t n = ff_law_size(Periodic, Degree);
    if (Degree < 0 || rT.size() != rY.size() || rT.size() <= n || !(T1 > T0))
        return false;
    std::vector<double> a(n * n, 0.0), b(n, 0.0), row(n);
    for (std::size_t i = 0; i < rT.size(); ++i) {
        ff_law_basis(Periodic, T0, T1, n, rT[i], row.data(), nullptr);
        for (std::size_t p = 0; p < n; ++p) {
            for (std::size_t q = 0; q < n; ++q)
                a[p * n + q] += row[p] * row[q];
            b[p] += row[p] * rY[i];
        }
    }
    if (!lsq_solve(a, b, n))
        return false;
    rOut.mPeriodic = Periodic;
    rOut.mT0 = T0;
    rOut.mT1 = T1;
    rOut.mCoef = std::move(b);
    return true;
}

PatchVec RollingBallFit::RailAt(double T) const {
    if (mRail.empty())
        return {0.0, 0.0, 0.0};
    const std::size_t n = mRail.size();
    if (mClosed) {
        T -= mLength * std::floor(T / mLength);
    } else if (T <= 0.0) {
        return ff_axpy(mRail[0], T, RailTangent(0.0));
    } else if (T >= mLength) {
        return ff_axpy(mRail[n - 1], T - mLength, RailTangent(mLength));
    }
    const std::size_t i = static_cast<std::size_t>(
        std::upper_bound(mRailT.begin(), mRailT.end(), T) - mRailT.begin());
    const std::size_t a = i == 0 ? 0 : i - 1;
    const std::size_t b = (a + 1) % n;
    const double t0 = mRailT[a], t1 = a + 1 < n ? mRailT[a + 1] : mLength;
    const double f = t1 > t0 ? (T - t0) / (t1 - t0) : 0.0;
    return ff_axpy(mRail[a], f, ff_sub(mRail[b], mRail[a]));
}

PatchVec RollingBallFit::RailTangent(double T) const {
    const std::size_t n = mRail.size();
    if (n < 2 || !(mLength > 0.0))
        return {0.0, 0.0, 0.0};
    // A central difference over a few segments: the rail's points are
    // averages, and one segment's direction is noisy.
    const double h = std::max(2.0 * mLength / static_cast<double>(n), 0.02 * mLength);
    double a = T - h, b = T + h;
    if (!mClosed) {
        a = std::clamp(a, 0.0, mLength);
        b = std::clamp(b, 0.0, mLength);
        if (!(b > a)) {
            a = T < 0.5 * mLength ? 0.0 : mLength - h;
            b = a + h;
        }
    }
    auto at = [&](double U) {
        if (mClosed)
            U -= mLength * std::floor(U / mLength);
        const std::size_t i = static_cast<std::size_t>(
            std::upper_bound(mRailT.begin(), mRailT.end(), U) - mRailT.begin());
        const std::size_t s = i == 0 ? 0 : std::min(i - 1, mClosed ? n - 1 : n - 2);
        const std::size_t e = (s + 1) % n;
        const double t0 = mRailT[s], t1 = s + 1 < n ? mRailT[s + 1] : mLength;
        const double f = t1 > t0 ? (U - t0) / (t1 - t0) : 0.0;
        return ff_axpy(mRail[s], f, ff_sub(mRail[e], mRail[s]));
    };
    return ff_unit(ff_sub(at(b), at(a)));
}

RollingBallFit rolling_ball_fit(const std::vector<PatchVec>& rPts,
                                const std::vector<PatchTri>& rTris, const ClosestPoint& rS1,
                                const ClosestPoint& rS2, double RMin, double RMax, double Tol) {
    RollingBallFit out;
    if (!(RMin > 0.0) || !(RMax > RMin) || !(Tol > 0.0))
        return out;
    std::vector<PatchVec> normals;
    std::vector<std::uint8_t> interior;
    ff_vertex_normals(rPts, rTris, normals, interior);
    std::vector<std::size_t> verts;
    for (std::size_t v = 0; v < rPts.size(); ++v)
        if (interior[v] && ff_norm(normals[v]) > 0.5)
            verts.push_back(v);
    if (verts.size() < 8)
        return out;
    const double inf = std::numeric_limits<double>::infinity();
    auto dist = [&](const ClosestPoint& rS, const PatchVec& rX) {
        PatchVec f{};
        return rS(rX, f) ? ff_norm(ff_sub(rX, f)) : inf;
    };
    auto misfit = [&](std::size_t V, int Side, double R) {
        const PatchVec c = ff_axpy(rPts[V], -Side * R, normals[V]);
        const double d1 = dist(rS1, c) - R, d2 = dist(rS2, c) - R;
        return d1 * d1 + d2 * d2;
    };
    // Which side the centres lie on, from a few vertices: a patch most of
    // whose vertices fit no ball is no blend, and is refused early.
    const std::vector<std::size_t> probe = ff_spread(verts, 32);
    int side = 0;
    std::size_t side_count = 0;
    double side_sum = inf;
    for (int s : {1, -1}) {
        std::size_t count = 0;
        double sum = 0.0;
        for (std::size_t v : probe) {
            double g = 0.0;
            ff_min_radius([&](double R) { return misfit(v, s, R); }, RMin, RMax, 24, 30, g);
            if (std::sqrt(0.5 * g) <= Tol)
                ++count;
            sum += g;
        }
        if (count > side_count || (count == side_count && sum < side_sum)) {
            side = s;
            side_count = count;
            side_sum = sum;
        }
    }
    if (side == 0 || 5 * side_count < 4 * probe.size())
        return out;
    // Every vertex's own radius.
    const std::vector<std::size_t> used = ff_spread(verts, 256);
    std::vector<std::size_t> kept;
    std::vector<double> radius;
    for (std::size_t v : used) {
        double g = 0.0;
        const double r =
            ff_min_radius([&](double R) { return misfit(v, side, R); }, RMin, RMax, 24, 40, g);
        if (std::sqrt(0.5 * g) <= Tol) {
            kept.push_back(v);
            radius.push_back(r);
        }
    }
    if (10 * kept.size() < 9 * used.size())
        return out;
    std::vector<double> sorted = radius;
    std::sort(sorted.begin(), sorted.end());
    const double r_med = sorted[sorted.size() / 2];
    // The rail. The boundary vertices near the first surface, chained along
    // the boundary edges, only order it: a segmentation boundary zigzags
    // between two rows of vertices. The rail itself is made of where the
    // vertices' balls touch the first surface, points of the contact curve,
    // averaged in bins along that order.
    const auto uses = ff_edge_uses(rTris);
    std::vector<std::uint8_t> on_rail(rPts.size(), 0);
    for (const auto& kv : uses)
        if (kv.second == 1)
            for (std::size_t v : {kv.first.first, kv.first.second})
                if (!on_rail[v] && dist(rS1, rPts[v]) <= 0.1 * r_med &&
                    dist(rS1, rPts[v]) < dist(rS2, rPts[v]))
                    on_rail[v] = 1;
    std::vector<std::vector<std::size_t>> adj(rPts.size());
    for (const auto& kv : uses)
        if (kv.second == 1 && on_rail[kv.first.first] && on_rail[kv.first.second]) {
            adj[kv.first.first].push_back(kv.first.second);
            adj[kv.first.second].push_back(kv.first.first);
        }
    // The component of the rail edges with the greatest length.
    std::vector<std::int64_t> comp(rPts.size(), -1);
    std::vector<double> comp_len;
    std::vector<std::size_t> comp_first;
    for (std::size_t s = 0; s < rPts.size(); ++s) {
        if (adj[s].empty() || comp[s] >= 0)
            continue;
        const std::int64_t id = static_cast<std::int64_t>(comp_len.size());
        comp_len.push_back(0.0);
        comp_first.push_back(s);
        std::vector<std::size_t> stack{s};
        comp[s] = id;
        while (!stack.empty()) {
            const std::size_t u = stack.back();
            stack.pop_back();
            for (std::size_t w : adj[u]) {
                if (u < w)
                    comp_len.back() += ff_norm(ff_sub(rPts[w], rPts[u]));
                if (comp[w] < 0) {
                    comp[w] = id;
                    stack.push_back(w);
                }
            }
        }
    }
    if (comp_len.empty())
        return out;
    const std::size_t best = static_cast<std::size_t>(
        std::max_element(comp_len.begin(), comp_len.end()) - comp_len.begin());
    bool closed = true;
    std::size_t start = comp_first[best];
    for (std::size_t v = 0; v < rPts.size(); ++v)
        if (comp[v] == static_cast<std::int64_t>(best) && adj[v].size() != 2) {
            closed = false;
            if (adj[v].size() == 1) {
                start = v;
                break;
            }
        }
    for (auto& a : adj)
        std::sort(a.begin(), a.end());
    std::vector<std::uint8_t> seen(rPts.size(), 0);
    std::vector<std::size_t> chain{start};
    seen[start] = 1;
    for (;;) {
        std::size_t next = rPts.size();
        for (std::size_t w : adj[chain.back()])
            if (!seen[w]) {
                next = w;
                break;
            }
        if (next == rPts.size())
            break;
        seen[next] = 1;
        chain.push_back(next);
    }
    closed = closed && chain.size() >= 4;
    if (chain.size() < 3)
        return out;
    // The nearest point of a polyline to a point, as its arc length.
    auto along = [&](const std::vector<PatchVec>& rLine, const std::vector<double>& rAt,
                     bool Closed, double Total, const PatchVec& rP) {
        const std::size_t nseg = Closed ? rLine.size() : rLine.size() - 1;
        double nearest = inf, at = 0.0;
        for (std::size_t s = 0; s < nseg; ++s) {
            const PatchVec& a = rLine[s];
            const PatchVec& b = rLine[(s + 1) % rLine.size()];
            const double f = ff_segment_fraction(rP, a, b);
            const double d = ff_norm(ff_sub(rP, ff_axpy(a, f, ff_sub(b, a))));
            if (d < nearest) {
                nearest = d;
                const double t1 = s + 1 < rLine.size() ? rAt[s + 1] : Total;
                at = rAt[s] + f * (t1 - rAt[s]);
            }
        }
        return at;
    };
    std::vector<PatchVec> coarse;
    std::vector<double> coarse_at;
    double coarse_len = 0.0;
    for (std::size_t i = 0; i < chain.size(); ++i) {
        if (i > 0)
            coarse_len += ff_norm(ff_sub(rPts[chain[i]], rPts[chain[i - 1]]));
        coarse.push_back(rPts[chain[i]]);
        coarse_at.push_back(coarse_len);
    }
    if (closed)
        coarse_len += ff_norm(ff_sub(rPts[chain.front()], rPts[chain.back()]));
    if (!(coarse_len > 0.0))
        return out;
    std::vector<PatchVec> centres(kept.size()), feet(kept.size());
    std::vector<double> order(kept.size());
    for (std::size_t i = 0; i < kept.size(); ++i) {
        centres[i] = ff_axpy(rPts[kept[i]], -side * radius[i], normals[kept[i]]);
        if (!rS1(centres[i], feet[i]))
            return out;
        order[i] = along(coarse, coarse_at, closed, coarse_len, feet[i]);
    }
    // Half as many bins as the chain has vertices smooths a two-row zigzag.
    const std::size_t nbin = std::max<std::size_t>(4, chain.size() / 2);
    std::vector<PatchVec> bin_sum(nbin, PatchVec{0.0, 0.0, 0.0});
    std::vector<std::size_t> bin_count(nbin, 0);
    for (std::size_t i = 0; i < kept.size(); ++i) {
        const std::size_t b = std::min(
            nbin - 1, static_cast<std::size_t>(order[i] / coarse_len * static_cast<double>(nbin)));
        bin_sum[b] = ff_axpy(bin_sum[b], 1.0, feet[i]);
        ++bin_count[b];
    }
    out.mClosed = closed;
    double len = 0.0;
    for (std::size_t b = 0; b < nbin; ++b) {
        if (bin_count[b] == 0)
            continue;
        const double w = 1.0 / static_cast<double>(bin_count[b]);
        const PatchVec p{bin_sum[b][0] * w, bin_sum[b][1] * w, bin_sum[b][2] * w};
        if (!out.mRail.empty())
            len += ff_norm(ff_sub(p, out.mRail.back()));
        out.mRail.push_back(p);
        out.mRailT.push_back(len);
    }
    if (out.mRail.size() < 3)
        return out;
    if (closed)
        len += ff_norm(ff_sub(out.mRail.front(), out.mRail.back()));
    out.mLength = len;
    if (!(len > 0.0))
        return out;
    // Each vertex's parameter: where its ball touches the rail, a point in
    // the vertex's own cross-section (the nearest rail point to the vertex
    // itself drifts along a rail whose radius varies).
    for (std::size_t i = 0; i < kept.size(); ++i) {
        out.mT.push_back(along(out.mRail, out.mRailT, closed, len, feet[i]));
        out.mCentres.push_back(centres[i]);
    }
    // The radius law: the lowest degree whose residual stays within a
    // quarter of the tolerance.
    const int max_degree = out.mClosed ? 6 : 8;
    bool found = false;
    for (int d = 0; d <= max_degree && !found; ++d) {
        FeatureLaw law;
        if (!feature_law_fit(out.mT, radius, out.mClosed, 0.0, len, d, law))
            break;
        double sum = 0.0, worst = 0.0;
        for (std::size_t i = 0; i < radius.size(); ++i) {
            const double e = std::fabs(radius[i] - law.At(out.mT[i]));
            sum += e * e;
            worst = std::max(worst, e);
        }
        const double rms = std::sqrt(sum / static_cast<double>(radius.size()));
        if (rms <= 0.25 * Tol) {
            found = true;
            out.mRadius = std::move(law);
            out.mRms = rms;
            out.mMax = worst;
        }
    }
    if (!found)
        return out;
    double lo = inf, hi = -inf;
    for (int i = 0; i <= 64; ++i) {
        const double r = out.mRadius.At(len * i / 64.0);
        lo = std::min(lo, r);
        hi = std::max(hi, r);
    }
    if (!(lo > 0.0))
        return out;
    out.mVariable = hi - lo > Tol;
    out.mSide = side;
    out.mOk = true;
    return out;
}

RollingBallSection rolling_ball_section(const RollingBallFit& rFit, const ClosestPoint& rS1,
                                        const ClosestPoint& rS2, double T, const PatchVec& rTangent,
                                        double Scale) {
    RollingBallSection out;
    if (!rFit.mOk || rFit.mT.empty())
        return out;
    const double r = rFit.mRadius.At(T);
    if (!(r > 0.0))
        return out;
    // The guess: the vertex centre nearest along the rail, moved with it.
    std::size_t nearest = 0;
    double gap = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < rFit.mT.size(); ++i) {
        double d = std::fabs(rFit.mT[i] - T);
        if (rFit.mClosed) {
            d = std::fmod(d, rFit.mLength);
            d = std::min(d, rFit.mLength - d);
        }
        if (d < gap) {
            gap = d;
            nearest = i;
        }
    }
    const PatchVec on_rail = rFit.RailAt(T);
    PatchVec c =
        ff_axpy(rFit.mCentres[nearest], 1.0, ff_sub(on_rail, rFit.RailAt(rFit.mT[nearest])));
    // The section's plane goes through the rail point, so consecutive balls
    // advance along the spine as the rail does.
    const PatchVec t = ff_unit(rTangent);
    c = ff_axpy(c, -ff_dot(ff_sub(c, on_rail), t), t);
    PatchVec f1{}, f2{};
    for (int it = 0; it < 60; ++it) {
        if (!rS1(c, f1) || !rS2(c, f2))
            return out;
        const PatchVec g1 = ff_sub(c, f1), g2 = ff_sub(c, f2);
        const double d1 = ff_norm(g1), d2 = ff_norm(g2);
        if (!(d1 > 0.0) || !(d2 > 0.0))
            return out;
        const PatchVec u1{g1[0] / d1, g1[1] / d1, g1[2] / d1};
        const PatchVec u2{g2[0] / d2, g2[1] / d2, g2[2] / d2};
        const PatchVec w1 = ff_axpy(u1, -ff_dot(u1, t), t), w2 = ff_axpy(u2, -ff_dot(u2, t), t);
        const double a11 = ff_dot(u1, w1), a12 = ff_dot(u1, w2), a21 = ff_dot(u2, w1),
                     a22 = ff_dot(u2, w2);
        const double det = a11 * a22 - a12 * a21;
        if (!(std::fabs(det) > 1e-12))
            return out;
        const double e1 = d1 - r, e2 = d2 - r;
        const double al = (-e1 * a22 + e2 * a12) / det, be = (-e2 * a11 + e1 * a21) / det;
        const PatchVec step = ff_axpy(ff_axpy({0.0, 0.0, 0.0}, al, w1), be, w2);
        c = ff_axpy(c, 1.0, step);
        if (ff_norm(step) <= 1e-13 * Scale)
            break;
    }
    if (!rS1(c, f1) || !rS2(c, f2))
        return out;
    const double tol = 1e-6 * Scale;
    if (std::fabs(ff_norm(ff_sub(c, f1)) - r) > tol || std::fabs(ff_norm(ff_sub(c, f2)) - r) > tol)
        return out;
    // A foot far off the cross-section belongs to another ball.
    const double off = std::fabs(rFit.mRadius.Slope(T)) * r + 0.5 * r;
    if (std::fabs(ff_dot(ff_sub(f1, c), t)) > off || std::fabs(ff_dot(ff_sub(f2, c), t)) > off)
        return out;
    out.mOk = true;
    out.mCentre = c;
    out.mFoot1 = f1;
    out.mFoot2 = f2;
    out.mRadius = r;
    return out;
}

PatchVec HelicalFit::At(double S, double Theta) const {
    const double v = mProfile.At(S);
    double rho = mRhoOfZeta ? v : S, zeta = mRhoOfZeta ? S : v;
    if (mTaper != 0.0 || mQ != 0.0) {
        const double g = Advance(Theta);
        rho += mTaper * g;
        zeta += g;
    } else {
        zeta += mK * Theta;
    }
    PatchVec p = ff_axpy(mOrigin, rho * std::cos(Theta), mE1);
    p = ff_axpy(p, rho * std::sin(Theta), mE2);
    return ff_axpy(p, zeta, mAxis);
}

double HelicalFit::Pitch() const {
    return 2.0 * kFfPi * mK;
}

double HelicalFit::PitchAt(double Theta) const {
    return 2.0 * kFfPi * (mK + 2.0 * mQ * (Theta - mThetaRef));
}

double HelicalFit::Advance(double Theta) const {
    return mK * Theta + mQ * (Theta - mThetaRef) * (Theta - mThetaRef);
}

namespace {

// Axial height, distance from the axis and angle of every vertex about the
// axis (rO, rA) in the frame (rE1, rE2).
void ff_cylindrical(const std::vector<PatchVec>& rPts, const PatchVec& rO, const PatchVec& rA,
                    const PatchVec& rE1, const PatchVec& rE2, std::vector<double>& rZ,
                    std::vector<double>& rRho, std::vector<double>& rPhi) {
    const std::size_t n = rPts.size();
    rZ.resize(n);
    rRho.resize(n);
    rPhi.resize(n);
    for (std::size_t v = 0; v < n; ++v) {
        const PatchVec d = ff_sub(rPts[v], rO);
        rZ[v] = ff_dot(d, rA);
        const PatchVec radial = ff_axpy(d, -rZ[v], rA);
        rRho[v] = ff_norm(radial);
        rPhi[v] = std::atan2(ff_dot(radial, rE2), ff_dot(radial, rE1));
    }
}

// The angles unwrapped continuously over the edges, breadth-first from the
// lowest vertex; false when they do not unwrap (the patch winds round the
// axis) or the patch is in more than one piece.
bool ff_unwrap(const std::vector<std::vector<std::size_t>>& rAdj, const std::vector<double>& rPhi,
               std::vector<double>& rTheta) {
    const std::size_t n = rAdj.size();
    rTheta.assign(n, 0.0);
    std::vector<std::uint8_t> seen(n, 0);
    std::size_t start = 0;
    while (start < n && rAdj[start].empty())
        ++start;
    if (start == n)
        return false;
    rTheta[start] = rPhi[start];
    seen[start] = 1;
    std::queue<std::size_t> todo;
    todo.push(start);
    while (!todo.empty()) {
        const std::size_t u = todo.front();
        todo.pop();
        for (std::size_t v : rAdj[u]) {
            const double step = ff_wrap(rPhi[v] - rPhi[u]);
            if (!seen[v]) {
                seen[v] = 1;
                rTheta[v] = rTheta[u] + step;
                todo.push(v);
            } else if (std::fabs(rTheta[v] - rTheta[u] - step) > 1e-6) {
                return false;
            }
        }
    }
    for (std::size_t v = 0; v < n; ++v)
        if (!seen[v] && !rAdj[v].empty())
            return false;
    return true;
}

// The profile points (s, y) of the vertices: (rho, zeta), or (zeta, rho).
struct FfProfileFit {
    bool mOk = false;
    FeatureLaw mLaw;
    double mS0 = 0.0, mS1 = 0.0, mRms = 0.0, mMax = 0.0;
};

FfProfileFit ff_profile_fit(const std::vector<double>& rS, const std::vector<double>& rY,
                            int Degree) {
    FfProfileFit out;
    if (rS.empty())
        return out;
    out.mS0 = *std::min_element(rS.begin(), rS.end());
    out.mS1 = *std::max_element(rS.begin(), rS.end());
    if (!(out.mS1 > out.mS0) || !feature_law_fit(rS, rY, false, out.mS0, out.mS1, Degree, out.mLaw))
        return out;
    double sum = 0.0;
    for (std::size_t i = 0; i < rS.size(); ++i) {
        const double slope = out.mLaw.Slope(rS[i]);
        const double e = std::fabs(rY[i] - out.mLaw.At(rS[i])) / std::sqrt(1.0 + slope * slope);
        sum += e * e;
        out.mMax = std::max(out.mMax, e);
    }
    out.mRms = std::sqrt(sum / static_cast<double>(rS.size()));
    out.mOk = std::isfinite(out.mRms);
    return out;
}

}  // namespace

HelicalFit helical_fit(const std::vector<PatchVec>& rPts, const std::vector<PatchTri>& rTris,
                       const PatchVec* pOrigin, const PatchVec* pAxis, double Tol, bool General,
                       double TaperSeed) {
    HelicalFit out;
    const std::size_t n = rPts.size();
    if (n < 8 || rTris.empty() || !(Tol > 0.0) || (pOrigin == nullptr) != (pAxis == nullptr))
        return out;
    std::vector<PatchVec> normals;
    std::vector<std::uint8_t> interior;
    ff_vertex_normals(rPts, rTris, normals, interior);
    std::vector<double> area(n, 0.0);
    std::vector<std::vector<std::size_t>> adj(n);
    for (const PatchTri& t : rTris) {
        const double a =
            0.5 * ff_norm(ff_cross(ff_sub(rPts[t[1]], rPts[t[0]]), ff_sub(rPts[t[2]], rPts[t[0]])));
        for (int k = 0; k < 3; ++k) {
            area[t[static_cast<std::size_t>(k)]] += a / 3.0;
            adj[t[static_cast<std::size_t>(k)]].push_back(t[static_cast<std::size_t>((k + 1) % 3)]);
            adj[t[static_cast<std::size_t>((k + 1) % 3)]].push_back(t[static_cast<std::size_t>(k)]);
        }
    }
    for (auto& a : adj) {
        std::sort(a.begin(), a.end());
        a.erase(std::unique(a.begin(), a.end()), a.end());
    }
    std::vector<std::size_t> live;
    double wsum = 0.0;
    PatchVec centre{0.0, 0.0, 0.0};
    for (std::size_t v = 0; v < n; ++v)
        if (!adj[v].empty()) {
            live.push_back(v);
            wsum += area[v];
            centre = ff_axpy(centre, area[v], rPts[v]);
        }
    if (live.size() < 8 || !(wsum > 0.0))
        return out;
    centre = {centre[0] / wsum, centre[1] / wsum, centre[2] / wsum};
    double scale = 0.0;
    for (std::size_t v : live)
        scale = std::max(scale, ff_norm(ff_sub(rPts[v], centre)));
    if (!(scale > 0.0))
        return out;
    // The screw: its axis, a point on it and its advance per radian.
    PatchVec axis{}, origin{};
    double k = 0.0;
    if (pAxis == nullptr) {
        double A[9] = {0}, B[9] = {0}, C[9] = {0};
        for (std::size_t v : live) {
            const PatchVec x{(rPts[v][0] - centre[0]) / scale, (rPts[v][1] - centre[1]) / scale,
                             (rPts[v][2] - centre[2]) / scale};
            const PatchVec& nv = normals[v];
            const PatchVec m = ff_cross(x, nv);
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j) {
                    A[i * 3 + j] +=
                        area[v] * m[static_cast<std::size_t>(i)] * m[static_cast<std::size_t>(j)];
                    B[i * 3 + j] +=
                        area[v] * m[static_cast<std::size_t>(i)] * nv[static_cast<std::size_t>(j)];
                    C[i * 3 + j] +=
                        area[v] * nv[static_cast<std::size_t>(i)] * nv[static_cast<std::size_t>(j)];
                }
        }
        // C^-1 by its adjugate.
        const double det = C[0] * (C[4] * C[8] - C[5] * C[7]) - C[1] * (C[3] * C[8] - C[5] * C[6]) +
                           C[2] * (C[3] * C[7] - C[4] * C[6]);
        const double ctr = C[0] + C[4] + C[8];
        if (!(std::fabs(det) > 1e-12 * ctr * ctr * ctr))
            return out;
        const double Ci[9] = {(C[4] * C[8] - C[5] * C[7]) / det, (C[2] * C[7] - C[1] * C[8]) / det,
                              (C[1] * C[5] - C[2] * C[4]) / det, (C[5] * C[6] - C[3] * C[8]) / det,
                              (C[0] * C[8] - C[2] * C[6]) / det, (C[2] * C[3] - C[0] * C[5]) / det,
                              (C[3] * C[7] - C[4] * C[6]) / det, (C[1] * C[6] - C[0] * C[7]) / det,
                              (C[0] * C[4] - C[1] * C[3]) / det};
        // S = A - B C^-1 B^T.
        double BCi[9] = {0}, S[9] = {0};
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                for (int l = 0; l < 3; ++l)
                    BCi[i * 3 + j] += B[i * 3 + l] * Ci[l * 3 + j];
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) {
                double s = A[i * 3 + j];
                for (int l = 0; l < 3; ++l)
                    s -= BCi[i * 3 + l] * B[j * 3 + l];
                S[i * 3 + j] = s;
            }
        for (int i = 0; i < 3; ++i)
            for (int j = i + 1; j < 3; ++j)
                S[i * 3 + j] = S[j * 3 + i] = 0.5 * (S[i * 3 + j] + S[j * 3 + i]);
        double w[3], vec[9];
        lsq_sym3_eigen(S, w, vec);
        // One screw only: a plane, a cylinder or a sphere is left unchanged by
        // more than one, and its two smallest eigenvalues both vanish.
        if (!(w[1] > 100.0 * std::max(w[0], 1e-15 * std::fabs(w[2]))))
            return out;
        const PatchVec c{vec[0], vec[1], vec[2]};
        PatchVec cbar{0.0, 0.0, 0.0};
        for (int i = 0; i < 3; ++i) {
            double s = 0.0;
            for (int l = 0; l < 3; ++l)  // (C^-1 B^T c)_i = sum_l Ci_il (B^T c)_l
                s += Ci[i * 3 + l] *
                     (B[0 * 3 + l] * c[0] + B[1 * 3 + l] * c[1] + B[2 * 3 + l] * c[2]);
            cbar[static_cast<std::size_t>(i)] = -s;
        }
        axis = c;
        k = scale * ff_dot(c, cbar);
        origin = ff_axpy(centre, scale, ff_cross(c, cbar));
    } else {
        origin = *pOrigin;
        axis = ff_unit(*pAxis);
        if (!(ff_norm(axis) > 0.0))
            return out;
        double num = 0.0, den = 0.0;
        for (std::size_t v : live) {
            const double na = ff_dot(normals[v], axis);
            num += area[v] * na * ff_dot(normals[v], ff_cross(axis, ff_sub(rPts[v], origin)));
            den += area[v] * na * na;
        }
        if (!(den > 1e-12 * wsum))
            return out;  // no normal leans along the axis: a cylinder
        k = -num / den;
    }
    // A screw with no advance is a rotation: a surface of revolution.
    if (!(std::fabs(k) * 2.0 * kFfPi > 1e-3 * scale))
        return out;
    // Move the origin to the foot of the patch's centre on the axis.
    origin = ff_axpy(origin, ff_dot(ff_sub(centre, origin), axis), axis);
    PatchVec e1{}, e2{};
    lsq_frame(axis.data(), e1.data(), e2.data());
    std::vector<double> z, rho, phi, theta0;
    ff_cylindrical(rPts, origin, axis, e1, e2, z, rho, phi);
    double rho_max = 0.0;
    for (std::size_t v : live)
        rho_max = std::max(rho_max, rho[v]);
    for (std::size_t v : live)
        if (!(rho[v] > 1e-9 * rho_max))
            return out;  // the patch touches the axis
    if (!ff_unwrap(adj, phi, theta0))
        return out;
    // The geometry of the vertices for a perturbed screw; the angles follow
    // the reference unwrapping to their nearest branch.
    const std::vector<std::size_t> sub = ff_spread(live, 400);
    // The general screw's reference angle and angular span: the pitch's rate
    // of change is taken about the patch's mean angle, and its LM step is
    // scaled by the span so the parameters move on a comparable scale.
    double theta_ref = 0.0, theta_span = 0.0;
    if (General) {
        double lo = std::numeric_limits<double>::infinity(), hi = -lo, sum = 0.0;
        for (std::size_t v : live) {
            sum += theta0[v];
            lo = std::min(lo, theta0[v]);
            hi = std::max(hi, theta0[v]);
        }
        theta_ref = sum / static_cast<double>(live.size());
        theta_span = hi - lo;
        if (!(theta_span > 0.0))
            return out;
    }
    struct Screw {
        PatchVec mO, mA, mE1, mE2;
        double mK;
        double mTaper = 0.0, mQ = 0.0;
    };
    auto screw_of = [&](const std::vector<double>& rP) {
        Screw s;
        s.mA = ff_unit(ff_axpy(ff_axpy(axis, rP[0], e1), rP[1], e2));
        s.mO = ff_axpy(ff_axpy(origin, rP[2], e1), rP[3], e2);
        s.mE1 = ff_unit(ff_axpy(e1, -ff_dot(e1, s.mA), s.mA));
        s.mE2 = ff_cross(s.mA, s.mE1);
        s.mK = rP[4];
        if (rP.size() > 5) {
            s.mTaper = rP[5];
            s.mQ = rP[6];
        }
        return s;
    };
    auto profile_points = [&](const Screw& rS, const std::vector<std::size_t>& rWhich,
                              bool RhoOfZeta, std::vector<double>& rSv, std::vector<double>& rYv,
                              std::vector<double>* pTheta) {
        rSv.clear();
        rYv.clear();
        if (pTheta)
            pTheta->clear();
        for (std::size_t v : rWhich) {
            const PatchVec d = ff_sub(rPts[v], rS.mO);
            const double zz = ff_dot(d, rS.mA);
            const PatchVec radial = ff_axpy(d, -zz, rS.mA);
            const double rr = ff_norm(radial);
            const double ph = std::atan2(ff_dot(radial, rS.mE2), ff_dot(radial, rS.mE1));
            const double th = theta0[v] + ff_wrap(ph - theta0[v]);
            double zeta, rho_u = rr;
            if (General) {
                // Unscrewed: the advance and the taper's outward move undone.
                const double g = rS.mK * th + rS.mQ * (th - theta_ref) * (th - theta_ref);
                zeta = zz - g;
                rho_u = rr - rS.mTaper * g;
            } else {
                zeta = zz - rS.mK * th;
            }
            rSv.push_back(RhoOfZeta ? zeta : rho_u);
            rYv.push_back(RhoOfZeta ? rho_u : zeta);
            if (pTheta)
                pTheta->push_back(th);
        }
    };
    std::vector<double> params{0.0, 0.0, 0.0, 0.0, k};
    if (General) {
        // The taper starts where a neighbouring cone puts it, measured along
        // the axis this fit settled on.
        params.push_back(TaperSeed);
        params.push_back(0.0);
    }
    std::vector<double> sv, yv;
    // The profile's orientation for the refinement: the one a quartic
    // explains better at the start.
    bool rho_of_zeta = false;
    {
        const Screw s0 = screw_of(params);
        profile_points(s0, sub, false, sv, yv, nullptr);
        const FfProfileFit a = ff_profile_fit(sv, yv, 4);
        profile_points(s0, sub, true, sv, yv, nullptr);
        const FfProfileFit b = ff_profile_fit(sv, yv, 4);
        if (!a.mOk && !b.mOk)
            return out;
        rho_of_zeta = !a.mOk || (b.mOk && b.mRms < a.mRms);
    }
    const LsqResidual residual = [&](const std::vector<double>& rP, double* pOut) {
        const Screw s = screw_of(rP);
        std::vector<double> ss, yy;
        profile_points(s, sub, rho_of_zeta, ss, yy, nullptr);
        const FfProfileFit f = ff_profile_fit(ss, yy, 4);
        for (std::size_t i = 0; i < sub.size(); ++i) {
            if (!f.mOk) {
                pOut[i] = 1e3 * scale;
                continue;
            }
            const double slope = f.mLaw.Slope(ss[i]);
            pOut[i] = (yy[i] - f.mLaw.At(ss[i])) / std::sqrt(1.0 + slope * slope);
        }
    };
    if (General)
        lsq_levenberg_marquardt(
            params, {1.0, 1.0, scale, scale, scale, 1.0, scale / (theta_span * theta_span)},
            sub.size(), residual, 80);
    else
        lsq_levenberg_marquardt(params, {1.0, 1.0, scale, scale, scale}, sub.size(), residual, 60);
    const Screw fin = screw_of(params);
    // The profile of the lowest degree within the tolerance, over every vertex.
    std::vector<double> th;
    for (int d = 1; d <= 8 && !out.mOk; ++d)
        for (bool roz : {false, true}) {
            profile_points(fin, live, roz, sv, yv, &th);
            const FfProfileFit f = ff_profile_fit(sv, yv, d);
            if (!f.mOk || f.mMax > 0.4 * Tol)
                continue;
            out.mOk = true;
            out.mRhoOfZeta = roz;
            out.mProfile = f.mLaw;
            out.mS0 = f.mS0;
            out.mS1 = f.mS1;
            out.mRms = f.mRms;
            out.mMax = f.mMax;
            out.mStraight = d == 1;
            break;
        }
    if (!out.mOk)
        return out;
    out.mOk = false;
    out.mOrigin = fin.mO;
    out.mAxis = fin.mA;
    out.mE1 = fin.mE1;
    out.mE2 = fin.mE2;
    out.mK = fin.mK;
    if (General) {
        out.mTaper = fin.mTaper;
        out.mQ = fin.mQ;
        out.mThetaRef = theta_ref;
    }
    out.mThetaMin = *std::min_element(th.begin(), th.end());
    out.mThetaMax = *std::max_element(th.begin(), th.end());
    // Half a turn at least: less of a helical surface is as well something else.
    out.mOk =
        out.mThetaMax - out.mThetaMin >= kFfPi && std::fabs(out.mK) * 2.0 * kFfPi > 1e-3 * scale;
    return out;
}

// ---------------------------------------------------------------------------
// Chamfers and conic blends, and the boundary of a vertex blend (roadmap 6.3).
// ---------------------------------------------------------------------------

namespace {

// A polyline's point at arc length T: wrapped on a closed line, continued
// along the end segments past an open one.
PatchVec ff_line_at(const std::vector<PatchVec>& rLine, const std::vector<double>& rAt, bool Closed,
                    double Length, double T) {
    const std::size_t n = rLine.size();
    if (n == 0)
        return {0.0, 0.0, 0.0};
    if (n == 1)
        return rLine[0];
    if (Closed) {
        T -= Length * std::floor(T / Length);
    } else if (T <= 0.0) {
        const double d = rAt[1] - rAt[0];
        return d > 0.0 ? ff_axpy(rLine[0], T / d, ff_sub(rLine[1], rLine[0])) : rLine[0];
    } else if (T >= Length) {
        const double d = rAt[n - 1] - rAt[n - 2];
        return d > 0.0 ? ff_axpy(rLine[n - 1], (T - Length) / d, ff_sub(rLine[n - 1], rLine[n - 2]))
                       : rLine[n - 1];
    }
    const std::size_t i =
        static_cast<std::size_t>(std::upper_bound(rAt.begin(), rAt.end(), T) - rAt.begin());
    const std::size_t a = i == 0 ? 0 : std::min(i - 1, Closed ? n - 1 : n - 2);
    const std::size_t b = (a + 1) % n;
    const double t0 = rAt[a], t1 = a + 1 < n ? rAt[a + 1] : Length;
    const double f = t1 > t0 ? (T - t0) / (t1 - t0) : 0.0;
    return ff_axpy(rLine[a], f, ff_sub(rLine[b], rLine[a]));
}

// The polyline's unit tangent at T, a central difference over a few segments.
PatchVec ff_line_tangent(const std::vector<PatchVec>& rLine, const std::vector<double>& rAt,
                         bool Closed, double Length, double T) {
    const std::size_t n = rLine.size();
    if (n < 2 || !(Length > 0.0))
        return {0.0, 0.0, 0.0};
    const double h = std::max(1.5 * Length / static_cast<double>(n), 0.01 * Length);
    double a = T - h, b = T + h;
    if (!Closed) {
        a = std::clamp(a, 0.0, Length);
        b = std::clamp(b, 0.0, Length);
        if (!(b > a)) {
            a = T < 0.5 * Length ? 0.0 : Length - h;
            b = a + h;
        }
    }
    return ff_unit(ff_sub(ff_line_at(rLine, rAt, Closed, Length, b),
                          ff_line_at(rLine, rAt, Closed, Length, a)));
}

// The arc length of a polyline's point nearest to rP.
double ff_line_along(const std::vector<PatchVec>& rLine, const std::vector<double>& rAt,
                     bool Closed, double Length, const PatchVec& rP) {
    const std::size_t nseg = Closed ? rLine.size() : rLine.size() - 1;
    double nearest = std::numeric_limits<double>::infinity(), at = 0.0;
    for (std::size_t s = 0; s < nseg; ++s) {
        const PatchVec& a = rLine[s];
        const PatchVec& b = rLine[(s + 1) % rLine.size()];
        const double f = ff_segment_fraction(rP, a, b);
        const double d = ff_norm(ff_sub(rP, ff_axpy(a, f, ff_sub(b, a))));
        if (d < nearest) {
            nearest = d;
            const double t1 = s + 1 < rLine.size() ? rAt[s + 1] : Length;
            at = rAt[s] + f * (t1 - rAt[s]);
        }
    }
    return at;
}

// Where the plane through rO normal to rN cuts a polyline (an open one
// continued along its end segments as far again as it is long), the cut
// nearest rO; false when it misses.
bool ff_plane_cut(const std::vector<PatchVec>& rLine, bool Closed, double Length,
                  const PatchVec& rO, const PatchVec& rN, PatchVec& rCut) {
    const std::size_t n = rLine.size();
    if (n < 2)
        return false;
    std::vector<PatchVec> line = rLine;
    if (!Closed) {
        const PatchVec d0 = ff_unit(ff_sub(rLine[0], rLine[1]));
        const PatchVec d1 = ff_unit(ff_sub(rLine[n - 1], rLine[n - 2]));
        line.insert(line.begin(), ff_axpy(rLine[0], Length, d0));
        line.push_back(ff_axpy(rLine[n - 1], Length, d1));
    }
    const std::size_t nseg = Closed ? line.size() : line.size() - 1;
    double best = std::numeric_limits<double>::infinity();
    bool found = false;
    for (std::size_t s = 0; s < nseg; ++s) {
        const PatchVec& a = line[s];
        const PatchVec& b = line[(s + 1) % line.size()];
        const double da = ff_dot(ff_sub(a, rO), rN), db = ff_dot(ff_sub(b, rO), rN);
        if ((da > 0.0 && db > 0.0) || (da < 0.0 && db < 0.0) || da == db)
            continue;
        const PatchVec c = ff_axpy(a, da / (da - db), ff_sub(b, a));
        const double d = ff_norm(ff_sub(c, rO));
        if (d < best) {
            best = d;
            rCut = c;
            found = true;
        }
    }
    return found;
}

// The rational quadratic (rR1, rA, rR2; W) at U.
PatchVec ff_conic_at(const PatchVec& rR1, const PatchVec& rA, const PatchVec& rR2, double W,
                     double U) {
    const double b0 = (1.0 - U) * (1.0 - U), b1 = 2.0 * W * U * (1.0 - U), b2 = U * U;
    const double den = b0 + b1 + b2;
    PatchVec p = ff_axpy(ff_axpy({0.0, 0.0, 0.0}, b0 / den, rR1), b1 / den, rA);
    return ff_axpy(p, b2 / den, rR2);
}

// The distance from rP to the conic, its parameter searched over [-0.25,
// 1.25] (the blend's section, and a little past both ends).
double ff_conic_distance(const PatchVec& rR1, const PatchVec& rA, const PatchVec& rR2, double W,
                         const PatchVec& rP) {
    constexpr int kScan = 48;
    auto dist = [&](double U) { return ff_norm(ff_sub(ff_conic_at(rR1, rA, rR2, W, U), rP)); };
    int arg = 0;
    double low = std::numeric_limits<double>::infinity();
    for (int i = 0; i <= kScan; ++i) {
        const double d = dist(-0.25 + 1.5 * i / kScan);
        if (d < low) {
            low = d;
            arg = i;
        }
    }
    double a = -0.25 + 1.5 * std::max(arg - 1, 0) / kScan;
    double b = -0.25 + 1.5 * std::min(arg + 1, kScan) / kScan;
    const double g = 0.5 * (std::sqrt(5.0) - 1.0);
    double x1 = b - g * (b - a), x2 = a + g * (b - a);
    double f1 = dist(x1), f2 = dist(x2);
    for (int it = 0; it < 40; ++it) {
        if (f1 <= f2) {
            b = x2;
            x2 = x1;
            f2 = f1;
            x1 = b - g * (b - a);
            f1 = dist(x1);
        } else {
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = a + g * (b - a);
            f2 = dist(x2);
        }
    }
    return std::min(low, dist(0.5 * (a + b)));
}

// The boundary edges rKeep selects, chained: the component of greatest
// length, from one of its ends (or anywhere on a loop).
bool ff_boundary_chain(const std::vector<PatchVec>& rPts,
                       const std::map<std::pair<std::size_t, std::size_t>, int>& rUses,
                       const std::function<bool(std::size_t, std::size_t)>& rKeep,
                       std::vector<std::size_t>& rChain, bool& rClosed) {
    const std::size_t n = rPts.size();
    std::vector<std::vector<std::size_t>> adj(n);
    for (const auto& kv : rUses)
        if (kv.second == 1 && rKeep(kv.first.first, kv.first.second)) {
            adj[kv.first.first].push_back(kv.first.second);
            adj[kv.first.second].push_back(kv.first.first);
        }
    std::vector<std::int64_t> comp(n, -1);
    std::vector<double> comp_len;
    std::vector<std::size_t> comp_first;
    for (std::size_t s = 0; s < n; ++s) {
        if (adj[s].empty() || comp[s] >= 0)
            continue;
        const std::int64_t id = static_cast<std::int64_t>(comp_len.size());
        comp_len.push_back(0.0);
        comp_first.push_back(s);
        std::vector<std::size_t> stack{s};
        comp[s] = id;
        while (!stack.empty()) {
            const std::size_t u = stack.back();
            stack.pop_back();
            for (std::size_t w : adj[u]) {
                if (u < w)
                    comp_len.back() += ff_norm(ff_sub(rPts[w], rPts[u]));
                if (comp[w] < 0) {
                    comp[w] = id;
                    stack.push_back(w);
                }
            }
        }
    }
    if (comp_len.empty())
        return false;
    const std::size_t best = static_cast<std::size_t>(
        std::max_element(comp_len.begin(), comp_len.end()) - comp_len.begin());
    bool closed = true;
    std::size_t start = comp_first[best];
    for (std::size_t v = 0; v < n; ++v)
        if (comp[v] == static_cast<std::int64_t>(best) && adj[v].size() != 2) {
            closed = false;
            if (adj[v].size() == 1) {
                start = v;
                break;
            }
        }
    for (auto& a : adj)
        std::sort(a.begin(), a.end());
    std::vector<std::uint8_t> seen(n, 0);
    rChain.assign(1, start);
    seen[start] = 1;
    for (;;) {
        std::size_t next = n;
        for (std::size_t w : adj[rChain.back()])
            if (!seen[w]) {
                next = w;
                break;
            }
        if (next == n)
            break;
        seen[next] = 1;
        rChain.push_back(next);
    }
    rClosed = closed && rChain.size() >= 4;
    return rChain.size() >= 3;
}

}  // namespace

PatchVec SectionBlendFit::Normal(double T) const {
    const PatchVec r1 = ff_line_at(mRail1, mRail1T, mClosed, mLength, T);
    PatchVec tan = ff_line_tangent(mRail1, mRail1T, mClosed, mLength, T);
    // The edge the two surfaces would meet along, through the rail point:
    // the direction both are level in, n1 x n2.
    PatchVec n1{}, n2{};
    if (mNormal1 && mNormal2 && mNormal1(r1, n1) && mNormal2(r1, n2)) {
        PatchVec e = ff_cross(n1, n2);
        if (ff_norm(e) > 0.1) {
            e = ff_unit(e);
            tan = ff_dot(e, tan) < 0.0 ? PatchVec{-e[0], -e[1], -e[2]} : e;
        }
    }
    return tan;
}

double SectionBlendFit::ParameterOf(const PatchVec& rX) const {
    // From the nearest rail point, the parameter whose section plane holds
    // rX: a secant iteration on (rX - R1(T)) . N(T).
    double t = ff_line_along(mRail1, mRail1T, mClosed, mLength, rX);
    auto f = [&](double T) {
        return ff_dot(ff_sub(rX, ff_line_at(mRail1, mRail1T, mClosed, mLength, T)), Normal(T));
    };
    const double h = 1e-3 * mLength;
    double f0 = f(t);
    for (int it = 0; it < 30 && std::fabs(f0) > 1e-13 * mLength; ++it) {
        const double d = (f(t + h) - f0) / h;
        if (!(std::fabs(d) > 0.0))
            break;
        const double step = std::clamp(-f0 / d, -0.1 * mLength, 0.1 * mLength);
        t += step;
        f0 = f(t);
    }
    return t;
}

bool SectionBlendFit::Section(double T, PatchVec& rR1, PatchVec& rR2, PatchVec& rApex) const {
    if (mRail1.size() < 2 || mRail2.size() < 2)
        return false;
    rR1 = ff_line_at(mRail1, mRail1T, mClosed, mLength, T);
    const PatchVec tan = Normal(T);
    if (!(ff_norm(tan) > 0.0) || !ff_plane_cut(mRail2, mClosed, mLength2, rR1, tan, rR2))
        return false;
    const PatchVec chord = ff_sub(rR2, rR1);
    if (!(ff_norm(chord) > 0.0))
        return false;
    if (mChamfer) {
        rApex = ff_axpy(rR1, 0.5, chord);
        return true;
    }
    // The tangent lines of the two surfaces in the section's plane, each
    // heading towards the other end, and where they meet.
    PatchVec n1{}, n2{};
    if (!mNormal1 || !mNormal2 || !mNormal1(rR1, n1) || !mNormal2(rR2, n2))
        return false;
    PatchVec d1 = ff_unit(ff_cross(tan, n1)), d2 = ff_unit(ff_cross(tan, n2));
    if (!(ff_norm(d1) > 0.0) || !(ff_norm(d2) > 0.0))
        return false;
    if (ff_dot(d1, chord) < 0.0)
        d1 = {-d1[0], -d1[1], -d1[2]};
    if (ff_dot(d2, chord) > 0.0)
        d2 = {-d2[0], -d2[1], -d2[2]};
    // rR1 + a d1 = rR2 + b d2, by least squares.
    const double a11 = 1.0, a12 = -ff_dot(d1, d2), a22 = 1.0;
    const double r1 = ff_dot(d1, chord), r2 = -ff_dot(d2, chord);
    const double det = a11 * a22 - a12 * a12;
    if (!(det > 1e-10))
        return false;  // parallel tangents: no apex
    const double a = (r1 * a22 - a12 * r2) / det, b = (a11 * r2 - a12 * r1) / det;
    if (!(a > 0.0) || !(b > 0.0))
        return false;
    const PatchVec p1 = ff_axpy(rR1, a, d1), p2 = ff_axpy(rR2, b, d2);
    rApex = {0.5 * (p1[0] + p2[0]), 0.5 * (p1[1] + p2[1]), 0.5 * (p1[2] + p2[2])};
    return true;
}

bool SectionBlendFit::At(double T, double U, PatchVec& rP) const {
    PatchVec r1{}, r2{}, apex{};
    if (!Section(T, r1, r2, apex))
        return false;
    if (mChamfer) {
        rP = ff_axpy(r1, U, ff_sub(r2, r1));
        return true;
    }
    const double w = mWeight.At(T);
    if (!(w > 0.0))
        return false;
    rP = ff_conic_at(r1, apex, r2, w, U);
    return true;
}

SectionBlendFit section_blend_fit(const std::vector<PatchVec>& rPts,
                                  const std::vector<PatchTri>& rTris, const ClosestPoint& rS1,
                                  const NormalAt& rN1, const ClosestPoint& rS2, const NormalAt& rN2,
                                  double Tol,
                                  const std::function<int(std::size_t, std::size_t)>& rEdgeSide) {
    SectionBlendFit out;
    if (!(Tol > 0.0) || rTris.empty())
        return out;
    std::vector<PatchVec> normals;
    std::vector<std::uint8_t> interior;
    ff_vertex_normals(rPts, rTris, normals, interior);
    const double inf = std::numeric_limits<double>::infinity();
    auto dist = [&](const ClosestPoint& rS, const PatchVec& rX, PatchVec* pFoot) {
        PatchVec f{};
        if (!rS(rX, f))
            return inf;
        if (pFoot)
            *pFoot = f;
        return ff_norm(ff_sub(rX, f));
    };
    // The rails: the boundary edges along each surface, chained -- those the
    // caller says it lies across, else those whose ends both lie on it.
    const auto uses = ff_edge_uses(rTris);
    std::vector<std::uint8_t> on1(rPts.size(), 0), on2(rPts.size(), 0);
    if (!rEdgeSide)
        for (const auto& kv : uses)
            if (kv.second == 1)
                for (std::size_t v : {kv.first.first, kv.first.second}) {
                    const double d1 = dist(rS1, rPts[v], nullptr), d2 = dist(rS2, rPts[v], nullptr);
                    if (d1 <= Tol && d1 < d2)
                        on1[v] = 1;
                    else if (d2 <= Tol && d2 < d1)
                        on2[v] = 1;
                }
    auto along_side = [&](int Side) {
        return [&, Side](std::size_t A, std::size_t B) {
            if (rEdgeSide)
                return rEdgeSide(A, B) == Side;
            const std::vector<std::uint8_t>& on = Side == 1 ? on1 : on2;
            return on[A] != 0 && on[B] != 0;
        };
    };
    std::vector<std::size_t> chain1, chain2;
    bool closed1 = false, closed2 = false;
    if (!ff_boundary_chain(rPts, uses, along_side(1), chain1, closed1) ||
        !ff_boundary_chain(rPts, uses, along_side(2), chain2, closed2) || closed1 != closed2)
        return out;
    out.mClosed = closed1;
    auto rail = [&](const std::vector<std::size_t>& rChain, const ClosestPoint& rS,
                    std::vector<PatchVec>& rLine, std::vector<double>& rAt, double& rLen) {
        rLen = 0.0;
        for (std::size_t v : rChain) {
            PatchVec f = rPts[v];
            if (!std::isfinite(dist(rS, rPts[v], &f)))
                return false;
            if (!rLine.empty())
                rLen += ff_norm(ff_sub(f, rLine.back()));
            rLine.push_back(f);
            rAt.push_back(rLen);
        }
        if (out.mClosed)
            rLen += ff_norm(ff_sub(rLine.front(), rLine.back()));
        return rLen > 0.0;
    };
    if (!rail(chain1, rS1, out.mRail1, out.mRail1T, out.mLength) ||
        !rail(chain2, rS2, out.mRail2, out.mRail2T, out.mLength2))
        return out;
    out.mNormal1 = rN1;
    out.mNormal2 = rN2;
    // Every interior vertex in its section: the chord's distance, and the
    // conic weight its barycentric coordinates give.
    std::vector<std::size_t> verts;
    for (std::size_t v = 0; v < rPts.size(); ++v)
        if (interior[v] && ff_norm(normals[v]) > 0.5)
            verts.push_back(v);
    verts = ff_spread(verts, 256);
    if (verts.size() < 6)
        return out;
    std::vector<double> at, chord_dist, weight_t, weight;
    std::vector<std::array<PatchVec, 4>> sections;  // r1, r2, apex, vertex
    SectionBlendFit probe = out;
    probe.mChamfer = false;
    for (std::size_t v : verts) {
        const PatchVec& x = rPts[v];
        const double t = out.ParameterOf(x);
        PatchVec r1{}, r2{}, apex{};
        out.mChamfer = true;
        if (!out.Section(t, r1, r2, apex))
            continue;
        const double f = ff_segment_fraction(x, r1, r2);
        at.push_back(t);
        chord_dist.push_back(ff_norm(ff_sub(x, ff_axpy(r1, f, ff_sub(r2, r1)))));
        if (!probe.Section(t, r1, r2, apex))
            continue;
        // Barycentric coordinates in (r1, apex, r2), in that triangle's plane.
        const PatchVec e1 = ff_sub(apex, r1), e2 = ff_sub(r2, r1), q = ff_sub(x, r1);
        const double g11 = ff_dot(e1, e1), g12 = ff_dot(e1, e2), g22 = ff_dot(e2, e2);
        const double det = g11 * g22 - g12 * g12;
        if (!(det > 0.0))
            continue;
        const double q1 = ff_dot(q, e1), q2 = ff_dot(q, e2);
        const double b1 = (q1 * g22 - q2 * g12) / det, b2 = (q2 * g11 - q1 * g12) / det;
        const double b0 = 1.0 - b1 - b2;
        if (!(b0 > 0.0) || !(b2 > 0.0) || b1 < 0.0)
            continue;
        weight_t.push_back(t);
        weight.push_back(b1 / (2.0 * std::sqrt(b0 * b2)));
        sections.push_back({r1, r2, apex, x});
    }
    if (10 * at.size() < 9 * verts.size())
        return out;
    // A chamfer: every vertex on its section's chord.
    double worst = 0.0, sum = 0.0;
    for (double d : chord_dist) {
        worst = std::max(worst, d);
        sum += d * d;
    }
    if (worst <= 0.4 * Tol) {
        out.mChamfer = true;
        out.mRms = std::sqrt(sum / static_cast<double>(chord_dist.size()));
        out.mMax = worst;
        out.mOk = true;
        return out;
    }
    // A conic: the weight law of the lowest degree whose conics pass near
    // every vertex.
    out.mChamfer = false;
    if (10 * weight.size() < 9 * verts.size())
        return out;
    const int max_degree = out.mClosed ? 6 : 8;
    for (int d = 0; d <= max_degree; ++d) {
        FeatureLaw law;
        if (!feature_law_fit(weight_t, weight, out.mClosed, 0.0, out.mLength, d, law))
            break;
        double w_lo = inf;
        for (int i = 0; i <= 64; ++i)
            w_lo = std::min(w_lo, law.At(out.mLength * i / 64.0));
        if (!(w_lo > 0.0))
            continue;
        worst = 0.0;
        sum = 0.0;
        for (std::size_t i = 0; i < sections.size(); ++i) {
            const auto& sec = sections[i];
            const double e = ff_conic_distance(sec[0], sec[2], sec[1], law.At(weight_t[i]), sec[3]);
            worst = std::max(worst, e);
            sum += e * e;
        }
        if (worst <= 0.4 * Tol) {
            out.mWeight = std::move(law);
            out.mRms = std::sqrt(sum / static_cast<double>(sections.size()));
            out.mMax = worst;
            out.mOk = true;
            return out;
        }
    }
    return out;
}

VertexBlendBoundary vertex_blend_boundary(
    const std::vector<PatchTri>& rTris,
    const std::function<std::int64_t(std::size_t, std::size_t)>& rAcross) {
    VertexBlendBoundary out;
    // The directed boundary edges, as the triangles wind them: one loop when
    // each boundary vertex starts exactly one.
    std::map<std::pair<std::size_t, std::size_t>, int> uses = ff_edge_uses(rTris);
    std::map<std::size_t, std::size_t> next;
    for (const PatchTri& t : rTris)
        for (int k = 0; k < 3; ++k) {
            const std::size_t a = t[static_cast<std::size_t>(k)];
            const std::size_t b = t[static_cast<std::size_t>((k + 1) % 3)];
            if (uses[{std::min(a, b), std::max(a, b)}] != 1)
                continue;
            if (!next.emplace(a, b).second)
                return out;  // two boundary edges leave one vertex: not a disc
        }
    if (next.size() < 3)
        return out;
    const std::size_t start = next.begin()->first;
    std::size_t v = start;
    do {
        out.mLoop.push_back(v);
        auto it = next.find(v);
        if (it == next.end() || out.mLoop.size() > next.size())
            return out;
        out.mAcross.push_back(rAcross(v, it->second));
        v = it->second;
    } while (v != start);
    if (out.mLoop.size() != next.size())
        return out;  // more than one loop: a hole
    // The runs: rotate so the loop starts where the neighbour changes.
    const std::size_t n = out.mLoop.size();
    std::size_t first = n;
    for (std::size_t i = 0; i < n; ++i)
        if (out.mAcross[i] != out.mAcross[(i + n - 1) % n]) {
            first = i;
            break;
        }
    if (first == n)
        return out;  // one neighbour all round
    std::rotate(out.mLoop.begin(), out.mLoop.begin() + static_cast<std::ptrdiff_t>(first),
                out.mLoop.end());
    std::rotate(out.mAcross.begin(), out.mAcross.begin() + static_cast<std::ptrdiff_t>(first),
                out.mAcross.end());
    for (std::size_t i = 0; i < n; ++i) {
        if (i > 0 && out.mAcross[i] == out.mAcross[i - 1])
            continue;
        if (out.mAcross[i] < 0)
            return out;  // an edge with nothing across
        if (std::find(out.mNeighbours.begin(), out.mNeighbours.end(), out.mAcross[i]) !=
            out.mNeighbours.end())
            return out;  // a neighbour met twice: its run is not contiguous
        out.mNeighbours.push_back(out.mAcross[i]);
    }
    out.mOk = out.mNeighbours.size() >= 3;
    return out;
}

JunctionPoint junction_closest(const std::vector<ClosestPoint>& rSurfaces,
                               const std::vector<NormalAt>& rNormals, const PatchVec& rStart,
                               double Scale, int MaxIterations) {
    JunctionPoint out;
    const std::size_t n = rSurfaces.size();
    if (n == 0 || rNormals.size() != n)
        return out;
    constexpr double kMu = 1e-6;
    // The largest distance to the surfaces from @p rX, and their feet and
    // normals there; false when a foot or a normal is missing.
    std::vector<PatchVec> foot(n), normal(n);
    auto measure = [&](const PatchVec& rX, double& rWorst) {
        rWorst = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            if (!rSurfaces[i](rX, foot[i]))
                return false;
            const PatchVec d{rX[0] - foot[i][0], rX[1] - foot[i][1], rX[2] - foot[i][2]};
            const double len = ff_norm(d);
            rWorst = std::max(rWorst, len);
            if (rNormals[i] && rNormals[i](rX, normal[i]))
                continue;
            if (!(len > 0.0))
                return false;
            normal[i] = {d[0] / len, d[1] / len, d[2] / len};
        }
        return true;
    };
    PatchVec p = rStart;
    double worst = 0.0;
    if (!measure(p, worst))
        return out;
    int it = 0;
    for (; it < MaxIterations; ++it) {
        std::vector<double> a(9, 0.0);
        double b[3];
        for (int r = 0; r < 3; ++r) {
            a[static_cast<std::size_t>(r * 3 + r)] = kMu;
            b[r] = kMu * p[static_cast<std::size_t>(r)];
        }
        for (std::size_t i = 0; i < n; ++i) {
            const PatchVec& nn = normal[i];
            const double off = ff_dot(nn, foot[i]);
            for (std::size_t r = 0; r < 3; ++r) {
                for (std::size_t c = 0; c < 3; ++c)
                    a[r * 3 + c] += nn[r] * nn[c];
                b[r] += nn[r] * off;
            }
        }
        if (!chol_factor(a, 3))
            break;
        chol_solve(a, 3, b);
        PatchVec step{b[0] - p[0], b[1] - p[1], b[2] - p[2]};
        bool moved = false;
        for (int half = 0; half <= 10; ++half) {
            const PatchVec q{p[0] + step[0], p[1] + step[1], p[2] + step[2]};
            double w = 0.0;
            if (measure(q, w) && w <= worst) {
                moved = ff_norm(step) > 1e-13 * Scale;
                p = q;
                worst = w;
                break;
            }
            step = {0.5 * step[0], 0.5 * step[1], 0.5 * step[2]};
        }
        if (!moved)
            break;
    }
    // The feet and normals at the point kept (a rejected trial overwrote them).
    if (!measure(p, worst))
        return out;
    out.mOk = true;
    out.mPoint = p;
    out.mMaxDistance = worst;
    out.mIterations = it;
    return out;
}

}  // namespace detail
}  // namespace meshioplusplus

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

// fit_primitives. See operations/fit_primitives.hpp for the contract.

// System includes
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/fit_primitives.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/parallel.hpp"

// Project includes (private, not installed)
#include "../detail/cell_region_mask.hpp"
#include "../detail/crease_edges.hpp"
#include "../detail/primitive_grow.hpp"
#include "../detail/lsq_fit.hpp"
#include "../detail/primitive_geom.hpp"
#include "../detail/surface_rings.hpp"
#include "../detail/math_utils.hpp"
#include "meshioplusplus/detail/space_filling.hpp"
#include "meshioplusplus/region.hpp"

namespace meshioplusplus {
namespace {

constexpr const char* kFpPrefix = "meshio++: fit_primitives: ";
constexpr double kFpPi = 3.14159265358979323846;
constexpr std::size_t kFpParams = 10;

// One patch's samples, shifted so the points' mean is the origin.
struct FpPatch {
    std::vector<double> mPts;      // 3 per point
    std::vector<double> mNormals;  // 3 per face, unit (or zero)
    std::vector<double> mCentres;  // 3 per face
    double mShift[3] = {0.0, 0.0, 0.0};
    double mDiag = 0.0;
    std::size_t NumPoints() const { return mPts.size() / 3; }
    std::size_t NumFaces() const { return mCentres.size() / 3; }
};

// A candidate: parameters (shifted frame) and the distances of the points.
struct FpFit {
    bool mOk = false;
    double mParams[kFpParams] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double mRms = std::numeric_limits<double>::quiet_NaN();
    double mMax = std::numeric_limits<double>::quiet_NaN();
    double mFraction = 0.0;
};

void fp_normalize(double* pV) {
    const double len = std::sqrt(detail::MathUtils::Dot(pV, pV));
    if (len > 0.0)
        for (int k = 0; k < 3; ++k)
            pV[k] /= len;
}

// Sign an axis so its largest-magnitude component (the first on a tie) is positive.
void fp_canonical_axis(double* pV) {
    int big = 0;
    for (int k = 1; k < 3; ++k)
        if (std::fabs(pV[k]) > std::fabs(pV[big]))
            big = k;
    if (pV[big] < 0.0)
        for (int k = 0; k < 3; ++k)
            pV[k] = -pV[k];
}

// Covariance (sum of outer products) of 3-vectors about a mean.
void fp_scatter(const std::vector<double>& rV, const double* pMean, double* pOut) {
    for (int k = 0; k < 9; ++k)
        pOut[k] = 0.0;
    const std::size_t n = rV.size() / 3;
    for (std::size_t i = 0; i < n; ++i) {
        const double d[3] = {rV[i * 3] - pMean[0], rV[i * 3 + 1] - pMean[1],
                             rV[i * 3 + 2] - pMean[2]};
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b)
                pOut[a * 3 + b] += d[a] * d[b];
    }
}

void fp_mean(const std::vector<double>& rV, double* pOut) {
    pOut[0] = pOut[1] = pOut[2] = 0.0;
    const std::size_t n = rV.size() / 3;
    for (std::size_t i = 0; i < n; ++i)
        for (int k = 0; k < 3; ++k)
            pOut[k] += rV[i * 3 + k];
    for (int k = 0; k < 3; ++k)
        pOut[k] /= static_cast<double>(n);
}

// Residual statistics of a set of signed distances.
void fp_stats(const std::vector<double>& rDist, double Tol, FpFit& rFit) {
    double sum = 0.0, mx = 0.0;
    std::size_t inside = 0;
    for (double d : rDist) {
        const double a = std::fabs(d);
        sum += d * d;
        mx = a > mx ? a : mx;
        inside += a <= Tol ? 1 : 0;
    }
    const double n = static_cast<double>(rDist.size());
    rFit.mRms = std::sqrt(sum / n);
    rFit.mMax = mx;
    rFit.mFraction = static_cast<double>(inside) / n;
}

// The axis direction cylinders and cones share: the direction the face
// normals vary least along; and the axis point nearest every normal line.
bool fp_axis(const FpPatch& rP, double* pDir, double* pPoint) {
    const std::size_t nf = rP.NumFaces();
    if (nf < 2)
        return false;
    double mean[3];
    fp_mean(rP.mNormals, mean);
    double cov[9];
    fp_scatter(rP.mNormals, mean, cov);
    double w[3], v[9];
    detail::lsq_sym3_eigen(cov, w, v);
    for (int k = 0; k < 3; ++k)
        pDir[k] = v[k];
    std::vector<double> a(9, 0.0), b(3, 0.0);
    for (std::size_t f = 0; f < nf; ++f) {
        const double* n = &rP.mNormals[f * 3];
        const double* c = &rP.mCentres[f * 3];
        for (int r = 0; r < 3; ++r)
            for (int s = 0; s < 3; ++s) {
                const double m = (r == s ? 1.0 : 0.0) - n[r] * n[s];
                a[static_cast<std::size_t>(r * 3 + s)] += m;
                b[static_cast<std::size_t>(r)] += m * c[s];
            }
    }
    const double wt = static_cast<double>(nf);
    for (int r = 0; r < 3; ++r)
        for (int s = 0; s < 3; ++s)
            a[static_cast<std::size_t>(r * 3 + s)] += wt * pDir[r] * pDir[s];
    if (!detail::lsq_solve(a, b, 3))
        return false;
    for (int k = 0; k < 3; ++k)
        pPoint[k] = b[static_cast<std::size_t>(k)];
    return true;
}

// The point and direction of an axis parameterised as an offset and a tilt
// in the frame (e1, e2, e3) about a base point.
void fp_axis_from(const double* pBase, const double* pE1, const double* pE2, const double* pE3,
                  double U, double V, double W, double Ax, double Ay, double* pPoint,
                  double* pDir) {
    for (int k = 0; k < 3; ++k) {
        pPoint[k] = pBase[k] + U * pE1[k] + V * pE2[k] + W * pE3[k];
        pDir[k] = pE3[k] + Ax * pE1[k] + Ay * pE2[k];
    }
    fp_normalize(pDir);
}

FpFit fp_plane(const FpPatch& rP, double Tol) {
    FpFit fit;
    const std::size_t n = rP.NumPoints();
    if (n < 3)
        return fit;
    const double origin[3] = {0.0, 0.0, 0.0};
    double cov[9], w[3], v[9];
    fp_scatter(rP.mPts, origin, cov);
    detail::lsq_sym3_eigen(cov, w, v);
    double nrm[3] = {v[0], v[1], v[2]};
    double total[3];
    fp_mean(rP.mNormals, total);
    if (detail::MathUtils::Dot(total, nrm) < 0.0)
        for (int k = 0; k < 3; ++k)
            nrm[k] = -nrm[k];
    std::vector<double> dist(n);
    for (std::size_t i = 0; i < n; ++i)
        dist[i] = detail::MathUtils::Dot(&rP.mPts[i * 3], nrm);
    fp_stats(dist, Tol, fit);
    for (int k = 0; k < 3; ++k)
        fit.mParams[3 + k] = nrm[k];
    fit.mOk = true;
    return fit;
}

FpFit fp_sphere(const FpPatch& rP, double Tol, int MaxIt) {
    FpFit fit;
    const std::size_t n = rP.NumPoints();
    if (n < 4)
        return fit;
    std::vector<double> a(16, 0.0), b(4, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const double* p = &rP.mPts[i * 3];
        const double row[4] = {2.0 * p[0], 2.0 * p[1], 2.0 * p[2], 1.0};
        const double rhs = detail::MathUtils::Dot(p, p);
        for (int r = 0; r < 4; ++r) {
            for (int s = 0; s < 4; ++s)
                a[static_cast<std::size_t>(r * 4 + s)] += row[r] * row[s];
            b[static_cast<std::size_t>(r)] += row[r] * rhs;
        }
    }
    if (!detail::lsq_solve(a, b, 4))
        return fit;
    const double r2 = b[3] + b[0] * b[0] + b[1] * b[1] + b[2] * b[2];
    if (!(r2 > 0.0))
        return fit;
    std::vector<double> theta = {b[0], b[1], b[2], std::sqrt(r2)};
    const std::vector<double> typical = {rP.mDiag, rP.mDiag, rP.mDiag, rP.mDiag};
    auto resid = [&](const std::vector<double>& rT, double* pOut) {
        for (std::size_t i = 0; i < n; ++i) {
            const double* p = &rP.mPts[i * 3];
            const double d[3] = {p[0] - rT[0], p[1] - rT[1], p[2] - rT[2]};
            pOut[i] = std::sqrt(detail::MathUtils::Dot(d, d)) - rT[3];
        }
    };
    detail::lsq_levenberg_marquardt(theta, typical, n, resid, MaxIt);
    if (!(theta[3] > 0.0) || !(theta[3] < 1e3 * rP.mDiag))
        return fit;
    std::vector<double> dist(n);
    resid(theta, dist.data());
    fp_stats(dist, Tol, fit);
    fit.mParams[0] = theta[0];
    fit.mParams[1] = theta[1];
    fit.mParams[2] = theta[2];
    fit.mParams[6] = theta[3];
    fit.mOk = true;
    return fit;
}

// Distance of a point to the line through pPoint along unit pDir.
double fp_line_distance(const double* pP, const double* pPoint, const double* pDir) {
    const double v[3] = {pP[0] - pPoint[0], pP[1] - pPoint[1], pP[2] - pPoint[2]};
    const double h = detail::MathUtils::Dot(v, pDir);
    const double r[3] = {v[0] - h * pDir[0], v[1] - h * pDir[1], v[2] - h * pDir[2]};
    return std::sqrt(detail::MathUtils::Dot(r, r));
}

// The axis point nearest the origin's projection: the foot of the origin
// (the shifted centroid) on the line.
void fp_foot_of_origin(const double* pPoint, const double* pDir, double* pOut) {
    const double h = -detail::MathUtils::Dot(pPoint, pDir);
    for (int k = 0; k < 3; ++k)
        pOut[k] = pPoint[k] + h * pDir[k];
}

FpFit fp_cylinder(const FpPatch& rP, double Tol, int MaxIt) {
    FpFit fit;
    const std::size_t n = rP.NumPoints();
    double dir[3], point[3];
    if (n < 5 || !fp_axis(rP, dir, point))
        return fit;
    double e1[3], e2[3];
    detail::lsq_frame(dir, e1, e2);
    double r0 = 0.0;
    for (std::size_t i = 0; i < n; ++i)
        r0 += fp_line_distance(&rP.mPts[i * 3], point, dir);
    r0 /= static_cast<double>(n);
    std::vector<double> theta = {0.0, 0.0, 0.0, 0.0, r0};
    const std::vector<double> typical = {rP.mDiag, rP.mDiag, 1.0, 1.0, rP.mDiag};
    auto resid = [&](const std::vector<double>& rT, double* pOut) {
        double pt[3], d[3];
        fp_axis_from(point, e1, e2, dir, rT[0], rT[1], 0.0, rT[2], rT[3], pt, d);
        for (std::size_t i = 0; i < n; ++i)
            pOut[i] = fp_line_distance(&rP.mPts[i * 3], pt, d) - rT[4];
    };
    detail::lsq_levenberg_marquardt(theta, typical, n, resid, MaxIt);
    if (!(theta[4] > 0.0) || !(theta[4] < 1e3 * rP.mDiag))
        return fit;
    std::vector<double> dist(n);
    resid(theta, dist.data());
    fp_stats(dist, Tol, fit);
    double pt[3], d[3], foot[3];
    fp_axis_from(point, e1, e2, dir, theta[0], theta[1], 0.0, theta[2], theta[3], pt, d);
    fp_canonical_axis(d);
    fp_foot_of_origin(pt, d, foot);
    for (int k = 0; k < 3; ++k) {
        fit.mParams[k] = foot[k];
        fit.mParams[3 + k] = d[k];
    }
    fit.mParams[6] = theta[4];
    fit.mOk = true;
    return fit;
}

FpFit fp_cone(const FpPatch& rP, double Tol, int MaxIt) {
    FpFit fit;
    const std::size_t n = rP.NumPoints();
    double dir[3], point[3];
    if (n < 6 || !fp_axis(rP, dir, point))
        return fit;
    // Radius against height: rho = a + b h.
    double sh = 0, sr = 0, shh = 0, shr = 0;
    for (std::size_t i = 0; i < n; ++i) {
        const double* p = &rP.mPts[i * 3];
        const double v[3] = {p[0] - point[0], p[1] - point[1], p[2] - point[2]};
        const double h = detail::MathUtils::Dot(v, dir);
        const double rho = fp_line_distance(p, point, dir);
        sh += h;
        sr += rho;
        shh += h * h;
        shr += h * rho;
    }
    const double nn = static_cast<double>(n);
    const double den = nn * shh - sh * sh;
    if (!(den > 0.0))
        return fit;
    const double slope = (nn * shr - sh * sr) / den;
    const double icpt = (sr - slope * sh) / nn;
    if (!(std::fabs(slope) > 1e-6))
        return fit;
    double axis[3] = {dir[0], dir[1], dir[2]};
    if (slope < 0.0)
        for (int k = 0; k < 3; ++k)
            axis[k] = -axis[k];
    const double h0 = -icpt / slope;
    double apex[3];
    for (int k = 0; k < 3; ++k)
        apex[k] = point[k] + h0 * dir[k];
    double e1[3], e2[3];
    detail::lsq_frame(axis, e1, e2);
    std::vector<double> theta = {0.0, 0.0, 0.0, 0.0, 0.0, std::atan(std::fabs(slope))};
    const std::vector<double> typical = {rP.mDiag, rP.mDiag, rP.mDiag, 1.0, 1.0, 1.0};
    auto resid = [&](const std::vector<double>& rT, double* pOut) {
        double a[3], d[3];
        fp_axis_from(apex, e1, e2, axis, rT[0], rT[1], rT[2], rT[3], rT[4], a, d);
        const double ca = std::cos(rT[5]), sa = std::sin(rT[5]);
        for (std::size_t i = 0; i < n; ++i) {
            const double* p = &rP.mPts[i * 3];
            const double v[3] = {p[0] - a[0], p[1] - a[1], p[2] - a[2]};
            const double h = detail::MathUtils::Dot(v, d);
            const double r[3] = {v[0] - h * d[0], v[1] - h * d[1], v[2] - h * d[2]};
            pOut[i] = std::sqrt(detail::MathUtils::Dot(r, r)) * ca - h * sa;
        }
    };
    detail::lsq_levenberg_marquardt(theta, typical, n, resid, MaxIt);
    if (!(theta[5] > 1e-3) || !(theta[5] < 0.5 * kFpPi - 1e-3))
        return fit;
    std::vector<double> dist(n);
    resid(theta, dist.data());
    fp_stats(dist, Tol, fit);
    double a[3], d[3];
    fp_axis_from(apex, e1, e2, axis, theta[0], theta[1], theta[2], theta[3], theta[4], a, d);
    for (int k = 0; k < 3; ++k) {
        fit.mParams[k] = a[k];
        fit.mParams[3 + k] = d[k];
    }
    fit.mParams[8] = theta[5] * 180.0 / kFpPi;
    fit.mOk = true;
    return fit;
}

// Tube centres for a minor radius: centroids moved back along their normals.
// Returns the circle's centre, axis, major radius and a misfit score.
bool fp_torus_guess(const FpPatch& rP, double Minor, double* pCentre, double* pAxis, double& rMajor,
                    double& rScore) {
    const std::size_t nf = rP.NumFaces();
    std::vector<double> q(nf * 3);
    for (std::size_t f = 0; f < nf; ++f)
        for (int k = 0; k < 3; ++k)
            q[f * 3 + k] = rP.mCentres[f * 3 + k] - Minor * rP.mNormals[f * 3 + k];
    double qm[3], cov[9], w[3], v[9];
    fp_mean(q, qm);
    fp_scatter(q, qm, cov);
    detail::lsq_sym3_eigen(cov, w, v);
    const double axis[3] = {v[0], v[1], v[2]};
    double e1[3], e2[3];
    detail::lsq_frame(axis, e1, e2);
    std::vector<double> a(9, 0.0), b(3, 0.0);
    for (std::size_t f = 0; f < nf; ++f) {
        const double d[3] = {q[f * 3] - qm[0], q[f * 3 + 1] - qm[1], q[f * 3 + 2] - qm[2]};
        const double u = detail::MathUtils::Dot(d, e1), t = detail::MathUtils::Dot(d, e2);
        const double row[3] = {2.0 * u, 2.0 * t, 1.0};
        const double rhs = u * u + t * t;
        for (int r = 0; r < 3; ++r) {
            for (int s = 0; s < 3; ++s)
                a[static_cast<std::size_t>(r * 3 + s)] += row[r] * row[s];
            b[static_cast<std::size_t>(r)] += row[r] * rhs;
        }
    }
    if (!detail::lsq_solve(a, b, 3))
        return false;
    const double r2 = b[2] + b[0] * b[0] + b[1] * b[1];
    if (!(r2 > 0.0))
        return false;
    rMajor = std::sqrt(r2);
    double score = 0.0;
    for (std::size_t f = 0; f < nf; ++f) {
        const double d[3] = {q[f * 3] - qm[0], q[f * 3 + 1] - qm[1], q[f * 3 + 2] - qm[2]};
        const double u = detail::MathUtils::Dot(d, e1) - b[0],
                     t = detail::MathUtils::Dot(d, e2) - b[1];
        const double radial = std::sqrt(u * u + t * t) - rMajor;
        const double off = detail::MathUtils::Dot(d, axis);
        score += radial * radial + off * off;
    }
    rScore = score;
    for (int k = 0; k < 3; ++k) {
        pCentre[k] = qm[k] + b[0] * e1[k] + b[1] * e2[k];
        pAxis[k] = axis[k];
    }
    return true;
}

FpFit fp_torus(const FpPatch& rP, double Tol, int MaxIt) {
    FpFit fit;
    const std::size_t n = rP.NumPoints();
    if (n < 7 || rP.NumFaces() < 3)
        return fit;
    bool found = false;
    double best = 0.0, centre[3] = {0, 0, 0}, axis[3] = {0, 0, 1}, major = 0.0, minor = 0.0;
    for (int sign = 0; sign < 2; ++sign)
        for (int j = 0; j < 48; ++j) {
            const double rho = rP.mDiag * std::pow(10.0, -3.0 + 3.0 * j / 47.0);
            const double s = sign == 0 ? rho : -rho;
            double c[3], ax[3], R, score;
            if (!fp_torus_guess(rP, s, c, ax, R, score))
                continue;
            if (!found || score < best) {
                found = true;
                best = score;
                for (int k = 0; k < 3; ++k) {
                    centre[k] = c[k];
                    axis[k] = ax[k];
                }
                major = R;
                minor = rho;
            }
        }
    if (!found)
        return fit;
    double e1[3], e2[3];
    detail::lsq_frame(axis, e1, e2);
    std::vector<double> theta = {0.0, 0.0, 0.0, 0.0, 0.0, major, minor};
    const std::vector<double> typical = {rP.mDiag, rP.mDiag, rP.mDiag, 1.0,
                                         1.0,      rP.mDiag, rP.mDiag};
    auto resid = [&](const std::vector<double>& rT, double* pOut) {
        double c[3], d[3];
        fp_axis_from(centre, e1, e2, axis, rT[0], rT[1], rT[2], rT[3], rT[4], c, d);
        for (std::size_t i = 0; i < n; ++i) {
            const double* p = &rP.mPts[i * 3];
            const double v[3] = {p[0] - c[0], p[1] - c[1], p[2] - c[2]};
            const double h = detail::MathUtils::Dot(v, d);
            const double r[3] = {v[0] - h * d[0], v[1] - h * d[1], v[2] - h * d[2]};
            const double radial = std::sqrt(detail::MathUtils::Dot(r, r)) - rT[5];
            pOut[i] = std::sqrt(radial * radial + h * h) - rT[6];
        }
    };
    detail::lsq_levenberg_marquardt(theta, typical, n, resid, MaxIt);
    if (!(theta[5] > 0.0) || !(theta[6] > 0.0) || !(theta[5] < 1e3 * rP.mDiag) ||
        !(theta[6] < 1e3 * rP.mDiag))
        return fit;
    std::vector<double> dist(n);
    resid(theta, dist.data());
    fp_stats(dist, Tol, fit);
    double c[3], d[3];
    fp_axis_from(centre, e1, e2, axis, theta[0], theta[1], theta[2], theta[3], theta[4], c, d);
    fp_canonical_axis(d);
    for (int k = 0; k < 3; ++k) {
        fit.mParams[k] = c[k];
        fit.mParams[3 + k] = d[k];
    }
    fit.mParams[6] = theta[5];
    fit.mParams[7] = theta[6];
    fit.mOk = true;
    return fit;
}

struct FpOutcome {
    std::int64_t mType = 0;
    double mParams[kFpParams];
    double mRms = std::numeric_limits<double>::quiet_NaN();
    double mMax = std::numeric_limits<double>::quiet_NaN();
    double mFraction = 0.0;
    double mTol = 0.0;
    detail::PrimitiveFormMetrics mForm{};
    // The points the fit used, and their RMS distance (roadmap 7.1.2).
    double mSigma = std::numeric_limits<double>::quiet_NaN();
    std::int64_t mNumUsed = 0;
    std::int64_t mNumRejected = 0;
    // The best attempt when no type was accepted (world frame): what sigma
    // trimming measures residuals against.
    std::int64_t mTryType = 0;
    double mTryParams[kFpParams] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
};

// How far the facets' centres lie from a fitted surface (the parameters in
// the patch's shifted frame), as an RMS. It tells apart the surfaces two rims
// of mesh points lie on alike: any two coaxial circles lie on a sphere, but
// the sphere bulges between the rims, away from the chords of a bore.
double fp_centre_misfit(std::int64_t Type, const double* pParams, const FpPatch& rP) {
    const std::size_t nf = rP.NumFaces();
    double sum = 0.0;
    for (std::size_t f = 0; f < nf; ++f) {
        const double d = detail::primitive_point_distance(Type, pParams, &rP.mCentres[f * 3]);
        sum += d * d;
    }
    return nf > 0 ? std::sqrt(sum / static_cast<double>(nf)) : 0.0;
}

FpOutcome fp_fit_patch(const FpPatch& rP, const FitPrimitivesOptions& rOptions) {
    FpOutcome out;
    const double nan = std::numeric_limits<double>::quiet_NaN();
    for (std::size_t k = 0; k < kFpParams; ++k)
        out.mParams[k] = nan;
    const double tol = rOptions.mRelative ? rOptions.mTolerance * rP.mDiag : rOptions.mTolerance;
    out.mTol = tol;
    if (rP.NumPoints() == 0 || !(rP.mDiag > 0.0))
        return out;
    const int max_it = static_cast<int>(std::min<std::int64_t>(rOptions.mMaxIterations, 100000));
    bool have_best = false;
    for (int t = 1; t <= 5; ++t) {
        if (!(rOptions.mTypes & (1u << t)))
            continue;
        FpFit fit;
        switch (t) {
            case 1:
                fit = fp_plane(rP, tol);
                break;
            case 2:
                fit = fp_sphere(rP, tol, max_it);
                break;
            case 3:
                fit = fp_cylinder(rP, tol, max_it);
                break;
            case 4:
                fit = fp_cone(rP, tol, max_it);
                break;
            default:
                fit = fp_torus(rP, tol, max_it);
                break;
        }
        if (!fit.mOk || !std::isfinite(fit.mRms))
            continue;
        if (fit.mRms <= tol && fit.mFraction >= rOptions.mMinInlierFraction) {
            // A sphere through two rims of points may as well be a cylinder,
            // cone or torus: a cylinder or cone the facets' centres lie nearer
            // to wins, and a torus (with a parameter to spare) when they lie
            // under half as far from it.
            int best_t = t;
            if (t == 2) {
                double best_misfit = fp_centre_misfit(2, fit.mParams, rP);
                for (int u = 3; u <= 5; ++u) {
                    if (!(rOptions.mTypes & (1u << u)))
                        continue;
                    const FpFit other = u == 3   ? fp_cylinder(rP, tol, max_it)
                                        : u == 4 ? fp_cone(rP, tol, max_it)
                                                 : fp_torus(rP, tol, max_it);
                    if (!other.mOk || !std::isfinite(other.mRms) || !(other.mRms <= tol) ||
                        other.mFraction < rOptions.mMinInlierFraction)
                        continue;
                    const double misfit = fp_centre_misfit(u, other.mParams, rP);
                    if (misfit < (u == 5 ? 0.5 : 1.0) * best_misfit) {
                        best_misfit = misfit;
                        best_t = u;
                        fit = other;
                    }
                }
            }
            out.mType = best_t;
            for (std::size_t k = 0; k < kFpParams; ++k)
                out.mParams[k] = fit.mParams[k];
            for (int k = 0; k < 3; ++k)
                out.mParams[k] += rP.mShift[k];  // back from the centred frame
            out.mRms = fit.mRms;
            out.mMax = fit.mMax;
            out.mFraction = fit.mFraction;
            return out;
        }
        if (!have_best || fit.mRms < out.mRms) {
            have_best = true;
            out.mRms = fit.mRms;
            out.mMax = fit.mMax;
            out.mFraction = fit.mFraction;
            out.mTryType = t;
            for (std::size_t k = 0; k < kFpParams; ++k)
                out.mTryParams[k] = fit.mParams[k];
            for (int k = 0; k < 3; ++k)
                out.mTryParams[k] += rP.mShift[k];
        }
    }
    return out;
}

// A copy of the patch with only the points @p rUsed marks, centred on their
// own mean as the fitters expect (its faces and diagonal kept).
FpPatch fp_subset(const FpPatch& rP, const std::vector<std::uint8_t>& rUsed) {
    FpPatch sub;
    sub.mDiag = rP.mDiag;
    for (std::size_t i = 0; i < rP.NumPoints(); ++i)
        if (rUsed[i])
            for (std::size_t k = 0; k < 3; ++k)
                sub.mPts.push_back(rP.mPts[i * 3 + k]);
    double mean[3] = {0.0, 0.0, 0.0};
    if (!sub.mPts.empty())
        fp_mean(sub.mPts, mean);
    for (std::size_t i = 0; i < sub.NumPoints(); ++i)
        for (std::size_t k = 0; k < 3; ++k)
            sub.mPts[i * 3 + k] -= mean[k];
    sub.mNormals = rP.mNormals;
    sub.mCentres = rP.mCentres;
    for (std::size_t f = 0; f < sub.NumFaces(); ++f)
        for (std::size_t k = 0; k < 3; ++k)
            sub.mCentres[f * 3 + k] -= mean[k];
    for (int k = 0; k < 3; ++k)
        sub.mShift[k] = rP.mShift[k] + mean[k];
    return sub;
}

// The points in Hilbert order over their bounding box (ties: point order).
std::vector<std::int64_t> fp_hilbert_order(const FpPatch& rP) {
    const std::size_t n = rP.NumPoints();
    double lo[3], hi[3];
    for (int k = 0; k < 3; ++k) {
        lo[k] = std::numeric_limits<double>::infinity();
        hi[k] = -lo[k];
    }
    for (std::size_t i = 0; i < n; ++i)
        for (int k = 0; k < 3; ++k) {
            lo[k] = std::min(lo[k], rP.mPts[i * 3 + static_cast<std::size_t>(k)]);
            hi[k] = std::max(hi[k], rP.mPts[i * 3 + static_cast<std::size_t>(k)]);
        }
    double ext = 0.0;
    for (int k = 0; k < 3; ++k)
        ext = std::max(ext, hi[k] - lo[k]);
    const double scale = ext > 0.0 ? static_cast<double>((1u << detail::sfc_bits) - 1) / ext : 0.0;
    std::vector<std::uint64_t> keys(n);
    for (std::size_t i = 0; i < n; ++i) {
        std::uint32_t q[3];
        for (int k = 0; k < 3; ++k)
            q[k] = static_cast<std::uint32_t>(
                (rP.mPts[i * 3 + static_cast<std::size_t>(k)] - lo[k]) * scale);
        keys[i] = detail::sfc_hilbert_key(q, detail::sfc_bits);
    }
    return detail::sfc_stable_argsort(keys);
}

// A least median of squares consensus (Rousseeuw 1984) for a plane (3 points)
// or a sphere (4), over a fixed schedule of minimal subsets drawn by stride
// from the Hilbert-ordered points, no random sampling; the inliers of the
// best subset: points within 2.5 robust standard deviations of it.
std::vector<std::uint8_t> fp_lmeds_inliers(const FpPatch& rP, std::int64_t Type) {
    const std::size_t n = rP.NumPoints();
    const std::size_t m = Type == 1 ? 3 : 4;
    std::vector<std::uint8_t> used(n, 1);
    if (n < 2 * m)
        return used;
    const std::vector<std::int64_t> order = fp_hilbert_order(rP);
    const std::size_t trials = std::min<std::size_t>(500, n);
    double best_med = std::numeric_limits<double>::infinity();
    double best[kFpParams] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<double> r2(n);
    for (std::size_t j = 0; j < trials; ++j) {
        std::size_t idx[4];
        for (std::size_t k = 0; k < m; ++k)
            idx[k] = static_cast<std::size_t>(order[(j * 7919 + k * (n / m) + k * j) % n]);
        double prm[kFpParams] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        const double* a = &rP.mPts[idx[0] * 3];
        const double* b = &rP.mPts[idx[1] * 3];
        const double* c = &rP.mPts[idx[2] * 3];
        if (Type == 1) {
            const double u[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
            const double v[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
            double nrm[3] = {u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2],
                             u[0] * v[1] - u[1] * v[0]};
            const double len = std::sqrt(detail::MathUtils::Dot(nrm, nrm));
            if (!(len > 0.0))
                continue;
            for (int k = 0; k < 3; ++k) {
                prm[k] = a[k];
                prm[3 + k] = nrm[k] / len;
            }
        } else {
            std::vector<double> mat(16, 0.0), rhs(4, 0.0);
            for (std::size_t r = 0; r < 4; ++r) {
                const double* p = &rP.mPts[idx[r] * 3];
                mat[r * 4 + 0] = 2.0 * p[0];
                mat[r * 4 + 1] = 2.0 * p[1];
                mat[r * 4 + 2] = 2.0 * p[2];
                mat[r * 4 + 3] = 1.0;
                rhs[r] = detail::MathUtils::Dot(p, p);
            }
            if (!detail::lsq_solve(mat, rhs, 4))
                continue;
            const double rad2 = rhs[3] + rhs[0] * rhs[0] + rhs[1] * rhs[1] + rhs[2] * rhs[2];
            if (!(rad2 > 0.0))
                continue;
            prm[0] = rhs[0];
            prm[1] = rhs[1];
            prm[2] = rhs[2];
            prm[6] = std::sqrt(rad2);
        }
        for (std::size_t i = 0; i < n; ++i) {
            const double d = detail::primitive_signed_distance(Type, prm, &rP.mPts[i * 3]);
            r2[i] = d * d;
        }
        std::vector<double> sorted = r2;
        std::nth_element(sorted.begin(), sorted.begin() + static_cast<std::ptrdiff_t>(n / 2),
                         sorted.end());
        const double med = sorted[n / 2];
        if (med < best_med) {
            best_med = med;
            std::copy(prm, prm + kFpParams, best);
        }
    }
    if (!std::isfinite(best_med))
        return used;
    // Rousseeuw's robust scale, with the small-sample factor.
    const double scale = 1.4826 * (1.0 + 5.0 / static_cast<double>(n - m)) * std::sqrt(best_med);
    // An exact consensus has no spread: keep what lies on it to rounding.
    const double cut = std::max(2.5 * scale, 1e-12 * rP.mDiag);
    std::size_t kept = 0;
    for (std::size_t i = 0; i < n; ++i) {
        const double d = std::fabs(detail::primitive_signed_distance(Type, best, &rP.mPts[i * 3]));
        used[i] = d <= cut ? 1 : 0;
        kept += used[i];
    }
    if (kept < 2 * m)
        std::fill(used.begin(), used.end(), std::uint8_t{1});
    return used;
}

// The sigma-trimmed fit of a patch (roadmap 7.1.2): refit on the points within
// K times the RMS distance of the points kept, until that set stops changing.
FpOutcome fp_fit_patch_trimmed(const FpPatch& rP, const FitPrimitivesOptions& rOptions, double K,
                               std::vector<std::uint8_t>& rUsed) {
    const std::size_t n = rP.NumPoints();
    std::vector<std::uint8_t>& used = rUsed;
    used.assign(n, 1);
    if (rOptions.mTypes == (1u << 1) || rOptions.mTypes == (1u << 2))
        used = fp_lmeds_inliers(rP, rOptions.mTypes == (1u << 1) ? 1 : 2);
    // The tolerance stays the whole patch's.
    FitPrimitivesOptions o = rOptions;
    if (o.mRelative) {
        o.mTolerance = rOptions.mTolerance * rP.mDiag;
        o.mRelative = false;
    }
    FpOutcome out;
    std::vector<double> dist(n);
    for (int round = 0; round < 20; ++round) {
        out = fp_fit_patch(fp_subset(rP, used), o);
        const std::int64_t type = out.mType ? out.mType : out.mTryType;
        const double* prm = out.mType ? out.mParams : out.mTryParams;
        if (type == 0)
            break;
        double sum = 0.0;
        std::size_t count = 0;
        for (std::size_t i = 0; i < n; ++i) {
            double w[3];
            for (std::size_t k = 0; k < 3; ++k)
                w[k] = rP.mPts[i * 3 + k] + rP.mShift[k];
            dist[i] = detail::primitive_signed_distance(type, prm, w);
            if (used[i]) {
                sum += dist[i] * dist[i];
                ++count;
            }
        }
        const double sigma = count ? std::sqrt(sum / static_cast<double>(count)) : 0.0;
        if (!(sigma > 0.0))
            break;
        std::vector<std::uint8_t> next(n);
        std::size_t kept = 0;
        for (std::size_t i = 0; i < n; ++i) {
            next[i] = std::fabs(dist[i]) <= K * sigma ? 1 : 0;
            kept += next[i];
        }
        if (next == used || kept < 6)
            break;
        used.swap(next);
    }
    out.mTol = o.mTolerance;
    std::int64_t nused = 0;
    for (std::uint8_t u : used)
        nused += u;
    out.mNumUsed = nused;
    out.mNumRejected = static_cast<std::int64_t>(n) - nused;
    return out;
}

// Seed mode (roadmap 7.1.1): the label marking the cells grown from the seeds
// (0) and the rest (-1), the cells themselves, and the absolute tolerance the
// growth ran at, which the patch is then fitted at.
std::vector<std::int64_t> fp_seed_labels(const Mesh& rMesh, const FitPrimitivesOptions& rOptions,
                                         std::size_t NCells, FitPrimitivesOptions& rFit,
                                         std::vector<std::int64_t>& rCells,
                                         std::vector<std::string>& rWarnings) {
    std::vector<std::uint8_t> selected(NCells, 1);
    if (!rOptions.mRegion.empty())
        selected = detail::cell_region_mask(rMesh, rOptions.mRegion, NCells, kFpPrefix);
    const detail::PrimitiveGrowSurface surf =
        detail::build_primitive_grow_surface(rMesh, selected, rOptions.mFeatureAngleDeg, kFpPrefix);
    std::vector<std::int64_t> label(NCells, -1);
    const std::size_t nf = surf.NumFaces();
    if (nf == 0) {
        rWarnings.emplace_back("the surface has no 2-D cells to grow from");
        return label;
    }
    detail::PrimitiveGrowInput in = surf.Input();
    detail::prepare_primitive_grow(in, nf, rOptions.mRelative ? rOptions.mTolerance : 1.0);
    if (!rOptions.mRelative)
        in.mTol = rOptions.mTolerance;
    const std::vector<std::vector<std::int64_t>> areas = detail::primitive_smooth_areas(in, nf);
    std::vector<std::int64_t> patch(nf, -1);
    std::vector<std::uint8_t> mark(nf, 0);
    const std::uint32_t types = rOptions.mTypes;
    const bool try_plane = (types & (1u << 1)) != 0;
    const std::uint32_t curved = types & ~(1u << 1);
    std::vector<std::int64_t> faces;
    for (std::size_t k = 0; k < rOptions.mSeeds.size(); ++k) {
        const std::int64_t f =
            detail::nearest_centroid_face(surf.mCentroid, rOptions.mSeeds[k].data());
        if (f < 0)
            continue;
        if (patch[static_cast<std::size_t>(f)] != -1)
            continue;  // an earlier seed's region holds it
        const std::vector<std::int64_t>& area =
            areas[static_cast<std::size_t>(in.mSmoothArea[static_cast<std::size_t>(f)])];
        const double area_total = detail::primitive_grow_area_of(in, area);
        // The seed's whole smooth area when one primitive explains it (step 1
        // of `segment`'s primitive method: a bore between sharp edges), else a
        // region grown exactly, as `segment` grows a CAD tessellation (a
        // fillet's tangent planes stay out), else within the whole tolerance,
        // as a scan's noise wants.
        detail::PrimitiveGrowRegion reg;
        bool whole_free = true;
        for (std::int64_t g : area)
            whole_free = whole_free && patch[static_cast<std::size_t>(g)] == -1;
        if (whole_free && detail::primitive_grow_fit(in, area, types).mType > 0)
            reg.mFaces = area;
        in.mGrowTol = detail::kPrimitiveGrowExactShare * in.mTol;
        if (reg.mFaces.empty())
            reg = detail::primitive_grow_from_seed(in, f, try_plane, curved, area_total, 0.0,
                                                   detail::kPrimitiveGrowExactShare, patch, mark);
        if (reg.mFaces.empty()) {
            in.mGrowTol = in.mTol;
            reg = detail::primitive_grow_from_seed(in, f, try_plane, curved, area_total, 0.0, 1.0,
                                                   patch, mark);
        }
        if (reg.mFaces.empty()) {
            rWarnings.push_back(
                "seed " + std::to_string(k) +
                ": no primitive of the asked types grows from the facet nearest it");
            continue;
        }
        for (std::int64_t g : reg.mFaces) {
            patch[static_cast<std::size_t>(g)] = 0;
            faces.push_back(g);
        }
    }
    std::sort(faces.begin(), faces.end());
    // Peel the rim: boundary faces whose normal turns away from the fitted
    // surface's, ring by ring.
    if (rOptions.mRimAngleDeg < 180.0 && !faces.empty()) {
        const double cos_rim = std::cos(rOptions.mRimAngleDeg * kFpPi / 180.0);
        std::vector<std::uint8_t> in_region(nf, 0);
        for (std::int64_t g : faces)
            in_region[static_cast<std::size_t>(g)] = 1;
        for (int pass = 0; pass < 8; ++pass) {
            const detail::PrimitiveFit fit = detail::primitive_grow_fit(in, faces, types);
            if (fit.mType == 0)
                break;
            std::vector<std::int64_t> kept;
            for (std::int64_t g : faces) {
                const std::size_t gs = static_cast<std::size_t>(g);
                const std::int64_t ring = surf.mRings.mStart[gs + 1] - surf.mRings.mStart[gs];
                const std::int64_t neighbours = surf.mAdjStart[gs + 1] - surf.mAdjStart[gs];
                bool rim = neighbours < ring;
                for (std::int64_t e = surf.mAdjStart[gs]; e < surf.mAdjStart[gs + 1] && !rim; ++e)
                    rim = !in_region[static_cast<std::size_t>(
                        surf.mAdj[static_cast<std::size_t>(e)].first)];
                double nrm[3];
                const bool peel =
                    rim &&
                    detail::primitive_normal(fit.mType, fit.mParams, &surf.mCentroid[gs * 3],
                                             nrm) &&
                    std::fabs(detail::MathUtils::Dot(nrm, &surf.mNormals[gs * 3])) < cos_rim;
                if (peel)
                    in_region[gs] = 0;
                else
                    kept.push_back(g);
            }
            if (kept.size() == faces.size() || kept.empty())
                break;
            faces.swap(kept);
        }
    }
    for (std::int64_t g : faces) {
        const std::int64_t c = surf.mRings.mCell[static_cast<std::size_t>(g)];
        label[static_cast<std::size_t>(c)] = 0;
        rCells.push_back(c);
    }
    std::sort(rCells.begin(), rCells.end());
    if (faces.empty())
        rWarnings.emplace_back("no seed grew a region");
    rFit.mTolerance = in.mTol;
    rFit.mRelative = false;
    return label;
}

}  // namespace

namespace detail {

PrimitiveFit fit_primitive_samples(const std::vector<double>& rPts,
                                   const std::vector<double>& rNormals,
                                   const std::vector<double>& rCentres, double Tol,
                                   std::uint32_t Types, double MinInlierFraction, int MaxIt) {
    PrimitiveFit out;
    const std::size_t n = rPts.size() / 3;
    if (n == 0 || !(Tol > 0.0))
        return out;
    FpPatch patch;
    patch.mPts = rPts;
    fp_mean(patch.mPts, patch.mShift);
    double lo[3], hi[3];
    for (int k = 0; k < 3; ++k)
        lo[k] = hi[k] = patch.mPts[static_cast<std::size_t>(k)];
    for (std::size_t i = 0; i < n; ++i)
        for (int k = 0; k < 3; ++k) {
            double& v = patch.mPts[i * 3 + static_cast<std::size_t>(k)];
            lo[k] = v < lo[k] ? v : lo[k];
            hi[k] = v > hi[k] ? v : hi[k];
            v -= patch.mShift[k];
        }
    const double ext[3] = {hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]};
    patch.mDiag = std::sqrt(detail::MathUtils::Dot(ext, ext));
    patch.mNormals = rNormals;
    patch.mCentres = rCentres;
    for (std::size_t f = 0; f < patch.mCentres.size() / 3; ++f)
        for (int k = 0; k < 3; ++k)
            patch.mCentres[f * 3 + static_cast<std::size_t>(k)] -= patch.mShift[k];
    FitPrimitivesOptions options;
    options.mTolerance = Tol;
    options.mRelative = false;
    options.mMinInlierFraction = MinInlierFraction;
    options.mTypes = Types;
    options.mMaxIterations = MaxIt;
    const FpOutcome o = fp_fit_patch(patch, options);
    out.mType = o.mType;
    for (std::size_t k = 0; k < kFpParams; ++k)
        out.mParams[k] = o.mParams[k];
    out.mRms = o.mRms;
    out.mMax = o.mMax;
    return out;
}

}  // namespace detail

FitPrimitivesResult fit_primitives(const Mesh& rMesh, const FitPrimitivesOptions& rOptions) {
    if (!(rOptions.mTolerance > 0.0) || !std::isfinite(rOptions.mTolerance))
        throw std::invalid_argument(std::string(kFpPrefix) +
                                    "the tolerance must be finite and > 0");
    if (!(rOptions.mMinInlierFraction >= 0.0 && rOptions.mMinInlierFraction <= 1.0))
        throw std::invalid_argument(std::string(kFpPrefix) +
                                    "min_inlier_fraction must lie in [0, 1]");
    if ((rOptions.mTypes & kFitAllTypes) == 0 || (rOptions.mTypes & ~kFitAllTypes) != 0)
        throw std::invalid_argument(std::string(kFpPrefix) +
                                    "types must name at least one of plane, sphere, cylinder, "
                                    "cone, torus and nothing else");
    if (rOptions.mMaxIterations < 1)
        throw std::invalid_argument(std::string(kFpPrefix) + "max_iterations must be >= 1");
    if (rOptions.mUsedPoints < FitUsedPoints::Auto || rOptions.mUsedPoints > FitUsedPoints::Sigma1)
        throw std::invalid_argument(std::string(kFpPrefix) + "unknown used_points");
    if (!(rOptions.mFeatureAngleDeg >= 0.0 && rOptions.mFeatureAngleDeg <= 180.0) ||
        !(rOptions.mRimAngleDeg >= 0.0 && rOptions.mRimAngleDeg <= 180.0))
        throw std::invalid_argument(std::string(kFpPrefix) +
                                    "feature_angle and rim_angle must lie in [0, 180] degrees");
    for (const std::array<double, 3>& rSeed : rOptions.mSeeds)
        for (double v : rSeed)
            if (!std::isfinite(v))
                throw std::invalid_argument(std::string(kFpPrefix) + "a seed must be finite");
    const bool seeded = !rOptions.mSeeds.empty();
    if (seeded && rOptions.mOutputRegion.empty())
        throw std::invalid_argument(std::string(kFpPrefix) + "the output region needs a name");
    const bool feature_mode = seeded || !rOptions.mRegion.empty();
    const FitUsedPoints used_mode =
        rOptions.mUsedPoints == FitUsedPoints::Auto
            ? (feature_mode ? FitUsedPoints::Sigma3 : FitUsedPoints::All)
            : rOptions.mUsedPoints;
    const double trim_k = used_mode == FitUsedPoints::Sigma3   ? 3.0
                          : used_mode == FitUsedPoints::Sigma2 ? 2.0
                                                               : 1.0;
    FitPrimitivesResult result;
    FitPrimitivesOptions fit_opts = rOptions;

    const std::vector<std::int64_t> bases = detail::block_bases(rMesh);
    const std::size_t ncells = static_cast<std::size_t>(detail::total_cells(bases));
    std::vector<std::int64_t> label(ncells, 0);
    std::vector<std::int64_t> feature_cells;
    if (seeded) {
        label = fp_seed_labels(rMesh, rOptions, ncells, fit_opts, feature_cells, result.mWarnings);
    } else if (!rOptions.mRegion.empty()) {
        const std::vector<std::uint8_t> mask =
            detail::cell_region_mask(rMesh, rOptions.mRegion, ncells, kFpPrefix);
        for (std::size_t c = 0; c < ncells; ++c)
            label[c] = mask[c] ? 0 : -1;
    } else if (!rOptions.mPatchData.empty()) {
        if (!rMesh.HasCellData(rOptions.mPatchData))
            throw std::invalid_argument(std::string(kFpPrefix) + "no cell data named '" +
                                        rOptions.mPatchData +
                                        "' (run segment first, or name the label array)");
        for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
            const NDArray& a = rMesh.CellData(rOptions.mPatchData, b);
            const std::size_t nc = rMesh.Cells(b).NumCells();
            if (a.Size() != nc)
                throw std::invalid_argument(std::string(kFpPrefix) + "cell data '" +
                                            rOptions.mPatchData + "' must hold one label per cell");
            if (a.Dtype() == DType::Float32 || a.Dtype() == DType::Float64)
                throw std::invalid_argument(std::string(kFpPrefix) + "cell data '" +
                                            rOptions.mPatchData + "' must be an integer label");
            const std::size_t base = static_cast<std::size_t>(bases[b]);
            for (std::size_t c = 0; c < nc; ++c)
                label[base + c] = detail::read_int(a, c);
        }
    }

    std::vector<std::uint8_t> all(ncells, 1);
    const detail::SurfaceRings rings = detail::surface_rings(rMesh, all, kFpPrefix);
    const std::vector<double> xyz = detail::surface_xyz(rMesh);
    const std::size_t nf = rings.NumFaces();

    // Patch rows: one per distinct non-negative label on a surface cell.
    std::map<std::int64_t, std::size_t> row_of;
    for (std::size_t f = 0; f < nf; ++f) {
        const std::int64_t l = label[static_cast<std::size_t>(rings.mCell[f])];
        if (l >= 0)
            row_of.emplace(l, 0);
    }
    std::vector<std::int64_t> ids;
    for (auto& kv : row_of) {
        kv.second = ids.size();
        ids.push_back(kv.first);
    }
    const std::size_t m = ids.size();
    std::vector<std::vector<std::size_t>> faces_of(m);
    for (std::size_t f = 0; f < nf; ++f) {
        const std::int64_t l = label[static_cast<std::size_t>(rings.mCell[f])];
        if (l >= 0)
            faces_of[row_of[l]].push_back(f);
    }

    std::vector<FpOutcome> outcome(m);
    parallel_for(m, [&](std::size_t r) {
        FpPatch patch;
        std::vector<std::int64_t> verts;
        for (std::size_t f : faces_of[r])
            for (std::int64_t k = rings.mStart[f]; k < rings.mStart[f + 1]; ++k)
                verts.push_back(rings.mNodes[static_cast<std::size_t>(k)]);
        std::sort(verts.begin(), verts.end());
        verts.erase(std::unique(verts.begin(), verts.end()), verts.end());
        patch.mPts.resize(verts.size() * 3);
        for (std::size_t i = 0; i < verts.size(); ++i)
            for (int k = 0; k < 3; ++k)
                patch.mPts[i * 3 + static_cast<std::size_t>(k)] =
                    xyz[static_cast<std::size_t>(verts[i]) * 3 + static_cast<std::size_t>(k)];
        if (!verts.empty()) {
            fp_mean(patch.mPts, patch.mShift);
            double lo[3], hi[3];
            for (int k = 0; k < 3; ++k)
                lo[k] = hi[k] = patch.mPts[static_cast<std::size_t>(k)];
            for (std::size_t i = 0; i < verts.size(); ++i)
                for (int k = 0; k < 3; ++k) {
                    double& v = patch.mPts[i * 3 + static_cast<std::size_t>(k)];
                    lo[k] = v < lo[k] ? v : lo[k];
                    hi[k] = v > hi[k] ? v : hi[k];
                    v -= patch.mShift[k];
                }
            const double ext[3] = {hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]};
            patch.mDiag = std::sqrt(detail::MathUtils::Dot(ext, ext));
        }
        for (std::size_t f : faces_of[r]) {
            const std::int64_t* ring = rings.mNodes.data() + rings.mStart[f];
            const std::size_t size =
                static_cast<std::size_t>(rings.mStart[f + 1] - rings.mStart[f]);
            double nrm[3];
            detail::ring_unit_normal(xyz.data(), ring, size, nrm);
            double c[3] = {0.0, 0.0, 0.0};
            for (std::size_t k = 0; k < size; ++k)
                for (std::size_t d = 0; d < 3; ++d)
                    c[d] += xyz[static_cast<std::size_t>(ring[k]) * 3 + d];
            for (std::size_t d = 0; d < 3; ++d) {
                patch.mNormals.push_back(nrm[d]);
                patch.mCentres.push_back(c[d] / static_cast<double>(size) - patch.mShift[d]);
            }
        }
        std::vector<std::uint8_t> used;
        if (used_mode == FitUsedPoints::All) {
            outcome[r] = fp_fit_patch(patch, fit_opts);
            used.assign(verts.size(), 1);
            outcome[r].mNumUsed = static_cast<std::int64_t>(verts.size());
        } else {
            outcome[r] = fp_fit_patch_trimmed(patch, fit_opts, trim_k, used);
        }
        if (!rOptions.mExtendedReport)
            return;
        // Sigma, form error and extremes over the points used, in the world
        // frame the parameters are in.
        std::vector<double> world;
        for (std::size_t i = 0; i < verts.size(); ++i)
            if (used[i])
                for (int k = 0; k < 3; ++k)
                    world.push_back(patch.mPts[i * 3 + static_cast<std::size_t>(k)] +
                                    patch.mShift[k]);
        const std::size_t nused = world.size() / 3;
        FpOutcome& o = outcome[r];
        if (o.mType > 0 && nused > 0) {
            double sum = 0.0;
            for (std::size_t i = 0; i < nused; ++i) {
                const double d =
                    detail::primitive_signed_distance(o.mType, o.mParams, &world[i * 3]);
                sum += d * d;
            }
            o.mSigma = std::sqrt(sum / static_cast<double>(nused));
        }
        o.mForm = detail::primitive_form_metrics(o.mType, o.mParams, world.data(), nused);
    });

    result.mNumPatches = static_cast<std::int64_t>(m);
    for (const FpOutcome& o : outcome) {
        switch (o.mType) {
            case 1:
                ++result.mNumPlane;
                break;
            case 2:
                ++result.mNumSphere;
                break;
            case 3:
                ++result.mNumCylinder;
                break;
            case 4:
                ++result.mNumCone;
                break;
            case 5:
                ++result.mNumTorus;
                break;
            default:
                ++result.mNumUnfitted;
                break;
        }
    }

    // Per global cell, its patch row (or none).
    std::vector<std::int64_t> cell_row(ncells, -1);
    for (std::size_t f = 0; f < nf; ++f) {
        const std::size_t c = static_cast<std::size_t>(rings.mCell[f]);
        if (label[c] >= 0)
            cell_row[c] = static_cast<std::int64_t>(row_of[label[c]]);
    }

    Mesh& out = result.mMesh;
    out = detail::clone_mesh(rMesh,
                             [](DataLocation, const std::string&, std::string&) { return true; });
    const double nan = std::numeric_limits<double>::quiet_NaN();
    std::vector<NDArray> type_blocks, rms_blocks, max_blocks;
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const std::size_t base = static_cast<std::size_t>(bases[b]);
        const std::size_t nc = rMesh.Cells(b).NumCells();
        NDArray t = NDArray::Uninit(DType::Int64, {nc});
        NDArray r = NDArray::Uninit(DType::Float64, {nc});
        NDArray x = NDArray::Uninit(DType::Float64, {nc});
        for (std::size_t c = 0; c < nc; ++c) {
            const std::int64_t row = cell_row[base + c];
            const FpOutcome* o = row >= 0 ? &outcome[static_cast<std::size_t>(row)] : nullptr;
            t.As<std::int64_t>()[c] = o ? o->mType : -1;
            r.As<double>()[c] = o ? o->mRms : nan;
            x.As<double>()[c] = o ? o->mMax : nan;
        }
        type_blocks.push_back(std::move(t));
        rms_blocks.push_back(std::move(r));
        max_blocks.push_back(std::move(x));
    }
    out.AddCellData(kFitTypeName, std::move(type_blocks));
    out.AddCellData(kFitRmsName, std::move(rms_blocks));
    out.AddCellData(kFitMaxName, std::move(max_blocks));

    NDArray pid = NDArray::Uninit(DType::Int64, {m});
    NDArray ptype = NDArray::Uninit(DType::Int64, {m});
    NDArray params = NDArray::Uninit(DType::Float64, {m, kFpParams});
    NDArray prms = NDArray::Uninit(DType::Float64, {m});
    NDArray pmax = NDArray::Uninit(DType::Float64, {m});
    NDArray pfrac = NDArray::Uninit(DType::Float64, {m});
    NDArray ptol = NDArray::Uninit(DType::Float64, {m});
    NDArray psigma = NDArray::Uninit(DType::Float64, {m});
    NDArray pnused = NDArray::Uninit(DType::Int64, {m});
    NDArray pnrej = NDArray::Uninit(DType::Int64, {m});
    NDArray pform = NDArray::Uninit(DType::Float64, {m});
    NDArray pdmin = NDArray::Uninit(DType::Float64, {m});
    NDArray pdmax = NDArray::Uninit(DType::Float64, {m});
    NDArray plen = NDArray::Uninit(DType::Float64, {m});
    NDArray parc = NDArray::Uninit(DType::Float64, {m});
    NDArray pext = NDArray::Uninit(DType::Float64, {m, 2});
    for (std::size_t r = 0; r < m; ++r) {
        const detail::PrimitiveFormMetrics& f = outcome[r].mForm;
        psigma.As<double>()[r] = outcome[r].mSigma;
        pnused.As<std::int64_t>()[r] = outcome[r].mNumUsed;
        pnrej.As<std::int64_t>()[r] = outcome[r].mNumRejected;
        pform.As<double>()[r] = f.mFormError;
        pdmin.As<double>()[r] = f.mDiameterMin;
        pdmax.As<double>()[r] = f.mDiameterMax;
        plen.As<double>()[r] = f.mLength;
        parc.As<double>()[r] = f.mArcDeg;
        pext.As<double>()[r * 2] = f.mExtent[0];
        pext.As<double>()[r * 2 + 1] = f.mExtent[1];
    }
    for (std::size_t r = 0; r < m; ++r) {
        pid.As<std::int64_t>()[r] = ids[r];
        ptype.As<std::int64_t>()[r] = outcome[r].mType;
        for (std::size_t k = 0; k < kFpParams; ++k)
            params.As<double>()[r * kFpParams + k] = outcome[r].mParams[k];
        prms.As<double>()[r] = outcome[r].mRms;
        pmax.As<double>()[r] = outcome[r].mMax;
        pfrac.As<double>()[r] = outcome[r].mFraction;
        ptol.As<double>()[r] = outcome[r].mTol;
    }
    out.AddFieldData("fit:patch_id", std::move(pid));
    out.AddFieldData("fit:patch_type", std::move(ptype));
    out.AddFieldData("fit:patch_params", std::move(params));
    out.AddFieldData("fit:patch_rms", std::move(prms));
    out.AddFieldData("fit:patch_max", std::move(pmax));
    out.AddFieldData("fit:patch_fraction", std::move(pfrac));
    out.AddFieldData("fit:patch_tolerance", std::move(ptol));
    if (rOptions.mExtendedReport) {
        out.AddFieldData(kFitPatchFormErrorName, std::move(pform));
        out.AddFieldData(kFitPatchDiameterMinName, std::move(pdmin));
        out.AddFieldData(kFitPatchDiameterMaxName, std::move(pdmax));
        out.AddFieldData(kFitPatchLengthName, std::move(plen));
        out.AddFieldData(kFitPatchArcName, std::move(parc));
        out.AddFieldData(kFitPatchExtentName, std::move(pext));
        out.AddFieldData(kFitPatchSigmaName, std::move(psigma));
        out.AddFieldData(kFitPatchNumUsedName, std::move(pnused));
        out.AddFieldData(kFitPatchNumRejectedName, std::move(pnrej));
    }
    if (seeded) {
        result.mNumFeatureCells = static_cast<std::int64_t>(feature_cells.size());
        for (std::size_t k = out.NumRegions(); k-- > 0;)
            if (out.Region(k).mKind == RegionKind::Cell &&
                out.Region(k).mName == rOptions.mOutputRegion)
                out.RemoveRegion(k);
        NDArray entries = NDArray::Uninit(DType::Int64, {feature_cells.size()});
        std::copy(feature_cells.begin(), feature_cells.end(), entries.As<std::int64_t>());
        meshioplusplus::Region region(rOptions.mOutputRegion, RegionKind::Cell, std::move(entries));
        region.Canonicalize();
        out.AddRegion(std::move(region));
    }
    return result;
}

}  // namespace meshioplusplus

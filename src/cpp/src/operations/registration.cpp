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
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/registration.hpp"
#include "meshioplusplus/detail/spatial_hash.hpp"
#include "meshioplusplus/detail/surface_distance.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/types.hpp"

namespace meshioplusplus {

namespace {

using RegVec = std::array<double, 3>;
using RegMat = std::array<double, 9>;  // row-major 3x3

RegVec reg_apply(const RegMat& rR, const RegVec& rT, const RegVec& rP) {
    return {rR[0] * rP[0] + rR[1] * rP[1] + rR[2] * rP[2] + rT[0],
            rR[3] * rP[0] + rR[4] * rP[1] + rR[5] * rP[2] + rT[1],
            rR[6] * rP[0] + rR[7] * rP[1] + rR[8] * rP[2] + rT[2]};
}

RegMat reg_mul(const RegMat& rA, const RegMat& rB) {
    RegMat c{};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                c[static_cast<std::size_t>(3 * i + j)] += rA[static_cast<std::size_t>(3 * i + k)] *
                                                          rB[static_cast<std::size_t>(3 * k + j)];
    return c;
}

double reg_dot(const RegVec& rA, const RegVec& rB) {
    return rA[0] * rB[0] + rA[1] * rB[1] + rA[2] * rB[2];
}

RegVec reg_cross(const RegVec& rA, const RegVec& rB) {
    return {rA[1] * rB[2] - rA[2] * rB[1], rA[2] * rB[0] - rA[0] * rB[2],
            rA[0] * rB[1] - rA[1] * rB[0]};
}

/// Eigen-decomposition of a symmetric N x N matrix (row-major) by cyclic
/// Jacobi: eigenvalues in `rW`, eigenvectors as the columns of `rV`.
template <std::size_t N>
void reg_jacobi(std::array<double, N * N> A, std::array<double, N>& rW,
                std::array<double, N * N>& rV) {
    rV.fill(0.0);
    for (std::size_t i = 0; i < N; ++i)
        rV[i * N + i] = 1.0;
    for (int sweep = 0; sweep < 100; ++sweep) {
        double off = 0.0, diag = 0.0;
        for (std::size_t i = 0; i < N; ++i) {
            diag += A[i * N + i] * A[i * N + i];
            for (std::size_t j = i + 1; j < N; ++j)
                off += A[i * N + j] * A[i * N + j];
        }
        if (off <= 1e-30 * diag || off == 0.0)
            break;
        for (std::size_t p = 0; p < N; ++p)
            for (std::size_t q = p + 1; q < N; ++q) {
                const double apq = A[p * N + q];
                if (apq == 0.0)
                    continue;
                const double theta = (A[q * N + q] - A[p * N + p]) / (2.0 * apq);
                const double t = (theta >= 0.0 ? 1.0 : -1.0) /
                                 (std::abs(theta) + std::sqrt(theta * theta + 1.0));
                const double c = 1.0 / std::sqrt(t * t + 1.0), s = t * c;
                for (std::size_t k = 0; k < N; ++k) {  // A <- A J
                    const double akp = A[k * N + p], akq = A[k * N + q];
                    A[k * N + p] = c * akp - s * akq;
                    A[k * N + q] = s * akp + c * akq;
                }
                for (std::size_t k = 0; k < N; ++k) {  // A <- J^T A
                    const double apk = A[p * N + k], aqk = A[q * N + k];
                    A[p * N + k] = c * apk - s * aqk;
                    A[q * N + k] = s * apk + c * aqk;
                }
                for (std::size_t k = 0; k < N; ++k) {
                    const double vkp = rV[k * N + p], vkq = rV[k * N + q];
                    rV[k * N + p] = c * vkp - s * vkq;
                    rV[k * N + q] = s * vkp + c * vkq;
                }
            }
    }
    for (std::size_t i = 0; i < N; ++i)
        rW[i] = A[i * N + i];
}

/// The rotation of a unit quaternion (w, x, y, z).
RegMat reg_quaternion(double W, double X, double Y, double Z) {
    return {W * W + X * X - Y * Y - Z * Z, 2 * (X * Y - W * Z),
            2 * (X * Z + W * Y),           2 * (X * Y + W * Z),
            W * W - X * X + Y * Y - Z * Z, 2 * (Y * Z - W * X),
            2 * (X * Z - W * Y),           2 * (Y * Z + W * X),
            W * W - X * X - Y * Y + Z * Z};
}

/// The rotation by the vector `rW` (axis times angle), exactly.
RegMat reg_rodrigues(const RegVec& rW) {
    const double angle = std::sqrt(reg_dot(rW, rW));
    if (!(angle > 0.0))
        return {1, 0, 0, 0, 1, 0, 0, 0, 1};
    const double h = 0.5 * angle, s = std::sin(h) / angle;
    return reg_quaternion(std::cos(h), rW[0] * s, rW[1] * s, rW[2] * s);
}

/// The rotation angle of a rotation matrix.
double reg_angle(const RegMat& rR) {
    const double c = std::clamp(0.5 * (rR[0] + rR[4] + rR[8] - 1.0), -1.0, 1.0);
    // Near zero acos loses half the digits; the skew part keeps them.
    const double s =
        0.5 * std::sqrt((rR[7] - rR[5]) * (rR[7] - rR[5]) + (rR[2] - rR[6]) * (rR[2] - rR[6]) +
                        (rR[3] - rR[1]) * (rR[3] - rR[1]));
    return std::atan2(s, c);
}

struct RegPairs {
    std::vector<RegVec> mP;  ///< moved source points
    std::vector<RegVec> mQ;  ///< their closest target points
    std::vector<RegVec> mN;  ///< the target's unit normals there (point-to-plane)
    double mSumSq = 0.0;
};

/// The target: a triangle soup with per-cell normals, or a cloud in a grid.
class RegTarget {
public:
    RegTarget(const Mesh& rTarget, bool NeedNormals) {
        const auto& dims = topological_dimension();
        bool surface = false;
        for (const auto cb : rTarget.CellRange()) {
            const auto it = dims.find(cb.Type());
            if (it != dims.end() && it->second >= 2 && cb.NumCells() > 0)
                surface = true;
        }
        const std::size_t n = rTarget.NumPoints();
        const std::size_t dim = rTarget.PointDim();
        mPoints.assign(n, RegVec{0.0, 0.0, 0.0});
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t k = 0; k < dim && k < 3; ++k) {
                const double v = detail::read_double(rTarget.Points(), i * dim + k);
                if (!std::isfinite(v))
                    throw std::invalid_argument("meshio++: register_rigid: target point " +
                                                std::to_string(i) + " has a non-finite coordinate");
                mPoints[i][k] = v;
            }
        if (n == 0)
            throw std::invalid_argument("meshio++: register_rigid: the target has no points");
        if (surface) {
            mSoup = detail::build_triangle_soup(rTarget, "");
            mQuery = detail::build_distance_query(mSoup, SurfaceDistanceOptions{});
            mSurface = true;
            if (NeedNormals) {
                std::int64_t cells = 0;
                for (const auto cb : rTarget.CellRange())
                    cells += static_cast<std::int64_t>(cb.NumCells());
                mCellNormal.assign(static_cast<std::size_t>(cells), RegVec{0.0, 0.0, 0.0});
                for (std::size_t t = 0; t < mSoup.NumTriangles(); ++t) {
                    RegVec& c = mCellNormal[static_cast<std::size_t>(mSoup.mSourceCell[t])];
                    for (std::size_t a = 0; a < 3; ++a)
                        c[a] += mQuery.mFaceNormal[t][a];
                }
                for (RegVec& c : mCellNormal) {
                    const double len = std::sqrt(reg_dot(c, c));
                    if (len > 0.0)
                        c = {c[0] / len, c[1] / len, c[2] / len};
                    else
                        c = {std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0};
                }
            }
            return;
        }
        if (NeedNormals) {
            if (!rTarget.HasPointData("normals"))
                throw std::invalid_argument(
                    "meshio++: register_rigid: point_to_plane against a point cloud needs its "
                    "'normals' point data; run estimate_normals on the target first");
            const NDArray& nrm = rTarget.PointData("normals");
            if (nrm.Size() != n * 3)
                throw std::invalid_argument(
                    "meshio++: register_rigid: the target's 'normals' must be (n, 3)");
            mPointNormal.resize(n);
            for (std::size_t i = 0; i < n; ++i)
                for (std::size_t a = 0; a < 3; ++a)
                    mPointNormal[i][a] = detail::read_double(nrm, i * 3 + a);
        }
        // A grid of about four points per cell over the cloud's surface
        // (half its box's surface area; a speed choice, never the answer).
        RegVec lo = mPoints[0], hi = mPoints[0];
        for (const RegVec& q : mPoints)
            for (std::size_t a = 0; a < 3; ++a) {
                lo[a] = std::min(lo[a], q[a]);
                hi[a] = std::max(hi[a], q[a]);
            }
        const RegVec e{hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]};
        const double area = e[0] * e[1] + e[1] * e[2] + e[2] * e[0];
        const double extent = std::max({e[0], e[1], e[2]});
        double cell = area > 0.0 ? std::sqrt(4.0 * area / static_cast<double>(n)) : extent;
        cell = std::max(cell, std::max(extent, 1.0) * 1e-9);
        mGrid = detail::SpatialGrid(cell);
        for (std::size_t i = 0; i < n; ++i)
            mGrid.Insert(mGrid.KeyOf(mPoints[i].data()), static_cast<std::int64_t>(i));
    }

    /// The closest target point (and normal) of each moved source point.
    void Match(const std::vector<RegVec>& rMoved, std::vector<RegVec>& rQ, std::vector<RegVec>& rN,
               bool NeedNormals) const {
        const std::size_t m = rMoved.size();
        rQ.assign(m, RegVec{0.0, 0.0, 0.0});
        rN.assign(NeedNormals ? m : 0, RegVec{0.0, 0.0, 0.0});
        if (mSurface) {
            const std::vector<detail::ClosestPointHit> hits =
                detail::query_closest_points(mQuery, rMoved);
            for (std::size_t i = 0; i < m; ++i) {
                rQ[i] = hits[i].mPoint;
                if (NeedNormals)
                    rN[i] = mCellNormal[static_cast<std::size_t>(hits[i].mSourceCell)];
            }
            return;
        }
        parallel_for(m, [&](std::size_t i) {
            const std::int64_t j = Nearest(rMoved[i]);
            rQ[i] = mPoints[static_cast<std::size_t>(j)];
            if (NeedNormals)
                rN[i] = mPointNormal[static_cast<std::size_t>(j)];
        });
    }

private:
    std::vector<RegVec> mPoints;
    std::vector<RegVec> mPointNormal;
    std::vector<RegVec> mCellNormal;
    detail::TriangleSoup mSoup;
    detail::DistanceQuery mQuery;
    detail::SpatialGrid mGrid{1.0};
    bool mSurface = false;

    /// The nearest cloud point, by expanding shells (ties to the lower id).
    std::int64_t Nearest(const RegVec& rP) const {
        const detail::GridKey k = mGrid.KeyOf(rP.data());
        const detail::GridKey& lo = mGrid.OccupiedLo();
        const detail::GridKey& hi = mGrid.OccupiedHi();
        const auto cap = [](std::int64_t v, std::int64_t a, std::int64_t b) {
            return std::max(v > a ? v - a : a - v, v > b ? v - b : b - v);
        };
        const auto start = [](std::int64_t v, std::int64_t a, std::int64_t b) {
            return v < a ? a - v : (v > b ? v - b : std::int64_t{0});
        };
        const std::int64_t rcap =
            std::max({cap(k.x, lo.x, hi.x), cap(k.y, lo.y, hi.y), cap(k.z, lo.z, hi.z)});
        const std::int64_t r0 =
            std::max({start(k.x, lo.x, hi.x), start(k.y, lo.y, hi.y), start(k.z, lo.z, hi.z)});
        const double cell = mGrid.CellSize();
        double best_d2 = std::numeric_limits<double>::infinity();
        std::int64_t best = -1;
        for (std::int64_t r = r0; r <= rcap; ++r) {
            if (best >= 0 && r >= 1) {
                const double ring = static_cast<double>(r - 1) * cell;
                if (ring * ring > best_d2)
                    break;
            }
            mGrid.ForEachInShell(k, r, [&](const std::vector<std::int64_t>& rIds) {
                for (const std::int64_t id : rIds) {
                    const RegVec& q = mPoints[static_cast<std::size_t>(id)];
                    const double dx = rP[0] - q[0], dy = rP[1] - q[1], dz = rP[2] - q[2];
                    const double d2 = dx * dx + dy * dy + dz * dz;
                    if (d2 < best_d2 || (d2 == best_d2 && id < best)) {
                        best_d2 = d2;
                        best = id;
                    }
                }
            });
        }
        return best;
    }
};

/// The pairs under the transform (R, t), within `MaxDistance` when positive,
/// with finite normals when they are needed.
RegPairs reg_pairs(const RegTarget& rTarget, const std::vector<RegVec>& rSource, const RegMat& rR,
                   const RegVec& rT, double MaxDistance, bool NeedNormals) {
    std::vector<RegVec> moved(rSource.size());
    parallel_for(rSource.size(), [&](std::size_t i) { moved[i] = reg_apply(rR, rT, rSource[i]); });
    std::vector<RegVec> q, nrm;
    rTarget.Match(moved, q, nrm, NeedNormals);
    RegPairs out;
    for (std::size_t i = 0; i < moved.size(); ++i) {
        const RegVec d{q[i][0] - moved[i][0], q[i][1] - moved[i][1], q[i][2] - moved[i][2]};
        const double d2 = reg_dot(d, d);
        if (MaxDistance > 0.0 && d2 > MaxDistance * MaxDistance)
            continue;
        if (NeedNormals &&
            !(std::isfinite(nrm[i][0]) && std::isfinite(nrm[i][1]) && std::isfinite(nrm[i][2])))
            continue;
        out.mP.push_back(moved[i]);
        out.mQ.push_back(q[i]);
        if (NeedNormals)
            out.mN.push_back(nrm[i]);
        out.mSumSq += d2;
    }
    return out;
}

/// Horn's closed-form best rotation and translation for point pairs.
void reg_point_to_point(const RegPairs& rPairs, RegMat& rR, RegVec& rT) {
    const double m = static_cast<double>(rPairs.mP.size());
    RegVec pc{0, 0, 0}, qc{0, 0, 0};
    for (std::size_t i = 0; i < rPairs.mP.size(); ++i)
        for (std::size_t a = 0; a < 3; ++a) {
            pc[a] += rPairs.mP[i][a];
            qc[a] += rPairs.mQ[i][a];
        }
    for (std::size_t a = 0; a < 3; ++a) {
        pc[a] /= m;
        qc[a] /= m;
    }
    RegMat s{};  // S_ab = sum (p - pc)_a (q - qc)_b
    for (std::size_t i = 0; i < rPairs.mP.size(); ++i)
        for (std::size_t a = 0; a < 3; ++a)
            for (std::size_t b = 0; b < 3; ++b)
                s[3 * a + b] += (rPairs.mP[i][a] - pc[a]) * (rPairs.mQ[i][b] - qc[b]);
    const double sxx = s[0], sxy = s[1], sxz = s[2], syx = s[3], syy = s[4], syz = s[5], szx = s[6],
                 szy = s[7], szz = s[8];
    const std::array<double, 16> n{sxx + syy + szz,  syz - szy, szx - sxz, sxy - syx, syz - szy,
                                   sxx - syy - szz,  sxy + syx, szx + sxz, szx - sxz, sxy + syx,
                                   -sxx + syy - szz, syz + szy, sxy - syx, szx + sxz, syz + szy,
                                   -sxx - syy + szz};
    std::array<double, 4> w;
    std::array<double, 16> v;
    reg_jacobi<4>(n, w, v);
    std::size_t best = 0;
    for (std::size_t j = 1; j < 4; ++j)
        if (w[j] > w[best])
            best = j;
    double qw = v[0 * 4 + best], qx = v[1 * 4 + best], qy = v[2 * 4 + best], qz = v[3 * 4 + best];
    const double len = std::sqrt(qw * qw + qx * qx + qy * qy + qz * qz);
    qw /= len;
    qx /= len;
    qy /= len;
    qz /= len;
    rR = reg_quaternion(qw, qx, qy, qz);
    const RegVec rp = reg_apply(rR, {0, 0, 0}, pc);
    rT = {qc[0] - rp[0], qc[1] - rp[1], qc[2] - rp[2]};
}

/// The linearised point-to-plane step, applied as an exact rotation about the
/// pairs' centroid.
void reg_point_to_plane(const RegPairs& rPairs, RegMat& rR, RegVec& rT) {
    const double m = static_cast<double>(rPairs.mP.size());
    RegVec c{0, 0, 0};
    for (const RegVec& p : rPairs.mP)
        for (std::size_t a = 0; a < 3; ++a)
            c[a] += p[a];
    for (std::size_t a = 0; a < 3; ++a)
        c[a] /= m;
    double scale = 0.0;
    for (const RegVec& p : rPairs.mP) {
        const RegVec d{p[0] - c[0], p[1] - c[1], p[2] - c[2]};
        scale += reg_dot(d, d);
    }
    scale = std::sqrt(scale / m);
    if (!(scale > 0.0))
        scale = 1.0;
    std::array<double, 36> a{};
    std::array<double, 6> g{};
    for (std::size_t i = 0; i < rPairs.mP.size(); ++i) {
        const RegVec p{(rPairs.mP[i][0] - c[0]) / scale, (rPairs.mP[i][1] - c[1]) / scale,
                       (rPairs.mP[i][2] - c[2]) / scale};
        const RegVec q{(rPairs.mQ[i][0] - c[0]) / scale, (rPairs.mQ[i][1] - c[1]) / scale,
                       (rPairs.mQ[i][2] - c[2]) / scale};
        const RegVec& nn = rPairs.mN[i];
        const RegVec x = reg_cross(p, nn);
        const std::array<double, 6> row{x[0], x[1], x[2], nn[0], nn[1], nn[2]};
        const double b = reg_dot({q[0] - p[0], q[1] - p[1], q[2] - p[2]}, nn);
        for (std::size_t r = 0; r < 6; ++r) {
            g[r] += row[r] * b;
            for (std::size_t s = 0; s < 6; ++s)
                a[r * 6 + s] += row[r] * row[s];
        }
    }
    // A pseudo-inverse: directions the planes do not constrain (sliding along
    // a flat or a cylindrical target) are left still.
    std::array<double, 6> w;
    std::array<double, 36> v;
    reg_jacobi<6>(a, w, v);
    const double wmax = *std::max_element(w.begin(), w.end());
    std::array<double, 6> x{};
    for (std::size_t k = 0; k < 6; ++k) {
        if (!(w[k] > 1e-12 * wmax))
            continue;
        double proj = 0.0;
        for (std::size_t r = 0; r < 6; ++r)
            proj += v[r * 6 + k] * g[r];
        for (std::size_t r = 0; r < 6; ++r)
            x[r] += v[r * 6 + k] * proj / w[k];
    }
    rR = reg_rodrigues({x[0], x[1], x[2]});
    // x -> R (x - c) + c + s dt
    const RegVec rc = reg_apply(rR, {0, 0, 0}, c);
    rT = {c[0] - rc[0] + scale * x[3], c[1] - rc[1] + scale * x[4], c[2] - rc[2] + scale * x[5]};
}

}  // namespace

RegistrationMethod registration_method_from_name(const std::string& rName) {
    if (rName == "point_to_point")
        return RegistrationMethod::PointToPoint;
    if (rName == "point_to_plane")
        return RegistrationMethod::PointToPlane;
    throw std::invalid_argument("meshio++: register_rigid: unknown method '" + rName +
                                "' (expected 'point_to_point' or 'point_to_plane')");
}

RegistrationResult register_rigid(const Mesh& rSource, const Mesh& rTarget,
                                  const RegistrationOptions& rOptions) {
    if (rOptions.mMaxIterations < 1)
        throw std::invalid_argument("meshio++: register_rigid: max_iterations must be at least 1");
    if (!(rOptions.mTolerance >= 0.0) || !(rOptions.mMaxDistance >= 0.0))
        throw std::invalid_argument(
            "meshio++: register_rigid: tolerance and max_distance must not be negative");
    const std::size_t n = rSource.NumPoints();
    if (n == 0)
        throw std::invalid_argument("meshio++: register_rigid: the source has no points");
    const std::size_t dim = rSource.PointDim();
    std::vector<RegVec> src(n, RegVec{0.0, 0.0, 0.0});
    RegVec lo{0, 0, 0}, hi{0, 0, 0};
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t k = 0; k < dim && k < 3; ++k) {
            const double v = detail::read_double(rSource.Points(), i * dim + k);
            if (!std::isfinite(v))
                throw std::invalid_argument("meshio++: register_rigid: source point " +
                                            std::to_string(i) + " has a non-finite coordinate");
            src[i][k] = v;
        }
    lo = hi = src[0];
    for (const RegVec& p : src)
        for (std::size_t a = 0; a < 3; ++a) {
            lo[a] = std::min(lo[a], p[a]);
            hi[a] = std::max(hi[a], p[a]);
        }
    const double extent =
        std::sqrt((hi[0] - lo[0]) * (hi[0] - lo[0]) + (hi[1] - lo[1]) * (hi[1] - lo[1]) +
                  (hi[2] - lo[2]) * (hi[2] - lo[2]));
    const std::array<double, 16>& m0 = rOptions.mInitial;
    for (const double v : m0)
        if (!std::isfinite(v))
            throw std::invalid_argument(
                "meshio++: register_rigid: the initial matrix is not finite");
    RegMat rot{m0[0], m0[1], m0[2], m0[4], m0[5], m0[6], m0[8], m0[9], m0[10]};
    RegVec tr{m0[3], m0[7], m0[11]};

    const bool plane = rOptions.mMethod == RegistrationMethod::PointToPlane;
    const RegTarget target(rTarget, plane);
    RegistrationResult result;
    for (std::int64_t it = 1; it <= rOptions.mMaxIterations; ++it) {
        const RegPairs pairs = reg_pairs(target, src, rot, tr, rOptions.mMaxDistance, plane);
        if (pairs.mP.empty())
            throw std::invalid_argument(
                "meshio++: register_rigid: no source point lies within max_distance of the "
                "target");
        if (it == 1)
            result.mInitialRms = std::sqrt(pairs.mSumSq / static_cast<double>(pairs.mP.size()));
        RegMat d_rot;
        RegVec d_tr;
        if (plane)
            reg_point_to_plane(pairs, d_rot, d_tr);
        else
            reg_point_to_point(pairs, d_rot, d_tr);
        // The step's motion: its rotation, and how far it moves the source's
        // centre (not the origin, which may lie far away).
        const RegVec centre = reg_apply(
            rot, tr, {0.5 * (lo[0] + hi[0]), 0.5 * (lo[1] + hi[1]), 0.5 * (lo[2] + hi[2])});
        const RegVec stepped = reg_apply(d_rot, d_tr, centre);
        const RegVec shift{stepped[0] - centre[0], stepped[1] - centre[1], stepped[2] - centre[2]};
        const double moved = std::sqrt(reg_dot(shift, shift));
        rot = reg_mul(d_rot, rot);
        tr = reg_apply(d_rot, d_tr, tr);
        result.mIterations = it;
        if (reg_angle(d_rot) <= rOptions.mTolerance &&
            moved <= rOptions.mTolerance * std::max(extent, 1e-300)) {
            result.mConverged = true;
            break;
        }
    }
    const RegPairs last = reg_pairs(target, src, rot, tr, rOptions.mMaxDistance, plane);
    result.mNumPairs = static_cast<std::int64_t>(last.mP.size());
    result.mRms =
        last.mP.empty() ? 0.0 : std::sqrt(last.mSumSq / static_cast<double>(last.mP.size()));
    result.mMatrix = {rot[0], rot[1], rot[2], tr[0], rot[3], rot[4], rot[5], tr[1],
                      rot[6], rot[7], rot[8], tr[2], 0.0,    0.0,    0.0,    1.0};
    return result;
}

}  // namespace meshioplusplus

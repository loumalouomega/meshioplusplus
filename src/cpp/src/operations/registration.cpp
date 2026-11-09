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
#include "../detail/bounding_box_utils.hpp"
#include "../detail/jacobi_eigen.hpp"
#include "../detail/lsq_fit.hpp"
#include "../detail/point_knn.hpp"
#include "meshioplusplus/detail/surface_distance.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/types.hpp"
#include "meshioplusplus/detail/geometry.hpp"

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

/// Eigen-decomposition of a symmetric N x N matrix (row-major) by cyclic
/// Jacobi: eigenvalues in `rW`, eigenvectors as the columns of `rV`.
template <std::size_t N>
void reg_jacobi(std::array<double, N * N> A, std::array<double, N>& rW,
                std::array<double, N * N>& rV) {
    detail::jacobi_eigen_sym(A.data(), N, rW.data(), rV.data(), 1e-30);
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
    const double angle = std::sqrt(detail::vec3_dot(rW, rW));
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
    RegTarget(const Mesh& rTarget, bool NeedNormals, const std::string& rRegion = "") {
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
        if (!rRegion.empty() && !surface)
            throw std::invalid_argument(
                "meshio++: register_rigid: target_region needs a target with surface cells");
        if (surface) {
            mSoup = detail::build_triangle_soup(rTarget, rRegion);
            if (mSoup.NumTriangles() == 0)
                throw std::invalid_argument(
                    "meshio++: register_rigid: the target region has no surface cells");
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
                    const double len = std::sqrt(detail::vec3_dot(c, c));
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
        const detail::BoundingBox3 box = detail::BoundingBoxUtils::FromPoints(mPoints);
        const double area = detail::BoundingBoxUtils::HalfSurfaceArea(box);
        const double extent = detail::BoundingBoxUtils::MaxExtent(box);
        double cell = area > 0.0 ? std::sqrt(4.0 * area / static_cast<double>(n)) : extent;
        cell = std::max(cell, std::max(extent, 1.0) * 1e-9);
        mKnn = detail::PointKnn(mPoints, cell);
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
            const std::int64_t j = mKnn.Nearest(rMoved[i]);
            rQ[i] = mPoints[static_cast<std::size_t>(j)];
            if (NeedNormals)
                rN[i] = mPointNormal[static_cast<std::size_t>(j)];
        });
    }

    bool IsSurface() const { return mSurface; }
    const detail::TriangleSoup& Soup() const { return mSoup; }
    const std::vector<RegVec>& Points() const { return mPoints; }
    bool HasPointNormals() const { return !mPointNormal.empty(); }

private:
    std::vector<RegVec> mPoints;
    std::vector<RegVec> mPointNormal;
    std::vector<RegVec> mCellNormal;
    detail::TriangleSoup mSoup;
    detail::DistanceQuery mQuery;
    detail::PointKnn mKnn;
    bool mSurface = false;
};

/// The pairs under the transform (R, t), within `MaxDistance` when positive
/// and within `Trim` times the median pair distance when `Trim` is positive,
/// with finite normals when they are needed.
RegPairs reg_pairs(const RegTarget& rTarget, const std::vector<RegVec>& rSource, const RegMat& rR,
                   const RegVec& rT, double MaxDistance, bool NeedNormals, double Trim = 0.0) {
    std::vector<RegVec> moved(rSource.size());
    parallel_for(rSource.size(), [&](std::size_t i) { moved[i] = reg_apply(rR, rT, rSource[i]); });
    std::vector<RegVec> q, nrm;
    rTarget.Match(moved, q, nrm, NeedNormals);
    double cut2 =
        MaxDistance > 0.0 ? MaxDistance * MaxDistance : std::numeric_limits<double>::infinity();
    if (Trim > 0.0 && !moved.empty()) {
        std::vector<double> d2s(moved.size());
        for (std::size_t i = 0; i < moved.size(); ++i) {
            const RegVec d{q[i][0] - moved[i][0], q[i][1] - moved[i][1], q[i][2] - moved[i][2]};
            d2s[i] = detail::vec3_dot(d, d);
        }
        const std::size_t mid = d2s.size() / 2;
        std::nth_element(d2s.begin(), d2s.begin() + static_cast<std::ptrdiff_t>(mid), d2s.end());
        const double med = std::sqrt(d2s[mid]);
        // A floor keeps the exact pairs of a perfect fit (median zero).
        const double cut = std::max(Trim * med, 1e-300);
        cut2 = std::min(cut2, cut * cut);
    }
    RegPairs out;
    for (std::size_t i = 0; i < moved.size(); ++i) {
        const RegVec d{q[i][0] - moved[i][0], q[i][1] - moved[i][1], q[i][2] - moved[i][2]};
        const double d2 = detail::vec3_dot(d, d);
        if (d2 > cut2)
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

/// Horn's closed-form best rotation and translation for point pairs, each
/// weighted by `pWeights[i]` when given (unweighted otherwise).
void reg_point_to_point(const RegPairs& rPairs, RegMat& rR, RegVec& rT,
                        const std::vector<double>* pWeights = nullptr) {
    double m = static_cast<double>(rPairs.mP.size());
    RegVec pc{0, 0, 0}, qc{0, 0, 0};
    if (pWeights != nullptr) {
        m = 0.0;
        for (std::size_t i = 0; i < rPairs.mP.size(); ++i) {
            const double w = (*pWeights)[i];
            m += w;
            for (std::size_t a = 0; a < 3; ++a) {
                pc[a] += w * rPairs.mP[i][a];
                qc[a] += w * rPairs.mQ[i][a];
            }
        }
    } else
        for (std::size_t i = 0; i < rPairs.mP.size(); ++i)
            for (std::size_t a = 0; a < 3; ++a) {
                pc[a] += rPairs.mP[i][a];
                qc[a] += rPairs.mQ[i][a];
            }
    for (std::size_t a = 0; a < 3; ++a) {
        pc[a] /= m;
        qc[a] /= m;
    }
    RegMat s{};  // S_ab = sum w (p - pc)_a (q - qc)_b
    for (std::size_t i = 0; i < rPairs.mP.size(); ++i) {
        const double w = pWeights != nullptr ? (*pWeights)[i] : 1.0;
        for (std::size_t a = 0; a < 3; ++a)
            for (std::size_t b = 0; b < 3; ++b) {
                const double t = (rPairs.mP[i][a] - pc[a]) * (rPairs.mQ[i][b] - qc[b]);
                s[3 * a + b] += pWeights != nullptr ? w * t : t;
            }
    }
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
        scale += detail::vec3_dot(d, d);
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
        const RegVec x = detail::vec3_cross(p, nn);
        const std::array<double, 6> row{x[0], x[1], x[2], nn[0], nn[1], nn[2]};
        const double b = detail::vec3_dot({q[0] - p[0], q[1] - p[1], q[2] - p[2]}, nn);
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

/// The ICP loop's settings, shared by `register_rigid`'s stages and
/// `find_symmetry_plane`.
struct RegIcpSettings {
    bool mPlane = false;
    std::int64_t mMaxIterations = 50;
    double mTolerance = 1e-10;
    double mMaxDistance = 0.0;
    double mTrim = 0.0;
};

/// What one ICP run ended with.
struct RegIcpOutcome {
    double mInitialRms = 0.0;
    double mRms = 0.0;
    std::int64_t mIterations = 0;
    std::int64_t mNumPairs = 0;
    bool mConverged = false;
    RegPairs mLast;
};

/// ICP from (rRot, rTr), updated in place; `rLo`/`rHi` bound the source.
RegIcpOutcome reg_icp(const RegTarget& rTarget, const std::vector<RegVec>& rSrc,
                      const RegIcpSettings& rSet, RegMat& rRot, RegVec& rTr) {
    const detail::BoundingBox3 box = detail::BoundingBoxUtils::FromPoints(rSrc);
    const RegVec lo = box.mMin, hi = box.mMax;
    const double extent = detail::BoundingBoxUtils::Diagonal(box);
    RegIcpOutcome out;
    for (std::int64_t it = 1; it <= rSet.mMaxIterations; ++it) {
        const RegPairs pairs =
            reg_pairs(rTarget, rSrc, rRot, rTr, rSet.mMaxDistance, rSet.mPlane, rSet.mTrim);
        if (pairs.mP.empty())
            throw std::invalid_argument(
                "meshio++: register_rigid: no source point lies within max_distance of the "
                "target");
        if (it == 1)
            out.mInitialRms = std::sqrt(pairs.mSumSq / static_cast<double>(pairs.mP.size()));
        RegMat d_rot;
        RegVec d_tr;
        if (rSet.mPlane)
            reg_point_to_plane(pairs, d_rot, d_tr);
        else
            reg_point_to_point(pairs, d_rot, d_tr);
        // The step's motion: its rotation, and how far it moves the source's
        // centre (not the origin, which may lie far away).
        const RegVec centre = reg_apply(
            rRot, rTr, {0.5 * (lo[0] + hi[0]), 0.5 * (lo[1] + hi[1]), 0.5 * (lo[2] + hi[2])});
        const RegVec stepped = reg_apply(d_rot, d_tr, centre);
        const RegVec shift{stepped[0] - centre[0], stepped[1] - centre[1], stepped[2] - centre[2]};
        const double moved = std::sqrt(detail::vec3_dot(shift, shift));
        rRot = reg_mul(d_rot, rRot);
        rTr = reg_apply(d_rot, d_tr, rTr);
        out.mIterations = it;
        if (reg_angle(d_rot) <= rSet.mTolerance &&
            moved <= rSet.mTolerance * std::max(extent, 1e-300)) {
            out.mConverged = true;
            break;
        }
    }
    out.mLast = reg_pairs(rTarget, rSrc, rRot, rTr, rSet.mMaxDistance, rSet.mPlane, rSet.mTrim);
    out.mNumPairs = static_cast<std::int64_t>(out.mLast.mP.size());
    out.mRms = out.mLast.mP.empty()
                   ? 0.0
                   : std::sqrt(out.mLast.mSumSq / static_cast<double>(out.mLast.mP.size()));
    return out;
}

/// A principal frame: the centroid and the axes (columns of `mAxes`, by
/// decreasing spread, right-handed).
struct RegFrame {
    RegVec mCentre{0, 0, 0};
    RegMat mAxes{1, 0, 0, 0, 1, 0, 0, 0, 1};
};

/// The frame of second moments `rM` (about the origin, total weight `W`,
/// first moment `rS`).
RegFrame reg_frame_from_moments(const RegMat& rM, const RegVec& rS, double W) {
    RegFrame f;
    if (!(W > 0.0))
        return f;
    f.mCentre = {rS[0] / W, rS[1] / W, rS[2] / W};
    std::array<double, 9> cov{};
    for (std::size_t a = 0; a < 3; ++a)
        for (std::size_t b = 0; b < 3; ++b)
            cov[3 * a + b] = rM[3 * a + b] / W - f.mCentre[a] * f.mCentre[b];
    std::array<double, 3> w{};
    std::array<double, 9> v{};
    detail::lsq_sym3_eigen(cov.data(), w.data(), v.data());  // ascending, rows
    // Columns by decreasing spread: rows 2, 1, 0 of `v`.
    for (std::size_t k = 0; k < 3; ++k)
        for (std::size_t a = 0; a < 3; ++a)
            f.mAxes[3 * a + k] = v[3 * (2 - k) + a];
    const RegVec e0{f.mAxes[0], f.mAxes[3], f.mAxes[6]}, e1{f.mAxes[1], f.mAxes[4], f.mAxes[7]};
    const RegVec e2 = detail::vec3_cross(e0, e1);
    f.mAxes[2] = e2[0];
    f.mAxes[5] = e2[1];
    f.mAxes[8] = e2[2];
    return f;
}

/// The principal frame of a triangle soup (area-weighted, exact per
/// triangle) or, with no triangles, of the points.
RegFrame reg_frame(const detail::TriangleSoup* pSoup, const std::vector<RegVec>& rPoints) {
    RegMat m{};
    RegVec s{0, 0, 0};
    double w = 0.0;
    if (pSoup != nullptr && pSoup->NumTriangles() > 0) {
        for (std::size_t t = 0; t < pSoup->NumTriangles(); ++t) {
            const RegVec& a = pSoup->mCorners[3 * t];
            const RegVec& b = pSoup->mCorners[3 * t + 1];
            const RegVec& c = pSoup->mCorners[3 * t + 2];
            const RegVec ab{b[0] - a[0], b[1] - a[1], b[2] - a[2]};
            const RegVec ac{c[0] - a[0], c[1] - a[1], c[2] - a[2]};
            const RegVec x = detail::vec3_cross(ab, ac);
            const double area = 0.5 * std::sqrt(detail::vec3_dot(x, x));
            const RegVec sum{a[0] + b[0] + c[0], a[1] + b[1] + c[1], a[2] + b[2] + c[2]};
            // The integral of x x^T over a triangle: A/12 (sum v v^T + s s^T).
            for (std::size_t i = 0; i < 3; ++i)
                for (std::size_t j = 0; j < 3; ++j)
                    m[3 * i + j] +=
                        area / 12.0 * (a[i] * a[j] + b[i] * b[j] + c[i] * c[j] + sum[i] * sum[j]);
            for (std::size_t i = 0; i < 3; ++i)
                s[i] += area * sum[i] / 3.0;
            w += area;
        }
        if (w > 0.0)
            return reg_frame_from_moments(m, s, w);
        m = {};
        s = {0, 0, 0};
    }
    for (const RegVec& p : rPoints) {
        for (std::size_t i = 0; i < 3; ++i) {
            s[i] += p[i];
            for (std::size_t j = 0; j < 3; ++j)
                m[3 * i + j] += p[i] * p[j];
        }
        w += 1.0;
    }
    return reg_frame_from_moments(m, s, w);
}

/// The 24 proper signed permutation matrices, in a fixed order (permutations
/// lexicographic, then the sign patterns in binary order).
std::vector<RegMat> reg_proper_axis_maps() {
    std::vector<RegMat> out;
    std::array<int, 3> perm{0, 1, 2};
    do {
        for (int signs = 0; signs < 8; ++signs) {
            RegMat p{};
            for (int r = 0; r < 3; ++r)
                p[static_cast<std::size_t>(3 * r + perm[static_cast<std::size_t>(r)])] =
                    (signs >> r & 1) != 0 ? -1.0 : 1.0;
            const double det = p[0] * (p[4] * p[8] - p[5] * p[7]) -
                               p[1] * (p[3] * p[8] - p[5] * p[6]) +
                               p[2] * (p[3] * p[7] - p[4] * p[6]);
            if (det > 0.0)
                out.push_back(p);
        }
    } while (std::next_permutation(perm.begin(), perm.end()));
    return out;
}

RegMat reg_transpose(const RegMat& rA) {
    return {rA[0], rA[3], rA[6], rA[1], rA[4], rA[7], rA[2], rA[5], rA[8]};
}

/// The RMS of the closest `Keep` share of the pairs' distances.
double reg_trimmed_rms(const RegPairs& rPairs, double Keep) {
    std::vector<double> d2(rPairs.mP.size());
    for (std::size_t i = 0; i < d2.size(); ++i) {
        const RegVec d{rPairs.mQ[i][0] - rPairs.mP[i][0], rPairs.mQ[i][1] - rPairs.mP[i][1],
                       rPairs.mQ[i][2] - rPairs.mP[i][2]};
        d2[i] = detail::vec3_dot(d, d);
    }
    if (d2.empty())
        return std::numeric_limits<double>::infinity();
    std::sort(d2.begin(), d2.end());
    const std::size_t k =
        std::max<std::size_t>(1, static_cast<std::size_t>(Keep * static_cast<double>(d2.size())));
    double sum = 0.0;
    for (std::size_t i = 0; i < k; ++i)
        sum += d2[i];
    return std::sqrt(sum / static_cast<double>(k));
}

/// Every `Stride`-th point, so at most `Max` of them.
std::vector<RegVec> reg_stride(const std::vector<RegVec>& rPts, std::size_t Max) {
    if (rPts.size() <= Max || Max == 0)
        return rPts;
    const std::size_t stride = (rPts.size() + Max - 1) / Max;
    std::vector<RegVec> out;
    out.reserve(Max);
    for (std::size_t i = 0; i < rPts.size(); i += stride)
        out.push_back(rPts[i]);
    return out;
}

/// The finite points of @p rMesh, `pWhat` naming it in errors.
std::vector<RegVec> reg_points(const Mesh& rMesh, const char* pWhat) {
    const std::size_t n = rMesh.NumPoints();
    if (n == 0)
        throw std::invalid_argument(std::string("meshio++: ") + pWhat + " has no points");
    const std::size_t dim = rMesh.PointDim();
    std::vector<RegVec> out(n, RegVec{0.0, 0.0, 0.0});
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t k = 0; k < dim && k < 3; ++k) {
            const double v = detail::read_double(rMesh.Points(), i * dim + k);
            if (!std::isfinite(v))
                throw std::invalid_argument(std::string("meshio++: ") + pWhat + " point " +
                                            std::to_string(i) + " has a non-finite coordinate");
            out[i][k] = v;
        }
    return out;
}

bool reg_has_surface(const Mesh& rMesh) {
    const auto& dims = topological_dimension();
    for (const auto cb : rMesh.CellRange()) {
        const auto it = dims.find(cb.Type());
        if (it != dims.end() && it->second == 2 && cb.NumCells() > 0)
            return true;
    }
    return false;
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

RegistrationInitial registration_initial_from_name(const std::string& rName) {
    if (rName == "matrix")
        return RegistrationInitial::Matrix;
    if (rName == "auto")
        return RegistrationInitial::Auto;
    throw std::invalid_argument("meshio++: register_rigid: unknown initial '" + rName +
                                "' (expected 'matrix' or 'auto')");
}

RegistrationResult register_rigid(const Mesh& rSource, const Mesh& rTarget,
                                  const RegistrationOptions& rOptions) {
    if (rOptions.mMaxIterations < 1)
        throw std::invalid_argument("meshio++: register_rigid: max_iterations must be at least 1");
    if (!(rOptions.mTolerance >= 0.0) || !(rOptions.mMaxDistance >= 0.0))
        throw std::invalid_argument(
            "meshio++: register_rigid: tolerance and max_distance must not be negative");
    if (!(rOptions.mAdaptiveTrim >= 0.0) || !std::isfinite(rOptions.mAdaptiveTrim))
        throw std::invalid_argument(
            "meshio++: register_rigid: adaptive_trim must be finite and not negative");
    if (!(rOptions.mSymmetricTolerance >= 0.0))
        throw std::invalid_argument(
            "meshio++: register_rigid: symmetric_tolerance must not be negative");
    const std::size_t n = rSource.NumPoints();
    if (n == 0)
        throw std::invalid_argument("meshio++: register_rigid: the source has no points");
    std::vector<RegVec> src = reg_points(rSource, "register_rigid: source");
    detail::TriangleSoup src_soup;
    const bool src_surface = reg_has_surface(rSource);
    if (!rOptions.mSourceRegion.empty()) {
        if (!src_surface)
            throw std::invalid_argument(
                "meshio++: register_rigid: source_region needs a source with surface cells");
        src_soup = detail::build_triangle_soup(rSource, rOptions.mSourceRegion);
        std::vector<std::uint8_t> used(n, 0);
        for (const auto& tri : src_soup.mVertices)
            for (const std::int64_t v : tri)
                used[static_cast<std::size_t>(v)] = 1;
        std::vector<RegVec> kept;
        for (std::size_t i = 0; i < n; ++i)
            if (used[i] != 0)
                kept.push_back(src[i]);
        if (kept.empty())
            throw std::invalid_argument(
                "meshio++: register_rigid: the source region has no surface cells");
        src = std::move(kept);
    } else if (src_surface && rOptions.mInitialMode == RegistrationInitial::Auto)
        src_soup = detail::build_triangle_soup(rSource, "");
    const std::array<double, 16>& m0 = rOptions.mInitial;
    for (const double v : m0)
        if (!std::isfinite(v))
            throw std::invalid_argument(
                "meshio++: register_rigid: the initial matrix is not finite");
    RegMat rot{m0[0], m0[1], m0[2], m0[4], m0[5], m0[6], m0[8], m0[9], m0[10]};
    RegVec tr{m0[3], m0[7], m0[11]};

    const bool plane = rOptions.mMethod == RegistrationMethod::PointToPlane;
    const RegTarget target(rTarget, plane, rOptions.mTargetRegion);
    RegistrationResult result;
    if (rOptions.mInitialMode == RegistrationInitial::Auto) {
        const RegFrame fs = reg_frame(&src_soup, src);
        const RegFrame ft =
            reg_frame(target.IsSurface() ? &target.Soup() : nullptr, target.Points());
        const std::vector<RegVec> coarse_src = reg_stride(src, 1000);
        RegIcpSettings coarse;
        coarse.mPlane = false;
        coarse.mMaxIterations = 30;
        coarse.mTolerance = 1e-6;
        coarse.mMaxDistance = rOptions.mMaxDistance;
        coarse.mTrim = rOptions.mAdaptiveTrim;
        const std::vector<RegMat> maps = reg_proper_axis_maps();
        std::vector<RegMat> rots(maps.size());
        std::vector<RegVec> trs(maps.size());
        std::size_t best = 0;
        for (std::size_t k = 0; k < maps.size(); ++k) {
            RegMat r = reg_mul(ft.mAxes, reg_mul(maps[k], reg_transpose(fs.mAxes)));
            const RegVec rc = reg_apply(r, {0, 0, 0}, fs.mCentre);
            RegVec t{ft.mCentre[0] - rc[0], ft.mCentre[1] - rc[1], ft.mCentre[2] - rc[2]};
            const RegIcpOutcome o = reg_icp(target, coarse_src, coarse, r, t);
            rots[k] = r;
            trs[k] = t;
            result.mCandidateRms.push_back(reg_trimmed_rms(o.mLast, 0.9));
            if (result.mCandidateRms[k] < result.mCandidateRms[best])
                best = k;
        }
        result.mChosenCandidate = static_cast<std::int64_t>(best);
        if (rOptions.mSymmetricCheck) {
            const double rb = result.mCandidateRms[best];
            const detail::BoundingBox3 sbox = detail::BoundingBoxUtils::FromPoints(src);
            const double floor = 1e-12 * detail::BoundingBoxUtils::Diagonal(sbox);
            for (std::size_t k = 0; k < maps.size(); ++k) {
                if (k == best)
                    continue;
                const double apart = reg_angle(reg_mul(rots[k], reg_transpose(rots[best])));
                if (apart > 5.0 * 3.14159265358979323846 / 180.0 &&
                    result.mCandidateRms[k] <= (1.0 + rOptions.mSymmetricTolerance) * rb + floor) {
                    result.mWarnings.push_back(
                        "candidate " + std::to_string(k) +
                        " fits within symmetric_tolerance of "
                        "the chosen candidate " +
                        std::to_string(best) + " at a pose " +
                        std::to_string(apart * 180.0 / 3.14159265358979323846) +
                        " degrees away: the part may be symmetric, and its pose undecided");
                    break;
                }
            }
        }
        rot = rots[best];
        tr = trs[best];
    }
    RegIcpSettings fine;
    fine.mPlane = plane;
    fine.mMaxIterations = rOptions.mMaxIterations;
    fine.mTolerance = rOptions.mTolerance;
    fine.mMaxDistance = rOptions.mMaxDistance;
    fine.mTrim = rOptions.mAdaptiveTrim;
    const RegIcpOutcome o = reg_icp(target, src, fine, rot, tr);
    result.mInitialRms = o.mInitialRms;
    result.mRms = o.mRms;
    result.mIterations = o.mIterations;
    result.mConverged = o.mConverged;
    result.mNumPairs = o.mNumPairs;
    result.mNumRejected = static_cast<std::int64_t>(src.size()) - o.mNumPairs;
    result.mMatrix = {rot[0], rot[1], rot[2], tr[0], rot[3], rot[4], rot[5], tr[1],
                      rot[6], rot[7], rot[8], tr[2], 0.0,    0.0,    0.0,    1.0};
    return result;
}

LandmarkResult register_landmarks(const std::vector<std::array<double, 3>>& rSource,
                                  const std::vector<std::array<double, 3>>& rTarget,
                                  const std::vector<double>& rWeights) {
    const char* prefix = "meshio++: register_landmarks: ";
    if (rSource.size() != rTarget.size())
        throw std::invalid_argument(std::string(prefix) + "the source has " +
                                    std::to_string(rSource.size()) + " points, the target " +
                                    std::to_string(rTarget.size()));
    if (rSource.size() < 3)
        throw std::invalid_argument(std::string(prefix) + "needs at least three pairs, got " +
                                    std::to_string(rSource.size()));
    if (!rWeights.empty() && rWeights.size() != rSource.size())
        throw std::invalid_argument(std::string(prefix) + "needs one weight per pair");
    RegPairs pairs;
    std::vector<double> w(rSource.size(), 1.0);
    std::size_t positive = 0;
    for (std::size_t i = 0; i < rSource.size(); ++i) {
        for (std::size_t a = 0; a < 3; ++a)
            if (!std::isfinite(rSource[i][a]) || !std::isfinite(rTarget[i][a]))
                throw std::invalid_argument(std::string(prefix) + "pair " + std::to_string(i) +
                                            " has a non-finite coordinate");
        if (!rWeights.empty()) {
            if (!std::isfinite(rWeights[i]) || rWeights[i] < 0.0)
                throw std::invalid_argument(std::string(prefix) + "weight " + std::to_string(i) +
                                            " must be finite and not negative");
            w[i] = rWeights[i];
        }
        positive += w[i] > 0.0 ? 1 : 0;
        pairs.mP.push_back(rSource[i]);
        pairs.mQ.push_back(rTarget[i]);
    }
    if (positive < 3)
        throw std::invalid_argument(std::string(prefix) +
                                    "needs at least three pairs of positive weight");
    // Collinear (or coincident) sources leave the rotation about their line free.
    {
        RegMat m{};
        RegVec s{0, 0, 0};
        double tw = 0.0;
        for (std::size_t i = 0; i < pairs.mP.size(); ++i) {
            if (!(w[i] > 0.0))
                continue;
            for (std::size_t a = 0; a < 3; ++a) {
                s[a] += w[i] * pairs.mP[i][a];
                for (std::size_t b = 0; b < 3; ++b)
                    m[3 * a + b] += w[i] * pairs.mP[i][a] * pairs.mP[i][b];
            }
            tw += w[i];
        }
        std::array<double, 9> cov{};
        for (std::size_t a = 0; a < 3; ++a)
            for (std::size_t b = 0; b < 3; ++b)
                cov[3 * a + b] = m[3 * a + b] / tw - s[a] * s[b] / (tw * tw);
        std::array<double, 3> ev{};
        std::array<double, 9> vec{};
        detail::lsq_sym3_eigen(cov.data(), ev.data(), vec.data());
        if (!(ev[1] > 1e-20 * std::max(ev[2], 1e-300)) || !(ev[2] > 0.0))
            throw std::invalid_argument(std::string(prefix) +
                                        "the source points are collinear: the rotation about "
                                        "their line is undetermined");
    }
    RegMat r;
    RegVec t;
    reg_point_to_point(pairs, r, t, &w);
    LandmarkResult out;
    double sum = 0.0, tw = 0.0;
    for (std::size_t i = 0; i < pairs.mP.size(); ++i) {
        const RegVec p = reg_apply(r, t, pairs.mP[i]);
        const RegVec d{p[0] - pairs.mQ[i][0], p[1] - pairs.mQ[i][1], p[2] - pairs.mQ[i][2]};
        const double d2 = detail::vec3_dot(d, d);
        out.mResiduals.push_back(std::sqrt(d2));
        sum += w[i] * d2;
        tw += w[i];
    }
    out.mRms = std::sqrt(sum / tw);
    out.mMatrix = {r[0], r[1], r[2], t[0], r[3], r[4], r[5], t[1],
                   r[6], r[7], r[8], t[2], 0.0,  0.0,  0.0,  1.0};
    return out;
}

SymmetryPlaneResult find_symmetry_plane(const Mesh& rMesh, const SymmetryPlaneOptions& rOptions) {
    const char* prefix = "meshio++: find_symmetry_plane: ";
    if (!(rOptions.mTolerance >= 0.0) || !std::isfinite(rOptions.mTolerance))
        throw std::invalid_argument(std::string(prefix) +
                                    "tolerance must be finite and not negative");
    if (!(rOptions.mMinShare >= 0.0 && rOptions.mMinShare <= 1.0))
        throw std::invalid_argument(std::string(prefix) + "min_share must lie in [0, 1]");
    if (rOptions.mMaxSamples < 1 || rOptions.mMaxIterations < 1)
        throw std::invalid_argument(std::string(prefix) +
                                    "max_samples and max_iterations must be at least 1");
    for (std::size_t a = 0; a < 3; ++a)
        if (!std::isfinite(rOptions.mSeedOrigin[a]) || !std::isfinite(rOptions.mSeedNormal[a]))
            throw std::invalid_argument(std::string(prefix) + "the seed plane is not finite");
    const std::vector<RegVec> pts = reg_points(rMesh, "find_symmetry_plane: the input");
    const bool plane = reg_has_surface(rMesh);
    const RegTarget target(rMesh, plane, "");
    const RegFrame frame = reg_frame(plane ? &target.Soup() : nullptr, pts);
    const detail::BoundingBox3 box = detail::BoundingBoxUtils::FromPoints(pts);
    const double diag = detail::BoundingBoxUtils::Diagonal(box);
    if (!(diag > 0.0))
        throw std::invalid_argument(std::string(prefix) + "the input has no extent");
    SymmetryPlaneResult result;
    result.mTolerance = rOptions.mTolerance > 0.0 ? rOptions.mTolerance : 0.005 * diag;

    struct Plane {
        RegVec mO, mN;
    };
    std::vector<Plane> cands;
    for (std::size_t k = 0; k < 3; ++k)
        cands.push_back({frame.mCentre, {frame.mAxes[k], frame.mAxes[3 + k], frame.mAxes[6 + k]}});
    const double seed_len = std::sqrt(detail::vec3_dot(rOptions.mSeedNormal, rOptions.mSeedNormal));
    if (seed_len > 0.0)
        cands.push_back({rOptions.mSeedOrigin,
                         {rOptions.mSeedNormal[0] / seed_len, rOptions.mSeedNormal[1] / seed_len,
                          rOptions.mSeedNormal[2] / seed_len}});
    const auto reflect = [](const std::vector<RegVec>& rIn, const Plane& rPl) {
        std::vector<RegVec> out(rIn.size());
        for (std::size_t i = 0; i < rIn.size(); ++i) {
            const RegVec& p = rIn[i];
            const double d = (p[0] - rPl.mO[0]) * rPl.mN[0] + (p[1] - rPl.mO[1]) * rPl.mN[1] +
                             (p[2] - rPl.mO[2]) * rPl.mN[2];
            out[i] = {p[0] - 2.0 * d * rPl.mN[0], p[1] - 2.0 * d * rPl.mN[1],
                      p[2] - 2.0 * d * rPl.mN[2]};
        }
        return out;
    };
    // The plane of the improper isometry `x -> R H x + t` (H the reflection
    // across `rPl`): the eigenvector of its linear part's -1 eigenvalue.
    const auto extract = [&](const Plane& rPl, const RegMat& rR, const RegVec& rT) {
        RegMat h{};
        for (std::size_t a = 0; a < 3; ++a)
            for (std::size_t b = 0; b < 3; ++b)
                h[3 * a + b] = (a == b ? 1.0 : 0.0) - 2.0 * rPl.mN[a] * rPl.mN[b];
        const RegMat l = reg_mul(rR, h);
        const double c = 2.0 * detail::vec3_dot(rPl.mO, rPl.mN);
        const RegVec rn = reg_apply(rR, {0, 0, 0}, rPl.mN);
        const RegVec b{rT[0] + c * rn[0], rT[1] + c * rn[1], rT[2] + c * rn[2]};
        std::array<double, 9> sym{};
        for (std::size_t a = 0; a < 3; ++a)
            for (std::size_t k = 0; k < 3; ++k)
                sym[3 * a + k] = 0.5 * (l[3 * a + k] + l[3 * k + a]);
        std::array<double, 3> ev{};
        std::array<double, 9> vec{};
        detail::lsq_sym3_eigen(sym.data(), ev.data(), vec.data());
        RegVec nn{vec[0], vec[1], vec[2]};  // smallest eigenvalue, signed positive
        const double d = 0.5 * detail::vec3_dot(b, nn);
        const double off = detail::vec3_dot(frame.mCentre, nn) - d;
        Plane out;
        out.mN = nn;
        out.mO = {frame.mCentre[0] - off * nn[0], frame.mCentre[1] - off * nn[1],
                  frame.mCentre[2] - off * nn[2]};
        return out;
    };
    RegIcpSettings set;
    set.mPlane = plane;
    set.mMaxIterations = rOptions.mMaxIterations;
    set.mTolerance = 1e-6;
    set.mTrim = 3.0;
    const std::vector<RegVec> coarse = reg_stride(pts, 500);
    std::size_t best = 0;
    for (std::size_t k = 0; k < cands.size(); ++k) {
        RegMat r{1, 0, 0, 0, 1, 0, 0, 0, 1};
        RegVec t{0, 0, 0};
        const RegIcpOutcome o = reg_icp(target, reflect(coarse, cands[k]), set, r, t);
        result.mCandidateRms.push_back(reg_trimmed_rms(o.mLast, 0.9));
        if (result.mCandidateRms[k] < result.mCandidateRms[best])
            best = k;
    }
    result.mChosenCandidate = static_cast<std::int64_t>(best);
    // Fine: from the chosen candidate on the larger sample, then twice more
    // from the plane it gives, so the recovered plane is its own fixed point.
    const std::vector<RegVec> fine =
        reg_stride(pts, static_cast<std::size_t>(rOptions.mMaxSamples));
    set.mTolerance = 1e-12;
    Plane cur = cands[best];
    for (int round = 0; round < 3; ++round) {
        RegMat r{1, 0, 0, 0, 1, 0, 0, 0, 1};
        RegVec t{0, 0, 0};
        reg_icp(target, reflect(fine, cur), set, r, t);
        cur = extract(cur, r, t);
    }
    // Measure the plane itself: its exact mirror image against the surface.
    const RegPairs meas = reg_pairs(target, reflect(fine, cur), {1, 0, 0, 0, 1, 0, 0, 0, 1},
                                    {0, 0, 0}, 0.0, false, 0.0);
    std::size_t within = 0;
    for (std::size_t i = 0; i < meas.mP.size(); ++i) {
        const RegVec d{meas.mQ[i][0] - meas.mP[i][0], meas.mQ[i][1] - meas.mP[i][1],
                       meas.mQ[i][2] - meas.mP[i][2]};
        within += detail::vec3_dot(d, d) <= result.mTolerance * result.mTolerance ? 1 : 0;
    }
    result.mRms =
        std::sqrt(meas.mSumSq / static_cast<double>(std::max<std::size_t>(1, meas.mP.size())));
    result.mShare =
        static_cast<double>(within) / static_cast<double>(std::max<std::size_t>(1, meas.mP.size()));
    result.mOrigin = cur.mO;
    result.mNormal = cur.mN;
    if (result.mShare < rOptions.mMinShare)
        throw std::invalid_argument(std::string(prefix) +
                                    "no symmetry plane: the best plane's share " +
                                    std::to_string(result.mShare) +
                                    " of samples within tolerance is below "
                                    "min_share " +
                                    std::to_string(rOptions.mMinShare));
    return result;
}

}  // namespace meshioplusplus

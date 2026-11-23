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
//                   SPDX-License-Identifier: MIT
//                   meshio++ default license: LICENSE
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//

// Streamlines of a point vector field. See detail/streamlines.hpp.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

// Project includes
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/detail/geometry.hpp"

// Project includes (private, not installed)
#include "streamlines.hpp"
#include "math_utils.hpp"

namespace meshioplusplus {
namespace detail {
namespace {

// How a cell of a given type name is cut into simplices, as corner indices into
// its connectivity row; `Corners` is how many leading nodes it reads.
struct SlCut {
    int mDim = 0;
    std::size_t mCorners = 0;
    std::vector<std::array<int, 4>> mSimplices;
};

bool sl_cut(const std::string& rType, SlCut& rCut) {
    auto starts = [&](const char* pPrefix) { return rType.rfind(pPrefix, 0) == 0; };
    if (starts("triangle")) {
        rCut = {2, 3, {{0, 1, 2, -1}}};
    } else if (starts("quad")) {
        rCut = {2, 4, {{0, 1, 2, -1}, {0, 2, 3, -1}}};
    } else if (starts("tetra")) {
        rCut = {3, 4, {{0, 1, 2, 3}}};
    } else if (starts("hexahedron")) {
        // Six tetrahedra around the diagonal 0-6.
        rCut = {
            3,
            8,
            {{0, 1, 2, 6}, {0, 2, 3, 6}, {0, 3, 7, 6}, {0, 7, 4, 6}, {0, 4, 5, 6}, {0, 5, 1, 6}}};
    } else if (starts("wedge")) {
        rCut = {3, 6, {{0, 1, 2, 3}, {1, 2, 5, 3}, {1, 5, 4, 3}}};
    } else if (starts("pyramid")) {
        rCut = {3, 5, {{0, 1, 2, 4}, {0, 2, 3, 4}}};
    } else {
        return false;
    }
    return true;
}

struct SlSimplex {
    std::array<std::int64_t, 4> mNode;
};

// The simplices of the mesh's highest supported dimension, and the (sorted,
// unique) names of the cell types that had to be left out.
struct SlSoup {
    int mDim = 0;
    std::vector<SlSimplex> mSimplices;
    std::vector<std::string> mSeen;
};

SlSoup sl_collect(const Mesh& rMesh) {
    SlSoup out;
    std::vector<SlSimplex> by_dim[2];
    for (const auto cb : rMesh.CellRange()) {
        const std::size_t n = cb.NumCells();
        const std::string& type = cb.Type();
        if (n > 0 && std::find(out.mSeen.begin(), out.mSeen.end(), type) == out.mSeen.end())
            out.mSeen.push_back(type);
        SlCut cut;
        if (n == 0 || cb.IsRagged() || cb.IsPolyhedron() || !sl_cut(type, cut))
            continue;
        const NDArray& conn = cb.Conn();
        const std::size_t cols = conn.Size() / n;
        if (cols < cut.mCorners)
            continue;
        std::vector<SlSimplex>& dest = by_dim[cut.mDim - 2];
        for (std::size_t r = 0; r < n; ++r)
            for (const std::array<int, 4>& s : cut.mSimplices) {
                SlSimplex simplex;
                simplex.mNode = {-1, -1, -1, -1};
                for (std::size_t k = 0; k < static_cast<std::size_t>(cut.mDim) + 1; ++k)
                    simplex.mNode[k] = read_int(conn, r * cols + static_cast<std::size_t>(s[k]));
                dest.push_back(simplex);
            }
    }
    std::sort(out.mSeen.begin(), out.mSeen.end());
    if (!by_dim[1].empty()) {
        out.mDim = 3;
        out.mSimplices = std::move(by_dim[1]);
    } else if (!by_dim[0].empty()) {
        out.mDim = 2;
        out.mSimplices = std::move(by_dim[0]);
    }
    return out;
}

// A uniform bucket grid over the simplices' boxes, filled in simplex order, and
// the point-in-simplex test and interpolation on top of it.
class SlField {
public:
    SlField(const SlSoup& rSoup, const double* pXyz, const double* pVec, double Diagonal,
            double Tolerance)
        : mrSoup(rSoup), mpXyz(pXyz), mpVec(pVec), mTol(Tolerance), mDim(rSoup.mDim) {
        const std::size_t count = rSoup.mSimplices.size();
        std::vector<Vec3> lo(count);
        std::vector<Vec3> hi(count);
        const std::size_t corners = static_cast<std::size_t>(mDim) + 1;
        for (std::size_t s = 0; s < count; ++s) {
            for (std::size_t c = 0; c < 3; ++c) {
                lo[s][c] = std::numeric_limits<double>::infinity();
                hi[s][c] = -lo[s][c];
            }
            for (std::size_t k = 0; k < corners; ++k) {
                const double* p = point(rSoup.mSimplices[s].mNode[k]);
                for (std::size_t c = 0; c < 3; ++c) {
                    lo[s][c] = std::min(lo[s][c], p[c]);
                    hi[s][c] = std::max(hi[s][c], p[c]);
                }
            }
            for (std::size_t c = 0; c < 3; ++c)
                if (!(lo[s][c] == lo[s][c]) || !(hi[s][c] == hi[s][c]))
                    lo[s][c] = hi[s][c] = 0.0;
        }
        for (std::size_t c = 0; c < 3; ++c) {
            mLo[c] = std::numeric_limits<double>::infinity();
            mHi[c] = -mLo[c];
        }
        for (std::size_t s = 0; s < count; ++s)
            for (std::size_t c = 0; c < 3; ++c) {
                mLo[c] = std::min(mLo[c], lo[s][c]);
                mHi[c] = std::max(mHi[c], hi[s][c]);
            }
        for (std::size_t c = 0; c < 3; ++c) {
            mLo[c] -= mTol;
            mHi[c] += mTol;
        }
        // About two simplices per bucket, over the axes that have any extent.
        const double flat = 1e-6 * Diagonal;
        double volume = 1.0;
        int live = 0;
        for (std::size_t c = 0; c < 3; ++c)
            if (mHi[c] - mLo[c] > flat + 2.0 * mTol) {
                volume *= mHi[c] - mLo[c];
                ++live;
            }
        const double buckets = std::max(1.0, static_cast<double>(count) / 2.0);
        double size = live > 0 ? std::pow(volume / buckets, 1.0 / live) : 1.0;
        for (std::size_t c = 0; c < 3; ++c) {
            const double extent = mHi[c] - mLo[c];
            std::size_t r = 1;
            if (extent > flat + 2.0 * mTol && size > 0.0)
                r = static_cast<std::size_t>(std::clamp(std::ceil(extent / size), 1.0, 256.0));
            mRes[c] = r;
            mCell[c] = extent / static_cast<double>(r);
        }
        mOffset.assign(mRes[0] * mRes[1] * mRes[2] + 1, 0);
        // Two passes, in simplex order: count, then fill.
        for (int pass = 0; pass < 2; ++pass) {
            std::vector<std::size_t> fill(mOffset.size() - 1, 0);
            for (std::size_t s = 0; s < count; ++s) {
                std::size_t a[3];
                std::size_t b[3];
                for (std::size_t c = 0; c < 3; ++c) {
                    a[c] = bucket(c, lo[s][c] - mTol);
                    b[c] = bucket(c, hi[s][c] + mTol);
                }
                for (std::size_t z = a[2]; z <= b[2]; ++z)
                    for (std::size_t y = a[1]; y <= b[1]; ++y)
                        for (std::size_t x = a[0]; x <= b[0]; ++x) {
                            const std::size_t id = x + mRes[0] * (y + mRes[1] * z);
                            if (pass == 0)
                                ++mOffset[id + 1];
                            else
                                mItems[mOffset[id] + fill[id]++] = s;
                        }
            }
            if (pass == 0) {
                for (std::size_t i = 1; i < mOffset.size(); ++i)
                    mOffset[i] += mOffset[i - 1];
                mItems.assign(mOffset.back(), 0);
            }
        }
    }

    /// The field at `pPoint` (a point of the plane, for a surface, is moved onto
    /// the triangle): the interpolated vector, and the simplex that holds it.
    bool Sample(double* pPoint, double* pVector, std::size_t* pSimplex = nullptr) const {
        for (std::size_t c = 0; c < 3; ++c)
            if (!(pPoint[c] > mLo[c]) || !(pPoint[c] < mHi[c]))
                return false;
        std::size_t at[3];
        for (std::size_t c = 0; c < 3; ++c)
            at[c] = bucket(c, pPoint[c]);
        const std::size_t id = at[0] + mRes[0] * (at[1] + mRes[1] * at[2]);
        double best = std::numeric_limits<double>::infinity();
        bool found = false;
        double best_vec[3] = {0, 0, 0};
        double best_point[3] = {0, 0, 0};
        std::size_t best_simplex = 0;
        for (std::size_t i = mOffset[id]; i < mOffset[id + 1]; ++i) {
            const std::size_t s = mItems[i];
            double lambda[4];
            double distance = 0.0;
            Vec3 normal = {0, 0, 0};
            if (!Locate(s, pPoint, lambda, distance, normal))
                continue;
            if (std::fabs(distance) >= best)
                continue;
            best = std::fabs(distance);
            found = true;
            best_simplex = s;
            for (std::size_t c = 0; c < 3; ++c) {
                double v = 0.0;
                for (std::size_t k = 0; k < static_cast<std::size_t>(mDim) + 1; ++k)
                    v += lambda[k] *
                         mpVec[3 * static_cast<std::size_t>(mrSoup.mSimplices[s].mNode[k]) + c];
                best_vec[c] = v;
                best_point[c] = pPoint[c] - (mDim == 2 ? distance * normal[c] : 0.0);
            }
            if (mDim == 2) {
                const double along = vec3_dot({best_vec[0], best_vec[1], best_vec[2]}, normal);
                for (std::size_t c = 0; c < 3; ++c)
                    best_vec[c] -= along * normal[c];
            }
            if (mDim == 3)
                break;  // simplices of a volume tile it: the first hit is the cell
        }
        if (!found)
            return false;
        for (std::size_t c = 0; c < 3; ++c) {
            pVector[c] = best_vec[c];
            pPoint[c] = best_point[c];
        }
        if (pSimplex != nullptr)
            *pSimplex = best_simplex;
        return true;
    }

    const double* point(std::int64_t Id) const { return mpXyz + 3 * static_cast<std::size_t>(Id); }

    double Measure(std::size_t Simplex) const {
        const SlSimplex& s = mrSoup.mSimplices[Simplex];
        if (mDim == 3) {
            const Vec3 a = detail::MathUtils::Sub(point(s.mNode[0]), point(s.mNode[3]));
            const Vec3 b = detail::MathUtils::Sub(point(s.mNode[1]), point(s.mNode[3]));
            const Vec3 c = detail::MathUtils::Sub(point(s.mNode[2]), point(s.mNode[3]));
            return std::fabs(vec3_dot(a, vec3_cross(b, c))) / 6.0;
        }
        const Vec3 a = detail::MathUtils::Sub(point(s.mNode[1]), point(s.mNode[0]));
        const Vec3 b = detail::MathUtils::Sub(point(s.mNode[2]), point(s.mNode[0]));
        const Vec3 n = vec3_cross(a, b);
        return 0.5 * std::sqrt(vec3_dot(n, n));
    }

    Vec3 Centroid(std::size_t Simplex) const {
        const SlSimplex& s = mrSoup.mSimplices[Simplex];
        const std::size_t corners = static_cast<std::size_t>(mDim) + 1;
        Vec3 out = {0, 0, 0};
        for (std::size_t k = 0; k < corners; ++k)
            for (std::size_t c = 0; c < 3; ++c)
                out[c] += point(s.mNode[k])[c] / static_cast<double>(corners);
        return out;
    }

private:
    std::size_t bucket(std::size_t Axis, double Coordinate) const {
        if (mRes[Axis] == 1 || !(mCell[Axis] > 0.0))
            return 0;
        const double t = (Coordinate - mLo[Axis]) / mCell[Axis];
        if (!(t > 0.0))
            return 0;
        return std::min(static_cast<std::size_t>(t), mRes[Axis] - 1);
    }

    // The barycentric coordinates of the point in a simplex; true when it is in
    // (a volume) or within `mTol` of (a surface). `rDistance` is the signed
    // distance to the plane of a triangle.
    bool Locate(std::size_t Simplex, const double* pPoint, double* pLambda, double& rDistance,
                Vec3& rNormal) const {
        const SlSimplex& s = mrSoup.mSimplices[Simplex];
        constexpr double kEps = 1e-9;
        if (mDim == 3) {
            const double* d = point(s.mNode[3]);
            const Vec3 a = detail::MathUtils::Sub(point(s.mNode[0]), d);
            const Vec3 b = detail::MathUtils::Sub(point(s.mNode[1]), d);
            const Vec3 c = detail::MathUtils::Sub(point(s.mNode[2]), d);
            const Vec3 p = detail::MathUtils::Sub(pPoint, d);
            const double det = vec3_dot(a, vec3_cross(b, c));
            if (!(std::fabs(det) > 0.0))
                return false;
            pLambda[0] = vec3_dot(p, vec3_cross(b, c)) / det;
            pLambda[1] = vec3_dot(a, vec3_cross(p, c)) / det;
            pLambda[2] = vec3_dot(a, vec3_cross(b, p)) / det;
            pLambda[3] = 1.0 - pLambda[0] - pLambda[1] - pLambda[2];
            for (int k = 0; k < 4; ++k)
                if (!(pLambda[k] >= -kEps))
                    return false;
            return true;
        }
        const double* origin = point(s.mNode[0]);
        const Vec3 u = detail::MathUtils::Sub(point(s.mNode[1]), origin);
        const Vec3 v = detail::MathUtils::Sub(point(s.mNode[2]), origin);
        const Vec3 w = detail::MathUtils::Sub(pPoint, origin);
        const Vec3 n = vec3_cross(u, v);
        const double area2 = vec3_dot(n, n);
        if (!(area2 > 0.0))
            return false;
        const double len = std::sqrt(area2);
        rNormal = {n[0] / len, n[1] / len, n[2] / len};
        rDistance = vec3_dot(w, rNormal);
        if (!(std::fabs(rDistance) <= mTol))
            return false;
        pLambda[1] = vec3_dot(vec3_cross(w, v), n) / area2;
        pLambda[2] = vec3_dot(vec3_cross(u, w), n) / area2;
        pLambda[0] = 1.0 - pLambda[1] - pLambda[2];
        for (int k = 0; k < 3; ++k)
            if (!(pLambda[k] >= -kEps))
                return false;
        return true;
    }

    const SlSoup& mrSoup;
    const double* mpXyz;
    const double* mpVec;
    double mTol;
    int mDim;
    Vec3 mLo = {0, 0, 0};
    Vec3 mHi = {0, 0, 0};
    Vec3 mCell = {1, 1, 1};
    std::array<std::size_t, 3> mRes = {1, 1, 1};
    std::vector<std::size_t> mOffset;
    std::vector<std::size_t> mItems;
};

// The unit direction of the field at a point, moved onto the surface for a
// triangle mesh; false where the field is gone (outside, zero or not finite).
bool sl_direction(const SlField& rField, double Floor, double* pPoint, double* pDir) {
    double v[3];
    if (!rField.Sample(pPoint, v))
        return false;
    const double mag = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    if (!std::isfinite(mag) || !(mag > Floor))
        return false;
    for (std::size_t c = 0; c < 3; ++c)
        pDir[c] = v[c] / mag;
    return true;
}

// One way from the seed: RK4 in arc length, stopping at the boundary, at a
// stagnation point or when the length is used up.
void sl_trace(const SlField& rField, double Floor, const Vec3& rSeed, double Sign, double Step,
              std::size_t MaxSteps, std::vector<Vec3>& rOut) {
    Vec3 p = rSeed;
    for (std::size_t it = 0; it < MaxSteps; ++it) {
        double k1[3];
        double k2[3];
        double k3[3];
        double k4[3];
        double q[3] = {p[0], p[1], p[2]};
        if (!sl_direction(rField, Floor, q, k1))
            return;
        for (std::size_t c = 0; c < 3; ++c)
            k1[c] *= Sign;
        auto stage = [&](const double* pK, double Scale, double* pOut) {
            double x[3] = {p[0] + Scale * Step * pK[0], p[1] + Scale * Step * pK[1],
                           p[2] + Scale * Step * pK[2]};
            if (!sl_direction(rField, Floor, x, pOut)) {
                for (std::size_t c = 0; c < 3; ++c)
                    pOut[c] = pK[c];  // the stage left the mesh: keep the last slope
                return;
            }
            for (std::size_t c = 0; c < 3; ++c)
                pOut[c] *= Sign;
        };
        stage(k1, 0.5, k2);
        stage(k2, 0.5, k3);
        stage(k3, 1.0, k4);
        double next[3];
        for (std::size_t c = 0; c < 3; ++c)
            next[c] = p[c] + Step / 6.0 * (k1[c] + 2.0 * k2[c] + 2.0 * k3[c] + k4[c]);
        double unused[3];
        if (!rField.Sample(next, unused))
            return;
        p = {next[0], next[1], next[2]};
        rOut.push_back(p);
    }
}

// What both tracers take from the mesh and the field before they seed: the
// simplices, the diagonal of the box around them, the fastest speed and the
// speed floor. `mOk` is false where there is nothing to trace (a zero-size mesh
// or a zero field).
struct SlSetup {
    SlSoup mSoup;
    double mDiag = 0.0;
    double mFloor = 0.0;
    double mMaxMag = 0.0;
    bool mOk = false;
};

SlSetup sl_setup(const Mesh& rMesh, const double* pXyz, const double* pVec) {
    SlSetup out;
    out.mSoup = sl_collect(rMesh);
    const SlSoup& soup = out.mSoup;
    if (soup.mSimplices.empty()) {
        std::string seen;
        for (const std::string& type : soup.mSeen)
            seen += (seen.empty() ? "" : ", ") + type;
        throw std::invalid_argument(
            "streamlines need triangle, quad, tetra, hexahedron, wedge or pyramid cells (this "
            "mesh has: " +
            (seen.empty() ? std::string("no cells") : seen) + ")");
    }
    // The box around the traced cells, and the scales taken from its diagonal.
    double lo[3] = {std::numeric_limits<double>::infinity(), 0, 0};
    double hi[3] = {0, 0, 0};
    bool first = true;
    for (const SlSimplex& s : soup.mSimplices)
        for (std::size_t k = 0; k < static_cast<std::size_t>(soup.mDim) + 1; ++k) {
            const double* p = pXyz + 3 * static_cast<std::size_t>(s.mNode[k]);
            for (std::size_t c = 0; c < 3; ++c) {
                lo[c] = first ? p[c] : std::min(lo[c], p[c]);
                hi[c] = first ? p[c] : std::max(hi[c], p[c]);
            }
            first = false;
        }
    out.mDiag = std::sqrt((hi[0] - lo[0]) * (hi[0] - lo[0]) + (hi[1] - lo[1]) * (hi[1] - lo[1]) +
                          (hi[2] - lo[2]) * (hi[2] - lo[2]));
    if (!(out.mDiag > 0.0) || !std::isfinite(out.mDiag))
        return out;
    // The slowest speed that still counts as moving, from the fastest corner.
    for (const SlSimplex& s : soup.mSimplices)
        for (std::size_t k = 0; k < static_cast<std::size_t>(soup.mDim) + 1; ++k) {
            const double* v = pVec + 3 * static_cast<std::size_t>(s.mNode[k]);
            const double m = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
            if (std::isfinite(m))
                out.mMaxMag = std::max(out.mMaxMag, m);
        }
    if (!(out.mMaxMag > 0.0))
        return out;
    out.mFloor = 1e-9 * out.mMaxMag;
    out.mOk = true;
    return out;
}

// The seeds: simplices ranked at equal measure along a golden-ratio sequence (a
// plain stride would line up with the rows of a structured mesh), each started
// at its centroid.
std::vector<Vec3> sl_seed_points(const SlField& rField, std::size_t Count, std::size_t Seeds) {
    std::vector<double> cumulative(Count);
    double total = 0.0;
    for (std::size_t s = 0; s < Count; ++s) {
        const double m = rField.Measure(s);
        total += std::isfinite(m) ? m : 0.0;
        cumulative[s] = total;
    }
    const std::size_t seeds = std::min(Seeds, Count);
    std::vector<Vec3> seed_points;
    for (std::size_t k = 0; k < seeds && total > 0.0; ++k) {
        double t = 0.5 + static_cast<double>(k) * 0.6180339887498949;
        t -= std::floor(t);
        const auto it = std::lower_bound(cumulative.begin(), cumulative.end(), t * total);
        const std::size_t s =
            std::min<std::size_t>(static_cast<std::size_t>(it - cumulative.begin()), Count - 1);
        seed_points.push_back(rField.Centroid(s));
    }
    return seed_points;
}

// The field's locator tolerance: half a render step on a surface, a whisker in
// a volume.
double sl_tolerance(int Dim, double Diag) {
    return Dim == 2 ? 0.5 * (0.005 * Diag) : 1e-9 * Diag;
}

struct SlTimePoint {
    Vec3 mP;
    Vec3 mV;  // the field's velocity at mP (not signed by the direction)
    double mT;
};

// Trace one way from the seed by true velocity (`dx/dt = sign * v`), with fixed
// RK4 or adaptive Dormand-Prince 5(4). The seed is not added. Returns why the
// line stopped.
StreamlineEnd sl_trace_timed(const SlField& rField, double Floor, const Vec3& rSeed, double Sign,
                             const StreamlineTraceOptions& rOptions, double StepTime,
                             double MaxLength, double ToleranceLength, std::size_t MaxSteps,
                             std::vector<SlTimePoint>& rOut) {
    // The velocity at x, signed by the direction; false outside the mesh.
    auto velocity = [&](const Vec3& x, Vec3& v) {
        double q[3] = {x[0], x[1], x[2]};
        double vv[3];
        if (!rField.Sample(q, vv))
            return false;
        v = {Sign * vv[0], Sign * vv[1], Sign * vv[2]};
        return true;
    };
    auto add = [](const Vec3& a, const Vec3& b, double s) {
        return Vec3{a[0] + s * b[0], a[1] + s * b[1], a[2] + s * b[2]};
    };
    auto norm = [](const Vec3& a) { return std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]); };

    Vec3 p = rSeed;
    Vec3 v;
    if (!velocity(p, v))
        return StreamlineEnd::LeftDomain;
    double t = 0.0;
    double length = 0.0;
    double h = StepTime;
    for (std::size_t it = 0; it < MaxSteps; ++it) {
        if (!(norm(v) > Floor))
            return StreamlineEnd::Stagnation;
        Vec3 next;
        Vec3 vnext;
        double used = h;
        if (!rOptions.mAdaptive) {
            // Fixed RK4. A stage that leaves the mesh keeps the last slope.
            Vec3 k1 = v;
            Vec3 k2;
            Vec3 k3;
            Vec3 k4;
            if (!velocity(add(p, k1, 0.5 * h), k2))
                k2 = k1;
            if (!velocity(add(p, k2, 0.5 * h), k3))
                k3 = k2;
            if (!velocity(add(p, k3, h), k4))
                k4 = k3;
            for (std::size_t c = 0; c < 3; ++c)
                next[c] = p[c] + h / 6.0 * (k1[c] + 2.0 * k2[c] + 2.0 * k3[c] + k4[c]);
            if (!velocity(next, vnext))
                return StreamlineEnd::LeftDomain;
        } else {
            // Dormand-Prince 5(4): retry the step with a smaller h until its
            // error estimate is within the tolerance; a stage outside the mesh
            // counts as a failure.
            bool accepted = false;
            for (int attempt = 0; attempt < 16 && !accepted; ++attempt) {
                Vec3 k[7];
                k[0] = v;
                bool ok = true;
                auto stage = [&](int i, std::initializer_list<std::pair<int, double>> terms) {
                    Vec3 x = p;
                    for (const auto& term : terms)
                        x = add(x, k[term.first], h * term.second);
                    if (!velocity(x, k[i]))
                        ok = false;
                };
                stage(1, {{0, 1.0 / 5.0}});
                if (ok)
                    stage(2, {{0, 3.0 / 40.0}, {1, 9.0 / 40.0}});
                if (ok)
                    stage(3, {{0, 44.0 / 45.0}, {1, -56.0 / 15.0}, {2, 32.0 / 9.0}});
                if (ok)
                    stage(4, {{0, 19372.0 / 6561.0},
                              {1, -25360.0 / 2187.0},
                              {2, 64448.0 / 6561.0},
                              {3, -212.0 / 729.0}});
                if (ok)
                    stage(5, {{0, 9017.0 / 3168.0},
                              {1, -355.0 / 33.0},
                              {2, 46732.0 / 5247.0},
                              {3, 49.0 / 176.0},
                              {4, -5103.0 / 18656.0}});
                if (ok) {
                    // The fifth-order solution, which is also where k7 is taken (FSAL).
                    Vec3 y = p;
                    const double b[6] = {35.0 / 384.0,     0.0,
                                         500.0 / 1113.0,   125.0 / 192.0,
                                         -2187.0 / 6784.0, 11.0 / 84.0};
                    for (int i = 0; i < 6; ++i)
                        y = add(y, k[i], h * b[i]);
                    if (!velocity(y, k[6])) {
                        ok = false;
                    } else {
                        const double e[7] = {
                            71.0 / 57600.0,      0.0,          -71.0 / 16695.0, 71.0 / 1920.0,
                            -17253.0 / 339200.0, 22.0 / 525.0, -1.0 / 40.0};
                        Vec3 err = {0, 0, 0};
                        for (int i = 0; i < 7; ++i)
                            err = add(err, k[i], h * e[i]);
                        const double en = norm(err);
                        if (en <= ToleranceLength) {
                            accepted = true;
                            next = y;
                            vnext = k[6];
                            used = h;
                            // Grow for the next step, within bounds of the first.
                            const double grow =
                                en > 0.0 ? 0.9 * std::pow(ToleranceLength / en, 0.2) : 5.0;
                            h = std::clamp(h * std::clamp(grow, 0.2, 5.0), StepTime * 1e-6,
                                           StepTime * 50.0);
                        } else {
                            h = std::max(
                                h * std::clamp(0.9 * std::pow(ToleranceLength / en, 0.2), 0.1, 0.9),
                                StepTime * 1e-6);
                        }
                    }
                }
                if (!ok)
                    h *= 0.5;
            }
            if (!accepted)
                return StreamlineEnd::LeftDomain;
        }
        length += norm(Vec3{next[0] - p[0], next[1] - p[1], next[2] - p[2]});
        t += used;
        p = next;
        v = vnext;
        SlTimePoint pt;
        pt.mP = p;
        pt.mV = {Sign * vnext[0], Sign * vnext[1], Sign * vnext[2]};  // the field's own velocity
        pt.mT = Sign * t;
        rOut.push_back(pt);
        if (length >= MaxLength)
            return StreamlineEnd::MaxLength;
    }
    return StreamlineEnd::MaxSteps;
}

}  // namespace

Streamlines trace_streamlines(const Mesh& rMesh, const double* pXyz, const double* pVec,
                              const StreamlineOptions& rOptions) {
    const SlSetup setup = sl_setup(rMesh, pXyz, pVec);
    const SlSoup& soup = setup.mSoup;
    Streamlines out;
    out.mDim = soup.mDim;
    out.mNumSimplices = soup.mSimplices.size();
    if (!setup.mOk)
        return out;
    const double diag = setup.mDiag;
    const double floor = setup.mFloor;
    const double step = 0.005 * diag;
    const std::size_t max_steps = static_cast<std::size_t>(
        std::clamp(std::ceil(rOptions.mLength * diag / step), 1.0, 2000.0));

    const SlField field(soup, pXyz, pVec, diag, sl_tolerance(soup.mDim, diag));
    const std::vector<Vec3> seed_points =
        sl_seed_points(field, soup.mSimplices.size(), rOptions.mSeeds);

    std::vector<std::vector<Vec3>> lines(seed_points.size());
    parallel_for(seed_points.size(), [&](std::size_t i) {
        std::vector<Vec3> back;
        std::vector<Vec3> forward;
        Vec3 seed = seed_points[i];
        double v[3];
        if (!field.Sample(seed.data(), v))
            return;
        sl_trace(field, floor, seed, -1.0, step, max_steps, back);
        sl_trace(field, floor, seed, +1.0, step, max_steps, forward);
        std::vector<Vec3>& line = lines[i];
        line.assign(back.rbegin(), back.rend());
        line.push_back(seed);
        line.insert(line.end(), forward.begin(), forward.end());
        if (line.size() < 2)
            line.clear();
    });
    for (const std::vector<Vec3>& line : lines) {
        if (line.empty())
            continue;
        for (const Vec3& p : line)
            out.mXyz.insert(out.mXyz.end(), p.begin(), p.end());
        out.mStart.push_back(out.mXyz.size() / 3);
    }
    return out;
}

StreamlineTraces trace_streamlines_timed(const Mesh& rMesh, const double* pXyz, const double* pVec,
                                         const StreamlineTraceOptions& rOptions) {
    const SlSetup setup = sl_setup(rMesh, pXyz, pVec);
    const SlSoup& soup = setup.mSoup;
    StreamlineTraces out;
    out.mDim = soup.mDim;
    if (!setup.mOk)
        return out;
    const double diag = setup.mDiag;
    const double floor = rOptions.mTerminalSpeed > 0.0 ? rOptions.mTerminalSpeed : setup.mFloor;
    const double step_length = rOptions.mStep > 0.0 ? rOptions.mStep : 0.005 * diag;
    const double step_time = step_length / setup.mMaxMag;
    const double max_length = rOptions.mMaxLength * diag;
    const std::size_t max_steps = rOptions.mMaxSteps > 0 ? rOptions.mMaxSteps : 20000;
    const double tolerance = rOptions.mTolerance * diag;

    const SlField field(soup, pXyz, pVec, diag, sl_tolerance(soup.mDim, diag));
    std::vector<Vec3> seeds;
    if (!rOptions.mSeedPoints.empty()) {
        for (std::size_t i = 0; i + 2 < rOptions.mSeedPoints.size(); i += 3)
            seeds.push_back({rOptions.mSeedPoints[i], rOptions.mSeedPoints[i + 1],
                             rOptions.mSeedPoints[i + 2]});
    } else {
        seeds = sl_seed_points(field, soup.mSimplices.size(), rOptions.mSeeds);
    }

    struct SlOne {
        bool mValid = false;
        std::vector<SlTimePoint> mBack;
        std::vector<SlTimePoint> mForward;
        Vec3 mSeed = {0, 0, 0};
        Vec3 mSeedVelocity = {0, 0, 0};
        StreamlineEnd mBackEnd = StreamlineEnd::LeftDomain;
        StreamlineEnd mForwardEnd = StreamlineEnd::LeftDomain;
    };
    std::vector<SlOne> traced(seeds.size());
    parallel_for(seeds.size(), [&](std::size_t i) {
        SlOne& one = traced[i];
        Vec3 seed = seeds[i];
        double v[3];
        if (!field.Sample(seed.data(), v))
            return;  // the seed (moved onto a surface) is outside the mesh
        one.mValid = true;
        one.mSeed = seed;
        one.mSeedVelocity = {v[0], v[1], v[2]};
        if (rOptions.mBackward)
            one.mBackEnd = sl_trace_timed(field, floor, seed, -1.0, rOptions, step_time, max_length,
                                          tolerance, max_steps, one.mBack);
        if (rOptions.mForward)
            one.mForwardEnd = sl_trace_timed(field, floor, seed, +1.0, rOptions, step_time,
                                             max_length, tolerance, max_steps, one.mForward);
    });
    for (std::size_t i = 0; i < traced.size(); ++i) {
        const SlOne& one = traced[i];
        if (!one.mValid) {
            ++out.mNumSeedsOutside;
            continue;
        }
        if (one.mBack.size() + one.mForward.size() < 1)
            continue;  // a point, not a line
        StreamlineTrace line;
        line.mSeedId = static_cast<std::int64_t>(i);
        line.mBackwardEnd = one.mBackEnd;
        line.mForwardEnd = one.mForwardEnd;
        auto push = [&](const Vec3& p, const Vec3& v, double t) {
            line.mXyz.insert(line.mXyz.end(), p.begin(), p.end());
            line.mVelocity.insert(line.mVelocity.end(), v.begin(), v.end());
            line.mTime.push_back(t);
        };
        for (auto it = one.mBack.rbegin(); it != one.mBack.rend(); ++it)
            push(it->mP, it->mV, it->mT);
        push(one.mSeed, one.mSeedVelocity, 0.0);
        for (const SlTimePoint& pt : one.mForward)
            push(pt.mP, pt.mV, pt.mT);
        out.mLines.push_back(std::move(line));
    }
    return out;
}

}  // namespace detail
}  // namespace meshioplusplus

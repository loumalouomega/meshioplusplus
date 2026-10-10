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

// Project includes (private, not installed)
#include "streamlines.hpp"

namespace meshioplusplus {
namespace detail {
namespace {

using Vec3 = std::array<double, 3>;

Vec3 sl_sub(const double* pA, const double* pB) {
    return {pA[0] - pB[0], pA[1] - pB[1], pA[2] - pB[2]};
}

Vec3 sl_cross(const Vec3& rA, const Vec3& rB) {
    return {rA[1] * rB[2] - rA[2] * rB[1], rA[2] * rB[0] - rA[0] * rB[2],
            rA[0] * rB[1] - rA[1] * rB[0]};
}

double sl_dot(const Vec3& rA, const Vec3& rB) {
    return rA[0] * rB[0] + rA[1] * rB[1] + rA[2] * rB[2];
}

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
                const double along = sl_dot({best_vec[0], best_vec[1], best_vec[2]}, normal);
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
            const Vec3 a = sl_sub(point(s.mNode[0]), point(s.mNode[3]));
            const Vec3 b = sl_sub(point(s.mNode[1]), point(s.mNode[3]));
            const Vec3 c = sl_sub(point(s.mNode[2]), point(s.mNode[3]));
            return std::fabs(sl_dot(a, sl_cross(b, c))) / 6.0;
        }
        const Vec3 a = sl_sub(point(s.mNode[1]), point(s.mNode[0]));
        const Vec3 b = sl_sub(point(s.mNode[2]), point(s.mNode[0]));
        const Vec3 n = sl_cross(a, b);
        return 0.5 * std::sqrt(sl_dot(n, n));
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
            const Vec3 a = sl_sub(point(s.mNode[0]), d);
            const Vec3 b = sl_sub(point(s.mNode[1]), d);
            const Vec3 c = sl_sub(point(s.mNode[2]), d);
            const Vec3 p = sl_sub(pPoint, d);
            const double det = sl_dot(a, sl_cross(b, c));
            if (!(std::fabs(det) > 0.0))
                return false;
            pLambda[0] = sl_dot(p, sl_cross(b, c)) / det;
            pLambda[1] = sl_dot(a, sl_cross(p, c)) / det;
            pLambda[2] = sl_dot(a, sl_cross(b, p)) / det;
            pLambda[3] = 1.0 - pLambda[0] - pLambda[1] - pLambda[2];
            for (int k = 0; k < 4; ++k)
                if (!(pLambda[k] >= -kEps))
                    return false;
            return true;
        }
        const double* origin = point(s.mNode[0]);
        const Vec3 u = sl_sub(point(s.mNode[1]), origin);
        const Vec3 v = sl_sub(point(s.mNode[2]), origin);
        const Vec3 w = sl_sub(pPoint, origin);
        const Vec3 n = sl_cross(u, v);
        const double area2 = sl_dot(n, n);
        if (!(area2 > 0.0))
            return false;
        const double len = std::sqrt(area2);
        rNormal = {n[0] / len, n[1] / len, n[2] / len};
        rDistance = sl_dot(w, rNormal);
        if (!(std::fabs(rDistance) <= mTol))
            return false;
        pLambda[1] = sl_dot(sl_cross(w, v), n) / area2;
        pLambda[2] = sl_dot(sl_cross(u, w), n) / area2;
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

}  // namespace

Streamlines trace_streamlines(const Mesh& rMesh, const double* pXyz, const double* pVec,
                              const StreamlineOptions& rOptions) {
    const SlSoup soup = sl_collect(rMesh);
    Streamlines out;
    if (soup.mSimplices.empty()) {
        std::string seen;
        for (const std::string& type : soup.mSeen)
            seen += (seen.empty() ? "" : ", ") + type;
        throw std::invalid_argument(
            "streamlines need triangle, quad, tetra, hexahedron, wedge or pyramid cells (this "
            "mesh has: " +
            (seen.empty() ? std::string("no cells") : seen) + ")");
    }
    out.mDim = soup.mDim;
    out.mNumSimplices = soup.mSimplices.size();

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
    const double diag =
        std::sqrt((hi[0] - lo[0]) * (hi[0] - lo[0]) + (hi[1] - lo[1]) * (hi[1] - lo[1]) +
                  (hi[2] - lo[2]) * (hi[2] - lo[2]));
    if (!(diag > 0.0) || !std::isfinite(diag))
        return out;
    const double step = 0.005 * diag;
    const std::size_t max_steps = static_cast<std::size_t>(
        std::clamp(std::ceil(rOptions.mLength * diag / step), 1.0, 2000.0));

    // The slowest speed that still counts as moving, from the fastest corner.
    double max_mag = 0.0;
    for (const SlSimplex& s : soup.mSimplices)
        for (std::size_t k = 0; k < static_cast<std::size_t>(soup.mDim) + 1; ++k) {
            const double* v = pVec + 3 * static_cast<std::size_t>(s.mNode[k]);
            const double m = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
            if (std::isfinite(m))
                max_mag = std::max(max_mag, m);
        }
    if (!(max_mag > 0.0))
        return out;
    const double floor = 1e-9 * max_mag;

    const SlField field(soup, pXyz, pVec, diag, soup.mDim == 2 ? 0.5 * step : 1e-9 * diag);

    // Seeds: simplices ranked at equal measure along a golden-ratio sequence
    // (a plain stride would line up with the rows of a structured mesh), each
    // started at its centroid.
    const std::size_t count = soup.mSimplices.size();
    std::vector<double> cumulative(count);
    double total = 0.0;
    for (std::size_t s = 0; s < count; ++s) {
        const double m = field.Measure(s);
        total += std::isfinite(m) ? m : 0.0;
        cumulative[s] = total;
    }
    const std::size_t seeds = std::min(rOptions.mSeeds, count);
    std::vector<Vec3> seed_points;
    for (std::size_t k = 0; k < seeds && total > 0.0; ++k) {
        double t = 0.5 + static_cast<double>(k) * 0.6180339887498949;
        t -= std::floor(t);
        const auto it = std::lower_bound(cumulative.begin(), cumulative.end(), t * total);
        const std::size_t s =
            std::min<std::size_t>(static_cast<std::size_t>(it - cumulative.begin()), count - 1);
        seed_points.push_back(field.Centroid(s));
    }

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

}  // namespace detail
}  // namespace meshioplusplus

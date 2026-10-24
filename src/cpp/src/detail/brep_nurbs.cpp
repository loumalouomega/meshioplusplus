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

// Kernel-free trimmed NURBS faces. See detail/brep_nurbs.hpp.

// System includes
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <utility>
#include <vector>

// Project includes (private, not installed)
#include "brep_nurbs.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

constexpr double kBnPi = 3.14159265358979323846;
constexpr double kBnTwoPi = 2.0 * kBnPi;

double bn_dot(const PatchVec& rA, const PatchVec& rB) {
    return rA[0] * rB[0] + rA[1] * rB[1] + rA[2] * rB[2];
}

double bn_dist(const PatchVec& rA, const PatchVec& rB) {
    const double d[3] = {rA[0] - rB[0], rA[1] - rB[1], rA[2] - rB[2]};
    return std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
}

// The knot span of `T` (clamped to the domain) of `N` poles of degree `P`.
std::size_t bn_span(int P, const std::vector<double>& rK, std::size_t N, double T) {
    const std::size_t p = static_cast<std::size_t>(P);
    if (T >= rK[N])
        return N - 1;
    if (T <= rK[p])
        return p;
    std::size_t lo = p, hi = N;
    while (hi - lo > 1) {
        const std::size_t mid = (lo + hi) / 2;
        if (T < rK[mid])
            hi = mid;
        else
            lo = mid;
    }
    return lo;
}

// The non-zero basis functions at `T` on `Span` (The NURBS Book, A2.2).
void bn_basis(std::size_t Span, double T, int P, const std::vector<double>& rK, double* pN) {
    double left[kBsplineMaxDegree + 2], right[kBsplineMaxDegree + 2];
    pN[0] = 1.0;
    for (int j = 1; j <= P; ++j) {
        left[j] = T - rK[Span + 1 - static_cast<std::size_t>(j)];
        right[j] = rK[Span + static_cast<std::size_t>(j)] - T;
        double saved = 0.0;
        for (int r = 0; r < j; ++r) {
            const double den = right[r + 1] + left[j - r];
            const double tmp = den != 0.0 ? pN[r] / den : 0.0;
            pN[r] = saved + right[r + 1] * tmp;
            saved = left[j - r] * tmp;
        }
        pN[j] = saved;
    }
}

// `x` moved by whole turns to within half a turn of `Near`.
double bn_unwrap(double X, double Near) {
    return X + kBnTwoPi * std::round((Near - X) / kBnTwoPi);
}

// A loop point: its position and its raw parameters (the plane's
// coordinates, a revolution's angle and height or latitude, a B-spline's
// own), and its parameter on the trim's edge curve.
struct BnSample {
    PatchVec mP{0.0, 0.0, 0.0};
    double mX[2] = {0.0, 0.0};
    double mT = 0.0;
};

// A trim being built: its samples along the loop; `pSrc` the given trim run
// backwards when `mRev`, or none for a seam (`mSeam`) or pole (`mPole`).
struct BnTrim {
    std::vector<BnSample> mS;
    const BrepNurbsTrim* pSrc = nullptr;
    bool mRev = false;
    std::int64_t mSeam = -1;
    bool mSeamAlong = true;
    bool mPole = false;
};

struct BnLoop {
    std::vector<BnTrim> mTrims;
    int mWind[2] = {0, 0};
};

// The point of a given trim at `T` along the loop.
PatchVec bn_at(const BnTrim& rT, double T) {
    return rT.pSrc->mAt(rT.mRev ? 1.0 - T : T);
}

// A loop's samples in order, each trim's first point included once.
template <class F>
void bn_each(const BnLoop& rL, F&& rF) {
    for (const BnTrim& t : rL.mTrims)
        for (std::size_t i = 0; i + 1 < t.mS.size(); ++i)
            rF(t.mS[i]);
}

// Twice the signed area of a loop in the given parameters (of `Map`).
template <class F>
double bn_area(const BnLoop& rL, F&& rMap) {
    std::vector<std::array<double, 2>> ring;
    bn_each(rL, [&](const BnSample& rS) { ring.push_back(rMap(rS)); });
    double a = 0.0;
    for (std::size_t i = 0; i < ring.size(); ++i) {
        const auto& p = ring[i];
        const auto& q = ring[(i + 1) % ring.size()];
        a += p[0] * q[1] - p[1] * q[0];
    }
    return a;
}

}  // namespace

UnitArc unit_arc(double A0, double A1) {
    UnitArc arc;
    arc.mA0 = A0;
    arc.mA1 = A1;
    const double span = A1 - A0;
    arc.mSegments =
        static_cast<std::size_t>(std::max(1.0, std::ceil(span / (0.5 * kBnPi) - 1e-12)));
    const double step = span / static_cast<double>(arc.mSegments);
    const double half = 0.5 * step;
    const double w = std::cos(half);
    arc.mKnots.assign(3, A0);
    for (std::size_t k = 1; k < arc.mSegments; ++k) {
        const double t = A0 + static_cast<double>(k) * step;
        arc.mKnots.push_back(t);
        arc.mKnots.push_back(t);
    }
    for (int i = 0; i < 3; ++i)
        arc.mKnots.push_back(A1);
    for (std::size_t k = 0; k <= arc.mSegments; ++k) {
        const double a = A0 + static_cast<double>(k) * step;
        arc.mPoles.push_back({std::cos(a), std::sin(a), 1.0});
        if (k < arc.mSegments) {
            const double m = a + half;
            arc.mPoles.push_back({std::cos(m) / w, std::sin(m) / w, w});
        }
    }
    return arc;
}

double unit_arc_param(const UnitArc& rArc, double A) {
    const double step = (rArc.mA1 - rArc.mA0) / static_cast<double>(rArc.mSegments);
    const double half = 0.5 * step;
    double k = std::floor((A - rArc.mA0) / step);
    k = std::min(std::max(k, 0.0), static_cast<double>(rArc.mSegments) - 1.0);
    const double psi = A - (rArc.mA0 + k * step + half);
    const double s = 0.5 * (std::tan(0.5 * psi) / std::tan(0.5 * half) + 1.0);
    return rArc.mA0 + (k + s) * step;
}

PatchVec brep_curve_point(const BrepCurve& rC, double T) {
    const std::size_t n = rC.mControlPoints.Shape()[0];
    const std::size_t span = bn_span(rC.mDegree, rC.mKnots, n, T);
    double nb[kBsplineMaxDegree + 1];
    bn_basis(span, T, rC.mDegree, rC.mKnots, nb);
    const double* xyz = rC.mControlPoints.As<double>();
    const double* w = rC.mWeights.As<double>();
    double acc[4] = {0.0, 0.0, 0.0, 0.0};
    for (int i = 0; i <= rC.mDegree; ++i) {
        const std::size_t k =
            span - static_cast<std::size_t>(rC.mDegree) + static_cast<std::size_t>(i);
        const double wk = rC.mRational ? w[k] : 1.0;
        for (int j = 0; j < 3; ++j)
            acc[j] += nb[i] * wk * xyz[k * 3 + static_cast<std::size_t>(j)];
        acc[3] += nb[i] * wk;
    }
    return {acc[0] / acc[3], acc[1] / acc[3], acc[2] / acc[3]};
}

PatchVec nurbs_patch_point(const NurbsPatch& rP, double U, double V) {
    const int p = rP.mDegrees[0], q = rP.mDegrees[1];
    const std::size_t nu = rP.mKnots[0].size() - static_cast<std::size_t>(p) - 1;
    const std::size_t nv = rP.mKnots[1].size() - static_cast<std::size_t>(q) - 1;
    const std::size_t su = bn_span(p, rP.mKnots[0], nu, U);
    const std::size_t sv = bn_span(q, rP.mKnots[1], nv, V);
    double bu[kBsplineMaxDegree + 1], bv[kBsplineMaxDegree + 1];
    bn_basis(su, U, p, rP.mKnots[0], bu);
    bn_basis(sv, V, q, rP.mKnots[1], bv);
    const double* xyz = rP.mControlPoints.As<double>();
    const double* w = rP.mWeights.As<double>();
    double acc[4] = {0.0, 0.0, 0.0, 0.0};
    for (int b = 0; b <= q; ++b) {
        const std::size_t j = sv - static_cast<std::size_t>(q) + static_cast<std::size_t>(b);
        for (int a = 0; a <= p; ++a) {
            const std::size_t i = su - static_cast<std::size_t>(p) + static_cast<std::size_t>(a);
            const std::size_t k = j * nu + i;
            const double c = bu[a] * bv[b] * w[k];
            for (int d = 0; d < 3; ++d)
                acc[d] += c * xyz[k * 3 + static_cast<std::size_t>(d)];
            acc[3] += c;
        }
    }
    return {acc[0] / acc[3], acc[1] / acc[3], acc[2] / acc[3]};
}

BrepCurve brep_curve_line(const PatchVec& rA, const PatchVec& rB) {
    return brep_curve_polyline({rA, rB});
}

BrepCurve brep_curve_polyline(const std::vector<PatchVec>& rPts) {
    BrepCurve c;
    const std::size_t n = rPts.size();
    c.mDegree = 1;
    c.mKnots.push_back(0.0);
    for (std::size_t i = 0; i < n; ++i)
        c.mKnots.push_back(static_cast<double>(i));
    c.mKnots.push_back(static_cast<double>(n - 1));
    c.mControlPoints = NDArray(DType::Float64, {n, 3});
    c.mWeights = NDArray(DType::Float64, {n});
    for (std::size_t i = 0; i < n; ++i) {
        for (int j = 0; j < 3; ++j)
            c.mControlPoints.As<double>()[i * 3 + static_cast<std::size_t>(j)] =
                rPts[i][static_cast<std::size_t>(j)];
        c.mWeights.As<double>()[i] = 1.0;
    }
    c.mActiveRange = {0.0, static_cast<double>(n - 1)};
    return c;
}

BrepCurve brep_curve_arc(const PatchVec& rCentre, const PatchVec& rX, const PatchVec& rY, double R,
                         double A0, double A1) {
    const UnitArc arc = unit_arc(std::min(A0, A1), std::max(A0, A1));
    BrepCurve c;
    c.mDegree = 2;
    c.mRational = true;
    c.mKnots = arc.mKnots;
    const std::size_t n = arc.mPoles.size();
    c.mControlPoints = NDArray(DType::Float64, {n, 3});
    c.mWeights = NDArray(DType::Float64, {n});
    for (std::size_t i = 0; i < n; ++i) {
        const auto& q = arc.mPoles[i];
        for (std::size_t j = 0; j < 3; ++j)
            c.mControlPoints.As<double>()[i * 3 + j] =
                rCentre[j] + R * (q[0] * rX[j] + q[1] * rY[j]);
        c.mWeights.As<double>()[i] = q[2];
    }
    c.mActiveRange = {arc.mA0, arc.mA1};
    return A1 < A0 ? brep_curve_reverse(c) : c;
}

BrepCurve brep_curve_reverse(const BrepCurve& rC) {
    BrepCurve c = rC;
    const std::size_t n = rC.mControlPoints.Shape()[0];
    const double a = rC.mKnots.front(), b = rC.mKnots.back();
    for (std::size_t i = 0; i < rC.mKnots.size(); ++i)
        c.mKnots[i] = a + b - rC.mKnots[rC.mKnots.size() - 1 - i];
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < 3; ++j)
            c.mControlPoints.As<double>()[i * 3 + j] =
                rC.mControlPoints.As<double>()[(n - 1 - i) * 3 + j];
        c.mWeights.As<double>()[i] = rC.mWeights.As<double>()[n - 1 - i];
    }
    if (rC.mIds.Size() == n)
        for (std::size_t i = 0; i < n; ++i)
            c.mIds.As<std::int64_t>()[i] = rC.mIds.As<std::int64_t>()[n - 1 - i];
    if (rC.mActiveRange.size() == 2)
        c.mActiveRange = {a + b - rC.mActiveRange[1], a + b - rC.mActiveRange[0]};
    return c;
}

BrepCurve brep_curve_bspline(const BsplineCurve& rC) {
    BrepCurve c;
    c.mDegree = rC.mDegree;
    c.mKnots = rC.mKnots;
    const std::size_t n = rC.mPoles.size();
    c.mControlPoints = NDArray(DType::Float64, {n, 3});
    c.mWeights = NDArray(DType::Float64, {n});
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < 3; ++j)
            c.mControlPoints.As<double>()[i * 3 + j] = rC.mPoles[i][j];
        c.mWeights.As<double>()[i] = 1.0;
    }
    c.mActiveRange = {rC.T0(), rC.T1()};
    return c;
}

bool brep_nurbs_face(const BrepNurbsSurface& rSurface,
                     const std::vector<std::vector<BrepNurbsTrim>>& rLoops, bool Sense, double Tol,
                     BrepNurbsFace& rOut, std::string& rWhy) {
    rOut = BrepNurbsFace{};
    const int kind = rSurface.mKind;
    const bool revolution = kind >= 2 && kind <= 5;
    const bool periodic[2] = {revolution, kind == 5};
    if (rLoops.empty()) {
        rWhy = "a face without loops";
        return false;
    }
    if (kind == 6 && rSurface.pBspline == nullptr) {
        rWhy = "a B-spline face without its surface";
        return false;
    }

    // The scale, for the tests of degenerate positions.
    double lo[3] = {1e300, 1e300, 1e300}, hi[3] = {-1e300, -1e300, -1e300};
    for (const auto& loop : rLoops)
        for (const BrepNurbsTrim& t : loop) {
            const auto take = [&](const PatchVec& rP) {
                for (std::size_t j = 0; j < 3; ++j) {
                    lo[j] = std::min(lo[j], rP[j]);
                    hi[j] = std::max(hi[j], rP[j]);
                }
            };
            if (t.mAt) {
                take(t.mAt(0.0));
                take(t.mAt(0.5));
                take(t.mAt(1.0));
            } else {
                for (const PatchVec& p : t.mPoints)
                    take(p);
            }
        }
    const double scale = std::max(
        1e-300, std::sqrt((hi[0] - lo[0]) * (hi[0] - lo[0]) + (hi[1] - lo[1]) * (hi[1] - lo[1]) +
                          (hi[2] - lo[2]) * (hi[2] - lo[2])));
    const double tiny = 1e-9 * scale;

    // The raw parameters of a point; `pNear` the previous sample's, for the
    // unwrapping and a B-spline's warm start.
    const BsplineSurface* bs = rSurface.pBspline;
    bool raw_ok = true;
    const auto raw = [&](const PatchVec& rP, const BnSample* pNear, double* pX) {
        const PatchVec d{rP[0] - rSurface.mOrigin[0], rP[1] - rSurface.mOrigin[1],
                         rP[2] - rSurface.mOrigin[2]};
        const double dx = bn_dot(d, rSurface.mX), dy = bn_dot(d, rSurface.mY),
                     dz = bn_dot(d, rSurface.mZ);
        switch (kind) {
            case 1:
                pX[0] = dx;
                pX[1] = dy;
                return;
            case 6: {
                double u = 0.0, v = 0.0;
                double best = std::numeric_limits<double>::infinity();
                if (pNear != nullptr) {
                    u = pNear->mX[0];
                    v = pNear->mX[1];
                    best = bspline_project(*bs, rP.data(), u, v);
                }
                if (!(best <= 4.0 * Tol + tiny)) {
                    // Global starts: the eight nearest nodes of a grid over the
                    // domain, each projected (a strongly curved face can bring
                    // another part of itself nearer to a grid node than the
                    // part the point lies on).
                    constexpr int kGrid = 24;
                    constexpr std::size_t kStarts = 8;
                    std::vector<std::pair<double, std::pair<double, double>>> nodes;
                    nodes.reserve((kGrid + 1) * (kGrid + 1));
                    for (int i = 0; i <= kGrid; ++i)
                        for (int j = 0; j <= kGrid; ++j) {
                            const double gu = bs->U0() + (bs->U1() - bs->U0()) * i / kGrid;
                            const double gv = bs->V0() + (bs->V1() - bs->V0()) * j / kGrid;
                            BsplineEval e;
                            bspline_eval(*bs, gu, gv, 0, e);
                            nodes.push_back({bn_dist({e.mS[0], e.mS[1], e.mS[2]}, rP), {gu, gv}});
                        }
                    const std::size_t starts = std::min(kStarts, nodes.size());
                    std::partial_sort(nodes.begin(),
                                      nodes.begin() + static_cast<std::ptrdiff_t>(starts),
                                      nodes.end());
                    for (std::size_t k = 0; k < starts; ++k) {
                        double u2 = nodes[k].second.first, v2 = nodes[k].second.second;
                        const double d2 = bspline_project(*bs, rP.data(), u2, v2);
                        if (d2 < best) {
                            best = d2;
                            u = u2;
                            v = v2;
                        }
                    }
                }
                pX[0] = u;
                pX[1] = v;
                return;
            }
            default:
                break;
        }
        const double rho = std::sqrt(dx * dx + dy * dy);
        if (!(rho > tiny)) {
            raw_ok = false;
            return;
        }
        pX[0] = std::atan2(dy, dx);
        if (kind == 3 || kind == 4) {
            pX[1] = dz;
        } else if (kind == 2) {
            pX[1] = std::atan2(dz, rho);
            if (std::fabs(pX[1]) > 0.5 * kBnPi - 1e-6)
                raw_ok = false;
        } else {
            pX[1] = std::atan2(dz, rho - rSurface.mRadius);
        }
        if (pNear != nullptr)
            for (int k = 0; k < 2; ++k)
                if (periodic[k])
                    pX[k] = bn_unwrap(pX[k], pNear->mX[k]);
    };

    // The loops' samples, against the surface's normal when it points in.
    std::vector<BnLoop> loops(rLoops.size());
    for (std::size_t li = 0; li < rLoops.size(); ++li) {
        const std::vector<BrepNurbsTrim>& src = rLoops[li];
        const std::size_t n = src.size();
        BnSample prev;
        bool have_prev = false;
        for (std::size_t k = 0; k < n; ++k) {
            BnTrim t;
            t.pSrc = &src[Sense ? k : n - 1 - k];
            t.mRev = !Sense;
            std::vector<std::pair<PatchVec, double>> pts;
            if (t.pSrc->mAt) {
                const int m = t.pSrc->mCircle ? 16 : 8;
                for (int i = 0; i <= m; ++i) {
                    const double s = static_cast<double>(i) / m;
                    pts.emplace_back(bn_at(t, s), s);
                }
            } else {
                const std::size_t m = t.pSrc->mPoints.size();
                if (m < 2) {
                    rWhy = "a trim with fewer than two points";
                    return false;
                }
                for (std::size_t i = 0; i < m; ++i)
                    pts.emplace_back(t.pSrc->mPoints[t.mRev ? m - 1 - i : i],
                                     static_cast<double>(i) / static_cast<double>(m - 1));
            }
            for (const auto& pt : pts) {
                BnSample s;
                s.mP = pt.first;
                s.mT = pt.second;
                raw(s.mP, have_prev ? &prev : nullptr, s.mX);
                if (!raw_ok) {
                    rWhy = "a loop runs through the surface's axis or pole";
                    return false;
                }
                t.mS.push_back(s);
                prev = s;
                have_prev = true;
            }
            loops[li].mTrims.push_back(std::move(t));
        }
        const BnSample& first = loops[li].mTrims.front().mS.front();
        const BnSample& last = loops[li].mTrims.back().mS.back();
        for (int k = 0; k < 2; ++k)
            if (periodic[k])
                loops[li].mWind[k] =
                    static_cast<int>(std::lround((last.mX[k] - first.mX[k]) / kBnTwoPi));
    }

    // Which direction the loops wind round, if any.
    int wrap = -1;
    std::vector<std::size_t> winding;
    for (std::size_t li = 0; li < loops.size(); ++li) {
        const int* w = loops[li].mWind;
        if (w[0] != 0 && w[1] != 0) {
            rWhy = "a loop winds round both directions of a torus";
            return false;
        }
        for (int k = 0; k < 2; ++k)
            if (w[k] != 0) {
                if (std::abs(w[k]) > 1) {
                    rWhy = "a loop winds round more than once";
                    return false;
                }
                if (wrap >= 0 && wrap != k) {
                    rWhy = "loops wind round different directions";
                    return false;
                }
                wrap = k;
                winding.push_back(li);
            }
    }

    // Whole turns: a loop moved by them, and a loop's mean parameter.
    const auto shift = [&](BnLoop& rL, int Dim, double By) {
        for (BnTrim& t : rL.mTrims)
            for (BnSample& s : t.mS)
                s.mX[Dim] += By;
    };
    const auto mean = [&](const BnLoop& rL, int Dim) {
        double sum = 0.0;
        std::size_t count = 0;
        bn_each(rL, [&](const BnSample& rS) {
            sum += rS.mX[Dim];
            ++count;
        });
        return sum / static_cast<double>(std::max<std::size_t>(count, 1));
    };
    const auto raw_xy = [](const BnSample& rS) {
        return std::array<double, 2>{rS.mX[0], rS.mX[1]};
    };

    // The composed outer loop of a face its loops wind round, and the loops
    // that are not part of it.
    std::vector<std::size_t> plain;
    for (std::size_t li = 0; li < loops.size(); ++li)
        if (std::find(winding.begin(), winding.end(), li) == winding.end())
            plain.push_back(li);
    // Across a winding, the loops' other parameter: the face lies to the
    // left of the loop that winds forwards, so on a torus the loop that winds
    // back goes within a turn on that side, and the holes between them.
    if (wrap >= 0 && periodic[1 - wrap] && winding.size() == 2) {
        const int d = wrap, e = 1 - wrap;
        BnLoop* pa = &loops[winding[0]];
        BnLoop* pb = &loops[winding[1]];
        if (pa->mWind[d] < 0)
            std::swap(pa, pb);
        // Forwards in the first parameter, the face is above; in the second,
        // to the right of it (counter-clockwise).
        const double base = mean(*pa, e);
        const double from = d == 0 ? base : base - kBnTwoPi;
        shift(*pb, e, -kBnTwoPi * std::floor((mean(*pb, e) - from) / kBnTwoPi));
        for (std::size_t li = 0; li < loops.size(); ++li)
            if (std::find(winding.begin(), winding.end(), li) == winding.end())
                shift(loops[li], e, -kBnTwoPi * std::floor((mean(loops[li], e) - from) / kBnTwoPi));
    }
    std::vector<BnLoop> out_loops;
    std::vector<bool> out_outer;
    // The left side of the strip a composed loop bounds, at the other
    // parameter `E`, for placing the plain loops in it.
    std::function<double(double)> strip_left;
    const auto seam_trim = [&](const double* pA, const double* pB, std::int64_t Seam, bool Along,
                               bool Pole) {
        BnTrim t;
        BnSample a, b;
        a.mX[0] = pA[0];
        a.mX[1] = pA[1];
        b.mX[0] = pB[0];
        b.mX[1] = pB[1];
        t.mS = {a, b};
        t.mSeam = Seam;
        t.mSeamAlong = Along;
        t.mPole = Pole;
        return t;
    };
    if (wrap >= 0) {
        const int d = wrap, e = 1 - wrap;
        if (winding.size() == 2) {
            BnLoop* pa = &loops[winding[0]];
            BnLoop* pb = &loops[winding[1]];
            if (pa->mWind[d] < 0)
                std::swap(pa, pb);
            if (pa->mWind[d] != 1 || pb->mWind[d] != -1) {
                rWhy = "two loops wind round the same way";
                return false;
            }
            const BnSample& a_end = pa->mTrims.back().mS.back();
            const BnSample& a_start = pa->mTrims.front().mS.front();
            const double b0 = pb->mTrims.front().mS.front().mX[d];
            shift(*pb, d, kBnTwoPi * std::round((a_end.mX[d] - b0) / kBnTwoPi));
            const BnSample& b_start = pb->mTrims.front().mS.front();
            const BnSample& b_end = pb->mTrims.back().mS.back();
            BnLoop composed;
            composed.mTrims = pa->mTrims;
            composed.mTrims.push_back(seam_trim(a_end.mX, b_start.mX, 0, true, false));
            composed.mTrims.insert(composed.mTrims.end(), pb->mTrims.begin(), pb->mTrims.end());
            composed.mTrims.push_back(seam_trim(b_end.mX, a_start.mX, 0, false, false));
            const double bx = b_end.mX[d], be = b_end.mX[e], ax = a_start.mX[d], ae = a_start.mX[e];
            strip_left = [=](double E) {
                return std::fabs(ae - be) > 0.0 ? bx + (ax - bx) * (E - be) / (ae - be)
                                                : std::min(ax, bx);
            };
            rOut.mSeams.emplace_back();
            rOut.mSeamHasCurve.push_back(true);
            out_loops.push_back(std::move(composed));
            out_outer.push_back(true);
        } else if (winding.size() == 1 && d == 0 && (kind == 2 || kind == 4)) {
            BnLoop& la = loops[winding[0]];
            const BnSample a_end = la.mTrims.back().mS.back();
            const BnSample a_start = la.mTrims.front().mS.front();
            // Through the pole (a sphere's) or the apex (a cone's) on the
            // side that leaves the face to the loops' left.
            std::vector<double> poles;
            if (kind == 2)
                poles = {0.5 * kBnPi, -0.5 * kBnPi};
            else
                poles = {0.0};
            bool found = false;
            for (double pole : poles) {
                BnLoop composed;
                composed.mTrims = la.mTrims;
                const double p1[2] = {a_end.mX[0], pole}, p2[2] = {a_start.mX[0], pole};
                composed.mTrims.push_back(seam_trim(a_end.mX, p1, 0, true, false));
                composed.mTrims.push_back(seam_trim(p1, p2, 1, true, true));
                composed.mTrims.push_back(seam_trim(p2, a_start.mX, 0, false, false));
                if (bn_area(composed, [](const BnSample& rS) {
                        return std::array<double, 2>{rS.mX[0], rS.mX[1]};
                    }) > 0.0) {
                    const double left = std::min(a_start.mX[0], a_end.mX[0]);
                    strip_left = [=](double) { return left; };
                    rOut.mSeams.emplace_back();
                    rOut.mSeamHasCurve.push_back(true);
                    rOut.mSeams.emplace_back();
                    rOut.mSeamHasCurve.push_back(false);
                    out_loops.push_back(std::move(composed));
                    out_outer.push_back(true);
                    found = true;
                    break;
                }
            }
            if (!found) {
                rWhy = "a loop round the axis bounds no side with a pole or apex";
                return false;
            }
        } else {
            rWhy = "loops wind round the surface in a way no seam closes";
            return false;
        }
        // The plain loops into the strip.
        for (std::size_t li : plain) {
            BnLoop& l = loops[li];
            double md = 0.0, me = 0.0;
            std::size_t count = 0;
            bn_each(l, [&](const BnSample& rS) {
                md += rS.mX[d];
                me += rS.mX[e];
                ++count;
            });
            md /= static_cast<double>(count);
            me /= static_cast<double>(count);
            const double left = strip_left(me);
            shift(l, d, -kBnTwoPi * std::floor((md - left) / kBnTwoPi));
            bool inside = true;
            bn_each(l, [&](const BnSample& rS) {
                const double x = strip_left(rS.mX[e]);
                inside = inside && rS.mX[d] >= x - 1e-12 && rS.mX[d] <= x + kBnTwoPi + 1e-12;
            });
            if (!inside) {
                rWhy = "a hole crosses the seam";
                return false;
            }
            out_loops.push_back(l);
            out_outer.push_back(false);
        }
    } else {
        // No winding: the one counter-clockwise loop is the outer one, and
        // each hole goes within half a turn of it.
        std::size_t outer = loops.size();
        for (std::size_t li = 0; li < loops.size(); ++li)
            if (bn_area(loops[li], raw_xy) > 0.0) {
                if (outer < loops.size()) {
                    rWhy = "the face has several outer loops";
                    return false;
                }
                outer = li;
            }
        if (outer == loops.size()) {
            rWhy = "the face has no outer loop";
            return false;
        }
        for (int k = 0; k < 2; ++k)
            if (periodic[k]) {
                const double at = mean(loops[outer], k);
                for (std::size_t li = 0; li < loops.size(); ++li)
                    if (li != outer)
                        shift(loops[li], k,
                              kBnTwoPi * std::round((at - mean(loops[li], k)) / kBnTwoPi));
            }
        for (BnLoop& l : loops) {
            out_loops.push_back(l);
            out_outer.push_back(false);
        }
    }

    // On a plane, a circular trim's centre in the raw parameters and its
    // angles, unwrapped along its samples: `{cu, cv, a0, a1}`.
    const auto circle_of = [&](const BnTrim& rT) {
        const PatchVec dc{rT.pSrc->mCentre[0] - rSurface.mOrigin[0],
                          rT.pSrc->mCentre[1] - rSurface.mOrigin[1],
                          rT.pSrc->mCentre[2] - rSurface.mOrigin[2]};
        const double cu = bn_dot(dc, rSurface.mX), cv = bn_dot(dc, rSurface.mY);
        double a0 = 0.0, a = 0.0;
        for (std::size_t i = 0; i < rT.mS.size(); ++i) {
            const double x = std::atan2(rT.mS[i].mX[1] - cv, rT.mS[i].mX[0] - cu);
            a = i == 0 ? x : bn_unwrap(x, a);
            if (i == 0)
                a0 = a;
        }
        return std::array<double, 4>{cu, cv, a0, a};
    };

    // The domain over every sample and, on a plane, every circular trim's
    // extreme points: the trim is the exact arc, which must lie within the
    // patch (evaluation clamps to it), not only its samples.
    double dlo[2] = {1e300, 1e300}, dhi[2] = {-1e300, -1e300};
    const auto take = [&](double U, double V) {
        dlo[0] = std::min(dlo[0], U);
        dhi[0] = std::max(dhi[0], U);
        dlo[1] = std::min(dlo[1], V);
        dhi[1] = std::max(dhi[1], V);
    };
    for (const BnLoop& l : out_loops)
        for (const BnTrim& t : l.mTrims) {
            for (const BnSample& s : t.mS)
                take(s.mX[0], s.mX[1]);
            if (t.pSrc == nullptr || kind != 1 || !t.pSrc->mCircle)
                continue;
            const auto c = circle_of(t);
            const double r = t.pSrc->mRadius, quarter = 0.5 * kBnPi;
            const auto k0 = static_cast<long long>(std::ceil(std::min(c[2], c[3]) / quarter));
            const auto k1 = static_cast<long long>(std::floor(std::max(c[2], c[3]) / quarter));
            for (long long k = k0; k <= k1; ++k) {
                // The point at angle k pi/2, exact on the axes.
                static constexpr double kCos[4] = {1.0, 0.0, -1.0, 0.0};
                static constexpr double kSin[4] = {0.0, 1.0, 0.0, -1.0};
                const auto q = static_cast<std::size_t>((k % 4 + 4) % 4);
                take(c[0] + r * kCos[q], c[1] + r * kSin[q]);
            }
        }
    for (int k = 0; k < 2; ++k)
        if (!(dhi[k] > dlo[k]) && kind != 6) {
            rWhy = "the loops span no area in the surface's parameters";
            return false;
        }

    // The patch, and the map from raw parameters to its own.
    NurbsPatch& patch = rOut.mPatch;
    UnitArc arc_u, arc_v;
    const auto put_grid = [&](std::size_t Nu, std::size_t Nv) {
        patch.mControlPoints = NDArray(DType::Float64, {Nu * Nv, 3});
        patch.mWeights = NDArray(DType::Float64, {Nu * Nv});
    };
    const auto set_pole = [&](std::size_t K, const PatchVec& rP, double W) {
        for (std::size_t j = 0; j < 3; ++j)
            patch.mControlPoints.As<double>()[K * 3 + j] = rP[j];
        patch.mWeights.As<double>()[K] = W;
    };
    const PatchVec& o = rSurface.mOrigin;
    const PatchVec& ax = rSurface.mX;
    const PatchVec& ay = rSurface.mY;
    const PatchVec& az = rSurface.mZ;
    if (kind == 1) {
        const double pad[2] = {1e-3 * (dhi[0] - dlo[0]) + tiny, 1e-3 * (dhi[1] - dlo[1]) + tiny};
        const double u0 = dlo[0] - pad[0], u1 = dhi[0] + pad[0];
        const double v0 = dlo[1] - pad[1], v1 = dhi[1] + pad[1];
        patch.mDegrees = {1, 1};
        patch.mKnots = {{u0, u0, u1, u1}, {v0, v0, v1, v1}};
        put_grid(2, 2);
        for (std::size_t j = 0; j < 2; ++j)
            for (std::size_t i = 0; i < 2; ++i) {
                const double u = i == 0 ? u0 : u1, v = j == 0 ? v0 : v1;
                set_pole(j * 2 + i,
                         {o[0] + u * ax[0] + v * ay[0], o[1] + u * ax[1] + v * ay[1],
                          o[2] + u * ax[2] + v * ay[2]},
                         1.0);
            }
    } else if (kind == 6) {
        patch.mDegrees = {bs->mDegreeU, bs->mDegreeV};
        patch.mKnots = {bs->mKnotsU, bs->mKnotsV};
        put_grid(bs->mNu, bs->mNv);
        for (std::size_t j = 0; j < bs->mNv; ++j)
            for (std::size_t i = 0; i < bs->mNu; ++i) {
                const double* p = &bs->mPoles[3 * (i + bs->mNu * j)];
                set_pole(j * bs->mNu + i, {p[0], p[1], p[2]}, 1.0);
            }
    } else {
        arc_u = unit_arc(dlo[0], dhi[0]);
        // The profile in the (radius, height) half-plane: its poles and weights.
        std::vector<std::array<double, 3>> profile;
        std::vector<double> profile_knots;
        int profile_degree = 1;
        if (kind == 3 || kind == 4) {
            const double t = kind == 3 ? 0.0 : rSurface.mTanHalfAngle;
            const double r0 = kind == 3 ? rSurface.mRadius : dlo[1] * t;
            const double r1 = kind == 3 ? rSurface.mRadius : dhi[1] * t;
            profile = {{r0, dlo[1], 1.0}, {r1, dhi[1], 1.0}};
            profile_knots = {dlo[1], dlo[1], dhi[1], dhi[1]};
        } else {
            arc_v = unit_arc(dlo[1], dhi[1]);
            profile_degree = 2;
            profile_knots = arc_v.mKnots;
            for (const auto& q : arc_v.mPoles) {
                if (kind == 2)
                    profile.push_back({rSurface.mRadius * q[0], rSurface.mRadius * q[1], q[2]});
                else
                    profile.push_back(
                        {rSurface.mRadius + rSurface.mMinor * q[0], rSurface.mMinor * q[1], q[2]});
            }
        }
        patch.mDegrees = {2, profile_degree};
        patch.mKnots = {arc_u.mKnots, profile_knots};
        const std::size_t nu = arc_u.mPoles.size(), nv = profile.size();
        put_grid(nu, nv);
        for (std::size_t j = 0; j < nv; ++j)
            for (std::size_t i = 0; i < nu; ++i) {
                const auto& c = arc_u.mPoles[i];
                const double rho = profile[j][0], z = profile[j][1];
                PatchVec p;
                for (std::size_t k = 0; k < 3; ++k)
                    p[k] = o[k] + rho * (c[0] * ax[k] + c[1] * ay[k]) + z * az[k];
                set_pole(j * nu + i, p, c[2] * profile[j][2]);
            }
    }
    const auto uv_of = [&](const BnSample& rS) -> std::array<double, 2> {
        if (!revolution)
            return {rS.mX[0], rS.mX[1]};
        const double u = unit_arc_param(arc_u, std::min(std::max(rS.mX[0], dlo[0]), dhi[0]));
        const double v = (kind == 3 || kind == 4)
                             ? rS.mX[1]
                             : unit_arc_param(arc_v, std::min(std::max(rS.mX[1], dlo[1]), dhi[1]));
        return {u, v};
    };

    // Trims with an edge curve refined until the surface between their
    // points follows it.
    const double follow = std::max(0.02 * Tol, 10.0 * tiny);
    for (BnLoop& l : out_loops)
        for (BnTrim& t : l.mTrims) {
            if (t.pSrc == nullptr || !t.pSrc->mAt || (kind == 1 && t.pSrc->mCircle))
                continue;
            for (int pass = 0; pass < 12 && t.mS.size() < 4096; ++pass) {
                std::vector<BnSample> next;
                bool changed = false;
                for (std::size_t i = 0; i + 1 < t.mS.size(); ++i) {
                    next.push_back(t.mS[i]);
                    const auto a = uv_of(t.mS[i]), b = uv_of(t.mS[i + 1]);
                    const PatchVec on =
                        nurbs_patch_point(patch, 0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1]));
                    BnSample m;
                    m.mT = 0.5 * (t.mS[i].mT + t.mS[i + 1].mT);
                    m.mP = bn_at(t, m.mT);
                    if (bn_dist(on, m.mP) <= follow)
                        continue;
                    raw(m.mP, &t.mS[i], m.mX);
                    if (!raw_ok) {
                        rWhy = "a trim runs through the surface's axis or pole";
                        return false;
                    }
                    next.push_back(m);
                    changed = true;
                }
                next.push_back(t.mS.back());
                t.mS = std::move(next);
                if (!changed)
                    break;
            }
        }

    // The orientation: a composed loop and a plain face's largest loop are
    // outer and run counter-clockwise; every other loop is a hole.
    std::vector<double> area(out_loops.size());
    for (std::size_t i = 0; i < out_loops.size(); ++i)
        area[i] = bn_area(out_loops[i], uv_of);
    if (wrap < 0) {
        std::size_t outer = 0;
        for (std::size_t i = 1; i < area.size(); ++i)
            if (std::fabs(area[i]) > std::fabs(area[outer]))
                outer = i;
        out_outer[outer] = true;
    }
    for (std::size_t i = 0; i < out_loops.size(); ++i)
        if (out_outer[i] ? !(area[i] > 0.0) : !(area[i] < 0.0)) {
            rWhy = "the loops are not one outer loop round its holes";
            return false;
        }

    // The seams' 3-D curves: the surface along their first trim.
    std::vector<bool> seam_done(rOut.mSeams.size(), false);
    for (const BnLoop& l : out_loops)
        for (const BnTrim& t : l.mTrims) {
            if (t.mSeam < 0 || !t.mSeamAlong || t.mPole)
                continue;
            const std::size_t k = static_cast<std::size_t>(t.mSeam);
            if (seam_done[k])
                continue;
            seam_done[k] = true;
            const auto a = uv_of(t.mS.front()), b = uv_of(t.mS.back());
            std::vector<PatchVec> pts;
            for (int i = 0; i <= 64; ++i) {
                const double s = static_cast<double>(i) / 64.0;
                pts.push_back(
                    nurbs_patch_point(patch, a[0] + (b[0] - a[0]) * s, a[1] + (b[1] - a[1]) * s));
            }
            CurveFitOptions co;
            co.mTolerance = std::max(0.05 * Tol, 10.0 * tiny);
            const CurveFitResult fit = bspline_curve_fit(pts, co);
            rOut.mSeams[k] = fit.mOk ? brep_curve_bspline(fit.mCurve) : brep_curve_polyline(pts);
        }

    // Every loop point within the tolerance of the surface at its parameters.
    for (const BnLoop& l : out_loops)
        for (const BnTrim& t : l.mTrims) {
            if (t.pSrc == nullptr)
                continue;
            for (const BnSample& s : t.mS) {
                const auto uv = uv_of(s);
                const double off = bn_dist(nurbs_patch_point(patch, uv[0], uv[1]), s.mP);
                if (off > 4.0 * Tol + 10.0 * tiny) {
                    rWhy = "a loop point lies " + std::to_string(off) + " off the " +
                           (kind == 6 ? "B-spline" : "analytic") +
                           " surface at its parameters (tolerance " + std::to_string(Tol) + ")";
                    return false;
                }
            }
        }

    // The face: trims as parameter-space curves along their loops.
    BrepFace& face = rOut.mFace;
    face.mRational = revolution;
    face.mSwappedNormal = !Sense;
    face.mTrimmed = true;
    std::int64_t trim_index = 0;
    for (std::size_t li = 0; li < out_loops.size(); ++li) {
        BrepLoop loop;
        loop.mType = out_outer[li] ? "outer" : "inner";
        for (const BnTrim& t : out_loops[li].mTrims) {
            BrepTrim bt;
            bt.mTrimIndex = trim_index++;
            bt.mDirection = true;
            if (t.pSrc != nullptr && kind == 1 && t.pSrc->mCircle) {
                // The exact arc, its angles unwrapped along the samples.
                const auto c = circle_of(t);
                bt.mCurve = brep_curve_arc({c[0], c[1], 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0},
                                           t.pSrc->mRadius, c[2], c[3]);
            } else {
                std::vector<PatchVec> uv;
                for (const BnSample& s : t.mS) {
                    const auto q = uv_of(s);
                    if (uv.empty() || q[0] != uv.back()[0] || q[1] != uv.back()[1])
                        uv.push_back({q[0], q[1], 0.0});
                }
                if (uv.size() < 2)
                    uv.push_back(uv.back());
                // A trim along a parameter line is that line.
                const double span_u = patch.mKnots[0].back() - patch.mKnots[0].front();
                const double span_v = patch.mKnots[1].back() - patch.mKnots[1].front();
                bool iso_u = true, iso_v = true;
                for (const PatchVec& q : uv) {
                    iso_u = iso_u && std::fabs(q[0] - uv.front()[0]) <= 1e-12 * span_u;
                    iso_v = iso_v && std::fabs(q[1] - uv.front()[1]) <= 1e-12 * span_v;
                }
                if (iso_u || iso_v)
                    uv = {uv.front(), uv.back()};
                bt.mCurve = brep_curve_polyline(uv);
            }
            BrepNurbsFace::Ref ref;
            if (t.pSrc != nullptr) {
                ref.mEdge = t.pSrc->mEdge;
                ref.mAlong = t.pSrc->mAlong != t.mRev;
            } else {
                ref.mSeam = t.mSeam;
                ref.mAlong = t.mSeamAlong;
            }
            rOut.mRefs.push_back(ref);
            loop.mTrims.push_back(std::move(bt));
        }
        face.mLoops.push_back(std::move(loop));
    }
    // The outer loop first.
    for (std::size_t li = 0; li < face.mLoops.size(); ++li)
        if (face.mLoops[li].mType == "outer" && li > 0) {
            std::size_t before = 0;
            for (std::size_t k = 0; k < li; ++k)
                before += face.mLoops[k].mTrims.size();
            const std::size_t count = face.mLoops[li].mTrims.size();
            std::rotate(face.mLoops.begin(), face.mLoops.begin() + static_cast<std::ptrdiff_t>(li),
                        face.mLoops.begin() + static_cast<std::ptrdiff_t>(li) + 1);
            std::rotate(rOut.mRefs.begin(),
                        rOut.mRefs.begin() + static_cast<std::ptrdiff_t>(before),
                        rOut.mRefs.begin() + static_cast<std::ptrdiff_t>(before + count));
            std::int64_t ti = 0;
            for (BrepLoop& l : face.mLoops)
                for (BrepTrim& t : l.mTrims)
                    t.mTrimIndex = ti++;
            break;
        }
    return true;
}

}  // namespace detail
}  // namespace meshioplusplus

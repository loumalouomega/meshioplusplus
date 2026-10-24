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
// The robust predicates of detail/predicates.hpp: a floating-point filter
// with Shewchuk's stage-A error bounds, then exact expansion arithmetic.
// Anonymous-namespace helpers are prefixed `pred_` (the amalgamation
// concatenates every translation unit into one).

// System includes
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

// Project includes
#include "predicates.hpp"

// The error-free transforms below are exact only when every product is
// rounded on its own: switch floating-point contraction (fused multiply-add)
// off for this file, whatever flags the including translation unit uses.
// MSVC does not contract under its default /fp:precise.
#if defined(__clang__)
#define MESHIOPLUSPLUS_PRED_NO_CONTRACT _Pragma("clang fp contract(off)")
#elif defined(__GNUC__)
#pragma GCC push_options
#pragma GCC optimize("fp-contract=off")
#define MESHIOPLUSPLUS_PRED_NO_CONTRACT
#else
#define MESHIOPLUSPLUS_PRED_NO_CONTRACT
#endif

namespace meshioplusplus {
namespace detail {

namespace {

static_assert(std::numeric_limits<double>::is_iec559 && std::numeric_limits<double>::digits == 53,
              "the robust predicates need IEEE-754 binary64");

/// Half an ulp of 1 (2^-53) and Dekker's splitter (2^27 + 1).
constexpr double kPredEps = 1.1102230246251565404236316680908203125e-16;
constexpr double kPredSplitter = 134217729.0;

/// Shewchuk's stage-A bounds, relative to each determinant's permanent.
constexpr double kPredCcwBound = (3.0 + 16.0 * kPredEps) * kPredEps;
constexpr double kPredIccBound = (10.0 + 96.0 * kPredEps) * kPredEps;
constexpr double kPredO3dBound = (7.0 + 56.0 * kPredEps) * kPredEps;
constexpr double kPredIspBound = (16.0 + 224.0 * kPredEps) * kPredEps;

/// An expansion: nonoverlapping components in increasing magnitude, zeros
/// eliminated; empty is zero.
using PredExpansion = std::vector<double>;

inline void pred_two_sum(double A, double B, double& rX, double& rY) {
    MESHIOPLUSPLUS_PRED_NO_CONTRACT
    rX = A + B;
    const double bv = rX - A;
    const double av = rX - bv;
    const double br = B - bv;
    const double ar = A - av;
    rY = ar + br;
}

inline void pred_two_diff(double A, double B, double& rX, double& rY) {
    MESHIOPLUSPLUS_PRED_NO_CONTRACT
    rX = A - B;
    const double bv = A - rX;
    const double av = rX + bv;
    const double br = bv - B;
    const double ar = A - av;
    rY = ar + br;
}

inline void pred_split(double A, double& rHi, double& rLo) {
    MESHIOPLUSPLUS_PRED_NO_CONTRACT
    const double c = kPredSplitter * A;
    const double abig = c - A;
    rHi = c - abig;
    rLo = A - rHi;
}

inline void pred_two_product(double A, double B, double& rX, double& rY) {
    MESHIOPLUSPLUS_PRED_NO_CONTRACT
    rX = A * B;
    double ahi = 0.0;
    double alo = 0.0;
    double bhi = 0.0;
    double blo = 0.0;
    pred_split(A, ahi, alo);
    pred_split(B, bhi, blo);
    const double err1 = rX - (ahi * bhi);
    const double err2 = err1 - (alo * bhi);
    const double err3 = err2 - (ahi * blo);
    rY = (alo * blo) - err3;
}

/// `A - B` exactly, as an expansion.
PredExpansion pred_diff(double A, double B) {
    double x = 0.0;
    double y = 0.0;
    pred_two_diff(A, B, x, y);
    PredExpansion e;
    if (y != 0.0)
        e.push_back(y);
    if (x != 0.0)
        e.push_back(x);
    return e;
}

/// `e + f` (Shewchuk's FAST-EXPANSION-SUM with zero elimination): merge by
/// magnitude, then a chain of two-sums.
PredExpansion pred_sum(const PredExpansion& rE, const PredExpansion& rF) {
    if (rE.empty())
        return rF;
    if (rF.empty())
        return rE;
    std::vector<double> g;
    g.reserve(rE.size() + rF.size());
    std::size_t i = 0;
    std::size_t j = 0;
    while (i < rE.size() && j < rF.size()) {
        if (std::fabs(rE[i]) <= std::fabs(rF[j]))
            g.push_back(rE[i++]);
        else
            g.push_back(rF[j++]);
    }
    while (i < rE.size())
        g.push_back(rE[i++]);
    while (j < rF.size())
        g.push_back(rF[j++]);
    PredExpansion h;
    h.reserve(g.size());
    double q = g[0];
    for (std::size_t k = 1; k < g.size(); ++k) {
        double qnew = 0.0;
        double hh = 0.0;
        pred_two_sum(q, g[k], qnew, hh);
        if (hh != 0.0)
            h.push_back(hh);
        q = qnew;
    }
    if (q != 0.0)
        h.push_back(q);
    return h;
}

PredExpansion pred_negate(PredExpansion E) {
    for (double& r_c : E)
        r_c = -r_c;
    return E;
}

PredExpansion pred_sub(const PredExpansion& rE, const PredExpansion& rF) {
    return pred_sum(rE, pred_negate(rF));
}

/// `e * b` (SCALE-EXPANSION with zero elimination).
PredExpansion pred_scale(const PredExpansion& rE, double B) {
    PredExpansion h;
    if (rE.empty() || B == 0.0)
        return h;
    h.reserve(2 * rE.size());
    double q = 0.0;
    double hh = 0.0;
    pred_two_product(rE[0], B, q, hh);
    if (hh != 0.0)
        h.push_back(hh);
    for (std::size_t i = 1; i < rE.size(); ++i) {
        double p1 = 0.0;
        double p0 = 0.0;
        pred_two_product(rE[i], B, p1, p0);
        double sum = 0.0;
        pred_two_sum(q, p0, sum, hh);
        if (hh != 0.0)
            h.push_back(hh);
        pred_two_sum(p1, sum, q, hh);
        if (hh != 0.0)
            h.push_back(hh);
    }
    if (q != 0.0)
        h.push_back(q);
    return h;
}

/// `e * f`, one scaling per component of `f`.
PredExpansion pred_mul(const PredExpansion& rE, const PredExpansion& rF) {
    PredExpansion acc;
    for (const double c : rF)
        acc = pred_sum(acc, pred_scale(rE, c));
    return acc;
}

/// The expansion's value rounded to a double: its sign is the sign of the
/// largest component, the exact sign of the expansion.
double pred_estimate(const PredExpansion& rE) {
    double s = 0.0;
    for (const double c : rE)
        s += c;
    if (s == 0.0 && !rE.empty())
        return rE.back();
    return s;
}

/// `a * d - b * c` of four expansions.
PredExpansion pred_minor(const PredExpansion& rA, const PredExpansion& rD, const PredExpansion& rB,
                         const PredExpansion& rC) {
    return pred_sub(pred_mul(rA, rD), pred_mul(rB, rC));
}

double pred_orient2d_exact(const double* pA, const double* pB, const double* pC) {
    const PredExpansion acx = pred_diff(pA[0], pC[0]);
    const PredExpansion acy = pred_diff(pA[1], pC[1]);
    const PredExpansion bcx = pred_diff(pB[0], pC[0]);
    const PredExpansion bcy = pred_diff(pB[1], pC[1]);
    return pred_estimate(pred_minor(acx, bcy, acy, bcx));
}

double pred_incircle_exact(const double* pA, const double* pB, const double* pC, const double* pD) {
    const PredExpansion adx = pred_diff(pA[0], pD[0]);
    const PredExpansion ady = pred_diff(pA[1], pD[1]);
    const PredExpansion bdx = pred_diff(pB[0], pD[0]);
    const PredExpansion bdy = pred_diff(pB[1], pD[1]);
    const PredExpansion cdx = pred_diff(pC[0], pD[0]);
    const PredExpansion cdy = pred_diff(pC[1], pD[1]);
    const PredExpansion alift = pred_sum(pred_mul(adx, adx), pred_mul(ady, ady));
    const PredExpansion blift = pred_sum(pred_mul(bdx, bdx), pred_mul(bdy, bdy));
    const PredExpansion clift = pred_sum(pred_mul(cdx, cdx), pred_mul(cdy, cdy));
    PredExpansion det = pred_mul(alift, pred_minor(bdx, cdy, cdx, bdy));
    det = pred_sum(det, pred_mul(blift, pred_minor(cdx, ady, adx, cdy)));
    det = pred_sum(det, pred_mul(clift, pred_minor(adx, bdy, bdx, ady)));
    return pred_estimate(det);
}

double pred_orient3d_exact(const double* pA, const double* pB, const double* pC, const double* pD) {
    const PredExpansion adx = pred_diff(pA[0], pD[0]);
    const PredExpansion ady = pred_diff(pA[1], pD[1]);
    const PredExpansion adz = pred_diff(pA[2], pD[2]);
    const PredExpansion bdx = pred_diff(pB[0], pD[0]);
    const PredExpansion bdy = pred_diff(pB[1], pD[1]);
    const PredExpansion bdz = pred_diff(pB[2], pD[2]);
    const PredExpansion cdx = pred_diff(pC[0], pD[0]);
    const PredExpansion cdy = pred_diff(pC[1], pD[1]);
    const PredExpansion cdz = pred_diff(pC[2], pD[2]);
    PredExpansion det = pred_mul(adz, pred_minor(bdx, cdy, cdx, bdy));
    det = pred_sum(det, pred_mul(bdz, pred_minor(cdx, ady, adx, cdy)));
    det = pred_sum(det, pred_mul(cdz, pred_minor(adx, bdy, bdx, ady)));
    return pred_estimate(det);
}

double pred_insphere_exact(const double* pA, const double* pB, const double* pC, const double* pD,
                           const double* pE) {
    PredExpansion d[4][3];
    const double* pts[4] = {pA, pB, pC, pD};
    for (int i = 0; i < 4; ++i)
        for (int k = 0; k < 3; ++k)
            d[i][k] = pred_diff(pts[i][k], pE[k]);
    // The 2x2 minors of the (x, y) columns, then the 3x3 cofactors and the
    // lifted determinant in the paper's arrangement.
    const auto xy = [&](int I, int J) { return pred_minor(d[I][0], d[J][1], d[J][0], d[I][1]); };
    const PredExpansion ab = xy(0, 1);
    const PredExpansion bc = xy(1, 2);
    const PredExpansion cd = xy(2, 3);
    const PredExpansion da = xy(3, 0);
    const PredExpansion ac = xy(0, 2);
    const PredExpansion bd = xy(1, 3);
    const PredExpansion abc =
        pred_sum(pred_sub(pred_mul(d[0][2], bc), pred_mul(d[1][2], ac)), pred_mul(d[2][2], ab));
    const PredExpansion bcd =
        pred_sum(pred_sub(pred_mul(d[1][2], cd), pred_mul(d[2][2], bd)), pred_mul(d[3][2], bc));
    const PredExpansion cda =
        pred_sum(pred_sum(pred_mul(d[2][2], da), pred_mul(d[3][2], ac)), pred_mul(d[0][2], cd));
    const PredExpansion dab =
        pred_sum(pred_sum(pred_mul(d[3][2], ab), pred_mul(d[0][2], bd)), pred_mul(d[1][2], da));
    PredExpansion lift[4];
    for (int i = 0; i < 4; ++i)
        lift[i] = pred_sum(pred_sum(pred_mul(d[i][0], d[i][0]), pred_mul(d[i][1], d[i][1])),
                           pred_mul(d[i][2], d[i][2]));
    const PredExpansion det = pred_sum(pred_sub(pred_mul(lift[3], abc), pred_mul(lift[2], dab)),
                                       pred_sub(pred_mul(lift[1], cda), pred_mul(lift[0], bcd)));
    return pred_estimate(det);
}

/// `(a - p) . (b - p)` and `(a - p) x (b - p)` exactly.
void pred_dot_cross_exact(const double* pA, const double* pB, const double* pP, PredExpansion& rDot,
                          PredExpansion& rCross) {
    const PredExpansion apx = pred_diff(pA[0], pP[0]);
    const PredExpansion apy = pred_diff(pA[1], pP[1]);
    const PredExpansion bpx = pred_diff(pB[0], pP[0]);
    const PredExpansion bpy = pred_diff(pB[1], pP[1]);
    rDot = pred_sum(pred_mul(apx, bpx), pred_mul(apy, bpy));
    rCross = pred_minor(apx, bpy, apy, bpx);
}

double pred_diametral_exact(const double* pA, const double* pB, const double* pP) {
    PredExpansion dot;
    PredExpansion cross;
    pred_dot_cross_exact(pA, pB, pP, dot, cross);
    return pred_estimate(dot);
}

double pred_lens_angle_exact(const double* pA, const double* pB, const double* pP) {
    PredExpansion dot;
    PredExpansion cross;
    pred_dot_cross_exact(pA, pB, pP, dot, cross);
    return pred_estimate(pred_sub(pred_scale(pred_mul(dot, dot), 3.0), pred_mul(cross, cross)));
}

}  // namespace

void predicates_two_product(double A, double B, double& rX, double& rY) {
    pred_two_product(A, B, rX, rY);
}

double orient2d(const double* pA, const double* pB, const double* pC) {
    MESHIOPLUSPLUS_PRED_NO_CONTRACT
    const double detleft = (pA[0] - pC[0]) * (pB[1] - pC[1]);
    const double detright = (pA[1] - pC[1]) * (pB[0] - pC[0]);
    const double det = detleft - detright;
    double detsum = 0.0;
    if (detleft > 0.0) {
        if (detright <= 0.0)
            return det;
        detsum = detleft + detright;
    } else if (detleft < 0.0) {
        if (detright >= 0.0)
            return det;
        detsum = -detleft - detright;
    } else {
        return det;
    }
    const double errbound = kPredCcwBound * detsum;
    if (det >= errbound || -det >= errbound)
        return det;
    return pred_orient2d_exact(pA, pB, pC);
}

double incircle(const double* pA, const double* pB, const double* pC, const double* pD) {
    MESHIOPLUSPLUS_PRED_NO_CONTRACT
    const double adx = pA[0] - pD[0];
    const double bdx = pB[0] - pD[0];
    const double cdx = pC[0] - pD[0];
    const double ady = pA[1] - pD[1];
    const double bdy = pB[1] - pD[1];
    const double cdy = pC[1] - pD[1];
    const double bdxcdy = bdx * cdy;
    const double cdxbdy = cdx * bdy;
    const double alift = adx * adx + ady * ady;
    const double cdxady = cdx * ady;
    const double adxcdy = adx * cdy;
    const double blift = bdx * bdx + bdy * bdy;
    const double adxbdy = adx * bdy;
    const double bdxady = bdx * ady;
    const double clift = cdx * cdx + cdy * cdy;
    const double det =
        alift * (bdxcdy - cdxbdy) + blift * (cdxady - adxcdy) + clift * (adxbdy - bdxady);
    const double permanent = (std::fabs(bdxcdy) + std::fabs(cdxbdy)) * alift +
                             (std::fabs(cdxady) + std::fabs(adxcdy)) * blift +
                             (std::fabs(adxbdy) + std::fabs(bdxady)) * clift;
    const double errbound = kPredIccBound * permanent;
    if (det > errbound || -det > errbound)
        return det;
    return pred_incircle_exact(pA, pB, pC, pD);
}

double orient3d(const double* pA, const double* pB, const double* pC, const double* pD) {
    MESHIOPLUSPLUS_PRED_NO_CONTRACT
    const double adx = pA[0] - pD[0];
    const double bdx = pB[0] - pD[0];
    const double cdx = pC[0] - pD[0];
    const double ady = pA[1] - pD[1];
    const double bdy = pB[1] - pD[1];
    const double cdy = pC[1] - pD[1];
    const double adz = pA[2] - pD[2];
    const double bdz = pB[2] - pD[2];
    const double cdz = pC[2] - pD[2];
    const double bdxcdy = bdx * cdy;
    const double cdxbdy = cdx * bdy;
    const double cdxady = cdx * ady;
    const double adxcdy = adx * cdy;
    const double adxbdy = adx * bdy;
    const double bdxady = bdx * ady;
    const double det = adz * (bdxcdy - cdxbdy) + bdz * (cdxady - adxcdy) + cdz * (adxbdy - bdxady);
    const double permanent = (std::fabs(bdxcdy) + std::fabs(cdxbdy)) * std::fabs(adz) +
                             (std::fabs(cdxady) + std::fabs(adxcdy)) * std::fabs(bdz) +
                             (std::fabs(adxbdy) + std::fabs(bdxady)) * std::fabs(cdz);
    const double errbound = kPredO3dBound * permanent;
    if (det > errbound || -det > errbound)
        return det;
    return pred_orient3d_exact(pA, pB, pC, pD);
}

double insphere(const double* pA, const double* pB, const double* pC, const double* pD,
                const double* pE) {
    MESHIOPLUSPLUS_PRED_NO_CONTRACT
    const double aex = pA[0] - pE[0];
    const double bex = pB[0] - pE[0];
    const double cex = pC[0] - pE[0];
    const double dex = pD[0] - pE[0];
    const double aey = pA[1] - pE[1];
    const double bey = pB[1] - pE[1];
    const double cey = pC[1] - pE[1];
    const double dey = pD[1] - pE[1];
    const double aez = pA[2] - pE[2];
    const double bez = pB[2] - pE[2];
    const double cez = pC[2] - pE[2];
    const double dez = pD[2] - pE[2];

    const double aexbey = aex * bey;
    const double bexaey = bex * aey;
    const double ab = aexbey - bexaey;
    const double bexcey = bex * cey;
    const double cexbey = cex * bey;
    const double bc = bexcey - cexbey;
    const double cexdey = cex * dey;
    const double dexcey = dex * cey;
    const double cd = cexdey - dexcey;
    const double dexaey = dex * aey;
    const double aexdey = aex * dey;
    const double da = dexaey - aexdey;
    const double aexcey = aex * cey;
    const double cexaey = cex * aey;
    const double ac = aexcey - cexaey;
    const double bexdey = bex * dey;
    const double dexbey = dex * bey;
    const double bd = bexdey - dexbey;

    const double abc = aez * bc - bez * ac + cez * ab;
    const double bcd = bez * cd - cez * bd + dez * bc;
    const double cda = cez * da + dez * ac + aez * cd;
    const double dab = dez * ab + aez * bd + bez * da;

    const double alift = aex * aex + aey * aey + aez * aez;
    const double blift = bex * bex + bey * bey + bez * bez;
    const double clift = cex * cex + cey * cey + cez * cez;
    const double dlift = dex * dex + dey * dey + dez * dez;

    const double det = (dlift * abc - clift * dab) + (blift * cda - alift * bcd);

    const double ab_p = std::fabs(aexbey) + std::fabs(bexaey);
    const double bc_p = std::fabs(bexcey) + std::fabs(cexbey);
    const double cd_p = std::fabs(cexdey) + std::fabs(dexcey);
    const double da_p = std::fabs(dexaey) + std::fabs(aexdey);
    const double ac_p = std::fabs(aexcey) + std::fabs(cexaey);
    const double bd_p = std::fabs(bexdey) + std::fabs(dexbey);
    const double aezp = std::fabs(aez);
    const double bezp = std::fabs(bez);
    const double cezp = std::fabs(cez);
    const double dezp = std::fabs(dez);
    const double permanent = (aezp * bc_p + bezp * ac_p + cezp * ab_p) * dlift +
                             (dezp * ab_p + aezp * bd_p + bezp * da_p) * clift +
                             (cezp * da_p + dezp * ac_p + aezp * cd_p) * blift +
                             (bezp * cd_p + cezp * bd_p + dezp * bc_p) * alift;
    // Twice the paper's bound: conservative, the exact stage decides the rest.
    const double errbound = 2.0 * kPredIspBound * permanent;
    if (det > errbound || -det > errbound)
        return det;
    return pred_insphere_exact(pA, pB, pC, pD, pE);
}

double diametral(const double* pA, const double* pB, const double* pP) {
    MESHIOPLUSPLUS_PRED_NO_CONTRACT
    // orient2d's shape (two products of differences, summed): its bound.
    const double xx = (pA[0] - pP[0]) * (pB[0] - pP[0]);
    const double yy = (pA[1] - pP[1]) * (pB[1] - pP[1]);
    const double dot = xx + yy;
    const double errbound = kPredCcwBound * (std::fabs(xx) + std::fabs(yy));
    if (dot > errbound || -dot > errbound)
        return dot;
    return pred_diametral_exact(pA, pB, pP);
}

double lens_angle(const double* pA, const double* pB, const double* pP) {
    MESHIOPLUSPLUS_PRED_NO_CONTRACT
    const double apx = pA[0] - pP[0];
    const double apy = pA[1] - pP[1];
    const double bpx = pB[0] - pP[0];
    const double bpy = pB[1] - pP[1];
    const double xx = apx * bpx;
    const double yy = apy * bpy;
    const double xy = apx * bpy;
    const double yx = apy * bpx;
    const double dot = xx + yy;
    const double cross = xy - yx;
    const double det = 3.0 * (dot * dot) - cross * cross;
    // The dot and cross products are each within (3 + 16 eps) eps of their
    // permanent P (orient2d's bound), so their squares are within about
    // 6 eps P^2, the weighted sum within 24 eps P^2 and its roundings
    // within 8 eps P^2: 64 eps P^2 bounds the error with room to spare.
    const double permanent = std::fabs(xx) + std::fabs(yy) + std::fabs(xy) + std::fabs(yx);
    const double errbound = 64.0 * kPredEps * (permanent * permanent);
    if (det > errbound || -det > errbound)
        return det;
    return pred_lens_angle_exact(pA, pB, pP);
}

}  // namespace detail
}  // namespace meshioplusplus

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC pop_options
#endif
#undef MESHIOPLUSPLUS_PRED_NO_CONTRACT

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
// Per-cell mesh quality metrics + inverted/degenerate detection. Every metric
// is evaluated on the cell's corner nodes (quadratic variants reduce to their
// linear parent) using only standard C++ math and the uniform mesh API, so it
// runs under every mesh backend. Metric formulas follow the common Verdict/VTK
// conventions; the ideal (unit) values are pinned by tests.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/quality.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/cell_faces.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {

namespace {

using detail::Vec3;

constexpr int NUM_METRICS = 11;
// Fixed metric order (index === position in QualityReport::mMetrics).
enum QualityMetric {
    QM_VOLUME = 0,
    QM_SCALED_JACOBIAN,
    QM_ASPECT_RATIO,
    QM_SKEWNESS,
    QM_MIN_ANGLE,
    QM_MAX_ANGLE,
    QM_WARPAGE,
    QM_MIN_DIHEDRAL,
    QM_MAX_DIHEDRAL,
    QM_INVERTED,
    QM_DEGENERATE,
};

const char* const QUALITY_METRIC_NAMES[NUM_METRICS] = {
    "quality:volume",   "quality:scaled_jacobian", "quality:aspect_ratio",
    "quality:skewness", "quality:min_angle",       "quality:max_angle",
    "quality:warpage",  "quality:min_dihedral",    "quality:max_dihedral",
    "quality:inverted", "quality:degenerate",
};

constexpr double QUALITY_NAN = std::numeric_limits<double>::quiet_NaN();
constexpr double QUALITY_RAD2DEG = 57.29577951308232087680;  // 180/pi

enum class QualityFamily { None, Tri, Quad, Tetra, Hex, Wedge, Pyr };

QualityFamily quality_family(CellType ct) {
    const int cc = detail::cell_corner_count(ct);
    const int dim = cell_type_dimension(ct);
    if (dim == 2) {
        if (cc == 3)
            return QualityFamily::Tri;
        if (cc == 4)
            return QualityFamily::Quad;
        return QualityFamily::None;
    }
    if (dim == 3) {
        // Only the types with a known face table (tetra(10), hex(20/27),
        // wedge(15), pyramid(13/14)); exotic high-order variants are skipped.
        if (!detail::skin_supported(ct))
            return QualityFamily::None;
        if (cc == 4)
            return QualityFamily::Tetra;
        if (cc == 8)
            return QualityFamily::Hex;
        if (cc == 6)
            return QualityFamily::Wedge;
        if (cc == 5)
            return QualityFamily::Pyr;
    }
    return QualityFamily::None;
}

// --- small geometry helpers -------------------------------------------------

double quality_edge_len(const Vec3& a, const Vec3& b) {
    return detail::vec3_norm(detail::vec3_sub(b, a));
}

// Interior angle (degrees) at the vertex from which edges u and v emanate.
double quality_angle_deg(const Vec3& u, const Vec3& v, double eps) {
    const double nu = detail::vec3_norm(u);
    const double nv = detail::vec3_norm(v);
    if (nu < eps || nv < eps)
        return QUALITY_NAN;
    double c = detail::vec3_dot(u, v) / (nu * nv);
    c = std::clamp(c, -1.0, 1.0);
    return std::acos(c) * QUALITY_RAD2DEG;
}

// det of the three unit edge vectors; NaN if any is degenerate.
double quality_det_unit3(const Vec3& a, const Vec3& b, const Vec3& c, double eps) {
    const double na = detail::vec3_norm(a);
    const double nb = detail::vec3_norm(b);
    const double nc = detail::vec3_norm(c);
    if (na < eps || nb < eps || nc < eps)
        return QUALITY_NAN;
    return detail::det3(detail::vec3_scale(a, 1.0 / na), detail::vec3_scale(b, 1.0 / nb),
                        detail::vec3_scale(c, 1.0 / nc));
}

// Signed volume via a cell-centroid / face-centroid fan over outward-wound
// faces (robust to non-planar faces). Positive for a well-oriented cell.
double quality_facefan_volume(const std::vector<Vec3>& rCorners,
                              const std::vector<detail::CellFaceDef>& rFaces) {
    Vec3 cc = {0.0, 0.0, 0.0};
    for (const Vec3& p : rCorners)
        cc = detail::vec3_add(cc, p);
    cc = detail::vec3_scale(cc, 1.0 / static_cast<double>(rCorners.size()));
    double vol = 0.0;
    for (const detail::CellFaceDef& f : rFaces) {
        Vec3 fc = {0.0, 0.0, 0.0};
        for (std::uint8_t k = 0; k < f.mNumCorners; ++k)
            fc = detail::vec3_add(fc, rCorners[f.mNodes[k]]);
        fc = detail::vec3_scale(fc, 1.0 / static_cast<double>(f.mNumCorners));
        for (std::uint8_t k = 0; k < f.mNumCorners; ++k) {
            const Vec3& a = rCorners[f.mNodes[k]];
            const Vec3& b = rCorners[f.mNodes[(k + 1) % f.mNumCorners]];
            vol += detail::triple_product(detail::vec3_sub(a, cc), detail::vec3_sub(b, cc),
                                          detail::vec3_sub(fc, cc)) /
                   6.0;
        }
    }
    return vol;
}

// Equiangle skewness over a set of quad-face corner angles (degrees). 0 = ideal.
double quality_equiangle_skew(double min_angle, double max_angle) {
    return std::max((max_angle - 90.0) / 90.0, (90.0 - min_angle) / 90.0);
}

// --- per-family metric evaluation (writes into m[NUM_METRICS], pre-NaN) -----

void quality_tri(const std::vector<Vec3>& c, bool is2d, double eps, double* m) {
    const Vec3& p0 = c[0];
    const Vec3& p1 = c[1];
    const Vec3& p2 = c[2];
    const double a = quality_edge_len(p0, p1);
    const double b = quality_edge_len(p1, p2);
    const double cc = quality_edge_len(p2, p0);

    double area;
    if (is2d)
        area = 0.5 * ((p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1]));
    else
        area = 0.5 * detail::vec3_norm(
                         detail::vec3_cross(detail::vec3_sub(p1, p0), detail::vec3_sub(p2, p0)));
    m[QM_VOLUME] = area;
    const double abs_area = std::fabs(area);
    const bool degenerate = abs_area < eps || a < eps || b < eps || cc < eps;
    m[QM_INVERTED] = (is2d && area < 0.0) ? 1.0 : 0.0;
    m[QM_DEGENERATE] = degenerate ? 1.0 : 0.0;
    if (degenerate)
        return;

    const double a0 = quality_angle_deg(detail::vec3_sub(p1, p0), detail::vec3_sub(p2, p0), eps);
    const double a1 = quality_angle_deg(detail::vec3_sub(p0, p1), detail::vec3_sub(p2, p1), eps);
    const double a2 = quality_angle_deg(detail::vec3_sub(p0, p2), detail::vec3_sub(p1, p2), eps);
    m[QM_MIN_ANGLE] = std::min({a0, a1, a2});
    m[QM_MAX_ANGLE] = std::max({a0, a1, a2});
    // Radius ratio R/(2r) = a*b*c*(a+b+c) / (16 A^2); 1 for an equilateral tri.
    m[QM_ASPECT_RATIO] = a * b * cc * (a + b + cc) / (16.0 * abs_area * abs_area);
}

void quality_quad(const std::vector<Vec3>& c, bool is2d, double eps, double* m) {
    const Vec3& p0 = c[0];
    const Vec3& p1 = c[1];
    const Vec3& p2 = c[2];
    const Vec3& p3 = c[3];
    const double l0 = quality_edge_len(p0, p1);
    const double l1 = quality_edge_len(p1, p2);
    const double l2 = quality_edge_len(p2, p3);
    const double l3 = quality_edge_len(p3, p0);

    double area;
    if (is2d) {
        area = 0.5 * ((p0[0] * p1[1] - p1[0] * p0[1]) + (p1[0] * p2[1] - p2[0] * p1[1]) +
                      (p2[0] * p3[1] - p3[0] * p2[1]) + (p3[0] * p0[1] - p0[0] * p3[1]));
    } else {
        Vec3 n = {0.0, 0.0, 0.0};
        const Vec3* pts[4] = {&p0, &p1, &p2, &p3};
        for (int i = 0; i < 4; ++i)
            n = detail::vec3_add(n, detail::vec3_cross(*pts[i], *pts[(i + 1) % 4]));
        area = 0.5 * detail::vec3_norm(n);
    }
    m[QM_VOLUME] = area;
    const double abs_area = std::fabs(area);
    const double lmin = std::min({l0, l1, l2, l3});
    const double lmax = std::max({l0, l1, l2, l3});
    const bool degenerate = abs_area < eps || lmin < eps;
    m[QM_INVERTED] = (is2d && area < 0.0) ? 1.0 : 0.0;
    m[QM_DEGENERATE] = degenerate ? 1.0 : 0.0;
    if (degenerate)
        return;

    m[QM_ASPECT_RATIO] = lmax / lmin;
    // Equiangle skewness over the four interior corner angles.
    const double t0 = quality_angle_deg(detail::vec3_sub(p1, p0), detail::vec3_sub(p3, p0), eps);
    const double t1 = quality_angle_deg(detail::vec3_sub(p2, p1), detail::vec3_sub(p0, p1), eps);
    const double t2 = quality_angle_deg(detail::vec3_sub(p3, p2), detail::vec3_sub(p1, p2), eps);
    const double t3 = quality_angle_deg(detail::vec3_sub(p0, p3), detail::vec3_sub(p2, p3), eps);
    m[QM_SKEWNESS] = quality_equiangle_skew(std::min({t0, t1, t2, t3}), std::max({t0, t1, t2, t3}));

    // Warpage: max out-of-plane angle over the two diagonal splits (degrees).
    auto split_warp = [&](const Vec3& a0, const Vec3& a1, const Vec3& a2, const Vec3& b0,
                          const Vec3& b1, const Vec3& b2) -> double {
        const Vec3 na = detail::vec3_cross(detail::vec3_sub(a1, a0), detail::vec3_sub(a2, a0));
        const Vec3 nb = detail::vec3_cross(detail::vec3_sub(b1, b0), detail::vec3_sub(b2, b0));
        const double lna = detail::vec3_norm(na);
        const double lnb = detail::vec3_norm(nb);
        if (lna < eps || lnb < eps)
            return QUALITY_NAN;
        double d = std::fabs(detail::vec3_dot(na, nb) / (lna * lnb));
        d = std::clamp(d, 0.0, 1.0);
        return std::acos(d) * QUALITY_RAD2DEG;
    };
    const double w0 = split_warp(p0, p1, p2, p2, p3, p0);  // diagonal 0-2
    const double w1 = split_warp(p1, p2, p3, p3, p0, p1);  // diagonal 1-3
    if (std::isfinite(w0) && std::isfinite(w1))
        m[QM_WARPAGE] = std::max(w0, w1);
}

// Solve the 3x3 system whose rows are a0,a1,a2 for x (Cramer). ok=false if singular.
Vec3 quality_solve3(const Vec3& a0, const Vec3& a1, const Vec3& a2, const Vec3& b, double eps,
                    bool& rOk) {
    auto det = [](const Vec3& r0, const Vec3& r1, const Vec3& r2) {
        return r0[0] * (r1[1] * r2[2] - r1[2] * r2[1]) - r0[1] * (r1[0] * r2[2] - r1[2] * r2[0]) +
               r0[2] * (r1[0] * r2[1] - r1[1] * r2[0]);
    };
    const double d = det(a0, a1, a2);
    if (std::fabs(d) < eps) {
        rOk = false;
        return {0.0, 0.0, 0.0};
    }
    rOk = true;
    Vec3 x;
    for (int k = 0; k < 3; ++k) {
        Vec3 c0 = a0, c1 = a1, c2 = a2;
        c0[k] = b[0];
        c1[k] = b[1];
        c2[k] = b[2];
        x[k] = det(c0, c1, c2) / d;
    }
    return x;
}

void quality_tetra(const std::vector<Vec3>& c, double eps, double* m) {
    const Vec3& p0 = c[0];
    const Vec3& p1 = c[1];
    const Vec3& p2 = c[2];
    const Vec3& p3 = c[3];
    const double vol = detail::triple_product(detail::vec3_sub(p1, p0), detail::vec3_sub(p2, p0),
                                              detail::vec3_sub(p3, p0)) /
                       6.0;
    m[QM_VOLUME] = vol;
    const double edges[6] = {quality_edge_len(p0, p1), quality_edge_len(p0, p2),
                             quality_edge_len(p0, p3), quality_edge_len(p1, p2),
                             quality_edge_len(p1, p3), quality_edge_len(p2, p3)};
    double emin = edges[0];
    for (double e : edges)
        emin = std::min(emin, e);

    // Scaled Jacobian: sqrt(2) * min corner det of unit outgoing edges.
    const double j0 = quality_det_unit3(detail::vec3_sub(p1, p0), detail::vec3_sub(p2, p0),
                                        detail::vec3_sub(p3, p0), eps);
    const double j1 = quality_det_unit3(detail::vec3_sub(p2, p1), detail::vec3_sub(p0, p1),
                                        detail::vec3_sub(p3, p1), eps);
    const double j2 = quality_det_unit3(detail::vec3_sub(p0, p2), detail::vec3_sub(p1, p2),
                                        detail::vec3_sub(p3, p2), eps);
    const double j3 = quality_det_unit3(detail::vec3_sub(p0, p3), detail::vec3_sub(p2, p3),
                                        detail::vec3_sub(p1, p3), eps);
    double sj = QUALITY_NAN;
    if (std::isfinite(j0) && std::isfinite(j1) && std::isfinite(j2) && std::isfinite(j3))
        sj = std::sqrt(2.0) * std::min({j0, j1, j2, j3});

    const double abs_vol = std::fabs(vol);
    const bool degenerate =
        abs_vol < eps || emin < eps || (std::isfinite(sj) && std::fabs(sj) < eps);
    m[QM_INVERTED] = vol < 0.0 ? 1.0 : 0.0;
    m[QM_DEGENERATE] = degenerate ? 1.0 : 0.0;
    if (degenerate)
        return;
    m[QM_SCALED_JACOBIAN] = sj;

    // Radius (aspect) ratio R / (3 r); 1 for a regular tetra.
    auto tri_area = [&](const Vec3& x, const Vec3& y, const Vec3& z) {
        return 0.5 * detail::vec3_norm(
                         detail::vec3_cross(detail::vec3_sub(y, x), detail::vec3_sub(z, x)));
    };
    const double face_sum =
        tri_area(p0, p1, p2) + tri_area(p0, p1, p3) + tri_area(p0, p2, p3) + tri_area(p1, p2, p3);
    if (face_sum > eps) {
        const double r = 3.0 * abs_vol / face_sum;
        // Circumcenter O: 2 (p_i - p0) . O = |p_i|^2 - |p0|^2, i = 1,2,3.
        const Vec3 a0 = detail::vec3_scale(detail::vec3_sub(p1, p0), 2.0);
        const Vec3 a1 = detail::vec3_scale(detail::vec3_sub(p2, p0), 2.0);
        const Vec3 a2 = detail::vec3_scale(detail::vec3_sub(p3, p0), 2.0);
        const Vec3 rhs = {detail::vec3_norm_sq(p1) - detail::vec3_norm_sq(p0),
                          detail::vec3_norm_sq(p2) - detail::vec3_norm_sq(p0),
                          detail::vec3_norm_sq(p3) - detail::vec3_norm_sq(p0)};
        bool ok = false;
        const Vec3 o = quality_solve3(a0, a1, a2, rhs, eps, ok);
        if (ok && r > eps) {
            const double R = detail::vec3_norm(detail::vec3_sub(o, p0));
            m[QM_ASPECT_RATIO] = R / (3.0 * r);
        }
    }

    // Min/max dihedral angle (degrees) from outward face normals.
    const int fverts[4][3] = {{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}};
    const int edge_faces[6][2] = {{0, 1}, {0, 2}, {1, 2}, {0, 3}, {1, 3}, {2, 3}};
    Vec3 centroid = detail::vec3_scale(
        detail::vec3_add(detail::vec3_add(p0, p1), detail::vec3_add(p2, p3)), 0.25);
    Vec3 fn[4];
    bool fn_ok = true;
    for (int f = 0; f < 4; ++f) {
        const Vec3& v0 = c[fverts[f][0]];
        const Vec3& v1 = c[fverts[f][1]];
        const Vec3& v2 = c[fverts[f][2]];
        Vec3 n = detail::vec3_cross(detail::vec3_sub(v1, v0), detail::vec3_sub(v2, v0));
        n = detail::vec3_normalize(n, eps);
        if (detail::vec3_norm(n) < 0.5) {
            fn_ok = false;
            break;
        }
        const Vec3 fc =
            detail::vec3_scale(detail::vec3_add(detail::vec3_add(v0, v1), v2), 1.0 / 3.0);
        if (detail::vec3_dot(n, detail::vec3_sub(fc, centroid)) < 0.0)
            n = detail::vec3_scale(n, -1.0);
        fn[f] = n;
    }
    if (fn_ok) {
        double dmin = 360.0, dmax = -360.0;
        for (int e = 0; e < 6; ++e) {
            const double dot =
                std::clamp(detail::vec3_dot(fn[edge_faces[e][0]], fn[edge_faces[e][1]]), -1.0, 1.0);
            const double dih = 180.0 - std::acos(dot) * QUALITY_RAD2DEG;
            dmin = std::min(dmin, dih);
            dmax = std::max(dmax, dih);
        }
        m[QM_MIN_DIHEDRAL] = dmin;
        m[QM_MAX_DIHEDRAL] = dmax;
    }
}

void quality_hex(const std::vector<Vec3>& c, const std::vector<detail::CellFaceDef>& rFaces,
                 double eps, double* m) {
    const double vol = quality_facefan_volume(c, rFaces);
    m[QM_VOLUME] = vol;

    // 12 edges.
    const int ep[12][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6},
                           {6, 7}, {7, 4}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};
    double emin = std::numeric_limits<double>::max();
    for (const auto& e : ep)
        emin = std::min(emin, quality_edge_len(c[e[0]], c[e[1]]));

    // Scaled Jacobian: min corner det of unit outgoing element edges.
    const int tj[8][3] = {{1, 3, 4}, {2, 0, 5}, {3, 1, 6}, {0, 2, 7},
                          {7, 5, 0}, {4, 6, 1}, {5, 7, 2}, {6, 4, 3}};
    double sj = std::numeric_limits<double>::max();
    bool sj_ok = true;
    for (int n = 0; n < 8; ++n) {
        const double d = quality_det_unit3(detail::vec3_sub(c[tj[n][0]], c[n]),
                                           detail::vec3_sub(c[tj[n][1]], c[n]),
                                           detail::vec3_sub(c[tj[n][2]], c[n]), eps);
        if (!std::isfinite(d)) {
            sj_ok = false;
            break;
        }
        sj = std::min(sj, d);
    }

    const double abs_vol = std::fabs(vol);
    const bool degenerate = abs_vol < eps || emin < eps || (sj_ok && std::fabs(sj) < eps) || !sj_ok;
    m[QM_INVERTED] = vol < 0.0 ? 1.0 : 0.0;
    m[QM_DEGENERATE] = degenerate ? 1.0 : 0.0;
    if (degenerate)
        return;
    m[QM_SCALED_JACOBIAN] = sj;

    // Equiangle skewness over the 24 face-corner angles.
    double amin = 360.0, amax = -360.0;
    for (const detail::CellFaceDef& f : rFaces) {
        for (std::uint8_t k = 0; k < f.mNumCorners; ++k) {
            const Vec3& pc = c[f.mNodes[k]];
            const Vec3& pn = c[f.mNodes[(k + 1) % f.mNumCorners]];
            const Vec3& pp = c[f.mNodes[(k + f.mNumCorners - 1) % f.mNumCorners]];
            const double ang =
                quality_angle_deg(detail::vec3_sub(pn, pc), detail::vec3_sub(pp, pc), eps);
            if (std::isfinite(ang)) {
                amin = std::min(amin, ang);
                amax = std::max(amax, ang);
            }
        }
    }
    if (amax > -360.0)
        m[QM_SKEWNESS] = quality_equiangle_skew(amin, amax);
}

void quality_wedge(const std::vector<Vec3>& c, const std::vector<detail::CellFaceDef>& rFaces,
                   double eps, double* m) {
    const double vol = quality_facefan_volume(c, rFaces);
    m[QM_VOLUME] = vol;

    const int ep[9][2] = {{0, 1}, {1, 2}, {2, 0}, {3, 4}, {4, 5}, {5, 3}, {0, 3}, {1, 4}, {2, 5}};
    double emin = std::numeric_limits<double>::max();
    for (const auto& e : ep)
        emin = std::min(emin, quality_edge_len(c[e[0]], c[e[1]]));

    // Scaled Jacobian at the 6 corners; (2/sqrt3) factor => 1 for a right prism.
    const int tj[6][3] = {{1, 2, 3}, {2, 0, 4}, {0, 1, 5}, {5, 4, 0}, {3, 5, 1}, {4, 3, 2}};
    double sj = std::numeric_limits<double>::max();
    bool sj_ok = true;
    for (int n = 0; n < 6; ++n) {
        const double d = quality_det_unit3(detail::vec3_sub(c[tj[n][0]], c[n]),
                                           detail::vec3_sub(c[tj[n][1]], c[n]),
                                           detail::vec3_sub(c[tj[n][2]], c[n]), eps);
        if (!std::isfinite(d)) {
            sj_ok = false;
            break;
        }
        sj = std::min(sj, d);
    }

    const double abs_vol = std::fabs(vol);
    const bool degenerate = abs_vol < eps || emin < eps || (sj_ok && std::fabs(sj) < eps) || !sj_ok;
    m[QM_INVERTED] = vol < 0.0 ? 1.0 : 0.0;
    m[QM_DEGENERATE] = degenerate ? 1.0 : 0.0;
    if (degenerate)
        return;
    m[QM_SCALED_JACOBIAN] = (2.0 / std::sqrt(3.0)) * sj;
}

void quality_pyramid(const std::vector<Vec3>& c, const std::vector<detail::CellFaceDef>& rFaces,
                     double eps, double* m) {
    const double vol = quality_facefan_volume(c, rFaces);
    m[QM_VOLUME] = vol;

    const int ep[8][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {0, 4}, {1, 4}, {2, 4}, {3, 4}};
    double emin = std::numeric_limits<double>::max();
    for (const auto& e : ep)
        emin = std::min(emin, quality_edge_len(c[e[0]], c[e[1]]));

    // Scaled Jacobian at the 4 base corners: two base edges + the apex edge.
    // sqrt(2) factor => 1 for a reference right square pyramid (J1).
    double sj = std::numeric_limits<double>::max();
    bool sj_ok = true;
    for (int i = 0; i < 4; ++i) {
        const int nxt = (i + 1) % 4;
        const int prv = (i + 3) % 4;
        const double d =
            quality_det_unit3(detail::vec3_sub(c[nxt], c[i]), detail::vec3_sub(c[prv], c[i]),
                              detail::vec3_sub(c[4], c[i]), eps);
        if (!std::isfinite(d)) {
            sj_ok = false;
            break;
        }
        sj = std::min(sj, d);
    }

    const double abs_vol = std::fabs(vol);
    const bool degenerate = abs_vol < eps || emin < eps || (sj_ok && std::fabs(sj) < eps) || !sj_ok;
    m[QM_INVERTED] = vol < 0.0 ? 1.0 : 0.0;
    m[QM_DEGENERATE] = degenerate ? 1.0 : 0.0;
    if (degenerate)
        return;
    m[QM_SCALED_JACOBIAN] = std::sqrt(2.0) * sj;
}

// Evaluate one cell's metrics into m (assumed pre-filled with NaN).
void quality_eval_cell(QualityFamily family, CellType ct, const std::vector<Vec3>& c, bool is2d,
                       double eps, double* m) {
    switch (family) {
        case QualityFamily::Tri:
            quality_tri(c, is2d, eps, m);
            break;
        case QualityFamily::Quad:
            quality_quad(c, is2d, eps, m);
            break;
        case QualityFamily::Tetra:
            quality_tetra(c, eps, m);
            break;
        case QualityFamily::Hex:
            quality_hex(c, detail::cell_faces(ct), eps, m);
            break;
        case QualityFamily::Wedge:
            quality_wedge(c, detail::cell_faces(ct), eps, m);
            break;
        case QualityFamily::Pyr:
            quality_pyramid(c, detail::cell_faces(ct), eps, m);
            break;
        default:
            break;
    }
}

using CellMetrics = std::array<double, NUM_METRICS>;

// A deep copy of an NDArray that always owns its buffer (safe even if the
// source is a view over foreign memory).
NDArray quality_owned_copy(const NDArray& rArr) {
    NDArray c = rArr;
    c.MakeOwned();
    return c;
}

// Backend-agnostic mesh clone through the uniform API (the KRATOS backend's
// Mesh is not copy-constructible, so we cannot rely on `Mesh out = in;`).
Mesh quality_clone_mesh(const Mesh& rMesh) {
    Mesh out;
    out.AssignPoints(quality_owned_copy(rMesh.Points()));
    for (const auto cb : rMesh.CellRange()) {
        if (cb.IsPolyhedron()) {
            std::vector<std::vector<std::vector<std::int64_t>>> cells(cb.NumCells());
            for (std::size_t c = 0; c < cb.NumCells(); ++c) {
                cells[c].resize(cb.NumFaces(c));
                for (std::size_t f = 0; f < cb.NumFaces(c); ++f) {
                    auto face = cb.Face(c, f);
                    cells[c][f].assign(face.first, face.first + face.second);
                }
            }
            out.AddPolyhedronBlock(std::string(cb.Type()), std::move(cells));
        } else if (cb.IsRagged()) {
            std::vector<std::vector<std::int64_t>> rows(cb.NumCells());
            for (std::size_t c = 0; c < cb.NumCells(); ++c)
                rows[c].assign(cb.Row(c), cb.Row(c) + cb.RowSize(c));
            out.AddPolygonBlock(std::string(cb.Type()), std::move(rows));
        } else {
            out.AddCellBlock(std::string(cb.Type()), quality_owned_copy(cb.Conn()));
        }
    }
    for (const std::string& name : rMesh.PointDataNames())
        out.AddPointData(name, quality_owned_copy(rMesh.PointData(name)));
    for (const std::string& name : rMesh.CellDataNames()) {
        std::vector<NDArray> blocks;
        for (std::size_t b = 0; b < rMesh.CellDataNumBlocks(name); ++b)
            blocks.push_back(quality_owned_copy(rMesh.CellData(name, b)));
        out.AddCellData(name, std::move(blocks));
    }
    for (const std::string& name : rMesh.FieldDataNames())
        out.AddFieldData(name, quality_owned_copy(rMesh.FieldData(name)));
    // Named regions pass through verbatim: attaching quality metrics adds
    // cell_data and renumbers nothing.
    for (std::size_t i = 0; i < rMesh.NumRegions(); ++i)
        out.AddRegion(rMesh.Region(i));
    // Property sets ride along too: these operations preserve the mesh's shape,
    // so dropping a deck's material data here would be new lossiness. They are
    // keyed by id, not by entity index, so there is nothing to remap.
    for (std::size_t i = 0; i < rMesh.NumPropertySets(); ++i)
        out.AddPropertySet(rMesh.GetPropertySet(i));
    return out;
}

}  // namespace

QualityReport compute_quality(const Mesh& rMesh) {
    const double eps = 1e-14;
    const NDArray& points = rMesh.Points();
    const std::size_t pdim = rMesh.PointDim();
    const bool is2d = pdim == 2;

    QualityReport rep;

    // --- per-block per-cell metric values (compute-bound → parallel) ---
    std::vector<std::vector<CellMetrics>> block_vals;
    block_vals.reserve(rMesh.NumCellBlocks());
    for (const auto cb : rMesh.CellRange()) {
        const CellType ct = cell_type_from_name(cb.Type());
        const std::size_t nc = cb.NumCells();
        rep.mNumCells += static_cast<std::int64_t>(nc);
        std::vector<CellMetrics> vals(nc);
        const QualityFamily family =
            (cb.IsRagged() || cb.IsPolyhedron()) ? QualityFamily::None : quality_family(ct);
        if (family == QualityFamily::None) {
            for (CellMetrics& v : vals)
                v.fill(QUALITY_NAN);
        } else {
            const NDArray& conn = cb.Conn();
            const std::size_t npc = cb.NodesPerCell();
            const int corner_count = detail::cell_corner_count(ct);
            parallel_for(nc, [&](std::size_t i) {
                CellMetrics& v = vals[i];
                v.fill(QUALITY_NAN);
                std::vector<Vec3> coords;
                detail::read_corner_coords(points, pdim, conn, i * npc,
                                           static_cast<std::size_t>(corner_count), coords);
                quality_eval_cell(family, ct, coords, is2d, eps, v.data());
            });
        }
        block_vals.push_back(std::move(vals));
    }

    // --- reduction: summaries + inverted/degenerate counts (serial) ---
    struct Acc {
        double mMin = std::numeric_limits<double>::infinity();
        double mMax = -std::numeric_limits<double>::infinity();
        double mSum = 0.0;
        std::int64_t mCount = 0;
    };
    std::array<Acc, NUM_METRICS> acc;
    for (const std::vector<CellMetrics>& vals : block_vals) {
        for (const CellMetrics& v : vals) {
            for (int mi = 0; mi < NUM_METRICS; ++mi) {
                const double x = v[mi];
                if (!std::isfinite(x))
                    continue;
                Acc& a = acc[mi];
                a.mMin = std::min(a.mMin, x);
                a.mMax = std::max(a.mMax, x);
                a.mSum += x;
                ++a.mCount;
            }
            if (v[QM_INVERTED] == 1.0)
                ++rep.mNumInverted;
            if (v[QM_DEGENERATE] == 1.0)
                ++rep.mNumDegenerate;
        }
    }

    // --- histograms (second pass now that min/max are known) ---
    std::array<QualityMetricSummary, NUM_METRICS> summ;
    for (int mi = 0; mi < NUM_METRICS; ++mi) {
        const Acc& a = acc[mi];
        QualityMetricSummary& s = summ[mi];
        s.mCount = a.mCount;
        if (a.mCount == 0)
            continue;
        s.mMin = a.mMin;
        s.mMax = a.mMax;
        s.mMean = a.mSum / static_cast<double>(a.mCount);
    }
    // Fixed-size chunks over the block-major global cell index, one partial
    // histogram per chunk, merged serially in chunk order (stats.cpp's chunked
    // reduction; bins are integer counts, so the merge is bit-identical to the
    // serial loop for any thread count).
    const int K = QualityMetricSummary::K;
    std::vector<std::size_t> starts(block_vals.size() + 1, 0);
    for (std::size_t bi = 0; bi < block_vals.size(); ++bi)
        starts[bi + 1] = starts[bi] + block_vals[bi].size();
    const std::size_t total_cells = starts.back();
    constexpr std::size_t hist_chunk = 4096;
    const std::size_t nchunks = (total_cells + hist_chunk - 1) / hist_chunk;
    using Hist = std::array<std::array<std::int64_t, QualityMetricSummary::K>, NUM_METRICS>;
    std::vector<Hist> partial(nchunks);
    parallel_for(
        nchunks,
        [&](std::size_t ci) {
            Hist& h = partial[ci];
            const std::size_t lo = ci * hist_chunk;
            const std::size_t hi = std::min(total_cells, lo + hist_chunk);
            std::size_t bi = static_cast<std::size_t>(
                std::upper_bound(starts.begin(), starts.end(), lo) - starts.begin() - 1);
            for (std::size_t g = lo; g < hi; ++g) {
                while (g >= starts[bi + 1])
                    ++bi;
                const CellMetrics& v = block_vals[bi][g - starts[bi]];
                for (int mi = 0; mi < NUM_METRICS; ++mi) {
                    const double x = v[mi];
                    if (!std::isfinite(x))
                        continue;
                    const QualityMetricSummary& s = summ[mi];
                    const double span = s.mMax - s.mMin;
                    int bin = 0;
                    if (span > 0.0)
                        bin = static_cast<int>((x - s.mMin) / span * K);
                    bin = std::clamp(bin, 0, K - 1);
                    ++h[mi][bin];
                }
            }
        },
        /*grain=*/1);
    for (const Hist& h : partial)
        for (int mi = 0; mi < NUM_METRICS; ++mi)
            for (int k = 0; k < K; ++k)
                summ[mi].mHistogram[k] += h[mi][k];

    // --- assemble the report ---
    for (int mi = 0; mi < NUM_METRICS; ++mi) {
        rep.mMetrics.emplace_back(QUALITY_METRIC_NAMES[mi], summ[mi]);
        std::vector<NDArray> arrays;
        arrays.reserve(block_vals.size());
        for (const std::vector<CellMetrics>& vals : block_vals) {
            NDArray a = NDArray::Uninit(DType::Float64, {vals.size(), 1});
            double* d = a.As<double>();
            parallel_for_bw(vals.size(), [&](std::size_t c) { d[c] = vals[c][mi]; });
            arrays.push_back(std::move(a));
        }
        rep.mCellArrays.emplace_back(QUALITY_METRIC_NAMES[mi], std::move(arrays));
    }

    return rep;
}

Mesh attach_quality(const Mesh& rMesh) {
    QualityReport rep = compute_quality(rMesh);
    Mesh out = quality_clone_mesh(rMesh);
    if (rMesh.NumCellBlocks() == 0)
        return out;
    for (std::pair<std::string, std::vector<NDArray>>& metric : rep.mCellArrays)
        out.AddCellData(metric.first, std::move(metric.second));
    return out;
}

}  // namespace meshioplusplus

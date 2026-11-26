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
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/polyhedron.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/parallel.hpp"

// Project includes (private, not installed)
#include "../detail/fem_kernels.hpp"
#include "../detail/math_utils.hpp"

namespace meshioplusplus {

namespace {

using detail::Vec3;

constexpr int NUM_METRICS = 22;
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
    // Cell measures: size is the unsigned area/volume, the rest are lengths.
    QM_SIZE,
    QM_DIAGONAL,
    QM_MIN_EDGE,
    QM_MAX_EDGE,
    // The Verdict metrics beyond the ones above (SAND2007-1751): each is set
    // only for the families it is defined for.
    QM_EDGE_RATIO,
    QM_CONDITION,
    QM_SHAPE,
    QM_ASPECT_GAMMA,
    QM_STRETCH,
    QM_TAPER,
    QM_DIAGONAL_RATIO,
};

const char* const QUALITY_METRIC_NAMES[NUM_METRICS] = {
    "quality:volume",         "quality:scaled_jacobian", "quality:aspect_ratio",
    "quality:skewness",       "quality:min_angle",       "quality:max_angle",
    "quality:warpage",        "quality:min_dihedral",    "quality:max_dihedral",
    "quality:inverted",       "quality:degenerate",      "quality:size",
    "quality:diagonal",       "quality:min_edge",        "quality:max_edge",
    "quality:edge_ratio",     "quality:condition",       "quality:shape",
    "quality:aspect_gamma",   "quality:stretch",         "quality:taper",
    "quality:diagonal_ratio",
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
    // The quad's normal: +z in 2-D, the sum of its edges' cross products
    // (twice its vector area) in 3-D.
    Vec3 normal = {0.0, 0.0, 1.0};
    if (is2d) {
        area = 0.5 * ((p0[0] * p1[1] - p1[0] * p0[1]) + (p1[0] * p2[1] - p2[0] * p1[1]) +
                      (p2[0] * p3[1] - p3[0] * p2[1]) + (p3[0] * p0[1] - p0[0] * p3[1]));
    } else {
        Vec3 n = {0.0, 0.0, 0.0};
        const Vec3* pts[4] = {&p0, &p1, &p2, &p3};
        for (int i = 0; i < 4; ++i)
            n = detail::vec3_add(n, detail::vec3_cross(*pts[i], *pts[(i + 1) % 4]));
        area = 0.5 * detail::vec3_norm(n);
        if (area > 0.0)
            normal = {n[0] / (2.0 * area), n[1] / (2.0 * area), n[2] / (2.0 * area)};
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
    // Scaled Jacobian: the least over the corners of (e1 x e2) . n / (|e1|
    // |e2|), e1 and e2 the corner's two edges, n the quad's normal (Verdict's
    // definition); 1 for a rectangle, negative at a reflex or inverted corner.
    {
        const Vec3* pts[4] = {&p0, &p1, &p2, &p3};
        double sj = std::numeric_limits<double>::infinity();
        for (int i = 0; i < 4; ++i) {
            const Vec3 e1 = detail::vec3_sub(*pts[(i + 1) % 4], *pts[i]);
            const Vec3 e2 = detail::vec3_sub(*pts[(i + 3) % 4], *pts[i]);
            const double l = detail::vec3_norm(e1) * detail::vec3_norm(e2);
            sj = std::min(sj, detail::vec3_dot(detail::vec3_cross(e1, e2), normal) / l);
        }
        m[QM_SCALED_JACOBIAN] = sj;
    }
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
    return detail::MathUtils::Solve3(a0, a1, a2, b, eps, rOk);
}

void quality_tetra(const std::vector<Vec3>& c, double eps, double* m) {
    const Vec3& p0 = c[0];
    const Vec3& p1 = c[1];
    const Vec3& p2 = c[2];
    const Vec3& p3 = c[3];
    const double vol = detail::tet_volume(p0, p1, p2, p3);
    m[QM_VOLUME] = vol;
    const double edges[6] = {quality_edge_len(p0, p1), quality_edge_len(p0, p2),
                             quality_edge_len(p0, p3), quality_edge_len(p1, p2),
                             quality_edge_len(p1, p3), quality_edge_len(p2, p3)};
    double emin = edges[0];
    for (double e : edges)
        emin = std::min(emin, e);

    // Scaled Jacobian: sqrt(2) * min corner det of unit outgoing edges.
    const double j0 = detail::det_unit3(detail::vec3_sub(p1, p0), detail::vec3_sub(p2, p0),
                                        detail::vec3_sub(p3, p0), eps);
    const double j1 = detail::det_unit3(detail::vec3_sub(p2, p1), detail::vec3_sub(p0, p1),
                                        detail::vec3_sub(p3, p1), eps);
    const double j2 = detail::det_unit3(detail::vec3_sub(p0, p2), detail::vec3_sub(p1, p2),
                                        detail::vec3_sub(p3, p2), eps);
    const double j3 = detail::det_unit3(detail::vec3_sub(p0, p3), detail::vec3_sub(p2, p3),
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
    double sj = std::numeric_limits<double>::max();
    bool sj_ok = true;
    for (int n = 0; n < 8; ++n) {
        const double d =
            detail::det_unit3(detail::vec3_sub(c[detail::kHexCornerEdges[n][0]], c[n]),
                              detail::vec3_sub(c[detail::kHexCornerEdges[n][1]], c[n]),
                              detail::vec3_sub(c[detail::kHexCornerEdges[n][2]], c[n]), eps);
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
    double sj = std::numeric_limits<double>::max();
    bool sj_ok = true;
    for (int n = 0; n < 6; ++n) {
        const double d =
            detail::det_unit3(detail::vec3_sub(c[detail::kWedgeCornerEdges[n][0]], c[n]),
                              detail::vec3_sub(c[detail::kWedgeCornerEdges[n][1]], c[n]),
                              detail::vec3_sub(c[detail::kWedgeCornerEdges[n][2]], c[n]), eps);
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
            detail::det_unit3(detail::vec3_sub(c[nxt], c[i]), detail::vec3_sub(c[prv], c[i]),
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

// The edges of each family, as corner-index pairs in VTK order.
constexpr int kQualityTriEdges[3][2] = {{0, 1}, {1, 2}, {2, 0}};
constexpr int kQualityQuadEdges[4][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
constexpr int kQualityTetraEdges[6][2] = {{0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3}};
constexpr int kQualityHexEdges[12][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6},
                                         {6, 7}, {7, 4}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};
constexpr int kQualityWedgeEdges[9][2] = {{0, 1}, {1, 2}, {2, 0}, {3, 4}, {4, 5},
                                          {5, 3}, {0, 3}, {1, 4}, {2, 5}};
constexpr int kQualityPyramidEdges[8][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0},
                                            {0, 4}, {1, 4}, {2, 4}, {3, 4}};

// Size, diagonal and the shortest and longest edge, written after the family's
// own metrics so a cell that the family call left early (degenerate) still
// reports its lengths. size is the unsigned area/volume the family computed.
/// The Frobenius condition number `|F| |F^-1| / 3` of the 3x3 matrix whose
/// columns are @p f1, @p f2, @p f3 (the inverse's rows are the columns' cross
/// products over the determinant). NaN when the determinant is zero.
double quality_cond3(const Vec3& f1, const Vec3& f2, const Vec3& f3) {
    const double det = detail::det3(f1, f2, f3);
    if (det == 0.0 || !std::isfinite(det))
        return QUALITY_NAN;
    const double fn2 =
        detail::vec3_norm_sq(f1) + detail::vec3_norm_sq(f2) + detail::vec3_norm_sq(f3);
    const double in2 = detail::vec3_norm_sq(detail::vec3_cross(f2, f3)) +
                       detail::vec3_norm_sq(detail::vec3_cross(f3, f1)) +
                       detail::vec3_norm_sq(detail::vec3_cross(f1, f2));
    return std::sqrt(fn2) * std::sqrt(in2) / (std::fabs(det) * 3.0);
}

// The Verdict metrics after the cell measures (m holds the edge lengths).
void quality_verdict_metrics(QualityFamily family, const std::vector<Vec3>& c, double* m) {
    const double emin = m[QM_MIN_EDGE];
    const double emax = m[QM_MAX_EDGE];
    if (emin > 0.0 && std::isfinite(emax))
        m[QM_EDGE_RATIO] = emax / emin;
    const double s3 = 1.0 / std::sqrt(3.0);
    const double s6 = 1.0 / std::sqrt(6.0);
    switch (family) {
        case QualityFamily::Tri: {
            // F = [e1 e2] W^-1 against the equilateral triangle; condition is
            // tr(M) / (2 sqrt(det M)) for M = F^T F, and shape its reciprocal.
            const Vec3 e1 = detail::vec3_sub(c[1], c[0]);
            const Vec3 e2 = detail::vec3_sub(c[2], c[0]);
            const Vec3 f2 =
                detail::vec3_add(detail::vec3_scale(e1, -s3), detail::vec3_scale(e2, 2.0 * s3));
            const double m11 = detail::vec3_norm_sq(e1);
            const double m22 = detail::vec3_norm_sq(f2);
            const double m12 = detail::vec3_dot(e1, f2);
            const double det = m11 * m22 - m12 * m12;
            if (det > 0.0 && std::isfinite(det)) {
                m[QM_CONDITION] = (m11 + m22) / (2.0 * std::sqrt(det));
                m[QM_SHAPE] = 2.0 * std::sqrt(det) / (m11 + m22);
            }
            break;
        }
        case QualityFamily::Quad: {
            double cond = 0.0;
            bool ok = true;
            for (int k = 0; k < 4 && ok; ++k) {
                const Vec3 a = detail::vec3_sub(c[(k + 1) % 4], c[k]);
                const Vec3 b = detail::vec3_sub(c[(k + 3) % 4], c[k]);
                const double area = detail::vec3_norm(detail::vec3_cross(a, b));
                if (!(area > 0.0)) {
                    ok = false;
                    break;
                }
                cond = std::max(cond,
                                (detail::vec3_norm_sq(a) + detail::vec3_norm_sq(b)) / (2.0 * area));
            }
            if (ok)
                m[QM_CONDITION] = cond;
            const double d = std::max(quality_edge_len(c[0], c[2]), quality_edge_len(c[1], c[3]));
            if (d > 0.0)
                m[QM_STRETCH] = std::sqrt(2.0) * emin / d;
            const Vec3 x1 =
                detail::vec3_add(detail::vec3_sub(c[1], c[0]), detail::vec3_sub(c[2], c[3]));
            const Vec3 x2 =
                detail::vec3_add(detail::vec3_sub(c[3], c[0]), detail::vec3_sub(c[2], c[1]));
            const Vec3 x12 =
                detail::vec3_add(detail::vec3_sub(c[0], c[1]), detail::vec3_sub(c[2], c[3]));
            const double den = std::min(detail::vec3_norm(x1), detail::vec3_norm(x2));
            if (den > 0.0)
                m[QM_TAPER] = detail::vec3_norm(x12) / den;
            break;
        }
        case QualityFamily::Tetra: {
            // F = [e1 e2 e3] W^-1 against the regular tetrahedron.
            const Vec3 e1 = detail::vec3_sub(c[1], c[0]);
            const Vec3 e2 = detail::vec3_sub(c[2], c[0]);
            const Vec3 e3 = detail::vec3_sub(c[3], c[0]);
            const Vec3 f2 =
                detail::vec3_add(detail::vec3_scale(e1, -s3), detail::vec3_scale(e2, 2.0 * s3));
            const Vec3 f3 = detail::vec3_add(
                detail::vec3_add(detail::vec3_scale(e1, -s6), detail::vec3_scale(e2, -s6)),
                detail::vec3_scale(e3, 3.0 * s6));
            m[QM_CONDITION] = quality_cond3(e1, f2, f3);
            const double det = detail::det3(e1, f2, f3);
            if (det == 0.0 || !std::isfinite(det)) {
                // degenerate: shape stays NaN
            } else if (det < 0.0) {
                m[QM_SHAPE] = 0.0;
            } else {
                const double fn2 =
                    detail::vec3_norm_sq(e1) + detail::vec3_norm_sq(f2) + detail::vec3_norm_sq(f3);
                m[QM_SHAPE] = 3.0 * std::cbrt(det * det) / fn2;
            }
            // aspect gamma: (root mean square edge)^3 / (6 sqrt(2) V).
            const double vol = std::fabs(detail::det3(e1, e2, e3)) / 6.0;
            double sum2 = 0.0;
            for (int a = 0; a < 4; ++a)
                for (int b = a + 1; b < 4; ++b)
                    sum2 += detail::vec3_norm_sq(detail::vec3_sub(c[b], c[a]));
            const double srms = std::sqrt(sum2 / 6.0);
            if (vol > 0.0)
                m[QM_ASPECT_GAMMA] = srms * srms * srms / (6.0 * std::sqrt(2.0) * vol);
            break;
        }
        case QualityFamily::Hex: {
            static const int kCorner[8][4] = {{0, 1, 3, 4}, {1, 2, 0, 5}, {2, 3, 1, 6},
                                              {3, 0, 2, 7}, {4, 7, 5, 0}, {5, 4, 6, 1},
                                              {6, 5, 7, 2}, {7, 6, 4, 3}};
            double cond = 0.0;
            bool ok = true;
            for (int k = 0; k < 8 && ok; ++k) {
                const int* q = kCorner[k];
                const double cc = quality_cond3(detail::vec3_sub(c[q[1]], c[q[0]]),
                                                detail::vec3_sub(c[q[2]], c[q[0]]),
                                                detail::vec3_sub(c[q[3]], c[q[0]]));
                if (!std::isfinite(cc))
                    ok = false;
                else
                    cond = std::max(cond, cc);
            }
            if (ok)
                m[QM_CONDITION] = cond;
            const double d[4] = {quality_edge_len(c[0], c[6]), quality_edge_len(c[1], c[7]),
                                 quality_edge_len(c[2], c[4]), quality_edge_len(c[3], c[5])};
            const double dmax = std::max(std::max(d[0], d[1]), std::max(d[2], d[3]));
            const double dmin = std::min(std::min(d[0], d[1]), std::min(d[2], d[3]));
            if (dmax > 0.0) {
                m[QM_DIAGONAL_RATIO] = dmin / dmax;
                m[QM_STRETCH] = std::sqrt(3.0) * emin / dmax;
            }
            break;
        }
        default:
            break;
    }
}

void quality_cell_measures(QualityFamily family, const std::vector<Vec3>& c, double* m) {
    const int(*edges)[2] = nullptr;
    int ne = 0;
    switch (family) {
        case QualityFamily::Tri:
            edges = kQualityTriEdges;
            ne = 3;
            break;
        case QualityFamily::Quad:
            edges = kQualityQuadEdges;
            ne = 4;
            break;
        case QualityFamily::Tetra:
            edges = kQualityTetraEdges;
            ne = 6;
            break;
        case QualityFamily::Hex:
            edges = kQualityHexEdges;
            ne = 12;
            break;
        case QualityFamily::Wedge:
            edges = kQualityWedgeEdges;
            ne = 9;
            break;
        case QualityFamily::Pyr:
            edges = kQualityPyramidEdges;
            ne = 8;
            break;
        default:
            return;
    }
    if (std::isfinite(m[QM_VOLUME]))
        m[QM_SIZE] = std::fabs(m[QM_VOLUME]);
    double emin = std::numeric_limits<double>::infinity();
    double emax = 0.0;
    for (int k = 0; k < ne; ++k) {
        const double l = quality_edge_len(c[edges[k][0]], c[edges[k][1]]);
        emin = std::min(emin, l);
        emax = std::max(emax, l);
    }
    m[QM_MIN_EDGE] = emin;
    m[QM_MAX_EDGE] = emax;
    double diag = 0.0;
    for (std::size_t a = 0; a < c.size(); ++a)
        for (std::size_t b = a + 1; b < c.size(); ++b)
            diag = std::max(diag, quality_edge_len(c[a], c[b]));
    m[QM_DIAGONAL] = diag;
    quality_verdict_metrics(family, c, m);
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
    quality_cell_measures(family, c, m);
}

using CellMetrics = std::array<double, NUM_METRICS>;

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
        if (cb.IsPolyhedron()) {
            // A REDUCED set, deliberately: volume, inverted and degenerate are
            // the three that are well defined for a cell bounded by arbitrary
            // polygons. Aspect ratio, skewness, warpage and the angle metrics
            // are all defined against a reference element a polyhedron does not
            // have; inventing plausible-looking numbers for them would be worse
            // than NaN, which is `compute_quality`'s standing convention for a
            // metric that does not apply.
            //
            // "Inverted" here means UNORIENTABLE, not negative-volume: a
            // polyhedron's stored winding is not a contract (doc/polyhedra.md),
            // so there is no convention for a cell to violate -- what a cell can
            // fail to be is a closed orientable surface.
            parallel_for(nc, [&](std::size_t i) {
                CellMetrics& v = vals[i];
                v.fill(QUALITY_NAN);
                // Per-thread scratch: cell_rings clears both on entry, so each
                // cell reuses the capacity the thread's largest cell grew.
                static thread_local detail::CellRings rings;
                static thread_local std::vector<Vec3> coords;
                if (!detail::cell_rings(cb, i, points, pdim, rings, coords))
                    return;
                const bool orientable = detail::orient_rings(rings, coords.data()) !=
                                        detail::RingOrientation::Unorientable;
                v[QM_INVERTED] = orientable ? 0.0 : 1.0;
                if (!orientable) {
                    v[QM_DEGENERATE] = 1.0;
                    return;
                }
                const detail::PolyMeasure pm = detail::poly_measure(rings, coords.data());
                v[QM_VOLUME] = std::abs(pm.mVolume);
                v[QM_SIZE] = v[QM_VOLUME];
                // Relative to the cell's own size, never an absolute epsilon --
                // the same rule grad_green_gauss_3d uses.
                const double scale = pm.mSurfaceArea * std::sqrt(pm.mSurfaceArea);
                v[QM_DEGENERATE] = (std::abs(pm.mVolume) > eps * scale) ? 0.0 : 1.0;
            });
            block_vals.push_back(std::move(vals));
            continue;
        }
        const QualityFamily family = cb.IsRagged() ? QualityFamily::None : quality_family(ct);
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
                // Per-thread scratch (at most 8 corners): read_corner_coords
                // clears it on entry.
                static thread_local std::vector<Vec3> coords;
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
    using Hist = std::array<std::array<std::int64_t, QualityMetricSummary::K>, NUM_METRICS>;
    const Hist hist = parallel_reduce(
        total_cells, hist_chunk, Hist{},
        [&](std::size_t lo, std::size_t hi) {
            Hist h{};
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
            return h;
        },
        [&](Hist acc, const Hist& h) {
            for (int mi = 0; mi < NUM_METRICS; ++mi)
                for (int k = 0; k < K; ++k)
                    acc[mi][k] += h[mi][k];
            return acc;
        });
    for (int mi = 0; mi < NUM_METRICS; ++mi)
        for (int k = 0; k < K; ++k)
            summ[mi].mHistogram[k] += hist[mi][k];

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
    Mesh out = detail::clone_mesh(rMesh);
    if (rMesh.NumCellBlocks() == 0)
        return out;
    for (std::pair<std::string, std::vector<NDArray>>& metric : rep.mCellArrays)
        out.AddCellData(metric.first, std::move(metric.second));
    return out;
}

}  // namespace meshioplusplus

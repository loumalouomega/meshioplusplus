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
#include <stdexcept>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/operations/auto_align.hpp"
#include "../detail/lsq_fit.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/surface_distance.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/operations/primitives.hpp"
#include "meshioplusplus/operations/transform.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {

namespace {

using AaVec = detail::Vec3;

constexpr const char* kAaPrefix = "meshio++: auto_align: ";
constexpr double kAaPi = 3.14159265358979323846;

/// The unit normals and areas of the voting facets.
struct AaFacets {
    std::vector<AaVec> mNormal;
    std::vector<double> mArea;
    std::size_t mNumReal = 0;  ///< facets of the input; a pseudo-facet may follow
    double mRealArea = 0.0;
};

/// The area-weighted unit-normal mean of the facets within `CosTol` of
/// `rDir`, or `rDir` itself when none is.
AaVec aa_refine(const AaFacets& rF, const AaVec& rDir, double CosTol) {
    AaVec sum{0.0, 0.0, 0.0};
    for (std::size_t t = 0; t < rF.mNormal.size(); ++t)
        if (detail::vec3_dot(rF.mNormal[t], rDir) >= CosTol)
            sum = detail::vec3_add(sum, detail::vec3_scale(rF.mNormal[t], rF.mArea[t]));
    const double len = detail::vec3_norm(sum);
    return len > 0.0 ? detail::vec3_scale(sum, 1.0 / len) : rDir;
}

/// `rV` with its component along the unit `rAxis` removed, normalised; zero
/// when nothing is left.
AaVec aa_orthogonal(const AaVec& rV, const AaVec& rAxis) {
    const AaVec o = detail::vec3_sub(rV, detail::vec3_scale(rAxis, detail::vec3_dot(rV, rAxis)));
    const double len = detail::vec3_norm(o);
    return len > 1e-12 ? detail::vec3_scale(o, 1.0 / len) : AaVec{0.0, 0.0, 0.0};
}

}  // namespace

AutoAlignResult auto_align(const Mesh& rMesh, const AutoAlignOptions& rOptions) {
    if (!(rOptions.mAngleToleranceDeg > 0.0 && rOptions.mAngleToleranceDeg < 45.0))
        throw std::invalid_argument(std::string(kAaPrefix) +
                                    "angle_tolerance must lie in (0, 45) degrees");
    if (rOptions.mSubdivisions < 0 || rOptions.mSubdivisions > 5)
        throw std::invalid_argument(std::string(kAaPrefix) + "subdivisions must lie in [0, 5]");
    const std::size_t n = rMesh.NumPoints();
    const std::size_t dim = rMesh.PointDim();
    for (std::size_t i = 0; i < n * dim; ++i)
        if (!std::isfinite(detail::read_double(rMesh.Points(), i)))
            throw std::invalid_argument(std::string(kAaPrefix) + "point " +
                                        std::to_string(i / std::max<std::size_t>(dim, 1)) +
                                        " has a non-finite coordinate");
    const detail::TriangleSoup soup = detail::build_triangle_soup(rMesh, rOptions.mRegion);
    AaFacets f;
    AaVec vector_area{0.0, 0.0, 0.0};
    for (std::size_t t = 0; t < soup.NumTriangles(); ++t) {
        const AaVec& a = soup.mCorners[3 * t];
        const AaVec x = detail::vec3_cross(detail::vec3_sub(soup.mCorners[3 * t + 1], a),
                                           detail::vec3_sub(soup.mCorners[3 * t + 2], a));
        const double len = detail::vec3_norm(x);
        if (!(len > 0.0))
            continue;
        f.mNormal.push_back(detail::vec3_scale(x, 1.0 / len));
        f.mArea.push_back(0.5 * len);
        f.mRealArea += 0.5 * len;
        vector_area = detail::vec3_add(vector_area, detail::vec3_scale(x, 0.5));
    }
    f.mNumReal = f.mNormal.size();
    if (f.mNumReal == 0)
        throw std::invalid_argument(std::string(kAaPrefix) +
                                    "the input has no surface cell of positive area");
    AutoAlignResult result;
    const double open = detail::vec3_norm(vector_area);
    if (open > 1e-3 * f.mRealArea) {
        f.mNormal.push_back(detail::vec3_scale(vector_area, -1.0 / open));
        f.mArea.push_back(open);
        result.mOpeningArea = open;
    }

    // The Gauss-sphere histogram.
    const Mesh ico = icosa_sphere(1.0, rOptions.mSubdivisions);
    const std::size_t nb = ico.NumPoints();
    std::vector<AaVec> bins(nb);
    for (std::size_t b = 0; b < nb; ++b) {
        for (std::size_t k = 0; k < 3; ++k)
            bins[b][k] = detail::read_double(ico.Points(), b * ico.PointDim() + k);
        bins[b] = detail::vec3_scale(bins[b], 1.0 / detail::vec3_norm(bins[b]));
    }
    std::vector<std::size_t> bin_of(f.mNormal.size(), 0);
    parallel_for(f.mNormal.size(), [&](std::size_t t) {
        std::size_t best = 0;
        double best_dot = -2.0;
        for (std::size_t b = 0; b < nb; ++b) {
            const double d = detail::vec3_dot(f.mNormal[t], bins[b]);
            if (d > best_dot) {
                best_dot = d;
                best = b;
            }
        }
        bin_of[t] = best;
    });
    std::vector<double> mass(nb, 0.0);
    for (std::size_t t = 0; t < f.mNormal.size(); ++t)
        mass[bin_of[t]] += f.mArea[t];
    const double tol = rOptions.mAngleToleranceDeg * kAaPi / 180.0;
    const double cos_tol = std::cos(tol), sin_tol = std::sin(tol);
    std::vector<double> score(nb, 0.0);
    for (std::size_t b = 0; b < nb; ++b)
        for (std::size_t c = 0; c < nb; ++c)
            if (detail::vec3_dot(bins[b], bins[c]) >= cos_tol)
                score[b] += mass[c];
    const auto best_bin = [&](const AaVec* pPerp) {
        std::size_t best = nb;
        for (std::size_t b = 0; b < nb; ++b) {
            if (pPerp != nullptr && std::abs(detail::vec3_dot(bins[b], *pPerp)) > sin_tol)
                continue;
            if (!(score[b] > 0.0))
                continue;
            if (best == nb || score[b] > score[best])
                best = b;
        }
        return best;
    };
    const auto second_dir = [&](const AaVec& rFirst) {
        const std::size_t b = best_bin(&rFirst);
        AaVec d{0.0, 0.0, 0.0};
        if (b < nb)
            d = aa_orthogonal(aa_refine(f, bins[b], cos_tol), rFirst);
        if (detail::vec3_norm(d) == 0.0) {
            result.mWarnings.push_back(
                "no surface faces perpendicular to the primary direction; the rotation about "
                "it is arbitrary");
            AaVec e2{0.0, 0.0, 0.0};
            detail::lsq_frame(rFirst.data(), d.data(), e2.data());
        }
        return d;
    };

    // The face pose.
    const std::size_t b1 = best_bin(nullptr);
    const AaVec d1 = aa_refine(f, bins[b1], cos_tol);
    const AaVec d2 = second_dir(d1);
    const AaVec d3 = detail::vec3_cross(d1, d2);
    double face_area = 0.0;
    for (std::size_t t = 0; t < f.mNumReal; ++t) {
        const AaVec& nn = f.mNormal[t];
        const double m =
            std::max({std::abs(detail::vec3_dot(nn, d1)), std::abs(detail::vec3_dot(nn, d2)),
                      std::abs(detail::vec3_dot(nn, d3))});
        face_area += m >= cos_tol ? f.mArea[t] : 0.0;
    }
    result.mFaceShare = face_area / f.mRealArea;

    // The axis pose.
    std::array<double, 9> scatter{};
    for (std::size_t t = 0; t < f.mNumReal; ++t)
        for (std::size_t i = 0; i < 3; ++i)
            for (std::size_t j = 0; j < 3; ++j)
                scatter[3 * i + j] += f.mArea[t] * f.mNormal[t][i] * f.mNormal[t][j];
    std::array<double, 3> ev{};
    std::array<double, 9> vec{};
    detail::lsq_sym3_eigen(scatter.data(), ev.data(), vec.data());
    const AaVec axis{vec[0], vec[1], vec[2]};
    double axis_area = 0.0;
    for (std::size_t t = 0; t < f.mNumReal; ++t) {
        const double c = std::abs(detail::vec3_dot(f.mNormal[t], axis));
        axis_area += (c <= sin_tol || c >= cos_tol) ? f.mArea[t] : 0.0;
    }
    result.mAxisShare = axis_area / f.mRealArea;

    AaVec up, xdir;
    if (result.mAxisShare > result.mFaceShare + 0.05) {
        result.mMode = "axis";
        double plus = 0.0, minus = 0.0;
        for (std::size_t t = 0; t < f.mNormal.size(); ++t) {
            const double c = detail::vec3_dot(f.mNormal[t], axis);
            plus += c >= cos_tol ? f.mArea[t] : 0.0;
            minus += c <= -cos_tol ? f.mArea[t] : 0.0;
        }
        // The end with the larger area facing along the axis goes down.
        up = plus > minus ? detail::vec3_scale(axis, -1.0) : axis;
        xdir = second_dir(up);
        result.mShare = result.mAxisShare;
        result.mPrimary = up;
    } else {
        result.mMode = "face";
        up = detail::vec3_scale(d1, -1.0);
        xdir = d2;
        result.mShare = result.mFaceShare;
        result.mPrimary = d1;
    }
    result.mSecondary = xdir;
    const AaVec ydir = detail::vec3_cross(up, xdir);
    std::array<double, 16> m{xdir[0], xdir[1], xdir[2], 0.0, ydir[0], ydir[1], ydir[2], 0.0,
                             up[0],   up[1],   up[2],   0.0, 0.0,     0.0,     0.0,     1.0};
    // Centre the box on the Z axis, its lowest point at z = 0.
    std::array<double, 3> lo{0, 0, 0}, hi{0, 0, 0};
    for (std::size_t i = 0; i < n; ++i) {
        AaVec p{0.0, 0.0, 0.0};
        for (std::size_t k = 0; k < dim && k < 3; ++k)
            p[k] = detail::read_double(rMesh.Points(), i * dim + k);
        for (std::size_t r = 0; r < 3; ++r) {
            const double v = m[4 * r] * p[0] + m[4 * r + 1] * p[1] + m[4 * r + 2] * p[2];
            lo[r] = i == 0 ? v : std::min(lo[r], v);
            hi[r] = i == 0 ? v : std::max(hi[r], v);
        }
    }
    m[3] = -0.5 * (lo[0] + hi[0]);
    m[7] = -0.5 * (lo[1] + hi[1]);
    m[11] = -lo[2];
    result.mMatrix = m;
    result.mMesh = transform(rMesh, transform_from_matrix(m.data()), rOptions.mRotateVectorData);
    return result;
}

}  // namespace meshioplusplus

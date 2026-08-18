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
#include <cmath>
#include <cstring>

// Project includes
#include "meshioplusplus/detail/decimate_common.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {
namespace detail {

DecimCsr decim_vertex_faces_csr(const DecimFaces& rFaces, std::size_t n) {
    DecimCsr csr;
    csr.mXadj.assign(n + 1, 0);
    const std::size_t nf = rFaces.mNumFaces;
    for (std::size_t i = 0; i < nf * 3; ++i)
        ++csr.mXadj[static_cast<std::size_t>(rFaces.mCorners[i]) + 1];
    for (std::size_t i = 0; i < n; ++i)
        csr.mXadj[i + 1] += csr.mXadj[i];
    csr.mAdj.resize(static_cast<std::size_t>(csr.mXadj[n]));
    std::vector<std::int64_t> cursor(csr.mXadj.begin(), csr.mXadj.end() - 1);
    for (std::size_t f = 0; f < nf; ++f)
        for (std::size_t k = 0; k < 3; ++k) {
            const std::size_t v = static_cast<std::size_t>(rFaces.mCorners[f * 3 + k]);
            csr.mAdj[static_cast<std::size_t>(cursor[v]++)] = static_cast<std::int64_t>(f);
        }
    return csr;
}

void decim_face_planes(const DecimFaces& rFaces, const std::vector<double>& rXyz,
                       std::vector<double>& rQuadK, std::vector<double>& rNormals) {
    const std::size_t nf = rFaces.mNumFaces;
    rQuadK.assign(nf * 10, 0.0);
    rNormals.assign(nf * 3, 0.0);
    parallel_for(nf, [&](std::size_t f) {
        const std::int64_t* c = rFaces.mCorners.data() + f * 3;
        const double* p0 = rXyz.data() + static_cast<std::size_t>(c[0]) * 3;
        const double* p1 = rXyz.data() + static_cast<std::size_t>(c[1]) * 3;
        const double* p2 = rXyz.data() + static_cast<std::size_t>(c[2]) * 3;
        const double e1x = p1[0] - p0[0];
        const double e1y = p1[1] - p0[1];
        const double e1z = p1[2] - p0[2];
        const double e2x = p2[0] - p0[0];
        const double e2y = p2[1] - p0[1];
        const double e2z = p2[2] - p0[2];
        const double nx = e1y * e2z - e1z * e2y;
        const double ny = e1z * e2x - e1x * e2z;
        const double nz = e1x * e2y - e1y * e2x;
        const double len2 = nx * nx + ny * ny + nz * nz;
        if (len2 == 0.0)
            return;  // degenerate: zero quadric, zero normal
        const double len = std::sqrt(len2);
        const double a = nx / len;
        const double b = ny / len;
        const double cc = nz / len;
        const double d = -(a * p0[0] + b * p0[1] + cc * p0[2]);
        const double w = 0.5 * len;
        const double wa = w * a;
        const double wb = w * b;
        const double wc = w * cc;
        const double wd = w * d;
        double* k = rQuadK.data() + f * 10;
        k[0] = wa * a;
        k[1] = wa * b;
        k[2] = wa * cc;
        k[3] = wa * d;
        k[4] = wb * b;
        k[5] = wb * cc;
        k[6] = wb * d;
        k[7] = wc * cc;
        k[8] = wc * d;
        k[9] = wd * d;
        double* nrm = rNormals.data() + f * 3;
        nrm[0] = a;
        nrm[1] = b;
        nrm[2] = cc;
    });
}

std::vector<double> decim_accumulate_quadrics(const DecimCsr& rCsr, std::size_t n,
                                              const std::vector<double>& rQuadK) {
    std::vector<double> q(n * 10, 0.0);
    parallel_for(n, [&](std::size_t v) {
        double* qv = q.data() + v * 10;
        for (std::int64_t k = rCsr.mXadj[v]; k < rCsr.mXadj[v + 1]; ++k) {
            const double* kf =
                rQuadK.data() +
                static_cast<std::size_t>(rCsr.mAdj[static_cast<std::size_t>(k)]) * 10;
            for (int i = 0; i < 10; ++i)
                qv[i] += kf[i];
        }
    });
    return q;
}

void decim_mark_features(const DecimCsr& rCsr, std::size_t n, const std::vector<double>& rNormals,
                         double CosThreshold, std::vector<std::uint8_t>& rPinned) {
    parallel_for(n, [&](std::size_t v) {
        const std::int64_t b = rCsr.mXadj[v];
        const std::int64_t e = rCsr.mXadj[v + 1];
        for (std::int64_t p = b; p < e; ++p) {
            const double* na = rNormals.data() +
                               static_cast<std::size_t>(rCsr.mAdj[static_cast<std::size_t>(p)]) * 3;
            if (na[0] == 0.0 && na[1] == 0.0 && na[2] == 0.0)
                continue;
            for (std::int64_t q = p + 1; q < e; ++q) {
                const double* nb =
                    rNormals.data() +
                    static_cast<std::size_t>(rCsr.mAdj[static_cast<std::size_t>(q)]) * 3;
                if (nb[0] == 0.0 && nb[1] == 0.0 && nb[2] == 0.0)
                    continue;
                if (na[0] * nb[0] + na[1] * nb[1] + na[2] * nb[2] < CosThreshold) {
                    rPinned[v] = 1;
                    return;
                }
            }
        }
    });
}

double decim_quadric_error(const double* q, double x, double y, double z) {
    return q[0] * x * x + q[4] * y * y + q[7] * z * z +
           2.0 * (q[1] * x * y + q[2] * x * z + q[5] * y * z) +
           2.0 * (q[3] * x + q[6] * y + q[8] * z) + q[9];
}

DecimPlaced decim_place(const DecimPlaceCtx& rCtx, std::int64_t a, std::int64_t b) {
    const std::vector<double>& xyz = *rCtx.mpXyz;
    const std::vector<double>& quads = *rCtx.mpQ;
    const double* xa = xyz.data() + static_cast<std::size_t>(a) * 3;
    const double* xb = xyz.data() + static_cast<std::size_t>(b) * 3;
    const double* qa = quads.data() + static_cast<std::size_t>(a) * 10;
    const double* qb = quads.data() + static_cast<std::size_t>(b) * 10;
    double q[10];
    for (int i = 0; i < 10; ++i)
        q[i] = qa[i] + qb[i];

    DecimPlaced out;
    const std::vector<std::uint8_t>& pinned = *rCtx.mpPinned;
    if (pinned[static_cast<std::size_t>(a)]) {
        out.mX[0] = xa[0];
        out.mX[1] = xa[1];
        out.mX[2] = xa[2];
        out.mErr = decim_quadric_error(q, out.mX[0], out.mX[1], out.mX[2]);
        return out;
    }
    if (pinned[static_cast<std::size_t>(b)]) {
        out.mX[0] = xb[0];
        out.mX[1] = xb[1];
        out.mX[2] = xb[2];
        out.mErr = decim_quadric_error(q, out.mX[0], out.mX[1], out.mX[2]);
        return out;
    }

    if (rCtx.mPlacement == DecimatePlacement::Endpoint) {
        const double err_a = decim_quadric_error(q, xa[0], xa[1], xa[2]);
        const double err_b = decim_quadric_error(q, xb[0], xb[1], xb[2]);
        if (err_b < err_a) {  // tie -> a, the lower id
            out.mX[0] = xb[0];
            out.mX[1] = xb[1];
            out.mX[2] = xb[2];
            out.mErr = err_b;
        } else {
            out.mX[0] = xa[0];
            out.mX[1] = xa[1];
            out.mX[2] = xa[2];
            out.mErr = err_a;
        }
        return out;
    }

    bool solved = false;
    if (rCtx.mPlacement == DecimatePlacement::Optimal) {
        // Minimize E: solve A x = -bvec with A the quadric's upper-left 3x3 and
        // bvec = (q3, q6, q8), via the cofactor (adjugate) form. On an exactly
        // planar patch A is singular by construction, so the midpoint fallback
        // is the hot path there -- expected, not a defect.
        const double c00 = q[4] * q[7] - q[5] * q[5];
        const double c01 = q[2] * q[5] - q[1] * q[7];
        const double c02 = q[1] * q[5] - q[2] * q[4];
        const double det = q[0] * c00 + q[1] * c01 + q[2] * c02;
        double scale = std::abs(q[0]);
        scale = std::max(scale, std::abs(q[1]));
        scale = std::max(scale, std::abs(q[2]));
        scale = std::max(scale, std::abs(q[4]));
        scale = std::max(scale, std::abs(q[5]));
        scale = std::max(scale, std::abs(q[7]));
        if (std::abs(det) > 1e-12 * (scale * scale * scale)) {
            const double c11 = q[0] * q[7] - q[2] * q[2];
            const double c12 = q[1] * q[2] - q[0] * q[5];
            const double c22 = q[0] * q[4] - q[1] * q[1];
            const double inv = 1.0 / det;
            out.mX[0] = -(c00 * q[3] + c01 * q[6] + c02 * q[8]) * inv;
            out.mX[1] = -(c01 * q[3] + c11 * q[6] + c12 * q[8]) * inv;
            out.mX[2] = -(c02 * q[3] + c12 * q[6] + c22 * q[8]) * inv;
            solved = true;
        }
    }
    if (!solved) {  // Midpoint, or Optimal's ill-conditioned fallback
        out.mX[0] = (xa[0] + xb[0]) * 0.5;
        out.mX[1] = (xa[1] + xb[1]) * 0.5;
        out.mX[2] = (xa[2] + xb[2]) * 0.5;
    }
    out.mErr = decim_quadric_error(q, out.mX[0], out.mX[1], out.mX[2]);
    return out;
}

void decim_face_normal(const std::vector<double>& rXyz, const std::int64_t* pCorners,
                       std::int64_t A, std::int64_t B, const double* pSub, double pOut[3]) {
    double p[3][3];
    for (int k = 0; k < 3; ++k) {
        const std::int64_t id = pCorners[k];
        if (pSub != nullptr && (id == A || id == B)) {
            p[k][0] = pSub[0];
            p[k][1] = pSub[1];
            p[k][2] = pSub[2];
        } else {
            const double* src = rXyz.data() + static_cast<std::size_t>(id) * 3;
            p[k][0] = src[0];
            p[k][1] = src[1];
            p[k][2] = src[2];
        }
    }
    const double e1x = p[1][0] - p[0][0];
    const double e1y = p[1][1] - p[0][1];
    const double e1z = p[1][2] - p[0][2];
    const double e2x = p[2][0] - p[0][0];
    const double e2y = p[2][1] - p[0][1];
    const double e2z = p[2][2] - p[0][2];
    pOut[0] = e1y * e2z - e1z * e2y;
    pOut[1] = e1z * e2x - e1x * e2z;
    pOut[2] = e1x * e2y - e1y * e2x;
}

std::vector<std::int64_t> decim_vertex_ring(const std::vector<std::int64_t>& rVFaces,
                                            const std::vector<std::int64_t>& rCorners,
                                            std::int64_t v) {
    std::vector<std::int64_t> ring;
    ring.reserve(rVFaces.size() * 2);
    for (std::int64_t f : rVFaces) {
        const std::int64_t* c = rCorners.data() + static_cast<std::size_t>(f) * 3;
        for (int k = 0; k < 3; ++k)
            if (c[k] != v)
                ring.push_back(c[k]);
    }
    std::sort(ring.begin(), ring.end());
    ring.erase(std::unique(ring.begin(), ring.end()), ring.end());
    return ring;
}

std::size_t decim_count_common(const std::vector<std::int64_t>& rA,
                               const std::vector<std::int64_t>& rB) {
    std::size_t i = 0;
    std::size_t j = 0;
    std::size_t common = 0;
    while (i < rA.size() && j < rB.size()) {
        if (rA[i] < rB[j])
            ++i;
        else if (rB[j] < rA[i])
            ++j;
        else {
            ++common;
            ++i;
            ++j;
        }
    }
    return common;
}

void decim_sorted_erase(std::vector<std::int64_t>& rVec, std::int64_t Value) {
    auto it = std::lower_bound(rVec.begin(), rVec.end(), Value);
    if (it != rVec.end() && *it == Value)
        rVec.erase(it);
}

void decim_sorted_insert(std::vector<std::int64_t>& rVec, std::int64_t Value) {
    rVec.insert(std::lower_bound(rVec.begin(), rVec.end(), Value), Value);
}

}  // namespace detail
}  // namespace meshioplusplus

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

// The crease test shared by feature_edges, decimate, decimate_volume and
// smooth. See crease_edges.hpp for the rule.

// System includes
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

// Project includes
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/parallel.hpp"

// Project includes (private, not installed)
#include "crease_edges.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

// One use of an undirected edge by one face ring.
struct CeEdgeUse {
    std::int64_t mLo;
    std::int64_t mHi;
    std::int64_t mFace;
    bool mForward;  // the ring walks it lo -> hi
};

bool ce_use_less(const CeEdgeUse& rA, const CeEdgeUse& rB) {
    if (rA.mLo != rB.mLo)
        return rA.mLo < rB.mLo;
    if (rA.mHi != rB.mHi)
        return rA.mHi < rB.mHi;
    if (rA.mFace != rB.mFace)
        return rA.mFace < rB.mFace;
    return rA.mForward < rB.mForward;
}

}  // namespace

void ring_unit_normal(const double* pXyz, const std::int64_t* pRing, std::size_t Size,
                      double* pOut) {
    Vec3 nrm = {0.0, 0.0, 0.0};
    for (std::size_t k = 0; k < Size; ++k) {
        const double* a = pXyz + static_cast<std::size_t>(pRing[k]) * 3;
        const double* b = pXyz + static_cast<std::size_t>(pRing[(k + 1) % Size]) * 3;
        nrm = vec3_add(nrm, vec3_cross(Vec3{a[0], a[1], a[2]}, Vec3{b[0], b[1], b[2]}));
    }
    const Vec3 unit = vec3_normalize(nrm);
    pOut[0] = unit[0];
    pOut[1] = unit[1];
    pOut[2] = unit[2];
}

std::vector<CreaseEdge> crease_edges(const std::vector<std::int64_t>& rStart,
                                     const std::vector<std::int64_t>& rNodes,
                                     const std::vector<double>& rUnitNormal,
                                     double FeatureAngleDeg) {
    const std::size_t nf = rStart.empty() ? 0 : rStart.size() - 1;

    // Per face, how many non-collapsed ring edges it contributes; then each
    // face writes its uses at its prefix offset.
    std::vector<std::uint64_t> count(nf, 0);
    parallel_for(nf, [&](std::size_t f) {
        const std::int64_t b = rStart[f];
        const std::int64_t s = rStart[f + 1] - b;
        std::uint64_t c = 0;
        for (std::int64_t k = 0; k < s; ++k)
            c += rNodes[static_cast<std::size_t>(b + k)] !=
                         rNodes[static_cast<std::size_t>(b + (k + 1) % s)]
                     ? 1
                     : 0;
        count[f] = s >= 2 ? c : 0;
    });
    std::vector<std::uint64_t> at(nf);
    const std::uint64_t nuses =
        parallel_exclusive_scan(count.data(), nf, at.data(), std::uint64_t{0});
    std::vector<CeEdgeUse> uses(nuses);
    parallel_for(nf, [&](std::size_t f) {
        if (count[f] == 0)
            return;
        const std::int64_t b = rStart[f];
        const std::int64_t s = rStart[f + 1] - b;
        std::uint64_t w = at[f];
        for (std::int64_t k = 0; k < s; ++k) {
            const std::int64_t u = rNodes[static_cast<std::size_t>(b + k)];
            const std::int64_t v = rNodes[static_cast<std::size_t>(b + (k + 1) % s)];
            if (u == v)
                continue;
            const std::int64_t face = static_cast<std::int64_t>(f);
            uses[w++] = u < v ? CeEdgeUse{u, v, face, true} : CeEdgeUse{v, u, face, false};
        }
    });
    parallel_sort(uses.begin(), uses.end(), ce_use_less);

    const double cos_thr = std::cos(FeatureAngleDeg * 3.14159265358979323846 / 180.0);
    const double nan = std::numeric_limits<double>::quiet_NaN();
    std::vector<CreaseEdge> out;
    for (std::size_t b = 0; b < uses.size();) {
        std::size_t e = b + 1;
        while (e < uses.size() && uses[e].mLo == uses[b].mLo && uses[e].mHi == uses[b].mHi)
            ++e;
        CreaseEdge edge;
        edge.mLo = uses[b].mLo;
        edge.mHi = uses[b].mHi;
        edge.mUses = static_cast<std::int32_t>(e - b);
        edge.mAngleDeg = nan;
        bool report = edge.mUses != 2;
        if (edge.mUses == 2) {
            const CeEdgeUse& ux = uses[b];
            const CeEdgeUse& uy = uses[b + 1];
            if (ux.mFace == uy.mFace) {
                b = e;
                continue;  // a ring revisiting its own edge: not a surface edge
            }
            edge.mInconsistent = ux.mForward == uy.mForward;
            const double* nx = &rUnitNormal[static_cast<std::size_t>(ux.mFace) * 3];
            const double* ny = &rUnitNormal[static_cast<std::size_t>(uy.mFace) * 3];
            const bool zx = nx[0] == 0.0 && nx[1] == 0.0 && nx[2] == 0.0;
            const bool zy = ny[0] == 0.0 && ny[1] == 0.0 && ny[2] == 0.0;
            if (!zx && !zy) {
                double dot = nx[0] * ny[0] + nx[1] * ny[1] + nx[2] * ny[2];
                if (edge.mInconsistent)
                    dot = -dot;
                edge.mSharp = dot < cos_thr;
                const double c = dot < -1.0 ? -1.0 : (dot > 1.0 ? 1.0 : dot);
                edge.mAngleDeg = std::acos(c) * (180.0 / 3.14159265358979323846);
            }
            report = edge.mInconsistent || edge.mSharp;
        }
        if (report)
            out.push_back(edge);
        b = e;
    }
    return out;
}

std::vector<CreaseEdge> crease_edges_triangles(const std::vector<std::int64_t>& rCorners,
                                               const std::vector<double>& rUnitNormal,
                                               double FeatureAngleDeg) {
    const std::size_t nf = rCorners.size() / 3;
    std::vector<std::int64_t> start(nf + 1);
    for (std::size_t f = 0; f <= nf; ++f)
        start[f] = static_cast<std::int64_t>(f * 3);
    return crease_edges(start, rCorners, rUnitNormal, FeatureAngleDeg);
}

void pin_crease_endpoints(const std::vector<CreaseEdge>& rEdges,
                          std::vector<std::uint8_t>& rPinned) {
    for (const CreaseEdge& e : rEdges) {
        if (!e.IsCrease())
            continue;
        rPinned[static_cast<std::size_t>(e.mLo)] = 1;
        rPinned[static_cast<std::size_t>(e.mHi)] = 1;
    }
}

}  // namespace detail
}  // namespace meshioplusplus

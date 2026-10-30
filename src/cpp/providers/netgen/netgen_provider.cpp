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

// The Netgen tetrahedral mesher as a provider: the nglib half of
// `tetrahedralize`. Compiled into the core in a monolithic build
// (-DMESHIOPLUSPLUS_WITH_NETGEN=ON) and into the meshioplusplus-netgen add-on
// in a split one. See providers.hpp and operations/tetrahedralize.cpp, which
// keeps the validation, the region bookkeeping, the native octree route and
// the statistics.

#ifdef MESHIOPLUSPLUS_HAS_NETGEN

// System includes
#include <algorithm>
#include <array>
#include <cstdint>
#include <mutex>
#include <stdexcept>
#include <string>
#include <vector>

// External includes
namespace nglib {
#include <nglib.h>
}  // namespace nglib

// Project includes
#include "meshioplusplus/providers.hpp"
#include "netgen_provider.hpp"

// The version CMake read (nglib.h carries none); "unknown" for a hand build.
#ifndef MESHIOPLUSPLUS_NETGEN_VERSION
#define MESHIOPLUSPLUS_NETGEN_VERSION "unknown"
#endif

namespace meshioplusplus {

namespace {

// nglib's state is global: one call at a time, and one Ng_Init.
std::mutex& netgen_mutex() {
    static std::mutex m;
    return m;
}

/// Deletes the mesh on every exit.
struct NgMeshGuard {
    nglib::Ng_Mesh* mMesh;
    explicit NgMeshGuard(nglib::Ng_Mesh* pMesh) : mMesh(pMesh) {}
    NgMeshGuard(const NgMeshGuard&) = delete;
    NgMeshGuard& operator=(const NgMeshGuard&) = delete;
    ~NgMeshGuard() { nglib::Ng_DeleteMesh(mMesh); }
};

TetMeshRegionResult netgen_mesh_region(const TetMeshRegionRequest& rRequest) {
    std::lock_guard<std::mutex> lock(netgen_mutex());
    static std::once_flag s_init;
    std::call_once(s_init, [] { nglib::Ng_Init(); });

    nglib::Ng_Mesh* p_m = nglib::Ng_NewMesh();
    const NgMeshGuard guard(p_m);
    for (const auto& r_p : rRequest.mPoints) {
        double p[3] = {r_p[0], r_p[1], r_p[2]};
        nglib::Ng_AddPoint(p_m, p);
    }
    for (const auto& r_t : rRequest.mTriangles) {
        int pi[3] = {static_cast<int>(r_t[0]) + 1, static_cast<int>(r_t[1]) + 1,
                     static_cast<int>(r_t[2]) + 1};
        nglib::Ng_AddSurfaceElement(p_m, nglib::NG_TRIG, pi);
    }
    nglib::Ng_Meshing_Parameters mp;
    if (rRequest.mMaxSize > 0.0)
        mp.maxh = rRequest.mMaxSize;
    mp.grading = rRequest.mGrading;
    mp.optsteps_3d = rRequest.mOptimizeSteps;
    mp.optvolmeshenable = rRequest.mOptimizeSteps > 0 ? 1 : 0;
    mp.optsurfmeshenable = 0;
    const nglib::Ng_Result status = nglib::Ng_GenerateVolumeMesh(p_m, &mp);
    const int ne = nglib::Ng_GetNE(p_m);
    const int nlocal = nglib::Ng_GetNP(p_m);
    if (status != nglib::NG_OK || ne <= 0)
        throw std::runtime_error("Netgen failed (status " +
                                 std::to_string(static_cast<int>(status)) +
                                 "); is the surface free of self-intersections "
                                 "(find_intersections)?");

    TetMeshRegionResult result;
    result.mPoints.resize(static_cast<std::size_t>(nlocal));
    for (int i = 1; i <= nlocal; ++i) {
        double p[3];
        nglib::Ng_GetPoint(p_m, i, p);
        result.mPoints[static_cast<std::size_t>(i - 1)] = {p[0], p[1], p[2]};
    }

    // The surface kept: its points where they were, its triangles as given.
    const std::size_t n_in = rRequest.mPoints.size();
    for (std::size_t i = 0; i < n_in && result.mBoundaryPreserved; ++i)
        if (i >= result.mPoints.size() || result.mPoints[i] != rRequest.mPoints[i])
            result.mBoundaryPreserved = false;
    const int n_tri = static_cast<int>(rRequest.mTriangles.size());
    if (nglib::Ng_GetNSE(p_m) != n_tri) {
        result.mBoundaryPreserved = false;
    } else {
        std::vector<std::array<std::int64_t, 3>> want, got;
        for (const auto& r_t : rRequest.mTriangles) {
            std::array<std::int64_t, 3> k = r_t;
            std::sort(k.begin(), k.end());
            want.push_back(k);
        }
        for (int e = 1; e <= n_tri; ++e) {
            int pi[NG_SURFACE_ELEMENT_MAXPOINTS];
            nglib::Ng_GetSurfaceElement(p_m, e, pi);
            std::array<std::int64_t, 3> k{};
            for (std::size_t j = 0; j < 3; ++j)
                k[j] = pi[j] >= 1 && pi[j] <= static_cast<int>(n_in) ? pi[j] - 1 : -1;
            std::sort(k.begin(), k.end());
            got.push_back(k);
        }
        std::sort(want.begin(), want.end());
        std::sort(got.begin(), got.end());
        if (want != got)
            result.mBoundaryPreserved = false;
    }

    result.mTets.reserve(static_cast<std::size_t>(ne));
    for (int e = 1; e <= ne; ++e) {
        int pi[NG_VOLUME_ELEMENT_MAXPOINTS];
        if (nglib::Ng_GetVolumeElement(p_m, e, pi) != nglib::NG_TET)
            throw std::runtime_error("Netgen returned an element that is not a tetrahedron");
        result.mTets.push_back({pi[0] - 1, pi[1] - 1, pi[2] - 1, pi[3] - 1});
    }
    return result;
}

}  // namespace

namespace providers {

TetMesherProvider make_netgen_provider() {
    TetMesherProvider provider;
    provider.mName = "netgen";
    provider.mVersion = MESHIOPLUSPLUS_NETGEN_VERSION;
    provider.mMeshRegion = netgen_mesh_region;
    return provider;
}

void register_netgen() {
    register_tet_mesher(make_netgen_provider());
}

}  // namespace providers

}  // namespace meshioplusplus

#endif  // MESHIOPLUSPLUS_HAS_NETGEN

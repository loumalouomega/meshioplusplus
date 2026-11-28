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

/**
 * @file providers/proj/proj_provider.cpp
 * @brief Optional PROJ coordinate reference system provider (MIT; roadmap
 * 9.10.2). Built only when MESHIOPLUSPLUS_WITH_PROJ=ON or as a plugin.
 */

#ifdef MESHIOPLUSPLUS_HAS_PROJ

// System includes
#include <array>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

// External includes
#include <proj.h>

// Project includes
#include "meshioplusplus/abi_version.hpp"
#include "meshioplusplus/providers.hpp"

namespace meshioplusplus {
namespace providers {

namespace {

struct ProjCtxDeleter {
    void operator()(PJ_CONTEXT* pCtx) const noexcept { proj_context_destroy(pCtx); }
};
struct ProjDeleter {
    void operator()(PJ* pPj) const noexcept { proj_destroy(pPj); }
};

std::string proj_error_text(PJ_CONTEXT* pCtx) {
    const int code = proj_context_errno(pCtx);
    const char* p_text = proj_context_errno_string(pCtx, code);
    return p_text ? std::string(p_text) : ("PROJ error " + std::to_string(code));
}

}  // namespace

CrsProvider make_proj_provider() {
    CrsProvider p;
    p.mName = "proj";
    p.mVersion = std::to_string(PROJ_VERSION_MAJOR) + "." + std::to_string(PROJ_VERSION_MINOR) +
                 "." + std::to_string(PROJ_VERSION_PATCH);
    p.mTransform = [](const CrsRequest& rReq) -> std::vector<std::array<double, 3>> {
        // One context per call: PROJ contexts are not thread-safe.
        std::unique_ptr<PJ_CONTEXT, ProjCtxDeleter> ctx(proj_context_create());
        if (!ctx)
            throw std::runtime_error("PROJ: proj_context_create failed");
        std::unique_ptr<PJ, ProjDeleter> raw(proj_create_crs_to_crs(
            ctx.get(), rReq.mSourceCrs.c_str(), rReq.mTargetCrs.c_str(), nullptr));
        if (!raw)
            throw std::invalid_argument("PROJ: cannot build a transformation from '" +
                                        rReq.mSourceCrs + "' to '" + rReq.mTargetCrs +
                                        "': " + proj_error_text(ctx.get()));
        // Easting/longitude first, as pyproj's always_xy=True.
        std::unique_ptr<PJ, ProjDeleter> pj(proj_normalize_for_visualization(ctx.get(), raw.get()));
        if (!pj)
            throw std::invalid_argument("PROJ: cannot normalise the axis order: " +
                                        proj_error_text(ctx.get()));
        std::vector<std::array<double, 3>> out = rReq.mPoints;
        if (out.empty())
            return out;
        constexpr std::size_t stride = sizeof(std::array<double, 3>);
        const std::size_t n = out.size();
        const std::size_t done =
            proj_trans_generic(pj.get(), PJ_FWD, &out[0][0], stride, n, &out[0][1], stride, n,
                               &out[0][2], stride, n, nullptr, 0, 0);
        if (done != n || proj_errno(pj.get()) != 0)
            throw std::invalid_argument("PROJ: transformation failed: " +
                                        proj_error_text(ctx.get()));
        return out;
    };
    return p;
}

}  // namespace providers
}  // namespace meshioplusplus

#ifdef MESHIOPLUSPLUS_PROJ_PLUGIN
extern "C" int meshioplusplus_plugin_abi_version(void) {
    return MESHIOPLUSPLUS_ABI_VERSION;
}
extern "C" void meshioplusplus_register_plugin(void) {
    meshioplusplus::register_crs_provider(meshioplusplus::providers::make_proj_provider());
}
#endif

#endif  // MESHIOPLUSPLUS_HAS_PROJ

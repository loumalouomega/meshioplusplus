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
 * @file providers/manifold/manifold_provider.cpp
 * @brief Optional Manifold surface boolean provider (Apache-2.0; roadmap 9.1.3).
 * Built only when MESHIOPLUSPLUS_WITH_MANIFOLD=ON or as a plugin.
 */

#ifdef MESHIOPLUSPLUS_HAS_MANIFOLD

#include <manifold/manifold.h>

#include <array>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include "meshioplusplus/abi_version.hpp"
#include "meshioplusplus/providers.hpp"

namespace meshioplusplus {
namespace providers {

BooleanProvider make_manifold_provider() {
    BooleanProvider p;
    p.mName = "manifold";
    p.mVersion = "3.x";
    p.mBoolean = [](const BooleanRequest& rReq) -> BooleanReply {
        auto to_manifold = [](const std::vector<std::array<double, 3>>& rPts,
                              const std::vector<std::array<std::int64_t, 3>>& rTris) {
            manifold::MeshGL64 gl;
            gl.numProp = 3;
            gl.vertProperties.reserve(rPts.size() * 3);
            for (const auto& pt : rPts) {
                gl.vertProperties.push_back(pt[0]);
                gl.vertProperties.push_back(pt[1]);
                gl.vertProperties.push_back(pt[2]);
            }
            gl.triVerts.reserve(rTris.size() * 3);
            for (const auto& t : rTris) {
                gl.triVerts.push_back(static_cast<std::uint64_t>(t[0]));
                gl.triVerts.push_back(static_cast<std::uint64_t>(t[1]));
                gl.triVerts.push_back(static_cast<std::uint64_t>(t[2]));
            }
            return manifold::Manifold(gl);
        };
        const manifold::Manifold a = to_manifold(rReq.mPointsA, rReq.mTrianglesA);
        const manifold::Manifold b = to_manifold(rReq.mPointsB, rReq.mTrianglesB);
        const manifold::OpType op = rReq.mOp == 1   ? manifold::OpType::Subtract
                                    : rReq.mOp == 2 ? manifold::OpType::Intersect
                                                    : manifold::OpType::Add;
        const manifold::Manifold res = a.Boolean(b, op);
        const manifold::MeshGL64 out_gl = res.GetMeshGL64();
        BooleanReply rep;
        const std::size_t npts = out_gl.vertProperties.size() / out_gl.numProp;
        rep.mPoints.resize(npts);
        for (std::size_t i = 0; i < npts; ++i)
            rep.mPoints[i] = {out_gl.vertProperties[i * out_gl.numProp + 0],
                              out_gl.vertProperties[i * out_gl.numProp + 1],
                              out_gl.vertProperties[i * out_gl.numProp + 2]};
        const std::size_t ntris = out_gl.triVerts.size() / 3;
        rep.mTriangles.resize(ntris);
        for (std::size_t i = 0; i < ntris; ++i)
            rep.mTriangles[i] = {static_cast<std::int64_t>(out_gl.triVerts[3 * i + 0]),
                                 static_cast<std::int64_t>(out_gl.triVerts[3 * i + 1]),
                                 static_cast<std::int64_t>(out_gl.triVerts[3 * i + 2])};
        return rep;
    };
    return p;
}

}  // namespace providers
}  // namespace meshioplusplus

#ifdef MESHIOPLUSPLUS_MANIFOLD_PLUGIN
extern "C" int meshioplusplus_plugin_abi_version(void) {
    return MESHIOPLUSPLUS_ABI_VERSION;
}
extern "C" void meshioplusplus_register_plugin(void) {
    meshioplusplus::register_boolean_provider(meshioplusplus::providers::make_manifold_provider());
}
#endif

#endif  // MESHIOPLUSPLUS_HAS_MANIFOLD

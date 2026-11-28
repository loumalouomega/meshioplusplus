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
 * @file providers/scotch/scotch_provider.cpp
 * @brief Optional Scotch graph partitioner provider (CeCILL-C; roadmap 9.7.2).
 * Built only when MESHIOPLUSPLUS_WITH_SCOTCH=ON or as a plugin.
 */

#ifdef MESHIOPLUSPLUS_HAS_SCOTCH

#include <cstdio>
#include <scotch.h>

#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include "meshioplusplus/abi_version.hpp"
#include "meshioplusplus/providers.hpp"

namespace meshioplusplus {
namespace providers {

GraphPartitionerProvider make_scotch_provider() {
    GraphPartitionerProvider p;
    p.mName = "scotch";
    p.mVersion = std::to_string(SCOTCH_VERSION) + "." + std::to_string(SCOTCH_RELEASE) + "." +
                 std::to_string(SCOTCH_PATCHLEVEL);
    p.mPartition = [](const GraphPartitionerRequest& rReq) -> std::vector<std::int64_t> {
        SCOTCH_randomSeed(static_cast<SCOTCH_Num>(rReq.mSeed));
        const SCOTCH_Num nvtxs = static_cast<SCOTCH_Num>(rReq.mNumVertices);
        std::vector<SCOTCH_Num> verttab(rReq.mNumVertices + 1);
        for (std::size_t i = 0; i <= rReq.mNumVertices; ++i)
            verttab[i] = static_cast<SCOTCH_Num>(rReq.mpXadj[i]);
        const std::size_t nedges = static_cast<std::size_t>(rReq.mpXadj[rReq.mNumVertices]);
        std::vector<SCOTCH_Num> edgetab(nedges);
        for (std::size_t i = 0; i < nedges; ++i)
            edgetab[i] = static_cast<SCOTCH_Num>(rReq.mpAdjncy[i]);
        std::vector<SCOTCH_Num> velotab;
        SCOTCH_Num* p_velotab = nullptr;
        if (rReq.mpWeights) {
            velotab.resize(rReq.mNumVertices);
            for (std::size_t i = 0; i < rReq.mNumVertices; ++i)
                velotab[i] = static_cast<SCOTCH_Num>(rReq.mpWeights[i]);
            p_velotab = velotab.data();
        }
        SCOTCH_Graph graph;
        if (SCOTCH_graphInit(&graph) != 0)
            throw std::runtime_error("SCOTCH_graphInit failed");
        if (SCOTCH_graphBuild(&graph, 0, nvtxs, verttab.data(), nullptr, p_velotab, nullptr,
                              static_cast<SCOTCH_Num>(nedges), edgetab.data(), nullptr) != 0) {
            SCOTCH_graphExit(&graph);
            throw std::runtime_error("SCOTCH_graphBuild failed");
        }
        SCOTCH_Strat strat;
        SCOTCH_stratInit(&strat);
        std::vector<SCOTCH_Num> part(rReq.mNumVertices, 0);
        const int rc =
            SCOTCH_graphPart(&graph, static_cast<SCOTCH_Num>(rReq.mNParts), &strat, part.data());
        SCOTCH_stratExit(&strat);
        SCOTCH_graphExit(&graph);
        if (rc != 0)
            throw std::runtime_error("SCOTCH_graphPart failed with code " + std::to_string(rc));
        std::vector<std::int64_t> out(rReq.mNumVertices);
        for (std::size_t i = 0; i < rReq.mNumVertices; ++i)
            out[i] = static_cast<std::int64_t>(part[i]);
        return out;
    };
    return p;
}

}  // namespace providers
}  // namespace meshioplusplus

#ifdef MESHIOPLUSPLUS_SCOTCH_PLUGIN
extern "C" int meshioplusplus_plugin_abi_version(void) {
    return MESHIOPLUSPLUS_ABI_VERSION;
}
extern "C" void meshioplusplus_register_plugin(void) {
    meshioplusplus::add_graph_partitioner(meshioplusplus::providers::make_scotch_provider());
}
#endif

#endif  // MESHIOPLUSPLUS_HAS_SCOTCH

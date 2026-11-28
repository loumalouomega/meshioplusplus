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
 * @file providers/metis/metis_provider.cpp
 * @brief Optional METIS graph partitioner provider (Apache-2.0 since METIS 5.1.1;
 * roadmap 9.7.2). Built only when MESHIOPLUSPLUS_WITH_METIS=ON or as a plugin.
 */

#ifdef MESHIOPLUSPLUS_HAS_METIS

#include <metis.h>

#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include "meshioplusplus/abi_version.hpp"
#include "meshioplusplus/providers.hpp"

namespace meshioplusplus {
namespace providers {

GraphPartitionerProvider make_metis_provider() {
    GraphPartitionerProvider p;
    p.mName = "metis";
    p.mVersion = std::to_string(METIS_VER_MAJOR) + "." + std::to_string(METIS_VER_MINOR) + "." +
                 std::to_string(METIS_VER_SUBMINOR);
    p.mPartition = [](const GraphPartitionerRequest& rReq) -> std::vector<std::int64_t> {
        idx_t nvtxs = static_cast<idx_t>(rReq.mNumVertices);
        idx_t ncon = 1;
        idx_t nparts = static_cast<idx_t>(rReq.mNParts);
        std::vector<idx_t> xadj(rReq.mNumVertices + 1);
        for (std::size_t i = 0; i <= rReq.mNumVertices; ++i)
            xadj[i] = static_cast<idx_t>(rReq.mpXadj[i]);
        const std::size_t nedges = static_cast<std::size_t>(rReq.mpXadj[rReq.mNumVertices]);
        std::vector<idx_t> adjncy(nedges);
        for (std::size_t i = 0; i < nedges; ++i)
            adjncy[i] = static_cast<idx_t>(rReq.mpAdjncy[i]);
        std::vector<idx_t> vwgt;
        idx_t* p_vwgt = nullptr;
        if (rReq.mpWeights) {
            vwgt.resize(rReq.mNumVertices);
            for (std::size_t i = 0; i < rReq.mNumVertices; ++i)
                vwgt[i] = static_cast<idx_t>(rReq.mpWeights[i]);
            p_vwgt = vwgt.data();
        }
        idx_t options[METIS_NOPTIONS];
        METIS_SetDefaultOptions(options);
        options[METIS_OPTION_SEED] = static_cast<idx_t>(rReq.mSeed);
        options[METIS_OPTION_UFACTOR] = static_cast<idx_t>(rReq.mImbalance * 1000.0);
        idx_t objval = 0;
        std::vector<idx_t> part(rReq.mNumVertices, 0);
        const int rc =
            METIS_PartGraphKway(&nvtxs, &ncon, xadj.data(), adjncy.data(), p_vwgt, nullptr, nullptr,
                                &nparts, nullptr, nullptr, options, &objval, part.data());
        if (rc != METIS_OK)
            throw std::runtime_error("METIS_PartGraphKway failed with code " + std::to_string(rc));
        std::vector<std::int64_t> out(rReq.mNumVertices);
        for (std::size_t i = 0; i < rReq.mNumVertices; ++i)
            out[i] = static_cast<std::int64_t>(part[i]);
        return out;
    };
    return p;
}

}  // namespace providers
}  // namespace meshioplusplus

#ifdef MESHIOPLUSPLUS_METIS_PLUGIN
extern "C" int meshioplusplus_plugin_abi_version(void) {
    return MESHIOPLUSPLUS_ABI_VERSION;
}
extern "C" void meshioplusplus_register_plugin(void) {
    meshioplusplus::add_graph_partitioner(meshioplusplus::providers::make_metis_provider());
}
#endif

#endif  // MESHIOPLUSPLUS_HAS_METIS

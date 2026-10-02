// SPDX-License-Identifier: MIT
#pragma once

#include "pugixml.hpp"
#include "meshioplusplus/detail/xdmf_common.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus::xdmfdetail {

// Private: pugixml must not appear in installed headers. Both XDMF writers use
// the same region encoding and heavy-data store, including empty named sets.
void xdmf_write_sets(pugi::xml_node grid, xdmfcommon::DataItemStore& rStore, const Mesh& rMesh);

}  // namespace meshioplusplus::xdmfdetail

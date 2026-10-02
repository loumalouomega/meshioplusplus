// SPDX-License-Identifier: MIT
#include "meshioplusplus/formats/vtp.hpp"
#include "../detail/vtk_xml_pieces.hpp"

namespace meshioplusplus {

Mesh read_vtp(const std::string& rPath, const ReadOptions& rOpts) {
    return detail::vtk_xml_read_pieces(rPath, rOpts, "PolyData", "VTP");
}

MeshMetadata read_vtp_metadata(const std::string& rPath, const ReadOptions&) {
    return detail::vtk_xml_pieces_metadata(rPath, "PolyData", "VTP");
}

}  // namespace meshioplusplus

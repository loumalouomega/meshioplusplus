// SPDX-License-Identifier: MIT
#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <vector>
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/read_options.hpp"

namespace meshioplusplus::detail {

// Core-private shared reader for serial PolyData and structured XML pieces.
Mesh vtk_xml_read_pieces(const std::string& rPath, const ReadOptions& rOpts, const char* pType,
                         const char* pFormat);
MeshMetadata vtk_xml_pieces_metadata(const std::string& rPath, const char* pType,
                                     const char* pFormat);

// Point dimensions, x-fastest. Legacy structured datasets also use lines/quads.
void vtk_structured_cells(const std::array<std::int64_t, 3>& rDims,
                          std::vector<std::int64_t>& rConn, std::vector<std::int64_t>& rOffsets,
                          std::vector<std::int64_t>& rTypes, bool VolumeOnly = false);

}  // namespace meshioplusplus::detail

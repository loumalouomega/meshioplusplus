// SPDX-License-Identifier: MIT
#pragma once

// Core-private: pugixml and payload views never cross an installed API.
// Implemented with the VTU reader, the original owner of this decoder.
#include <cstddef>
#include <optional>
#include <string>
#include <string_view>
#include "pugixml.hpp"
#include "meshioplusplus/detail/file_source.hpp"
#include "meshioplusplus/detail/vtu_binary.hpp"
#include "meshioplusplus/ndarray.hpp"

namespace meshioplusplus::detail {

struct VtuContext {
    VtkCodec mCodec = VtkCodec::None;
    std::size_t mHeaderSize = 4;
    bool mBigEndian = false;
    const unsigned char* mRaw = nullptr;
    std::size_t mRawLen = 0;
    const char* mBase64 = nullptr;
    std::size_t mBase64Len = 0;
};

// Owns both XML and raw file bytes. Contexts borrow it for the duration of a
// read; every decoded NDArray owns its storage independently of this source.
struct VtuSource {
    pugi::xml_document mDoc;
    std::optional<FileSource> mFile;
    std::string_view mBytes;
    std::size_t mRawStart = 0;
    std::size_t mRawStop = 0;
    bool mIsRaw = false;
};

void vtu_load(const std::string& rPath, unsigned int ParseOptions, VtuSource& rSource,
              const char* pType = "UnstructuredGrid", const char* pFormat = "VTU");
VtuContext vtk_xml_read_context(const VtuSource& rSource, const char* pFormat);
NDArray vtu_read_data_array(const pugi::xml_node& rDa, const VtuContext& rCtx, int& rNumComponents);

}  // namespace meshioplusplus::detail

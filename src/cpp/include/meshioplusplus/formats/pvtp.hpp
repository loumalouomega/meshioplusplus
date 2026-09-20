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
#pragma once

/**
 * @file formats/pvtp.hpp
 * @brief VTK XML parallel polygonal data (`.pvtp`): `formats/pvtu.hpp`'s index
 * over `.vtp` pieces (v15.0.0).
 *
 * `<VTKFile type="PPolyData"><PPolyData GhostLevel="N">...<Piece Source=
 * "stem/stem_0000.vtp"/>`. Everything in `formats/pvtu.hpp` applies -- carving
 * by `partition:part`, the declaration check, ghosts, merging with one region
 * per piece -- with `.vtp` in place of `.vtu`.
 */

// System includes
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/detail/vtu_binary.hpp"
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/formats/pvtu.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/read_options.hpp"

namespace meshioplusplus {

/// @copydoc write_pvtu
MESHIOPLUSPLUS_API void write_pvtp(const std::string& rPath, const Mesh& rMesh, bool binary = true,
                                   bool zlib = true);

/// @copydoc write_pvtu_codec
MESHIOPLUSPLUS_API void write_pvtp_codec(const std::string& rPath, const Mesh& rMesh, bool binary,
                                         detail::VtkCodec codec,
                                         const std::string& rPartKey = kPvtuPartKey);

/// @copydoc write_pvtu_pieces_codec
MESHIOPLUSPLUS_API void write_pvtp_pieces_codec(const std::string& rPath,
                                                const std::vector<const Mesh*>& rPieces,
                                                bool binary, detail::VtkCodec codec);

/// @copydoc read_pvtu
MESHIOPLUSPLUS_API Mesh read_pvtp(const std::string& rPath, const ReadOptions& rOpts = {});

/// @copydoc read_pvtu_metadata
MESHIOPLUSPLUS_API MeshMetadata read_pvtp_metadata(const std::string& rPath,
                                                   const ReadOptions& rOpts = {});

}  // namespace meshioplusplus

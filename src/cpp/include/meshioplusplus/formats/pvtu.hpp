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
 * @file formats/pvtu.hpp
 * @brief VTK XML parallel unstructured grid (`.pvtu`): an index that declares
 * the arrays and names one `.vtu` piece per part (v14.1.0, roadmap §1.1).
 *
 * A `.pvtu` is `<VTKFile type="PUnstructuredGrid"><PUnstructuredGrid
 * GhostLevel="N"><PPointData/><PCellData/><PPoints/><Piece Source="..."/>...`.
 * The index carries no geometry, only the *declarations* of the arrays every
 * piece holds (name, type, component count) and the relative path of each piece.
 * `.pvtp` (`formats/pvtp.hpp`) is the same over `.vtp` pieces.
 *
 * ### Write
 *
 * `write_pvtu_pieces_codec` is the primitive: it takes already-carved pieces --
 * `partition`'s output, halo layers included -- and writes each as a `.vtu` in a
 * sibling directory named after the index's stem (`<stem>/<stem>_0000.vtu`,
 * zero-padded like a `{step}` pattern, relative `Source=` paths). **Every piece
 * must declare identical arrays** (name, type, component count); that is
 * validated before anything is created, so a refusal leaves nothing on disk.
 * A piece with no cells is legal and written (an idle rank): skipping it would
 * renumber the parts.
 *
 * `write_pvtu_codec` takes one mesh and carves it by the integer `cell_data`
 * array named by its part key (`partition:part`, what `partition_labels`
 * produces); without that array the mesh is one piece.
 *
 * A `partition:ghost` array (`partition`'s halo tag: 0 owned, L reached at layer
 * L) becomes `vtkGhostType` on cells (`DUPLICATECELL` for any layer) and on
 * points (`DUPLICATEPOINT` when no owned cell of the piece uses the point), and
 * sets the index's `GhostLevel` to the deepest layer. `partition:ghost` itself is
 * written too, since `vtkGhostType` collapses layer 2 onto layer 1. An array the
 * caller already named `vtkGhostType` is passed through unchanged.
 *
 * ### Read
 *
 * `read_pvtu` reads every piece and combines them with `merge()` (no welding: a
 * partition's interface points appear once per piece; `clean` with `weld=true`
 * fuses them), with one `RegionKind::Cell` region per piece (`piece_0`,
 * `piece_1`, ...). `ReadOptions::mPiece` keeps one piece instead. A file with a
 * single piece reads as that piece, with no region. Piece paths are resolved
 * against the index's own directory.
 *
 * Ghost cells are kept by default -- a reader must not silently discard data --
 * and `GhostPolicy::Drop` removes every cell with a `vtkGhostType` bit set (and
 * the points only they used) before merging, which reconstructs the partition
 * of unity a halo'd `partition` broke.
 */

// System includes
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/detail/vtu_binary.hpp"
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/read_options.hpp"

namespace meshioplusplus {

/// What `read_pvtu`/`read_pvtp`/`read_pvd` do with ghost cells.
enum class GhostPolicy {
    Keep,  ///< Leave `vtkGhostType` cells in the mesh (the default).
    Drop,  ///< Remove cells with a ghost bit set, and the points only they used.
};

/**
 * Options specific to the parallel-index readers.
 *
 * A struct of its own, passed as a defaulted trailing parameter, rather than a
 * new `ReadOptions` member: growing `ReadOptions` is an ABI-tier-A layout change
 * that also reaches the aggregates embedding it.
 */
struct PvtuReadOptions {
    GhostPolicy mGhosts = GhostPolicy::Keep;
};

/// The integer `cell_data` array `write_pvtu*` carves a mesh by (`partition_labels`'s).
inline constexpr const char* kPvtuPartKey = "partition:part";

/**
 * @brief Write a mesh as a `.pvtu` index plus one `.vtu` piece per part.
 * @param rPath the output `.pvtu` path; pieces land in `<stem>/` beside it.
 * @param rMesh the mesh; carved by `partition:part` when that array exists.
 * @param binary base64-encode each piece's arrays instead of writing text.
 * @param zlib compress each piece's binary blocks. Ignored when @p binary is false.
 * @throws WriteError on an unwritable path or a non-integer part array.
 */
MESHIOPLUSPLUS_API void write_pvtu(const std::string& rPath, const Mesh& rMesh, bool binary = true,
                                   bool zlib = true);

/**
 * @brief `write_pvtu` with an explicit block codec and part key.
 * @param rPartKey the integer `cell_data` array to carve by; empty (or absent
 *        from the mesh) writes one piece.
 */
MESHIOPLUSPLUS_API void write_pvtu_codec(const std::string& rPath, const Mesh& rMesh, bool binary,
                                         detail::VtkCodec codec,
                                         const std::string& rPartKey = kPvtuPartKey);

/**
 * @brief Write already-carved pieces (e.g. `partition`'s) as an index plus one file each.
 * @throws WriteError on no pieces, or when the pieces do not declare identical
 *         arrays (the message names the array and both pieces).
 */
MESHIOPLUSPLUS_API void write_pvtu_pieces_codec(const std::string& rPath,
                                                const std::vector<const Mesh*>& rPieces,
                                                bool binary, detail::VtkCodec codec);

/**
 * @brief Read a `.pvtu`: every piece merged (one region each), or one piece.
 * @param rOpts selective-read options, forwarded to each piece; `mPieceSet`
 *        selects one piece (`ResolvePiece`).
 * @param rGhost what to do with ghost cells.
 * @throws ReadError on an unparsable index, a missing piece file, or a piece
 *         that is not `.vtu`/`.vtp`.
 */
MESHIOPLUSPLUS_API Mesh read_pvtu(const std::string& rPath, const ReadOptions& rOpts = {},
                                  const PvtuReadOptions& rGhost = {});

/// The sum of the pieces' own metadata, in the order a real read would produce.
MESHIOPLUSPLUS_API MeshMetadata read_pvtu_metadata(const std::string& rPath,
                                                   const ReadOptions& rOpts = {});

}  // namespace meshioplusplus

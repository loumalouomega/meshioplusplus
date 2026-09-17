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
 * @file formats/vtm.hpp
 * @brief VTK XML MultiBlock (`.vtm`): an index file plus one `.vtu` piece per
 * cell block (v11.6.0, roadmap §1 tier B4, part 3 of 3).
 *
 * A `.vtm` file is `<VTKFile type="vtkMultiBlockDataSet"><vtkMultiBlockDataSet>
 * <Block index="0"><DataSet index="i" name="..." file="stem/stem_i.vtu"/>...
 * </Block></vtkMultiBlockDataSet></VTKFile>` -- the index never carries
 * geometry itself, only a list of piece files, one per meshio++ CellBlock.
 *
 * ### Write
 *
 * `write_vtm` creates a directory next to the index (named after the index's
 * own stem) and writes one `.vtu` piece per `CellBlock`: the piece carries
 * that block's cells and cell_data, and the mesh's full point_data, then is
 * pruned with `clean(remove_orphans=true)` so each piece is self-contained
 * (only the points that block actually references). Points are therefore
 * duplicated across pieces sharing a boundary -- documented, not a bug: a
 * `.vtm` piece is a standalone `.vtu` by design, readable on its own.  Each
 * piece is named `"block_<i>"` in the index's `name=` attribute.
 *
 * ### Read
 *
 * `read_vtm` parses the index, reads every `DataSet` piece with `read_vtu`,
 * and combines them with `operations/merge.hpp`'s `merge()` (no welding, so
 * no piece's own points move or fuse with another's -- multiblock pieces are
 * pre-separated by construction, not a set of coincident-point fragments to
 * weld back together). Each piece's cells become one `RegionKind::Cell`
 * region in the merged mesh, named from the index's `name=` attribute (via
 * `merge()`'s own per-input cell index map, so the region is correct
 * regardless of whether same-typed blocks from different pieces ended up
 * consolidated into one output CellBlock).
 *
 * ### Deliberately not supported
 *
 * Nested `<Block>` elements are read structurally (every `DataSet` anywhere
 * under `<vtkMultiBlockDataSet>` is collected, in document order) but nesting
 * itself is not reproduced in the meshio++ mesh -- there is nothing in the
 * uniform API to hold a block hierarchy, only a flat list of named regions.
 * `vtkPolyData` pieces (as opposed to `vtkUnstructuredGrid`) are read via the
 * extension of their own `file=` attribute, so a hand-written `.vtm` mixing
 * `.vtu` and `.vtp` pieces round-trips on read even though this writer always
 * emits `.vtu`.
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/detail/vtu_binary.hpp"
#include "meshioplusplus/read_options.hpp"

namespace meshioplusplus {

/**
 * @brief Write a mesh as VTK XML MultiBlock: an index plus one `.vtu` piece
 *        per cell block.
 * @param rPath the output `.vtm` path; pieces land in a sibling directory
 *        named after its stem (`<stem>/`).
 * @param rMesh the mesh; each `CellBlock` becomes one piece.
 * @param binary base64-encode each piece's arrays instead of writing text.
 * @param zlib compress each piece's binary blocks. Ignored when @p binary is
 *        false.
 * @throws WriteError if the index or a piece file cannot be opened.
 */
MESHIOPLUSPLUS_API void write_vtm(const std::string& rPath, const Mesh& rMesh, bool binary = true,
                                  bool zlib = true);

/**
 * @brief Write a VTK XML MultiBlock with an explicit per-piece block codec.
 * @param rPath the output path.
 * @param rMesh the mesh.
 * @param binary base64-encode each piece's arrays instead of writing text.
 * @param codec the block compressor; `None` writes uncompressed base64.
 * @throws WriteError as `write_vtm`, and when @p codec is not in this build.
 */
MESHIOPLUSPLUS_API void write_vtm_codec(const std::string& rPath, const Mesh& rMesh, bool binary,
                                        detail::VtkCodec codec);

/**
 * @brief Read a VTK XML MultiBlock file.
 * @param rPath the input `.vtm` path.
 * @param rOpts selective-read options, forwarded to each piece's `read_vtu`;
 *        `mPointsOnly` and `mDataArrays` apply.
 * @return the pieces merged into one mesh (no welding), with one
 *         `RegionKind::Cell` region per piece.
 * @throws ReadError if the index cannot be parsed, or a piece cannot be read
 *         (a non-`.vtu`/`.vtp` piece, in particular).
 */
MESHIOPLUSPLUS_API Mesh read_vtm(const std::string& rPath, const ReadOptions& rOpts = {});

/**
 * @brief Summarize a VTK XML MultiBlock file: the sum of its pieces' own
 *        metadata, read without materializing any piece's arrays.
 *
 * `mCellBlocks` reflects the same first-seen, same-type-consolidation order
 * `read_vtm`'s own `merge()` call would produce, so a caller can rely on it
 * agreeing with a real read.
 */
MESHIOPLUSPLUS_API MeshMetadata read_vtm_metadata(const std::string& rPath,
                                                  const ReadOptions& rOpts = {});

}  // namespace meshioplusplus

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
 * @file vtp.hpp
 * @brief VTK XML PolyData (.vtp) C++ reader/writer.
 *
 * The same VTK-XML container as VTU (shared `<DataArray>` machinery in
 * `detail/vtk_xml.hpp` + base64/zlib framing in `detail/vtu_binary.hpp`),
 * with a `<PolyData>` grid holding `<Verts>/<Lines>/<Polys>/<Strips>`
 * connectivity+offsets sections instead of `<Cells>`. Only surface cells are
 * representable: `vertex` (Verts), `line` (Lines), and
 * `triangle`/`quad`/`polygon` (Polys). Cell data follows VTK's canonical
 * PolyData cell order — Verts, then Lines, then Polys, then Strips — in
 * both directions. Triangle strips, poly-vertex/poly-line rows, multiple
 * pieces, appended data, and lzma compression raise (Python fallback).
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
 * @brief Write a mesh as VTK XML PolyData (.vtp).
 *
 * Writable cell types: `vertex` -> Verts, `line` -> Lines, and
 * `triangle`/`quad`/`polygon` (rectangular or ragged polygon blocks) ->
 * Polys. Blocks are emitted grouped in VTK's canonical PolyData cell order
 * (Verts, Lines, Polys), and cell_data is reordered in lockstep. Points are
 * always written with three components (z padded with 0 for 2D meshes).
 *
 * @param rPath filesystem path of the output `.vtp` file
 * @param rMesh the mesh to write
 * @param binary true for base64 "binary" DataArrays, false for ASCII
 * @param zlib compress binary DataArrays with zlib (requires a zlib build;
 *             ignored for ASCII)
 * @throws WriteError if the file cannot be opened or the mesh contains a
 *         cell type PolyData cannot hold (volume or quadratic cells,
 *         polyhedra)
 */
MESHIOPLUSPLUS_API void write_vtp(const std::string& rPath, const Mesh& rMesh, bool binary, bool zlib);

/**
 * @brief Write a `.vtp` choosing the block-compression codec explicitly.
 *
 * `write_vtp(..., zlib)` is exactly this with `Zlib`/`None`, and zlib remains
 * the default everywhere -- nothing changes for existing callers or files.
 *
 * `VtkCodec::LZ4` writes `vtkLZ4DataCompressor` in LZ4's raw block format,
 * which VTK and ParaView read. `VtkCodec::ZSTD` writes
 * `vtkZSTDDataCompressor`, which is a **meshio++ extension**: VTK ships no
 * ZSTD compressor, so such a file is readable by meshio++ but not by ParaView.
 *
 * @throws WriteError naming the CMake option when the codec was not compiled
 *         into this build.
 */
MESHIOPLUSPLUS_API void write_vtp_codec(const std::string& rPath, const Mesh& rMesh, bool binary,
                     detail::VtkCodec codec);

/**
 * @brief Read a VTK XML PolyData (.vtp) file.
 *
 * Verts rows become `vertex` cells, Lines rows `line` cells, and Polys rows
 * `triangle`/`quad`/`polygon` cells (grouped by row size); cell_data is
 * split per block in VTK's canonical Verts/Lines/Polys order.
 *
 * @param rPath filesystem path of the `.vtp` file
 * @param rOpts optional narrowing of what is materialized; the default reads
 *        everything. Unwanted `<PointData>`/`<CellData>` arrays are skipped
 *        before their payload is decoded.
 * @return the read Mesh
 * @throws ReadError on malformed XML, a non-PolyData file, triangle strips,
 *         poly-vertex/poly-line rows, multiple pieces, appended data, or
 *         lzma compression (all deferred to the Python reader)
 */
MESHIOPLUSPLUS_API Mesh read_vtp(const std::string& rPath, const ReadOptions& rOpts = {});

/**
 * @brief Summarize a `.vtp` without decoding its heavy arrays.
 *
 * PolyData has no `types` array -- cell types follow from each section's
 * per-cell size -- so only the sections' `offsets` (one value per cell) are
 * decoded and the connectivity is skipped entirely.
 *
 * @param rPath filesystem path to summarize
 * @return the summary; `mHasBBox` is false (see `read_vtu_metadata`)
 * @throws ReadError on the same unsupported constructs as `read_vtp`
 */
MESHIOPLUSPLUS_API MeshMetadata read_vtp_metadata(const std::string& rPath, const ReadOptions& rOpts = {});

}  // namespace meshioplusplus

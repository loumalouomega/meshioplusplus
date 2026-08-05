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
 * @file formats/vti.hpp
 * @brief VTK XML ImageData (`.vti`): a regular lattice whose geometry is three
 * attributes rather than a point array.
 *
 * ### Why this format and not a raw dump
 *
 * `compute_sdf` and `voxelize` produce a lattice, and a lattice written as an
 * explicit `.vtu` spends the overwhelming majority of its bytes re-stating an
 * index formula. ImageData states it instead: `Origin`, `Spacing` and
 * `WholeExtent` **are** the grid header, which makes `.vti` the one format in
 * meshio++ that round-trips a generated grid's geometry exactly. (No format
 * persists arbitrary `field_data`, so the `sdf:*` keys do not survive any write
 * -- here they do not need to, because the geometry itself carries the same
 * information and `detail::lattice_from_mesh` recovers it.)
 *
 * The container is the same VTK XML this repo already reads and writes, so the
 * `<DataArray>` codec (`detail/vtk_xml.hpp`) and the base64 + block-compression
 * framing (`detail/vtu_binary.hpp`) are reused **verbatim** rather than
 * reimplemented -- that is an argument about test surface, not about typing.
 *
 * ### The mesh side of the deal
 *
 * A `Mesh` has no implicit geometry, so:
 *
 * - **`read_vti` expands** the extent into explicit points and one `hexahedron`
 *   cell block, through `detail/grid_lattice.hpp` -- the same numbering `grid()`
 *   and `voxelize()` produce, which is what makes `read_vti(write_vti(m)) == m`
 *   an identity rather than a coincidence.
 * - **`write_vti` requires a lattice.** A mesh that is not one has no `Origin`/
 *   `Spacing`/`WholeExtent` to write and raises `WriteError` by name. That
 *   includes a *partial* lattice (`voxelize`'s `surface`/`inside` fills, or an
 *   octree): ImageData cannot express a hole, and silently filling one in would
 *   write a different mesh than the caller handed over.
 *
 * ### Deliberately not supported (both raise, so a shim falls back to Python)
 *
 * - `<AppendedData>` -- the VTU C++ reader declines it too, for the same reason.
 * - More than one `<Piece>`, or a piece whose `Extent` is not the `WholeExtent`.
 * - `header_type="UInt64"` is supported on read (the header size is honoured);
 *   the writer always emits the default `UInt32`, as the VTU writer does.
 * - lzma, and any codec this build was compiled without -- by name.
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/detail/vtu_binary.hpp"
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/read_options.hpp"

namespace meshioplusplus {

/**
 * @brief Write a mesh as VTK XML ImageData.
 * @param rPath the output path.
 * @param rMesh the mesh; must be a dense lattice (see the file docs).
 * @param binary base64-encode the arrays instead of writing them as text.
 * @param zlib compress the binary blocks. Ignored when @p binary is false.
 * @throws WriteError when @p rMesh is not a dense lattice, or when zlib was
 *         requested and this build has none.
 */
MESHIOPLUSPLUS_API void write_vti(const std::string& rPath, const Mesh& rMesh, bool binary = true,
                                  bool zlib = true);

/**
 * @brief Write a mesh as VTK XML ImageData with an explicit block codec.
 * @param rPath the output path.
 * @param rMesh the mesh; must be a dense lattice.
 * @param binary base64-encode the arrays instead of writing them as text.
 * @param codec the block compressor; `None` writes uncompressed base64.
 * @throws WriteError as `write_vti`, and when @p codec is not in this build.
 */
MESHIOPLUSPLUS_API void write_vti_codec(const std::string& rPath, const Mesh& rMesh, bool binary,
                                        detail::VtkCodec codec);

/**
 * @brief Read a VTK XML ImageData file, expanding it into explicit cells.
 * @param rPath the input path.
 * @param rOpts selective-read options; `mPointsOnly` and `mDataArrays` apply.
 * @return a mesh with one `hexahedron` block in `detail/grid_lattice.hpp`'s
 *         numbering, or a point-only mesh when the extent has no cells.
 * @throws ReadError on a construct the C++ reader declines (see the file docs).
 */
MESHIOPLUSPLUS_API Mesh read_vti(const std::string& rPath, const ReadOptions& rOpts = {});

/**
 * @brief Summarize a VTK XML ImageData file without decoding its arrays.
 *
 * Unusually cheap even by metadata-reader standards: `WholeExtent` gives both the
 * point and the cell count outright, so nothing at all is decoded.
 */
MESHIOPLUSPLUS_API MeshMetadata read_vti_metadata(const std::string& rPath,
                                                  const ReadOptions& rOpts = {});

}  // namespace meshioplusplus

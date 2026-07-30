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
 * @file cgns.hpp
 * @brief CGNS (.cgns) C++ reader/writer, a genuine CGNS/SIDS-compliant
 *        unstructured-mesh subset stored in HDF5 (the ADF-over-HDF5 mapping
 *        every real CGNS tool uses), readable by cgnslib/ParaView/VTK.
 *
 * Rewritten in v9.8.0 from a private tetrahedra-only encoding that wrote no
 * CGNS node attributes at all and (for any non-tetra mesh) an HDF5 file its
 * own reader rejected on read-back. **The `" data"` leading-space dataset
 * name is not an ad hoc convention** — `#define D_DATA " data"` in cgnslib's
 * own `ADFH.c` is the real ADF-over-HDF5 mapping every node's payload uses;
 * the previous docs here and in doc/formats/cgns.md claiming otherwise were
 * wrong. See doc/formats/cgns.md for the full on-disk layout, the supported/
 * rejected cell-type table, and what CI can and cannot verify.
 *
 * On-disk tree (unstructured, one zone):
 * `/` (root node attrs + `" format"`/`" hdf5version"`) -> `CGNSLibraryVersion`
 * -> `Base` (`CGNSBase_t`, `" data"` = `[CellDim, PhysDim]`) -> `Zone1`
 * (`Zone_t`, `" data"` shape `(3,1)` = `[NVertex, NCell, 0]`) -> `ZoneType`
 * (`"Unstructured"`) + `GridCoordinates` (`CoordinateX`/`Y`/`Z`, the last
 * only when `PhysDim == 3`) + one `Elements_t` section per meshio++ cell
 * block (`ElementRange`/`ElementConnectivity`, 1-based, row-major/
 * element-major unlike MED's Fortran order). Every node carries CGNS's four
 * standard attributes (`name`/`label`/`type`/`flags`) and the file/every
 * group is created with HDF5 link+attribute creation-order tracking, both
 * load-bearing: cgnslib's `has_child`/`has_data` iterate creation order with
 * no name-order fallback.
 *
 * Since v9.9.0 the zone also carries `FlowSolution_t` nodes for point/cell
 * data — `GridLocation` `Vertex` for `point_data`, `CellCenter` for
 * `cell_data`, one `DataArray_t` per scalar. CGNS has no component concept
 * (no `NumberOfComponents` in the SIDS), so a k-component array is split into
 * `<name>_0..<name>_{k-1}` and re-joined on read; `cell_data` needs every
 * block at the zone's `CellDim`, and solutions are read only for a
 * single-zone file. See doc/formats/cgns.md's "Data mapping".
 *
 * A file with no `CGNSBase_t` node (the pre-v9.8.0 layout, or upstream
 * meshio's own writer) is still read via a legacy fallback path that
 * reproduces the old reader's exact behavior and messages, so no existing
 * file becomes unreadable.
 */

#ifdef MESHIOPLUSPLUS_HAS_HDF5

// System includes
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Write `mesh` as a CGNS/SIDS-compliant unstructured HDF5 file.
 *
 * Every cell block becomes its own `Elements_t` section (see the file
 * comment for the full tree); `CellDim`/`PhysDim` and the zone's vertex/cell
 * counts are derived from the mesh. Only the fixed-node-count types listed
 * in doc/formats/cgns.md's type table are supported — anything else
 * (`polygon`/`polyhedron*`, ragged blocks, or a type whose CGNS node
 * ordering meshio++ does not yet implement) throws by name rather than
 * writing a file with a guessed or silently dropped section.
 *
 * @param rPath filesystem path to write
 * @param rMesh the mesh to write — every cell block is emitted, in mesh
 *        order, as its own section (multiple blocks of the same type stay
 *        separate sections, unlike MED)
 * @param gzip_level HDF5 gzip compression level applied to every dataset;
 *        write-only, negative disables compression (the default)
 * @throws WriteError for an unsupported/unverified-ordering cell type, a
 *         ragged block, or `PointDim() > 3`
 */
MESHIOPLUSPLUS_API void write_cgns(const std::string& rPath, const Mesh& rMesh, int gzip_level);

/**
 * @brief Read a CGNS/HDF5 file: either one written by @ref write_cgns (or
 *        any other CGNS/SIDS-compliant unstructured-mesh writer covering
 *        the same type table), or a pre-v9.8.0 meshio++/legacy file.
 *
 * The spec path is selected by the presence of a root child node whose
 * `label` attribute is `"CGNSBase_t"` — the discriminator is structural,
 * never the file extension or a version number, and a file with neither
 * shape is a `ReadError` naming both. Every `Unstructured` zone under the
 * first `CGNSBase_t` is concatenated (a `Structured` zone is rejected rather
 * than silently dropped); `Elements_t` sections are found by `label`, never
 * by name, and ordered ascending by `ElementRange[0]` (which reproduces
 * `write_cgns`'s own block order exactly). The point array's column count
 * follows however many of `CoordinateX/Y/Z` the file actually has (at least
 * 2) rather than always `(n,3)`, so a 2-D-authored file round-trips its
 * point shape.
 *
 * @param rPath filesystem path to read
 * @return the read Mesh (points + one cell block per section, in
 *         `ElementRange`-sorted order; no point_data/cell_data/field_data)
 * @throws ReadError if the file has neither a spec `CGNSBase_t` node nor the
 *         legacy `Base/Zone1/GridElements` layout, a zone is `Structured`, a
 *         section has a `MIXED`/`NGON_n`/`NFACE_n`/unknown/unverified-
 *         ordering `ElementType`, or a section's connectivity size doesn't
 *         match its declared element count
 */
MESHIOPLUSPLUS_API Mesh read_cgns(const std::string& rPath);

}  // namespace meshioplusplus

#endif  // MESHIOPLUSPLUS_HAS_HDF5

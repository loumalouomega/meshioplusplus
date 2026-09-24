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
 * @file ansys.hpp
 * @brief Ansys/Fluent mesh (.msh) C++ reader/writer.
 *
 * Not to be confused with the unrelated Ansys MAPDL "coded database" format
 * handled by ansysinp.hpp. This is the Fluent `.msh` format: fully
 * parenthesis-nested "Scheme-like" sections `(<index> ...)`, where the index
 * may be a bare decimal (ASCII payload) or prefixed `20`/`30` for a binary
 * payload (`20xx` = float32 nodes / int32 cells, `30xx` = float64 / int64).
 * All connectivity and zone-header integers in both ASCII and binary bodies
 * are **hexadecimal** — the format's defining quirk. Section `10` gives node
 * blocks, `13` the faces (nodes plus the cells `c0 c1` on either side) and `12`
 * the cell zones. The reader rebuilds each cell from its faces, keeps boundary
 * faces as surface cells, and records each cell's zone in
 * `cell_data["ansys:zone"]` with one named cell region per zone; a legacy
 * meshio file (cell sections with connectivity bodies) is read as cells only.
 *
 * See doc/formats/ansys.md for the full section grammar and the
 * element-type/face-type code tables.
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Write `mesh` as a Fluent .msh file in Fluent's face-based layout.
 *
 * The cells are the blocks of the highest dimension (2 or 3), in cell zones
 * from `cell_data["ansys:zone"]` (else one per block). Every face is written
 * once with its nodes and the cells `c0 c1` on either side: 3-D faces from
 * `detail/face_mesh.hpp` with the right-hand normal into `c0`, 2-D edges with
 * `c0` on the left. Interior faces form one zone; a boundary face is in the
 * wall zone of the lower-dimensional block (or `ansys:zone`) whose facet
 * matches it, the rest in a default wall. Zone names (`45` sections) come from
 * the regions. Polygons and polyhedra are written as element type 7.
 *
 * @param rPath filesystem path to write
 * @param rMesh the mesh to write
 * @param binary write the node section as `3010` (float64) and the face
 *        sections as `2013` (int32), or everything ASCII (`false`)
 * @throws WriteError if the points are not 2-D or 3-D, the mesh has no 2-D or
 *         3-D cells, or a binary id does not fit 32 bits
 */
MESHIOPLUSPLUS_API void write_ansys(const std::string& rPath, const Mesh& rMesh, bool binary);

/**
 * @brief Read a Fluent .msh file.
 *
 * Parses `(0 ...)`/`(1 ...)`/`(2 ...)` header/comment/dimension sections
 * (bracket-skipped), `(<pfx>10 ...)` node sections (ascii one point per
 * line, binary a raw float32/float64 block), and `(<pfx>12 ...)` cell
 * sections (dead zones -> no cells; `mixed` zones structurally skipped, body
 * not decoded). All hexadecimal header/body integers are converted; the
 * result is one flat `Mesh.mCells` list with the first point-zone's `first`
 * index subtracted from every connectivity array.
 *
 * @param rPath filesystem path to read
 * @return the read Mesh (point_data/cell_data/field_data always empty)
 * @throws ReadError if the file can't be opened, a section header is
 *         malformed or truncated, a cell zone's `element-type` isn't one of
 *         the known volume codes (0/1/.../6), or a face (`13`) section
 *         carries a data body — **any** real face section always defers the
 *         whole file to the Python fallback, so files with real boundary
 *         face zones (a common real-world case) are never handled here;
 *         binary "mixed" faces additionally raise unconditionally even in
 *         the Python path
 * @note point_data/cell_data/field_data are never produced — this format
 *       carries no per-node/per-cell field values, only geometry and zone
 *       structure
 */
MESHIOPLUSPLUS_API Mesh read_ansys(const std::string& rPath);

}  // namespace meshioplusplus

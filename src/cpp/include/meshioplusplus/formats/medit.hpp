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
 * @file medit.hpp
 * @brief Medit / GMF (.mesh) ASCII C++ reader/writer.
 *
 * Medit (INRIA "libMeshb" GMF) is a keyword-section format. This header
 * covers only the **ASCII** `.mesh` variant: whitespace/`#`-comment
 * tokenized, `MeshVersionFormatted`/`Dimension` header, then `Vertices`
 * (`x1 x2 [x3] ref` rows, 1-based) and element keyword sections
 * (`Edges`/`Triangles`/`Quadrilaterals`/`Tetrahedra`/`Prisms`/`Pyramids`/
 * `Hexahedra`/`Hexaedra`), each a count followed by that many
 * `<node ids> ref` rows. The trailing per-row integer becomes
 * `point_data["medit:ref"]`/`cell_data["medit:ref"]`; Medit stores at most
 * one such integer column, so only the first int-dtype point/cell data
 * array is kept on write (extras are dropped with a warning).
 *
 * The binary `.meshb` GMF variant (position-indexed records, endianness
 * flip via a leading magic code, version-dependent int/float widths) is
 * `read_medit_binary` / `write_medit_binary`; `read_medit` and `write_medit`
 * dispatch on the `"b"` filename suffix.
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Write a Mesh to an ASCII Medit (.mesh) file.
 *
 * Emits `MeshVersionFormatted 2` (float64 coordinates) and `Dimension`,
 * then one keyword section per cell type present (`Vertices` is always
 * written; `Edges`/`Triangles`/`Quadrilaterals`/`Tetrahedra`/`Prisms`/
 * `Pyramids`/`Hexahedra` as applicable), each row 1-based node ids
 * followed by a trailing reference integer. At most one int-dtype
 * point_data array and one int-dtype cell_data array are used to populate
 * the `ref` columns (Medit's single-reference-column limitation); if
 * `point_data`/`cell_data` contain more than one integer candidate, the
 * first is used and the rest are dropped with a warning. Ends with `End`.
 *
 * @param path filesystem path to the .mesh file to create/overwrite
 * @param mesh the mesh to write
 * @throws WriteError on an unsupported cell type
 * @note reads `point_data`/`cell_data` key `"medit:ref"` if present
 *       (preferred over other int-dtype arrays for the ref column)
 */
MESHIOPLUSPLUS_API void write_medit_ascii(const std::string& rPath, const Mesh& rMesh);

/**
 * @brief Read an ASCII Medit (.mesh) file into a Mesh.
 *
 * Parses `MeshVersionFormatted` (0/1 -> float32 coords, 2 -> float64) and
 * `Dimension`, then every recognized element keyword section, converting
 * 1-based node ids to 0-based. Sections such as `Corners`, `Normals`,
 * `NormalAtVertices`, `SubDomainFromMesh`, `VertexOnGeometricVertex`/
 * `Edge`, `EdgeOnGeometricEdge`, `Identifier`, `Geometry`,
 * `RequiredVertices`, `TangentAtVertices`, `Tangents`, `Ridges` are
 * recognized only enough to be token-skipped, and are otherwise discarded.
 *
 * @param path filesystem path to the .mesh file to read
 * @return the read Mesh, with `point_data["medit:ref"]` and
 *         `cell_data["medit:ref"]` (one array per cell block) populated
 *         from each row's trailing reference integer
 * @throws ReadError on a malformed file or unrecognized required section
 */
MESHIOPLUSPLUS_API Mesh read_medit_ascii(const std::string& rPath);

/**
 * @brief Read a binary Medit (`.meshb`, libMeshb GMF) file.
 *
 * A leading code `1` (native byte order) or `16777216` (the opposite one), the
 * version 1-4 that fixes the integer, float and position widths, then the
 * keyword sections, each `[keyword][next position][count][records]`. Vertices,
 * the element sections and the `Corners` ... `RequiredTriangles` id lists are
 * read as the ASCII reader reads them (`medit:ref`, regions); every other
 * section is skipped through its next-section position (a warning names it).
 * @throws ReadError on a malformed or truncated file
 */
MESHIOPLUSPLUS_API Mesh read_medit_binary(const std::string& rPath);

/**
 * @brief Write a binary Medit (`.meshb`) file: little-endian, 64-bit float
 *        coordinates, version 3 (32-bit integers) or 4 when an index or
 *        reference needs 64 bits, one section per cell block and the id lists.
 *        The Python twin writes the same bytes.
 * @throws WriteError for points that are not 2-D or 3-D
 */
MESHIOPLUSPLUS_API void write_medit_binary(const std::string& rPath, const Mesh& rMesh);

/// `read_medit_binary` for a path ending in `b`, else `read_medit_ascii`.
MESHIOPLUSPLUS_API Mesh read_medit(const std::string& rPath);

/// `write_medit_binary` for a path ending in `b`, else `write_medit_ascii`.
MESHIOPLUSPLUS_API void write_medit(const std::string& rPath, const Mesh& rMesh);

}  // namespace meshioplusplus

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
 * `Hexahedra`/`Hexaedra`, and the high-order `EdgesP2`/`EdgesP3`/
 * `TrianglesP2`/`TrianglesP3`/`TetrahedraP2` as `line3`/`line4`/
 * `triangle6`/`triangle10`/`tetra10`), each a count followed by that many
 * `<node ids> ref` rows. libMeshb defines no node order for the high-order
 * sections: a file states its own in a `<Keyword>Ordering` section (one
 * integer tuple per node), which the readers apply and the writers emit
 * (doc/formats/medit.md). The trailing per-row integer becomes
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
 * `Pyramids`/`Hexahedra` and the high-order sections, each after its
 * `<Keyword>Ordering`, as applicable), each row 1-based node ids
 * followed by a trailing reference integer. At most one int-dtype
 * point_data array and one int-dtype cell_data array are used to populate
 * the `ref` columns (Medit's single-reference-column limitation); if
 * `point_data`/`cell_data` contain more than one integer candidate, the
 * first is used and the rest are dropped with a warning. Ends with `End`.
 *
 * A cell type Medit has no keyword for is skipped with a warning.
 *
 * @param path filesystem path to the .mesh file to create/overwrite
 * @param mesh the mesh to write
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
 * `TangentAtVertices`, `Tangents` are recognized only enough to be
 * token-skipped, and are otherwise discarded; `Corners`,
 * `RequiredVertices`, `Ridges`, `RequiredEdges` and `RequiredTriangles`
 * become regions. A high-order section is put in meshio's node order by
 * the file's `<Keyword>Ordering` (meshio's order without one); an ordering
 * that is not a permutation of meshio++'s tuples throws.
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

/**
 * @brief Read a Medit solution: ascii `.sol`, or binary `.solb` when the path
 *        ends in `b` (libMeshb's `SolAt*` keywords, MMG's size maps and metrics).
 *
 * The result has no geometry: `NumPoints()` rows of zero-width points carrying
 * the `SolAtVertices` fields as `point_data["medit:sol:<i>"]` (`(n,)` for a
 * scalar, `(n, components)` otherwise; a symmetric matrix keeps the file's
 * lower-triangle order) and the dimension as `field_data["medit:dimension"]`.
 * Element sections need their mesh: they are skipped with a warning, and
 * `read_medit_sol_into` attaches them.
 *
 * @throws ReadError on a malformed or truncated file.
 */
MESHIOPLUSPLUS_API Mesh read_medit_sol(const std::string& rPath);

/**
 * @brief Attach the fields of a Medit solution to `rMesh`: `SolAtVertices`
 *        fields as `point_data["medit:sol:<i>"]`, `SolAt<Elements>` fields as
 *        `cell_data["medit:sol:<i>"]` (NaN on the blocks of other types; a name
 *        two types give different shapes gets a `:<cell type>` suffix).
 *
 * @throws ReadError when the file is malformed, or its dimension or row counts
 *         disagree with the mesh.
 */
MESHIOPLUSPLUS_API void read_medit_sol_into(const std::string& rPath, Mesh& rMesh);

/**
 * @brief Write a mesh's floating-point point and cell data as a Medit solution:
 *        ascii `.sol`, or little-endian binary `.solb` when the path ends in `b`.
 *
 * Each array becomes one field of the type its component count names in the
 * mesh's dimension (1: scalar, `dim`: vector, `dim*(dim+1)/2`: symmetric
 * matrix, `dim*dim`: full matrix); integer arrays and other widths are skipped
 * with a warning. A geometry-less mesh takes its dimension from
 * `field_data["medit:dimension"]`.
 *
 * @throws WriteError when there is nothing to write or no 2-D/3-D dimension.
 */
MESHIOPLUSPLUS_API void write_medit_sol(const std::string& rPath, const Mesh& rMesh);

/// `read_medit_binary` for a path ending in `b`, else `read_medit_ascii`.
MESHIOPLUSPLUS_API Mesh read_medit(const std::string& rPath);

/// `write_medit_binary` for a path ending in `b`, else `write_medit_ascii`.
MESHIOPLUSPLUS_API void write_medit(const std::string& rPath, const Mesh& rMesh);

}  // namespace meshioplusplus

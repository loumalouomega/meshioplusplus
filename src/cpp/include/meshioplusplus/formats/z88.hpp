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
 * @file z88.hpp
 * @brief Z88 / Z88Aurora structure file (`z88i1.txt`) reader/writer, with the
 *        `z88o2.txt` displacements and `z88o3.txt` stresses.
 *
 * Z88's files have fixed names, so dispatch is by basename (`is_z88_filename`),
 * ahead of the `.txt` extension that belongs to the xyz reader.
 *
 * `z88i1.txt`: a header line whose first three integers are the dimension, the
 * node and the element count (Z88OS v15 writes `ndim nnodes nelem ndof kflag`;
 * Z88 <= V13 and Z88Aurora V1 add the material count and more flags, and put
 * material lines after the elements, which are skipped); one line per node,
 * `id ndof x y [z]`; two lines per element, `id type` then its node ids (the
 * count is fixed by the type). `KFLAG = 1` means cylindrical coordinates
 * (r, phi in degrees, z), converted to Cartesian.
 *
 * Element types: 1 hex8, 10 hex20, 16 tet10, 17 tet4, 7/8/20/23 quad8, 3/14/15/
 * 18/24 tri6, 6 tri3, 2/4/5/9/13/25 two-node bars and beams. The cubic and
 * layered types 11/12 (12-node quads), 19 (16-node plate), 21 (16-node shell)
 * and 22 (12-node shell) keep their corners (quad, quad, hexahedron, wedge)
 * with a warning. The Z88 type is kept as the `z88:type` cell data. Hexahedra
 * list their faces top first (`"z88"` tables of `detail/node_order.hpp`), and
 * tet10's last three mid-edge nodes run 2-4, 3-4, 1-4.
 *
 * Results next to the structure file are attached on read: `z88o2.txt` as
 * `point_data["U"]` (2, 3 or 6 columns, the widest node's; NaN where a node has
 * fewer), and the solid and plane-stress blocks of `z88o3.txt` as the
 * element's mean stress: `cell_data["SIG"]` (XX YY ZZ XY YZ ZX in 3-D, XX YY XY
 * in 2-D) and `cell_data["SIGV"]` when the file has the equivalent stress. Other
 * element families' stress blocks (beams, plates, shells, tori) are skipped.
 * See doc/formats/z88.md.
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Whether @p rPath's basename is one of Z88's fixed file names
 *        (`z88i1.txt`, `z88structure.txt`, `z88o2.txt`, `z88o3.txt`, any case).
 */
MESHIOPLUSPLUS_API bool is_z88_filename(const std::string& rPath);

/**
 * @brief Read a Z88 structure file, with its results when present.
 *
 * @param rPath the `z88i1.txt` (or `z88structure.txt`) to read; a `z88o2.txt`
 *        or `z88o3.txt` path reads the structure file next to it
 * @param Results whether to attach `z88o2.txt`/`z88o3.txt` from the same
 *        directory
 * @return the mesh
 * @throws ReadError if the file can't be read, is truncated, a field is
 *         malformed, an element names an undefined node or has an unknown type
 */
MESHIOPLUSPLUS_API Mesh read_z88(const std::string& rPath, bool Results = true);

/**
 * @brief Write `rMesh` as a Z88OS v15 structure file (`z88i1.txt`).
 *
 * The element type comes from the `z88:type` cell data when it is compatible
 * with the cell, else from the cell type (hexahedron 1, hexahedron20 10,
 * tetra 17, tetra10 16, triangle6 14, quad8 7, triangle 6 (2-D only), line 4 in
 * 3-D and 9 in 2-D). Other cells, regions and data arrays are dropped with a
 * warning and a provenance note. The dimension is 2 when every z coordinate is
 * zero and no cell is 3-D.
 *
 * @param rPath filesystem path to write
 * @param rMesh the mesh to write
 * @param Stubs also write empty `z88i2.txt` (no constraints) and `z88i5.txt`
 *        (no surface loads) next to it, so the deck is complete
 * @throws WriteError if nothing is left to write
 */
MESHIOPLUSPLUS_API void write_z88(const std::string& rPath, const Mesh& rMesh, bool Stubs = false);

}  // namespace meshioplusplus

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
 * @file mfem.hpp
 * @brief MFEM mesh (`.mesh`) and grid function (`.gf`) C++ reader/writer.
 *
 * MFEM's own ASCII mesh, the one GLVis opens (`MFEM mesh v1.0`, `v1.2`, `v1.3`):
 * `dimension`, `elements` (`attribute geometry vertex...`), `boundary` (the same
 * for faces or edges), `vertices`, and in v1.3 named `attribute_sets` and
 * `bdr_attribute_sets`. A curved mesh replaces the vertex coordinates by a
 * `nodes` grid function.
 *
 *  - Elements and boundary elements are separate cell blocks, one per type, the
 *    elements first. Their attribute is the `mfem:attribute` cell data and an
 *    `attribute_<n>` (elements) or `boundary_<n>` (boundary) cell region tagged
 *    `n`; an attribute set is a cell region of its own name.
 *  - An `H1` order-2 `nodes` space (`H1_<d>D_P2`, legacy `Quadratic`) gives
 *    `line3`/`triangle6`/`quad9`/`tetra10`/`wedge18`/`hexahedron27` cells. The
 *    degrees of freedom are numbered the way MFEM numbers them -- vertices, then
 *    edges and faces in order of first appearance, then element interiors -- and
 *    the points are those degrees of freedom, in that order.
 *  - Higher orders (and Bernstein or serendipity spaces) keep only the vertices,
 *    with a warning. `L2_T1_<d>D_P1` nodes (periodic meshes) give each element
 *    its own points.
 *  - Every geometry, the prism included, is in meshio++'s node order (MFEM's own
 *    VTK export reverses prisms for classic VTK's winding; meshio++ does not).
 *
 * Non-conforming (`MFEM NC mesh`, v1.1 `vertex_parents`), NURBS and INLINE
 * meshes are refused. A parallel rank file is read as its local part. See
 * doc/formats/mfem.md.
 */

// System includes
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// One grid function (`.gf`) to read onto an MFEM mesh, and the data name it gets.
struct MfemGridFunction {
    std::string mName;
    std::string mPath;
};

/**
 * @brief Read an MFEM `.mesh` file.
 * @param rPath filesystem path to read
 * @return the mesh
 * @throws ReadError if the file can't be read, is not a conforming MFEM mesh, a
 *         section is malformed or truncated, or the `nodes` space is unsupported
 */
MESHIOPLUSPLUS_API Mesh read_mfem(const std::string& rPath);

/**
 * @brief Read an MFEM `.mesh` file and grid functions defined on it.
 *
 * An `H1` order-1 or order-2 grid function becomes point data (an order-2 field
 * on a linear mesh makes the cells quadratic, with the new nodes at edge, face
 * and cell centres); an `L2` order-0 one becomes cell data (NaN on boundary
 * cells). Any other space is skipped with a warning.
 *
 * @param rPath filesystem path of the mesh
 * @param rGridFunctions the `.gf` files, each with the data name it gets
 * @return the mesh with the grid functions as data
 * @throws ReadError as `read_mfem(rPath)`, or if a grid function does not match
 *         the mesh
 */
MESHIOPLUSPLUS_API Mesh read_mfem(const std::string& rPath,
                                  const std::vector<MfemGridFunction>& rGridFunctions);

/**
 * @brief Write `rMesh` as an MFEM `.mesh` file.
 *
 * The cells of the highest dimension are the elements and those one dimension
 * lower the boundary; side regions add boundary elements. Attributes come from
 * `mfem:attribute`, else from cell regions, else 1; named regions become v1.3
 * attribute sets. Quadratic cells (the serendipity ones completed with their
 * face and cell centres) are written as an `H1_<d>D_P2` `nodes` space; a mesh
 * with pyramids is written linear. Data arrays are dropped with a warning (see
 * the overload).
 *
 * @param rPath filesystem path to write
 * @param rMesh the mesh to write
 * @throws WriteError for a mesh with no elements MFEM can hold
 */
MESHIOPLUSPLUS_API void write_mfem(const std::string& rPath, const Mesh& rMesh);

/**
 * @brief `write_mfem(rPath, rMesh)`, and, when `GridFunctions` is true, one
 * `<stem>.<name>.gf` next to it per data array: point data as `H1` at the mesh's
 * order, cell data as `L2` order 0 (element cells only).
 */
MESHIOPLUSPLUS_API void write_mfem(const std::string& rPath, const Mesh& rMesh, bool GridFunctions);

}  // namespace meshioplusplus

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
 *  - Higher orders give VTK Lagrange cells (`H1` and the legacy `Cubic`, since
 *    v16.11.0); Bernstein and serendipity spaces are evaluated at the Lagrange
 *    nodes. `L2_T1_<d>D_P1` nodes (periodic meshes) give each element its own
 *    points.
 *  - Every geometry, the prism included, is in meshio++'s node order (MFEM's own
 *    VTK export reverses prisms for classic VTK's winding; meshio++ does not).
 *
 * A non-conforming mesh (`MFEM NC mesh`) is read as its leaf elements, and a
 * NURBS mesh (`MFEM NURBS mesh` v1.0/v1.1) as its knot-span elements, each a
 * Lagrange cell of the highest knot-vector order through the rational patch
 * geometry; the overloads that take a `NurbsInfo` (since v16.40.0) also return
 * the patches themselves, and write them back. INLINE meshes and NURBS meshes
 * with non-conforming patches are refused. A parallel rank file is read with
 * its siblings. See doc/formats/mfem.md.
 */

// System includes
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/nurbs.hpp"
#include "meshioplusplus/read_options.hpp"

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
 * An `H1` grid function becomes point data (an order-2 field on a linear mesh
 * makes the cells quadratic, with the new nodes at edge, face and cell centres;
 * order 3 and up makes them VTK Lagrange cells, since v16.11.0); an `L2`
 * order-0 one becomes cell data (NaN on boundary cells). Any other space is
 * skipped with a warning.
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
 * @brief `read_mfem(rPath, rGridFunctions)`, honouring `ReadOptions::mPiece`
 * for a parallel mesh (since v16.11.0).
 *
 * A rank file of a parallel mesh (`<prefix>.000000`, ... as `ParMesh::Save`
 * writes them) is read with its sibling ranks as one mesh: the vertices the
 * communication groups share are merged, and `cell_data["partition:part"]`
 * holds each cell's rank. With `mPieceSet`, rank `mPiece` alone is read (still
 * labelled). A grid function of a parallel mesh is named by one of its rank
 * files, `<name>.000000`, ...; each rank reads its own.
 *
 * @param rPath filesystem path of the mesh (or of one rank file)
 * @param rGridFunctions the `.gf` files, each with the data name it gets
 * @param rOptions `mPiece`/`mPieceSet` select one rank of a parallel mesh
 * @return the mesh with the grid functions as data
 * @throws ReadError as the other overloads, or for a piece out of range
 */
MESHIOPLUSPLUS_API Mesh read_mfem(const std::string& rPath,
                                  const std::vector<MfemGridFunction>& rGridFunctions,
                                  const ReadOptions& rOptions);

/**
 * @brief `read_mfem(rPath, rGridFunctions, rOptions)`, and the NURBS patches of a
 * NURBS mesh in @p rInfo (since v16.40.0).
 *
 * For an `MFEM NURBS mesh`, @p rInfo receives what the sampling into Lagrange
 * cells drops: per patch its degrees, knot vectors (oriented along the patch's
 * own directions), control net and weights, and the control values of each
 * NURBS grid function read with it (`mFieldNames` are their data names); the
 * file's patch topology; and the patch and knot span of each element cell
 * (the first block) and the boundary patch of each boundary cell. Any other
 * mesh leaves @p rInfo empty. @p rInfo is replaced on success and unchanged if
 * the read fails.
 *
 * @param rPath filesystem path of the mesh
 * @param rGridFunctions the `.gf` files, each with the data name it gets
 * @param rOptions `mPiece`/`mPieceSet` select one rank of a parallel mesh
 * @param rInfo receives the NURBS side channel
 * @return the mesh with the grid functions as data
 * @throws ReadError as the other overloads
 */
MESHIOPLUSPLUS_API Mesh read_mfem(const std::string& rPath,
                                  const std::vector<MfemGridFunction>& rGridFunctions,
                                  const ReadOptions& rOptions, NurbsInfo& rInfo);

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

/**
 * @brief Write the NURBS patches of @p rInfo as an `MFEM NURBS mesh v1.0`
 * (since v16.40.0); an empty @p rInfo writes as `write_mfem(rPath, rMesh,
 * GridFunctions)`.
 *
 * The file is MFEM's global form (`knotvectors`, `weights` and a `NURBS<p>`
 * nodes space with Cartesian control points), so `read_mfem` of it gives back
 * @p rMesh and @p rInfo exactly; the patches, their attributes and the
 * boundary come from @p rInfo. @p rMesh must still be the sampling of the
 * patches -- its points, its element cells and the point data of the NURBS
 * fields equal what reading the file would give, within 1e-9 of their
 * magnitude -- so a mesh moved or edited after the read is refused rather than
 * written with stale patches. With `GridFunctions`, each field of @p rInfo
 * whose point data the mesh still has is written as a NURBS grid function
 * `<stem>.<name>.gf`, element cell data as `L2` order 0; other point data,
 * field data and regions other than the attributes are dropped with a warning.
 * Spacing formulas of a v1.1 file are not kept (its knots are).
 *
 * @throws WriteError for an inconsistent @p rInfo, or a mesh that no longer
 *         matches it ("the mesh no longer matches its NURBS side channel")
 */
MESHIOPLUSPLUS_API void write_mfem(const std::string& rPath, const Mesh& rMesh, bool GridFunctions,
                                   const NurbsInfo& rInfo);

}  // namespace meshioplusplus

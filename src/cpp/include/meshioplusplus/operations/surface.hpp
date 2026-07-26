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
 * @file operations/surface.hpp
 * @brief Surface/boundary extraction — the general form of boundary-skin
 * extraction, part of the dependency-free mesh-operations layer.
 *
 * A *facet* (a face of a 3D cell, or an edge of a 2D cell) is on the boundary
 * when its sorted corner-node key occurs exactly once across all cells of the
 * chosen dimension — facets shared by two cells are interior and cancel out
 * (the face-hashing algorithm of Kratos Multiphysics' `SkinDetectionProcess`).
 *
 * `extract_surface` chooses the dimension automatically: if the mesh has any
 * 3D volume cells it returns their boundary faces (`triangle`/`quad` and the
 * quadratic `triangle6`/`quad8`/`quad9`); otherwise, for a surface mesh, it
 * returns the boundary edges (`line`/`line3`). Per-cell-type facet topology
 * comes from `detail/cell_faces.hpp` (3D) and `detail/cell_edges.hpp` (2D).
 *
 * `extract_skin` (declared in `skin.hpp`) is the volume-only, `linearize`-able
 * special case; both share the implementation in `operations/surface.cpp`.
 */

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Extract the boundary of a mesh's highest-dimension cells as a new mesh.
 *
 * Picks the dimension automatically: 3D volume cells present → boundary faces
 * (`triangle`, `triangle6`, `quad`, `quad8`, `quad9`); otherwise 2D surface
 * cells → boundary edges (`line`, `line3`). Facets whose sorted corner-node
 * key occurs exactly once are kept, in cell-block enumeration order,
 * partitioned into the fixed canonical block order for the chosen dimension.
 * Points are *compacted* (only referenced nodes, in ascending original-index
 * order) and `point_data` rows are subset the same way; `cell_data`,
 * `field_data`, sets, and `info` are dropped.
 *
 * Unsupported same-dimension blocks (ragged/polyhedron blocks, Lagrange and
 * other very-high-order types) are skipped with a warning; lower-dimension
 * blocks are ignored.
 *
 * @param rMesh the mesh to extract the boundary of
 * @param recordParentIds when true, attaches an `Int64` cell-data array
 *        `"surface:parent_cell"` giving, for each output facet, the global
 *        input-cell index (counted block-major over every input cell of every
 *        block) of the unique cell that owns it
 * @return a new mesh holding the boundary facets
 * @throws std::invalid_argument if `rMesh` has no supported 2D/3D cell block
 */
MESHIOPLUSPLUS_API Mesh extract_surface(const Mesh& rMesh, bool recordParentIds = false);

/**
 * @brief Whether a mesh has at least one cell block `extract_surface` supports.
 * @param rMesh the mesh to test
 * @return true when `extract_surface(rMesh)` would have input to work on
 */
MESHIOPLUSPLUS_API bool has_surface_extractable_cells(const Mesh& rMesh);

/**
 * @brief Copy a mesh's cell data onto the boundary facets extracted from it.
 *
 * `extract_surface` drops cell data -- a facet is not a cell, so there is no
 * general answer to what its value should be. With `recordParentIds` it does
 * record which input cell owned each facet, and for a *viewer* the useful
 * answer is exactly that owner's value: colouring a solid by its per-cell
 * material or tag is the common case, and without this the array simply
 * vanishes between the mesh and the renderer.
 *
 * `rSurface` must have been produced by `extract_surface(rSource, true)`; the
 * `"surface:parent_cell"` array is what this reads. If it is absent this is a
 * no-op. Arrays whose block layout does not match `rSource` are skipped rather
 * than guessed at.
 *
 * @param rSource the mesh the boundary was extracted from
 * @param rSurface the extracted boundary, modified in place
 */
MESHIOPLUSPLUS_API void gather_cell_data_onto_surface(const Mesh& rSource, Mesh& rSurface);

namespace detail {

/**
 * @brief The shared core of `extract_surface` and `extract_skin`, with the
 *        `linearize` and `recordParentIds` flags exposed orthogonally.
 *
 * The two public entry points each pin one of the flags -- `extract_surface`
 * always keeps the input order (`linearize=false`) and `extract_skin` never
 * records provenance -- so neither can produce a *linearized* boundary that
 * still says which cell each facet came from. The SVG/TikZ writers need
 * exactly that combination to colour a projected volume mesh by cell data.
 *
 * This is deliberately not part of the stable API: `extract_surface` and
 * `extract_skin` are bound flat into the C, Fortran and WASM surfaces, and no
 * binding consumer needs the combination, so widening those signatures would
 * cost four ABIs to serve one internal caller.
 *
 * @param rMesh the mesh to extract the boundary of
 * @param forceFaceMode true selects volume-only face mode (`extract_skin`);
 *        false auto-selects faces or edges by the mesh's max dimension
 * @param linearize when true, emit linear facets for higher-order input
 * @param recordParentIds when true, attach `"surface:parent_cell"` (see
 *        `extract_surface`)
 * @param pOpName the caller's name, used in error and warning messages
 * @return a new mesh holding the boundary facets
 */
MESHIOPLUSPLUS_API Mesh surface_extract(const Mesh& rMesh, bool forceFaceMode, bool linearize, bool recordParentIds,
                     const char* pOpName);

}  // namespace detail

}  // namespace meshioplusplus

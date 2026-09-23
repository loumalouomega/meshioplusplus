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
 * @file elmer.hpp
 * @brief Elmer mesh directory (ElmerSolver's native mesh) C++ reader/writer.
 *
 * An Elmer mesh is a *directory*, not a file: `mesh.header` (counts),
 * `mesh.nodes` (`id part x y z`), `mesh.elements` (`id body type nodes…`),
 * `mesh.boundary` (`id boundary parent1 parent2 type nodes…`) and, optionally,
 * `mesh.names` (`$ name = id` lines naming bodies and boundaries). Element type
 * codes are family × 100 + node count (`101`, `202`/`203`/`204`, `303`/`306`/
 * `310`, `404`/`408`/`409`, `504`/`510`, `605`/`613`, `706`/`715`/`718`,
 * `808`/`820`/`827`); `820` and `827` go through the `"elmer"` tables of
 * `detail/node_order.hpp`, every other code is in meshio++'s own order.
 *
 * Reading:
 *  - bulk elements and boundary elements become separate cell blocks, one per
 *    type code; body ids become `Cell` regions (tag = the id) named from
 *    `mesh.names` or `body_<id>`, boundary ids likewise (`boundary_<id>`). A name
 *    used for both a body and a boundary is prefixed `body:`/`boundary:`.
 *    Boundary elements stay cells (not `Side` regions): Elmer allows boundary
 *    elements with no parent (point conditions) and between two bodies.
 *  - a partitioned mesh (`partitioning.N/part.n.*`, as ElmerGrid writes it) is
 *    merged into one mesh, shared nodes and halo elements once each, with the
 *    0-based part of every cell in `cell_data["partition:part"]`. The path may be
 *    the `partitioning.N` directory itself, or a mesh directory holding exactly
 *    one; a serial mesh next to exactly one `partitioning.N` is read serially
 *    and labelled from it. `ReadOptions::mPiece` reads one part alone.
 *  - `rPath` may also name the `mesh.header` file, standing for its directory.
 *  - ElmerGrid's binary `mesh.*.bin` variant is not read.
 *
 * Writing (always serial, ASCII): the cells of the highest dimension are the
 * bulk elements, every lower-dimensional cell a boundary element, and every
 * `Side` region's facets further boundary elements; parents are regenerated from
 * the mesh. Body and boundary ids come from the `Cell`/`Side` regions (their tag
 * when positive and unused, otherwise the next free id), cells in no region get
 * a fresh id per block. See doc/formats/elmer.md.
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/read_options.hpp"

namespace meshioplusplus {

/**
 * @brief Read an Elmer mesh directory.
 * @param rPath the mesh directory, a `partitioning.N` directory, or `mesh.header`
 * @param rOptions `mPiece`/`mPieceSet` select one part of a partitioned mesh
 * @return the mesh, bodies and boundaries as `Cell` regions
 * @throws ReadError if the directory holds no Elmer mesh, a file is malformed, a
 *         type code has no meshio++ cell type, or an element names an undefined
 *         node
 */
MESHIOPLUSPLUS_API Mesh read_elmer(const std::string& rPath, const ReadOptions& rOptions = {});

/**
 * @brief Write `rMesh` as a serial Elmer mesh directory, created if needed.
 *
 * Writes `mesh.header`, `mesh.nodes`, `mesh.elements`, `mesh.boundary` and
 * `mesh.names` (the last always, carrying the provenance block and the names of
 * the regions written as bodies and boundaries). Point regions and data arrays
 * have no Elmer equivalent and are dropped with a warning; an existing
 * `partitioning.*` directory is left in place, with a warning.
 *
 * @param rPath the directory to write
 * @param rMesh the mesh to write
 * @throws WriteError for a cell type with no Elmer code, or a mesh with no cells
 */
MESHIOPLUSPLUS_API void write_elmer(const std::string& rPath, const Mesh& rMesh);

}  // namespace meshioplusplus

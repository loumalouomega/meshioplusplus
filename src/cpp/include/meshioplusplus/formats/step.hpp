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
 * @file formats/step.hpp
 * @brief STEP (ISO 10303-21, AP214) writer: a closed surface as a faceted
 * B-rep solid a CAD tool opens -- the CAD hand-off of the skin-to-CAD work.
 *
 * The writer is write-only: reading STEP means a CAD kernel (roadmap §6.2.1; see doc/spikes.md),
 * which the core deliberately does not carry.
 *
 * Each surface facet becomes an `ADVANCED_FACE` on a `PLANE`, bounded by an
 * `EDGE_LOOP` of `ORIENTED_EDGE`s over `EDGE_CURVE`s on `LINE`s; an edge two
 * facets share is one `EDGE_CURVE`, so the shell is topologically closed with
 * no sewing. Every edge-connected piece of the surface is one
 * `CLOSED_SHELL` inside one `MANIFOLD_SOLID_BREP`, collected in an
 * `ADVANCED_BREP_SHAPE_REPRESENTATION` under the usual AP214 product
 * structure, in millimetres.
 *
 *  - **Input**: `triangle`, `quad` and `polygon` blocks (quadratic surface
 *    cells contribute their corners), or -- when the mesh has volume cells --
 *    the skin of those (`extract_skin`), as the `stl` writer takes it. Lines
 *    and vertices are dropped with a warning.
 *  - **Planarity**: a facet of four or more corners further than `1e-9` of
 *    its size from its own plane is fanned into triangles from its first
 *    corner, so every face lies exactly on its `PLANE`.
 *  - **Closed and oriented**: an open, non-manifold or inconsistently wound
 *    edge throws `WriteError`, naming the first few as point-id pairs (run
 *    `repair` first). A piece wound inwards (negative signed volume) is
 *    written reversed, with a warning, so every solid has positive volume.
 *  - **Bytes**: coordinates are written with `%.17g` (with the decimal point
 *    STEP requires) and the header time stamp is fixed, so a file is
 *    reproducible and the numpy twin writes the same bytes.
 *
 * One face per facet: merging coplanar facets into one face is not done.
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Write a closed surface (or a volume mesh's skin) as a STEP faceted solid.
 * @param rPath the `.step` / `.stp` file to create.
 * @param rMesh the mesh to write.
 * @throws WriteError on an open, non-manifold or inconsistently wound surface,
 *         a degenerate facet, a non-finite coordinate, or no facets at all.
 */
MESHIOPLUSPLUS_API void write_step(const std::string& rPath, const Mesh& rMesh);

}  // namespace meshioplusplus

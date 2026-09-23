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
 * @file febio.hpp
 * @brief FEBio input file (`.feb`) C++ reader/writer: the mesh, not the model.
 *
 * A `.feb` file is XML with a `febio_spec` root whose `version` attribute says
 * where the mesh lives: `<Geometry>` in spec 2.5, `<Mesh>` (plus
 * `<MeshDomains>`) in 3.0 and 4.0. All three are read; 4.0 is written.
 *
 *  - `<Nodes>` give the points (ids may be sparse); each `<Elements type=
 *    name=>` block one cell block and a `Cell` region of its name. hex27 goes
 *    through the `"febio"` table of `detail/node_order.hpp`; every other type is
 *    in meshio++'s (VTK's) order.
 *  - `<NodeSet>` is a `Point` region, `<ElementSet>` a `Cell` region, and a
 *    `<Surface>` a `Side` region when every facet is a face of a solid element
 *    (otherwise its facets become a cell block and a `Cell` region). `<Edge>`
 *    and `<DiscreteSet>` become line blocks with a `Cell` region.
 *  - `<MeshData>` `<NodeData>`/`<ElementData>` become point/cell data, NaN
 *    outside the set they are defined on.
 *  - Materials, loads, boundary conditions, steps and every other section are
 *    not read. The spec-2.5/3.0 `<Part>`/`<Instance>` form is refused.
 *
 * See doc/formats/febio.md.
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/read_options.hpp"

namespace meshioplusplus {

/**
 * @brief Read the mesh of an FEBio `.feb` file (spec 2.5, 3.0 or 4.0).
 * @param rPath filesystem path to read
 * @param rOptions `mLenient` downgrades element types with no meshio++ cell
 *        type (tet5, tet15) to their linear or quadratic form with a warning
 * @return the mesh, with sets, surfaces and domains as regions
 * @throws ReadError on malformed XML, an unsupported spec version, the `<Part>`
 *         form, an unknown element type, or a reference to an undefined node,
 *         element or set
 */
MESHIOPLUSPLUS_API Mesh read_febio(const std::string& rPath, const ReadOptions& rOptions = {});

/**
 * @brief Write `rMesh` as an FEBio spec-4.0 `.feb` mesh.
 *
 * Writes `<Module>`, a placeholder `<Material>` per domain (FEBio refuses a
 * domain whose material is undefined), `<Mesh>` and `<MeshDomains>`. A 2-D block
 * whose cells are all faces of solids is written as a `<Surface>`, and a line
 * block whose cells are all edges of other cells as an `<Edge>`; every other
 * block is an `<Elements>` block with a solid, shell or beam domain. Point
 * regions become `<NodeSet>`s, cell regions `<ElementSet>`s (or name the block
 * they cover exactly), side regions `<Surface>`s. Vertex cells and data arrays
 * are dropped with a warning.
 *
 * @param rPath filesystem path to write
 * @param rMesh the mesh to write
 * @throws WriteError for a cell type with no FEBio element type
 */
MESHIOPLUSPLUS_API void write_febio(const std::string& rPath, const Mesh& rMesh);

}  // namespace meshioplusplus

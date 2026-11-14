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
 * @file formats/threemf.hpp
 * @brief 3MF, the 3D Manufacturing Format's core specification (roadmap 2.4.5):
 * what PrusaSlicer, Cura and Bambu Studio exchange.
 *
 * The reader opens the ZIP package (`detail/zip_archive.hpp`), follows
 * `_rels/.rels` to the model part (`3D/3dmodel.model` by default) and instantiates
 * every `<build>` `<item>`: an object's `<mesh>` (`<vertices>`, `<triangles>`) or its
 * `<components>`, each with its `transform` (the 3x4 row-vector matrix the
 * specification defines) composed down the tree. All triangles form one block,
 * each build item a cell region named after its object (`object_<id>` without a
 * name). The model's `unit` becomes `meshio:length_unit`. Base materials
 * (`<basematerials>`, the triangle's `pid`/`p1` or its object's `pid`/`pindex`)
 * become the `uint8` RGBA cell data `3mf:color`. A model that requires an extension
 * (`requiredextensions`) is refused by name: the materials, production, slice and
 * beam-lattice extensions are out of scope.
 *
 * The writer emits `[Content_Types].xml`, `_rels/.rels` and `3D/3dmodel.model`: one
 * object per cell region (the triangles in no region form one more object, `mesh`),
 * each with the vertices it uses, a build item each, the unit from
 * `meshio:length_unit` when it is one of 3MF's (millimetre otherwise) and base
 * materials from `3mf:color`. Triangles only. The Python twin is
 * `meshioplusplus.threemf`, which writes the same bytes.
 * @see doc/formats/3mf.md
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/** @brief Read a 3MF package's build. */
MESHIOPLUSPLUS_API Mesh read_3mf(const std::string& rPath);

/** @brief Write a triangle mesh as a 3MF package. */
MESHIOPLUSPLUS_API void write_3mf(const std::string& rPath, const Mesh& rMesh);

}  // namespace meshioplusplus

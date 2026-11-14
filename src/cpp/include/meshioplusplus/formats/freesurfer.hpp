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
 * @file formats/freesurfer.hpp
 * @brief FreeSurfer binary surfaces, with `curv` overlays and `annot` parcellations
 * (roadmap 2.4.4).
 *
 * A triangle surface (`lh.pial`, `rh.white`, ...): the magic `0xFFFFFE`, a
 * "created by" line ended by two newlines, big-endian `int32` vertex and face
 * counts, `float32` coordinates and `int32` faces; whatever follows (the volume
 * geometry tags FreeSurfer appends) is kept verbatim in
 * `field_data["freesurfer:tail"]` and the creator line in
 * `field_data["freesurfer:created_by"]` (bytes), and both are written back. Quad
 * surfaces (`0xFFFFFF`/`0xFFFFFD` magic) are refused by name.
 *
 * Overlays (`FreesurferReadOptions`): each new-format `curv` file (`0xFFFFFF`
 * magic, vertex and face counts, values per vertex, `float32` values) becomes point
 * data named after the file with its `lh.`/`rh.` prefix dropped (`lh.thickness` is
 * `thickness`); an `annot` file becomes the Int32 point data `freesurfer:annot`
 * (each vertex's structure index, -1 where its colour matches no entry) and one
 * point region per colour-table entry (tag = the index), from either colour-table
 * layout.
 *
 * The writer emits a triangle surface from one triangle block. The Python twin is
 * `meshioplusplus.freesurfer`; nibabel reads what it writes.
 * @see doc/formats/freesurfer.md
 */

// System includes
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/** @brief The overlay files to read with a FreeSurfer surface. */
struct FreesurferReadOptions {
    std::vector<std::string> mCurv;  ///< curv files, each one point data array
    std::string mAnnot;              ///< an annot file ("" = none)
};

/** @brief Read a FreeSurfer triangle surface and its overlays. */
MESHIOPLUSPLUS_API Mesh read_freesurfer(const std::string& rPath,
                                        const FreesurferReadOptions& rOptions);

/** @brief Read a FreeSurfer triangle surface. */
MESHIOPLUSPLUS_API Mesh read_freesurfer(const std::string& rPath);

/** @brief Write one triangle block as a FreeSurfer triangle surface. */
MESHIOPLUSPLUS_API void write_freesurfer(const std::string& rPath, const Mesh& rMesh);

}  // namespace meshioplusplus

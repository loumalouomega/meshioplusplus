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
 * @file formats/gifti.hpp
 * @brief GIFTI (`.gii`) surfaces: neuroimaging meshes and their per-vertex data
 * (roadmap 2.4.4).
 *
 * The reader takes every `DataArray`: `ASCII`, `Base64Binary` and
 * `GZipBase64Binary` encodings (the zlib stream inflated with window 15), either
 * `Endian`, `RowMajorOrder` or `ColumnMajorOrder`, and the NIfTI integer and float
 * types. `NIFTI_INTENT_POINTSET` is the points (its first `CoordinateSystemTransformMatrix`
 * in `field_data["gifti:transform"]`, 4x4, and the space names as bytes in
 * `gifti:data_space` and `gifti:transformed_space`), `NIFTI_INTENT_TRIANGLE` the
 * triangle block, `NIFTI_INTENT_LABEL` an integer point data array whose
 * `LabelTable` entries become point regions (tag = the label key), and every other
 * array with one row per point is point data named by its `Name` metadata (else
 * `gifti:<intent>_<index>`). External-file encodings are refused by name.
 *
 * The writer emits one triangle block as `POINTSET` (float32) and `TRIANGLE`
 * (int32), then each point data array (float32, or int32 for an integer array;
 * `NIFTI_INTENT_LABEL` for the one named `gifti:label` or any integer array a
 * tagged point region labels), little-endian `Base64Binary`, row-major, with a
 * `LabelTable` from the tagged point regions. The Python twin is
 * `meshioplusplus.gifti`; nibabel reads what it writes.
 * @see doc/formats/gifti.md
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/** @brief Read a GIFTI surface and its per-vertex arrays. */
MESHIOPLUSPLUS_API Mesh read_gifti(const std::string& rPath);

/** @brief Write a triangle surface and its point data as GIFTI. */
MESHIOPLUSPLUS_API void write_gifti(const std::string& rPath, const Mesh& rMesh);

}  // namespace meshioplusplus

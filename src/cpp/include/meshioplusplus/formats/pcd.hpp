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
 * @file pcd.hpp
 * @brief Point Cloud Library PCD (v0.7) reader/writer.
 *
 * A PCD file is a text header (`VERSION`, `FIELDS`, `SIZE`, `TYPE`, `COUNT`, `WIDTH`,
 * `HEIGHT`, `VIEWPOINT`, `POINTS`, `DATA`) followed by the points as ASCII rows
 * (`DATA ascii`), little-endian packed records (`binary`) or a struct-of-arrays block
 * inside an LZF stream (`binary_compressed`, two `uint32` sizes then the stream).
 *
 * **Mesh mapping.** The points plus exactly one `vertex` block (what
 * `subsample_points` emits). Point data is keyed by field name, with three
 * conventions: `normal_x/y/z` -> `"normals"` (n, 3); PCL's `rgb`/`rgba`, a `uint32`
 * `0x00RRGGBB`/`0xAARRGGBB` that lives in a `float32` slot (unpacked by bit-cast, never
 * by value) -> `"rgb"` (n, 3) / `"rgba"` (n, 4) `uint8`; any other field keeps its own
 * name and dtype, with `COUNT > 1` giving (n, count). `_` padding fields are skipped.
 * `x`/`y`/`z` keep the file's precision (all `F4` -> float32, otherwise float64).
 *
 * **Organised clouds** (`HEIGHT > 1`) are kept whole, NaN rows included, with
 * `WIDTH`/`HEIGHT` in `field_data["pcd:width"]`/`["pcd:height"]` (restored on write
 * when `WIDTH * HEIGHT` equals the point count). A `VIEWPOINT` that is not the identity
 * `0 0 0 1 0 0 0` goes to `field_data["pcd:viewpoint"]`; it is recorded, never applied.
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/** @brief The `DATA` encoding written by `write_pcd`. */
enum class PcdData {
    Ascii,             ///< `DATA ascii`
    Binary,            ///< `DATA binary` (the default)
    BinaryCompressed,  ///< `DATA binary_compressed` (LZF over struct-of-arrays)
};

/** @brief Options for `read_pcd`. */
struct PcdReadOptions {
    /// Drop points whose x, y or z is not finite (an organised cloud's invalid
    /// returns) and, when any is dropped, the organisation with them.
    bool mDropInvalid = false;
};

/**
 * @brief Read a PCD file.
 * @param rPath filesystem path
 * @param rOptions read options
 * @return points, one `vertex` block, the point data described above
 * @throws ReadError on a malformed header, an unsupported field type, or truncated or
 *         inconsistent data
 */
MESHIOPLUSPLUS_API Mesh read_pcd(const std::string& rPath,
                                 const PcdReadOptions& rOptions = PcdReadOptions());

/**
 * @brief Write a PCD file.
 *
 * Writes the points and the point data; cells other than `vertex`, cell data and
 * unrelated field data are dropped with a warning and a provenance note. Fields are
 * `x y z`, then the point data in name order. `normals`, `curvature` and `intensity` are
 * stored in the points' precision.
 *
 * @param rPath filesystem path
 * @param rMesh the mesh
 * @param data the `DATA` encoding
 * @param float64_points write `x y z` (and normals/curvature/intensity) as float64
 *        instead of PCL's float32 (a lossy narrowing is recorded in the provenance)
 * @throws WriteError if the file cannot be created
 */
MESHIOPLUSPLUS_API void write_pcd(const std::string& rPath, const Mesh& rMesh,
                                  PcdData data = PcdData::Binary, bool float64_points = false);

}  // namespace meshioplusplus

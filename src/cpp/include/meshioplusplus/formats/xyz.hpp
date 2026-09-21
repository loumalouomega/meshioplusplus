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
 * @file xyz.hpp
 * @brief Headerless ASCII point clouds (`.xyz`, `.xyzn`, `.xyzrgb`, `.asc`, `.pts`,
 *        `.txt`).
 *
 * One point per line, columns separated by whitespace, commas or semicolons; blank
 * lines and `#` / `//` comments are skipped. XYZ is a convention, not a specification,
 * so the columns are resolved in this order: an explicit column list; a header comment
 * naming them (`# x y z nx ny nz`, what the writer emits); the column count and the
 * file extension (3 = xyz, 4 = xyz + scalar, 6 = by extension or by value range, `.pts`
 * with 7 = xyz + intensity + rgb); anything still ambiguous is an error asking for the
 * list. Chemistry XYZ (atom count, comment line, `element x y z`) shares the extension
 * and is refused by name. `.pts` files start with a point count that is validated.
 *
 * **Mesh mapping.** float64 points plus one `vertex` block. Columns `nx ny nz` ->
 * `"normals"` (n, 3); `r g b [a]` -> `"rgb"` / `"rgba"` `uint8` (byte values, or unit
 * floats scaled by 255); every other named column is a float64 scalar under its own name.
 */

// System includes
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/** @brief Options for `read_xyz`. */
struct XyzReadOptions {
    /// Column names, one per column: `x y z nx ny nz r g b a`, any other name is a
    /// scalar, `_` skips the column. Empty = infer.
    std::vector<std::string> mColumns;
    /// Column delimiter. Empty = detect (`;`, then `,`, else whitespace).
    std::string mDelimiter;
};

/**
 * @brief Read a headerless ASCII point cloud.
 * @throws ReadError on ragged or non-numeric rows, an ambiguous column layout, a
 *         chemistry XYZ file, or a `.pts` count mismatch
 */
MESHIOPLUSPLUS_API Mesh read_xyz(const std::string& rPath,
                                 const XyzReadOptions& rOptions = XyzReadOptions());

/**
 * @brief Write a point cloud as `# x y z ...` plus one row per point.
 *
 * Cells other than `vertex` and all cell data are dropped with a warning. `normals`
 * becomes `nx ny nz`, `rgb`/`rgba` become `r g b [a]`, other arrays `name` or
 * `name_0 .. name_k`.
 *
 * @param rFloatFormat a `printf` floating conversion without the `%` (e.g. `.16e`);
 *        empty = `.9g` for float32 columns and `.17g` for the rest
 */
MESHIOPLUSPLUS_API void write_xyz(const std::string& rPath, const Mesh& rMesh,
                                  const std::string& rFloatFormat = "");

}  // namespace meshioplusplus

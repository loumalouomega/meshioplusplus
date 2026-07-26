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
 * @file obj_off.hpp
 * @brief Wavefront OBJ (.obj) and Geomview OFF (.off) ASCII C++
 *        readers/writers — two distinct surface formats sharing one
 *        header.
 *
 * **OBJ**: a line-oriented format with `v x y z` (points), `vn`/`vt`
 * (vertex normals / texture coordinates, stored raw with arbitrary column
 * count), `s` (smooth-shading toggle, ignored), `f i1[/t1[/n1]] ...` (faces
 * — only the leading vertex index of each `i/t/n` token is used; a run of
 * faces stays in one cell block until the per-face vertex count changes or
 * a new `g` line appears), and `g <name>` (starts a new group, incrementing
 * a running group-id counter starting at -1; only the numeric id is kept,
 * not the name string). Faces are grouped by vertex count into `triangle`
 * (3), `quad` (4), or `polygon` (else). Blank trailing groups are dropped.
 *
 * **OFF**: a minimal format — a literal `"OFF"` first line, a
 * `<nverts> <nfaces> <nedges>` header (edge count parsed but discarded),
 * `nverts` coordinate rows, then `nfaces` rows each `n i0 i1 ... i(n-1)`.
 * As with OBJ, faces are grouped by vertex count into `triangle` (3),
 * `quad` (4), or `polygon` (else) cell blocks, with a run of same-count
 * faces staying in one block until the count changes; a leading count
 * below 3 is a hard `ReadError`. OFF carries no point_data, cell_data, or
 * field_data at all.
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Write a Mesh to a Geomview OFF (.off) file.
 *
 * Emits the `"OFF"` header line, `<nverts> <nfaces> 0` (edge count always
 * 0), vertex coordinate rows, then one `n i0 ... i(n-1)` row per face,
 * for every `triangle`/`quad`/`polygon` cell block in mesh order (a
 * `polygon` block must be rectangular — same vertex count for every cell
 * in the block). Any other cell type is skipped with a warning.
 *
 * @param rPath filesystem path to the .off file to create/overwrite
 * @param rMesh the mesh to write
 */
MESHIOPLUSPLUS_API void write_off(const std::string& rPath, const Mesh& rMesh);

/**
 * @brief Read a Geomview OFF (.off) file into a Mesh.
 *
 * Validates the `"OFF"` first line, reads the vertex/face/edge counts
 * (edge count discarded), then `nverts` coordinate rows and `nfaces` face
 * rows. Faces are grouped by vertex count into `triangle`/`quad`/`polygon`
 * cell blocks, exactly like the OBJ reader above.
 *
 * @param rPath filesystem path to the .off file to read
 * @return the read Mesh (points + triangle/quad/polygon cell blocks only —
 *         no point_data/cell_data/field_data)
 * @throws ReadError if the first line isn't `"OFF"`, or any face row's
 *         leading vertex count is below 3
 */
MESHIOPLUSPLUS_API Mesh read_off(const std::string& rPath);

/**
 * @brief Write a Mesh to a Wavefront OBJ (.obj) file.
 *
 * Emits `v x y z` rows, `vn`/`vt` rows from `point_data["obj:vn"]`/
 * `["obj:vt"]` if present, and `f` face rows grouped by cell block
 * (triangle/quad/polygon), 1-based indices. Group (`g`) lines are emitted
 * per distinct value found in `cell_data["obj:group_ids"]`, if present.
 *
 * @param rPath filesystem path to the .obj file to create/overwrite
 * @param rMesh the mesh to write
 * @throws WriteError on an unsupported cell type
 * @note reads `point_data["obj:vn"]`, `point_data["obj:vt"]`,
 *       `cell_data["obj:group_ids"]` if present
 */
MESHIOPLUSPLUS_API void write_obj(const std::string& rPath, const Mesh& rMesh);

/**
 * @brief Read a Wavefront OBJ (.obj) file into a Mesh.
 *
 * Parses `v` (points), `vn`/`vt` (stored raw, arbitrary column count),
 * and `f` (faces; only the leading vertex index of each `i[/t[/n]]` token
 * is kept — the `/vt`/`/vn` index references are discarded entirely).
 * Faces are grouped by vertex count into `triangle`/`quad`/`polygon` cell
 * blocks; a run of same-count faces breaks into a new block whenever the
 * count changes **or** a `g` line appears, even if the count didn't
 * change. `g <name>` increments a running group-id counter (starting at
 * -1 for faces before any `g` line); only the id survives, not the name.
 * Empty trailing groups are dropped after the full file is scanned.
 *
 * @param rPath filesystem path to the .obj file to read
 * @return the read Mesh, with `point_data["obj:vn"]` (if any `vn` lines
 *         were seen), `point_data["obj:vt"]` (if any `vt` lines were
 *         seen), and `cell_data["obj:group_ids"]` (one int array per cell
 *         block, the originating group id, `-1` if before the first `g`)
 * @throws ReadError on a malformed file
 */
MESHIOPLUSPLUS_API Mesh read_obj(const std::string& rPath);

}  // namespace meshioplusplus

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
 * @file ensight.hpp
 * @brief EnSight Gold (.case/.geo) C++ reader/writer.
 *
 * EnSight Gold stores a dataset as a small `.case` index file plus a
 * geometry file (conventionally `.geo`) and, optionally, one file per
 * `VARIABLE` entry. The `.case` `FORMAT`/`GEOMETRY` sections (the file must
 * declare `type: ensight gold`) and the Gold geometry file in both ASCII and
 * C-binary form (the leading `"C Binary"` 80-char record selects binary;
 * `"Fortran Binary"` is rejected) are read; since v11.3.0 (roadmap §1 tier
 * B1) so are `TIME` (a transient file's step times) and `VARIABLE` (per-node/
 * per-element scalar/vector `point_data`/`cell_data`) — see `read_ensight`'s
 * own doc comment. Binary files use 32-bit ints/floats in the writing
 * machine's byte order; the reader auto-detects a foreign byte order from
 * the plausibility of the part-number/node-count records and byte-swaps
 * accordingly.
 *
 * Element keywords `point`, `bar2/3`, `tria3/6`, `quad4/8`, `tetra4/10`,
 * `pyramid5/13`, `penta6/15`, `hexa8/20` map to the corresponding meshio
 * cell types; ragged `nsided`/`nfaced` sections are read into polygon /
 * polyhedron blocks (grouped by node count, the openfoam convention). Node
 * ordering matches meshio for every type except `penta15`, which differs
 * from meshio's `wedge15` by the involution
 * `{0,2,1,3,5,4, 8,7,6, 11,10,9, 12,14,13}` (the same map VTK's EnSight
 * readers apply). Per the Gold specification connectivity is **positional**
 * (1-based index into the part's coordinate list); `node id given/ignore`
 * id arrays are present in the file but skipped.
 *
 * Multi-part files are concatenated into one point array; every per-part
 * element section becomes its own cell block, and when the file has two or
 * more parts the owning part number is recorded as the integer cell_data
 * field `"ensight:part"`. The writer emits a single part (`node id assign`,
 * `element id assign`) and drops point/cell/field data (mesh-only scope) --
 * the roadmap's variable-*reading* item leaves variable-*writing* to §1.
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/read_options.hpp"

namespace meshioplusplus {

/**
 * @brief Write a mesh as an EnSight Gold `.case` + `.geo` sibling pair.
 *
 * `rPath` may name either sibling (`.case` or `.geo`); the other is derived
 * from the shared stem and both files are written. The geometry is a single
 * part (id 1, description "Mesh") with `node id assign` / `element id
 * assign` numbering; coordinates are always emitted with three components
 * (z padded with 0 for 2D meshes). ASCII floats use `%12.5e` (EnSight caps
 * usable precision at 6 significant digits); binary emits 32-bit
 * ints/floats in host byte order with 80-char string records.
 *
 * @param rPath filesystem path to either the `.case` or `.geo` sibling
 * @param rMesh the mesh to write
 * @param binary true for C-binary geometry, false for ASCII
 * @throws WriteError if a file cannot be opened, the mesh contains a cell
 *         type without an EnSight keyword, a ragged (polygon/polyhedron)
 *         block is present (not written in v1), points have more than three
 *         components, or a binary mesh exceeds 32-bit counts
 * @note point/cell/field data are not written (geometry-only scope).
 */
MESHIOPLUSPLUS_API void write_ensight(const std::string& rPath, const Mesh& rMesh, bool binary);

/**
 * @brief Read an EnSight Gold `.case` file (or a Gold geometry file directly).
 *
 * A `.case` path is parsed for its `GEOMETRY`/`model:` entry (resolved
 * relative to the case file's directory; transient wildcard names are
 * rejected); any other path is treated as a Gold geometry file. ASCII and
 * C-binary geometries are both handled, including foreign-endian binaries.
 *
 * @param rPath filesystem path to a `.case` file or a Gold `.geo` file
 * @return the read Mesh (parts concatenated; one cell block per per-part
 *         element section; `"ensight:part"` cell_data when the file has
 *         two or more parts)
 * @throws ReadError on malformed input, a non-Gold case file, Fortran-binary
 *         geometry, an unknown element keyword, or out-of-range connectivity
 */
MESHIOPLUSPLUS_API Mesh read_ensight(const std::string& rPath);

/**
 * @brief `read_ensight` with read options — reads a `.case` file's
 *        `VARIABLE` sections and selects one step of a transient one.
 *
 * A `scalar per node:`/`vector per node:`/`scalar per element:`/`vector per
 * element:` entry becomes `point_data`/`cell_data` under its own name, read
 * from the file its (possibly `mTimeStep`-templated) filename names. A
 * variable file mirrors the geometry file's own `part`/section structure
 * exactly (see the file doc comment), which is how a value lands on the
 * right point or cell block with no name matching against the geometry at
 * all -- and why a variable file whose part sequence does not match the
 * geometry's own is a `ReadError`, not a best-effort guess. `mTimeStep`
 * (0-based, negative counts from the end) resolves against the case file's
 * `TIME` `time values:` (only the *first* `time set:` is honoured when a
 * file has more than one) and templates every `*`-bearing filename with
 * `filename start number:` + step × `filename increment:`, zero-padded to
 * the `*` run's own width. A file with no `TIME` section has exactly one
 * step; a non-default `mTimeStep` against it is refused rather than
 * silently answering step 0. `mPointsOnly`/`mMetadataOnly`/`mDataArrays`
 * narrow which variables are read, same as every other format.
 *
 * A bare geometry file (no `.case`) never reaches any of this: there is no
 * `VARIABLE`/`TIME` section to read, so `rOptions`' data-selecting fields are
 * silently inert on that path, exactly as before this overload existed.
 *
 * @param rPath filesystem path to a `.case` file or a Gold `.geo` file
 * @param rOptions read options
 * @return the read Mesh, as the plain overload, plus point_data/cell_data
 *         for every wanted `VARIABLE` entry
 * @throws ReadError as the plain overload, plus a `mTimeStep` out of range,
 *         a `VARIABLE` kind other than the four listed above being silently
 *         skipped (not an error), or a variable file whose structure does
 *         not match the geometry's
 */
MESHIOPLUSPLUS_API Mesh read_ensight(const std::string& rPath, const ReadOptions& rOptions);

/**
 * @brief Summarize an EnSight `.case` file's available time steps.
 *
 * Reads only the `.case` file's `TIME` section for `mTimeValues`; the mesh
 * shape (`mNumPoints`/`mCellBlocks`) still needs a full geometry read (no
 * native header-only shape scan, unlike CGNS/Gmsh 4.1), so
 * `mFellBackToFullRead` is always `true` here -- the same shape Exodus's own
 * metadata override has. A bare geometry file (no `.case`, hence no `TIME`
 * section to read) throws, which lets `registry_read_metadata`'s fallback
 * take over and answer from a plain full read instead.
 *
 * @param rPath filesystem path to a `.case` file
 * @param rOptions unused (metadata carries no timestep of its own to select)
 * @return the file's time values, with the geometry's shape from a full read
 * @throws ReadError when `rPath` is not a `.case` file, or as `read_ensight`.
 */
MESHIOPLUSPLUS_API MeshMetadata read_ensight_metadata(const std::string& rPath,
                                                      const ReadOptions& rOptions);

}  // namespace meshioplusplus

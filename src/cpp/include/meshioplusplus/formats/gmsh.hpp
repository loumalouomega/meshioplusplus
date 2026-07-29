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
 * @file gmsh.hpp
 * @brief Gmsh mesh format (.msh, versions 2.2 and 4.1) C++ reader/writer.
 *
 * `$MeshFormat` (`version filetype datasize`; `filetype` 0=ascii, 1=binary,
 * with a 4-byte endianness-detection integer `1` for binary) is read first
 * and picks the reader: `"2"`/`"2.2"` -> the 2.2 path, `"4"`/`"4.1"` -> the
 * 4.1 path. The C++ reader (@ref read_gmsh) handles **versions 2.2 and 4.1
 * only** — version 4.0 (whose `$Entities` layout differs: 6 bounding-box
 * doubles even for points, and no bounding-entity lists) and `$Periodic`
 * records always throw and defer to the Python reader (see
 * doc/formats/gmsh.md#quirks-limitations).
 *
 * **Version 2.2**: `$PhysicalNames`; `$Nodes` (ascii `id x y z` rows, or
 * binary `(int32 id, 3xdouble)`); `$Elements` (ascii `id type ntags
 * tag1..tagN node1..nodeK` per line, binary per-block `elem_type num_elems
 * num_tags` header then flat int32 rows). The first two element tags become
 * `cell_data["gmsh:physical"]`/`cell_data["gmsh:geometrical"]`.
 *
 * **Version 4.1** restructures node/element blocks: `$Nodes` header
 * `numEntityBlocks numNodes minNodeTag maxNodeTag`, per block `entityDim
 * entityTag parametric numNodesInBlock` followed by a node-tag list then a
 * matching coordinate list (tags may be sparse/out of order, requiring a
 * tag->index remap); `$Elements` likewise groups per entity block, with rows
 * of `elementTag node1..nodeK`. `$Entities` maps each `(entityDim, entityTag)`
 * to its physical tags and to the entities bounding it — it is the **only**
 * place 4.1 records physical-group membership, so `cell_data["gmsh:physical"]`
 * (and therefore every named region) comes from there rather than from the
 * element rows the way 2.2 does. `point_data["gmsh:dim_tags"]` (an `(N,2)`
 * `(entity_dim, entity_tag)` array) and the bounding entities (@ref GmshInfo,
 * surfacing in Python as `cell_sets["gmsh:bounding_entities"]`) are
 * v4.1-only concepts.
 *
 * Five element types need a node-order permutation between Gmsh and
 * meshio++ (`tetra10`, `hexahedron20`, `hexahedron27`, `wedge15`,
 * `pyramid13` — see doc/formats/gmsh.md for the exact permutation arrays);
 * everything else uses natural order. The C++ type table covers a curated
 * subset up through roughly `hexahedron125`/`tetra286` — not the full
 * ~110-entry Python table — so a file referencing a higher-order type
 * outside that subset falls back to Python transparently. `field_data` maps
 * from `$PhysicalNames` as `[phys_num, phys_dim]`; `mesh.gmsh_periodic` (a
 * mesh-level attribute, not a data-dict key) is only ever populated by the
 * Python reader.
 */

// System includes
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/read_options.hpp"

namespace meshioplusplus {

/**
 * @brief Gmsh data that has no place on a `Mesh` (the `MedInfo`/`ExodusInfo`
 *        side-channel pattern).
 *
 * A reader overload fills one, a writer overload consumes one, and the shared
 * registry -- and therefore every flat binding (WASM, C API, Fortran, Julia,
 * R) -- passes none: a documented gap, not a silent loss.
 */
struct MESHIOPLUSPLUS_API GmshInfo {
    /**
     * @brief Format 4.1 `$Entities` bounding-entity tags, one entry per cell
     *        block (the tags of the entity that block belongs to).
     *
     * These are the boundary of the block's geometric entity, one dimension
     * down, and the tags are **signed** -- the sign carries the boundary's
     * orientation. That is why they cannot be a @ref Region (whose entries are
     * sorted non-negative indices) and why they need this channel at all. On
     * the Python side they surface as `cell_sets["gmsh:bounding_entities"]`.
     *
     * Empty for dimension-0 blocks (points have no boundary) and for files
     * with no `$Entities` section at all.
     */
    std::vector<std::vector<std::int32_t>> mBoundingEntities;
};

/**
 * @brief Write `mesh` to `path` as a Gmsh 2.2 .msh file (ascii or binary).
 *
 * Emits `$MeshFormat` (version "2.2"), `$PhysicalNames` (from
 * `field_data`), `$Nodes`, and `$Elements` with `gmsh:physical`/
 * `gmsh:geometrical` as the first two element tags. Applies the gmsh <->
 * meshio++ node-order permutation for `tetra10`/`hexahedron20`/
 * `hexahedron27`/`wedge15`/`pyramid13`.
 *
 * @param rPath filesystem path to write
 * @param rMesh the mesh to write
 * @param binary write node/element bodies as binary (`true`, with the
 *        endianness-detection integer) or ASCII (`false`)
 * @throws WriteError if a cell block's type has no Gmsh type-code mapping
 * @note reads/writes `cell_data["gmsh:physical"]`/`cell_data["gmsh:geometrical"]`
 *       and `field_data` (as `$PhysicalNames`)
 * @note the shim only attempts this C++ path when `float_fmt == ".16e"` and
 *       `mesh.gmsh_periodic` is unset
 */
MESHIOPLUSPLUS_API void write_gmsh22(const std::string& rPath, const Mesh& rMesh, bool binary);

/**
 * @brief Write `mesh` to `path` as a Gmsh 4.1 .msh file (ascii or binary).
 *
 * When the mesh carries `point_data["gmsh:dim_tags"]`, a full `$Entities`
 * section is emitted (one entity per unique `(dim, tag)` node pair, its
 * physical tag taken from the cell block living on it) and `$Nodes` is split
 * into one block per entity — which is what makes physical-group *membership*
 * survive a 4.1 round-trip. Without `gmsh:dim_tags` there is no entity
 * structure to describe, so no `$Entities` is written and `$Nodes` is a single
 * block, byte for byte as before.
 *
 * Entity bounding boxes are written as zeros (as the Python reference writer
 * does); gmsh recomputes them on load.
 *
 * @param rPath filesystem path to write
 * @param rMesh the mesh to write
 * @param binary write node/element bodies as binary (`true`) or ASCII
 *        (`false`)
 * @throws WriteError if a cell block's type has no Gmsh type-code mapping, an
 *         entity dimension is outside 0..3, or two cell blocks claim the same
 *         `(dim, entity tag)` — that entity would have no single physical tag
 * @note this overload writes no bounding entities; use the @ref GmshInfo
 *       overload to carry them
 * @note the shim only attempts this C++ path when `float_fmt == ".16e"` and
 *       the mesh has no `gmsh_periodic`
 */
MESHIOPLUSPLUS_API void write_gmsh41(const std::string& rPath, const Mesh& rMesh, bool binary);

/**
 * @brief Write a Gmsh 4.1 .msh file, including the `$Entities` bounding
 *        entities.
 *
 * Identical to @ref write_gmsh41(const std::string&, const Mesh&, bool) except
 * that each entity's bounding-entity list is taken from
 * `rInfo.mBoundingEntities` (indexed by cell block). Entries beyond the block
 * count, and dimension-0 entities, are ignored.
 */
MESHIOPLUSPLUS_API void write_gmsh41(const std::string& rPath, const Mesh& rMesh, bool binary,
                                     const GmshInfo& rInfo);

/**
 * @brief Read a Gmsh .msh file (versions 2.2 and 4.1 only).
 *
 * Dispatches on the `$MeshFormat` version string; parses `$PhysicalNames`,
 * `$Nodes`, `$Elements` (applying the gmsh <-> meshio++ node-order
 * permutation where needed), and, for 4.1, `$Entities`/per-entity node and
 * element blocks with node-tag->index remapping.
 *
 * @param rPath filesystem path to read
 * @return the read Mesh, with `cell_data["gmsh:physical"]`/
 *         `cell_data["gmsh:geometrical"]` (2.2: the first two element tags;
 *         4.1: the block's entity tag and its `$Entities` physical tag, the
 *         latter present only when the file tags any entity at all, and 0 for
 *         the untagged blocks), `point_data["gmsh:dim_tags"]` (v4.1 only),
 *         `field_data` from `$PhysicalNames`, and one `Cell` region per named
 *         physical group
 * @throws ReadError for anything not handled by the C++ path — version not
 *         2.2/4.1 (4.0's `$Entities` layout differs), `$Periodic` records,
 *         a Gmsh element type outside the curated type-code subset, or
 *         parametric nodes — so the Python reader can take over
 * @note the C++ reader never populates `mesh.gmsh_periodic`; only the
 *       Python fallback does, for files containing `$Periodic`
 * @note this overload discards the `$Entities` bounding-entity tags; use the
 *       @ref GmshInfo overload to keep them
 */
MESHIOPLUSPLUS_API Mesh read_gmsh(const std::string& rPath, const ReadOptions& rOpts = {});

/**
 * @brief Read a Gmsh .msh file, keeping the data that has no place on a `Mesh`.
 *
 * Identical to @ref read_gmsh(const std::string&, const ReadOptions&) except
 * that `rInfo` receives the format 4.1 `$Entities` bounding-entity tags.
 *
 * @param rPath filesystem path to read
 * @param rInfo receives the side-channel data (cleared of nothing — append-only
 *        into whatever the caller passes)
 * @param rOpts selective-read options
 * @return the read Mesh
 */
MESHIOPLUSPLUS_API Mesh read_gmsh(const std::string& rPath, GmshInfo& rInfo,
                                  const ReadOptions& rOpts = {});

/**
 * @brief Summarize a `.msh` without reading its node coordinates or connectivity.
 *
 * **Format 4.1 only.** 4.1 groups elements into typed blocks whose headers carry
 * the type and count, so the summary walks block headers and skips each
 * payload -- by exact byte arithmetic for binary, line counts for ascii.
 * Format 2.2 stores a type per element, so there is no cheap path to have;
 * `read_gmsh_metadata` throws for it and `registry_read_metadata` falls back to
 * a full read (reporting `mFellBackToFullRead`).
 *
 * @param rPath filesystem path to summarize
 * @return the summary; `mHasBBox` is false (reading it would defeat the point)
 * @throws ReadError for format 2.2, `$Entities`/`$Periodic`, or an unsupported
 *         element type -- exactly what `read_gmsh` rejects
 */
MESHIOPLUSPLUS_API MeshMetadata read_gmsh_metadata(const std::string& rPath, const ReadOptions& rOpts = {});

}  // namespace meshioplusplus

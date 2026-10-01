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
 * @file mdpa.hpp
 * @brief Kratos Multiphysics MDPA (`.mdpa`) ascii mesh C++ reader/writer.
 *
 * MDPA is a line-oriented, block-structured ascii format: every section is
 * `Begin <Block> [args]` ... `End <Block>`, `//` starts a comment (to the end
 * of the line, and a header may carry one with no separating space, as in
 * `Begin Elements Element3D3N// group`), and tabs are as good as spaces.
 * The blocks this reader/writer understands:
 *
 *  - `Begin Nodes` — `id x y z` rows, or bare `x y z` rows (auto-detected by
 *    column count; a bare row takes its position as its id). Ids may be
 *    arbitrary — gapped and non-monotonic both read, as a real Kratos deck
 *    left by a SubModelPart extraction or an entity removal needs — and
 *    connectivity, `NodalData` and `SubModelPartNodes` resolve through a
 *    file-id → row map built lazily on the first id that is not `row + 1`
 *    (`abaqus.cpp`'s `mPointIds` pattern). Points come back in **file order**,
 *    never sorted by id; a duplicate id is a `ReadError`. **Original ids
 *    survive a write**, too (see #kMdpaIdName): whenever they were not already
 *    the trivial `1..n` a fresh write would produce anyway, they are attached
 *    as ordinary `point_data`/`cell_data`, which `write_mdpa` reads back.
 *  - `Begin Elements <KratosName>` / `Begin Conditions <KratosName>` —
 *    `id property_id n1 n2 ...` rows. The Kratos entity name resolves to a
 *    meshio cell type through `backends/kratos_names.hpp`
 *    (`cell_type_from_kratos_name`), with a longest-suffix fallback so
 *    application-specific names such as `SmallDisplacementElement3D4N` still
 *    resolve (via their `Element3D4N` suffix). A new cell block is started
 *    whenever the type differs from the previous one, exactly as the Python
 *    reference does, so block order follows the file. Property ids become
 *    Int64 `cell_data["gmsh:physical"]` — the name the Python reference uses.
 *    Element and condition ids (each their own independent 1-based counter, in
 *    file order across every block of that kind) are preserved the same way as
 *    node ids -- see #kMdpaIdName.
 *  - `Begin ModelPartData` — `KEY value` pairs, kept as one-element Float64
 *    `field_data` entries.
 *  - `Begin Properties <id>` — the material data. Its body has no place on the
 *    `Mesh` (`NDArray` has no string dtype, so `CONSTITUTIVE_LAW
 *    LinearElastic3DLaw` is not expressible as `field_data`), so it travels in
 *    the #MdpaInfo side channel — the `MedInfo`/`ExodusInfo` pattern. Reading
 *    without one parses the block and drops it with a warning; the flat
 *    bindings (C API, Fortran, Julia, R, WASM) do exactly that, a documented
 *    gap rather than a silent loss.
 *  - `Begin NodalData <VAR>[n]` — sparse per-node values → `point_data`
 *    (Float64, NaN where a node is not listed). A leading `0`/`1` "fixed"
 *    column is detected and stored as Int64 `"<VAR>_fixed_status"` (`-1` where
 *    unspecified); a zero-component variable becomes an Int64 membership flag.
 *  - `Begin ElementalData` / `Begin ConditionalData` — the same, per entity id
 *    → `cell_data` (one array per cell block, the repo-wide convention).
 *  - `Begin SubModelPart <Name>` — `SubModelPartNodes` become a
 *    `RegionKind::Point` region and `SubModelPartElements`/`Conditions` a
 *    `RegionKind::Cell` one (global, block-major cell indices), both named
 *    after the sub-model-part. Nested parts are flattened to a `parent/child`
 *    path name, which round-trips as a single name.
 *
 * Kratos orders the nodes of `hexahedron20` and `hexahedron27` differently
 * from VTK/meshio; the permutation is applied on read and undone on write, so
 * an MDPA file stays valid for Kratos.
 *
 * @note cell_data key produced/consumed: `"gmsh:physical"` (the Kratos
 *       property id of each element/condition); `point_data`/
 *       `cell_data[kMdpaIdName]` (`"mdpa:id"`) for original node/entity ids,
 *       when they were not already the trivial `1..n` renumbering — see
 *       #kMdpaIdName.
 *
 * ## The blocks the `Mesh` cannot hold
 *
 * MDPA also carries content with no place on a `Mesh`: top-level `Begin Table`
 * blocks, `Begin Geometries`, `Begin Mesh <id>`, a non-numeric
 * `ModelPartData` value, `SubModelPartData`/`SubModelPartTables` bodies, and
 * blocks such as `Begin Constraints` that no tool but Kratos interprets. The
 * **#MdpaInfo overloads read and write every one of them**
 * (`read_mdpa(path, info)` / `write_mdpa(path, mesh, info)`, and
 * `mio_read_with_info` / `mio_write_with_info` on the flat ABI): tables,
 * geometries, mesh blocks, text data and sub-model-part data are parsed into
 * #MdpaInfo, and any other top-level block is kept verbatim as an
 * #MdpaRawBlock, so it round-trips byte for byte.
 *
 * Without an #MdpaInfo there is nowhere to put that content, so rather than
 * silently dropping it the info-less overloads (and so the registry, `mio_read`
 * and the native CLI) **throw `ReadError` naming the construct**.
 * **`ReadOptions::mLenient` downgrades every one of those to a warning plus a
 * skip**, and what was skipped is recorded in `MdpaInfo::mSkippedConstructs`.
 * Geometry and constraint membership is carried in #MdpaSubModelPart as raw
 * file ids, not mesh regions. Two things still throw even under `mLenient`,
 * because skipping them would return a mesh that is quietly wrong rather than
 * merely incomplete:
 *
 *  - a duplicate node id (two coordinate rows claiming one id is
 *    unrepresentable, not merely incomplete);
 *  - a malformed row, an unknown entity name, or connectivity naming a node
 *    the `Nodes` block does not define.
 *
 * The writer emits the mesh-level blocks (`ModelPartData` from scalar
 * `field_data`, `Properties`, `Nodes`, `Elements`/`Conditions`,
 * `NodalData`/`ElementalData`/`ConditionalData`, `SubModelPart`s from named
 * regions) and, given an #MdpaInfo, its `Table`, `Geometries`, `Mesh` and raw
 * blocks. `RegionKind::Side` regions are dropped with a warning (MDPA has no
 * facet-set concept).
 */

// System includes
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/properties.hpp"
#include "meshioplusplus/read_options.hpp"

namespace meshioplusplus {

/** @brief One `Begin Properties <id>` block (see #PropertySet). */
using MdpaProperties = PropertySet;

/** @brief One `KEY value` entry of a properties block (see #PropertyValue). */
using MdpaProperty = PropertyValue;

/**
 * @brief The original-id carrier: `point_data`/`cell_data["mdpa:id"]`.
 *
 * The uniform mesh API has no id-translation layer -- `Mesh::Points()`/`Conn()`
 * are dense 0-based arrays where "point index `i`" *is* row `i` -- so there is
 * nowhere on the `Mesh` itself to remember a file's original node/element/
 * condition numbering. `read_mdpa` attaches it as ordinary data instead: a
 * point_data array of Int64 node ids (one per point, in read order) and/or a
 * per-block cell_data array of Int64 element/condition ids, **only when those
 * ids were not already the trivial `1..n` renumbering `write_mdpa` would
 * produce anyway** -- so a sequential (or id-less) deck is completely
 * unaffected and a re-write of it stays byte-identical to before. `write_mdpa`
 * honours the array when present (falling back to the old renumbering when it
 * is absent, the wrong length, or the wrong dtype) and throws `WriteError` on
 * a duplicate value, since writing one would silently produce an invalid file.
 * Riding as plain data means it survives (and is renumbered by) the ordinary
 * mesh operations for free -- crop/split/etc. carry `point_data`/`cell_data`
 * through their existing row-selection machinery with no MDPA-specific code.
 */
inline constexpr const char* kMdpaIdName = "mdpa:id";

/**
 * @brief The Kratos entity spelling of one cell block.
 *
 * The name is kept rather than re-derived because MDPA's is
 * application-specific: `SmallDisplacementElement3D4N` resolves to `tetra` via
 * the longest-suffix fallback, but `kratos_element_name(Tetrahedra4)` only ever
 * gives back the canonical `Element3D4N`, so a round trip through the C++ core
 * silently downgraded every deck's element names before this existed.
 *
 * `mIsCondition` is recorded rather than inferred from the name: the writer's
 * fallback rule ("a block whose default Kratos element name is 2-D is a
 * Condition") is a heuristic about cell types, and guessing from an
 * application-specific name would be a second, different guess.
 */
struct MdpaEntityName {
    /** @brief The name as it appeared, or empty to derive it from the cell type. */
    std::string mName;
    /** @brief Whether the block was a `Conditions` block rather than `Elements`. */
    bool mIsCondition = false;
};

/**
 * @brief One `Begin Geometries <Name>` run: connectivity with no property id.
 *
 * Kratos geometries (the entities CAD/IGA applications and `CreateNewGeometry`
 * use) are not elements, so they are kept apart from the mesh's cell blocks:
 * an operation would otherwise treat them as cells. A new block starts
 * wherever the name changes, so every name the file spelled survives.
 */
struct MdpaGeometryBlock {
    /** @brief The Kratos geometry name as read (`Triangle3D3`); empty to derive it. */
    std::string mName;
    /** @brief The meshio cell type (`"triangle"`). */
    std::string mType;
    /** @brief Int64 `(n, k)`: 0-based point rows, in meshio node order. */
    NDArray mConn;
    /** @brief The file's geometry ids, one per row. */
    std::vector<std::int64_t> mIds;
};

/**
 * @brief One `Begin Mesh <id>` block (a Kratos multi-level mesh).
 *
 * Node members are point rows, like region entries; element and condition
 * members are the file's own ids, written back verbatim -- the Python
 * reference's rule, since a Mesh block names entities rather than rows.
 */
struct MdpaMeshBlock {
    /** @brief The mesh id from the header (never 0: Kratos reserves it). */
    std::int64_t mId = 0;
    /** @brief The `MeshData` entries, in file order. */
    std::vector<PropertyValue> mData;
    /** @brief `MeshNodes`, as 0-based point rows. */
    std::vector<std::int64_t> mNodes;
    /** @brief `MeshElements`, as file element ids. */
    std::vector<std::int64_t> mElementIds;
    /** @brief `MeshConditions`, as file condition ids. */
    std::vector<std::int64_t> mConditionIds;
};

/**
 * @brief The non-membership content of one `Begin SubModelPart`.
 *
 * Membership (nodes, elements, conditions) is on the mesh as regions, which
 * operations remap; this holds only what a region cannot: the part's
 * `SubModelPartData` entries and table, geometry and constraint ids. These ids
 * stay in file order and are not remapped by mesh operations; constraints
 * themselves remain opaque #MdpaRawBlock content.
 */
struct MdpaSubModelPart {
    /** @brief The hierarchical name (`"parent/child"`), as the region is named. */
    std::string mName;
    /** @brief The `SubModelPartData` entries, in file order. */
    std::vector<PropertyValue> mData;
    /** @brief The `SubModelPartTables` ids, in file order. */
    std::vector<std::int64_t> mTables;
    /** @brief The `SubModelPartGeometries` file ids, in file order. */
    std::vector<std::int64_t> mGeometryIds;
    /** @brief The `SubModelPartConstraints` file ids, in file order. */
    std::vector<std::int64_t> mConstraintIds;
};

/**
 * @brief A top-level block kept verbatim, such as `Begin Constraints`.
 *
 * What only Kratos interprets -- master-slave constraints, and any block a
 * future Kratos adds -- is carried as text so a round trip loses nothing.
 */
struct MdpaRawBlock {
    /** @brief The header line with any comment removed (`"Begin Constraints ..."`). */
    std::string mHeader;
    /** @brief The body lines, verbatim, each ending in `'\n'`. */
    std::string mBody;
    /** @brief The terminating line (`"End Constraints"`). */
    std::string mEnd;
};

/**
 * @brief MDPA content the C++ `Mesh` cannot hold, carried alongside it.
 *
 * The `MedInfo`/`ExodusInfo`/`OpenFoamInfo` pattern: a reader overload fills
 * one, a writer overload consumes one, and the registry passes none. The flat
 * bindings reach it through `mio_read_with_info` / `mio_write_with_info`.
 *
 * Round-tripping a real deck is exactly `read_mdpa(in, info)` followed by
 * `write_mdpa(out, mesh, info)`; both halves are additive, so `read_mdpa(in)`
 * and `write_mdpa(out, mesh)` behave as they always did.
 */
struct MdpaInfo {
    /** @brief Every `Begin Properties <id>` block, in file order. */
    std::vector<MdpaProperties> mProperties;
    /** @brief One entry per cell block, in mesh block order. */
    std::vector<MdpaEntityName> mEntityNames;
    /**
     * @brief What `ReadOptions::mLenient` skipped, in the order it was met.
     *
     * Empty after a strict read by construction -- a strict read either has
     * nothing to skip or throws. Reporting it is what keeps "lenient" from
     * meaning "silently lossy".
     */
    std::vector<std::string> mSkippedConstructs;
    /**
     * @brief The non-numeric `ModelPartData` entries (`mText`), in file order.
     *
     * Numeric entries stay one-element Float64 `field_data`, as they always
     * were; `NDArray` has no string dtype for the rest.
     */
    std::vector<PropertyValue> mModelPartData;
    /** @brief Every top-level `Begin Table` (`mIsTable`; `mKey` = header arguments). */
    std::vector<PropertyValue> mTables;
    /** @brief Every `Begin Geometries` run, in file order. */
    std::vector<MdpaGeometryBlock> mGeometries;
    /** @brief Every `Begin Mesh <id>` block, in file order. */
    std::vector<MdpaMeshBlock> mMeshBlocks;
    /** @brief Sub-model-parts with data or table/geometry/constraint ids, in file order. */
    std::vector<MdpaSubModelPart> mSubModelParts;
    /** @brief Every other top-level block, verbatim, in file order. */
    std::vector<MdpaRawBlock> mRawBlocks;
};

/**
 * @brief Write a mesh to a Kratos MDPA file.
 *
 * Emits, in order: `ModelPartData` (one line per one-element `field_data`
 * array), the `Properties` blocks, `Nodes` (always three coordinates —
 * 2-D meshes are padded with `z = 0`), then one `Elements` or `Conditions`
 * block per cell block, `NodalData`/`ElementalData`/`ConditionalData` for the
 * remaining data arrays and one `SubModelPart` per named region.
 *
 * With no #MdpaInfo the properties have no bodies, but one **empty** block is
 * emitted per distinct id the entity rows actually reference, ascending. That
 * is a fix, not a change of policy: the rows have always written their
 * `gmsh:physical` value as the property id while the header was hard-coded to a
 * single `Properties 0`, so a tagged mesh produced a file referencing
 * undeclared properties, which Kratos's own `ModelPartIO` rejects. A mesh whose
 * ids are all 0 — every mesh with no `gmsh:physical` — still emits exactly
 * `Begin Properties 0` / `End Properties` and so is byte-identical to before.
 *
 * A block is written as `Conditions` when the default Kratos *element* name
 * for its cell type is a 2-D one (`Element2D4N`, ...) and as `Elements`
 * otherwise — the rule the Python reference applies for a mesh with no
 * physical tags, which is what keeps a quad mesh writing as
 * `SurfaceCondition3D4N`. Node/element/condition ids default to `row + 1`
 * (elements and conditions each their own 1-based counter over every block of
 * that kind, in mesh order) **unless the mesh carries #kMdpaIdName**, in which
 * case those original ids are written back instead — see #kMdpaIdName for the
 * exact contract, including what counts as present/valid and the duplicate-id
 * `WriteError`. Every place an entity or node is referenced elsewhere in the
 * file (`NodalData`/`ElementalData`/`ConditionalData` row keys, `SubModelPart`
 * node/element/condition lists) uses the same resolved id, so the file is
 * always internally consistent whichever numbering was actually used. The
 * property id of a cell is its `cell_data["gmsh:physical"]` value when that
 * array exists, else 0.
 *
 * @param rPath filesystem path to write
 * @param rMesh the mesh to write
 * @throws WriteError on an unopenable output path, a ragged/polyhedron cell
 *         block (MDPA has no such entity), a cell type with no Kratos name, or
 *         a duplicate value in #kMdpaIdName (which would silently produce an
 *         invalid Kratos deck)
 * @note reads `cell_data["gmsh:physical"]` for the per-entity property id;
 *       `point_data["<VAR>_fixed_status"]` for the `NodalData` fixed column;
 *       `point_data`/`cell_data[kMdpaIdName]` for original node/entity ids.
 */
MESHIOPLUSPLUS_API void write_mdpa(const std::string& rPath, const Mesh& rMesh);

/**
 * @brief Write a mesh to MDPA, restoring the content `read_mdpa` set aside.
 *
 * Identical to the two-argument form except that `rInfo` supplies the
 * `Properties` bodies, the per-block Kratos entity names, the text
 * `ModelPartData` entries, and the `Table`, `Geometries`, `Mesh`,
 * sub-model-part data and raw blocks, written in that file position: tables
 * after `Properties`, geometries after the entity blocks, sub-model-part data
 * inside its `SubModelPart`, then `Mesh` and raw blocks last. A block with no
 * entry in `MdpaInfo::mEntityNames`, or one whose `mName` is empty, falls back
 * to the derived name exactly as the two-argument form does, so a partially
 * filled #MdpaInfo is legal.
 *
 * @param rPath filesystem path to write
 * @param rMesh the mesh to write
 * @param rInfo the side-channel content
 * @throws WriteError as the two-argument form, and on a geometry row naming a
 *         point the mesh does not have
 */
MESHIOPLUSPLUS_API void write_mdpa(const std::string& rPath, const Mesh& rMesh,
                                   const MdpaInfo& rInfo);

/**
 * @brief Read a Kratos MDPA mesh file.
 *
 * Parses the blocks listed in the file-level documentation, producing points,
 * cell blocks with `gmsh:physical` property ids, `point_data`/`cell_data` from
 * the `*Data` blocks, `field_data` from `ModelPartData` and `Point`/`Cell`
 * regions from the `SubModelPart`s.
 *
 * @param rPath filesystem path to read
 * @return the read Mesh
 * @throws ReadError on a malformed or unterminated block, on connectivity
 *         referring to a node that does not exist, and on every construct
 *         only an #MdpaInfo can hold (see "The blocks the `Mesh` cannot hold"
 *         above), since this overload has nowhere to put it
 * @note cell_data key produced: `"gmsh:physical"`.
 */
MESHIOPLUSPLUS_API Mesh read_mdpa(const std::string& rPath);

/**
 * @brief Read a Kratos MDPA mesh file honouring @p rOptions.
 *
 * Only `ReadOptions::mLenient` changes anything here; the narrowing options are
 * applied by the caller (`registry_read`, the Python `_apply_read_filter`)
 * after the read, as they are for every format with no native selective path.
 *
 * @param rPath filesystem path to read
 * @param rOptions read options; `mLenient` downgrades the unsupported-construct
 *        rejections listed above to a warning plus a skip
 * @return the read Mesh
 */
MESHIOPLUSPLUS_API Mesh read_mdpa(const std::string& rPath, const ReadOptions& rOptions);

/**
 * @brief Read a Kratos MDPA mesh file, keeping what the `Mesh` cannot hold.
 *
 * The overload a round trip needs: `rInfo` comes back carrying the `Properties`
 * bodies, the per-block Kratos entity names, every block listed under "The
 * blocks the `Mesh` cannot hold" and (under `mLenient`) the list of skipped
 * constructs, all of which `write_mdpa(path, mesh, info)` puts back, including
 * geometry and constraint membership of nested sub-model-parts.
 *
 * @param rPath filesystem path to read
 * @param rInfo out: the side-channel content; cleared first
 * @param rOptions read options (see the two-argument overload)
 * @return the read Mesh
 */
MESHIOPLUSPLUS_API Mesh read_mdpa(const std::string& rPath, MdpaInfo& rInfo,
                                  const ReadOptions& rOptions = {});

}  // namespace meshioplusplus

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
 *  - `Begin Nodes` — `id x y z` rows. Node ids must be `1..n` in order (see
 *    the limitations below); connectivity is 1-based into them.
 *  - `Begin Elements <KratosName>` / `Begin Conditions <KratosName>` —
 *    `id property_id n1 n2 ...` rows. The Kratos entity name resolves to a
 *    meshio cell type through `backends/kratos_names.hpp`
 *    (`cell_type_from_kratos_name`), with a longest-suffix fallback so
 *    application-specific names such as `SmallDisplacementElement3D4N` still
 *    resolve (via their `Element3D4N` suffix). A new cell block is started
 *    whenever the type differs from the previous one, exactly as the Python
 *    reference does, so block order follows the file. Property ids become
 *    Int64 `cell_data["gmsh:physical"]` — the name the Python reference uses.
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
 *       property id of each element/condition).
 *
 * ## Limitations (deliberate, and reported by throwing)
 *
 * The C++ `Mesh` has no place for MDPA's non-mesh metadata, so rather than
 * silently dropping it the reader **throws `ReadError` naming the construct**
 * — which lets the Python shim fall back to the pure-Python reference
 * (`meshioplusplus/mdpa/_mdpa.py`), whose `mesh.misc_data` does carry it:
 *
 *  - `Begin Table` (top-level), `Begin Geometries`, `Begin Mesh <id>` and
 *    `Begin Constraints` blocks;
 *  - a non-numeric `ModelPartData` value, and non-empty `SubModelPartData` /
 *    `SubModelPartTables` / `SubModelPartGeometries` / `SubModelPartConstraints`
 *    sub-blocks;
 *  - any unrecognized `Begin <Block>`.
 *
 * **`ReadOptions::mLenient` downgrades every one of those to a warning plus a
 * skip**, which is what makes a production deck readable where there is no
 * Python to fall back to (the C API, Fortran, Julia, R, WASM, the native CLI).
 * What was skipped is recorded in `MdpaInfo::mSkippedConstructs`. Two things
 * are deliberately *not* on that list and still throw even under `mLenient`,
 * because skipping them would return a mesh that is quietly wrong rather than
 * merely incomplete:
 *
 *  - node ids that are not `1..n` in ascending order (the format allows
 *    arbitrary ids; honouring them would need a renumbering the Python
 *    reference does not do);
 *  - a malformed row, an unknown entity name, or connectivity naming a node
 *    that does not exist.
 *
 * The writer emits the mesh-level blocks (`ModelPartData` from scalar
 * `field_data`, `Properties`, `Nodes`, `Elements`/`Conditions`,
 * `NodalData`/`ElementalData`/`ConditionalData`, `SubModelPart`s from named
 * regions); it never writes `Tables`, `Geometries` or `Mesh` blocks, and
 * `RegionKind::Side` regions are dropped with a warning (MDPA has no facet-set
 * concept).
 */

// System includes
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
 * @brief MDPA content the C++ `Mesh` cannot hold, carried alongside it.
 *
 * The `MedInfo`/`ExodusInfo`/`OpenFoamInfo` pattern: a reader overload fills
 * one, a writer overload consumes one, and the registry (and therefore every
 * flat binding) passes none — a documented gap, not a silent loss.
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
 * `SurfaceCondition3D4N`. Element and condition ids are two independent
 * 1-based counters. The property id of a cell is its
 * `cell_data["gmsh:physical"]` value when that array exists, else 0.
 *
 * @param rPath filesystem path to write
 * @param rMesh the mesh to write
 * @throws WriteError on an unopenable output path, a ragged/polyhedron cell
 *         block (MDPA has no such entity), or a cell type with no Kratos name
 * @note reads `cell_data["gmsh:physical"]` for the per-entity property id;
 *       `point_data["<VAR>_fixed_status"]` for the `NodalData` fixed column.
 */
MESHIOPLUSPLUS_API void write_mdpa(const std::string& rPath, const Mesh& rMesh);

/**
 * @brief Write a mesh to MDPA, restoring the content `read_mdpa` set aside.
 *
 * Identical to the two-argument form except that `rInfo` supplies the
 * `Properties` bodies and the per-block Kratos entity names. A block with no
 * entry in `MdpaInfo::mEntityNames`, or one whose `mName` is empty, falls back
 * to the derived name exactly as the two-argument form does, so a partially
 * filled #MdpaInfo is legal.
 *
 * @param rPath filesystem path to write
 * @param rMesh the mesh to write
 * @param rInfo the side-channel content (properties, entity names)
 * @throws WriteError as the two-argument form
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
 *         referring to a node that does not exist, and — by design — on every
 *         construct listed under "Limitations" above, so that a caller with a
 *         Python fallback can defer to the richer reference reader
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
 * bodies, the per-block Kratos entity names and (under `mLenient`) the list of
 * skipped constructs, all of which `write_mdpa(path, mesh, info)` puts back.
 *
 * @param rPath filesystem path to read
 * @param rInfo out: the side-channel content; cleared first
 * @param rOptions read options (see the two-argument overload)
 * @return the read Mesh
 */
MESHIOPLUSPLUS_API Mesh read_mdpa(const std::string& rPath, MdpaInfo& rInfo,
                                  const ReadOptions& rOptions = {});

}  // namespace meshioplusplus

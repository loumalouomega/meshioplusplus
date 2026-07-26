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
 *  - `Begin Properties <id>` — accepted only when empty (see limitations).
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
 *  - `Begin Table` (top-level or inside `Properties`), `Begin Geometries`,
 *    `Begin Mesh <id>` and `Begin Constraints` blocks;
 *  - a non-empty `Begin Properties <id>` body, a non-numeric `ModelPartData`
 *    value, and non-empty `SubModelPartData` / `SubModelPartTables` /
 *    `SubModelPartGeometries` / `SubModelPartConstraints` sub-blocks;
 *  - node ids that are not `1..n` in ascending order (the format allows
 *    arbitrary ids; honouring them would need a renumbering the Python
 *    reference does not do, and ignoring them silently produces a wrong mesh);
 *  - any unrecognized `Begin <Block>`.
 *
 * The writer emits the mesh-level blocks only (`ModelPartData` from scalar
 * `field_data`, an empty `Properties 0`, `Nodes`, `Elements`/`Conditions`,
 * `NodalData`/`ElementalData`/`ConditionalData`, `SubModelPart`s from named
 * regions); it never writes `Tables`, `Geometries`, `Mesh` blocks or
 * per-property data, and `RegionKind::Side` regions are dropped with a warning
 * (MDPA has no facet-set concept).
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Write a mesh to a Kratos MDPA file.
 *
 * Emits, in order: `ModelPartData` (one line per one-element `field_data`
 * array), an empty `Properties 0`, `Nodes` (always three coordinates —
 * 2-D meshes are padded with `z = 0`), then one `Elements` or `Conditions`
 * block per cell block, `NodalData`/`ElementalData`/`ConditionalData` for the
 * remaining data arrays and one `SubModelPart` per named region.
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

}  // namespace meshioplusplus

# Kratos / MDPA (`.mdpa`)

The [Kratos Multiphysics](https://github.com/KratosMultiphysics/Kratos/wiki/Input-data) model-part data format: block-structured ASCII (`Begin ... / End ...`). This is the largest and most feature-rich format meshio++ supports — it has no C++ implementation.

| | |
|---|---|
| **Format name** | `mdpa` |
| **Extensions** | `.mdpa` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshioplusplus

mesh = meshioplusplus.read("model.mdpa")
meshioplusplus.mdpa.write("out.mdpa", mesh, float_fmt=".16e", binary=False)
```

- **`float_fmt`** — coordinate format.
- **`binary`** — MDPA is ASCII-only; passing `binary=True` raises `WriteError` unconditionally.

## File structure

A single pass over `Begin <X> ... End <X>` blocks:

- **`ModelPartData`** — `key value` pairs (`//` comments stripped) → `field_data[key]` (parsed as float when possible, else kept as a string).
- **`Nodes`** — rows of either `id x y z` or bare `x y z` (auto-detected by column count).
- **`Elements <KratosType>`** / **`Conditions <KratosType>`** — the header's Kratos type name is matched by exact match first, then by longest-substring match (to avoid e.g. `"Line"` ambiguously matching inside `"Line3D2"`); each row is `id property_id node_ids...`. If the type can't be resolved from the header at all, it's inferred purely from node count. Property ids become `gmsh:physical`/`gmsh:geometrical`-style tags (MDPA reuses gmsh's tag-key convention here).
- **`Geometries <Type>`** — like Elements/Conditions but with **no property id column**; stored separately as `mesh.geometries_block` (**not** part of `mesh.cells`), a non-standard mesh-level attribute.
- **`Table <id> <var1> <var2> ...`** — rows until `End Table`; a malformed header (too few parts, non-integer id, no variables) is warned-and- skipped; stored as `field_data[f"table_{id}"] = {"variables": [...], "data": ndarray}`.
- **`Properties <id>`** — key/value pairs (auto-typed float → int-if-integer → else string) plus any nested `Table` blocks, stored under `field_data[f"properties_{id}"]`.
- **`NodalData <VAR[n]>`** / **`ElementalData`/`ConditionalData <VAR>`** — values per entity; missing entities are densified with `NaN`. An optional leading "fixed" flag column (`0`/`1`) is heuristically detected (only treated as a flag if the value is exactly 0/1 **and** more numeric values follow on the same row) — if any fixed-status is seen, a parallel `{var}_fixed_status` array is produced (sentinel `-1` = "not specified"). Scalar (0-component) variables are treated as boolean-by-presence: listed ids get `1`, unlisted get `0`.
- **`SubModelPart <Name>`** (nestable, joined with `/` for a hierarchical key e.g. `"Outer/Inner"`) — sub-blocks `SubModelPartData`, `SubModelPartTables`, `SubModelPartNodes` (0-based, validated against the point count), `SubModelPartElements`/`Conditions` (raw 1-based ids kept **unconverted**, explicitly to preserve exact round-trip values). Not implemented: `SubModelPartGeometries`, `Constraints` sub-blocks.
- **`Mesh <id> [name]`** — an alternate/coarser mesh representation; `id=0` is invalid per Kratos convention and skipped with a warning. Sub-blocks: `MeshData`, `MeshNodes` (0-based, validated), `MeshElements`/`Conditions` (raw 1-based ids kept as-is).

All of the round-trip-only bookkeeping above (element/condition/geometry id maps, SubModelPart hierarchy, alternate-mesh data) accumulates in `mesh.misc_data` — a **non-standard mesh attribute** specific to this format.

## Cell types & node ordering

The Kratos type tables (Geometries/Elements/Conditions) are large; a representative slice:

| Kratos | meshio++ | Kratos | meshio++ |
|---|---|---|---|
| `Line2D2`, `Element3D2N` | `line` | `Tetrahedra3D4`, `Element3D4N` | `tetra` |
| `Line2D3`, `LineElement3D3N` | `line3` | `Tetrahedra3D10` | `tetra10` |
| `Triangle2D3`, `Element3D3N` | `triangle` | `Hexahedra3D8`, `Element3D8N` | `hexahedron` |
| `Triangle2D6`, `Element2D6N` | `triangle6` | `Hexahedra3D20` | `hexahedron20` |
| `Quadrilateral2D4`, `Element2D4N` | `quad` | `Hexahedra3D27` | `hexahedron27` |
| `Quadrilateral2D8` | `quad8` | `Prism3D6`, `Element3D6N` | `wedge` |
| `Quadrilateral2D9` | `quad9` | `Element3D5N` | `pyramid` |
| `Point2D`, `Element2D1N` | `vertex` | `Element3D13N`/`15N` | `wedge15` |

**Quadratic hexahedron node-order permutation** (the format's key gotcha, applied only for `hexahedron20`/`hexahedron27`) — read applies the argsort of the Kratos-order array below; write applies the array itself directly (a true inverse pair):

```
hex20 kratos order: [0,1,2,3,4,5,6,7,8,11,10,9,16,19,18,17,12,13,14,15]
hex27 kratos order: [0,1,2,3,4,5,6,7,8,11,10,9,16,19,18,17,12,15,14,13,
                      20,23,21,24,22,25,26]
```

All other cell types are assumed to already share meshio++'s VTK-style ordering (no permutation applied).

## Data mapping

MDPA has an unusually rich set of data keys, several structured differently from every other format meshio++ supports:

- `field_data[key]` — `ModelPartData` scalars.
- `field_data[f"table_{id}"]` — `{"variables": [...], "data": ndarray}`.
- `field_data[f"properties_{id}"]` — a dict, possibly containing nested `table_<id>` entries.
- `point_data[VAR]` — Kratos variable names verbatim (e.g. `TEMPERATURE`, `DISPLACEMENT` as an `(n,3)` array from a `DISPLACEMENT[3]` header).
- `point_data[f"{VAR}_fixed_status"]` — sentinel `-1`/`0`/`1`.
- **`cell_data[<meshio_type>]["gmsh:physical"]`/`["gmsh:geometrical"]`/`[VAR]`** — unlike every other meshio++ format, MDPA's cell_data is nested **by cell type name** as an inner dict (`{cell_type: {var: array}}`), not the usual flat `{var: [array_per_block]}` convention.
- `mesh.misc_data` — non-standard attribute: `reader_element_ids_info`, `reader_condition_ids_info`, `mdpa_geometry_ids_info`, `submodelpart_info`, `meshes`.
- `mesh.geometries_block` — non-standard attribute, a list of `CellBlock`s from `Begin Geometries`.

## Quirks & limitations

- The `cell_data` nested-by-type structure (`{cell_type: {var: array}}`) is a genuine structural departure from meshio++'s usual flat convention — code consuming MDPA-read meshes needs to account for this specifically.
- The h20/h27 permutation tables are applied **directly** on write (not their inverse) and via **argsort** on read — this is intentional and correct (the two operations really are exact inverses of each other), but worth internalizing since it looks asymmetric at first glance.
- `SubModelPartElements`/`Conditions` and `MeshElements`/`Conditions` store **raw, unconverted 1-based ids** rather than remapped local indices — this assumes element/condition ids are stable across a read→write cycle (true unless entities are reordered in between).
- Malformed rows in almost every block type (bad `Table` headers, data-row/variable-count mismatches, out-of-range `SubModelPartNodes` entries) are warned-and-skipped rather than raising — MDPA parsing is deliberately lenient/best-effort given how varied real Kratos input decks are.
- Writing `ElementalData`/`ConditionalData` omits any entity that was entirely `NaN` (never had data) rather than writing `NaN` literally.
- No C++ implementation exists for this format at all.

## Notes

- `tests/python/meshes/mdpa/test_small_cube.mdpa` — a small unit-cube tet mesh.
- `tests/python/meshes/mdpa/test_submodelpart.mdpa` — a 2D quad mesh with a nested `SubModelPart` containing Nodes/Elements/Conditions/empty-Geometries/ empty-Constraints sub-blocks.
- Additional `tests/python/input/mdpa/test_*.mdpa` fixtures target node-order permutation edge cases, minimal/degenerate geometries, hierarchical SubModelParts, and varied Table layouts.
- `tests/python/test_mdpa.py` also builds many MDPA snippets **inline** (not as files) covering nearly every block type, including a `test_roundtrip_all_blocks` exercising almost all of them at once with a NaN-aware comparison helper.

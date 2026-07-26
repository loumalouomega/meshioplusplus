# Kratos / MDPA (`.mdpa`)

The [Kratos Multiphysics](https://github.com/KratosMultiphysics/Kratos/wiki/Input-data) model-part data format: block-structured ASCII (`Begin ... / End ...`). This is the largest and most feature-rich format meshio++ supports. There are **two implementations**: the pure-Python reference below — which the Python API always uses for reading — and a C++ core reader/writer covering the mesh-level blocks, which is what makes `.mdpa` reachable from the C API, Fortran, Julia, R, WebAssembly and the native CLI (see [C++ core](#c-core) below).

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

## C++ core

`meshioplusplus::read_mdpa` / `write_mdpa` (`src/cpp/src/formats/mdpa.cpp`) implement the mesh-level part of the format against the uniform mesh API, and are registered in the shared dispatch registry — so `.mdpa` now works from the [C API](../c_api.md), [Fortran](../fortran.md), [Julia](../julia.md), [R](../r.md), [WebAssembly](../wasm.md) and the native CLI. The Kratos entity-name tables are the ones already shared with the [KRATOS mesh backend](../cpp_backends.md) (`backends/kratos_names.hpp`), extended with a longest-suffix fallback so application-specific names such as `SmallDisplacementElement3D4N` resolve through their `Element3D4N` suffix.

What the C++ core maps:

| MDPA | C++ `Mesh` |
|---|---|
| `Nodes` | `points` (always 3 columns) |
| `Elements` / `Conditions` | cell blocks, in file order, plus Int64 `cell_data["gmsh:physical"]` (the property id) |
| `ModelPartData` | one-element Float64 `field_data` entries (numeric values only) |
| `NodalData` | `point_data` (+ `"<VAR>_fixed_status"`) |
| `ElementalData` / `ConditionalData` | `cell_data`, **one array per cell block** — the repo-wide convention, *not* the reference reader's nested-by-cell-type layout |
| `SubModelPart` | a `Point` and/or `Cell` [named region](../regions.md); nested parts flatten to a `parent/child` name |

Everything the C++ `Mesh` cannot hold makes the reader **throw `ReadError` naming the construct** rather than dropping it silently: `Table`, `Geometries`, `Mesh <id>` and `Constraints` blocks, a non-empty `Properties` body, a non-numeric `ModelPartData` value, non-empty `SubModelPartData`/`SubModelPartTables`, node ids that are not `1..n` in order, and any unrecognized block. The writer emits `ModelPartData`, an empty `Properties 0`, `Nodes`, `Elements`/`Conditions`, the `*Data` blocks and one `SubModelPart` per named region; `Side` regions are dropped with a warning (MDPA has no facet-set concept).

**The Python `meshioplusplus.mdpa.read` deliberately does not use it.** Only the reference reader produces `mesh.misc_data`, `mesh.geometries_block` and the nested-by-cell-type `cell_data` this page documents, so preferring the C++ reader would silently change the Python API's output; reach it explicitly with `meshioplusplus._core.mdpa_read(path)` when you want the standard layout. `meshioplusplus.mdpa.write` *does* use the C++ writer, for real file paths and meshes that carry none of the MDPA extras (`misc_data`, `geometries_block`, `field_data`), falling back to the reference writer otherwise.

## Notes

- `tests/python/meshes/mdpa/test_small_cube.mdpa` — a small unit-cube tet mesh.
- `tests/python/meshes/mdpa/test_submodelpart.mdpa` — a 2D quad mesh with a nested `SubModelPart` containing Nodes/Elements/Conditions/empty-Geometries/ empty-Constraints sub-blocks.
- Additional `tests/python/input/mdpa/test_*.mdpa` fixtures target node-order permutation edge cases, minimal/degenerate geometries, hierarchical SubModelParts, and varied Table layouts.
- `tests/python/test_mdpa.py` also builds many MDPA snippets **inline** (not as files) covering nearly every block type, including a `test_roundtrip_all_blocks` exercising almost all of them at once with a NaN-aware comparison helper.

## Properties, entity names and lenient reads (v9.1.0)

The C++ reader used to throw on a non-empty `Begin Properties` body — which
essentially every production deck has — so the C++ path handled geometry-only
files. Three changes fix that.

### `MdpaInfo`: the side channel

Properties bodies and per-block Kratos entity names have no place on the `Mesh`:
`NDArray` has ten numeric dtypes and no string one, so `CONSTITUTIVE_LAW
LinearElastic3DLaw` has no `field_data` representation at all. They travel in a
typed side-channel struct instead — the `MedInfo`/`ExodusInfo` pattern:

```cpp
meshioplusplus::MdpaInfo info;
meshioplusplus::Mesh mesh = meshioplusplus::read_mdpa("model.mdpa", info);
// ... operate on mesh ...
meshioplusplus::write_mdpa("out.mdpa", mesh, info);   // properties and names restored
```

`read_mdpa(path)` and `write_mdpa(path, mesh)` are unchanged. Without an
`MdpaInfo` the properties body is parsed and dropped with a warning, which is
what the registry — and therefore the C API, Fortran, Julia, R, WASM and the
native CLI — does. That is a documented flat-ABI gap, not a silent loss.

Values are typed where they can be and verbatim where they cannot:

| In the file | In `PropertyValue` |
|---|---|
| `DENSITY 7850.0` | `mValues`, Float64 `{1}` |
| `Begin Table 4 T E ... End Table` | `mValues`, Float64 `(n, k)`, `mIsTable`, `mKey` = the header arguments |
| `CONSTITUTIVE_LAW LinearElastic3DLaw` | `mText`, verbatim |
| `LOCAL_AXES [3] (1.0, 0.0, 0.0)` | `mText`, verbatim |

The text fallback is not a gap to close later: it is what makes an unrecognized
value **lossless**, since it is re-emitted byte for byte, and it is what the
pure-Python reference does (`float()` fails, keep the raw string).

### Entity names

`MdpaInfo::mEntityNames` holds one `{mName, mIsCondition}` per cell block. It
matters because the derivation is lossy in one direction only:
`SmallDisplacementElement3D4N` resolves to `tetra` through the longest-suffix
fallback, but `tetra` only ever derives back to the canonical `Element3D4N`.
Without the info, a round trip renames every application element.

The reader's block-splitting key includes the entity name, so two adjacent
`SmallDisplacementElement3D4N` and `TotalLagrangianElement3D4N` blocks — both
`tetra` — stay separate rather than collapsing onto one name.

### `Properties` declarations on write

The writer emitted a single hard-coded `Begin Properties 0` while the entity rows
wrote their `gmsh:physical` value as the property id, so a tagged mesh produced a
file referencing undeclared properties, which Kratos's own `ModelPartIO` rejects.
Without an `MdpaInfo` the writer now emits one **empty** block per distinct id
the rows actually reference, ascending. A mesh whose ids are all 0 — every mesh
with no `gmsh:physical` — still emits exactly the same two lines as before.

### `--lenient`

`ReadOptions::mLenient` (`--lenient` on the native CLI) downgrades the remaining
rejections — `Table`, `Geometries`, `Mesh`, `Constraints`, non-empty
`SubModelPartData`/`Tables`/`Geometries`/`Constraints`, a non-numeric
`ModelPartData` value — to a warning plus a skip, recorded in
`MdpaInfo::mSkippedConstructs`. What still throws, even under `mLenient`:
non-sequential node ids, a malformed row, an unknown entity name, and
connectivity naming a node that does not exist — skipping any of those would
return a mesh that is quietly wrong rather than merely incomplete.

`meshioplusplus.mdpa.read` remains the pure-Python reference reader and is
unaffected by all of the above; it already carries this content in
`mesh.misc_data` and `field_data["properties_<id>"]`.

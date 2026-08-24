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
- **`Nodes`** — rows of either `id x y z` or bare `x y z` (auto-detected by column count). Ids may be **arbitrary** — gapped and non-monotonic both read, which is what a real Kratos deck left by a SubModelPart extraction or an entity removal looks like — and a bare row takes its **position** as its id. Points come back in **file order**, never sorted by id; connectivity, `NodalData`, `SubModelPartNodes` and `MeshNodes` all resolve through the resulting file-id → row map. A **duplicate id is a `ReadError`** (two coordinate rows claiming one id is unrepresentable, not merely incomplete), as is connectivity naming an id the block does not define.
- **`Elements <KratosType>`** / **`Conditions <KratosType>`** — the header's Kratos type name is matched by exact match first, then by longest-substring match (to avoid e.g. `"Line"` ambiguously matching inside `"Line3D2"`); each row is `id property_id node_ids...`. If the type can't be resolved from the header at all, it's inferred purely from node count. Property ids become `gmsh:physical`/`gmsh:geometrical`-style tags (MDPA reuses gmsh's tag-key convention here).
- **`Geometries <Type>`** — like Elements/Conditions but with **no property id column**; stored separately as `mesh.geometries_block` (**not** part of `mesh.cells`), a non-standard mesh-level attribute.
- **`Table <id> <var1> <var2> ...`** — rows until `End Table`; a malformed header (too few parts, non-integer id, no variables) is warned-and- skipped; stored as `field_data[f"table_{id}"] = {"variables": [...], "data": ndarray}`.
- **`Properties <id>`** — key/value pairs (auto-typed float → int-if-integer → else string) plus any nested `Table` blocks, stored under `field_data[f"properties_{id}"]`.
- **`NodalData <VAR[n]>`** / **`ElementalData`/`ConditionalData <VAR>`** — values per entity; missing entities are densified with `NaN`. An optional leading "fixed" flag column (`0`/`1`) is heuristically detected (only treated as a flag if the value is exactly 0/1 **and** more numeric values follow on the same row) — if any fixed-status is seen, a parallel `{var}_fixed_status` array is produced (sentinel `-1` = "not specified"). Scalar (0-component) variables are treated as boolean-by-presence: listed ids get `1`, unlisted get `0`.
- **`SubModelPart <Name>`** (nestable, joined with `/` for a hierarchical key e.g. `"Outer/Inner"`) — sub-blocks `SubModelPartData`, `SubModelPartTables`, `SubModelPartNodes` (0-based rows, resolved through the node-id map; an id the file never defined is dropped), `SubModelPartElements`/`Conditions` (raw 1-based ids kept **unconverted**, explicitly to preserve exact round-trip values). Not implemented: `SubModelPartGeometries`, `Constraints` sub-blocks.
- **`Mesh <id> [name]`** — an alternate/coarser mesh representation; `id=0` is invalid per Kratos convention and skipped with a warning. Sub-blocks: `MeshData`, `MeshNodes` (0-based rows, resolved through the node-id map), `MeshElements`/`Conditions` (raw 1-based ids kept as-is).

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
- `SubModelPartElements`/`Conditions` and `MeshElements`/`Conditions` (`Begin Mesh`, not `SubModelPart`) store **raw, unconverted 1-based ids** internally (in `misc_data`), but the *written* reference is resolved through `reader_element_ids_info`/`reader_condition_ids_info` for `SubModelPart` (see below) — `Begin Mesh` blocks are the one place that still writes the raw stored id verbatim (a deliberate, narrower exception; see [Original ids preserved on write](#original-ids-preserved-on-write-v9-14-0)). The *node* lists resolve through the node-id map either way.
- **Original ids are preserved on write since v9.14.0** — see [Original ids preserved on write](#original-ids-preserved-on-write-v9-14-0) below.
- The reference reader resolves connectivity **as it goes**, so a *gapped* deck whose `Elements`/`Conditions`/`Geometries` block precedes its `Nodes` block falls back to "row = id − 1" and warns; the C++ reader defers resolution to the end and is order-independent. Real decks always put `Nodes` first.
- Malformed rows in almost every block type (bad `Table` headers, data-row/variable-count mismatches, `SubModelPartNodes` entries naming an undefined node) are warned-and-skipped rather than raising — MDPA parsing is deliberately lenient/best-effort given how varied real Kratos input decks are.
- Writing `ElementalData`/`ConditionalData` omits any entity that was entirely `NaN` (never had data) rather than writing `NaN` literally.

## C++ core

`meshioplusplus::read_mdpa` / `write_mdpa` (`src/cpp/src/formats/mdpa.cpp`) implement the mesh-level part of the format against the uniform mesh API, and are registered in the shared dispatch registry — so `.mdpa` now works from the [C API](../c_api.md), [Fortran](../fortran.md), [Julia](../julia.md), [R](../r.md), [WebAssembly](../wasm.md) and the native CLI. The Kratos entity-name tables are the ones already shared with the [KRATOS mesh backend](../cpp_backends.md) (`backends/kratos_names.hpp`), extended with a longest-suffix fallback so application-specific names such as `SmallDisplacementElement3D4N` resolve through their `Element3D4N` suffix.

What the C++ core maps:

| MDPA | C++ `Mesh` |
|---|---|
| `Nodes` | `points` (always 3 columns), in file order, with a file-id → row map for arbitrary ids; original ids ride along as `point_data["mdpa:id"]` when non-trivial |
| `Elements` / `Conditions` | cell blocks, in file order, plus Int64 `cell_data["gmsh:physical"]` (the property id) and `cell_data["mdpa:id"]` (the original entity id, when non-trivial) |
| `ModelPartData` | one-element Float64 `field_data` entries (numeric values only) |
| `NodalData` | `point_data` (+ `"<VAR>_fixed_status"`) |
| `ElementalData` / `ConditionalData` | `cell_data`, **one array per cell block** — the repo-wide convention, *not* the reference reader's nested-by-cell-type layout |
| `SubModelPart` | a `Point` and/or `Cell` [named region](../regions.md); nested parts flatten to a `parent/child` name |

Everything the C++ `Mesh` cannot hold makes the reader **throw `ReadError` naming the construct** rather than dropping it silently: `Table`, `Geometries`, `Mesh <id>` and `Constraints` blocks, a non-empty `Properties` body, a non-numeric `ModelPartData` value, non-empty `SubModelPartData`/`SubModelPartTables`, a duplicate node id, connectivity naming a node the `Nodes` block does not define, and any unrecognized block. The writer emits `ModelPartData`, an empty `Properties 0`, `Nodes`, `Elements`/`Conditions`, the `*Data` blocks and one `SubModelPart` per named region; `Side` regions are dropped with a warning (MDPA has no facet-set concept). Node/element/condition ids default to `row + 1`/independent 1-based counters **unless the mesh carries `point_data`/`cell_data["mdpa:id"]`**, in which case the original ids are written back — see [Original ids preserved on write](#original-ids-preserved-on-write-v9-14-0).

**The Python `meshioplusplus.mdpa.read` deliberately does not use it.** Only the reference reader produces `mesh.misc_data`, `mesh.geometries_block` and the nested-by-cell-type `cell_data` this page documents, so preferring the C++ reader would silently change the Python API's output; reach it explicitly with `meshioplusplus._core.mdpa_read(path)` when you want the standard layout. `meshioplusplus.mdpa.write` *does* use the C++ writer, for real file paths and meshes that carry none of the MDPA extras (`misc_data`, `geometries_block`, `field_data`), falling back to the reference writer otherwise.

## Notes

- `tests/python/meshes/mdpa/test_small_cube.mdpa` — a small unit-cube tet mesh.
- `tests/python/meshes/mdpa/test_submodelpart.mdpa` — a 2D quad mesh with a nested `SubModelPart` containing Nodes/Elements/Conditions/empty-Geometries/ empty-Constraints sub-blocks.
- Additional `tests/python/input/mdpa/test_*.mdpa` fixtures target node-order permutation edge cases, minimal/degenerate geometries, hierarchical SubModelParts, and varied Table layouts.
- `tests/python/test_mdpa.py` also builds many MDPA snippets **inline** (not as files) covering nearly every block type, including a `test_roundtrip_all_blocks` exercising almost all of them at once with a NaN-aware comparison helper.

## Properties, entity names and lenient reads (v9.1.0)

The C++ reader used to throw on a non-empty `Begin Properties` body — which essentially every production deck has — so the C++ path handled geometry-only files. Three changes fix that.

### `MdpaInfo`: the side channel

Properties bodies and per-block Kratos entity names have no place on the `Mesh`: `NDArray` has ten numeric dtypes and no string one, so `CONSTITUTIVE_LAW LinearElastic3DLaw` has no `field_data` representation at all. They travel in a typed side-channel struct instead — the `MedInfo`/`ExodusInfo` pattern:

```cpp
meshioplusplus::MdpaInfo info;
meshioplusplus::Mesh mesh = meshioplusplus::read_mdpa("model.mdpa", info);
// ... operate on mesh ...
meshioplusplus::write_mdpa("out.mdpa", mesh, info);   // properties and names restored
```

`read_mdpa(path)` and `write_mdpa(path, mesh)` are unchanged.

**Since v9.2.0 you rarely need an `MdpaInfo` for the properties.** The bodies are carried on the `Mesh` itself, through the uniform API's `AddPropertySet` / `GetPropertySet` (see [the C++ API](../cpp_api.md)), so a plain

```cpp
meshioplusplus::Mesh mesh = meshioplusplus::registry_read("model.mdpa", "mdpa", {});
meshioplusplus::write_mdpa("out.mdpa", mesh);      // material data preserved
```

round-trips them — and under the KRATOS backend `mesh.GetModelPart()` hands back `Properties` blocks with real values rather than bare ids, which is what makes `to_model_part`'s "apply property" overload transfer anything. Through v9.1.0 the bodies rode this side channel *only*, and nothing reachable from `registry_readers()` could ask for one, so every registry-based consumer — which is all of them — got the ids and no material data.

`MdpaInfo` keeps two jobs the mesh channel deliberately does not do: it preserves the **file order** of the blocks (the mesh canonicalizes them to ascending id), and it carries `mEntityNames` and `mSkippedConstructs`. When both are supplied, the `MdpaInfo` wins.

Property sets do not cross a boundary where the mesh is materialized in the host language — the Python numpy `Mesh`, WASM's JS objects, the C API's `mio_mesh` accessors — since a `PropertySet` has no numpy or embind analogue. File-to-file paths (`mio_read`/`mio_write`, the WASM `convert`, the native CLI) keep them, because the `Mesh` never leaves the core.

Values are typed where they can be and verbatim where they cannot:

| In the file | In `PropertyValue` |
|---|---|
| `DENSITY 7850.0` | `mValues`, Float64 `{1}` |
| `Begin Table 4 T E ... End Table` | `mValues`, Float64 `(n, k)`, `mIsTable`, `mKey` = the header arguments |
| `CONSTITUTIVE_LAW LinearElastic3DLaw` | `mText`, verbatim |
| `LOCAL_AXES [3] (1.0, 0.0, 0.0)` | `mText`, verbatim |

The text fallback is not a gap to close later: it is what makes an unrecognized value **lossless**, since it is re-emitted byte for byte, and it is what the pure-Python reference does (`float()` fails, keep the raw string).

### Entity names

`MdpaInfo::mEntityNames` holds one `{mName, mIsCondition}` per cell block. It matters because the derivation is lossy in one direction only: `SmallDisplacementElement3D4N` resolves to `tetra` through the longest-suffix fallback, but `tetra` only ever derives back to the canonical `Element3D4N`. Without the info, a round trip renames every application element.

The reader's block-splitting key includes the entity name, so two adjacent `SmallDisplacementElement3D4N` and `TotalLagrangianElement3D4N` blocks — both `tetra` — stay separate rather than collapsing onto one name.

### `Properties` declarations on write

The writer emitted a single hard-coded `Begin Properties 0` while the entity rows wrote their `gmsh:physical` value as the property id, so a tagged mesh produced a file referencing undeclared properties, which Kratos's own `ModelPartIO` rejects. Without an `MdpaInfo` the writer emits the mesh's own property sets when it has any (v9.2.0), plus one **empty** block per referenced id no set covers; with none at all, one empty block per distinct id the rows reference, ascending. A mesh whose ids are all 0 — every mesh with no `gmsh:physical` — still emits exactly the same two lines as before.

### `--lenient`

`ReadOptions::mLenient` (`--lenient` on the native CLI) downgrades the remaining rejections — `Table`, `Geometries`, `Mesh`, `Constraints`, non-empty `SubModelPartData`/`Tables`/`Geometries`/`Constraints`, a non-numeric `ModelPartData` value — to a warning plus a skip, recorded in `MdpaInfo::mSkippedConstructs`. What still throws, even under `mLenient`: a duplicate node id, a malformed row, an unknown entity name, and connectivity naming a node that does not exist — skipping any of those would return a mesh that is quietly wrong rather than merely incomplete.

Note that **arbitrary node ids are not on that list and need no `mLenient`**: accepting them is strictly more *correct*, not more lenient, so a plain strict read handles a gapped deck. Both readers agree on one exactly — points in file order, connectivity resolved to the same rows, `point_data` keyed by the real file id — which `tests/python/test_mdpa.py::test_cpp_and_python_agree_on_gapped_ids` pins against a deck kept textually identical to the gtest suite's. The one remaining divergence is a *mixed* `Nodes` block (some rows with ids, some without): the C++ reader accepts it under the "a bare row takes its position" rule, while the reference reader's `np.loadtxt` is rectangular and rejects it.

`meshioplusplus.mdpa.read` remains the pure-Python reference reader and is unaffected by all of the above; it already carries this content in `mesh.misc_data` and `field_data["properties_<id>"]`.

## Original ids preserved on write (v9.14.0)

The read-side arbitrary-id support above (v9.13.0) closed the *read* half of roadmap `doc/roadmap.md`'s MDPA section; v9.14.0 closes the *write* half: node/element/condition ids read from a gapped or non-monotonic deck now survive a re-write, instead of both writers unconditionally renumbering to `1..n`.

**The carrier is ordinary data, not a new `Mesh` slot.** The uniform mesh API has no id-translation layer — `Mesh::Points()`/`Conn()` are dense 0-based arrays where "point index `i`" *is* row `i` — so there is nowhere on the `Mesh` itself to remember a file's original numbering. Both readers instead attach it as:

- `point_data["mdpa:id"]` — Int64, one entry per point, in read (row) order.
- `cell_data["mdpa:id"]` — Int64, one array per cell block (C++) / one array per meshio cell type (the Python reference's existing nesting), the original element/condition id.

**Attached only when it matters.** Node ids get the array exactly when `node_ids_dense` (the read-side lazy map) ever flipped to `false` — i.e. the file's ids were not already `1..n`. Elements and conditions get it exactly when either kind's own independent 1-based file-order counter ever disagreed with a row's actual id. A sequential (or id-less) deck therefore picks up no `mdpa:id` at all, and a re-write of it takes the *exact* old code path — byte-identical output, not merely equivalent. This is the same "only when it matters" contract the read-side feature established, applied to the write side.

**Both writers honour the array when present**, falling back to the old renumbering when it is absent, the wrong length, or (for `cell_data`) missing from any one block/type — never partially honouring it, since that would risk assigning the same id to two different entities. A **duplicate value** in either array is a hard `WriteError` — two nodes (or two elements, or two conditions; elements and conditions have independent Kratos id namespaces, so an element and a condition may legitimately share a numeric id) claiming one id would silently produce an ambiguous file, which is unrepresentable rather than merely incomplete, the same class of error the reader's duplicate-id `ReadError` already covers.

**Every place a node or entity is referenced elsewhere in the file resolves through the same written id**, not a bare `row + 1`/original id — this is what actually makes the feature correct rather than half-working: connectivity, `NodalData`/`ElementalData`/`ConditionalData` row keys, and `SubModelPart` node/element/condition lists. Getting only the `Nodes` block itself right while leaving connectivity on `+ 1` was the first, wrong implementation of this feature during development — caught by `test_cpp_and_python_agree_on_gapped_ids`'s round-trip check failing with `"connectivity refers to node id 1, which the file's Nodes block does not define"` the moment a re-read was attempted, which is why every gtest/pytest case here re-reads the written file rather than only inspecting it.

**A real correctness bug surfaced and got fixed along the way.** The Python reference reader has always stored `SubModelPartElements`/`Conditions` membership as the *original*, unconverted file ids (`misc_data["submodelpart_info"][name]["elements_raw"/"conditions_raw"]`). Before v9.14.0, the writer re-emitted those raw ids **verbatim** — which was already wrong whenever a plain read→write renumbered entities (the pre-v9.14.0 universal case) or reclassified one across the Elements/Conditions boundary (`_compute_blocks_name`'s dimension heuristic can turn a `Condition` into an `Element`, as happens to `tests/python/input/mdpa/test_submodelparts_hierarchical.mdpa`'s lone `LineCondition3D2N`): the emitted `SubModelPartConditions` entry could name an id that does not exist anywhere in the very file being written, or that now belongs to a different block kind — a real, reproducible corrupt-file bug that predates this feature and was invisible because the round-trip test covering that fixture only ever compared parsed-back `misc_data`, never validated that the written file's own cross-references resolved. The fix: each raw id is resolved through `misc_data["reader_element_ids_info"]`/`["reader_condition_ids_info"]` (original id → `(meshio_type, local_idx)`, captured at read time) and then through the writer's own `mdpa_written_entity_ids` (that same key → the id actually written this time, whether preserved or freshly renumbered), so the emitted reference is always real. An id that fails to resolve (only possible if the mesh was edited between read and write) is warned about and dropped rather than emitted dangling; a mesh with no reader info at all — never read via this module — falls back to writing the raw id verbatim, matching the pre-fix behaviour for that case. `Begin Mesh`'s `MeshElements`/`Conditions` sub-blocks are a narrower, **deliberately unfixed** instance of the same historical shortcut: an existing test (`test_roundtrip_all_blocks`) asserts they carry the raw ids verbatim, specifically so two independent reads' `misc_data` compare equal, so that one path keeps its documented, tested behaviour rather than being changed as a side effect of this work; `MeshNodes` (a node reference, not an entity reference) *does* resolve through the preserved node ids, since nothing tests the old `+ 1` behaviour there.

**C++ writer duplicate/fallback rules**, mirrored exactly by the Python reference writer:

- Node ids: honoured when `HasPointData("mdpa:id")` and its size matches `NumPoints()`; a duplicate value throws.
- Entity ids: honoured when `HasCellData("mdpa:id")` and `CellDataNumBlocks` equals the block count, **and** every individual block's array length matches that block's cell count (checked up front; any mismatch disables preservation for the whole write rather than partially trusting it). A duplicate value throws, checked separately within elements and within conditions.

**What is still out of scope**: geometries and `Begin Mesh` blocks are Python-reference-only, non-standard mesh attributes (`mesh.geometries_block`, `mesh.misc_data["meshes"]`); their entity-id lists (as opposed to node references, which are fixed) keep their pre-existing behaviour rather than gaining the same resolution machinery `SubModelPart` got. This closes roadmap `doc/roadmap.md`'s MDPA section in full.

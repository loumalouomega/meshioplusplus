# MED / Salome (`.med`)

The [MED](https://docs.salome-platform.org/latest/dev/MEDCoupling/developer/med-file.html) format (Salome/Code-Aster), stored in HDF5. This is the most structurally involved format meshio++ supports.

| | |
|---|---|
| **Format name** | `med` |
| **Extensions** | `.med` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | `h5py` |

## Reading & writing

```python
import meshioplusplus

mesh = meshioplusplus.read("mesh.med")
meshioplusplus.med.write("out.med", mesh, med_version="4.1.0")
```

- **`med_version`** — the `MAJ.MIN.REL` triple written to `INFOS_GENERALES` (e.g. `"4.1.0"`, `"4.0.0"`, `"3.0.0"`); default `"4.1.0"`. An unparsable string falls back to `4, 1, 0`.

`meshioplusplus.med` also exposes two standalone multi-mesh functions with no single-mesh equivalent:

```python
meshes, mesh_names = meshioplusplus.med.read_med_multi("multi.med")
meshioplusplus.med.write_med_multi("out.med", meshes, mesh_names=["fluid", "solid"])
```

- **`read_med_multi(filename, **kwargs)`** — reads every mesh under `ENS_MAA`, returning `(list[Mesh], list[str])`.
- **`write_med_multi(filename, meshes, mesh_names=None, med_version="4.1.0", **kwargs)`** — writes several meshes into one file. Missing names default to `mesh_<i>`; duplicates are de-duplicated with a numeric suffix. Field names that collide across meshes are disambiguated with an `@<mesh_name>` suffix. Forces HDF5 link-creation-order tracking for the whole file (required by medfile/Salome/mdump to read the result), restoring the previous `h5py` global config afterward.

**Note on native acceleration**: when built with `MESHIO_WITH_HDF5`, the C++ core (`meshioplusplus._core.med_read`/`med_write`) handles the **mesh-representation** part of MED exactly — points, point/cell tags, families with `GRO` group names, **named regions** derived from those families (and, since v9.8.0, from `cell_data["gmsh:physical"]` directly — see the gmsh-physical-group bridging entry below), optional `NUM` global numbering, an `INFOS_GENERALES` version check, the mesh-level metadata attributes, node-orientation permutations, and `POG`/`POG2` ragged polygons — and `meshioplusplus.med.read`/`write` use it by default, falling back to the Python/h5py implementation (as when HDF5 is absent) for the constructs the C++ path deliberately does **not** replicate byte-for-byte. Since v8.7.0 that no longer includes ordinary `CHA` **fields**: the C++ path writes/reads the **single-timestep common case** directly (plain `point_data`/`cell_data` arrays, no units, no component names, no multiple timesteps, no MED-4.1 optimization bitmask). Python is still used for: units (`med:field_units`), multi-timestep metadata (`med:step_meta`) or name-encoding (`"Name[idx] - pdt"`), non-default profiles/ELGA, and **multi-mesh** files (`read_med_multi`/`write_med_multi`, always Python). See [Quirks & limitations](#quirks-limitations).

## File structure

HDF5 groups, in write order:

```
INFOS_GENERALES                       # attrs MAJ/MIN/REL from med_version
ENS_MAA/<mesh_name>                   # mesh_name defaults to "mesh"
  (attrs DIM, ESP = points.shape[1]; REP=0; UNT/UNI = mesh.unit_time/unit_coords;
   SRT=1; NOM=<16-char-padded axis names>; DES = mesh.description or
   "Mesh created with meshio++"; TYP=0)
  -0000000000000000001-0000000000000000001    # the (single) time-step group
    (attrs CGT=1, NDT=-1, NOR=-1, PDT=-1.0)
    NOE                                # nodes
      COO                              # Fortran-order-flattened coordinates
      FAM                              # optional: per-point family/tag id
    MAI                                # mailles (cells)
      <MED type>/                      # one group per cell block, e.g. "HE8"
        NOD                            # Fortran-order 1-based connectivity
        FAM                            # optional: per-cell family/tag id
FAS/<mesh_name>
  FAMILLE_ZERO                         # attr NUM=0, always present
  NOEUD/                               # optional: point-tag family info
    FAM_<id>_<name1>_<name2>.../GRO/NOM  # 80-byte NUL-padded names
  ELEME/                               # optional: cell-tag family info, same layout
CHA/<field_name>/                      # fields (name may carry a bitmask, see below)
  (attrs MAI=mesh_name, TYP=6, NCO=n_components, NOM=<16-char-padded component names>)
  <step key>                            # one group per (NDT, NOR) time step
    NOE | NOE.<MEDTYPE> | MAI.<MEDTYPE>       # support: nodal / ELNO / ELEM
      MED_NO_PROFILE_INTERNAL/                # (or a real profile name)
        CO                                    # Fortran-order-flattened values
```

Point/cell coordinate and connectivity arrays are stored **Fortran-ordered** (column-major); the C++ core flattens/unflattens explicitly to match, since C++ has no native Fortran-order array type.

**Multi-timestep field names**: a field written at several times is stored as one `CHA` group per *base* name, with a step group per `(NDT, NOR)` pair; on read, all but the first step are surfaced as separate `point_data`/`cell_data` keys named `"{base_name}[{NDT}] - {PDT:g}"` (parsed back by `_parse_med_field_name`, a regex `(.+)\[(\d+)\]\s*-\s*([0-9.eE+-]+)$`).

## Cell types

| meshio++ | MED | meshio++ | MED |
|---|---|---|---|
| `vertex` | `PO1` | `tetra` | `TE4` |
| `line` | `SE2` | `tetra10` | `T10` |
| `line3` | `SE3` | `hexahedron` | `HE8` |
| `line4` | `SE4` | `hexahedron20` | `H20` |
| `triangle` | `TR3` | `pyramid` | `PY5` |
| `triangle6` | `TR6` | `pyramid13` | `P13` |
| `triangle7` | `TR7` | `wedge` | `PE6` |
| `quad` | `QU4` | `wedge15` | `P15` |
| `quad8` | `QU8` | `polygon` | `POG` |
| `quad9` | `QU9` | `polygon2` | `POG2` |

`polygon`/`polygon2` (both mapped to MED's `MED_POLYGON`/`MED_POLYGON2`, entity `MED_CELL`) support **ragged** cell blocks — a Voronoi-style mesh mixing 4-gons through 7-gons in one block reads back as a Python `list` of per-polygon node arrays rather than a rectangular ndarray (see [`_mesh.py`'s `CellBlock`](mdpa.md) uniform-vs-ragged detection: a block is still stored as an ndarray when every polygon in it happens to have the same vertex count).

**Node-orientation permutation** (`_med_node_perm`, linear 3D types only — applied identically on read and write, since the permutation is a fixed involution-pair):

```
tetra:      [0, 1, 3, 2]
pyramid:    [0, 3, 2, 1, 4]
wedge:      [3, 4, 5, 0, 1, 2]
hexahedron: [4, 5, 6, 7, 0, 1, 2, 3]
```

Quadratic 3D types (`tetra10`, `hexahedron20`, `pyramid13`, `wedge15`) share the same meshio++↔MED orientation difference, but no corners+midpoints permutation is implemented for them yet — they're read and written **unconverted** and may come out mis-oriented; a warning is emitted the first time one is encountered.

## Data mapping

- `point_data["point_tags"]` — per-point family/tag id.
- `cell_data["cell_tags"]` — per-cell-block family/tag id array.
- `mesh.point_tags` / `mesh.cell_tags` — mesh-level attributes (not `point_data`/`cell_data`), holding `{set_id: [subset_name, ...]}` read from `FAS/NOEUD`/`FAS/ELEME`.
- `mesh.point_tag_groups` / `mesh.cell_tag_groups` — mesh-level attributes, `{set_id: "FAM_<id>"}` short link names; always present (as a dict, possibly empty) after any Python `read()`, regardless of whether the source file had a `FAS` section at all.
- `mesh.mesh_name` / `mesh.description` / `mesh.unit_time` / `mesh.unit_coords` — mesh-level metadata attributes read from/written to `ENS_MAA`'s `NOM` (mesh group name)/`DES`/`UNT`/`UNI`. All default to `""`/`"mesh"` when absent; `description` defaults to `"Mesh created with meshio++"` on write if unset. Values round-trip through `latin-1` and are stripped of surrounding whitespace and NUL padding on read (MED files from other tools may fixed-width-pad these attributes).
- `field_data["med:nom"]` — list of component-name-lists, one per field, in field-iteration order (point_data fields, then cell_data fields).
- `field_data["med:field_units"]` / `field_data["med:step_meta"]` — dicts (not arrays) carrying per-field physical units and per-step `(NDT, NOR, PDT)` metadata; `tests/python/helpers.py::write_read`'s generic field_data comparison explicitly skips these three `med:*` keys since they aren't array-like.
- Arbitrary named point/cell data → `CHA` fields.
- **Gmsh physical-group bridging** (write-only, unconditional, native in C++ since v9.8.0): if `cell_data["gmsh:physical"]` is present, each distinct physical id becomes an element family (negative id, per MED convention) even when no `cell_tags`/`cell_sets` were set explicitly — named via `field_data` if a matching Gmsh physical-group name exists, else `f"FAM_{fid}"`/`group_{id}`; an id already covered by a named `Cell` region (e.g. one Gmsh 4.1 `$Entities` already attached on read) is not also added under the fallback name. This is what lets a Gmsh-imported mesh round-trip its physical groups through a `.med` write without the caller doing anything extra — including through the flat bindings (WASM/C API/Fortran), which have no Python to fall back to and used to hit an unconditional `WriteError` here. See also [`_pick_best_format`](gmsh.md) which prefers the gmsh writer over other `.msh`-extension candidates specifically when it detects `cell_tags`/ `point_tags`/`med:*` markers headed the other way.
- MED 4.1 **bitmask** attributes (`LEN`/`LGC`/`LNA`/`LAA`/etc., via `med/_med41.py`'s `FieldBitmaskWriter`) are written on every field to record which entity/geometry types are present across time steps, as a single 32-bit integer per attribute rather than a list of strings — required for medfile/Salome/mdump compatibility with MED ≥4.1.

## Named regions

Since v9.6.0, `FAS`/`GRO` family group names attach as `mesh.regions` — one `Point` region per `NOEUD` group name and one `Cell` region (global block-major indices) per `ELEME` group name — in addition to (not instead of) the existing `point_tags`/`cell_tags`/`point_data["point_tags"]`/`cell_data["cell_tags"]` representation, which is unaffected and still what makes a MED→MED round-trip byte-identical. This is what promotes MED out of the regions `PHASE_2` deferred table (see `doc/regions.md`); `mesh.point_sets`/`mesh.cell_sets` (the compat views over `.regions`) are populated exactly as before.

- **Read**: a family with `N` group names contributes its points/cells to all `N` regions (a multi-group family, MED's mechanism for "this entity belongs to several groups at once"). A region's `dim`/`tag` are left unspecified (`-1`) — a group name may be spread across several family ids that disagree on dimension, so there is no single honest value to report.
- **Write precedence**: when the mesh already carries `point_data["point_tags"]` (resp. `cell_data["cell_tags"]`), that native data is written exactly as before and any `Region` of the same kind is ignored. Only when the mesh carries **no** native tags of its own — the common case for a mesh converted from another format, e.g. Abaqus — are families synthesized from `Point`/`Cell` regions: one family per unique combination of region names an entity belongs to, node families numbered positive from `+1`, element families negative from `-1` (the same MED convention `_ensure_med_families` already used for the Python writer). This closes a real gap: before v9.6.0, a regions-only mesh (no native tags) written through the C++ path silently produced a file with **no** `FAS` groups at all, because nothing about the write ever raised to trigger the Python bridging. Since v9.8.0 the same cell-side synthesis **also** folds in `cell_data["gmsh:physical"]` (every Gmsh-sourced mesh's raw per-cell tag column) into the same per-cell name sets as any `Cell` region — an id resolves to a name via `field_data` (first element only, no dimension disambiguation) when one exists, else `f"group_{id}"`; an id that already resolves to an existing `Cell` region name is not added again. This is the C++ port of `_ensure_med_families`'s cell-side bridging in `_med.py`, so the two stay byte-compatible.
- **`Side` regions have no MED equivalent** (a facet is not a node or an element) and are dropped with a warning; they never prevent `Point`/`Cell` regions in the same mesh from being written.
- **The `.` in a group name is not treated specially by this reader/writer.** Kratos's own MED application (`MedApplication`) uses `.` as a SubModelPart nesting separator when *it* writes MED family names, and meshio++'s KRATOS mesh backend correspondingly warn-skips a region name containing `.` when materializing SubModelParts (`BuildSubModelPartsFromRegions`, which reserves `.` for `mdpa`'s `/`-joined nesting translated at that layer). A MED file carrying such a name round-trips its region unaffected through every other backend; only KRATOS-backend SubModelPart materialization is affected, and only for that one name.

## Global numbering

Since v9.6.0, the optional `NUM` datasets Salome/Code_Aster/Kratos write (global node/element ids) are carried as `point_data["med:num"]`/`cell_data["med:num"]`, following the repository's `<format>:<thing>` data convention rather than the `MedInfo` side channel (which the flat bindings — WASM, C API, Fortran — cannot see). Cell `NUM` is carried only when **every** cell block has it; a partial NUM array is not a mesh-wide numbering and is dropped with a warning rather than fabricated (no `iota` filler, unlike some MED consumers). A mesh with no `med:num` writes exactly as before — no `NUM` dataset appears anywhere.

## MED version check

Since v9.6.0, both readers check `INFOS_GENERALES`'s `MAJ` attribute and reject (`ReadError`) a file written by a MED major version newer than 4, with a message naming the file's actual `MAJ.MIN.REL` — instead of failing later with an unrelated structural error. Older majors are unaffected (the repository's own `cylinder.med` fixture is MED 3.0.0 and keeps reading); a file with no `INFOS_GENERALES` at all is also unaffected. The check is deliberately duplicated in the C++ reader and the Python fallback with the same message, so a file rejected by one is rejected by the other rather than silently misread by whichever path is tried second.

## Lenient reads

Since v9.9.0 the C++ reader accepts `ReadOptions` (MED is in `registry_readers_ex()`, which is what makes the options reach WASM, the C API, Fortran, Julia, R and both CLIs with no per-binding code), and `ReadOptions::mLenient` — the mechanism `mdpa` established in v9.1.0 — opens the whole set of constructs that previously threw `"handled by Python fallback"`.

**Strict reads are unchanged.** Every construct below still throws by default, so `meshioplusplus.med.read` still falls back to the pure-Python reference and the Python surface is byte-identical to v9.8.0's. What changes is that a Python-less binding can now get through a real Salome/Code_Aster file at all, instead of failing on sight of it.

Under `mLenient`, constructs that can be *described* are read into `MedInfo` rather than merely skipped:

| `MedInfo` field | carries |
| --- | --- |
| `mFieldUnits` | each field's `(UNI, UNT)` where non-empty |
| `mStepMeta` | each field's `(NDT, NOR, PDT)` where not the write-side default |
| `mFieldTimeValues` | every step's `PDT`, in step order |
| `mSkippedConstructs` | one entry per construct that was skipped, for diagnostics |

Constructs that have no representation at all drop **that one field** with a `log::warn` and keep the rest of the file: a named `PFL` profile, an `ELNO`/`ELGA` support, and a field mixing nodal and cell supports. `ELNO`/`ELGA` is *structurally* impossible rather than merely unimplemented — the uniform mesh API's `cell_data` is always `(n,)` or `(n,k)`, never a per-node-within-cell 3-D shape — so no amount of work in this format closes it.

`MedInfo` is a side channel the registry drops (`ReadFn` has no info slot), so the flat bindings get the lenient *read* but not the recorded metadata — the same documented gap `ExodusInfo` and `GmshInfo` have.

### Selecting a time step

`ReadOptions::mTimeStep` selects one step of a multi-step `CHA` field (0-based, negative counts from the end — the `ResolveTimeStep` contract, see [selective reads](../selective_read.md#reading-one-time-step)). A multi-step field used to fail the read outright.

A non-default step is honoured **without** `mLenient`, deliberately: it is a request the Python shim never makes, so no Python behaviour depends on it. Step order comes from the `(NDT, NOR)` subgroup names, which MED zero-pads, so name order *is* step order.

## Quirks & limitations

- **Two supports for cell data**: `ELEM` (one value per cell, exactly 1 Gauss point) and `ELNO` (one value per node-per-cell, "defined at every node"); which one is used is decided by shape (`ndim <= 2` → ELEM, `shape[1] == num_nodes_per_cell[type]` → ELNO, else `ELGA`). **`ELGA` (general Gauss-point data at unknown points) is silently skipped on write** — there's no representation for arbitrary Gauss-point layouts.
- **A field's `NOM` carries 16 characters per component** (v9.9.0). Both writers previously wrote a single blank 16-char slot whatever the component count, which deviates from MED's convention for any k>1 field. When no explicit component names are supplied (`field_data["med:nom"]`, a Python-only convention), MED's own default spelling `V1..Vk` is generated. A scalar field keeps the single blank 16-char slot, so its bytes are unchanged. Consequence: a k>1 field written by either writer now reads back with `med:nom` populated where it previously came back empty.
- **A field's row count must equal its entity count** (v9.9.0). `write_med` raises a named `WriteError` when a `point_data` array's rows disagree with `NumPoints()` (or a `cell_data` block's with that block's cell count), naming the array and both counts. There was previously no write-side check at all: a mis-shaped array -- most often an `(n,k)` vector flattened to `(nk,)` by a caller that lost its component count -- wrote `NBR = nk` against `n` points and produced a file this very reader rejects, so the failure surfaced far from its cause. Multi-component fields themselves round-trip normally.
- Family names longer than 80 bytes (after `latin-1` encoding) raise `WriteError` rather than silently truncating.
- A family with no groups omits the `GRO` dataset entirely (rather than writing an empty one); an all-default/no-tags mesh likewise omits `FAS` family groups it doesn't need.
- **Same-type cell blocks are consolidated on write** (since v9.8.0; before that, rejected up front with `WriteError("MED files cannot have two sections of the same cell type.")`): MED has no way to represent two `MAI` sections of one type, but MSH 4.1's canonical structure is one cell block per *entity*, so real Gmsh-4.1-sourced meshes routinely carry several same-type blocks. `write_med` now groups blocks by type (first-seen order) and writes one section per type, concatenating the contributing blocks' connectivity/`FAM`/`NUM` — mirroring the Python reference's own write-time merge (`_med.py`'s `cells_by_type`), so this was never a Python-vs-C++ divergence, only a C++-vs-flat-bindings one. Blocks of the same type must still agree on node count to be merged; a disagreement is a `WriteError` naming the type. `FAM`/`NUM` are written for a consolidated section only when the source data (native `cell_tags`/`med:num`, or the region-synthesized families above) covers *every* contributing block — a partial array is dropped with a `log::warn` rather than partially written.
- Re-writing a field under an already-used name appends a new support group under that field's most recent timestep, rather than creating a distinct field.
- `FAS` (the families group) may live either under the mesh's own time-step group or at the top level (`f["FAS"][mesh_name]`) — both readers check the nested location first, then fall back to top-level.
- **C++ vs Python split (default path):** the C++ core handles points, point/ cell tags, families (with `GRO` group names) — additionally attached as named regions, see [Named regions](#named-regions) — optional `NUM` global numbering (see [Global numbering](#global-numbering)), an `INFOS_GENERALES` version check, the mesh-level metadata attributes (`mesh_name`/`description`/`unit_time`/`unit_coords`/ `point_tag_groups`/`cell_tag_groups`), the node-orientation permutations, `POG`/`POG2` ragged polygons, and (since v8.7.0) **ordinary `CHA` fields — the single-timestep common case only** (plain `point_data`/`cell_data` arrays, one `NOE` or one-or-more `MAI.<type>` support subgroups, fixed `ndt=1`/`nor=-1`/`pdt=0.0`, blank `UNI`/`UNT`/`NOM`, and **no MED-4.1 optimization bitmask** — a deliberate scope cut: our own reader never reads the bitmask, so its absence costs nothing for a meshio++ round-trip, only for interop with tools like Salome/MEDCoupling that use it). Field values are always resolved by **cell type** (`meshio_to_med()`), never a fixed block index — MED reorders `MAI` blocks alphabetically by type code on read, so block *position* does not survive a round-trip even though block *identity* does. The C++ reader **declines** (throws, deferring the whole file to Python) rather than silently drop information whenever a field's `UNI`/`UNT` is non-empty or a timestep's `NDT`/`NOR`/`PDT` isn't the write-side default — this is checked against the file's actual attribute values, not the group's name, since a hand-edited file can carry meaningful metadata under the same default-looking timestep key. The **enhanced-feature guard lives in the Python shim** (`meshioplusplus/med/__init__.py`), not the C++ core: `med:field_units`/`med:step_meta` are Python-only dict-valued `field_data` conventions that cannot survive the Python→C++ mesh conversion at all (`med_write`'s pybind wrapper always uses lenient field_data coercion, silently dropping any non-numeric entry before the C++ `Mesh` ever sees it), so a C++-side check for them would be dead code under every binding; the shim also defers when an array name uses the `"Name[idx] - pdt"` multi-timestep encoding, which the C++ writer has no notion of. It still **raises by default** (so `meshioplusplus.med.read`/`write` fall back to Python) for: any of those three enhanced-field cases, non-default **profiles**/ELGA (ELNO/ELGA cell data is structurally unreachable from the C++ core regardless — the uniform mesh API's `cell_data` is always `(n,)`/`(n,k)`, never a 3-D per-node-within-cell shape), and **multi-mesh** files. Since v9.9.0 everything on that list except multi-mesh is reachable under `ReadOptions::mLenient` (see [Lenient reads](#lenient-reads)), which is what a Python-less binding needs; the strict default is unchanged precisely so the shim keeps deferring and the Python surface does not move. `read_med_multi`/`write_med_multi` are always Python. The `gmsh:physical`→family bridging and same-type block consolidation used to be on this deferred-to-Python list too; both are native in the C++ core since v9.8.0 (see the two entries above) and no longer force a fallback — a `.msh`→`.med` conversion through the shared registry (WASM/C API/Fortran, which have no Python fallback at all) now succeeds for real Gmsh-4.1 output.
- **A cell block with no `FAM` array reads as family 0** (v9.9.0), instead of failing the whole file. MED spells "belongs to no family" as id **0**, so this is the file's own meaning rather than a guess, and a cell_data array must cover every block for the uniform mesh API to hold it at all. The Python reference was **also** wrong here, differently: it appended nothing for such a block, leaving `cell_tags` *shorter* than `mesh.cells`. Both now zero-fill identically.
- Ragged `polygon`/`polygon2` blocks (mixed vertex counts) round-trip through the C++ core as `POG`/`POG2` (CSR `NOD` + `INN` offset arrays); they cross the C++↔Python boundary as a copied list of arrays (ragged data cannot be zero-copy). `polyhedron*` blocks round-trip as `MED_POLYHEDRON` (`POE`) since v9.19.0, which needs **three** 1-based arrays where a polygon needs two: `NOD` (every face's nodes, flat), `INN` (face → start in `NOD`) and `IND` (cell → start in the face list). MED holds one section per type inside a `MAI` group, so every `polyhedron<N>` block is canonicalised to a single `POE` on write and regrouped by unique node count on read — the same bucketing the OpenFOAM, EnSight and CGNS readers use. (Before v9.19.0 this line claimed polyhedra were "Python-only for MED"; that was simply wrong — **neither** path had a `POE` entry.)

## Notes

- `tests/python/meshes/med/box.med` (Code_Aster 13.6) — single hexahedron, 8 points (sum 12), point_data `resu____DEPL` (displacement, shape `(8,3)`), cell_data `resu____EPSI_ELNO`/`resu____SIEF_ELNO` (ELNO strain/stress, shape `(1,8,6)`), `resu____ENEL_ELNO`/`resu____ENEL_ELEM` (energy, both supports).
- `tests/python/meshes/med/cylinder.med` (Salome 9.2.2, version downgraded to 3.0.0) — mixed cell types `{pyramid:18, quad:18, line:17, tetra:63, triangle:4}`, point tags summing to 52 with named families like `{2:["Side"], 3:["Side","Top"], 4:["Top"]}`, cell tags e.g. `{-6:["Top circle"], -9:["A","B"], ...}`.
- `tests/python/meshes/med/input_code_aster.med` (~4.8 MB) and `tests/python/meshes/med/voronoi_hex.med` (~15 KB, ragged Voronoi polygons) — larger fixtures covering the multi-mesh/polygon/metadata read paths above.
- Originally ported from upstream meshio; the multi-mesh, ragged-polygon, MED-4.1 bitmask, node-orientation, and gmsh-family-bridging enhancements were contributed by [Simvia's meshlane fork](https://github.com/simvia-tech/meshlane) and brought back into this repository (see `CHANGELOG.md`).

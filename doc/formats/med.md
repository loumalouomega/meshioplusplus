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

**Note on native acceleration**: when built with `MESHIO_WITH_HDF5`, the C++ core (`meshioplusplus._core.med_read`/`med_write`) handles the **mesh-representation** part of MED exactly — points, point/cell tags, families with `GRO` group names, the mesh-level metadata attributes, node-orientation permutations, and `POG`/`POG2` ragged polygons — and `meshioplusplus.med.read`/`write` use it by default, falling back to the Python/h5py implementation (as when HDF5 is absent) for the constructs the C++ path deliberately does **not** replicate byte-for-byte. Since v8.7.0 that no longer includes ordinary `CHA` **fields**: the C++ path writes/reads the **single-timestep common case** directly (plain `point_data`/`cell_data` arrays, no units, no component names, no multiple timesteps, no MED-4.1 optimization bitmask). Python is still used for: units (`med:field_units`), multi-timestep metadata (`med:step_meta`) or name-encoding (`"Name[idx] - pdt"`), the `gmsh:physical`→family bridging, non-default profiles/ELGA, and **multi-mesh** files (`read_med_multi`/`write_med_multi`, always Python). See [Quirks & limitations](#quirks-limitations).

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
| `triangle` | `TR3` | `hexahedron20` | `H20` |
| `triangle6` | `TR6` | `pyramid` | `PY5` |
| `triangle7` | `TR7` | `pyramid13` | `P13` |
| `quad` | `QU4` | `wedge` | `PE6` |
| `quad8` | `QU8` | `wedge15` | `P15` |
| `quad9` | `QU9` | `polygon` | `POG` |
| | | `polygon2` | `POG2` |

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
- **Gmsh physical-group bridging** (write-only, unconditional): if `cell_data["gmsh:physical"]` is present, each distinct physical id becomes an element family (negative id, per MED convention) even when no `cell_tags`/`cell_sets` were set explicitly — named via `field_data` if a matching Gmsh physical-group name exists, else `f"FAM_{fid}"`/`group_{id}`. This is what lets a Gmsh-imported mesh round-trip its physical groups through a `.med` write without the caller doing anything extra; see also [`_pick_best_format`](gmsh.md) which prefers the gmsh writer over other `.msh`-extension candidates specifically when it detects `cell_tags`/ `point_tags`/`med:*` markers headed the other way.
- MED 4.1 **bitmask** attributes (`LEN`/`LGC`/`LNA`/`LAA`/etc., via `med/_med41.py`'s `FieldBitmaskWriter`) are written on every field to record which entity/geometry types are present across time steps, as a single 32-bit integer per attribute rather than a list of strings — required for medfile/Salome/mdump compatibility with MED ≥4.1.

## Quirks & limitations

- **Two supports for cell data**: `ELEM` (one value per cell, exactly 1 Gauss point) and `ELNO` (one value per node-per-cell, "defined at every node"); which one is used is decided by shape (`ndim <= 2` → ELEM, `shape[1] == num_nodes_per_cell[type]` → ELNO, else `ELGA`). **`ELGA` (general Gauss-point data at unknown points) is silently skipped on write** — there's no representation for arbitrary Gauss-point layouts.
- Family names longer than 80 bytes (after `latin-1` encoding) raise `WriteError` rather than silently truncating.
- A family with no groups omits the `GRO` dataset entirely (rather than writing an empty one); an all-default/no-tags mesh likewise omits `FAS` family groups it doesn't need.
- Writing a mesh with two cell blocks of the same type is rejected up front (`WriteError`) — MED cannot represent two blocks of one type.
- Re-writing a field under an already-used name appends a new support group under that field's most recent timestep, rather than creating a distinct field.
- `FAS` (the families group) may live either under the mesh's own time-step group or at the top level (`f["FAS"][mesh_name]`) — both readers check the nested location first, then fall back to top-level.
- **C++ vs Python split (default path):** the C++ core handles points, point/ cell tags, families (with `GRO` group names), the mesh-level metadata attributes (`mesh_name`/`description`/`unit_time`/`unit_coords`/ `point_tag_groups`/`cell_tag_groups`), the node-orientation permutations, `POG`/`POG2` ragged polygons, and (since v8.7.0) **ordinary `CHA` fields — the single-timestep common case only** (plain `point_data`/`cell_data` arrays, one `NOE` or one-or-more `MAI.<type>` support subgroups, fixed `ndt=1`/`nor=-1`/`pdt=0.0`, blank `UNI`/`UNT`/`NOM`, and **no MED-4.1 optimization bitmask** — a deliberate scope cut: our own reader never reads the bitmask, so its absence costs nothing for a meshio++ round-trip, only for interop with tools like Salome/MEDCoupling that use it). Field values are always resolved by **cell type** (`meshio_to_med()`), never a fixed block index — MED reorders `MAI` blocks alphabetically by type code on read, so block *position* does not survive a round-trip even though block *identity* does. The C++ reader **declines** (throws, deferring the whole file to Python) rather than silently drop information whenever a field's `UNI`/`UNT` is non-empty or a timestep's `NDT`/`NOR`/`PDT` isn't the write-side default — this is checked against the file's actual attribute values, not the group's name, since a hand-edited file can carry meaningful metadata under the same default-looking timestep key. The **enhanced-feature guard lives in the Python shim** (`meshioplusplus/med/__init__.py`), not the C++ core: `med:field_units`/`med:step_meta` are Python-only dict-valued `field_data` conventions that cannot survive the Python→C++ mesh conversion at all (`med_write`'s pybind wrapper always uses lenient field_data coercion, silently dropping any non-numeric entry before the C++ `Mesh` ever sees it), so a C++-side check for them would be dead code under every binding; the shim also defers when an array name uses the `"Name[idx] - pdt"` multi-timestep encoding, which the C++ writer has no notion of. It still **raises** (so `meshioplusplus.med.read`/`write` fall back to Python) for: any of those three enhanced-field cases, the `gmsh:physical`→family **bridging** on write, non-default **profiles**/ELGA (ELNO/ELGA cell data is structurally unreachable from the C++ core regardless — the uniform mesh API's `cell_data` is always `(n,)`/`(n,k)`, never a 3-D per-node-within-cell shape), and **multi-mesh** files. `read_med_multi`/`write_med_multi` are always Python.
- Ragged `polygon`/`polygon2` blocks (mixed vertex counts) round-trip through the C++ core as `POG`/`POG2` (CSR `NOD` + `INN` offset arrays); they cross the C++↔Python boundary as a copied list of arrays (ragged data cannot be zero-copy). `polyhedron*` blocks remain Python-only for MED.

## Notes

- `tests/python/meshes/med/box.med` (Code_Aster 13.6) — single hexahedron, 8 points (sum 12), point_data `resu____DEPL` (displacement, shape `(8,3)`), cell_data `resu____EPSI_ELNO`/`resu____SIEF_ELNO` (ELNO strain/stress, shape `(1,8,6)`), `resu____ENEL_ELNO`/`resu____ENEL_ELEM` (energy, both supports).
- `tests/python/meshes/med/cylinder.med` (Salome 9.2.2, version downgraded to 3.0.0) — mixed cell types `{pyramid:18, quad:18, line:17, tetra:63, triangle:4}`, point tags summing to 52 with named families like `{2:["Side"], 3:["Side","Top"], 4:["Top"]}`, cell tags e.g. `{-6:["Top circle"], -9:["A","B"], ...}`.
- `tests/python/meshes/med/input_code_aster.med` (~4.8 MB) and `tests/python/meshes/med/voronoi_hex.med` (~15 KB, ragged Voronoi polygons) — larger fixtures covering the multi-mesh/polygon/metadata read paths above.
- Originally ported from upstream meshio; the multi-mesh, ragged-polygon, MED-4.1 bitmask, node-orientation, and gmsh-family-bridging enhancements were contributed by [Simvia's meshlane fork](https://github.com/simvia-tech/meshlane) and brought back into this repository (see `CHANGELOG.md`).

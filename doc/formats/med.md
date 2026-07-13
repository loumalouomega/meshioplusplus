# MED / Salome (`.med`)

The [MED](https://docs.salome-platform.org/latest/dev/MEDCoupling/developer/med-file.html)
format (Salome/Code-Aster), stored in HDF5. This is the most structurally
involved format meshio supports.

| | |
|---|---|
| **Format name** | `med` |
| **Extensions** | `.med` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | `h5py` (or a C++ build with HDF5) |

## Reading & writing

```python
import meshio

mesh = meshio.read("mesh.med")
meshio.med.write("out.med", mesh)
```

`write` takes no keyword arguments — MED does not support compression, unlike
most of the other HDF5-backed formats.

## File structure

HDF5 groups, in write order:

```
INFOS_GENERALES                       # attrs MAJ=3, MIN=0, REL=0
ENS_MAA/<mesh_name>                   # mesh_name is always literally "mesh"
  (attrs DIM, ESP = points.shape[1]; REP=0; UNT/UNI="";
   SRT=1; NOM=<16-char-padded axis names>; DES="Mesh created with meshio"; TYP=0)
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
CHA/<field_name>/                      # fields
  (attrs MAI=mesh_name, TYP=6, NCO=n_components, NOM=<16-char-padded component names>)
  0000000000000000000100000000000000000001   # time-step group (NDT=1, NOR=1, PDT=0.0)
    NOE | NOE.<MEDTYPE> | MAI.<MEDTYPE>       # support: nodal / ELNO / ELEM
      MED_NO_PROFILE_INTERNAL/                # (or a real profile name)
        CO                                    # Fortran-order-flattened values
```

Point/cell coordinate and connectivity arrays are stored **Fortran-ordered**
(column-major); the C++ core flattens/unflattens explicitly to match, since
C++ has no native Fortran-order array type.

## Cell types

| meshio | MED | meshio | MED |
|---|---|---|---|
| `vertex` | `PO1` | `tetra` | `TE4` |
| `line` | `SE2` | `tetra10` | `T10` |
| `line3` | `SE3` | `hexahedron` | `HE8` |
| `triangle` | `TR3` | `hexahedron20` | `H20` |
| `triangle6` | `TR6` | `pyramid` | `PY5` |
| `quad` | `QU4` | `pyramid13` | `P13` |
| `quad8` | `QU8` | `wedge` | `PE6` |
| | | `wedge15` | `P15` |

## Data mapping

- `point_data["point_tags"]` — per-point family/tag id.
- `cell_data["cell_tags"]` — per-cell-block family/tag id array.
- `mesh.point_tags` / `mesh.cell_tags` — mesh-level attributes (not
  `point_data`/`cell_data`), holding `{set_id: [subset_name, ...]}` read from
  `FAS/NOEUD`/`FAS/ELEME`. In the C++ core these travel through a dedicated
  out-of-band `MedInfo` struct rather than through `Mesh` itself.
- `field_data["med:nom"]` — list of component-name-lists, one per field, in
  field-iteration order (point_data fields, then cell_data fields).
- Arbitrary named point/cell data → `CHA` fields; multi-timestep fields get a
  `"{name}[i] - {t:g}"` suffix on read.

## Quirks & limitations

- **Two supports for cell data**: `ELEM` (one value per cell, exactly 1
  Gauss point) and `ELNO` (one value per node-per-cell, "defined at every
  node"); which one is used is decided by shape (`ndim &lt;= 2` → ELEM,
  `shape[1] == num_nodes_per_cell[type]` → ELNO, else `ELGA`).
  **`ELGA` (general Gauss-point data at unknown points) is silently skipped
  on write** in both Python and C++ — there's no representation for
  arbitrary Gauss-point layouts.
- **Non-default profiles** (MED's `PROFILS` mechanism restricting field data
  to a subset of points/cells) are fully supported by the Python reader
  (missing entries filled with `NaN`), but the **C++ reader explicitly
  refuses** any file with a `PROFILS` group or a non-`MED_NO_PROFILE_INTERNAL`
  profile, forcing the Python fallback.
- **Partial cell tags/fields** (not every cell block has a `FAM` or a field
  entry) are also explicitly rejected by the C++ path, since the fixed-size
  `Mesh` representation has no way to express "this block is missing" — 
  Python handles it via per-block `None`/empty defaults.
- Family names are 80-byte, NUL-padded, decoded character-by-character then
  stripped.
- `FAS` (the families group) may live either under the mesh's own time-step
  group or at the top level (`f["FAS"][mesh_name]`) — both readers check the
  nested location first, then fall back to top-level.
- Writing a mesh with two cell blocks of the same type is rejected up front
  (`WriteError`) — MED cannot represent two blocks of one type.
- Re-writing a field under an already-used name appends a new support group
  under that field's most recent timestep, rather than creating a distinct
  field.

## Notes

- `tests/meshes/med/box.med` (Code_Aster 13.6) — single hexahedron, 8 points
  (sum 12), point_data `resu____DEPL` (displacement, shape `(8,3)`),
  cell_data `resu____EPSI_ELNO`/`resu____SIEF_ELNO` (ELNO strain/stress,
  shape `(1,8,6)`), `resu____ENEL_ELNO`/`resu____ENEL_ELEM` (energy, both
  supports).
- `tests/meshes/med/cylinder.med` (Salome 9.2.2, version downgraded to
  3.0.0) — mixed cell types `{pyramid:18, quad:18, line:17, tetra:63,
  triangle:4}`, point tags summing to 52 with named families like
  `{2:["Side"], 3:["Side","Top"], 4:["Top"]}`, cell tags e.g.
  `{-6:["Top circle"], -9:["A","B"], ...}`.
- Read/written through the C++ core when built with `MESHIO_WITH_HDF5`,
  otherwise through `h5py` — see
  [native acceleration](../formats.md#native-acceleration-and-fallbacks).

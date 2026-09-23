# I-DEAS Universal File — UNV (`.unv`, `.uff`)

The [I-DEAS Universal File](https://www.ceas3.uc.edu/sdrluff/) is the ASCII interchange format of SDRC I-DEAS / Siemens NX, written and read by Salome (SMESH), Code_Aster, gmsh, ElmerGrid, OpenFOAM's `ideasUnvToFoam` and the experimental-modal-analysis tools (test rigs, pyuff). A file is a sequence of **datasets**, each opened by a `-1` line and a dataset-number line and closed by another `-1` line. Meshes usually carry the `.unv` extension; test-lab result files the `.uff` one. Both are the same format.

| | |
|---|---|
| **Format name** | `unv` |
| **Extensions** | `.unv`, `.uff` (also recognised by content) |
| **Read / Write** | ✓ / ✓ |
| **Steps** | ✓ (read): the steps of its results, see [Results and steps](#results-and-steps) |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshioplusplus

mesh = meshioplusplus.read("model.unv")                   # mesh, groups, first step of results
last = meshioplusplus.read("modes.unv", time_step=-1)     # the last step
meshioplusplus.read_metadata("frf.uff")["time_values"]    # every step's time / frequency
for freq, m in meshioplusplus.read_sequence("frf.uff"):   # one mesh per frequency line
    ...

meshioplusplus.unv.write("out.unv", mesh)                 # nodes, elements, groups, results
meshioplusplus.unv.write("out.unv", mesh, code_aster=True)
```

`read` takes the sequence-engine keywords (`time_step`, `points_only`, `arrays`). `write` takes two optional keywords:

- `code_aster=False`: when `True`, results are written as the legacy datasets **55** (`point_data`) and **56** (`cell_data`) in single precision (`6E13.5`, as the SDRL records read), instead of dataset 2414.
- `node_dataset=2411`: the node dataset to emit; `781` is also accepted.

## Datasets

| Dataset | Read | Written | Maps to |
|---|---|---|---|
| **2411** nodes (double precision, `D` exponents), **781**, legacy **15** | ✓ | 2411 or 781 | `points` |
| **2412** elements, legacy **780** | ✓ | 2412 | cell blocks, `cell_data["unv:pid"]` / `["unv:mid"]` |
| **2467**, **2477**, **2452**, **2435** permanent groups (entity quadruples); **2417**, **2429**, **2430**, **2432** (entity pairs) | ✓ | 2467 | point and cell [regions](../regions.md) |
| **164** units | ✓ | when the mesh carries it | `field_data["unv:units"]`, `["unv:unit_factors"]` |
| **2420** coordinate systems | ✓ (applied to nodes) | — | node coordinates |
| **2414** analysis data at nodes (location 1) / on elements (location 2) | ✓ | ✓ (default) | `point_data` / `cell_data`, one step each |
| **55** data at nodes, **56** data on elements | ✓ | with `code_aster=True` | `point_data` / `cell_data`, one step each |
| **58** / **58b** function at nodal DOF (ASCII / binary) | ✓ | — | `point_data`, one step per abscissa value |

Any other dataset is skipped. Data at nodes on elements (dataset 57, 2414 location 3) and data at points (2414 location 5) are skipped with a warning.

### Nodes and coordinate systems

Record 1 of a 2411 node names its definition coordinate system. When a 2420 dataset defines that system as **Cartesian**, the coordinates are moved into the global system as `x = M·x_local + o`, where `M` is the matrix's first three rows and `o` its fourth: the convention of Salome's reader, the only widely used consumer of the dataset. A **cylindrical** or **spherical** system is warned about and its nodes are left in local coordinates; a system number above 1 with no 2420 definition is taken as global, with a warning. Units (164) are recorded, not applied: `unv:unit_factors` holds the length, force, temperature and temperature-offset factors *to SI* as the file states them. The units code is read from columns 1–10 of record 1 (`I10,20A1,I10`), since some writers run the description into it (`5mm (milli-newton)`); since v16.4.0 Python no longer refuses such a record.

### Elements

Record 1 of a 2412 element is `label, FE descriptor, physical property, material property, colour, node count` (the legacy 780 has eight fields: the property numbers are its fourth and sixth). Beam descriptors (11, 21-25) carry an extra record (orientation node and cross sections) before the node labels; it is skipped on read and written as `0 1 1`. The physical and material property numbers become the integer `cell_data` `unv:pid` and `unv:mid`, and are written back from them (`unv:mid` falls back to `unv:pid`, `unv:pid` to 1).

The FE descriptor *and the node count* select the type; an element whose node count does not match its descriptor is skipped with a warning. This matters in practice: gmsh writes its `quad9` and `hexahedron27` under the linear descriptors 94 and 115 in its own node order.

| Descriptors | meshio++ type | Written as |
|---|---|---|
| 11, 21-25 (rods and beams, 2 or 3 nodes) | `line` / `line3` | 21 / 24 |
| 41, 51, 61, 74, 81, 91 | `triangle` | 91 |
| 42, 52, 62, 72, 82, 92 | `triangle6` | 92 |
| 44, 54, 64, 71, 84, 94, 122 | `quad` | 94 |
| 45, 55, 65, 75, 85, 95 | `quad8` (9 nodes: `quad9`) | 95 |
| 111 / 118 | `tetra` / `tetra10` | 111 / 118 |
| 112 / 113 | `wedge` / `wedge15` | 112 / 113 |
| 115 / 116 | `hexahedron` / `hexahedron20` | 115 / 116 |
| 119, 312 / 114 | `pyramid` / `pyramid13` | 312 / 114 |

Linear elements share VTK's node order. Parabolic elements list their nodes ring by ring with each mid-side node **between** its two corners: a quadratic hexahedron is the bottom ring (corner, mid, corner, mid, …), then the four vertical mid-edges, then the top ring. The permutation below (`meshio[perm[i]] = unv[i]`) was pinned against gmsh, which writes the same mesh as `.unv` and `.msh`, and matches Salome's reader. It lives in the shared [node-ordering registry](../node_ordering.md) under `"unv"`:

| Type | Permutation |
|---|---|
| `line3` | `[0, 2, 1]` |
| `triangle6` | `[0, 3, 1, 4, 2, 5]` |
| `quad8` / `quad9` | `[0, 4, 1, 5, 2, 6, 3, 7]` (+ `8`) |
| `tetra10` | `[0, 4, 1, 5, 2, 6, 7, 8, 9, 3]` |
| `pyramid13` | `[0, 5, 1, 6, 2, 7, 3, 8, 9, 10, 11, 12, 4]` |
| `wedge15` | `[0, 6, 1, 7, 2, 8, 12, 13, 14, 3, 9, 4, 10, 5, 11]` |
| `hexahedron20` | `[0, 8, 1, 9, 2, 10, 3, 11, 16, 17, 18, 19, 4, 12, 5, 13, 6, 14, 7, 15]` |

### Groups

A permanent group lists entities by type: **7 is a node, 8 a finite element**; other entity types (loads, restraint sets, …) are ignored. Each group becomes up to two [regions](../regions.md) of the same name, tagged with the group number: a `point` region of its nodes and a `cell` region of its elements (`dim` = the highest dimension among them). A group with neither is kept as an empty `cell` region. Members the file names but did not define are dropped with a warning. Groups therefore also appear through `mesh.point_sets` / `mesh.cell_sets`, and every flat binding (C, Fortran, Julia, R, WASM, the native CLI) sees them through the mesh's regions.

On write, every `point` and `cell` region becomes a 2467 group; a point and a cell region sharing a name are one group. The group number is the region's `tag` when it is positive and not already taken, otherwise the next free number. UNV has no facet group, so `side` regions are dropped with a warning (and a provenance note).

The property numbers `unv:pid` / `unv:mid` stay cell data rather than regions: as regions they would be written back as extra groups.

## Results and steps

Every result block (2414, 55, 56) and every abscissa value of a 58 function belongs to a **step**. Blocks with the same step merge (a displacement and a stress of the same mode are one step); steps are listed in the order the file first names them. `time_step` picks one (`0` is the first, negative counts from the end, out of range is an error), `read_metadata(...)["time_values"]` lists their values, and `read_sequence` walks them. A plain `read` returns the first step.

| Analysis type (2414 record 9 / 55 record 6) | Step number | Step value (`meshio:time`) |
|---|---|---|
| 1 static, 0 unknown, others | load set | the record's time when non-zero, else the load set |
| 2 normal mode | mode number | frequency |
| 3, 7 complex eigenvalue | mode number | \|Im λ\| / 2π |
| 4 transient, 9 static non-linear | time-step number | time |
| 5 frequency response | frequency number | frequency |
| 6 buckling | mode number | eigenvalue |
| 58 functions | sample index (1-based) | the abscissa (analysis 5 when the abscissa is a frequency, 4 when it is a time) |

When a file has more than one step, or a non-zero analysis type, the step carries `field_data` `meshio:time` (float), `unv:analysis` and `unv:step` (integers). A single step of an unknown analysis (what older meshio++ versions wrote) reads without them, exactly as before. On write, those three field data describe the one step written (records 9, 10 and 12 of 2414, 6-8 of 55/56).

Result arrays:

- The array is named after the dataset (2414 record 2, 55/56 ID line 1). A blank or `NONE` name falls back to the result type (`displacement`, `stress`, `temperature`, …). Names repeated within one step get `_2`, `_3`.
- Entities the block does not list are **NaN**, and NaN rows are not written back.
- Symmetric tensors are stored `Sxx Sxy Syy Sxz Syz Szz` and read into meshio++'s `xx yy zz xy yz zx`; general tensors (9 values) are stored column by column and read row by row. Both are reordered back on write.
- Complex data (data types 5, 6) becomes `<name>_real` and `<name>_imag`. Integer data (type 1) reads as float64.
- 2414 values are written in double precision, one line per entity (`E20.12`); `code_aster=True` writes single precision `6E13.5`.

### Functions (58 / 58b)

A test-lab file holds one dataset 58 per measured function: an FRF, a coherence, a spectrum … at one response node and direction, over a frequency (or time) axis. Functions of the same kind (function type, load case, reference node and axis) and the same abscissa grid are one group; each abscissa value is a step, and each group gives `point_data` arrays at that step:

- `<name>` with 3 components for the translational directions ±X/±Y/±Z, `<name>_rot` for the rotations, and a scalar `<name>` (or `<name>_scalar` beside a vector) for scalar points. `<name>` is the function type (`frf`, `coherence`, `auto_spectrum`, `psd`, …), qualified with the reference (`frf_ref1z`) only when two kinds of the same type coexist.
- Values are multiplied by the sign of the response and of the reference direction, so a function measured towards `-Z` lands in the `+Z` component with its sign flipped.
- Nodes without a function are NaN; a response node missing from the geometry drops its function with a warning; a file with functions but no nodes (15, 781 or 2411) is an error.
- Functions of the same kind on a different abscissa grid become their own steps, with a warning.

Dataset **58b** is the binary variant: 11 ASCII header records, then the raw values. The byte order must be as declared (1 little endian, 2 big endian) and the floating-point format IEEE 754. The size of the binary block is taken from record 7 (count, precision, spacing, complex or real), not from the header line: pyuff, among others, declares half the size for complex data.

### Code_Aster comments

Code_Aster's UNV writer ends integer records with a `%` comment: the dataset line (`55   %VALEURS AUX NOEUDS`) and every node label of a 55/57 record (`1     % NOEUD N1`). An integer record stops at the first token starting with `%` (since v16.4.0; both engines refused these files before).

## Compatibility with files written by meshio++ ≤ 15.5

Older meshio++ versions wrote three things against the SDRL records, and read them back consistently:

- Group entity types were swapped (8 for nodes, 7 for elements). Such a file now reads its node groups as element groups and the other way round.
- `hexahedron20` and `wedge15` put the vertical mid-edges after the top ring. Such a file now reads with those mid-nodes misplaced.
- `code_aster=True` wrote a padded 55 header (a line of eight zeros where record 7 belongs) and element data as dataset 57. That layout is recognised by its invalid record 7 and still read, with a warning.

Symmetric tensors were also written unpermuted, so a 6-component array from an older file reads with its `xy`/`yy`/`zz`/`xz` components moved.

## Notes

- The record layouts follow the UC-SDRL dataset pages (the 55, 58, 2411, 2412, 2414, 2420, 2429 and 2467 texts are quoted in [pyuff](https://pypi.org/project/pyuff/)'s `get_structure_*` functions). Node orders and the group layouts were checked against gmsh's writer and Salome's reader.
- The fixtures under `tests/python/meshes/unv/` come from three independent writers, regenerated by `tools/gen_unv_fixtures.py`: gmsh 4.15.2 (quadratic tetrahedra, wedges and hexahedra, each with its `.msh` twin), pyuff 2.5.6 (units, 55 mode shapes, 58 and 58b frequency-response functions) and hand-written record layouts (2414 modes and transients, 2420, 15/780, the older meshio++ layout).
- Not read: dataset 57 and 2414 data at nodes on elements (element-node results), 2414 data at points, rotation of results between coordinate systems, the full dataset-58 DSP metadata (axis labels, units exponents).

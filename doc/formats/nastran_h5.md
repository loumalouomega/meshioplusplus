# MSC Nastran HDF5 results (`.h5`)

The result database [MSC Nastran](https://hexagon.com/products/product-groups/computer-aided-engineering-software/msc-nastran) writes when the bulk data carries `MDLPRM,HDF5,...` (MSC Nastran 2016 and later). Patran and post-processing scripts read it instead of the older OP2. It holds the model and the results in one self-describing HDF5 file of compound tables, so reading it needs no Nastran installation: `.h5` → `.vtu` / VTKHDF / any meshio++ format, with every subcase, mode, time or frequency as a step.

| | |
|---|---|
| **Format name** | `nastran_h5` |
| **Extensions** | `.h5` (GiD's longer `.post.h5` stays `gid`) |
| **Read / Write** | ✓ / — ([read-only by design](../conformance.md#nastran-h5): MSC Nastran writes it, and no tool reads one written elsewhere) |
| **Extra dependencies** | HDF5 (the C++ core); `h5py` for the Python reader |

## Reading

```python
import meshioplusplus

mesh = meshioplusplus.read("job.h5")                   # the first result domain
mode3 = meshioplusplus.read("job.h5", time_step=2)     # the third one

mesh.point_data["DISPLACEMENT"]          # (n, 3) translations
mesh.point_data["DISPLACEMENT_ROT"]      # (n, 3) rotations
mode3.point_data["EIGENVECTOR"]          # a mode shape
mesh.cell_data["STRESS:X"]               # per block, NaN where the card has no such output
mesh.cell_data["STRESS:X@corner"]        # (cells, corners): the corner values of solids and shells
mesh.cell_data["STRESS:X1@ply"]          # (cells, plies): composite ply stresses
mesh.field_data["nastran:subcase"], mesh.field_data["meshio:time"]

meshioplusplus.nastran_h5.time_values("job.h5")   # every step's TIME_FREQ_EIGR
```

```bash
meshioplusplus convert job.h5 'mode_{step}.vtu'    # one .vtu per step
meshioplusplus convert job.h5 job.pvd              # or a ParaView collection
```

Both engines read the whole format: the C++ core, and a Python reference reader on `h5py` that the core falls back to. `.h5` is read-only: `meshioplusplus.write(..., file_format="nastran_h5")` is an error.

## Which files

Only MSC Nastran's schema is read. The reader refuses, with a `ReadError` naming the reason and never a guess:

- a `.h5` without `/NASTRAN/INPUT/NODE/GRID`, i.e. any other HDF5 file under that extension;
- a file whose `/NASTRAN` `VERSION` attribute does not start with `msc`. Other vendors write incompatible HDF5 layouts under similar names, and a file that reads without error but means something else is worse than no file. A file without the attribute is read.

The schema grows with every MSC release. The reader relies on nothing specific to one release: it reads the tables and members it knows by name, and treats a table it does not know as described below. It was checked against MSC Nastran 2020 output (`SCHEMA` 20200).

## The model

- **Points** are the `GRID` rows, in file order, in the basic coordinate system; see [Coordinate systems](#coordinate-systems). `SPOINT`/`EPOINT` scalar points are not points.
- **Cells.** Each `/NASTRAN/INPUT/ELEMENT/<card>` table with a cell type gives up to two cell blocks, in card-name order:

| Card | Cell type | Quadratic |
|---|---|---|
| `CTRIA3`, `CTRIAR` | `triangle` | |
| `CTRIA6` | `triangle` | `triangle6` |
| `CQUAD4`, `CQUADR`, `CSHEAR` | `quad` | |
| `CQUAD8` | `quad` | `quad8` |
| `CQUAD` | `quad` | `quad9` |
| `CTETRA` | `tetra` | `tetra10` |
| `CPYRAM` | `pyramid` | `pyramid13` |
| `CPENTA` | `wedge` | `wedge15` |
| `CHEXA` | `hexahedron` | `hexahedron20` |
| `CROD`, `CTUBE`, `CONROD`, `CBAR`, `CBEAM`, `CBUSH`, `CVISC`, `PLOTEL` | `line` | |
| `CONM2` | `vertex` | |

- **Order.** MSC stores the node list of a solid at its full quadratic width, zero-padded (a linear `CHEXA` has 20 slots, 12 of them 0). An element with every mid-side node is quadratic, one with none is linear, and one with only some (a `CTRIA6` listing `6 4 60 0 65 64`) is read as the linear type with a warning. A card's linear elements come first, then its quadratic ones.
- **Node order.** `CHEXA`/`CPENTA` list their mid-side nodes bottom ring, vertical edges, top ring; meshio++ lists bottom, top, vertical. They are permuted on read (`connectivity[k] = G[perm[k]]`, hexahedron20 `0..11, 16..19, 12..15`, wedge15 `0..8, 12, 13, 14, 9, 10, 11`), the same tables as the [bulk Nastran](./nastran.md) reader. `CTETRA`, `CPYRAM`, `CTRIA6`, `CQUAD8` need none. The tests build each quadratic card from its edge midpoints and check that every mid-node lies on its edge.
- **Springs and dampers** (since v16.12.0): `CELAS1`, `CELAS2`, `CDAMP1` and `CDAMP2` are lines between their two GRIDs, or vertices when grounded (one of them 0).
- **Skipped cards.** Scalar-point elements (`CELAS3/4`, `CDAMP3/4`, `CMASS*`) and any card with no cell type (`CHBDYE`, `CBUSH1D`...) are skipped with one warning listing them. An element that connects a SPOINT is skipped too; one that names a GRID the file does not define is an error.
- **Cell data.** `nastran:eid` (the element id) and `nastran:pid` (the property id, -1 for a card without one, such as `CONM2` or `CONROD`), integer, per block.
- **Regions.** One [cell region](../regions.md) per property id, named after the property card that defines it (`PSHELL_4`, `PSOLID_2`, or `PID_<n>` when no card does; a flat table such as `PSHELL` wins over a grouped one such as `PCOMP/IDENTITY` for an id both define), with the property id as its tag and the highest cell dimension among its cells as its dimension.

## Coordinate systems

Since v16.10.0 a GRID whose `CP` is not 0 is moved from its local system to the basic one, and every nodal vector result of a GRID whose `CD` is not 0 (the system its results are output in) is rotated to basic at that point. The systems are the `CORD1R`/`CORD1C`/`CORD1S` and `CORD2R`/`CORD2C`/`CORD2S` tables of `/NASTRAN/INPUT/COORDINATE_SYSTEM`. They are resolved in any order: a `CORD1*` through its three GRIDs (which may lie in another local system) and a `CORD2*` through its reference system. Cylindrical `(r, θ, z)` and spherical `(r, θ, φ)` coordinates, in degrees, are converted to the system's Cartesian axes first. A vector at a point of a cylindrical or spherical system is rotated by the local axes at that point, so the rotation differs from point to point.

What is rotated: every nodal table with `X, Y, Z` members (and its `RX, RY, RZ` rotations) and all of `GRID_FORCE`. Scalar tables, element results (in the element's own system) and complex tables stay as written. `point_data["nastran:cp"]` and `["nastran:cd"]` still record the system ids when any is non-zero, and an info line says how many GRIDs moved. A system that cannot be resolved (a missing `CORD*` card, or a cycle) leaves its GRIDs and their results as written, with a warning naming the ids. The bulk-data [`nastran`](./nastran.md) reader does not apply `CP`.

## Results and steps

Every result table has an INDEX twin, `/INDEX/NASTRAN/RESULT/<same path>`, whose rows give the position and length of the table's rows for each **domain**. A domain is one row of `/NASTRAN/RESULT/DOMAINS`: a subcase, and within it a mode, a time step or a frequency. The domains any nodal or elemental INDEX references are the steps of the [sequence engine](../sequences.md), in `DOMAINS` order. `time_step` picks one (`0` is the first, negative counts from the end, out of range is an error), `read_metadata(...).time_values` lists them, and `read_sequence` walks them. Only the chosen domain's rows are read, through HDF5 hyperslabs. A `DOMAINS` row no table uses (the eigenvalue summary's own domains, for instance) is not a step. A file with no results is one step with no field data.

| Field data | Meaning |
|---|---|
| `meshio:time` | the domain's `TIME_FREQ_EIGR`, as written: the time of a transient run, the frequency of a frequency response, the eigenvalue of a mode (for a real mode that is λ = ω², not a frequency; for a buckling mode the load factor) |
| `nastran:domain` | the `DOMAINS` id |
| `nastran:subcase`, `nastran:step`, `nastran:mode` | the domain's subcase, step and mode number |
| `nastran:analysis` | the domain's `ANALYSIS` code, as written |
| `nastran:eigi` | the imaginary part of a complex eigenvalue |

`meshio:time` is not converted to a frequency: the `ANALYSIS` codes are not documented publicly, and the one that looks like "modes" also carries buckling load factors in a SOL 105 file, so no conversion could be applied safely. Rigid-body modes all have eigenvalue 0, so a `.pvd` written from them repeats `timestep="0.0"`; write one file per step (`'mode_{step}.vtu'`) or select them by `nastran:mode`.

A result table can cover only some domains: in a SOL 105 file `DISPLACEMENT` is in the static domain and `EIGENVECTOR` in the buckling ones. Each step has the arrays of the tables that have rows for it.

### Nodal results

Each `/NASTRAN/RESULT/NODAL/<T>` table with at most one row per node and domain becomes point data. Members are grouped by name:

| Members | Point data |
|---|---|
| `X, Y, Z` | `<T>` `(n, 3)` |
| `RX, RY, RZ` | `<T>_ROT` `(n, 3)` |
| `XR, YR, ZR` / `XI, YI, ZI` (a `<B>_CPLX` table) | `<B>_real` / `<B>_imag` `(n, 3)` |
| `RXR, RYR, RZR` / `RXI, RYI, RZI` | `<B>_ROT_real` / `<B>_ROT_imag` |
| a lone `VALUE` (`TEMPERATURE`) | `<T>` `(n,)` |
| any other float member `M` (`KINETIC_ENERGY`'s `KET1`...) | `<T>:<M>` `(n,)` |

A point the table has no row for is NaN. MSC writes rows for SPOINTs into the same tables (a `EIGENVECTOR` of 40 GRIDs and 3 SPOINTs has 43 rows per mode); those rows are dropped. `GRID_FORCE`, which has one row per (node, element), is read as described under [Several values per element](#several-values-per-element).

### Element results

Each `/NASTRAN/RESULT/ELEMENTAL/<G>/<T>` table (`<G>` is `STRESS`, `STRAIN`, `ELEMENT_FORCE`, `ENERGY`...) with at most one row per element and domain becomes cell data, one array per float member `M`, named `<G>:<M>` and shared by every table of the group: `STRESS:X` holds the `X` of the `HEXA`, `PENTA` and `TETRA` tables alike, `STRESS:X1` the fibre-1 stress of the shell tables. A member that is an array (a solid's `X(9)`: the centre and then each corner; a corner-output shell's `(5)`; a beam's 11 stations) contributes its **first** entry to `<G>:<M>`, which is the centre for solids and `CEN` shells and end A for a beam; its other entries are read as described next. A cell whose card has no such member is NaN. Rows are joined to cells by `EID` (`ID` for the `ENERGY` tables).

### Several values per element

Since v16.10.0 the values beyond the centre are read too. They are added beside `<G>:<M>`, which keeps the centre value, as `(cells, columns)` cell data, NaN-padded, with `field_data["nastran:layout:<name>"] = [columns, 1]`:

| Array | Source | Column |
|---|---|---|
| `<G>:<M>@corner` | the corner entries of a solid's or corner-output shell's array member | the GRID's position in the cell's connectivity |
| `<G>:<M>@ply` | the per-ply rows of a `_COMP` table | `PLY - 1` |
| `<G>:<M>@station` | a `BEAM` table's 11 stations, and the `SD` rows of `BARS` | the station; stations not output stay NaN |
| `GRID_FORCE:<M>` (`F1`...`M3`) | the element rows of `GRID_FORCE` | the GRID's position in the cell's connectivity |

The `GRID_FORCE` rows that belong to no element (`*TOTALS*`, `APP-LOAD`, `F-OF-SPC`...) become point data `GRID_FORCE:<label>:<M>`, the label being the row's `ELNAME` without blanks and asterisks (`GRID_FORCE:TOTALS:F1`). A cell whose card has no such values is NaN in every column. A table with repeated rows of another kind (`ENERGY`'s summary row with ID 100000000) is still skipped with a warning.

`CONM2` elements are not matched by element results: they have none, and MSC accepts a `CONM2` sharing its id with a structural element. Any other id shared by two cards gets a warning, and the results go to the first.

## What is not read

- **Everything outside `/NASTRAN/INPUT/NODE`, `/ELEMENT`, `/PROPERTY` and `/NASTRAN/RESULT/NODAL`, `/ELEMENTAL`**: loads, constraints, materials, the eigenvalue summary, optimization and aerodynamic results.
- **Other vendors' HDF5 schemas**, refused as above until a file is in hand. OP2 files are read by [`nastran_op2`](./nastran_op2.md), with the same mesh and member names; OP4 and punch files are not read.

## Verification

The test suite reads six files MSC Nastran 2020 wrote (SOL 101, 103, 105, 107, 108 and 159, from the [pyNastran](https://github.com/SteveDoyle2/pyNastran) repository, BSD 3-Clause, credited in `tests/python/meshes/nastran_h5/README.md`). Every translation and rotation vector of every domain — 38 tables, the SOL 103 eigenvectors included — equals pyNastran 1.4.1's own reading of the same file bit for bit, and a `.vtu` sequence converted from the SOL 103 file, read back by ParaView, carries the same eigenvectors. The two engines agree exactly on every array of every step of every file, and every solid has a positive volume. The corner, ply, station and grid point force values match pyNastran's on the same files. The coordinate systems are checked on a probe derived from a real file, with GRIDs in a tilted `CORD2C`, a `CORD2S` defined in it and a `CORD1R` through GRIDs of the `CORD2C` (`tools/gen_nastran_cord_reference.py`): points match pyNastran's to 1e-16 and displacements to 4e-16, compared with upstream pyNastran's rotation matrices, because 1.4.1 composes cylindrical and spherical ones the wrong way round.

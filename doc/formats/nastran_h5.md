# MSC Nastran HDF5 results (`.h5`)

The result database [MSC Nastran](https://hexagon.com/products/product-groups/computer-aided-engineering-software/msc-nastran) writes when the bulk data carries `MDLPRM,HDF5,...` (MSC Nastran 2016 and later). Patran and post-processing scripts read it instead of the older OP2. It holds the model and the results in one self-describing HDF5 file of compound tables, so reading it needs no Nastran installation: `.h5` → `.vtu` / VTKHDF / any meshio++ format, with every subcase, mode, time or frequency as a step.

| | |
|---|---|
| **Format name** | `nastran_h5` |
| **Extensions** | `.h5` (GiD's longer `.post.h5` stays `gid`) |
| **Read / Write** | ✓ / — (read-only) |
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

- **Points** are the `GRID` rows, in file order. Their coordinates are taken **as written**. A GRID whose `CP` is not 0 is given in a local coordinate system, which is not applied: a warning names how many, and `point_data["nastran:cp"]` keeps the system id (only when some CP is non-zero). A non-zero `CD` (the system the results are output in) gets a warning and `point_data["nastran:cd"]` the same way. `SPOINT`/`EPOINT` scalar points are not points.
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
- **Skipped cards.** Scalar-point elements (`CELAS*`, `CDAMP*`, `CMASS*`) and any card with no cell type (`CHBDYE`, `CBUSH1D`...) are skipped with one warning listing them. An element that connects a SPOINT is skipped too; one that names a GRID the file does not define is an error.
- **Cell data.** `nastran:eid` (the element id) and `nastran:pid` (the property id, -1 for a card without one, such as `CONM2` or `CONROD`), integer, per block.
- **Regions.** One [cell region](../regions.md) per property id, named after the property card that defines it (`PSHELL_4`, `PSOLID_2`, or `PID_<n>` when no card does; a flat table such as `PSHELL` wins over a grouped one such as `PCOMP/IDENTITY` for an id both define), with the property id as its tag and the highest cell dimension among its cells as its dimension.

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

A point the table has no row for is NaN. MSC writes rows for SPOINTs into the same tables (a `EIGENVECTOR` of 40 GRIDs and 3 SPOINTs has 43 rows per mode); those rows are dropped. `GRID_FORCE`, which has one row per (node, element), is skipped with a warning.

### Element results

Each `/NASTRAN/RESULT/ELEMENTAL/<G>/<T>` table (`<G>` is `STRESS`, `STRAIN`, `ELEMENT_FORCE`, `ENERGY`...) with at most one row per element and domain becomes cell data, one array per float member `M`, named `<G>:<M>` and shared by every table of the group: `STRESS:X` holds the `X` of the `HEXA`, `PENTA` and `TETRA` tables alike, `STRESS:X1` the fibre-1 stress of the shell tables. A member that is an array (a solid's `X(9)`: the centre and then each corner; a corner-output shell's `(5)`; a beam's 11 stations) contributes its **first** entry, which is the centre for solids and `CEN` shells and end A for a beam. A cell whose card has no such member is NaN. Rows are joined to cells by `EID` (`ID` for the `ENERGY` tables). Per-ply `_COMP` tables, with one row per ply, are skipped with a warning.

`CONM2` elements are not matched by element results: they have none, and MSC accepts a `CONM2` sharing its id with a structural element. Any other id shared by two cards gets a warning, and the results go to the first.

## What is not read

- **Coordinate systems.** `CP`/`CD` are recorded, not applied (above).
- **Everything outside `/NASTRAN/INPUT/NODE`, `/ELEMENT`, `/PROPERTY` and `/NASTRAN/RESULT/NODAL`, `/ELEMENTAL`**: loads, constraints, materials, the eigenvalue summary, optimization and aerodynamic results.
- **Per-ply, per-station and per-grid values** beyond the first entry of an array member.
- **Other vendors' HDF5 schemas**, refused as above until a file is in hand. OP2 files are read by [`nastran_op2`](./nastran_op2.md), with the same mesh and member names; OP4 and punch files are not read.

## Verification

The test suite reads six files MSC Nastran 2020 wrote (SOL 101, 103, 105, 107, 108 and 159, from the [pyNastran](https://github.com/SteveDoyle2/pyNastran) repository, BSD 3-Clause, credited in `tests/python/meshes/nastran_h5/README.md`). Every translation and rotation vector of every domain — 38 tables, the SOL 103 eigenvectors included — equals pyNastran 1.4.1's own reading of the same file bit for bit, and a `.vtu` sequence converted from the SOL 103 file, read back by ParaView, carries the same eigenvectors. The two engines agree exactly on every array of every step of every file, and every solid has a positive volume.

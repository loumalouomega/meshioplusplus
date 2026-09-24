# Nastran OP2 results (`.op2`)

The binary result file [MSC Nastran](https://hexagon.com/products/product-groups/computer-aided-engineering-software/msc-nastran) and [Simcenter Nastran](https://plm.sw.siemens.com/en-US/simcenter/mechanical-simulation/nastran/) (NX Nastran) write with `PARAM,POST,-1` (or `-2`), and that Femap, Patran and most post-processors read. MSC runs that write HDF5 are read by [`nastran_h5`](./nastran_h5.md); OP2 is the route for Simcenter/NX Nastran and for MSC decks without `MDLPRM,HDF5`. meshio++ reads it as a sequence, `.op2` → `.vtu` / VTKHDF / any meshio++ format, with every subcase, mode, time or frequency as a step.

| | |
|---|---|
| **Format name** | `nastran_op2` |
| **Extensions** | `.op2` |
| **Read / Write** | ✓ / — (read-only) |
| **Extra dependencies** | none |

## Reading

```python
import meshioplusplus

mesh = meshioplusplus.read("job.op2")                 # the first step
mode3 = meshioplusplus.read("job.op2", time_step=2)   # the third one

mesh.point_data["DISPLACEMENT"]        # (n, 3) translations
mesh.point_data["DISPLACEMENT_ROT"]    # (n, 3) rotations
mode3.point_data["EIGENVECTOR"]        # a mode shape
mesh.cell_data["STRESS:VON_MISES"]     # solids' centre von Mises, NaN elsewhere
mesh.cell_data["STRESS:VON_MISES1"]    # shells' fibre-1 von Mises
mesh.field_data["nastran:subcase"], mesh.field_data["meshio:time"]

meshioplusplus.nastran_op2.time_values("job.op2")
```

```bash
meshioplusplus convert job.op2 'mode_{step}.vtu'
```

Both engines read the whole format: the C++ core, and a Python reference reader the core falls back to. `.op2` is read-only.

## Which files

An OP2 is a Fortran unformatted file of named tables. meshio++ reads MSC and NX/Simcenter files, 32-bit (4-byte words) and 64-bit (8-byte words, still with 4-byte record markers), in either byte order. The table names of 64-bit NX files, written four characters per word, are read as the 32-bit ones. A formatted (text) OP2 is refused. A file is recognised by its `.op2` extension, or by content: the `PARAM,POST,-1` header (a date and the tape code record) is sniffed in any word size and byte order.

## The model

The mesh is built as the [`nastran_h5`](./nastran_h5.md) reader builds it, from the same card table (`detail/nastran_model.hpp`):

- **Points** are the `GEOM1` GRIDs, in file order. Coordinates are taken as written: a GRID with `CP != 0` (or `CD != 0`) gets a warning and `point_data["nastran:cp"]` (`["nastran:cd"]`). NX's 32-bit GRID records with double-precision coordinates are read as such.
- **Cells.** The `GEOM2` element records with a cell type: `CQUAD4`, `CQUADR`, `CSHEAR`, `CQUAD8`, `CQUAD`, `CTRIA3`, `CTRIAR`, `CTRIA6`, `CTETRA`, `CPYRAM`, `CPENTA`, `CHEXA` (linear or quadratic by their mid-side nodes), `CROD`, `CTUBE`, `CONROD`, `CBAR`, `CBEAM`, `CBUSH`, `CVISC`, `PLOTEL` as lines and `CONM2` as vertices, in card-name order; NX 2019 and later's `CQUAD4`/`CTRIA3` records (keys 15401 and 15301) are read as those cards. A card whose record length fits more than one layout (MSC and NX `CQUAD4`s differ by a word) takes the first whose entries validate: ids increasing, corner nodes that are GRIDs. Records with no cell type (`CELAS*`, `CDAMP*`, `CMASS*`, `CGAP`...) are named in one warning.
- **Cell data.** `nastran:eid` and `nastran:pid` (-1 for `CONM2`, `CONROD`, `PLOTEL`).
- **Regions.** One cell region per property id, named `<PTYPE>_<pid>` after its `EPT` record (`PSHELL`, `PSOLID`, `PCOMP`, `PCOMPG`, `PBAR`, `PBARL`, `PBEAM`, `PBEAML`, `PROD`, `PSHEAR`, `PTUBE`, `PBUSH`...), `PID_<pid>` when none defines it. NX writes each `PCOMP` a second time as a `PSHELL` with material ids from 100000000: that twin is ignored, and a `PCOMPG` wins over a `PCOMP` of the same id, as pyNastran reads them.
- **Without geometry.** A file with no `GEOM1` GRIDs (written without the model, or with `PARAM,POST,-2` defaults) takes its mesh from the input deck beside it: `<stem>.bdf`, `.dat`, `.nas` or `.blk`, read by the [bulk Nastran](./nastran.md) reader, with `nastran:eid` added. Without one, the `ReadError` lists the paths it looked for. The bulk reader does not follow `INCLUDE` cards.

## Results and steps

A results table alternates a 146-word header (approach, table and element-type codes, subcase, mode or time, width) with the data it describes. Every distinct (subcase, analysis code, word 5 of the header: load set, mode number, time or frequency) of the tables below is a step, in the order the file first names it: `time_step` picks one, `read_metadata(...).time_values` lists them and `read_sequence` walks them.

| Field data | Meaning |
|---|---|
| `meshio:time` | the time (transient) or frequency, the eigenvalue of a mode (real modes: λ = ω²) or the post-buckling eigenvalue, 0 for a static subcase; as written, not converted |
| `nastran:subcase`, `nastran:analysis` | the subcase and the approach's analysis code (1 statics, 2 modes, 5 frequency, 6 transient, 8 post-buckling...) |
| `nastran:mode` | the mode number (analysis 2, 8 and 9), else 0 |

### Nodal results

Real SORT1 tables with 8 words per node (`OUG*`, `BOUG*`, `OQG*`, `OQMG*`, `OPG*`) become point data by their table code, with the `nastran_h5` names: `DISPLACEMENT`, `EIGENVECTOR`, `VELOCITY`, `ACCELERATION`, `SPC_FORCE`, `MPC_FORCE` (the `OQMG` tables, which share the SPC forces' code) and `APPLIED_LOAD`, each `(n, 3)` with a `_ROT` twin; a thermal displacement table is `TEMPERATURE` `(n,)`. A point with no row is NaN, and scalar-point rows are dropped. When two tables give the same point a value in the same step, the first one read keeps it.

### Element results

Real SORT1 stress (`OES*`) and strain (`OSTR*`) tables become cell data `STRESS:<M>` or `STRAIN:<M>`, the **centre** value per element, NaN on cells without it:

| Element type (code) | Members |
|---|---|
| `CROD` (1), `CONROD` (10) | `A`, `MSA`, `T`, `MST` |
| `CTUBE` (3) | `AS`, `MSA`, `TS`, `MST` |
| `CSHEAR` (4) | `TMAX`, `TAVG`, `MS` |
| `CBAR` (34) | `X1A`...`X4A`, `AX`, `MAXA`, `MINA`, `MST`, `X1B`...`X4B`, `MAXB`, `MINB`, `MSC` |
| `CQUAD4` (33), `CTRIA3` (74), and the centre of `CQUAD8` (64), `CTRIAR` (70), `CTRIA6` (75), `CQUADR` (82), corner `CQUAD4` (144) | per fibre `k` = 1, 2: `FD<k>`, `X<k>`, `Y<k>`, `TXY<k>`, `ANGLE<k>`, `MAJOR<k>`, `MINOR<k>`, `VON_MISES<k>` (or `MAX_SHEAR<k>`) |
| `CTETRA` (39), `CHEXA` (67), `CPENTA` (68), `CPYRAM` (255) | `X`, `Y`, `Z`, `TXY`, `TYZ`, `TZX`, `PRINCIPAL_A`, `PRINCIPAL_B`, `PRINCIPAL_C` (in the file's A, B, C order, not sorted), `PRESSURE`, `VON_MISES` (or `OCT_SHEAR`) |

The members common to [`nastran_h5`](./nastran_h5.md) (`X1`, `TXY1`, `X`, `TZX`, `A`, `TMAX`...) have its names, so an OP2 and an HDF5 file of the same run give the same arrays; the derived values the HDF5 schema does not store have their own names. Whether the last plate and solid column is von Mises or the maximum (octahedral) shear, and whether a table holds stress or strain, come from the header's stress code. Elements share results by id; `CONM2` elements, which have none, are left out.

**Not read**, and named in one warning: complex, random and SORT2 tables, the corner and per-ply values (the `nastran_h5` reader keeps only the centre too: see the [roadmap](../roadmap.md)), every other element type (`CBEAM`, springs, composites, NX's newer solids...), and every other table (`OEF` forces, `OGPFB` grid point forces, `ONRGY` energies, contact and bolt tables...).

## Verification

The test suite reads eleven files from the [pyNastran](https://github.com/SteveDoyle2/pyNastran) test models (BSD 3-Clause, credited in `tests/python/meshes/nastran_op2/README.md`): MSC and NX output, 32- and 64-bit, SOL 101 statics (among them `static_solid_shell_bar`, a solid, shell and bar model), modes, buckling, transient heat transfer, an NX 2412 file, and one with its geometry tables removed next to its deck. On every step of every file the nodal vectors and the centre stress and strain equal pyNastran 1.4.1's reading to float precision, and the two engines agree exactly; the same holds for 40 further pyNastran models checked outside the suite. No file from an MSC or NX version newer than those samples, and no big-endian file Nastran wrote, has been read.

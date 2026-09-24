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
mesh.cell_data["STRESS:X@corner"]      # (cells, nodes): solid corner values
mesh.cell_data["GRID_FORCE:F1"]        # (cells, nodes): grid point forces
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

- **Points** are the `GEOM1` GRIDs, in file order, in the basic coordinate system (since v16.10.0): a GRID with `CP != 0` is moved from its local system, which the `GEOM1` `CORD1R/C/S` and `CORD2R/C/S` records define, and `point_data["nastran:cp"]` (`["nastran:cd"]`) records the system ids. The rules are those of [`nastran_h5`](./nastran_h5.md#coordinate-systems). NX's 32-bit GRID records with double-precision coordinates are read as such. A mesh taken from a sibling deck keeps the deck reader's raw coordinates.
- **Cells.** The `GEOM2` element records with a cell type: `CQUAD4`, `CQUADR`, `CSHEAR`, `CQUAD8`, `CQUAD`, `CTRIA3`, `CTRIAR`, `CTRIA6`, `CTETRA`, `CPYRAM`, `CPENTA`, `CHEXA` (linear or quadratic by their mid-side nodes), `CROD`, `CTUBE`, `CONROD`, `CBAR`, `CBEAM`, `CBUSH`, `CVISC`, `PLOTEL` as lines, `CONM2` as vertices, and (since v16.12.0) the springs and dampers `CELAS1`, `CELAS2`, `CDAMP1`, `CDAMP2` as lines between their two GRIDs or vertices when grounded (one of them 0), in card-name order; NX 2019 and later's `CQUAD4`/`CTRIA3` records (keys 15401 and 15301) are read as those cards. A card whose record length fits more than one layout (MSC and NX `CQUAD4`s differ by a word) takes the first whose entries validate: ids increasing, corner nodes that are GRIDs. Records with no cell type (`CELAS3/4` and `CDAMP3/4`, which join scalar points, `CMASS*`, `CGAP`...) are named in one warning, and a spring with a scalar point at either end is skipped.
- **Cell data.** `nastran:eid` and `nastran:pid` (-1 for `CONM2`, `CONROD`, `PLOTEL`).
- **Regions.** One cell region per property id, named `<PTYPE>_<pid>` after its `EPT` record (`PSHELL`, `PSOLID`, `PCOMP`, `PCOMPG`, `PBAR`, `PBARL`, `PBEAM`, `PBEAML`, `PROD`, `PSHEAR`, `PTUBE`, `PBUSH`...), `PID_<pid>` when none defines it. NX writes each `PCOMP` a second time as a `PSHELL` with material ids from 100000000: that twin is ignored, and a `PCOMPG` wins over a `PCOMP` of the same id, as pyNastran reads them.
- **Without geometry.** A file with no `GEOM1` GRIDs (written without the model, or with `PARAM,POST,-2` defaults) takes its mesh from the input deck beside it: `<stem>.bdf`, `.dat`, `.nas` or `.blk`, read by the [bulk Nastran](./nastran.md) reader, with `nastran:eid` added. The bulk reader does not follow `INCLUDE` cards. Without a deck (since v16.12.0) the points come from the basic grid point table: `BGPDTS` rows (output system and basic coordinates, in internal order) named by the `EQEXINS` table, or NX's newer `BGPDT`, whose rows carry the GRID id; scalar points are left out, and elements still come from `GEOM2` when the file has it. Their output systems are unknown without `GEOM1`, so results in a non-basic `CD` stay in it, with a warning. With neither, the `ReadError` lists the paths it looked for. A fluid (acoustic) GRID's `CD = -1` reads as the basic system.
- **`PARAM,POST,-2`** files carry no header: the first table's name opens them. They read as a `POST,-1` file does, and are recognised without the `.op2` extension by that name (a one-word 2, an eight-character upper-case name, a one-word -1).

## Results and steps

A results table alternates a 146-word header (approach, table and element-type codes, subcase, mode or time, width) with the data it describes. Every distinct (subcase, analysis code, word 5 of the header: load set, mode number, time or frequency) of the tables below is a step, in the order the file first names it; grid point force and energy tables join a step the other tables made (their word 5 is 0 in a static subcase), and an energy table makes its own only when no table names its subcase and analysis. `time_step` picks one, `read_metadata(...).time_values` lists them and `read_sequence` walks them.

| Field data | Meaning |
|---|---|
| `meshio:time` | the time (transient) or frequency, the eigenvalue of a mode (real modes: λ = ω²) or the post-buckling eigenvalue, 0 for a static subcase; as written, not converted |
| `nastran:subcase`, `nastran:analysis` | the subcase and the approach's analysis code (1 statics, 2 modes, 5 frequency, 6 transient, 8 post-buckling...) |
| `nastran:mode` | the mode number (analysis 2, 8 and 9), else 0 |
| `nastran:eigi` | a complex mode's imaginary eigenvalue (analysis 9), with `meshio:time` its real part, both from the `CLAMA` table (since v16.12.0) |

### Nodal results

Real SORT1 tables with 8 words per node (`OUG*`, `BOUG*`, `OQG*`, `OQMG*`, `OPG*`) become point data by their table code, with the `nastran_h5` names: `DISPLACEMENT`, `EIGENVECTOR`, `VELOCITY`, `ACCELERATION`, `SPC_FORCE`, `MPC_FORCE` (the `OQMG` tables, which share the SPC forces' code) and `APPLIED_LOAD`, each `(n, 3)` with a `_ROT` twin; a thermal displacement table is `TEMPERATURE` `(n,)`. A point with no row is NaN, and scalar-point rows are dropped. When two tables give the same point a value in the same step, the first one read keeps it. Vectors of a GRID with `CD != 0` are rotated from that system to basic (since v16.10.0), except in the `BOUG*` tables, which Nastran already writes in basic.

Since v16.12.0:

- **Complex tables** (frequency response, complex modes; 14 words per node) give `<name>_real` and `<name>_imag` (`DISPLACEMENT_real`, `DISPLACEMENT_ROT_imag`...), the `nastran_h5` names. Magnitude and phase (format 3) are converted, the phase in degrees. Both parts are rotated to basic.
- **Random tables** (`OUGPSD`, `OUGRMS`, `OQGNO`...: the table code's hundreds name the quantity) give `<name>_PSD`, `_ATO`, `_RMS`, `_NO` or `_CRM`. They are spectral quantities, not vectors, so they stay in the GRID's output system.
- **SORT2 tables** (one GRID or element over every step, each row led by its time, frequency or mode) are pivoted: every row joins the step its first word names, so a SORT2 table and its SORT1 twin give the same arrays.

### Element results

Stress (`OES*`) and strain (`OSTR*`) tables become cell data `STRESS:<M>` or `STRAIN:<M>`, the **centre** value per element (a beam's first station), NaN on cells without it:

| Element type (code) | Members |
|---|---|
| `CROD` (1), `CONROD` (10) | `A`, `MSA`, `T`, `MST` |
| `CTUBE` (3) | `AS`, `MSA`, `TS`, `MST` |
| `CSHEAR` (4) | `TMAX`, `TAVG`, `MS` |
| `CBAR` (34) | `X1A`...`X4A`, `AX`, `MAXA`, `MINA`, `MST`, `X1B`...`X4B`, `MAXB`, `MINB`, `MSC` |
| `CQUAD4` (33), `CTRIA3` (74), and the centre of `CQUAD8` (64), `CTRIAR` (70), `CTRIA6` (75), `CQUADR` (82), corner `CQUAD4` (144) | per fibre `k` = 1, 2: `FD<k>`, `X<k>`, `Y<k>`, `TXY<k>`, `ANGLE<k>`, `MAJOR<k>`, `MINOR<k>`, `VON_MISES<k>` (or `MAX_SHEAR<k>`) |
| `CTETRA` (39), `CHEXA` (67), `CPENTA` (68), `CPYRAM` (255) | `X`, `Y`, `Z`, `TXY`, `TYZ`, `TZX`, `PRINCIPAL_A`, `PRINCIPAL_B`, `PRINCIPAL_C` (in the file's A, B, C order, not sorted), `PRESSURE`, `VON_MISES` (or `OCT_SHEAR`) |
| `CBEAM` (2) | `SD`, `XC`, `XD`, `XE`, `XF`, `MAX`, `MIN`, `MST`, `MSC` |
| `CBAR` stations (100) | `SD`, `XC`, `XD`, `XE`, `XF`, `AX`, `MAX`, `MIN`, `MS` |
| composite `CQUAD4` (95), `CQUAD8` (96), `CTRIA3` (97), `CTRIA6` (98), NX `CQUADR` (232), `CTRIAR` (233) | per ply: `X1`, `Y1`, `T1`, `L1`, `L2`, `ANGLE`, `MAJOR`, `MINOR`, `VON_MISES` (or `MAX_SHEAR`) |
| `CELAS1`...`CELAS4` (11-14) | `S` |
| `CBUSH` (102) | `TX`, `TY`, `TZ`, `RX`, `RY`, `RZ` |
| NX `CHEXA` (300), `CPENTA` (301), `CTETRA` (302), `CPYRAM` (303) | corners only (no centre value): `X`, `Y`, `Z`, `TXY`, `TYZ`, `TZX`, `VON_MISES` as `<M>@corner` |

Since v16.12.0 the other element tables are read too:

- **Complex** stresses and strains (rods, springs, `CSHEAR`, `CBAR`, `CBEAM`, centroidal and corner shells, solids, composites, `CBUSH`, and NX's von Mises variants `OESVM1`/`OSTRVM1`) give `<M>_real` and `<M>_imag` for each complex member; the real ones (fibre distances, von Mises of the `VM` tables) keep their names.
- **Random** stresses and strains (`OESPSD`, `OESRMS`...) have their own magnitude-only layouts: rods `A`, `T`; `CSHEAR` `TMAX`, `TAVG`; `CBAR` `X1A`...`AX`...`X4B`; `CBEAM` stations; shells `FD<k>`, `X<k>`, `Y<k>`, `TXY<k>` (and `VON_MISES<k>`); solids `X`...`TZX` (and `VON_MISES`); composites per ply; `CBUSH`. The names carry the random suffix (`STRESS:X1_RMS`). No file with random element tables was available: these layouts follow pyNastran's reader and are not checked against one.
- **Element forces** (`OEF*`) become `ELEMENT_FORCE:<M>`, real, complex (`_real`, `_imag`) or random (suffixed): rods, `CTUBE`, `CONROD` and `CVISC` `AF`, `TRQ`; springs and dampers `F`; `CBAR` `BM1A`, `BM2A`, `BM1B`, `BM2B`, `TS1`, `TS2`, `AF`, `TRQ`; `CBAR` stations (100) `SD`, `BM1`, `BM2`, `TS1`, `TS2`, `AF`, `TRQ` per station; `CSHEAR` `F41`...`F14`, `KF1`, `S12`...`KF4`, `S41`; shells `MX`, `MY`, `MXY`, `BMX`, `BMY`, `BMXY`, `TX`, `TY` (centre, and `@corner` for corner output); `CBUSH` `FX`...`MZ`; `CGAP` `FX`, `SFY`, `SFZ`, `U`, `V`, `W`, `SV`, `SW`; `CBEAM` `SD`, `BM1`, `BM2`, `TS1`, `TS2`, `AF`, `TTRQ`, `WTRQ` per station. The complex `CSHEAR` layout is the forces' real then imaginary parts, then the kick forces' and shear flows' (pyNastran pairs them otherwise).
- **Heat fluxes** (an `OEF` of a heat-transfer run) become `HEAT_FLUX:<M>`: line, shell and solid elements `XGRAD`, `YGRAD`, `ZGRAD`, `XFLUX`, `YFLUX`, `ZFLUX`; `CHBDYE/G/P` `FAPPLIED`, `FREECONV`, `FORCECONV`, `FRAD`, `FTOTAL`; `CONV` `FREECONV`, `FREECONVK`. `CHBDY*` elements have no cell, so their fluxes are not kept.
- **Energies**: strain energy (`ONRGY*`) becomes `ENERGY:ENERGY`, `ENERGY:PCT`, `ENERGY:DEN`, and kinetic energy (`OEKE*`) `KINETIC_ENERGY:<M>` likewise; complex energy has `ENERGY_real`/`_imag`. The total row of each set is no element and is dropped.

The members common to [`nastran_h5`](./nastran_h5.md) (`X1`, `TXY1`, `X`, `TZX`, `A`, `TMAX`...) have its names, so an OP2 and an HDF5 file of the same run give the same arrays; the derived values the HDF5 schema does not store have their own names. Whether the last plate and solid column is von Mises or the maximum (octahedral) shear, and whether a table holds stress or strain, come from the header's stress code. Elements share results by id; `CONM2` elements, which have none, are left out.

Since v16.10.0 the values beyond the centre are read too, with the arrays and `nastran:layout:<name>` field data of [`nastran_h5`](./nastran_h5.md#several-values-per-element): `<name>@corner` for the corners of solids (39, 67, 68, 255) and corner-output shells (64, 70, 75, 82, 144), `<name>@ply` for the composite plies (column `PLY - 1`) and `<name>@station` for `CBEAM`'s 11 stations and `CBAR`'s station entries (stations not output stay NaN). The grid point force tables (`OGPFB*`) give `GRID_FORCE:F1`...`M3` per element node and `GRID_FORCE:<label>:<M>` point data for the rows of no element (`*TOTALS*`, `APP-LOAD`, `F-OF-SPC`..., named without blanks and asterisks: `GRID_FORCE:TOTALS:F1`), rotated from each GRID's `CD` to basic. MSC writes a different header word 5 in `OGPFB1` than in the other tables of the subcase, so grid point forces join the step the subcase's other tables made instead of starting their own.

**Not read**, and named in one warning: complex grid point forces, heat-transfer temperature rates and loads, nonlinear stress tables (`OESNLX*`), energy loss (`OEDE`), every other element type, and every other table (contact and bolt tables...).

## Verification

The test suite reads eleven files from the [pyNastran](https://github.com/SteveDoyle2/pyNastran) test models (BSD 3-Clause, credited in `tests/python/meshes/nastran_op2/README.md`): MSC and NX output, 32- and 64-bit, SOL 101 statics (among them `static_solid_shell_bar`, a solid, shell and bar model), modes, buckling, transient heat transfer, an NX 2412 file, and one with its geometry tables removed next to its deck. On every step of every file the nodal vectors, the centre stress and strain, and the corner, ply, station and grid point force values equal pyNastran 1.4.1's reading to float precision, and the two engines agree exactly. The nodal vectors and centre values also match on 40 further pyNastran models checked outside the suite (before the corner, ply, station and grid point force values were read). No file from an MSC or NX version newer than those samples, and no big-endian file Nastran wrote, has been read.

Since v16.12.0 three more pyNastran models are in the suite: `freq_elements2` (frequency response, magnitude and phase, every element family's complex forces, stresses and strains, NX's von Mises variants), `modes_complex_elements` (complex modes) and `test_vba` (random: PSD, RMS and NO displacements, accelerations and SPC forces, and acoustic GRIDs). Their complex, random and SORT2 nodal values, and every step's element forces, heat fluxes and strain energies across the suite, equal pyNastran's to float precision (pyNastran keeps complex values in single precision). Left out of that comparison where pyNastran misreads: complex `CSHEAR` forces, `CBAR` station forces (it files them with `CBAR` 34's), complex `CTRIA3` from `OESVM1` and complex composite von Mises. Derived files check the rest: `time_thermal_elements_sort2_only.op2` (the SORT2 tables alone read as their SORT1 twins), `static_elements_bgpdt.op2` and `sol401_tstep1_bgpdt.op2` (`GEOM1` removed: the `BGPDTS`/`EQEXINS` and NX `BGPDT` points equal the GRIDs), and a header-less copy of `static_solid_shell_bar.op2` for `PARAM,POST,-2`. No sample has random element tables, complex SORT2 nodal tables, `CBUSH` or `CGAP` forces, `CBUSH` real stresses or NX solids 300-303.

The coordinate systems are checked on `static_solid_shell_bar_cord.op2`, derived from `static_solid_shell_bar.op2` by pointing every GRID's `CP` and `CD` at the file's own `CORD2R`, `CORD2C` and `CORD2S` systems (`tools/gen_nastran_cord_reference.py`): points match pyNastran's to 1e-12 and displacements to the file's single precision. They are compared with upstream pyNastran's rotation matrices, because 1.4.1 composes cylindrical and spherical ones the wrong way round.

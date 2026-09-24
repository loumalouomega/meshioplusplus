# OpenRadioss time-history files (`<run>T01`, `T02`, …)

The time-history file an OpenRadioss (or Radioss) run writes (`/TFILE` in the engine deck): at each output time, the global energies and momenta, the variables of every part and subset, and the variables each `/TH/<kind>` group of the starter deck asks for (nodes, bricks, shells, springs, rigid bodies, interfaces, sections…). It holds numbers per entity, not a mesh; the model is read from the starter deck by [`radioss`](./radioss.md), and its animation files by [`radioss_anim`](./radioss_anim.md). New in v16.12.0.

| | |
|---|---|
| **Format name** | `radioss_th` |
| **File names** | a stem, then `T` and two digits, no extension (`crashT01`; a restart run writes `crashT02`…); also recognised by content, its title record |
| **Read / Write** | ✓ / — (read-only) |
| **Extra dependencies** | — |

## Reading

```python
import meshioplusplus

first = meshioplusplus.read("crashT01")                  # the first output
last = meshioplusplus.read("crashT01", time_step=-1)     # the last one

last.field_data["meshio:time"]
last.field_data["radioss_th:global:kinetic_energy"]      # one value
last.field_data["radioss_th:part:1:IE"]                  # part 1's internal energy
last.field_data["radioss_th:node:1"]                     # TH group 1 (nodes): (nodes, variables)
last.field_data["radioss_th:node:1:ids"]                 # its node ids
last.field_data["radioss_th:node:1:variables"]           # its variable codes

meshioplusplus.radioss_th.time_values("crashT01")
for time, mesh in meshioplusplus.read_sequence("crashT01"):
    ...
```

```bash
meshioplusplus info crashT01
```

Both engines read the whole format: the C++ core, and a Python reference reader the core falls back to.

## The file

A big-endian Fortran unformatted file (each record framed by its length). First the descriptions:

- a title record (the file version, then 80 characters) and the run's date;
- from version 3051 on, three more records, the last the mass, length and time unit factors;
- the counts of parts, materials, properties, subsets, TH groups and global variables, and the global variables' codes;
- each part (id, title, its variable codes), material and property (id, title), subset (id, child subsets, parts, variable codes, title), and TH group (id, type, entity count, variable count, title; each entity's id and title; the variable codes).

Titles are 40 characters before version 3041, 80 up to 4020 and 100 from 4021. Then one output per time: the time, the global variables, the part variables, the subset variables and one record per TH group (entities × variables).

## Steps and data

The steps are the outputs, in file order. A run stopped while writing leaves an incomplete last output, which is not read. Each step has no points (an empty `vertex` block) and holds field data only:

- **`meshio:time`**: the output's time.
- **Global variables**, one value each, named by their code: `radioss_th:global:internal_energy`, `kinetic_energy`, `x_momentum`, `y_momentum`, `z_momentum`, `mass`, `time_step`, `rotation_energy`, `external_work`, `spring_energy`, `contact_energy`, `hourglass_energy`, `elastic_contact_energy`, `frictional_contact_energy`, `damping_contact_energy`, `plastic_work`, `added_mass`, `percentage_added_mass`, `inlet_mass`, `outlet_mass`, `inlet_energy`, `outlet_energy` (codes 1 to 22; others `var<code>`).
- **Part and subset variables**, one value each: `radioss_th:part:<id>:<name>` and `radioss_th:subset:<id>:<name>`, the name Radioss's keyword for the code (`IE`, `KE`, `XMOM`, `YMOM`, `ZMOM`, `MASS`, `HE`, `TURBKE`, `XCG`…`ZCG`, `XXMOM`…`ZZMOM`, `IXX`…`IZX`, `RIE`, `KERB`, `RKERB`, `RKE`, `ERODED`, `HEAT`, `VX`, `VY`, `VZ`, `PW`; others `var<code>`).
- **TH groups**: `radioss_th:<kind>:<id>`, an (entities × variables) array, with `:ids` (the entity ids) and `:variables` (the variable codes, in column order). The kind comes from the group's type: `node` (0), `brick` (1), `quad` (2), `shell` (3), `truss` (4), `beam` (5), `spring` (6), `sh3n` (7), `sphcel` (51), `inter` (101), `rwall` (102), `rbody` (103), `sectio` (104), `monvol` (107), `accel` (108); others `type<n>`.
- **`radioss_th:unit_factors`**: the mass, length and time factors, in files that have them.

`arrays` keeps only the named field data; `points_only` keeps only the time.

## Not read

- The titles (the run's, the parts', the groups' and their entities') and the date.
- The names of the TH group variables. The file stores their codes; the names (`DX`, `SX`, `FNX`…) are in Radioss's own tables, per group kind, and `th_to_csv` reads them from a `<run>T01_TITLES` file the engine does not write by default.

## Verification

The layout follows OpenRadioss's `th_to_csv` converter (MIT). The test suite reads a T01 that OpenRadioss wrote for a meshio++ deck (two bricks, one end held, the other pushed; see `tests/python/meshes/radioss_th/column/README.md`), and every value of its 8 outputs equals `th_to_csv`'s reading to the 7 digits it prints; the two engines agree exactly. Outside the suite, the 81 T01 files of OpenRadioss's QA tests read the same in both engines and match `th_to_csv` on every file it reads correctly (it miscounts files with a single output). Newer layouts (version 3051 and later, with unit factors and longer titles) are tested on synthetic files only: every OpenRadioss file seen is version 3040.

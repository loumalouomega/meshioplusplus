# OpenRadioss time-history files (`<run>T01`, `T02`, …)

The time-history file an OpenRadioss (or Radioss) run writes (`/TFILE` in the engine deck): at each output time, the global energies and momenta, the variables of every part and subset, and the variables each `/TH/<kind>` group of the starter deck asks for (nodes, bricks, shells, springs, rigid bodies, interfaces, sections…). It holds numbers per entity, not a mesh; the model is read from the starter deck by [`radioss`](./radioss.md), and its animation files by [`radioss_anim`](./radioss_anim.md). Read since v16.12.0; [written](#writing) too.

| | |
|---|---|
| **Format name** | `radioss_th` |
| **File names** | a stem, then `T` and two digits, no extension (`crashT01`; a restart run writes `crashT02`…); also recognised by content, its title record |
| **Read / Write** | ✓ / ✓ (field data only; see [Writing](#writing)) |
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

## Writing

```python
import meshioplusplus

meshioplusplus.write("copyT01", mesh)                               # one output, named by its file name
meshioplusplus.write("history.bin", mesh, file_format="radioss_th")  # or by the format

with meshioplusplus.radioss_th.SeriesWriter("runT01") as w:         # one file, an output per step
    for time, mesh in steps:
        w.write(time, mesh)
meshioplusplus.write_sequence("runT01", steps)                      # the same
```

The writer is the reader's inverse: it takes the `radioss_th:*` field data back to the variable codes, ids and groups, and ignores points and cells, which a time-history file does not hold. `meshio:time` is the output's time (0 when absent). The names are read as the reader writes them:

- `radioss_th:global:<name>`, `radioss_th:part:<id>:<name>` and `radioss_th:subset:<id>:<name>` are one value each; a name is looked up in the reader's own tables, and `var<code>` is the code itself.
- `radioss_th:<kind>:<id>` is a TH group, an (entities × variables) array with its `:ids` and `:variables`. Without them the ids are 1 to n and the variable codes 1 to m. A kind not in the list above is `type<n>`.
- `radioss_th:unit_factors` (three values) makes the file version 4021, with 100-character titles and the unit factors; without it the file is version 3040, the version OpenRadioss writes, with 40-character titles.

Values are stored as big-endian `float32` and ids and codes as `int32`, so a read of a written file returns the values rounded to single precision (exactly, for a file read from OpenRadioss). Anything else in `field_data` has no place in the file and is dropped with a warning that names it. A group array that is not 2-D, an `:ids` or `:variables` of the wrong length, a value that is not an integer in 32-bit range and `radioss_th:unit_factors` without three values are a `WriteError`.

The bytes depend on the names and values only: parts, subsets and TH groups are written sorted by id (groups by type, then id) and their variables by code, whatever the order in `field_data`, so a read of a written file lists the same data but not necessarily in the original file order. The titles are not the run's: the title record holds the meshio++ credit line (the [provenance](../provenance.md) tag, the title record is a bounded slot), the date record `meshio++`, and the parts, subsets, groups and entities are titled `PART <id>`, `SUBSET <id>`, `TH <kind> <id>` and `<kind> <id>`. No materials, properties, child subsets or subset parts are written (the reader reads none of them either), and a subset carries only its variables.

A series writes the descriptions with its first step and one output per step. Every step must have the first step's variables (the same names, ids and variable codes); a step that does not is a `WriteError` naming it, and `{step}` in the path writes one file per step instead. A step's time is the argument of `write` (or the pair's first item in `write_sequence`), not its `meshio:time`. One step is `write`'s bytes. Both engines write the same bytes: the C++ `write_radioss_th` and `RadiossThSeriesWriter` (bound as `_core.radioss_th_write` and `_core.RadiossThSeriesWriter`), and the Python twin the core falls back to when it declines the file. The native `sequence_to_timeseries` and pipeline fan-in take the format too. The framing is the shared Fortran record framer (`detail/fortran_records.hpp` and `_fortran_records.py`), whose reader tolerates a truncated last record here, as a run killed mid-write leaves one.

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

The writer is checked against the hand-framed layout the reader's fixtures were first built with (record by record, framed without the library's framer), by reading its files back through both readers, and by rewriting every output of the OpenRadioss T01 above and reading it back value for value. Files written with version 3040 are the layout `th_to_csv` read in the checks above; `th_to_csv` is not installed in the test environment, so a written file is not run through it by the suite.

# LS-DYNA binary output (`binout`)

The binary file [LS-DYNA](https://lsdyna.ansys.com/) writes in place of its ASCII databases when `*DATABASE_<name>` asks for `BINARY = 2` or `3`: `nodout` (node histories), `glstat` (global energies), `matsum` (per-part energies), `elout` (element histories), `rcforc`, `secforc`, `rbdout` and the others, all in one file (`binout`, continued in `binout%001`...; an MPP run writes one per processor, `binout0000`...). It is an LSDA file (Livermore Software Data Archival), the container lasso-python and LS-PrePost read. New in v16.12.0.

| | |
|---|---|
| **Format name** | `lsdyna_binout` |
| **Extensions** | none: found by its file name, `binout` or `binout0000`... (any case), or by its LSDA header |
| **Read / Write** | ✓ / ✓ (one file, no continuation files; see [Writing](#writing) for what is left, tracked in the [roadmap](../roadmap.md#_3-reader-and-writer-parity), §3.2.3) |
| **Extra dependencies** | none |

## Reading

```python
import meshioplusplus

mesh = meshioplusplus.read("run/binout")                 # the first output
last = meshioplusplus.read("run/binout", time_step=-1)   # the last one

last.points                                  # the nodout nodes where they are now
last.point_data["lsdyna:nid"]                # their ids
last.point_data["displacement"]              # (n, 3); also rotation, velocity, ...
last.field_data["binout:glstat:kinetic_energy"]
last.field_data["binout:matsum:internal_energy"]   # one value per matsum id
last.field_data["binout:matsum:ids"]

meshioplusplus.lsdyna_binout.time_values("run/binout")
```

```bash
meshioplusplus convert run/binout 'nodes_{step}.vtu'   # one point cloud per output
meshioplusplus info run/binout
```

Both engines read the whole format: the C++ core, and a Python reference reader the core falls back to. The writer is [below](#writing).

## The file

An LSDA file starts with an 8-byte header: its size, the sizes of its length, offset, command and type fields, and its byte order. Records are a length and a command; a chain of symbol tables (reached through offset records) holds directory (`CD`) and variable (`VARIABLE`: type, offset, count) entries, and each variable's values sit in a `DATA` record. The continuation files `binout%001`... beside the file are read with it. Variable types are 1-, 2-, 4- and 8-byte integers (signed and unsigned), 4- and 8-byte reals, and links to other variables; 1-byte variables are text (titles, dates) and are not read.

Each database is a directory holding `metadata` (`ids`, `title`, `legend`...) and one folder per output time, `d000001`, `d000002`..., each with its `time` and the variables of that output. `elout` and `jntforc` nest one level deeper (`elout/beam`, `elout/shell`...).

## Steps and data

The steps are `nodout`'s outputs or, in a file without `nodout`, those of the first database with any (in path order). `time_step` picks one; `read_metadata(...).time_values` lists the times; `read_sequence` walks them.

- **Points.** With `nodout`: its nodes (`metadata/ids`, as `point_data["lsdyna:nid"]`) at their current coordinates (`x_coordinate`, `y_coordinate`, `z_coordinate`), one `vertex` cell each. Without it, no points.
- **Point data.** `displacement`, `rotation`, `velocity`, `rotational_velocity`, `acceleration`, `rotational_acceleration` from `nodout`'s `x_`/`y_`/`z_` and `rx_`/`ry_`/`rz_` variables; any other per-node `nodout` variable under its own name.
- **Field data.** `meshio:time`; every other `nodout` variable as `binout:nodout:<variable>`; and every other database's **latest output at or before the step's time** (their times compared in single precision, as the file stores them) as `binout:<database>:<variable>` (`binout:glstat:kinetic_energy`, `binout:elout/beam:axial`), with that output's own time as `binout:<database>:time` and the database's metadata ids as `binout:<database>:ids`. The databases write at their own intervals (`glstat` every second, `nodout` every tenth of one), so a step can carry an older `glstat`. Integers stay integers.

## Writing

```python
import meshioplusplus

meshioplusplus.write("binout", mesh)                                   # one output, named by its file name
meshioplusplus.write("run.bin", mesh, file_format="lsdyna_binout")     # or by the format

with meshioplusplus.lsdyna_binout.SeriesWriter("binout") as w:         # one file, an output per step
    for time, mesh in steps:
        w.write(time, mesh)
meshioplusplus.write_sequence("binout", steps, file_format="lsdyna_binout")   # the same
```

The writer is the reader's inverse. It writes one LSDA file (8-byte lengths and offsets, little-endian, the DATA records first and the one symbol table behind them), with reals stored as `float32` and integers as `int32`, as LS-DYNA writes them. `meshio:time` is the output's time (0 when absent).

- **`nodout`** comes from the mesh. `/nodout/metadata/ids` is `point_data["lsdyna:nid"]` (else 1 to n, and an id outside the 32-bit range is a `WriteError`). Each output is a folder `/nodout/d000001`, `d000002`... holding `time`, the points as `x_coordinate`, `y_coordinate`, `z_coordinate` (a 2-D point set gets `z = 0`), the vectors `displacement`, `rotation`, `velocity`, `rotational_velocity`, `acceleration` and `rotational_acceleration` (each `(n, 3)`, as `x_displacement`, `rx_displacement`...), every other one-dimensional per-point array under its own name, and every `binout:nodout:<variable>` field data as a variable of its own. A per-point array is read back as point data only for two or more nodes: with one node it comes back as field data, and a mesh with one node drops such an array with a warning.
- **Every other database** comes from the field data `binout:<database>:<variable>` (a nested database as `binout:elout/beam:axial`): its folders hold `time` (`binout:<database>:time`, else the step's time) and the variables, and `binout:<database>:ids` is `/<database>/metadata/ids`. A database name that is a step folder (`d000001`) or `metadata`, or a name with an empty part, or a variable name that is empty, over 255 characters, not ASCII or holding a `/` is a `WriteError`.
- A mesh with no point is a binout with no `nodout`: its first database's outputs are the steps, as the reader takes them. A mesh with no point and no `binout:*` field data is a `WriteError`.
- **Not written**, and named in one warning: cells (the points are a `vertex` cloud on the way back), regions, cell data, field data that is not `meshio:time` or under `binout:`, and per-point arrays that are not one value or three per point.

A series writes `nodout`'s ids once and an output folder per step. A database gets a new folder only when its own `binout:<database>:time` changes, which is the inverse of what `read` attaches (the latest output at or before the step's time): a database written every second beside a `nodout` written every tenth keeps its cadence, and `write_sequence(read_sequence(...))` round trips. A database may be absent from a step and appear later; once it has appeared its variables must match its first output's. Every step must have the first step's number of points and ids and `nodout` variables, or the step is a `WriteError` naming it (and `{step}` in the path writes one file per step instead). A step's time is the argument of `write` (or the pair's first item in `write_sequence`), not its `meshio:time`. One step is `write`'s bytes. The file is assembled when the series closes. Both engines write the same bytes: the C++ `write_lsdyna_binout` and `LsdynaBinoutSeriesWriter` and their Python twins, which the test suite compares, and a hand-built file spells out.

Not yet: continuation files (`binout%001`), 4-byte-offset and big-endian files, the databases' text metadata (titles, legends, dates) and links; they are the rest of [roadmap 3.2.3](../roadmap.md#_3-reader-and-writer-parity).

lasso-python 2.0.4 reads the written files (`Binout`): the databases, `nodout`'s ids, times and coordinates, and the `glstat` and `elout/beam` outputs of a series, each database with its own number of outputs. They were not opened in LS-PrePost.

## Not read

- The databases' text metadata (titles, legends, dates, revisions).
- An MPP run's per-processor files together: each `binout0000`... reads on its own.

## Verification

The layout follows lasso-python's `lsda_py3` (BSD 3-Clause), itself LSTC's. The test suite reads a binout LS-DYNA wrote (`glstat` and `rwforc`, from Ansys' [example data](https://github.com/ansys/example-data), MIT) and a subset of another (a Hybrid III dummy's `nodout`, `elout/beam`, `glstat` and `matsum`, rewritten with lasso-python's LSDA writer from the 15 MB original), and every value equals lasso-python's `Binout` reading; the two engines agree exactly. The full 15 MB original (1001 `nodout` outputs, 13 databases) and qd-cae's `swforc` binout were checked outside the suite.

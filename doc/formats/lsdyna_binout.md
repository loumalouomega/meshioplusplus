# LS-DYNA binary output (`binout`)

The binary file [LS-DYNA](https://lsdyna.ansys.com/) writes in place of its ASCII databases when `*DATABASE_<name>` asks for `BINARY = 2` or `3`: `nodout` (node histories), `glstat` (global energies), `matsum` (per-part energies), `elout` (element histories), `rcforc`, `secforc`, `rbdout` and the others, all in one file (`binout`, continued in `binout%001`...; an MPP run writes one per processor, `binout0000`...). It is an LSDA file (Livermore Software Data Archival), the container lasso-python and LS-PrePost read. New in v16.12.0.

| | |
|---|---|
| **Format name** | `lsdyna_binout` |
| **Extensions** | none: found by its file name, `binout` or `binout0000`... (any case), or by its LSDA header |
| **Read / Write** | ✓ / — (read-only) |
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

Both engines read the whole format: the C++ core, and a Python reference reader the core falls back to.

## The file

An LSDA file starts with an 8-byte header: its size, the sizes of its length, offset, command and type fields, and its byte order. Records are a length and a command; a chain of symbol tables (reached through offset records) holds directory (`CD`) and variable (`VARIABLE`: type, offset, count) entries, and each variable's values sit in a `DATA` record. The continuation files `binout%001`... beside the file are read with it. Variable types are 1-, 2-, 4- and 8-byte integers (signed and unsigned), 4- and 8-byte reals, and links to other variables; 1-byte variables are text (titles, dates) and are not read.

Each database is a directory holding `metadata` (`ids`, `title`, `legend`...) and one folder per output time, `d000001`, `d000002`..., each with its `time` and the variables of that output. `elout` and `jntforc` nest one level deeper (`elout/beam`, `elout/shell`...).

## Steps and data

The steps are `nodout`'s outputs or, in a file without `nodout`, those of the first database with any (in path order). `time_step` picks one; `read_metadata(...).time_values` lists the times; `read_sequence` walks them.

- **Points.** With `nodout`: its nodes (`metadata/ids`, as `point_data["lsdyna:nid"]`) at their current coordinates (`x_coordinate`, `y_coordinate`, `z_coordinate`), one `vertex` cell each. Without it, no points.
- **Point data.** `displacement`, `rotation`, `velocity`, `rotational_velocity`, `acceleration`, `rotational_acceleration` from `nodout`'s `x_`/`y_`/`z_` and `rx_`/`ry_`/`rz_` variables; any other per-node `nodout` variable under its own name.
- **Field data.** `meshio:time`; every other `nodout` variable as `binout:nodout:<variable>`; and every other database's **latest output at or before the step's time** (their times compared in single precision, as the file stores them) as `binout:<database>:<variable>` (`binout:glstat:kinetic_energy`, `binout:elout/beam:axial`), with that output's own time as `binout:<database>:time` and the database's metadata ids as `binout:<database>:ids`. The databases write at their own intervals (`glstat` every second, `nodout` every tenth of one), so a step can carry an older `glstat`. Integers stay integers.

## Not read

- The databases' text metadata (titles, legends, dates, revisions).
- An MPP run's per-processor files together: each `binout0000`... reads on its own.

## Verification

The layout follows lasso-python's `lsda_py3` (BSD 3-Clause), itself LSTC's. The test suite reads a binout LS-DYNA wrote (`glstat` and `rwforc`, from Ansys' [example data](https://github.com/ansys/example-data), MIT) and a subset of another (a Hybrid III dummy's `nodout`, `elout/beam`, `glstat` and `matsum`, rewritten with lasso-python's LSDA writer from the 15 MB original), and every value equals lasso-python's `Binout` reading; the two engines agree exactly. The full 15 MB original (1001 `nodout` outputs, 13 databases) and qd-cae's `swforc` binout were checked outside the suite.

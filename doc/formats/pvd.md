# PVD — ParaView collection (`.pvd`)

A time-indexed list of VTK XML files: the on-disk face of the [sequence engine](/sequences) for ParaView, where a `.pvd` gives the time slider its steps (v14.1.0, roadmap §1.1). An entry names a serial or parallel XML file — `.vtu`, `.vtp`, `.vtm`, `.pvtu` or `.pvtp`, never legacy `.vtk` — so `.pvd` → [`.pvtu`](./pvtu.md) → `.vtu` is the ordinary layout of a partitioned transient run. The format is documented on the [ParaView wiki](https://www.paraview.org/Wiki/ParaView/Data_formats) rather than in VTK.

| | |
|---|---|
| **Format name** | `pvd` |
| **Extensions** | `.pvd` |
| **Read / Write** | ✓ / ✓ |
| **Time steps** | ✓ — `time_step=`, `read_metadata()["time_values"]`, [`read_sequence`](/sequences) and [`write_sequence`](/sequences) (fan-in) |
| **Extra dependencies** | — (zlib for compressed pieces, as VTU) |

## Reading & writing

```python
import meshioplusplus

# many steps: one .vtu per step next to the index, the time of each step in the index
meshioplusplus.write_sequence("run.pvd", ((t, mesh_at(t)) for t in times))

mesh = meshioplusplus.read("run.pvd")                    # step 0
last = meshioplusplus.read("run.pvd", time_step=-1)      # the final state
meta = meshioplusplus.read_metadata("run.pvd")
meta["time_values"]                                      # every step's time, off the index alone

for t, mesh in meshioplusplus.read_sequence("run.pvd"):  # one mesh alive at a time
    ...

# one step
meshioplusplus.pvd.write("case.pvd", mesh)               # time from field_data["meshio:time"], else 0
```

`meshioplusplus.read`/`write` dispatch on the extension, so `.pvd` needs no explicit format name.

## File structure

Writing `run.pvd` creates a *sibling directory* named after the index's own stem, holding one `.vtu` per step, numbered with zero padding (width 4):

```
run.pvd
run/
  run_0000.vtu
  run_0001.vtu
  run_0002.vtu
```

```xml
<?xml version="1.0"?>
<VTKFile type="Collection" version="1.0" byte_order="LittleEndian">
<Collection>
<DataSet timestep="0" part="0" file="run/run_0000.vtu"/>
<DataSet timestep="0.5" part="0" file="run/run_0001.vtu"/>
<DataSet timestep="2" part="0" file="run/run_0002.vtu"/>
</Collection>
</VTKFile>
```

`file=` is always relative, written with forward slashes and XML-escaped. Times are written as the shortest spelling that reads back exactly (`0.1`, not `0.10000000000000001`), independent of the process locale.

## Two axes: time and part

A `<DataSet>` carries `timestep`, `part` and (optionally) `group` and `name`. The two index axes map onto the two selectors meshio++ already has:

| `.pvd` attribute | meshio++ selector | default |
|---|---|---|
| `timestep=` | `time_step=` (`ReadOptions::mTimeStep`) | step `0`, the earliest |
| `part=` | `piece=` (`ReadOptions::mPiece`) | unset: every part of the chosen step, merged |
| `group=` / `name=` | region name | — |

- **Steps** are the distinct `timestep` values in ascending order, whatever order the entries are listed in and however unevenly they are spaced. A missing `timestep` is `0`, as ParaView reads it. `time_step=-1` is the last step; out of range names the step count.
- **Parts** of a step are its entries, ordered by `part` and then by document order. `read()` merges them (no welding) with one `RegionKind::Cell` region per entry, named from `name=`, else `group/part_<p>`, else `part_<p>`. `piece=k` keeps one of them alone, with no region.
- A step with a **single entry** is returned as that file reads: a step that is one `.pvtu` keeps its own `piece_<i>` regions and is not wrapped in a second all-cells region.
- The chosen step's time is attached as `field_data["meshio:time"]`, so fan-out and fan-in close on time.
- `read_metadata` reports every step's time from the index alone — no piece is opened for `time_values` — and summarizes step 0's pieces for the rest, which is what a default `read` returns.

```python
meshioplusplus.read("run.pvd", time_step=3)               # step 3, all its parts merged
meshioplusplus.read("run.pvd", time_step=3, piece=2)      # step 3, part 2 alone
```

## Composition: `.pvd` → `.pvtu` → `.vtu`

The normal layout for a partitioned transient run is one `.pvtu` per step, each naming its per-part `.vtu` files, under one `.pvd`:

```python
for k, (t, mesh) in enumerate(steps):
    meshioplusplus.pvtu.write_pieces(f"step{k}.pvtu", meshioplusplus.partition(mesh, 4, ghost_layers=1))
# then a .pvd whose entries are timestep="t" file="step{k}.pvtu"
```

Both levels go through the same piece reader, so they nest with no special case; `ghosts="drop"` (`meshioplusplus.pvd.read(..., ghosts="drop")`) reaches every child. A parallel index does not nest (a `.pvtu` names `.vtu` files only), so no cycle is expressible.

## Quirks & limitations

- **A plain `write` is one step.** `meshioplusplus.write("x.pvd", mesh)` writes a one-entry collection. Many steps go through `write_sequence`, which streams them: one mesh alive at a time, and the index is rewritten after **every** step, so a run that is killed leaves a collection ParaView opens covering every finished step.
- **Step pieces are `.vtu`.** Whatever the mesh, each step is written as a `.vtu` (polygonal data included). The `PvdSeriesWriter` class does not take a piece format.
- **A `TimeValue` field-data array is not read as a time.** The index's `timestep` is the only time source; a piece's own `<FieldData>` is neither written by the `.vtu` writer nor read by the C++ `.vtu` reader.
- **Non-finite times are refused** on write; a `timestep` that is not a number is refused on read.
- An entry whose `file=` is missing, does not exist (naming the attribute and the path — an absolute path from another machine is not searched for elsewhere) or is not one of the five XML extensions raises `ReadError`. An empty collection reads as an empty mesh.
- **Windows-style backslashes** in a `file=` attribute are not translated.
- Opening a `.pvd` in ParaView is not covered by the automated tests: vanilla VTK has no collection reader (`vtkPVDReader` ships with ParaView). What is covered is the index's structure, and that every file it names opens in VTK's own readers.

## Notes

`PvdSeriesWriter` (C++, `formats/pvd.hpp`) is what the native sequence driver (`sequence_to_timeseries`) writes a `.pvd` fan-in through, and so what the C API's `mio_sequence_to_timeseries` and the wrappers built on it reach; that route is covered by the C++ tests, not by tests of each wrapper. The Python `meshioplusplus.pvd.SeriesWriter` is its twin and the one `write_sequence` uses. `write_sequence` to a `.pvd` accepts `Encoding` (ASCII or binary pieces) and rejects `Codec` and `FloatFormat`, the rule every transient writer applies to an option it cannot honour.

## See also

- [Sequences](/sequences) — the engine `.pvd` is a reader and writer for.
- [`.pvtu`](./pvtu.md) / [`.pvtp`](./pvtp.md) — the parallel indices a step may name.
- [VTU](./vtu.md) — what each step is; [VTM](./vtm.md) — the multiblock index; [VTKHDF](./vtkhdf.md) — a transient dataset in one file.

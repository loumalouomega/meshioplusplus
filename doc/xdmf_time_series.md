# XDMF Time Series

XDMF is the only format in meshio++ with built-in support for temporal (time series) data. The mesh topology is written once; field data is written per time step.

Requires `h5py` when using the default `data_format="HDF"`.

---

## Writing a time series

```python
import meshioplusplus

with meshioplusplus.xdmf.TimeSeriesWriter("simulation.xdmf") as writer:
    writer.write_points_cells(points, cells)
    for t, phi in time_steps:
        writer.write_data(t, point_data={"phi": phi})
```

### `TimeSeriesWriter(filename, data_format="HDF")`

| Parameter | Default | Description |
|-----------|---------|-------------|
| `filename` | — | Path to the `.xdmf` file |
| `data_format` | `"HDF"` | `"HDF"` (companion `.h5`), `"XML"` (inline), or `"Binary"` (separate `.bin` files) |

Must be used as a context manager (`with` statement). The `.xdmf` file is written on `__exit__`.

### `writer.write_points_cells(points, cells)`

Write the shared mesh topology. Must be called before `write_data`.

### `writer.write_data(t, point_data=None, cell_data=None)`

Write field data for one time step `t` (a float). Both `point_data` and `cell_data` are dicts of `str -> numpy array`.

---

## Reading a time series

```python
with meshioplusplus.xdmf.TimeSeriesReader("simulation.xdmf") as reader:
    points, cells = reader.read_points_cells()
    for k in range(reader.num_steps):
        t, point_data, cell_data = reader.read_data(k)
```

### `TimeSeriesReader(filename)`

Parses the XDMF file on construction. Only XDMF version 3 is supported for time series.

### `reader.num_steps`

Total number of time steps stored in the file.

### `reader.read_points_cells()`

Returns `(points, cells)` — the shared mesh topology as a numpy array and a list of `CellBlock`.

### `reader.read_data(k)`

Returns `(t, point_data, cell_data)` for time step index `k`.

---

## C++ (`XdmfTimeSeriesWriter`)

The write side also exists natively, in [`meshioplusplus/formats/xdmf_time_series.hpp`](https://github.com/loumalouomega/meshioplusplus/blob/main/src/cpp/include/meshioplusplus/formats/xdmf_time_series.hpp). This is the one writer that is **not** in the shared format registry (`registry.cpp`): a series is a stateful multi-call object, so there is no `(path, mesh)` call for the registry to name.

```cpp
#include "meshioplusplus/formats/xdmf_time_series.hpp"

meshioplusplus::XdmfTimeSeriesWriter w("simulation.xdmf");  // "HDF" by default
w.WritePointsCells(mesh);            // the static grid, once
for (int k = 0; k < nsteps; ++k) {
    solve(mesh);
    w.WriteData(k * dt, mesh);       // only point_data/cell_data are consumed
}
w.Finalize();                        // the destructor would do this too
```

`WritePointsCells` takes the whole `Mesh` and uses its points and cells; `WriteData` takes a `Mesh` and uses its `point_data`/`cell_data`, so a solver can pass the same object it is updating in place. Arrays are written in the uniform API's sorted-name order, which makes a series byte-deterministic across mesh backends. `Finalize()` is idempotent and is called by the destructor; an explicit call exists so a write failure surfaces as an exception rather than being swallowed during unwinding.

::: warning Deleting the output
Because the destructor writes the `.xdmf`, removing that file while the writer is still alive **recreates it**:

```cpp
{
    XdmfTimeSeriesWriter w(path, "XML");
    w.WritePointsCells(m);
    w.WriteData(0.0, {u});
    std::filesystem::remove(path);   // "clean up"
}   // <-- the destructor finalizes here, recreating path
```

The second half is what makes this bite: with `XdmfSeriesMode::Append` the *next* run continues a series you believed deleted, so the symptom appears one run later as a wrong step count rather than at the delete. Destroy the writer before removing its output — end the scope, or call `Finalize()` explicitly and delete afterwards. This matters most in a test suite that writes, asserts and cleans up inside a single scope.
:::

Two differences from the Python writer, both deliberate:

- the `"HDF"` companion is a **sibling of the `.xdmf`** (`<path minus extension>.h5`), not `<stem>.h5` in the process's working directory — the Python spelling breaks the moment the output path has a directory component;
- `"HDF"` on a build without HDF5 support throws by name from the constructor (naming `-DMESHIOPLUSPLUS_WITH_HDF5=ON`) rather than the symbol being absent.

Reading is the ordinary `read_xdmf`: it resolves a temporal collection structurally (it never runs an XInclude/XPointer pass), takes the geometry from the static grid and the attributes from the step [`ReadOptions::mTimeStep`](selective_read.md) selects — `0` first, `-1` last, out of range an error naming the count. `read_metadata` reports every step's `<Time Value>` in `mTimeValues` from the XML alone, without touching a payload.

## C API

The same handle, flat:

```c
mio_xdmf_series* s = mio_xdmf_series_create("simulation.xdmf", "HDF", -1);
mio_xdmf_series_write_points_cells(s, mesh);
for (k = 0; k < nsteps; ++k)
    mio_xdmf_series_write_data(s, k * dt, mesh);
mio_xdmf_series_free(s);   /* finalizes if mio_xdmf_series_finalize wasn't called */
```

`mio_xdmf_series_num_steps` reports the count written so far. `mio_xdmf_series_create` returns `NULL` (with `mio_last_error()` set) for an unknown data format or for `"HDF"` without HDF5 support. The three languages riding this C API — Fortran, Julia and R — each expose the same handle in their own idiom, below.

One rule is shared by all three, and is why `finalize` exists as a call of its own rather than only as part of the free: **a write failure during the implicit finalize cannot be reported** — from a Fortran `free()`, a Julia `close`, or an R garbage-collection finalizer alike. Call `finalize` explicitly when you want to see one. The other shared rule is that the `.xdmf` light data is buffered until then, so a series is only readable once it has been finalized.

## Fortran

`type(mio_xdmf_series)`, with the module's usual conventions: type-bound procedures, optional `stat`/`errmsg` on every fallible one, and handles freed explicitly (no finalizer) exactly like `type(mio_mesh)`. See [Fortran](./fortran.md#transient-time-series-xdmf-writing).

```fortran
type(mio_xdmf_series) :: series

call series%create('simulation.xdmf')      ! "HDF" by default
call series%write_points_cells(m)          ! the static grid, once
do k = 0, nsteps - 1
    call solve(m)
    call series%write_data(k*dt, m)        ! point_data/cell_data only
end do
print *, series%num_steps()
call series%finalize()                     ! free() would do this too
call series%free()
```

`create` takes optional `data_format` and `gzip_level`. An unknown format, or `'HDF'` on a build without HDF5, fails through `stat`/`errmsg` rather than aborting.

## Julia

`XdmfSeries`, released by a GC finalizer with a deterministic, idempotent `close` — the [`Mesh`](./julia.md#memory-management) rule. See [Julia](./julia.md#transient-time-series-xdmf-writing).

```julia
s = XdmfSeries("simulation.xdmf")          # "HDF" by default
write_points_cells!(s, mesh)
for k in 0:nsteps-1
    solve!(mesh)
    write_data!(s, k * dt, mesh)
end
num_steps(s)
finalize!(s)                               # close(s) would do this too
close(s)
```

There is also a `do`-block form, `XdmfSeries(path) do s ... end`, which closes the series even if the body throws. The exported name is `finalize!`, not `finalize`: `Base.finalize` runs an object's GC finalizer and means something quite different.

## R

`mio_xdmf_series()`, an external pointer with its **own** tag — so a `mio_mesh` and a series are never accepted for one another — and a registered finalizer. See [R](./r.md#transient-time-series-xdmf-writing).

```r
s <- mio_xdmf_series("simulation.xdmf")        # "HDF" by default
mio_xdmf_series_write_points_cells(s, mesh)
for (k in 0:9) {
  mio_xdmf_series_write_data(s, k * 0.1, mesh)
}
mio_xdmf_series_num_steps(s)
mio_xdmf_series_finalize(s)                    # release() would do this too
mio_xdmf_series_release(s)
```

`mio_xdmf_series_release()` is the idempotent deterministic free and `mio_xdmf_series_is_open()` the predicate. Release the handle before the directory it writes into goes away: a finalizer running later would try to write the `.xdmf` into a path that no longer exists, and cannot report it.

## Python (the C++ writer, explicitly)

The C++ writer is reachable from Python as `_core.XdmfTimeSeriesWriter`. It is **not** wired underneath `meshioplusplus.xdmf.TimeSeriesWriter`, which keeps its own documented behaviour untouched — the two are not drop-in equivalents, so the C++ one is exposed *additionally* and by name (the same choice made for [`.mdpa`](./formats/mdpa.md#c-core)).

```python
from meshioplusplus import _core

with _core.XdmfTimeSeriesWriter("simulation.xdmf") as w:   # "HDF" by default
    w.write_points_cells(mesh)                             # the static grid, once
    for k in range(nsteps):
        solve(mesh)
        w.write_data(k * dt, mesh)      # only point_data/cell_data are used
# the .xdmf is written on __exit__
```

`XdmfTimeSeriesWriter(path, data_format="HDF", gzip_level=-1)` has `write_points_cells(mesh)`, `write_data(time, mesh)`, `finalize()`, the read-only properties `num_steps` and `finalized`, and `__enter__`/`__exit__` (which finalizes, as the Python writer's does — it never suppresses the body's exception).

Which one to reach for:

| | `meshioplusplus.xdmf.TimeSeriesWriter` | `_core.XdmfTimeSeriesWriter` |
|---|---|---|
| Arguments | raw `points, cells` / `point_data=`, `cell_data=` | whole `Mesh` objects |
| `"HDF"` companion | `<stem>.h5` **relative to the CWD** | `<path minus extension>.h5`, a **sibling** of the `.xdmf` |
| Array order | dict insertion order | the uniform API's sorted-name order |
| Needs `h5py` | yes, for `"HDF"` | no — the core's own HDF5 |
| Available without a compiled core | yes | no |

Use the C++ one when you already hold a `Mesh` (a solver updating one in place is the motivating case), when the output path has a directory component, or when you want the write to cost no Python per step. Use the Python one when you have loose numpy arrays, want `h5py` interop, or are running against a build with no compiled core.

Reading is unchanged either way: `meshioplusplus.xdmf.TimeSeriesReader`, or the ordinary `meshioplusplus.read` (which resolves the temporal collection structurally and materializes one step), with `meshioplusplus.read_metadata(...)["time_values"]` to enumerate the steps.

`"HDF"` on a build without HDF5 raises `meshioplusplus.WriteError` from the constructor, naming `-DMESHIOPLUSPLUS_WITH_HDF5=ON` — a clean Python exception, not a missing symbol.

## WASM

The same handle, in JavaScript — the **one stateful** binding in [`@meshioplusplus/wasm`](./wasm.md), since every other one is a pure function over a mesh object:

```javascript
const w = m.createXdmfTimeSeriesWriter('/series.xdmf');  // 'HDF' by default
w.writePointsCells(mesh);
for (let k = 0; k < nsteps; ++k) w.writeData(k * dt, stepMesh(k));
w.close();                                  // the files appear HERE

const xdmf = m.FS.readFile('/series.xdmf');
const h5 = m.FS.readFile('/series.h5');     // 'HDF' writes TWO files
```

The object has `writePointsCells(mesh)`, `writeData(time, mesh)`, `finalize()`, `numSteps()`, `finalized()` and `close()`. **`'HDF'` puts two files in the virtual filesystem** — the `.xdmf` and its sibling `.h5` — and nothing at all until `finalize()`/`close()` runs; copy both out of `FS`. `'HDF'` works in the shipped artifact, which links a wasm32 HDF5.

Under the ergonomic wrapper the raw binding is an opaque integer handle plus seven free functions rather than an embind `class_`, so that no live C++ object (and no Emscripten `.delete()`) reaches JS and every error arrives as a readable `Error`; see [WASM § Transient XDMF](./wasm.md#transient-time-series-xdmf).

---

## Notes

- The mesh topology is stored once in the XDMF file and referenced by each time step using XInclude.
- With `data_format="HDF"`, all numerical data goes into a companion `<stem>.h5` file. Both files must be present to read.
- `data_format="XML"` embeds all data directly into the XML, which avoids external files but produces large `.xdmf` files.
- `data_format="Binary"` writes one `.bin` file per data array; useful when HDF5 is not available.

## Crash resilience: `Flush()` and append mode

By default the light data (`.xdmf`) is written once, at `Finalize()`: the temporal collection element has to enclose every step. That is fine for a script that runs to completion and wrong for a solve that is killed, hits a node failure, or is simply still running — those leave heavy data on disk and no readable `.xdmf` at all.

`Flush()` writes the document as it currently stands without finalizing, so the file opens in ParaView and covers every step written so far:

```cpp
meshioplusplus::XdmfTimeSeriesWriter w("run.xdmf");
w.WritePointsCells(mesh);
for (int k = 0; k < nsteps; ++k) {
    solve(mesh);
    w.WriteData(k * dt, mesh);
    if (k % 10 == 0)
        w.Flush();          // a checkpoint every ten steps
}
w.Finalize();
```

The document goes to a sibling temp file and is then `rename`d over the target, so a crash *during* a flush cannot truncate the previous one; heavy data is flushed first, so the `.xdmf` never names a dataset that is not on disk yet.

`SetAutoFlush(true)` does it after every `WriteData`. It is **off by default** because a flush re-serializes the whole document, so flushing every step is quadratic in the step count — and for `"XML"`, whose heavy data lives *in* that document, quadratic in the data volume too. Pick a cadence that matches what you can afford to lose.

### Continuing a series

`XdmfSeriesMode::Append` continues the temporal collection already at the path instead of overwriting it, which is what a restarted analysis needs:

```cpp
meshioplusplus::XdmfTimeSeriesWriter w("run.xdmf", "HDF", -1,
                                       meshioplusplus::XdmfSeriesMode::Append);
// w.NumSteps() is already the count the earlier run wrote.
```

Appending to a path with **no file yet is not an error** — it is simply a fresh series, so a restartable solver can pass `Append` unconditionally instead of probing the filesystem. `WritePointsCells` still rejects a second call: the static grid is exactly the thing that must not be duplicated, and an appended series reuses the one already in the file.

The heavy-data counter resumes past what is already on disk by **scanning the container** — the `.h5` root group, or probing `<base>N.bin` — rather than trusting the document, because a mis-resumed counter would silently overwrite `data0` instead of failing.

The point and cell counts are recovered from the document too, out of `<Topology NumberOfElements>` and the geometry `<DataItem>`'s `Dimensions`, so an appended series validates the [`NamedArray` overload](#writing-a-step-from-solver-arrays) exactly as a fresh one does. **In v9.1.0 it did not**: the counts stayed at 0 and that overload rejected every array with `expected 0 (0 x 1)`, so a resumed solver had to fall back to passing a whole `Mesh` per step. If a foreign document declares neither, meshio++ warns once and skips the length check rather than failing — the `<DataItem>` carries its own `Dimensions`, so the output is still valid.

The collection and the static grid are located **structurally**, by the same resolver `read_xdmf` uses: the first `GridType="Collection" CollectionType="Temporal"` under `<Domain>`, and any `GridType="Uniform"` sibling. A grid named something other than `mesh` is therefore recognized (v9.1.0 matched the literal name and would append a second static grid), and a non-version-3 document is rejected rather than silently continued.

## Writing a step from solver arrays

After `WritePointsCells` the geometry and connectivity are fixed, and each step produces nodal and elemental *values*, not a new mesh. The `NamedArray` overload takes them directly:

```cpp
meshioplusplus::XdmfTimeSeriesWriter::NamedArray u;
u.mName = "displacement";
u.mNumComponents = 3;                     // row-major, NumPoints * 3 values
u.mValues = solver.NodalDisplacements();
w.WriteData(t, {u});
```

Under the KRATOS mesh backend this is not merely tidier: passing a `Mesh` per step re-stages the whole `ModelPart` every output step when only the values changed.

Row counts are validated against the static grid and a mismatch throws naming the array — deliberately stricter than the `Mesh` overload, which reads the counts off the mesh and cannot get them wrong. Arrays are emitted in the order given, unlike the `Mesh` overload's sorted-name order.

### On other surfaces

| Surface | Flush | Append | Solver arrays |
|---|---|---|---|
| C++ | `Flush()`, `SetAutoFlush()` | `XdmfSeriesMode::Append` | `WriteData(t, point, cell)` |
| C API | `mio_xdmf_series_flush` | `mio_xdmf_series_create_ex` + `mio_xdmf_series_opts.mode` | `mio_xdmf_series_write_data_arrays` + `mio_named_array` |
| Fortran | `s%flush()` | `s%create(..., mode='append')` | — (see below) |
| Julia | `flush!(s)` | `XdmfSeries(path; mode=:append)` | `write_data!(s, t, Dict(...))` |
| R | `mio_xdmf_series_flush()` | `mio_xdmf_series(..., mode = "append")` | — |
| WASM | `w.flush()` | `{ mode: 'append' }` | `w.writeDataArrays(t, {...})` |
| Python | `w.flush()`, `w.auto_flush` | `mode="append"` | `w.write_data_arrays(t, {...})` |

**Fortran has no solver-array overload**, deliberately: an array of derived types holding interop pointers is a poor fit for Fortran, and a Fortran solver already holds an `mio_mesh` handle it can `add_point_data` into before `write_data`.

## See also

- [Sequences](sequences.md) — the driver that fans a directory of single-step files **in** to a series through this writer, and fans a series **out** to one file per step. It streams (one mesh alive at a time), so a 500-step run needs no more memory than a one-step one.

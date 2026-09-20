# VTKHDF Time Series

VTKHDF stores a transient case as **one static grid plus one appended step per solve**: the geometry is written once, every field is a chunked, extendable dataset, and the `Steps` group is a table of offsets into those flat arrays. This is the on-disk form of the [sequence engine](./sequences.md): a run that produces `out_0000.vtu … out_0500.vtu` becomes one file with random access per step, and it opens in ParaView with a working time slider. See also the [format page](./formats/vtkhdf.md).

Two writers produce the **same layout**: a pure-Python one (`meshioplusplus.vtkhdf.TimeSeriesWriter`, h5py) and a C++ one (`meshioplusplus._core.VtkhdfTimeSeriesWriter`, `VtkhdfTimeSeriesWriter` in C++). Neither is buffered.

---

## Writing a time series

```python
import meshioplusplus

with meshioplusplus.vtkhdf.TimeSeriesWriter("simulation.vtkhdf") as writer:
    writer.write_points_cells(points, cells)
    for t, u in time_steps:
        writer.write_data(t, point_data={"u": u})
```

### `TimeSeriesWriter(filename, compression=None, compression_opts=4, mode="truncate")`

| Parameter | Default | Description |
|-----------|---------|-------------|
| `filename` | — | Path to the `.vtkhdf` file (a path, not a buffer) |
| `compression` | `None` | `"gzip"` or `None`; applies to the appended datasets |
| `compression_opts` | `4` | gzip level 0-9 |
| `mode` | `"truncate"` | `"truncate"` starts a fresh series; `"append"` continues the one at `filename`. A path that does not exist yet is just a fresh series, so a solver can pass `"append"` unconditionally |

Must be used as a context manager; `__exit__` finalizes.

### `writer.write_points_cells(points, cells)`

Write the shared grid. Once, before `write_data`. A polyhedral grid is supported and bumps the declared version to 2.5. When continuing an existing series this writes nothing and only checks that the point and cell counts match, so a driver can call it unconditionally.

### `writer.write_data(t, point_data=None, cell_data=None, field_data=None)`

Append one step at simulation time `t`. `cell_data` values may be one array or a per-block list (concatenated in block order). Arrays keep their dtype. `field_data` may vary in length from step to step, but not in component count. The reserved `meshio:time` field is skipped, since `t` already records it.

### The set of array names is fixed at the first step

A later step that introduces a new name, drops one, changes a dtype or a component count, or has the wrong number of rows raises `WriteError` naming the array. VTK indexes `PointDataOffsets/<name>` by step, so a short table would be a silent truncation; a refused step leaves the file exactly as it was.

### Durability

Every step lands in the file as it is written, and `NSteps` is rewritten after each — so a run that is killed leaves a file ParaView opens, covering every completed step. That is the difference from the [XDMF writer](./xdmf_time_series.md), whose `.xdmf` light data appears only when the writer finalizes: here a flush is a cheap `H5Fflush`, not a re-serialization of a growing document, so auto-flush defaults to **on** (`writer.auto_flush = False` turns it off).

Durable means the file is complete and readable once the writer closes **or the process dies** — a killed run needs no recovery step. It does not mean a *second process* can read it while the writer is still open: HDF5 takes a file lock for writing, so a concurrent reader (a ParaView session tailing a live run) needs `HDF5_USE_FILE_LOCKING=FALSE`, or to wait. Within one Python process the two engines can also collide for the same reason, because h5py and the C++ core each bring their own libhdf5; the C++ reader then declines and the shim falls back to h5py.

---

## Reading a time series

```python
with meshioplusplus.vtkhdf.TimeSeriesReader("simulation.vtkhdf") as reader:
    reader.num_steps                    # 12, from Steps/Values
    reader.times                        # [0.0, 0.5, ...]
    points, cells = reader.read_points_cells()
    t, point_data, cell_data = reader.read_data(k)
```

`read_data(k)` accepts a negative `k`. The generic entry points work too: `meshioplusplus.read(path, time_step=k)`, `read_metadata(path)["time_values"]`, and the [sequence engine](./sequences.md)'s `read_sequence`/`TimeSeries`.

---

## The C++ writer

```cpp
#include "meshioplusplus/formats/vtkhdf_time_series.hpp"

meshioplusplus::VtkhdfTimeSeriesWriter w("out.vtkhdf");
w.WritePointsCells(mesh);
for (int k = 0; k < nsteps; ++k) {
    solve(mesh);
    w.WriteData(k * dt, mesh);   // point_data, cell_data and field_data are consumed
}
w.Finalize();                    // the destructor would do this too
```

`WriteData(time, const std::vector<NamedArray>& point, const std::vector<NamedArray>& cell = {})` takes raw solver arrays with no `Mesh` in between. The class is move-only, and a moved-from writer answers its observers as for a finished series while `WritePointsCells` and both `WriteData` overloads throw. `VtkhdfSeriesMode::Append` continues an existing series.

From Python the same class is `meshioplusplus._core.VtkhdfTimeSeriesWriter`, taking whole `Mesh` objects:

```python
with meshioplusplus._core.VtkhdfTimeSeriesWriter("out.vtkhdf") as w:
    w.write_points_cells(mesh)
    w.write_data(0.5, mesh)
    w.write_data_arrays(1.0, {"u": u})   # name -> array, no Mesh
```

It is deliberately **not** wired underneath the Python class: the Python writer's documented API takes raw arrays where this one takes a `Mesh`. The sequence engine prefers the C++ writer when it is available and falls back to the Python one otherwise.

---

## On-disk layout

```
VTKHDF/                    Type="UnstructuredGrid"  Version=[2,0]
  NumberOfPoints           [n]              length 1, written once
  NumberOfCells            [m]
  NumberOfConnectivityIds  [c]
  Points, Connectivity, Offsets, Types      written once
  PointData/u              (N*n, ...)       chunked, extendable, one step's rows appended per step
  CellData/p               (N*m, ...)
  FieldData/gain           concatenated over steps
  Steps/                   NSteps=N
    Values                 (N,)
    PartOffsets, PointOffsets, CellOffsets, ConnectivityIdOffsets   all zero: every step points at the same geometry
    NumberOfParts          (N,) all 1
    PointDataOffsets/u     (N,) = [0, n, 2n, …]      CellDataOffsets/p (N,)
    FieldDataOffsets/gain  (N,)                       FieldDataSizes/gain (N, 2) = [ncomp, ntuples]
```

With polyhedra, `FaceConnectivityOffsets`, `FaceOffsetsOffsets` and `PolyhedronToFaceIdOffsets` are added, all zero.

---

## Notes

The two writers are checked against each other on this layout (dataset names, shapes, dtypes and values), and against `vtkHDFReader` itself: `TIME_STEPS` matches the step times, the geometry is identical at every step, and every array is exact — including a transient polyhedral mesh.

A transient *partitioned* file written by other tools is read (see the [format page](./formats/vtkhdf.md#partitions)), but these writers produce a single static grid; a partitioned transient case can be assembled step by step through the sequence engine or read back per piece with `piece=`.

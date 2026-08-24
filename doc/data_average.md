# Point ↔ cell data averaging

`point_data_to_cell_data` and `cell_data_to_point_data` move data between locations without changing the geometry. They are [data operations](/data_operations), not file formats.

- **point → cell**: each cell's value is the mean over its own nodes.
- **cell → point**: each point's value is the mean over the cells incident to it, optionally weighted by each cell's |measure|.

![Point -> cell -> point averaging round-trip, same colour scale both sides](/images/data_average_roundtrip.png)

```python
import meshioplusplus as mp

mesh = mp.read("part.vtu")

# Every point_data array becomes a cell_data array of the same name.
cells = mp.point_data_to_cell_data(mesh)

# Only the named arrays, written under a suffixed name.
cells = mp.point_data_to_cell_data(mesh, keys=["T", "p"], suffix="_c")

# Cell -> point, weighted by cell area/volume so large cells count for more.
points = mp.cell_data_to_point_data(mesh, keys=["stress"], weighted=True)
```

## Weighting

`weighted=False` (the default) counts every incident cell equally. `weighted=True` weights each contribution by the cell's |measure| — length for a 1D cell, area for a 2D cell, volume for a 3D one. On a mesh with cells of very different sizes the two differ substantially; on a uniform mesh they coincide.

Cells whose measure cannot be computed — ragged polygon blocks; a polyhedron block's volume **is** computable since v9.16.0 — fall back to a unit weight, with one warning per call.

## Components and data types

Both directions work **component-wise**, so scalar, vector and tensor arrays are handled identically. A scalar array stays 1-D; a vector array keeps its `(n, ncomp)` shape.

The output is **always `float64`**, whatever the input dtype: the mean of an `int32` field is not an `int32`.

## Edge cases

- A point touched by no cell — or by no cell carrying a finite value — yields `NaN`. So does a cell whose nodes all carry non-finite values.
- Non-finite values never contribute to an average. See the [NaN policy](/data_operations#nan-policy).
- Ragged polygon blocks average over the row's nodes; polyhedron blocks average over the **distinct** nodes across all of a cell's faces.

::: tip Determinism
The cell → point direction is a scatter, and floating-point addition is not associative, so its accumulation pass runs **serially** on purpose: the result must not depend on the thread count. Only the final divide is parallelised.
:::

## CLI

```bash
meshioplusplus data to-cell  in.vtu out.vtu --keys T,p --target-suffix _c
meshioplusplus data to-point in.vtu out.vtu --keys stress --weighted
```

Omit `--keys` to convert every array at the source location. See the [CLI reference](/cli#meshioplusplus-data).

## Other languages

- **C API** — `mio_data_point_to_cell(mesh, names, count, suffix)` and `mio_data_cell_to_point(mesh, names, count, weight, suffix)`, where `weight` is `MIO_WEIGHT_UNIFORM` or `MIO_WEIGHT_MEASURE`. Pass `count == 0` for "every array". See the [C API reference](/c_api).
- **Fortran** — `m%data_point_to_cell()` and `m%data_cell_to_point(["stress"], weight=MIO_WEIGHT_MEASURE)`. See the [Fortran reference](/fortran).
- **WebAssembly / JavaScript** — `dataPointToCell(mesh, ["T"], "_c")` and `dataCellToPoint(mesh, ["stress"], "measure", "")`. See the [WebAssembly reference](/wasm).

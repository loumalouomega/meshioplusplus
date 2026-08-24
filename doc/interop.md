# Interoperability

Convert a meshio++ `Mesh` to and from the in-memory objects of the tools you reach for next — **without a file round-trip** — sharing the underlying numpy buffers wherever the target accepts them as they are.

```python
import meshioplusplus as mio

mesh = mio.read("bracket.msh")
grid = mio.to_pyvista(mesh)          # a pyvista.UnstructuredGrid, buffers shared
grid.plot(scalars="T")
```

Five targets ship today: **PyVista**, **trimesh**, **Apache Arrow / Parquet**, **pandas** and **polars**. Open3D and DOLFINx are [Phase 2](#phase-2-open3d-and-dolfinx). The same story continues on the [GPU handoff page](./gpu) — DLPack export and CuPy transfer — and the [ML data handling page](./ml) — graphs, feature matrices, dataset export and PyTorch/JAX tensors — both built on this module's payload conventions.

Nothing here is part of the C++ core, which stays dependency-free. This is pure Python over the numpy the readers already return: no format is added, no binding is touched, and every third-party import is lazy.

## Installation

| Extra | Brings | For |
|---|---|---|
| `meshioplusplus[pyvista]` | [PyVista](https://pyvista.org/) (BSD-3-Clause) | `to_pyvista` / `from_pyvista` |
| `meshioplusplus[trimesh]` | [trimesh](https://trimesh.org/) (MIT) | `to_trimesh` / `from_trimesh` |
| `meshioplusplus[arrow]` | [pyarrow](https://arrow.apache.org/) (Apache-2.0) | `to_arrow` / `write_parquet` / … |
| `meshioplusplus[pandas]` | [pandas](https://pandas.pydata.org/) (BSD-3-Clause) | `to_pandas` |
| `meshioplusplus[polars]` | [polars](https://pola.rs/) (MIT) | `to_polars` |
| `meshioplusplus[interop]` | all five | everything on this page |

These are **not** in `meshioplusplus[all]`. `[all]` means "the optional dependencies the *formats* need" — h5py and netCDF4 — and quietly adding a PyVista-sized dependency there would surprise anyone who installs it to read an XDMF file. `[viewer]`, `[kahip]` and `[codecs]` stand apart for the same reason.

Each entry point imports its target *inside* the function and raises a named error when it is missing, so an install without the extra behaves exactly as before:

```
ImportError: meshio++: to_pyvista: pyvista is not installed;
             install it with `pip install meshioplusplus[pyvista]`
```

`has_pyvista()`, `has_trimesh()`, `has_arrow()`, `has_pandas()`, `has_polars()`, `has_open3d()` and `has_dolfinx()` answer the same question without raising.

## Architecture: the pure payload layer

`src/python/meshioplusplus/_interop.py` is split exactly the way [`_viewer.py`](./viewer) is, and for the same reason. The bulk of it is a **pure payload layer** — `_to_vtk_payload`, `_to_triangles_payload`, `_to_table_payload` — which imports no third-party library at all, does not mutate its input, and returns plain numpy plus frozen dataclasses. The public `to_*` / `from_*` functions are thin wrappers over it.

That split is what makes the genuinely subtle parts testable with none of the optional libraries installed: block-major cell indexing, which arrays are shared and which are copied, and the routing of non-triangle input through existing operations. `tests/python/test_interop.py` runs that half in every CI leg; the gated half runs in the dedicated `interop` job.

Two things are **reused rather than transcribed**: the meshio↔VTK cell-type map and node orderings come from `_vtk_common.py` (the same tables the VTU/VTK writers use), and the global block-major cell index comes from `_regions.block_bases` — the Python owner of the core's `detail/cell_index.hpp` numbering. Nothing here re-derives a block offset.

## The zero-copy contract

A buffer is **shared** when the target accepts the array as-is: contiguous, supported dtype, right shape. Every `to_*` takes `zero_copy_only` — with one deliberate exception: [`to_polars`](#pandas-and-polars) takes none, because polars always copies into its own Arrow-backed buffers, and a flag that could only ever raise would be dishonest:

- `zero_copy_only=False` (default) — copies happen, and each is recorded in a `notes` list surfaced as a warning. Nothing is lost silently.
- `zero_copy_only=True` — any step that would copy raises `InteropError` (a `ValueError`), naming the array and the reason. pyarrow's own `zero_copy_only` is the precedent.

```python
grid = mio.to_pyvista(mesh, zero_copy_only=True)
# InteropError: meshio++: to_pyvista: connectivity cannot be shared with the
#               target: connectivity dtype int32 is not int64 (VTK's vtkIdType)
```

### What the flag governs

It governs arrays that **exist in the `Mesh`**: `points`, each `point_data` array, each `cell_data` array, and — for a single-block mesh — connectivity.

It does **not** govern *derived* arrays. VTK's `offsets` and `celltypes` have no meshio++ counterpart at all, and a multi-block mesh has no single connectivity array to share, so those are always constructed. Making them fatal would mean the flag rejected every mixed-type mesh, which is the normal case rather than an edge case. Derived arrays are listed on the payload's `derived` field but are not warnings — they are not a loss, just a difference in data model.

### What forces a copy

| Cause | Where |
|---|---|
| connectivity is not `int64` (VTK's `vtkIdType`) | PyVista |
| a `wedge` block — meshio's node order (the gmsh Prism) differs from VTK's | PyVista |
| points are not `(P, 3)`; a 2-D mesh is padded with zeros | both |
| points are not `float64` | trimesh |
| any array that is not C-contiguous, or is `bool` (widened to `uint8`) | all |
| `cell_data` spanning more than one block (concatenated block-major) | derived |
| the `x`/`y`/`z` columns — strided slices of a row-major `points` array | derived |

Arrays that cannot be represented at all — `object` dtype, complex, strings — are **dropped with a note**, never an exception.

### Lifetime

meshio++'s readers hand back capsule-backed numpy owned by the C++ core, so a wrapper that shares a buffer without holding a reference is a use-after-free rather than a wrong answer. Every returned wrapper carries a `_meshioplusplus_refs` list holding each shared array, so it stays valid after the source mesh is garbage-collected. (VTK's own `numpy_to_vtk(deep=0)` also attaches a `_numpy_reference`, but this does not rely on that.) A test converts, deletes the source `Mesh`, forces `gc.collect()`, churns the allocator, and re-reads the values.

## PyVista

```python
grid = mio.to_pyvista(mesh)
mesh2 = mio.from_pyvista(grid)
```

PyVista ships its own `from_meshio`, but it targets the upstream `meshio` package and does not recognize `meshioplusplus` — hence these functions.

The blocks are turned into the VTK 9 `connectivity` / `offsets` / `celltypes` triple in block order, so **mixed-type meshes are the normal case**. `point_data` and `cell_data` land on `grid.point_data` / `grid.cell_data`, the latter concatenated block-major.

| meshio++ | PyVista | Notes |
|---|---|---|
| `points` | `grid.points` | shared for `float32`/`float64` `(P, 3)`; a 2-D mesh is padded |
| cell blocks | `connectivity`/`offsets`/`celltypes` | in block order; `wedge` is reordered |
| ragged `polygon` blocks | `VTK_POLYGON` with per-row offsets | supported |
| `polyhedron`, `custom`, unmapped types | — | **dropped**, warned by name |
| `point_data` / `cell_data` | `grid.point_data` / `grid.cell_data` | shared where dtype/contiguity allow |
| `Point` / `Cell` regions | `region:<name>` int8 masks + a sidecar | see below |
| `Side` regions | — | **dropped**: VTK has no `(cell, local facet)` concept |
| `field_data`, `mesh.info`, `gmsh_periodic` | — | not carried |

### Regions

`Point` and `Cell` regions become **int8 0/1 mask arrays** named `region:<name>` on `grid.point_data` / `grid.cell_data`, which is directly colourable in ParaView. A mask alone cannot express a region's `dim` and `tag`, so a JSON sidecar rides in `grid.field_data["meshioplusplus:regions"]`:

```json
[{"name": "inlet", "kind": "point", "dim": 2, "tag": 11},
 {"name": "steel", "kind": "cell",  "dim": 3, "tag": 7}]
```

`from_pyvista` reconstructs regions exactly when the sidecar is present, so a gmsh physical group's integer tag survives the round-trip. A grid built by something else carries the masks but not the sidecar; its regions come back with `dim = tag = -1` and a warning.

A cell region entry landing in a **dropped** block is removed from the mask, with a note — the global block-major index advances over skipped blocks, so nothing silently lands on the wrong cell.

`Side` regions are the one asymmetry: VTK has no `(cell, local facet)` concept, so they are dropped with a warning naming them, and `from_pyvista` cannot invent them back.

## trimesh

```python
tm = mio.to_trimesh(mesh)
mesh2 = mio.from_trimesh(tm)
```

trimesh is **triangles only**. Non-triangle input is routed through meshio++'s own operations rather than a reimplementation, in this order, each warned about:

1. volume cells → [`extract_surface`](./extract_surface)
2. higher-order 2-D cells → [`convert_cells("linearize")`](./convert_cells)
3. quads and polygons → `convert_cells("simplexify")`
4. anything still not a triangle (`line`, `vertex`) → dropped with a note

`payload.ops` records exactly which of these ran.

| meshio++ | trimesh | Notes |
|---|---|---|
| `points` | `tm.vertices` | shared for `float64` `(P, 3)` |
| triangle blocks | `tm.faces` | shared for a single `int64` block |
| `point_data` | `tm.vertex_attributes` | |
| `cell_data` | `tm.face_attributes` | for the surviving triangle blocks |
| `color`/`rgb`/`*_color`/`*_rgb` | `trimesh.visual.ColorVisuals` | matched by **name**, never inferred from a `[0, 1]` range |
| `Point` / `Cell` regions | `region:<name>` masks + `tm.metadata` sidecar | only what the composed operations preserve |
| `Side` regions | — | dropped |

`trimesh.Trimesh(..., process=False)` is used deliberately: the default `process=True` merges coincident vertices, which silently renumbers the mesh and would invalidate every index-based attribute.

**Regions follow the composed operations' documented behaviour and get no second policy here.** In practice: a pure-triangle input keeps its `Point`/`Cell` regions; anything routed through `extract_surface` loses them, because `extract_surface` drops regions (its output cells and points are newly created and have no correspondence with the input's). `convert_cells` remaps them.

## Arrow and Parquet

```python
table = mio.to_arrow(mesh, location="cell")   # a pyarrow.Table
mio.write_parquet(mesh, "data.parquet", location="point")
arrays = mio.read_parquet("data.parquet")     # {name: ndarray}
```

::: warning This is not a mesh format.
It exports `point_data` / `cell_data` **for analytics** — pandas, polars, DuckDB — and does not round-trip geometry. It is deliberately *not* registered in meshio++'s format registry, so `meshioplusplus convert mesh.vtu out.parquet` does not work and will not be made to.
:::

For pandas and polars specifically there is also a [direct path](#pandas-and-polars) — `to_pandas` / `to_polars` — with no pyarrow or Parquet detour; Arrow remains the route when you need an on-disk file, a self-describing schema, or shared buffers.

**Point table** — `x`/`y`/`z` (as many as the mesh has) plus one column per `point_data` array.

**Cell table** — `block` (the block index), `cell_type`, `cell` (the **global block-major** cell index, the same numbering regions and [`partition_labels`](./partition) use), plus one column per `cell_data` array, concatenated block-major.

### Multi-component arrays

A `(N, 3)` array becomes a single Arrow **`fixed_size_list<double>[3]`** column, not three `v_0` / `v_1` / `v_2` columns. Flattening loses the shape, and `from_arrow` restores `(N, 3)` exactly:

```python
arrays = mio.from_arrow(mio.to_arrow(mesh))
arrays["v"].shape        # (N, 3)
```

### Schema metadata

The table is self-describing — the metadata survives a Parquet round-trip:

| Key | Value |
|---|---|
| `meshioplusplus:version` | the meshio++ version that wrote it |
| `meshioplusplus:location` | `point` or `cell` |
| `meshioplusplus:num_points` / `:num_cells` | the mesh's counts |
| `meshioplusplus:cell_types` | JSON `[[type, count], …]` in block order |
| `meshioplusplus:regions` | JSON `[{name, kind, dim, tag}, …]` |

```python
arrays, metadata = mio.read_parquet("data.parquet", return_metadata=True)
metadata["meshioplusplus:cell_types"]   # '[["triangle",1],["quad",1],["tetra",2]]'
```

Arrow's zero-copy is genuine for numeric buffers, and the implementation asserts it rather than assuming it — note that `ChunkedArray.combine_chunks()` copies even for a single chunk, so `from_arrow` takes the chunk directly.

### CLI

```bash
meshioplusplus data export mesh.vtu out.parquet --location cell
```

A sub-verb of the [`data` group](./cli#meshioplusplus-data). It exists in the **Python CLI only** — the native C++ CLI has no counterpart, since pyarrow is a Python library. Same caveat as above: it exports data, not a mesh.

## pandas and polars

```python
df = mio.to_pandas(mesh, location="cell")     # a pandas.DataFrame
pf = mio.to_polars(mesh, location="point")    # a polars.DataFrame
```

The same tabular export as [Arrow](#arrow-and-parquet) — same columns, same point/cell tables, same "not a mesh format" caveat — but straight to a frame, with **no pyarrow dependency and no Parquet detour**. Both are thin wrappers over the same `_to_table_payload` seam `to_arrow` uses.

### Multi-component arrays, revisited

The Arrow rule above — a `(N, 3)` array stays one `fixed_size_list` column, never `v_0`/`v_1`/`v_2` — holds for polars too: `to_polars` maps it to a **`pl.Array(Float64, 3)`** column, so the shape survives structurally.

pandas is the one target that genuinely cannot comply: a DataFrame column is one-dimensional, full stop. So `to_pandas` — and only `to_pandas` — flattens a `(N, k)` array into suffixed columns `v_0` … `v_{k-1}`, which is also the layout a training script wants. The grouping is recorded so nothing is lost:

```python
df = mio.to_pandas(mesh, location="point")
list(df.columns)                              # [x, y, z, T, v_0, v_1, v_2]
df.attrs["meshioplusplus:components"]         # '{"v":["v_0","v_1","v_2"]}'
```

If an expanded name collides with an existing column (a scalar `v_1` next to a vector `v`), the first-produced column wins and the later one is dropped with a warning — pandas cannot hold two columns of one name.

### Metadata

`df.attrs` carries all the [schema metadata](#schema-metadata) keys plus `meshioplusplus:components` above. Note pandas drops `attrs` on many operations, so read it before transforming the frame. polars has **no** attrs/metadata slot at all; the shape survives structurally, and `to_arrow` remains the self-describing-table path.

### Copy semantics

`to_pandas` takes `zero_copy_only` like every other `to_*`, and sharing is measured, never assumed. A multi-component array **always** copies here (the suffix expansion slices a row-major array), so under `zero_copy_only=True` it raises, pointing at `to_arrow` as the shared-buffer path for vector data. `to_polars` takes no `zero_copy_only` — polars always copies into its own Arrow-backed buffers, so the frame is independent of the mesh by construction.

### No `from_pandas` / `from_polars`

Deliberate, not an omission: the "from" direction of a table target returns plain arrays, not a `Mesh` (a table never carried the geometry), and a frame is one `{c: df[c].to_numpy() for c in df.columns}` away from the dict `from_arrow` already returns.

## Phase 2: Open3D and DOLFINx

`has_open3d()` and `has_dolfinx()` return `False`, and `to_open3d` / `to_dolfinx` raise `NotImplementedError` naming the phase. Both slot into the same pure-payload seam; the constraints are recorded here so the design does not have to be rediscovered.

**Open3D** holds a `TriangleMesh` (surface) or a `TetraMesh` (tets), so `_to_triangles_payload` covers the first and a `_to_tetra_payload` would cover the second. Two things make it unlike the targets above: `Vector3dVector` typically **copies**, so the zero-copy contract has to be stated per structure rather than globally, and the wheel is ~400 MB, which is a lot to put behind an extra people might install by reflex.

**DOLFINx** is the harder one. It accepts **single-cell-type meshes only**, needs a `ufl`/`basix` domain object alongside the arrays, is built around an MPI communicator rather than a process-local object, and — the actual work — requires a **node-ordering permutation** from VTK/gmsh ordering to DOLFINx/basix ordering. `dolfinx.io.gmshio`'s `cell_perm_gmsh` is the reference implementation. A `_to_dolfinx_payload` would own that permutation, keeping it testable without DOLFINx installed, exactly as the VTK payload is testable without VTK. DOLFINx is conda/apt-only, so it can never be a pip extra and its CI leg cannot be a plain `pip install`.

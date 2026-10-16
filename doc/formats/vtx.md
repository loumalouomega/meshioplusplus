# DOLFINx VTX (ADIOS2 `.bp`)

The output of DOLFINx's `VTXWriter` (`dolfinx.io.VTXWriter`): an ADIOS2 BP4 or BP5 *directory* (`u.bp/` holding `md.idx`, `md.0`, `data.0`…) with a `vtk.xml` schema attribute and, per time step, the mesh and the functions written. It is what ParaView opens with its ADIOS2 VTX reader. New in v16.13.0; written too since the writer joined.

| | |
|---|---|
| **Format name** | `vtx` |
| **Extensions** | `.bp` (a directory; one without the extension is found by its `md.idx` and `data.0`) |
| **Read / Write** | ✓ / ✓ ([conformance](../conformance.md#vtx)) |
| **Time steps** | ✓ — `time_step=`, `read_metadata()["time_values"]`, [`read_sequence`](/sequences) |
| **Extra dependencies** | ADIOS2: a core built with `MESHIOPLUSPLUS_WITH_ADIOS2=ON`, or the `adios2` Python package (`pip install meshioplusplus[adios2]`) |

## Reading

```python
import meshioplusplus

mesh = meshioplusplus.read("u.bp")                    # step 0
last = meshioplusplus.read("u.bp", time_step=-1)      # the last step
last.field_data["meshio:time"]                        # its time
meshioplusplus.vtx.time_values("u.bp")                # every step's time

for t, mesh in meshioplusplus.read_sequence("u.bp"):  # one mesh alive at a time
    ...

# a parallel run: weld each rank's ghost points onto their owners
serial = meshioplusplus.vtx.read("u.bp", ghosts="drop")
```

```bash
meshioplusplus info u.bp
meshioplusplus convert u.bp u.vtkhdf                  # every step
meshioplusplus convert 'run_*.bp' run.pvd             # a glob keeps .bp directories
```

ADIOS2 is a heavy dependency for one consumer, so it is never a default: the PyPI wheels have no ADIOS2 in their core, and read `.bp` through the `adios2` package when it is installed. A core built with `MESHIOPLUSPLUS_WITH_ADIOS2=ON` (see [installation](/installation)) reads it natively, which is also what the C, Fortran, Julia and R bindings and the native CLI need; the WASM build has no ADIOS2. The two engines give the same meshes. A process never holds both: a core built with ADIOS2 does not fall back to the Python reader, whose `adios2` package brings its own copy of the ADIOS2 libraries.

## Writing

```python
meshioplusplus.write("u.bp", mesh)                                   # one step
meshioplusplus.write_sequence("u.bp", steps, file_format="vtx")     # a step per mesh
```

```bash
meshioplusplus convert model.vtu u.bp
```

The writer produces `VTXWriter`'s layout for one rank (BP5): the `vtk.xml` schema, then per step `step` (the time, `field_data["meshio:time"]`), `NumberOfNodes` and `NumberOfEntities`, the rank-local arrays `geometry`, `connectivity` and one variable per point or cell array, and the single value `types`. The mesh is written at every step (`VTXMeshPolicy.update`). Cells are written as the `VTK_LAGRANGE_*` types DOLFINx uses wherever the node order allows (every degree-1 cell, `line3`, `triangle6`, `quad9`, `tetra10`), as plain VTK types otherwise. Like DOLFINx's files, a file holds **one cell type**: a mesh with several is a `WriteError` (write each type to its own file), as is a mesh with no cells. A series needs the first step's points, cells and array names at every step. Arrays keep their dtype; those of more than two dimensions, or whose rows are not the points' or cells', are dropped with a warning, and `geometry`, `connectivity`, `types`, `step`, `TIME`, `NumberOfNodes` and `NumberOfEntities` are reserved names. The ghost arrays are written only when the mesh carries them. A core built with ADIOS2 writes through its own library, otherwise the `adios2` package does; both write the same variables and values (the byte layout of a BP5 directory is ADIOS2's).

## The layout

`vtk.xml` is the VTK XML model with variable names in place of data: `Points` names the `geometry` variable (points × 3), `Cells` the `connectivity` (cells × (1 + nodes), each row prefixed with its node count, VTK legacy style) and `types` variables, `PointData`/`CellData` one variable per array, and a `TIME` array the `step` variable. Every step holds `step`; the mesh variables are written at every step (`VTXMeshPolicy.update`) or only at the first (`VTXMeshPolicy.reuse`). Every MPI rank writes one block of each variable, with its own point numbering and its ghost points.

## Mapping

- **Steps** are the file's ADIOS2 steps; a step's time is its `step` value (the step index where there is none), attached as `field_data["meshio:time"]`. A step without mesh variables takes the mesh of the last step that has one.
- **Ranks** are concatenated in block order, each rank's connectivity offset by the points before it, as ParaView merges them. `ghosts="drop"` (`ReadOptions::mGhosts = GhostPolicy::Drop`) instead welds every point onto the owned point with the same `vtkOriginalPointIds` and removes `vtkGhostType` and `vtkOriginalPointIds`; the result is the serial run's mesh (checked on a two-rank run of the fixture below).
- **Cells**: DOLFINx writes `VTK_LAGRANGE_*` types for every degree. A Lagrange cell whose node order is a linear or quadratic VTK cell's is read as that cell: every degree-1 cell (`line`, `triangle`, `quad`, `tetra`, `hexahedron`, `wedge`, `pyramid`), and the 3-node curve (`line3`), 6-node triangle (`triangle6`), 9-node quadrilateral (`quad9`) and 10-node tetrahedron (`tetra10`). Other degrees stay `VTK_LAGRANGE_*`, as the VTU reader reads them.
- **Data**: `PointData` arrays are point data and `CellData` arrays (DG0 functions) cell data, named as in the schema; a vector is padded to three components and a tensor to nine, as DOLFINx writes them; a complex function is two arrays, `<name>_real` and `<name>_imag`. `vtkOriginalPointIds` and `vtkGhostType` are kept as point data. `arrays=` and `points_only=` narrow the data.

## Quirks of DOLFINx 0.11 output

- The schema's `NumberOfCells` names a variable DOLFINx writes as `NumberOfEntities` when the mesh comes from a function space; the counts are taken from the `connectivity` blocks, never from those variables.
- In the first step of a `reuse` file DOLFINx writes `vtkOriginalPointIds` and `vtkGhostType` a second time, as one-value blocks whose payload is not what their count says, and ADIOS2 then misreads the next rank's block too. A block that matches no rank's size is never read, and a ghost array whose blocks do not match the ranks one to one is taken from the nearest step where they do (DOLFINx rewrites both, unchanged, every step). ParaView 6.1 appends the stray blocks to the arrays (26 ids for 25 points).
- ParaView 6.1 aborts on every step after the first of a `reuse` file; meshio++ reads them all.

## Not read

- `adios4dolfinx` checkpoints and Fides output are other ADIOS2 layouts, with no `vtk.xml`; such a `.bp` is a `ReadError` that says so.
- Ranks: the writer writes one rank; DOLFINx reads its own meshes from XDMF, which meshio++ also writes.

## Validation

The fixtures are real DOLFINx 0.11.0 runs (ADIOS2 2.12.1, conda-forge, September 2026): the heat equation (P1 triangles, three steps, `reuse`, BP5), on one rank and on two; a P2 vector and a DG0 field on tetrahedra (`update`, BP4); a mesh-only quadrilateral file; and a complex field in DOLFINx's layout. Every step both engines read matches ParaView 6.1's `ADIOS2VTXReader` (points, connectivity, cell types and arrays, bit for bit) where ParaView can read it, and the raw ADIOS2 arrays elsewhere; see `tests/python/meshes/vtx/README.md`.

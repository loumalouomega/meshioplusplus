# Polyhedral refinement (subdivide)

`subdivide(mesh, record_parent_ids=False)` splits every eligible 3D cell into
one polyhedral child per face, connected to a new interior point. It is a mesh
**operation**, not a file format, uses only standard C++, and runs under every
mesh backend.

`refine` and `decimate` both raise by name on a polyhedron block, pointing at
[`convert_cells(mode="simplexify")`](/convert_cells) — both are built on fixed
same-type subdivision templates, and an arbitrary polyhedron has none.
`subdivide` is the answer for polyhedral refinement specifically.

```python
import meshioplusplus as mp

mesh = mp.read("bracket.msh")
out = mp.subdivide(mesh, record_parent_ids=True)
mp.write("bracket_subdivided.vtu", out)
```

## No per-type template table

Unlike `refine`, which is driven by a fixed table per cell type
(`detail/refine_templates.hpp`), `subdivide` needs none. It goes through
`detail::cell_rings` — the same uniform face-ring abstraction `gradient` and
`compute_quality` already use — which treats a **tabulated** type (`tetra`,
`hexahedron`, `wedge`, `pyramid`, and their quadratic variants, reduced to
corners) and an **existing polyhedron block** identically. So the same code
handles every 3D cell type the mesh already supports.

## The construction

For each eligible 3D cell:

1. `detail::cell_rings` gives the cell's faces as global node-id rings,
   uniformly whether the cell is tabulated or already a polyhedron.
2. `detail::orient_rings` repairs winding so every face points outward; a cell
   whose faces are not a closed orientable surface **raises** (`ValueError` in
   Python), naming the cell — the same guard `convert_cells(mode="simplexify")`
   uses for its own polyhedron branch.
3. One new point is added per cell: the plain arithmetic mean of the cell's
   own corner node coordinates. This is deliberately **not** the volume
   centroid (`poly_measure`'s `mCentroid`, a different, volume-weighted
   point) — the corner average is what makes the children's total volume
   *literally the same sum* as the parent's own volume computation.
4. For each face, one polyhedron child is emitted whose boundary is that face
   **unchanged** (original winding, so a neighbouring cell across it still
   sees the identical face) plus one new triangle per face edge, connecting
   that edge to the cell's new interior point.

This makes the result **automatically conforming**: no closure, no 2:1
balance, no `refine:hanging`/`refine:entity` analogue is needed, because a
shared face between two input cells is never touched — only the interior of
each cell is subdivided.

## Output structure

**One polyhedron output block per input 3D block**, with genuinely mixed cell
shapes inside — `AddPolyhedronBlock` stores cells as ragged CSR with no
constraint that they share a node or face count, so there is no need to group
children by distinct node count into `polyhedronN` blocks the way some
*readers* (CGNS's `NFACE_n`, OpenFOAM, EnSight) do for their own
format-compatibility reasons. A single hexahedron already demonstrates this: a
wedge cell's 5 faces (2 triangular, 3 quadrilateral) subdivide into 2
tetrahedron-shaped children and 3 square-pyramid-shaped children, all in one
coherent output block.

Non-3D blocks (2D/1D boundary markers) and 3D blocks with no `cell_faces` row
(the full-Lagrange family — `hexahedron64`, `tetra20`, the `VTK_LAGRANGE_*`
types) pass through unchanged.

## Volume conservation

Since every child is a fan from the parent's own corner average, the
children's total volume equals the parent's — checkable near-exactly, not just
by the divergence theorem's abstract equality:

```python
before = mp.compute_stats(mesh)["signed_volume"]
out = mp.subdivide(mesh)
after = mp.compute_stats(out)["signed_volume"]
assert abs(after - before) < 1e-9 * abs(before)
```

## Regions

Point and Cell regions (and so `point_sets`/`cell_sets`) survive: a parent's
children occupy a contiguous run in the output block, so the carry uses
`CellMapKind::FirstChild` — the same shape `convert_cells` already uses for
its own one-to-many splits. Named **Side** regions do not survive: `FirstChild`
drops them unconditionally, since a child's facets have no correspondence with
the parent's — the same limitation `convert_cells(mode="simplexify")` already
has and documents, not a new gap this operation introduces.

There is no point map to request (unlike `convert_cells`): `subdivide` never
prunes or renumbers an original point, and new apex points are by
construction unreferenced by any existing Point region.

## Naming

`record_parent_ids=True` attaches an Int64 `subdivide:parent_cell` cell_data
array recording, per output cell, the index of the input cell it came from
**within its own block** (blocks correspond 1:1 with the input).

## No numpy fallback

This operation is **C++-core only**, with no pure-Python reference
implementation at all — unlike most operations in this repo, which fall back
to numpy when the compiled core is unavailable. The winding repair
(`orient_rings`) is a discrete branch on the sign of an enclosed volume, and a
second, independently written implementation of that branch could land on the
opposite side for a near-degenerate cell and then diverge macroscopically
rather than by round-off — the same reasoning `_convert_cells.py`'s polyhedron
branch and `_smooth.py`'s inversion guard already document. A pure-Python
build (`meshioplusplus._core` unavailable) raises `NotImplementedError` naming
the reason, for any input, tabulated or not.

## CLI

```sh
meshioplusplus subdivide IN OUT [--record-parent-ids]
```

Both CLIs produce byte-identical output (there being only one implementation).
See the [CLI reference](/cli#meshioplusplus).

## Other languages

```c
mio_subdivide_result* r = mio_subdivide(mesh, /*record_parent_ids=*/1);
const mio_mesh* out = mio_subdivide_result_mesh(r);
/* ... */
mio_mesh* owned = mio_subdivide_result_take_mesh(r);
mio_subdivide_result_free(r);
```

```fortran
type(mio_mesh) :: s
s = m%subdivide(record_parent_ids=.true., stat=st)
```

```julia
s = subdivide(mesh; record_parent_ids=true)
s.mesh, s.cell_maps
```

```r
s <- mio_subdivide(mesh, record_parent_ids = TRUE)
s$mesh; s$cell_maps
```

```js
const out = await m.subdivide(mesh, true);
```

On the flat ABIs (C, Fortran, Julia, R, WASM), `SubdivideResult` carries no
point map at all — a deliberate divergence from `convert_cells`'s shape, since
subdivide never prunes or renumbers a point. Fortran/WASM additionally skip
exposing the per-block cell maps (a documented flat-ABI gap, the same one
`convert_cells`'s own Fortran/WASM bindings already have).

This operation is also reachable as a `Subdivide` step in the
[settings pipeline](/pipeline) and in the browser viewer's `convertSurfaceOps`
chain (`{op: 'subdivide', recordParentIds: true}`).

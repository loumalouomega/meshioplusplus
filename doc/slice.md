# Slicing / cross-sections

`meshioplusplus.slice(mesh, origin, normal)` computes the **planar cross-section** of a mesh — the actual intersection of the mesh with a plane, one topological dimension below the cut cells. It is a mesh **operation**, not a file format, uses only standard C++/numpy, and runs under every mesh backend.

- a **3D volume** mesh cut by the plane yields a **surface** mesh (`triangle`/`quad` section faces);
- a **2D surface** mesh cut by the plane yields a **line** mesh (segments).

This is different from [`crop`](/crop) in plane mode, which keeps *whole cells* on one side of the plane: `slice` computes the intersection itself and lowers the dimension.

Its data-driven sibling is [`isosurface`](/isosurface), which cuts where a scalar field equals an isovalue rather than where the distance to a plane is zero. The two share one marching-tetrahedra cutter, so everything below about watertightness, winding, degeneracy and determinism holds identically there.

```python
import meshioplusplus

vol = meshioplusplus.read("part.vtu")            # a tetra / hex / wedge mesh

# the cross-section at z = 0.5
section = meshioplusplus.slice(vol, origin=(0, 0, 0.5), normal=(0, 0, 1))

# record which input cell each section face was cut from
section = meshioplusplus.slice(
    vol, origin=(0, 0, 0.5), normal=(0, 0, 1), record_parent_ids=True
)

meshioplusplus.write("section.vtu", section)
```

The plane is `origin` (a point on it) plus `normal` (its normal — non-zero, need not be unit length; a zero normal raises). The result is a brand-new mesh whose points are all cut points, so `point_sets`/`cell_sets`, `mesh.info` and `gmsh_periodic` are **not** carried; each section cell inherits its parent cell's `cell_data`, and interpolated `point_data` is promoted to Float64.

## How it works — marching tetrahedra

The input is first simplexified with [`convert_cells`](/convert_cells) (every 3D cell becomes a tetrahedron, every 2D cell a triangle), so each cell's cross-section is a well-defined convex primitive — no angular point-sorting and no non-convex-hex ambiguity. Per simplex the plane crosses 0, 3 or 4 edges of a tetra (a `triangle` or `quad`) or 0 or 2 edges of a triangle (a `line`), from a fixed case table keyed by the vertices' above/below-plane pattern.

The honest consequence: a hexahedron or wedge section is the **union of its simplices' sections** — correctly located and shaped, just triangulated rather than returned as a single polygon.

The signed distance of node *i* to the plane is `d_i = dot(x_i − origin, normal)`. An edge `(i, j)` is crossed when `d_i` and `d_j` are on opposite sides; the crossing point is `x_i + t·(x_j − x_i)` with `t = d_i / (d_i − d_j)`, and every `point_data` array is interpolated at the same `t`.

## Degeneracy rule

A node exactly on the plane (`d_i == 0`) is classified on the **positive** side (`d_i >= 0`). This one rule makes the sign mask *total*, which is what keeps the section well-defined at the awkward cases:

- a plane **grazing a shared internal face** is emitted by exactly one of the two incident cells (the one whose off-plane apex is below), so the face appears **once**, never doubled;
- any section primitive whose crossings collapse to a point or a line (area or length below a bbox-relative tolerance) is dropped, so there are no zero-area/zero-length section cells.

## Watertight sections

A crossing point that lies on an edge shared by two simplices is emitted as a **single** output node: cut points are keyed by the sorted pair of simplexified-node ids and deduped, so the section is watertight (a crossing on a shared edge is one node, referenced by both faces). The one documented exception is a **grazing** configuration (the plane passing exactly through mesh vertices/edges), where coincident-but-distinct nodes can remain — the section is still geometrically correct, just not vertex-minimal there.

## Winding

Every section face is wound so that its Newell normal points toward the `+normal` side of the plane — a globally consistent orientation, pinned by a Newell-normal invariant test (`tests/cpp/test_slice.cpp`).

## Provenance

With `record_parent_ids=True`, an Int64 `slice:parent_cell` `cell_data` array records, per section cell, the **global (block-major)** index of the input cell it was cut from — mirroring `surface:parent_cell` from [surface extraction](/extract_surface). The parent's `cell_data` row is replicated onto each section cell regardless of this flag.

## Determinism

Output is byte-identical across the three mesh backends, across thread counts, and across the C++-core / numpy-fallback boundary (pinned by `tests/python/test_slice.py::test_cpp_matches_python`): the signed-distance and coordinate/`point_data` passes are parallel, but the cut-and-dedup sweep that hands out node ids is **serial** over a fixed `(block, cell, ring-edge)` traversal, and the crossing `t` is computed from the sorted endpoint pair so both sides of a shared edge agree bit-for-bit.

## CLI

```sh
meshioplusplus slice IN OUT --origin 0,0,0.5 --normal 0,0,1
meshioplusplus slice part.msh section.vtu --normal=0,0,-1 --record-parent-ids
```

Negative components use the `--normal=…` / `--origin=…` form (a bare `--normal -1,0,0` would be read as a new option). Available in both the Python CLI and the native `meshioplusplus` binary — see the [CLI reference](/cli#meshioplusplus-slice).

## Other languages

The operation is exposed on every binding surface:

```c
/* C API */
double origin[3] = {0, 0, 0.5}, normal[3] = {0, 0, 1};
mio_mesh* section = mio_slice(mesh, origin, normal, /*record_parent_ids=*/1);
```

```fortran
! Fortran: type-bound procedure (like m%transform / m%crop_plane)
section = mesh%slice([0.d0, 0.d0, 0.5d0], [0.d0, 0.d0, 1.d0], &
                     record_parent_ids=.true.)
```

```js
// WASM
const section = m.slice(mesh, [0, 0, 0.5], [0, 0, 1], true);
```

As everywhere on the flat bindings, `point_sets`/`cell_sets` never cross — and here they are not carried at all, since the section is entirely new topology. The browser [viewer](/viewer)'s planar-section feature runs this operation.

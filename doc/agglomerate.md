# Polyhedral coarsening (agglomerate)

`meshioplusplus.agglomerate(mesh, target_group_size=8, merge_coplanar_faces=False, coplanar_angle=1.0, min_sphericity=0.0)` merges groups of cells into single larger polyhedral cells — the many-to-one counterpart to [`subdivide`](/subdivide). It is a mesh **operation**, not a file format, uses only standard C++, and runs under every mesh backend.

`decimate` raises by name on a polyhedron, pointing at `convert_cells(mode="simplexify")` — its fixed-template QEM edge collapse has no analogue for merging arbitrary polyhedral cells. `agglomerate` is a genuinely different algorithm: greedy seed-and-grow over the mesh's shared-face dual, then one polyhedron per group whose faces are exactly that group's *external* boundary.

```python
import meshioplusplus as mp

mesh = mp.read("bracket.msh")
coarse = mp.agglomerate(mesh, target_group_size=8)
mp.write("bracket_coarse.vtu", coarse)
```

## The construction

![Greedy seed-and-grow over shared faces on a grid of unequal cells, showing which neighbour each group absorbs and which faces are dropped](/diagrams/agglomerate_grow.svg)

1. `detail::build_global_faces(mesh)` re-expresses the mesh's volume cells as a globally deduplicated face list with owner/neighbour pairing — the same machinery the OpenFOAM writer and the CGNS `NGON_n`/`NFACE_n` writer share. A mesh with any **non-manifold** face (used by three or more cells) is **refused**: the owner/neighbour classification below is only well-defined on a 2-manifold face, and guessing would silently misclassify a boundary.
2. **Greedy seed-and-grow**, serial and deterministic: cells are seeded in ascending order; a group absorbs its unclaimed face-neighbour with the largest *accumulated* shared-face area (summed over every face the group's current members share with that neighbour) until it reaches `target_group_size`, or no unclaimed neighbour remains — a short group at a mesh boundary or pocket is expected, not an error.
3. Each group emits **one** polyhedron cell: walk every member's faces. A face whose other side is also in the group is internal and dropped (this happens from both sides, so it is never emitted twice); every other face is the group's own external boundary and is kept, wound exactly as the member's own local orientation already records it — no new orientation logic is needed.

This is deliberately **face-adjacency**, never node-adjacency (the kind `partition`'s ghost layers and `gradient`'s stencil use): merging on shared-node adjacency could fuse two cells touching only at a single pinch-point vertex, producing a non-manifold union.

## Compact groups: the sphericity gate

Greedy growth by shared area favours compact groups, but on a stretched or graded mesh a group can still grow a long arm along the dual. `min_sphericity` (default `0`, off) refuses to absorb a candidate when the union would be less round than that. The measure is the sphericity `π^(1/3) (6V)^(2/3) / A`: 1 for a ball, about 0.81 for a cube, lower for anything elongated or flat. The group's volume `V` and external area `A` are kept incrementally and exactly: absorbing a cell adds its volume, adds its area and subtracts twice the area it shares with the group. A refused candidate is skipped for that group only and can still seed or join another group, so the grouping stays deterministic. Groups may end up smaller than `target_group_size`. The report's `num_rejected` counts the refusals.

```python
coarse, report = mp.agglomerate(mesh, target_group_size=8, min_sphericity=0.7, return_report=True)
report["num_rejected"]
```

## Fusing coplanar faces

Without merging, a coarse cell keeps every fine face on its boundary: a flat wall that eight hexahedra tiled stays as eight quads. `merge_coplanar_faces=True` fuses such faces into one polygon. Faces are grouped by the pair of sides they separate (two groups, or a group and the mesh boundary), so both neighbours get the identical fused polygon (reversed) and the mesh stays conforming. Within a pair, edge-adjacent faces whose normals lie within `coplanar_angle` degrees (default 1) of the patch's first face fuse when the patch's outline is one simple loop, with no hole and no pinch point. Every polyhedron a fusion touches is checked to remain closed; if a fusion would break that, the patch keeps its original faces. Volume is conserved either way, since the fused polygon covers exactly the faces it replaces. The report's `num_faces_merged` counts the faces removed (fused faces minus the polygons that replace them, per group side).

Fusing drops the edge-interior nodes of the patch from the face lists. Points are still never renumbered, so the result is valid; [`clean`](/clean) with `remove_orphans=True` compacts the point set if a node is no longer referenced at all.

## Output structure

Non-volume blocks (2D/1D boundary markers, and any 3D block with no face table) pass through unchanged, in their original relative position. Every volume cell is consumed into **one** new `polyhedron` block, emitted at the position the *first* original volume block occupied — so a mesh whose volume cells already form one contiguous run keeps its overall block order. Cells inside that block have whatever face/node count their own group boundary produces; `AddPolyhedronBlock` stores genuinely ragged CSR with no same-shape constraint, so there is no grouping-by-size step to do — the same simplification `subdivide` already established.

## Points and data

Points are **not** compacted: a group can leave an interior node unreferenced — the same never-prune-or-renumber precedent `subdivide` set for its own orphan case. [`clean`](/clean) with `remove_orphans=True` is the documented follow-up for a caller who wants a minimal point set. `point_data` and `field_data` therefore need no remapping at all.

`cell_data` for a pass-through block is copied verbatim. For the merged block, each group's row is its **first member's** row (ascending order within the group) — the same keep-first convention `clean`'s point weld already uses, since there is no single principled value for an array like a material tag once several cells with (possibly) different values merge. An array whose block count does not match the input mesh is dropped with a warning rather than guessed at.

## Volume conservation

Since every surviving face is a face the input mesh already had (just regrouped), a merge conserves volume **exactly** — checkable to an identity, not a tolerance:

```python
before = mp.compute_stats(mesh)["signed_volume"]
out = mp.agglomerate(mesh, target_group_size=2)
after = mp.compute_stats(out)["signed_volume"]
assert after == before
```

## Regions

Point and Cell regions survive: the carry uses `CellMapKind::Global`, a single flat map (input global cell → output global cell). Named **Side** regions survive too (since v16.26.0): a facet of a merged cell moves to the output polyhedron's face containing it, matched by its nodes. With `merge_coplanar_faces=True` that face may be the merged polygon, and several facets then name the same face. See [how operations treat regions](/regions#how-operations-treat-regions).

There is no point map to request: `agglomerate` never prunes or renumbers an original point.

## No numpy fallback

This operation is **C++-core only**, with no pure-Python reference implementation — the same reasoning `subdivide` already documents: the emit step depends transitively on a winding repair (a discrete branch on the sign of an enclosed volume) that a second, independently written implementation could disagree with near-degenerate cells. A pure-Python build (`meshioplusplus._core` unavailable) raises `NotImplementedError` naming the reason, for any input.

## CLI

```sh
meshioplusplus agglomerate IN OUT [--target-group-size N] [--merge-coplanar-faces]
                             [--coplanar-angle DEG] [--min-sphericity S] [--json]
```

`--json` prints `{cells_in, cells_out, num_faces_merged, num_rejected}`.

Both CLIs produce byte-identical output (there being only one implementation). See the [CLI reference](/cli#meshioplusplus).

## Other languages

```c
mio_agglomerate_result* r = mio_agglomerate(mesh, /*target_group_size=*/8);
const mio_mesh* out = mio_agglomerate_result_mesh(r);
/* ... */
mio_mesh* owned = mio_agglomerate_result_take_mesh(r);
mio_agglomerate_result_free(r);
```

```c
mio_agglomerate_opts opts;
mio_agglomerate_opts_init(&opts);           /* groups of 8, no merging, no gate */
opts.merge_coplanar_faces = 1;
opts.min_sphericity = 0.7;
int64_t merged = 0, rejected = 0;
mio_agglomerate_result* r = mio_agglomerate_ex(mesh, &opts, &merged, &rejected);
```

```fortran
type(mio_mesh) :: c
c = m%agglomerate(target_group_size=8_int64, stat=st, merge_coplanar_faces=.true., &
                  min_sphericity=0.7_real64, num_faces_merged=merged, num_rejected=rejected)
```

```julia
a = agglomerate(mesh; target_group_size=8, merge_coplanar_faces=true, min_sphericity=0.7)
a.mesh, a.cell_map, a.num_faces_merged, a.num_rejected
```

```r
a <- mio_agglomerate(mesh, target_group_size = 8, merge_coplanar_faces = TRUE, min_sphericity = 0.7)
a$mesh; a$cell_map; a$num_faces_merged; a$num_rejected
```

```js
const out = await m.agglomerate(mesh, 8);
const r = await m.agglomerate(mesh, 8, true, { mergeCoplanarFaces: true, minSphericity: 0.7 });
```

On the flat ABIs, `AgglomerateResult` carries a single **flat** cell map (unlike `SubdivideResult`'s per-block one) — an output cell's index is a function of which group it joined, not which input block it came from. The C API's `mio_agglomerate_result_cell_map` accordingly takes **no `block` parameter**, and there is correspondingly no `mio_agglomerate_result_num_cell_maps` (there is exactly one array). WASM's `agglomerate(mesh, targetGroupSize, returnMaps, options)` returns `{mesh, cellMap, numFacesMerged, numRejected}` (the same flat map) when `returnMaps` is set.

This operation is also reachable as an `Agglomerate` step in the [settings pipeline](/pipeline) (keys `TargetGroupSize`, `MergeCoplanarFaces`, `CoplanarAngle`, `MinSphericity`) and in the browser viewer's `convertSurfaceOps` chain (`{op: 'agglomerate', targetGroupSize: 8, mergeCoplanarFaces: true}`).

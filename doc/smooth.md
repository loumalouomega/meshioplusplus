# Smoothing (Laplacian / Taubin / ODT)

`meshioplusplus.smooth(mesh, method="taubin", iterations=10)` improves **element shape** by relaxing each free node toward the centroid of its edge neighbours (`method="odt"` instead moves each free interior tet vertex toward the volume-weighted average of its incident tets' circumcenters — see [ODT smoothing](#odt-smoothing) below). It is a mesh **operation** (like [refine](/refine) and [clean](/clean)), not a file format, and uses only standard C++/numpy, so it runs under every mesh backend.

It is a **pure coordinate move**. Connectivity, cell counts, `cell_data`, `field_data` and `point_data` *values* all come through unchanged; only the point coordinates differ, and the points array keeps its input dtype. That makes `smooth` the geometric counterpart to [mesh quality](/mesh_quality), which only measures the thing this operation fixes.

![Element shape repaired: a jittered hexahedron block before and after smoothing, coloured by scaled Jacobian](/images/smooth_quality.png)

```python
import meshioplusplus

mesh = meshioplusplus.read("rough.msh")

# the default: Taubin, 10 iterations, boundary and features pinned
better = meshioplusplus.smooth(mesh)

# Laplacian is stronger per pass but shrinks the mesh
strong = meshioplusplus.smooth(mesh, method="laplacian", iterations=5)

# ask what actually happened
out, report = meshioplusplus.smooth(mesh, iterations=20, return_report=True)
print(report)   # {'num_nodes_moved': 27, 'max_displacement': 1.012, ...}

meshioplusplus.write("better.vtu", better)
```

## Two methods

Both are driven by the same edge-neighbour centroid displacement `L(i) = mean(x[j] for j adjacent to i) - x[i]`.

| method | per iteration | shrinks? |
|---|---|---|
| `laplacian` | one `x <- x + lambda*L(x)` pass (`lambda` defaults to 0.5) | **yes**, unavoidably |
| `taubin` | a `+lambda` pass (0.33) then a `-mu` pass (−0.34) | no |

**Laplacian** is a low-pass filter with no pass band: every spatial frequency is attenuated, including the zeroth. Run to convergence, a closed surface collapses toward a point. It is the right choice when a handful of passes is enough and volume loss does not matter.

**Taubin** follows each smoothing pass with a deliberate *un*-smoothing pass whose larger magnitude overshoots, so the pair has a **pass band** below `k = 1/lambda + 1/mu` where features are preserved, and a stop band above it where noise is attenuated. With the defaults that cutoff is ≈ 0.089, and hundreds of iterations leave the enclosed volume essentially unchanged. This is why it is the default.

The difference is not subtle — the same noisy grid, 40 iterations of each, with the boundary deliberately left free:

![Laplacian collapses the mesh over 40 iterations; Taubin denoises it while holding its size](/images/smooth_shrinkage.png)

`lambda` must lie in `(0, 1)`, and Taubin additionally requires `mu < -lambda < 0` — otherwise the pass pair amplifies instead of filtering, and `smooth` raises rather than diverging quietly. Passing a **negative `lambda`** selects that method's own default (0.5 or 0.33); this is why the parameter's default is `-1.0` rather than a number.

::: tip Python naming
The Python parameter is `lambda_`, with a trailing underscore, because `lambda` is a Python keyword. Every other surface spells it `lambda`.
:::

## ODT smoothing

`method="odt"` is optimal-Delaunay-triangulation smoothing: rather than the edge-neighbour centroid `laplacian`/`taubin` use, each free **interior** tet vertex moves toward the volume-weighted average of its incident tets' circumcenters. It is closed-form, not iterative-within-a-pass — one pass computes the target directly from the current mesh, the same Jacobi update `laplacian`/`taubin` already use.

**Tet-only.** `smooth` otherwise works on any mesh; `method="odt"` raises by name on any non-`tetra` block, rather than silently pinning or ignoring it. This is what closes the "ODT" half of `doc/roadmap.md`'s "Volumetric CVD/ODT" bullet — honestly, as ODT *smoothing* on an existing tet mesh's fixed connectivity, not ODT *remeshing* (which would need to change the connectivity itself; see [`remesh_volume`](/remesh_volume) for the operation that actually generates a fresh tet mesh).

```python
mesh = meshioplusplus.read("tangled.vtu")  # a tetra mesh

better, report = meshioplusplus.smooth(mesh, method="odt", iterations=10, return_report=True)
```

`lambda`'s negative-sentinel-means-"this method's own default" convention (above) extends to ODT: its default is `0.9`, a near-full step toward the circumcenter average, damped rather than exactly `1.0` for the same reason the other two methods are not exactly `1.0`. `mu` is silently ignored by ODT, matching `laplacian`'s own silent ignoring of it — deliberately *not* a named-error option the way `RemeshOptions.max_anisotropy` is under the wrong metric, since making only the newest method strict while its sibling stays lenient would be a worse inconsistency than either rule applied uniformly.

**Boundary pinning still applies**: with `fix_boundary=True` (the default), a tet mesh's boundary vertices are pinned exactly as under `laplacian`/`taubin`, so ODT only ever relaxes the interior. The inversion guard (below) applies unchanged too — an unguarded ODT step can invert a tet near a concave boundary, which is why the numpy fallback (below) refuses `method="odt"` unconditionally rather than only under `guard_inversion=True`.

**This is C++-core only, with no numpy fallback at all** — unlike `laplacian`/`taubin`, whose continuous part the fallback does implement (just not the inversion guard). ODT's target is a genuinely different closed-form computation (a circumcenter, not a centroid), and pairing it with the same "raise unless the guard is off" fallback posture would still leave an unguarded ODT step meaningless on a concave mesh; `meshioplusplus.smooth(mesh, method="odt")` raises `NotImplementedError` naming the reason when `_core` is unavailable, for any input.

## The neighbour graph is edges, not cliques

Two nodes are neighbours **iff they share a cell edge** — the same `cell_refine_edges` table [refine](/refine) uses, so the two cannot drift. Ragged polygon rows and polyhedron faces contribute their ring edges.

This matters. [Reordering](/reorder) builds the *clique* graph (every pair of nodes sharing a cell) because that is the sparsity pattern of the FEM matrix whose bandwidth it exists to reduce. Smoothing on that graph would pull each hexahedron corner toward its three face diagonals and its body diagonal as hard as toward its three real edge neighbours, bevelling cubes toward spheres. On the edge graph, a structured hexahedron block is a fixed point.

## What is pinned

A node never moves if it is:

- **on the boundary**, when `fix_boundary=True` (the default). Boundary is the once-used-facet test: a face used by exactly one 3D cell, or — in a pure surface mesh — an edge used by exactly one 2D cell. Without this the domain's own shape changes: a closed surface contracts and an open domain's walls creep inward.
- **a feature node**, when `preserve_features=True` (the default). These are boundary nodes where two incident boundary facet *normals* differ by more than `feature_angle` degrees (default 30; 0° is coplanar, 90° is the edge of a box — the `vtkFeatureEdges` convention). Without it, unfeatured boundary smoothing rounds every geometric corner off.
- **named in `frozen`** — an index array, or the name of a `mesh.point_sets` entry. Unioned with the other pins, never subtracted from them.
- **referenced only by blocks whose edge topology is unknown** — the higher-order family (`tetra10`, `quad8`, …), the VTK-Lagrange types, `custom`, `vertex`. A `tetra10`'s mid-edge nodes have no meaningful centroid: moving them would silently flatten a deliberately curved mesh, and moving only the corners would strand the mid-nodes off their edges. An unknown neighbourhood gives no defined target, so the node holds still — and says so once, via a warning.

## The inversion guard

Smoothing can in principle push a node through a neighbouring face. With `guard_inversion=True` (the default), a node's candidate position is committed only if no incident cell's signed measure would change from valid to inverted; otherwise the node holds still **for that pass** and the event is counted in `num_skipped_inversion`.

The rule is **"do no harm", not "preserve the sign"**:

- a currently-**valid** cell may not be made invalid — the guard's whole purpose;
- a cell that arrives **already inverted** imposes no constraint at all. Blocking its sign change would mean the guard prevents smoothing from *repairing* a tangle, which is one of the main reasons to reach for this operation. (An early version did exactly that: 5 inverted cells in, 4 still inverted out with the guard on, versus 0 with it off.)
- an exactly **degenerate** cell likewise imposes nothing, so a sliver's nodes are not pinned forever.

Three measures are used, by cell family:

| cells | measure |
|---|---|
| `tetra`, `hexahedron`, `wedge`, `pyramid` | signed volume by the outward face fan |
| 2D cells in a 2D mesh (`PointDim == 2`) | signed shoelace area |
| 2D cells in a 3D mesh | **normal flip**: `dot(n_before, n_after) <= 0` |

The third has no counterpart in [`compute_quality`](/mesh_quality), which correctly declines to call a triangle floating in 3D "inverted" — a facet in space has no intrinsic orientation. Smoothing needs only the weaker, purely relative question *did this facet just fold back over itself?*, so it asks that one, and a triangle surface mesh is guarded rather than left free to fold.

Polyhedron cells **do** have a signed measure since v9.16.0 (the corner-average fan over the cell's own faces), so the guard covers them — and their boundary is detected too, so `fix_boundary` pins their outer nodes. Cells with no signed measure at all (`line`, `custom`) are **skipped, not pinned**: a missing safety check is not a missing target, and pinning there would freeze an entire beam-element curve for no geometric reason.

The guard costs a second CSR (node → incident cell) plus a corner table, both built only when it is enabled.

## What changes

- **Points move. Nothing else does.** Point and cell counts, block structure, connectivity, `cell_data`, `field_data` and `point_data` *values* are all carried through unchanged, and the points array keeps its dtype (iteration runs in `double` regardless, with a single cast on write-back, so a Float32 mesh does not accumulate one rounding per pass).
- **`point_sets` / `cell_sets` pass through untouched** — unlike [refine](/refine)/[crop](/crop)/[split](/split), whose shims must remap them, this operation never adds, removes or renumbers a node or a cell.
- `mesh.info` and `gmsh_periodic` are not carried through (as for every other operation).
- The report carries `num_nodes_moved` (nodes whose **net** displacement exceeded a scale-relative tolerance), `max_displacement` (measured against the input, not the previous pass) and `num_skipped_inversion` (a count of `(node, pass)` rejection *events*, so one node rejected in five passes counts five).

## Determinism

Every pass is **Jacobi**: all new positions are computed from the previous pass's positions and committed together, so no node ever observes a half-updated neighbour and the result cannot depend on traversal order. The update loop contains no hashing, no sorting and no floating-point reduction; neighbour sums run in ascending node id, the boundary set comes from a serial dedup pass over a parallel-filled buffer (the same phase split [`extract_surface`](/extract_surface) uses), and the summary counters are folded serially out of per-node slots.

Output is byte-identical across the MESHIO/NATIVE/KRATOS backends and across thread counts — verified on a 52k-node mesh.

The numpy fallback reproduces the continuous part of the algorithm to within round-off, but **does not implement the inversion guard**: that is a discrete branch on the sign of a cell measure, and near a degenerate cell two implementations could legitimately land on opposite sides and then diverge macroscopically. Rather than ship a second, subtly different guard, the fallback raises if one is requested — pass `guard_inversion=False` to use it.

## CLI

```bash
meshioplusplus smooth rough.msh better.vtu
meshioplusplus smooth rough.msh better.vtu --method laplacian --iterations 5
meshioplusplus smooth rough.msh better.vtu --feature-angle 45 --no-guard-inversion
meshioplusplus smooth rough.msh better.vtu --lambda 0.4 --mu=-0.45
meshioplusplus smooth tangled.vtu better.vtu --method odt --iterations 10
```

Negative values need the `--mu=-0.45` form, as for [transform](/transform) and [crop](/crop). See the [CLI reference](/cli).

## Other languages

- **C API** — `mio_smooth(mesh, method, iterations, lambda, mu, fix_boundary, preserve_features, feature_angle, guard_inversion, &nodes_moved, &max_displacement, &skipped_inversion)` returns a plain `mio_mesh*`; the three counter out-params are individually optional. `frozen` is not exposed on the flat ABI. `method` is passed straight through, so `"odt"` needs no C API surface of its own beyond the string. See the [C API reference](/c_api).
- **Fortran** — `mesh%smooth(method, iterations [, lambda, mu, fix_boundary, preserve_features, feature_angle, guard_inversion, nodes_moved, max_displacement, skipped_inversion])`. See the [Fortran reference](/fortran).
- **Julia** — `smooth(mesh; method="taubin", iterations=10, ...)`. See the [Julia reference](/julia).
- **R** — `mio_smooth(mesh, method="taubin", iterations=10L, ...)`. See the [R reference](/r).
- **WebAssembly / JavaScript** — `smooth(mesh, method, iterations, lambda, mu, fixBoundary, preserveFeatures, featureAngle, guardInversion)`, returning `{mesh, numNodesMoved, maxDisplacement, numSkippedInversion}`. See the [WebAssembly reference](/wasm).

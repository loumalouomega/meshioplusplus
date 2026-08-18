# Conservative interpolation (mass-preserving field transfer)

`meshioplusplus.conservative_interpolate(source, target)` mass-preservingly
transfers a **source** mesh's data arrays onto a **target** mesh: over the
region the two meshes share, `sum(target value * target measure)` equals
`sum(source value * source measure)`. This is the property
[interpolate](/interpolate)'s `"barycentric"` mode does not have — a
pointwise sampler has no notion of "how much of the source region a target
sample stands for". CFD/FEM solver coupling and successive-remeshing
workflows need exactly this guarantee; `interpolate` remains the right tool
for genuinely pointwise sampling (probing a field at coordinates,
visualisation resampling, ...).

It is a mesh **operation**, not a file format, deliberately a **separate
operation** from `interpolate` rather than a third method — a genuinely
different algorithm (measure-weighted geometric overlap, not point sampling)
with different guarantees, mirroring the [decimate](/decimate) /
[decimate_volume](/decimate_volume) split.

The result is a **new mesh that is a copy of the target**: geometry,
connectivity, the target's own data arrays and its `point_sets`/`cell_sets`
are preserved exactly, with the requested source arrays remapped onto it.
Unlike `interpolate`, an unset `arrays` transfers **every source `point_data`
AND `cell_data` array** — there is one algorithm regardless of location, so
`interpolate`'s "cell_data only when named" special case does not apply here.
Output arrays are always Float64 (a measure-weighted mean is not integral).

![A coarse cell_data field mass-preservingly remapped onto a finer, offset target mesh](/images/conservative_interpolate_transfer.png)

```python
import meshioplusplus

coarse = meshioplusplus.read("solution.vtu")   # carries cell_data "pressure"
fine = meshioplusplus.read("remeshed.vtu")

out = meshioplusplus.conservative_interpolate(coarse, fine)

# name the arrays; uncovered target cells get -1 instead of 0
out = meshioplusplus.conservative_interpolate(
    coarse, fine, arrays=["pressure"], default_value=-1.0
)

meshioplusplus.write("mapped.vtu", out)
```

## Algorithm

Both meshes are first reduced to triangles (2D) or tetrahedra (3D) via
[convert_cells](/convert_cells)`(mode="simplexify")` — the same call
`interpolate`'s `"barycentric"` mode already makes. This is what lets the
operation accept **ragged and polyhedron blocks for free**: Simplexify
already fans them into simplices via a shipped, tested path, so no general
polygon/polyhedron clipper is needed. `source` and `target` must share the
same maximum topological dimension; there is no cross-dimension remap.

For every pair of overlapping simplices (found via a bucket-grid spatial
hash over bounding boxes, `detail/spatial_hash.hpp`, the same broad-phase
`interpolate`'s barycentric mode already uses), the **exact** overlap
measure is computed:

- **2D**: a Sutherland-Hodgman convex polygon clip of the target triangle
  against the source triangle's three edges, measured with the resulting
  polygon's area (the same xy-plane assumption `interpolate`'s barycentric
  triangle path documents);
- **3D**: since both operands are always tetrahedra, a bounded (not
  general-polyhedron) convex-polytope clip — the source tetrahedron's four
  faces clipped against the target tetrahedron's four half-spaces in turn,
  each cut capped with a fan-triangulated, angle-sorted polygon — measured
  via the closed-triangulated-surface divergence theorem.

A target cell's value is the overlap-measure-weighted mean of every source
cell it intersects. A target cell whose covered fraction (against its own
unclipped measure) falls below a small relative tolerance is filled with
`default_value`; the count is reported through one aggregated warning per
call.

`cell_data` is transferred by this algorithm directly. `point_data` is
transferred by **composition**:
[point_data_to_cell_data](/data_operations) lumps the source array onto a
cell proxy (the unweighted mean of a cell's own corner values), the same
overlap algorithm remaps that proxy onto the target's cells exactly, and
[cell_data_to_point_data](/data_operations)`(weight="measure")` distributes
the result back onto the target's points. Conservation is **exact only for
the middle step** — the two lumping steps that sandwich it are each
already-documented approximations of their own (an unweighted corner mean,
and a measure-weighted point mean) — so the overall `point_data` path is a
layered approximation, not exact nodal/FEM conservation. No dual-cell /
control-volume machinery is built.

This operation deliberately reports **no** integral/conservation diagnostic
of its own — measuring how well conservation held on a given mesh is the job
of the roadmap's separate field-integration companion, not duplicated here.

There is deliberately no `extrapolate` flag (unlike `interpolate`): a silent
nearest-source-cell fallback for uncovered cells would break the exact
conservation guarantee for precisely the boundary cells a caller is most
likely to reach for it, on the one operation whose entire purpose is that
guarantee.

## Arrays and locations

`arrays=None` (the default) transfers **every source `point_data` and
`cell_data` array**, in sorted name order. A name present in both source
locations transfers both; a name present in neither raises.

## Name conflicts

If a transferred name already exists on the target, `on_conflict` decides:
`"error"` (default) raises, `"overwrite"` replaces the target's array, and
`"suffix"` writes to `name + "_interp"` — raising if that name is taken too,
so nothing is ever silently clobbered.

## Determinism

Output is byte-identical across the three mesh backends and across thread
counts (the same idioms `interpolate`/`data_average` already use):

- independent per-target-simplex work runs in parallel, with each simplex's
  own candidate accumulation in a fixed (deduplicated, ascending
  source-simplex index) order;
- the many-to-few scatter from simplices back onto the target's original
  cells runs **serially**, in ascending target-simplex index — floating-point
  addition is not associative, so this cannot be parallelised without making
  the result thread-count-dependent (the same rule
  `operations/data_average.cpp`'s cell-to-point averaging documents);
- the 3D clip kernel recentres on the source tetrahedron's own corner average
  before any arithmetic (the numerical-stability lesson `detail/polyhedron.hpp`'s
  `poly_measure` already documents), and its angle-sorted capping ring breaks
  an exact tie by ascending internal build index rather than insertion order.

**There is no pure-Python fallback.** The 3D clip kernel is a discrete-branch
geometric algorithm (half-space in/out classification, cutting-plane chord
deduplication, angle-sorted cap triangulation) of the same class
[subdivide](/subdivide)/[agglomerate](/agglomerate)/[decimate_volume](/decimate_volume)
already document as unsafe to give a second, independently-written
implementation — near a degenerate or tangent overlap the two could silently
disagree. `meshioplusplus.conservative_interpolate` raises
`NotImplementedError` by name when the compiled `_core` extension is
unavailable.

## CLI

```sh
meshioplusplus conservative-interpolate SOURCE TARGET OUT
meshioplusplus conservative-interpolate coarse.vtu fine.vtu out.vtu --arrays pressure
meshioplusplus conservative-interpolate a.msh b.msh out.vtu --on-conflict suffix
meshioplusplus conservative-interpolate a.msh b.msh out.vtu --default-value=-1
```

Available in both the Python CLI and the native `meshioplusplus` binary.

## Other languages

The operation is exposed on every binding surface:

```c
/* C API: arrays as char** + count (NULL / count <= 0 = all point_data AND cell_data) */
mio_mesh* out = mio_conservative_interpolate(source, target, NULL, 0, 0.0, "error");
```

```fortran
! Fortran: module-level function (two-mesh input, like mio_merge/mio_interpolate)
mapped = mio_conservative_interpolate(coarse, fine, stat=st)
```

```julia
# Julia
out = conservative_interpolate(coarse, fine; arrays=["pressure"])
```

```r
# R
out <- mio_conservative_interpolate(coarse, fine, arrays = "pressure")
```

```js
// WASM: arrays as a JS string array ([] = all point_data AND cell_data)
const out = m.conservativeInterpolate(coarse, fine, ['pressure'], 0, 'error');
```

The flat surfaces share the Python surface's semantics; as everywhere else on
the flat bindings, `point_sets`/`cell_sets` never cross (the target's sets
are re-attached by the Python shim only).

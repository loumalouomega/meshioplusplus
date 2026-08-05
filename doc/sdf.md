# Signed distance to a surface

How far a point is from a triangle skin, and which side of it the point is on.

```python
import meshioplusplus as mp

skin = mp.read("bunny.stl")

# distances at arbitrary points
d = mp.sample_distance(skin, [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]])

# or attached to a mesh, as ordinary data
grid = mp.grid([64, 64, 64], origin=(-1, -1, -1), spacing=(0.05,) * 3)
field = mp.distance_to_surface(grid, skin)
mp.write("field.vtu", field)              # carries sdf:distance as point_data
```

Negative is inside, by the usual convention. Because the result is ordinary
`point_data`, every writer carries it, `--color-by sdf:distance` colours it, and
`isosurface(field, "sdf:distance", [0.0])` reconstructs the surface.

## The sign is the hard part

An unsigned distance is a minimisation with no interesting failure mode. A
*signed* one has to decide which side of the surface a query is on, and the
standard mistake — taking the sign from the **nearest triangle's** normal — is
exactly right on convex geometry and exactly wrong wherever two nearly-opposite
faces meet.

The fix is the angle-weighted **pseudonormal** (Bærentzen & Aanæs): the sign comes
from the normal of the nearest *feature*, which is the triangle only when the
closest point lies in its interior:

| nearest feature | normal used |
|---|---|
| face | the triangle's own normal |
| edge | the sum of the two incident faces' unit normals |
| vertex | the angle-weighted sum of the incident faces' unit normals |

::: tip Where this actually matters
A *reentrant* corner does **not** expose the difference — at a concave edge both
incident faces happen to give the correct sign. The failure needs two incident
faces whose normals are nearly opposite, which is a **sharp spike**, not a notch.
`tests/cpp/test_surface_distance.cpp` uses a sliver prism whose tip subtends about
1.7°, and demonstrates the naive method getting it wrong rather than merely
asserting the right answer.
:::

## Signs on offer

| `sign` | behaviour | cost |
|---|---|---|
| `pseudonormal` (default) | exact for a watertight, consistently wound surface | free — it reuses the nearest-triangle query |
| `winding-number` | robust to holes, self-intersection and per-component orientation flips | **O(triangles) per query** |
| `unsigned` | no sign at all; the only meaningful mode on an open sheet | free |

`winding-number` has no acceleration structure short of a fast-multipole
expansion, so `max_winding_work` (default 2e9 on `n_queries × n_triangles`)
refuses by name rather than silently running for an hour. Pair it with a band.

## Is your surface even closed?

```python
mp.surface_watertight_check(skin)
# {'boundary_edges': 0, 'non_manifold_edges': 0, 'inconsistent_pairs': 0,
#  'degenerate_triangles': 0, 'watertight': True}
```

The counts are reported separately because they need different fixes — "your STL
has 12 boundary edges" is actionable, "not watertight" is not. `watertight_check`
is `"warn"` by default, and can be `"off"` or `"error"`.

## The band

`band=r` clamps distances beyond `r` and stops the search there, which is what
makes a narrow-band field cheap. Clamped values are marked in a **mandatory**
Int64 `sdf:band` array (`1` = computed, `0` = clamped): a clamped value is
byte-indistinguishable from a computed one, so without that array a caller cannot
tell the difference.

Inside the band, a banded run is **byte-identical** to a full-field one — the band
is an optimization, not a second answer, and a test pins that.

## The accelerator cannot change the answer

The nearest-triangle search runs over a bucket grid, and every candidate
comparison is totally ordered on `(squared distance, triangle id)`. That makes the
accelerator *provably unobservable*: the same query at any bucket size returns
byte-identical distances **and** the same nearest cell.

Two things follow, and both are load-bearing rather than incidental:

- `grid_cell_size` is a public option mainly so a test can prove the invariance
  (`test_the_bucket_size_does_not_change_the_answer`). It is also a legitimate
  tuning knob.
- The numpy reference implementation has **no accelerator at all** — a brute-force
  scan in ascending triangle order is the same computation, which is why the two
  are bit-identical.

The bucket-size heuristic was retuned once, after measurement, without
re-validating a single distance: sizing buckets by the mean triangle alone made a
64³ inside-fill of the 112k-triangle Stanford bunny take 19 seconds, because a
query far from the surface had to expand through hundreds of empty shells. Adding
the domain's own extent to the rule cut that to 2.9 seconds with **identical**
output.

## Where byte-parity stops

`sign="winding-number"` sums one `atan2` per triangle and compares the total
against a threshold. `atan2` is not correctly rounded, and its last-ulp behaviour
differs between libm and numpy, so the two implementations can genuinely disagree
on a near-tangent point. The numpy reference therefore **raises**
`NotImplementedError` rather than hiding the difference behind a tolerance — the
same call `_smooth.py` makes for its inversion guard and `_partition.py` for ghost
layers. The C++ path is still tested, against analytic fields and against the
pseudonormal answer.

Everything else here is `+ - * /`, comparisons and one correctly-rounded `sqrt`,
so it is bit-identical by construction.

## What the surface may be

`triangle` blocks are taken as they are; `quad` and rectangular `polygon` blocks
are fanned with the same deterministic fan `convert_cells(simplexify)` uses, so
the two cannot disagree about which diagonal a quad is split on. A volume or
polyhedron block is refused by name, pointing at `extract_surface`; a higher-order
block is refused pointing at `linearize`.

`sdf:closest_cell` (opt-in) names a cell of the mesh **you passed in**, not of the
triangulation you never saw.

## Not yet implemented

`compute_sdf` — generate a grid *and* fill it in one call, including an adaptive
octree — is declared but throws by name. Compose `voxelize` with
`distance_to_surface` in the meantime. Its option and result layouts are final
from v9.24.0 so that adding the implementation is a pure `.cpp` change; see
[`doc/abi.md`](abi.md) for why that matters.

## See also

- [Regular grids](voxelize.md) — where to evaluate the field.
- [Isosurfaces](isosurface.md) — contour it back to a surface.

# Decimation (QEM edge collapse)

`meshioplusplus.decimate(mesh, ratio=0.25)` reduces a **surface** mesh's face count by greedy quadric-error-metric edge collapse, preserving shape, boundaries and features. It is the resolution-*reducing* inverse of [refine](/refine), completing the pair, and — like the other mesh operations — uses only standard C++/numpy, so it runs under every mesh backend.

```python
import meshioplusplus

mesh = meshioplusplus.read("scan.stl")

# keep 25% of the faces
coarse = meshioplusplus.decimate(mesh, ratio=0.25)

# or stop at an absolute face count / an error budget
coarse = meshioplusplus.decimate(mesh, target_faces=5000)
coarse = meshioplusplus.decimate(mesh, max_error=1e-6)

coarse, report = meshioplusplus.decimate(mesh, ratio=0.25, return_report=True)
print(report["faces_removed"], report["collapses_rejected"])
```

Both CLIs expose it as a verb:

```bash
meshioplusplus decimate scan.stl coarse.stl --ratio 0.25
meshioplusplus decimate scan.stl coarse.stl --target-faces 5000 --placement midpoint
```

## The algorithm in a paragraph

Every vertex accumulates a 4×4 symmetric **quadric**: the sum, over its incident triangles, of the squared-distance-to-plane form of each triangle's plane, weighted by area (Garland & Heckbert, SIGGRAPH '97). Collapsing an edge merges its two endpoints into one surviving vertex whose *error* is the summed endpoint quadric evaluated at the surviving position — geometrically, the weighted squared distance to all the planes the two vertices have ever represented. The cheapest edge is collapsed, its quadrics are summed onto the survivor, the survivor's edges are re-scored, and the loop repeats until the stopping criterion is met. Cheap collapses are the ones that don't change the shape; expensive ones are at corners and creases — so error-ordered greedy collapse removes detail exactly where the surface is flat.

## Scope: surfaces in, triangles out

The operating type is `triangle`. `quad` and rectangular `polygon` blocks are **triangulated first** (via [`convert_cells`](/convert_cells)' simplexify), so the output is all-triangle even when the input was not; the block structure stays 1:1 with the input (a block may come back with zero faces).

Everything else raises **by name** rather than guessing:

- a mesh with **3D volume cells** points at `extract_surface` — QEM decimation of a *volume* mesh (tetra-collapse validity, no boundary-shape objective) is a different and much harder problem, deliberately out of scope rather than silently skinned;
- **higher-order** cells (`triangle6`, …) point at `linearize`;
- **ragged** polygon/polyhedron blocks;
- **`line`/`vertex`** blocks mixed in with the surface (their nodes would dangle after the collapse).

## Stopping criteria

Exactly one of the three must be given:

| criterion | meaning |
|---|---|
| `ratio` | fraction of the (triangulated) faces to KEEP, in `(0, 1]` |
| `target_faces` | absolute face count to stop at |
| `max_error` | collapse only while the cheapest candidate's quadric error is at most this (squared mesh units) |

A collapse removes one or two faces, so `ratio`/`target_faces` land **within one collapse** of the request. When pinning leaves no collapsible edge before the target is reached, the run warns and returns normally — the report tells the truth.

## Placement, and how data follows

`placement=` decides where the surviving vertex goes:

- **`optimal`** (default) — the minimizer of the summed quadric (a 3×3 solve), falling back to the edge midpoint when the system is ill-conditioned (`|det| ≤ 1e-12 · max|A|³`). On an *exactly planar* patch the system is singular by construction — the quadric is flat along the plane — so flat geometry effectively runs on midpoints; that is expected, not a defect.
- **`midpoint`** — always the edge midpoint.
- **`endpoint`** — the endpoint with the lower quadric error (tie → the lower vertex id).

Float-kind `point_data` at the survivor is blended between the two endpoints at the parameter `t` obtained by **projecting the placed point onto the edge, clamped to `[0, 1]`**. That `t` is exact for `midpoint`/`endpoint` and an *approximation* for `optimal` — the optimal point need not lie on the edge. Blended arrays keep their float dtype (the run accumulates in Float64 with one cast at the end). **Integer** `point_data` takes the survivor's own value — a blended material id is meaningless. Each surviving face keeps its own `cell_data` row; `field_data` passes through; `point_sets`/`cell_sets` are remapped through the returned maps in the Python shim.

## What is pinned

Mirroring [`smooth`](/smooth)'s vocabulary and defaults, a pinned vertex never moves and is never removed:

- **`preserve_boundary`** (default on) — boundary vertices, by the classic once-used-edge test (an edge used by exactly one triangle is boundary). An edge between two pinned vertices never enters the queue, which is what keeps an open patch's **outline exactly intact**; a collapse toward a pinned vertex keeps that vertex's own position regardless of `placement`.
- **`preserve_features`** (default on, `feature_angle=30`) — vertices where two incident face **normals** differ by more than the angle (the `vtkFeatureEdges` convention), which is what keeps a cube's corners and creases sharp. Conservative v1: crease vertices are fully pinned rather than allowed to slide along the crease.
- **`frozen`** — an optional index array, or the name of one of `mesh.point_sets`.

## Validity guards

Each guard **rejects the individual collapse** (counted in the report's `collapses_rejected`) rather than aborting:

- the **link condition** — the vertices adjacent to both endpoints must be exactly the opposite vertices of the faces on the edge, and an edge used by more than two faces is non-manifold and never collapsible. This is what guarantees the output cannot become non-manifold or change topology (no pinched surfaces, no merged sheets).
- **normal-flip rejection** — a surviving incident face whose normal would turn by 90° or more (`n_before · n_after ≤ 0`) rejects the collapse. Like `smooth`'s inversion guard this is strictly *do no harm*: a face that is already degenerate imposes no constraint, and a survivor collapsing to zero area is caught by the same test.

## Determinism

Setup is parallel with a fixed floating-point order (each vertex sums its quadric in ascending incident-face index); the greedy loop is **serial** — QEM is inherently sequential — driven by a priority queue with the total order *(error, lower endpoint id, higher endpoint id)* and lazy version-stamped deletion. Output is **byte-identical** across the three mesh backends, across thread counts, and across the C++-core/numpy-fallback boundary: the pure Python fallback is an expression-for-expression twin, pinned by `tests/python/test_decimate.py::test_cpp_matches_python`.

## Other language surfaces

- **C API** — `mio_decimate(...)` returning an opaque `mio_decimate_result` (mesh borrow/take, zero-copy point/cell maps, counter getters). The `frozen` mask is not exposed across the C ABI (a documented flat-ABI gap, like `mio_smooth`'s).
- **Fortran** — `coarse = m%decimate(ratio=0.25_real64, faces_removed=n, ...)` with optional counter out-args and a 1-based `point_map`.
- **WASM** — `decimate(mesh, ratio, targetFaces, maxError, placement, ...)` → `{mesh, facesRemoved, pointsRemoved, collapsesRejected, maxErrorApplied}`, also available as a `convertSurfaceOps` pipeline op (`{op: "decimate"}`, defaulting to `ratio: 0.5`).
- **CLI** — the `decimate` verb in both the Python and the native CLI (see [CLI](/cli)).

The returned maps make the result composable: `point_map` sends every input point to its **survivor's** output index (collapsed points map to the survivor, not −1 — usable for remapping external per-point arrays), and the per-input-block `cell_maps` send each input cell to its first surviving triangle (−1 when none survived).

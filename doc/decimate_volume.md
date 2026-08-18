# Volume decimation (QEM tet-edge collapse)

`meshioplusplus.decimate_volume(mesh, ratio=0.25)` reduces a **tetrahedral**
mesh's cell count by greedy quadric-error-metric tet-edge collapse. It is the
volume-mesh sibling of [`decimate`](/decimate) — a **separate operation**, not
a mode on it: `decimate` keeps raising by name on any 3D volume block,
pointing here, and this operation never touches `decimate`'s own machinery.
Like the other mesh operations, it uses only standard C++, so it runs under
every mesh backend.

```python
import meshioplusplus

mesh = meshioplusplus.read("solid.vtu")

# keep 25% of the tets
coarse = meshioplusplus.decimate_volume(mesh, ratio=0.25)

# or stop at an absolute tet count / an error budget
coarse = meshioplusplus.decimate_volume(mesh, target_cells=5000)
coarse = meshioplusplus.decimate_volume(mesh, max_error=1e-6)

coarse, report = meshioplusplus.decimate_volume(
    mesh, ratio=0.25, return_report=True
)
print(report["tets_removed"], report["collapses_rejected"])
```

Both CLIs expose it as a verb:

```bash
meshioplusplus decimate-volume solid.vtu coarse.vtu --ratio 0.25
meshioplusplus decimate-volume solid.vtu coarse.vtu --target-cells 5000 --placement midpoint
```

## Boundary vertices participate

Unlike `decimate`, whose `preserve_boundary` defaults **on** (boundary
vertices never move), `decimate_volume`'s `preserve_boundary` defaults
**off** — boundary vertices are scored and collapsed by a real quadric-error
objective, exactly like every other vertex. Pass `preserve_boundary=True` to
reproduce `decimate`'s own pinned-boundary behaviour instead.

## The algorithm in a paragraph

Every vertex accumulates a Garland–Heckbert plane **quadric** — but only from
its incident **boundary triangles** (the mesh's own outer skin, found via the
mesh's owner/neighbour face dual). A purely interior vertex, by construction,
touches no boundary triangle at all, so its quadric is the exact zero matrix,
not merely small. For a candidate edge `(a, b)` the **combined quadric**
`Q(a) + Q(b)` feeds the *same* placement solve `decimate` uses:

- if at least one endpoint touches the boundary, `Q_ab` is non-degenerate (or
  the solve's own ill-conditioning bound routes it to the midpoint) —
  identical to `decimate`'s own code path;
- if **both** endpoints are purely interior, `Q_ab` is exactly zero, which
  that *same* bound classifies as degenerate and falls back to the
  **midpoint** — no interior-specific placement logic is needed at all.

Scoring needs one explicit split, though: a quadric error evaluated against an
exact-zero quadric is identically zero for *every* interior-interior edge,
which would tie the whole interior of the mesh and fall through to id order
rather than a useful ordering. Every candidate therefore carries a **regime**:
`0` (boundary-touching, scored by real quadric error) or `1` (purely
interior, scored by squared edge length) — regime-0 collapses are always
considered ahead of regime-1 ones in the greedy queue.

## Scope: tets in, tets out

The operating type is `tetra`, exclusively. Everything else raises **by
name** rather than guessing:

- a non-tetra **3D volume cell** (hexahedron, wedge, pyramid) or a
  **polyhedron/ragged** block points at `convert_cells(mode="simplexify")`;
- any **non-3D block** mixed in with the tets (its nodes would dangle after
  the collapse) points at dropping it first, e.g. via `split`.

The block structure stays 1:1 with the input (a block may come back with
zero cells).

## Stopping criteria

Exactly one of the three must be given:

| criterion | meaning |
|---|---|
| `ratio` | fraction of the tets to KEEP, in `(0, 1]` |
| `target_cells` | absolute tet count to stop at |
| `max_error` | collapse only while the cheapest **boundary-touching** candidate's quadric error is at most this (squared mesh units) |

`max_error` is only strictly meaningful for regime-0 (boundary-touching)
entries — real quadric-error units. Compared directly against a regime-1
entry's squared edge length it is a documented rough tool for mixed meshes,
not a claimed exact criterion; `ratio`/`target_cells` are regime-agnostic and
remain exactly well-defined either way. A collapse removes every tet sharing
the edge (one or more), so `ratio`/`target_cells` land **within one
collapse** of the request. When pinning leaves no collapsible edge before the
target is reached, the run warns and returns normally.

## Placement, and how data follows

`placement=` decides where the surviving vertex goes — `optimal` (default,
the quadric minimizer, midpoint when ill-conditioned or purely interior),
`midpoint`, or `endpoint` (the endpoint with the lower quadric error) —
exactly `decimate`'s own three modes, reused unchanged. A boundary vertex's
quadric is built from boundary-triangle planes only, so its minimizer
naturally stays near the surface.

Float-kind `point_data` at the survivor is blended between the two endpoints
at the clamped edge-projection parameter, exactly as `decimate` does.
**Integer** `point_data` keeps the survivor's own value. Each surviving tet
keeps its own `cell_data` row; `field_data` passes through; named regions
(and so `point_sets`/`cell_sets`) are carried by the C++ core itself — Point
and Cell regions survive, named Side regions do not (a removed tet's local
facet numbering has no correspondence to the survivor's).

## What is pinned

Mirroring `decimate`'s vocabulary, a pinned vertex never moves and is never
removed:

- **`preserve_boundary`** (default **off**, see above) — every boundary
  vertex, by the once-used-face test on the mesh's own outer skin.
- **`preserve_features`** (default on, `feature_angle=30`) — boundary
  vertices whose incident **boundary-triangle** normals differ by more than
  the angle, keeping corners and creases of the outer surface sharp.
- **`frozen`** — an optional index array, or the name of one of
  `mesh.point_sets`.

## Validity guards

Tet-only, so simpler than the general 3D case: a tet's *other* incident cells
are exactly the tets shared by both endpoints. Each guard **rejects the
individual collapse** (counted in the report's `collapses_rejected`) rather
than aborting:

- the **vertex-link condition** — the exact set of nodes adjacent to both
  endpoints (via their own incident tets) must equal the "opposite corners"
  of every shared tet — pure integer set equality, no floating point;
- the **duplicate-tet guard** — no surviving tet incident to one endpoint
  alone may share its "opposite triangle" with a surviving tet incident to
  the other endpoint alone (which would otherwise produce two tets with
  identical corners);
- **tet-inversion rejection** (`smooth`'s "do no harm") — a surviving tet
  whose signed volume is non-zero must not change sign under the candidate
  placement; an already-degenerate tet imposes no constraint.

**Boundary-touching collapses additionally run `decimate`'s own** ring/
shared-face link condition and normal-flip check, over the mesh's own outer
skin — so the boundary surface cannot tear, pinch or fold independently of
the interior guards above.

## Determinism

Setup is parallel with a fixed floating-point order, matching `decimate`'s
own; the greedy loop is **serial**, driven by a priority queue with lazy
version-stamped deletion. **This operation is C++-core only, with no numpy
fallback at all** — like [`subdivide`](/subdivide)/[`agglomerate`](/agglomerate)
and unlike `decimate` itself. The volume-specific validity guards sit on top
of the mesh's own outer skin, whose incremental maintenance through hundreds
of collapses is exactly the kind of intricate bookkeeping a second,
independently written implementation would risk silently diverging from.
Calling `decimate_volume` without the compiled `meshioplusplus._core`
extension raises `NotImplementedError` naming the reason.

## Other language surfaces

- **C API** — `mio_decimate_volume(...)` returning an opaque
  `mio_decimate_volume_result` (mesh borrow/take, zero-copy point/cell maps,
  counter getters) — its own result type, distinct from `mio_decimate_result`.
  The `frozen` mask is not exposed across the C ABI (a documented flat-ABI
  gap, like `mio_decimate`'s own).
- **CLI** — the `decimate-volume` verb in both the Python and the native CLI
  (see [CLI](/cli)); `--preserve-boundary` is an opt-**in** flag here (the
  opposite of `decimate`'s opt-out `--no-preserve-boundary`).
- **Pipeline** — a `DecimateVolume` step (`{Ratio, TargetCells, MaxError,
  Placement, PreserveBoundary, PreserveFeatures, FeatureAngle}`), dispatched
  generically off the shared op table.
- **MCP** — the `decimate_volume` tool.

The returned maps make the result composable: `point_map` sends every input
point to its **survivor's** output index (collapsed points map to the
survivor, not −1), and the per-input-block `cell_maps` send each input cell
to its own output index, or −1 when it did not survive.

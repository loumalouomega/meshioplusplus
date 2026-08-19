# Surface remeshing (ACVD clustering)

`meshioplusplus.remesh(mesh, num_clusters)` replaces a **surface** mesh's
triangulation with a new, near-uniformly-sized, well-shaped one at a
caller-chosen vertex count, by approximated centroidal Voronoi diagram (ACVD)
clustering. It is the one resolution-changing operation in this repo that
does not work on the input's own triangulation: [`refine`](/refine)
subdivides the input's cells, [`decimate`](/decimate) collapses its edges,
[`subdivide`](/subdivide)/[`agglomerate`](/agglomerate) restructure it into
polyhedra, [`smooth`](/smooth) moves its points — all of them inherit the
input's element shapes. QEM decimation in particular can only *remove*
elements, so a badly-shaped input stays badly-shaped at every target count.
`remesh` instead partitions the surface into `num_clusters` compact regions
and builds the **dual** of that partition, so output element quality is a
property of the clustering rather than of the input. Like the other mesh
operations, it uses only standard C++, so it runs under every mesh backend.

```python
import meshioplusplus

mesh = meshioplusplus.read("bracket.stl")

out = meshioplusplus.remesh(mesh, num_clusters=5000)

# the feature-preserving metric keeps sharp edges and corners
out = meshioplusplus.remesh(mesh, num_clusters=5000, metric="quadric")

out, report = meshioplusplus.remesh(mesh, num_clusters=5000, return_report=True)
print(report["num_clusters"], report["num_isolated_clusters"])
```

Both CLIs expose it as a verb:

```bash
meshioplusplus remesh bracket.stl out.vtu --num-clusters 5000
meshioplusplus remesh bracket.stl out.vtu --num-clusters 5000 --metric quadric
```

## Attribution and licence reasoning

The reference method is Valette's ACVD, but ACVD itself is
[CeCILL-B](https://github.com/valette/ACVD) and hard-bound to VTK, both
incompatible with this project's MIT/dependency-free posture — its source is
never read or vendored here. The isotropic clustering engine below is instead
derived from [`pyvista/pyacvd`](https://github.com/pyvista/pyacvd)
(`src/clustering.cpp`/`src/pyacvd/clustering.py`, nanobind), **MIT License,
Copyright (c) 2017-2024 The PyVista Developers** — itself an independent
implementation of the published research of S. Valette and J.-M. Chassery,
not of ACVD's own code, and whose algorithm is plain-array C++ with no VTK at
all. Deriving from it is therefore a clean MIT→MIT chain, requiring only a
preserved copyright notice (see `CITATION.cff`) — not the CeCILL-B
interface/website credit obligations ACVD's own licence would impose.

`metric="quadric"` (the feature-preserving variant, "ACVDQ" in the
literature) has **no pyacvd counterpart** — pyacvd implements only the
isotropic method. It is a fresh synthesis over this project's own
pre-existing Garland-Heckbert quadric machinery
(`detail/decimate_common.hpp`, shared with `decimate`/`decimate_volume`),
built from the general "metric-dependent discrete Voronoi diagram" concept
(Valette, Chassery & Prost, IEEE TVCG 2008) — never from ACVD's own
`vtkQEMetricForClustering.h`, which this project does not read.

## The algorithm in a paragraph

Each input vertex is an "item" carrying a weight (a third of its incident
triangle area, so the weights sum to the surface area) and a weighted
position. A cluster holds only two accumulators — `sgamma = sum(w * x)` and
`srho = sum(w)` — which is what makes testing a candidate move O(1) and the
whole method fast:

- **Seeding** grows `num_clusters` regions breadth-first from the lowest
  unassigned vertex, each absorbing neighbours until it reaches its share of
  the remaining area. Deliberately RNG-free (unlike pyacvd), so a given mesh
  always yields the same clustering.
- **Energy minimisation** sweeps the edges whose endpoints lie in different
  clusters and, for each, compares the current badness against moving either
  endpoint across, accepting whichever move lowers the total. A cluster
  unmodified in the previous sweep is skipped, so late sweeps are cheap.
- **Repair** is not optional: seeding can leave vertices unassigned and
  minimisation can split a cluster into disconnected pieces. Unassigned
  vertices are grown into from their neighbours; a split cluster keeps its
  largest component and the rest are re-grown, alternating with further
  minimisation up to `max_repair_passes`.
- **The dual** places one output vertex at each cluster's representative
  point and emits one triangle for every input triangle whose three vertices
  lie in three *distinct* clusters, deduplicated. Winding is fixed by
  comparing each output triangle's normal against the mean input normal of
  the clusters it joins.

### Metric: isotropic vs quadric

Both metrics maximize cluster compactness — the acceptance test is always a
"badness to minimize" derived from the same `-|sgamma|^2/srho` form — but
`metric="quadric"` **adds** a Garland-Heckbert quadric flatness term on top,
and changes the representative point:

- **`isotropic`** (default) — representative point = area-weighted centroid.
  Fast and robust; rounds sharp edges and corners toward the interior, since
  a corner cluster's centroid is pulled off the true surface.
- **`quadric`** — each cluster also accumulates a 10-entry quadric (summed
  from its members' per-vertex quadrics), and the representative point is
  that quadric's own minimizer (falling back to the isotropic centroid when
  the solve is ill-conditioned — an exactly flat region, the common case). A
  flat cluster's quadric is near-singular in-plane and degenerates to the
  centroid; a crease or corner's quadric is not, and pins the point onto the
  feature instead of smearing it toward a centroid.

The quadric term is **added to**, not substituted for, the isotropic
compactness term during clustering itself — a pure-quadric objective has no
compactness pull at all (a cluster confined to a common plane can grow
arbitrarily thin and snake-like along a low-curvature direction without its
quadric error rising), which was measured to repeatedly produce disconnected,
non-manifold clusters even after every repair pass. The additive combination
keeps the compactness guarantee while still letting the quadric term decide
where it actually matters — near a real feature.

## Subdivision matters more than iteration count

The clustering is a *discrete* approximation of a centroidal Voronoi diagram
whose resolution is the number of items per cluster, so a target near the
input's own vertex count gives poor results until the input is refined.
`subdivide` applies that many uniform [`refine`](/refine) passes first (each
multiplying the triangle count by four); leaving it unset (the default)
auto-picks the smallest count reaching `subsample_ratio` (default 10) items
per cluster, capped at `max_subdivide` (default 4). `subdivide=0` disables
subdivision entirely.

## Scope: surface in, triangles out

The operating type is `triangle`; `quad` and rectangular `polygon` blocks are
triangulated first via `convert_cells(mode="simplexify")`. Everything else
raises **by name**, mirroring `decimate`'s own scope exactly:

- a **3D volume cell** points at `extract_surface` (the volumetric CVD/ODT
  problem is genuinely different, not a silent skin);
- a **higher-order** cell points at `convert_cells(mode="linearize")`;
- a **ragged** polygon/polyhedron block, or a `line`/`vertex` block, raises
  outright.

## What does NOT survive

The output is a brand-new mesh with new points and new connectivity, so —
unlike every other operation in this layer — there is **no meaningful point
or cell map**, and `point_data`, `cell_data` and named regions are
**dropped**. Transfer a field onto the result with
[`interpolate`](/interpolate) or
[`conservative_interpolate`](/conservative_interpolate), which is exactly the
composition those operations exist for. `field_data` passes through
verbatim.

## Boundaries are not specially protected

An open surface's boundary vertices are ordinary items with no extra pinning
or dual-edge insertion. `remesh` runs on an open mesh without error (matching
the reference method's own default), but the outline near an open edge is
not guaranteed preserved — boundary/manifoldness protection remains a
documented, still-open item on `doc/roadmap.md`.

## Topology is not preserved

Two surface sheets closer together than a cluster can merge, and the genus
can change. This is inherent to dualising a discrete Voronoi partition, not
an implementation limit — a target vertex count far below the feature scale
will visibly simplify topology.

## Output manifoldness is best-effort

The dual of a discrete Voronoi partition need not be 2-manifold. Repair
removes the common cause (disconnected clusters), and the report's
`num_isolated_clusters` names what could not be fixed — check it rather than
assuming; a non-zero value means the output may be non-manifold near those
clusters.

## Determinism, and no numpy fallback

Seeding is RNG-free, the energy sweep visits edges in a fixed order and is
**serial** (the objective is inherently sequential — every accepted move
changes the accumulators the next test reads), and the setup passes fill
disjoint slots in parallel with a fixed floating-point order. Output is
therefore byte-identical across the three mesh backends and thread counts —
but there is deliberately **no numpy twin**, like
[`subdivide`](/subdivide)/[`agglomerate`](/agglomerate)/[`decimate_volume`](/decimate_volume)/[`conservative_interpolate`](/conservative_interpolate):
a single near-tie decided differently by an independent implementation
diverges into a genuinely different clustering, not a last-ulp difference.
Calling `remesh` without the compiled `meshioplusplus._core` extension raises
`NotImplementedError` naming the reason, for any input.

## Other language surfaces

- **C API** — `mio_remesh(mesh, num_clusters, subdivide, subsample_ratio,
  max_subdivide, max_iterations, max_repair_passes, metric,
  &num_clusters_out, &num_iterations, &subdivide_applied,
  &num_isolated_clusters)` → a plain `mio_mesh*` (the counter out-params are
  all nullable, `mio_smooth`'s shape — no opaque result handle, since there
  are no maps to hand back).
- **Fortran** — `m%remesh(num_clusters, ...)`, the same optional-argument
  shape as `m%estimate_error`.
- **Julia** — `remesh(mesh, num_clusters; subdivide=-1, metric=:isotropic,
  ...)` → a `(; mesh, num_clusters, num_iterations, subdivide_applied,
  num_isolated_clusters)` NamedTuple.
- **R** — `mio_remesh(mesh, num_clusters, ...)` → a named list of the same
  five fields (counters as `double`, R having no native int64).
- **WASM** — `remesh(mesh, numClusters, subdivide, subsampleRatio,
  maxSubdivide, maxIterations, maxRepairPasses, metric)`, and reachable as a
  `{op: 'remesh', numClusters: ...}` `convertSurfaceOps`/pipeline step.
- **CLI** — the `remesh` verb in both the Python and the native CLI (see
  [CLI](/cli)): `--num-clusters` (required), `--subdivide`,
  `--subsample-ratio`, `--max-subdivide`, `--iterations`, `--repair-passes`,
  `--metric isotropic|quadric`.
- **Pipeline** — a `Remesh` step (`{NumClusters, Subdivide, SubsampleRatio,
  MaxSubdivide, MaxIterations, MaxRepairPasses, Metric}`), dispatched
  generically off the shared op table — the one step whose output has no
  correspondence to its input, exactly like a fresh mesh from `Voxelize`.
- **MCP** — the `remesh` tool.

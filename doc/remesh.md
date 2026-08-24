# Surface remeshing (ACVD clustering)

`meshioplusplus.remesh(mesh, num_clusters)` replaces a **surface** mesh's triangulation with a new, near-uniformly-sized, well-shaped one at a caller-chosen vertex count, by approximated centroidal Voronoi diagram (ACVD) clustering. It is the one resolution-changing operation in this repo that does not work on the input's own triangulation: [`refine`](/refine) subdivides the input's cells, [`decimate`](/decimate) collapses its edges, [`subdivide`](/subdivide)/[`agglomerate`](/agglomerate) restructure it into polyhedra, [`smooth`](/smooth) moves its points — all of them inherit the input's element shapes. QEM decimation in particular can only *remove* elements, so a badly-shaped input stays badly-shaped at every target count. `remesh` instead partitions the surface into `num_clusters` compact regions and builds the **dual** of that partition, so output element quality is a property of the clustering rather than of the input. Like the other mesh operations, it uses only standard C++, so it runs under every mesh backend.

```python
import meshioplusplus

mesh = meshioplusplus.read("bracket.stl")

out = meshioplusplus.remesh(mesh, num_clusters=5000)

# the feature-preserving metric keeps sharp edges and corners
out = meshioplusplus.remesh(mesh, num_clusters=5000, metric="quadric")

# the anisotropic metric elongates elements along low-curvature directions
out = meshioplusplus.remesh(mesh, num_clusters=5000, metric="anisotropic", max_anisotropy=4.0)

out, report = meshioplusplus.remesh(mesh, num_clusters=5000, return_report=True)
print(report["num_clusters"], report["num_isolated_clusters"])
```

Both CLIs expose it as a verb:

```bash
meshioplusplus remesh bracket.stl out.vtu --num-clusters 5000
meshioplusplus remesh bracket.stl out.vtu --num-clusters 5000 --metric quadric
meshioplusplus remesh bracket.stl out.vtu --num-clusters 5000 --metric anisotropic --max-anisotropy 4
```

## Attribution and licence reasoning

The reference method is Valette's ACVD, but ACVD itself is [CeCILL-B](https://github.com/valette/ACVD) and hard-bound to VTK, both incompatible with this project's MIT/dependency-free posture — its source is never read or vendored here. The isotropic clustering engine below is instead derived from [`pyvista/pyacvd`](https://github.com/pyvista/pyacvd) (`src/clustering.cpp`/`src/pyacvd/clustering.py`, nanobind), **MIT License, Copyright (c) 2017-2024 The PyVista Developers** — itself an independent implementation of the published research of S. Valette and J.-M. Chassery, not of ACVD's own code, and whose algorithm is plain-array C++ with no VTK at all. Deriving from it is therefore a clean MIT→MIT chain, requiring only a preserved copyright notice (see `CITATION.cff`) — not the CeCILL-B interface/website credit obligations ACVD's own licence would impose.

`metric="quadric"` (the feature-preserving variant, "ACVDQ" in the literature) has **no pyacvd counterpart** — pyacvd implements only the isotropic method. It is a fresh synthesis over this project's own pre-existing Garland-Heckbert quadric machinery (`detail/decimate_common.hpp`, shared with `decimate`/`decimate_volume`), built from the general "metric-dependent discrete Voronoi diagram" concept (Valette, Chassery & Prost, IEEE TVCG 2008) — never from ACVD's own `vtkQEMetricForClustering.h`, which this project does not read.

`metric="anisotropic"` has no pyacvd counterpart either and cites the same TVCG 2008 paper for the same reason: a metric-weighted Voronoi diagram is the paper's whole subject, and per-vertex curvature-tensor metrics are its headline anisotropic application, not something specific to ACVD's own implementation. It is built entirely from this project's own osculating-paraboloid curvature fit (already present for `gradation`, see below) and the same quadric-accumulator machinery `metric="quadric"` uses — see "Curvature-anisotropic remeshing" below for how the two connect.

## The algorithm in a paragraph

Each input vertex is an "item" carrying a weight (a third of its incident triangle area, so the weights sum to the surface area) and a weighted position. A cluster holds only two accumulators — `sgamma = sum(w * x)` and `srho = sum(w)` — which is what makes testing a candidate move O(1) and the whole method fast:

- **Seeding** grows `num_clusters` regions breadth-first from the lowest unassigned vertex, each absorbing neighbours until it reaches its share of the remaining area. Deliberately RNG-free (unlike pyacvd), so a given mesh always yields the same clustering.
- **Energy minimisation** sweeps the edges whose endpoints lie in different clusters and, for each, compares the current badness against moving either endpoint across, accepting whichever move lowers the total. A cluster unmodified in the previous sweep is skipped, so late sweeps are cheap.
- **Repair** is not optional: seeding can leave vertices unassigned and minimisation can split a cluster into disconnected pieces. Unassigned vertices are grown into from their neighbours; a split cluster keeps its largest component and the rest are re-grown, alternating with further minimisation up to `max_repair_passes`.
- **The dual** places one output vertex at each cluster's representative point and emits one triangle for every input triangle whose three vertices lie in three *distinct* clusters, deduplicated. Winding is fixed by comparing each output triangle's normal against the mean input normal of the clusters it joins.

### Metric: isotropic, quadric or anisotropic

All three metrics maximize cluster compactness — the acceptance test is always a "badness to minimize" derived from the same `-|sgamma|^2/srho` form — but `metric="quadric"` **adds** a Garland-Heckbert quadric flatness term on top and `metric="anisotropic"` **replaces** it with a curvature-tensor error term (both change the representative point too):

- **`isotropic`** (default) — representative point = area-weighted centroid. Fast and robust; rounds sharp edges and corners toward the interior, since a corner cluster's centroid is pulled off the true surface.
- **`quadric`** — each cluster also accumulates a 10-entry quadric (summed from its members' per-vertex quadrics), and the representative point is that quadric's own minimizer (falling back to the isotropic centroid when the solve is ill-conditioned — an exactly flat region, the common case). A flat cluster's quadric is near-singular in-plane and degenerates to the centroid; a crease or corner's quadric is not, and pins the point onto the feature instead of smearing it toward a centroid.
- **`anisotropic`** — clusters shaped by a local curvature tensor rather than isotropic distance, so elongated features (a fillet, a pipe, a rib) are meshed with elongated elements at a fraction of the vertex count an isotropic run would need. See "Curvature-anisotropic remeshing" below.

The quadric term is **added to**, not substituted for, the isotropic compactness term during clustering itself — a pure-quadric objective has no compactness pull at all (a cluster confined to a common plane can grow arbitrarily thin and snake-like along a low-curvature direction without its quadric error rising), which was measured to repeatedly produce disconnected, non-manifold clusters even after every repair pass. The additive combination keeps the compactness guarantee while still letting the quadric term decide where it actually matters — near a real feature. The anisotropic term does not need this stabilizer: unlike a flat quadric, a curvature tensor is always SPD (full rank), so it already penalizes spread in every direction — just anisotropically — and is used pure.

## Subdivision matters more than iteration count

The clustering is a *discrete* approximation of a centroidal Voronoi diagram whose resolution is the number of items per cluster, so a target near the input's own vertex count gives poor results until the input is refined. `subdivide` applies that many uniform [`refine`](/refine) passes first (each multiplying the triangle count by four); leaving it unset (the default) auto-picks the smallest count reaching `subsample_ratio` (default 10) items per cluster, capped at `max_subdivide` (default 4). `subdivide=0` disables subdivision entirely.

## Scope: surface in, triangles out

The operating type is `triangle`; `quad` and rectangular `polygon` blocks are triangulated first via `convert_cells(mode="simplexify")`. Everything else raises **by name**, mirroring `decimate`'s own scope exactly:

- a **3D volume cell** points at `extract_surface` (the volumetric CVD/ODT problem is genuinely different, not a silent skin);
- a **higher-order** cell points at `convert_cells(mode="linearize")`;
- a **ragged** polygon/polyhedron block, or a `line`/`vertex` block, raises outright.

## What does NOT survive

The output is a brand-new mesh with new points and new connectivity, so — unlike every other operation in this layer — there is **no meaningful point or cell map**, and `point_data`, `cell_data` and named regions are **dropped**. Transfer a field onto the result with [`interpolate`](/interpolate) or [`conservative_interpolate`](/conservative_interpolate), which is exactly the composition those operations exist for. `field_data` passes through verbatim.

## Curvature gradation

`gradation` (the exponent `gamma` in the item weight `area * kappa**gamma`) concentrates clusters where the surface bends more sharply, instead of spreading them by area alone. `kappa` is a per-vertex curvature magnitude from a local osculating-paraboloid fit over the vertex's 1-ring: fit `h = a*u**2 + b*u*w + c*w**2` in a local tangent frame, then `kappa = max(|kappa1|, |kappa2|)` from the eigenvalues of `[[2a, b], [b, 2c]]`. A near-singular fit (low valence, near-collinear neighbourhood) falls back to `kappa = 0`. `gradation = 0.0` (the default) means `kappa**0 == 1` identically, so curvature is not even computed and every weight reproduces plain area weighting byte-for-byte — every pre-existing test and every closed-mesh example in this repo is unaffected. Positive `gradation` pulls resolution toward high-curvature regions (a bump, a fillet, a sharp bend); it applies identically under both `metric="isotropic"` and `metric="quadric"` — gradation is orthogonal to the metric choice, not a third one.

## Curvature-anisotropic remeshing

`metric="anisotropic"` shapes each cluster with a per-vertex curvature *tensor* instead of a scalar distance, so a cluster is free to be long in a low-curvature direction and short in a high-curvature one — a fillet, a pipe or a rib is then meshed with elongated triangles at a fraction of the vertex count an isotropic (or quadric) run needs to resolve the same feature. `max_anisotropy` (default `4.0`) caps how elongated a single vertex's metric is allowed to be.

**Reuses the gradation fit, does not duplicate it.** The same osculating-paraboloid fit described above (`h = a*u**2 + b*u*w + c*w**2` in a local tangent frame) already yields the two principal curvatures `kappa1`/`kappa2` as the eigenvalues of `[[2a, b], [b, 2c]]`; gradation keeps only their magnitude. The anisotropic metric additionally keeps the corresponding **eigen*vectors***, i.e. the two principal *directions* in the tangent plane, mapped back into world coordinates. Neither the sign of an eigenvalue nor a consistent choice between the two eigenvectors matters here — only the outer product `e ⊗ e` ever enters the metric, and `(-e) ⊗ (-e) == e ⊗ e`, so there is no tie-break to get wrong or to keep consistent between neighbouring vertices.

**Metric construction, per vertex:**

1. Target edge lengths along each principal direction are `~ 1/sqrt(|kappa|)` (short across curvature, long along it); their ratio is clamped to `max_anisotropy` before anything else is computed, which is what keeps a single near-degenerate vertex from producing an unusably long or thin cluster.
2. The eigenvalue along the surface **normal** is pinned to the tighter (higher-curvature) of the two in-plane eigenvalues — not left at some independent value — which is what keeps the resulting ellipsoidal metric ball tangent to the true surface rather than letting a cluster drift off it in the normal direction.
3. The resulting 3x3 symmetric positive-definite tensor is then uniformly rescaled so its determinant is exactly `1`. This is a **reshape, never a rescale**: the metric only ever changes an item's relative in-cluster distances, never the overall element size, which is what keeps `max_anisotropy` orthogonal to both `num_clusters` (global element count) and `gradation` (density) rather than interacting with either.
4. A perfectly flat vertex (zero curvature in both directions, the degenerate case the ratio clamp already guards) yields **exactly** the identity tensor — `metric="anisotropic"` on an exactly flat mesh is therefore bit-for-bit identical to `metric="isotropic"`.

**Packs into the existing quadric accumulator, not a new one.** With a per-vertex SPD metric `M`, the anisotropic clustering energy `E = sum(w * (x - v)^T M (x - v))` expands algebraically into `v^T A v + 2 b.v + c` with `A = sum(w*M)`, `b = -sum(w*M*x)`, `c = sum(w * x^T M x)` — which is *exactly* the same 10-double Garland-Heckbert quadric layout (`[aa,ab,ac,ad,bb,bc,bd,cc,cd,dd]`) `metric="quadric"` already accumulates per cluster, minimised by the same solver. Choosing `metric="anisotropic"` therefore needs no new per-vertex or per-cluster storage at all — it changes only which 10 doubles a vertex contributes, never how they are accumulated, minimised or read back.

**The badness term is used pure, unlike `metric="quadric"`'s additive stabilizer.** `metric="quadric"`'s clustering acceptance test *adds* the quadric error to the isotropic compactness term (see above), because a flat quadric is rank-1 and does not penalise a cluster snaking along its own null direction. A curvature-tensor metric has no such null direction — it is SPD by construction (the normal-pinning step above guarantees full rank) — so it already penalises spread in *every* direction, just anisotropically, up to the `max_anisotropy` ratio. Using it as an additive stabilizer on top of isotropic compactness would only double-count that penalty, so the anisotropic badness is the pure quadric-error form.

`max_anisotropy` set to a non-default value under any metric other than `"anisotropic"` is a **named error**, not a silently-ignored option — the same posture `write_options.hpp` and selective `refine`'s mutually-exclusive selectors already take.

## Boundaries and output manifoldness

`preserve_boundary` (default `True`, and so free on the closed meshes most examples and every pre-existing test use) detects the input's open boundary, if any, via a plain edge-use-count pass — a triangle edge used by exactly one face is a boundary edge. Boundary vertices are then seeded **before** the interior BFS, so a cluster's boundary segment stays anchored to a contiguous stretch of the outline rather than being absorbed piecemeal by whichever interior cluster reaches it first. The dual gains a companion pass: every boundary edge whose two endpoints land in different clusters emits one `line` cell between those clusters' representative points, so the output gains a **second, optional `line` cell block** carrying the boundary polyline alongside the `triangle` block. Without this, a boundary-adjacent triangle simply has no "third neighbour" across the missing side and is silently dropped from the triangle dual, leaving no coherent output boundary at all.

This is a clean-room design achieving the roadmap's stated boundary goals (built from this project's own existing conventions — the pinning idiom `decimate`/`smooth` already use, the dual-edge dedup idiom the triangle dual already uses), not a reproduction of ACVD's own boundary-fixing algorithm, whose CeCILL-B source this project does not read (see the licence reasoning above).

## Topology is not preserved

Two surface sheets closer together than a cluster can merge, and the genus can change. This is inherent to dualising a discrete Voronoi partition, not an implementation limit — a target vertex count far below the feature scale will visibly simplify topology.

## Output manifoldness is best-effort

The dual of a discrete Voronoi partition need not be 2-manifold, for two distinct reasons, reported separately. Disconnected clusters are `num_isolated_clusters`' concern, unchanged since the operation's first release. Non-manifold **output vertices** — a vertex whose incident dual-triangle fan does not form a single loop (interior) or open chain (boundary), i.e. a "bowtie" — are `num_non_manifold_vertices`'s concern; both are folded into the same repair loop (regrow the implicated clusters, minimise again, up to `max_repair_passes`) rather than two separate loops. Either counter non-zero means a pathological input still produced non-manifold output — check both rather than assuming.

## Determinism, and no numpy fallback

Seeding is RNG-free, the energy sweep visits edges in a fixed order and is **serial** (the objective is inherently sequential — every accepted move changes the accumulators the next test reads), and the setup passes fill disjoint slots in parallel with a fixed floating-point order. Output is therefore byte-identical across the three mesh backends and thread counts — but there is deliberately **no numpy twin**, like [`subdivide`](/subdivide)/[`agglomerate`](/agglomerate)/[`decimate_volume`](/decimate_volume)/[`conservative_interpolate`](/conservative_interpolate): a single near-tie decided differently by an independent implementation diverges into a genuinely different clustering, not a last-ulp difference. Calling `remesh` without the compiled `meshioplusplus._core` extension raises `NotImplementedError` naming the reason, for any input.

## Other language surfaces

- **C API** — a flat `mio_remesh(mesh, num_clusters, subdivide, subsample_ratio, max_subdivide, max_iterations, max_repair_passes, metric, gradation, preserve_boundary, &num_clusters_out, &num_iterations, &subdivide_applied, &num_isolated_clusters, &num_non_manifold_vertices)` → a plain `mio_mesh*` (the counter out-params are all nullable, `mio_smooth`'s shape — no opaque result handle, since there are no maps to hand back), plus the `mio_refine_ex`/`mio_refine_opts` growth-path pattern applied a second time: `mio_remesh_ex(mesh, const mio_remesh_opts* opts, mio_remesh_report* report)`, an append-only `reserved`-tailed options struct (`mio_remesh_opts_init` fills the defaults, including `max_anisotropy`) and a matching `mio_remesh_report` struct — introduced specifically because `mio_remesh`'s flat signature had already grown once and a `max_anisotropy` parameter would have broken every C/Fortran/Julia/R caller a second time. Plain `mio_remesh` is unchanged and delegates to `mio_remesh_ex` with `max_anisotropy` left at its default.
- **Fortran** — `m%remesh(num_clusters, ..., gradation=..., preserve_boundary=..., max_anisotropy=..., num_non_manifold_vertices=...)`, routed internally through `mio_remesh_ex`/`mio_remesh_opts` (a `bind(c)` mirror type with a runtime layout guard, the same discipline `mio_refine_opts`' Fortran mirror already uses).
- **Julia** — `remesh(mesh, num_clusters; subdivide=-1, metric=:isotropic, gradation=0.0, preserve_boundary=true, max_anisotropy=4.0, ...)` → a `(; mesh, num_clusters, num_iterations, subdivide_applied, num_isolated_clusters, num_non_manifold_vertices)` NamedTuple, also routed through `mio_remesh_ex` with a mirrored, layout-guarded options struct.
- **R** — `mio_remesh(mesh, num_clusters, ..., gradation = 0.0, preserve_boundary = TRUE, max_anisotropy = 4.0)` → a named list of the same six fields (counters as `double`, R having no native int64); builds a real `mio_remesh_opts` from the installed C header directly, needing no manual struct mirror.
- **WASM** — `remesh(mesh, numClusters, subdivide, subsampleRatio, maxSubdivide, maxIterations, maxRepairPasses, metric, gradation, preserveBoundary, maxAnisotropy)`, and reachable as a `{op: 'remesh', numClusters: ...}` `convertSurfaceOps`/pipeline step (`Gradation`/`PreserveBoundary`/`MaxAnisotropy` join that step's param list). WASM binds directly against the C++ `remesh()` overload and needed no `_ex` growth path of its own.
- **CLI** — the `remesh` verb in both the Python and the native CLI (see [CLI](/cli)): `--num-clusters` (required), `--subdivide`, `--subsample-ratio`, `--max-subdivide`, `--iterations`, `--repair-passes`, `--metric isotropic|quadric|anisotropic`, `--gradation`, `--no-preserve-boundary`, `--max-anisotropy`.
- **Pipeline** — a `Remesh` step (`{NumClusters, Subdivide, SubsampleRatio, MaxSubdivide, MaxIterations, MaxRepairPasses, Metric, Gradation, PreserveBoundary, MaxAnisotropy}`), dispatched generically off the shared op table — the one step whose output has no correspondence to its input, exactly like a fresh mesh from `Voxelize`.
- **MCP** — the `remesh` tool.

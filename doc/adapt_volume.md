# Size-field tetrahedral remeshing (`adapt_volume`)

`adapt_volume` turns a valid but over-refined tetrahedral mesh, typically the graded output of [`octree_mesh`](octree_mesh.md), into one with fewer and better-shaped tetrahedra. It works towards an isotropic size field with a gradation and a bound on how far the boundary may move, and keeps corners, ridges, materials, regions and boundary tags. It is native, MIT-licensed, has no dependency and is deterministic: an alternative to post-processing with the LGPL [MMG backend](mmg.md)'s `adapt`.

```python
import meshioplusplus as mp

tets = mp.octree_mesh(mp.read("part.stl"), min_size=0.5)
out, report = mp.adapt_volume(tets, hmin=0.5, hgrad=1.5, return_report=True)
print(report["num_tets_before"], report["num_tets_after"], report["max_boundary_deviation"])

# the same as an opt-in octree_mesh post-stage
out = mp.octree_mesh(mp.read("part.stl"), min_size=0.5, post_adapt=True, adapt_hgrad=1.5)
```

```bash
meshioplusplus adapt-volume tets.vtu adapted.vtu --hmin 0.5 --hgrad 1.5 --json
meshioplusplus octree-mesh part.stl adapted.vtu --min-size 0.5 --post-adapt --adapt-hgrad 1.5
```

## Where it sits

[`optimize_volume`](optimize_volume.md) raises a tetrahedral mesh's worst quality with flips and relocation but keeps its vertex set, so it cannot coarsen. [`decimate_volume`](decimate_volume.md) coarsens by quadric shape error with no size field, gradation or Hausdorff bound. [`remesh_volume`](remesh_volume.md) discards the input for a fresh lattice. `adapt_volume` changes the vertex set and the connectivity together, towards a target size, under the guards a remesher needs: it is the volume counterpart of MMG's `mmg3d` isotropic mode, written from the published algorithms.

## The size field

Every vertex gets a target edge length `h`:

- on the boundary, by `boundary_sizing`:
  - `"input"` (the default): the smallest of the mean length of its boundary edges (the input boundary keeps its own resolution), the `size_array` point data where given, and, with `curvature_sizing`, `sqrt(8 hausd R)` for the local radius of curvature `R` (the chord of a circle of radius `R` with sagitta `hausd`);
  - `"geometry"`: `hmax`, bounded by `size_array` and by `sqrt(8 hausd R)`, as MMG sizes a boundary from its geometry, so a flat or gently curved boundary coarsens too and the Hausdorff and normal guards keep its shape;
- inside, `hmax` (and `size_array` where given);
- then graded, `h(v) = min(h(v), h(u) + (hgrad - 1) |uv|)` along every edge, by a Dijkstra sweep from the smallest sizes, and clamped to `[hmin, hmax]`.

`hmin` therefore also raises the boundary sizes: under `"input"` sizing it is the lever that matches `octree_mesh`'s `min_size` when the stuffed boundary has edges shorter than the skin needs. The default keeps the specified rule (the boundary's resolution is the caller's); `"geometry"` is what to use for the fewest tetrahedra (see [measured](#measured)). By default `hmin` is the smallest boundary size (no clamp), `hmax` an eighth of the bounding-box diagonal and `hausd` 1e-3 of it (`hausd_relative` reads `hausd` as a fraction of the diagonal). Sizes vary geometrically along an edge, so its metric length is `|ab| (hb - ha) / (ha hb ln(hb / ha))`, `|ab| / ha` when the sizes agree (Dapogny, Dobrzynski and Frey 2014), and a new vertex takes the interpolated size.

## The cycle

Each iteration runs four passes, each serial in a fixed order:

1. **Collapse** the edges shorter than 4/5 in the metric, shortest first (a queue with version stamps): the half-edge collapse `a -> b` moves `a`'s tetrahedra onto `b`. Both directions are tried, the less constrained end going first.
2. **Split** the edges longer than 4/3, longest first, at the metric midpoint. A boundary edge's new vertex goes on its ridge chain (by arc length) or is projected onto its patch of the input boundary when the guards allow, and stays at the straight midpoint otherwise. The 4/5 and 4/3 thresholds are Botsch and Kobbelt's (2004).
3. **Swap**: the 2-3 flip, edge removal (3-2 and larger rings, up to 7) and multi-face removal of [`optimize_volume`](optimize_volume.md), run on each material class alone, so an interface face and every boundary face stay. Each move strictly raises the least dihedral angle of the tetrahedra it replaces.
4. **Smooth**: inside, `optimize_volume`'s guarded ODT relocation (Alliez et al. 2005), or with `quality_metric="amips"` a descent of the star's conformal AMIPS energy (Hu et al. 2020), then the worst-cell smoothing for tetrahedra under 40 degrees. On the boundary, a surface vertex slides half way to the mean of its boundary neighbours in the tangent plane and back onto its patch, and a ridge vertex half way to the middle of its two chain neighbours along the chain.

The iterations stop when the topological edits (collapses, splits, swaps) of one iteration change fewer than 1% of the tetrahedra, at `max_iterations` (10), or once the least dihedral angle reaches `target_min_dihedral`. A polish of swaps and worst-cell smoothing (at most three rounds) ends the run.

## The guards

Every edit runs its preconditions, simulates the edit, checks the postconditions on the simulated star and only then commits (Jiang et al. 2022's precondition / edit / postcondition structure; a refused edit changes nothing). A collapse of `a` onto `b` is refused, and counted by the first guard that refuses it, when:

| Guard | Refused when | Report counter |
|---|---|---|
| Class | `a` is a corner or frozen; a surface vertex would leave its patch (the edge is not an edge of two boundary faces of `a`'s patch); a ridge vertex would leave its chain | `num_rejected_by_class` |
| Link | the link condition fails: the vertices, edges or triangles common to the links of `a` and `b` are not exactly those of the link of `ab` (a common triangle would duplicate a tetrahedron); on the boundary, the vertices sharing a boundary face with both ends are not exactly the far corners of the faces through `ab`, or two boundary faces would merge | `num_rejected_by_link` |
| Inversion | a tetrahedron is not strictly positive by the exact `orient3d` predicate after the move, or flatter than `1e-14` of the diagonal cubed | `num_rejected_by_inversion` |
| Length | an edge the collapse creates is longer than 4/3 (which also keeps collapses and splits from undoing each other) | `num_rejected_by_length` |
| Quality | the least quality of the new star is below `min(old star's least, quality_min)` | `num_rejected_by_quality` |
| Normal | a boundary face's normal would turn by more than `max_normal_deviation` (30 degrees) | `num_rejected_by_normal` |
| Hausdorff | a sample of a new boundary face lies farther than `0.9 hausd` from its patch of the input boundary, or a point of the input boundary the old faces covered is no longer within `max(0.9 hausd, its old distance)` of the new faces and the unchanged ones around them | `num_rejected_by_hausdorff` |

Splits check the children's orientation and quality, and a projected boundary vertex the normals and both halves of the envelope too. Smoothing moves keep every tetrahedron positive and the star's least quality, and a boundary move also the normals and the envelope. Swaps compare the worst new tetrahedron with the worst old one. `quality_min` defaults to 10 degrees for `min_dihedral` and 0.2 for `mean_ratio` and `amips` (the mean ratio `12 (3 V)^(2/3) / sum(l^2)`, 1 for the regular tetrahedron; the conformal AMIPS energy is its reciprocal).

The forward half of the Hausdorff guard samples each new boundary face on a barycentric grid of spacing at most `hausd` (2 to 6 subdivisions) against the input boundary's bounding-volume hierarchy, filtered to the face's own patch; the reverse half indexes the input boundary's vertices, edge midpoints and centroids in a uniform grid. The 10% margin covers the face between its samples.

## Classification

The input's boundary faces are the faces of one tetrahedron and the faces between two tetrahedra of different material class (an interface is an internal boundary). A class is the input block, the cell regions a tetrahedron belongs to and its values of the material arrays (`octree_mesh:material`, `medit:ref`, `gmsh:physical`, `gmsh:geometrical`, `cell_tags`, or `material_arrays`). A boundary edge is a ridge when it does not have exactly two boundary faces, when its faces differ in side-region tags or classes, or (with `preserve_features`) when their normals turn by more than `feature_angle` (30 degrees). Patches are the faces connected across non-ridge edges. Ridges form chains between corners, where the ridge graph branches, ends, turns by more than `corner_angle` (30 degrees) or changes the pair of patches it separates; a closed loop with no corner (a cylinder's rim) is cut into two chains at two corners, as `feature_edges`' graph does. Corners never move; ridge vertices move only along their chain, surface vertices only within their patch, interior vertices freely. `preserve_boundary` freezes every boundary vertex, so the boundary triangles are left bit for bit (an interior vertex may still collapse onto one), and `frozen` pins any vertex.

## Guarantees

- No tetrahedron is inverted or flat (exact `orient3d`), except one that was already flat in the input and that no edit could remove; a tangled input, two tetrahedra on the same side of their shared face, is refused by name (run [`untangle`](untangle.md) first).
- The output is conforming, its outer boundary closed, and every material class keeps its side of every interface, with the interface topology unchanged.
- The least quality never drops below `min(the input's least, quality_min)`, and swaps and smoothing never lower it.
- `max_boundary_deviation`, the symmetric Hausdorff distance between the output's boundary and the input's, measured with [`hausdorff_distance`](hausdorff.md) after the run, is at most `hausd` on the tests' sphere, box, cylinder and nested shells; the guard is sampled, so this is a measured bound, not a proof.
- With `preserve_boundary` the boundary vertices and triangles are bit for bit those of the input and the volume is preserved to rounding.
- The output is byte-identical across the mesh backends, the parallel backends and thread counts: every pass is serial in a fixed order with total-order keys, and only per-element measurements run in parallel. `tools/bench_ops.sh --hash` checks this on the `adapt_volume` row ([benchmarks](benchmarks.md#determinism-check)); on this change the digests agreed on SEQ and OpenMP at 1, 4 and 8 threads, for `adapt_volume`, `optimize_volume` and `octree_mesh_tet` at tiers S and M (TBB was not installed).

## What survives

Float `point_data` is interpolated on splits and keeps the survivor's value on collapses; smoothing moves points and not their data (as in `optimize_volume`). Integer `point_data` keeps the value of the input point each vertex comes from. The material arrays are kept by class; any other `cell_data` (`octree:level`, `octree_mesh:shape`, quality arrays) is dropped with a warning (recompute quality with `attach_quality`). Cell regions are rebuilt from the classes, side regions from the boundary faces' tags (a split face's halves keep its tags), and point regions keep their surviving members; a split vertex joins a point region when both ends of its edge are members and it has the kind of one of them, so the midpoint of an edge between two corners does not join a region of corners. `field_data` and property sets pass through. The output has one `tetra` block per input block.

## Report

`num_tets_before`, `num_tets_after`, `num_points_before`, `num_points_after`; the edits `num_collapsed`, `num_split`, `num_swapped`, `num_smoothed`; the refusals `num_rejected_by_link`, `_inversion`, `_quality`, `_length`, `_class`, `_hausdorff`, `_normal`; `min_dihedral_before/after`, `max_dihedral_before/after`, `mean_quality_before/after`; `max_boundary_deviation`; the `hmin`, `hmax` and `hausd` used; `volume_before`, `volume_after` and `region_volume_change` per class; the classification counts `num_corners`, `num_ridge_vertices`, `num_chains`, `num_patches`; `iterations`; `timings` per stage (`setup`, `collapse`, `split`, `swap`, `smooth`, `polish`, `output`); `warnings`.

## Measured

`tools/adapt_volume_bench.py` meshes each skin with `octree_mesh` at `min_size` (a fraction of the bounding-box diagonal), then remeshes the tetrahedra with `adapt_volume` at `hmin` = `min_size`, `hmax` = the octree's `max_size`, `hgrad` 1.3 and the default `hausd` (1e-3 of the diagonal), once with each `boundary_sizing`, and runs MMG 5.8's `mmg3d_O3` on the same tetrahedra with the same four parameters (as a separate program on a Medit file: nothing links it). tetgen and gmsh were not installed and are skipped. The skins are the unit sphere (`icosa_sphere(1, 5)`), the unit cube, a cylinder of radius 0.5 and height 2, `example/data/Bunny.stl`, and the Happy Buddha and XYZ RGB dragon scans from common-3d-test-models (`--fetch`). One run each, serial, on a 4-core x86-64 container, Release build; the raw rows, with the scaled-Jacobian minimum and dihedral histograms, are `benchmark/adapt_volume.csv` (input sizes) and `benchmark/adapt_volume_geometry.csv` (geometry sizes).

"Under 20°" is the share of tetrahedra whose least dihedral angle is below 20 degrees. "Hausdorff to skin" is the symmetric [`hausdorff_distance`](hausdorff.md) from the output's boundary to the *original skin*, so it includes `octree_mesh`'s own distance from the skin (its first row per case); `adapt_volume`'s bound is on its distance from *its input's* boundary, which its report gives as `max_boundary_deviation`. The volume error is against the `octree_mesh` mesh. On the sphere, `octree_mesh` meshes the icosphere's own resolution at both `min_size` values the script tries, so only one is shown.

| Case | `min_size` | Method | Tetrahedra | Points | Dihedral range | Under 20° | Hausdorff to skin | Volume error | Seconds |
|---|---|---|---:|---:|---|---:|---:|---:|---:|
| sphere | 0.02 | octree_mesh | 18,552 | 3,871 | 35.3–125.3 | 0.00% | 0.00323 | 0.000% | 0.6 |
| sphere | 0.02 | adapt_volume, input sizes | 14,975 | 3,276 | 38.1–135.1 | 0.00% | 0.00338 | 0.022% | 5.0 |
| sphere | 0.02 | adapt_volume, geometry sizes | 7,215 | 1,861 | 34.9–138.4 | 0.00% | 0.00385 | 0.096% | 3.9 |
| sphere | 0.02 | mmg3d | 12,094 | 2,665 | 25.3–136.5 | 0.00% | 0.00325 | 0.003% | 1.2 |
| box | 0.03 | octree_mesh | 35,702 | 8,359 | 14.9–149.8 | 0.89% | 2.78e-16 | 0.000% | 0.7 |
| box | 0.03 | adapt_volume, input sizes | 39,800 | 8,684 | 36.4–138.2 | 0.00% | 2.78e-16 | 0.000% | 12.8 |
| box | 0.03 | adapt_volume, geometry sizes | 1,283 | 328 | 30.6–138.9 | 0.00% | 2.48e-16 | 0.000% | 1.8 |
| box | 0.03 | mmg3d | 1,399 | 343 | 25.6–136.6 | 0.00% | 1e-15 | 0.000% | 1.2 |
| box | 0.015 | octree_mesh | 75,995 | 18,219 | 14.9–149.8 | 1.41% | 2.78e-16 | 0.000% | 1.8 |
| box | 0.015 | adapt_volume, input sizes | 96,749 | 21,056 | 33.6–142.7 | 0.00% | 2.48e-16 | 0.000% | 34.3 |
| box | 0.015 | adapt_volume, geometry sizes | 1,291 | 327 | 26.3–136.9 | 0.00% | 2e-16 | 0.000% | 2.9 |
| box | 0.015 | mmg3d | 1,350 | 334 | 20.0–148.4 | 0.07% | 1.01e-15 | 0.000% | 1.6 |
| cylinder | 0.02 | octree_mesh | 63,758 | 13,797 | 11.6–160.3 | 1.66% | 0.00297 | 0.000% | 1.4 |
| cylinder | 0.02 | adapt_volume, input sizes | 54,631 | 11,453 | 35.4–139.5 | 0.00% | 0.0029 | 0.099% | 14.8 |
| cylinder | 0.02 | adapt_volume, geometry sizes | 9,887 | 2,347 | 27.1–142.3 | 0.00% | 0.00313 | 0.202% | 6.2 |
| cylinder | 0.02 | mmg3d | 16,943 | 3,831 | 16.0–149.0 | 0.04% | 0.0022 | 0.030% | 2.0 |
| cylinder | 0.01 | octree_mesh | 204,598 | 42,829 | 7.9–167.1 | 2.87% | 0.00297 | 0.000% | 4.5 |
| cylinder | 0.01 | adapt_volume, input sizes | 164,220 | 33,340 | 31.3–143.9 | 0.00% | 0.00294 | 0.082% | 52.6 |
| cylinder | 0.01 | adapt_volume, geometry sizes | 10,736 | 2,534 | 30.1–136.2 | 0.00% | 0.00336 | 0.204% | 11.1 |
| cylinder | 0.01 | mmg3d | 23,218 | 5,050 | 19.8–145.6 | 0.00% | 0.00215 | 0.033% | 4.5 |
| bunny | 0.012 | octree_mesh | 242,914 | 49,423 | 0.3–179.5 | 2.05% | 0.584 | 0.000% | 7.8 |
| bunny | 0.012 | adapt_volume, input sizes | 165,202 | 33,408 | 2.4–169.5 | 0.04% | 0.587 | 0.046% | 54.4 |
| bunny | 0.012 | adapt_volume, geometry sizes | 45,849 | 10,384 | 4.0–169.8 | 0.19% | 0.587 | 0.223% | 28.5 |
| bunny | 0.012 | mmg3d | 75,251 | 16,354 | 1.1–175.8 | 0.10% | 0.418 | 0.046% | 6.3 |
| bunny | 0.008 | octree_mesh | 629,287 | 125,848 | 0.2–179.6 | 2.07% | 0.313 | 0.000% | 19.4 |
| bunny | 0.008 | adapt_volume, input sizes | 405,123 | 80,444 | 0.3–179.3 | 0.01% | 0.313 | 0.023% | 133.9 |
| bunny | 0.008 | adapt_volume, geometry sizes | 62,203 | 13,704 | 3.8–167.0 | 0.13% | 0.313 | 0.225% | 48.1 |
| bunny | 0.008 | mmg3d | 112,033 | 23,522 | 3.0–170.9 | 0.07% | 0.277 | 0.092% | 13.7 |
| buddha | 0.008 | octree_mesh | 525,773 | 108,850 | 0.1–179.9 | 2.87% | 0.00209 | 0.000% | 29.1 |
| buddha | 0.008 | adapt_volume, input sizes | 376,153 | 78,800 | 1.1–178.5 | 0.65% | 0.00216 | 0.039% | 146.4 |
| buddha | 0.008 | adapt_volume, geometry sizes | 96,771 | 24,293 | 0.9–178.5 | 3.98% | 0.00209 | 0.210% | 64.8 |
| buddha | 0.008 | mmg3d | 182,096 | 40,891 | 0.1–179.7 | 2.09% | 0.00209 | 0.060% | 22.3 |
| dragon | 0.008 | octree_mesh | 195,220 | 42,010 | 0.1–179.8 | 9.26% | 3.07 | 0.000% | 18.3 |
| dragon | 0.008 | adapt_volume, input sizes | 163,314 | 37,757 | 2.1–177.2 | 4.04% | 3.07 | 0.141% | 76.2 |
| dragon | 0.008 | adapt_volume, geometry sizes | 59,796 | 16,654 | 0.7–178.7 | 15.75% | 3.07 | 0.274% | 31.6 |
| dragon | 0.008 | mmg3d | 82,985 | 21,387 | 0.5–179.2 | 9.01% | 3.07 | 0.007% | 19.3 |

![The octree_mesh sphere after adapt_volume with geometry sizes, cut open, each tetrahedron coloured by its least dihedral angle](/images/adapt_volume_sphere.png)

**Against the targets.** The goal was about 1.3 times MMG's tetrahedron count with a least dihedral angle close to MMG's.

- With `boundary_sizing="geometry"` the count is 0.46 to 0.96 times MMG's on every case. The least dihedral angle is above MMG's on every case (on the scans barely: 0.9 against 0.1 degrees on the Buddha, 0.7 against 0.5 on the dragon). On the two scans, though, more tetrahedra end below 20 degrees than with MMG: 4.0% against 2.1% (Buddha) and 15.8% against 9.0% (dragon). Collapses there are accepted down to the quality floor, `min(the input's least, quality_min)`, which a scan's sub-degree octree slivers put at the input's least angle. The volume changes by 0.1–0.27%, against MMG's 0.003–0.09%: the boundary moves within `hausd`, and MMG's moves less.
- With the default `"input"` sizing, the boundary keeps the octree's resolution. The count is then 1.2 times MMG's on the sphere, about 2 times on the scans, 2.2 to 7 times on the bunny and cylinder, and 28 to 72 times on the box, whose flat faces MMG coarsens to a few hundred triangles. On the analytic cases its least angle is the highest of the three (31–38 degrees), and its volume error is below the geometry runs'.
- `max_boundary_deviation` stayed within `hausd` on every run: on the sphere 0.00142 (input) and 0.00273 (geometry) against 0.00346, on the cylinder 0.00220 and 0.00225 against 0.00245, and 0 on the box.

**Runtime.** `adapt_volume` takes 1.5 to 4.5 times MMG's time with geometry sizes and 4 to 21 times with input sizes; most of it is the envelope and link checks of refused collapses. Every pass is serial, for determinism; `mmg3d` is serial too.

**Gradation.** At `hgrad` 1.3 the size grows by at most 30% of the distance, slower than `octree_mesh`'s 2:1 grading, so under input sizing the interior gains tetrahedra; a looser gradation removes them. On the unit sphere at `min_size=0.04` (13,560 tetrahedra from `octree_mesh`, least angle 27.9 degrees):

| `hgrad` | Input sizes: tetrahedra | Least angle | Geometry sizes: tetrahedra | Least angle |
|---:|---:|---:|---:|---:|
| 1.3 | 15,449 | 34.3 | 6,960 | 27.0 |
| 1.6 | 8,111 | 36.5 | 5,201 | 27.4 |
| 2.0 | 6,451 | 34.7 | 3,893 | 27.8 |
| 3.0 | 4,696 | 30.8 | 3,532 | 27.2 |

`max_boundary_deviation` was 0.0026–0.0031 against `hausd` 0.00346 in every row.

The whole benchmark at `hgrad` 2 with geometry sizes, against `mmg3d` at `hgrad` 2 (`benchmark/adapt_volume_geometry_hgrad2.csv`; its seconds were taken while a build ran, so they are left out), gives 0.44 to 0.96 times MMG's tetrahedra, from 4,304 on the sphere (MMG 5,590) to 46,407 on the Buddha (MMG 87,948). The least dihedral angle is above MMG's on every case, and on the scans fewer tetrahedra fall under 10 degrees than with MMG (0.74% against 1.16% on the Buddha, 2.88% against 4.14% on the dragon). The volume still changes more than with MMG (0.10–0.28% against 0.003–0.11%).

**Limits.**

- The Hausdorff guard samples faces, so the bound is measured on these cases, not proven.
- About 2% of edges still churn between collapse and split per iteration once smoothing moves vertices. The stop rule counts only topological edits, so a run usually takes all `max_iterations`.
- R, Octave and WASM have no `post_adapt`; they call `adapt_volume` on `octree_mesh`'s output instead.
- No pass runs in parallel.

## Everywhere else

| Surface | Call |
| --- | --- |
| C++ | `AdaptVolumeResult r = adapt_volume(mesh, opts);` with `AdaptVolumeOptions` (`mHmin`, `mHmax`, `mHgrad`, `mHausd`, `mHausdRelative`, `mSizeArray`, `mCurvatureSizing`, `mFeatureAngle`, `mCornerAngle`, `mMaxNormalDeviation`, `mPreserveBoundary`, `mPreserveFeatures`, `mQualityMetric`, `mQualityMin`, `mMaxIterations`, `mTargetMinDihedral`, `mFrozen`, `mMaterialArrays`, `mEnableCollapse/Split/Swap/Smooth`, `mBoundarySizing`); `OctreeMeshOptions::mPostAdapt` and `mAdapt` |
| Python | `adapt_volume(mesh, *, hmin=0.0, hmax=0.0, hgrad=1.3, hausd=0.0, hausd_relative=False, size_array=None, curvature_sizing=False, feature_angle=30.0, corner_angle=30.0, max_normal_deviation=30.0, preserve_boundary=False, preserve_features=True, quality_metric="min_dihedral", quality_min=-1.0, max_iterations=10, target_min_dihedral=0.0, frozen=None, material_arrays=None, collapse=True, split=True, swap=True, smooth=True, boundary_sizing="input", return_report=False)`; `octree_mesh(..., post_adapt=False, adapt_hmin=, adapt_hmax=, adapt_hgrad=, adapt_hausd=, adapt_max_iterations=, adapt_preserve_boundary=, adapt_boundary_sizing=, adapt_options=None)` |
| Both CLIs | `meshioplusplus adapt-volume IN OUT [--hmin H] [--hmax H] [--hgrad G] [--hausd D] [--preserve-boundary] [--max-iterations N] [--boundary-sizing input|geometry] ... [--json]`; `octree-mesh ... --post-adapt [--adapt-*]` ([CLI](cli.md#meshioplusplus-adapt-volume)) |
| Pipeline | `{"Op": "AdaptVolume", "Hmin": 0.5, "Hgrad": 1.5}`; `OctreeMesh` takes `PostAdapt` and `Adapt*` ([pipelines](pipeline.md)) |
| C | `mio_adapt_volume_opts_init(&o); o.hgrad = 1.5; mio_mesh* out = mio_adapt_volume(mesh, &o, &report);`; `mio_octree_mesh_opts.post_adapt` and `.adapt`, `mio_octree_mesh_adapt_report` ([C API](c_api.md)) |
| Fortran | `out = m%adapt_volume(hgrad=1.5d0, report=rep, stat=ierr)` ([Fortran](fortran.md)) |
| Julia | `r = adapt_volume(mesh; hgrad=1.5)`; `octree_mesh(skin; post_adapt=true, adapt_hgrad=1.5)` ([Julia](julia.md)) |
| R | `r <- mio_adapt_volume(mesh, hgrad = 1.5)` ([R](r.md)) |
| Octave / MATLAB | `r = mio.adapt_volume(mesh, 'hgrad', 1.5)` ([Octave](octave.md)) |
| WASM | `m.adaptVolume(mesh, {hgrad: 1.5})` ([WASM](wasm.md)) |
| MCP | tool `adapt_volume(input_path, output_path, hmin, hmax, hgrad, hausd, ...)`; `octree_mesh` takes `post_adapt` and `adapt_*` ([MCP](mcp.md)) |

R, Octave and WASM have no `post_adapt` on `octree_mesh`: call `adapt_volume` on its mesh, which is the same operation. Fortran, Julia, R, Octave and WASM do not take `frozen` or the material arrays; Fortran, R, Octave and WASM not the curvature, corner and normal options either, which keep their defaults.

C++-core only, with no numpy fallback, like `optimize_volume` and `remesh_volume`: every edit is a discrete branch on exact predicates and near-tie comparisons, so a second implementation could diverge into a different mesh. `MESHIOPLUSPLUS_ADAPT_VOLUME_CHECK=1` runs the invariant checker after every pass and throws on a violation; the fuzz target `meshioplusplus_fuzz_adapt_volume` ([fuzzing](fuzzing.md)) feeds it ulp-jittered, flat, sliver and offset tetrahedra.

## Attribution

Implemented from published descriptions only: Dapogny, Dobrzynski and Frey, "Three-dimensional adaptive domain remeshing, implicit domain meshing, and applications to free and moving boundary problems", *J. Comput. Phys.* 262 (2014) — the metric edge length, the split / collapse / swap / move cycle and Hausdorff-bounded surface changes; Klingner and Shewchuk, "Aggressive Tetrahedral Mesh Improvement", *16th International Meshing Roundtable* (2007), and Klingner's thesis, UCB/EECS-2008-145 — the schedule and quality measures; Freitag and Ollivier-Gooch, "Tetrahedral mesh improvement using swapping and smoothing", *IJNME* 40(21) (1997); Shewchuk, "Two Discrete Optimization Algorithms for the Topological Improvement of Tetrahedral Meshes" (2002) — edge and multi-face removal, through `optimize_volume`; Hu et al., "Fast Tetrahedral Meshing in the Wild", *ACM TOG* 39(4) (2020) — conformal AMIPS and envelope-checked edits; Jiang et al., "Declarative Specification for Unstructured Mesh Editing Algorithms", *ACM TOG* 41(6) (2022) — precondition, edit, postcondition and rollback; Botsch and Kobbelt, "A Remeshing Approach to Multiresolution Modeling" (2004) — the 4/5 and 4/3 thresholds; Alliez et al., "Variational Tetrahedral Meshing", SIGGRAPH 2005 — ODT. The exact `orient3d` is Shewchuk's public-domain predicate, already in the core. No code of MMG, CGAL, TetGen, Gmsh/HXT, Netgen, Mesquite, TetWild, fTetWild or Stellar is read or vendored.

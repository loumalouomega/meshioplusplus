# Hausdorff distance

`hausdorff_distance(a, b)` measures how far apart two surfaces are at their worst: the largest distance from a point of one to the other. It is the number a remeshing, a decimation, a format round trip or a surrogate's reconstructed geometry is checked against — "no point moved further than this".

```python
import meshioplusplus as mp

before = mp.read("part.stl")
after = mp.decimate(before, ratio=0.25)

r = mp.hausdorff_distance(before, after)
r["distance"]                  # max(a_to_b, b_to_a)
r["a_to_b"], r["b_to_a"]       # the one-sided maxima
r["mean_a_to_b"], r["rms_a_to_b"]
r["worst_point_a"]             # the sample of A farthest from B

r = mp.hausdorff_distance(before, after, face_samples=4)   # a tighter estimate
```

## What is measured

Each surface is sampled, and every sample's unsigned distance to the *other* surface is found with the [signed-distance](/sdf) kernel — the same bucket-grid nearest-triangle search, so the two can never disagree. The one-sided distances are the largest sample distance each way, and the Hausdorff distance is the larger of the two. The mean and root-mean-square of each side's sample distances are reported alongside, as are the sample counts and the worst sample of each side.

## Sampling, and why vertices alone give a lower bound

The farthest point of a surface from another lies somewhere on the first surface, not necessarily at a vertex. With `face_samples=0` (the default) only the vertices of each surface are sampled: exact when the farthest point is a vertex — a surface against a refinement of itself, a rigid motion — and a lower bound otherwise. `face_samples=s` also samples the centroid of each of the `s * s` sub-triangles every triangle splits into, which tightens the bound as `s` grows; a flat square against a tent of height 1 over it is `0` at the vertices, `0.47` at `s = 1`, `0.59` at `s = 4` and `0.65` at `s = 8`, approaching the true `0.71` at the centre. The mean and RMS are over the samples, not area-weighted.

## Inputs

Each mesh is a surface — triangles, quads and polygons — or a volume mesh, whose linear [skin](/extract_skin) is compared. `region_a` / `region_b` restrict a surface to a named cell region (a region on a volume mesh is refused: extract its surface first). A mesh with no surface triangles is refused rather than reported as distance zero. The result is deterministic: samples are generated and reduced in a fixed order, and the worst sample is the first to attain the maximum.

## CLI

```sh
meshioplusplus hausdorff A B [--face-samples S] [--region-a NAME] [--region-b NAME] \
    [--input-format-a FMT] [--input-format-b FMT] [--max D] [--json]
```

`--max D` makes it a gate: the exit status is 1 when the distance exceeds `D`, so a CI job can assert on a geometry change the way it asserts on a test.

## Other languages

```c
mio_hausdorff_opts opts;
mio_hausdorff_opts_init(&opts);
opts.face_samples = 4;
mio_hausdorff_report r;
mio_hausdorff_distance(a, b, &opts, &r);   /* r.distance, r.a_to_b, ... */
```

```fortran
d = a%hausdorff_distance(b, face_samples=4_int64, a_to_b=ab)
```

```julia
r = hausdorff_distance(a, b; face_samples=4)
```

```r
r <- mio_hausdorff_distance(a, b, face_samples = 4)
```

```js
const r = m.hausdorffDistance(a, b, 4);
```

Two meshes do not fit a single-mesh [pipeline](/pipeline) step, so there is none; the MCP tool is `hausdorff` (with an optional `max_distance` that adds `passed`).

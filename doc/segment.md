# Segment

`segment(mesh)` splits a surface into smooth patches: regions of cells bounded by the surface's feature edges, each a candidate for one face of a CAD model. It is the first step of turning a triangle skin (an STL, a scan, the output of [`extract_skin`](/extract_skin) or [`isosurface`](/isosurface)) into a solid a CAD tool can edit — [`fit_primitives`](/fit_primitives) measures each patch and the [`step`](/formats/step) writer hands the skin to CAD; see [skin to CAD](/cad_fit). A segmented surface is useful on its own too: as boundary-condition patches, or as input for [`decimate`](/decimate) and [`remesh`](/remesh).

```python
import meshioplusplus as mp

skin = mp.read("bracket.stl")

seg = mp.segment(skin)                                   # bounded by 30-degree creases
seg = mp.segment(skin, feature_angle=20)                  # stricter creases
seg = mp.segment(skin, curvature_tolerance=0.5)           # split a fillet from its planes
seg = mp.segment(skin, normal_tolerance=45)               # bound each patch's normal spread
seg = mp.segment(skin, min_patch_cells=10)                # absorb slivers
seg = mp.segment(skin, method="primitive")                # one patch per plane, cylinder, fillet ...
seg, report = mp.segment(skin, return_report=True)
report["num_patches"], report["num_merged"], report["num_feature_edges"]

seg.cell_data["segment:patch"]                           # the patch of every cell
seg.field_data["segment:adjacency"]                      # rows (a, b, shared edges)
[r.name for r in seg.regions]                            # patch:0, patch:1, ...
```

## How patches grow

Two cells are neighbours when they share an edge that exactly those two cells use; an edge used once (open), three or more times (non-manifold) or by a pair whose dihedral angle exceeds `feature_angle` is never crossed. A patch grows from its seed — the lowest-numbered cell not yet in a patch — breadth first, across every neighbour edge that passes the tests below. The dihedral angle is measured after reorienting a pair that walks its edge the same way, so inconsistent winding does not split a patch (it is the crease test of [`feature_edges`](/feature_edges)).

| Option | Default | A neighbour joins when |
|---|---|---|
| `feature_angle` | 30 | the dihedral angle across the shared edge is at most this many degrees |
| `normal_tolerance` | 0 (off) | its normal, orientation-corrected, lies within this many degrees of the seed's |
| `curvature_tolerance` | 0 (off) | its curvature estimate differs by at most this much (1/length) from the mean estimate of the cells the patch holds so far |

The curvature estimate of a cell is the largest, over its smooth paired edges, of the edge's dihedral angle in radians divided by the distance between the two cells' centroids: about the larger principal curvature, `1/r` on a cylinder or sphere of radius `r` and `0` on a plane. A fillet tangent to its planes has no crease to stop a patch, but its curvature differs from theirs: `curvature_tolerance` is the test that separates them. `normal_tolerance` instead bounds how far a patch may turn, which splits a smooth closed surface (a cylinder, a sphere) into pieces.

After growing, every patch with fewer than `min_patch_cells` cells is merged into the neighbouring patch it shares the most edges with (ties to the lowest patch id), in one pass in ascending patch order; a patch with no neighbour stays. `region="name"` restricts the segmentation to the cells of a named cell region.

## Primitive method

A triangle skin exported from CAD meets a fillet's tangent planes at smooth edges, so angles alone give one patch for both, and no primitive fits that patch; a finely tessellated cylinder is smooth everywhere. `method="primitive"` (C++ `SegmentMethod::Primitive`) ends a patch where its surface ends instead, the reverse-engineering approach of [stlToSolid](https://github.com/Crypto69/stlToSolid) and [stl2curves](https://github.com/Duvdivan/stl2curves). With `t` the `tolerance` (default `1e-3`) times the bounding-box diagonal:

1. Each **smooth area** (cells connected across edges below `feature_angle`) that one primitive explains — [`fit_primitives`](/fit_primitives)' fit with every point within `t` — is one patch: a whole bore, a whole flat face.
2. Any other is **grown** from seeds, largest cell first. From each seed a plane grows across smooth edges while every corner of the next cell lies within 2% of `t` of it (refitted as it doubles), and a sphere, cylinder, cone or torus fitted within `t` to the seed's neighbourhood grows the same way; the curved region wins when it covers more than 1.5 times the plane's area. A plane must hold two cells and 2% of its smooth area, a curved region 1% of it and three points per parameter of its primitive (12 for a sphere, 21 for a torus), and a curved fit's radius must not exceed the diagonal. A region is then trimmed to the cells its own refit explains within 2% of `t`, and kept when that refit's RMS is too.
3. Cells no primitive claims are grown by angles as above; `normal_tolerance` and `curvature_tolerance` apply to these only. An area the kept regions explain less than half of is not split at all, as a freeform surface fits small spheres and tori anywhere.

The strict 2% bound is what a CAD tessellation meets, its points lying exactly on its surfaces, and what keeps a plane from taking the first strip of a fillet tangent to it; a scan wants a `tolerance` some fifty times its noise. The method is deterministic (serial growth, ties to the lowest cell) and has no numpy twin: without the core the Python function raises `NotImplementedError`. [`fit_cad`](/fit_cad) runs it on a skin that has no `segment:patch` yet.

## Refining the patches

`segment` draws its boundaries once: segmented by angles, a fillet tangent to a wall shares the wall's patch, and a cylinder cut at a smooth edge stays in two. [`refine_patches`](/refine_patches) takes the `segment:patch` labels as they come out and refines them by variational shape approximation — boundaries move to where one surface ends, a primitive is split out of a patch no primitive explains, and neighbours one primitive explains merge — then rewrites `segment:patch`, `segment:adjacency` and the patch regions in the same layout, so [`fit_primitives`](/fit_primitives) and [`fit_cad`](/fit_cad) read the result unchanged. Pass it the `feature_angle` the surface was segmented with.

## Output

The input, unchanged — points, connectivity, block order, data and existing regions are bit-identical — plus:

| Output | Type | Content |
|---|---|---|
| `segment:patch` (cell data) | `int64` | the patch of each cell, numbered from 0 in the order of each patch's lowest cell; `-1` on lines, vertices and cells outside `region` |
| `segment:adjacency` (field data) | `int64`, `(m, 3)` | one row `(a, b, edges)` per pair of patches `a < b` sharing paired edges, ascending |
| `<prefix>:<k>` (cell regions) | | one region per patch; `prefix` defaults to `patch` |

Cell regions already named `<prefix>:<integer>` — an earlier run's — are replaced. The report counts the patches, the patches merged, and the feature, open and non-manifold edges of the surface.

Surface cells only: triangles, quads, polygons and their quadratic variants (which contribute their corners). A volume block is refused; take its skin first with [`extract_skin`](/extract_skin).

## Determinism

Normals and the sorted edge table are built in parallel; the growth itself is a serial walk (FIFO queue, neighbours in ascending cell order), because which patch claims a cell depends on the order the walk reaches it — the serial phase is the price of an answer that does not depend on the thread count. The result is byte-identical on every backend, thread count and in the numpy twin, which the tests pin.

## CLI

```sh
meshioplusplus segment IN OUT [--angle DEG] [--normal-tolerance DEG] \
    [--curvature-tolerance K] [--min-patch-cells N] [--region NAME] \
    [--prefix NAME] [--method angle|primitive] [--tolerance T] [--json] [--quiet]
```

## Other languages

```c
mio_segment_opts opts;
mio_segment_opts_init(&opts);
opts.curvature_tolerance = 0.5;
opts.method = MIO_SEGMENT_PRIMITIVE;                    /* or MIO_SEGMENT_ANGLE */
opts.tolerance = 1e-3;
mio_segment_report report;
mio_mesh* seg = mio_segment(mesh, &opts, &report);
```

```fortran
seg = m%segment(curvature_tolerance=0.5_real64, num_patches=npatch)
seg = m%segment(method='primitive', tolerance=1e-3_real64)
```

```julia
r = segment(mesh; curvature_tolerance=0.5)
r = segment(mesh; method="primitive", tolerance=1e-3)
r.mesh, r.num_patches
```

```r
r <- mio_segment(mesh, curvature_tolerance = 0.5)
r <- mio_segment(mesh, method = "primitive", tolerance = 1e-3)
```

```js
const r = m.segment(mesh, 30, 0, 0.5);
r.mesh.cell_data['segment:patch'];
const p = m.segment(mesh, 30, 0, 0, 1, '', 'patch', 'primitive', 1e-3);
```

It is also a settings-pipeline step, `{"Op": "Segment", "FeatureAngle": 30, "NormalTolerance": 0, "CurvatureTolerance": 0, "MinPatchCells": 1, "Region": "", "Prefix": "patch", "Method": "angle", "Tolerance": 0.001}` — see [pipelines](/pipeline) — and the MCP tool `segment`.

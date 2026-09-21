# Normals

`compute_normals(mesh, point_normals=…, cell_normals=…, split_angle=…)` attaches the point and cell normals of a surface, and can split vertices at creases so that every point carries exactly one normal — what a renderer, a shading model or a [glTF export](/formats/gltf) needs, and what a point-cloud learner reads from a `.pcd` or `.xyz` file.

```python
import meshioplusplus as mp

mesh = mp.read("part.stl")

out = mp.compute_normals(mesh)                                # point_data["normals"], one smooth normal per point
out = mp.compute_normals(mesh, cell_normals=True)             # plus cell_data["normals"]
out = mp.compute_normals(mesh, split_angle=30)                # crease at 30 degrees: points are duplicated
out = mp.compute_normals(mesh, split_angle=30, record_parent_ids=True)

out, report = mp.compute_normals(mesh, split_angle=30, return_report=True)
report["num_added_points"]                # points the split appended
report["quality"]["inconsistent_pairs"]   # nonzero: some normals average faces that disagree about "out"

mp.write("part.xyz", mp.compute_normals(mesh))               # the normals become the xyz normal columns
```

## A vertex normal belongs to a patch, not to a position

At the edge of a cube one position has three normals. Averaging them gives one normal that is right for none of the faces, which is what a smooth-shaded render of a hard-edged part looks like. With `split_angle=None` (the default) `compute_normals` returns that average: one normal per point, the angle- or area-weighted mean of the incident unit face normals.

With a `split_angle` the corners around a vertex are grouped into *smooth fans* and every fan beyond the first gets its own copy of the point, so the mesh carries one normal per point and a crease renders as a crease.

## How the fans are found

Two corners of one vertex belong to the same fan when the triangles that own them are joined along an edge through that vertex. An edge joins its two triangles only when all of these hold:

* exactly two triangles use it;
* they walk it in opposite directions, so they agree about which side is out;
* neither is degenerate;
* the angle between their face normals does not exceed `split_angle`.

So a boundary edge, a non-manifold edge and a pair wound the same way always cut. Triangles fanned from one polygon are always joined, so a non-planar quad or polygon never splits along its own diagonals.

The joins are resolved with a union-find whose root is always the smallest corner index of its set, and the groups are numbered by ascending root, so the partition and the numbering depend only on the mesh and the angle — not on the order edges happen to be visited in. A split angle of `180` keeps every proper manifold pair joined: on a closed, consistently wound surface it is the unsplit result, bit for bit.

## The split layout

The layout is the one [`repair`](/repair) uses for its bowtie copies. The original points keep their indices; the copies are appended after them. `point_data` is gathered by row (a copy carries its source's values), a copy joins its source's point regions, and `record_parent_ids=True` attaches `normals:parent_point`, the input point each output point came from. The numbering of cells is untouched, so cell data, cell regions and side regions ride through verbatim.

A group whose incident faces sum to no direction (only degenerate triangles) is not worth a copy of its point: it folds into the point's first defined group.

## It never reorients

Two triangles that disagree about which side is out cannot both be right, and averaging them gives a wrong normal that looks plausible. `compute_normals` never repairs its input. `report["quality"]` carries the input's edge defect counts, `inconsistent_pairs` first among them; when it is nonzero and no split was asked for, a warning says so. The fix is [`repair(mesh, fix_orientation=True)`](/repair) first. A split always cuts at such an edge.

## Weighting

`weight="angle"` (the default) weights each incident face by the angle it subtends at the vertex — the geometrically correct choice (Baerentzen and Aanaes), and the one the [signed distance](/sdf)'s pseudonormals use. `weight="area"` weights by face area: not exact at a vertex, but free of `acos`, so the numpy twin reproduces it exactly on every platform. The cell normals do not depend on the weight.

## Arrays

| Array | Location | Shape | Written when |
|---|---|---|---|
| `normals` | `point_data` | `(n, 3)`, unit, `float64` | `point_normals=True` (default) |
| `normals` | `cell_data` | `(cells, 3)`, unit, `float64` | `cell_normals=True` |
| `normals:parent_point` | `point_data` | `(n,)`, `int64` | `record_parent_ids=True` |

The name `normals` is the one the [PCD](/formats/pcd) and [XYZ](/formats/xyz) readers and writers already use, so `compute_normals` followed by a write to either emits the normal columns with no further step. An existing array of that name is replaced.

A cell normal is the unit vector area of the cell: the sum of its fan of triangle cross products (Newell's normal), which does not depend on which corner the fan starts from and is exact for a non-planar polygon. Points that no selected triangle touches, and points whose faces sum to nothing, are `NaN` and counted in the report as `num_isolated` and `num_undefined`.

## What it accepts

The input is a surface: triangles, quads and polygons (ragged or rectangular), in 2-D or 3-D space; 2-D points are padded with `z = 0`, so a flat mesh gets `(0, 0, 1)`. Lines and vertices are ignored and get `NaN`. A volume or polyhedron block is refused by name pointing at [`extract_surface`](/extract_surface), and a higher-order surface block pointing at [`convert_cells`](/convert_cells)'s `linearize` — the operation does not silently reduce a volume to its skin.

`region="name"` restricts the operation to one named cell region; cells outside it keep their original connectivity and every point they alone touch stays `NaN`.

## Results across the two engines

The compiled core and the numpy fallback (`_normals.py`) are pinned against each other on cubes, spheres, an open cylinder and a mixed polygon/triangle mesh, for both weights and for no split, `0`, `30` and `180` degree splits. The point normals, cell normals, split partition, output points, connectivity and report are byte-identical — the fallback calls `math.acos` per corner rather than numpy's vectorised `arccos`, which may differ from libm in the last bit, and the core is built with `-ffp-contract=off`.

## CLI

```sh
meshioplusplus normals IN OUT \
    [--cell] [--no-point] [--weight angle|area] \
    [--split-angle DEG] [--record-parent-ids] \
    [--region NAME] [--quiet]
```

Both CLIs produce byte-identical files. See the [CLI reference](/cli).

## Other languages

```c
mio_normals_opts opts;
mio_normals_opts_init(&opts);
opts.split = 1;              /* split_angle defaults to 30 degrees */
mio_normals_report report;
mio_mesh* out = mio_compute_normals(mesh, &opts, &report);
report.num_added_points;
```

```fortran
type(mio_mesh) :: out
integer(int64) :: nadd
out = m%normals(split_angle=30.0_real64, num_added_points=nadd)
```

```julia
r = compute_normals(mesh; split_angle=30, record_parent_ids=true)
r.mesh, r.num_added_points, r.quality.inconsistent_pairs
```

```r
r <- mio_compute_normals(mesh, split_angle = 30)
r$mesh; r$num_added_points
```

```js
const r = await m.computeNormals(mesh, true, true, 'angle', 30, true, '');
r.mesh.point_data['normals'];
r.numAddedPoints;
```

In WASM a negative `splitAngle` (the default) means no split. The operation is also a settings-pipeline step: `{"Op": "Normals", "PointNormals": true, "CellNormals": false, "Weight": "angle", "SplitAngle": 30, "RecordParentIds": false, "Region": ""}`, where an absent `SplitAngle` means no split — see [pipelines](/pipeline). The MCP server exposes it as `normals`.

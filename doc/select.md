# Select

`select(mesh, name, kind)` writes one named [region](/regions): the points, cells or cell sides that pass a box, a half-space, a ball, a data or quality threshold, a connected component, or a smooth patch grown from a seed point. It turns a boundary-condition patch, a contact surface or a material zone into a line of a script, reproducible on the next remesh, rather than a pick in a GUI; the region is what a solver deck names and what [`find_interface`](/region_adjacency) selects by.

![A plate with a bolt hole: the hole's wall and the top face grown from seed points, and one outer side picked by a box, each a named side region](/images/select.png)

```python
import meshioplusplus as mp

m = mp.read("bracket.vtu")

m = mp.select(m, "inlet", "side", box=(-1, -1, -1, 0.0, 1, 1))          # boundary faces in a box
m = mp.select(m, "top", "side", seed=(0.5, 0.5, 2.0))                   # the smooth face nearest a point
m = mp.select(m, "rim", "point", sphere=((0, 0, 0), 0.1))               # points in a ball
m = mp.select(m, "upper", plane=((0, 0, 1), (0, 0, 1)), mode="centroid")  # cells above z = 1
m = mp.select(m, "steel", where=("material", "==", 2))                  # cells by cell_data
m = mp.select(m, "slivers", quality=("min_dihedral", "<", 10))          # cells by quality
m = mp.select(m, "part1", component=1)                                  # the second loose part
m, report = mp.select(m, "wall", "side", seed=(0, 0, 0), feature_angle=20, return_report=True)
report["num_selected"], report["num_candidates"], report["num_grown"]
```

The result is the input, unchanged — points, connectivity, block order, data and every other region are bit-identical — plus the region. A region of the same name and kind is replaced, whatever its dimension or tag, so rerunning a script leaves one region; `tag=` sets its integer tag (a gmsh physical tag, an Exodus side-set id).

## What can be selected

| `kind` | Candidates | Entries |
|---|---|---|
| `"cell"` (default) | every cell | global, block-major cell indices |
| `"point"` | every point | point indices |
| `"side"` | the facets of the cells of the mesh's highest dimension: faces of the 3-D cells, or edges of the 2-D cells when there are none; with `boundary_only=True` (the default) only facets owned by one such cell | `(cell, local facet)` pairs in the numbering [regions](/regions) gives |

A side region's dimension is its facets' (2 for faces, 1 for edges); a point or cell region's is `-1`. Polygons name their edges and polyhedra their faces in storage order; cells whose type has no facet table (the VTK Lagrange types) contribute no sides.

## Criteria

Every criterion given must pass, and at least one is required. Each is evaluated over all candidates on its own and the results are intersected — a patch grows over the whole surface even when a box also restricts the selection.

| Criterion | Argument | A candidate passes when |
|---|---|---|
| Box, half-space, ball | `box=(lo3, hi3)` or six numbers, `plane=(point3, normal3)`, `sphere=(center3, radius)` | it is inside all the shapes given: `lo <= p <= hi`, `(p - point) . normal >= 0`, `\|p - center\| <= radius`. A point tests itself; a cell or side reduces its nodes by `mode`: all inside (`"all"`, default), any inside (`"any"`), or its centroid, the mean of its distinct nodes (`"centroid"`) |
| Data | `where=(array, op, value)` | the comparison holds on a scalar `point_data` array (points) or `cell_data` array (cells; a side takes its owner cell's value). `op` is one of `<`, `<=`, `>`, `>=`, `==`, `!=`, the vocabulary [`crop`](/crop) and [`refine`](/refine) share; a NaN never matches |
| Quality | `quality=(metric, op, value)` | the comparison holds on a [`compute_quality`](/mesh_quality) metric (`scaled_jacobian`, `min_dihedral`, `aspect_ratio`, ...), computed on the fly. Cells and sides only; a metric that does not apply to a cell type is NaN there, so the cell never matches |
| Component | `component=k` | the cell is in the k-th node-connected component, numbered as [`split`](/split) with `by="component"` numbers its pieces (each component's lowest cell, ascending). A point passes when a cell of the component uses it, a side when its owner is in it |
| Grow | `seed=(x, y, z)` | the facet is in the patch grown from the facet nearest the seed (by centroid; ties to the lowest index) |

The point data and cell data are not converted into each other: a `where` on a point region refuses a `cell_data` name, and on a cell region a `point_data` name, saying which conversion (`data to-cell` or `data to-point`, see [data averaging](/data_average)) to run.

### Growing a patch

The patch grows breadth first from the seed facet across every edge that exactly two facets use and whose dihedral angle is at most `feature_angle` (default 30 degrees) — the edges [`feature_edges`](/feature_edges) would not draw, by the same crease test. With `normal_tolerance > 0` a facet also joins only when its normal lies within that many degrees of the seed's, which splits a smooth closed surface (a sphere, a cylinder) into a cap.

The facets grown over depend on the kind:

- `"cell"`: the 2-D cells of the mesh. A volume mesh has none; select its sides instead.
- `"side"`: the candidate faces (the boundary faces, or every face with `boundary_only=False`). A mesh whose highest dimension is 2 has edges for sides, which have no normal; select its cells instead.
- `"point"`: the 2-D cells, or the boundary faces of the 3-D cells when there are none; the patch's nodes are selected.

Unlike [`segment`](/segment), which partitions the whole surface into patches, `select` grows one patch from a point you name — the one a boundary condition goes on.

## Feeding the interface operations

Cell regions are what [`find_interface`](/region_adjacency), `region_adjacency` and `contact_pairs` select by, and side regions what `split_interface` takes:

```python
m = mp.select(m, "left", plane=((1, 0, 0), (-1, 0, 0)), mode="centroid")
m = mp.select(m, "right", plane=((1, 0, 0), (1, 0, 0)), mode="centroid")
faces, report = mp.find_interface(m, "left", "right", return_report=True)
```

## Determinism and engines

Per-candidate tests run in parallel into disjoint slots; the growth is a serial breadth-first walk with neighbours in ascending order, and the component labels are [`split`](/split)'s own. The region is the same on every mesh backend and thread count.

`select` is C++-core only: the growth walks the core's crease test and the quality criterion its metrics, and a numpy twin of either would be a second implementation to keep byte-identical (the [`remesh_volume`](/remesh_volume) precedent). Without the compiled core the Python function raises `NotImplementedError`.

## CLI

```sh
meshioplusplus select IN OUT --name NAME [--kind cell|point|side] [--tag N] \
    [--box XMIN,YMIN,ZMIN,XMAX,YMAX,ZMAX] [--plane PX,PY,PZ,NX,NY,NZ] [--sphere CX,CY,CZ,R] \
    [--where 'ARRAY OP VALUE'] [--quality 'METRIC OP VALUE'] [--component K] \
    [--seed X,Y,Z] [--feature-angle DEG] [--normal-tolerance DEG] \
    [--mode all|any|centroid] [--all-facets] [--json] [--quiet]
```

Both CLIs take the same flags; `--where` and `--quality` parse `NAME OP VALUE` as `crop --where` does.

## Other languages

```c
mio_select_opts opts;
mio_select_opts_init(&opts);           /* an all-zero struct is NOT the default */
opts.name = "top";
opts.kind = MIO_REGION_SIDE;
opts.use_seed = 1;
opts.seed[0] = 0.5; opts.seed[1] = 0.5; opts.seed[2] = 2.0;
opts.quality_metric = "scaled_jacobian";
opts.quality_compare = MIO_REFINE_GT;
opts.quality_value = 0.2;
mio_select_report report;
mio_mesh* out = mio_select(mesh, &opts, &report);
```

```fortran
out = m%select_region('top', kind='side', seed=[0.5_real64, 0.5_real64, 2.0_real64], &
                      num_selected=n)
out = m%select_region('slivers', quality_metric='min_dihedral', quality_op='<', &
                      quality_value=10.0_real64)
```

```julia
r = select_region(mesh, "top"; kind="side", seed=(0.5, 0.5, 2.0))
r = select_region(mesh, "steel"; where_array="material", where_op="==", where_value=2)
r.mesh, r.num_selected
```

```r
r <- mio_select(mesh, "top", kind = "side", seed = c(0.5, 0.5, 2))
r <- mio_select(mesh, "part", component = 2)
```

```js
const r = m.selectRegion(mesh, 'top', { kind: 'side', seed: [0.5, 0.5, 2] });
const s = m.selectRegion(mesh, 'slivers', { quality: ['min_dihedral', '<', 10] });
```

Fortran's method is `select_region`, because Fortran reads `select` as the start of a construct, and Julia's is too, so it does not clash with `DataFrames.select`. In Fortran, Julia and R `component` is 1-based, like every other index there; in C, Python, WebAssembly, the CLIs and the pipeline it is 0-based.

It is also a settings-pipeline step, `{"Op": "Select", "Name": "top", "Kind": "side", "Seed": [0.5, 0.5, 2.0]}`, with the keys `Name`, `Kind`, `Tag`, `Mode`, `BoundaryOnly`, `Bbox`, `Point` + `Normal`, `Center` + `Radius`, `Where` + `Compare` + `Value`, `Quality` + `QualityCompare` + `QualityValue`, `Component`, `Seed`, `FeatureAngle` and `NormalTolerance` — see [pipelines](/pipeline) — and the MCP tool `select`.

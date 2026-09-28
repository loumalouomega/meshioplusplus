# Feature edges

`feature_edges(mesh, feature_angle=30)` returns the edges of a surface that matter to anyone meshing, loading or inspecting it — the sharp creases, the open boundary, the non-manifold junctions and the pairs of faces that disagree about which side is out — as a mesh of `line` cells. It is the crease test [`decimate`](/decimate), [`decimate_volume`](/decimate_volume) and [`smooth`](/smooth) pin nodes with, exposed for inspection and for picking the edges a boundary condition lives on.

```python
import meshioplusplus as mp

mesh = mp.read("part.stl")

edges = mp.feature_edges(mesh)                          # 30 degrees, every category
edges = mp.feature_edges(mesh, feature_angle=60)        # only the sharper creases
edges = mp.feature_edges(mesh, boundary=False)          # creases only, not the open rim
edges, report = mp.feature_edges(mesh, return_report=True)
report["num_feature"], report["num_boundary"], report["num_non_manifold"], report["num_inconsistent"]

mp.write("part_edges.vtp", edges)                       # view it on top of the surface
```

## Which edges

Every edge of the surface's faces is classified from the faces that use it:

| Category | `feature:kind` | Rule |
|---|---|---|
| non-manifold | 3 | used by three or more faces |
| boundary | 2 | used by one face |
| inconsistent | 4 | used by two faces that walk it the same way (they disagree about "out") |
| feature | 1 | used by two faces whose normals are more than `feature_angle` degrees apart |

An edge in several categories is labelled by the first selected one in the order of the table. For an inconsistent pair the angle is measured after reorienting one face, so a flat but inconsistently wound pair is inconsistent and not sharp. Only the two faces that share an edge are compared: a smooth surface that is merely coarsely tessellated has no feature edges however few triangles it has. A face with no area gives no angle, so its edges are never sharp.

The edges of a quad or polygon are its own ring edges; a non-planar quad never reports its diagonal. Higher-order cells contribute their corners.

## Volume meshes

When the mesh has volume cells, the surface examined is their skin: the faces used by one selected cell (the [global face table](/polyhedra) orients each outward). `region="name"` restricts to a named cell region, of surface cells or of volume cells whose skin is taken — so the skin of a sub-volume includes its interface with the rest of the mesh.

## Output

The input's points, verbatim and not compacted, so a line's node ids are the input's own; one `line` block in ascending `(low, high)` endpoint order; and two cell-data arrays:

| Array | Type | Content |
|---|---|---|
| `feature:kind` | `int64` | 1 feature, 2 boundary, 3 non-manifold, 4 inconsistent |
| `feature:angle` | `float64` | the dihedral angle in degrees; `NaN` where it is undefined |

Point data, field data and point regions ride through; cell and side regions name cells the edge mesh does not have and are dropped with a warning. The report counts every category over the whole surface, whether or not it was selected for output.

## How decimate and smooth use it

Since v16.23.0 the three operations that preserve features pin the endpoints of every sharp or non-manifold edge this test reports (open edges stay their own `preserve_boundary` / `fix_boundary` decision). Before, each carried a per-vertex test that compared *every* pair of faces around a vertex, adjacent or not, which pinned every vertex of a coarsely tessellated sphere and let a decimation do nothing. `smooth` on a surface mesh also used to look only at the boundary polyline, in the XY plane; it now also holds the creases of the surface itself, and measures the polyline's corners in 3-D.

## CLI

```sh
meshioplusplus feature-edges IN OUT [--angle DEG] \
    [--no-feature] [--no-boundary] [--no-non-manifold] [--no-inconsistent] \
    [--region NAME] [--json] [--quiet]
```

## Other languages

```c
mio_feature_edges_opts opts;
mio_feature_edges_opts_init(&opts);
opts.feature_angle = 45.0;
mio_feature_edges_report report;
mio_mesh* edges = mio_feature_edges(mesh, &opts, &report);
```

```fortran
edges = m%feature_edges(feature_angle=45.0_real64, num_feature=nfeat)
```

```julia
r = feature_edges(mesh; feature_angle=45)
r.mesh, r.num_feature
```

```r
r <- mio_feature_edges(mesh, feature_angle = 45)
```

```js
const r = m.featureEdges(mesh, 45);
r.mesh.cell_data['feature:kind'];
```

It is also a settings-pipeline step, `{"Op": "FeatureEdges", "FeatureAngle": 30, "Feature": true, "Boundary": true, "NonManifold": true, "Inconsistent": true, "Region": ""}` — see [pipelines](/pipeline) — and the MCP tool `feature_edges`.

# Boundary layers — `boundary_layer`

`boundary_layer` grows prismatic layers from a surface (roadmap 6.9; [design](./boundary_layer_plan.md)): `num_layers` layers of cells along the surface's normals, wedges from triangles and hexahedra from quads, the first `first_height` thick and each next `growth` times the one before. The surface is any manifold, consistently oriented triangle and quad mesh — a closed one from [`mesh_cad_surface`](./mesh_cad.md) or [`extract_surface`](./extract_surface.md), or an open one.

```python
import meshioplusplus as mp

skin = mp.mesh_cad_surface(mp.read("wing.cad.json"), 0.05)
layers, report = mp.boundary_layer(skin, 1e-3, num_layers=8, growth=1.25, return_report=True)
print(report["num_wedges"], report["min_scaled_jacobian"], report["num_shrunk"])
```

![Layers grown out of the all-quad unit cube (hexahedra) and the capped cylinder (wedges) and into the quad-dominant sphere, cut open, coloured by layer](/images/boundary_layer.png)

## How it works

1. **Normals.** Each point grows along its *most normal* normal, the unit vector `n` maximising the least `n · n_f` over its faces (Aubry and Löhner 2009), found from the angle-weighted mean by steps toward the face that sees it least (100, shrinking as `1 / (k + 2)`). Its offsets are divided by that least cosine (at most four times), so every face's layers keep their thickness at a ridge or a corner: the cube's cap lies exactly on the planes the total thickness off its faces. A point no direction sees all its faces from (a fold past 90 degrees, cosine under 0.05) does not grow.
2. **Collisions.** A ray from each point along its direction against the surface (a bounding-volume hierarchy, Möller–Trumbore): the point's layers take at most `collision_fraction` (0.45) of the distance to what it hits, so layers growing toward each other from both sides of a gap stop before they meet. The share a point keeps is smoothed `smooth_iterations` (5) times over its neighbours, each pass taking the lesser of a point's own share and its mean with its neighbours', so it is never raised and the layers taper rather than step.
3. **Validity.** Every cell is judged by its least corner Jacobian over its corner edges' lengths (1 for a right prism). The points of a cell that folds have their layers halved, at most `max_shrink` (6) times, and then collapsed: a collapsed point stays on the surface, and the cells around it become pyramids (one collapsed corner of a wedge) and tetrahedra (two), its quads split into triangles in every layer so the faces the cells share stay conforming. With `collapse="layers"` a point instead drops its top layer, so it grows fewer layers at their full thickness, and the cells above its last layer end in pyramids and tetrahedra on it (R. V. Garimella and M. S. Shephard, Int. J. Numer. Meth. Engng. 49, 2000); a point with none left is collapsed. No cell is written folded.
4. **Fronts.** Once every cell is valid, each cell is tested against the surface's faces and the last layer's (`front_check`, on by default), the advancing-layer front test of Pirzadeh (1996) done in volume: a face that crosses one of the cell's faces, by the exact triangle–triangle kernel of `detail/tri_tri` (Shewchuk's `orient3d` and `orient2d` decide every case, in the spirit of Guigue and Devillers 2003), or lies inside the cell, meets it, over a bounding-volume hierarchy. A face whose column shares a surface point with the cell is its neighbour and is not tested. The points of a cell that meets a face, and of a last-layer face it meets, give way as a folding cell's do, and the layers are judged again. A ray misses what lies between its point's neighbours; the front test does not: under a lid of a third of its width hung 0.5 above it, a plate's corners' rays miss the lid and its layers would reach 1.0, through the lid; with the test they stop at 0.25. On the cube and the sphere it finds nothing, so the meshes are the test-off ones byte for byte.
5. **Smoothed normals.** `normal_smooth_iterations` Laplacian passes (Jacobi, so the order does not matter) move each direction toward its neighbours' mean, a move kept only while every face around the point still sees it at more than half the least cosine its most normal normal had; the stretch then follows the new least cosine. Off by default.

## Output

The cells come first (`wedge`, `hexahedron`, `pyramid`, `tetra`, the blocks that are not empty), then the surface (the *wall*) and the last layer (the *cap*) as `triangle` and `quad` faces, split like their columns. The faces of the wall and the cap are oriented with their normals along the growth. The cell data `boundary_layer:layer` holds each cell's layer (1 to `num_layers`), 0 for the wall and `num_layers + 1` for the cap, and the Cell regions `boundary_layer:cells`, `boundary_layer:wall` and `boundary_layer:cap` name them. Every cell data array of the surface's faces (`cad:face` among them) is carried to the cells and faces grown from them, its point data to every copy of a point, a Cell region of the surface to everything grown from its faces and a Point region to every copy of its points; side regions are dropped with a warning.

The volume the cap closes off is filled by [`tetrahedralize`](./tetrahedralize.md), which keeps the cap's triangles, so the tetrahedra conform to the last layer's faces (`tests/cpp/test_tetrahedralize.cpp`, the sphere's layers grown inward), or by [`octree_mesh`](./octree_mesh.md), fitted to the cap rather than conforming to it.

## Options and report

| option | meaning |
|---|---|
| `first_height` | the first layer's thickness (required, positive) |
| `num_layers` | 1 to 1000 (3) |
| `growth` | each layer over the one before, 1 to 10 (1.2) |
| `reverse` | grow against the faces' normals (right-hand rule) rather than along them (false) |
| `collision_fraction` | the share of the distance to the surface a point's ray hits that its layers may take, `(0, 1]` (0.45) |
| `smooth_iterations` | passes smoothing that share over neighbours (5) |
| `max_shrink` | halvings of a point's layers before they collapse (6) |
| `collapse` | how a point whose cells fold or cross gives way: `points`, its layers halved and then collapsed, or `layers`, its top layer dropped (`points`) |
| `normal_smooth_iterations` | Laplacian passes over the directions, 0 to 1000 (0) |
| `front_check` | test every cell against the surface's and the last layer's faces (true) |

The report has `num_cells`, `num_wedges`, `num_hexahedra`, `num_pyramids`, `num_tetrahedra`, `num_shrunk` (points whose layers are thinner than asked), `num_collapsed` (of those, points that do not grow) and `num_not_visible` (of those, points no direction sees all their faces from), `num_fewer_layers` (with `collapse="layers"`, the points that grow some layers but fewer than asked), `num_front_crossings` (the cell and face pairs the front test found, over every round), `min_scaled_jacobian` over the cells, `min_thickness` and `max_thickness` (the total along the face normals, over points), `cap_closed` (every cap edge shared by two cap faces) and `warnings`. A surface with no triangle or quad, another 2-D or 3-D cell type (lines and vertices are skipped), an edge used by more than two faces or faces that disagree on an edge's direction is refused by name.

## Bindings

| surface | entry |
|---|---|
| Python | `meshioplusplus.boundary_layer(mesh, first_height, *, num_layers=3, growth=1.2, reverse=False, collision_fraction=0.45, smooth_iterations=5, max_shrink=6, return_report=False)`; no numpy twin |
| C++ | `boundary_layer(const Mesh&, const BoundaryLayerOptions&)` → `BoundaryLayerResult` (`operations/boundary_layer.hpp`) |
| C | `mio_boundary_layer(mesh, opts, report)` with `mio_boundary_layer_opts` (64 bytes, through `mio_boundary_layer_opts_init`) and `mio_boundary_layer_report` (128 bytes) |
| Fortran | `mio_boundary_layer(mesh, first_height, num_layers, growth, reverse, collision_fraction, smooth_iterations, max_shrink, report, stat, errmsg)` |
| Julia | `boundary_layer(mesh; first_height, num_layers, growth, reverse, collision_fraction, smooth_iterations, max_shrink)` |
| R | `mio_boundary_layer(mesh, first_height, num_layers, growth, reverse, collision_fraction, smooth_iterations, max_shrink)` |
| WASM | `boundaryLayer(mesh, {firstHeight, numLayers, growth, reverse, collisionFraction, smoothIterations, maxShrink})` |
| CLI | `boundary-layer INFILE OUTFILE --first-height H [--num-layers N] [--growth G] [--reverse] [--collision-fraction F] [--smooth N] [--max-shrink N] [--json]`, both CLIs |
| MCP | `boundary_layer(input_path, output_path, first_height, ...)` |
| Pipeline | `BoundaryLayer` with `FirstHeight` (required), `NumLayers`, `Growth`, `Reverse`, `CollisionFraction`, `SmoothIterations`, `MaxShrink`, `Collapse`, `NormalSmoothIterations`, `FrontCheck`; the step reports `NumCells`, `NumShrunk`, `NumCollapsed`, `NumFewerLayers`, `NumFrontCrossings` and `MinScaledJacobian` |

## Limitations

- Layers give way where they would collide (the ray, then the front test), by halving or by dropping layers, not by bending around what they meet.
- The inside test of the front check takes each valid cell as convex, which a prism whose side faces twist is only nearly; a face that touches a cell is counted as meeting it.
- Linear cells only, by decision: a curved surface is linearized first. Layers of curved cells (`wedge15`, `hexahedron20`) need their validity proved, by a Bézier bound on a prism's or a hexahedron's Jacobian that the core does not have ([decisions](/roadmap#non-goals-and-decisions-taken)).

## References

- A. Aubry and R. Löhner, "On the 'most normal' normal", *Communications in Numerical Methods in Engineering* 25 (2009).
- S. Pirzadeh, "Three-dimensional unstructured viscous grids by the advancing-layers method", *AIAA Journal* 34 (1996).
- R. V. Garimella and M. S. Shephard, "Boundary layer mesh generation for viscous flow simulations", *International Journal for Numerical Methods in Engineering* 49 (2000).
- T. Möller and B. Trumbore, "Fast, minimum storage ray-triangle intersection", *Journal of Graphics Tools* 2 (1997).

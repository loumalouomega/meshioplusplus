# Sweeps — `extrude` and `revolve`

`extrude` sweeps a surface into a volume, or a curve into a surface, along a vector (roadmap 3.1.2); `revolve` sweeps it around an axis (3.1.3). The input's top-dimensional cells are the *profile*: `triangle` and `quad` cells sweep into `wedge` and `hexahedron` cells, or, when the mesh has none, `line` cells into `quad` cells. Both run on the layer sweep [`boundary_layer`](./boundary_layer.md) grows its layers with (`detail/layer_sweep`), and their output follows the [parametric primitives](./primitives.md): the swept cells only, with named side regions.

```python
import meshioplusplus as mp

slab = mp.extrude(mp.rectangle(dims=(8, 4)), (0.0, 0.0, 0.5), num_layers=6, ratio=1.3)
ring, report = mp.revolve(mp.annulus(), axis=(1.0, 0.0, 0.0), axis_point=(0.0, 2.0, 0.0), return_report=True)
```

## `extrude`

The points are copied to `num_layers` levels along `direction`, the whole displacement of the last level: level `k` lies at `x + t_k d`, with `t_k = k / L` by default, `t_k = (r^k - 1) / (r^L - 1)` for a geometric `ratio = r` (each layer `r` times as thick as the one below), or the explicit `offsets` (fractions of `direction`, increasing, positive, the last 1; not with a ratio other than 1). Uniform levels are computed as `x + (d k) / L`, the way `box` computes its coordinates, so extruding `rectangle` gives `box`'s points and cells bit for bit, and `bottom`, `top` and the four carried side regions equal its `zmin`, `zmax`, `xmin`, `xmax`, `ymin` and `ymax`.

Each face is turned to point along the direction before it is swept (`num_flipped` counts them), so every cell has a positive volume whatever the faces' winding; a face, or a curve's line, parallel to the direction is refused by its index. A curve's quads face `(b - a) x d` for a line from `a` to `b`: a counter-clockwise circle extruded along +z gives an outward wall.

## `revolve`

The points are turned about the axis through `axis_point` along `axis` by `angle` degrees (the right-hand rule; a negative angle turns the other way) in `num_layers` equal steps. Angles are degrees here, as in [`transform`](./transform.md), where the `annulus` primitive takes radians. A full turn (`|angle| = 360`) closes on itself: its last level is the profile again, so there is neither a `start` nor an `end`. It needs at least 3 steps, and no step may reach 180 degrees. The default axis, `y` through the origin, turns an `(r, z)` profile drawn in the xy plane into the solid it describes.

A point on the axis (within `tolerance` of the profile's bounding-box diagonal) stays where it is, so the cells around it collapse: a triangle with one corner on the axis gives a pyramid per step, with two a tetrahedron, a quad with an edge on the axis a wedge, a quad with one corner on it a pyramid and a wedge (split along its first diagonal, as `boundary_layer` splits a collapsing column), and a line with one end on it a triangle — the poles of a revolved curve. A rectangle touching the axis revolved in `L` steps is a polygonal cylinder whose volume is exactly `L/2 sin(2π/L) R² H`, which the tests check.

Each face is turned to point along its motion. A face whose corners lie on both sides of the axis is refused (it crosses it), as is one that moves within its own plane (an axis normal to the profile's plane), and a profile with faces on both sides of the axis cannot turn more than 180 degrees, since each side would sweep through the other; a smaller turn of it is allowed, with a warning. A line lying on the axis would sweep to nothing and is refused.

## Output

The cells come layer by layer (every cell of layer 1, then layer 2, ...), in blocks `wedge`, `hexahedron`, `pyramid`, `tetra` for a volume and `quad`, `triangle` for a surface, empty blocks left out. The points are the profile's, then each level's own copies, level by level (a point on the axis has none).

Generated `Side` regions, of dimension 2 for a volume and 1 for a surface:

| region | `extrude` | `revolve` |
|---|---|---|
| `bottom` / `start` | the profile, facing away from the sweep | the profile, unless the turn is full |
| `top` / `end` | the last level | the last level, unless the turn is full |
| `side` | what the profile's open border sweeps to: the edges used by one face, or a curve's two ends | the same; nothing for a closed profile or a curve's ends on the axis |

An input region of the same name wins, with a warning.

Carried from the input: point data per copy, cell data per layer (an array whose dtype or width differs between the profile's blocks is dropped, with a warning), field data. A `Cell` region of profile cells names their columns; a `Point` region every copy of its points; a `Side` region of profile facets (the edges of 2-D cells) the faces those edges sweep to; and a `Cell` region of lower-dimensional cells — a gmsh mesh's boundary lines, a curve's end vertices — the faces (or edges) those cells sweep to, as a `Side` region. Every region's dimension grows by one and its tag is kept. The lower-dimensional cells themselves are not swept, and their cell data is dropped with a warning. With `record_layers`, an Int32 `extrude:layer` (`revolve:layer`) cell array holds each cell's layer, 1 to `L`.

Higher-order, polygon, polyhedron and volume cells are refused by name; [`convert_cells`](./convert_cells.md) linearizes or triangulates a profile first.

## Options and report

| option | `extrude` | `revolve` |
|---|---|---|
| `direction` | the whole displacement (`(0, 0, 1)`) | |
| `axis`, `axis_point` | | the axis (`(0, 1, 0)` through the origin) |
| `angle` | | degrees, `0 < |angle| <= 360` (360) |
| `num_layers` | layers (1) | steps (16) |
| `ratio`, `offsets` | the grading (1, none) | |
| `tolerance` | | on the axis, a fraction of the bounding-box diagonal (1e-9) |
| `record_layers` | attach `extrude:layer` (false) | attach `revolve:layer` (false) |
| `max_cells` | refuse more profile cells times layers (20 000 000; non-positive disables) | the same |

The report has `num_cells`, `num_hexahedra`, `num_wedges`, `num_pyramids`, `num_tetrahedra`, `num_quads` and `num_triangles` (a curve's), `num_flipped` (profile cells turned to face the sweep), `min_scaled_jacobian` (the least corner Jacobian over the corner edges' lengths of a volume cell, 1 for a right prism; 0 for a surface) and `warnings`; `revolve` adds `num_on_axis` (profile points on the axis) and `full`.

## Bindings

| surface | entry |
|---|---|
| Python | `meshioplusplus.extrude(mesh, direction=(0, 0, 1), *, num_layers=1, ratio=1.0, offsets=None, record_layers=False, max_cells=20_000_000, return_report=False)` and `meshioplusplus.revolve(mesh, axis=(0, 1, 0), *, axis_point=(0, 0, 0), angle=360.0, num_layers=16, tolerance=1e-9, record_layers=False, max_cells=20_000_000, return_report=False)`; no numpy twin |
| C++ | `extrude(const Mesh&, const ExtrudeOptions&)` and `revolve(const Mesh&, const RevolveOptions&)` → `SweepResult` (`operations/extrude.hpp`) |
| C | `mio_extrude(mesh, opts, report)` with `mio_extrude_opts` (96 bytes; `offsets` borrowed for the call) and `mio_revolve(mesh, opts, report)` with `mio_revolve_opts` (112 bytes), both reporting into `mio_sweep_report` (128 bytes), with `mio_extrude_warnings` and `mio_revolve_warnings` |
| Fortran | `mio_extrude(mesh, direction, num_layers, ratio, offsets, record_layers, max_cells, report, stat, errmsg)` and `mio_revolve(mesh, axis, axis_point, angle, num_layers, tolerance, record_layers, max_cells, report, stat, errmsg)` with `mio_sweep_report` |
| Julia | `extrude(mesh; direction, num_layers, ratio, offsets, record_layers, max_cells)` and `revolve(mesh; axis, axis_point, angle, num_layers, tolerance, record_layers, max_cells)` |
| R | `mio_extrude(mesh, direction, num_layers, ratio, offsets, record_layers, max_cells)` and `mio_revolve(mesh, axis, axis_point, angle, num_layers, tolerance, record_layers, max_cells)` |
| WASM | `extrude(mesh, {direction, numLayers, ratio, offsets, recordLayers, maxCells})` and `revolve(mesh, {axis, axisPoint, angle, numLayers, tolerance, recordLayers, maxCells})` |
| CLI | `extrude INFILE OUTFILE [--direction=x,y,z] [--num-layers N] [--ratio R] [--offsets a,b,..] [--record-layers] [--max-cells N] [--json]` and `revolve INFILE OUTFILE [--axis=x,y,z] [--axis-point=x,y,z] [--angle DEG] [--num-layers N] [--tolerance T] [--record-layers] [--max-cells N] [--json]`, both CLIs |
| MCP | `extrude(input_path, output_path, direction, num_layers, ratio, offsets, record_layers, max_cells)` and `revolve(input_path, output_path, axis, axis_point, angle, num_layers, tolerance, record_layers, max_cells)` |
| Pipeline | `Extrude` with `Direction`, `NumLayers`, `Ratio`, `Offsets`, `RecordLayers`, `MaxCells`, and `Revolve` with `Axis`, `AxisPoint`, `Angle`, `NumLayers`, `Tolerance`, `RecordLayers`, `MaxCells`; they report `NumCells`, `MinScaledJacobian` and `NumFlipped` (`Extrude`) or `NumOnAxis` (`Revolve`) |

## Limitations

- Linear profiles only: a quadratic profile is linearized first, and the swept cells are linear (`wedge15` and `hexahedron20` sweeps would need the curved validity bound [`boundary_layer`](./boundary_layer.md#limitations) also lacks).
- `revolve` places the levels' points on circles, but the cells between them are straight-sided: a coarse turn is a polygonal solid, not a curved one.
- A curve's side cannot be told apart in a profile on both sides of the axis, so only a face profile is refused for a turn of more than 180 degrees through itself.

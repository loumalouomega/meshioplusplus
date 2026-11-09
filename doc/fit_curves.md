# Lines, arcs and circles from sections

`fit_curves(sections)` turns chained sections into a sketch: each chain of a chained line mesh ([`slice(chain=True)`](./slice.md#chained-sections) or `sections`) is split into the runs a line or a circular arc explains within a tolerance, in its section plane's 2-D frame. The sketch is what a drawing of a scanned profile, a DXF file or a CAD profile solid is made of (roadmap 7.6.2).

```python
import meshioplusplus as mp

scan = mp.read("part_scan.stl")
profile = mp.sections(scan, direction=(0, 0, 1), offsets=[2.0])
sketch, report = mp.fit_curves(profile, tolerance=0.01, snap_angle=1.0, return_report=True)
for piece in report["pieces"]:
    print(piece["type"], piece["start"], piece["end"], piece["radius"])
```

## The frame

Each chain is planar. Its plane is `normal` through the chain's points when given, else the least-spread direction of the points of all the chains cut at one `section:offset` (all chains without one). The frame's in-plane axes complete the normal deterministically, and its origin is the normal times the plane's offset, so the sketches of parallel sections share one 2-D frame. The report lists each frame as a row-major 4 × 4 whose columns are the two in-plane axes, the normal and the origin, so `x = F (u, v, 0, 1)`; a chain lying off its plane by more than the tolerance is warned about.

## How a chain is split

1. **Points with no shape are skipped.** A point on the straight segment between its neighbours (where a section cut a triangle's diagonal, the chord's midpoint) is left out of the fits; it still counts in the residuals and its segments in the pieces.
2. **A whole circle.** A closed chain whose points all lie within the tolerance of one circle is one circle.
3. **Split.** Otherwise the run is split recursively (a closed chain first at its first point and the point farthest from it): a run one line explains within the tolerance (a total least-squares line) is a line; failing that, a run one circle explains (Kåsa's algebraic circle, refined by Gauss–Newton on the geometric distance) is an arc when it spans at least `min_arc_points` points and `min_arc_sweep` degrees and no chord of it is longer than four times its median chord (so an arc never bridges a straight side); failing both, the run is split at its point farthest from its chord.
4. **Merge.** Neighbouring pieces whose union one line or one arc explains are merged, left to right, until none is; a closed chain is then started at its last piece's start and fitted again, kept when that needs fewer pieces (the seam).
5. **Snap.** With `snap_angle` > 0, a line within that many degrees of the frame's axes or their 45-degree diagonals is turned onto it about its midpoint.
6. **Corners.** Consecutive pieces' shared end moves to their intersection (line-line, line-arc or arc-arc, the nearest to the data point) when one lies within ten tolerances of it; open ends are their data points projected onto the piece.

Residuals are measured at the chain's points, not between them. A slot sliced through a tessellated prism (the probe) gives two lines and two arcs with radii to 1e-6; a plate with a hole gives four lines and a circle.

## What comes back

The input with the Int64 cell data `sketch:piece` (each segment's piece) and the field data `sketch:type` (1 line, 2 arc, 3 circle), `sketch:params` (ten numbers a piece: start x, y, end x, y, centre x, y, radius, sweep in degrees, residual, points), `sketch:chain`, `sketch:frame` and `sketch:frames` (the frames, sixteen numbers each), so a format or a later pipeline step reads the sketch. The report lists the `pieces` (`type`, `start`, `end`, `centre`, `radius`, `sweep` — signed, counter-clockwise positive, 360 for a circle — `residual`, `num_points`, `chain`, `frame`), the `frames`, `num_lines`, `num_arcs`, `num_circles`, `max_residual`, the `tolerance` applied (0 takes 1e-3 of the bounding-box diagonal) and `warnings`; `meshioplusplus._fit_curves.sketch_json(report)` serialises the sketch.

Core only, no numpy twin. Deterministic: no sampling, every tie to the lower index.

## Everywhere else

| Surface | Call |
| --- | --- |
| Python | `fit_curves(sections, *, tolerance=0.0, snap_angle=0.0, normal=None, min_arc_points=4, min_arc_sweep=5.0, return_report=False)` |
| Both CLIs | `meshioplusplus fit-curves IN [OUT] [--tolerance T] [--snap-angle DEG] [--normal x,y,z] [--min-arc-points N] [--min-arc-sweep DEG] [--json]` |
| Pipeline | `{"Op": "FitCurves", "Tolerance": 0, "SnapAngle": 0, "Normal": [0, 0, 1], "MinArcPoints": 4, "MinArcSweep": 5}` after a chained `Slice` or `Sections`; counters `NumLines`, `NumArcs`, `NumCircles`, `MaxResidual` |
| C | `mio_fit_curves_opts_init(&o); mio_mesh* sk = mio_fit_curves(sections, &o, &report);` then `mio_fit_curves_piece(i, &piece)` (`mio_sketch_piece`, 104 bytes; `MIO_SKETCH_LINE`, `_ARC`, `_CIRCLE`), `mio_fit_curves_frame(i, frame)` and `mio_fit_curves_warnings`; `mio_fit_curves_opts` 88 bytes, `mio_fit_curves_report` 96 |
| Fortran | `sk = mio_fit_curves(sections [, tolerance, snap_angle, normal, min_arc_points, min_arc_sweep, report, stat, errmsg])`, `call mio_fit_curves_piece(i, piece)` (1-based) |
| Julia | `fit_curves(sections; tolerance, snap_angle, normal, min_arc_points, min_arc_sweep)` → `(; mesh, pieces, frames, ...)`, each piece's end as `stop` |
| R | `mio_fit_curves(sections, tolerance, snap_angle, normal, min_arc_points, min_arc_sweep)`, `pieces` a data frame |
| Octave / MATLAB | `mio.fit_curves(sections, ...)`, `pieces` a struct array |
| WASM | `m.fitCurves(sections, {tolerance, snapAngle, normal, minArcPoints, minArcSweep})` |
| MCP | tool `fit_curves(input_path, output_path, tolerance, snap_angle, normal, min_arc_points, min_arc_sweep)` |

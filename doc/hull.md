# Convex hull and planar Delaunay

`convex_hull(points)` and `delaunay_2d(points)` turn scattered points into geometry: a closed hull around a scan, the footprint of a part, or a terrain triangulated from survey points. Both read only the points (a mesh's cells are ignored; an `(n, 2)` or `(n, 3)` array works too), and every decision is an exact geometric predicate (Shewchuk's `orient2d`, `orient3d` and `incircle`), so coplanar, collinear and cocircular points are handled without a tolerance and the result is the same on every platform.

```python
import meshioplusplus as mp

scan = mp.read("part_scan.xyz")
hull, report = mp.convex_hull(scan, return_report=True)    # closed triangle surface
report["volume"], report["area"], report["num_vertices"]

footprint = mp.convex_hull(scan, dimension=2)              # x, y outline: a line loop

survey = mp.read("survey.pcd")
terrain = mp.delaunay_2d(survey)                           # triangulated in x, y; z kept
mp.write("terrain.vtu", terrain)
```

## `convex_hull`

- **3-D** (quickhull, Barber, Dobkin and Huhdanpaa, ACM TOMS 22(4), 1996): a closed `triangle` surface wound outward. A point is outside a face only when `orient3d` says so exactly, and faces that turn out coplanar are merged into one convex facet and re-triangulated from its corners, so a point that lies on a facet or an edge is never a vertex: the hull of a 6 × 6 × 6 lattice has 8 vertices and 12 triangles.
- **2-D** (Andrew's monotone chain over `x, y` with exact `orient2d`): a counter-clockwise closed `line` loop through the strict corners.
- **Dimension.** `dimension=0` (the default) takes 2 for points with two coordinates and for 3-D points whose `z` are all equal, else 3. `dimension=2` on a 3-D cloud takes the outline of its `x, y` projection; the vertices keep their own `z`.
- **Output.** The hull's vertices alone, in ascending input order, with their point data gathered and `hull:point_id` (Int64) naming the input point each one is. Field data rides along; cell data and regions do not (the cells were ignored).
- **Refused by name:** fewer than three points, a non-finite coordinate, collinear points in 2-D, coplanar points in 3-D (naming `dimension=2`).

The report holds `dimension`, `num_vertices`, `num_facets` (triangles, or edges in 2-D), `area` (the enclosed area in 2-D, the surface area in 3-D), `volume` (3-D) and `perimeter` (2-D).

## `delaunay_2d`

The Delaunay triangulation of the points' `x, y` coordinates: no point lies strictly inside any triangle's circumcircle. It is the constrained triangulation of `detail/cdt2d` (Lawson flips over Hilbert-ordered insertion) with the convex hull's edges as its only constraints, and those are Delaunay edges, so the result is the Delaunay triangulation; among cocircular points (a grid) the choice is deterministic.

- **Output.** Every input point, in input order, with its coordinates (a terrain keeps its `z`) and point data, and one counter-clockwise `triangle` block covering the convex hull. Point ids match the input, so the triangulation can be laid over the original data.
- **Duplicates.** A point that coincides in `x, y` with an earlier one is left out of the triangles (`num_duplicates`); it stays in the point list.
- **Range.** The predicates stay exact while no product inside them overflows, so each coordinate must be zero or have a magnitude in [2⁻¹²⁸, 2¹²⁸] (rescale with [`transform`](./transform.md) first otherwise).
- **Refused by name:** fewer than three distinct points, collinear points, a non-finite or out-of-range coordinate.

The report holds `num_triangles`, `num_duplicates` and `num_boundary` (the points on the hull's boundary, its corners and the points on its edges): `num_triangles = 2 n − num_boundary − 2` for `n` distinct points.

The triangulation covers the convex hull: a concave footprint or a hole in the survey is bridged by long, thin triangles. Trim them afterwards (by edge length, with [`select`](./select.md), or by a quality measure), or mesh the outline with [`triangulate_polygon_2d`](./mesh_cad.md) when the boundary is known.

## Speed

Both are sequential. On one core: a million normally distributed points give their hull in about 0.2 s, 200 000 points all on a sphere (every one a vertex) in about 2 s, and 200 000 uniform points their Delaunay triangulation in about 0.5 s.

There is no numpy twin (the [`geodesic_distance`](./geodesic.md) precedent): floating-point numpy cannot reproduce exact predicates on degenerate input, which is exactly where a twin would disagree. The Python functions raise `NotImplementedError` without the compiled core.

## Everywhere else

| Surface | Call |
| --- | --- |
| Python | `convex_hull(points, *, dimension=0, return_report=False)`, `delaunay_2d(points, *, return_report=False)` |
| Pipeline | `{"Op": "ConvexHull", "Dimension": 0}` (counters `Dimension`, `NumVertices`, `NumFacets`, `Area`, `Volume`); `{"Op": "Delaunay2D"}` (counters `NumTriangles`, `NumDuplicates`) |
| Both CLIs | `meshioplusplus hull IN OUT [--dimension 0\|2\|3] [--json]`, `meshioplusplus delaunay IN OUT [--json]` |
| C | `mio_mesh* h = mio_convex_hull(mesh, 0, &report);` (`mio_hull_report`), `mio_mesh* t = mio_delaunay_2d(mesh, &num_triangles, &num_duplicates, &num_boundary);` |
| Fortran | `h = m%convex_hull(dimension=0, num_vertices=nv, volume=v, stat=ierr)`, `t = m%delaunay_2d(num_triangles=nt, stat=ierr)` |
| Julia | `convex_hull(mesh; dimension=0)`, `delaunay_2d(mesh)`; named tuples with `mesh` and the counts |
| R | `mio_convex_hull(mesh, dimension = 0)`, `mio_delaunay_2d(mesh)`; lists with `mesh` and the counts |
| WASM | `m.convexHull(mesh, {dimension})`, `m.delaunay2d(mesh)` |
| MCP | tools `convex_hull(input_path, output_path, dimension)` and `delaunay_2d(input_path, output_path)` |

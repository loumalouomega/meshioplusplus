# Silhouette outlines

`silhouette(mesh, normal)` is the outline of a part seen along a direction: its surface projected onto a plane and unioned, as closed chains with its holes (roadmap 7.6.3). With a clearance it is the outline of a nest, a jig or a foam insert that holds the part loosely or by a press fit, and with `extrude=True` it is that nest's prism. The outline is a chained line mesh, so [`fit_curves`](./fit_curves.md) turns it into lines and arcs and the [DXF](./formats/dxf.md) writer sends it to a laser cutter at 1:1.

<!--pytest-codeblocks:skip-->

```python
import meshioplusplus

outline, report = meshioplusplus.silhouette(part, (0, 0, 1), return_report=True)
report["num_outer"], report["num_holes"], report["area"]

nest = meshioplusplus.silhouette(part, (0, 0, 1), clearance=0.5)  # 0.5 mm of play
sketch = meshioplusplus.fit_curves(nest, tolerance=0.01)
meshioplusplus.write("nest.dxf", sketch)

insert = meshioplusplus.silhouette(part, (0, 0, 1), clearance=-0.2, extrude=True)
meshioplusplus.write("insert.stl", insert)
```

## Options

| Option | Default | Meaning |
|---|---|---|
| `normal` | `(0, 0, 1)` | the viewing direction; the outline lies in the plane normal to it |
| `clearance` | 0 | > 0 grows the outline by this distance (a loose fit), < 0 shrinks it (a press fit) |
| `round_angle` | 180 | a corner the offset opens is rounded with an arc when it turns by less than this many degrees, else kept sharp |
| `arc_segments` | 64 | polygon sides a full turn of arc takes (at least 4) |
| `extrude` | false | return the closed triangle surface of the prism over the outline instead |
| `region` | none | outline only this cell region of a surface |
| `merge_tolerance` | 1e-9 | projected points closer than this times the projected box's diagonal are one point |

The input is a surface (triangles, quads and polygons, fanned as `convert_cells` fans them) or a volume mesh, outlined through its skin. Line and point cells carry no area and are left out; a mesh with nothing else is refused by name.

## How it is computed

The surface's points are projected onto the plane through the origin normal to `normal`, in the frame `lsq_fit`'s `lsq_frame` completes it to, and the union of the projected triangles is found by the shared kernel `detail/planar_union.hpp`:

1. **Candidate edges.** An edge two triangles share is inside the union when the two lie on either side of it in the plane, which is when the surface faces the same way on both. Every other edge may lie on the outline: the surface's boundary edges and its contour generators, where it turns from facing the viewer to facing away. Triangles seen edge-on project to nothing and are left out.
2. **Crossings.** The candidate edges are split where they meet: at the endpoint one touches another at, at the ends of a collinear overlap, and at the point where two cross. The meetings are found with exact orientation tests (`cdt2d_find_crossings`, Shewchuk's predicates), the crossing points computed in double, and the search is repeated until no edges meet.
3. **Coverage.** The split edges constrain a planar Delaunay triangulation (`detail/cdt2d.hpp`). A triangle is covered or not as a whole, and so are its neighbours across an unconstrained edge, so each group of triangles joined that way is decided once, by its largest triangle's centroid tested exactly against the projected triangles found through a bucket grid.
4. **Loops.** The covered triangles' edges whose other side is uncovered are walked with the covered side on the left, turning round each vertex through covered triangles only, so two islands touching at a corner stay two loops. Outer boundaries run counter-clockwise about the normal and holes clockwise. Each loop starts at its lowest point, and loops are ordered by it.

Points closer than the merge tolerance are merged, so two points of the surface that project almost onto each other (the top and bottom edges of a vertical wall) are one. Everything is sequential and the output depends only on the input.

## Clearance

A clearance `c` offsets the outline by the Minkowski sum (`c > 0`) or difference (`c < 0`) with a disk of radius `|c|`, through the same kernel. The region is unioned with, or has subtracted from it, a strip along its boundary on the offset side: a rectangle `|c|` wide on each edge, and at each corner the offset opens (a convex corner when growing, a concave one when shrinking) either an arc or a mitre.

- **An arc** when the corner turns by less than `round_angle`: a polygon whose sides are tangent to the true arc, with `arc_segments` sides a full turn, so the offset is never less than the clearance anywhere: a loose fit always clears, a press fit always grips.
- **A mitre** otherwise: the two offset edges extended until they meet, bevelled when the mitre would be longer than four times the clearance. The default, 180, rounds every corner, which is the true offset; a smaller angle keeps the corners of a machined part sharp.

A corner the offset does not open needs nothing, since the rectangles meet there. A clearance that shrinks the outline to nothing returns no loops, with a warning.

## Output

The outline is a line mesh in the plane at the part's lowest point along the normal, one closed chain per loop in walking order, with the cell data that [`fit_curves`](./fit_curves.md) and the [chained sections](./slice.md#chained-sections) use: the Int64 `section:chain` (the loop) and `section:closed` (1), and the Float64 `section:offset` (the plane's height). The Int64 `silhouette:hole` is 1 on the segments of a hole. With `extrude=True` the result is instead one `triangle` block: the outline's triangulation at the lowest and the highest height along the normal and a wall of two triangles on each outline edge, closed and wound outward. Both carry the field data `silhouette:normal` (the unit normal) and `silhouette:extent` (the lowest and highest heights), and keep the input's `meshio:length_unit`.

The report gives `num_loops`, `num_outer`, `num_holes`, `num_triangles` (the projected triangles of non-zero area), `area` (of the outline's region), `perimeter` (the loops' total length), `min_height`, `max_height` and `warnings`.

The silhouette of a sphere along any axis is one circle of its radius within the tessellation's sag, and growing a unit square by 0.1 with mitred corners gives an area of exactly 1.44 and shrinking it 0.64: `tests/cpp/test_silhouette.cpp` and `tests/python/test_silhouette.py` pin both, and a stress test against a Monte-Carlo estimate of the projected area agrees within its noise over random directions.

## Surfaces

C++ `silhouette(mesh, SilhouetteOptions)` (`operations/silhouette.hpp`), Python `meshioplusplus.silhouette` (C++ core only, no numpy twin), both CLIs (`silhouette`), the `Silhouette` [pipeline](./pipeline.md) step, C `mio_silhouette`, Fortran `mio_silhouette`, Julia `silhouette`, R `mio_silhouette`, Octave/MATLAB `mio.silhouette`, WASM `silhouette` and the MCP tool `silhouette`.

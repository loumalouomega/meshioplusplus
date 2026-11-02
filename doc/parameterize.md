# UV parameterization

`parameterize(mesh)` gives a triangle surface a UV atlas: two texture coordinates per point, in the unit square, as `point_data["uv"]`. It cuts the surface into charts, flattens each without flipping a triangle, and packs the charts side by side. The OBJ, glTF and USD writers turn `uv` into texture coordinates, so the result can be textured, baked or painted in any renderer or DCC tool. Uses: texture baking for digital twins, the [Blender](./blender.md) hand-off, and a 2-D domain for fields that live on a surface.

```python
import meshioplusplus as mp

bunny = mp.read("bunny.stl")
out, report = mp.parameterize(bunny, return_report=True)
report["num_charts"], report["num_flipped"], report["stretch"]
mp.write("bunny.glb", out)    # TEXCOORD_0
mp.write("bunny.obj", out)    # vt, faces as v/vt
mp.write("bunny.usda", out)   # primvars:st
```

## How it works

1. **Charts.** Each edge-connected piece of the surface that is a topological disc (one boundary loop, Euler characteristic 1) is a chart. A closed or holed piece of genus 0 is cut into discs along a cut graph; anything else (a torus, a piece that is not manifold) is halved along the direction its normals spread most. Both cuts are [`fit_cad`](./fit_cad.md)'s. A piece above 8192 triangles is halved before any flattening is tried, since its conformal map is slow and rarely stretches little enough to keep.
2. **Flattening.** Least-squares conformal maps (Lévy, Petitjean, Ray and Maillot, "Least squares conformal maps for automatic texture atlas generation", SIGGRAPH 2002), else mean-value coordinates on a circle (Floater, "Mean value coordinates", CAGD 20(1), 2003), each accepted only when no triangle flips. A chart that flips anyway, or whose L2 stretch exceeds `max_stretch`, is halved again, up to `max_depth` times. Each chart is scaled so its area in the plane equals its area on the surface (one texel density everywhere), mirrored if needed so its triangles keep their winding (counter-clockwise in `uv` seen from the side the normals point to), and turned to its principal axes.
3. **Packing.** The charts' bounding boxes go onto shelves, tallest first, `padding` times the square root of the surface's area apart, and the atlas is scaled uniformly into the unit square, so `uv` keeps the aspect ratio of the surface.

A point on a seam belongs to two or more charts, so it is **split**: the first copy keeps its index (the chart of the first triangle that uses it), the others are appended, and `point_data["parameterize:parent_point"]` names the input point of every output point. Point regions, `line` and `vertex` cells keep pointing at the right points, every point array is copied to the copies, cell data and regions ride along, and `cell_data["parameterize:chart"]` holds each triangle's chart (`-1` on the other blocks). A point in no triangle gets `NaN` coordinates. Weld the copies back with [`clean`](./clean.md) when only the geometry matters.

## What it gives

Measured with the defaults:

| Surface | Charts | Flipped | Stretch | Packing | Seam points | Time |
| --- | --- | --- | --- | --- | --- | --- |
| disk, 384 triangles | 1 | 0 | 1.000 | 0.78 | 0 | 0.01 s |
| icosphere, 1280 triangles | 2 | 0 | 1.081 | 0.40 | 40 | 0.02 s |
| torus, 576 triangles | 4 | 0 | 1.042 | 0.40 | 77 | 0.01 s |
| Stanford bunny, 112k triangles | 47 | 0 | 1.109 | 0.35 | 2936 | 11 s |

`stretch` is Sander, Snyder, Gortler and Hoppe's area-weighted L2 stretch ("Texture mapping progressive meshes", SIGGRAPH 2001) at each chart's own scale: 1 for a map that keeps lengths, as a flat disc's does. `packing_efficiency` is the fraction of the unit square the triangles cover; shelf packing of bounding boxes leaves the corners of round charts empty, so 0.3 to 0.5 is usual. A lower `max_stretch` gives more, smaller and less distorted charts.

## Texture coordinates in the formats

`uv`, Float64 `(n, 2)` point data, is the convention; three formats map it onto their own:

| Format | Written | Read |
| --- | --- | --- |
| [OBJ](./formats/obj.md) | one `vt u v` per point, every face corner `v/vt` (both engines) | `f v/vt` corners: one `uv` per point, a point whose corners name different `vt` rows split as above, with `obj:parent_point` |
| [glTF](./formats/gltf.md) | `TEXCOORD_0` on the triangles, `v` flipped to `1 - v` (glTF's origin is the image's top left); `uv` is then not a `_UV` custom attribute (both engines, the same bytes) | `TEXCOORD_0` as `uv`, flipped back |
| [USD](./formats/usd.md) | `primvars:st`, `texCoord2f[]`, vertex interpolation | `primvars:st` as `uv` |

## Scope

The surface is the mesh's `triangle` blocks; `line` and `vertex` blocks ride along. A volume cell is refused naming [`extract_surface`](./extract_surface.md); a quad, a polygon or a higher-order cell naming `convert_cells(mode="simplexify")` or `linearize`. A negative `padding` or `max_depth`, or a `max_stretch` below 1, is refused. Everything is serial and in input order, so the atlas is the same on every backend and thread count. There is no numpy twin (the [`geodesic_distance`](./geodesic.md) precedent): the conformal map is an iterative solve.

## Everywhere else

| Surface | Call |
| --- | --- |
| Python | `parameterize(mesh, *, padding=0.02, max_depth=8, max_stretch=2.0, return_report=False)`; report `num_charts`, `num_flipped`, `stretch`, `packing_efficiency`, `num_seam_points` |
| C++ | `ParameterizeResult parameterize(const Mesh&, const ParameterizeOptions&)` in `operations/parameterize.hpp` |
| Pipeline | `{"Op": "Parameterize", "Padding": 0.02, "MaxDepth": 8, "MaxStretch": 2}`; counters `NumCharts`, `NumFlipped`, `Stretch`, `PackingEfficiency` |
| Both CLIs | `meshioplusplus parameterize IN OUT [--padding F] [--max-depth N] [--max-stretch S] [--json]` |
| C | `mio_parameterize_opts_init(&o); mio_mesh* out = mio_parameterize(mesh, &o, &report);` |
| Fortran | `out = m%parameterize(max_stretch=1.5_real64, num_charts=nc, stat=ierr)` |
| Julia | `r = parameterize(mesh; max_stretch=1.5)`; `r.mesh`, `r.num_charts`, `r.stretch` |
| R | `r <- mio_parameterize(mesh, max_stretch = 1.5)`; `r$mesh`, `r$num_charts`, `r$stretch` |
| WASM | `m.parameterize(mesh, {padding, maxDepth, maxStretch})` returns `{mesh, numCharts, numFlipped, stretch, packingEfficiency, numSeamPoints}` |
| MCP | tool `parameterize(input_path, output_path, padding, max_depth, max_stretch)` |

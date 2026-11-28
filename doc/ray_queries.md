# Ray casting and point location

Two batch queries on a mesh, both returning arrays. `cast_rays` shoots rays at a **surface** and reports the nearest hit of each; `locate_points` finds the **cell** of a volume mesh that contains each point and the barycentric weights inside it. They are the lookups [`thickness`](/thickness), [`distance_to_surface`](/sdf) and [`interpolate`](/interpolate) already do internally, as operations of their own.

```python
import numpy as np
import meshioplusplus as mio

sphere = mio.uv_sphere(1.0, 64, 32)
r = mio.cast_rays(sphere, origins=[[3, 0, 0], [0, 3, 0]], directions=[-1, 0, 0])
print(r["hit"], r["distance"], r["cell"])

box = mio.box(dims=(2, 2, 2), size=(2.0, 2.0, 2.0))
loc = mio.locate_points(box, [[0.5, 0.5, 0.5], [9, 9, 9]])
print(loc["cell"], loc["nodes"][0], loc["barycentric"][0])
```

```sh
meshioplusplus cast-rays surface.stl origins.ply hits.vtu --direction 0,0,-1
meshioplusplus locate volume.vtu points.ply located.vtu
```

## cast_rays

`cast_rays(mesh, origins, directions, max_distance=0, region=None, all_hits=False)` shoots every ray `origin + t * direction` at the **triangle, quad and rectangular polygon cells** of `mesh` (a quad or polygon is fanned into triangles) and reports the nearest hit (or, with `all_hits=True`, every crossing sorted by ascending `t`) with `0 < t <= max_distance`. The direction need not be a unit vector (it is normalised, so `t` is a length), and one direction may serve every ray. A volume mesh is refused by name, pointing at `extract_surface`; a higher-order block at `linearize`. `region` restricts the surface to a named cell region.

| Key | Meaning |
| --- | --- |
| `hit` | whether the ray hit (nearest-hit mode) |
| `distance` | the distance to the hit, `NaN` on a miss; length `m` in `all_hits=True` mode |
| `point` | the hit point, `NaN` on a miss; shape `(m, 3)` in `all_hits=True` mode |
| `normal` | the unit triangle normal on its winding side `(b - a) x (c - a)`, `NaN` on a miss; shape `(m, 3)` in `all_hits=True` mode |
| `cell` | the global (block-major) index of the **input cell** hit, `-1` on a miss; length `m` in `all_hits=True` mode |
| `ray_offset` | CSR offsets of length `n_rays + 1` slicing the flat hit arrays per ray (present when `all_hits=True`) |
| `num_hits` | how many rays hit (nearest-hit mode), or total crossings (`all_hits=True`) |

A ray whose origin is on the surface does not hit it at `t = 0`. A ray through a shared edge or vertex hits once, at the lowest-numbered triangle: the traversal's barycentric tolerance of `1e-10` closes the seams, and `all_hits=True` deduplicates crossings within `1e-9 * max(1, |t|)` at the same vertex or edge. The rays are traced in parallel over a bounding-volume hierarchy and the result is identical on every backend and thread count.

**Not offered yet:** a reusable query handle that keeps the hierarchy between batches (roadmap 9.3.2, on request); each call builds its hierarchy.

## locate_points

`locate_points(mesh, points, tolerance=1e-12)` reduces `mesh` to simplices exactly as `interpolate`'s barycentric mode does (`convert_cells(Simplexify)`: a hexahedron, prism, pyramid, polyhedron or quadrilateral is split; a higher-order cell is linearised), looks each point up in a grid of the simplices' bounding boxes and takes the **first containing simplex in ascending order**, so a point on a shared face belongs to the lowest-numbered cell. A simplex contains a point when every barycentric weight is at least `-tolerance`. The search is the one `interpolate` uses (`detail/point_locator.hpp`), which is why `sum(w[k] * f[nodes[k]])` equals `interpolate`'s value of a point array `f` bit for bit.

| Key | Meaning |
| --- | --- |
| `cell` | the global index of the **input cell** that contains the point, `-1` outside |
| `nodes` | four point ids of the containing simplex (the input mesh's own ids; a triangle's fourth is `-1`), `-1` outside |
| `barycentric` | four weights with respect to `nodes` (a triangle's fourth is `0`), `0` outside |
| `num_located` | how many points were inside |

Volume meshes (tetrahedra after simplexification) and **planar** triangle meshes are located. A triangle's containment test is made in the xy-plane, but the grid lookup that precedes it is three-dimensional, so a point must lie in the mesh's plane, within a grid cell. The lookup also means `tolerance` widens the barycentric test only: a point outside the mesh by more than the rounding that keeps it in the same grid cell as a simplex is reported outside, whatever the tolerance. A surface in space is not located by containment: [`distance_to_surface`](/sdf) projects onto it and `cast_rays` intersects it.

## Surfaces

The operations are in C++ (`operations/cast_rays.hpp`, `operations/locate_points.hpp`), Python (core only, no numpy fallback), both CLIs (`cast-rays` and `locate`, which write the hit and located points as a vertex mesh), the [MCP server](/mcp) (`cast_rays`, `locate_points`), the [C API](/c_api) (`mio_cast_rays`, `mio_locate_points`, caller-allocated output arrays), [Fortran](/fortran) (`mesh%cast_rays`, `mesh%locate_points`), [Julia](/julia), [R](/r) (`mio_cast_rays`, `mio_locate_points`), [Octave / MATLAB](/octave) (`mio.cast_rays`, `mio.locate_points`) and [WebAssembly](/wasm) (`castRays`, `locatePoints`).

Indices follow each binding's convention: 0-based and `-1` for none in C, C++, Python and WebAssembly; 1-based and `0` for none in Fortran, Julia, R and Octave.

# Offset, thicken and hollow

`offset_surface(mesh, distance)` moves a surface along its normals so every face's plane moves by `distance`; `thicken(surface, thickness, side=...)` turns an open surface into a closed shell; `hollow(solid, wall)` gives a closed surface an inner wall (roadmap 7.8.4). The surface is triangles and quads; another cell type is refused by name.

<!--pytest-codeblocks:skip-->

```python
import meshioplusplus

grown, report = meshioplusplus.offset_surface(part, 0.5, return_report=True)
report["max_move"], report["num_folded"]
shell = meshioplusplus.thicken(scan, 2.0, side="inward")  # the scan stays the outside
cup = meshioplusplus.hollow(solid, 1.5)
```

## The direction and the distance

Each point moves along its offset direction: the most normal normal of Aubry and Löhner (2009), the angle-weighted mean of its faces' normals stepped towards the face that sees it least, the direction [`boundary_layer`](./boundary_layer.md) grows its layers along (one function, `detail/offset_normals.hpp`, computes both, so the first layer of a boundary layer and an offset by its height agree bit for bit). It moves by the distance over the least cosine between that direction and a face around it, at most four times the distance: both faces of a sharp edge keep the distance, so a cube offset by `d` is the cube grown by `d`, its corners moved by `sqrt(3) * d` (`tests/cpp/test_offset_surface.cpp`). `smooth_iterations` smooths the directions as `boundary_layer`'s `normal_smooth_iterations` does.

A positive distance moves along the faces' normals (their winding; [`repair`](./repair.md) orients a surface outward, [`flip_normals`](./flip_normals.md) the other way). A point that no direction sees from all its faces within a cosine of 0.05 (a saddle sharper than about 87 degrees) moves by the distance alone and is counted in `num_not_visible`; a point whose stretch reached its bound in `num_clamped`. The offset of a concave region by more than its radius folds: each face whose normal turns over is counted in `num_folded`, with a warning. There is no collision guard: the [`compute_sdf`](./sdf.md) + [`isosurface`](./isosurface.md) recipe is the robust, remeshing alternative for large offsets.

## Thicken and hollow

| Operation | Layers | Rim |
|---|---|---|
| `thicken(side="inward")` (default) | the surface, and its offset by `-thickness` reversed | quads join each boundary edge across |
| `thicken(side="outward")` | the offset by `+thickness`, and the surface reversed | the same |
| `thicken(side="both")` | the offsets by `+thickness / 2` and `-thickness / 2`, the second reversed | the same |
| `hollow(wall)` | the surface, and its offset by `-wall` reversed | none; an open surface is refused |

The points are the upper layer's then the lower's, the point data duplicated; the cells are one `triangle` and one `quad` block, and the Int64 cell data `thicken:side` is 0 on the upper layer, 1 on the lower and 2 on the rim. The input's cell data and regions are not carried onto the shell (a warning names them). A thickened hemisphere is closed and consistently wound, its volume that of the two half shells to 2%, and [`thickness`](./thickness.md) reads the requested wall to 1% (its 2nd percentile).

The report (`return_report=True`): `num_not_visible`, `num_clamped`, `num_folded`, `num_rim_quads` (thicken), `max_move` and `warnings`.

## Surfaces

C++ `offset_surface`, `thicken`, `hollow` and `thicken_side_from_name` (`operations/offset_surface.hpp`); Python `meshioplusplus.offset_surface`, `thicken` and `hollow` (C++ core only); both CLIs (`offset-surface IN OUT --distance D`, `thicken IN OUT --thickness T [--side inward|outward|both]`, `hollow IN OUT --wall W`, each with `[--smooth N] [--json]`); the `OffsetSurface` (`Distance`), `Thicken` (`Thickness`, `Side`) and `Hollow` (`Wall`) [pipeline](./pipeline.md) steps, each with `SmoothIterations`; C `mio_offset_surface`, `mio_thicken` (the `mio_thicken_side` enum) and `mio_hollow` (`mio_offset_report`, 64 bytes; `mio_offset_warnings`); Fortran, Julia, R, Octave/MATLAB, WASM `offsetSurface`, `thicken` and `hollow`; and the MCP tools `offset_surface`, `thicken` and `hollow`.

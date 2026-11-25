# Point-cloud conditioning

`filter_points` removes outliers from a point cloud and thins it. Scans read from PLY, PCD, XYZ or LAS carry noise and far more points than a reconstruction or a registration needs; this is the step between reading one and [`estimate_normals`](/point_normals) or [`reconstruct_surface`](/reconstruct_surface).

```python
import meshioplusplus as mio

scan = mio.read("scan.ply")
clean, report = mio.filter_points(
    scan, statistical_k=16, std_ratio=2.0, voxel_size=0.5, return_report=True
)
print(report["num_input"], len(clean.points), report["num_removed_statistical"])
```

```sh
meshioplusplus filter-points scan.ply clean.ply --outliers statistical --voxel 0.5
```

## The four steps

Up to four steps run **in this order**, each seeing the points the one before kept. A step is off at its default of `0`, and at least one must be on.

| Step | Options | Rule |
| --- | --- | --- |
| Statistical outliers | `statistical_k`, `std_ratio` | each point's mean distance to its `statistical_k` nearest other points is compared with the mean `mu` and the population standard deviation `sigma` of those means over all points; a point above `mu + std_ratio * sigma` goes (Rusu et al., *Towards 3D point cloud based object maps for household environments*, Robotics and Autonomous Systems 2008) |
| Radius outliers | `radius`, `min_neighbors` | a point with fewer than `min_neighbors` **other** points within `radius` goes |
| Voxel thinning | `voxel_size` | the points are binned into cubes of that side anchored at the bounding box's minimum corner, and one **real point** per occupied cube survives: the one nearest the mean of its cube's points, ties to the lower id. The output is a subset of the input, with every array intact; the cubes are visited in ascending (x, y, z) key order, which is the step's output order |
| Poisson-disk thinning | `min_distance` | greedy dart throwing in index order, with no random numbers (after Bridson, *Fast Poisson disk sampling in arbitrary dimensions*, SIGGRAPH Sketches 2007): a point is kept when no already-kept point lies closer than `min_distance`. The kept points are at least that far apart, and every dropped point is within it of a kept one |

All searches are exact and deterministic (ties to the lower id), the per-point work runs in parallel and every reduction is serial in input order: the result is identical on every backend and thread count.

## What the result is

A **point cloud**: the surviving points, every per-point array gathered with them (dtype kept), one `vertex` block, and `point_map` giving each output row's input point. Field data is kept. Point regions follow the surviving points (a group that loses every point stays, empty). Cells and Cell and Side regions are dropped, with a warning.

| Report key | Meaning |
| --- | --- |
| `point_map` | each output point's input point (0-based; 1-based in Fortran, Julia, R and Octave) |
| `num_input` | the input point count |
| `num_removed_statistical`, `num_removed_radius`, `num_removed_voxel`, `num_removed_poisson` | what each step removed |
| `warnings` | cells or regions dropped, a step skipped |

A negative or non-finite parameter, a non-finite coordinate, no step enabled, or a voxel size so small that the cloud's extent divided by it exceeds 10^15 throws by name.

## Not offered yet

The voxel step keeps a real point rather than the cube's mean point, so the result never invents a point. A step that averages the coordinates (and the float arrays) of each cube is the remainder of roadmap 9.9, together with moving `subsample_points`' count-driven `grid` and `farthest` methods into the core; the Python `subsample_points` keeps them meanwhile.

## Surfaces

C++ `filter_points` (`operations/filter_points.hpp`), Python (core only, no numpy fallback), both CLIs (`filter-points`), the [pipeline](/pipeline) step `FilterPoints`, the [MCP server](/mcp) tool `filter_points`, the [C API](/c_api) (`mio_filter_points`), [Fortran](/fortran) (`mesh%filter_points`), [Julia](/julia), [R](/r) (`mio_filter_points`), [Octave / MATLAB](/octave) (`mio.filter_points`) and [WebAssembly](/wasm) (`filterPoints`).

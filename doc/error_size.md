# Sizes from an error estimate — `error_size`

`error_size` turns a per-cell error indicator into a per-point size field. [`estimate_error`](./error.md) measures where a field is badly resolved and writes `error:zz`, one value per cell; [`mesh_cad_surface`](./mesh_cad.md)'s `background` reads a size per point. `error_size` joins the two, so an adaptation loop re-meshes from the CAD model instead of moving points off it.

```python
import meshioplusplus as mp

model = mp.read("part.cad.json")
mesh = mp.mesh_cad_surface(model, 0.2)
mesh.point_data["f"] = solve(mesh)  # a field on the mesh
mesh = mp.estimate_error(mesh, "f")  # cell_data error:zz
sized, report = mp.error_size(mesh, target_error=1e-3, min_size=0.01, return_report=True)
adapted = mp.mesh_cad_surface(model, 0.2, background=sized)  # reads point_data "size"
```

## How a size is chosen

The procedure is the adaptive one of Zienkiewicz and Zhu (1987). With `N` cells evaluated and a target global error `eta`, the error each cell may carry when the error is spread evenly is `e = eta / sqrt(N)`. A cell of size `h_K` (its diameter, the largest distance between two of its nodes) whose indicator is `eta_K` asks for `h_K (e / eta_K)^(1/p)`, where `p` is the rate at which the indicator falls with `h`: the element order for the recovered-gradient error, 1 for linear elements. A cell whose indicator is zero asks for no change of its own and is bounded by `max_size` alone.

The size is clamped to `[min_size, max_size]` (0 is no bound) and written per point as the least of its cells' sizes: the finer request wins, as every size source of `mesh_cad_surface` refines. A point with no evaluated cell (a cell whose indicator or size is not finite is skipped) reads `NaN`, which `background` treats as no size there.

Geometry, connectivity, regions, property sets and every existing array pass through; only the size array is added. The result is byte-identical across mesh backends and thread counts.

On the unit sphere with a `tanh` layer across it (notebook 42), the indicator of the layer's field on a mesh at `h = 0.1`, asked for 1.5 times that mesh's global error, re-meshes the sphere with 3,137 points at an L2 interpolation error of 1.54e-2, where a uniform mesh at `h = 0.05` needs 11,273 points for 1.59e-2.

![The sphere meshed uniformly at h = 0.05 and sized by error_size from the layer's error indicator, coloured by the layer's field](/images/mesh_cad_error_size.png)

## Options and report

| option | meaning |
|---|---|
| `error_array` | the one-component `cell_data` indicator (`error:zz`) |
| `target_error` | the target global error `eta`, in the indicator's norm (required, positive) |
| `order` | the rate `p` (1) |
| `min_size`, `max_size` | the bounds on the size written (0: none) |
| `output` | the `point_data` array written (`size`, what `background` reads by default) |
| `overwrite` | allow replacing an existing array of that name |

The report has `num_cells` (evaluated), `num_skipped` (cells whose indicator or size is not finite), `num_clamped` (cell sizes the bounds clamped), `min_size` and `max_size` (over the points written). A missing indicator, one that is `point_data` or has more than one component, a target error or order that is not finite and positive, a negative bound, `min_size > max_size` or an output-name collision is refused by name.

## Bindings

| surface | entry |
|---|---|
| Python | `meshioplusplus.error_size(mesh, error_array="error:zz", *, target_error, order=1.0, min_size=0.0, max_size=0.0, output="size", overwrite=False, return_report=False)`; a numpy twin answers when the core is absent |
| C++ | `error_size(const Mesh&, const ErrorSizeOptions&)` → `ErrorSizeResult` (`operations/error_size.hpp`) |
| C | `mio_error_size(mesh, &opts, &report)` with `mio_error_size_opts` (80 bytes, through `mio_error_size_opts_init`) and `mio_error_size_report` (64 bytes) |
| Fortran | `m%error_size(target_error, [error_array, order, min_size, max_size, output, overwrite, …])` |
| Julia | `error_size(mesh; target_error, error_array, order, min_size, max_size, output, overwrite)` |
| R | `mio_error_size(mesh, target_error, error_array, order, min_size, max_size, output, overwrite)` |
| WASM | `errorSize(mesh, {targetError, errorArray, order, minSize, maxSize, output, overwrite})` |
| CLI | `data error-size INFILE OUTFILE --target-error E [--array NAME] [--order P] [--min-size H] [--max-size H] [--output NAME] [--overwrite] [--json]`, both CLIs |
| MCP | `error_size(input_path, output_path, target_error, ...)` |
| Pipeline | `ErrorSize` after an `EstimateError` step, with `TargetError`, `Array`, `Order`, `MinSize`, `MaxSize`, `Output` |

## References

- O. C. Zienkiewicz and J. Z. Zhu, "A simple error estimator and adaptive procedure for practical engineering analysis", *International Journal for Numerical Methods in Engineering* 24 (1987).

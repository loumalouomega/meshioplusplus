# Colorize

`colorize(mesh, field)` bakes a scalar field into 8-bit colours through a colormap, so the colouring travels with the mesh into files and tools that show colours but not fields: MeshLab, a PLY viewer, a glTF scene, a 3-D print slicer. A deviation map is the usual input; any point or cell field works.

```python
import meshioplusplus as mp

measured = mp.deviation(scan, nominal)
coloured = mp.colorize(measured, "deviation:distance", cmap="coolwarm",
                       symmetric=True, steps=10)
mp.write("deviation.ply", coloured)                        # vertex colours
mp.write("deviation.glb", coloured, color_by="rgb")        # the same colours in glTF
```

## The colours

The field's values `v` map to `t = (v - vmin) / (vmax - vmin)`, clamped to `[0, 1]` (a degenerate range maps to the middle of the map), then through one of the built-in 256-entry colormaps (`viridis`, `coolwarm`, `turbo`, `magma`, `inferno`, `plasma`, `grey`, `synthwave` and their `_r` reverses) that the [renderer](/tui) and the [glTF](/formats/gltf) writer share. The range is the field's finite minimum and maximum unless `vmin` and `vmax` give it; `symmetric=True` takes `+-max |v|` so zero sits mid-map, as a signed deviation wants. A vector field is reduced to its magnitude unless `component` picks one entry. A non-finite value takes `nan_color` (`#rrggbb`, default `#808080`).

**Bands.** With `steps` of 2 or more, `t` is replaced by the centre of the equal band it falls in, so the map reads as `steps` flat colours at fixed steps: a tolerance band diagram rather than a gradient. The last band holds `t = 1`.

## What comes back

The input, with three UInt8 arrays named `red`, `green` and `blue` beside the field, on its location: point data for a point field, cell data for a cell field. Points, cells and every other array are unchanged. Those are the names PLY's vertex and face colour properties take, and the PLY writer writes them as `uchar` properties under them; the glTF writer takes them as `COLOR_0` with `color_by="rgb"`. The report gives the `vmin` and `vmax` used and the `location`.

## Core only

The colouring runs in the compiled core and has no numpy twin. Without the core, the Python entry point raises `NotImplementedError` naming it.

## CLI

```sh
meshioplusplus colorize IN OUT --field F [--cmap NAME] [--vmin V --vmax V] [--symmetric]
                        [--steps N] [--nan-color #rrggbb] [--component K] [--json]
```

## Pipeline

```json
{ "Op": "Colorize", "Field": "deviation:distance", "Cmap": "coolwarm", "VMin": 0, "VMax": 0, "Symmetric": true, "Steps": 10, "NanColor": "#808080", "Component": -1 }
```

`VMin >= VMax` takes the automatic range; the step reports `VMin` and `VMax`. See [pipelines](/pipeline).

## Bindings

C `mio_colorize(mesh, &opts, &report)` with `mio_colorize_opts` (96 bytes, `mio_colorize_opts_init`) and `mio_colorize_report` (56 bytes, `location` 0 point or 1 cell); Fortran `mio_colorize(mesh, field[, cmap, vmin, vmax, symmetric, steps, nan_color, component, report])`, Julia `colorize(mesh, field; ...)`, R `mio_colorize(mesh, field, ...)` and Octave/MATLAB `mio.colorize(mesh, field, ...)` (1-based `component`), WASM `colorize(mesh, {field, cmap, vmin, vmax, symmetric, steps, nanColor, component})` (0-based) and the MCP tool `colorize`.

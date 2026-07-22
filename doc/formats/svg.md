# SVG (`.svg`)

[Scalable Vector Graphics](https://www.w3.org/TR/SVG/) output — a **write-only** visualization format that draws the mesh edges. Flat 2D meshes draw directly; genuinely 3D meshes are **rendered**: the boundary skin of any volume cells is extracted (see [`extract_skin`](../extract_skin.md)) and projected through an orthographic camera, painted back-to-front.

| | |
|---|---|
| **Format name** | `svg` |
| **Extensions** | `.svg` |
| **Read / Write** | — / ✓ |
| **Extra dependencies** | — |

## Reading & writing

There is no reader — `register_format` is called with `read=None`. Full write signature:

```python
import meshioplusplus

meshioplusplus.svg.write(
    "out.svg", mesh,
    float_fmt=".3f",
    stroke_width=None,
    image_width=100,
    fill="#c8c5bd",
    stroke="#000080",
    azimuth=45.0,
    elevation=35.264389682754654,
    roll=0.0,
    color_by=None,
    component=None,
    cmap="viridis",
    vmin=None,
    vmax=None,
    nan_color="#808080",
    colorbar=False,
)
```

- **`float_fmt`** — coordinate number format.
- **`stroke_width`** — edge stroke width; if `None` (default), auto-computed as 1% of the mesh's on-canvas width.
- **`image_width`** — output SVG width in user units, default `100` (not `None`) — deliberately non-trivial because some SVG viewers (e.g. `eog`) mis-render images whose natural width is close to `1` unit.
- **`fill`** / **`stroke`** — cell fill and edge colours, defaulted to match ParaView's default rendering colours (per an inline source comment).
- **`azimuth`** / **`elevation`** / **`roll`** — orthographic camera angles in degrees, used only for genuinely 3D input. The view direction is `(cos el·cos az, cos el·sin az, sin el)` with z as the up reference; `roll` rotates the image in-screen. The defaults (`45`, `atan(1/√2) ≈ 35.264°`) look down the `(1,1,1)` diagonal — the classic CAD isometric view.
- **`color_by`** … **`colorbar`** — data-driven colouring, see below. With `color_by` unset (the default) the output is byte-identical to previous releases.

### Data-driven colouring

`color_by` names a `point_data` or `cell_data` array; each drawn face then gets its own `fill` attribute (overriding the document-level `path { fill: … }` rule) instead of the flat `fill` colour.

```python
meshioplusplus.svg.write(
    "quality.svg",
    meshioplusplus.attach_quality(mesh),
    color_by="quality:scaled_jacobian",
    cmap="viridis",
    colorbar=True,
)
```

![a bracket coloured by scaled Jacobian](/images/color_by_quality.svg)

- **Point vs cell.** The name is looked up in `point_data` first, then `cell_data`; present in neither is an error naming what *is* available. **Point data** colours a face by the **mean of its corner values**; **cell data** by **its owning cell's value**.
- **Volume meshes.** On the projected path the drawn faces are skin facets, not input cells. When colouring by cell data the skin is therefore extracted **with provenance** (the `surface:parent_cell` array from [`extract_surface`](../extract_surface.md)), so each facet gets the value of the volume cell it actually bounds. A 2D surface mesh needs no indirection — its faces *are* its cells.
- **`component`** picks one component of a multi-component array; without it a multi-component row reduces to its **magnitude**. A 1-D scalar array stays scalar.
- **`cmap`** is one of the built-in colormaps — `viridis`, `coolwarm`, `turbo` — baked into the core as 256-entry lookup tables. **No extra dependency**: matplotlib is needed only to *regenerate* the tables (`tools/gen_colormaps.py`), never to use them.
- **`vmin` / `vmax`** set the mapped range. Left unset, the range is the **finite minimum and maximum among the drawn faces** — so the visible figure spans the whole colorbar. Note this differs from ParaView, which ranges over the whole array including cells that are never drawn. Values outside the range are **clamped**. `vmin > vmax` is an error; `vmin == vmax` (which an all-constant array produces) maps everything to the middle of the colormap rather than dividing by zero.
- **`nan_color`** is used verbatim, without going through the colormap, for any face whose value is NaN or infinite. Such values are also excluded from the auto range.
- **`colorbar`** appends a 32-swatch vertical gradient with `vmin`/`vmax` labels to the right of the figure. It widens **only the `viewBox`** — the scaling factor, the stroke width and every mesh coordinate are computed from the unmodified width, so turning it on moves nothing already on the canvas. The swatches are `<rect>`, not `<path>`, so the document-level `path` rule does not stroke them.

### 3D input

A mesh whose points have a non-zero z extent takes the 3D rendering path: if it contains supported volume cells, `extract_skin(mesh, linearize=True)` runs first (higher-order faces collapse to corner triangles/quads); a 3D *shell* mesh (only surface cells) is projected as-is (`triangle6`/`quad8`/`quad9` are corner-linearized). Faces are then sorted back-to-front by view-space centroid depth (painter's algorithm, stable sort) and emitted with the same `<path>` templates as the flat path. There is no backface culling and no shading in v1 (a fixed fill colour; shading is future work).

## File structure

A single `<svg>` root containing one `<path>` element per drawable cell (**not** `<polygon>`) — chosen deliberately: the comment in the source notes that `svgo` (a common SVG optimizer) converts `<polygon>`s to `<path>`s but drops style information when it does so, so meshio++ emits paths directly to sidestep that.

Path `d` templates (space-separated coordinate pairs, `float_fmt`-formatted):

| cell type | path template |
|---|---|
| `line` | `M x0 y0L x1 y1` (open, no closing `Z`) |
| `triangle` | `M x0 y0L x1 y1L x2 y2Z` |
| `quad` | `M x0 y0L x1 y1L x2 y2L x3 y3Z` |

If `points.shape[1] == 3` and every z coordinate is `~0` (`atol=1e-14`), the mesh is treated as flat 2D and drawn exactly as in previous releases (byte-identical); otherwise the 3D projected path above applies. The y-coordinate is flipped (`max_y + min_y - y`) to convert from the mesh/math convention (y-up) to SVG's screen convention (y-down) — in the 3D path the flip applies to the projected coordinates.

## Cell types

`line`, `triangle`, `quad` (plus, on the 3D path, corner-linearized `triangle6`/`quad8`/`quad9` and the volume types accepted by [`extract_skin`](../extract_skin.md)). Any other cell block present in the mesh is **silently dropped — no warning at all** (unlike most other meshio++ writers' warn-and-skip convention for unsupported cell types).

## Data mapping

No data array is written to the file. One array can be *read* to drive the face colours via `color_by` (see above); everything else — and all of `field_data` — is ignored, so only geometry and connectivity affect an uncoloured figure.

## Quirks & limitations

- No diagonal/winding correction on `quad` cells — a "crossed" (non-convex, bowtie) node ordering renders incorrectly with no error raised.
- Unsupported cells vanish from the output silently.
- The painter's algorithm sorts whole faces by centroid depth — mutually intersecting faces (which a closed skin never has) can stack in the wrong order; there is no per-pixel depth test.
- Write-only; there is no way to read an SVG back into a `Mesh`.
- Colouring is a **Python + C++-direct + CLI** feature. The C API, Fortran and WebAssembly surfaces reach this writer through the shared registry, whose `(path, mesh)` writer entries structurally cannot carry parameters, so they always emit the fixed default styling — a documented gap of the same kind as the point/cell-set gaps in `diff`/`merge`/`split`.
- Only *one* array can drive the colours; there is no multi-field or per-block colouring.
- The colorbar's labels are the range endpoints only — there are no intermediate ticks.

## Notes

- Backed by the **C++ core** (`write_svg`) with a pure-Python fallback: `meshioplusplus.svg.write` uses the C++ writer for real file paths and falls back to Python for file-object/buffer targets or on any error. Registered in the shared dispatch registry, so it is also reachable from the WASM, C API, and Fortran flat bindings (write-only, fixed default styling, **the default isometric camera and no colouring** — per-call options are exposed only through the Python `write`, the direct C++ call, and the CLI).
- The C++ writer and the pure-Python reference are **byte-identical**, on the flat path, the projected 3D path and the coloured path alike. `tests/test_svg.py` asserts exact `read_bytes()` equality across a matrix of colouring options and additionally pins the provenance mapping, the corner averaging, the auto range, clamping, `nan_color` and the colorbar's append-only geometry. `src/cpp/tests/test_svg_tikz.cpp` covers the C++ writer directly, including the golden `fill="#440154"` spelling and the four invalid-option errors.
- Colouring shares its resolution layer with the TikZ writer: `src/cpp/include/meshioplusplus/detail/face_color.hpp` and its twin `src/python/meshioplusplus/_facecolor.py`, over the colormap tables in `detail/colormap.{hpp,cpp}` / `_colormap.py`.

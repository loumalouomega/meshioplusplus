# TikZ (`.tikz`)

[TikZ/PGF](https://tikz.dev/) output — a **write-only** visualization format that draws the mesh cells as a LaTeX figure. Flat 2D meshes draw directly; genuinely 3D meshes are **rendered**: the boundary skin of any volume cells is extracted (see [`extract_skin`](../extract_skin.md)) and projected through an orthographic camera, drawn back-to-front. By default it emits a standalone, directly `pdflatex`-compilable document; it can also emit a bare `tikzpicture` snippet for `\input` into a larger LaTeX document. It is the LaTeX counterpart to the [SVG](./svg.md) writer (the project logo — the Stanford bunny — is produced through this machinery, see `logo/gen_logo_tikz.py`).

| | |
|---|---|
| **Format name** | `tikz` |
| **Extensions** | `.tikz` |
| **Read / Write** | — / ✓ |
| **Extra dependencies** | — |

## Reading & writing

There is no reader — `register_format` is called with `read=None`. Full write signature:

```python
import meshioplusplus

meshioplusplus.tikz.write(
    "out.tikz", mesh,
    float_fmt=".6f",
    standalone=True,
    line_width=None,
    fill="gray!30",
    draw="black",
    scale=None,
    azimuth=45.0,
    elevation=35.264389682754654,
    roll=0.0,
)
```

- **`float_fmt`** — coordinate number format.
- **`standalone`** — if `True` (default), wrap the `tikzpicture` in a full `\documentclass{standalone}` + `\usepackage{tikz}` document that compiles directly with `pdflatex`. If `False`, emit only the `\begin{tikzpicture}…\end{tikzpicture}` environment for `\input` into an existing document.
- **`line_width`** — TikZ line width for edges, e.g. `"0.4pt"`; if `None` (default), TikZ's own default width is used. When set it is applied both on the `tikzpicture` options and on each `\draw`.
- **`fill`** — [xcolor](https://ctan.org/pkg/xcolor) fill spec for the filled faces (triangles/quads), e.g. `"gray!30"`, `"blue!20"`.
- **`draw`** — xcolor spec for the edge stroke.
- **`scale`** — optional `\begin{tikzpicture}[scale=…]` factor; if `None` (default), coordinates are emitted verbatim and no `scale` key is added.
- **`azimuth`** / **`elevation`** / **`roll`** — orthographic camera angles in degrees, used only for genuinely 3D input; same semantics and defaults (the classic CAD isometric view) as the [SVG writer](./svg.md).

### 3D input

A mesh whose points have a non-zero z extent takes the 3D rendering path: supported volume cells are skin-extracted first (`extract_skin(mesh, linearize=True)`); a 3D shell mesh is projected as-is (`triangle6`/`quad8`/`quad9` corner-linearized). Faces are sorted back-to-front by view-space centroid depth (painter's algorithm) and drawn with the same `\draw` templates as the flat path — nearer filled faces cover the hidden edges of farther ones. No backface culling, no shading in v1.

## File structure

One `\draw` command per drawable cell inside a single `tikzpicture` environment. Each cell's vertices are emitted as `(x,y)` coordinate pairs (`float_fmt`-formatted) joined by TikZ's `--` path operator:

| cell type | `\draw` template |
|---|---|
| `line` | `\draw[draw=…] (x0,y0) -- (x1,y1);` (open, no `cycle`) |
| `triangle` | `\draw[fill=…, draw=…] (x0,y0) -- (x1,y1) -- (x2,y2) -- cycle;` |
| `quad` | `\draw[fill=…, draw=…] (x0,y0) -- (x1,y1) -- (x2,y2) -- (x3,y3) -- cycle;` |

If `points.shape[1] == 3` and every z coordinate is `~0` (`atol=1e-14`), the mesh is treated as flat 2D and drawn exactly as in previous releases (byte-identical); otherwise the 3D projected path above applies.

Unlike the SVG writer, the y-coordinate is **not** flipped — TikZ/PGF already uses the math convention (y grows upward), so mesh (or projected) coordinates map straight onto the canvas.

## Cell types

`line`, `triangle`, `quad` (plus, on the 3D path, corner-linearized `triangle6`/`quad8`/`quad9` and the volume types accepted by [`extract_skin`](../extract_skin.md)). Any other cell block present in the mesh is **silently dropped** (matching the SVG writer's behaviour).

## Data mapping

None consumed — `point_data`/`cell_data`/`field_data` are ignored entirely; only geometry and cell connectivity affect the output.

## Quirks & limitations

- No winding correction on `quad` cells — a "crossed" (bowtie) node ordering renders incorrectly with no error raised.
- Unsupported cells vanish from the output silently.
- The painter's algorithm sorts whole faces by centroid depth — mutually intersecting faces (which a closed skin never has) can stack in the wrong order; there is no per-pixel depth test.
- Write-only; there is no way to read a TikZ figure back into a `Mesh`.

## Notes

- Backed by the **C++ core** (`write_tikz`) with a pure-Python fallback: `meshioplusplus.tikz.write` uses the C++ writer for real file paths and falls back to Python for file-object/buffer targets or on any error. The C++ writer is byte-for-byte identical to the Python reference — **including the 3D projected path** (the camera arithmetic in `detail/projection.hpp` and `_projection.py` is kept expression-for-expression identical for this reason). Registered in the shared dispatch registry, so it is also reachable from the WASM, C API, and Fortran flat bindings (write-only, fixed default styling and the default isometric camera; the flat surface always emits the standalone document).
- `tests/test_tikz.py` checks the document/`tikzpicture` wrappers and `\draw` count, cross-checks the C++ and Python writers are byte-identical (2D **and** 3D), and covers the `standalone=False` snippet, volume-skin rendering, and camera angles. `cpp/tests/test_svg_tikz.cpp` covers the C++ writer directly (standalone vs snippet, filled faces vs open lines, `\draw` count, style/scale options, the projected 3D paths).

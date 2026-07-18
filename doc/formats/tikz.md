# TikZ (`.tikz`)

[TikZ/PGF](https://tikz.dev/) output for 2D meshes — a **write-only** visualization format that draws the mesh cells as a LaTeX figure. By default it emits a standalone, directly `pdflatex`-compilable document; it can also emit a bare `tikzpicture` snippet for `\input` into a larger LaTeX document. It is the LaTeX counterpart to the [SVG](./svg.md) writer.

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
)
```

- **`float_fmt`** — coordinate number format.
- **`standalone`** — if `True` (default), wrap the `tikzpicture` in a full `\documentclass{standalone}` + `\usepackage{tikz}` document that compiles directly with `pdflatex`. If `False`, emit only the `\begin{tikzpicture}…\end{tikzpicture}` environment for `\input` into an existing document.
- **`line_width`** — TikZ line width for edges, e.g. `"0.4pt"`; if `None` (default), TikZ's own default width is used. When set it is applied both on the `tikzpicture` options and on each `\draw`.
- **`fill`** — [xcolor](https://ctan.org/pkg/xcolor) fill spec for the filled faces (triangles/quads), e.g. `"gray!30"`, `"blue!20"`.
- **`draw`** — xcolor spec for the edge stroke.
- **`scale`** — optional `\begin{tikzpicture}[scale=…]` factor; if `None` (default), coordinates are emitted verbatim and no `scale` key is added.

## File structure

One `\draw` command per drawable cell inside a single `tikzpicture` environment. Each cell's vertices are emitted as `(x,y)` coordinate pairs (`float_fmt`-formatted) joined by TikZ's `--` path operator:

| cell type | `\draw` template |
|---|---|
| `line` | `\draw[draw=…] (x0,y0) -- (x1,y1);` (open, no `cycle`) |
| `triangle` | `\draw[fill=…, draw=…] (x0,y0) -- (x1,y1) -- (x2,y2) -- cycle;` |
| `quad` | `\draw[fill=…, draw=…] (x0,y0) -- (x1,y1) -- (x2,y2) -- (x3,y3) -- cycle;` |

Points must be flat 2D: if `points.shape[1] == 3`, every z coordinate must be `~0` (`atol=1e-14`), else `WriteError`. Only the first two columns are used.

Unlike the SVG writer, the y-coordinate is **not** flipped — TikZ/PGF already uses the math convention (y grows upward), so mesh coordinates map straight onto the canvas.

## Cell types

`line`, `triangle`, `quad` only. Any other cell block present in the mesh is **silently dropped** (matching the SVG writer's behaviour).

## Data mapping

None consumed — `point_data`/`cell_data`/`field_data` are ignored entirely; only geometry and cell connectivity affect the output.

## Quirks & limitations

- No winding correction on `quad` cells — a "crossed" (bowtie) node ordering renders incorrectly with no error raised.
- Non-`line`/`triangle`/`quad` cells vanish from the output silently.
- Write-only; there is no way to read a TikZ figure back into a `Mesh`.
- No C++ implementation — this format is Python-only.

## Notes

- No reference test fixture (there being no reader to test against one); `tests/test_tikz.py` asserts that writing does not raise, checks for the document/`tikzpicture` wrappers, and covers the `standalone=False` snippet and the non-flat-3D `WriteError` path.
- Implemented in pure Python (no C++ core path).

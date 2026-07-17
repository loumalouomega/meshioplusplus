# SVG (`.svg`)

[Scalable Vector Graphics](https://www.w3.org/TR/SVG/) output for 2D meshes — a **write-only** visualization format that draws the mesh edges.

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
)
```

- **`float_fmt`** — coordinate number format.
- **`stroke_width`** — edge stroke width; if `None` (default), auto-computed as 1% of the mesh's on-canvas width.
- **`image_width`** — output SVG width in user units, default `100` (not `None`) — deliberately non-trivial because some SVG viewers (e.g. `eog`) mis-render images whose natural width is close to `1` unit.
- **`fill`** / **`stroke`** — cell fill and edge colours, defaulted to match ParaView's default rendering colours (per an inline source comment).

## File structure

A single `<svg>` root containing one `<path>` element per drawable cell (**not** `<polygon>`) — chosen deliberately: the comment in the source notes that `svgo` (a common SVG optimizer) converts `<polygon>`s to `<path>`s but drops style information when it does so, so meshio++ emits paths directly to sidestep that.

Path `d` templates (space-separated coordinate pairs, `float_fmt`-formatted):

| cell type | path template |
|---|---|
| `line` | `M x0 y0L x1 y1` (open, no closing `Z`) |
| `triangle` | `M x0 y0L x1 y1L x2 y2Z` |
| `quad` | `M x0 y0L x1 y1L x2 y2L x3 y3Z` |

Points must be flat 2D: if `points.shape[1] == 3`, every z coordinate must be `~0` (`atol=1e-14`), else `WriteError`. The y-coordinate is flipped (`max_y + min_y - y`) to convert from the mesh/math convention (y-up) to SVG's screen convention (y-down).

## Cell types

`line`, `triangle`, `quad` only. Any other cell block present in the mesh is **silently dropped — no warning at all** (unlike most other meshio++ writers' warn-and-skip convention for unsupported cell types).

## Data mapping

None consumed — point_data/cell_data/field_data are ignored entirely; only geometry and cell connectivity affect the output.

## Quirks & limitations

- No diagonal/winding correction on `quad` cells — a "crossed" (non-convex, bowtie) node ordering renders incorrectly with no error raised.
- Non-`line`/`triangle`/`quad` cells vanish from the output silently.
- Write-only; there is no way to read an SVG back into a `Mesh`.
- No C++ implementation — this format is Python-only.

## Notes

- No reference test fixture (there being no reader to test against one); `tests/test_svg.py` only asserts that writing a mesh doesn't raise.
- Implemented in pure Python (no C++ core path).

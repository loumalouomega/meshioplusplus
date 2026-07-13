# SVG (`.svg`)

[Scalable Vector Graphics](https://www.w3.org/TR/SVG/) output for 2D meshes — a
**write-only** visualization format that draws the mesh edges.

| | |
|---|---|
| **Format name** | `svg` |
| **Extensions** | `.svg` |
| **Read / Write** | — / ✓ |
| **Extra dependencies** | — |

## Writing

```python
import meshio

meshio.svg.write("out.svg", mesh,
    float_fmt=".3f",
    stroke_width=None,
    image_width=100,
    fill="#c8c5bd",
    stroke="#000080",
)
```

- **`float_fmt`** — coordinate format.
- **`stroke_width`** — edge width (auto if `None`).
- **`image_width`** — output width in pixels.
- **`fill`** / **`stroke`** — cell fill and edge colours.

## Cell types

2D `triangle` / `quad` meshes.

## Notes

- Write-only, 2D only. Intended for quick visual inspection, not data exchange.
- Implemented in pure Python (no C++ core path).

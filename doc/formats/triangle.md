# Triangle (`.node` / `.ele` / `.poly`)

The file formats of Shewchuk's [Triangle](https://www.cs.cmu.edu/~quake/triangle.html) 2D mesh generator — the planar analogue of [TetGen](./tetgen.md): a `.node`/`.ele` sibling pair (points + triangles) and the `.poly` PSLG format (points + boundary segments).

| | |
|---|---|
| **Format name** | `triangle` |
| **Extensions** | `.node`, `.ele` (shared with tetgen), `.poly` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshioplusplus

mesh = meshioplusplus.read("mesh.node")                     # 2D pair -> triangle (see dispatch below)
mesh = meshioplusplus.triangle.read("mesh.node")            # explicit
meshioplusplus.triangle.write("out.node", mesh)             # writes out.node + out.ele
meshioplusplus.triangle.write("out.poly", mesh)             # writes the PSLG
meshioplusplus.write("out.node", mesh, file_format="triangle")
```

## Extension dispatch (vs tetgen)

`.node`/`.ele` are registered to **tetgen first** (3D pairs stay the default everywhere, including the WASM/C/Fortran flat bindings). Reading through the generic `meshioplusplus.read()`:

- a **3D** `.node`/`.ele` pair reads via tetgen, exactly as before;
- a **2D** pair makes tetgen raise (`Need 3D/2D points` mismatch) and the extension dispatcher falls through to `triangle` automatically — no explicit format needed.

*Writing* to `.node`/`.ele` without an explicit format still picks tetgen; pass `file_format="triangle"` (or call `meshioplusplus.triangle.write`) for 2D meshes. Only `.poly` defaults to `triangle`. The flat bindings (WASM, C API, Fortran) have a single-winner extension map, so reading a 2D `.node` there always needs `format="triangle"`.

## File structure

**`.node`**: header `npoints 2 nattrs nbmarkers`, rows `idx x y attrs.. markers..`. The node index base (0 or 1) is auto-detected from the first row and indices must be exactly consecutive (as in tetgen). A lone `.node` file (no `.ele` sibling) reads as a point cloud.

**`.ele`**: header `ntriangles 3|6 nattrs`, rows `idx n0..n(k-1) attrs..` — 3-node files become `triangle` blocks, 6-node files `triangle6`.

**`.poly`**: a vertex section (inline, or the sibling `.node` when the vertex count is 0), a segment section (`idx a b [marker]` → one `line` cell block), a hole section, and optional regional attributes. Holes and regions are skipped with a warning (they describe meshing directives, not mesh entities).

Write: a `.node`/`.ele` path writes the pair (all `triangle`/`triangle6` blocks, which must share one type; other cell types are skipped with a warning); a `.poly` path writes inline vertices, the `line` blocks as segments, and zero holes. Points must be 2D.

## Data mapping

Mirrors tetgen's naming:

- `point_data["triangle:attr{k}"]` — vertex attribute columns.
- `point_data["triangle:ref"]` / `"triangle:ref2"` / ... — boundary marker columns.
- `cell_data["triangle:ref"]` / ... — element attribute columns (`.ele`) or segment markers (`.poly`).

## Quirks & limitations

- The format spans multiple files and **cannot** be read from or written to a buffer.
- `.ele` attribute columns are read as `float64` (Triangle's regional attributes may be non-integer), unlike tetgen's integer parsing.
- Mixing `triangle` and `triangle6` blocks in one mesh raises `WriteError` (a `.ele` file has a single nodes-per-triangle header).
- `.poly` files with both segments and triangles are not a thing in Triangle itself — triangles are never written to `.poly` (use the `.node`/`.ele` pair).

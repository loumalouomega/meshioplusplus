# FreeFem++ mesh (`.msh`)

The [FreeFem++](https://freefem.org/) `.msh` mesh format (as handled by
[FEconv](https://github.com/victorsndvg/FEconv)). ASCII, with a volume-element
block and a boundary-element block, each entity carrying an integer
region/boundary label.

| | |
|---|---|
| **Format name** | `freefem` |
| **Extensions** | `.msh` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshioplusplus

# `.msh` is shared with ansys and gmsh — be explicit on read:
mesh = meshioplusplus.read("mesh.msh", file_format="freefem")
meshioplusplus.freefem.write("out.msh", mesh)
```

`write` takes no keyword arguments.

## File structure

```
nver n_el1 n_el2
x y [z] ref        # nver rows
v0 v1 ... ref      # n_el1 rows (volume element)
v0 v1 ... ref      # n_el2 rows (boundary element)
```

The header is 3 integers: vertex count, then the two element-block counts.
The spatial dimension is **inferred** from the first vertex row's token count
minus one (must resolve to 2 or 3) — there is no explicit dimension field.

- **2D**: `n_el1` rows are `triangle` (3 nodes + ref), `n_el2` rows are
  `line` (2 nodes + ref, the boundary edges).
- **3D**: `n_el1` rows are `tetra` (4 nodes + ref), `n_el2` rows are
  `triangle` (3 nodes + ref, the boundary faces).

All connectivity is 1-based. Blank lines are ignored; the reader consumes
tokens from non-blank lines only (`_nonblank_lines`).

## Cell types

Linear only, dimension-dependent as above: `{triangle, line}` for 2D meshes,
`{tetra, triangle}` for 3D meshes.

## Data mapping

- `point_data["freefem:ref"]` — the per-vertex label.
- `cell_data["freefem:ref"]` — the per-element label, one array per cell
  block; defaults to zero if not present when writing.

## Quirks & limitations

- The `.msh` extension is shared with [`ansys`](./ansys.md) and
  [`gmsh`](./gmsh.md). On extension-based auto-detection, the three formats
  are tried in registration order and `freefem` is attempted **last** — pass
  `file_format="freefem"` explicitly to select it reliably (see
  `test_explicit_file_format`).
- The writer only ever emits the two cell types appropriate for the mesh's
  dimension; any other cell type triggers a Python `warn` (and, in the C++
  writer, a hard `WriteError` that forces the Python fallback, which then
  performs the warn-and-skip).
- Points are written with full Python `repr()` precision in the Python writer
  vs. `%.16e` in the C++ writer — same effective precision, different string
  form (does not affect round-trip correctness).

## Notes

- Fully handled by the C++ core.
- No reference fixture exists under `tests/meshes/freefem/`; tests round-trip
  `tri_mesh_2d`, `tri_mesh`, and `tet_mesh`.

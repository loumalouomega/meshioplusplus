# XYZ (`.xyz`, `.xyzn`, `.xyzrgb`, `.asc`, `.pts`, `.txt`)

Headerless ASCII point clouds: one point per line, columns separated by whitespace, commas or semicolons; blank lines and `#` / `//` comments are skipped. XYZ is a convention, not a specification, so the column meaning is resolved in this order: an explicit `columns=` list; a header comment naming the columns (`# x y z nx ny nz`, what the writer emits); the column count and the file extension (3 = xyz, 4 = xyz + scalar, 6 = by extension or by value range, `.pts` with 7 = xyz + intensity + rgb); anything still ambiguous is an error asking for the list.

| | |
|---|---|
| **Format name** | `xyz` |
| **Extensions** | `.xyz`, `.xyzn`, `.xyzrgb`, `.asc`, `.pts`, `.txt` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshioplusplus

mesh = meshioplusplus.read("cloud.xyz")
meshioplusplus.xyz.write("out.xyz", mesh)
mesh = meshioplusplus.xyz.read("odd.xyz", columns=["x", "y", "z", "_", "p"])
```

- **`columns`** — one name per column: `x y z nx ny nz r g b a`, any other name is a scalar, `_` skips the column. Empty = infer.
- **`delimiter`** — column delimiter. Empty = detect (`;`, then `,`, else whitespace).
- **`float_fmt`** (write) — a format spec such as `".16e"`; by default float32 columns use `.9g` and the rest `.17g`.

## Mesh mapping

Float64 points plus one `vertex` block — what [`subsample_points`](../point_budgets.md) produces. Columns `nx ny nz` → `"normals"` (n, 3); `r g b [a]` → `"rgb"` / `"rgba"` `uint8` (byte values, or unit floats scaled by 255); every other named column is a float64 scalar under its own name.

## Quirks & limitations

- **Chemistry XYZ shares the extension** (atom count on line 1, comment on line 2, element symbols in column 1). It is detected and refused with a message naming it — meshio++ has no reader for it.
- **Six columns are ambiguous by count alone** (`xyz` + normals *or* RGB): resolved by value range (unit vectors are normals, byte-valued integers are colours), overridable by extension (`.xyzn` / `.xyzrgb`) or by `columns=`.
- **`.pts` files start with a point count** that is validated against the rows found.
- Only `vertex` cells are writable (anything else is skipped with a warning); all cell data is dropped. Points below 3D are padded with zeros.

## Notes

- The C++ core and the Python reference are byte-identical (same header comment, same float formatting).

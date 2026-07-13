# TetGen (`.node` / `.ele`)

The [TetGen](https://wias-berlin.de/software/tetgen/fformats.html) mesh format:
a pair of sibling files sharing a stem — `<name>.node` (points, attributes,
boundary markers) and `<name>.ele` (tetrahedra, region attributes).

| | |
|---|---|
| **Format name** | `tetgen` |
| **Extensions** | `.node`, `.ele` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshio

mesh = meshio.read("mesh.node")   # reads the .node/.ele pair
meshio.tetgen.write("out.node", mesh, float_fmt=".16e")
```

- **`float_fmt`** — coordinate format.

Either path (`.node` or `.ele`) selects the pair.

## File structure

`.node`: header `npoints dim nattrs nbmarkers` (dim must be 3), then
`idx x y z attrs… markers…`. `.ele`: header `ntets 4 nattrs`, then
`idx n0 n1 n2 n3 attrs…`. The node index base is auto-detected.

## Cell types

`tetra` only.

## Data mapping

- Point attributes → `point_data["tetgen:attr{k}"]`; boundary markers →
  `point_data["tetgen:ref"]` (etc.).
- Region attributes → `cell_data["tetgen:ref"]` (etc.).

## Notes

- The format spans two files and **cannot** be read from or written to a buffer.
- Fully handled by the C++ core (Python fallback for buffers / non-default
  `float_fmt`).

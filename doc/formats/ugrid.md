# AFLR UGRID (`.ugrid`)

The [AFLR UGRID](https://www.simcenter.msstate.edu/software/documentation/ug_io/3d_grid_file_type_ugrid.html)
format: a binary/ASCII surface+volume mesh whose byte layout is encoded in the
penultimate filename suffix.

| | |
|---|---|
| **Format name** | `ugrid` |
| **Extensions** | `.ugrid` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshio

mesh = meshio.read("sphere.b8.ugrid")   # flavour from the suffix
meshio.ugrid.write("out.lb8.ugrid", mesh)
```

The penultimate suffix selects the flavour: `ascii` plus binary variants
`b8l/b8/b4` (big-endian), `lb8l/lb8/lb4` (little-endian), and `r8/r4/lr8/lr4`
(Fortran-record). These encode `{big,little}` endianness, `{4,8}`-byte ints and
`{4,8}`-byte floats.

## File structure

A 7-integer header (point / triangle / quad / tetra / pyramid / wedge / hex
counts), then points, surface connectivity (1-based), surface boundary tags, and
volume elements. Fortran-record flavours wrap each record in length markers.

## Cell types

`triangle`, `quad` (surface), `tetra`, `pyramid`, `wedge`, `hexahedron` (volume),
with the pyramid node reorder.

## Data mapping

- Surface boundary tags (and zero-filled volume tags) → `cell_data["ugrid:ref"]`.

## Notes

- Fully handled by the C++ core across every flavour (host-relative byte
  swapping).

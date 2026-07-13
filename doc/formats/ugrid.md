# AFLR UGRID (`.ugrid`)

The [AFLR UGRID](https://www.simcenter.msstate.edu/software/documentation/ug_io/3d_grid_file_type_ugrid.html)
format: a binary/ASCII surface+volume mesh whose byte layout is encoded in
the **penultimate** filename suffix (e.g. `foo.lb8.ugrid` → flavour `lb8`).

| | |
|---|---|
| **Format name** | `ugrid` |
| **Extensions** | `.ugrid` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshio

mesh = meshio.read("sphere.b8.ugrid")   # flavour taken from the suffix
meshio.ugrid.write("out.lb8.ugrid", mesh)
```

No kwargs — the flavour is entirely inferred from the filename.

## File-flavour table

| suffix | kind | float dtype | int dtype |
|---|---|---|---|
| *(none)* | ascii | native | native |
| `b8l` | C, big-endian | `>f8` | `>i8` |
| `b8` | C, big-endian | `>f8` | `>i4` |
| `b4` | C, big-endian | `>f4` | `>i4` |
| `lb8l` | C, little-endian | `<f8` | `<i8` |
| `lb8` | C, little-endian | `<f8` | `<i4` |
| `lb4` | C, little-endian | `<f4` | `<i4` |
| `r8` | Fortran-record, big-endian | `>f8` | `>i4` |
| `r4` | Fortran-record, big-endian | `>f4` | `>i4` |
| `lr8` | Fortran-record, little-endian | `<f8` | `<i4` |
| `lr4` | Fortran-record, little-endian | `<f4` | `<i4` |

"C" flavours have no extra framing; "Fortran-record" flavours (`r*`/`lr*`)
wrap the whole body in Fortran unformatted-record markers: a leading integer
byte-count, the payload, then a trailing integer repeating the same
byte-count. There are exactly **2** such records total: record 1 is the
7-integer header; record 2 is everything else (points, connectivity,
boundary tags) combined.

## File structure

Header (always 7 integers): `num_points, num_triangle, num_quad, num_tetra,
num_pyramid, num_wedge, num_hexahedron`. Body, in fixed order:

1. Points: `num_points × 3` floats.
2. Triangle connectivity (1-based, `-1`'d on read) — only if
   `num_triangle > 0`.
3. Quad connectivity (1-based) — only if `num_quad > 0`.
4. Triangle boundary tags — only if `num_triangle > 0`.
5. Quad boundary tags — only if `num_quad > 0`.
6. Tetra connectivity.
7. Pyramid connectivity — **reordered** via `[1,0,3,4,2]` on read
   (`[1,0,4,2,3]` on write) — see below.
8. Wedge connectivity.
9. Hexahedron connectivity.

Volume elements (6-9) get **zero-filled** boundary tags — UGRID has no
per-volume-element tag concept, so meshio synthesizes zeros for uniformity
with the surface tags.

## Cell types

`triangle`, `quad` (surface); `tetra`, `pyramid`, `wedge`, `hexahedron`
(volume). Node ordering matches meshio's convention for every type **except
pyramids**, which need the permutation above — this is the one place in the
format where getting the order wrong would silently produce inverted-volume
elements rather than an outright error, so it's specifically covered by a
signed-volume regression test.

## Data mapping

- `cell_data["ugrid:ref"]` — one array per written cell block, in the fixed
  type order (triangle, quad get real boundary tags; tetra/pyramid/wedge/
  hexahedron get all-zero tags).

## Quirks & limitations

- The writer enforces **at most one cell block per known type**
  (`ValueError` otherwise) and skips unknown types with a warning.
- Fortran-record byte-count markers are **written but not validated on
  re-read** — a corrupted record-length byte wouldn't be caught by the
  reader.
- Boundary tags on write come from the first integer `cell_data` array found
  (warns if more than one candidate exists); defaults to all-`1` if none.
- 4-byte-float variants (`b4`/`lb4`/`r4`/`lr4`) inherently have lower
  round-trip precision than the 8-byte variants — reflected in the test
  suite's looser tolerance for those flavours.

## Notes

- `tests/meshes/ugrid/pyra_cube.ugrid` (ascii) — a unit cube split into 6
  pyramids; used to verify the pyramid permutation via a signed-volume sum
  (≈1.0).
- `tests/meshes/ugrid/sphere_mixed.1.lb8.ugrid` (little-endian binary) —
  3270 points, 864 triangles, 3024 wedges, 9072 tets, boundary tag counts
  `{1:432, 2:216, 3:216}`.
- `tests/meshes/ugrid/hch_strct.4.lb8.ugrid` — 306 points, 12 tri, 178 quad,
  96 wedge, 144 hex; also used to check total surface area, catching any
  connectivity-order regression.
- Fully handled by the C++ core across every flavour (host-relative byte
  swapping).

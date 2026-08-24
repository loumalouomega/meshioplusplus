# VTI — VTK XML ImageData (`.vti`)

A regular lattice whose geometry is **three attributes** rather than a point array. It is the natural home for a generated grid — `grid`, `voxelize` with `fill="all"`, and `compute_sdf` with `structure="voxel"` all produce exactly the shape this format stores — and the only format in meshio++ that round-trips such a grid's header. See the [VTK XML file formats documentation](https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html).

| | |
|---|---|
| **Format name** | `vti` |
| **Extensions** | `.vti` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — (zlib for compressed output, as VTU) |

## Reading & writing

```python
import meshioplusplus

field = meshioplusplus.compute_sdf(skin, resolution=(128, 128, 128))
meshioplusplus.vti.write("field.vti", field,
    binary=True,            # base64-encode the arrays instead of writing text
    compression="zlib",     # "zlib", "lz4", "zstd", or None
    header_type=None,       # "UInt32" (default) or "UInt64"
)

back = meshioplusplus.vti.read("field.vti")
```

`meshioplusplus.read`/`write` dispatch on the extension, so `.vti` needs no explicit format name.

## File structure

```xml
<?xml version="1.0"?>
<VTKFile type="ImageData" version="0.1" byte_order="LittleEndian"
         compressor="vtkZLibDataCompressor">
<ImageData WholeExtent="0 16 0 16 0 16"
           Origin="-0.17320508075688773 -0.17320508075688773 -0.17320508075688773"
           Spacing="0.084150635094610967 0.084150635094610967 0.084150635094610967">
<Piece Extent="0 16 0 16 0 16">
<PointData>
  <DataArray type="Float64" Name="sdf:distance" format="binary">…</DataArray>
</PointData>
<CellData/>
</Piece>
</ImageData>
</VTKFile>
```

There is **no `<Points>` section**: `Origin`, `Spacing` and `WholeExtent` *are* the geometry. `WholeExtent` counts **points**, so `0 16` on an axis means 16 cells and 17 planes.

## The mesh side of the deal

A `Mesh` has no implicit geometry, so the two directions are not symmetric:

- **`read` expands** the extent into explicit points and one `hexahedron` block, through the same numbering [`grid`](/voxelize) produces (x fastest, then y, then z). That is what makes `read(write(m)) == m` an identity rather than a coincidence.
- **`write` requires a lattice.** A mesh that is not one has no extent to write and raises `WriteError` by name.

::: warning A *partial* grid cannot be written either
`voxelize`'s `surface`/`inside` fills and `compute_sdf`'s octree produce meshes made of the right cells in the right places — and ImageData still has no way to express their holes. Writing one would silently fill them in, so it is refused. Write those as `.vtu`, which stores the cells explicitly.
:::

The recovery is exact rather than a fit: `grid` writes every coordinate as `origin + index * spacing` evaluated independently, so all points on one plane carry bit-identical values and an exact sort-and-unique recovers the planes with no tolerance. Only the *spacing* needs one — consecutive gaps of `origin + i*h` differ in the last bits — so uniformity is checked to a relative `1e-9`, which rejects a genuinely graded mesh by orders of magnitude rather than by ulps.

## Data mapping

| meshio++ | VTI |
|---|---|
| `point_data` | `<PointData>`, one `<DataArray>` per array, sorted by name |
| `cell_data` | `<CellData>` — a lattice has one block, so one array each |
| `field_data` | **not written** (see below) |
| `point_sets` / `cell_sets` | not written |

Both `Origin` and `Spacing` are written with `%.17g` — the round-trip width for a double. The stream default of six significant digits would lose about ten digits of a real origin, placing the grid ~1e-7 off its own points with nothing downstream to flag it.

## Why this format exists for `compute_sdf`

A generated grid carries an `sdf:*` `field_data` header describing itself. **No file format persists arbitrary `field_data`** — not gmsh, not MED, not VTU — so a grid written anywhere and read back has lost it. `.vti` does not need to persist it: its three attributes are the same information, and reading them back reconstructs the identical mesh.

## Quirks & limitations

- **`<AppendedData>` is not supported** in either implementation, and raises. The VTU C++ reader declines it too, for the same reason.
- **One `<Piece>` only**, and its `Extent` must equal the `WholeExtent`. A partial piece's arrays are sized to the *piece*, so reading them against the whole extent would be silently misaligned.
- **A non-identity `Direction`** (a rotated lattice) raises: an axis-aligned hexahedron grid cannot express it without baking the rotation into the coordinates, which is a different mesh from the one the file describes.
- **lzma is rejected** by both readers. Python has the module; declining it is a deliberate parity choice, so that the two readers accept the same files.
- `header_type="UInt64"` is honoured on read; the writer always emits the default `UInt32`, as the VTU writer does.
- The writer always declares `LittleEndian`; the readers honour `byte_order`.

## Notes

The C++ reader/writer share `detail/vtk_xml.hpp`'s `<DataArray>` codec and `detail/vtu_binary.hpp`'s base64 + 32 KiB block framing with VTU and VTP verbatim, so `--codec zlib|lz4|zstd` works here exactly as it does there.

Testing does **not** rely on the round trip through this library's own writer, which a consistently wrong `Origin` or a transposed extent would survive unchanged. `tests/cpp/test_vti.cpp` and `tests/python/test_vti.py` additionally assert on the raw written bytes and read hand-written files whose expected geometry is stated independently — a non-zero `Origin`, an extent not starting at zero, an inverted extent, a foreign `Piece`.

## See also

- [Signed distance fields](/sdf) — what usually produces a `.vti`.
- [Regular grids](/voxelize) — the numbering this format shares.
- [VTU](/formats/vtu) — the explicit-cell sibling, for a partial grid.

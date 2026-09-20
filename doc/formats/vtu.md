# VTU — VTK XML UnstructuredGrid (`.vtu`)

The [serial VTK XML](https://vtk.org/Wiki/VTK_XML_Formats) UnstructuredGrid format: an XML container whose `DataArray` payloads can be inline ASCII, inline base64 binary, or appended raw/base64 binary, optionally block-compressed.

| | |
|---|---|
| **Format name** | `vtu` |
| **Extensions** | `.vtu` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshioplusplus

mesh = meshioplusplus.read("mesh.vtu")
meshioplusplus.vtu.write("out.vtu", mesh,
    binary=True,
    compression="zlib",  # "zlib", "lzma", or None
    header_type=None,    # "UInt32" or "UInt64"
)
```

- **`binary`** — base64-encoded binary DataArrays (`True`) or ASCII.
- **`compression`** — block compression filter for binary data (`vtkZLibDataCompressor`/`vtkLZMADataCompressor`, or none).
- **`header_type`** — integer type used for the binary block header/sizes (default `UInt32`).

## File structure

```xml
<VTKFile type="UnstructuredGrid" byte_order="..." [compressor="..."] [header_type="..."]>
  <UnstructuredGrid>
    <Piece NumberOfPoints=".." NumberOfCells="..">
      <Points><DataArray .../></Points>
      <Cells>
        <DataArray Name="connectivity" .../>
        <DataArray Name="offsets" .../>
        <DataArray Name="types" .../>
        <!-- polyhedron only: -->
        <DataArray Name="faces" .../><DataArray Name="faceoffsets" .../>
      </Cells>
      <PointData><DataArray Name=".." .../></PointData>
      <CellData><DataArray Name=".." .../></CellData>
    </Piece>
    <FieldData><DataArray Name=".." .../></FieldData>
  </UnstructuredGrid>
  <AppendedData encoding="base64">_<blob></AppendedData>
</VTKFile>
```

**Binary encoding scheme** (matching VTK's own convention exactly, so files round-trip byte-for-byte with other VTK tools):

- Uncompressed: `base64(header[header_type: total_nbytes] + raw_bytes)`.
- Compressed: `base64(header[nblocks, blocksize=32768, last_block_size, csize_0..csize_{n-1}])`, followed by a **separate** base64 blob of `concat(compressed_block_0..n-1)`. Header fields use the file's declared `header_type` dtype throughout.

## Cell types

The full VTK cell set, including VTK Lagrange high-order cells (`VTK_LAGRANGE_*`). See [VTK](./vtk.md) for the shared numeric type-code table.

## Data mapping

`<PointData>`/`<CellData>` map generically to `point_data`/`cell_data`; `cell_sets` round-trip as extra data arrays with an info-level message (VTU has no native set concept).

`<FieldData>` ↔ `mesh.field_data` (v15.0.0): mesh-level arrays travel in a `<FieldData>` element on the grid, before the `<Piece>`, exactly where VTK's own writers put it, one `<DataArray>` per name with the `NumberOfTuples` VTK requires (and `NumberOfComponents` for a two-or-more-dimensional array; a higher rank is flattened to `(tuples, components)`, and a rank-0 scalar reads back as a length-1 array). Reading accepts it on the grid and inside a `<Piece>` (the piece's overriding the grid's); an array of a non-numeric type (`type="String"`) has no meshio++ dtype and is skipped with a warning rather than failing the read. A mesh without field data writes no `<FieldData>` element at all, so every existing file is byte-identical. A `TimeValue` array is VTK's "time in field data" convention: ParaView reads it as the dataset's time step (verified against ParaView 6.1.1), and a [`.pvd`](./pvd.md) entry with no `timestep=` takes its step time from it. A value that is not a numeric array has no VTK type and is not written (Python warns).

## Quirks & limitations

- **Raw/appended binary without valid XML**: when appended binary data contains raw bytes that break XML parsing, the Python reader falls back to a regex-based manual split of the file into header/data/footer before continuing — **the C++ reader does not implement this path at all** and raises on any `<AppendedData>` section, forcing the Python fallback.
- **lzma compression is Python-only** — the C++ reader/writer explicitly reject it.
- **Polyhedron cells are read and written by both engines**, including meshes that mix polyhedra with other cell types in one file (the C++ reader and writer carry VTU's `faces`/`faceoffsets` streams, with `-1` marking a non-polyhedral cell; the Python writer has accepted mixing since v9.19.0), as an OpenFOAM-derived mesh always requires.
- **Multi-`<Piece>` files**: the Python reader merges all pieces (concatenating points/cells/point_data across them); the C++ reader only supports a single `<Piece>` and throws otherwise.
- A `header_type` other than the default (`None`, meaning `UInt32`) always forces the Python path.
- 2D points are auto-padded to 3D on write (warning in Python; silent in C++).
- Byte order: the Python writer records the system's native byte order in the `byte_order` attribute; the C++ writer always declares `LittleEndian`.

## Notes

- `tests/python/meshes/vtu/00_raw_binary.vtu`, `01_raw_binary_int64.vtu` (uses an Int64 header type), `02_raw_compressed.vtu` (zlib-compressed appended data) — each a 162-point, 64-cell `tetra` mesh, all exercised via the raw-binary fallback path described above.
- The C++ core handles ASCII, uncompressed binary, and **zlib** binary (when built with `MESHIO_WITH_ZLIB`; otherwise the Python stdlib handles zlib too) — see [native acceleration](../formats.md#native-acceleration-and-fallbacks).

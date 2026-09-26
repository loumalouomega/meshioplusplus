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
    appended=False,      # True: raw binary in one <AppendedData> section
)
```

- **`binary`** — base64-encoded binary DataArrays (`True`) or ASCII.
- **`compression`** — block compression filter for binary data (`vtkZLibDataCompressor`/`vtkLZMADataCompressor`, or none).
- **`header_type`** — integer type used for the binary block header/sizes (default `UInt32`).
- **`appended`** (since v16.21.0) — write every array as `format="appended"` into one `<AppendedData encoding="raw">` section at the end of the file: the same header and (compressed) blocks as `binary`, as raw bytes rather than base64, so the file is about a quarter smaller and is read without a base64 decode. VTK's own writers default to this layout. It needs `binary=True`, and either engine writes it (the C++ core for `None`/`zlib`/`lz4`/`zstd`, the Python reference for `lzma` too). The same encoding is `WriteEncoding::RawAppended` in C++ (`write_vtu_appended`), `MIO_ENCODING_RAW_APPENDED` in C, `"raw_appended"` for the [pipeline](../pipeline.md), WASM and the MCP `convert` tool's `mode`, and `--appended` in both CLIs; every other format refuses it by name.

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

- **Appended data** (since v16.6.0 in both engines): `<AppendedData encoding="raw">` is cut out of the file as bytes before the XML is parsed (a raw payload may hold any byte, `<` included), and every `format="appended"` array is read at its `offset` into it; `encoding="base64"` payloads are decoded the same way, whether each array's header was encoded together with its body or apart from it (VTK's appended writer pads them separately). Offsets padded with blanks (`offset="    1660"`, as FEconv writes them) are accepted. Before v16.6.0 the C++ reader refused any `<AppendedData>` and the Python fallback failed on padded offsets and on base64 appended data.
- **Byte order**: `byte_order="BigEndian"` is honoured for headers and array bodies, inline or appended, in both engines; arrays are returned in native order.
- **lzma compression is Python-only** — the C++ reader/writer explicitly reject it.
- **Polyhedron cells are read and written by both engines**, including meshes that mix polyhedra with other cell types in one file (the C++ reader and writer carry VTU's `faces`/`faceoffsets` streams, with `-1` marking a non-polyhedral cell; the Python writer has accepted mixing since v9.19.0), as an OpenFOAM-derived mesh always requires. Before v16.6.0 the Python writer (which writes every polyhedral file: the C++ writer is not reached for ragged blocks) put each polyhedron's `faceoffsets` entry one short when an ordinary cell preceded it, and VTK 9.4+ (ParaView 6) refused such files.
- **Multi-`<Piece>` files** are merged by both engines (since v16.6.0; the C++ reader used to refuse them and the Python one kept only the last piece's cells): points are concatenated, each piece's node ids, offsets and polyhedron face streams shifted accordingly, and point and cell data kept when every piece has the array; an array missing from some piece is dropped with a warning. A piece that declares points but has no `<Points>` is a `ReadError`.
- A `header_type` other than the default (`None`, meaning `UInt32`) always forces the Python path.
- 2D points are auto-padded to 3D on write (warning in Python; silent in C++).
- Byte order: the Python writer records the system's native byte order in the `byte_order` attribute; the C++ writer always declares `LittleEndian`.

## Notes

- `tests/python/meshes/vtu/00_raw_binary.vtu`, `01_raw_binary_int64.vtu` (uses an Int64 header type), `02_raw_compressed.vtu` (zlib-compressed appended data) — each a 162-point, 64-cell `tetra` mesh, all read through the raw appended path described above. `raw_bigendian.vtu`, `raw_zlib.vtu`, `base64_appended.vtu`, `binary_bigendian.vtu` and `two_pieces_polyhedron.vtu` (generated by `tools/gen_feconv_quirk_fixtures.py`) reproduce the framings of FEconv's samples: the same pyramid-on-a-polyhedral-cube in each, read identically by both engines.
- The C++ core handles ASCII, uncompressed binary, and **zlib** binary (when built with `MESHIO_WITH_ZLIB`; otherwise the Python stdlib handles zlib too) — see [native acceleration](../formats.md#native-acceleration-and-fallbacks).

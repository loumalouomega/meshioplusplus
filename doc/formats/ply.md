# PLY (`.ply`)

The [Polygon File Format](https://en.wikipedia.org/wiki/PLY_(file_format)) (aka Stanford triangle format): a header describing named element groups with typed properties, followed by ASCII or (little/big-endian) binary data.

| | |
|---|---|
| **Format name** | `ply` |
| **Extensions** | `.ply` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshioplusplus

mesh = meshioplusplus.read("bunny.ply")
meshioplusplus.ply.write("out.ply", mesh, binary=True)
```

- **`binary`** — binary (`True`, the default) or ASCII.
- **`skin`** — when `True` (the default) and the mesh contains supported 3D volume cells (tetra, hexahedron, wedge, pyramid and their common higher-order variants), the boundary skin is extracted first (see [`extract_skin`](../extract_skin.md)) and only that skin is written as triangles/quads, with the vertex table compacted to the referenced points. Pre-existing surface blocks are dropped with a warning. `skin=False` restores the legacy behavior (volume cells are skipped).

## File structure

```
ply
format ascii 1.0 | format binary_little_endian 1.0 | format binary_big_endian 1.0
comment ...
element vertex <N>
property <type> <name>          # repeated, x/y/[z] conventionally first
[element face <M>
 property list <count_type> <index_type> vertex_indices
 [property <type> <name> ...]]  # extra per-face scalar properties (Python reader only)
end_header
<vertex data, N rows>
<face data, M rows: count + indices [+ extra properties]>
```

Endianness is read directly from the `format` line (`ascii`/ `binary_little_endian`/`binary_big_endian`), unlike VTU/VTK which use a separate byte-order attribute.

## Cell types

Faces are grouped by vertex count: 1→`vertex`, 2→`line`, 3→`triangle`, 4→`quad`, else→`polygon`.

## Data mapping

- `point_data[<name>]` — any vertex property beyond `x`/`y`/`z` (e.g. `nx,ny,nz` normals, `confidence`, `intensity`).
- `cell_data[<name>]` — any face property beyond `vertex_indices` (**Python reader only** — see limitations).

## Quirks & limitations

- **Ragged face-list parsing in binary mode**: since each face row's length is only known by reading its own leading count, the reader first walks the whole buffer computing per-row byte offsets, then groups **consecutive constant-length runs** into separate cell blocks — this differs from most other formats, where all same-typed cells end up in one block regardless of their position in the file; here, position matters.
- **Extra face properties beyond the index list are Python-only** — the C++ reader explicitly rejects any face `property` beyond the single index list, as well as list-typed *vertex* properties; both force a Python fallback.
- 64-bit integer cell data is **silently downcast to int32** on write (PLY has no 64-bit integer property type), with a warning; the writer also requires all cell blocks share one dtype.
- Multi-dimensional point-data is not writable (skipped per-key with a warning in Python; silently filtered in C++).
- Only `vertex`/`line`/`triangle`/`quad`/`polygon` are writable cell types. 3D volume cells write their extracted boundary skin by default (`skin=True`); with `skin=False` (or for other unsupported types) they are skipped — the Python reference writer warns per skipped block, while the default C++ path skips silently.
- `obj_info` header lines (a Meshlab convention) are skipped without being parsed.
- The uchar-property naming: `uchar` is mapped to a signed 1-byte type in the *count-field* parsing path (matching a common "uchar-as-count" convention some tools use) even though `uchar` normally means unsigned — this is deliberate, not an oversight, and only affects the binary list-count field, not general vertex/face data typing.

## Notes

- `tests/meshes/ply/bun_zipper_res4.ply` — the decimated Stanford bunny, 948 triangle cells, point-sum ≈34.14584.
- `tests/meshes/ply/tet.ply` — a tiny 4-cell tetrahedron surface, point-sum 6 (exercises extra face properties, forcing the Python path).
- The C++ core handles ASCII and both binary endiannesses for the common case (index-only faces, non-list vertex properties).

# VTS — VTK XML StructuredGrid (`.vts`)

A lattice with explicit points but *implicit* connectivity: `nx * ny * nz` hexahedra over a `WholeExtent` corner grid, the same topology `.vti` states as three attributes — except `.vts` writes the points out, exactly as `.vtu` does, so a structured mesh's points need not sit on an even grid at all (v11.6.0, roadmap [§1](../roadmap.md#_1-wasm-parity) tier B4). See the [VTK XML file formats documentation](https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html).

| | |
|---|---|
| **Format name** | `vts` |
| **Extensions** | `.vts` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — (zlib for compressed output, as VTU) |

## Reading & writing

```python
import meshioplusplus

mesh = meshioplusplus.grid([16, 16, 16])
meshioplusplus.vts.write("field.vts", mesh,
    binary=True,            # base64-encode the arrays instead of writing text
    compression="zlib",     # "zlib", "lz4", "zstd", or None
    header_type=None,       # "UInt32" (default) or "UInt64"
)

back = meshioplusplus.vts.read("field.vts")
```

`meshioplusplus.read`/`write` dispatch on the extension, so `.vts` needs no explicit format name.

## File structure

```xml
<?xml version="1.0"?>
<VTKFile type="StructuredGrid" version="0.1" byte_order="LittleEndian">
<StructuredGrid WholeExtent="0 2 0 2 0 2">
<Piece Extent="0 2 0 2 0 2">
<Points>
  <DataArray type="Float64" Name="Points" NumberOfComponents="3" format="ascii">…</DataArray>
</Points>
<PointData/>
<CellData/>
</Piece>
</StructuredGrid>
</VTKFile>
```

`WholeExtent` counts **points**, exactly as `.vti`'s does: `0 2` on an axis means 2 cells and 3 planes. Unlike `.vti`, the points themselves are an explicit `<DataArray>` — connectivity is still never written, only implied by extent + index order (x fastest, then y, then z, the same numbering [`grid`](/voxelize) produces).

## The mesh side of the deal

The two directions are not symmetric, and not for the same reason `.vti`'s aren't:

- **`read` never requires a lattice.** It reads `<Points>` verbatim and builds `hexahedron` connectivity from `WholeExtent` alone via the index formula — a genuinely curved structured mesh (VTK allows one) reads correctly, since nothing here checks that the points tile a regular grid.
- **`write` does require one**, exactly as `.vti`'s writer does (`detail::lattice_from_mesh`/`meshioplusplus._grid.lattice_from_mesh`): the writer has no way to recover which `(nx, ny, nz)` a mesh's points were meant to tile without one. Once confirmed, the mesh's own points are written **unchanged** — there is nothing to recompute, unlike `.vti`'s `Origin`/`Spacing` attributes.

::: warning A *partial* grid cannot be written either
Same as `.vti`: `voxelize`'s `surface`/`inside` fills and `compute_sdf`'s octree are the right cells in the right places, but `lattice_from_mesh` still cannot recover an extent from them. Write those as `.vtu`, which stores the cells explicitly.
:::

## Data mapping

| meshio++ | VTS |
|---|---|
| `points` | `<Points>`, one `Float64` `<DataArray>` |
| `point_data` | `<PointData>`, one `<DataArray>` per array, sorted by name |
| `cell_data` | `<CellData>` — a lattice has one block, so one array each |
| `field_data` | **not written** |
| `point_sets` / `cell_sets` | not written |

## Quirks & limitations

Identical to [`.vti`'s](./vti.md#quirks-limitations): `<AppendedData>` is declined by both readers; one `<Piece>` only, whose `Extent` must equal the `WholeExtent`; lzma is rejected by both readers (a deliberate parity choice, not a capability gap — Python has the module); `header_type="UInt64"` is honoured on read, the writer always emits the default `UInt32`.

**Degenerate (2-D/1-D) extents are not expanded to quad/line/vertex cells.** VTK itself allows a `WholeExtent` with one or more axes at zero cells; this reader, like `.vti`'s, only ever emits `hexahedron` connectivity, so a degenerate extent reads back as a points-only mesh. Documented, not silent — a genuinely 2-D/1-D structured file is a follow-up, not implemented here.

## Notes

Reuses the same `detail/vtk_xml.hpp`/`detail/vtu_binary.hpp` `<DataArray>` codec machinery `.vti`/`.vtu` do — `--codec zlib|lz4|zstd` works here exactly as it does there. The pure-Python fallback (`meshioplusplus.vts._vts`) reuses `meshioplusplus._grid._lattice_py` for the hexahedron index formula and `VtuReader`'s bound decode methods, the same reuse pattern `.vti`'s Python twin uses.

Testing does not rely on the round trip through this library's own writer alone: `tests/cpp/test_vts.cpp`/`tests/python/test_vts.py` also read a hand-written file whose points are deliberately *not* on a regular grid — the identity this format exists for over `.vti` — and assert the C++ and Python engines produce cross-readable bytes.

## See also

- [VTI](./vti.md) — the implicit-geometry sibling, for a mesh that *is* a uniform lattice.
- [VTU](./vtu.md) — the fully-explicit-cell sibling, for a partial grid or non-hexahedron topology.
- [Regular grids](/voxelize) — the numbering this format shares.

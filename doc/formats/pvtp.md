# PVTP — VTK XML parallel polygonal data (`.pvtp`)

The [`.pvtu`](./pvtu.md) index over `.vtp` pieces: it declares the arrays every piece holds and names one polygonal-data piece per part, for surface meshes, line sets and point sets that are split across ranks (v15.0.0). Everything on the [`.pvtu` page](./pvtu.md) applies unchanged — carving by `partition:part`, the declaration check, ghost cells, merging with one region per piece, `piece=` — with `.vtp` in place of `.vtu`. See the [VTK XML file formats documentation](https://docs.vtk.org/en/latest/vtk_file_formats/vtkxml_file_format.html).

| | |
|---|---|
| **Format name** | `pvtp` |
| **Extensions** | `.pvtp` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — (zlib for compressed pieces, as VTP) |

## Reading & writing

```python
import meshioplusplus

surface = meshioplusplus.read("wing.vtp")
surface.cell_data["partition:part"] = meshioplusplus.partition_labels(surface, 4)
meshioplusplus.pvtp.write("wing.pvtp", surface)       # wing.pvtp + wing/wing_0000.vtp ... wing_0003.vtp

merged = meshioplusplus.pvtp.read("wing.pvtp")
one = meshioplusplus.pvtp.read("wing.pvtp", piece=2)
meshioplusplus.pvtp.write_pieces("wing.pvtp", meshioplusplus.partition(surface, 4, ghost_layers=1))
```

## File structure

```xml
<?xml version="1.0"?>
<VTKFile type="PPolyData" version="1.0" byte_order="LittleEndian">
<PPolyData GhostLevel="0">
<PPointData>
<PDataArray type="Float64" Name="u"/>
</PPointData>
<PPoints>
<PDataArray type="Float64" Name="Points" NumberOfComponents="3"/>
</PPoints>
<Piece Source="wing/wing_0000.vtp"/>
<Piece Source="wing/wing_0001.vtp"/>
</PPolyData>
</VTKFile>
```

## Quirks & limitations

The pieces are `.vtp` files, so they carry the cell types PolyData can: vertices, lines, triangles, quads and polygons — see [VTP](./vtp.md#quirks-limitations). Every other limitation is the [`.pvtu`](./pvtu.md#quirks-limitations) one.

## See also

- [`.pvtu`](./pvtu.md) — the full description of the index, its ghost cells and its declarations.
- [VTP](./vtp.md) — what every piece actually is.
- [partition](/partition) — the operation whose output this is.

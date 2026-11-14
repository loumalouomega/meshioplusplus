# PVTI — VTK XML parallel image-data index (`.pvti`)

The parallel index over `.vti` pieces: it declares the whole grid's `WholeExtent` and names one structured piece per part, each with the `Extent` it covers inside it (read only). Reading places every piece by its `Extent`, so the points of a plane two pieces share become one layer of points (nothing is welded: the placement is the weld) and the cells are the whole lattice's, with one `piece_<i>` cell region per piece, or one piece with `piece=`. The pieces must cover `WholeExtent`; a missing piece file, a piece whose point count disagrees with its `Extent`, or a gap is a `ReadError`. See the [VTK XML file formats documentation](https://docs.vtk.org/en/latest/vtk_file_formats/vtkxml_file_format.html) and the [`.pvtu` page](./pvtu.md) for the unstructured counterpart.

| | |
|---|---|
| **Format name** | `pvti` |
| **Extensions** | `.pvti` |
| **Read / Write** | ✓ / — (see [conformance](./conformance.md#pvti)) |
| **Extra dependencies** | — (zlib for compressed pieces, as `.vti`) |

## Reading

```python
import meshioplusplus

mesh = meshioplusplus.read("grid.pvti")
one = meshioplusplus.pvti.read("grid.pvti", piece=2)
```

Python only: the native core does not read it, and the native CLI does not know the extension.
To write one, write the whole lattice as a `.vti`; there is no writer (a general mesh is not a lattice, so the write side cannot meet the round-trip conformance contract).

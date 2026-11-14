# PVTS — VTK XML parallel structured-grid index (`.pvts`)

The parallel index over `.vts` pieces: it declares the whole grid's `WholeExtent` and names one structured piece per part, each with the `Extent` it covers inside it (read only). Reading places every piece by its `Extent`, so the points of a plane two pieces share become one layer of points (nothing is welded: the placement is the weld) and the cells are the whole lattice's, with one `piece_<i>` cell region per piece, or one piece with `piece=`. The pieces must cover `WholeExtent`; a missing piece file, a piece whose point count disagrees with its `Extent`, or a gap is a `ReadError`. See the [VTK XML file formats documentation](https://docs.vtk.org/en/latest/vtk_file_formats/vtkxml_file_format.html) and the [`.pvtu` page](./pvtu.md) for the unstructured counterpart.

| | |
|---|---|
| **Format name** | `pvts` |
| **Extensions** | `.pvts` |
| **Read / Write** | ✓ / — (see [conformance](./conformance.md#pvts)) |
| **Extra dependencies** | — (zlib for compressed pieces, as `.vts`) |

## Reading

```python
import meshioplusplus

mesh = meshioplusplus.read("grid.pvts")
one = meshioplusplus.pvts.read("grid.pvts", piece=2)
```

Python only: the native core does not read it, and the native CLI does not know the extension.
To write one, write the whole lattice as a `.vts`; there is no writer (a general mesh is not a lattice, so the write side cannot meet the round-trip conformance contract).

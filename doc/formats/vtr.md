# VTR — VTK XML RectilinearGrid (`.vtr`)

A lattice whose per-axis point coordinates are three independent, only *monotonic* 1-D arrays — the genuinely more general sibling of `.vti`'s uniform `Origin`/`Spacing`: a graded grid, finer near a wall and coarser far from it, is a RectilinearGrid, never an ImageData (v11.6.0, roadmap [§1](../roadmap.md#_1-wasm-parity) tier B4). See the [VTK XML file formats documentation](https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html).

| | |
|---|---|
| **Format name** | `vtr` |
| **Extensions** | `.vtr` |
| **Read / Write** | ✓ / ✓ (write needs a **uniform** lattice — see below) |
| **Extra dependencies** | — (zlib for compressed output, as VTU) |

## Reading & writing

```python
import meshioplusplus

mesh = meshioplusplus.grid([16, 16, 16])
meshioplusplus.vtr.write("field.vtr", mesh,
    binary=True,            # base64-encode the arrays instead of writing text
    compression="zlib",     # "zlib", "lz4", "zstd", or None
    header_type=None,       # "UInt32" (default) or "UInt64"
)

back = meshioplusplus.vtr.read("field.vtr")
```

`meshioplusplus.read`/`write` dispatch on the extension, so `.vtr` needs no explicit format name.

## File structure

```xml
<?xml version="1.0"?>
<VTKFile type="RectilinearGrid" version="0.1" byte_order="LittleEndian">
<RectilinearGrid WholeExtent="0 2 0 2 0 2">
<Piece Extent="0 2 0 2 0 2">
<Coordinates>
  <DataArray type="Float64" Name="x_coordinates" format="ascii">0 1 4</DataArray>
  <DataArray type="Float64" Name="y_coordinates" format="ascii">0 1 2</DataArray>
  <DataArray type="Float64" Name="z_coordinates" format="ascii">0 1 2</DataArray>
</Coordinates>
<PointData/>
<CellData/>
</Piece>
</RectilinearGrid>
</VTKFile>
```

`WholeExtent` counts **points**, as `.vti`'s/`.vts`'s do. Point `(i, j, k)` sits at `(x_coordinates[i], y_coordinates[j], z_coordinates[k])` — a **tensor product**, not three independent offsets — so the coordinate arrays alone (never explicit `<Points>`) fully determine the geometry.

## The mesh side of the deal

The two directions are asymmetric, and more so than `.vti`'s or `.vts`'s:

- **`read` is fully general.** It builds points from the tensor product of the file's own three coordinate arrays, with **no uniformity check at all** — a genuinely graded grid reads correctly, and nothing here would even notice if it were not.
- **`write` requires a *uniform* lattice**, exactly as `.vti`'s writer does (`detail::lattice_from_mesh`/`meshioplusplus._grid.lattice_from_mesh`).

::: warning A genuinely graded mesh cannot be written as `.vtr` today
Recovering three arbitrary per-axis coordinate arrays from an unstructured point set needs a detector `lattice_from_mesh` does not implement — it only ever recovers one `Origin`/`Spacing` pair. A mesh whose axis spacing is not uniform (and is therefore exactly the kind of geometry `.vtr` exists to express) raises `WriteError`, the same message a non-lattice mesh gets. This is a **documented follow-up**, not a silent gap: see `doc/roadmap.md`. A *partial* grid (`voxelize`'s `surface`/`inside` fills, or an octree) cannot be written either, for the same reason `.vti`'s can't.
:::

## Data mapping

| meshio++ | VTR |
|---|---|
| `points` | `<Coordinates>`, three `Float64` `<DataArray>`s (`x_coordinates`/`y_coordinates`/`z_coordinates`), never written back individually |
| `point_data` | `<PointData>`, one `<DataArray>` per array, sorted by name |
| `cell_data` | `<CellData>` — a lattice has one block, so one array each |
| `field_data` | **not written** |
| `point_sets` / `cell_sets` | not written |

## Quirks & limitations

Identical to [`.vti`'s](./vti.md#quirks-limitations): `<AppendedData>` is declined by both readers; one `<Piece>` only, whose `Extent` must equal the `WholeExtent`; lzma is rejected by both readers; `header_type="UInt64"` is honoured on read, the writer always emits the default `UInt32`.

**Degenerate (2-D/1-D) extents are not expanded to quad/line/vertex cells**, the same remainder `.vts` has — this reader, like `.vti`'s and `.vts`'s, only ever emits `hexahedron` connectivity.

**The writer's uniform-lattice restriction is the format's main current limitation** — see the warning above.

## Notes

Reuses the same `detail/vtk_xml.hpp`/`detail/vtu_binary.hpp` `<DataArray>` codec machinery `.vti`/`.vts`/`.vtu` do. The pure-Python fallback (`meshioplusplus.vtr._vtr`) reuses `VtuReader`'s bound decode methods, the same pattern `.vti`'s and `.vts`'s Python twins use; its own point-construction is a `numpy.meshgrid` tensor product rather than a call into `meshioplusplus._grid`, since that module's own lattice builder assumes uniform spacing.

Testing exercises a hand-written, genuinely non-uniform fixture (`tests/cpp/test_vtr.cpp`/`tests/python/test_vtr.py`) — the identity this format exists for over `.vti`/`.vts` — in addition to the usual round-trip-through-this-library's-own-writer coverage.

## See also

- [VTI](./vti.md) — the uniform-spacing sibling, simpler and smaller on disk when it applies.
- [VTS](./vts.md) — the explicit-points sibling, for a mesh that is curved rather than merely graded.
- [VTU](./vtu.md) — the fully-explicit-cell sibling, for a partial grid or non-hexahedron topology.

# VTM — VTK XML MultiBlock (`.vtm`)

An index file plus one `.vtu` piece per cell block: unlike `.vti`/`.vts`/`.vtr`, `.vtm` carries no geometry of its own at all, and there is no lattice restriction — any mesh with one or more cell blocks round-trips (v11.6.0, tier B4, part 3 of 3). See the [VTK XML file formats documentation](https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html).

| | |
|---|---|
| **Format name** | `vtm` |
| **Extensions** | `.vtm` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — (zlib for compressed pieces, as VTU) |

## Reading & writing

```python
import meshioplusplus

mesh = meshioplusplus.read("assembly.med")  # any mesh with 2+ cell blocks
meshioplusplus.vtm.write("assembly.vtm", mesh,
    binary=True,            # base64-encode each piece's arrays instead of writing text
    compression="zlib",     # "zlib", "lz4", "zstd", or None
    header_type=None,       # "UInt32" (default) or "UInt64"
)

back = meshioplusplus.vtm.read("assembly.vtm")
```

`meshioplusplus.read`/`write` dispatch on the extension, so `.vtm` needs no explicit format name.

## File structure

Writing `assembly.vtm` creates a *sibling directory* named after the index's own stem, holding one piece per cell block:

```
assembly.vtm
assembly/
  assembly_0.vtu
  assembly_1.vtu
```

```xml
<?xml version="1.0"?>
<VTKFile type="vtkMultiBlockDataSet" version="1.0" byte_order="LittleEndian">
<vtkMultiBlockDataSet>
<Block index="0">
<DataSet index="0" name="block_0" file="assembly/assembly_0.vtu"/>
<DataSet index="1" name="block_1" file="assembly/assembly_1.vtu"/>
</Block>
</vtkMultiBlockDataSet>
</VTKFile>
```

Each piece is a standalone, independently readable `.vtu` file — open `assembly/assembly_0.vtu` on its own in ParaView, or with `meshioplusplus.read()`, and it works exactly as any other `.vtu` would.

## The mesh side of the deal

- **`write`** carves the mesh into one piece per `CellBlock`: the piece carries that block's cells and cell_data, and the mesh's full point_data, then is pruned with `clean(remove_orphans=true)` so it only keeps the points that block actually references. A point shared by two blocks is therefore **duplicated** across their two pieces — each piece is a standalone `.vtu` by design, not a fragment that depends on its siblings. `field_data` is not carried onto pieces: it describes the whole mesh, not one block, and `.vtm` has no per-block place to put it.
- **`read`** parses the index, reads every piece with the best available reader for its own extension (`.vtu` or `.vtp`), and combines them with [`merge()`](/merge) using `weld=false` — multiblock pieces are pre-separated by construction, not coincident-point fragments to weld back together. Cell blocks of the same meshio++ type from *different* pieces consolidate into one output block, exactly as `merge()` always does. Each piece becomes one `RegionKind::Cell` region in the merged mesh, named from the index's `name=` attribute — correct even when same-typed pieces consolidate, since it comes from `merge()`'s own per-input cell index map rather than a naive block-count assumption.

```python
back = meshioplusplus.vtm.read("assembly.vtm")
for region in back.regions:
    print(region.name, region.kind, len(region.entries))  # one "cell" region per piece
```

Nested `<Block>` elements are read structurally (every `<DataSet>` anywhere under `<vtkMultiBlockDataSet>` is collected, in document order) but the nesting itself is not reproduced — there is nothing in the uniform mesh API to hold a block hierarchy, only the flat list of named regions above. A hand-written `.vtm` mixing `.vtu` and `.vtp` pieces reads correctly even though this writer always emits `.vtu`.

## Data mapping

| meshio++ | VTM |
|---|---|
| `points` | duplicated across pieces, once per block that references them |
| `point_data` | copied onto every piece (the full array, then pruned to that piece's own points) |
| `cell_data` | each piece gets its own block's slice |
| `field_data` | **not written** (no per-block place to put it) |
| named cell regions | recovered from `name=`, one per piece, not round-tripped on write |

## Quirks & limitations

Only `.vtu`/`.vtp` pieces are read; any other `file=` extension raises `ReadError` by name. A piece containing polyhedron cells is declined by the underlying `.vtu` writer, the same restriction [`.vtu`'s own page](./vtu.md#quirks-limitations) documents — a mesh with a polyhedron block cannot be written as `.vtm` either.

Regions are not round-tripped symmetrically: `write` does not look for existing regions to name pieces after (every piece is written as `"block_<i>"`), while `read` always attaches one region per piece. Writing a mesh that already carries `.vtm`-shaped regions and reading the result back therefore does not reproduce the original region names.

## Notes

Pieces reuse the exact same `.vtu` writer/reader (`write_vtu_codec`/`read_vtu`) every other VTK XML piece file in this library does — `--codec zlib|lz4|zstd` works here exactly as it does for `.vtu`. The pure-Python fallback (`meshioplusplus.vtm._vtm`) delegates each piece to the public `meshioplusplus.vtu`/`meshioplusplus.vtp` readers/writers (so pieces still get C++ acceleration when available) and reuses the internal numpy `merge()` implementation for combining them.

Testing exercises the roadmap's own stated probe — write two (in practice three, to also cover same-type block consolidation) blocks to `.vtm`, read them back as two meshes — in `tests/cpp/test_vtm.cpp`/`tests/python/test_vtm.py`, plus a check that each piece is independently readable on its own.

## See also

- [VTU](./vtu.md) — the format every `.vtm` piece actually is.
- [VTS](./vts.md) / [VTR](./vtr.md) — the lattice-restricted siblings, for a single structured block.
- [merge](/merge) — the operation `.vtm`'s reader is built on.

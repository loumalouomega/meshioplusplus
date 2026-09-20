# PVTU — VTK XML parallel unstructured grid (`.pvtu`)

The index of a partitioned dataset: it *declares* the arrays every piece holds and names one `.vtu` piece per part, so a mesh split across ranks opens in ParaView as one dataset (v14.1.0, roadmap §1.1). This is the on-disk face of [`partition`](/partition): `partition` → `.pvtu` → read-merge returns the original mesh, up to point ordering. See the [VTK XML file formats documentation](https://docs.vtk.org/en/latest/vtk_file_formats/vtkxml_file_format.html) and [`.pvtp`](./pvtp.md) for the `.vtp` twin.

| | |
|---|---|
| **Format name** | `pvtu` |
| **Extensions** | `.pvtu` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — (zlib for compressed pieces, as VTU) |

## Reading & writing

```python
import meshioplusplus

mesh = meshioplusplus.read("case.vtu")
labels = meshioplusplus.partition_labels(mesh, 4)      # one Int64 array per cell block
mesh.cell_data["partition:part"] = labels
meshioplusplus.pvtu.write("case.pvtu", mesh,           # carved by partition:part
    binary=True,            # base64-encode each piece's arrays instead of writing text
    compression="zlib",     # "zlib", "lz4", "zstd", or None
    header_type=None,       # "UInt32" (default) or "UInt64"
    part_key="partition:part",  # the integer cell_data array to carve by
)

merged = meshioplusplus.pvtu.read("case.pvtu")         # every piece merged, one region each
second = meshioplusplus.pvtu.read("case.pvtu", piece=1)  # or one piece alone
```

`meshioplusplus.read`/`write` dispatch on the extension, so `.pvtu` needs no explicit format name, and `read(..., piece=k)` works through the generic call (see [selective reads](/selective_read#picking-a-piece)).

To write the pieces `partition` returns — with their halo layers — hand them over as a list:

```python
pieces = meshioplusplus.partition(mesh, 4, ghost_layers=1)
meshioplusplus.pvtu.write_pieces("case.pvtu", pieces)
```

From the command line, a `partition` output path ending in `.pvtu` (or `-o pvtu`) with no `{part}` token writes one index over every part rather than one file per part:

```bash
meshioplusplus partition -n 4 --ghost-layers 1 case.vtu case.pvtu
```

## File structure

Writing `case.pvtu` creates a *sibling directory* named after the index's own stem, holding one piece per part, numbered with zero padding (the width `{step}` patterns use: at least four digits):

```
case.pvtu
case/
  case_0000.vtu
  case_0001.vtu
  case_0002.vtu
  case_0003.vtu
```

```xml
<?xml version="1.0"?>
<VTKFile type="PUnstructuredGrid" version="1.0" byte_order="LittleEndian">
<PUnstructuredGrid GhostLevel="1">
<PPointData>
<PDataArray type="Float64" Name="u"/>
<PDataArray type="UInt8" Name="vtkGhostType"/>
</PPointData>
<PCellData>
<PDataArray type="Float64" Name="c"/>
<PDataArray type="Int64" Name="partition:ghost"/>
<PDataArray type="UInt8" Name="vtkGhostType"/>
</PCellData>
<PPoints>
<PDataArray type="Float64" Name="Points" NumberOfComponents="3"/>
</PPoints>
<Piece Source="case/case_0000.vtu"/>
<Piece Source="case/case_0001.vtu"/>
</PUnstructuredGrid>
</VTKFile>
```

The index carries no geometry and no data, only the *declarations* of the arrays (`PPoints`, `PPointData`, `PCellData`: name, type, component count) and one `<Piece Source=>` per part. Each piece is a standalone, independently readable `.vtu` file. `Source=` is always relative, always written with forward slashes, and XML-escaped.

## The mesh side of the deal

- **`write`** carves the mesh by the integer `cell_data` array named by `part_key` (default `partition:part`, what `partition_labels` and `partition --labels-only` produce): one piece per part id `0..max`, each pruned to the points its own cells reference, so a point on a part boundary is **duplicated** into every piece that touches it. A part that owns no cells is still written as an (empty) piece, because skipping it would renumber the parts after it. Without the array the whole mesh is one piece. A non-integer or negative part array is refused.
- **`write_pieces`** is the primitive under it: already-carved pieces go in as they are. This is how `partition`'s output is written with its halo layers, since a single mesh cannot express a cell that is a ghost of several parts.
- **Every piece must declare identical arrays** — name, type and component count, for `Points`, point data and cell data. That is checked **before anything is created**, so a refusal leaves no directory and no half-written index; the error names the array, the piece and both declarations. A piece with no cells at all may omit its cell arrays (an idle rank allocated none).
- **`read`** parses the index, reads every piece with the best available `.vtu`/`.vtp` reader, and combines them with [`merge()`](/merge) using `weld=false`. Each piece becomes one `RegionKind::Cell` region named `piece_0`, `piece_1`, …. A point on a part boundary therefore appears once per piece; `clean` with `weld=true` fuses them:

```python
merged = meshioplusplus.read("case.pvtu")
whole = meshioplusplus.clean(merged, weld=True)   # the original mesh, up to point ordering
```

- **`piece=k`** keeps only piece `k` (negative counts from the end; out of range names the piece count) and attaches no region. A file with a single piece reads as that piece, with no region.
- Piece paths are resolved against the **index's own directory**, so a partitioned case can be moved as a folder. A path with spaces, `&` or other XML-escaped characters is read as written. An *absolute* path written on another machine that does not exist here is an error naming the attribute and the path — it is deliberately not searched for elsewhere, because a fallback could read the wrong file.
- A `.pvtu` names `.vtu` pieces (and, leniently, `.vtp`). It does not nest: another `.pvtu` as a piece is refused by name, so no cycle is expressible. The layer above is [`.pvd`](./pvd.md), which may name `.pvtu` files.

## Ghost cells

`partition(..., ghost_layers=N)` grows each piece by `N` layers of neighbouring parts' cells and tags them `partition:ghost` (0 = owned, `L` = reached at layer `L`). Writing translates that into VTK's own vocabulary:

- `vtkGhostType` on **cells** (`UInt8`): `0` for an owned cell, `DUPLICATECELL` (1) for a cell in any halo layer.
- `vtkGhostType` on **points**: `DUPLICATEPOINT` (1) when no owned cell of the piece uses the point, else `0`.
- `GhostLevel="N"` on the index, the deepest layer present.
- `partition:ghost` itself is written too, as an ordinary `Int64` cell array: `vtkGhostType` collapses layer 2 onto layer 1, so keeping both makes meshio++ → meshio++ exact, and ParaView simply ignores the extra array.

An array you already named `vtkGhostType` is passed through unchanged (all of VTK's bits survive, `REFINEDCELL` = 8 and `HIDDENCELL` = 32 included), and nothing is fabricated when there is neither. The caller's meshes are never modified.

On read the ghost cells are **kept by default** — a reader must not silently discard data, and keeping is what makes read → write round-trip. `ghosts="drop"` removes every cell with any `vtkGhostType` bit set, and the points only those cells used, from each piece *before* merging, and then removes the ghost arrays:

```python
kept = meshioplusplus.pvtu.read("case.pvtu")                    # more cells than the original: the halo is real data
whole = meshioplusplus.pvtu.read("case.pvtu", ghosts="drop")    # the partition of unity again
```

Merging does not deduplicate ghost cells by itself; `ghosts="drop"` is the precise tool for that. From C++ the same choice is `PvtuReadOptions::mGhosts` (`GhostPolicy::Keep` / `Drop`), a defaulted trailing parameter of `read_pvtu` rather than a `ReadOptions` member, so `ReadOptions` (and the ABI) does not change.

## Data mapping

| meshio++ | PVTU |
|---|---|
| `points` | duplicated into every piece that references them |
| `point_data` | copied onto every piece and declared under `PPointData` |
| `cell_data` | each piece gets its own cells' slice, declared under `PCellData` |
| `partition:part` | the carving key; kept as an ordinary array in each piece |
| `partition:ghost` | kept, and translated to `vtkGhostType` + `GhostLevel` |
| `field_data` | **not written** (the `.vtu` writer does not emit it) |
| named cell regions | recovered as `piece_<i>` on read, not round-tripped on write |

## Quirks & limitations

- **`field_data` is not carried.** The `.vtu` writer emits none, so it cannot travel through a piece, and the C++ `.vtu` reader does not read a piece's `<FieldData>` (the Python reference reader does). Both are properties of `.vtu`, not of the index.
- **The Python reference handles rectangular cell blocks only** when carving, deriving ghost arrays or dropping ghosts; a polygon or polyhedron block raises `NotImplementedError` there, and the C++ core (the normal path) handles them.
- **Regions are not round-tripped symmetrically**: `write` does not look at existing regions, while `read` always attaches one per piece.
- The index's own `byte_order` and `header_type` attributes are not used for reading: each piece's own header is authoritative, so an index and its pieces may disagree.
- Opening the index in ParaView is not covered by the automated tests; what is covered is VTK's own `vtkXMLPUnstructuredGridReader` / `vtkXMLPPolyDataReader` reading every file either engine writes, and both engines reading what `vtkXMLPUnstructuredGridWriter` writes.

## Notes

The pieces reuse the exact `.vtu` writer/reader every other VTK XML piece file uses — `compression="zlib|lz4|zstd"` works here exactly as it does for `.vtu`. The pure-Python reference (`meshioplusplus.pvtu._pvtu`, over `meshioplusplus._pvtk_index`) delegates each piece to the public `meshioplusplus.vtu`/`vtp` readers/writers, so pieces still get C++ acceleration when the core is present; both engines are tested against each other and against VTK.

## See also

- [`.pvtp`](./pvtp.md) — the same index over `.vtp` pieces.
- [`.pvd`](./pvd.md) — a time-indexed collection, whose entries may be `.pvtu` files.
- [VTU](./vtu.md) — what every piece actually is; [VTM](./vtm.md) — the index that carves by cell *block* instead of by part.
- [partition](/partition) — the operation whose output this is, halo layers included.
- [merge](/merge) — the operation the reader is built on; [selective reads](/selective_read#picking-a-piece) — `piece=`.

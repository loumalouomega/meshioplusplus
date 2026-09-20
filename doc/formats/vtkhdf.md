# VTKHDF — Kitware's HDF5-based VTK format (`.vtkhdf`)

The one ParaView-native format that puts time, partitions and fields in a single file: geometry and connectivity once, every field as a flat array, and small offset tables saying which slice of each array belongs to which step and which piece. That is the shape of a physics-ML dataset — one trajectory is one file, with random access per step. See the [VTKHDF documentation](https://docs.vtk.org/en/latest/vtk_file_formats/vtkhdf_file_format/index.html).

| | |
|---|---|
| **Format name** | `vtkhdf` |
| **Extensions** | `.vtkhdf`; `.hdf` is also registered (see [the `.hdf` extension](#the-hdf-extension)) |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | HDF5 (C++ core built with `MESHIOPLUSPLUS_WITH_HDF5`, otherwise `h5py`) |

## Reading & writing

```python
import meshioplusplus

meshioplusplus.vtkhdf.write("case.vtkhdf", mesh,
    compression="gzip",        # "gzip" (default) or None; datasets under 4 KiB are never compressed
    compression_opts=4,        # gzip level 0-9
    dataset_type="UnstructuredGrid",  # or "PolyData", "PartitionedDataSetCollection", "MultiBlockDataSet"
    version=None,              # None = the oldest version that covers the content; or (major, minor)
)

mesh = meshioplusplus.read("case.vtkhdf")                  # every piece merged, step 0
last = meshioplusplus.read("run.vtkhdf", time_step=-1)     # the last step of a transient file
part = meshioplusplus.read("case.vtkhdf", piece=1)         # one partition, not the merged mesh
```

`meshioplusplus.read`/`write` dispatch on the extension, so `.vtkhdf` needs no explicit format name. The reader also takes `points_only` and `arrays` (narrowing what is materialized) and, on `meshioplusplus.vtkhdf.read` itself, `lenient` (see [PolyData](#polydata)).

A transient file is best written through the [time-series writers](../vtkhdf_time_series.md), or through the [sequence engine](../sequences.md) (`convert 'out_*.vtu' series.vtkhdf`), which fans a set of per-step files into one.

## Supported types

| `Type` | Read | Write | Notes |
|---|---|---|---|
| `UnstructuredGrid` | ✓ | ✓ | polyhedra, partitions, time |
| `PolyData` | ✓ | ✓ | `Vertices`/`Lines`/`Polygons`/`Strips`; never chosen automatically on write |
| `PartitionedDataSetCollection` | ✓ | ✓ | blocks are `UnstructuredGrid`/`PolyData`; no time |
| `MultiBlockDataSet` | ✓ | ✓ | nesting is flattened; no time |
| `ImageData`, `OverlappingAMR`, `HyperTreeGrid`, `Table`, `RectilinearGrid`, `StructuredGrid` | — | — | refused by name |

## The mesh side of the deal

| meshio++ | VTKHDF |
|---|---|
| `points` | `Points`, `(n, 3)`; a 2-D mesh is zero-padded. The stored dtype is kept on read (VTK writes `float32`) and the mesh's own dtype on write |
| cell blocks | `Connectivity`/`Offsets`/`Types`, with VTK's ids and node order. The linear wedge is the one type whose node order differs from meshio's and is permuted both ways |
| `point_data` | `PointData/<name>`, `(n,)` or `(n, k)`; higher rank is folded to `(n, prod)`; a `bool` array is stored as `uint8` |
| `cell_data` | `CellData/<name>`, one flat array over all cells (concatenated in block order); split back per block on read |
| `field_data` | `FieldData/<name>`; a non-numeric entry is skipped with a warning on write |
| named cell regions | one per partition or composite block on read (`piece_<i>`, or the block's name); the block layout on write of a composite type |
| provenance | the `VTKHDF` group's `meshioplusplus:provenance` attribute (plain UTF-8 text) |

`Type` is written as a **fixed-length ASCII** string attribute and `Version` as a two-element `int64` array. `Type` as HDF5's default variable-length UTF-8 has tripped readers, so it is never written that way; both forms are accepted on read.

## Versions

The version written is the **oldest that covers the content**, never a fixed maximum: `2.0` for a plain mesh (the version that added time), `2.1` for a composite, `2.5` once a polyhedron is present. Pin one with `version=(major, minor)`; a version too old for the content is a `WriteError` naming the feature (`polyhedral cells needs VTKHDF 2.5 or newer`). VTK master documents `2.8`, and this build knows minors up to 8.

The reader accepts major versions 1 and 2 and **dispatches on features, not on the declared minor**: if `FaceConnectivity` exists it decodes polyhedra whatever `Version` claims, and if `Steps` exists it honours it. Files in the wild declare sloppy versions, so the minor only produces a good message — a newer minor than 8 is read with a warning, a different major is refused.

## Partitions

A file with several partitions (VTK's `vtkPartitionedDataSet`, or `partition` output written piece by piece) has one entry per piece in every `NumberOf*` dataset, with the big arrays concatenated. `Offsets` holds `NumberOfCells + 1` entries **per piece**, and connectivity holds **piece-local** point ids.

By default the pieces are **merged** into one mesh, in order, and one cell region per piece is attached, named `piece_0`, `piece_1`, …; the region's `tag` is the piece's position. Pass `piece=k` (negative counts from the end) to read **one piece alone**, in which case no regions are attached. An out-of-range piece is an error naming the piece count. A format whose reader has no pieces refuses `piece` rather than returning the merged mesh.

Polyhedra are bucketed into `polyhedron<N>` blocks (N = the cell's unique node count) **per piece**, so a piece's cells always occupy one contiguous range of the merged mesh, and reading `piece=k` alone reproduces region `k` of the merged read exactly. Runs of any other cell type join across piece boundaries.

## Composites

A `PartitionedDataSetCollection` or `MultiBlockDataSet` keeps each block as a flat group directly under `VTKHDF`, with a hierarchy of soft links in an `Assembly` group. On read the blocks are merged in Assembly order — never re-welded, never regrouped by type — and each becomes one cell region named after its Assembly link (a nested block keeps its own name; nesting is flattened). A block that is itself partitioned yields one region per partition, `<block>/piece_<j>`. Pass `piece=k` to select across the flattened `(block, partition)` list.

On write the mesh is carved into blocks along its **cell regions** when they partition the cells (disjoint, covering every cell), and along its cell blocks otherwise (named `block_<i>`, with a warning if regions were present but unusable). Each block is pruned to the points it references. Blocks are written in `(tag, name)` order, which is what a read composite's regions carry, so a read→write round trip restores the order; regions without tags come out in name order. Both engines use the same rule, because the C++ mesh stores regions sorted by name and cannot remember an insertion order.

Three measured constraints of `vtkHDFReader` shape the writer, all confirmed against VTK 9.7 rather than taken from the prose specification:

- the `VTKHDF` group, the `Assembly` group and its children must **track link creation order** — a reader given an untracked one aborts the whole process, so the writer creates them with `H5P_CRT_ORDER_TRACKED | H5P_CRT_ORDER_INDEXED`;
- a `PartitionedDataSetCollection` block needs its **`Index`** attribute — without it the reader returns empty blocks;
- a composite's **root may hold only blocks and the `Assembly`**: a root `FieldData` group is read as a phantom block and the read fails. `field_data` is therefore replicated onto every block's own `FieldData` group and read back from the first block.

## Time

The `Steps` group (`NSteps` attribute, `Values`, and per-step offset tables) turns the flat arrays into a time series. `time_step=k` selects a step (negative counts from the end; out of range is an error naming the step count) and attaches that step's time as `field_data["meshio:time"]`, so a fanned-out step carries its time like every other format. A static file has exactly one step. `read_metadata(path)["time_values"]` lists every step's time without touching the geometry.

A static mesh is written **once**: every step's `PartOffsets`, `PointOffsets`, `CellOffsets` and `ConnectivityIdOffsets` are zero, while `PointData/<name>` gains `NumberOfPoints` rows per step and `PointDataOffsets/<name>` counts `k * NumberOfPoints`. `vtkHDFReader` plays this layout: `TIME_STEPS` is right, the geometry is identical at every step and each array is exact. (VTK's own writer repeats the geometry for every step instead.)

Offsets are read tolerantly: the data offset tables may hold `NSteps` or `NSteps + 1` entries (VTK writes the latter), `CellOffsets` may be `(N,)` or `(N, 1)` for an unstructured grid and is `(N, 4)` for `PolyData`, and `Values` may be `float32` or `float64`. For polyhedra the face-level starts are derived from the `NumberOf*` counters rather than from `FaceConnectivityOffsets` and its siblings, because VTK's writer leaves the last piece out of those tables for partitioned steps.

**A transient composite is not supported.** A `PartitionedDataSetCollection`/`MultiBlockDataSet` (or a block) carrying a `Steps` group with more than one step is refused by name, because its layout could not be measured: `vtkHDFWriter` in VTK 9.7 crashes on it. Only `UnstructuredGrid` and `PolyData` files are transient.

## Polyhedra

A polyhedral cell is stored with four arrays: for cell `i`, `PolyhedronToFaces[PolyhedronOffsets[i] : PolyhedronOffsets[i+1]]` lists its face ids, and face `f` is `FaceConnectivity[FaceOffsets[f] : FaceOffsets[f+1]]`. `PolyhedronOffsets` has `NumberOfCells + 1` entries and covers **every** cell — a non-polyhedral cell is a zero-length span — which is how a mixed mesh (an OpenFOAM mesh always mixes hexahedra, polyhedra and boundary faces) works. The cell's own `Connectivity` row is its **sorted unique node set**, exactly what the [`.vtu`](./vtu.md) writer emits.

VTKHDF's face table may be shared between neighbouring polyhedra; meshio++ owns faces per cell, so on read the sharing is expanded and not recorded, and on write faces are appended without de-duplication (correct, at some size cost). See [polyhedra](../polyhedra.md).

## PolyData

`Vertices`, `Lines`, `Polygons` and `Strips` each hold their own connectivity set. On read, `Polygons` rows of 3, 4 or more nodes become `triangle`, `quad` and `polygon`; a `Vertices` row of one node is a `vertex` and a `Lines` row of two nodes a `line`. VTK's cell order — and the order of `CellData` — is **piece-major** (each partition's cells in category order), while the topology arrays are category-major; the reader interleaves accordingly, which is what keeps `cell_data` aligned.

Constructs with no meshio++ type — **poly-vertex, poly-line and triangle strips** — are a `ReadError` naming the construct. `lenient=True` (`ReadOptions::mLenient`) skips them with a warning instead, and drops their `cell_data` rows with them. Triangulating a strip is a deliberate non-goal: it changes the cell count, which is [`convert_cells`](../convert_cells.md)' job.

On write, `dataset_type="PolyData"` regroups the blocks into VTK's canonical order (`Vertices`, `Lines`, `Polygons`, `Strips`) and reorders `cell_data` to match, so a mesh whose blocks are in another order does not round-trip block-for-block. Any block that is not a vertex, line, triangle, quad or polygon is a `WriteError` naming it and pointing at `UnstructuredGrid` — deliberately stricter than the [`.vtp`](./vtp.md) writer.

## The `.hdf` extension

`.hdf` is a generic HDF5 extension, so it is registered but not sniffed: a file with no `/VTKHDF` group is refused with a message telling the caller to name the format explicitly (`meshioplusplus.read(path, file_format="vtkhdf")`, or the right format for the file). The content sniffer deliberately does not claim generic HDF5 magic.

## Quirks & limitations

- `PolyData`, composites and polyhedral files fall back to a full read for `read_metadata` (their block structure needs the arrays); a polyhedron-free `UnstructuredGrid` is summarized from `Steps/Values`, the counters and `Types` alone.
- `compress`/`decompress` do **not** support `vtkhdf`: they rewrite a file in place from a single merged read, which would silently drop every step but the first from a transient file. `convert` to a new file with `--compression gzip|none` (MCP: `compression`) is supported.
- A partitioned write from `partition` is not produced directly: write each piece as a composite block (`dataset_type="PartitionedDataSetCollection"` with one cell region per piece) or use the [time-series writers](../vtkhdf_time_series.md) for a transient case.
- `NumberOfPoints`, `Points` and the other partition datasets are written for a single piece; a plain mesh has one partition.

## C++ core and the Python reference

Both engines read and write everything on this page and are cross-checked: a sweep of 28 files (VTK-written and meshio++-written) across 14 read-option combinations gives identical results from the two readers, and every file the C++ writer produces reads back through the Python reference. The C++ path is used whenever the core is built with HDF5 and the file is a path; the `h5py` reference is used otherwise, or when the core declines a file (`MESHIOPLUSPLUS_STRICT_CORE=1` turns a decline into an error).

The C++ API is `write_vtkhdf(path, mesh, gzip_level, type, version)`, `read_vtkhdf(path[, ReadOptions])` and `read_vtkhdf_metadata` in `meshioplusplus/formats/vtkhdf.hpp`, plus `VtkhdfTimeSeriesWriter` in `vtkhdf_time_series.hpp`. `ReadOptions` gained `mPiece`/`mPieceSet` for this format (ABI 14); see [selective reads](../selective_read.md).

## See also

- [VTU](./vtu.md) — the XML sibling whose cell-type tables, node order and polyhedron convention this format reuses.
- [VTM](./vtm.md) — the multi-file equivalent of a composite.
- [VTKHDF time series](../vtkhdf_time_series.md) — the stateful writers.
- [Sequences](../sequences.md) — how a set of per-step files becomes one transient `.vtkhdf`.
- [Regions](../regions.md), [polyhedra](../polyhedra.md), [selective reads](../selective_read.md).

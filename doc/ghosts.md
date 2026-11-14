# Ghosts, owners and global ids

A decision record (roadmap 10.4): the names meshio++ gives to the arrays a partitioned or numbered mesh carries, which formats read and write them, and what an operation does with them. The names live once in `src/cpp/src/detail/ghost_names.hpp` and its Python twin `meshioplusplus/_ghosts.py`; `tests/python/test_ghost_names.py` keeps the two in step. meshio++ stays serial (see [partition](./partition.md#meshio-is-serial)): these arrays describe the pieces a rank holds, they never imply a communicator.

## The names

| Array | Where | Dtype | Meaning |
|---|---|---|---|
| `partition:ghost` | `cell_data` | Int64 | `0` for a cell its piece owns, `L` for a halo cell reached at layer `L`. The cross-format name: an operation that cares about ghosts looks for this one. |
| `partition:part` | `cell_data` | Int64 | The 0-based piece a cell belongs to. Written by `partition_labels`, the `Partition` pipeline step and the partitioned readers (Elmer, MFEM); writers that split by piece (pvtu, pvtp, Elmer, OpenFOAM, Tecplot, VTX) carve by it. |
| `partition:owner` | `cell_data`, `point_data` | Int64 | Reserved: the owning piece of a halo cell or a shared point. No reader or operation writes it yet; it is named here so the Nemesis reader (roadmap 2.3.8.2) and the merge helper (2.3.8.1) agree on it. |
| `partition:original_point_id`, `partition:original_cell_id` | `point_data`, `cell_data` | Int64 | 0-based indices into the mesh a piece was cut from (`partition(..., record_ids=True)`). They survive a `.pvtu` round trip, so pieces can be put back together exactly (below). |
| `global:node_id`, `global:cell_id` | `point_data`, `cell_data` | Int64 | A file's own global numbering, 1-based or not, as the file states it. Exodus reads `node_num_map` and `elem_num_map` into them and writes them back. |
| `vtkGhostType` | `cell_data`, `point_data` | UInt8 | VTK's bitfield (`DUPLICATECELL = 1`, `DUPLICATEPOINT = 1`), kept verbatim by the VTK readers and always written as UInt8. |
| `vtkOriginalPointIds` | `point_data` | Int64 | VTK's global point ids; the VTX reader welds ghost points through them. |
| `gmsh:ghost` | — | — | Reserved for Gmsh's own flag; today the Gmsh reader turns a kept `$GhostElements` section into `partition:ghost` (`1` for a ghost) directly. |

A format-native array (`vtkGhostType`, a vendor's own id array such as `med:num`, `mdpa:id` or `dpf:node_id`) is never renamed: the same-format round trip stays byte-identical, and the `global:` names are used where a format's ids are global ids in this sense.

## Which formats carry them

| Format | Reads | Writes |
|---|---|---|
| `pvtu`, `pvtp`, `pvd` | `vtkGhostType` of each piece; `ghosts="drop"` removes the halo (and `partition:ghost` when it is then all zero) | `vtkGhostType` derived from `partition:ghost`, with `GhostLevel`; `partition:ghost` itself too, so the layer depth survives |
| `vtu`, `vtp` | `vtkGhostType` and every other array as stored | as stored; `vtkGhostType` forced to UInt8 |
| `vtx` | `vtkGhostType`, `vtkOriginalPointIds`; `ghosts="drop"` welds ghost points onto their owners | both, from `partition:part` |
| `gmsh` | `$GhostElements` → `partition:ghost` | the section, kept verbatim |
| `exodus` | `node_num_map` → `global:node_id`, `elem_num_map` → `global:cell_id` (an identity map is not kept) | the maps, from those arrays; they are not also written as variables |
| `elmer`, `mfem` | `partition:part` for a partitioned directory or set of files | Elmer: a `partitioning.N` directory from `partition:part` |
| `cgns` | none: SIDS has no ghost flag for unstructured zones, and `Rind` (structured only) is stripped | — |

## What operations do with them

- `partition` writes `partition:ghost` (with `ghost_layers > 0`) and the `partition:original_*` ids (with `record_ids=True`). Pieces overlap when ghosted, so `partition_labels` refuses `ghost_layers`.
- Data-only operations and `transform` pass them through unchanged; subset operations (`crop`, `split`, `clean`, `extract_*`) carry them as they carry any `cell_data`/`point_data` row, so an id keeps meaning the cell or point it named.
- `merge` concatenates them; pieces of one partition merged this way still hold their ids, which is what reassembly needs.
- Reading a `.pvtu` with `ghosts="drop"` and ordering by the `partition:original_*` ids gives back the mesh `partition` cut, exactly (`tests/python/test_pvtu.py::test_ghosted_pieces_rebuild_the_serial_mesh_from_their_ids`). The general helper that merges pieces by global id, whatever produced them, is roadmap 2.3.8.1, and takes these names as its inputs.

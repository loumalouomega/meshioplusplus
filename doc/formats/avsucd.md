# AVS-UCD (`.avs`)

The [AVS Unstructured Cell Data](https://lanl.github.io/LaGriT/pages/docs/read_avs.html) format: an ASCII header of counts, node coordinates, cells (with a per-cell material id), then optional node- and cell-data sections.

| | |
|---|---|
| **Format name** | `avsucd` |
| **Extensions** | `.avs` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshioplusplus

mesh = meshioplusplus.read("mesh.avs")
meshioplusplus.avsucd.write("out.avs", mesh)
```

`write` takes no keyword arguments.

## File structure

`#`-comment lines are skipped. First non-comment line: 5 integers `num_nodes num_cells num_node_data num_cell_data 0` (the trailing field is unused).

- **Nodes**: `num_nodes` rows of `id x y z`. `id` is an **arbitrary integer** (not necessarily sequential or 1-based); an id→index map is built while reading.
- **Cells**: `num_cells` rows of `id material_id avsucd_type_name node_id0 node_id1 ...`. `material_id` (also an arbitrary integer) becomes `cell_data["avsucd:material"]`; node ids are resolved through the same id→index map, so cells may reference nodes by their original file ids. Node count per row is simply "everything after the first 3 fields" — no fixed per-type count table is consulted on read.
- **Node-data / cell-data** (present only if the corresponding header count is nonzero): a header line `n_arrays size_0 size_1 ... size_{n-1}` (component count per array), then `n_arrays` label lines of the form `"<name>, real"` (the `", real"` suffix is discarded — data is always treated as float), then `num_entities` rows of `id v0 v1 ... `, again resolved through the id map.

Write re-numbers everything sequentially from 1, in the same section order, with `%.14e` (Python) / `%.14e` and `%.17g` (C++: data vs. points respectively) float precision.

## Cell types

| AVS-UCD | meshio++ | AVS-UCD | meshio++ |
|---|---|---|---|
| `pt` | `vertex` | `tet` | `tetra` |
| `line` | `line` | `pyr` | `pyramid` |
| `tri` | `triangle` | `prism` | `wedge` |
| `quad` | `quad` | `hex` | `hexahedron` |

Node-order permutation, meshio++ → AVS-UCD:

| type | permutation |
|---|---|
| `tetra` | `[0, 1, 3, 2]` |
| `pyramid` | `[4, 0, 1, 2, 3]` |
| `wedge` | `[3, 4, 5, 0, 1, 2]` |
| `hexahedron` | `[4, 5, 6, 7, 0, 1, 2, 3]` |

AVS-UCD → meshio++ is the same table for `tetra`/`wedge`/`hexahedron` (all involutions), but **not** for `pyramid`, which uses the (functionally correct but textually different) inverse `[1, 2, 3, 4, 0]`.

## Data mapping

- `cell_data["avsucd:material"]` — the cell record's material id, one array per block.
- Any other point_data/cell_data name comes directly from the data section's label lines.

## Quirks & limitations

- Node and cell ids in the file are **arbitrary integers**; both read and write maintain explicit id↔index maps so files with sparse, reordered, or non-contiguous numbering are handled correctly (write always renumbers sequentially from 1).
- On write, `avsucd:material` is chosen as the **first** integer-typed cell_data array found; if others exist they're dropped, with a warning in the Python writer (`"AVS-UCD can only write one cell data array... Skipping ..."`) but silently in the C++ writer.
- Cell-data arrays spanning multiple cell blocks are read as one flat array then re-split via cumulative block-length offsets — this assumes the blocks are contiguous in the order they were originally read.
- 2D points are promoted to 3D on write with a warning, and the Python writer does this **in place** on the caller's `Mesh.points` — a user-visible side effect worth being aware of.
- Data-array label cleanup (`strip()` + replace spaces with `_`) is not reversible — a name with meaningful internal spaces is altered irrecoverably on read.

## Notes

- Fully handled by the C++ core.
- No reference fixture exists under `tests/meshes/avsucd/`; tests round-trip `empty_mesh`, `tri_mesh`, `quad_mesh`, `tri_quad_mesh`, `tet_mesh`, `hex_mesh`, plus a data variant naming `"avsucd:material"` explicitly alongside a scalar and a 3-vector array.

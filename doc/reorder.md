# Mesh reordering / renumbering

`meshioplusplus.reorder(mesh, method="rcm")` renumbers a mesh's nodes and elements to **reduce sparse-matrix bandwidth** and **improve cache locality**. It is a *pure permutation*: no nodes or cells are added or removed, the geometry and all data are preserved, and the applied permutations are returned so you can remap external per-node / per-cell arrays. It is a mesh **operation** (like [surface extraction](/extract_surface) and [quality metrics](/mesh_quality)), not a file format, and uses only standard C++/numpy, so it runs under every mesh backend.

![Node-pair sparsity before and after RCM reordering — the diagonal band is the point](/images/reorder_bandwidth.png)

```python
import meshioplusplus

mesh = meshioplusplus.read("part.vtu")
print(meshioplusplus.compute_bandwidth(mesh))          # before

out = meshioplusplus.reorder(mesh, method="rcm")
print(meshioplusplus.compute_bandwidth(out))           # after (smaller)

meshioplusplus.write("part_reordered.vtu", out)

# also return the permutations (old index -> new index)
out, node_perm, cell_perms = meshioplusplus.reorder(
    mesh, method="rcm", return_permutation=True
)
external_dof = external_dof[node_perm]                  # remap your own arrays
```

## Methods

| `method` | strategy | uses | good for |
|---|---|---|---|
| `"rcm"` (default) | Reverse Cuthill–McKee on the node-adjacency graph | connectivity | reducing matrix bandwidth / profile |
| `"morton"` | Morton (Z-order) space-filling curve | coordinates | cache locality, connectivity-free |
| `"hilbert"` | Hilbert space-filling curve | coordinates | cache locality (better than Morton) |

- **RCM** builds the node graph (nodes sharing a cell are neighbours) as a CSR structure, finds a pseudo-peripheral start node per connected component (via the rooted-level-structure heuristic), runs a Cuthill–McKee BFS visiting neighbours in increasing degree order, and reverses the result. Disconnected components and isolated/orphan nodes are handled.
- **Morton / Hilbert** quantize the node coordinates to 21 bits per axis over the mesh bounding box and sort the nodes along the curve. These optimise *spatial* locality (spatially close nodes get nearby indices); they do not target bandwidth the way RCM does.

## What gets permuted

- **Points** are reordered by the node permutation.
- **Connectivity** of every cell block is remapped to the new node indices, and each block's cells are reordered by their minimum new node index (deriving a per-block element permutation).
- **`point_data`** rows follow their nodes; **`cell_data`** rows follow their cells.
- **`field_data`** is preserved unchanged (mesh-global).
- **`point_sets` / `cell_sets`** indices are remapped to the new numbering (done in the Python layer).

`mesh.info` and `gmsh_periodic` are not carried through the operation.

## Return value

By default `reorder` returns the renumbered `Mesh`. With `return_permutation=True` it returns `(mesh, node_permutation, cell_permutations)`:

- `node_permutation` — an `int64` array, `old index → new index` (`new_points[node_permutation[i]] == old_points[i]`), zero-copy from the C++ core where possible.
- `cell_permutations` — one `int64` array per cell block, `old → new`.

## Bandwidth

`meshioplusplus.compute_bandwidth(mesh)` returns the connectivity bandwidth: the maximum over all cells of `(max node index − min node index)` within a cell. It is a convenient before/after measure of an RCM reordering.

## CLI

```bash
meshioplusplus reorder input.vtu output.vtu --method rcm       # default
meshioplusplus reorder input.vtu output.vtu --method morton
meshioplusplus reorder input.vtu output.vtu --method hilbert
meshioplusplus reorder input.vtu output.vtu --method rcm --report
```

`--report` prints the connectivity bandwidth before and after. See the [CLI reference](/cli).

## Other languages

The operation is exposed across every binding surface:

- **C API** — `mio_reorder(mesh, method)` returns a `mio_reorder_result` (borrow the mesh with `mio_reorder_result_mesh`, read the permutations with `mio_reorder_result_node_perm` / `..._cell_perm`, or take ownership of the mesh with `mio_reorder_result_take_mesh`); `mio_compute_bandwidth(mesh)`. See the [C API reference](/c_api).
- **Fortran** — `mesh%reorder(method, node_perm=...)` and `mesh%compute_bandwidth()`. See the [Fortran reference](/fortran).
- **WebAssembly / JavaScript** — `reorder(mesh, method)` returns `{mesh, nodePermutation, cellPermutations}`; `computeBandwidth(mesh)`. See the [WebAssembly reference](/wasm).

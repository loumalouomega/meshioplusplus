# Adjacency graphs and graph files

`adjacency(mesh, kind)` returns a mesh's node or cell graph in compressed-sparse-row form, `(indptr, indices)`, and `write_graph` writes it as a METIS or Scotch graph file. The graphs are the ones the core already builds for `reorder`, `smooth` and `partition`, exposed: the sparsity pattern and preallocation of an assembled finite-element matrix, the input of an external partitioner, and the CSR a graph network or a sparse solver takes directly.

```python
import meshioplusplus as mp
import scipy.sparse as sp

mesh = mp.read("bracket.vtu")

indptr, indices = mp.adjacency(mesh)              # nodes joined by a cell edge
indptr[-1] // 2                                    # the number of edges
indices[indptr[7]:indptr[8]]                       # the neighbours of node 7, ascending

indptr, indices = mp.adjacency(mesh, "clique")    # nodes sharing a cell: the FEM pattern
pattern = sp.csr_matrix((np.ones(len(indices)), indices, indptr))

indptr, indices = mp.adjacency(mesh, "cell")      # cells sharing a facet: the dual graph
mp.write_graph("bracket.graph", mesh, "cell")      # METIS, for gpmetis or KaHIP
mp.write_graph("bracket.grf", mesh, "cell")        # Scotch, for gpart
```

## The CSR contract

The neighbours of vertex `i` are `indices[indptr[i]:indptr[i + 1]]`. Both arrays are int64 and 0-based; every row is sorted ascending, with no duplicate and no self-loop, and every undirected edge is stored in both directions, so `indptr[-1]` is twice the edge count and the matrix the pair describes is symmetric. The graph is built through the uniform mesh API, so every mesh backend gives the same one, and its fill passes are serial with independent per-row sorts, so it is byte-identical on every parallel backend and thread count. The input mesh is never modified.

## Kinds

| kind | vertices | edge between two vertices when | consumer |
| --- | --- | --- | --- |
| `"node"` (default) | points | they are the ends of a cell edge | graph networks, smoothing, Laplacians |
| `"clique"` | points | they belong to one cell | the sparsity pattern (without the diagonal) and preallocation of an assembled FEM matrix; bandwidth reordering |
| `"cell"` | cells, numbered block-major | they share a facet: a face when the mesh holds any 3-D cell, else an edge | METIS, Scotch, KaHIP and ParMETIS; `partition`'s KaHIP method cuts this graph |

- **`"node"`** takes its edges from the linear types' edge tables (line, triangle, quad, tetra, pyramid, wedge, hexahedron), the rings of polygon rows (ragged rows and uniform `polygon*` blocks alike) and the face rings of polyhedra. The higher-order family, the VTK Lagrange types and `custom` have no tabulated edges and contribute none, so their nodes are isolated unless another block joins them; linearize with [`convert_cells`](/convert_cells) to keep their corners. It is [`edge_index`](/ml#graphs-for-gnns-edge-index)'s node graph in CSR form, except that `adjacency` also follows polyhedron faces, which `edge_index` skips with a warning.
- **`"clique"`** contributes every cell type, higher-order and polyhedral cells included: a `tetra10`'s ten nodes are all coupled, as they are in its element matrix.
- **`"cell"`** numbers cells in the global block-major order that regions and [`partition_labels`](/partition) use. Cells of lower dimension than the mesh (the line blocks beside a surface, the boundary triangles beside a volume) are isolated vertices. Polyhedra take part through their own faces and polygons through their ring edges, so a hexahedron and a polyhedron that meet on a face are neighbours. A facet shared by more than two cells joins every later cell to the first one, as the KaHIP graph always has.

Connectivity entries outside the points are ignored, as the shared adjacency kernels ignore them.

## Graph files

`write_graph(path, graph, kind="node", *, file_format=None)` writes a `Mesh`'s `kind` graph, or any `(indptr, indices)` pair, in one of two plain-text dialects; `file_format` defaults to the extension's.

| `file_format` | extension | layout |
| --- | --- | --- |
| `"metis"` | `.graph`, `.metis` | a header `n m` (vertices, undirected edges), then one line per vertex listing its **1-based** neighbours; an isolated vertex is an empty line. METIS, ParMETIS, KaHIP and mt-metis read it. |
| `"scotch"` | `.grf` | the version `0`, then `n 2m`, then `0 000` (base 0, no labels or weights), then one line per vertex: its degree and its **0-based** neighbours. |

Both describe an undirected simple graph without weights, so the writer refuses an asymmetric graph, a self-loop, a neighbour listed twice and a neighbour outside the vertices, naming the row. It returns `{"format", "num_vertices", "num_edges", "max_degree"}`. The files are graph exports, not mesh formats: they are not in the [format registry](/formats), `read` does not read them, and a partition computed from one is attached back through `partition_labels`' numbering.

## Engines

The C++ core builds the graph and writes the files; the Python twins reproduce both byte for byte (`tests/python/test_adjacency.py` holds them to it) for a pure-Python install, except the cell graph of a polyhedral mesh, which raises `NotImplementedError` naming the core.

## CLI

```sh
meshioplusplus graph IN OUT.graph [--kind node|clique|cell] [--format metis|scotch] [--json]
```

Both CLIs write the same bytes; `--json` reports `output`, `kind`, `format`, `num_vertices`, `num_edges` and `max_degree`. The MCP server's `graph` tool takes the same arguments.

## Other languages

The flat bindings return the same graph; Fortran, Julia and R shift both arrays to **1-based** (METIS's Fortran numbering and the layout of a Julia `SparseMatrixCSC`'s `colptr`/`rowval`), C and WebAssembly keep them 0-based.

```c
mio_adjacency* g = mio_adjacency_create(mesh, "cell");   /* NULL kind: "node" */
int64_t rows, nnz;
const int64_t* indptr = mio_adjacency_indptr(g, &rows);  /* rows = n + 1 */
const int64_t* indices = mio_adjacency_indices(g, &nnz);
mio_adjacency_free(g);
mio_write_graph(mesh, "cell", "part.graph", NULL, &nv, &ne, &max_degree);
```

```fortran
call m%adjacency(indptr, indices, kind='clique')   ! 1-based
call m%write_graph('part.graph', kind='cell')
```

```julia
g = adjacency(mesh; kind="cell")   # (; indptr, indices), 1-based
write_graph("part.graph", mesh; kind="cell")
```

```r
g <- mio_adjacency(mesh, "cell")   # list(indptr, indices), 1-based
mio_write_graph(mesh, "part.graph", kind = "cell")
```

```js
const { indptr, indices } = m.adjacency(mesh, 'clique');   // Int32Array, 0-based
m.writeGraph(mesh, '/part.graph', 'cell');                  // into MEMFS
```

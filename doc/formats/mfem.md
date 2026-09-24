# MFEM mesh and grid functions (`.mesh`, `.gf`)

[MFEM](https://mfem.org)'s own ASCII mesh, the file its examples write and [GLVis](https://glvis.org) opens, and the grid functions (`.gf`) that hold a field on it. A mesh lists `elements` and `boundary` elements by an integer **attribute** and a geometry code, then either `vertices` with their coordinates or, for a curved high-order mesh, a `nodes` grid function that gives the geometry.

| | |
|---|---|
| **Format name** | `mfem` |
| **Extensions** | `.mesh`, shared with [Medit](./medit.md): a file whose first line names an MFEM mesh is read as one |
| **Read / Write** | ✓ / ✓ (grid functions both ways) |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshioplusplus

mesh = meshioplusplus.read("star-q2.mesh")                          # an order-2 mesh
mesh = meshioplusplus.mfem.read("refined.mesh", {"u": "sol.gf"})     # with a field
mesh = meshioplusplus.mfem.read("refined.mesh", ["sol.gf"])          # named by its stem
mesh = meshioplusplus.read("fichera-q3.mesh")                       # order 3: VTK Lagrange cells
mesh = meshioplusplus.mfem.read("run.mesh.000000", {"u": "u.000000"})  # a parallel run, merged
meshioplusplus.mfem.write("out.mesh", mesh, grid_functions=True)     # + out.<name>.gf
```

Both engines (the C++ core and the pure-Python reference) read the same meshes and write the same bytes. In C++ the fields are `read_mfem(path, {{"u", "sol.gf"}})` and `write_mfem(path, mesh, true)`. The [MCP server](../mcp.md)'s `convert` takes the same `grid_functions` and `write_grid_functions` options.

**`.mesh` is shared with Medit.** The Python `read` tries Medit and then MFEM. The native resolver (CLI, C, Fortran, Julia, R, WASM) looks at the first line of an existing `.mesh` file and hands an `MFEM mesh …` file to `mfem`. It is the one content-aware extension default. Writing a `.mesh` without a format still writes Medit: pass `file_format="mfem"` (CLI `--output-format mfem`).

## Reading

| Version | Read as |
|---|---|
| `MFEM mesh v1.0` | a conforming mesh |
| `MFEM mesh v1.2` | the same, up to `mfem_mesh_end`; the rank file of a parallel mesh is read with its siblings (see [Parallel meshes](#parallel-meshes)) |
| `MFEM mesh v1.3` | v1.2 plus `attribute_sets` and `bdr_attribute_sets` |
| `MFEM mesh v1.1` | a conforming mesh; one with the legacy `vertex_parents` and `coarse_elements` sections is its leaf mesh, which is what it lists |
| `MFEM NC mesh v1.0`, `v1.1` | its leaf elements (see [Non-conforming meshes](#non-conforming-meshes)) |
| NURBS, INLINE | refused with a `ReadError` naming why |

- Elements and boundary elements are separate cell blocks, one per type, the elements first. Geometries 0–7 map to `vertex`, `line`, `triangle`, `quad`, `tetra`, `hexahedron`, `wedge` and `pyramid`. Every geometry, the prism included, is in meshio++'s node order.
- The attribute is the `mfem:attribute` cell data and an `attribute_<n>` (element) or `boundary_<n>` (boundary) `cell` region tagged `n`. A v1.3 attribute set is a `cell` region of its own name, holding the cells of all its attributes.
- `#` lines are comments. Boundary elements stay cells, like [Elmer](./elmer.md)'s, because MFEM allows boundary elements between two elements and interior ones, which a side region cannot express.

## Curved meshes

The `nodes` grid function names its space by a finite element collection:

| Collection | Read as |
|---|---|
| `H1_<d>D_P1`, `Linear`, `H1Pos_…_P1` | linear cells, vertices from the nodes |
| `H1_<d>D_P2`, `H1@…_P2`, `Quadratic` | **order-2 cells**: `line3`, `triangle6`, `quad9`, `tetra10`, `wedge18`, `hexahedron27` |
| `H1_<d>D_P<p>`, `H1@G`/`H1@U…_P<p>` with p ≥ 3, and the legacy `Cubic` (v16.11.0) | **VTK Lagrange cells** of order p (`VTK_LAGRANGE_CURVE`, `…_TRIANGLE`, `…_QUADRILATERAL`, `…_TETRAHEDRON`, `…_HEXAHEDRON`, `…_WEDGE`), their nodes on the equispaced lattice VTK uses |
| `L2_T1_<d>D_P1` (periodic meshes) | linear cells, **each element with its own points** |
| Bernstein (`H1Pos`) or serendipity (`H1Ser`) at order 2+, and other bases | **the vertices only**, with a warning |
| anything else | a `ReadError` |

An order-2 mesh's points are MFEM's degrees of freedom, **in MFEM's order**: the vertices, then the edges in order of first appearance (elements in file order, each element's edges in MFEM's local order), then the quadrilateral faces in 3-D (numbered over all faces the same way), then the element interiors. Each cell's nodes are placed by the vertices of the edge, face or cell they sit on, so edge and face orientation never matters at order 2. A mesh with pyramids is read at order 1, since meshio++ has no curved pyramid.

**Arbitrary order** (v16.11.0). MFEM's order-p nodes sit on Gauss–Lobatto points (`H1@G`, the default), closed-uniform points (`H1@U`) or, for `Cubic`, the equispaced points of MFEM's legacy cubic elements. meshio++ numbers MFEM's degrees of freedom the way MFEM does (vertices, edges, faces, interiors; each edge and face in the orientation of the element that first holds it), places every one at its reference position, and evaluates each cell's nodal basis at the VTK Lagrange points by interpolation, so the curved geometry is kept exactly (an order-p polynomial through order-p points) rather than resampled. Nodes shared by neighbouring cells are shared points. The result is what MFEM's own high-order VTK output (`ParaViewDataCollection` with `SetHighOrderOutput`) holds, point for point. meshio++ has no fixed type for these cells: they are the ragged `VTK_LAGRANGE_*` types, which [VTU](./vtu.md) and [VTK](./vtk.md) write and read.

## Non-conforming meshes

An `MFEM NC mesh` (v16.11.0) stores a refinement tree: root elements (`rank attr geom ref_type` and their children or vertices, `geom = -1` for an unused slot), `vertex_parents` (a vertex made on the edge between two others; v1.1 adds the position along it) and the `coordinates` of the root vertices. meshio++ reads the **leaf elements** and boundary elements, placing each new vertex from its parents; the hanging nodes are left as they are (a leaf's edge may end in the middle of its neighbour's edge), with a warning. Elements of other ranks (a parallel NC file) are dropped with a warning. Grid functions on an NC mesh are skipped with a warning: their degrees of freedom follow MFEM's space-filling-curve ordering of the leaves, which is not read.

## Parallel meshes

A parallel MFEM run saves one file per rank, `<prefix>.000000`, `<prefix>.000001`, … (v16.11.0). Opening any one of them reads **the whole mesh**: its siblings are found beside it, and shared vertices and higher-order nodes are merged, so every cell appears once, and `cell_data["partition:part"]` holds each cell's rank.

- `ParMesh::ParPrint` files (`MFEM mesh v1.2`, with `communication_groups` and the shared vertices, edges and faces of each group) are merged through their groups.
- `ParMesh::Save` files are plain serial meshes whose boundary also lists the faces shared with other ranks. Their boundary vertices are merged by position (to `1e-9` of the mesh size), and a boundary face that two ranks both list, being interior, is dropped.
- `piece=k` (`mfem.read(path, piece=k)`, `ReadOptions::mPiece` in C++) reads rank k alone.
- A rank's grid function is named by one of its files: `{"u": "sol.000000"}` reads `sol.000000`, `sol.000001`, … with the mesh.

Non-conforming parallel meshes are refused.

## Grid functions

| Space | Read as |
|---|---|
| `H1` order 1 or 2 | `point_data` (vectors as `(n, vdim)`); both `Ordering: 0` (byNODES) and `1` (byVDIM) |
| `H1` order 3 and up (v16.11.0) | `point_data` at the VTK Lagrange nodes, interpolated as the geometry is |
| `L2` order 0 | `cell_data`, element cells only; NaN on the boundary cells |
| anything else | skipped, with a warning |

An order-2 field on a linear mesh makes the cells quadratic, with the new nodes at the edge, face and cell centres, and an order-p field on a lower-order mesh makes them VTK Lagrange cells of order p. An order-1 field on an order-2 mesh is interpolated linearly to the extra nodes, which is exact. A field whose size does not match the mesh is a `ReadError`.

**A field must be saved with the mesh it lives on.** MFEM re-marks triangles and tetrahedra when it loads a mesh (it reorders their vertices for refinement), and a grid function numbers its degrees of freedom on the re-marked mesh. A `.gf` saved in a session therefore matches the mesh *printed* in that session (`mesh.Print`), not always the original file of a triangle or tetrahedral mesh; pair them as MFEM's own examples do.

## Writing

- The cells of the highest dimension are the **elements**, and those one dimension lower the **boundary**. Each facet of a `side` region becomes one more boundary element. Lower-dimensional cells are dropped with a warning.
- The attribute of each cell is `mfem:attribute`, else its cell region's tag (or a fresh number), else 1. A cell region other than `attribute_<n>`/`boundary_<n>` becomes an attribute set, which makes the file `MFEM mesh v1.3`; otherwise it is v1.0.
- Quadratic cells make an `H1_<d>D_P2` `nodes` space. `quad8`, `hexahedron20` and `wedge15` are completed with their face and body centres from the serendipity map, which keeps them curved. Linear cells in a quadratic mesh get their extra nodes at the centres.
- VTK Lagrange cells of order p (v16.11.0) make an `H1_<d>D_P<p>` space on Gauss–Lobatto points, their nodal values interpolated from the VTK lattice; MFEM reads such a file back to the same geometry. Lower-order cells in it are raised to order p exactly; serendipity cells are placed from their corners, with a warning. A mesh with pyramids is written with corners only, with a warning.
- Point coordinates and field values are written with 17 significant digits.
- `grid_functions=True` writes one `<stem>.<name>.gf` per data array: point data as `H1` at the mesh's order and cell data as `L2` order 0. Without it, and for field data always, **data is dropped** with a warning and a provenance note. Point regions are dropped the same way (MFEM has no node sets).

## Verification

Reading and writing were checked against MFEM 4.10 itself through PyMFEM, with MFEM's own sample meshes (BSD-3, the test fixtures):

- For every element of 14 samples (orders 1–3, triangles to pyramids, curved and mixed), every node meshio++ reads sits where MFEM's element transformation puts it, to `0.0`.
- MFEM reads the files meshio++ writes to the same geometry and field values, to `1.7e-16`.
- Arbitrary order (v16.11.0): on 14 meshes of orders 3–5 (MFEM's `fichera-q3`, `escher-p3`, `toroid-wedge` and `rt-2d-p4-tri`, and one-element reference meshes of every shape curved and warped by MFEM, Gauss–Lobatto and closed-uniform), every Lagrange node and field value equals MFEM's own high-order VTK output to `1.1e-13`, and MFEM reads the order-p files meshio++ writes to `1.4e-13`.
- Non-conforming (v16.11.0): the leaf and boundary elements of `amr-quad` match the ones MFEM builds, attribute and corners.
- Parallel (v16.11.0): meshes MFEM saved over 3 and 4 ranks with MPI, in both layouts, read to MFEM's own per-rank high-order output, with exactly the serial mesh's vertex count once merged.
- `tools/gen_mfem_fixtures.py` and `tools/gen_mfem_parallel_fixtures.py` freeze MFEM's evaluation, so the tests check every node and field value without MFEM.

## Notes

- **MFEM's own VTK export reverses prisms.** Its `PrismMap` exists because classic VTK winds the wedge the other way round from MFEM and meshio++. A prism read from a `.mesh` and one from an MFEM-written `.vtk` therefore differ in node order, not in shape.
- **NURBS meshes** (`MFEM NURBS mesh`) are still refused: their geometry is a set of patches with global control points and knot vectors, which would have to be evaluated to Lagrange cells. That stays in the [roadmap](../roadmap.md).

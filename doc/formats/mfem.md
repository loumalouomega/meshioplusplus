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
meshioplusplus.mfem.write("out.mesh", mesh, grid_functions=True)     # + out.<name>.gf
```

Both engines (the C++ core and the pure-Python reference) read the same meshes and write the same bytes. In C++ the fields are `read_mfem(path, {{"u", "sol.gf"}})` and `write_mfem(path, mesh, true)`. The [MCP server](../mcp.md)'s `convert` takes the same `grid_functions` and `write_grid_functions` options.

**`.mesh` is shared with Medit.** The Python `read` tries Medit and then MFEM. The native resolver (CLI, C, Fortran, Julia, R, WASM) looks at the first line of an existing `.mesh` file and hands an `MFEM mesh …` file to `mfem`. It is the one content-aware extension default. Writing a `.mesh` without a format still writes Medit: pass `file_format="mfem"` (CLI `--output-format mfem`).

## Reading

| Version | Read as |
|---|---|
| `MFEM mesh v1.0` | a conforming mesh |
| `MFEM mesh v1.2` | the same, up to `mfem_mesh_end`; a rank file of a parallel mesh is read as its local part, with a warning that its communication groups are ignored |
| `MFEM mesh v1.3` | v1.2 plus `attribute_sets` and `bdr_attribute_sets` |
| `MFEM mesh v1.1` | a conforming mesh unless it has a `vertex_parents` section |
| `MFEM NC mesh`, `vertex_parents`, NURBS, INLINE | refused with a `ReadError` naming why |

- Elements and boundary elements are separate cell blocks, one per type, the elements first. Geometries 0–7 map to `vertex`, `line`, `triangle`, `quad`, `tetra`, `hexahedron`, `wedge` and `pyramid`. Every geometry, the prism included, is in meshio++'s node order.
- The attribute is the `mfem:attribute` cell data and an `attribute_<n>` (element) or `boundary_<n>` (boundary) `cell` region tagged `n`. A v1.3 attribute set is a `cell` region of its own name, holding the cells of all its attributes.
- `#` lines are comments. Boundary elements stay cells, like [Elmer](./elmer.md)'s, because MFEM allows boundary elements between two elements and interior ones, which a side region cannot express.

## Curved meshes

The `nodes` grid function names its space by a finite element collection:

| Collection | Read as |
|---|---|
| `H1_<d>D_P1`, `Linear`, `H1Pos_…_P1` | linear cells, vertices from the nodes |
| `H1_<d>D_P2`, `H1@…_P2`, `Quadratic` | **order-2 cells**: `line3`, `triangle6`, `quad9`, `tetra10`, `wedge18`, `hexahedron27` |
| `L2_T1_<d>D_P1` (periodic meshes) | linear cells, **each element with its own points** |
| higher orders, `Cubic`, Bernstein (`H1Pos`) or serendipity (`H1Ser`) at order 2+ | **the vertices only**, with a warning; see the note below |
| anything else | a `ReadError` |

An order-2 mesh's points are MFEM's degrees of freedom, **in MFEM's order**: the vertices, then the edges in order of first appearance (elements in file order, each element's edges in MFEM's local order), then the quadrilateral faces in 3-D (numbered over all faces the same way), then the element interiors. Each cell's nodes are placed by the vertices of the edge, face or cell they sit on, so edge and face orientation never matters at order 2. A mesh with pyramids is read at order 1, since MFEM has no order-2 pyramid meshio++ can hold.

## Grid functions

| Space | Read as |
|---|---|
| `H1` order 1 or 2 | `point_data` (vectors as `(n, vdim)`); both `Ordering: 0` (byNODES) and `1` (byVDIM) |
| `L2` order 0 | `cell_data`, element cells only; NaN on the boundary cells |
| anything else | skipped, with a warning |

An order-2 field on a linear mesh makes the cells quadratic, with the new nodes at the edge, face and cell centres. An order-1 field on an order-2 mesh is interpolated linearly to the extra nodes, which is exact. A field whose size does not match the mesh is a `ReadError`.

## Writing

- The cells of the highest dimension are the **elements**, and those one dimension lower the **boundary**. Each facet of a `side` region becomes one more boundary element. Lower-dimensional cells are dropped with a warning.
- The attribute of each cell is `mfem:attribute`, else its cell region's tag (or a fresh number), else 1. A cell region other than `attribute_<n>`/`boundary_<n>` becomes an attribute set, which makes the file `MFEM mesh v1.3`; otherwise it is v1.0.
- Quadratic cells make an `H1_<d>D_P2` `nodes` space. `quad8`, `hexahedron20` and `wedge15` are completed with their face and body centres from the serendipity map, which keeps them curved. Linear cells in a quadratic mesh get their extra nodes at the centres. A mesh with pyramids is written linear, with a warning.
- Point coordinates and field values are written with 17 significant digits.
- `grid_functions=True` writes one `<stem>.<name>.gf` per data array: point data as `H1` at the mesh's order and cell data as `L2` order 0. Without it, and for field data always, **data is dropped** with a warning and a provenance note. Point regions are dropped the same way (MFEM has no node sets).

## Verification

Reading and writing were checked against MFEM 4.10 itself through PyMFEM, with MFEM's own sample meshes (BSD-3, the test fixtures):

- For every element of 14 samples (orders 1–3, triangles to pyramids, curved and mixed), every node meshio++ reads sits where MFEM's element transformation puts it, to `0.0`.
- MFEM reads the files meshio++ writes to the same geometry and field values, to `1.7e-16`.
- `tools/gen_mfem_fixtures.py` freezes MFEM's evaluation, so the tests check every node and field value without MFEM.

## Notes

- **MFEM's own VTK export reverses prisms.** Its `PrismMap` exists because classic VTK winds the wedge the other way round from MFEM and meshio++. A prism read from a `.mesh` and one from an MFEM-written `.vtk` therefore differ in node order, not in shape.
- **Arbitrary order.** meshio++'s fixed cell types stop at order 2 for these shapes. An order-3+ MFEM mesh (Gauss–Lobatto points) keeps only its vertices, and its curvature is lost. That, non-conforming meshes, and parallel meshes as a whole stay in the [roadmap](../roadmap.md).

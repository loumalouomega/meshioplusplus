# Polyhedra and ragged cells

Most cell types have a fixed node count: a `tetra` always has four nodes, a `hexahedron` always eight. A **ragged** block does not — its cells vary in size, which is what a Voronoi polygon or a general polyhedral CFD cell needs. meshio++ stores two ragged shapes, and this page is the single place their rules are stated, because five surfaces and four language bindings all have to agree on them.

## The two shapes

**1-level ragged** — a `polygon*` block. Each cell is one list of node ids of whatever length. MED `POG`, EnSight `nsided` and VTP `Polys` produce these.

**2-level ragged** — a `polyhedron*` block. Each cell is a list of **faces**, each face a list of node ids. There is no separate cell-node list: a cell's nodes are whatever its faces reference. OpenFOAM, EnSight `nfaced` and VTU's `VTK_POLYHEDRON` produce these.

The type name encodes the node count where the format supplies one — `polyhedron4` has four distinct nodes, `polyhedron12` twelve — which is how readers group cells into blocks. See [Cell types](/cell_types).

::: warning Check the type name, not `IsRagged()`
A block whose cells all *happen* to have the same node count may be stored rectangularly and is then not ragged at all, even though its type is `polygon`. Code that must reject polyhedra (a format writer, say) has to test the **type name**, not `IsRagged()`.
:::

## The flat CSR vocabulary

Nested lists do not cross a flat ABI, so every non-C++ surface carries ragged blocks as flat CSR arrays. There is **one** set of names:

| Array | Length | Meaning |
| --- | --- | --- |
| `nodes` (JS: `data`) | `num_nodes` | every face's node ids concatenated, 0-based |
| `face_offsets` (JS: `faceOffsets`, or `rowOffsets` when 1-level) | `num_faces + 1` | each face's start index into `nodes` |
| `cell_offsets` (JS: `cellOffsets`) | `num_cells + 1` | each cell's start index into the **face list** |

So face `f` of cell `c` is

```
nodes[ face_offsets[ cell_offsets[c] + f ] .. face_offsets[ cell_offsets[c] + f + 1 ] )
```

A 1-level block has exactly one face per cell, so `face_offsets` *is* the row-offsets array and there is no `cell_offsets` at all — the C API returns `NULL` and JS omits the key, rather than synthesising a `0,1,2,…` identity that a caller could mistake for information.

::: tip A naming skew worth knowing
The `NATIVE` mesh backend stores the same two arrays under the opposite names: its `mRowOffsets` is this page's `face_offsets`, and its `mFaceOffsets` is this page's `cell_offsets`. The CSR vocabulary above is what every *binding* uses; `native_mesh.hpp`'s field names are internal and predate it.
:::

![A two-level ragged block as the three CSR arrays cell_offsets, face_offsets and nodes, with the indexing of one face traced through them](/diagrams/polyhedra_csr.svg)

## Winding

Faces **should** be wound so that the right-hand normal points *out* of the cell. That is what OpenFOAM's `polyMesh` gives on read, what VTK's type 42 asks for, and what meshio++ writes.

meshio++ does **not require** it. Real meshes arrive with faces wound inconsistently — including, as it happens, this repository's own long-standing polyhedral test fixture — and rejecting them would be less useful than fixing them. Every geometric kernel therefore repairs the winding per cell before measuring, and reports the face set as unorientable only when it is genuinely not a closed orientable surface (some undirected edge is not used exactly twice).

The practical consequence: a cell volume never comes back negative because someone's faces were inside-out.

## Non-planar faces

A face with four or more corners that are not coplanar does not bound a unique volume — the answer depends on how you triangulate it. meshio++'s answer is the **corner-average fan**: the face is the union of triangles from each of its edges to the face's corner average.

That is chosen, not arbitrary:

- It is **independent of which corner the face's node list starts at**. A fan about the first node is not, and OpenFOAM's `faces` file, MED's `INN` and VTU's `faces` stream all pick that start node arbitrarily — so two readers of the same mesh could otherwise legitimately disagree about a volume.
- It is the same surface [`gradient`](/gradient)'s Green–Gauss integration already uses, for an independent reason (the corner average is the only apex whose value is known exactly for a linear field).
- It is OpenFOAM's own cell-volume decomposition, so a round-tripped case reports the volumes its solver would.

## Support across the surfaces

| Surface | Build a ragged block | Read one back |
| --- | --- | --- |
| C++ | `AddPolygonBlock` / `AddPolyhedronBlock` | `CellView::Row`/`RowSize`, `NumFaces`/`Face` |
| Python | `Mesh(..., cells=[("polyhedron4", [[[...]]])])` | `CellBlock.data` (a nested list) |
| [C API](/c_api) | `mio_mesh_add_polygon_block` / `_polyhedron_block` | `mio_poly_conn_create` + the three array accessors |
| [Fortran](/fortran) | `m%add_polygon_block` / `m%add_polyhedron_block` | `m%get_polygon_block` / `m%get_polyhedron_block` (CSR, 1-based) |
| [Julia](/julia) | `add_polygon_block!` / `add_polyhedron_block!` | `polygon_block` / `polyhedron_block` (nested vectors, 1-based) |
| [R](/r) | `mio_add_polygon_block()` / `mio_add_polyhedron_block()` | `mio_polygon_block()` / `mio_polyhedron_block()` (nested lists, 1-based) |
| [WASM](/wasm) | `{type, data, rowOffsets}` / `{type, data, faceOffsets, cellOffsets}` on the mesh object | the same shape back |

Each language gets its **natural** shape over one flat ABI: Fortran has no ragged array type so it gets the CSR triple; Julia has nested vectors and R has lists, so those are what they hand back. This is the same policy that already gives all three column-major arrays and 1-based indices.

### The C API's snapshot is not a borrow

Every other getter in the C API is a zero-copy borrow that dies at the next mutating call ([rule 3](/c_api)). `mio_poly_conn` is an **owning snapshot** instead, and deliberately so: the `MESHIO` mesh backend stores ragged blocks as nested vectors, so there is no offsets array inside the mesh to point at, and the `KRATOS` backend rebuilds its blocks lazily. The snapshot is therefore *safer* than a borrow — it stays valid across mutating calls — and must be released with `mio_poly_conn_free()`.

Julia's `polygon_block`/`polyhedron_block` copy out of that snapshot and free it, so there is no `_ptr` form for ragged blocks: a `MeshBorrow` carries the mesh's mutation generation as its guard, which is the wrong guard here.

## Format support

Reading polyhedra: **MED** (`POE`), **EnSight** (`nfaced`), **OpenFOAM**, **VTU** (type 42) and **CGNS** (`NGON_n`+`NFACE_n`). CGNS additionally reads ADF containers and the CGNS 3.x section layout on a build with the optional [cgnslib backend](/formats/cgns). Reading jagged polygons: **MED** (`POG`/`POG2`), **EnSight** (`nsided`), **VTP**, **OpenFOAM**.

Writing polyhedra: **MED** (`POE`), **EnSight** (`nfaced`) and **VTU** (type 42) since v9.19.0; **OpenFOAM** since v9.20.0, where they are the native cell shape; **CGNS** (`NGON_n`+`NFACE_n`) since v9.21.0. Writing jagged polygons: **MED** (`POG`/`POG2`), **EnSight** (`nsided`), **VTP**, **CGNS** (`NGON_n`).

See [Formats](/formats) for the current table and each format page for what it can express. A writer that cannot represent a ragged block fails naming the format rather than silently dropping cells.

## Operations

**Structural** operations carry ragged blocks through unchanged: `crop`, `split`, `clean`, `merge`, `transform`, `reorder`, `diff` and the data operations.

**Geometric** operations go through the kernel (`detail/polyhedron.hpp`) as of v9.16.0:

| Operation | On a polyhedron |
| --- | --- |
| [`stats`](/stats) | volume, area and centroid, via the corner-average fan. A cell that is not a closed orientable surface is excluded and warned about, not silently counted as zero |
| [`gradient`](/gradient) | fully supported — Green–Gauss needs only faces, so a polyhedron goes through the same code as a hexahedron and recovers a linear field exactly |
| [`compute_quality`](/mesh_quality) | a **reduced** set: `volume`, `inverted`, `degenerate`. Every metric defined against a reference element stays `NaN` — see below |
| [`extract_surface`](/extract_surface) / `extract_skin` | fully supported, including faces of any arity: a pentagon becomes a row of a ragged `polygon` block |
| [`smooth`](/smooth) | fully supported: the inversion guard and `fix_boundary` both cover them |
| [`partition`](/partition) | fully supported, both methods: SFC works off centroids, and the KaHIP dual graph connects cells sharing a face |
| [`clean`](/clean) | degenerate and duplicate polyhedra are dropped; the duplicate key is the sorted set of sorted faces |
| [`data_average`](/data_average) | `Measure` weighting works, instead of falling back to a unit weight |

### Why `quality` reports so little

Aspect ratio, skewness, warpage and the angle metrics are all defined against a *reference element* — an ideal tetrahedron, an ideal hexahedron. A polyhedron has no reference element, so those numbers would be invented rather than measured. meshio++ reports `NaN`, which is its standing convention for a metric that does not apply, rather than something plausible-looking.

"Inverted" also means something different here: since winding is repaired rather than required, there is no convention for a cell to violate. What a polyhedron *can* fail to be is a closed orientable surface, and that is what the flag reports.

### Not yet supported

- `refine` and `decimate` raise by name — both are built on fixed subdivision templates, and an arbitrary polyhedron has none. Both messages point at `convert_cells(simplexify)`, which since v9.17.0 decomposes a polyhedron into tetrahedra using the *same* fan the kernel measures — so the decomposition conserves volume exactly. That also makes [`slice`](/slice), [`isosurface`](/isosurface) and `interpolate --barycentric` work on polyhedra, since all three go through it.
- The pure-numpy `convert_cells` fallback is rectangular-only, so on a build with no compiled core a polyhedral simplexify raises `NotImplementedError`.

One shared facet key (`detail::FacetKey`, variable arity with a four-entry inline fast path) backs `extract_surface`, `smooth`'s boundary marking and `partition`'s dual graph. That matters beyond tidiness: with a separate key per cell shape, a hexahedron and a polyhedron meeting on a face would not recognise each other, and both would report that interior face as boundary.

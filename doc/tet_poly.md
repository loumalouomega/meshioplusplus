# Tetrahedral ⇄ polyhedral conversion (tet_to_poly / poly_to_tet)

`meshioplusplus.tet_to_poly(mesh)` converts a tetrahedral mesh into a polyhedral one. It builds the **median dual**, the tet-to-poly conversion CFD meshers use: one polyhedron per vertex, typically four to six times fewer cells than tetrahedra, with many more neighbours each. `meshioplusplus.poly_to_tet(mesh)` goes the other way and decomposes a polyhedral (or mixed) mesh into tetrahedra. Both are mesh **operations**, not file formats. They use only standard C++ and run under every mesh backend.

![The polyhedral dual of a tetrahedral ball, cut away and shrunk, coloured by the number of faces per cell](/images/tet_to_poly.png)

Polyhedral meshes themselves are first-class in meshio++: see [polyhedra and ragged cells](/polyhedra) for how they are stored and which formats read and write them.

```python
import meshioplusplus as mp

tets = mp.read("domain.msh")
poly = mp.tet_to_poly(tets)  # one polyhedron per vertex
mp.write("domain_poly.vtu", poly)

back = mp.poly_to_tet(poly)  # tetrahedra again (not the original ones)
```

## tet_to_poly: the median dual

Every vertex `v` of the tetrahedral mesh becomes one polyhedral cell, bounded by:

- **One interior face per edge `(v, w)` leaving the vertex.** It is the ring of the tetrahedron centroids and triangle centroids around that edge, wound so its normal points from `v` to `w`. The cell of `w` gets the same ring reversed, so the result is conforming. Where the edge reaches the boundary the ring is open, and it closes through the edge midpoint.
- **Boundary faces, where the vertex lies on the boundary.** Each boundary triangle `(v, a, b)` contributes the planar quadrilateral `v, mid(v,a), centroid, mid(v,b)`, its corner of the triangle. Around `v` these quadrilaterals are fused into one polygon per **run** of consecutive triangles.

A run breaks where two consecutive boundary triangles' normals differ by more than `feature_angle` degrees (default 1), where the side changes (the mesh boundary versus an interface with a particular other block), or where the input triangle block changes. A run that goes all the way round drops `v`, which then lies inside the polygon. Any other run keeps `v` as a corner, so feature edges and corners survive as edges and vertices of the dual.

With the default angle every boundary face is planar, and the dual's total volume equals the tetrahedral mesh's up to rounding. A larger angle gives fewer boundary faces, but a fused face over a curved boundary is not planar, and the volume then differs slightly.

```python
poly, report = mp.tet_to_poly(tets, feature_angle=30.0, return_report=True)
report["num_boundary_faces"], report["num_new_points"], report["num_split_vertices"]
```

### Interfaces and pinched vertices

A dual cell is built per **(vertex, `tetra` block, face-connected component of that block's tetrahedra around the vertex)**. A vertex on the interface between two blocks therefore yields one cell on each side, so a material interface stays sharp, and the two cells meet on a conforming interface face. A pinched vertex, where tetrahedra touch only at the vertex, yields one cell per lobe. `num_split_vertices` counts the vertices that gave more than one cell.

### Output structure

- **Points.** The input points keep their indices. An interior vertex ends up unreferenced; [`clean`](/clean) with `remove_orphans=True` compacts the point set. New points follow in a fixed order: the midpoints of the edges that reach a boundary (in sorted edge order), the centroid of every triangle face, then the centroid of every tetrahedron. `point_data` is extended to the new points by averaging the nodes each one is built from.
- **Blocks.** Each `tetra` block becomes one `polyhedron` block, in the same position. Each `triangle` block becomes a ragged `polygon` block holding the dual boundary faces over its triangles; on an interface, the faces from the lower-index block's side are used. A triangle that is not on a boundary is dropped with a warning, and so is every other block (lines, vertices, and so on).
- **`cell_data`.** A polyhedron's value is the |volume|-weighted mean over its tetrahedra for floating-point arrays, and the value of its largest tetrahedron (ties: the lowest index) for integer arrays. Polygon blocks do the same with triangle areas. An array that does not have one row per cell of every block is dropped with a warning. `field_data` passes through.
- **`point_data_to_cell_data=True`** also copies each vertex's `point_data` onto its dual cells and boundary polygons as `cell_data`, the natural map from a nodal field to a cell-centred one. A name already used by a `cell_data` array is skipped with a warning.
- **`record_vertex_ids=True`** attaches an Int64 `dual:vertex` `cell_data` array holding each output cell's source vertex. The report's `vertex_map` holds the same information as one flat array over the output cells.
- **Regions.** Point regions survive, since points are not renumbered. Cell and side regions are kept as empty groups with a warning, because a dual cell belongs to a vertex, not to a tetrahedron.

### What it refuses

A volume block that is not linear `tetra` raises `ValueError` naming the type. For `tetra10`, linearize with [`convert_cells`](/convert_cells) first. For hexahedra, wedges, pyramids or existing polyhedra, decompose with `poly_to_tet` first. It also refuses a mesh with no tetrahedra, with points that are not 3D, or with a face shared by three or more tetrahedra (non-manifold).

## poly_to_tet

`poly_to_tet` is [`convert_cells(mode="simplexify")`](/convert_cells) with one addition: a polyhedron that already is a tetrahedron (four triangular faces over four distinct nodes, closed and orientable) becomes that one positively oriented `tetra`, with no new points. Every other cell takes simplexify's own path:

- a polyhedron uses the volume-conserving corner-average fan, one tetrahedron per face edge, adding the cell's and each face's corner average as points ([polyhedra](/polyhedra#non-planar-faces));
- a hexahedron, wedge or pyramid uses the fixed templates;
- a polygon uses a triangle fan.

A tetrahedral mesh stored as `polyhedron4` cells, as some readers produce, therefore comes back exactly, whatever the winding of its faces. `poly_to_tet(tet_to_poly(mesh))` is a valid, volume-conserving tetrahedral mesh, but it is not `mesh`: the dual has forgotten the original tetrahedra.

```python
tets, report = mp.poly_to_tet(poly, return_report=True)
report["num_exact_tets"], report["num_fanned_cells"]
```

`record_parent_ids=True` attaches the `convert:parent_cell` array, and the report's `point_map` and `cell_maps` have `convert_cells`' shape.

Neither operation has a pure-numpy fallback. Both depend on winding and orientation decisions that a second implementation could resolve differently near degenerate cells, so without the compiled core they raise `NotImplementedError`.

## Other surfaces

| Surface | `tet_to_poly` | `poly_to_tet` |
| --- | --- | --- |
| CLI | `meshioplusplus tet-to-poly IN OUT [--feature-angle DEG] [--point-data-to-cell-data] [--record-vertex-ids] [--json]` | `meshioplusplus poly-to-tet IN OUT [--record-parent-ids] [--json]` |
| [Pipeline](/pipeline) | `TetToPoly` | `PolyToTet` |
| [C API](/c_api) | `mio_tet_to_poly` (`mio_tet_to_poly_opts`) | `mio_poly_to_tet` |
| [Fortran](/fortran) | `m%tet_to_poly(...)` | `m%poly_to_tet(...)` |
| [Julia](/julia) | `tet_to_poly(mesh; ...)` | `poly_to_tet(mesh; ...)` |
| [R](/r) | `mio_tet_to_poly()` | `mio_poly_to_tet()` |
| [WASM](/wasm) | `tetToPoly(mesh, options, returnMaps)` | `polyToTet(mesh, recordParentIds, returnMaps)` |
| [MCP](/mcp) | `tet_to_poly` | `poly_to_tet` |

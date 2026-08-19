# Surface / boundary extraction

`meshioplusplus.extract_surface(mesh, record_parent_ids=False)` derives the **boundary of a mesh's highest-dimension cells** as a new `Mesh`. It is the general form of [skin extraction](./extract_skin.md): the same face-hashing algorithm (a facet whose sorted corner-node key occurs **exactly once** is on the boundary), but it picks the dimension automatically and also handles 2D surface meshes.

![The mesh's boundary, extracted as a new surface mesh](/images/extract_surface.png)

```python
import meshioplusplus

vol = meshioplusplus.read("part.msh")           # tetra/hexa/wedge/pyramid mesh
surf = meshioplusplus.extract_surface(vol)      # triangle/quad surface mesh
meshioplusplus.write("surface.stl", surf)

sheet = meshioplusplus.read("sheet.vtu")        # triangle/quad surface mesh
edges = meshioplusplus.extract_surface(sheet)   # line/line3 boundary loop
```

- If the mesh has any **3D volume cells**, the boundary **faces** are returned (`triangle`/`quad` and the quadratic `triangle6`/`quad8`/`quad9`).
- Otherwise, for a **2D surface mesh**, the boundary **edges** are returned (`line`/`line3`).

## Supported cell types

**Volume → faces** — `tetra`, `tetra10`, `hexahedron`, `hexahedron20`, `hexahedron27`, `wedge`, `wedge15`, `pyramid`, `pyramid13`, `pyramid14` (same face tables as [skin extraction](./extract_skin.md#supported-cell-types)).

**Surface → edges**

| surface type | boundary edges |
|---|---|
| `triangle` | 3 × `line` |
| `triangle6` | 3 × `line3` |
| `quad` | 4 × `line` |
| `quad8` / `quad9` | 4 × `line3` |

Facet matching uses the **corner nodes only**, so mixed linear/quadratic meshes still cancel shared facets correctly.

## Semantics

- **Output block order** is fixed: `triangle`, `triangle6`, `quad`, `quad8`, `quad9` for faces; `line`, `line3` for edges (only non-empty blocks appear), in cell-block enumeration order within each block.
- **Points are compacted** to the referenced subset in ascending original-index order; `point_data` rows are subset the same way. `field_data`, sets, and `info` are dropped.
- **`record_parent_ids=True`** attaches an `int64` `cell_data` array `"surface:parent_cell"` giving, for each output facet, the global index of the unique input cell that owns it (counted block-major over every input cell of every block). Point compaction does not affect it — it is per-output-cell.
- **Polyhedron blocks are supported** since v9.16.0: a cell's own faces are hashed like any others, so a polyhedron and a hexahedron meeting on a face cancel each other out, and a face of any arity is emitted (one with more than four corners becomes a row of a ragged `polygon` block).
- **Unsupported same-dimension blocks** (ragged polygon blocks, Lagrange and very-high-order types) are skipped with a warning. If no supported 2D/3D block exists, `ValueError` is raised.

## Relationship to `extract_skin`

`extract_surface` and [`extract_skin`](./extract_skin.md) share one implementation (`src/cpp/src/operations/surface.cpp`). Use `extract_skin` when you specifically want the volume skin with the optional `linearize` flag (what the STL/PLY/SVG/TikZ writers use); use `extract_surface` for the auto-dimension behavior, 2D-edge extraction, and parent-cell provenance.

## Cross-language

Available in every binding surface: Python (`extract_surface`), the C API (`mio_extract_surface`), Fortran (`mesh%extract_surface()`), and WASM (`extractSurface`), plus the CLI verb `meshioplusplus extract-surface`.

## Implementation

The extractor runs over the uniform mesh API (so it works under all three mesh backends) with a parallel facet-key phase and a serial dedup/emit phase that keeps output byte-identical across backends and thread counts. Face topology comes from `src/cpp/include/meshioplusplus/detail/cell_faces.hpp` and edge topology from `detail/cell_edges.hpp`; a pure-numpy twin (`meshioplusplus._surface._extract_surface_py`) produces identical output, asserted by `tests/python/test_extract_surface.py`.

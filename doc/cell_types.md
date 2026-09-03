# Cell Types

meshio++ uses its own canonical type names. Every format reader maps native element names to these; every writer maps them back.

Node ordering follows the VTK convention where available. The figures below are generated from the code's own topology tables (`_convert_cells._ELEVATE` for mid-edge nodes, `_skin._CELL_FACES` and `_refine_templates.QUAD_FACES` for face and body centres), so they cannot drift from what the readers, writers and operations actually assume; corner nodes are drawn solid, mid-edge nodes blue, face centres orange and the body centre violet. One convention is worth stating: `pyramid14` is `pyramid13` plus a base-centre node 13, the order every consumer in the core (`detail/cell_faces.hpp`, the CGNS table) uses, while the gmsh reader applies no permutation to gmsh's own edge-lexicographic type 14, a gap recorded in the [roadmap](./roadmap.md#9-format-completeness).

## 0-D

| Type | Nodes |
|------|-------|
| `vertex` | 1 |

## 1-D (line elements)

![Node ordering of vertex, line and line3](/diagrams/cell_types_line.svg)

| Type | Nodes | Description |
|------|-------|-------------|
| `line` | 2 | Linear line |
| `line3` | 3 | Quadratic line |
| `line4` | 4 | Cubic line |
| `line5` | 5 | Quartic line |
| `line6` | 6 | |
| `line7` | 7 | |
| `line8` | 8 | |
| `line9` | 9 | |
| `line10` | 10 | |
| `line11` | 11 | |

## 2-D surface elements

### Triangles

![Node ordering of triangle and triangle6](/diagrams/cell_types_triangle.svg)

| Type | Nodes |
|------|-------|
| `triangle` | 3 |
| `triangle6` | 6 |
| `triangle10` | 10 |
| `triangle15` | 15 |
| `triangle21` | 21 |
| `triangle28` | 28 |
| `triangle36` | 36 |
| `triangle45` | 45 |
| `triangle55` | 55 |
| `triangle66` | 66 |

### Quadrilaterals

![Node ordering of quad, quad8 and quad9](/diagrams/cell_types_quad.svg)

| Type | Nodes |
|------|-------|
| `quad` | 4 |
| `quad8` | 8 |
| `quad9` | 9 |
| `quad16` | 16 |
| `quad25` | 25 |
| `quad36` | 36 |
| `quad49` | 49 |
| `quad64` | 64 |
| `quad81` | 81 |
| `quad100` | 100 |
| `quad121` | 121 |

### Arbitrary polygons

| Type | Nodes |
|------|-------|
| `polygon` | variable |

For `polygon` cells, all cells within a `CellBlock` must have the same number of nodes (the array shape determines this). Multiple `CellBlock` entries with type `polygon` but different node counts are permitted.

## 3-D volume elements

### Tetrahedra

![Node ordering of tetra and tetra10](/diagrams/cell_types_tetra.svg)

| Type | Nodes |
|------|-------|
| `tetra` | 4 |
| `tetra10` | 10 |
| `tetra20` | 20 |
| `tetra35` | 35 |
| `tetra56` | 56 |
| `tetra84` | 84 |
| `tetra120` | 120 |
| `tetra165` | 165 |
| `tetra220` | 220 |
| `tetra286` | 286 |

### Hexahedra

![Node ordering of hexahedron, hexahedron20, hexahedron24 and hexahedron27](/diagrams/cell_types_hexahedron.svg)

| Type | Nodes |
|------|-------|
| `hexahedron` | 8 |
| `hexahedron20` | 20 |
| `hexahedron24` | 24 |
| `hexahedron27` | 27 |
| `hexahedron64` | 64 |
| `hexahedron125` | 125 |
| `hexahedron216` | 216 |
| `hexahedron343` | 343 |
| `hexahedron512` | 512 |
| `hexahedron729` | 729 |
| `hexahedron1000` | 1000 |
| `hexahedron1331` | 1331 |

### Wedges (prisms)

![Node ordering of wedge, wedge15 and wedge18](/diagrams/cell_types_wedge.svg)

| Type | Nodes |
|------|-------|
| `wedge` | 6 |
| `wedge15` | 15 |
| `wedge18` | 18 |
| `wedge40` | 40 |
| `wedge75` | 75 |
| `wedge126` | 126 |
| `wedge196` | 196 |
| `wedge288` | 288 |
| `wedge405` | 405 |
| `wedge550` | 550 |

### Pyramids

![Node ordering of pyramid, pyramid13 and pyramid14](/diagrams/cell_types_pyramid.svg)

| Type | Nodes |
|------|-------|
| `pyramid` | 5 |
| `pyramid13` | 13 |
| `pyramid14` | 14 |

### Arbitrary polyhedra

| Type | Description |
|------|-------------|
| `polyhedron4` | Tetrahedron as polyhedron |
| `polyhedron5` | Pyramid as polyhedron |
| `polyhedronnN` | N-faced polyhedron |

Polyhedron cells store their face connectivity as a list of lists (ragged), not a rectangular numpy array. The type name encodes the number of faces: `polyhedron4` has 4 faces.

## VTK Lagrange types

Higher-order VTK Lagrange elements (arbitrary polynomial order, controlled at runtime):

| Type |
|------|
| `VTK_LAGRANGE_CURVE` |
| `VTK_LAGRANGE_TRIANGLE` |
| `VTK_LAGRANGE_QUADRILATERAL` |
| `VTK_LAGRANGE_TETRAHEDRON` |
| `VTK_LAGRANGE_HEXAHEDRON` |
| `VTK_LAGRANGE_WEDGE` |
| `VTK_LAGRANGE_PYRAMID` |

These are read and written by the VTK/VTU readers/writers only.

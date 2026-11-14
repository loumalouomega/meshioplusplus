# 3MF — 3D Manufacturing Format (`.3mf`)

The [3MF](https://3mf.io) core specification: the package slicers exchange, with units, several objects and colour where STL has none (roadmap 2.4.5). Consumers: PrusaSlicer, Cura, Bambu Studio. Both engines read and write it (`formats/threemf.cpp` on the shared [ZIP container](../common_utilities.md), `meshioplusplus.threemf`), and they write the same bytes. The format is named `3mf`; its Python package is `meshioplusplus.threemf`, since a module name cannot start with a digit.

| | |
|---|---|
| **Format name** | `3mf` |
| **Extensions** | `.3mf` |
| **Read / Write** | ✓ / ✓ (triangles) |
| **Extra dependencies** | — (zlib for deflated members) |

## Reading

```python
import meshioplusplus

part = meshioplusplus.read("plate.3mf")
```

- The package's `_rels/.rels` names the model part (`3D/3dmodel.model` when it does not).
- Every `<build>` `<item>` is instantiated: an object's `<mesh>` (`<vertices>`, `<triangles>`) or its `<components>`, with each `transform` (the specification's 3x4 row-vector matrix) composed down the tree. An object built twice appears twice.
- All triangles are one block; each build item is a cell region named after its object (`object_<id>` when it has no name; a repeated name gets `_<k>`).
- The model's `unit` (micron, millimeter, centimeter, inch, foot, meter; millimetre by default) is `field_data["meshio:length_unit"]` in metres.
- Base materials (`<basematerials>`, a triangle's `pid`/`p1` or its object's `pid`/`pindex`) are the `uint8` RGBA cell data `3mf:color`; a triangle with no material is `0 0 0 0`.
- A model that lists `requiredextensions` is refused by name: the materials, production, slice and beam-lattice extensions are out of scope.

## Writing

```python
meshioplusplus.write("plate.3mf", mesh)
```

`[Content_Types].xml`, `_rels/.rels` and `3D/3dmodel.model`, deflated when that helps and dated 1980-01-01 so the bytes depend on the mesh alone. Each cell region is an object with the vertices it uses (a triangle in several regions goes to the first, in region order; the triangles in none form one more object, `mesh`), each object a build item. The unit comes from `meshio:length_unit` when it is one of 3MF's, millimetre otherwise; `3mf:color` becomes base materials. Triangles only: other cells are refused by name.

## Verification

lib3mf 2.5 (the 3MF Consortium's library) reads what meshio++ writes, and meshio++ reads what lib3mf writes back (`tests/python/test_threemf.py`). PrusaSlicer has not been tried.

# Gmsh (`.msh`)

The [Gmsh](https://gmsh.info/doc/texinfo/gmsh.html#File-formats) mesh format, supporting file versions **2.2**, **4.0** and **4.1**, in ASCII and binary.

| | |
|---|---|
| **Format name** | `gmsh` (writes v4.1), `gmsh22` (writes v2.2) |
| **Extensions** | `.msh` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshioplusplus

mesh = meshioplusplus.read("mesh.msh")          # version auto-detected
meshioplusplus.gmsh.write("out.msh", mesh,
    fmt_version="4.1",  # "2.2", "4.0" (write not supported), or "4.1"
    binary=True,
    float_fmt=".16e",
)
```

- **`fmt_version`** — output MSH version (write supports `"2.2"` and `"4.1"`).
- **`binary`** — write the node/element/data bodies in binary (`True`) or ASCII (`False`).
- **`float_fmt`** — ASCII coordinate format string.

Via the generic dispatch, `file_format="gmsh"` writes v4.1 and `file_format="gmsh22"` writes v2.2.

## File structure

`$MeshFormat` is read first: `version filetype datasize` (`filetype` 0=ascii, 1=binary; if binary, a 4-byte integer `1` follows to detect endianness). The version string's major component picks the reader: `"2"`/`"2.2"` → the 2.2 reader, `"4.0"` → the 4.0 reader, `"4"`/`"4.1"` → the 4.1 reader. `$Comments` blocks before `$MeshFormat` are skipped.

**Version 2.2**: `$PhysicalNames`; `$Nodes` (ascii `id x y z` rows, or binary `(int32 id, 3×double)` structs); `$Elements` (ascii: `id type ntags tag1..tagN node1..nodeK` per line; binary: a per-block header `elem_type num_elems num_tags` then flat `int32` rows of `[element_id, tags…, nodes…]`). The first two element tags are `gmsh:physical`/`gmsh:geometrical`; further tags produce a warning. Optional `$Periodic`, `$NodeData`/`$ElementData` (post-processing views).

**Version 4.0** adds `$Entities`: per-dimension counts, then per-entity `tag, 6 doubles (bbox, discarded), num_physicals, physicals[]` (plus, for dim&gt;0, a discarded BREP-bounding-entity list). `$Nodes` and `$Elements` are grouped into per-entity blocks; critically, elements reference nodes **by node tag, not by array index**, requiring an inverse tag→index remap before connectivity can be built.

**Version 4.1** (the default write target) restructures the block headers: `$Entities` starts with `numPoints numCurves numSurfaces numVolumes` (4 `size_t`s); each entity also records its `numBoundingXxx` BREP-boundary count (populating `cell_sets["gmsh:bounding_entities"]`). `$Nodes` header: `numEntityBlocks numNodes minNodeTag maxNodeTag`; per block `entityDim entityTag parametric numNodesInBlock`, then a **node-tag list**, then a matching **coordinate list** — tags may be sparse or out of order (handled via `np.unique(..., return_inverse=True)`). `$Elements` header: `numEntityBlocks numElements minElementTag maxElementTag`; per block `entityDim entityTag elementType numElementsInBlock`, then rows of `elementTag node1..nodeK`.

`$NodeData`/`$ElementData` (used across versions): a block of string tags (the first is the data name), real tags (only "time" is meaningfully used), integer tags `[timestep, num_components, num_items]`, then `num_items × (1+num_components)` values — the leading column is the 1-based node/element index, discarded after use.

## Cell types & node ordering

`$Elements` type codes map to meshio++ types via a large table (`_gmsh_to_meshio_type` in `common.py`) covering linear through very-high-order elements — types 1 through 110, spanning `line`…`line11`, `triangle`…`triangle66`, `quad`…`quad121`, `tetra`…`tetra286`, `hexahedron`…`hexahedron1000`, `wedge`…`wedge550`, and `pyramid`/`pyramid13`/`pyramid14`.

Five element types need a node-order permutation between Gmsh and meshio++ (everything else uses natural order):

| type | gmsh → meshio++ | meshio++ → gmsh |
|---|---|---|
| `tetra10` | `[0,1,2,3,4,5,6,7,9,8]` (self-inverse) | same |
| `hexahedron20` | `[0,1,2,3,4,5,6,7,8,11,13,9,16,18,19,17,10,12,14,15]` | `[0,1,2,3,4,5,6,7,8,11,16,9,17,10,18,19,12,15,13,14]` |
| `hexahedron27` | hex20 permutation + `[22,23,21,24,20,25,26]` | hex20 inverse + `[24,22,20,21,23,25,26]` |
| `wedge15` | `[0,1,2,3,4,5,6,9,7,12,14,13,8,10,11]` | `[0,1,2,3,4,5,6,8,12,7,13,14,9,11,10]` |
| `pyramid13` | `[0,1,2,3,4,5,8,10,6,7,9,11,12]` | `[0,1,2,3,4,5,8,9,6,10,7,11,12]` |

## Data mapping

- `cell_data["gmsh:physical"]`, `cell_data["gmsh:geometrical"]` — the first two element tags (also recognized as `"cell_tags"` on write).
- `point_data["gmsh:dim_tags"]` — v4.1 only, an `(N, 2)` int array of `(entity_dim, entity_tag)` per node.
- `cell_sets["gmsh:bounding_entities"]` — v4.1 only.
- `field_data[name] = [phys_num, phys_dim]` — from `$PhysicalNames`.
- Arbitrary `point_data`/`cell_data` from `$NodeData`/`$ElementData`.
- `mesh.gmsh_periodic` — a mesh-level attribute (not a data-dict key) holding `[dim, (slave_tag, master_tag), affine_or_None, node_pairs]` per periodic relation, from `$Periodic`.

## Quirks & limitations

- Version strings are normalized: `"2"` → 2.2, `"4"` → 4.1.
- Gmsh can't distinguish a `(n,)` shape from `(n,1)` for post-processing data; the reader squeezes single-component arrays to 1D.
- Elements in v4.0/4.1 are addressed by **node tag**, not array position — the most structurally distinctive quirk of this format relative to nearly every other one meshio++ supports.
- v4.1 write requires `gmsh:dim_tags` in `point_data` to emit more than one cell type; without it, only a single cell type can be written (`WriteError` otherwise).
- `$Periodic` record layout differs across all three versions (e.g. v4.0 binary uses a *negative* node count as a sentinel meaning "an affine transform follows", then reads a fixed 16 floats for that transform).
- The C++ type table covers up through `hexahedron125`/`tetra286`(sic — the exact upper bound is a curated subset, not the full ~110-entry Python table); a file referencing a higher-order type outside that subset falls back to Python transparently.
- The C++ shim always tries the C++ reader first, falling back to Python on any exception. On write, C++ is only attempted for `float_fmt == ".16e"`, no `gmsh_periodic`, and (`fmt_version == "2.2"`) or (`"4.1"` with no `gmsh:dim_tags`) — meaning v4.0 write, and any v4.1 write carrying `gmsh:dim_tags` or periodic data, always go through Python.

## Notes

- `tests/python/meshes/msh/insulated-2.2.msh` and `insulated-4.1.msh` — the same mesh (Gmsh 4.2.2, `-format msh2` for the 2.2 variant): 111 triangles (2 Physical Surfaces) + 21 lines (1 Physical Line). Used to check point sums, cell counts, and `gmsh:physical`/`gmsh:geometrical`/`cell_sets` consistency, including the v4.1 `$Entities`-bearing variant.

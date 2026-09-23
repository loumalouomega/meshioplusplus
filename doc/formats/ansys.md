# Ansys Fluent mesh (`.msh`)

The ANSYS Fluent `.msh` mesh format, also written by TGrid, GAMBIT and ANSYS Meshing: parenthesis-nested "Scheme-like" sections, in ASCII, binary, or a mix of both within one file.

| | |
|---|---|
| **Format name** | `ansys` |
| **Extensions** | `.msh` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshioplusplus

mesh = meshioplusplus.read("mesh.msh", file_format="ansys")
meshioplusplus.ansys.write("out.msh", mesh, binary=True)
```

- **`binary`** — write binary (`True`) or ASCII section bodies.

## File structure

Every section is `(<index> ...)`; the file is read as bytes, since ASCII and binary sections may be mixed in one file. The index may be a bare decimal (ASCII payload) or prefixed `20`/`30` (binary payload: `20xx` = float32/int32, `30xx` = float64/int64). Every integer in a header or an ASCII body is **hexadecimal**.

- `(0 ...)` comment, `(1 "...")` header: skipped (brackets inside quotes are ignored). `(2 <dim>)` gives the dimension (2 or 3).
- **Nodes** — `(<pfx>10 (zone-id first last type ND) (...))`. The zones may come in any order (TGrid lists them out of order); each is placed by its id range, and ids are rebased so the smallest `first` becomes point 0.
- **Faces** — `(<pfx>13 (zone-id first last bc-type face-type) (...))`. A row is `n0 .. nk c0 c1`: the face's nodes and the two cells it separates, `0` meaning none. Face types 2, 3 and 4 have that many nodes; in a mixed (`0`) or polygonal (`5`) zone each row leads with its node count, in ASCII and binary alike.
- **Cells** — `(<pfx>12 (zone-id first last zone-type element-type))`. A Fluent cell zone normally only declares a range of cell ids (a mixed zone may add the list of each cell's type); `zone-type 0` is a dead zone.
- `(39 (id type name) ...)` and `(45 (id type name) ...)` name the zones; their id is **decimal**, unlike every other header.
- A section is a declaration when no body follows its header, including a binary one closed by `End of Binary Section NNNN)`. Unknown sections are skipped; an unknown binary one up to its `End of Binary Section` marker.

### Cells from faces

The cells are rebuilt from their faces. The right-hand normal of `n0 .. nk` points into `c0` (in 2-D, walking `n0 → n1` leaves `c0` on the left), so each face is outward for `c1` and reversed for `c0`. Each cell's outward faces, taken in file order, go through the kernel the OpenFOAM reader also uses (`formats/face_cells_common.hpp`, `_face_cells.py`): four triangles make a `tetra`, a quad and four triangles a `pyramid`, two triangles and three quads a `wedge`, six quads a `hexahedron`, all in meshio++'s node order with a positive volume (the orientation is checked geometrically). Anything else is a `polyhedron<N>` of its outward faces. In 2-D a cell's edges are chained into one counter-clockwise ring: a `triangle`, a `quad` or a `polygon<N>`. A declared cell that no face names is skipped with a warning.

Boundary faces are kept as surface cells (`triangle`, `quad`, `polygon<N>`, or `line` in 2-D), wound outward from the domain. Faces of an `interior` zone (bc-type 2) with a cell on both sides are dropped; an interior-zone face with one side empty is kept (FEconv writes every face into one interior zone).

### Legacy layout

A cell section whose element type has a fixed node count **and** a body is meshio's own layout: the body is the connectivity, `npc` hexadecimal ids per cell. That is what the writer produces, and such a file reads as it always did, cells only, with no zone data or regions. Fluent itself never writes it and cannot read it.

Write emits, in order:

```
(1 "meshio++ VERSION")
(2 DIM)
(10 (0 1 N_POINTS 0))                          -- node-count declaration
(12 (0 1 TOTAL_CELLS 0))                       -- cell-count declaration
(10|3010 (1 1 N 1 DIM)( DATA ))                -- node block
(12|2012|3012 (1 FIRST LAST 1 ANSYS_TYPE)( DATA ))   -- per cell block
```

## Cell types

Element-type codes of a `12` section:

| code | meshio++ type | nodes |
|---|---|---|
| 0 | mixed | — |
| 1 | `triangle` | 3 |
| 2 | `tetra` | 4 |
| 3 | `quad` | 4 |
| 4 | `hexahedron` | 8 |
| 5 | `pyramid` | 5 |
| 6 | `wedge` | 6 |
| 7 | polyhedron | — |

On read the type comes from the faces, not from this code. meshio++ → Ansys codes on write: `triangle:1, tetra:2, quad:3, hexahedron:4, pyramid:5, wedge:6` (no writer support for mixed or polyhedral zones).

## Data mapping

- `cell_data["ansys:zone"]` — the zone id of every cell, volume and boundary alike (face-based files only).
- Regions — one cell region per zone, named from its `(39 ...)`/`(45 ...)` declaration (`zone_<id>` when there is none), with the zone id as `tag` and the cells' dimension as `dim`. See [named regions](../regions.md).
- No point data or field data.

## Quirks & limitations

- Before v16.4.0 the reader returned the faces only, never a volume cell, and the C++ reader handed every file with face sections to Python. It now rebuilds the cells in both engines, which give identical meshes, zones and regions. **Interior faces are no longer returned** as `triangle`/`quad` blocks.
- Out-of-order node zones, a missing blank in `(13(`, GAMBIT's boundary rows with `c0 = 0`, bodyless binary declarations and binary mixed face zones are all accepted; before v16.4.0 each of them failed.
- A leading node count is assumed for polygonal rows in a mixed face zone (the samples only have fixed-size rows there).
- Hanging-node trees (`(58 ...)`/`(59 ...)`), periodic shadows and cell-tree data are skipped.
- The writer still emits the legacy layout; a face-based writer that Fluent reads is on the [roadmap](../roadmap.md).
- 2D/3D validity (`dim in {2,3}`) is only checked on write.

## Notes

- `tests/python/meshes/ansys/cells3d.msh`, `tgrid2d.msh`, `gambit2d.msh` and `binary3d.msh`, generated by `tools/gen_feconv_quirk_fixtures.py`, hold cells known only through their faces: two tetrahedra sharing a face, a hexahedron, a wedge, a pyramid and a polyhedron; 2-D TGrid and GAMBIT variants; binary sections. Tests check positive volumes, outward boundary and polyhedron faces, zones, regions and C++/Python parity. The 29 Fluent meshes of FEconv's `examples/` (GPL, not committed) read identically in both engines, all cells positively oriented; `ansys_mesh.msh` gives the same 224 tetrahedra as its I-DEAS UNV twin (`tests/python/test_feconv_examples.py`, with `MESHIOPLUSPLUS_FECONV_DIR` set).
- The legacy write-read round trip is tested for synthetic meshes (`empty_mesh`, `tri_mesh`, `tri_mesh_2d`, `quad_mesh`, `tri_quad_mesh`, `tet_mesh`, `hex_mesh`, `pyramid_mesh`, `wedge_mesh`), ASCII and binary.
- `.msh` is shared with [`gmsh`](./gmsh.md) and [`freefem`](./freefem.md); on auto-detection `ansys` is tried first. Pass `file_format` to disambiguate.

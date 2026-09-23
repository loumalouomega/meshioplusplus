# COMSOL mesh (`.mphtxt`, `.mphbin`)

[COMSOL Multiphysics](https://www.comsol.com)' native mesh files: text `.mphtxt` and its binary twin `.mphbin`, the documented way into and out of COMSOL that keeps domain ids and selections.

| | |
|---|---|
| **Format names** | `mphtxt`, `mphbin` ([its own page](./mphbin.md)) |
| **Extensions** | `.mphtxt`, `.mphbin` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshioplusplus

mesh = meshioplusplus.read("model.mphtxt")      # or "model.mphbin"
mesh.cell_data["mphtxt:geom"]                   # the geometric entity of every element
[(r.name, r.dim) for r in mesh.regions]         # COMSOL Selections

meshioplusplus.mphtxt.write("out.mphtxt", mesh)
meshioplusplus.mphbin.write("out.mphbin", mesh)
```

Neither writer takes keyword arguments. Both engines read the same mesh and write the same bytes.

## File structure

Both files serialise the same sequence of integers, doubles and strings ([file structure](https://doc.comsol.com/6.1/doc/com.comsol.help.comsol/comsol_api_fileformats.53.18.html)):

1. The file version `0 1`.
2. The tags, one per object (`mesh1`, `mesh1_sel1`, …), and the type names (`obj`), one per object.
3. The objects, each opened by `0 0 1` and a class name. meshio++ reads `Mesh` and `Selection`; the first object of any other class (a `Geom3` geometry, for instance) stops the read with a warning, since its length cannot be known. A file without a Mesh object is refused.

A **`Mesh` object** holds its version, the space dimension, the vertex count, the lowest vertex index (0 in current files, 1 in some old ones) and the coordinates, then per element type: its name, nodes per element, the element count and connectivity, and one **geometric entity index** per element. Version 4 (what COMSOL writes today) stops there. Older versions add, per type, the elements' parameter values and up/down pairs. A parameter record holds a varying number of values per node — one on an edge, three on a 3-D boundary — so the reader finds the record length by where the rest of the element record fits. A file may give no entity indices at all; its elements then get COMSOL's default, domain 1 or entity 0.

A **`Selection` object** holds its version, a label, the tag of the Mesh it selects from, a dimension and a list of entities.

**Text.** Values are whitespace-separated and `#` starts a comment to the end of the line. A string is its length, one blank and exactly that many characters, which may include blanks and `#`: `12 Copper Piece`.

**Binary.** Integers are little-endian int32, doubles float64, and a string is an int32 length followed by one int32 code point per character. There is no comment, so no provenance tag either.

Several Mesh objects are merged into one mesh, their points one after the other; each object's cells become a `cell` region named by its tag.

## Cell types & node order

| COMSOL | meshio++ | COMSOL | meshio++ |
|---|---|---|---|
| `vtx` | `vertex` | | |
| `edg` | `line` | `edg2` | `line3` |
| `tri` | `triangle` | `tri2` | `triangle6` |
| `quad` | `quad` | `quad2` | `quad9` |
| `tet` | `tetra` | `tet2` | `tetra10` |
| `pyr` | `pyramid` | `pyr2` | `pyramid14` |
| `prism` | `wedge` | `prism2` | `wedge18` |
| `hex` | `hexahedron` | `hex2` | `hexahedron27` |

COMSOL numbers an element's corners in tensor order (x fastest: a quad is `0 1 3 2` as a ring), then every other node of the element's quadratic lattice in lexicographic (z, y, x) order ([Mesh Element Numbering Conventions](https://doc.comsol.com/5.5/doc/com.comsol.help.comsol/comsol_api_mesh.40.33.html)). A `tri2` is corners, then mid-edges 0–1, 0–2, 1–2; a `hex2` puts its bottom-face centre ninth in lattice order. The permutations live in the [node-ordering registry](../node_ordering.md) under the key `mphtxt`. They agree with the tables of [AWS Palace](https://github.com/awslabs/palace)'s COMSOL reader composed with the gmsh ones. On every quadratic element of the published COMSOL files checked (below), each mid-edge node sits on its edge's midpoint. Linear prisms keep their order: the base normal points to the top, as in meshio++.

Before v16.1.0 meshio++ used the identity for every quadratic type, which misplaced the mid-edge nodes of `tri2`, `tet2`, `quad2`, `prism2` and `hex2`, and did not know `pyr2`. Files written by those releases with quadratic elements read back differently now.

## Selections

A Selection becomes a `cell` region named by its label, with the Selection's dimension. Its entries are the elements of that dimension, in the Mesh it names, whose entity index it lists. Duplicate labels get a ` (2)` suffix and a warning.

On write:

- **Entity indices.** They come from `cell_data["mphtxt:geom"]` when the mesh has it. Otherwise the writer numbers them per dimension: each pairwise-disjoint cell region whose cells share a dimension becomes the next entity, and the cells in none share one more. Domains (cells of the space dimension) count from 1, boundaries, edges and points from 0. Regions are taken in `(kind, name, dim, tag)` order.
- **Selections.** Every cell region that is exactly a union of whole entities of one dimension is written as a Selection; any other region is dropped with a warning. Without `mphtxt:geom`, every disjoint region is its own entity, so it always survives, and so does a region that unions several.
- **What is dropped.** Point and side regions are dropped with a warning, as are data arrays other than `mphtxt:geom`.

The writer writes a single Mesh object, so the regions of a multi-object file (which are not unions of whole entities there) are dropped on write.

## Data mapping

- `cell_data["mphtxt:geom"]` — each element's geometric entity index, as in the file: domains from 1; points, edges and boundaries from 0.
- `mesh.regions` — Selections (above), plus one region per Mesh object when there are several.

## Quirks & limitations

- Element parameter values and up/down pairs (Mesh versions before 4) are read past and never written.
- A Selection of a Mesh object the file does not contain is skipped with a warning.
- COMSOL's geometry objects, sectionwise result exports and spreadsheet data are not read.
- `mphbin` writes int32 values; a mesh with more than 2³¹−1 vertices or elements cannot be written.

## Notes

- The readers were run over the COMSOL files published with [deal.II](https://github.com/dealii/dealii/tree/master/tests/grid/grids/comsol), [FEconv](https://github.com/victorsndvg/FEconv/tree/master/examples/test/mphtxt) and [Wolfram's FEMAddOns](https://github.com/WolframResearch/FEMAddOns/tree/master/Resources/Meshes/Comsol) — 40 files, versions 2 and 4, linear and quadratic, with selections and several Mesh objects. Both engines read every mesh identically. One file holds only geometry and is refused, as it should be. These files are not redistributed here.
- The fixtures in `tests/python/meshes/comsol/` are generated by `tools/gen_comsol_fixtures.py`, which spells out COMSOL's numbering from the lattice rule rather than from the registry: two domains with Selections (also as `.mphbin`), every quadratic type, a version 2 Mesh with parameters and up/down pairs, and two Mesh objects.

# libMesh (`.xda`, `.xdr`)

libMesh's native mesh file (`XdrIO`), also the mesh MOOSE writes with `--mesh-only` and in its checkpoints. `.xda` is ASCII, `.xdr` the same stream of values in big-endian XDR. MOOSE users are also served by [Exodus](./exodus.md); this format is for meshes that exist only in libMesh's own file.

| | |
|---|---|
| **Format name** | `libmesh` |
| **Extensions** | `.xda`, `.xdr` (also recognised by content: the `libMesh-` version string, as text or as an XDR string) |
| **Read / Write** | ✓ / — |
| **Extra dependencies** | — |

## Reading

```python
import meshioplusplus

mesh = meshioplusplus.read("mesh.xdr")   # or meshioplusplus.libmesh.read(...)
```

`read` takes no options. The encoding comes from the content, not the extension. Both engines (the C++ core and the pure-Python reference) read the same meshes.

## The stream

A version string (`libMesh-0.7.0+` to `libMesh-1.8.0`), the element and node counts, four "inline or separate file" flags (boundary conditions, subdomain, processor and p-level ids), then:

| Section | Since | Read as |
|---|---|---|
| per-field integer sizes | 0.9.2 | the connectivity integer width (4 or 8 bytes in XDR; header integers are 8 bytes from 1.3.0 on) |
| extra integer names, elemsets | 1.8.0 | skipped |
| subdomain id → name map | 0.9.2 | region names |
| one connectivity block per refinement level | | cells (see below) |
| coordinates | | points |
| node unique ids | 0.9.6 | skipped |
| side sets (with a name map from 0.9.2) | | `side` regions |
| node sets | 0.9.2 | `point` regions |
| edge and shell-face sets | 1.1.0 | skipped with a warning |

Legacy pre-`libMesh` files (`DEAL 003`, `LIBM 0`) are refused, as libMesh itself refuses them. Compressed files (`.xda.gz`, `.xdr.bz2`) must be decompressed first.

## Cells

- Only **active** elements (the leaves of the refinement tree) become cells. When the file holds refinement levels, `libmesh:level` records each cell's level.
- The subdomain id is the `libmesh:subdomain` cell data and a `cell` region per subdomain, named by the file's map, else `subdomain_<id>` (tag = id). An inline p-level becomes `libmesh:p_level`.
- Node ids that no element uses (libMesh writes their coordinates as NaN) are dropped; the original ids are then kept as `libmesh:id` point data.

| libMesh type | meshio++ type |
|---|---|
| EDGE2 / EDGE3 / EDGE4 | `line` / `line3` / `line4` |
| TRI3, TRISHELL3, TRI3SUBDIVISION / TRI6 / TRI7 | `triangle` / `triangle6` / `triangle7` |
| QUAD4, QUADSHELL4 / QUAD8, QUADSHELL8 / QUAD9, QUADSHELL9 | `quad` / `quad8` / `quad9` |
| TET4 / TET10 / TET14 | `tetra` / `tetra10` / `tetra10` (face nodes dropped) |
| HEX8 / HEX20 / HEX27 | `hexahedron` / `hexahedron20` / `hexahedron27` |
| PRISM6 / PRISM15 / PRISM18 / PRISM20, PRISM21 | `wedge` / `wedge15` / `wedge18` / `wedge18` (extra nodes dropped) |
| PYRAMID5 / PYRAMID13 / PYRAMID14 / PYRAMID18 | `pyramid` / `pyramid13` / `pyramid14` / `pyramid14` (extra nodes dropped) |
| NODEELEM | `vertex` |
| infinite elements (INF*) | skipped with a warning |
| C0POLYGON, C0POLYHEDRON, REMOTEELEM | `ReadError` |

Dropping extra nodes is logged and recorded as a provenance note; the dropped nodes stay as points.

## Node order

HEX20, HEX27, PRISM15 and PRISM18 list their vertical mid-edge nodes before the top ring (and HEX27 its face centres bottom, y−, x+, y+, x−, top), so they use the `"libmesh"` tables of the [node-ordering registry](../node_ordering.md), taken from libMesh's own VTK connectivity (`cell_hex27.C` …). Every other type is already in meshio++'s order.

## Boundaries

A side set entry names an element and one of its sides in libMesh's numbering (`Hex8::side_nodes_map` …). It becomes a `side` region entry through the facet with the same corner nodes. libMesh keeps boundary conditions on the coarse elements, so a side of a refined element is carried down to every active descendant whose side lies on it. The sides of line elements (their end points) have no side-region form and are skipped with a warning.

## Validation

The fixtures under `tests/python/meshes/libmesh/` are written by `tools/gen_libmesh_fixtures.py` in libMesh's own layout, since libMesh (LGPL) is not a dependency. Outside the repository the reader was run on libMesh's `reference_elements/` and `tests/meshes/` samples (both encodings for every element type, and an AMR mesh with side sets): every cell has a positive volume and every higher-order node sits where meshio++'s tables put it. Three of libMesh's hand-written reference files (`one_pyramid13/14/18.xda`) announce inline subdomain ids they do not contain; libMesh itself cannot read them either, and meshio++ refuses them.

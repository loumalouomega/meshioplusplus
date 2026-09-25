# libMesh (`.xda`, `.xdr`)

libMesh's native mesh file (`XdrIO`), also the mesh MOOSE writes with `--mesh-only` and in its checkpoints. `.xda` is ASCII, `.xdr` the same stream of values in big-endian XDR. MOOSE users are also served by [Exodus](./exodus.md); this format is for meshes that exist only in libMesh's own file.

| | |
|---|---|
| **Format name** | `libmesh` |
| **Extensions** | `.xda`, `.xdr`, and each with `.gz` or `.bz2` (also recognised by content: the `libMesh-` version string, as text or as an XDR string) |
| **Read / Write** | ✓ / ✓ (writer v16.11.0) |
| **Extra dependencies** | — (the native core uses zlib for gzip and, when built with `MESHIOPLUSPLUS_WITH_BZIP2`, libbz2 for bzip2; without them Python does it) |

## Reading

```python
import meshioplusplus

mesh = meshioplusplus.read("mesh.xdr")   # or meshioplusplus.libmesh.read(...)
```

`read` takes no options. The encoding comes from the content, not the extension, and gzip or bzip2 compression is inflated whatever the file is called (libMesh itself compresses only ASCII files: its `.xdr.gz` and `.xdr.bz2` are plain XDR). Both engines (the C++ core and the pure-Python reference) read the same meshes; the core inflates gzip through zlib and, since v16.12.0, bzip2 through libbz2 when built with `-DMESHIOPLUSPLUS_WITH_BZIP2=ON` (`_core.__has_bzip2__`); a build without it leaves bzip2 to the Python reader, so its native CLI and C, Fortran, Julia, R and WASM bindings read gzip only.

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
| edge and shell-face sets | 1.1.0 | `line` cells and `cell` regions (see [Boundaries](#boundaries)) |

Legacy pre-`libMesh` files (`DEAL 003`, `LIBM 0`) are refused, as libMesh itself refuses them.

## Cells

- Only **active** elements (the leaves of the refinement tree) become cells. When the file holds refinement levels, `libmesh:level` records each cell's level.
- **The refinement tree** (v16.12.0) is kept too, as field data, so the writer can write it back: `libmesh:tree`, one row per element of the file in file order (level by level), `cell, parent, type, subdomain, p_level` (the cell is −1 for a refined element, the parent −1 at level 0, the type libMesh's `ElemType` number), and `libmesh:tree:nodes`, the element's nodes as point indices in libMesh's order, −1 past them.
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
| C0POLYGON, C0POLYHEDRON, REMOTEELEM | `ReadError` (see below) |

Dropping extra nodes is logged and recorded as a provenance note; the dropped nodes stay as points.

No `.xda`/`.xdr` file can hold a C0POLYGON or C0POLYHEDRON: `XdrIO` writes no node count per element and takes it from the element type, which for these two types has none, so libMesh cannot write them either (`XdrIO::pack_element` asserts the fixed count). They remain a `ReadError`.

## Node order

HEX20, HEX27, PRISM15 and PRISM18 list their vertical mid-edge nodes before the top ring (and HEX27 its face centres bottom, y−, x+, y+, x−, top), so they use the `"libmesh"` tables of the [node-ordering registry](../node_ordering.md), taken from libMesh's own VTK connectivity (`cell_hex27.C` …). Every other type is already in meshio++'s order.

## Boundaries

A side set entry names an element and one of its sides in libMesh's numbering (`Hex8::side_nodes_map` …). It becomes a `side` region entry through the facet with the same corner nodes. libMesh keeps boundary conditions on the coarse elements, so a side of a refined element is carried down to every active descendant whose side lies on it. The sides of line elements (their end points) have no side-region form and are skipped with a warning.

Edge sets and shell-face sets (libMesh 1.1.0 on) share the side sets' id space and name map:

- An **edge set** entry names an element and one of its edges (`Hex8::edge_nodes_map` …). Each named edge becomes a `line` cell (a `line3` with the element's mid-edge node when the element is quadratic), shared by every entry naming it, in a `cell` region `<name>:edge` (tag = id). These cells are not libMesh elements: their `libmesh:subdomain` (and `libmesh:level`, `libmesh:p_level`) is −1. An edge of a refined element spans its coarse edge.
- A **shell-face set** entry names a 2-D element and its face 0 or 1. It becomes a `cell` region `<name>:shellface0` or `<name>:shellface1` on the element (on every active descendant of a refined one). libMesh gives the two faces no geometric meaning, so the names keep its numbers.

`<name>` is the side-set name map's entry for the id, else `boundary_<id>`.

## Writing

```python
meshioplusplus.write("mesh.xdr", mesh)   # XDR; "mesh.xda" for ASCII
meshioplusplus.write("mesh.xda.gz", mesh)  # gzip-compressed ASCII
```

`write` takes no options and emits the `libMesh-1.8.0` layout that `XdrIO::write` produces: 8-byte ids, inline subdomain ids, no unique ids, no processor ids. The encoding comes from the extension (`.xdr` is XDR, anything else ASCII). A trailing `.gz` or `.bz2` compresses an ASCII file, as libMesh does; a `.xdr.gz` or `.xdr.bz2` is written as plain XDR, because that is what libMesh writes and all it reads under those names (its `Xdr` class opens XDR files with plain stdio). The native writer compresses (v16.12.0: gzip through zlib, bzip2 through libbz2 when built with it); a build without them has the core write the plain stream and Python compress it. gzip has no timestamp (the header of Python's `gzip.compress(data, mtime=0)`), so the bytes are reproducible, and both engines write the same bytes.

- **Cells.** Every cell whose type is in the table above (the plain type, not the shell or subdivision variant) becomes an element.
- **The refinement tree.** When `libmesh:tree` still describes the cells (each written cell is exactly one of its leaves, with the same type and nodes; every row valid, levels in order), the whole tree is written back, level by level, the leaves' subdomains and p-levels from the cells and their ancestors' from the tree (v16.12.0). Otherwise, or without it, each cell is a level-0 element and the mesh is written flat, with a warning when a tree was there but no longer matches. `polygon`, `polyhedron`, Lagrange and other cells libMesh has no type for are dropped with a warning and a provenance note. The `"libmesh"` node-ordering tables are applied in reverse.
- **Subdomains.** `libmesh:subdomain` when present, else the first `cell` region that holds the cell (its tag, or a fresh id past the largest tag), else 0. Ids must fit libMesh's `subdomain_id_type` (0 to 65534). Region names other than `subdomain_<id>` go into the subdomain name map.
- **Boundaries.** `side` regions become side sets, each entry matched to the libMesh side with the same corner nodes. When the tree is written, libMesh keeps boundary ids on level-0 elements only, so each set is lifted: a level-0 element's side (edge, shell face) is in the set when every active piece of it is; entries that cover only part of one are dropped with a warning. An edge cell spanning a whole level-0 edge, as the reader makes from a refined element's edge set, goes on that edge. The `<name>:edge` and `<name>:shellface<k>` cell regions the reader makes become edge and shell-face sets again; their line cells are not written as elements. `point` regions become node sets. A region's tag is its id; untagged regions get fresh ids (one id space for side, edge and shell-face sets, another for node sets). Names other than `boundary_<id>`/`nodeset_<id>` go into the name maps. A side or edge that no written element has is dropped with a warning.
- **Node ids.** Points are numbered in order, unless `libmesh:id` holds distinct non-negative ids (as the reader leaves it when the file had unused ids): then those ids are kept and the unused ones written as NaN in XDR, as libMesh writes them, and as 0 in ASCII, because libMesh's ASCII reader cannot parse the `nan` its own writer prints (it only loads the nodes elements use, so the value is never looked at).
- **p-levels.** `libmesh:p_level` is written inline. `libmesh:level` and other data arrays are not written.

## Validation

The fixtures under `tests/python/meshes/libmesh/` are written by `tools/gen_libmesh_fixtures.py` in libMesh's own layout, since libMesh (LGPL) is not a dependency, except `tree_quad.xd?` and `tree_hex20.xdr`, which libMesh 1.8.0 wrote (from `tools/libmesh_tree_driver.cpp`): a quad mesh refined into three levels with a side set, a node set and p-levels, and a refined HEX20 with an edge set on a level-0 edge. meshio++ writes them back with the same tree and sets, and ASCII line for line as libMesh did (its node-set rows aside, which libMesh orders its own way). For v16.12.0, libMesh read the files meshio++ writes back from these and from libMesh-refined hexahedra with the same elements at every level, active elements, sets and measure, and read its `.xda.gz`, `.xda.bz2`, `.xdr.gz` and `.xdr.bz2`. For v16.11.0 libMesh 1.8.0 was built from source outside the repository: it reads every fixture and every file the writer produces (from these fixtures, from libMesh-written meshes with edge and shell-face sets, gzip and bzip2 files, and from Gmsh, Exodus and Abaqus meshes) with the same elements, subdomains, side, edge, shell-face and node sets and names, and the same measure; meshio++ reads libMesh's own files the same way. Outside the repository the reader was run on libMesh's `reference_elements/` and `tests/meshes/` samples (both encodings for every element type, and an AMR mesh with side sets): every cell has a positive volume and every higher-order node sits where meshio++'s tables put it. Three of libMesh's hand-written reference files (`one_pyramid13/14/18.xda`) announce inline subdomain ids they do not contain; libMesh itself cannot read them either, and meshio++ refuses them.

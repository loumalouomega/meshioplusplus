<!--pytest-codeblocks:skipfile-->
# libMesh `.xda`/`.xdr` fixtures

Written by `tools/gen_libmesh_fixtures.py`, value by value in the layout of libMesh's `XdrIO::write` and `xdr_cxx.C`, each once as ASCII (`.xda`) and once as XDR (`.xdr`). libMesh itself (LGPL-2.1) is not a dependency and none of its sample meshes is copied.

| File | Version | What it exercises |
|---|---|---|
| `hex27.*` | `libMesh-1.3.0`, 8-byte fields | two HEX27 in subdomains 1 (`left`) and 2, unique ids, side sets `inlet` (1) and 2 on the end faces, node set `wall` (5) |
| `one_hex.*` | `libMesh-0.7.0+` | one HEX8 and six side sets in the old 4-byte header layout |
| `amr_quad.*` | `libMesh-1.8.0` | a QUAD4 refined once: only the four children are cells; inline p-levels; an extra element integer; an unused node id whose coordinates are NaN; a side set on the parent's bottom edge that reaches the two bottom children |
| `tet14_prism18.*` | `libMesh-1.3.0`, 4-byte fields | TET14 and PRISM18, read as `tetra10` and `wedge18` |
| `edges_shell.*` | `libMesh-1.8.0` | a HEX20 (subdomain 1) with a QUAD8 shell (`skin`, 2) on its top face; edge sets `axis` (11) and 12, the latter naming one edge twice (from the hex and from the shell), so it reads as one `line3`; shell faces `front` (20, face 0) and 21 (face 1); a side set `bottom` (1) |
| `hex27.xda.gz`, `hex27.xdr.bz2` | as `hex27.*` | the same meshes gzip- and bzip2-compressed (gzip without a timestamp) |
| `tree_quad.xda`, `tree_quad.xdr` | `libMesh-1.8.0`, **written by libMesh** | a 2×2 QUAD4 square refined uniformly, then one child again (24 elements over three levels, 19 active), side sets from `build_square` (0 named `bottom`), a node set `corner` (7), p-level 1 on every third element, subdomain `plate` |
| `tree_hex20.xdr` | `libMesh-1.8.0`, **written by libMesh** | a HEX20 refined uniformly (9 elements), with an edge set `axis` (11) on the level-0 element's edge 0 |

libMesh 1.8.0, built from source outside the repository, reads every fixture here except the bzip2 one (that build had no bzip2) with the sets and names listed; the ASCII `amr_quad.xda` loses its boundary data there, because libMesh's ASCII reader cannot parse the `nan` coordinates libMesh writes for unused ids (the XDR twin reads fully). Outside the repository the reader was also run on libMesh's own sample meshes (`reference_elements/`, `tests/meshes/xdrio_elements/`, `mesh_assign_test_mesh.xda` and the `.xda.gz` meshes): every cell has a positive volume and every mid-edge and face node sits where meshio++'s tables put it.

The `tree_*` files were written by libMesh 1.8.0 itself, running `tools/libmesh_tree_driver.cpp` (libMesh's LGPL covers the program, not the meshes it writes); they are libMesh's own refinement trees, so the writer is checked against them.

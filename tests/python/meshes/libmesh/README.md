<!--pytest-codeblocks:skipfile-->
# libMesh `.xda`/`.xdr` fixtures

Written by `tools/gen_libmesh_fixtures.py`, value by value in the layout of libMesh's `XdrIO::write` and `xdr_cxx.C`, each once as ASCII (`.xda`) and once as XDR (`.xdr`). libMesh itself (LGPL-2.1) is not a dependency and none of its sample meshes is copied.

| File | Version | What it exercises |
|---|---|---|
| `hex27.*` | `libMesh-1.3.0`, 8-byte fields | two HEX27 in subdomains 1 (`left`) and 2, unique ids, side sets `inlet` (1) and 2 on the end faces, node set `wall` (5) |
| `one_hex.*` | `libMesh-0.7.0+` | one HEX8 and six side sets in the old 4-byte header layout |
| `amr_quad.*` | `libMesh-1.8.0` | a QUAD4 refined once: only the four children are cells; inline p-levels; an extra element integer; an unused node id whose coordinates are NaN; a side set on the parent's bottom edge that reaches the two bottom children |
| `tet14_prism18.*` | `libMesh-1.3.0`, 4-byte fields | TET14 and PRISM18, read as `tetra10` and `wedge18` |

Outside the repository the reader was also run on libMesh's own sample meshes (`reference_elements/`, `tests/meshes/xdrio_elements/`, `mesh_assign_test_mesh.xda` and the `.xda.gz` meshes): every cell has a positive volume and every mid-edge and face node sits where meshio++'s tables put it.

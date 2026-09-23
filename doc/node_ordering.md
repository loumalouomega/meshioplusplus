# Node ordering

meshio++ stores every cell in the [VTK node order](./cell_types.md). A format that numbers the nodes of an element differently needs a permutation on the way in and its inverse on the way out. Since v16.0.0 those permutations live in one registry, keyed by *(format, cell type)*:

- C++: `detail/node_order.hpp`, `node_order(format, cell_type)`;
- Python: `meshioplusplus._node_order`, `node_order(fmt, cell_type)`, `to_meshio(...)` and `from_meshio(...)`.

A type with no entry uses the identity. Every entry holds **both directions as gather tables**:

- `to_meshio[k]` is the file slot meshio++ node `k` comes from on read: `meshio[k] = file[to_meshio[k]]`;
- `from_meshio[j]` is the meshio++ node written to file slot `j`: `file[j] = meshio[from_meshio[j]]`.

Most tables are their own inverse, but not all. MED's `hexahedron27` is not, and a writer that reused its read table would scramble the face centres, so a format never has to know.

## What is in the registry

| Format | Cell types with a permutation | Pinned against |
|---|---|---|
| `med` | `tetra`, `pyramid`, `wedge`, `hexahedron` and their quadratic forms up to `wedge18`/`hexahedron27` | MEDCoupling's `CellModel.cxx` edge and face tables (orientation and mid-edges); `hexahedron27`/`wedge18` (v16.0.0) against a file written by MED-fichier itself and Code_Aster's MED reader |
| `code_aster` | `wedge15`, `wedge18`, `hexahedron20`, `hexahedron27` | Code_Aster's gmsh reader (`inigms.F90`) exactly, and its MED reader (`lrmtyp.F90`) up to a symmetry of the reference cell; see [Code_Aster](./formats/code_aster.md#node-order) |
| `elmer` | `hexahedron20`, `hexahedron27` | ElmerSolver's `elements.def` reference coordinates and its own VTU writer's `Elmer2VtkIndexes` permutation; gmsh meshes converted by ElmerGrid read back node for node as meshio++'s gmsh reader reads them (v16.2.0); see [Elmer](./formats/elmer.md#node-order) |
| `febio` (`.feb` and `.xplt`) | `hexahedron27` | the shape functions of `FEHex27` in FEBio's `FECore/FESolidElementShape.cpp` (the mid-height face centres run y−, x+, y+, x−); FEBio's own parsers and plot files, read against FEBio 4.12 (v16.2.0); see [FEBio](./formats/febio.md#node-order) |
| `flux` | `tetra`, `tetra10`, `pyramid`, `wedge`, `wedge15`, `hexahedron`, `hexahedron20` | FEconv's FLUX samples: every solid is VTK's element mirrored (base face clockwise); read through these tables all have positive Jacobians and mid-edge nodes at edge midpoints, and equal their I-DEAS UNV twins row for row. `wedge15` has no sample and follows the same rule |
| `frd` | `hexahedron20`, `wedge15`, `line3` | `ccx` 2.23 output for the same `.inp` |
| `patran` | `hexahedron20`, `wedge15` | the Patran Reference Manual's Element Library (the vertical mid-edges come before the top ring), with fixtures written from its edge lists (v16.5.0); see [Patran](./formats/patran.md#node-order) |
| `mphtxt` (also `mphbin`) | `quad`, `pyramid`, `hexahedron`, `triangle6`, `quad9`, `tetra10`, `pyramid14`, `wedge18`, `hexahedron27` | COMSOL's "Mesh Element Numbering Conventions" (corners in tensor order, then the quadratic lattice in lexicographic order), real COMSOL files (deal.II, FEconv, Wolfram FEMAddOns), and AWS Palace's COMSOL-to-gmsh tables composed with the gmsh ones (v16.1.0) |
| `unv` | `line3`, `triangle6`, `quad8`, `quad9`, `tetra10`, `pyramid13`, `wedge15`, `hexahedron20` | gmsh's `.unv`/`.msh` twins and Salome's SMESH driver |

The gmsh, CGNS, GiD, Exodus and Kratos tables still live in their own readers. They move here when those formats are next touched ([roadmap §1.15](./roadmap.md)).

## Self-test

`tests/cpp/test_node_order.cpp` and `tests/python/test_node_order.py` check every entry the same way:

- it is a permutation of its cell type's node count, and its two directions are inverses;
- a reference element (mid-edge nodes on midpoints, face and body centres on centroids) written through `from_meshio` and read back through `to_meshio` is again a valid, positively oriented element;
- the C++ and Python tables are identical;
- the `code_aster` tables equal Code_Aster's gmsh reader composed with the gmsh tables, and agree with its MED reader composed with the `med` tables.

<!--pytest-codeblocks:skipfile-->
# DOLFINx VTX (`.bp`) fixtures

Written by `tools/gen_vtx_fixtures.py` with DOLFINx 0.11.0, ADIOS2 2.12.1 and MPICH from conda-forge (September 2026):

| Directory | What it is |
|---|---|
| `heat.bp` | the heat equation, P1 on a 4 x 4 triangle mesh, backward Euler, three written steps (t = 0, 0.01, 0.02); `VTXMeshPolicy.reuse` (the mesh is in step 0 only); BP5 |
| `elastic.bp` | a P2 vector field (`displacement`) and a DG0 field (`density`) on a 2 x 2 x 1 tetrahedral mesh, two steps; `VTXMeshPolicy.update`; BP4 |
| `quads.bp` | a mesh-only file (no function) of a 2 x 3 quadrilateral mesh; BP5 |
| `heat_np2.bp` | `heat.bp`'s run on two MPI ranks: two blocks per step, each with its own point numbering and ghost points |
| `plain.bp` | an ADIOS2 file with no `vtk.xml` schema, which the reader refuses |
| `complex.bp` | two triangles with a complex P1 field as `p_real`/`p_imag`, written with the `adios2` package in DOLFINx's layout (a complex DOLFINx build is a separate PETSc flavour) |

`<name>.reference.npz` is ParaView 6.1's `ADIOS2VTXReader` output after `MergeBlocks`, recorded by `tools/gen_vtx_reference.py`: per step `k`, `k<k>/points`, `k<k>/connectivity`, `k<k>/offsets`, `k<k>/types` and one `k<k>/point/<name>` or `k<k>/cell/<name>` per array, plus `times` and `steps`. ParaView 6.1 aborts on every step after the first of a `reuse` file (`heat`, `heat_np2`), so those steps are not in the reference; `<name>.raw.npz` holds every step's function arrays as the `adios2` package reads them (the rank blocks concatenated, `k<k>/<name>`), their oracle.

DOLFINx 0.11 writes `vtkOriginalPointIds` and `vtkGhostType` a second time in the first step of a `reuse` file, as one-value blocks whose payload is not what their count says; ParaView appends them to the arrays (26 ids for 25 points in `heat`), meshio++ skips them.

The files are ours (MIT); DOLFINx's licence covers the program, not the files it writes.

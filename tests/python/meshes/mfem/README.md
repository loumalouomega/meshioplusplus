<!--pytest-codeblocks:skipfile-->
# MFEM reference fixtures

**The meshes were written by MFEM, not by meshio++.** They come unmodified from the `data/` directory of [MFEM](https://github.com/mfem/mfem) (BSD 3-Clause, notice below; commit `45799ebe`). The grid functions and `reference.npz` were written by MFEM itself through PyMFEM 4.10, by `tools/gen_mfem_fixtures.py`.

| File | What it exercises |
|---|---|
| `star-q2.mesh` | order-2 `H1_2D_P2` nodes on quadrilaterals (`quad9`), boundary edges (`line3`) |
| `escher-p2.mesh` | order-2 tetrahedra (`tetra10`) and their curved triangle faces |
| `fichera-q2.mesh` | order-2 hexahedra (`hexahedron27`): edge, face and interior dofs |
| `fichera-mixed-p2.mesh` | order-2 tetrahedra, hexahedra and prisms in one mesh (mixed face types) |
| `compass.mesh` | `MFEM mesh v1.3` with `attribute_sets` and `bdr_attribute_sets` |
| `tinyzoo-3d.mesh` | one hexahedron, prism, pyramid and tetrahedron |
| `periodic-square.mesh` | discontinuous `L2_T1_2D_P1` nodes (a periodic mesh) |
| `escher-p3.mesh` | legacy `Cubic` nodes, read at the vertices only |
| `amr-quad.mesh` | `MFEM NC mesh v1.0`, which must be refused |
| `star-q2.u.gf`, `star-q2.v.gf` | `H1` order-2 scalar and 2-vector (byNODES) fields |
| `fichera-q2.w.gf` | an `H1` order-2 3-vector field stored byVDIM |
| `compass.t.gf`, `compass.q.gf`, `compass.e.gf` | `H1` order 1, `H1` order 2 on a linear mesh, and `L2` order 0 |

`reference.npz` freezes MFEM's own evaluation: for every vertex, edge, quadrilateral face and element of six of the meshes, keyed by its sorted global vertex ids, the point MFEM's element transformation gives at the entity's reference centre and the value of each grid function there. `test_mfem.py` checks every node of every meshio++ cell against it, so the test needs no MFEM.

## MFEM licence

```text
BSD 3-Clause License

Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```

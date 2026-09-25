<!--pytest-codeblocks:skipfile-->
# Abaqus results-file (`.fil`) fixtures

No free solver writes `.fil`, and no Abaqus licence was available, so the fixtures come from two places.

| File | Source | What it exercises |
|---|---|---|
| `model.fil`, `model_le.fil`, `model_be.fil` | `tools/gen_abaqus_fil_fixtures.py`, record by record from the Abaqus Analysis User's Guide "Results file output format" | one model and two increments in ASCII and in little- and big-endian binary: 1900/1990 continuation, a 1940 long set label, 1932 set continuation, a user element, stresses at integration points and at two shell section points, invariants at the centroid, nodal forces at element nodes, averaged nodal values, `U` and `RF`; every value is a closed-form function repeated in the generator |
| `extras.fil`, `extras_le.fil` | the same generator | the records no public file carries: a contact surface (1501/1502) and contact output (1503-1521), element matrices (1001-1031, one split over two records), a frequency step's two modes (1980), a modal dynamic increment (301/302/304/308, 1999 energies, a rebar and a whole-element record) and an Explicit increment |
| `cjekel/TenBarArea.{fil,dat,inp}` | [Introduction-to-Python-Numerical-Analysis-for-Engineers-and-Scientist](https://github.com/cjekel/Introduction-to-Python-Numerical-Analysis-for-Engineers-and-Scientist), unmodified (MIT, Charles Jekel, `cjekel/LICENSE`) | a real binary run of Abaqus 6.14 (T2D2 trusses) with its `.dat` printout: `U` and `S11` are checked against it |
| `bertoldi/job-strip-angle-buckle-45.fil` | [inverse-design-textile-metamaterials](https://github.com/bertoldi-collab/inverse-design-textile-metamaterials), unmodified (MIT, Bertoldi Group, `bertoldi/LICENSE`; Git LFS) | a real binary buckling run of Abaqus 6.23 (S4R): five modes (1980) |
| `pybaqus/*.fil` | [pybaqus](https://github.com/cristobaltapia/pybaqus) `tests/abaqus/fil/` (commit `0bb36285`), unmodified | real Abaqus 2023 ASCII output (release 6.23-1): C3D8, CPS4, CPE3, and a 2-D model with gaps in its numbering |

The synthetic binary files use 8-byte integers and 4-byte Fortran record markers around each 512-word block. The pybaqus files were read by pybaqus as well as by meshio++; both agree on every node's `U` and every integration point's `S` (see `tests/python/test_abaqus_fil.py`). `TenBarArea` is the roadmap's "matches the `.dat` printout" check, on trusses.

## Licence

pybaqus is MIT-licensed: Copyright (c) 2020 Cristóbal Tapia Camú

```text
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

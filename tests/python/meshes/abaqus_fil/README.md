<!--pytest-codeblocks:skipfile-->
# Abaqus results-file (`.fil`) fixtures

No free solver writes `.fil`, and no Abaqus licence was available, so the fixtures come from two places.

| File | Source | What it exercises |
|---|---|---|
| `model.fil`, `model_le.fil`, `model_be.fil` | `tools/gen_abaqus_fil_fixtures.py`, record by record from the Abaqus Analysis User's Guide "Results file output format" | one model and two increments in ASCII and in little- and big-endian binary: 1900/1990 continuation, a 1940 long set label, 1932 set continuation, a user element, stresses at integration points and at two shell section points, invariants at the centroid, nodal forces at element nodes, averaged nodal values, `U` and `RF`; every value is a closed-form function repeated in the generator |
| `pybaqus/*.fil` | [pybaqus](https://github.com/cristobaltapia/pybaqus) `tests/abaqus/fil/` (commit `0bb36285`), unmodified | real Abaqus 2023 ASCII output (release 6.23-1): C3D8, CPS4, CPE3, and a 2-D model with gaps in its numbering |

The synthetic binary files use 8-byte integers and 4-byte Fortran record markers around each 512-word block. The pybaqus files were read by pybaqus as well as by meshio++; both agree on every node's `U` and every integration point's `S` (see `tests/python/test_abaqus_fil.py`). The "matches the `.dat` printout" check of the roadmap needs a real Abaqus run and is not done.

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

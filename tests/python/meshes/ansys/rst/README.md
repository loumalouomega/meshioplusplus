<!--pytest-codeblocks:skipfile-->
# Ansys MAPDL result file fixtures

**Written by MAPDL, not by meshio++.** They are the external-validation layer of the `ansys_rst` reader: without files the solver wrote, the tests could only prove that meshio++ agrees with itself.

They come unmodified from the [pymapdl-reader](https://github.com/ansys/pymapdl-reader) repository (MIT, notice below; commit `488ad8b`):

| File | Upstream path | What it exercises |
|---|---|---|
| `file.rst` | `ansys/mapdl/reader/examples/file.rst` | a modal analysis of 40 SOLID186 bricks: six modes, frequencies as set times |
| `hex_201.rst` | `tests/testfiles/hex_201.rst` | bit- and windowed-sparse records, node and element components |
| `shell181_2021R1.rst` | `tests/testfiles/shell181_2021R1.rst` | the 2021R1 element-type index (`mapFlag`), int16 element records, rotations |
| `beam44.rst` | `tests/testfiles/beam44.rst` | a release 13 file, whose shorter headers carry 32-bit pointers only |
| `cyc_stress.rst` | `tests/testfiles/rst/cyc_stress.rst` | nodes with rotated coordinate systems |
| `cyclic_v182.rst` | `tests/testfiles/cyclic_reader/cyclic_v182.rst` | a cyclic-symmetry model: the base sector with a warning |
| `file.rth` | `tests/testfiles/file.rth` | a thermal result file: `TEMP` |

`pymapdl_reference.npz` freezes pymapdl-reader's reading of every fixture (its VTK grid, each set's time, node numbers and nodal DOF solution), and `test_ansys_rst.py` compares both meshio++ engines against it: the same cells as sets of node coordinates, bit-identical set times, and nodal solutions to 1e-12. `tools/gen_ansys_rst_reference.py` regenerates it. None of these files holds a result set that covers only some nodes, which pymapdl-reader misreads (it sizes a record's doubles by its length in 4-byte words); the synthetic file `test_ansys_rst.py` writes covers that case instead.

## pymapdl-reader licence

MIT License

Copyright (c) 2021 - 2026 Synopsys, Inc. and ANSYS, Inc. All rights reserved.

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

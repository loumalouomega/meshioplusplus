<!--pytest-codeblocks:skipfile-->
# Ansys MAPDL `.cdb` reference fixtures

**Written by MAPDL, Workbench and HyperMesh, not by meshio++.** They are the external-validation layer of the `ansysInp` reader: without decks other tools wrote, the tests could only prove that meshio++ agrees with itself.

They come unmodified from the [mapdl-archive](https://github.com/akaszynski/mapdl-archive) repository (MIT, notice below; commit `d378d17f`), which inherited most of them from pymapdl-reader. The CRLF line endings two of them carry are kept (`-text` in `.gitattributes`).

| File | Upstream path | What it exercises |
|---|---|---|
| `all_solid_cells.cdb` | `tests/test_data/all_solid_cells.cdb` | SOLID186 hexahedron plus its degenerate wedge, pyramid and tetrahedron; CRLF |
| `etblock.cdb` | `tests/test_data/etblock.cdb` | element types from an `ETBLOCK` table instead of `ET` lines |
| `workbench_193.cdb` | `tests/test_data/workbench_193.cdb` | Workbench's mixed-width node format `(1i7,2i9,6e21.13)`; nodes only |
| `hypermesh.cdb` | `tests/test_data/hypermesh.cdb` | HyperMesh output: `(3i8,6e16.9)`, `(19i8)`, padded `EBLOCK` header; CRLF |
| `ErnoRadiation.cdb` | `tests/test_data/ErnoRadiation.cdb` | SURF152 quads with an extra orientation node; node and element components |
| `mixed_missing_midside.cdb` | `tests/test_data/mixed_missing_midside.cdb` | missing midside nodes (node `0`, and rows cut short) on tetrahedra, pyramids and triangles |
| `HexBeam.cdb` | `src/mapdl_archive/examples/HexBeam.cdb` | SOLID186 hexahedra with node and element components |
| `sector.cdb` | `src/mapdl_archive/examples/sector.cdb` | degenerate SOLID185 wedges; a `CMBLOCK` header with a trailing `!` comment |

`mapdl_archive_reference.npz` freezes mapdl-archive 0.4.2's reading of every fixture (its VTK grid and its component sizes), and `test_ansysInp.py` compares both meshio++ engines against it, cell by cell, as sets of node coordinates. `tools/gen_ansys_cdb_reference.py` regenerates it. In `mixed_missing_midside.cdb`, mapdl-archive places three missing triangle midside nodes at the origin, while meshio++ places them at their edge midpoints; the test allows exactly those three cells to differ.

## mapdl-archive licence

MIT License

Copyright (c) 2023-2024 Alex Kaszynski. All rights reserved.

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

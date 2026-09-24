<!--pytest-codeblocks:skipfile-->
# MSC Nastran HDF5 reference fixtures

**Written by MSC Nastran 2020 (`VERSION=msc20200`, `SCHEMA=20200`), not by meshio++.** They are the external-validation layer of the `nastran_h5` reader: without files a real solver wrote, the tests could only prove that meshio++ agrees with itself.

They come unmodified from the [pyNastran](https://github.com/SteveDoyle2/pyNastran) repository (BSD 3-Clause, notice below). Each one runs the same small model through a different solution sequence:

| File | Upstream path | SOL | What it exercises |
|---|---|---|---|
| `static_elements.h5` | `models/elements/static_elements.h5` | 101 | one static domain, every element card, element results |
| `modes_elements.h5` | `models/elements/modes_elements.h5` | 103 | three modes; SPOINT rows mixed into `EIGENVECTOR` |
| `freq_elements.h5` | `models/elements/freq_elements.h5` | 108 | `_CPLX` tables, five frequencies |
| `modes_complex_elements.h5` | `models/elements/modes_complex_elements.h5` | 107 | complex eigenvectors |
| `buckling_solid_shell_bar.h5` | `models/sol_101_elements/buckling_solid_shell_bar.h5` | 105 | tables covering different domain sets (static + buckling) |
| `time_thermal_elements.h5` | `models/elements/time_thermal_elements.h5` | 159 | scalar `TEMPERATURE`, a card with no cell type (`CHBDYE`) |

The eigenvectors and displacements were checked against pyNastran 1.4.1's own `.h5` reader (`pyNastran/dev/h5`): identical, bit for bit, on every GRID of every domain. `pynastran_reference.npz` freezes that reading (node ids and six components per nodal vector table and domain, 38 tables), so `test_nastran_h5.py` repeats the check without pyNastran installed; `tools/gen_nastran_h5_reference.py` regenerates it. The `.op2` files next to them upstream come from a 2017 run of a different model revision, so they are *not* a reference for these files.

`cord_reference.npz` puts GRIDs of `static_elements.h5` in local systems (a tilted CORD2C, a CORD2S defined in it, a CORD1R through three GRIDs): the coordinates, CP, CD and systems `test_nastran_h5.py` writes into a copy of the file, and pyNastran's basic positions and displacements (`tools/gen_nastran_cord_reference.py`).

## pyNastran licence

Copyright (c) 2011-2026 Steven Doyle.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

  a. Redistributions of source code must retain the above copyright notice,
     this list of conditions and the following disclaimer.
  b. Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the following disclaimer in the
     documentation and/or other materials provided with the distribution.
  c. Neither the name of the pyNastran developers nor the names of any
     contributors may be used to endorse or promote products derived from
     this software without specific prior written permission.


THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS
BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
THE POSSIBILITY OF SUCH DAMAGE.

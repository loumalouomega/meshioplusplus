<!--pytest-codeblocks:skipfile-->
# Nastran OP2 reference fixtures

**Written by MSC Nastran and Simcenter (NX) Nastran, not by meshio++,** except `solid_bending_no_geom.op2` and `static_solid_shell_bar_cord.op2` (below). They are the external-validation layer of the `nastran_op2` reader.

They come unmodified from the [pyNastran](https://github.com/SteveDoyle2/pyNastran) repository's `models/` folder at commit `cb663b5d` (BSD 3-Clause, notice below), copied by `tools/gen_nastran_op2_reference.py`:

| File | Upstream path | Writer | What it exercises |
|---|---|---|---|
| `static_solid_shell_bar.op2` | `sol_101_elements/static_solid_shell_bar.op2` | NX 8.5 | SOL 101: solids, shells, bars, rods; displacement, SPC/MPC forces, loads, stress and strain |
| `mode_solid_shell_bar.op2` | `sol_101_elements/mode_solid_shell_bar.op2` | MSC | SOL 103: three modes |
| `buckling_solid_shell_bar.op2` | `sol_101_elements/buckling_solid_shell_bar.op2` | MSC | SOL 105: a static and four buckling steps in one subcase |
| `solid_bending.op2`, `solid_bending.bdf` | `solid_bending/` | MSC | a tetrahedral solid in bending, and its deck |
| `d173.op2` | `msc/64_bit/d173.op2` | MSC, 64-bit | 8-byte words, CONM2 masses |
| `sol401_tstep1.op2` | `nx/sol401/sol401_tstep1.op2` | NX 2412 | double-precision GRID records, NX's newer CQUAD4 key |
| `time_thermal_elements.op2` | `elements/time_thermal_elements.op2` | MSC | SOL 159: temperatures over time |
| `time_thermal_elements_sort2_nx.op2` | `elements/time_thermal_elements_sort2_nx.op2` | NX | SORT2 tables (skipped) beside SORT1 ones |
| `static_elements.op2` | `elements/static_elements.op2` | MSC | every element card, partial mid-side nodes, PCOMP/PCOMPG, springs and dampers |
| `ctetra10.op2` | `unit/pload4/ctetra10.op2` | MSC | a ten-node tetrahedron |
| `freq_elements2.op2` | `elements/freq_elements2.op2` | NX 10.1 | SOL 111 frequency response, magnitude/phase: complex nodal results, forces, stresses and strains of every element family, NX's von Mises variants |
| `modes_complex_elements.op2` | `elements/modes_complex_elements.op2` | NX 10.1 | SOL 107 complex modes, CLAMA |
| `test_vba.op2` | `nx/test_vba/test_vba.op2` | NX 2206 | random response: PSD, RMS and NO tables; acoustic GRIDs (CD = -1) |

`solid_bending_no_geom.op2` is `solid_bending.op2` with its `GEOM*` and `EPT*` tables removed (every other byte unchanged), and `solid_bending_no_geom.bdf` a copy of its deck: the reader's sibling-deck route.

`static_elements_bgpdt.op2` and `sol401_tstep1_bgpdt.op2` are their files with the `GEOM1` table removed (every other byte unchanged), and `time_thermal_elements_sort2_only.op2` is `time_thermal_elements_sort2_nx.op2` without its SORT1 result tables (`OUGV1`, `OPG1`, `OEF1X`): the basic-grid-point-table and SORT2 routes. `tools/gen_nastran_op2_reference.py` writes all three.

`static_solid_shell_bar_cord.op2` is `static_solid_shell_bar.op2` with the CP and CD words of its GRID record set to the file's own CORD2R/C/S systems (every other byte unchanged), and `cord_reference.npz` pyNastran's basic positions and subcase-1 displacements of it; `tools/gen_nastran_cord_reference.py` writes both (see its docstring for why the displacements use upstream pyNastran's rotation matrices rather than 1.4.1's).

`pynastran_reference.npz` freezes pyNastran 1.4.1's reading of every file: per step the nodal vectors by GRID id, the centre stress and strain by element id, and the multi-valued results (composite plies, plate and solid corners, CBEAM stations, grid point forces) with a `|cols` array naming each value's ply, GRID or station distance, under meshio++'s names. `test_nastran_op2.py` compares against it without pyNastran installed; `tools/gen_nastran_op2_reference.py` regenerates it.

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

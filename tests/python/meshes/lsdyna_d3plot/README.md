<!--pytest-codeblocks:skipfile-->
# LS-DYNA d3plot reference fixtures

Two kinds of family, each in its own folder with its `d3plot` and numbered members:

- **`lasso/`**: [lasso-python](https://github.com/open-lasso-python/lasso-python)'s own test data (`test/test_data/` at commit `723b9303`, BSD 3-Clause, notice below), real LS-DYNA output copied unmodified:

| Folder | What it exercises |
|---|---|
| `simple_d3plot` | 4696 shells, three layers of stress and history variables, element deletion flags, the state in `d3plot01` |
| `d3plot_solid_int` | 16 solids with 8 integration points and 16 shells with 5 layers, plastic strain, deletion flags, 22 states over 23 files |
| `d3plot_beamip` | one beam with 4 integration points |
| `d3plot_node_temperature` | nodal temperatures, 23 states |
| `order_d3plot` | members numbered up to `d3plot100`, read in numeric order |

- **`generated/`**: families `tools/gen_d3plot_reference.py` writes with lasso-python's `D3plot.write_d3plot`: a hexahedron, a wedge and a tetrahedron (degenerate 8-node solids), five quads and a triangle (a degenerate quad), two beams and three parts with titles, over four states in which the plate bends, plastic strain grows, three shells and the tetrahedron are deleted. `shell_solid` holds its states in `d3plot01`, `shell_solid_split` one per file, `shell_solid_double` is the same model in double precision.

`lasso_reference.npz` freezes lasso-python 2.0.4's reading of every family and state (coordinates, ids, times, displacement, velocity, temperature, solid and shell stress, plastic strain and deletion flags). `test_lsdyna_d3plot.py` compares against it without lasso-python installed; `tools/gen_d3plot_reference.py` regenerates it (and the `generated/` families).

## lasso-python licence

Copyright 2022 lasso-open-source

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors
   may be used to endorse or promote products derived from this software without
   specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

<!--pytest-codeblocks:skipfile-->
# Nastran / OptiStruct reference decks

`cylinder.fem` (and `cylinder_cells_first.fem`, the same deck with the elements first) is a HyperMesh-generated Optistruct (Nastran-like) mesh file. It contains the same mesh information as that in `../med/cylinder.med`.

| File | Origin | What it exercises |
|---|---|---|
| `cylinder.fem`, `cylinder_cells_first.fem` | HyperMesh 2017.3 export | a component (`$HMMOVE` + `$HMNAME COMP` after the elements), CBAR/CTRIA3/CQUAD4/CPYRA/CTETRA |
| `optistruct_mixed.fem` | `tools/gen_optistruct_fixture.py` (written by hand, not by meshio++) | small/large/free fields, continuations, CHEXA20/CTETRA10 node order, components by `$HMMOVE` and by property id, SET cards, skipped optimization/contact cards |
| `composite_plate_2022.fem` | pyNastran `models/optistruct/composite_strain_bug/composite_plate_2022.fem`, unmodified (BSD 3-Clause, notice below) | a real OptiStruct deck: free field, a component known only by its property id, `SET` cards named by `$HMSET`, PCOMP/PLY cards skipped quietly |

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

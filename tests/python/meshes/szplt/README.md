<!--pytest-codeblocks:skipfile-->
# Tecplot SZL (`.szplt`) fixtures

Written by `tools/gen_szplt_fixtures.cpp` through TecIO 2018.3 (the TecIO source Tecplot distributes, as vendored by SU2 at commit `bc15466`; SZL file version 105). Every data set is described once and written twice: as `<name>.szplt` through TecIO's writer and as `<name>.dat`, an ASCII Tecplot twin written by the generator itself (TecIO's current writer API writes only SZL).

| Data set | What it holds |
|---|---|
| `fe_mixed` | two FE zones, two tetrahedra and two bricks; `P` nodal, `Q` cell-centred |
| `ordered` | an ordered 3 x 2 x 2 zone; `P` nodal, `Q` cell-centred |
| `transient` | three quad zones on strand 1 at t = 0, 0.5, 1; the later zones share `X`, `Y` and the connectivity of the first, `R` is passive |
| `triangles` | two-dimensional triangles with single-precision data |

The files are ours (MIT). TecIO itself is not in this repository.

<!--pytest-codeblocks:skipfile-->
# OpenRadioss starter-deck fixtures

Written by `tools/gen_radioss_fixtures.py` from the card layouts of OpenRadioss's `hm_cfg_files`. OpenRadioss's own QA decks are CC BY-NC and are not copied.

| File | What it exercises |
|---|---|
| `deck_0000.rad` + `deck_include.inc` | format 2019: every element card the reader reads (degenerate `/BRICK` rows, `/BRIC20`, `/TETRA10`, an inverted `/TETRA4`, `/PENTA6`, `/SHELL` quad and triangle, `/SH3N`, `/BEAM`, `/TRUSS`, a single-node `/SPRING`), parts, a nested `/SUBSET`, every group form the reader resolves and a `/GRBRIC/BOX` it does not, a `/SURF/SEG`, a comma-separated `/NODE` line, an `#include` ended by `#enddata`, and cards after `/END` |
| `old_0000.rad` | format 44: 8- and 16-column fields |
| `gmsh_hex20_0000.rad`, `gmsh_tet10_0000.rad` (+ `.msh` twins) | written by gmsh 4.15 (its own Radioss exporter) for a small box, with the same meshes as gmsh `.msh` files: read back, every cell has the nodes of its `.msh` twin, so the `/BRIC20` table agrees with gmsh's `getVertexRAD` on a file gmsh wrote. Regenerated with `tools/gen_radioss_fixtures.py --gmsh` |

Outside the repository the reader was run on the 81 starter decks of OpenRadioss's `qa-tests/` (the shipped `_0000.rad` files, some without their `#include` files): every one reads, and every cell has a positive volume.

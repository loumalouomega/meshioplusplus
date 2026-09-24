<!--pytest-codeblocks:skipfile-->
# A real OpenRadioss run

`column_0000.rad` (starter) and `column_0001.rad` (engine) are a meshio++ deck, written by hand for these tests: two 1 mm steel bricks stacked in z, the bottom nodes held, the top nodes pushed down at 1 mm/ms, run for 0.004 ms (units kg, mm, ms). The engine asks for a time-history file every 0.0005 ms and an animation file every 0.002 ms, with displacements, velocities, von Mises stress and the brick stress tensors.

OpenRadioss (release `latest-20260728`, `starter_linux64_gf` and `engine_linux64_gf`) ran the deck and wrote:

| File | What it is |
|---|---|
| `columnT01` | the time-history file: 22 global variables, part 1's IE, KE and MASS, TH group 1 (nodes 1, 9 and 12: DX DY DZ VX VY VZ) and TH group 2 (bricks 1 and 2: SX SY SZ), 8 outputs |
| `columnA001`, `columnA002` | the animation files at t = 0 and t ≈ 0.002 ms |
| `columnT01_th_to_csv.csv` | OpenRadioss's own `th_to_csv` reading of `columnT01`: one row per output, the time then every value in file order (its garbled header dropped) |

The deck is ours (MIT); OpenRadioss's licence covers the program, not the files it writes. A single-element-thick column needs `/DT` 0.5: at the default scale the run is unstable.

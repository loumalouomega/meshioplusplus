<!--pytest-codeblocks:skipfile-->
# Real Patran neutral files

Written by P3/PATRAN and PATRAN themselves, redistributed unmodified (renamed) under their projects' permissive licences, whose texts sit beside them.

| File | Written by | Source | Licence |
|---|---|---|---|
| `warp3d_test_25.out` | P3/PATRAN 3.0 (1996) | WARP3D `example_problems_threads/test_25_patran.out`, <https://github.com/rhdodds/warp3d> | University of Illinois/NCSA Open Source License, Copyright (c) 2011 University of Illinois at Urbana-Champaign ([LICENSE.warp3d](LICENSE.warp3d)) |
| `warp3d_test_28.out` | PATRAN 2.5 | WARP3D `example_problems_threads/test_28_patran.out` | same |
| `tahoe_square.pat` | P3/PATRAN 3.0 | Tahoe `benchmark_XML/level.0/geometry/square.1.patran.neutral`, <https://github.com/samanseifi/Tahoe> | BSD-3-Clause, Copyright (c) 2014 Regents of the University of Colorado ([LICENSE.tahoe](LICENSE.tahoe)) |

`tests/python/test_patran.py` reads each with both engines and checks the cell orientation and, for `tahoe_square.pat`, the named components (packet 21).

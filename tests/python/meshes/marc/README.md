<!--pytest-codeblocks:skipfile-->
# MSC Marc fixtures

Written by `tools/gen_marc_fixtures.py` from the layouts of Marc Volume C (program input) and Volume D (PLDUMP2000, the post file), in the styles Marc Mentat writes. No permissively licensed Marc deck or post file exists: the real files the readers were checked on are not copied (see below).

| File | What it exercises |
|---|---|
| `hex20.dat` | two 20-node bricks (type 21) in fixed 5/10-column fields: non-sequential node numbers, a node defined twice, element and node sets with `TO`/`BY` ranges, `C` continuation, set names, `EXCEPT` and `AND`, a data line of another option that starts with a set name, and a history section after `END OPTION` |
| `hex20_extended.dat` | the same model after `EXTENDED`: 10/20-column fields, Mentat's reals with an exponent and no `E` (`5.000000000000000-1`), a negative real touching the node number before it |
| `mixed_free.dat` | free (comma) format: a brick, a brick with collapsed nodes (a wedge), a 10-node tetrahedron, a shell quad, a beam, an element of a type meshio++ does not read (skipped with a warning), and a face set (read as field data) |
| `plane_quad8.dat` | 8-node plane-strain quads (type 27) with two coordinates per node |
| `include_main.dat` + `include_nodes.inc` + `include/sets.inc` | free format, `INCLUDE`d twice (the second relative to the first included file): a Herrmann brick (84) and quad (80) whose pressure nodes have no coordinates, a composite brick (149), a 3-node rebar line (168, middle node second), an element set, an edge set and a face set |
| `remesh.t19` | the model of `results.t19` for its first increment; the second remeshes (flag in block 517, model blocks repeated after block 519) into three bricks, with an edge set (type 12) and a face set (type 13) |
| `results.t19` | a formatted post file of two 8-node bricks: two increments, displacements and reaction forces, the equivalent stress and the stress tensor (codes 17, 311..316) at eight integration points, node and element sets, and a block the reader skips |

Outside the repository the readers were run on:

- the twelve Marc Mentat 2020 decks of DAMASK's element-library tests (`tests/integration/resources/Marc_element_lib/`, AGPL-3.0, not copied): element types 6, 7, 11, 21, 27, 54, 57, 117, 125, 127, 134 and 136 in extended format. Every cell is positively oriented, every mid-edge node sits at its edge's midpoint, and the grain sets partition the elements;
- a `.t19` written by Marc for the FEDES project's small beam example (`models/Example-3D-beam-small/Model1-Input-Marc.t19` in github.com/nightly/fedes, not copied): mesh, sets, loads and stress tensor read as the file describes them, and its von Mises stress (code 17) matches its stress tensor's, which fixes the tensor's component order.

The C++ and Python readers agree bit for bit on all of them.

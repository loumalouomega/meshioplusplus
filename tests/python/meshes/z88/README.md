<!--pytest-codeblocks:skipfile-->
# Z88 fixtures

Z88's files have fixed names, so each deck is a directory. The structure files are written by `tools/gen_z88_fixtures.py` from the Z88 manual's layouts; Z88OS (GPL-2) is not a dependency and none of its example decks is copied.

| Directory | What it exercises |
|---|---|
| `cantilever/` | five hex20 (type 10), a Z88OS v15 header, and the inputs Z88R needs (`z88i2.txt` clamps x = 0 and loads x = 100, `z88mat.txt`, `51.txt`, `z88elp.txt`, `z88int.txt`, `z88i5.txt`) |
| `tets/` | a tet10 (type 16) and a tet4 (type 17) |
| `plate_v13/` | a 2-D quad8 (type 7) and tri6 (type 14) with the Z88 <= V13 / Aurora V1 header and a material line after the elements |
| `polar/` | two tri6 in cylindrical input (`KFLAG = 1`) |

`cantilever/z88o2.txt` and `cantilever/z88o3.txt` are **Z88OS V15's own output** for that deck: the prebuilt `bin/unix64/z88r` of [Z88OS](https://github.com/LSCAD/Z88OS) run as `z88r -c -choly`, with the `z88.dyn`, `z88.fcd` and `z88man.txt` of its distribution (not committed). The tip deflection, 1.849 mm, is within 3 % of Euler-Bernoulli's 1.905 mm for this coarse mesh. Z88R also solved the deck meshio++ writes back from the read mesh, and gave identical displacements.

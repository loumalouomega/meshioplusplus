<!--pytest-codeblocks:skipfile-->
# Z88 fixtures

Z88's files have fixed names, so each deck is a directory. The structure files are written by `tools/gen_z88_fixtures.py` from the Z88 manual's layouts; Z88OS (GPL-2) is not a dependency and none of its example decks is copied.

| Directory | What it exercises |
|---|---|
| `cantilever/` | five hex20 (type 10), a Z88OS v15 header, and the inputs Z88R needs (`z88i2.txt` clamps x = 0 and loads x = 100, `z88mat.txt`, `51.txt`, `z88elp.txt`, `z88int.txt`, `z88i5.txt`) |
| `tets/` | a tet10 (type 16) and a tet4 (type 17) |
| `plate_v13/` | a 2-D quad8 (type 7) and tri6 (type 14) with the Z88 <= V13 / Aurora V1 header and a material line after the elements |
| `polar/` | two tri6 in cylindrical input (`KFLAG = 1`) |
| `frame/` | a 3-D portal of three beams (type 2) and a diagonal truss (type 4), clamped feet and a sway load; beam section values in `z88elp.txt` |
| `plate/` | four quad8 plates (type 20) in a 2-D file, clamped at x = 0, a point load at a corner |
| `shell/` | four tri6 shells (type 24) flat in z = 0 of a 3-D file, clamped at x = 0 and loaded at x = 200 |
| `torus/` | two axisymmetric quad8 (type 8), held axially at z = 0 and pushed outwards |

The `z88o2.txt`, `z88o3.txt` and `z88o4.txt` of `cantilever/`, `frame/`, `plate/`, `shell/` and `torus/` are **Z88OS V15's own output** for those decks: the prebuilt `bin/unix64/z88r` of [Z88OS](https://github.com/LSCAD/Z88OS) run as `z88r -c -choly`, with the `z88.dyn` and `z88.fcd` of its distribution and a `z88man.txt` with von Mises stresses (`ISFLAG 1`) and the beam, plate or shell flag the deck needs (`IBFLAG`, `IPFLAG`, `IHFLAG`); none of these is committed. The cantilever's tip deflection, 1.849 mm, is within 3 % of Euler-Bernoulli's 1.905 mm for this coarse mesh. Z88R also solved the deck meshio++ writes back from the read mesh, and gave identical displacements.

# PLOT3D — multi-block structured grids (`.p3d`, `.x`, `.g`)

The NASA PLOT3D format: multi-block structured grids, with `.q` solution and `.f` function files, what OVERFLOW, CFL3D and the NASA Turbulence Modeling Resource exchange (roadmap 2.4.2). Both engines read and write it (`formats/plot3d.cpp`, `meshioplusplus.plot3d`).

| | |
|---|---|
| **Format name** | `plot3d` |
| **Extensions** | `.p3d`, `.p3dfmt`, `.x`, `.g` (`.xyz` stays the point cloud [`xyz`](./xyz.md); pass `file_format="plot3d"`) |
| **Read / Write** | ✓ / ✓ (lattice blocks only; see [conformance](./conformance.md#plot3d)) |
| **Extra dependencies** | — |

## Reading

```python
import meshioplusplus

grid = meshioplusplus.read("wing.x")
grid = meshioplusplus.plot3d.read("wing.x", solution="wing.q", function="wing.f")
grid = meshioplusplus.plot3d.read("odd.x", precision="single", multiblock=False)
```

A PLOT3D file does not say which variant it is, so the reader infers it: ASCII, Fortran sequential unformatted (record markers) or raw binary (a C stream); single or multi-block (a leading block count); 2-D (`ni nj`, `x y`) or 3-D (`ni nj nk`, `x y z`); single or double precision; with or without an `iblank` array; either byte order. Every combination is checked against the file's exact size (a Fortran file's record sizes). When the size fits several readings, the one whose coordinates are finite and of sane magnitude (no subnormals, which a misread precision produces) wins. A file that still reads several ways is a `ReadError` naming each reading; the read options `precision` (`"single"`/`"double"`), `multiblock`, `iblank` and `dimensions` (2/3) pick one (`Plot3dReadOptions` in C++).

- Each grid block is one cell block of hexahedra (quads in 2-D, lines when only one dimension exceeds one), its points in PLOT3D's order (`i` fastest) after the previous block's, and a cell region `block_<b>`.
- The block dimensions are `field_data["plot3d:dims"]` (`(nblocks, 3)`, `nk = 1` in 2-D); a 2-D grid has `z = 0`.
- `iblank` is the Int32 point data `plot3d:iblank` (overset grids use negative values, so it is not a `uint8`).
- `solution=` reads a `.q` file (`[nblocks] dims`, then per block the four reference values and the `dim + 2` conserved variables): `density`, `momentum` (a 3-vector, `0` for `z` in 2-D) and `energy` point data, and the reference values (Mach number, angle of attack, Reynolds number, time) per block in `field_data["plot3d:q_reference"]`. Its encoding is inferred the same way and its dimensions must match the grid's.
- `function=` reads a `.f` file (`[nblocks]`, per block its dimensions and `nvar`, then the values) as `plot3d:f`, `(npoints, nvar)`.

## Writing

```python
meshioplusplus.write("out.x", grid)
```

The writer emits the Fortran unformatted multi-block double-precision 3-D form (little-endian, 4-byte record markers), with `iblank` when the mesh carries `plot3d:iblank`. It needs `field_data["plot3d:dims"]` whose blocks hold exactly the mesh's points, and refuses any other mesh by name: a general mesh is not a set of lattices. A grid read from that form writes back byte-equal; a 2-D grid is written as 3-D with `nk = 1`. Solutions and function files are not written.

## Verification

The variants are generated in the tests (`tests/python/test_plot3d.py`, `tests/cpp/test_plot3d.cpp`) from the layouts above, in both engines; no NASA-written file has been checked yet (the TMR grids are the planned fixture).

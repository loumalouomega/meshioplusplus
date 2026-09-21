# CalculiX results (`.frd`)

The result file of [CalculiX](http://www.dhondt.de/): the ASCII `.frd` that `ccx` writes and its post-processor `cgx` reads. CalculiX takes Abaqus-style `.inp` input, which meshio++ already writes, so reading its results closes the loop for an open-source solver: mesh → `.inp` → `ccx` → `.frd` → `.vtu` / VTKHDF, with no proprietary tool in the chain.

| | |
|---|---|
| **Format name** | `frd` |
| **Extensions** | `.frd` |
| **Read / Write** | ✓ / — (read-only) |
| **Extra dependencies** | — |

## Reading

```python
import meshioplusplus

mesh = meshioplusplus.read("beam.frd")                 # the first increment
last = meshioplusplus.read("beam.frd", time_step=-1)   # the last one

mesh.point_data["DISP"]       # (n, 3) displacements
mesh.point_data["STRESS"]     # (n, 6) xx yy zz xy yz zx
mesh.field_data["meshio:time"]

# von Mises and principal values beside every stress / strain tensor
mesh = meshioplusplus.frd.read("beam.frd", derived=True)
mesh.point_data["STRESS_mises"], mesh.point_data["STRESS_principal"]

meshioplusplus.write("beam.vtu", mesh)                # or: meshioplusplus convert beam.frd beam.vtu
```

Both engines (the C++ core and the Python reference) read the whole format. Buffers are read by the Python reader. `.frd` is read-only: `meshioplusplus.write(..., file_format="frd")` is an error.

## What is read

The mesh is **the mesh CalculiX wrote, not the mesh you gave it**. `ccx` expands shells and beams into solids before it writes the results, so an `S8` shell comes back as a `hexahedron20` and a `B32` beam as another one, with more nodes than the `.inp` had. Trusses (`T3D2`, `T3D3`) and solids come back as themselves.

| FRD type | Cell type | Nodes | Written by |
|---|---|---|---|
| 1 | `hexahedron` | 8 | `C3D8`, and the expanded `S4`, `B31` |
| 2 | `wedge` | 6 | `C3D6`, and the expanded `S3` |
| 3 | `tetra` | 4 | `C3D4` |
| 4 | `hexahedron20` | 20 | `C3D20`, and the expanded `S8`, `B32` |
| 5 | `wedge15` | 15 | `C3D15`, and the expanded `S6` |
| 6 | `tetra10` | 10 | `C3D10` |
| 7 | `triangle` | 3 | cgx only |
| 8 | `triangle6` | 6 | cgx only |
| 9 | `quad` | 4 | cgx only |
| 10 | `quad8` | 8 | cgx only |
| 11 | `line` | 2 | `T3D2` |
| 12 | `line3` | 3 | `T3D3` |

Types 7 to 10 are the shells of the cgx manual; `ccx` 2.23 never writes them, so they are tested against a hand-written file rather than solver output. An element of any other type is skipped with a warning.

Consecutive elements of one type form one cell block. The element group and material numbers are the integer cell data `frd:group` and `frd:material`.

### Node order

The he20, pe15 and be3 layouts differ from the Abaqus order that meshio++ uses, and are permuted on read (`connectivity[k] = frd_nodes[perm[k]]`):

| FRD type | Permutation | What differs |
|---|---|---|
| 4 (`hexahedron20`) | `0..11, 16..19, 12..15` | the four vertical mid-edge nodes come before the top ring |
| 5 (`wedge15`) | `0..8, 12, 13, 14, 9, 10, 11` | the three vertical mid-edge nodes come before the top ring |
| 12 (`line3`) | `0, 2, 1` | end, mid, end instead of end, end, mid |

The tables were confirmed against `ccx` 2.23: the `.frd` of each single-element deck in the test suite lists its nodes in the permuted order, and after the permutation every mid-node lies on its edge. The Python tests check exactly that. A quadratic **beam** is the exception: `ccx` expands it into a `hexahedron20` whose "mid-edge" nodes are the far end of the beam, which is its own convention and no geometry to check against.

## Results and steps

Each `100C` block of the file is one result (`DISP`, `STRESS`, `TOSTRAIN`, `NDTEMP`, `FORC`, `FLUX`, `ERROR` and whatever else the `-4` record names) at one **increment**. The increments are the steps of the [sequence engine](../sequences.md): `time_step` picks one (`0` is the first, negative counts from the end, out of range is an error), `read_metadata(...).time_values` lists them, and `read_sequence` walks them.

| Field data | Meaning |
|---|---|
| `meshio:time` | the increment's value: the time of a transient run, the load factor of a static one, the **frequency** of a frequency step (cycles per time unit, as `ccx` prints it, not the eigenvalue) |
| `frd:step` | the step counter of the `100C` record (the mode number of a frequency step) |
| `frd:analysis` | the analysis type: 0 static, 1 time step, 2 frequency, 3 buckling, 4 user |

Increments are told apart by their frame id as well as their value, because the modes of a frequency step can share a frequency (the two bending modes of a square beam do).

- **Point data.** One array per result, named as the file names it: one component is `(n,)`, three `(n, 3)`, six `(n, 6)`, more `(n, c)`. A node the block has no value for is NaN. A repeated name within one increment gets a `_2` suffix. Components the file marks as calculated (the `ALL` entry of `DISP`) are not data and are dropped.
- **Symmetric tensors** keep the file's order `xx yy zz xy yz zx` (`SXX SYY SZZ SXY SYZ SZX`), which is also VTK's, so a `.vtu` written from them displays correctly. This is the convention for a six-component symmetric tensor in meshio++; see [mesh data model](../mesh_data_model.md).
- **`points_only`** skips every result block and `arrays` parses only the blocks it names, so reading one field of a large file does not parse the rest.

### Derived fields

`derived=True` adds two arrays beside each six-component `STRESS`, `TOSTRAIN`, `MESTRAIN` or `ZZSTR`:

| Array | Meaning |
|---|---|
| `<NAME>_mises` | `sqrt(½((xx−yy)² + (yy−zz)² + (zz−xx)²) + 3(xy² + yz² + zx²))`, the formula of [ccx2paraview](https://github.com/calculix/ccx2paraview), also for strains |
| `<NAME>_principal` | `(n, 3)`, the eigenvalues in ascending order (min, mid, max) |

Both are NaN where any of the six components is. `derived` is an argument of the format's own reader (`meshioplusplus.frd.read`, `_core.frd_read`, `read_frd(path, options, FrdReadOptions{})` in C++): the generic `read`, the CLIs, the MCP server and the C, Fortran, Julia, R and WASM surfaces do not carry it. Where it is needed there, read the tensor and derive with your own tool for now.

## Record layout

A `.frd` is a stream of fixed-column records keyed by their first columns:

| Key | Record |
|---|---|
| `1C`, `1U` | header and user text (ignored), `1PSTEP` and friends (ignored) |
| `2C` | node block: `-1`, node id, three `E12.5` |
| `3C` | element block: `-1`, element id, type, group, material; then `-2` lines of node ids |
| `100C` | one result block: the frame id, the value, the node count, the analysis type, the step counter and the format flag; then `-4` (name, component count, type), one `-5` per component, `-1` rows and `-3` |
| `9999` | end |

The format flag is 0 for the **short** layout (`I5` ids), 1 for the **long** one (`I10` ids) and 2 for binary. `ccx` writes the long layout. **Values are `E12.5` with no separator**, so negatives run together (`7-1.18144E-06` is a node id and a value) and every field is sliced by column. A result with more than six components continues on `-2` lines with a blank node field.

## Quirks & limitations

- **Binary `.frd` (flag 2) is not read.** It is an error that says so.
- **The `.dat` tabular file is not read.** It orders tensor components differently (`sxx syy szz sxy sxz syz`) and holds integration-point values rather than nodal ones.
- **A second `2C` or `3C` block** (some post-processor exports repeat them) is ignored with a warning; the first wins.
- **The whole file is parsed before a step is chosen**, and `read_metadata` reads it in full (`fell_back_to_full_read` is true): there is no header-only path.
- **No node ids are kept.** Points are numbered in the order of the node block; results are matched by id.

## See also

- [Abaqus](./abaqus.md) — the `.inp` that `ccx` reads.
- [Sequences](../sequences.md) — the steps of a multi-increment file.
- [PVD](./pvd.md) and [VTKHDF](./vtkhdf.md) — the time-series outputs a `.frd` converts to.

# Ansys MAPDL results (`.rst` / `.rth`)

The `.rst` file is the binary results database Ansys MAPDL writes (Ansys Mechanical writes the same file through MAPDL). `.rth` is its thermal counterpart. It holds the model once and then one **result set** per saved load step, substep or mode. meshio++ reads it one set at a time, so an analysis converts to a `.vtu` sequence. The model's input is usually a [`.cdb`](./ansysinp.md) deck.

| | |
|---|---|
| **Format name** | `ansys_rst` |
| **Extensions** | `.rst`, `.rth` (also recognised by content: the standard header of MAPDL file 12) |
| **Read / Write** | ✓ / — (read-only) |
| **Extra dependencies** | — |

## Reading

```python
import meshioplusplus

mesh = meshioplusplus.read("file.rst")                # the first result set
mode = meshioplusplus.read("file.rst", time_step=3)   # the fourth
meshioplusplus.ansys_rst.time_values("file.rst")      # every set's time or frequency
for time, mesh in meshioplusplus.read_sequence("file.rst"):
    ...
```

```bash
meshioplusplus convert file.rst 'mode_{step}.vtu'   # one .vtu per result set
meshioplusplus info file.rst                        # blocks, regions, arrays, times
```

`read` also takes `points_only`/`arrays` (read only the named arrays) and `lenient` (skip elements whose type has no meshio++ cell). Both engines read every set alike.

## File structure

The file is a sequence of Fortran-style records: an `i32` length in 4-byte words, an `i32` flags word, the data and a trailer word. Records point to each other by word offsets from the start of the file, split into low and high 32-bit halves since release 15. The layout is Ansys's `fdresu.inc`, as the open reader [pymapdl-reader](https://github.com/ansys/pymapdl-reader) (MIT) documents it:

- the **standard header** (file number 12, release, title) and the **results header**: counts, the analysis type and pointers to the tables below;
- the **nodal equivalence table** (the node numbers in solution order), the **data-set index** (one pointer per set), and the **time** and **load-step** tables;
- the **geometry**: element types (routine number, key options, node count), nodes (number, coordinates, rotation angles), elements (material, type, real, section, number, nodes) and components;
- per set, a **solution header** (its DOFs and pointers relative to itself) and the **nodal solution**.

Integer, `int16`, `float32` and `float64` records are read, dense or bit- and windowed-sparse. zlib-compressed records (written with `/FCOMP,RST,1` or higher) are refused, naming `/FCOMP,RST,0`. A file must be little-endian, as MAPDL writes on every current platform.

## Mapping

| `.rst` | meshio++ |
|---|---|
| nodes, element types, elements | points and cells, exactly as a [`.cdb`](./ansysinp.md#cell-types) deck's: degenerate shapes resolved, missing midside nodes created at their edge midpoints |
| element routine, type, `MAT`, `REAL`, `SECNUM` | `cell_data["ansys:element"]`, `"ansys:type"`, `"ansys:mat"`, `"ansys:real"`, `"ansys:secnum"` |
| node and element components | `point` and `cell` regions |
| the set's time, or its frequency in a modal analysis | `field_data["meshio:time"]` |
| its load step, substep and cumulative number | `field_data["ansys:load_step"]`, `"ansys:substep"`, `"ansys:cumulative"` |
| `UX UY UZ`, `ROTX ROTY ROTZ`, `AX AY AZ`, `VX VY VZ` | point data `U`, `ROT`, `A`, `V`: 3-component vectors, a missing component 0 |
| other DOFs (`TEMP`, `PRES`, `VOLT`, `MAG`, `CURR` ...) | scalar point data under the DOF's label |

**Rotated nodes.** MAPDL stores a node's solution in its nodal coordinate system, whose axes are the global ones turned about Z by THXY, then about the new X by THYZ, then about the newest Y by THZX. meshio++ rotates each vector DOF to the global axes, `v_global = Rz(THXY) Rx(THYZ) Ry(THZX) v_nodal`. pymapdl-reader applies the three turns in the opposite order, which agrees with this only for nodes rotated about a single axis, as every node in the test fixtures is.

**Nodes without a solution** are NaN: a set can hold only some nodes (with element birth and death, for example), and a record then lists which. MAPDL's undefined value (2^100) is NaN too.

## Not read

- **Element results** (stresses, strains, nodal forces, element energies) and reaction forces. The nodal solution is read; element records are the next step (see the [roadmap](../roadmap.md)).
- **A partial file of a distributed solve** (`file0.rst`, `file1.rst` ...), whose header counts fewer nodes than the model: a `ReadError` asks for the combined file, which MAPDL writes by default or with `RESCOMBINE`.
- **Cyclic-symmetry expansion**: a cyclic model's base sector is read, with a warning; its other sectors are not generated.
- A results file whose geometry changes between sets (rezoning) is read with the geometry its results header points to, for every set.

## Notes

- **Validation.** The test fixtures in `tests/python/meshes/ansys/rst/` are result files from pymapdl-reader's repository, written by MAPDL releases 13 to 2021 R1: static, modal and thermal, solids, shells and beams, sparse records, rotated nodes and a cyclic model. Both engines read the cells, set times and nodal solutions pymapdl-reader reads, to 1e-12 (bit-identical in practice), frozen in `pymapdl_reference.npz` by `tools/gen_ansys_rst_reference.py`. Across all 23 files in that repository that are not distributed partials, the cells agree except for elements meshio++ skips under `lenient` (FOLLW201 follower elements, TARGE170 pilot nodes) and 740 CONTA174 contact elements of one file, which pymapdl-reader leaves empty and meshio++ reads as quadrilaterals whose midside nodes lie on their edges. The nodal solutions agree except on result sets that hold only some nodes, where pymapdl-reader reads past the record (it sizes a record's doubles by its length in 4-byte words). A synthetic file in the tests pins the three-axis rotation, a partial set and the undefined value.

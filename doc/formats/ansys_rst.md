# Ansys MAPDL results (`.rst` / `.rth`)

The `.rst` file is the binary results database Ansys MAPDL writes (Ansys Mechanical writes the same file through MAPDL). `.rth` is its thermal counterpart. It holds the model once and then one **result set** per saved load step, substep or mode. meshio++ reads it one set at a time, so an analysis converts to a `.vtu` sequence: the nodal solution, the reaction forces and the element results (nodal stresses, strains and forces). The files of a distributed solve are read together, and a static cyclic-symmetry model can be expanded to its full rotor. The model's input is usually a [`.cdb`](./ansysinp.md) deck.

| | |
|---|---|
| **Format name** | `ansys_rst`; `ansys_rst_cyclic` for the full rotor of a cyclic model (by name only) |
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
mesh.point_data["S"]      # stresses averaged at the corner nodes
mesh.cell_data["S"][0]    # per element node: (cells, nodes, 6), NaN at midside nodes

meshioplusplus.read("file0.rst")                      # a distributed solve: file0..fileN.rst
rotor = meshioplusplus.read("sector.rst", file_format="ansys_rst_cyclic")
```

```bash
meshioplusplus convert file.rst 'mode_{step}.vtu'   # one .vtu per result set
meshioplusplus info file.rst                        # blocks, regions, arrays, times
meshioplusplus convert --input-format ansys_rst_cyclic sector.rst rotor.vtu
```

`read` also takes `points_only`/`arrays` (read only the named arrays) and `lenient` (skip elements whose type has no meshio++ cell). Both engines read every set alike.

## File structure

The file is a sequence of Fortran-style records: an `i32` length in 4-byte words, an `i32` flags word, the data and a trailer word. Records point to each other by word offsets from the start of the file, split into low and high 32-bit halves since release 15. The layout is Ansys's `fdresu.inc`, as the open reader [pymapdl-reader](https://github.com/ansys/pymapdl-reader) (MIT) documents it:

- the **standard header** (file number 12, release, title) and the **results header**: counts, the analysis type and pointers to the tables below;
- the **nodal equivalence table** (the node numbers in solution order), the **data-set index** (one pointer per set), and the **time** and **load-step** tables;
- the **geometry**: element types (routine number, key options, node count), nodes (number, coordinates, rotation angles), elements (material, type, real, section, number, nodes) and components;
- per set, a **solution header** (its DOFs and pointers relative to itself), the **nodal solution**, the **reaction forces** (a table of node and DOF positions and their values) and the **element solution**: for each element a table of pointers, one per kind of record (`EMS ENF ENS ENG EGR EEL EPL ECR ETH EUL ...`), relative to itself.

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
| reaction forces of `UX UY UZ` and `ROTX ROTY ROTZ` | point data `RF` and `RMOM`, rotated to the global axes; NaN at nodes with none |
| reactions of other DOFs (`TEMP` ...) | scalar point data `RF_<DOF>` (`RF_TEMP` is the reaction heat flow) |
| element nodal stresses (`ENS`) | `S`: cell data per element node `(cells, nodes, 6)` and point data `(points, 6)` averaged at the corner nodes, `xx yy zz xy yz xz` |
| elastic, plastic, creep and thermal strains (`EEL`, `EPL`, `ECR`, `ETH`) | `EPEL`, `EPPL`, `EPCR`, `EPTH`, laid out as `S` (engineering shear strains, as MAPDL stores them) |
| the top surface of a layered shell (SHELL181/281 with `KEYOPT(8)=0`) | `S@top`, `EPEL@top` ... beside the bottom surface's `S`, `EPEL` ... |
| element nodal forces (`ENF`) | cell data `ENF`, `(cells, nodes, DOFs)` in the set's DOF order |

**Rotated nodes.** MAPDL stores a node's solution in its nodal coordinate system, whose axes are the global ones turned about Z by THXY, then about the new X by THYZ, then about the newest Y by THZX. meshio++ rotates each vector DOF to the global axes, `v_global = Rz(THXY) Rx(THYZ) Ry(THZX) v_nodal`. pymapdl-reader applies the three turns in the opposite order, which agrees with this only for nodes rotated about a single axis, as every node in the test fixtures is.

**Nodes without a solution** are NaN: a set can hold only some nodes (with element birth and death, for example), and a record then lists which. MAPDL's undefined value (2^100) is NaN too.

## Element results

An element's stresses and strains are stored at its corner nodes (`nodstr` of them, the corners of a quadratic element), in the element's coordinate system. meshio++ rotates each tensor to the global axes by the element's Euler angles (`EUL`, one set for the element or one per node, MAPDL's 3-1-2 convention) and places it on the cell's node the record's node became, degenerate shapes included. Midside nodes carry no value (NaN).

- **Cell data** keeps each element's own values, discontinuous across elements: `cell_data["S"]` has one `(cells, nodes, 6)` array per block.
- **Point data** averages them at each corner node over the elements that have a value there, as pymapdl-reader's `nodal_stress` and `nodal_elastic_strain` do (their seventh item, the equivalent strain, is not kept). Midside nodes are NaN.
- A negative entry `-n` in an element's pointer table is a record of `n` zeros that MAPDL did not write (unloaded elements of a cyclic model do this): it reads as zeros and counts in the average.
- **Layered shells.** SHELL181 and SHELL281 with `KEYOPT(8)=0` store the bottom surface, then the top: the bottom is `S`, `EPEL` ..., the top `S@top`, `EPEL@top` ... (NaN in cells of other types).
- Before release 14.5 a node's `ENS` record holds 11 items (the six components, three principal stresses, the intensity and the equivalent stress); the first six are read.
- **Nodal forces** (`ENF`) are per element node, for every node the element has forces at, in the set's DOF order. At each node they sum to minus the reaction force or the applied load.
- Line and point elements (beams, links, masses) carry no element results in meshio++: their records hold section forces and stresses that are not tensors at nodes.

`arrays=[...]` narrows what is read by name (`S`, `S@top`, `EPEL`, `ENF`, `RF` ...), and `points_only` skips all results.

## Distributed solves

A distributed MAPDL solve can leave one results file per process, `<job>0.rst` to `<job>N.rst`, each with the nodes and elements of its part of the model. Reading the main file, `<job>0.rst`, reads the others beside it: their models are merged by node number (the nodes on the partition interfaces appear in several files), each element keeps the results of the file that holds it, and reactions at shared nodes are summed. A partial file other than the main one is refused, naming the main one; so is a main file whose partial files do not cover the global node list (one is missing). MAPDL writes the combined file itself by default, or with `RESCOMBINE`; both give the same mesh and results.

## Cyclic symmetry

A cyclic-symmetry model holds one sector (the base sector, element numbers up to `csEls`) and, before release 18.2, a duplicate sector. `ansys_rst` reads the file as stored, with a warning. `ansys_rst_cyclic` reads a **static** cyclic analysis as the full rotor: the base sector's cells and points are repeated round the cyclic axis once per sector (global Z, or the Z axis of the local coordinate system the model names), and every vector and tensor result is rotated with its sector; `cell_data["ansys:sector"]` numbers the copies and `field_data["ansys:sectors"]` counts them. Coincident nodes on the sector boundaries are not merged (`meshioplusplus clean` merges them). Regions hold their members in every sector. A modal cyclic analysis is refused: expanding its harmonic indices needs the paired modes combined, which is not implemented.

## Not read

- **zlib-compressed records** (`/FCOMP,RST,1` or higher) are refused, naming `/FCOMP,RST,0`: no compressed file was available to check the record framing against.
- **Other element records**: energies (`ENG`), miscellaneous items (`EMS`, `EMN`), fluxes and gradients (`EFX`, `EGR`), nonlinear data (`ENL`), state variables (`ESV`), contact data (`ECT`).
- **Modal cyclic expansion** (above), and the full rotor of any analysis other than a static one.
- A results file whose geometry changes between sets (rezoning) is read with the geometry its results header points to, for every set.

## Notes

- **Validation.** The test fixtures in `tests/python/meshes/ansys/rst/` are result files from pymapdl-reader's repository, written by MAPDL releases 13 to 2021 R1: static, modal and thermal, solids, shells and beams, sparse records, rotated nodes, two cyclic models and a distributed solve with its combined file. Both engines read the cells, set times and nodal solutions pymapdl-reader reads, to 1e-12 (bit-identical in practice), frozen in `pymapdl_reference.npz` by `tools/gen_ansys_rst_reference.py`. The averaged stresses and elastic strains match pymapdl-reader's `nodal_stress` and `nodal_elastic_strain` exactly on every solid model (a static SOLID186 beam, a modal one, the distributed solve, the cyclic models), and the reactions match `nodal_reaction_forces` once pymapdl-reader's values, which it leaves in the nodal coordinate systems, are rotated. Its element results of SHELL181/281 models are not compared: pymapdl-reader corrupts its heap reading them; nor are the stresses of the release 13 file, which it does not read correctly (its `element_stress` refuses files before 14.5; meshio++'s reading of them is checked by the principal stresses matching the eigenvalues of the six components). The distributed solve merges to the combined file MAPDL wrote, and the full rotors of both static cyclic models match pymapdl-reader's to 1e-15. The element nodal forces of a static SOLID186 model sum at every node to minus its reaction or its applied load. A synthetic file pins a rotated element, a layered shell, an all-zero record and the reactions of a rotated node. Across all 23 files in that repository that are not distributed partials, the cells agree except for elements meshio++ skips under `lenient` (FOLLW201 follower elements, TARGE170 pilot nodes) and 740 CONTA174 contact elements of one file, which pymapdl-reader leaves empty and meshio++ reads as quadrilaterals whose midside nodes lie on their edges. The nodal solutions agree except on result sets that hold only some nodes, where pymapdl-reader reads past the record (it sizes a record's doubles by its length in 4-byte words). A synthetic file in the tests pins the three-axis rotation, a partial set and the undefined value.

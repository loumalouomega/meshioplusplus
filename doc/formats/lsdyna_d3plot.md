# LS-DYNA state database (`d3plot`)

The binary result database [LS-DYNA](https://lsdyna.ansys.com/) writes for post-processing, and what LS-PrePost, lasso-python and most public crash datasets start from. A run writes the family `d3plot`, `d3plot01`, `d3plot02`...: the geometry and then one state per output time. meshio++ reads the family as a sequence, `d3plot` → `.vtu` / VTKHDF / any meshio++ format, with every state as a step. The input deck (`.k`) is a different format: see [LS-DYNA keyword decks](./lsdyna.md).

| | |
|---|---|
| **Format name** | `lsdyna_d3plot` |
| **Extensions** | none: found by its file name, `d3plot` (any case), or by its control block |
| **Read / Write** | ✓ / — (read-only) |
| **Extra dependencies** | none |

## Reading

```python
import meshioplusplus

mesh = meshioplusplus.read("run/d3plot")                  # the first state
last = meshioplusplus.read("run/d3plot", time_step=-1)    # the last one

last.point_data["displacement"]             # (n, 3), relative to the initial geometry
last.cell_data["effective_plastic_strain"]  # per block: (cells, points)
last.cell_data["stress"]                    # per block: (cells, points * 6)
last.field_data["lsdyna_d3plot:layout:stress"]   # [points, 6]
last.cell_data["lsdyna:alive"]              # int8: 1 active, 0 deleted

meshioplusplus.lsdyna_d3plot.time_values("run/d3plot")   # the time of every state
```

```bash
meshioplusplus convert run/d3plot 'state_{step}.vtu'    # one .vtu per state
meshioplusplus convert run/d3plot run.pvd               # or a ParaView collection
```

Both engines read the whole format: the C++ core, and a Python reference reader the core falls back to. `d3plot` is read-only.

## Which files

- **The family.** Open the base file, `d3plot`. The files named after it with a number (`d3plot01`, `d3plot02`... `d3plot100`) are found beside it and taken in numeric order. Together they are one stream of states: a state never straddles two files, and zero padding at the end of a file is skipped. As in lasso-python, a file holds as many whole states as fit before its last non-zero word (lasso-python counts bytes, which loses a state of a big-endian file), and that word is the end mark (-999999) LS-DYNA writes after the last state. A numbered member opened on its own is refused with a message naming the base file, and a sequence glob such as `run/d3plot*` keeps the base file only.
- **Precision and byte order.** Single (4-byte words) and double precision (8-byte words) runs, in either byte order, are told apart by the control block: its file-type word and its dimension word must be plausible in exactly one reading.
- **File types.** `d3plot` and `d3eigv` (the eigenvector database, whose states are the modes) are read. `d3part` and interface-force files, femzip-compressed files (`d3plot.fz`, or the femzip mark in the header), two-dimensional databases (`NDIM = 2`), CFD and multi-solver data and adaptive-remeshing runs (`NADAPT != 0`, whose geometry changes between states) are refused with a `ReadError` that names them.

## The model

- **Points** are the initial coordinates, in the file's node order. The user node ids are `point_data["lsdyna:nid"]` (the internal numbers 1..n when the file has no numbering section).
- **Cells.** One cell block per element family and cell type, in the order solids, thick shells, beams, shells, and within a family in order of first appearance:

| Family | Cell types |
|---|---|
| 8-node solids | `hexahedron`, and `wedge`, `pyramid`, `tetra` from repeated nodes (the [keyword reader's](./lsdyna.md) collapse rules, `detail/degenerate_solid.hpp`) |
| 10-node solids (`NEL8 < 0`) | `tetra10`: the 8-node connectivity plus the two extra nodes, in LS-DYNA's order, which is meshio++'s |
| thick shells | as the 8-node solids |
| beams | `line` (the orientation node is dropped) |
| shells | `quad`, `triangle` when `n4 = n3`, `quad8` for the 8-node shells (`NEL48`) |

- **Cell data.** `lsdyna:eid` (the user element id) and `lsdyna:part` (the user part id), integer, per block.
- **Regions.** One [cell region](../regions.md) per part with elements, tagged with the user part id and named by the part title (the `90001` titles section), or `Part <id>` without one.
- **Part ids.** With a numbering section that lists them, the element's part index maps through its first part-id array, as lasso-python reads it; otherwise through the titles section's ids when it lists every part; otherwise the part number itself.

## States and data

Every state is a step of the [sequence engine](../sequences.md): `time_step` picks one (`0` is the first, negative counts from the end), `read_metadata(...).time_values` lists the times, and `read_sequence` walks them. Only the chosen state is read.

| Field data | Meaning |
|---|---|
| `meshio:time` | the state's time |
| `lsdyna:state` | its index in the family |
| `global_kinetic_energy`, `global_internal_energy`, `global_total_energy`, `global_velocity` | the global variables |
| `part_internal_energy`, `part_kinetic_energy`, `part_velocity`, `part_mass`, `part_hourglass_energy` | one value (or vector) per part, when the file has them |

**Point data**, whichever the header flags: `displacement` (the current coordinates minus the initial ones; the points themselves stay where the geometry put them, as in [FEBio plot files](./xplt.md)), `velocity`, `acceleration`, `temperature` (3 components with temperature layers), `heat_flux`, `mass_scaling`, `temperature_gradient`, `residual_forces`, `residual_moments`.

**Cell data.** A variable written per integration point (solids: 1 or 8) or through-thickness layer (shells and thick shells: `MAXINT`; beams: their integration points) is an array of shape `(cells, points * components)`, point-major, with `field_data["lsdyna_d3plot:layout:<name>"] = [points, components]`: the same convention as the [Abaqus `.fil`](./abaqus_fil.md) reader. A family with fewer points than the widest one, or without the variable, is NaN-padded. Point `k` means integration point `k` for a solid and layer `k` for a shell (for three layers: mid, inner, outer).

| Cell data | Families |
|---|---|
| `stress` (6: xx, yy, zz, xy, yz, zx), `effective_plastic_strain` | solids, thick shells, shells; beams' plastic strain per integration point too |
| `history_variables` | solids (`NEIPH`), shells and thick shells (`NEIPS`), beams (`NEIPB`) |
| `strain` (6, per integration point) | solids (`ISTRN`) |
| `strain_inner`, `strain_outer` (6 each) | shells and thick shells (`ISTRN`) |
| `plastic_strain_tensor`, `thermal_strain_tensor` | solids and shells, when `IDTDT` flags them |
| `thickness`, `internal_energy`, `shell_bending_moment`, `shell_shear_force`, `shell_normal_force`, `shell_element_variables` | shells |
| `beam_axial_force`, `beam_shear_force`, `beam_bending_moment`, `beam_torsion_moment`, `beam_axial_stress`, `beam_shear_stress`, `beam_axial_strain` | beams |
| `thermal_variables` | solids (`NT3D`) |

Shells in rigid parts (material type 20, listed in the material-type section of `NDIM = 5` and `7` files) have no state data: their values are NaN.

**Deletion.** The deletion flags are an int8 mask, `lsdyna:alive`: 1 for an active element or node, 0 for a deleted one. With `MDLOPT = 2` (element deletion) it is cell data; with `MDLOPT = 1` (node deletion) it is point data. Deleted cells stay in the mesh, so every state has the same cells; a consumer that wants them gone filters on the mask.

## What is not read

- **Rigid bodies, rigid roads, SPH particles and airbag particles.** Their geometry and state words are skipped (their sizes are computed so that everything after them still reads), with a warning for SPH and airbags.
- **Rigid-wall forces**, which follow the per-part values in the globals.
- **The sibling files** `d3part`, `d3thdt` and `binout` (an LSDA container), and interface-force files.
- **20-, 27-node and cubic solids** (`NEL20`, `NEL27`, `NEL21P`, `NEL15T`, `NEL20T`, `NEL40P`, `NEL64`): refused by name.

## Verification

The word offsets and the state layout follow the LS-DYNA database manual, as [lasso-python](https://github.com/open-lasso-python/lasso-python) (BSD 3-Clause) reads it. The test suite reads lasso-python's five real test families (a shell model with deletion flags, a shell and solid model with 22 states, a beam with integration points, a thermal model and a family split over `d3plot` ... `d3plot100`) and three shell + solid + beam families written by lasso-python's own `D3plot.write_d3plot` (in single and double precision, in one file and in one file per state), whose shells and a solid are deleted state by state. On every state of every family the displacement, the solid and shell stress, the plastic strain and the deletion flags equal lasso-python 2.0.4's reading, and the two engines agree exactly on every array.

Two deviations from lasso-python, both toward the manual: the two extra nodes of ten-node solids are read once, after the 8-node connectivity (lasso-python reads a second copy, which its writer also writes and which is skipped here when present), and the displacement is computed in double precision. No file LS-DYNA itself wrote with element deletion events, rigid bodies, SPH or airbags has been read; those paths are checked only against the manual and lasso-python.

# OpenRadioss animation files (`<run>A001`, `A002`, …)

The animation files an OpenRadioss (or Radioss) run writes, one per output time: the model as it is at that time and the results the run was asked to animate (`/ANIM` in the engine deck). Each file is one state, so a run's files together are a transient [sequence](../sequences.md). The starter deck that defines the model is read by [`radioss`](./radioss.md).

| | |
|---|---|
| **Format name** | `radioss_anim` |
| **File names** | a stem, then `A` and three or more digits (`crashA001`, `crashA042`); also recognised by content, the magic number `0x542C` |
| **Read / Write** | ✓ / ✓ ([conformance](../conformance.md#radioss-anim)) |
| **Extra dependencies** | — |

## Reading

```python
import meshioplusplus

state = meshioplusplus.read("crashA010")                  # or meshioplusplus.radioss_anim.read(...)
for time, mesh in meshioplusplus.read_sequence("crashA*"):
    ...                                                    # every state, in time order
```

```bash
meshioplusplus convert 'crashA*' crash.vtkhdf              # the whole run as one time series
```

`read` takes no options. Both engines (the C++ core and the pure-Python reference) read the same meshes. The state's time is `field_data["meshio:time"]` and `read_metadata(...).time_values`, so a sequence over a run's files takes each file's own time rather than the number in its name.

## Writing

```python
meshioplusplus.write("runA001", mesh)                              # one state
meshioplusplus.write_sequence("runA001", steps, file_format="radioss_anim")   # runA001, runA002, ...
```

```bash
meshioplusplus convert model.vtu runA001                             # the file name decides the format
```

Both engines write the same bytes. A mesh is one state, its time `field_data["meshio:time"]` (0 when absent); a series is one file per step, numbered from the number in the first name (`runA001`, `runA002`, …; a name that is not animation-shaped is a stem, `run` giving `runA001`). Vertices become SPH particles, lines 1-D elements, triangles and quads facets (a triangle repeats its last node), and tetrahedra, pyramids, wedges and hexahedra 8-node bricks with repeated nodes, the patterns the reader collapses again. Elements are grouped by family and part, so the blocks come back in the file's order. Nodal arrays of one or three components are scalars and vectors; cell arrays of one component are scalars, six components (xx yy zz xy yz zx) tensors (a facet keeps xx, yy and xy), nine components 1-D force and moment sets. `radioss:part` or the cell regions name the parts; `radioss:alive`, `radioss:element_id`, `radioss:node_id`, `radioss:mass`, `radioss:material` and `radioss:property` are written when present (a missing partner array is written as zeros or a running number). Everything numeric is stored as a 32-bit float, so values come back to that precision. The **time-history lists** are regions: each item of the file (the 0-based index of a node, or of an element of one family, and a label of up to 49 characters) is an entry of the region `radioss_th:<label>`, a `point` region for the nodes and a `cell` region for the elements (a label that names elements of several families is one region), and the writer sets flag 6 and writes the regions named `radioss_th:` back as the four lists; they are not parts. The indices are read as the file's own order within the family (the fixtures' generator and the layout in the OpenRadioss reader agree on a 0-based position; no run of OpenRadioss's `anim_to_vtk`, which skips the lists, checks it). Other cell types and arrays of other shapes are dropped with a warning; the older layouts are not written. The **nodal normals** (three 16-bit words per node, between the facet part table and the scalar names) are skipped on reading and written as zeros: their scale and sign are not documented in any source available here and `anim_to_vtk` does not use them, so a writer that invented an encoding would be unchecked.

## The file

A big-endian binary file with the magic `0x542C`: the time, three titles and ten flags, then

| Section | When | Holds |
|---|---|---|
| 2-D | always | the nodes (coordinates, nodal scalars and vectors) and the 4-node facets of shells, with element scalars and 2-D tensors |
| 3-D | flag 3 | 8-node bricks with element scalars and tensors |
| 1-D | flag 4 | 2-node elements (beams, springs, rigid-body links) with element scalars and force/moment sets |
| masses, ids | flags 1, 2 | element and nodal masses; node and element ids |
| hierarchy | flag 5 | each part's subset, material and property; subset, material and property names |
| time history | flag 6 | the nodes and elements saved for time history, as the `radioss_th:<label>` regions |
| SPH | flag 8 | particles with scalars and tensors |

The layout is the one OpenRadioss's `anim_to_vtk` converter reads (MIT); no code is copied. Older layouts (magic numbers `0x5426` to `0x542B`, from Radioss before OpenRadioss) are refused: neither `anim_to_vtk` nor OpenRadioss's own readers document them, and no such file has been found. Every file OpenRadioss writes, including the test suite's own run (`tests/python/meshes/radioss_th/column`), is `0x542C`.

## The mesh

- **Points** are the state's coordinates (the deformed shape), as 64-bit floats.
- **Cells**, in the file's family order: 1-D elements as `line`, facets as `quad` (or `triangle` when a node repeats), bricks as `hexahedron` or the `tetra`, `pyramid` or `wedge` a brick with repeated nodes stands for (as the [starter reader](./radioss.md) collapses them), SPH particles as `vertex`. One block per type, in order of first appearance.
- **Point data**: every nodal scalar and vector under its name (`Velocity`, `Displacement`, `Contact Forces`…), `radioss:node_id` and `radioss:mass` when the file has them.
- **Cell data**: every element scalar under its name, the union over the families, NaN on the cells of a family that does not have it (a 2-D and a 3-D `Plastic Strain` are one array). Tensors have six components `xx yy zz xy yz zx`; a 2-D tensor's `zz`, `yz` and `zx` are NaN. A 1-D force/moment set has nine components, NaN off the 1-D cells. `radioss:part`, and `radioss:alive` (1 while the element exists, 0 once it is deleted), always; `radioss:element_id`, `radioss:mass`, `radioss:material`, `radioss:property` when every family has them.
- **Regions**: each part is a `cell` region named by its title (tag = part id). A 1-D and a 2-D part may share an id.

`anim_to_vtk` writes a facet's alive flag, 0xFF, as its "erosion status" 0 and the bricks' 1 as 1; meshio++ reads every non-zero flag as alive. `anim_to_vtk` also places a 3-D tensor's fifth and sixth values at `xz` and `yz`, where its own reading code documents them as `yz` and `zx`; meshio++ keeps the file's order.

## Validation

The fixtures under `tests/python/meshes/radioss_anim/` are written by `tools/gen_radioss_anim_fixtures.py` (two states of a small model with every section), and `anim_to_vtk/` holds what `anim_to_vtk` makes of them; the tests compare every array with it. Outside the repository both engines were also run on three animation files OpenRadioss wrote (Kitware's `openradioss-to-vtkhdf` test data, bricks, shells, beams and SPH): they agree with each other and with `anim_to_vtk` to its printed precision.

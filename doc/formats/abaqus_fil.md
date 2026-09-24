# Abaqus results file (`.fil`)

The results file Abaqus/Standard (and Explicit, with `*FILE OUTPUT`) writes on request with `*NODE FILE` and `*EL FILE`: the one route into Abaqus results that needs neither the ODB API nor an Abaqus installation. For `.odb` see the [roadmap](../roadmap.md) (an exporter script is planned). The input deck itself is the [`abaqus`](./abaqus.md) format.

| | |
|---|---|
| **Format name** | `abaqus_fil` |
| **Extensions** | `.fil` (also recognised by content, ASCII and binary) |
| **Read / Write** | ✓ / — |
| **Extra dependencies** | — |

## Reading

```python
import meshioplusplus

mesh = meshioplusplus.read("job.fil")                        # the first increment
last = meshioplusplus.abaqus_fil.read("job.fil", time_step=-1)
times = meshioplusplus.abaqus_fil.time_values("job.fil")
series = meshioplusplus.read_sequence("job.fil")             # every increment
```

`read(filename, points_only=False, arrays=None, time_step=0)`. Both engines (the C++ core and the pure-Python reference) read the same meshes.

## Records

A `.fil` file is a sequence of records `[length, key, attributes…]` of 8-byte words.

- **Binary** (the default): 512-word blocks, each one Fortran record, in the writing machine's byte order (both orders and 4- or 8-byte markers are recognised, and files without markers are read too). A record may cross block boundaries; the last record of an increment is padded to the end of its block.
- **ASCII** (`*FILE FORMAT, ASCII`): each record starts with `*`; items are `I` + a two-digit digit count + the digits, `D` (or `E`) + a 22-character real, `A` + 8 characters, run together on 80-column lines. Items are read by width, so a `*` inside a text item does not start a record.

| Key | Read as |
|---|---|
| 1900 / 1990 | an element and its continuation: a cell, its label kept as `cell_data["abaqus:id"]` |
| 1901 | a node, its label kept as `point_data["abaqus:id"]` |
| 1931 / 1932, 1933 / 1934 | a node set / element set and its continuation: a `point` / `cell` region |
| 1940 | a label cross-reference: a set name longer than 8 characters is written as a number, resolved here (assembly sets keep Abaqus's `ASSEMBLY_INSTANCE_SET` spelling) |
| 2000 / 2001 | the start and end of an increment: a step |
| 1911 | an output request: nodal or element records follow |
| 1 | the element header (element, integration or node point, section point, location) of the element records after it |
| 1921, 1922, 1902 | release, heading, active degrees of freedom: skipped |
| contact (15xx), element matrices (10xx), modal generalised (301–310), substructure | skipped with a warning |

Element types map through the `.inp` table and, for the ones it does not list, by family and node count (`C3D…`, `CPE…`/`CPS…`/`CAX…`, shells `S…`, membranes `M3D…`, trusses and beams); their nodes are in meshio++'s order. Types with no meshio++ cell (user elements, connectors, springs) are skipped with a warning.

## Steps and data

Every increment is a step: `time_step` selects one (0 = first, negative counts from the end), and `field_data` carries `meshio:time` (the total time), `abaqus:step`, `abaqus:increment`, `abaqus:step_time` and `abaqus:procedure`. The increment's records become data named by the Abaqus identifier (`U`, `V`, `A`, `RF`, `CF`, `COORD`, `NT`, … for nodes; `S`, `E`, `SINV`, `PE`, `LE`, `PEEQ`, `MISES`, `ENER`, `NFORC`, … for elements; `key_<n>` for a key without one):

| Where | Becomes |
|---|---|
| nodal record | `point_data`, NaN for nodes without a value |
| element record at the integration points (location 0) | `cell_data` of shape `(cells, points × components)`, point-major |
| at the centroid (1) or for the whole element (5) | `cell_data` of shape `(cells, components)` |
| at the element nodes (2) | `cell_data` of shape `(cells, nodes × components)`, point-major in the cell's node order |
| averaged at the nodes (4) | `point_data` |
| rebar (3) | skipped |

Every array is rectangular and has the same width in every cell block, so every writer can hold it (VTU, XDMF, …): the widest block sets the number of points and components, and the other blocks are padded with NaN. Per-point data is flattened point-major, column `point × components + component`; its `(points, components)` is `field_data["abaqus:layout:<name>"]`, so `a.reshape(len(a), *layout)` restores it. A single column drops its axis. Section points above 1 (shell and beam layers; continuum elements write 0) get `@sp<k>` appended to the name: `S` holds the first layer and the solids, `S@sp5` the fifth layer. A block with no value for a field holds NaN.

Components keep the file's order, as every reader does ([mesh data model](../mesh_data_model.md)): a solid's `S` is 11, 22, 33, 12, **13, 23**, where meshio++'s six-component convention is `xx yy zz xy yz zx`. Von Mises does not depend on the order of the shear components; for principal values, swap the last two first. Von Mises is also `SINV`'s first component when the run wrote invariants; otherwise compute it from `S` with the [tensor-invariants operation](../tensor_invariants.md) on a `(cells, 6)` array, such as the mean over the points.

## Validation

No free solver writes `.fil` and no Abaqus licence was available. The fixtures under `tests/python/meshes/abaqus_fil/` are real Abaqus 2023 ASCII output from pybaqus's test suite (MIT), whose `U` and integration-point `S` match pybaqus's own reading on all 236 values, and synthetic files written record by record from the Abaqus guide in ASCII and both binary byte orders. A 16 MB binary file from AbaqusFilFile-Translator's example (not committed) reads with its nine increments. The roadmap's check against a run's `.dat` printout still needs a real Abaqus run.

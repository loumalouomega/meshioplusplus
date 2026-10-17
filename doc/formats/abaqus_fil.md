# Abaqus results file (`.fil`)

The results file Abaqus/Standard (and Explicit, with `*FILE OUTPUT`) writes on request with `*NODE FILE` and `*EL FILE`: the one route into Abaqus results that needs neither the ODB API nor an Abaqus installation. For `.odb`, run the exporter script in Abaqus's own Python ([Abaqus `.odb` route](../routes/abaqus_odb.md), v16.13.0). The input deck itself is the [`abaqus`](./abaqus.md) format.

| | |
|---|---|
| **Format name** | `abaqus_fil` |
| **Extensions** | `.fil` (also recognised by content, ASCII and binary) |
| **Read / Write** | ✓ / ✓ ([writing](#writing); the remaining gaps are in the [roadmap](../roadmap.md#_3-reader-and-writer-parity), §3.2.6) |
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
| 1980 | a mode of an eigenvalue step: it opens a step of its own (see below) |
| 1999 | total energies: `field_data["abaqus:ALLKE"]`, `ALLSE`, `ALLWK`, … in the increment (Explicit's layout when the procedure is Explicit) |
| 301–310 | a mode-based dynamic step's generalized displacements, velocities, accelerations, base motions, phases and energies per mode: `field_data["abaqus:GU"]`, `GV`, `GA`, `BM` (the motion type and six components; the base name is dropped), `GPU`, `GPV`, `GPA`, `SNE`, `KE`, `T` |
| 1501 / 1502 | a contact surface and its facets: a `side` region (facets of rigid surfaces and unknown elements are left out, with a warning) |
| 1503 / 1504 / 1511–1592 | contact output: the node header, then `point_data["CSTRESS"]`, `CDSTRESS`, `CDISP`, `CFN`, `CFS`, `CAREA`, `CMN`, … at the slave nodes |
| 1001 / 1002 / 1011 / 1012 / 1021 / 1022 / 1031 | element matrices (see below) |
| substructure, contour integrals, cavity radiation, section output | skipped with a warning |

Element types map through the `.inp` table and, for the ones it does not list, by family and node count (`C3D…`, `CPE…`/`CPS…`/`CAX…`, shells `S…`, membranes `M3D…`, trusses and beams); their nodes are in meshio++'s order. Types with no meshio++ cell (user elements, connectors, springs) are skipped with a warning.

## Steps and data

Every increment is a step: `time_step` selects one (0 = first, negative counts from the end), and `field_data` carries `meshio:time` (the total time), `abaqus:step`, `abaqus:increment`, `abaqus:step_time` and `abaqus:procedure`. Since v16.12.0 an eigenfrequency or buckling step, which writes one increment of several eigenvectors, is **one step per mode**: each starts at its 1980 record and carries `abaqus:mode`, `abaqus:eigenvalue` and, for a frequency step, `abaqus:generalized_mass`, `abaqus:composite_damping`, `abaqus:participation_factor` and `abaqus:effective_mass` (one value per degree of freedom); the modes share their increment's total time. The increment's records become data named by the Abaqus identifier (`U`, `V`, `A`, `RF`, `CF`, `COORD`, `NT`, … for nodes; `S`, `E`, `SINV`, `PE`, `LE`, `PEEQ`, `MISES`, `ENER`, `NFORC`, … for elements; `key_<n>` for a key without one):

| Where | Becomes |
|---|---|
| nodal record | `point_data`, NaN for nodes without a value |
| element record at the integration points (location 0) | `cell_data` of shape `(cells, points × components)`, point-major |
| at the centroid (1) or for the whole element (5) | `cell_data` of shape `(cells, components)` |
| at the element nodes (2) | `cell_data` of shape `(cells, nodes × components)`, point-major in the cell's node order |
| averaged at the nodes (4) | `point_data` |
| rebar (3) | like integration-point data, named `<name>@rebar:<rebar>` (v16.12.0) |

Every array is rectangular and has the same width in every cell block, so every writer can hold it (VTU, XDMF, …): the widest block sets the number of points and components, and the other blocks are padded with NaN. Per-point data is flattened point-major, column `point × components + component`; its `(points, components)` is `field_data["abaqus:layout:<name>"]`, so `a.reshape(len(a), *layout)` restores it. A single column drops its axis. Section points above 1 (shell and beam layers; continuum elements write 0) get `@sp<k>` appended to the name: `S` holds the first layer and the solids, `S@sp5` the fifth layer. A block with no value for a field holds NaN.

Keys that differ between Abaqus/Standard and Abaqus/Explicit follow the increment's procedure (17, 21 and 74 are Explicit): key 79 is `RATIO` or `ERV`, and the total energies follow Explicit's attribute list. Explicit-only element records (`CKE`, `CKLE`, `CKLS`, `CKSTAT`, `CKEMAG`, `EMSF`, `EDT`, `CFAILST`, `CDMG`, `CDIF`, `CDIM`, `CDIP`) have their names.

**Element matrices** (`*ELEMENT MATRIX OUTPUT`, v16.12.0) are field data, read from the whole file: the values of every element one after the other in `abaqus:stiffness`, `abaqus:mass` and `abaqus:load` (a record split into continuation records is joined), and per element a row `(element label, offset, count, flag)` in `abaqus:<kind>:index`, the flag being 1 for a symmetric matrix (its upper triangle by columns, as the file stores it), 0 for a nonsymmetric one (full, by columns) and the load case for a load vector. The degrees of freedom of each element (record 1002) are `abaqus:matrix_dofs` with the same index.

Components keep the file's order, as every reader does ([mesh data model](../mesh_data_model.md)): a solid's `S` is 11, 22, 33, 12, **13, 23**, where meshio++'s six-component convention is `xx yy zz xy yz zx`. Von Mises does not depend on the order of the shear components; for principal values, swap the last two first. Von Mises is also `SINV`'s first component when the run wrote invariants; otherwise compute it from `S` with the [tensor-invariants operation](../tensor_invariants.md) on a `(cells, 6)` array, such as the mean over the points.

## Writing

```python
import meshioplusplus

meshioplusplus.write("job.fil", mesh)                       # ASCII, one increment
meshioplusplus.write("job.fil", mesh, binary=True)          # 512-word blocks, little-endian, 8-byte words

with meshioplusplus.abaqus_fil.SeriesWriter("job.fil") as w:   # the model once, an increment per step
    for time, mesh in steps:
        w.write(time, mesh)
meshioplusplus.write_sequence("job.fil", steps, binary=True)   # the same
```

`meshioplusplus.abaqus_fil.write(filename, mesh, binary=False)`. Both engines write the same bytes; the C++ core writes the file and the Python twin takes over when it declines. The writer is the reader's inverse:

- **The model** is written once: a 1921 header, the active degrees of freedom (1902), the nodes (1901, labels from `point_data["abaqus:id"]`, else 1 to n), the elements (1900 and its 1990 continuations, labels from `cell_data["abaqus:id"]`, else numbered block by block) and the `point` and `cell` regions as node sets (1931/1932) and element sets (1933/1934). A region name longer than 8 characters goes through a 1940 label cross-reference, as Abaqus does. Labels must be positive and unique, and a region or cell that refers to a point the mesh does not have is a `WriteError`.
- **Cell types** become Abaqus element types the results file can name: `line` T3D2, `line3` T3D3, `triangle` S3, `quad` S4R, `triangle6` STRI65, `quad8` S8R, `quad9` S9R5, `tetra` C3D4, `tetra10` C3D10, `pyramid` C3D5, `pyramid13` C3D13, `wedge` C3D6, `wedge15` C3D15, `hexahedron` C3D8R and `hexahedron20` C3D20R. A block with no Abaqus type is left out with a warning.
- **An increment** is a 2000 record (total time `meshio:time`, else the step's time; `abaqus:step_time`, `abaqus:procedure`, `abaqus:step` and `abaqus:increment` when present), its results and a 2001 record. Per-point arrays named by an Abaqus identifier (`U`, `V`, `A`, `RF`, `NT`, …, or `key_<n>`) are nodal records; per-cell arrays (`S`, `E`, `PEEQ`, …) are element records, at the integration points (width `points × components`, from `abaqus:layout:<name>`), at section points (`@sp<k>`), at the element nodes, at the centroid or for the whole element, and at rebar (`@rebar:<name>`): the inverse of the table under [Steps and data](#steps-and-data). A trailing NaN is trimmed from a record and an inner NaN written as 0, so a node or element with no value has no record. A value that is not finite elsewhere, or outside the range of a `D22.15` item, is a `WriteError`.
- **Modes**: `abaqus:mode`, `abaqus:eigenvalue` and the frequency step's mass and participation fields become a 1980 record before the mode's nodal records.
- **Total energies** (record 1999) are the `abaqus:ALLKE`, `ALLSE`, `ALLWK`, `ALLPD`, `ALLIE`, `ETOTAL` ... field data of the step, positional by its procedure (Explicit's names for procedures 17, 21 and 74): those present up to the last one, an absent one in between written as 0.
- **Element matrices** (1001, 1002, 1011/1012, 1021/1022, 1031): per row of each `abaqus:matrix_dofs`, `stiffness`, `mass` and `load` array's `:index` table (element, offset, length, flag) a 1001 record and the values; the flag picks the symmetric or the full record.
- **Contact surfaces**: a `side` region is a 1501 surface with a 1502 record per (element, face number), the face number from the cell's Abaqus face order.
- **Not written**, and named in one warning: arrays and regions with no place in the file, contact output (1503 to 1592) and the mode-based dynamics records (301 to 310).
- **ASCII** puts each record on 80-column lines starting with `*` (`I`, `D` with a `D` exponent, `A`, a blank line after 2001); **binary** writes 512-word blocks with 4-byte Fortran markers and pads the 2001 record to the block end.

A series writes the model from the first step and an increment per step; a step whose points or cells differ from the first step's is a `WriteError` (`{step}` in the path writes one file per step instead). The step's time is the argument of `write` (or the pair's first item in `write_sequence`). One step is `write`'s bytes.

pybaqus 0.2.17 reads the ASCII files (it does not read binary `.fil` files); this reader reads both layouts back. They were not opened in Abaqus/CAE.

**Byte order and markers.** `big_endian=True` writes the binary layout big-endian and `marker_size=8` 8-byte record markers (both go through `write`, `SeriesWriter`, `write_sequence` and the C++ `AbaqusFilWriteOptions`); the readers read every combination. Not yet: contact output, the mode-based dynamics records and 4-byte words (the rest of [roadmap 3.2.6](../roadmap.md#_3-reader-and-writer-parity)).

## Validation

No free solver writes `.fil` and no Abaqus licence was available. The fixtures under `tests/python/meshes/abaqus_fil/` are real Abaqus 2023 ASCII output from pybaqus's test suite (MIT), whose `U` and integration-point `S` match pybaqus's own reading on all 236 values, and synthetic files written record by record from the Abaqus guide in ASCII and both binary byte orders. A 16 MB binary file from AbaqusFilFile-Translator's example (not committed) reads with its nine increments.

v16.12.0: a real binary run with its printout, TenBarArea (Abaqus 6.14, T2D2 trusses, MIT), reads with the `U` of every node and `S11` of every element that its `.dat` prints (to the printed digits), and a real buckling run (Abaqus 6.23, S4R, MIT) with its five modes and eigenvalues; both are committed. The records no public file carries (frequency modes with their masses, modal dynamics, energies, contact surfaces and output, element matrices, rebar, Explicit's energies) are checked on the synthetic `extras.fil`, written from the guide. A run with shells or `C3D20R` checked against its `.dat` still needs a licence.

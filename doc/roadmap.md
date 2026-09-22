# meshio++ roadmap

Status at time of writing: **v15.6.0** — 53 core formats plus four Python-only physics-ML ones, thirty-nine mesh operations + six data operations, six language surfaces (Python / C / Fortran / Julia / R / WASM), two viewers plus a browser dataset manager, a Blender add-on, a ParaView plugin, an MCP server, a settings-driven pipeline engine, a dataset-manifest layer with a PhysicsNeMo adapter, and a versioned ABI (`MESHIOPLUSPLUS_ABI_VERSION` 15).

This document lists what is *not* built. Nothing here duplicates shipped functionality; where a feature partially exists, the shipped half is named and the gap is stated explicitly. Release history lives in [`CHANGELOG.md`](https://github.com/loumalouomega/meshioplusplus/blob/main/CHANGELOG.md), not here.

Effort key: **S** = days, **M** = a couple of weeks, **L** = a month or more, **XL** = a project in its own right. An item names a **probe** — the failing test that proves the gap — wherever one is cheap.

## How this file works

- **Sections are ordered, and the order is the recommendation.** Each section states what belongs in it; an item that sounds exciting does not move up for that reason. An empty section is removed, not kept as a placeholder.
- **A closed item is removed, not struck through.** Its history is the `CHANGELOG.md` entry and the feature's own `doc/` page; a partly closed item is narrowed to what remains (the `AGENTS.md` change checklist rule).
- **Defect-shaped items go in a *Correctness debts* section at the top, regardless of size** — behaviour that loses data, mis-orients cells, does not terminate on valid input or fails silently is not a feature request, even when the fix and the feature are the same work. None is open today, so there is no such section; the first defect found opens one as §1 and renumbers the rest.
- **An item estimated from a doc, a `.d.ts` or a changelog alone says so** ("verify first") and names its probe; the code has repeatedly been more or less capable than its description.
- **[Non-goals](#non-goals-and-decisions-taken) record decisions already taken**, with their reasons, so they are not re-proposed as gaps.

## The map

![The roadmap at a glance: open items grouped by section, shaded by effort, with dependency arrows and the items that need a design pass or a research spike first](/diagrams/roadmap_map.svg)

---

## 1. Format reach

*Admission: a format a simulation or physics-ML workflow actually exchanges, or the missing half of a shipped one. "Exchanges" means a file that crosses a tool boundary in a real pipeline (mesher → solver, solver → post-processor, solver → training set), not a format that merely exists. Every item names the consumer on the other side of the file.*

Sizes: **S** ≈ days, **M** ≈ two weeks, **L** ≈ a month or more, including tests, docs page and CLI wiring. Sizes assume the shared infrastructure in §1.21 exists; the first format that needs a component pays for it.

Link legend: unmarked links were opened or returned by a search while this section was written (September 2026); links marked † are quoted from memory and must be checked before relying on them. Vendor documentation portals move often — when a link is dead, search the document title.

---

### A. FEM interchange, solver inputs and result files — native implementations

Python meshio covers none of §1.1–§1.18 except OptiStruct `.fem` as plain Nastran bulk; these are greenfield. Cross-validation partners are named in each block.

### 1.1 MSC Nastran HDF5 results `.h5` (read) — **S–M**

Result database written by MSC Nastran (since 2016, `MDLPRM,HDF5`); on the other side is Patran or a post script. Self-describing; MSC ships the schema with each release.

- **Why.** HDF5 is already linked; far cheaper than OP2 for the same results.
- **Scope (read).** Mesh from `/NASTRAN/INPUT`, nodal and element results from `/NASTRAN/RESULT`, subcases/modes/frequencies as a sequence.
- **Structure.** Compound-type tables: `/NASTRAN/INPUT/NODE/GRID`, `/NASTRAN/INPUT/ELEMENT/*`, `/NASTRAN/RESULT/NODAL/DISPLACEMENT` (and `_CPLX`), `/NASTRAN/RESULT/ELEMENTAL/…`, `/NASTRAN/RESULT/DOMAINS` (subcase, step, mode, frequency…), `/INDEX/…` giving row ranges per domain.
- **Mapping.** GRID → points; ELEMENT tables → cells, regions by property id; NODAL → point data; ELEMENTAL → cell data; DOMAINS rows → sequence entries; `_CPLX` → real/imaginary pair.
- **Pitfalls.** Result rows must be joined to DOMAINS through `DOMAIN_ID`/INDEX; other Nastran vendors write *incompatible* HDF5 schemas — detect and refuse rather than misread; schema grows every release.
- **Out of scope.** OP2/OP4/punch (see §1.17); non-MSC HDF5 dialects until a file is in hand.
- **Done when.** A SOL 103 `.h5` converts to a `.vtu` sequence with eigenvectors per mode matching pyNastran.
- **References.** [The Nastran HDF5 result database (Hexagon docs, login may be required)](https://help.hexagonmi.com/bundle/MSC_Nastran_2021/page/Nastran_Combined_Book/refman/interother/TOC.The.Nastran.HDF5.Result.xhtml) · [MSC blog: HDF5 in Nastran and Patran](https://simulatemore.mscsoftware.com/hdf5-a-useful-enhancement-for-msc-nastran-and-patran/) · [pyNastran issue with access pattern](https://github.com/SteveDoyle2/pyNastran/issues/308) · [NASTRAN_CoFE: open-source solver that *writes* MSC-format `.h5` — a licence-free test-file generator](https://github.com/vtpasquale/NASTRAN_CoFE) · [vcti-dataspace-nastran-h5](https://pypi.org/project/vcti-dataspace-nastran-h5/)

### 1.2 Code_Aster `.mail` mesh (read / write) — **S**

Native Code_Aster ASCII mesh; on the other side is Code_Aster/Salome-Meca.

- **Why.** Cheap: node ordering is MED's, and the MED writer already ships.
- **Structure.** Free-format, 80-column lines, blank and comma separators. `TITRE`; `COOR_2D`/`COOR_3D`; element keywords `POI1`, `SEG2/3/4`, `TRIA3/6/7`, `QUAD4/8/9`, `TETRA4/10`, `PENTA6/15/18`, `PYRAM5/13`, `HEXA8/20/27`; `GROUP_NO`/`GROUP_MA`; each block ends at `FINSF`, the file at `FIN`. Entities are *named* (`N1`, `M1`), not numbered.
- **Mapping.** Keywords → cell types through the MED permutations; `GROUP_MA` → regions/cell sets; `GROUP_NO` → node sets.
- **Pitfalls.** Names are limited to 8 characters in the classic format (MED allows far longer group names — truncate with collision checks on write); a keyword stays active until `FINSF`.
- **Out of scope.** `.comm` command files; `.rmed` results (already MED).
- **Done when.** A `HEXA20` `.mail` round-trips through MED and back with identical connectivity and groups.
- **References.** [U3.01.00, mesh file description (EN, v14)](https://code-aster.org/V2/doc/v14/en/man_u/u3/u3.01.00.pdf) · [same, v12](https://code-aster.org/doc/v12/en/man_u/u3/u3.01.00.pdf) · [U4.21.01 `LIRE_MAILLAGE`](https://code-aster.org/V2/doc/v12/fr/man_u/u4/u4.21.01.pdf) · [U7.01.21 MED reading, group-name limits](https://code-aster.org/V2/doc/v12/fr/man_u/u7/u7.01.21.pdf) · [U1.03.00 principles (EN)](https://biba1632.gitlab.io/code-aster-manuals/docs/user/u1.03.00.html)

### 1.3 Altair OptiStruct `.fem` (read; Nastran reader extension) — **S**

HyperMesh/OptiStruct bulk-data dialect; on the other side is OptiStruct.

- **Why.** A test matrix and a small extension for the shipped Nastran reader, not a new format — but it carries component names that plain bulk data lacks.
- **Structure.** Nastran small/large/free field, an I/O-options section before `BEGIN BULK`, HyperMesh comment cards (`$HMNAME COMP`, `$HMMOVE`, `$HWCOLOR`) holding component names and membership, optimization cards (`DESVAR`, `DRESP1`, `DTPL`…), `CONTACT`/`TIE`.
- **Mapping.** `$HMNAME COMP` + `$HMMOVE` → named regions; everything else through the Nastran path.
- **Pitfalls.** The `$HM` comments are load-bearing; optimization and contact cards must be skipped with a warning, not an error.
- **Done when.** A HyperMesh-exported `.fem` reads with components as named regions and unknown optimization cards skipped.
- **References.** [OptiStruct help](https://help.altair.com/hwsolvers/os/index.htm) † · [pyNastran (has an OptiStruct mode)](https://github.com/SteveDoyle2/pyNastran)

### 1.4 COMSOL `.mphtxt` / `.mphbin` mesh (read / write) — **M**

COMSOL's native mesh exchange; on the other side is COMSOL Multiphysics.

- **Why.** Large multiphysics user base; the documented route in and out of COMSOL that keeps domain ids.
- **Structure.** Version line (`0 1`), tags, types, then serialised objects (`0 0 1`, class name). `Mesh` object (version 4): `sdim`, number of vertices, lowest vertex index, coordinates, then per element type (`vtx`, `edg`, `tri`, `quad`, `tet`, `pyr`, `prism`, `hex`, and second-order `edg2`, `tri2`, `tet2`, `hex2`…): nodes per element, element count, connectivity, geometric-entity index per element. Optional `Selection` objects (label, dimension, entity list). `.mphbin` is the same serialisation in binary.
- **Mapping.** Domain-level elements → cells, entity index → region; boundary/edge/vertex elements → face/edge/node sets; Selections → named sets.
- **Pitfalls.** Node ordering is tensor-product, **not VTK**, for quad, hex, prism and pyramid; entity indices are 0-based for points/edges/boundaries but **1-based for domains**; the lowest vertex index may be 0 or 1.
- **Out of scope.** Geometry objects; spreadsheet/sectionwise result exports.
- **Done when.** A mixed tet/prism file with two domains converts to `.vtu` with correct node order and domain ids, and a written file imports into COMSOL.
- **References.** [File structure (v6.1)](https://doc.comsol.com/6.1/doc/com.comsol.help.comsol/comsol_api_fileformats.53.18.html) · [Serializable classes](https://doc.comsol.com/6.1/doc/com.comsol.help.comsol/comsol_api_fileformats.53.23.html) · [`Mesh` class fields](https://doc.comsol.com/5.6/doc/com.comsol.help.comsol/comsol_api_fileformats.50.30.html) · [`Selection` class](https://doc.comsol.com/6.1/doc/com.comsol.help.comsol/comsol_api_fileformats.53.39.html) · [Mesh Import and Export Guide (PDF)](https://www.comsol.com/model/download/1273541/COMSOL_MeshImportExportGuide.pdf) · [FEconv](https://github.com/victorsndvg/FEconv) † and [ElmerGrid](https://www.nic.funet.fi/index/elmer/doc/ElmerGridManual.pdf) as readers to compare against

### 1.5 FEBio `.feb` input and `.xplt` results (read / write `.feb`; read `.xplt`) — **M**

Biomechanics FEM; on the other side is FEBio / FEBio Studio.

- **Structure.** `.feb`: XML, `febio_spec` 2.5/3.0/4.0 — `<Mesh>` with `<Nodes>`, `<Elements type= name=>`, `<NodeSet>`, `<Surface>`, `<ElementSet>`, `<DiscreteSet>`; `<MeshDomains>`; `<MeshData>`. `.xplt`: tagged chunk binary (root tag `0x00464542`), header, dictionary of nodal/domain/surface variables, mesh section, one state section per step; optional compression.
- **Mapping.** `<Elements>` blocks → regions; sets and surfaces → node/face/element sets; `.xplt` dictionary → point/cell data; states → sequence.
- **Pitfalls.** Tag layout changes between spec 2.5, 3.0 and 4.0 — dispatch on the version attribute; the published binary spec lags the code (version tag `0x0004`+ in practice, node ids added in 3.3), so FEBio Studio's reader is the working reference.
- **Out of scope.** Materials, steps, loads.
- **Done when.** A `.feb` mesh round-trips and its `.xplt` converts to a `.vtu` sequence matching FEBio Studio.
- **References.** [FEBio manual: plotfile section](https://help.febio.org/FEBio/FEBio_um_2_9/FEBio_um_2-9-Subsection-3.17.2.html) · [Binary database spec discussion, version 0x0004](https://forums.febio.org/forum/febio-forums/users-forum/1172-xplt-binary-file-spec-log-files) · [Spec PDF location](https://forums.febio.org/forum/febio-forums/users-forum/1102-outputting-custom-variables) · [FEBio Studio (reference `xplt` reader, MIT)](https://github.com/febiosoftware/FEBioStudio) · [FEBio](https://github.com/febiosoftware/FEBio/wiki) · [febio-python](https://github.com/Nobregaigor/febio-python) · [interFEBio](https://github.com/andresutrera/interFEBio)

### 1.6 Elmer mesh directory (read / write) — **S–M**

ElmerSolver's native mesh; on the other side is Elmer.

- **Structure.** A *directory*: `mesh.header`, `mesh.nodes`, `mesh.elements`, `mesh.boundary`, optional `mesh.names`; partitioned meshes under `partitioning.N/part.n.{header,nodes,elements,boundary,shared}`. Element codes = family × 100 + node count: 101; 202/203; 303/306; 404/408/409; 504/510; 605/613; 706/715; 808/820/827. Boundary lines: id, boundary number, two parent elements, type, nodes.
- **Mapping.** Body ids → regions; boundary numbers → face sets; `mesh.names` → names; partitions ↔ `partition`.
- **Pitfalls.** The path is a directory (shared API work with OpenFOAM); parent-element references must be regenerated on write; known third-party export bugs with quadratic tets make a good regression test.
- **Out of scope.** `.sif` files; results (Elmer already writes VTU).
- **Done when.** A quadratic-tet mesh directory round-trips and ElmerSolver runs on the written copy.
- **References.** [ElmerSolver manual, appendix on mesh file format](https://www.nic.funet.fi/index/elmer/doc/ElmerSolverManual.pdf) · [ElmerGrid manual](https://www.nic.funet.fi/index/elmer/doc/ElmerGridManual.pdf) · [Forum: `mesh.boundary` line layout](https://www.elmerfem.org/forum/viewtopic.php?t=7761) · [elmerfem source](https://github.com/ElmerCSC/elmerfem) †

### 1.7 MSC Patran neutral file `.pat` / `.out` (read / write) — **S–M**

Legacy neutral file still exported by Patran, Cubit and ANSA and imported by Fluent and Feko.

- **Structure.** "Packets": a header card `I2,8I8` (`IT, ID, IV, KC, N1–N5`) followed by `KC` data cards. 25 title, 26 summary, 01 node, 02 element (shape in `IV`: 2 bar, 3 tri, 4 quad, 5 tet, 7 wedge, 8 hex), 03 material, 04 element properties, 06–08 loads and displacements, 21 named component, 99 end.
- **Mapping.** Packet 01 → points; 02 → cells; 21 → regions/sets; property id → region fallback.
- **Pitfalls.** Fixed-width cards; element shape lives in the header, not the data; ids have gaps; only the text form is supported by Patran itself.
- **Out of scope.** Loads/BC packets; result files `.nod`/`.els`/`.dis` until a consumer asks.
- **Done when.** A Cubit-exported hex/tet file round-trips with named components as sets.
- **References.** [Interface to PATRAN 2 Neutral File guide (PDF)](https://help-be.hexagonmi.com/bundle/Patran_2021.1_Interface_to_PATRAN_2_Neutral_File_Preference_Guide/raw/resource/enus/Patran_2021.1_Interface_to_PATRAN_2_Neutral_File_Preference_Guide.pdf) · [Patran reference: the neutral file, packet layout](http://web.mscsoftware.com/training_videos/patran/reverb3/Basic%20Functions/formats_topics.13.3.html) · [Fluent's supported packets](https://www.afs.enea.it/project/neptunius/docs/fluent/html/ug/node176.htm) · [Feko's supported packets](https://help.altair.com/2023/feko/topics/feko/user_guide/appendix/editfeko_cards/geometry/card_in_import_patran_neutral_file_feko_r.htm)

### 1.8 Femap neutral file `.neu` (read; mesh write) — **M**

Femap's documented interchange and the only open route to a Femap model.

- **Structure.** Data blocks bracketed by `-1` lines, block id on the next line: 100 header (version), 402 properties, 403 nodes, 404 elements, 405 coordinate systems, 408 groups, 413 layers, 450 output sets, 451 output vectors (1051 in newer versions), 601 materials. Element topology ids (0 line2, 2 tri3, 3 tri6, 4 quad4, 5 quad8, 6 tet4, 7 wedge6, 8 brick8, 10 tet10, 11 wedge15, 12 brick20) with zero-padded 20-slot node arrays.
- **Mapping.** 403/404 → mesh; property id → region; 408 → sets; 450/451 → point/cell data and a sequence.
- **Pitfalls.** Block layouts depend on the Femap version in block 100 — dispatch on it and keep one sample per version.
- **Out of scope.** The binary `.modfem` database (§1.20).
- **Done when.** A brick20 model reads identically from files written by two Femap versions.
- **References.** [Femap neutral file format (document mirror)](https://vdocuments.mx/femap-neutral-file-format.html) · the authoritative copy is `neutral.pdf` in the Femap installation's `pdf/` folder

### 1.9 MFEM `.mesh` and `.gf` (read / write) — **M**

MFEM/GLVis native mesh and grid functions; on the other side is a high-order research code.

- **Structure.** `MFEM mesh v1.0`–`v1.3`, `MFEM NC mesh v1.0`: `dimension`; `elements` (attribute, geometry type 0 point … 7 pyramid, vertices); `boundary`; `vertices`; optional `nodes` block (a `FiniteElementSpace` + `GridFunction`: collection, order, `VDim`, `Ordering`) for curved meshes; v1.3 attribute sets; parallel format is per-rank with `communication_groups`.
- **Mapping.** Attribute → region; boundary attribute → face set; `nodes`/`.gf` of H1 order p → Lagrange cells of order p, else corner-only with a warning.
- **Pitfalls.** High-order nodal ordering vs VTK Lagrange ordering; `Ordering` byNODES vs byVDIM; NC meshes.
- **Out of scope.** NC hierarchy; non-H1 spaces; data collections.
- **Done when.** An order-2 mesh with a `.gf` field converts to an order-2 Lagrange `.vtu` matching GLVis.
- **References.** [Mesh formats](https://mfem.org/mesh-formats/) † · [Mesh v1.0 format](https://mfem.org/mesh-format-v1.0/) · [v1.x format discussion](https://github.com/mfem/mfem/issues/205) · [MFEM source and `data/` samples](https://github.com/mfem/mfem) †

### 1.10 libMesh `.xda` / `.xdr` (read) — **S–M**

libMesh/MOOSE-adjacent native mesh. `.xda` is ASCII, `.xdr` the same stream in big-endian XDR. Header: version string, element and node counts, boundary-condition counts, field sizes; element blocks by refinement level; subdomain ids; side/node-set ids and names.

- **Pitfalls.** XDR endianness; AMR level and p-level fields; `TET14`, `PYRAMID14/18`, `PRISM20/21` have no VTK counterpart (drop extra nodes with a warning).
- **Done when.** A `HEX27` file reads with correct ordering and subdomains as regions. MOOSE users are already served by Exodus.
- **References.** [libMesh source (`XdrIO`)](https://github.com/libMesh/libmesh) † · [libMesh docs](https://libmesh.github.io/) †

### 1.11 Z88 / Z88Aurora (read / write mesh; read results) — **S**

Open-source teaching/SME FEM. `z88i1.txt`: line 1 = dimension, nodes, elements, DOF, flags; node lines; two-line element records with Z88 element numbers (e.g. 1 hex8, 10 hex20, 16 tet10, 17 tet4). `z88o2.txt` displacements, `z88o3.txt` stresses.

- **Done when.** A hex20 structure file round-trips and displacements attach as point data.
- **References.** [Z88 home and manuals](https://z88.de/) †

### 1.12 Abaqus `.fil` results file (read) — **M**

The only route into Abaqus results that needs neither the ODB API nor an Abaqus install; written on request (`*NODE FILE`, `*EL FILE`, …).

- **Structure.** Sequential records of 8-byte words: `[length, key, attributes…]`. Binary by default; ASCII (`*FILE FORMAT, ASCII`) writes `*` record starts with `I`/`D`/`A`-tagged items on 80-column lines. Keys: 1900 element definition (+1990 continuation), 1901 node, 1911 output request (0 element / 1 nodal / 2 modal / 3 set energy; set name), 1921 release and counts, 1922 heading, 1931/1933 sets, 1940 label cross-reference, 2000 increment start, 2001 increment end; data keys such as 101 `U`, 104 `RF`, 11 `S`, 21 `E`; record 1 = element header (element, integration point, section point, location).
- **Mapping.** 1900/1901 → mesh; sets → sets; nodal records → point data; element records → cell data per integration point; 2000/2001 → sequence.
- **Pitfalls.** Labels over 8 characters become integers resolved through record 1940; keys differ between Standard and Explicit; Explicit writes degenerate bricks as wedges/tets; binary files carry Fortran record markers; location (integration point / centroid / nodal average) comes from the element header.
- **Out of scope.** Writing; substructure records.
- **Done when.** A Standard static case converts to `.vtu` with `U` and von Mises from `S` matching the `.dat` printout.
- **References.** [Results file output format, record keys (2016 guide mirror)](https://ceae-server.colorado.edu/v2016/books/usb/pt02ch05s01afi01.html) · [Records written for any file output request (2017 mirror)](https://abaqus-docs.mit.edu/2017/English/SIMACAEOUTRefMap/simaout-c-recordswrittenforanyfileoutputrequest.htm) · [ASCII item encoding](http://abaqusdocs.eait.uq.edu.au/v6.10ef/books/usb/pt02ch04s01aus36.html) · [User post-processing of results files, overview](https://ceae-server.colorado.edu/v2016/books/exa/ch15s01abo02.html) · [FJOIN example: which header records to keep](https://classes.engineering.wustl.edu/2009/spring/mase5513/abaqus/docs/v6.6/books/exa/ch12s01aex131.html)

### 1.13 Radioss / OpenRadioss Starter `_0000.rad` (read) — **S–M**

A fourth keyword family; on the other side is an open explicit solver. OpenRadioss reads LS-DYNA `.k` natively, and `.k` shipped in v15.2.0, so check whether that already covers the demand; if not, build on the card tokenizer in `detail/keyword_card.hpp`.

- **Structure.** `#RADIOSS STARTER`, `/BEGIN` (version, units), then `/KEYWORD/option/id` blocks with 10-column fields: `/NODE`, `/BRICK`, `/TETRA4`, `/TETRA10`, `/BRIC20`, `/SHELL`, `/SH3N`, `/QUAD`, `/BEAM`, `/TRUSS`, `/SPRING`, `/PART`, `/SUBSET`, `/GRNOD`, `/GRBRIC`, `/GRSHEL`, `/SURF`, `#include`.
- **Mapping.** `/PART` → regions (elements are grouped under their part id); `/GR*` and `/SURF` → sets; `/SUBSET` → region hierarchy (flattened).
- **Results.** `A001…` animation files are converted by OpenRadioss's own `anim_to_vtk` (and community anim → d3plot / VTKHDF converters); document that route rather than parsing ANIM natively. `T01` time history → `th_to_csv`.
- **Done when.** An OpenRadioss example deck reads with parts as regions and correct element totals.
- **References.** [Radioss reference guide](https://help.altair.com/hwsolvers/rad/index.htm) † · [OpenRadioss](https://github.com/OpenRadioss/OpenRadioss) · [`tools/anim_to_vtk`, `th_to_csv`](https://github.com/OpenRadioss/OpenRadioss/tree/main/tools) · [OpenRadioss Tools repo (converters, GUI)](https://github.com/OpenRadioss/Tools) · [ANIM structure discussion](https://github.com/orgs/OpenRadioss/discussions/552) · [Vortex-Radioss anim → d3plot](https://github.com/Vortex-CAE/Vortex-Radioss) †

### 1.14 MSC Marc input `.dat` (read) and `.t19` (read, later) — **M**

Marc's input deck, and the formatted ASCII variant of its post file.

- **Structure.** Parameter section (`title`, `sizing`, `elements`, `end`), model definition (`connectivity`, `coordinates`, `define element|node set`, …, `end option`), history definition. Fixed 5/10-column, extended-precision and comma free formats. Marc element type numbers (7 = hex8; verify 21 hex20, 134 tet4, 127 tet10, 75 shell against Volume B before coding) with Marc node orderings.
- **Mapping.** Element type → cell + permutation; `define` sets → sets/regions.
- **`.t19`.** Post codes and record layout are described in Volume C; a native reader is feasible and is the only licence-free route to Marc results. The binary `.t16` stays on the PyPost route (§1.20).
- **Done when.** A hex20 deck reads with correct node order and element sets as regions.
- **References.** [Marc 2024.2 Volume C: Program Input (PDF)](https://documentation-be.hexagon.com/bundle/Marc_2024.2-Volume_C_Program_Input/raw/resource/enus/Marc_2024.2-Volume_C_Program_Input.pdf) · [Volume C, older, with post codes (mirror)](http://www.sd.ruhr-uni-bochum.de/downloads/links/marc_manuals/online_documentation_marc_2005/volc.pdf) · [Volume A: theory and user information (mirror)](http://www.sd.ruhr-uni-bochum.de/downloads/links/marc_manuals/online_documentation_marc_2003/vola.pdf) · Volume B (element library) sits beside them on the same mirror

### 1.15 ANSYS `.cdb` archive (read / write) — **S–M**; `.rst` / `.rth` results (read) — **L**

`CDWRITE` archives and MAPDL binary results; on the other side is Ansys Mechanical/MAPDL.

- **`.cdb`.** Blocked ASCII: `NBLOCK` with a literal Fortran format line (e.g. `(3i9,6e21.13e3)`), `EBLOCK` (solid and non-solid layouts; 19 fields then continuation), `CMBLOCK` (components, negative values = ranges), `ET`/`ETBLOCK`, `TYPE`/`MAT`/`REAL`. Element *type number* → ANSYS element name → topology; degenerate shapes by repeated nodes.
- **`.rst`.** Fortran-style blocked binary: 100-integer standard header, results header, nodal and element equivalence tables, data-set index, per-set solution header with *relative* pointers, nodal solution, element solution records (`ENS`, `EEL`, …). Compressed/sparse records exist (`/FCOMP,RST,0` avoids them).
- **Mapping.** Components → sets; element type/material → regions; data-set index → sequence; nodal/element solution → point/cell data.
- **Pitfalls.** Parse the format line, never assume widths; 32- vs 64-bit pointers; distributed-solve file layouts.
- **Done when.** A `.cdb` with mixed solids and components round-trips; a modal `.rst` reads mode shapes matching pymapdl-reader.
- **References.** [Coded database file commands (`NBLOCK`, `EBLOCK`, `CMBLOCK`)](https://ansyshelp.ansys.com/public/Views/Secured/corp/v251/en/ans_prog/Hlp_P_INT3_3.html) · [Guide to Interfacing with ANSYS: format of binary data files (old PDF mirror)](http://www.free-motor.org/DATA/p_int_2188.pdf) · [pymapdl-reader (MIT; C readers built from Ansys headers)](https://github.com/ansys/pymapdl-reader) · [mapdl-archive](https://github.com/akaszynski/mapdl-archive) † · [PyDPF, the vendor route](https://dpf.docs.pyansys.com/) †

### 1.16 LS-DYNA `d3plot` family (read) — **L**

The state database behind most public crash datasets; on the other side is LS-PrePost — and physics-ML pipelines that start from d3plot.

- **Structure.** Word-addressed binary: control words (title, `NDIM`, `NUMNP`, `NEL8/2/4/T`, `NV3D/2D/1D`, `NGLBV`, `IT/IU/IV/IA`, `MATTYP`, `MAXINT`, `MDLOPT`, `NARBS`, `EXTRA`…), geometry, user ids, then states (time, globals, nodal, element, deletion arrays). Families split across `d3plot`, `d3plot01`, … on 512-word boundaries. Siblings `d3part`, `d3eigv`, `d3thdt`; `binout` is a separate LSDA container.
- **Mapping.** Parts → regions; states → sequence; deletion arrays → a mask field; globals → field data.
- **Pitfalls.** Single/double precision and endianness must be sniffed; element connectivity is degenerate like the keyword file (the collapse rules are in `formats/lsdyna.cpp`); `MAXINT` encodes integration points and flags; adaptive runs restart the geometry.
- **Out of scope.** Interface force files, CPM/airbag particles, multi-solver data.
- **Done when.** A shell + solid family converts to a sequence with displacement, plastic strain and deletion matching lasso-python.
- **References.** [LS-DYNA Database Binary Output Files (2014 revision)](https://www.dynasupport.com/manuals/additional/ls-dyna-database-manual-2014) · [lasso-python docs](https://open-lasso-python.github.io/lasso-python/dyna/) · [lasso-python source (MIT)](https://github.com/open-lasso-python/lasso-python/releases) · [dynareadout (C, d3plot + binout + key files)](https://github.com/PucklaJ/dynareadout)

### 1.17 Nastran OP2 (read) — **M–L**

Fortran-record binary of named data blocks (`GEOM1–4`, `EPT`, `MPT`, `OUGV1`, `OES1X`, `OSTR1`, `OEF1X`, `OQG1`, `OGPFB1`, `ONRGY1`…), MSC and NX/Simcenter dialects, 32/64-bit. Only after §1.1 shows demand beyond HDF5.

- **Mapping.** `OUG*` → point data; `OES`/`OSTR`/`OEF` → cell data; subcases/modes → sequence.
- **Pitfalls.** The table zoo is enormous — scope to displacements, eigenvectors, solid/shell stress and strain, SPC forces; refuse the rest by name.
- **Done when.** A SOL 101 file's displacement and stress match pyNastran.
- **References.** [pyNastran (BSD; the reference implementation)](https://github.com/SteveDoyle2/pyNastran) · [pyNastran docs](https://pynastran-git.readthedocs.io/) †

### 1.18 Tecplot binary `.plt` (read) — **M**

Documented in the Data Format Guide's appendix: magic `#!TDV112`, header (title, variables, zone records marked `299.0`), end-of-header marker `357.0`, then zone data. Reuses the ASCII reader's zone → region/sequence mapping.

- **Out of scope.** `.szplt` — TecIO only (§1.20).
- **Done when.** An FE `.plt` written by `preplot` from an ASCII file converts identically to that ASCII file.
- **References.** [Data Format Guide, binary file format appendix](https://tecplot.azureedge.net/products/360/current/360-data-format.html) · [`tecio`, a pure-Python `.plt` reader/writer (GPL — behaviour reference only, do not copy code)](https://pypi.org/project/tecio/) · [Tecplot file types](https://tecplot.com/2016/09/16/tecplot-data-file-types-dat-plt-szplt/)

---

### B. Formats that need a vendor runtime or a heavy optional dependency

### 1.19 DOLFINx ADIOS2 `.bp` (read; optional ADIOS2 build) — **M**

`VTXWriter` output: an ADIOS2 BP4/BP5 directory with a `vtk.xml` schema attribute and step-wise `geometry`, `connectivity`, `types`, `NumberOfNodes`, `NumberOfEntities` plus field variables. ADIOS2 is a heavy dependency for one consumer, so this is an opt-in build flag, never a default. `adios4dolfinx` checkpoints and Fides output are different layouts — out of scope.

- **Done when.** A DOLFINx demo's `.bp` reads as a sequence matching ParaView's VTX reader.
- **References.** [DOLFINx](https://github.com/FEniCS/dolfinx) † · [ADIOS2 documentation](https://adios2.readthedocs.io/) † · [VTK ADIOS2 module (VTX schema reader)](https://docs.vtk.org/en/latest/modules/vtk-modules/IO/ADIOS2/README.html) · [VTX reader changes in VTK 9.4](https://docs.vtk.org/en/latest/release_details/9.4/adios2-vtx-reader-changes.html) · [adios4dolfinx](https://github.com/jorgensd/adios4dolfinx) †

### 1.20 Vendor-runtime routes (exporter scripts and optional plugins, no native parser) — **M** in total

For these, the deliverable is a documented, tested route — a script that runs inside the vendor's own Python and writes VTKHDF/XDMF, or a plugin compiled against an SDK the user already owns. meshio++ never redistributes vendor libraries.

- **Abaqus `.odb`.** No public specification; readable only through the ODB API on a licensed install, version-locked (`abaqus upgrade -odb`). Ship an `abaqus python` exporter script (instances → regions, steps/frames → sequence, field outputs by position) and document it; an optional C++ plugin against the ODB API comes only if a user asks. Where no licence is available, §1.12 (`.fil`) is the fallback. — [ODB2VTK (reference exporter)](https://github.com/Arris-Composites/ODB2VTK) · [abqpy](https://github.com/haiiliin/abqpy) †
- **MSC Marc `.t16`.** Read through PyPost (`py_post`), shipped with Marc/Mentat. Ship an exporter script; §1.14 covers `.t19` natively.
- **Femap `.modfem`.** COM API on Windows only. The route is "export a neutral file" → §1.8.
- **Tecplot `.szplt`.** Undocumented; TecIO only. Optional link against a user-installed TecIO, or "re-save as `.plt`" → §1.18. — [TecIO](https://tecplot.com/products/tecio-library/) †
- **ANSYS results beyond §1.15.** Document the PyDPF export route.
- **Done when.** Each route has a docs page, a script under `contrib/`, and one manual test recorded with the vendor version used.

---

### 1.21 Shared infrastructure worth building once

- **Fixed-width card tokenizer, remaining half: format-line parsing.** `detail/keyword_card.hpp` (v15.2.0, built for LS-DYNA) already splits a card into standard, long, I10 and free layouts per line and parses real fields with Fortran spellings; Nastran/OptiStruct, Radioss, Marc, Patran, UNV and ANSYS `.cdb` should adopt it instead of growing their own (the `.frd` reader of v15.3.0 slices its own columns: the layout is fixed by the record key, not by a format line). What it lacks is parsing a Fortran format string such as `(3i9,6e20.13)`, which `.cdb` needs.
- **Node-ordering permutation registry** keyed by (format, element type), with a self-test that maps a reference element through every table and checks a positive Jacobian — UNV, COMSOL, Femap, Marc, libMesh, MFEM high-order, ANSYS, MED/`.mail`. The `.frd` he20/pe15/be3 tables (v15.3.0) and the UNV sandwich tables (v15.6.0, pinned against gmsh and Salome) are constants inside `formats/frd.cpp`/`formats/unv.cpp` and their Python twins and should move into it.
- **Fortran unformatted-record reader** with endianness and 4/8-byte marker sniffing — OP2, ANSYS `.rst`, Abaqus `.fil`, EnSight Fortran-binary. (`d3plot` is word-addressed, not record-framed, but shares the sniffing.)
- **HDF5 schema helper, remaining half: typed compound-table reads** for Nastran `.h5` (§1.1). The rest of it — fixed-length string attributes, integer-array attributes, hyperslab row reads, chunked appendable datasets, creation-ordered groups and soft links — shipped in `detail/hdf5_util.hpp` with VTKHDF.
- **Directory-as-format support** in the path and CLI layer — Elmer, OpenFOAM, ADIOS2 `.bp`.
- **Optional-plugin mechanism and a `contrib/` script convention** — §1.19, §1.20.

### 1.22 Suggested order

1. **First FEM wave (build tokenizer + permutation registry):** §1.2 `.mail` → §1.3 OptiStruct.
2. **HDF5 results:** §1.1 Nastran `.h5`.
3. **Second FEM wave:** §1.4 COMSOL → §1.6 Elmer → §1.5 FEBio → §1.15 ANSYS `.cdb`.
4. **Legacy and research reach:** §1.7 Patran → §1.8 Femap → §1.12 Abaqus `.fil` → §1.9 MFEM → §1.10 libMesh → §1.11 Z88 → §1.13 Radioss → §1.14 Marc.
5. **Heavy binaries, on demand:** §1.16 d3plot (pull forward if crash-dataset users appear) → §1.15 `.rst` → §1.17 OP2 → §1.18 `.plt`.
6. **Routes:** §1.20 scripts as users ask; §1.19 when a FEniCSx user asks.

### Considered, not queued

- **Graphics, CAD, VFX or lab-specific rather than simulation interchange:** 3MF, X3D, DXF, Alembic, OpenVDB, Silo, LAS/LAZ, E57, OBJ materials (`.mtl`).
- **Closed FEM databases with no route worth shipping:** Altair `.h3d` (undisclosed format, SDK not redistributable — Altair's own advice for OpenRadioss is ANIM → VTK, see the [OpenRadioss discussion](https://github.com/orgs/OpenRadioss/discussions/552)); HyperMesh `.hm`; Abaqus `.sim`.
- **No native format to implement:** deal.II (reads and writes UCD, VTK, Gmsh, Abaqus, Exodus — all reachable today); MOOSE (Exodus).
- **Nastran OP4 matrices and punch `.pch`:** matrix and card dumps, not mesh/field exchange; pyNastran serves them.
- **Non-MSC Nastran HDF5 dialects:** incompatible schemas; revisit with a file in hand.
- **Partitioning niche:** gmsh `$PartitionedEntities` (Kratos partitions through MDPA).
- **Operators and tooling, not formats** (listed only because they came up in the same review): ICP registration; mesh duals; cylindrical/spherical coordinate transforms; shell completions.

Revisit any of them when a consumer asks with a file in hand — a real deck or result file, the tool version that wrote it, and the tool that needs the converted output.

### Open verification items before coding

- Marc element numbers other than 7, and their node orderings, against Volume B (§1.14).
- Every link marked †.

---

## 2. Spack package upkeep

*Admission: a package that already exists but lives in someone else's repository and has fallen behind the release. Not a feature, and ordered straight after format reach because it is how HPC users get meshio++ at all ([issue #3](https://github.com/loumalouomega/meshioplusplus/issues/3): "Interesting for HPC"; the maintainer's follow-up: done, but it has to be kept current with every release).*

Both recipes are upstream in [`spack/spack-packages`](https://github.com/spack/spack-packages) — [`meshioplusplus`](https://github.com/spack/spack-packages/blob/develop/repos/spack_repo/builtin/packages/meshioplusplus/package.py) (a `CMakePackage`: C API, Fortran, installable C++ API, CLI) and [`py-meshioplusplus`](https://github.com/spack/spack-packages/blob/develop/repos/spack_repo/builtin/packages/py_meshioplusplus/package.py) (a `PythonPackage`), added in [PR #5624](https://github.com/spack/spack-packages/pull/5624). Nothing in this repository builds or tests them, so nothing notices when they drift.

- **Verify first, then bump the versions.** Read on 2026-09-21, the newest tagged `version()` in both recipes is **9.10.0** (plus `master`) against **15.2.0** here, and the `meshioplusplus` recipe's `url` still points at the `v9.10.0` tarball. Probe: `spack versions meshioplusplus py-meshioplusplus` against `git tag`. Add a `version(..., sha256=...)` line per release worth keeping (`spack checksum` computes it); `spack install py-meshioplusplus@15.2.0` and `meshioplusplus@15.2.0` must then concretize and build. **S**
- **Audit the recipes against six major versions of drift.** Every `when="@X:"` guard, the `depends_on` floors (`py-scikit-build-core@0.8:`, `py-pybind11@2.11:`, `python@3.8:`, `cmake@3.15:`) and the `conflicts("%gcc@:9")` must be re-checked against `pyproject.toml` and `CMakeLists.txt` at the new version, and the docstring's "~40 unstructured mesh formats" is now 52. The `+cxx_api` libraries install as `libmeshioplusplus_core_<backend>.so.<abi>` and `MESHIOPLUSPLUS_ABI_VERSION` is now 15, while the C library keeps `SOVERSION 0` ([ABI policy](./abi.md)) — confirm the recipe finds both, and that a `find_package(meshioplusplus X.Y.Z EXACT CONFIG)` consumer against a Spack-installed `+cxx_api` still resolves ([C++ API](./cpp_api.md)). **S**
- **Variants that lag the CMake options.** The recipes expose `hdf5`, `netcdf`, `zlib`, `zstd`, `lz4`, `kahip`, `fortran`, `cli`, `parallel`, `mesh_backend`, `cxx_api` and `cxx_api_backends`. `CMakeLists.txt` also defines `MESHIOPLUSPLUS_WITH_CGNSLIB`, `MESHIOPLUSPLUS_WITH_GIDPOST`, `MESHIOPLUSPLUS_WITH_POLYSCOPE`, `MESHIOPLUSPLUS_WITH_EIGEN` and `MESHIOPLUSPLUS_WITH_JSON`. Verify which of them post-date 9.10.0 and which have a Spack package to depend on (`cgns` does); expose the ones an HPC build would choose, leave the rest at their defaults, and keep every variant named as the Conan option and vcpkg feature are ([C API](./c_api.md#package-managers-conan-vcpkg-spack)). **S–M**
- **Test the matrix, not one install.** Build both recipes with `spack install --test=root` for the default variants, `+fortran`, `+cxx_api` (each `cxx_api_backends`), `parallel=openmp|tbb|kokkos` and `mesh_backend=kratos`, and run `spack style` and `spack audit` before opening the PR. None of this is in CI here today; a smoke job that installs from a `spack-packages` checkout on a schedule would catch the next drift without a person remembering to look. **S–M**
- **Make the bump a release step.** [Installation → Spack](./installation.md#spack) and the C API page both say a new release "needs no action here" because a `version(...)` line is added upstream after each tag; the drift above is what that sentence produces when it is nobody's job. Add the upstream PR to the release checklist in `AGENTS.md` (after the tag: `spack checksum`, one `version()` line per recipe, `spack style`) and replace the sentence in both docs. Opening the PR from the release workflow needs a token on a fork and is a follow-up only if the manual step is skipped again. **S**
- **Done when.** Both recipes list the current release, `spack install py-meshioplusplus@<current> +hdf5 +netcdf +zlib` and `spack install meshioplusplus@<current> +fortran +cxx_api` succeed on a clean Spack, the docs no longer claim the step is automatic, and issue #3 is closed.
- **References.** [Spack packaging guide](https://spack.readthedocs.io/en/latest/packaging_guide_creation.html) † · [`spack checksum`](https://spack.readthedocs.io/en/latest/command_index.html#spack-checksum) † · [PR #5624, the original submission](https://github.com/spack/spack-packages/pull/5624) · [Installation → Spack](./installation.md#spack)

---

## 3. Quality of implementation

*Admission: work that makes every other item safer to land. None of it is a feature, so none of it competes for the same attention — it can run in parallel with everything.*

- **A sanitizer CI leg** — ASan and UBSan over the existing `cpp-tests` job. No workflow passes `-fsanitize` today, and a fuzzer that finds a crash without one reports a symptom rather than the out-of-range read behind it. The precondition for the next item. **S**
- **Fuzzing the readers** (libFuzzer, then OSS-Fuzz if the project is accepted). 43 mostly hand-rolled parsers are reachable from a C ABI, a browser, a VS Code extension and an MCP server — untrusted input reaches them by design. One fuzz target per `registry_readers()` entry, seeded from `tests/python/meshes/`. The highest-value non-feature item in this document. **M**
- **A format conformance matrix** — one canonical mesh written to and read back from every writable format, asserting per format what survives (points, each cell type, point/cell/field data and their dtypes, each region kind) against a declared expectation. `tests/python/test_region_roundtrip.py` already does this for regions over Gmsh/Abaqus/MED, and `tests/cpp/test_sequence.cpp`'s `WriteSupportsTimeAgreesWithReality` is the registry-iterating shape to generalise it to. The declared expectations become a lossiness column in the [format table](./formats.md), which today has only Read/Write/dependencies, with lossiness scattered across its notes and fifty-five per-format quirks sections. The canonical mesh should be a primitive constructor from [§6](#_6-operations). **M**
- **Property-based testing** (Hypothesis) over the invariants the docs already articulate: partition-of-unity, volume conservation, conformity, byte-identical determinism, map composition. **M**
- **A benchmark harness that covers what ships, with a CI leg.** The suite exists (`benchmark/`, up to ~1M synthetic tets in Python, 257k in the C++ backend benchmark, which is off by default) but no CI job runs any of it, so a performance regression is found by a user. It is also narrow: `benchmark/bench.py` times 6 format labels of the 43 the core reads, and [benchmarks](./benchmarks.md) has no numbers for any operation or for the parallel backends — `src/cpp/benchmark/bench_backends.cpp` compares mesh backends only. Widen `bench.py` to every registry format; add a `bench_ops.cpp` (`extract_surface`, `smooth`, `refine`, `merge`, `clean`, `compute_sdf`, `decimate`, `partition`, `reorder`) over a size sweep and SEQ/OpenMP/TBB/Kokkos; add a 10M+ cell tier; run it on a schedule that records rather than gates. Every [§4](#_4-performance) item is gated on this showing its before/after, and it decides whether the scale items in [§8](#_8-long-run-spike-first) matter at all. **S–M**
- **Test and install the ParaView plugin.** `tools/paraview-meshioplusplus-plugin.py` ships as a reader and writer, but nothing tests it, and the `data_files` entry that would install it is commented out in `pyproject.toml`, so [its page](./paraview_plugin.md) describes a plugin path nothing writes. A `pvpython` smoke step plus the install fix. **S**
- **Finish the fallback narrowing.** The per-format shims route every decline through `core_declined` (`_fallback.py`), but two halves of the same defect remain. First, **44 broad `except Exception` handlers in 36 package-root operation files** (`_clean.py`, `_data_average.py`, `_curvature.py`, …) still wrap the C++ core and, unlike a format fallback, silently substitute a *different algorithm*; each operation needs its own ruling on which exceptions mean "unsupported" and which mean "bad input", with `_error.py`'s `(ValueError, TypeError)` re-raise as the model. Second, the C++ core has one `ReadError` for both "malformed file" and "construct I deliberately decline", and about a hundred `std::stoi`/`stoll`/`.at()` call sites in the format readers leak `ValueError`/`IndexError` instead (`tests/python/meshes/tecplot/quad_zone_space.tec` makes the C++ Tecplot reader throw `std::stoull`, and now logs a warning before the Python twin reads it); a `ReadError` subclass, or a per-entry-point wrapper in `_core.cpp` (not the global translator, which would remap the operations' `std::invalid_argument`), makes the "recognised decline" contract total. Running the suite once under `MESHIOPLUSPLUS_STRICT_CORE=1` sizes both. **M**

---

## 4. Performance

*Admission: a measured or code-verified slowdown in a path a user hits, with the shape of the fix named. Nothing here is scheduled before the [§3](#_3-quality-of-implementation) harness can show its before/after.*

Two findings frame the section. First, **the serial phases below are deliberate**: each is documented in the code as a determinism pin, not an oversight — output is byte-identical across parallel backends and thread counts, and the reference-file tests enforce it — so every fix must keep that guarantee and prove it with a SEQ-versus-OpenMP diff, not assert it. Second, **every parallel item is conditional on the backend**: a SEQ build (and the `stl` fallback without TBB) runs `parallel_for` sequentially, so each change must also show that SEQ does not get slower — a parallel sort is O(n log n) where the hash map it replaces is O(n).

**Measured regressions** (verify first: `benchmark/results.csv` is one run on one machine). Gmsh binary write at 0.68× legacy meshio and MED read at 0.81–0.93×, both reported as parity or better in [benchmarks](./benchmarks.md) until this change; XDMF (HDF5) write at 0.92–0.97×; XDMF read at 1.0× on a single-block mesh against 10× on mixed topology, which suggests the post-read conversions — `xdmf.cpp` has no `parallel_for` at all. Re-measure on the widened harness, then profile. **S each to size**

**Text I/O.**

- **A shared tokenizer and number path.** 29 readers split each line with `std::istringstream` and `>>` into a `std::vector<std::string>` — one stream and N heap allocations per line (`su2.cpp`, `vtk_read.cpp`, `mdpa.cpp`, `avsucd.cpp`, `tecplot.cpp`, `unv.cpp`, `flac3d.cpp` and 22 more) — while `gmsh.cpp`'s `GmshCursor`, a `string_view` cursor calling `strtod` straight on the buffer, is the in-repo model to copy. Built on the shipped `detail/fast_number.hpp` (locale-independent `parse_double`/`snprintf_c`) and migrated reader by reader against the reference files; every one of those streams is already pinned to the classic locale by `detail/classic_stream.hpp`, so what remains here is performance only. The worst cases: XDMF ASCII `DataItem`s (a `std::string` and a dtype switch per scalar), VTU ASCII arrays (`push_back` with no `reserve`, then a second pass with a dtype switch per element) and `mdpa.cpp`, which materializes one `std::string` per line of the whole file before parsing. **M**
- **Hoist the per-element dtype switch.** `detail::dispatch_dtype` (`detail/value_io.hpp`) exists to move a `DType` switch out of a hot loop, and is used in 14 source files against nearly 500 per-element `read_double`/`read_int`/`read_point` call sites. In I/O the hot ones are the VTU ASCII writer and `vtu_to_int64` (`detail/vtk_xml.cpp`) and the Exodus reader's index shift, coordinate transpose and `column_stack` (`exodus.cpp`), which are also fully serial — the exact treatment MED already had, which took it from 0.2–0.6× to parity. **S each**
- **Parallel row formatting in ASCII writers.** `abaqus.cpp` and `ansysinp.cpp` format rows in parallel into one string per row and stream them in order, with byte-identical output; `vtu`, `vtp`, `vtk`, `medit` and `tecplot` still format serially with a locale-aware `ostream << int` per index, and the OpenFOAM writer formats every coordinate through `ostream << std::setprecision(16)`, the slowest formatting route in the standard library. The reference files are the gate, not the claim — `setprecision(16)`'s exact output has to be reproduced. **S per writer**

**Binary I/O (VTU).**

- **`b64decode` is serial and pushes one byte at a time** (`detail/vtu_binary.cpp`), while `b64encode` beside it is parallel; decode is the VTU read hot path. Its inverse table is guarded by a hand-rolled `static bool init`, which is a data race as soon as decoding is parallel or the GIL is released — replace it with a magic static in the same change. A branchless resize-then-index loop first, then parallel chunks after a whitespace pre-scan (chunks are not independent without one). **S**
- **The VTU binary read copies each payload five or six times**: `vtu_strip` copies the base64 text and then `substr`s it, `vtu_parse_binary` decodes into a `std::vector` and `memcpy`s that into the array, `vtu_decode_uncompressed` builds another vector to strip the header, and the codec returns one `std::vector` per 32 KiB block which is then copied again (on write, the blocks are concatenated with no `reserve`). Decode straight into the destination array; independent `<DataArray>`s can then decode in parallel. It roughly halves peak memory too. **M**
- **Raw `<AppendedData>` is not supported** in either direction. It is the VTU encoding with no base64 at all, so files meshio++ writes could skip both items above; a format-reach item as much as a performance one. **M**

**Memory and allocation.**

- **Small, mechanical:** five readers slurp the file by hand instead of going through `detail/file_source.hpp` (`ply`, `medit`, `ansys`, `wkt`, `stl`), and several accumulate without `reserve` although the count is in the header (`gid_read.cpp`, `openfoam.cpp`, `ansysinp.cpp`, `unv.cpp`). **S**
- **Per-cell heap allocations in every polyhedral reader.** `AddPolygonBlock`/`AddPolyhedronBlock` take nested `std::vector`s, so the CGNS, EnSight, MED, OpenFOAM, UNV, FLAC3D and VTU readers allocate one vector per cell and one per face; the VTU path additionally sorts each cell's nodes and buckets cells through a `std::map`. A CSR ingestion overload — flat node ids plus face and cell offsets, the shape the WASM binding already crosses with — removes all of it at once. It changes the installed API, so it pairs with an ABI bump. **M–L**

**Serial phases inside parallel operations.**

- **One shared, deterministic facet and edge table.** Seven places fill facet or edge keys in parallel and then deduplicate them through a single-threaded `unordered_map` (`surface.cpp`, `smooth.cpp`, `partition.cpp`, `convert_cells.cpp`, `refine.cpp`, `detail/marching.cpp`, `detail/surface_distance.cpp`) — on a 10M-tet mesh, ~40M hash operations on one core while the rest idle; it is the dominant phase of `extract_surface`, `extract_skin`, `smooth`, `partition`, `refine` and `elevate`. Nine independent key and hash types exist for it, and `detail/face_mesh.hpp` has only four consumers. `refine.cpp` already states the property a parallel version needs — numbering is a pure function of (block, cell, slot) — so a parallel sort of (key, slot) followed by a segmented first-occurrence scan reproduces it exactly there; the other six sites need that argument made before the claim. `parallel.hpp` has only `parallel_for`, so the first step is a `parallel_sort`/`parallel_reduce` primitive; the second is caching the table so an N-step pipeline stops rebuilding it N times. **L**
- **`optimize_volume` rebuilds everything on every sweep** — a fresh `Mesh` with new point and connectivity arrays, `smooth`'s node adjacency and boundary hash, and three single-threaded `unordered_map`s for the 2-3 and 3-2 flips — with one `parallel_for` in the whole file. Hoist the adjacency out of the sweep and maintain the face and edge maps across flips. **S–M**
- **Welding is single-threaded.** `clean` keeps its own point grid rather than `detail::SpatialGrid` (its header says so), reads each coordinate through the dtype switch, and deduplicates cells with a `std::string` built per cell as the hash key, alongside three per-cell vectors; `merge` scans 27 buckets per point on one core. Move to `SpatialGrid`, fixed-size keys and a chunked scan merged keep-first. **M**
- **Distance-kernel construction is serial** — the triangle soup, edge table, pseudonormals and grid insertion (`detail/surface_distance.cpp`) — while its queries are parallel; `compute_sdf`, `shrinkwrap`, `voxelize`, `compute_curvature` and `repair` all pay it. A BVH is deliberately *not* the fix: `detail/surface_distance.hpp` rejects one because its visit order would make the equidistant-triangle tiebreak observable. **M**
- **`reorder` sorts an index vector through an indirect comparator** after a serial bounding-box scan, and `partition`'s space-filling-curve path does the same; a radix sort on the 64-bit curve key replaces both. **S**
- **Hoist the dtype switch in operations**, as in I/O: `refine` (29 call sites), `interpolate` (21), `convert_cells` (13), `repair` and `diff` (7 each), `partition` (6). **S each**
- **Operations with no parallel phase at all**: `remesh` (serial clustering), `remesh_volume` (serial cut), `agglomerate` (a `std::set` and an `unordered_map` allocated per seed), `split`, `undo_green`, and `hessian` (two `gradient` passes plus deep copies). `compute_quality` (chunked partials merged in chunk order) and `sobolev_deform` (a gather-form sparse product) are the in-repo models. **L in total**
- **Radius and k-nearest search in the core.** `proximity_graph` is numpy-only — its own docstring measures 200k points at 6.4 s for a radius graph and 23 s for k=16 on one core — while `detail::SpatialGrid` already has the needed primitives. Expose the neighbour search as a core operation the Python layer calls, keeping graph assembly in Python per [Non-goals](#non-goals-and-decisions-taken). **M**

**Boundaries and startup.**

- **The Python bindings never release the GIL.** Nothing in `bindings/python` uses `gil_scoped_release` or a `call_guard`, so a multi-gigabyte read or a 10M-cell operation blocks every other Python thread, and no caller can convert files in a thread pool. Release after the numpy→`Mesh` conversion and re-acquire before the `Mesh`→numpy one, once the `b64decode` table race above is fixed, auditing that no released region touches a Python object. **S–M**
- **`import meshioplusplus` takes ~180–210 ms, ~160–190 ms of it the CLI.** The package `__init__` imports `_cli`, which imports every verb module, and `_common.py` imports `rich` for library use; every MCP process and CLI call pays it. Load `_cli` lazily through a module `__getattr__` and move `rich` into the CLI. The heavy optional dependencies (`h5py`, `netCDF4`, `torch`, `pxr`, `vtk`) are already imported inside functions. **S**
- **A declined C++ read costs a full parse before the Python one starts.** `vtu_read.cpp` builds the whole XML DOM, base64 bodies included, before rejecting lzma, appended data or multiple pieces; `gmsh.cpp` rejects `$Periodic` only after `$Nodes` and `$Elements` are parsed; with an ambiguous extension one `.msh` can be parsed up to six times. A cheap pre-flight — attribute and section-header scans — before the expensive parse, with only early rejections falling back. Every decline is now logged at DEBUG (and raises under `MESHIOPLUSPLUS_STRICT_CORE=1`), which is what shows how often it happens; the durable fix is [§5](#_5-core-parity-across-surfaces)'s core parity. **S–M**
- **The MCP server re-reads the input file on every tool call** — 65 call sites in `mcp/_tools.py` go through an uncached `_load()`, so an agent's info → clean → decimate → convert parses one file four times. A bounded cache keyed on (path, `mtime_ns`, size); anything weaker manufactures a stale-read bug. **S**
- **pybind11 per-call overheads**, together: the `Mesh`→numpy conversion re-imports the `Mesh` class on every call, the contiguity check does a Python attribute lookup per array, operations clone connectivity they never change (`transform`), and polygon/polyhedron blocks cross the boundary one node id and one face at a time where the WASM binding already uses a CSR triple. **S–M**
- **Flat-binding accessor copies**, together: R copies the points twice and shifts connectivity to 1-based with a scalar loop, and Julia's safe accessors `copy` the borrowed view; both are documented, and both are fixable behind the same accessor names. **S**
- **The browser viewer round-trips every operation result through a VTP file** — written into MEMFS, copied out with `.slice()`, then parsed again by vtk.js on the main thread. A typed-array mesh channel between worker and renderer removes three full passes. **L**

**Deliberately not**, and recorded so it is not re-proposed:

- **Explicit SIMD intrinsics or `-march` flags** — portability across wheels, WASM and the release binaries is worth more than the scalar kernels cost; revisit only if the harness shows a kernel dominating.
- **Kokkos device execution** — every `parallel_for` body captures host pointers (`parallel.hpp`); the GPU route is the DLPack/CuPy handoff ([GPU handoff](./gpu.md)).
- **A BVH in place of the uniform grid** — the tiebreak argument above.
- **Tuning the pure-Python fallback readers** — the fix is making the C++ path accept the file ([§5](#_5-core-parity-across-surfaces)), not a faster fallback; `_decimate.py`'s heap-based twin is deleted once `decimate` accepts its inputs, not optimised.

*Recommended posture:* the harness and the measured regressions first; then the small isolated wins — `b64decode`, the lazy CLI import, `optimize_volume` — then the shared facet table, the largest total win; the boundary items as their consumers ask.

---

## 5. Core parity across surfaces

*Admission: something the Python layer can do that the C++ core cannot, or that the core can do and a binding cannot reach.* A construct that forces the Python fallback is not "slower from C" — it is **unreadable** from C, Fortran, Julia, R, WASM and the native CLI, none of which has a fallback. Ordered by this project's own consumers, Kratos first.

- **MDPA beyond mesh-level blocks.** The C++ core reads and writes `Nodes`/`Elements`/`Conditions`/`SubModelPart`s, but `Begin Table`, `Begin Geometries`, `Begin Mesh <id>`, `Begin Constraints` and non-numeric `ModelPartData` throw — or, under a lenient read, are skipped and listed in `MdpaInfo`, which no flat binding exposes ([MDPA](./formats/mdpa.md#c-core)). **M**
- **Gmsh `$Periodic` and format 4.0 in the C++ core.** `$Periodic`, both directions: a periodic 4.1 file is unreadable from every flat binding today ([Gmsh](./formats/gmsh.md)). The C++ reader also accepts only versions 2.2 and 4.1 (`gmsh.cpp`), so a 4.0 file, which the Python reader reads, is unreadable from them too. Pairs with periodic node matching in [§6](#_6-operations). **S–M**
- **VTK-family constructs the C++ readers refuse**: multi-`<Piece>` `.vtu` (the Python reader merges pieces) and legacy `.vtk` structured points, structured grid and rectilinear grid ([VTU](./formats/vtu.md), [VTK](./formats/vtk.md)). **S–M**
- **XDMF 2 and XPath references** — the C++ core implements XDMF 3 only, and `Reference="XML"` `DataItem`s not at all ([XDMF](./formats/xdmf.md)). **M**
- **MED multi-mesh files and profiles**, which are Python-only and not reachable even under a lenient C++ read ([MED](./formats/med.md)). **M**
- **Netgen extras** — periodic `identifications`, `materials`/`bcnames`/`cd2names`/`cd3names`, `edgesegmentsgi2` and the `.vol.gz` container ([Netgen](./formats/netgen.md)). **S–M**
- **An Exodus writer that carries sets and steps.** It writes element blocks but no node sets or side sets, so only element-block regions round-trip, and it writes one step per file; a multi-step writer is a stateful object of the `XdmfTimeSeriesWriter` shape ([Exodus](./formats/exodus.md)). **M**
- **Sets → regions, phase 2.** `ansysInp` sets still travel in a side-channel struct (UNV groups became regions in v15.6.0), XDMF `Sets` are not mapped, and VTU/VTP need a documented region convention ([Named regions](./regions.md)). **M**
- **Named Side regions surviving operations.** `subdivide`, `agglomerate`, `undo_green` and `convert_cells(simplexify)` drop them through their parent-cell remap, and the cutters (`slice`, `isosurface`, `extract_surface`, `extract_skin`) drop them outright ([Named regions](./regions.md)); for a Kratos model, Side regions are where the boundary conditions live. **M**
- **A structured pipeline report on the flat ABI.** C, Fortran, Julia and R receive status plus `mio_last_error()` only; a caller-buffer JSON accessor is recorded as a follow-up, as are the v2 multi-mesh steps (`Inputs:` for `Merge`/`Interpolate`/`UndoGreen`, an `Output.Pattern` for `Split`/partition) ([pipelines](./pipeline.md)). **S–M**
- **PCD `binary_compressed` on the flat ABI.** The C++ API and Python write it (`write_pcd(..., PcdData::BinaryCompressed)`, `data="binary_compressed"`), but `WriteOptions`/`mio_write_opts` only carry the VTK block codecs, so C, Fortran, Julia, R and WASM can read it and cannot write it; the fix is an appended `MIO_CODEC_LZF` enumerator routed to `write_pcd` ([PCD](./formats/pcd.md)). **S**
- **glTF write options on the flat ABI.** `mio_write("x.glb")` and `writeMesh` write the defaults; the colour field, colormap and range, split angle, up axis and unit scale are reachable from Python, C++, both CLIs and MCP only, so C, Fortran, Julia, R and WASM cannot colour a `.glb` or change its axis ([glTF](./formats/gltf.md)). The fix is a `mio_gltf_opts` struct with a reserved tail, following the `mio_curvature_opts` shape. **S**
- **Point/cell sets beyond regions in the core**, so the `convert -s/-d` sets↔data conversions work in the native CLI and flat bindings ([Julia](./julia.md)). **S–M**

---

## 6. Operations

*Admission: a new operation, or a public face for machinery that already exists privately inside one.*

**Generation.** Almost every operation transforms a mesh you already have. The exceptions all start from something else — `grid` from a lattice (`detail/grid_lattice.hpp`), `voxelize`/`compute_sdf` from a surface's bounding box, `remesh_volume` from a closed surface — and nothing builds a shape from parameters, sweeps one, or triangulates a domain.

- **Primitive constructors** — `box`, `sphere`, `cylinder`, `disk`, in their own `operations/primitives.hpp` beside `grid`. Dependency-free, and it removes the fixture-file dependency from tests, docs, notebooks, the browser demo and the MCP server; it is also the canonical mesh the [§3](#_3-quality-of-implementation) conformance matrix needs. Highest leverage per line of code in this document. **S**
- **`extrude`** — 2-D → 3-D sweep (triangle → wedge, quad → hexahedron) with `nlayers` and per-layer offsets, carrying regions to side and cap regions. The most-requested generation primitive. **M**
- **`revolve`** — `extrude`'s rotational sibling around an axis, sharing its layer machinery; degenerate cells on the axis are the only new work. **M**
- **Delaunay / constrained 2-D meshing** — genuinely useful, but robust geometric predicates are where dependency-free stops paying. Better as an optional Triangle or Gmsh backend, off by default, following the KaHIP pattern. **L**

**Analysis and editing.**

- **Feature edges as a line mesh.** The feature-angle crease test exists three times (`decimate`, `smooth`, `remesh` each carry `mFeatureAngleDeg`) and only ever pins nodes; one public op emitting `line` cells serves inspection, boundary-condition picking and those three in one place. **S**
- **Hausdorff distance** between two meshes — a symmetric max-reduction over the shipped `distance_to_surface`, returning the scalar the remesh/decimate tests and the conformance matrix want to assert on. **S**
- **Periodic node-pair matching** — given two boundary regions and a transform, return the matched node pairs. `$Periodic` already round-trips as metadata and `proximity_graph` already does minimum-image search; Kratos periodic conditions are the consumer. **S–M**
- **A quality gate** — pass/fail thresholds over the metrics and histograms `compute_quality` already produces, a `check` CLI verb that exits non-zero, and uniform `--json` output (today only a handful of verbs take it — not `info`, `quality`, `diff` or `convert`). What a CI pipeline over meshes actually scripts. **S–M**
- **Region set algebra** — union, intersection, difference, rename and retag of named regions. Regions are a first-class layer with add/replace/enumerate only; every consumer that builds boundary conditions from them hand-rolls this. **S**
- **Time-axis resampling of sequences.** The sequence engine reads time values and drives N→N, fan-in and fan-out, but never resamples; aligning two solvers' timelines (or a solver and a surrogate's) is the missing step before a pairwise `diff` or a training pair. **S–M**
- **`agglomerate` follow-ups** — coplanar boundary-face merging (fusing adjacent group faces on one plane into a single polygon) and a shape-quality absorption gate, both deferred when it shipped ([agglomerate](./agglomerate.md)). **S–M**

---

## 7. Ecosystem reach

*Admission: getting what exists to the people who would use it.*

- **Registration and distribution.** All pure logistics and all blocking real adoption; the effort is small per target but the calendar time is set by each registry's review, so start early. **S each, calendar-bound**

| Target | Status today |
| --- | --- |
| conda-forge | Not started. |
| CRAN | `R CMD check --as-cran` already runs in CI, so the gate exists; only submission remains ([R](./r.md)). |
| Julia General | Not registered, so `Pkg.add("MeshioPlusPlus")` does not work yet ([Julia](./julia.md)). |
| Blender Extensions Platform | Extension zips are built and attached to every release; only the listing remains ([Blender](./blender.md)). |
| ConanCenter and the vcpkg registry | Recipes are self-hosted and CI-validated under `packages/`; neither is submitted ([C API](./c_api.md#package-managers-conan-vcpkg-spack)). |
| ParaView | The plugin exists; installing it with the wheel is the [§3](#_3-quality-of-implementation) item. |

- **Rust bindings** over the C API — the next language by scientific adoption after Julia and R, and the ABI/`SOVERSION` work makes it cheap. **M**
- **Interop phase 2** — the Open3D and DOLFINx bridges exist as named stubs that raise `NotImplementedError` ([interop](./interop.md)); the pinned-memory staging for CuPy is wired in C++ (v8.5.0) but not from Python ([GPU handoff](./gpu.md)). **S–M each**
- **A CLI chatbot** — `meshioplusplus chat`, a natural-language entry point from the terminal. The MCP server already exposes the whole Python surface to an agent, but only to a host that speaks MCP (Claude Desktop, an IDE); a user with an LLM API key cannot ask the CLI "why does this file fail to convert" or "clean this mesh and tell me what changed". A thin client over existing work, not new mesh functionality. **S–M**
  - **Reuse, not reimplementation.** The tool-calling loop drives `mcp/_tools.py`'s `TOOL_REGISTRY` in-process — never a second copy of its dispatch, path sandbox or result sanitizing — so every tool added there reaches the chatbot for free, under the parity guarantee `test_every_operation_has_a_tool` already gives MCP. The pattern is proven, not hypothetical: `mcp/_http.py` already dispatches the dashboard's `/api/tools/<name>` through the same registry.
  - **v1 scope.** A stateless-per-turn loop (history + tool schemas out, tool calls executed, results fed back until a plain-text reply) over files under the CWD or `--root`, reusing the MCP sandbox; context lives only in the terminal session — no persistence, no server.
  - **Packaging.** A `chat` extra carrying the model provider's SDK, kept out of `[all]` and imported only by the new `_cli/_chat.py`, so `_tools.py`, the default CI matrix and every other surface are untouched; without the extra or an API key the verb fails by name, following the `meshioplusplus-mcp` entry-point precedent. Provider-agnostic wiring and persisted history are follow-ups, not v1.
  - **Docs.** A `doc/chat.md` page modelled on [MCP](./mcp.md), plus the README and tool-table updates the `AGENTS.md` checklist requires.

---

## 8. Long run (spike first)

*Admission: work whose shape is unknown until an investigation writes it down. Findings before code.*

**Scale.** Memory-mapped reads ([mmap](./mmap.md), via `ReadOptions` and the C ABI but not Python's `read()`) roughly halve the peak footprint of a large read, and the XDMF series appender and Python's chunked `write_dataset` write a series or a dataset without holding it — but nothing writes one mesh larger than memory, and no operation streams. Run the [§3](#_3-quality-of-implementation) benchmark tier first; it decides whether either item below matters.

- **Streaming / chunked writes of one mesh**, the counterpart to selective and memory-mapped reads. **L**
- **Out-of-core operations** for the ops that are already block-local. **XL**

**NURBS and higher-order geometry.** The data model is strictly linear/Lagrange polytopes: a `CellBlock` is a cell-type string plus a node-index array. NURBS is a genuinely different object — control points, weights, knot vectors, and a parametric mapping — and CAD/IGA formats (STEP, IGES, Rhino 3dm, `.iga`) express geometry that no current cell type can hold. This is the most architecturally invasive item on the list and should be approached as a research spike, not a feature.

- **Spike: how far can the current model stretch?** Higher-order Lagrange cells already exist (`hexahedron27`, VTK-Lagrange types); a rational Bézier/NURBS patch needs *weights* and a *knot vector*, which have nowhere to live. Determine whether a side-channel struct (the `MedInfo`/`GmshInfo` precedent) suffices, or whether the `Mesh` needs a genuine second entity kind. Write the finding up before committing. **M**
- **Read-only CAD ingestion first**: a NURBS surface tessellated to a triangle mesh at a requested tolerance, with the parametric data carried out-of-band. This delivers most of the practical value (getting CAD into the mesh world) without touching the data model, and is the natural first release. **L**
- **A real IGA data model** — patches, control nets, weights, knots, trimming curves — plus formats and evaluation. This is XL, likely a separate library or a major version, and should only be attempted if the spike shows real demand.
- **Dependency reality**: robust STEP/IGES parsing effectively means OpenCASCADE, which is a heavyweight LGPL dependency. If ingestion goes ahead, it must follow the KaHIP/Polyscope pattern — strictly optional, off by default, never in the core, licence implications documented. **Findings before code.**
- **No C++ tessellator exists yet** — `tessellate` is Python-only — so read-only CAD ingestion reaching every surface also carries that port.

*Recommended posture: spike and document; do not schedule implementation until the spike says what shape it takes.*

---

## Non-goals and decisions taken

Recorded so they are not re-proposed as gaps.

- **MPI in the library** — none planned ([C++ API](./cpp_api.md)); `partition`'s ghost layers produce the halo an MPI assembly in the owning application needs.
- **Solver-coupled physics-ML** — assembled solver residuals, adjoints and Sobolev training, co-simulation, active-learning *labeling*, adaptive remeshing driven by a surrogate, and MPI model-part gathering. Every one needs a live solver (an assembly routine, its tangent, its communicator) and meshio++ has no notion of a discrete system; they belong in the application that owns the solver, the division [Symbolic and physics](physicsnemo/symbolic_and_physics.md) describes.
- **The Python-only layers stay Python** — `pmsh`, `zarr`, `cae` and `usd`, and the physics-ML surface (`tessellate`, grids, point budgets, proximity graphs, datasets, training) are export targets and tooling for a training pipeline, registered in Python rather than the shared C++ registry. A core kernel they call (the neighbour search in [§4](#_4-performance)) does not change that.
- **Polyscope in the release CLI binaries** — excluded deliberately, so `view`/`screenshot` there report the build flag rather than opening a window ([viewer](./viewer.md)).
- **General meshing algorithms as operations** — hex and hex-dominant meshing, boolean/CSG, boundary-layer inflation, quadrangulation and geodesic distance are each a library in their own right; the answer is an optional backend (the KaHIP pattern), not a native implementation.
- **KaHIP in the WASM build** — no Emscripten port, and a graph partitioner would bloat every consumer's bundle; `"auto"` resolving to the SFC method is the answer, and `"kahip"` throwing by name is the contract.
- **zstd/lz4 codecs, memory mapping and the Kokkos backend under Emscripten** — the codecs have no Emscripten port and are compiled out, CMake refuses the Kokkos backend under Emscripten, and there is nothing to map inside MEMFS.
- **Single-file output for DOLFIN, TetGen and EnSight from WASM**, or lifting DOLFIN's simplicial restriction — both are facts of the formats; `writeMesh`'s written-paths return value (v11.2.0) is the fix for the bookkeeping they cause.
- **Polyhedron blocks as a WASM gap** — `vtu`, `ensight`, `cgns`, `med` and `openfoam` all write them from WASM; a consumer that cannot is constrained by its own data model.

---

## Suggested sequencing

Open work only; what shipped is in `CHANGELOG.md`.

1. **Format reach ([§1](#_1-format-reach))** — the first FEM wave (`.mail`, OptiStruct).
2. **Spack package upkeep ([§2](#_2-spack-package-upkeep))** — a small, mechanical catch-up (recipes, variants, one release-checklist line) that unblocks HPC users on the current release; then a checklist step, so it stays current.
3. **Sanitizer leg, then fuzzing ([§3](#_3-quality-of-implementation))** — a parallel track from day one; it does not compete for the same attention as features.
4. **Performance ([§4](#_4-performance))** — the harness and the measured regressions first, then the isolated wins (`b64decode`, the lazy CLI import, `optimize_volume`), then the shared facet table.
5. **Primitive constructors ([§6](#_6-operations))** — a few days, and a prerequisite of the conformance matrix, every demo surface and `extrude`/`revolve`.
6. **Core parity ([§5](#_5-core-parity-across-surfaces))** — MDPA first, then Side-region survival and sets → regions, then the rest by consumer demand.
7. **Registration ([§7](#_7-ecosystem-reach))** — calendar-bound, so start the submissions early and let them run alongside everything else.
8. **Long-run spikes ([§8](#_8-long-run-spike-first))** — the benchmark tier decides the scale items; the NURBS spike is scheduled independently of the rest.

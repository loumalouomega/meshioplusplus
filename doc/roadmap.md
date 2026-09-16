# meshio++ roadmap

Status at time of writing: **v11.1.0** — 43 core formats plus four Python-only physics-ML ones, thirty-eight mesh operations + five data operations, six language surfaces (Python / C / Fortran / Julia / R / WASM), two viewers plus a browser dataset manager, a Blender add-on, a ParaView plugin, an MCP server, a settings-driven pipeline engine, a dataset-manifest layer with a PhysicsNeMo adapter, and a versioned ABI (`MESHIOPLUSPLUS_ABI_VERSION` 12).

This document lists what is *not* built. Nothing here duplicates shipped functionality; where a feature partially exists, the shipped half is named and the gap is stated explicitly. Release history lives in [`CHANGELOG.md`](https://github.com/loumalouomega/meshioplusplus/blob/main/CHANGELOG.md), not here.

Effort key: **S** = days, **M** = a couple of weeks, **L** = a month or more, **XL** = a project in its own right. An item names a **probe** — the failing test that proves the gap — wherever one is cheap.

## How this file works

- **Sections are ordered, and the order is the recommendation.** Each section states what belongs in it; an item that sounds exciting does not move up for that reason. An empty section is removed, not kept as a placeholder.
- **A closed item is removed, not struck through.** Its history is the `CHANGELOG.md` entry and the feature's own `doc/` page; a partly closed item is narrowed to what remains (the `AGENTS.md` change checklist rule).
- **Defect-shaped items go to [§2](#_2-correctness-debts) regardless of size** — behaviour that loses data, mis-orients cells or fails silently is not a feature request, even when the fix and the feature are the same work.
- **An item estimated from a doc, a `.d.ts` or a changelog alone says so** ("verify first") and names its probe; the code has repeatedly been more or less capable than its description.
- **[Non-goals](#non-goals-and-decisions-taken) record decisions already taken**, with their reasons, so they are not re-proposed as gaps.

## The map

![The roadmap at a glance: open items grouped by section, shaded by effort, with dependency arrows and the items that need a design pass or a research spike first](/diagrams/roadmap_map.svg)

---

## 1. WASM parity

**The gap.** The WASM build is not a reduced port: since v8.0.0 it ships every C++ format (43 readable / 46 writable keys, HDF5, netCDF, CGNS and GiD included), both a sequential and a threaded (`_mt`) variant, and 78 bound I/O and operation functions (plus the XDMF series handle). What it lacks is measured rather than guessed. Its main downstream consumer, the Kratos MDPA Preview VS Code extension (on `@meshioplusplus/wasm` 10.21.1), records a list of *non-goals and known constraints* that are really meshio++ gaps, and every one was re-checked against v11.1.0: no format reader, writer or the WASM loader changed in between, so all of them still hold. Each item below names the downstream constraint it would lift and its **probe** — a `tests/wasm/smoke.mjs` step that fails before the fix. Like that extension's own rule, an item is not estimated from `index.d.ts` alone: `.d.ts`-level green has been wrong before (a `.med` export was type-correct and unreadable; a mismatched `locateFile` produced a `LinkError` naming neither the file nor the variant).

**Defects** — behaviour that loses data or fails opaquely; first, regardless of size. Two core-writer defects a WASM consumer is the first to hit (DOLFIN multi-block cell data, the MED 4.1 bitmask) are owned by [§2](#_2-correctness-debts), since fixing them fixes every surface.

- **Load failures name nothing.** `loadMeshioPlusPlus()` (`src/wasm/src/index.mjs`) instantiates the chosen variant with no guard, so a `locateFile` returning the sequential binary to the `_mt` glue (or a missing `.wasm`) surfaces as whatever Emscripten throws — historically a bare `LinkError`. Catch at instantiation and rethrow naming the variant, the filename Emscripten requested and the URL `locateFile` returned. *Probe:* a smoke step with a deliberately mismatched `locateFile`, asserting the message names both. **S**
- **Integer dtypes are widened to Float64 at the boundary**, so an in-memory `decimate(mesh)` blends integer tags as floats (the C++ core keeps the surviving vertex's row for integer point data — the boundary destroyed the dtype before it got there), `partition:ghost`/`refine:*`/material ids come back as doubles, and an Int64 id past 2^53 silently rounds. Carry the dtype across in both directions, following the `*_components` sibling-map precedent (`Int32Array`/`BigInt64Array`/`Uint8Array` where the source was integer). The path-based `convert`/`convertSurfaceOps`/`runPipeline` already avoid the boundary and stay unaffected. *Lifts:* the downstream "decimate blends every field including integer tags as float64" exclusion. *Probe:* decimate a surface carrying an Int32 cell tag and an integer point field; assert the returned typed-array class and exact values. **M**

**A parity guard** — a Python test that extracts the `emscripten::function` names from `bindings/wasm/js_bindings.cpp` and checks them against the `loadMeshioPlusPlus()` wrapper keys, the `MeshioPlusPlusModule` members in `index.d.ts` and the operation-shaped `mio_*` functions of the C header, with named exemptions — the `test_every_operation_has_a_tool` pattern applied to the flat surfaces. Nothing checks this today, and the doc/type drift corrected when this section was written (an `extractSkin` still described as unexposed, polyhedron writing described as impossible, wrong format keys, five functions missing from the wrapper's JSDoc) is what that absence produces. Can start in parallel with everything else. **S**

**The JS surface over C++ that already exists** — no core change, no ABI change.

- **Index maps.** The C++ result structs of `crop`, `split`, `convert_cells`, `refine`, `partition`, `clean`, `decimate` and `merge` all return `point_map`/`cell_maps`, and the WASM binding discards them; `clean`, `decimate` and `extractSkin` have not even a `record*` flag to fall back on. Return them on request (`returnMaps: true`, since they cost a copy each). *Lifts:* part of the downstream "adopting meshio++'s returned mesh" exclusion — the lost original entity ids; the other two losses it names (Kratos entity kinds, property ids) have no meshio++ counterpart and stay the consumer's business. *Probe:* `clean` with welding, asserting a composed map reproduces the input coordinates. **M**
- **Side channels.** None of the `*Info` structs (`formats/*.hpp`) crosses the boundary, because the registry's `ReadFn`/`WriteFn` carry none. Return them from `readMeshSelective` as an optional `info` object and accept them on write: `OpenFoamInfo` (patch names and types — the shared registry writer currently puts every boundary face into one `defaultFaces` patch), `MedInfo` (field units, and the constructs a lenient read dropped), `MdpaInfo` (skipped constructs, properties bodies), `AnsysInfo`/`UnvInfo` sets. *Lifts:* the downstream extension re-deriving OpenFOAM patch names itself after every write, and its lenient-read diagnostic that cannot say what was dropped. *Probe:* write an OpenFOAM case with two named patches and read the `boundary` file back. **M**
- **Per-format write options** matching the C API's `mio_write_ex`/`mio_write_opts`: binary vs ASCII, codec and gzip level, float format, VTK 4.2 vs 5.1. **S–M**
- **Report what a write produced.** `writeMesh`/`convert` return `void`, yet tetgen, triangle, ensight, dolfin, openfoam, gid (ascii) and HDF-backed xdmf write sibling files or directories into MEMFS that the caller must know to copy out. Return the list of written paths. *Lifts:* the downstream extension's own companion-file bookkeeping, which it rebuilds per format. *Probe:* `writeMesh` to `.xdmf` and `.foam`, asserting the returned list against `FS.readdir`. **S**
- **`frozen` masks for `smooth` and `decimate`**, exposed in the C++ and Python APIs only. Land them across every flat binding (C, Fortran, Julia, R, WASM) in one change, since each doc page calls the gap shared. **S**
- **A direct `decimateVolume` binding.** Today it is reachable only as a `convertSurfaceOps`/`runPipeline` step. **S**
- **A stateful sequence reader** (open / step count / time / read step), which C, Fortran, Julia and R have and WASM does not — it has only `sequenceEntries` and the two conversion helpers. **S–M**

**Core format gaps a browser consumer hits first** — C++ work that reaches every surface at once, listed here because WASM has no Python reference implementation to fall back on.

- **In-file timelines.** Only XDMF, Exodus and GiD report real time values; each fix below is independent and closes a row of the downstream extension's transient audit (which pins today's behaviour, so each lands as a visible diff there). MED: a `read_med_metadata` filling `mTimeValues` (today metadata falls back to a full default read, which rejects a multi-step field — the sequence layer then closes for free, see [sequences](./sequences.md)). CGNS: `BaseIterativeData`/`ZoneIterativeData` and one solution per step, instead of every `FlowSolution` loading into the same name. Tecplot: every `ZONE` with its `SOLUTIONTIME`/`STRANDID`, instead of stopping at the first. Gmsh: `$NodeData`/`$ElementData` time tags, instead of keeping the first section per name. EnSight: a transient `model:` wildcard plus `VARIABLE`/`TIME`, instead of throwing on `*`. *Probe (each):* a two-step fixture whose `readMetadata(path).timeValues` has length 2. **L in total, S–M each**
- **OpenFOAM beyond one polyMesh.** Read `cellZones`/`faceZones`/`pointZones` as cell/side/point regions, and write them back — the writer currently *deletes* existing zone files as stale. Select a region in multi-region cases (`constant/<region>/polyMesh`), reconstruct decomposed `processor*/` cases, and read time-directory `vol*Field`/`point*Field` dictionaries onto the mesh. *Lifts:* the downstream extension parsing time directories natively in TypeScript, and its zone and multi-region/decomposed diagnostics. *Probe:* a case with one `cellZones` entry round-tripped through `readMesh`/`writeMesh`. **M–L**
- **Named groups surviving a Gmsh 4.1 export.** `$PhysicalNames` rows are emitted only for regions with a tag ≥ 0, and membership only exists where `gmsh:dim_tags` does, so groups read from Abaqus, MED or MDPA (tag −1) vanish in a `.msh` export. Allocate physical tags for untagged regions and write the `$Entities` that carry their membership. *Lifts:* the downstream "SubModelParts do not survive a gmsh export" constraint. *Probe:* an Abaqus `ELSET` round-tripped through `.msh`. **M**
- **VTK XML structured and multiblock files.** `.vti` is the only structured XML format; there is no `.vts`/`.vtr` reader or writer and no `.vtm` at all. A `.vtm` writer is an index plus one piece file per block, with the pieces reported by the written-paths list above — the same index-plus-pieces machinery as `.pvd`/`.pvtu` in [§6](#_6-format-reach), so one implementation serves both. *Lifts:* the downstream extension's read-only multiblock and structured-grid support. *Probe:* write two blocks to `.vtm`, read them back as two meshes. **M**

**Deliberately not**, and recorded so it is not re-proposed:

- **KaHIP in the WASM build** — no Emscripten port, and a graph partitioner would bloat every consumer's bundle; `"auto"` resolving to the SFC method is the answer, and `"kahip"` throwing by name is the contract.
- **zstd/lz4 codecs, memory mapping and the Kokkos backend under Emscripten** — the codecs have no Emscripten port and are compiled out, CMake refuses the Kokkos backend under Emscripten, and there is nothing to map inside MEMFS.
- **The Python-only layers** — see [Non-goals](#non-goals-and-decisions-taken); for a browser in particular they are export targets for a training pipeline, not interchange it needs.
- **Single-file output for DOLFIN, TetGen and EnSight**, or lifting DOLFIN's simplicial restriction — both are facts of the formats; the written-paths list above is the fix for the bookkeeping they cause.
- **Polyhedron blocks as a WASM gap** — `vtu`, `ensight`, `cgns`, `med` and `openfoam` all write them from WASM; a consumer that cannot is constrained by its own data model.

*Recommended posture:* the two opaque-failure defects and the parity guard first — small, and the guard stops the documentation from drifting again; then the JS-surface items, which need no core change and no ABI bump; then the core format gaps one at a time, each against its probe, rather than as a WASM release of their own.

---

## 2. Correctness debts

*Admission: behaviour that loses data, mis-orients cells, or drops something without saying so. First regardless of size, and each is independently shippable.*

- **MED quadratic 3-D node ordering is not converted.** `tetra10`, `hexahedron20`, `pyramid13` and `wedge15` are read and written with a warning and no meshio++↔MED permutation in either engine (`src/cpp/src/formats/med.cpp`, `src/python/meshioplusplus/med/_med.py`), so they may be mis-oriented for Salome, code_aster and code_saturne — the tools the format exists for. The linear 3-D permutations beside them are the template. *Probe:* a Salome-written `tetra10`/`hexahedron20` fixture whose cell volumes are all positive after a read. **S–M**
- **Gmsh `pyramid14` (and `wedge18`) node ordering is not permuted.** Both engines' reorder tables stop at `pyramid13`/`wedge15`, while [cell types](./cell_types.md) documents `pyramid14` as `pyramid13` plus a base-centre node and gmsh's own type 14 is edge-lexicographic. *Probe:* a gmsh-generated second-order pyramid whose base-centre node lies on the base plane after a read. **S, verify first**
- **The `vtk51` registry alias writes a VTK 4.2 file.** `src/python/meshioplusplus/vtk/__init__.py` maps both `vtk42` and `vtk51` to the 4.2 writer — documented as a quirk in [VTK](./formats/vtk.md) rather than fixed. *Probe:* `write(..., file_format="vtk51")` and assert the `# vtk DataFile Version 5.1` header. **S**
- **Writers that drop data without a warning.** The C++ AVS-UCD writer keeps only the first integer cell-data array, silently (the Python writer warns); the SVG and TikZ writers skip unsupported cell blocks with no warning at all ([SVG](./formats/svg.md), [TikZ](./formats/tikz.md)). Every other writer follows warn-and-skip, and provenance's conversion assumptions sweep those warnings. *Probe:* each writer under `pytest.warns`. **S**
- **DOLFIN writes every cell-data block to the same sibling file.** `write_dolfin` loops over `CellDataNumBlocks` but names each file `<stem>_<name>.xml` without the block, so on a multi-block mesh the last block overwrites the earlier ones. Write the one simplicial block the mesh file actually holds, and pin it. *Probe:* two triangle blocks with distinct cell data; read the sibling back. **S, verify first**
- **The C++ MED writer emits no MED 4.1 field bitmask attributes** (`LEN`/`LGC`/`LNA`/…), which the Python writer does; every flat binding and the native CLI use the C++ writer with no Python to defer to. The source calls this an interoperability gap with Salome/MEDCoupling, while [MED](./formats/med.md) calls the attributes required for medfile/mdump — establish which is true with `mdump` before sizing the fix. *Probe:* a C++-written, field-carrying `.med` opened by medfile. **S–M, verify first**

---

## 3. Quality of implementation

*Admission: work that makes every other item safer to land. None of it is a feature, so none of it competes for the same attention — it can run in parallel with everything.*

- **A sanitizer CI leg** — ASan and UBSan over the existing `cpp-tests` job. No workflow passes `-fsanitize` today, and a fuzzer that finds a crash without one reports a symptom rather than the out-of-range read behind it. The precondition for the next item. **S**
- **Fuzzing the readers** (libFuzzer, then OSS-Fuzz if the project is accepted). 43 mostly hand-rolled parsers are reachable from a C ABI, a browser, a VS Code extension and an MCP server — untrusted input reaches them by design. One fuzz target per `registry_readers()` entry, seeded from `tests/python/meshes/`. The highest-value non-feature item in this document. **M**
- **A format conformance matrix** — one canonical mesh written to and read back from every writable format, asserting per format what survives (points, each cell type, point/cell/field data and their dtypes, each region kind) against a declared expectation. `tests/python/test_region_roundtrip.py` already does this for regions over Gmsh/Abaqus/MED, and `tests/cpp/test_sequence.cpp`'s `WriteSupportsTimeAgreesWithReality` is the registry-iterating shape to generalise it to. The declared expectations become a lossiness column in the [format table](./formats.md), which today has only Read/Write/dependencies, with lossiness scattered across its notes and fifty-five per-format quirks sections. The canonical mesh should be a primitive constructor from [§5](#_5-operations). **M**
- **Property-based testing** (Hypothesis) over the invariants the docs already articulate: partition-of-unity, volume conservation, conformity, byte-identical determinism, map composition. **M**
- **A large-mesh benchmark tier with a CI leg.** The suite exists (`benchmark/`, up to ~1M synthetic tets in Python, 257k in the C++ backend benchmark, which is off by default) but no CI job runs any of it, so a performance regression is found by a user. Add a 10M+ cell tier and a scheduled job that records rather than gates; it also decides whether the two scale items in [§8](#_8-long-run-spike-first) matter at all. **S**
- **Test and install the ParaView plugin.** `tools/paraview-meshioplusplus-plugin.py` ships as a reader and writer, but nothing tests it, and the `data_files` entry that would install it is commented out in `pyproject.toml`, so [its page](./paraview_plugin.md) describes a plugin path nothing writes. A `pvpython` smoke step plus the install fix. **S**

---

## 4. Core parity across surfaces

*Admission: something the Python layer can do that the C++ core cannot, or that the core can do and a binding cannot reach.* A construct that forces the Python fallback is not "slower from C" — it is **unreadable** from C, Fortran, Julia, R, WASM and the native CLI, none of which has a fallback. Ordered by this project's own consumers, Kratos first.

- **MDPA beyond mesh-level blocks.** The C++ core reads and writes `Nodes`/`Elements`/`Conditions`/`SubModelPart`s, but `Begin Table`, `Begin Geometries`, `Begin Mesh <id>`, `Begin Constraints` and non-numeric `ModelPartData` throw — or, under a lenient read, are skipped and listed in `MdpaInfo`, which no flat binding exposes ([MDPA](./formats/mdpa.md#c-core)). **M**
- **Gmsh `$Periodic` in the C++ core**, both directions: a periodic 4.1 file is unreadable from every flat binding today ([Gmsh](./formats/gmsh.md)). Pairs with periodic node matching in [§5](#_5-operations). **S–M**
- **VTK-family constructs the C++ readers refuse**: multi-`<Piece>` `.vtu` (the Python reader merges pieces) and legacy `.vtk` structured points, structured grid and rectilinear grid ([VTU](./formats/vtu.md), [VTK](./formats/vtk.md)). **S–M**
- **XDMF 2 and XPath references** — the C++ core implements XDMF 3 only, and `Reference="XML"` `DataItem`s not at all ([XDMF](./formats/xdmf.md)). **M**
- **MED multi-mesh files and profiles**, which are Python-only and not reachable even under a lenient C++ read ([MED](./formats/med.md)). **M**
- **Netgen extras** — periodic `identifications`, `materials`/`bcnames`/`cd2names`/`cd3names`, `edgesegmentsgi2` and the `.vol.gz` container ([Netgen](./formats/netgen.md)). **S–M**
- **An Exodus writer that carries sets and steps.** It writes element blocks but no node sets or side sets, so only element-block regions round-trip, and it writes one step per file; a multi-step writer is a stateful object of the `XdmfTimeSeriesWriter` shape ([Exodus](./formats/exodus.md)). **M**
- **Sets → regions, phase 2.** `ansysInp` and `unv` sets still travel in side-channel structs, XDMF `Sets` are not mapped, and VTU/VTP need a documented region convention ([Named regions](./regions.md)). This closes, for every surface at once, the side-channel gap [§1](#_1-wasm-parity) names for WASM. **M**
- **Named Side regions surviving operations.** `subdivide`, `agglomerate`, `undo_green` and `convert_cells(simplexify)` drop them through their parent-cell remap, and the cutters (`slice`, `isosurface`, `extract_surface`, `extract_skin`) drop them outright ([Named regions](./regions.md)); for a Kratos model, Side regions are where the boundary conditions live. **M**
- **A structured pipeline report on the flat ABI.** C, Fortran, Julia and R receive status plus `mio_last_error()` only; a caller-buffer JSON accessor is recorded as a follow-up, as are the v2 multi-mesh steps (`Inputs:` for `Merge`/`Interpolate`/`UndoGreen`, an `Output.Pattern` for `Split`/partition) ([pipelines](./pipeline.md)). **S–M**
- **Point/cell sets beyond regions in the core**, so the `convert -s/-d` sets↔data conversions work in the native CLI and flat bindings ([Julia](./julia.md)). **S–M**

---

## 5. Operations

*Admission: a new operation, or a public face for machinery that already exists privately inside one.*

**Generation.** Almost every operation transforms a mesh you already have. The exceptions all start from something else — `grid` from a lattice (`detail/grid_lattice.hpp`), `voxelize`/`compute_sdf` from a surface's bounding box, `remesh_volume` from a closed surface — and nothing builds a shape from parameters, sweeps one, or triangulates a domain.

- **Primitive constructors** — `box`, `sphere`, `cylinder`, `disk`, in their own `operations/primitives.hpp` beside `grid`. Dependency-free, and it removes the fixture-file dependency from tests, docs, notebooks, the browser demo and the MCP server; it is also the canonical mesh the [§3](#_3-quality-of-implementation) conformance matrix needs. Highest leverage per line of code in this document. **S**
- **`extrude`** — 2-D → 3-D sweep (triangle → wedge, quad → hexahedron) with `nlayers` and per-layer offsets, carrying regions to side and cap regions. The most-requested generation primitive. **M**
- **`revolve`** — `extrude`'s rotational sibling around an axis, sharing its layer machinery; degenerate cells on the axis are the only new work. **M**
- **Delaunay / constrained 2-D meshing** — genuinely useful, but robust geometric predicates are where dependency-free stops paying. Better as an optional Triangle or Gmsh backend, off by default, following the KaHIP pattern. **L**

**Analysis and editing.**

- **`compute_normals`** as attachable point and cell data. Angle-weighted pseudonormals are already computed inside `distance_to_surface` and never exposed; a public op is the prerequisite for glTF export ([§6](#_6-format-reach)) and for shading in the viewers. **S**
- **Feature edges as a line mesh.** The feature-angle crease test exists three times (`decimate`, `smooth`, `remesh` each carry `mFeatureAngleDeg`) and only ever pins nodes; one public op emitting `line` cells serves inspection, boundary-condition picking and those three in one place. **S**
- **Hausdorff distance** between two meshes — a symmetric max-reduction over the shipped `distance_to_surface`, returning the scalar the remesh/decimate tests and the conformance matrix want to assert on. **S**
- **Periodic node-pair matching** — given two boundary regions and a transform, return the matched node pairs. `$Periodic` already round-trips as metadata and `proximity_graph` already does minimum-image search; Kratos periodic conditions are the consumer. **S–M**
- **A quality gate** — pass/fail thresholds over the metrics and histograms `compute_quality` already produces, a `check` CLI verb that exits non-zero, and uniform `--json` output (today only a handful of verbs take it — not `info`, `quality`, `diff` or `convert`). What a CI pipeline over meshes actually scripts. **S–M**
- **Region set algebra** — union, intersection, difference, rename and retag of named regions. Regions are a first-class layer with add/replace/enumerate only; every consumer that builds boundary conditions from them hand-rolls this. **S**
- **Time-axis resampling of sequences.** The sequence engine reads time values and drives N→N, fan-in and fan-out, but never resamples; aligning two solvers' timelines (or a solver and a surrogate's) is the missing step before a pairwise `diff` or a training pair. **S–M**
- **`agglomerate` follow-ups** — coplanar boundary-face merging (fusing adjacent group faces on one plane into a single polygon) and a shape-quality absorption gate, both deferred when it shipped ([agglomerate](./agglomerate.md)). **S–M**

---

## 6. Format reach

*Admission: a format a simulation or physics-ML workflow actually exchanges, or the missing half of a shipped one.*

- **VTKHDF** — Kitware's HDF5-based VTK file format, read natively by ParaView and designed for large, transient and partitioned unstructured data, which the XML formats handle with one file per piece and step. HDF5 is already linked into every build that has MED/CGNS, so this is a reader and writer, not a new dependency. **M**
- **`.pvd` and `.pvtu`/`.pvtp`** — the on-disk ParaView face of the sequence engine (a `.pvd` is a time-indexed collection) and of `partition` (a `.pvtu` is one piece per part plus an index). Shares the index-plus-pieces machinery of the `.vtm` item in [§1](#_1-wasm-parity). **S–M**
- **Point-cloud formats** `.xyz` and `.pcd` — `select_points`, `subsample_points` and `proximity_graph` form a point-cloud path with no point-cloud file at either end of it. **S**
- **LS-DYNA keyword input `.k`** — the one major solver-input keyword format missing beside Abaqus, Nastran and ANSYS; the same ASCII reader shape, with `*PART` as regions. **M**
- **CalculiX `.frd` results (read)** — CalculiX takes Abaqus-style `.inp` input, which meshio++ already writes; reading its results file closes the loop for an open-source solver. **S–M**
- **glTF / `.glb` (write)** — a web-native skin export for the browser viewer, dashboards and anything downstream of the Blender path; needs `compute_normals` first. **S–M**
- **Half-gaps of shipped formats** — OpenFOAM *binary* write (binary read exists) **S**; SU2 multizone (`NZONE`) **S**; Tecplot's non-transient multiple zones and EnSight variables are the non-timeline halves of the [§1](#_1-wasm-parity) timeline item and land with it.

*Considered, not queued:* 3MF, X3D, DXF, Alembic, OpenVDB and Silo (graphics, VFX or lab-specific rather than simulation interchange); Nastran OP2 and Tecplot binary (large specifications better served by pyNastran and Tecplot's own tools); LAS; OBJ materials (`.mtl`); gmsh `$PartitionedEntities` (Kratos partitions through MDPA); ICP registration; mesh duals; cylindrical/spherical coordinate transforms; shell completions. Revisit any of them when a consumer asks with a file in hand.

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

Recorded so they are not re-proposed as gaps. Surface-specific decisions stay with their section (WASM's are in [§1](#_1-wasm-parity)).

- **MPI in the library** — none planned ([C++ API](./cpp_api.md)); `partition`'s ghost layers produce the halo an MPI assembly in the owning application needs.
- **Solver-coupled physics-ML** — assembled solver residuals, adjoints and Sobolev training, co-simulation, active-learning *labeling*, adaptive remeshing driven by a surrogate, and MPI model-part gathering. Every one needs a live solver (an assembly routine, its tangent, its communicator) and meshio++ has no notion of a discrete system; they belong in the application that owns the solver, the division [Symbolic and physics](physicsnemo/symbolic_and_physics.md) describes.
- **The Python-only layers stay Python** — `pmsh`, `zarr`, `cae` and `usd`, and the physics-ML surface (`tessellate`, grids, point budgets, proximity graphs, datasets, training) are export targets and tooling for a training pipeline, registered in Python rather than the shared C++ registry.
- **Polyscope in the release CLI binaries** — excluded deliberately, so `view`/`screenshot` there report the build flag rather than opening a window ([viewer](./viewer.md)).
- **General meshing algorithms as operations** — hex and hex-dominant meshing, boolean/CSG, boundary-layer inflation, quadrangulation and geodesic distance are each a library in their own right; the answer is an optional backend (the KaHIP pattern), not a native implementation.

---

## Suggested sequencing

Open work only; what shipped is in `CHANGELOG.md`.

1. **WASM parity ([§1](#_1-wasm-parity))** — the loader and dtype defects and the parity guard first; the JS-surface items next (no core change, no ABI bump); its core format gaps one at a time, each against its probe.
2. **Correctness debts ([§2](#_2-correctness-debts))** — small, independent, and any of them can be picked up between larger items; the MED and gmsh ordering pair first, since a mis-oriented cell is invisible until a solver rejects it.
3. **Sanitizer leg, then fuzzing ([§3](#_3-quality-of-implementation))** — a parallel track from day one; it does not compete for the same attention as features.
4. **Primitive constructors ([§5](#_5-operations))** — a few days, and a prerequisite of the conformance matrix, every demo surface and `extrude`/`revolve`.
5. **Core parity ([§4](#_4-core-parity-across-surfaces))** — MDPA first, then Side-region survival and sets → regions, then the rest by consumer demand.
6. **Format reach ([§6](#_6-format-reach))** — `.pvd`/`.pvtu` together with §1's `.vtm`, VTKHDF after them; `compute_normals` before glTF.
7. **Registration ([§7](#_7-ecosystem-reach))** — calendar-bound, so start the submissions early and let them run alongside everything else.
8. **Long-run spikes ([§8](#_8-long-run-spike-first))** — the benchmark tier decides the scale items; the NURBS spike is scheduled independently of the rest.

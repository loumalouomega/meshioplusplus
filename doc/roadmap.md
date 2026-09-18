# meshio++ roadmap

Status at time of writing: **v12.0.0** — 46 core formats plus four Python-only physics-ML ones, thirty-eight mesh operations + five data operations, six language surfaces (Python / C / Fortran / Julia / R / WASM), two viewers plus a browser dataset manager, a Blender add-on, a ParaView plugin, an MCP server, a settings-driven pipeline engine, a dataset-manifest layer with a PhysicsNeMo adapter, and a versioned ABI (`MESHIOPLUSPLUS_ABI_VERSION` 13).

This document lists what is *not* built. Nothing here duplicates shipped functionality; where a feature partially exists, the shipped half is named and the gap is stated explicitly. Release history lives in [`CHANGELOG.md`](https://github.com/loumalouomega/meshioplusplus/blob/main/CHANGELOG.md), not here.

Effort key: **S** = days, **M** = a couple of weeks, **L** = a month or more, **XL** = a project in its own right. An item names a **probe** — the failing test that proves the gap — wherever one is cheap.

## How this file works

- **Sections are ordered, and the order is the recommendation.** Each section states what belongs in it; an item that sounds exciting does not move up for that reason. An empty section is removed, not kept as a placeholder.
- **A closed item is removed, not struck through.** Its history is the `CHANGELOG.md` entry and the feature's own `doc/` page; a partly closed item is narrowed to what remains (the `AGENTS.md` change checklist rule).
- **Defect-shaped items go to [§1](#_1-correctness-debts) regardless of size** — behaviour that loses data, mis-orients cells, does not terminate on valid input or fails silently is not a feature request, even when the fix and the feature are the same work.
- **An item estimated from a doc, a `.d.ts` or a changelog alone says so** ("verify first") and names its probe; the code has repeatedly been more or less capable than its description.
- **[Non-goals](#non-goals-and-decisions-taken) record decisions already taken**, with their reasons, so they are not re-proposed as gaps.

## The map

![The roadmap at a glance: open items grouped by section, shaded by effort, with dependency arrows and the items that need a design pass or a research spike first](/diagrams/roadmap_map.svg)

---

## 1. Correctness debts

*Admission: behaviour that loses data, mis-orients cells, does not terminate on valid input, or drops something without saying so. First regardless of size, and each is independently shippable.*

- **MED quadratic 3-D node ordering is not converted.** `tetra10`, `hexahedron20`, `pyramid13` and `wedge15` are read and written with a warning and no meshio++↔MED permutation in either engine (`src/cpp/src/formats/med.cpp`, `src/python/meshioplusplus/med/_med.py`), so they may be mis-oriented for Salome, code_aster and code_saturne — the tools the format exists for. The linear 3-D permutations beside them are the template. *Probe:* a Salome-written `tetra10`/`hexahedron20` fixture whose cell volumes are all positive after a read. **S–M**
- **Gmsh `pyramid14` (and `wedge18`) node ordering is not permuted.** Both engines' reorder tables stop at `pyramid13`/`wedge15`, while [cell types](./cell_types.md) documents `pyramid14` as `pyramid13` plus a base-centre node and gmsh's own type 14 is edge-lexicographic. *Probe:* a gmsh-generated second-order pyramid whose base-centre node lies on the base plane after a read. **S, verify first**
- **Writers that drop data without a warning.** The C++ AVS-UCD writer keeps only the first integer cell-data array, silently (the Python writer warns); the SVG and TikZ writers skip unsupported cell blocks with no warning at all ([SVG](./formats/svg.md), [TikZ](./formats/tikz.md)). Every other writer follows warn-and-skip, and provenance's conversion assumptions sweep those warnings. *Probe:* each writer under `pytest.warns`. **S**
- **The C++ MED writer emits no MED 4.1 field bitmask attributes** (`LEN`/`LGC`/`LNA`/…), which the Python writer does; every flat binding and the native CLI use the C++ writer with no Python to defer to. The source calls this an interoperability gap with Salome/MEDCoupling, while [MED](./formats/med.md) calls the attributes required for medfile/mdump — establish which is true with `mdump` before sizing the fix. *Probe:* a C++-written, field-carrying `.med` opened by medfile. **S–M, verify first**
- **Number parsing follows the process locale.** 31 C++ readers parse with `strtod`/`strtoll` — chosen over `std::from_chars` because its floating-point overload is an Emscripten/libc++ hazard (`gid_read.cpp`) — and no first-party code pins `LC_NUMERIC`, so a host that adopts the environment's locale — Qt applications call `setlocale(LC_ALL, "")` at startup, as does any Python code that runs `locale.setlocale(locale.LC_ALL, "")` — may read `1.5` as `1` under a comma-decimal locale. A `detail/fast_number.hpp` that uses `from_chars`/`to_chars` where `<version>` reports `__cpp_lib_to_chars` and a C-locale `strtod_l` otherwise — the feature-test pattern `detail/format_compat.hpp` already uses — fixes it and is the base of the [§3](#_3-performance) text-I/O items. *Probe:* read an ASCII fixture from Python after `locale.setlocale(locale.LC_ALL, "de_DE.UTF-8")`. **S–M, verify first**
- **Reader fallback prints to stdout and swallows the reason.** `read()` `print`s the error for every ambiguous-extension candidate that fails (a `.msh` tries `ansys`, `gmsh` and `freefem` in turn, `_helpers.py`), which pollutes a CLI pipeline's output, and over forty format packages wrap the C++ reader in `except Exception: pass`, so nobody can see why the fast path declined. Route both through `warnings`/logging and let only a recognised "not handled here" exception fall back. **S**

---

## 2. Quality of implementation

*Admission: work that makes every other item safer to land. None of it is a feature, so none of it competes for the same attention — it can run in parallel with everything.*

- **A sanitizer CI leg** — ASan and UBSan over the existing `cpp-tests` job. No workflow passes `-fsanitize` today, and a fuzzer that finds a crash without one reports a symptom rather than the out-of-range read behind it. The precondition for the next item. **S**
- **Fuzzing the readers** (libFuzzer, then OSS-Fuzz if the project is accepted). 43 mostly hand-rolled parsers are reachable from a C ABI, a browser, a VS Code extension and an MCP server — untrusted input reaches them by design. One fuzz target per `registry_readers()` entry, seeded from `tests/python/meshes/`. The highest-value non-feature item in this document. **M**
- **A format conformance matrix** — one canonical mesh written to and read back from every writable format, asserting per format what survives (points, each cell type, point/cell/field data and their dtypes, each region kind) against a declared expectation. `tests/python/test_region_roundtrip.py` already does this for regions over Gmsh/Abaqus/MED, and `tests/cpp/test_sequence.cpp`'s `WriteSupportsTimeAgreesWithReality` is the registry-iterating shape to generalise it to. The declared expectations become a lossiness column in the [format table](./formats.md), which today has only Read/Write/dependencies, with lossiness scattered across its notes and fifty-five per-format quirks sections. The canonical mesh should be a primitive constructor from [§5](#_5-operations). **M**
- **Property-based testing** (Hypothesis) over the invariants the docs already articulate: partition-of-unity, volume conservation, conformity, byte-identical determinism, map composition. **M**
- **A benchmark harness that covers what ships, with a CI leg.** The suite exists (`benchmark/`, up to ~1M synthetic tets in Python, 257k in the C++ backend benchmark, which is off by default) but no CI job runs any of it, so a performance regression is found by a user. It is also narrow: `benchmark/bench.py` times 6 format labels of the 43 the core reads, and [benchmarks](./benchmarks.md) has no numbers for any operation or for the parallel backends — `src/cpp/benchmark/bench_backends.cpp` compares mesh backends only. Widen `bench.py` to every registry format; add a `bench_ops.cpp` (`extract_surface`, `smooth`, `refine`, `merge`, `clean`, `compute_sdf`, `decimate`, `partition`, `reorder`) over a size sweep and SEQ/OpenMP/TBB/Kokkos; add a 10M+ cell tier; run it on a schedule that records rather than gates. Every [§3](#_3-performance) item is gated on this showing its before/after, and it decides whether the scale items in [§8](#_8-long-run-spike-first) matter at all. **S–M**
- **Test and install the ParaView plugin.** `tools/paraview-meshioplusplus-plugin.py` ships as a reader and writer, but nothing tests it, and the `data_files` entry that would install it is commented out in `pyproject.toml`, so [its page](./paraview_plugin.md) describes a plugin path nothing writes. A `pvpython` smoke step plus the install fix. **S**

---

## 3. Performance

*Admission: a measured or code-verified slowdown in a path a user hits, with the shape of the fix named. Nothing here is scheduled before the [§2](#_2-quality-of-implementation) harness can show its before/after.*

Two findings frame the section. First, **the serial phases below are deliberate**: each is documented in the code as a determinism pin, not an oversight — output is byte-identical across parallel backends and thread counts, and the reference-file tests enforce it — so every fix must keep that guarantee and prove it with a SEQ-versus-OpenMP diff, not assert it. Second, **every parallel item is conditional on the backend**: a SEQ build (and the `stl` fallback without TBB) runs `parallel_for` sequentially, so each change must also show that SEQ does not get slower — a parallel sort is O(n log n) where the hash map it replaces is O(n).

**Measured regressions** (verify first: `benchmark/results.csv` is one run on one machine). Gmsh binary write at 0.68× legacy meshio and MED read at 0.81–0.93×, both reported as parity or better in [benchmarks](./benchmarks.md) until this change; XDMF (HDF5) write at 0.92–0.97×; XDMF read at 1.0× on a single-block mesh against 10× on mixed topology, which suggests the post-read conversions — `xdmf.cpp` has no `parallel_for` at all. Re-measure on the widened harness, then profile. **S each to size**

**Text I/O.**

- **A shared tokenizer and number path.** 29 readers split each line with `std::istringstream` and `>>` into a `std::vector<std::string>` — one stream and N heap allocations per line (`su2.cpp`, `vtk_read.cpp`, `mdpa.cpp`, `avsucd.cpp`, `tecplot.cpp`, `unv.cpp`, `flac3d.cpp` and 22 more) — while `gmsh.cpp`'s `GmshCursor`, a `string_view` cursor calling `strtod` straight on the buffer, is the in-repo model to copy. Built on the `fast_number.hpp` from [§1](#_1-correctness-debts) and migrated reader by reader against the reference files. The worst cases: XDMF ASCII `DataItem`s (a `std::string` and a dtype switch per scalar), VTU ASCII arrays (`push_back` with no `reserve`, then a second pass with a dtype switch per element) and `mdpa.cpp`, which materializes one `std::string` per line of the whole file before parsing. **M**
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
- **A declined C++ read costs a full parse before the Python one starts.** `vtu_read.cpp` builds the whole XML DOM, base64 bodies included, before rejecting lzma, appended data or multiple pieces; `gmsh.cpp` rejects `$Periodic` only after `$Nodes` and `$Elements` are parsed; with an ambiguous extension one `.msh` can be parsed up to six times. A cheap pre-flight — attribute and section-header scans — before the expensive parse, with only early rejections falling back. Pairs with the [§1](#_1-correctness-debts) swallowing item; the durable fix is [§4](#_4-core-parity-across-surfaces)'s core parity. **S–M**
- **The MCP server re-reads the input file on every tool call** — 65 call sites in `mcp/_tools.py` go through an uncached `_load()`, so an agent's info → clean → decimate → convert parses one file four times. A bounded cache keyed on (path, `mtime_ns`, size); anything weaker manufactures a stale-read bug. **S**
- **pybind11 per-call overheads**, together: the `Mesh`→numpy conversion re-imports the `Mesh` class on every call, the contiguity check does a Python attribute lookup per array, operations clone connectivity they never change (`transform`), and polygon/polyhedron blocks cross the boundary one node id and one face at a time where the WASM binding already uses a CSR triple. **S–M**
- **Flat-binding accessor copies**, together: R copies the points twice and shifts connectivity to 1-based with a scalar loop, and Julia's safe accessors `copy` the borrowed view; both are documented, and both are fixable behind the same accessor names. **S**
- **The browser viewer round-trips every operation result through a VTP file** — written into MEMFS, copied out with `.slice()`, then parsed again by vtk.js on the main thread. A typed-array mesh channel between worker and renderer removes three full passes. **L**

**Deliberately not**, and recorded so it is not re-proposed:

- **Explicit SIMD intrinsics or `-march` flags** — portability across wheels, WASM and the release binaries is worth more than the scalar kernels cost; revisit only if the harness shows a kernel dominating.
- **Kokkos device execution** — every `parallel_for` body captures host pointers (`parallel.hpp`); the GPU route is the DLPack/CuPy handoff ([GPU handoff](./gpu.md)).
- **A BVH in place of the uniform grid** — the tiebreak argument above.
- **Tuning the pure-Python fallback readers** — the fix is making the C++ path accept the file ([§4](#_4-core-parity-across-surfaces)), not a faster fallback; `_decimate.py`'s heap-based twin is deleted once `decimate` accepts its inputs, not optimised.

*Recommended posture:* the harness and the measured regressions first; then the small isolated wins — the polyhedron rescan in [§1](#_1-correctness-debts), `b64decode`, the lazy CLI import, `optimize_volume` — then the shared facet table, the largest total win; the boundary items as their consumers ask.

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
- **Sets → regions, phase 2.** `ansysInp` and `unv` sets still travel in side-channel structs, XDMF `Sets` are not mapped, and VTU/VTP need a documented region convention ([Named regions](./regions.md)). **M**
- **Named Side regions surviving operations.** `subdivide`, `agglomerate`, `undo_green` and `convert_cells(simplexify)` drop them through their parent-cell remap, and the cutters (`slice`, `isosurface`, `extract_surface`, `extract_skin`) drop them outright ([Named regions](./regions.md)); for a Kratos model, Side regions are where the boundary conditions live. **M**
- **A structured pipeline report on the flat ABI.** C, Fortran, Julia and R receive status plus `mio_last_error()` only; a caller-buffer JSON accessor is recorded as a follow-up, as are the v2 multi-mesh steps (`Inputs:` for `Merge`/`Interpolate`/`UndoGreen`, an `Output.Pattern` for `Split`/partition) ([pipelines](./pipeline.md)). **S–M**
- **Point/cell sets beyond regions in the core**, so the `convert -s/-d` sets↔data conversions work in the native CLI and flat bindings ([Julia](./julia.md)). **S–M**

---

## 5. Operations

*Admission: a new operation, or a public face for machinery that already exists privately inside one.*

**Generation.** Almost every operation transforms a mesh you already have. The exceptions all start from something else — `grid` from a lattice (`detail/grid_lattice.hpp`), `voxelize`/`compute_sdf` from a surface's bounding box, `remesh_volume` from a closed surface — and nothing builds a shape from parameters, sweeps one, or triangulates a domain.

- **Primitive constructors** — `box`, `sphere`, `cylinder`, `disk`, in their own `operations/primitives.hpp` beside `grid`. Dependency-free, and it removes the fixture-file dependency from tests, docs, notebooks, the browser demo and the MCP server; it is also the canonical mesh the [§2](#_2-quality-of-implementation) conformance matrix needs. Highest leverage per line of code in this document. **S**
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
- **`.pvd` and `.pvtu`/`.pvtp`** — the on-disk ParaView face of the sequence engine (a `.pvd` is a time-indexed collection) and of `partition` (a `.pvtu` is one piece per part plus an index). Shares the index-plus-pieces machinery the shipped [`.vtm` writer](./formats/vtm.md) already has (v11.6.0). **S–M**
- **Point-cloud formats** `.xyz` and `.pcd` — `select_points`, `subsample_points` and `proximity_graph` form a point-cloud path with no point-cloud file at either end of it. **S**
- **LS-DYNA keyword input `.k`** — the one major solver-input keyword format missing beside Abaqus, Nastran and ANSYS; the same ASCII reader shape, with `*PART` as regions. **M**
- **CalculiX `.frd` results (read)** — CalculiX takes Abaqus-style `.inp` input, which meshio++ already writes; reading its results file closes the loop for an open-source solver. **S–M**
- **glTF / `.glb` (write)** — a web-native skin export for the browser viewer, dashboards and anything downstream of the Blender path; needs `compute_normals` first. **S–M**
- **Half-gaps of shipped formats** — OpenFOAM *binary* write (binary read exists) **S**; SU2 multizone (`NZONE`) **S**; Tecplot's non-transient multiple zones **S**; EnSight variable *write* (reading already shipped) **S**.

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
| ParaView | The plugin exists; installing it with the wheel is the [§2](#_2-quality-of-implementation) item. |

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

**Scale.** Memory-mapped reads ([mmap](./mmap.md), via `ReadOptions` and the C ABI but not Python's `read()`) roughly halve the peak footprint of a large read, and the XDMF series appender and Python's chunked `write_dataset` write a series or a dataset without holding it — but nothing writes one mesh larger than memory, and no operation streams. Run the [§2](#_2-quality-of-implementation) benchmark tier first; it decides whether either item below matters.

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
- **The Python-only layers stay Python** — `pmsh`, `zarr`, `cae` and `usd`, and the physics-ML surface (`tessellate`, grids, point budgets, proximity graphs, datasets, training) are export targets and tooling for a training pipeline, registered in Python rather than the shared C++ registry. A core kernel they call (the neighbour search in [§3](#_3-performance)) does not change that.
- **Polyscope in the release CLI binaries** — excluded deliberately, so `view`/`screenshot` there report the build flag rather than opening a window ([viewer](./viewer.md)).
- **General meshing algorithms as operations** — hex and hex-dominant meshing, boolean/CSG, boundary-layer inflation, quadrangulation and geodesic distance are each a library in their own right; the answer is an optional backend (the KaHIP pattern), not a native implementation.
- **KaHIP in the WASM build** — no Emscripten port, and a graph partitioner would bloat every consumer's bundle; `"auto"` resolving to the SFC method is the answer, and `"kahip"` throwing by name is the contract.
- **zstd/lz4 codecs, memory mapping and the Kokkos backend under Emscripten** — the codecs have no Emscripten port and are compiled out, CMake refuses the Kokkos backend under Emscripten, and there is nothing to map inside MEMFS.
- **Single-file output for DOLFIN, TetGen and EnSight from WASM**, or lifting DOLFIN's simplicial restriction — both are facts of the formats; `writeMesh`'s written-paths return value (v11.2.0) is the fix for the bookkeeping they cause.
- **Polyhedron blocks as a WASM gap** — `vtu`, `ensight`, `cgns`, `med` and `openfoam` all write them from WASM; a consumer that cannot is constrained by its own data model.

---

## Suggested sequencing

Open work only; what shipped is in `CHANGELOG.md`.

1. **Correctness debts ([§1](#_1-correctness-debts))** — small, independent, and any of them can be picked up between larger items; the MED and gmsh ordering pair first, since a mis-oriented cell is invisible until a solver rejects it.
2. **Sanitizer leg, then fuzzing ([§2](#_2-quality-of-implementation))** — a parallel track from day one; it does not compete for the same attention as features.
3. **Performance ([§3](#_3-performance))** — the harness and the measured regressions first, then the isolated wins (the polyhedron rescan, `b64decode`, the lazy CLI import, `optimize_volume`), then the shared facet table.
4. **Primitive constructors ([§5](#_5-operations))** — a few days, and a prerequisite of the conformance matrix, every demo surface and `extrude`/`revolve`.
5. **Core parity ([§4](#_4-core-parity-across-surfaces))** — MDPA first, then Side-region survival and sets → regions, then the rest by consumer demand.
6. **Format reach ([§6](#_6-format-reach))** — `.pvd`/`.pvtu`, VTKHDF after them; `compute_normals` before glTF.
7. **Registration ([§7](#_7-ecosystem-reach))** — calendar-bound, so start the submissions early and let them run alongside everything else.
8. **Long-run spikes ([§8](#_8-long-run-spike-first))** — the benchmark tier decides the scale items; the NURBS spike is scheduled independently of the rest.

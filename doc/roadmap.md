# meshio++ roadmap

Status at time of writing: **v16.14.0** — 76 core formats plus four Python-only physics-ML ones, thirty-nine mesh operations + six data operations, six language surfaces (Python / C / Fortran / Julia / R / WASM), two viewers plus a browser dataset manager, a Blender add-on, a ParaView plugin, an MCP server, a settings-driven pipeline engine, a dataset-manifest layer with a PhysicsNeMo adapter, and a versioned ABI (`MESHIOPLUSPLUS_ABI_VERSION` 17), sanitizer and fuzzing gates, and a format conformance matrix.

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

Sizes: **S** ≈ days, **M** ≈ two weeks, **L** ≈ a month or more, including tests, docs page and CLI wiring. Sizes assume the shared infrastructure in §1.1 exists; the first format that needs a component pays for it.

Link legend: unmarked links were opened or returned by a search while this section was written (September 2026); links marked † are quoted from memory and must be checked before relying on them. Vendor documentation portals move often — when a link is dead, search the document title.

---

Every format this section queued is now read: the native FEM interchange, solver input and result-file readers of v16.5.0 to v16.12.0 (Tecplot, Femap, MFEM, Patran, Abaqus `.fil`, Radioss, Z88, libMesh, Marc, Ansys `.rst`, LS-DYNA and OP2), DOLFINx's ADIOS2 `.bp` ([`vtx`](./formats/vtx.md)) and Tecplot SZL through TecIO ([`szplt`](./formats/szplt.md)) of v16.13.0, and the [vendor-runtime routes](./vendor_routes.md) for the files only a vendor's software reads (Abaqus `.odb`, Marc `.t16`, Femap `.modfem`, Ansys results beyond the native reader). What they still lack is a check against a run of the vendor tool, listed under [Awaiting a licensed run](#awaiting-a-licensed-run), and the shared components below.

### 1.1 Shared infrastructure worth building once

- **Fixed-width card tokenizer, remaining adoption.** `detail/keyword_card.hpp` (v15.2.0, built for LS-DYNA) splits a card into standard, long, I10 and free layouts per line, parses real fields with Fortran spellings and, since v16.3.0, parses a Fortran format line such as `(3i9,6e21.13e3)` into the fields `split_fixed` cuts (the ANSYS `.cdb` reader's). The Nastran/OptiStruct (v16.1.0) and Patran (v16.5.0) readers use it too, and the Radioss (v16.7.0) and Marc (v16.8.0) readers its field converters (their column widths follow the deck: Radioss's input version, Marc's `EXTENDED`); UNV should adopt it instead of growing its own (the `.frd` reader of v15.3.0 slices its own columns: the layout is fixed by the record key, not by a format line).
- **Node-ordering permutation registry, remaining half.** `detail/node_order.hpp` and `_node_order.py` (v16.0.0, built for Code_Aster `.mail`) hold the MED, Code_Aster, `.frd`, UNV, COMSOL, Elmer, FEBio, FLUX, Patran, libMesh, Radioss and Z88 tables in both directions, with a self-test that maps a reference element through every table and checks midpoints and a positive Jacobian ([node ordering](./node_ordering.md)). What remains: the gmsh, CGNS, GiD, Exodus and Kratos tables still live inside their readers and should move in when those formats are next touched (Marc, v16.8.0, needs no entry: every Marc element type read is in meshio++'s node order but the 3-node rebar lines, whose swap lives in the reader; ANSYS's degenerate shapes and Femap's 20-slot brick layout are slot maps, not permutations, and live in their readers; MFEM needs none, since its order-2 nodes are placed by the entity they sit on).
- **Fortran unformatted-record reader, remaining adoption.** `detail/fortran_records.hpp` (v16.7.0, built for Abaqus `.fil`) sniffs 4/8-byte markers in either byte order and splits a file into its records, and `detail/binary_stream.hpp` reads values in a chosen byte order (libMesh `.xdr`). The OP2 reader (v16.9.0) builds on it, and the Python twins share `_fortran_records.py`; EnSight Fortran binary should too. The ANSYS `.rst` reader addresses records by pointer and keeps its own view, and `d3plot` (v16.9.0) is word-addressed, not record-framed.
- **Directory-as-format support, remaining half.** Reading sniffs a directory by its files since v16.2.0 (Elmer, OpenFOAM, and a DOLFINx `.bp` since v16.13.0), and sequence globs keep suffixed directories (`.pmsh`, `.zarr`, and `.bp` in both engines). What remains: an Elmer directory has no suffix, so a glob cannot find it and it must be listed explicitly.

### 1.2 Reader and writer parity

*Admission: a format the [format table](./formats.md) lists as read-only, or as writing less than it reads.* The table is the probe: every row whose Read and Write columns differ is a candidate, and today those are the rows below. A result file has nothing to write back (the solver is the only producer, and `frd`'s page says so), so parity is not automatically the goal; each row needs a decision, written into its `doc/formats/<format>.md` page, before any code.

- **Triage first.** Sort every row into one of three outcomes: *write it* (a file a downstream tool reads), *stay read-only by design* (say why, and list it under [Non-goals](#non-goals-and-decisions-taken)), or *narrow the gap* (a writer that keeps less than the reader reads). The candidates by kind: **input decks and meshes** with no writer — `radioss` (`.rad` starter deck) and `marc` (`.dat` deck); **partial writers** — `femap` (no results) and `z88` (no results, `side` regions or other data); **result files** — `abaqus_fil`, `ansys_rst`, `ansys_rst_cyclic`, `frd`, `lsdyna_binout`, `lsdyna_d3plot`, `marc_t19`, `nastran_h5`, `nastran_op2`, `radioss_anim`, `radioss_th`, `xplt`, plus `vtx` and `szplt`, whose only reader is a heavy or vendor library; and the **export-only** `gltf`, `svg` and `tikz`, which are a surface, not a mesh, and have no reader to pair. **S**
- **Writers for the decks.** A `radioss` starter-deck writer and a `marc` deck writer, each in both engines, checked by reading the file back and by the vendor's own reader where one is at hand (OpenRadioss's starter is open; Marc is under [Awaiting a licensed run](#awaiting-a-licensed-run)). Node-order tables already exist in both directions ([node ordering](./node_ordering.md)); the new work is the card layouts, regions as sets, and a declaration in `tests/python/conformance_spec.py`. **M each**
- **Narrow the partial writers.** `femap` writes the mesh and groups but not the results it reads (output vectors), and `z88` writes the structure file, its input files and Z88Aurora's sets but drops `side` regions, other data arrays and the results it reads ([Femap](./formats/femap.md), [Z88](./formats/z88.md)). Each is either closed or recorded on its page as the deliberate limit. **S–M each**
- **Result-file writers, only with a consumer.** A synthetic `.frd`, `.op2` or d3plot is worth writing when a downstream tool reads it (a post-processor fed from a surrogate's output is the case to look for), not to fill a cell of the table. Until then they are recorded as read-only by design. **M–L each**
- **Guard against drift.** Extend the conformance matrix, or add a test over the registry, so a new format that reads without writing must carry an explicit `read-only` reason, and the table's reader/writer columns are checked against the registry rather than kept by hand. Probe: register a reader with no writer and no exemption; the test must fail. **S**
- **Done when.** Every row where Read and Write differ either is closed or names its reason on its format page, the table and the registry agree by test, and the [C](./c_api.md), [WASM](./wasm.md), [Fortran](./fortran.md), [Julia](./julia.md) and [R](./r.md) docs list the new writers.

### 1.3 Suggested order

1. **Reader and writer parity ([§1.2](#_1-2-reader-and-writer-parity)):** the triage and the drift guard first, then the deck writers.
2. **Checks against the vendor tools:** each item under [Awaiting a licensed run](#awaiting-a-licensed-run) as soon as a licence, or a file the tool wrote, is in hand.

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

- Every link marked †. The Femap `feFileWriteNeutral` argument list and the Tecplot macro's `$!READDATASET` form in `contrib/` were taken from the vendor references, not from a run.
- Ansys `.rst` records flagged both bit-sparse and windowed-sparse (three files of Ansys' DPF example data): their layout is neither form's, and pymapdl-reader misreads them too.
- MFEM NURBS meshes with non-conforming patches (`patch_topology` with hanging knots) are refused; no user has asked.

#### Awaiting a licensed run

Checks that need the vendor tool itself, or a file it wrote that no public source has. Each reader works from the vendor's documentation and every public file found (September 2026); what remains is one probe per line.

The [vendor-runtime routes](./vendor_routes.md) of v16.13.0 are tested against stand-ins of the vendor APIs only. Each needs one recorded run, with the vendor version written into its page; what to install for it:

- **Abaqus `.odb`.** Abaqus 2020 or later (its `abaqus python`, whose bundled numpy the script needs; 2023 and before run Python 2.7, 2024 and later Python 3). Run `contrib/abaqus_odb/odb_to_vtu.py` on a job with shells, C3D20R and an assembly-level set, and compare its nodal `U` and centroidal `S` against the same job's `.fil` read natively ([route](./routes/abaqus_odb.md)).
- **Marc `.t16`.** Marc/Mentat 2020 or later with PyPost (`py_post`, run with Mentat's Python). Run `contrib/marc_t16/t16_to_vtu.py` on a job that writes both `.t16` and `.t19`, and compare against the `.t19` read natively; confirm that post file position 0 holds the model ([route](./routes/marc_t16.md)).
- **Ansys results via PyDPF.** A DPF Server (Ansys 2021 R1 or later, or the standalone DPF Server 2024 R2 or later) and `pip install ansys-dpf-core pyvista meshioplusplus[h5py]`. Run `contrib/ansys_dpf/rst_to_vtkhdf.py` on a structural `.rst` and compare against the native `.rst` reader; then on a `/FCOMP,RST,1` file, which the native reader refuses ([route](./routes/ansys_dpf.md)).
- **Femap `.modfem`.** Femap 2301 or later on Windows, and `pip install pywin32`. Run `contrib/femap/export_neutral.py`, check `feFileWriteNeutral`'s argument list against that version's API help, and read the neutral file back ([route](./routes/femap_modfem.md)).
- **Tecplot `.szplt`.** Tecplot 360 2023 or later, and `pip install pytecplot`. Save a data set as `.szplt` and as `.plt` from Tecplot itself and read both (the TecIO-backed reader against the native one); run `contrib/tecplot/szplt_to_plt.py` and the macro `contrib/tecplot/szplt_to_plt.mcr`; and point `MESHIOPLUSPLUS_TECIO_LIBRARY` at Tecplot's own shared `libtecio` for the ctypes reader ([route](./routes/tecplot_szplt.md)).

The readers:

- **Tecplot.** A `.plt` written by `preplot` from the ASCII fixtures (`#!TDV112`), read identically to its `.dat`; binary files of versions 100–105 and 109–111, of which no sample was found (the VisIt suite covers 71, 75, 106–108).
- **Femap.** Element (404) and group (408) records written by Femap 10 and by 12 or later (the public 2401 files hold none); the slot layout of the 13-node pyramid (topology 19); a file meshio++ writes, imported into Femap.
- **Patran.** A neutral file exported by ANSA or HyperMesh.
- **Abaqus `.fil`.** A binary run with shells and C3D20R, checked against its `.dat` printout (the truss run TenBarArea is checked).
- **Radioss.** Animation files of the layouts before OpenRadioss (magic numbers `0x5426` to `0x542B`): neither OpenRadioss's tools nor its readers document them.
- **Marc.** A run of Marc itself; post files before revision 9 (Marc K7 and earlier, Volume D 2000's PLDUMP layout); Volume A's edge and face numbering, to map edge and face sets (kept as `(cell, number)` field data) onto facets.
- **Ansys `.rst`.** A file written with `/FCOMP,RST,1` or higher (zlib records): `detail/zlib_inflate` is ready, the record framing needs a sample.
- **LS-DYNA.** A family from LS-DYNA with adaptive remeshing (`NADAPT`), the `d3thdt` time-history database, the cubic and bubble solids (21P, 15T, 20T, 40P, 64), airbag particles and rigid roads (read from the manual and lasso-python only), a real `d3part`, and the whole reader against LS-PrePost; whether LS-DYNA writes a ten-node solid's extra nodes once or twice.
- **Nastran OP2.** Random element tables, complex SORT2 nodal tables, CBUSH and CGAP forces, CBUSH real stresses, and NX's solid types 300–303 (read from pyNastran's layouts; no sample).

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

The sanitizer leg, the reader fuzzing, the format conformance matrix, the property tests, the benchmark harness, the ParaView plugin test and the fallback narrowing all shipped in v16.14.0 ([fuzzing and sanitizers](./fuzzing.md), [format conformance](./conformance.md), [benchmarks](./benchmarks.md)). What remains:

- **OSS-Fuzz submission.** The harness already takes its format from its binary's name, and `tools/fuzz/oss-fuzz/` holds a draft `project.yaml` and `build.sh`; what is left is the submission to `google/oss-fuzz` and the maintainer contact it needs. Calendar-bound, like the registries in [§7](#_7-ecosystem-reach). **S**
- **Fuzz the library-backed readers.** The campaign skips the formats whose bytes go to HDF5, netCDF, ADIOS2 or TecIO (`tools/fuzz/not_fuzzed.txt`): those libraries are not instrumented in the build, and their own parsers are not this project's to fix. meshio++'s code that walks the objects they return (dataset shapes, attribute types, link targets) is reachable all the same; fuzzing it means building the libraries with the sanitizers, or fuzzing at the object level behind a stub. Probe: an HDF5 file whose connectivity dataset's shape disagrees with its declared element count. **M**
- **A benchmark trend store.** The `benchmark` workflow records CSVs as artifacts and never gates; comparing runs means downloading them. A store that keeps one row set per run (a branch or a release asset) and a page that plots it would make a slow drift visible without anyone looking for it. **S**

---

## 4. Performance

*Admission: a measured or code-verified slowdown in a path a user hits, with the shape of the fix named. Every item shows its before/after on the benchmark harness (`benchmark/bench.py` for I/O, `tools/bench_ops.sh` for operations; [benchmarks](./benchmarks.md)).*

Two findings frame the section. First, **the serial phases below are deliberate**: each is documented in the code as a determinism pin, not an oversight — output is byte-identical across parallel backends and thread counts, which repeated-run tests and the C++-versus-numpy byte comparisons enforce — so every fix must keep that guarantee and prove it with a SEQ-versus-OpenMP diff, not assert it. Second, **every parallel item is conditional on the backend**: a SEQ build (and the `stl` fallback without TBB) runs `parallel_for` sequentially, so each change must also show that SEQ does not get slower — a parallel sort is O(n log n) where the hash map it replaces is O(n). SEQ is not a corner case: the Linux wheels and the native CLI release binaries are built SEQ, and WASM ships a SEQ build beside the threaded one, so an algorithmic win reaches every user while a parallel one reaches source, conda and threaded-WASM builds.

**Landing an item.** Each item names only what is particular to it. All of them:

1. Record the before on SEQ, OpenMP and TBB at 1, 4 and 8 threads: `tools/bench_ops.sh before.csv "SEQ OPENMP TBB" "1 4 8" -- --hash --ops <op> --tier M --tier L`, or `benchmark/bench.py --sizes M,L --formats <format>`.
2. Show the output is byte-identical across backends and thread counts, and to the digests before the change: `meshioplusplus_bench_ops --hash`, with `BASELINE=before.csv` for the second sweep ([benchmarks](./benchmarks.md#determinism-check)). The rows each item needs exist since v16.15.0; an item whose operation has none adds one first.
3. Show SEQ is not slower.
4. Classify the change by the [ABI policy](./abi.md): the body of an exported, non-inline function is free; an inline or template body in an installed header is Tier B; a new function in an installed header is additive and goes in the [ABI review](./abi_reviews.md); a changed signature, or the layout of an installed type, is Tier A. Regenerate the single header.
5. Correct every code comment that describes the old algorithm (the items name the stale ones already found), and update the numbers in [benchmarks](./benchmarks.md).

**Text I/O.**

- **The tokenizer's remaining readers.** Since v16.20.0 the core-private `detail/text_cursor.hpp` splits a mapped file into `string_view` lines and blank-separated tokens, parses bounded tokens (`parse_double_token`, `parse_int_token`, and `strtoll_token`/`strtoull_token`/`parse_double_prefix` for the lenient C-library semantics), and offers `TextStream`, a drop-in for a classic-locale `istringstream >>` pinned against it token for token by `test_text_cursor.cpp`. `mdpa`, `su2`, `avsucd` and `tecplot` parse views of the mapped file; XDMF's `DataItem` and VTU's ASCII arrays parse straight into the typed array; and the streams of `obj`, `stl`, `ply`, `gmsh`, `openfoam`, `flac3d`, `netgen`, `permas`, `triangle`, `tetgen`, `freefem`, `flux`, `abaqus`, `dex`, `ip`, `wkt`, `vti`, `vtr` and `vts` are `TextStream`s. What remains is performance only: every stream is pinned to the classic locale. **S per reader**
  - **Streams.** `ensight` (six), `vtk_read`, `off`, `mfm`, `patran`, `radioss`, `ansys` and one in `tecplot` still build a `std::istringstream`.
  - **Allocating tokenizers.** `gid_split`, `pcd.cpp`, `xyz.cpp`, `marc.cpp`, `lsdyna.cpp`, `femap.cpp`, `radioss.cpp`, `ansysinp.cpp`'s `ans_commas` and `keyword_card`'s `split_card`/`split_fixed` (which six readers share) return a `std::string` per token, as do `medit.cpp`'s `Tokenizer::next` and `unv_real`. `abaqus`, `ansysinp`, `dex`, `ip`, `flux`, `marc`, `netgen`, `permas`, `nastran`, `xyz` and `lsdyna` still hold the whole file as a `std::vector<std::string>` of lines; `split_lines` over a `FileSource` replaces it.
  - **Fold in.** Eight private cursors predate it, each in an anonymous namespace: `GmshCursor`, `VtkCursor` (`vtk_read.cpp`), `EnsightAsciiCursor`, `UnvLineReader`, `MfLexer` (`mfem.cpp`), `MdpaCursor`, `GidLineCursor` and `FnCursor` (`femap.cpp`).
  - **Gate.** Each reader's fixtures and its C++ and Python tests; `test_text_cursor.cpp`; `test_no_locale_sensitive_number_io.py`; `test_gmsh.py`'s byte comparison of the C++ and Python ASCII writers. Measure with `bench.py --sizes M,L --formats <format>`.
  - **Watch.** `GmshCursor::next_int` parses a double, so it accepts `2.0` and `1e3`; a strict integer path would reject files read today (`tests/python/meshes/gmsh/indented.msh`), so keep a lenient variant. The views end where the file ends: parse through the bounded token functions, never `parse_double` over a view, which reads up to the first character that cannot continue a number and relies on a terminator a mapped `FileSource` does not have.
- **Hoist the per-element dtype switch in the remaining formats.** Since v16.19.0 the VTK XML and Exodus loops this item named switch on the dtype once per array, through the private `detail/typed_view.hpp` (`Int64View`, `DoubleView`: an array as `read_int`/`read_double` read it, zero-copy when the dtype matches). About 390 per-element `read_double`/`read_int` calls remain in the formats, the most in `gmsh` (35), `cgns` (32), `med` (20), `libmesh`, `febio`, `gid`, `z88` and `pcd`; `write_double`/`write_int` dispatch per element too. **S each**
  - **Watch.** Each `dispatch_dtype` site instantiates its body for every dtype: watch the WASM artifact size `wasm.yml` reports; a view converts once instead. The output is bit-identical, since `read_double` is exactly a `static_cast<double>`.

**Memory and allocation.**

- **Polyhedral readers still building nested vectors.** Every backend stores ragged blocks as CSR since v16.16.0 and takes it whole through the CSR `AddPolygonBlock`/`AddPolyhedronBlock` overloads, which the VTK XML reconstruction, OpenFOAM, CGNS, MED, the C API and the WASM binding use; `ensight` (`nsided`/`nfaced`), `cgns_mll`, `tecplot`, `vtkhdf` and `ansys` still build one vector per cell and per face and let the nested overload convert, and `reorder` rebuilds ragged blocks the same way when it rewrites them (`compute_normals` builds CSR since v16.18.0). Move each onto the CSR overloads. **S each**
  - **Watch.** WASM's `Int32Array` offsets overflow past 2³¹ entries.

**Serial phases inside parallel operations.**

- **Waiting for the next ABI bump.** Since v16.18.0 every operation this section listed runs its serial phases in parallel or in gather form, byte-identical to before: the shared facet table (with `group_slots`' counting sort chunked), welding and duplicate cells on every block kind, the distance kernel's construction, marching, `split`, `hessian` and the cell-to-point average, `compute_normals`, `remesh_volume`, `remesh`'s setup passes, `agglomerate`, `undo_green`, a core neighbour search for `proximity_graph`, and the per-element dtype switch in nine operations. `test_op_goldens.cpp` pins 28 cases to digests taken from the previous implementation. What is left needs an installed layout to change: **M**
  - **Caching the facet table.** An N-step pipeline still rebuilds the facet table N times; caching it on the `Mesh` needs the `Mesh` layout, which `test_abi_layout.cpp` pins.
  - **The distance grid's map.** `build_distance_query` groups its (bucket, triangle) pairs in parallel, then moves each bucket into `SpatialGrid`'s `unordered_map` one at a time; a sorted bucket table would change `DistanceQuery`'s installed layout. A BVH is deliberately *not* the fix: `detail/surface_distance.hpp` rejects one because its visit order would make the equidistant-triangle tiebreak observable.
- **Small leftovers.** `partition`'s dual graph keeps its own map (KaHIP builds only); `gradient`'s and `diff`'s per-cell readers (`grad_load_cell`, which reads through the installed `cell_rings`, and `diff_cell_nodes`) still switch on the dtype per element; and writes (`write_double`, `write_int`) dispatch per value in `refine` and `convert_cells`. **S**
- **Serial by design**, so not re-proposed: the keep-first weld loop (chains make it order-dependent), the union sweeps of `split` and `compute_normals` (cheap, and their partitions are order-independent anyway), and `remesh`'s energy sweep (every accepted move changes what the next test reads; it has no numpy twin because a near-tie decided differently yields a different clustering).

**Boundaries and startup.**

- **The Python bindings never release the GIL.** Nothing in `bindings/python` uses `gil_scoped_release` or a `call_guard`, so a multi-gigabyte read or a 10M-cell operation blocks every other Python thread, and no caller can convert files in a thread pool. It matters most where the core is SEQ — the Linux wheels — since threads are then the only parallelism; `_sequence.py` runs a process pool today because the operations release no GIL. **S–M**
  - **Where.** There is no central wrapper: 252 `m.def` lambdas each convert inline (`py_to_mesh` 141 times, `mesh_to_py` 121). A write converts, then calls the core: release in between, and re-acquire before its `PyMeshRefs` die. A read returns `mesh_to_py(read_x(path))`, so each must be split — `Mesh m; { release; m = read_x(path); } return mesh_to_py(std::move(m));` — with `core_read_options`, which takes Python objects, evaluated first. A small `nogil` helper inside each lambda is the practical form; `py::call_guard<py::gil_scoped_release>` fits only lambdas that touch no Python object, such as `sniff_format`.
  - **Keep the GIL**, or a process-wide mutex, around the libraries that are not thread-safe: HDF5 (`hdf5_util`, CGNS, VTKHDF), netCDF (Exodus), the CGNS MLL and ADIOS2. The one callback into Python, `gid_write_series`'s step pull, re-acquires it.
  - **Other races** a released GIL exposes: `std::localtime` in `h5m.cpp` (use `localtime_r`/`localtime_s`), and a `getenv` per write in `provenance.cpp`, racing with writes to `os.environ`. Logging and provenance are already thread-safe (a mutex, and `thread_local` state on both sides).
  - **Watch.** Each OpenMP region runs `omp_get_max_threads()` threads, so N Python threads oversubscribe the machine N times: document it, or expose a thread cap. Another thread can still write into input arrays the core is viewing — a data race, not a use-after-free, since `PyMeshRefs` keeps them alive.
  - **Gate.** A new thread-pool test (results equal to serial calls), and a heartbeat thread counting ticks during one long read, which shows directly whether the GIL was released.
- **Lazy submodules for `import meshioplusplus`.** `python -X importtime` bills ~260 ms to `_cli` only because it is the first name in the package's `from . import (_cli, …)`, so it pays for numpy and everything else the package imports anyway. With `_cli` stubbed (Python 3.14, September 2026): numpy ~90 ms, the `_core` extension ~58 ms, the ~320 other package modules ~79 ms, and `_cli` itself ~7 ms; `rich`, `xml.sax.saxutils` (with `urllib`, `http.client` and `ssl`) and `importlib.metadata` left the import path in v16.15.0, and `test_import_footprint.py` keeps them out. **M**
  - **Lazy `_cli`** saves only ~12 ms, and changes behaviour: format registration follows import order, and without `_cli` first, `.msh` resolves to ansys, freefem, gmsh instead of ansys, gmsh, freefem, so every Gmsh file would try FreeFEM first. Pin the order explicitly before making it lazy (`test_helpers.py` documents today's). `_cli` is in `__all__`, so `from meshioplusplus import *` would trigger a module `__getattr__`.
  - **The real lever** is loading operation and format submodules on first use through a module `__getattr__`, for which the registration-order pin is also the precondition. The CLI's entry point imports the whole package regardless, so these gains reach library and MCP users (the latter once per process), not a single CLI call.
  - **Gate.** Extend `test_import_footprint.py` with `meshioplusplus._cli` and the lazily loaded submodules. Measure with `SKBUILD_EDITABLE_SKIP` set in an editable install, whose `editable.rebuild` otherwise adds a build check to every import.
- **pybind11 per-call overheads.** `mesh_to_py` imports the `Mesh` class on every call (and `regions_to_py` the `Region` class), `ensure_contiguous` reads each array's byte order through a Python attribute lookup, and polygon and polyhedron blocks cross one node id and one face vector at a time (`ragged_cellblock_from_py`), and come back as one list append per polygon or face (`ragged_data_to_py`), where the WASM binding already crosses a CSR triple (`mesh_to_val`). **S–M**
  - **Class cache.** Through `py::gil_safe_call_once_and_store` or a deliberately leaked handle: a plain `static py::object` crashes at interpreter shutdown. `dtype::byteorder()` exists in the installed pybind11; check it against the `pybind11>=2.11` floor.
  - **Clones are semantics.** `transform` clones the connectivity and data it does not change — as 17 operation files do through `clone_mesh`, `clone_geometry` or `quality_clone_mesh` — but that clone is what makes the result a new mesh that aliases nothing, and the Python twin copies too. Dropping it is an API decision, not an optimisation. `transform`'s binding also refuses ragged meshes (`allow_ragged=false`), so those always run the Python twin.
- **R copies the points twice.** `R_mio_points` (`mesh.c`) copies them through `mio_r_copy_as_real`, then `memcpy`s the result into a second matrix, and `conn_matrix` converts the connectivity into one allocation and shifts it to 1-based into another with a scalar loop; `shape_data` in the same file already copies once and sets the `dim` attribute. Fixable behind the same accessor names ([R](./r.md) documents R as copy-only). Julia needs nothing: `points` is one copy and `connectivity` one fused copy-and-shift loop, the minimum for an owning result. **S**
- **The browser viewer round-trips every operation result through a VTP file.** `renderPipeline` (`src/viewer/src/worker/worker.ts`) writes the result to `/out.vtp` in MEMFS as binary, zlib-compressed, base64 VTP; `take()` copies it out with `FS.readFile` and again with a `.slice()` (`readFile` already returns a fresh copy); and `Renderer.load` parses it on the main thread with vtk.js's `XMLPolyDataReader` — an XML parse, a base64 decode and an inflate. Every apply also re-reads the original file and replays the whole pipeline. **L**
  - **A cheap first step (S).** Write the VTP uncompressed or as raw appended data (vtk.js reads it), and drop the redundant `.slice()` behind a guard for a future WASMFS build, whose `readFile` may return a view onto the heap.
  - **The typed-array channel.** A variant of `convertSurfaceOps` that returns `mesh_to_val(surface)` instead of writing a file reuses the typed-array `Mesh` the npm library already returns (`points` as a `Float64Array`, `Int32Array` cells per block, `*_components` for multi-component data). The reason `worker/protocol.ts` and `js_bindings.cpp` give for VTP — that flat meshes cannot carry multi-component arrays — has been stale since `*_components` arrived in v9.9.0.
  - **Watch.** The offline embedded page and the dataset page still need a VTP; the picker needs `surface:parent_cell`; vtk.js wants its own cell layout (a count, then the ids), not per-block arrays; and `mesh_to_val` narrows connectivity to `Int32` and widens points to `Float64`.
  - **Gate.** `tests/viewer/*.spec.mjs` (Playwright: `npm run test:e2e` in `src/viewer`) and `tests/wasm/smoke.mjs`. There is no viewer benchmark: time `applyOps` in `ops.spec.mjs` until the viewer's status changes.

**Deliberately not**, and recorded so it is not re-proposed:

- **Explicit SIMD intrinsics or `-march` flags** — portability across wheels, WASM and the release binaries is worth more than the scalar kernels cost; revisit only if the harness shows a kernel dominating.
- **Kokkos device execution** — every `parallel_for` body captures host pointers (`parallel.hpp`); the GPU route is the DLPack/CuPy handoff ([GPU handoff](./gpu.md)).
- **A BVH in place of the uniform grid** — the tiebreak argument above.
- **Tuning the pure-Python fallback readers** — the fix is making the C++ path accept the file ([§5](#_5-core-parity-across-surfaces)), not a faster fallback; `_decimate.py`'s heap-based twin is deleted once `decimate` accepts its inputs, not optimised.

*Recommended posture:* the boundary items as their consumers ask, and the tokenizer's remaining readers alongside other work in them; the operations' layout-bound leftovers with the next ABI bump. Every item is judged by the determinism check (`--hash`, v16.15.0).

---

## 5. Core parity across surfaces

*Admission: something the Python layer can do that the C++ core cannot, or that the core can do and a binding cannot reach.* A construct that forces the Python fallback is not "slower from C" — it is **unreadable** from C, Fortran, Julia, R, WASM and the native CLI, none of which has a fallback. Ordered by this project's own consumers, Kratos first.

- **MDPA beyond mesh-level blocks.** The C++ core reads and writes `Nodes`/`Elements`/`Conditions`/`SubModelPart`s, but `Begin Table`, `Begin Geometries`, `Begin Mesh <id>`, `Begin Constraints` and non-numeric `ModelPartData` throw — or, under a lenient read, are skipped and listed in `MdpaInfo`, which no flat binding exposes ([MDPA](./formats/mdpa.md#c-core)). **M**
- **Gmsh `$Periodic` and format 4.0 in the C++ core.** `$Periodic`, both directions: a periodic 4.1 file is unreadable from every flat binding today ([Gmsh](./formats/gmsh.md)). The C++ reader also accepts only versions 2.2 and 4.1 (`gmsh.cpp`), so a 4.0 file, which the Python reader reads, is unreadable from them too. Pairs with periodic node matching in [§6](#_6-operations). **S–M**
- **VTK-family constructs the C++ readers refuse**: `<AppendedData>` and multiple `<Piece>`s in `.vtp`, `.vts`, `.vtr` and `.vti` (the `.vtu` reader handles both since v16.6.0; its appended-data decoder in `vtu_read.cpp` is the one to share), and legacy `.vtk` structured points, structured grid and rectilinear grid ([VTU](./formats/vtu.md), [VTK](./formats/vtk.md)). **S–M**
- **XDMF 2 and XPath references** — the C++ core implements XDMF 3 only, and `Reference="XML"` `DataItem`s not at all ([XDMF](./formats/xdmf.md)). **M**
- **MED multi-mesh files and profiles**, which are Python-only and not reachable even under a lenient C++ read ([MED](./formats/med.md)). **M**
- **Netgen extras** — periodic `identifications`, `materials`/`bcnames`/`cd2names`/`cd3names`, `edgesegmentsgi2` and the `.vol.gz` container ([Netgen](./formats/netgen.md)). **S–M**
- **An Exodus writer that carries sets and steps.** It writes element blocks but no node sets or side sets, so only element-block regions round-trip, and it writes one step per file; a multi-step writer is a stateful object of the `XdmfTimeSeriesWriter` shape ([Exodus](./formats/exodus.md)). **M**
- **Sets → regions, phase 2.** XDMF `Sets` are not mapped (UNV groups became regions in v15.6.0 and Ansys components in v16.3.0), and VTU/VTP need a documented region convention ([Named regions](./regions.md)). **M**
- **Named Side regions surviving operations.** `subdivide`, `agglomerate`, `undo_green` and `convert_cells(simplexify)` drop them through their parent-cell remap, and the cutters (`slice`, `isosurface`, `extract_surface`, `extract_skin`) drop them outright ([Named regions](./regions.md)); for a Kratos model, Side regions are where the boundary conditions live. **M**
- **A structured pipeline report on the flat ABI.** C, Fortran, Julia and R receive status plus `mio_last_error()` only; a caller-buffer JSON accessor is recorded as a follow-up, as are the v2 multi-mesh steps (`Inputs:` for `Merge`/`Interpolate`/`UndoGreen`, an `Output.Pattern` for `Split`/partition) ([pipelines](./pipeline.md)). **S–M**
- **PCD `binary_compressed` on the flat ABI.** The C++ API and Python write it (`write_pcd(..., PcdData::BinaryCompressed)`, `data="binary_compressed"`), but `WriteOptions`/`mio_write_opts` only carry the VTK block codecs, so C, Fortran, Julia, R and WASM can read it and cannot write it; the fix is an appended `MIO_CODEC_LZF` enumerator routed to `write_pcd` ([PCD](./formats/pcd.md)). **S**
- **glTF write options on the flat ABI.** `mio_write("x.glb")` and `writeMesh` write the defaults; the colour field, colormap and range, split angle, up axis and unit scale are reachable from Python, C++, both CLIs and MCP only, so C, Fortran, Julia, R and WASM cannot colour a `.glb` or change its axis ([glTF](./formats/gltf.md)). The fix is a `mio_gltf_opts` struct with a reserved tail, following the `mio_curvature_opts` shape. **S**
- **Point/cell sets beyond regions in the core**, so the `convert -s/-d` sets↔data conversions work in the native CLI and flat bindings ([Julia](./julia.md)). **S–M**

---

## 6. Operations

*Admission: a new operation, or a public face for machinery that already exists privately inside one.*

**Generation.** Almost every operation transforms a mesh you already have. The exceptions all start from something else — `grid` from a lattice (`detail/grid_lattice.hpp`), `voxelize`/`compute_sdf` from a surface's bounding box, `remesh_volume` from a closed surface — and nothing builds a shape from parameters, sweeps one, or triangulates a domain.

- **Primitive constructors** — `box`, `sphere`, `cylinder`, `disk`, in their own `operations/primitives.hpp` beside `grid`. Dependency-free, and it removes the fixture-file dependency from tests, docs, notebooks, the browser demo and the MCP server; and it could replace the test-local canonical mesh of the [format conformance matrix](./conformance.md). Highest leverage per line of code in this document. **S**
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
| ParaView | The plugin installs with the wheel (`share/meshioplusplus/paraview/`) and is tested against conda-forge ParaView in CI ([ParaView plugin](./paraview_plugin.md)); no listing on a plugin index yet. |

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

**Scale.** Memory-mapped reads ([mmap](./mmap.md), via `ReadOptions` and the C ABI but not Python's `read()`) roughly halve the peak footprint of a large read, and the XDMF series appender and Python's chunked `write_dataset` write a series or a dataset without holding it — but nothing writes one mesh larger than memory, and no operation streams. Run the benchmark harness's `--tier XL` (about 10M cells) first; it decides whether either item below matters.

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

1. **Format reach ([§1](#_1-format-reach))** — reader/writer parity (triage and the drift guard first), the shared components as formats need them, and the checks against the vendor tools as licences appear.
2. **Spack package upkeep ([§2](#_2-spack-package-upkeep))** — a small, mechanical catch-up (recipes, variants, one release-checklist line) that unblocks HPC users on the current release; then a checklist step, so it stays current.
3. **The OSS-Fuzz submission ([§3](#_3-quality-of-implementation))** — calendar-bound, so it runs alongside everything else.
4. **Performance ([§4](#_4-performance))** — the boundaries, and the tokenizer's remaining readers as they are touched; the operations' layout-bound leftovers with the next ABI bump.
5. **Primitive constructors ([§6](#_6-operations))** — a few days, and a prerequisite of the conformance matrix, every demo surface and `extrude`/`revolve`.
6. **Core parity ([§5](#_5-core-parity-across-surfaces))** — MDPA first, then Side-region survival and sets → regions, then the rest by consumer demand.
7. **Registration ([§7](#_7-ecosystem-reach))** — calendar-bound, so start the submissions early and let them run alongside everything else.
8. **Long-run spikes ([§8](#_8-long-run-spike-first))** — the benchmark tier decides the scale items; the NURBS spike is scheduled independently of the rest.

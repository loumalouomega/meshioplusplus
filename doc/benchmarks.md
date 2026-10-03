# Benchmarks

How does the meshio++ C++ core compare with the original pure-Python [meshio](https://github.com/nschloe/meshio)? The [`benchmark/`](https://github.com/loumalouomega/meshioplusplus/tree/main/benchmark) folder times read and write conversions on the formats that **both** libraries support, on the same in-memory mesh.

Both libraries expose an identical `Mesh` / `read` / `write` API, so the harness ([`benchmark/bench.py`](https://github.com/loumalouomega/meshioplusplus/blob/main/benchmark/bench.py)) hands one geometry to each and times it. The legacy pure-Python meshio is imported from source (it needs no build); meshio++ is the installed package. The headline input is the bundled **`example.msh`** — a real Gmsh mesh of a mechanical bracket (~52k nodes, ~293k cells, mixed triangles + tetrahedra). Reproduce everything with [`benchmark/01_benchmark.ipynb`](https://github.com/loumalouomega/meshioplusplus/blob/main/benchmark/01_benchmark.ipynb).

## Where the C++ core helps (and where it doesn't)

meshio++ moves the parsing/serialising hot loops into C++, so the win is largest exactly where pure-Python is slowest — **text/ASCII formats**:

- **VTU binary + zlib** — the zlib block compression parallelises across cores (the C++ core defaults to an **OpenMP** backend with dynamic scheduling, which also load-balances hybrid P+E-core CPUs), so this is the biggest win: **~16× write**, ~2.3× read.
- **VTU ASCII** — the C++ number formatter and parser are several times faster than the Python/numpy text path (~7× write, ~5× read).
- **XDMF (HDF5)** — reads are much faster on mixed-topology meshes (~10× on the bracket). On a single-block mesh both libraries spend most of their time in gzip, which HDF5 runs on one thread; since v16.15.0 meshio++ writes gzip datasets in ~1 MiB row chunks compressed in parallel, and inflates any row-chunked gzip dataset in parallel on read, so that time now divides across cores (a SEQ build, as in the Linux wheels, still runs it on one).
- **MED (HDF5)** — with the Eigen-backed Fortran↔C transpose fused with the node reorder, MED writes faster (~1.1× in `benchmark/results.csv`). Reads measured 0.81–0.93× there; v16.15.0 removed the two zero-fills of every dataset the reader immediately overwrites and made the family-to-region pass one sweep over the entities instead of one per family, which brings the single-block read level with pure Python.
- **Gmsh binary** — the reader decodes straight from the slurped buffer into an owning array that is *moved* into the cell block (no copy), so reads are **~1.7×** faster; the writer buffers each block into one `write`. `benchmark/results.csv` (July 2026) shows the write at 0.68×, but re-measuring on v16.14.0 against meshio 5.3.5 put it at about 1.9×: the recorded figure does not reproduce, and the notebook's next run replaces it.
- **VTK binary** — writes at parity (fused gather+byte-swap, one `write`). Reads now **beat** pure-Python both on single-cell-type meshes (~1.45×, the connectivity `NDArray` is moved straight into the cell block) and on mixed-topology meshes (**~1.1×** on the bracket, up from ~0.4× originally): reader output buffers skip the zero-fill they immediately overwrite, and the per-type block copy is chunked across threads so its first-touch page faults are serviced concurrently. Endianness conversion uses single-instruction `bswap` intrinsics throughout.

For **plain binary dumps**, pure-Python meshio streams the whole array through numpy's `fromfile`/`tofile` at C speed — a high bar — but with the redundant connectivity copies removed (the reader now *adopts* the byte-swapped buffer as the cell array in the common single-type case), meshio++'s VTK/Gmsh reads now land **at or above parity** there. HDF5 (MED, XDMF) is even-to-faster. MDPA is the ~1× control: a C++ reader exists, but Python's `read()` deliberately keeps calling the pure-Python one (see [formats](./formats.md#when-the-native-path-declines)), so this benchmark measures the same code on both sides of the comparison.

These formats were previously *slower* in meshio++ (VTK/Gmsh binary and MED read all landed at 0.2–0.6×); the current numbers reflect an optimisation pass — bulk-buffered binary I/O (one `write` per section, fused gather+byte-swap; bulk `memcpy` decode on read), **zero-copy cell reconstruction** (the connectivity buffer is reshaped and moved into the cell block, not copied), a real parallel backend (OpenMP by default), thread-capping for the memory-bandwidth-bound loops, and Eigen for the MED transpose. Output stays byte-identical throughout (the round-trip and reference-file tests are the gate).

This is the honest shape of it: meshio++ is a large win for text and compute-bound formats (ASCII, zlib) and now at or above parity on the binary and HDF5 formats too — including single-cell-type binary *reads*, which the zero-copy reconstruction brought level with (or past) numpy's vectorised `fromfile`. Mixed-topology binary reads, which can't adopt the buffer directly, land just under parity.

### v16.18–v16.21: I/O and operation parallelism

Four releases moved more of the remaining serial work into C++ and onto more cores, all byte-identical to what came before (`test_io_baseline.py` pins the writers; `bench_ops --hash` pins the operations across backends and thread counts):

- **A shared text tokenizer** (v16.21.0, `detail/text_cursor.hpp`): `mdpa`, `su2`, `avsucd` and `tecplot` parse `string_view` tokens of the mapped file instead of allocating a `std::string` per line and token — on `bench.py`'s M meshes, su2 reads 3.5× faster, avsucd 2.9×, tecplot 2.1×, and the C++ mdpa reader (used by every flat binding but not by Python's own `read()`, see [formats](./formats.md#when-the-native-path-declines)) 3.2×. VTU ASCII arrays and XDMF `DataItem`s parse straight into their typed array the same way.
- **Nineteen more readers through `FileSource`** (v16.20.0; see [memory-mapped reading](./mmap.md#coverage)), and ASCII writers (VTU, VTP, legacy VTK, Medit, Tecplot, OpenFOAM points, Abaqus, Ansys `.cdb`, LS-DYNA) formatting rows in parallel chunks instead of serially or through a per-row heap string.
- **Raw appended VTU** (v16.21.0, `appended=True`/`--appended`; see [VTU](./formats/vtu.md)): on a 216,000-hexahedron grid, 2.6× faster to write and 6.9× faster to read than inline binary, and 1.6×/1.4× with zlib.
- **Operation parallelism** (v16.18.0–v16.19.0): the distance kernel (`shrinkwrap`, `remesh_volume`, `sample_distance`/`distance_to_surface`), welding (`clean`, `merge`), the shared facet table's counting sort (`extract_surface`, `isosurface`, `slice`), `cell_data_to_point_data`/`gradient`/`hessian`'s gather, `compute_normals`, `remesh_volume` and `remesh`'s setup passes, and the `proximity_graph` pair search (`neighbor_pairs`, in the C++ core since v16.18.0 — see [proximity graphs](./proximity_graphs.md)) all moved from a serial or hash-map-based pass to a parallel one.

## Real mesh (`example.msh`)

Read/write time (log scale) and speedup on the actual bracket mesh:

![read/write timings on example.msh](/benchmarks/benchmark_times.svg)

![speedup on example.msh](/benchmarks/benchmark_speedup.svg)

Speedup = *legacy time / meshio++ time*. Bars in the shaded region mean meshio++ is faster; to the left of the dashed line the pure-Python numpy path wins.

## Does the speedup grow with mesh size?

Both libraries are O(n), so the relative speedup settles to a per-format constant on non-trivial meshes — but the edges behave differently by format:

- **Text formats (VTU ASCII)** — the write speedup *climbs* out of the small-mesh regime as fixed per-call overheads amortise, then plateaus. A large real mesh realises the full speedup; a tiny one does not.
- **Compressed binary (VTU + zlib)** — the OpenMP-parallel zlib compression *grows* with size as there is more work to spread across cores.
- **Plain binary (VTK/Gmsh)** — writes track parity; single-cell-type reads now match or beat numpy (the connectivity buffer is adopted with no copy).

![speedup vs mesh size](/benchmarks/benchmark_scaling.svg)

::: tip Parallel backend
The C++ core parallelises with a compile-time backend (`AUTO` → OpenMP by default). Memory-bandwidth-bound loops (byte-swap, transpose, gather) are thread-capped because they saturate bandwidth after a few threads; compute-bound loops (zlib, base64) use all cores. Check the active backend with `python -c "import meshioplusplus._core as c; print(c.__parallel_backend__)"` — if it prints `stl` without TBB linked, `parallel_for` runs sequentially.
:::

## Reproducing

```sh
uv pip install --python .venv matplotlib jupyter nbconvert ipykernel
cd benchmark
../.venv/bin/jupyter nbconvert --to notebook --execute --inplace 01_benchmark.ipynb
```

The notebook records the machine, library versions, and the inputs (the bundled `example.msh` bracket plus a synthetic tetrahedral cube and a size sweep), runs the harness, writes `results.csv`, and regenerates the plots above. Numbers are single-machine and indicative — the *shape* of the result is the point, not the exact factors.

## Every format

`benchmark/bench.py` also times a write and a read of **every** format meshio++ both writes and reads back, each fed the largest input its [conformance declaration](./conformance.md) says it keeps: the synthetic tetrahedral cube for volume formats, its surface for surface formats (STL, OBJ, PLY, …), its points for point clouds.

```sh
python benchmark/bench.py --sizes S,M --out results_all.csv            # every format
python benchmark/bench.py --sizes L --formats vtu,gmsh22,xdmf,med      # a subset
```

The sizes are 6·(n−1)³ tetrahedra for n = 16, 36 and 56 points per edge (about 20k, 250k and 1M). Legacy meshio is optional here: with `MESHIO_LEGACY_SRC` pointing at a source checkout's `src`, or `meshio` installed, the curated comparison above fills its legacy columns; without it only the meshio++ columns are written.

## Operations

`meshioplusplus_bench_ops` (built with `-DMESHIOPLUSPLUS_BUILD_BENCHMARKS=ON`, `src/cpp/benchmark/bench_ops.cpp`) times the operations on the same cube — a structured Kuhn tetrahedralisation, every tet positively oriented so its surface is consistently wound — or on its surface: `extract_surface`, `extract_skin`, `smooth` (the surface) and `smooth_volume` (a jittered cube), `refine`, `elevate` (`convert_cells` to quadratic), `merge` (welding a mesh onto itself), `clean` and `clean_weld` (two unwelded copies of the cube), `compute_sdf` (48³ grid), `decimate` (half the surface), `partition` (8 parts), `reorder` (RCM) and `reorder_hilbert`, `optimize_volume` (the jittered cube), `agglomerate`, `split` (by component), `undo_green` (a green-closed refinement of every seventh cell), `hessian`, `compute_normals`, `compute_curvature`, `shrinkwrap` (an inflated copy of the surface onto the surface), `voxelize`, and since v16.19.0 `remesh` (a quarter of the surface's points as clusters), `remesh_volume` (1.5 lattice steps), `sample_distance` and `distance_to_surface` (the jittered cube's points against the surface), `isosurface` (two levels) and `slice` of a smooth point field, `linearize` (`convert_cells` back from quadratic), `interpolate` (barycentric, onto the jittered cube), `repair`, `sobolev_deform` (a surface displacement field) and `gradient` (least squares, at the points). Tiers S, M and L are about 10k, 160k and 750k tetrahedra; `--tier XL` (about 10M) is opt-in. Output is one CSV row per tier and operation: `backend,threads,op,cells,median_s,runs,digest`.

The parallel backend is a compile-time choice, so `tools/bench_ops.sh` configures one tree per backend (SEQ, OpenMP, TBB; one that fails to configure is skipped) and sweeps `OMP_NUM_THREADS`, which the benchmark also applies to TBB:

```sh
CMAKE_BUILD_PARALLEL_LEVEL=4 tools/bench_ops.sh benchmark/results_ops.csv "SEQ OPENMP TBB" "1 2 4 8" -- --tier S --tier M
```

### Determinism check

Every operation's output is meant to be byte-identical across parallel backends and thread counts; the serial phases in several operations exist to guarantee that. `--hash` makes each row carry a 64-bit digest of its result — points, connectivity, point/cell/field data and regions in block and sorted-name order, plus the result's index maps and counters — taken from the untimed warmup run, and `tools/bench_ops.sh` then fails if one row's digest differs between any two backends or thread counts. With `BASELINE=<earlier.csv>` it also fails if a digest differs from that earlier sweep, which is how a change that must not alter output proves it:

```sh
CMAKE_BUILD_PARALLEL_LEVEL=4 tools/bench_ops.sh before.csv "SEQ OPENMP TBB" "1 4 8" -- --hash --tier S --tier M --runs 1
# ... change the code ...
BASELINE=before.csv CMAKE_BUILD_PARALLEL_LEVEL=4 tools/bench_ops.sh after.csv "SEQ OPENMP TBB" "1 4 8" -- --hash --tier S --tier M --runs 1
```

The `determinism` job in `ci.yml` runs the first form on the small tier at 1 and 4 threads on every pull request. The digest is also available to the C++ tests (`src/cpp/benchmark/mesh_digest.hpp`), where golden digests pin an operation's output across a rewrite.

## Clang-Tidy performance audit

`tools/performance-tidy.sh` runs the installed Clang-Tidy's built-in [`performance-*` checks](https://clang.llvm.org/extra/clang-tidy/checks/list.html) as an advisory audit, independently of `.clang-tidy`'s include-hygiene gate. It never applies fixes. Findings are candidates, not measured speedups: a cheap `CellView` value is intentional, clones must retain their ownership semantics, and changing a public enum's width, signature or installed inline body needs the [ABI policy](./abi.md), even when the tool offers a replacement.

```sh
# Configure a separate Clang/SEQ tree, then audit its core, C API and native CLI.
tools/performance-tidy.sh --jobs 4
tools/performance-tidy.sh --mesh-backend NATIVE --jobs 4
tools/performance-tidy.sh --mesh-backend KRATOS --jobs 4

# Reuse a configured tree (including Python or optional libraries when enabled).
tools/performance-tidy.sh --build-dir build/performance-tidy-meshio --out build/performance-tidy-meshio/before.json

# A selected-TU scan is explicitly marked as selected, not whole-tree coverage.
tools/performance-tidy.sh --build-dir build/performance-tidy-meshio --file src/cpp/src/formats/ansys.cpp
```

The default build tree is `build/performance-tidy-<mesh-backend>`; CMake module scanning is disabled so GCC-only dependency-scanner flags cannot contaminate Clang-Tidy's command line. Set `CLANG_TIDY`, `CC`, `CXX` or `PYTHON` to choose executables. To cover Python, OpenMP/TBB, WASM or optional dependencies, configure their real compilation databases separately with `-DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCMAKE_CXX_SCAN_FOR_MODULES=OFF`, and pass `--build-dir`; the default configuration does not enable Python, WASM or opt-in libraries. Check the report's actual `compile_defines`, not just `WITH_*=ON` requests, to see which conditional paths were included.

The JSON report records the tool version, enabled checks, explicit Clang-Tidy profile, database hash, CMake settings, actual compile definitions for the scanned commands, scanned translation units and deduplicated source-location/check findings, with the TUs that reported each one. Installed-header findings are labelled for ABI review rather than automatically rejected or accepted. Each TU also has a raw log beside the report. A non-zero tool exit or parsing error makes the scan fail and marks its report incomplete; a new scan invalidates an earlier successful report before configuration or tool startup. A database containing no first-party sources, a stale command or a requested source without a compile command is an error, never an invitation to guess flags. Third-party and generated single-header sources are excluded.

For clangd's editing feedback, opt in locally through a matching fragment in your user `clangd/config.yaml` (for example, `~/.config/clangd/config.yaml` on Linux), substituting your checkout and build paths. This does not change the repository's include-only `.clang-tidy`; clangd's default fast-check filter may skip expensive checks, so use the standalone runner for the exhaustive inventory.

```yaml
If:
  PathMatch: '/path/to/meshioplusplus/(src/cpp|bindings)/.*'
CompileFlags:
  CompilationDatabase: /path/to/meshioplusplus/build/performance-tidy-meshio
Diagnostics:
  ClangTidy:
    Add: 'performance-*'
```

The initial audit (2026-10-03, Clang-Tidy 23.1.1, 22 checks, SEQ, HDF5/netCDF/zlib enabled) completed without parsing failures:

| Mesh backend | First-party TUs | Unique findings | Additional boundary coverage |
| --- | ---: | ---: | --- |
| MESHIO | 189 | 522 | C API, native CLI, Python/pybind11 |
| NATIVE | 188 | 285 | C API, native CLI |
| KRATOS | 188 | 286 | C API, native CLI |

The MESHIO inventory groups as follows; a location can have more than one check, so these are diagnostic counts rather than distinct edits:

| Candidate class | Diagnostics | Triage |
| --- | ---: | --- |
| Value parameters, copied initializations and moves | 287 | Review ownership first; many are pybind11 refcount suggestions, small metadata copies or necessary owned clones. |
| Vector capacity planning | 59 | First batch: the Fluent writer's two mesh-sized lists; small metadata lists and unchecked reader counts are deferred. |
| String concatenation, find/character overloads and view conversion | 91 | Many are error paths or header parsing; no reader throughput win established yet. |
| Enum width | 85 | Installed enum widths are ABI-sensitive; private enums still need evidence before narrowing. |

These counts are a baseline for triage, not a zero-warning target. Public enum-width and nested-vector ingestion-signature suggestions require a separate ABI decision; pybind11 value-parameter suggestions require lifetime/refcount review. Many remaining reservations concern small metadata vectors or untrusted reader counts, and many string-concatenation suggestions are error paths. Prioritize measured hot paths and do not reserve directly from an unchecked file count. Optional CGNS MLL, ADIOS2, TecIO, KaHIP, bzip2, zstd/lz4 and viewer paths, other parallel configurations and WASM remain separate coverage work. For each accepted implementation-only batch, preserve the before report, record timing/allocation evidence, run native and Python parity tests, and apply the [determinism check](#determinism-check) to operation changes. The include-hygiene CI policy remains unchanged.

### First batch: Fluent cell-zone list capacity

`write_ansys` now reserves the mesh-derived size of its cell-type list and, for a mixed-type cell zone, the outer row list. This removes repeated capacity growth without changing ingestion, traversal order or serial/parallel work. The two `performance-inefficient-vector-operation` findings disappear on every mesh backend; complete post-change SEQ scans report 520/283/284 findings for MESHIO/NATIVE/KRATOS, with no new diagnostics. This is an exported non-inline function-body change, with no installed declaration, inline body, signature or layout changes: C++ ABI 22 is unchanged, and the single header is regenerated.

On an AMD Ryzen 7 255 (16 logical CPUs), Clang 23.1.1 Release/SEQ, `benchmark/bench.py --sizes M,L --formats ansys --repeats 7` gives the following binary-writer medians. These differences are small enough that the change is treated as an allocation reduction, not an established throughput speedup; no material SEQ slowdown was observed.

| Input | Cells | Before | After |
| --- | ---: | ---: | ---: |
| M tetrahedral grid | 257,250 | 288.8 ms | 288.6 ms |
| L tetrahedral grid | 998,250 | 1123.0 ms | 1104.1 ms |

Counting ordinary global `operator new` requests around a warmed-up direct native write, with the same output pathname before and after, gives the following reductions on SEQ/libstdc++; ASCII and binary have identical reductions. These are cumulative allocation requests saved per write, **not peak RSS measurements**. The mixed inputs split alternate squares of 420×420 and 816×816 grids into triangles, keep the other squares as quads, and join both types into zone 7.

| Input | Cells | Allocation calls eliminated | Requested bytes eliminated |
| --- | ---: | ---: | ---: |
| M uniform zone | 257,250 | 18 | 1,068,148 (1.02 MiB) |
| L uniform zone | 998,250 | 20 | 4,395,604 (4.19 MiB) |
| M mixed zone | 264,600 | 38 | 21,951,300 (20.93 MiB) |
| L mixed zone | 998,784 | 40 | 30,754,276 (29.33 MiB) |

The [raw auxiliary sweep](https://github.com/loumalouomega/meshioplusplus/blob/main/benchmark/performance_tidy_fluent.csv) preserves 112 before/after rows: both M/L inputs, both zone shapes, ASCII/binary, SEQ at 1 thread, and OpenMP/TBB at 1/4/8 threads. Every SHA-256 agrees across backends, thread counts, repeats and the pre-change baseline (`5726a4a2`). SEQ/TBB use Clang 23.1.1, OpenMP uses GCC 16.2.1; TBB is capped with `tbb::global_control`, not merely `OMP_NUM_THREADS`. Its three-repeat writer timings had background build/test activity and the allocation interposer loaded (counting disabled during timing), so they are retained for transparency, not used to claim a speedup; the seven-repeat uninstrumented harness timings above are the primary SEQ comparison.

Native C++ tests cover ASCII/binary output and an explicit triangle/quad zone on MESHIO, NATIVE and KRATOS; Python tests compare that mixed zone byte-for-byte with the reference writer and retain the existing I/O byte baselines. The new test joins cell blocks using `ansys:zone`, because default block-per-zone inputs do not exercise the mixed-type row list. Run backend CTest suites serially or give them separate temporary directories: `mt::temp_path` is process-local, so concurrent suites sharing one `TMPDIR` can overwrite each other's fixtures.

## In CI

The weekly `benchmark` workflow (also runnable by hand) runs both: every format at size M, and the operations for SEQ, OpenMP and TBB at 1, 2 and 4 threads, with `--hash`. It uploads the CSVs as artifacts and prints them in the job summary. Successful default-branch runs also publish immutable records with commit, run/attempt, machine, compiler, dependencies and benchmark parameters to `benchmark-data`; the [benchmark trends page](./benchmark_trends.md) plots that history through the existing Pages deployment. Its timings never fail a build — a hosted runner is noisy, so they are for trends across runs, not for gating one change — but a digest that differs between backends does. Every [performance](./roadmap.md#_3-performance) item on the roadmap is expected to show its before and after with these tools.

## Mesh-backend benchmarks

The C++ core's [mesh backend](cpp_backends.md) (MESHIO / NATIVE / KRATOS) is an exclusive compile-time choice, so `benchmark/bench_backends.sh` builds one benchmark binary per backend (`src/cpp/benchmark/bench_backends.cpp`, enabled with `-DMESHIOPLUSPLUS_BUILD_BENCHMARKS=ON`) and collates a CSV (`benchmark/results_backends.csv`). Method mirrors the Python harness: warmup + median of 5 (`std::chrono`), a synthetic structured tet cube (default 6·35³ = 257k tets over 46k shared points), and four kinds of rows:

- **ingest** — building the mesh through the uniform ingestion API (the reader side's cost);
- **traverse** — a full writer-side accessor sweep;
- **to_modelpart** (KRATOS only) — the one-time `GetModelPart()` materialization: Nodes, Elements/Conditions, variables, and the automatic tag SubModelParts;
- **write/read** per format — full file round-trips (gmsh 4.1 binary, vtu binary+zlib, vtk binary, medit ASCII, su2).

Representative single-machine numbers (257k tets):

| op | meshio | native | kratos |
| -- | ------ | ------ | ------ |
| ingest | 0.7 ms | 1.0 ms | 0.8 ms |
| traverse | 0.6 ms | 0.7 ms | 0.8 ms |
| to_modelpart | — | — | 65 ms |
| gmsh 4.1 binary write / read | 17 / 3.6 ms | 17 / 8.4 ms | 18 / 2.8 ms |
| vtu (binary+zlib) write / read | 28 / 19 ms | 19 / 38 ms | 25 / 16 ms |
| medit ASCII write / read | 79 / 59 ms | 70 / 63 ms | 89 / 63 ms |

The takeaway: because ingestion is move-based for canonical (Float64/Int64) arrays and the KRATOS backend materializes its ModelPart lazily, **format I/O costs the same under every backend** (differences above are run-to-run noise on parse-bound paths); the only real extra is the explicit, one-time `to_modelpart` conversion — the O(n) entity-creation pass any Kratos exchange has to pay. Reproduce with:

```sh
./benchmark/bench_backends.sh          # optional: grid size, e.g. `... 50`
```

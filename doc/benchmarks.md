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

## In CI

The weekly `benchmark` workflow (also runnable by hand) runs both: every format at size M, and the operations for SEQ, OpenMP and TBB at 1, 2 and 4 threads, with `--hash`. It uploads the CSVs as artifacts and prints them in the job summary. Its timings never fail a build — a hosted runner is noisy, so they are for trends across runs, not for gating one change — but a digest that differs between backends does. Every [performance](./roadmap.md#_4-performance) item on the roadmap is expected to show its before and after with these tools.

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

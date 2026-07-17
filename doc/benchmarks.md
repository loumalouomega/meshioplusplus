# Benchmarks

How does the meshio++ C++ core compare with the original pure-Python [meshio](https://github.com/nschloe/meshio)? The [`benchmark/`](https://github.com/loumalouomega/meshioplusplus/tree/main/benchmark) folder times read and write conversions on the formats that **both** libraries support, on the same in-memory mesh.

Both libraries expose an identical `Mesh` / `read` / `write` API, so the harness ([`benchmark/bench.py`](https://github.com/loumalouomega/meshioplusplus/blob/main/benchmark/bench.py)) hands one geometry to each and times it. The legacy pure-Python meshio is imported from source (it needs no build); meshio++ is the installed package. The headline input is the bundled **`example.msh`** — a real Gmsh mesh of a mechanical bracket (~52k nodes, ~293k cells, mixed triangles + tetrahedra). Reproduce everything with [`benchmark/01_benchmark.ipynb`](https://github.com/loumalouomega/meshioplusplus/blob/main/benchmark/01_benchmark.ipynb).

## Where the C++ core helps (and where it doesn't)

meshio++ moves the parsing/serialising hot loops into C++, so the win is largest exactly where pure-Python is slowest — **text/ASCII formats**:

- **VTU binary + zlib** — the zlib block compression parallelises across cores (the C++ core defaults to an **OpenMP** backend with dynamic scheduling, which also load-balances hybrid P+E-core CPUs), so this is the biggest win: **~16× write**, ~2.3× read.
- **VTU ASCII** — the C++ number formatter and parser are several times faster than the Python/numpy text path (~7× write, ~5× read).
- **XDMF read** — much faster on mixed-topology meshes (~10× on the bracket), roughly even on a single-block mesh.
- **MED (HDF5)** — with the Eigen-backed Fortran↔C transpose fused with the node reorder, MED is now at parity or better (~1.2× write, ~1.0× read).
- **Gmsh binary** — the writer buffers each block into one `write`, and the reader decodes straight from the slurped buffer into an owning array that is *moved* into the cell block (no copy): ~1.4× write, **~1.8× read**.
- **VTK binary** — writes at parity (fused gather+byte-swap, one `write`). Reads now **beat** pure-Python both on single-cell-type meshes (~1.45×, the connectivity `NDArray` is moved straight into the cell block) and on mixed-topology meshes (**~1.1×** on the bracket, up from ~0.4× originally): reader output buffers skip the zero-fill they immediately overwrite, and the per-type block copy is chunked across threads so its first-touch page faults are serviced concurrently. Endianness conversion uses single-instruction `bswap` intrinsics throughout.

For **plain binary dumps**, pure-Python meshio streams the whole array through numpy's `fromfile`/`tofile` at C speed — a high bar — but with the redundant connectivity copies removed (the reader now *adopts* the byte-swapped buffer as the cell array in the common single-type case), meshio++'s VTK/Gmsh reads now land **at or above parity** there. HDF5 (MED, XDMF) is even-to-faster. A Python-only format (MDPA) is the ~1× control.

These formats were previously *slower* in meshio++ (VTK/Gmsh binary and MED read all landed at 0.2–0.6×); the current numbers reflect an optimisation pass — bulk-buffered binary I/O (one `write` per section, fused gather+byte-swap; bulk `memcpy` decode on read), **zero-copy cell reconstruction** (the connectivity buffer is reshaped and moved into the cell block, not copied), a real parallel backend (OpenMP by default), thread-capping for the memory-bandwidth-bound loops, and Eigen for the MED transpose. Output stays byte-identical throughout (the round-trip and reference-file tests are the gate).

This is the honest shape of it: meshio++ is a large win for text and compute-bound formats (ASCII, zlib) and now at or above parity on the binary and HDF5 formats too — including single-cell-type binary *reads*, which the zero-copy reconstruction brought level with (or past) numpy's vectorised `fromfile`. Mixed-topology binary reads, which can't adopt the buffer directly, land just under parity.

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

## Mesh-backend benchmarks

The C++ core's [mesh backend](cpp_backends.md) (MESHIO / NATIVE / KRATOS) is an exclusive compile-time choice, so `benchmark/bench_backends.sh` builds one benchmark binary per backend (`cpp/benchmark/bench_backends.cpp`, enabled with `-DMESHIOPLUSPLUS_BUILD_BENCHMARKS=ON`) and collates a CSV (`benchmark/results_backends.csv`). Method mirrors the Python harness: warmup + median of 5 (`std::chrono`), a synthetic structured tet cube (default 6·35³ = 257k tets over 46k shared points), and four kinds of rows:

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

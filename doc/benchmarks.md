# Benchmarks

How does the meshio++ C++ core compare with the original pure-Python
[meshio](https://github.com/nschloe/meshio)? The
[`benchmark/`](https://github.com/<org>/meshioplusplus/tree/main/benchmark)
folder times read and write conversions on the formats that **both** libraries
support, on the same in-memory mesh.

Both libraries expose an identical `Mesh` / `read` / `write` API, so the harness
([`benchmark/bench.py`](https://github.com/<org>/meshioplusplus/blob/main/benchmark/bench.py))
hands one geometry to each and times it. The legacy pure-Python meshio is
imported from source (it needs no build); meshio++ is the installed package. The
headline input is the bundled **`example.msh`** — a real Gmsh mesh of a
mechanical bracket (~52k nodes, ~293k cells, mixed triangles + tetrahedra).
Reproduce everything with
[`benchmark/01_benchmark.ipynb`](https://github.com/<org>/meshioplusplus/blob/main/benchmark/01_benchmark.ipynb).

## Where the C++ core helps (and where it doesn't)

meshio++ moves the parsing/serialising hot loops into C++, so the win is
largest exactly where pure-Python is slowest — **text/ASCII formats**:

- **VTU binary + zlib** — the zlib block compression parallelises across cores
  (the C++ core defaults to an **OpenMP** backend), so this is now the biggest
  win: **~12× write** on the bracket (~16× on a larger cube), ~2× read.
- **VTU ASCII** — the C++ number formatter and parser are several times faster
  than the Python/numpy text path (~7× write, ~5× read).
- **XDMF read** — much faster on mixed-topology meshes (~10× on the bracket),
  roughly even on a single-block mesh.
- **MED (HDF5)** — with the Eigen-backed Fortran↔C transpose fused with the
  node reorder, MED is now at parity or better (~1.2× write, ~1.0× read).
- **Gmsh binary** — the writer buffers each block into one `write` instead of a
  stream call per scalar: ~1.3× write; reads land ~0.8×.

For **plain binary dumps** (legacy VTK binary) pure-Python meshio streams the
whole array through numpy's `fromfile`/`tofile` at C speed, which is hard to
beat — meshio++ writes are now at parity (~0.9–1.1×) and reads land ~0.7–0.8×
(MED read is at parity, ~0.9–1.2×), since numpy's single-pass `fromfile` remains
the ceiling for the reader's parse + byte-swap + cell reconstruction. A
Python-only format (MDPA) is the ~1× control.

These formats were previously *slower* in meshio++ (VTK/Gmsh binary and MED read
all landed at 0.2–0.6×); the current numbers reflect an optimisation pass —
bulk-buffered binary I/O (one `write` per section, fused gather+byte-swap; bulk
`memcpy` decode and an identity-remap fast path on read), passing the int64
connectivity buffer straight to cell reconstruction (no copy), a real parallel
backend (OpenMP by default), thread-capping for the memory-bandwidth-bound
loops, and Eigen for the MED transpose. Output stays byte-identical throughout
(the round-trip and reference-file tests are the gate).

This is the honest shape of it: meshio++ is a large win for text and
compute-bound formats (ASCII, zlib) and now at least at parity on the binary and
HDF5 formats, with plain-binary *reads* the one place numpy's vectorised I/O
still leads.

## Real mesh (`example.msh`)

Read/write time (log scale) and speedup on the actual bracket mesh:

![read/write timings on example.msh](/benchmarks/benchmark_times.svg)

![speedup on example.msh](/benchmarks/benchmark_speedup.svg)

Speedup = *legacy time / meshio++ time*. Bars in the shaded region mean meshio++
is faster; to the left of the dashed line the pure-Python numpy path wins.

## Does the speedup grow with mesh size?

Both libraries are O(n), so the relative speedup settles to a per-format
constant on non-trivial meshes — but the edges behave differently by format:

- **Text formats (VTU ASCII)** — the write speedup *climbs* out of the
  small-mesh regime as fixed per-call overheads amortise, then plateaus. A large
  real mesh realises the full speedup; a tiny one does not.
- **Compressed binary (VTU + zlib)** — the OpenMP-parallel zlib compression
  *grows* with size as there is more work to spread across cores.
- **Plain binary (VTK)** — writes track parity; reads stay a little below,
  since numpy's fully vectorised `fromfile` is a hard single-pass baseline.

![speedup vs mesh size](/benchmarks/benchmark_scaling.svg)

::: tip Parallel backend
The C++ core parallelises with a compile-time backend (`AUTO` → OpenMP by
default). Memory-bandwidth-bound loops (byte-swap, transpose, gather) are
thread-capped because they saturate bandwidth after a few threads; compute-bound
loops (zlib, base64) use all cores. Check the active backend with
`python -c "import meshioplusplus._core as c; print(c.__parallel_backend__)"` —
if it prints `stl` without TBB linked, `parallel_for` runs sequentially.
:::

## Reproducing

```sh
uv pip install --python .venv matplotlib jupyter nbconvert ipykernel
cd benchmark
../.venv/bin/jupyter nbconvert --to notebook --execute --inplace 01_benchmark.ipynb
```

The notebook records the machine, library versions, and the inputs (the bundled
`example.msh` bracket plus a synthetic tetrahedral cube and a size sweep), runs
the harness, writes `results.csv`, and regenerates the plots above. Numbers are
single-machine and indicative — the *shape* of the result is the point, not the
exact factors.

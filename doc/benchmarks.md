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

- **VTU ASCII** — the C++ number formatter and parser are several times faster
  than the Python/numpy text path (~8× write, ~5× read on the bracket).
- **VTU binary + zlib** — a solid, consistent speedup (~1.7× each way).
- **XDMF read** — much faster on mixed-topology meshes (~10× on the bracket),
  though roughly even on a single-block mesh.

For **plain binary dumps** (legacy VTK binary, single-type Gmsh binary) the
pure-Python meshio already reads and writes with numpy's `fromfile`/`tofile` at
C speed, so there is little left to win — and the pybind11 boundary can make
meshio++ *slightly slower* on the very simplest layouts. **HDF5** formats
(XDMF-HDF, MED) go through the same HDF5 C library either way, so they are
roughly even. A Python-only format (MDPA) is the ~1× control — neither library
has a C++ path for it.

This is the honest shape of it: meshio++ is a large win for text-heavy formats
and complex readers, roughly neutral for numpy-friendly binary blobs, and never
changes the file contents — the fallbacks guarantee byte-compatible output.

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
  small-mesh regime (~5× → ~8×) as fixed per-call overheads amortise, then
  plateaus. A large real mesh realises the full speedup; a tiny one does not.
- **Compressed binary (VTU + zlib)** — roughly flat (~1.7×).
- **Plain binary (VTK)** — pure-Python meshio streams it through numpy's
  fully vectorised `fromfile`/`tofile`, which only pulls *further* ahead as the
  mesh grows, so meshio++'s ratio drifts a little lower with size here.

![speedup vs mesh size](/benchmarks/benchmark_scaling.svg)

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

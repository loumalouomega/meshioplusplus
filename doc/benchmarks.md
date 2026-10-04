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

### OFF, Medit, XYZ and UNV: second reader batch

The [112-row matrix](https://github.com/loumalouomega/meshioplusplus/blob/main/benchmark/text_io_batch2.csv) covers these four readers on the same M/L harness inputs, compiler settings, SEQ/OpenMP/TBB modules and thread limits as the MFM experiment. Before modules come from the preserved original archive, with the four source files verified against its source commit; after modules and source snapshots are separately archived under `build/text-io-batch2/`. Timing medians use seven uninstrumented repeats after warmup. Allocation counts are three identical warmed native-dispatch reads in separate instrumented processes; they include Python dispatch/binding overhead and are cumulative ordinary `operator new` requests, not peak memory. Digests include geometry, point/cell/field data and regions and agree across both stages and all seven configurations for every format/size pair.

| SEQ reader | M read before → after | L read before → after | M allocation calls before → after | L allocation calls before → after |
| --- | ---: | ---: | ---: | ---: |
| OFF | 3.6 → 0.9 ms | 9.0 → 2.5 ms | 22,097 → 47 | 54,498 → 48 |
| Medit ASCII | 38.0 → 33.8 ms | 148.7 → 132.5 ms | 140,017 → 49 | 526,889 → 41 |
| XYZ | 1.8 → 1.0 ms | 4.5 → 2.5 ms | 58,455 → 29,487 | 144,621 → 72,693 |
| UNV | 95.5 → 99.6 ms | 387.2 → 394.5 ms | 5,377,168 → 5,237,200 | 20,778,542 → 20,251,694 |

OFF replaces both file extraction and its count-line stream with `TextStream` over a `FileSource`. Medit retains its comment-aware tokenizer but returns views and parses bounded prefixes. XYZ retains comment/column/delimiter/PTS and molecular-file rules while viewing its source lines and numeric tokens; header names remain owned and the output arrays are copied. UNV normalizes short Fortran real fields into a terminated stack buffer instead of an owned string; long fields still have an owned fallback. None changes writer code, public API or installed headers, so C++ ABI 22 is unchanged. The remaining private-cursor consolidation and dtype-switch work are not part of this batch.

The [32-row SEQ confirmation](https://github.com/loumalouomega/meshioplusplus/blob/main/benchmark/text_io_batch2_seq_confirmation.csv) records before→after and after→before orders rather than replacing the original sweep. OFF, Medit and XYZ retain read gains in both orders. UNV is an allocation-only result, not a speedup: M reads are 94.43→95.04 and 98.51→96.21 ms; L reads are 388.35→389.39 and 385.29→391.07 ms, with no material slowdown. Writer medians fluctuate despite identical implementations, notably UNV M (210.55→214.77 and 179.73→211.94 ms); no writer speedup is claimed. Keep those noisy observations in the evidence. Direct native suites run on MESHIO/NATIVE/KRATOS; fresh-process Python tests exercise both forced buffered and mmap-attempt paths, CRLF/no-final-newline and owned results after unlinking and replacing the source.

### EnSight, VTK, Ansys, Patran and Tecplot: stream batch

The [140-row matrix](https://github.com/loumalouomega/meshioplusplus/blob/main/benchmark/text_io_batch3.csv) and [40-row SEQ confirmation](https://github.com/loumalouomega/meshioplusplus/blob/main/benchmark/text_io_batch3_seq_confirmation.csv) preserve matching file/mesh digests before and after on all seven configurations. The driver uses the harness's M/L inputs with `binary=False` for EnSight, VTK and Ansys, so it actually measures their ASCII paths; Tecplot and Patran use their usual text writers. Archived modules/source/compile commands and raw logs are in `build/text-io-batch3/`. The timing/count methodology is the same as the second batch, including separate instrumentation and both SEQ measurement orders.

| SEQ reader | M read before → after | L read before → after | M allocation calls before → after | L allocation calls before → after |
| --- | ---: | ---: | ---: | ---: |
| EnSight ASCII | 16.3 → 16.9 ms | 60.5 → 58.9 ms | 57 → 54 | 49 → 46 |
| VTK ASCII | 20.0 → 20.7 ms | 76.2 → 75.1 ms | 67 → 61 | 59 → 53 |
| Ansys ASCII | 408.2 → 416.4 ms | 2127.8 → 2164.5 ms | 6,615,646 → 6,615,645 | 25,592,351 → 25,592,350 |
| Patran mesh | 292.3 → 287.7 ms | 1353.3 → 1330.9 ms | 5,611,802 → 5,611,802 | 21,721,422 → 21,721,422 |
| Tecplot | 57.6 → 57.2 ms | 455.8 → 443.2 ms | 257,406 → 257,406 | 998,402 → 998,402 |

These are small stream-overhead cleanups, not bulk-reader speedups. EnSight's case/time/id records, VTK's header-token splitter, Ansys's zone-name record and Patran's result header now use `TextStream`. Patran's mesh-only harness does **not** exercise result headers; its existing direct C++ and native/reference Python text/binary-result tests cover that migration, with no result-reader performance claim. Tecplot's allocating stream helper had no call sites: removing it has no runtime effect, and its counts remain identical. Confirmation SEQ read medians do not show a material slowdown; unchanged writer timings remain noisy (including EnSight and Patran L) and are retained rather than interpreted as a gain. All changes are core-private/non-inline implementations; API, installed headers and ABI 22 are unchanged. The remaining allocating tokenizers and private-cursor consolidation remain open.

### Femap, PCD and Radioss engine records: fourth batch

Femap's comma fields now view its source records, with bounded real parsing and owned titles/results. PCD's header fields view the source; ASCII data, header integers and viewpoint numbers use a terminated stack buffer for short tokens and an owned long-token fallback, retaining the existing C parser's range/errno and embedded-NUL semantics. Radioss engine records use shared blank-separated views and bounded whole-token floating parsing. Radioss fixed/comma fields and include-expanded source lines are **not** migrated in this batch. No writer, installed header or public API changes; ABI 22 remains unchanged.

The [84-row matrix](https://github.com/loumalouomega/meshioplusplus/blob/main/benchmark/text_io_batch4.csv) measures Femap, **ASCII** PCD (`binary=False`) and Radioss **starter** meshes on the harness's M/L inputs, with the same seven configurations and separate timing/allocation instrumentation. All before/after file and full parsed-data digests agree. Archives are under `build/text-io-batch4-ascii/`. The starter harness does not exercise engine fields: Radioss's direct native/native-reference engine tests and the new long-token/ignored-nonnumeric test cover that path, but no engine performance gain is claimed.

| SEQ reader | M read before → after | L read before → after | M allocation calls before → after | L allocation calls before → after |
| --- | ---: | ---: | ---: | ---: |
| Femap | 328.1 → 211.7 ms | 1482.7 → 968.5 ms | 8,480,531 → 6,197,003 | 32,859,077 → 24,013,789 |
| PCD ASCII | 1.0 → 1.1 ms | 2.7 → 2.6 ms | 106 → 96 | 106 → 96 |
| Radioss starter | 232.6 → 250.9 ms | 1140.5 → 1141.0 ms | 5,267,765 → 5,267,765 | 20,742,055 → 20,742,055 |

Femap cumulative requested bytes drop from 1,713,613,273 to 973,491,843 for M and 6,683,154,139 to 3,808,758,163 for L (about 43%); these are not peak-memory measurements. The [24-row SEQ confirmation](https://github.com/loumalouomega/meshioplusplus/blob/main/benchmark/text_io_batch4_seq_confirmation.csv) retains both orders: Femap M 328.75→193.21 and 319.86→190.42 ms, L 1506.31→978.33 and 1487.66→974.36 ms. Radioss's initial M slowdown does not reproduce (232.27→231.49 and 229.37→231.69 ms). PCD is a small allocation result, not a speedup claim: M 1.03→0.98 and 1.04→0.98 ms; L 2.51→2.69 and 2.59→2.55 ms, with the small one-order regression retained as timing variability rather than hidden. Native tests run on all three mesh backends, and the fresh-process ownership/prefix/error tests also pass on the original archived SEQ module, proving they pin existing behavior rather than broaden parsing.

### DEX, IP, FLUX and PERMAS: line-view batch

The [112-row matrix](https://github.com/loumalouomega/meshioplusplus/blob/main/benchmark/text_io_batch5.csv) and [32-row both-order confirmation](https://github.com/loumalouomega/meshioplusplus/blob/main/benchmark/text_io_batch5_seq_confirmation.csv) retain identical before/after file and full mesh digests on all seven configurations. All four readers now view `FileSource` lines instead of owning every line; FLUX and PERMAS also view blank-separated body tokens. IP still owns each parenthesis-normalized record and its column names; DEX owns header metadata. Existing numeric-prefix and stream-extraction behavior, writer bytes, public API and installed headers are unchanged, leaving ABI 22 unchanged. DEX requires a nodal field: the harness now supplies `benchmark:field = arange(npoints)` only for that input, without mutating the shared meshes or changing conformance.

| SEQ reader | M read before → after | L read before → after | M allocation calls before → after | L allocation calls before → after |
| --- | ---: | ---: | ---: | ---: |
| DEX | 9.3 → 8.6 ms | 37.7 → 33.2 ms | 233,124 → 186,567 | 877,949 → 702,411 |
| IP | 1.5 → 1.1 ms | 3.6 → 2.6 ms | 39,108 → 26,093 | 101,030 → 67,375 |
| FLUX | 159.2 → 66.9 ms | 554.9 → 261.9 ms | 985,651 → 257,419 | 3,825,865 → 998,433 |
| PERMAS | 66.7 → 37.7 ms | 268.7 → 153.6 ms | 2,470,195 → 1,987,474 | 9,547,803 → 7,690,314 |

The first FLUX candidate allocated a temporary token vector per record, increasing M/L allocation calls to 2,455,382/9,511,276 despite a read gain. Its original raw matrix is preserved in `build/text-io-batch5/after/`; the published matrix uses `after-refined/`, which reuses that vector and reduces requested bytes to 217,824,419/868,814,051 from 433,964,810/1,728,718,042. This is cumulative allocation traffic, not peak RSS. Before source/settings and all archived modules remain in the same tree. Confirmation reproduces all four reader gains, including FLUX M 163.71→65.01 and 162.14→65.65 ms, L 545.95→259.36 and 553.04→262.12 ms. Writer timing fluctuations are retained, including FLUX L 830.57→885.25 ms in the reversed order; writer code is unchanged and no writer performance gain is claimed. Focused native tests cover MESHIO/NATIVE/KRATOS, and Python checks cover source ownership, CRLF/no-final-newline and forced buffered/mmap-attempt reads.

### Abaqus, Ansys coded databases, Marc, Nastran and Netgen: deck-line batch

These five readers now view their mapped/buffered source lines. Marc retains each included file in stable `deque<FileSource>` storage through deck parsing and keeps its formatted post-file source alive with the post reader; returned names, sets and arrays remain owned. Netgen views its inflated gzip buffer instead of copying it into a string stream, retaining concatenated-member handling and count bounds. Abaqus data-section rows and Nastran flattened logical cards remain owned; their tokenizers are not closed by this batch. No installed headers/API/writer changes; ABI 22 is unchanged.

The [140-row matrix](https://github.com/loumalouomega/meshioplusplus/blob/main/benchmark/text_io_batch6.csv) has identical before/after file and mesh digests across all seven configurations. The source/module/settings archives are in `build/text-io-batch6-lines/`; the published after archive is `after-refined/`. The first Netgen candidate copied each view into a temporary string before stripping it, leaving M allocation calls essentially unchanged; its raw evidence remains in `after/`, and the refinement strips the view directly. The harness uses normal text/deck writers, not gzip, includes or Marc post results; those paths are covered by native/native-reference fixtures and fresh-process source-ownership tests, including a nested Marc include.

| SEQ reader | M read before → after | L read before → after | M allocation calls before → after | L allocation calls before → after |
| --- | ---: | ---: | ---: | ---: |
| Abaqus | 93.3 → 89.5 ms | 369.1 → 322.9 ms | 2,711,303 → 2,407,418 | 10,444,757 → 9,270,901 |
| Ansys coded database | 179.8 → 166.1 ms | 718.0 → 597.4 ms | 5,143,990 → 4,840,073 | 19,900,568 → 18,726,691 |
| Marc deck | 156.2 → 154.9 ms | 628.8 → 567.4 ms | 4,242,675 → 3,938,765 | 16,405,217 → 15,231,347 |
| Nastran | 212.8 → 206.4 ms | 871.7 → 798.6 ms | 7,338,703 → 6,988,158 | 28,345,645 → 26,996,182 |
| Netgen | 67.6 → 56.9 ms | 267.9 → 213.3 ms | 2,220,753 → 1,916,854 | 8,568,399 → 7,394,529 |

The [40-row confirmation](https://github.com/loumalouomega/meshioplusplus/blob/main/benchmark/text_io_batch6_seq_confirmation.csv) retains both SEQ orders and reader gains. Netgen M reads are 65.47→55.38 and 65.70→54.95 ms, L 258.96→210.18 and 259.36→215.44 ms. All readers reduce cumulative requested bytes as well as allocation calls; these are not peak-memory measurements. Unchanged writer timings remain noisy, including Abaqus L 133.99→163.92 ms in the reversed order; no writer speedup is claimed. The five-format/ownership/gate tests pass on SEQ/OpenMP/TBB and direct native suites run on MESHIO/NATIVE/KRATOS.

### LS-DYNA, Radioss card fields and GiD quoted tokens

The [84-row matrix](https://github.com/loumalouomega/meshioplusplus/blob/main/benchmark/text_io_batch7.csv) and [24-row both-order confirmation](https://github.com/loumalouomega/meshioplusplus/blob/main/benchmark/text_io_batch7_seq_confirmation.csv) match file/full-mesh digests across both stages and all seven configurations. This batch uses the immediately preceding rebuilt modules as its baseline, with source snapshots and implementation diffs at both stages in `build/text-io-batch7/`, because Radioss's engine fields were already migrated. Compiler commands and settings match. LS-DYNA views lines and block records over each include's source, explicitly retaining the original split's final empty record. Radioss fixed/comma fields view the stable expanded deck lines; numeric parsing still calls the existing owned card parser. GiD views quoted/unquoted tokens over each record and parses short numbers from terminated stack buffers; metadata names and long-number fallbacks remain owned. No writer, installed-header or API changes; ABI 22 is unchanged.

| SEQ reader | M read before → after | L read before → after | M allocation calls before → after | L allocation calls before → after |
| --- | ---: | ---: | ---: | ---: |
| LS-DYNA | 198.2 → 188.6 ms | 934.8 → 1000.3 ms | 3,449,240 → 2,841,425 | 13,372,224 → 11,024,489 |
| Radioss starter | 230.3 → 241.3 ms | 1139.8 → 1098.1 ms | 5,267,759 → 5,267,759 | 20,742,052 → 20,742,052 |
| GiD ASCII | 126.9 → 137.7 ms | 561.4 → 576.4 ms | 1,823,380 → 1,823,380 | 7,043,186 → 7,043,186 |

The initial timing regressions are retained rather than replaced. Both-order confirmation reproduces no material reader slowdown: LS-DYNA M 192.61→173.98 and 193.61→172.06 ms, L 947.52→890.69 and 938.45→877.37 ms; Radioss M 231.70→220.86 and 233.06→219.69 ms, L 1158.20→1099.74 and 1144.19→1098.79 ms; GiD M 129.59→124.70 and 129.95→122.82 ms, L 610.72→556.40 and 605.44→555.03 ms. Cumulative requested bytes decrease for every reader/configuration, even where allocation call counts stay identical because the numeric strings used small-string storage. They are not peak-memory measurements. Unchanged writer timing noise is retained, with no writer gain claimed. Ownership/prefix tests also pass on the baseline SEQ module; direct native tests cover all three mesh backends, and Python suites include GiD ASCII/binary/HDF5 paths and LS-DYNA/Radioss include/metadata behavior.

### Shared-card views and bounded Fortran numbers

The [196-row matrix](https://github.com/loumalouomega/meshioplusplus/blob/main/benchmark/text_io_batch8.csv) and [56-row both-order confirmation](https://github.com/loumalouomega/meshioplusplus/blob/main/benchmark/text_io_batch8_seq_confirmation.csv) preserve file/full-mesh digests across both stages and all seven configurations. Six readers (LS-DYNA, Ansys coded databases, Marc, Radioss, Patran and UNV) use core-private card/token views and short terminated number buffers; Nastran also uses the bounded numeric counterpart while keeping its reassembled logical cards owned. The installed `keyword_card.hpp` owning API is unchanged. Direct tests compare values, signed zero, error messages, D/implicit exponents, range limits, embedded NULs and stack/owned-buffer boundaries against that API. No writer/public-header/signature/default changes; ABI 22 is unchanged. Source/module/settings archives are in `build/text-io-batch8-cards/`, against the preceding rebuilt modules.

| SEQ reader | M read before → after | L read before → after | M allocation calls before → after | L allocation calls before → after |
| --- | ---: | ---: | ---: | ---: |
| LS-DYNA | 171.6 → 160.0 ms | 891.7 → 916.4 ms | 2,841,428 → 2,328,212 | 11,024,492 → 9,029,996 |
| Ansys coded database | 152.4 → 146.6 ms | 603.5 → 563.7 ms | 4,840,073 → 4,700,105 | 18,726,691 → 18,199,843 |
| Marc deck | 147.1 → 124.7 ms | 581.9 → 473.4 ms | 3,938,765 → 2,806,066 | 15,231,347 → 10,860,888 |
| Radioss starter | 221.2 → 245.5 ms | 1210.2 → 1198.6 ms | 5,267,765 → 4,882,853 | 20,742,055 → 19,246,183 |
| Patran neutral | 288.2 → 294.0 ms | 1414.9 → 1821.8 ms | 5,611,802 → 5,611,802 | 21,721,422 → 21,721,422 |
| Nastran | 197.9 → 213.2 ms | 839.1 → 895.8 ms | 6,988,158 → 6,754,878 | 26,996,182 → 26,055,382 |
| UNV | 91.2 → 103.1 ms | 374.7 → 400.9 ms | 5,237,200 → 5,237,200 | 20,251,694 → 20,251,694 |

Initial regressions are retained, including Patran L's reader 1414.88→1821.78 ms and unchanged writer 1221.52→1413.31 ms. Confirmation reproduces reader gains for LS-DYNA, Ansys, Marc and Patran, with no material confirmed SEQ slowdown for the other three. Patran M reads are 293.46→263.01 and 289.81→264.15 ms, L 1330.39→1219.53 and 1274.99→1190.54 ms. Marc L reads are 582.27→457.77 and 566.54→465.27 ms. Nastran's small increases (M 193.94→197.97, L 784.05→796.93 ms in the reversed order) and all writer fluctuations remain in the confirmation; Nastran/Radioss/UNV have no speedup claim. UNV's benchmark takes the already migrated free-field mesh path, so unchanged allocation counts/bytes are expected; fixed-record fixtures cover its migrated paths. Patran's call counts stay identical while requested bytes fall from 703,906,827→457,390,491 (M) and 2,736,711,739→1,783,419,723 (L), primarily smaller token-vector storage. These are cumulative requests, not peak RSS. Marc post results and other non-mesh records are covered by fixtures, not those mesh timing rows.

The initial confirmation launch accidentally permitted an editable-install rebuild, which failed on the full `/tmp` filesystem before any confirmation timing row was collected. Its failure log remains in the archive; recovery explicitly disables editable rebuilds and uses project-local temporary storage with the same archived modules and measurement paths. The focused gates passed 908 Python tests (four optional skips) per SEQ/OpenMP/TBB module and 118 native tests per MESHIO/NATIVE/KRATOS backend. The remaining allocating Abaqus/Netgen tokenizers and Nastran logical-card storage are separate work, not claimed as closed here.

### Shared byte and record cursors

Eight format-private adapters now share the core-private `TextCursor`/`RecordCursor` implementations: Gmsh, VTK legacy, EnSight ASCII, UNV, MFEM, MDPA, GiD and Femap. They retain comments, quoting, binary positioning, diagnostics and any deliberately owned normalized names/tokens. Gmsh still reads integer identifiers through the floating path (`2.0` and `1e3`), and tests cover hexadecimal coordinates, signed zero, prefix position/failure behavior, long fields, embedded NULs and source replacement. Installed headers, writer code/bytes and public APIs are unchanged; ABI 22 is unchanged. The [224-row matrix](https://github.com/loumalouomega/meshioplusplus/blob/main/benchmark/text_io_batch9.csv) and [64-row both-order confirmation](https://github.com/loumalouomega/meshioplusplus/blob/main/benchmark/text_io_batch9_seq_confirmation.csv) match file/full exposed-mesh digests across both stages and all seven configurations, against the preceding rebuilt modules with matching compiler commands. Source/module/settings archives are in `build/text-io-batch9-cursors/`.

| SEQ reader | M read before → after | L read before → after |
| --- | ---: | ---: |
| Gmsh 4.1 ASCII | 22.2 → 20.5 ms | 89.1 → 74.8 ms |
| VTK legacy ASCII | 19.9 → 16.6 ms | 78.2 → 60.0 ms |
| EnSight ASCII | 16.4 → 13.6 ms | 62.7 → 48.9 ms |
| UNV | 95.4 → 102.8 ms | 404.9 → 407.1 ms |
| MFEM | 113.5 → 129.1 ms | 854.5 → 861.9 ms |
| MDPA native | 31.8 → 36.0 ms | 125.3 → 130.7 ms |
| GiD ASCII | 126.4 → 140.9 ms | 694.5 → 662.1 ms |
| Femap neutral | 204.7 → 215.2 ms | 1078.1 → 970.4 ms |

Allocation calls and cumulative requested bytes are identical for every format/configuration: this is consolidation/bounded prefix parsing, not an allocation-cleanup claim. Both SEQ orders reproduce Gmsh, VTK and EnSight reader gains: Gmsh L 82.81→69.99 and 81.28→69.03 ms, VTK L 76.04→54.05 and 75.36→53.93 ms, EnSight L 58.35→44.40 and 58.66→44.82 ms. The other adapters are neutral within the retained fluctuations, not speedups; GiD's reversed-order increases (M 118.42→122.72, L 533.69→550.47 ms) oppose the first order's modest gains, while Femap L rises 915.60→922.10 and 918.92→923.10 ms. Initial regressions and unchanged-writer noise remain published, including Gmsh L writes 297.61→390.37 ms and VTK L confirmation writes 158.98→201.77 ms.

Gmsh 4.1 uses a directly tested single-type tetrahedral input, without changing its mixed-cell conformance rejection. MDPA is explicitly read through `_core.mdpa_read` at both stages: its public Python shim deliberately uses the Python reference and would not measure the native cursor. The CSV identifies that reader; its digest covers the standard mesh exposed by the binding, while native MDPA tests cover property/side-channel behavior. The initial input-selection and Python-MDPA digest failures remain in the archive, before any complete matrix was recorded. ASCII selection is explicit for Gmsh, VTK and EnSight; binary/other grammar paths are covered by fixtures, not those timing rows. Focused gates passed 1,147 Python tests with one expected failure per SEQ/OpenMP/TBB module, 240 native tests per MESHIO/NATIVE/KRATOS backend and 174 ownership/compatibility tests on the archived baseline SEQ module.

### Remaining deck-token views

Abaqus comma fields/data-section rows, Netgen whitespace fields, Nastran logical-card chunks/fields and FLAC3D ASCII rows now use bounded views. Nested include sources and inflated gzip buffers remain live through parsing. Ordinary Nastran large fields are adjacent source views; tolerated mixed free/fixed continuations retain stable owned concatenations in a deque rather than extending a view across commas. Duplicates, trailing empty fields, numeric-prefix leniency, names/sets, connectivity and writers are unchanged. The [112-row matrix](https://github.com/loumalouomega/meshioplusplus/blob/main/benchmark/text_io_batch10.csv) and [32-row both-order confirmation](https://github.com/loumalouomega/meshioplusplus/blob/main/benchmark/text_io_batch10_seq_confirmation.csv) preserve file/full-mesh digests across stages and all seven configurations. They compare the preceding rebuilt modules with matching compiler commands/settings; archives are in `build/text-io-batch10-tokens/`. No installed headers/public APIs changed; ABI 22 is unchanged.

| SEQ reader | M read before → after | L read before → after | M allocation calls before → after | L allocation calls before → after |
| --- | ---: | ---: | ---: | ---: |
| Abaqus | 86.0 → 52.1 ms | 335.7 → 198.3 ms | 2,407,418 → 1,613,016 | 9,270,901 → 6,220,739 |
| Nastran | 233.7 → 196.1 ms | 916.9 → 723.0 ms | 6,754,878 → 6,171,676 | 26,055,382 → 23,828,820 |
| Netgen | 59.0 → 38.7 ms | 232.8 → 141.4 ms | 1,916,854 → 1,426,326 | 7,394,529 → 5,518,201 |
| FLAC3D ASCII | 70.2 → 41.6 ms | 285.3 → 153.7 ms | 2,220,781 → 818,521 | 8,568,429 → 3,170,489 |

Both SEQ orders reproduce reader gains: Abaqus L 320.06→179.13 and 314.77→177.55 ms, Nastran L 822.18→655.01 and 811.48→662.90 ms, Netgen L 219.87→136.29 and 219.46→136.82 ms, FLAC3D L 277.26→147.94 and 272.31→147.93 ms. All four reduce allocation calls and cumulative requested bytes in every configuration. SEQ L requested bytes fall from 1,068,508,549→657,547,999 (Abaqus), 2,432,190,351→1,703,386,391 (Nastran), 792,157,678→473,771,897 (Netgen) and 782,706,085→201,687,829 (FLAC3D); these are cumulative requests, not peak RSS. Unchanged-writer noise remains published, including reversed-order Abaqus L writes 131.75→164.42 ms and Netgen L 265.34→289.25 ms; no writer speedup is claimed.

Focused tests pass 346 Python cases per SEQ/OpenMP/TBB module and 60 native cases per MESHIO/NATIVE/KRATOS backend. The archived baseline also passes 187 ownership/compatibility cases, including long native prefixes, source replacement, 300 mixed-format Nastran joins, nested Abaqus includes and Netgen gzip. The native Abaqus include-local/global id behavior predates this batch and differs from the Python reference for the cross-include set in that test; its expected data is pinned directly rather than claiming new reference compatibility. Gzip, binary FLAC3D and richer set/property paths are fixture-covered but not measured by the ASCII mesh rows.

### Source-copy removal and temporary ownership guards

DEX numeric tokens now view the source and use bounded prefix parsing; D/d-exponent normalization uses a short terminated stack buffer or an owned fallback for long tokens. Gmsh time/metadata headers, OBJ trimmed rows and OpenFOAM vector/field text retain existing source views rather than first constructing strings. PLY's ASCII face path names its row owner explicitly. This is **not a dangling-read crash fix**: the preceding `TextStream(std::string&&)` already moved temporary strings into `mOwned`. That owning overload is retained and directly tested with short and long temporaries. A const owning temporary, which cannot move into the owner, is now rejected; the genuinely borrowing `TextCursor` also rejects owning rvalues. Public APIs, installed headers, writer code/bytes and ABI 22 are unchanged.

The [140-row matrix](https://github.com/loumalouomega/meshioplusplus/blob/main/benchmark/text_io_batch11.csv) and [40-row both-order confirmation](https://github.com/loumalouomega/meshioplusplus/blob/main/benchmark/text_io_batch11_seq_confirmation.csv) preserve file/full-mesh digests against the preceding rebuilt modules and across all seven configurations, with matching compiler commands/settings. Archives, the corrected ownership audit and the retained test-key failure/recovery are in `build/text-io-batch11-lifetimes/`. Gmsh timings select ASCII; OBJ, OpenFOAM and DEX also use text paths. PLY timing uses its **default binary** writer/reader, not an ASCII performance claim; ASCII faces are exercised by its fixtures and source-replacement tests. OpenFOAM field-copy removal is fixture-covered, not timed by the geometry-only mesh input.

| SEQ reader | M read before → after | L read before → after | M allocation calls before → after | L allocation calls before → after |
| --- | ---: | ---: | ---: | ---: |
| Gmsh 4.1 ASCII | 18.4 → 21.4 ms | 70.2 → 76.6 ms | 56 → 56 | 45 → 45 |
| OBJ | 2.60 → 2.60 ms | 6.33 → 6.39 ms | 64,100 → 44,180 | 161,286 → 108,986 |
| PLY binary | 2.01 → 1.97 ms | 4.59 → 4.97 ms | 132,377 → 132,377 | 326,779 → 326,779 |
| OpenFOAM geometry | 212.7 → 225.4 ms | 874.5 → 880.6 ms | 8,565,369 → 8,565,369 | 33,499,026 → 33,499,026 |
| DEX | 8.98 → 7.89 ms | 36.75 → 30.09 ms | 186,570 → 140,042 | 702,414 → 526,926 |

DEX gains reproduce in both SEQ orders (M 8.57→6.99 and 8.46→7.00 ms, L 35.21→29.82 and 34.89→29.30 ms). OBJ L reads are 6.07→5.77 and 6.07→5.74 ms; other copy-removal paths are neutral or fluctuating, not reader speedup claims. The initial Gmsh slowdown is much smaller in confirmation (M 19.47→20.20 and 19.34→19.72 ms, L 69.01→68.90 and 69.12→69.51 ms); PLY M changes direction between orders (2.02→1.85 and 1.99→2.28 ms). All observed timing rows remain published. In particular, unchanged OBJ M writes increase in both orders (6.14→6.52 and 5.78→6.64 ms), an unfavorable observation requiring follow-up rather than a broad no-writer-regression claim. DEX and OBJ reduce cumulative allocation requests in all configurations; SEQ L bytes fall 52,300,098→46,859,989 and 10,825,536→9,407,320 respectively. Gmsh, PLY binary and OpenFOAM geometry allocations/bytes are unchanged.

The focused gates pass 615 Python tests with one optional skip per SEQ/OpenMP/TBB module, 78 selected native cases per MESHIO/NATIVE/KRATOS backend (one external `checkMesh` test skipped), and 201 source-ownership cases on the archived baseline SEQ module. Seven standalone cursor parity/lifetime tests pass AddressSanitizer and UndefinedBehaviorSanitizer. DEX long exponents, signed zero and embedded-NUL prefixes preserve native semantics; source replacement checks both mapped and buffered reads in fresh processes.

### Dtype-hoisted OFF, IP, FLUX and PERMAS writers

Roadmap §3.1.2's first writer batch. The OFF, IP, FLUX and PERMAS writers build one core-private `detail::DoubleView`/`detail::Int64View` per point, connectivity or reference array before the row loop, instead of calling `read_double`/`read_int` (a `dispatch_dtype` switch) per element. The view points into the array when its dtype already is float64/int64 and otherwise makes one converted copy, element for element as `read_double`/`read_int` convert it, so the written bytes cannot change. Readers are untouched. The change lives in non-inline `.cpp` bodies and a private header: no installed header, signature or layout changes, and C++ ABI 22 is unchanged.

The [112-row matrix](https://github.com/loumalouomega/meshioplusplus/blob/main/benchmark/text_io_batch12_float64_int64.csv) and [32-row two-round confirmation](https://github.com/loumalouomega/meshioplusplus/blob/main/benchmark/text_io_batch12_float64_int64_seq_confirmation.csv) use the harness's canonical float64/int64 inputs, the zero-copy case, through strict Python dispatch. File and full-mesh digests agree across stages and all seven configurations. Reader and writer allocation calls and requested bytes are identical before and after in every configuration: canonical inputs allocate no view storage.

| SEQ writer | M write, round 1 | M write, round 2 | L write, round 1 | L write, round 2 |
| --- | ---: | ---: | ---: | ---: |
| OFF | 6.64 → 6.29 ms | 6.98 → 6.88 ms | 13.07 → 13.87 ms | 14.23 → 13.21 ms |
| IP | 4.09 → 4.25 ms | 4.60 → 4.51 ms | 9.92 → 10.50 ms | 10.56 → 9.87 ms |
| FLUX | 229.32 → 224.40 ms | 231.05 → 225.72 ms | 850.81 → 860.19 ms | 895.45 → 863.10 ms |
| PERMAS | 73.23 → 74.44 ms | 81.27 → 74.05 ms | 287.84 → 300.61 ms | 302.76 → 275.20 ms |

These timings are **neutral**: every L row changes direction between rounds, and the writers spend their time formatting numbers, not dispatching on the dtype. No writer speedup is claimed for canonical inputs. The unchanged readers fluctuate the same way (PERMAS L 150.38→158.57 and 153.02→155.99 ms; FLUX L 293.26→301.37 and 389.78→287.40 ms), and all rows remain published. The hoist pays off for non-canonical inputs (float32 points, int32 connectivity), where the per-element switch is replaced by one parallel conversion. `test_text_io_dtypes.cpp` and `test_text_io_dtypes.py` pin byte identity across all ten dtypes and in 2-D and 3-D, and check that a caller's arrays are not modified.

### Dtype-hoisted Gmsh, FEBio, Z88, Patran, Femap and MDPA writers

Roadmap §3.1.2's second writer batch, in the same shape as the [first](#dtype-hoisted-off-ip-flux-and-permas-writers): the Gmsh 2.2/4.1 (ASCII and binary), FEBio, Z88, Patran, Femap and MDPA writers read their points, connectivity, per-cell tags/types/properties and data arrays through one private `detail::DoubleView`/`detail::Int64View` per array, instead of a `read_double`/`read_int` dtype switch per element. Z88's constraint and surface-load arrays, MDPA's per-value formatter (now a small array-level helper) and Femap's mesh fingerprint use views too. Readers are untouched. The change is in non-inline `.cpp` bodies and a private header, so installed headers and C++ ABI 22 are unchanged.

The [168-row matrix](https://github.com/loumalouomega/meshioplusplus/blob/main/benchmark/text_io_batch13_float64_int64.csv) (SEQ at one thread; OpenMP and TBB at 1/4/8) and the [48-row two-round confirmation](https://github.com/loumalouomega/meshioplusplus/blob/main/benchmark/text_io_batch13_float64_int64_seq_confirmation.csv) use the canonical float64/int64 benchmark inputs, where a view is zero-copy. Gmsh timing selects its ASCII 4.1 path; the binary paths and the 2.2 writer are covered by the dtype tests. Written-file and parsed-mesh SHA-256 values agree across stages, backends and thread counts, and reader and writer allocation calls and requested bytes are identical before and after in every configuration: canonical inputs allocate no view storage.

| SEQ writer | M write, round 1 | M write, round 2 | L write, round 1 | L write, round 2 |
| --- | ---: | ---: | ---: | ---: |
| Gmsh 4.1 ASCII | 101.3 → 99.9 ms | 100.3 → 91.7 ms | 387.5 → 386.1 ms | 383.9 → 384.4 ms |
| FEBio | 79.4 → 75.7 ms | 81.3 → 76.4 ms | 312.2 → 300.6 ms | 310.0 → 305.2 ms |
| Z88 | 170.1 → 175.6 ms | 172.0 → 172.7 ms | 703.3 → 687.6 ms | 674.2 → 675.6 ms |
| Patran | 472.5 → 482.2 ms | 480.7 → 479.4 ms | 1859.0 → 1892.4 ms | 1873.5 → 1831.7 ms |
| Femap | 201.4 → 224.2 ms | 208.3 → 210.5 ms | 799.4 → 846.0 ms | 862.7 → 833.0 ms |
| MDPA | 136.9 → 138.2 ms | 133.6 → 127.4 ms | 514.7 → 525.9 ms | 486.3 → 514.2 ms |

These timings are **neutral**: FEBio writes are the only ones that improve in all four cells (about 1.5–5%); the others change direction between rounds, and Femap's and MDPA's first-round increases are not reproduced in round 2 but are retained. The writers spend their time formatting numbers, so no canonical-input speedup is claimed. The unchanged readers fluctuate by more than the writers do (Femap L 1850.3 → 1472.8 and 1805.0 → 1784.9 ms; Z88 L 639.2 → 761.4 and 750.1 → 764.5 ms), which bounds how much any one cell can be trusted. All rows remain published. The hoist matters for non-canonical inputs (float32 points, int32 connectivity), where the per-element switch becomes one parallel conversion.

`test_text_io_dtypes.cpp` and `test_text_io_dtypes.py` now pin byte identity against canonical storage across all ten dtypes, in 2-D and 3-D, for every writer in both batches, and check that the caller's arrays are not modified. NATIVE and KRATOS ingest arrays as canonical float64/int64, so on those backends the check exercises the same path through converted copies.

Measured on this session's container (4 cores, GCC 13.3 with `-O3 -DNDEBUG -ffp-contract=off`, SEQ/OpenMP/TBB modules built from the same tree), seven repeats after one warmup, TBB capped with `tbb::global_control` and OpenMP with `OMP_NUM_THREADS`; absolute times are not comparable with the earlier batches' machine. The first matrix hashed MDPA's parsed `properties_0` dictionary through its pointer bytes, which varied between processes while every written file was identical; the digest now hashes the dictionary's text and MDPA was re-measured, so the published MDPA rows come from that second run and the other formats from the first.

### Dtype-hoisted CGNS, MED, libMesh, PCD, GiD and glTF

Roadmap §3.1.2's third batch. These formats read their points, connectivity, offsets, tags and data through one private `detail::DoubleView`/`detail::Int64View` per array instead of a `read_double`/`read_int` dtype switch per element, in the writers and in the per-value loops of the CGNS, MED and PCD readers (MED's polyhedron and polygon offset tables, CGNS coordinate and connectivity columns, PCD's coordinate and normal columns). libMesh builds its per-block `p_level` views on first use, so only the blocks of written cells are touched; PCD's ASCII writer holds one view per output column; glTF's colour, attribute, normal and region loops are hoisted, and its per-cell `gltf_cell_ids` keeps its single read because a view per cell would cost more than the switch. The changes are in non-inline `.cpp` bodies and a private header, so installed headers and C++ ABI 22 are unchanged.

The [140-row matrix](https://github.com/loumalouomega/meshioplusplus/blob/main/benchmark/text_io_batch14_float64_int64.csv) (SEQ at one thread; OpenMP and TBB at 1/4/8) and the [40-row two-round confirmation](https://github.com/loumalouomega/meshioplusplus/blob/main/benchmark/text_io_batch14_float64_int64_seq_confirmation.csv) use the canonical float64/int64 inputs through strict Python dispatch, on the HDF5-enabled modules. CGNS and MED outputs embed per-write HDF5 metadata that varies between processes, so their `file_sha256` is the digest of every dataset's name, dtype, shape and bytes; every other format hashes the file bytes. All digests, and the parsed-mesh digests, agree across stages, backends and thread counts, and reader and writer allocation calls and requested bytes are identical before and after in every configuration. glTF is write-only, so the round-trip harness cannot time it; it is covered by the dtype tests only.

| SEQ | M write, rounds 1 / 2 | L write, rounds 1 / 2 | M read, rounds 1 / 2 | L read, rounds 1 / 2 |
| --- | ---: | ---: | ---: | ---: |
| CGNS | 51.1 → 53.1 / 52.6 → 43.6 ms | 211.3 → 170.7 / 220.4 → 220.9 ms | 14.6 → 14.7 / 14.7 → 11.8 ms | 57.9 → 56.6 / 77.6 → 47.6 ms |
| MED | 8.5 → 8.5 / 8.5 → 8.4 ms | 35.7 → 35.5 / 59.7 → 33.7 ms | 2.4 → 2.4 / 2.4 → 2.6 ms | 7.1 → 7.6 / 28.4 → 8.1 ms |
| libMesh | 101.1 → 95.9 / 93.0 → 94.8 ms | 376.7 → 352.1 / 408.7 → 413.3 ms | 137.2 → 157.1 / 173.7 → 161.4 ms | 620.9 → 708.8 / 633.8 → 716.7 ms |
| PCD | 0.5 → 0.7 / 0.5 → 0.7 ms | 0.9 → 0.8 / 1.0 → 0.7 ms | 0.2 → 0.2 / 0.2 → 0.2 ms | 0.5 → 0.5 / 0.5 → 0.4 ms |
| GiD | 120.0 → 119.1 / 124.0 → 129.0 ms | 482.8 → 491.5 / 470.9 → 495.6 ms | 201.9 → 194.3 / 194.4 → 198.3 ms | 902.5 → 898.7 / 872.0 → 843.5 ms |

These timings are **neutral**, and two observations are unfavorable and retained rather than explained away. libMesh L **reads** are about 12% slower in both orders (620.9 → 708.8 and 633.8 → 716.7 ms) although the libMesh reader code is untouched in this batch, so this is not attributable to the hoist and remains open for follow-up (code layout or another cause has not been ruled out). PCD M **writes** rise by 0.2 ms in both rounds (0.5 → 0.7 ms), while L changes direction; at sub-millisecond scale this is retained as an observation, and no PCD speedup is claimed. The CGNS and MED rows include 2–4x outliers in single cells (MED L reads 28.4 ms, CGNS L reads 77.6 ms) that appear only in the second round's before run; their cause was not isolated, and they are the reason no CGNS or MED claim is made. The hoist pays off for non-canonical inputs, where the per-element switch becomes one parallel conversion.

A first measurement of this batch caught two allocation regressions before publication, which were fixed rather than documented as expected: the PCD reader made a converted double copy of every float32 coordinate column (+3 allocation calls, up to +436 KB requested), and the CGNS reader's `std::deque` of coordinate views cost two allocations. PCD now converts each column in place with one dtype dispatch per column and CGNS uses a fixed array; the published matrix is the rerun, and the same harness change (dataset digests for the HDF5 outputs) applies to both stages. `test_text_io_dtypes.cpp` and `test_text_io_dtypes.py` pin libMesh, GiD and glTF byte for byte against canonical storage across all ten dtypes, and PCD, CGNS and MED, which record each array's dtype in the file, by round-tripped values; every writer is also checked to leave the caller's arrays untouched.

Measured in the same container and with the same build settings as the [previous batch](#dtype-hoisted-gmsh-febio-z88-patran-femap-and-mdpa-writers), with HDF5 enabled and NetCDF disabled, seven repeats after one warmup.

### Dtype-hoisted OBJ, STL, PLY, TetGen, UGRID, DOLFIN, FreeFEM, AVS-UCD, WKT, Triangle, SVG and TikZ writers

Roadmap §3.1.2's fourth batch, the simple text writers. They read points, connectivity, tags and data through one private `detail::DoubleView`/`detail::Int64View` per array instead of a `read_double`/`read_int` dtype switch per element. Arrays that were looked up by name or block inside a row loop (Triangle, TetGen, AVS-UCD, FreeFEM) get one view before the loop; an array absent for a block keeps its original fallback value, so no view is made for it. The PLY writer holds one view per scalar point property, of the kind its dtype is written as; the PLY reader's per-column loop and scalar header read are unchanged. The changes are in non-inline `.cpp` bodies and a private header, so installed headers and C++ ABI 22 are unchanged.

The [252-row matrix](https://github.com/loumalouomega/meshioplusplus/blob/main/benchmark/text_io_batch15_float64_int64.csv) (SEQ at one thread; OpenMP and TBB at 1/4/8) and the [72-row two-round confirmation](https://github.com/loumalouomega/meshioplusplus/blob/main/benchmark/text_io_batch15_float64_int64_seq_confirmation.csv) cover the nine formats the round-trip harness can time, on the canonical float64/int64 inputs. Triangle writes 2-D points only and the harness inputs are 3-D, and SVG and TikZ are write-only, so those three are covered by the dtype tests alone. Written-file and parsed-mesh SHA-256 values agree across stages, backends and thread counts, and reader and writer allocation calls and requested bytes are identical before and after in every configuration. The before modules are built from `d545e00`, and the after modules from the same tree plus the batch's diff.

| SEQ | M write, rounds 1 / 2 | L write, rounds 1 / 2 | M read, rounds 1 / 2 | L read, rounds 1 / 2 |
| --- | ---: | ---: | ---: | ---: |
| AVS-UCD | 145.33 → 143.51 / 146.99 → 148.27 ms | 568.63 → 554.35 / 558.62 → 539.13 ms | 85.16 → 87.23 / 80.03 → 87.66 ms | 337.33 → 334.77 / 328.72 → 335.81 ms |
| DOLFIN | 240.56 → 239.58 / 239.59 → 235.01 ms | 937.75 → 934.54 / 926.88 → 936.87 ms | 211.29 → 209.22 / 211.74 → 205.76 ms | 836.79 → 819.79 / 816.98 → 803.70 ms |
| FreeFEM | 119.80 → 122.17 / 128.09 → 117.68 ms | 468.92 → 476.48 / 466.23 → 467.15 ms | 65.68 → 64.26 / 66.08 → 65.03 ms | 265.27 → 256.29 / 263.18 → 254.60 ms |
| OBJ | 8.16 → 8.30 / 8.07 → 7.97 ms | 19.71 → 20.58 / 18.95 → 19.71 ms | 5.06 → 5.23 / 5.06 → 5.29 ms | 12.63 → 13.62 / 12.47 → 13.10 ms |
| PLY | 2.21 → 2.16 / 2.30 → 2.10 ms | 4.82 → 5.29 / 5.09 → 4.85 ms | 4.01 → 3.95 / 3.84 → 3.94 ms | 9.62 → 9.85 / 9.40 → 9.87 ms |
| STL | 47.31 → 49.84 / 54.56 → 47.65 ms | 129.59 → 120.00 / 119.78 → 118.93 ms | 22.93 → 23.02 / 23.00 → 22.81 ms | 60.44 → 57.93 / 56.45 → 57.31 ms |
| TetGen | 125.03 → 122.04 / 129.69 → 118.69 ms | 480.00 → 473.45 / 483.45 → 469.04 ms | 110.68 → 116.20 / 116.88 → 114.66 ms | 518.14 → 518.21 / 500.29 → 502.56 ms |
| UGRID | 96.07 → 96.28 / 101.01 → 95.72 ms | 377.36 → 362.46 / 370.52 → 363.86 ms | 59.16 → 59.42 / 62.89 → 60.25 ms | 237.63 → 239.00 / 236.66 → 240.33 ms |
| WKT | 44.50 → 42.76 / 44.07 → 43.32 ms | 107.75 → 109.45 / 109.27 → 109.30 ms | 25.23 → 23.49 / 23.44 → 23.21 ms | 59.75 → 60.75 / 59.57 → 58.55 ms |

These timings are **neutral**. Most writer cells move by a few percent in either direction, and the writers spend their time formatting numbers, so no canonical-input speedup is claimed. Two retained observations point the same way in both rounds. OBJ L writes are about 4% slower (19.71 → 20.58 and 18.95 → 19.71 ms); this is the writer follow-up promised by the [OBJ M-write observation](#source-copy-removal-and-temporary-ownership-guards), whose M writes are now neutral (8.16 → 8.30 and 8.07 → 7.97 ms), so that earlier increase is not reproduced. The unchanged OBJ reader moves by a similar amount in the same direction (L 12.63 → 13.62 and 12.47 → 13.10 ms), which suggests a shift common to the module rather than a writer effect, but this was not isolated. The hoist pays off for non-canonical inputs, where the per-element switch becomes one parallel conversion.

`test_text_io_dtypes.cpp` and `test_text_io_dtypes.py` extend the earlier batches to these writers. Every file a writer produces, including companions such as TetGen's `.ele` and DOLFIN's mesh-function files, is compared with canonical storage across all ten dtypes. UGRID (integer labels), AVS-UCD (integer materials) and DOLFIN (`float` versus `int` mesh functions) choose their output by whether an array is float or integer, so they are compared with the same class in float64 or int64 rather than always float64/int64. PLY records property dtypes in its header and is checked by round-tripped values.

## Every format

`benchmark/bench.py` also times a write and a read of **every** format meshio++ both writes and reads back, each fed the largest input its [conformance declaration](./conformance.md) says it keeps: the synthetic tetrahedral cube for volume formats, its surface for surface formats (STL, OBJ, PLY, …), its points for point clouds.

Explicitly tested single-type exceptions also participate: MFM and Gmsh 4.1 reject the mixed conformance mesh but round-trip the harness's tetrahedral cube. This does not change their conformance declarations or imply mixed-cell support. Use `--formats mfm,gmsh` to measure those supported reader/writer paths.

```sh
python benchmark/bench.py --sizes S,M --out results_all.csv            # every format
python benchmark/bench.py --sizes L --formats vtu,gmsh22,xdmf,med      # a subset
python benchmark/bench.py --sizes M,L --formats mfm --repeats 7      # single-type MFM
```

The sizes are 6·(n−1)³ tetrahedra for n = 16, 36 and 56 points per edge (about 20k, 250k and 1M). Legacy meshio is optional here: with `MESHIO_LEGACY_SRC` pointing at a source checkout's `src`, or `meshio` installed, the curated comparison above fills its legacy columns; without it only the meshio++ columns are written.

### Text I/O: MFM baseline

The first prerequisite for roadmap §3.1.1 was MFM benchmark coverage. Its original native reader stored the remaining file in owning token strings, including the discarded reference arrays. The [before/after evidence](https://github.com/loumalouomega/meshioplusplus/blob/main/benchmark/text_io_mfm.csv) preserves all 14 baseline rows and adds 14 after rows on MESHIO: SEQ at one thread and OpenMP/TBB at 1/4/8 threads. Written-file SHA-256 values and parsed-mesh SHA-256 values (points, block types/connectivity and `mfm:ref`, including array shapes and dtypes) agree across stages and all configurations; parsed arrays also match the generated input, with the writer's default subdomain value of one.

Recorded on 2026-10-03 against native source commit `3b0791205bd68cdd50d3127d1aa9d08047ba9900`, AMD Ryzen 7 255, Release/MESHIO, HDF5/netCDF/zlib enabled. SEQ/TBB use Clang 23.1.1 and OpenMP uses GCC 16.2.1, all with `-O3 -DNDEBUG -ffp-contract=off`. Timings use the harness's one warmup and median of seven runs, without an allocation interposer or concurrent project builds/tests. TBB is explicitly capped with `tbb::global_control`; OpenMP uses `OMP_NUM_THREADS` with dynamic adjustment disabled. Compare each configuration with its own later result, not compiler-to-compiler timings.

| SEQ input | Cells | Read median | Write median | Reader allocation calls | Cumulative requested bytes |
| --- | ---: | ---: | ---: | ---: | ---: |
| M tetrahedral grid | 257,250 | 166.6 ms | 188.2 ms | 140,031 | 551,509,791 (525.96 MiB) |
| L tetrahedral grid | 998,250 | 576.3 ms | 523.6 ms | 526,913 | 2,203,755,809 (2101.67 MiB) |

Allocation counting is separate: ordinary global `operator new` requests around a warmed-up direct native read, with three identical count/byte results per configuration and size. These are cumulative requests, **not peak RSS**; pathname lengths can affect small allocations, so the same input path is retained before and after. The local `build/text-io-mfm/{before,after}/` archives preserve copied native modules with SHA-256 checks, source/helper snapshots, actual compile commands, CMake caches and raw JSON/logs; editable auto-rebuilding is disabled and each measurement explicitly loads its archived module. The after source is anchored at `2d39e4020d2cac0a8372e170eb9ce92296741614` plus the preserved MFM implementation diff; compiler versions and actual MFM compile commands match the baseline.

### MFM: bounded token views

The reader now uses `FileSource`, `TextStream` for the header and blank-separated token views for the body. Bounded prefix helpers retain its lenient C-library number semantics; discarded reference tokens and trailing tokens remain ignored. Mesh arrays own their storage after the source is released. The writer is unchanged. This changes an exported non-inline implementation only: no installed headers, signatures, layouts or default arguments change, C++ ABI 22 stays unchanged, and the single header is regenerated.

| SEQ input | Original read median | Token-view read median | Reader allocation calls after | Cumulative requested bytes after |
| --- | ---: | ---: | ---: | ---: |
| M tetrahedral grid | 166.6 ms | 39.1 ms (4.27×) | 57 | 279,846,415 (266.88 MiB) |
| L tetrahedral grid | 576.3 ms | 248.7 ms (2.32×) | 59 | 1,117,887,823 (1066.10 MiB) |

The unchanged writer showed noise in the initial sequential sweep (SEQ M/L 188.2/523.6 → 198.4/591.6 ms), so a separate [uninstrumented SEQ confirmation](https://github.com/loumalouomega/meshioplusplus/blob/main/benchmark/text_io_mfm_seq_confirmation.csv) retains both measurement orders, before→after and after→before, with seven repetitions each. M writer medians are 192.1→195.8 and 193.3→194.1 ms; L medians are 536.0→533.9 and 540.9→534.9 ms. No material writer slowdown reproduces. Confirmation read medians are 166.6→38.5 and 171.1→38.6 ms for M, and 742.6→154.8 and 680.3→153.8 ms for L; the original sweep remains in the evidence rather than being replaced by more favorable timings. No parallel speedup or peak-RSS reduction is claimed.

Direct native tests on MESHIO/NATIVE/KRATOS and Python tests cover every supported linear type, native/reference geometry parity, mixed-type write rejection, blank/CRLF input, missing final newline, empty meshes, section truncation, existing lenient prefixes and ownership after removing/overwriting the input. The shared tokenizer and locale guards, plus Gmsh's native/reference ASCII byte comparisons, remain the gate. MFM is removed from the roadmap's remaining stream-reader list; the other readers and dtype-hoisting work remain open.

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

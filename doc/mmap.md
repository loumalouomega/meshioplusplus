# Memory-mapped reading

Readers that consume a whole file at once can memory-map it instead of copying it into a heap buffer. This is on by default for large files and needs no code change.

## What it actually buys

Measured against one bulk `read()` into a pre-sized buffer — which is what these readers already did — mapping saves **one full-file `memcpy` and the second copy of the file in RAM**. It does *not* make the page faults go away.

So this is a **memory-footprint feature more than a throughput one**. On a multi-GB mesh it roughly halves peak RSS during the read. The clear exception is a reader that was building an `ostringstream` and calling `.str()` — that pays for two extra copies, and mapping removes both.

If you were expecting a large wall-clock win, this is not that feature; selective reads (see [Selective reads](selective_read.md)) are where the read-time savings are.

## When it engages

| `MmapMode` | Behaviour |
| --- | --- |
| `Auto` (default) | map regular files ≥ 16 MiB |
| `On` | map whenever possible |
| `Off` | never map |

Mapping is **advisory**. A pipe, a zero-length file, a platform without support, WASM, or any failure in the mapping call itself falls back silently to a buffered read. It never turns a readable file into an error.

`MESHIOPLUSPLUS_MMAP_THRESHOLD` (bytes) overrides the `Auto` threshold, mainly for benchmarking; `0` maps everything.

### One deliberate exception

A mapped buffer behaves as if NUL-terminated only because the kernel zero-fills the remainder of the final page — which the C string functions some parsers use (`strtod`) quietly rely on. That slack does **not** exist when the file size is an exact multiple of the page size, so those files are never mapped, even under `On`. It costs a buffered read roughly once in 4096 files and buys back a hard guarantee.

## Coverage

Every whole-file reader now goes through `FileSource`: **gmsh** (which additionally honours an explicit `ReadOptions::mMmap`), **vtk**, **ensight**, **ugrid**'s ASCII branch and **openfoam**.

`openfoam` gained the most: it previously slurped via `ostringstream` + `.str()`, paying for **two** extra full-file copies on top of the read, both of which are now gone.

Streaming readers are unaffected by design — `ugrid`'s binary branch reads each section straight into its destination array with no whole-file intermediate, so there is nothing to map, and the HDF5/netCDF-backed formats manage their own I/O.

## Safety

Nothing in a returned `Mesh` ever points into the mapping: every reader copies what it parses into owning storage through the uniform mesh API. The mapping is therefore a function-local RAII object released when the reader returns — no keep-alive, no shared ownership, no capsule.

A test reads a file, then deletes **and overwrites** it, and re-validates the mesh, so a violation of that rule fails loudly in CI rather than becoming a use-after-free downstream.

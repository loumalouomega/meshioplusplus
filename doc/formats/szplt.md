# Tecplot SZL (`.szplt`)

Tecplot's subzone-loadable format, the default output of Tecplot 360 since 2016 and of TecIO-based solvers (SU2, FUN3D…). It is undocumented, and TecIO, Tecplot's own library, is its only reader, so meshio++ reads it through a TecIO the user already has. New in v16.13.0; written too since the writer joined.

| | |
|---|---|
| **Format name** | `szplt` |
| **Extensions** | `.szplt` (also found by content, `#!SZPLT`) |
| **Read / Write** | ✓ / ✓ ([conformance](../conformance.md#szplt)) |
| **Time steps** | ✓ — as [Tecplot](./tecplot.md): the distinct solution times |
| **Extra dependencies** | TecIO: a core built with `MESHIOPLUSPLUS_WITH_TECIO=ON`, or a shared TecIO named by `MESHIOPLUSPLUS_TECIO_LIBRARY` |

## Reading

```python
import meshioplusplus

mesh = meshioplusplus.read("flow.szplt")
last = meshioplusplus.read("flow.szplt", time_step=-1)
meshioplusplus.szplt.time_values("flow.szplt")
```

A `.szplt` reads exactly as the same data saved as `.plt` or `.dat` does: TecIO's reader fills the zone model the ASCII and binary Tecplot readers share, so zones, regions, point and cell data, variable sharing, passive variables, shared connectivity and time steps follow [Tecplot](./tecplot.md#data-mapping). FEPOLYGON and FEPOLYHEDRON zones are refused by name (save such a file as `.plt`).

## Writing

```python
meshioplusplus.write("flow.szplt", mesh)                                   # a static data set
meshioplusplus.write_sequence("flow.szplt", steps, file_format="szplt")   # a strand of zones per step
```

```bash
meshioplusplus convert model.vtu model.szplt
```

TecIO writes the file, in the layout the [ASCII and `.plt` writers](./tecplot.md#reading-writing) use: one zone per cell block (named after a cell region that is exactly the block, else `block_<i>`), X, Y, (Z) and the point data as nodal variables, cell data as cell-centred variables (a variable a block has no values for is passive in its zone), vector and tensor arrays split into `<name>_<k>`, and doubles throughout. A wedge or pyramid is a degenerate brick and comes back as a hexahedron, as in `.plt`; vertices and polygon or polyhedron blocks are skipped with a warning. TecIO shares a variable between zones only together with their connectivity, so the blocks of one mesh cannot share a point array as they do in the other Tecplot writers: a mesh of one block is one zone with every point, and a mesh of several gives each zone the points its cells use (a point several blocks use is repeated, one no cell uses is dropped), as TecIO's own files have it. A series is a strand (`STRANDID` 1) of zones per step at the step's time, the later steps sharing the first step's coordinates and connectivity, so every step needs the first step's points and cells and the same variables. The core built with TecIO and the Python twin through `MESHIOPLUSPLUS_TECIO_LIBRARY` write the same zones; they are never loaded together, since the twin needs a second, shared TecIO.

## Getting TecIO

meshio++ never redistributes TecIO. Two ways to point it at yours:

- **A native build** with `-DMESHIOPLUSPLUS_WITH_TECIO=ON -DTECIO_ROOT=<dir>`, where `<dir>` holds `TECIO.h` (with `tecio_Exports.h` and `StandardIntegralTypes.h`) and `libtecio` — a static `libtecio.a` also needs `pthread`, which the build adds. This is what the C, Fortran, Julia and R bindings and the native CLI read through. TecIO is the source Tecplot distributes on its [TecIO page](https://tecplot.com/products/tecio-library/) (free, a registration form), or the copy a Tecplot 360 install carries; `build/configure.sh --with-tecio <dir>` sets both flags.
- **The Python reader** loads a *shared* TecIO through ctypes: set `MESHIOPLUSPLUS_TECIO_LIBRARY` to its path (Tecplot 360 ships `libtecio.so` / `tecio.dll` in its `bin` directory), or have it on the library search path. TecIO's own CMake builds a static library with hidden symbols; for a shared one, build it with `add_library(tecio SHARED …)` and without `-fvisibility=hidden`.

Without either, reading raises an `ImportError` naming both and the third way: save the file as `.plt` in Tecplot ([Tecplot route](../routes/tecplot_szplt.md)).

## Validation

The fixtures (`tests/python/meshes/szplt`) are written by TecIO 2018.3 — the TecIO source SU2 vendors — each beside an ASCII twin holding the same data: two FE zones of tetrahedra and bricks with cell-centred data, an ordered zone, a transient file sharing variables and connectivity with passive variables, and single-precision triangles. Both engines read every `.szplt` identically to its twin. A `.szplt` written by Tecplot 360 itself is listed in the [roadmap](../roadmap.md#awaiting-a-licensed-run).

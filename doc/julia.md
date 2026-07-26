# Julia

meshio++ ships a Julia package, `MeshioPlusPlus`, layered on the [C API](/c_api) via `ccall` — the same way the [Fortran](/fortran) module is, and aimed at the same HPC audience:

```julia
import MeshioPlusPlus as mio
using MeshioPlusPlus

m = mio.read("bracket.msh")
println(m)                     # Mesh(9231 points, 42145 cells in 3 blocks)
surf = extract_surface(m)
mio.write(surf, "bracket_surface.vtu")
```

::: warning This binding is not MIT
`bindings/julia/` is released under the **GNU General Public License, version 3 (GPL-3.0)**, not the MIT licence covering the rest of meshio++. GPL-3.0 is a **copyleft** license, not a permission-required one: anyone — including a company — may use, modify or sell it commercially with **no permission needed**, but *distributing* it or a modified version of it must be under GPL-3.0 too, with source available. Purely private/internal use that is never distributed carries no obligation.

The C API and C++ core this binding calls are unaffected and stay MIT. Calling that stable, non-GPL C ABI via `ccall`/`dlopen` at runtime is the standard "linking exception" case — it does not require the C library to also be GPL. See [`bindings/julia/LICENSE`](https://github.com/loumalouomega/meshioplusplus/blob/master/bindings/julia/LICENSE).
:::

## Building and installing

The package binds the **installed** C library; no C++ is compiled by it.

```sh
./build/configure.sh --c-api --build
cmake --install build/cpp-release --prefix /opt/meshioplusplus
```

Then either put `/opt/meshioplusplus/lib` on the loader path, or point the package straight at the library:

```julia
ENV["MESHIOPLUSPLUS_LIB"] = "/opt/meshioplusplus/lib/libmeshioplusplus.so"

using Pkg
Pkg.develop(path="bindings/julia/MeshioPlusPlus")
```

`MESHIOPLUSPLUS_LIB` is checked first, then the standard loader path (`Libdl.find_library`). Failing both raises an error naming these two build commands.

::: warning A build with HDF5 can fail to load from Julia on Debian/Ubuntu
If the C API was built with `-DMESHIOPLUSPLUS_WITH_HDF5=ON`, loading it from a Julia process can fail with `libcurl.so.4: version 'CURL_OPENSSL_4' not found`. This is a real Debian/Ubuntu + Julia interaction, not a meshio++ bug: `libhdf5-dev` there pulls in `libhdf5_openmpi`, which links the *system* `libcurl`, while Julia bundles its own `libcurl` — and depending on the process's library-resolution order (an IJulia kernel subprocess is more likely to hit it than an interactive `julia` invocation), the two can conflict. Build with `-DMESHIOPLUSPLUS_WITH_HDF5=OFF -DMESHIOPLUSPLUS_WITH_NETCDF=OFF` if this binding is all you need from the library, as the `julia` CI job does.
:::

::: tip No registration yet, and no JLL
GPL-3.0 **is** OSI-approved, so the package is eligible for Julia's General registry — but registering it is a separate follow-up step, not done yet, so `Pkg.add("MeshioPlusPlus")` will not work until then; install by path or by URL in the meantime.

There is also deliberately no [JLL](https://docs.binarybuilder.org/stable/jll/): shipping a binary artifact through BinaryBuilder/Yggdrasil is a real distribution step and belongs in a follow-up rather than being faked here.
:::

## Array layout and 1-based indexing

The C core stores points row-major as `(num_points, dim)` and connectivity as `(num_cells, nodes_per_cell)`. Because Julia is column-major, the **same memory** is naturally

```julia
points(m)                      # (dim, num_points)
connectivity(m, 1)             # (nodes_per_cell, num_cells)
point_data(m, "displacement")  # (components, num_points)
```

so nothing is ever transposed — exactly the reasoning the Fortran module documents. Node indices are **1-based** here, and the ±1 shift happens inside the **copying** accessors, where a copy is made anyway:

| accessor | copies? | node indices |
|---|---|---|
| `points(m)` | yes | n/a — coordinates carry no indices |
| `points_ptr(m)` | **no** (zero-copy borrow) | n/a |
| `connectivity(m, i)` | yes | **1-based** |
| `connectivity_ptr(m, i)` | **no** | **0-based** — the ABI's own |

A borrow cannot be shifted without copying it, which is the whole point of a borrow — hence the two names. `point_data_ptr`, `cell_data_ptr` and `field_data_ptr` follow the same pattern.

Index maps and permutations (`refine`, `convert_cells`, `decimate`, `partition`, `reorder`) are copies, so they are 1-based too, and the C API's `-1` "pruned / absent" sentinel becomes **`0`** — never a valid 1-based index. That is verbatim the Fortran rule; the bindings agree deliberately.

`partition_labels` is the exception: those are part **ids**, not indices, so they stay in `0:nparts-1`.

## The borrow window

`points_ptr` and friends return a `MeshBorrow`, a real zero-copy view of the buffer inside the C++ core (`unsafe_wrap`). The C API's rule is that such a pointer stays valid until the next **mutating** `mio_mesh_*` call on that mesh, or until it is freed; read-only accessors never invalidate it.

`MeshBorrow` enforces exactly that rather than trusting the caller:

```julia
b = points_ptr(m)
b[1, 2]                                # fine
add_point_data!(m, "T", temperatures)  # a MUTATING call ends the window
b[1, 2]                                # BorrowError -- not a stale read
```

It also holds a reference to the owning `Mesh`, so the mesh cannot be garbage-collected while a view of it is alive. `parent(b)` returns the raw `Array` for hot loops; that escape hatch is **unchecked** and valid only inside the window.

## Memory management

A `Mesh` releases its handle through a **finalizer** — the one real difference from the Fortran module, where meshes are freed explicitly with `call m%free()`. `close(m)` releases one deterministically and is idempotent.

Operations producing an opaque C result (`split`, `partition`, `reorder`, `refine`, `decimate`, `convert_cells`) always **transfer ownership** of the mesh out of that result rather than handing back a borrow into it, so a piece stays valid after the result is gone:

```julia
for (key, piece) in mio.split(m; by="type")
    mio.write(piece, "part_$key.vtu")   # still valid; the result is long gone
end
```

## Names that would shadow `Base`

`read`, `write`, `convert`, `merge`, `split` and `diff` are **not exported**, because they collide with `Base`. Call them qualified:

```julia
import MeshioPlusPlus as mio
m = mio.read("bracket.msh")
```

Everything else — the accessors, the setters, `regions`, and the remaining operations — is exported. `close` and `isopen` are proper `Base` method extensions.

## Error handling

Every failure raises a `MeshioError` carrying the C API's own thread-local message; a status code never reaches the caller.

```julia
julia> mio.read("/nope.vtu")
ERROR: MeshioError (read error): meshio++: cannot open file '/nope.vtu'
```

Using a borrow outside its window raises `BorrowError` instead.

## Named regions

```julia
add_region!(m, "inlet", :point, [1, 3, 5])
add_region!(m, "solid", :cell, [1, 2]; dim=3, tag=17)
add_region!(m, "wall", :side, Int64[1 2; 0 2])   # (cell, facet) pairs

for r in regions(m)
    println(r.name, " ", r.kind, " ", size(r.entries, 2), " entries")
end
```

Entries are 1-based, with one exception documented in [`doc/regions.md`](/regions): for a `:side` region the second row is a **facet ordinal within the cell type**, not a mesh index, so it is passed through unshifted — the same rule Fortran's `partition_labels` follows.

## Selective reads and time steps

`ReadOptions` narrows what a read materializes, and since v8.6.0 also picks which time step
of a multi-step file to decode:

```julia
using MeshioPlusPlus

m = MeshioPlusPlus.read("big.vtu"; options = ReadOptions(points_only = true))

# 0 (the default) is the first step; negative counts from the end.
last = MeshioPlusPlus.read("run.exo"; options = ReadOptions(time_step = -1))

meta = read_metadata("run.exo")
meta.time_values          # [0.0, 0.5, 1.0] -- how many steps `time_step` may name
```

Out of range is an error naming the available count, never a silent clamp.
`meta.time_values` is empty for a format with no time concept, so `length(...)` is always
safe. Honoured by `exodus`; see [Selective reads](/selective_read).

## Transient (time-series) XDMF writing

`XdmfSeries` is the write half of the above, and the one writer `write` cannot express: a
series is a **stateful** multi-call object, so it gets its own handle rather than a keyword.
The grid goes out once and each solve appends a cheap step. See
[XDMF time series](/xdmf_time_series).

```julia
using MeshioPlusPlus

s = XdmfSeries("simulation.xdmf")          # "HDF" by default
write_points_cells!(s, mesh)               # the static grid, once
for k in 0:nsteps-1
    solve!(mesh)
    write_data!(s, k * dt, mesh)           # point_data/cell_data only
end
num_steps(s)                               # 10
finalize!(s)                               # close(s) would do this too
close(s)
```

The function form closes the series afterwards even if the body throws:

```julia
XdmfSeries("simulation.xdmf"; data_format = "XML") do s
    write_points_cells!(s, mesh)
    write_data!(s, 0.0, mesh)
end
```

`data_format` is `"HDF"` (the default; needs an HDF5-enabled library), `"XML"` (everything
inline in the `.xdmf`) or `"Binary"`; `gzip_level` applies to `"HDF"` datasets only and is
negative (no compression) by default. An unknown format, or `"HDF"` against a library built
without HDF5, throws a `MeshioError` from the constructor.

Two things worth knowing before reading the result back:

- the `.xdmf` light data is **buffered until the series is finalized**, so the file is only
  readable after `finalize!(s)` (or `close(s)`, which finalizes first);
- `close` swallows a write failure during that implicit finalize — call `finalize!`
  explicitly to see one.

Like `Mesh`, the handle is released by a GC finalizer and `close` is the deterministic,
idempotent form. The name is `finalize!` rather than `finalize` because `Base.finalize`
runs an object's GC finalizer and means something quite different. Reading a finished
series back is the ordinary `MeshioPlusPlus.read(path; options = ReadOptions(time_step = k))`.

## Documented gaps

These are gaps in the **C ABI**, shared with the [Fortran](/fortran) and [R](/r) bindings; the Julia package invents no workaround for any of them:

- **point / cell sets beyond regions** never reach the C++ core at all;
- the **`frozen` pin mask** of `smooth` and `decimate`;
- **per-cell-type counts** in `stats` — use `cell_block_types` with `cell_block_info`;
- **ragged block connectivity** (polygons and polyhedra of varying size): `cell_block_info(m, i).is_ragged` reports them, and `connectivity` then throws rather than returning something wrong;
- the combined **`data_manage`** — `data_drop` / `data_keep` / `data_rename` compose to the same effect;
- **Exodus provenance strings** (`qa_records` / `info_records`): they ride the `ExodusInfo` side channel, which like `MedInfo` does not cross the flat ABI, so `mesh.info` has no counterpart here. Geometry, data, regions and time steps are unaffected.

One further limitation is a *format* one rather than an ABI one, and easy to trip over: **gmsh does not currently round-trip named regions.** The writer emits the `$PhysicalNames` entry but does not attach the physical tag to an entity in `$Elements`, so a reader finds nothing to rebuild the group from. This is pre-existing meshio++ behaviour, reproducible from Python; `abaqus` round-trips regions correctly.

## Tests

```sh
MESHIOPLUSPLUS_LIB=/opt/meshioplusplus/lib/libmeshioplusplus.so \
  julia --project=bindings/julia/MeshioPlusPlus -e 'using Pkg; Pkg.test()'
```

The suite uses the same deliberately non-square fixture as [`tests/fortran/test_fortran_api.f90`](https://github.com/loumalouomega/meshioplusplus/blob/master/tests/fortran/test_fortran_api.f90) — 5 points × 3 dims, 2 tetrahedra × 4 nodes, 3-component vector data — so a transposed mapping or a missed shift cannot cancel out and pass anyway. It pins the column-major identity, the 1-based/0-based accessor pair, the borrow window, regions, and every operation.

## v9.1.0 additions

- `ReadOptions(; lenient=true)` — see [`doc/selective_read.md`](selective_read.md).
- XDMF series: `flush!(s)`, `finalized(s)`, `XdmfSeries(path; mode=:append,
  auto_flush=false)`, and `write_data!(s, t, Dict("u" => values))` for writing a
  step from raw arrays with no `Mesh` in between. An `n x k` matrix is
  transposed on the way out, since Julia is column-major and the C ABI expects
  `k` components per entity row-major.

`flush!` is named with a bang for the same reason `finalize!` is: `Base.flush`
means "flush this IO stream". `MdpaInfo` is not exposed (as for every flat
binding).

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
`bindings/julia/` is released under the **PolyForm Noncommercial License 1.0.0**, not the MIT licence covering the rest of meshio++. Noncommercial use — personal, academic, research, teaching — is free; **use inside a private, proprietary or otherwise commercial project, including internal use at a company, requires the prior explicit written permission** of the copyright holder, Vicente Mataix Ferrándiz (<tote1989@gmail.com>).

The C API and C++ core this binding calls are unaffected and stay MIT. See [`bindings/julia/LICENSE`](https://github.com/loumalouomega/meshioplusplus/blob/master/bindings/julia/LICENSE).
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

::: tip No registry, and no JLL — for now
A non-OSI licence makes the package **ineligible for Julia's General registry**, which requires one, so `Pkg.add("MeshioPlusPlus")` will not work; install by path, by URL, or from a private registry.

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

## Documented gaps

These are gaps in the **C ABI**, shared with the [Fortran](/fortran) and [R](/r) bindings; the Julia package invents no workaround for any of them:

- **point / cell sets beyond regions** never reach the C++ core at all;
- the **`frozen` pin mask** of `smooth` and `decimate`;
- **per-cell-type counts** in `stats` — use `cell_block_types` with `cell_block_info`;
- **ragged block connectivity** (polygons and polyhedra of varying size): `cell_block_info(m, i).is_ragged` reports them, and `connectivity` then throws rather than returning something wrong;
- the combined **`data_manage`** — `data_drop` / `data_keep` / `data_rename` compose to the same effect.

One further limitation is a *format* one rather than an ABI one, and easy to trip over: **gmsh does not currently round-trip named regions.** The writer emits the `$PhysicalNames` entry but does not attach the physical tag to an entity in `$Elements`, so a reader finds nothing to rebuild the group from. This is pre-existing meshio++ behaviour, reproducible from Python; `abaqus` round-trips regions correctly.

## Tests

```sh
MESHIOPLUSPLUS_LIB=/opt/meshioplusplus/lib/libmeshioplusplus.so \
  julia --project=bindings/julia/MeshioPlusPlus -e 'using Pkg; Pkg.test()'
```

The suite uses the same deliberately non-square fixture as [`tests/fortran/test_fortran_api.f90`](https://github.com/loumalouomega/meshioplusplus/blob/master/tests/fortran/test_fortran_api.f90) — 5 points × 3 dims, 2 tetrahedra × 4 nodes, 3-component vector data — so a transposed mapping or a missed shift cannot cancel out and pass anyway. It pins the column-major identity, the 1-based/0-based accessor pair, the borrow window, regions, and every operation.

# meshio++ for Julia

Julia bindings for [meshio++](https://github.com/loumalouomega/meshioplusplus):
read and write ~40 mesh formats used in scientific computing and FEM, and run
meshio++'s operations layer (surface extraction, quality, refine, decimate,
partition, smooth, slice, isosurface, …) on the results.

## Licence — read this first

**This binding is not MIT like the rest of meshio++.** It is released under the
**GNU General Public License, version 3 (GPL-3.0)** — see [`LICENSE`](LICENSE)
for the full, unmodified text.

GPL-3.0 is a **copyleft** license, not a permission-required one:

* anyone — including a company — may use, modify, and even sell this binding
  commercially, **with no permission needed**;
* the condition is on *conveying* (distributing) it: if you distribute this
  binding or a modified version of it, you must do so under GPL-3.0 too, with
  source available to whoever receives it;
* purely private/internal use that is never distributed outside your
  organization carries **no obligation** under GPL-3.0 at all.

The meshio++ C API and C++ core this binding calls are unaffected and remain
MIT under the repository root `LICENSE`. Only the Julia sources in this
directory are GPL-3.0. Calling a stable, non-GPL C ABI via `ccall`/`dlopen` at
runtime — as this binding does — is the standard "linking exception" case: it
does not require the C library on the other side of that ABI to also be
GPL-licensed.

Because GPL-3.0 **is** OSI-approved, this package is eligible for Julia's
General registry (registration is a separate follow-up step, not done yet);
until then, install by path or URL — see below.

## Prerequisite: build and install the C library

This package binds the installed C API, exactly as the Fortran module does. No
C++ is compiled here.

```sh
./build/configure.sh --c-api --build
cmake --install build/cpp-release --prefix /your/prefix
```

Then either put `/your/prefix/lib` on the loader path, or point the package
straight at the library:

```julia
ENV["MESHIOPLUSPLUS_LIB"] = "/your/prefix/lib/libmeshioplusplus.so"
```

`MESHIOPLUSPLUS_LIB` is checked first, then the standard loader path. A
failure to find the library raises an error naming both build commands.

There is deliberately **no JLL**: shipping a binary artifact through
BinaryBuilder/Yggdrasil is a real distribution step, and it is a follow-up
rather than something to fake here.

## Install

```julia
using Pkg
Pkg.develop(path="bindings/julia/MeshioPlusPlus")   # from a checkout
```

## Usage

```julia
import MeshioPlusPlus as mio
using MeshioPlusPlus

m = mio.read("bracket.msh")
println(m)                      # Mesh(9231 points, 42145 cells in 3 blocks)

pts  = points(m)                # (dim, num_points) copy
conn = connectivity(m, 1)       # (nodes_per_cell, num_cells), 1-based copy

surf = extract_surface(m)
mio.write(surf, "bracket_surface.vtu")
```

`read`, `write`, `convert`, `merge`, `split` and `diff` are **not exported**,
because they would shadow the `Base` functions of the same name; call them
qualified as above. Everything else is exported.

## The two conventions

### Column-major — no transpose anywhere

Julia is column-major and the C core row-major, so

```julia
points(m)          # (dim, num_points)
connectivity(m, i) # (nodes_per_cell, num_cells)
point_data(m, "displacement")  # (components, num_points)
```

are the **same memory** as the C API's row-major `(num_points, dim)`,
`(num_cells, nodes_per_cell)` and `(num_points, components)`. Nothing is
transposed on either side. This is the same reasoning the Fortran module
documents.

### 1-based indexing — only where a copy happens

Connectivity is 0-based at the ABI. The ±1 shift happens inside the **copying**
accessors, where a copy is being made anyway:

| accessor | copies? | node indices |
|---|---|---|
| `connectivity(m, i)` | yes | **1-based** |
| `connectivity_ptr(m, i)` | **no** (zero-copy borrow) | **0-based** — the ABI's own |
| `points(m)` | yes | n/a (coordinates carry no indices) |
| `points_ptr(m)` | **no** | n/a |

A borrow cannot be shifted without copying it, which is the whole point of a
borrow — hence the two names.

Index maps and permutations (`refine`, `convert_cells`, `decimate`,
`partition`, `reorder`) are copies, so they are 1-based too, and the C API's
`-1` "pruned / absent" sentinel becomes **`0`** — never a valid 1-based index.
That is exactly the Fortran module's rule; the bindings agree deliberately.

`partition_labels` returns part **ids**, not indices, so its values stay in
`0:nparts-1` and are **not** shifted.

Region entries are 1-based too, with one exception: for a `:side` region the
second row is a **facet ordinal within the cell type**, not a mesh index, so it
is passed through unshifted.

## The borrow window

`points_ptr`, `connectivity_ptr`, `point_data_ptr`, `cell_data_ptr` and
`field_data_ptr` return a `MeshBorrow` — a real zero-copy view of the buffer
inside the C++ core, via `unsafe_wrap`.

The C API's rule is that a borrowed pointer stays valid until the next
**mutating** `mio_mesh_*` call on that mesh, or until the mesh is freed;
read-only accessors never invalidate it. `MeshBorrow` enforces exactly that:

```julia
b = points_ptr(m)
b[1, 2]                                  # fine
add_point_data!(m, "T", temperatures)    # a MUTATING call ends the window
b[1, 2]                                  # BorrowError -- not a stale read
```

It also holds a reference to the owning `Mesh`, so the mesh cannot be
garbage-collected while a view of it is alive. `parent(b)` gives the raw
`Array` for hot loops; that escape hatch is **unchecked** and valid only inside
the window.

## Memory management

A `Mesh` releases its handle through a **finalizer** — unlike the Fortran
module, where meshes are freed explicitly with `call m%free()`. Call
`close(m)` for deterministic release; it is idempotent.

Operations that produce an opaque C result (`split`, `partition`, `reorder`,
`refine`, `decimate`, `convert_cells`) always **transfer ownership** of the
mesh out of that result rather than handing back a borrow into it, so a piece
stays valid after the result is freed:

```julia
for (key, piece) in mio.split(m; by="type")
    mio.write(piece, "part_$key.vtu")    # still valid; the result is long gone
end
```

## Errors

Every failure raises a `MeshioError` carrying the C API's own thread-local
message; a status code never reaches the caller.

```julia
julia> mio.read("/nope.vtu")
ERROR: MeshioError (read error): meshio++: cannot open file '/nope.vtu'
```

Using a borrow outside its window raises `BorrowError` instead.

## Documented gaps

These are gaps in the **C ABI**, shared with the Fortran and R bindings, not
things this package chose to leave out. It invents no workaround for any of
them:

* **point / cell sets beyond regions** never reach the C++ core at all;
* the **`frozen` pin mask** of `smooth` and `decimate`;
* **per-cell-type counts** in the statistics report — use `cell_block_types`
  with `cell_block_info` instead;
* **ragged block connectivity** (polygons and polyhedra of varying size).
  `cell_block_info(m, i).is_ragged` reports them; `connectivity` then throws
  rather than returning something wrong;
* the combined **`data_manage`** — the `data_drop` / `data_keep` /
  `data_rename` primitives compose to the same effect.

One further limitation is not a C-ABI gap but a format one, worth knowing
because it is easy to trip over: **gmsh does not currently round-trip named
regions.** The writer emits the `$PhysicalNames` entry, but does not attach
the physical tag to an entity in `$Elements`, so a reader finds nothing to
reconstruct the group from. This is a pre-existing meshio++ behaviour,
reproducible from Python with no Julia involved; `abaqus` round-trips regions
correctly and is what the test suite uses.

## Tests

```sh
MESHIOPLUSPLUS_LIB=/your/prefix/lib/libmeshioplusplus.so \
  julia --project=bindings/julia/MeshioPlusPlus -e 'using Pkg; Pkg.test()'
```

# meshio++ for R

R bindings for [meshio++](https://github.com/loumalouomega/meshioplusplus):
read and write ~40 mesh formats used in scientific computing and FEM, and run
meshio++'s operations layer (surface extraction, quality, refine, decimate,
partition, smooth, slice, isosurface, …) on the results.

MIT-licensed, like the rest of meshio++. (The sibling Julia binding is *not* —
see `bindings/julia/LICENSE`.)

## Prerequisite: build and install the C library

This package binds the installed C API, exactly as the Fortran module does.
No C++ is compiled here.

```sh
./build/configure.sh --c-api --build
cmake --install build/cpp-release --prefix /your/prefix
```

Then make it discoverable, either through pkg-config:

```sh
export PKG_CONFIG_PATH=/your/prefix/lib/pkgconfig
```

or by naming the prefix directly:

```sh
export MESHIOPLUSPLUS_HOME=/your/prefix
```

At run time the shared library must also be on the loader path
(`LD_LIBRARY_PATH` on Linux, `DYLD_LIBRARY_PATH` on macOS). On Windows,
`configure` is not run: set `MESHIOPLUSPLUS_HOME` and put the DLL on `PATH`.

## Install

```sh
R CMD INSTALL bindings/r/meshioplusplus
```

## Usage

```r
library(meshioplusplus)

m <- mio_read("bracket.msh")
m
#> <mio_mesh: 9231 points, 42145 cells in 3 block(s)>

pts  <- mio_points(m)          # dim x num_points
conn <- mio_connectivity(m, 1) # nodes_per_cell x num_cells, 1-based

surf <- mio_extract_surface(m)
mio_write(surf, "bracket_surface.vtu")
mio_release(m)
```

Every exported function is prefixed `mio_`, which keeps the package clear of
base R names such as `points()`, `stats()`, `split()` and `merge()`.

## Column-major — no transpose anywhere

R matrices are column-major and the C core row-major, so

```r
mio_points(m)                        # dim x num_points
mio_connectivity(m, i)               # nodes_per_cell x num_cells
mio_point_data(m, "displacement")    # components x num_points
```

are the **same layout** as the C API's row-major `(num_points, dim)`,
`(num_cells, nodes_per_cell)` and `(num_points, components)`. Nothing is
transposed on either side — one `memcpy` and the shape is already right. This
is the same reasoning the Fortran and Julia bindings document.

The one deliberate exception is `mio_transform()`, whose 4×4 affine matrix *is*
transposed on the way out: that is a mathematical object, not a mesh array.

## Copy-only — R has no borrow

R vectors are R-managed, so the C API's zero-copy borrow **does not survive
into R** without ALTREP machinery that is out of scope here. **Every accessor
copies.** That is a real difference from the Julia and Fortran bindings, which
do hand out live views, and it is stated here rather than papered over.

Accordingly there is no `_ptr` accessor at all. The 0-based reader is named
`mio_connectivity_raw()`, deliberately *not* `mio_connectivity_ptr()`, so
nobody reads it as a borrow:

| accessor | node indices |
|---|---|
| `mio_connectivity(m, i)` | **1-based** |
| `mio_connectivity_raw(m, i)` | **0-based** — the ABI's own |

Both copy. The ±1 shift happens in that copy, which is where the Fortran and
Julia bindings put it too.

## 64-bit integers

R has no native 64-bit integer type. The `int64` arrays the C API returns —
connectivity, region entries, index maps, permutations — therefore come back as
`double`. That is exact to 2^53, far beyond any mesh this library will meet,
and it avoids a hard dependency on **bit64** for a theoretical case. The dtype
a data array was actually stored as is reported in its `"dtype"` attribute:

```r
attr(mio_point_data(m, "temperature"), "dtype")
#> [1] "float64"
```

## Indexing

Cell-block indices, data-name indices, connectivity and region entries are all
**1-based**. Index maps and permutations (`mio_refine()`,
`mio_convert_cells()`, `mio_decimate()`, `mio_partition()`, `mio_reorder()`)
are 1-based too, with the C API's `-1` "pruned / absent" sentinel becoming
**`0`** — never a valid 1-based index. That is exactly the Fortran and Julia
rule; the three bindings agree deliberately.

`mio_partition_labels()` is the exception: those are part **ids**, not indices,
so they stay in `0:(nparts - 1)`.

Region entries are 1-based with one further exception: for a `"side"` region
the second row is a **facet ordinal within the cell type**, not a mesh index,
so it is passed through unshifted.

## Memory management

A mesh handle is an external pointer with a registered finalizer, so it is
released when garbage-collected. `mio_release()` frees it immediately and is
idempotent; using a released handle is a clean R error, not a crash.

Operations that produce an opaque C result (`mio_split()`, `mio_partition()`,
`mio_reorder()`, `mio_refine()`, `mio_decimate()`, `mio_convert_cells()`)
always **transfer ownership** of the mesh out of that result rather than
returning a borrow into it, so a piece stays valid after the result is gone.

## Errors

Every failure becomes an R condition carrying the C API's own thread-local
message; a status code never reaches R.

```r
mio_read("/nope.vtu")
#> Error: meshio++: cannot open file '/nope.vtu'
```

## Documented gaps

These are gaps in the **C ABI**, shared with the Fortran and Julia bindings,
not things this package chose to leave out. It invents no workaround for any
of them:

* **point / cell sets beyond regions** never reach the C++ core at all;
* the **`frozen` pin mask** of `mio_smooth()` and `mio_decimate()`;
* **per-cell-type counts** in `mio_stats()` — use `mio_cell_block_types()` with
  `mio_cell_block_info()` instead;
* ~~ragged block connectivity~~ — **closed in v9.15.0**. `mio_polygon_block()`
  and `mio_polyhedron_block()` read them as nested 1-based lists;
  `mio_add_polygon_block()` and `mio_add_polyhedron_block()` build them.
  `mio_connectivity()` still raises, since a ragged block has no matrix;
* the combined **`data_manage`** — `mio_data_drop()` / `mio_data_keep()` /
  `mio_data_rename()` compose to the same effect.

One further limitation is not a C-ABI gap but a format one: **gmsh does not
currently round-trip named regions.** The writer emits the `$PhysicalNames`
entry but does not attach the physical tag to an entity in `$Elements`, so a
reader finds nothing to reconstruct the group from. This is a pre-existing
meshio++ behaviour, reproducible from Python; `abaqus` round-trips regions
correctly and is what the test suite uses.

A second limitation is specific to this binding: **the data setters always
write `Float64`.** `mio_add_point_data()`, `mio_append_cell_data()` and
`mio_add_field_data()` copy through `REALSXP` regardless of the R vector's
storage mode — R has no integer type reaching the C ABI here, the write-side
twin of the 64-bit-integer read limitation above. In practice this means
`mio_split(by = "region")` needs an integer cell-data tag that must come from
a *read* file (a gmsh physical group, an Exodus block id, an
`mio_isosurface()`-produced `iso:index`, …) rather than one built fresh in R —
see `example/r/03_mesh_operations.ipynb`, which hits exactly this and uses
`by = "type"` instead.

## Why plain `.Call`, not Rcpp

Zero dependencies, no C++ toolchain in `R CMD check`, and it matches the flat
C ABI the rest of meshio++'s bindings sit on. The whole surface is scalars,
vectors and matrices, where `Rf_allocVector` / `Rf_allocMatrix` plus a `memcpy`
is all that is needed; Rcpp would buy convenience this surface does not need.

## Tests

```sh
R CMD build bindings/r/meshioplusplus
R CMD check --as-cran meshioplusplus_8.3.0.tar.gz
```

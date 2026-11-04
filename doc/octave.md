# Octave / MATLAB

meshio++ ships a binding for GNU Octave and MATLAB, `bindings/octave/meshioplusplus`, layered on the [C API](/c_api) like the [Fortran](/fortran), [Julia](/julia) and [R](/r) bindings:

```matlab
m = mio.read('bracket.msh');
m
%   <mio.MeshHandle: 9231 points, 42145 cells in 3 block(s)>

surf = mio.extract_surface(m);
mio.write(surf, 'bracket_surface.vtu');
p = mio.points(surf);              % num_points x 3
t = mio.connectivity(surf, 1);     % num_cells x 3, 1-based
trisurf(t, p(:, 1), p(:, 2), p(:, 3));
```

Every function lives in the `mio` package namespace, so the binding never shadows a built-in (`points`, `merge`, `split`, `stats`, ...). The binding is a port of the R binding, function for function: R's `mio_foo(...)` is `mio.foo(...)` here, with the same arguments and defaults. Each function has help text (`help mio.read`). The [R page](/r) and the per-operation pages describe what each function computes. This page covers what differs in Octave and MATLAB.

The package is MIT-licensed, like the rest of meshio++.

## Octave and MATLAB

Development and CI use **GNU Octave** (8 or later). The sources are written for both runtimes:

- The native part is plain C over the **MEX API** (`mex.h`), which both runtimes implement. Octave's C++ `oct.h` interface is never used. One gateway, `meshioplusplus_mex`, exports every entry point.
- The `.m` files keep to the language subset the two share: `%` comments, `~=`, `end`, single-quoted strings, `fprintf`, `onCleanup`. `tests/python/test_octave_binding_parity.py` rejects Octave-only syntax (`#`, `!=`, `endfunction`, `++`, double-quoted strings, ...).
- Meshes and the other native objects are `classdef` handle classes, which both runtimes support.

MATLAB is not exercised in CI, because hosted runners do not provide it. A MATLAB-specific problem should be reported as a bug.

## Building and installing

The binding uses the **installed** C library and compiles only its own C shim:

```sh
./build/configure.sh --c-api --build
cmake --install build/cpp-release --prefix /opt/meshioplusplus
```

The library is found the same way the R package finds it:

1. `MESHIOPLUSPLUS_CFLAGS` / `MESHIOPLUSPLUS_LIBS`, for a fully manual setup.
2. `MESHIOPLUSPLUS_HOME`, an install prefix.
3. `pkg-config meshioplusplus`, through the relocatable `.pc` the C API installs (`export PKG_CONFIG_PATH=/opt/meshioplusplus/lib/pkgconfig`).

In **Octave**, install the package from a tarball:

```sh
tar -C bindings/octave -czf meshioplusplus.tar.gz meshioplusplus
octave --eval "pkg install meshioplusplus.tar.gz"
```

Then load it with `pkg load meshioplusplus`.

In **MATLAB** (or in Octave without `pkg`), add the package to the path and compile the gateway once:

```matlab
addpath('/path/to/meshioplusplus/bindings/octave/meshioplusplus/inst');
mio_build();                                    % or mio_build('prefix', '/opt/meshioplusplus')
```

`mio_build` calls `mkoctfile --mex` in Octave and `mex -R2018a` in MATLAB, and writes `meshioplusplus_mex.<mexext>` next to itself. At run time the shared library must be on the loader path: `LD_LIBRARY_PATH` on Linux, `DYLD_LIBRARY_PATH` on macOS and `PATH` on Windows. For MATLAB, set it before MATLAB starts.

If no C library is found, the build fails with a message that says how to build and install one, rather than with a bare compiler error.

## Calling convention

Arguments that R requires are positional. Optional arguments can follow positionally in R's order, or be given as name/value pairs:

```matlab
mio.connectivity(m, 2)                               % block 2
mio.connectivity(m, 'block', 2)                      % the same
mio.read('a.msh', 'gmsh')                            % format, positionally
mio.read('a.msh', 'format', 'gmsh', 'lenient', true)
mio.write(m, 'out.vtu', 'encoding', 'binary', 'codec', 'zstd')
```

The first character argument that names an option starts the name/value part. A scalar struct of options in place of all the optional arguments is accepted too. A misspelt option name is an error. The defaults are R's, translated: R's `NULL` is `[]`, `character(0)` is `{}`, and `TRUE`/`FALSE` are `true`/`false`. As in R, `'arrays', []` reads every array and `'arrays', {}` reads none.

A function that R returns `invisible(mesh)` from (the builders `mio.set_points`, `mio.add_cell_block`, ...) returns the mesh, so calls can chain. One R function changes shape: `mio.with_provenance` takes a function handle, `mio.with_provenance(@() ..., ...)`, where R takes an expression.

## Array layout: one entity per row

MATLAB puts one point, cell or entry per **row**. That is the shape `patch`, `trisurf` and `scatter3` expect, so every per-entity array crosses that way:

| Value | Shape | Class |
| --- | --- | --- |
| `mio.points(m)` | `num_points x dim` | `double` |
| `mio.connectivity(m, b)` | `num_cells x nodes_per_cell`, 1-based | `double` |
| `mio.point_data(m, name)` | `num_points x components...` (a scalar field is a column) | stored class |
| `mio.cell_data(m, name, b)` | `num_cells x components...` | stored class |
| region `entries` | `num_entries x stride` | `double` |
| index maps, id lists | `n x 1` columns | `double` |
| `bbox` | `2 x 3`: the minimum corner row, then the maximum corner row | `double` |

Here R differs: R keeps the C API's row-major memory as a column-major `(dim x n)` matrix and never transposes. The Octave binding copies and transposes in the same pass, so it costs nothing extra over R's copy. Setters take the same row orientation. Mathematical matrices keep their meaning: a 4x4 affine for `mio.transform` or `mio.match_periodic_nodes` is written as it reads, `T(1:3, 4)` being the translation.

Ragged blocks come back as cell arrays. `mio.polygon_block` returns one 1-based **row** vector per cell. `mio.polyhedron_block` returns, per cell, a cell array of face rows. The matching builders `mio.add_polygon_block` and `mio.add_polyhedron_block` take the same shapes.

## Types

MATLAB has every integer width, so data arrays keep their **stored class**: an `int32` material id comes back as `int32`, a `float32` field as `single`. A setter keeps the class it is given, and `logical` is stored as `uint8`. R, which lacks most of these types, converts everything to `double`. Points, connectivity, indices and counts are `double` here, exact to 2^53. Integer results that R returns as `integer(1)` come back as `double` scalars.

## Indexing

Connectivity, cell blocks, region indices, point and cell ids, and index maps are **1-based**. In an index map, the C API's "absent" sentinel `-1` becomes `0`, which is never a valid 1-based index. This is the rule the Fortran, Julia and R bindings follow. A side region's second column is a facet ordinal within the cell type, not a mesh index, so it is not shifted. `mio.connectivity_raw` returns the C API's 0-based numbering, as a copy.

## Handles and memory

A mesh is a `mio.MeshHandle`. Format side channels (`mio.FormatInfo`), sequences (`mio.SequenceHandle`), and XDMF and Exodus series writers (`mio.XdmfSeries`, `mio.ExodusSeries`) are handles too. All derive from `mio.Handle < handle`:

- The native object is freed when the last MATLAB reference goes away (`delete`), or earlier by `mio.release(m)` or `m.release()`. Releasing is idempotent.
- A handle holds only an opaque id into a registry inside the MEX. Each id carries a generation counter, so a released handle, a handle loaded back from a `.mat` file or a forged id is **rejected with an error, never dereferenced**. Passing the wrong kind (a series where a mesh is expected) is an error too. This is the counterpart of R's external-pointer tags.
- The MEX stays locked (`mexLock`) while any handle is live, so `clear mex` cannot unload it from under them.
- Handles have reference semantics: `b = a` shares one mesh. Operations always return new meshes and leave their inputs untouched.

Every accessor copies, as in R: the C API's zero-copy borrows do not survive into runtime-managed arrays. To work on a whole mesh as plain data, copy it into a struct and back:

```matlab
s = m.toStruct();                       % points, cells, point_data, cell_data, field_data, regions
m2 = mio.MeshHandle.fromStruct(s);
```

Data names such as `gmsh:physical` are not valid field names, so the struct's data sections are struct arrays of `name`/`data` pairs.

## Error handling

A failing call raises a MATLAB error whose identifier starts with `meshioplusplus:` and whose message is the C library's own (`mio.last_error()` returns the same text):

```matlab
try
    mio.read('missing.vtu');
catch err
    disp(err.identifier);   % meshioplusplus:read
    disp(err.message);
end
```

A `mio_status` code never reaches user code. Temporary memory is MEX-managed and is released if a call raises, and native objects a call has allocated are freed before it raises.

## Named regions, selective reads and series

These follow the R binding exactly, apart from the layout rules above:

- `mio.regions(m)` returns an `n x 1` struct array with the fields `name`, `kind`, `dim`, `tag` and `entries`.
- `mio.add_region`, `mio.edit_regions`, `mio.remove_region` and `mio.match_periodic_nodes` edit and pair regions. See [regions](./regions.md).
- `mio.read(path, 'arrays', {...}, 'points_only', true, 'time_step', k, 'piece', p)` and `mio.read_metadata(path)` do [selective reads](./selective_read.md).
- `mio.xdmf_series`, `mio.exodus_series` and `mio.sequence` return handle objects whose `_write_*`, `_flush`, `_finalize`, `_step` and similar functions mirror R's. See [XDMF time series](./xdmf_time_series.md), [Exodus](./formats/exodus.md#stateful-series-writing) and [sequences](./sequences.md).
- Pipeline reports (`mio.pipeline_run_json_report`, ...) return the JSON text as a char row vector, which `jsondecode` parses.

## Software rendering

`f = mio.render(mesh, ...)` rasterises a mesh with the core's software renderer. Every option is a name/value pair named after its `mio_render_opts` field, exactly as in R. The frame is a struct:

- `f.rgba` is an `H x W x 4` `uint8` array with straight alpha. Row 1 is the top of the picture and column 1 its left edge, so `image(f.rgba(:, :, 1:3))` shows it upright and `imwrite(f.rgba(:, :, 1:3), 'view.png', 'Alpha', f.rgba(:, :, 4))` saves it.
- `f.cell_ids` is `H x W` in the same orientation, holding **1-based** cell ids with `0` where no cell is drawn. R returns 0-based ids, with `-1` for none.
- `f.range` is `[vmin vmax]`, or `[]` when no field is mapped.
- `f.notes` is an `n x 1` cellstr.

`mio.frame_array(f)` is `double(f.rgba) / 255`. `mio.render_png` returns the encoded PNG as a `uint8` column, `mio.write_snapshot` writes a file, and `mio.render_text` returns terminal art as a char row (UTF-8 decoded). See [terminal rendering](/tui) for every option.

## Where the shapes differ from R

Lists follow one rule. A named list becomes a 1x1 struct. A list of same-shaped records (regions, partitions, data_info entries) becomes an `n x 1` struct array. An unnamed list becomes an `n x 1` cell, and a character vector an `n x 1` cellstr. An R data frame becomes a struct of `n x 1` columns. The few other differences:

- `[pieces, keys] = mio.split(m, ...)` returns an `n x 1` cell of meshes and an `n x 1` cellstr of keys, where R returns a named list.
- A single coordinate (an origin, a spacing, a bbox corner, a centroid, a worst point) is a `1 x 3` row, like one row of `mio.points`. Statistics triples (min, mean, max) are `3 x 1` columns.
- Warning texts that R returns as one string with newlines come back as an `n x 1` cellstr.
- `mio.line` names R's `end` option `'stop'`, because `end` is reserved.
- Index arguments (`frozen` nodes, `sources`, refined `cells`) must be positive integers. A `0` or a fraction raises "is 1-based here", where R would pass it on to the library.

Octave prints `warning: function .../+mio/diff.m shadows a built-in function` (and the same for `merge` and `pipe`) when the package is added to the path. Package functions cannot shadow anything, so the built-in `diff`, `merge` and `pipe` keep working. `warning('off', 'Octave:shadowed-function')` silences the message.

## Tests

`bindings/octave/meshioplusplus/tests` holds the suite. The files port the R testthat files, and the fixtures are the same five-point, two-tetra mesh, asymmetric in every axis, so a transposed mapping or a missed 1-based shift cannot cancel out:

```sh
cd bindings/octave/meshioplusplus/tests
octave --no-gui --eval "exit(run_tests() > 0)"          # or run_tests('render')
```

`run_tests` works in MATLAB too (`matlab -batch "exit(run_tests() > 0)"`). Tests that need an optional dependency (MED and HDF5, netCDF, OpenCASCADE, ...) or an environment-supplied fixture are skipped. The `octave` CI job builds the binding both ways (`pkg install` and `mio_build`) and runs the suite after each. The static Python test `tests/python/test_octave_binding_parity.py` fails when:

- an R export lacks its `inst/+mio/<name>.m` wrapper,
- an R `.Call` entry point lacks its MEX command, or the two have different argument counts,
- a wrapper calls an unregistered command, or
- a `.m` file uses Octave-only syntax.

## Extending the binding

A new R entry point gets its Octave counterpart in the same change:

- Add a handler `static mxArray *O_mio_x(const mxArray **a)` to the `src/` file that mirrors the R source file, with its arguments in R's `.Call` order, and add the handler to that file's table.
- Add an `inst/+mio/x.m` wrapper that ports the R function's defaults and normalisation through `mio.internal.opts` / `choice` / `block`.
- Add a test in `tests/test_*.m`.

`src/mio_oct.h` documents the helper layer and the rules above.

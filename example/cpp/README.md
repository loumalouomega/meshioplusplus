# Examples (C++)

The C++ counterpart of [`../python/`](../python/): the same tour of
meshio++, called directly against the **C++ core** instead of the Python
bindings, on the same bundled bracket geometry. No Python is involved --
these notebooks run on the [xeus-cpp](https://github.com/compiler-research/xeus-cpp)
Jupyter kernel (a `clang-repl`-based C++ interpreter) against
[`meshioplusplus.hpp`](../../src/single_include/meshioplusplus/meshioplusplus.hpp),
the committed single-header amalgamation -- `#include` it,
`#define MESHIOPLUSPLUS_IMPLEMENTATION` once, and every format and operation
in the core is available with no separate build step.

| Notebook | What it shows |
|----------|---------------|
| [`01_read_and_visualize.ipynb`](01_read_and_visualize.ipynb) | Read the mesh, inspect it, render the full part and a cropped interior view, and colour a corner by quality/height. |
| [`02_convert_and_inspect.ipynb`](02_convert_and_inspect.ipynb) | Convert to VTU/VTK/XDMF/Gmsh/PLY, compare file sizes, verify every round trip. |
| [`03_mesh_operations.ipynb`](03_mesh_operations.ipynb) | The same operations tour as [`../python/03_mesh_operations.ipynb`](../python/03_mesh_operations.ipynb): surface/skin extraction, quality, reorder, diff, sniff, transform, clean, crop, merge, split, stats, convert_cells, refine, partition, smooth, interpolate, the five data operations, and selective reads. |

## Rendering without PyVista

There is no PyVista/VTK in C++, so instead of rasterizing a GL view, these
notebooks lean on meshio++'s own **SVG writer**: it already does the camera
projection and data-driven colouring PyVista would
([`doc/formats/svg.md#data-driven-colouring`](../../doc/formats/svg.md#data-driven-colouring)),
needs no extra dependency, and keeps "meshio++'s core has no rendering
dependency" true. [`mio_notebook.hpp`](mio_notebook.hpp) is the small,
notebook-only helper that makes this ergonomic:

- `render(mesh, color_by, cmap, colorbar, vmin, vmax, azimuth, elevation, width)` wraps `write_svg`, reads the result back, and returns a `SvgImage`.
- `bar_chart(...)` / `histogram_chart(...)` are hand-rolled SVG primitives (no plotting library) for the handful of numeric charts the Python notebooks draw with matplotlib, coloured from meshio++'s own `detail::colormap_lookup` tables so the palette matches the mesh renders.
- `hex_block(n, spacing)` builds the small synthetic hexahedron-grid fixtures several operation demos use (`convert_cells`, `refine`, `smooth`, `interpolate`).
- An `SvgImage`/`mime_bundle_repr` pair plugs into xeus-cpp's rich-display mechanism; cells call `xcpp::display(render(...))` explicitly (see "xeus-cpp quirks" below for why).

Only meshio++'s **3 built-in colormaps** are available in C++: `viridis`,
`coolwarm`, `turbo` (`detail::colormap_names()`) -- no `cividis`/`plasma`/`magma`.

Renders **crop to a legible region** before drawing wherever the source mesh
has more than a few thousand cells. An SVG's size scales linearly with facet
count (unlike a raster PNG), so naively rendering the bracket's full ~58k-facet
skin makes a multi-MB file whose individual elements are sub-pixel -- the
exact issue [`tools/gen_doc_images.py`](../../tools/gen_doc_images.py) already
works around for the project's static doc figures. `03_mesh_operations.ipynb`
computes a reusable `vizCorner` crop early and reuses (or re-derives, e.g.
centred on a deformation rather than a domain corner) it throughout.

## `example.vtu`, not `example/example.msh`

The input here is [`example.vtu`](example.vtu), not the original
`example/example.msh`. That was originally because Gmsh 4.1's `$Entities`
section (present in that file) was a gap in meshio++'s **C++** reader, which
the Python bindings papered over with a fallback a C++-only kernel does not
have. **Since v9.7.0 the C++ reader handles `$Entities`**, so the `.msh` could
now be read directly; the `.vtu` is kept because these notebooks are committed
with outputs and re-deriving them is a separate change. It was pre-converted
once, from the repo root, with the Python bindings:

```sh
.venv/bin/python3 -c "
import numpy as np
import meshioplusplus as mp
src = mp.read('example/example.msh')
tri = np.concatenate([cb.data for cb in src.cells if cb.type == 'triangle'])
tet = np.concatenate([cb.data for cb in src.cells if cb.type == 'tetra'])
mesh = mp.Mesh(src.points, [('triangle', tri), ('tetra', tet)])
mp.write('example/cpp/example.vtu', mesh, binary=False)
"
```

ASCII, not binary: the amalgamation's default build has no
`MESHIOPLUSPLUS_HAS_ZLIB`, and VTU's binary writer defaults to zlib
compression (see "xeus-cpp / amalgamation quirks" below).

## Setting up the kernel

The xeus-cpp kernel is a heavy, compiler-adjacent dependency (it bundles a
`clang-repl`) that only ships via conda-forge -- there is no pip wheel. Install
it into an isolated [micromamba](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html)
environment (kept out of the repo and out of `.venv`; `.micromamba/` is
gitignored):

```sh
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj -C /tmp/mmbin bin/micromamba
export MAMBA_ROOT_PREFIX=$PWD/.micromamba
/tmp/mmbin/bin/micromamba create -y -n xcpp -c conda-forge xeus-cpp jupyterlab notebook nbconvert
```

This registers several kernels (`xcpp17`/`xcpp20`/`xcpp23`, `xc11`/`xc17`/`xc23`,
...); the notebooks here pin `xcpp20` to match the core's C++20 baseline
(`CMakeLists.txt`'s `CMAKE_CXX_STANDARD`).

## Running them

The notebooks are committed **with their outputs** so they render on GitHub
without any setup, same convention as `example/`. To re-run:

```sh
export MAMBA_ROOT_PREFIX=$PWD/.micromamba
/tmp/mmbin/bin/micromamba run -n xcpp jupyter nbconvert --to notebook --execute --inplace example/cpp/*.ipynb
```

nbconvert sets the kernel's working directory to the notebook's own folder,
which is why the cells use paths relative to `example/cpp/`
(`../../src/single_include/...`, `mio_notebook.hpp`, `example.vtu`) rather than
absolute ones.

## xeus-cpp / amalgamation quirks worth knowing before editing these notebooks

- **Cells share one persistent global namespace**, like a Python kernel's --
  a `for` loop or `NDArray`/`Mesh`/etc. declared at a cell's top level (not
  inside `{ }`) is visible, and *must not collide*, in every later cell. This
  bit twice while writing `03_mesh_operations.ipynb`: a local named `dup`
  collided with libc's `dup(int)` (pulled in transitively), and `cs`/`fs`
  collided with `namespace fs = std::filesystem;`, which the amalgamation
  itself declares. When in doubt, scope scratch variables in `{ }` or give
  them a distinctive name.
- **No auto-print of a cell's trailing expression** (unlike classic
  xeus-cling notebooks) -- this build of xeus-cpp (0.10.0) does not evaluate
  and display a bare trailing value. Use `std::cout` for text and
  `xcpp::display(...)` (from `xcpp/xdisplay.hpp`, included by
  `mio_notebook.hpp`) for rich (e.g. SVG) output.
- **A capture-default lambda (`[&]`/`[=]`) cannot appear directly in a
  cell's top-level code** -- clang-repl rejects it as "non-local lambda
  expression cannot have a capture-default", a restriction real translation-unit
  compilation doesn't have. Wrap such logic in a captureless lambda, an
  ordinary loop, or a `{ }` block that is itself inside a function.
- **Confirm a fast way to catch both of the above before spending a slow
  nbconvert cycle on them**: concatenate every code cell's source into one
  `int main() { ... }` (top-level `#include`/`#define`/`using` lines stay
  outside `main`) and compile with a real `g++`/`clang++` -- it reproduces
  the REPL's shared-global-scope semantics closely enough to catch
  redefinitions in seconds, though it can't check anything that depends on
  the real `libxeus`/`libclang-repl` runtime (e.g. `xcpp::display` itself
  won't link standalone).
- Optional deps (`zlib`, `HDF5`, ...) are **compiled out** by default in the
  amalgamation (no `-DMESHIOPLUSPLUS_HAS_*`), unlike the Python wheel. The
  registry's default writer table bakes in zlib-compressed VTU/VTP
  (`registry.cpp`), so call `write_vtu`/`write_vtp` directly with an explicit
  `zlib=false` rather than going through `registry_writers()` for those two
  formats.

# Examples (Julia)

The Julia counterpart of [`../python/`](../python/) and [`../cpp/`](../cpp/): the same tour of meshio++,
called through the **`MeshioPlusPlus` Julia package** (`bindings/julia/`) instead of Python or the C++
core directly. This binding rides on the same flat C API as [Fortran](../../doc/fortran.md) and
[R](../r/) — see [`doc/julia.md`](../../doc/julia.md) for the package itself.

| Notebook | What it shows |
|----------|---------------|
| [`01_read_and_visualize.ipynb`](01_read_and_visualize.ipynb) | Read the mesh, inspect it, render the full part and a cropped interior view, and chart element quality/height for a corner. |
| [`02_convert_and_inspect.ipynb`](02_convert_and_inspect.ipynb) | Convert to VTU/VTK/XDMF/Gmsh/PLY, compare file sizes, verify every round trip. |
| [`03_mesh_operations.ipynb`](03_mesh_operations.ipynb) | The same operations tour as [`../cpp/03_mesh_operations.ipynb`](../cpp/03_mesh_operations.ipynb): surface/skin extraction, quality, reorder, diff, sniff, transform, clean, crop, merge, split, stats, convert_cells, refine, decimate, partition, smooth, interpolate, slice, isosurface, the five data operations, and selective reads. |

## Rendering without PyVista — and without per-call SVG colouring

There is no PyVista/VTK in Julia, so like the C++ notebooks these lean on meshio++'s own **SVG writer**.
Unlike the C++ notebooks, though, this binding reaches that writer only through the **flat C API**'s
`mio.write(mesh, path; format="svg")` — and `CLAUDE.md` documents that data-driven colouring is a
**flat-ABI gap**: `registry.cpp`'s `(path, mesh)` writer lambdas cannot carry per-call parameters (`color_by`,
`cmap`, `colorbar`, camera angles), so every render here uses the writer's fixed default styling. Where the
C++ notebook colours a mesh by a quality metric or a field, these notebooks show the same information as a
small chart instead — [`mio_notebook.jl`](mio_notebook.jl)'s `bar_chart`/`histogram_chart`, hand-rolled SVG
primitives with no plotting-library dependency (mirroring `mio_notebook.hpp`'s approach, translated to
Julia), since meshio++'s own colormap tables aren't reachable from the C API either.

- `render(mesh)` wraps `mio.write(..., format="svg")`, reads the result back, and returns an `SvgImage`
  that Jupyter displays via its `image/svg+xml` MIME hook.
- `bar_chart(...)` / `histogram_chart(...)` are the small numeric charts standing in for colour.
- `hex_block(n; spacing=1.0)` builds the small synthetic hexahedron-grid fixtures several operation demos
  use (`convert_cells`, `refine`, `smooth`, `interpolate`, `slice`, `isosurface`), mirroring the C++ and
  Python notebooks' equivalents.

## `example.vtu`, not `example/example.msh`

The input here is [`example.vtu`](example.vtu), copied verbatim from [`../cpp/example.vtu`](../cpp/example.vtu)
(both are Git-LFS-tracked, so this costs no meaningful repo size) rather than reading `example/example.msh`
directly: this binding sits on the same **C++ reader** as the C++ notebooks, and Gmsh 4.1's `$Entities`
section (present in the original file) is a documented gap there — Python's fallback isn't available with
no Python involved. See [`../cpp/README.md`](../cpp/README.md) for how it was produced.

## Setting up the kernel

[IJulia](https://github.com/JuliaLang/IJulia.jl) is the standard, officially-maintained Jupyter kernel for
Julia:

```sh
julia -e 'using Pkg; Pkg.add("IJulia")'
```

This registers a `julia-1.x` kernelspec (`~/.local/share/jupyter/kernels/`, or wherever `jupyter --paths`
points) usable from any Jupyter frontend, including one from an entirely different Python environment (no
`--project` is baked into the kernel; each notebook's first cell `Pkg.activate`s
[`bindings/julia/MeshioPlusPlus`](../../bindings/julia/MeshioPlusPlus) itself, so the kernel's own default
project doesn't matter).

## Building and installing the C API first

This binding consumes the **installed** C library — no separate build step for the notebooks themselves,
but the library has to exist first:

```sh
cmake -S . -B build-capi -DCMAKE_BUILD_TYPE=Release -DMESHIOPLUSPLUS_BUILD_PYTHON=OFF \
  -DMESHIOPLUSPLUS_BUILD_C_API=ON -DMESHIOPLUSPLUS_WITH_HDF5=OFF -DMESHIOPLUSPLUS_WITH_NETCDF=OFF
cmake --build build-capi
cmake --install build-capi --prefix inst-nohdf5
```

::: warning Build with `-DMESHIOPLUSPLUS_WITH_HDF5=OFF` for Julia specifically
On Debian/Ubuntu, a C API built **with** HDF5 links `libhdf5_openmpi`, which links the *system* `libcurl`.
Julia bundles its own `libcurl`, and an IJulia kernel subprocess's `dlopen` of that library can fail with
`libcurl.so.4: version 'CURL_OPENSSL_4' not found` — a real Debian/Ubuntu + Julia interaction (confirmed by
building the same C API with HDF5 on vs. off), not a meshio++ bug. See [`doc/julia.md`](../../doc/julia.md).
None of these notebooks need HDF5-backed formats.
:::

## Running them

The notebooks are committed **with their outputs** so they render on GitHub without any setup. To re-run:

```sh
jupyter nbconvert --to notebook --execute --inplace example/julia/*.ipynb
```

nbconvert sets the kernel's working directory to the notebook's own folder, which is why the cells use
paths relative to `example/julia/` (`example.vtu`, `mio_notebook.jl`, `../../bindings/julia/MeshioPlusPlus`)
rather than absolute ones, and why the first cell of every notebook sets
`ENV["MESHIOPLUSPLUS_LIB"]` relative to `@__DIR__` before `Pkg.activate`ing and `using`-ing the package.

## A correctness note from writing these notebooks

The smoothing demo in `03_mesh_operations.ipynb` caught a real bug: this binding's `smooth()` originally
defaulted `mu` to `-0.53` (both here and in the [R binding](../r/)), an invented value never checked
against the actual default (`-0.34`, per Fortran's `meshioplusplus.f90` and the Python bindings). A `-0.53`
pass-band is wide enough to reliably make a tangled mesh **worse**, not better — caught by cross-checking
the same synthetic jittered-hex demo against the Python bindings, which recovered cleanly while Julia and R
did not. Both bindings' default is now `-0.34`, matching Fortran/Python, and the demo recovers as expected.

# Examples (R)

The R counterpart of [`../python/`](../python/) and [`../cpp/`](../cpp/): the same tour of meshio++,
called through the **`meshioplusplus` R package** (`bindings/r/`) instead of Python or the C++ core
directly. This binding rides on the same flat C API as [Fortran](../../doc/fortran.md) and
[Julia](../julia/) — see [`doc/r.md`](../../doc/r.md) for the package itself.

| Notebook | What it shows |
|----------|---------------|
| [`01_read_and_visualize.ipynb`](01_read_and_visualize.ipynb) | Read the mesh, inspect it, render the full part and a cropped interior view, and chart element quality/height for a corner. |
| [`02_convert_and_inspect.ipynb`](02_convert_and_inspect.ipynb) | Convert to VTU/VTK/XDMF/Gmsh/PLY, compare file sizes, verify every round trip. |
| [`03_mesh_operations.ipynb`](03_mesh_operations.ipynb) | The same operations tour as [`../cpp/03_mesh_operations.ipynb`](../cpp/03_mesh_operations.ipynb): surface/skin extraction, quality, reorder, diff, sniff, transform, clean, crop, merge, split, stats, convert_cells, refine, decimate, partition, smooth, interpolate, slice, isosurface, the five data operations, and selective reads. |

## Rendering without PyVista — and without per-call SVG colouring

There is no PyVista/VTK in R, so like the C++ notebooks these lean on meshio++'s own **SVG writer**. Unlike
the C++ notebooks, this binding reaches that writer only through the **flat C API**'s
`mio_write(mesh, path, format = "svg")` — and `CLAUDE.md` documents that data-driven colouring is a
**flat-ABI gap**: `registry.cpp`'s `(path, mesh)` writer lambdas cannot carry per-call parameters
(`color_by`, `cmap`, `colorbar`, camera angles), so every render here uses the writer's fixed default
styling. Where the C++ notebook colours a mesh by a quality metric or a field, these notebooks show the
same information as a small chart instead — [`mio_notebook.R`](mio_notebook.R)'s
`mio_bar_chart`/`mio_histogram`, plain **base R graphics** (`barplot`/`hist`) that IRkernel automatically
captures as an inline PNG, needing no extra package at all.

- `render(mesh)` wraps `mio_write(..., format = "svg")`, reads the result back, and displays it via
  `IRdisplay::display_svg()`.
- `mio_bar_chart(...)` / `mio_histogram(...)` are the small numeric charts standing in for colour.
- `hex_block(n, spacing = 1.0)` builds the small synthetic hexahedron-grid fixtures several operation demos
  use (`convert_cells`, `refine`, `smooth`, `interpolate`, `slice`, `isosurface`), mirroring the C++, Python
  and Julia notebooks' equivalents.

## `example.vtu`, not `example/example.msh`

The input here is [`example.vtu`](example.vtu), copied verbatim from [`../cpp/example.vtu`](../cpp/example.vtu)
(both are Git-LFS-tracked, so this costs no meaningful repo size) rather than reading `example/example.msh`
directly: this binding sits on the same **C++ reader** as the C++ notebooks, and Gmsh 4.1's `$Entities`
section (present in the original file) is a documented gap there. See
[`../cpp/README.md`](../cpp/README.md) for how it was produced.

## Setting up the kernel

[IRkernel](https://irkernel.github.io/) is the standard, officially-maintained Jupyter kernel for R. It and
the meshio++ R package are installed into an isolated
[micromamba](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html) environment, the same
pattern [`../cpp/README.md`](../cpp/README.md) uses for its xeus-cpp kernel (kept out of the repo and out
of `.venv`; `.micromamba/` is gitignored):

```sh
export MAMBA_ROOT_PREFIX=$PWD/.micromamba
micromamba create -y -n r-notebooks -c conda-forge r-base r-irkernel c-compiler
```

`c-compiler` is required even though the meshio++ R package's own build doesn't need one beyond what
`R CMD INSTALL` already invokes — conda-forge's `r-base` does not pull in a compiler by default, and
without one `R CMD INSTALL` fails looking for `x86_64-conda-linux-gnu-cc`.

Then build and install the R package into that environment (see the C API step below first):

```sh
PKG_CONFIG_PATH=$PWD/inst-nohdf5/lib/pkgconfig \
  .micromamba/envs/r-notebooks/lib/R/bin/R CMD INSTALL bindings/r/meshioplusplus
.micromamba/envs/r-notebooks/lib/R/bin/Rscript -e \
  'IRkernel::installspec(name = "ir-meshioplusplus", displayname = "R (meshioplusplus)")'
```

`IRkernel::installspec` writes a `kernel.json` with no `env` block by default; add
`"env": {"LD_LIBRARY_PATH": "<repo>/inst-nohdf5/lib"}` to the generated
`~/.local/share/jupyter/kernels/ir-meshioplusplus/kernel.json` so the kernel finds
`libmeshioplusplus.so` at runtime regardless of the invoking shell's environment.

## Building and installing the C API first

```sh
cmake -S . -B build-capi -DCMAKE_BUILD_TYPE=Release -DMESHIOPLUSPLUS_BUILD_PYTHON=OFF \
  -DMESHIOPLUSPLUS_BUILD_C_API=ON -DMESHIOPLUSPLUS_WITH_HDF5=OFF -DMESHIOPLUSPLUS_WITH_NETCDF=OFF
cmake --build build-capi
cmake --install build-capi --prefix inst-nohdf5
```

HDF5 is off here to match the sibling [Julia notebooks](../julia/), which genuinely need it off (see that
folder's README); R itself has no such conflict, but reusing one prefix for both keeps the file-size
comparisons in `02_convert_and_inspect.ipynb` comparable (XDMF's data format follows the build — HDF5 when
available, inline XML otherwise — so a mixed setup would make an unrelated file-size difference look like a
per-language one).

## Running them

The notebooks are committed **with their outputs** so they render on GitHub without any setup. To re-run:

```sh
jupyter nbconvert --to notebook --execute --inplace example/r/*.ipynb
```

nbconvert sets the kernel's working directory to the notebook's own folder, which is why the cells use
paths relative to `example/r/` (`example.vtu`, `mio_notebook.R`) rather than absolute ones.

## Two correctness notes from writing these notebooks

- The smoothing demo in `03_mesh_operations.ipynb` caught a real bug: this binding's `mio_smooth()`
  originally defaulted `mu` to `-0.53` (both here and in the [Julia binding](../julia/)), an invented value
  never checked against the actual default (`-0.34`, per Fortran's `meshioplusplus.f90` and the Python
  bindings). A `-0.53` pass-band is wide enough to reliably make a tangled mesh **worse**, not better —
  caught by cross-checking the same synthetic jittered-hex demo against the Python bindings, which
  recovered cleanly while Julia and R did not. Both bindings' default is now `-0.34`, matching
  Fortran/Python, and the demo recovers as expected.
- The split demo surfaced a binding-specific gap, not a bug: `mio_split(by = "region")` needs a genuinely
  **integer** cell-data tag, but every one of this package's data setters
  (`mio_add_point_data`/`mio_append_cell_data`/`mio_add_field_data`) writes `Float64` regardless of the R
  vector's storage mode — R has no integer type reaching the C ABI here, the write-side twin of the
  64-bit-integer read limitation `doc/r.md` already documents. A region tag built fresh in R therefore
  cannot drive `by = "region"`; the notebook uses `by = "type"` instead, and the isosurface demo's
  `iso:index` split works because that tag is produced by the C++ core itself, never passed through an R
  setter. See [`doc/r.md`](../../doc/r.md#documented-gaps).

# Examples (Octave/MATLAB)

The Octave/MATLAB counterpart of [`../r/`](../r/), [`../python/`](../python/) and [`../cpp/`](../cpp/): the same tour of meshio++, called through the **Octave/MATLAB binding** (`bindings/octave/meshioplusplus`) instead of Python or the C++ core directly. The binding is a port of the R package, function for function: R's `mio_foo(...)` is `mio.foo(...)` here. It rides on the same flat C API as [Fortran](../../doc/fortran.md), [Julia](../julia/) and [R](../r/); see [`doc/octave.md`](../../doc/octave.md) for the binding itself.

| Notebook | What it shows |
|----------|---------------|
| [`01_read_and_visualize.ipynb`](01_read_and_visualize.ipynb) | Read the mesh, inspect it, render the full part and a cropped interior view, colour a corner by element quality and by height, and chart both. |
| [`02_convert_and_inspect.ipynb`](02_convert_and_inspect.ipynb) | Convert to VTU/VTK/XDMF/Gmsh/PLY, compare file sizes, verify every round trip. |
| [`03_mesh_operations.ipynb`](03_mesh_operations.ipynb) | The same operations tour as [`../r/03_mesh_operations.ipynb`](../r/03_mesh_operations.ipynb): surface/skin extraction, quality, reorder, diff, sniff, transform, clean, crop, merge, split (including an `int32` region tag, which R cannot build), stats, convert_cells, refine, decimate, partition, smooth, interpolate, slice (also as a `trisurf`), isosurface, gradient, the five data operations, sequences and selective reads. |
| [`04_terminal_rendering.ipynb`](04_terminal_rendering.ipynb) | The software rasterizer with every option reached through the C API (`mio.render`, `mio.render_png`, `mio.render_text`, `mio.write_snapshot`): colouring by data, the `cell_ids` pick buffer, contour lines, arrows, a warp, streamlines on the skin and through a cut-away, the terminal form and `write_snapshot`. |

## Rendering

There is no PyVista/VTK in Octave or MATLAB. The R notebooks draw through meshio++'s SVG writer, which the flat C API reaches only with its fixed default styling (`doc/formats/svg.md` documents that as a flat-ABI gap). These notebooks draw through the **software rasterizer** instead: `mio.render` returns a frame whose `rgba` field is an `H x W x 4` `uint8` array with row 1 at the top, so `image(f.rgba(:, :, 1:3)); axis image off` shows it upright, and every option (`color_by`, `cmap`, `colorbar`, `cutaway`, ...) is available. The octave kernel captures each figure as an inline PNG. Charts use plain `barh`, `hist` and `trisurf`.

The helpers are notebook-only, one function per file so that MATLAB finds them too:

- [`mio_notebook.m`](mio_notebook.m) is the setup script every notebook runs first. It compiles the MEX gateway with `mio_build` when it is missing, silences Octave's gnuplot notice and sets a white figure background.
- [`render.m`](render.m) `render(mesh, ...)` draws a mesh with smooth shading and feature edges on white, through [`mio_frame_image.m`](mio_frame_image.m), and prints the frame's colour notes when a field is mapped. Extra name/value pairs go to `mio.render`.
- [`mio_frame_image.m`](mio_frame_image.m) `mio_frame_image(mesh, ...)` renders a frame and shows it in a figure the size of the frame.
- [`mio_bar_chart.m`](mio_bar_chart.m) and [`mio_histogram.m`](mio_histogram.m) are the small numeric charts.
- [`hex_block.m`](hex_block.m) `hex_block(n, spacing)` builds the synthetic hexahedron-grid fixture several operation demos use, mirroring the C++, Python, Julia and R notebooks.
- [`mio_nb_temp_path.m`](mio_nb_temp_path.m) names a fresh scratch file.

Each code cell draws at most one figure. The octave kernel exports all of a cell's figures at the size of the first one and in no fixed order, so one figure per cell keeps the outputs stable.

## `example.vtu`

The input is [`example.vtu`](example.vtu), the same file as [`../r/example.vtu`](../r/example.vtu) and [`../cpp/example.vtu`](../cpp/example.vtu), copied verbatim. It is Git-LFS-tracked like those copies; run `git lfs pull` after cloning. See [`../cpp/README.md`](../cpp/README.md) for how it was produced.

## Setting up

### Octave, the kernel and inline graphics

```sh
sudo apt-get install -y octave octave-dev gnuplot-nox fonts-freefont-otf
pip install octave_kernel
```

- [`octave_kernel`](https://github.com/Calysto/octave_kernel) registers a Jupyter kernel named `octave`. It runs `octave` from `PATH` (or `OCTAVE_EXECUTABLE`).
- `octave-dev` provides `mkoctfile`, which `mio_build` uses to compile the MEX gateway.
- Without a display, Octave's only graphics toolkit is **gnuplot**, so inline plots need the `gnuplot` binary: without it `available_graphics_toolkits` is empty and nothing can be drawn. `gnuplot-nox` is enough. The kernel prints each figure with `print`, so no X server is needed.
- `fonts-freefont-otf` provides `FreeSans.otf`, the font Octave's text renderer is built against. A `--no-install-recommends` install of Octave leaves it out, and then every axes fails with `ft_text_renderer: invalid bounding box, cannot render`. `fonts-freefont-ttf` is not a substitute.

With an X server (or `xvfb-run`), the `qt` and `fltk` toolkits work as well; the notebooks do not depend on which one draws.

### The C library

```sh
./build/configure.sh --c-api --build
cmake --install build/cpp-release --prefix "$PWD/inst-capi"
```

The committed outputs were produced against a C library built with HDF5 on, so the XDMF file in `02_convert_and_inspect.ipynb` keeps its heavy data in a `.h5` file next to it (the R and Julia notebooks were run with HDF5 off and show inline XML). The size comparison counts both files, so it holds either way.

### The binding

The notebooks add the binding to the path themselves, relative to this folder:

```matlab
warning('off', 'Octave:shadowed-function');   % a known false positive, see doc/octave.md
addpath(fullfile(pwd, '..', '..', 'bindings', 'octave', 'meshioplusplus', 'inst'));
mio_notebook
```

`mio_notebook` builds `inst/meshioplusplus_mex.mex` with `mio_build` on first use (`inst/*.mex` is gitignored). `mio_build` finds the C library through `MESHIOPLUSPLUS_HOME`, `MESHIOPLUSPLUS_CFLAGS`/`MESHIOPLUSPLUS_LIBS` or `pkg-config`, so set one of them for the first run, or build once by hand:

```sh
MESHIOPLUSPLUS_HOME=$PWD/inst-capi octave --no-gui --eval \
  "addpath('bindings/octave/meshioplusplus/inst'); mio_build()"
```

At run time `libmeshioplusplus.so` must be on the loader path. Export `LD_LIBRARY_PATH` before starting Jupyter (`DYLD_LIBRARY_PATH` on macOS, `PATH` on Windows), or add an `env` block to the kernel spec (`jupyter kernelspec list` shows where `octave/kernel.json` is):

```json
"env": {"LD_LIBRARY_PATH": "/path/to/meshioplusplus/inst-capi/lib"}
```

## Running them

The notebooks are committed **with their outputs** so they render on GitHub without any setup. To re-run:

```sh
export LD_LIBRARY_PATH=$PWD/inst-capi/lib MESHIOPLUSPLUS_HOME=$PWD/inst-capi
jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.kernel_name=octave \
  --ExecutePreprocessor.timeout=600 example/octave/*.ipynb
```

nbconvert sets the kernel's working directory to the notebook's own folder, which is why the cells use paths relative to `example/octave/`. The octave kernel reports an Octave error as cell output rather than failing the run, so check the executed notebooks for `error:` before committing them.

## In MATLAB

The cells and helpers keep to the language subset Octave and MATLAB share (`%` comments, single quotes, `~=`), so they run in MATLAB too: open them with the [MATLAB kernel for Jupyter](https://github.com/mathworks/jupyter-matlab-proxy), or paste the cells into a live script. Compile the gateway once with `mio_build` (it calls `mex` there) and set the library path before MATLAB starts. MATLAB's default colormap is parula, not Octave's viridis, which changes the look of the `cell_ids` picture in `04_terminal_rendering.ipynb`; the rasterizer's own colour maps are the same everywhere.

## Differences from the R notebooks

- Arrays put one entity per row: `mio.points(m)` is `num_points x 3` where R's is `3 x num_points`, so R's `mio_points(m)[3, ]` is `p = mio.points(m); p(:, 3)` here, and the bump, shuffle and field computations are column operations with no loop.
- Renders are coloured by the metric or field the C++ notebook colours by, because the rasterizer takes that option and the R notebooks' SVG route does not. The R notebooks' charts are kept alongside.
- `mio.split(..., 'by', 'region')` works with a tag built here: MATLAB has integer classes and the setters keep them, where every R setter writes `Float64`. `03_mesh_operations.ipynb` splits the surface by an `int32` tag next to the R notebook's split by type.
- The jittered block in the smoothing demo uses `rng(3)`, whose stream differs from R's `set.seed(3)`, so its numbers differ from the R notebook's; the outcome (a tangled block, minimum scaled Jacobian below zero, recovered to about 0.9) does not.

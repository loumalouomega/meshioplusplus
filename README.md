<p align="center">
  <a href="https://github.com/loumalouomega/meshioplusplus"><img alt="meshio++" src="https://raw.githubusercontent.com/loumalouomega/meshioplusplus/main/logo/logo-with-text.svg" width="60%"></a>
  <p align="center">I/O for mesh files.</p>
</p>

[![PyPi Version](https://img.shields.io/pypi/v/meshioplusplus.svg?style=flat-square)](https://pypi.org/project/meshioplusplus/) [![npm Version](https://img.shields.io/npm/v/%40meshioplusplus%2Fwasm.svg?style=flat-square)](https://www.npmjs.com/package/@meshioplusplus/wasm) [![PyPI pyversions](https://img.shields.io/pypi/pyversions/meshioplusplus.svg?style=flat-square)](https://pypi.org/project/meshioplusplus/) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21384760.svg?style=flat-square)](https://doi.org/10.5281/zenodo.21384760)

[![GitHub stars](https://img.shields.io/github/stars/loumalouomega/meshioplusplus.svg?style=flat-square&logo=github&label=Stars&logoColor=white)](https://github.com/loumalouomega/meshioplusplus) [![PyPi downloads](https://img.shields.io/pypi/dm/meshioplusplus.svg?style=flat-square)](https://pypistats.org/packages/meshioplusplus)

[![gh-actions](https://img.shields.io/github/actions/workflow/status/loumalouomega/meshioplusplus/ci.yml?branch=main&style=flat-square)](https://github.com/loumalouomega/meshioplusplus/actions?query=workflow%3Aci) [![codecov](https://img.shields.io/codecov/c/github/loumalouomega/meshioplusplus.svg?style=flat-square)](https://app.codecov.io/gh/loumalouomega/meshioplusplus) [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg?style=flat-square)](https://github.com/psf/black)

There are various mesh formats available for representing unstructured meshes. meshio++ can read and write all of the following and smoothly converts between them:

> [Abaqus](https://help.3ds.com/2024/english/dssimulia_established/SIMACAEMODRefMap/simamod-c-inputsyntax.htm) (`.inp`),
> ANSYS msh (`.msh`),
> [Ansys/APDL coded database](https://www.ansys.com) (`.cdb`, `.inp`),
> [AVS-UCD](https://lanl.github.io/LaGriT/pages/docs/read_avs.html) (`.avs`),
> [CGNS](https://cgns.github.io/) (`.cgns`),
> [DOLFIN XML](https://manpages.ubuntu.com/manpages/jammy/en/man1/dolfin-convert.1.html) (`.xml`),
> [COMSOL](https://www.comsol.com) (`.mphtxt`),
> [Exodus](https://nschloe.github.io/meshio/exodus.pdf) (`.e`, `.exo`),
> [EnSight Gold](https://vis.lbl.gov/archive/NERSC/Software/ensight/doc/OnlineHelp/UM-C11.pdf) (geometry, `.case`/`.geo`),
> [FLAC3D](https://www.itascacg.com/software/flac3d) (`.f3grid`),
> [FLUX](https://www.altair.com/flux/) (mesh `.pf3`, field `.dex`),
> [FreeFem++](https://freefem.org/) (`.msh`),
> [H5M](https://www.mcs.anl.gov/~fathom/moab-docs/h5mmain.html) (`.h5m`),
> [HMF](https://loumalouomega.github.io/meshioplusplus/formats/hmf) (`.hmf`, experimental, meshio++-specific),
> [I-deas Universal / UNV](https://www.ceas3.uc.edu/sdrluff/) (`.unv`),
> [ANSYS Fluent interpolation](https://github.com/victorsndvg/FEconv) (`.ip`),
> [Kratos/MDPA](https://github.com/KratosMultiphysics/Kratos/wiki/Input-data) (`.mdpa`),
> [Medit](https://people.sc.fsu.edu/~jburkardt/data/medit/medit.html) (`.mesh`, `.meshb`),
> [MED/Salome](https://docs.salome-platform.org/latest/dev/MEDCoupling/developer/med-file.html) (`.med`),
> [Modulef](https://github.com/victorsndvg/FEconv) (mesh `.mfm`, field `.mff`),
> [Nastran](https://help.autodesk.com/view/NSTRN/2019/ENU/?guid=GUID-42B54ACB-FBE3-47CA-B8FE-475E7AD91A00) (bulk data, `.bdf`, `.fem`, `.nas`),
> [Netgen](https://github.com/ngsolve/netgen) (`.vol`, `.vol.gz`),
> [Neuroglancer precomputed format](https://github.com/google/neuroglancer/tree/master/src/datasource/precomputed#mesh-representation-of-segmented-object-surfaces),
> [Gmsh](https://gmsh.info/doc/texinfo/gmsh.html#File-formats) (format versions 2.2, 4.0, and 4.1, `.msh`),
> [OBJ](https://en.wikipedia.org/wiki/Wavefront_.obj_file) (`.obj`),
> [OFF](https://segeval.cs.princeton.edu/public/off_format.html) (`.off`),
> [OpenFOAM polyMesh](https://www.openfoam.com/) (`.foam`, read-only),
> [PERMAS](https://www.intes.de) (`.post`, `.post.gz`, `.dato`, `.dato.gz`),
> [PLY](<https://en.wikipedia.org/wiki/PLY_(file_format)>) (`.ply`),
> [STL](<https://en.wikipedia.org/wiki/STL_(file_format)>) (`.stl`),
> [Tecplot .dat](http://paulbourke.net/dataformats/tp/),
> [TetGen .node/.ele](https://wias-berlin.de/software/tetgen/fformats.html),
> [Triangle .node/.ele/.poly](https://www.cs.cmu.edu/~quake/triangle.html),
> [SVG](https://www.w3.org/TR/SVG/) (2D output only) (`.svg`),
> [TikZ](https://tikz.dev/) (2D LaTeX output only) (`.tikz`),
> [SU2](https://su2code.github.io/docs_v7/Mesh-File/) (`.su2`),
> [UGRID](https://www.simcenter.msstate.edu/software/documentation/ug_io/3d_grid_file_type_ugrid.html) (`.ugrid`),
> [VTK](https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf) (`.vtk`),
> [VTP](https://docs.vtk.org/en/latest/vtk_file_formats/vtkxml_file_format.html) (`.vtp`),
> [VTU](https://vtk.org/Wiki/VTK_XML_Formats) (`.vtu`),
> [WKT](https://en.wikipedia.org/wiki/Well-known_text_representation_of_geometry) ([TIN](https://en.wikipedia.org/wiki/Triangulated_irregular_network)) (`.wkt`),
> [XDMF](https://xdmf.org/index.php/XDMF_Model_and_Format) (`.xdmf`, `.xmf`).

meshio++ ships a **C++20 core** (built with pybind11 + scikit-build-core) that reads and writes most formats with zero-copy numpy at the I/O boundary, plus optional HDF5/netCDF acceleration and a **selectable parallel backend** (`AUTO` by default — prefers OpenMP, then STL+TBB, then sequential; override with `-DMESHIOPLUSPLUS_PARALLEL_BACKEND=...`). Every format has a pure-Python fallback, so behaviour and file compatibility are identical whether or not the native libraries are present. For a standalone C++ build use `build/configure.sh` (Linux/macOS) or `build/configure.bat` (Windows). Full docs (install, data model, per-format options, CLI) live at [the documentation site](https://loumalouomega.github.io/meshioplusplus/) (sources under [`doc/`](https://github.com/loumalouomega/meshioplusplus/tree/main/doc)).

Install with

```
pip install meshioplusplus[all]
```

(`[all]` pulls in all optional dependencies. By default, meshio++ only uses numpy.) You can then use the command-line tool

<!--pytest-codeblocks:skip-->

```sh
meshioplusplus convert    input.msh output.vtk   # convert between two formats

meshioplusplus info       input.xdmf             # show some info about the mesh

meshioplusplus compress   input.vtu              # compress the mesh file
meshioplusplus decompress input.vtu              # decompress the mesh file

meshioplusplus binary     input.msh              # convert to binary format
meshioplusplus ascii      input.msh              # convert to ASCII format
```

with any of the supported formats.

In Python, simply do

<!--pytest-codeblocks:skip-->

```python
import meshioplusplus

mesh = meshioplusplus.read(
    filename,  # string, os.PathLike, or a buffer/open file
    # file_format="stl",  # optional if filename is a path; inferred from extension
    # see meshioplusplus convert --help for all possible formats
)
# mesh.points, mesh.cells, mesh.cells_dict, ...

# mesh.vtk.read() is also possible
```

to read a mesh. To write, do

```python
import meshioplusplus

# two triangles and one quad
points = [
    [0.0, 0.0],
    [1.0, 0.0],
    [0.0, 1.0],
    [1.0, 1.0],
    [2.0, 0.0],
    [2.0, 1.0],
]
cells = [
    ("triangle", [[0, 1, 2], [1, 3, 2]]),
    ("quad", [[1, 4, 5, 3]]),
]

mesh = meshioplusplus.Mesh(
    points,
    cells,
    # Optionally provide extra data on points, cells, etc.
    point_data={"T": [0.3, -1.2, 0.5, 0.7, 0.0, -3.0]},
    # Each item in cell data must match the cells array
    cell_data={"a": [[0.1, 0.2], [0.4]]},
)
mesh.write(
    "foo.vtk",  # str, os.PathLike, or buffer/open file
    # file_format="vtk",  # optional if first argument is a path; inferred from extension
)

# Alternative with the same options
meshioplusplus.write_points_cells("foo.vtk", points, cells)
```

For both input and output, you can optionally specify the exact `file_format` (in case you would like to enforce ASCII over binary VTK, for example).

#### Time series

The [XDMF format](https://xdmf.org/index.php/XDMF_Model_and_Format) supports time series with a shared mesh. You can write times series data using meshio++ with

<!--pytest-codeblocks:skip-->

```python
with meshioplusplus.xdmf.TimeSeriesWriter(filename) as writer:
    writer.write_points_cells(points, cells)
    for t in [0.0, 0.1, 0.21]:
        writer.write_data(t, point_data={"phi": data})
```

and read it with

<!--pytest-codeblocks:skip-->

```python
with meshioplusplus.xdmf.TimeSeriesReader(filename) as reader:
    points, cells = reader.read_points_cells()
    for k in range(reader.num_steps):
        t, point_data, cell_data = reader.read_data(k)
```

### ParaView plugin

<img alt="gmsh paraview" src="https://nschloe.github.io/meshio/gmsh-paraview.png" width="60%">
*A Gmsh file opened with ParaView.*

If you have downloaded a binary version of ParaView, you may proceed as follows.

- Install meshio++ for the Python major version that ParaView uses (check `pvpython --version`)
- Open ParaView
- Find the file `paraview-meshioplusplus-plugin.py` of your meshio++ installation (on Linux: `~/.local/share/paraview-5.9/plugins/`) and load it under _Tools / Manage Plugins / Load New_
- _Optional:_ Activate _Auto Load_

You can now open all meshio++-supported files in ParaView.

### Benchmarks

How much does the C++ core help? The [`benchmark/`](https://github.com/loumalouomega/meshioplusplus/tree/main/benchmark) folder times read/write conversions against the original pure-Python [meshio](https://github.com/nschloe/meshio) on the formats both support (same in-memory mesh, same machine). The headline input is the bundled [`example.msh`](https://github.com/loumalouomega/meshioplusplus/blob/main/example/example.msh) — a real Gmsh bracket (~52k nodes, ~293k cells).

<img alt="meshio vs meshio++ speedup on example.msh" src="https://raw.githubusercontent.com/loumalouomega/meshioplusplus/main/benchmark/plots/benchmark_speedup.svg" width="85%">

meshio++'s biggest wins are the parallel and text paths: **VTU binary+zlib ~16× write** (the zlib blocks run across cores via an OpenMP backend with dynamic scheduling — hybrid P+E-core CPUs load-balance too), **VTU ASCII ~7× write / ~5× read**, and mixed-topology **XDMF read ~10×**. The binary and HDF5 formats that used to be *slower* — VTK/Gmsh binary, UGRID, and MED — are now at or above parity after an optimisation pass (bulk-buffered binary I/O, single-instruction `bswap` endianness conversion, a real parallel backend, an Eigen-backed MED transpose, **zero-copy cell reconstruction** that moves the connectivity buffer straight into the mesh, and uninitialised reader buffers + thread-parallel block copies so nothing is written twice); binary **reads** now match or beat numpy's `fromfile` — Gmsh ~1.7×, single-type VTK ~1.45×, and even mixed-topology VTK ~1.1×. Output stays byte-identical throughout.

The speedup is per-element: text/parallel formats climb out of the small-mesh regime and plateau (large meshes realise the full speedup):

<img alt="speedup vs mesh size" src="https://raw.githubusercontent.com/loumalouomega/meshioplusplus/main/benchmark/plots/benchmark_scaling.svg" width="85%">

Full methodology and a reproducible notebook are on the [Benchmarks](https://loumalouomega.github.io/meshioplusplus/benchmarks) doc page (source: [`benchmark/01_benchmark.ipynb`](https://github.com/loumalouomega/meshioplusplus/blob/main/benchmark/01_benchmark.ipynb)).

### Installation

meshio++ is [available from the Python Package Index](https://pypi.org/project/meshioplusplus/), so simply run

```
pip install meshioplusplus
```

to install.

Additional dependencies (`netcdf4`, `h5py`) are required for some of the output formats and can be pulled in by

```
pip install meshioplusplus[all]
```

For JavaScript / browser use, the C++ core also ships as a WebAssembly npm package covering 29 of the formats above:

```
npm install @meshioplusplus/wasm
```

See the [WebAssembly / JavaScript](https://loumalouomega.github.io/meshioplusplus/wasm) doc page for usage and the format-support table.

### C / Fortran API

For HPC codes written in C or Fortran, the C++ core also builds as an installable shared library (`libmeshioplusplus`, pure-C99 header, pkg-config + `find_package` support) with a modern OO Fortran 2008 module on top:

```
./build/configure.sh --fortran --tests --build     # --c-api for the C API alone
cmake --install build/cpp-release --prefix /opt/meshioplusplus
```

```c
mio_mesh* m = mio_read("in.msh", NULL);
printf("%lld points\n", (long long)mio_mesh_num_points(m));
mio_write("out.vtu", m, NULL);
mio_mesh_free(m);
```

```fortran
use meshioplusplus
type(mio_mesh) :: m
call m%read("in.msh")
call m%write("out.vtu")
call m%free()
```

The C API is also packaged for **Conan** (root [`conanfile.py`](conanfile.py)) and **vcpkg** (overlay port under [`ports/meshioplusplus/`](ports/meshioplusplus)), both driving the same install/`find_package` path:

```
conan create . -o meshioplusplus/*:with_hdf5=True
vcpkg install meshioplusplus --overlay-ports=ports
```

Full mesh access (build meshes from raw arrays, zero-copy readback) is covered on the [C API](https://loumalouomega.github.io/meshioplusplus/c_api) and [Fortran](https://loumalouomega.github.io/meshioplusplus/fortran) doc pages.

### Single-header C++

The whole C++ core is also amalgamated into one self-contained, [STB](https://github.com/nothings/stb)-style header — [`single_include/meshioplusplus/meshioplusplus.hpp`](single_include/meshioplusplus/meshioplusplus.hpp) — with pugixml bundled and no external dependencies by default. Drop it in, no CMake or linking required:

```cpp
// in exactly ONE .cpp:
#define MESHIOPLUSPLUS_IMPLEMENTATION
#include "meshioplusplus/meshioplusplus.hpp"
// elsewhere: just #include it (declarations only)
```

```
g++ -std=c++20 -I single_include main.cpp
```

It is generated by `./tools/amalgamate.sh` and kept in sync by CI. See the [single-header](https://loumalouomega.github.io/meshioplusplus/single_header) doc page (optional HDF5/netCDF/zlib formats via `MESHIOPLUSPLUS_HAS_*` macros).

### C++ mesh backends

Standalone C++ builds (no Python) can swap the in-memory mesh structure at compile time via `MESHIOPLUSPLUS_MESH_BACKEND` — every format works identically under each backend:

- **MESHIO** (default; the Python extension and PyPI wheels always use it) — mirrors the Python `meshio.Mesh`;
- **NATIVE** — the fastest pure-C++ structure (canonical Float64/Int64 storage, cell-type enum, CSR ragged blocks); the WebAssembly build uses it;
- **KRATOS** — a [Kratos Multiphysics](https://github.com/KratosMultiphysics/Kratos)-style `ModelPart` (Nodes/Elements/Conditions/SubModelParts) plus a header-only templated bridge that populates a real `Kratos::ModelPart` with no Kratos build dependency.

```
./build/configure.sh --mesh-backend NATIVE --tests --build
```

See the [C++ mesh backends](https://loumalouomega.github.io/meshioplusplus/cpp_backends) doc page.

### Testing

To run the meshio++ unit tests, check out this repository, install it with the test extras, and type

```
pytest tests/
```

### License

meshio++ is published under the [MIT license](https://en.wikipedia.org/wiki/MIT_License).

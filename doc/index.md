---
layout: home

hero:
  name: meshio++
  text: I/O and operations for many mesh formats
  tagline: One unified mesh data model, 43 file formats, 39 mesh and data operations, a fast C++ core with pure-Python fallbacks, and six language surfaces over it.
  image:
    src: /logo-icon.svg
    alt: meshio++
  actions:
    - theme: brand
      text: Quickstart
      link: /quickstart
    - theme: alt
      text: Architecture
      link: /architecture
    - theme: alt
      text: Supported formats
      link: /formats
    - theme: alt
      text: GitHub
      link: https://github.com/loumalouomega/meshioplusplus

features:
  - icon: 🗂️
    title: 43 formats
    details: Read 43 and write 46 mesh formats — VTK, VTU, XDMF, Gmsh, MED, Exodus, CGNS, Abaqus, Nastran, UNV, OpenFOAM, GiD, COMSOL, FLUX and more — through a single API, with per-format options where the format has them.
    link: /formats
  - icon: 🧩
    title: One data model
    details: A single Mesh (points, cell blocks, point and cell data, field data, named regions of points, cells or sides) bridges every format, so conversion is one call and a group survives the round trip.
    link: /mesh_data_model
  - icon: ⚡
    title: Fast C++ core, Python fallbacks
    details: A C++20 core behind a pybind11 extension with zero-copy numpy at the I/O boundary, optional HDF5, netCDF and compression codecs, and a pure-Python reference implementation for every format.
    link: /architecture
  - icon: 🛠️
    title: 34 mesh + 5 data operations
    details: Quality, skin and surface extraction, reordering, cleaning, cropping, splitting, partitioning, refinement and coarsening, decimation, remeshing, smoothing, slicing, isosurfaces, gradients, interpolation and more, on every surface.
    link: /architecture#the-operations-layer
  - icon: 🌐
    title: Six language surfaces
    details: Python, a pure-C99 API with a versioned ABI, Fortran, Julia, R and a WebAssembly package, plus an installable C++ API and a single-header build, all over the same core and the same format registry.
    link: /c_api
  - icon: ⌨️
    title: Two CLIs and a pipeline engine
    details: A Python CLI and a Python-free native binary with the same verbs, a declarative settings.json pipeline that every surface runs, and a sequence driver for transient, multi-file datasets.
    link: /cli
  - icon: 👁️
    title: Viewers, Blender, ParaView
    details: A browser viewer and dataset manager on the WebAssembly build, a Polyscope desktop viewer, a Blender 4.2+ add-on and a ParaView reader plugin.
    link: /viewer
  - icon: 🤖
    title: Built for agents and ML
    details: An MCP server exposing every operation as a stateless tool, zero-copy interop with PyVista, trimesh, Arrow and pandas, DLPack and CuPy handoff, graph and feature-matrix export, and a PhysicsNeMo adapter.
    link: /mcp
---

## What is meshio++?

meshio++ reads and writes unstructured mesh files. It supports 43 formats through one unified [data model](./mesh_data_model.md), so you can convert between any of them and run the same operations on the result, from the command line or from Python:

```python
import meshioplusplus

mesh = meshioplusplus.read("input.msh")   # format inferred from the extension
mesh.write("output.vtu")
```

```sh
meshioplusplus convert input.msh output.vtu
meshioplusplus info input.xdmf
```

Under the Python API sits one C++ core with six language surfaces on top of it, and the tools — CLIs, viewers, an MCP server, a Blender add-on — are built on those surfaces:

![meshio++ architecture: one C++ core, six language surfaces, and the tools built on them](/diagrams/architecture.svg)

See the [Quickstart](./quickstart.md) to get going, the [Architecture](./architecture.md) page for the map, the [Supported formats](./formats.md) table for the full list and per-format options, and the [CLI reference](./cli.md) for the command-line tools.

---
layout: home

hero:
  name: meshio++
  text: I/O for many mesh formats
  tagline: One unified mesh data model, 35+ file formats, a fast C++ core with pure-Python fallbacks.
  image:
    src: /logo-icon.svg
    alt: meshio++
  actions:
    - theme: brand
      text: Quickstart
      link: /quickstart
    - theme: alt
      text: Supported formats
      link: /formats
    - theme: alt
      text: GitHub
      link: https://github.com/loumalouomega/meshioplusplus

features:
  - title: 35+ formats
    details: Read and write VTK, VTU, XDMF, Gmsh, MED, Exodus, CGNS, Abaqus, Nastran, UNV, COMSOL, FLUX, and many more — all through a single API.
  - title: Unified data model
    details: A single Mesh object (points, cells, point/cell data, field data, sets) bridges every format, so conversion is one call.
  - title: Fast C++ core
    details: A pybind11 extension with zero-copy numpy at the I/O boundary and optional HDF5/netCDF acceleration — with pure-Python fallbacks everywhere.
  - title: Swappable mesh backends
    details: Standalone C++ builds choose the in-memory structure at compile time — the meshio-mirroring default, a fastest-possible native layout (used by the WASM build), or a Kratos Multiphysics-style ModelPart with a header-only bridge.
---

## What is meshio++?

meshio++ reads and writes unstructured mesh files. It supports over 35 formats and provides one unified [data model](./mesh_data_model.md) so you can convert between any of them, from the command line or from Python:

```python
import meshioplusplus

mesh = meshioplusplus.read("input.msh")   # format inferred from the extension
mesh.write("output.vtu")
```

```sh
meshioplusplus convert input.msh output.vtu
meshioplusplus info input.xdmf
```

See the [Quickstart](./quickstart.md) to get going, the [Supported formats](./formats.md) table for the full list and per-format options, and the [CLI reference](./cli.md) for the command-line tools.

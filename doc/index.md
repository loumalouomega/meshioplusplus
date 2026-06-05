# meshio Documentation

meshio is a Python library for reading and writing unstructured mesh files. It supports over 30 formats and provides a single unified data model so you can convert between any of them.

## Contents

- [Installation](installation.md) — install options, optional dependencies
- [Quickstart](quickstart.md) — reading, writing, and converting meshes in a few lines
- [Mesh data model](mesh_data_model.md) — `Mesh`, `CellBlock`, points, cells, data fields, sets
- [Cell types](cell_types.md) — all supported element types and their node counts
- [Supported formats](formats.md) — full table of formats, extensions, read/write support, and format-specific options
- [CLI reference](cli.md) — `meshio convert`, `info`, `ascii`, `binary`, `compress`, `decompress`
- [XDMF time series](xdmf_time_series.md) — writing and reading temporal simulation data
- [Extending meshio](extending.md) — registering custom formats, adding a new built-in format
- [ParaView plugin](paraview_plugin.md) — loading meshio-supported files directly in ParaView

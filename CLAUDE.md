# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this project is

meshio is a Python library for reading and writing many mesh file formats used in scientific computing and FEM (Finite Element Method). It provides a unified `Mesh` data structure that all format readers/writers convert to and from.

## Commands

**Install for development:**
```bash
pip install -e ".[all]"
```

**Run all tests:**
```bash
pytest tests/
```

**Run a single test file:**
```bash
pytest tests/test_gmsh.py
```

**Run a specific test:**
```bash
pytest tests/test_gmsh.py::test_gmsh22
```

**Lint:**
```bash
black --check .
flake8 .
isort --check .
```

**Format:**
```bash
isort .
black .
```

**Via tox:**
```bash
tox           # runs pytest
tox -e lint   # lint check
tox -e fmt    # reformat
```

**CLI usage:**
```bash
meshio convert input.msh output.vtk
meshio info mesh.msh
meshio ascii mesh.msh         # convert to ASCII
meshio binary mesh.msh        # convert to binary
```

## Architecture

**Core data model** (`src/meshio/_mesh.py`):
- `Mesh`: central class holding `points` (numpy array), `cells` (list of `CellBlock`), `point_data`, `cell_data`, `field_data`, `point_sets`, `cell_sets`, and optionally `gmsh_periodic`
- `CellBlock`: holds a cell type string (e.g. `"triangle"`, `"tetra10"`) and a numpy array of node indices
- `_common.py`: `num_nodes_per_cell` dict maps cell type names to node counts; this is the canonical registry of supported element types

**Format registration** (`src/meshio/_helpers.py`):
- Each format module calls `register_format(name, extensions, reader_fn, writer_map)` at module import time
- `read(filename)` auto-detects format from file extension via `extension_to_filetypes`; `write(filename, mesh)` does the same
- Format modules are imported in `__init__.py`, which triggers registration

**Format module layout** (e.g. `src/meshio/gmsh/`):
- `__init__.py` exports `read`, `write`, and any type-mapping helpers
- `main.py` (or `_<format>.py`) implements `read(filename)` and `write(filename, mesh)`, calls `register_format` at the bottom
- Multi-version formats (Gmsh, VTK) dispatch to version-specific submodules

**Cell type naming convention**: meshio uses its own type names (`triangle`, `tetra10`, `hexahedron20`, etc.). Each format module has a type mapping between its native names and meshio names (e.g. `_gmsh_to_meshio_type`, `_meshio_to_gmsh_type` in `gmsh/common.py`).

**Tests** (`tests/`):
- `helpers.py` contains pre-built `Mesh` fixtures (`tri_mesh`, `tet_mesh`, `hex_mesh`, etc.) and `write_read()` — the standard round-trip test pattern
- Per-format test files (`test_gmsh.py`, `test_vtu.py`, etc.) parametrize over mesh types and call `helpers.write_read(tmp_path, writer, reader, mesh, atol)`
- `tests/input/` and `tests/meshes/` contain real mesh files used in read-only tests

**CLI** (`src/meshio/_cli/`): subcommands (`convert`, `info`, `ascii`, `binary`, `compress`, `decompress`) each live in their own module and are registered in `_main.py`.

## Adding a new format

1. Create `src/meshio/<format>/` with `__init__.py`, implement `read()` and `write()`, call `register_format()` at module level
2. Import the module in `src/meshio/__init__.py`
3. Add type mappings in a `common.py` if the format uses different cell type names
4. Add `tests/test_<format>.py` using `helpers.write_read()`

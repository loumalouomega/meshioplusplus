# Supported Formats

## Format table

| Format name | Extensions | Read | Write | Extra dependencies |
|-------------|-----------|------|-------|--------------------|
| `abaqus` | `.inp` | ✓ | ✓ | — |
| `ansys` | `.msh` | ✓ | ✓ | — |
| `avsucd` | `.avs` | ✓ | ✓ | — |
| `cgns` | `.cgns` | ✓ | ✓ | `h5py` |
| `dolfin-xml` | `.xml` | ✓ | ✓ | — |
| `exodus` | `.e`, `.exo`, `.ex2` | ✓ | ✓ | `netCDF4` |
| `flac3d` | `.f3grid` | ✓ | ✓ | — |
| `gmsh` / `gmsh22` | `.msh` | ✓ | ✓ | — |
| `h5m` | `.h5m` | ✓ | ✓ | `h5py` |
| `hmf` | `.hmf` | ✓ | ✓ | `h5py` |
| `mdpa` | `.mdpa` | ✓ | ✓ | — |
| `med` | `.med` | ✓ | ✓ | `h5py` |
| `medit` | `.mesh`, `.meshb` | ✓ | ✓ | — |
| `nastran` | `.bdf`, `.fem`, `.nas` | ✓ | ✓ | — |
| `netgen` | `.vol`, `.vol.gz` | ✓ | ✓ | — |
| `neuroglancer` | (no extension) | ✓ | ✓ | — |
| `obj` | `.obj` | ✓ | ✓ | — |
| `off` | `.off` | ✓ | ✓ | — |
| `permas` | `.post`, `.post.gz`, `.dato`, `.dato.gz` | ✓ | ✓ | — |
| `ply` | `.ply` | ✓ | ✓ | — |
| `stl` | `.stl` | ✓ | ✓ | — |
| `su2` | `.su2` | ✓ | ✓ | — |
| `svg` | `.svg` | — | ✓ | — |
| `tecplot` | `.dat`, `.tec` | ✓ | ✓ | — |
| `tetgen` | `.ele` / `.node` | ✓ | ✓ | — |
| `ugrid` | `.ugrid` | ✓ | ✓ | — |
| `vtk` / `vtk42` / `vtk51` | `.vtk` | ✓ | ✓ | — |
| `vtu` | `.vtu` | ✓ | ✓ | — |
| `wkt` | `.wkt` | ✓ | ✓ | — |
| `xdmf` | `.xdmf`, `.xmf` | ✓ | ✓ | `h5py` (for HDF data) |

**Note on `.msh`:** Both `ansys` and `gmsh` use `.msh`. When reading, meshio tries `ansys` first, then `gmsh`. Specify `file_format` explicitly to avoid ambiguity.

**Note on `tetgen`:** The format spans two files (`.node` + `.ele`). It cannot be read from or written to a buffer.

**Note on `svg`:** Write-only, 2D meshes only.

---

## Format-specific write options

All writers are called as `meshio.write(filename, mesh, file_format=..., **kwargs)` or `mesh.write(filename, **kwargs)`. The `**kwargs` depend on the format.

### Gmsh (`.msh`)

```python
meshio.gmsh.write(filename, mesh,
    fmt_version="4.1",   # "2.2", "4.0", or "4.1"
    binary=True,
    float_fmt=".16e",
)
```

Use `file_format="gmsh22"` to write version 2.2 via the generic `meshio.write`.

### VTU (`.vtu`)

```python
meshio.vtu.write(filename, mesh,
    binary=True,
    compression="zlib",   # "zlib", "lzma", or None
    header_type=None,     # "UInt32" or "UInt64"
)
```

### VTK (`.vtk`)

```python
meshio.vtk.write(filename, mesh,
    binary=True,
    # For version selection use file_format="vtk42" or "vtk51"
)
```

`file_format="vtk"` writes VTK 5.1. `file_format="vtk42"` or `"vtk51"` select specific versions.

### XDMF (`.xdmf`, `.xmf`)

```python
meshio.xdmf.write(filename, mesh,
    data_format="HDF",        # "HDF", "XML", or "Binary"
    compression="gzip",       # h5py compression filter (HDF only)
    compression_opts=4,       # compression level
)
```

With `data_format="HDF"`, meshio writes a companion `.h5` file alongside the `.xdmf`. With `"XML"`, all data is embedded in the XML. With `"Binary"`, data is written to separate `.bin` files.

### Medit (`.mesh`)

```python
meshio.medit.write(filename, mesh,
    float_fmt=".16e",
)
```

### PLY (`.ply`)

```python
meshio.ply.write(filename, mesh,
    binary=True,
)
```

### STL (`.stl`)

```python
meshio.stl.write(filename, mesh,
    binary=False,
)
```

### MED (`.med`)

```python
meshio.med.write(filename, mesh,
    compression="gzip",
    compression_opts=4,
)
```

### CGNS (`.cgns`)

```python
meshio.cgns.write(filename, mesh,
    compression="gzip",
    compression_opts=4,
)
```

### Nastran (`.bdf`)

```python
meshio.nastran.write(filename, mesh,
    point_format="fixed-large",   # or "fixed-small", "free"
    cell_format="fixed-small",
)
```

### FLAC3D (`.f3grid`)

```python
meshio.flac3d.write(filename, mesh,
    float_fmt=".16e",
    binary=False,
)
```

### SU2 (`.su2`)

`meshio.su2.write(filename, mesh)` — no extra options.

### AVS-UCD (`.avs`)

`meshio.avsucd.write(filename, mesh)` — no extra options.

### Abaqus (`.inp`)

`meshio.abaqus.write(filename, mesh)` — no extra options.

### DOLFIN-XML (`.xml`)

`meshio.dolfin.write(filename, mesh)` — no extra options.

---

## CLI format names

When using `meshio convert -o <format>`, use one of the format names from the first column of the table above (e.g. `gmsh`, `gmsh22`, `vtk`, `vtk42`, `vtu`, `xdmf`, …).

# Supported Formats

## Format table

Each format name links to a detailed reference page (structure, options, data
mapping, and the C++ vs Python behaviour).

| Format name | Extensions | Read | Write | Extra dependencies |
|-------------|-----------|------|-------|--------------------|
| [`abaqus`](./formats/abaqus.md) | `.inp` | ✓ | ✓ | — |
| [`ansys`](./formats/ansys.md) | `.msh` | ✓ | ✓ | — |
| [`avsucd`](./formats/avsucd.md) | `.avs` | ✓ | ✓ | — |
| [`cgns`](./formats/cgns.md) | `.cgns` | ✓ | ✓ | `h5py` |
| [`dolfin-xml`](./formats/dolfin.md) | `.xml` | ✓ | ✓ | — |
| [`exodus`](./formats/exodus.md) | `.e`, `.exo`, `.ex2` | ✓ | ✓ | `netCDF4` |
| [`flac3d`](./formats/flac3d.md) | `.f3grid` | ✓ | ✓ | — |
| [`flux`](./formats/flux.md) | `.pf3` | ✓ | ✓ | — |
| [`freefem`](./formats/freefem.md) | `.msh` | ✓ | ✓ | — |
| [`gmsh` / `gmsh22`](./formats/gmsh.md) | `.msh` | ✓ | ✓ | — |
| [`h5m`](./formats/h5m.md) | `.h5m` | ✓ | ✓ | `h5py` |
| [`hmf`](./formats/hmf.md) | `.hmf` | ✓ | ✓ | `h5py` |
| [`mdpa`](./formats/mdpa.md) | `.mdpa` | ✓ | ✓ | — |
| [`med`](./formats/med.md) | `.med` | ✓ | ✓ | `h5py` |
| [`medit`](./formats/medit.md) | `.mesh`, `.meshb` | ✓ | ✓ | — |
| [`mfm`](./formats/mfm.md) | `.mfm` | ✓ | ✓ | — |
| [`mphtxt`](./formats/mphtxt.md) | `.mphtxt` | ✓ | ✓ | — |
| [`nastran`](./formats/nastran.md) | `.bdf`, `.fem`, `.nas` | ✓ | ✓ | — |
| [`netgen`](./formats/netgen.md) | `.vol`, `.vol.gz` | ✓ | ✓ | — |
| [`neuroglancer`](./formats/neuroglancer.md) | (no extension) | ✓ | ✓ | — |
| [`obj`](./formats/obj.md) | `.obj` | ✓ | ✓ | — |
| [`off`](./formats/off.md) | `.off` | ✓ | ✓ | — |
| [`permas`](./formats/permas.md) | `.post`, `.post.gz`, `.dato`, `.dato.gz` | ✓ | ✓ | — |
| [`ply`](./formats/ply.md) | `.ply` | ✓ | ✓ | — |
| [`stl`](./formats/stl.md) | `.stl` | ✓ | ✓ | — |
| [`su2`](./formats/su2.md) | `.su2` | ✓ | ✓ | — |
| [`svg`](./formats/svg.md) | `.svg` | — | ✓ | — |
| [`tecplot`](./formats/tecplot.md) | `.dat`, `.tec` | ✓ | ✓ | — |
| [`tetgen`](./formats/tetgen.md) | `.ele` / `.node` | ✓ | ✓ | — |
| [`ugrid`](./formats/ugrid.md) | `.ugrid` | ✓ | ✓ | — |
| [`unv`](./formats/unv.md) | `.unv` | ✓ | ✓ | — |
| [`vtk` / `vtk42` / `vtk51`](./formats/vtk.md) | `.vtk` | ✓ | ✓ | — |
| [`vtu`](./formats/vtu.md) | `.vtu` | ✓ | ✓ | — |
| [`wkt`](./formats/wkt.md) | `.wkt` | ✓ | ✓ | — |
| [`xdmf`](./formats/xdmf.md) | `.xdmf`, `.xmf` | ✓ | ✓ | `h5py` (for HDF data) |

**Note on `.msh`:** `ansys`, `freefem`, and `gmsh` all use `.msh`. When reading, meshio tries them in registration order and uses the first that parses the file. Specify `file_format` explicitly (e.g. `file_format="freefem"`) to avoid ambiguity.

**Note on `tetgen`:** The format spans two files (`.node` + `.ele`). It cannot be read from or written to a buffer.

**Note on `svg`:** Write-only, 2D meshes only.

**Note on `mfm`:** Single element type per file (non-hybrid), linear elements only.

**Note on FEconv-derived formats (`unv`, `mfm`, `freefem`, `mphtxt`, `flux`):** These readers/writers were implemented against the [FEconv](https://github.com/victorsndvg/FEconv) format documentation. `unv` handles the parabolic mid-node "sandwich" ordering and maps permanent groups (dataset 2467) to `point_sets`/`cell_sets`; `mphtxt` and `flux` round-trip per-element region references as `cell_data` (`mphtxt:geom`, `pf3:ref`). Node orderings for higher-order elements round-trip losslessly but may differ from the originating tool's internal ordering for some element types.

---

## Native acceleration and fallbacks

meshio ships a C++ core (`meshio._core`, built with pybind11 + scikit-build-core). Most formats read and write through the C++ core with zero-copy numpy at the I/O boundary; each has a pure-Python fallback that is used automatically when the C++ path can't handle a file or when the extension was built without an optional dependency:

- **HDF5** (`cgns`, `h5m`, `hmf`, `med`, and XDMF `data_format="HDF"`) — C++ when built with `MESHIO_WITH_HDF5`, otherwise `h5py`.
- **netCDF** (`exodus`) — C++ when built with `MESHIO_WITH_NETCDF`, otherwise `netCDF4`.
- **zlib** (VTU zlib compression) — C++ when built with `MESHIO_WITH_ZLIB`, otherwise the Python stdlib.

Behaviour and file compatibility are identical either way; the native paths are only faster. Install the optional runtime deps with `pip install meshio[all]`.

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

`meshio.med.write(filename, mesh)` — no extra options (MED does not compress).

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

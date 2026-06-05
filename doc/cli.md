# CLI Reference

The `meshio` command-line tool is installed alongside the Python package.

```
meshio --version
meshio --help
meshio <subcommand> --help
```

---

## meshio convert

Convert a mesh file from one format to another.

```
meshio convert [options] INFILE OUTFILE
```

| Option | Short | Description |
|--------|-------|-------------|
| `--input-format FORMAT` | `-i` | Force input format (skip extension detection) |
| `--output-format FORMAT` | `-o` | Force output format |
| `--ascii` | `-a` | Write ASCII variant (default: binary where available) |
| `--float-format FMT` | `-f` | Float format string for ASCII output (default: `.16e`) |
| `--sets-to-int-data` | `-s` | Convert point/cell sets to integer data arrays |
| `--int-data-to-sets` | `-d` | Convert integer data arrays to point/cell sets |

**Examples:**

```sh
meshio convert mesh.msh mesh.vtu
meshio convert -i gmsh -o vtk mesh.msh mesh.vtk
meshio convert --ascii mesh.msh mesh.vtu
meshio convert --sets-to-int-data mesh.inp mesh.xdmf
```

---

## meshio info

Print a summary of a mesh file.

```
meshio info [options] INFILE
```

| Option | Short | Description |
|--------|-------|-------------|
| `--input-format FORMAT` | `-i` | Force input format |

Output includes: number of points, cell blocks and their types/counts, point/cell sets, point/cell data names, field data names. It also warns if cells reference nonexistent points or if there are unused points.

**Example:**

```sh
meshio info mesh.msh
```

---

## meshio compress

Compress the data in a mesh file (formats that support compression, e.g. VTU).

```
meshio compress [options] INFILE
```

| Option | Short | Description |
|--------|-------|-------------|
| `--input-format FORMAT` | `-i` | Force input format |

---

## meshio decompress

Decompress the data in a mesh file.

```
meshio decompress [options] INFILE
```

| Option | Short | Description |
|--------|-------|-------------|
| `--input-format FORMAT` | `-i` | Force input format |

---

## meshio ascii

Convert a mesh file to its ASCII representation (in-place).

```
meshio ascii [options] INFILE
```

| Option | Short | Description |
|--------|-------|-------------|
| `--input-format FORMAT` | `-i` | Force input format |

---

## meshio binary

Convert a mesh file to its binary representation (in-place).

```
meshio binary [options] INFILE
```

| Option | Short | Description |
|--------|-------|-------------|
| `--input-format FORMAT` | `-i` | Force input format |

---

## Format names

The `--input-format` and `--output-format` options accept any of the registered format names. The full list is shown by `meshio convert --help`. Common values:

`abaqus`, `ansys`, `avsucd`, `cgns`, `dolfin-xml`, `exodus`, `flac3d`, `gmsh`, `gmsh22`, `h5m`, `hmf`, `mdpa`, `med`, `medit`, `nastran`, `netgen`, `obj`, `off`, `permas`, `ply`, `stl`, `su2`, `svg`, `tecplot`, `tetgen`, `ugrid`, `vtk`, `vtk42`, `vtk51`, `vtu`, `wkt`, `xdmf`

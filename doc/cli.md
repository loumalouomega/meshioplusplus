# CLI Reference

The `meshioplusplus` command-line tool is installed alongside the Python package.

```
meshioplusplus --version
meshioplusplus --help
meshioplusplus <subcommand> --help
```

---

## meshioplusplus convert

Convert a mesh file from one format to another.

```
meshioplusplus convert [options] INFILE OUTFILE
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
meshioplusplus convert mesh.msh mesh.vtu
meshioplusplus convert -i gmsh -o vtk mesh.msh mesh.vtk
meshioplusplus convert --ascii mesh.msh mesh.vtu
meshioplusplus convert --sets-to-int-data mesh.inp mesh.xdmf
```

---

## meshioplusplus info

Print a summary of a mesh file.

```
meshioplusplus info [options] INFILE
```

| Option | Short | Description |
|--------|-------|-------------|
| `--input-format FORMAT` | `-i` | Force input format |

Output includes: number of points, cell blocks and their types/counts, point/cell sets, point/cell data names, field data names. It also warns if cells reference nonexistent points or if there are unused points.

**Example:**

```sh
meshioplusplus info mesh.msh
```

---

## meshioplusplus compress

Compress the data in a mesh file (formats that support compression, e.g. VTU).

```
meshioplusplus compress [options] INFILE
```

| Option | Short | Description |
|--------|-------|-------------|
| `--input-format FORMAT` | `-i` | Force input format |

---

## meshioplusplus decompress

Decompress the data in a mesh file.

```
meshioplusplus decompress [options] INFILE
```

| Option | Short | Description |
|--------|-------|-------------|
| `--input-format FORMAT` | `-i` | Force input format |

---

## meshioplusplus ascii

Convert a mesh file to its ASCII representation (in-place).

```
meshioplusplus ascii [options] INFILE
```

| Option | Short | Description |
|--------|-------|-------------|
| `--input-format FORMAT` | `-i` | Force input format |

---

## meshioplusplus binary

Convert a mesh file to its binary representation (in-place).

```
meshioplusplus binary [options] INFILE
```

| Option | Short | Description |
|--------|-------|-------------|
| `--input-format FORMAT` | `-i` | Force input format |

---

## Format names

The `--input-format` and `--output-format` options accept any of the registered format names. The full list is shown by `meshioplusplus convert --help`. Common values:

`abaqus`, `ansys`, `avsucd`, `cgns`, `dolfin-xml`, `exodus`, `flac3d`, `gmsh`, `gmsh22`, `h5m`, `hmf`, `mdpa`, `med`, `medit`, `nastran`, `netgen`, `obj`, `off`, `permas`, `ply`, `stl`, `su2`, `svg`, `tecplot`, `tetgen`, `ugrid`, `vtk`, `vtk42`, `vtk51`, `vtu`, `wkt`, `xdmf`

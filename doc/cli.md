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
meshioplusplus convert mesh.msh skin.stl   # volume mesh -> boundary-skin STL
```

Converting a 3D volume mesh to STL or PLY writes its extracted boundary
skin (the writers' default — see [Skin extraction](./extract_skin.md));
converting to SVG or TikZ renders it with the default isometric camera.

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

## meshioplusplus quality

Print a per-cell [mesh quality](./mesh_quality.md) report (min/mean/max and
counts of inverted/degenerate cells).

```
meshioplusplus quality [options] INFILE
```

| Option | Short | Description |
|--------|-------|-------------|
| `--input-format FORMAT` | `-i` | Force input format |
| `--output FILE` | `-o` | Also write the metrics into `FILE` as `cell_data` |

**Examples:**

```sh
meshioplusplus quality part.vtu
meshioplusplus quality part.vtu -o part_quality.vtu
```

---

## meshioplusplus extract-surface

Extract the [boundary surface/edges](./extract_surface.md) of a mesh (volume →
faces, 2D surface → edges) and write it out.

```
meshioplusplus extract-surface [options] INFILE OUTFILE
```

| Option | Short | Description |
|--------|-------|-------------|
| `--input-format FORMAT` | `-i` | Force input format |
| `--output-format FORMAT` | `-o` | Force output format |
| `--parent-ids` | `-p` | Record each facet's parent cell id as `cell_data` |

**Examples:**

```sh
meshioplusplus extract-surface part.vtu surface.stl
meshioplusplus extract-surface --parent-ids part.vtu surface.vtu
```

---

## meshioplusplus reorder

[Renumber](./reorder.md) a mesh's nodes/elements to reduce matrix bandwidth
(RCM) or improve cache locality (Morton / Hilbert), as a pure permutation.

```
meshioplusplus reorder [options] INFILE OUTFILE
```

| Option | Short | Description |
|--------|-------|-------------|
| `--method METHOD` | `-m` | `rcm` (default), `morton`, or `hilbert` |
| `--report` | `-r` | Print the connectivity bandwidth before and after |
| `--input-format FORMAT` | `-i` | Force input format |
| `--output-format FORMAT` | `-o` | Force output format |

**Examples:**

```sh
meshioplusplus reorder part.vtu reordered.vtu
meshioplusplus reorder part.vtu reordered.vtu --method hilbert
meshioplusplus reorder part.vtu reordered.vtu --method rcm --report
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

# CLI Reference

The `meshioplusplus` command-line tool is installed alongside the Python package.

```
meshioplusplus --version
meshioplusplus --help
meshioplusplus <subcommand> --help
```

The same verbs are also available as a **standalone native C++ binary** that needs no Python interpreter and no pybind11 extension — see [Native CLI (C API)](./c_api.md#native-command-line-binary) for how to build it (`build/configure.sh --cli --build`), or download a ready-to-run, statically-linked build for Linux/macOS/Windows from the [GitHub Releases](https://github.com/loumalouomega/meshioplusplus/releases) page. It shares the [C API](./c_api.md)'s flat-surface limitations: point/cell **sets** and `convert -s/-d` are unavailable there (they live only in the Python `Mesh`), and there is no Python fallback for formats whose C++ reader raises. Everything below otherwise applies identically to both.

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

## meshioplusplus diff

[Compare](./diff.md) two meshes and report whether they are equivalent within a
tolerance. The **exit code is nonzero when the meshes differ** and zero when they
are equal, so it drops straight into CI / shell scripts / Makefiles.

```
meshioplusplus diff [options] INFILE_A INFILE_B
```

| Option | Short | Description |
|--------|-------|-------------|
| `--atol ATOL` | | Absolute tolerance in `abs_err <= atol + rtol*|expected|` (default `1e-12`) |
| `--rtol RTOL` | | Relative tolerance (default `1e-9`) |
| `--unordered` | | Match points by spatial proximity (tolerant to a shuffled node order) |
| `--exact` | | Only a bitwise-identical result passes (tolerated drift exits nonzero) |
| `--quiet` | `-q` | Print nothing; communicate equality only via the exit code |
| `--input-format-a FORMAT` | | Force the format of the first file |
| `--input-format-b FORMAT` | | Force the format of the second file |

**Examples:**

```sh
meshioplusplus diff a.vtu b.vtu
meshioplusplus diff a.vtu b.vtu --atol 1e-8 --rtol 1e-6
meshioplusplus diff a.msh b.vtu --unordered
meshioplusplus diff expected.vtu actual.vtu --quiet || echo "regression!"
```

---

## meshioplusplus merge

Merge two or more mesh files into one.

```
meshioplusplus merge [options] FILE... OUTFILE
```

Takes two or more input meshes followed by the output file.

| Option | Short | Description |
|--------|-------|-------------|
| `--input-format FORMAT` | `-i` | Force input format (applied to every input) |
| `--output-format FORMAT` | `-o` | Force output format |
| `--weld` | | Merge coincident nodes within `--atol` |
| `--atol ATOL` | | Coincidence tolerance for `--weld` (default `1e-8`) |
| `--data-policy POLICY` | | `intersection` (default, keep only data keys present in every input) or `fill` (keep every key, filling missing rows with NaN) |
| `--drop-duplicate-cells` | | With `--weld`, drop cells that become identical after welding |
| `--no-source-tag` | | Do not add the per-cell `source_mesh_id` tag |
| `--quiet` | `-q` | Do not print the merge summary |

Prints a summary of points/cells in and out (and points welded, with `--weld`) unless `--quiet` is given.

**Examples:**

```sh
meshioplusplus merge a.vtu b.vtu merged.vtu
meshioplusplus merge a.vtu b.vtu c.vtu merged.vtu --weld --atol 1e-6
meshioplusplus merge a.vtu b.vtu merged.vtu --data-policy fill
```

---

## meshioplusplus transform

Apply an affine transform to a mesh's point coordinates (see [transform](/transform)).

```
meshioplusplus transform [options] INFILE OUTFILE
```

| Option | Description |
|--------|-------------|
| `--translate x,y,z` | Translation |
| `--scale sx,sy,sz` | Per-axis scale (or a single scalar for uniform) |
| `--rotate axis,deg` | Rotation; `axis` is `x`/`y`/`z` or `nx,ny,nz` (angle in degrees) |
| `--matrix m00,...,m33` | A row-major 4×4 affine matrix (16 values) |
| `--scale-units FACTOR` | Uniform unit-scale factor (e.g. `0.001`) |
| `--rotate-data` | Also rotate vector/tensor `point_data` by the transform |
| `--input-format` / `--output-format` (`-i`/`-o`) | Force input/output format |

Give exactly one transform source. Values starting with `-` need the `=` form
(`--translate=-1,0,0`).

**Examples:**

```sh
meshioplusplus transform in.vtu out.vtu --translate 1,2,3
meshioplusplus transform in.vtu out.vtu --rotate z,90
meshioplusplus transform in.vtu out.vtu --scale-units 0.001
```

---

## meshioplusplus clean

Weld / prune / de-dup a mesh in one pass (see [clean](/clean)).

```
meshioplusplus clean [options] INFILE OUTFILE
```

| Option | Description |
|--------|-------------|
| `--weld` | Fuse coincident points within `--atol` |
| `--atol ATOL` | Weld tolerance (default `1e-8`) |
| `--remove-orphans` | Drop unused points |
| `--drop-degenerate` | Drop degenerate cells |
| `--drop-duplicates` | Drop exact-duplicate cells |
| `--input-format` / `--output-format` (`-i`/`-o`) | Force input/output format |

With no step flags, the default set runs (remove-orphans + drop-degenerate +
drop-duplicates, **no** weld). A removal summary is printed.

**Examples:**

```sh
meshioplusplus clean in.vtu out.vtu
meshioplusplus clean in.vtu out.vtu --weld --atol 1e-6
```

---

## meshioplusplus crop

Extract the part of a mesh inside a bounding box or half-space (see [crop](/crop)).

```
meshioplusplus crop [options] INFILE OUTFILE
```

| Option | Description |
|--------|-------------|
| `--bbox xmin,ymin,zmin,xmax,ymax,zmax` | Axis-aligned bounding box |
| `--plane px,py,pz,nx,ny,nz` | Half-space (point + normal), keep `(p−point)·normal ≥ 0` |
| `--mode all\|any` | Keep a cell if ALL (default) or ANY node is inside |
| `--record-ids` | Attach original point/cell ids as data arrays |
| `--input-format` / `--output-format` (`-i`/`-o`) | Force input/output format |

Give exactly one of `--bbox`/`--plane`. Negative values need the `=` form
(`--bbox=-1,-1,-1,1,1,1`).

**Examples:**

```sh
meshioplusplus crop in.vtu out.vtu --bbox 0,0,0,1,1,1
meshioplusplus crop in.vtu out.vtu --plane 0.5,0,0,1,0,0 --mode any
```

---

## meshioplusplus split

Partition a mesh into several files by type, region, or connected component (see [split](/split)).

```
meshioplusplus split [options] INFILE OUTPATTERN
```

`OUTPATTERN` must contain `{key}`, replaced by each piece's key.

| Option | Description |
|--------|-------------|
| `--by type\|region\|component` | Split criterion (default `type`) |
| `--tag NAME` | For `--by region`: the integer `cell_data` name to split on |
| `--input-format` / `--output-format` (`-i`/`-o`) | Force input/output format |

Prints how many pieces were produced and their sizes.

**Examples:**

```sh
meshioplusplus split in.vtu 'out_{key}.vtu' --by type
meshioplusplus split in.vtu 'out_{key}.vtu' --by component
```

---

## meshioplusplus stats

Print geometric statistics of a mesh (see [stats](/stats)).

```
meshioplusplus stats [options] INFILE
```

| Option | Description |
|--------|-------------|
| `--json` | Emit the statistics as JSON |
| `--input-format` (`-i`) | Force input format |

Prints the bounding box, extent, centroid, per-cell-type counts, total area,
signed/unsigned volume, and inverted-cell count. This complements `info` (which
is topological) with geometric measures.

**Examples:**

```sh
meshioplusplus stats mesh.vtu
meshioplusplus stats mesh.vtu --json
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

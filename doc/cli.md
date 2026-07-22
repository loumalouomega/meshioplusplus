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

**Data-driven colouring** (SVG/TikZ output only):

| Option | Description |
|--------|-------------|
| `--color-by NAME` | `point_data` or `cell_data` array to colour the faces by |
| `--component I` | Component of a multi-component array (default: its magnitude) |
| `--cmap NAME` | `viridis` (default), `coolwarm` or `turbo` |
| `--vmin V` / `--vmax V` | Colour range (default: the drawn faces' finite range) |
| `--nan-color C` | Colour for NaN/infinite values (default: `#808080` / `gray`) |
| `--colorbar` | Append a gradient bar with min/max labels |

Point data colours a face by the mean of its corner values, cell data by its
owning cell's value — for a volume mesh, found through the skin's
`surface:parent_cell` provenance. `--color-by` with any other output format is
an error, as is any of the modifier flags without `--color-by`. See the
[SVG](./formats/svg.md#data-driven-colouring) and [TikZ](./formats/tikz.md)
format pages for the full semantics.

**Examples:**

```sh
meshioplusplus convert mesh.msh mesh.vtu
meshioplusplus convert -i gmsh -o vtk mesh.msh mesh.vtk
meshioplusplus convert --ascii mesh.msh mesh.vtu
meshioplusplus convert --sets-to-int-data mesh.inp mesh.xdmf
meshioplusplus convert mesh.msh skin.stl   # volume mesh -> boundary-skin STL

# colour a vector figure by a field
meshioplusplus convert mesh.vtu figure.svg --color-by temperature --colorbar
meshioplusplus convert mesh.vtu figure.tikz --color-by damage --cmap coolwarm \
    --vmin 0 --vmax 1
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

## meshioplusplus slice

Compute the planar cross-section of a mesh — the intersection with a plane, one
dimension below the cut cells (a volume mesh → a `triangle`/`quad` surface, a 2D
surface → a `line` mesh). See [slice](/slice).

```
meshioplusplus slice [options] INFILE OUTFILE
```

| Option | Description |
|--------|-------------|
| `--origin x,y,z` | A point on the cutting plane (default `0,0,0`) |
| `--normal x,y,z` | The plane normal (default `0,0,1`; non-zero) |
| `--record-parent-ids` | Attach `slice:parent_cell` (the input cell each section cell was cut from) |
| `--input-format` / `--output-format` (`-i`/`-o`) | Force input/output format |

Negative components need the `=` form (`--normal=0,0,-1`). Unlike `crop`, which
keeps whole cells on one side, `slice` computes the intersection itself.

**Examples:**

```sh
meshioplusplus slice in.vtu section.vtu --origin 0,0,0.5 --normal 0,0,1
meshioplusplus slice part.msh section.vtu --normal=0,0,-1 --record-parent-ids
```

---

## meshioplusplus isosurface

Compute the level set(s) of a scalar `point_data` field — the data-driven
sibling of `slice`, and like it one dimension below the cut cells (a volume mesh
→ a `triangle`/`quad` surface, a 2D surface → a `line` contour). See
[isosurface](/isosurface).

```
meshioplusplus isosurface [options] INFILE OUTFILE
```

| Option | Description |
|--------|-------------|
| `--array NAME` | The `point_data` array to contour (**required**) |
| `--values v1,v2,…` | The isovalues (**required**); sorted ascending, duplicates dropped |
| `--component I` | Component of a multi-component array; the row magnitude by default |
| `--record-parent-ids` | Attach `iso:parent_cell` (the input cell each contour cell was cut from) |
| `--input-format` / `--output-format` (`-i`/`-o`) | Force input/output format |

Negative isovalues need the `=` form (`--values=-1.5`). A `cell_data` name is
rejected: cell data is piecewise constant and has no level set — convert it with
`meshioplusplus data to-point` first. Every contour cell is tagged with
`iso:value` (Float64) and `iso:index` (Int64, the ordinal — the integer tag
`split --by region --tag …` needs).

**Examples:**

```sh
meshioplusplus isosurface part.vtu shell.vtu --array T --values 350
meshioplusplus isosurface part.vtu shells.vtu --array T --values 300,350,400
meshioplusplus isosurface part.vtu shell.vtu --array v --values=-1.5 --component 2

# one file per contour, via the integer ordinal tag
meshioplusplus split shells.vtu 'contour_{key}.vtu' --by region --tag iso:index
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

## meshioplusplus convert-cells

Convert a mesh's element representation (see [convert_cells](/convert_cells)).
Distinct from `convert`, which changes the *file format*.

```
meshioplusplus convert-cells [options] INFILE OUTFILE
```

| Option | Description |
|--------|-------------|
| `--mode linearize\|simplexify\|elevate` | Conversion to perform (default `linearize`) |
| `--record-parent-ids` | Attach `convert:parent_cell` cell_data of the source cell indices |
| `--input-format` / `--output-format` (`-i`/`-o`) | Force input/output format |

`linearize` drops higher-order nodes (`tetra10` → `tetra`) and prunes the points
that become unreferenced; `simplexify` decomposes cells into simplices of the
same dimension (`hexahedron` → 6 `tetra`); `elevate` promotes linear cells to
serendipity quadratic (`triangle` → `triangle6`), adding a node per unique edge.
A polyhedron block under `simplexify`, and `quad9`/`hexahedron27` under
`elevate`, are errors.

**Examples:**

```sh
meshioplusplus convert-cells in.msh out.vtu --mode linearize
meshioplusplus convert-cells in.msh out.vtu --mode simplexify --record-parent-ids
meshioplusplus convert-cells in.msh out.vtu --mode elevate
```

---

## meshioplusplus refine

Uniformly refine a mesh, subdividing every cell into same-type children (see
[refine](/refine)).

```
meshioplusplus refine [options] INFILE OUTFILE
```

| Option | Description |
|--------|-------------|
| `--levels N` | How many times to subdivide (default `1`) |
| `--record-parent-ids` | Attach `refine:parent_cell` cell_data of the original cell indices |
| `--input-format` / `--output-format` (`-i`/`-o`) | Force input/output format |

One level splits a `triangle`/`quad` into 4 and a `tetra`/`wedge`/`hexahedron`
into 8, inserting nodes at edge, quad-face and body midpoints. Those nodes are
shared between neighbouring cells, so the result has no hanging nodes, and
`point_data` is interpolated onto them. Higher-order cells, `pyramid`, and ragged
polygon/polyhedron blocks have no same-type subdivision and are errors —
`convert-cells --mode linearize` (or `--mode simplexify`) first.

Note the cell count grows as `4^levels` (2D) or `8^levels` (3D), so `--levels 3`
is already a 512× increase on a volume mesh.

**Examples:**

```sh
meshioplusplus refine in.msh out.vtu
meshioplusplus refine in.msh out.vtu --levels 2
meshioplusplus refine coarse.vtu fine.vtu --levels 2 --record-parent-ids
```

---

## meshioplusplus smooth

Relax point coordinates toward their edge-neighbour centroids to improve element
shape (see [smoothing](/smooth)).

```
meshioplusplus smooth [options] INFILE OUTFILE
```

| Option | Description |
|--------|-------------|
| `--method taubin\|laplacian` | Smoothing operator (default `taubin`, which does not shrink; `laplacian` is stronger per pass but shrinks) |
| `--iterations N` | How many iterations to run; for `taubin` one iteration is two passes (default `10`) |
| `--lambda L` | Relaxation factor of the smoothing pass, in (0, 1); the default depends on `--method` (`0.5` laplacian, `0.33` taubin) |
| `--mu=M` | `taubin` only: un-shrinking factor, must satisfy `mu < -lambda < 0` (default `-0.34`; the negative value needs the `--mu=` form) |
| `--no-fix-boundary` | Let boundary nodes move too (they are pinned by default) |
| `--no-preserve-features` | Do not pin feature nodes (sharp corners/creases are kept by default) |
| `--feature-angle A` | Angle in degrees between boundary facet normals above which their shared nodes are pinned as features (default `30`) |
| `--no-guard-inversion` | Do not reject moves that would invert an incident cell |
| `--quiet` (`-q`) | Suppress the summary output |
| `--input-format` / `--output-format` (`-i`/`-o`) | Force input/output format |

Only the point coordinates move: connectivity, `cell_data`, `field_data` and
`point_data` values come through unchanged, the point and cell counts are
unchanged, and the points array keeps its input dtype. Neighbours are the nodes
joined by an actual cell *edge*, so a structured hex block is a fixed point.
Nodes of blocks whose edge topology is unknown — the higher-order family, the
VTK-Lagrange types, `custom` — are pinned rather than guessed at. Unless
suppressed, the summary reports the number of nodes moved, the largest net
displacement, and how many moves the inversion guard rejected.

**Examples:**

```sh
meshioplusplus smooth noisy.vtu smooth.vtu
meshioplusplus smooth noisy.vtu smooth.vtu --iterations 40
meshioplusplus smooth in.msh out.vtu --method laplacian --lambda 0.4
meshioplusplus smooth in.msh out.vtu --mu=-0.4 --feature-angle 45
meshioplusplus smooth in.msh out.vtu --no-fix-boundary --no-guard-inversion -q
```

---

## meshioplusplus interpolate

Sample data arrays from a SOURCE mesh onto a TARGET mesh (see
[interpolation](/interpolate)). The output is a copy of the target — geometry,
connectivity and its own data preserved exactly — with the requested source
arrays sampled onto it.

```
meshioplusplus interpolate [options] SOURCE TARGET OUTFILE
```

| Option | Description |
|--------|-------------|
| `--method nearest\|barycentric` | Nearest source point (default; dtype-preserving) or linear interpolation in a simplexified source (exact on a linear field; Float64) |
| `--arrays a,b` | Comma-separated source array names to transfer (default: every source `point_data` array; `cell_data` transfers only when named) |
| `--extrapolate` | `barycentric` only: give a target point outside the source domain the nearest source point's value instead of the default value |
| `--default-value=V` | `barycentric` only: fill value for target points outside the source domain (default `0`; negative values need the `--default-value=` form) |
| `--on-conflict error\|overwrite\|suffix` | What to do when a transferred name already exists on the target (default `error`; `suffix` writes to `NAME_interp`) |
| `--quiet` (`-q`) | Suppress the transfer summary |
| `--input-format` / `--output-format` (`-i`/`-o`) | Force input (both files) / output format |

Source `point_data` is sampled at the target's points, source `cell_data` at
the target's cell centroids — always by nearest source-cell centroid, whatever
the method. Under `barycentric` the source is simplexified first, so on a
quad/hex source the result is the simplex-linear interpolant, and triangle
sources are evaluated in the xy-plane (use `nearest` for a curved surface
embedded in 3D).

**Examples:**

```sh
meshioplusplus interpolate coarse.vtu fine.vtu out.vtu
meshioplusplus interpolate coarse.vtu fine.vtu out.vtu --method barycentric
meshioplusplus interpolate a.msh b.msh out.vtu --arrays T,v --on-conflict suffix
meshioplusplus interpolate a.msh b.msh out.vtu --method barycentric --extrapolate
```

---

## meshioplusplus partition

Decompose a mesh into N balanced parts for domain decomposition (see
[partitioning](/partition)) — the count-driven complement to `split`.

```
meshioplusplus partition [options] INFILE OUTPATTERN
```

`OUTPATTERN` must contain `{part}` (e.g. `out_{part}.vtu`), expanded once per
piece — except with `--labels-only`, where it is a single plain path.

| Option | Description |
|--------|-------------|
| `--nparts N` (`-n`) | Number of parts (required, `>= 1`) |
| `--method M` | `sfc` (Hilbert curve cut, always available), `kahip` (KaHIP `kaffpa`; needs a KaHIP-enabled build, fails by name otherwise), or `auto` (default: kahip when available, else sfc) |
| `--imbalance F` | KaHIP only: allowed imbalance fraction in (0, 1) (default `0.03`) |
| `--mode M` | KaHIP only: `fast` / `eco` / `strong` (default `eco`) |
| `--seed N` | KaHIP only: random seed (deterministic per seed; default `0`) |
| `--weights NAME` | Scalar `cell_data` array of per-cell weights to balance instead of the cell count |
| `--record-ids` | Attach `partition:original_point_id` / `partition:original_cell_id` arrays to every piece |
| `--labels-only` | Write the input once with the Int64 `partition:part` cell_data attached instead of writing pieces |
| `--ghost-layers N` | Reserved; only `0` is supported |
| `--input-format` / `--output-format` (`-i`/`-o`) | Force input/output format |

Every piece keeps the input's cell-block structure 1:1 (empty blocks included,
unlike `split`), so the pieces recombine into the input — each cell lands in
exactly one piece.

**Examples:**

```sh
meshioplusplus partition domain.msh 'domain_{part}.vtu' --nparts 4
meshioplusplus partition domain.msh 'domain_{part}.vtu' -n 16 --method kahip --mode strong
meshioplusplus partition domain.msh labelled.vtu --nparts 4 --labels-only
meshioplusplus partition domain.msh 'p_{part}.vtu' -n 8 --weights cost --record-ids
```

---

## meshioplusplus data

A nested group of nine verbs operating on a mesh's `point_data` / `cell_data` /
`field_data` arrays (see [data operations](/data_operations)). **The geometry is
never modified** by any of them — points, connectivity, block order and block
types come through bit-identical.

```
meshioplusplus data <subcommand> [options]
```

| Subcommand | Description |
|------------|-------------|
| `info` | Summarize every data array (see [data summary](/data_info)) |
| `rename` | Rename data arrays (see [array management](/data_manage)) |
| `drop` | Drop data arrays by name |
| `keep` | Keep only the named data arrays |
| `to-cell` | Average `point_data` onto the cells (see [averaging](/data_average)) |
| `to-point` | Average `cell_data` onto the points |
| `calc` | Derive an array from an expression (see [expressions](/data_calc)) |
| `clamp` | Clamp values into a range (see [conditioning](/data_condition)) |
| `normalize` | Rescale values to a target range |

Every verb takes `--input-format` (`-i`), and every verb but `info` takes an
`OUTFILE` and `--output-format` (`-o`).

::: warning Colons in names
Data names routinely contain colons (`gmsh:physical`). `data rename` therefore
splits its `OLD:NEW` value on the **last** colon — `--point gmsh:physical:tag`
renames `gmsh:physical` to `tag`. `data calc` splits `NAME = EXPR` on the
**first** `=`. `drop`/`keep` take a comma-separated name list with no prefix, so
colons there are unambiguous. The Python CLI and the native binary implement
identical rules.
:::

### data info

| Option | Description |
|--------|-------------|
| `--json` | Emit the summary as JSON |

Prints location, name, dtype, component count, entry count, min/max/mean and
NaN/inf counts for every array. Read-only.

### data rename / drop / keep

| Option | Description |
|--------|-------------|
| `--point`, `--cell`, `--field` | `OLD:NEW` for `rename` (repeatable); a comma-separated name list for `drop`/`keep` |
| `--ignore-missing` | Skip names that do not exist instead of failing (`drop`/`keep`) |

For `keep`, a location that is not named at all is left untouched; naming it
with an empty list drops everything there.

### data to-cell / to-point

| Option | Description |
|--------|-------------|
| `--keys` | Comma-separated names to convert (default: all at the source location) |
| `--target-suffix` | Append this to each output name (default: keep the same name) |
| `--weighted` | `to-point` only: weight by cell measure (area/volume) instead of counting cells equally |

The output is always `float64` — the mean of an integer field is not an integer.

### data calc

| Option | Description |
|--------|-------------|
| `--point`, `--cell`, `--field` | `NAME = EXPRESSION` (repeatable) |
| `--overwrite` | Allow replacing an array that already exists |

The expression grammar accepts `+ - * /`, unary minus, parentheses, numeric
literals, array names, and `abs`/`sqrt`/`min`/`max`/`norm` — nothing else is
evaluated.

### data clamp / normalize

| Option | Description |
|--------|-------------|
| `--point`, `--cell`, `--field` | Comma-separated names (default: all at that location) |
| `--min`, `--max` | `clamp` only: the bounds (both required) |
| `--to LO,HI` | `normalize` only: target range (default `0,1`) |
| `--zero-mean` | `normalize` only: standardize to zero mean / unit std instead |
| `--magnitude` | Condition by row magnitude instead of per component |
| `--nan` | `ignore` (default), `replace` or `fail` |
| `--nan-value` | Replacement used with `--nan replace` |
| `--suffix` | Store as `NAME+SUFFIX` instead of replacing in place |

**Examples:**

```sh
meshioplusplus data info mesh.vtu
meshioplusplus data info mesh.vtu --json

meshioplusplus data rename in.vtu out.vtu --point T:temperature
meshioplusplus data drop   in.vtu out.vtu --point a,b --cell c
meshioplusplus data keep   in.vtu out.vtu --point T,p --cell mat

meshioplusplus data to-cell  in.vtu out.vtu --keys T,p --target-suffix _c
meshioplusplus data to-point in.vtu out.vtu --keys stress --weighted

meshioplusplus data calc in.vtu out.vtu --point "speed = norm(velocity)"
meshioplusplus data calc in.vtu out.vtu --cell  "dp = p_new - p_old"

meshioplusplus data clamp     in.vtu out.vtu --point T --min 0 --max 100
meshioplusplus data normalize in.vtu out.vtu --cell damage --to 0,1
meshioplusplus data normalize in.vtu out.vtu --point T --zero-mean
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

## Selective reads and fast summaries

`info --fast` summarizes a file from its header instead of loading it, and `convert` can
narrow what it reads:

```bash
meshioplusplus info --fast big.vtu
meshioplusplus convert --points-only in.vtu out.vtu     # geometry, no data arrays
meshioplusplus convert --arrays u,p in.vtu out.vtu      # only these data arrays
```

`--points-only` keeps connectivity — it narrows data, not topology. `arrays` with an empty
list keeps no arrays; omitting the flag keeps every array.

Formats without a header-only path are read in full and `info --fast` says so explicitly
(`no header-only path for this format; the file was read in full`) rather than implying a
saving that did not happen. See [Selective reads](selective_read.md).

`--points-only`/`--arrays` are rejected alongside `-s`/`-d`, which convert exactly the data
arrays that were skipped.

## Compression codecs

```bash
meshioplusplus compress --codec lz4 mesh.vtu
```

`--codec zlib|lz4|zstd` selects the VTK XML block codec for `.vtu`/`.vtp`. zlib is the
default; `lz4` stays ParaView-readable, `zstd` is a meshio++ extension that ParaView cannot
read. The flag is **rejected** for formats with no block codec rather than silently ignored.
See [Compression codecs](codecs.md).

Both CLIs — the Python one and the native `meshioplusplus` binary — accept these identically.

## `view` / `screenshot` — the native viewer

```sh
meshioplusplus view part.msh
meshioplusplus view part.msh --kind surface --color-by material
meshioplusplus screenshot part.msh out.png --size 1600x1200
```

Options: `--input-format/-i`, `--kind {auto,surface,volume,curve,points}`,
`--color-by NAME`, `--name NAME`; `screenshot` adds `--size WIDTHxHEIGHT` and
`--transparent`.

These mirror the Python CLI's verbs, but in the **native binary** they are only
functional in a build configured with
[Polyscope](https://polyscope.run):

```sh
git submodule update --init --recursive     # Polyscope vendors its own submodules
build/configure.sh --cli --with-polyscope --build
```

They are listed in `--help` in every build; without the flag they report it
rather than silently not existing. The **prebuilt release binaries do not
include the viewer** — they are deliberately dependency-free single files, and
Polyscope needs OpenGL, GLFW and X11. Use the Python CLI (`pip install
meshioplusplus[viewer]`) or the [browser viewer](/viewer) if you would rather
not build from source.

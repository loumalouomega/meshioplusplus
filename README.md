<p align="center">
  <a href="https://github.com/loumalouomega/meshioplusplus"><img alt="meshio++" src="https://raw.githubusercontent.com/loumalouomega/meshioplusplus/main/doc/logo/logo-with-text.svg" width="60%"></a>
  <p align="center">I/O for mesh files.</p>
</p>


[![PyPi Version](https://img.shields.io/pypi/v/meshioplusplus.svg?style=flat-square)](https://pypi.org/project/meshioplusplus/) [![npm Version](https://img.shields.io/npm/v/%40meshioplusplus%2Fwasm.svg?style=flat-square)](https://www.npmjs.com/package/@meshioplusplus/wasm) [![PyPI pyversions](https://img.shields.io/pypi/pyversions/meshioplusplus.svg?style=flat-square)](https://pypi.org/project/meshioplusplus/) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21384760.svg?style=flat-square)](https://doi.org/10.5281/zenodo.21384760)

[![C++][c++-image]][c++standard] [![Python](https://img.shields.io/badge/Python-3.9%2B-3776ab.svg?style=flat-square&logo=python&logoColor=white)](https://pypi.org/project/meshioplusplus/) [![C](https://img.shields.io/badge/C-99-a8b9cc.svg?style=flat-square&logo=c&logoColor=white)](doc/c_api.md) [![Fortran](https://img.shields.io/badge/Fortran-2008-734f96.svg?style=flat-square&logo=fortran&logoColor=white)](doc/fortran.md) [![WebAssembly](https://img.shields.io/badge/WebAssembly-npm-654ff0.svg?style=flat-square&logo=webassembly&logoColor=white)](https://www.npmjs.com/package/@meshioplusplus/wasm) [![TypeScript](https://img.shields.io/badge/TypeScript-viewer-3178c6.svg?style=flat-square&logo=typescript&logoColor=white)](src/viewer/)

[![GitHub stars](https://img.shields.io/github/stars/loumalouomega/meshioplusplus.svg?style=flat-square&logo=github&label=Stars&logoColor=white)](https://github.com/loumalouomega/meshioplusplus) [![PyPi downloads](https://img.shields.io/pypi/dm/meshioplusplus.svg?style=flat-square)](https://pypistats.org/packages/meshioplusplus)

[![gh-actions](https://img.shields.io/github/actions/workflow/status/loumalouomega/meshioplusplus/ci.yml?branch=master&style=flat-square)](https://github.com/loumalouomega/meshioplusplus/actions?query=workflow%3Aci) [![codecov](https://img.shields.io/codecov/c/github/loumalouomega/meshioplusplus.svg?style=flat-square)](https://app.codecov.io/gh/loumalouomega/meshioplusplus) [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg?style=flat-square)](https://github.com/psf/black)

[c++-image]: https://img.shields.io/badge/C++-20-blue.svg?style=flat&logo=c%2B%2B
[c++standard]: https://isocpp.org/std/the-standard

There are various mesh formats available for representing unstructured meshes. meshio++ can read and write all of the following and smoothly converts between them:

> [Abaqus](https://help.3ds.com/2024/english/dssimulia_established/SIMACAEMODRefMap/simamod-c-inputsyntax.htm) (`.inp`),
> ANSYS msh (`.msh`),
> [Ansys/APDL coded database](https://www.ansys.com) (`.cdb`, `.inp`),
> [AVS-UCD](https://lanl.github.io/LaGriT/pages/docs/read_avs.html) (`.avs`),
> [CGNS](https://cgns.github.io/) (`.cgns`),
> [DOLFIN XML](https://manpages.ubuntu.com/manpages/jammy/en/man1/dolfin-convert.1.html) (`.xml`),
> [COMSOL](https://www.comsol.com) (`.mphtxt`),
> [Exodus](https://nschloe.github.io/meshio/exodus.pdf) (`.e`, `.exo`),
> [EnSight Gold](https://vis.lbl.gov/archive/NERSC/Software/ensight/doc/OnlineHelp/UM-C11.pdf) (geometry, `.case`/`.geo`),
> [FLAC3D](https://www.itascacg.com/software/flac3d) (`.f3grid`),
> [FLUX](https://www.altair.com/flux/) (mesh `.pf3`, field `.dex`),
> [FreeFem++](https://freefem.org/) (`.msh`),
> [H5M](https://www.mcs.anl.gov/~fathom/moab-docs/h5mmain.html) (`.h5m`),
> [HMF](https://loumalouomega.github.io/meshioplusplus/formats/hmf) (`.hmf`, experimental, meshio++-specific),
> [I-deas Universal / UNV](https://www.ceas3.uc.edu/sdrluff/) (`.unv`),
> [ANSYS Fluent interpolation](https://github.com/victorsndvg/FEconv) (`.ip`),
> [Kratos/MDPA](https://github.com/KratosMultiphysics/Kratos/wiki/Input-data) (`.mdpa`),
> [Medit](https://people.sc.fsu.edu/~jburkardt/data/medit/medit.html) (`.mesh`, `.meshb`),
> [MED/Salome](https://docs.salome-platform.org/latest/dev/MEDCoupling/developer/med-file.html) (`.med`),
> [Modulef](https://github.com/victorsndvg/FEconv) (mesh `.mfm`, field `.mff`),
> [Nastran](https://help.autodesk.com/view/NSTRN/2019/ENU/?guid=GUID-42B54ACB-FBE3-47CA-B8FE-475E7AD91A00) (bulk data, `.bdf`, `.fem`, `.nas`),
> [Netgen](https://github.com/ngsolve/netgen) (`.vol`, `.vol.gz`),
> [Neuroglancer precomputed format](https://github.com/google/neuroglancer/tree/master/src/datasource/precomputed#mesh-representation-of-segmented-object-surfaces),
> [Gmsh](https://gmsh.info/doc/texinfo/gmsh.html#File-formats) (format versions 2.2, 4.0, and 4.1, `.msh`),
> [OBJ](https://en.wikipedia.org/wiki/Wavefront_.obj_file) (`.obj`),
> [OFF](https://segeval.cs.princeton.edu/public/off_format.html) (`.off`),
> [OpenFOAM polyMesh](https://www.openfoam.com/) (`.foam`, read-only),
> [PERMAS](https://www.intes.de) (`.post`, `.post.gz`, `.dato`, `.dato.gz`),
> [PLY](<https://en.wikipedia.org/wiki/PLY_(file_format)>) (`.ply`),
> [STL](<https://en.wikipedia.org/wiki/STL_(file_format)>) (`.stl`),
> [Tecplot .dat](http://paulbourke.net/dataformats/tp/),
> [TetGen .node/.ele](https://wias-berlin.de/software/tetgen/fformats.html),
> [Triangle .node/.ele/.poly](https://www.cs.cmu.edu/~quake/triangle.html),
> [SVG](https://www.w3.org/TR/SVG/) (output only; 2D direct, 3D via skin projection) (`.svg`),
> [TikZ](https://tikz.dev/) (LaTeX output only; 2D direct, 3D via skin projection) (`.tikz`),
> [SU2](https://su2code.github.io/docs_v7/Mesh-File/) (`.su2`),
> [UGRID](https://www.simcenter.msstate.edu/software/documentation/ug_io/3d_grid_file_type_ugrid.html) (`.ugrid`),
> [VTK](https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf) (`.vtk`),
> [VTP](https://docs.vtk.org/en/latest/vtk_file_formats/vtkxml_file_format.html) (`.vtp`),
> [VTU](https://vtk.org/Wiki/VTK_XML_Formats) (`.vtu`),
> [WKT](https://en.wikipedia.org/wiki/Well-known_text_representation_of_geometry) ([TIN](https://en.wikipedia.org/wiki/Triangulated_irregular_network)) (`.wkt`),
> [XDMF](https://xdmf.org/index.php/XDMF_Model_and_Format) (`.xdmf`, `.xmf`).

<p align="center">
  <img alt="" src="https://raw.githubusercontent.com/loumalouomega/meshioplusplus/main/doc/logo/logo-icon-square.png" width="64">
</p>

meshio++ ships a **C++20 core** (built with pybind11 + scikit-build-core) that reads and writes most formats with zero-copy numpy at the I/O boundary, plus optional HDF5/netCDF acceleration and a **selectable parallel backend** (`AUTO` by default — prefers OpenMP, then STL+TBB, then sequential; override with `-DMESHIOPLUSPLUS_PARALLEL_BACKEND=...`). Every format has a pure-Python fallback, so behaviour and file compatibility are identical whether or not the native libraries are present. For a standalone C++ build use `build/configure.sh` (Linux/macOS) or `build/configure.bat` (Windows). Full docs (install, data model, per-format options, CLI) live at [the documentation site](https://loumalouomega.github.io/meshioplusplus/) (sources under [`doc/`](https://github.com/loumalouomega/meshioplusplus/tree/main/doc)).

Install with

```
pip install meshioplusplus[all]
```

(`[all]` pulls in all optional dependencies. By default, meshio++ only uses numpy.) You can then use the command-line tool

<!--pytest-codeblocks:skip-->

```sh
meshioplusplus convert    input.msh output.vtk   # convert between two formats

meshioplusplus info       input.xdmf             # show some info about the mesh

meshioplusplus compress   input.vtu              # compress the mesh file
meshioplusplus decompress input.vtu              # decompress the mesh file

meshioplusplus binary     input.msh              # convert to binary format
meshioplusplus ascii      input.msh              # convert to ASCII format

meshioplusplus merge      a.vtu b.vtu out.vtu    # merge meshes (optional --weld)

meshioplusplus transform  in.vtu out.vtu --translate 1,2,3   # affine transform
meshioplusplus clean      in.vtu out.vtu --weld              # weld / prune / de-dup
meshioplusplus crop       in.vtu out.vtu --bbox 0,0,0,1,1,1  # subset by region
meshioplusplus split      in.vtu 'out_{key}.vtu' --by type   # split by criterion
meshioplusplus stats      mesh.vtu                           # geometric statistics
meshioplusplus convert-cells in.msh out.vtu --mode simplexify  # hexes -> tetra
meshioplusplus refine     in.vtu out.vtu --levels 2          # uniform subdivision
meshioplusplus partition  in.vtu 'out_{part}.vtu' --nparts 4 # N balanced parts
meshioplusplus smooth     in.vtu out.vtu --iterations 20     # relax node positions
meshioplusplus interpolate src.vtu tgt.vtu out.vtu           # transfer fields across meshes
meshioplusplus slice      in.vtu out.vtu --normal 0,0,1      # planar cross-section

meshioplusplus data info  mesh.vtu                           # summarize data arrays
meshioplusplus data calc  in.vtu out.vtu --point "s = norm(v)"   # derive a field
meshioplusplus data to-cell  in.vtu out.vtu --keys T         # point -> cell average
meshioplusplus data normalize in.vtu out.vtu --cell damage --to 0,1
```

with any of the supported formats.

The same verbs are available as a **standalone C++ binary** that needs no Python: grab a ready-to-run, statically-linked build for Linux/macOS/Windows from the [GitHub Releases](https://github.com/loumalouomega/meshioplusplus/releases) page, or build it yourself with `build/configure.sh --cli --build` (or `-DMESHIOPLUSPLUS_BUILD_CLI=ON`). It links only the C++ core, so point/cell *sets* (and `convert -s/-d`) — which live only in the Python `Mesh` — are unavailable there; use the Python CLI for those.

In Python, simply do

<!--pytest-codeblocks:skip-->

```python
import meshioplusplus

mesh = meshioplusplus.read(
    filename,  # string, os.PathLike, or a buffer/open file
    # file_format="stl",  # optional if filename is a path; inferred from extension
    # see meshioplusplus convert --help for all possible formats
)
# mesh.points, mesh.cells, mesh.cells_dict, ...

# mesh.vtk.read() is also possible
```

to read a mesh. To write, do

```python
import meshioplusplus

# two triangles and one quad
points = [
    [0.0, 0.0],
    [1.0, 0.0],
    [0.0, 1.0],
    [1.0, 1.0],
    [2.0, 0.0],
    [2.0, 1.0],
]
cells = [
    ("triangle", [[0, 1, 2], [1, 3, 2]]),
    ("quad", [[1, 4, 5, 3]]),
]

mesh = meshioplusplus.Mesh(
    points,
    cells,
    # Optionally provide extra data on points, cells, etc.
    point_data={"T": [0.3, -1.2, 0.5, 0.7, 0.0, -3.0]},
    # Each item in cell data must match the cells array
    cell_data={"a": [[0.1, 0.2], [0.4]]},
)
mesh.write(
    "foo.vtk",  # str, os.PathLike, or buffer/open file
    # file_format="vtk",  # optional if first argument is a path; inferred from extension
)

# Alternative with the same options
meshioplusplus.write_points_cells("foo.vtk", points, cells)
```

For both input and output, you can optionally specify the exact `file_format` (in case you would like to enforce ASCII over binary VTK, for example).

#### Skin extraction

`meshioplusplus.extract_skin` derives the boundary surface of a 3D volume mesh (the [Kratos `SkinDetectionProcess`](https://github.com/KratosMultiphysics/Kratos) face-hashing algorithm — faces occurring exactly once are boundary; points are compacted, `point_data` follows):

<!--pytest-codeblocks:skip-->

```python
vol = meshioplusplus.read("part.msh")     # tetra/hexa/wedge/pyramid mesh
skin = meshioplusplus.extract_skin(vol)   # triangle/quad/... surface mesh
```

The **STL and PLY writers do this automatically** for volume meshes (pass `skin=False` for the legacy drop-volume-cells behavior), and the **SVG/TikZ writers render 3D meshes** by projecting the skin through an orthographic camera (`azimuth`/`elevation`/`roll` in degrees, default the classic CAD isometric view) with painter's-algorithm depth ordering — that is exactly how the Stanford-bunny logo above is drawn.

#### Publication-quality vector figures

The SVG and TikZ writers can colour each face by a data array, turning them into figures you can drop straight into a paper — resolution-independent, and with **no extra dependency**: the colormaps are built into the core.

<!--pytest-codeblocks:skip-->

```python
annotated = meshioplusplus.attach_quality(mesh)
meshioplusplus.write(
    "quality.svg", annotated,
    color_by="quality:scaled_jacobian",   # or any point_data / cell_data name
    cmap="viridis",                       # viridis / coolwarm / turbo
    colorbar=True,
)
```

<img alt="a bracket coloured by scaled Jacobian" src="https://raw.githubusercontent.com/loumalouomega/meshioplusplus/main/doc/public/images/color_by_quality.svg" width="85%">

*The bundled bracket coloured by element quality — the same figure `tools/gen_doc_images.py` regenerates.*

**Point data** colours a face by the mean of its corner values, **cell data** by its owning cell's value — for a volume mesh, tracked through the extracted skin's parent-cell provenance, so a per-cell material or metric lands on the right facet. Multi-component arrays reduce to a `component` or to their magnitude; `vmin`/`vmax` set the range (default: the drawn faces' finite range), and non-finite values take `nan_color`. From the command line:

<!--pytest-codeblocks:skip-->

```sh
meshioplusplus convert mesh.vtu figure.svg --color-by temperature --colorbar
```

Colouring is available from Python, from C++ directly, and from both CLIs; the flat C/Fortran/WebAssembly bindings reach these writers through the shared registry and always emit the default styling.

#### Surface extraction

`meshioplusplus.extract_surface` is the general form of skin extraction: it picks the dimension automatically (a volume mesh → boundary faces, a 2D surface mesh → boundary edges) and can record each facet's parent cell id (`record_parent_ids=True`). See the [surface extraction docs](https://meshioplusplus.readthedocs.io) (`doc/extract_surface.md`).

<!--pytest-codeblocks:skip-->

```python
surf = meshioplusplus.extract_surface(vol)                  # faces (or edges for a 2D mesh)
edges = meshioplusplus.extract_surface(sheet, record_parent_ids=True)
```

#### Mesh quality

`meshioplusplus.compute_quality` scores every cell on a set of geometric quality metrics (area/volume, scaled Jacobian, aspect ratio, skewness, interior/dihedral angles, warpage) and flags inverted/degenerate cells; `attach_quality` writes them back as `cell_data`. See `doc/mesh_quality.md`.

<!--pytest-codeblocks:skip-->

```python
report = meshioplusplus.compute_quality(mesh)
print(report["num_inverted"], "inverted cells")
annotated = meshioplusplus.attach_quality(mesh)   # metrics as cell_data
```

#### Reordering / renumbering

`meshioplusplus.reorder` renumbers nodes and elements to reduce sparse-matrix bandwidth (Reverse Cuthill–McKee) or improve cache locality (Morton / Hilbert space-filling curves). It is a pure permutation — geometry and all data preserved — and returns the applied node/cell permutations so external arrays can be remapped. `compute_bandwidth` measures the before/after connectivity bandwidth. See `doc/reorder.md`.

<!--pytest-codeblocks:skip-->

```python
out = meshioplusplus.reorder(mesh, method="rcm")            # "morton" / "hilbert" too
out, node_perm, cell_perms = meshioplusplus.reorder(mesh, return_permutation=True)
print(meshioplusplus.compute_bandwidth(mesh), "->", meshioplusplus.compute_bandwidth(out))
```

#### Comparison (diff)

`meshioplusplus.diff` compares two meshes and reports whether they are equivalent within a tolerance (`abs_err <= atol + rtol*|expected|`), with a structured breakdown (points, cells, data, named sets) and an overall verdict (`identical` / `equal within tolerance` / `different`); `meshes_equal` is the boolean wrapper for test suites. An optional `unordered=True` mode matches nodes by spatial proximity, so a shuffled node order still compares equal. See `doc/diff.md`.

<!--pytest-codeblocks:skip-->

```python
assert meshioplusplus.meshes_equal(a, b, atol=1e-8)         # ideal in a regression test
report = meshioplusplus.diff(a, b, unordered=True)          # tolerant to shuffled node order
print(report["verdict"])
```

The `meshioplusplus diff a.vtu b.vtu` CLI verb sets a nonzero exit code when meshes differ, for direct use in CI / Makefiles.

#### Merge / combine

`meshioplusplus.merge` combines two or more meshes into one: it concatenates points (offsetting connectivity so indices stay valid), merges cell blocks by type, concatenates data (per a configurable `data_policy`), and tags each cell's origin. With `weld=True` it fuses coincident nodes across inputs within `atol` using a spatial hash (never O(N²)) — the standard way to stitch adjacent blocks into a watertight mesh. Overlapping set / field-data names are namespaced by source id. See `doc/merge.md`.

<!--pytest-codeblocks:skip-->

```python
combined = meshioplusplus.merge([a, b, c])                 # concatenate
welded = meshioplusplus.merge([left, right], weld=True, atol=1e-8)  # fuse the shared interface
```

#### Editing (transform / clean / crop / split) and statistics

A bundle of dependency-free mesh-editing utilities:

- **`meshioplusplus.transform`** — apply an affine transform (translate / scale / rotate / 4×4 matrix / unit-scale) to the points; connectivity and data are carried through. See `doc/transform.md`.
- **`meshioplusplus.clean`** — weld coincident points (spatial hash), drop degenerate and duplicate cells, and remove orphaned points, in one toggleable pass. See `doc/clean.md`.
- **`meshioplusplus.crop`** — extract the part of a mesh inside a bounding box or half-space, pruning unused points (`mode="all"`/`"any"`). See `doc/crop.md`.
- **`meshioplusplus.split`** — partition a mesh into several by cell type, connected component (flood-fill), or region (`cell_sets` / integer tag). See `doc/split.md`.
- **`meshioplusplus.compute_stats`** — geometric statistics (bounding box, centroid, per-type counts, area, signed/unsigned volume, inverted cells) — the geometric complement to `info`. See `doc/stats.md`.

<!--pytest-codeblocks:skip-->

```python
out = meshioplusplus.transform(mesh, rotate=("z", 90))
out = meshioplusplus.clean(mesh, weld=True, atol=1e-8)
sub = meshioplusplus.crop(mesh, bbox=[0, 0, 0, 1, 1, 1])
pieces = meshioplusplus.split(mesh, by="type")             # {"triangle": ..., ...}
s = meshioplusplus.compute_stats(mesh)                     # dict of measures
```

#### Cell conversion (linearize / simplexify / elevate)

**`meshioplusplus.convert_cells`** converts a mesh's *element representation* — which cell types it is built from — while leaving the object it describes intact. See `doc/convert_cells.md`.

- **`mode="linearize"`** — every higher-order cell becomes its linear base (`tetra10` → `tetra`, `hexahedron27` → `hexahedron`), keeping the corner connectivity verbatim and pruning the nodes that become unreferenced.
- **`mode="simplexify"`** — every cell is decomposed into simplices of the same topological dimension (`quad` → 2 `triangle`, `hexahedron` → 6 `tetra`, `wedge` → 3, `pyramid` → 2, an n-gon into an (n−2)-triangle fan). No points are added, each parent's `cell_data` is replicated to its children, and every emitted simplex is positively oriented with volume conserved.
- **`mode="elevate"`** — every linear cell is promoted to its serendipity quadratic counterpart (`triangle` → `triangle6`, `hexahedron` → `hexahedron20`), adding one node per unique edge at the edge midpoint with `point_data` set to the endpoint mean.

<!--pytest-codeblocks:skip-->

```python
linear = meshioplusplus.convert_cells(mesh, mode="linearize")
tets = meshioplusplus.convert_cells(mesh, mode="simplexify")   # hexes -> tetra
quadratic = meshioplusplus.convert_cells(mesh, mode="elevate")
```

Each mode is idempotent on cells it does not apply to, so it is safe on a mixed-order mesh, and output is byte-identical across mesh backends and thread counts.

#### Refinement

**`meshioplusplus.refine`** subdivides every cell into congruent children of the *same* cell type, increasing a mesh's resolution: `line` → 2, `triangle` → 4, `quad` → 4, `tetra` → 8, `wedge` → 8, `hexahedron` → 8, with `levels=n` applying the templates `n` times. See `doc/refine.md`.

New nodes sit at the midpoints of the parent's edges, quad faces and (hexahedron only) body, and carry the mean of that entity's corner values for every `point_data` array — so a linear field is interpolated exactly. Mid-edge and quad-face-centre nodes are **shared** between every cell touching the entity, so the refined mesh has no hanging nodes; each parent's `cell_data` row is replicated to its children.

<!--pytest-codeblocks:skip-->

```python
fine = meshioplusplus.refine(mesh)                      # one level
finer = meshioplusplus.refine(mesh, levels=2)           # 64x the cells in 3D
tagged = meshioplusplus.refine(mesh, record_parent_ids=True)
```

Children inherit the parent's orientation (zero newly-inverted cells for a well-oriented input), and volume is conserved — exactly for `tetra` always, and for `wedge`/`hexahedron` when the parent is affine. Higher-order cells, `pyramid`, and ragged blocks have no same-type subdivision and raise by name.

#### Partitioning

**`meshioplusplus.partition`** decomposes a mesh into exactly N balanced pieces for domain decomposition — the count-driven complement to the criterion-driven `split`. See `doc/partition.md`.

- **SFC** (the default fallback, always available, dependency-free): cells are cut into contiguous ranges along a Hilbert space-filling curve of their centroids — equal-weight part sizes differ by at most one cell, `weights=<cell_data>` balances a per-cell cost instead, and the assignment is deterministic and byte-identical across mesh backends and thread counts.
- **KaHIP** (the optional quality path): the shared-face dual graph goes through [KaHIP](https://github.com/KaHIP/KaHIP)'s serial `kaffpa()`, which actively minimizes the edge cut. Configure `imbalance` (default 3%), `mode` (`fast`/`eco`/`strong`, default `eco` — eco/strong carry the quality) and `seed`. **KaHIP is MIT-licensed like meshio++ itself, so enabling it changes nothing about licensing**; it is bring-your-own (`-DMESHIOPLUSPLUS_WITH_KAHIP=ON` + `KAHIP_ROOT`, Conan `with_kahip`, vcpkg feature `kahip` — next to the HDF5/zstd-style optional deps), links only the serial interface (no MPI), and `pip install meshioplusplus[kahip]` gives pure-Python installs the same quality path via the MIT `kahip` wheel. Requesting it where absent fails by name — never a silent downgrade.

<!--pytest-codeblocks:skip-->

```python
pieces = meshioplusplus.partition(mesh, 4)                    # list of 4 meshes
labels = meshioplusplus.partition_labels(mesh, 4)             # per-block Int64 part ids
quality = meshioplusplus.partition(mesh, 16, method="kahip", mode="strong")
```

Pieces keep the input's block structure 1:1, so they recombine into the input: every cell lands in exactly one piece (`ghost_layers` is reserved for halo growth and raises for now).

#### Smoothing

**`meshioplusplus.smooth`** relaxes point coordinates toward their edge-neighbour centroids to improve element shape, leaving topology and every data value alone: **only the points move**. See `doc/smooth.md`.

Both operators are driven by the same centroid displacement. **Laplacian** (`x <- x + lambda*L(x)`) smooths strongly per pass but shrinks — over 40 iterations on a jittered 8×8 quad grid it contracts the bounding box by 57%. **Taubin** (the default) follows each `+lambda` pass with a larger-magnitude `-mu` pass that deliberately un-shrinks, leaving the same grid 3.6% smaller. Neighbours are the nodes joined by an actual cell *edge*, not the element clique, so a structured hex block is a fixed point rather than being bevelled toward a sphere. Boundary nodes, feature nodes (incident boundary facet normals differing by more than `feature_angle`), an optional `frozen` mask, and the nodes of blocks whose edge topology is unknown are all pinned by default, and the inversion guard rejects any move that would turn a valid cell inverted.

<!--pytest-codeblocks:skip-->

```python
relaxed = meshioplusplus.smooth(mesh)                          # 10 Taubin iterations
harder = meshioplusplus.smooth(mesh, iterations=40)            # shrink-free even so
lap = meshioplusplus.smooth(mesh, method="laplacian", lambda_=0.4)  # note the underscore
out, report = meshioplusplus.smooth(mesh, return_report=True)  # nodes moved, max displacement
```

`lambda_` carries a trailing underscore because `lambda` is a Python keyword, and a negative value means "this method's own default" (0.5 Laplacian, 0.33 Taubin). Point and cell counts, connectivity, `cell_data`, `field_data`, `point_data` values and the points array's dtype all come through unchanged, and output is byte-identical across mesh backends and thread counts.

#### Field transfer (interpolation)

**`meshioplusplus.interpolate`** samples a **source** mesh's data arrays onto a **target** mesh — the first cross-mesh operation that transfers data (diff compares, merge concatenates). The result is a copy of the target, its own geometry/data/sets preserved exactly, with the requested source arrays sampled on: source `point_data` at the target's points, source `cell_data` by nearest source-cell centroid. `method="nearest"` (default) copies the nearest source point's value bit-for-bit; `method="barycentric"` simplexifies the source first and interpolates linearly — exact on a linear field — with `default_value`/`extrapolate` deciding what happens outside the source domain. Both search grids are bucket-grid spatial hashes (never O(N²)), and output is byte-identical across backends, thread counts and the C++/numpy boundary. See `doc/interpolate.md`.

<!--pytest-codeblocks:skip-->

```python
coarse = meshioplusplus.read("solution.vtu")   # carries point_data "T"
fine = meshioplusplus.read("remeshed.vtu")
mapped = meshioplusplus.interpolate(coarse, fine, method="barycentric")
meshioplusplus.write("mapped.vtu", mapped)
```

#### Slicing / cross-sections

**`meshioplusplus.slice`** computes the planar cross-section of a mesh — the actual intersection of the mesh with a plane, one topological dimension below the cut cells: a 3D volume mesh yields a `triangle`/`quad` surface, a 2D surface mesh a `line` mesh. Unlike `crop` (plane mode), which keeps whole cells on one side, `slice` computes the intersection and lowers the dimension. It uses robust marching tetrahedra (the input is simplexified first, so each cell's cross-section is a well-defined convex primitive), deduping crossing points on shared edges into single nodes so the section is watertight, and winding every face consistently toward the `+normal` side. Each section cell inherits its parent's `cell_data`; `record_parent_ids=True` attaches `slice:parent_cell`, and `point_data` is interpolated at the cut. Output is byte-identical across backends, thread counts and the C++/numpy boundary. See `doc/slice.md`.

<!--pytest-codeblocks:skip-->

```python
vol = meshioplusplus.read("part.vtu")                          # a tetra/hex/wedge mesh
section = meshioplusplus.slice(vol, origin=(0, 0, 0.5), normal=(0, 0, 1))
meshioplusplus.write("section.vtu", section)                  # a triangle/quad surface at z=0.5
```

(`slice` shadows the Python built-in only as a module attribute — `meshioplusplus.slice` is intended.)

These operations are exposed across every binding surface (Python, C API, Fortran, WASM) and as the CLI verbs `meshioplusplus quality`, `meshioplusplus extract-surface`, `meshioplusplus reorder`, `meshioplusplus diff`, `meshioplusplus merge`, `meshioplusplus transform`, `meshioplusplus clean`, `meshioplusplus crop`, `meshioplusplus slice`, `meshioplusplus split`, `meshioplusplus stats`, `meshioplusplus convert-cells`, `meshioplusplus refine`, `meshioplusplus partition`, `meshioplusplus smooth`, and `meshioplusplus interpolate`.

#### Data operations (rename / average / calc / condition / summarize)

A second bundle operates on the **data arrays** a mesh carries (`point_data` / `cell_data` / `field_data`) rather than on its geometry, which none of them ever modifies:

- **`meshioplusplus.data_rename` / `data_drop` / `data_keep`** — rewrite which arrays a mesh carries and under what names; values, dtypes and shapes are copied verbatim. See `doc/data_manage.md`.
- **`meshioplusplus.point_data_to_cell_data` / `cell_data_to_point_data`** — move data between locations by averaging, optionally weighted by cell area/volume. See `doc/data_average.md`.
- **`meshioplusplus.data_calc`** — derive a new array from an elementwise expression (`+ - * /`, parentheses, `abs`/`sqrt`/`min`/`max`/`norm`) evaluated by a hand-written parser — no external parser library, no arbitrary-code path. See `doc/data_calc.md`.
- **`meshioplusplus.data_condition`** — clamp, normalize to a target range, or standardize to zero mean / unit standard deviation, per component or by row magnitude. See `doc/data_condition.md`.
- **`meshioplusplus.data_info`** — a read-only per-array summary (dtype, shape, components, min/max/mean, NaN/inf counts) — the data-side complement to `info` and `compute_stats`. See `doc/data_info.md`.

<!--pytest-codeblocks:skip-->

```python
out = meshioplusplus.data_calc(mesh, "norm(velocity)", location="point", output="speed")
out = meshioplusplus.point_data_to_cell_data(out, keys=["speed"], suffix="_c")
out = meshioplusplus.data_condition(out, "cell", ["speed_c"], mode="normalize")
out = meshioplusplus.data_rename(out, "point", "T", "temperature")
arrays = meshioplusplus.data_info(out)                     # list of per-array dicts
```

These are likewise exposed across every binding surface, and as the nine CLI verbs under the `meshioplusplus data` group (`info`, `rename`, `drop`, `keep`, `to-cell`, `to-point`, `calc`, `clamp`, `normalize`). See `doc/data_operations.md`.

#### Time series

The [XDMF format](https://xdmf.org/index.php/XDMF_Model_and_Format) supports time series with a shared mesh. You can write times series data using meshio++ with

<!--pytest-codeblocks:skip-->

```python
with meshioplusplus.xdmf.TimeSeriesWriter(filename) as writer:
    writer.write_points_cells(points, cells)
    for t in [0.0, 0.1, 0.21]:
        writer.write_data(t, point_data={"phi": data})
```

and read it with

<!--pytest-codeblocks:skip-->

```python
with meshioplusplus.xdmf.TimeSeriesReader(filename) as reader:
    points, cells = reader.read_points_cells()
    for k in range(reader.num_steps):
        t, point_data, cell_data = reader.read_data(k)
```

### Interactive viewer

[![Try it in your browser](https://img.shields.io/badge/try%20it-in%20your%20browser-4c9ffe?logo=webassembly&logoColor=white)](https://loumalouomega.github.io/meshioplusplus/viewer/)

<img alt="the meshio++ browser viewer" src="https://loumalouomega.github.io/meshioplusplus/viewer/browser-viewer.png" width="85%">

*The browser viewer — reading, rendering and converting entirely client-side.*

One call, two backends:

```python
import meshioplusplus

mesh = meshioplusplus.read("part.msh")
meshioplusplus.view(mesh)                      # pick a backend automatically
meshioplusplus.view(mesh, backend="polyscope") # a native desktop window
meshioplusplus.view(mesh, backend="browser")   # vtk.js, in a browser or notebook
```

The **desktop** backend is [Polyscope](https://polyscope.run), an optional Python-only extra (`pip install meshioplusplus[viewer]`). It draws solids you can slice into, colours by any point or cell array, and renders headless screenshots for CI and docs:

```python
meshioplusplus.screenshot(mesh, "part.png", color_by="temperature")
```

<img alt="the example bracket in Polyscope, coloured by scaled Jacobian" src="https://loumalouomega.github.io/meshioplusplus/viewer/desktop-viewer.png" width="70%">

*The bundled `example.msh` bracket coloured by element quality — this image is generated by `screenshot()` itself.*

The **[browser](https://loumalouomega.github.io/meshioplusplus/viewer/)** backend needs nothing extra. The same app is hosted as a live demo: drag in any supported format, colour by point or cell data, and convert and download to another format — all client-side, with no server and no upload. Since it runs the WebAssembly build, every format meshio++ reads works there too.

From the command line:

```sh
meshioplusplus view part.msh
meshioplusplus screenshot part.msh part.png --size 1600 1200
```

See [the viewer docs](https://loumalouomega.github.io/meshioplusplus/viewer.html) for how volume meshes are handled and what each backend can and cannot do.

### ParaView plugin

<img alt="gmsh paraview" src="https://nschloe.github.io/meshio/gmsh-paraview.png" width="60%">
*A Gmsh file opened with ParaView.*

If you have downloaded a binary version of ParaView, you may proceed as follows.

- Install meshio++ for the Python major version that ParaView uses (check `pvpython --version`)
- Open ParaView
- Find the file `paraview-meshioplusplus-plugin.py` of your meshio++ installation (on Linux: `~/.local/share/paraview-5.9/plugins/`) and load it under _Tools / Manage Plugins / Load New_
- _Optional:_ Activate _Auto Load_

You can now open all meshio++-supported files in ParaView.

### Benchmarks

How much does the C++ core help? The [`benchmark/`](https://github.com/loumalouomega/meshioplusplus/tree/main/benchmark) folder times read/write conversions against the original pure-Python [meshio](https://github.com/nschloe/meshio) on the formats both support (same in-memory mesh, same machine). The headline input is the bundled [`example.msh`](https://github.com/loumalouomega/meshioplusplus/blob/main/example/example.msh) — a real Gmsh bracket (~52k nodes, ~293k cells).

<img alt="meshio vs meshio++ speedup on example.msh" src="https://raw.githubusercontent.com/loumalouomega/meshioplusplus/main/benchmark/plots/benchmark_speedup.svg" width="85%">

meshio++'s biggest wins are the parallel and text paths: **VTU binary+zlib ~16× write** (the zlib blocks run across cores via an OpenMP backend with dynamic scheduling — hybrid P+E-core CPUs load-balance too), **VTU ASCII ~7× write / ~5× read**, and mixed-topology **XDMF read ~10×**. The binary and HDF5 formats that used to be *slower* — VTK/Gmsh binary, UGRID, and MED — are now at or above parity after an optimisation pass (bulk-buffered binary I/O, single-instruction `bswap` endianness conversion, a real parallel backend, an Eigen-backed MED transpose, **zero-copy cell reconstruction** that moves the connectivity buffer straight into the mesh, and uninitialised reader buffers + thread-parallel block copies so nothing is written twice); binary **reads** now match or beat numpy's `fromfile` — Gmsh ~1.7×, single-type VTK ~1.45×, and even mixed-topology VTK ~1.1×. Output stays byte-identical throughout.

The speedup is per-element: text/parallel formats climb out of the small-mesh regime and plateau (large meshes realise the full speedup):

<img alt="speedup vs mesh size" src="https://raw.githubusercontent.com/loumalouomega/meshioplusplus/main/benchmark/plots/benchmark_scaling.svg" width="85%">

Full methodology and a reproducible notebook are on the [Benchmarks](https://loumalouomega.github.io/meshioplusplus/benchmarks) doc page (source: [`benchmark/01_benchmark.ipynb`](https://github.com/loumalouomega/meshioplusplus/blob/main/benchmark/01_benchmark.ipynb)).

### Reading only what you need

```python
import meshioplusplus

mesh = meshioplusplus.read("big.vtu", points_only=True)   # geometry, no data arrays
mesh = meshioplusplus.read("big.vtu", arrays=["u", "p"])  # only these arrays
meta = meshioplusplus.read_metadata("big.vtu")            # counts/names, no heavy arrays
```

VTU, VTP, XDMF and Gmsh skip the unwanted array bodies outright; other formats are read in full and filtered, and `meta["fell_back_to_full_read"]` says which happened. Large files can also be memory-mapped (automatic above 16 MiB), which roughly halves peak memory during a read. See [selective reads](doc/selective_read.md) and [memory-mapped reading](doc/mmap.md).

VTK XML output can additionally use **lz4** (ParaView-readable) or **zstd** (a meshio++ extension) instead of zlib, when built with `-DMESHIOPLUSPLUS_WITH_LZ4=ON` / `-DMESHIOPLUSPLUS_WITH_ZSTD=ON`. zlib remains the default. See [compression codecs](doc/codecs.md).

### Installation

meshio++ is [available from the Python Package Index](https://pypi.org/project/meshioplusplus/), so simply run

```
pip install meshioplusplus
```

to install.

Additional dependencies (`netcdf4`, `h5py`) are required for some of the output formats and can be pulled in by

```
pip install meshioplusplus[all]
```

For JavaScript / browser use, the C++ core also ships as a WebAssembly npm package covering 29 of the formats above:

```
npm install @meshioplusplus/wasm
```

See the [WebAssembly / JavaScript](https://loumalouomega.github.io/meshioplusplus/wasm) doc page for usage and the format-support table.

### C / Fortran API

For HPC codes written in C or Fortran, the C++ core also builds as an installable shared library (`libmeshioplusplus`, pure-C99 header, pkg-config + `find_package` support) with a modern OO Fortran 2008 module on top:

```
./build/configure.sh --fortran --tests --build     # --c-api for the C API alone
cmake --install build/cpp-release --prefix /opt/meshioplusplus
```

```c
mio_mesh* m = mio_read("in.msh", NULL);
printf("%lld points\n", (long long)mio_mesh_num_points(m));
mio_write("out.vtu", m, NULL);
mio_mesh_free(m);
```

```fortran
use meshioplusplus
type(mio_mesh) :: m
call m%read("in.msh")
call m%write("out.vtu")
call m%free()
```

The C API is also packaged for **Conan** (root [`conanfile.py`](conanfile.py)) and **vcpkg** (overlay port under [`packages/vcpkg/meshioplusplus/`](packages/vcpkg/meshioplusplus)), both driving the same install/`find_package` path:

```
conan create . -o meshioplusplus/*:with_hdf5=True
vcpkg install meshioplusplus --overlay-ports=ports
```

Full mesh access (build meshes from raw arrays, zero-copy readback) is covered on the [C API](https://loumalouomega.github.io/meshioplusplus/c_api) and [Fortran](https://loumalouomega.github.io/meshioplusplus/fortran) doc pages.

### Single-header C++

The whole C++ core is also amalgamated into one self-contained, [STB](https://github.com/nothings/stb)-style header — [`src/single_include/meshioplusplus/meshioplusplus.hpp`](src/single_include/meshioplusplus/meshioplusplus.hpp) — with pugixml bundled and no external dependencies by default. Drop it in, no CMake or linking required:

```cpp
// in exactly ONE .cpp:
#define MESHIOPLUSPLUS_IMPLEMENTATION
#include "meshioplusplus/meshioplusplus.hpp"
// elsewhere: just #include it (declarations only)
```

```
g++ -std=c++20 -I src/single_include main.cpp
```

It is generated by `./tools/amalgamate.sh` and kept in sync by CI. See the [single-header](https://loumalouomega.github.io/meshioplusplus/single_header) doc page (optional HDF5/netCDF/zlib formats via `MESHIOPLUSPLUS_HAS_*` macros).

[`example/cpp/`](example/cpp/) is the C++ counterpart of [`example/python/`](example/python/): the same tour of meshio++, called directly against this single header instead of the Python bindings, on the [xeus-cpp](https://github.com/compiler-research/xeus-cpp) Jupyter kernel — no PyVista either, renders go through meshio++'s own SVG writer.

### C++ mesh backends

Standalone C++ builds (no Python) can swap the in-memory mesh structure at compile time via `MESHIOPLUSPLUS_MESH_BACKEND` — every format works identically under each backend:

- **MESHIO** (default; the Python extension and PyPI wheels always use it) — mirrors the Python `meshio.Mesh`;
- **NATIVE** — the fastest pure-C++ structure (canonical Float64/Int64 storage, cell-type enum, CSR ragged blocks); the WebAssembly build uses it;
- **KRATOS** — a [Kratos Multiphysics](https://github.com/KratosMultiphysics/Kratos)-style `ModelPart` (Nodes/Elements/Conditions/SubModelParts) plus a header-only templated bridge that populates a real `Kratos::ModelPart` with no Kratos build dependency.

```
./build/configure.sh --mesh-backend NATIVE --tests --build
```

See the [C++ mesh backends](https://loumalouomega.github.io/meshioplusplus/cpp_backends) doc page.

### Testing

To run the meshio++ unit tests, check out this repository, install it with the test extras, and type

```
pytest tests/
```

### License

meshio++ is published under the [MIT license](https://en.wikipedia.org/wiki/MIT_License).

### Acknowledgements

**meshio++ is a fork of [meshio](https://github.com/nschloe/meshio)** by [Nico Schlömer](https://github.com/nschloe) and its contributors (MIT). meshio is where the `Mesh` data model, the cell-type naming, and the great majority of the format readers and writers come from; this fork adds a C++20 core, an operations layer, and the C/Fortran/WebAssembly bindings on top of that foundation. The original copyright is retained in [`LICENSE`](LICENSE).

Some code and test fixtures also come from [Simvia's `meshlane` fork](https://github.com/simvia-tech/meshlane) of meshio (MIT), credited per-change in [`CITATION.cff`](CITATION.cff) and [`CHANGELOG.md`](CHANGELOG.md).

The C++ core is dependency-free by design. Everything below is either optional, bundled, or confined to one binding or tool.

#### Runtime dependencies (Python)

| Project | Used for | License |
| --- | --- | --- |
| [NumPy](https://numpy.org/) | the array type the whole data model is built on | BSD-3-Clause |
| [Rich](https://github.com/Textualize/rich) | CLI output formatting | MIT |

#### Optional dependencies

| Project | Extra | Used for | License |
| --- | --- | --- | --- |
| [Polyscope](https://polyscope.run/) | `[viewer]` | the desktop viewer and `screenshot()` | MIT |
| [h5py](https://www.h5py.org/) | `[all]` | the Python fallback for CGNS, H5M, MED, XDMF | BSD-3-Clause |
| [netCDF4](https://unidata.github.io/netcdf4-python/) | `[all]` | the Python fallback for Exodus | MIT |
| [KaHIP](https://kahip.github.io/) | `[kahip]` / CMake | the quality graph-partitioning backend | MIT |
| [zstandard](https://github.com/indygreg/python-zstandard), [lz4](https://github.com/python-lz4/python-lz4) | `[codecs]` | the Python fallback for the optional VTK block codecs | BSD-3-Clause / BSD-2-Clause |

#### Bundled and build-time

| Project | Used for | License |
| --- | --- | --- |
| [pugixml](https://pugixml.org/) | XML parsing in the C++ core (vendored in `src/cpp/third_party/`) | MIT |
| [Eigen](https://eigen.tuxfamily.org/) | the MED Fortran↔C transpose (git submodule, optional) | MPL-2.0 |
| [pybind11](https://github.com/pybind/pybind11) | the Python bindings | BSD-3-Clause |
| [scikit-build-core](https://github.com/scikit-build/scikit-build-core) | the CMake-driven build backend | Apache-2.0 |
| [Emscripten](https://emscripten.org/) | the WebAssembly build | MIT / NCSA |
| [GoogleTest](https://github.com/google/googletest) | the C++ test suite | BSD-3-Clause |
| [zlib](https://zlib.net/), [Zstandard](https://facebook.github.io/zstd/), [LZ4](https://lz4.org/) | optional compression codecs | zlib / BSD-3-Clause / BSD-2-Clause |
| [HDF5](https://www.hdfgroup.org/solutions/hdf5/), [netCDF](https://www.unidata.ucar.edu/software/netcdf/) | optional native paths for HDF5/netCDF-backed formats | BSD-3-Clause / MIT-like |

#### Browser viewer (isolated under [`src/viewer/`](src/viewer/))

| Project | Used for | License |
| --- | --- | --- |
| [vtk.js](https://kitware.github.io/vtk-js/) | all rendering, colour maps and the scalar bar | BSD-3-Clause |
| [Vite](https://vite.dev/) | the app build | MIT |
| [vite-plugin-singlefile](https://github.com/richardtallent/vite-plugin-singlefile) | the self-contained wheel-bundled build | MIT |
| [Playwright](https://playwright.dev/) | the end-to-end tests and the documentation screenshots | Apache-2.0 |

#### Also

[Kratos Multiphysics](https://github.com/KratosMultiphysics/Kratos) (BSD-3-Clause) is the source of the skin-detection algorithm, the EnSight writer logic, the KaHIP partitioning approach and the `FindKaHIP.cmake` module, all these 3 implementation are from the original author of this project as well. Also thanks to its `ModelPart` design informs the KRATOS mesh backend. [VTK](https://vtk.org/) and [Verdict](https://github.com/sandialabs/verdict) (both BSD-3-Clause) define the mesh-quality formulas. The documentation is built with [VitePress](https://vitepress.dev/) (MIT) and [Doxygen](https://www.doxygen.nl/) (GPL-2.0, used as a tool only). The logo renders the [Stanford Bunny](https://www.thingiverse.com/thing:88208) ("Stanford Bunny — Digitized!" by MakerBot, CC-BY).

Thank you to all of them.

# MMG backend

[MMG](https://www.mmgtools.org/) remeshes simplicial meshes: `mmg2d` planar triangles, `mmgs` triangle surfaces in 3-D and `mmg3d` tetrahedra. meshio++ links it as an optional backend for four operations, adapting a mesh to a size or an anisotropic metric, embedding a level set in it, meshing a planar domain from its boundary and moving it by a displacement. MMG's result cannot be reproduced by a second implementation, so there is no numpy twin, as for [`remesh_volume`](./remesh_volume.md): a build without MMG raises `NotImplementedError` (Python) or throws `Unsupported` (C++), naming the flag.

The file route needs no link at all: meshio++ reads and writes MMG's native files, Medit [`.mesh`](./formats/medit.md) and [`.sol`](./formats/medit_sol.md), so `mmg3d_O3`, `mmgs_O3` and `mmg2d_O3` can be driven through files from any build.

## adapt

`adapt` remeshes a simplicial mesh to a size or a metric (roadmap 3.2.1).

```python
import numpy as np
import meshioplusplus as mio

square = mio.rectangle(dims=(16, 16), cell_type="triangle")
out = mio.adapt(square, hsiz=0.05)                       # a constant edge size
out = mio.adapt(square, size=0.02 + 0.2 * square.points[:, 0])  # a size field

# An anisotropic metric from an analytic field's Hessian
square.point_data["u"] = np.tanh(20 * (square.points[:, 0] - 0.5))
h = mio.hessian(square, "u", location="point")
m = mio.metric_from_hessian(h, "u:hessian", error=1e-3, hmin=1e-3, hmax=0.2)
out, report = mio.adapt(m, metric="metric", return_report=True)
```

- **Which MMG.** A mesh with tetrahedra goes to `MMG3D_mmg3dlib`, with wedges as MMG prisms, triangles and quads as boundary faces and lines as ridges. A planar triangle mesh (2-D points, or 3-D points on one z, which comes back) goes to `MMG2D_mmg2dlib`, with quads passing and lines as boundary edges. Any other triangle surface goes to `MMGS_mmgslib`, with lines as ridges. Hexahedra, pyramids, higher-order and ragged blocks are refused with a pointer to [`convert_cells`](./convert_cells.md), and so are prisms without tetrahedra and quads without triangles.
- **The size.** Give at most one of `hsiz` (a constant edge size), `size` (a size per point: a `point_data` name or an array) and `metric` (a symmetric tensor per point, in the [file order](#symmetric-tensors)). With none, MMG sizes from the mesh itself. [`metric_from_hessian`](#metric_from_hessian) makes a metric from [`hessian`](./hessian.md)'s output; a 3-D metric given to mmg2d keeps its xy block.
- **Parameters.** `hmin` and `hmax` bound the edge size, `hausd` how far the boundary may move, `hgrad` and `hgradreq` the gradation. `angle_detection` is the ridge angle in degrees (45, MMG's default; `None` turns it off). `optim`, `noinsert`, `noswap`, `nomove`, `nosurf`, `nreg`, `opnbdy` and `mem` are MMG's switches of those names; mmgs has no `nosurf` or `opnbdy`. `verbose` defaults to -1, since MMG prints to stdout. A negative or absent value leaves MMG's default.
- **Per-region sizes.** `local={"inlet": {"hmax": 0.01}}` sets `hmin`, `hmax` and `hausd` for the entities of a cell or side region, through `Set_localParameter` (tetrahedra and boundary triangles in mmg3d, triangles in mmgs, triangles and edges in mmg2d).

### Entities

Each element and facet MMG is given carries a reference standing for its regions and its `medit:ref`, so they come back under their names:

- **Cell regions** come back as the cells that carry their reference: an element's region holds the new elements of its subdomain, a facet cell's region the new facets on it. A region that ends up empty stays, empty.
- **Side regions** become boundary facets (triangles in 3-D, edges in 2-D) with their own references. They come back as `(cell, local facet)` pairs on the new elements, through the facet's first owner.
- **A `medit:ref` cell array** comes back with the values of its cells; a `medit:ref` point array is kept as the vertices' references.
- **Required entities.** `frozen` (point ids, a boolean mask or a point region's name) and the regions named in `required` become required vertices, elements and facets, which MMG neither moves nor removes. The regions `medit:RequiredVertices`, `medit:RequiredEdges` and `medit:RequiredTriangles` are required without being named. `medit:Corners` become MMG corners and `medit:Ridges` ridges.
- **Point regions** survive when every point is required: the points are matched by their unchanged coordinates. Any other point region is dropped with a warning.
- **Facet cells.** Boundary triangles (mmg3d) and lines come back when the input had some. MMG's new untagged facets come with them, since MMG returns every boundary facet. MMG gives the facets it creates on an interface the reference of an element, which meshio++ reads as untagged rather than as a region's.

### What does not survive

As for [`remesh`](./remesh.md), the points are new, so other `point_data` and `cell_data` are dropped with a warning; compose with [`interpolate`](./interpolate.md). `field_data` carries through.

### Report

`return_report=True` returns, beside the mesh:

- `status`: `"success"` for `MMG5_SUCCESS`, or `"low_failure"` for `MMG5_LOWFAILURE`, a valid mesh MMG could not finish adapting (a warning, not an error). `MMG5_STRONGFAILURE` raises `RuntimeError`.
- `library` (`"mmg2d"`, `"mmgs"` or `"mmg3d"`), and the points and elements in and out.
- `quality`: the `quality:aspect_ratio` summary of the result (1 is ideal; `min`, `max`, `mean`, `count` and a ten-bin `histogram`).
- `warnings`.

## metric_from_hessian

`metric_from_hessian(mesh, hessian, error=1e-2, hmin=1e-3, hmax=1.0, output="metric")` turns a row-major Hessian (`(n, 9)`, or `(n, 4)` in 2-D) into the metric `M = c/error |H|`. The Hessian is symmetrised, its eigenvalues are taken in absolute value and scaled by `c = (1/2)(d/(d+1))²` for P1 interpolation in `d` dimensions (Frey and Alauzet, *Anisotropic mesh adaptation for CFD computations*, CMAME 194, 2005), and the result is clamped to `[1/hmax², 1/hmin²]`. It is written as `point_data[output]`, six components (three for a 2-D Hessian) in the file order. It needs no MMG and is in every build.

## Other language surfaces

- **C API:** `mio_adapt(mesh, opts, report)` with `mio_adapt_opts` (which embeds the shared `mio_mmg_params`, `mio_mmg_params_init`) and `mio_mmg_report`, and `mio_metric_from_hessian`. The structs are append-only with reserved tails ([C API](./c_api.md)).
- **Fortran:** `mio_adapt(mesh, hsiz=, size_field=, metric=, ..., frozen=, required=, report=)` (1-based `frozen`) and `mio_metric_from_hessian`.
- **Julia:** `adapt(mesh; hsiz, size, metric, ..., frozen, required, local_params)` (1-based `frozen`) and `metric_from_hessian`.
- **R:** `mio_adapt(mesh, hsiz, size, metric, ..., frozen, required)` (1-based `frozen`) and `mio_metric_from_hessian`.
- **CLI:** `meshioplusplus adapt IN OUT --hsiz 0.05` (or `--size NAME`, `--metric NAME`), with `--hmin`, `--hmax`, `--hausd`, `--hgrad`, `--hgradreq`, `--angle-detection`, `--optim`, `--noinsert`, `--noswap`, `--nomove`, `--nosurf`, `--nreg`, `--opnbdy`, `--frozen REGION` and repeatable `--required REGION`, in the Python and native CLIs alike.
- **Pipeline:** an `Adapt` step with `Hsiz`, `Size`, `Metric`, `Hmin`, `Hmax`, `Hausd`, `Hgrad`, `Hgradreq`, `AngleDetection`, the switches, `Frozen` (a point region) and `Required`; it reports `NumPointsIn/Out`, `NumElementsIn/Out`, `LowFailure`, `QualityMin` and `QualityMean` ([pipelines](./pipeline.md)).
- **MCP:** the `adapt` tool, which reports `mmg_available` (an error with `mmg_available: false` without MMG), and `metric_from_hessian` ([MCP](./mcp.md)).
- **WebAssembly:** only `hasMmg()`, always `false`, until roadmap 3.2.8.

## Building with MMG

MMG is bring-your-own: meshio++'s CMake finds an install through the package MMG installs (`mmg_DIR`, or `CMAKE_PREFIX_PATH`), never vendors it and never downloads it. `build/build-mmg.sh` builds the pinned MMG 5.8.0 as shared libraries.

```sh
build/build-mmg.sh ~/mmg-install                 # MMG alone
build/build-mmg.sh ~/mmg-install --with-elas     # also ELAS, for move (GPL-3.0, below)

build/configure.sh --tests --with-mmg --mmg-root ~/mmg-install --build
CMAKE_ARGS="-DMESHIOPLUSPLUS_WITH_MMG=ON -Dmmg_DIR=$HOME/mmg-install/lib/cmake/mmg" pip install .
export LD_LIBRARY_PATH=$HOME/mmg-install/lib:$LD_LIBRARY_PATH   # MMG and ELAS carry no RPATH
```

- **Linked dynamically only.** CMake links the shared `Mmg::libmmg_so`; an install with only the static library is refused with a status message, and the operations then raise.
- **What CMake reads from the install.** The width of `MMG5_int` from `mmg/common/mmgcmakedefines.h` (checked again by a `static_assert` when the core compiles), and whether MMG was built with ELAS, from `libmmg`'s reference to ELAS's `LS_elastic` (MMG's headers do not record it).
- **Never in the published binaries.** The PyPI wheels, the release CLI binaries and both WebAssembly artifacts have no MMG; the browser port is [roadmap 3.2.8](./roadmap.md).

Ask the build what it has:

```python
from meshioplusplus import _core
_core.__has_mmg__, _core.__has_mmg_elas__
```

The C API has `mio_has_mmg()`, `mio_has_mmg_elas()` and `mio_mmg_index_bits()`; Fortran `mio_has_mmg()`, `mio_has_mmg_elas()` and `mio_mmg_index_bits()`; Julia `has_mmg()` and `has_mmg_elas()`; R `mio_has_mmg()` and `mio_has_mmg_elas()`; WebAssembly `hasMmg()`, always `false`; C++ `mmg_available()`, `mmg_has_elas()` and `mmg_index_bits()` in `operations/mmg.hpp`.

## Licence

MMG is **LGPL-3.0-or-later**. meshio++ is MIT, and a build without MMG contains no MMG code. A build with it links `libmmg` dynamically, so a recipient can replace MMG without rebuilding meshio++; whoever redistributes such a build owes the LGPL-3.0's obligations for MMG (its licence text, its source or a written offer of it, and the freedom to relink against a modified `libmmg`). This is why MMG is in none of the published binaries.

**ELAS is GPL-3.0.** `move` needs an MMG built with [ELAS](https://github.com/ISCDtoolbox/LinearElasticity) and its [Commons](https://github.com/ISCDtoolbox/Commons), both GPL-3.0. An MMG linked to them, and a meshio++ linked to that MMG, form a GPL-3.0 work when distributed. `build/build-mmg.sh` builds ELAS only when asked (`--with-elas`); a build without it keeps MMG's LGPL terms and refuses `move` by name. See [third-party notices](https://github.com/loumalouomega/meshioplusplus/blob/master/THIRD_PARTY_NOTICES.md#mmg).

## Embedding MMG

Three hazards of calling MMG as a library, read in the MMG 5.8.0 source, are handled once for every operation (`src/cpp/src/detail/mmg_guard.hpp`):

- **Signal handlers.** Every MMG library entry point (`MMG3D_mmg3dlib`, `MMG3D_mmg3dls`, `MMG3D_mmg3dmov`, `MMGS_mmgslib`, `MMGS_mmgsls`, `MMG2D_mmg2dlib`, `MMG2D_mmg2dmesh`, `MMG2D_mmg2dls`, `MMG2D_mmg2dmov`) calls `signal()` for `SIGABRT`, `SIGFPE`, `SIGILL`, `SIGSEGV`, `SIGTERM` and `SIGINT`, installing a handler that prints a message and calls `exit(EXIT_FAILURE)`, and never puts the old ones back. meshio++ saves the six handlers before each call and restores them after it, on every path including errors, so Python's Ctrl-C handler survives. **During** a call MMG's handlers stand: a Ctrl-C, or a crash inside MMG, ends the whole process rather than raising, so an MMG call cannot be interrupted.
- **Index width.** `MMG5_int` is `int32_t` unless MMG was configured with `-DMMG5_INT=int64_t`. A mesh beyond the range is refused by name, before any data reaches MMG, rather than truncated; `mmg_index_bits()` reports the width.
- **Threads.** Whether two meshes can be adapted at once in one process is not established, so every MMG call is serialised behind one mutex. ParMmg, MMG's MPI layer, is out of scope ([non-goals](./roadmap.md#non-goals-and-decisions-taken)).

## Symmetric tensors

A metric is a symmetric tensor per point, and MMG has two orders for it. Its files ([`.sol`](./formats/medit_sol.md)) store the lower triangle row by row, `m11 m21 m22 m31 m32 m33` (in 2-D `m11 m21 m22`); its API (`MMG3D_Set_tensorSol`) takes the upper triangle row by row, `m11 m12 m13 m22 m23 m33`, and MMG's reader swaps the third and fourth values (`inout.c`). meshio++'s metric arrays are in the **file order**, so a metric read from a `.sol` is used as it is, and the operations swap for the API.

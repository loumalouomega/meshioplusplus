# Field-aligned quads — `remesh_quads`

`remesh_quads` remeshes a triangle surface into quads that follow a smooth cross field, so they line up with the surface's principal directions and stay square everywhere, a sphere's poles included. It hands the surface to a quad remesher, an installed program run as a separate process: **QuadriFlow** (J. Huang, Y. Zhou, M. Niessner, J. R. Shewchuk and L. Guibas, *QuadriFlow: a scalable and robust method for quadrangulation*, SGP 2018; MIT) or **Instant Meshes** (W. Jakob, M. Tarini, D. Panozzo and O. Sorkine-Hornung, *Instant field-aligned meshes*, SIGGRAPH Asia 2015; BSD-3). The result is a new surface near the input, not on its points; the report says how near.

```python
import meshioplusplus as mp

skin = mp.mesh_cad_surface(mp.read("sphere.cad.json"), 0.1, edges=False)
quads, report = mp.remesh_quads(skin, return_report=True)
print(report["num_quads"], report["min_scaled_jacobian"], report["hausdorff_to_surface"])
```

## Why not recombine

[`mesh_cad_surface(quads="all")`](./mesh_cad.md) recombines the mesher's own triangles into quads by a maximum matching. That keeps the triangulation's points and conforms to the CAD faces, but the quads inherit the triangles' directions: where those swirl, at a sphere's poles, the quads are poor. A field-aligned remesher places new points along a cross field smoothed over the whole surface, so its quads are square wherever the field is smooth and leave the field's few singularities, eight on a sphere, as the only irregular vertices.

## Measured

On the unit sphere meshed by `mesh_cad_surface` at `h = 0.1` (5,540 triangles), QuadriFlow at its default target of half the triangles makes 2,602 quads, no triangles, whose least scaled Jacobian is 0.368 and mean 0.987, within 0.0037 of the input surface (seed 0; seed 1 gives 2,571 quads, 0.405). Recombining the same triangles (`quads="all"`) makes 2,778 quads whose least scaled Jacobian is 0.116, at the poles. The worst quad depends on the seed, since it sits at one of the cross field's singularities: on the sphere of the CAD JSON fixture (5,446 triangles) seeds 0, 1 and 2 give 0.25, 0.58 and 0.48, the mean 0.987 every time. So the roadmap's probe, at least 0.3, holds for most seeds but not all; `tests/cpp/test_remesh_quads.cpp` and `tests/python/test_remesh_quads.py` pin what holds for every seed tried (all quads, a mean of at least 0.95, a worst of at least 0.2 and above recombination's), and CI's external-meshers job builds QuadriFlow from its pinned commit to run them.

## The remeshers

A quad remesher is a [provider](./cpp_api.md#providers-and-plugins) registered by name (`QuadRemesherProvider`, `add_quad_remesher`, `quad_remesher_names()`). Every native build registers the two programs when it finds them, on `PATH` or named by `MESHIOPLUSPLUS_QUADRIFLOW` and `MESHIOPLUSPLUS_INSTANT_MESHES` (`MESHIOPLUSPLUS_NO_EXTERNAL_MESHERS=1` turns them off, with the [external tetrahedral meshers](./tetrahedralize.md#external-meshers-tetgen-and-gmsh)); `remesh_quads_backends()` lists what this install has, and `backend="auto"` takes QuadriFlow, then Instant Meshes.

- **`quadriflow`** gets the surface as a full-precision OBJ and runs `quadriflow -i in.obj -o out.obj -f <target_faces>` with `-sharp`, `-boundary`, `-adaptive` and `-seed` as asked. It always returns quads only. It needs only Eigen and Boost to build (`cmake -DBUILD_FREE_LICENSE=ON`); no Linux distribution packages it.
- **`instant_meshes`** runs the `Instant Meshes` program (also found as `instant-meshes` or `InstantMeshes`) in batch mode, `-d -o out.obj -f <target_faces>`, deterministic, with `-D` (quad-dominant, triangles allowed) for `dominant`, `-b` for `boundary` and `-c 30` for `sharp`. It builds with its GUI library, so a build needs OpenGL and X headers. It is not exercised by CI, nor was it run while this was written: its build failed here on its bundled TBB with GCC 13.

The program's OBJ is read back by its `v` and `f` records: quads, and the triangles a quad-dominant run leaves (a larger face is refused). Every face must name a point and every point be finite, or the call fails naming the program; a program that fails quotes the last lines of its output. Each run works in a fresh temporary directory, removed afterwards (`MESHIOPLUSPLUS_KEEP_TEMP=1` keeps it). They never run in WebAssembly.

Running them as programs keeps meshio++ MIT: neither is linked, vendored or shipped. A linked add-on (the remeshers' solvers extracted into a library behind a `meshioplusplus-quadriflow` wheel, as the Netgen and MMG add-ons are built) would bring them to the wheels and the browser; it was weighed against this route and not taken, since it is most of the work and the programs already run everywhere else.

## Options and report

| option | meaning |
|---|---|
| `target_faces` | the number of faces asked for; 0 takes half the input's triangles (the programs treat it as a target) |
| `backend` | `auto` (QuadriFlow, else Instant Meshes) or a remesher's name |
| `sharp` | align with and keep the sharp edges |
| `boundary` | keep the boundary of an open surface (true) |
| `adaptive` | let the face size follow the curvature (QuadriFlow) |
| `dominant` | allow triangles where quads fit badly (Instant Meshes) |
| `seed` | the program's random seed, so a run repeats |
| `max_cells` | refuse a target or result of more faces (50,000,000) |

The input is linear triangles: lines and vertices are skipped, any other cell type is refused by name. The report has `num_quads`, `num_triangles` (left among the quads, with a warning), `num_points`, `min_scaled_jacobian` and `mean_scaled_jacobian` (Verdict's, the measure [`mesh_cad_surface`](./mesh_cad.md) reports for its quads), `hausdorff_to_surface` (the two-sided distance to the input, sampled), `backend` and `warnings`. The input's length unit carries over; its data and regions do not, since the points are new.

## Bindings

| surface | entry |
|---|---|
| Python | `meshioplusplus.remesh_quads(mesh, *, target_faces=0, backend="auto", sharp=False, boundary=True, adaptive=False, dominant=False, seed=0, max_cells=50_000_000, return_report=False)`, `remesh_quads_backends()`; no numpy twin |
| C++ | `remesh_quads(const Mesh&, const RemeshQuadsOptions&)` → `RemeshQuadsResult`, `remesh_quads_backends()` (`operations/remesh_quads.hpp`); `QuadRemesherProvider`, `add_quad_remesher`, `quad_remesher(name)`, `quad_remesher_names()` (`providers.hpp`) |
| C | `mio_remesh_quads(mesh, &opts, &report)` with `mio_remesh_quads_opts` (through `mio_remesh_quads_opts_init`) and `mio_remesh_quads_report` (`backend_name`); `mio_quad_remesher_count()`, `mio_quad_remesher_name(i)` |
| Fortran | `mio_remesh_quads(mesh, ...)`, `mio_quad_remesher_count()`, `mio_quad_remesher_name(i)`, `mio_remesh_quads_backend_name(report)` |
| Julia | `remesh_quads(mesh; target_faces, backend, sharp, boundary, adaptive, dominant, seed, max_cells)`, `remesh_quads_backends()` |
| R | `mio_remesh_quads(mesh, ...)`, `mio_remesh_quads_backends()` |
| WASM | `remeshQuads(mesh, {targetFaces, backend, sharp, boundary, adaptive, dominant, seed, maxCells})`, which always throws there: no program runs in WebAssembly |
| CLI | `remesh-quads INFILE OUTFILE [--faces N] [--backend NAME] [--sharp] [--no-boundary] [--adaptive] [--dominant] [--seed S] [--max-cells N] [--json]`, both CLIs; `meshioplusplus --version` lists the remeshers |
| MCP | `remesh_quads(input_path, output_path, ...)`, reporting `backends_available` |
| Pipeline | `RemeshQuads` with `TargetFaces`, `Backend`, `Sharp`, `Boundary`, `Adaptive`, `Dominant`, `Seed`, `MaxCells` |

# Tetrahedra keeping a skin — `tetrahedralize`

`tetrahedralize` fills a closed triangle surface with tetrahedra and keeps the surface's triangles exactly: every input triangle is a facet of the result, every input point is a point of it at its index. It is the conforming volume route for a [`mesh_cad_surface`](./mesh_cad.md) result or a [boundary layer](./boundary_layer.md)'s cap, where [`octree_mesh`](./octree_mesh.md) fits its tetrahedra to the skin within a measured distance instead. The tetrahedra come from an optional backend, Netgen's nglib, in builds configured with `-DMESHIOPLUSPLUS_WITH_NETGEN=ON`; without it the operation throws `Unsupported` naming the flag.

```python
import meshioplusplus as mp

model = mp.read("part.cad.json")
skin = mp.mesh_cad_surface(model, 0.2, min_angle=25)
volume, report = mp.tetrahedralize(skin, max_size=0.2, return_report=True)
print(report["num_tetrahedra"], report["boundary_preserved"], report["min_dihedral_deg"])
```

## Input

The input is a closed, manifold, consistently oriented triangle surface. Line and vertex blocks are skipped (a `mesh_cad_surface` result's edges); a quad, polygon, higher-order or volume block is refused by name, so an all-quad or curved surface is meshed from its linear triangle result, and a quad cap is split first. An open, non-manifold, misoriented or degenerate surface is refused with the counts of the watertightness check, the same counts `mesh_cad_surface` reports. A self-intersecting surface passes that check and fails in the backend, reported as an error rather than meshed wrong.

## Regions

The surface's connected shells are oriented outward (a shell whose triangles face in is turned), and each shell bounds a region: the inside of the shell minus the shells directly inside it, found by the winding number of one of their points. A box inside a box gives two regions, the shell between them and the inner box. Each region is meshed on its own, from its shell and its inner shells turned inward, so neighbouring regions share the inner shell's triangles exactly and the tetrahedra conform across it.

The result holds the input's points first, unchanged and at their indices, then the new interior points; a `tetra` block (positively oriented, region by region) and the surface's `triangle` block as given. The cell data `tetrahedralize:region` names each tetrahedron's region (-1 on the triangles), and Cell regions `region:<k>` hold them. When the surface carries `face_array` (`cad:face`, one value per CAD face), that array is carried to the triangles (-1 on the tetrahedra) and each value becomes a Side region `face:<v>` on the facets the triangles are.

## The backend

Netgen (J. Schöberl, 1997) meshes the volume by an advancing front with Delaunay insertion and then optimises it: its volume mesher keeps the surface mesh it is given, which `tetrahedralize` checks (`boundary_preserved`, with `hausdorff_to_surface` zero when it holds), and `optimize_steps` passes of its volume optimisation follow. `max_size` bounds the tetrahedra's edges, `max_volume` adds the edge of a regular tetrahedron of that volume, `(6 √2 V)^(1/3)`, to the bound, and `grading` sets how fast the size grows away from the surface; with neither bound the surface's triangles size the inside. nglib keeps global state, so calls into it are serialised behind one lock and a region is meshed in one thread: the result is deterministic, the same bytes run after run.

Netgen is LGPL-2.1. It is linked dynamically, never vendored, and never in the published wheels or binaries ([third-party notices](https://github.com/loumalouomega/meshioplusplus/blob/main/THIRD_PARTY_NOTICES.md)); `build/configure.sh --with-netgen` (or `--netgen-root <prefix>`) enables it against a distribution's `libnglib-dev` or an install of its own. The WebAssembly build has no Netgen: `tetrahedralize` throws there by name and `hasNetgen()` is false.

## Measured

On the unit sphere meshed by `mesh_cad_surface` at `h = 0.2` and 25 degrees, `tetrahedralize` makes 3,688 tetrahedra whose boundary is the skin itself, with dihedral angles of 8.6 to 164.4 degrees; their volume is the skin's enclosed volume to rounding, so `4/3 π` is met to the skin's own chordal error. A box inside a box gives two regions that share the inner box's triangles. `tests/cpp/test_tetrahedralize.cpp` and `tests/python/test_tetrahedralize.py` pin these, and CI's `netgen` job runs them against Ubuntu's Netgen.

[`optimize_volume`](./optimize_volume.md) (ODT relocation with 2-3 and 3-2 flips) runs on the result unchanged; what it reaches on Netgen's output is recorded below. The moves the literature's strongest schedules add (edge removal, multi-face removal; Klingner and Shewchuk 2007) are not in `optimize_volume`.

## Options and report

| option | meaning |
|---|---|
| `backend` | `auto` (the first built in) or `netgen` |
| `max_size` | the largest edge asked of the tetrahedra (0: none) |
| `max_volume` | the largest volume asked of a tetrahedron (0: none) |
| `grading` | how fast the size grows away from the surface, in `[0, 1]` (0.3) |
| `optimize_steps` | passes of the backend's volume optimisation, 0 to 100 (3) |
| `max_cells` | refuse a result of more tetrahedra (50,000,000) |
| `face_array` | the cell data naming each triangle's CAD face (`cad:face`) |

The report has `num_tetrahedra`, `num_points` (every point), `num_new_points` (of those, inside the surface), `num_regions`, `num_boundary_triangles`, `boundary_preserved`, `hausdorff_to_surface`, `min_dihedral_deg` and `max_dihedral_deg`, `volume` (the tetrahedra's total), `backend` (the one that ran) and `warnings`.

## Other backends

Geogram (BSD-3) was the first candidate and was rejected: its `mesh_tetrahedralize`, the constrained tetrahedralisation of a closed surface, creates its Delaunay kernel as `Delaunay::create(3, "tetgen")`, and TetGen is AGPL-3.0, which this repository excludes with CGAL's `Mesh_3` and Gmsh (GPL). Two candidates stay recorded for a second backend under the same contract (`TetrahedralizeBackend`):

- **fTetWild** (MPL-2.0; Hu et al. 2020) meshes within an envelope around the surface, so it survives dirty input and does not keep the triangles exactly: it would report `boundary_preserved` false and a Hausdorff distance. It depends on Geogram's Delaunay kernel (without TetGen) and libigl.
- **The exact constrained Delaunay tetrahedralisation** of Diazzi, Panozzo, Vaxman and Attene (2023), on Attene's indirect predicates (Attene 2020), recovers every input triangle without approximation and meshed all valid Thingi10K models in its authors' tests. Its reference code is LGPL only when built with `-DLGPL=ON` (its default links GPL components), and it ships as an executable, so it would be a file route like the [vendor routes](./vendor_routes.md), not a link.

A native constrained tetrahedralisation stays a non-goal: the Delaunay kernel is short, but boundary recovery with Steiner points is where the robust libraries spend their years ([non-goals](./roadmap.md#non-goals-and-decisions-taken)).

## Bindings

| surface | entry |
|---|---|
| Python | `meshioplusplus.tetrahedralize(mesh, *, backend="auto", max_size=0.0, max_volume=0.0, grading=0.3, optimize_steps=3, max_cells=50_000_000, face_array="cad:face", return_report=False)`; `meshioplusplus._core.__has_netgen__`; no numpy twin |
| C++ | `tetrahedralize(const Mesh&, const TetrahedralizeOptions&)` → `TetrahedralizeResult`, `tetrahedralize_has_netgen()` (`operations/tetrahedralize.hpp`) |
| C | `mio_tetrahedralize(mesh, &opts, &report)` with `mio_tetrahedralize_opts` (through `mio_tetrahedralize_opts_init`) and `mio_tetrahedralize_report`; `mio_has_netgen()` |
| Fortran | `mio_tetrahedralize(...)`, see [Fortran](./fortran.md) |
| Julia | `tetrahedralize(mesh; backend, max_size, max_volume, grading, optimize_steps, max_cells, face_array)` |
| R | `mio_tetrahedralize(mesh, ...)` |
| WASM | `tetrahedralize(mesh, {...})` throws by name; `hasNetgen()` is `false` |
| CLI | `tetrahedralize INFILE OUTFILE [--max-size H] [--max-volume V] [--grading G] [--optimize-steps N] [--max-cells N] [--face-array NAME] [--backend auto\|netgen] [--json]`, both CLIs |
| MCP | `tetrahedralize(input_path, output_path, ...)` |
| Pipeline | `Tetrahedralize` with `Backend`, `MaxSize`, `MaxVolume`, `Grading`, `OptimizeSteps`, `MaxCells`, `FaceArray` |

## References

- J. Schöberl, "NETGEN: an advancing front 2D/3D-mesh generator based on abstract rules", *Computing and Visualization in Science* 1 (1997).
- B. M. Klingner and J. R. Shewchuk, "Aggressive tetrahedral mesh improvement", *Proceedings of the 16th International Meshing Roundtable* (2007).
- Y. Hu, T. Schneider, B. Wang, D. Zorin and D. Panozzo, "Fast tetrahedral meshing in the wild", *ACM Transactions on Graphics* 39 (2020).
- L. Diazzi, D. Panozzo, A. Vaxman and M. Attene, "Constrained Delaunay tetrahedrization: a robust and practical approach", *ACM Transactions on Graphics* 42 (2023).
- M. Attene, "Indirect predicates for geometric constructions", *Computer-Aided Design* 126 (2020).
- Q. Zhou and A. Jacobson, "Thingi10K: a dataset of 10,000 3D-printing models", arXiv:1605.04797 (2016).

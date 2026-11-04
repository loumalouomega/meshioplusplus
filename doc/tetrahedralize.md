# Tetrahedra from a skin — `tetrahedralize`

`tetrahedralize` fills a closed triangle surface with tetrahedra: a [`mesh_cad_surface`](./mesh_cad.md) result, a [boundary layer](./boundary_layer.md)'s cap, any watertight skin. It has two backends. **Netgen** keeps the surface's triangles exactly: every input triangle is a facet of the result and every input point a point of it at its index, so the volume conforms to the skin and to whatever else shares it. The **octree** backend is native and always built: [`octree_mesh`](./octree_mesh.md)'s isosurface stuffing, whose boundary is fitted to the skin within a measured distance rather than kept, with dihedral angles bounded where the stuffing's guarantee holds. `backend="auto"`, the default, takes Netgen when the build has it (`-DMESHIOPLUSPLUS_WITH_NETGEN=ON`) and the octree backend otherwise, with a warning saying so; `backend="netgen"` in a build without it throws `Unsupported` naming the flag. So `tetrahedralize` works in every wheel, CLI and WebAssembly build, and Netgen is the upgrade that keeps the skin. The installed **TetGen** and **Gmsh** programs are two more backends that keep it (`backend="tetgen"`, `"gmsh"`): `tetrahedralize` writes each region's surface, runs the program as a separate process and reads the tetrahedra back ([external meshers](#external-meshers-tetgen-and-gmsh)).

```python
import meshioplusplus as mp

model = mp.read("part.cad.json")
skin = mp.mesh_cad_surface(model, 0.2, min_angle=25)
volume, report = mp.tetrahedralize(skin, max_size=0.2, return_report=True)
print(report["num_tetrahedra"], report["boundary_preserved"], report["min_dihedral_deg"])
```

## Input

The input is a closed, manifold, consistently oriented triangle surface. Line and vertex blocks are skipped (a `mesh_cad_surface` result's edges); a quad, polygon, higher-order or volume block is refused by name, so an all-quad or curved surface is meshed from its linear triangle result, and a quad cap is split first. An open, non-manifold, misoriented or degenerate surface is refused with the counts of the watertightness check, the same counts `mesh_cad_surface` reports. A self-intersecting surface passes that check and fails in the backend, reported as an error rather than meshed wrong; [`find_intersections`](/find_intersections) finds the crossing triangles beforehand.

## Regions

The surface's connected shells are oriented outward (a shell whose triangles face in is turned), and each shell bounds a region: the inside of the shell minus the shells directly inside it, found by the winding number of one of their points. A box inside a box gives two regions, the shell between them and the inner box. Each region is meshed on its own, from its shell and its inner shells turned inward. With Netgen neighbouring regions share the inner shell's triangles exactly and the tetrahedra conform across it; the octree backend fits each region to its own shells, so two regions do not share points across their interface, which it warns.

The result holds the points, a `tetra` block (positively oriented, region by region) and a `triangle` block. With Netgen the input's points come first, unchanged and at their indices, then the new interior points, and the triangles are the surface's as given; with the octree backend every point is new and the triangles are the tetrahedra's boundary facets, outward. The cell data `tetrahedralize:region` names each tetrahedron's region (-1 on the triangles), and Cell regions `region:<k>` hold them. When the surface carries `face_array` (`cad:face`, one value per CAD face), that array is carried to the triangles (-1 on the tetrahedra; an octree boundary facet takes the value of the surface triangle nearest its centroid) and each value becomes a Side region `face:<v>` on the facets the triangles are.

## Netgen: the skin kept

Netgen (J. Schöberl, 1997) meshes the volume by an advancing front with Delaunay insertion and then optimises it: its volume mesher keeps the surface mesh it is given, which `tetrahedralize` checks (`boundary_preserved`, with `hausdorff_to_surface` zero when it holds), and `optimize_steps` passes of its volume optimisation follow. `max_size` bounds the tetrahedra's edges, `max_volume` adds the edge of a regular tetrahedron of that volume, `(6 √2 V)^(1/3)`, to the bound, and `grading` sets how fast the size grows away from the surface; with neither bound the surface's triangles size the inside. nglib keeps global state, so calls into it are serialised behind one lock and a region is meshed in one thread: the result is deterministic, the same bytes run after run.

Netgen is LGPL-2.1. Natively it is linked dynamically, never vendored, and not in the base wheel or the CLI binaries ([third-party notices](https://github.com/loumalouomega/meshioplusplus/blob/main/THIRD_PARTY_NOTICES.md)): the PyPI route is the `meshioplusplus-netgen` add-on wheel (`pip install "meshioplusplus[netgen]"`, `has_netgen()` reports it), and the CLI and C API can load it as a [plugin](./cpp_api.md#providers-and-plugins); `build/configure.sh --with-netgen` (or `--netgen-root <prefix>`) enables it against a distribution's `libnglib-dev` or an install of its own. In the browser the [full WebAssembly artifact](./wasm.md#lean-and-full-builds) links Netgen 6.2.2401 statically (`hasNetgen()` is true there), built from source by `build/build-wasm-deps.sh --full`; the lean artifact has none, so `tetrahedralize` there is the octree backend.

## External meshers: TetGen and Gmsh

A tetrahedral mesher is a [provider](./cpp_api.md#providers-and-plugins), and besides the default one (Netgen) any number can be registered by name and chosen with `backend="<name>"` (C++ `TetrahedralizeBackend::Named` and `TetrahedralizeOptions::mProvider`; `tetrahedralize_backends()` lists them). Two come with every native build and register themselves when their programs are installed, found on `PATH` or named by `MESHIOPLUSPLUS_TETGEN` and `MESHIOPLUSPLUS_GMSH` (`MESHIOPLUSPLUS_NO_EXTERNAL_MESHERS=1` turns both off):

- **`tetgen`** writes the region's surface as a TetGen [`.smesh`](./formats/tetgen.md) (its inner shells as holes, a point inside each) and runs `tetgen -pYQ`: a constrained tetrahedralisation that adds no point on the boundary (`-Y`), with quality (`-q`) when `optimize_steps` is positive and a volume bound (`-a`) from `max_size`, the volume of a regular tetrahedron of that edge. It reads back `.1.node` and `.1.ele`.
- **`gmsh`** writes the surface as a full-precision ASCII STL, one `solid` per shell, with a `.geo` that makes each shell a surface loop of one volume, and runs `gmsh -3` (`-clmax` from `max_size`; `optimize_steps` 0 turns its optimisation off). It reads back the MSH 4.1 output's tetrahedra.

The core checks what either returns, as it checks Netgen and any provider since: every input point back, matched to the program's copy of it within 1e-9 of the region's diagonal (a text round trip rounds the last digit) and then given the input's exact coordinates, and the faces of exactly one tetrahedron the input's triangles. A program that moves the boundary is refused naming it, so `boundary_preserved` is true and `hausdorff_to_surface` zero whenever the call returns. `grading` is Netgen's alone. The programs run in a fresh temporary directory, removed afterwards (`MESHIOPLUSPLUS_KEEP_TEMP=1` keeps it), with no shell, their output in a log whose last lines a failure quotes. They never run in WebAssembly, where the names are not registered.

Running them as programs, never linking them, is what their licences leave open: TetGen is AGPL-3.0 and Gmsh GPL-2.0-or-later, and meshio++ neither ships nor links either (the [decision](./roadmap.md#non-goals-and-decisions-taken) on them stands for linking). On Ubuntu `apt install tetgen gmsh` provides both, and CI's `netgen` job runs the tests against them.

## Octree: the skin fitted

The octree backend meshes each region with [`octree_mesh`](./octree_mesh.md) under uniform skin grading at the surface's own size (its mean edge length, or `max_size` when that is smaller): a body-centred cubic lattice stuffed to the skin, its vertices warped onto it and the result fitted, untangled and optimised for `optimize_steps` rounds (F. Labelle and J. R. Shewchuk, *Isosurface stuffing*, 2007). Where the report of `octree_mesh` says its guarantee holds the dihedral angles lie between 10.7 and 164.8 degrees by construction; the boundary is within `hausdorff_to_surface` of the skin, the two-sided distance between them, and `boundary_preserved` is false. The inside grows to `max_size` or an eighth of the region's diagonal; `grading` is Netgen's alone. It is native, deterministic across parallel backends and thread counts, and in every build.

## Measured

On the unit sphere meshed by `mesh_cad_surface` at `h = 0.2` and 25 degrees, `tetrahedralize` makes 3,145 tetrahedra whose boundary is the skin itself, with dihedral angles of 8.5 to 164.7 degrees (3,680 and 8.9 to 164.4 with `max_size=0.2`); their volume, 4.1500, is the skin's enclosed volume to rounding, 0.9% under `4/3 π`, the skin's own chordal error. A box inside a box gives two regions that share the inner box's triangles. The octree backend fills the same sphere with 15,480 tetrahedra whose dihedral angles lie between 21.4 and 134.3 degrees, its boundary within 0.0097 (0.05 `h`) of the skin and its volume, 4.1252, 0.6% under the skin's; the unit cube's skin it fills exactly, with angles of 45 to 98 degrees. The two trade as their names say: Netgen keeps the skin with a fifth of the tetrahedra and wider angle extremes, the octree backend bounds the angles and moves the boundary. `tests/cpp/test_tetrahedralize.cpp` and `tests/python/test_tetrahedralize.py` pin these, and CI's `netgen` job runs the Netgen tests against Ubuntu's Netgen.

[`optimize_volume`](./optimize_volume.md) takes the tetrahedra once the triangle block is dropped. Its ODT energy objective (ODT relocation with 2-3 and 3-2 flips, Alliez et al. 2005; `objective="energy"`, the only one up to v18.0.0) was measured on Netgen's output (October 2026): ODT relocation lowers an energy summed over the cells, not the worst one, so it trades extremes for the average.

| skin (`h`) | Netgen, `optimize_steps` | Netgen's dihedrals | slivers under 15° | after `optimize_volume` | slivers |
|---|---:|---|---:|---|---:|
| sphere (0.2) | 0 | 0.84 – 177.5 | 138 | 0.18 – 179.6 | 70 |
| sphere (0.2) | 3 | 8.5 – 164.7 | 5 | 1.65 – 176.0 | 39 |
| sphere (0.1) | 0 | 0.05 – 179.8 | 1,111 | 0.05 – 179.8 | 518 |
| sphere (0.1) | 3 | 5.2 – 167.7 | 48 | 2.0 – 174.1 | 336 |
| capped cylinder (0.2) | 0 | 1.7 – 176.3 | 133 | 6.6 – 169.3 | 23 |
| capped cylinder (0.2) | 3 | 13.0 – 160.3 | 2 | 2.0 – 176.8 | 6 |
| unit cube (0.2) | 3 | 23.5 – 137.3 | 0 | 16.2 – 149.3 | 0 |

On Netgen's raw output the energy objective halves the slivers; after Netgen's three passes it makes the extremes worse on every skin. The default objective since, `worst` (guarded relocation, edge removal, multi-face removal and worst-cell smoothing; Klingner and Shewchuk 2007, Freitag and Ollivier-Gooch 1997), never lowers the least dihedral angle and raises it past Netgen's own passes: 17.4 to 32.4 degrees on the sphere and 13.3 to 19.6 on the capped cylinder, with no tet left under 15 degrees ([measured](./optimize_volume.md#measured)).

![The unit sphere's skin at h = 0.2 filled by tetrahedralize's octree backend and by Netgen, each cut open at x = 0 and coloured by each tetrahedron's smallest dihedral angle, beside the histograms of those angles](/images/tetrahedralize_sphere.png)

## Options and report

| option | meaning |
|---|---|
| `backend` | `auto` (Netgen when built in, else the octree backend, with a warning), `netgen`, `octree`, or a registered mesher's name (`tetgen`, `gmsh`; `tetrahedralize_backends()`) |
| `max_size` | the largest edge asked of the tetrahedra (0: none) |
| `max_volume` | the largest volume asked of a tetrahedron (0: none) |
| `grading` | how fast the size grows away from the surface, in `[0, 1]` (0.3; Netgen) |
| `optimize_steps` | passes of the backend's volume optimisation, 0 to 100 (3) |
| `max_cells` | refuse a result of more tetrahedra (50,000,000) |
| `face_array` | the cell data naming each triangle's CAD face (`cad:face`) |

The report has `num_tetrahedra`, `num_points` (every point), `num_new_points` (of those, not the input's: every point with the octree backend), `num_regions`, `num_boundary_triangles`, `boundary_preserved`, `hausdorff_to_surface`, `min_dihedral_deg` and `max_dihedral_deg`, `volume` (the tetrahedra's total), `backend` (the one that ran) and `warnings`.

## Other backends

Geogram (BSD-3) was the first candidate and was rejected: its `mesh_tetrahedralize`, the constrained tetrahedralisation of a closed surface, creates its Delaunay kernel as `Delaunay::create(3, "tetgen")`, and TetGen is AGPL-3.0, which this repository excludes from what it links, with CGAL's `Mesh_3` and Gmsh (GPL); TetGen and Gmsh are reachable as [external meshers](#external-meshers-tetgen-and-gmsh) instead. Two candidates stay recorded for another backend under the same contract (a named `TetMesherProvider`):

- **fTetWild** (MPL-2.0; Hu et al. 2020) meshes within an envelope around the surface, so it survives dirty input and does not keep the triangles exactly: it would report `boundary_preserved` false and a Hausdorff distance. It depends on Geogram's Delaunay kernel (without TetGen) and libigl.
- **The exact constrained Delaunay tetrahedralisation** of Diazzi, Panozzo, Vaxman and Attene (2023), on Attene's indirect predicates (Attene 2020), recovers every input triangle without approximation and meshed all valid Thingi10K models in its authors' tests. Its reference code is LGPL only when built with `-DLGPL=ON` (its default links GPL components), and it ships as an executable, so it would be a file route like the [vendor routes](./vendor_routes.md), not a link.

A native constrained tetrahedralisation stays a non-goal: the Delaunay kernel is short, but boundary recovery with Steiner points is where the robust libraries spend their years ([non-goals](./roadmap.md#non-goals-and-decisions-taken)).

## Bindings

| surface | entry |
|---|---|
| Python | `meshioplusplus.tetrahedralize(mesh, *, backend="auto", max_size=0.0, max_volume=0.0, grading=0.3, optimize_steps=3, max_cells=50_000_000, face_array="cad:face", return_report=False)`; `has_netgen()`, `tetrahedralize_backends()`; no numpy twin |
| C++ | `tetrahedralize(const Mesh&, const TetrahedralizeOptions&)` → `TetrahedralizeResult`, `tetrahedralize_has_netgen()`, `tetrahedralize_backends()` (`operations/tetrahedralize.hpp`); `add_tet_mesher`, `tet_mesher(name)`, `tet_mesher_names()` (`providers.hpp`) |
| C | `mio_tetrahedralize(mesh, &opts, &report)` with `mio_tetrahedralize_opts` (through `mio_tetrahedralize_opts_init`; `MIO_TETRAHEDRALIZE_BACKEND_NAMED` with `provider`) and `mio_tetrahedralize_report` (`backend_name`); `mio_has_netgen()`, `mio_tet_mesher_count()`, `mio_tet_mesher_name(i)` |
| Fortran | `mio_tetrahedralize(...)`, `mio_tet_mesher_count()`, `mio_tet_mesher_name(i)`, `mio_tetrahedralize_backend_name(report)`, see [Fortran](./fortran.md) |
| Julia | `tetrahedralize(mesh; backend, max_size, max_volume, grading, optimize_steps, max_cells, face_array)`, `tetrahedralize_backends()` |
| R | `mio_tetrahedralize(mesh, ...)`, `mio_tetrahedralize_backends()` |
| WASM | `tetrahedralize(mesh, {backend, maxSize, maxVolume, grading, optimizeSteps, maxCells, faceArray})`: Netgen in the full artifact (`hasNetgen()` is `true`), the octree backend in the lean one, where `backend: "netgen"` throws naming the full artifact; no external mesher |
| CLI | `tetrahedralize INFILE OUTFILE [--max-size H] [--max-volume V] [--grading G] [--optimize-steps N] [--max-cells N] [--face-array NAME] [--backend NAME] [--json]`, both CLIs; `meshioplusplus --version` lists the backends |
| MCP | `tetrahedralize(input_path, output_path, ...)`, reporting `netgen_available` and `backends_available` |
| Pipeline | `Tetrahedralize` with `Backend`, `MaxSize`, `MaxVolume`, `Grading`, `OptimizeSteps`, `MaxCells`, `FaceArray` |

## References

- J. Schöberl, "NETGEN: an advancing front 2D/3D-mesh generator based on abstract rules", *Computing and Visualization in Science* 1 (1997).
- F. Labelle and J. R. Shewchuk, "Isosurface stuffing: fast tetrahedral meshes with good dihedral angles", *ACM Transactions on Graphics* 26 (2007).
- B. M. Klingner and J. R. Shewchuk, "Aggressive tetrahedral mesh improvement", *Proceedings of the 16th International Meshing Roundtable* (2007).
- P. Alliez, D. Cohen-Steiner, M. Yvinec and M. Desbrun, "Variational tetrahedral meshing", *ACM Transactions on Graphics* 24 (2005).
- L. A. Freitag and C. Ollivier-Gooch, "Tetrahedral mesh improvement using swapping and smoothing", *International Journal for Numerical Methods in Engineering* 40 (1997).
- Y. Hu, T. Schneider, B. Wang, D. Zorin and D. Panozzo, "Fast tetrahedral meshing in the wild", *ACM Transactions on Graphics* 39 (2020).
- L. Diazzi, D. Panozzo, A. Vaxman and M. Attene, "Constrained Delaunay tetrahedrization: a robust and practical approach", *ACM Transactions on Graphics* 42 (2023).
- M. Attene, "Indirect predicates for geometric constructions", *Computer-Aided Design* 126 (2020).
- Q. Zhou and A. Jacobson, "Thingi10K: a dataset of 10,000 3D-printing models", arXiv:1605.04797 (2016).

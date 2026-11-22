# Tetrahedra from a skin — `tetrahedralize`

`tetrahedralize` fills a closed triangle surface with tetrahedra: a [`mesh_cad_surface`](./mesh_cad.md) result, a [boundary layer](./boundary_layer.md)'s cap, any watertight skin. It has two backends. **Netgen** keeps the surface's triangles exactly: every input triangle is a facet of the result and every input point a point of it at its index, so the volume conforms to the skin and to whatever else shares it. The **octree** backend is native and always built: [`octree_mesh`](./octree_mesh.md)'s isosurface stuffing, whose boundary is fitted to the skin within a measured distance rather than kept, with dihedral angles bounded where the stuffing's guarantee holds. `backend="auto"`, the default, takes Netgen when the build has it (`-DMESHIOPLUSPLUS_WITH_NETGEN=ON`) and the octree backend otherwise, with a warning saying so; `backend="netgen"` in a build without it throws `Unsupported` naming the flag. So `tetrahedralize` works in every wheel, CLI and WebAssembly build, and Netgen is the upgrade that keeps the skin. The installed **TetGen** and **Gmsh** programs are two more backends that keep it (`backend="tetgen"`, `"gmsh"`): `tetrahedralize` writes each region's surface, runs the program as a separate process and reads the tetrahedra back ([external meshers](#external-meshers-tetgen-gmsh-and-ftetwild)). The installed **fTetWild** is reached the same way (`backend="ftetwild"`), but fits its boundary within an envelope around the skin, as the octree backend does, rather than keeping it.

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

## External meshers: TetGen, Gmsh and fTetWild

A tetrahedral mesher is a [provider](./cpp_api.md#providers-and-plugins), and besides the default one (Netgen) any number can be registered by name and chosen with `backend="<name>"` (C++ `TetrahedralizeBackend::Named` and `TetrahedralizeOptions::mProvider`; `tetrahedralize_backends()` lists them). Three come with every native build and register themselves when their programs are installed, found on `PATH` or named by `MESHIOPLUSPLUS_TETGEN`, `MESHIOPLUSPLUS_GMSH` and `MESHIOPLUSPLUS_FTETWILD` (`MESHIOPLUSPLUS_NO_EXTERNAL_MESHERS=1` turns them all off):

- **`tetgen`** writes the region's surface as a TetGen [`.smesh`](./formats/tetgen.md) (its inner shells as holes, a point inside each) and runs `tetgen -pYQ`: a constrained tetrahedralisation that adds no point on the boundary (`-Y`), with quality (`-q`) when `optimize_steps` is positive and a volume bound (`-a`) from `max_size`, the volume of a regular tetrahedron of that edge. It reads back `.1.node` and `.1.ele`.
- **`gmsh`** writes the surface as a full-precision ASCII STL, one `solid` per shell, with a `.geo` that makes each shell a surface loop of one volume, and runs `gmsh -3` (`-clmax` from `max_size`; `optimize_steps` 0 turns its optimisation off). It reads back the MSH 4.1 output's tetrahedra.
- **`ftetwild`** writes the surface as a full-precision OBJ (its inner shells facing in, so fTetWild's winding-number filter makes them holes) and runs fTetWild's `FloatTetwild_bin -i in.obj -o out.msh --level 6 --no-color` (Y. Hu, T. Schneider, B. Wang, D. Zorin and D. Panozzo, *Fast Tetrahedral Meshing in the Wild*, ACM TOG 39(4), 2020), with `-a` (the ideal edge length) from `max_size` and `--max-its 0` when `optimize_steps` is 0; otherwise fTetWild's own defaults hold (an ideal edge of 0.05 of the bounding box diagonal, an envelope of 1e-3 of it, its own optimisation schedule). It reads back the MSH 2.2 output's tetrahedra. The program is found as `FloatTetwild_bin`, `fTetWild` or `ftetwild`; build it from [fTetWild](https://github.com/loumalouomega/fTetWild) with CMake (`cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build --target FloatTetwild_bin`).

The core checks what TetGen or Gmsh returns, as it checks Netgen and any provider since: every input point back, matched to the program's copy of it within 1e-9 of the region's diagonal (a text round trip rounds the last digit) and then given the input's exact coordinates, and the faces of exactly one tetrahedron the input's triangles. A program that moves the boundary is refused naming it, so `boundary_preserved` is true and `hausdorff_to_surface` zero whenever the call returns. fTetWild is the exception by design: it does not keep the skin's triangles, and its provider says so (`TetMesherProvider::mKeepsBoundary` false), so its result is taken like the octree backend's: every point is new, the triangle block is the tetrahedra's boundary facets (each given the `face_array` value of the nearest skin triangle), `boundary_preserved` is false and `hausdorff_to_surface` is the sampled two-sided distance between that boundary and the skin; a surface of several regions gives each region its own mesh, which do not share points across their interfaces (a warning says so). fTetWild survives dirty input that the others refuse, but `tetrahedralize` still asks for a closed, oriented skin first. `grading` is Netgen's alone. The programs run in a fresh temporary directory, removed afterwards (`MESHIOPLUSPLUS_KEEP_TEMP=1` keeps it), with no shell, their output in a log whose last lines a failure quotes. They never run in WebAssembly, where the names are not registered.

Running them as programs, never linking them, is what their licences leave open: TetGen is AGPL-3.0 and Gmsh GPL-2.0-or-later, and meshio++ neither ships nor links either (the [decision](./roadmap.md#non-goals-and-decisions-taken) on them stands for linking). fTetWild is MPL-2.0, which would allow linking, but it runs as a program too, so its Geogram, libigl and TBB stay out of meshio++'s build. On Ubuntu `apt install tetgen gmsh` provides the first two; CI's `netgen` job runs the tests against them and against fTetWild built from source.

## Octree: the skin fitted

The octree backend meshes each region with [`octree_mesh`](./octree_mesh.md) under uniform skin grading at the surface's own size (its mean edge length, or `max_size` when that is smaller): a body-centred cubic lattice stuffed to the skin, its vertices warped onto it and the result fitted, untangled and optimised for `optimize_steps` rounds (F. Labelle and J. R. Shewchuk, *Isosurface stuffing*, 2007). Where the report of `octree_mesh` says its guarantee holds the dihedral angles lie between 10.7 and 164.8 degrees by construction; the boundary is within `hausdorff_to_surface` of the skin, the two-sided distance between them, and `boundary_preserved` is false. The inside grows to `max_size` or an eighth of the region's diagonal; `grading` is Netgen's alone. It is native, deterministic across parallel backends and thread counts, and in every build.

## Measured

On the unit sphere meshed by `mesh_cad_surface` at `h = 0.2` and 25 degrees, `tetrahedralize` makes 3,145 tetrahedra whose boundary is the skin itself, with dihedral angles of 8.5 to 164.7 degrees (3,680 and 8.9 to 164.4 with `max_size=0.2`); their volume, 4.1500, is the skin's enclosed volume to rounding, 0.9% under `4/3 π`, the skin's own chordal error. A box inside a box gives two regions that share the inner box's triangles. The octree backend fills the same sphere with 15,480 tetrahedra whose dihedral angles lie between 21.4 and 134.3 degrees, its boundary within 0.0097 (0.05 `h`) of the skin and its volume, 4.1252, 0.6% under the skin's; the unit cube's skin it fills exactly, with angles of 45 to 98 degrees. The two trade as their names say: Netgen keeps the skin with a fifth of the tetrahedra and wider angle extremes, the octree backend bounds the angles and moves the boundary. `tests/cpp/test_tetrahedralize.cpp` and `tests/python/test_tetrahedralize.py` pin these, and CI's `netgen` job runs the Netgen tests against Ubuntu's Netgen. [Backend comparison](#backend-comparison) measures all the backends on four skins, a scanned bunny among them.

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

## Backend comparison

[`72_tetrahedralize_backends.ipynb`](https://github.com/loumalouomega/meshioplusplus/blob/main/example/python/72_tetrahedralize_backends.ipynb) meshes four skins with every backend, the native octree backend plain and with `post_adapt=True` ([`adapt_volume`](./adapt_volume.md) on its tetrahedra with every `adapt_*` keyword at its default, v20.1.0; the same `octree_mesh` call as the plain backend, so the two differ only by the remeshing pass), the octree backend again at a count matched to the others (`octree (matched)`, below), Netgen, TetGen, Gmsh and fTetWild, and measures each result with [`compute_quality`](./mesh_quality.md) and the report (November 2026, Ubuntu 24.04 with its packaged Netgen, TetGen 1.5.0 and Gmsh, and fTetWild built from [its fork](https://github.com/loumalouomega/fTetWild) with TBB; `h` is `max_size`, times are the best of three). The skins are the unit sphere (1,280 triangles, `h` 0.2), the unit cube (300, 0.2), a cylinder of radius 0.5 and height 1 (800, 0.2) and the Stanford bunny scan of `example/data`, remeshed to 2,996 triangles by `remesh` and oriented by `repair` (bounding box about 108 units, `h` 5.97), a complex skin with thin ears and concave folds. At equal `h` the plain octree's count is set by the skin's own edge length and comes out two to four times the others', so `octree (matched)` is the same `octree_mesh` call with its skin size and largest leaf both scaled by one factor, found by a few secant steps on the count, until the count is within 15% of the median of Netgen, TetGen and Gmsh on that skin: a factor of 1.44 gives 5,112 tetrahedra on the sphere against a target of 5,193, 1.24 gives 960 on the cube (941), 1.47 gives 1,880 on the cylinder (1,677) and 1.39 gives 11,907 on the bunny (11,566). fTetWild runs at its defaults with `max_size = h` as its ideal edge. "SJ" is the scaled Jacobian, "AR" the aspect ratio, "slivers" the tetrahedra whose smallest dihedral angle is under 15 degrees, and the distance is the report's `hausdorff_to_surface` (`hausdorff_to_skin` for the adapted and matched octree); "kept" says whether the skin's triangles are the result's boundary. No backend produced an inverted tetrahedron.

![Counts, time and slivers by backend and skin](/images/tetrahedralize_backends_counts.png)

| skin | backend | tets | time s | kept | distance | volume error % | dihedral min - max | slivers | SJ < 0.2 | SJ min / mean | AR max / mean |
|---|---|---:|---:|:-:|---:|---:|---|---:|---:|---|---|
| sphere | octree | 15,480 | 0.20 | no | 0.0072 | -1.50 | 31.9 - 121.5 | 0 | 0 | 0.327 / 0.848 | 2.22 / 1.11 |
| sphere | octree (matched) | 5,112 | 0.09 | no | 0.0112 | -2.19 | 35.7 - 117.8 | 0 | 0 | 0.450 / 0.815 | 2.57 / 1.13 |
| sphere | octree + adapt | 13,433 | 2.82 | no | 0.0076 | -1.60 | 36.7 - 131.4 | 0 | 0 | 0.225 / 0.664 | 2.61 / 1.24 |
| sphere | netgen | 4,070 | 0.17 | yes | 0 | -0.86 | 22.3 - 142.8 | 0 | 0 | 0.217 / 0.611 | 3.05 / 1.29 |
| sphere | tetgen | 8,307 | 0.10 | yes | 0 | -0.86 | 9.0 - 164.6 | 86 | 48 | 0.117 / 0.637 | 5.53 / 1.34 |
| sphere | gmsh | 5,193 | 0.13 | yes | 0 | -0.86 | 13.6 - 155.5 | 10 | 0 | 0.205 / 0.598 | 3.33 / 1.34 |
| sphere | ftetwild | 5,671 | 0.47 | no | 0.0020 | -0.93 | 25.1 - 139.5 | 0 | 0 | 0.310 / 0.712 | 2.09 / 1.21 |
| cube | octree | 1,800 | 0.04 | no | 0 | 0 | 45.0 - 90.0 | 0 | 0 | 0.577 / 0.821 | 1.43 / 1.18 |
| cube | octree (matched) | 960 | 0.03 | no | 0 | 0 | 45.0 - 90.8 | 0 | 0 | 0.575 / 0.794 | 1.43 / 1.20 |
| cube | octree + adapt | 1,800 | 0.07 | no | 0 | 0 | 45.0 - 100.5 | 0 | 0 | 0.577 / 0.804 | 1.53 / 1.17 |
| cube | netgen | 941 | 0.09 | yes | 0 | 0 | 27.0 - 139.2 | 0 | 0 | 0.394 / 0.657 | 1.97 / 1.22 |
| cube | tetgen | 1,436 | 0.02 | yes | 0 | 0 | 9.2 - 164.5 | 26 | 16 | 0.118 / 0.578 | 5.80 / 1.45 |
| cube | gmsh | 641 | 0.07 | yes | 0 | 0 | 9.7 - 159.8 | 3 | 11 | 0.099 / 0.548 | 4.25 / 1.49 |
| cube | ftetwild | 2,489 | 1.53 | no | 0.0007 | 0 | 19.9 - 146.2 | 0 | 0 | 0.223 / 0.588 | 3.89 / 1.34 |
| cylinder | octree | 3,744 | 0.07 | no | 0.0054 | -0.88 | 27.1 - 127.8 | 0 | 0 | 0.360 / 0.819 | 1.85 / 1.13 |
| cylinder | octree (matched) | 1,880 | 0.05 | no | 0.0133 | -1.32 | 25.1 - 132.6 | 0 | 0 | 0.250 / 0.662 | 2.25 / 1.28 |
| cylinder | octree + adapt | 3,602 | 0.35 | no | 0.0054 | -0.93 | 38.3 - 134.9 | 0 | 0 | 0.326 / 0.717 | 2.36 / 1.21 |
| cylinder | netgen | 1,677 | 0.14 | yes | 0 | -0.42 | 5.9 - 170.2 | 36 | 126 | 0.065 / 0.449 | 45.99 / 1.89 |
| cylinder | tetgen | 1,966 | 0.04 | yes | 0 | -0.42 | 2.0 - 174.2 | 128 | 196 | 0.019 / 0.432 | 95.48 / 2.22 |
| cylinder | gmsh | 1,503 | 0.09 | yes | 0 | -0.42 | 1.2 - 175.3 | 112 | 219 | 0.027 / 0.355 | 150.76 / 2.66 |
| cylinder | ftetwild | 2,168 | 2.21 | no | 0.0009 | -0.47 | 11.7 - 142.5 | 7 | 5 | 0.127 / 0.595 | 4.69 / 1.33 |
| bunny | octree | 30,786 | 0.64 | no | 2.5368 | -0.50 | 0.2 - 179.7 | 91 | 147 | 0.002 / 0.815 | 50,405 / 3.94 |
| bunny | octree (matched) | 11,907 | 0.26 | no | 3.5069 | -1.12 | 0.4 - 179.0 | 143 | 243 | 0.002 / 0.755 | 6,225 / 2.34 |
| bunny | octree + adapt | 29,061 | 7.72 | no | 2.5368 | -0.59 | 8.6 - 176.7 | 28 | 92 | 0.002 / 0.636 | 444.34 / 1.33 |
| bunny | netgen | 11,566 | 1.14 | yes | 0 | 0 | 0.0 - 179.9 | 775 | 1,322 | 0.000 / 0.469 | 973,746 / 119.6 |
| bunny | tetgen | 17,330 | 0.15 | yes | 0 | 0 | 0.0 - 179.1 | 586 | 803 | 0.000 / 0.595 | 75,532 / 9.59 |
| bunny | gmsh | 9,335 | 0.28 | yes | 0 | 0 | 0.0 - 177.3 | 220 | 370 | 0.000 / 0.527 | 79,834 / 12.32 |
| bunny | ftetwild | 18,818 | 8.84 | no | 0.1355 | -0.09 | 11.0 - 155.7 | 14 | 22 | 0.110 / 0.704 | 7.24 / 1.20 |

![Scaled Jacobian and dihedral-angle distributions](/images/tetrahedralize_backends_quality.png)

What the numbers say:

- **Skin kept or fitted.** The three backends that keep the skin conform to it exactly (distance 0, the same volume error, the skin's own chordal error); the octree backends trade that for the best element shapes on the primitives: plain octree has no sliver and a scaled Jacobian of at least 0.33 on the sphere, cube and cylinder, and on the cube it reproduces the skin exactly with angles of 45 to 90 degrees. Of the three that keep the skin Netgen has the best shapes on the sphere, cube and cylinder (its optimisation passes), TetGen's constrained tetrahedralisation is the fastest and leaves the most slivers on the primitives, Gmsh sits between them.
- **At a comparable count** the octree still has the best shapes on the primitives: matched to Gmsh's count on the sphere (5,112 against 5,193) its scaled Jacobian floor rises to 0.45 and its dihedral angles stay within 35.7 and 117.8 degrees, against 0.21 and 13.6 to 155.5 for Gmsh, and on the cube it matches Netgen's count (960 against 941) with angles of 45 to 91 degrees against Netgen's 27 to 139. What the coarser lattice costs is the boundary: the distance to the skin grows by 1.4 to 2.5 times (0.0072 to 0.0112 on the sphere, 2.54 to 3.51 on the bunny) and the volume error with it (-2.2 percent on the sphere, -1.1 on the bunny), and on the bunny the slivers rise from 91 to 143, since the lattice resolves the ears more coarsely. The plain octree's extra elements buy a closer boundary, not better shapes.
- **fTetWild** fits the skin within its envelope (1e-3 of the bounding box diagonal by default): the distance is 0.002 on the sphere, under 0.001 on the cube and cylinder and 0.14 on the bunny (0.13 percent of its size, against the octree's 2.54), and the volume within 0.1 percent of the skin's on the bunny. It has the best worst elements of every backend on the bunny, by far: smallest dihedral angle 11 degrees, 14 slivers and a worst aspect ratio of 7 where the meshers that keep the skin reach 0 degrees, 220 to 775 slivers and aspect ratios of 75,000 to 970,000, since it is free to move the boundary off the scan's creases. On the sphere and the cylinder its worst shapes sit between the octree's and Netgen's (scaled Jacobian floor 0.31 and 0.13, against Netgen's 0.22 and 0.07), on the cube below Netgen's (0.22 against 0.39), with more tetrahedra than Netgen (5,671 against 4,070 on the sphere, 2,489 against 941 on the cube), and it is the slowest backend on the primitives (0.5 to 2.2 s) and on the bunny (8.8 s), where its envelope checks dominate. Its counts move by about one percent from run to run, since its parallel passes are not deterministic.
- **Count against quality.** Fewer elements is not better by itself: Gmsh and Netgen use 30 to 60 percent fewer tetrahedra than TetGen on the sphere and cube, with fewer bad ones. The plain octree's interior is the finest at equal `h` (3.8 times Netgen's count on the sphere) and so the most expensive for a solver.
- **The adapted octree** never has more elements than the plain one on the same lattice (0.87 of the plain count on the sphere, 1.00 on the cube, whose lattice is already minimal, 0.96 on the cylinder, 0.94 on the bunny) and leaves the boundary where it was (the distance to the skin is unchanged to the fourth digit). What it buys is shape at the extremes: the smallest dihedral angle rises from 27.1 to 38.3 degrees on the cylinder, 31.9 to 36.7 on the sphere and 0.2 to 8.6 on the bunny, whose slivers fall from 91 to 28 and worst aspect ratio from 50,405 to 444. The mean scaled Jacobian falls (0.85 to 0.66 on the sphere), since the lattice's regular cells are the best a mesh can have and any remeshing moves away from them, and the call costs 2 to 13 times as long (3.2 s against 0.22 s on the sphere). It coarsens only a little because the plain octree already sits at the target size; a coarser target is what reduces the count, and `adapt_boundary_sizing="geometry"` can coarsen a flat or gently curved boundary further (0.45 of the plain count on the sphere in a trial) but refined the cylinder's (1.47), so it is not a default.
- **The complex skin** separates the backends. The scan's thin ears and creases leave slivers everywhere: with the skin kept every mesher makes flat tetrahedra (dihedral angles up to 177-180 degrees, 220 to 775 slivers, Netgen the most); the plain octree has 91 slivers and a scaled Jacobian floor of 0.002, but moves the boundary by 2.54 units (2.4 percent of the size) and loses 0.5 percent of the volume. A skin like this is better prepared first ([`repair`](./repair.md), `remesh`, a coarser `h` near the ears) when the skin has to be kept.

![Cut-away views of the sphere's and the bunny's tetrahedra, coloured by scaled Jacobian](/images/tetrahedralize_backends_cutaway.png)

The figures are the notebook's matplotlib fallback, since the machine that executed it had no off-screen GL; the notebook uses PyVista renders where GL is available. TetGen, Gmsh and fTetWild appear only when their programs are installed and Netgen when built in or added, so a run elsewhere has the columns it can. The timings are a single machine's wall clock and the sizes come from one run each, so read the ordering, not the third digit.

## Options and report

| option | meaning |
|---|---|
| `backend` | `auto` (Netgen when built in, else the octree backend, with a warning), `netgen`, `octree`, or a registered mesher's name (`tetgen`, `gmsh`, `ftetwild`; `tetrahedralize_backends()`) |
| `max_size` | the largest edge asked of the tetrahedra (0: none) |
| `max_volume` | the largest volume asked of a tetrahedron (0: none) |
| `grading` | how fast the size grows away from the surface, in `[0, 1]` (0.3; Netgen) |
| `optimize_steps` | passes of the backend's volume optimisation, 0 to 100 (3) |
| `max_cells` | refuse a result of more tetrahedra (50,000,000) |
| `face_array` | the cell data naming each triangle's CAD face (`cad:face`) |

The report has `num_tetrahedra`, `num_points` (every point), `num_new_points` (of those, not the input's: every point with the octree and fTetWild backends), `num_regions`, `num_boundary_triangles`, `boundary_preserved`, `hausdorff_to_surface`, `min_dihedral_deg` and `max_dihedral_deg`, `volume` (the tetrahedra's total), `backend` (the one that ran) and `warnings`.

## Other backends

Geogram (BSD-3) was the first candidate and was rejected: its `mesh_tetrahedralize`, the constrained tetrahedralisation of a closed surface, creates its Delaunay kernel as `Delaunay::create(3, "tetgen")`, and TetGen is AGPL-3.0, which this repository excludes from what it links, with CGAL's `Mesh_3` and Gmsh (GPL); TetGen and Gmsh are reachable as [external meshers](#external-meshers-tetgen-gmsh-and-ftetwild) instead, and so is fTetWild (MPL-2.0; Hu et al. 2020), the candidate recorded here before, which meshes within an envelope around the surface and reports `boundary_preserved` false and a Hausdorff distance. One candidate stays recorded for another backend under the same contract (a named `TetMesherProvider`):
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

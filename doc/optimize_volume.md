# ODT remeshing (`optimize_volume`)

`optimize_volume` raises a tetrahedral mesh's worst element quality by *ODT remeshing* — relocating vertices **and** changing connectivity. By default (`objective="worst"`) it raises the least dihedral angle and never lowers it; `objective="energy"` is the ODT energy with scaled-Jacobian flips it ran alone up to v18.0.0 (see [two objectives](#two-objectives)). It is the genuine "ODT remeshing" the roadmap's volumetric bullet named, and the missing third member of a trio whose other two each do only half the job:

- [`smooth`](smooth.md#odt-smoothing) with `method="odt"` is ODT *smoothing*: it moves each free interior tet vertex to the volume-weighted circumcenter average of its incident tets, on the mesh's **fixed** connectivity. It cannot fix a badly *connected* tetrahedralization — a sliver that no vertex motion removes.
- [`remesh_volume`](remesh_volume.md) *generates* a fresh tet mesh from a signed-distance lattice (surface-in, volume-out); it discards the input's tets entirely rather than improving them.

`optimize_volume` alternates the ODT vertex relocation above with quality-improving topological **flips** (2-3 and 3-2), so both the vertex positions and the connectivity change — which is what makes it *remeshing* rather than *smoothing*, and the resolution-preserving quality optimiser [`refine`](refine.md)/[`decimate_volume`](decimate_volume.md) are not (`refine` subdivides the input's own possibly-bad cells; `decimate_volume` can only remove them).

## Predicate-free by design

The roadmap deliberately rejects the literal Delaunay/ODT-remeshing method because a robust 3D Delaunay kernel needs in-sphere/orientation predicates — exactly where this project's dependency-free posture stops paying (Geogram is a large dependency; TetGen and CGAL's `Mesh_3` are AGPL/GPL). `optimize_volume` needs **none** of that. A flip is applied only when, using the *signed volume* of the candidate tets alone (`detail::cell_volume_from_corners`/`detail::det3`, never an in-sphere test):

1. every new tet is non-degenerate and the local configuration is convex (a pure signed-volume test), and
2. the **minimum** quality (scaled Jacobian) over the new tets strictly exceeds the minimum over the tets it replaces, by `min_improvement`.

Criterion 2 is Freitag & Ollivier-Gooch's local-mesh-improvement-by-swapping rule ("Tetrahedral mesh improvement using swapping and smoothing", 1997). Because every accepted flip strictly raises a bounded quantity — the worst incident quality — the process is monotone in worst quality and therefore terminates, with no Delaunay optimality argument and no predicate.

A **2-3 flip** replaces the two tets sharing an interior triangular face with three tets around the edge joining their two apexes; a **3-2 flip** is its inverse (three tets around an interior edge become two capping the ring). The relocation half is delegated verbatim to `smooth`'s `method="odt"`, reusing all of its boundary/feature/frozen pinning and inversion guard.

![A 2-3 flip replaces two tetrahedra sharing a face with three around the new edge between their apexes, and the 3-2 flip is its inverse](/diagrams/flips_23_32.svg)

## Two objectives

The rule above, with the scaled Jacobian, lowers an energy summed over the cells; an ODT relocation moves a vertex to where its star's energy is lower even when one of its tets gets worse. Measured on Netgen's tetrahedra, that **lowers the least dihedral angle** on every skin once Netgen's own passes have run (the `energy` column below). It is kept as `objective="energy"` (C++ `OptimizeVolumeObjective::Energy`), byte for byte what the operation did up to v18.0.0.

The default, `objective="worst"`, targets the worst cell instead (Freitag and Ollivier-Gooch 1997; Klingner and Shewchuk 2007). The quality is the least **dihedral angle** of a tet, and every step is accepted only when it raises the least angle of the tets it touches, so the least angle of the mesh never drops. Each sweep runs:

1. **Guarded ODT relocation.** `smooth`'s ODT pass, then the guard: a vertex whose star's least angle fell below what it was before the pass goes back, in ascending vertex order, repeatedly, until no moved vertex's star is worse than before. Every tet that changed has a moved vertex in it, so none ends below the least angle the pass started from.
2. **Topological moves**, each over interior faces and edges only:
   - the **2-3 flip**, by the least dihedral angle;
   - **edge removal**: the *k* tets around an interior edge, whose far vertices form a closed ring, are replaced by the 2 (*k* − 2) tets over the triangulation of the ring that maximises their least angle, found by the dynamic program over polygon triangulations; *k* = 3 is the 3-2 flip, and `max_ring_size` (default 7) bounds *k*;
   - **multi-face removal**, its inverse: the faces sandwiched between two vertices *a* and *b* — each with the tets (face, *a*) and (face, *b*) — grown from an interior face into a disc of at most `max_ring_size` boundary vertices, become the tets around the new edge *ab*.

   A move is considered only where the worst tet it replaces is below 40 degrees, and applied only when every new tet is positively oriented, the new tets fill exactly the old ones' volume, and their least angle beats the old least angle by `min_improvement` degrees.
3. **Worst-cell smoothing.** The free interior vertices of tets below 40 degrees, worst first, each climb the least dihedral angle of their star: the search direction is the shortest vector in the convex hull of the gradients of the angles within a small band of the least (the active set), and a backtracking line search accepts a step only when the star's least angle rises. Boundary and frozen vertices stay.

The sweeps stop at a fixed point or after `max_iterations`. The report adds `num_edge_removals` (rings of four or more; three-tet rings count as `num_32_flips`), `num_multi_face_removals`, `num_worst_smoothed`, and `min_dihedral_before`/`min_dihedral_after` (both objectives report them); under `worst`, `min_dihedral_after >= min_dihedral_before` always.

### Measured

On Netgen's tetrahedra, from the `netgen-mesher` Python package's volume mesher (`maxh = 0.2`, with and without its three optimisation passes) filling two of meshio++'s own skins — `icosa_sphere(1, 3)` and the triangulated skin of a wedge `cylinder(radius=1, height=2)` — the least dihedral angle in and out, and the tets under 15 degrees (`tools/optvol_probe.py`, October 2026):

| skin | Netgen passes | tets | least angle in | `energy` | `worst` | under 15° in / `energy` / `worst` |
|---|---:|---:|---:|---:|---:|---|
| sphere | 0 | 4,990 | 0.97 | 0.11 | 31.23 | 372 / 52 / 0 |
| sphere | 3 | 3,703 | 17.37 | 12.49 | 32.35 | 0 / 1 / 0 |
| capped cylinder | 0 | 15,344 | 1.63 | 0.03 | 20.93 | 300 / 115 / 0 |
| capped cylinder | 3 | 11,439 | 13.30 | 8.32 | 19.55 | 2 / 56 / 0 |

The `ball` primitive's tetrahedra (`ball(n=6, cell_type="tetra")`, 10,368 of them, flattened into slivers along the sphere) go from a least angle of 0.08 degrees and 372 tets under 15 to 30.2 degrees and none under `worst`, and to 2.6 degrees and 60 under `energy` (`example/python/50_optimize_volume.ipynb`):

![The ball after optimize_volume's worst objective, cut open and coloured by each tetrahedron's least dihedral angle, every one above 30 degrees](/images/optimize_volume_worst.png)

On the cylinder the worst tets out have two to four of their vertices on the boundary, which smoothing does not move; on the sphere it passes 30 degrees, the figure Klingner and Shewchuk report. `tests/python/test_optimize_volume_netgen.py` pins the guarantee on Netgen's output (it skips without `netgen-mesher`), and `tests/cpp/test_optimize_volume.cpp` on fixtures the core builds alone.

## Boundary is invariant

A 2-3 flip acts only on an *interior* face (shared by two tets) and a 3-2 flip, an edge removal or a multi-face removal only on *interior* edges and faces, so none ever touches a boundary face; worst-cell smoothing never moves a boundary vertex. Combined with `preserve_boundary` pinning boundary vertices during relocation, the output's boundary surface is **byte-identical** to the input's: watertight in ⇒ watertight out, with none of the coincident-edge risk `remesh_volume`'s surface warp carries.

## What survives

The **point set is invariant** — relocation moves points, flips only reconnect them, and no point is ever added or removed — so `point_data`, `field_data` and named **Point** regions carry through unchanged. `cell_data` has no correspondence across a flip (a 2-3 flip replaces two cells with three) and is **dropped with a warning**, as are named **Cell**/**Side** regions; the output is a single `tetra` block. This is more generous than `remesh_volume` (which drops point data too, having a genuinely new point set), and it is honest: what is preserved is exactly what a preserved point set can preserve. To keep `cell_data` or a field defined on the boundary, compose with [`interpolate`](interpolate.md) or [`conservative_interpolate`](conservative_interpolate.md).

## Scope and determinism

Tet-only: a non-`tetra` 3D block raises pointing at `convert_cells(mode="simplexify")`, a non-3D block alongside the tets raises pointing at `split`, and a ragged/polyhedron or tet-less mesh raises by name — the same tet-only scope `smooth`'s ODT method and `decimate_volume` enforce. Per-tet quality and each sweep's face/edge adjacency are built in parallel into disjoint slots; the relocation guard, the move-application loops and the worst-cell smoothing are serial in a fixed order, so output is byte-identical across the three mesh backends and across thread counts.

Like [`remesh_volume`](remesh_volume.md)/[`subdivide`](subdivide.md)/[`agglomerate`](agglomerate.md)/[`decimate_volume`](decimate_volume.md), this operation is **C++-core only, with no numpy fallback** — the flip acceptance is a discrete branch on the sign of a volume and on a near-tie quality comparison, so a second independent implementation could land on the other side of such a tie and diverge into a different connectivity. The Python function raises `NotImplementedError` by name when the compiled core is unavailable.

## Usage

::: code-group

```python [Python]
import meshioplusplus as mio

mesh = mio.read("volume.vtu")          # a tetrahedral mesh
out, report = mio.optimize_volume(mesh, max_iterations=10, return_report=True)
print(report["num_flips"], report["min_dihedral_before"], report["min_dihedral_after"])
old = mio.optimize_volume(mesh, objective="energy")   # the ODT energy, as up to v18.0.0
mio.write("optimized.vtu", out)
```

```bash [CLI]
meshioplusplus optimize-volume volume.vtu optimized.vtu
# --max-iterations N --no-relocate --no-flip --no-preserve-boundary --min-improvement E
# --objective worst|energy --max-ring-size K
```

:::

Bindings mirror the operation on every surface: Python `optimize_volume`, C `mio_optimize_volume` and `mio_optimize_volume_ex` (`mio_optimize_volume_opts` with `objective` and `max_ring_size`, `mio_optimize_volume_report` with every counter; the flat entry point keeps its signature and runs the default objective), Fortran `m%optimize_volume`, Julia `optimize_volume`, R `mio_optimize_volume`, WASM `optimizeVolume`, the `optimize-volume` CLI verb in both CLIs, an `OptimizeVolume` settings-pipeline step, and an `optimize_volume` MCP tool. The relocation pin mask (`frozen`) is exposed on the Python and C++ APIs only — a documented flat-ABI gap of its own, distinct from `smooth`'s (whose `frozen` reaches every flat binding: C `mio_smooth_ex`, Fortran, Julia, R, and WASM).

## Attribution

The flip-acceptance rule (apply a local transformation only when it strictly improves the worst incident element quality) and the non-smooth worst-angle smoothing are the published method of Freitag & Ollivier-Gooch, "Tetrahedral mesh improvement using swapping and smoothing" (*Int. J. Numer. Methods Eng.*, 1997); edge removal by the dynamic program over ring triangulations and multi-face removal are those of Shewchuk, "Two discrete optimization algorithms for the topological improvement of tetrahedral meshes" (2002), and Klingner & Shewchuk, "Aggressive tetrahedral mesh improvement" (*16th International Meshing Roundtable*, 2007) — all implemented from the descriptions only. No external mesh library, and in particular no in-sphere/orientation predicate kernel, is read or vendored. The ODT vertex relocation is this repo's own `smooth` `method="odt"` (Alliez et al. 2005; Chen & Xu).

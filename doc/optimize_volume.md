# ODT remeshing (`optimize_volume`)

`optimize_volume` raises a tetrahedral mesh's worst element quality by *ODT remeshing* — relocating vertices **and** changing connectivity. It is the genuine "ODT remeshing" the roadmap's volumetric bullet named, and the missing third member of a trio whose other two each do only half the job:

- [`smooth`](smooth.md#odt-smoothing) with `method="odt"` is ODT *smoothing*: it moves each free interior tet vertex to the volume-weighted circumcenter average of its incident tets, on the mesh's **fixed** connectivity. It cannot fix a badly *connected* tetrahedralization — a sliver that no vertex motion removes.
- [`remesh_volume`](remesh_volume.md) *generates* a fresh tet mesh from a signed-distance lattice (surface-in, volume-out); it discards the input's tets entirely rather than improving them.

`optimize_volume` alternates the ODT vertex relocation above with quality-improving topological **flips** (2-3 and 3-2), so both the vertex positions and the connectivity change — which is what makes it *remeshing* rather than *smoothing*, and the resolution-preserving quality optimiser [`refine`](refine.md)/[`decimate_volume`](decimate_volume.md) are not (`refine` subdivides the input's own possibly-bad cells; `decimate_volume` can only remove them).

## Predicate-free by design

The roadmap deliberately rejects the literal Delaunay/ODT-remeshing method because a robust 3D Delaunay kernel needs in-sphere/orientation predicates — exactly where this project's dependency-free posture stops paying (Geogram is a large dependency; TetGen and CGAL's `Mesh_3` are AGPL/GPL). `optimize_volume` needs **none** of that. A flip is applied only when, using the *signed volume* of the candidate tets alone (`detail::cell_volume_from_corners`/`detail::det3`, never an in-sphere test):

1. every new tet is non-degenerate and the local configuration is convex (a pure signed-volume test), and
2. the **minimum** quality (scaled Jacobian) over the new tets strictly exceeds the minimum over the tets it replaces, by `min_improvement`.

Criterion 2 is Freitag & Ollivier-Gooch's local-mesh-improvement-by-swapping rule ("Tetrahedral mesh improvement using swapping and smoothing", 1997). Because every accepted flip strictly raises a bounded quantity — the worst incident quality — the process is monotone in worst quality and therefore terminates, with no Delaunay optimality argument and no predicate.

A **2-3 flip** replaces the two tets sharing an interior triangular face with three tets around the edge joining their two apexes; a **3-2 flip** is its inverse (three tets around an interior edge become two capping the ring). The relocation half is delegated verbatim to `smooth`'s `method="odt"`, reusing all of its boundary/feature/frozen pinning and inversion guard.

## Boundary is invariant

A 2-3 flip acts only on an *interior* face (shared by two tets) and a 3-2 flip only on an *interior* edge (all its incident faces interior), so neither ever touches a boundary face. Combined with `preserve_boundary` pinning boundary vertices during relocation, the output's boundary surface is **byte-identical** to the input's: watertight in ⇒ watertight out, with none of the coincident-edge risk `remesh_volume`'s surface warp carries.

## What survives

The **point set is invariant** — relocation moves points, flips only reconnect them, and no point is ever added or removed — so `point_data`, `field_data` and named **Point** regions carry through unchanged. `cell_data` has no correspondence across a flip (a 2-3 flip replaces two cells with three) and is **dropped with a warning**, as are named **Cell**/**Side** regions; the output is a single `tetra` block. This is more generous than `remesh_volume` (which drops point data too, having a genuinely new point set), and it is honest: what is preserved is exactly what a preserved point set can preserve. To keep `cell_data` or a field defined on the boundary, compose with [`interpolate`](interpolate.md) or [`conservative_interpolate`](conservative_interpolate.md).

## Scope and determinism

Tet-only: a non-`tetra` 3D block raises pointing at `convert_cells(mode="simplexify")`, a non-3D block alongside the tets raises pointing at `split`, and a ragged/polyhedron or tet-less mesh raises by name — the same tet-only scope `smooth`'s ODT method and `decimate_volume` enforce. Per-tet quality and each sweep's face/edge adjacency are built in parallel into disjoint slots; the flip-application loop is serial in a fixed order, so output is byte-identical across the three mesh backends and across thread counts.

Like [`remesh_volume`](remesh_volume.md)/[`subdivide`](subdivide.md)/[`agglomerate`](agglomerate.md)/[`decimate_volume`](decimate_volume.md), this operation is **C++-core only, with no numpy fallback** — the flip acceptance is a discrete branch on the sign of a volume and on a near-tie quality comparison, so a second independent implementation could land on the other side of such a tie and diverge into a different connectivity. The Python function raises `NotImplementedError` by name when the compiled core is unavailable.

## Usage

::: code-group

```python [Python]
import meshioplusplus as mio

mesh = mio.read("volume.vtu")          # a tetrahedral mesh
out, report = mio.optimize_volume(mesh, max_iterations=10, return_report=True)
print(report["num_flips"], report["min_quality_before"], report["min_quality_after"])
mio.write("optimized.vtu", out)
```

```bash [CLI]
meshioplusplus optimize-volume volume.vtu optimized.vtu
# --max-iterations N --no-relocate --no-flip --no-preserve-boundary --min-improvement E
```

:::

Bindings mirror the operation on every surface: Python `optimize_volume`, C `mio_optimize_volume`, Fortran `m%optimize_volume`, Julia `optimize_volume`, R `mio_optimize_volume`, WASM `optimizeVolume`, the `optimize-volume` CLI verb in both CLIs, an `OptimizeVolume` settings-pipeline step, and an `optimize_volume` MCP tool. The relocation pin mask (`frozen`) is exposed on the Python and C++ APIs only — a documented flat-ABI gap shared with `smooth`.

## Attribution

The flip-acceptance rule (apply a local transformation only when it strictly improves the worst incident element quality) is the published method of Freitag & Ollivier-Gooch, "Tetrahedral mesh improvement using swapping and smoothing" (*Int. J. Numer. Methods Eng.*, 1997) — implemented from the description only. No external mesh library, and in particular no in-sphere/orientation predicate kernel, is read or vendored. The ODT vertex relocation is this repo's own `smooth` `method="odt"` (Alliez et al. 2005; Chen & Xu).

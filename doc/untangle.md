# Untangling

`untangle(mesh)` repairs inverted tetrahedra and hexahedra by moving nodes until every cell corner has a positive signed volume, keeping the connectivity. Tangled meshes come out of large morphs ([`morph`](./morph.md) with `allow_inverted`), out of explicit-dynamics result files whose deformed configuration folded, and out of hand edits.

```python
import meshioplusplus as mp

mesh = mp.read("deformed.vtu")
fixed = mp.untangle(mesh)                                 # boundary fixed, only nodes near tangles move
fixed = mp.untangle(mesh, fix_boundary=False, rings=3)    # let more nodes move
better = mp.untangle(mesh, smooth=True)                    # untangle and smooth everything

fixed, report = mp.untangle(mesh, return_report=True)
report["num_inverted_before"], report["num_inverted_after"], report["min_quality_after"]
```

## What it minimises

The untangling and smoothing functional of Escobar, Rodríguez, Montenegro, Montero and González-Yuste ("Simultaneous untangling and smoothing of tetrahedral meshes", CMAME 192, 2003). Every corner simplex `c`, a tetrahedron itself or one of a hexahedron's eight corners with its three edges (Knupp's corner Jacobians), contributes

```
|S|² / (3 h(det S)^(2/3)),      h(s) = (s + sqrt(s² + 4 δ²)) / 2
```

with `S` the corner's edge matrix against the ideal one (a regular tetrahedron, a cube corner), in units of the mesh's mean edge length. `h` is a regularised signed volume: positive even when `det S` is negative, so the functional is finite across an inversion and falls as the corner turns right side out; `δ` comes from the most negative corner of each sweep (Escobar et al.'s rule, with `ε = 1e-3`) and is 0 once none is negative, where the term is the inverse mean ratio, 1 for an ideal corner. A sweep visits the free nodes in ascending order (Gauss–Seidel), each taking up to four steps of steepest descent from a central-difference gradient and a backtracking line search that accepts only a decrease; it is serial, so the result is the same on every backend.

**What moves.** By default only the nodes within `rings` edge rings (2) of an inverted cell, recomputed each sweep, and it stops as soon as no cell is inverted: a repair, not a remesh. `smooth=True` moves every free node until the functional settles (a relative change below 1e-6 per sweep), smoothing as it untangles. The boundary nodes (of faces used by one cell) stay fixed with `fix_boundary` (the default), as do the `frozen` nodes and the nodes of any volume cell it does not measure (`wedge`, `pyramid`, polyhedra, higher order).

"Inverted" here means a corner with a non-positive signed volume, which is stricter for a hexahedron than [`check_quality`](./mesh_quality.md)'s whole-cell volume: a hexahedron can have a positive volume and a folded corner.

## Against `smooth` and `optimize_volume`

Measured, as the roadmap asked, on a tetrahedral cube (10 × 10 × 10 hexahedra split into 6 000 tetrahedra) with a fifth of its interior nodes jittered by 0.08, which inverts 486 tetrahedra:

| Operation | Inverted after | Connectivity | Hexahedra |
| --- | --- | --- | --- |
| [`smooth`](./smooth.md) | 19 (its inversion guard keeps a valid mesh valid; it repairs a tangle only incidentally) | kept | yes |
| [`optimize_volume`](./optimize_volume.md) | 0 | changed (it flips edges) | refused |
| `untangle` | 0, in 2 sweeps, 0.1 s | kept | yes |

On the roadmap's probe, one interior node of a tetrahedral cube pushed through a face (6 tetrahedra inverted), all three end valid; `untangle` takes one sweep, moving only the nodes within two rings of the tangle, and lifts the least corner scaled Jacobian from -0.24 to 0.30.

## Report

`num_inverted_before`, `num_inverted_after`, `min_quality_before` and `min_quality_after` (the least corner scaled Jacobian, the determinant over the edge lengths' product: 1 for a cube corner, 0.71 for a regular tetrahedron, negative when inverted), `sweeps` and `moves`. Python warns when cells are still inverted after `max_sweeps`; free more nodes then (`fix_boundary=False`, more `rings`), since a tangle whose fix needs a boundary node to move cannot be undone with the boundary fixed.

There is no numpy twin (the [`geodesic_distance`](./geodesic.md) precedent).

## Everywhere else

| Surface | Call |
| --- | --- |
| Python | `untangle(mesh, *, fix_boundary=True, frozen=None, rings=2, smooth=False, max_sweeps=100, return_report=False)` |
| Pipeline | `{"Op": "Untangle", "FixBoundary": true, "Rings": 2, "Smooth": false, "MaxSweeps": 100}`; counters `NumInvertedBefore`, `NumInvertedAfter`, `MinQualityAfter`, `Sweeps` |
| Both CLIs | `meshioplusplus untangle IN OUT [--no-fix-boundary] [--rings N] [--smooth] [--max-sweeps N] [--json]` |
| C | `mio_untangle_opts_init(&o); mio_mesh* out = mio_untangle(mesh, &o, &report);` |
| Fortran | `out = m%untangle(rings=3, num_inverted_after=n, stat=ierr)` |
| Julia | `r = untangle(mesh; rings=3)`; `r.mesh`, `r.num_inverted_after`, ... |
| R | `r <- mio_untangle(mesh, rings = 3)`; `r$mesh`, `r$num_inverted_after`, ... |
| WASM | `m.untangle(mesh, {fixBoundary, frozen, rings, smooth, maxSweeps})` |
| MCP | tool `untangle(input_path, output_path, fix_boundary, frozen, rings, smooth, max_sweeps)` |

# Mesh morphing

`morph(mesh, displacement)` moves a mesh's boundary (or any chosen points) by a prescribed displacement and carries the motion smoothly into the rest of the mesh, keeping its connectivity. Uses: shape optimisation, where a design surface moves and the volume mesh around it must follow, and fluid–structure interaction, where the fluid mesh moves with the structure.

```python
import meshioplusplus as mp

mesh = mp.read("duct.vtu")                       # carries point_data["d"] on its moving wall
out = mp.morph(mesh, "d")                        # boundary takes d, the interior follows
out = mp.morph(mesh, d, steps=4)                 # d an (n, 3) array; large motion in 4 increments
out = mp.morph(mesh, "d", fix_boundary=False, constrained="wall")   # only the wall's points pinned

out, report = mp.morph(mesh, "d", return_report=True)
report["num_inverted"], report["max_displacement"], report["converged"]
```

## What it solves

The **constrained** points take the displacement as given: the boundary points (`fix_boundary`, the points of facets used by one cell; give the boundary that stays put a zero displacement), the points `constrained` names (a `point_set`, a nonzero point-data array, or a boolean mask), and any point no surface or volume cell uses. Every other point moves by the **harmonic extension** of those values: per component, `K u = 0` at the free points with `u = d` at the constrained ones, `K` the P1 stiffness of the mesh's top-dimensional cells (hexahedra, wedges, pyramids, quads and polygons decomposed into simplices for the assembly, as [`convert_cells(mode="simplexify")`](./convert_cells.md) does; the mesh itself keeps its cells), solved by Jacobi-preconditioned conjugate gradients.

**Stiffening.** Each simplex's stiffness is scaled by `(mean measure / its measure)^stiffening` (Jacobian-based stiffening; Tezduyar et al. 1992, Stein, Tezduyar and Benney 2003), so small cells, typically the boundary layer next to a moving wall, move nearly rigidly and keep their shape, while large cells far away absorb the deformation. `stiffening=0` is the plain Laplace extension. **Steps.** With `steps` > 1 the displacement is applied in equal increments, the stiffness reassembled on the moved mesh each time: what lets a large motion through.

Rotating the inner circle of an annulus (radii 0.3 and 1.5, 48 × 12 quads split into triangles):

| Rotation | `stiffening` 0 | `stiffening` 1 | `stiffening` 1, `steps` 4 |
| --- | --- | --- | --- |
| 30° | valid | valid | valid |
| 60° | valid | valid | valid |
| 90° | 48 cells inverted | 48 cells inverted | valid |

## The inversion guard

Afterwards every linear volume cell's signed volume and every surface cell's orientation (its Newell normal) is compared with the input's. A cell that turned inside out is counted (`num_inverted`), and unless `allow_inverted=True` the call raises, naming the remedies: more `steps`, a higher `stiffening`. It never returns an inverted mesh silently.

The roadmap's probe, lifting a hexahedral box's top face by 0.1 (a bump, highest at its centre), inverts no cell; a rigid motion of the whole boundary moves every point rigidly (the harmonic extension of a linear field is that field).

## Scope

- Linear cells: `tetra`, `hexahedron`, `wedge`, `pyramid`, `triangle`, `quad`, `polygon`; `line` and `vertex` blocks ride along. A polyhedron is refused naming `poly_to_tet`, a higher-order cell naming `linearize`. The extension runs over the highest-dimensional cells (a volume mesh's boundary faces, if present as cells, add no stiffness).
- The displacement is `(n, dim)` or `(n, 3)`; only its constrained rows are read. A 2-D mesh moves in its plane.
- Points move; connectivity, data, regions and the points' dtype-independent values pass through. A constrained point lands exactly on its target, whatever the steps.

Compared with [`sobolev_deform`](./sobolev_deform.md), which filters a displacement the caller has for every point and pins some points at zero, `morph` takes a displacement for some points only and makes up the rest. There is no numpy twin (the [`geodesic_distance`](./geodesic.md) precedent).

## Report

`num_constrained`, `num_free`, `num_inverted`, `iterations` (conjugate-gradient iterations, all solves), `max_displacement` and `converged` (false when a solve hit its cap, with a warning).

## Everywhere else

| Surface | Call |
| --- | --- |
| Python | `morph(mesh, displacement, *, fix_boundary=True, constrained=None, stiffening=1.0, steps=1, allow_inverted=False, tolerance=1e-10, max_iterations=0, return_report=False)` |
| Pipeline | `{"Op": "Morph", "Array": "d", "FixBoundary": true, "ConstrainedArray": "", "Stiffening": 1, "Steps": 1, "AllowInverted": false}`; counters `NumConstrained`, `NumInverted`, `MaxDisplacement`, `Converged` |
| Both CLIs | `meshioplusplus morph IN OUT --array d [--no-fix-boundary] [--constrained-array NAME] [--stiffening S] [--steps N] [--allow-inverted] [--json]` |
| C | `mio_morph_opts_init(&o); o.array = "d"; o.steps = 4; mio_mesh* out = mio_morph(mesh, &o, &report);` |
| Fortran | `out = m%morph("d", steps=4, num_inverted=ninv, stat=ierr)` |
| Julia | `r = morph(mesh, "d"; steps=4)`; `r.mesh`, `r.num_inverted`, ... |
| R | `r <- mio_morph(mesh, "d", steps = 4)`; `r$mesh`, `r$num_inverted`, ... |
| WASM | `m.morph(mesh, {array: 'd', steps: 4, ...})` returns `{mesh, numConstrained, numFree, numInverted, iterations, maxDisplacement, converged}` |
| MCP | tool `morph(input_path, output_path, array, fix_boundary, constrained, stiffening, steps, allow_inverted)` |

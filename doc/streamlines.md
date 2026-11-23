# Streamlines

`streamlines(mesh, "U", seeds=…)` follows a point **vector field** through the mesh and returns the paths as a new mesh of line segments: the points of every line, one `line` cell between consecutive points, and the data a post-processor wants along a path. It is a mesh **operation**, not a file format.

```python
import meshioplusplus as mp

mesh = mp.read("flow.vtu")                                   # carries a 3-component point array "U"
lines = mp.streamlines(mesh, "U", seeds=200)                 # about 200 seeds, both directions
lines = mp.streamlines(mesh, "U", seed_points=[[0.1, 0.5, 0.5]],
                       direction="forward", integrator="rk45")
t = lines.point_data["IntegrationTime"]                      # the field's time since the seed
mp.write("lines.vtp", lines)
```

## What it traces

The field is the **piecewise-linear interpolation** of the point array over the mesh's cells, each split into simplices: triangles and quads (two triangles) for a surface, tetrahedra, hexahedra (six), wedges (three) and pyramids (two) for a volume; a higher-order cell contributes its corner nodes. A volume mesh is traced in its **volume**, a surface mesh on its **surface** (the field is projected onto the triangle it is on and the point stays on it); the dimension of the highest cells decides which. A 2-component array is the field `(u, v, 0)`.

Each line solves `dx/dt = v(x)` from its seed, so **`IntegrationTime` is the field's own time** (and a faster field covers more distance per unit of it), unlike the unit-speed arc length of [`render`](/viewer)'s streamlines, which keep their own fixed tracer and are unchanged.

| `integrator` | Method |
| --- | --- |
| `"rk4"` (default) | Fixed-step classical Runge–Kutta; the step is `step` as a length at the fastest speed (default `0.005` of the mesh's diagonal) |
| `"rk45"` | Dormand–Prince 5(4) with step control: a step is retried smaller until its error estimate is within `tolerance` (a fraction of the diagonal, default `1e-6`); a stage that leaves the mesh counts as a failure |

## The result

| Array | Where | Meaning |
| --- | --- | --- |
| `IntegrationTime` | point data | the time since the seed; **negative** behind it, increasing along every line |
| `Speed` | point data | `|v|` at the point |
| `SeedId` | point data and cell data | the index of the line's seed (explicit seeds in the order given, or the ranked simplex centroids) |
| the field, under its own name | point data | `v` at the point |

`return_report=True` adds `num_lines`, `num_seeds_outside`, `seed_ids`, `dim` and, per line, `end_backward` and `end_forward`: the index of why that end stopped, in `STOP_REASONS`:

| Reason | Meaning |
| --- | --- |
| `left_domain` | the next point is outside the mesh (or, adaptive, no step within the tolerance stays inside) |
| `stagnation` | the speed fell to `terminal_speed` or below (default `1e-9` of the fastest speed) |
| `max_length` | the path reached `max_length`, a fraction of the diagonal (default `0.5`) in that direction |
| `max_steps` | the step limit (default 20000) was used up, as in a very slow region |

## Seeds

Without `seed_points`, `seeds` simplices are ranked at **equal measure along a golden-ratio sequence** (a plain stride would line up with the rows of a structured mesh) and each is started at its centroid: the choice `render`'s streamlines make, so a count gives the same lines in both. Explicit `seed_points` (shape `(n, 3)`, or `(n, 2)` for a planar mesh) are used instead; a seed that is not inside the mesh makes no line and is counted in `num_seeds_outside`. A line of a single point (both ends stopped at once) is dropped.

## Determinism, and what it is not

The lines are traced in parallel into per-seed buffers joined in seed order, so the result does not depend on the parallel backend or the thread count. The locator is a uniform bucket grid, so the cost is linear in the number of simplices and the steps.

- **C++ core only**, with no numpy fallback: the adaptive step control branches on rounded error estimates, which a second implementation could decide differently. Without the compiled `_core` the function raises `NotImplementedError`.
- The piecewise-linear field is **first order**: a streamline of a smooth field on a coarse mesh follows the mesh's interpolation, not the underlying field, and the integrator's tolerance says nothing about that error. A rigid rotation `(-y, x)` is linear, so its circles are exact to the integrator's accuracy; the tests pin them.
- It is a *tracer*, not a particle physics model: no inertia, no diffusion, no time-dependent field (trace each step of a [sequence](/sequences) in turn).

## CLI

```bash
meshioplusplus streamlines flow.vtu lines.vtu --array U --seeds 200
meshioplusplus streamlines flow.vtu lines.vtu --array U --seed-points "0.1,0.5,0.5;0.9,0.5,0.5" \
    --direction forward --integrator rk45 --max-length 1.0 --json
```

See the [CLI reference](/cli#meshioplusplus-streamlines); both CLIs print the same JSON (`output`, `num_lines`, `num_points`, `num_seeds_outside`, `dim`). The MCP server carries it as the `streamlines` tool.

## Other languages

The C++ API is `streamlines(mesh, StreamlinesOptions)` in `operations/streamlines.hpp` (`StreamlinesResult` returns the mesh, the counts and the per-line seed ids and end reasons). It is **not yet in the flat bindings** (C, Fortran, Julia, R, Octave, WebAssembly): see the [roadmap](/roadmap#math-on-meshes).

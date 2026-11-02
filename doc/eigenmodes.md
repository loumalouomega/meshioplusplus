# Laplacian eigenmodes

`laplacian_eigenmodes(surface, num_modes)` computes the lowest eigenpairs of a triangle surface's Laplace–Beltrami operator: the surface's own Fourier basis, smooth functions ordered by how fast they vary across it. Uses: spectral positional encodings for graph networks (each point gets the values of the first modes as features, independent of how the mesh is numbered; see [machine learning](./ml.md)), shape descriptors (the eigenvalues do not change under rigid motion or remeshing), and the first vibration modes of a membrane.

```python
import meshioplusplus as mp

skin = mp.read("bracket_skin.stl")
out = mp.laplacian_eigenmodes(skin, 16)
phi = out.point_data["laplacian:modes"]           # (n, 16), one column per mode
lam = out.field_data["laplacian:eigenvalues"]     # (16,), ascending

out, report = mp.laplacian_eigenmodes(skin, 32, return_report=True)
report["iterations"], report["max_residual"], report["converged"]
```

## What it solves

Over the surface's linear triangles, the generalised eigenproblem `K φ = λ M φ`, with `K` the cotangent stiffness and `M` the lumped mass (a third of each triangle's area at each corner): the operators [`geodesic_distance`](./geodesic.md) uses, from the same code. A boundary carries the natural (Neumann) condition. The modes are `M`-orthonormal (`φᵢᵀ M φⱼ = δᵢⱼ`), so `Σ M φᵢ f` is the `i`-th coefficient of a field `f`, and each is signed so that its largest-magnitude entry is positive. Each connected piece has a zero eigenvalue with a constant mode, so the count of zeros is the count of pieces.

On the unit sphere the eigenvalues are `l (l + 1)`, each `2 l + 1` times; on an icosahedral sphere (the roadmap's probe):

| Points | `λ₁ … λ₃` (2) | `λ₄ … λ₈` (6) | `λ₉ … λ₁₅` (12) | Time, 16 modes |
| --- | --- | --- | --- | --- |
| 162 | 1.99991 | 5.86450 | 11.32 – 11.34 | 0.05 s |
| 642 | 1.99999 | 5.96586 | 11.83 | 0.24 s |
| 2 562 | 2.00000 | 5.99145 | 11.957 – 11.958 | 1.1 s |
| 10 242 | 2.00000 | 5.99786 | 11.989 – 11.990 | 8 s |

The first two clusters stay exactly degenerate (the icosahedron's symmetry keeps them so); the seven-fold `l = 3` cluster splits into four and three, as the icosahedral group requires, and closes under refinement. Times are on 4 cores.

## How

Block inverse iteration with a shift, and Rayleigh–Ritz. A block of `max(2 num_modes, num_modes + 8)` vectors is repeatedly multiplied by `(K + s M)⁻¹ M` (`s` one over the surface's area, which keeps the system positive definite and the zero eigenvalues reachable), each solve Jacobi-preconditioned conjugate gradients warm-started from the previous iterate and only as accurate as the block already is; the block is `M`-orthonormalised and the eigenpairs of `K` within it taken. Because a whole block iterates, a repeated eigenvalue comes out with all its copies, which a single-vector Lanczos iteration would miss. It stops when every wanted pair's residual `‖K φ − λ M φ‖` (in the `M⁻¹` norm, relative to `λ + s`) is below `tolerance` (1e-8, so the eigenvalues are accurate to about 1e-14 relative), or after `max_iterations` (300) with a warning. The convergence rate is the ratio of the last wanted eigenvalue to the first one outside the block, so it takes 20 to 40 iterations whatever the mesh size.

Everything is deterministic: the start block is a fixed hash of the point ids, and the solves and products sum in fixed chunks. Within a repeated eigenvalue the basis returned is one of many, the same on every run, backend and thread count; compare subspaces, not columns, across meshes. There is no numpy twin (the [`geodesic_distance`](./geodesic.md) precedent).

## Scope

- **The surface is the mesh's `triangle` blocks**: a volume cell is refused naming `extract_surface`, a quad or polygon naming `convert_cells(mode="simplexify")`, a higher-order cell naming `linearize`. Line and vertex blocks ride along; a point in no triangle gets `NaN` modes. A zero-area triangle adds nothing (`num_degenerate`).
- `num_modes` is at least 1 and at most the surface's point count.
- Everything else passes through: points, connectivity, every other data array, regions.

## Report

`eigenvalues`, `iterations`, `max_residual` (the largest relative residual of the wanted pairs), `converged` and `num_degenerate`.

## Everywhere else

| Surface | Call |
| --- | --- |
| Python | `laplacian_eigenmodes(mesh, num_modes=10, *, tolerance=1e-8, max_iterations=0, return_report=False)` |
| Pipeline | `{"Op": "LaplacianEigenmodes", "NumModes": 10, "Tolerance": 1e-8, "MaxIterations": 0}`; counters `Iterations`, `MaxResidual`, `Converged` |
| Both CLIs | `meshioplusplus eigenmodes IN OUT [-k N] [--tolerance T] [--max-iterations N] [--json]` |
| C | `mio_eigenmode_opts_init(&o); o.num_modes = 16; mio_mesh* out = mio_laplacian_eigenmodes(mesh, &o, &report);` |
| Fortran | `out = m%laplacian_eigenmodes(num_modes=16, converged=ok, stat=ierr)` |
| Julia | `r = laplacian_eigenmodes(mesh; num_modes=16)`; `r.mesh`, `r.converged` |
| R | `r <- mio_laplacian_eigenmodes(mesh, num_modes = 16)`; `r$mesh`, `r$converged` |
| WASM | `m.laplacianEigenmodes(mesh, {numModes, tolerance, maxIterations})` returns `{mesh, eigenvalues, iterations, maxResidual, converged}` |
| MCP | tool `laplacian_eigenmodes(input_path, output_path, num_modes, tolerance, max_iterations)` |

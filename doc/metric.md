# Metric toolbox

`prepare_metric(mesh, arrays)` combines and conditions anisotropic metric fields, with no MMG needed: it intersects the metrics of several fields so one mesh resolves all of them, normalises the result to a target vertex count, clamps its sizes and grades it. The output feeds [`adapt`](./mmg.md) (or any metric-driven mesher, through a Medit `.sol`).

```python
import meshioplusplus as mp

mesh = mp.read("flow.vtu")                                  # carries "p" and "T"
for name in ("p", "T"):
    mesh = mp.hessian(mesh, name, location="point")
    mesh = mp.metric_from_hessian(mesh, f"{name}:hessian", output=f"m_{name}")

m = mp.prepare_metric(mesh, ["m_p", "m_T"], target_vertices=20000,
                      hmin=1e-3, hmax=0.2, gradation=1.5)
out = mp.adapt(m, metric="metric")                          # needs the MMG backend
```

## The stages

A metric field is `point_data` with six components `m11 m21 m22 m31 m32 m33` per point in 3-D (three, `m11 m21 m22`, in 2-D): the lower-triangle order [`metric_from_hessian`](./mmg.md#metric_from_hessian) writes and a Medit `.sol` holds. A metric `M` asks for edges of unit length measured in `M`: size `1/sqrt(λ)` along each eigenvector. Per point, in order:

1. **Intersection** of the fields, folded left to right by simultaneous reduction (Frey and George, *Mesh Generation*, 2008, §10.6): with `C = A^(-1/2) B A^(-1/2) = Q diag(μ) Qᵀ`, `A ∩ B = A^(1/2) Q diag(max(1, μ)) Qᵀ A^(1/2)`. It is the smallest metric that asks, in every direction, for at least what each field does (its unit ball is the largest ellipsoid inside both), so a mesh that follows it resolves every field's features at once.
2. **Normalisation** (`target_vertices` > 0): the continuous-mesh L^p normalisation of Loseille and Alauzet ("Continuous mesh framework part II", SIAM J. Numer. Anal. 49, 2011), `M ← N^(2/d) (∫ det(M)^(p/(2p+d)))^(−2/d) det(M)^(−1/(2p+d)) M`, `d` the metric's dimension and the integral over the mesh's `d`-dimensional cells by lumped point measures. It makes the metric's complexity `∫ sqrt(det M)`, about the vertex count a mesher produces from it, equal `N` (before clamping). `norm` is `p`: 1 or 2 spend vertices where the error is (2 is the common choice for the L² interpolation error), `0` means p = ∞, a uniform scaling that keeps the metric's shape.
3. **Clamping** (`hmin`, `hmax`): every eigenvalue into `[1/hmax², 1/hmin²]`, so every size asked lies in `[hmin, hmax]`.
4. **Gradation** (`gradation` = β > 1): Alauzet's metric-space gradation ("Size gradation control of anisotropic meshes", Finite Elem. Anal. Des. 46, 2010). Along every mesh edge `pq`, of length `ℓ` measured in `M_p`, `M_q ← M_q ∩ (1 + ℓ ln β)^(−2) M_p` and the same the other way, swept over the edges in order until no metric changes by more than 1e-3 (or `max_sweeps`): sizes grow by at most a factor β per unit length, so a tiny feature does not sit in a sea of huge cells. It only ever refines. The result is clamped again.

## Results

The roadmap's probe: two fronts, `tanh(40 (x − 0.5))` and `tanh(40 (y − 0.4))` on a unit square, each field's Hessian metric, intersected. In every direction at every point the intersection asks for at least what each metric does (checked over 64 random directions), equals front 1's metric across front 1 and front 2's across front 2, and asks there for sizes over 50 times finer than far from both. Normalising a uniform anisotropic metric to 3 000 vertices gives a complexity of 3 000 for p = 1, 2 and ∞ and keeps its 16 : 1 anisotropy.

On a planar mesh (2-D points, or every `z` equal) with no volume cell, a 3-D metric, which is what `hessian` and `metric_from_hessian` produce there, is normalised over the triangles and quads by its x–y block (`d` = 2). A 3-D metric on a curved surface (MMG's `mmgs` case) can be intersected, clamped and graded, but not normalised, which needs cells of the metric's dimension and says so. Everything else passes through: points, cells, the input fields, every other array. Serial and deterministic. There is no numpy twin (the `metric_from_hessian` precedent).

## Report

`complexity_before` (of the intersection) and `complexity_after` (of the result), `num_graded` (points the gradation refined) and `sweeps`.

## Everywhere else

| Surface | Call |
| --- | --- |
| Python | `prepare_metric(mesh, arrays, *, output="metric", target_vertices=0.0, norm=2.0, hmin=0.0, hmax=0.0, gradation=0.0, max_sweeps=100, return_report=False)` |
| Pipeline | `{"Op": "PrepareMetric", "Arrays": ["m_p", "m_T"], "Output": "metric", "TargetVertices": 20000, "Norm": 2, "Hmin": 0.001, "Hmax": 0.2, "Gradation": 1.5, "MaxSweeps": 100}`; counters `ComplexityBefore`, `ComplexityAfter`, `NumGraded` |
| Both CLIs | `meshioplusplus metric IN OUT --arrays m_p,m_T [--output NAME] [--target-vertices N] [--norm P] [--hmin H] [--hmax H] [--gradation B] [--max-sweeps N] [--json]` |
| C | `mio_metric_opts_init(&o); o.arrays = names; o.num_arrays = 2; o.target_vertices = 2e4; mio_mesh* out = mio_prepare_metric(mesh, &o, &report);` |
| Fortran | `out = m%prepare_metric([character(4) :: "m_p", "m_T"], target_vertices=2d4, stat=ierr)` |
| Julia | `r = prepare_metric(mesh, ["m_p", "m_T"]; target_vertices=2e4)`; `r.mesh`, `r.complexity_after`, ... |
| R | `r <- mio_prepare_metric(mesh, c("m_p", "m_T"), target_vertices = 2e4)`; `r$mesh`, ... |
| WASM | `m.prepareMetric(mesh, {arrays: ['m_p', 'm_T'], targetVertices: 2e4, ...})` returns `{mesh, complexityBefore, complexityAfter, numGraded, sweeps}` |
| MCP | tool `prepare_metric(input_path, output_path, arrays, output, target_vertices, norm, hmin, hmax, gradation)` |

# Mesh operators as sparse matrices

`mesh_operator` returns the mesh's own **P1 finite-element operators** as sparse matrices: the cotangent Laplacian, the lumped and the consistent mass, and the per-triangle gradient of a triangle surface. They are the operators [`geodesic_distance`](/geodesic), [`laplacian_eigenmodes`](/eigenmodes) and [`sobolev_deform`](/sobolev_deform) solve with, returned as data so an application can build its own systems on them. Nothing is solved and no physics is assembled: the solve of a user's system stays the application's.

```python
import meshioplusplus as mp

mesh = mp.icosa_sphere(1.0, 3)
K = mp.mesh_operator(mesh, "laplacian")         # scipy.sparse.csr_matrix
M = mp.mesh_operator(mesh, "mass_lumped")
G = mp.mesh_operator(mesh, "gradient")          # (3 T, n)

f = mesh.points[:, 0]
grad = (G @ f).reshape(-1, 3)                   # the gradient of f on every triangle
lap = (M.diagonal() ** -1) * (K @ f)            # M^-1 K f, the (positive) Laplacian of f
```

| `kind` | Matrix | Shape |
| --- | --- | --- |
| `laplacian` | `K_ij = ∫ ∇φᵢ·∇φⱼ dA`: symmetric, positive semidefinite, rows summing to 0, the natural (Neumann) boundary, no constraint | `n × n` |
| `mass_lumped` | `Mᵢᵢ = Σ A_t / 3` over the triangles at `i` | `n × n`, diagonal |
| `mass_consistent` | `Mᵢⱼ = ∫ φᵢφⱼ dA`: `A_t / 12 (1 + δᵢⱼ)` per triangle | `n × n` |
| `gradient` | row `3 t + k` is component `k` of the gradient of the P1 field on triangle `t`, so `G f` laid out `(t, k)` is the per-triangle gradient; `K = Gᵀ diag(A) G` | `3 T × n` |

`n` is the point count and `T` the number of triangles of the `triangle` blocks, in input order. The rows of a zero-area triangle are empty in `gradient`, and such triangles contribute nothing to the others (`return_report=True` gives `num_degenerate`). Columns are ascending and unique within a row, and the matrices are the same on every backend.

## Smoothing a field

`smooth_field(mesh, array, length, *, output=None, tolerance=1e-12, max_iterations=0, return_report=False)` solves the screened Poisson problem `(M + l² K) g = M f` over the lumped mass and the cotangent stiffness above, and adds `point_data["<array>:smooth"]` (or `output`). Features shorter than `l` are flattened, longer ones pass nearly unchanged, the `M`-weighted mean is preserved exactly, and `length=0` returns the field; the points do not move (that is [`sobolev_deform`](/sobolev_deform)'s job). A vector array is smoothed component by component, a point in no triangle keeps its value, and each solve is Jacobi-preconditioned conjugate gradients, serial and the same on every backend. Like `mesh_operator` it is **C++-core only**, on triangle surfaces, with the verb `smooth-field` and the MCP tool `smooth_field`.

```bash
meshioplusplus smooth-field noisy.vtu smooth.vtu --array curvature --length 0.2
```

The **field Laplacian** needs no operation of its own: `M⁻¹ K f` is `(K @ f) / M.diagonal()` from the matrices above (positive for a function with a maximum where the curvature is negative).

## Scope

The surface is the mesh's `triangle` blocks, as for `geodesic_distance`; a volume cell (`extract_surface` makes one into triangles), a quad or polygon (`convert_cells(mode="simplexify")`) or a higher-order surface cell (`linearize`) is refused by name, and line and vertex blocks are ignored. **Tetrahedra are not offered**: they wait on the shared simplex kernels of roadmap 3.4.4.

## Python

`mesh_operator(mesh, kind="laplacian", *, sparse=None, return_report=False)` is **C++-core only** (no numpy fallback, the `geodesic_distance` precedent). It returns a `scipy.sparse.csr_matrix` when scipy imports, else the tuple `(indptr, indices, values, shape)`; `sparse=True` or `False` forces either.

## CLI and MCP

```bash
meshioplusplus mesh-operator sphere.vtu laplacian.npz --kind laplacian --json
```

writes `indptr`, `indices`, `values` and `shape` to the `.npz` file (read it back with `scipy.sparse.csr_matrix((v, i, p), shape=shape)`); `--json` prints the shape and the entry count. The verb is Python-only, and the MCP server carries it as the `mesh_operator` tool. The C API and the flat bindings (Fortran, Julia, R, Octave, WebAssembly) do not carry `mesh_operator` or `smooth_field` yet: a CSR result handle is its own design (roadmap 15.2).

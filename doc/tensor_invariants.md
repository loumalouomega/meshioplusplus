# Tensor invariants

`tensor_invariants` derives von Mises, principal, hydrostatic and deviatoric fields (and, opt-in, trace, determinant, the three invariants, the inverse and the eigenvectors) from a symmetric or general 3x3 tensor data array. It is a [data operation](/data_operations), not a file format.

Before this operation the only route to these was the CalculiX [`.frd`](/formats/frd) reader's format-specific `derived=` option, which the generic `read`, the CLIs, MCP and the flat bindings did not carry. `.frd`'s `derived=` now calls this operation internally and keeps its existing output names.

## Input layout

An array is treated as one tensor per row, selected by its trailing component count:

- **6 components**: a symmetric tensor `xx yy zz xy yz zx`, the order used throughout this codebase (see [the mesh data model](/mesh_data_model)).
- **9 components**: a general 3x3 tensor, row-major (`xx xy xz yx yy yz zx zy zz`), the layout [`gradient`](/gradient) and [`hessian`](/hessian) produce. `mises` and `principal` are only defined for a symmetric tensor, so both use the symmetric part `0.5*(T + T^T)`.

Any other component count is an error when the array is named explicitly, and is silently skipped when `keys` is omitted (every 6- or 9-component array at the location is processed).

```python
import meshioplusplus as mp

mesh = mp.read("result.vtu")  # carries a "stress" point_data array, 6 components

mesh = mp.tensor_invariants(mesh, "point", ["stress"])
# stress_mises, stress_principal, stress_hydrostatic, stress_deviatoric

# Only the ones you need:
mesh = mp.tensor_invariants(mesh, "point", ["stress"], outputs=["mises"])
```

## Outputs

Each output is written as `prefix + name + "_" + output + suffix`:

| Output | Formula | Shape |
| --- | --- | --- |
| `mises` | `sqrt(0.5*((xx-yy)^2+(yy-zz)^2+(zz-xx)^2+6*(xy^2+yz^2+zx^2)))`, symmetric part for a 9-component input | `(n,)` |
| `principal` | eigenvalues of the symmetric part, ascending (min, mid, max) | `(n, 3)` |
| `hydrostatic` | `(xx+yy+zz)/3`, the mean of the tensor's own (unsymmetrized) diagonal | `(n,)` |
| `deviatoric` | the input with `hydrostatic` subtracted from its three diagonal entries, off-diagonal entries unchanged | same shape as the input |
| `trace` (opt-in) | `xx+yy+zz` of the tensor as given | `(n,)` |
| `determinant` (opt-in) | `det A` of the 3x3, as the triple product of its rows | `(n,)` |
| `invariants` (opt-in) | `I1 = trace`, `I2 = (xx·yy−xy·yx) + (yy·zz−yz·zy) + (xx·zz−xz·zx)` (the sum of the principal 2x2 minors, `½(tr² − tr A²)`), `I3 = det` | `(n, 3)` |
| `eigenvectors` (opt-in) | the unit eigenvectors of the symmetric part, as the row-major 3x3 whose **column** `k` (components `k`, `k + 3`, `k + 6`) belongs to the k-th (ascending) `principal` value, so `A·v = λ·v` for `v = (e[k], e[k+3], e[k+6])` | `(n, 9)` |
| `inverse` (opt-in) | `adj(A) / det(A)`; a row whose determinant is exactly zero is `NaN` | same shape as the input |

`eigenvectors` is opt-in too, but uses the symmetric part of a 9-component input, as `mises` and `principal` do. Each vector's sign is fixed (its first component of magnitude above `1e-8` is positive), so the result does not depend on how the solver's rotations landed. A repeated eigenvalue has no unique eigenvectors: the pair returned is one orthonormal basis of its eigenspace, and another build may return a different one (the values are not affected). Do not compare the vectors of such rows across engines.

The four algebraic outputs are **opt-in**: a default call, and `outputs=None`, still writes only the first four, so existing results do not gain arrays. They read the tensor **as given**, like `hydrostatic` and `deviatoric` and unlike `mises` and `principal`: a 9-component input is not symmetrised, so the determinant of `[[1,2,3],[0,1,4],[5,6,0]]` is 1. A 6-component input is symmetric by construction and its `inverse` has six components, in the same order. There is no warning channel, so a singular row is `NaN` and nothing more; a nearly singular one is inverted as it stands, and conditioning is the caller's to check.

Products of two arrays (`dot`, `cross`, `matmul`) are not offered: see the [roadmap](/roadmap#math-on-meshes).

A row with any non-finite input component produces `NaN` in every requested output for that row.

## Locations

`location="point"` or `"cell"`; `"field"` is rejected — a field array has no per-row correspondence to reduce. For `cell_data`, one output array is written per cell block, matching the mesh's own block structure.

## Data types

Every output is `float64`, like the rest of the arithmetic-producing data operations ([`data_calc`](/data_calc), `point_data_to_cell_data`/`cell_data_to_point_data`).

## CLI

```bash
meshioplusplus data invariants result.vtu out.vtu --names stress
meshioplusplus data invariants result.vtu out.vtu --names stress --outputs mises,hydrostatic
meshioplusplus data invariants result.vtu out.vtu --location cell --names stress --prefix p_
```

Omit `--names` to process every 6- or 9-component array at `--location` (default `point`). `--no-overwrite` fails instead of silently replacing an existing array of the target name. See the [CLI reference](/cli#meshioplusplus-data).

## Other languages

- **C API** — `mio_tensor_invariants(mesh, location, names, count, outputs, prefix, suffix, overwrite)` with the `mio_tensor_invariant` bitmask (`MIO_TINV_MISES`, `_PRINCIPAL`, `_HYDROSTATIC`, `_DEVIATORIC`, `_ALL` (the four of those), and the opt-in `_TRACE`, `_DETERMINANT`, `_INVARIANTS`, `_INVERSE`, `_EIGENVECTORS`). See the [C API reference](/c_api).
- **Fortran** — `m%tensor_invariants(MIO_DATA_POINT, ["stress"], outputs=ior(MIO_TINV_MISES, MIO_TINV_PRINCIPAL))`. See the [Fortran reference](/fortran).
- **Julia** — `tensor_invariants(m, :point, ["stress"]; outputs=[:mises, :principal])`. See the [Julia reference](/julia).
- **R** — `mio_tensor_invariants(mesh, "point", "stress", outputs = c("mises", "principal"))`. See the [R reference](/r).
- **WebAssembly / JavaScript** — `tensorInvariants(mesh, "point", ["stress"], ["mises", "principal"])`. See the [WebAssembly reference](/wasm).

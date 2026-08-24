# Interpolation (field transfer)

`meshioplusplus.interpolate(source, target)` samples a **source** mesh's data arrays onto a **target** mesh — the repo's first cross-mesh operation that transfers data ([diff](/diff) compares two meshes, [merge](/merge) concatenates them; neither resamples). It is a mesh **operation**, not a file format, uses only standard C++/numpy, and runs under every mesh backend.

The result is a **new mesh that is a copy of the target**: geometry, connectivity, the target's own data arrays and its `point_sets`/`cell_sets` are preserved exactly, with the requested source arrays sampled onto it. Source `point_data` is sampled at the target's points; source `cell_data` at the target's cell centroids. Source sets are *not* transferred, and `mesh.info`/`gmsh_periodic` are not carried.

![A coarse source field transferred onto a finer, offset target mesh](/images/interpolate_transfer.png)

```python
import meshioplusplus

coarse = meshioplusplus.read("solution.vtu")   # carries point_data "T"
fine = meshioplusplus.read("remeshed.vtu")

# nearest source point (robust, any cell type, dtype-preserving)
out = meshioplusplus.interpolate(coarse, fine)

# linear interpolation, exact on a linear field
out = meshioplusplus.interpolate(coarse, fine, method="barycentric")

# name the arrays; points outside the source keep -1 instead of erroring
out = meshioplusplus.interpolate(
    coarse, fine, method="barycentric", arrays=["T"], default_value=-1.0
)

meshioplusplus.write("mapped.vtu", out)
```

## Two methods

| method | value at a target sample | dtype | needs |
|---|---|---|---|
| `nearest` (default) | the nearest source **point**'s value, copied bit-for-bit | preserved | source points only |
| `barycentric` | linear interpolation in the containing source **simplex** | Float64 | triangle/tetra cells (after simplexify) |

**Nearest** is piecewise constant and never fails: every target sample finds *some* nearest source point (tie-break: smallest squared distance, then lowest source point index). Values are copied as raw bytes, so integer tags — including `int64` payloads beyond 2⁵³ — survive exactly.

**Barycentric** first simplexifies the source with [convert_cells](/convert_cells) so containment and weights are exact and uniform. Two caveats follow directly from that:

- on a quad/hex source the result is the **simplex-linear** interpolant, not true bi-/trilinear — exact on a linear field either way, and the difference on smooth data is second order;
- only the source's maximum topological dimension is sampled (tetrahedra when any exist, else triangles), and triangle barycentrics are computed **in the xy-plane**. For a curved surface embedded in 3D, use `nearest`.

`cell_data` always transfers by **nearest source-cell centroid**, whatever the method — there is no meaningful barycentric interpolation of a piecewise constant field.

## Arrays and locations

`arrays=None` (the default) transfers **every source `point_data` array**, in sorted name order. `cell_data` transfers only when named explicitly in `arrays`; a name present in both source locations transfers both, and a name present in neither raises. Cross-location transfer (point → cell and back) is deliberately out of scope — compose with the [data operations](/data_operations): `data to-cell` / `to-point`.

## Outside the domain

Under `barycentric`, a target point contained in no source simplex receives `default_value` (in every component), and one warning names how many points that was. Pass `extrapolate=True` to fall back to the nearest source point's value instead. `nearest` has no outside — the nearest point always exists.

## Name conflicts

If a transferred name already exists on the target, `on_conflict` decides: `"error"` (default) raises, `"overwrite"` replaces the target's array, and `"suffix"` writes to `name + "_interp"` — raising if that name is taken too, so nothing is ever silently clobbered.

## Determinism

The output is byte-identical across the three mesh backends, across thread counts, and across the C++-core / numpy-fallback boundary (pinned by `tests/python/test_interpolate.py::test_cpp_matches_python`):

- both search grids (source points, simplex bounding boxes) are bucket-grid spatial hashes (`detail/spatial_hash.hpp`, shared with [merge](/merge)'s weld) built by a parallel key fill and a **serial** ascending insert — no O(N²) anywhere;
- the grid cell size is derived without `cbrt`: the smallest integer `R` with `R·R·R ≥ n`, then one division of the largest bbox extent — an integer loop plus one IEEE division, bit-reproducible in numpy;
- the nearest search expands Chebyshev shells and stops once `((r−1)·cell)²` exceeds the best squared distance; ties break to the lowest source index;
- barycentric containment takes the **first containing simplex in ascending simplex index** (weights ≥ −10⁻¹² — a dimensionless volume ratio, hence scale-invariant), and weights are applied corner by corner in connectivity order, never through `np.sum`.

Every target value is computed independently from the immutable source, so the parallel loops write disjoint slots and thread count cannot affect the result.

## CLI

```sh
meshioplusplus interpolate SOURCE TARGET OUT
meshioplusplus interpolate coarse.vtu fine.vtu out.vtu --method barycentric
meshioplusplus interpolate a.msh b.msh out.vtu --arrays T,v --on-conflict suffix
meshioplusplus interpolate a.msh b.msh out.vtu --extrapolate --default-value=-1
```

Available in both the Python CLI and the native `meshioplusplus` binary — see the [CLI reference](/cli#meshioplusplus-interpolate).

## Other languages

The operation is exposed on every binding surface:

```c
/* C API: arrays as char** + count (NULL / count <= 0 = all point_data) */
mio_mesh* out = mio_interpolate(source, target, "barycentric",
                                NULL, 0, /*extrapolate=*/0, 0.0, "error");
```

```fortran
! Fortran: module-level function (two-mesh input, like mio_merge)
mapped = mio_interpolate(coarse, fine, method='barycentric', stat=st)
```

```js
// WASM: arrays as a JS string array ([] = all point_data)
const out = m.interpolate(coarse, fine, 'barycentric', ['T'], false, 0, 'error');
```

The flat surfaces share the Python surface's semantics; as everywhere else on the flat bindings, `point_sets`/`cell_sets` never cross (the target's sets are re-attached by the Python shim only).

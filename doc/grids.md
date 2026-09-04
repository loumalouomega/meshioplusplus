---
title: Mesh and regular grids
description: Sampling a mesh's fields onto a voxel grid and scattering a prediction back, for convolutional and superresolution models.
---

# Mesh and regular grids

A convolutional model does not take a mesh. It takes a dense rectangular array, and getting a field from an unstructured mesh onto one — and the prediction back off it — is the whole of the data path a grid-shaped architecture needs. This page is that path: a lattice as a value, the two transfers, and the one metric that says whether a super-resolved field has the right small-scale content.

Everything here is pure Python over machinery that already exists. Sampling is [`interpolate`](interpolate.md) onto the lattice's points, so it inherits that operation's exactness for a linear field; the column contract is [`feature_matrix`](ml.md)'s, so a grid's channels and a graph's node features are named by one rule rather than two. Nothing in the C++ core, the WASM build or any binding changed.

```python
import meshioplusplus as mio

mesh = mio.read("case_0042.vtu")

spec  = mio.GridSpec.from_mesh(mesh, resolution=(64, 64, 64))
array = mio.sample_grid(mesh, spec, fields=["T", "vel"])

array.values      # (C, D, H, W) float64 — channels first, then z, y, x
array.channels    # ('T', 'vel_0', 'vel_1', 'vel_2') — the contract
array.coverage    # fraction of grid points inside the mesh

fine   = mio.resample_grid(array.values, spec, spec.upscale(2))   # the baseline
result = mio.scatter_grid(mio.GridArray(fine, spec.upscale(2), array.channels), mesh)
mio.write("predicted.vtu", result)
```

## The array layout

A sampled grid is **`(C, D, H, W)`**: channels first, then **z, y, x**. That is `torch.nn.Conv3d`'s layout and VTK ImageData's, and it is what [`grid`](voxelize.md)'s x-fastest point numbering already produces, so a flat per-point array reshapes into it with no transpose:

```
values[c, k, j, i]      is channel c at lattice point (i, j, k)
```

`GridSpec` carries both spellings and they are deliberately different words:

| | order | means |
|---|---|---|
| `spec.dims` | **world** `(nx, ny, nz)` | **cell** counts — the same spelling `grid`, `voxelize` and `compute_sdf` use |
| `spec.shape` | **tensor** `(D, H, W)` | **point** counts, `(nz+1, ny+1, nx+1)` — the array's own shape |

This is the single most likely thing to get wrong, and the failure is silent: a transposed grid trains, converges and predicts plausible numbers in the wrong places. Every schema this module writes records `"layout": "channels_first_zyx"` verbatim so a checkpoint read months later says which way round it is.

## The shared-box rule

Take the box from the **fine** mesh:

```python
coarse = mio.GridSpec.from_mesh(fine_mesh, resolution=(32, 32, 32))
fine   = coarse.upscale(2)
```

Then every fine node lies inside the coarse grid, and scattering the prediction back is interpolation everywhere rather than extrapolation at the boundary. `upscale` recomputes the spacing from the box rather than dividing it, so the two specs cover the same domain exactly — `spec.same_bounds(other)` checks it and `spec.scaling_factor(other)` returns the integer ratio an SRResNet needs, or `None` when the pair does not have one.

Sampling happens at lattice **points**, not cell centres, for the same reason: the points reach the box's faces and the centres do not.

## Coverage, and why it is reported

A grid around a concave domain is mostly outside it. Points outside the source take `fill_value` (or their nearest source value with `extrapolate=True`), and a model trained on such a grid learns the fill as if it were physics.

`array.coverage` is the fraction of lattice points inside the source. It is free when `extrapolate` is off — a constant probe channel rides the same interpolation call, and barycentric weights summing to one make it an exact mask — and costs one extra pass when `extrapolate` is on, since "inside" then has to be established separately. Pass `fill_value=float("nan")` to make the fill visible to every downstream reduction rather than plausible.

## Trilinear, and why it is not `interpolate`

[`interpolate`](interpolate.md)'s `"barycentric"` mode would also evaluate a field on a lattice, and both are exact for a **linear** field, but they are not the same interpolant and the difference is not small. On a unit cell carrying `u = x·y`:

| point | exact | `interpolate_grid` | `interpolate(..., "barycentric")` |
|---|---|---|---|
| (0.5, 0.5, 0.5) | 0.25 | **0.25** | 0.5 |
| (0.25, 0.75, 0.5) | 0.1875 | **0.1875** | 0.25 |

Trilinear reproduces the field exactly — `x·y` *is* trilinear — while the simplex decomposition is off by a factor of two at the centre, the chosen diagonal showing through. It is also what `torch.nn.functional.interpolate` does and what every superresolution baseline is reported against. Adding it as a third `interpolate` method would additionally break that function's byte-identical C++/numpy twin contract, for an algorithm that only makes sense on a dense lattice.

Queries outside the box are **clamped** to it rather than extrapolated: a trilinear form extrapolates as a product of linears and diverges quickly, and under the shared-box rule an outside query is a numerical excursion at a face rather than a real request.

## Caching grids

Write them as **`.vti`**. VTK ImageData stores the lattice as origin/spacing/extent and regenerates the points arithmetically, so `GridArray.from_mesh` recovers the spec exactly, with no tolerance and no sidecar:

```python
mio.write("coarse.vti", array.to_mesh())
back = mio.GridArray.from_mesh(mio.read("coarse.vti"))
assert back.spec == array.spec
```

A grid written to `.vtu` or `.msh` is a lossy round trip for this purpose — the points survive as coordinates but the *lattice* does not — and reading one back raises by name rather than accepting a warped grid. `to_mesh()` produces an ordinary `hexahedron` mesh, so [`view`](viewer.md), [`crop`](crop.md), [`isosurface`](isosurface.md) and the browser viewer all work on a grid with no new code; that is the reason the grid is a mesh here rather than a bespoke container.

## Power spectra

The honest measure of whether a super-resolved or generated field carries the right *small-scale content*, as opposed to a plausible-looking one with the high wavenumbers smoothed away. A pointwise error cannot see the difference — two single-mode fields of the same amplitude have identical RMS and completely different spectra — so a superresolution result reported only as an error norm has not been evaluated at the scales it exists to reproduce.

```python
ps = mio.power_spectrum(array.channel("T"), spec)
ps.wavenumber, ps.power, ps.counts     # cycles per unit length
```

Three properties worth knowing. The power sums **exactly** to `mean(field**2)` (Parseval), which holds only because every mode is binned — including the ones past the isotropic Nyquist that live in the corners of the box; truncating the tail would break it, so `counts` is reported instead and a caller cuts it themselves. A vector field's components are summed, giving the energy spectrum. And an **anisotropic** lattice raises: averaging over shells is only meaningful when a shell means the same thing on every axis.

## Two-dimensional operators

FNO, AFNO and 2-D U-Nets take `(C, H, W)`. A lattice always has at least two points on every axis — one cell has two corners — so the thin axis of a planar problem still arrives with a plane at each face:

```python
plane = mio.squeeze_grid(array.values, axis=2, index=0)   # world z, the lo face
back  = mio.expand_grid(plane, axis=2)
```

`axis` is a **world** axis (0=x, 1=y, 2=z), like every other three-vector here, not the tensor axis it maps to. `index` is required whenever the axis is longer than one, because silently taking the first plane would discard the other without saying so.

## Command line

```bash
meshioplusplus grid-sample   case.vtu grid.vti --resolution 64,64,64
meshioplusplus grid-resample grid.vti fine.vti --factor 2
meshioplusplus grid-scatter  fine.vti case.vtu predicted.vtu
meshioplusplus grid-spectrum grid.vti --field T --json
```

`grid-sample` takes the same six-field lattice vocabulary as [`voxelize`](voxelize.md) and [`compute_sdf`](sdf.md) — `--resolution` or `--cell-size`, plus `--bounds`, `--padding`, `--padding-relative` and `--max-cells` — because "the grid around this mesh" should mean one thing across the library.

## MCP

Four tools: `grid_sample`, `grid_scatter`, `grid_resample` and `grid_power_spectrum`, all path-in/path-out like every other tool, so an agent can drive the whole data path without holding a mesh. See [the MCP server](mcp.md).

## Pairing a coarse grid with a fine one

`upscale` is right for resampling and **wrong** for a superresolution pair, so there are two methods and they are named apart:

| | multiplies | on a 4×4×4-cell grid at ×2 | use for |
|---|---|---|---|
| `spec.upscale(s)` | **cells** | 9×9×9 points; all 125 coarse points nest | `resample_grid` |
| `spec.upscale_samples(s)` | **samples** | 10×10×10 points; only the 8 corners are shared | a model pair |

A convolutional upsampler multiplies sample counts — `SRResNet(scaling_factor=2)` maps `(B, C, 5, 5, 5)` to `(B, C, 10, 10, 10)` — so a pair built with `upscale` fails as a shape error deep inside a loss. Both preserve the box exactly; `spec.sample_scaling_factor(other)` returns the single integer such a model must be parametrized by, or `None` when the pair does not have one.

See [paired cases](datasets.md#paired-cases) for describing the two sides in a manifest and [`grid_sample_pair`](physicsnemo.md#grid-samples) for producing them.

## What is not here

An `srresnet` family in the training spec is the next item on [the roadmap](roadmap.md). `cell_data` is deliberately refused in both directions: a piecewise-constant field has no value at a point, so convert it with `cell_data_to_point_data` (CLI `data to-point`) first, which makes the approximation explicit rather than hiding it inside the sampler.

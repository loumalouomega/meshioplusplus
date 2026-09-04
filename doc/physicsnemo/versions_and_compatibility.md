---
title: Versions and compatibility
description: The physicsnemo release this documentation describes, the optional extras and what each pulls in, what changed between 2.1 and 2.2, and how to upgrade without being fooled by green tests.
---

# Versions and compatibility

This section describes **physicsnemo 2.2.0**, which is what meshio++'s adapter is exercised against. Upstream's 2.2.1 is a fix-only release (Python 3.14 packaging and annotation introspection) and the 2.3.0 changelog was still empty when this was written, so 2.2 is the current surface. Where a page says an API "does not exist" or "raises", it means in 2.2.

::: warning This page ages faster than the rest
Every other page here describes ideas that stay true. This one describes a snapshot of a moving library, and it was last checked against the versions in the table below. Treat it as a baseline to compare against, not as a current fact.
:::

## What was installed when this was written

None of these is a *requirement* — meshio++ declares no ML dependency at all — but if something behaves differently for you, this is the baseline.

| Package | Version | Role |
|---|---|---|
| `nvidia-physicsnemo` | 2.2.0 | everything ML |
| `torch` | 2.13.0 (CUDA 13.0 build) | everything ML |
| `tensordict` | 0.14.0 | the mesh tensorclass and datapipe samples; 2.2 requires `tensordict[zarr]>=0.14` |
| `warp-lang` | 1.16.0 | Warp kernels: ball query, SDF, remeshing, the FFD deformer; 2.2 requires `>=1.14` |
| `zarr` | 3.3.0 | mesh and datapipe Zarr I/O |
| `torch-geometric`, `torch_scatter` | 2.8.0, 2.1.2 | the graph models (`PYG_AVAILABLE` upstream checks both) |
| `onnx`, `onnxscript`, `onnxruntime` | 1.22.0, 0.7.1, 1.29.0 | ONNX export (torch's exporter needs `onnxscript`) and CPU inference |
| `gpytorch` | 1.15.2 | the GP uncertainty head (`uq-extras`) |
| `pyvista` | 0.48.4 | mesh rendering |
| `cupy-cuda13x` | 14.1.1 | GPU array interchange (upstream's `cu13` extra pins `<14`; 14.1.1 works) |

## The extras, and what each pulls in

`pip install "nvidia-physicsnemo[<extra>]"`. Read from the installed package's metadata:

| Extra | Pulls in | You want it for |
|---|---|---|
| `cu12` / `cu13` | `cuml`, `cupy`, `nvidia-dali`, `pylibraft` for that CUDA major, plus torch/torchvision | GPU-side data loading and clustering kernels; pick the one matching your driver |
| `mesh-extras` | `matplotlib`, `pyvista`, `vtk` | `physicsnemo.mesh` visualization and pyvista conversion |
| `datapipes-extras` | `dask`, `netcdf4`, `tfrecord`, `xarray`, `zarr` | the climate and Zarr readers |
| `gnns` | `torch-geometric`, `torch-scatter`, `torch-sparse`, `torch-cluster`, plus `pyvista`, `vtk`, `stl`, `scipy`, `mlflow`, `wandb` | every graph model; GraphCast needs `torch-sparse` |
| `model-extras`, `nn-extras`, `utils-extras` | `scipy`, `stl`, `vtk`, `mlflow`, `wandb`, `line-profiler` | logging and profiling helpers, STL I/O |
| `sym` | `sympy` | `physicsnemo.sym` (the module itself is bundled; only SymPy is extra) |
| `uq-extras` | `gpytorch` | the variational GP heads |
| `natten-cu12` / `natten-cu13` | `natten` | neighborhood attention layers |
| `transformer-engine-cu12` / `-cu13` | `transformer-engine[core,pytorch]` | `use_te=True` on Transolver-family models (fp8 attention) |

The core dependencies (always installed) include `torch`, `tensordict[zarr]`, `warp-lang`, `h5py`, `onnx`, `hydra-core`/`omegaconf`, `s3fs`/`fsspec`, `timm`, `nvtx`, `treelib` and `einops`.

**One packaging trap worth knowing before it costs you an afternoon:** the graph models import `torch_scatter` directly, and its prebuilt wheels lag torch releases by weeks to months. A perfectly ordinary `pip install torch` can therefore leave you with a torch for which no `torch_scatter` wheel exists yet, and the failure appears deep inside a model's forward pass rather than at install time. Check the wheel index for your torch version before pinning it.

## What changed from 2.1 to 2.2

The lesson that recurs across this table is worth stating on its own: **a layout change that keeps shapes identical passes every test whose fixture is symmetric.**

| Change in 2.2 | Why it matters |
|---|---|
| Mesh-calculus gradients are derivative-first `(N, D, C)` from *every* backend (2.1's least-squares backend was channel-major) | the canonical silent break — a symmetric canary gradient cannot see it, so anything consuming a gradient needs a canary that is asymmetric *and* non-square |
| `ShardTensor` moved from `physicsnemo.distributed` to `physicsnemo.domain_parallel` | an import path, loud and easy |
| `remesh` is Warp-backed (`pyacvd` dropped) and its count targets output *vertices*, not cells; it raises for anything but a surface in 3-D | a target count that used to mean cells now means vertices — same type, different meaning, no error |
| `Mesh.save`/`load` gained a Zarr backend (`mesh.io.to_zarr`/`from_zarr`) | additive; removes the reason to install the curator for Zarr alone |
| GeoTransolver and FLARE promoted out of `experimental` | with `ShardTensor` support and activation checkpointing |
| `nn.functional.signed_distance_field` returns a 3-tuple `(sdf, hit_points, hit_faces)` | a tuple where a tensor used to be, so it fails loudly |
| `fill_interior` gained exact-boundary 2-D filling; its `n = 3` still raises `NotImplementedError` | tetrahedral filling of a closed surface is still not upstream's job |
| `export_to_onnx_stream` no longer runs the model twice, but still exposes no `dynamic_axes` | export through `torch.onnx.export` directly if an axis must stay dynamic |
| `datapipes/protocols.py` rewritten (`_PrefetchResult` became `HostPayload`) | only affects code subclassing the protocol internals |
| `integrate` gained `nan_policy`; `integrate_cell_data`/`integrate_point_data` deprecated in favour of `integrate(...)` | a deprecation with a straightforward replacement |
| The legacy diffusion modules (`samplers.legacy_deterministic_sampler`, `metrics.legacy_losses`, `preconditioners.legacy`) now warn they will be deprecated | new work should target the protocol API |
| Added: `shrink_and_perturb_`, the mesh deformers (`sobolev_deform`, `shrinkwrap`, RBF, FFD), the fixed-topology energies, the grid divergence/curl/Laplacian functionals, `farthest_point_sampling`, FSDP2 checkpoint support | all additive |
| `poisson_sample_indices_fixed` removed | |

## How to check what you have

```python
import physicsnemo, torch, tensordict
print(physicsnemo.__version__, torch.__version__, tensordict.__version__)
```

## Upgrading: compare skip counts, not pass/fail

Two things make a green test run after an upgrade *insufficient* evidence, and both apply to any project with optional ML dependencies:

1. **Optional-dependency gates skip green.** Graph, ONNX, GP-head and visualization tests each self-skip when their import fails. If an upgrade renames the symbol a gate probes, whole classes of test turn into skips and the run stays green. Record the skip count before and after; it must not grow.
2. **Symmetric fixtures hide layout flips.** A unit cube hides every length scale and a field with a symmetric gradient hides a transpose. Any new fixture wants a non-unit extent and an asymmetric canary for exactly this reason.

The upstream [release notes](https://docs.nvidia.com/physicsnemo/latest/release-notes/index.html) and repository changelog list every change per release; the table above is the subset that touches a mesh pipeline.

Back to [PhysicsNeMo basics](./overview.md), or on to [PhysicsNeMo integration](../physicsnemo.md) for what meshio++ implements.

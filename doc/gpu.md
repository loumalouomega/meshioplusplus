# GPU handoff (DLPack / CuPy)

Move a meshio++ `Mesh`'s arrays to and from device memory — **without a file round-trip** — through the standard array-exchange protocols, so meshio++ output feeds CuPy, the RAPIDS stack, PyTorch, JAX and Numba directly, and device arrays can be written back out through every existing writer.

Be precise about what this does and does not avoid: **a host→device move is always a transfer across the bus (PCIe/NVLink) — there is no such thing as a free one.** What this feature removes is (a) the file round-trip people do today, and (b) every *extra* copy on either side of that one transfer, by handing buffers over through protocols the consumer adopts in place. On this page "zero-copy" appears only where it is literally true — sharing host buffers into the payload, and on-device adoption of the returned CuPy arrays — never for the transfer itself. The interop layer's [`zero_copy_only` contract](./interop#the-zero-copy-contract) is a **host** buffer-sharing contract and stays exactly that.

```python
import meshioplusplus as mio

mesh = mio.read("bracket.msh")
gpu = mio.to_cupy(mesh)                 # ONE host→device transfer per array
gpu.point_data["T"] *= 2.0              # any CuPy/RAPIDS kernel from here
mio.write("scaled.vtu", mio.from_cupy(gpu))   # one device→host copy, then any writer
```

Nothing here is part of the C++ core, which stays dependency-free. This is pure Python (`src/python/meshioplusplus/_gpu.py`) over the numpy the readers already return; every third-party import is lazy and named.

## Installation

There is **deliberately no pip extra** for this feature. CuPy has no single installable distribution — its wheels are CUDA-version-specific — so a `meshioplusplus[gpu]` extra pinning any one of them would fail to resolve or install the wrong binary for most users. Install the wheel matching your toolkit yourself:

| Your toolkit | Install |
|---|---|
| CUDA 13.x | `pip install cupy-cuda13x` |
| CUDA 12.x | `pip install cupy-cuda12x` |
| CUDA 11.x | `pip install cupy-cuda11x` |
| ROCm | the matching `cupy-rocm-*` wheel (untested by meshio++) |

The install error names this recipe rather than a pip extra. Two predicates answer the availability question without raising, and they are **not** the same question: `has_cupy()` (importable — true on a machine with no GPU at all, since CuPy fails only at first use) and `has_cuda_device()` (importable *and* at least one device actually answers).

`to_dlpack` and `from_cupy`'s host paths need only numpy ≥ 1.23 (where `np.from_dlpack` and `ndarray.__dlpack__` exist) — no CuPy, no GPU.

## Architecture: the pure payload layer

`_gpu.py` is split exactly the way [`_interop.py`](./interop) is: the bulk is the pure builder `_to_device_payload`, which imports no third-party library, mutates nothing, and returns plain numpy in a frozen `DevicePayload` — which is what makes everything subtle here testable in the default CI matrix on a machine with no GPU. `to_cupy` / `from_cupy` are thin gated wrappers.

A `DevicePayload` carries:

| Field | Contents |
|---|---|
| `points` | the `(n, dim)` coordinates |
| `cells` | one `(type, connectivity)` pair per rectangular cell block |
| `block_bases` | prefix sums over the payload's own blocks — the block-major cell numbering, from `_regions.block_bases`, never re-derived |
| `point_data` / `cell_data` | name → array / name → per-block tuple of arrays |
| `field_data` | host-side passthrough, never transferred |
| `regions` | `(name, kind, dim, tag, entries)` index arrays — see below |
| `device` | the DLPack `(device_type, device_id)`: `(1, 0)` host, `(2, i)` CUDA |
| `dropped`, `shared`, `notes`, `derived`, `refs` | bookkeeping (see below) |

`cell_data` stays **per-block** — the mesh-faithful shape, since blocks may disagree on component shape — and `block_bases` rides along so a consumer that wants the flat block-major array concatenates once (`cp.concatenate`) without re-deriving an offset. Ragged `polygon`/`polyhedron` blocks cannot be a DLPack tensor and are dropped with a note, listed in `dropped`.

Every lossy or copying step lands in `notes`, surfaced as a warning — the interop layer's standing convention. `derived` records data-model facts that are not losses (the device transfer itself, region renumbering after a ragged drop) and is never warned.

## DLPack: the primary export (host too)

```python
payload = mio.to_dlpack(mesh)           # host numpy, nothing transferred
t = torch.from_dlpack(payload.points)   # adopted in place
a = np.from_dlpack(payload.point_data["T"])   # ditto — genuinely zero-copy
```

`to_dlpack(mesh)` returns a host `DevicePayload` whose every array is a numpy ndarray — and a numpy ndarray natively implements `__dlpack__` / `__dlpack_device__` (reporting `(kDLCPU, 0)`). The payload is the container; **each array speaks the protocol itself**, so hand any of them straight to `np.from_dlpack`, `torch.from_dlpack`, `jax.dlpack`, `cupy.from_dlpack`, … DLPack was chosen as the primary because it is framework-neutral *and* covers host arrays — one export path serves CPU and GPU consumers, and the host half genuinely exercises the machinery in CI with no GPU involved.

This is a host buffer-sharing export, squarely inside the [zero-copy contract](./interop#the-zero-copy-contract): canonical float64/int64 arrays are shared as-is, anything that would copy is noted, and `to_dlpack(mesh, zero_copy_only=True)` raises rather than copying silently. Lifetime works the way numpy's own DLPack does: the exported capsule keeps the exporting array alive, which keeps the reader's C++-owned buffer alive.

## PyTorch / JAX convenience: `to_torch`, `to_jax`

When the consumer is known, the per-array `torch.from_dlpack` calls above are already wrapped for you:

```python
t = mio.to_torch(mesh)                 # torch tensors, host, zero-copy
t = mio.to_torch(mesh, device="cuda")  # + one bus transfer per array
j = mio.to_jax(mesh)                   # jax arrays, JAX's default device
```

Both return the same `DevicePayload` shape with every array adopted for you — `to_torch` on the host genuinely zero-copy, `device=` an honest recorded transfer; `to_jax` on JAX's default device with a documented x64 caveat. There is deliberately no `[torch]`/`[jax]` pip extra (the CuPy reasoning). The full story — including the PyG `Data` recipe with `edge_index`/`feature_matrix` — lives on the [ML data handling](./ml) page.

## dtypes: explicit opt-in downcasts only

meshio++'s canonical storage is float64/int64 and that is what you get by default. Most GPU pipelines want fp32, so both exports take explicit opt-ins — never silent, each cast recorded in the warned `notes`:

- `float32=True` — every float array (`points`, float `point_data` / `cell_data`) is cast to float32; integer arrays are untouched.
- `int32=True` — connectivity, `block_bases` and region entries are cast to int32 **after a range check**: an index outside the int32 range raises, naming the array, rather than wrapping. int32 is the floor — there is no int16 option.

Combining either with `zero_copy_only=True` raises immediately: an explicit downcast is a copy by definition.

## CuPy: the device transfer

```python
gpu = mio.to_cupy(mesh)                       # requires a working CUDA device
gpu = mio.to_cupy(mesh, float32=True, int32=True)   # kernel-friendly dtypes
back = mio.from_cupy(gpu)                     # one device→host copy per array
```

`to_cupy` builds the same payload and moves it: **one bus transfer per array, never called "zero-copy"** — and there is deliberately no `zero_copy_only` parameter here, because it would be a lie on a path that is 100% copies. The returned CuPy arrays expose both DLPack and `__cuda_array_interface__` (CAI), so everything downstream (RAPIDS, Numba, torch via DLPack) adopts them on-device without further copies.

`from_cupy` accepts a `DevicePayload` or any duck-typed object with the same attributes, and each array may be anything speaking DLPack or CAI — a CuPy array, a torch tensor, a Numba device array:

- a numpy array passes through;
- a **host** DLPack exporter (`kDLCPU` — a torch CPU tensor, a JAX host array) is adopted zero-copy via `np.from_dlpack`, no CuPy needed;
- a device-resident array comes back with **one deliberate device→host copy** (its own `.get()` when it has one, else CuPy adopts it on-device via CAI/DLPack — zero-copy — and copies back once). That copy is the documented cost of returning to the CPU world, not a trick, and it is noted.

The result is an ordinary `Mesh`: every existing writer, operation and binding works on it unchanged.

### Regions on device

`Point` and `Cell` regions are index arrays and `Side` regions are `(cell, facet)` pairs — all three ride along **as device index arrays**, carrying `name`/`kind`/`dim`/`tag`. This is deliberately *not* the mask-plus-JSON-sidecar convention the [PyVista/trimesh export](./interop#regions) uses: that convention exists only because those targets have no structured slot to put a region in. `DevicePayload` is meshio++'s own type, so the index arrays stay what they are — lossless (`Side` pairs survive, which no other interop target manages), smaller than masks, and directly usable as gather indices in a kernel. `cell` and `side` entries are numbered in the payload's own block-major order (which differs from the mesh's only when a ragged block was dropped — then they are renumbered at export, and entries in dropped blocks are removed with a note).

### Pinned staging

```python
gpu = mio.to_cupy(mesh, pinned=True)                  # synchronous, staged
stream = cupy.cuda.Stream(non_blocking=True)
gpu = mio.to_cupy(mesh, pinned=True, stream=stream)   # async DMA
stream.synchronize()                                  # ...before using gpu
```

`pinned=True` stages each array through page-locked (pinned) host memory, which is what allows the DMA engine to run a true async copy. With the default `stream=None` the call synchronizes before returning and the staging buffers are released. With a caller-supplied stream the copies are enqueued asynchronously and the staging buffers are kept alive on `payload.refs` — the caller **must synchronize the stream** before using the arrays or dropping the payload.

Honestly scoped: pinned staging pays off for **repeated** transfers of large buffers; for a one-shot transfer of a small mesh the extra host-to-pinned copy can cost as much as it saves, and meshio++ has not benchmarked a crossover point — measure on your own meshes before turning it on. The bigger win — reading files **directly into** pinned memory, skipping the staging copy entirely — is Phase 2 (below).

## Testing and CI

The pure payload builder and the DLPack **host** round-trip (export, re-adopt with `np.from_dlpack`, assert values and that sharing really held) run in the default CI matrix and the `interop` job — no GPU involved, and that half is real coverage of the export machinery, not a mock.

**The device path is not covered by public CI.** GitHub-hosted runners have no CUDA device, so the gated tests in `tests/python/test_gpu.py` — guarded by `pytest.importorskip("cupy")` *plus* an actual `getDeviceCount() > 0` check — skip there, and the CuPy-only lines read as uncovered on the Codecov patch check by design. To run the device suite, use any machine with a CUDA device and a matching CuPy wheel (a self-hosted or cloud GPU runner does fine):

```bash
pip install -e . cupy-cuda13x pytest    # or the wheel matching your toolkit
python -m pytest tests/python/test_gpu.py -v      # gated tests now run
```

## Phase 2: reading directly into pinned memory

By default the readers allocate ordinary pageable memory, so a pinned transfer stages through one extra host copy. **The C++ half of removing that copy is in place since v8.5.0**: every owning `NDArray` buffer is allocated through an optional process-global hook —

```cpp
#include <meshioplusplus/ndarray.hpp>

// meshioplusplus::BufferAllocator: plain C callbacks {alloc, free, user}.
meshioplusplus::set_buffer_allocator(my_pinned_allocator);  // nullptr = default heap
```

— consulted only at allocation time. Each buffer keeps its own `shared_ptr` reference to the allocator it was born with, so arrays outlive any later `set_buffer_allocator` call: uninstalling the hook never orphans live buffers, and a buffer allocated pinned is freed through the pinned pool no matter when it dies. Two consequences worth knowing: a *copy* of an array allocates through the allocator current at copy time (pinned-ness does not propagate through copies), and views are unaffected (they own nothing). Content, zero-copy and determinism contracts are unchanged; the capsule handoff to numpy needs no change — the buffer's own deleter fires when numpy drops the last reference.

**The Python/CuPy wiring is the remaining follow-up**: a `meshioplusplus._gpu.pinned_reads()` context manager installing a `cupy.cuda.alloc_pinned_memory`-backed allocator (the C++ shim must hold the `py::function` pair and a `ptr -> PinnedMemoryPointer` map, take the GIL in both callbacks, and guard the free path with `Py_IsInitialized()` for interpreter-shutdown stragglers), plus a `to_cupy(pinned=True)` fast path that skips the staging copy for arrays the hook allocated. Recorded the way [Open3D and DOLFINx's constraints](./interop#phase-2-open3d-and-dolfinx) are: known, bounded, and waiting on a real need rather than speculation.

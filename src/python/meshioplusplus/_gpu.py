"""GPU handoff at the I/O boundary — DLPack export, CuPy transfer.

Move a :class:`~meshioplusplus._mesh.Mesh`'s arrays to and from device memory
**without a file round-trip**, through the standard exchange protocols:

* **DLPack** (the primary export) — :func:`to_dlpack` returns a payload of host
  numpy arrays, every one of which natively speaks ``__dlpack__`` /
  ``__dlpack_device__`` (numpy ≥ 1.23). DLPack is framework-neutral (CuPy,
  PyTorch, JAX, TensorFlow, Numba, the RAPIDS stack) and — the property this
  module leans on — it also covers **host** arrays (device type ``kDLCPU``), so
  the same export path serves CPU consumers and is fully testable with no GPU.
* **``__cuda_array_interface__``** (CAI, the older CUDA-specific protocol) —
  **consumed** on ingest (:func:`from_cupy` accepts any object exposing it) and
  exposed for free on every CuPy array :func:`to_cupy` returns.

Be precise about what is and is not avoided: **a host→device move is always a
transfer across the bus (PCIe/NVLink) — there is no free one.** What this
module removes is (a) the file round-trip, and (b) every *extra* copy on either
side of that one transfer, by handing buffers over through protocols the
consumer adopts in place. The word "zero-copy" is used here only for the two
places it is literally true: sharing host buffers into the payload
(:func:`to_dlpack`, governed by the same ``zero_copy_only`` contract as
:mod:`meshioplusplus._interop` — a **host** buffer-sharing contract), and
on-device adoption of the returned CuPy arrays via DLPack/CAI. The device
transfer itself is never called that. Symmetrically, :func:`from_cupy` performs
**one deliberate device→host copy per array** and says so.

Nothing here is part of the C++ core. This is pure Python over the numpy the
readers already return; every third-party import is lazy and named.

There is deliberately **no pip extra** for CuPy: its wheels are
CUDA-version-specific (``cupy-cuda13x``, ``cupy-cuda12x``, ``cupy-rocm-*``), so
a ``meshioplusplus[gpu]`` extra pinning any one of them would break for most
users. Install the wheel matching your toolkit; the install error names the
recipe. See ``doc/gpu.md``.

Architecture
------------
Split exactly like :mod:`meshioplusplus._interop` and for the same reason: the
bulk is the **pure payload builder** :func:`_to_device_payload`, which imports
no third-party library, does not mutate its input, and returns plain numpy in a
frozen dataclass — testable in the default CI matrix on a machine with no GPU.
:func:`to_cupy` / :func:`from_cupy` are thin gated wrappers. Block-major cell
numbering comes from :func:`~meshioplusplus._regions.block_bases` (the Python
owner of ``detail/cell_index.hpp``) — never re-derived — and ``block_bases``
rides in the payload so a GPU consumer can concatenate per-block arrays
block-major without re-deriving an offset either.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, NamedTuple

import numpy as np

from ._interop import (
    InteropError,
    _copy_note,
    _derived_note,
    _emit,
    _importable,
    _is_native,
    _prepare_data_array,
)
from ._regions import block_bases

#: DLPack device types this module cares about. ``kDLCPU`` is what every host
#: array reports; CUDA is 2 and ROCm 10 (recorded as CuPy reports them —
#: nothing here special-cases ROCm, and it is documented as untested).
_KDLCPU = 1

_I32 = np.iinfo(np.int32)


class GpuCellBlock(NamedTuple):
    """One rectangular cell block: the meshio++ type name + its connectivity."""

    type: str
    data: Any


class GpuRegion(NamedTuple):
    """A named region as plain index arrays.

    Unlike the PyVista/trimesh export, regions here are **not** flattened into
    mask arrays with a JSON sidecar — that convention exists only because those
    targets have no structured slot for ``(name, kind, dim, tag, entries)``.
    This payload is our own type, so the index arrays ride along unchanged:
    lossless (``dim``/``tag`` carried; ``side`` ``(cell, facet)`` pairs
    survive, which no other interop target manages), smaller than a mask, and
    directly usable as gather indices in a kernel. ``cell`` and ``side``
    entries are numbered in the **payload's own** block-major order.
    """

    name: str
    kind: str
    dim: int
    tag: int
    entries: Any


@dataclass(frozen=True)
class DevicePayload:
    """A mesh's arrays, ready to hand to a device (or already there).

    ``device`` is the DLPack ``(device_type, device_id)`` pair: ``(1, 0)`` for
    a host payload from :func:`to_dlpack`, ``(2, i)`` for CUDA device ``i``
    from :func:`to_cupy` (ROCm reports 10). ``field_data`` is a host-side
    passthrough and is never transferred. ``refs`` is the keep-alive for
    asynchronous pinned staging buffers — see :func:`to_cupy`.
    """

    points: Any
    cells: tuple  # of GpuCellBlock, ragged blocks dropped
    block_bases: Any  # prefix sums over the payload's OWN cells
    point_data: dict  # name -> array
    cell_data: dict  # name -> tuple of per-block arrays, aligned with cells
    field_data: dict  # host passthrough, never transferred
    regions: tuple  # of GpuRegion
    device: tuple  # DLPack (device_type, device_id)
    dropped: tuple  # ((block_index, cell_type), ...) blocks that could not ride
    shared: dict  # label -> bool, host buffer sharing only
    notes: tuple  # every lossy or copying step, warned by the public wrappers
    derived: tuple  # constructions that are a data-model fact, never warned
    refs: tuple  # keep-alive: pinned staging buffers of an async transfer


# --------------------------------------------------------------------------- #
# pure payload layer (no third-party imports)                                 #
# --------------------------------------------------------------------------- #
def _maybe_float32(a, label, float32, notes, shared):
    """Apply the explicit ``float32=True`` opt-in downcast, recorded, never silent."""
    if not float32 or a is None or a.dtype.kind != "f" or a.dtype == np.float32:
        return a
    notes.append(f"{label} was downcast {a.dtype.name}→float32 (float32=True)")
    shared[label] = False
    return a.astype(np.float32)


def _checked_int_cast(a, label, notes, shared, *, op):
    """int64 → int32 with a range check — int32 is the floor, there is no int16."""
    if a.dtype == np.int32:
        return a
    if a.size and (int(a.max()) > _I32.max or int(a.min()) < _I32.min):
        raise InteropError(
            f"meshio++: {op}: {label} cannot be downcast to int32: values span "
            f"[{int(a.min())}, {int(a.max())}], outside the int32 range"
        )
    notes.append(f"{label} was downcast {a.dtype.name}→int32 (int32=True)")
    shared[label] = False
    return a.astype(np.int32)


def _per_block_cell_data(mesh, name, kept, notes, shared, zero_copy_only, float32, op):
    """``cell_data[name]`` for the kept blocks, one prepared array per block.

    Deliberately per-block (the mesh-faithful shape — blocks may disagree on
    component shape), unlike ``_interop._flatten_cell_data``: the payload
    carries ``block_bases`` so a consumer that wants the flat block-major array
    concatenates once, without re-deriving an offset. The drop policies mirror
    ``_flatten_cell_data`` — refusing beats guessing.
    """
    blocks = mesh.cell_data[name]
    if len(blocks) != len(mesh.cells):
        notes.append(
            f"cell_data '{name}' has {len(blocks)} block(s) for "
            f"{len(mesh.cells)} cell block(s); dropped"
        )
        return None
    out = []
    for b in kept:
        blk = blocks[b]
        if blk is None:
            notes.append(f"cell_data '{name}' has an empty block; dropped")
            return None
        label = f"cell_data '{name}' block {b}"
        a = _prepare_data_array(op, label, blk, notes, shared, zero_copy_only)
        if a is None:
            return None
        if a.shape[0] != len(mesh.cells[b]):
            notes.append(
                f"cell_data '{name}' block {b} has {a.shape[0]} row(s) for the "
                f"{len(mesh.cells[b])}-cell '{mesh.cells[b].type}' block; dropped"
            )
            return None
        out.append(_maybe_float32(a, label, float32, notes, shared))
    return tuple(out)


def _remap_cell_entries(region, old_bases, kept, new_bases, notes):
    """Renumber a cell/side region into the payload's own block-major numbering.

    Only reached when a ragged block was dropped. Entries whose cell lies in a
    dropped block are removed with a note (the ``_region_payload`` policy);
    ``side`` entries keep their facet column untouched.
    """
    entries = np.asarray(region.entries, dtype=np.int64)
    side = region.kind == "side"
    cells_col = entries[:, 0] if side else entries
    blocks = np.searchsorted(old_bases, cells_col, side="right") - 1
    keep_mask = np.isin(blocks, kept)
    if not bool(keep_mask.all()):
        notes.append(
            f"{region.kind} region '{region.name}': "
            f"{int((~keep_mask).sum())} entrie(s) lie in dropped cell blocks"
        )
    entries = entries[keep_mask]
    blocks = blocks[keep_mask]
    cells_col = cells_col[keep_mask]
    new_cells = new_bases[np.searchsorted(kept, blocks)] + (
        cells_col - old_bases[blocks]
    )
    if side:
        return np.column_stack([new_cells, entries[:, 1]]).astype(np.int64)
    return new_cells.astype(np.int64)


def _to_device_payload(
    mesh, *, float32=False, int32=False, zero_copy_only=False, op="to_dlpack"
):
    """The pure builder: a HOST :class:`DevicePayload`, plain numpy throughout."""
    if zero_copy_only and (float32 or int32):
        raise InteropError(
            f"meshio++: {op}: an explicit dtype downcast is a copy by "
            "definition; drop float32/int32 or zero_copy_only"
        )
    notes, derived, shared = [], [], {}

    points = _prepare_data_array(
        op, "points", mesh.points, notes, shared, zero_copy_only
    )
    if points is None:
        raise ValueError(
            f"meshio++: {op}: points cannot be exported "
            f"(dtype {np.asarray(mesh.points).dtype})"
        )
    points = _maybe_float32(points, "points", float32, notes, shared)

    cells, kept, dropped = [], [], []
    for b, cb in enumerate(mesh.cells):
        if not isinstance(cb.data, np.ndarray) or cb.data.dtype == object:
            dropped.append((b, cb.type))
            notes.append(
                f"cell block {b} ('{cb.type}') is ragged and cannot be a "
                "DLPack tensor; dropped"
            )
            continue
        label = f"cells[{b}] ('{cb.type}')"
        conn = cb.data
        if _is_native(conn):
            shared[label] = True
        else:
            _copy_note(op, label, "not C-contiguous", notes, shared, zero_copy_only)
            conn = np.ascontiguousarray(conn)
        if int32:
            conn = _checked_int_cast(conn, label, notes, shared, op=op)
        kept.append(b)
        cells.append(GpuCellBlock(cb.type, conn))
    kept = np.asarray(kept, dtype=np.int64)

    old_bases = block_bases(mesh.cells)
    new_bases = block_bases([mesh.cells[b] for b in kept])
    renumbered = len(kept) != len(mesh.cells)
    bases = new_bases
    if int32:
        bases = _checked_int_cast(bases, "block_bases", notes, shared, op=op)

    point_data = {}
    for name, raw in mesh.point_data.items():
        label = f"point_data '{name}'"
        a = _prepare_data_array(op, label, raw, notes, shared, zero_copy_only)
        if a is None:
            continue
        point_data[name] = _maybe_float32(a, label, float32, notes, shared)

    cell_data = {}
    for name in mesh.cell_data:
        got = _per_block_cell_data(
            mesh, name, kept, notes, shared, zero_copy_only, float32, op
        )
        if got is not None:
            cell_data[name] = got

    regions = []
    for region in getattr(mesh, "regions", []):
        label = f"region '{region.name}'"
        entries = np.asarray(region.entries, dtype=np.int64)
        if renumbered and region.kind in ("cell", "side"):
            entries = _remap_cell_entries(region, old_bases, kept, new_bases, notes)
            _derived_note(
                f"{label} entries",
                "renumbered into the payload's kept-block cell numbering",
                derived,
            )
            shared[label] = False
        else:
            shared[label] = entries is region.entries
        if int32:
            entries = _checked_int_cast(entries, label, notes, shared, op=op)
        regions.append(
            GpuRegion(region.name, region.kind, region.dim, region.tag, entries)
        )

    return DevicePayload(
        points=points,
        cells=tuple(cells),
        block_bases=bases,
        point_data=point_data,
        cell_data=cell_data,
        field_data=dict(mesh.field_data),
        regions=tuple(regions),
        device=(_KDLCPU, 0),
        dropped=tuple(dropped),
        shared=shared,
        notes=tuple(notes),
        derived=tuple(derived),
        refs=(),
    )


def _dlpack_device(obj):
    """``obj.__dlpack_device__()`` as a tuple, or None when unavailable."""
    try:
        return tuple(obj.__dlpack_device__())
    except Exception:
        return None


def _adopt_host(obj, label, notes):
    """One array back to host numpy — the consumption ladder.

    Everything except a device-resident array with no ``.get()`` works without
    CuPy installed: numpy passes through, any DLPack host (``kDLCPU``) exporter
    (a torch CPU tensor, a JAX host array) is adopted zero-copy via
    ``np.from_dlpack``, and a CuPy-like object copies itself back with its own
    ``.get()``. A device-resident CAI/DLPack object without ``.get()`` (a torch
    CUDA tensor, a Numba device array) needs CuPy to adopt it on device —
    zero-copy — before the one deliberate device→host copy.
    """
    if isinstance(obj, np.ndarray):
        return obj
    dev = _dlpack_device(obj)
    on_device = hasattr(obj, "__cuda_array_interface__") or (
        dev is not None and dev[0] != _KDLCPU
    )
    if on_device:
        notes.append(f"{label} was copied device→host")
        if hasattr(obj, "get"):
            return np.asarray(obj.get())
        cp = _require_cupy("from_cupy")
        return cp.asnumpy(cp.asarray(obj))
    if hasattr(obj, "__dlpack__") and hasattr(np, "from_dlpack"):
        try:
            return np.from_dlpack(obj)
        except (TypeError, RuntimeError, BufferError):
            pass
    return np.asarray(obj)


# --------------------------------------------------------------------------- #
# availability + the named install error                                      #
# --------------------------------------------------------------------------- #
def has_cupy() -> bool:
    """Whether CuPy is importable (which does NOT imply a device is present)."""
    return _importable("cupy")


def has_cuda_device() -> bool:
    """Whether CuPy is importable AND at least one device actually answers.

    CuPy imports fine on a machine with no GPU and fails only at first use, so
    checking the import alone is not an availability check.
    """
    if not _importable("cupy"):
        return False
    try:
        import cupy

        return int(cupy.cuda.runtime.getDeviceCount()) > 0
    except Exception:
        return False


def _require_cupy(op):
    """Import CuPy or raise naming the wheel recipe — not a pip extra.

    Deliberately not ``_interop._require``: that helper's message promises
    ``pip install meshioplusplus[extra]``, and no such extra exists (or should)
    for CuPy — the wheels are CUDA-version-specific.
    """
    if _importable("cupy"):
        import cupy

        return cupy
    raise ImportError(
        f"meshio++: {op}: cupy is not installed. There is deliberately no pip "
        "extra for it: CuPy wheels are CUDA-version-specific, so install the "
        "wheel matching your toolkit — e.g. `pip install cupy-cuda13x` "
        "(CUDA 13.x), `pip install cupy-cuda12x` (CUDA 12.x), "
        "`cupy-cuda11x` (CUDA 11.x), or a `cupy-rocm-*` wheel for ROCm. "
        "See doc/gpu.md."
    )


def _require_device(op):
    cp = _require_cupy(op)
    try:
        count = int(cp.cuda.runtime.getDeviceCount())
    except Exception:
        count = 0
    if count <= 0:
        raise RuntimeError(
            f"meshio++: {op}: cupy is installed but no CUDA device is "
            "available (cupy imports fine on a machine with no GPU and fails "
            "only at first use)"
        )
    return cp


# --------------------------------------------------------------------------- #
# public API                                                                  #
# --------------------------------------------------------------------------- #
def to_dlpack(mesh, *, float32=False, int32=False, zero_copy_only=False):
    """Export the mesh's arrays as a host :class:`DevicePayload`.

    Every array in the payload is a numpy ndarray and therefore natively
    implements ``__dlpack__`` / ``__dlpack_device__`` (reporting
    ``(kDLCPU, 0)``) — hand any of them to ``np.from_dlpack``,
    ``torch.from_dlpack``, ``jax.dlpack``, ``cupy.from_dlpack``, … directly.
    The payload is the container; each array speaks the protocol itself.

    This is a **host** buffer-sharing export, squarely inside the interop
    layer's ``zero_copy_only`` contract: canonical arrays are shared as-is, and
    anything that would copy either lands in ``notes`` (warned) or raises
    :class:`~meshioplusplus._interop.InteropError` under
    ``zero_copy_only=True``. ``float32=True`` / ``int32=True`` are explicit
    opt-in downcasts (recorded, range-checked for indices) and are refused
    outright in combination with ``zero_copy_only=True`` — a downcast is a copy
    by definition.
    """
    payload = _to_device_payload(
        mesh,
        float32=float32,
        int32=int32,
        zero_copy_only=zero_copy_only,
        op="to_dlpack",
    )
    _emit("to_dlpack", payload.notes)
    return payload


def to_cupy(mesh, *, float32=False, int32=False, pinned=False, stream=None):
    """Transfer the mesh's arrays to the current CUDA device via CuPy.

    Returns a :class:`DevicePayload` whose arrays are ``cupy.ndarray``s (each
    exposing both DLPack and ``__cuda_array_interface__``). **The host→device
    move is one bus transfer per array — it is never "zero-copy"**; what this
    call avoids is the file round-trip and any extra staging copy beyond it.
    There is no ``zero_copy_only`` parameter here on purpose: it would be a lie
    on a path that is 100% copies.

    ``pinned=True`` stages each array through page-locked (pinned) host memory,
    which is what enables true async DMA. With the default ``stream=None`` the
    call synchronizes before returning and the staging buffers are released;
    with a caller-supplied ``cupy.cuda.Stream`` the copies are enqueued
    asynchronously and the staging buffers are kept alive on ``payload.refs`` —
    **the caller must synchronize the stream before using the arrays or
    dropping the payload's refs**.

    Requires a working CUDA device, not just an importable CuPy — see
    :func:`has_cuda_device`.
    """
    cp = _require_device("to_cupy")
    host = _to_device_payload(
        mesh, float32=float32, int32=int32, zero_copy_only=False, op="to_cupy"
    )
    payload = _to_device(host, cp, pinned=pinned, stream=stream)
    _emit("to_cupy", payload.notes)
    return payload


def from_cupy(payload):
    """Rebuild a normal :class:`~meshioplusplus._mesh.Mesh` from a device payload.

    Accepts a :class:`DevicePayload` or any duck-typed object with the same
    attributes; each array may be a numpy/cupy array or anything exposing
    DLPack or ``__cuda_array_interface__`` (torch tensors, Numba device
    arrays). Device-resident arrays are brought back with **one deliberate
    device→host copy each** (recorded in the warned notes — that is the
    documented cost, not a trick); host DLPack exporters are adopted zero-copy.
    The result is an ordinary ``Mesh``: every existing writer, operation and
    binding works on it unchanged.
    """
    from ._mesh import Mesh
    from ._regions import Region

    notes = []
    points = _adopt_host(payload.points, "points", notes)
    cells = []
    for b, blk in enumerate(getattr(payload, "cells", ())):
        ctype, data = blk[0], blk[1]
        cells.append((ctype, _adopt_host(data, f"cells[{b}] ('{ctype}')", notes)))
    point_data = {
        name: _adopt_host(a, f"point_data '{name}'", notes)
        for name, a in dict(getattr(payload, "point_data", {})).items()
    }
    cell_data = {
        name: [
            _adopt_host(a, f"cell_data '{name}' block {b}", notes)
            for b, a in enumerate(blocks)
        ]
        for name, blocks in dict(getattr(payload, "cell_data", {})).items()
    }
    regions = [
        Region(
            r[0],
            r[1],
            _adopt_host(r[4], f"region '{r[0]}'", notes),
            dim=int(r[2]),
            tag=int(r[3]),
        )
        for r in getattr(payload, "regions", ())
    ]
    mesh = Mesh(
        points,
        cells,
        point_data=point_data,
        cell_data=cell_data,
        field_data=dict(getattr(payload, "field_data", {})),
        regions=regions,
    )
    _emit("from_cupy", notes)
    return mesh


# --------------------------------------------------------------------------- #
# device transfer (CuPy only past this line)                                  #
# --------------------------------------------------------------------------- #
def _pinned_transfer(host, cp, stream):
    """host → pinned staging → device. Returns (device array, staging or None)."""
    if host.size == 0:
        return cp.asarray(host), None
    mem = cp.cuda.alloc_pinned_memory(host.nbytes)
    staging = np.frombuffer(mem, host.dtype, host.size).reshape(host.shape)
    staging[...] = host
    dev = cp.empty(host.shape, host.dtype)
    dev.set(staging, stream=stream)
    return dev, staging


def _to_device(host, cp, *, pinned, stream):
    """Move a host payload's arrays to the device; metadata carries over."""
    staging_refs = []

    def move(a):
        if pinned:
            dev, staging = _pinned_transfer(a, cp, stream)
            if stream is not None and staging is not None:
                staging_refs.append(staging)
            return dev
        if stream is not None:
            with stream:
                return cp.asarray(a)
        return cp.asarray(a)

    points = move(host.points)
    cells = tuple(GpuCellBlock(blk.type, move(blk.data)) for blk in host.cells)
    bases = move(host.block_bases)
    point_data = {name: move(a) for name, a in host.point_data.items()}
    cell_data = {
        name: tuple(move(a) for a in blocks) for name, blocks in host.cell_data.items()
    }
    regions = tuple(
        GpuRegion(r.name, r.kind, r.dim, r.tag, move(r.entries)) for r in host.regions
    )
    if pinned and stream is None:
        # dev.set() on pinned memory is asynchronous on the current stream;
        # synchronize so the staging buffers can be released here and now.
        cp.cuda.get_current_stream().synchronize()
    derived = host.derived + (
        "arrays were copied host→device (the one bus transfer this call "
        "exists to perform)",
    )
    return DevicePayload(
        points=points,
        cells=cells,
        block_bases=bases,
        point_data=point_data,
        cell_data=cell_data,
        field_data=host.field_data,
        regions=regions,
        device=_dlpack_device(points) or (2, 0),
        dropped=host.dropped,
        shared=host.shared,
        notes=host.notes,
        derived=derived,
        refs=tuple(staging_refs),
    )


# --------------------------------------------------------------------------- #
# PyTorch / JAX tensor handoff                                                #
# --------------------------------------------------------------------------- #
def has_torch() -> bool:
    """Whether PyTorch is importable. There is deliberately no pip extra."""
    return _importable("torch")


def has_jax() -> bool:
    """Whether JAX is importable. There is deliberately no pip extra."""
    return _importable("jax")


def _require_framework(
    op,
    module,
    hint,
    doc="doc/ml.md",
    qualifier="(pick the wheel matching your accelerator)",
):
    """Import a heavyweight optional dependency or raise naming the install command.

    The CuPy precedent, not ``_interop._require``: there is deliberately no
    ``meshioplusplus[torch]``/``[jax]`` extra — torch's default Linux wheel
    bundles CUDA at ~900 MB and JAX's accelerator story lives in its own
    extras (``jax[cuda12]``, …), so an extra pinning either would surprise far
    more people than it would help (the Open3D ~400 MB-wheel reasoning).
    ``doc`` names the page the error points at — the PhysicsNeMo adapter
    reuses this rather than growing a second no-extra error shape.

    ``qualifier`` is the parenthetical after the install hint. It defaults to
    the accelerator wording every ML caller wants, so their messages are
    unchanged; the Blender bridge overrides it, since ``bpy``'s constraint is a
    CPython pin rather than a GPU.
    """
    if _importable(module):
        return __import__(module)
    raise ImportError(
        f"meshio++: {op}: {module} is not installed. There is deliberately no "
        f"pip extra for it; install it directly with `{hint}` {qualifier}. "
        f"See {doc}."
    )


def _adopted_payload(host, adopt, extra_derived):
    """A ``DevicePayload`` with every array passed through ``adopt``."""
    points = adopt(host.points)
    cells = tuple(GpuCellBlock(b.type, adopt(b.data)) for b in host.cells)
    payload = DevicePayload(
        points=points,
        cells=cells,
        block_bases=adopt(host.block_bases),
        point_data={name: adopt(a) for name, a in host.point_data.items()},
        cell_data={
            name: tuple(adopt(a) for a in blocks)
            for name, blocks in host.cell_data.items()
        },
        field_data=host.field_data,
        regions=tuple(
            GpuRegion(r.name, r.kind, r.dim, r.tag, adopt(r.entries))
            for r in host.regions
        ),
        device=_dlpack_device(points) or host.device,
        dropped=host.dropped,
        shared=host.shared,
        notes=host.notes,
        derived=host.derived + extra_derived,
        refs=(),
    )
    return payload


def to_torch(mesh, *, device=None, float32=False, int32=False):
    """The mesh's arrays as PyTorch tensors, adopted over DLPack.

    Returns a :class:`DevicePayload` whose arrays are ``torch.Tensor``s. With
    ``device=None`` (default) every tensor **adopts the host numpy buffer
    zero-copy** via ``torch.from_dlpack`` — nothing is transferred, and the
    tensors share memory with the mesh. Passing ``device=`` (``"cuda"``, a
    ``torch.device``, …) then moves each tensor with ``.to(device)`` — per the
    GPU honesty rule that is **one bus transfer per array**, recorded in the
    payload's ``derived`` list exactly like :func:`to_cupy`'s transfer.

    ``float32=True`` / ``int32=True`` are the same explicit opt-in downcasts
    :func:`to_dlpack` takes (applied on the host side, before adoption).

    Raises
    ------
    ImportError
        when torch is not installed (naming ``pip install torch`` — there is
        deliberately no pip extra, see doc/ml.md).
    """
    torch = _require_framework("to_torch", "torch", "pip install torch")

    host = _to_device_payload(
        mesh, float32=float32, int32=int32, zero_copy_only=False, op="to_torch"
    )

    if device is None:
        extra = ("tensors adopt the host buffers via DLPack — zero-copy",)

        def adopt(a):
            return torch.from_dlpack(a)

    else:
        extra = (
            f"tensors were moved to {device!r} with .to() (the one bus "
            "transfer this call exists to perform)",
        )

        def adopt(a):
            return torch.from_dlpack(a).to(device)

    payload = _adopted_payload(host, adopt, extra)
    _emit("to_torch", payload.notes)
    return payload


def to_jax(mesh, *, float32=False, int32=False):
    """The mesh's arrays as JAX arrays, adopted over DLPack.

    Returns a :class:`DevicePayload` whose arrays are ``jax.Array``s, placed
    on JAX's **default device** — on a GPU machine that placement is itself a
    transfer, which is JAX's committed behaviour rather than this function's
    choice. There is no ``device=`` parameter: device placement belongs to
    JAX (``jax.default_device``, ``jax.device_put``).

    One JAX-specific honesty note: under JAX's default configuration
    (``jax_enable_x64=False``) 64-bit arrays cannot exist, so meshio++'s
    canonical float64/int64 arrays are adopted via a fallback that lets JAX
    apply its own 32-bit demotion — a copy, counted and warned once. Pass
    ``float32=True``/``int32=True`` to make the narrowing explicit on the
    meshio++ side, or enable x64 in JAX for lossless adoption.

    Raises
    ------
    ImportError
        when jax is not installed (naming ``pip install jax`` — there is
        deliberately no pip extra, see doc/ml.md).
    """
    _require_framework("to_jax", "jax", "pip install jax")
    from jax import dlpack as jax_dlpack
    from jax import numpy as jnp

    host = _to_device_payload(
        mesh, float32=float32, int32=int32, zero_copy_only=False, op="to_jax"
    )

    fallbacks = []

    def adopt(a):
        try:
            return jax_dlpack.from_dlpack(a)
        except Exception:
            fallbacks.append(1)
            return jnp.asarray(a)

    payload = _adopted_payload(
        host,
        adopt,
        ("arrays were placed on JAX's default device (JAX's own placement)",),
    )
    if fallbacks:
        payload = DevicePayload(
            **{
                **{f: getattr(payload, f) for f in DevicePayload.__dataclass_fields__},
                "notes": payload.notes
                + (
                    f"{len(fallbacks)} array(s) could not be adopted over "
                    "DLPack and were converted by jnp.asarray (a copy; under "
                    "jax_enable_x64=False JAX demotes 64-bit dtypes)",
                ),
            }
        )
    _emit("to_jax", payload.notes)
    return payload

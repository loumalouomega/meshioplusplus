"""Tests for the GPU handoff layer (DLPack export, CuPy transfer).

Two halves, deliberately — the ``test_interop.py`` pattern.

The **always-on** half runs in the default CI matrix on machines with no GPU
and no CuPy: the pure payload builder, and the DLPack **host** round-trip —
which genuinely exercises the export machinery, because DLPack covers host
arrays (``kDLCPU``) and numpy ≥ 1.23 speaks it natively.

The **gated** half needs a real CUDA device. ``pytest.importorskip("cupy")``
alone is NOT enough: CuPy imports fine on a machine with no GPU and fails only
at first use, so the gate also checks ``getDeviceCount() > 0``. Public CI has
no GPU runner — this half is documented as not covered there (doc/gpu.md).
"""

from __future__ import annotations

import gc
from types import SimpleNamespace

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _gpu
from meshioplusplus._interop import InteropError
from meshioplusplus._mesh import Mesh
from meshioplusplus._regions import Region, block_bases

_HAS_DLPACK = hasattr(np, "from_dlpack") and hasattr(np.ndarray, "__dlpack__")
needs_dlpack = pytest.mark.skipif(not _HAS_DLPACK, reason="numpy < 1.23: no DLPack")


def _device_available():
    try:
        import cupy

        return int(cupy.cuda.runtime.getDeviceCount()) > 0
    except Exception:
        return False


needs_gpu = pytest.mark.skipif(
    not _device_available(), reason="no CUDA device (cupy missing or no GPU)"
)


# --------------------------------------------------------------------------- #
# fixtures                                                                    #
# --------------------------------------------------------------------------- #
def _mixed_mesh():
    """triangle + quad + tetra, with point/cell/field data and all three region kinds."""
    points = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    cells = [
        ("triangle", np.array([[0, 1, 2]], dtype=np.int64)),
        ("quad", np.array([[0, 1, 2, 3]], dtype=np.int64)),
        ("tetra", np.array([[0, 1, 2, 4], [0, 2, 3, 4]], dtype=np.int64)),
    ]
    return Mesh(
        points,
        cells,
        point_data={
            "T": np.arange(5, dtype=np.float64),
            "v": np.arange(15, dtype=np.float64).reshape(5, 3),
        },
        cell_data={
            "mat": [
                np.array([10], dtype=np.int64),
                np.array([20], dtype=np.int64),
                np.array([30, 31], dtype=np.int64),
            ]
        },
        field_data={"time": np.array([1.5])},
        regions=[
            Region("inlet", "point", [0, 1], dim=0, tag=5),
            Region("steel", "cell", [1, 3], dim=3, tag=7),
            Region("wall", "side", [[0, 1], [2, 0]], dim=2, tag=9),
        ],
    )


def _triangle_mesh():
    """A single-block, int64/float64 triangle mesh — the shareable-everything case."""
    points = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]]
    )
    return Mesh(
        points,
        [("triangle", np.array([[0, 1, 2], [0, 2, 3]], dtype=np.int64))],
        point_data={
            "T": np.arange(4, dtype=np.float64),
            "v": np.arange(12, dtype=np.float64).reshape(4, 3),
        },
        cell_data={"mat": [np.array([1, 2], dtype=np.int64)]},
        regions=[
            Region("inlet", "point", [0, 1], dim=0, tag=5),
            Region("skin", "cell", [1], dim=2, tag=9),
        ],
    )


def _ragged_mesh():
    """triangle + ragged polygon + quad; regions reference kept AND dropped cells.

    Global block-major numbering: triangle cell = 0, polygon cells = 1 and 2
    (ragged, dropped by the payload), quad cell = 3. After the drop the kept
    numbering is triangle = 0, quad = 1.
    """
    points = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.5, 1.5, 0.0],
        ]
    )
    cells = [
        ("triangle", np.array([[0, 1, 2]], dtype=np.int64)),
        ("polygon", [[0, 1, 2, 3, 4], [0, 1, 2]]),
        ("quad", np.array([[0, 1, 2, 3]], dtype=np.int64)),
    ]
    return Mesh(
        points,
        cells,
        cell_data={
            "mat": [
                np.array([10], dtype=np.int64),
                np.array([20, 21], dtype=np.int64),
                np.array([30], dtype=np.int64),
            ]
        },
        regions=[
            Region("mixed", "cell", [0, 2, 3], dim=2, tag=1),
            Region("wall", "side", [[3, 0]], dim=1, tag=2),
        ],
    )


# --------------------------------------------------------------------------- #
# always-on: the pure payload builder                                         #
# --------------------------------------------------------------------------- #
def test_device_payload_mixed_types():
    mesh = _mixed_mesh()
    payload = _gpu._to_device_payload(mesh)

    assert payload.device == (1, 0)  # kDLCPU: a host payload
    assert payload.dropped == ()
    assert [blk.type for blk in payload.cells] == ["triangle", "quad", "tetra"]
    # canonical dtypes, untouched by default
    assert payload.points.dtype == np.float64
    assert all(blk.data.dtype == np.int64 for blk in payload.cells)
    # block-major numbering comes from _regions.block_bases, never re-derived
    assert np.array_equal(payload.block_bases, block_bases(mesh.cells))
    # everything canonical is shared, not copied
    assert np.shares_memory(payload.points, mesh.points)
    assert np.shares_memory(payload.cells[2].data, mesh.cells[2].data)
    assert np.shares_memory(payload.point_data["v"], mesh.point_data["v"])
    assert payload.notes == ()
    # field_data is a host passthrough
    assert np.array_equal(payload.field_data["time"], [1.5])


def test_device_payload_cell_data_is_per_block():
    mesh = _mixed_mesh()
    payload = _gpu._to_device_payload(mesh)

    blocks = payload.cell_data["mat"]
    assert len(blocks) == len(payload.cells)
    for got, expected in zip(blocks, mesh.cell_data["mat"]):
        assert np.array_equal(got, expected)
        assert np.shares_memory(got, expected)


def test_device_payload_regions_are_index_arrays():
    mesh = _mixed_mesh()
    payload = _gpu._to_device_payload(mesh)

    by_name = {r.name: r for r in payload.regions}
    assert set(by_name) == {"inlet", "steel", "wall"}
    inlet, steel, wall = by_name["inlet"], by_name["steel"], by_name["wall"]
    assert (inlet.kind, inlet.dim, inlet.tag) == ("point", 0, 5)
    assert np.array_equal(inlet.entries, [0, 1])
    assert (steel.kind, steel.dim, steel.tag) == ("cell", 3, 7)
    assert np.array_equal(steel.entries, [1, 3])
    # Side (cell, facet) pairs ride along — no other interop target keeps them.
    assert wall.kind == "side"
    assert wall.entries.shape == (2, 2)
    assert np.array_equal(wall.entries, [[0, 1], [2, 0]])
    # index arrays are shared with the mesh's own regions, not rebuilt
    assert np.shares_memory(steel.entries, mesh.regions[1].entries)


def test_device_payload_drops_ragged_blocks_and_remaps_regions():
    mesh = _ragged_mesh()
    payload = _gpu._to_device_payload(mesh)

    assert payload.dropped == ((1, "polygon"),)
    assert any("ragged" in n for n in payload.notes)
    assert [blk.type for blk in payload.cells] == ["triangle", "quad"]
    # bases are over the payload's OWN cells (1 triangle + 1 quad)
    assert np.array_equal(payload.block_bases, [0, 1, 2])
    # cell_data follows the kept blocks
    assert [a.tolist() for a in payload.cell_data["mat"]] == [[10], [30]]

    by_name = {r.name: r for r in payload.regions}
    # global cells [0, 2, 3]: 0 kept -> 0, 2 in the dropped block, 3 -> 1
    assert np.array_equal(by_name["mixed"].entries, [0, 1])
    assert any("dropped cell blocks" in n for n in payload.notes)
    # side entry (3, 0): cell renumbered, facet column untouched
    assert np.array_equal(by_name["wall"].entries, [[1, 0]])
    assert any("renumbered" in d for d in payload.derived)


def test_device_payload_float32_is_opt_in_and_noted():
    mesh = _mixed_mesh()
    assert _gpu._to_device_payload(mesh).points.dtype == np.float64  # default

    payload = _gpu._to_device_payload(mesh, float32=True)
    assert payload.points.dtype == np.float32
    assert payload.point_data["T"].dtype == np.float32
    # integer arrays are untouched by float32
    assert payload.cells[0].data.dtype == np.int64
    assert payload.cell_data["mat"][0].dtype == np.int64
    assert any("float32=True" in n for n in payload.notes)
    assert payload.shared["points"] is False


def test_device_payload_int32_range_check():
    mesh = _mixed_mesh()
    payload = _gpu._to_device_payload(mesh, int32=True)
    assert all(blk.data.dtype == np.int32 for blk in payload.cells)
    assert payload.block_bases.dtype == np.int32
    assert all(r.entries.dtype == np.int32 for r in payload.regions)
    assert any("int32=True" in n for n in payload.notes)

    # an index outside the int32 range refuses by name rather than wrapping
    big = Mesh(
        np.zeros((3, 3)),
        [("triangle", np.array([[0, 1, 2**31]], dtype=np.int64))],
    )
    with pytest.raises(InteropError, match=r"cells\[0\].*int32"):
        _gpu._to_device_payload(big, int32=True)


def test_device_payload_data_policy_matches_interop():
    mesh = _triangle_mesh()
    mesh.point_data["mask"] = np.array([True, False, True, False])
    mesh.point_data["obj"] = np.array([None, 1, "x", 2.0], dtype=object)
    payload = _gpu._to_device_payload(mesh)

    assert payload.point_data["mask"].dtype == np.uint8  # bool widened, noted
    assert any("uint8" in n for n in payload.notes)
    assert "obj" not in payload.point_data  # dropped, noted
    assert any("obj" in n and "dropped" in n for n in payload.notes)


# --------------------------------------------------------------------------- #
# always-on: DLPack host round trip                                           #
# --------------------------------------------------------------------------- #
@needs_dlpack
def test_to_dlpack_arrays_speak_dlpack():
    payload = meshioplusplus.to_dlpack(_mixed_mesh())
    for a in (payload.points, payload.cells[0].data, payload.point_data["T"]):
        assert hasattr(a, "__dlpack__")
        assert tuple(a.__dlpack_device__()) == (1, 0)  # kDLCPU


@needs_dlpack
def test_to_dlpack_host_round_trip_shares():
    mesh = _mixed_mesh()
    payload = meshioplusplus.to_dlpack(mesh)
    adopted = np.from_dlpack(payload.points)
    assert np.array_equal(adopted, mesh.points)
    # host adoption is genuinely zero-copy — measured, not assumed
    assert np.shares_memory(adopted, mesh.points)


def test_to_dlpack_zero_copy_only():
    mesh = _triangle_mesh()
    mesh.point_data["T"] = np.arange(8, dtype=np.float64)[::2]  # non-contiguous
    with pytest.raises(InteropError, match="point_data 'T'"):
        meshioplusplus.to_dlpack(mesh, zero_copy_only=True)

    # an explicit downcast is a copy by definition
    with pytest.raises(InteropError, match="downcast is a copy"):
        meshioplusplus.to_dlpack(_triangle_mesh(), float32=True, zero_copy_only=True)


@needs_dlpack
def test_from_cupy_round_trips_a_host_payload():
    mesh = _mixed_mesh()
    back = meshioplusplus.from_cupy(meshioplusplus.to_dlpack(mesh))
    assert meshioplusplus.meshes_equal(mesh, back)
    assert back.regions == mesh.regions  # incl. dim/tag and the side pairs
    assert np.array_equal(back.field_data["time"], mesh.field_data["time"])


class _DLPackOnly:
    """An object exposing ONLY the DLPack protocol (delegating to numpy)."""

    def __init__(self, a):
        self._a = a

    def __dlpack__(self, **kwargs):
        return self._a.__dlpack__(**kwargs)

    def __dlpack_device__(self):
        return self._a.__dlpack_device__()


@needs_dlpack
def test_from_cupy_adopts_foreign_dlpack_objects():
    mesh = _triangle_mesh()
    payload = SimpleNamespace(
        points=_DLPackOnly(mesh.points),
        cells=(("triangle", _DLPackOnly(mesh.cells[0].data)),),
        point_data={"T": _DLPackOnly(mesh.point_data["T"])},
        cell_data={},
        field_data={},
        regions=(),
    )
    back = meshioplusplus.from_cupy(payload)
    # a kDLCPU exporter is adopted zero-copy, no CuPy involved
    assert np.shares_memory(back.points, mesh.points)
    assert np.shares_memory(back.point_data["T"], mesh.point_data["T"])
    assert np.array_equal(back.cells[0].data, mesh.cells[0].data)


class _FakeDeviceArray:
    """Claims to live on a CUDA device (CAI present), offers no way back."""

    __cuda_array_interface__ = {"shape": (1,), "typestr": "<f8", "data": (0, False)}


def test_from_cupy_device_array_without_cupy_names_the_wheel(monkeypatch):
    monkeypatch.setattr(_gpu, "_importable", lambda module: False)
    payload = SimpleNamespace(
        points=_FakeDeviceArray(),
        cells=(),
        point_data={},
        cell_data={},
        field_data={},
        regions=(),
    )
    with pytest.raises(ImportError, match="cupy-cuda12x") as excinfo:
        meshioplusplus.from_cupy(payload)
    # the wheel recipe, deliberately NOT a broken pip extra
    assert "meshioplusplus[" not in str(excinfo.value)


def test_to_cupy_without_cupy_names_the_wheel(monkeypatch):
    monkeypatch.setattr(_gpu, "_importable", lambda module: False)
    with pytest.raises(ImportError, match="cupy-cuda12x") as excinfo:
        meshioplusplus.to_cupy(_triangle_mesh())
    assert "meshioplusplus[" not in str(excinfo.value)


def test_predicates_are_booleans():
    assert isinstance(meshioplusplus.has_cupy(), bool)
    assert isinstance(meshioplusplus.has_cuda_device(), bool)
    assert isinstance(meshioplusplus.has_torch(), bool)
    assert isinstance(meshioplusplus.has_jax(), bool)
    if not meshioplusplus.has_cupy():
        assert meshioplusplus.has_cuda_device() is False


def test_public_api_is_exported():
    for name in (
        "to_dlpack",
        "to_cupy",
        "from_cupy",
        "has_cupy",
        "has_cuda_device",
        "to_torch",
        "to_jax",
        "has_torch",
        "has_jax",
    ):
        assert name in meshioplusplus.__all__
        assert hasattr(meshioplusplus, name)


# --------------------------------------------------------------------------- #
# PyTorch / JAX tensor handoff                                                #
# --------------------------------------------------------------------------- #
def test_require_framework_names_the_install_command(monkeypatch):
    monkeypatch.setattr(_gpu, "_importable", lambda module: False)
    with pytest.raises(ImportError, match="pip install torch") as excinfo:
        meshioplusplus.to_torch(_triangle_mesh())
    # No pip extra exists for torch/jax, deliberately — the CuPy precedent.
    assert "meshioplusplus[" not in str(excinfo.value)
    with pytest.raises(ImportError, match="pip install jax"):
        meshioplusplus.to_jax(_triangle_mesh())


def test_to_torch_adopts_host_buffers():
    torch = pytest.importorskip("torch")
    mesh = _triangle_mesh()
    payload = meshioplusplus.to_torch(mesh)

    assert isinstance(payload.points, torch.Tensor)
    assert np.array_equal(payload.points.numpy(), mesh.points)
    # Host adoption is genuinely zero-copy — measured, not assumed.
    assert np.shares_memory(payload.points.numpy(), mesh.points)
    assert np.shares_memory(payload.point_data["T"].numpy(), mesh.point_data["T"])
    assert payload.cells[0].type == "triangle"
    assert np.array_equal(payload.cells[0].data.numpy(), mesh.cells[0].data)


def test_to_jax_adopts_arrays():
    jnp = pytest.importorskip("jax.numpy")
    mesh = _triangle_mesh()
    payload = meshioplusplus.to_jax(mesh)

    assert isinstance(payload.points, jnp.ndarray)
    assert np.allclose(np.asarray(payload.points), mesh.points)
    assert np.allclose(np.asarray(payload.point_data["T"]), mesh.point_data["T"])


@needs_dlpack
def test_payload_outlives_the_source_mesh(tmp_path):
    """Readers hand back capsule-backed numpy owned by the C++ core — the
    payload must keep it alive, or this is a use-after-free rather than a
    wrong answer."""
    p = tmp_path / "m.vtu"
    meshioplusplus.write(p, _triangle_mesh())

    def build():
        mesh = meshioplusplus.read(p)
        return meshioplusplus.to_dlpack(mesh), np.array(mesh.points)

    payload, expected = build()
    gc.collect()
    for _ in range(3):
        churn = np.random.rand(1_000_000)
        del churn
    gc.collect()

    assert np.array_equal(payload.points, expected)
    assert np.array_equal(np.from_dlpack(payload.points), expected)


# --------------------------------------------------------------------------- #
# gated: a real CUDA device (NOT covered by public CI — doc/gpu.md)           #
# --------------------------------------------------------------------------- #
@needs_gpu
def test_cupy_payload_is_on_device():
    import cupy as cp

    payload = meshioplusplus.to_cupy(_mixed_mesh())
    assert payload.device[0] != 1  # not kDLCPU: CUDA (2) or ROCm (10)
    for a in (payload.points, payload.cells[0].data, payload.point_data["T"]):
        assert isinstance(a, cp.ndarray)
        # the arrays speak both protocols: DLPack and the older CAI
        assert hasattr(a, "__dlpack__")
        assert hasattr(a, "__cuda_array_interface__")


@needs_gpu
def test_cupy_round_trip(tmp_path):
    mesh = _mixed_mesh()
    back = meshioplusplus.from_cupy(meshioplusplus.to_cupy(mesh))
    assert meshioplusplus.meshes_equal(mesh, back)
    assert back.regions == mesh.regions
    # end-to-end: device arrays feed every existing writer unchanged
    meshioplusplus.write(tmp_path / "roundtrip.vtu", back)
    again = meshioplusplus.read(tmp_path / "roundtrip.vtu")
    assert meshioplusplus.meshes_equal(mesh, again)


@needs_gpu
def test_cupy_dtype_opt_ins():
    import cupy as cp

    payload = meshioplusplus.to_cupy(_mixed_mesh(), float32=True, int32=True)
    assert payload.points.dtype == cp.float32
    assert payload.cells[0].data.dtype == cp.int32
    assert any("float32=True" in n for n in payload.notes)
    assert any("int32=True" in n for n in payload.notes)


@needs_gpu
def test_cupy_pinned_equivalence():
    import cupy as cp

    mesh = _mixed_mesh()
    plain = meshioplusplus.to_cupy(mesh)
    pinned = meshioplusplus.to_cupy(mesh, pinned=True)
    stream = cp.cuda.Stream(non_blocking=True)
    asynchronous = meshioplusplus.to_cupy(mesh, pinned=True, stream=stream)
    assert asynchronous.refs  # staging buffers stay alive until the sync
    stream.synchronize()

    for a, b, c in (
        (plain.points, pinned.points, asynchronous.points),
        (plain.point_data["v"], pinned.point_data["v"], asynchronous.point_data["v"]),
        (plain.cells[2].data, pinned.cells[2].data, asynchronous.cells[2].data),
    ):
        assert np.array_equal(cp.asnumpy(a), cp.asnumpy(b))
        assert np.array_equal(cp.asnumpy(a), cp.asnumpy(c))


@needs_gpu
def test_from_cupy_consumes_raw_cai_object():
    import cupy as cp

    class _CaiOnly:
        def __init__(self, a):
            self.__cuda_array_interface__ = a.__cuda_array_interface__
            self._keep = a

    mesh = _triangle_mesh()
    device = meshioplusplus.to_cupy(mesh)
    payload = SimpleNamespace(
        points=_CaiOnly(device.points),
        cells=(("triangle", _CaiOnly(device.cells[0].data)),),
        point_data={},
        cell_data={},
        field_data={},
        regions=(),
    )
    back = meshioplusplus.from_cupy(payload)
    assert np.array_equal(back.points, cp.asnumpy(device.points))
    assert np.array_equal(back.cells[0].data, mesh.cells[0].data)


@needs_gpu
def test_has_cuda_device_is_true_here():
    assert meshioplusplus.has_cupy()
    assert meshioplusplus.has_cuda_device()

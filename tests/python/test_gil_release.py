"""The C++ core releases the GIL around its reads, writes and operations
(roadmap 3.4.1).

Two gates: a thread pool gives the same results as serial calls, and a
heartbeat thread keeps ticking during one long native read or operation --
which it cannot do while the call holds the GIL.
"""

import os
import threading
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pytest

import meshioplusplus

from . import helpers

try:
    from meshioplusplus import _core
except ImportError:
    _core = None

needs_core = pytest.mark.skipif(_core is None, reason="needs the C++ core")

_WORKERS = 8


def _assert_same_mesh(a, b):
    np.testing.assert_array_equal(a.points, b.points)
    assert len(a.cells) == len(b.cells)
    for block_a, block_b in zip(a.cells, b.cells):
        assert block_a.type == block_b.type
        np.testing.assert_array_equal(block_a.data, block_b.data)
    assert sorted(a.point_data) == sorted(b.point_data)
    for name in a.point_data:
        np.testing.assert_array_equal(a.point_data[name], b.point_data[name])


def _lattice(n, seed):
    mesh = meshioplusplus.grid((n, n, n), spacing=(1.0, 0.5, 0.25))
    rng = np.random.default_rng(seed)
    mesh.point_data["u"] = rng.random(len(mesh.points))
    return mesh


# (extension, mesh factory) pairs over formats whose native paths release the
# GIL; vtkhdf joins when the build has HDF5, to show the paths that keep it
# are still correct when called from many threads.
def _cases():
    cases = [
        (".vtu", lambda i: _lattice(6 + i % 3, i)),
        (".vtk", lambda i: _lattice(5 + i % 3, i)),
        (".msh", lambda i: _lattice(4 + i % 3, i)),
        (".mdpa", lambda i: _lattice(4 + i % 2, i)),
        (".ply", lambda i: helpers.tri_mesh),
        (".stl", lambda i: helpers.tri_mesh),
        (".obj", lambda i: helpers.tri_mesh),
    ]
    if _core is not None and _core.__has_hdf5__:
        cases.append((".vtkhdf", lambda i: _lattice(4 + i % 2, i)))
    return cases


@needs_core
@pytest.mark.parametrize("extension,factory", _cases())
def test_thread_pool_round_trips_equal_serial(tmp_path, extension, factory):
    count = 2 * _WORKERS
    meshes = [factory(i) for i in range(count)]

    def round_trip(i, tag):
        path = str(tmp_path / f"{tag}_{i}{extension}")
        meshioplusplus.write(path, meshes[i])
        return meshioplusplus.read(path)

    serial = [round_trip(i, "serial") for i in range(count)]
    with ThreadPoolExecutor(max_workers=_WORKERS) as pool:
        threaded = list(pool.map(lambda i: round_trip(i, "threaded"), range(count)))
    for a, b in zip(serial, threaded):
        _assert_same_mesh(a, b)
    if extension in (".vtu", ".vtk", ".msh", ".ply", ".stl", ".obj"):
        # Writers are deterministic, so the files themselves agree too.
        for i in range(count):
            with open(tmp_path / f"serial_{i}{extension}", "rb") as f:
                expected = f.read()
            with open(tmp_path / f"threaded_{i}{extension}", "rb") as f:
                assert f.read() == expected


@needs_core
def test_thread_pool_operations_equal_serial():
    meshes = [_lattice(5 + i % 4, i) for i in range(2 * _WORKERS)]

    def run(mesh):
        skin = meshioplusplus.extract_skin(mesh)
        fine = meshioplusplus.refine(mesh)
        moved = meshioplusplus.transform(mesh, translate=(1.0, 2.0, 3.0))
        quality = meshioplusplus.compute_quality(mesh)
        return skin, fine, moved, quality

    serial = [run(m) for m in meshes]
    with ThreadPoolExecutor(max_workers=_WORKERS) as pool:
        threaded = list(pool.map(run, meshes))
    for expected, got in zip(serial, threaded):
        for a, b in zip(expected[:3], got[:3]):
            _assert_same_mesh(a, b)
        assert expected[3]["num_cells"] == got[3]["num_cells"]
        for name, arrays in expected[3]["cell_arrays"].items():
            for a, b in zip(arrays, got[3]["cell_arrays"][name]):
                np.testing.assert_array_equal(a, b)


def _ticks_during(call):
    """Run ``call`` while a daemon thread counts 1 ms naps; return
    ``(ticks, elapsed)``. A call that holds the GIL starves the counter."""
    ticks = 0
    started = threading.Event()
    stop = threading.Event()

    def heartbeat():
        nonlocal ticks
        started.set()
        while not stop.is_set():
            time.sleep(0.001)
            ticks += 1

    thread = threading.Thread(target=heartbeat, daemon=True)
    thread.start()
    started.wait()
    time.sleep(0.01)
    before = ticks
    t0 = time.perf_counter()
    call()
    elapsed = time.perf_counter() - t0
    during = ticks - before
    stop.set()
    thread.join()
    return during, elapsed


def _assert_heartbeat(call):
    during, elapsed = _ticks_during(call)
    if elapsed < 0.05:
        pytest.skip(f"the call took {elapsed * 1e3:.0f} ms, too short to judge")
    # Unimpeded, the counter ticks about once per (1 ms nap + wake-up). Holding
    # the GIL leaves it at most a tick or two; a fifth of the ideal rate is far
    # from both, and robust to a loaded machine.
    expected = elapsed / 0.002
    assert during >= max(5, 0.2 * expected), (
        f"heartbeat ticked {during} times in {elapsed * 1e3:.0f} ms: "
        "the native call held the GIL"
    )


@needs_core
def test_heartbeat_runs_during_a_native_read(tmp_path):
    path = str(tmp_path / "big.vtu")
    mesh = _lattice(int(os.environ.get("MESHIOPLUSPLUS_GIL_TEST_N", "70")), 0)
    _core.vtu_write(path, mesh, True, False)
    _assert_heartbeat(lambda: _core.vtu_read(path))


@needs_core
def test_heartbeat_runs_during_a_native_operation():
    mesh = _lattice(int(os.environ.get("MESHIOPLUSPLUS_GIL_TEST_N", "70")), 0)
    _assert_heartbeat(lambda: _core.compute_quality(mesh))

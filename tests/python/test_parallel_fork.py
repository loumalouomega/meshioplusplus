"""The THREADS parallel backend across ``fork`` (roadmap 5.2).

The Linux wheels build the core's ``std::thread`` pool. A process that used it
and then forks (``multiprocessing``'s default start method on Linux) hands the
child a copy of the pool's memory without its threads; the core rebuilds the
pool in the child. Skipped unless the core was built with THREADS.
"""

import multiprocessing
import sys

import numpy as np
import pytest

import meshioplusplus

pytestmark = [
    pytest.mark.skipif(
        meshioplusplus._core.__parallel_backend__ != "threads",
        reason="needs the THREADS parallel backend",
    ),
    pytest.mark.skipif(sys.platform == "win32", reason="fork is POSIX-only"),
]


def _cube(n):
    """A structured grid of hexahedra, large enough for parallel loops."""
    x = np.linspace(0.0, 1.0, n + 1)
    pts = np.stack(np.meshgrid(x, x, x, indexing="ij"), axis=-1).reshape(-1, 3)

    def idx(i, j, k):
        return (i * (n + 1) + j) * (n + 1) + k

    i, j, k = np.meshgrid(np.arange(n), np.arange(n), np.arange(n), indexing="ij")
    i, j, k = i.ravel(), j.ravel(), k.ravel()
    hexes = np.stack(
        [
            idx(i, j, k),
            idx(i + 1, j, k),
            idx(i + 1, j + 1, k),
            idx(i, j + 1, k),
            idx(i, j, k + 1),
            idx(i + 1, j, k + 1),
            idx(i + 1, j + 1, k + 1),
            idx(i, j + 1, k + 1),
        ],
        axis=1,
    )
    return meshioplusplus.Mesh(pts, [("hexahedron", hexes)])


def _surface_size(n, queue):
    queue.put(len(meshioplusplus.extract_surface(_cube(n)).cells[0].data))


def test_the_pool_works_in_a_forked_child():
    n = 24
    expected = len(meshioplusplus.extract_surface(_cube(n)).cells[0].data)
    assert expected == 6 * n * n
    ctx = multiprocessing.get_context("fork")
    queue = ctx.Queue()
    child = ctx.Process(target=_surface_size, args=(n, queue))
    child.start()
    child.join(timeout=120)
    assert child.exitcode == 0
    assert queue.get(timeout=10) == expected

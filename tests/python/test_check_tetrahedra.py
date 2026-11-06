"""check_tetrahedra: inverted, degenerate and overlapping tetrahedra."""

import sys

import numpy as np
import pytest

import meshioplusplus as mio

try:
    from meshioplusplus import _core  # noqa: F401

    HAVE_CORE = True
except ImportError:  # pragma: no cover
    HAVE_CORE = False

needs_core = pytest.mark.skipif(not HAVE_CORE, reason="needs the compiled core")

CORNER = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)


def _mesh(points, tets, kind="tetra"):
    return mio.Mesh(np.asarray(points, dtype=float), [(kind, np.asarray(tets))])


@pytest.fixture(params=["core", "numpy"])
def engine(request, monkeypatch):
    # Both the compiled operation and its numpy twin must give the same answers.
    if request.param == "numpy":
        monkeypatch.setitem(sys.modules, "meshioplusplus._core", None)
    elif not HAVE_CORE:
        pytest.skip("needs the compiled core")
    return request.param


def test_a_positive_tetrahedron_is_valid(engine):
    report = mio.check_tetrahedra(_mesh(CORNER, [[0, 1, 2, 3]]))
    assert report["num_tetrahedra"] == 1
    assert report["num_inverted"] == report["num_degenerate"] == 0
    assert report["num_overlaps"] == 0 and report["valid"] is True
    assert report["min_volume"] == pytest.approx(1 / 6)
    assert report["max_volume"] == pytest.approx(1 / 6)
    assert report["warnings"] == []


def test_swapping_two_nodes_inverts_it(engine):
    report = mio.check_tetrahedra(_mesh(CORNER, [[0, 2, 1, 3]]))
    assert report["num_inverted"] == 1 and report["valid"] is False
    assert report["min_volume"] == pytest.approx(-1 / 6)


def test_a_flat_tetrahedron_is_degenerate(engine):
    flat = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]]
    report = mio.check_tetrahedra(_mesh(flat, [[0, 1, 2, 3]]))
    assert report["num_degenerate"] == 1 and report["num_inverted"] == 0
    assert report["valid"] is False


def test_the_factor_sets_the_tolerance(engine):
    thin = _mesh([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0.2, 0.2, 1e-4]], [[0, 1, 2, 3]])
    assert mio.check_tetrahedra(thin)["valid"] is True
    assert mio.check_tetrahedra(thin, zero_volume_factor=1e-3)["num_degenerate"] == 1


def test_overlaps_across_a_shared_face(engine):
    apart = _mesh(np.vstack([CORNER, [[0, 0, -1]]]), [[0, 1, 2, 3], [0, 2, 1, 4]])
    assert mio.check_tetrahedra(apart)["num_overlaps"] == 0
    same_side = _mesh(
        np.vstack([CORNER, [[0.2, 0.2, 0.5]]]), [[0, 1, 2, 3], [0, 1, 2, 4]]
    )
    report = mio.check_tetrahedra(same_side)
    assert report["num_overlaps"] == 1 and report["valid"] is False


def test_a_face_used_by_three_tetrahedra_is_not_examined(engine):
    # Three tetrahedra on one face: not a pair, so not counted.
    points = np.vstack([CORNER, [[0.2, 0.2, 0.5], [0.3, 0.3, 0.4]]])
    report = mio.check_tetrahedra(
        _mesh(points, [[0, 1, 2, 3], [0, 1, 2, 4], [0, 1, 2, 5]])
    )
    assert report["num_overlaps"] == 0


def test_quadratic_tetrahedra_are_examined_by_their_corners(engine):
    data = np.arange(10).reshape(1, 10)
    points = np.vstack([CORNER, np.zeros((6, 3)) + 0.25])
    report = mio.check_tetrahedra(_mesh(points, data, kind="tetra10"))
    assert report["num_tetrahedra"] == 1


def test_a_mesh_without_tetrahedra_is_reported(engine):
    quad = _mesh([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]], [[0, 1, 2, 3]], "quad")
    report = mio.check_tetrahedra(quad)
    assert report["num_tetrahedra"] == 0 and report["valid"] is True
    assert report["warnings"] and "no tetrahedra" in report["warnings"][0]


@pytest.mark.parametrize("factor", [-1.0, float("nan"), float("inf")])
def test_a_bad_factor_is_refused(factor):
    with pytest.raises(ValueError, match="zero_volume_factor"):
        mio.check_tetrahedra(_mesh(CORNER, [[0, 1, 2, 3]]), zero_volume_factor=factor)


@needs_core
def test_the_core_and_the_numpy_twin_agree(monkeypatch):
    rng = np.random.default_rng(3)
    points = rng.random((40, 3))
    tets = rng.integers(0, 40, size=(60, 4))
    tets = tets[[len(set(t)) == 4 for t in tets]]
    mesh = _mesh(points, tets)
    native = mio.check_tetrahedra(mesh)
    monkeypatch.setitem(sys.modules, "meshioplusplus._core", None)
    twin = mio.check_tetrahedra(mesh)
    for key in ("num_tetrahedra", "num_inverted", "num_degenerate", "num_overlaps"):
        assert native[key] == twin[key]
    assert native["min_volume"] == pytest.approx(twin["min_volume"])
    assert native["max_volume"] == pytest.approx(twin["max_volume"])


@needs_core
def test_an_octree_mesh_is_valid():
    out = mio.octree_mesh(mio.icosa_sphere(1.0, 3), max_size=0.5, min_size=0.1)
    report = mio.check_tetrahedra(out)
    assert report["num_tetrahedra"] > 0
    assert report["num_inverted"] == 0 and report["num_overlaps"] == 0

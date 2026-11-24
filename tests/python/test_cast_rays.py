"""Rays and point location (roadmap 9.3): the Python surface of ``cast_rays``
and ``locate_points``, their CLI verbs and the agreement with ``interpolate``."""

import json

import numpy as np
import pytest

import meshioplusplus as mio

_core = pytest.importorskip("meshioplusplus._core")
pytestmark = pytest.mark.skipif(
    not hasattr(_core, "cast_rays"), reason="needs the compiled core"
)


def _cube():
    pts = np.array(
        [
            [0, 0, 0],
            [1, 0, 0],
            [1, 1, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 0, 1],
            [1, 1, 1],
            [0, 1, 1],
        ],
        dtype=float,
    )
    tris = [
        [0, 2, 1], [0, 3, 2], [4, 5, 6], [4, 6, 7], [0, 1, 5], [0, 5, 4],
        [1, 2, 6], [1, 6, 5], [2, 3, 7], [2, 7, 6], [3, 0, 4], [3, 4, 7],
    ]  # fmt: skip
    return mio.Mesh(pts, [("triangle", tris)])


def test_rays_hit_the_faces_of_a_cube():
    r = mio.cast_rays(
        _cube(),
        [[-1, 0.5, 0.5], [2, 0.5, 0.5], [0.5, 0.5, 5]],
        [[1, 0, 0], [-1, 0, 0], [0, 0, -1]],
    )
    assert r["hit"].tolist() == [True, True, True]
    assert r["distance"] == pytest.approx([1.0, 1.0, 4.0])
    assert r["normal"] == pytest.approx(np.array([[-1, 0, 0], [1, 0, 0], [0, 0, 1]]))
    assert r["num_hits"] == 3
    assert r["point"][2] == pytest.approx([0.5, 0.5, 1.0])


def test_misses_are_nan_and_the_distance_limit_applies():
    r = mio.cast_rays(_cube(), [[-1, 5, 5]], [1, 0, 0])
    assert not r["hit"][0] and np.isnan(r["distance"][0]) and r["cell"][0] == -1
    near = mio.cast_rays(_cube(), [[-1, 0.5, 0.5]], [1, 0, 0], max_distance=0.5)
    assert near["num_hits"] == 0


def test_cast_rays_refusals():
    with pytest.raises(ValueError):
        mio.cast_rays(_cube(), [[0, 0, 0]], [0, 0, 0])
    with pytest.raises(ValueError):
        mio.cast_rays(_cube(), [[0, 0, 0]], [1, 0, 0], region="nope")
    tet = mio.Mesh(
        [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], [("tetra", [[0, 1, 2, 3]])]
    )
    with pytest.raises(ValueError):
        mio.cast_rays(tet, [[2, 2, 2]], [-1, -1, -1])


def test_locate_points_in_hexahedra():
    box = mio.box(dims=(2, 2, 2), size=(2.0, 2.0, 2.0))
    r = mio.locate_points(box, [[0.5, 0.5, 0.5], [1.5, 1.5, 1.5], [5, 5, 5]])
    assert r["cell"].tolist() == [0, 7, -1]
    assert r["num_located"] == 2
    assert (r["nodes"][2] == -1).all() and (r["barycentric"][2] == 0).all()
    p = np.array([0.5, 0.5, 0.5])
    pts = np.asarray(box.points)
    assert r["barycentric"][0] @ pts[r["nodes"][0]] == pytest.approx(p)
    assert r["barycentric"][0].sum() == pytest.approx(1.0)


def test_locate_points_agrees_with_interpolate():
    box = mio.box(dims=(3, 3, 3), size=(1.0, 1.0, 1.0))
    pts = np.asarray(box.points)
    f = np.sin(3 * pts[:, 0]) + pts[:, 1] ** 2 - 0.5 * pts[:, 2]
    box.point_data["f"] = f
    q = np.random.default_rng(3).uniform(0.05, 0.95, size=(40, 3))
    got = mio.interpolate(
        box,
        mio.Mesh(q, [("vertex", np.arange(40).reshape(-1, 1))]),
        method="barycentric",
    ).point_data["f"]
    r = mio.locate_points(box, q)
    assert (r["cell"] >= 0).all()
    mine = np.array(
        [
            sum(w * f[n] for w, n in zip(wi, ni) if n >= 0)
            for wi, ni in zip(r["barycentric"], r["nodes"])
        ]
    )
    assert np.array_equal(mine, got)


def test_locate_points_refusals():
    box = mio.box(dims=(1, 1, 1))
    with pytest.raises(ValueError):
        mio.locate_points(box, [[0.5, 0.5, 0.5]], tolerance=-1.0)
    line = mio.Mesh([[0, 0, 0], [1, 0, 0]], [("line", [[0, 1]])])
    with pytest.raises(ValueError):
        mio.locate_points(line, [[0.5, 0, 0]])


def test_cli_verbs(tmp_path, capsys):
    from meshioplusplus._cli._main import main

    s = str(tmp_path / "cube.vtu")
    mio.write(s, _cube())
    rays = str(tmp_path / "rays.vtu")
    mio.write(rays, mio.Mesh([[-1, 0.5, 0.5], [-1, 9, 9]], []))
    out = str(tmp_path / "hits.vtu")
    assert main(["cast-rays", s, rays, out, "--direction", "1,0,0", "--json"]) in (
        0,
        None,
    )
    rep = json.loads(capsys.readouterr().out)
    assert rep["num_rays"] == 2 and rep["num_hits"] == 1
    assert mio.read(out).point_data["ray:distance"][0] == pytest.approx(1.0)

    box = str(tmp_path / "box.vtu")
    mio.write(box, mio.box(dims=(2, 2, 2), size=(2.0, 2.0, 2.0)))
    pts = str(tmp_path / "pts.vtu")
    mio.write(pts, mio.Mesh([[0.5, 0.5, 0.5], [9, 9, 9]], []))
    out2 = str(tmp_path / "located.vtu")
    assert main(["locate", box, pts, out2, "--json"]) in (0, None)
    rep = json.loads(capsys.readouterr().out)
    assert rep["num_points"] == 2 and rep["num_located"] == 1

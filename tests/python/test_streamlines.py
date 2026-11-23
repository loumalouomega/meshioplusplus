"""Streamlines of a point vector field as a mesh of line segments (roadmap 15.4)."""

import numpy as np
import pytest

import meshioplusplus as mp
from meshioplusplus import _streamlines

pytest.importorskip("meshioplusplus._core")


def _grid(n=20, field=lambda x, y: (-y, x, 0.0)):
    xs = np.linspace(-1.0, 1.0, n + 1)
    pts = np.array([[x, y, 0.0] for y in xs for x in xs])
    tris = []
    for j in range(n):
        for i in range(n):
            a = j * (n + 1) + i
            tris.append([a, a + 1, a + n + 2])
            tris.append([a, a + n + 2, a + n + 1])
    mesh = mp.Mesh(pts, [("triangle", np.array(tris))])
    mesh.point_data["u"] = np.array([field(x, y) for x, y, _ in pts])
    return mesh


def test_a_rotation_traces_circles_of_constant_radius():
    out, rep = mp.streamlines(
        _grid(),
        "u",
        seed_points=[[0.5, 0.0, 0.0]],
        direction="forward",
        max_length=1.0,
        return_report=True,
    )
    assert rep["num_lines"] == 1 and rep["dim"] == 2
    r = np.linalg.norm(out.points[:, :2], axis=1)
    np.testing.assert_allclose(r, 0.5, atol=2e-3)
    t = out.point_data["IntegrationTime"]
    assert t[-1] > 5.0 and np.all(np.diff(t) > 0)
    np.testing.assert_allclose(out.point_data["Speed"], r, atol=2e-3)
    assert _streamlines.STOP_REASONS[rep["end_forward"][0]] == "max_length"


def test_both_directions_have_negative_and_positive_time():
    out = mp.streamlines(_grid(), "u", seed_points=[[0.5, 0, 0]], max_length=0.3)
    t = out.point_data["IntegrationTime"]
    assert t[0] < 0.0 < t[-1]
    assert len(out.cells) == 1 and out.cells[0].type == "line"
    assert len(out.cells[0].data) == len(out.points) - 1
    assert list(out.cell_data["SeedId"][0]) == [0] * len(out.cells[0].data)


def test_the_integrators_agree():
    kw = dict(seed_points=[[0.4, 0.1, 0.0]], direction="forward", max_length=0.4)
    a = mp.streamlines(_grid(), "u", integrator="rk4", **kw)
    b = mp.streamlines(_grid(), "u", integrator="rk45", **kw)
    ra, rb = np.linalg.norm(a.points[-1, :2]), np.linalg.norm(b.points[-1, :2])
    assert ra == pytest.approx(rb, abs=2e-3)
    assert rb == pytest.approx(np.hypot(0.4, 0.1), abs=2e-3)


def test_stop_reasons():
    uniform = _grid(8, lambda x, y: (1.0, 0.0, 0.0))
    _, rep = mp.streamlines(
        uniform, "u", seed_points=[[0, 0, 0]], max_length=5.0, return_report=True
    )
    names = _streamlines.STOP_REASONS
    assert (
        names[rep["end_backward"][0]] == names[rep["end_forward"][0]] == "left_domain"
    )
    out, rep = mp.streamlines(
        uniform,
        "u",
        seed_points=[[0, 0, 0]],
        max_length=5.0,
        max_steps=3,
        return_report=True,
    )
    assert names[rep["end_forward"][0]] == "max_steps" and len(out.points) == 7
    _, rep = mp.streamlines(
        uniform, "u", seed_points=[[0, 0, 0]], max_length=0.1, return_report=True
    )
    assert names[rep["end_forward"][0]] == "max_length"
    slow = _grid(8, lambda x, y: (0.5 * (1.0 - x), 0.0, 0.0))
    _, rep = mp.streamlines(
        slow,
        "u",
        seed_points=[[-0.5, 0, 0]],
        direction="forward",
        max_length=5.0,
        terminal_speed=0.05,
        return_report=True,
    )
    assert names[rep["end_forward"][0]] == "stagnation"


def test_a_seed_outside_the_mesh_makes_no_line():
    _, rep = mp.streamlines(
        _grid(),
        "u",
        seed_points=[[0.5, 0, 0], [5, 5, 0]],
        max_length=0.1,
        return_report=True,
    )
    assert rep["num_lines"] == 1 and rep["num_seeds_outside"] == 1
    assert list(rep["seed_ids"]) == [0]


def test_seed_counts_and_two_dimensional_seed_points():
    out, rep = mp.streamlines(
        _grid(), "u", seeds=12, max_length=0.2, return_report=True
    )
    assert 6 < rep["num_lines"] <= 12
    assert out.point_data["u"].shape[1] == 3
    out2 = mp.streamlines(_grid(), "u", seed_points=[[0.5, 0.0]], max_length=0.1)
    assert len(out2.points) > 2


def test_a_volume_is_traced_in_its_volume():
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
        float,
    )
    tets = [
        [0, 1, 2, 6],
        [0, 2, 3, 6],
        [0, 3, 7, 6],
        [0, 7, 4, 6],
        [0, 4, 5, 6],
        [0, 5, 1, 6],
    ]
    mesh = mp.Mesh(pts, [("tetra", np.array(tets))])
    mesh.point_data["u"] = np.tile([0.0, 0.0, 2.0], (8, 1))
    out, rep = mp.streamlines(
        mesh, "u", seed_points=[[0.3, 0.4, 0.5]], max_length=2.0, return_report=True
    )
    assert rep["dim"] == 3
    np.testing.assert_allclose(out.points[:, 0], 0.3, atol=1e-9)
    np.testing.assert_allclose(out.points[:, 1], 0.4, atol=1e-9)
    t = out.point_data["IntegrationTime"]
    assert t[-1] - t[0] == pytest.approx((out.points[-1, 2] - out.points[0, 2]) / 2.0)


def test_a_two_component_field_is_padded():
    mesh = _grid()
    mesh.point_data["u2"] = mesh.point_data["u"][:, :2].copy()
    out = mp.streamlines(mesh, "u2", seed_points=[[0.5, 0, 0]], max_length=0.2)
    assert out.point_data["u2"].shape[1] == 2


@pytest.mark.parametrize(
    "kwargs,needle",
    [
        (dict(array="nothing"), "nothing"),
        (dict(array="u", direction="up"), "direction"),
        (dict(array="u", integrator="euler"), "integrator"),
        (dict(array="u", max_length=0.0), "positive"),
        (dict(array="u", seed_points=[[1, 2, 3, 4]]), "shape"),
    ],
)
def test_refusals(kwargs, needle):
    array = kwargs.pop("array")
    with pytest.raises(ValueError, match=needle):
        mp.streamlines(_grid(), array, **kwargs)


def test_a_scalar_field_and_a_mesh_without_cells_are_refused():
    mesh = _grid()
    mesh.point_data["t"] = np.ones(len(mesh.points))
    with pytest.raises(ValueError, match="2 or 3"):
        mp.streamlines(mesh, "t")
    lines = mp.Mesh(np.array([[0, 0, 0], [1, 0, 0.0]]), [("line", np.array([[0, 1]]))])
    lines.point_data["u"] = np.array([[1.0, 0, 0], [1.0, 0, 0]])
    with pytest.raises(ValueError, match="cells"):
        mp.streamlines(lines, "u")


def test_the_result_is_reproducible():
    a = mp.streamlines(_grid(), "u", seeds=30, integrator="rk45")
    b = mp.streamlines(_grid(), "u", seeds=30, integrator="rk45")
    np.testing.assert_array_equal(a.points, b.points)
    np.testing.assert_array_equal(
        a.point_data["IntegrationTime"], b.point_data["IntegrationTime"]
    )


def test_the_cli_verb(tmp_path, capsys):
    from meshioplusplus._cli._main import main

    src = str(tmp_path / "flow.vtu")
    mp.write(src, _grid())
    out = str(tmp_path / "lines.vtu")
    rc = main(
        [
            "streamlines",
            src,
            out,
            "--array",
            "u",
            "--seeds",
            "8",
            "--max-length",
            "0.2",
            "--json",
        ]
    )
    assert rc in (0, None)
    import json

    report = json.loads(capsys.readouterr().out)
    assert report["num_lines"] > 0 and report["dim"] == 2
    back = mp.read(out)
    assert "IntegrationTime" in back.point_data

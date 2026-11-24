"""Mass, centre of mass and inertia (roadmap 9.4): closed forms on a box, the
parallel-axis theorem, shells, density arrays, regions and the refusals."""

import numpy as np
import pytest

import meshioplusplus as mio

_core = pytest.importorskip("meshioplusplus._core")
pytestmark = pytest.mark.skipif(
    not hasattr(_core, "mass_properties"), reason="needs the compiled core"
)


def _box_points(x=0.0, y=0.0, z=0.0):
    dx = [0, 2, 2, 0, 0, 2, 2, 0]
    dy = [0, 0, 3, 3, 0, 0, 3, 3]
    dz = [0, 0, 0, 0, 4, 4, 4, 4]
    return np.array([[x + a, y + b, z + c] for a, b, c in zip(dx, dy, dz)])


def _hex_box(x=0.0, y=0.0, z=0.0):
    return mio.Mesh(_box_points(x, y, z), [("hexahedron", [[0, 1, 2, 3, 4, 5, 6, 7]])])


def _tri_box(flip=False):
    quads = [
        [0, 3, 2, 1],
        [4, 5, 6, 7],
        [0, 1, 5, 4],
        [3, 7, 6, 2],
        [0, 4, 7, 3],
        [1, 2, 6, 5],
    ]
    tris = []
    for q in quads:
        tris += [[q[0], q[1], q[2]], [q[0], q[2], q[3]]]
    if flip:
        tris = [[t[0], t[2], t[1]] for t in tris]
    return mio.Mesh(_box_points(), [("triangle", tris)])


def _check_box(report, rho, origin=(0.0, 0.0, 0.0)):
    t = report["total"]
    m = 24.0 * rho
    assert t["measure"] == pytest.approx(24.0)
    assert t["mass"] == pytest.approx(m)
    assert t["centre"] == pytest.approx(
        [origin[0] + 1.0, origin[1] + 1.5, origin[2] + 2.0], rel=1e-9
    )
    assert t["inertia"][:3] == pytest.approx(
        [m * 25 / 12, m * 20 / 12, m * 13 / 12], rel=1e-9
    )
    assert t["inertia"][3:] == pytest.approx([0, 0, 0], abs=1e-8 * m)
    assert t["principal_moments"][0] == pytest.approx(m * 13 / 12, rel=1e-9)
    assert t["principal_moments"][2] == pytest.approx(m * 25 / 12, rel=1e-9)
    assert t["gyration_radii"][0] == pytest.approx((13 / 12) ** 0.5, rel=1e-9)


def test_hexahedron_box():
    r = mio.mass_properties(_hex_box(), density=2.0)
    assert r["kind"] == "solid"
    _check_box(r, 2.0)


def test_far_from_origin_keeps_accuracy():
    r = mio.mass_properties(_hex_box(1e6, -2e6, 3e6), density=7850.0)
    _check_box(r, 7850.0, (1e6, -2e6, 3e6))


@pytest.mark.parametrize("flip", [False, True])
def test_closed_surface_equals_the_solid(flip):
    _check_box(mio.mass_properties(_tri_box(flip), density=2.0), 2.0)


def test_parallel_axis_about_the_origin():
    mesh = _hex_box(1.0, 2.0, 3.0)
    c = mio.mass_properties(mesh, density=1.5)["total"]
    g = mio.mass_properties(mesh, density=1.5, about="origin")["total"]
    m, d = c["mass"], np.array(c["centre"])
    d2 = d @ d
    ic = np.array(c["inertia"])
    expect = ic + m * np.array(
        [
            d2 - d[0] ** 2,
            d2 - d[1] ** 2,
            d2 - d[2] ** 2,
            -d[0] * d[1],
            -d[1] * d[2],
            -d[2] * d[0],
        ]
    )
    assert np.array(g["inertia"]) == pytest.approx(expect, rel=1e-9, abs=1e-9)


def test_principal_axes_diagonalise_the_tensor():
    t = mio.mass_properties(_hex_box(), density=1.0)["total"]
    xx, yy, zz, xy, yz, zx = t["inertia"]
    inertia = np.array([[xx, xy, zx], [xy, yy, yz], [zx, yz, zz]])
    axes = np.array(t["principal_axes"]).reshape(3, 3)
    diag = axes.T @ inertia @ axes
    assert np.diag(diag) == pytest.approx(t["principal_moments"], rel=1e-9)


def test_tetrahedron():
    mesh = mio.Mesh(
        [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], [("tetra", [[0, 1, 2, 3]])]
    )
    t = mio.mass_properties(mesh, density=6.0, about="origin")["total"]
    assert t["measure"] == pytest.approx(1 / 6)
    assert t["mass"] == pytest.approx(1.0)
    assert t["centre"] == pytest.approx([0.25, 0.25, 0.25])
    assert t["inertia"][0] == pytest.approx(0.2)
    assert t["inertia"][3] == pytest.approx(-0.05)


def test_shell_of_a_plate():
    mesh = mio.Mesh(
        [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]], [("quad", [[0, 1, 2, 3]])]
    )
    r = mio.mass_properties(mesh, density=8.0, thickness=0.5)
    assert r["kind"] == "shell"
    t = r["total"]
    assert t["measure"] == pytest.approx(1.0)
    assert t["mass"] == pytest.approx(4.0)
    assert t["inertia"][2] == pytest.approx(4.0 / 6.0)
    assert t["inertia"][0] == pytest.approx(4.0 / 12.0)


def test_open_surface_without_thickness_is_refused():
    mesh = mio.Mesh(
        [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]], [("quad", [[0, 1, 2, 3]])]
    )
    with pytest.raises(ValueError, match="closed"):
        mio.mass_properties(mesh)


def test_density_array_and_regions():
    pts = [[x, y, z] for z in (0, 1, 2) for y in (0, 1) for x in (0, 1)]
    cells = [[0, 1, 3, 2, 4, 5, 7, 6], [4, 5, 7, 6, 8, 9, 11, 10]]
    mesh = mio.Mesh(
        pts,
        [("hexahedron", cells)],
        cell_data={"rho": [np.array([1.0, 3.0])]},
    )
    mesh.regions.append(mio.Region("top", "cell", [1]))
    r = mio.mass_properties(mesh, density_array="rho")
    assert r["total"]["mass"] == pytest.approx(4.0)
    assert r["total"]["centre"][2] == pytest.approx((0.5 + 3 * 1.5) / 4)
    assert [e["name"] for e in r["regions"]] == ["top"]
    assert r["regions"][0]["mass"] == pytest.approx(3.0)
    only = mio.mass_properties(mesh, region="top")
    assert only["total"]["mass"] == pytest.approx(1.0)
    assert only["regions"] == []


def test_refusals():
    box = _hex_box()
    with pytest.raises(ValueError):
        mio.mass_properties(box, density_array="nope")
    with pytest.raises(ValueError):
        mio.mass_properties(box, region="nope")
    with pytest.raises(ValueError):
        mio.mass_properties(box, thickness=-1.0)
    with pytest.raises(ValueError):
        mio.mass_properties(box, density=float("nan"))
    with pytest.raises(ValueError):
        mio.mass_properties(box, about="moon")
    line = mio.Mesh([[0, 0, 0], [1, 0, 0]], [("line", [[0, 1]])])
    with pytest.raises(ValueError):
        mio.mass_properties(line)


def test_cli_json(tmp_path, capsys):
    from meshioplusplus._cli._main import main

    path = str(tmp_path / "box.vtu")
    mio.write(path, _hex_box())
    assert main(["mass", path, "--density", "2", "--json"]) in (0, None)
    import json

    out = json.loads(capsys.readouterr().out)
    assert out["total"]["mass"] == pytest.approx(48.0)

"""Octree / quadtree leaves around a skin (``octree``).

The C++/numpy parity tests assert **raw bytes**, as in test_voxelize.py.
"""

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import octree, voxelize


def _sphere(radius=1.0, nu=24, nv=12):
    pts = [[0, 0, -radius]]
    for j in range(1, nv):
        th = np.pi * j / nv
        for i in range(nu):
            ph = 2 * np.pi * i / nu
            pts.append(
                [
                    radius * np.sin(th) * np.cos(ph),
                    radius * np.sin(th) * np.sin(ph),
                    -radius * np.cos(th),
                ]
            )
    pts.append([0, 0, radius])
    north = len(pts) - 1

    def at(j, i):
        return 1 + (j - 1) * nu + i % nu

    tris = []
    for i in range(nu):
        tris.append([0, at(1, i + 1), at(1, i)])
        tris.append([north, at(nv - 1, i), at(nv - 1, i + 1)])
    for j in range(1, nv - 1):
        for i in range(nu):
            tris.append([at(j, i), at(j, i + 1), at(j + 1, i + 1)])
            tris.append([at(j, i), at(j + 1, i + 1), at(j + 1, i)])
    return meshioplusplus.Mesh(
        np.array(pts, dtype=float), [("triangle", np.array(tris, dtype=np.int64))]
    )


def _circle(radius=1.0, n=64, close=True):
    t = 2 * np.pi * np.arange(n) / n
    pts = np.column_stack([radius * np.cos(t), radius * np.sin(t), np.zeros(n)])
    lines = [[i, i + 1] for i in range(n - 1)]
    if close:
        lines.append([n - 1, 0])
    return meshioplusplus.Mesh(pts, [("line", np.array(lines, dtype=np.int64))])


BOX = [-1.2, -1.2, -1.2, 1.2, 1.2, 1.2]


def _volumes(mesh, report, dim):
    level = np.asarray(mesh.cell_data["octree:level"][0])
    edge = report["root_size"] / 2.0**level
    return float(np.sum(edge**dim))


def test_uniform_tree_is_the_voxelize_lattice():
    sphere = _sphere(nu=16, nv=8)
    b = [-1, -1, -1, 1, 1, 1]
    tree, rep = octree(sphere, min_level=3, max_level=3, bounds=b, return_report=True)
    vox = voxelize(sphere, resolution=(8, 8, 8), bounds=b)
    assert rep["num_leaves"] == 512 and rep["num_hanging"] == 0
    assert np.asarray(tree.points).tobytes() == np.asarray(vox.points).tobytes()


def test_cut_leaves_sit_at_max_level():
    tree = octree(_sphere(), bounds=BOX, max_level=4)
    level = np.asarray(tree.cell_data["octree:level"][0])
    cut = np.asarray(tree.cell_data["octree:cut"][0])
    assert np.all(level[cut == 1] == 4)
    assert level.min() < 4


def test_inside_volume_converges_to_the_ball():
    errs = []
    for max_level in (3, 4, 5):
        tree, rep = octree(
            _sphere(nu=48, nv=24),
            bounds=BOX,
            max_level=max_level,
            fill="inside",
            return_report=True,
        )
        errs.append(abs(_volumes(tree, rep, 3) - 4 / 3 * np.pi) / (4 / 3 * np.pi))
    assert errs[0] > errs[1] > errs[2]
    assert errs[2] < 0.1


def test_quadtree_of_a_circle_converges_to_its_area():
    tree, rep = octree(
        _circle(n=256),
        bounds=[-1.2, -1.2, 0, 1.2, 1.2, 0],
        max_level=6,
        fill="inside",
        return_report=True,
    )
    assert rep["dimension"] == 2
    assert tree.cells[0].type == "quad"
    assert abs(_volumes(tree, rep, 2) - np.pi) < 0.05
    assert rep["num_hanging"] > 0


def test_face_balance_holds():
    tree = octree(_sphere(), bounds=BOX, max_level=4, balance="face")
    pts = np.asarray(tree.points)
    conn = np.asarray(tree.cells[0].data)
    lo, hi = pts[conn[:, 0]], pts[conn[:, 6]]
    level = np.asarray(tree.cell_data["octree:level"][0])
    eps = 1e-12
    for a in range(len(conn)):
        touch = (np.abs(hi[a] - lo) < eps) | (np.abs(hi - lo[a]) < eps)[:, :]
        overlap = (lo[a] < hi - eps) & (lo < hi[a] - eps)
        face = (touch.sum(axis=1) == 1) & (overlap.sum(axis=1) == 2)
        assert np.all(np.abs(level[face] - level[a]) <= 1)


def test_hanging_flag_is_set_where_levels_meet():
    tree, rep = octree(_sphere(), bounds=BOX, max_level=3, return_report=True)
    h = np.asarray(tree.point_data["octree:hanging"])
    assert h.sum() == rep["num_hanging"] > 0


def test_open_polyline_is_refused_for_inside():
    with pytest.raises(ValueError, match="does not close"):
        octree(_circle(close=False), max_level=3, fill="inside")


def test_budget_is_enforced():
    with pytest.raises(ValueError, match="max_cells"):
        octree(_sphere(), bounds=BOX, max_level=5, max_cells=100)


def test_bad_arguments():
    with pytest.raises(ValueError):
        octree(_sphere(), min_level=3, max_level=2)
    with pytest.raises(ValueError):
        octree(_sphere(), balance="bogus")
    with pytest.raises(ValueError):
        octree(_sphere(), fill="bogus")


@pytest.mark.parametrize(
    "case",
    [
        dict(skin="sphere", bounds=BOX, max_level=3),
        dict(skin="sphere", bounds=BOX, max_level=3, balance="full"),
        dict(skin="sphere", bounds=BOX, max_level=3, fill="surface"),
        dict(
            skin="sphere", bounds=BOX, max_level=3, fill="inside", attach_distance=True
        ),
        dict(
            skin="sphere", bounds=BOX, max_level=3, sample_distance=0.2, balance="none"
        ),
        dict(
            skin="circle",
            bounds=[-1.2, -1.2, 0, 1.2, 1.2, 0],
            max_level=4,
            fill="inside",
        ),
        dict(skin="circle", max_level=4, root_cells=(2, 2, 1), padding=0.1),
    ],
)
def test_cpp_matches_python(case):
    core = pytest.importorskip("meshioplusplus._core")
    from meshioplusplus._octree import _octree_py

    case = dict(case)
    skin = _sphere() if case.pop("skin") == "sphere" else _circle()
    kw = dict(
        root_cells=list(case.get("root_cells", (1, 1, 1))),
        min_level=0,
        max_level=case["max_level"],
        bounds=case.get("bounds"),
        padding=case.get("padding", 0.0),
        padding_relative=0.0,
        dimension=0,
        sample_distance=case.get("sample_distance", 0.0),
        balance=case.get("balance", "face"),
        fill=case.get("fill", "all"),
        attach_distance=case.get("attach_distance", False),
        even_odd=False,
        max_cells=20000000,
    )
    got = core.octree(skin, watertight_check="off", sign="pseudonormal", **kw)
    want, rep = _octree_py(skin, sign="pseudonormal", **kw)
    assert got["num_leaves"] == rep["num_leaves"]
    assert got["num_hanging"] == rep["num_hanging"]
    a, b = got["mesh"], want
    assert np.asarray(a.points).tobytes() == np.asarray(b.points).tobytes()
    assert (
        np.asarray(a.cells[0].data).tobytes() == np.asarray(b.cells[0].data).tobytes()
    )
    for name in b.cell_data:
        if name == "sdf:distance":
            # A distance is a float reduction the two sides evaluate in different
            # orders; the structure it decides (fill, leaves) is pinned by bytes.
            np.testing.assert_allclose(
                a.cell_data[name][0], b.cell_data[name][0], rtol=1e-9, atol=1e-12
            )
            continue
        assert (
            np.asarray(a.cell_data[name][0]).tobytes()
            == np.asarray(b.cell_data[name][0]).tobytes()
        )
    assert (
        np.asarray(a.point_data["octree:hanging"]).tobytes()
        == np.asarray(b.point_data["octree:hanging"]).tobytes()
    )

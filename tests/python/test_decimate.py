"""Surface decimation (``meshioplusplus.decimate``)."""

import numpy as np
import pytest

import meshioplusplus as mp
from meshioplusplus._decimate import _decimate_py


def _grid(n=4):
    """A regular (n+1)x(n+1) triangulated grid over the unit square (z = 0)."""
    pts = np.array(
        [[i / n, j / n, 0.0] for j in range(n + 1) for i in range(n + 1)], dtype=float
    )

    def v(i, j):
        return j * (n + 1) + i

    cells = []
    for j in range(n):
        for i in range(n):
            cells.append([v(i, j), v(i + 1, j), v(i + 1, j + 1)])
            cells.append([v(i, j), v(i + 1, j + 1), v(i, j + 1)])
    return mp.Mesh(pts, [("triangle", np.array(cells, dtype=np.int64))])


def _sphere(levels=3):
    """An octahedron refined ``levels`` times, pushed onto the unit sphere.

    Same recipe as tests/cpp/test_decimate.cpp's ``sphere_mesh``.
    """
    octa = mp.Mesh(
        np.array(
            [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]],
            dtype=float,
        ),
        [
            (
                "triangle",
                np.array(
                    [
                        [0, 2, 4],
                        [2, 1, 4],
                        [1, 3, 4],
                        [3, 0, 4],
                        [2, 0, 5],
                        [1, 2, 5],
                        [3, 1, 5],
                        [0, 3, 5],
                    ],
                    dtype=np.int64,
                ),
            )
        ],
    )
    fine = mp.refine(octa, levels=levels)
    pts = np.asarray(fine.points, dtype=float)
    pts = pts / np.linalg.norm(pts, axis=1)[:, None]
    return mp.Mesh(pts, [("triangle", np.asarray(fine.cells[0].data, dtype=np.int64))])


def _total_faces(mesh):
    return sum(len(cb.data) for cb in mesh.cells)


def _has_point(mesh, xyz):
    return bool((np.asarray(mesh.points) == np.asarray(xyz)).all(axis=1).any())


def _edge_use_counts(mesh):
    counts = {}
    for cb in mesh.cells:
        for row in np.asarray(cb.data):
            for k in range(3):
                a, b = int(row[k]), int(row[(k + 1) % 3])
                key = (a, b) if a < b else (b, a)
                counts[key] = counts.get(key, 0) + 1
    return counts


def test_planar_grid_keeps_outline():
    grid = _grid(4)  # 32 faces, 16 boundary + 9 interior points
    out, report = mp.decimate(grid, target_faces=1, return_report=True)
    assert _total_faces(out) <= 16
    assert report["faces_removed"] == 32 - _total_faces(out)
    for k in range(5):
        t = k / 4
        for p in [(t, 0, 0), (t, 1, 0), (0, t, 0), (1, t, 0)]:
            assert _has_point(out, p)
    # The sheet's area is unchanged: the outline is intact and everything
    # stays in the z = 0 plane.
    assert np.allclose(np.asarray(out.points)[:, 2], 0.0)
    assert mp.compute_stats(out)["total_area"] == pytest.approx(1.0, abs=1e-9)
    assert max(_edge_use_counts(out).values()) <= 2


def test_target_lands_within_one_collapse():
    grid = _grid(4)
    for target in (30, 24, 21, 20, 16):
        out = mp.decimate(grid, target_faces=target)
        assert target - 1 <= _total_faces(out) <= target


def test_ratio_matches_face_target():
    grid = _grid(4)
    by_ratio = mp.decimate(grid, ratio=0.625)  # 20 of 32
    by_faces = mp.decimate(grid, target_faces=20)
    assert _total_faces(by_ratio) == _total_faces(by_faces)


def test_sphere_keeps_shape_and_orientation():
    sphere = _sphere(3)  # 512 faces
    area_in = mp.compute_stats(sphere)["total_area"]
    out = mp.decimate(sphere, ratio=0.25)
    assert 127 <= _total_faces(out) <= 128

    stats = mp.compute_stats(out)
    assert stats["total_area"] == pytest.approx(area_in, rel=0.25)
    assert np.allclose(stats["bbox_min"], [-1, -1, -1], atol=0.25)
    assert np.allclose(stats["bbox_max"], [1, 1, 1], atol=0.25)

    # No flipped normals: the decimated sphere stays star-shaped around the
    # origin.
    pts = np.asarray(out.points)
    tri = pts[np.asarray(out.cells[0].data)]
    normals = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    centroids = tri.mean(axis=1)
    assert (np.einsum("ij,ij->i", normals, centroids) > 0).all()
    assert max(_edge_use_counts(out).values()) <= 2


def test_max_error_stops_correctly():
    # Flat grid: every collapse is error-free, so a tiny budget still decimates
    # down to the pinned outline.
    out, report = mp.decimate(_grid(4), max_error=1e-9, return_report=True)
    assert _total_faces(out) <= 16
    assert report["max_error_applied"] <= 1e-9

    # Strictly convex bump: every collapse has real error, so nothing moves.
    n = 4
    pts = np.array(
        [
            [i / n, j / n, (i / n) ** 2 + (j / n) ** 2]
            for j in range(n + 1)
            for i in range(n + 1)
        ]
    )
    bump = mp.Mesh(pts, [("triangle", np.asarray(_grid(4).cells[0].data))])
    out, report = mp.decimate(
        bump, max_error=1e-12, preserve_features=False, return_report=True
    )
    assert report["faces_removed"] == 0
    assert _total_faces(out) == 32


def test_cube_corners_and_creases_survive():
    cube = mp.extract_skin(
        mp.Mesh(
            np.array(
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
            ),
            [("hexahedron", np.array([[0, 1, 2, 3, 4, 5, 6, 7]]))],
        ),
        linearize=True,
    )
    fine = mp.refine(cube, levels=2)
    out = mp.decimate(fine, ratio=0.25)
    assert _total_faces(out) < 2 * _total_faces(fine)

    # Every input vertex on a cube edge (two coordinates in {0, 1}) is a
    # feature vertex and must survive at its exact position.
    for p in np.asarray(fine.points):
        if sum(1 for c in p if c in (0.0, 1.0)) >= 2:
            assert _has_point(out, p)


def test_frozen_vertex_survives():
    grid = _grid(4)
    centre = 2 * 5 + 2  # (0.5, 0.5)
    out = mp.decimate(grid, target_faces=1, frozen=[centre])
    assert _has_point(out, (0.5, 0.5, 0.0))

    grid.point_sets = {"keep": np.array([centre])}
    out = mp.decimate(grid, target_faces=1, frozen="keep")
    assert _has_point(out, (0.5, 0.5, 0.0))
    assert "keep" in out.point_sets and len(out.point_sets["keep"]) == 1

    with pytest.raises(ValueError, match="no point_set named"):
        mp.decimate(grid, target_faces=1, frozen="nope")
    with pytest.raises(ValueError, match="out of range"):
        mp.decimate(grid, target_faces=1, frozen=[99])


def test_maps_and_sets():
    grid = _grid(4)
    grid.point_sets = {"corners": np.array([0, 4, 20, 24])}
    grid.cell_sets = {"first": [np.array([0, 1])]}
    out, report = mp.decimate(grid, target_faces=16, return_report=True)

    # Pinned corners survive, so their (remapped, deduped) set has 4 entries.
    assert len(out.point_sets["corners"]) == 4
    assert (out.point_sets["corners"] < len(out.points)).all()
    assert all(
        idx is None or (np.asarray(idx) < len(out.cells[b].data)).all()
        for b, idx in enumerate(out.cell_sets["first"])
    )
    assert report["points_removed"] == len(grid.points) - len(out.points)


def test_quads_are_triangulated_blockwise():
    # All seven points of the strip are boundary, so nothing collapses; the
    # quad block still comes back triangulated and blocks stay 1:1.
    strip = mp.Mesh(
        np.array(
            [
                [0, 0, 0],
                [1, 0, 0],
                [2, 0, 0],
                [3, 1, 0],
                [2, 1, 0],
                [1, 1, 0],
                [0, 1, 0],
            ],
            dtype=float,
        ),
        [
            ("triangle", np.array([[0, 1, 5], [0, 5, 6]])),
            ("quad", np.array([[1, 2, 4, 5]])),
            ("triangle", np.array([[2, 3, 4]])),
        ],
    )
    out = mp.decimate(strip, target_faces=1)
    assert [cb.type for cb in out.cells] == ["triangle", "triangle", "triangle"]
    assert [len(cb.data) for cb in out.cells] == [2, 2, 1]
    assert len(out.points) == 7


def test_point_data_blend_and_integer_survivor():
    sphere = _sphere(2)
    n = len(sphere.points)
    sphere.point_data["T"] = np.asarray(sphere.points)[:, 0].astype(np.float32)
    sphere.point_data["id"] = (np.arange(n) % 7).astype(np.int32)
    sphere.cell_data["tag"] = [
        (np.arange(len(sphere.cells[0].data)) % 3).astype(np.int32)
    ]
    out = mp.decimate(sphere, ratio=0.5, preserve_features=False)

    t = np.asarray(out.point_data["T"])
    assert t.dtype == np.float32
    assert (t >= -1.0 - 1e-6).all() and (t <= 1.0 + 1e-6).all()
    ids = np.asarray(out.point_data["id"])
    assert ids.dtype.kind in "iu"
    assert set(np.unique(ids)) <= set(range(7))
    tags = np.asarray(out.cell_data["tag"][0])
    assert len(tags) == len(out.cells[0].data)
    assert "convert:parent_cell" not in out.cell_data


def test_raises_by_name():
    tet = mp.Mesh(
        np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float),
        [("tetra", np.array([[0, 1, 2, 3]]))],
    )
    with pytest.raises(ValueError, match="extract_surface"):
        mp.decimate(tet, ratio=0.5)

    tri6 = mp.Mesh(
        np.array(
            [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0.5, 0, 0], [0.5, 0.5, 0], [0, 0.5, 0]],
            dtype=float,
        ),
        [("triangle6", np.array([[0, 1, 2, 3, 4, 5]]))],
    )
    with pytest.raises(ValueError, match="linearize"):
        mp.decimate(tri6, ratio=0.5)

    mixed = mp.Mesh(
        np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]], dtype=float),
        [("triangle", np.array([[0, 1, 2]])), ("line", np.array([[0, 1]]))],
    )
    with pytest.raises(ValueError, match="line"):
        mp.decimate(mixed, ratio=0.5)

    lines = mp.Mesh(
        np.array([[0, 0, 0], [1, 0, 0]], dtype=float),
        [("line", np.array([[0, 1]]))],
    )
    with pytest.raises(ValueError, match="surface"):
        mp.decimate(lines, ratio=0.5)


@pytest.mark.parametrize(
    "kwargs",
    [
        {},  # no criterion
        {"ratio": 0.5, "target_faces": 4},  # two criteria
        {"ratio": 1.5},  # out of range
        {"ratio": 0.5, "placement": "nearest"},  # unknown placement
    ],
)
def test_invalid_options_raise(kwargs):
    with pytest.raises(ValueError):
        mp.decimate(_grid(2), **kwargs)


def test_determinism_two_runs():
    sphere = _sphere(2)
    a = mp.decimate(sphere, ratio=0.4, preserve_features=False)
    b = mp.decimate(sphere, ratio=0.4, preserve_features=False)
    assert a.points.tobytes() == b.points.tobytes()
    for x, y in zip(a.cells, b.cells):
        assert np.asarray(x.data).tobytes() == np.asarray(y.data).tobytes()


@pytest.mark.parametrize("placement", ["optimal", "midpoint", "endpoint"])
@pytest.mark.parametrize(
    "criterion",
    [{"target_ratio": 0.4}, {"target_faces": 40}, {"max_error": 1e-3}],
)
def test_cpp_matches_python(placement, criterion):
    core = pytest.importorskip("meshioplusplus._core")
    for mesh in (_grid(4), _sphere(2)):
        n = len(mesh.points)
        mesh.point_data["T"] = np.asarray(mesh.points)[:, 0].astype(np.float32)
        mesh.point_data["v"] = np.asarray(mesh.points).astype(np.float64)
        mesh.point_data["id"] = (np.arange(n) % 5).astype(np.int32)
        mesh.cell_data["tag"] = [
            (np.arange(len(mesh.cells[0].data)) % 3).astype(np.int32)
        ]

        args = dict(target_ratio=-1.0, target_faces=-1, max_error=-1.0)
        args.update(criterion)
        got = core.decimate(
            mesh,
            args["target_ratio"],
            args["target_faces"],
            args["max_error"],
            placement,
            True,
            False,  # preserve_features off: the coarse sphere pins itself
            30.0,
            None,
        )
        ref_mesh, ref_pm, ref_cm, ref_report = _decimate_py(
            mesh,
            args["target_ratio"],
            args["target_faces"],
            args["max_error"],
            placement,
            True,
            False,
            30.0,
            None,
        )

        out = got["mesh"]
        assert out.points.dtype == ref_mesh.points.dtype
        assert np.array_equal(out.points, ref_mesh.points)
        assert [cb.type for cb in out.cells] == [cb.type for cb in ref_mesh.cells]
        for x, y in zip(out.cells, ref_mesh.cells):
            dx, dy = np.asarray(x.data), np.asarray(y.data)
            assert dx.dtype == dy.dtype and np.array_equal(dx, dy)
        assert set(out.point_data) == set(ref_mesh.point_data)
        for k in out.point_data:
            x, y = np.asarray(out.point_data[k]), np.asarray(ref_mesh.point_data[k])
            assert x.dtype == y.dtype and np.array_equal(x, y)
        assert set(out.cell_data) == set(ref_mesh.cell_data)
        for k in out.cell_data:
            for xb, yb in zip(out.cell_data[k], ref_mesh.cell_data[k]):
                xb, yb = np.asarray(xb), np.asarray(yb)
                assert xb.dtype == yb.dtype and np.array_equal(xb, yb)
        assert np.array_equal(np.asarray(got["point_map"]), ref_pm)
        assert len(got["cell_maps"]) == len(ref_cm)
        for xm, ym in zip(got["cell_maps"], ref_cm):
            assert np.array_equal(np.asarray(xm), ym)
        assert got["faces_removed"] == ref_report["faces_removed"]
        assert got["points_removed"] == ref_report["points_removed"]
        assert got["collapses_rejected"] == ref_report["collapses_rejected"]
        assert got["max_error_applied"] == ref_report["max_error_applied"]


def test_report_shape():
    out, report = mp.decimate(_grid(2), ratio=0.5, return_report=True)
    assert set(report) == {
        "faces_removed",
        "points_removed",
        "collapses_rejected",
        "max_error_applied",
    }
    assert isinstance(out, mp.Mesh)

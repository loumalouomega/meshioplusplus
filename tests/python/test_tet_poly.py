"""Tests for tetrahedral <-> polyhedral conversion (``tet_to_poly`` and
``poly_to_tet``)."""

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import poly_to_tet, tet_to_poly

try:
    from meshioplusplus import _core
except ImportError:  # pragma: no cover - a pure-Python build
    _core = None

needs_core = pytest.mark.skipif(
    _core is None, reason="tet_to_poly/poly_to_tet have no pure-Python fallback"
)


def _box(n=2):
    return meshioplusplus.box(dims=(n, n, n), cell_type="tetra")


def _volume(mesh):
    return meshioplusplus.compute_stats(mesh)["signed_volume"]


@needs_core
def test_dual_of_a_box_has_one_cell_per_vertex_and_its_volume():
    mesh = _box(3)
    out, report = tet_to_poly(mesh, return_report=True)
    assert [c.type for c in out.cells] == ["polyhedron"]
    assert len(out.cells[0].data) == len(mesh.points)
    assert _volume(out) == pytest.approx(1.0, abs=1e-12)
    # 6 * 4 face vertices x 1 + 12 * 2 edge vertices x 2 + 8 corners x 3
    assert report["num_boundary_faces"] == 96
    assert report["num_split_vertices"] == 0
    assert len(out.points) == len(mesh.points) + report["num_new_points"]
    np.testing.assert_array_equal(
        np.sort(report["vertex_map"]), np.arange(len(mesh.points))
    )
    # the input points keep their indices
    np.testing.assert_array_equal(out.points[: len(mesh.points)], mesh.points)


@needs_core
def test_a_curved_mesh_with_inverted_tetrahedra_converts():
    # The tetrahedral ball has a few inverted (tangled) tets near its surface;
    # the dual orients every tet consistently with its neighbours, so it still
    # yields one closed cell per vertex and the same signed volume.
    mesh = meshioplusplus.ball(n=6, cell_type="tetra")
    out, report = tet_to_poly(mesh, return_report=True)
    assert len(out.cells[0].data) == len(mesh.points)
    assert report["num_split_vertices"] == 0
    assert _volume(out) == pytest.approx(_volume(mesh), rel=1e-9)
    coarse, smooth = tet_to_poly(mesh, feature_angle=30.0, return_report=True)
    assert smooth["num_boundary_faces"] < report["num_boundary_faces"]


@needs_core
def test_feature_angle_fuses_boundary_runs():
    mesh = _box(2)
    _, sharp = tet_to_poly(mesh, return_report=True)
    _, smooth = tet_to_poly(mesh, feature_angle=95.0, return_report=True)
    assert sharp["num_boundary_faces"] == 54
    assert smooth["num_boundary_faces"] == 26


@needs_core
def test_data_rules_and_regions():
    pts = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0.5, 0.5, 0.5]], float)
    mesh = meshioplusplus.Mesh(
        pts,
        [("tetra", [[0, 1, 2, 4], [0, 2, 3, 4]]), ("triangle", [[0, 1, 2], [0, 2, 3]])],
        point_data={"x": pts[:, 0].copy()},
        cell_data={
            "w": [np.array([1.0, 3.0]), np.array([10.0, 20.0])],
            "id": [np.array([7, 9]), np.array([1, 2])],
        },
    )
    out = tet_to_poly(mesh, point_data_to_cell_data=True, record_vertex_ids=True)
    assert [c.type for c in out.cells] == ["polyhedron", "polygon"]
    vid = out.cell_data["dual:vertex"][0]
    w = out.cell_data["w"][0]
    ids = out.cell_data["id"][0]
    for v, wv, iv in zip(vid, w, ids):
        if v == 1:
            assert (wv, iv) == (1.0, 7)
        elif v == 3:
            assert (wv, iv) == (3.0, 9)
        else:
            assert (wv, iv) == (2.0, 7)
    np.testing.assert_array_equal(out.cell_data["x"][0], pts[vid, 0])
    # the polygons over the base carry the triangle block's data by area
    for v, wv in zip(out.cell_data["dual:vertex"][1], out.cell_data["w"][1]):
        assert wv == {1: 10.0, 3: 20.0}.get(int(v), 15.0)
    assert len(out.point_data["x"]) == len(out.points)


@needs_core
def test_two_tetra_blocks_split_the_interface_vertices():
    box = _box(2)
    tets = box.cells[0].data
    cx = box.points[tets].mean(axis=1)[:, 0]
    mesh = meshioplusplus.Mesh(
        box.points, [("tetra", tets[cx < 0.5]), ("tetra", tets[cx >= 0.5])]
    )
    out, report = tet_to_poly(mesh, return_report=True)
    assert len(out.cells) == 2
    assert report["num_split_vertices"] == 9
    assert sum(len(c.data) for c in out.cells) == 27 + 9
    assert _volume(out) == pytest.approx(1.0, abs=1e-12)


@needs_core
def test_dual_survives_a_vtu_round_trip(tmp_path):
    out = tet_to_poly(_box(2))
    path = str(tmp_path / "dual.vtu")
    meshioplusplus.write(path, out)
    back = meshioplusplus.read(path)
    assert sum(len(c.data) for c in back.cells) == 27
    assert _volume(back) == pytest.approx(1.0, abs=1e-12)


@needs_core
@pytest.mark.parametrize(
    "mesh, match",
    [
        (meshioplusplus.box(dims=(1, 1, 1)), "not linear tetra"),
        (meshioplusplus.box(dims=(1, 1, 1), cell_type="tetra", order=2), "linearize"),
        (meshioplusplus.icosa_sphere(1.0, 1), "no tetra block"),
    ],
)
def test_refusals(mesh, match):
    with pytest.raises(ValueError, match=match):
        tet_to_poly(mesh)


@needs_core
def test_feature_angle_out_of_range():
    with pytest.raises(ValueError, match="mFeatureAngleDeg"):
        tet_to_poly(_box(1), feature_angle=180.0)


@needs_core
def test_poly_to_tet_recovers_tetrahedra_stored_as_polyhedra_exactly():
    pts = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]], float)
    mesh = meshioplusplus.Mesh(
        pts,
        [
            (
                "polyhedron4",
                [
                    [[0, 2, 1], [0, 1, 3], [1, 2, 3], [0, 3, 2]],
                    [[1, 2, 3], [1, 4, 2], [2, 4, 3], [3, 4, 1]],  # inward-wound
                ],
            )
        ],
    )
    out, report = poly_to_tet(mesh, return_report=True)
    assert report["num_exact_tets"] == 2
    assert report["num_fanned_cells"] == 0
    assert len(out.points) == 5
    assert [c.type for c in out.cells] == ["tetra"]
    assert sorted(out.cells[0].data[0]) == [0, 1, 2, 3]
    p = out.points[out.cells[0].data]
    signed = np.einsum(
        "ij,ij->i", np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), p[:, 3] - p[:, 0]
    )
    assert np.all(signed > 0)


@needs_core
def test_round_trip_through_the_dual_conserves_volume():
    mesh = _box(2)
    back, report = poly_to_tet(tet_to_poly(mesh), return_report=True)
    assert report["num_exact_tets"] == 0
    assert report["num_fanned_cells"] == 27
    assert [c.type for c in back.cells] == ["tetra"]
    assert _volume(back) == pytest.approx(1.0, abs=1e-12)


@needs_core
def test_poly_to_tet_matches_simplexify_on_non_tetrahedral_polyhedra():
    dual = tet_to_poly(_box(1))
    a = poly_to_tet(dual, record_parent_ids=True)
    b = meshioplusplus.convert_cells(dual, mode="simplexify", record_parent_ids=True)
    assert meshioplusplus.meshes_equal(a, b)


def test_no_core_raises(monkeypatch):
    import builtins

    real_import = builtins.__import__

    def fake_import(name, globals=None, locals=None, fromlist=(), level=0):
        if fromlist and "_core" in fromlist and level == 1:
            raise ImportError("no core")
        return real_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    with pytest.raises(NotImplementedError, match="C\\+\\+-core only"):
        tet_to_poly(None)
    with pytest.raises(NotImplementedError, match="C\\+\\+-core only"):
        poly_to_tet(None)

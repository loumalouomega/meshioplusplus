"""Tests for ``compute_normals`` -- point and cell normals of a surface, with
optional vertex splitting at creases. Fixtures mirror
``tests/cpp/test_compute_normals.cpp``; the numpy twin in ``_normals.py`` is
pinned against the compiled core."""

import numpy as np
import pytest

import meshioplusplus as mio
from meshioplusplus import _normals as nrm
from meshioplusplus import compute_normals

from .test_curvature import icosphere, open_cylinder

try:
    from meshioplusplus import _core
except ImportError:  # pragma: no cover - a pure-Python build
    _core = None

needs_core = pytest.mark.skipif(_core is None, reason="needs the compiled core")


def cube_quads():
    pts = [
        [0, 0, 0],
        [1, 0, 0],
        [1, 1, 0],
        [0, 1, 0],
        [0, 0, 1],
        [1, 0, 1],
        [1, 1, 1],
        [0, 1, 1],
    ]
    cells = [
        [0, 3, 2, 1],
        [4, 5, 6, 7],
        [0, 1, 5, 4],
        [3, 7, 6, 2],
        [0, 4, 7, 3],
        [1, 2, 6, 5],
    ]
    return mio.Mesh(np.array(pts, float), [("quad", np.array(cells))])


def _twin(mesh, **kw):
    return nrm._compute_normals_py(
        mesh,
        kw.get("point_normals", True),
        kw.get("cell_normals", False),
        kw.get("weight", "angle"),
        kw.get("split_angle"),
        kw.get("record_parent_ids", False),
        kw.get("region", ""),
    )


def _same(a, b):
    """The compiled core's result and the twin's, array for array."""
    assert np.array_equal(a.points, b.points)
    for ca, cb in zip(a.cells, b.cells):
        assert ca.type == cb.type
        assert np.array_equal(np.asarray(ca.data), np.asarray(cb.data))
    assert sorted(a.point_data) == sorted(b.point_data)
    for k in a.point_data:
        np.testing.assert_array_equal(a.point_data[k], b.point_data[k], err_msg=k)
    assert sorted(a.cell_data) == sorted(b.cell_data)
    for k in a.cell_data:
        for x, y in zip(a.cell_data[k], b.cell_data[k]):
            np.testing.assert_array_equal(x, y, err_msg=k)


def test_flat_square_has_up_normals():
    m = mio.Mesh(
        np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]], float),
        [("triangle", np.array([[0, 1, 2], [0, 2, 3]]))],
    )
    out = compute_normals(m, cell_normals=True)
    np.testing.assert_array_equal(out.point_data["normals"], [[0, 0, 1]] * 4)
    np.testing.assert_array_equal(out.cell_data["normals"][0], [[0, 0, 1]] * 2)
    assert "normals" not in m.point_data  # the input is never modified


def test_two_dimensional_points_get_an_up_normal():
    m = mio.Mesh(
        np.array([[0, 0], [1, 0], [1, 1]], float), [("triangle", np.array([[0, 1, 2]]))]
    )
    out = compute_normals(m)
    assert out.point_data["normals"].shape == (3, 3)
    np.testing.assert_array_equal(out.point_data["normals"][:, 2], 1.0)


def test_splitting_a_cube_gives_axis_aligned_normals():
    out, rep = compute_normals(
        cube_quads(),
        split_angle=30,
        cell_normals=True,
        record_parent_ids=True,
        return_report=True,
    )
    assert len(out.points) == 24
    assert rep["num_added_points"] == 16 and rep["num_split_points"] == 8
    n = out.point_data["normals"]
    assert np.all(np.sort(np.abs(n), axis=1)[:, :2] == 0.0)
    assert np.all(np.abs(n).max(axis=1) == 1.0)
    # every cell's corners sit on points carrying the cell's own normal
    for c, row in enumerate(out.cells[0].data):
        np.testing.assert_array_equal(n[row], [out.cell_data["normals"][0][c]] * 4)
    parent = out.point_data["normals:parent_point"]
    np.testing.assert_array_equal(parent[:8], np.arange(8))
    np.testing.assert_array_equal(out.points[8:], cube_quads().points[parent[8:]])


def test_a_full_split_angle_on_a_closed_manifold_is_the_unsplit_result():
    smooth = compute_normals(cube_quads())
    wide = compute_normals(cube_quads(), split_angle=180)
    assert len(wide.points) == 8
    np.testing.assert_array_equal(
        wide.point_data["normals"], smooth.point_data["normals"]
    )


@pytest.mark.parametrize("weight", ["angle", "area"])
def test_a_sphere_has_radial_normals(weight):
    m = icosphere(3, 2.0)
    out = compute_normals(m, weight=weight)
    np.testing.assert_allclose(
        out.point_data["normals"],
        m.points / np.linalg.norm(m.points, axis=1)[:, None],
        atol=1e-2,
    )


def test_a_sphere_below_the_split_angle_is_not_split():
    out, rep = compute_normals(icosphere(3), split_angle=30, return_report=True)
    assert rep["num_added_points"] == 0
    assert len(out.points) == len(icosphere(3).points)


def test_an_inconsistent_pair_is_reported_and_always_split():
    m = mio.Mesh(
        np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]], float),
        [("triangle", np.array([[0, 1, 2], [3, 1, 2]]))],
    )
    _out, rep = compute_normals(m, split_angle=180, return_report=True)
    assert rep["quality"]["inconsistent_pairs"] == 1
    assert rep["num_added_points"] == 2


def test_lines_and_vertices_are_ignored():
    m = mio.Mesh(
        np.array([[0, 0, 0], [1, 0, 0]], float), [("line", np.array([[0, 1]]))]
    )
    out, rep = compute_normals(m, cell_normals=True, return_report=True)
    assert np.isnan(out.point_data["normals"]).all()
    assert np.isnan(out.cell_data["normals"][0]).all()
    assert rep["num_isolated"] == 2


def test_volumes_and_higher_order_surfaces_are_refused_by_name():
    tet = mio.Mesh(
        np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], float),
        [("tetra", np.array([[0, 1, 2, 3]]))],
    )
    with pytest.raises(ValueError, match="extract_surface"):
        compute_normals(tet)
    t6 = mio.Mesh(
        np.array(
            [[0, 0, 0], [2, 0, 0], [0, 2, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]], float
        ),
        [("triangle6", np.array([[0, 1, 2, 3, 4, 5]]))],
    )
    with pytest.raises(ValueError, match="linearize"):
        compute_normals(t6)


def test_bad_arguments_are_rejected():
    m = cube_quads()
    with pytest.raises(ValueError, match="split angle"):
        compute_normals(m, split_angle=181)
    with pytest.raises(ValueError, match="split angle"):
        compute_normals(m, split_angle=-1)
    with pytest.raises(ValueError, match="weight"):
        compute_normals(m, weight="bogus")
    with pytest.raises(ValueError, match="no cell region named 'nowhere'"):
        compute_normals(m, region="nowhere")


def test_a_region_restricts_the_surface():
    m = mio.Mesh(
        np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]], float),
        [("triangle", np.array([[0, 1, 2], [0, 2, 3]]))],
        regions=[mio.Region("first", "cell", [0])],
    )
    out, rep = compute_normals(m, region="first", return_report=True)
    n = out.point_data["normals"]
    np.testing.assert_array_equal(n[:3, 2], 1.0)
    assert np.isnan(n[3]).all()
    assert rep["num_isolated"] == 1


def test_splitting_carries_data_and_regions_along():
    m = cube_quads()
    m.point_data["T"] = 10.0 * np.arange(8)
    m.cell_data["id"] = [np.arange(6.0)]
    m.regions = [mio.Region("corner", "point", [0]), mio.Region("top", "cell", [1])]
    out = compute_normals(m, split_angle=30, record_parent_ids=True)
    parent = out.point_data["normals:parent_point"]
    np.testing.assert_array_equal(out.point_data["T"], 10.0 * parent)
    by_name = {r.name: r for r in out.regions}
    assert len(by_name["corner"].entries) == 3  # the point and its two copies
    assert list(by_name["top"].entries) == [1]
    np.testing.assert_array_equal(out.cell_data["id"][0], np.arange(6.0))


def test_an_existing_normals_array_is_replaced():
    m = cube_quads()
    m.point_data["normals"] = np.full((8, 3), 7.0)
    out = compute_normals(m)
    assert not (out.point_data["normals"] == 7.0).any()


def test_normals_reach_the_point_cloud_formats(tmp_path):
    out = compute_normals(cube_quads(), split_angle=30)
    path = tmp_path / "cube.xyz"
    mio.write(
        str(path),
        mio.Mesh(
            out.points,
            [("vertex", np.arange(24)[:, None])],
            point_data={"normals": out.point_data["normals"]},
        ),
    )
    rows = [ln.split() for ln in path.read_text().split("\n") if ln and ln[0] != "#"]
    assert len(rows[0]) == 6


# --------------------------------------------------------------------------- #
# the numpy twin against the compiled core                                     #
# --------------------------------------------------------------------------- #
def _twin_cases():
    ragged = mio.Mesh(
        np.array(
            [[0, 0, 0], [1, 0, 0], [1.5, 1, 0], [0.5, 2, 0], [-0.5, 1, 0], [0, 0, 1]],
            float,
        ),
        [("polygon", [[0, 1, 2, 3, 4]]), ("triangle", np.array([[0, 5, 4]]))],
    )
    return {
        "cube": cube_quads(),
        "sphere": icosphere(2, 1.5),
        "cylinder": open_cylinder(),
        "ragged": ragged,
    }


@needs_core
@pytest.mark.parametrize("name", ["cube", "sphere", "cylinder", "ragged"])
@pytest.mark.parametrize("weight", ["angle", "area"])
@pytest.mark.parametrize("split", [None, 30.0, 180.0, 0.0])
def test_the_twin_matches_the_core(name, weight, split):
    m = _twin_cases()[name]
    kw = dict(
        weight=weight, split_angle=split, cell_normals=True, record_parent_ids=True
    )
    core, core_rep = compute_normals(m, return_report=True, **kw)
    twin, twin_rep = _twin(m, **kw)
    _same(core, twin)
    assert core_rep == twin_rep


@needs_core
def test_the_twin_matches_the_core_on_a_region_and_with_data():
    m = _twin_cases()["cube"]
    m.point_data["T"] = np.arange(8.0)
    m.regions = [mio.Region("top", "cell", [1, 2, 3]), mio.Region("p", "point", [0, 6])]
    kw = dict(split_angle=30.0, region="top", record_parent_ids=True)
    core, core_rep = compute_normals(m, return_report=True, **kw)
    twin, twin_rep = _twin(m, **kw)
    _same(core, twin)
    assert core_rep == twin_rep
    assert {r.name: r.entries.tolist() for r in core.regions} == {
        r.name: r.entries.tolist() for r in twin.regions
    }


# --------------------------------------------------------------------------- #
# the CLI and the pipeline step                                                #
# --------------------------------------------------------------------------- #
def test_cli_normals_verb(tmp_path):
    src = tmp_path / "cube.vtu"
    dst = tmp_path / "cube_n.vtu"
    mio.write(str(src), cube_quads())
    mio._cli.main(
        ["normals", str(src), str(dst), "--split-angle", "30", "--cell", "-q"]
    )
    out = mio.read(str(dst))
    assert len(out.points) == 24
    assert "normals" in out.point_data and "normals" in out.cell_data

    with pytest.raises(ValueError, match="split angle"):
        mio._cli.main(["normals", str(src), str(dst), "--split-angle", "400", "-q"])


def _pipeline_settings(tmp_path, step):
    src = tmp_path / "cube.vtu"
    mio.write(str(src), cube_quads())
    return {
        "Version": 1,
        "Input": {"Path": str(src)},
        "Operations": [dict(Op="Normals", **step)],
        "Output": {"Path": str(tmp_path / "out.vtu")},
    }


def test_pipeline_normals_step(tmp_path):
    settings = _pipeline_settings(tmp_path, {"SplitAngle": 30, "CellNormals": True})
    report = mio.run_pipeline(settings)
    assert report["steps"][0]["op"] == "Normals"
    assert report["steps"][0]["NumAddedPoints"] == 16
    out = mio.read(settings["Output"]["Path"])
    assert len(out.points) == 24 and "normals" in out.point_data

    # an absent SplitAngle means one smooth normal per point
    smooth = _pipeline_settings(tmp_path, {})
    mio.run_pipeline(smooth)
    assert len(mio.read(smooth["Output"]["Path"]).points) == 8

    with pytest.raises(ValueError):
        mio.run_pipeline(_pipeline_settings(tmp_path, {"Bogus": 1}))


@needs_core
def test_pipeline_normals_step_matches_the_c_engine(tmp_path):
    import json

    if not hasattr(_core, "run_pipeline_json"):
        pytest.skip("this build carries no C++ pipeline engine")
    settings = _pipeline_settings(tmp_path, {"SplitAngle": 30, "RecordParentIds": True})
    mio.run_pipeline(settings)
    py = mio.read(settings["Output"]["Path"])
    settings["Output"]["Path"] = str(tmp_path / "out_cpp.vtu")
    try:
        _core.run_pipeline_json(json.dumps(settings))
    except RuntimeError as exc:
        if "no JSON parser" in str(exc):
            pytest.skip("this build carries no JSON parser for the C++ engine")
        raise
    cpp = mio.read(settings["Output"]["Path"])
    np.testing.assert_array_equal(py.points, cpp.points)
    np.testing.assert_array_equal(py.point_data["normals"], cpp.point_data["normals"])

"""Tests for the mesh comparison ("diff") operation."""

import copy

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import Mesh, diff, meshes_equal
from meshioplusplus._cli import main as cli_main

from . import helpers


def _tri():
    # Deepcopy shared fixtures before any mutation (some Python writers mutate
    # their input; see the shared-fixture caveat).
    return copy.deepcopy(helpers.tri_mesh)


def test_identical():
    a = _tri()
    b = _tri()
    report = diff(a, b)
    assert report["verdict"] == "identical"
    assert meshes_equal(a, b) is True
    assert report["point_data"]["only_in_a"] == []
    assert not report["block_count_mismatch"]


def test_perturb_below_tolerance():
    a = _tri()
    b = _tri()
    b.points = np.array(b.points, dtype=float)
    b.points[2, 1] += 1e-10
    report = diff(a, b, atol=1e-8, rtol=0.0)
    assert report["verdict"] == "equal within tolerance"
    assert meshes_equal(a, b, atol=1e-8, rtol=0.0) is True
    assert report["points"]["max_abs_error"] > 0.0
    assert report["points"]["num_exceeding"] == 0


def test_perturb_above_tolerance_worst_index():
    a = _tri()
    b = _tri()
    ncols = np.asarray(a.points).shape[1]
    b.points = np.array(b.points, dtype=float)
    b.points[2, 1] += 1.0
    report = diff(a, b, atol=1e-8, rtol=1e-6)
    assert report["verdict"] == "different"
    assert not meshes_equal(a, b, atol=1e-8, rtol=1e-6)
    assert report["points"]["num_exceeding"] >= 1
    # Worst offender is point 2, coordinate 1 -> flat index 2*ncols + 1.
    assert report["points"]["worst_index"] == 2 * ncols + 1


def test_different_cell_count():
    a = _tri()
    b = Mesh(np.asarray(a.points), [("triangle", np.asarray(a.cells[0].data)[:1])])
    report = diff(a, b)
    assert report["verdict"] == "different"
    assert report["blocks"][0]["count_mismatch"]
    assert report["blocks"][0]["count_a"] != report["blocks"][0]["count_b"]


def test_different_cell_type():
    a = Mesh(
        [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]],
        [("quad", [[0, 1, 2, 3]])],
    )
    b = Mesh(
        [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]],
        [("triangle", [[0, 1, 2]])],
    )
    report = diff(a, b)
    assert report["verdict"] == "different"
    assert report["blocks"][0]["type_mismatch"]
    assert report["blocks"][0]["type_a"] == "quad"
    assert report["blocks"][0]["type_b"] == "triangle"


def test_different_connectivity():
    a = _tri()
    conn = np.asarray(a.cells[0].data).copy()
    conn[-1] = conn[-1][::-1]  # reverse the last triangle's node order
    b = Mesh(np.asarray(a.points), [("triangle", conn)])
    report = diff(a, b)
    assert report["verdict"] == "different"
    assert report["blocks"][0]["conn_mismatch_count"] == 1
    assert report["blocks"][0]["first_mismatches"] == [len(conn) - 1]


def test_added_point_data_key():
    a = _tri()
    b = _tri()
    a.point_data["temperature"] = np.arange(len(a.points), dtype=float)
    report = diff(a, b)
    assert report["verdict"] == "different"
    assert report["point_data"]["only_in_a"] == ["temperature"]
    assert report["point_data"]["only_in_b"] == []


def test_removed_cell_data_key():
    a = _tri()
    b = _tri()
    b.cell_data["region"] = [np.zeros(len(a.cells[0].data), dtype=np.int32)]
    report = diff(a, b)
    assert report["verdict"] == "different"
    assert report["cell_data"]["only_in_b"] == ["region"]


def test_point_data_values_within_tolerance():
    a = _tri()
    b = _tri()
    n = len(a.points)
    a.point_data["t"] = np.arange(n, dtype=float)
    b.point_data["t"] = np.arange(n, dtype=float) + 1e-10
    report = diff(a, b, atol=1e-8, rtol=0.0)
    assert report["verdict"] == "equal within tolerance"
    shared = report["point_data"]["shared"]
    assert len(shared) == 1
    assert shared[0]["name"] == "t"
    assert shared[0]["max_abs_error"] > 0.0


@pytest.mark.parametrize(
    "writer_fmt,ext,atol",
    [
        ("vtu", ".vtu", 1e-12),
        ("vtk", ".vtk", 1e-6),
    ],
)
def test_roundtrip_regression(tmp_path, writer_fmt, ext, atol):
    # The headline CI use case: write to another format, read back, assert equal.
    orig = copy.deepcopy(helpers.tet_mesh)
    path = tmp_path / f"mesh{ext}"
    meshioplusplus.write(str(path), copy.deepcopy(orig), file_format=writer_fmt)
    back = meshioplusplus.read(str(path), file_format=writer_fmt)
    assert meshes_equal(orig, back, atol=atol), diff(orig, back, atol=atol)


def test_unordered_shuffled_nodes():
    a = _tri()
    pts = np.asarray(a.points, dtype=float)
    conn = np.asarray(a.cells[0].data)
    n = len(pts)
    rng = np.random.default_rng(0)
    perm = rng.permutation(n)  # old -> new
    inv = np.empty(n, dtype=np.int64)
    inv[perm] = np.arange(n)
    pts_b = pts[inv]  # new node i holds old node inv[i]
    conn_b = perm[conn]  # remap connectivity to new node numbering
    b = Mesh(pts_b, [("triangle", conn_b)])

    # Ordered comparison sees the shuffle as different...
    assert diff(a, b)["verdict"] == "different"
    # ...unordered proximity matching sees them as equal.
    report = diff(a, b, unordered=True)
    assert report["verdict"] == "identical"
    assert not report["correspondence_failed"]
    assert report["blocks"][0]["conn_mismatch_count"] == 0
    assert meshes_equal(a, b, unordered=True)


def test_unordered_correspondence_fails():
    a = _tri()
    pts = np.asarray(a.points, dtype=float) + 100.0  # move far away
    b = Mesh(pts, [("triangle", np.asarray(a.cells[0].data))])
    report = diff(a, b, unordered=True, atol=1e-6, rtol=1e-6)
    assert report["correspondence_failed"]
    assert report["verdict"] == "different"


def test_sets_reported():
    a = _tri()
    b = _tri()
    a.point_sets = {"left": np.array([0, 3])}
    b.cell_sets = {"all": [np.array([0, 1])]}
    report = diff(a, b)
    assert report["verdict"] == "different"
    assert report["sets"]["point_sets"]["only_in_a"] == ["left"]
    assert report["sets"]["cell_sets"]["only_in_b"] == ["all"]
    # meshes_equal ignores sets.
    a2 = _tri()
    a2.point_sets = {"left": np.array([0, 3])}
    assert meshes_equal(a2, _tri())


def test_cpp_matches_python():
    from meshioplusplus import _diff

    a = copy.deepcopy(helpers.tet_mesh)
    b = copy.deepcopy(helpers.tet_mesh)
    b.points = np.array(b.points, dtype=float)
    b.points[1, 0] += 1e-3
    cpp = _diff.diff(a, b, atol=1e-8, rtol=1e-6)
    py = _diff._diff_py(a, b, 1e-8, 1e-6, False, 10)
    assert cpp["verdict"] == py["verdict"]
    assert cpp["points"]["worst_index"] == py["points"]["worst_index"]
    assert cpp["points"]["num_exceeding"] == py["points"]["num_exceeding"]
    np.testing.assert_allclose(
        cpp["points"]["max_abs_error"], py["points"]["max_abs_error"]
    )


def test_cli_exit_codes(tmp_path):
    a = copy.deepcopy(helpers.tri_mesh)
    pa = tmp_path / "a.vtu"
    pb = tmp_path / "b.vtu"
    meshioplusplus.write(str(pa), copy.deepcopy(a))
    b = copy.deepcopy(a)
    b.points = np.array(b.points, dtype=float)
    b.points[2, 1] += 1e-6
    meshioplusplus.write(str(pb), b)

    assert cli_main(["diff", str(pa), str(pa)]) == 0
    assert cli_main(["diff", str(pa), str(pb), "--atol", "1e-4"]) == 0
    assert cli_main(["diff", str(pa), str(pb), "--atol", "1e-9"]) == 1
    # --exact fails even on tolerated drift.
    assert cli_main(["diff", str(pa), str(pb), "--atol", "1e-4", "--exact"]) == 1


def test_cli_quiet(tmp_path, capsys):
    a = copy.deepcopy(helpers.tri_mesh)
    pa = tmp_path / "a.vtu"
    meshioplusplus.write(str(pa), copy.deepcopy(a))
    rc = cli_main(["diff", str(pa), str(pa), "--quiet"])
    captured = capsys.readouterr()
    assert rc == 0
    assert captured.out == ""

"""Tests for the mesh merge/weld operation (meshioplusplus.merge)."""

import copy

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import Mesh, merge, meshes_equal
from meshioplusplus._cli import main as cli_main


def _block_a():
    # unit square split into 2 triangles
    return Mesh(
        [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]],
        [("triangle", [[0, 1, 2], [0, 2, 3]])],
    )


def _block_b():
    # adjacent square sharing the x=1 edge (points (1,0,0), (1,1,0) coincide)
    return Mesh(
        [[1, 0, 0], [2, 0, 0], [2, 1, 0], [1, 1, 0]],
        [("triangle", [[0, 1, 2], [0, 2, 3]])],
    )


def _ncells(mesh):
    return sum(len(cb.data) for cb in mesh.cells)


# --------------------------------------------------------------------------- #
# concatenation                                                               #
# --------------------------------------------------------------------------- #
def test_concatenation_counts():
    a, b = _block_a(), _block_b()
    m = merge([a, b])
    assert len(m.points) == len(a.points) + len(b.points)
    assert _ncells(m) == _ncells(a) + _ncells(b)
    # same-type blocks merged
    assert len(m.cells) == 1
    # connectivity indices remain valid
    for cb in m.cells:
        conn = np.asarray(cb.data)
        assert conn.min() >= 0
        assert conn.max() < len(m.points)


def test_offset_correctness():
    a, b = _block_a(), _block_b()
    m = merge([a, b])
    # B's first cell (local {0,1,2}) is offset by 4 -> {4,5,6}; node 4 must be
    # B's original point 0 = (1,0,0).
    conn = np.asarray(m.cells[0].data)
    assert list(conn[2]) == [4, 5, 6]
    np.testing.assert_allclose(m.points[4], [1, 0, 0])


def test_source_tag():
    a, b = _block_a(), _block_b()
    m = merge([a, b])
    assert "source_mesh_id" in m.cell_data
    tag = np.concatenate([np.asarray(x).ravel() for x in m.cell_data["source_mesh_id"]])
    assert len(tag) == _ncells(m)
    assert list(tag) == [0, 0, 1, 1]


def test_no_source_tag():
    m = merge([_block_a(), _block_b()], source_tag=False)
    assert "source_mesh_id" not in m.cell_data


# --------------------------------------------------------------------------- #
# welding                                                                     #
# --------------------------------------------------------------------------- #
def test_weld_shared_interface():
    m = merge([_block_a(), _block_b()], weld=True, atol=1e-9)
    # 8 input points, 2 coincident pairs -> 6 output points
    assert len(m.points) == 6
    assert _ncells(m) == 4
    for cb in m.cells:
        assert np.asarray(cb.data).max() < 6


def test_weld_tolerance():
    a = _block_a()
    # B shifted +1e-3 in x: near but not exactly coincident with A's right edge
    b = Mesh(
        [[1.001, 0, 0], [2, 0, 0], [2, 1, 0], [1.001, 1, 0]],
        [("triangle", [[0, 1, 2], [0, 2, 3]])],
    )
    # below tolerance -> nothing welds
    assert len(merge([a, b], weld=True, atol=1e-6).points) == 8
    # above tolerance -> the 2 near-coincident pairs weld
    assert len(merge([a, b], weld=True, atol=1e-2).points) == 6


def test_drop_duplicate_cells():
    a = _block_a()
    m = merge([a, copy.deepcopy(a)], weld=True, atol=1e-9, drop_duplicate_cells=True)
    assert len(m.points) == 4
    assert _ncells(m) == 2


# --------------------------------------------------------------------------- #
# data policy                                                                 #
# --------------------------------------------------------------------------- #
def test_data_policy_intersection():
    a = _block_a()
    a.point_data = {
        "T": np.array([1.0, 2.0, 3.0, 4.0]),
        "onlyA": np.array([9, 8, 7, 6]),
    }
    b = _block_b()
    b.point_data = {"T": np.array([5.0, 6.0, 7.0, 8.0])}
    m = merge([a, b])  # intersection is the default
    assert "T" in m.point_data
    assert "onlyA" not in m.point_data
    assert len(np.asarray(m.point_data["T"])) == 8


def test_data_policy_fill():
    a = _block_a()
    a.point_data = {
        "T": np.array([1.0, 2.0, 3.0, 4.0]),
        "onlyA": np.array([9.0, 8.0, 7.0, 6.0]),
    }
    b = _block_b()
    b.point_data = {"T": np.array([5.0, 6.0, 7.0, 8.0])}
    m = merge([a, b], data_policy="fill")
    assert "onlyA" in m.point_data
    oa = np.asarray(m.point_data["onlyA"])
    assert np.isnan(oa[4:]).all()
    assert np.allclose(oa[:4], [9, 8, 7, 6])


def test_concatenated_data_length():
    a = _block_a()
    a.point_data = {"T": np.array([1.0, 2.0, 3.0, 4.0])}
    b = _block_b()
    b.point_data = {"T": np.array([5.0, 6.0, 7.0, 8.0])}
    m = merge([a, b])
    np.testing.assert_allclose(np.asarray(m.point_data["T"]), [1, 2, 3, 4, 5, 6, 7, 8])


def test_bad_data_policy():
    with pytest.raises(ValueError):
        merge([_block_a(), _block_b()], data_policy="bogus")


# --------------------------------------------------------------------------- #
# sets (namespaced by source id, remapped)                                    #
# --------------------------------------------------------------------------- #
def test_point_sets_namespaced_and_remapped():
    a = _block_a()
    a.point_sets = {"left": np.array([0, 3])}
    b = _block_b()
    b.point_sets = {"left": np.array([0, 3])}
    m = merge([a, b])
    assert set(m.point_sets.keys()) == {"0:left", "1:left"}
    # B's local indices 0,3 -> global 4,7 in the concatenation
    assert list(m.point_sets["1:left"]) == [4, 7]


def test_cell_sets_remapped():
    a = _block_a()
    a.cell_sets = {"g": [np.array([0])]}
    b = _block_b()
    b.cell_sets = {"g": [np.array([1])]}
    m = merge([a, b])
    assert set(m.cell_sets.keys()) == {"0:g", "1:g"}
    # A's cell 0 stays at output cell 0
    assert list(m.cell_sets["0:g"][0]) == [0]
    # B's cell 1 -> output cell 3
    assert list(m.cell_sets["1:g"][0]) == [3]


def test_field_data_namespaced():
    a = _block_a()
    a.field_data = {"mat": np.array([1, 2]), "onlyA": np.array([9])}
    b = _block_b()
    b.field_data = {"mat": np.array([3, 4])}
    m = merge([a, b])
    # colliding key namespaced, unique key kept
    assert "0:mat" in m.field_data and "1:mat" in m.field_data
    assert "onlyA" in m.field_data


# --------------------------------------------------------------------------- #
# C++ core vs pure-python reference                                           #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("weld", [False, True])
@pytest.mark.parametrize("policy", ["intersection", "fill"])
def test_cpp_matches_python(weld, policy):
    from meshioplusplus import _merge

    a = _block_a()
    a.point_data = {"T": np.array([1.0, 2.0, 3.0, 4.0])}
    a.cell_data = {"p": [np.array([10, 11])]}
    b = _block_b()
    b.point_data = {"T": np.array([5.0, 6.0, 7.0, 8.0])}
    b.cell_data = {"p": [np.array([20, 21])]}

    cpp = merge([a, b], weld=weld, atol=1e-9, data_policy=policy)
    py, _, _ = _merge._merge_py([a, b], weld, 1e-9, True, policy, False)
    assert meshes_equal(cpp, py, atol=1e-12), meshioplusplus.diff(cpp, py, atol=1e-12)


# --------------------------------------------------------------------------- #
# round-trip                                                                   #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("weld", [False, True])
def test_roundtrip_write_read(tmp_path, weld):
    m = merge([_block_a(), _block_b()], weld=weld, atol=1e-9)
    path = tmp_path / "merged.vtu"
    meshioplusplus.write(str(path), copy.deepcopy(m))
    back = meshioplusplus.read(str(path))
    assert meshes_equal(m, back, atol=1e-10), meshioplusplus.diff(m, back)


# --------------------------------------------------------------------------- #
# CLI                                                                          #
# --------------------------------------------------------------------------- #
def test_cli_merge(tmp_path, capsys):
    pa = tmp_path / "a.vtu"
    pb = tmp_path / "b.vtu"
    out = tmp_path / "out.vtu"
    meshioplusplus.write(str(pa), _block_a())
    meshioplusplus.write(str(pb), _block_b())

    rc = cli_main(["merge", str(pa), str(pb), str(out)])
    assert rc == 0
    captured = capsys.readouterr()
    assert "merged 2 meshes" in captured.out
    read_back = meshioplusplus.read(str(out))
    assert len(read_back.points) == 8

    # weld path + summary line
    outw = tmp_path / "outw.vtu"
    rc = cli_main(["merge", str(pa), str(pb), str(outw), "--weld", "--atol", "1e-8"])
    assert rc == 0
    captured = capsys.readouterr()
    assert "points welded: 2" in captured.out
    assert len(meshioplusplus.read(str(outw)).points) == 6


def test_cli_quiet(tmp_path, capsys):
    pa = tmp_path / "a.vtu"
    pb = tmp_path / "b.vtu"
    out = tmp_path / "out.vtu"
    meshioplusplus.write(str(pa), _block_a())
    meshioplusplus.write(str(pb), _block_b())
    rc = cli_main(["merge", str(pa), str(pb), str(out), "--quiet"])
    captured = capsys.readouterr()
    assert rc == 0
    assert captured.out == ""


def test_cli_too_few_files(tmp_path):
    pa = tmp_path / "a.vtu"
    meshioplusplus.write(str(pa), _block_a())
    with pytest.raises(SystemExit):
        cli_main(["merge", str(pa)])

"""Tests for ``check_conformity`` (roadmap 3.2.1).

The probes are meshes whose defects are known by construction: the hanging
nodes ``refine(closure="balanced")`` records in ``refine:hanging``, the
conforming closures, the green parents ``undo_green`` restores beside red
children, and hand-built T-junctions, overlaps and unwelded faces.
"""

import numpy as np
import pytest

import meshioplusplus as mp

pytest.importorskip("meshioplusplus._core")


def _flags(mesh):
    return np.asarray(mesh.point_data["conformity:hanging"]).reshape(-1)


def _t_junction():
    # A coarse quad beside two fine ones: the fine pair's shared node, 7, sits
    # on the coarse quad's edge.
    pts = np.array(
        [[0, 0], [1, 0], [1, 1], [0, 1], [2, 0], [2, 0.5], [2, 1], [1, 0.5]],
        dtype=float,
    )
    return mp.Mesh(
        pts, [("quad", np.array([[0, 1, 2, 3], [1, 4, 5, 7], [7, 5, 6, 2]]))]
    )


def test_a_grid_conforms():
    out, report = mp.check_conformity(mp.box((3, 3, 3)), return_report=True)
    assert report["dimension"] == 3
    assert report["conforming"] is True
    assert report["interface_mode"] == "conforming"
    assert report["num_boundary_facets"] == 6 * 9
    assert report["warnings"] == []
    assert not _flags(out).any()


def test_a_t_junction():
    out, report = mp.check_conformity(_t_junction(), return_report=True)
    assert report["dimension"] == 2
    assert report["num_t_junctions"] == 1
    assert report["num_partial_overlaps"] == 2
    assert report["interface_mode"] == "proximity"
    assert list(np.nonzero(_flags(out))[0]) == [7]
    assert _flags(out)[7] == 2


def test_the_mesh_passes_through_and_flags_are_optional():
    mesh = _t_junction()
    mesh.point_data["marker"] = np.arange(8, dtype=float)
    out = mp.check_conformity(mesh)
    assert np.asarray(out.points).tobytes() == np.asarray(mesh.points).tobytes()
    assert np.array_equal(out.cells[0].data, mesh.cells[0].data)
    assert np.array_equal(out.point_data["marker"], mesh.point_data["marker"])
    bare = mp.check_conformity(mesh, attach_flags=False)
    assert "conformity:hanging" not in bare.point_data


def test_unwelded_faces_suggest_weld():
    a = mp.box((1, 1, 1))
    b = mp.transform(a, translate=(1.0, 0.0, 0.0))
    both = mp.merge([a, b], weld=False)
    _, report = mp.check_conformity(both, return_report=True)
    assert report["num_unwelded_facets"] == 1
    assert report["num_hanging_nodes"] == 0
    assert report["interface_mode"] == "weld"


def test_quadratic_cells_conform():
    tets = mp.convert_cells(mp.box((2, 2, 2), cell_type="tetra"), mode="elevate")
    _, report = mp.check_conformity(tets, return_report=True)
    assert report["conforming"] is True


def _hex_grid(n):
    return mp.box((n, n, n), (n, n, n))


def test_balanced_refine_flags_exactly_refine_hanging():
    out = mp.refine(_hex_grid(4), cells=[21], closure="balanced")
    checked, report = mp.check_conformity(out, return_report=True)
    want = np.asarray(out.point_data["refine:hanging"]).reshape(-1)
    assert report["num_hanging_nodes"] == int(want.sum()) > 0
    assert np.array_equal(_flags(checked) != 0, want != 0)
    assert report["interface_mode"] == "proximity"


@pytest.mark.parametrize("closure", ["redgreen", "propagate"])
def test_conforming_closures_report_nothing(closure):
    out = mp.refine(_hex_grid(3), cells=[13], closure=closure)
    _, report = mp.check_conformity(out, return_report=True)
    assert report["conforming"] is True, report


def test_undo_green_leaves_hanging_nodes():
    coarse = mp.rectangle((4, 4), (4, 4), cell_type="triangle", split="up")
    fine = mp.refine(coarse, cells=[8], record_hierarchy=True, record_levels=True)
    _, before = mp.check_conformity(fine, return_report=True)
    assert before["conforming"] is True
    undone = mp.undo_green(coarse, fine)
    _, after = mp.check_conformity(undone, return_report=True)
    # The restored green parents meet red children whose mid-edge nodes now
    # hang on their edges; undo_green's orphan points are not counted.
    assert after["num_t_junctions"] > 0
    assert after["interface_mode"] == "proximity"


def test_lines_are_not_examined():
    line = mp.Mesh(np.array([[0.0, 0, 0], [1, 0, 0]]), [("line", np.array([[0, 1]]))])
    _, report = mp.check_conformity(line, return_report=True)
    assert report["dimension"] == 0
    assert len(report["warnings"]) == 1


@pytest.mark.parametrize("tolerance", [-1.0, float("inf"), float("nan")])
def test_a_bad_tolerance_is_refused(tolerance):
    with pytest.raises(ValueError, match="tolerance"):
        mp.check_conformity(_t_junction(), tolerance=tolerance)

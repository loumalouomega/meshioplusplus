"""Side regions surviving operations (v16.26.0): the C++ core against the numpy
twins in ``_side_carry.py``, and the geometric invariants a carry must keep."""

import numpy as np
import pytest

import meshioplusplus as mio
from meshioplusplus import Region
from meshioplusplus._facets import facet_nodes

_core = pytest.importorskip("meshioplusplus._core")


def _declined(*args, **kwargs):
    raise RuntimeError("declined for the test")


def _force_twin(monkeypatch, name):
    monkeypatch.setattr(_core, name, _declined)


def _hex_pair():
    """Two unit hexahedra side by side, with a side region on their floor."""
    pts = np.array(
        [(x, y, z) for z in (0, 1) for y in (0, 1) for x in (0, 1, 2)], dtype=float
    )

    def p(x, y, z):
        return z * 6 + y * 3 + x

    hexes = [
        [p(x, 0, 0), p(x + 1, 0, 0), p(x + 1, 1, 0), p(x, 1, 0)]
        + [p(x, 0, 1), p(x + 1, 0, 1), p(x + 1, 1, 1), p(x, 1, 1)]
        for x in (0, 1)
    ]
    mesh = mio.Mesh(pts, [("hexahedron", np.array(hexes))])
    # Face 4 of a hexahedron is its z = 0 face (cell_faces order).
    mesh.regions = [Region("floor", "side", np.array([[0, 4], [1, 4]]))]
    return mesh


def _tri_strip():
    """A 4 x 3 grid of triangles with a side region on its y = 0 edge."""
    nx, ny = 4, 3
    pts = np.array([(x, y, 0.0) for y in range(ny + 1) for x in range(nx + 1)], float)

    def p(x, y):
        return y * (nx + 1) + x

    tris = []
    for y in range(ny):
        for x in range(nx):
            a, b, c, d = p(x, y), p(x + 1, y), p(x + 1, y + 1), p(x, y + 1)
            tris += [[a, b, c], [a, c, d]]
    mesh = mio.Mesh(pts, [("triangle", np.array(tris))])
    mesh.regions = [Region("bottom", "side", np.array([[2 * x, 0] for x in range(nx)]))]
    return mesh


def _side(mesh, name):
    return next(r for r in mesh.regions if r.name == name and r.kind == "side")


def _measure(mesh, name, plane_axis):
    """Total length/area of a side region, asserting every facet lies on the
    coordinate plane ``plane_axis == 0``."""
    total = 0.0
    for cell, facet in _side(mesh, name).entries.tolist():
        _, nodes = facet_nodes(mesh, cell, facet)
        pts = np.asarray(mesh.points)[nodes]
        assert np.allclose(pts[:, plane_axis], 0.0), pts
        if len(nodes) == 2:
            total += float(np.linalg.norm(pts[1] - pts[0]))
        else:
            corners = pts[
                : 4 if len(nodes) >= 8 else (3 if len(nodes) in (3, 6) else 4)
            ]
            if len(corners) == 3:
                total += 0.5 * np.linalg.norm(
                    np.cross(corners[1] - corners[0], corners[2] - corners[0])
                )
            else:
                total += 0.5 * np.linalg.norm(
                    np.cross(corners[2] - corners[0], corners[3] - corners[1])
                )
    return total


@pytest.mark.parametrize("mode", ["simplexify", "elevate", "linearize"])
def test_convert_cells_carries_the_floor_and_the_twin_agrees(mode, monkeypatch):
    mesh = _hex_pair()
    core = mio.convert_cells(mesh, mode)
    assert _measure(core, "floor", 2) == pytest.approx(2.0)
    _force_twin(monkeypatch, "convert_cells")
    twin = mio.convert_cells(mesh, mode)
    np.testing.assert_array_equal(
        _side(core, "floor").entries, _side(twin, "floor").entries
    )


def test_refine_carries_the_floor_and_the_twin_agrees(monkeypatch):
    mesh = _hex_pair()
    core = mio.refine(mesh)
    # Each unit face becomes four quarter faces.
    assert len(_side(core, "floor").entries) == 8
    assert _measure(core, "floor", 2) == pytest.approx(2.0)
    _force_twin(monkeypatch, "refine")
    twin = mio.refine(mesh)
    np.testing.assert_array_equal(
        _side(core, "floor").entries, _side(twin, "floor").entries
    )


def test_decimate_keeps_the_boundary_edges_and_the_twin_agrees(monkeypatch):
    mesh = _tri_strip()
    core = mio.decimate(mesh, ratio=0.5)
    assert _measure(core, "bottom", 1) == pytest.approx(4.0)
    _force_twin(monkeypatch, "decimate")
    twin = mio.decimate(mesh, ratio=0.5)
    np.testing.assert_array_equal(
        _side(core, "bottom").entries, _side(twin, "bottom").entries
    )


def test_subdivide_and_agglomerate_keep_the_floor():
    mesh = _hex_pair()
    assert _measure(mio.subdivide(mesh), "floor", 2) == pytest.approx(2.0)
    assert _measure(
        mio.agglomerate(mesh, target_group_size=2), "floor", 2
    ) == pytest.approx(2.0)


def test_repair_finds_a_flipped_edge_by_its_nodes():
    mesh = _tri_strip()
    data = mesh.cells[0].data.copy()
    data[0] = data[0][::-1]  # (a, b, c) -> (c, b, a): edge (a, b) is now edge 1
    mesh = mio.Mesh(mesh.points, [("triangle", data)])
    mesh.regions = [
        Region("bottom", "side", np.array([[0, 1]] + [[2 * x, 0] for x in (1, 2, 3)]))
    ]
    out = mio.repair(mesh)
    assert len(_side(out, "bottom").entries) == 4
    assert _measure(out, "bottom", 1) == pytest.approx(4.0)


def test_refine_then_undo_green_restores_the_side_region():
    coarse = _tri_strip()
    fine = mio.refine(
        coarse, cells=[0, 2, 4], record_hierarchy=True, record_levels=True
    )
    back = mio.undo_green(coarse, fine)
    # The red-refined cells keep their halves; the green ones merge back.
    assert _measure(back, "bottom", 1) == pytest.approx(4.0)
    assert _measure(fine, "bottom", 1) == pytest.approx(4.0)


@pytest.mark.parametrize("op", ["extract_surface", "extract_skin"])
def test_a_side_region_becomes_a_cell_region_of_the_surface(op, monkeypatch):
    mesh = _hex_pair()
    mesh.regions = list(mesh.regions) + [
        Region("corner", "point", np.array([0])),
        Region("left", "cell", np.array([0])),
    ]
    core = getattr(mio, op)(mesh)
    floor = next(r for r in core.regions if r.name == "floor")
    assert floor.kind == "cell" and len(floor.entries) == 2
    for c in floor.entries.tolist():
        # Both floor cells lie on z = 0.
        pts = np.asarray(core.points)[core.cells[0].data[c]]
        assert np.allclose(pts[:, 2], 0.0)
    assert not any(r.name == "left" for r in core.regions)
    corner = next(r for r in core.regions if r.name == "corner")
    np.testing.assert_allclose(np.asarray(core.points)[corner.entries[0]], [0, 0, 0])
    _force_twin(monkeypatch, op)
    twin = getattr(mio, op)(mesh)
    assert sorted((r.name, r.kind, r.entries.tolist()) for r in core.regions) == sorted(
        (r.name, r.kind, r.entries.tolist()) for r in twin.regions
    )


def test_an_interior_facet_is_lost_with_a_warning(capfd):
    mesh = _hex_pair()
    # Face 2 of hexahedron 0 is x = 1, shared with hexahedron 1: not boundary.
    shared = next(
        f
        for f in range(6)
        if np.allclose(np.asarray(mesh.points)[facet_nodes(mesh, 0, f)[1]][:, 0], 1.0)
    )
    mesh.regions = [Region("wall", "side", np.array([[0, shared]]))]
    out = mio.extract_surface(mesh)
    assert next(r for r in out.regions if r.name == "wall").entries.size == 0
    assert "lost 1" in capfd.readouterr().err

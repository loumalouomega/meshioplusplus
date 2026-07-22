"""Tests for the planar cross-section operation ``meshioplusplus.slice``."""

import numpy as np
import pytest

import meshioplusplus as mp
from meshioplusplus._slice import _slice_py


# --------------------------------------------------------------------------- #
# geometry builders with hand-computable answers                              #
# --------------------------------------------------------------------------- #
def _hex_cube(n=2):
    """The unit cube [0,1]^3 as an n x n x n grid of hexahedra, f = z."""
    xs = np.linspace(0.0, 1.0, n + 1)
    pts = np.array([[x, y, z] for z in xs for y in xs for x in xs], dtype=np.float64)

    def idx(i, j, k):
        return (k * (n + 1) + j) * (n + 1) + i

    cells = [
        [
            idx(i, j, k),
            idx(i + 1, j, k),
            idx(i + 1, j + 1, k),
            idx(i, j + 1, k),
            idx(i, j, k + 1),
            idx(i + 1, j, k + 1),
            idx(i + 1, j + 1, k + 1),
            idx(i, j + 1, k + 1),
        ]
        for k in range(n)
        for j in range(n)
        for i in range(n)
    ]
    m = mp.Mesh(pts, [("hexahedron", np.array(cells))])
    m.point_data["f"] = pts[:, 2].copy()
    m.cell_data["mat"] = [np.arange(len(cells), dtype=np.int32)]
    return m


def _tet_cube(n=2):
    """The same unit cube, but each hexahedron split into 6 tetrahedra."""
    hexes = _hex_cube(n)
    tets = mp.convert_cells(hexes, mode="simplexify")
    return tets


def _quad_sheet():
    """A 2x1 sheet of quads in the xy-plane, g = x."""
    pts = np.array([[0, 0], [1, 0], [1, 1], [0, 1], [2, 0], [2, 1]], dtype=np.float64)
    m = mp.Mesh(pts, [("quad", np.array([[0, 1, 2, 3], [1, 4, 5, 2]]))])
    m.point_data["g"] = pts[:, 0].copy()
    return m


def _section_area(mesh):
    area = 0.0
    for cb in mesh.cells:
        for cell in np.asarray(cb.data):
            p = mesh.points[cell]
            if len(cell) == 3:
                area += 0.5 * np.linalg.norm(np.cross(p[1] - p[0], p[2] - p[0]))
            elif len(cell) == 4:
                area += 0.5 * np.linalg.norm(np.cross(p[2] - p[0], p[3] - p[1]))
    return area


# --------------------------------------------------------------------------- #
# behaviour                                                                    #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("builder", [_hex_cube, _tet_cube])
def test_midplane_area_and_field(builder):
    mesh = builder(2)
    sec = mp.slice(mesh, origin=(0, 0, 0.5), normal=(0, 0, 1))
    # The mid-plane cross-section of the unit cube is a unit square.
    assert _section_area(sec) == pytest.approx(1.0, abs=1e-12)
    # A volume mesh sections into surface faces at z = 0.5.
    assert all(cb.type in ("triangle", "quad") for cb in sec.cells)
    assert np.allclose(sec.points[:, 2], 0.5)
    # The linear field f = z is interpolated exactly to 0.5 at the cut.
    assert np.allclose(sec.point_data["f"], 0.5)


@pytest.mark.parametrize("builder", [_hex_cube, _tet_cube])
def test_section_faces_normals_align_with_plane(builder):
    sec = mp.slice(builder(2), origin=(0, 0, 0.5), normal=(0, 0, 1))
    for cb in sec.cells:
        for cell in np.asarray(cb.data):
            p = sec.points[cell]
            nrm = np.cross(p[1] - p[0], p[2] - p[0])
            # every face is wound with its Newell normal on the +normal side
            assert nrm[2] > 1e-12


def test_missing_plane_is_empty():
    sec = mp.slice(_hex_cube(2), origin=(0, 0, 5.0), normal=(0, 0, 1))
    assert len(sec.cells) == 0
    assert len(sec.points) == 0


def test_2d_mesh_yields_line_segments():
    sec = mp.slice(_quad_sheet(), origin=(0, 0.5, 0), normal=(0, 1, 0))
    assert [cb.type for cb in sec.cells] == ["line"]
    total = 0.0
    for cb in sec.cells:
        for cell in np.asarray(cb.data):
            p = sec.points[cell]
            total += np.linalg.norm(p[1] - p[0])
    # the plane y = 0.5 crosses both unit-width quads -> total length 2.0
    assert total == pytest.approx(2.0, abs=1e-12)
    # g = x is interpolated exactly along the cut
    assert sec.point_data["g"].min() == pytest.approx(0.0)
    assert sec.point_data["g"].max() == pytest.approx(2.0)


def test_grazing_shared_face_not_double_emitted():
    # Two hexes stacked in z sharing the internal face at z = 0.5; the plane
    # lies exactly on that shared face.
    pts = np.array(
        [[x, y, z] for z in (0.0, 0.5, 1.0) for y in (0.0, 1.0) for x in (0.0, 1.0)],
        dtype=np.float64,
    )
    idx = lambda i, j, k: (k * 2 + j) * 2 + i  # noqa: E731
    cells = [
        [
            idx(0, 0, k),
            idx(1, 0, k),
            idx(1, 1, k),
            idx(0, 1, k),
            idx(0, 0, k + 1),
            idx(1, 0, k + 1),
            idx(1, 1, k + 1),
            idx(0, 1, k + 1),
        ]
        for k in range(2)
    ]
    mesh = mp.Mesh(pts, [("hexahedron", np.array(cells))])
    sec = mp.slice(mesh, origin=(0, 0, 0.5), normal=(0, 0, 1))
    # The grazing face is emitted once (area 1.0), not twice (area 2.0).
    assert _section_area(sec) == pytest.approx(1.0, abs=1e-12)


def test_watertight_shared_edges_are_deduped():
    # A clean interior cut (z = 0.4 lies strictly inside the bottom hex layer):
    # a crossing point on an edge shared by two simplices must be a single output
    # node, so no two distinct output nodes share a location -- the section is
    # watertight (grazing configs are the documented exception, tested above).
    sec = mp.slice(_hex_cube(2), origin=(0, 0, 0.4), normal=(0, 0, 1))
    assert len(sec.points) > 0
    assert len(np.unique(sec.points, axis=0)) == len(sec.points)
    # Every connectivity index refers to a real node (no dangling ids).
    for cb in sec.cells:
        assert np.asarray(cb.data).max() < len(sec.points)


def test_parent_cell_provenance():
    # z = 0.25 cleanly cuts the bottom hex layer (cells 0..3 in the 2x2x2 grid).
    mesh = _hex_cube(2)
    sec = mp.slice(mesh, origin=(0, 0, 0.25), normal=(0, 0, 1), record_parent_ids=True)
    assert "slice:parent_cell" in sec.cell_data
    parents = np.concatenate(
        [np.asarray(b) for b in sec.cell_data["slice:parent_cell"]]
    )
    # Only the bottom layer of hexes straddles z = 0.25.
    assert set(parents.tolist()) == {0, 1, 2, 3}
    # cell_data is replicated from the parent cell (same block order as parents).
    mat = np.concatenate([np.asarray(b) for b in sec.cell_data["mat"]])
    assert np.array_equal(mat, parents.astype(mat.dtype))


def test_zero_normal_raises():
    with pytest.raises(ValueError):
        mp.slice(_hex_cube(1), origin=(0, 0, 0), normal=(0, 0, 0))


def test_sets_not_carried():
    mesh = _hex_cube(2)
    mesh.point_sets = {"a": np.array([0, 1, 2])}
    sec = mp.slice(mesh, origin=(0, 0, 0.5), normal=(0, 0, 1))
    assert not getattr(sec, "point_sets", None)


# --------------------------------------------------------------------------- #
# C++ / numpy byte-parity and determinism                                     #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize(
    "origin,normal,rec",
    [
        ((0, 0, 0.5), (0, 0, 1), True),
        ((0, 0, 0.4), (0, 0, 1), False),
        ((0.5, 0.5, 0.5), (1, 1, 1), True),
        ((0, 0, 5.0), (0, 0, 1), True),
    ],
)
def test_cpp_matches_python(origin, normal, rec):
    core = pytest.importorskip("meshioplusplus._core")
    mesh = _hex_cube(3)
    mesh.point_data["v"] = mesh.points.copy()  # a vector field too
    got = core.slice(mesh, list(map(float, origin)), list(map(float, normal)), rec)
    ref = _slice_py(mesh, origin, normal, rec)

    assert got.points.dtype == ref.points.dtype
    assert np.array_equal(got.points, ref.points)
    assert [cb.type for cb in got.cells] == [cb.type for cb in ref.cells]
    for a, b in zip(got.cells, ref.cells):
        da, db = np.asarray(a.data), np.asarray(b.data)
        assert da.dtype == db.dtype and np.array_equal(da, db)
    assert set(got.point_data) == set(ref.point_data)
    for k in got.point_data:
        x, y = np.asarray(got.point_data[k]), np.asarray(ref.point_data[k])
        assert x.dtype == y.dtype and np.array_equal(x, y)
    assert set(got.cell_data) == set(ref.cell_data)
    for k in got.cell_data:
        for xb, yb in zip(got.cell_data[k], ref.cell_data[k]):
            xb, yb = np.asarray(xb), np.asarray(yb)
            assert xb.dtype == yb.dtype and np.array_equal(xb, yb)


def test_determinism_two_runs():
    mesh = _hex_cube(3)
    a = mp.slice(mesh, origin=(0, 0, 0.4), normal=(0, 0, 1), record_parent_ids=True)
    b = mp.slice(mesh, origin=(0, 0, 0.4), normal=(0, 0, 1), record_parent_ids=True)
    assert a.points.tobytes() == b.points.tobytes()
    for x, y in zip(a.cells, b.cells):
        assert np.asarray(x.data).tobytes() == np.asarray(y.data).tobytes()

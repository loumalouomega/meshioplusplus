"""Tests for the isosurface / contour operation ``meshioplusplus.isosurface``."""

import numpy as np
import pytest

import meshioplusplus as mp
from meshioplusplus._isosurface import _isosurface_py


# --------------------------------------------------------------------------- #
# geometry builders with hand-computable answers                              #
# --------------------------------------------------------------------------- #
def _hex_cube(n=2, lo=0.0, hi=1.0):
    """``[lo,hi]^3`` as an n x n x n grid of hexahedra, carrying two fields."""
    xs = np.linspace(lo, hi, n + 1)
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
    m.point_data["fx"] = pts[:, 0].copy()  # a linear field
    m.point_data["r2"] = (pts * pts).sum(axis=1)  # a radial field
    m.point_data["g"] = 2.0 * pts[:, 1] + 1.0  # a second linear field
    m.cell_data["mat"] = [np.arange(len(cells), dtype=np.int32)]
    return m


def _tet_cube(n=2):
    """The same unit cube, but each hexahedron split into 6 tetrahedra."""
    return mp.convert_cells(_hex_cube(n), mode="simplexify")


def _quad_sheet():
    """A 2x1 sheet of quads in the xy-plane, ``h = x``."""
    pts = np.array([[0, 0], [1, 0], [1, 1], [0, 1], [2, 0], [2, 1]], dtype=np.float64)
    m = mp.Mesh(pts, [("quad", np.array([[0, 1, 2, 3], [1, 4, 5, 2]]))])
    m.point_data["h"] = pts[:, 0].copy()
    return m


def _area(mesh):
    area = 0.0
    for cb in mesh.cells:
        for cell in np.asarray(cb.data):
            p = mesh.points[cell]
            if len(cell) == 3:
                area += 0.5 * np.linalg.norm(np.cross(p[1] - p[0], p[2] - p[0]))
            elif len(cell) == 4:
                area += 0.5 * np.linalg.norm(np.cross(p[2] - p[0], p[3] - p[1]))
    return area


def _flat(mesh, name):
    return np.concatenate([np.asarray(b) for b in mesh.cell_data[name]])


# --------------------------------------------------------------------------- #
# behaviour                                                                    #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("builder", [_hex_cube, _tet_cube])
def test_linear_field_level_set_is_planar(builder):
    mesh = builder(2)
    iso = mp.isosurface(mesh, "fx", 0.5)
    # f = x is linear, so its 0.5 level set is exactly the plane x = 0.5 ...
    assert np.allclose(iso.points[:, 0], 0.5, atol=1e-12)
    # ... whose cross-section of the unit cube is a unit square.
    assert _area(iso) == pytest.approx(1.0, abs=1e-12)
    assert all(cb.type in ("triangle", "quad") for cb in iso.cells)


def test_sphere_area_matches_analytic():
    # f = x^2 + y^2 + z^2 on [-1,1]^3: the r^2 level set is a full sphere of
    # radius r, whose area is 4*pi*r^2 up to the discretization error.
    r = 0.6
    iso = mp.isosurface(_hex_cube(12, lo=-1.0, hi=1.0), "r2", r * r)
    assert _area(iso) == pytest.approx(4.0 * np.pi * r * r, rel=0.06)


@pytest.mark.parametrize("builder", [_hex_cube, _tet_cube])
def test_contoured_field_is_exactly_the_isovalue(builder):
    iso = mp.isosurface(builder(2), "fx", 0.5)
    # Not `approx`: the contoured array is written exactly (interpolation only
    # reaches the isovalue to within round-off).
    assert np.all(np.asarray(iso.point_data["fx"]) == 0.5)


def test_other_arrays_are_interpolated_linearly():
    mesh = _hex_cube(2)
    iso = mp.isosurface(mesh, "fx", 0.5)
    # g = 2y + 1 is linear, so interpolation at the cut is exact.
    assert np.allclose(iso.point_data["g"], 2.0 * iso.points[:, 1] + 1.0, atol=1e-12)


def test_magnitude_reduction_is_not_forced_exact():
    mesh = _hex_cube(2)
    mesh.point_data["v"] = np.asarray(mesh.points).copy()
    iso = mp.isosurface(mesh, "v", 1.0)  # |v| = 1, no component given
    assert len(iso.cells) > 0
    got = np.linalg.norm(np.asarray(iso.point_data["v"]), axis=1)
    # |lerp(v)| != lerp(|v|) mathematically, so this is close but not exact.
    assert np.allclose(got, 1.0, atol=5e-2)


def test_component_selection_contours_that_component():
    mesh = _hex_cube(2)
    mesh.point_data["v"] = np.asarray(mesh.points).copy()
    iso = mp.isosurface(mesh, "v", 0.5, component=1)
    assert np.all(np.asarray(iso.point_data["v"])[:, 1] == 0.5)
    assert np.allclose(iso.points[:, 1], 0.5, atol=1e-12)


def test_multiple_isovalues_are_tagged_in_ascending_order():
    iso = mp.isosurface(_hex_cube(2), "fx", [0.75, 0.25])
    values = _flat(iso, "iso:value")
    indices = _flat(iso, "iso:index")
    assert values.dtype == np.float64 and indices.dtype == np.int64
    # Sorted ascending on the way in, so index 0 is the smaller value.
    assert dict(zip(indices.tolist(), values.tolist())) == {0: 0.25, 1: 0.75}
    # Each contour is the plane at its own isovalue.
    for value in (0.25, 0.75):
        pts = iso.points[
            np.unique(
                np.concatenate(
                    [
                        np.asarray(cb.data)[np.asarray(b) == value].ravel()
                        for cb, b in zip(iso.cells, iso.cell_data["iso:value"])
                    ]
                )
            )
        ]
        assert np.allclose(pts[:, 0], value, atol=1e-12)


def test_duplicate_isovalues_are_cut_once():
    once = mp.isosurface(_hex_cube(2), "fx", [0.5])
    twice = mp.isosurface(_hex_cube(2), "fx", [0.5, 0.5])
    assert _area(twice) == pytest.approx(_area(once), abs=1e-12)


def test_iso_index_round_trips_through_split():
    iso = mp.isosurface(_hex_cube(2), "fx", [0.25, 0.75])
    pieces = mp.split(iso, by="region", tag="iso:index")
    assert len(pieces) == 2
    for piece in pieces.values():
        assert len(np.unique(_flat(piece, "iso:value"))) == 1


def test_out_of_range_isovalue_is_empty_not_an_error():
    mesh = _hex_cube(2)
    # One value inside, one outside: only the inside one contributes.
    iso = mp.isosurface(mesh, "fx", [0.5, 9.0])
    assert len(np.unique(_flat(iso, "iso:value"))) == 1
    # All outside -> an empty mesh (still not an error).
    empty = mp.isosurface(mesh, "fx", [-3.0, 9.0])
    assert len(empty.cells) == 0
    assert len(empty.points) == 0


def test_2d_mesh_yields_line_contours():
    iso = mp.isosurface(_quad_sheet(), "h", 0.5)
    assert [cb.type for cb in iso.cells] == ["line"]
    total = 0.0
    for cb in iso.cells:
        for cell in np.asarray(cb.data):
            p = iso.points[cell]
            total += np.linalg.norm(p[1] - p[0])
    # h = x crosses the left quad only, over its full unit height.
    assert total == pytest.approx(1.0, abs=1e-12)
    assert np.all(np.asarray(iso.point_data["h"]) == 0.5)


def test_cell_data_field_raises_by_name():
    mesh = _hex_cube(2)
    with pytest.raises(ValueError) as excinfo:
        mp.isosurface(mesh, "mat", 1.0)
    message = str(excinfo.value)
    assert "mat" in message
    # The error must name the conversion that makes the field contourable.
    assert "to-point" in message or "cell_data_to_point_data" in message


def test_unknown_array_raises_listing_what_exists():
    with pytest.raises(ValueError) as excinfo:
        mp.isosurface(_hex_cube(2), "nope", 0.5)
    assert "fx" in str(excinfo.value)


def test_bad_arguments_raise():
    mesh = _hex_cube(1)
    with pytest.raises(ValueError):
        mp.isosurface(mesh, "fx", [])
    with pytest.raises(ValueError):
        mp.isosurface(mesh, "fx", [float("nan")])
    with pytest.raises(ValueError):
        mp.isosurface(mesh, "fx", 0.5, component=7)


def test_plateau_at_the_isovalue_is_emitted_once():
    # Two hexes stacked in z sharing the internal face at z = 0.5, with a field
    # equal to the isovalue exactly on that shared face. `d >= 0` is the
    # positive side, so the sign mask is total and the face is emitted once.
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
    mesh.point_data["fz"] = pts[:, 2].copy()
    iso = mp.isosurface(mesh, "fz", 0.5)
    assert _area(iso) == pytest.approx(1.0, abs=1e-12)


def test_watertight_shared_edges_are_deduped():
    iso = mp.isosurface(_hex_cube(2), "fx", 0.4)
    assert len(iso.points) > 0
    # A crossing on an edge shared by two simplices is a single output node.
    assert len(np.unique(iso.points, axis=0)) == len(iso.points)
    for cb in iso.cells:
        assert np.asarray(cb.data).max() < len(iso.points)


def test_faces_are_wound_toward_increasing_field():
    mesh = _hex_cube(2)
    iso = mp.isosurface(mesh, "fx", 0.5)
    for cb in iso.cells:
        for cell in np.asarray(cb.data):
            p = iso.points[cell]
            nrm = np.cross(p[1] - p[0], p[2] - p[0])
            # f = x increases along +x, so every normal points that way.
            assert nrm[0] > 1e-12


def test_parent_cell_provenance_and_replicated_cell_data():
    mesh = _hex_cube(2)
    iso = mp.isosurface(mesh, "fx", 0.25, record_parent_ids=True)
    assert "iso:parent_cell" in iso.cell_data
    parents = _flat(iso, "iso:parent_cell")
    # x = 0.25 straddles only the i = 0 column of the 2x2x2 hex grid.
    assert set(parents.tolist()) == {0, 2, 4, 6}
    mat = _flat(iso, "mat")
    assert np.array_equal(mat, parents.astype(mat.dtype))


def test_sets_not_carried():
    mesh = _hex_cube(2)
    mesh.point_sets = {"a": np.array([0, 1, 2])}
    iso = mp.isosurface(mesh, "fx", 0.5)
    assert not getattr(iso, "point_sets", None)


# --------------------------------------------------------------------------- #
# C++ / numpy byte-parity and determinism                                     #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize(
    "array,isovalues,component,rec",
    [
        ("fx", [0.5], None, True),
        ("fx", [0.25, 0.75], None, False),
        ("r2", [1.0], None, True),
        ("v", [0.5], 0, True),  # a vector array, one component
        ("v", [1.0], None, False),  # the same array by magnitude
        ("fx", [9.0], None, True),  # nothing crossed
    ],
)
def test_cpp_matches_python(array, isovalues, component, rec):
    core = pytest.importorskip("meshioplusplus._core")
    mesh = _hex_cube(3)
    mesh.point_data["v"] = np.asarray(mesh.points).copy()
    got = core.isosurface(
        mesh,
        array,
        [float(v) for v in isovalues],
        -1 if component is None else component,
        rec,
    )
    ref = _isosurface_py(mesh, array, isovalues, component, rec)

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
    a = mp.isosurface(mesh, "r2", [0.5, 1.5], record_parent_ids=True)
    b = mp.isosurface(mesh, "r2", [0.5, 1.5], record_parent_ids=True)
    assert a.points.tobytes() == b.points.tobytes()
    for x, y in zip(a.cells, b.cells):
        assert np.asarray(x.data).tobytes() == np.asarray(y.data).tobytes()

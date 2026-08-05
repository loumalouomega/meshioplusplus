"""Regular grids and signed distance (``grid``, ``voxelize``, ``sample_distance``).

The C++/numpy parity tests here assert **raw bytes**, not closeness: the two
implementations are supposed to be the same computation, and a tolerance would
hide exactly the drift the pairing exists to catch.
"""

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import (
    distance_to_surface,
    grid,
    sample_distance,
    surface_watertight_check,
    voxelize,
)


def _prism(footprint, lo, hi):
    """A closed prism over a counter-clockwise xy polygon, wound outward."""
    n = len(footprint)
    pts = [[x, y, lo] for x, y in footprint] + [[x, y, hi] for x, y in footprint]
    tris = []
    for k in range(1, n - 1):
        tris.append([0, k + 1, k])
        tris.append([n, n + k, n + k + 1])
    for i in range(n):
        j = (i + 1) % n
        tris.append([i, j, n + j])
        tris.append([i, n + j, n + i])
    return meshioplusplus.Mesh(
        np.array(pts, dtype=float), [("triangle", np.array(tris, dtype=np.int64))]
    )


def _box(half=1.0):
    return _prism(
        [(-half, -half), (half, -half), (half, half), (-half, half)], -half, half
    )


def _unit_cube():
    return _prism([(0, 0), (1, 0), (1, 1), (0, 1)], 0, 1)


def _l_solid():
    return _prism([(0, 0), (2, 0), (2, 1), (1, 1), (1, 2), (0, 2)], 0, 1)


def _l_contains(p):
    if not (0 < p[0] < 2 and 0 < p[1] < 2 and 0 < p[2] < 1):
        return False
    return not (p[0] > 1 and p[1] > 1)


def _box_sdf(p, half):
    q = np.abs(p) - half
    return np.linalg.norm(np.maximum(q, 0.0)) + min(np.max(q), 0.0)


# --------------------------------------------------------------------------- #
# grid                                                                         #
# --------------------------------------------------------------------------- #
def test_grid_counts_and_type():
    m = grid([3, 4, 5])
    assert len(m.points) == 4 * 5 * 6
    assert len(m.cells) == 1
    assert m.cells[0].type == "hexahedron"
    assert len(m.cells[0].data) == 3 * 4 * 5


def test_grid_cells_are_perfect_hexahedra():
    # The free structural oracle: a right parallelepiped has scaled Jacobian
    # exactly 1, so a transposed axis shows up as an inverted cell.
    m = grid([3, 3, 3], spacing=(2.0, 0.5, 1.5))
    q = meshioplusplus.attach_quality(m)
    sj = np.asarray(q.cell_data["quality:scaled_jacobian"][0]).reshape(-1)
    assert np.all(sj == 1.0)
    st = meshioplusplus.compute_stats(m)
    assert st["num_inverted"] == 0
    assert st["signed_volume"] == pytest.approx(27 * 2.0 * 0.5 * 1.5)


def test_grid_coordinates_are_origin_plus_index_times_spacing():
    m = grid([4, 1, 1], origin=(-1.0, 2.0, 0.5), spacing=(0.25, 3.0, 7.0))
    pts = np.asarray(m.points)
    for i in range(5):
        assert pts[i, 0] == -1.0 + i * 0.25


def test_grid_empty_is_a_mesh_not_a_throw():
    m = grid([0, 0, 0])
    assert len(m.points) == 0
    assert len(m.cells) == 0


def test_grid_bad_arguments_raise_by_name():
    with pytest.raises(ValueError, match="negative"):
        grid([-1, 1, 1])
    with pytest.raises(ValueError, match="spacing"):
        grid([1, 1, 1], spacing=(0.0, 1.0, 1.0))
    with pytest.raises(ValueError, match="max_cells"):
        grid([100, 100, 100], max_cells=1000)


@pytest.mark.parametrize(
    "dims,origin,spacing",
    [
        ([2, 2, 2], (0.0, 0.0, 0.0), (1.0, 1.0, 1.0)),
        ([3, 1, 4], (0.5, -1.0, 2.0), (0.25, 3.0, 0.5)),
        ([1, 1, 1], (0.0, 0.0, 0.0), (1.0, 1.0, 1.0)),
    ],
)
def test_grid_cpp_matches_python(dims, origin, spacing):
    core = pytest.importorskip("meshioplusplus._core")
    from meshioplusplus._grid import _grid_py

    got = core.grid(dims, list(origin), list(spacing), 20000000)
    want = _grid_py(np.array(dims), np.array(origin), np.array(spacing))
    assert np.asarray(got.points).tobytes() == np.asarray(want.points).tobytes()
    assert (
        np.asarray(got.cells[0].data).tobytes()
        == np.asarray(want.cells[0].data).tobytes()
    )


# --------------------------------------------------------------------------- #
# voxelize                                                                     #
# --------------------------------------------------------------------------- #
def test_voxelize_all_covers_the_bounding_box():
    out, report = voxelize(_unit_cube(), resolution=[4, 4, 4], return_report=True)
    assert len(out.cells[0].data) == 64
    assert report["num_occupied"] == 64
    assert report["dims"] == [4, 4, 4]
    assert report["spacing"] == pytest.approx([0.25, 0.25, 0.25])
    # The whole point of the design: the result is an ordinary mesh.
    assert meshioplusplus.compute_stats(out)["signed_volume"] == pytest.approx(1.0)


def test_voxelize_inside_keeps_the_interior():
    out, report = voxelize(
        _unit_cube(),
        resolution=[5, 5, 5],
        bounds=[-0.5, -0.5, -0.5, 1.5, 1.5, 1.5],
        fill="inside",
        return_report=True,
    )
    # Cell size 0.4; centres at -0.3, 0.1, 0.5, 0.9, 1.3 -- three per axis inside.
    assert report["num_occupied"] == 27
    assert len(out.cells[0].data) == 27


def test_voxelize_surface_is_a_shell():
    out, report = voxelize(
        _unit_cube(), resolution=[6, 6, 6], fill="surface", return_report=True
    )
    assert 0 < report["num_occupied"] < 216
    assert len(out.cells[0].data) == report["num_occupied"]


def test_voxelize_requires_exactly_one_of_resolution_and_cell_size():
    with pytest.raises(ValueError, match="exactly one"):
        voxelize(_unit_cube())
    with pytest.raises(ValueError, match="exactly one"):
        voxelize(_unit_cube(), resolution=[2, 2, 2], cell_size=0.5)


def test_voxelize_cell_size_covers_the_box():
    _out, report = voxelize(_unit_cube(), cell_size=0.3, return_report=True)
    assert report["dims"][0] == 4  # ceil(1.0 / 0.3)
    assert report["spacing"][0] == pytest.approx(0.3)


def test_voxelize_budget_refuses_by_name():
    with pytest.raises(ValueError, match="max_cells"):
        voxelize(_unit_cube(), resolution=[100, 100, 100], max_cells=1000)


def test_voxelize_unknown_fill_raises_by_name():
    with pytest.raises(ValueError, match="unknown fill"):
        voxelize(_unit_cube(), resolution=[2, 2, 2], fill="solid")


def test_voxelize_occupancy_is_attached_on_request():
    out = voxelize(_unit_cube(), resolution=[3, 3, 3], attach_occupancy=True)
    assert "voxel:occupancy" in out.cell_data
    assert np.all(np.asarray(out.cell_data["voxel:occupancy"][0]) == 1)


@pytest.mark.parametrize(
    "fill,kwargs",
    [
        ("all", dict(resolution=[4, 4, 4])),
        ("surface", dict(resolution=[6, 6, 6])),
        (
            "inside",
            dict(resolution=[5, 5, 5], bounds=[-0.5, -0.5, -0.5, 1.5, 1.5, 1.5]),
        ),
    ],
)
def test_cpp_matches_python(fill, kwargs):
    """Byte-identity across the C++-core/numpy-fallback boundary.

    This pins **parity only** -- both sides can be identically wrong -- so it never
    substitutes for the analytic tests above and in tests/cpp/.
    """
    core = pytest.importorskip("meshioplusplus._core")
    from meshioplusplus._voxelize import _voxelize_py

    mesh = _unit_cube()
    got = core.voxelize(
        mesh,
        kwargs.get("resolution"),
        None,
        kwargs.get("bounds"),
        0.0,
        0.0,
        fill,
        False,
        20000000,
        "pseudonormal",
        "off",
    )
    want, report = _voxelize_py(
        mesh,
        kwargs.get("resolution"),
        None,
        kwargs.get("bounds"),
        0.0,
        0.0,
        fill,
        False,
        20000000,
        "pseudonormal",
    )
    assert got["num_occupied"] == report["num_occupied"]
    assert np.asarray(got["mesh"].points).tobytes() == np.asarray(want.points).tobytes()
    a = np.asarray(got["mesh"].cells[0].data)
    b = np.asarray(want.cells[0].data)
    assert a.dtype == b.dtype and a.shape == b.shape
    assert a.tobytes() == b.tobytes()


# --------------------------------------------------------------------------- #
# distance                                                                     #
# --------------------------------------------------------------------------- #
def test_box_distance_is_exact():
    # A box's surface is represented exactly by triangles, so there is no
    # tessellation error to hide behind.
    box = _box(1.0)
    rng = np.random.default_rng(0)
    q = rng.uniform(-2.5, 2.5, size=(300, 3))
    got = np.asarray(sample_distance(box, q))
    want = np.array([_box_sdf(p, 1.0) for p in q])
    assert np.allclose(got, want, atol=1e-12)


def test_sign_is_correct_throughout_the_l_solid():
    solid = _l_solid()
    rng = np.random.default_rng(1)
    q = rng.uniform(-0.5, 2.5, size=(400, 3))
    got = np.asarray(sample_distance(solid, q))
    for p, d in zip(q, got):
        assert (d < 0.0) == _l_contains(p), f"sign wrong at {p}"


def test_unsigned_is_never_negative():
    got = np.asarray(
        sample_distance(
            _l_solid(),
            np.random.default_rng(2).uniform(-1, 3, size=(100, 3)),
            sign="unsigned",
        )
    )
    assert np.all(got >= 0.0)


def test_the_bucket_size_does_not_change_the_answer():
    """The reason every candidate comparison is totally ordered."""
    box = _box(1.0)
    q = np.random.default_rng(3).uniform(-2.0, 2.0, size=(200, 3))
    reference = np.asarray(sample_distance(box, q))
    for cell in (0.1, 0.37, 1.0, 4.0):
        got = np.asarray(sample_distance(box, q, grid_cell_size=cell))
        assert (
            got.tobytes() == reference.tobytes()
        ), f"cell size {cell} changed the answer"


def test_watertight_check_counts_the_defects():
    assert surface_watertight_check(_box(1.0))["watertight"]
    sheet = meshioplusplus.Mesh(
        np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0]]),
        [("triangle", np.array([[0, 1, 2]], dtype=np.int64))],
    )
    q = surface_watertight_check(sheet)
    assert not q["watertight"]
    assert q["boundary_edges"] == 3


def test_distance_to_surface_attaches_point_data():
    box = _box(1.0)
    query = grid([4, 4, 4], origin=(-2.0, -2.0, -2.0), spacing=(1.0, 1.0, 1.0))
    out, report = distance_to_surface(
        query, box, return_report=True, record_inside=True
    )
    assert "sdf:distance" in out.point_data
    assert "sdf:inside" in out.point_data
    assert len(out.point_data["sdf:distance"]) == len(query.points)
    assert report["quality"]["watertight"]
    # The free cross-subsystem check: the zero level set of the field is the box.
    contour = meshioplusplus.isosurface(out, "sdf:distance", [0.0])
    assert sum(len(cb.data) for cb in contour.cells) > 0


def test_distance_to_surface_cell_location():
    box = _box(1.0)
    query = grid([3, 3, 3], origin=(-1.5, -1.5, -1.5), spacing=(1.0, 1.0, 1.0))
    out = distance_to_surface(query, box, location="center")
    assert "sdf:distance" in out.cell_data
    assert len(out.cell_data["sdf:distance"][0]) == 27


def test_a_banded_run_agrees_exactly_inside_the_band():
    box = _box(1.0)
    query = grid([6, 6, 6], origin=(-2.0, -2.0, -2.0), spacing=(4.0 / 6.0,) * 3)
    full = distance_to_surface(query, box)
    banded = distance_to_surface(query, box, band=0.75)
    flag = np.asarray(banded.point_data["sdf:band"]).reshape(-1)
    a = np.asarray(full.point_data["sdf:distance"])
    b = np.asarray(banded.point_data["sdf:distance"])
    assert np.any(flag == 1) and np.any(flag == 0)
    assert a[flag == 1].tobytes() == b[flag == 1].tobytes()


def test_a_volume_mesh_is_refused_by_name():
    tet = meshioplusplus.Mesh(
        np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]),
        [("tetra", np.array([[0, 1, 2, 3]], dtype=np.int64))],
    )
    with pytest.raises(ValueError, match="extract_surface"):
        surface_watertight_check(tet)


def test_bad_option_names_raise_by_name():
    box = _box(1.0)
    q = np.zeros((1, 3))
    with pytest.raises(ValueError, match="unknown sign"):
        sample_distance(box, q, sign="magic")
    with pytest.raises(ValueError, match="unknown weight"):
        sample_distance(box, q, weight="magic")


@pytest.mark.parametrize("weight", ["angle", "area"])
@pytest.mark.parametrize("surface", ["box", "l"])
def test_distance_cpp_matches_python(weight, surface):
    """Byte-identity of the distance field across the C++/numpy boundary.

    The numpy reference deliberately has no accelerator: the bucket grid is
    provably unable to change the answer (see
    ``test_the_bucket_size_does_not_change_the_answer``), so a brute-force scan in
    ascending triangle order is the same computation.
    """
    core = pytest.importorskip("meshioplusplus._core")
    from meshioplusplus._sdf import _sample_py

    mesh = _box(1.0) if surface == "box" else _l_solid()
    q = np.random.default_rng(4).uniform(-2.5, 2.5, size=(200, 3))
    got = np.asarray(
        core.sample_distance(
            mesh, q, "pseudonormal", weight, 0.0, "off", "", 0.0, 2.0e9
        )
    )
    want, _cells, _band = _sample_py(mesh, q, "pseudonormal", weight, 0.0)
    assert got.dtype == want.dtype and got.shape == want.shape
    assert got.tobytes() == want.tobytes()


def test_the_numpy_reference_refuses_the_winding_number():
    """It sums one atan2 per triangle and thresholds the total; atan2 is not
    correctly rounded, so the two implementations could genuinely disagree.
    Refusing is honest where a tolerance would not be."""
    from meshioplusplus._sdf import _sample_py

    with pytest.raises(NotImplementedError, match="winding-number"):
        _sample_py(_box(1.0), np.zeros((1, 3)), "winding-number", "angle", 0.0)

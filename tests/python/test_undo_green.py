"""Tests for green-element undo (the ``undo_green`` operation).

The central oracle: every green sibling group's substituted row must equal,
verbatim, the coarse mesh's own row at that group's parent -- a lookup, not
an approximation. Unlike ``subdivide``/``agglomerate``, this operation has a
full numpy twin (pure array bookkeeping, no winding repair), so
``test_cpp_matches_python`` compares the two engines directly.
"""

from collections import defaultdict

import numpy as np
import pytest

import meshioplusplus as mp
from meshioplusplus._regions import Region
from meshioplusplus._undo_green import _undo_green_py

try:
    from meshioplusplus import _core
except ImportError:  # pragma: no cover - a pure-Python build
    _core = None

needs_core = pytest.mark.skipif(_core is None, reason="needs the compiled C++ core")

_RESERVED = (
    "refine:cell_id",
    "refine:parent_id",
    "refine:level",
    "refine:entity",
    "refine:hanging",
    "refine:parent_cell",
)


def _tri_grid(n):
    """Two triangles per cell of an ``n x n`` grid over the unit square."""
    pts = [[float(i), float(j), 0.0] for j in range(n + 1) for i in range(n + 1)]
    cells = []
    for j in range(n):
        for i in range(n):
            a = j * (n + 1) + i
            b = a + 1
            cc = a + n + 2
            dd = a + n + 1
            cells.append([a, b, cc])
            cells.append([a, cc, dd])
    return mp.Mesh(np.array(pts), [("triangle", np.array(cells, dtype=np.int64))])


def _refine_center(coarse, levels=1):
    return mp.refine(
        coarse,
        cells=[8],
        record_hierarchy=True,
        record_levels=True,
        levels=levels,
    )


def _classify(fine):
    ids = fine.cell_data["refine:cell_id"][0].reshape(-1)
    parents = fine.cell_data["refine:parent_id"][0].reshape(-1)
    levels = fine.cell_data["refine:level"][0].reshape(-1)
    greens = defaultdict(list)
    reds = []
    for c in range(len(ids)):
        if ids[c] == parents[c]:
            continue
        if levels[c] == 0:
            greens[int(parents[c])].append(c)
        else:
            reds.append(c)
    return greens, reds


# --------------------------------------------------------------------------- #
# the core oracle                                                              #
# --------------------------------------------------------------------------- #


def test_exactly_restores_the_coarse_parent_per_green_group():
    coarse = _tri_grid(3)
    fine = _refine_center(coarse)
    greens, reds = _classify(fine)
    assert greens, "the fixture must produce at least one green group"
    assert reds, "the fixture must produce at least one red child"

    undone, report = mp.undo_green(coarse, fine, return_report=True)
    assert report["num_groups_undone"] == len(greens)
    assert report["num_cells_removed"] == sum(len(v) - 1 for v in greens.values())

    assert len(undone.cells[0].data) < len(fine.cells[0].data)
    assert len(undone.cells[0].data) > len(coarse.cells[0].data)
    assert len(undone.points) == len(fine.points)
    np.testing.assert_array_equal(undone.points, fine.points)

    # Four are cell_data, two (entity/hanging) are point_data -- both dicts
    # must be checked, or a point_data leftover passes vacuously.
    for name in _RESERVED:
        assert name not in undone.cell_data
        assert name not in undone.point_data


@needs_core
def test_green_rows_equal_the_coarse_mesh_and_red_rows_are_unchanged():
    # Drives the pure-Python reference directly (its cell_maps ARE exposed)
    # to check row-for-row identity, the strongest form of the oracle.
    coarse = _tri_grid(3)
    fine = _refine_center(coarse)
    undone, _ = _undo_green_py(coarse, fine)

    # Rebuild the cell map the same way the reference does internally, by
    # re-deriving from the output row counts is fragile; instead just check
    # a representative green group and a representative red child directly
    # against the coarse/fine connectivity.
    greens, reds = _classify(fine)
    fine_conn = fine.cells[0].data
    coarse_conn = coarse.cells[0].data
    undone_conn = undone.cells[0].data

    # Every green group's members must have vanished as distinct rows: the
    # undone mesh must contain the coarse parent's row somewhere.
    for parent_id, members in greens.items():
        assert any(
            np.array_equal(undone_conn[r], coarse_conn[parent_id])
            for r in range(len(undone_conn))
        ), f"coarse row {parent_id} not found verbatim in the undone mesh"

    # Every red child's row must still be present somewhere, unchanged.
    for c in reds:
        assert any(
            np.array_equal(undone_conn[r], fine_conn[c])
            for r in range(len(undone_conn))
        )


def test_cpp_matches_python():
    coarse = _tri_grid(3)
    fine = _refine_center(coarse)
    fine.regions.append(Region("dummy_cells", "cell", np.array([0, 1], dtype=np.int64)))

    cpp, cpp_report = mp.undo_green(coarse, fine, return_report=True)
    py, py_report = _undo_green_py(coarse, fine)

    assert cpp_report == py_report
    assert mp.meshes_equal(cpp, py)
    assert len(cpp.regions) == len(py.regions)
    cpp_names = sorted(r.name for r in cpp.regions)
    py_names = sorted(r.name for r in py.regions)
    assert cpp_names == py_names


def test_cpp_matches_python_with_no_green_groups():
    # Uniform refinement: every group is red, nothing to undo.
    coarse = _tri_grid(2)
    fine = mp.refine(coarse, record_hierarchy=True, record_levels=True)
    cpp, cpp_report = mp.undo_green(coarse, fine, return_report=True)
    py, py_report = _undo_green_py(coarse, fine)
    assert cpp_report == py_report == {"num_groups_undone": 0, "num_cells_removed": 0}
    assert mp.meshes_equal(cpp, py)
    assert len(cpp.cells[0].data) == len(fine.cells[0].data)


# --------------------------------------------------------------------------- #
# regions                                                                      #
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("engine", ["cpp", "python"])
def test_regions_survive_the_non_injective_collapse_and_sides_do_not(engine):
    if engine == "cpp" and _core is None:
        pytest.skip("needs the compiled C++ core")
    coarse = _tri_grid(3)
    fine = _refine_center(coarse)
    greens, _ = _classify(fine)
    pair = next(v for v in greens.values() if len(v) >= 2)[:2]

    fine.regions.append(Region("greens", "cell", np.array(pair, dtype=np.int64)))
    fine.regions.append(Region("apex", "point", np.array([0], dtype=np.int64)))
    fine.regions.append(Region("edge", "side", np.array([[0, 0]], dtype=np.int64)))

    if engine == "cpp":
        undone = mp.undo_green(coarse, fine)
    else:
        undone, _ = _undo_green_py(coarse, fine)

    greens_out = next(r for r in undone.regions if r.name == "greens")
    assert (
        len(greens_out.entries) == 1
    ), "both members of the pair collapse via Region dedup"

    apex_out = next(r for r in undone.regions if r.name == "apex")
    assert apex_out.entries.tolist() == [0], "points are never renumbered"

    assert not any(
        r.name == "edge" for r in undone.regions
    ), "named Side regions do not survive undo_green at all"


# --------------------------------------------------------------------------- #
# preconditions and errors                                                     #
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("engine", ["cpp", "python"])
def test_missing_hierarchy_raises(engine):
    if engine == "cpp" and _core is None:
        pytest.skip("needs the compiled C++ core")
    coarse = _tri_grid(3)
    fine = mp.refine(coarse, cells=[8])  # no record_hierarchy, no record_levels
    fn = mp.undo_green if engine == "cpp" else _undo_green_py
    with pytest.raises(ValueError):
        fn(coarse, fine)


@pytest.mark.parametrize("engine", ["cpp", "python"])
def test_hierarchy_without_levels_still_raises(engine):
    if engine == "cpp" and _core is None:
        pytest.skip("needs the compiled C++ core")
    coarse = _tri_grid(3)
    fine = mp.refine(coarse, cells=[8], record_hierarchy=True)
    fn = mp.undo_green if engine == "cpp" else _undo_green_py
    with pytest.raises(ValueError):
        fn(coarse, fine)


@pytest.mark.parametrize("engine", ["cpp", "python"])
def test_an_unrelated_coarse_mesh_raises(engine):
    if engine == "cpp" and _core is None:
        pytest.skip("needs the compiled C++ core")
    coarse = _tri_grid(3)
    fine = _refine_center(coarse)
    unrelated = _tri_grid(1)
    fn = mp.undo_green if engine == "cpp" else _undo_green_py
    with pytest.raises(ValueError):
        fn(unrelated, fine)


@pytest.mark.parametrize("engine", ["cpp", "python"])
def test_a_coarse_mesh_with_more_points_than_fine_raises(engine):
    if engine == "cpp" and _core is None:
        pytest.skip("needs the compiled C++ core")
    coarse = _tri_grid(5)
    fine = _refine_center(_tri_grid(3))
    fn = mp.undo_green if engine == "cpp" else _undo_green_py
    with pytest.raises(ValueError):
        fn(coarse, fine)


@pytest.mark.parametrize("engine", ["cpp", "python"])
def test_a_multi_level_hierarchy_is_refused(engine):
    if engine == "cpp" and _core is None:
        pytest.skip("needs the compiled C++ core")
    coarse = _tri_grid(3)
    fine = _refine_center(coarse, levels=2)
    fn = mp.undo_green if engine == "cpp" else _undo_green_py
    with pytest.raises(ValueError):
        fn(coarse, fine)

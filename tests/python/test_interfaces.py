"""Tests for the roadmap §5.2 interface/contact operations."""

import numpy as np
import pytest

import meshioplusplus as mio
from meshioplusplus import _interfaces as interfaces
from meshioplusplus._regions import Region

try:
    from meshioplusplus import _core
except ImportError:  # pragma: no cover - pure-Python build
    _core = None

needs_core = pytest.mark.skipif(_core is None, reason="needs the compiled core")


def two_tets():
    return mio.Mesh(
        np.array(
            [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 0, -1]],
            dtype=float,
        ),
        [("tetra", np.array([[0, 1, 2, 3], [0, 2, 1, 4]], dtype=np.int64))],
        regions=[
            Region("upper", "cell", [0], dim=3),
            Region("lower", "cell", [1], dim=3),
            Region("points", "point", [0, 1]),
        ],
        point_data={"temperature": np.arange(5.0)},
        field_data={"case": np.array([7], dtype=np.int64)},
    )


def test_region_adjacency_returns_shared_facet_and_provenance():
    out = mio.region_adjacency(two_tets())
    assert [(block.type, len(block.data)) for block in out.cells] == [("triangle", 1)]
    assert out.cell_data["interface:region_a"][0].tolist() == [0]
    assert out.cell_data["interface:region_b"][0].tolist() == [1]
    # Default selection uses the canonical cell-region ordering: lower, upper.
    assert out.cell_data["interface:parent_cell_a"][0].tolist() == [1]
    assert out.cell_data["interface:parent_cell_b"][0].tolist() == [0]
    assert out.cell_data["interface:shared_count"][0].tolist() == [2]
    assert out.cell_data["interface:measure"][0].tolist() == pytest.approx([0.5])
    assert [(r.name, r.kind, r.entries.tolist()) for r in out.regions] == [
        ("points", "point", [0, 1]),
        ("adjacency:0:1", "cell", [0]),
    ]
    source = two_tets()
    np.testing.assert_array_equal(out.points, source.points)
    np.testing.assert_array_equal(out.point_data["temperature"], np.arange(5.0))
    np.testing.assert_array_equal(out.field_data["case"], [7])


def test_region_adjacency_matches_the_numpy_reference():
    mesh = two_tets()
    native = mio.region_adjacency(mesh)
    reference = interfaces._region_adjacency_py(mesh)
    assert [(b.type, b.data.tolist()) for b in native.cells] == [
        (b.type, b.data.tolist()) for b in reference.cells
    ]
    assert native.cell_data.keys() == reference.cell_data.keys()
    for name in native.cell_data:
        for actual, expected in zip(native.cell_data[name], reference.cell_data[name]):
            np.testing.assert_array_equal(actual, expected)


def test_region_adjacency_selector_order_and_2d_measure():
    mesh = mio.Mesh(
        np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0], [0, 1, 0], [1, 1, 0], [2, 1, 0]]),
        [("quad", np.array([[0, 1, 4, 3], [1, 2, 5, 4]], dtype=np.int64))],
        regions=[
            Region("left", "cell", [0], dim=2),
            Region("right", "cell", [1], dim=2),
        ],
    )
    out = mio.region_adjacency(mesh, ["right", "left"])
    assert out.cells[0].type == "line"
    assert out.cells[0].data.tolist() == [[4, 1]]
    assert out.cell_data["interface:region_a"][0].tolist() == [0]
    assert out.cell_data["interface:region_b"][0].tolist() == [1]
    assert out.cell_data["interface:measure"][0].tolist() == pytest.approx([1.0])


def test_region_adjacency_uses_blocks_if_cell_regions_are_absent():
    mesh = mio.Mesh(
        np.array(
            [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 0, -1]],
            dtype=float,
        ),
        [
            ("tetra", np.array([[0, 1, 2, 3]], dtype=np.int64)),
            ("tetra", np.array([[0, 2, 1, 4]], dtype=np.int64)),
        ],
    )
    out = mio.region_adjacency(mesh)
    assert [(block.type, len(block.data)) for block in out.cells] == [("triangle", 1)]
    assert out.regions[-1].name == "adjacency:0:1"


@pytest.mark.parametrize(
    "regions",
    [["upper"], ["upper", "upper"], ["missing", "lower"]],
)
def test_region_adjacency_rejects_invalid_selections(regions):
    with pytest.raises(ValueError):
        mio.region_adjacency(two_tets(), regions)


def test_empty_selector_list_is_the_default_selection():
    assert len(mio.region_adjacency(two_tets(), []).cells[0].data) == 1


def test_region_adjacency_can_disambiguate_same_name_regions():
    mesh = two_tets()
    mesh.regions.append(Region("upper", "cell", [0], dim=3, tag=9))
    with pytest.raises(ValueError, match="ambiguous"):
        mio.region_adjacency(mesh, ["upper", "lower"])
    out = mio.region_adjacency(
        mesh,
        [{"name": "upper", "dim": 3, "tag": 9}, {"name": "lower", "dim": 3}],
    )
    assert len(out.cells[0].data) == 1


@needs_core
def test_region_adjacency_python_twin_matches_core():
    mesh = two_tets()
    core = mio.region_adjacency(mesh)
    twin = interfaces._region_adjacency_py(mesh)
    assert [(b.type, b.data.tolist()) for b in core.cells] == [
        (b.type, b.data.tolist()) for b in twin.cells
    ]
    for name in core.cell_data:
        for a, b in zip(core.cell_data[name], twin.cell_data[name]):
            np.testing.assert_array_equal(a, b)


def test_find_interface_returns_facet_mesh_side_regions_and_report():
    mesh = two_tets()
    out, report = mio.find_interface(mesh, "upper", "lower", return_report=True)
    assert [(block.type, block.data.tolist()) for block in out.cells] == [
        ("triangle", [[0, 2, 1]])
    ]
    assert report["num_pairs"] == 1
    assert report["area"] == pytest.approx(0.5)
    assert report["max_gap"] == 0.0
    assert report["unmatched_a"] == report["unmatched_b"] == 3
    assert report["side_a"].kind == report["side_b"].kind == "side"
    assert report["side_a"].entries.tolist() == [[0, 3]]
    assert report["side_b"].entries.tolist() == [[1, 3]]
    assert out.cell_data["interface:partner_cell"][0].tolist() == [1]


def test_find_interface_proximity_and_python_twin_match():
    a = mio.Mesh(
        np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1.0]]),
        [("tetra", np.array([[0, 1, 2, 3]], dtype=np.int64))],
        regions=[Region("part", "cell", [0], dim=3)],
    )
    b = mio.Mesh(
        np.array([[0, 0, 0], [0, 1, 0], [1, 0, 0], [0, 0, -1.0]]),
        [("tetra", np.array([[0, 1, 2, 3]], dtype=np.int64))],
        regions=[Region("part", "cell", [0], dim=3)],
    )
    kwargs = dict(mesh_b=b, mode="proximity", gap_tolerance=1e-10, angle_tolerance=5)
    native, report = mio.find_interface(a, "part", "part", return_report=True, **kwargs)
    twin = interfaces._find_interface_py(a, "part", "part", **kwargs)
    assert report["num_pairs"] == 1
    assert report["area"] == pytest.approx(0.5)
    assert [(block.type, block.data.tolist()) for block in native.cells] == [
        (block.type, block.data.tolist()) for block in twin["mesh"].cells
    ]


def stacked_grids(n_lower, n_upper):
    """Two separately meshed blocks that touch on the plane z = 1."""
    lower = mio.grid([n_lower] * 2 + [2], spacing=(1 / n_lower, 1 / n_lower, 0.5))
    upper = mio.grid(
        [n_upper] * 2 + [2],
        origin=(0.0, 0.0, 1.0),
        spacing=(1 / n_upper, 1 / n_upper, 0.5),
    )
    both = mio.merge([lower, upper], source_tag=False)
    n_lo = len(lower.cells[0].data)
    both.regions = [
        Region("lower", "cell", np.arange(n_lo), dim=3),
        Region(
            "upper", "cell", np.arange(n_lo, n_lo + len(upper.cells[0].data)), dim=3
        ),
    ]
    return both


def test_find_interface_between_separately_meshed_parts():
    # The roadmap probe: both modes recover the shared plane, whose area is exact.
    nonmatching = stacked_grids(4, 6)
    conforming = mio.find_interface(nonmatching, "lower", "upper", return_report=True)[
        1
    ]
    assert conforming["num_pairs"] == 0

    kwargs = dict(mode="proximity", gap_tolerance=1e-9, angle_tolerance=5.0)
    interface, report = mio.find_interface(
        nonmatching, "lower", "upper", return_report=True, **kwargs
    )
    assert report["num_pairs"] == 16
    assert report["area"] == pytest.approx(1.0)
    assert report["max_gap"] == pytest.approx(0.0, abs=1e-9)
    twin = interfaces._find_interface_py(nonmatching, "lower", "upper", **kwargs)
    assert twin["report"]["num_pairs"] == 16
    assert [b.data.tolist() for b in interface.cells] == [
        b.data.tolist() for b in twin["mesh"].cells
    ]

    # With the finer part as master the same plane is recovered facet by facet.
    finer = mio.find_interface(
        nonmatching, "lower", "upper", master="b", return_report=True, **kwargs
    )[1]
    assert finer["num_pairs"] == 36
    assert finer["area"] == pytest.approx(1.0)

    # Meshed alike and welded, the conforming mode finds the same plane exactly.
    welded = mio.merge(
        [
            mio.grid([4, 4, 2], spacing=(0.25, 0.25, 0.5)),
            mio.grid([4, 4, 2], origin=(0, 0, 1.0), spacing=(0.25, 0.25, 0.5)),
        ],
        weld=True,
        source_tag=False,
    )
    welded.regions = [
        Region("lower", "cell", np.arange(32), dim=3),
        Region("upper", "cell", np.arange(32, 64), dim=3),
    ]
    same = mio.find_interface(welded, "lower", "upper", return_report=True)[1]
    assert same["num_pairs"] == 16
    assert same["area"] == pytest.approx(1.0)


def test_find_interface_gap_sign_follows_master_side():
    a = mio.Mesh(
        np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1.0]]),
        [("tetra", np.array([[0, 1, 2, 3]], dtype=np.int64))],
        regions=[Region("part", "cell", [0], dim=3)],
    )
    b = mio.Mesh(
        np.array([[0, 0, -0.02], [0, 1, -0.02], [1, 0, -0.02], [0, 0, -1.0]]),
        [("tetra", np.array([[0, 1, 2, 3]], dtype=np.int64))],
        regions=[Region("part", "cell", [0], dim=3)],
    )
    kwargs = dict(mesh_b=b, mode="proximity", gap_tolerance=0.03, angle_tolerance=5)
    master_a = mio.find_interface(a, "part", "part", master="a", **kwargs)
    master_b = mio.find_interface(a, "part", "part", master="b", **kwargs)
    gap_a = master_a.cell_data["interface:gap"][0][0]
    gap_b = master_b.cell_data["interface:gap"][0][0]
    assert abs(gap_a) == pytest.approx(0.02)
    assert abs(gap_b) == pytest.approx(0.02)
    # Measured along each master's own outward normal, a separation is negative.
    assert gap_a < 0 and gap_b < 0


def test_conforming_two_mesh_mode_does_not_match_equal_ids_at_different_positions():
    a = mio.Mesh(
        np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1.0]]),
        [("tetra", np.array([[0, 1, 2, 3]], dtype=np.int64))],
        regions=[Region("part", "cell", [0])],
    )
    b = mio.Mesh(
        np.array([[0, 0, 2], [1, 0, 2], [0, 1, 2], [0, 0, 3.0]]),
        [("tetra", np.array([[0, 1, 2, 3]], dtype=np.int64))],
        regions=[Region("part", "cell", [0])],
    )
    mesh, report = mio.find_interface(a, "part", "part", mesh_b=b, return_report=True)
    assert report["num_pairs"] == 0
    assert not mesh.cells


def test_contact_pairs_returns_projection_and_matches_numpy_twin():
    mesh = two_tets()
    native = mio.contact_pairs(mesh, "points", "lower")
    twin = interfaces._contact_pairs_py(mesh, "points", mesh, "lower")
    assert native["slave_point"].tolist() == [0, 1]
    assert native["master_cell"].tolist() == [1, 1]
    assert native["unmatched"].size == 0
    assert native["local_coordinates"].shape == (2, 3)
    assert native["gap"].tolist() == pytest.approx([0.0, 0.0])
    for key in native:
        np.testing.assert_allclose(native[key], twin[key])


def test_contact_pairs_reports_polygon_fan_subfacet_for_local_coordinates():
    mesh = mio.Mesh(
        np.array(
            [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0.25, 0.75, 0.1]],
            dtype=np.float64,
        ),
        [("quad", np.array([[0, 1, 2, 3]], dtype=np.int64))],
        regions=[Region("master", "cell", [0]), Region("slave", "point", [4])],
    )
    native = mio.contact_pairs(mesh, "slave", "master", tolerance=0.2)
    twin = interfaces._contact_pairs_py(mesh, "slave", mesh, "master", tolerance=0.2)
    assert native["master_subfacet"].tolist() == [1]
    assert native["local_coordinates"].tolist()[0] == pytest.approx([0.25, 0.25, 0.5])
    for key in native:
        np.testing.assert_allclose(native[key], twin[key])


def test_split_interface_carries_point_data_regions_and_cohesive_cells():
    mesh = two_tets()
    side = Region("cut", "side", [[0, 3]])
    native, report = mio.split_interface(
        mesh, side, add_cohesive=True, return_report=True
    )
    twin, twin_report = interfaces._split_interface_py(mesh, side, add_cohesive=True)
    assert (
        report == twin_report == {"num_duplicated_points": 3, "num_cohesive_cells": 1}
    )
    assert native.points.shape[0] == 8
    assert native.cells[-1].type == "wedge"
    assert native.point_data["temperature"].tolist() == [0, 1, 2, 3, 4, 0, 1, 2]
    assert [r.entries.tolist() for r in native.regions if r.name == "points"] == [
        [0, 1, 5, 6]
    ]
    assert [(block.type, block.data.tolist()) for block in native.cells] == [
        (block.type, block.data.tolist()) for block in twin.cells
    ]
    np.testing.assert_array_equal(
        native.point_data["temperature"], twin.point_data["temperature"]
    )


def test_split_interface_deduplicates_cohesive_cells_when_both_sides_are_selected():
    mesh = two_tets()
    side = Region("cut", "side", [[0, 3], [1, 3]])
    out, report = mio.split_interface(mesh, side, add_cohesive=True, return_report=True)
    assert report["num_cohesive_cells"] == 1
    assert out.cells[-1].type == "wedge"


def test_split_interface_2d_cohesive_line_uses_trace_b_cell_data():
    mesh = mio.Mesh(
        np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0], [0, 1, 0], [1, 1, 0], [2, 1, 0]]),
        [("quad", np.array([[0, 1, 4, 3], [1, 2, 5, 4]], dtype=np.int64))],
    )
    out, report = mio.split_interface(
        mesh, Region("cut", "side", [[0, 1]]), add_cohesive=True, return_report=True
    )
    assert report == {"num_duplicated_points": 2, "num_cohesive_cells": 1}
    assert out.cells[-1].type == "line"
    assert out.cells[-1].data.shape == (1, 2)
    assert out.cell_data["cohesive:trace_b"][-1].shape == (1, 2)
    assert (
        out.cell_data["cohesive:trace_b"][-1][0].tolist()
        != out.cells[-1].data[0].tolist()
    )

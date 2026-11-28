"""Tests for Roadmap Tier 9 Wave 2 operations (9.2, 9.5, 9.6.1, 9.7.1, 9.7.3)."""

import numpy as np
import pytest

import meshioplusplus as mio


def test_subdivide_surface_loop_and_catmull_clark():
    vol = mio.box((1, 1, 1), cell_type="tetra")
    skin = mio.extract_surface(vol)
    skin.point_data["x"] = skin.points[:, 0].copy()

    out, rep = mio.subdivide_surface(
        skin,
        scheme="loop",
        levels=2,
        crease_angle=0.0,
        record_parent_ids=True,
        return_report=True,
    )
    assert out.cells[0].type == "triangle"
    assert len(out.cells[0].data) == len(skin.cells[0].data) * 16
    assert "subdivide:parent_cell" in out.cell_data
    assert "x" in out.point_data

    qskin = mio.extract_surface(mio.box((1, 1, 1), cell_type="hexahedron"))
    qout, qrep = mio.subdivide_surface(
        qskin, scheme="catmull_clark", levels=1, crease_angle=60.0, return_report=True
    )
    assert qout.cells[0].type == "quad"
    assert len(qout.cells[0].data) == 24
    assert qrep["num_crease_edges"] == 12

    with pytest.raises(ValueError, match="triangles only"):
        mio.subdivide_surface(qskin, scheme="loop")


def test_convert_cells_order3_and_full_lagrange():
    tet = mio.box((1, 1, 1), cell_type="tetra")
    tet.point_data["f"] = (
        2.0 * tet.points[:, 0] - tet.points[:, 1] + 3.0 * tet.points[:, 2]
    )

    t20 = mio.convert_cells(tet, "elevate", order=3)
    assert t20.cells[0].type == "tetra20"
    expected = 2.0 * t20.points[:, 0] - t20.points[:, 1] + 3.0 * t20.points[:, 2]
    np.testing.assert_allclose(t20.point_data["f"], expected, atol=1e-12)

    t10 = mio.convert_cells(t20, "reduce", order=2)
    assert t10.cells[0].type == "tetra10"

    t4 = mio.convert_cells(t10, "reduce", order=1)
    assert t4.cells[0].type == "tetra"
    assert len(t4.points) == len(tet.points)

    hex_m = mio.box((1, 1, 1), cell_type="hexahedron")
    h27 = mio.convert_cells(hex_m, "elevate", order=2, full=True)
    assert h27.cells[0].type == "hexahedron27"
    assert len(h27.points) == 27


def test_refine_undo_green_first_and_single_mesh_undo():
    base = mio.box((2, 2, 2), cell_type="tetra")
    p1 = mio.refine(
        base, levels=1, cells=[0], record_levels=True, record_hierarchy=True
    )
    undone, rep = mio.undo_green(p1, return_report=True)
    assert rep["num_groups_undone"] > 0
    assert len(undone.cells[0].data) < len(p1.cells[0].data)

    cur = p1
    for _ in range(3):
        cur = mio.refine(
            cur,
            levels=1,
            cells=[0],
            record_levels=True,
            record_hierarchy=True,
            undo_green_first=True,
        )
    qm = mio.attach_quality(cur)
    assert float(np.min(qm.cell_data["quality:scaled_jacobian"][0])) > 0.02


def test_partition_connected_parts_and_unknown_provider():
    grid = mio.box((4, 4, 1), cell_type="hexahedron")
    lbls = mio.partition_labels(grid, 4, method="sfc", connected_parts=True)
    assert len(lbls) == 1
    assert len(lbls[0]) == len(grid.cells[0].data)

    with pytest.raises(ValueError, match="unknown or unavailable method"):
        mio.partition_labels(grid, 2, method="nonexistent_provider_xyz")

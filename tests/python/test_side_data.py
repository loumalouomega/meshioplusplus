"""Tests for roadmap §8.4 facet/edge data (sides_to_cells / cells_to_sides)."""

import numpy as np
import pytest

import meshioplusplus as mp


def _hex_with_side_set():
    m = mp.grid((2, 1, 1))
    # Two boundary quad faces on cell 0 (local facet 0) and cell 1 (local facet 0)
    sides = np.array([[0, 0], [1, 0]], dtype=np.int64)
    m.regions = [mp.Region("bottom", "side", sides)]
    m.field_data["exodus:sset:bottom:pressure"] = np.array(
        [10.0, 20.0], dtype=np.float64
    )
    return m


def test_sides_to_cells_and_cells_to_sides_roundtrip():
    m = _hex_with_side_set()
    mat, report = mp.sides_to_cells(m, "bottom", return_report=True)
    assert report["num_facets"] == 2
    assert len(mat.cells) == 2
    assert mat.cells[1].type == "quad"
    assert mat.cells[1].data.shape == (2, 4)
    assert "meshio:side_of" in mat.cell_data
    np.testing.assert_array_equal(mat.cell_data["meshio:side_of"][1], [[0, 0], [1, 0]])
    assert "pressure" in mat.cell_data
    np.testing.assert_allclose(mat.cell_data["pressure"][1], [10.0, 20.0])

    # Modify facet cell_data and convert back to Side region + exodus:sset:*
    mat.cell_data["pressure"][1] = np.array([30.0, 40.0], dtype=np.float64)
    back, back_rep = mp.cells_to_sides(mat, drop_facet_blocks=True, return_report=True)
    assert back_rep["num_sides"] == 2
    assert len(back.cells) == 1
    assert "meshio:side_of" not in back.cell_data
    side_regs = [r for r in back.regions if r.kind == "side" and r.name == "bottom"]
    assert len(side_regs) == 1
    np.testing.assert_array_equal(side_regs[0].entries, [[0, 0], [1, 0]])
    np.testing.assert_allclose(
        back.field_data["exodus:sset:bottom:pressure"], [30.0, 40.0]
    )


def test_exodus_side_set_variable_roundtrip_via_sides_to_cells(tmp_path):
    """Roadmap 8.4 Done when: Exodus side-set variable reads as cell data on a
    quad block that sides_to_cells produces, and writes back."""
    m = _hex_with_side_set()
    exo1 = tmp_path / "in.exo"
    exo2 = tmp_path / "out.exo"
    mp.write(exo1, m)

    read_m = mp.read(exo1)
    mat = mp.sides_to_cells(read_m)
    assert mat.cells[-1].type == "quad"
    np.testing.assert_allclose(mat.cell_data["pressure"][-1], [10.0, 20.0])

    back = mp.cells_to_sides(mat, drop_facet_blocks=True)
    mp.write(exo2, back)
    final_m = mp.read(exo2)
    np.testing.assert_allclose(
        final_m.field_data["exodus:sset:bottom:pressure"], [10.0, 20.0]
    )


def test_sides_to_cells_unknown_region_raises():
    m = _hex_with_side_set()
    with pytest.raises(ValueError, match="no side region named"):
        mp.sides_to_cells(m, "missing")

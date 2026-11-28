"""Tests for roadmap §8.6 structured-grid annotations (meshio:structured:<block>)."""

import numpy as np
import pytest

import meshioplusplus as mp
from meshioplusplus.vti import _vti as py_vti
from meshioplusplus.vtr import _vtr as py_vtr
from meshioplusplus.vts import _vts as py_vts


def _make_curvilinear_vts(path):
    """Write a genuinely curved 3x3x2 StructuredGrid XML file directly."""
    ni, nj, nk = 3, 3, 2
    pts = []
    for k in range(nk):
        for j in range(nj):
            for i in range(ni):
                # Non-linear warp so geometric uniform lattice detection rejects it.
                x = float(i) + 0.15 * np.sin(j)
                y = float(j) + 0.1 * (i * i)
                z = float(k) + 0.05 * i * j
                pts.append(f"{x:.12g} {y:.12g} {z:.12g}")
    xml = (
        '<?xml version="1.0"?>\n'
        '<VTKFile type="StructuredGrid" version="0.1" byte_order="LittleEndian">\n'
        f'<StructuredGrid WholeExtent="0 {ni - 1} 0 {nj - 1} 0 {nk - 1}">\n'
        f'<Piece Extent="0 {ni - 1} 0 {nj - 1} 0 {nk - 1}">\n'
        "<Points>\n"
        '<DataArray type="Float64" Name="Points" NumberOfComponents="3" format="ascii">\n'
        + "\n".join(pts)
        + "\n</DataArray>\n</Points>\n</Piece>\n</StructuredGrid>\n</VTKFile>\n"
    )
    path.write_text(xml)


def _make_graded_vtr(path):
    """Write a non-uniform (graded) 4x3x2 RectilinearGrid XML file directly."""
    xs = [0.0, 0.1, 0.4, 1.0]
    ys = [0.0, 0.25, 1.0]
    zs = [0.0, 1.0]
    xml = (
        '<?xml version="1.0"?>\n'
        '<VTKFile type="RectilinearGrid" version="0.1" byte_order="LittleEndian">\n'
        '<RectilinearGrid WholeExtent="0 3 0 2 0 1">\n'
        '<Piece Extent="0 3 0 2 0 1">\n'
        "<Coordinates>\n"
        '<DataArray type="Float64" Name="x_coordinates" format="ascii">\n'
        + "\n".join(str(v) for v in xs)
        + '\n</DataArray>\n<DataArray type="Float64" Name="y_coordinates" format="ascii">\n'
        + "\n".join(str(v) for v in ys)
        + '\n</DataArray>\n<DataArray type="Float64" Name="z_coordinates" format="ascii">\n'
        + "\n".join(str(v) for v in zs)
        + "\n</DataArray>\n</Coordinates>\n</Piece>\n</RectilinearGrid>\n</VTKFile>\n"
    )
    path.write_text(xml)
    return xs, ys, zs


@pytest.mark.parametrize("engine", ["core", "py"])
def test_curvilinear_vts_roundtrip_and_smooth_invalidation(tmp_path, engine):
    a_vts = tmp_path / "a.vts"
    b_vts = tmp_path / "b.vts"
    _make_curvilinear_vts(a_vts)

    if engine == "core":
        m = mp.read(a_vts)
    else:
        m = py_vts.read(str(a_vts))

    assert "meshio:structured:0" in m.field_data
    np.testing.assert_array_equal(m.field_data["meshio:structured:0"], [3, 3, 2])

    # Round-trip curvilinear grid through vts writer
    if engine == "core":
        mp.write(b_vts, m)
        m2 = mp.read(b_vts)
    else:
        py_vts.write(str(b_vts), m)
        m2 = py_vts.read(str(b_vts))
    np.testing.assert_allclose(m2.points, m.points)
    np.testing.assert_array_equal(m2.cells[0].data, m.cells[0].data)

    # After smooth, structured annotation is stripped and write_vts refuses by name.
    smoothed = mp.smooth(m, iterations=1)
    assert "meshio:structured:0" not in smoothed.field_data
    with pytest.raises(mp.WriteError, match="StructuredGrid"):
        if engine == "core":
            mp.write(tmp_path / "c.vts", smoothed)
        else:
            py_vts.write(str(tmp_path / "c.vts"), smoothed)


@pytest.mark.parametrize("engine", ["core", "py"])
def test_graded_vtr_roundtrip(tmp_path, engine):
    a_vtr = tmp_path / "a.vtr"
    b_vtr = tmp_path / "b.vtr"
    _make_graded_vtr(a_vtr)

    if engine == "core":
        m = mp.read(a_vtr)
        mp.write(b_vtr, m)
        m2 = mp.read(b_vtr)
    else:
        m = py_vtr.read(str(a_vtr))
        py_vtr.write(str(b_vtr), m)
        m2 = py_vtr.read(str(b_vtr))

    np.testing.assert_array_equal(m.field_data["meshio:structured:0"], [4, 3, 2])
    np.testing.assert_allclose(m2.points, m.points)
    np.testing.assert_array_equal(m2.cells[0].data, m.cells[0].data)


@pytest.mark.parametrize("engine", ["core", "py"])
def test_vti_sets_structured_annotation(tmp_path, engine):
    g = mp.grid((2, 3, 2))
    out = tmp_path / f"g_{engine}.vti"
    if engine == "core":
        mp.write(out, g)
        back = mp.read(out)
    else:
        py_vti.write(str(out), g)
        back = py_vti.read(str(out))
    np.testing.assert_array_equal(back.field_data["meshio:structured:0"], [3, 4, 3])


def test_to_grid_and_from_grid_roundtrip():
    g = mp.grid((2, 3, 2), origin=(1.0, 2.0, 3.0), spacing=(0.5, 0.5, 0.5))
    g.point_data["T"] = np.arange(len(g.points), dtype=np.float64)
    ga = mp.to_grid(g, fields=["T"])
    assert ga.values.shape == (1, 3, 4, 3)
    m = mp.from_grid(ga)
    assert "meshio:structured:0" in m.field_data
    np.testing.assert_array_equal(m.field_data["meshio:structured:0"], [3, 4, 3])
    np.testing.assert_allclose(m.point_data["T"].reshape(-1), g.point_data["T"])

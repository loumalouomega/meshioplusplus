"""Tests for roadmap 8.3: integration-point data convention, ops, and MED ELGA."""

from __future__ import annotations

import numpy as np
import pytest

import meshioplusplus as mio


def _unit_hex_mesh() -> mio.Mesh:
    points = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
        ],
        dtype=np.float64,
    )
    cells = [("hexahedron", np.array([[0, 1, 2, 3, 4, 5, 6, 7]], dtype=np.int64))]
    return mio.Mesh(points, cells)


def _unit_tet_mesh() -> mio.Mesh:
    points = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ],
        dtype=np.float64,
    )
    cells = [("tetra", np.array([[0, 1, 2, 3]], dtype=np.int64))]
    return mio.Mesh(points, cells)


def _two_quad_mesh() -> mio.Mesh:
    points = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
            [2.0, 1.0, 0.0],
        ],
        dtype=np.float64,
    )
    cells = [
        (
            "quad",
            np.array([[0, 1, 4, 3], [1, 2, 5, 4]], dtype=np.int64),
        )
    ]
    return mio.Mesh(points, cells)


@pytest.mark.parametrize("factory", [_unit_hex_mesh, _unit_tet_mesh, _two_quad_mesh])
def test_linear_field_points_to_ip_and_back(factory):
    mesh = factory()
    x = mesh.points[:, 0]
    y = mesh.points[:, 1]
    z = mesh.points[:, 2]
    linear_scalar = 2.0 * x - 3.0 * y + 5.0 * z + 7.0
    linear_vec = np.column_stack([linear_scalar, -x + 4.0 * y + 1.0])
    mesh.point_data["temp"] = linear_scalar
    mesh.point_data["stress"] = linear_vec

    ip_mesh, rep1 = mio.points_to_ip(
        mesh, rule_id=2, drop_source=True, return_report=True
    )
    assert set(rep1["converted"]) == {"temp", "stress"}
    assert "temp" not in ip_mesh.point_data
    assert "meshio:ip:temp" in ip_mesh.field_data
    assert "meshio:ip:stress" in ip_mesh.field_data

    layout_s = np.asarray(ip_mesh.field_data["meshio:ip:stress"]).ravel()
    assert int(layout_s[1]) == 2
    assert int(layout_s[2]) == 2

    # Extrapolating back recovers the exact linear nodal field.
    rec_mesh, rep2 = mio.ip_to_points(ip_mesh, drop_source=False, return_report=True)
    assert set(rep2["converted"]) == {"temp", "stress"}
    np.testing.assert_allclose(rec_mesh.point_data["temp"], linear_scalar, atol=1e-11)
    np.testing.assert_allclose(rec_mesh.point_data["stress"], linear_vec, atol=1e-11)

    # Cell average matches the exact field value at each cell centroid.
    cell_mesh = mio.ip_to_cells(ip_mesh)
    assert "meshio:ip:temp" not in cell_mesh.field_data
    centroids = np.array(
        [mesh.points[cell].mean(axis=0) for cell in mesh.cells[0].data],
        dtype=np.float64,
    )
    expected_c = (
        2.0 * centroids[:, 0] - 3.0 * centroids[:, 1] + 5.0 * centroids[:, 2] + 7.0
    )
    np.testing.assert_allclose(cell_mesh.cell_data["temp"][0], expected_c, atol=1e-11)


def test_element_nodal_and_reduced_integration():
    mesh = _unit_hex_mesh()
    vals = np.linspace(1.0, 8.0, 8, dtype=np.float64)
    mesh.point_data["u"] = vals

    # rule_id=0: element-nodal round-trip
    elno = mio.points_to_ip(mesh, rule_id=0, drop_source=True)
    assert list(np.asarray(elno.field_data["meshio:ip:u"]).ravel()) == [8, 1, 0]
    back = mio.ip_to_points(elno)
    np.testing.assert_allclose(back.point_data["u"], vals, atol=1e-12)

    # rule_id=1: 1-point centroid rule preserves constant fields exactly
    mesh_const = _unit_hex_mesh()
    mesh_const.point_data["c"] = np.full(8, 42.5, dtype=np.float64)
    ip1 = mio.points_to_ip(mesh_const, rule_id=1, drop_source=True)
    assert list(np.asarray(ip1.field_data["meshio:ip:c"]).ravel()) == [1, 1, 1]
    back1 = mio.ip_to_points(ip1)
    np.testing.assert_allclose(back1.point_data["c"], 42.5, atol=1e-12)


def test_data_info_reports_integration_points():
    mesh = _unit_tet_mesh()
    mesh.point_data["s"] = np.array([1.0, 2.0, 3.0, 4.0])
    ip_mesh = mio.points_to_ip(mesh, rule_id=2, drop_source=True)

    info = mio.data_info(ip_mesh)
    cell_entries = [
        e for e in info if e["location"] == "cell_data" and e["name"] == "s"
    ]
    assert len(cell_entries) == 1
    assert cell_entries[0]["integration_points"] == {
        "points": 4,
        "width": 1,
        "rule_id": 2,
    }


def test_med_elga_native_roundtrip(tmp_path):
    mesh = _unit_hex_mesh()
    x, y, z = mesh.points[:, 0], mesh.points[:, 1], mesh.points[:, 2]
    stress = np.column_stack([10.0 * x + 1.0, 20.0 * y + 2.0, 30.0 * z + 3.0])
    mesh.point_data["S"] = stress
    ip_mesh = mio.points_to_ip(mesh, rule_id=2, drop_source=True)

    med_path = tmp_path / "elga.med"
    mio.write(med_path, ip_mesh, file_format="med")

    # Verify NGA == 8 was written in the HDF5 MED file.
    h5py = pytest.importorskip("h5py")
    with h5py.File(med_path, "r") as f:
        ts_name = next(iter(f["CHA"]["S"].keys()))
        mai = f["CHA"]["S"][ts_name]["MAI.HE8"]
        pfl_name = next(iter(mai.keys()))
        assert int(mai[pfl_name].attrs["NGA"]) == 8

    # Read back natively and verify meshio:ip:S and exact nodal recovery.
    reloaded = mio.read(med_path, file_format="med")
    assert "meshio:ip:S" in reloaded.field_data
    layout = np.asarray(reloaded.field_data["meshio:ip:S"]).ravel()
    assert int(layout[0]) == 8
    assert int(layout[1]) == 3
    np.testing.assert_allclose(
        reloaded.cell_data["S"][0], ip_mesh.cell_data["S"][0], atol=1e-12
    )

    recovered = mio.ip_to_points(reloaded)
    np.testing.assert_allclose(recovered.point_data["S"], stress, atol=1e-11)

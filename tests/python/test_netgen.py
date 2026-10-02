import pathlib

import numpy as np
import pytest

import meshioplusplus

from . import helpers

test_set = [
    helpers.empty_mesh,
    helpers.line_mesh,
    helpers.tri_mesh_2d,
    helpers.tri_mesh,
    helpers.triangle6_mesh,
    helpers.quad_mesh,
    helpers.quad8_mesh,
    helpers.tri_quad_mesh,
    helpers.tet_mesh,
    helpers.tet10_mesh,
    helpers.hex_mesh,
    helpers.hex20_mesh,
    helpers.pyramid_mesh,
    helpers.wedge_mesh,
    helpers.add_cell_data(helpers.tri_mesh, [("netgen:index", (), int)]),
]


this_dir = pathlib.Path(__file__).resolve().parent
netgen_mesh_directory = this_dir / "meshes" / "netgen"

PERIODIC_1D = "periodic_1d.vol"
PERIODIC_2D = "periodic_2d.vol"
PERIODIC_3D = "periodic_3d.vol"

netgen_meshes = [PERIODIC_1D, PERIODIC_2D, PERIODIC_3D]


@pytest.mark.parametrize("mesh", test_set)
@pytest.mark.parametrize("suffix", [".vol", ".vol.gz"])
def test(mesh, suffix, tmp_path):
    helpers.write_read(
        tmp_path,
        meshioplusplus.netgen.write,
        meshioplusplus.netgen.read,
        mesh,
        1.0e-13,
        extension=suffix,
    )


expected_identifications = {
    PERIODIC_1D: np.array([[1, 51, 1]]),
    PERIODIC_2D: np.array(
        [[2, 1, 4], [3, 4, 4], [9, 17, 4], [10, 18, 4], [11, 19, 4], [12, 20, 4]]
    ),
    PERIODIC_3D: np.array(
        [
            [1, 3, 1],
            [2, 5, 1],
            [4, 7, 1],
            [6, 8, 1],
            [9, 11, 1],
            [10, 12, 1],
            [15, 13, 1],
            [16, 14, 1],
            [21, 19, 1],
            [22, 20, 1],
            [25, 23, 1],
            [26, 24, 1],
            [38, 54, 1],
            [39, 55, 1],
            [40, 56, 1],
        ]
    ),
}

expected_identificationtypes = {
    PERIODIC_1D: np.array([[2]]),
    PERIODIC_2D: np.array([[1, 1, 1, 2]]),
    PERIODIC_3D: np.array([[2]]),
}

expected_field_data = {
    PERIODIC_1D: {},
    PERIODIC_2D: {"outer": [3, 1], "periodic": [4, 1]},
    PERIODIC_3D: {"outer": [6, 2], "default": [3, 2]},
}


@pytest.mark.parametrize("netgen_mesh", [PERIODIC_1D, PERIODIC_2D, PERIODIC_3D])
def test_advanced(netgen_mesh, tmp_path):
    mesh = meshioplusplus.read(str(netgen_mesh_directory / netgen_mesh))

    p = tmp_path / f"{netgen_mesh}_out.vol"
    mesh.write(p)
    mesh_out = meshioplusplus.read(p)

    assert np.all(
        mesh.info["netgen:identifications"] == expected_identifications[netgen_mesh]
    )
    assert np.all(
        mesh.info["netgen:identifications"] == mesh_out.info["netgen:identifications"]
    )
    assert np.all(
        mesh.info["netgen:identificationtypes"]
        == expected_identificationtypes[netgen_mesh]
    )
    assert np.all(
        mesh.info["netgen:identificationtypes"]
        == mesh_out.info["netgen:identificationtypes"]
    )
    for kk, vv in mesh.field_data.items():
        assert np.all(vv == expected_field_data[netgen_mesh][kk])
        assert np.all(vv == mesh_out.field_data[kk])
    for cd, cd_out in zip(
        mesh.cell_data["netgen:index"], mesh_out.cell_data["netgen:index"]
    ):
        assert np.all(cd == cd_out)


@pytest.mark.parametrize("suffix", [".vol", ".vol.gz"])
def test_extras_use_native_core(tmp_path, monkeypatch, suffix):
    core = pytest.importorskip("meshioplusplus._core")
    if suffix == ".vol.gz" and not core.__has_zlib__:
        pytest.skip("gzip requires native zlib")
    mesh = meshioplusplus.Mesh(
        np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
        [("triangle", np.array([[0, 1, 2]]))],
        field_data={"wall": [2, 2], "edge": [1, 1]},
        cell_data={"netgen:index": [np.array([2])]},
        info={
            "netgen:identifications": np.array([[1, 2, 1]]),
            "netgen:identificationtypes": np.array([[2]]),
        },
    )
    path = tmp_path / ("extras" + suffix)
    monkeypatch.setenv("MESHIOPLUSPLUS_STRICT_CORE", "1")
    meshioplusplus.write(path, mesh)
    native = meshioplusplus.read(path)
    reference = meshioplusplus.netgen._netgen.read(path)
    np.testing.assert_array_equal(native.points, reference.points)
    for key in mesh.field_data:
        np.testing.assert_array_equal(native.field_data[key], reference.field_data[key])
    for key in mesh.info:
        np.testing.assert_array_equal(native.info[key], reference.info[key])


def test_native_two_line_edges_and_auxiliary_sections(tmp_path):
    core = pytest.importorskip("meshioplusplus._core")
    path = tmp_path / "two-line.vol"
    path.write_text(
        "mesh3d\ndimension\n3\nsurf1 surf2 p1 p2\nedgesegmentsgi2\n1\n"
        "7 0 1 2\n0 0 0 0\npoints\n2\n0 0 0\n1 0 0\n"
        "face_colours\n1\n1 0 0\nsingular_points\n1\n1\nendmesh\n"
    )
    native = core.netgen_read(str(path))
    reference = meshioplusplus.netgen._netgen.read(path)
    np.testing.assert_array_equal(native.cells[0].data, reference.cells[0].data)
    np.testing.assert_array_equal(native.cell_data["netgen:index"][0], [7])


@pytest.mark.parametrize(
    "section", ["identifications", "identificationtypes", "materials"]
)
def test_native_extras_reject_truncated_sections(tmp_path, section):
    core = pytest.importorskip("meshioplusplus._core")
    path = tmp_path / "truncated.vol"
    path.write_text(f"mesh3d\ndimension\n3\n{section}\n1\n")
    with pytest.raises(Exception, match="Netgen"):
        core.netgen_read(str(path))

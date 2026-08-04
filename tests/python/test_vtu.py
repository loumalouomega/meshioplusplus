import pathlib
import tempfile

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus.vtu import _vtu

from . import helpers

test_set = [
    # helpers.empty_mesh,
    helpers.line_mesh,
    helpers.tri_mesh,
    helpers.tri_mesh_one_cell,
    helpers.triangle6_mesh,
    helpers.quad_mesh,
    helpers.quad8_mesh,
    helpers.tri_quad_mesh,
    helpers.polygon_mesh,
    helpers.polygon_mesh_one_cell,
    helpers.polygon2_mesh,
    helpers.tet_mesh,
    helpers.tet10_mesh,
    helpers.hex_mesh,
    helpers.hex20_mesh,
    helpers.pyramid_mesh,
    helpers.wedge_mesh,
    helpers.polyhedron_mesh,
    helpers.lagrange_high_order_mesh,
    helpers.add_point_data(helpers.tri_mesh, 1),
    helpers.add_point_data(helpers.tri_mesh, 2),
    helpers.add_point_data(helpers.tri_mesh, 3),
    helpers.add_cell_data(helpers.tri_mesh, [("a", (), np.float64)]),
    helpers.add_cell_data(helpers.tri_quad_mesh, [("a", (), np.float64)]),
    helpers.add_cell_data(helpers.tri_mesh, [("a", (2,), np.float32)]),
    helpers.add_cell_data(helpers.tri_mesh, [("b", (3,), np.float64)]),
    helpers.add_cell_data(helpers.polygon_mesh, [("a", (), np.float32)]),
    helpers.add_cell_data(helpers.polyhedron_mesh, [("a", (2,), np.float32)]),
]


@pytest.mark.parametrize("mesh", test_set)
@pytest.mark.parametrize(
    "data_type", [(False, None), (True, None), (True, "lzma"), (True, "zlib")]
)
def test(mesh, data_type, tmp_path):
    binary, compression = data_type

    def writer(*args, **kwargs):
        return meshioplusplus.vtu.write(
            *args, binary=binary, compression=compression, **kwargs
        )

    # ASCII files are only meant for debugging, VTK stores only 11 digits
    # <https://gitlab.kitware.com/vtk/vtk/-/issues/17038#note_264052>
    tol = 1.0e-15 if binary else 1.0e-10
    helpers.write_read(tmp_path, writer, meshioplusplus.vtu.read, mesh, tol)


def test_generic_io(tmp_path):
    helpers.generic_io(tmp_path / "test.vtu")
    # With additional, insignificant suffix:
    helpers.generic_io(tmp_path / "test.0.vtu")


@pytest.mark.parametrize(
    "filename, ref_cells, ref_num_cells, ref_num_pnt",
    [
        ("00_raw_binary.vtu", "tetra", 162, 64),
        ("01_raw_binary_int64.vtu", "tetra", 162, 64),
        ("02_raw_compressed.vtu", "tetra", 162, 64),
    ],
)
def test_read_from_file(filename, ref_cells, ref_num_cells, ref_num_pnt):
    this_dir = pathlib.Path(__file__).resolve().parent
    filename = this_dir / "meshes" / "vtu" / filename

    mesh = meshioplusplus.read(filename)
    assert len(mesh.cells) == 1
    assert ref_cells == mesh.cells[0].type
    assert len(mesh.cells[0].data) == ref_num_cells
    assert len(mesh.points) == ref_num_pnt


# --- malformed-input / error-path coverage ---


def test_vtu_wrong_root_tag_raises(tmp_path):
    # A well-formed XML document that is not a VTKFile is rejected on both the
    # C++ and Python paths, so the ReadError surfaces through the shim.
    p = tmp_path / "bad.vtu"
    p.write_text('<?xml version="1.0"?><NotVTK></NotVTK>')
    with pytest.raises(meshioplusplus.ReadError):
        meshioplusplus.read(p, file_format="vtu")


def test_vtu_not_xml_raises(tmp_path):
    p = tmp_path / "bad.vtu"
    p.write_text("this is not xml at all")
    with pytest.raises(meshioplusplus.ReadError):
        meshioplusplus.read(p, file_format="vtu")


def test_polyhedra_mix_with_other_cell_types_both_engines():
    """VTU's `faceoffsets` carries -1 for a non-polyhedral cell, which is
    exactly how the format expresses a mixed mesh -- and an OpenFOAM mesh always
    mixes. The C++ and Python paths must agree, because Windows CI builds native
    paths off and runs the reference implementation.
    """
    pts = np.array(
        [
            [0, 0, 0],
            [1, 0, 0],
            [1, 1, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 0, 1],
            [1, 1, 1],
            [0, 1, 1],
        ],
        float,
    )
    tet_faces = [
        np.array([0, 2, 1]),
        np.array([0, 1, 3]),
        np.array([1, 2, 3]),
        np.array([2, 0, 3]),
    ]
    mesh = meshioplusplus.Mesh(
        pts,
        [
            ("triangle", np.array([[4, 5, 6]])),
            ("polyhedron4", [tet_faces]),
        ],
    )

    for writer, reader, tag in (
        (_vtu.write, _vtu.read, "python"),
        (meshioplusplus.vtu.write, meshioplusplus.vtu.read, "shim"),
    ):
        with tempfile.TemporaryDirectory() as d:
            p = pathlib.Path(d) / "mixed.vtu"
            writer(p, mesh)
            out = reader(p)
            types = [c.type for c in out.cells]
            assert any(t.startswith("polyhedron") for t in types), (tag, types)
            assert "triangle" in types, (tag, types)

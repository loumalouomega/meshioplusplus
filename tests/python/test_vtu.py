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


def _wedge_faces(b):
    return [
        np.array([b + 0, b + 2, b + 1]),
        np.array([b + 3, b + 4, b + 5]),
        np.array([b + 0, b + 1, b + 4, b + 3]),
        np.array([b + 1, b + 2, b + 5, b + 4]),
        np.array([b + 2, b + 0, b + 3, b + 5]),
    ]


def _tet_faces(b):
    return [
        np.array([b + 0, b + 2, b + 1]),
        np.array([b + 0, b + 1, b + 3]),
        np.array([b + 1, b + 2, b + 3]),
        np.array([b + 2, b + 0, b + 3]),
    ]


def test_polyhedra_with_alternating_node_counts_keep_their_cell_data():
    # A polyhedral run that mixes node counts is bucketed into one block per count, so a
    # bucket's cells are not contiguous in the file. cell_data has to follow the cells.
    pts = np.random.default_rng(0).random((16, 3))
    mesh = meshioplusplus.Mesh(
        pts,
        [
            ("polyhedron6", [_wedge_faces(0)]),
            ("polyhedron4", [_tet_faces(6)]),
            ("polyhedron6", [_wedge_faces(10)]),
        ],
        cell_data={"tag": [np.array([1.0]), np.array([2.0]), np.array([3.0])]},
    )

    for writer, reader, label in (
        (_vtu.write, _vtu.read, "python"),
        (meshioplusplus.vtu.write, meshioplusplus.vtu.read, "shim"),
    ):
        with tempfile.TemporaryDirectory() as d:
            p = pathlib.Path(d) / "alt.vtu"
            writer(p, mesh)
            out = reader(p)
        by_type = {}
        for cell, tag in zip(out.cells, out.cell_data["tag"]):
            by_type[cell.type] = np.asarray(tag).tolist()
        assert by_type == {"polyhedron6": [1.0, 3.0], "polyhedron4": [2.0]}, (
            label,
            by_type,
        )


# tools/gen_feconv_quirk_fixtures.py: one pyramid on one polyhedral cube, in
# every framing FEconv's samples use (and base64 appended, which they do not).
_QUIRKS = pathlib.Path(__file__).resolve().parent / "meshes" / "vtu"


def _core_vtu_read(path):
    from meshioplusplus import _core

    return _core.vtu_read(str(path))


@pytest.mark.parametrize("reader", [_core_vtu_read, _vtu.read], ids=["core", "python"])
@pytest.mark.parametrize(
    "name", ["raw_bigendian", "raw_zlib", "base64_appended", "binary_bigendian"]
)
def test_appended_and_bigendian_framings(reader, name):
    from meshioplusplus import _core

    if (
        reader is _core_vtu_read
        and name == "raw_zlib"
        and not getattr(_core, "__has_zlib__", False)
    ):
        pytest.skip("this build has no zlib")
    mesh = reader(_QUIRKS / f"{name}.vtu")
    assert [c.type for c in mesh.cells] == ["pyramid", "polyhedron8"]
    np.testing.assert_array_equal(mesh.cells[0].data, [[4, 5, 6, 7, 8]])
    faces = [np.asarray(f).tolist() for f in mesh.cells[1].data[0]]
    assert faces[0] == [3, 0, 4, 7] and faces[-1] == [4, 5, 6, 7]
    np.testing.assert_array_equal(mesh.points[8], [0.5, 0.5, 2.0])
    np.testing.assert_array_equal(mesh.point_data["heat"], np.arange(9.0))
    assert [np.asarray(d).tolist() for d in mesh.cell_data["material"]] == [[7], [9]]
    assert mesh.cell_data["material"][0].dtype.isnative


@pytest.mark.parametrize("reader", [_core_vtu_read, _vtu.read], ids=["core", "python"])
def test_two_pieces_with_polyhedra_merge(reader):
    mesh = reader(_QUIRKS / "two_pieces_polyhedron.vtu")
    assert len(mesh.points) == 18
    assert [c.type for c in mesh.cells] == ["pyramid", "polyhedron8"] * 2
    # The second piece's node ids follow the first piece's nine points.
    np.testing.assert_array_equal(mesh.cells[2].data, [[13, 14, 15, 16, 17]])
    faces = [np.asarray(f).tolist() for f in mesh.cells[3].data[0]]
    assert faces[0] == [12, 9, 13, 16]
    # Cell data in both pieces survives; point data in only one is dropped.
    assert [np.asarray(d).tolist() for d in mesh.cell_data["material"]] == [
        [7],
        [9],
    ] * 2
    assert "heat" not in (mesh.point_data or {})


def test_two_pieces_metadata_matches_the_read():
    from meshioplusplus import _core

    meta = _core.read_metadata(str(_QUIRKS / "two_pieces_polyhedron.vtu"), "vtu")
    assert meta["num_points"] == 18
    assert meta["cell_data_names"] == ["material"]
    assert meta["point_data_names"] == []


def test_faceoffsets_after_an_ordinary_block_end_at_the_stream(tmp_path):
    # faceoffsets are END offsets into `faces`. After a non-polyhedral block
    # (-1 entries) the next polyhedron's offset used to be rebased on that -1,
    # one short, and VTK (9.4+, ParaView 6) refused the whole file.
    import re

    pts = np.random.default_rng(1).random((20, 3))
    mesh = meshioplusplus.Mesh(
        pts,
        [
            ("tetra", [[0, 1, 2, 3]]),
            ("polyhedron4", [_tet_faces(4)]),
            ("hexahedron", [list(range(8, 16))]),
            ("polyhedron4", [_tet_faces(16)]),
        ],
    )
    p = tmp_path / "mixed.vtu"
    _vtu.write(p, mesh, binary=False)
    text = p.read_text()

    def array(name):
        m = re.search(rf'Name="{name}"[^>]*>(.*?)</DataArray>', text, re.S)
        return [int(v) for v in m.group(1).split()]

    faces, faceoffsets = array("faces"), array("faceoffsets")
    assert faceoffsets[0] == -1 and faceoffsets[2] == -1
    assert faceoffsets[3] == len(faces)
    assert faceoffsets[1] == 1 + 4 * 4  # one tetrahedron: count + 4 x (3 + 1)
    out = _vtu.read(p)
    assert [c.type for c in out.cells] == [c.type for c in mesh.cells]

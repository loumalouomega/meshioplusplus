import pathlib
import tempfile

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus.vtu import _vtu

from . import bezier_vtu_files, helpers

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
    helpers.bezier_mesh,
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


@pytest.mark.parametrize(
    "mesh",
    [
        helpers.tri_mesh,
        helpers.polyhedron_mesh,
        helpers.add_point_data(helpers.tri_mesh, 3),
        helpers.add_cell_data(helpers.tri_quad_mesh, [("a", (), np.float64)]),
    ],
)
@pytest.mark.parametrize("compression", [None, "zlib", "lzma"])
@pytest.mark.parametrize("engine", ["core", "python"])
def test_raw_appended_round_trip(mesh, compression, engine, tmp_path):
    # v16.21.0: every array in one raw <AppendedData> section. The Python
    # reference writes it too (lzma only there), and both engines read either.
    if (
        engine == "core"
        and compression == "lzma"
        and not getattr(meshioplusplus._core, "__has_lzma__", False)
    ):
        pytest.skip("this build has no lzma")
    writer = meshioplusplus.vtu.write if engine == "core" else _vtu.write

    def write(filename, m):
        writer(filename, m, binary=True, compression=compression, appended=True)
        text = pathlib.Path(filename).read_bytes()
        assert b'<AppendedData encoding="raw">' in text
        assert b'format="appended"' in text
        assert b'format="binary"' not in text

    helpers.write_read(tmp_path, write, meshioplusplus.vtu.read, mesh, 1.0e-15)
    helpers.write_read(tmp_path, write, _vtu.read, mesh, 1.0e-15)


def test_raw_appended_needs_binary(tmp_path):
    for writer in (meshioplusplus.vtu.write, _vtu.write):
        with pytest.raises(meshioplusplus.WriteError, match="needs binary=True"):
            writer(tmp_path / "a.vtu", helpers.tri_mesh, binary=False, appended=True)


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


def test_a_compressed_block_needs_zlib():
    # raw_zlib.vtu holds zlib-compressed appended data. The core reads it only
    # where zlib is built in and refuses it with a ReadError elsewhere.
    from meshioplusplus import _core

    path = _QUIRKS / "raw_zlib.vtu"
    if getattr(_core, "__has_zlib__", False):
        mesh = _core.vtu_read(str(path))
        assert [c.type for c in mesh.cells] == ["pyramid", "polyhedron8"]
    else:
        with pytest.raises(meshioplusplus.ReadError, match="zlib"):
            _core.vtu_read(str(path))


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


def test_uint64_above_int64_max_survives_an_ascii_round_trip(tmp_path):
    """Both C++ halves went through int64: the writer printed a UInt64 value
    above INT64_MAX as a negative number, and the reader clamped one."""
    from meshioplusplus import _core

    values = np.array([0, 2**63 + 5, 2**64 - 1], dtype=np.uint64)
    mesh = meshioplusplus.Mesh(
        np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
        [("triangle", np.array([[0, 1, 2]]))],
        point_data={"u": values},
    )
    cpp_path = tmp_path / "cpp.vtu"
    _core.vtu_write_codec(str(cpp_path), mesh, False, "none")
    assert str(2**64 - 1) in cpp_path.read_text()
    py_path = tmp_path / "py.vtu"
    _vtu.write(py_path, mesh, binary=False)

    for path in (cpp_path, py_path):
        for read in (lambda p: _core.vtu_read(str(p)), _vtu.read):
            got = np.asarray(read(path).point_data["u"])
            assert got.dtype == np.uint64
            np.testing.assert_array_equal(got, values)


# --- Bezier cells and their weights (v16.40.0) -------------------------------

_BEZIER_VTU = """<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian">
<UnstructuredGrid><Piece NumberOfPoints="9" NumberOfCells="1">
<Points><DataArray type="Float64" NumberOfComponents="3" format="ascii">
1 0 0 2 0 0 0 2 0 0 1 0 1.5 0 0 2 2 0 0 1.5 0 1 1 0 1.5 1.5 0
</DataArray></Points>
<Cells>
<DataArray type="Int64" Name="connectivity" format="ascii">0 1 2 3 4 5 6 7 8</DataArray>
<DataArray type="Int64" Name="offsets" format="ascii">9</DataArray>
<DataArray type="UInt8" Name="types" format="ascii">77</DataArray>
</Cells>
<PointData RationalWeights="w">
<DataArray type="Float64" Name="RationalWeights" format="ascii">0 0 0 0 0 0 0 0 0</DataArray>
<DataArray type="Float64" Name="w" format="ascii">
1 1 1 1 1 0.70710678118654757 1 0.70710678118654757 0.70710678118654757
</DataArray>
</PointData>
<CellData HigherOrderDegrees="deg">
<DataArray type="Float64" Name="deg" NumberOfComponents="3" format="ascii">2 2 0</DataArray>
</CellData>
</Piece></UnstructuredGrid></VTKFile>
"""


@pytest.mark.parametrize("reader", [_core_vtu_read, _vtu.read], ids=["core", "python"])
def test_bezier_marked_arrays_read_under_their_vtk_names(reader, tmp_path):
    # VTK marks the weights and degrees of Bezier cells with a section
    # attribute naming the array; both engines read the marked array under the
    # attribute's own name, and an unmarked array of that name gives way.
    p = tmp_path / "bezier.vtu"
    p.write_text(_BEZIER_VTU)
    mesh = reader(p)
    assert [c.type for c in mesh.cells] == ["VTK_BEZIER_QUADRILATERAL"]
    np.testing.assert_array_equal(mesh.cells[0].data, [np.arange(9)])
    assert sorted(mesh.point_data) == ["RationalWeights"]
    np.testing.assert_allclose(
        mesh.point_data["RationalWeights"],
        helpers.bezier_mesh.point_data["RationalWeights"],
    )
    assert sorted(mesh.cell_data) == ["HigherOrderDegrees"]
    np.testing.assert_array_equal(
        np.asarray(mesh.cell_data["HigherOrderDegrees"][0]).reshape(-1), [2, 2, 0]
    )


@pytest.mark.parametrize("engine", ["core", "python"])
def test_bezier_weights_are_marked_on_write(engine, tmp_path):
    p = tmp_path / "bezier.vtu"
    writer = meshioplusplus.vtu.write if engine == "core" else _vtu.write
    writer(p, helpers.bezier_mesh, binary=False)
    text = p.read_text()
    assert 'RationalWeights="RationalWeights"' in text
    assert 'HigherOrderDegrees="HigherOrderDegrees"' in text
    for read in (_core_vtu_read, _vtu.read):
        mesh = read(p)
        assert [c.type for c in mesh.cells] == [
            c.type for c in helpers.bezier_mesh.cells
        ]
        np.testing.assert_allclose(
            mesh.point_data["RationalWeights"],
            helpers.bezier_mesh.point_data["RationalWeights"],
        )
    # A mesh with neither array writes no marker.
    writer(p, helpers.tri_mesh, binary=False)
    assert "RationalWeights" not in p.read_text()


def test_vtk_reads_the_weights_we_mark(tmp_path):
    vtk = pytest.importorskip("vtk")
    p = tmp_path / "bezier.vtu"
    meshioplusplus.vtu.write(p, helpers.bezier_mesh)
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(p))
    reader.Update()
    grid = reader.GetOutput()
    assert grid.GetCellType(0) == 77
    weights = grid.GetPointData().GetRationalWeights()
    assert weights is not None and weights.GetName() == "RationalWeights"
    assert abs(weights.GetValue(5) - np.sqrt(0.5)) < 1e-15
    assert grid.GetCellData().GetHigherOrderDegrees() is not None


@pytest.mark.parametrize("reader", [_core_vtu_read, _vtu.read], ids=["core", "python"])
@pytest.mark.parametrize("encoding", ["ascii", "binary"])
def test_bezier_files_vtk_wrote(reader, encoding, tmp_path):
    # Written by VTK 9.7.1 (tools/gen_bezier_vtu_fixtures.py): the marked
    # arrays are named "Weights" and "Degrees" there.
    def vtk_file(name):
        p = tmp_path / f"{name}_{encoding}.vtu"
        p.write_bytes(getattr(bezier_vtu_files, f"{name}_{encoding}".upper()))
        return p

    quad = reader(vtk_file("quarter_annulus"))
    assert [c.type for c in quad.cells] == ["VTK_BEZIER_QUADRILATERAL"]
    np.testing.assert_allclose(quad.points, helpers.bezier_mesh.points)
    np.testing.assert_array_equal(
        quad.point_data["RationalWeights"],
        helpers.bezier_mesh.point_data["RationalWeights"],
    )
    np.testing.assert_array_equal(
        np.asarray(quad.cell_data["HigherOrderDegrees"][0]).reshape(-1), [2, 2, 0]
    )
    hexa = reader(vtk_file("anisotropic_hex"))
    assert [c.type for c in hexa.cells] == ["VTK_BEZIER_HEXAHEDRON"]
    assert hexa.cells[0].data.shape == (1, 12)
    np.testing.assert_array_equal(
        hexa.point_data["RationalWeights"], [1.0] * 8 + [2.0] * 4
    )
    np.testing.assert_array_equal(
        np.asarray(hexa.cell_data["HigherOrderDegrees"][0]).reshape(-1), [2, 1, 1]
    )


# --- the edge order of higher-order hexahedra across versions (VTK 9.1) -----


def _hex27(cell_type):
    pts = np.array([[i, i % 3, i % 5] for i in range(27)], dtype=float)
    return meshioplusplus.Mesh(pts, [(cell_type, np.arange(27).reshape(1, 27))])


@pytest.mark.parametrize("engine", ["core", "python"])
@pytest.mark.parametrize(
    "cell_type", ["VTK_LAGRANGE_HEXAHEDRON", "VTK_BEZIER_HEXAHEDRON"]
)
def test_higher_order_hexahedra_are_written_as_vtk_reads_them(
    engine, cell_type, tmp_path
):
    # VTK 9.1 swapped the nodes of edges 10 and 11 and reads a file older than
    # XML 2.2 the old way; meshio++ holds the new order, so it writes 2.2.
    vtk = pytest.importorskip("vtk")
    p = tmp_path / "hex.vtu"
    writer = meshioplusplus.vtu.write if engine == "core" else _vtu.write
    writer(p, _hex27(cell_type), binary=False)
    assert 'version="2.2"' in p.read_text()
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(p))
    reader.Update()
    cell = reader.GetOutput().GetCell(0)
    assert [cell.GetPointId(i) for i in range(27)] == list(range(27))
    for read in (_core_vtu_read, _vtu.read):
        np.testing.assert_array_equal(read(p).cells[0].data, [np.arange(27)])
    # Without such cells the version (and so every byte) is as before.
    writer(p, helpers.tri_mesh, binary=False)
    assert 'version="0.1"' in p.read_text()


@pytest.mark.parametrize("reader", [_core_vtu_read, _vtu.read], ids=["core", "python"])
@pytest.mark.parametrize("version", ["0.1", "1.0", "2.2", "2.3"])
def test_an_older_file_has_its_hexahedron_edges_swapped(reader, version, tmp_path):
    p = tmp_path / "hex.vtu"
    meshioplusplus.vtu.write(p, _hex27("VTK_LAGRANGE_HEXAHEDRON"), binary=False)
    p.write_text(p.read_text().replace('version="2.2"', f'version="{version}"'))
    want = np.arange(27)
    if version in ("0.1", "1.0"):
        want[[18, 19]] = want[[19, 18]]
    np.testing.assert_array_equal(reader(p).cells[0].data, [want])

import copy
import pathlib

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus.ensight import _ensight

from . import helpers

test_set = [
    helpers.line_mesh,
    helpers.tri_mesh,
    helpers.tri_mesh_2d,
    helpers.quad_mesh,
    helpers.tri_quad_mesh,
    helpers.tet_mesh,
    helpers.tet10_mesh,
    helpers.hex_mesh,
    helpers.hex20_mesh,
    helpers.pyramid_mesh,
    helpers.wedge_mesh,
    helpers.wedge15_mesh,
]


@pytest.mark.parametrize("mesh", test_set)
@pytest.mark.parametrize("binary", [False, True])
def test(mesh, binary, tmp_path):
    def writer(filename, mesh):
        meshioplusplus.ensight.write(filename, mesh, binary=binary)

    helpers.write_read(
        tmp_path,
        writer,
        meshioplusplus.ensight.read,
        mesh,
        1.0e-6 if binary else 1.0e-5,
        extension=".case",
    )


@pytest.mark.parametrize("mesh", [helpers.tet10_mesh, helpers.wedge15_mesh])
@pytest.mark.parametrize("binary", [False, True])
def test_cross_compat(mesh, binary, tmp_path):
    # shim (C++ where available) write -> pure-Python read
    p = tmp_path / "cross.case"
    meshioplusplus.ensight.write(p, copy.deepcopy(mesh), binary=binary)
    out = _ensight.read(p)
    assert np.allclose(mesh.points, out.points, atol=1.0e-5, rtol=0.0)
    assert np.array_equal(mesh.cells[0].data, out.cells[0].data)

    # pure-Python write -> shim (C++ where available) read
    _ensight.write(p, copy.deepcopy(mesh), binary=binary)
    out = meshioplusplus.ensight.read(p)
    assert np.allclose(mesh.points, out.points, atol=1.0e-5, rtol=0.0)
    assert np.array_equal(mesh.cells[0].data, out.cells[0].data)


def test_geo_path_write(tmp_path):
    # Writing to the .geo sibling produces both files and both are readable.
    p = tmp_path / "mesh.geo"
    meshioplusplus.ensight.write(p, copy.deepcopy(helpers.tet_mesh), binary=False)
    assert (tmp_path / "mesh.case").is_file()
    assert (tmp_path / "mesh.geo").is_file()
    for path in (tmp_path / "mesh.case", tmp_path / "mesh.geo"):
        out = meshioplusplus.ensight.read(path)
        assert np.allclose(helpers.tet_mesh.points, out.points, atol=1.0e-5, rtol=0.0)
        assert np.array_equal(helpers.tet_mesh.cells[0].data, out.cells[0].data)


@pytest.mark.parametrize("reader", [meshioplusplus.ensight.read, _ensight.read])
def test_reference_file(reader):
    # Multi-part, node/element ids given (skipped: Gold connectivity is
    # positional), tetra4 + tria3 + nsided, "ensight:part" tagging.
    this_dir = pathlib.Path(__file__).resolve().parent
    mesh = reader(this_dir / "meshes" / "ensight" / "simple.case")

    assert mesh.points.shape == (9, 3)
    ref_points = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [2.0, 0.5, 0.0],
        ]
    )
    assert np.allclose(mesh.points, ref_points, atol=1.0e-12, rtol=0.0)

    assert [c.type for c in mesh.cells] == ["tetra", "triangle", "polygon"]
    assert np.array_equal(mesh.cells[0].data, [[0, 1, 2, 3]])
    assert np.array_equal(mesh.cells[1].data, [[4, 5, 6], [4, 6, 7]])
    assert np.array_equal(np.asarray(mesh.cells[2].data), [[4, 5, 8, 6]])

    part = mesh.cell_data["ensight:part"]
    assert np.array_equal(part[0], [1])
    assert np.array_equal(part[1], [2, 2])
    assert np.array_equal(part[2], [2])


@pytest.mark.parametrize("reader", [meshioplusplus.ensight.read, _ensight.read])
def test_byteswapped_binary(reader, tmp_path):
    # A foreign-endian C-binary geometry must be detected and byte-swapped.
    swapped = ">" if np.little_endian else "<"

    def str80(s):
        return s.encode().ljust(80, b"\0")

    points = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    )
    buf = b"".join(
        [
            str80("C Binary"),
            str80("byteswapped test"),
            str80("written for meshio++"),
            str80("node id assign"),
            str80("element id assign"),
            str80("part"),
            np.array([1], dtype=swapped + "i4").tobytes(),
            str80("Mesh"),
            str80("coordinates"),
            np.array([4], dtype=swapped + "i4").tobytes(),
            np.ascontiguousarray(points.T, dtype=swapped + "f4").tobytes(),
            str80("tetra4"),
            np.array([1], dtype=swapped + "i4").tobytes(),
            np.array([1, 2, 3, 4], dtype=swapped + "i4").tobytes(),
        ]
    )
    p = tmp_path / "swapped.geo"
    p.write_bytes(buf)

    mesh = reader(p)
    assert np.allclose(mesh.points, points, atol=1.0e-6, rtol=0.0)
    assert np.array_equal(mesh.cells[0].data, [[0, 1, 2, 3]])


def _core_write(path, mesh, binary=True, fortran=False):
    # Straight through the C++ core: the shim would fall back to the Python
    # twin if the core declined, which would hide a broken native writer.
    meshioplusplus._core.ensight_write(str(path), mesh, binary, fortran)


_ENGINES = {"core": _core_write, "python": _ensight.write}


def _ragged_mesh():
    """Triangles, a jagged ``polygon`` block (3-, 5- and 4-gons) and an
    ``nfaced`` block mixing node counts (tetra, pyramid, tetra, hexahedron)."""
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
            [2.0, 0.5, 0.0],
            [2.0, 0.5, 1.0],
        ]
    )
    polygons = [[0, 1, 2], [1, 8, 2, 6, 5], [0, 3, 7, 4]]
    tet = [[0, 2, 1], [0, 1, 4], [1, 2, 4], [0, 4, 2]]
    pyramid = [[0, 1, 2, 3], [0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]]
    hexa = [
        [0, 3, 2, 1],
        [4, 5, 6, 7],
        [0, 1, 5, 4],
        [2, 3, 7, 6],
        [0, 4, 7, 3],
        [1, 2, 6, 5],
    ]
    cells = [
        ("triangle", np.array([[0, 1, 2], [0, 2, 3]])),
        ("polygon", [np.array(p) for p in polygons]),
        (
            "polyhedron",
            [[np.array(f) for f in cell] for cell in (tet, pyramid, tet, hexa)],
        ),
    ]
    return meshioplusplus.Mesh(points, cells)


def _polygon_rows(block):
    return [list(map(int, row)) for row in block.data]


def _polyhedron_cells(block):
    return [[list(map(int, face)) for face in cell] for cell in block.data]


@pytest.mark.parametrize("engine", sorted(_ENGINES))
@pytest.mark.parametrize("binary", [False, True])
def test_nsided_nfaced_round_trip(engine, binary, tmp_path):
    mesh = _ragged_mesh()
    case = tmp_path / "ragged.case"
    _ENGINES[engine](case, copy.deepcopy(mesh), binary=binary)
    for reader in (meshioplusplus.ensight.read, _ensight.read):
        back = reader(case)
        np.testing.assert_allclose(back.points, mesh.points, atol=1.0e-5)
        # nfaced cells come back grouped by node count, in first-seen order.
        assert [c.type for c in back.cells] == [
            "triangle",
            "polygon",
            "polyhedron4",
            "polyhedron5",
            "polyhedron8",
        ]
        assert np.array_equal(back.cells[0].data, mesh.cells[0].data)
        assert _polygon_rows(back.cells[1]) == _polygon_rows(mesh.cells[1])
        tet, pyramid, tet2, hexa = _polyhedron_cells(mesh.cells[2])
        assert _polyhedron_cells(back.cells[2]) == [tet, tet2]
        assert _polyhedron_cells(back.cells[3]) == [pyramid]
        assert _polyhedron_cells(back.cells[4]) == [hexa]


@pytest.mark.parametrize("engine", sorted(_ENGINES))
def test_nsided_uniform_polygon_block_round_trips(engine, tmp_path):
    # A block whose polygons all have one node count is a rectangular array.
    mesh = copy.deepcopy(helpers.polygon_mesh)
    assert isinstance(mesh.cells[3].data, np.ndarray)
    case = tmp_path / "poly.case"
    _ENGINES[engine](case, mesh, binary=False)
    back = _ensight.read(case)
    assert [c.type for c in back.cells] == ["triangle", "quad", "polygon", "polygon"]
    assert _polygon_rows(back.cells[2]) == _polygon_rows(mesh.cells[2])
    assert _polygon_rows(back.cells[3]) == _polygon_rows(mesh.cells[3])


def _byte_parity_meshes():
    mixed = _ragged_mesh()
    empty = copy.deepcopy(helpers.tri_mesh)
    empty.cells.append(meshioplusplus.CellBlock("polygon", []))
    empty.cells.append(meshioplusplus.CellBlock("polyhedron", []))
    return {
        "mixed": mixed,
        "polygon_mesh": helpers.polygon_mesh,
        "polygon_mesh_one_cell": helpers.polygon_mesh_one_cell,
        "polygon2_mesh": helpers.polygon2_mesh,
        "polyhedron_mesh": helpers.polyhedron_mesh,
        "empty_ragged_blocks": empty,
    }


@pytest.mark.parametrize("name", sorted(_byte_parity_meshes()))
@pytest.mark.parametrize(
    "binary, fortran", [(False, False), (True, False), (True, True)]
)
def test_engines_write_the_same_nsided_nfaced_bytes(name, binary, fortran, tmp_path):
    mesh = _byte_parity_meshes()[name]
    a, b = tmp_path / "a.case", tmp_path / "b.case"
    _core_write(a, copy.deepcopy(mesh), binary=binary, fortran=fortran)
    _ensight.write(b, copy.deepcopy(mesh), binary=binary, fortran=fortran)
    assert (tmp_path / "a.geo").read_bytes() == (tmp_path / "b.geo").read_bytes()


def test_shim_writes_ragged_blocks_through_the_core(tmp_path):
    # The public writer reaches the native nsided/nfaced writer (no silent
    # fall-back to the Python twin): under strict mode a decline would raise.
    from meshioplusplus import _fallback

    _fallback.set_strict_core(True)
    try:
        case = tmp_path / "shim.case"
        meshioplusplus.ensight.write(case, _ragged_mesh(), binary=True)
    finally:
        _fallback.set_strict_core(None)
    assert [c.type for c in _ensight.read(case).cells][:2] == ["triangle", "polygon"]


def test_wildcard_model_raises(tmp_path):
    p = tmp_path / "transient.case"
    p.write_text(
        "FORMAT\ntype: ensight gold\n\nGEOMETRY\nmodel: 1 transient.****.geo\n"
    )
    with pytest.raises(meshioplusplus.ReadError):
        _ensight.read(p)


def test_read_metadata_and_time_step_reads_a_transient_case(tmp_path):
    """roadmap §1 tier B1: read_ensight_metadata reports the .case file's
    TIME time_values, and time_step selects one step's VARIABLE files -- one
    file per step (EnSight's own convention), so the fixture is two
    per-node scalar files, written by hand (plain text, no external tool)."""
    mesh = copy.deepcopy(helpers.tri_mesh)
    path = tmp_path / "transient.case"
    meshioplusplus.ensight.write(path, mesh, binary=False)

    with open(path, "a") as f:
        f.write(
            "TIME\n"
            "time set:              1\n"
            "number of steps:       2\n"
            "filename start number: 0\n"
            "filename increment:    1\n"
            "time values:\n"
            "0.0\n"
            "2.5\n"
            "VARIABLE\n"
            "scalar per node:    1  pressure  pressure.****.scl\n"
        )

    n = len(mesh.points)
    for step, base in enumerate([10.0, 11.0]):
        vpath = tmp_path / f"pressure.{step:04d}.scl"
        lines = ["pressure", "part", "         1", "coordinates"]
        lines += [str(base + i * 10.0) for i in range(n)]
        vpath.write_text("\n".join(lines) + "\n")

    meta = meshioplusplus.read_metadata(path, "ensight")
    assert meta["time_values"] == [0.0, 2.5]
    assert meta["fell_back_to_full_read"] is True

    mesh0 = meshioplusplus.ensight.read(path, time_step=0)
    assert mesh0.point_data["pressure"][0] == 10.0
    mesh1 = meshioplusplus.ensight.read(path, time_step=1)
    assert mesh1.point_data["pressure"][0] == 11.0
    mesh_last = meshioplusplus.ensight.read(path, time_step=-1)
    assert mesh_last.point_data["pressure"][0] == 11.0

    with pytest.raises(meshioplusplus.ReadError):
        meshioplusplus.ensight.read(path, time_step=5)


@pytest.mark.parametrize("binary", [False, True])
def test_variable_write_round_trips_scalar_vector_tensor_and_constant(binary, tmp_path):
    """roadmap §1.1: write_ensight now writes a VARIABLE section -- scalar,
    vector (2-component padded to 3) and tensor symm point_data, tensor
    asym cell_data, and a field_data scalar as 'constant per case'."""
    mesh = copy.deepcopy(helpers.tri_mesh)
    mesh.point_data["temp"] = np.array([1.0, 2.0, 3.0, 4.0])
    mesh.point_data["vel2d"] = np.array(
        [[1.0, 10.0], [2.0, 20.0], [3.0, 30.0], [4.0, 40.0]]
    )
    mesh.point_data["sig"] = np.arange(24, dtype="f8").reshape(4, 6)
    mesh.cell_data["eps"] = [np.arange(18, dtype="f8").reshape(2, 9)]
    mesh.field_data["gravity"] = np.array([9.81])

    path = tmp_path / "vars.case"
    meshioplusplus.ensight.write(path, mesh, binary=binary)
    out = meshioplusplus.ensight.read(path)

    atol = 1.0e-6 if binary else 1.0e-5
    assert np.allclose(out.point_data["temp"], mesh.point_data["temp"], atol=atol)
    assert out.point_data["vel2d"].shape == (4, 3)
    assert np.allclose(
        out.point_data["vel2d"][:, :2], mesh.point_data["vel2d"], atol=atol
    )
    assert np.allclose(out.point_data["vel2d"][:, 2], 0.0, atol=atol)
    assert np.allclose(out.point_data["sig"], mesh.point_data["sig"], atol=atol)
    assert np.allclose(out.cell_data["eps"][0], mesh.cell_data["eps"][0], atol=atol)
    assert out.field_data["gravity"][0] == pytest.approx(9.81, abs=atol)


def test_variable_write_skips_unsupported_component_counts(tmp_path):
    mesh = copy.deepcopy(helpers.tri_mesh)
    mesh.point_data["weird"] = np.arange(20, dtype="f8").reshape(4, 5)
    mesh.point_data["ok"] = np.array([1.0, 2.0, 3.0, 4.0])

    path = tmp_path / "skip.case"
    meshioplusplus.ensight.write(path, mesh, binary=False)
    out = meshioplusplus.ensight.read(path)

    assert "weird" not in out.point_data
    assert np.allclose(out.point_data["ok"], mesh.point_data["ok"])


# --------------------------------------------------------------------------- #
# Fortran binary, and binary variable files other readers can read            #
# --------------------------------------------------------------------------- #


def _data_mesh():
    pts = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]], dtype=float)
    return meshioplusplus.Mesh(
        pts,
        [("tetra", np.array([[0, 1, 2, 3]])), ("triangle", np.array([[1, 2, 4]]))],
        point_data={"T": np.arange(1.0, 6.0), "V": np.arange(15.0).reshape(5, 3)},
        cell_data={"S": [np.array([7.0]), np.array([8.0])]},
    )


def _reframe(path, width, order):
    """Rewrite a little-endian Fortran-binary file with `width`-byte markers in
    byte order `order`, swapping the numeric records to match."""
    import struct

    from meshioplusplus._fortran_records import fortran_records, sniff_fortran_records

    data = path.read_bytes()
    records = fortran_records(data, sniff_fortran_records(data), "test")
    fmt = order + ("i" if width == 4 else "q")
    out = b""
    for off, size in records:
        payload = data[off : off + size]
        if size != 80 and order == ">":
            payload = np.frombuffer(payload, "<u4").byteswap().tobytes()
        out += struct.pack(fmt, size) + payload + struct.pack(fmt, size)
    path.write_bytes(out)


@pytest.mark.parametrize("engine", ["core", "python"])
def test_fortran_binary_round_trips(engine, tmp_path):
    mesh = _data_mesh()
    case = tmp_path / "m.case"
    if engine == "core":
        meshioplusplus._core.ensight_write(str(case), mesh, True, True)
    else:
        _ensight.write(case, mesh, binary=True, fortran=True)
    geo = (tmp_path / "m.geo").read_bytes()
    assert geo[:4] == np.int32(80).tobytes() and geo[4:18] == b"Fortran Binary"
    for back in (meshioplusplus.read(case), _ensight.read(case)):
        np.testing.assert_allclose(back.points, mesh.points)
        assert [c.type for c in back.cells] == ["tetra", "triangle"]
    if engine == "core":  # the Python writer is geometry-only
        back = meshioplusplus.read(case)
        np.testing.assert_allclose(back.point_data["V"], mesh.point_data["V"])
        np.testing.assert_allclose(back.cell_data["S"][1], [8.0])


@pytest.mark.parametrize("width, order", [(4, ">"), (8, "<"), (8, ">")])
def test_fortran_binary_reads_any_marker_width_and_byte_order(width, order, tmp_path):
    mesh = _data_mesh()
    case = tmp_path / "m.case"
    meshioplusplus._core.ensight_write(str(case), mesh, True, True)
    for f in tmp_path.iterdir():
        if f.suffix != ".case":
            _reframe(f, width, order)
    back = meshioplusplus.read(case)
    np.testing.assert_allclose(back.points, mesh.points)
    np.testing.assert_allclose(back.point_data["T"], mesh.point_data["T"])
    np.testing.assert_allclose(_ensight.read(case).points, mesh.points)


@pytest.mark.parametrize("fortran", [False, True])
def test_vtk_reads_binary_variables(fortran, tmp_path):
    # A variable file starts with its description: until v16.17.0 the writer
    # put a "C Binary" record first, and VTK/ParaView read no variables.
    vtk = pytest.importorskip("vtk")
    from vtk.util.numpy_support import vtk_to_numpy

    mesh = _data_mesh()
    case = tmp_path / "m.case"
    meshioplusplus._core.ensight_write(str(case), mesh, True, fortran)
    reader = vtk.vtkEnSightGoldBinaryReader()
    reader.SetCaseFileName(str(case))
    reader.Update()
    block = reader.GetOutput().GetBlock(0)
    np.testing.assert_allclose(vtk_to_numpy(block.GetPoints().GetData()), mesh.points)
    np.testing.assert_allclose(
        vtk_to_numpy(block.GetPointData().GetArray("T")), mesh.point_data["T"]
    )
    np.testing.assert_allclose(
        vtk_to_numpy(block.GetPointData().GetArray("V")), mesh.point_data["V"]
    )
    np.testing.assert_allclose(
        vtk_to_numpy(block.GetCellData().GetArray("S")), [7.0, 8.0]
    )


def test_variable_files_written_before_v16_17_still_read(tmp_path):
    mesh = _data_mesh()
    case = tmp_path / "m.case"
    meshioplusplus._core.ensight_write(str(case), mesh, True, False)
    for name in ("T.scl", "V.vec", "S.escl"):
        f = tmp_path / name
        f.write_bytes(b"C Binary".ljust(80, b"\0") + f.read_bytes())
    back = meshioplusplus.read(case)
    np.testing.assert_allclose(back.point_data["V"], mesh.point_data["V"])
    np.testing.assert_allclose(back.cell_data["S"][0], [7.0])


def test_fortran_needs_binary(tmp_path):
    with pytest.raises(meshioplusplus.WriteError):
        _ensight.write(tmp_path / "m.case", _data_mesh(), binary=False, fortran=True)
    with pytest.raises(meshioplusplus.WriteError):
        meshioplusplus.ensight.write(tmp_path / "m.case", _data_mesh(), False, True)


# --------------------------------------------------------------------------- #
# Roadmap 3.4.5: the Python twin carries VARIABLE reading and writing         #
# --------------------------------------------------------------------------- #


def _variable_mesh():
    mesh = copy.deepcopy(_data_mesh())
    mesh.point_data["sig"] = np.arange(30, dtype="f8").reshape(5, 6)
    mesh.point_data["vel2d"] = np.arange(10, dtype="f8").reshape(5, 2)
    mesh.cell_data["eps"] = [np.arange(9, dtype="f8").reshape(1, 9)] * 2
    mesh.field_data["gravity"] = np.array([9.81])
    return mesh


@pytest.mark.parametrize("mode", ["ascii", "binary", "fortran"])
def test_engines_write_the_same_variable_bytes(mode, tmp_path):
    from meshioplusplus import _core

    mesh = _variable_mesh()
    binary, fortran = mode != "ascii", mode == "fortran"
    a, b = tmp_path / "core", tmp_path / "twin"
    a.mkdir()
    b.mkdir()
    _core.ensight_write(str(a / "m.case"), mesh, binary, fortran)
    _ensight.write(b / "m.case", mesh, binary=binary, fortran=fortran)
    names = sorted(p.name for p in a.iterdir())
    assert names == sorted(p.name for p in b.iterdir())
    assert "sig.tsym" in names and "eps.etasym" in names
    for name in names:
        assert (a / name).read_bytes() == (b / name).read_bytes(), name


@pytest.mark.parametrize("mode", ["ascii", "binary", "fortran"])
def test_engines_read_each_others_variables(mode, tmp_path):
    from meshioplusplus import _core

    mesh = _variable_mesh()
    binary, fortran = mode != "ascii", mode == "fortran"
    path = tmp_path / "m.case"
    _core.ensight_write(str(path), mesh, binary, fortran)
    native, twin = _core.ensight_read(str(path), 0), _ensight.read(path)
    assert sorted(native.point_data) == sorted(twin.point_data)
    for name in native.point_data:
        np.testing.assert_array_equal(native.point_data[name], twin.point_data[name])
    for name in native.cell_data:
        for x, y in zip(native.cell_data[name], twin.cell_data[name]):
            np.testing.assert_array_equal(x, y)
    assert native.field_data["gravity"][0] == twin.field_data["gravity"][0]


def test_python_twin_resolves_time_steps(tmp_path):
    mesh = copy.deepcopy(helpers.tri_mesh)
    path = tmp_path / "transient.case"
    _ensight.write(path, mesh, binary=False)
    with open(path, "a") as f:
        f.write(
            "TIME\ntime set: 1\nnumber of steps: 2\nfilename start number: 0\n"
            "filename increment: 1\ntime values:\n0.0\n2.5\n"
            "VARIABLE\nscalar per node: 1 pressure pressure.****.scl\n"
        )
    n = len(mesh.points)
    for step, base in enumerate([10.0, 11.0]):
        lines = ["pressure", "part", "         1", "coordinates"]
        lines += [str(base + i) for i in range(n)]
        (tmp_path / f"pressure.{step:04d}.scl").write_text("\n".join(lines) + "\n")
    assert _ensight.read(path, 0).point_data["pressure"][0] == 10.0
    assert _ensight.read(path, -1).point_data["pressure"][0] == 11.0
    with pytest.raises(meshioplusplus.ReadError):
        _ensight.read(path, 5)

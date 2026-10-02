"""Independent multi-piece XML and legacy structured VTK fixtures."""

import base64
import importlib
import xml.etree.ElementTree as ET
import zlib

import numpy as np
import pytest

import meshioplusplus as mio
from meshioplusplus import _fallback

TYPES = {
    "vtp": "PolyData",
    "vts": "StructuredGrid",
    "vtr": "RectilinearGrid",
    "vti": "ImageData",
}


def xml_fixture(fmt, missing=False, incompatible=False, local_regions=True):
    root = ET.Element("VTKFile", type=TYPES[fmt], byte_order="LittleEndian")
    attrs = {} if fmt == "vtp" else {"WholeExtent": "2 4 0 1 0 1"}
    grid = ET.SubElement(root, TYPES[fmt], attrs)

    def array(parent, name, data, dtype="Int64", components=None):
        attrs = {
            "Name": name,
            "type": dtype,
            "format": "ascii",
            "NumberOfTuples": str(np.asarray(data).size // (components or 1)),
        }
        if components:
            attrs["NumberOfComponents"] = str(components)
        da = ET.SubElement(parent, "DataArray", attrs)
        da.text = " ".join(str(x) for x in np.asarray(data).reshape(-1))
        return da

    array(ET.SubElement(grid, "FieldData"), "units", [42])
    expected = []
    for part in range(2):
        points = np.array(
            [
                [x, y, z]
                for z in range(2 if fmt != "vtp" else 1)
                for y in range(2)
                for x in (2 + part, 3 + part)
            ],
            dtype=float,
        )
        expected.append(points)
        attrs = (
            {"NumberOfPoints": "4", "NumberOfPolys": "1"}
            if fmt == "vtp"
            else {"Extent": f"{2 + part} {3 + part} 0 1 0 1"}
        )
        piece = ET.SubElement(grid, "Piece", attrs)
        if fmt in ("vts", "vtp"):
            array(ET.SubElement(piece, "Points"), "Points", points, "Float64", 3)
        if fmt == "vtp":
            polys = ET.SubElement(piece, "Polys")
            array(polys, "connectivity", [0, 1, 3, 2])
            array(polys, "offsets", [4])
        if fmt == "vtr":
            coordinates = ET.SubElement(piece, "Coordinates")
            for axis, values in enumerate(([2 + part, 3 + part], [0, 1], [0, 1])):
                array(coordinates, f"axis_{axis}", values, "Float64")
        pd = ET.SubElement(piece, "PointData")
        array(pd, "tag", np.arange(len(points)) + 10 * part)
        if not missing or part == 0:
            array(
                pd,
                "optional",
                np.ones(len(points)),
                "Float32" if incompatible and part else "Float64",
            )
        cd = ET.SubElement(piece, "CellData")
        array(cd, "material", [20 + part])
        array(cd, "vtkGhostType", [part], "UInt8")
        if local_regions:
            fd = ET.SubElement(piece, "FieldData")
            array(fd, "region:point:corner", [0])
            array(fd, "region:cell:body", [0])
            array(fd, "region:side:wall", [[0, 1]], components=2)
            array(fd, "region:point:empty", [])
    return ET.tostring(root), np.concatenate(expected)


@pytest.mark.parametrize("fmt", TYPES)
@pytest.mark.parametrize(
    "missing,incompatible", [(False, False), (True, False), (False, True)]
)
def test_xml_piece_assembly(tmp_path, monkeypatch, fmt, missing, incompatible):
    data, points = xml_fixture(fmt, missing, incompatible)
    path = tmp_path / f"pieces.{fmt}"
    path.write_bytes(data)
    core = pytest.importorskip("meshioplusplus._core")
    readers = [
        getattr(core, fmt + "_read"),
        importlib.import_module(f"meshioplusplus.{fmt}._{fmt}").read,
    ]
    for read in readers:
        mesh = read(str(path))
        np.testing.assert_array_equal(mesh.points, points)
        assert sum(len(cb.data) for cb in mesh.cells) == 2
        np.testing.assert_array_equal(
            np.concatenate(mesh.cell_data["material"]), [20, 21]
        )
        np.testing.assert_array_equal(
            np.concatenate(mesh.cell_data["vtkGhostType"]), [0, 1]
        )
        assert ("optional" in mesh.point_data) == (not missing and not incompatible)
        regions = {(r.kind, r.name): r.entries for r in mesh.regions}
        np.testing.assert_array_equal(regions["point", "corner"], [0, len(points) // 2])
        np.testing.assert_array_equal(regions["cell", "body"], [0, 1])
        np.testing.assert_array_equal(regions["side", "wall"], [[0, 1], [1, 1]])
        assert regions["point", "empty"].size == 0
        np.testing.assert_array_equal(mesh.field_data["units"], [42])
        assert mesh.cells[0].data[-1].min() >= len(points) // 2
    monkeypatch.setattr(_fallback, "_strict", True)
    assert len(mio.read(path).points) == len(points)
    metadata = mio.read_metadata(path)
    assert metadata["num_points"] == len(points)
    assert ("optional" in metadata["point_data_names"]) == (
        not missing and not incompatible
    )
    bare = mio.read(path, points_only=True)
    assert bare.point_data == {} and bare.cell_data == {} and bare.field_data == {}
    assert len(bare.regions) == 4
    selected = mio.read(path, arrays=["material"])
    assert selected.point_data == {} and selected.field_data == {}
    assert list(selected.cell_data) == ["material"]


def legacy_fixture(dataset, dims, binary=False, version="4.2"):
    parts = [
        f"# vtk DataFile Version {version}\nindependent structured fixture\n{'BINARY' if binary else 'ASCII'}\nDATASET {dataset}\nDIMENSIONS {' '.join(map(str, dims))}\n".encode()
    ]

    def values(data, dtype="f8"):
        arr = np.asarray(data, dtype=(">" if binary else "=") + dtype)
        parts.append(
            arr.tobytes()
            if binary
            else " ".join(str(x) for x in arr.reshape(-1)).encode()
        )
        parts.append(b"\n")

    axes = [np.arange(n, dtype=float) * (k + 1) + k for k, n in enumerate(dims)]
    points = np.array([[x, y, z] for z in axes[2] for y in axes[1] for x in axes[0]])
    if dataset == "STRUCTURED_POINTS":
        parts.append(b"ORIGIN 0 1 2\nASPECT_RATIO 1 2 3\n")
    elif dataset == "STRUCTURED_GRID":
        parts.append(f"POINTS {len(points)} double\n".encode())
        values(points)
    else:
        for name, axis in zip("XYZ", axes):
            parts.append(f"{name}_COORDINATES {len(axis)} double\n".encode())
            values(axis)
    nc = int(np.prod([max(n - 1, 1) for n in dims]))
    parts.append(
        f"POINT_DATA {len(points)}\nSCALARS tag int 1\nLOOKUP_TABLE default\n".encode()
    )
    values(np.arange(len(points)), "i4")
    parts.append(b"VECTORS velocity double\n")
    values(np.arange(len(points) * 3).reshape(-1, 3))
    parts.append(f"CELL_DATA {nc}\nTENSORS stress double\n".encode())
    values(np.arange(nc * 9).reshape(-1, 3, 3))
    return b"".join(parts), points


@pytest.mark.parametrize(
    "dataset", ["STRUCTURED_POINTS", "STRUCTURED_GRID", "RECTILINEAR_GRID"]
)
@pytest.mark.parametrize(
    "dims", [(3, 2, 2), (3, 2, 1), (1, 3, 2), (1, 1, 4), (1, 1, 1)]
)
@pytest.mark.parametrize("binary", [False, True])
@pytest.mark.parametrize("version", ["4.2", "5.1"])
def test_legacy_structured(tmp_path, monkeypatch, dataset, dims, binary, version):
    raw, points = legacy_fixture(dataset, dims, binary, version)
    path = tmp_path / "structured.vtk"
    path.write_bytes(raw)
    core = pytest.importorskip("meshioplusplus._core")
    reference = importlib.import_module("meshioplusplus.vtk._main").read(path)
    native = core.vtk_read(str(path))
    np.testing.assert_array_equal(native.points, points)
    np.testing.assert_array_equal(native.points, reference.points)
    assert [cb.type for cb in native.cells] == [cb.type for cb in reference.cells]
    for a, b in zip(native.cells, reference.cells):
        np.testing.assert_array_equal(a.data, b.data)
    for name in ("tag", "velocity"):
        np.testing.assert_array_equal(
            native.point_data[name], reference.point_data[name]
        )
    np.testing.assert_array_equal(
        native.cell_data["stress"][0], reference.cell_data["stress"][0]
    )
    monkeypatch.setattr(_fallback, "_strict", True)
    assert len(mio.read(path).points) == len(points)


def encode_xml(data, encoding, header, order, compressed):
    root = ET.fromstring(data)
    bo = ">" if order == "BigEndian" else "<"
    hd = np.dtype(bo + ("u8" if header == "UInt64" else "u4"))
    root.set("byte_order", order)
    root.set("header_type", header)
    if compressed:
        root.set("compressor", "vtkZLibDataCompressor")
    payload = bytearray()
    dtypes = {"Int64": "i8", "UInt8": "u1", "Float64": "f8", "Float32": "f4"}
    for da in root.iter("DataArray"):
        raw = np.array(
            (da.text or "").split(), dtype=bo + dtypes[da.get("type")]
        ).tobytes()
        if compressed:
            body = zlib.compress(raw) if raw else b""
            h = np.array(
                [int(bool(raw)), max(len(raw), 1), len(raw)]
                + ([len(body)] if raw else []),
                dtype=hd,
            ).tobytes()
        else:
            h, body = np.array([len(raw)], dtype=hd).tobytes(), raw
        encoded = base64.b64encode(h) + base64.b64encode(body)
        if encoding == "binary":
            da.set("format", "binary")
            da.text = (encoded if compressed else base64.b64encode(h + body)).decode()
        else:
            da.set("format", "appended")
            da.set("offset", str(len(payload)))
            da.text = None
            payload.extend(h + body if encoding == "raw" else encoded)
    if encoding == "binary":
        return ET.tostring(root)
    prefix = ET.tostring(root).removesuffix(b"</VTKFile>")
    return (
        prefix
        + f'<AppendedData encoding="{encoding}">_'.encode()
        + payload
        + b"</AppendedData></VTKFile>"
    )


@pytest.mark.parametrize("fmt", TYPES)
@pytest.mark.parametrize("encoding", ["binary", "base64", "raw"])
@pytest.mark.parametrize(
    "header,order", [("UInt32", "LittleEndian"), ("UInt64", "BigEndian")]
)
@pytest.mark.parametrize("compressed", [False, True])
def test_xml_piece_framing(tmp_path, fmt, encoding, header, order, compressed):
    core = pytest.importorskip("meshioplusplus._core")
    if compressed and not getattr(core, "__has_zlib__", False):
        pytest.skip("native build lacks zlib")
    data, points = xml_fixture(fmt)
    path = tmp_path / f"framed.{fmt}"
    path.write_bytes(encode_xml(data, encoding, header, order, compressed))
    for read in (
        getattr(core, fmt + "_read"),
        importlib.import_module(f"meshioplusplus.{fmt}._{fmt}").read,
    ):
        mesh = read(str(path))
        np.testing.assert_array_equal(mesh.points, points)
        np.testing.assert_array_equal(
            np.concatenate(mesh.cell_data["material"]), [20, 21]
        )


@pytest.mark.parametrize("fmt", ["vts", "vtr", "vti"])
@pytest.mark.parametrize(
    "extent",
    [
        "3 2 0 1 0 1",
        "0 5 0 1 0 1",
        "2 x 0 1 0 1",
        "2 3 0 1",
        "-9223372036854775808 9223372036854775807 0 1 0 1",
    ],
)
def test_bad_piece_extents(tmp_path, fmt, extent):
    data, _ = xml_fixture(fmt)
    root = ET.fromstring(data)
    root.find(TYPES[fmt]).find("Piece").set("Extent", extent)
    path = tmp_path / f"bad.{fmt}"
    path.write_bytes(ET.tostring(root))
    core = pytest.importorskip("meshioplusplus._core")
    for read in (
        getattr(core, fmt + "_read"),
        importlib.import_module(f"meshioplusplus.{fmt}._{fmt}").read,
    ):
        with pytest.raises(mio.ReadError):
            read(str(path))


def test_vtp_mixed_sections_and_empty_piece(tmp_path):
    raw, points = xml_fixture("vtp")
    root = ET.fromstring(raw)
    grid = root.find("PolyData")
    for piece in grid.findall("Piece"):
        for section, conn, offsets in (("Verts", "0", "1"), ("Lines", "0 1", "2")):
            node = ET.SubElement(piece, section)
            for name, text in (("connectivity", conn), ("offsets", offsets)):
                ET.SubElement(
                    node, "DataArray", type="Int64", Name=name, format="ascii"
                ).text = text
        for da in piece.find("CellData"):
            da.text = " ".join([da.text] * 3)
    # Supply empty arrays too so every piece has a compatible data schema.
    empty = ET.SubElement(grid, "Piece", NumberOfPoints="0")
    for section in ("PointData", "CellData"):
        holder = ET.SubElement(empty, section)
        for da in grid.find("Piece").find(section):
            ET.SubElement(holder, "DataArray", da.attrib)
    path = tmp_path / "mixed.vtp"
    path.write_bytes(ET.tostring(root))
    core = pytest.importorskip("meshioplusplus._core")
    reference = importlib.import_module("meshioplusplus.vtp._vtp").read(path)
    mesh = core.vtp_read(str(path))
    np.testing.assert_array_equal(mesh.points, points)
    assert [cb.type for cb in mesh.cells] == ["vertex", "line", "quad"] * 2
    for a, b in zip(mesh.cells, reference.cells):
        np.testing.assert_array_equal(a.data, b.data)
    np.testing.assert_array_equal(
        np.concatenate(mesh.cell_data["material"]), [20] * 3 + [21] * 3
    )
    regions = {(r.kind, r.name): r.entries for r in mesh.regions}
    np.testing.assert_array_equal(regions["cell", "body"], [0, 3])
    np.testing.assert_array_equal(regions["side", "wall"], [[0, 1], [3, 1]])


@pytest.mark.parametrize(
    "dataset", ["STRUCTURED_POINTS", "STRUCTURED_GRID", "RECTILINEAR_GRID"]
)
@pytest.mark.parametrize("mutation", ["dimensions", "points", "data", "truncate"])
def test_malformed_legacy_structured(tmp_path, dataset, mutation):
    raw, _ = legacy_fixture(dataset, (3, 2, 2), binary=True)
    if mutation == "dimensions":
        raw = raw.replace(b"DIMENSIONS 3 2 2", b"DIMENSIONS -3 2 2")
    elif mutation == "points":
        if dataset == "STRUCTURED_GRID":
            raw = raw.replace(b"POINTS 12 double", b"POINTS 13 double")
        elif dataset == "RECTILINEAR_GRID":
            raw = raw.replace(b"X_COORDINATES 3 double", b"X_COORDINATES 4 double")
        else:
            raw = raw.replace(b"ASPECT_RATIO 1 2 3", b"ASPECT_RATIO 1 2")
    elif mutation == "data":
        raw = raw.replace(b"POINT_DATA 12", b"POINT_DATA 13")
    else:
        raw = raw[:-8]
    path = tmp_path / "bad.vtk"
    path.write_bytes(raw)
    core = pytest.importorskip("meshioplusplus._core")
    with pytest.raises(mio.ReadError):
        core.vtk_read(str(path))


@pytest.mark.parametrize("fmt", ["vts", "vtr", "vti"])
@pytest.mark.parametrize("all_empty", [False, True])
def test_empty_structured_pieces(tmp_path, fmt, all_empty):
    raw, points = xml_fixture(fmt, local_regions=False)
    root = ET.fromstring(raw)
    grid = root.find(TYPES[fmt])
    if all_empty:
        for piece in grid.findall("Piece"):
            grid.remove(piece)
        grid.set("WholeExtent", "0 -1 0 -1 0 -1")
    empty = ET.SubElement(grid, "Piece", Extent="0 -1 0 -1 0 -1")
    if not all_empty:
        for section in ("PointData", "CellData"):
            holder = ET.SubElement(empty, section)
            for da in grid.find("Piece").find(section):
                ET.SubElement(holder, "DataArray", da.attrib)
    path = tmp_path / f"empty.{fmt}"
    path.write_bytes(ET.tostring(root))
    core = pytest.importorskip("meshioplusplus._core")
    for read in (
        getattr(core, fmt + "_read"),
        importlib.import_module(f"meshioplusplus.{fmt}._{fmt}").read,
    ):
        mesh = read(str(path))
        assert len(mesh.points) == (0 if all_empty else len(points))
        assert sum(len(cb.data) for cb in mesh.cells) == (0 if all_empty else 2)
        if not all_empty:
            assert "tag" in mesh.point_data and "material" in mesh.cell_data


@pytest.mark.parametrize("fmt", TYPES)
@pytest.mark.parametrize("encoding", ["binary", "raw", "base64"])
def test_vtk_produced_multiple_pieces(tmp_path, fmt, encoding):
    vtk = pytest.importorskip("vtk")
    from vtk.util.numpy_support import numpy_to_vtk

    image = vtk.vtkImageData()
    image.SetDimensions(3, 2, 2)
    if fmt == "vti":
        source = image
    elif fmt == "vtr":
        source = vtk.vtkRectilinearGrid()
        source.SetDimensions(3, 2, 2)
        source.SetXCoordinates(
            numpy_to_vtk(np.array([0.0, 0.25, 1.0], dtype=np.float64), deep=True)
        )
        source.SetYCoordinates(numpy_to_vtk(np.array([0.0, 1.0]), deep=True))
        source.SetZCoordinates(numpy_to_vtk(np.array([0.0, 1.0]), deep=True))
    else:
        if fmt == "vts":
            source = vtk.vtkStructuredGrid()
            source.SetDimensions(3, 2, 2)
            points = vtk.vtkPoints()
            points.SetData(
                numpy_to_vtk(
                    np.array(
                        [image.GetPoint(i) for i in range(image.GetNumberOfPoints())]
                    ),
                    deep=True,
                )
            )
            source.SetPoints(points)
        else:
            surface = vtk.vtkDataSetSurfaceFilter()
            surface.SetInputData(image)
            surface.Update()
            source = surface.GetOutput()
    array = numpy_to_vtk(
        np.arange(source.GetNumberOfPoints(), dtype=np.int32), deep=True
    )
    array.SetName("tag")
    source.GetPointData().AddArray(array)
    writers = {
        "vtp": vtk.vtkXMLPolyDataWriter,
        "vts": vtk.vtkXMLStructuredGridWriter,
        "vtr": vtk.vtkXMLRectilinearGridWriter,
        "vti": vtk.vtkXMLImageDataWriter,
    }
    writer = writers[fmt]()
    writer.SetInputData(source)
    writer.SetNumberOfPieces(2)
    writer.SetHeaderTypeToUInt64()
    writer.SetCompressorTypeToZLib()
    if encoding == "binary":
        writer.SetDataModeToBinary()
    else:
        writer.SetDataModeToAppended()
        writer.SetEncodeAppendedData(encoding == "base64")
    path = tmp_path / f"producer.{fmt}"
    writer.SetFileName(str(path))
    assert writer.Write() == 1
    core = pytest.importorskip("meshioplusplus._core")
    native = getattr(core, fmt + "_read")(str(path))
    reference = importlib.import_module(f"meshioplusplus.{fmt}._{fmt}").read(path)
    np.testing.assert_array_equal(native.points, reference.points)
    np.testing.assert_array_equal(native.point_data["tag"], reference.point_data["tag"])
    # SetInputData has no streaming source, so VTK emits the complete input twice.
    assert len(native.points) == 2 * source.GetNumberOfPoints()
    assert sum(len(cb.data) for cb in native.cells) == 2 * source.GetNumberOfCells()


@pytest.mark.parametrize("fmt", TYPES)
def test_duplicate_array_cannot_replace_missing_piece_data(tmp_path, fmt):
    raw, _ = xml_fixture(fmt, missing=True)
    root = ET.fromstring(raw)
    pd = root.find(TYPES[fmt]).find("Piece").find("PointData")
    duplicate = ET.fromstring(ET.tostring(pd.find("DataArray[@Name='optional']")))
    pd.append(duplicate)
    path = tmp_path / f"duplicate.{fmt}"
    path.write_bytes(ET.tostring(root))
    core = pytest.importorskip("meshioplusplus._core")
    for read in (
        getattr(core, fmt + "_read"),
        importlib.import_module(f"meshioplusplus.{fmt}._{fmt}").read,
    ):
        with pytest.raises(mio.ReadError, match="duplicate"):
            read(str(path))


@pytest.mark.parametrize(
    "attribute,value",
    [("NumberOfPoints", "-4"), ("NumberOfPolys", "2"), ("NumberOfStrips", "1")],
)
def test_bad_polydata_counts(tmp_path, attribute, value):
    raw, _ = xml_fixture("vtp")
    root = ET.fromstring(raw)
    root.find("PolyData").find("Piece").set(attribute, value)
    path = tmp_path / "bad-count.vtp"
    path.write_bytes(ET.tostring(root))
    core = pytest.importorskip("meshioplusplus._core")
    for read in (
        core.vtp_read,
        importlib.import_module("meshioplusplus.vtp._vtp").read,
    ):
        with pytest.raises(mio.ReadError):
            read(str(path))

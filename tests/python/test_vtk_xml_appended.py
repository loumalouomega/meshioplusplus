"""Independent framing oracles: no meshio++ writer produces these payloads."""

import base64
import importlib
import zlib

import numpy as np
import pytest

import meshioplusplus as mio

FORMATS = {
    "vtp": "PolyData",
    "vts": "StructuredGrid",
    "vtr": "RectilinearGrid",
    "vti": "ImageData",
}


def appended_fixture(
    fmt, encoding="raw", header="UInt32", order="LittleEndian", codec=None, fields=False
):
    """One cube or surface quad, with multi-component and integer fields."""
    bo = "<" if order == "LittleEndian" else ">"
    hd = np.dtype(bo + ("u4" if header == "UInt32" else "u8"))
    payload = bytearray()

    def array(name, values, vtk_type, nc=0):
        dtype = {"Float64": "f8", "Int64": "i8", "Int32": "i4"}[vtk_type]
        raw = np.asarray(values, dtype=bo + dtype).tobytes()
        if codec:
            # Multiple blocks, independently compressed, including a short last block.
            chunks = [raw[k : k + 32] for k in range(0, len(raw), 32)]
            if codec == "lz4":
                lz4 = pytest.importorskip("lz4.block")
                blocks = [lz4.compress(c, store_size=False) for c in chunks]
            elif codec == "zstd":
                zstd = pytest.importorskip("zstandard")
                blocks = [zstd.ZstdCompressor().compress(c) for c in chunks]
            else:
                blocks = [zlib.compress(c) for c in chunks]
            h = np.array(
                [len(blocks), 32, len(chunks[-1]) if chunks else 0]
                + [len(b) for b in blocks],
                dtype=hd,
            ).tobytes()
            body = b"".join(blocks)
        else:
            h, body = np.array([len(raw)], dtype=hd).tobytes(), raw
        offset = len(payload)
        # Padded header and body encoded separately, as VTK does.
        payload.extend(
            h + body
            if encoding == "raw"
            else base64.b64encode(h) + base64.b64encode(body)
        )
        comps = f' NumberOfComponents="{nc}"' if nc else ""
        return f'<DataArray Name="{name}" type="{vtk_type}"{comps} format="appended" offset="{offset}"/>'

    points = np.array(
        [[x, y, z] for z in range(2) for y in range(2) for x in range(2)], dtype=float
    )
    if fmt == "vtp":
        points = points[:4]
    n = len(points)
    pd = array("vec", np.arange(n * 3).reshape(n, 3) * 0.5, "Float64", 3)
    pd += array("tag", np.arange(n), "Int64")
    cd = array("mat", [19], "Int32")
    fd = ""
    if fields:
        fd = (
            "<FieldData>"
            + array("units", [42], "Int64")
            + array("region:point:left", [0, 2], "Int64")
            + "</FieldData>"
        )
    geometry = ""
    if fmt in ("vtp", "vts"):
        geometry = "<Points>" + array("Points", points, "Float64", 3) + "</Points>"
    if fmt == "vtp":
        geometry += (
            "<Polys>"
            + array("connectivity", [0, 1, 3, 2], "Int64")
            + array("offsets", [4], "Int64")
            + "</Polys>"
        )
        # An explicitly empty section checks zero-block compressed arrays too.
        geometry += (
            "<Verts>"
            + array("connectivity", [], "Int64")
            + array("offsets", [], "Int64")
            + "</Verts>"
        )
        attrs = 'NumberOfPoints="4" NumberOfVerts="0" NumberOfLines="0" NumberOfStrips="0" NumberOfPolys="1"'
        grid_attrs = ""
    else:
        attrs = 'Extent="0 1 0 1 0 1"'
        grid_attrs = ' WholeExtent="0 1 0 1 0 1"'
        if fmt == "vti":
            grid_attrs += ' Origin="0 0 0" Spacing="1 1 1"'
        if fmt == "vtr":
            geometry = (
                "<Coordinates>"
                + "".join(
                    array(f"{axis}_coordinates", [0, 1], "Float64") for axis in "xyz"
                )
                + "</Coordinates>"
            )
    compressors = {
        "zlib": "vtkZLibDataCompressor",
        "lz4": "vtkLZ4DataCompressor",
        "zstd": "vtkZSTDDataCompressor",
    }
    compressor = f' compressor="{compressors[codec]}"' if codec else ""
    dataset = FORMATS[fmt]
    xml = f'<VTKFile type="{dataset}" version="1.0" byte_order="{order}" header_type="{header}"{compressor}><{dataset}{grid_attrs}>{fd}<Piece {attrs}><PointData>{pd}</PointData><CellData>{cd}</CellData>{geometry}</Piece></{dataset}><AppendedData encoding="{encoding}">_'
    return xml.encode() + payload + b"</AppendedData></VTKFile>", points


def readers(fmt):
    core = pytest.importorskip("meshioplusplus._core")
    reference = importlib.import_module(f"meshioplusplus.{fmt}._{fmt}").read
    return [lambda path: getattr(core, fmt + "_read")(str(path)), reference]


def assert_mesh(mesh, points):
    np.testing.assert_array_equal(mesh.points, points)
    np.testing.assert_array_equal(
        mesh.point_data["tag"].reshape(-1), np.arange(len(points))
    )
    np.testing.assert_array_equal(
        mesh.point_data["vec"], np.arange(len(points) * 3).reshape(-1, 3) * 0.5
    )
    np.testing.assert_array_equal(mesh.cell_data["mat"][0].reshape(-1), [19])
    np.testing.assert_array_equal(
        mesh.cells[0].data,
        [[0, 1, 3, 2]] if len(points) == 4 else [[0, 1, 3, 2, 4, 5, 7, 6]],
    )


@pytest.mark.parametrize("fmt", FORMATS)
@pytest.mark.parametrize("encoding", ["raw", "base64"])
@pytest.mark.parametrize("header", ["UInt32", "UInt64"])
@pytest.mark.parametrize("order", ["LittleEndian", "BigEndian"])
@pytest.mark.parametrize("codec", [None, "zlib"])
def test_independent_payloads(fmt, encoding, header, order, codec, tmp_path):
    core = pytest.importorskip("meshioplusplus._core")
    if codec and not core.__has_zlib__:
        pytest.skip("native build lacks zlib")
    blob, points = appended_fixture(fmt, encoding, header, order, codec)
    path = tmp_path / f"fixture.{fmt}"
    path.write_bytes(blob)
    for reader in readers(fmt):
        assert_mesh(reader(path), points)
    meta = core.read_metadata(str(path), fmt)
    assert meta["num_points"] == len(points)
    assert meta["num_cells"] == 1
    assert meta["point_data_names"] == ["tag", "vec"]
    assert meta["cell_data_names"] == ["mat"]
    # Ownership: deleting the input cannot invalidate decoded buffers.
    mesh = readers(fmt)[0](path)
    path.unlink()
    assert_mesh(mesh, points)


@pytest.mark.parametrize("fmt", FORMATS)
@pytest.mark.parametrize("encoding", ["raw", "base64"])
def test_selective_and_strict_dispatch(fmt, encoding, tmp_path, monkeypatch):
    from meshioplusplus import _fallback

    monkeypatch.setenv("MESHIOPLUSPLUS_STRICT_CORE", "1")
    monkeypatch.setattr(_fallback, "_strict", True)
    blob, points = appended_fixture(fmt, encoding)
    path = tmp_path / f"fixture.{fmt}"
    path.write_bytes(blob)
    mesh = mio.read(path, arrays=["tag"])
    assert list(mesh.point_data) == ["tag"]
    assert not mesh.cell_data
    mesh = mio.read(path, points_only=True)
    np.testing.assert_array_equal(mesh.points, points)
    assert len(mesh.cells) == 1
    assert not mesh.point_data and not mesh.cell_data


@pytest.mark.parametrize("fmt", FORMATS)
@pytest.mark.parametrize("codec", ["lz4", "zstd"])
@pytest.mark.parametrize("encoding", ["raw", "base64"])
def test_optional_codecs_or_explicit_native_decline(fmt, codec, encoding, tmp_path):
    core = pytest.importorskip("meshioplusplus._core")
    blob, points = appended_fixture(fmt, encoding, "UInt64", "BigEndian", codec)
    path = tmp_path / f"codec.{fmt}"
    path.write_bytes(blob)
    native, reference = readers(fmt)
    assert_mesh(reference(path), points)
    if getattr(core, f"__has_{codec}__", False):
        assert_mesh(native(path), points)
    else:
        with pytest.raises(mio.ReadError, match=codec.upper()):
            native(path)
        with pytest.raises(mio.ReadError, match=codec.upper()):
            core.read_metadata(str(path), fmt)


@pytest.mark.parametrize("fmt", FORMATS)
@pytest.mark.parametrize(
    "codec, compressor",
    [("lz4", "vtkLZ4DataCompressor"), ("zstd", "vtkZSTDDataCompressor")],
)
def test_native_missing_codec_fails_before_payload_decode(
    fmt, codec, compressor, tmp_path
):
    core = pytest.importorskip("meshioplusplus._core")
    if getattr(core, f"__has_{codec}__", False):
        pytest.skip("native codec is enabled")
    blob, _ = appended_fixture(fmt, codec="zlib")
    blob = blob.replace(b"vtkZLibDataCompressor", compressor.encode())
    path = tmp_path / f"missing_codec.{fmt}"
    path.write_bytes(blob)
    with pytest.raises(mio.ReadError, match=codec.upper()):
        readers(fmt)[0](path)


@pytest.mark.parametrize("encoding", ["raw", "base64"])
def test_vtp_appended_field_data_and_regions(encoding, tmp_path):
    blob, points = appended_fixture("vtp", encoding, fields=True)
    path = tmp_path / "regions.vtp"
    path.write_bytes(blob)
    for reader in readers("vtp"):
        mesh = reader(path)
        assert_mesh(mesh, points)
        np.testing.assert_array_equal(mesh.field_data["units"], [42])
        assert [(r.name, r.kind, list(r.entries)) for r in mesh.regions] == [
            ("left", "point", [0, 2])
        ]
    mesh = mio.read(path, points_only=True)
    assert not mesh.point_data and not mesh.field_data
    assert [(r.name, list(r.entries)) for r in mesh.regions] == [("left", [0, 2])]


@pytest.mark.parametrize("fmt", FORMATS)
@pytest.mark.parametrize(
    "bad", [b"-1", b"nope", b"3junk", b"18446744073709551616", b"999999999"]
)
def test_invalid_offsets(fmt, bad, tmp_path):
    blob, _ = appended_fixture(fmt)
    blob = blob.replace(b'offset="0"', b'offset="' + bad + b'"', 1)
    path = tmp_path / f"bad.{fmt}"
    path.write_bytes(blob)
    for reader in readers(fmt):
        with pytest.raises(mio.ReadError):
            reader(path)


@pytest.mark.parametrize("fmt", FORMATS)
@pytest.mark.parametrize(
    "bad", ["marker", "truncated", "size", "header", "encoding", "missing"]
)
def test_malformed_payloads(fmt, bad, tmp_path):
    blob, _ = appended_fixture(fmt)
    start = blob.index(b">_", blob.index(b"<AppendedData")) + 2
    if bad == "marker":
        blob = blob[: start - 1] + b"!" + blob[start:]
    elif bad == "truncated":
        blob = blob[:start] + b"\x00\x00" + b"</AppendedData></VTKFile>"
    elif bad == "size":
        blob = blob[:start] + b"\xff" * 4 + blob[start + 4 :]
    elif bad == "header":
        blob = blob.replace(b'header_type="UInt32"', b'header_type="Int32"')
    elif bad == "encoding":
        blob, _ = appended_fixture(fmt, "base64")
        blob = blob.replace(b'encoding="base64"', b'encoding="unknown"')
    else:
        blob = blob[: blob.index(b"<AppendedData")] + b"</VTKFile>"
    path = tmp_path / f"bad.{fmt}"
    path.write_bytes(blob)
    for reader in readers(fmt):
        with pytest.raises(mio.ReadError):
            reader(path)


@pytest.mark.parametrize("fmt", FORMATS)
@pytest.mark.parametrize("nc", ["0", "7"])
def test_inconsistent_components_fail_instead_of_ignoring_reshape(fmt, nc, tmp_path):
    from meshioplusplus._exceptions import CorruptionError

    blob, _ = appended_fixture(fmt)
    blob = blob.replace(
        b'NumberOfComponents="3"', f'NumberOfComponents="{nc}"'.encode(), 1
    )
    path = tmp_path / f"bad_components.{fmt}"
    path.write_bytes(blob)
    for reader in readers(fmt):
        with pytest.raises((mio.ReadError, CorruptionError)):
            reader(path)


@pytest.mark.parametrize("fmt", FORMATS)
@pytest.mark.parametrize("encoding", ["raw", "base64"])
def test_external_vtk_writer(fmt, encoding, tmp_path):
    vtk = pytest.importorskip("vtk")
    from vtk.util.numpy_support import numpy_to_vtk

    points = np.array(
        [[x, y, z] for z in range(2) for y in range(2) for x in range(2)], dtype=float
    )
    if fmt == "vtp":
        points = points[:4]
        data = vtk.vtkPolyData()
        cells = vtk.vtkCellArray()
        cells.InsertNextCell(4, [0, 1, 3, 2])
        data.SetPolys(cells)
    else:
        data = {
            "vts": vtk.vtkStructuredGrid,
            "vtr": vtk.vtkRectilinearGrid,
            "vti": vtk.vtkImageData,
        }[fmt]()
        data.SetDimensions(2, 2, 2)
        if fmt == "vtr":
            for setter in (
                data.SetXCoordinates,
                data.SetYCoordinates,
                data.SetZCoordinates,
            ):
                setter(numpy_to_vtk(np.array([0.0, 1.0])))
    if fmt in ("vtp", "vts"):
        vp = vtk.vtkPoints()
        vp.SetData(numpy_to_vtk(points))
        data.SetPoints(vp)
    for name, values in (
        ("tag", np.arange(len(points), dtype=np.int64)),
        ("vec", np.arange(len(points) * 3).reshape(-1, 3) * 0.5),
    ):
        arr = numpy_to_vtk(values)
        arr.SetName(name)
        data.GetPointData().AddArray(arr)
    arr = numpy_to_vtk(np.array([19], dtype=np.int32))
    arr.SetName("mat")
    data.GetCellData().AddArray(arr)
    path = tmp_path / f"vtk.{fmt}"
    writer = {
        "vtp": vtk.vtkXMLPolyDataWriter,
        "vts": vtk.vtkXMLStructuredGridWriter,
        "vtr": vtk.vtkXMLRectilinearGridWriter,
        "vti": vtk.vtkXMLImageDataWriter,
    }[fmt]()
    writer.SetInputData(data)
    writer.SetFileName(str(path))
    writer.SetDataModeToAppended()
    writer.SetEncodeAppendedData(encoding == "base64")
    writer.SetHeaderTypeToUInt64()
    writer.SetCompressorTypeToZLib()
    assert writer.Write() == 1
    for reader in readers(fmt):
        assert_mesh(reader(path), points)

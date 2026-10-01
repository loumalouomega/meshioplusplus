"""Independent multi-piece XML and legacy structured VTK fixtures."""

import importlib
import xml.etree.ElementTree as ET

import numpy as np
import pytest

import meshioplusplus as mio

TYPES = {"vtp": "PolyData", "vts": "StructuredGrid", "vtr": "RectilinearGrid", "vti": "ImageData"}


def xml_fixture(fmt, missing=False, incompatible=False, local_regions=True):
    root = ET.Element("VTKFile", type=TYPES[fmt], byte_order="LittleEndian")
    attrs = {} if fmt == "vtp" else {"WholeExtent": "2 4 0 1 0 1"}
    grid = ET.SubElement(root, TYPES[fmt], attrs)

    def array(parent, name, data, dtype="Int64", components=None):
        attrs = {"Name": name, "type": dtype, "format": "ascii"}
        if components:
            attrs["NumberOfComponents"] = str(components)
        da = ET.SubElement(parent, "DataArray", attrs)
        da.text = " ".join(str(x) for x in np.asarray(data).reshape(-1))
        return da

    array(ET.SubElement(grid, "FieldData"), "units", [42])
    expected = []
    for part in range(2):
        points = np.array([[x, y, z] for z in range(2 if fmt != "vtp" else 1) for y in range(2) for x in (2 + part, 3 + part)], dtype=float)
        expected.append(points)
        attrs = {"NumberOfPoints": "4", "NumberOfPolys": "1"} if fmt == "vtp" else {"Extent": f"{2 + part} {3 + part} 0 1 0 1"}
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
            array(pd, "optional", np.ones(len(points)), "Float32" if incompatible and part else "Float64")
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
@pytest.mark.parametrize("missing,incompatible", [(False, False), (True, False), (False, True)])
def test_xml_piece_assembly(tmp_path, monkeypatch, fmt, missing, incompatible):
    data, points = xml_fixture(fmt, missing, incompatible)
    path = tmp_path / f"pieces.{fmt}"
    path.write_bytes(data)
    core = pytest.importorskip("meshioplusplus._core")
    readers = [getattr(core, fmt + "_read"), importlib.import_module(f"meshioplusplus.{fmt}._{fmt}").read]
    for read in readers:
        mesh = read(str(path))
        np.testing.assert_array_equal(mesh.points, points)
        assert sum(len(cb.data) for cb in mesh.cells) == 2
        np.testing.assert_array_equal(np.concatenate(mesh.cell_data["material"]), [20, 21])
        np.testing.assert_array_equal(np.concatenate(mesh.cell_data["vtkGhostType"]), [0, 1])
        assert ("optional" in mesh.point_data) == (not missing and not incompatible)
        regions = {(r.kind, r.name): r.entries for r in mesh.regions}
        np.testing.assert_array_equal(regions["point", "corner"], [0, len(points) // 2])
        np.testing.assert_array_equal(regions["cell", "body"], [0, 1])
        np.testing.assert_array_equal(regions["side", "wall"], [[0, 1], [1, 1]])
        assert regions["point", "empty"].size == 0
        np.testing.assert_array_equal(mesh.field_data["units"], [42])
        assert mesh.cells[0].data[-1].min() >= len(points) // 2
    monkeypatch.setenv("MESHIOPLUSPLUS_STRICT_CORE", "1")
    assert len(mio.read(path).points) == len(points)
    metadata = mio.read_metadata(path)
    assert metadata["num_points"] == len(points)
    assert ("optional" in metadata["point_data_names"]) == (not missing and not incompatible)


def legacy_fixture(dataset, dims, binary=False, version="4.2"):
    parts = [f"# vtk DataFile Version {version}\nindependent structured fixture\n{'BINARY' if binary else 'ASCII'}\nDATASET {dataset}\nDIMENSIONS {' '.join(map(str, dims))}\n".encode()]

    def values(data, dtype="f8"):
        arr = np.asarray(data, dtype=(">" if binary else "=") + dtype)
        parts.append(arr.tobytes() if binary else " ".join(str(x) for x in arr.reshape(-1)).encode())
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
    parts.append(f"POINT_DATA {len(points)}\nSCALARS tag int 1\nLOOKUP_TABLE default\n".encode())
    values(np.arange(len(points)), "i4")
    parts.append(b"VECTORS velocity double\n")
    values(np.arange(len(points) * 3).reshape(-1, 3))
    parts.append(f"CELL_DATA {nc}\nTENSORS stress double\n".encode())
    values(np.arange(nc * 9).reshape(-1, 3, 3))
    return b"".join(parts), points


@pytest.mark.parametrize("dataset", ["STRUCTURED_POINTS", "STRUCTURED_GRID", "RECTILINEAR_GRID"])
@pytest.mark.parametrize("dims", [(3, 2, 2), (3, 2, 1), (1, 3, 2), (1, 1, 4)])
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
        np.testing.assert_array_equal(native.point_data[name], reference.point_data[name])
    np.testing.assert_array_equal(native.cell_data["stress"][0], reference.cell_data["stress"][0])
    monkeypatch.setenv("MESHIOPLUSPLUS_STRICT_CORE", "1")
    assert len(mio.read(path).points) == len(points)

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


@pytest.mark.parametrize("writer", [meshioplusplus.ensight.write, _ensight.write])
def test_ragged_write_raises(writer, tmp_path):
    # The C++ path raises too; the shim then falls back to the Python writer,
    # which raises the WriteError seen here.
    with pytest.raises(meshioplusplus.WriteError):
        writer(tmp_path / "poly.case", copy.deepcopy(helpers.polygon_mesh))


def test_wildcard_model_raises(tmp_path):
    p = tmp_path / "transient.case"
    p.write_text(
        "FORMAT\ntype: ensight gold\n\nGEOMETRY\nmodel: 1 transient.****.geo\n"
    )
    with pytest.raises(meshioplusplus.ReadError):
        _ensight.read(p)

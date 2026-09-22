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

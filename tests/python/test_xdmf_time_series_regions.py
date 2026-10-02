"""Shared named regions in temporal XDMF, independent of the field steps."""

import xml.etree.ElementTree as ET

import numpy as np
import pytest

import meshioplusplus as mio
from meshioplusplus import _core, _fallback


def region_mesh(mixed=True):
    points = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 1.0, 0.0],
        ]
    )
    cells = [("tetra", np.array([[0, 1, 2, 3]], dtype=np.int64))]
    if mixed:
        cells += [
            ("triangle", np.array([[1, 4, 2]], dtype=np.int64)),
            ("line", np.array([[0, 4]], dtype=np.int64)),
        ]
    regions = [
        mio.Region('anchors & "rim"', "point", [0, 4], dim=0, tag=7),
        mio.Region("empty points", "point", [], dim=0, tag=8),
        mio.Region("volume", "cell", [0], dim=3, tag=21),
        mio.Region("empty cells", "cell", [], dim=3, tag=22),
        mio.Region(
            "wall", "side", [[0, 1], [1, 2]] if mixed else [[0, 1]], dim=2, tag=9
        ),
        mio.Region("empty sides", "side", np.empty((0, 2), dtype=np.int64)),
    ]
    return mio.Mesh(points, cells, regions=regions)


def assert_regions(actual, expected):
    actual = {(r.kind, r.name, r.dim, r.tag): r.entries for r in actual}
    expected = {(r.kind, r.name, r.dim, r.tag): r.entries for r in expected}
    assert actual.keys() == expected.keys()
    for key in expected:
        np.testing.assert_array_equal(actual[key], expected[key])


def write_series(path, fmt, native, mesh):
    if fmt == "HDF":
        pytest.importorskip("h5py")
    if native:
        if not hasattr(_core, "XdmfTimeSeriesWriter"):
            pytest.skip("native series writer unavailable")
        if fmt == "HDF" and not getattr(_core, "__has_hdf5__", False):
            pytest.skip("native HDF5 unavailable")
        with _core.XdmfTimeSeriesWriter(str(path), fmt) as writer:
            writer.write_points_cells(mesh)
            for k in range(3):
                step = mio.Mesh(
                    mesh.points,
                    mesh.cells,
                    point_data={"phi": np.full(len(mesh.points), float(k))},
                    cell_data={
                        "rho": [
                            np.array([float(k + b)]) for b in range(len(mesh.cells))
                        ]
                    },
                )
                writer.write_data(k * 0.5, step)
    else:
        with mio.xdmf.TimeSeriesWriter(path, fmt) as writer:
            writer.write_points_cells(mesh.points, mesh.cells, regions=mesh.regions)
            for k in range(3):
                writer.write_data(
                    k * 0.5,
                    point_data={"phi": np.full(len(mesh.points), float(k))},
                    cell_data={
                        "rho": [
                            np.array([float(k + b)]) for b in range(len(mesh.cells))
                        ]
                    },
                )


@pytest.mark.parametrize("fmt", ["XML", "Binary", "HDF"])
@pytest.mark.parametrize("native", [False, True])
@pytest.mark.parametrize("mixed", [False, True])
def test_series_regions_roundtrip(tmp_path, monkeypatch, fmt, native, mixed):
    # Preserve the Python writer's documented CWD-relative heavy-data contract.
    monkeypatch.chdir(tmp_path)
    mesh = region_mesh(mixed)
    original = [cb.data.copy() for cb in mesh.cells]
    path = tmp_path / "series.xdmf"
    write_series(path, fmt, native, mesh)
    for block, expected in zip(mesh.cells, original):
        np.testing.assert_array_equal(block.data, expected)
    root = ET.parse(path).getroot()
    static = root.find("Domain/Grid[@Name='mesh']")
    sets = static.findall("Set")
    assert len(sets) == (7 if mixed else 6)
    assert [s.get("SetType") for s in sets if s.get("Name") == "wall"] == (
        ["Face", "Edge"] if mixed else ["Face"]
    )
    collection = root.find("Domain/Grid[@CollectionType='Temporal']")
    for step in collection.findall("Grid"):
        assert not step.findall("Set")
        include = step.find("{http://www.w3.org/2003/XInclude}include")
        assert "self::Set" in include.get("xpointer")
    # Reading must also work with a CWD unrelated to the series directory.
    other = tmp_path / "other"
    other.mkdir()
    monkeypatch.chdir(other)
    with mio.xdmf.TimeSeriesReader(path) as reader:
        assert reader.regions == []
        points, cells = reader.read_points_cells()
        np.testing.assert_array_equal(points, mesh.points)
        for a, b in zip(cells, mesh.cells):
            np.testing.assert_array_equal(a.data, b.data)
        assert_regions(reader.regions, mesh.regions)
        for k in range(3):
            t, pd, cd = reader.read_data(k)
            assert t == k * 0.5
            np.testing.assert_array_equal(
                pd["phi"], np.full(len(mesh.points), float(k))
            )
            np.testing.assert_array_equal(
                np.concatenate(cd["rho"]), np.arange(len(mesh.cells)) + k
            )
            assert_regions(reader.regions, mesh.regions)
        reader.read_points_cells()
        assert_regions(reader.regions, mesh.regions)
    if hasattr(_core, "xdmf_read") and (
        fmt != "HDF" or getattr(_core, "__has_hdf5__", False)
    ):
        monkeypatch.setattr(_fallback, "_strict", True)
        for k in (0, 1, -1):
            got = mio.read(path, time_step=k)
            assert_regions(got.regions, mesh.regions)
            assert_regions(
                mio.read(path, time_step=k, points_only=True).regions, mesh.regions
            )
            assert_regions(
                mio.read(path, time_step=k, arrays=["phi"]).regions, mesh.regions
            )
        metadata = mio.read_metadata(path)
        assert len(metadata["regions"]) == 6
        if not mixed:
            np.testing.assert_array_equal(metadata["time_values"], [0.0, 0.5, 1.0])


def test_python_legacy_arguments_and_dict_cells(tmp_path):
    mesh = region_mesh(False)
    path = tmp_path / "legacy.xdmf"
    with mio.xdmf.TimeSeriesWriter(path, "XML") as writer:
        writer.write_points_cells(mesh.points, {cb.type: cb.data for cb in mesh.cells})
        writer.write_data(0.0)
    with mio.xdmf.TimeSeriesReader(path) as reader:
        assert len(reader.read_points_cells()) == 2
        assert reader.regions == []
        assert len(reader.read_data(0)) == 3


def test_python_mixed_vertex_topology_preserves_regions_and_input(tmp_path):
    mesh = region_mesh(False)
    mesh.cells.append(mio.CellBlock("vertex", np.array([[4]], dtype=np.int64)))
    mesh.regions.append(mio.Region("vertex", "cell", [1], dim=0, tag=31))
    vertex = mesh.cells[-1].data.copy()
    path = tmp_path / "vertex.xdmf"
    write_series(path, "XML", False, mesh)
    np.testing.assert_array_equal(mesh.cells[-1].data, vertex)
    with mio.xdmf.TimeSeriesReader(path) as reader:
        _, cells = reader.read_points_cells()
        np.testing.assert_array_equal(cells[-1].data, vertex)
        assert_regions(reader.regions, mesh.regions)


@pytest.mark.parametrize("fmt", ["XML", "Binary", "HDF"])
def test_shared_regions_survive_flush_append_and_array_steps(tmp_path, fmt):
    if not hasattr(_core, "XdmfTimeSeriesWriter"):
        pytest.skip("native writer unavailable")
    if fmt == "HDF" and not getattr(_core, "__has_hdf5__", False):
        pytest.skip("native HDF5 unavailable")
    mesh = region_mesh(False)
    path = tmp_path / "append.xdmf"
    with _core.XdmfTimeSeriesWriter(str(path), fmt) as writer:
        writer.write_points_cells(mesh)
        writer.write_data_arrays(
            0.0, {"phi": np.arange(len(mesh.points), dtype=float)}, {}
        )
        writer.flush()
        assert_regions(_core.xdmf_read(str(path)).regions, mesh.regions)
    with _core.XdmfTimeSeriesWriter(str(path), fmt, mode="append") as writer:
        writer.write_data_arrays(1.0, {"phi": np.ones(len(mesh.points))}, {})
    got = _core.xdmf_read(str(path), time_step=-1)
    assert_regions(got.regions, mesh.regions)
    assert (
        len(ET.parse(path).getroot().find("Domain/Grid[@Name='mesh']").findall("Set"))
        == 6
    )


@pytest.mark.parametrize("malformed", ["missing", "length"])
def test_python_reader_rejects_malformed_shared_side_set(tmp_path, malformed):
    mesh = region_mesh(False)
    path = tmp_path / "bad.xdmf"
    write_series(path, "XML", False, mesh)
    root = ET.parse(path).getroot()
    side = root.find("Domain/Grid[@Name='mesh']/Set[@Name='wall']")
    local = side.findall("DataItem")[1]
    if malformed == "missing":
        side.remove(local)
    else:
        local.set("Dimensions", "2")
        local.text = "1 2"
    ET.ElementTree(root).write(path)
    with mio.xdmf.TimeSeriesReader(path) as reader:
        with pytest.raises(mio.ReadError):
            reader.read_points_cells()

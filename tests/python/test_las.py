"""LAS (roadmap 2.4.6.1): files written by laspy, read in both engines; ours read by
laspy."""

import numpy as np
import pytest

import meshioplusplus as mio
from meshioplusplus import _fallback
from meshioplusplus.las import _las

laspy = pytest.importorskip("laspy")


def _laspy_file(path, fmt, version):
    rng = np.random.default_rng(fmt)
    n = 50
    las = laspy.create(point_format=fmt, file_version=version)
    las.header.offsets = [100.0, 200.0, 5.0]
    las.header.scales = [0.01, 0.01, 0.001]
    las.x = rng.uniform(100, 110, n)
    las.y = rng.uniform(200, 210, n)
    las.z = rng.uniform(5, 6, n)
    las.intensity = rng.integers(0, 60000, n)
    las.return_number = rng.integers(1, 4, n)
    las.number_of_returns = np.full(n, 3)
    las.classification = rng.integers(0, 20, n)
    las.point_source_id = rng.integers(0, 1000, n)
    if "gps_time" in las.point_format.dimension_names:
        las.gps_time = rng.uniform(0, 1e6, n)
    if "red" in las.point_format.dimension_names:
        las.red = rng.integers(0, 65535, n)
        las.green = rng.integers(0, 65535, n)
        las.blue = rng.integers(0, 65535, n)
    las.write(path)
    return laspy.read(path)


@pytest.fixture(params=["core", "python"])
def reader(request, monkeypatch):
    if request.param == "core":
        monkeypatch.setattr(_fallback, "_strict", True)
        return mio.las.read
    return _las.read


@pytest.mark.parametrize(
    "fmt, version",
    [
        (0, "1.2"),
        (1, "1.2"),
        (3, "1.2"),
        (2, "1.3"),
        (6, "1.4"),
        (7, "1.4"),
        (8, "1.4"),
    ],
)
def test_laspy_files_read(tmp_path, reader, fmt, version):
    path = tmp_path / f"f{fmt}.las"
    ref = _laspy_file(path, fmt, version)
    mesh = reader(path)
    assert np.allclose(mesh.points, np.column_stack([ref.x, ref.y, ref.z]))
    assert mesh.cells[0].type == "vertex" and len(mesh.cells[0].data) == 50
    for name in (
        "intensity",
        "return_number",
        "number_of_returns",
        "classification",
        "point_source_id",
    ):
        assert np.array_equal(mesh.point_data[name], np.asarray(ref[name])), name
    if "gps_time" in ref.point_format.dimension_names:
        assert np.allclose(mesh.point_data["gps_time"], ref.gps_time)
    if "red" in ref.point_format.dimension_names:
        assert np.array_equal(mesh.point_data["blue"], np.asarray(ref.blue))
    assert mesh.field_data["las:point_format"].tolist() == [fmt]
    every = reader(path, every=10)
    assert len(every.points) == 5
    picked = reader(path, arrays=["intensity"])
    assert set(picked.point_data) == {"intensity"}


@pytest.mark.parametrize("fmt", [6, 7, 8])
def test_we_write_what_laspy_reads_with_equal_bytes(tmp_path, fmt):
    src = tmp_path / "src.las"
    ref = _laspy_file(src, fmt, "1.4")
    mesh = mio.las.read(src)
    a, b = tmp_path / "a.las", tmp_path / "b.las"
    mio.las.write(a, mesh)
    _las.write(b, mesh)
    assert a.read_bytes() == b.read_bytes()
    back = laspy.read(a)
    assert back.header.point_format.id == fmt
    assert np.allclose(back.x, ref.x) and np.allclose(back.z, ref.z)
    assert np.array_equal(
        np.asarray(back.classification), np.asarray(ref.classification)
    )
    assert np.array_equal(np.asarray(back.return_number), np.asarray(ref.return_number))


def test_any_point_set_writes_and_errors(tmp_path, reader):
    mesh = mio.Mesh(
        np.array([[0.5, 1.5, 2.5], [3.25, 4.0, 5.0]]), [("line", np.array([[0, 1]]))]
    )
    path = tmp_path / "p.las"
    mio.write(path, mesh)
    back = reader(path)
    assert np.allclose(back.points, mesh.points, atol=1e-3)
    bad = tmp_path / "bad.las"
    bad.write_bytes(b"not a las file" * 30)
    with pytest.raises(mio.ReadError, match="LASF"):
        reader(bad)

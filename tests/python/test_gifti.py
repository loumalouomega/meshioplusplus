"""GIFTI (roadmap 2.4.4): both engines, and nibabel both ways."""

import numpy as np
import pytest

import meshioplusplus as mio
from meshioplusplus import _fallback
from meshioplusplus._regions import Region
from meshioplusplus.gifti import _gifti

nib = pytest.importorskip("nibabel")


def _surface():
    pts = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
    tri = np.array([[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]])
    mesh = mio.Mesh(
        pts,
        [("triangle", tri)],
        point_data={
            "thickness": np.array([1.5, 2.5, 3.5, 4.5]),
            "gifti:label": np.array([1, 1, 2, 2], dtype=np.int32),
        },
        regions=[
            Region("cortex", "point", [0, 1], tag=1),
            Region("wall", "point", [2, 3], tag=2),
        ],
    )
    return mesh


@pytest.fixture(params=["core", "python"])
def engine(request, monkeypatch):
    if request.param == "core":
        monkeypatch.setattr(_fallback, "_strict", True)
        return mio.gifti.read, mio.gifti.write
    return _gifti.read, _gifti.write


def test_round_trip_and_equal_bytes(tmp_path, engine):
    read, write = engine
    path = tmp_path / "s.gii"
    write(path, _surface())
    back = read(path)
    assert np.allclose(back.points, _surface().points)
    assert np.asarray(back.cells[0].data).tolist() == _surface().cells[0].data.tolist()
    assert np.allclose(back.point_data["thickness"], [1.5, 2.5, 3.5, 4.5])
    regions = {r.name: (np.asarray(r.entries).tolist(), r.tag) for r in back.regions}
    assert regions == {"cortex": ([0, 1], 1), "wall": ([2, 3], 2)}
    other = tmp_path / "o.gii"
    (_gifti.write if write is mio.gifti.write else mio.gifti.write)(other, _surface())
    assert other.read_bytes() == path.read_bytes()


def test_nibabel_reads_what_we_write(tmp_path):
    path = tmp_path / "s.gii"
    mio.write(path, _surface())
    img = nib.load(path)
    coords, faces = img.agg_data(("pointset", "triangle"))
    assert np.allclose(coords, _surface().points)
    assert np.asarray(faces).tolist() == _surface().cells[0].data.tolist()
    names = [da.meta.get("Name") for da in img.darrays]
    assert "thickness" in names
    assert img.labeltable.get_labels_as_dict() == {1: "cortex", 2: "wall"}


@pytest.mark.parametrize("encoding", ["ASCII", "Base64Binary", "GZipBase64Binary"])
@pytest.mark.parametrize("endian", ["little", "big"])
def test_we_read_what_nibabel_writes(tmp_path, engine, encoding, endian):
    read, _ = engine
    gi = nib.gifti
    order = "<" if endian == "little" else ">"
    pts = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]], dtype=np.float32)
    img = gi.GiftiImage()
    kw = dict(encoding=encoding, endian=endian)
    img.add_gifti_data_array(
        gi.GiftiDataArray(
            pts.astype(order + "f4"), intent="NIFTI_INTENT_POINTSET", **kw
        )
    )
    img.add_gifti_data_array(
        gi.GiftiDataArray(
            np.array([[0, 1, 2]], dtype=order + "i4"),
            intent="NIFTI_INTENT_TRIANGLE",
            **kw,
        )
    )
    shape = gi.GiftiDataArray(
        np.array([7.0, 8.0, 9.0], dtype=order + "f4"), intent="NIFTI_INTENT_SHAPE", **kw
    )
    shape.meta["Name"] = "curv"
    img.add_gifti_data_array(shape)
    path = tmp_path / f"n_{encoding}_{endian}.gii"
    nib.save(img, path)
    mesh = read(path)
    assert np.allclose(mesh.points, pts)
    assert np.asarray(mesh.cells[0].data).tolist() == [[0, 1, 2]]
    assert np.allclose(mesh.point_data["curv"], [7, 8, 9])


def test_errors(tmp_path, engine):
    read, write = engine
    bad = tmp_path / "b.gii"
    bad.write_text("<GIFTI></GIFTI>")
    with pytest.raises(mio.ReadError, match="POINTSET"):
        read(bad)
    quads = mio.Mesh(np.zeros((4, 3)), [("quad", np.array([[0, 1, 2, 3]]))])
    with pytest.raises(mio.WriteError, match="triangle"):
        write(tmp_path / "q.gii", quads)

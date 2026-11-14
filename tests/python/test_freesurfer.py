"""FreeSurfer surfaces, curv and annot (roadmap 2.4.4): both engines, nibabel both ways."""

import numpy as np
import pytest

import meshioplusplus as mio
from meshioplusplus import _fallback
from meshioplusplus.freesurfer import _freesurfer

nib = pytest.importorskip("nibabel")
fsio = nib.freesurfer.io

PTS = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=np.float32)
TRI = np.array([[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]], dtype=np.int32)


@pytest.fixture(params=["core", "python"])
def engine(request, monkeypatch):
    if request.param == "core":
        monkeypatch.setattr(_fallback, "_strict", True)
        return mio.freesurfer.read, mio.freesurfer.write
    return _freesurfer.read, _freesurfer.write


@pytest.fixture
def subject(tmp_path):
    fsio.write_geometry(tmp_path / "lh.pial", PTS, TRI, create_stamp="created by test")
    fsio.write_morph_data(tmp_path / "lh.thickness", np.array([1, 2, 3, 4], np.float32))
    ctab = np.array([[25, 5, 25, 0, 0], [220, 20, 10, 0, 0]], dtype=np.int64)
    fsio.write_annot(
        tmp_path / "lh.aparc.annot",
        np.array([0, 0, 1, 1]),
        ctab,
        [b"bankssts", b"frontal"],
    )
    return tmp_path


def test_surface_curv_and_annot(subject, engine):
    read, _ = engine
    mesh = read(
        subject / "lh.pial",
        curv=[subject / "lh.thickness"],
        annot=subject / "lh.aparc.annot",
    )
    assert np.allclose(mesh.points, PTS)
    assert np.asarray(mesh.cells[0].data).tolist() == TRI.tolist()
    assert np.allclose(mesh.point_data["thickness"], [1, 2, 3, 4])
    assert np.asarray(mesh.point_data["freesurfer:annot"]).tolist() == [0, 0, 1, 1]
    regions = {r.name: (np.asarray(r.entries).tolist(), r.tag) for r in mesh.regions}
    assert regions == {"bankssts": ([0, 1], 0), "frontal": ([2, 3], 1)}
    assert bytes(mesh.field_data["freesurfer:created_by"]) == b"created by test"


def test_write_round_trips_through_nibabel_with_equal_bytes(subject, engine, tmp_path):
    read, write = engine
    mesh = read(subject / "lh.pial")
    out = tmp_path / "out.pial"
    write(out, mesh)
    assert out.read_bytes() == (subject / "lh.pial").read_bytes()
    coords, faces = fsio.read_geometry(out)
    assert np.allclose(coords, PTS) and faces.tolist() == TRI.tolist()


def test_errors(tmp_path, engine):
    read, write = engine
    quad = tmp_path / "q.pial"
    quad.write_bytes(b"\xff\xff\xff" + b"\0" * 12)
    with pytest.raises(mio.ReadError, match="quad surface"):
        read(quad)
    with pytest.raises(mio.WriteError, match="triangle"):
        write(
            tmp_path / "x.pial",
            mio.Mesh(np.zeros((4, 3)), [("quad", np.array([[0, 1, 2, 3]]))]),
        )

"""3MF (roadmap 2.4.5): packages built here, both engines, lib3mf both ways."""

import zipfile

import numpy as np
import pytest

import meshioplusplus as mio
from meshioplusplus import _core, _fallback
from meshioplusplus._regions import Region
from meshioplusplus.threemf import _threemf

HAS_ZLIB = getattr(_core, "__has_zlib__", False)

MODEL = """<?xml version="1.0" encoding="UTF-8"?>
<model unit="inch" xmlns="http://schemas.microsoft.com/3dmanufacturing/core/2015/02">
<resources>
<basematerials id="1"><base name="red" displaycolor="#FF0000"/>
<base name="blue" displaycolor="#0000FF80"/></basematerials>
<object id="2" type="model" name="tri" pid="1" pindex="0"><mesh><vertices>
<vertex x="0" y="0" z="0"/><vertex x="1" y="0" z="0"/><vertex x="0" y="1" z="0"/>
</vertices><triangles><triangle v1="0" v2="1" v3="2"/>
<triangle v1="0" v2="2" v3="1" pid="1" p1="1"/></triangles></mesh></object>
<object id="3" type="model"><components>
<component objectid="2" transform="1 0 0 0 1 0 0 0 1 10 0 0"/></components></object>
</resources>
<build><item objectid="2"/><item objectid="3" transform="2 0 0 0 2 0 0 0 2 0 0 5"/></build>
</model>
"""


def _package(path, model=MODEL, rels=True, method=zipfile.ZIP_DEFLATED):
    with zipfile.ZipFile(path, "w", method) as z:
        z.writestr("[Content_Types].xml", "<Types/>")
        if rels:
            z.writestr(
                "_rels/.rels",
                '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/'
                'relationships"><Relationship Target="/3D/3dmodel.model" Id="r" '
                f'Type="{_threemf.MODEL_TYPE}"/></Relationships>',
            )
        z.writestr("3D/3dmodel.model", model)


@pytest.fixture(params=["core", "python"])
def engine(request, monkeypatch):
    if request.param == "core":
        if not HAS_ZLIB:
            pytest.skip("the core reads deflated packages only with zlib")
        monkeypatch.setattr(_fallback, "_strict", True)
        return mio.threemf.read, mio.threemf.write
    return _threemf.read, _threemf.write


def test_build_components_transforms_and_materials(tmp_path, engine):
    read, _ = engine
    path = tmp_path / "a.3mf"
    _package(path)
    mesh = read(path)
    assert len(mesh.points) == 6
    # The second item: the component shifts x by 10, then the item scales by 2 and
    # lifts z by 5.
    assert np.allclose(mesh.points[3:], [[20, 0, 5], [22, 0, 5], [20, 2, 5]])
    assert np.asarray(mesh.cells[0].data).tolist() == [
        [0, 1, 2],
        [0, 2, 1],
        [3, 4, 5],
        [3, 5, 4],
    ]
    assert {r.name: np.asarray(r.entries).tolist() for r in mesh.regions} == {
        "tri": [0, 1],
        "object_3": [2, 3],
    }
    colours = np.asarray(mesh.cell_data["3mf:color"][0]).tolist()
    assert colours[0] == [255, 0, 0, 255] and colours[1] == [0, 0, 255, 128]
    assert np.allclose(mesh.field_data["meshio:length_unit"], [0.0254])


def test_required_extensions_are_refused(tmp_path, engine):
    read, _ = engine
    path = tmp_path / "x.3mf"
    _package(
        path, MODEL.replace("<model ", '<model requiredextensions="p" '), rels=False
    )
    with pytest.raises(mio.ReadError, match="requires the extensions"):
        read(path)


def test_write_round_trip_with_equal_bytes(tmp_path, engine):
    read, write = engine
    src = tmp_path / "a.3mf"
    _package(src)
    mesh = mio.threemf.read(src)
    out = tmp_path / "out.3mf"
    write(out, mesh)
    other = tmp_path / "other.3mf"
    (_threemf.write if write is mio.threemf.write else mio.threemf.write)(other, mesh)
    if HAS_ZLIB:  # without zlib the core stores while Python deflates
        assert out.read_bytes() == other.read_bytes()
    back = read(out)
    # Objects are written in region-name order, so the points come back regrouped.
    assert np.allclose(np.sort(back.points, axis=0), np.sort(mesh.points, axis=0))
    assert {r.name for r in back.regions} == {"tri", "object_3"}
    assert np.array_equal(
        back.cell_data["3mf:color"][0], mesh.cell_data["3mf:color"][0]
    )
    assert np.allclose(back.field_data["meshio:length_unit"], [0.0254])


def test_non_triangles_are_refused(tmp_path, engine):
    _, write = engine
    quads = mio.Mesh(np.zeros((4, 3)), [("quad", np.array([[0, 1, 2, 3]]))])
    with pytest.raises(mio.WriteError, match="triangles"):
        write(tmp_path / "q.3mf", quads)


def test_lib3mf_reads_what_we_write_and_the_reverse(tmp_path):
    lib3mf = pytest.importorskip("lib3mf")
    w = lib3mf.get_wrapper()
    mesh = mio.Mesh(
        np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float),
        [("triangle", np.array([[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]]))],
        regions=[Region("half", "cell", [0, 1])],
    )
    ours = tmp_path / "ours.3mf"
    mio.write(ours, mesh)
    model = w.CreateModel()
    model.QueryReader("3mf").ReadFromFile(str(ours))
    it = model.GetMeshObjects()
    counts = []
    while it.MoveNext():
        obj = it.GetCurrentMeshObject()
        counts.append((obj.GetName(), obj.GetVertexCount(), obj.GetTriangleCount()))
    assert sorted(counts) == [("half", 4, 2), ("mesh", 4, 2)]
    theirs = tmp_path / "theirs.3mf"
    model.QueryWriter("3mf").WriteToFile(str(theirs))
    for read in (mio.threemf.read, _threemf.read):
        back = read(theirs)
        assert len(back.cells[0].data) == 4
        assert {r.name for r in back.regions} == {"half", "mesh"}


def test_core_reads_a_stored_package_without_zlib(tmp_path, monkeypatch):
    monkeypatch.setattr(_fallback, "_strict", True)
    path = tmp_path / "s.3mf"
    _package(path, method=zipfile.ZIP_STORED)
    assert len(mio.threemf.read(path).points) == 6

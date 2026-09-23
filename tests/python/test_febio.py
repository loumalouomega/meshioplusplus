"""FEBio input files (``.feb``): both engines, the fixtures written by
``tools/gen_febio_fixtures.py`` in FEBio's own node numbering (spec 2.5, 3.0
and 4.0), and the spec-4.0 writer."""

import pathlib

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _core
from meshioplusplus.febio import _febio as py_febio

from .test_node_order import _is_valid

MESHES = pathlib.Path(__file__).parent / "meshes" / "febio"
FIXTURES = sorted(MESHES.glob("*.feb"))


@pytest.fixture(params=["core", "python"])
def engine(request):
    """Both engines behind the same read/write signatures."""
    if request.param == "core":
        return meshioplusplus.febio
    return py_febio


def _regions(mesh):
    return sorted(
        (r.kind, r.name, r.dim, r.tag, tuple(np.asarray(r.entries).ravel()))
        for r in mesh.regions
    )


def _same(a, b):
    np.testing.assert_array_equal(a.points, b.points)
    assert [c.type for c in a.cells] == [c.type for c in b.cells]
    for x, y in zip(a.cells, b.cells):
        np.testing.assert_array_equal(x.data, y.data)
    assert _regions(a) == _regions(b)
    assert sorted(a.point_data) == sorted(b.point_data)
    for name in a.point_data:
        np.testing.assert_array_equal(a.point_data[name], b.point_data[name])
    assert sorted(a.cell_data) == sorted(b.cell_data)
    for name in a.cell_data:
        for x, y in zip(a.cell_data[name], b.cell_data[name]):
            np.testing.assert_array_equal(x, y)


@pytest.mark.parametrize("path", FIXTURES, ids=[p.name for p in FIXTURES])
def test_engines_read_the_same_mesh(path):
    _same(meshioplusplus.febio.read(path), py_febio.read(path))


def test_every_solid_is_valid_in_meshio_order():
    """Mid-edge nodes on midpoints, hex27 face centres on their faces, positive
    volumes: the FEBio node order was undone correctly."""
    mesh = meshioplusplus.read(MESHES / "all_elements_v40.feb")
    solids = [b for b in mesh.cells if b.type.rstrip("0123456789") != "line"]
    solids = [b for b in solids if b.type not in ("quad", "triangle")]
    assert [b.type for b in solids] == [
        "tetra",
        "tetra10",
        "wedge",
        "wedge15",
        "pyramid",
        "pyramid13",
        "hexahedron",
        "hexahedron20",
        "hexahedron27",
    ]
    for block in solids:
        for row in block.data:
            assert _is_valid(mesh.points[row], block.type), block.type


def test_what_each_mesh_section_becomes(engine):
    mesh = engine.read(MESHES / "all_elements_v40.feb")
    got = {(r.name, r.kind): r for r in mesh.regions}
    # Every <Elements> block is a cell region tagged with its domain's material.
    assert got[("part_hex27", "cell")].tag == 9
    assert got[("part_quad4", "cell")].dim == 2
    # Named <Nodes> blocks are node sets, as in FEBio; ranges expand.
    assert len(got[("SolidNodes", "point")].entries) == 108
    assert got[("first_ten", "point")].entries.tolist() == list(range(10))
    # A surface on solid faces is a side region; one off them its own block.
    assert got[("hex_base", "side")].entries.tolist() == [[6, 4]]
    floating = got[("floating", "cell")]
    assert mesh.cells[10].type == "triangle" and floating.entries.tolist() == [10]
    assert [mesh.cells[k].type for k in (11, 12)] == ["line", "line"]
    assert got[("hex_edge", "cell")].entries.tolist() == [11]
    assert got[("springs", "cell")].entries.tolist() == [12]
    # MeshData: NaN outside the set the array is defined on.
    temperature = mesh.point_data["temperature"]
    hex8 = mesh.cells[6].data[0]
    assert temperature[hex8[:4]].tolist() == [10, 20, 30, 40]
    assert np.isnan(np.delete(temperature, hex8[:4])).all()
    assert mesh.cell_data["fiber"][6].tolist() == [[1, 0, 0]]
    assert np.isnan(mesh.cell_data["fiber"][0]).all()


@pytest.mark.parametrize("version", ["v25", "v30"])
def test_older_specs(engine, version):
    mesh = engine.read(MESHES / f"block_{version}.feb")
    got = {r.name: r for r in mesh.regions}
    assert (got["left"].tag, got["right"].tag) == (1, 2)
    assert got["fixed"].entries.tolist() == [0, 1, 2, 3]
    assert got["both"].entries.tolist() == [0, 1]
    assert got["top"].kind == "side" and len(got["top"].entries) == 2
    assert [a.tolist() for a in mesh.cell_data["thickness"]] == [[0.5], [1.0]]


def test_writers_write_the_same_bytes(tmp_path):
    for path in FIXTURES:
        mesh = meshioplusplus.read(path)
        mesh.point_data.clear()
        mesh.cell_data.clear()
        _core.febio_write(str(tmp_path / "cpp.feb"), mesh)
        py_febio.write(tmp_path / "py.feb", mesh)
        assert (tmp_path / "cpp.feb").read_bytes() == (tmp_path / "py.feb").read_bytes()


@pytest.mark.parametrize("path", FIXTURES, ids=[p.name for p in FIXTURES])
def test_round_trip(engine, path, tmp_path):
    mesh = meshioplusplus.read(path)
    mesh.point_data.clear()
    mesh.cell_data.clear()
    engine.write(tmp_path / "out.feb", mesh)
    back = engine.read(tmp_path / "out.feb")
    np.testing.assert_allclose(back.points, mesh.points)
    assert [(c.type, c.data.tolist()) for c in back.cells] == [
        (c.type, c.data.tolist()) for c in mesh.cells
    ]
    # Names and memberships survive; the tags (material ids) are renumbered.
    before = {(r.name, r.kind): np.asarray(r.entries).tolist() for r in mesh.regions}
    after = {(r.name, r.kind): np.asarray(r.entries).tolist() for r in back.regions}
    assert after == before


def test_what_the_writer_makes_of_foreign_blocks(tmp_path):
    """Boundary triangles on a tet become a <Surface> (read back as a side
    region); lines on its edges an <Edge>; a lone two-node line a
    <DiscreteSet>; vertices are dropped."""
    points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [3, 3, 3]], float)
    mesh = meshioplusplus.Mesh(
        points,
        [
            ("tetra", np.array([[0, 1, 2, 3]])),
            ("triangle", np.array([[0, 2, 1]])),
            ("line", np.array([[0, 1]])),
            ("line", np.array([[3, 4]])),
            ("vertex", np.array([[4]])),
        ],
    )
    meshioplusplus.write(tmp_path / "out.feb", mesh)
    text = (tmp_path / "out.feb").read_text()
    assert '<Surface name="Surface2">' in text
    assert '<Edge name="Edge3">' in text
    assert '<DiscreteSet name="DiscreteSet4">' in text
    assert "<delem>4,5</delem>" in text
    assert '<SolidDomain name="Part1" mat="Part1"/>' in text
    back = meshioplusplus.read(tmp_path / "out.feb")
    assert [c.type for c in back.cells] == ["tetra", "line", "line"]
    assert {(r.name, r.kind) for r in back.regions} == {
        ("Part1", "cell"),
        ("Surface2", "side"),
        ("Edge3", "cell"),
        ("DiscreteSet4", "cell"),
    }


def test_errors(engine, tmp_path):
    bad = tmp_path / "bad.feb"
    bad.write_text('<febio_spec version="1.2"><Geometry/></febio_spec>')
    with pytest.raises(meshioplusplus.ReadError, match="1.2"):
        engine.read(bad)
    bad.write_text(
        '<febio_spec version="4.0"><Mesh><Nodes><node id="1">0,0,0</node></Nodes>'
        '<Elements type="tet4" name="p"><elem id="1">1,2,3,4</elem></Elements>'
        "</Mesh></febio_spec>"
    )
    with pytest.raises(meshioplusplus.ReadError, match="undefined node 2"):
        engine.read(bad)
    bad.write_text(
        '<febio_spec version="4.0"><Mesh><Nodes>'
        + "".join(f'<node id="{i}">{i},0,0</node>' for i in range(1, 6))
        + '</Nodes><Elements type="tet5" name="p"><elem id="1">1,2,3,4,5</elem>'
        "</Elements></Mesh></febio_spec>"
    )
    with pytest.raises(meshioplusplus.ReadError, match="tet5"):
        engine.read(bad)
    assert engine.read(bad, lenient=True).cells[0].type == "tetra"
    bad.write_text(
        '<febio_spec version="3.0"><Geometry><Part name="a"/></Geometry></febio_spec>'
    )
    with pytest.raises(meshioplusplus.ReadError, match="Part"):
        engine.read(bad)


def test_mesh_from_another_file(engine, tmp_path):
    (tmp_path / "driver.feb").write_text(
        '<febio_spec version="4.0"><Module type="solid"/>'
        f'<Mesh from="{MESHES / "block_v30.feb"}"/></febio_spec>'
    )
    # A 4.0 driver pulling the <Mesh> of a 3.0 file, as FEBio's `from=` does.
    mesh = engine.read(tmp_path / "driver.feb")
    assert [c.type for c in mesh.cells] == ["hexahedron", "hexahedron"]


def test_sniffed_and_registered():
    assert meshioplusplus.sniff_format(MESHES / "block_v25.feb") == "febio"
    out = meshioplusplus.formats()
    assert out["extensions"][".feb"] == ["febio"]

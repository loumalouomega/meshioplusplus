"""COMSOL native meshes, text (``.mphtxt``) and binary (``.mphbin``): both
engines, the fixtures written by ``tools/gen_comsol_fixtures.py`` in COMSOL's
own node numbering (not by meshio++), Selections as regions, Mesh versions 2
and 4, and the writer's entity indices and Selections.

The readers were also run over the COMSOL files published with deal.II, FEconv
and Wolfram's FEMAddOns (not redistributed here): both engines agree on every
one and put every quadratic mid-edge node on its edge's midpoint.
"""

import io
import pathlib

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _core
from meshioplusplus.mphtxt import _mphtxt as py_mphtxt

from . import helpers
from .test_node_order import _PYR, _PYR_E, _SHAPES, _is_valid

MESHES = pathlib.Path(__file__).resolve().parent / "meshes" / "comsol"
FIXTURES = [
    "two_domains.mphtxt",
    "quadratic.mphtxt",
    "legacy_v2.mphtxt",
    "two_objects.mphtxt",
]

MESHES_TO_TEST = [
    helpers.line_mesh,
    helpers.tri_mesh,
    helpers.tri_mesh_2d,
    helpers.triangle6_mesh,
    helpers.quad_mesh,
    helpers.tet_mesh,
    helpers.tet10_mesh,
    helpers.hex_mesh,
    helpers.wedge_mesh,
    helpers.tri_quad_mesh,
]


@pytest.mark.parametrize("mesh", MESHES_TO_TEST)
def test_io(mesh, tmp_path):
    helpers.write_read(
        tmp_path, meshioplusplus.mphtxt.write, meshioplusplus.mphtxt.read, mesh, 1.0e-12
    )


@pytest.mark.parametrize("mesh", MESHES_TO_TEST)
def test_io_binary(mesh, tmp_path):
    helpers.write_read(
        tmp_path,
        meshioplusplus.mphbin.write,
        meshioplusplus.mphbin.read,
        mesh,
        0.0,
        extension=".mphbin",
    )


def test_generic_io(tmp_path):
    helpers.generic_io(tmp_path / "test.mphtxt")
    helpers.generic_io(tmp_path / "test.0.mphtxt")
    helpers.generic_io(tmp_path / "test.mphbin")


ENGINES = {
    "core": (_core.mphtxt_read, _core.mphbin_read),
    "python": (py_mphtxt.read, py_mphtxt.read_binary),
}


@pytest.fixture(params=sorted(ENGINES))
def engine(request):
    return ENGINES[request.param]


def _regions(mesh):
    return {(r.kind, r.name, r.dim, r.tag): r.entries.tolist() for r in mesh.regions}


def test_two_domains(engine):
    # the done-when of roadmap §1.2: a mixed tet/prism file with two domains,
    # read with meshio++'s node order and the domain ids
    for read, name in zip(engine, ("two_domains.mphtxt", "two_domains.mphbin")):
        mesh = read(str(MESHES / name))
        assert [c.type for c in mesh.cells] == [
            "vertex",
            "line",
            "triangle",
            "quad",
            "tetra",
            "wedge",
        ]
        geom = [g.tolist() for g in mesh.cell_data["mphtxt:geom"]]
        assert geom == [[0, 3], [0, 1], [0, 2, 1], [3], [1], [2]]
        # quad: COMSOL's tensor order becomes a ring
        assert mesh.cells[3].data.tolist() == [[1, 2, 5, 4]]
        assert _regions(mesh) == {
            ("cell", "Lower Part", 3, -1): [8],
            ("cell", "Part #2", 3, -1): [9],  # '#' inside a string is no comment
            ("cell", "caps", 2, -1): [4, 6],
        }
        for block in mesh.cells[4:]:
            assert _is_valid(mesh.points[block.data[0]], block.type)


def test_quadratic_node_order(engine):
    mesh = engine[0](str(MESHES / "quadratic.mphtxt"))
    shapes = dict(_SHAPES)
    shapes["pyramid14"] = (_PYR, _PYR_E, [(0, 1, 2, 3)], False)
    types = [c.type for c in mesh.cells]
    assert types == [
        "tetra10",
        "pyramid14",
        "wedge18",
        "hexahedron27",
        "triangle6",
        "quad9",
        "line3",
    ]
    for block in mesh.cells:
        corners, edges, faces, body = shapes[block.type]
        x = mesh.points[block.data[0]]
        c = x[: len(corners)]
        k = len(corners)
        for a, b in edges:
            np.testing.assert_allclose(x[k], (c[a] + c[b]) / 2)
            k += 1
        for f in faces:
            np.testing.assert_allclose(x[k], c[list(f)].mean(0))
            k += 1
        if body:
            np.testing.assert_allclose(x[k], c.mean(0))
        if block.type not in ("triangle6", "quad9", "line3", "pyramid14"):
            assert _is_valid(x, block.type)


def test_legacy_version_2(engine):
    # parameters with three values per node on faces and one on edges, up/down
    # pairs, lowest vertex index 1
    mesh = engine[0](str(MESHES / "legacy_v2.mphtxt"))
    assert [c.type for c in mesh.cells] == ["vertex", "line3", "quad9", "hexahedron27"]
    assert [g.tolist() for g in mesh.cell_data["mphtxt:geom"]] == [
        [0],
        [0],
        [4, 5],
        [1],
    ]
    x = mesh.points[mesh.cells[3].data[0]]
    np.testing.assert_allclose(x[26], [0.5, 0.5, 0.5])  # body centre last


def test_several_mesh_objects(engine):
    mesh = engine[0](str(MESHES / "two_objects.mphtxt"))
    assert len(mesh.points) == 7
    assert mesh.cells[1].data.tolist() == [[3, 4, 5], [4, 6, 5]]
    assert _regions(mesh) == {
        ("cell", "mesh1", -1, -1): [0],
        ("cell", "mesh2", -1, -1): [1, 2],
        ("cell", "right", 2, -1): [2],  # a selection of the second object
    }


@pytest.mark.parametrize("name", FIXTURES)
def test_engines_agree(name):
    a = _core.mphtxt_read(str(MESHES / name))
    b = py_mphtxt.read(str(MESHES / name))
    assert np.array_equal(a.points, b.points)
    assert [(c.type, c.data.tolist()) for c in a.cells] == [
        (c.type, c.data.tolist()) for c in b.cells
    ]
    assert [g.tolist() for g in a.cell_data["mphtxt:geom"]] == [
        g.tolist() for g in b.cell_data["mphtxt:geom"]
    ]
    assert _regions(a) == _regions(b)


def test_text_and_binary_agree():
    a = _core.mphtxt_read(str(MESHES / "two_domains.mphtxt"))
    b = _core.mphbin_read(str(MESHES / "two_domains.mphbin"))
    assert np.array_equal(a.points, b.points)
    assert _regions(a) == _regions(b)


@pytest.mark.parametrize("name", FIXTURES)
def test_writers_emit_the_same_bytes_and_round_trip(name, tmp_path):
    mesh = _core.mphtxt_read(str(MESHES / name))
    for ext, core_write, py_write, read in [
        (".mphtxt", _core.mphtxt_write, py_mphtxt.write, _core.mphtxt_read),
        (".mphbin", _core.mphbin_write, py_mphtxt.write_binary, _core.mphbin_read),
    ]:
        c, p = tmp_path / ("c" + ext), tmp_path / ("p" + ext)
        core_write(str(c), mesh)
        py_write(str(p), mesh)
        assert c.read_bytes() == p.read_bytes()
        back = read(str(c))
        assert [(x.type, x.data.tolist()) for x in back.cells] == [
            (x.type, x.data.tolist()) for x in mesh.cells
        ]
        if name != "two_objects.mphtxt":  # object regions are not entity unions
            assert _regions(back) == _regions(mesh)


@pytest.mark.parametrize(
    "write", [_core.mphtxt_write, py_mphtxt.write], ids=["core", "python"]
)
def test_writer_derives_entities_from_regions(write, tmp_path, capfd):
    tet = helpers.tet_mesh
    Region = meshioplusplus.Region
    mesh = meshioplusplus.Mesh(
        tet.points,
        tet.cells,
        regions=[
            Region("a", "cell", [0], dim=3),
            Region("b", "cell", [1], dim=3),
            Region("both", "cell", [0, 1], dim=3),  # a union of whole entities
            Region("pts", "point", [0]),
        ],
    )
    path = tmp_path / "out.mphtxt"
    write(str(path), mesh)
    back = meshioplusplus.mphtxt.read(str(path))
    # domains count from 1 in region order
    assert back.cell_data["mphtxt:geom"][0].tolist() == [1, 2]
    assert _regions(back) == {
        ("cell", "a", 3, -1): [0],
        ("cell", "b", 3, -1): [1],
        ("cell", "both", 3, -1): [0, 1],
    }
    assert "'pts' dropped" in " ".join(capfd.readouterr().err.split())

    # a region splitting an entity cannot be a Selection
    mesh.regions = [Region("half", "cell", [0], dim=3)]
    mesh.cell_data["mphtxt:geom"] = [np.array([4, 4])]
    write(str(path), mesh)
    assert meshioplusplus.mphtxt.read(str(path)).regions == []
    assert "'half' is not a union" in " ".join(capfd.readouterr().err.split())


def test_no_regions_default_entities(tmp_path):
    # without regions or mphtxt:geom: domains are 1, lower dimensions 0
    mesh = meshioplusplus.Mesh(
        helpers.tet_mesh.points,
        [("tetra", helpers.tet_mesh.cells[0].data), ("triangle", [[0, 1, 2]])],
    )
    path = tmp_path / "out.mphtxt"
    meshioplusplus.mphtxt.write(path, mesh)
    geom = meshioplusplus.mphtxt.read(path).cell_data["mphtxt:geom"]
    assert [g.tolist() for g in geom] == [[1, 1], [0]]
    assert "4 # version" in path.read_text()


def test_buffers():
    text = (MESHES / "two_domains.mphtxt").read_text()
    mesh = py_mphtxt.read(io.StringIO(text))
    assert len(mesh.regions) == 3
    out = io.StringIO()
    py_mphtxt.write(out, mesh)
    assert py_mphtxt.read(io.StringIO(out.getvalue())).regions == mesh.regions
    raw = (MESHES / "two_domains.mphbin").read_bytes()
    binary = py_mphtxt.read_binary(io.BytesIO(raw))
    assert binary.regions == mesh.regions


@pytest.mark.parametrize(
    "read", [_core.mphtxt_read, py_mphtxt.read], ids=["core", "python"]
)
def test_errors(read, tmp_path, capfd):
    path = tmp_path / "bad.mphtxt"
    path.write_text("0 2\n0\n0\n")
    with pytest.raises(meshioplusplus.ReadError, match="version"):
        read(str(path))
    path.write_text("0 1\n1 5 geom1\n1 3 obj\n0 0 1\n5 Geom3\n")
    with pytest.raises(meshioplusplus.ReadError, match="no Mesh object"):
        read(str(path))
    assert "'Geom3'" in capfd.readouterr().err
    path.write_text(
        "0 1\n1 5 mesh1\n1 3 obj\n0 0 1\n4 Mesh\n4\n2\n3\n0\n0 0\n1 0\n0 1\n"
        "1\n3 tri\n3\n1\n0 1 7\n1\n1\n"
    )
    with pytest.raises(meshioplusplus.ReadError, match="vertex 7"):
        read(str(path))


def test_mphbin_is_not_given_a_provenance_slot(tmp_path):
    from meshioplusplus import _provenance

    path = tmp_path / "out.mphbin"
    with _provenance.scope(_provenance.Mode.REQUIRED):
        with pytest.raises(meshioplusplus.WriteError):
            py_mphtxt.write_binary(str(path), helpers.tri_mesh)

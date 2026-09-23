"""Code_Aster native mesh (``.mail``): both engines, the fixtures written by
``tools/gen_code_aster_fixtures.py`` in Code_Aster's own node numbering, and a
round trip through MED."""

import pathlib

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _core
from meshioplusplus.code_aster import _code_aster as py_code_aster

from .test_node_order import _is_valid

MESHES = pathlib.Path(__file__).parent / "meshes" / "code_aster"
FIXTURES = sorted(MESHES.glob("*.mail"))

# meshio++ (VTK) mid-edge corner pairs, for the geometry checks.
_EDGES = {
    "line3": [(0, 1)],
    "triangle6": [(0, 1), (1, 2), (2, 0)],
    "triangle7": [(0, 1), (1, 2), (2, 0)],
    "quad8": [(0, 1), (1, 2), (2, 3), (3, 0)],
    "quad9": [(0, 1), (1, 2), (2, 3), (3, 0)],
    "tetra10": [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)],
    "pyramid13": [(0, 1), (1, 2), (2, 3), (3, 0), (0, 4), (1, 4), (2, 4), (3, 4)],
    "wedge15": [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (0, 3), (1, 4), (2, 5)],
    "wedge18": [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (0, 3), (1, 4), (2, 5)],
    "hexahedron20": [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4)]
    + [(0, 4), (1, 5), (2, 6), (3, 7)],
    "hexahedron27": [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4)]
    + [(0, 4), (1, 5), (2, 6), (3, 7)],
}


_CORNERS = {
    "line3": 2,
    "triangle6": 3,
    "triangle7": 3,
    "quad8": 4,
    "quad9": 4,
    "tetra10": 4,
    "pyramid13": 5,
    "wedge15": 6,
    "wedge18": 6,
    "hexahedron20": 8,
    "hexahedron27": 8,
}


@pytest.fixture(params=["core", "python"])
def engine(request):
    """Both engines behind the same read/write signatures."""
    if request.param == "core":
        return meshioplusplus.code_aster
    return py_code_aster


def _blocks(mesh):
    return [(b.type, np.asarray(b.data)) for b in mesh.cells]


def _regions(mesh):
    return sorted(
        (r.kind, r.name, r.dim, tuple(np.asarray(r.entries).ravel()))
        for r in mesh.regions
    )


def _same(a, b):
    np.testing.assert_array_equal(a.points, b.points)
    assert [t for t, _ in _blocks(a)] == [t for t, _ in _blocks(b)]
    for (_, x), (_, y) in zip(_blocks(a), _blocks(b)):
        np.testing.assert_array_equal(x, y)
    assert _regions(a) == _regions(b)


def _mixed_mesh():
    """Every cell type the format holds, one cell each, plus regions."""
    rows = []
    points = []
    for k, (keyword, (cell_type, count)) in enumerate(py_code_aster._TYPES.items()):
        base = len(points)
        points += [[float(k), float(j), 0.5 * j] for j in range(count)]
        rows.append((cell_type, np.arange(base, base + count)[None]))
    mesh = meshioplusplus.Mesh(np.array(points), rows)
    mesh.regions = [
        meshioplusplus.Region("SOLIDS", "cell", np.arange(10, 20), dim=3),
        meshioplusplus.Region("FIRST", "point", [0, 1, 2]),
        meshioplusplus.Region("FIRST", "cell", [0]),
    ]
    return mesh


@pytest.mark.parametrize("path", FIXTURES, ids=[p.name for p in FIXTURES])
def test_engines_agree_on_every_fixture(path):
    _same(meshioplusplus.code_aster.read(path), py_code_aster.read(path))


@pytest.mark.parametrize("path", FIXTURES, ids=[p.name for p in FIXTURES])
def test_fixture_cells_are_valid_in_meshio_order(engine, path):
    """The fixtures list nodes in Code_Aster's numbering; after reading, every
    mid-edge node sits on its meshio++ edge's midpoint and the volume cells are
    valid, positively oriented elements."""
    mesh = engine.read(path)
    points = mesh.points
    if points.shape[1] == 2:
        points = np.hstack([points, np.zeros((len(points), 1))])
    checked = 0
    for cell_type, data in _blocks(mesh):
        edges = _EDGES.get(cell_type, [])
        corners = _CORNERS.get(cell_type, 0)
        for row in data:
            x = points[row]
            for k, (a, b) in enumerate(edges):
                np.testing.assert_allclose(
                    x[corners + k], (x[a] + x[b]) / 2, err_msg=cell_type
                )
                checked += 1
            if meshioplusplus.Mesh(points, [(cell_type, [row])]).cells[0].dim == 3:
                assert _is_valid(x, cell_type), cell_type
    assert checked > 0


def test_hexa20_block_groups():
    mesh = meshioplusplus.read(MESHES / "hexa20_block.mail")
    regions = {(r.kind, r.name): r for r in mesh.regions}
    assert set(regions) == {("cell", "VOLUME"), ("cell", "X0"), ("point", "X0")}
    np.testing.assert_array_equal(regions[("cell", "VOLUME")].entries, [0, 1])
    assert regions[("cell", "VOLUME")].dim == 3
    assert regions[("cell", "X0")].dim == 2
    np.testing.assert_allclose(mesh.points[regions[("point", "X0")].entries][:, 0], 0.0)


def test_plate_2d_is_two_dimensional_and_warns_on_a_repeated_member(engine, capfd):
    mesh = engine.read(MESHES / "plate_2d.mail")
    assert "more than once" in capfd.readouterr().err
    assert mesh.points.shape == (11, 2)
    corners = [r for r in mesh.regions if r.name == "CORNERS"][0]
    np.testing.assert_array_equal(corners.entries, [0, 4])


def test_engines_write_the_same_bytes(tmp_path):
    mesh = _mixed_mesh()
    _core.code_aster_write(str(tmp_path / "cpp.mail"), mesh)
    py_code_aster.write(tmp_path / "py.mail", mesh)
    cpp = (tmp_path / "cpp.mail").read_bytes()
    assert cpp == (tmp_path / "py.mail").read_bytes()
    assert b"Written by meshio++" in cpp
    assert all(len(line) <= 80 for line in cpp.split(b"\n"))


def test_round_trip_keeps_cells_and_regions(engine, tmp_path):
    mesh = _mixed_mesh()
    engine.write(tmp_path / "m.mail", mesh)
    back = engine.read(tmp_path / "m.mail")
    np.testing.assert_array_equal(back.points, mesh.points)
    assert [b.type for b in back.cells] == [b.type for b in mesh.cells]
    for x, y in zip(back.cells, mesh.cells):
        np.testing.assert_array_equal(x.data, y.data)
    assert [
        (r.kind, r.name, list(r.entries))
        for r in sorted(back.regions, key=lambda r: r.key)
    ] == [
        (r.kind, r.name, list(r.entries))
        for r in sorted(mesh.regions, key=lambda r: r.key)
    ]


def test_hexa20_through_med_and_back(tmp_path):
    """The roadmap's done-when: a HEXA20 .mail through MED and back keeps its
    connectivity and groups."""
    pytest.importorskip("h5py")
    mesh = meshioplusplus.read(MESHES / "hexa20_block.mail")
    meshioplusplus.write(tmp_path / "a.mail", mesh)
    meshioplusplus.write(tmp_path / "b.med", mesh)
    via_med = meshioplusplus.read(tmp_path / "b.med")
    via_med.point_data.clear()
    via_med.cell_data.clear()
    meshioplusplus.write(tmp_path / "c.mail", via_med)

    def body(p):
        return [line for line in p.read_text().splitlines() if not line.startswith("%")]

    assert body(tmp_path / "a.mail") == body(tmp_path / "c.mail")


def test_group_names_are_sanitised_the_same_way(tmp_path, capfd):
    mesh = _mixed_mesh()
    mesh.regions = [
        meshioplusplus.Region("a group name longer than twenty-four", "cell", [0]),
        meshioplusplus.Region("a group name longer than twenty-four!", "cell", [1]),
        meshioplusplus.Region("tête", "point", [0]),
    ]
    py_code_aster.write(tmp_path / "py.mail", mesh)
    assert "is written as" in capfd.readouterr().err
    _core.code_aster_write(str(tmp_path / "cpp.mail"), mesh)
    assert (tmp_path / "cpp.mail").read_bytes() == (tmp_path / "py.mail").read_bytes()
    names = sorted(r.name for r in meshioplusplus.read(tmp_path / "cpp.mail").regions)
    assert names == ["a_group_name_longer_th_1", "a_group_name_longer_than", "t__te"]


def test_side_regions_and_data_are_dropped_with_a_warning(engine, tmp_path, capfd):
    mesh = _mixed_mesh()
    mesh.point_data["T"] = np.zeros(len(mesh.points))
    mesh.regions.append(meshioplusplus.Region("WALL", "side", [[10, 0]], dim=2))
    engine.write(tmp_path / "m.mail", mesh)
    err = capfd.readouterr().err
    assert "side region" in err and "no data arrays" in err
    text = (tmp_path / "m.mail").read_text()
    assert "WALL" not in text
    assert "data-dropped" in text and "regions-dropped" in text


@pytest.mark.parametrize(
    "body, match",
    [
        ("COOR_3D\n N1 0 0 0\nFINSF\nSEG2\n M1 N1 N9\nFINSF\nFIN\n", "N9"),
        ("COOR_3D\n N1 0 0 0\n N1 1 0 0\nFINSF\nFIN\n", "defined twice"),
        ("COOR_3D\n N1 0 0 0\nFINSF\nGROUP_NO\n G N2\nFINSF\nFIN\n", "N2"),
        ("COOR_3D\n N1 0 0 inf\nFINSF\nFIN\n", "inf"),
        ("COOR_3D\n N1 0 0 0\n", "ends inside a block"),
        ("COOR_3D\n N1 0 0 0\nFINSF\nSEG2\n M1 N1\nFINSF\nFIN\n", "fewer than 2"),
        ("COOR_3D\n N1 0 0 0\nFINSF\n12 13\nFIN\n", "expected a keyword"),
    ],
)
def test_errors_name_the_culprit(engine, tmp_path, body, match):
    path = tmp_path / "bad.mail"
    path.write_text(body)
    with pytest.raises(meshioplusplus.ReadError, match=match):
        engine.read(path)


def test_unsupported_keywords_are_skipped(engine, tmp_path, capfd):
    path = tmp_path / "m.mail"
    path.write_text(
        "COOR_3D\n N1 0 0 0\n N2 1 0 0\nFINSF\n"
        "SEG22\n M9 N1 N2 N1 N2\nFINSF\nSYS_COOR\n anything\nFINSF\n"
        "SEG2\n M1 N1 N2\nFINSF\nFIN\n"
    )
    mesh = engine.read(path)
    assert "SEG22" in capfd.readouterr().err
    assert [b.type for b in mesh.cells] == ["line"]


def test_columns_past_80_are_ignored(engine, tmp_path, capfd):
    path = tmp_path / "m.mail"
    path.write_text(
        "COOR_3D\n N1 0 0 0\n N2 1 0 0\n N3 0 1 0\nFINSF\n"
        "SEG2\n M1 N1 N2\n" + " " * 80 + " M2 N2 N3\nFINSF\nFIN\n"
    )
    mesh = engine.read(path)
    assert "80 columns" in capfd.readouterr().err
    assert len(mesh.cells[0].data) == 1


def test_extension_and_sniffing(tmp_path):
    assert meshioplusplus._helpers._filetypes_from_path(pathlib.Path("x.mail")) == [
        "code_aster"
    ]
    path = tmp_path / "no_extension"
    path.write_bytes((MESHES / "plate_2d.mail").read_bytes())
    assert meshioplusplus.sniff_format(path) == "code_aster"

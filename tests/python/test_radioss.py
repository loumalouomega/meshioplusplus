"""OpenRadioss starter decks (``*_0000.rad``): both engines, the decks written by
``tools/gen_radioss_fixtures.py``, groups, subsets, surfaces and includes."""

import pathlib
import shutil

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _core
from meshioplusplus._facets import facet_nodes
from meshioplusplus.radioss import _radioss as py_radioss

from .test_node_order import _is_valid

MESHES = pathlib.Path(__file__).parent / "meshes" / "radioss"
DECKS = sorted(MESHES.glob("*_0000.rad"))


@pytest.fixture(params=["core", "python"])
def read(request):
    if request.param == "core":
        return meshioplusplus.radioss.read
    return py_radioss.read


def _regions(mesh):
    return sorted(
        (r.kind, r.name, r.dim, r.tag, tuple(np.asarray(r.entries).ravel()))
        for r in mesh.regions
    )


@pytest.mark.parametrize("path", DECKS, ids=[p.name for p in DECKS])
def test_engines_agree_on_every_deck(path):
    a = meshioplusplus.radioss.read(path)
    b = py_radioss.read(path)
    np.testing.assert_array_equal(a.points, b.points)
    assert [c.type for c in a.cells] == [c.type for c in b.cells]
    for x, y in zip(a.cells, b.cells):
        np.testing.assert_array_equal(x.data, y.data)
    for name in ("radioss:part", "radioss:property", "radioss:material"):
        for x, y in zip(a.cell_data[name], b.cell_data[name]):
            np.testing.assert_array_equal(x, y)
    assert _regions(a) == _regions(b)


def test_every_card(read):
    mesh = read(MESHES / "deck_0000.rad")
    counts = {c.type: len(c.data) for c in mesh.cells}
    assert counts == {
        "line": 4,  # beam, truss, spring, the included truss (the single-node spring goes)
        "hexahedron": 1,
        "wedge": 3,  # two degenerate bricks and the /PENTA6
        "tetra": 2,  # a degenerate brick and the inverted /TETRA4
        "pyramid": 1,
        "hexahedron20": 1,
        "tetra10": 1,
        "quad": 1,
        "triangle": 2,  # the /SHELL with n3 == n4 and the /SH3N
    }
    for block in mesh.cells:
        for row in np.asarray(block.data):
            assert _is_valid(np.asarray(mesh.points)[row], block.type), block.type
    # Node 99 after /END is not read; node 1 is the comma-separated line.
    assert len(mesh.points) == 88
    np.testing.assert_array_equal(mesh.points[0], [0.0, 0.0, 0.0])
    assert int(mesh.field_data["radioss:version"]) == 2019


def _region(mesh, kind, name):
    (r,) = [r for r in mesh.regions if r.kind == kind and r.name == name]
    return r


def _types(mesh, cells):
    types = [c.type for c in mesh.cells for _ in range(len(c.data))]
    return sorted(types[c] for c in cells)


def test_parts_and_subsets(read):
    mesh = read(MESHES / "deck_0000.rad")
    bricks = _region(mesh, "cell", "bricks")
    assert bricks.tag == 1 and bricks.dim == 3
    assert _types(mesh, bricks.entries) == [
        "hexahedron",
        "pyramid",
        "tetra",
        "wedge",
        "wedge",
    ]
    assert _region(mesh, "cell", "Part 4").tag == 4  # blank title
    assert _types(mesh, _region(mesh, "cell", "included truss").entries) == ["line"]
    # Subset 1 holds part 1 and, through subset 2, parts 2 and 3.
    assert len(_region(mesh, "cell", "solid parts").entries) == 12
    assert len(_region(mesh, "cell", "quadratic subset").entries) == 7
    parts = np.concatenate(mesh.cell_data["radioss:part"])
    props = np.concatenate(mesh.cell_data["radioss:property"])
    assert set(zip(parts.tolist(), props.tolist())) == {
        (1, 11),
        (2, 12),
        (3, 13),
        (4, 14),
        (5, 15),
    }


def test_groups(read):
    mesh = read(MESHES / "deck_0000.rad")
    assert len(_region(mesh, "point", "hex nodes").entries) == 8
    assert len(_region(mesh, "point", "quadratic nodes").entries) == 20 + 10 + 4 + 6
    # ids 1 2 3 4 then -3: the hexahedron, a wedge and the degenerate tetra.
    assert _types(mesh, _region(mesh, "cell", "some bricks").entries) == [
        "hexahedron",
        "tetra",
        "wedge",
    ]
    assert _types(mesh, _region(mesh, "cell", "part 2 solids").entries) == [
        "hexahedron20",
        "tetra",
        "tetra10",
        "wedge",
    ]
    union = _region(mesh, "cell", "union")
    assert len(union.entries) == 7
    assert len(_region(mesh, "cell", "subset solids").entries) == 9
    assert len(_region(mesh, "cell", "quad shell").entries) == 1
    assert len(_region(mesh, "cell", "in a box").entries) == 0
    assert len(_region(mesh, "point", "skin nodes").entries) == 8 + 3


def test_surface_segments(read):
    mesh = read(MESHES / "deck_0000.rad")
    skin = _region(mesh, "side", "skin")
    assert skin.dim == 2 and len(skin.entries) == 3
    node_counts = sorted(
        len(facet_nodes(mesh, int(c), int(f))[1]) for c, f in skin.entries
    )
    assert node_counts == [4, 4, 6]  # two quad faces and a tetra10 face


def test_old_format_columns(read):
    mesh = read(MESHES / "old_0000.rad")
    assert [c.type for c in mesh.cells] == ["hexahedron"]
    assert _is_valid(np.asarray(mesh.points)[mesh.cells[0].data[0]], "hexahedron")
    assert int(mesh.field_data["radioss:version"]) == 44
    assert len(_region(mesh, "point", "bottom").entries) == 4


def test_errors(read, tmp_path):
    engine = tmp_path / "run_0001.rad"
    engine.write_text("#RADIOSS ENGINE\n/RUN/run/1\n10.0\n/END\n")
    with pytest.raises(meshioplusplus.ReadError, match="engine"):
        read(engine)
    deck = tmp_path / "bad_0000.rad"
    deck.write_text(
        "#RADIOSS STARTER\n/BEGIN\nbad\n      2019         0\n\n\n"
        "/NODE\n         1                 0.0                 0.0                 0.0\n"
        "/TRUSS/1\n         1         1         2\n/END\n"
    )
    with pytest.raises(meshioplusplus.ReadError, match="undefined node 2"):
        read(deck)


def test_sniff_and_dispatch(tmp_path):
    path = tmp_path / "deck.inp2"
    shutil.copy(MESHES / "old_0000.rad", path)
    assert meshioplusplus.sniff_format(path) == "radioss"
    assert _core.sniff_format(str(path)) == "radioss"
    assert meshioplusplus.read(MESHES / "old_0000.rad").cells[0].type == "hexahedron"

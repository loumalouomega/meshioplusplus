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


@pytest.mark.parametrize("name", ["gmsh_hex20", "gmsh_tet10"])
def test_gmsh_written_decks_match_their_msh_twins(read, name):
    """gmsh writes Radioss's own node order (``getVertexRAD``); read back, every
    cell has the nodes, in meshio++'s order, of the same cell in gmsh's .msh."""
    rad = read(MESHES / f"{name}_0000.rad")
    msh = meshioplusplus.read(MESHES / f"{name}.msh")
    (block,) = rad.cells
    (twin,) = [c for c in msh.cells if c.type == block.type]
    corners = 8 if block.type == "hexahedron20" else 4

    def by_corners(mesh, rows):
        pts = np.round(np.asarray(mesh.points), 4)  # .rad keeps 6 digits
        return {tuple(sorted(map(tuple, pts[row[:corners]]))): pts[row] for row in rows}

    ours = by_corners(rad, np.asarray(block.data))
    theirs = by_corners(msh, np.asarray(twin.data))
    assert len(ours) == len(theirs) == len(block.data)
    for key, nodes in ours.items():
        np.testing.assert_allclose(nodes, theirs[key], atol=1e-4)
    for row in np.asarray(block.data):
        assert _is_valid(np.asarray(rad.points)[row], block.type)


def _features_deck():
    """Two bricks (part 1) side by side on x, two shells (part 2) on their top,
    in millimetres with metre work units, with every box, generator and surface
    form the reader resolves."""
    i10 = lambda *v: "".join(f"{x:10d}" for x in v)  # noqa: E731
    r20 = lambda *v: "".join(f"{x:20.6f}" for x in v)  # noqa: E731
    node = {}
    lines = ["#RADIOSS STARTER", "/BEGIN", "features", i10(2022, 0)]
    lines += [f"{'kg':>20}{'mm':>20}{'s':>20}", f"{'kg':>20}{'m':>20}{'s':>20}"]
    lines.append("/NODE")
    for k in range(2):
        for j in range(2):
            for i in range(3):
                nid = 1 + i + 3 * j + 6 * k
                node[(i, j, k)] = nid
                lines.append(f"{nid:10d}" + r20(1000.0 * i, 1000.0 * j, 1000.0 * k))
    lines.append("/BRICK/1")
    for i in range(2):
        corners = [(i, 0, 0), (i + 1, 0, 0), (i + 1, 1, 0), (i, 1, 0)]
        ids = [node[c] for c in corners] + [node[(a, b, 1)] for a, b, _ in corners]
        lines.append(i10(101 + i, *ids))
    lines.append("/SHELL/2")
    for i in range(2):
        lines.append(
            i10(201 + i, node[(i, 0, 1)], node[(i + 1, 0, 1)], node[(i + 1, 1, 1)])
            + i10(node[(i, 1, 1)])
        )
    lines += ["/PART/1", "solid", i10(1, 5), "/PART/2", "skin", i10(2, 6)]
    lines += ["/BOX/RECTA/1", "near half", i10(0, 0, 0)]
    lines += [r20(-1.0, -1.0, -1.0), r20(1001.0, 1001.0, 1001.0)]
    # A sphere round node 1 (its centre given as a node), 10 mm across.
    lines += ["/BOX/SPHER/2", "corner", f"{1:10d}{'':20}" + r20(10.0), r20(0, 0, 0)]
    lines += ["/BOX/CYLIN/3", "x axis", i10(0, 0) + f"{'':10}" + r20(10.0)]
    lines += [r20(-1.0, 0.0, 0.0), r20(2001.0, 0.0, 0.0)]
    lines += ["/BOX/BOX/4", "near half but the corner", i10(1, -2)]
    lines += ["/GRNOD/BOX/10", "in box 1", i10(1)]
    lines += ["/GRNOD/BOX/11", "in box 4", i10(4)]
    lines += ["/GRNOD/BOX/12", "on the axis", i10(3)]
    lines += ["/GRBRIC/BOX/13", "bricks inside", i10(1)]
    lines += ["/GRBRIC/BOX2/14", "bricks touching", i10(1)]
    lines += ["/GRNOD/GENE/15", "ranges", i10(1, 3, 7, 8)]
    lines += ["/GRNOD/GEN_INCR/16", "every fifth", i10(1, 12, 5)]
    lines += ["/GRSHEL/GENE/17", "second shell", i10(202, 202)]
    lines += ["/GRNOD/NODE/18/1", "unit suffix", i10(5)]
    lines += ["/SURF/PART/EXT/20", "solid outside", i10(1)]
    lines += ["/SURF/PART/21", "skin shells", i10(2)]
    lines += ["/SURF/PART/ALL/22", "solid all", i10(1)]
    lines += ["/SURF/GRBRIC/EXT/23", "first brick outside", i10(13)]
    lines += ["/SURF/GRBRIC/FREE/24", "first brick free", i10(13)]
    lines += ["/SURF/SURF/25", "second brick free", i10(20, -24)]
    lines += ["/SURF/MAT/26", "material 6", i10(6)]
    lines.append("/END")
    return "\n".join(lines) + "\n"


def test_units_boxes_generators_and_surfaces(read, tmp_path):
    deck = tmp_path / "features_0000.rad"
    deck.write_text(_features_deck())
    mesh = read(deck)
    # Millimetres in, metres of work: lengths scaled by 1e-3.
    assert float(mesh.field_data["radioss:length_scale"]) == pytest.approx(1e-3)
    assert np.asarray(mesh.points).max() == pytest.approx(2.0)
    points = np.asarray(mesh.points)

    def ids_of(kind, name):
        return np.asarray(_region(mesh, kind, name).entries)

    assert len(ids_of("point", "in box 1")) == 8
    near = ids_of("point", "in box 4")
    assert len(near) == 7 and not np.any(np.all(points[near] == 0.0, axis=1))
    np.testing.assert_allclose(points[ids_of("point", "on the axis")][:, 1:], 0.0)
    assert len(ids_of("point", "on the axis")) == 3
    assert len(ids_of("cell", "bricks inside")) == 1
    assert len(ids_of("cell", "bricks touching")) == 2
    assert len(ids_of("point", "ranges")) == 5
    assert len(ids_of("point", "every fifth")) == 3
    assert len(ids_of("cell", "second shell")) == 1
    suffix = _region(mesh, "point", "unit suffix")
    assert suffix.tag == 18 and len(suffix.entries) == 1

    def faces(name):
        return len(np.asarray(_region(mesh, "side", name).entries))

    assert faces("solid outside") == 10  # 12 faces, the shared one twice
    assert faces("skin shells") == 2
    assert faces("solid all") == 12
    assert faces("first brick outside") == 6  # external to the group
    assert faces("first brick free") == 5  # the face shared with brick 102 is not
    assert faces("second brick free") == 5
    assert faces("material 6") == 2


def _more_deck():
    """The two bricks and shells of the features deck, plus a /UNIT'd node block
    in metres, a fixed skew turned 45 degrees about z, a box aligned with it,
    box surfaces and the two analytical surfaces."""
    i10 = lambda *v: "".join(f"{x:10d}" for x in v)  # noqa: E731
    r20 = lambda *v: "".join(f"{x:20.6f}" for x in v)  # noqa: E731
    lines = _features_deck().splitlines()[:-1]  # without /END
    lines += ["/UNIT/7", "metres", f"{'kg':>20}{'m':>20}{'s':>20}"]
    lines += ["/NODE/7", f"{13:10d}" + r20(0.5, 0.5, 2.0)]
    lines += ["/SKEW/FIX/3", "turned", r20(0, 0, 0), r20(1, 1, 0), r20(-1, 1, 0)]
    lines += ["/BOX/RECTA/5", "along the skew", i10(0, 0, 3)]
    lines += [r20(100.0, -300.0, -100.0), r20(0.0, 2500.0, 1500.0)]
    lines += ["/GRNOD/BOX/30", "in the skewed box", i10(5)]
    lines += ["/SURF/BOX/31", "shells in box 1", i10(1)]
    lines += ["/SURF/BOX2/32", "shells touching box 1", i10(1)]
    lines += ["/SURF/BOX/EXT/33", "outside in box 1", i10(1)]
    lines += ["/SURF/BOX/ALL/34", "all in box 1", i10(1)]
    lines += ["/SURF/PLANE/35", "floor", r20(0, 0, 1000.0), r20(0, 0, 2000.0)]
    lines += ["/SURF/ELLIPS/36", "egg", i10(3, 4), r20(1000.0, 0, 0)]
    lines += [r20(100.0, 200.0, 300.0)]
    lines.append("/END")
    return "\n".join(lines) + "\n"


def test_units_skews_box_surfaces_and_engine_deck(read, tmp_path):
    deck = tmp_path / "more_0000.rad"
    deck.write_text(_more_deck())
    (tmp_path / "more_0001.rad").write_text(
        "#RADIOSS ENGINE\n/RUN/more/1\n                10.0\n/ANIM/DT\n0.0 0.5\n"
        "/ANIM/ELEM/SIGX\n/TFILE/4\n0.01\n/END\n"
    )
    mesh = read(deck)
    points = np.asarray(mesh.points)
    # /NODE/7 is in metres (/UNIT/7), the rest in millimetres
    np.testing.assert_allclose(points[-1], [0.5, 0.5, 2.0])
    assert points[:-1].max() == pytest.approx(2.0)
    # the box along the skew: the nodes with y >= x and x + y <= 2.5
    inside = np.asarray(_region(mesh, "point", "in the skewed box").entries)
    assert len(inside) == 6
    assert np.all(points[inside][:, 1] >= points[inside][:, 0])

    def faces(name):
        return len(np.asarray(_region(mesh, "side", name).entries))

    assert faces("shells in box 1") == 1
    assert faces("shells touching box 1") == 2
    assert faces("outside in box 1") == 6  # a shell and five faces of brick 101
    assert faces("all in box 1") == 8  # and brick 101's inner face, from both sides
    np.testing.assert_allclose(
        mesh.field_data["radioss:surf_plane:35"], [0, 0, 1, 0, 0, 2]
    )
    egg = mesh.field_data["radioss:surf_ellips:36"]
    s = 0.5**0.5
    np.testing.assert_allclose(
        egg, [4, 1, 0, 0, 0.1, 0.2, 0.3, s, s, 0, -s, s, 0, 0, 0, 1], atol=1e-15
    )
    # the engine deck beside the starter
    fd = mesh.field_data
    np.testing.assert_array_equal(fd["radioss:engine:RUN/more/1"], [10.0])
    np.testing.assert_array_equal(fd["radioss:engine:ANIM/DT"], [0.0, 0.5])
    assert fd["radioss:engine:ANIM/ELEM/SIGX"].size == 0
    np.testing.assert_array_equal(fd["radioss:engine:TFILE/4"], [0.01])
    engine = read(tmp_path / "more_0001.rad")
    assert len(engine.points) == 0 and not engine.cells
    np.testing.assert_array_equal(
        engine.field_data["radioss:engine:ANIM/DT"], [0.0, 0.5]
    )


def test_version_4_decks_title_their_keywords(read, tmp_path):
    """A 4.x deck (the version on the #RADIOSS STARTER line, no /BEGIN): 8 and
    16 column fields and titles in the keyword, as OpenRadioss's BAR2V41BD02."""
    i8 = lambda *v: "".join(f"{x:8d}" for x in v)  # noqa: E731
    r16 = lambda *v: "".join(f"{x:16.9f}" for x in v)  # noqa: E731
    lines = ["#RADIOSS STARTER      41QUAD4X                0", "/TITLE", "a 4.1 deck"]
    lines.append("/NODE")
    for k, (x, y) in enumerate([(0, 0), (1, 0), (1, 1), (0, 1), (2, 0), (2, 1)], 1):
        lines.append(i8(k) + r16(0.0, x, y))
    lines += ["/QUAD/1", i8(1, 1, 2, 3, 4), i8(2, 2, 5, 6, 3)]
    lines += ["/PART/1/COPPER", i8(1, 1)]
    lines += ["/GRNOD/NODE/2/FIXED_NODES", i8(1, 4)]
    lines += ["/SKEW/FIX/1/turned", r16(0, 1, 0), r16(0, 0, 1)]
    lines.append("/END")
    deck = tmp_path / "old41_0000.rad"
    deck.write_text("\n".join(lines) + "\n")
    mesh = read(deck)
    assert int(mesh.field_data["radioss:version"]) == 41
    assert [(c.type, len(c.data)) for c in mesh.cells] == [("quad", 2)]
    assert np.asarray(mesh.points)[4].tolist() == [0.0, 2.0, 0.0]
    assert sorted(_region(mesh, "point", "FIXED_NODES").entries.tolist()) == [0, 3]
    assert len(_region(mesh, "cell", "COPPER").entries) == 2


# --------------------------------------------------------------------------- #
# Writing                                                                     #
# --------------------------------------------------------------------------- #

WRITERS = {
    "core": lambda p, m, stubs=False: _core.radioss_write(str(p), m, stubs),
    "python": lambda p, m, stubs=False: py_radioss.write(p, m, stubs=stubs),
}


def _every_card_mesh():
    """One cell of every type the writer has a card for, a node group and a
    surface on a solid face and on a shell."""
    from meshioplusplus._regions import Region

    c = np.array(
        [
            [0, 0, 0],
            [1, 0, 0],
            [1, 1, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 0, 1],
            [1, 1, 1],
            [0, 1, 1],
        ],
        dtype=float,
    )
    edges = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4)]
    edges += [(0, 4), (1, 5), (2, 6), (3, 7)]
    h20 = np.vstack([c] + [(c[a] + c[b])[None] / 2 for a, b in edges])
    t = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
    t_edges = [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)]
    t10 = np.vstack([t] + [(t[a] + t[b])[None] / 2 for a, b in t_edges])
    w = c[[0, 1, 3, 4, 5, 7]]
    p = np.vstack([c[:4], [[0.5, 0.5, 1.0]]])
    q = c[:4] + [0, 0, -1]
    blocks = [
        ("hexahedron20", h20),
        ("tetra10", t10 + [2, 0, 0]),
        ("hexahedron", c + [4, 0, 0]),
        ("tetra", t + [6, 0, 0]),
        ("wedge", w + [8, 0, 0]),
        ("pyramid", p + [10, 0, 0]),
        ("quad", q),
        ("triangle", q[:3] + [2, 0, 0]),
        ("line", q[:2] + [4, 0, 0]),
    ]
    pts, cells = [], []
    for ctype, x in blocks:
        cells.append((ctype, np.arange(len(x))[None] + sum(len(y) for y in pts)))
        pts.append(x)
    mesh = meshioplusplus.Mesh(np.vstack(pts), cells)
    mesh.regions = [
        Region("fixed", "point", np.array([0, 1, 2]), tag=4),
        Region("skin", "side", np.array([[2, 0], [6, 0], [7, 2]]), tag=7),
        Region("solids", "cell", np.array([0, 1, 2, 3, 4, 5]), tag=3),
    ]
    return mesh


def _same_mesh(a, b):
    np.testing.assert_array_equal(a.points, b.points)
    assert [c.type for c in a.cells] == [c.type for c in b.cells]
    for x, y in zip(a.cells, b.cells):
        np.testing.assert_array_equal(x.data, y.data)


@pytest.mark.parametrize("engine", sorted(WRITERS))
def test_every_card_round_trips(engine, tmp_path):
    mesh = _every_card_mesh()
    out = tmp_path / "every_0000.rad"
    WRITERS[engine](out, mesh)
    for back in (meshioplusplus.radioss.read(out), py_radioss.read(out)):
        _same_mesh(mesh, back)
        by_name = {(r.kind, r.name): r for r in back.regions}
        np.testing.assert_array_equal(by_name[("point", "fixed")].entries, [0, 1, 2])
        assert by_name[("point", "fixed")].tag == 4
        # The shell's own face, and the hexahedron's first face.
        sides = sorted(
            map(tuple, np.asarray(by_name[("side", "skin")].entries).tolist())
        )
        assert sides == [(2, 0), (6, 0), (7, 0)]
        # A region of one element family is a part.
        assert by_name[("cell", "solids")].tag == 3


@pytest.mark.parametrize("path", DECKS, ids=[p.name for p in DECKS])
def test_decks_round_trip_and_the_writers_agree(path, tmp_path):
    mesh = py_radioss.read(path)
    written = {}
    for engine, writer in WRITERS.items():
        for stubs in (False, True):
            out = tmp_path / engine / str(stubs) / path.name
            out.parent.mkdir(parents=True)
            writer(out, mesh, stubs)
            written[(engine, stubs)] = out.read_bytes()
    assert written[("core", False)] == written[("python", False)]
    assert written[("core", True)] == written[("python", True)]
    back = py_radioss.read(tmp_path / "core" / "False" / path.name)
    _same_mesh(mesh, back)
    for name in ("radioss:part", "radioss:property", "radioss:material"):
        for x, y in zip(mesh.cell_data[name], back.cell_data[name]):
            np.testing.assert_array_equal(x, y)
    assert _regions(mesh) == _regions(back)


def test_titles_stay_one_title_line(tmp_path):
    from meshioplusplus._regions import Region

    mesh = _every_card_mesh()
    names = ["#not a comment", "/not a keyword", "two\nlines", "x" * 150]
    mesh.regions = [
        Region(n, "point", np.array([k]), tag=k + 1) for k, n in enumerate(names)
    ]
    for engine, writer in WRITERS.items():
        out = tmp_path / f"{engine}_0000.rad"
        writer(out, mesh)
        back = {
            r.tag: r.name for r in py_radioss.read(out).regions if r.kind == "point"
        }
        assert back == {
            1: "#not a comment",
            2: "/not a keyword",
            3: "two lines",
            4: "x" * 100,
        }


def _stub_props(text):
    """{property id: (type, first field)} and {part id: property id} of a deck."""
    lines = text.split("\n")
    props, parts = {}, {}
    for k, line in enumerate(lines):
        if line.startswith("/PROP/"):
            kind, pid = line.split("/")[2:4]
            props[int(pid)] = (kind, lines[k + 2][:10].strip())
        if line.startswith("/PART/"):
            parts[int(line[6:])] = int(lines[k + 2][:10])
    return props, parts


def test_stub_properties_fit_each_card(tmp_path):
    # The every-card mesh: its "solids" region (every solid) is one part, whose
    # bricks mix /BRICK and /BRIC20 (Isolid 0, which the starter turns to 16
    # for the /BRIC20); the shells and the truss have their own types.
    mesh = _every_card_mesh()
    # A /BRIC20 on its own gets Isolid 16.
    alone = meshioplusplus.Mesh(mesh.points, [mesh.cells[0]])
    for engine, writer in WRITERS.items():
        out = tmp_path / f"{engine}_0000.rad"
        writer(out, mesh, True)
        text = out.read_text()
        assert "/MAT/LAW1/1" in text
        props, parts = _stub_props(text)
        assert sorted(props[p] for p in set(parts.values())) == [
            ("SHELL", "0"),
            ("SHELL", "0"),
            ("SOLID", "0"),
            ("TRUSS", ""),
        ]
        writer(out, alone, True)
        props, _ = _stub_props(out.read_text())
        assert list(props.values()) == [("SOLID", "16")]
    # A part of solids and shells has no one stub property.
    mixed = _every_card_mesh()
    mixed.cell_data["radioss:part"] = [
        np.ones(len(c.data), dtype=int) for c in mixed.cells
    ]
    for writer in WRITERS.values():
        with pytest.raises(meshioplusplus.WriteError, match="no one stub property"):
            writer(tmp_path / "m_0000.rad", mixed, True)


def test_write_dispatch_and_refusals(tmp_path):
    mesh = _every_card_mesh()
    meshioplusplus.write(tmp_path / "d_0000.rad", mesh)
    assert "/BRIC20/" in (tmp_path / "d_0000.rad").read_text()
    empty = meshioplusplus.Mesh(np.zeros((0, 3)), [])
    for writer in WRITERS.values():
        with pytest.raises(meshioplusplus.WriteError, match="needs nodes"):
            writer(tmp_path / "e_0000.rad", empty)


def _openradioss():
    import os

    root = os.environ.get("MESHIOPLUSPLUS_OPENRADIOSS")
    if not root:
        pytest.skip(
            "set MESHIOPLUSPLUS_OPENRADIOSS to an OpenRadioss install to run its starter"
        )
    return pathlib.Path(root)


def _run_starter(root, deck):
    import os
    import subprocess

    env = dict(os.environ)
    env.update(
        OPENRADIOSS_PATH=str(root),
        RAD_CFG_PATH=str(root / "hm_cfg_files"),
        OMP_NUM_THREADS="1",
        LD_LIBRARY_PATH=os.pathsep.join(
            [
                str(root / "extlib" / "hm_reader" / "linux64"),
                env.get("LD_LIBRARY_PATH", ""),
            ]
        ),
    )
    exe = root / "exec" / "starter_linux64_gf"
    run = subprocess.run(
        [str(exe), "-i", deck.name, "-np", "1"],
        cwd=deck.parent,
        env=env,
        capture_output=True,
        text=True,
        timeout=600,
    )
    report = (deck.parent / (deck.stem + ".out")).read_text(errors="replace")
    return run.returncode, report


@pytest.mark.parametrize("path", DECKS, ids=[p.name for p in DECKS])
def test_openradioss_starter_accepts_written_decks(path, tmp_path):
    # Opt-in: MESHIOPLUSPLUS_OPENRADIOSS names an OpenRadioss install (its
    # exec/, hm_cfg_files/ and extlib/); v16.17.0 was checked with the
    # latest-20260728 build.
    root = _openradioss()
    deck = tmp_path / path.name
    meshioplusplus.radioss.write(deck, py_radioss.read(path), stubs=True)
    code, report = _run_starter(root, deck)
    assert code == 0, report[-3000:]
    assert "ERROR TERMINATION" not in report

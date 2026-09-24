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

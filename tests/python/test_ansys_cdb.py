"""Ansys MAPDL coded database (``.cdb``): both engines against the decks
mapdl-archive ships (and its frozen reading of them), degenerate shapes, missing
midside nodes, components as regions, and the writer's MAPDL layouts."""

import collections
import pathlib
import re

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _core
from meshioplusplus.ansysInp import _ansysInp as py_ansys

from .test_node_order import _SHAPES, _is_valid, _reference

MESHES = pathlib.Path(__file__).parent / "meshes" / "ansys"
FIXTURES = sorted(MESHES.glob("*.cdb"))

# VTK cell type ids mapdl-archive's grid uses, as meshio++ types.
_VTK = {
    1: "vertex",
    3: "line",
    5: "triangle",
    9: "quad",
    10: "tetra",
    12: "hexahedron",
    13: "wedge",
    14: "pyramid",
    21: "line3",
    22: "triangle6",
    23: "quad8",
    24: "tetra10",
    25: "hexahedron20",
    26: "wedge15",
    27: "pyramid13",
}


@pytest.fixture(params=["core", "python"])
def engine(request):
    """Both engines behind the same read/write signatures."""
    if request.param == "core":
        return meshioplusplus.ansysInp
    return py_ansys


def _regions(mesh):
    return sorted(
        (r.kind, r.name, tuple(np.asarray(r.entries).ravel())) for r in mesh.regions
    )


def _same(a, b):
    np.testing.assert_array_equal(a.points, b.points)
    assert [c.type for c in a.cells] == [c.type for c in b.cells]
    for x, y in zip(a.cells, b.cells):
        np.testing.assert_array_equal(x.data, y.data)
    assert sorted(a.cell_data) == sorted(b.cell_data)
    for name in a.cell_data:
        for x, y in zip(a.cell_data[name], b.cell_data[name]):
            np.testing.assert_array_equal(x, y)
    assert _regions(a) == _regions(b)


def _cell_keys(points, cell_type, rows):
    """A multiset of cells, each as its type and its sorted node coordinates."""
    out = collections.Counter()
    for row in rows:
        coords = tuple(sorted(tuple(np.round(p, 9)) for p in points[row]))
        out[(cell_type, coords)] += 1
    return out


@pytest.mark.parametrize("path", FIXTURES, ids=[p.name for p in FIXTURES])
def test_engines_agree_on_every_fixture(path):
    _same(
        meshioplusplus.ansysInp.read(path, lenient=True),
        py_ansys.read(path, lenient=True),
    )


@pytest.mark.parametrize("path", FIXTURES, ids=[p.name for p in FIXTURES])
def test_cells_and_components_match_mapdl_archive(engine, path):
    """Every cell mapdl-archive builds, as the same type over the same node
    coordinates; every component with the same size."""
    ref = np.load(MESHES / "mapdl_archive_reference.npz")
    key = path.stem
    mesh = engine.read(path, lenient=True)

    points, offsets = ref[f"{key}/points"], ref[f"{key}/offsets"]
    conn = ref[f"{key}/connectivity"]
    expected = collections.Counter()
    for c, vtk in enumerate(ref[f"{key}/celltypes"]):
        expected += _cell_keys(
            points, _VTK[int(vtk)], [conn[offsets[c] : offsets[c + 1]]]
        )
    got = collections.Counter()
    for block in mesh.cells:
        got += _cell_keys(mesh.points, block.type, block.data)
    missing, extra = expected - got, got - expected
    if key == "mixed_missing_midside":
        # mapdl-archive puts three missing triangle midsides at the origin;
        # meshio++ puts them at their edge midpoints.
        assert sum(missing.values()) == sum(extra.values()) == 3
        assert {t for t, _ in missing} == {t for t, _ in extra} == {"triangle6"}
    else:
        assert not missing and not extra

    for kind, region_kind in (("node", "point"), ("elem", "cell")):
        sizes = dict(
            zip(ref[f"{key}/{kind}_components"].tolist(), ref[f"{key}/{kind}_sizes"])
        )
        ours = {r.name: len(r.entries) for r in mesh.regions if r.kind == region_kind}
        assert ours == sizes


@pytest.mark.parametrize("name", ["all_solid_cells", "sector", "HexBeam"])
def test_solids_are_positively_oriented(engine, name):
    mesh = engine.read(MESHES / f"{name}.cdb")
    checked = 0
    for block in mesh.cells:
        corners = len(_SHAPES[block.type][0])
        for row in block.data:
            # Only the corners: real decks put midsides on curved edges.
            assert _is_valid(
                mesh.points[row[:corners]], block.type.rstrip("0123456789")
            )
            checked += 1
    assert checked > 0


def test_degenerate_solid186_shapes(engine):
    mesh = engine.read(MESHES / "all_solid_cells.cdb")
    assert [b.type for b in mesh.cells] == [
        "hexahedron20",
        "wedge15",
        "pyramid13",
        "tetra10",
    ]
    for data in mesh.cell_data["ansys:element"]:
        np.testing.assert_array_equal(data, [186])


def test_missing_midside_nodes_sit_on_their_edges(engine, capfd):
    mesh = engine.read(MESHES / "mixed_missing_midside.cdb")
    found = re.search(r"(\d+) missing midside node\(s\)", capfd.readouterr().err)
    n_file = len(mesh.points) - int(found.group(1))
    tets = [b for b in mesh.cells if b.type == "tetra10"][0]
    # SOLID92 rows written with nine nodes (the last midside left off) still
    # make quadratic tetrahedra.
    assert len(tets.data) == 236
    added = 0
    edges = [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)]
    for row in tets.data:
        for k, (a, b) in enumerate(edges):
            if row[4 + k] >= n_file:
                x = mesh.points[row]
                np.testing.assert_allclose(x[4 + k], (x[a] + x[b]) / 2)
                added += 1
    assert added > 0


def test_workbench_mixed_width_node_format(engine):
    mesh = engine.read(MESHES / "workbench_193.cdb")
    assert mesh.points.shape == (3, 3)
    assert not mesh.cells


def test_components_are_regions_and_sets(tmp_path):
    mesh = meshioplusplus.read(MESHES / "HexBeam.cdb")
    assert sorted(mesh.point_sets) == ["NCOMP2", "NODE_SELECTION"]
    assert sorted(mesh.cell_sets) == ["ECOMP1", "ECOMP2"]
    assert sum(len(s) for s in mesh.cell_sets["ECOMP1"]) == 22


def _mixed_mesh():
    """Every cell type the writer holds, one cell each (reference geometry),
    plus point and cell regions."""
    types = [
        "hexahedron",
        "hexahedron20",
        "wedge",
        "wedge15",
        "pyramid",
        "pyramid13",
        "tetra",
        "tetra10",
        "quad",
        "quad8",
        "triangle",
        "triangle6",
    ]
    points, cells = [], []
    for k, cell_type in enumerate(types):
        shape = _reference("triangle6")[:3] if cell_type == "triangle" else None
        x = (_reference(cell_type) if shape is None else shape) + [3.0 * k, 0.0, 0.0]
        cells.append((cell_type, np.arange(len(points), len(points) + len(x))[None]))
        points += x.tolist()
    base = len(points)
    points += [[0.0, 5.0, 0.0], [1.0, 5.0, 0.0], [0.5, 5.0, 0.0]]
    cells += [
        ("line", [[base, base + 1]]),
        ("line3", [[base, base + 1, base + 2]]),
        ("vertex", [[base]]),
    ]
    mesh = meshioplusplus.Mesh(np.array(points), cells)
    mesh.regions = [
        meshioplusplus.Region("SOLIDS", "cell", np.arange(0, 8)),
        meshioplusplus.Region("FIRST", "point", [0, 1, 2, 3, 7]),
        meshioplusplus.Region("FIRST", "cell", [0, 14]),
    ]
    return mesh


def test_engines_write_the_same_bytes(tmp_path):
    mesh = _mixed_mesh()
    _core.ansysinp_write(str(tmp_path / "cpp.cdb"), mesh)
    py_ansys.write(tmp_path / "py.cdb", mesh)
    cpp = (tmp_path / "cpp.cdb").read_bytes()
    assert cpp == (tmp_path / "py.cdb").read_bytes()
    assert b"Written by meshio++" in cpp


def test_round_trip_keeps_cells_and_regions(engine, tmp_path):
    mesh = _mixed_mesh()
    engine.write(tmp_path / "m.cdb", mesh)
    back = engine.read(tmp_path / "m.cdb")
    np.testing.assert_array_equal(back.points, mesh.points)
    assert [b.type for b in back.cells] == [b.type for b in mesh.cells]
    for x, y in zip(back.cells, mesh.cells):
        np.testing.assert_array_equal(x.data, y.data)
    assert _regions(back) == _regions(mesh)
    for block in back.cells:
        if block.type in _SHAPES and len(_SHAPES[block.type][0]) > 3:
            assert _is_valid(back.points[block.data[0]], block.type), block.type


def test_writer_uses_mapdl_layouts(engine, tmp_path):
    """Wedges, pyramids and tetrahedra as degenerate bricks under SOLID185/186
    (the native SOLID187/285 for tetrahedra), triangles as degenerate quads."""
    engine.write(tmp_path / "m.cdb", _mixed_mesh())
    lines = (tmp_path / "m.cdb").read_text().splitlines()
    ets = sorted(line for line in lines if line.startswith("ET,"))
    assert ets == sorted(
        f"ET,{k},{r}"
        for k, r in enumerate([185, 186, 285, 187, 181, 281, 188, 189, 21], start=1)
    )
    start = lines.index("(19i10)", lines.index("(3i9,6e21.13e3)") + 1) + 1
    rows = []
    k = start
    while lines[k].strip() != "-1":
        head = [int(lines[k][c : c + 10]) for c in range(0, len(lines[k]), 10)]
        nodes = head[11:]
        k += 1
        if head[8] > 8:
            nodes += [int(lines[k][c : c + 10]) for c in range(0, len(lines[k]), 10)]
            k += 1
        assert len(nodes) == head[8]
        rows.append(nodes)
    wedge = rows[2]
    assert len(wedge) == 8 and wedge[2] == wedge[3] and wedge[6] == wedge[7]
    tetra10 = rows[7]
    assert len(tetra10) == 10  # SOLID187, not a degenerate SOLID186
    tri = rows[10]
    assert len(tri) == 4 and tri[2] == tri[3]


def test_element_attributes_survive_a_round_trip(engine, tmp_path):
    mesh = engine.read(MESHES / "HexBeam.cdb")
    engine.write(tmp_path / "m.cdb", mesh)
    back = engine.read(tmp_path / "m.cdb")
    for name in (
        "ansys:element",
        "ansys:type",
        "ansys:mat",
        "ansys:real",
        "ansys:secnum",
    ):
        for x, y in zip(mesh.cell_data[name], back.cell_data[name]):
            np.testing.assert_array_equal(x, y)
    assert _regions(back) == _regions(mesh)


def test_element_type_that_does_not_fit_falls_back(engine, tmp_path, capfd):
    mesh = _mixed_mesh()
    mesh.cell_data["ansys:element"] = [
        np.array([187 if b.type == "hexahedron" else 0]) for b in mesh.cells
    ]
    engine.write(tmp_path / "m.cdb", mesh)
    assert "do not fit" in capfd.readouterr().err
    back = engine.read(tmp_path / "m.cdb")
    assert back.cell_data["ansys:element"][0].tolist() == [185]


def test_side_regions_are_dropped_with_a_warning(engine, tmp_path, capfd):
    mesh = _mixed_mesh()
    mesh.regions.append(meshioplusplus.Region("WALL", "side", [[0, 0]], dim=2))
    engine.write(tmp_path / "m.cdb", mesh)
    assert "side region" in capfd.readouterr().err
    text = (tmp_path / "m.cdb").read_text()
    assert "WALL" not in text and "regions-dropped" in text


_DECK = """\
/PREP7
{et}
NBLOCK,6,SOLID
(3i9,6e21.13e3)
        1        0        0 0.0000000000000E+000 0.0000000000000E+000 0.0000000000000E+000
        2        0        0 1.0000000000000E+000 0.0000000000000E+000 0.0000000000000E+000
        3        0        0 1.0000000000000E+000 1.0000000000000E+000 0.0000000000000E+000
        4        0        0 0.0000000000000E+000 1.0000000000000E+000 0.0000000000000E+000
N,R5.3,LOC,       -1,
EBLOCK,19,SOLID
(19i9)
{row}
       -1
{tail}FINISH
"""
_QUAD_ROW = "        1        1        1        1        0        0        0        0        4        0        1        1        2        3        4"  # noqa: E501


@pytest.mark.parametrize(
    "et, cell_type",
    [
        ("ET,1,181", "quad"),
        ("ET,1,SHELL181", "quad"),
        ("et, 1, plane182", "quad"),
        ("ET,1,200\nKEYOP,1,1,6", "quad"),  # MESH200 KEYOPT(1) = 6: 4-node area
        ("ET,1,200\nKEYOPT,1,1,2", "line3"),  # a set third node: quadratic
        ("ET,1,21", "vertex"),
    ],
)
def test_element_types_by_number_name_and_keyopt(engine, tmp_path, et, cell_type):
    path = tmp_path / "m.cdb"
    path.write_text(_DECK.format(et=et, row=_QUAD_ROW, tail=""))
    assert [b.type for b in engine.read(path).cells] == [cell_type]


def test_unknown_element_type_fails_or_is_skipped(engine, tmp_path, capfd):
    path = tmp_path / "m.cdb"
    path.write_text(_DECK.format(et="ET,1,999", row=_QUAD_ROW, tail=""))
    with pytest.raises(meshioplusplus.ReadError, match="element type 999"):
        engine.read(path)
    assert not engine.read(path, lenient=True).cells
    assert "skipped" in capfd.readouterr().err


@pytest.mark.parametrize(
    "et, row, tail, match",
    [
        ("", _QUAD_ROW, "", "no ET or ETBLOCK"),
        (
            "ET,1,181",
            _QUAD_ROW.replace("        4\n", "        9\n")[:-1] + "9",
            "",
            "node 9",
        ),
        ("ET,1,181", _QUAD_ROW, "CMBLOCK,C,NODE,1\n(8i10)\n        -3\n", "range end"),
        (
            "ET,1,181",
            _QUAD_ROW.replace("        2        3", "        x        3"),
            "",
            "bad integer",
        ),
    ],
)
def test_errors_name_the_culprit(engine, tmp_path, et, row, tail, match):
    path = tmp_path / "bad.cdb"
    path.write_text(_DECK.format(et=et, row=row, tail=tail))
    with pytest.raises(meshioplusplus.ReadError, match=match):
        engine.read(path)


def test_short_component_block_and_comments(engine, tmp_path):
    """A CMBLOCK whose header count overstates its entries ends at the next
    command; `!` comments on command lines are ignored."""
    tail = "CMBLOCK,TOP  ,NODE,       5  ! two ranges\n(8i10)\n         3        -4\n"
    path = tmp_path / "m.cdb"
    path.write_text(_DECK.format(et="ET,1,181 ! shell", row=_QUAD_ROW, tail=tail))
    mesh = engine.read(path)
    assert mesh.point_sets["TOP"].tolist() == [2, 3]


def test_non_solid_eblock_is_skipped_with_a_warning(engine, tmp_path, capfd):
    path = tmp_path / "m.cdb"
    text = _DECK.format(et="ET,1,181", row=_QUAD_ROW, tail="")
    text = text.replace("EBLOCK,19,SOLID", "EBLOCK,10,,1")
    path.write_text(text)
    mesh = engine.read(path)
    assert not mesh.cells
    assert "non-solid EBLOCK" in capfd.readouterr().err


def test_extension_and_registry(tmp_path):
    assert "ansysInp" in meshioplusplus._helpers._filetypes_from_path(
        pathlib.Path("x.cdb")
    )
    mesh = meshioplusplus.read(MESHES / "sector.cdb")
    meshioplusplus.write(tmp_path / "m.cdb", mesh)
    back = meshioplusplus.read(tmp_path / "m.cdb")
    assert _regions(back) == _regions(mesh)

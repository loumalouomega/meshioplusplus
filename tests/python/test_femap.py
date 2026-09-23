"""Femap neutral file (``.neu``): both engines, real Femap and EMSolution files
(``tools/gen_femap_fixtures.py``), output sets as steps, and round trips."""

import pathlib

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _core
from meshioplusplus.femap import _femap as py_femap

from .test_node_order import _is_valid

MESHES = pathlib.Path(__file__).parent / "meshes" / "femap"
FIXTURES = sorted(MESHES.glob("*.neu"))


@pytest.fixture(params=["core", "python"])
def engine(request):
    """Both engines behind the same read/write signatures."""
    if request.param == "core":
        return meshioplusplus.femap
    return py_femap


def _blocks(mesh):
    return [(b.type, np.asarray(b.data)) for b in mesh.cells]


def _regions(mesh):
    return sorted(
        (r.kind, r.name, r.dim, r.tag, tuple(np.asarray(r.entries).ravel()))
        for r in mesh.regions
    )


def _same(a, b):
    np.testing.assert_array_equal(a.points, b.points)
    assert [t for t, _ in _blocks(a)] == [t for t, _ in _blocks(b)]
    for (_, x), (_, y) in zip(_blocks(a), _blocks(b)):
        np.testing.assert_array_equal(x, y)
    assert _regions(a) == _regions(b)
    assert sorted(a.point_data) == sorted(b.point_data)
    for k in a.point_data:
        np.testing.assert_array_equal(a.point_data[k], b.point_data[k])
    assert sorted(a.cell_data) == sorted(b.cell_data)
    for k in a.cell_data:
        for x, y in zip(a.cell_data[k], b.cell_data[k]):
            np.testing.assert_array_equal(x, y)
    assert sorted(a.field_data) == sorted(b.field_data)


@pytest.mark.parametrize("path", FIXTURES, ids=[p.name for p in FIXTURES])
def test_engines_agree_on_every_fixture(path):
    steps = len(meshioplusplus.femap.time_values(path))
    assert steps == len(py_femap.time_values(path))
    for step in range(max(steps, 1)):
        _same(
            meshioplusplus.femap.read(path, time_step=step),
            py_femap.read(path, time_step=step),
        )


@pytest.mark.parametrize(
    "name, cell_type",
    [
        ("A342", "tetra10"),
        ("A352", "wedge15"),
        ("A362", "hexahedron20"),
        ("mesh_sample", None),
    ],
)
def test_cells_are_valid_meshio_cells(engine, name, cell_type):
    """Femap's degenerate-brick node slots, read from real files, give valid,
    positively oriented meshio++ cells (mid-edge nodes on their edges)."""
    mesh = engine.read(MESHES / f"{name}.neu")
    types = [t for t, _ in _blocks(mesh)]
    if cell_type:
        assert types == [cell_type]
    else:
        assert types == ["hexahedron", "pyramid", "tetra", "wedge"]
    for t, data in _blocks(mesh):
        for row in data:
            assert _is_valid(mesh.points[row], t), t


def test_properties_and_groups_become_regions(engine):
    mesh = engine.read(MESHES / "MC361.neu")
    regions = {(r.kind, r.name): r for r in mesh.regions}
    assert ("cell", "EALL") in regions  # the property, by its 402 title
    assert regions[("cell", "EALL")].tag == 1
    assert ("point", "XMIN NSET") in regions
    xmin = regions[("point", "XMIN NSET")]
    np.testing.assert_allclose(mesh.points[xmin.entries, 0], 0.0)
    assert len(xmin.entries) == 9
    assert [list(np.unique(p)) for p in mesh.cell_data["femap:property"]] == [[1]]
    assert [list(np.unique(t)) for t in mesh.cell_data["femap:type"]] == [[25]]


def test_output_sets_are_steps(engine):
    path = MESHES / "ems_results_1051.neu"
    times = engine.time_values(path)
    assert len(times) == 13 and times[1] == pytest.approx(0.00166667)
    last = engine.read(path, time_step=-1)
    assert float(last.field_data["meshio:time"]) == pytest.approx(times[-1])
    assert int(last.field_data["femap:set"]) == 13
    with pytest.raises(meshioplusplus.ReadError):
        engine.read(path, time_step=13)


def test_both_result_encodings_give_the_same_data(engine):
    """EMSolution wrote the same 13 output sets as 451 records and as 1051
    ranges; both read to the same arrays."""
    for step in (0, 5, 12):
        a = engine.read(MESHES / "ems_results_451.neu", time_step=step)
        b = engine.read(MESHES / "ems_results_1051.neu", time_step=step)
        assert sorted(a.point_data) == sorted(b.point_data) and a.point_data
        for k in a.point_data:
            np.testing.assert_array_equal(a.point_data[k], b.point_data[k])
        for k in a.cell_data:
            for x, y in zip(a.cell_data[k], b.cell_data[k]):
                np.testing.assert_array_equal(x, y)


@pytest.mark.parametrize("name", ["v2020_results", "v82_results", "mystran_results"])
def test_result_blocks_of_three_writers(engine, name):
    """Output blocks written by Femap 2020.1 (1051, seven-line output set
    tails), Femap 8.2 (451) and MYSTRAN (451, a line before the first block).
    The 8.2 and MYSTRAN values match femap_neutral_parser's reading of them."""
    mesh = engine.read(MESHES / f"{name}.neu", time_step=0)
    assert len(engine.time_values(MESHES / f"{name}.neu")) == 2
    assert [t for t, _ in _blocks(mesh)] == ["line"]
    if name == "mystran_results":
        t1 = mesh.point_data["T1  translation"]
    else:
        t1 = mesh.point_data["T1 Translation"]
    assert t1.shape == (12,)
    assert not np.isnan(t1).all()


def test_v2020_and_v82_hold_the_same_results(engine):
    """The same Nastran run exported by Femap 2020.1 (1051) and 8.2 (451)."""
    for step in (0, 1):
        a = engine.read(MESHES / "v2020_results.neu", time_step=step)
        b = engine.read(MESHES / "v82_results.neu", time_step=step)
        common = set(a.point_data) & set(b.point_data)
        assert "Total Translation" in common
        for k in common:
            np.testing.assert_allclose(a.point_data[k], b.point_data[k], rtol=1e-6)


def test_points_only_and_array_selection(engine):
    path = MESHES / "ems_results_451.neu"
    assert not engine.read(path, points_only=True).point_data
    full = engine.read(path)
    name = sorted(full.point_data)[0]
    picked = engine.read(path, arrays=[name])
    assert list(picked.point_data) == [name] and not [
        k for k in picked.cell_data if not k.startswith("femap:")
    ]


@pytest.mark.parametrize("path", FIXTURES, ids=[p.name for p in FIXTURES])
def test_engines_write_the_same_bytes(path, tmp_path):
    mesh = meshioplusplus.femap.read(path)
    _core.femap_write(str(tmp_path / "cpp.neu"), mesh)
    py_femap.write(tmp_path / "py.neu", mesh)
    assert (tmp_path / "cpp.neu").read_bytes() == (tmp_path / "py.neu").read_bytes()


@pytest.mark.parametrize(
    "name", ["A342", "A352", "A362", "A731", "B741", "MC361", "mesh_sample"]
)
def test_round_trip(engine, name, tmp_path):
    mesh = engine.read(MESHES / f"{name}.neu")
    engine.write(tmp_path / "out.neu", mesh)
    back = engine.read(tmp_path / "out.neu")
    np.testing.assert_array_equal(back.points, mesh.points)
    for (t, x), (u, y) in zip(_blocks(mesh), _blocks(back)):
        assert t == u
        np.testing.assert_array_equal(x, y)
    assert _regions(back) == _regions(mesh)
    for k in ("femap:property", "femap:type"):
        for x, y in zip(mesh.cell_data[k], back.cell_data[k]):
            np.testing.assert_array_equal(x, y)


def test_writer_makes_groups_and_drops_what_it_cannot_hold(engine, tmp_path, capfd):
    points = np.array(
        [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]], dtype=float
    )
    mesh = meshioplusplus.Mesh(
        points, [("tetra", [[0, 1, 2, 3]]), ("polygon", [[0, 1, 4, 2, 3]])]
    )
    mesh.point_data["T"] = np.zeros(5)
    mesh.regions = [
        meshioplusplus.Region("fixed", "point", [0, 1]),
        meshioplusplus.Region("fixed", "cell", [0], dim=3),
        meshioplusplus.Region("wall", "side", [[0, 1]], dim=2),
    ]
    engine.write(tmp_path / "m.neu", mesh)
    err = " ".join(capfd.readouterr().err.split())
    assert "polygon" in err and "side region" in err and "dropped" in err
    back = engine.read(tmp_path / "m.neu")
    assert [t for t, _ in _blocks(back)] == ["tetra"]
    regions = {(r.kind, r.name): r for r in back.regions}
    assert ("point", "fixed") in regions and ("cell", "fixed") in regions
    assert regions[("point", "fixed")].tag == regions[("cell", "fixed")].tag
    assert ("cell", "property_1") in regions


def _neu(*blocks):
    out = []
    for bid, lines in blocks:
        out.append(
            f"   -1\n{bid:6d}\n" + "".join(line + "\n" for line in lines) + "   -1\n"
        )
    return "".join(out)


_NODES = (
    403,
    [
        "1,0,0,1,46,0,0,0,0,0,0,0.,0.,0.,0,",
        "2,0,0,1,46,0,0,0,0,0,0,1.,0.,0.,0,",
    ],
)


def _bar(n2=2, topology=0):
    return (
        404,
        [
            f"7,124,1,1,{topology},1,0,0,0,0,0,0,0,",
            f"1,{n2},0,0,0,0,0,0,0,0,",
            "0,0,0,0,0,0,0,0,0,0,",
            "0.,0.,0.,",
            "0.,0.,0.,",
            "0.,0.,0.,",
            "0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,",
        ],
    )


@pytest.mark.parametrize(
    "blocks, match",
    [
        ([(100, ["<NULL>", "8.2,"])], "no nodes"),
        ([_NODES, _bar(n2=9)], "undefined node 9"),
        ([_NODES, _bar(n2=0)], "no node in slot 1"),
        ([_NODES, _NODES], "defined twice"),
        ([(403, ["1,0,0,1,46,0,0,0,0,0,0,0.,x,0.,0,"])], "bad coordinate"),
        ([(403, ["1,0,0,1,46,0.,0.,0.,"])], "fields"),
        (
            [_NODES, (404, ["7,124,1,1,0,1,0,0,", "1,2,0,0,0,0,0,0,0,0,"])],
            "ends inside",
        ),
    ],
)
def test_errors_name_the_culprit(engine, tmp_path, blocks, match):
    path = tmp_path / "bad.neu"
    path.write_text(_neu(*blocks))
    with pytest.raises(meshioplusplus.ReadError, match=match):
        engine.read(path)


def test_skipped_topologies_and_element_node_lists(engine, tmp_path, capfd):
    rigid = (
        404,
        [
            "8,124,1,29,13,1,0,0,0,0,0,0,0,",
            "1,0,0,0,0,0,0,0,0,0,",
            "0,0,0,0,0,0,0,0,0,0,",
            "0.,0.,0.,",
            "0.,0.,0.,",
            "0.,0.,0.,",
            "0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,",
            "2,0,1.,1,1,1,0,0,0,",
            "-1,",
        ],
    )
    path = tmp_path / "m.neu"
    path.write_text(_neu((100, ["<NULL>", "9.3,"]), _NODES, rigid, _bar()))
    mesh = engine.read(path)
    assert "rigid" in " ".join(capfd.readouterr().err.split())
    assert [(t, len(d)) for t, d in _blocks(mesh)] == [("line", 1)]


def test_extension_and_sniffing(tmp_path):
    assert meshioplusplus._helpers._filetypes_from_path(pathlib.Path("x.neu")) == [
        "femap"
    ]
    for p in (MESHES / "A361.neu", MESHES / "mystran_results.neu"):
        if p.exists():
            assert meshioplusplus.sniff_format(p) == "femap"
    path = tmp_path / "no_extension"
    path.write_bytes((MESHES / "mystran_results.neu").read_bytes())
    assert meshioplusplus.sniff_format(path) == "femap"

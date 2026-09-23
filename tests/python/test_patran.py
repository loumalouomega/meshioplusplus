"""MSC Patran 2 neutral file (``.pat``/``.out``): both engines, the fixtures
written by ``tools/gen_patran_fixtures.py`` in Patran's own node numbering, and
round trips."""

import pathlib

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _core
from meshioplusplus.patran import _patran as py_patran

from .test_node_order import _is_valid

MESHES = pathlib.Path(__file__).parent / "meshes" / "patran"
FIXTURES = sorted(MESHES.glob("*.pat"))


@pytest.fixture(params=["core", "python"])
def engine(request):
    """Both engines behind the same read/write signatures."""
    if request.param == "core":
        return meshioplusplus.patran
    return py_patran


def _blocks(mesh):
    return [(b.type, np.asarray(b.data)) for b in mesh.cells]


def _regions(mesh):
    return sorted(
        (r.kind, r.name, r.dim, r.tag, tuple(np.asarray(r.entries).ravel()))
        for r in mesh.regions
    )


def _same(a, b, exact=True):
    if exact:
        np.testing.assert_array_equal(a.points, b.points)
    else:
        np.testing.assert_allclose(a.points, b.points, rtol=1e-9, atol=1e-12)
    assert [t for t, _ in _blocks(a)] == [t for t, _ in _blocks(b)]
    for (_, x), (_, y) in zip(_blocks(a), _blocks(b)):
        np.testing.assert_array_equal(x, y)
    assert _regions(a) == _regions(b)


@pytest.mark.parametrize("path", FIXTURES, ids=[p.name for p in FIXTURES])
def test_engines_agree_on_every_fixture(path):
    a = meshioplusplus.patran.read(path)
    b = py_patran.read(path)
    _same(a, b)
    for x, y in zip(a.cell_data["patran:property"], b.cell_data["patran:property"]):
        np.testing.assert_array_equal(x, y)


def test_quadratic_cells_are_valid_in_meshio_order(engine):
    """The fixture lists mid-edge nodes in Patran's order (vertical edges before
    the top ring for hex20/wedge15); after reading, every cell is a valid,
    positively oriented meshio++ element."""
    mesh = engine.read(MESHES / "quadratic.pat")
    types = [t for t, _ in _blocks(mesh)]
    assert types == [
        "hexahedron20",
        "wedge15",
        "tetra10",
        "pyramid13",
        "quad8",
        "triangle6",
        "line3",
    ]
    for cell_type, data in _blocks(mesh):
        for row in data:
            assert _is_valid(mesh.points[row], cell_type), cell_type


def test_mixed_linear_components_and_property_fallback(engine):
    mesh = engine.read(MESHES / "mixed_linear.pat")
    assert [(t, len(d)) for t, d in _blocks(mesh)] == [
        ("hexahedron", 1),
        ("wedge", 1),
        ("tetra", 1),
        ("pyramid", 1),
        ("quad", 2),
        ("line", 1),
    ]
    for cell_type, data in _blocks(mesh):
        if cell_type != "line":
            for row in data:
                assert _is_valid(mesh.points[row], cell_type), cell_type
    regions = {(r.kind, r.name): r for r in mesh.regions}
    assert set(regions) == {
        ("cell", "SOLIDS"),
        ("point", "FIXED"),
        ("cell", "MIXED"),
        ("point", "MIXED"),
        ("cell", "property_2"),
        ("cell", "property_3"),
    }
    np.testing.assert_array_equal(regions[("cell", "SOLIDS")].entries, [0, 1, 2, 3])
    assert regions[("cell", "SOLIDS")].tag == 1 and regions[("cell", "SOLIDS")].dim == 3
    np.testing.assert_array_equal(regions[("point", "FIXED")].entries, [0, 1, 2, 3])
    np.testing.assert_array_equal(regions[("cell", "MIXED")].entries, [4])
    np.testing.assert_array_equal(regions[("point", "MIXED")].entries, [0, 1])
    # the second quad and the bar are in no component
    np.testing.assert_array_equal(regions[("cell", "property_2")].entries, [5])
    np.testing.assert_array_equal(regions[("cell", "property_3")].entries, [6])
    assert regions[("cell", "property_3")].tag == 3
    assert [list(v) for v in mesh.cell_data["patran:property"]] == [
        [1],
        [1],
        [1],
        [1],
        [2, 2],
        [3],
    ]


@pytest.mark.parametrize("path", FIXTURES, ids=[p.name for p in FIXTURES])
def test_engines_write_the_same_bytes(path, tmp_path):
    mesh = meshioplusplus.read(path)
    _core.patran_write(str(tmp_path / "cpp.pat"), mesh)
    py_patran.write(tmp_path / "py.pat", mesh)
    assert (tmp_path / "cpp.pat").read_bytes() == (tmp_path / "py.pat").read_bytes()


@pytest.mark.parametrize("path", FIXTURES, ids=[p.name for p in FIXTURES])
def test_round_trip_keeps_cells_regions_and_properties(engine, path, tmp_path):
    mesh = engine.read(path)
    engine.write(tmp_path / "out.pat", mesh)
    back = engine.read(tmp_path / "out.pat")
    _same(mesh, back, exact=False)
    for x, y in zip(
        mesh.cell_data["patran:property"], back.cell_data["patran:property"]
    ):
        np.testing.assert_array_equal(x, y)


def test_writer_merges_point_and_cell_regions_and_drops_the_rest(
    engine, tmp_path, capfd
):
    points = np.array(
        [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1.5]], dtype=float
    )
    mesh = meshioplusplus.Mesh(
        points, [("tetra", [[0, 1, 2, 3]]), ("vertex", [[4]]), ("line", [[0, 4]])]
    )
    mesh.point_data["T"] = np.zeros(len(points))
    mesh.regions = [
        meshioplusplus.Region("a_rather_long_component", "cell", [0, 2], dim=3),
        meshioplusplus.Region("a_rather_long_component", "point", [0, 4]),
        meshioplusplus.Region("WALL", "side", [[0, 1]], dim=2),
    ]
    engine.write(tmp_path / "m.pat", mesh)
    err = " ".join(capfd.readouterr().err.split())
    assert "vertex" in err and "side region" in err and "no data arrays" in err
    text = (tmp_path / "m.pat").read_text()
    assert "cells-dropped" not in text  # the title holds the tag only
    back = engine.read(tmp_path / "m.pat")
    assert [t for t, _ in _blocks(back)] == ["tetra", "line"]
    regions = {(r.kind, r.name): r for r in back.regions}
    assert set(regions) == {("cell", "a_rather_l"), ("point", "a_rather_l")} or set(
        regions
    ) == {("cell", "a_rather_lon"), ("point", "a_rather_lon")}
    cell = next(r for r in back.regions if r.kind == "cell")
    point = next(r for r in back.regions if r.kind == "point")
    assert cell.tag == point.tag
    # the line was cell 2 of the written mesh; after dropping the vertex it is 1
    np.testing.assert_array_equal(cell.entries, [0, 1])
    np.testing.assert_array_equal(point.entries, [0, 4])


def test_title_carries_the_provenance_tag(engine, tmp_path):
    mesh = meshioplusplus.read(MESHES / "mixed_linear.pat")
    engine.write(tmp_path / "m.pat", mesh)
    lines = (tmp_path / "m.pat").read_text().splitlines()
    assert lines[0].startswith("25")
    assert lines[1].startswith("Written by meshio++ v")


def _pat(*cards):
    return "".join(c + "\n" for c in cards)


_NODE = [
    " 1       1       0       2       0       0       0       0       0",
    " 0.000000000E+00 0.000000000E+00 0.000000000E+00",
    "1G       6       0       0  000000",
]


@pytest.mark.parametrize(
    "body, match",
    [
        (
            _pat(
                *_NODE,
                " 2       7       2       2       0       0       0       0       0",
                "       2       0       1       0",
                "       1       9",
            ),
            "undefined node 9",
        ),
        (_pat(*_NODE, *_NODE), "defined twice"),
        (
            _pat(
                " 1       1       0       2       0       0       0       0       0",
                " 0.0000000x0E+00",
                "1G",
            ),
            "invalid real",
        ),
        (
            _pat(" 1       1       0       5       0       0       0       0       0"),
            "file ends first",
        ),
        (
            _pat(
                *_NODE,
                " 2       7       2       1       0       0       0       0       0",
                "       2       0       1       0",
            ),
            "has 1 data cards",
        ),
    ],
)
def test_errors_name_the_culprit(engine, tmp_path, body, match):
    path = tmp_path / "bad.pat"
    path.write_text(body)
    with pytest.raises(meshioplusplus.ReadError, match=match):
        engine.read(path)


def test_unknown_shapes_and_component_types_are_skipped(engine, tmp_path, capfd):
    body = _pat(
        *_NODE,
        " 1       2       0       2       0       0       0       0       0",
        " 1.000000000E+00 0.000000000E+00 0.000000000E+00",
        "1G       6       0       0  000000",
        # a 4-node bar: no meshio++ equivalent
        " 2       5       2       2       0       0       0       0       0",
        "       4       0       1       0",
        "       1       2       1       2",
        " 2       6       2       2       0       0       0       0       0",
        "       2       0       1       0",
        "       1       2",
        # a component listing a coordinate frame (19) and a missing node
        "21       3       6       2       0       0       0       0       0",
        "C1",
        "      19       1       6       6       5      77",
        "99       0       0       1       0       0       0       0       0",
    )
    path = tmp_path / "m.pat"
    path.write_text(body)
    mesh = engine.read(path)
    err = " ".join(capfd.readouterr().err.split())
    assert "shape 2 with 4 nodes" in err
    assert "type 19" in err and "undefined" in err
    assert [t for t, _ in _blocks(mesh)] == ["line"]
    regions = {(r.kind, r.name): r for r in mesh.regions}
    np.testing.assert_array_equal(regions[("cell", "C1")].entries, [0])


def test_crlf_and_missing_end_packet(engine, tmp_path, capfd):
    body = _pat(*_NODE).replace("\n", "\r\n")
    path = tmp_path / "m.pat"
    path.write_bytes(body.encode())
    mesh = engine.read(path)
    assert "no end packet" in " ".join(capfd.readouterr().err.split())
    assert mesh.points.shape == (1, 3)


def test_extension_and_sniffing(tmp_path):
    for ext in (".pat", ".out"):
        assert meshioplusplus._helpers._filetypes_from_path(
            pathlib.Path("x" + ext)
        ) == ["patran"]
    path = tmp_path / "no_extension"
    path.write_bytes((MESHES / "mixed_linear.pat").read_bytes())
    assert meshioplusplus.sniff_format(path) == "patran"
    (tmp_path / "not_patran").write_text(
        "25 this is not a header card at all, just text\n"
    )
    assert meshioplusplus.sniff_format(tmp_path / "not_patran") != "patran"

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


REAL = sorted((MESHES / "real").glob("*.out")) + sorted((MESHES / "real").glob("*.pat"))


@pytest.mark.parametrize("path", REAL, ids=[p.name for p in REAL])
def test_real_exports(engine, path):
    """Files written by P3/PATRAN 3.0 and PATRAN 2.5 (WARP3D's and Tahoe's
    examples, see meshes/patran/real/README.md): both engines agree and every
    cell is positively oriented."""
    mesh = engine.read(path)
    _same(mesh, py_patran.read(path))
    assert mesh.cells
    for cell_type, data in _blocks(mesh):
        for row in data:
            assert _is_valid(mesh.points[row], cell_type), (path.name, cell_type)


def test_real_export_components_are_regions(engine):
    # Tahoe's square: 4 x 4 nodes on [-1, 1]^2 (P3/PATRAN 3.0), with named
    # components (packet 21) for the edges, two corners and the whole square.
    mesh = engine.read(MESHES / "real" / "tahoe_square.pat")
    regions = {(r.kind, r.name): np.asarray(r.entries) for r in mesh.regions}
    assert sorted(regions) == [
        ("cell", "SQUARE"),
        ("point", "BOTTOM"),
        ("point", "CORNER1"),
        ("point", "CORNER2"),
        ("point", "SQUARE"),
        ("point", "TOP"),
    ]
    assert len(regions[("cell", "SQUARE")]) == 9
    y = mesh.points[:, 1]
    np.testing.assert_array_equal(y[regions[("point", "BOTTOM")]], [y.min()] * 4)
    np.testing.assert_array_equal(y[regions[("point", "TOP")]], [y.max()] * 4)


def test_quad9_and_triangle7(engine, tmp_path):
    # Patran's QUAD9 and TRI7 list corners, mid-edges, then the centre: the
    # meshio++ order. Checked on a P3/PATRAN 3.0 export outside the repository
    # (Tahoe's cyl.out: every QUAD9 mid-edge node within 3 % of its edge's
    # midpoint on the curved surface).
    xy = [[0, 0], [2, 0], [2, 2], [0, 2], [1, 0], [2, 1], [1, 2], [0, 1], [1, 1]]
    xy += [[4, 0], [3, 1], [3, 0], [3, 1.5], [8 / 3, 2 / 3]]  # 12 is unused
    points = np.array([p + [0.0] for p in xy], dtype=float)
    mesh = meshioplusplus.Mesh(
        points,
        [
            ("quad9", np.array([[0, 1, 2, 3, 4, 5, 6, 7, 8]])),
            ("triangle7", np.array([[1, 9, 2, 11, 10, 5, 13]])),
        ],
    )
    assert _is_valid(points[mesh.cells[0].data[0]], "quad9")
    path = tmp_path / "mesh.pat"
    engine.write(path, mesh)
    # Packet 02 headers (I2,8I8): the shape code IV is the third field.
    headers = [line.split() for line in path.read_text().splitlines()]
    headers = [h for h in headers if len(h) == 9 and all(v.isdigit() for v in h)]
    assert [int(h[2]) for h in headers if h[0] == "2"] == [4, 3]
    back = engine.read(path)
    assert [b.type for b in back.cells] == ["quad9", "triangle7"]
    for x, y in zip(back.cells, mesh.cells):
        np.testing.assert_array_equal(x.data, y.data)


# --- loads, boundary conditions and result files ---------------------------------------------


def test_loads_and_boundary_conditions(engine):
    """mixed_linear.pat carries a uniform pressure on face 6 of the hex
    (packet 06) and a node temperature (packet 10)."""
    mesh = engine.read(MESHES / "mixed_linear.pat")
    np.testing.assert_array_equal(mesh.point_data["patran:temperature:1"][0], 300.0)
    assert np.isnan(mesh.point_data["patran:temperature:1"][1:]).all()
    flags = mesh.field_data["patran:distributed_load"]
    values = mesh.field_data["patran:distributed_load_values"]
    # cell 0 (the hex), set 1, surface load at the centroid, component 1, face 6
    np.testing.assert_array_equal(flags, [[0, 1, 1, 1, 0, 1] + [0] * 13 + [6]])
    assert values[0, 0] == 1.0 and np.isnan(values[0, 1:]).all()


def _synthetic_loads(mesh):
    n = len(mesh.points)
    ncells = sum(len(c.data) for c in mesh.cells)
    force = np.full((n, 6), np.nan)
    force[0, [0, 2]] = [1.5, -2.0]
    force[3] = np.arange(6)
    frame = np.zeros(n, dtype=np.int64)
    frame[3] = 5
    fixed = np.full((n, 6), np.nan)
    fixed[1, :3] = 0.0
    temps = np.full(ncells, np.nan)
    temps[1] = 25.0
    offsets = np.cumsum([0] + [len(c.data) for c in mesh.cells])
    mesh.point_data["patran:force:3"] = force
    mesh.point_data["patran:force_frame:3"] = frame
    mesh.point_data["patran:displacement:1"] = fixed
    mesh.cell_data["patran:element_temperature:7"] = [
        temps[offsets[b] : offsets[b + 1]] for b in range(len(mesh.cells))
    ]
    flags = np.zeros((2, 20), dtype=np.int64)
    flags[0] = [0, 4, 1, 1, 0, 1, 0, 0, 0, 0, 0] + [1, 1] + [0] * 6 + [2]
    flags[1] = [1, 4, 1, 0, 1, 1, 1, 0, 0, 0, 0] + [1, 1, 1] + [0] * 5 + [1]
    values = np.full((2, 54), np.nan)
    values[0, 0] = 10.0
    values[1, :6] = np.arange(6) + 0.5  # 2 components at 3 nodes
    mesh.field_data["patran:distributed_load"] = flags
    mesh.field_data["patran:distributed_load_values"] = values
    return mesh


@pytest.mark.parametrize("writer", ["core", "python"])
def test_loads_round_trip(writer, engine, tmp_path):
    mesh = _synthetic_loads(meshioplusplus.read(MESHES / "mixed_linear.pat"))
    w = meshioplusplus.patran if writer == "core" else py_patran
    w.write(tmp_path / "m.pat", mesh)
    back = engine.read(tmp_path / "m.pat")
    for key in ("patran:force:3", "patran:force_frame:3", "patran:displacement:1"):
        np.testing.assert_array_equal(back.point_data[key], mesh.point_data[key])
    for a, b in zip(
        back.cell_data["patran:element_temperature:7"],
        mesh.cell_data["patran:element_temperature:7"],
    ):
        np.testing.assert_array_equal(a, b)
    for key in ("patran:distributed_load", "patran:distributed_load_values"):
        np.testing.assert_array_equal(back.field_data[key], mesh.field_data[key])


def test_loads_write_the_same_bytes(tmp_path):
    mesh = _synthetic_loads(meshioplusplus.read(MESHES / "mixed_linear.pat"))
    meshioplusplus.patran.write(tmp_path / "a.pat", mesh)
    py_patran.write(tmp_path / "b.pat", mesh)
    assert (tmp_path / "a.pat").read_bytes() == (tmp_path / "b.pat").read_bytes()


WARP3D = MESHES / "real" / "warp3d_ssy"
WARP3D_RESULTS = {
    k: f"{WARP3D}.{k}" for k in ("wnfr00001", "wefe00001", "wnbd00001", "webs00001")
}


def test_warp3d_result_files(engine):
    """WARP3D's Patran 2.5 result files, text (wnfr nodal, wefe element) and
    Fortran binary (wnbd nodal, webs element), read onto its neutral file."""
    mesh = engine.read(f"{WARP3D}.out", WARP3D_RESULTS)
    assert [(c.type, len(c.data)) for c in mesh.cells] == [("hexahedron", 40)]
    assert mesh.point_data["wnfr00001"].shape == (164, 3)
    assert mesh.point_data["wnbd00001"].shape == (164, 3)
    assert mesh.cell_data["wefe00001"][0].shape == (40, 22)
    assert mesh.cell_data["webs00001"][0].shape == (40, 26)
    # the first records of the files, as printed
    np.testing.assert_array_equal(
        mesh.point_data["wnfr00001"][0], [0.390875e-04, -0.113240e-02, -0.139326e-03]
    )
    np.testing.assert_allclose(
        mesh.cell_data["wefe00001"][0][0, :3],
        [-0.583119e-03, 0.386872e-02, -0.113706e-03],
    )
    assert np.isfinite(mesh.point_data["wnbd00001"]).all()
    assert np.isfinite(mesh.cell_data["webs00001"][0]).all()
    # packet 08: WARP3D's constraints, two or three components per node
    fixed = mesh.point_data["patran:displacement:1"]
    assert (fixed[~np.isnan(fixed)] == 0).all() and (~np.isnan(fixed)).any()


def test_warp3d_results_engines_agree():
    a = meshioplusplus.patran.read(f"{WARP3D}.out", WARP3D_RESULTS)
    b = py_patran.read(f"{WARP3D}.out", WARP3D_RESULTS)
    for key in a.point_data:
        np.testing.assert_array_equal(a.point_data[key], b.point_data[key])
    for key in a.cell_data:
        for x, y in zip(a.cell_data[key], b.cell_data[key]):
            np.testing.assert_array_equal(x, y)


def test_result_file_errors(engine, tmp_path, capfd):
    bad = tmp_path / "bad.nod"
    bad.write_text("title\n     1     1  0.0  0     3\nsub\nsub\n       1 1.0\n")
    with pytest.raises(meshioplusplus.ReadError, match="ends inside"):
        engine.read(MESHES / "mixed_linear.pat", {"r": str(bad)})
    stray = tmp_path / "stray.nod"
    stray.write_text("title\n     1     1  0.0  0     1\nsub\nsub\n  999999 1.0\n")
    mesh = engine.read(MESHES / "mixed_linear.pat", {"r": str(stray)})
    assert np.isnan(mesh.point_data["r"]).all()
    assert "undefined node" in capfd.readouterr().err

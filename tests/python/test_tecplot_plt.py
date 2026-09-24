"""Binary Tecplot (.plt), ordered and face-based (FEPOLYGON/FEPOLYHEDRON)
zones: both engines, checked against the ASCII twin of every TecIO-written
fixture."""

import pathlib
import struct

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import ReadError
from meshioplusplus.tecplot import _tecplot

HERE = pathlib.Path(__file__).resolve().parent
PLT = HERE / "meshes" / "tecplot" / "plt"
NAMES = ["fe_mixed", "fe_2d", "ordered", "transient", "poly_2d", "poly_3d"]

try:
    from meshioplusplus import _core

    HAS_CORE = hasattr(_core, "tecplot_read")
except ImportError:  # pragma: no cover
    HAS_CORE = False

ENGINES = ["python"] + (["core"] if HAS_CORE else [])


def _read(engine, path, step=0):
    if engine == "core":
        return _core.tecplot_read(str(path), step)
    return _tecplot.read(str(path), step)


def _same(a, b):
    np.testing.assert_array_equal(a.points, b.points)
    assert [c.type for c in a.cells] == [c.type for c in b.cells]
    for ca, cb in zip(a.cells, b.cells):
        assert _cell_lists(ca) == _cell_lists(cb)
    assert set(a.point_data) == set(b.point_data)
    for k in a.point_data:
        np.testing.assert_array_equal(a.point_data[k], b.point_data[k])
    assert set(a.cell_data) == set(b.cell_data)
    for k in a.cell_data:
        for x, y in zip(a.cell_data[k], b.cell_data[k]):
            np.testing.assert_array_equal(x, y)
    key = lambda r: (r.name, r.kind, r.tag, tuple(r.entries))  # noqa: E731
    assert sorted(map(key, a.regions)) == sorted(map(key, b.regions))


def _cell_lists(block):
    """A block's cells as nested lists, ragged or not."""
    if block.type.startswith("polyhedron"):
        return [[list(map(int, f)) for f in cell] for cell in block.data]
    return [list(map(int, row)) for row in block.data]


def _volume(points, faces):
    """Signed volume of a polyhedron (positive when its faces are outward)."""
    v = 0.0
    for face in faces:
        q = points[np.asarray(face)]
        c = q.mean(axis=0)
        for i in range(len(q)):
            v += np.dot(c, np.cross(q[i], q[(i + 1) % len(q)])) / 6.0
    return v


def _steps(name):
    return max(1, len(_tecplot.time_values(PLT / f"{name}.dat")))


@pytest.mark.parametrize("engine", ENGINES)
@pytest.mark.parametrize("name", NAMES)
def test_plt_reads_as_its_ascii_twin(engine, name):
    for step in range(_steps(name)):
        _same(
            _read(engine, PLT / f"{name}.plt", step),
            _read(engine, PLT / f"{name}.dat", step),
        )


@pytest.mark.skipif(not HAS_CORE, reason="needs the C++ core")
@pytest.mark.parametrize("ext", [".plt", ".dat"])
@pytest.mark.parametrize("name", NAMES)
def test_engines_agree(name, ext):
    for step in range(_steps(name)):
        _same(
            _read("core", PLT / f"{name}{ext}", step),
            _read("python", PLT / f"{name}{ext}", step),
        )


def test_dispatch_sniff_and_metadata():
    path = PLT / "transient.plt"
    assert meshioplusplus.sniff_format(str(path)) == "tecplot"
    mesh = meshioplusplus.read(path, time_step=-1)
    assert [c.type for c in mesh.cells] == ["quad"]
    meta = meshioplusplus.read_metadata(str(path))
    assert meta["format"] == "tecplot"
    assert meta["time_values"] == [0.0, 0.5, 1.25]
    assert meta["num_points"] == 6
    # The ASCII twin under a .dat name goes to Tecplot, not Marc.
    assert len(meshioplusplus.read(PLT / "transient.dat").points) == 6


@pytest.mark.parametrize("engine", ENGINES)
def test_ordered_zones(engine):
    mesh = _read(engine, PLT / "ordered.plt")
    assert [(c.type, len(c.data)) for c in mesh.cells] == [
        ("line", 4),
        ("quad", 6),
        ("hexahedron", 4),
    ]
    assert [r.name for r in mesh.regions] == ["i_line", "ij_plane", "ijk_block"]
    p = mesh.points
    for block in mesh.cells[1:]:
        for c in block.data:
            a, b, d = p[c[0]], p[c[1]], p[c[3]]
            n = np.cross(b - a, d - a)
            if block.type == "quad":
                assert n[2] > 0
            else:
                assert np.dot(n, p[c[4]] - a) > 0
    # the IJ zone's first cell is (i, j) = (0, 0): its nodes are 0, 1, 5, 4 of the zone
    np.testing.assert_array_equal(mesh.cells[1].data[0] - 5, [0, 1, 5, 4])
    assert len(mesh.cell_data["C"][2]) == 4  # (3-1)*(3-1)*(2-1), ghosts dropped


@pytest.mark.parametrize("engine", ENGINES)
def test_sharing_and_passive(engine):
    mesh = _read(engine, PLT / "fe_2d.plt")
    assert [c.type for c in mesh.cells] == ["triangle", "quad", "quad", "line"]
    # zones 2 and 3 share X/Y with zone 1 but carry their own U: own points
    assert len(mesh.points) == 6 + 6 + 6 + 3
    np.testing.assert_array_equal(mesh.points[6:12], mesh.points[:6])
    # passive Q is NaN in the "quads" zone; the line zone stores Q at its nodes
    assert np.isnan(mesh.cell_data["Q"][1]).all()
    assert np.isnan(mesh.cell_data["Q"][3]).all()
    assert not np.isnan(mesh.point_data["Q"][18:]).any()
    assert np.isnan(mesh.point_data["Q"][:18]).all()
    # "quads_again" shares the connectivity of "quads"
    np.testing.assert_array_equal(mesh.cells[2].data - 12, mesh.cells[1].data - 6)


# --- hand-made files: data formats, byte order, refusals ---------------------------------


def _s(text, bo):
    return b"".join(struct.pack(bo + "i", ord(c)) for c in text) + struct.pack(
        bo + "i", 0
    )


def _write_plt(path, formats, zone_type=2, version=b"112", bo="<", values=None):
    """One FE triangle (or an FEPOLYGON header) with X, Y and a third variable in
    the given data formats."""
    i = lambda *v: struct.pack(bo + f"{len(v)}i", *v)  # noqa: E731
    out = (
        b"#!TDV"
        + version
        + i(1, 0)
        + _s("t", bo)
        + i(3)
        + _s("X", bo)
        + _s("Y", bo)
        + _s("V", bo)
    )
    out += (
        struct.pack(bo + "f", 299.0)
        + _s("z", bo)
        + i(-1, -1)
        + struct.pack(bo + "d", 0.0)
    )
    out += i(-1, zone_type, 0, 0, 0, 3)
    if zone_type in (6, 7):
        out += i(3, 6, 0, 0)
    out += i(1, 0, 0, 0, 0)
    out += struct.pack(bo + "f", 357.0) + struct.pack(bo + "f", 299.0)
    out += i(*formats) + i(0, 0, -1)
    values = values or [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [3.0, -2.0, 7.0]]
    out += struct.pack(bo + "6d", 0, 1, 0, 1, 0, 1)
    codes = {1: "f", 2: "d", 3: "i", 4: "h", 5: "B"}
    for fmt, vals in zip(formats, values):
        if fmt in codes:
            cast = int if fmt >= 3 else float
            out += struct.pack(bo + f"{len(vals)}" + codes[fmt], *map(cast, vals))
    out += i(0, 1, 2)
    path.write_bytes(out)
    return path


@pytest.mark.parametrize("engine", ENGINES)
@pytest.mark.parametrize("bo", ["<", ">"])
@pytest.mark.parametrize("fmt", [1, 2, 3, 4, 5])
def test_data_formats_and_byte_order(tmp_path, engine, bo, fmt):
    v = [3, 254, 7] if fmt == 5 else [3, -2, 7]  # a byte is unsigned
    path = _write_plt(
        tmp_path / "f.plt", [2, 1, fmt], bo=bo, values=[[0, 1, 0], [0, 0, 1], v]
    )
    mesh = _read(engine, path)
    np.testing.assert_array_equal(mesh.points, [[0, 0], [1, 0], [0, 1]])
    np.testing.assert_array_equal(mesh.point_data["V"], v)
    np.testing.assert_array_equal(mesh.cells[0].data, [[0, 1, 2]])


@pytest.mark.parametrize("engine", ENGINES)
def test_refusals(tmp_path, engine):
    with pytest.raises(ReadError, match="BIT"):
        _read(engine, _write_plt(tmp_path / "bit.plt", [2, 2, 6]))
    with pytest.raises(ReadError, match="version"):
        _read(engine, _write_plt(tmp_path / "old.plt", [2, 2, 2], version=b"102"))
    good = _write_plt(tmp_path / "good.plt", [2, 2, 2]).read_bytes()
    (tmp_path / "cut.plt").write_bytes(good[:-10])
    with pytest.raises(ReadError, match="truncated"):
        _read(engine, tmp_path / "cut.plt")


# --- face-based zones ---------------------------------------------------------------------


@pytest.mark.parametrize("engine", ENGINES)
def test_polygon_zones(engine):
    mesh = _read(engine, PLT / "poly_2d.plt")
    assert [(c.type, len(c.data)) for c in mesh.cells] == [("polygon", 3)] * 2
    # every ring counter-clockwise, starting at its first edge in the file
    assert _cell_lists(mesh.cells[0]) == [[0, 1, 4, 5], [4, 1, 2, 3, 6], [5, 4, 6, 7]]
    # the second zone shares X/Y and the face map but owns U: its own points
    assert _cell_lists(mesh.cells[1]) == [
        [v + 8 for v in r] for r in _cell_lists(mesh.cells[0])
    ]
    np.testing.assert_array_equal(mesh.points[8:], mesh.points[:8])
    assert [len(q) for q in mesh.cell_data["Q"]] == [3, 3]
    assert [r.name for r in mesh.regions] == ["polygons", "polygons_again"]


@pytest.mark.parametrize("engine", ENGINES)
def test_polyhedron_zones(engine):
    mesh = _read(engine, PLT / "poly_3d.plt")
    # one block per node count, in first-seen order: the zone's cells are
    # hex, tet, prism, hex -- so the second hex joins the first block
    assert [(c.type, len(c.data)) for c in mesh.cells] == [
        ("polyhedron8", 2),
        ("polyhedron4", 1),
        ("polyhedron10", 1),
        ("polyhedron8", 1),
    ]
    vols = [_volume(mesh.points, cell) for c in mesh.cells for cell in c.data]
    np.testing.assert_allclose(vols, [1, 1, 1 / 6, 2.5 * np.sin(0.4 * np.pi), 1])
    # cell-centred values follow their cells into the blocks
    dat = _read(engine, PLT / "poly_3d.dat")
    p = dat.cell_data["P"]
    assert [x.tolist() for x in p] == [x.tolist() for x in mesh.cell_data["P"]]
    zone_regions = {r.name: r.entries.tolist() for r in mesh.regions}
    assert zone_regions == {"polyhedra": [0, 1, 2, 3], "below": [4]}
    assert [x.tolist() for x in mesh.cell_data["tecplot:zone"]] == [
        [0, 0],
        [0],
        [0],
        [1],
    ]


@pytest.mark.skipif(not HAS_CORE, reason="needs the C++ core")
def test_polyhedron_metadata():
    meta = meshioplusplus.read_metadata(str(PLT / "poly_3d.plt"))
    assert [(b["type"], b["num_cells"]) for b in meta["cell_blocks"]] == [
        ("polyhedron8", 2),
        ("polyhedron4", 1),
        ("polyhedron10", 1),
        ("polyhedron8", 1),
    ]


@pytest.mark.parametrize("writer", ENGINES)
@pytest.mark.parametrize("reader", ENGINES)
@pytest.mark.parametrize("name", ["poly_2d", "poly_3d"])
def test_face_based_zones_round_trip(tmp_path, writer, reader, name):
    mesh = _read("python", PLT / f"{name}.plt")
    out = tmp_path / "out.dat"
    if writer == "core":
        _core.tecplot_write(str(out), mesh)
    else:
        _tecplot.write(str(out), mesh)
    text = out.read_text()
    assert ("FEPOLYHEDRON" if name == "poly_3d" else "FEPOLYGON") in text
    back = _read(reader, out)
    np.testing.assert_allclose(back.points, mesh.points)
    assert [c.type for c in back.cells] == [c.type for c in mesh.cells]
    for a, b in zip(back.cells, mesh.cells):
        assert _cell_lists(a) == _cell_lists(b)
    for a, b in zip(
        back.cell_data["P" if name == "poly_3d" else "Q"],
        mesh.cell_data["P" if name == "poly_3d" else "Q"],
    ):
        np.testing.assert_allclose(a, b)


@pytest.mark.parametrize("engine", ENGINES)
def test_face_based_refusals(tmp_path, engine):
    head = 'VARIABLES = "X" "Y"\n'
    square = "0 1 1 0\n0 0 1 1\n"
    # an open ring: three edges of a square
    open_ring = (
        head
        + "ZONE ZONETYPE=FEPOLYGON, NODES=4, ELEMENTS=1, FACES=3\n"
        + square
        + "1 2\n2 3\n3 4\n1 1 1\n0 0 0\n"
    )
    (tmp_path / "open.dat").write_text(open_ring)
    with pytest.raises(ReadError, match="closed polygon"):
        _read(engine, tmp_path / "open.dat")
    bad_node = (
        head
        + "ZONE ZONETYPE=FEPOLYGON, NODES=4, ELEMENTS=1, FACES=4\n"
        + square
        + "1 2\n2 3\n3 4\n4 9\n1 1 1 1\n0 0 0 0\n"
    )
    (tmp_path / "node.dat").write_text(bad_node)
    with pytest.raises(ReadError, match="out of range"):
        _read(engine, tmp_path / "node.dat")
    no_faces = head + "ZONE ZONETYPE=FEPOLYHEDRON, NODES=4, ELEMENTS=1\n" + square
    (tmp_path / "nofaces.dat").write_text(no_faces)
    with pytest.raises(ReadError, match="FACES"):
        _read(engine, tmp_path / "nofaces.dat")
    point = (
        head + "ZONE ZONETYPE=FEPOLYGON, NODES=4, ELEMENTS=1, FACES=4, "
        "DATAPACKING=POINT\n0 0\n1 0\n1 1\n0 1\n1 2\n2 3\n3 4\n4 1\n1 1 1 1\n0 0 0 0\n"
    )
    (tmp_path / "point.dat").write_text(point)
    with pytest.raises(ReadError, match="BLOCK"):
        _read(engine, tmp_path / "point.dat")


# --- ASCII ordered forms ------------------------------------------------------------------


@pytest.mark.parametrize("engine", ENGINES)
def test_ascii_ordered_forms(tmp_path, engine):
    path = tmp_path / "o.dat"
    path.write_text(
        'VARIABLES = "X" "Y" "P"\n'
        # old F=POINT form, J-ordered (I=1): a polyline along J
        'ZONE T="j", I=1, J=3, F=POINT\n'
        "0 0 1\n0 1 2\n0 2 3\n"
        # new form with POINT packing
        'ZONE T="ij", I=2, J=2, DATAPACKING=POINT\n'
        "0 0 4\n1 0 5\n0 1 6\n1 1 7\n"
        # VARSHARELIST without a zone number: the previous zone
        'ZONE T="ij2", I=2, J=2, VARSHARELIST=([1-2])\n'
        "8 9 10 11\n"
    )
    mesh = _read(engine, path)
    assert [(c.type, len(c.data)) for c in mesh.cells] == [
        ("line", 2),
        ("quad", 1),
        ("quad", 1),
    ]
    np.testing.assert_array_equal(mesh.cells[0].data, [[0, 1], [1, 2]])
    np.testing.assert_array_equal(mesh.cells[1].data, [[3, 4, 6, 5]])
    np.testing.assert_array_equal(
        mesh.point_data["P"], [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    )
    np.testing.assert_array_equal(mesh.points[7:], mesh.points[3:7])

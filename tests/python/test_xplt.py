"""FEBio plot files (``.xplt``): both engines, on plot files FEBio 4.12 wrote
for the models in ``tests/python/meshes/febio/xplt/`` (see
``tools/gen_febio_fixtures.py``) and on byte-swapped and truncated copies."""

import pathlib
import struct
import zlib

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _core
from meshioplusplus.xplt import _xplt as py_xplt

from .test_node_order import _is_valid

XPLT = pathlib.Path(__file__).parent / "meshes" / "febio" / "xplt"
FIXTURES = sorted(XPLT.glob("*.xplt"))


@pytest.fixture(params=["core", "python"])
def engine(request):
    """Both engines behind the same read signature."""
    if request.param == "core":
        return meshioplusplus.xplt
    return py_xplt


def _same(a, b):
    np.testing.assert_array_equal(a.points, b.points)
    assert [(c.type, c.data.tolist()) for c in a.cells] == [
        (c.type, c.data.tolist()) for c in b.cells
    ]
    regions = [
        sorted(
            (r.kind, r.name, r.dim, r.tag, np.asarray(r.entries).tolist())
            for r in m.regions
        )
        for m in (a, b)
    ]
    assert regions[0] == regions[1]
    for x, y in ((a.point_data, b.point_data), (a.field_data, b.field_data)):
        assert sorted(x) == sorted(y)
        for name in x:
            np.testing.assert_array_equal(x[name], y[name])
    assert sorted(a.cell_data) == sorted(b.cell_data)
    for name in a.cell_data:
        for x, y in zip(a.cell_data[name], b.cell_data[name]):
            np.testing.assert_array_equal(x, y)


@pytest.mark.parametrize("path", FIXTURES, ids=[p.name for p in FIXTURES])
def test_engines_read_every_state_alike(path):
    times = meshioplusplus.xplt.time_values(path)
    assert times == py_xplt.time_values(path)
    for k in range(len(times)):
        _same(
            meshioplusplus.xplt.read(path, time_step=k), py_xplt.read(path, time_step=k)
        )


def test_the_block_model(engine):
    mesh = engine.read(XPLT / "xplt_blocks.xplt", time_step=-1)
    # Two hexes, the shell, and the data surface "end" as a facet block.
    assert [c.type for c in mesh.cells] == ["hexahedron", "hexahedron", "quad", "quad"]
    assert [a.tolist() for a in mesh.cell_data["xplt:surface"]] == [[0], [0], [0], [1]]
    got = {(r.name, r.kind): r for r in mesh.regions}
    # Domains: solid ones by name, the nameless shell by its element set.
    assert got[("left", "cell")].tag == 1 and got[("right", "cell")].tag == 2
    assert got[("skin", "cell")].dim == 2
    assert got[("end", "side")].entries.tolist() == [[1, 1]]  # right hex, face x+
    fixed = got[("fixed", "point")].entries
    assert np.allclose(mesh.points[fixed][:, 0], 0)
    # The fixed nodes do not move; the loaded end does.
    u = mesh.point_data["displacement"]
    assert np.abs(u[fixed]).max() == 0
    assert u[:, 0].max() > 0
    assert mesh.cell_data["stress"][0].shape == (1, 6)
    assert mesh.point_data["nodal stress"].shape == (12, 6)
    assert mesh.field_data["meshio:time"].tolist() == [pytest.approx(1.0)]
    assert mesh.field_data["xplt:step"].tolist() == [3]
    # FEBio 4.12 writes the pressure-loaded surface's traction record empty:
    # a variable on no surface gives no array.
    assert "surface traction" not in mesh.cell_data
    assert got[("end", "cell")].entries.tolist() == [3]


def _facet_areas(points, rows):
    """Areas of planar quads as two triangles."""
    p = points[rows]
    a = np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
    b = np.cross(p[:, 2] - p[:, 0], p[:, 3] - p[:, 0])
    return 0.5 * (np.linalg.norm(a, axis=1) + np.linalg.norm(b, axis=1))


def test_surface_variables_ride_on_facet_blocks(engine):
    path = XPLT / "xplt_surface.xplt"
    for k in range(4):
        mesh = engine.read(path, time_step=k)
        # The bar, then one facet block per data surface: "end", "top".
        assert [(c.type, len(c.data)) for c in mesh.cells] == [
            ("hexahedron", 2),
            ("quad", 1),
            ("quad", 2),
        ]
        assert [a.tolist() for a in mesh.cell_data["xplt:surface"]] == [
            [0, 0],
            [1],
            [2, 2],
        ]
        got = {(r.name, r.kind): r for r in mesh.regions}
        assert got[("end", "cell")].entries.tolist() == [2]
        assert got[("top", "cell")].entries.tolist() == [3, 4]
        assert got[("end", "side")].entries.tolist() == [[1, 1]]
        # FEBio's facet areas are those of the deformed facets.
        deformed = mesh.points + mesh.point_data["displacement"]
        area = mesh.cell_data["facet area"]
        # FEBio 4.12 plots a variable on one surface only: "facet area" is
        # asked for on both, but only "top" has it.
        assert np.isnan(area[0]).all() and np.isnan(area[1]).all()
        np.testing.assert_allclose(
            area[2], _facet_areas(deformed, mesh.cells[2].data), rtol=1e-6
        )
        # A region variable holds one value per surface, on each facet.
        total = mesh.cell_data["surface area"]
        assert np.isnan(total[1]).all()
        np.testing.assert_allclose(total[2], area[2].sum(), rtol=1e-6)
    assert area[2][0] > 1.0  # the top was stretched


def test_edge_variables_ride_on_line_blocks(engine):
    # FEBio's one edge variable, "edge contact gap" (a value per edge node), on
    # the edge of an edge-to-surface sliding contact: three free nodes above
    # the bar, whose tip is bent up towards them.
    path = XPLT / "xplt_edge.xplt"
    gaps = []
    for k in range(4):
        mesh = engine.read(path, time_step=k)
        # The bar, the contact surface "top", then the edge "ridge".
        assert [(c.type, len(c.data)) for c in mesh.cells] == [
            ("hexahedron", 2),
            ("quad", 2),
            ("line", 2),
        ]
        np.testing.assert_array_equal(mesh.cells[2].data, [[12, 13], [13, 14]])
        assert [a.tolist() for a in mesh.cell_data["xplt:edge"]] == [
            [0, 0],
            [0, 0],
            [1, 1],
        ]
        got = {(r.name, r.kind): r for r in mesh.regions}
        assert got[("ridge", "cell")].entries.tolist() == [4, 5]
        assert got[("ridge", "cell")].dim == 1
        gap = mesh.point_data["edge contact gap"]
        assert np.isnan(gap[:12]).all()
        gaps.append(gap[12:])
    # No gap before the first solve; then the free nodes stay 0.25 above the
    # undeformed top face, the bar rising towards them, most at x=1.5.
    np.testing.assert_array_equal(gaps[0], 0)
    for g in gaps[1:]:
        assert (g < 0).all() and (-g <= 0.25).all()
        assert g[0] < g[1] < g[2]
    assert (np.diff(np.array(gaps[1:]), axis=0) > 0).all()


def test_remeshed_states_carry_their_own_mesh(engine):
    path = XPLT / "xplt_remesh.xplt"
    sizes = []
    for k, t in enumerate(py_xplt.time_values(path)):
        mesh = engine.read(path, time_step=k)
        sizes.append((len(mesh.points), len(mesh.cells[0].data)))
        assert mesh.cell_data["stress"][0].shape == (len(mesh.cells[0].data), 6)
        # The pulled end moved 0.1 t along x, however refined the mesh is.
        end = np.isclose(mesh.points[:, 0], 2)
        np.testing.assert_allclose(
            mesh.point_data["displacement"][end, 0], 0.1 * t, atol=1e-6
        )
    # FEBio's log: 12, 45, 78 and 111 nodes after each hex_refine pass.
    assert [n for n, _ in sizes] == [12, 45, 78, 111]
    assert sizes[1][1] == 16
    # Metadata describes the first mesh.
    info = meshioplusplus.read_metadata(path)
    assert info["time_values"] == pytest.approx([0, 1 / 3, 2 / 3, 1], abs=1e-6)
    assert len(list(meshioplusplus.read_sequence(path))) == 4


def test_compressed_file_reads_like_the_uncompressed_one():
    plain = meshioplusplus.xplt.read(XPLT / "xplt_blocks.xplt", time_step=2)
    packed = meshioplusplus.xplt.read(XPLT / "xplt_blocks_z.xplt", time_step=2)
    _same(plain, packed)


def test_hex27_is_in_meshio_order(engine):
    mesh = engine.read(XPLT / "xplt_hex27.xplt", time_step=-1)
    assert mesh.cells[0].type == "hexahedron27"
    row = mesh.cells[0].data[0]
    assert _is_valid(mesh.points[row], "hexahedron27")
    # Pushed down by 0.1 on top: the top face centre moved -0.1 in z.
    assert mesh.point_data["displacement"][row[25], 2] == pytest.approx(-0.1, abs=1e-6)


def test_sequences_and_metadata(tmp_path):
    path = XPLT / "xplt_blocks_z.xplt"
    assert meshioplusplus.read_metadata(path)["time_values"] == pytest.approx(
        [0, 1 / 3, 2 / 3, 1], abs=1e-6
    )
    steps = list(meshioplusplus.read_sequence(path))
    assert len(steps) == 4
    meshioplusplus.write_sequence(
        str(tmp_path / "step_{step}.vtu"), meshioplusplus.read_sequence(path)
    )
    written = sorted(tmp_path.glob("step_*.vtu"))
    assert len(written) == 4
    last = meshioplusplus.read(written[-1])
    np.testing.assert_allclose(
        last.point_data["displacement"], steps[-1][1].point_data["displacement"]
    )


def test_selective_reads(engine):
    path = XPLT / "xplt_blocks.xplt"
    only = engine.read(path, arrays=["stress"])
    assert list(only.cell_data) == ["stress"] and not only.point_data
    bare = engine.read(path, points_only=True)
    assert not bare.cell_data and not bare.point_data
    with pytest.raises(meshioplusplus.ReadError, match="out of range"):
        engine.read(path, time_step=4)


def _swapped(data):
    """The same (uncompressed) plot file in the other byte order: every 4-byte
    word swapped, except the bytes of strings and 64-byte names."""
    words = bytearray(data)

    def walk(start, stop):
        pos = start
        while pos + 8 <= stop:
            ident, size = struct.unpack_from("<II", data, pos)
            struct.pack_into(">II", words, pos, ident, size)
            a, b = pos + 8, pos + 8 + size
            if ident in _BRANCHES:
                walk(a, b)
            elif ident in _STRINGS:
                n = struct.unpack_from("<i", data, a)[0]
                struct.pack_into(">i", words, a, n)
            elif ident in _FIXED:
                pass
            else:
                for p in range(a, b, 4):
                    words[p : p + 4] = data[p : p + 4][::-1]
            pos = b

    struct.pack_into(">I", words, 0, py_xplt.MAGIC)
    walk(4, len(data))
    return bytes(words)


_BRANCHES = {
    py_xplt.ROOT,
    py_xplt.HEADER,
    py_xplt.DICTIONARY,
    py_xplt.DIC_ITEM,
    *py_xplt.DIC_GROUPS,
    py_xplt.MESH,
    py_xplt.NODE_SECTION,
    py_xplt.NODE_HEADER,
    py_xplt.DOMAIN_SECTION,
    py_xplt.DOMAIN,
    py_xplt.DOMAIN_HDR,
    py_xplt.DOM_ELEM_LIST,
    py_xplt.SURFACE_SECTION,
    py_xplt.SURFACE,
    py_xplt.SURFACE_HDR,
    py_xplt.FACE_LIST,
    py_xplt.NODESET_SECTION,
    py_xplt.NODESET,
    py_xplt.NODESET_HDR,
    py_xplt.PARTS_SECTION,
    py_xplt.PART,
    py_xplt.ELEMENTSET_SECTION,
    py_xplt.ELEMENTSET,
    py_xplt.ELEMENTSET_HDR,
    py_xplt.FACETSET_SECTION,
    py_xplt.FACETSET,
    py_xplt.FACETSET_HDR,
    py_xplt.FACETSET_LIST,
    py_xplt.STATE,
    py_xplt.STATE_HEADER,
    py_xplt.STATE_DATA,
    py_xplt.STATE_VARIABLE,
    # VAR_DATA holds (region id, size, floats) sub-chunks: walk it like a branch.
    py_xplt.STATE_VAR_DATA,
    *py_xplt.DATA_GROUPS,
    0x01048000,
    0x01048100,
    0x01048101,
    0x01048200,
    0x02030000,
}
_STRINGS = {
    0x01010005,
    0x01010006,
    py_xplt.DOM_NAME,
    py_xplt.SURFACE_NAME,
    py_xplt.NODESET_NAME,
    py_xplt.ELEMENTSET_NAME,
    py_xplt.FACETSET_NAME,
    0x01048104,
}
_FIXED = {py_xplt.DIC_ITEM_NAME, 0x01020006, 0x01020007, py_xplt.PART_NAME}


def test_byte_swapped_file_reads_the_same(engine, tmp_path):
    source = XPLT / "xplt_blocks.xplt"
    swapped = tmp_path / "swapped.xplt"
    swapped.write_bytes(_swapped(source.read_bytes()))
    assert meshioplusplus.sniff_format(swapped) == "xplt"
    _same(engine.read(swapped, time_step=-1), engine.read(source, time_step=-1))


def test_truncated_last_state_is_dropped(engine, tmp_path, capfd):
    source = (XPLT / "xplt_blocks_z.xplt").read_bytes()
    cut = tmp_path / "cut.xplt"
    cut.write_bytes(source[:-10])
    assert len(engine.time_values(cut)) == 3
    assert "truncated" in capfd.readouterr().err
    mesh = engine.read(cut, time_step=-1)
    assert mesh.field_data["xplt:step"].tolist() == [2]


def test_refuses_what_it_cannot_read(engine, tmp_path):
    junk = tmp_path / "junk.xplt"
    junk.write_bytes(b"not a plot file")
    with pytest.raises(meshioplusplus.ReadError, match="magic"):
        engine.read(junk)
    old = bytearray((XPLT / "xplt_hex27.xplt").read_bytes())
    # HDR_VERSION is the first leaf of the header: magic, root and header ids
    # and sizes (4 x 4 bytes), then its own id and size.
    version_at = 4 + 8 + 8 + 8
    assert struct.unpack_from("<I", old, version_at)[0] == 0x35
    struct.pack_into("<I", old, version_at, 0x08)
    junk.write_bytes(bytes(old))
    with pytest.raises(meshioplusplus.ReadError, match="0x0008"):
        engine.read(junk)


def test_zlib_streams_are_one_per_chunk():
    """The compressed fixture really is one zlib stream per state: inflating the
    first stream after the raw mesh gives a whole state chunk."""
    data = (XPLT / "xplt_blocks_z.xplt").read_bytes()
    pos = 4
    for _ in range(2):  # the raw root and mesh
        pos += 8 + struct.unpack_from("<I", data, pos + 4)[0]
    stream = zlib.decompressobj()
    chunk = stream.decompress(data[pos:])
    assert stream.eof and struct.unpack_from("<I", chunk, 0)[0] == py_xplt.STATE


def test_registered_and_sniffed():
    out = meshioplusplus.formats()
    assert out["extensions"][".xplt"] == ["xplt"]
    assert "xplt" in out["readable"] and "xplt" in out["writable"]
    for path in FIXTURES:
        assert meshioplusplus.sniff_format(path) == "xplt"
    assert hasattr(_core, "xplt_read")


# ---------------------------------------------------------------------------
# The writer: both engines write the same bytes, and either reader gets the
# mesh back (to the 32-bit floats the file holds).
# ---------------------------------------------------------------------------

from meshioplusplus.xplt import _write as py_write  # noqa: E402


@pytest.fixture(params=["core", "python"])
def writer(request):
    if request.param == "core":
        return lambda path, mesh: _core.xplt_write(str(path), mesh)
    return py_write.write


def _scratch_mesh():
    pts = np.array(
        [
            [0, 0, 0],
            [1, 0, 0],
            [1, 1, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 0, 1],
            [1, 1, 1],
            [0, 1, 1],
            [2, 0, 0],
            [2, 1, 0],
            [2, 0.5, 1],
        ],
        dtype=float,
    )
    cells = [
        ("hexahedron", np.array([[0, 1, 2, 3, 4, 5, 6, 7]])),
        ("tetra", np.array([[1, 8, 9, 5]])),
        ("triangle", np.array([[0, 1, 2]])),
        ("line", np.array([[0, 1], [1, 2]])),
    ]
    mesh = meshioplusplus.Mesh(pts, cells)
    mesh.point_data["displacement"] = np.arange(33.0).reshape(11, 3)
    mesh.point_data["temperature"] = np.arange(11.0)
    mesh.cell_data["pressure"] = [
        np.array([1.0]),
        np.array([2.0]),
        np.array([3.0]),
        np.array([4.0, 5.0]),
    ]
    mesh.cell_data["stress"] = [np.arange(6.0)[None] + 10 * k for k in range(3)] + [
        np.arange(12.0).reshape(2, 6)
    ]
    mesh.field_data["energy"] = np.array([4.5])
    mesh.field_data["meshio:time"] = np.array([0.5])
    mesh.regions = [
        meshioplusplus.Region("steel", "cell", np.array([0]), 3, 7),
        meshioplusplus.Region("pair", "cell", np.array([0, 1]), 3),
        meshioplusplus.Region("base", "point", np.array([0, 1, 2, 3])),
    ]
    return mesh


@pytest.mark.parametrize("path", FIXTURES, ids=[p.name for p in FIXTURES])
def test_the_engines_write_the_same_bytes(path, tmp_path):
    mesh = meshioplusplus.read(path)
    _core.xplt_write(str(tmp_path / "c.xplt"), mesh)
    py_write.write(tmp_path / "p.xplt", mesh)
    assert (tmp_path / "c.xplt").read_bytes() == (tmp_path / "p.xplt").read_bytes()


def test_a_scratch_mesh_survives_a_round_trip(writer, engine, tmp_path):
    mesh = _scratch_mesh()
    out = tmp_path / "scratch.xplt"
    writer(out, mesh)
    again = engine.read(out)
    np.testing.assert_allclose(again.points, mesh.points)
    assert [(c.type, c.data.tolist()) for c in again.cells] == [
        (c.type, c.data.tolist()) for c in mesh.cells
    ]
    np.testing.assert_allclose(
        again.point_data["displacement"], mesh.point_data["displacement"]
    )
    np.testing.assert_allclose(
        again.point_data["temperature"], mesh.point_data["temperature"]
    )
    for name in ("pressure", "stress"):
        for got, want in zip(again.cell_data[name], mesh.cell_data[name]):
            np.testing.assert_allclose(got, want)
    np.testing.assert_allclose(again.field_data["energy"], [4.5])
    assert float(again.field_data["meshio:time"][0]) == 0.5
    regions = {(r.kind, r.name): (r.tag, r.entries.tolist()) for r in again.regions}
    assert regions[("cell", "steel")] == (7, [0])
    assert regions[("cell", "pair")] == (-1, [0, 1])
    assert regions[("point", "base")] == (-1, [0, 1, 2, 3])


@pytest.mark.parametrize("path", FIXTURES, ids=[p.name for p in FIXTURES])
def test_a_febio_file_survives_a_round_trip(path, writer, engine, tmp_path):
    original = engine.read(path)
    out = tmp_path / "again.xplt"
    writer(out, original)
    again = engine.read(out)
    np.testing.assert_allclose(again.points, original.points)
    domains = [i for i, c in enumerate(original.cells)]
    assert [(c.type, c.data.tolist()) for c in again.cells[: len(domains)]] == [
        (c.type, c.data.tolist()) for c in original.cells[: len(domains)]
    ]
    for name, value in original.point_data.items():
        np.testing.assert_allclose(
            again.point_data[name], value, rtol=1e-6, equal_nan=True
        )


def test_vertices_and_unknown_cells_are_skipped(writer, engine, tmp_path):
    mesh = _scratch_mesh()
    mesh.cells.append(meshioplusplus.CellBlock("vertex", np.array([[0], [3]])))
    for name in mesh.cell_data:
        mesh.cell_data[name].append(np.zeros((2,) + mesh.cell_data[name][0].shape[1:]))
    out = tmp_path / "skip.xplt"
    writer(out, mesh)
    assert [c.type for c in engine.read(out).cells] == [
        "hexahedron",
        "tetra",
        "triangle",
        "line",
    ]


def test_a_series_is_a_state_per_step(tmp_path):
    steps = []
    for k in range(3):
        mesh = _scratch_mesh()
        mesh.point_data["temperature"] = mesh.point_data["temperature"] + k
        steps.append((0.5 * k, mesh))
    meshioplusplus.write_sequence(tmp_path / "run.xplt", iter(steps))
    assert meshioplusplus.xplt.time_values(tmp_path / "run.xplt") == [0.0, 0.5, 1.0]
    back = list(meshioplusplus.read_sequence(str(tmp_path / "run.xplt")))
    for (_, got), (_, mesh) in zip(back, steps):
        np.testing.assert_allclose(
            got.point_data["temperature"], mesh.point_data["temperature"]
        )
    # The Python twin writes the same bytes.
    with py_write.SeriesWriter(tmp_path / "py.xplt") as w:
        for time, mesh in steps:
            w.write(time, mesh)
    assert (tmp_path / "py.xplt").read_bytes() == (tmp_path / "run.xplt").read_bytes()


def test_a_series_refuses_a_different_mesh(tmp_path):
    other = _scratch_mesh()
    other.point_data["extra"] = np.zeros(11)
    for make in (meshioplusplus.xplt.SeriesWriter, py_write.SeriesWriter):
        with pytest.raises(meshioplusplus.WriteError):
            with make(tmp_path / "bad.xplt") as w:
                w.write(0.0, _scratch_mesh())
                w.write(1.0, other)

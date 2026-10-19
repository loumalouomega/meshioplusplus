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


# ---------------------------------------------------------------------------
# Surfaces, edges and facet sets: the inverse of what the reader builds.
# ---------------------------------------------------------------------------

SURFACE_FIXTURES = [XPLT / "xplt_surface.xplt", XPLT / "xplt_edge.xplt"]


def _lfs_pointer(path):
    return path.read_bytes()[:20].startswith(b"version https://git")


def _regions_of(mesh):
    return sorted(
        (r.kind, r.name, r.dim, r.tag, np.asarray(r.entries).tolist())
        for r in mesh.regions
    )


def _chunk_ids(path):
    """The top-level chunk ids and the mesh chunk's sections."""
    data = path.read_bytes()
    found = set()

    def walk(a, b):
        while a + 8 <= b:
            ident, size = struct.unpack_from("<II", data, a)
            found.add(ident)
            if ident == 0x01040000:
                walk(a + 8, a + 8 + size)
            a += 8 + size

    walk(4, len(data))
    return found


@pytest.mark.parametrize(
    "path", SURFACE_FIXTURES, ids=[p.name for p in SURFACE_FIXTURES]
)
def test_febio_surfaces_and_edges_survive_a_round_trip(path, writer, engine, tmp_path):
    if _lfs_pointer(path):
        pytest.skip("the FEBio fixture is a Git-LFS pointer")
    original = engine.read(path)
    assert "xplt:surface" in original.cell_data or "xplt:edge" in original.cell_data
    out = tmp_path / "again.xplt"
    writer(out, original)
    again = engine.read(out)
    np.testing.assert_array_equal(again.points, original.points)
    assert [(c.type, c.data.tolist()) for c in again.cells] == [
        (c.type, c.data.tolist()) for c in original.cells
    ]
    assert _regions_of(again) == _regions_of(original)
    for key in ("xplt:surface", "xplt:edge"):
        assert (key in again.cell_data) == (key in original.cell_data)
        for got, want in zip(
            again.cell_data.get(key, []), original.cell_data.get(key, [])
        ):
            np.testing.assert_array_equal(got, want)
    # The sections are really there, not folded into domains.
    ids = _chunk_ids(out)
    assert (0x01043000 in ids) == ("xplt:surface" in original.cell_data)
    assert (0x01048000 in ids) == ("xplt:edge" in original.cell_data)


def _surface_mesh():
    mesh = _scratch_mesh()
    mesh.cells[2] = meshioplusplus.CellBlock("triangle", np.array([[1, 8, 9]]))
    # block 2 (the triangle) is a data surface, block 3 (the lines) an edge
    mesh.cell_data["xplt:surface"] = [
        np.zeros(1, dtype=np.int64),
        np.zeros(1, dtype=np.int64),
        np.full(1, 7, dtype=np.int64),
        np.zeros(2, dtype=np.int64),
    ]
    mesh.cell_data["xplt:edge"] = [
        np.zeros(1, dtype=np.int64),
        np.zeros(1, dtype=np.int64),
        np.zeros(1, dtype=np.int64),
        np.full(2, 4, dtype=np.int64),
    ]
    del mesh.cell_data["pressure"], mesh.cell_data["stress"]
    mesh.regions += [
        meshioplusplus.Region("skin", "cell", np.array([2]), 2),
        meshioplusplus.Region("rim", "cell", np.array([3, 4]), 1),
        # the -x face of the hexahedron and the bottom face of the tetra
        meshioplusplus.Region("walls", "side", np.array([[0, 0], [0, 4], [1, 3]]), 2),
    ]
    return mesh


def test_a_surface_an_edge_and_a_facet_set_are_written(writer, engine, tmp_path):
    mesh = _surface_mesh()
    out = tmp_path / "facets.xplt"
    writer(out, mesh)
    ids = _chunk_ids(out)
    assert {0x01043000, 0x01047000, 0x01048000} <= ids
    again = engine.read(out)
    assert [(c.type, c.data.tolist()) for c in again.cells] == [
        (c.type, c.data.tolist()) for c in mesh.cells
    ]
    regions = {(r.kind, r.name): r for r in again.regions}
    assert regions[("cell", "skin")].entries.tolist() == [2]
    assert regions[("cell", "rim")].entries.tolist() == [3, 4]
    assert regions[("cell", "rim")].dim == 1
    # the surface "skin" lies on the tetra's face 3, which a side
    # region of the same name stands for; the facet set is its own region
    assert regions[("side", "skin")].entries.tolist() == [[1, 3]]
    assert sorted(regions[("side", "walls")].entries.tolist()) == sorted(
        mesh.regions[-1].entries.tolist()
    )
    assert again.cell_data["xplt:surface"][2].tolist() == [1]
    assert again.cell_data["xplt:edge"][3].tolist() == [1, 1]
    # neither marker is a variable, and the surface block is not a domain
    assert not any("xplt" in name for name in again.point_data)
    assert [r.name for r in again.regions if r.kind == "cell"].count("skin") == 1


def test_surfaces_need_the_marker_to_be_surfaces(writer, engine, tmp_path):
    mesh = _surface_mesh()
    del mesh.cell_data["xplt:surface"], mesh.cell_data["xplt:edge"]
    out = tmp_path / "plain.xplt"
    writer(out, mesh)
    assert not {0x01043000, 0x01048000} & _chunk_ids(out)
    assert 0x01047000 in _chunk_ids(out)
    again = engine.read(out)
    assert "xplt:surface" not in again.cell_data


def test_the_engines_write_the_same_surface_bytes(tmp_path):
    mesh = _surface_mesh()
    _core.xplt_write(str(tmp_path / "c.xplt"), mesh)
    py_write.write(tmp_path / "p.xplt", mesh)
    assert (tmp_path / "c.xplt").read_bytes() == (tmp_path / "p.xplt").read_bytes()
    with py_write.SeriesWriter(tmp_path / "ps.xplt") as w:
        w.write(0.0, mesh)
        w.write(1.0, mesh)
    meshioplusplus.write_sequence(
        tmp_path / "cs.xplt", iter([(0.0, mesh), (1.0, mesh)])
    )
    assert (tmp_path / "ps.xplt").read_bytes() == (tmp_path / "cs.xplt").read_bytes()


# ---------------------------------------------------------------------------
# Surface and edge variables, and compression.
# ---------------------------------------------------------------------------


def _dictionary_groups(path):
    """The dictionary's group ids and variable names, by group."""
    data = path.read_bytes()
    groups = {}
    root = struct.unpack_from("<II", data, 4)
    dictionary = None
    a = 12
    while a + 8 <= 12 + root[1]:
        ident, size = struct.unpack_from("<II", data, a)
        if ident == 0x01020000:
            dictionary = (a + 8, a + 8 + size)
        a += 8 + size
    a, end = dictionary
    while a + 8 <= end:
        gid, gsize = struct.unpack_from("<II", data, a)
        names, b = [], a + 8
        while b + 8 <= a + 8 + gsize:
            iid, isize = struct.unpack_from("<II", data, b)
            c = b + 8
            while c + 8 <= b + 8 + isize:
                cid, csize = struct.unpack_from("<II", data, c)
                if cid == 0x01020004:
                    names.append(
                        data[c + 8 : c + 8 + csize].split(b"\0", 1)[0].decode()
                    )
                c += 8 + csize
            b += 8 + isize
        groups[gid] = names
        a += 8 + gsize
    return groups


def _equal_cell_data(got, want, names):
    for name in names:
        assert name in got.cell_data, name
        for x, y in zip(got.cell_data[name], want.cell_data[name]):
            np.testing.assert_allclose(x, y, rtol=1e-6, equal_nan=True)


def test_febio_surface_variables_survive_a_round_trip(writer, engine, tmp_path):
    path = XPLT / "xplt_surface.xplt"
    if _lfs_pointer(path):
        pytest.skip("the FEBio fixture is a Git-LFS pointer")
    for step in range(4):
        original = engine.read(path, time_step=step)
        out = tmp_path / f"again{step}.xplt"
        writer(out, original)
        again = engine.read(out)
        _equal_cell_data(again, original, ["facet area", "surface area"])
        # they are surface variables, not domain variables
        groups = _dictionary_groups(out)
        assert "facet area" in groups[0x01025000]
        assert "facet area" not in groups.get(0x01024000, [])


def test_the_engines_write_the_same_bytes_with_surface_variables(tmp_path):
    original = meshioplusplus.read(XPLT / "xplt_surface.xplt", time_step=2)
    _core.xplt_write(str(tmp_path / "c.xplt"), original)
    py_write.write(tmp_path / "p.xplt", original)
    assert (tmp_path / "c.xplt").read_bytes() == (tmp_path / "p.xplt").read_bytes()


def _variable_mesh():
    """The surface mesh with a vector variable on the triangle surface, a scalar
    one on the edge and a scalar one shared with the domains."""
    mesh = _surface_mesh()
    nan = np.nan
    mesh.cell_data["traction"] = [
        np.full((1, 3), nan),
        np.full((1, 3), nan),
        np.array([[1.0, 2.0, 3.0]]),
        np.full((2, 3), nan),
    ]
    mesh.cell_data["gap"] = [
        np.full(1, nan),
        np.full(1, nan),
        np.full(1, nan),
        np.array([0.5, 0.25]),
    ]
    mesh.cell_data["shared"] = [
        np.array([1.0]),
        np.array([2.0]),
        np.array([3.0]),
        np.array([4.0, 5.0]),
    ]
    return mesh


def test_surface_and_edge_variables_are_written(writer, engine, tmp_path):
    mesh = _variable_mesh()
    out = tmp_path / "vars.xplt"
    writer(out, mesh)
    groups = _dictionary_groups(out)
    # "shared" has values on the domains, so it is a domain variable only
    assert groups[0x01025000] == ["traction"]
    assert groups[0x01026000] == ["gap"]
    assert groups[0x01024000] == ["shared"]
    again = engine.read(out)
    _equal_cell_data(again, mesh, ["traction", "gap"])
    np.testing.assert_allclose(
        np.concatenate(again.cell_data["shared"][:2]), [1.0, 2.0]
    )
    # the surface variable is NaN off its surface, the gap off its edge
    assert np.isnan(again.cell_data["traction"][0]).all()
    assert np.isnan(again.cell_data["gap"][2]).all()


def test_a_variable_that_only_the_domains_hold_stays_a_domain_variable(
    writer, engine, tmp_path
):
    mesh = _surface_mesh()
    mesh.cell_data["load"] = [
        np.array([1.0]),
        np.array([2.0]),
        np.full(1, np.nan),
        np.full(2, np.nan),
    ]
    out = tmp_path / "load.xplt"
    writer(out, mesh)
    groups = _dictionary_groups(out)
    assert groups[0x01024000] == ["load"]
    assert 0x01025000 not in groups and 0x01026000 not in groups


def test_the_engines_write_the_same_bytes_with_variables_and_series(tmp_path):
    mesh = _variable_mesh()
    _core.xplt_write(str(tmp_path / "c.xplt"), mesh)
    py_write.write(tmp_path / "p.xplt", mesh)
    assert (tmp_path / "c.xplt").read_bytes() == (tmp_path / "p.xplt").read_bytes()
    with py_write.SeriesWriter(tmp_path / "ps.xplt") as w:
        w.write(0.0, mesh)
        w.write(1.0, mesh)
    meshioplusplus.write_sequence(
        tmp_path / "cs.xplt", iter([(0.0, mesh), (1.0, mesh)])
    )
    assert (tmp_path / "ps.xplt").read_bytes() == (tmp_path / "cs.xplt").read_bytes()


def _core_has_zlib(tmp_path):
    try:
        _core.xplt_write(str(tmp_path / "probe.xplt"), _scratch_mesh(), True)
    except Exception:
        return False
    return True


def test_a_compressed_file_reads_like_the_plain_one(writer, engine, tmp_path):
    mesh = _variable_mesh()
    plain = tmp_path / "plain.xplt"
    packed = tmp_path / "packed.xplt"
    writer(plain, mesh)
    py_write.write(packed, mesh, compress=True)
    assert packed.read_bytes()[:4] == plain.read_bytes()[:4]
    # the header's compression flag is set and the state is one zlib stream
    data = packed.read_bytes()
    assert struct.unpack_from("<I", data, 40)[0] == 1
    pos = 12 + struct.unpack_from("<I", data, 8)[0]  # after the root chunk
    pos += 8 + struct.unpack_from("<I", data, pos + 4)[0]  # after the mesh chunk
    state = zlib.decompress(data[pos:])
    assert struct.unpack_from("<I", state, 0)[0] == 0x02000000
    _same(engine.read(packed), engine.read(plain))


def test_the_engines_write_the_same_compressed_bytes(tmp_path):
    if not _core_has_zlib(tmp_path):
        pytest.skip("the C++ core was built without zlib")
    mesh = _variable_mesh()
    _core.xplt_write(str(tmp_path / "c.xplt"), mesh, True)
    py_write.write(tmp_path / "p.xplt", mesh, compress=True)
    assert (tmp_path / "c.xplt").read_bytes() == (tmp_path / "p.xplt").read_bytes()
    with py_write.SeriesWriter(tmp_path / "ps.xplt", compress=True) as w:
        w.write(0.0, mesh)
        w.write(1.0, mesh)
    with meshioplusplus.xplt.SeriesWriter(tmp_path / "cs.xplt", compress=True) as w:
        w.write(0.0, mesh)
        w.write(1.0, mesh)
    assert (tmp_path / "ps.xplt").read_bytes() == (tmp_path / "cs.xplt").read_bytes()
    assert meshioplusplus.xplt.time_values(tmp_path / "cs.xplt") == [0.0, 1.0]


@pytest.mark.parametrize("path", FIXTURES, ids=[p.name for p in FIXTURES])
def test_a_febio_file_survives_a_compressed_round_trip(path, tmp_path):
    original = meshioplusplus.xplt.read(path)
    plain, packed = tmp_path / "plain.xplt", tmp_path / "packed.xplt"
    meshioplusplus.xplt.write(plain, original)
    meshioplusplus.xplt.write(packed, original, compress=True)
    _same(meshioplusplus.xplt.read(packed), meshioplusplus.xplt.read(plain))
    _same(py_xplt.read(packed), meshioplusplus.xplt.read(plain))


def test_compression_is_off_by_default(tmp_path):
    meshioplusplus.write(tmp_path / "a.xplt", _scratch_mesh())
    data = (tmp_path / "a.xplt").read_bytes()
    assert struct.unpack_from("<I", data, 40)[0] == 0


def _small_mesh():
    """A remesh of `_scratch_mesh`: fewer points and cells, the same variables."""
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
        ],
        dtype=float,
    )
    mesh = meshioplusplus.Mesh(
        pts, [("hexahedron", np.array([[0, 1, 2, 3, 4, 5, 6, 7]]))]
    )
    mesh.point_data["displacement"] = np.arange(24.0).reshape(8, 3)
    mesh.point_data["temperature"] = np.arange(8.0)
    mesh.cell_data["pressure"] = [np.array([9.0])]
    mesh.cell_data["stress"] = [np.arange(6.0)[None]]
    mesh.field_data["energy"] = np.array([7.5])
    return mesh


@pytest.mark.parametrize("compress", [False, True])
def test_a_remeshed_series_writes_a_mesh_section_where_the_mesh_changes(
    tmp_path, compress
):
    if compress and not meshioplusplus._core.__has_zlib__:
        pytest.skip("no zlib")
    steps = [(0.0, _scratch_mesh()), (1.0, _small_mesh()), (2.0, _small_mesh())]
    kwargs = {"compress": compress}
    with meshioplusplus.xplt.SeriesWriter(tmp_path / "core.xplt", **kwargs) as w:
        for time, mesh in steps:
            w.write(time, mesh)
    with py_write.SeriesWriter(tmp_path / "py.xplt", **kwargs) as w:
        for time, mesh in steps:
            w.write(time, mesh)
    assert (tmp_path / "core.xplt").read_bytes() == (tmp_path / "py.xplt").read_bytes()
    for engine_read in (meshioplusplus.xplt.read, py_xplt.read):
        for k, (_, mesh) in enumerate(steps):
            got = engine_read(tmp_path / "core.xplt", time_step=k)
            assert len(got.points) == len(mesh.points)
            assert [c.type for c in got.cells][0] == "hexahedron"
            np.testing.assert_allclose(
                got.point_data["temperature"], mesh.point_data["temperature"]
            )
    # The mesh section appears once per change: the third step reuses the second.
    plain = (tmp_path / "core.xplt").read_bytes()
    if not compress:
        assert plain.count(bytes.fromhex("00000401")) == 2


def test_a_remesh_with_the_same_counts_but_new_connectivity_is_a_new_mesh(tmp_path):
    first = _small_mesh()
    second = _small_mesh()
    second.cells[0].data[0] = [1, 0, 3, 2, 5, 4, 7, 6]
    for make in (meshioplusplus.xplt.SeriesWriter, py_write.SeriesWriter):
        out = tmp_path / f"{make.__module__}.xplt"
        with make(out) as w:
            w.write(0.0, first)
            w.write(1.0, second)
        got = py_xplt.read(out, time_step=1)
        assert got.cells[0].data[0].tolist() == [1, 0, 3, 2, 5, 4, 7, 6]


def _surface_node_mesh():
    """The surface mesh with a point array that only the triangle surface's nodes
    carry and one that only the edge's nodes carry, marked as the reader marks them."""
    mesh = _surface_mesh()
    contact = np.full(11, np.nan)
    contact[[1, 8, 9]] = [0.5, 1.5, 2.5]
    rim = np.full((11, 3), np.nan)
    rim[[0, 1, 2]] = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    mesh.point_data["contact"] = contact
    mesh.point_data["rim"] = rim
    mesh.field_data["xplt:surface_node:contact"] = np.array([1])
    mesh.field_data["xplt:edge_node:rim"] = np.array([1])
    return mesh


def test_surface_node_variables_return_to_their_group(writer, engine, tmp_path):
    mesh = _surface_node_mesh()
    out = tmp_path / "nodal.xplt"
    writer(out, mesh)
    groups = _dictionary_groups(out)
    assert groups[0x01025000] == ["contact"]
    assert groups[0x01026000] == ["rim"]
    # not nodal variables
    assert "contact" not in groups.get(0x01023000, [])
    again = engine.read(out)
    for name in ("contact", "rim"):
        np.testing.assert_allclose(
            again.point_data[name], mesh.point_data[name], equal_nan=True
        )
    # the markers come back, so a second round trip keeps the groups
    assert "xplt:surface_node:contact" in again.field_data
    assert "xplt:edge_node:rim" in again.field_data
    out2 = tmp_path / "nodal2.xplt"
    writer(out2, again)
    assert _dictionary_groups(out2) == groups


def test_the_engines_write_the_same_bytes_with_surface_node_variables(tmp_path):
    mesh = _surface_node_mesh()
    _core.xplt_write(str(tmp_path / "c.xplt"), mesh)
    py_write.write(tmp_path / "p.xplt", mesh)
    assert (tmp_path / "c.xplt").read_bytes() == (tmp_path / "p.xplt").read_bytes()


def test_a_marked_array_without_facets_stays_nodal(writer, engine, tmp_path):
    mesh = _scratch_mesh()
    mesh.field_data["xplt:surface_node:temperature"] = np.array([1])
    out = tmp_path / "plain.xplt"
    writer(out, mesh)
    assert "temperature" in _dictionary_groups(out)[0x01023000]

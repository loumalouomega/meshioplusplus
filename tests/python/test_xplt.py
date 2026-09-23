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
    assert [c.type for c in mesh.cells] == ["hexahedron", "hexahedron", "quad"]
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
    # The surface traction is surface data: not read.
    assert "surface traction" not in mesh.point_data


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
    assert "xplt" in out["readable"] and "xplt" not in out["writable"]
    for path in FIXTURES:
        assert meshioplusplus.sniff_format(path) == "xplt"
    assert hasattr(_core, "xplt_read")

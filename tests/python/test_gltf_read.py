"""Tests for the glTF 2.0 reader (``.glb`` / ``.gltf``).

Hand-built files exercise what the writer never emits (strips, fans, loops,
strides, sparse and normalized accessors, node hierarchies, ``data:`` URIs); the
writer's own output checks the round trip. The C++ core and the Python twin are
compared on every case.
"""

import base64
import json
import struct

import numpy as np
import pytest

import meshioplusplus as mio
from meshioplusplus import ReadError
from meshioplusplus.gltf import _gltf_read

CUBE_POINTS = np.array(
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
CUBE_QUADS = np.array(
    [[0, 3, 2, 1], [4, 5, 6, 7], [0, 1, 5, 4], [3, 7, 6, 2], [0, 4, 7, 3], [1, 2, 6, 5]]
)


def both(path, **kwargs):
    """Read with the dispatcher (core) and the Python twin; they must agree."""
    core = mio.gltf.read(str(path), **kwargs)
    twin = _gltf_read.read(str(path), **kwargs)
    assert core.points.shape == twin.points.shape
    np.testing.assert_allclose(core.points, twin.points, atol=1e-12)
    assert [b.type for b in core.cells] == [b.type for b in twin.cells]
    for a, b in zip(core.cells, twin.cells):
        np.testing.assert_array_equal(a.data, b.data)
    assert sorted(core.point_data) == sorted(twin.point_data)
    for name in core.point_data:
        np.testing.assert_allclose(core.point_data[name], twin.point_data[name])
    assert sorted(r.name for r in core.regions) == sorted(r.name for r in twin.regions)
    for r in core.regions:
        t = next(x for x in twin.regions if x.name == r.name)
        np.testing.assert_array_equal(r.entries, t.entries)
    return core


def gltf_file(tmp_path, binary, body, name="m.gltf"):
    uri = "data:application/octet-stream;base64," + base64.b64encode(binary).decode()
    doc = {"asset": {"version": "2.0"}}
    doc.update(body)
    doc["buffers"] = [{"byteLength": len(binary), "uri": uri}]
    path = tmp_path / name
    path.write_text(json.dumps(doc))
    return path


def test_round_trip_of_the_writers_output(tmp_path):
    mesh = mio.Mesh(CUBE_POINTS.copy(), [("quad", CUBE_QUADS.copy())])
    for suffix in (".glb", ".gltf"):
        path = tmp_path / ("cube" + suffix)
        mio.write(str(path), mesh)
        back = both(path, weld=True)
        assert len(back.points) == 8
        assert back.cells[0].type == "triangle" and len(back.cells[0].data) == 12
        np.testing.assert_allclose(
            np.sort(back.points, axis=0), np.sort(CUBE_POINTS, axis=0), atol=1e-6
        )
        # unwelded, the creased cube splits into 24 vertices
        assert len(both(path).points) == 24
        assert mio.read(str(path)).cells[0].type == "triangle"


@pytest.mark.parametrize("up", ["z", "y", "x"])
def test_the_writers_root_transform_is_undone(tmp_path, up):
    pts = CUBE_POINTS * [2.0, 3.0, 5.0] + [100.0, -50.0, 7.0]
    mesh = mio.Mesh(pts, [("quad", CUBE_QUADS.copy())])
    path = tmp_path / "t.glb"
    mio.gltf.write(str(path), mesh, up_axis=up, scale=0.5)
    back = both(path, weld=True)
    # scale is a unit change to metres, so the source coordinates come back scaled by 1
    np.testing.assert_allclose(
        np.sort(back.points, axis=0), np.sort(pts, axis=0), rtol=1e-6, atol=1e-4
    )


def test_regions_fields_and_colour(tmp_path):
    pts = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]], dtype=float)
    mesh = mio.Mesh(
        pts,
        [("triangle", np.array([[0, 1, 2], [1, 3, 2]]))],
        point_data={
            "temperature": np.array([1.0, 2.0, 3.0, 4.0]),
            "v": np.arange(8.0).reshape(4, 2),
        },
        regions=[mio.Region("left", "cell", [0]), mio.Region("right", "cell", [1])],
    )
    path = tmp_path / "f.glb"
    mio.gltf.write(str(path), mesh, color_by="temperature")
    back = both(path, weld=True)
    assert {r.name for r in back.regions} == {"left", "right"}
    assert set(back.point_data) >= {"temperature", "v", "gltf:color"}
    assert back.point_data["v"].shape == (4, 2)
    assert back.point_data["gltf:color"].shape == (4, 3)
    only = mio.gltf.read(str(path), arrays=["temperature"])
    assert list(only.point_data) == ["temperature"]
    assert not mio.gltf.read(str(path), points_only=True).point_data


def test_strips_fans_loops_strides_and_sparse(tmp_path):
    binary = b""
    for p in ([0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]):
        binary += struct.pack("<3f", *p) + b"\0\0\0\0"  # stride 16
    idx_off = len(binary)
    binary += struct.pack("<4H", 0, 1, 2, 3)
    sp_idx = len(binary)
    binary += struct.pack("<H", 3) + b"\0\0"
    sp_val = len(binary)
    binary += struct.pack("<f", 7.0)

    def body(mode):
        return {
            "scenes": [{"nodes": [0]}],
            "nodes": [{"name": "body", "mesh": 0}],
            "meshes": [
                {
                    "primitives": [
                        {
                            "attributes": {"POSITION": 0, "_S": 2},
                            "indices": 1,
                            "mode": mode,
                        }
                    ]
                }
            ],
            "accessors": [
                {"bufferView": 0, "componentType": 5126, "count": 4, "type": "VEC3"},
                {"bufferView": 1, "componentType": 5123, "count": 4, "type": "SCALAR"},
                {
                    "componentType": 5126,
                    "count": 4,
                    "type": "SCALAR",
                    "name": "s",
                    "sparse": {
                        "count": 1,
                        "indices": {"bufferView": 2, "componentType": 5123},
                        "values": {"bufferView": 3},
                    },
                },
            ],
            "bufferViews": [
                {"buffer": 0, "byteLength": idx_off, "byteStride": 16},
                {"buffer": 0, "byteOffset": idx_off, "byteLength": 8},
                {"buffer": 0, "byteOffset": sp_idx, "byteLength": 2},
                {"buffer": 0, "byteOffset": sp_val, "byteLength": 4},
            ],
        }

    expected = {
        5: ("triangle", 2),
        6: ("triangle", 2),
        2: ("line", 4),
        3: ("line", 3),
        0: ("vertex", 4),
    }
    for mode, (cell_type, count) in expected.items():
        path = gltf_file(tmp_path, binary, body(mode), f"m{mode}.gltf")
        mesh = both(path)
        assert [b.type for b in mesh.cells] == [cell_type]
        assert len(mesh.cells[0].data) == count
        np.testing.assert_array_equal(mesh.point_data["s"], [0, 0, 0, 7])
        assert [r.name for r in mesh.regions] == ["body"]
    strip = both(gltf_file(tmp_path, binary, body(5), "strip.gltf"))
    np.testing.assert_array_equal(strip.cells[0].data, [[0, 1, 2], [2, 1, 3]])
    fan = both(gltf_file(tmp_path, binary, body(6), "fan.gltf"))
    np.testing.assert_array_equal(fan.cells[0].data, [[0, 1, 2], [0, 2, 3]])


def test_node_hierarchy_and_normalized_quantized_positions(tmp_path):
    binary = bytes([0, 0, 0, 127, 0, 0, 0, 127, 0]) + b"\0\0\0"
    path = gltf_file(
        tmp_path,
        binary,
        {
            "nodes": [
                {"children": [1], "translation": [10, 0, 0]},
                {"mesh": 0, "scale": [2, 2, 2]},
            ],
            "meshes": [{"primitives": [{"attributes": {"POSITION": 0}}]}],
            "accessors": [
                {
                    "bufferView": 0,
                    "componentType": 5120,
                    "normalized": True,
                    "count": 3,
                    "type": "VEC3",
                }
            ],
            "bufferViews": [{"buffer": 0, "byteLength": 9}],
        },
    )
    mesh = both(path)
    np.testing.assert_allclose(
        mesh.points, [[10, 0, 0], [12, 0, 0], [10, 2, 0]], atol=1e-12
    )
    assert not mesh.regions


def test_a_glb_with_external_buffer_alongside(tmp_path):
    mesh = mio.Mesh(CUBE_POINTS.copy(), [("quad", CUBE_QUADS.copy())])
    path = tmp_path / "ext.gltf"
    mio.write(str(path), mesh)
    assert (tmp_path / "ext.bin").exists()
    assert len(both(path).points) == 24


def _bad(tmp_path, doc_extra, name="bad.gltf", binary=b"\0" * 12):
    return gltf_file(tmp_path, binary, doc_extra, name)


BASE = {
    "nodes": [{"mesh": 0}],
    "accessors": [{"bufferView": 0, "componentType": 5126, "count": 1, "type": "VEC3"}],
    "bufferViews": [{"buffer": 0, "byteLength": 12}],
}


def _with_prim(prim):
    doc = dict(BASE)
    doc["meshes"] = [{"primitives": [dict({"attributes": {"POSITION": 0}}, **prim)]}]
    return doc


@pytest.mark.parametrize(
    "prim, needle",
    [
        (
            {"extensions": {"KHR_draco_mesh_compression": {}}},
            "KHR_draco_mesh_compression",
        ),
        ({"indices": 0}, "indices"),
        ({"mode": 9}, "mode"),
    ],
)
def test_unsupported_and_malformed_primitives(tmp_path, prim, needle):
    path = _bad(tmp_path, _with_prim(prim))
    for reader in (lambda p: mio.gltf.read(p), lambda p: _gltf_read.read(p)):
        with pytest.raises(ReadError) as err:
            reader(str(path))
        assert needle in str(err.value)


def test_container_and_file_errors(tmp_path):
    p = tmp_path / "a.gltf"
    p.write_text("not json")
    with pytest.raises(ReadError):
        mio.gltf.read(str(p))
    p.write_text(
        json.dumps(
            {
                "asset": {"version": "2.0"},
                "buffers": [{"byteLength": 4, "uri": "../x.bin"}],
            }
        )
    )
    with pytest.raises(ReadError):
        mio.gltf.read(str(p))
    with pytest.raises(ReadError):
        _gltf_read.read(str(p))
    (tmp_path / "short.glb").write_bytes(b"glTF\x02\x00\x00\x00")
    with pytest.raises(ReadError):
        mio.gltf.read(str(tmp_path / "short.glb"))
    # a short buffer, a node cycle and a missing POSITION accessor
    doc = _with_prim({})
    doc["accessors"][0]["count"] = 2
    with pytest.raises(ReadError):
        mio.gltf.read(str(_bad(tmp_path, doc, "short.gltf")))
    cycle = {"nodes": [{"children": [0]}], "scenes": [{"nodes": [0]}]}
    with pytest.raises(ReadError):
        mio.gltf.read(str(_bad(tmp_path, cycle, "cycle.gltf")))

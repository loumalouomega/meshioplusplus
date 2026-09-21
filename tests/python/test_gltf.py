"""Tests for the glTF 2.0 writer (``.glb`` / ``.gltf`` + ``.bin``).

The files are decoded here with ``struct`` and ``json`` -- nothing in the writer
is trusted to check itself -- and the C++ core and the Python reference are
pinned against each other byte for byte. The Khronos glTF-Validator (the
``gltf-validator`` npm package) runs when ``node`` and the package are present;
set ``MESHIOPLUSPLUS_REQUIRE_GLTF_VALIDATOR=1`` to make its absence a failure
(the CI job does), and ``MESHIOPLUSPLUS_GLTF_VALIDATOR`` to the directory whose
``node_modules`` holds it.
"""

import json
import os
import shutil
import struct
import subprocess

import numpy as np
import pytest

import meshioplusplus as mio
from meshioplusplus import _colormap
from meshioplusplus.gltf import _gltf

try:
    from meshioplusplus import _core
except ImportError:  # pragma: no cover - a pure-Python build
    _core = None

needs_core = pytest.mark.skipif(_core is None, reason="needs the compiled core")

HERE = os.path.dirname(os.path.abspath(__file__))
VALIDATE = os.path.join(HERE, "gltf_validate.mjs")


# --------------------------------------------------------------------------- #
# fixtures                                                                     #
# --------------------------------------------------------------------------- #
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
    [
        [0, 3, 2, 1],
        [4, 5, 6, 7],
        [0, 1, 5, 4],
        [3, 7, 6, 2],
        [0, 4, 7, 3],
        [1, 2, 6, 5],
    ]
)


def cube_quads(**kw):
    return mio.Mesh(CUBE_POINTS.copy(), [("quad", CUBE_QUADS.copy())], **kw)


def cube_hex(**kw):
    return mio.Mesh(
        CUBE_POINTS.copy(), [("hexahedron", np.array([[0, 1, 2, 3, 4, 5, 6, 7]]))], **kw
    )


def tet_mesh(**kw):
    points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
    return mio.Mesh(points, [("tetra", np.array([[0, 1, 2, 3]]))], **kw)


def two_tets(**kw):
    points = np.array(
        [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]], dtype=float
    )
    return mio.Mesh(points, [("tetra", np.array([[0, 1, 2, 3], [1, 2, 3, 4]]))], **kw)


# --------------------------------------------------------------------------- #
# decoding                                                                     #
# --------------------------------------------------------------------------- #
def read_glb(path):
    """(json, bin) of a GLB, asserting the container rules on the way."""
    data = open(path, "rb").read()
    magic, version, length = struct.unpack_from("<III", data, 0)
    assert magic == 0x46546C67 and version == 2
    assert length == len(data), "the header length is the file length"
    off = 12
    chunks = []
    while off < len(data):
        clen, ctype = struct.unpack_from("<II", data, off)
        assert clen % 4 == 0, "chunks are 4-byte aligned"
        chunks.append((ctype, data[off + 8 : off + 8 + clen]))
        off += 8 + clen
    assert off == len(data)
    assert chunks[0][0] == 0x4E4F534A, "the JSON chunk comes first"
    raw = chunks[0][1]
    assert raw.rstrip(b" ") == raw.rstrip(), "the JSON chunk is padded with spaces"
    doc = json.loads(raw)
    bin_chunk = b""
    if len(chunks) > 1:
        assert chunks[1][0] == 0x004E4942 and len(chunks) == 2
        bin_chunk = chunks[1][1]
        assert (
            doc["buffers"][0]["byteLength"]
            <= len(bin_chunk)
            <= (doc["buffers"][0]["byteLength"] + 3)
        )
    return doc, bin_chunk


def accessor(doc, blob, index):
    acc = doc["accessors"][index]
    view = doc["bufferViews"][acc["bufferView"]]
    ncomp = {"SCALAR": 1, "VEC2": 2, "VEC3": 3, "VEC4": 4}[acc["type"]]
    dtype = {5125: "<u4", 5126: "<f4"}[acc["componentType"]]
    raw = blob[view["byteOffset"] : view["byteOffset"] + view["byteLength"]]
    arr = np.frombuffer(raw, dtype=dtype)
    assert arr.size == acc["count"] * ncomp
    return arr.reshape(-1, ncomp) if ncomp > 1 else arr


def primitives(doc):
    return [(m["name"], p) for m in doc["meshes"] for p in m["primitives"]]


def world_positions(doc, blob, prim):
    """Positions through the root node's TRS, as the viewer would apply them."""
    pos = accessor(doc, blob, prim["attributes"]["POSITION"]).astype(np.float64)
    root = doc["nodes"][0]
    s = root.get("scale", [1, 1, 1])
    pos = pos * np.array(s)
    x, y, z, w = root.get("rotation", [0, 0, 0, 1])
    rot = np.array(
        [
            [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
            [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
            [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)],
        ]
    )
    return pos @ rot.T + np.array(root["translation"])


def total_vertices(doc, blob):
    return sum(
        accessor(doc, blob, p["attributes"]["POSITION"]).shape[0]
        for _n, p in primitives(doc)
    )


# --------------------------------------------------------------------------- #
# the container                                                                #
# --------------------------------------------------------------------------- #
def test_glb_container_rules(tmp_path):
    path = tmp_path / "cube.glb"
    mio.write(str(path), cube_quads())
    doc, blob = read_glb(path)
    assert doc["asset"]["version"] == "2.0"
    assert doc["asset"]["generator"].startswith("Written by meshio++ v")
    assert doc["buffers"][0].get("uri") is None
    # every bufferView is inside the buffer and 4-byte aligned
    for view in doc["bufferViews"]:
        assert view["byteOffset"] % 4 == 0 and view["byteLength"] % 4 == 0
        assert (
            view["byteOffset"] + view["byteLength"] <= doc["buffers"][0]["byteLength"]
        )


def test_gltf_container_writes_a_sidecar(tmp_path):
    path = tmp_path / "my cube.gltf"
    mio.write(str(path), cube_quads())
    doc = json.loads(path.read_bytes())
    sidecar = tmp_path / "my cube.bin"
    assert sidecar.exists()
    assert doc["buffers"][0]["uri"] == "my%20cube.bin"
    assert doc["buffers"][0]["byteLength"] == sidecar.stat().st_size
    blob = sidecar.read_bytes()
    prim = doc["meshes"][0]["primitives"][0]
    assert accessor(doc, blob, prim["indices"]).size == 36

    # the container follows the suffix unless told otherwise
    forced = tmp_path / "forced.dat"
    mio.gltf.write(str(forced), cube_quads(), container="glb")
    read_glb(forced)


def test_the_format_is_write_only():
    formats = mio.formats()
    assert "gltf" in formats["writable"] and "gltf" not in formats["readable"]
    assert formats["extensions"][".glb"] == ["gltf"]
    assert formats["extensions"][".gltf"] == ["gltf"]
    with pytest.raises(Exception):
        mio.read("nothing.glb")


# --------------------------------------------------------------------------- #
# geometry                                                                     #
# --------------------------------------------------------------------------- #
def test_a_cube_of_quads(tmp_path):
    path = tmp_path / "cube.glb"
    mio.write(str(path), cube_quads())
    doc, blob = read_glb(path)
    ((_name, prim),) = primitives(doc)
    assert prim["mode"] == 4
    idx = accessor(doc, blob, prim["indices"])
    assert idx.size == 12 * 3
    pos = accessor(doc, blob, prim["attributes"]["POSITION"])
    normal = accessor(doc, blob, prim["attributes"]["NORMAL"])
    # a crease at every edge: 6 faces x 4 corners
    assert pos.shape == (24, 3) and normal.shape == (24, 3)
    np.testing.assert_allclose(np.linalg.norm(normal, axis=1), 1.0, atol=1e-6)
    assert np.all(np.abs(normal).max(axis=1) == 1.0)

    # POSITION carries min/max, and they are exact
    acc = doc["accessors"][prim["attributes"]["POSITION"]]
    assert acc["min"] == [float(x) for x in pos.min(axis=0)]
    assert acc["max"] == [float(x) for x in pos.max(axis=0)]

    # each triangle's winding agrees with its stored normal (outward)
    tri = pos[idx.reshape(-1, 3)]
    face = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    face /= np.linalg.norm(face, axis=1)[:, None]
    np.testing.assert_allclose(face, normal[idx.reshape(-1, 3)[:, 0]], atol=1e-6)


def test_the_source_coordinates_come_back_through_the_root_transform(tmp_path):
    points = CUBE_POINTS * 3.0 + np.array([1000.0, -2000.0, 500.0])
    mesh = mio.Mesh(points, [("quad", CUBE_QUADS)])
    path = tmp_path / "far.glb"
    mio.write(str(path), mesh)
    doc, blob = read_glb(path)
    ((_n, prim),) = primitives(doc)
    # z-up source: a -90 degree rotation about x, and small local coordinates
    assert doc["nodes"][0]["rotation"] == [
        -0.70710678118654757,
        0,
        0,
        0.70710678118654757,
    ]
    local = accessor(doc, blob, prim["attributes"]["POSITION"])
    assert np.abs(local).max() <= 1.5 + 1e-6, "recentred: the coordinates stay small"
    world = world_positions(doc, blob, prim)
    # the y-up world is (x, z, -y) of the source
    src = np.column_stack([points[:, 0], points[:, 2], -points[:, 1]])
    got = {tuple(np.round(w, 4)) for w in world}
    assert got == {tuple(np.round(s, 4)) for s in src}

    # recentre off keeps the coordinates and needs no translation
    path2 = tmp_path / "raw.glb"
    mio.gltf.write(str(path2), mesh, recenter=False)
    doc2, blob2 = read_glb(path2)
    raw = accessor(doc2, blob2, primitives(doc2)[0][1]["attributes"]["POSITION"])
    assert raw.max() >= 1000.0 or raw.min() <= -2000.0
    assert doc2["nodes"][0]["translation"] == [0, 0, 0]


@pytest.mark.parametrize("up", ["x", "y", "z"])
def test_up_axis_and_scale(tmp_path, up):
    mesh = cube_quads()
    path = tmp_path / f"up_{up}.glb"
    mio.gltf.write(str(path), mesh, up_axis=up, scale=0.001)
    doc, blob = read_glb(path)
    root = doc["nodes"][0]
    assert root["scale"] == [0.001] * 3
    ((_n, prim),) = primitives(doc)
    world = world_positions(doc, blob, prim)
    src = CUBE_POINTS * 0.001
    perm = {
        "z": lambda p: (p[0], p[2], -p[1]),
        "y": lambda p: tuple(p),
        "x": lambda p: (-p[1], p[0], p[2]),
    }[up]
    assert {tuple(np.round(w, 9)) for w in world} == {
        tuple(np.round(perm(s), 9)) for s in src
    }
    assert ("rotation" in root) == (up != "y")


def test_a_flat_mesh_is_y_up_with_no_rotation(tmp_path):
    mesh = mio.Mesh(
        np.array([[0, 0], [1, 0], [1, 1], [0, 1]], float),
        [("triangle", np.array([[0, 1, 2], [0, 2, 3]]))],
    )
    path = tmp_path / "flat.glb"
    mio.write(str(path), mesh)
    doc, blob = read_glb(path)
    assert "rotation" not in doc["nodes"][0]
    ((_n, prim),) = primitives(doc)
    normal = accessor(doc, blob, prim["attributes"]["NORMAL"])
    np.testing.assert_array_equal(normal, [[0, 0, 1]] * 4)  # smooth: no crease


def test_volume_blocks_export_their_skin(tmp_path):
    for mesh, verts in ((tet_mesh(), 12), (cube_hex(), 24)):
        path = tmp_path / "v.glb"
        mio.write(str(path), mesh)
        doc, blob = read_glb(path)
        assert total_vertices(doc, blob) == verts


def test_interior_faces_are_not_exported(tmp_path):
    path = tmp_path / "two.glb"
    mio.write(str(path), two_tets())
    doc, blob = read_glb(path)
    ((_n, prim),) = primitives(doc)
    # 8 facets minus the shared one twice over = 6 boundary triangles
    assert accessor(doc, blob, prim["indices"]).size == 6 * 3


def test_a_smoother_split_angle_shares_vertices(tmp_path):
    mesh = mio.Mesh(*_sphere_arrays())
    counts = {}
    for angle in (0.0, 30.0, 180.0):
        path = tmp_path / f"s{int(angle)}.glb"
        mio.gltf.write(str(path), mesh, split_angle=angle)
        doc, blob = read_glb(path)
        counts[angle] = total_vertices(doc, blob)
    assert counts[180.0] == len(mesh.points)  # smooth: one vertex per point
    assert counts[30.0] == counts[180.0]  # a sphere has no crease at 30 degrees
    assert counts[0.0] > counts[180.0]  # every face on its own


def _sphere_arrays():
    from .test_curvature import icosphere

    m = icosphere(1, 1.0)
    return m.points, [("triangle", m.cells[0].data)]


def test_normals_can_be_switched_off(tmp_path):
    path = tmp_path / "n.glb"
    mio.gltf.write(str(path), cube_quads(), normals=False)
    doc, blob = read_glb(path)
    ((_n, prim),) = primitives(doc)
    assert "NORMAL" not in prim["attributes"]
    assert accessor(doc, blob, prim["attributes"]["POSITION"]).shape[0] == 8


# --------------------------------------------------------------------------- #
# regions                                                                      #
# --------------------------------------------------------------------------- #
def test_one_node_per_region_and_a_boundary_patch_wins_over_the_skin(tmp_path):
    points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
    mesh = mio.Mesh(
        points,
        [("tetra", np.array([[0, 1, 2, 3]])), ("triangle", np.array([[0, 1, 2]]))],
        regions=[
            mio.Region("fluid", "cell", [0]),
            mio.Region("wall", "cell", [1]),
        ],
    )
    path = tmp_path / "r.glb"
    mio.write(str(path), mesh)
    doc, blob = read_glb(path)
    names = [m["name"] for m in doc["meshes"]]
    assert names == ["fluid", "wall"]
    assert [n["name"] for n in doc["nodes"][1:]] == ["fluid", "wall"]
    tri_counts = {
        m["name"]: accessor(doc, blob, m["primitives"][0]["indices"]).size // 3
        for m in doc["meshes"]
    }
    # the skin has 4 facets; the wall cell replaces one of them
    assert tri_counts == {"fluid": 3, "wall": 1}


def test_a_cell_in_several_regions_goes_to_the_smallest_and_leftovers_are_unassigned(
    tmp_path,
):
    mesh = mio.Mesh(
        np.array([[0, 0], [1, 0], [2, 0], [0, 1], [1, 1], [2, 1]], float),
        [
            ("quad", np.array([[0, 1, 4, 3], [1, 2, 5, 4]])),
            ("triangle", np.array([[0, 1, 4]])),
        ],
        regions=[
            mio.Region("everything", "cell", [0, 1]),
            mio.Region("left", "cell", [0]),
        ],
    )
    path = tmp_path / "r2.glb"
    mio.write(str(path), mesh)
    doc, blob = read_glb(path)
    counts = {
        m["name"]: accessor(doc, blob, m["primitives"][0]["indices"]).size // 3
        for m in doc["meshes"]
    }
    assert counts == {"everything": 2, "left": 2, "unassigned": 1}
    assert [m["name"] for m in doc["meshes"]] == ["everything", "left", "unassigned"]

    flat = tmp_path / "one.glb"
    mio.gltf.write(str(flat), mesh, by_region=False)
    doc1, _b = read_glb(flat)
    assert [m["name"] for m in doc1["meshes"]] == ["mesh"]


# --------------------------------------------------------------------------- #
# lines and points                                                             #
# --------------------------------------------------------------------------- #
def test_lines_and_points_and_point_clouds(tmp_path):
    lines = mio.Mesh(
        np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0]], float),
        [("line", np.array([[0, 1], [1, 2]]))],
    )
    p = tmp_path / "l.glb"
    mio.write(str(p), lines)
    doc, blob = read_glb(p)
    ((_n, prim),) = primitives(doc)
    assert prim["mode"] == 1
    assert accessor(doc, blob, prim["indices"]).size == 4
    assert "NORMAL" not in prim["attributes"]

    # a cell-less mesh -- a point cloud -- is one POINTS primitive, normals kept
    cloud = mio.Mesh(
        np.array([[0, 0, 0], [1, 0, 0], [0, 2, 0]], float),
        [],
        point_data={"normals": np.array([[0, 0, 2.0], [0, 0, 1.0], [0, 3.0, 0]])},
    )
    p2 = tmp_path / "c.glb"
    mio.write(str(p2), cloud)
    doc2, blob2 = read_glb(p2)
    ((_n2, prim2),) = primitives(doc2)
    assert prim2["mode"] == 0 and "indices" not in prim2
    normal = accessor(doc2, blob2, prim2["attributes"]["NORMAL"])
    np.testing.assert_array_equal(normal, [[0, 0, 1], [0, 0, 1], [0, 1, 0]])

    # a zero normal row: NORMAL is omitted, not invalid
    cloud.point_data["normals"][1] = 0.0
    p3 = tmp_path / "c2.glb"
    mio.write(str(p3), cloud)
    doc3, _b3 = read_glb(p3)
    assert "NORMAL" not in primitives(doc3)[0][1]["attributes"]

    verts = mio.Mesh(
        np.array([[0, 0, 0], [1, 1, 1]], float), [("vertex", np.array([[0], [1]]))]
    )
    p4 = tmp_path / "v.glb"
    mio.write(str(p4), verts)
    assert primitives(read_glb(p4)[0])[0][1]["mode"] == 0


def test_an_empty_mesh_is_a_valid_empty_scene(tmp_path):
    mesh = mio.Mesh(np.zeros((0, 3)), [])
    path = tmp_path / "empty.glb"
    mio.write(str(path), mesh)
    doc, blob = read_glb(path)
    assert doc["scenes"] == [{}] and "nodes" not in doc and blob == b""


# --------------------------------------------------------------------------- #
# fields and colour                                                            #
# --------------------------------------------------------------------------- #
def test_point_data_is_exported_as_raw_custom_attributes(tmp_path):
    mesh = cube_quads(
        point_data={
            "temperature": np.arange(8.0),
            "vel": np.arange(24.0).reshape(8, 3),
            "count": np.arange(8, dtype=np.int64),
            "T": np.ones(8),
            "t": np.ones(8) * 2,
            "bad": np.array([np.nan] + [1.0] * 7),
            "tensor": np.zeros((8, 6)),
        }
    )
    path = tmp_path / "f.glb"
    mio.write(str(path), mesh)
    doc, blob = read_glb(path)
    ((_n, prim),) = primitives(doc)
    attrs = prim["attributes"]
    assert {"_TEMPERATURE", "_VEL", "_COUNT", "_T", "_T_2"} <= set(attrs)
    assert "_BAD" not in attrs and "_TENSOR" not in attrs
    pos = accessor(doc, blob, attrs["POSITION"])
    temp = accessor(doc, blob, attrs["_TEMPERATURE"])
    # the value at a vertex is the value at its point
    for v in range(len(pos)):
        p = int(np.flatnonzero((CUBE_POINTS == pos[v] + 0.5).all(axis=1))[0])
        assert temp[v] == p
    assert doc["accessors"][attrs["_TEMPERATURE"]]["name"] == "temperature"
    assert doc["accessors"][attrs["_VEL"]]["type"] == "VEC3"

    off = tmp_path / "off.glb"
    mio.gltf.write(str(off), mesh, fields=False)
    assert not any(
        k.startswith("_") for k in primitives(read_glb(off)[0])[0][1]["attributes"]
    )


def _linear_color(t, cmap="viridis"):
    table = _colormap.colormap_table(cmap)
    rgb = _colormap.colormap_lookup(table, t)
    lin = np.array(_colormap.SRGB_TO_LINEAR_BITS, dtype=np.uint32).view(np.float32)
    return [float(lin[c]) for c in rgb]


def test_colouring_by_a_point_field(tmp_path):
    mesh = cube_quads(point_data={"T": np.arange(8.0)})
    path = tmp_path / "c.glb"
    mio.gltf.write(str(path), mesh, color_by="T")
    doc, blob = read_glb(path)
    ((_n, prim),) = primitives(doc)
    assert "KHR_materials_unlit" in doc["extensionsUsed"]
    assert "KHR_materials_unlit" in doc["materials"][0]["extensions"]
    assert doc["materials"][0]["pbrMetallicRoughness"]["baseColorFactor"] == [
        1,
        1,
        1,
        1,
    ]
    color = accessor(doc, blob, prim["attributes"]["COLOR_0"])
    pos = accessor(doc, blob, prim["attributes"]["POSITION"])
    for v in range(len(pos)):
        p = int(np.flatnonzero((CUBE_POINTS == pos[v] + 0.5).all(axis=1))[0])
        np.testing.assert_array_equal(
            color[v], np.array(_linear_color(p / 7.0), dtype=np.float32)
        )
    # the raw values ride along
    assert "_T" in prim["attributes"]

    lit = tmp_path / "lit.glb"
    mio.gltf.write(str(lit), mesh, color_by="T", unlit=False)
    assert "extensionsUsed" not in read_glb(lit)[0]


def test_colouring_by_a_cell_field_splits_vertices_per_cell(tmp_path):
    mesh = mio.Mesh(
        np.array([[0, 0], [1, 0], [2, 0], [0, 1], [1, 1], [2, 1]], float),
        [("quad", np.array([[0, 1, 4, 3], [1, 2, 5, 4]]))],
        cell_data={"c": [np.array([0.0, 1.0])]},
    )
    path = tmp_path / "cc.glb"
    mio.gltf.write(str(path), mesh, color_by="c", vmin=0.0, vmax=1.0)
    doc, blob = read_glb(path)
    ((_n, prim),) = primitives(doc)
    color = accessor(doc, blob, prim["attributes"]["COLOR_0"])
    assert color.shape[0] == 8, "the shared edge carries one colour per cell"
    np.testing.assert_array_equal(color[0], np.array(_linear_color(0.0), np.float32))
    np.testing.assert_array_equal(color[-1], np.array(_linear_color(1.0), np.float32))


def test_a_non_finite_value_takes_the_nan_colour(tmp_path):
    mesh = cube_quads(point_data={"T": np.array([np.nan] + [1.0] * 7)})
    path = tmp_path / "nan.glb"
    mio.gltf.write(str(path), mesh, color_by="T", nan_color="#ff0000", vmin=0, vmax=2)
    doc, blob = read_glb(path)
    ((_n, prim),) = primitives(doc)
    color = accessor(doc, blob, prim["attributes"]["COLOR_0"])
    pos = accessor(doc, blob, prim["attributes"]["POSITION"])
    red = np.array([1.0, 0.0, 0.0], np.float32)
    for v in range(len(pos)):
        if (pos[v] + 0.5 == 0).all():
            np.testing.assert_array_equal(color[v], red)


def test_colour_option_errors():
    m = cube_quads(point_data={"T": np.arange(8.0)})
    for kw, match in (
        (dict(color_by="nope"), "neither point_data nor cell_data"),
        (dict(color_by="T", cmap="nope"), "unknown colormap"),
        (dict(color_by="T", nan_color="red"), "#rrggbb"),
        (dict(color_by="T", vmin=2, vmax=1), "vmin"),
        (dict(color_by="T", component=3), "out of range"),
    ):
        with pytest.raises(ValueError, match=match):
            mio.gltf.write("x.glb", m, **kw)


def test_option_errors():
    m = cube_quads()
    for kw, match in (
        (dict(split_angle=181), "split_angle"),
        (dict(split_angle=-1), "split_angle"),
        (dict(scale=0), "scale"),
        (dict(up_axis="w"), "up axis"),
        (dict(container="zip"), "container"),
        (dict(normal_weight="huge"), "weight"),
    ):
        with pytest.raises(ValueError, match=match):
            mio.gltf.write("x.glb", m, **kw)


def test_a_non_finite_coordinate_is_refused(tmp_path):
    m = cube_quads()
    m.points[3, 1] = np.inf
    with pytest.raises(Exception, match="non-finite"):
        mio.write(str(tmp_path / "bad.glb"), m)


def test_writing_to_a_buffer(tmp_path):
    import io

    buf = io.BytesIO()
    mio.gltf.write(buf, cube_quads())
    path = tmp_path / "b.glb"
    path.write_bytes(buf.getvalue())
    read_glb(path)
    with pytest.raises(ValueError, match="binary"):
        mio.gltf.write(io.BytesIO(), cube_quads(), container="gltf")


# --------------------------------------------------------------------------- #
# the C++ core and the Python reference write the same bytes                   #
# --------------------------------------------------------------------------- #
def _parity_cases():
    from .test_curvature import icosphere, open_cylinder

    regioned = mio.Mesh(
        np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], float),
        [("tetra", np.array([[0, 1, 2, 3]])), ("triangle", np.array([[0, 1, 2]]))],
        point_data={"T": np.arange(4.0), "v": np.arange(12.0).reshape(4, 3)},
        cell_data={"c": [np.array([3.0]), np.array([5.0])]},
        regions=[mio.Region("fluid", "cell", [0]), mio.Region("wall", "cell", [1])],
    )
    ragged = mio.Mesh(
        np.array(
            [[0, 0, 0], [1, 0, 0], [1.5, 1, 0], [0.5, 2, 0], [-0.5, 1, 0], [0, 0, 1]],
            float,
        ),
        [("polygon", [[0, 1, 2, 3, 4]]), ("triangle", np.array([[0, 5, 4]]))],
    )
    lines = mio.Mesh(
        np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [5, 5, 5]], float),
        [("line", np.array([[0, 1], [1, 2]])), ("vertex", np.array([[3]]))],
    )
    cloud = mio.Mesh(
        np.random.default_rng(3).random((20, 3)),
        [],
        point_data={
            "normals": np.tile([0.0, 0.0, 1.0], (20, 1)),
            "p": np.linspace(0, 1, 20),
        },
    )
    return {
        "cube": cube_quads(point_data={"T": np.arange(8.0)}),
        "hex": cube_hex(point_data={"T": np.arange(8.0)}),
        "two_tets": two_tets(),
        "sphere": icosphere(2, 1.5),
        "cylinder": open_cylinder(),
        "regioned": regioned,
        "ragged": ragged,
        "lines": lines,
        "cloud": cloud,
    }


_OPTION_SETS = [
    {},
    {"normal_weight": "area", "split_angle": 60.0},
    {"normals": False, "fields": False, "by_region": False, "recenter": False},
    {"up_axis": "x", "scale": 0.001, "split_angle": 180.0},
    {"color_by": "T", "cmap": "turbo"},
    {"color_by": "c", "unlit": False, "vmin": 0.0, "vmax": 10.0},
]


@needs_core
@pytest.mark.parametrize("name", list(_parity_cases()))
@pytest.mark.parametrize("opts", _OPTION_SETS, ids=lambda o: ",".join(o) or "default")
@pytest.mark.parametrize("suffix", [".glb", ".gltf"])
def test_the_core_and_the_reference_write_the_same_bytes(tmp_path, name, opts, suffix):
    mesh = _parity_cases()[name]
    if ("color_by" in opts) and opts["color_by"] not in list(mesh.point_data) + list(
        mesh.cell_data
    ):
        pytest.skip("this fixture has no such array")
    a, b = tmp_path / f"core{suffix}", tmp_path / f"py{suffix}"
    _core.gltf_write(
        str(a),
        mesh,
        **{
            "container": "auto",
            "up_axis": opts.get("up_axis", "auto"),
            "normal_weight": opts.get("normal_weight", "angle"),
            "normals": opts.get("normals", True),
            "fields": opts.get("fields", True),
            "recenter": opts.get("recenter", True),
            "by_region": opts.get("by_region", True),
            "unlit": opts.get("unlit", True),
            "split_angle": opts.get("split_angle", 30.0),
            "scale": opts.get("scale", 1.0),
            "color_by": opts.get("color_by", ""),
            "cmap": opts.get("cmap", "viridis"),
            "vmin": opts.get("vmin"),
            "vmax": opts.get("vmax"),
        },
    )
    _gltf.write(str(b), mesh, **opts)
    ja, jb = a.read_bytes(), b.read_bytes()
    if suffix == ".gltf":
        # the sidecar's basename differs only by the file's own name
        ja, jb = ja.replace(b"core.bin", b"X.bin"), jb.replace(b"py.bin", b"X.bin")
        assert a.with_suffix(".bin").read_bytes() == b.with_suffix(".bin").read_bytes()
    assert ja == jb


# --------------------------------------------------------------------------- #
# the tables and provenance                                                    #
# --------------------------------------------------------------------------- #
@needs_core
def test_the_srgb_table_matches_the_core():
    assert list(_core.srgb_to_linear_table()) == list(_colormap.SRGB_TO_LINEAR_BITS)
    lin = np.array(_colormap.SRGB_TO_LINEAR_BITS, dtype=np.uint32).view(np.float32)
    assert lin[0] == 0.0 and lin[255] == 1.0
    assert abs(float(lin[128]) - 0.21586050) < 1e-7
    assert np.all(np.diff(lin) > 0)


def test_generator_carries_the_provenance_tag(tmp_path):
    from meshioplusplus import _provenance

    path = tmp_path / "p.glb"
    mio.write(str(path), cube_quads())
    generator = read_glb(path)[0]["asset"]["generator"]
    assert generator == _provenance.lines(_provenance.SlotTier.BLOCK)[0]
    assert generator.startswith("Written by meshio++ v")


# --------------------------------------------------------------------------- #
# external checks                                                              #
# --------------------------------------------------------------------------- #
def _validator_available():
    if shutil.which("node") is None:
        return False
    root = os.environ.get("MESHIOPLUSPLUS_GLTF_VALIDATOR", os.getcwd())
    return os.path.isdir(os.path.join(root, "node_modules", "gltf-validator"))


def _require_validator():
    if _validator_available():
        return
    if os.environ.get("MESHIOPLUSPLUS_REQUIRE_GLTF_VALIDATOR") == "1":
        pytest.fail("the glTF-Validator is required but node/gltf-validator is missing")
    pytest.skip("needs node and the gltf-validator npm package")


def _validate(path):
    proc = subprocess.run(
        ["node", VALIDATE, str(path)], capture_output=True, text=True, timeout=120
    )
    report = json.loads(proc.stdout.strip().splitlines()[-1]) if proc.stdout else {}
    assert (
        proc.returncode == 0 and report.get("numErrors") == 0
    ), f"{path.name}: {report.get('messages') or proc.stderr}"
    assert report["numWarnings"] == 0, report["messages"]
    return report


@pytest.mark.parametrize("name", list(_parity_cases()))
@pytest.mark.parametrize(
    "opts", _OPTION_SETS[:1] + _OPTION_SETS[3:5], ids=lambda o: ",".join(o) or "default"
)
@pytest.mark.parametrize("suffix", [".glb", ".gltf"])
def test_the_khronos_validator_accepts_the_output(tmp_path, name, opts, suffix):
    _require_validator()
    mesh = _parity_cases()[name]
    if ("color_by" in opts) and opts["color_by"] not in list(mesh.point_data) + list(
        mesh.cell_data
    ):
        pytest.skip("this fixture has no such array")
    path = tmp_path / f"m{suffix}"
    mio.gltf.write(str(path), mesh, **opts)
    _validate(path)


def test_the_khronos_validator_accepts_an_empty_scene(tmp_path):
    _require_validator()
    path = tmp_path / "empty.glb"
    mio.write(str(path), mio.Mesh(np.zeros((0, 3)), []))
    _validate(path)


def test_trimesh_loads_the_output(tmp_path):
    trimesh = pytest.importorskip("trimesh")
    path = tmp_path / "t.glb"
    mio.write(str(path), cube_hex())
    scene = trimesh.load(str(path), force="scene")
    geoms = list(scene.geometry.values())
    assert sum(len(g.faces) for g in geoms) == 12
    assert scene.bounds is not None
    # trimesh applies the node transform: y-up world of the z-up source
    lo, hi = scene.bounds
    np.testing.assert_allclose(hi - lo, [1, 1, 1], atol=1e-5)

    # a tall z-up part stands up in the y-up world
    tall = mio.Mesh(CUBE_POINTS * np.array([1.0, 1.0, 5.0]), [("quad", CUBE_QUADS)])
    path2 = tmp_path / "tall.glb"
    mio.write(str(path2), tall)
    lo, hi = trimesh.load(str(path2), force="scene").bounds
    np.testing.assert_allclose(hi - lo, [1, 5, 1], atol=1e-5)


# --------------------------------------------------------------------------- #
# the convert CLI                                                              #
# --------------------------------------------------------------------------- #
def test_cli_convert_to_glb_with_colour(tmp_path):
    src = tmp_path / "in.vtu"
    mio.write(str(src), cube_quads(point_data={"T": np.arange(8.0)}))
    dst = tmp_path / "out.glb"
    mio._cli.main(
        [
            "convert",
            str(src),
            str(dst),
            "--color-by",
            "T",
            "--cmap",
            "coolwarm",
            "--vmin",
            "0",
            "--vmax",
            "7",
            "--split-angle",
            "60",
            "--up-axis",
            "y",
        ]
    )
    doc, _blob = read_glb(dst)
    ((_n, prim),) = primitives(doc)
    assert "COLOR_0" in prim["attributes"]
    assert "rotation" not in doc["nodes"][0]

    plain = tmp_path / "plain.glb"
    mio._cli.main(["convert", str(src), str(plain)])
    assert "COLOR_0" not in primitives(read_glb(plain)[0])[0][1]["attributes"]


def test_cli_rejects_flags_that_do_not_apply(tmp_path):
    src = tmp_path / "in.vtu"
    mio.write(str(src), cube_quads(point_data={"T": np.arange(8.0)}))
    for extra, match in (
        (["--colorbar", "--color-by", "T"], "colorbar"),
        (["--split-angle", "10"], "glTF"),
    ):
        out = tmp_path / ("x.glb" if "--colorbar" in extra else "x.vtu")
        with pytest.raises(ValueError, match=match):
            mio._cli.main(["convert", str(src), str(out), *extra])

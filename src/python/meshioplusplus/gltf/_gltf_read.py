"""glTF 2.0 reader -- the Python twin of ``formats/gltf_read.cpp``.

A general reader, not one limited to meshio++'s own output: ``.glb`` and
``.gltf`` (external or ``data:`` buffers), every accessor layout (component
type, ``normalized``, ``byteStride``, sparse), node transforms applied down the
hierarchy, and the writer's root node (up axis, scale, recentring) undone.
``POINTS`` become ``vertex`` cells, ``LINES``/``LINE_STRIP``/``LINE_LOOP``
``line`` cells and ``TRIANGLES``/``TRIANGLE_STRIP``/``TRIANGLE_FAN``
``triangle`` cells, in that block order reversed (triangle, line, vertex).
A node name becomes a cell region, ``_NAME`` attributes become ``point_data``,
``COLOR_0`` becomes ``gltf:color`` and a ``POINTS`` primitive's ``NORMAL``
becomes ``normals``. Draco and meshopt are refused by name.
"""

from __future__ import annotations

import base64
import binascii
import json
import math
import os
import struct
from urllib.parse import unquote

import numpy as np

from .._common import warn
from .._exceptions import ReadError
from .._mesh import CellBlock, Mesh
from .._regions import Region

_MAX_COUNT = 1 << 27
_MAX_ZEROS = 1 << 26
_MAX_INT = 1 << 40
_HALF_SQRT2 = 0.70710678118654757

_DTYPES = {
    5120: np.dtype("<i1"),
    5121: np.dtype("<u1"),
    5122: np.dtype("<i2"),
    5123: np.dtype("<u2"),
    5125: np.dtype("<u4"),
    5126: np.dtype("<f4"),
}
_TYPES = {
    "SCALAR": 1,
    "VEC2": 2,
    "VEC3": 3,
    "VEC4": 4,
    "MAT2": 4,
    "MAT3": 9,
    "MAT4": 16,
}


def _fail(message):
    raise ReadError("gltf: " + message)


# ---- the JSON accessors ---------------------------------------------------- #


def _object(value, what):
    if not isinstance(value, dict):
        _fail(f"{what} is not an object")
    return value


def _uint(value, what):
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or float(value) != int(value)
        or value < 0
        or value > _MAX_INT
    ):
        _fail(f"{what} is not a non-negative integer")
    return int(value)


def _uint_or(obj, key, default, what):
    if key not in obj:
        return default
    return _uint(obj[key], f"{what}.{key}")


def _list(obj, key):
    value = obj.get(key) if isinstance(obj, dict) else None
    if value is None:
        return []
    if not isinstance(value, list):
        _fail(f"'{key}' is not an array")
    return value


def _numbers(value, count, what):
    if not isinstance(value, list) or len(value) != count:
        _fail(f"{what} is not an array of {count} numbers")
    out = []
    for v in value:
        if (
            isinstance(v, bool)
            or not isinstance(v, (int, float))
            or not math.isfinite(v)
        ):
            _fail(f"{what} holds a value that is not a finite number")
        out.append(float(v))
    return out


# ---- the container --------------------------------------------------------- #


def _load(filename):
    directory = (
        os.path.dirname(os.fspath(filename)) if not hasattr(filename, "read") else ""
    )
    if hasattr(filename, "read"):
        data = filename.read()
    else:
        try:
            with open(filename, "rb") as f:
                data = f.read()
        except OSError as exc:
            raise ReadError(f"gltf: cannot open '{filename}'") from exc
    if isinstance(data, str):
        data = data.encode("utf-8")
    binary = None
    if data[:4] == b"glTF":
        if len(data) < 20:
            _fail("the GLB is truncated")
        version, total = struct.unpack_from("<II", data, 4)
        if version != 2:
            _fail("only GLB version 2 is supported")
        if total > len(data):
            _fail("the GLB header declares more bytes than the file holds")
        pos = 12
        json_bytes = None
        while pos + 8 <= total:
            length, kind = struct.unpack_from("<II", data, pos)
            pos += 8
            if length > total - pos:
                _fail("a GLB chunk runs past the end of the file")
            if json_bytes is None:
                if kind != 0x4E4F534A:
                    _fail("the first GLB chunk is not JSON")
                json_bytes = data[pos : pos + length]
            elif kind == 0x004E4942 and binary is None:
                binary = data[pos : pos + length]
            pos += length + (-length) % 4
        if json_bytes is None:
            _fail("the GLB has no JSON chunk")
        text = json_bytes.rstrip(b" \x00")
    else:
        name = os.fspath(filename) if not hasattr(filename, "read") else ""
        if name.lower().endswith(".glb"):
            _fail(f"'{name}' has the .glb suffix but not the GLB magic")
        text = data
    try:
        root = json.loads(text.decode("utf-8-sig"))
    except (ValueError, UnicodeDecodeError) as exc:
        raise ReadError(f"gltf: invalid JSON ({exc})") from exc
    if not isinstance(root, dict):
        _fail("the document is not a JSON object")
    asset = root.get("asset")
    version = asset.get("version") if isinstance(asset, dict) else None
    if not isinstance(version, str) or not version.startswith("2."):
        _fail("not a glTF 2.x document (asset.version)")
    return root, binary, directory


def _uri_bytes(uri, directory):
    if uri.startswith("data:"):
        comma = uri.find(",")
        header = "" if comma < 0 else uri[5:comma]
        if not header.endswith(";base64"):
            _fail("a data: URI that is not base64")
        try:
            return base64.b64decode(uri[comma + 1 :] + "===", validate=False)
        except (binascii.Error, ValueError) as exc:
            raise ReadError("gltf: a data: URI is not valid base64") from exc
    path = unquote(uri)
    if not path or path[0] in "/\\" or "://" in path or "\\" in path or ":" in path:
        _fail(f"the buffer URI '{uri}' is not a relative file path")
    if any(part == ".." for part in path.split("/")):
        _fail(f"the buffer URI '{uri}' leaves the file's directory")
    try:
        with open(os.path.join(directory, path), "rb") as f:
            return f.read()
    except OSError as exc:
        raise ReadError(f"gltf: cannot open '{path}'") from exc


def _buffers(root, binary, directory):
    out = []
    for i, b in enumerate(_list(root, "buffers")):
        b = _object(b, f"buffers[{i}]")
        length = _uint_or(b, "byteLength", 0, "buffer")
        uri = b.get("uri")
        if uri is None:
            if i != 0 or binary is None:
                _fail(f"buffer {i} has no uri and no GLB BIN chunk")
            data = binary
        else:
            if not isinstance(uri, str):
                _fail("a buffer uri is not a string")
            data = _uri_bytes(uri, directory)
        if len(data) < length:
            _fail(f"buffer {i} is shorter than its byteLength")
        out.append(data)
    return out


# ---- accessors ------------------------------------------------------------- #


def _span(ctx, view, byte_offset, elem_size, count, use_view_stride):
    views = _list(ctx.root, "bufferViews")
    if view >= len(views):
        _fail("a bufferView index is out of range")
    v = _object(views[view], f"bufferViews[{view}]")
    ext = v.get("extensions")
    if isinstance(ext, dict) and "EXT_meshopt_compression" in ext:
        _fail("EXT_meshopt_compression is not supported")
    buffer = _uint_or(v, "buffer", -1, "bufferView")
    if buffer < 0 or buffer >= len(ctx.buffers):
        _fail("a bufferView names a buffer that does not exist")
    offset = _uint_or(v, "byteOffset", 0, "bufferView")
    length = _uint_or(v, "byteLength", 0, "bufferView")
    if offset + length > len(ctx.buffers[buffer]):
        _fail("a bufferView runs past the end of its buffer")
    stride = elem_size
    if use_view_stride:
        s = _uint_or(v, "byteStride", 0, "bufferView")
        if s:
            if s < elem_size or s > (1 << 20):
                _fail("a bufferView byteStride is out of range")
            stride = s
    if count > 0 and byte_offset + (count - 1) * stride + elem_size > length:
        _fail("an accessor runs past the end of its bufferView")
    return ctx.buffers[buffer], offset + byte_offset, stride


def _normalize(values, component_type):
    if component_type == 5120:
        return np.maximum(values / 127.0, -1.0)
    if component_type == 5121:
        return values / 255.0
    if component_type == 5122:
        return np.maximum(values / 32767.0, -1.0)
    if component_type == 5123:
        return values / 65535.0
    return values


def _strided(buf, offset, stride, count, ncomp, dtype):
    if count == 0:
        return np.zeros((0, ncomp), dtype=np.float64)
    arr = np.ndarray(
        shape=(count, ncomp),
        dtype=dtype,
        buffer=buf,
        offset=offset,
        strides=(stride, dtype.itemsize),
    )
    return arr.astype(np.float64)


class _Accessor:
    __slots__ = ("count", "ncomp", "component_type", "name", "values")


def _read_accessor(ctx, index):
    accessors = _list(ctx.root, "accessors")
    if index >= len(accessors):
        _fail("an accessor index is out of range")
    a = _object(accessors[index], f"accessors[{index}]")
    out = _Accessor()
    count = _uint_or(a, "count", 0, "accessor")
    if count > _MAX_COUNT:
        _fail("an accessor count is too large")
    out.count = count
    kind = a.get("type")
    if not isinstance(kind, str):
        _fail("an accessor has no type")
    if kind not in _TYPES:
        _fail(f"unknown accessor type '{kind}'")
    out.ncomp = _TYPES[kind]
    out.component_type = _uint_or(a, "componentType", 0, "accessor")
    if out.component_type not in _DTYPES:
        _fail(f"unknown componentType {out.component_type}")
    dtype = _DTYPES[out.component_type]
    normalized = a.get("normalized") is True
    name = a.get("name")
    out.name = name if isinstance(name, str) else ""
    if kind.startswith("MAT"):
        _fail("matrix accessors are not supported as vertex data")

    n = count * out.ncomp
    if "bufferView" not in a and n > _MAX_ZEROS:
        _fail("an accessor without a bufferView is too large")
    values = np.zeros((count, out.ncomp), dtype=np.float64)
    if "bufferView" in a:
        elem = dtype.itemsize * out.ncomp
        buf, offset, stride = _span(
            ctx,
            _uint(a["bufferView"], "accessor.bufferView"),
            _uint_or(a, "byteOffset", 0, "accessor"),
            elem,
            count,
            True,
        )
        values = _strided(buf, offset, stride, count, out.ncomp, dtype)
        if normalized:
            values = _normalize(values, out.component_type)
    sparse = a.get("sparse")
    if sparse is not None:
        sparse = _object(sparse, "sparse")
        scount = _uint_or(sparse, "count", 0, "sparse")
        if scount > count:
            _fail("a sparse accessor has more entries than elements")
        idx = _object(sparse.get("indices"), "sparse.indices")
        val = _object(sparse.get("values"), "sparse.values")
        itype = _uint_or(idx, "componentType", 0, "sparse.indices")
        if itype not in (5121, 5123, 5125):
            _fail("sparse indices must be unsigned")
        idtype = _DTYPES[itype]
        ib, ioff, _ = _span(
            ctx,
            _uint_or(idx, "bufferView", -1, "sparse.indices"),
            _uint_or(idx, "byteOffset", 0, "sparse.indices"),
            idtype.itemsize,
            scount,
            False,
        )
        elem = dtype.itemsize * out.ncomp
        vb, voff, _ = _span(
            ctx,
            _uint_or(val, "bufferView", -1, "sparse.values"),
            _uint_or(val, "byteOffset", 0, "sparse.values"),
            elem,
            scount,
            False,
        )
        at = _strided(ib, ioff, idtype.itemsize, scount, 1, idtype)[:, 0]
        if at.size and at.max() >= count:
            _fail("a sparse index is out of range")
        sv = _strided(vb, voff, elem, scount, out.ncomp, dtype)
        if normalized:
            sv = _normalize(sv, out.component_type)
        values = values.copy()
        values[at.astype(np.int64)] = sv
    out.values = values
    return out


# ---- transforms ------------------------------------------------------------ #


def _node_matrix(node):
    """Row-major 4x4 matrix of a node (the transpose of glTF's column-major)."""
    if "matrix" in node:
        v = _numbers(node["matrix"], 16, "node.matrix")
        return np.array(v, dtype=np.float64).reshape(4, 4).T
    t = [0.0, 0.0, 0.0]
    s = [1.0, 1.0, 1.0]
    q = [0.0, 0.0, 0.0, 1.0]
    if "translation" in node:
        t = _numbers(node["translation"], 3, "node.translation")
    if "scale" in node:
        s = _numbers(node["scale"], 3, "node.scale")
    if "rotation" in node:
        q = _numbers(node["rotation"], 4, "node.rotation")
    norm = math.sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3])
    if not norm > 0.0:
        _fail("a node rotation is the zero quaternion")
    x, y, z, w = (c / norm for c in q)
    m = np.eye(4)
    m[0, 0] = (1 - 2 * (y * y + z * z)) * s[0]
    m[1, 0] = (2 * (x * y + z * w)) * s[0]
    m[2, 0] = (2 * (x * z - y * w)) * s[0]
    m[0, 1] = (2 * (x * y - z * w)) * s[1]
    m[1, 1] = (1 - 2 * (x * x + z * z)) * s[1]
    m[2, 1] = (2 * (y * z + x * w)) * s[1]
    m[0, 2] = (2 * (x * z + y * w)) * s[2]
    m[1, 2] = (2 * (y * z - x * w)) * s[2]
    m[2, 2] = (1 - 2 * (x * x + y * y)) * s[2]
    m[0, 3], m[1, 3], m[2, 3] = t
    return m


def _writer_root(node):
    """The source-space offset of the writer's root node ``meshio++``, or None."""
    if node.get("name") != "meshio++" or "mesh" in node or "matrix" in node:
        return None
    t = [0.0, 0.0, 0.0]
    s = 1.0
    up = 2  # 0 = Z up, 1 = X up, 2 = Y up (no rotation)
    if "rotation" in node:
        q = _numbers(node["rotation"], 4, "node.rotation")

        def close(a, b, c, d):
            tol = 1e-9
            same = all(abs(x - y) < tol for x, y in zip(q, (a, b, c, d)))
            flip = all(abs(x + y) < tol for x, y in zip(q, (a, b, c, d)))
            return same or flip

        if close(-_HALF_SQRT2, 0, 0, _HALF_SQRT2):
            up = 0
        elif close(0, 0, _HALF_SQRT2, _HALF_SQRT2):
            up = 1
        elif not close(0, 0, 0, 1):
            return None
    if "scale" in node:
        n = _numbers(node["scale"], 3, "node.scale")
        if not n[0] > 0.0 or n[1] != n[0] or n[2] != n[0]:
            return None
        s = n[0]
    if "translation" in node:
        t = _numbers(node["translation"], 3, "node.translation")
    if up == 0:
        offset = [t[0] / s, -t[2] / s, t[1] / s]
    elif up == 1:
        offset = [t[1] / s, -t[0] / s, t[2] / s]
    else:
        offset = [t[0] / s, t[1] / s, t[2] / s]
    return [o + 0.0 for o in offset]


# ---- the mesh -------------------------------------------------------------- #


class _State:
    def __init__(self, want_data, arrays):
        self.points = []  # list of (n, 3) arrays
        self.npoints = 0
        self.cells = ([], [], [])  # triangle, line, vertex: lists of int64 rows
        self.cell_counts = [0, 0, 0]
        self.regions = {}  # name -> ([tri ids], [line ids], [vertex ids])
        self.fields = {}  # name -> [ncomp, [chunks], dropped]
        self.warned = set()
        self.want_data = want_data
        self.arrays = arrays

    def wants(self, name):
        return self.arrays is None or name in self.arrays


def _field_name(attr, accessor, mode):
    if attr == "COLOR_0":
        return "gltf:color"
    if attr == "NORMAL":
        return "normals" if mode == 0 else ""
    if len(attr) > 1 and attr[0] == "_":
        return accessor.name or attr[1:]
    return ""


def _primitive(ctx, prim, world, offset, region, state):
    ext = prim.get("extensions")
    if isinstance(ext, dict) and "KHR_draco_mesh_compression" in ext:
        _fail("KHR_draco_mesh_compression is not supported")
    mode = _uint_or(prim, "mode", 4, "primitive")
    if mode > 6:
        _fail(f"unknown primitive mode {mode}")
    attrs = _object(prim.get("attributes"), "primitive.attributes")
    if "POSITION" not in attrs:
        if "noposition" not in state.warned:
            state.warned.add("noposition")
            warn("gltf: a primitive without POSITION is skipped.")
        return
    pos = _read_accessor(ctx, _uint(attrs["POSITION"], "POSITION"))
    if pos.ncomp != 3:
        _fail("POSITION is not VEC3")
    nv = pos.count
    base = state.npoints
    if base + nv > _MAX_COUNT:
        _fail("the file holds too many vertices")

    if "indices" in prim:
        acc = _read_accessor(ctx, _uint(prim["indices"], "indices"))
        if acc.ncomp != 1 or acc.component_type not in (5121, 5123, 5125):
            _fail("indices must be unsigned SCALAR")
        flat = acc.values[:, 0]
        if flat.size and flat.max() >= nv:
            _fail("an index is out of range")
        index = flat.astype(np.int64)
    else:
        index = np.arange(nv, dtype=np.int64)

    p = pos.values @ world[:3, :3].T + world[:3, 3] + np.asarray(offset)
    if not np.all(np.isfinite(p)):
        _fail("a vertex position is not finite")
    state.points.append(p + 0.0)
    state.npoints += nv

    local = {}
    if state.want_data:
        for key, ref in attrs.items():
            if key == "POSITION" or isinstance(ref, bool):
                continue
            if not isinstance(ref, (int, float)):
                continue
            if not key or (key[0] != "_" and key not in ("COLOR_0", "NORMAL")):
                continue
            acc = _read_accessor(ctx, _uint(ref, key))
            if acc.count != nv:
                _fail(f"attribute {key} does not have one value per vertex")
            name = _field_name(key, acc, mode)
            if not name or acc.ncomp > 4 or not state.wants(name):
                continue
            if name in local:
                if "dup:" + name not in state.warned:
                    state.warned.add("dup:" + name)
                    warn(f"gltf: two attributes map to '{name}'; the first is kept.")
                continue
            local[name] = acc.values
    for name, values in local.items():
        entry = state.fields.get(name)
        if entry is None:
            entry = [values.shape[1], [], False]
            if base:
                entry[1].append(np.full((base, values.shape[1]), np.nan))
            state.fields[name] = entry
        elif entry[0] != values.shape[1] and not entry[2]:
            entry[2] = True
            warn(
                f"gltf: attribute '{name}' changes its component count between "
                "primitives; dropped."
            )
    for name, entry in state.fields.items():
        values = local.get(name)
        if values is not None and values.shape[1] == entry[0]:
            entry[1].append(values)
        else:
            entry[1].append(np.full((nv, entry[0]), np.nan))

    n = len(index)
    kind = 0
    cells = np.zeros((0, 3), dtype=np.int64)
    if mode == 0:
        kind, cells = 2, index.reshape(-1, 1)
    elif mode == 1:
        if n % 2:
            _fail("a LINES primitive has an odd number of indices")
        kind, cells = 1, index.reshape(-1, 2)
    elif mode in (2, 3):
        kind = 1
        cells = (
            np.stack([index[:-1], index[1:]], axis=1)
            if n > 1
            else np.zeros((0, 2), np.int64)
        )
        if mode == 2 and n > 2:
            cells = np.vstack([cells, [[index[-1], index[0]]]])
    elif mode == 4:
        if n % 3:
            _fail("a TRIANGLES primitive has a count that is not a multiple of 3")
        cells = index.reshape(-1, 3)
    elif mode == 5:
        if n >= 3:
            a, b, c = index[:-2], index[1:-1], index[2:]
            odd = (np.arange(n - 2) % 2) == 1
            first = np.where(odd, b, a)
            second = np.where(odd, a, b)
            keep = (a != b) & (b != c) & (a != c)
            cells = np.stack([first, second, c], axis=1)[keep]
    else:
        if n >= 3:
            a, b, c = np.full(n - 2, index[0]), index[1:-1], index[2:]
            keep = (a != b) & (b != c) & (a != c)
            cells = np.stack([a, b, c], axis=1)[keep]
    cells = cells.astype(np.int64) + base
    state.cells[kind].append(cells)
    first_id = state.cell_counts[kind]
    state.cell_counts[kind] += len(cells)
    if region and region not in ("unassigned", "mesh"):
        entry = state.regions.setdefault(region, ([], [], []))
        entry[kind].append(np.arange(first_id, first_id + len(cells), dtype=np.int64))


def _walk(ctx, node_index, parent, offset, visited, state):
    nodes = _list(ctx.root, "nodes")
    if node_index >= len(nodes):
        _fail("a node index is out of range")
    if node_index in visited:
        _fail("the node graph is not a forest (a node is reached twice)")
    visited.add(node_index)
    node = _object(nodes[node_index], f"nodes[{node_index}]")
    world = parent @ _node_matrix(node)
    if "mesh" in node:
        meshes = _list(ctx.root, "meshes")
        mi = _uint(node["mesh"], "node.mesh")
        if mi >= len(meshes):
            _fail("a mesh index is out of range")
        m = _object(meshes[mi], f"meshes[{mi}]")
        name = node.get("name")
        region = name if isinstance(name, str) else ""
        for prim in _list(m, "primitives"):
            prim = _object(prim, "primitive")
            if "targets" in prim and "morph" not in state.warned:
                state.warned.add("morph")
                warn("gltf: morph targets are skipped.")
            _primitive(ctx, prim, world, offset, region, state)
    children = node.get("children")
    if children is not None:
        if not isinstance(children, list):
            _fail("node.children is not an array")
        for c in children:
            _walk(ctx, _uint(c, "node.children"), world, offset, visited, state)


class _Ctx:
    def __init__(self, root, buffers):
        self.root = root
        self.buffers = buffers


def read(filename, points_only=False, arrays=None, weld=False):
    """Read a glTF 2.0 / GLB file (the Python twin of ``read_gltf``)."""
    root, binary, directory = _load(filename)
    for ext in _list(root, "extensionsRequired"):
        if ext in ("KHR_draco_mesh_compression", "EXT_meshopt_compression"):
            _fail(f"{ext} is not supported")
    ctx = _Ctx(root, _buffers(root, binary, directory))
    state = _State(not points_only and arrays != [], arrays)
    if arrays is not None:
        state.arrays = set(arrays)

    skipped = [k for k in ("animations", "skins") if _list(root, k)]
    if _list(root, "textures") or _list(root, "images"):
        skipped.append("textures")
    if skipped:
        warn("gltf: " + ", ".join(skipped) + " are skipped.")

    nodes = _list(root, "nodes")
    scenes = _list(root, "scenes")
    if scenes:
        si = _uint_or(root, "scene", 0, "document")
        if si >= len(scenes):
            _fail("the default scene index is out of range")
        roots = [
            _uint(n, "scene.nodes")
            for n in _list(_object(scenes[si], "scene"), "nodes")
        ]
    else:
        children = set()
        for n in nodes:
            for c in _list(n, "children"):
                children.add(_uint(c, "node.children"))
        roots = [i for i in range(len(nodes)) if i not in children]

    visited = set()
    zero = [0.0, 0.0, 0.0]
    for r in roots:
        if r >= len(nodes):
            _fail("a scene names a node that does not exist")
        node = _object(nodes[r], "node")
        offset = _writer_root(node) if len(roots) == 1 else None
        if offset is not None:
            visited.add(r)
            for c in _list(node, "children"):
                _walk(ctx, _uint(c, "node.children"), np.eye(4), offset, visited, state)
        else:
            _walk(ctx, r, np.eye(4), zero, visited, state)

    n = state.npoints
    points = np.vstack(state.points) if state.points else np.zeros((0, 3))
    fields = {}
    for name, (ncomp, chunks, dropped) in state.fields.items():
        if not dropped:
            fields[name] = np.vstack(chunks) if chunks else np.zeros((0, ncomp))
    cells = []
    for k in range(3):
        cells.append(
            np.vstack(state.cells[k])
            if state.cells[k]
            else np.zeros((0, (3, 2, 1)[k]), dtype=np.int64)
        )

    if weld and n:
        _, first, inverse = np.unique(
            points + 0.0, axis=0, return_index=True, return_inverse=True
        )
        inverse = np.asarray(inverse).reshape(-1)
        # Keep first-occurrence order, as the core does.
        order = np.argsort(first, kind="stable")
        rank = np.empty_like(order)
        rank[order] = np.arange(len(order))
        remap = rank[inverse]
        keep = first[order]
        points = points[keep]
        fields = {k: v[keep] for k, v in fields.items()}
        cells = [remap[c] if len(c) else c for c in cells]
        n = len(keep)

    blocks = []
    offsets = [0, 0, 0]
    total = 0
    for k, name in enumerate(("triangle", "line", "vertex")):
        offsets[k] = total
        if len(cells[k]):
            total += len(cells[k])
            blocks.append(CellBlock(name, cells[k].astype(np.int64)))
    point_data = {}
    for name, values in fields.items():
        point_data[name] = values[:, 0] if values.shape[1] == 1 else values
    regions = []
    for name in sorted(state.regions):
        ids = []
        for k in range(3):
            for chunk in state.regions[name][k]:
                ids.append(chunk + offsets[k])
        entries = np.concatenate(ids) if ids else np.zeros(0, dtype=np.int64)
        regions.append(Region(name, "cell", entries))
    return Mesh(
        points,
        blocks,
        point_data=point_data or None,
        regions=regions or None,
    )

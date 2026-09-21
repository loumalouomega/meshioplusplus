"""
I/O for the Point Cloud Library's PCD format (v0.7), cf.
<https://pointclouds.org/documentation/tutorials/pcd_file_format.html>.

A PCD file is a short text header (``VERSION``, ``FIELDS``, ``SIZE``, ``TYPE``,
``COUNT``, ``WIDTH``, ``HEIGHT``, ``VIEWPOINT``, ``POINTS``, ``DATA``) followed by
the points as ASCII rows, packed little-endian records (``binary``) or a
field-by-field (struct-of-arrays) block wrapped in an LZF stream
(``binary_compressed``). The mesh is the point set plus a single ``vertex``
block, exactly what ``subsample_points`` produces.
"""

import struct

import numpy as np

from .. import _provenance
from .._common import warn
from .._exceptions import ReadError, WriteError
from .._files import open_file
from .._mesh import Mesh
from . import _lzf

_NP_TYPES = {
    ("F", 4): "<f4",
    ("F", 8): "<f8",
    ("I", 1): "i1",
    ("I", 2): "<i2",
    ("I", 4): "<i4",
    ("I", 8): "<i8",
    ("U", 1): "u1",
    ("U", 2): "<u2",
    ("U", 4): "<u4",
    ("U", 8): "<u8",
}
_IDENTITY_VIEWPOINT = (0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0)
_DATA_MODES = ("ascii", "binary", "binary_compressed")
_POINT_DTYPES = {"float32": ("F", 4), "float64": ("F", 8)}


def read(filename, drop_invalid=False):
    with open_file(filename, "rb") as f:
        raw = f.read()
    if isinstance(raw, str):
        raw = raw.encode("latin-1")
    return _parse(bytes(raw), drop_invalid)


def _parse_header(raw):
    header = {}
    pos = 0
    while True:
        end = raw.find(b"\n", pos)
        if end < 0:
            raise ReadError("PCD: no DATA line found in the header")
        line = raw[pos:end].decode("latin-1").strip()
        pos = end + 1
        if not line or line.startswith("#"):
            continue
        key, _, rest = line.partition(" ")
        key = key.upper()
        header[key] = rest.split()
        if key == "DATA":
            return header, pos


def _parse(raw, drop_invalid):
    header, body = _parse_header(raw)
    if "FIELDS" not in header:
        raise ReadError("PCD: the header has no FIELDS line")
    names = header["FIELDS"]
    nfields = len(names)
    sizes = [int(s) for s in header.get("SIZE", [])]
    types = [t.upper() for t in header.get("TYPE", [])]
    counts = [int(c) for c in header.get("COUNT", ["1"] * nfields)]
    if not (len(sizes) == len(types) == len(counts) == nfields):
        raise ReadError("PCD: FIELDS, SIZE, TYPE and COUNT disagree on the field count")
    for t, s in zip(types, sizes):
        if (t, s) not in _NP_TYPES:
            raise ReadError(f"PCD: unsupported field type {t}{s}")

    try:
        width = int(header["WIDTH"][0]) if "WIDTH" in header else None
        height = int(header["HEIGHT"][0]) if "HEIGHT" in header else 1
        npoints = int(header["POINTS"][0]) if "POINTS" in header else None
    except (ValueError, IndexError):
        raise ReadError("PCD: malformed WIDTH/HEIGHT/POINTS in the header")
    if npoints is None:
        if width is None:
            raise ReadError("PCD: the header has neither POINTS nor WIDTH")
        npoints = width * height
    if width is None:
        width = npoints
    if width * height != npoints:
        warn(
            f"PCD: WIDTH*HEIGHT ({width}*{height}) != POINTS ({npoints}); using POINTS"
        )

    viewpoint = np.array(_IDENTITY_VIEWPOINT)
    if "VIEWPOINT" in header:
        try:
            vp = [float(v) for v in header["VIEWPOINT"]]
        except ValueError:
            raise ReadError("PCD: malformed VIEWPOINT")
        if len(vp) != 7:
            raise ReadError("PCD: VIEWPOINT needs 7 values (tx ty tz qw qx qy qz)")
        viewpoint = np.array(vp)

    mode = (header["DATA"][0] if header["DATA"] else "").lower()
    if mode not in _DATA_MODES:
        raise ReadError(f"PCD: unknown DATA mode '{mode}'")
    payload = raw[body:]

    dtypes = [np.dtype(_NP_TYPES[(t, s)]) for t, s in zip(types, sizes)]
    if mode == "ascii":
        columns = _read_ascii(payload, npoints, dtypes, counts)
    elif mode == "binary":
        columns = _read_binary(payload, npoints, dtypes, counts)
    else:
        columns = _read_compressed(payload, npoints, dtypes, counts)

    fields = dict(zip(names, zip(columns, types, sizes)))
    return _build_mesh(fields, names, npoints, width, height, viewpoint, drop_invalid)


def _read_ascii(payload, npoints, dtypes, counts):
    tokens = payload.split()
    per_row = sum(counts)
    if len(tokens) != npoints * per_row:
        raise ReadError(
            f"PCD: expected {npoints * per_row} values, found {len(tokens)}"
        )
    table = np.array(tokens).reshape(npoints, per_row) if npoints else None
    columns = []
    col = 0
    try:
        for dt, count in zip(dtypes, counts):
            if npoints == 0:
                columns.append(np.empty((0, count), dtype=dt.newbyteorder("=")))
                continue
            block = table[:, col : col + count]
            col += count
            if dt.kind == "f":
                values = block.astype(np.float64).astype(dt.newbyteorder("="))
            else:
                values = block.astype(np.int64 if dt.kind == "i" else np.uint64)
                values = values.astype(dt.newbyteorder("="))
            columns.append(values)
    except (ValueError, OverflowError):
        raise ReadError("PCD: non-numeric value in the ASCII data") from None
    return columns


def _read_binary(payload, npoints, dtypes, counts):
    dtype = np.dtype(
        [(f"f{i}", dt, (c,)) for i, (dt, c) in enumerate(zip(dtypes, counts))]
    )
    if len(payload) < npoints * dtype.itemsize:
        raise ReadError("PCD: binary data is shorter than the header declares")
    rows = np.frombuffer(payload, dtype=dtype, count=npoints)
    return [rows[f"f{i}"].astype(dt.newbyteorder("=")) for i, dt in enumerate(dtypes)]


def _read_compressed(payload, npoints, dtypes, counts):
    if len(payload) < 8:
        raise ReadError("PCD: binary_compressed data is missing its size prefix")
    comp, uncomp = struct.unpack("<II", payload[:8])
    if 8 + comp > len(payload):
        raise ReadError("PCD: binary_compressed data is shorter than its size prefix")
    expected = npoints * sum(dt.itemsize * c for dt, c in zip(dtypes, counts))
    if uncomp != expected:
        raise ReadError(
            f"PCD: binary_compressed declares {uncomp} bytes, the header implies {expected}"
        )
    raw = _lzf.decompress(payload[8 : 8 + comp], uncomp)
    columns = []
    offset = 0
    for dt, count in zip(dtypes, counts):
        block = np.frombuffer(raw, dtype=dt, count=npoints * count, offset=offset)
        columns.append(block.reshape(npoints, count).astype(dt.newbyteorder("=")))
        offset += npoints * count * dt.itemsize
    return columns


def _unpack_colour(column, type_, alpha):
    packed = np.ascontiguousarray(column[:, 0])
    if type_ == "F":
        packed = packed.astype(np.float32).view(np.uint32)  # bit-cast, never by value
    else:
        packed = packed.astype(np.uint32)
    channels = [(packed >> 16) & 255, (packed >> 8) & 255, packed & 255]
    if alpha:
        channels.append((packed >> 24) & 255)
    return np.stack(channels, axis=1).astype(np.uint8)


def _build_mesh(fields, names, npoints, width, height, viewpoint, drop_invalid):
    for axis in "xyz":
        if axis not in fields or fields[axis][0].shape[1] != 1:
            raise ReadError(f"PCD: the cloud has no scalar '{axis}' field")
    xyz = [fields[a][0][:, 0] for a in "xyz"]
    single = all(fields[a][0].dtype == np.float32 for a in "xyz")
    points = np.stack(xyz, axis=1).astype(np.float32 if single else np.float64)

    point_data = {}
    consumed = {"x", "y", "z"}
    normal_axes = ["normal_x", "normal_y", "normal_z"]
    if all(a in fields and fields[a][0].shape[1] == 1 for a in normal_axes):
        parts = [fields[a][0][:, 0] for a in normal_axes]
        single_n = all(p.dtype == np.float32 for p in parts)
        point_data["normals"] = np.stack(parts, axis=1).astype(
            np.float32 if single_n else np.float64
        )
        consumed.update(normal_axes)
    for name in names:
        if name in consumed or name == "_" or name in point_data:
            continue
        column, type_, size = fields[name]
        if name in ("rgb", "rgba") and column.shape[1] == 1 and size == 4:
            point_data[name] = _unpack_colour(column, type_, name == "rgba")
        elif column.shape[1] == 1:
            point_data[name] = np.ascontiguousarray(column[:, 0])
        else:
            point_data[name] = np.ascontiguousarray(column)

    field_data = {}
    organised = height > 1 and width * height == npoints
    if drop_invalid:
        keep = np.isfinite(points).all(axis=1)
        if not keep.all():
            points = points[keep]
            point_data = {k: v[keep] for k, v in point_data.items()}
            organised = False
    if organised:
        field_data["pcd:width"] = np.array([width], dtype=np.int64)
        field_data["pcd:height"] = np.array([height], dtype=np.int64)
    if tuple(viewpoint) != _IDENTITY_VIEWPOINT:
        field_data["pcd:viewpoint"] = viewpoint

    n = len(points)
    cells = [("vertex", np.arange(n, dtype=np.int64).reshape(-1, 1))]
    return Mesh(points, cells, point_data=point_data, field_data=field_data)


# ---------------------------------------------------------------------------
# writer


def _unique(name, used):
    name = "".join("_" if c.isspace() else c for c in name) or "field"
    candidate, k = name, 2
    while candidate in used:
        candidate = f"{name}_{k}"
        k += 1
    used.add(candidate)
    return candidate


def _pack_colour(values, alpha):
    channels = np.clip(np.rint(values.astype(np.float64)), 0, 255).astype(np.uint32)
    packed = (channels[:, 0] << 16) | (channels[:, 1] << 8) | channels[:, 2]
    if alpha:
        packed |= channels[:, 3] << 24
    return packed.astype(np.uint32)


def _typed(array):
    """(type letter, size, array) for a generic point-data array."""
    kind = array.dtype.kind
    if kind == "b":
        return "U", 1, array.astype(np.uint8)
    if kind == "f":
        size = 8 if array.dtype.itemsize > 4 else 4
        return "F", size, array.astype(np.float64 if size == 8 else np.float32)
    if kind in "iu":
        return kind.upper(), array.dtype.itemsize, array
    return None


def _collect_fields(points, point_data, point_type):
    """The ordered field list: (name, type, size, (n, count) array)."""
    n = len(points)
    ptype = _POINT_DTYPES[point_type]
    np_point = np.float64 if ptype[1] == 8 else np.float32
    used = {"x", "y", "z"}
    fields = [
        (a, "F", ptype[1], points[:, i : i + 1].astype(np_point))
        for i, a in enumerate("xyz")
    ]
    dropped = []
    for name in sorted(point_data):
        array = np.asarray(point_data[name])
        if name == "normals" and array.shape == (n, 3):
            for axis, i in zip("xyz", range(3)):
                field = _unique(f"normal_{axis}", used)
                fields.append(
                    (field, "F", ptype[1], array[:, i : i + 1].astype(np_point))
                )
            continue
        if name in ("rgb", "rgba") and array.dtype.kind in "iuf":
            width = 4 if name == "rgba" else 3
            if array.shape == (n, width):
                packed = _pack_colour(array, name == "rgba")
                if name == "rgb":
                    packed = packed.view(np.float32)
                    fields.append((_unique(name, used), "F", 4, packed.reshape(n, 1)))
                else:
                    fields.append((_unique(name, used), "U", 4, packed.reshape(n, 1)))
                continue
        if name in ("curvature", "intensity") and array.dtype.kind == "f":
            array = array.astype(np_point)
        typed = _typed(array)
        if typed is None:
            dropped.append(name)
            continue
        type_, size, array = typed
        count = int(np.prod(array.shape[1:])) if array.ndim > 1 else 1
        fields.append((_unique(name, used), type_, size, array.reshape(n, count)))
    return fields, dropped


def _format_ascii_value(value, type_, size):
    if type_ == "F":
        if value != value:
            return "nan"
        return format(float(value), ".9g" if size == 4 else ".17g")
    return str(int(value))


def _ascii_body(fields, n):
    columns = []
    for _, type_, size, array in fields:
        for j in range(array.shape[1]):
            columns.append([_format_ascii_value(v, type_, size) for v in array[:, j]])
    lines = [" ".join(col[i] for col in columns) for i in range(n)]
    return ("\n".join(lines) + "\n").encode() if n else b""


def write(filename, mesh, binary=True, data=None, point_dtype="float32"):
    if data is None:
        data = "binary" if binary else "ascii"
    if data not in _DATA_MODES:
        raise WriteError(f"PCD: data must be one of {_DATA_MODES}, got '{data}'")
    if point_dtype not in _POINT_DTYPES:
        raise WriteError("PCD: point_dtype must be 'float32' or 'float64'")

    points = np.asarray(mesh.points)
    if points.shape[1] < 3:
        warn("PCD requires 3D points; padding with zeros.")
        _provenance.note("point-padding", "points padded with zero coordinates to 3D")
        pad = np.zeros((len(points), 3 - points.shape[1]), dtype=points.dtype)
        points = np.concatenate([points, pad], axis=1)
    if point_dtype == "float32" and points.dtype != np.float32:
        narrowed = points.astype(np.float32)
        if not np.array_equal(
            narrowed.astype(np.float64), points.astype(np.float64), equal_nan=True
        ):
            _provenance.note(
                "dtype",
                "float64 coordinates written as float32 (point_dtype='float32')",
            )

    skipped = [c.type for c in mesh.cells if c.type != "vertex"]
    if skipped:
        string = ", ".join(skipped)
        warn(f"PCD holds points only. Skipping {string} cells.")
        _provenance.note(
            "cells-dropped", f"cell block(s) of type {string} have no PCD equivalent"
        )
    if mesh.cell_data:
        _provenance.note("data-dropped", "cell data has no PCD equivalent")
    extra = sorted(k for k in mesh.field_data if not k.startswith("pcd:"))
    if extra:
        _provenance.note(
            "data-dropped", "field data has no PCD equivalent: " + ", ".join(extra)
        )

    n = len(points)
    fields, dropped = _collect_fields(points[:, :3], mesh.point_data, point_dtype)
    if dropped:
        string = ", ".join(dropped)
        warn(f"PCD cannot store point data {string}; skipping.")
        _provenance.note(
            "data-dropped", f"point data not representable in PCD: {string}"
        )

    width, height = n, 1
    fd = mesh.field_data
    if "pcd:width" in fd and "pcd:height" in fd:
        w = int(np.asarray(fd["pcd:width"]).ravel()[0])
        h = int(np.asarray(fd["pcd:height"]).ravel()[0])
        if w * h == n and n:
            width, height = w, h
    viewpoint = _IDENTITY_VIEWPOINT
    if "pcd:viewpoint" in fd:
        vp = np.asarray(fd["pcd:viewpoint"], dtype=np.float64).ravel()
        if vp.size == 7:
            viewpoint = tuple(vp)

    head = [
        "# .PCD v0.7 - Point Cloud Data file format\n",
        _provenance.render_lines(_provenance.SlotTier.BLOCK, "# "),
        "VERSION 0.7\n",
        "FIELDS " + " ".join(f[0] for f in fields) + "\n",
        "SIZE " + " ".join(str(f[2]) for f in fields) + "\n",
        "TYPE " + " ".join(f[1] for f in fields) + "\n",
        "COUNT " + " ".join(str(f[3].shape[1]) for f in fields) + "\n",
        f"WIDTH {width}\n",
        f"HEIGHT {height}\n",
        "VIEWPOINT " + " ".join(format(float(v), ".17g") for v in viewpoint) + "\n",
        f"POINTS {n}\n",
        f"DATA {data}\n",
    ]
    if data == "ascii":
        body = _ascii_body(fields, n)
    elif data == "binary":
        dtype = np.dtype(
            [
                (f"f{i}", _NP_TYPES[(t, s)], (a.shape[1],))
                for i, (_, t, s, a) in enumerate(fields)
            ]
        )
        rows = np.empty(n, dtype=dtype)
        for i, (_, _, _, a) in enumerate(fields):
            rows[f"f{i}"] = a
        body = rows.tobytes()
    else:
        blocks = b"".join(
            np.ascontiguousarray(a, dtype=_NP_TYPES[(t, s)]).tobytes()
            for _, t, s, a in fields
        )
        packed = _lzf.compress(blocks)
        body = struct.pack("<II", len(packed), len(blocks)) + packed

    with open_file(filename, "wb") as fh:
        fh.write("".join(head).encode())
        fh.write(body)

"""ASPRS LAS point clouds: the Python reference of ``formats/las.cpp`` (roadmap
2.4.6.1). See that header for the fields; both engines write the same bytes."""

from __future__ import annotations

import struct

import numpy as np

from .._common import warn
from .._exceptions import ReadError, WriteError
from .._mesh import Mesh

RECORD = (20, 28, 26, 34, 57, 63, 30, 36, 38, 59, 67)


def _fields(fmt):
    """``[(name, dtype, extract(records))]``; ``records`` is ``(n, length)`` uint8."""

    def u8(off):
        return lambda r: r[:, off]

    def bits(off, shift, width):
        return lambda r: (r[:, off] >> shift) & ((1 << width) - 1)

    def typed(off, dt):
        size = np.dtype(dt).itemsize
        return (
            lambda r: np.ascontiguousarray(r[:, off : off + size]).view(dt).reshape(-1)
        )

    f = [("intensity", np.uint16, typed(12, "<u2"))]
    if fmt < 6:
        f += [
            ("return_number", np.uint8, bits(14, 0, 3)),
            ("number_of_returns", np.uint8, bits(14, 3, 3)),
            ("scan_direction_flag", np.uint8, bits(14, 6, 1)),
            ("edge_of_flight_line", np.uint8, bits(14, 7, 1)),
            ("classification", np.uint8, bits(15, 0, 5)),
            ("synthetic", np.uint8, bits(15, 5, 1)),
            ("key_point", np.uint8, bits(15, 6, 1)),
            ("withheld", np.uint8, bits(15, 7, 1)),
            ("scan_angle_rank", np.int8, typed(16, "i1")),
            ("user_data", np.uint8, u8(17)),
            ("point_source_id", np.uint16, typed(18, "<u2")),
        ]
    else:
        f += [
            ("return_number", np.uint8, bits(14, 0, 4)),
            ("number_of_returns", np.uint8, bits(14, 4, 4)),
            ("synthetic", np.uint8, bits(15, 0, 1)),
            ("key_point", np.uint8, bits(15, 1, 1)),
            ("withheld", np.uint8, bits(15, 2, 1)),
            ("overlap", np.uint8, bits(15, 3, 1)),
            ("scanner_channel", np.uint8, bits(15, 4, 2)),
            ("scan_direction_flag", np.uint8, bits(15, 6, 1)),
            ("edge_of_flight_line", np.uint8, bits(15, 7, 1)),
            ("classification", np.uint8, u8(16)),
            ("user_data", np.uint8, u8(17)),
            ("scan_angle", np.int16, typed(18, "<i2")),
            ("point_source_id", np.uint16, typed(20, "<u2")),
        ]
    gps, rgb, nir = {
        1: (20, 0, 0),
        4: (20, 0, 0),
        2: (0, 20, 0),
        3: (20, 28, 0),
        5: (20, 28, 0),
        6: (22, 0, 0),
        9: (22, 0, 0),
        7: (22, 30, 0),
        8: (22, 30, 36),
        10: (22, 30, 36),
    }.get(fmt, (0, 0, 0))
    if gps:
        f.append(("gps_time", np.float64, typed(gps, "<f8")))
    if rgb:
        f += [
            ("red", np.uint16, typed(rgb, "<u2")),
            ("green", np.uint16, typed(rgb + 2, "<u2")),
            ("blue", np.uint16, typed(rgb + 4, "<u2")),
        ]
    if nir:
        f.append(("nir", np.uint16, typed(nir, "<u2")))
    return f


def read(filename, every=1, arrays=None, points_only=False):
    try:
        with open(filename, "rb") as fh:
            d = fh.read()
    except OSError:
        raise ReadError(f"meshio++: las: could not open file: {filename}") from None
    if len(d) < 227 or d[:4] != b"LASF":
        raise ReadError(
            f"meshio++: las: '{filename}' is not a LAS file (no LASF signature)"
        )
    major, minor = d[24], d[25]
    if major != 1 or minor > 4:
        raise ReadError(
            f"meshio++: las: version {major}.{minor} is not read (1.0 to 1.4)"
        )
    header, offset, nvlr = struct.unpack_from("<HII", d, 94)
    if d[104] & 0xC0:
        raise ReadError(
            f"meshio++: las: '{filename}' is LAZ-compressed; meshio++ reads uncompressed "
            "LAS (LAZ needs LASzip)"
        )
    fmt = d[104] & 0x3F
    if fmt > 10:
        raise ReadError(f"meshio++: las: point format {fmt} is not read")
    (record,) = struct.unpack_from("<H", d, 105)
    if record < RECORD[fmt]:
        raise ReadError(
            f"meshio++: las: point record length {record} is shorter than format {fmt} needs"
        )
    (count,) = struct.unpack_from("<I", d, 107)
    if minor >= 4 and header >= 375 and count == 0:
        (count,) = struct.unpack_from("<Q", d, 247)
    scale = np.array(struct.unpack_from("<3d", d, 131))
    origin = np.array(struct.unpack_from("<3d", d, 155))
    if header > len(d) or offset > len(d) or offset < header:
        raise ReadError(
            "meshio++: las: the header size or point offset is out of range"
        )
    if count > (len(d) - offset) // record:
        raise ReadError(
            f"meshio++: las: the header promises {count} points but the file holds "
            f"{(len(d) - offset) // record}"
        )
    at = header
    for _ in range(nvlr):
        if at + 54 > offset:
            raise ReadError("meshio++: las: a VLR runs past the point data")
        at += 54 + struct.unpack_from("<H", d, at + 20)[0]
    if at > offset:
        raise ReadError("meshio++: las: the VLRs run past the point data")
    vlrs = d[header:at]
    every = max(int(every), 1)
    records = np.frombuffer(d, dtype=np.uint8, count=count * record, offset=offset)
    records = records.reshape(count, record)[::every]
    n = len(records)
    xyz = np.ascontiguousarray(records[:, :12]).view("<i4").reshape(n, 3)
    points = xyz * scale + origin
    point_data = {}
    if not points_only:
        for name, dtype, get in _fields(fmt):
            if arrays is not None and name not in arrays:
                continue
            point_data[name] = np.asarray(get(records)).astype(dtype)
    field_data = {
        "las:scale": scale,
        "las:offset": origin,
        "las:version": np.array([major, minor], dtype=np.int64),
        "las:point_format": np.array([fmt], dtype=np.int64),
    }
    if vlrs:
        field_data["las:vlrs"] = np.frombuffer(vlrs, dtype=np.uint8).copy()
    if minor >= 4 and header >= 375 and struct.unpack_from("<I", d, 243)[0] > 0:
        warn(f"las: '{filename}' has extended VLRs, which are not kept")
    cells = [("vertex", np.arange(n, dtype=np.int64).reshape(-1, 1))] if n else []
    return Mesh(points, cells, point_data=point_data, field_data=field_data)


def write(filename, mesh):
    pts = np.asarray(mesh.points, dtype=np.float64)
    n = len(pts)
    full = np.zeros((n, 3))
    full[:, : min(3, pts.shape[1])] = pts[:, :3]

    def has(name):
        return name in mesh.point_data and np.asarray(mesh.point_data[name]).size == n

    rgb = has("red") and has("green") and has("blue")
    fmt = (8 if has("nir") else 7) if rgb else 6
    record = RECORD[fmt]
    lo = full.min(axis=0) if n else np.zeros(3)
    hi = full.max(axis=0) if n else np.zeros(3)
    scale = np.array([0.001, 0.001, 0.001])
    origin = np.floor(lo)
    s = mesh.field_data.get("las:scale")
    if s is not None and np.asarray(s).size == 3:
        scale = np.asarray(s, dtype=np.float64).reshape(3)
    o = mesh.field_data.get("las:offset")
    if o is not None and np.asarray(o).size == 3:
        origin = np.asarray(o, dtype=np.float64).reshape(3)
    if not np.all(np.isfinite(scale)) or np.any(scale <= 0):
        raise WriteError("meshio++: las: las:scale must be three positive numbers")
    vlrs = b""
    v = mesh.field_data.get("las:vlrs")
    if v is not None:
        vlrs = np.asarray(v, dtype=np.uint8).tobytes()
    nvlr, at = 0, 0
    while at + 54 <= len(vlrs):
        at += 54 + struct.unpack_from("<H", vlrs, at + 20)[0]
        nvlr += 1

    def field(name, default):
        if has(name):
            return np.asarray(mesh.point_data[name], dtype=np.float64).reshape(-1)
        return np.full(n, float(default))

    q = np.round((full - origin) / scale)
    if np.any(q < np.iinfo(np.int32).min) or np.any(q > np.iinfo(np.int32).max):
        raise WriteError(
            "meshio++: las: a coordinate does not fit the int32 grid of las:scale and "
            "las:offset"
        )
    ret = field("return_number", 1).astype(np.int64).astype(np.uint64) & 15
    nret = field("number_of_returns", 1).astype(np.int64).astype(np.uint64) & 15

    def flag(name, shift, mask):
        return (field(name, 0).astype(np.int64).astype(np.uint64) & mask) << shift

    flags = (
        flag("synthetic", 0, 1)
        | flag("key_point", 1, 1)
        | flag("withheld", 2, 1)
        | flag("overlap", 3, 1)
        | flag("scanner_channel", 4, 3)
        | flag("scan_direction_flag", 6, 1)
        | flag("edge_of_flight_line", 7, 1)
    )
    angle = (
        field("scan_angle", 0)
        if has("scan_angle")
        else np.round(field("scan_angle_rank", 0) / 0.006)
    )
    dtype = [
        ("xyz", "<i4", 3),
        ("intensity", "<u2"),
        ("returns", "u1"),
        ("flags", "u1"),
        ("classification", "u1"),
        ("user_data", "u1"),
        ("scan_angle", "<i2"),
        ("point_source_id", "<u2"),
        ("gps_time", "<f8"),
    ]
    if fmt >= 7:
        dtype += [("red", "<u2"), ("green", "<u2"), ("blue", "<u2")]
    if fmt == 8:
        dtype += [("nir", "<u2")]
    rec = np.zeros(n, dtype=dtype)
    rec["xyz"] = q.astype(np.int32)

    def as_int(values, kind):
        return values.astype(np.int64).astype(kind)

    rec["intensity"] = as_int(field("intensity", 0), np.uint16)
    rec["returns"] = (ret | (nret << 4)).astype(np.uint8)
    rec["flags"] = flags.astype(np.uint8)
    rec["classification"] = as_int(field("classification", 0), np.uint64) & 255
    rec["user_data"] = as_int(field("user_data", 0), np.uint64) & 255
    rec["scan_angle"] = as_int(angle, np.int16)
    rec["point_source_id"] = as_int(field("point_source_id", 0), np.uint16)
    rec["gps_time"] = field("gps_time", 0)
    for name in ("red", "green", "blue", "nir"):
        if name in rec.dtype.names:
            rec[name] = as_int(field(name, 0), np.uint16)
    by_return = np.bincount(ret[ret >= 1].astype(np.int64) - 1, minlength=15)[:15]
    ident = b"meshio++".ljust(32, b"\0")
    head = (
        b"LASF"
        + struct.pack("<HH", 0, 16)
        + b"\0" * 16
        + bytes([1, 4])
        + ident
        + ident
        + struct.pack("<HHHII", 0, 0, 375, 375 + len(vlrs), nvlr)
        + struct.pack("<BH", fmt, record)
        + struct.pack("<6I", 0, 0, 0, 0, 0, 0)
        + struct.pack("<3d", *scale)
        + struct.pack("<3d", *origin)
        + struct.pack("<6d", hi[0], lo[0], hi[1], lo[1], hi[2], lo[2])
        + struct.pack("<QQIQ", 0, 0, 0, n)
        + struct.pack("<15Q", *[int(c) for c in by_return])
    )
    with open(filename, "wb") as f:
        f.write(head + vlrs + rec.tobytes())

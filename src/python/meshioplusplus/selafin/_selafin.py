"""TELEMAC Selafin results: the Python reference of ``formats/selafin.cpp`` (roadmap
2.4.3.4). See that header for the layout."""

from __future__ import annotations

import struct

import numpy as np

from .._exceptions import ReadError, WriteError
from .._fortran_records import fortran_record, fortran_records
from .._mesh import Mesh

UNIT_PREFIX = "selafin:unit:"
_LAYOUT = (4, ">")


def _bytes(text):
    return np.frombuffer(text.encode("latin-1"), dtype=np.uint8).copy()


def _trim(raw):
    return raw.decode("latin-1").rstrip(" \x00").lstrip(" ")


class _Records:
    def __init__(self, path):
        try:
            with open(path, "rb") as f:
                self.data = f.read()
        except OSError:
            raise ReadError(f"meshio++: selafin: could not open file: {path}") from None
        self.spans = fortran_records(self.data, _LAYOUT, "selafin", lenient=True)
        self.at = 0

    def next(self, what):
        if self.at >= len(self.spans):
            raise ReadError(f"meshio++: selafin: the file ends before its {what}")
        o, n = self.spans[self.at]
        self.at += 1
        return self.data[o : o + n]


def _header(rec):
    title = rec.next("title")
    if len(title) != 80:
        raise ReadError(
            "meshio++: selafin: the title record is not 80 characters (not a Selafin "
            "file, or not big-endian)"
        )
    h = {"title": _trim(title)}
    nbv = rec.next("variable counts")
    if len(nbv) != 8:
        raise ReadError("meshio++: selafin: the NBV1 NBV2 record is malformed")
    nvar = sum(struct.unpack(">2i", nbv))
    if nvar < 0 or nvar > len(rec.spans):
        raise ReadError(f"meshio++: selafin: variable count {nvar} is out of range")
    names, units = [], []
    for _ in range(nvar):
        r = rec.next("variable name")
        if len(r) != 32:
            raise ReadError("meshio++: selafin: a variable record is not 32 characters")
        names.append(_trim(r[:16]))
        units.append(_trim(r[16:]))
    iparam = rec.next("IPARAM")
    if len(iparam) != 40:
        raise ReadError("meshio++: selafin: IPARAM is not ten integers")
    iparam = struct.unpack(">10i", iparam)
    if iparam[9] == 1:
        rec.next("date")
    counts = rec.next("NELEM NPOIN NDP")
    if len(counts) != 16:
        raise ReadError("meshio++: selafin: the NELEM NPOIN NDP 1 record is malformed")
    nelem, npoin, ndp, _ = struct.unpack(">4i", counts)
    limit = len(rec.data)
    for what, n in (("element", nelem), ("point", npoin)):
        if n < 0 or n > limit:
            raise ReadError(
                f"meshio++: selafin: {what} count {n} exceeds what the file holds ({limit})"
            )
    if ndp < 0 or ndp > 64:
        raise ReadError(
            f"meshio++: selafin: node-per-element count {ndp} is out of range"
        )
    ikle = rec.next("IKLE")
    if len(ikle) != nelem * ndp * 4:
        raise ReadError("meshio++: selafin: IKLE does not hold NELEM * NDP integers")
    ikle = np.frombuffer(ikle, dtype=">i4").astype(np.int64)
    if len(ikle) and (ikle.min() < 1 or ikle.max() > npoin):
        raise ReadError("meshio++: selafin: IKLE names a point out of range")
    ipobo = rec.next("IPOBO")
    if len(ipobo) != npoin * 4:
        raise ReadError("meshio++: selafin: IPOBO does not hold NPOIN integers")
    x = rec.next("X")
    if npoin > 0 and len(x) == npoin * 8:
        w = 8
    elif len(x) == npoin * 4:
        w = 4
    else:
        raise ReadError("meshio++: selafin: X does not hold NPOIN reals")
    y = rec.next("Y")
    if len(y) != len(x):
        raise ReadError("meshio++: selafin: Y does not match X")
    dt = ">f4" if w == 4 else ">f8"
    h.update(
        names=names,
        units=units,
        nplan=iparam[6],
        nelem=nelem,
        npoin=npoin,
        ndp=ndp,
        w=w,
        ikle=(ikle - 1).reshape(nelem, ndp),
        ipobo=np.frombuffer(ipobo, dtype=">i4").astype(np.int32),
        x=np.frombuffer(x, dtype=dt).astype(np.float64),
        y=np.frombuffer(y, dtype=dt).astype(np.float64),
        first=rec.at,
        dtype=dt,
    )
    times = []
    r = rec.at
    while r < len(rec.spans) and rec.spans[r][1] == w:
        if not all(
            r + 1 + v < len(rec.spans) and rec.spans[r + 1 + v][1] == npoin * w
            for v in range(nvar)
        ):
            break
        o = rec.spans[r][0]
        times.append(float(np.frombuffer(rec.data[o : o + w], dtype=dt)[0]))
        r += 1 + nvar
    h["times"] = times
    return h


def _resolve(times, time_step):
    n = len(times)
    if n == 0:
        if time_step in (0, -1):
            return 0
        raise ReadError(
            f"meshio++: time step {time_step} requested, but this file carries no time steps"
        )
    k = time_step + n if time_step < 0 else time_step
    if not 0 <= k < n:
        raise ReadError(
            f"meshio++: time step {time_step} is out of range: this file has {n} "
            f"step{'' if n == 1 else 's'}"
        )
    return k


def time_values(filename):
    return _header(_Records(filename))["times"]


def read(filename, time_step=0):
    rec = _Records(filename)
    h = _header(rec)
    step = _resolve(h["times"], time_step)
    points = np.zeros((h["npoin"], 3))
    points[:, 0], points[:, 1] = h["x"], h["y"]
    ndp, nplan = h["ndp"], h["nplan"]
    kind = {3: "triangle", 4: "tetra" if nplan > 1 else "quad", 6: "wedge"}.get(ndp)
    if kind is None:
        raise ReadError(
            f"meshio++: selafin: elements of {ndp} nodes are not a TELEMAC element"
        )
    values = []
    if h["times"]:
        base = h["first"] + step * (1 + len(h["names"])) + 1
        for v in range(len(h["names"])):
            o, n = rec.spans[base + v]
            values.append(
                np.frombuffer(rec.data[o : o + n], dtype=h["dtype"]).astype(np.float64)
            )
        if nplan > 1:
            for v, name in enumerate(h["names"]):
                if name.startswith("ELEVATION Z") or name.startswith("COTE Z"):
                    points[:, 2] = values[v]
                    break
    field_data = {"selafin:title": _bytes(h["title"])}
    if h["times"]:
        field_data["meshio:time"] = np.array([h["times"][step]])
    point_data = {}
    for name, unit, value in zip(h["names"], h["units"], values):
        point_data[name] = value
        field_data[UNIT_PREFIX + name] = _bytes(unit)
    point_data["selafin:ipobo"] = h["ipobo"]
    return Mesh(
        points,
        [(kind, h["ikle"])] if h["nelem"] else [],
        point_data=point_data,
        field_data=field_data,
    )


def _padded(text, width):
    return text[:width].ljust(width).encode("latin-1")


def write(filename, mesh):
    if len(mesh.cells) != 1 or mesh.cells[0].type not in ("triangle", "quad"):
        raise WriteError(
            "meshio++: selafin: the writer takes one triangle or quad block (a 2-D "
            "TELEMAC mesh)"
        )
    conn = np.asarray(mesh.cells[0].data, dtype=np.int64)
    npoin = len(mesh.points)
    nelem, ndp = conn.shape
    names = [
        n
        for n in sorted(mesh.point_data)
        if n != "selafin:ipobo"
        and np.asarray(mesh.point_data[n]).size == npoin
        and np.asarray(mesh.point_data[n]).dtype in (np.float64, np.float32)
    ]
    title = "meshio++"
    if "selafin:title" in mesh.field_data:
        title = np.asarray(mesh.field_data["selafin:title"], np.uint8).tobytes()
        title = title.decode("latin-1")
    out = [fortran_record(_padded(title, 72) + b"SERAFIND", _LAYOUT)]
    out.append(fortran_record(struct.pack(">2i", len(names), 0), _LAYOUT))
    for name in names:
        unit = mesh.field_data.get(UNIT_PREFIX + name)
        unit = (
            ""
            if unit is None
            else np.asarray(unit, np.uint8).tobytes().decode("latin-1")
        )
        out.append(fortran_record(_padded(name, 16) + _padded(unit, 16), _LAYOUT))
    out.append(
        fortran_record(struct.pack(">10i", 1, 0, 0, 0, 0, 0, 0, 0, 0, 0), _LAYOUT)
    )
    out.append(fortran_record(struct.pack(">4i", nelem, npoin, ndp, 1), _LAYOUT))
    out.append(fortran_record((conn + 1).astype(">i4").tobytes(), _LAYOUT))
    ipobo = mesh.point_data.get("selafin:ipobo")
    if ipobo is None or np.asarray(ipobo).size != npoin:
        ipobo = np.zeros(npoin)
    out.append(fortran_record(np.asarray(ipobo).astype(">i4").tobytes(), _LAYOUT))
    pts = np.asarray(mesh.points, dtype=np.float64)
    for c in range(2):
        col = pts[:, c] if c < pts.shape[1] else np.zeros(npoin)
        out.append(fortran_record(col.astype(">f8").tobytes(), _LAYOUT))
    time = mesh.field_data.get("meshio:time")
    time = (
        float(np.asarray(time).reshape(-1)[0])
        if time is not None and np.size(time)
        else 0.0
    )
    out.append(fortran_record(struct.pack(">d", time), _LAYOUT))
    for name in names:
        out.append(
            fortran_record(
                np.asarray(mesh.point_data[name]).astype(">f8").tobytes(), _LAYOUT
            )
        )
    with open(filename, "wb") as f:
        f.write(b"".join(out))

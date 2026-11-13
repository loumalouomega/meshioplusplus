"""FreeSurfer surfaces, curv overlays and annot parcellations: the Python reference of
``formats/freesurfer.cpp`` (roadmap 2.4.4)."""

from __future__ import annotations

import os
import struct

import numpy as np

from .._exceptions import ReadError, WriteError
from .._mesh import Mesh
from .._regions import Region


class _Cursor:
    def __init__(self, path, what):
        try:
            with open(path, "rb") as f:
                self.data = f.read()
        except OSError:
            raise ReadError(
                f"meshio++: freesurfer: could not open file: {path}"
            ) from None
        self.at = 0
        self.what = f"{what} '{path}'"

    def need(self, n):
        if len(self.data) - self.at < n:
            raise ReadError(f"meshio++: freesurfer: {self.what} is truncated")

    def i32(self, n=None):
        count = 1 if n is None else n
        self.need(4 * count)
        out = np.frombuffer(self.data, dtype=">i4", count=count, offset=self.at)
        self.at += 4 * count
        return int(out[0]) if n is None else out.astype(np.int64)

    def f32(self, n):
        self.need(4 * n)
        out = np.frombuffer(self.data, dtype=">f4", count=n, offset=self.at)
        self.at += 4 * n
        return out

    def magic(self):
        self.need(3)
        b = self.data[self.at : self.at + 3]
        self.at += 3
        return (b[0] << 16) | (b[1] << 8) | b[2]

    def string(self):
        n = self.i32()
        if n < 0:
            raise ReadError(
                f"meshio++: freesurfer: a negative string length in {self.what}"
            )
        self.need(n)
        s = self.data[self.at : self.at + n]
        self.at += n
        return s.split(b"\0", 1)[0].decode("latin-1")


def _overlay_name(path):
    base = os.path.basename(str(path))
    if base.startswith(("lh.", "rh.")):
        base = base[3:]
    return base or "curv"


def _add_curv(mesh, path):
    c = _Cursor(path, "curv file")
    if c.magic() != 0xFFFFFF:
        raise ReadError(
            f"meshio++: freesurfer: '{path}' is not a new-format curv file (old-format "
            "curv is not read)"
        )
    nv = c.i32()
    c.i32()
    per = c.i32()
    if nv != len(mesh.points):
        raise ReadError(
            f"meshio++: freesurfer: '{path}' has {nv} values but the surface has "
            f"{len(mesh.points)} vertices"
        )
    if per != 1:
        raise ReadError(
            f"meshio++: freesurfer: '{path}' has {per} values per vertex; meshio++ reads one"
        )
    mesh.point_data[_overlay_name(path)] = c.f32(nv).astype(np.float32)


def _add_annot(mesh, path):
    c = _Cursor(path, "annot file")
    nv = c.i32()
    if nv != len(mesh.points):
        raise ReadError(
            f"meshio++: freesurfer: '{path}' labels {nv} vertices but the surface has "
            f"{len(mesh.points)}"
        )
    pairs = c.i32(2 * nv).reshape(-1, 2) if nv else np.zeros((0, 2), np.int64)
    if len(pairs) and (pairs[:, 0].min() < 0 or pairs[:, 0].max() >= nv):
        raise ReadError(f"meshio++: freesurfer: '{path}' names a vertex out of range")
    colour = np.zeros(nv, dtype=np.int64)
    colour[pairs[:, 0]] = pairs[:, 1]
    entries = []
    if c.at < len(c.data) and c.i32() == 1:
        n = c.i32()

        def rgb():
            r, g, b, _ = (c.i32() for _ in range(4))
            return r + g * 256 + b * 65536

        if n > 0:
            c.string()
            if n > len(c.data):
                raise ReadError(
                    "meshio++: freesurfer: colour-table entry count out of range"
                )
            for k in range(n):
                name = c.string()
                entries.append((k, name, rgb()))
        else:
            if -n != 2:
                raise ReadError(
                    f"meshio++: freesurfer: colour table version {-n} is not read"
                )
            c.i32()
            c.string()
            count = c.i32()
            if count < 0 or count > len(c.data):
                raise ReadError(
                    "meshio++: freesurfer: colour-table entry count out of range"
                )
            for _ in range(count):
                index = c.i32()
                name = c.string()
                entries.append((index, name, rgb()))
    by_colour = {}
    for index, _, col in entries:
        by_colour.setdefault(col, index)
    labels = np.array([by_colour.get(int(v), -1) for v in colour], dtype=np.int32)
    for index, name, _ in entries:
        mesh.regions.append(
            Region(name, "point", np.flatnonzero(labels == index), tag=index)
        )
    mesh.point_data["freesurfer:annot"] = labels


def read(filename, curv=None, annot=None):
    c = _Cursor(filename, "surface")
    magic = c.magic()
    if magic in (0xFFFFFF, 0xFFFFFD):
        raise ReadError(
            f"meshio++: freesurfer: '{filename}' is a quad surface; meshio++ reads "
            "triangle surfaces"
        )
    if magic != 0xFFFFFE:
        raise ReadError(
            f"meshio++: freesurfer: '{filename}' is not a FreeSurfer surface"
        )
    end = c.data.find(b"\n\n", c.at)
    if end < 0:
        raise ReadError(f"meshio++: freesurfer: '{filename}' has no creator line")
    created = c.data[c.at : end]
    c.at = end + 2
    nv, nf = c.i32(), c.i32()
    limit = len(c.data) // 12
    for what, n in (("vertex", nv), ("face", nf)):
        if n < 0 or n > limit:
            raise ReadError(
                f"meshio++: freesurfer: {what} count {n} exceeds what the file holds"
            )
    points = c.f32(3 * nv).astype(np.float64).reshape(nv, 3)
    faces = c.i32(3 * nf).reshape(nf, 3) if nf else np.zeros((0, 3), np.int64)
    if faces.size and (faces.min() < 0 or faces.max() >= nv):
        raise ReadError("meshio++: freesurfer: a face names a vertex out of range")
    field_data = {"freesurfer:created_by": np.frombuffer(created, np.uint8).copy()}
    if c.at < len(c.data):
        field_data["freesurfer:tail"] = np.frombuffer(c.data[c.at :], np.uint8).copy()
    mesh = Mesh(points, [("triangle", faces)] if nf else [], field_data=field_data)
    for path in curv or []:
        _add_curv(mesh, path)
    if annot:
        _add_annot(mesh, annot)
    return mesh


def write(filename, mesh):
    if len(mesh.cells) != 1 or mesh.cells[0].type != "triangle":
        raise WriteError(
            "meshio++: freesurfer: a FreeSurfer surface is one triangle block"
        )
    tri = np.asarray(mesh.cells[0].data, dtype=np.int64)
    created = mesh.field_data.get("freesurfer:created_by")
    created = (
        b"created by meshio++"
        if created is None
        else np.asarray(created, np.uint8).tobytes()
    )
    pts = np.asarray(mesh.points, dtype=np.float64)
    full = np.zeros((len(pts), 3))
    full[:, : min(3, pts.shape[1])] = pts[:, :3]
    out = [
        b"\xff\xff\xfe",
        created + b"\n\n",
        struct.pack(">2i", len(pts), len(tri)),
        full.astype(">f4").tobytes(),
        tri.astype(">i4").tobytes(),
    ]
    tail = mesh.field_data.get("freesurfer:tail")
    if tail is not None:
        out.append(np.asarray(tail, np.uint8).tobytes())
    with open(filename, "wb") as f:
        f.write(b"".join(out))

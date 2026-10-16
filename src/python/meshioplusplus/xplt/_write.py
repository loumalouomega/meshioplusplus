"""Write FEBio plot files (``.xplt``).

The pure-Python twin of ``write_xplt`` in ``src/cpp/src/formats/xplt.cpp``: both
engines write the same bytes. The file is the layout ``_xplt.py`` reads and
FEBio Studio opens (plot version 0x35, little-endian, uncompressed): the header
and the dictionary of variables, the mesh (nodes with ids, one domain per cell
block, node sets, element sets, parts) and one state per step.

Variables: point data are nodal variables, cell data are per-element domain
variables (``FMT_ITEM``) and numeric field data are global variables. Facet
surfaces, edges, side regions and vertices are not written.
"""

from __future__ import annotations

import struct

import numpy as np

from .._common import warn
from .._exceptions import WriteError
from .._files import is_buffer
from .._node_order import from_meshio
from .._regions import canonical_regions
from ._xplt import MAGIC

__all__ = ["write", "SeriesWriter"]

_VERSION = 0x35
_SOFTWARE = "meshio++"
# meshio++ cell type -> (plot element type, nodes); FEBioPlotFile.h's Elem_Type.
_ELEM = {
    "hexahedron": 0,
    "wedge": 1,
    "tetra": 2,
    "quad": 3,
    "triangle": 4,
    "line": 5,
    "hexahedron20": 6,
    "tetra10": 7,
    "hexahedron27": 9,
    "triangle6": 10,
    "quad8": 11,
    "quad9": 12,
    "wedge15": 13,
    "triangle10": 15,
    "pyramid": 16,
    "pyramid13": 18,
    "line3": 19,
}
_GROUPS = (
    # (dictionary tag, state tag, storage format)
    ("global", 0x01021000, 0x02020100, 3),
    ("node", 0x01023000, 0x02020300, 0),
    ("domain", 0x01024000, 0x02020400, 1),
)
# values per entry -> Var_Type; anything else is an ARRAY of that size (6).
_TYPE_OF_WIDTH = {1: 0, 3: 1, 6: 2, 9: 5}


def _u32(v):
    return struct.pack("<I", v & 0xFFFFFFFF)


def _chunk(ident, payload=b""):
    return _u32(ident) + _u32(len(payload)) + payload


def _string(text):
    raw = text.encode("latin-1", "replace")
    return struct.pack("<i", len(raw)) + raw


def _name64(text):
    raw = text.encode("latin-1", "replace")[:63]
    return raw + b"\0" * (64 - len(raw))


def _f32(values):
    return np.asarray(values, dtype="<f4").tobytes()


def _i32(values):
    return np.asarray(values, dtype="<i4").tobytes()


class _Var:
    def __init__(self, name, width):
        self.name = name
        self.width = width

    def key(self):
        return (self.name, self.width)


class _Model:
    """The mesh part of a plot file and the layout of its variables."""

    def __init__(self, mesh):
        points = np.asarray(mesh.points, dtype=np.float64)
        n = len(points)
        points = points.reshape(n, -1)
        self.points = np.zeros((n, 3))
        self.points[:, : min(3, points.shape[1])] = points[:, :3]
        self.n = n
        self.blocks = []  # (block index, cell type, conn in file order, first cell)
        skipped = {}
        self.new_index = {}  # global cell number -> element id (1-based)
        base = 0
        eid = 1
        for b, block in enumerate(mesh.cells):
            data = np.asarray(block.data, dtype=np.int64)
            if block.type not in _ELEM:
                if len(data):
                    skipped[block.type] = skipped.get(block.type, 0) + len(data)
                base += len(data)
                continue
            if len(data) and (data.min() < 0 or data.max() >= n):
                raise WriteError(
                    f"FEBio .xplt: a '{block.type}' cell refers to a point the mesh "
                    "does not have"
                )
            conn = from_meshio("febio", block.type, data)
            for r in range(len(data)):
                self.new_index[base + r] = eid
                eid += 1
            self.blocks.append((b, block.type, conn, base))
            base += len(data)
        self.ncells = base
        for ctype, count in skipped.items():
            warn(f"FEBio .xplt: {count} '{ctype}' cell(s) are not written")
        self.regions = canonical_regions(mesh)
        self._names(mesh)
        self.vars = {"global": [], "node": [], "domain": []}
        self._variables(mesh)
        self.signature = (
            n,
            [(t, c.shape, s) for _, t, c, s in self.blocks],
            {g: [v.key() for v in vs] for g, vs in self.vars.items()},
        )

    def _names(self, mesh):
        """Per domain: (name, part id)."""
        self.domains = []
        by_range = {}
        for r in self.regions:
            if r.kind == "cell":
                by_range.setdefault(tuple(np.asarray(r.entries).tolist()), r)
        for k, (_, ctype, conn, first) in enumerate(self.blocks):
            count = len(conn)
            match = by_range.get(tuple(range(first, first + count)))
            if match is not None and count:
                part = match.tag if match.tag >= 0 else k + 1
                self.domains.append((match.name, part, f"Part{k + 1}"))
            else:
                self.domains.append(("", k + 1, f"Part{k + 1}"))

    def _variables(self, mesh):
        for name in sorted(mesh.point_data):
            a = np.asarray(mesh.point_data[name])
            width = self._width(a, self.n, f"point data '{name}'")
            if width:
                self.vars["node"].append(_Var(name, width))
        for name in sorted(mesh.cell_data):
            if name.startswith("xplt:"):
                continue
            width = None
            for b, _, conn, _ in self.blocks:
                a = np.asarray(mesh.cell_data[name][b])
                w = self._width(a, len(conn), f"cell data '{name}'")
                if w is None or (width is not None and w != width):
                    width = None
                    break
                width = w
            else:
                if width:
                    self.vars["domain"].append(_Var(name, width))
                continue
            warn(f"FEBio .xplt: cell data '{name}' is not written")
        for name in sorted(mesh.field_data):
            if name == "meshio:time" or name.startswith("xplt:"):
                continue
            a = np.asarray(mesh.field_data[name])
            if a.dtype.kind in "fiub" and a.size:
                self.vars["global"].append(_Var(name, int(a.size)))
            else:
                warn(f"FEBio .xplt: field data '{name}' is not written")

    @staticmethod
    def _width(a, rows, what):
        if a.dtype.kind not in "fiub":
            warn(f"FEBio .xplt: {what} is not written")
            return None
        if rows == 0:
            return 1 if a.ndim <= 1 else int(a.shape[1]) if a.ndim == 2 else None
        if a.size % rows or a.ndim > 2:
            warn(f"FEBio .xplt: {what} is not written")
            return None
        return int(a.size // rows)

    # -- the mesh chunks -----------------------------------------------------------
    def header_and_dictionary(self):
        header = _chunk(
            0x01010000,
            _chunk(0x01010001, _u32(_VERSION))
            + _chunk(0x01010004, _u32(0))
            + _chunk(0x01010006, _string(_SOFTWARE)),
        )
        dictionary = b""
        for group, tag, _, fmt in _GROUPS:
            items = b""
            for v in self.vars[group]:
                vtype = _TYPE_OF_WIDTH.get(v.width, 6)
                asize = v.width if vtype == 6 else 0
                items += _chunk(
                    0x01020001,
                    _chunk(0x01020002, _u32(vtype))
                    + _chunk(0x01020003, _u32(fmt))
                    + _chunk(0x01020005, _u32(asize))
                    + _chunk(0x01020004, _name64(v.name)),
                )
            if items:
                dictionary += _chunk(tag, items)
        return _chunk(0x01000000, header + _chunk(0x01020000, dictionary))

    def mesh_chunk(self):
        coords = b"".join(_u32(i + 1) + _f32(self.points[i]) for i in range(self.n))
        nodes = _chunk(
            0x01041000,
            _chunk(
                0x01041100,
                _chunk(0x01041101, _u32(self.n)) + _chunk(0x01041102, _u32(3)),
            )
            + _chunk(0x01041200, coords),
        )
        domains = b""
        for (_, ctype, conn, first), (name, part, _) in zip(self.blocks, self.domains):
            head = (
                _chunk(0x01042102, _u32(_ELEM[ctype]))
                + _chunk(0x01042103, _u32(part))
                + _chunk(0x01032104, _u32(len(conn)))
            )
            if name:
                head += _chunk(0x01032105, _string(name))
            elements = b"".join(
                _chunk(0x01042201, _i32([self.new_index[first + r]] + row.tolist()))
                for r, row in enumerate(conn)
            )
            domains += _chunk(
                0x01042100, _chunk(0x01042101, head) + _chunk(0x01042200, elements)
            )
        body = nodes + _chunk(0x01042000, domains)

        def sets(kind, section, item, hdr, ident, size, name, lst, members):
            out = b""
            count = 0
            for r in self.regions:
                if r.kind != kind:
                    continue
                ids = members(np.asarray(r.entries).tolist())
                count += 1
                out += _chunk(
                    item,
                    _chunk(
                        hdr,
                        _chunk(ident, _u32(count))
                        + _chunk(size, _u32(len(ids)))
                        + _chunk(name, _string(r.name)),
                    )
                    + _chunk(lst, _i32(ids)),
                )
            return _chunk(section, out) if out else b""

        body += sets(
            "point",
            0x01044000,
            0x01044100,
            0x01044101,
            0x01044102,
            0x01044104,
            0x01044103,
            0x01044200,
            lambda e: [i for i in e if 0 <= i < self.n],
        )
        body += sets(
            "cell",
            0x01046000,
            0x01046100,
            0x01046101,
            0x01046102,
            0x01046104,
            0x01046103,
            0x01046200,
            lambda e: [self.new_index[i] for i in e if i in self.new_index],
        )
        parts = b"".join(
            _chunk(
                0x01045100,
                _chunk(0x01045101, _u32(part)) + _chunk(0x01045102, _name64(pname)),
            )
            for (_, part, pname) in self.domains
        )
        if parts:
            body += _chunk(0x01045000, parts)
        return _chunk(0x01040000, body)

    # -- one state -------------------------------------------------------------------
    def state_chunk(self, mesh, time, status=0):
        header = _chunk(
            0x02010000,
            _chunk(0x02010002, _f32([time])) + _chunk(0x02010003, _u32(status)),
        )
        data = b""
        for group, _, tag, _ in _GROUPS:
            variables = b""
            for k, v in enumerate(self.vars[group]):
                if group == "global":
                    regions = [(0, np.asarray(mesh.field_data[v.name]).reshape(-1))]
                elif group == "node":
                    regions = [(0, np.asarray(mesh.point_data[v.name]).reshape(-1))]
                else:
                    regions = [
                        (d + 1, np.asarray(mesh.cell_data[v.name][b]).reshape(-1))
                        for d, (b, _, _, _) in enumerate(self.blocks)
                    ]
                payload = b"".join(
                    _chunk(rid, _f32(values.astype(np.float64)))
                    for rid, values in regions
                )
                variables += _chunk(
                    0x02020001,
                    _chunk(0x02020002, _u32(k + 1)) + _chunk(0x02020003, payload),
                )
            if variables:
                data += _chunk(tag, variables)
        return _chunk(0x02000000, header + _chunk(0x02020000, data))


def _time_of(mesh):
    fd = getattr(mesh, "field_data", {}) or {}
    if "meshio:time" in fd:
        return float(np.asarray(fd["meshio:time"]).reshape(-1)[0])
    return 0.0


def _status_of(mesh):
    fd = getattr(mesh, "field_data", {}) or {}
    if "xplt:status" in fd:
        return int(np.asarray(fd["xplt:status"]).reshape(-1)[0])
    return 0


def _assemble(model, states):
    return (
        _u32(MAGIC)
        + model.header_and_dictionary()
        + model.mesh_chunk()
        + b"".join(states)
    )


def write(filename, mesh):
    """Write one FEBio plot file holding one state: the mesh and its data at
    ``field_data["meshio:time"]`` (0 when absent)."""
    if is_buffer(filename, "w"):
        raise WriteError("FEBio .xplt: cannot be written to a stream")
    model = _Model(mesh)
    state = model.state_chunk(mesh, _time_of(mesh), _status_of(mesh))
    with open(filename, "wb") as fh:
        fh.write(_assemble(model, [state]))


class SeriesWriter:
    """One plot file with a state per step; the mesh and the variables are those
    of the first step, which every later step must share."""

    def __init__(self, filename):
        self._path = filename
        if is_buffer(filename, "w"):
            raise WriteError("FEBio .xplt: cannot be written to a stream")
        self._model = None
        self._states = []
        self._closed = False

    def write(self, time, mesh):
        if self._model is None:
            self._model = _Model(mesh)
        else:
            other = _Model(mesh)
            if other.signature != self._model.signature:
                raise WriteError(
                    f"FEBio .xplt: series: step {len(self._states)}'s points, cells "
                    "or variables differ from the first step's; a plot file holds "
                    "one mesh -- write one file per step with '{step}'"
                )
        self._states.append(
            self._model.state_chunk(mesh, float(time), _status_of(mesh))
        )

    def close(self):
        if self._closed or self._model is None:
            self._closed = True
            return
        self._closed = True
        with open(self._path, "wb") as fh:
            fh.write(_assemble(self._model, self._states))

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False

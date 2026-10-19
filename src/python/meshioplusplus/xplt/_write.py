"""Write FEBio plot files (``.xplt``).

The pure-Python twin of ``write_xplt`` in ``src/cpp/src/formats/xplt.cpp``: both
engines write the same bytes. The file is the layout ``_xplt.py`` reads and
FEBio Studio opens (plot version 0x35, little-endian, uncompressed): the header
and the dictionary of variables, the mesh (nodes with ids, one domain per cell
block, node sets, element sets, parts) and one state per step.

Variables: point data are nodal variables, cell data are per-element domain
variables (``FMT_ITEM``) and numeric field data are global variables.

Surfaces, edges and facet sets are the inverse of what the reader builds. A
block of ``triangle``/``quad``-family (``line``/``line3``) cells carrying a
positive ``xplt:surface`` (``xplt:edge``) cell-data id is not a domain: the
blocks of one id are one surface (edge) named after the cell region holding
exactly its cells, and the id is not a variable. Every ``side`` region that no
written surface is named after becomes a facet set of the faces it names.

Surface and edge variables are the cell data that has finite values on the cells
of a surface (edge): one per-facet (``FMT_ITEM``) variable of that name, whose
region ``k`` is surface (edge) ``k`` and which is left out of a surface holding
only NaN. An array with a finite value on a domain cell is a domain variable
instead (a name is one or the other, as on reading). Point-valued surface
variables the reader averages onto the points come back as nodal variables, not
as these.

Compression (``compress=True``) writes every state as its own zlib stream, as
FEBio does, and sets the header's compression flag. Vertices are not written.
"""

from __future__ import annotations

import struct
import zlib

import numpy as np

from .._common import warn
from .._exceptions import WriteError
from .._facets import facet_nodes
from .._files import is_buffer
from .._node_order import from_meshio
from .._regions import canonical_regions
from ._xplt import MAGIC

__all__ = ["write", "SeriesWriter"]

_VERSION = 0x35
_LEVEL = 6  # zlib's default; the C++ twin deflates at the same level
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
# Facet node counts FEBio's surfaces take, and the cell types that hold them.
_SURFACE_NODES = {
    "triangle": 3,
    "quad": 4,
    "triangle6": 6,
    "triangle7": 7,
    "quad8": 8,
    "quad9": 9,
    "triangle10": 10,
}
_EDGE_NODES = {"line": 2, "line3": 3}
_FACET_WIDTH = {3, 4, 6, 7, 8, 9, 10}
_FACETSET_PAD = 10  # a facet set's leaves are padded to this many nodes
_EDGE_PAD = 3
_GROUPS = (
    # (dictionary tag, state tag, storage format)
    ("global", 0x01021000, 0x02020100, 3),
    ("node", 0x01023000, 0x02020300, 0),
    ("domain", 0x01024000, 0x02020400, 1),
    ("surface", 0x01025000, 0x02020500, 1),
    ("edge", 0x01026000, 0x02020600, 1),
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


def _marker(mesh, key, b, count):
    """The surface (edge) id of block ``b``: its ``key`` cell data when every
    cell of the block carries the same positive integer, else 0."""
    arrays = mesh.cell_data.get(key)
    if arrays is None or b >= len(arrays) or not count:
        return 0
    a = np.asarray(arrays[b]).reshape(-1)
    if a.dtype.kind not in "iu" or len(a) != count or (a != a[0]).any():
        return 0
    return max(int(a[0]), 0)


class _Var:
    def __init__(self, name, width, nodal=False):
        self.name = name
        self.width = width
        # a surface or edge variable stored per node (FMT_NODE), not per facet
        self.nodal = nodal

    def key(self):
        return (self.name, self.width, self.nodal)


class _Model:
    """The mesh part of a plot file and the layout of its variables."""

    def __init__(self, mesh, compress=False):
        self.compress = bool(compress)
        points = np.asarray(mesh.points, dtype=np.float64)
        n = len(points)
        points = points.reshape(n, -1)
        self.points = np.zeros((n, 3))
        self.points[:, : min(3, points.shape[1])] = points[:, :3]
        self.n = n
        self.blocks = []  # (block index, cell type, conn in file order, first cell)
        skipped = {}
        self.new_index = {}  # global cell number -> element id (1-based)
        carriers = {"surface": {}, "edge": {}}  # id -> [(first cell, rows)]
        self.carrier_cells = set()
        base = 0
        eid = 1
        for b, block in enumerate(mesh.cells):
            data = np.asarray(block.data, dtype=np.int64)
            carrier = None
            if block.type in _SURFACE_NODES:
                ident = _marker(mesh, "xplt:surface", b, len(data))
                carrier = ("surface", ident) if ident else None
            elif block.type in _EDGE_NODES:
                ident = _marker(mesh, "xplt:edge", b, len(data))
                carrier = ("edge", ident) if ident else None
            if carrier is not None:
                if data.min() < 0 or data.max() >= n:
                    raise WriteError(
                        f"FEBio .xplt: a '{block.type}' cell refers to a point the "
                        "mesh does not have"
                    )
                carriers[carrier[0]].setdefault(carrier[1], []).append((base, data, b))
                self.carrier_cells.update(range(base, base + len(data)))
                base += len(data)
                continue
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
        self.regions = canonical_regions(mesh)
        self.surfaces = self._carriers(carriers["surface"], "Surface")
        self.edges = self._carriers(carriers["edge"], "Edge")
        # per surface (edge), in id order: the cell blocks that hold its facets
        self.sources = {
            kind: [[b for _, _, b in carriers[kind][i]] for i in sorted(carriers[kind])]
            for kind in ("surface", "edge")
        }
        self.facetsets = self._facetsets(mesh)
        for ctype, count in skipped.items():
            warn(f"FEBio .xplt: {count} '{ctype}' cell(s) are not written")
        self._names(mesh)
        self.vars = {g: [] for g, _, _, _ in _GROUPS}
        self._variables(mesh)
        # Points, cells, facets and sets: a change writes a new mesh section
        # into a series. The dictionary (variable names and widths) cannot change.
        self.mesh_signature = (
            n,
            [(t, c.shape, s, c.tobytes()) for _, t, c, s in self.blocks],
            [(nm, [(nn, tuple(f)) for nn, f in fs]) for nm, fs in self.surfaces],
            [(nm, [(nn, tuple(f)) for nn, f in fs]) for nm, fs in self.edges],
            [(nm, [(nn, tuple(f)) for nn, f in fs]) for nm, fs in self.facetsets],
        )
        self.dictionary_signature = {
            g: [v.key() for v in vs] for g, vs in self.vars.items()
        }

    def _carriers(self, by_id, fallback):
        """Per surface (edge), in id order: (name, [(nodes per facet, nodes)])."""
        by_cells = {}
        for r in self.regions:
            if r.kind == "cell":
                by_cells.setdefault(tuple(sorted(np.asarray(r.entries).tolist())), r)
        out = []
        for k, ident in enumerate(sorted(by_id), start=1):
            cells, facets = [], []
            for first, data, _ in by_id[ident]:
                cells.extend(range(first, first + len(data)))
                facets.extend((data.shape[1], row) for row in data.tolist())
            match = by_cells.get(tuple(sorted(cells)))
            out.append((match.name if match else f"{fallback}{k}", facets))
        return out

    def _facetsets(self, mesh):
        """Per side region no written surface is named after: its faces."""
        taken = {name for name, _ in self.surfaces}
        out = []
        for r in self.regions:
            if r.kind != "side" or r.name in taken:
                continue
            facets, dropped = [], 0
            for cell, facet in np.asarray(r.entries).reshape(-1, 2).tolist():
                found = facet_nodes(mesh, cell, facet)
                if found is None or len(found[1]) not in _FACET_WIDTH:
                    dropped += 1
                    continue
                facets.append((len(found[1]), found[1]))
            if dropped:
                warn(
                    f"FEBio .xplt: {dropped} edge or missing facet(s) of the side "
                    f"region '{r.name}' are not written"
                )
            out.append((r.name, facets))
        return out

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
            if not width:
                continue
            # The reader marks a point array that came from surface (edge) node
            # values; it goes back to that group when the mesh still has such facets.
            group = next(
                (
                    g
                    for g in ("surface", "edge")
                    if f"xplt:{g}_node:{name}" in mesh.field_data and self.sources[g]
                ),
                "node",
            )
            self.vars[group].append(_Var(name, width, group != "node"))
        for name in sorted(mesh.cell_data):
            if name.startswith("xplt:"):
                continue
            arrays = mesh.cell_data[name]
            has_domain = self._has_values(arrays, [b for b, _, _, _ in self.blocks])
            # a variable with a value on a domain is a domain variable
            carrier = {
                g: 0 if has_domain else self._carrier_width(mesh, name, g)
                for g in ("surface", "edge")
            }
            written = False
            for group, width in carrier.items():
                if width:
                    self.vars[group].append(_Var(name, width))
                    written = True
            width = None
            for b, _, conn, _ in self.blocks:
                a = np.asarray(arrays[b])
                w = self._width(a, len(conn), f"cell data '{name}'", not written)
                if w is None or (width is not None and w != width):
                    width = None
                    break
                width = w
            else:
                if width and not written:
                    self.vars["domain"].append(_Var(name, width))
                continue
            if not written:
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
    def _has_values(arrays, blocks):
        """Whether any of the given blocks holds a finite (or non-float) value."""
        for b in blocks:
            if b >= len(arrays):
                continue
            a = np.asarray(arrays[b])
            if a.dtype.kind == "f":
                if np.isfinite(a).any():
                    return True
            elif a.dtype.kind in "iub" and a.size:
                return True
        return False

    def _carrier_width(self, mesh, name, group):
        """Values per facet of the cell data ``name`` on the surfaces (edges), or
        0 when no facet has a finite value or the blocks disagree."""
        arrays = mesh.cell_data[name]
        blocks = [b for srcs in self.sources[group] for b in srcs]
        if not blocks or not self._has_values(arrays, blocks):
            return 0
        width = None
        for b in blocks:
            if b >= len(arrays):
                return 0
            count = len(mesh.cells[b].data)
            a = np.asarray(arrays[b])
            if a.dtype.kind not in "fiub" or a.ndim > 2 or a.size % count:
                return 0
            w = a.size // count
            if width not in (None, w):
                return 0
            width = w
        return width or 0

    @staticmethod
    def _width(a, rows, what, loud=True):
        if a.dtype.kind not in "fiub":
            if loud:
                warn(f"FEBio .xplt: {what} is not written")
            return None
        if rows == 0:
            return 1 if a.ndim <= 1 else int(a.shape[1]) if a.ndim == 2 else None
        if a.size % rows or a.ndim > 2:
            if loud:
                warn(f"FEBio .xplt: {what} is not written")
            return None
        return int(a.size // rows)

    # -- the mesh chunks -----------------------------------------------------------
    def header_and_dictionary(self):
        header = _chunk(
            0x01010000,
            _chunk(0x01010001, _u32(_VERSION))
            + _chunk(0x01010004, _u32(1 if self.compress else 0))
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
                    + _chunk(0x01020003, _u32(0 if v.nodal else fmt))
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
                entries = np.asarray(r.entries).tolist()
                ids = members(entries)
                if kind == "cell" and entries and not ids:
                    if all(i in self.carrier_cells for i in entries):
                        continue  # a surface or edge, written as one
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
        # (section, item, header, id, size, name), (widest, list, leaf), padding
        body += self._facet_section(
            self.surfaces,
            (0x01043000, 0x01043100, 0x01043101, 0x01043102, 0x01043103, 0x01043104),
            (0x01043105, 0x01043200, 0x01043201),
            None,
        )
        body += self._facet_section(
            self.facetsets,
            (0x01047000, 0x01047100, 0x01047101, 0x01047102, 0x01047104, 0x01047103),
            (0x01047105, 0x01047200, 0x01047201),
            _FACETSET_PAD,
        )
        body += self._facet_section(
            self.edges,
            (0x01048000, 0x01048100, 0x01048101, 0x01048102, 0x01048103, 0x01048104),
            (0x01048105, 0x01048200, 0x01048201),
            _EDGE_PAD,
        )
        return _chunk(0x01040000, body)

    @staticmethod
    def _facet_section(items, head_tags, list_tags, pad):
        """A section of surfaces, facet sets or edges: FEBio's header (id, size,
        name, widest facet) and a list of ``[facet id, node count, nodes]``
        leaves, padded to ``pad`` nodes (the widest facet when ``None``)."""
        if not items:
            return b""
        section, item, header, t_id, t_size, t_name = head_tags
        t_widest, t_list, t_leaf = list_tags
        out = b""
        for k, (name, facets) in enumerate(items, start=1):
            width = pad or max((nn for nn, _ in facets), default=0)
            head = (
                _chunk(t_id, _u32(k))
                + _chunk(t_size, _u32(len(facets)))
                + _chunk(t_name, _string(name))
                + _chunk(t_widest, _u32(width))
            )
            leaves = b"".join(
                _chunk(t_leaf, _i32([f + 1, nn] + list(nodes) + [0] * (width - nn)))
                for f, (nn, nodes) in enumerate(facets)
            )
            out += _chunk(item, _chunk(header, head) + _chunk(t_list, leaves))
        return _chunk(section, out)

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
                elif group == "domain":
                    regions = [
                        (d + 1, np.asarray(mesh.cell_data[v.name][b]).reshape(-1))
                        for d, (b, _, _, _) in enumerate(self.blocks)
                    ]
                elif v.nodal:
                    # region k is surface (edge) k: the point values at its nodes,
                    # each node once, in the order its facets first list them
                    all_values = np.asarray(
                        mesh.point_data[v.name], dtype=np.float64
                    ).reshape(self.n, v.width)
                    facets = self.surfaces if group == "surface" else self.edges
                    regions = []
                    for s, (_, items) in enumerate(facets):
                        nodes = list(dict.fromkeys(p for _, f in items for p in f))
                        if not nodes:
                            continue
                        values = all_values[nodes].reshape(-1)
                        if not np.isnan(values).all():
                            regions.append((s + 1, values))
                else:
                    # region k is surface (edge) k: its facets' values in order,
                    # none for a surface that has no value
                    regions = []
                    for s, srcs in enumerate(self.sources[group]):
                        values = np.concatenate(
                            [
                                np.asarray(
                                    mesh.cell_data[v.name][b], dtype=np.float64
                                ).reshape(-1)
                                for b in srcs
                            ]
                        )
                        if not np.isnan(values).all():
                            regions.append((s + 1, values))
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


def _pack(model, state):
    """A state as the file stores it: with compression, its own zlib stream."""
    return zlib.compress(state, _LEVEL) if model.compress else state


def _assemble(model, states, mesh=None):
    return (
        _u32(MAGIC)
        + model.header_and_dictionary()
        + (model.mesh_chunk() if mesh is None else mesh)
        + b"".join(states)
    )


def write(filename, mesh, compress=False):
    """Write one FEBio plot file holding one state: the mesh and its data at
    ``field_data["meshio:time"]`` (0 when absent). ``compress=True`` stores the
    state as a zlib stream."""
    if is_buffer(filename, "w"):
        raise WriteError("FEBio .xplt: cannot be written to a stream")
    model = _Model(mesh, compress)
    state = _pack(model, model.state_chunk(mesh, _time_of(mesh), _status_of(mesh)))
    with open(filename, "wb") as fh:
        fh.write(_assemble(model, [state]))


class SeriesWriter:
    """One plot file with a state per step. The variables are those of the first
    step, which every later step must share; a step whose mesh differs from the
    previous one (a remeshed run) writes a new mesh section before its state."""

    def __init__(self, filename, compress=False):
        self._path = filename
        self._compress = bool(compress)
        if is_buffer(filename, "w"):
            raise WriteError("FEBio .xplt: cannot be written to a stream")
        self._model = None  # the first step's: it fixes the header and dictionary
        self._mesh_signature = None
        self._first_mesh = b""
        self._states = []
        self._closed = False

    def write(self, time, mesh):
        step = _Model(mesh, self._compress)
        chunks = b""
        if self._model is None:
            self._model = step
            self._first_mesh = step.mesh_chunk()
            self._mesh_signature = step.mesh_signature
        else:
            if step.dictionary_signature != self._model.dictionary_signature:
                raise WriteError(
                    f"FEBio .xplt: series: step {len(self._states)}'s variables "
                    "differ from the first step's; a plot file holds one variable "
                    "dictionary -- write one file per step with '{step}'"
                )
            if step.mesh_signature != self._mesh_signature:
                chunks = _pack(step, step.mesh_chunk())
                self._mesh_signature = step.mesh_signature
        state = step.state_chunk(mesh, float(time), _status_of(mesh))
        self._states.append(chunks + _pack(step, state))

    def close(self):
        if self._closed or self._model is None:
            self._closed = True
            return
        self._closed = True
        with open(self._path, "wb") as fh:
            fh.write(_assemble(self._model, self._states, self._first_mesh))

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False

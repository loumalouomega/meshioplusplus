"""I/O for libMesh ``.xda`` (ASCII) / ``.xdr`` (XDR binary) meshes.

The pure-Python twin of ``src/cpp/src/formats/libmesh.cpp``: both engines read
the same meshes.

Both encodings carry one value stream (libMesh's ``XdrIO``): a version string,
the element and node counts, four "inline or not" flags, per-field integer sizes
(0.9.2+), subdomain names, one connectivity block per refinement level, the
coordinates, then the side sets, node sets (0.9.2+) and edge and shell-face sets
(1.1.0+). ``.xdr`` is big-endian XDR.

Only active (leaf) elements become cells. The subdomain id is the
``libmesh:subdomain`` cell data and a cell region per subdomain; side sets
become side regions (carried down to refined children), node sets point
regions. HEX20/HEX27/PRISM15/PRISM18 use the ``"libmesh"`` node-order tables.
"""

import math
import struct

import numpy as np

from .. import _provenance
from .._common import warn
from .._exceptions import ReadError
from .._facets import FacetIndex
from .._files import open_file
from .._mesh import Mesh, topological_dimension
from .._node_order import node_order
from .._regions import Region

__all__ = ["read"]

# libMesh ElemType -> (file nodes, meshio++ type or None, nodes kept, shape, name).
_TYPES = [
    (2, "line", 2, "edge", "EDGE2"),
    (3, "line3", 3, "edge", "EDGE3"),
    (4, "line4", 4, "edge", "EDGE4"),
    (3, "triangle", 3, "tri", "TRI3"),
    (6, "triangle6", 6, "tri", "TRI6"),
    (4, "quad", 4, "quad", "QUAD4"),
    (8, "quad8", 8, "quad", "QUAD8"),
    (9, "quad9", 9, "quad", "QUAD9"),
    (4, "tetra", 4, "tet", "TET4"),
    (10, "tetra10", 10, "tet", "TET10"),
    (8, "hexahedron", 8, "hex", "HEX8"),
    (20, "hexahedron20", 20, "hex", "HEX20"),
    (27, "hexahedron27", 27, "hex", "HEX27"),
    (6, "wedge", 6, "prism", "PRISM6"),
    (15, "wedge15", 15, "prism", "PRISM15"),
    (18, "wedge18", 18, "prism", "PRISM18"),
    (5, "pyramid", 5, "pyramid", "PYRAMID5"),
    (13, "pyramid13", 13, "pyramid", "PYRAMID13"),
    (14, "pyramid14", 14, "pyramid", "PYRAMID14"),
    (2, None, 0, None, "INFEDGE2"),
    (4, None, 0, None, "INFQUAD4"),
    (6, None, 0, None, "INFQUAD6"),
    (8, None, 0, None, "INFHEX8"),
    (16, None, 0, None, "INFHEX16"),
    (18, None, 0, None, "INFHEX18"),
    (6, None, 0, None, "INFPRISM6"),
    (12, None, 0, None, "INFPRISM12"),
    (1, "vertex", 1, "point", "NODEELEM"),
    (0, None, 0, None, "REMOTEELEM"),
    (3, "triangle", 3, "tri", "TRI3SUBDIVISION"),
    (3, "triangle", 3, "tri", "TRISHELL3"),
    (4, "quad", 4, "quad", "QUADSHELL4"),
    (8, "quad8", 8, "quad", "QUADSHELL8"),
    (7, "triangle7", 7, "tri", "TRI7"),
    (14, "tetra10", 10, "tet", "TET14"),
    (20, "wedge18", 18, "prism", "PRISM20"),
    (21, "wedge18", 18, "prism", "PRISM21"),
    (18, "pyramid14", 14, "pyramid", "PYRAMID18"),
    (9, "quad9", 9, "quad", "QUADSHELL9"),
    (0, None, 0, None, "C0POLYGON"),
    (0, None, 0, None, "C0POLYHEDRON"),
]

# Corner nodes of each side in libMesh's side numbering (Hex8::side_nodes_map ...).
_SIDES = {
    "tri": [(0, 1), (1, 2), (2, 0)],
    "quad": [(0, 1), (1, 2), (2, 3), (3, 0)],
    "tet": [(0, 2, 1), (0, 1, 3), (1, 2, 3), (2, 0, 3)],
    "hex": [
        (0, 3, 2, 1),
        (0, 1, 5, 4),
        (1, 2, 6, 5),
        (2, 3, 7, 6),
        (3, 0, 4, 7),
        (4, 5, 6, 7),
    ],
    "prism": [(0, 2, 1), (0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5), (3, 4, 5)],
    "pyramid": [(0, 1, 4), (1, 2, 4), (2, 3, 4), (3, 0, 4), (0, 3, 2, 1)],
}


class _Stream:
    """libMesh's ``Xdr`` in READ (ASCII) or DECODE (XDR) mode (xdr_cxx.C)."""

    def __init__(self, data, xdr):
        self.data = data
        self.xdr = xdr
        self.pos = 0

    def fail(self, what):
        if self.xdr:
            where = f"byte {self.pos}"
        else:
            where = "line %d" % (self.data.count(b"\n", 0, self.pos) + 1)
        raise ReadError(f"libMesh: {what} ({where})")

    # -- XDR ------------------------------------------------------------------
    def _take(self, n):
        if self.pos + n > len(self.data):
            raise ReadError(
                f"libMesh: file is truncated (needs {n} more bytes at offset "
                f"{self.pos} of {len(self.data)})"
            )
        chunk = self.data[self.pos : self.pos + n]
        self.pos += n
        return chunk

    def _bin_int(self, width):
        return struct.unpack(">Q" if width == 8 else ">I", self._take(width))[0]

    # -- ASCII ----------------------------------------------------------------
    def _line_end(self):
        eol = self.data.find(b"\n", self.pos)
        return len(self.data) if eol < 0 else eol

    def _skip_line(self):
        eol = self._line_end()
        self.pos = eol + 1 if eol < len(self.data) else eol

    def _token(self):
        data, n = self.data, len(self.data)
        pos = self.pos
        while pos < n and data[pos : pos + 1].isspace():
            pos += 1
        if pos >= n:
            self.pos = pos
            self.fail("file ends early")
        start = pos
        while pos < n and not data[pos : pos + 1].isspace():
            pos += 1
        self.pos = pos
        return data[start:pos].decode("latin-1")

    def _int(self):
        tok = self._token()
        try:
            return int(tok)
        except ValueError:
            self.fail(f"bad integer '{tok}'")

    # -- the value kinds --------------------------------------------------------
    def string(self):
        if self.xdr:
            n = self._bin_int(4)
            s = self._take(n)
            self._take((4 - n % 4) % 4)
            return s.split(b"\0", 1)[0].decode("latin-1")
        eol = self._line_end()
        s = self.data[self.pos : eol].decode("latin-1")
        self.pos = eol + 1 if eol < len(self.data) else eol
        s = s.split("\t", 1)[0]
        return s[:-1] if s.endswith("\r") else s

    def scalar(self, width):
        if self.xdr:
            return self._bin_int(width)
        v = self._int()
        self._skip_line()
        return v

    def stream_int(self, width):
        return self._bin_int(width) if self.xdr else self._int()

    def stream_ints(self, width, count):
        if self.xdr:
            fmt = ">%d%s" % (count, "Q" if width == 8 else "I")
            return list(struct.unpack(fmt, self._take(width * count)))
        return [self._int() for _ in range(count)]

    def stream_reals(self, count):
        if self.xdr:
            return np.frombuffer(self._take(8 * count), dtype=">f8").astype(np.float64)
        out = np.empty(count, dtype=np.float64)
        for k in range(count):
            tok = self._token()
            try:
                out[k] = float(tok)
            except ValueError:
                if "nan" in tok.lower():
                    out[k] = math.nan
                else:
                    self.fail(f"bad real '{tok}'")
        return out

    def int_vector(self, width):
        n = self.scalar(4)
        out = [self.stream_int(width) for _ in range(n)]
        if not self.xdr:
            self._skip_line()
        return out

    def string_vector(self):
        n = self.scalar(4)
        out = [self.string() if self.xdr else self._token() for _ in range(n)]
        if not self.xdr:
            self._skip_line()
        return out


def _has_any(version, versions):
    return any(v in version for v in versions)


def _name_map(io, hw):
    n = io.scalar(hw)
    if n == 0:
        return {}
    ids = io.int_vector(hw)
    names = io.string_vector()
    return dict(zip(ids, names))


def _parse(data, xdr):
    io = _Stream(data, xdr)
    f = {
        "elements": [],
        "coords": np.zeros(0),
        "subdomain_names": {},
        "sideset_names": {},
        "nodeset_names": {},
        "sides": [],
        "nodesets": [],
        "skipped_edge_bcs": 0,
        "inline_p": False,
    }
    version = io.string()
    if "libMesh" not in version:
        raise ReadError(
            f"libMesh: '{version}' is a legacy (pre-libMesh) mesh file, which libMesh "
            "itself no longer reads"
        )
    v092 = _has_any(version, ("0.9.2", "0.9.6", "1.1.0", "1.3.0", "1.8.0"))
    v096 = _has_any(version, ("0.9.6", "1.1.0", "1.3.0", "1.8.0"))
    v110 = _has_any(version, ("1.1.0", "1.3.0", "1.8.0"))
    v130 = _has_any(version, ("1.3.0", "1.8.0"))
    v180 = "1.8.0" in version
    hw = 8 if v130 else 4

    n_elem = io.scalar(hw)
    n_nodes = io.scalar(hw)
    bc_file = io.string()
    sid_file = io.string()
    pid_file = io.string()
    pl_file = io.string()

    sizes = [8, 0, 0, 0, 0, 0, 0, 0]  # type uid pid sid p eid side bid
    if v092:
        sizes = [io.scalar(hw) for _ in range(8)]
    field_width = sizes[0]
    # Pre-1.3.0 files were written with 32-bit connectivity (XdrIO::read).
    tw = 4 if (not v130 or field_width == 4) else 8
    read_uid = v092 and sizes[1] != 0

    n_elem_ints = n_node_ints = 0
    if v180:
        io.scalar(hw)  # extra integer size
        n_node_ints = len(io.string_vector())
        n_elem_ints = len(io.string_vector())
        codes = io.int_vector(tw)
        for _ in codes:
            io.int_vector(tw)
    if v092:
        f["subdomain_names"] = _name_map(io, hw)
    if n_elem == 0:
        return f

    read_pid = pid_file == "."
    read_sid = sid_file == "."
    read_p = pl_file == "."
    f["inline_p"] = read_p
    elements = f["elements"]
    at_level = done = 0
    level = -1
    for e in range(n_elem):
        if done == at_level:
            at_level = io.scalar(tw)
            done = 0
            level += 1
        done += 1
        code = io.stream_int(tw)
        info = _TYPES[code] if code < len(_TYPES) else None
        if info is None or info[0] == 0:
            name = f" ({info[4]})" if info else ""
            io.fail(f"element {e} has unsupported type {code}{name}")
        if read_uid:
            io.stream_int(tw)
        parent = -1
        if level > 0:
            parent = io.stream_int(tw)
            if parent < 0 or parent >= e:
                io.fail(
                    f"element {e} names parent {parent}, which is not an earlier element"
                )
        if read_pid:
            io.stream_int(tw)
        sid = io.stream_int(tw) if read_sid else 0
        plev = io.stream_int(tw) if read_p else 0
        nodes = io.stream_ints(tw, info[0])
        for n in nodes:
            if n < 0 or n >= n_nodes:
                io.fail(f"element {e} names node {n} of {n_nodes}")
        if n_elem_ints:
            io.stream_ints(tw, n_elem_ints)
        elements.append((info, parent, sid, plev, level, nodes))

    f["coords"] = io.stream_reals(3 * n_nodes)
    if v096 and io.scalar(4) != 0:  # "presence of unique ids"
        io.stream_ints(8 if field_width == 8 else 4, n_nodes)
    if n_node_ints:
        io.stream_ints(tw, n_node_ints * n_nodes)

    if bc_file == "n/a":
        return f

    def triples():
        names = _name_map(io, hw) if v092 else {}
        n = io.scalar(hw)
        values = io.stream_ints(tw, 3 * n)
        return names, [tuple(values[k : k + 3]) for k in range(0, 3 * n, 3)]

    f["sideset_names"], f["sides"] = triples()
    if v092:
        f["nodeset_names"] = _name_map(io, hw)
        n = io.scalar(hw)
        values = io.stream_ints(tw, 2 * n)
        f["nodesets"] = [tuple(values[k : k + 2]) for k in range(0, 2 * n, 2)]
    if v110:
        f["skipped_edge_bcs"] += len(triples()[1])  # edge
        f["skipped_edge_bcs"] += len(triples()[1])  # shell face
    return f


def _on_side(corners, q):
    """Whether point ``q`` lies on the flat side whose corners are ``corners``."""
    p = [np.asarray(c, dtype=float) for c in corners]
    q = np.asarray(q, dtype=float)
    scale = max(np.linalg.norm(p[(k + 1) % len(p)] - p[k]) for k in range(len(p)))
    tol = 1e-6 * (scale if scale > 0 else 1.0)
    if len(p) == 2:
        d = p[1] - p[0]
        rel = q - p[0]
        len2 = float(d @ d)
        if len2 == 0.0:
            return False
        t = float(rel @ d) / len2
        off = np.cross(rel, d)
        return -1e-6 <= t <= 1 + 1e-6 and math.sqrt(float(off @ off) / len2) <= tol

    def in_triangle(a, b, c):
        n = np.cross(b - a, c - a)
        n2 = float(n @ n)
        if n2 == 0.0:
            return False
        if abs(float((q - a) @ n)) / math.sqrt(n2) > tol:
            return False
        eps = -1e-6 * n2
        return (
            float(np.cross(b - a, q - a) @ n) >= eps
            and float(np.cross(c - b, q - b) @ n) >= eps
            and float(np.cross(a - c, q - c) @ n) >= eps
        )

    if len(p) == 3:
        return in_triangle(p[0], p[1], p[2])
    return in_triangle(p[0], p[1], p[2]) or in_triangle(p[0], p[2], p[3])


def read(filename):
    with open_file(filename, "rb") as fh:
        data = fh.read()
    if isinstance(data, str):
        data = data.encode("latin-1")
    # XDR starts with the version string's 4-byte big-endian length.
    xdr = False
    if len(data) >= 8 and not data.startswith(b"libMesh"):
        length = struct.unpack(">I", data[:4])[0]
        xdr = 0 < length < 256 and data[4:11] == b"libMesh"
    f = _parse(data, xdr)

    elements = f["elements"]
    ne = len(elements)
    active = [True] * ne
    children = [[] for _ in range(ne)]
    max_level = 0
    for e, el in enumerate(elements):
        max_level = max(max_level, el[4])
        if el[1] >= 0:
            active[el[1]] = False
            children[el[1]].append(e)

    coords = f["coords"].reshape(-1, 3)
    nn = len(coords)
    used = np.zeros(nn, dtype=bool)
    for el in elements:
        used[el[5]] = True
    kept = np.flatnonzero(used)
    node_index = np.full(nn, -1, dtype=np.int64)
    node_index[kept] = np.arange(len(kept))

    order_of_types = []
    by_type = {}
    skipped = {}
    dropped = {}
    for e, el in enumerate(elements):
        if not active[e]:
            continue
        info = el[0]
        if info[1] is None:
            skipped[info[4]] = skipped.get(info[4], 0) + 1
            continue
        if info[1] not in by_type:
            by_type[info[1]] = []
            order_of_types.append(info[1])
        by_type[info[1]].append(e)
        if info[2] < info[0]:
            dropped[info[4]] = dropped.get(info[4], 0) + 1
    for name in sorted(skipped):
        warn(
            f"libMesh: skipping {skipped[name]} {name} element(s) (no meshio++ "
            "equivalent)"
        )
    for name in sorted(dropped):
        warn(
            f"libMesh: {dropped[name]} {name} element(s) keep only the nodes of the "
            "nearest meshio++ type"
        )
        _provenance.note(
            "high-order-dropped",
            f"{dropped[name]} {name} element(s) lost their extra nodes",
        )

    cells = []
    cell_of = [-1] * ne
    cell_dim = []
    sid_blocks, level_blocks, p_blocks = [], [], []
    for cell_type in order_of_types:
        members = by_type[cell_type]
        order = node_order("libmesh", cell_type)
        k = len(order.to_meshio) if order else elements[members[0]][0][2]
        src = list(order.to_meshio) if order else list(range(k))
        conn = np.empty((len(members), k), dtype=np.int64)
        for r, e in enumerate(members):
            nodes = elements[e][5]
            conn[r] = node_index[[nodes[j] for j in src]]
            cell_of[e] = len(cell_dim)
            cell_dim.append(topological_dimension[cell_type])
        cells.append((cell_type, conn))
        sid_blocks.append(np.array([elements[e][2] for e in members], dtype=np.int64))
        level_blocks.append(np.array([elements[e][4] for e in members], dtype=np.int64))
        p_blocks.append(np.array([elements[e][3] for e in members], dtype=np.int64))

    mesh = Mesh(coords[kept].astype(np.float64), cells)
    if len(kept) != nn:
        mesh.point_data["libmesh:id"] = kept.astype(np.int64)
    if not cells:
        return mesh
    mesh.cell_data["libmesh:subdomain"] = sid_blocks
    if max_level > 0:
        mesh.cell_data["libmesh:level"] = level_blocks
    if f["inline_p"]:
        mesh.cell_data["libmesh:p_level"] = p_blocks

    regions = []
    by_sid = {}
    sid_dim = {}
    for e, el in enumerate(elements):
        if cell_of[e] < 0:
            continue
        by_sid.setdefault(el[2], []).append(cell_of[e])
        sid_dim[el[2]] = max(sid_dim.get(el[2], -1), cell_dim[cell_of[e]])
    for sid in sorted(by_sid):
        name = f["subdomain_names"].get(sid, f"subdomain_{sid}")
        regions.append(
            Region(
                name, "cell", np.array(by_sid[sid], dtype=np.int64), sid_dim[sid], sid
            )
        )

    if f["sides"]:
        facets = FacetIndex(mesh)
        by_id = {}
        id_dim = {}
        unmatched = 0
        point_sides = 0

        def add(bid, e, s):
            nonlocal unmatched
            el = elements[e]
            corners = [node_index[el[5][k]] for k in _SIDES[el[0][3]][s]]
            hit = facets.find(corners)
            cell = cell_of[e]
            owner = None
            if hit is not None and hit.first[0] == cell:
                owner = hit.first
            elif hit is not None and hit.second is not None and hit.second[0] == cell:
                owner = hit.second
            if owner is None:
                unmatched += 1
                return
            by_id.setdefault(bid, set()).add((owner[0], owner[1]))
            id_dim[bid] = max(id_dim.get(bid, -1), cell_dim[cell] - 1)

        for elem, side, bid in f["sides"]:
            if elem < 0 or elem >= ne:
                unmatched += 1
                continue
            el = elements[elem]
            if el[0][3] == "edge":
                point_sides += 1
                continue
            sides = _SIDES.get(el[0][3], [])
            if side < 0 or side >= len(sides) or el[0][1] is None:
                unmatched += 1
                continue
            if active[elem]:
                add(bid, elem, side)
                continue
            side_pts = [coords[el[5][k]] for k in sides[side]]
            stack = list(children[elem])
            while stack:
                c = stack.pop()
                if not active[c]:
                    stack.extend(children[c])
                    continue
                ch = elements[c]
                if ch[0][1] is None:
                    continue
                for cs, local in enumerate(_SIDES.get(ch[0][3], [])):
                    if all(_on_side(side_pts, coords[ch[5][k]]) for k in local):
                        add(bid, c, cs)
        if unmatched:
            warn(
                f"libMesh: {unmatched} boundary side(s) match no cell facet and are "
                "skipped"
            )
        if point_sides:
            warn(
                f"libMesh: {point_sides} boundary side(s) of line elements (end points) "
                "have no side region form and are skipped"
            )
        for bid in sorted(by_id):
            name = f["sideset_names"].get(bid, f"boundary_{bid}")
            entries = np.array(sorted(by_id[bid]), dtype=np.int64).reshape(-1, 2)
            regions.append(Region(name, "side", entries, id_dim[bid], bid))

    by_nodeset = {}
    for node, bid in f["nodesets"]:
        if 0 <= node < nn and node_index[node] >= 0:
            by_nodeset.setdefault(bid, []).append(int(node_index[node]))
    for bid in sorted(by_nodeset):
        name = f["nodeset_names"].get(bid, f"nodeset_{bid}")
        regions.append(
            Region(name, "point", np.array(by_nodeset[bid], dtype=np.int64), -1, bid)
        )
    if f["skipped_edge_bcs"]:
        warn(
            f"libMesh: {f['skipped_edge_bcs']} edge/shell-face boundary condition(s) "
            "skipped"
        )
    mesh.regions = regions
    return mesh

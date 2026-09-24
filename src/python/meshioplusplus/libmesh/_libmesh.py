"""I/O for libMesh ``.xda`` (ASCII) / ``.xdr`` (XDR binary) meshes.

The pure-Python twin of ``src/cpp/src/formats/libmesh.cpp``: both engines read
the same meshes.

Both encodings carry one value stream (libMesh's ``XdrIO``): a version string,
the element and node counts, four "inline or not" flags, per-field integer sizes
(0.9.2+), subdomain names, one connectivity block per refinement level, the
coordinates, then the side sets, node sets (0.9.2+) and edge and shell-face sets
(1.1.0+). ``.xdr`` is big-endian XDR, optionally gzip- or bzip2-compressed.

Only active (leaf) elements become cells. The subdomain id is the
``libmesh:subdomain`` cell data and a cell region per subdomain; side sets
become side regions (carried down to refined children), node sets point
regions, edge sets ``line``/``line3`` cells in a ``<name>:edge`` cell region and
shell-face sets ``<name>:shellface<k>`` cell regions. HEX20/HEX27/PRISM15/PRISM18
use the ``"libmesh"`` node-order tables. ``write`` is the reverse, in the
libMesh-1.8.0 layout.
"""

import bz2
import gzip
import math
import struct

import numpy as np

from .. import _provenance
from .._common import warn
from .._exceptions import ReadError, WriteError
from .._facets import FacetIndex, facet_nodes
from .._files import open_file
from .._mesh import Mesh, topological_dimension
from .._node_order import node_order
from .._regions import Region

__all__ = ["read", "write"]

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

# Corner nodes of each edge, libMesh's edge numbering (Hex8::edge_nodes_map ...).
# A quadratic element's mid-edge node for edge k is node ``vertices + k``.
_EDGES = {
    "tri": [(0, 1), (1, 2), (2, 0)],
    "quad": [(0, 1), (1, 2), (2, 3), (3, 0)],
    "tet": [(0, 1), (1, 2), (0, 2), (0, 3), (1, 3), (2, 3)],
    "hex": [
        (0, 1),
        (1, 2),
        (2, 3),
        (0, 3),
        (0, 4),
        (1, 5),
        (2, 6),
        (3, 7),
        (4, 5),
        (5, 6),
        (6, 7),
        (4, 7),
    ],
    "prism": [(0, 1), (1, 2), (0, 2), (0, 3), (1, 4), (2, 5), (3, 4), (4, 5), (3, 5)],
    "pyramid": [(0, 1), (1, 2), (2, 3), (0, 3), (0, 4), (1, 4), (2, 4), (3, 4)],
}
_VERTICES = {
    "point": 1,
    "edge": 2,
    "tri": 3,
    "quad": 4,
    "tet": 4,
    "pyramid": 5,
    "prism": 6,
    "hex": 8,
}

# Region-name suffixes for the sets that are not sides.
_EDGE_SUFFIX = ":edge"
_SHELLFACE_SUFFIX = ":shellface"


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
        "edges": [],
        "shellfaces": [],
        "nodesets": [],
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

    # Side, edge and shell-face sets share one id space, and each repeats the
    # sideset name map.
    def triples():
        if v092:
            for k, v in _name_map(io, hw).items():
                f["sideset_names"].setdefault(k, v)
        n = io.scalar(hw)
        values = io.stream_ints(tw, 3 * n)
        return [tuple(values[k : k + 3]) for k in range(0, 3 * n, 3)]

    f["sides"] = triples()
    if v092:
        f["nodeset_names"] = _name_map(io, hw)
        n = io.scalar(hw)
        values = io.stream_ints(tw, 2 * n)
        f["nodesets"] = [tuple(values[k : k + 2]) for k in range(0, 2 * n, 2)]
    if v110:
        f["edges"] = triples()
        f["shellfaces"] = triples()
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
    # libMesh writes `.xda.gz`/`.xdr.gz` through gzip and `.bz2` through bzip2.
    if data[:2] == b"\x1f\x8b":
        data = gzip.decompress(data)
    elif data[:3] == b"BZh":
        data = bz2.decompress(data)
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

    # Edge sets -> line cells (the element's edge corners, plus its mid-edge
    # node when it is quadratic), shared by every set naming that edge. They
    # are not libMesh elements: their subdomain (and level, p-level) is -1.
    by_edge_set = {}
    if cells and f["edges"]:
        seen = {}
        rows = ([], [])  # line, line3
        members = []
        unmatched = 0
        for elem, k, bid in f["edges"]:
            if elem < 0 or elem >= ne:
                unmatched += 1
                continue
            el = elements[elem]
            table = _EDGES.get(el[0][3], [])
            if k < 0 or k >= len(table):
                unmatched += 1
                continue
            nv = _VERTICES[el[0][3]]
            row = [int(node_index[el[5][c]]) for c in table[k]]
            if len(el[5]) >= nv + len(table):
                row.append(int(node_index[el[5][nv + k]]))
            key = tuple(sorted(row))
            kind = 1 if len(row) == 3 else 0
            if key not in seen:
                seen[key] = (kind, len(rows[kind]))
                rows[kind].append(row)
            members.append((bid, seen[key]))
        first = [0, 0]
        for kind in (0, 1):
            if not rows[kind]:
                continue
            first[kind] = len(cell_dim)
            n = len(rows[kind])
            cells.append(("line3" if kind else "line", np.array(rows[kind], np.int64)))
            for blocks in (sid_blocks, level_blocks, p_blocks):
                blocks.append(np.full(n, -1, dtype=np.int64))
            cell_dim.extend([1] * n)
        for bid, (kind, row) in members:
            by_edge_set.setdefault(bid, set()).add(first[kind] + row)
        if unmatched:
            warn(
                f"libMesh: {unmatched} edge boundary condition(s) name no element edge "
                "and are skipped"
            )

    mesh = Mesh(coords[kept].astype(np.float64), cells)
    if len(kept) != nn:
        mesh.point_data["libmesh:id"] = kept.astype(np.int64)
    if not cells:
        return mesh
    mesh.cell_data["libmesh:subdomain"] = sid_blocks
    if max_level > 0:
        mesh.cell_data["libmesh:level"] = level_blocks
        # The refinement tree, every element in file order, so the writer can
        # write it back: `libmesh:tree` (cell or -1, parent row or -1, libMesh
        # type, subdomain, p-level) and `libmesh:tree:nodes` (its nodes as
        # point indices, in libMesh's order, -1 past them).
        code_of = {id(t): k for k, t in enumerate(_TYPES)}
        width = max(len(el[5]) for el in elements)
        tree_nodes = np.full((ne, width), -1, dtype=np.int64)
        for e, el in enumerate(elements):
            tree_nodes[e, : len(el[5])] = node_index[list(el[5])]
        mesh.field_data["libmesh:tree"] = np.array(
            [
                [cell_of[e], el[1], code_of[id(el[0])], el[2], el[3]]
                for e, el in enumerate(elements)
            ],
            dtype=np.int64,
        ).reshape(-1, 5)
        mesh.field_data["libmesh:tree:nodes"] = tree_nodes
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

    def set_name(bid):
        return f["sideset_names"].get(bid, f"boundary_{bid}")

    for bid in sorted(by_edge_set):
        entries = np.array(sorted(by_edge_set[bid]), dtype=np.int64)
        regions.append(Region(set_name(bid) + _EDGE_SUFFIX, "cell", entries, 1, bid))

    # Shell-face sets -> cell regions `<name>:shellface<k>` on the 2-D cells (a
    # refined element's face is carried to all its active descendants).
    if f["shellfaces"]:
        by_face = {}
        unmatched = 0
        for elem, face, bid in f["shellfaces"]:
            if elem < 0 or elem >= ne or face not in (0, 1):
                unmatched += 1
                continue
            found = by_face.setdefault((bid, face), set())
            stack = [elem]
            while stack:
                e = stack.pop()
                if not active[e]:
                    stack.extend(children[e])
                    continue
                if cell_of[e] >= 0 and cell_dim[cell_of[e]] == 2:
                    found.add(cell_of[e])
                else:
                    unmatched += 1
        if unmatched:
            warn(
                f"libMesh: {unmatched} shell-face boundary condition(s) name no 2-D "
                "element and are skipped"
            )
        for bid, face in sorted(by_face):
            found = by_face[(bid, face)]
            if found:
                regions.append(
                    Region(
                        f"{set_name(bid)}{_SHELLFACE_SUFFIX}{face}",
                        "cell",
                        np.array(sorted(found), dtype=np.int64),
                        2,
                        bid,
                    )
                )

    by_nodeset = {}
    for node, bid in f["nodesets"]:
        if 0 <= node < nn and node_index[node] >= 0:
            by_nodeset.setdefault(bid, []).append(int(node_index[node]))
    for bid in sorted(by_nodeset):
        name = f["nodeset_names"].get(bid, f"nodeset_{bid}")
        regions.append(
            Region(name, "point", np.array(by_nodeset[bid], dtype=np.int64), -1, bid)
        )
    mesh.regions = regions
    return mesh


# --- writing -------------------------------------------------------------------------

# meshio++ type -> libMesh ElemType, for the types libMesh stores as they are.
_CODES = {
    "line": 0,
    "line3": 1,
    "line4": 2,
    "triangle": 3,
    "triangle6": 4,
    "quad": 5,
    "quad8": 6,
    "quad9": 7,
    "tetra": 8,
    "tetra10": 9,
    "hexahedron": 10,
    "hexahedron20": 11,
    "hexahedron27": 12,
    "wedge": 13,
    "wedge15": 14,
    "wedge18": 15,
    "pyramid": 16,
    "pyramid13": 17,
    "pyramid14": 18,
    "vertex": 27,
    "triangle7": 33,
}


class _Out:
    """libMesh's ``Xdr`` in WRITE (ASCII) or ENCODE (XDR) mode, ``_Stream`` reversed."""

    def __init__(self, xdr):
        self.xdr = xdr
        self.out = bytearray() if xdr else []

    def _int(self, v, width):
        self.out += struct.pack(">Q" if width == 8 else ">I", v % (1 << (8 * width)))

    def _str(self, s):
        b = s.encode()
        self.out += struct.pack(">I", len(b)) + b + b"\0" * ((4 - len(b) % 4) % 4)

    def _comment(self, comment):
        self.out.append(f"\t {comment}\n" if comment else "\n")

    def string(self, s, comment=""):
        if self.xdr:
            self._str(s)
            return
        self.out.append(s)
        self._comment(comment)

    def scalar(self, v, comment, width=8):
        if self.xdr:
            self._int(v, width)
            return
        self.out.append(str(v))
        self._comment(comment)

    def int_vector(self, values, comment=""):
        if self.xdr:
            self._int(len(values), 4)
            for v in values:
                self._int(v, 8)
            return
        self.scalar(len(values), "# vector length")
        self.out.append("".join(f"{v}\t " for v in values))
        self._comment(comment)

    def string_vector(self, values, comment=""):
        if self.xdr:
            self._int(len(values), 4)
            for v in values:
                self._str(v)
            return
        self.scalar(len(values), "# vector length")
        self.out.append("".join(f"{v}\t " for v in values))
        self._comment(comment)

    def ints(self, values):
        if self.xdr:
            for v in values:
                self._int(v, 8)
            return
        self.out.append(" ".join(str(v) for v in values) + "\n")

    def reals(self, values):
        if self.xdr:
            self.out += np.asarray(values, dtype=">f8").tobytes()
            return
        parts = [f"{v:.17e}" for v in values.tolist()]
        for k in range(0, len(parts), 3):
            self.out.append(" ".join(parts[k : k + 3]) + "\n")

    def data(self):
        return bytes(self.out) if self.xdr else "".join(self.out).encode()


def _tree(mesh, elems, elem_of, ncells, npts):
    """The refinement tree (``libmesh:tree``, as the reader keeps it) when it
    still describes these cells: each written cell is exactly one leaf of the
    tree, with the same type and nodes. Returns (rows, row of each written
    element), a row being (cell, parent, code, subdomain, p-level, level,
    nodes); None (with a warning) when it does not match, or is absent."""
    t = mesh.field_data.get("libmesh:tree")
    tn = mesh.field_data.get("libmesh:tree:nodes")
    if t is None or tn is None:
        return None
    t = np.asarray(t)
    tn = np.asarray(tn)
    ok = (
        t.ndim == 2
        and t.shape[1] == 5
        and len(t) > 0
        and tn.ndim == 2
        and len(tn) == len(t)
        and tn.shape[1] > 0
    )
    rows = []
    has_child = [False] * (len(t) if ok else 0)
    last_level = 0
    for r in range(len(t) if ok else 0):
        cell, parent, code, sd, plev = (int(x) for x in t[r])
        info = _TYPES[code] if 0 <= code < len(_TYPES) else None
        nodes = []
        for v in tn[r].tolist():
            if v < 0:
                break
            nodes.append(int(v))
        ok = (
            info is not None
            and info[0] > 0
            and -1 <= parent < r
            and 0 <= sd <= 65534
            and len(nodes) == info[0]
            and all(v < npts for v in nodes)
        )
        if not ok:
            break
        level = 0
        if parent >= 0:
            level = rows[parent][5] + 1
            has_child[parent] = True
        ok = level >= last_level
        last_level = level
        rows.append((cell, parent, code, sd, plev, level, nodes))
    named = set()
    row_of_elem = [-1] * len(elems)
    for r in range(len(rows) if ok else 0):
        cell, _, code, _, _, _, nodes = rows[r]
        if has_child[r]:
            ok = cell == -1
        else:
            ok = 0 <= cell < ncells and elem_of[cell] >= 0 and cell not in named
            if ok:
                named.add(cell)
                i = int(elem_of[cell])
                ok = elems[i][0] == code and list(elems[i][2]) == nodes
                row_of_elem[i] = r
        if not ok:
            break
    ok = ok and all(r >= 0 for r in row_of_elem)
    if not ok:
        warn(
            "libMesh: `libmesh:tree` no longer matches the cells; the mesh is written "
            "flat, as its active cells"
        )
        _provenance.note(
            "refinement-tree-dropped", "libmesh:tree does not match the cells"
        )
        return None
    return rows, row_of_elem


def _lift(tree, points, sides, edges, shellfaces):
    """libMesh keeps boundary ids on level-0 elements: a set's entries are
    lifted to each level-0 side (edge, shell face) whose active pieces are all
    in the set; the rest are dropped with a warning."""
    rows, row_of_elem = tree
    nt = len(rows)
    kids = [[] for _ in range(nt)]
    leaf = [True] * nt
    for r, row in enumerate(rows):
        if row[1] >= 0:
            kids[row[1]].append(r)
            leaf[row[1]] = False
    pts = np.zeros((len(points), 3))
    pts[:, : min(points.shape[1], 3)] = points[:, :3]

    def leaves(root):
        out, stack = [], [root]
        while stack:
            r = stack.pop()
            if leaf[r]:
                out.append(r)
            else:
                stack.extend(reversed(kids[r]))
        return out

    def table_of(r, use_edges):
        shape = _TYPES[rows[r][2]][3]
        return (
            [list(e) for e in _EDGES.get(shape, [])]
            if use_edges
            else _SIDES.get(shape, [])
        )

    def pieces_on(root, on, use_edges):
        out = []
        for c in leaves(root):
            for k, local in enumerate(table_of(c, use_edges)):
                if all(_on_side(on, pts[rows[c][6][v]]) for v in local):
                    out.append((c, k))
        return out

    roots = [r for r in range(nt) if rows[r][1] < 0]
    lost = 0

    def lift(values, use_edges):
        nonlocal lost
        want = {}
        for e, k, bid in values:
            want.setdefault(bid, set()).add((row_of_elem[e], k))
        out = set()
        for bid in sorted(want):
            pieces = want[bid]
            covered = set()
            for r0 in roots:
                for k, local in enumerate(table_of(r0, use_edges)):
                    on = [pts[rows[r0][6][v]] for v in local]
                    d = pieces_on(r0, on, use_edges)
                    if d and all(x in pieces for x in d):
                        out.add((r0, k, bid))
                        covered.update(d)
            lost += sum(1 for x in pieces if x not in covered)
        return out

    sides = lift(sides, False)
    edges = lift(edges, True)
    want = {}
    for e, face, bid in shellfaces:
        want.setdefault((bid, face), set()).add(row_of_elem[e])
    lifted = set()
    for bid, face in sorted(want):
        cells = want[(bid, face)]
        covered = set()
        for r0 in roots:
            under = leaves(r0)
            if all(c in cells for c in under):
                lifted.add((r0, face, bid))
                covered.update(under)
        lost += sum(1 for c in cells if c not in covered)
    if lost:
        warn(
            f"libMesh: {lost} boundary set entries cover only part of a level-0 "
            "element's side, edge or face and are dropped (libMesh keeps boundary ids "
            "on level-0 elements)"
        )
        _provenance.note(
            "regions-dropped", f"{lost} boundary entries are not whole level-0 sides"
        )
    return sides, edges, lifted


def _shellface(name):
    """A shell-face region's ``(face, base name)``, else ``None``."""
    for k in (0, 1):
        suffix = f"{_SHELLFACE_SUFFIX}{k}"
        if len(name) > len(suffix) and name.endswith(suffix):
            return k, name[: -len(suffix)]
    return None


def write(filename, mesh):
    """Write a libMesh ``.xda`` (ASCII) or ``.xdr`` (XDR) mesh, libMesh-1.8.0.

    The encoding comes from the extension; a trailing ``.gz``/``.bz2``
    compresses an ASCII file (a ``.xdr.gz``/``.xdr.bz2`` is plain XDR, as
    libMesh writes it). See ``doc/formats/libmesh.md``.
    """
    name = str(filename).lower()
    compress = None
    for suffix, module in ((".gz", gzip), (".bz2", bz2)):
        if name.endswith(suffix):
            compress = module
            name = name[: -len(suffix)]
    xdr = name.endswith(".xdr")
    # libMesh's XDR files ignore the suffix (its `Xdr` opens them with plain
    # stdio): a `.xdr.gz`/`.xdr.bz2` is plain XDR, which libMesh then reads.
    if xdr:
        compress = None

    blocks = list(mesh.cells)
    starts = [0]
    block_dim = []
    for block in blocks:
        starts.append(starts[-1] + len(block.data))
        block_dim.append(topological_dimension.get(block.type, -1))
    ncells = starts[-1]
    cell_block = np.repeat(np.arange(len(blocks)), np.diff(starts)).astype(np.int64)

    def dim_of(cell):
        return block_dim[cell_block[cell]]

    # Region roles: edge sets and shell faces (by name suffix), side sets, node
    # sets, and the remaining cell regions as subdomains.
    regions = sorted(getattr(mesh, "regions", []) or [], key=lambda r: r.key)
    edge_sets, shell_sets, side_sets, node_sets, subdomains = [], [], [], [], []
    placeholder = np.zeros(ncells, dtype=bool)
    for reg in regions:
        entries = np.asarray(reg.entries, dtype=np.int64)
        if reg.kind == "side":
            side_sets.append(reg)
            continue
        if reg.kind == "point":
            node_sets.append(reg)
            continue
        cells = entries.ravel()
        if not np.all((cells >= 0) & (cells < ncells)):
            continue
        dims = {dim_of(c) for c in cells}
        if len(cells) and reg.name.endswith(_EDGE_SUFFIX) and dims == {1}:
            if len(reg.name) > len(_EDGE_SUFFIX):
                edge_sets.append(reg)
                placeholder[cells] = True
                continue
        if len(cells) and _shellface(reg.name) and dims == {2}:
            shell_sets.append(reg)
            continue
        subdomains.append(reg)

    # Elements: every cell libMesh has a type for, in cell order.
    elems = []  # (code, cell, nodes in libMesh order)
    elem_of = np.full(ncells, -1, dtype=np.int64)
    dropped = {}
    for b, block in enumerate(blocks):
        code = _CODES.get(block.type) if isinstance(block.data, np.ndarray) else None
        if code is None:
            if len(block.data):
                dropped[block.type] = dropped.get(block.type, 0) + len(block.data)
            continue
        conn = np.asarray(block.data, dtype=np.int64)
        order = node_order("libmesh", block.type)
        if order is not None:
            conn = conn[:, list(order.from_meshio)]
        for r, row in enumerate(conn.tolist()):
            g = starts[b] + r
            if placeholder[g]:
                continue
            elem_of[g] = len(elems)
            elems.append((code, g, row))
    for t in sorted(dropped):
        warn(
            f"libMesh: {dropped[t]} '{t}' cell(s) have no libMesh element type and are "
            "dropped"
        )
        _provenance.note(
            "cells-dropped", f"{dropped[t]} '{t}' cell(s) have no libMesh element type"
        )

    # Node ids: `libmesh:id` when it is a valid numbering, else the point index.
    points = np.asarray(mesh.points, dtype=np.float64)
    if points.ndim == 1:
        points = points.reshape(-1, 1)
    npts = len(points)
    node_id = np.arange(npts, dtype=np.int64)
    max_node_id = npts
    ids = mesh.point_data.get("libmesh:id")
    if ids is not None:
        ids = np.asarray(ids).ravel()
        if (
            len(ids) == npts
            and np.all(ids >= 0)
            and len(np.unique(ids)) == npts
            and np.issubdtype(ids.dtype, np.integer)
        ):
            node_id = ids.astype(np.int64)
            max_node_id = int(ids.max()) + 1 if npts else 0

    # Subdomain ids: `libmesh:subdomain`, else the first cell region holding
    # the cell (its tag, or a fresh id), else 0.
    sid = np.zeros(ncells, dtype=np.int64)
    subdomain_names = {}
    if "libmesh:subdomain" in mesh.cell_data:
        for b, arr in enumerate(mesh.cell_data["libmesh:subdomain"]):
            sid[starts[b] : starts[b + 1]] = np.asarray(arr).ravel().astype(np.int64)
        for reg in subdomains:
            if reg.tag >= 0 and reg.name != f"subdomain_{reg.tag}":
                subdomain_names.setdefault(reg.tag, reg.name)
    else:
        nxt = max([0] + [reg.tag + 1 for reg in subdomains])
        assigned = np.zeros(ncells, dtype=bool)
        for reg in subdomains:
            if reg.tag >= 0:
                sd = reg.tag
            else:
                sd = nxt
                nxt += 1
            for c in np.asarray(reg.entries, dtype=np.int64).ravel():
                if not assigned[c]:
                    assigned[c] = True
                    sid[c] = sd
            if reg.name != f"subdomain_{sd}":
                subdomain_names.setdefault(sd, reg.name)
    for _, cell, _ in elems:
        if not 0 <= sid[cell] <= 65534:
            raise WriteError(
                f"libMesh: subdomain id {sid[cell]} is outside libMesh's 0..65534"
            )
    p_level = mesh.cell_data.get("libmesh:p_level")
    write_p = p_level is not None

    # Boundary ids: one id space for side, edge and shell-face sets.
    next_bid = max([0] + [r.tag + 1 for r in side_sets + edge_sets + shell_sets])
    bids = {}
    sideset_names = {}

    def boundary_id(reg, base):
        nonlocal next_bid
        if id(reg) in bids:
            return bids[id(reg)]
        if reg.tag >= 0:
            bid = reg.tag
        else:
            bid = next_bid
            next_bid += 1
        bids[id(reg)] = bid
        if base != f"boundary_{bid}":
            sideset_names.setdefault(bid, base)
        return bid

    # Side sets: (element, libMesh side, id), matched by the facet's corners.
    sides = set()
    sides_lost = 0
    for reg in side_sets:
        bid = boundary_id(reg, reg.name)
        for cell, facet in np.asarray(reg.entries, dtype=np.int64).reshape(-1, 2):
            hit = None
            if 0 <= cell < ncells and elem_of[cell] >= 0:
                hit = facet_nodes(mesh, int(cell), int(facet))
            if hit is None:
                sides_lost += 1
                continue
            ftype, fnodes = hit
            corners = (
                2
                if topological_dimension.get(ftype) == 1
                else (3 if ftype.startswith("triangle") else 4)
            )
            key = sorted(fnodes[:corners])
            e = int(elem_of[cell])
            code, _, nodes = elems[e]
            for s, local in enumerate(_SIDES.get(_TYPES[code][3], [])):
                if sorted(nodes[k] for k in local) == key:
                    sides.add((e, s, bid))
                    break
            else:
                sides_lost += 1

    # Edge sets: each line cell on the first element holding that edge.
    edges = set()
    edges_lost = 0
    # Edge cells no active element holds (a refined element's whole edge):
    # (corner, corner, id), matched against the tree's level-0 elements.
    edges_unmatched = []
    if edge_sets:
        owner = {}
        for e, (code, _, nodes) in enumerate(elems):
            for k, (a, b) in enumerate(_EDGES.get(_TYPES[code][3], [])):
                owner.setdefault(tuple(sorted((nodes[a], nodes[b]))), (e, k))
        for reg in edge_sets:
            bid = boundary_id(reg, reg.name[: -len(_EDGE_SUFFIX)])
            for cell in np.asarray(reg.entries, dtype=np.int64).ravel():
                b = cell_block[cell]
                row = blocks[b].data[cell - starts[b]]
                key = tuple(sorted((int(row[0]), int(row[1]))))
                hit = owner.get(key)
                if hit is None:
                    edges_unmatched.append((key, bid))
                    continue
                edges.add((hit[0], hit[1], bid))

    # Shell faces.
    shellfaces = set()
    for reg in shell_sets:
        face, base = _shellface(reg.name)
        bid = boundary_id(reg, base)
        for cell in np.asarray(reg.entries, dtype=np.int64).ravel():
            if elem_of[cell] >= 0:
                shellfaces.add((int(elem_of[cell]), face, bid))

    # Node sets.
    next_nid = max([0] + [r.tag + 1 for r in node_sets])
    nodesets = set()
    nodeset_names = {}
    for reg in node_sets:
        if reg.tag >= 0:
            nid = reg.tag
        else:
            nid = next_nid
            next_nid += 1
        if reg.name != f"nodeset_{nid}":
            nodeset_names.setdefault(nid, reg.name)
        for p in np.asarray(reg.entries, dtype=np.int64).ravel():
            if 0 <= p < npts:
                nodesets.add((int(node_id[p]), nid))
    tree = _tree(mesh, elems, elem_of, ncells, npts)
    if tree is not None:
        sides, edges, shellfaces = _lift(tree, points, sides, edges, shellfaces)
        # Whole edges of refined level-0 elements.
        root_edge = {}
        for r0, row in enumerate(tree[0]):
            if row[1] >= 0:
                break
            for k, (a, b) in enumerate(_EDGES.get(_TYPES[row[2]][3], [])):
                root_edge.setdefault(tuple(sorted((row[6][a], row[6][b]))), (r0, k))
        still = []
        for key, bid in edges_unmatched:
            hit = root_edge.get(key)
            if hit is None:
                still.append((key, bid))
            else:
                edges.add((hit[0], hit[1], bid))
        edges_unmatched = still
    edges_lost += len(edges_unmatched)
    if sides_lost or edges_lost:
        warn(
            f"libMesh: {sides_lost} side and {edges_lost} edge set entries match no "
            "written element and are dropped"
        )
        _provenance.note(
            "regions-dropped",
            f"{sides_lost + edges_lost} side/edge set entries match no libMesh element",
        )
    bcs = bool(sides or edges or shellfaces or nodesets)

    # The stream, as XdrIO::write lays it out (libMesh-1.8.0, 8-byte ids).
    io = _Out(xdr)
    io.string("libMesh-1.8.0")
    io.scalar(len(tree[0]) if tree is not None else len(elems), "# number of elements")
    io.scalar(max_node_id, "# number of nodes")
    io.string("." if bcs else "n/a", "# boundary condition specification file")
    io.string(".", "# subdomain id specification file")
    io.string("n/a", "# processor id specification file")
    io.string("." if write_p else "n/a", "# p-level specification file")
    io.scalar(8, "# type size")
    io.scalar(0, "# uid size")
    io.scalar(0, "# pid size")
    io.scalar(8, "# sid size")
    io.scalar(8 if write_p else 0, "# p-level size")
    for label in ("eid", "side", "bid"):
        io.scalar(8 if bcs else 0, f"# {label} size")
    io.scalar(0, "# extra integer size")
    io.string_vector([], "# node integer names")
    io.string_vector([], "# elem integer names")
    io.int_vector([], "# elemset codes")

    def name_map(names, comment):
        io.scalar(len(names), comment)
        if names:
            io.int_vector(sorted(names))
            io.string_vector([names[k] for k in sorted(names)])

    name_map(subdomain_names, "# subdomain id to name map")
    p_flat = None
    if write_p:
        p_flat = np.concatenate(
            [np.asarray(a).ravel().astype(np.int64) for a in p_level]
        )
    legend = "p_level " if write_p else ""
    if tree is not None:
        # Level by level, each element after its parent; the leaves' subdomain
        # and p-level from the cells, their ancestors' from the tree.
        rows = tree[0]
        for r, (cell, parent, code, sd, plev, level, nodes) in enumerate(rows):
            if r == 0 or level != rows[r - 1][5]:
                n = sum(1 for row in rows[r:] if row[5] == level)
                io.scalar(
                    n,
                    f"# n_elem at level {level}, [ type "
                    + ("parent " if level else "")
                    + f"sid {legend}(n0 ... nN-1) ]",
                )
            rec = [code] + ([parent] if level else [])
            if cell >= 0:
                rec.append(int(sid[cell]))
                if write_p:
                    rec.append(int(p_flat[cell]))
            else:
                rec.append(sd)
                if write_p:
                    rec.append(plev)
            rec += node_id[nodes].tolist()
            io.ints(rec)
    elif elems:
        io.scalar(
            len(elems), f"# n_elem at level 0, [ type sid {legend}(n0 ... nN-1) ]"
        )
    for code, cell, nodes in elems if tree is None else []:
        rec = [code, int(sid[cell])]
        if write_p:
            rec.append(int(p_flat[cell]))
        rec += node_id[nodes].tolist()
        io.ints(rec)

    # Unused node ids: NaN in XDR, as libMesh writes them; 0 in ASCII, where
    # libMesh's `>>` cannot read back the `nan` it writes (it only loads nodes
    # elements use, so the value is never looked at).
    coords = np.full((max_node_id, 3), math.nan if xdr else 0.0)
    pd = min(points.shape[1], 3)
    coords[node_id, :] = 0.0
    coords[node_id, :pd] = points[:, :pd]
    io.reals(coords.ravel())
    io.scalar(0, "# presence of unique ids", width=4)

    def triples(values, comment):
        name_map(sideset_names, "# sideset id to name map")
        io.scalar(len(values), comment)
        for t in sorted(values):
            io.ints(list(t))

    triples(sides, "# number of side boundary conditions")
    name_map(nodeset_names, "# nodeset id to name map")
    io.scalar(len(nodesets), "# number of nodesets")
    for t in sorted(nodesets):
        io.ints(list(t))
    triples(edges, "# number of edge boundary conditions")
    triples(shellfaces, "# number of shellface boundary conditions")

    data = io.data()
    if compress is not None:
        data = gzip.compress(data, mtime=0) if compress is gzip else bz2.compress(data)
    with open_file(filename, "wb") as fh:
        fh.write(data)

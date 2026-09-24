"""FEBio plot files (``.xplt``): the pure-Python reference reader.

The twin of ``src/cpp/src/formats/xplt.cpp``. An ``.xplt`` file is a tree of
tagged chunks (``u32 id, u32 size, payload``) after a 4-byte magic: a root with
the header and the dictionary of variables, the mesh, then one state per saved
time. With compression on, every chunk after the first mesh is its own zlib
stream. The layout follows FEBio's writer (``FEBioPlot/FEBioPlotFile.cpp``) and
FEBio Studio's reader (``XPLTLib/xpltReader3.cpp``); see
``doc/formats/xplt.md``.
"""

from __future__ import annotations

import struct
import zlib

import numpy as np

from .._common import warn
from .._exceptions import ReadError
from .._facets import FacetIndex
from .._mesh import CellBlock, Mesh
from .._node_order import node_order
from .._regions import Region

MAGIC = 0x00464542

ROOT = 0x01000000
HEADER = 0x01010000
HDR_VERSION = 0x01010001
HDR_COMPRESSION = 0x01010004
DICTIONARY = 0x01020000
DIC_ITEM = 0x01020001
DIC_ITEM_TYPE = 0x01020002
DIC_ITEM_FMT = 0x01020003
DIC_ITEM_NAME = 0x01020004
DIC_ITEM_ARRAYSIZE = 0x01020005
DIC_GROUPS = {
    0x01021000: "global",
    0x01023000: "node",
    0x01024000: "domain",
    0x01025000: "surface",
    0x01026000: "edge",
}
MESH = 0x01040000
NODE_SECTION = 0x01041000
NODE_HEADER = 0x01041100
NODE_SIZE = 0x01041101
NODE_COORDS = 0x01041200
DOMAIN_SECTION = 0x01042000
DOMAIN = 0x01042100
DOMAIN_HDR = 0x01042101
DOM_ELEM_TYPE = 0x01042102
DOM_PART_ID = 0x01042103
DOM_NAME = 0x01032105
DOM_ELEM_LIST = 0x01042200
ELEMENT = 0x01042201
SURFACE_SECTION = 0x01043000
SURFACE = 0x01043100
SURFACE_HDR = 0x01043101
SURFACE_NAME = 0x01043104
FACE_LIST = 0x01043200
FACE = 0x01043201
NODESET_SECTION = 0x01044000
NODESET = 0x01044100
NODESET_HDR = 0x01044101
NODESET_NAME = 0x01044103
NODESET_LIST = 0x01044200
PARTS_SECTION = 0x01045000
PART = 0x01045100
PART_ID = 0x01045101
PART_NAME = 0x01045102
ELEMENTSET_SECTION = 0x01046000
ELEMENTSET = 0x01046100
ELEMENTSET_HDR = 0x01046101
ELEMENTSET_NAME = 0x01046103
ELEMENTSET_LIST = 0x01046200
FACETSET_SECTION = 0x01047000
FACETSET = 0x01047100
FACETSET_HDR = 0x01047101
FACETSET_NAME = 0x01047103
FACETSET_LIST = 0x01047200
FACET = 0x01047201
EDGE_SECTION = 0x01048000
EDGE = 0x01048100
EDGE_HDR = 0x01048101
EDGE_NAME = 0x01048104
EDGE_LIST = 0x01048200
LINE = 0x01048201
STATE = 0x02000000
STATE_HEADER = 0x02010000
STATE_TIME = 0x02010002
STATE_STATUS = 0x02010003
STATE_DATA = 0x02020000
STATE_VARIABLE = 0x02020001
STATE_VAR_ID = 0x02020002
STATE_VAR_DATA = 0x02020003
DATA_GROUPS = {
    0x02020100: "global",
    0x02020300: "node",
    0x02020400: "domain",
    0x02020500: "surface",
    0x02020600: "edge",
}

# PLT element type -> (meshio++ type, nodes in the file); FEBioPlotFile.h's Elem_Type.
_ELEM = {
    0: ("hexahedron", 8),
    1: ("wedge", 6),
    2: ("tetra", 4),
    3: ("quad", 4),
    4: ("triangle", 3),
    5: ("line", 2),
    6: ("hexahedron20", 20),
    7: ("tetra10", 10),
    8: ("tet15", 15),
    9: ("hexahedron27", 27),
    10: ("triangle6", 6),
    11: ("quad8", 8),
    12: ("quad9", 9),
    13: ("wedge15", 15),
    14: ("tet20", 20),
    15: ("triangle10", 10),
    16: ("pyramid", 5),
    17: ("tet5", 5),
    18: ("pyramid13", 13),
    19: ("line3", 3),
}
# FEBio types with no meshio++ cell, and what lenient keeps of them.
_LOSSY = {"tet15": "tetra10", "tet5": "tetra", "tet20": None}
_NODES = {
    "hexahedron": 8,
    "wedge": 6,
    "tetra": 4,
    "quad": 4,
    "triangle": 3,
    "line": 2,
    "hexahedron20": 20,
    "tetra10": 10,
    "hexahedron27": 27,
    "triangle6": 6,
    "quad8": 8,
    "quad9": 9,
    "wedge15": 15,
    "triangle10": 10,
    "pyramid": 5,
    "pyramid13": 13,
    "line3": 3,
}
_CORNERS = {3: 3, 4: 4, 6: 3, 7: 3, 8: 4, 9: 4, 10: 3}
_FACET_TYPES = {3: "triangle", 4: "quad", 6: "triangle6", 7: "triangle7", 8: "quad8"}
_FACET_TYPES[9] = "quad9"
_FACET_TYPES[10] = "triangle10"
# Var_Type -> float count per value (ARRAY types multiply the array size).
_WIDTH = {0: 1, 1: 3, 2: 6, 3: 3, 4: 21, 5: 9}
_FMT = {0: "node", 1: "item", 2: "mult", 3: "region", 4: "matpoints"}


def _fail(what):
    raise ReadError(f"FEBio .xplt: {what}")


class _Reader:
    def __init__(self, data, swap):
        self.data = data
        self.swap = swap
        self.end = ">" if swap else "<"

    def u32(self, pos):
        return struct.unpack_from(self.end + "I", self.data, pos)[0]

    def i32s(self, start, stop):
        return np.frombuffer(self.data[start:stop], dtype=np.dtype(self.end + "i4"))

    def f32s(self, start, stop):
        return np.frombuffer(self.data[start:stop], dtype=np.dtype(self.end + "f4"))

    def chunks(self, start, stop):
        """(id, payload start, payload end) of every chunk in [start, stop)."""
        pos = start
        while pos + 8 <= stop:
            ident, size = self.u32(pos), self.u32(pos + 4)
            if pos + 8 + size > stop:
                _fail("a chunk runs past its parent (truncated file?)")
            yield ident, pos + 8, pos + 8 + size
            pos += 8 + size

    def child(self, start, stop, ident):
        for i, a, b in self.chunks(start, stop):
            if i == ident:
                return a, b
        return None

    def int_of(self, start, stop, ident, default=None):
        found = self.child(start, stop, ident)
        return (
            default
            if found is None
            else struct.unpack_from(self.end + "i", self.data, found[0])[0]
        )

    def string(self, start, stop):
        """A length-prefixed string (i32 length, then bytes)."""
        if stop - start < 4:
            return ""
        n = struct.unpack_from(self.end + "i", self.data, start)[0]
        return (
            bytes(self.data[start + 4 : start + 4 + n]).decode("latin-1").rstrip("\0")
        )

    def fixed(self, start, stop):
        """A 64-byte NUL-padded string."""
        return bytes(self.data[start:stop]).split(b"\0", 1)[0].decode("latin-1")


def _top_level(data):
    """The top-level chunks as (reader, id, start, stop), inflating compressed ones."""
    if len(data) < 4:
        _fail("file too short")
    magic = struct.unpack_from("<I", data, 0)[0]
    if magic == MAGIC:
        swap = False
    elif struct.unpack_from(">I", data, 0)[0] == MAGIC:
        swap = True
    else:
        _fail("not an FEBio plot file (bad magic)")
    raw = _Reader(data, swap)
    out = []
    pos = 4
    compressed = False
    while pos < len(data):
        if len(out) >= 2 and compressed:
            stream = zlib.decompressobj()
            try:
                chunk = stream.decompress(data[pos:])
            except zlib.error:
                chunk = b""
            if not stream.eof:
                warn("FEBio .xplt: the last state is truncated and was not read")
                break
            pos = len(data) - len(stream.unused_data)
            reader = _Reader(chunk, swap)
            if len(chunk) < 8:
                break
            ident, size = reader.u32(0), reader.u32(4)
            out.append((reader, ident, 8, min(8 + size, len(chunk))))
            continue
        if pos + 8 > len(data):
            warn("FEBio .xplt: trailing bytes after the last state were ignored")
            break
        ident, size = raw.u32(pos), raw.u32(pos + 4)
        if pos + 8 + size > len(data):
            if ident == STATE:
                warn("FEBio .xplt: the last state is truncated and was not read")
                break
            _fail("the file is truncated")
        out.append((raw, ident, pos + 8, pos + 8 + size))
        if ident == ROOT:
            header = raw.child(pos + 8, pos + 8 + size, HEADER)
            if header is not None:
                compressed = raw.int_of(*header, HDR_COMPRESSION, 0) != 0
        pos += 8 + size
    return out


def _dictionary(reader, start, stop):
    items = {g: [] for g in DIC_GROUPS.values()}
    for gid, a, b in reader.chunks(start, stop):
        group = DIC_GROUPS.get(gid)
        if group is None:
            continue
        for iid, c, d in reader.chunks(a, b):
            if iid != DIC_ITEM:
                continue
            vtype = reader.int_of(c, d, DIC_ITEM_TYPE, 0)
            fmt = reader.int_of(c, d, DIC_ITEM_FMT, 0)
            asize = reader.int_of(c, d, DIC_ITEM_ARRAYSIZE, 0)
            found = reader.child(c, d, DIC_ITEM_NAME)
            name = reader.fixed(*found) if found else ""
            if "=" in name:
                name = name.split("=", 1)[1]
            if vtype == 6:
                width = max(asize, 1)
            elif vtype == 7:
                width = 3 * max(asize, 1)
            else:
                width = _WIDTH.get(vtype, 1)
            items[group].append((name, _FMT.get(fmt, "item"), width))
    return items


def _facets(reader, a, b, list_tag, leaf_tag):
    """Facets of a surface or facet set: (node count, 0-based nodes)."""
    out = []
    found = reader.child(a, b, list_tag)
    if found is None:
        return out
    for fid, c, d in reader.chunks(*found):
        if fid != leaf_tag:
            continue
        values = reader.i32s(c, d)
        nn = int(values[1])
        out.append((nn, [int(v) for v in values[2 : 2 + nn]]))
    return out


def _mesh_section(reader, start, stop, version):
    mesh = {
        "nodes": None,
        "domains": [],
        "surfaces": [],
        "edges": [],
        "nodesets": [],
        "elemsets": [],
    }
    mesh["parts"] = {}
    for sid, a, b in reader.chunks(start, stop):
        if sid == NODE_SECTION:
            head = reader.child(a, b, NODE_HEADER)
            n = reader.int_of(*head, NODE_SIZE, 0) if head else 0
            coords = reader.child(a, b, NODE_COORDS)
            if coords is None:
                _fail("the mesh has no node coordinates")
            words = (coords[1] - coords[0]) // 4
            if words == 4 * n:
                xyz = reader.f32s(*coords).reshape(n, 4)[:, 1:]
            elif words == 3 * n:
                xyz = reader.f32s(*coords).reshape(n, 3)
            else:
                _fail("the node coordinates do not match the node count")
            mesh["nodes"] = xyz.astype(np.float64)
        elif sid == DOMAIN_SECTION:
            for did, c, d in reader.chunks(a, b):
                if did != DOMAIN:
                    continue
                hdr = reader.child(c, d, DOMAIN_HDR)
                etype = reader.int_of(*hdr, DOM_ELEM_TYPE, -1)
                part = reader.int_of(*hdr, DOM_PART_ID, -1)
                named = reader.child(*hdr, DOM_NAME)
                name = reader.string(*named) if named else ""
                if etype not in _ELEM:
                    _fail(f"unknown element type {etype}")
                cell_type, nn = _ELEM[etype]
                ids, conn = [], []
                listed = reader.child(c, d, DOM_ELEM_LIST)
                for eid, e, f in reader.chunks(*listed) if listed else []:
                    if eid != ELEMENT:
                        continue
                    values = reader.i32s(e, f)
                    ids.append(int(values[0]))
                    conn.append(values[1 : 1 + nn])
                conn = np.array(conn, dtype=np.int64).reshape(-1, nn)
                mesh["domains"].append((cell_type, part, name, ids, conn))
        elif sid == SURFACE_SECTION:
            for sid2, c, d in reader.chunks(a, b):
                if sid2 != SURFACE:
                    continue
                hdr = reader.child(c, d, SURFACE_HDR)
                named = reader.child(*hdr, SURFACE_NAME) if hdr else None
                name = reader.string(*named) if named else ""
                mesh["surfaces"].append(
                    (name, _facets(reader, c, d, FACE_LIST, FACE), False)
                )
        elif sid == FACETSET_SECTION:
            for sid2, c, d in reader.chunks(a, b):
                if sid2 != FACETSET:
                    continue
                hdr = reader.child(c, d, FACETSET_HDR)
                named = reader.child(*hdr, FACETSET_NAME) if hdr else None
                name = reader.string(*named) if named else ""
                mesh["surfaces"].append(
                    (name, _facets(reader, c, d, FACETSET_LIST, FACET), True)
                )
        elif sid == EDGE_SECTION:
            # an edge (FEBio 3.5+): named segments [edge id, nodes, n1 .. n3]
            for sid2, c, d in reader.chunks(a, b):
                if sid2 != EDGE:
                    continue
                hdr = reader.child(c, d, EDGE_HDR)
                named = reader.child(*hdr, EDGE_NAME) if hdr else None
                name = reader.string(*named) if named else ""
                mesh["edges"].append((name, _facets(reader, c, d, EDGE_LIST, LINE)))
        elif sid in (NODESET_SECTION, ELEMENTSET_SECTION):
            item, hdr_tag, name_tag, list_tag = (
                (NODESET, NODESET_HDR, NODESET_NAME, NODESET_LIST)
                if sid == NODESET_SECTION
                else (ELEMENTSET, ELEMENTSET_HDR, ELEMENTSET_NAME, ELEMENTSET_LIST)
            )
            key = "nodesets" if sid == NODESET_SECTION else "elemsets"
            for sid2, c, d in reader.chunks(a, b):
                if sid2 != item:
                    continue
                hdr = reader.child(c, d, hdr_tag)
                named = reader.child(*hdr, name_tag) if hdr else None
                name = reader.string(*named) if named else ""
                listed = reader.child(c, d, list_tag)
                values = reader.i32s(*listed).tolist() if listed else []
                mesh[key].append((name, values))
        elif sid == PARTS_SECTION:
            for pid, c, d in reader.chunks(a, b):
                if pid != PART:
                    continue
                ident = reader.int_of(c, d, PART_ID, -1)
                named = reader.child(c, d, PART_NAME)
                mesh["parts"][ident] = reader.fixed(*named) if named else ""
    if mesh["nodes"] is None:
        _fail("the file has no mesh")
    return mesh


def _parse(filename):
    with open(filename, "rb") as fh:
        data = memoryview(fh.read())
    top = _top_level(data)
    if not top or top[0][1] != ROOT:
        _fail("the file does not start with its root section")
    reader, _, a, b = top[0]
    header = reader.child(a, b, HEADER)
    version = reader.int_of(*header, HDR_VERSION, 0) if header else 0
    if not 0x0030 <= version <= 0x00FF:
        _fail(
            f"plot file version 0x{version:04x} is not supported "
            "(FEBio 3 and 4 write 0x0030 and later)"
        )
    found = reader.child(a, b, DICTIONARY)
    dictionary = (
        _dictionary(reader, *found) if found else {g: [] for g in DIC_GROUPS.values()}
    )
    # A remeshed run writes a mesh section before the states that use it.
    meshes, states, state_mesh = [], [], []
    for t in top[1:]:
        if t[1] == MESH:
            meshes.append(t)
        elif t[1] == STATE:
            states.append(t)
            state_mesh.append(max(len(meshes) - 1, 0))
    if not meshes:
        _fail("the file has no mesh")
    return version, dictionary, (meshes, state_mesh), states


def _mesh_of(parsed, index):
    version, _, (meshes, state_mesh), states = parsed
    mreader, _, c, d = meshes[state_mesh[index] if states else 0]
    return _mesh_section(mreader, c, d, version)


def _state_time(state):
    reader, _, a, b = state
    hdr = reader.child(a, b, STATE_HEADER)
    found = reader.child(*hdr, STATE_TIME) if hdr else None
    time = float(reader.f32s(found[0], found[0] + 4)[0]) if found else 0.0
    status = reader.int_of(*hdr, STATE_STATUS, None) if hdr else None
    return time, status


def time_values(filename):
    """The times of the plot file's states."""
    return [_state_time(s)[0] for s in _parse(filename)[3]]


def read(filename, points_only=False, arrays=None, time_step=0, lenient=False):
    parsed = _parse(str(filename))
    _, dictionary, _, states = parsed
    # the step first: a remeshed run reads the mesh that step uses
    index = 0
    if not states:
        if time_step not in (0, -1):
            _fail(f"time step {time_step} is out of range: the file has no states")
    else:
        n = len(states)
        index = time_step + n if time_step < 0 else time_step
        if not 0 <= index < n:
            raise ReadError(
                f"time step {time_step} is out of range: the file has {n} step(s)"
            )
    raw = _mesh_of(parsed, index)

    # -- mesh ---------------------------------------------------------------------
    cells = []
    file_conn = []  # connectivity in the file's (FEBio's) node order, per domain
    elem_index = {}
    base = 0
    regions = {}
    lossy = 0
    for k, (cell_type, part, name, ids, conn) in enumerate(raw["domains"]):
        if cell_type in _LOSSY:
            keep = _LOSSY[cell_type]
            if not lenient or keep is None:
                hint = " (read with lenient to downgrade it)" if keep else ""
                _fail(f"element type {cell_type} has no meshio++ cell type{hint}")
            lossy += 1
            cell_type = keep
        file_conn.append(conn)
        kept = conn[:, : _NODES[cell_type]]
        order = node_order("febio", cell_type)
        if order is not None:
            kept = kept[:, list(order.to_meshio)]
        cells.append(CellBlock(cell_type, np.ascontiguousarray(kept)))
        for r, ident in enumerate(ids):
            elem_index[ident] = base + r
        # FEBio 4 names only solid domains; a nameless one takes the name of
        # the element set holding exactly its elements (FEBio writes one per
        # <Elements> block), else its part's.
        if not name:
            mine = set(ids)
            name = next((n for n, m in raw["elemsets"] if set(m) == mine), "")
        region = name or raw["parts"].get(part, "") or f"Part{k + 1}"
        entry = regions.setdefault(("cell", region), [part, _dim(cell_type), []])
        entry[1] = max(entry[1], _dim(cell_type))
        entry[2].extend(range(base, base + len(ids)))
        base += len(ids)
    if lossy:
        warn(f"FEBio .xplt: {lossy} domain(s) downgraded to a meshio++ cell type")
    mesh = Mesh(raw["nodes"], cells)

    for name, members in raw["nodesets"]:
        regions.setdefault(("point", name), [-1, -1, []])[2].extend(members)
    for name, members in raw["elemsets"]:
        entry = regions.setdefault(("cell", name), [-1, -1, []])
        entry[2].extend(elem_index[i] for i in members if i in elem_index)

    # Surfaces: a data surface is a block of facet cells per facet type (its
    # variables live there) plus a side region when its facets lie on the cells;
    # a facet set is only the side region, or facet blocks when some facet lies
    # on no cell.
    faces = None
    surface_cells = []  # per data surface, the global cell of each facet
    extra = []  # (type, rows, data surface id or 0)
    seen_sides = set()
    surface_id = 0
    for name, facets, facet_set in raw["surfaces"]:
        if not facet_set:
            surface_id += 1
        if not facets:
            if not facet_set:
                surface_cells.append([])
            continue
        if faces is None:
            faces = FacetIndex(mesh, surface_edges=False)
        sides = []
        for nn, nodes in facets:
            hit = faces.find(nodes[: _CORNERS.get(nn, nn)])
            if hit is None:
                sides = None
                break
            sides.append(hit.first)
        if sides is not None and name not in seen_sides:
            seen_sides.add(name)
            regions.setdefault(("side", name), [-1, 2, []])[2].extend(sides)
        if sides is not None and facet_set:
            continue
        by_type = {}
        for f, (nn, nodes) in enumerate(facets):
            cell_type = _FACET_TYPES.get(nn)
            if cell_type is None or len(nodes) != nn:
                continue
            by_type.setdefault(cell_type, []).append(f)
        cells_of = [-1] * len(facets)
        entry = regions.setdefault(("cell", name), [-1, 2, []])
        for cell_type, members in by_type.items():
            rows = []
            for f in members:
                rows.append(facets[f][1])
                cells_of[f] = base
                entry[2].append(base)
                base += 1
            extra.append((cell_type, rows, 0 if facet_set else surface_id, 0))
        if not facet_set:
            surface_cells.append(cells_of)
    # Edges: a block of line cells per segment type, a cell region each, and
    # each segment's global cell for the edge's data.
    edge_cells = []
    for e, (name, segments) in enumerate(raw["edges"], start=1):
        by_type = {}
        for f, (nn, nodes) in enumerate(segments):
            cell_type = {2: "line", 3: "line3"}.get(nn)
            if cell_type is None or len(nodes) != nn:
                continue
            by_type.setdefault(cell_type, []).append(f)
        cells_of = [-1] * len(segments)
        entry = regions.setdefault(("cell", name), [-1, 1, []])
        entry[1] = max(entry[1], 1)
        for cell_type, members in by_type.items():
            rows = []
            for f in members:
                rows.append(segments[f][1])
                cells_of[f] = base
                entry[2].append(base)
                base += 1
            extra.append((cell_type, rows, 0, e))
        edge_cells.append(cells_of)
    n_domain_blocks = len(mesh.cells)
    for cell_type, rows, _, _ in extra:
        mesh.cells.append(CellBlock(cell_type, np.array(rows, dtype=np.int64)))
    for key, pos in (("xplt:surface", 2), ("xplt:edge", 3)):
        wants_ids = not points_only and (arrays is None or key in arrays)
        if wants_ids and any(e[pos] > 0 for e in extra):
            ids = [
                np.zeros(len(c.data), dtype=np.int64)
                for c in mesh.cells[:n_domain_blocks]
            ]
            ids += [np.full(len(e[1]), e[pos], dtype=np.int64) for e in extra]
            mesh.cell_data[key] = ids

    kind_order = {"point": 0, "cell": 1, "side": 2}
    out = []
    for (kind, name), (tag, dim, entries) in regions.items():
        arr = np.array(entries, dtype=np.int64)
        if kind == "side":
            arr = arr.reshape(-1, 2)
        out.append(Region(name, kind, arr, dim=dim, tag=tag))
    out.sort(key=lambda r: (kind_order[r.kind], r.name, r.dim, r.tag))
    mesh.regions = out

    # -- the chosen state -------------------------------------------------------------
    if not states:
        return mesh
    mesh.time_values = [_state_time(s)[0] for s in states]
    time, status = _state_time(states[index])
    mesh.field_data["meshio:time"] = np.array([time], dtype=np.float64)
    mesh.field_data["xplt:step"] = np.array([index], dtype=np.int64)
    if status is not None:
        mesh.field_data["xplt:status"] = np.array([status], dtype=np.int64)
    if points_only:
        return mesh
    wanted = None if arrays is None else set(arrays)
    data_surfaces = [
        facets for _, facets, facet_set in raw["surfaces"] if not facet_set
    ]
    surfaces = list(zip(data_surfaces, surface_cells))
    edges = list(zip([segments for _, segments in raw["edges"]], edge_cells))
    _read_state(states[index], dictionary, mesh, file_conn, wanted, surfaces, edges)
    return mesh


def _dim(cell_type):
    family = cell_type.rstrip("0123456789")
    return {"line": 1, "triangle": 2, "quad": 2}.get(family, 3)


def _read_state(state, dictionary, mesh, file_conn, wanted, surfaces=(), edges=()):
    reader, _, a, b = state
    found = reader.child(a, b, STATE_DATA)
    if found is None:
        return
    n_points = len(mesh.points)
    sizes = [len(c.data) for c in mesh.cells]
    skipped = []
    point_sums = {}
    for gid, c, d in reader.chunks(*found):
        group = DATA_GROUPS.get(gid)
        if group is None:
            continue
        items = dictionary[group]
        for vid, e, f in reader.chunks(c, d):
            if vid != STATE_VARIABLE:
                continue
            var = reader.int_of(e, f, STATE_VAR_ID, 0)
            if not 1 <= var <= len(items):
                continue
            name, fmt, width = items[var - 1]
            if wanted is not None and name not in wanted:
                continue
            data = reader.child(e, f, STATE_VAR_DATA)
            if data is None:
                continue
            regions = [
                (rid, reader.f32s(g, h).astype(np.float64))
                for rid, g, h in reader.chunks(*data)
            ]
            if group in ("surface", "edge"):
                # region k is data surface (or edge) k: per facet or segment
                # (item) or one value (region) -> its cells; per node (node,
                # first seen over its facets or segments) or per facet node
                # (mult) -> averaged at the points
                carriers = edges if group == "edge" else surfaces
                if fmt in ("item", "region"):
                    bases = np.concatenate([[0], np.cumsum(sizes)])
                    blocks = [np.full((n, width), np.nan) for n in sizes]
                    landed = False
                    for rid, values in regions:
                        k = rid - 1
                        if not 0 <= k < len(carriers):
                            continue
                        cells = carriers[k][1]
                        need = width if fmt == "region" else len(cells) * width
                        if len(values) < need:
                            continue
                        for f, cell in enumerate(cells):
                            if cell < 0:
                                continue
                            landed = True
                            bk = int(np.searchsorted(bases, cell, side="right")) - 1
                            row = (
                                values[:width]
                                if fmt == "region"
                                else values[f * width : (f + 1) * width]
                            )
                            blocks[bk][cell - bases[bk]] = row
                    # a variable on no surface (edge) of this mesh gives no array
                    if landed:
                        mesh.cell_data[name] = [
                            blk[:, 0].copy() if width == 1 else blk for blk in blocks
                        ]
                elif fmt in ("node", "mult"):
                    for k, (facets, _) in enumerate(carriers):
                        for rid, values in regions:
                            if rid != k + 1:
                                continue
                            flat = [p for _, nodes in facets for p in nodes]
                            nodes = flat if fmt == "mult" else list(dict.fromkeys(flat))
                            if not nodes or len(values) < len(nodes) * width:
                                continue
                            total, count = point_sums.setdefault(
                                name, (np.zeros((n_points, width)), np.zeros(n_points))
                            )
                            rows = values[: len(nodes) * width].reshape(
                                len(nodes), width
                            )
                            keep = [i for i, p in enumerate(nodes) if p < n_points]
                            np.add.at(total, [nodes[i] for i in keep], rows[keep])
                            np.add.at(count, [nodes[i] for i in keep], 1)
                elif name not in skipped:
                    skipped.append(name)
            elif group == "global":
                for _, values in regions[:1]:
                    mesh.field_data[name] = values
            elif group == "node":
                for _, values in regions[:1]:
                    if len(values) < n_points * width:
                        continue
                    values = values[: n_points * width].reshape(n_points, width)
                    mesh.point_data[name] = values[:, 0] if width == 1 else values
            elif fmt in ("item", "region"):
                blocks = [np.full((n, width), np.nan) for n in sizes]
                for rid, values in regions:
                    k = rid - 1
                    if not 0 <= k < len(file_conn):
                        continue
                    ne = len(file_conn[k])
                    if len(values) < (width if fmt == "region" else ne * width):
                        continue
                    if fmt == "region":
                        blocks[k][:] = values[:width]
                    else:
                        blocks[k][:] = values[: ne * width].reshape(ne, width)
                mesh.cell_data[name] = [
                    blk[:, 0].copy() if width == 1 else blk for blk in blocks
                ]
            elif fmt in ("node", "mult"):
                total, count = point_sums.setdefault(
                    name, (np.zeros((n_points, width)), np.zeros(n_points))
                )
                for rid, values in regions:
                    k = rid - 1
                    if not 0 <= k < len(file_conn):
                        continue
                    conn = file_conn[k]
                    if fmt == "node":
                        local = list(dict.fromkeys(conn.ravel().tolist()))
                        if len(values) < len(local) * width:
                            continue
                        rows = values[: len(local) * width].reshape(len(local), width)
                        np.add.at(total, local, rows)
                        np.add.at(count, local, 1)
                    else:
                        flat = conn.ravel()
                        if len(values) < len(flat) * width:
                            continue
                        rows = values[: len(flat) * width].reshape(len(flat), width)
                        np.add.at(total, flat, rows)
                        np.add.at(count, flat, 1)
            else:
                if name not in skipped:
                    skipped.append(name)
    for name, (total, count) in point_sums.items():
        if name in mesh.point_data:
            warn(
                f"FEBio .xplt: nodal '{name}' kept; the element-node '{name}' was dropped"
            )
            continue
        with np.errstate(invalid="ignore", divide="ignore"):
            values = total / count[:, None]
        mesh.point_data[name] = values[:, 0] if values.shape[1] == 1 else values
    if skipped:
        warn(
            "FEBio .xplt: material-point variables are not read: " + ", ".join(skipped)
        )

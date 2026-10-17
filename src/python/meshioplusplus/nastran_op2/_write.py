"""Nastran OP2 result files, written.

The twin of the writing half of ``src/cpp/src/formats/nastran_op2.cpp``: both
engines write the same bytes. The file is the reader's inverse for the core of
the format: 32-bit little-endian words (4-byte record markers), the MSC
dialect, SORT1 real tables. It holds the model (``GEOM1`` GRIDs, ``GEOM2``
element cards, ``EPT`` property cards), the nodal tables (``OUGV1``,
``OQG1``, ``OQMG1``, ``OPG1``) and the centre stress and strain tables
(``OES1X1``, ``OSTR1X``) of shells, solids and rods, one block per step.
"""

import re
import struct

import numpy as np

from .. import _provenance
from .._common import warn
from .._exceptions import WriteError
from .._files import is_buffer
from .._fortran_records import fortran_record
from ..nastran._model import CARD_TYPE

TIME_KEY = "meshio:time"

# Fortran records with 4-byte little-endian markers.
_LAYOUT = (4, "<")

# The device code of every approach and row id (print).
_DEVICE = 1

_VERSION = "XXXXXXXX"
_TAPE_ID = "NASTRAN FORT TAPE ID CODE - "

# Nastran numbers the quadratic hexahedron's and wedge's mid-side nodes bottom,
# vertical, top; meshio++ (VTK) bottom, top, vertical: G[perm[k]] = conn[k].
_HEXA20 = list(range(12)) + [16, 17, 18, 19, 12, 13, 14, 15]
_PENTA15 = list(range(9)) + [12, 13, 14, 9, 10, 11]

# card -> ((k1, k2, k3) key, record words, node count of the linear card,
# first node word, node permutation to Nastran's order or None)
_CARDS = {
    "CDAMP1": ((201, 2, 69), 6, 2, 2, None),
    "CELAS1": ((601, 6, 73), 6, 2, 2, None),
    "CONM2": ((1501, 15, 64), 13, 1, 1, None),
    "CONROD": ((1601, 16, 47), 8, 2, 1, None),
    "CHEXA": ((7308, 73, 253), 22, 8, 2, _HEXA20),
    "CPENTA": ((4108, 41, 280), 17, 6, 2, _PENTA15),
    "CPYRAM": ((17200, 172, 1000), 16, 5, 2, None),
    "CQUAD": ((9108, 91, 507), 11, 4, 2, None),
    "CQUAD4": ((2958, 51, 177), 14, 4, 2, None),
    "CQUAD8": ((4701, 47, 326), 17, 4, 2, None),
    "CROD": ((3001, 30, 48), 4, 2, 2, None),
    "CTETRA": ((5508, 55, 217), 12, 4, 2, None),
    "CTRIA3": ((5959, 59, 282), 13, 3, 2, None),
    "CTRIA6": ((4801, 48, 327), 14, 3, 2, None),
    "CTUBE": ((3701, 37, 49), 4, 2, 2, None),
}

# meshio++ cell type -> (card, nodes written); a line takes CROD, CTUBE or
# CONROD by its property.
_TYPES = {
    "vertex": ("CONM2", 1),
    "line": (None, 2),
    "triangle": ("CTRIA3", 3),
    "triangle6": ("CTRIA6", 6),
    "quad": ("CQUAD4", 4),
    "quad8": ("CQUAD8", 8),
    "quad9": ("CQUAD", 9),
    "tetra": ("CTETRA", 4),
    "tetra10": ("CTETRA", 10),
    "pyramid": ("CPYRAM", 5),
    "pyramid13": ("CPYRAM", 13),
    "wedge": ("CPENTA", 6),
    "wedge15": ("CPENTA", 15),
    "hexahedron": ("CHEXA", 8),
    "hexahedron20": ("CHEXA", 20),
}

# The property id of a cell block with no 'nastran:pid': one per family (shells,
# solids, rods), as ids are shared by every property card; a vertex has none.
_DEFAULT_PID = {
    "vertex": -1,
    "line": 3,
    "triangle": 1,
    "triangle6": 1,
    "quad": 1,
    "quad8": 1,
    "quad9": 1,
    "tetra": 2,
    "tetra10": 2,
    "pyramid": 2,
    "pyramid13": 2,
    "wedge": 2,
    "wedge15": 2,
    "hexahedron": 2,
    "hexahedron20": 2,
}

# Property cards: ((k1, k2, k3) key, record words)
_PSHELL = ((2302, 23, 283), 11)
_PSOLID = ((2402, 24, 281), 7)
_PROD = ((902, 9, 29), 6)
_PTUBE = ((1602, 16, 30), 5)
_PELAS = ((302, 3, 46), 4)
_PDAMP = ((202, 2, 45), 2)

# The writer's card for each element type number of ``nastran:card`` on a line or a
# vertex, the cards each of those cell types may hold, and the cards it does not
# write yet (CBAR, CBEAM, CBUSH, CELAS2, CDAMP2, CVISC, PLOTEL, CSHEAR): their rows
# are left out rather than written as another card.
_LINE_CARDS = {
    30: "CROD",
    16: "CONROD",
    37: "CTUBE",
    6: "CELAS1",
    2: "CDAMP1",
    15: "CONM2",
}
_BLOCK_CARDS = {
    "vertex": ("CONM2", "CELAS1", "CDAMP1"),
    "line": ("CROD", "CONROD", "CTUBE", "CELAS1", "CDAMP1"),
}
_CARD_OF_TYPE = {number: name for name, number in CARD_TYPE.items()}
# A cell of another type whose card this writer leaves out: CSHEAR is a quad.
_UNWRITTEN_CELLS = {31: "CSHEAR"}

# The property and material values the field data can set, as
# ``field_data["nastran:<CARD>_<id>:<MEMBER>"]``: card -> {member: (word, default)}.
# ``MID*`` members are integers; the others are floats. Words are those of the
# ``EPT`` and ``MPT`` records; the member names are the ``nastran_h5`` columns.
PROPERTY_VALUES = {
    "PSHELL": {"MID1": (1, 1), "T": (2, 1.0), "NSM": (7, 0.0)},
    "PROD": {
        "MID": (1, 1),
        "A": (2, 1.0),
        "J": (3, 0.0),
        "C": (4, 0.0),
        "NSM": (5, 0.0),
    },
    "PTUBE": {"MID": (1, 1), "OD": (2, 1.0), "T": (3, 0.0), "NSM": (4, 0.0)},
    "PSOLID": {"MID": (1, 1)},
    "PELAS": {"K": (1, 1.0), "GE": (2, 0.0), "S": (3, 0.0)},
    "PDAMP": {"B": (1, 1.0)},
    "MAT1": {
        "E": (1, 1.0),
        "G": (2, 0.0),
        "NU": (3, 0.3),
        "RHO": (4, 0.0),
        "A": (5, 0.0),
        "TREF": (6, 0.0),
        "GE": (7, 0.0),
    },
}
_VALUE_KEY = re.compile(r"^nastran:([A-Z0-9]+)_([0-9]+):([A-Z0-9]+)$")


def property_values(mesh):
    """``{(card, id): {member: value}}`` from the ``nastran:<CARD>_<id>:<MEMBER>``
    field data; keys that name no card or member here are left alone."""
    out = {}
    for name, array in mesh.field_data.items():
        m = _VALUE_KEY.match(name)
        if m is None or m.group(1) not in PROPERTY_VALUES:
            continue
        card, ident, member = m.group(1), int(m.group(2)), m.group(3)
        if member not in PROPERTY_VALUES[card]:
            continue
        values = np.asarray(array).ravel()
        if values.size != 1:
            _wfail(f"'{name}' must hold one value")
        out.setdefault((card, ident), {})[member] = (
            int(values[0]) if member.startswith("MID") else float(values[0])
        )
    return out


def material_ids(properties, values):
    """The MAT1 ids: those the properties name (MID default 1) and those given."""
    ids = set()
    for ptype in ("PSHELL", "PSOLID", "PROD", "PTUBE"):  # the families with a material
        for pid in properties.get(ptype, ()):
            given = values.get((ptype, pid), {})
            ids.add(given.get("MID", given.get("MID1", 1)))
    ids.update(i for (card, i) in values if card == "MAT1")
    return sorted(ids) or []


# Nodal tables in file order: (array, table, table code, thermal flag)
_NODAL = (
    ("DISPLACEMENT", "OUGV1", 1, 0),
    ("EIGENVECTOR", "OUGV1", 7, 0),
    ("VELOCITY", "OUGV1", 10, 0),
    ("ACCELERATION", "OUGV1", 11, 0),
    ("TEMPERATURE", "OUGV1", 1, 1),
    ("SPC_FORCE", "OQG1", 3, 0),
    ("MPC_FORCE", "OQMG1", 3, 0),
    ("APPLIED_LOAD", "OPG1", 2, 0),
)
_NODAL_TABLES = ("OUGV1", "OQG1", "OQMG1", "OPG1")
_RESULT_TABLES = _NODAL_TABLES + ("OEF1X", "OES1X1", "OSTR1X")

# Element result families in file order: (element type, card, members). The
# members are the words after the element id of a row (rods, shells), or the
# values of the centre block (solids).
_ROD = ["A", "MSA", "T", "MST"]
_TUBE = ["AS", "MSA", "TS", "MST"]
_SOLID = [
    "X",
    "Y",
    "Z",
    "TXY",
    "TYZ",
    "TZX",
    "PRINCIPAL_A",
    "PRINCIPAL_B",
    "PRINCIPAL_C",
]
_SOLID_WORD = {
    # member -> word of the 21-word node block (word 0 is the GRID)
    "X": 1,
    "TXY": 2,
    "PRINCIPAL_A": 3,
    "PRESSURE": 7,
    "VON_MISES": 8,
    "OCT_SHEAR": 8,
    "Y": 9,
    "TYZ": 10,
    "PRINCIPAL_B": 11,
    "Z": 15,
    "TZX": 16,
    "PRINCIPAL_C": 17,
}
_PLATE = [
    "FD1", "X1", "Y1", "TXY1", "ANGLE1", "MAJOR1", "MINOR1", "VM1",
    "FD2", "X2", "Y2", "TXY2", "ANGLE2", "MAJOR2", "MINOR2", "VM2",
]  # fmt: skip
_FAMILIES = (
    (1, "CROD", "rod"),
    (3, "CTUBE", "rod"),
    (10, "CONROD", "rod"),
    (33, "CQUAD4", "plate"),
    (39, "CTETRA", "solid"),
    (67, "CHEXA", "solid"),
    (68, "CPENTA", "solid"),
    (74, "CTRIA3", "plate"),
    (255, "CPYRAM", "solid"),
)
_CORNERS = {"CTETRA": 4, "CHEXA": 8, "CPENTA": 6, "CPYRAM": 5}

# Element force families (``OEF1X``, table code 4): (element type, card, members) of
# the rods (axial force and torque) and the plates (the eight membrane and bending
# forces), as ``ELEMENT_FORCE:<member>`` cell data.
_SHELL_FORCE = ["MX", "MY", "MXY", "BMX", "BMY", "BMXY", "TX", "TY"]
_FORCE_FAMILIES = (
    (1, "CROD", ["AF", "TRQ"]),
    (3, "CTUBE", ["AF", "TRQ"]),
    (10, "CONROD", ["AF", "TRQ"]),
    (33, "CQUAD4", _SHELL_FORCE),
    (74, "CTRIA3", _SHELL_FORCE),
)

# Arrays and field data that mean something without being results.
_QUIET_POINT = ("nastran:cp", "nastran:cd", "nastran:nid")
_QUIET_CELL = ("nastran:eid", "nastran:pid", "nastran:card", "nastran:c1", "nastran:c2")

_PTYPE = re.compile(r"^([A-Z][A-Z0-9]*)_([0-9]+)$")


def _wfail(message):
    raise WriteError(f"Nastran OP2: {message}")


# -- framing ---------------------------------------------------------------------------


def _marker(value):
    return fortran_record(struct.pack("<i", value), _LAYOUT)


def _record(payload):
    return _marker(len(payload) // 4) + fortran_record(payload, _LAYOUT)


def _text(text):
    """A name of up to 8 characters, in two words."""
    return text[:8].ljust(8).encode("ascii", "replace")


def _table(name, records):
    """A named table: its name, the 7-word header, then every record behind a
    marker triple, and the closing triple and 0."""
    out = [
        _record(_text(name)),
        _marker(-1),
        _record(struct.pack("<7i", 101, 0, 0, 0, 0, 0, 0)),
    ]
    k = -2
    out += [_marker(k), _marker(1), _marker(0), _record(_text(name))]
    for payload in records:
        k -= 1
        out += [_marker(k), _marker(1), _marker(0), _record(payload)]
    k -= 1
    out += [_marker(k), _marker(1), _marker(0), _marker(0)]
    return b"".join(out)


def _file_header():
    return b"".join(
        [
            _marker(3),
            fortran_record(struct.pack("<3i", 1, 1, 0), _LAYOUT),
            _marker(7),
            fortran_record(_TAPE_ID.encode("ascii"), _LAYOUT),
            _record(_text(_VERSION)),
            _marker(-1),
            _marker(0),
        ]
    )


def _key(k1, k2, k3):
    return struct.pack("<3i", k1, k2, k3)


# -- words -----------------------------------------------------------------------------


def _f32(values):
    with np.errstate(over="ignore"):  # an overflow is an infinity, as in C++
        return np.ascontiguousarray(values, dtype="<f4").view("<u4")


def _i32(values):
    return np.ascontiguousarray(np.asarray(values).astype("<i4")).view("<u4")


def _bytes(words):
    return np.ascontiguousarray(words, dtype="<u4").tobytes()


# -- the model -------------------------------------------------------------------------


def _signature(mesh):
    """What every step of a series must share: the points' count and the
    cell blocks' types and shapes."""
    return (len(mesh.points), [(c.type, np.asarray(c.data).shape) for c in mesh.cells])


def _node_ids(mesh, n):
    """The GRID id of each point: ``point_data["nastran:nid"]``, else 1 to n."""
    given = mesh.point_data.get("nastran:nid")
    if given is None:
        return np.arange(1, n + 1, dtype=np.int64)
    given = np.asarray(given)
    if given.size != n or given.dtype.kind not in "iu":
        _wfail("'nastran:nid' must hold one integer per point")
    ids = given.astype(np.int64).reshape(n)
    if n and ids.min() < 1:
        _wfail("GRID ids must be positive")
    if len(np.unique(ids)) != n:
        _wfail("GRID ids are not unique")
    return ids


def _components(mesh, b, n, ctype):
    """The ``(n, 2)`` component numbers of block ``b``: ``nastran:c1`` and
    ``nastran:c2``, else 1 at each end. A vertex's second end is ground: 0."""
    out = np.ones((n, 2), dtype=np.int64)
    for k, name in enumerate(("nastran:c1", "nastran:c2")):
        if name in mesh.cell_data:
            values = np.asarray(mesh.cell_data[name][b], dtype=np.int64)
            if values.size != n:
                _wfail(f"'{name}' must have one value per cell")
            out[:, k] = values
    if ctype == "vertex":
        out[:, 1] = 0
    return out


def _written(codes, ctype, card):
    """The card of each row of a block that carries ``nastran:card``, the rows to
    keep, and the names of the cards left out. A line or vertex takes the card its
    number names, and a mismatch with its cell type is an error; any other cell
    keeps the card of its type, unless its number is one this writer leaves out."""
    cards = np.empty(len(codes), dtype=object)
    keep = np.ones(len(codes), dtype=bool)
    left = set()
    for i, code in enumerate(np.asarray(codes, dtype=np.int64).tolist()):
        if ctype in _BLOCK_CARDS:
            name = _LINE_CARDS.get(code)
            if name is None:
                keep[i] = False
                left.add(_CARD_OF_TYPE.get(code, f"card {code}"))
            elif name not in _BLOCK_CARDS[ctype]:
                _wfail(f"a '{ctype}' block cannot hold {name} cards")
            else:
                cards[i] = name
        elif code in _UNWRITTEN_CELLS:
            keep[i] = False
            left.add(_UNWRITTEN_CELLS[code])
        else:
            cards[i] = card
    return cards, keep, left


class _Model:
    """The mesh as the file holds it: GRIDs, element cards and properties."""

    def __init__(self, mesh, per_card_ids=False):
        self.npoints = len(mesh.points)
        points = np.asarray(mesh.points, dtype=np.float64).reshape(self.npoints, -1)
        self.points = np.zeros((self.npoints, 3))
        self.points[:, : min(3, points.shape[1])] = points[:, :3]
        self.node_ids = _node_ids(mesh, self.npoints)
        self.signature = _signature(mesh)
        self.values = property_values(mesh)
        ptube = set()
        for region in getattr(mesh, "regions", None) or []:
            m = _PTYPE.match(region.name)
            if m and m.group(1) == "PTUBE":
                ptube.add(int(m.group(2)))
        sequential = 1
        have_eid = "nastran:eid" in mesh.cell_data
        have_pid = "nastran:pid" in mesh.cell_data
        self.blocks = []  # per written block: dict
        self.skipped = []
        self.unwritten = set()  # element cards the file does not hold yet
        have_card = "nastran:card" in mesh.cell_data
        for b, block in enumerate(mesh.cells):
            conn = np.asarray(block.data)
            spec = _TYPES.get(block.type)
            if spec is None:
                if len(conn):
                    self.skipped.append(block.type)
                continue
            n = len(conn)
            card, nodes = spec
            if n and (conn.ndim != 2 or conn.shape[1] != nodes):
                _wfail(f"a '{block.type}' block needs {nodes} nodes per cell")
            if n and (conn.min() < 0 or conn.max() >= self.npoints):
                _wfail(
                    f"a '{block.type}' cell refers to a point the mesh does not have"
                )
            if have_eid:
                eids = np.asarray(mesh.cell_data["nastran:eid"][b], dtype=np.int64)
                if eids.size != n:
                    _wfail("'nastran:eid' must have one value per cell")
            else:
                eids = np.arange(sequential, sequential + n, dtype=np.int64)
            sequential += n
            if have_pid:
                pids = np.asarray(mesh.cell_data["nastran:pid"][b], dtype=np.int64)
                if pids.size != n:
                    _wfail("'nastran:pid' must have one value per cell")
            else:
                pids = np.full(n, _DEFAULT_PID[block.type], dtype=np.int64)
            comps = _components(mesh, b, n, block.type)
            if have_card:
                cards, keep, left = _written(
                    mesh.cell_data["nastran:card"][b], block.type, card
                )
                self.unwritten.update(left)
            else:
                keep = np.ones(n, dtype=bool)
                if card is None:
                    cards = np.where(
                        pids <= 0,
                        "CONROD",
                        np.where(np.isin(pids, sorted(ptube)), "CTUBE", "CROD"),
                    )
                else:
                    cards = np.full(n, card, dtype=object)
            if not keep.all():
                conn, eids, pids = conn[keep], eids[keep], pids[keep]
                comps, cards = comps[keep], cards[keep]
                n = len(eids)
            if not n:
                continue
            self.blocks.append(
                {
                    "index": b,
                    "type": block.type,
                    "conn": conn.astype(np.int64),
                    "nodes": nodes,
                    "eids": eids,
                    "pids": pids,
                    "cards": cards,
                    "comps": comps,
                }
            )
        for blk in self.blocks:
            if np.any(blk["eids"] <= 0):
                _wfail("element ids must be positive")
        every = (
            np.concatenate([blk["eids"] for blk in self.blocks])
            if self.blocks
            else np.zeros(0, dtype=np.int64)
        )
        if per_card_ids:
            # A card's table is keyed by element id; two cards may share one.
            keyed = [
                (str(c), int(e))
                for blk in self.blocks
                for c, e in zip(blk["cards"], blk["eids"])
            ]
            if len(set(keyed)) != len(keyed):
                _wfail("element ids are not unique")
        elif len(np.unique(every)) != len(every):
            _wfail("element ids are not unique")
        # Cards in name order, each one's rows by element id.
        self.cards = {}  # card -> ([(block, rows)], row order by element id)
        for name in sorted(_CARDS):
            parts = []
            for blk in self.blocks:
                rows = np.nonzero(blk["cards"] == name)[0]
                if len(rows):
                    parts.append((blk, rows))
            if not parts:
                continue
            eid = np.concatenate([blk["eids"][r] for blk, r in parts])
            order = np.argsort(eid, kind="stable")
            self.cards[name] = (parts, order)
        # Property ids by family.
        self.properties = {
            "PSHELL": set(),
            "PSOLID": set(),
            "PROD": set(),
            "PTUBE": set(),
            "PELAS": set(),
            "PDAMP": set(),
        }
        for blk in self.blocks:
            for name, ptype in (
                ("CTRIA3", "PSHELL"),
                ("CTRIA6", "PSHELL"),
                ("CQUAD4", "PSHELL"),
                ("CQUAD8", "PSHELL"),
                ("CQUAD", "PSHELL"),
                ("CTETRA", "PSOLID"),
                ("CPENTA", "PSOLID"),
                ("CPYRAM", "PSOLID"),
                ("CHEXA", "PSOLID"),
                ("CROD", "PROD"),
                ("CTUBE", "PTUBE"),
                ("CELAS1", "PELAS"),
                ("CDAMP1", "PDAMP"),
            ):
                rows = blk["cards"] == name
                if rows.any():
                    self.properties[ptype].update(blk["pids"][rows].tolist())
        used = {}
        for ptype, pids in self.properties.items():
            for p in pids:
                if p in used:
                    _wfail(
                        f"property id {p} is used by {used[p]} and {ptype} elements: "
                        "shells, solids and rods need ids of their own"
                    )
                used[p] = ptype
        for p in used:
            if p <= 0:
                _wfail("property ids must be positive")

    def element_rows(self, name):
        """(eids, pids, [(block, rows)], order) of a card; the first two by id."""
        parts, order = self.cards[name]
        eid = np.concatenate([blk["eids"][r] for blk, r in parts])[order]
        pid = np.concatenate([blk["pids"][r] for blk, r in parts])[order]
        return eid, pid, parts, order

    # -- geometry ---------------------------------------------------------------------

    def geom1(self):
        n = self.npoints
        words = np.zeros((n, 8), dtype="<u4")
        order = np.argsort(self.node_ids, kind="stable")  # GEOM1 holds GRIDs by id
        words[:, 0] = _i32(self.node_ids[order])
        words[:, 2:5] = _f32(self.points[order])
        records = [_key(4501, 45, 1) + _bytes(words), _key(65535, 65535, 65535)]
        return _table("GEOM1", records)

    def geom2(self):
        records = []
        for name in sorted(self.cards):
            key, size, linear, first, perm = _CARDS[name]
            parts, order = self.cards[name]
            eid = np.concatenate([blk["eids"][r] for blk, r in parts])
            pid = np.concatenate([blk["pids"][r] for blk, r in parts])
            comp = np.concatenate([blk["comps"][r] for blk, r in parts])
            width = max(blk["nodes"] for blk, _ in parts)
            nodes = np.zeros((len(eid), max(width, 1)), dtype=np.int64)
            at = 0
            for blk, rows in parts:
                conn = self.node_ids[blk["conn"][rows]]
                if blk["nodes"] > linear and perm is not None:
                    placed = np.zeros_like(conn)
                    for i, p in enumerate(perm):
                        placed[:, p] = conn[:, i]
                    conn = placed
                nodes[at : at + len(rows), : conn.shape[1]] = conn
                at += len(rows)
            eid, pid, nodes, comp = eid[order], pid[order], nodes[order], comp[order]
            words = np.zeros((len(eid), size), dtype="<u4")
            words[:, 0] = _i32(eid)
            if name == "CONM2":
                words[:, 1] = _i32(nodes[:, 0])
            elif name == "CONROD":
                words[:, 1] = _i32(nodes[:, 0])
                words[:, 2] = _i32(nodes[:, 1])
                words[:, 3] = _i32(np.ones(len(eid)))  # MID
                words[:, 4] = _f32(np.ones(len(eid)))  # A
            else:
                words[:, 1] = _i32(pid)
                words[:, 2 : 2 + nodes.shape[1]] = _i32(nodes)
                if name in ("CELAS1", "CDAMP1"):
                    words[:, 4:6] = _i32(comp)
            records.append(_key(*key) + _bytes(words))
        records.append(_key(65535, 65535, 65535))
        return _table("GEOM2", records)

    def ept(self):
        records = []
        for ptype, (key, size) in (
            ("PROD", _PROD),
            ("PSHELL", _PSHELL),
            ("PSOLID", _PSOLID),
            ("PTUBE", _PTUBE),
            ("PELAS", _PELAS),
            ("PDAMP", _PDAMP),
        ):
            pids = sorted(self.properties[ptype])
            if not pids:
                continue
            words = np.zeros((len(pids), size), dtype="<u4")
            words[:, 0] = _i32(pids)
            for k, pid in enumerate(pids):
                given = self.values.get((ptype, pid), {})
                for member, (word, default) in PROPERTY_VALUES[ptype].items():
                    value = given.get(member, default)
                    if member.startswith("MID"):
                        words[k, word] = _i32(np.array([value]))[0]
                    else:
                        words[k, word] = _f32(np.array([value]))[0]
            if ptype == "PSOLID":
                words[:, 6] = 0x43454D53  # FCTN 'SMEC'
            records.append(_key(*key) + _bytes(words))
        records.append(_key(65535, 65535, 65535))
        return _table("EPT", records)

    def mpt(self):
        """A MAT1 per material the properties name: unit modulus unless the field
        data gives its values."""
        ids = material_ids(self.properties, self.values) or [1]
        words = np.zeros((len(ids), 12), dtype="<u4")
        for k, mid in enumerate(ids):
            words[k, 0] = _i32(np.array([mid]))[0]
            given = self.values.get(("MAT1", mid), {})
            for member, (word, default) in PROPERTY_VALUES["MAT1"].items():
                words[k, word] = _f32(np.array([given.get(member, default)]))[0]
        return _table(
            "MPT", [_key(103, 1, 77) + _bytes(words), _key(65535, 65535, 65535)]
        )


# -- steps -----------------------------------------------------------------------------


def _scalar(mesh, name, default):
    array = mesh.field_data.get(name)
    if array is None:
        return default
    values = np.asarray(array).ravel()
    return values[0] if values.size else default


class _Header:
    """The step words every result header carries."""

    def __init__(self, mesh, time):
        self.time = float(np.float32(time))
        self.subcase = int(_scalar(mesh, "nastran:subcase", 1))
        default = 1 if self.time == 0.0 else 6
        self.analysis = int(_scalar(mesh, "nastran:analysis", default))
        self.mode = int(_scalar(mesh, "nastran:mode", 1))
        if self.subcase <= 0 or self.analysis <= 0 or self.analysis > 99:
            _wfail("'nastran:subcase' and 'nastran:analysis' must be positive codes")
        w5, w6 = 1, 0.0
        if self.analysis in (5, 6, 10, 12):
            w5 = struct.unpack("<i", struct.pack("<f", self.time))[0]
            w6 = 0.0
            self.w5f = self.time
        elif self.analysis in (2, 8, 9):
            w5, w6 = self.mode, self.time
            self.w5f = 0.0
        else:
            self.w5f = 0.0
        self.w5 = w5
        self.w6 = w6

    def key(self):
        return (self.subcase, self.analysis, self.w5)

    def words(self, tcode, etype, num_wide, s_code=0, thermal=0, format_code=1):
        w = np.zeros(146, dtype="<u4")
        w[50:] = 0x20202020  # the title, subtitle and label: blanks
        w[0] = self.analysis * 10 + _DEVICE
        w[1] = tcode
        w[2] = etype
        w[3] = self.subcase
        w[4] = np.uint32(self.w5 & 0xFFFFFFFF)
        w[5:6] = _f32(np.array([self.w6]))
        w[8] = format_code
        w[9] = num_wide
        w[10] = s_code
        w[22] = thermal
        return _bytes(w)


def _step_tables(model, mesh, time, taken):
    """The step's header and the records of its result tables: ``{table:
    [payloads]}``, and the ``(table, table code, element type)`` of every block."""
    header = _Header(mesh, time)
    tables = {}
    layout = []
    # Nodal tables.
    for name, table, tcode, thermal in _NODAL:
        base = mesh.point_data.get(name)
        if base is None:
            continue
        taken.add(name)
        base = np.asarray(base, dtype=np.float64)
        ncomp = 1 if thermal else 3
        if base.size != model.npoints * ncomp:
            _wfail(f"'{name}' must have {ncomp} component(s) per point")
        base = base.reshape(model.npoints, ncomp)
        rot = None
        if not thermal:
            rot = mesh.point_data.get(name + "_ROT")
            if rot is not None:
                taken.add(name + "_ROT")
                rot = np.asarray(rot, dtype=np.float64)
                if rot.size != model.npoints * 3:
                    _wfail(f"'{name}_ROT' must have 3 components per point")
                rot = rot.reshape(model.npoints, 3)
        keep = np.nonzero(~np.all(np.isnan(base), axis=1))[0]
        if not len(keep):
            continue
        keep = keep[np.argsort(model.node_ids[keep], kind="stable")]  # rows by id
        words = np.zeros((len(keep), 8), dtype="<u4")
        words[:, 0] = _i32(model.node_ids[keep] * 10 + _DEVICE)
        words[:, 1] = 1
        words[:, 2 : 2 + ncomp] = _f32(base[keep])
        if rot is not None:
            words[:, 5:8] = _f32(rot[keep])
        layout.append((table, tcode, 0))
        tables.setdefault(table, []).extend(
            [header.words(tcode, 0, 8, 0, thermal), _bytes(words)]
        )
    # Complex nodal tables: real and imaginary parts (``<name>_real``/``_imag`` and
    # their ``_ROT_`` rotations), table code 1000 + the real one, format code 2.
    for name, table, tcode, thermal in _NODAL:
        if thermal or name + "_real" not in mesh.point_data:
            continue
        if name + "_imag" not in mesh.point_data:
            continue
        parts = []
        for suffix in ("_real", "_imag", "_ROT_real", "_ROT_imag"):
            array = mesh.point_data.get(name + suffix)
            if array is None:
                parts.append(None)
                continue
            taken.add(name + suffix)
            array = np.asarray(array, dtype=np.float64)
            if array.size != model.npoints * 3:
                _wfail(f"'{name}{suffix}' must have 3 components per point")
            parts.append(array.reshape(model.npoints, 3))
        real, imag, rot_real, rot_imag = parts
        shown = np.isnan(real).all(axis=1) & np.isnan(imag).all(axis=1)
        keep = np.nonzero(~shown)[0]
        if not len(keep):
            continue
        keep = keep[np.argsort(model.node_ids[keep], kind="stable")]  # rows by id
        words = np.zeros((len(keep), 14), dtype="<u4")
        words[:, 0] = _i32(model.node_ids[keep] * 10 + _DEVICE)
        words[:, 1] = 1
        words[:, 2:5] = _f32(real[keep])
        words[:, 8:11] = _f32(imag[keep])
        if rot_real is not None:
            words[:, 5:8] = _f32(rot_real[keep])
        if rot_imag is not None:
            words[:, 11:14] = _f32(rot_imag[keep])
        layout.append((table, 1000 + tcode, 0))
        tables.setdefault(table, []).extend(
            [header.words(1000 + tcode, 0, 14, 0, 0, 2), _bytes(words)]
        )
    # Element tables.
    _force_tables(model, mesh, header, tables, taken, layout)
    for group, table in (("STRESS", "OES1X1"), ("STRAIN", "OSTR1X")):
        _element_tables(model, mesh, header, group, table, tables, taken, layout)
    return header, tables, layout


def _gather(mesh, name, parts, order, width=None):
    """A cell-data array of the cells of a card, in the card's row order:
    ``(n,)``, or ``(n, width)`` for a per-corner array; ``None`` when no block
    of the card has it."""
    arrays = mesh.cell_data.get(name)
    if arrays is None:
        return None
    shape = (sum(len(rows) for _, rows in parts),) + (() if width is None else (width,))
    out = np.full(shape, np.nan)
    at = 0
    for blk, rows in parts:
        a = np.asarray(arrays[blk["index"]], dtype=np.float64)
        n = len(blk["eids"])
        if width is None:
            if n and (a.size < n or a.size % n):
                _wfail(f"'{name}' must have one value per cell")
            out[at : at + len(rows)] = a.reshape(n, -1)[rows, 0]
        elif a.ndim == 2 and a.shape[0] == n:
            k = min(width, a.shape[1])
            out[at : at + len(rows), :k] = a[rows][:, :k]
        at += len(rows)
    return out[order]


def _force_tables(model, mesh, header, tables, taken, layout):
    """``ELEMENT_FORCE:<member>`` of the rods and plates, one ``OEF1X`` block each."""
    for etype, card, members in _FORCE_FAMILIES:
        if card not in model.cards:
            continue
        eid, _pid, parts, order = model.element_rows(card)
        columns = {}
        for m in members:
            column = _gather(mesh, "ELEMENT_FORCE:" + m, parts, order)
            if column is not None:
                taken.add("ELEMENT_FORCE:" + m)
                columns[m] = column
        if not columns:
            continue
        shown = np.zeros(len(eid), dtype=bool)
        for column in columns.values():
            shown |= ~np.isnan(column)
        keep = np.nonzero(shown)[0]
        if not len(keep):
            continue
        wide = 1 + len(members)
        words = np.zeros((len(keep), wide), dtype="<u4")
        words[:, 0] = _i32(eid[keep] * 10 + _DEVICE)
        for j, m in enumerate(members):
            if m in columns:
                words[:, 1 + j] = _f32(np.nan_to_num(columns[m][keep], nan=0.0))
        layout.append(("OEF1X", 4, etype))
        tables.setdefault("OEF1X", []).extend(
            [header.words(4, etype, wide, 0), _bytes(words)]
        )


def _element_tables(model, mesh, header, group, table, tables, taken, layout):
    prefix = group + ":"
    for etype, card, kind in _FAMILIES:
        if card not in model.cards:
            continue
        eid, _pid, parts, order = model.element_rows(card)
        count = len(eid)
        if kind == "rod":
            members = _TUBE if card == "CTUBE" else _ROD
        elif kind == "plate":
            shear = any(
                prefix + f"MAX_SHEAR{k}" in mesh.cell_data
                and prefix + f"VON_MISES{k}" not in mesh.cell_data
                for k in (1, 2)
            )
            vm = "MAX_SHEAR" if shear else "VON_MISES"
            members = [m.replace("VM", vm) for m in _PLATE]
        else:
            octa = (
                prefix + "OCT_SHEAR" in mesh.cell_data
                and prefix + "VON_MISES" not in mesh.cell_data
            )
            members = _SOLID + ["PRESSURE", "OCT_SHEAR" if octa else "VON_MISES"]
        columns = {}  # member -> (count,) values; absent members are not here
        for m in members:
            column = _gather(mesh, prefix + m, parts, order)
            if column is not None:
                taken.add(prefix + m)
                columns[m] = column
        if kind == "solid":
            ncorner = _CORNERS[card]
            corners = {}
            for m in members:
                matrix = _gather(mesh, prefix + m + "@corner", parts, order, ncorner)
                if matrix is not None:
                    taken.add(prefix + m + "@corner")
                    corners[m] = matrix
        shown = np.zeros(count, dtype=bool)
        for column in columns.values():
            shown |= ~np.isnan(column)
        keep = np.nonzero(shown)[0]
        if not len(keep):
            continue

        def value(m, rows=keep):
            """The member's values at the kept rows: 0 where absent or NaN."""
            if m not in columns:
                return np.zeros(len(rows))
            return np.nan_to_num(columns[m][rows], nan=0.0)

        bit0 = 0 if members[-1] in ("OCT_SHEAR", "MAX_SHEAR2") else 1
        s_code = (10 if group == "STRAIN" else 0) | bit0  # strain: bits 8 and 2
        if kind == "solid":
            wide = 4 + 21 * (1 + _CORNERS[card])
        else:
            wide = 1 + len(members)
        words = np.zeros((len(keep), wide), dtype="<u4")
        words[:, 0] = _i32(eid[keep] * 10 + _DEVICE)
        if kind != "solid":
            for j, m in enumerate(members):
                words[:, 1 + j] = _f32(value(m))
        else:
            ncorner = _CORNERS[card]
            ids = np.concatenate(
                [model.node_ids[blk["conn"][r][:, :ncorner]] for blk, r in parts]
            )
            ids = ids[order][keep]
            words[:, 2] = 0x44495247  # 'GRID'
            words[:, 3] = ncorner
            for k in range(1 + ncorner):
                base = 4 + 21 * k
                words[:, base] = 0x2F4E4543 if k == 0 else _i32(ids[:, k - 1])
                for m in members:
                    v = value(m)
                    if k and m in corners:
                        c = corners[m][keep, k - 1]
                        v = np.where(np.isnan(c), v, c)
                    words[:, base + _SOLID_WORD[m]] = _f32(v)
        layout.append((table, 5, etype))
        tables.setdefault(table, []).extend(
            [header.words(5, etype, wide, s_code), _bytes(words)]
        )


def _dropped(model, mesh, taken):
    """Names of arrays that have no place in the file."""
    names = []
    for name in mesh.point_data:
        if name not in taken and name not in _QUIET_POINT:
            names.append(name)
    for name in mesh.cell_data:
        if name not in taken and name not in _QUIET_CELL:
            names.append(name)
    for name in mesh.field_data:
        if name != TIME_KEY and not name.startswith("nastran:"):
            names.append(name)
    return sorted(set(names))


def _warn_dropped(model, mesh, taken):
    names = _dropped(model, mesh, taken)
    if names:
        warn(
            "Nastran OP2: arrays without a place in the file are not written: "
            + ", ".join(names)
        )
    if model.skipped:
        warn(
            "Nastran OP2: cell blocks with no element card are not written: "
            + ", ".join(sorted(set(model.skipped)))
        )
    if model.unwritten:
        warn(
            "Nastran OP2: element cards not written yet, their elements are left out: "
            + ", ".join(sorted(model.unwritten))
        )


def _time_of(mesh):
    array = mesh.field_data.get(TIME_KEY)
    if array is None:
        return 0.0
    values = np.asarray(array, dtype=np.float64).ravel()
    return float(values[0]) if values.size else 0.0


def _assemble(model, steps):
    """The whole file: header, model tables, then every result table with all
    steps' blocks, step by step."""
    out = [_file_header(), model.geom1(), model.geom2(), model.ept()]
    if any(model.properties.values()):
        out.append(model.mpt())
    for table in _RESULT_TABLES:
        records = []
        for tables in steps:
            records.extend(tables.get(table, []))
        if records:
            out.append(_table(table, records))
    out.append(_marker(0))  # the end of the file
    return b"".join(out)


@_provenance.slotless_writer
def write(filename, mesh):
    model = _Model(mesh)
    taken = set()
    _header, tables, _layout = _step_tables(model, mesh, _time_of(mesh), taken)
    _warn_dropped(model, mesh, taken)
    data = _assemble(model, [tables])
    if is_buffer(filename, "w"):
        filename.write(data)
    else:
        with open(filename, "wb") as fh:
            fh.write(data)


class SeriesWriter:
    """Write a time series into one OP2 file: the model once, and every
    table's block per step. The tables hold every step, so the file is
    assembled when the writer closes. A step's time is the argument, not its
    ``meshio:time``; its subcase, analysis and mode come from its
    ``nastran:subcase``, ``nastran:analysis`` and ``nastran:mode``. Every step
    must have the first step's points count, cell blocks and result arrays, and
    a subcase, analysis and mode or time of its own. Use as a context manager,
    or call :meth:`close`."""

    def __init__(self, filename):
        self._filename = filename
        self._model = None
        self._keys = []
        self._count = 0
        self._layout = None
        self._tables = []

    def write(self, time, mesh):
        taken = set()
        index = self._count
        if self._model is None:
            self._model = _Model(mesh)
        elif self._model.signature != _signature(mesh):
            _wfail(
                f"series: step {index}'s points or cells differ from the first "
                "step's; an OP2 file holds one model -- write one file per step "
                "with '{step}'"
            )
        header, tables, layout = _step_tables(self._model, mesh, float(time), taken)
        if self._layout is None:
            self._layout = layout
            _warn_dropped(self._model, mesh, taken)
        elif layout != self._layout:
            _wfail(
                f"series: step {index}'s result arrays differ from the first "
                "step's; a result table holds one block per step"
            )
        # A step with no result table is no step of the file: its key is free.
        if layout and header.key() in self._keys:
            _wfail(
                f"series: step {index} has the subcase, analysis and load set, "
                "mode or time of an earlier step; give it a 'nastran:subcase' "
                "or a different time"
            )
        if layout:
            self._keys.append(header.key())
        self._count += 1
        self._tables.append(tables)

    def close(self):
        if self._filename is None:
            return
        filename, self._filename = self._filename, None
        if not self._count:
            return
        data = _assemble(self._model, self._tables)
        if is_buffer(filename, "w"):
            filename.write(data)
        else:
            with open(filename, "wb") as fh:
                fh.write(data)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        if exc[0] is None:
            self.close()
        else:
            self._filename = None
        return False

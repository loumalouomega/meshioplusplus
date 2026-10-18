"""Nastran OP2 result files, written.

The twin of the writing half of ``src/cpp/src/formats/nastran_op2.cpp``: both
engines write the same bytes. The file is the reader's inverse for the core of
the format: 32-bit little-endian words (4-byte record markers), the MSC
dialect, SORT1 real tables. It holds the model (``GEOM1`` GRIDs, ``GEOM2``
element cards, ``EPT`` property cards), the nodal tables (``OUGV1``,
``OQG1``, ``OQMG1``, ``OPG1``) and the centre stress and strain tables
(``OES1X1``, ``OSTR1X``) of shells, solids and rods, one block per step.
"""

import math
import re
import struct

import numpy as np

from .. import _provenance
from .._common import warn
from .._exceptions import WriteError
from .._files import is_buffer
from .._fortran_records import fortran_record
from ..nastran._model import CARD_TYPE, pbeam_station_slots

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
    "CBAR": ((2408, 24, 180), 16, 2, 2, None),
    "CBEAM": ((5408, 54, 261), 18, 2, 2, None),
    "CBUSH": ((2608, 26, 60), 14, 2, 2, None),
    "CDAMP2": ((301, 3, 70), 6, 2, 2, None),
    "CELAS2": ((701, 7, 74), 8, 2, 2, None),
    "CGAP": ((1908, 19, 104), 9, 2, 2, None),
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
    "CSHEAR": ((3101, 31, 61), 6, 4, 2, None),
    "CTETRA": ((5508, 55, 217), 12, 4, 2, None),
    "CTRIA3": ((5959, 59, 282), 13, 3, 2, None),
    "CTRIA6": ((4801, 48, 327), 14, 3, 2, None),
    "CTUBE": ((3701, 37, 49), 4, 2, 2, None),
    "CVISC": ((3901, 39, 50), 4, 2, 2, None),
    "PLOTEL": ((5201, 52, 11), 3, 2, 1, None),
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
_PBAR = ((52, 20, 181), 19)
# A PBEAM is a header (pid, mid, segments, ccf, a zero) and 12 slots of 16 words, the first 11
# holding the stations (see ``pbeam_station_slots``) and the last one a spare:
# SO, XXB, A, I1, I2, I12, J and NSM, then the outline points C1, C2, D1 and D2 as x, y.
_PBEAM = ((5402, 54, 262), 197)
# A PBUSH is the MSC layout of 23 words: the id, six stiffnesses K, six damping B, six GE
# and the stress and strain coefficients SA, ST, EA and ET.
_PBUSH = ((1402, 14, 37), 23)
_BUSH_ARRAYS = (("K", 1), ("B", 7), ("GE", 13))
_BEAM_STATION_MEMBERS = ("SO", "XXB", "A", "I1", "I2", "I12", "J", "NSM")
_BEAM_POINT_MEMBERS = ("C1", "C2", "D1", "D2")
_PELAS = ((302, 3, 46), 4)
_PDAMP = ((202, 2, 45), 2)
_PVISC = ((1802, 18, 31), 3)
_PGAP = ((2102, 21, 121), 11)
_PSHEAR = ((1002, 10, 42), 6)
_PCOMP = ((2706, 27, 287), 8)

# The writer's card for each element type number of ``nastran:card`` on a line or a
# vertex, the cards each of those cell types may hold, and the cards it does not
# write yet: their rows are left out rather than written as another card.
_LINE_CARDS = {
    30: "CROD",
    16: "CONROD",
    37: "CTUBE",
    6: "CELAS1",
    2: "CDAMP1",
    15: "CONM2",
    24: "CBAR",
    54: "CBEAM",
    26: "CBUSH",
    3: "CDAMP2",
    7: "CELAS2",
    39: "CVISC",
    52: "PLOTEL",
    19: "CGAP",
}
_BLOCK_CARDS = {
    "vertex": ("CONM2", "CELAS1", "CDAMP1", "CELAS2", "CDAMP2"),
    "line": (
        "CROD",
        "CONROD",
        "CTUBE",
        "CELAS1",
        "CDAMP1",
        "CBAR",
        "CBEAM",
        "CBUSH",
        "CELAS2",
        "CDAMP2",
        "CVISC",
        "PLOTEL",
        "CGAP",
    ),
}
_CARD_OF_TYPE = {number: name for name, number in CARD_TYPE.items()}
# The element type number of a CSHEAR, which only a quad holds.
_CSHEAR = 31

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
    "PBAR": {
        "MID": (1, 1),
        "A": (2, 1.0),
        "I1": (3, 0.0),
        "I2": (4, 0.0),
        "J": (5, 0.0),
        "NSM": (6, 0.0),
        "C1": (8, 0.0),
        "C2": (9, 0.0),
        "D1": (10, 0.0),
        "D2": (11, 0.0),
        "E1": (12, 0.0),
        "E2": (13, 0.0),
        "F1": (14, 0.0),
        "F2": (15, 0.0),
        "K1": (16, 0.0),
        "K2": (17, 0.0),
        "I12": (18, 0.0),
    },
    "PELAS": {"K": (1, 1.0), "GE": (2, 0.0), "S": (3, 0.0)},
    "PDAMP": {"B": (1, 1.0)},
    "PVISC": {"CE": (1, 1.0), "CR": (2, 0.0)},
    "PGAP": {
        "U0": (1, 0.0),
        "F0": (2, 0.0),
        "KA": (3, 1.0e8),
        "KB": (4, 1.0e-6),
        "KT": (5, 0.0),
        "MU1": (6, 0.0),
        "MU2": (7, 0.0),
        "TMAX": (8, 0.0),
        "MAR": (9, 100.0),
        "TRMIN": (10, 0.001),
    },
    "PSHEAR": {
        "MID": (1, 1),
        "T": (2, 1.0),
        "NSM": (3, 0.0),
        "F1": (4, 0.0),
        "F2": (5, 0.0),
    },
    # The header words of a PCOMP (its layer count is negative for a symmetric layup,
    # ``SYM``), and ``PLY``: the plies as ``(n, 4)`` rows of MID, T, THETA and SOUT.
    "PCOMP": {
        "Z0": (2, 0.0),
        "NSM": (3, 0.0),
        "SB": (4, 0.0),
        "FT": (5, 0),
        "TREF": (6, 0.0),
        "GE": (7, 0.0),
        "SYM": (None, 0),
        "PLY": (None, None),
    },
    "MAT1": {
        "E": (1, 1.0),
        "G": (2, 0.0),
        "NU": (3, 0.3),
        "RHO": (4, 0.0),
        "A": (5, 0.0),
        "TREF": (6, 0.0),
        "GE": (7, 0.0),
    },
    # The anisotropic shell material (17 words): the stiffness matrix, the thermal
    # expansion and the stress limits; MCSID is an integer.
    "MAT2": {
        "G11": (1, 1.0),
        "G12": (2, 0.0),
        "G13": (3, 0.0),
        "G22": (4, 1.0),
        "G23": (5, 0.0),
        "G33": (6, 1.0),
        "RHO": (7, 0.0),
        "A1": (8, 0.0),
        "A2": (9, 0.0),
        "A3": (10, 0.0),
        "TREF": (11, 0.0),
        "GE": (12, 0.0),
        "ST": (13, 0.0),
        "SC": (14, 0.0),
        "SS": (15, 0.0),
        "MCSID": (16, 0),
    },
    # The orthotropic ply material (19 words).
    "MAT8": {
        "E1": (1, 1.0),
        "E2": (2, 1.0),
        "NU12": (3, 0.3),
        "G12": (4, 0.0),
        "G1Z": (5, 0.0),
        "G2Z": (6, 0.0),
        "RHO": (7, 0.0),
        "A1": (8, 0.0),
        "A2": (9, 0.0),
        "TREF": (10, 0.0),
        "XT": (11, 0.0),
        "XC": (12, 0.0),
        "YT": (13, 0.0),
        "YC": (14, 0.0),
        "S": (15, 0.0),
        "GE": (16, 0.0),
        "F12": (17, 0.0),
        "STRN": (18, 0.0),
    },
}
# The material cards, with their MPT keys and record words.
_MATERIALS = (
    ("MAT1", (103, 1, 77), 12),
    ("MAT2", (203, 2, 78), 17),
    ("MAT8", (2503, 25, 288), 19),
)
_VALUE_KEY = re.compile(r"^nastran:([A-Z0-9]+)_([0-9]+):([A-Z0-9]+)$")
# The members that are integers besides the ``MID`` ones.
_INT_MEMBERS = ("FT", "SYM", "MCSID")


def is_int_member(member):
    """Whether a property or material member holds an integer."""
    return member.startswith("MID") or member in _INT_MEMBERS


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
        if member == "PLY":
            plies = np.asarray(array, dtype=np.float64)
            if plies.ndim != 2 or plies.shape[1] != 4 or not len(plies):
                _wfail(f"'{name}' must hold (plies, 4) rows of MID, T, THETA and SOUT")
            out.setdefault((card, ident), {})[member] = plies
            continue
        values = np.asarray(array).ravel()
        if values.size != 1:
            _wfail(f"'{name}' must hold one value")
        out.setdefault((card, ident), {})[member] = (
            int(values[0]) if is_int_member(member) else float(values[0])
        )
    return out


def material_ids(properties, values, beams=None):
    """The MAT1 ids: those the properties name (MID default 1) and those given. A PBEAM's
    MID is in the beams' field data (``{pid: {member: array}}``), not in the values."""
    ids = set()
    for ptype in (
        "PSHELL",
        "PSOLID",
        "PROD",
        "PTUBE",
        "PBAR",
        "PSHEAR",
    ):  # with a material
        for pid in properties.get(ptype, ()):
            given = values.get((ptype, pid), {})
            ids.add(given.get("MID", given.get("MID1", 1)))
    for pid in properties.get("PBEAM", ()):
        mid = (beams or {}).get(pid, {}).get("MID")
        ids.add(1 if mid is None else int(np.asarray(mid).ravel()[0]))
    for pid in properties.get("PCOMP", ()):
        plies = np.asarray(values[("PCOMP", pid)]["PLY"]).reshape(-1, 4)
        ids.update(int(mid) for mid in plies[:, 0])
    ids.update(material_cards(values))
    return sorted(ids) or []


def material_cards(values):
    """``{id: card}`` of the materials the field data gives (``nastran:MAT1_<id>:*``,
    ``MAT2`` or ``MAT8``); an id given as two cards is an error."""
    cards = {}
    for card, ident in sorted(values):
        if card not in ("MAT1", "MAT2", "MAT8"):
            continue
        if cards.setdefault(ident, card) != card:
            _wfail(
                f"material id {ident} is given as both {cards[ident]} and {card} "
                "field data"
            )
    return cards


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
_QUIET_CELL = (
    "nastran:eid",
    "nastran:pid",
    "nastran:card",
    "nastran:c1",
    "nastran:c2",
    "nastran:x",
    "nastran:g0",
    "nastran:pin",
    "nastran:w",
    "nastran:k",
    "nastran:ge",
    "nastran:s",
    "nastran:b",
)

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


def _raw_bar_fields(mesh, b, n):
    """The per-cell values of block ``b`` that the card fields are made of, as a dict
    of arrays (``None`` where the mesh has none): ``x`` ``(n, 3)``, ``g0``, ``pin``
    ``(n, 2)``, ``w`` ``(n, 6)``, ``k``, ``ge``, ``s`` and ``b``."""
    widths = {"x": 3, "g0": 1, "pin": 2, "w": 6, "k": 1, "ge": 1, "s": 1, "b": 1}
    names = {
        "x": "nastran:x",
        "g0": "nastran:g0",
        "pin": "nastran:pin",
        "w": "nastran:w",
        "k": "nastran:k",
        "ge": "nastran:ge",
        "s": "nastran:s",
        "b": "nastran:b",
    }
    out = {}
    for key, name in names.items():
        if name not in mesh.cell_data:
            out[key] = None
            continue
        values = np.asarray(mesh.cell_data[name][b], dtype=np.float64)
        if values.size != n * widths[key]:
            _wfail(f"'{name}' must have {widths[key]} value(s) per cell")
        out[key] = values.reshape(n, widths[key])
    return out


def _bar_matrix(raw, cards):
    """The ``(n, 12)`` card fields of the rows, by their cards: a CBAR's orientation
    X1..X3, G0, pins PA and PB and offsets W1A..W3A and W1B..W3B, and a CBUSH's X1..X3
    and G0; a CELAS2's K, GE and S in columns 0 to 2; a CDAMP2's B in column 0. Each row
    takes only its own card's."""
    n = len(cards)
    out = np.zeros((n, 12))
    oriented = np.isin(cards, ["CBAR", "CBEAM", "CBUSH", "CGAP"])
    if oriented.any():
        for key, lo, hi in (("x", 0, 3), ("g0", 3, 4)):
            if raw[key] is not None:
                out[oriented, lo:hi] = raw[key][oriented]
    bars = np.isin(cards, ["CBAR", "CBEAM"])  # a CBUSH has no pins or offsets
    if bars.any():
        for key, lo, hi in (("pin", 4, 6), ("w", 6, 12)):
            if raw[key] is not None:
                out[bars, lo:hi] = raw[key][bars]
    springs = np.asarray(cards == "CELAS2")
    if springs.any():
        for key, col in (("k", 0), ("ge", 1), ("s", 2)):
            if raw[key] is not None:
                out[springs, col] = raw[key][springs, 0]
    dampers = np.asarray(cards == "CDAMP2")
    if dampers.any() and raw["b"] is not None:
        out[dampers, 0] = raw["b"][dampers, 0]
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
        elif code == _CSHEAR and ctype == "quad":
            cards[i] = "CSHEAR"
        elif code == _CSHEAR:
            keep[i] = False
            left.add("CSHEAR")
        else:
            cards[i] = card
    return cards, keep, left


# Degrees to radians: the reader's constant, so both engines round the same way.
_DEGREE = 3.14159265358979323846 / 180.0


def _cord_system(mesh, cid):
    """System ``cid``'s type, origin A and axes X, Y, Z (basic), from its
    ``nastran:CORD2_<cid>:*`` field data."""
    prefix = f"nastran:CORD2_{cid}:"
    if any(prefix + m not in mesh.field_data for m in ("TYPE", "A", "X", "Y", "Z")):
        _wfail(
            f"GRIDs are in coordinate system {cid}, which has no '{prefix}*' field data"
        )
    kind = np.asarray(mesh.field_data[prefix + "TYPE"]).ravel()
    if kind.size != 1:
        _wfail(f"'{prefix}TYPE' must have one value")
    if int(kind[0]) not in (1, 2, 3):
        _wfail(f"'{prefix}TYPE' must be 1, 2 or 3")
    vectors = []
    for m in ("A", "X", "Y", "Z"):
        value = np.asarray(mesh.field_data[prefix + m], dtype=np.float64).ravel()
        if value.size != 3:
            _wfail(f"'{prefix}{m}' must have 3 values")
        vectors.append([float(v) for v in value])
    return (int(kind[0]), *vectors)


def _to_local(system, point):
    """The basic ``point`` in ``system``'s local coordinates: the inverse of the reader's
    transform, in the operation order of the C++ core."""
    ctype, o, x, y, z = system
    d0, d1, d2 = point[0] - o[0], point[1] - o[1], point[2] - o[2]
    xl = d0 * x[0] + d1 * x[1] + d2 * x[2]
    yl = d0 * y[0] + d1 * y[1] + d2 * y[2]
    zl = d0 * z[0] + d1 * z[1] + d2 * z[2]
    if ctype == 1:
        return [xl, yl, zl]
    if ctype == 2:
        return [math.sqrt(xl * xl + yl * yl), math.atan2(yl, xl) / _DEGREE, zl]
    r = math.sqrt(xl * xl + yl * yl + zl * zl)
    theta = 0.0
    if r > 0.0:
        theta = math.acos(min(1.0, max(-1.0, zl / r))) / _DEGREE
    return [r, theta, math.atan2(yl, xl) / _DEGREE]


def _local_points(mesh, points):
    """Each point's coordinates in its own system (basic for CP 0), each point's CP and
    the systems the CPs name, by id: what GEOM1 holds. Output systems (CD) are not
    written, so a non-zero CD is reported."""
    n = len(points)
    cp = np.zeros(n, dtype=np.int64)
    if "nastran:cp" in mesh.point_data:
        cp = np.asarray(mesh.point_data["nastran:cp"], dtype=np.int64).reshape(-1)
        if cp.size != n:
            _wfail("'nastran:cp' must hold one integer per point")
    if np.any(cp < 0):
        _wfail("coordinate system ids must be positive")
    cords = {}
    for cid in sorted(set(cp[cp != 0].tolist())):
        cords[cid] = _cord_system(mesh, cid)
    local = np.array(points, dtype=np.float64)
    for i in np.flatnonzero(cp):
        local[i] = _to_local(cords[int(cp[i])], points[i])
    if "nastran:cd" in mesh.point_data:
        named = int(np.count_nonzero(np.asarray(mesh.point_data["nastran:cd"]) > 0))
        if named:
            warn(
                f"Nastran OP2: {named} GRID(s) have an output system (nastran:cd) that is "
                "not written: their CD stays 0 and their results stay in basic coordinates"
            )
    return cp, local, cords


def _bush_values(mesh):
    """The PBUSH field data by id: ``{pid: {member: array}}`` from ``nastran:PBUSH_<id>:*``."""
    out = {}
    for name, value in mesh.field_data.items():
        match = re.fullmatch(r"nastran:PBUSH_([0-9]+):([A-Z]+)", name)
        if match:
            out.setdefault(int(match.group(1)), {})[match.group(2)] = np.asarray(value)
    return out


def _beam_values(mesh):
    """The PBEAM field data by id: ``{pid: {member: array}}`` from ``nastran:PBEAM_<id>:*``."""
    out = {}
    for name, value in mesh.field_data.items():
        match = re.fullmatch(r"nastran:PBEAM_([0-9]+):([A-Z0-9]+)", name)
        if match:
            out.setdefault(int(match.group(1)), {})[match.group(2)] = np.asarray(value)
    return out


class _Model:
    """The mesh as the file holds it: GRIDs, element cards and properties."""

    def __init__(self, mesh, per_card_ids=False):
        self.npoints = len(mesh.points)
        points = np.asarray(mesh.points, dtype=np.float64).reshape(self.npoints, -1)
        self.points = np.zeros((self.npoints, 3))
        self.points[:, : min(3, points.shape[1])] = points[:, :3]
        self.node_ids = _node_ids(mesh, self.npoints)
        self.cp, self.local, self.cords = _local_points(mesh, self.points)
        self.signature = _signature(mesh)
        self.values = property_values(mesh)
        self.beams = _beam_values(mesh)
        self.bushes = _bush_values(mesh)
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
            raw_bar = _raw_bar_fields(mesh, b, n)
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
                raw_bar = {
                    key: (None if value is None else value[keep])
                    for key, value in raw_bar.items()
                }
                n = len(eids)
            if not n:
                continue
            if np.any(cards == "CBAR") and not {"nastran:x", "nastran:g0"} & set(
                mesh.cell_data
            ):
                _wfail(
                    "a CBAR needs its orientation: nastran:x, or nastran:g0 for a G0"
                )
            for card, key in (("CELAS2", "nastran:k"), ("CDAMP2", "nastran:b")):
                if np.any(cards == card) and key not in mesh.cell_data:
                    _wfail(f"a {card} needs its value: {key}")
            if np.any(cards == "CBEAM"):
                if raw_bar["x"] is None and raw_bar["g0"] is None:
                    _wfail(
                        "a CBEAM needs its orientation: nastran:x, or nastran:g0 for a G0"
                    )
            for card in ("CBUSH", "CGAP"):
                if (
                    np.any(cards == card)
                    and raw_bar["x"] is None
                    and raw_bar["g0"] is None
                ):
                    _wfail(
                        f"a {card} needs its orientation: nastran:x, or nastran:g0 for a G0"
                    )
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
                    "bar": _bar_matrix(raw_bar, cards),
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
            "PBAR": set(),
            "PBEAM": set(),
            "PBUSH": set(),
            "PELAS": set(),
            "PDAMP": set(),
            "PVISC": set(),
            "PGAP": set(),
            "PSHEAR": set(),
            "PCOMP": set(),
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
                ("CBAR", "PBAR"),
                ("CBEAM", "PBEAM"),
                ("CBUSH", "PBUSH"),
                ("CELAS1", "PELAS"),
                ("CDAMP1", "PDAMP"),
                ("CVISC", "PVISC"),
                ("CGAP", "PGAP"),
                ("CSHEAR", "PSHEAR"),
            ):
                rows = blk["cards"] == name
                if rows.any():
                    self.properties[ptype].update(blk["pids"][rows].tolist())
        for pid in sorted(self.properties["PSHELL"]):
            if "PLY" in self.values.get(("PCOMP", pid), {}):
                self.properties["PSHELL"].discard(pid)
                self.properties["PCOMP"].add(pid)
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
        material_cards(self.values)  # refuses a material id given as two cards

    def element_rows(self, name):
        """(eids, pids, [(block, rows)], order) of a card; the first two by id."""
        parts, order = self.cards[name]
        eid = np.concatenate([blk["eids"][r] for blk, r in parts])[order]
        pid = np.concatenate([blk["pids"][r] for blk, r in parts])[order]
        return eid, pid, parts, order

    # -- geometry ---------------------------------------------------------------------

    def geom1(self):
        records = []
        # The coordinate systems the GRIDs are in: CORD2C, then CORD2R, then CORD2S, by id.
        for ctype, key in ((2, (2001, 20, 9)), (1, (2101, 21, 8)), (3, (2201, 22, 10))):
            ids = [cid for cid in sorted(self.cords) if self.cords[cid][0] == ctype]
            if not ids:
                continue
            systems = [self.cords[cid] for cid in ids]
            origin = np.array([s[1] for s in systems], dtype=np.float64)
            x_axis = np.array([s[2] for s in systems], dtype=np.float64)
            z_axis = np.array([s[4] for s in systems], dtype=np.float64)
            words = np.zeros((len(ids), 13), dtype="<u4")
            words[:, 0] = _i32(ids)
            words[:, 1] = ctype
            words[:, 2] = 2
            words[:, 3] = 0  # RID 0: the points are in basic coordinates
            words[:, 4:7] = _f32(origin)
            words[:, 7:10] = _f32(origin + z_axis)  # B: a point on +z
            words[:, 10:13] = _f32(origin + x_axis)  # C: a point in the xz plane
            records.append(_key(*key) + _bytes(words))
        n = self.npoints
        words = np.zeros((n, 8), dtype="<u4")
        order = np.argsort(self.node_ids, kind="stable")  # GEOM1 holds GRIDs by id
        words[:, 0] = _i32(self.node_ids[order])
        words[:, 1] = _i32(self.cp[order])  # CP: the system the coordinates are in
        words[:, 2:5] = _f32(self.local[order])
        records += [_key(4501, 45, 1) + _bytes(words), _key(65535, 65535, 65535)]
        return _table("GEOM1", records)

    def geom2(self):
        records = []
        for name in sorted(self.cards):
            key, size, linear, first, perm = _CARDS[name]
            parts, order = self.cards[name]
            eid = np.concatenate([blk["eids"][r] for blk, r in parts])
            pid = np.concatenate([blk["pids"][r] for blk, r in parts])
            comp = np.concatenate([blk["comps"][r] for blk, r in parts])
            bar = np.concatenate([blk["bar"][r] for blk, r in parts])
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
            eid, pid, nodes = eid[order], pid[order], nodes[order]
            comp, bar = comp[order], bar[order]
            words = np.zeros((len(eid), size), dtype="<u4")
            words[:, 0] = _i32(eid)
            if name == "CONM2":
                words[:, 1] = _i32(nodes[:, 0])
            elif name == "CGAP":
                # X1..X3 at 4 to 6, or the G0 at 4 (flag 2); the flag at 7 (1 for a
                # vector); the coordinate system CID (8) stays 0.
                g0 = bar[:, 3].astype(np.int64)
                flag = np.where(g0 != 0, 2, 1)
                words[:, 1] = _i32(pid)
                words[:, 2:4] = _i32(nodes)
                words[:, 4:7] = _f32(bar[:, 0:3])
                words[flag == 2, 4] = _i32(g0[flag == 2])
                words[flag == 2, 5:7] = 0
                words[:, 7] = _i32(flag)
            elif name == "PLOTEL":
                words[:, 1] = _i32(nodes[:, 0])  # no property: the nodes at 1 and 2
                words[:, 2] = _i32(nodes[:, 1])
            elif name == "CONROD":
                words[:, 1] = _i32(nodes[:, 0])
                words[:, 2] = _i32(nodes[:, 1])
                words[:, 3] = _i32(np.ones(len(eid)))  # MID
                words[:, 4] = _f32(np.ones(len(eid)))  # A
            else:
                words[:, 1] = _i32(pid)
                words[:, 2 : 2 + nodes.shape[1]] = _i32(nodes)
                if name in ("CELAS1", "CDAMP1", "CELAS2", "CDAMP2"):
                    words[:, 4:6] = _i32(comp)
                if name == "CELAS2":
                    # K at 1, GE at 6 and S at 7 (the reader's "bar" columns 0 to 2)
                    words[:, 1] = _f32(bar[:, 0])
                    words[:, 6] = _f32(bar[:, 1])
                    words[:, 7] = _f32(bar[:, 2])
                if name == "CDAMP2":
                    words[:, 1] = _f32(bar[:, 0])
                if name == "CBAR":
                    # X1..X3 at 4 to 6, or the G0 at 4 (flag 2); the flag at 7
                    g0 = bar[:, 3].astype(np.int64)
                    flag = np.where(g0 != 0, 2, 0)
                    words[:, 4:7] = _f32(bar[:, 0:3])
                    words[flag == 2, 4] = _i32(g0[flag == 2])
                    words[flag == 2, 5:7] = 0
                    words[:, 7] = _i32(flag)
                    words[:, 8:10] = _i32(bar[:, 4:6].astype(np.int64))
                    words[:, 10:16] = _f32(bar[:, 6:12])
                if name == "CBUSH":
                    # X1..X3 at 4 to 6, or the G0 at 4 (flag 2); the flag at 7; the
                    # element system 0 and the stress coefficient 0.5; output system -1.
                    g0 = bar[:, 3].astype(np.int64)
                    flag = np.where(g0 != 0, 2, 0)
                    words[:, 4:7] = _f32(bar[:, 0:3])
                    words[flag == 2, 4] = _i32(g0[flag == 2])
                    words[flag == 2, 5:7] = 0
                    words[:, 7] = _i32(flag)
                    words[:, 9] = _f32(np.full(len(eid), 0.5))
                    words[:, 10] = _i32(np.full(len(eid), -1, dtype=np.int64))
                if name == "CBEAM":
                    # X1..X3 at 6 to 8 (flag 1: global) or the G0 at 6 with zero 7 and 8
                    # (flag 2), the flag at 9, PA and PB at 10 and 11, the offsets at 12
                    # to 17; the scalar ends SA and SB stay zero.
                    g0 = bar[:, 3].astype(np.int64)
                    flag = np.where(g0 != 0, 2, 1)
                    words[:, 6:9] = _f32(bar[:, 0:3])
                    words[flag == 2, 6] = _i32(g0[flag == 2])
                    words[flag == 2, 7:9] = 0
                    words[:, 9] = _i32(flag)
                    words[:, 10:12] = _i32(bar[:, 4:6].astype(np.int64))
                    words[:, 12:18] = _f32(bar[:, 6:12])
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
            ("PBAR", _PBAR),
            ("PELAS", _PELAS),
            ("PDAMP", _PDAMP),
            ("PVISC", _PVISC),
            ("PGAP", _PGAP),
            ("PSHEAR", _PSHEAR),
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
        if self.properties["PCOMP"]:
            records.append(_key(*_PCOMP[0]) + _bytes(self._pcomp_words()))
        if self.properties["PBEAM"]:
            records.append(_key(*_PBEAM[0]) + _bytes(self._pbeam_words()))
        if self.properties["PBUSH"]:
            records.append(_key(*_PBUSH[0]) + _bytes(self._pbush_words()))
        records.append(_key(65535, 65535, 65535))
        return _table("EPT", records)

    def _pcomp_words(self):
        """The PCOMP records: the header words (pid, layer count, Z0, NSM, SB, FT, TREF,
        GE) and the plies (MID, T, THETA, SOUT), one entry per id."""
        out = []
        for pid in sorted(self.properties["PCOMP"]):
            given = self.values[("PCOMP", pid)]
            plies = given["PLY"]
            count = len(plies)
            head = np.zeros(8, dtype="<u4")
            head[0] = _i32(np.array([pid]))[0]
            head[1] = _i32(np.array([-count if given.get("SYM", 0) else count]))[0]
            for k, member in ((2, "Z0"), (3, "NSM"), (4, "SB"), (6, "TREF"), (7, "GE")):
                head[k] = _f32(np.array([given.get(member, 0.0)]))[0]
            head[5] = _i32(np.array([int(given.get("FT", 0))]))[0]
            body = np.zeros((count, 4), dtype="<u4")
            body[:, 0] = _i32(plies[:, 0].astype(np.int64))
            body[:, 1] = _f32(plies[:, 1])
            body[:, 2] = _f32(plies[:, 2])
            body[:, 3] = _i32(plies[:, 3].astype(np.int64))
            out.append(np.concatenate([head, body.ravel()]))
        return np.concatenate(out)

    def _pbeam_words(self):
        """The PBEAM records, one per id: the header (pid, MID, segments, CCF and a zero),
        then the stations in slots of 16 words, from the id's field data. An id with none
        gets placeholders: one segment of constant section, with A, I1, I2 and J of 1.
        """
        out = []
        for pid in sorted(self.properties["PBEAM"]):
            given = self.beams.get(pid, {})
            mid = int(np.asarray(given.get("MID", [1])).ravel()[0])
            nseg = int(np.asarray(given.get("NSEG", [1])).ravel()[0])
            ccf = int(np.asarray(given.get("CCF", [1])).ravel()[0])
            n = nseg + 1
            if not 1 <= nseg <= 10:
                _wfail(f"PBEAM {pid} has {nseg} segments; a PBEAM holds 1 to 10")
            defaults = {
                "SO": np.ones(n),
                "XXB": np.linspace(0.0, 1.0, n),
                "A": np.ones(n),
                "I1": np.ones(n),
                "I2": np.ones(n),
                "I12": np.zeros(n),
                "J": np.ones(n),
                "NSM": np.zeros(n),
            }
            stations = np.zeros((n, 16))
            for column, member in enumerate(_BEAM_STATION_MEMBERS):
                values = np.asarray(
                    given.get(member, defaults[member]), dtype=np.float64
                ).ravel()
                if values.size != n:
                    _wfail(f"'nastran:PBEAM_{pid}:{member}' must have {n} value(s)")
                stations[:, column] = values
            for first, member in zip((8, 10, 12, 14), _BEAM_POINT_MEMBERS):
                values = np.asarray(
                    given.get(member, np.zeros((n, 2))), dtype=np.float64
                ).reshape(-1, 2)
                if values.shape != (n, 2):
                    _wfail(f"'nastran:PBEAM_{pid}:{member}' must have {n} x 2 value(s)")
                stations[:, first : first + 2] = values
            words = np.zeros(_PBEAM[1], dtype="<u4")
            words[:4] = _i32([pid, mid, nseg, ccf])
            slots = np.zeros((11, 16))
            slots[pbeam_station_slots(n)] = stations
            words[5 : 5 + 16 * 11] = _f32(slots).ravel()
            out.append(words)
        return np.concatenate(out)

    def _pbush_words(self):
        """The PBUSH records, one per id: the id, K, B and GE (six each) and SA, ST, EA and
        ET, from the id's field data. An id with none gets K of 1, B of 0 and GE of 0, and
        stress and strain coefficients of 1. A GE given as one value is written to all six.
        """
        out = []
        for pid in sorted(self.properties["PBUSH"]):
            given = self.bushes.get(pid, {})
            words = np.zeros(_PBUSH[1], dtype="<u4")
            words[0] = _i32(np.array([pid]))[0]
            for member, first in _BUSH_ARRAYS:
                default = {"K": 1.0, "B": 0.0, "GE": 0.0}[member]
                values = np.asarray(
                    given.get(member, [default]), dtype=np.float64
                ).ravel()
                if values.size == 1:
                    values = np.full(6, values[0])
                if values.size != 6:
                    _wfail(f"'nastran:PBUSH_{pid}:{member}' must have 1 or 6 value(s)")
                words[first : first + 6] = _f32(values)
            for word, member in zip(range(19, 23), ("SA", "ST", "EA", "ET")):
                value = np.asarray(given.get(member, [1.0]), dtype=np.float64).ravel()
                if value.size != 1:
                    _wfail(f"'nastran:PBUSH_{pid}:{member}' must have 1 value")
                words[word] = _f32(value)[0]
            out.append(words)
        return np.concatenate(out)

    def mpt(self):
        """A MAT1 per material the properties name: unit modulus unless the field
        data gives its values. A material given as ``nastran:MAT2_<id>:*`` or
        ``nastran:MAT8_<id>:*`` is written as that card instead."""
        ids = material_ids(self.properties, self.values, self.beams) or [1]
        given_cards = material_cards(self.values)
        records = []
        for card, key, size in _MATERIALS:
            mids = [mid for mid in ids if given_cards.get(mid, "MAT1") == card]
            if not mids:
                continue
            words = np.zeros((len(mids), size), dtype="<u4")
            for k, mid in enumerate(mids):
                words[k, 0] = _i32(np.array([mid]))[0]
                given = self.values.get((card, mid), {})
                for member, (word, default) in PROPERTY_VALUES[card].items():
                    value = np.array([given.get(member, default)])
                    words[k, word] = (
                        _i32(value)[0] if is_int_member(member) else _f32(value)[0]
                    )
            records.append(_key(*key) + _bytes(words))
        records.append(_key(65535, 65535, 65535))
        return _table("MPT", records)


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

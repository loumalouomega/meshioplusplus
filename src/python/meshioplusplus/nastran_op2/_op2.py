"""I/O for Nastran OP2 result files (``.op2``), read-only.

The pure-Python twin of ``src/cpp/src/formats/nastran_op2.cpp``: both engines read
the same meshes.

An OP2 is a Fortran unformatted file (4-byte record markers; 4- or 8-byte words
in 64-bit runs) of named tables. Each table is its name, a header record, then
records separated by marker triples ``[-k, 1, 0]`` and closed by a ``0``. The
geometry tables hold card records keyed by three integers: ``GEOM1`` the GRIDs,
``GEOM2`` the elements, ``EPT`` the properties. A results table (``OUGV1``,
``OQG1``, ``OES1X1``, ``OSTR1X``...) alternates a 146-word header (table code,
subcase, mode or time, element type, width) with the data it describes. The
record layouts follow the MSC/NX DMAP manuals as pyNastran reads them. See
``doc/formats/nastran_op2.md``.
"""

import os

import numpy as np

from .._common import warn
from .._exceptions import ReadError
from .._fortran_records import fortran_records, sniff_fortran_records
from .._mesh import Mesh
from ..nastran._model import (
    CARDS,
    CoordCard,
    CoordSystems,
    add_cells,
    apply_frames,
    rotate_to_basic,
)

__all__ = ["read", "time_values"]

TIME_KEY = "meshio:time"
WHO = "Nastran OP2"
_NAN = float("nan")


def _fail(message):
    raise ReadError(f"{WHO}: {message}")


# -- framing ---------------------------------------------------------------------------


class _Stream:
    """The file's Fortran records, read as markers (one word) and data records."""

    def __init__(self, data):
        layout = sniff_fortran_records(data)
        if layout is None:
            _fail("the file is not Fortran unformatted (no record framing)")
        self.data = data
        self.blocks = fortran_records(data, layout, WHO)
        self.order = layout[1]
        first = self.blocks[0][1]
        if first not in (4, 8):
            _fail("the first record is not a one-word marker")
        self.ws = first
        self.itype = np.dtype(self.order + ("i4" if first == 4 else "i8"))
        self.ftype = np.dtype(self.order + ("f4" if first == 4 else "f8"))
        self.i = 0

    def done(self):
        return self.i >= len(self.blocks)

    def _marker_value(self, index):
        if index >= len(self.blocks):
            _fail("the file is truncated")
        off, n = self.blocks[index]
        if n != self.ws:
            _fail(f"expected a marker at byte {off - 4}, found a {n}-byte record")
        return int(np.frombuffer(self.data, self.itype, 1, off)[0])

    def peek(self):
        return self._marker_value(self.i)

    def marker(self):
        v = self._marker_value(self.i)
        self.i += 1
        return v

    def block(self):
        if self.done():
            _fail("the file is truncated")
        off, n = self.blocks[self.i]
        self.i += 1
        return self.data[off : off + n]

    def record(self):
        """A marker ``n > 0`` and its data, joined with any continuation blocks."""
        n = self.marker()
        if n <= 0:
            _fail(f"expected a record, found marker {n}")
        parts = [self.block()]
        while not self.done():
            off, size = self.blocks[self.i]
            if size != self.ws or self._marker_value(self.i) <= 0:
                break
            self.i += 1
            parts.append(self.block())
        return b"".join(parts)


def _text(raw, ws):
    """A Nastran string: 8-byte words hold 4 characters then 4 blanks in NX
    64-bit files (interlaced) and 8 characters in MSC ones."""
    if ws == 8 and len(raw) % 8 == 0 and len(raw) >= 16:
        halves = [raw[k + 4 : k + 8] for k in range(0, len(raw), 8)]
        if all(h == b"    " for h in halves) and raw[8:12] != b"    ":
            raw = b"".join(raw[k : k + 4] for k in range(0, len(raw), 8))
    return raw.decode("latin-1")


def _parse(data):
    """(version text, [(table name, [(mark, record bytes)])])."""
    s = _Stream(data)
    version = ""
    if s.peek() == 3:  # PARAM,POST,-1 header: date, tape code, version
        s.marker()
        s.block()
        s.marker()
        s.block()
        version = _text(s.record(), s.ws).strip()
        if s.marker() != -1 or s.marker() != 0:
            _fail("the file header is not closed by the markers -1, 0")
    tables = []
    while not s.done():
        m = s.peek()
        if m == 0:
            s.marker()
            continue
        if m < 0:
            _fail(f"expected a table name, found marker {m}")
        name = _text(s.record(), s.ws)[:8].strip()
        records = []
        mark = 0
        while True:
            if s.done():
                _fail(f"table {name} is truncated")
            m = s.peek()
            if m == 0:
                s.marker()
                break
            if m == -1:
                s.marker()
                mark = -1
                continue
            if m < 0:
                triple = (s.marker(), s.marker(), s.marker())
                if triple[1] != 1 or triple[2] not in (0, 1):
                    _fail(f"table {name}: malformed marker triple {triple}")
                mark = m
                continue
            records.append((mark, s.record()))
        tables.append((name, records))
    return s, version, tables


# -- geometry --------------------------------------------------------------------------

# (k1, k2) of a GEOM2 element record -> (card, candidate entry sizes in words,
# word index of the first node, node count, pid word or None). Where versions
# wrote different sizes (NX and MSC CQUAD4: 14 and 15 words) the first size
# whose entries validate (increasing ids, corners that are GRIDs) wins, in
# pyNastran's order.
_ELEMENTS = {
    (2408, 24): ("CBAR", (16,), 2, 2, 1),
    (5408, 54): ("CBEAM", (18,), 2, 2, 1),
    (2608, 26): ("CBUSH", (14,), 2, 2, 1),
    (1501, 15): ("CONM2", (13,), 1, 1, None),
    (1601, 16): ("CONROD", (8,), 1, 2, None),
    (4108, 41): ("CPENTA", (17,), 2, 15, 1),
    (17200, 172): ("CPYRAM", (16,), 2, 13, 1),
    (9108, 91): ("CQUAD", (11,), 2, 9, 1),
    (2958, 51): ("CQUAD4", (14, 15), 2, 4, 1),
    (4701, 47): ("CQUAD8", (17, 16, 18), 2, 8, 1),
    (8009, 80): ("CQUADR", (14, 15), 2, 4, 1),
    (3001, 30): ("CROD", (4,), 2, 2, 1),
    (3101, 31): ("CSHEAR", (6,), 2, 4, 1),
    (7308, 73): ("CHEXA", (22,), 2, 20, 1),
    (5508, 55): ("CTETRA", (12,), 2, 10, 1),
    (5959, 59): ("CTRIA3", (13, 14), 2, 3, 1),
    (4801, 48): ("CTRIA6", (13, 14, 15), 2, 6, 1),
    (9200, 92): ("CTRIAR", (13, 14), 2, 3, 1),
    (3701, 37): ("CTUBE", (4,), 2, 2, 1),
    (3901, 39): ("CVISC", (4,), 2, 2, 1),
    (5201, 52): ("PLOTEL", (3,), 1, 2, None),
    # NX 2019 and later write CQUAD4/CTRIA3 as these (pyNastran's CQUADRN/CTRIARN)
    (15401, 154): ("CQUAD4", (14, 15), 2, 4, 1),
    (15301, 153): ("CTRIA3", (13, 14), 2, 3, 1),
}

# Element records without a cell type, named in the warning.
_NAMED_ELEMENT_KEYS = {
    (601, 6): "CELAS1",
    (701, 7): "CELAS2",
    (801, 8): "CELAS3",
    (901, 9): "CELAS4",
    (201, 2): "CDAMP1",
    (301, 3): "CDAMP2",
    (401, 4): "CDAMP3",
    (501, 5): "CDAMP4",
    (1001, 10): "CMASS1",
    (1101, 11): "CMASS2",
    (1201, 12): "CMASS3",
    (1301, 13): "CMASS4",
    (1401, 14): "CONM1",
    (1908, 19): "CGAP",
    (5608, 56): "CBUSH1D",
    (6108, 61): "CTRIAX6",
    (9008, 90): "CQUADX",
    (10108, 101): "CTRIAX",
}

# (k1, k2) of an EPT record -> (card, entry size in words, or a walker name)
_PROPERTIES = {
    (2302, 23): ("PSHELL", 11),
    (2402, 24): ("PSOLID", 7),
    (902, 9): ("PROD", 6),
    (1002, 10): ("PSHEAR", 6),
    (1602, 16): ("PTUBE", 5),
    (52, 20): ("PBAR", 19),
    (5402, 54): ("PBEAM", 197),
    (302, 3): ("PELAS", 4),
    (1802, 18): ("PVISC", 3),
    (202, 2): ("PDAMP", 2),
    (402, 4): ("PMASS", 2),
    (2102, 21): ("PGAP", 11),
    (4706, 47): ("PLSOLID", 7),
    (4606, 46): ("PLPLANE", 11),
    (1402, 14): ("PBUSH", (18, 23, 24, 27)),
    (2706, 27): ("PCOMP", "pcomp"),
    (15006, 150): ("PCOMPG", "pcompg"),
    (9102, 91): ("PBARL", "terminated"),
    (9202, 92): ("PBEAML", "terminated"),
}


def _ints(s, raw):
    return np.frombuffer(raw, s.itype).astype(np.int64)


def _floats(s, raw):
    return np.frombuffer(raw, s.ftype).astype(np.float64)


def _geometry_records(s, tables, prefixes):
    """``(key, int words)`` of every card record of the tables named ``prefixes``."""
    out = []
    for name, records in tables:
        if not name.startswith(prefixes):
            continue
        for mark, raw in records:
            if mark > -3 or len(raw) < 3 * s.ws or len(raw) % s.ws:
                continue
            words = _ints(s, raw)
            out.append(((int(words[0]), int(words[1]), int(words[2])), raw))
    return out


def _valid_elements(rows, first, corners_count, grids):
    """Entries whose ids increase and whose corner nodes are GRIDs."""
    if rows.shape[0] == 0:
        return False
    eid = rows[:, 0]
    if np.any(eid <= 0) or np.any(np.diff(eid) <= 0):
        return False
    corners = rows[:, first : first + corners_count]
    if np.any(corners <= 0):
        return False
    if grids is not None:
        return bool(np.all(np.isin(corners, grids)))
    return True


def _grid_rows(s, raw, k3):
    """(ids, cp, xyz, cd) of a GRID record, or ``None``.

    Entries are ``id, cp, x, y, z, cd, ps, seid``: 8 words, the coordinates in
    the file's precision; or, in 32-bit files, 11 words with the coordinates as
    doubles (NX writes them under the key (4501, 45, 1120001)).
    """
    body = raw[3 * s.ws :]
    nwords = len(body) // s.ws
    layouts = [8, 11] if k3 != 1120001 else [11, 8]
    for size in layouts:
        if size == 11 and s.ws != 4:
            continue
        if nwords == 0 or nwords % size:
            continue
        n = nwords // size
        i = s.order + ("i4" if s.ws == 4 else "i8")
        if size == 8:
            f = s.order + ("f4" if s.ws == 4 else "f8")
            dt = np.dtype(
                [
                    ("id", i),
                    ("cp", i),
                    ("x", f, (3,)),
                    ("cd", i),
                    ("ps", i),
                    ("seid", i),
                ]
            )
        else:
            dt = np.dtype(
                [
                    ("id", i),
                    ("cp", i),
                    ("x", s.order + "f8", (3,)),
                    ("cd", i),
                    ("ps", i),
                    ("seid", i),
                ]
            )
        rows = np.frombuffer(body, dt, n)
        gid = rows["id"].astype(np.int64)
        if np.any(gid <= 0) or np.any(np.diff(gid) <= 0) or np.any(rows["cd"] < 0):
            continue
        return (
            gid.tolist(),
            rows["cp"].astype(np.int64).tolist(),
            rows["x"].astype(np.float64).tolist(),
            rows["cd"].astype(np.int64).tolist(),
        )
    return None


# CORD1C/R/S and CORD2C/R/S: key -> (type, defined by GRIDs)
_CORDS = {
    (1701, 17): (2, True),
    (1801, 18): (1, True),
    (1901, 19): (3, True),
    (2001, 20): (2, False),
    (2101, 21): (1, False),
    (2201, 22): (3, False),
}


def _cord_rows(s, raw, ctype, by_grids, out):
    """CORD1 entries ``cid, type, 1|2, g1, g2, g3`` (6 words); CORD2 entries
    ``cid, type, 2, rid, a1..c3`` (13 words in the file's precision) or, in
    32-bit files, 22 words with the coordinates as doubles (NX's GEOM1N)."""
    body = raw[3 * s.ws :]
    nwords = len(body) // s.ws
    for size in [6] if by_grids else [13, 22]:
        if (size == 22 and s.ws != 4) or nwords == 0 or nwords % size:
            continue
        n = nwords // size
        i = s.order + ("i4" if s.ws == 4 else "i8")
        if by_grids:
            rows = (
                np.frombuffer(body, np.dtype(i), n * 6).reshape(n, 6).astype(np.int64)
            )
            cards = [CoordCard(r[0], ctype, grids=r[3:6].tolist()) for r in rows]
            # The two flag words vary (pyNastran has seen 1 or 2 for a CORD1C).
            ok = all(r[0] > 0 and min(r[3:6]) > 0 for r in rows.tolist())
        else:
            f = s.order + (("f4" if s.ws == 4 else "f8") if size == 13 else "f8")
            dt = np.dtype([("h", i, (4,)), ("abc", f, (9,))])
            rows = np.frombuffer(body, dt, n)
            head = rows["h"].astype(np.int64).tolist()
            abc = rows["abc"].astype(np.float64).tolist()
            cards = [CoordCard(h[0], ctype, rid=h[3], abc=a) for h, a in zip(head, abc)]
            ok = all(h[0] > 0 and h[1] == ctype and h[3] >= 0 for h in head)
        if ok:
            out.extend(cards)
            return True
    return False


def _read_grids(s, tables):
    """GRID ids, coordinates, CP and CD, and the SPOINT/EPOINT ids."""
    ids, xyz, cp, cd = [], [], [], []
    spoints = set()
    cords = []
    for (k1, k2, _k3), raw in _geometry_records(s, tables, ("GEOM1",)):
        words = _ints(s, raw)[3:]
        if (k1, k2) == (4501, 45):
            grids = _grid_rows(s, raw, _k3)
            if grids is None:
                warn(
                    f"{WHO}: a GRID record of {len(words)} words matches no GRID layout; skipped"
                )
                continue
            gid, gcp, gxyz, gcd = grids
            ids.extend(gid)
            cp.extend(gcp)
            xyz.extend(gxyz)
            cd.extend(gcd)
            continue
        cord = _CORDS.get((k1, k2))
        if cord is not None and not _cord_rows(s, raw, *cord, cords):
            warn(
                f"{WHO}: a coordinate system record ({k1}, {k2}) of {len(words)} words "
                "matches no layout; skipped"
            )
    for (k1, k2, _k3), raw in _geometry_records(s, tables, ("GEOM2", "GEOM1")):
        if (k1, k2) in ((5551, 49), (707, 7)):  # SPOINT, EPOINT
            spoints.update(_ints(s, raw)[3:].tolist())
    return ids, xyz, cp, cd, spoints, cords


def _read_elements(s, tables, grids):
    cards = {}
    skipped = []
    grid_array = np.array(sorted(grids), dtype=np.int64) if grids else None
    for key, raw in _geometry_records(s, tables, ("GEOM2",)):
        k12 = key[:2]
        if k12 == (65535, 65535) or k12 in ((5551, 49), (707, 7)):
            continue
        spec = _ELEMENTS.get(k12)
        if spec is None:
            name = _NAMED_ELEMENT_KEYS.get(k12, f"({key[0]},{key[1]},{key[2]})")
            if name not in skipped:
                skipped.append(name)
            continue
        card, sizes, first, count, pid_word = spec
        words = _ints(s, raw)[3:]
        rows = None
        for size in sizes:
            if len(words) % size:
                continue
            candidate = words.reshape(-1, size)
            if _valid_elements(candidate, first, CARDS[card][1], grid_array):
                rows = candidate
                break
        if rows is None:
            warn(
                f"{WHO}: a {card} record of {len(words)} words matches no known layout; skipped"
            )
            continue
        entry = cards.setdefault(card, [[], [], []])
        entry[0].append(rows[:, 0])
        entry[1].append(
            rows[:, pid_word] if pid_word is not None else np.full(len(rows), -1)
        )
        entry[2].append(rows[:, first : first + count])
    if skipped:
        warn(f"{WHO}: skipped element records with no cell type: " + ", ".join(skipped))
    out = []
    for card in sorted(cards):
        eid, pid, g = cards[card]
        out.append((card, np.concatenate(eid), np.concatenate(pid), np.concatenate(g)))
    return out


def _entries_terminated(words, start):
    """End (exclusive) of an entry that runs to a ``-1`` word, or ``None``."""
    k = start
    while k < len(words):
        if words[k] == -1:
            return k + 1
        k += 1
    return None


# A pid named by several records keeps the first of these, else its first record.
_PROPERTY_PRIORITY = ("PCOMPG", "PCOMP")


def _read_properties(s, tables):
    """Property id -> property card name.

    NX also writes each PCOMP as a PSHELL whose material ids are 100000000 and
    up; those PSHELL entries are ignored, and a PCOMPG wins over a PCOMP of the
    same id, as pyNastran reads them.
    """
    found = {}
    for (k1, k2, _k3), raw in _geometry_records(s, tables, ("EPT",)):
        spec = _PROPERTIES.get((k1, k2))
        if spec is None:
            continue
        card, layout = spec
        words = _ints(s, raw)[3:]
        pids = []
        if isinstance(layout, int):
            if len(words) % layout == 0:
                rows = words.reshape(-1, layout)
                if card == "PSHELL":
                    rows = rows[rows[:, 1] < 100000000]
                pids = rows[:, 0].tolist()
        elif isinstance(layout, tuple):
            for size in layout:
                if len(words) % size == 0:
                    ids = words[::size]
                    if np.all(ids > 0) and np.all(np.diff(ids) > 0):
                        pids = ids.tolist()
                        break
        elif layout == "pcomp":
            k = 0
            while k + 8 <= len(words):
                nlayers = abs(int(words[k + 1]))
                if words[k] <= 0 or nlayers == 0 or k + 8 + 4 * nlayers > len(words):
                    break
                pids.append(int(words[k]))
                k += 8 + 4 * nlayers
        elif layout == "pcompg":
            k = 0
            while k + 8 <= len(words):
                if words[k] <= 0:
                    break
                pid = int(words[k])
                k += 8
                while k + 5 <= len(words) and not np.all(words[k : k + 5] == -1):
                    k += 5
                if k + 5 > len(words):
                    break
                pids.append(pid)
                k += 5
        else:  # an entry from its pid to a -1 word (PBARL, PBEAML)
            k = 0
            while k + 6 <= len(words):
                if words[k] <= 0 or (pids and words[k] <= pids[-1]):
                    break
                end = _entries_terminated(words, k + 6)
                if end is None:
                    break
                pids.append(int(words[k]))
                k = end
        for p in pids:
            found.setdefault(int(p), []).append(card)
    ptype = {}
    for pid, cards in found.items():
        chosen = next((c for c in _PROPERTY_PRIORITY if c in cards), cards[0])
        ptype[pid] = chosen
    return ptype


# -- results ---------------------------------------------------------------------------

_NODAL_PREFIXES = ("OUG", "BOUG", "OQG", "OQMG", "OPG")
_ELEMENT_PREFIXES = ("OES", "OSTR")

# table code -> point data name
_NODAL_NAMES = {
    1: "DISPLACEMENT",
    7: "EIGENVECTOR",
    10: "VELOCITY",
    11: "ACCELERATION",
    3: "SPC_FORCE",
    39: "MPC_FORCE",
    2: "APPLIED_LOAD",
}

_ELEMENT_TYPE_NAMES = {
    1: "CROD",
    2: "CBEAM",
    3: "CTUBE",
    4: "CSHEAR",
    10: "CONROD",
    11: "CELAS1",
    12: "CELAS2",
    13: "CELAS3",
    14: "CELAS4",
    33: "CQUAD4",
    34: "CBAR",
    39: "CTETRA",
    64: "CQUAD8",
    67: "CHEXA",
    68: "CPENTA",
    70: "CTRIAR",
    74: "CTRIA3",
    75: "CTRIA6",
    82: "CQUADR",
    95: "CQUAD4 composite",
    96: "CQUAD8 composite",
    97: "CTRIA3 composite",
    98: "CTRIA6 composite",
    100: "CBAR stations",
    102: "CBUSH",
    144: "CQUAD4 corner",
    255: "CPYRAM",
}

_PLATE_NAMES = ["FD1", "X1", "Y1", "TXY1", "ANGLE1", "MAJOR1", "MINOR1", None,
                "FD2", "X2", "Y2", "TXY2", "ANGLE2", "MAJOR2", "MINOR2", None]  # fmt: skip
_SOLID_NAMES = {  # word in the 21-word block -> name
    1: "X",
    2: "TXY",
    3: "PRINCIPAL_A",  # the principal values in the file's A, B, C order (unsorted)
    7: "PRESSURE",
    8: None,  # von Mises or octahedral shear, per s_code
    9: "Y",
    10: "TYZ",
    11: "PRINCIPAL_B",
    15: "Z",
    16: "TZX",
    17: "PRINCIPAL_C",
}
_BAR_NAMES = ["X1A", "X2A", "X3A", "X4A", "AX", "MAXA", "MINA", "MST",
              "X1B", "X2B", "X3B", "X4B", "MAXB", "MINB", "MSC"]  # fmt: skip
_SOLID_NODES = {39: 5, 67: 9, 68: 7, 255: 6}  # centre + corners
_CORNER_PLATE_NODES = {64: 5, 70: 4, 75: 4, 82: 5, 144: 5}


def _element_layout(etype, num_wide, s_code):
    """[(word, name)] of the element's centre values, or ``None``."""
    vm = "VON_MISES" if s_code & 1 else "MAX_SHEAR"
    if etype in (1, 10) and num_wide == 5:
        return [(1, "A"), (2, "MSA"), (3, "T"), (4, "MST")]
    if etype == 3 and num_wide == 5:
        return [(1, "AS"), (2, "MSA"), (3, "TS"), (4, "MST")]
    if etype == 4 and num_wide == 4:
        return [(1, "TMAX"), (2, "TAVG"), (3, "MS")]
    if etype == 34 and num_wide == 16:
        return [(1 + k, n) for k, n in enumerate(_BAR_NAMES)]
    if etype in (33, 74) and num_wide == 17:
        return [
            (1 + k, n if n else vm + str(1 + k // 8))
            for k, n in enumerate(_PLATE_NAMES)
        ]
    if etype in _CORNER_PLATE_NODES and num_wide == 2 + 17 * _CORNER_PLATE_NODES[etype]:
        return [
            (3 + k, n if n else vm + str(1 + k // 8))
            for k, n in enumerate(_PLATE_NAMES)
        ]
    if etype in _SOLID_NODES and num_wide == 4 + 21 * _SOLID_NODES[etype]:
        octa = "VON_MISES" if s_code & 1 else "OCT_SHEAR"
        return [(4 + k, n if n else octa) for k, n in sorted(_SOLID_NAMES.items())]
    return None


def _time_of(analysis, w5_int, w5_float, w6_float):
    if analysis in (5, 6, 10, 12):
        return w5_float
    if analysis in (2, 8, 9):
        return w6_float
    return 0.0


class _File:
    def __init__(self, filename):
        path = os.fspath(filename)
        with open(path, "rb") as f:
            data = f.read()
        if not data:
            _fail(f"'{path}' is empty")
        self.path = path
        self.stream, self.version, self.tables = _parse(data)
        self._scan_results()

    def _scan_results(self):
        """Every result block of a supported table: (step key, table, header, data)."""
        s = self.stream
        ws = s.ws
        self.steps = []  # dicts with key, subcase, analysis, mode, time
        step_index = {}
        self.blocks = []  # (step, kind, name info, data bytes)
        skipped = []

        def skip(reason):
            if reason not in skipped:
                skipped.append(reason)

        for name, records in self.tables:
            nodal = name.startswith(_NODAL_PREFIXES)
            elemental = name.startswith(_ELEMENT_PREFIXES)
            if not (nodal or elemental):
                if name.startswith("O"):
                    skip(name)
                continue
            header = None
            for mark, raw in records:
                if mark > -3:
                    continue
                if len(raw) == 146 * ws:
                    header = raw
                    continue
                if header is None:
                    continue
                h = _ints(s, header)
                f = _floats(s, header)
                approach, tcode, etype, subcase = (int(v) for v in h[:4])
                device = approach % 10
                analysis = (approach - device) // 10
                table_code = tcode % 1000
                sort_code = tcode // 1000
                num_wide = int(h[9])
                s_code = int(h[10])
                thermal = int(h[22])
                if sort_code & 1:
                    skip(f"{name} (complex)")
                    continue
                if sort_code & 4:
                    skip(f"{name} (random)")
                    continue
                if sort_code & 2:
                    skip(f"{name} (SORT2)")
                    continue
                if nodal:
                    base = _NODAL_NAMES.get(table_code)
                    # MPC forces share the SPC forces' table code; the name tells.
                    if name.startswith("OQMG") and table_code in (3, 39):
                        base = "MPC_FORCE"
                    if base is None:
                        skip(f"{name} (table code {table_code})")
                        continue
                    if num_wide != 8:
                        skip(f"{name} ({num_wide} words per node)")
                        continue
                    if thermal == 1:
                        if table_code != 1:
                            skip(f"{name} (thermal table code {table_code})")
                            continue
                        base = "TEMPERATURE"
                    info = ("nodal", base, name.startswith("BOUG"))
                else:
                    if table_code != 5:
                        skip(f"{name} (table code {table_code})")
                        continue
                    layout = _element_layout(etype, num_wide, s_code)
                    if layout is None:
                        ename = _ELEMENT_TYPE_NAMES.get(etype, f"type {etype}")
                        skip(f"{name} {ename}")
                        continue
                    group = "STRAIN" if s_code & 8 else "STRESS"
                    info = ("element", group, layout, num_wide)
                w5 = int(h[4])
                key = (subcase, analysis, w5)
                if key not in step_index:
                    step_index[key] = len(self.steps)
                    self.steps.append(
                        {
                            "subcase": subcase,
                            "analysis": analysis,
                            "mode": w5 if analysis in (2, 8, 9) else 0,
                            "time": _time_of(analysis, w5, float(f[4]), float(f[5])),
                        }
                    )
                self.blocks.append((step_index[key], info, raw))
        self.skipped = skipped


def time_values(filename):
    """The time of every step, in file order: frequency or time, the eigenvalue of a
    mode, 0 for a static subcase."""
    return [st["time"] for st in _File(filename).steps]


def _sibling_deck(path):
    stem, _ = os.path.splitext(path)
    tried = []
    for ext in (".bdf", ".dat", ".nas", ".blk", ".BDF", ".DAT", ".NAS", ".BLK"):
        candidate = stem + ext
        tried.append(candidate)
        if os.path.isfile(candidate):
            return candidate, tried
    return None, tried


def _mesh_from_deck(path):
    from ..nastran._nastran import read as read_deck

    mesh = read_deck(path)
    grid_index = {int(g): i for i, g in enumerate(mesh.points_id.tolist())}
    cell_index = {}
    offsets, sizes = [], []
    base = 0
    for block_ids in mesh.cells_id:
        offsets.append(base)
        sizes.append(len(block_ids))
        for i, e in enumerate(block_ids.tolist()):
            cell_index.setdefault(int(e), base + i)
        base += len(block_ids)
    if mesh.cells:
        mesh.cell_data["nastran:eid"] = [
            np.asarray(b, dtype=np.int64) for b in mesh.cells_id
        ]
    del mesh.points_id
    del mesh.cells_id
    return mesh, grid_index, cell_index, offsets, sizes


def _build_mesh(nf):
    s = nf.stream
    ids, xyz, cp, cd, spoints, cords = _read_grids(s, nf.tables)
    if not ids:
        deck, tried = _sibling_deck(nf.path)
        if deck is None:
            _fail(
                "the file has no GEOM1 GRID records (rerun with PARAM,POST,-1 or "
                "provide the input deck beside it); looked for " + ", ".join(tried)
            )
        return _mesh_from_deck(deck) + (CoordSystems(), None)
    grid_index = {}
    for i, g in enumerate(ids):
        if g in grid_index:
            _fail(f"GRID {g} is defined twice")
        grid_index[g] = i
    points = np.asarray(xyz, dtype=np.float64).reshape(len(ids), 3)
    point_data = {}
    systems = apply_frames(points, point_data, cords, ids, cp, cd, WHO)
    cards = _read_elements(s, nf.tables, set(ids))
    cells, cell_data, regions, cell_index, offsets, sizes = add_cells(
        cards, grid_index, spoints, _read_properties(s, nf.tables), WHO
    )
    mesh = Mesh(points, cells, point_data=point_data, cell_data=cell_data)
    mesh.regions = regions
    return mesh, grid_index, cell_index, offsets, sizes, systems, cd


def read(filename, points_only=False, arrays=None, time_step=0):
    nf = _File(filename)
    mesh, grid_index, cell_index, offsets, sizes, systems, cd = _build_mesh(nf)
    mesh.time_values = [st["time"] for st in nf.steps]
    n = len(nf.steps)
    if n == 0:
        if time_step not in (0, -1):
            raise ReadError(
                f"time step {time_step} is out of range: the file has no steps"
            )
        if nf.skipped:
            warn(f"{WHO}: result tables not read: " + ", ".join(nf.skipped))
        return mesh
    index = time_step + n if time_step < 0 else time_step
    if not 0 <= index < n:
        raise ReadError(
            f"time step {time_step} is out of range: the file has {n} step(s)"
        )
    step = nf.steps[index]
    mesh.field_data[TIME_KEY] = np.array([step["time"]], dtype=np.float64)
    for key in ("subcase", "analysis", "mode"):
        mesh.field_data["nastran:" + key] = np.array([step[key]], dtype=np.int64)
    if nf.skipped:
        warn(f"{WHO}: result tables not read: " + ", ".join(nf.skipped))
    if points_only:
        return mesh

    def wants(name):
        return arrays is None or name in arrays

    s = nf.stream
    npts = len(mesh.points)
    ncells = sum(sizes)
    point_arrays = {}
    cell_arrays = {}
    for step_id, info, raw in nf.blocks:
        if step_id != index:
            continue
        ints = _ints(s, raw)
        floats = _floats(s, raw)
        if info[0] == "nodal":
            base = info[1]
            if len(ints) % 8:
                _fail(f"a {base} record is not a whole number of rows")
            rows_i = ints.reshape(-1, 8)
            rows_f = floats.reshape(-1, 8)
            nid = rows_i[:, 0] // 10
            targets = [(r, grid_index.get(int(g))) for r, g in enumerate(nid.tolist())]
            outputs = (
                [(base, 2, 1)]
                if base == "TEMPERATURE"
                else [(base, 2, 3), (base + "_ROT", 5, 3)]
            )
            for name, c0, nc in outputs:
                if not wants(name):
                    continue
                values = point_arrays.get(name)
                if values is None:
                    values = point_arrays[name] = np.full((npts, nc), _NAN)
                for r, p in targets:
                    if p is not None and np.isnan(values[p, 0]):
                        v = rows_f[r, c0 : c0 + nc].tolist()
                        # Results are in the GRID's output system (CD) unless the table is BOUG*.
                        if nc == 3 and not info[2] and cd is not None:
                            rotate_to_basic(systems, cd[p], mesh.points[p].tolist(), v)
                        values[p] = v
        else:
            _, group, layout, num_wide = info
            if len(ints) % num_wide:
                _fail(f"a {group} record is not a whole number of elements")
            rows_i = ints.reshape(-1, num_wide)
            rows_f = floats.reshape(-1, num_wide)
            eids = (rows_i[:, 0] // 10).tolist()
            targets = [(r, cell_index.get(int(e))) for r, e in enumerate(eids)]
            for word, member in layout:
                name = f"{group}:{member}"
                if not wants(name):
                    continue
                values = cell_arrays.get(name)
                if values is None:
                    values = cell_arrays[name] = np.full(ncells, _NAN)
                for r, c in targets:
                    if c is not None and np.isnan(values[c]):
                        values[c] = rows_f[r, word]
    for name, values in point_arrays.items():
        mesh.point_data[name] = values[:, 0] if values.shape[1] == 1 else values
    for name, values in cell_arrays.items():
        mesh.cell_data[name] = [
            values[offsets[b] : offsets[b] + sizes[b]] for b in range(len(sizes))
        ]
    return mesh

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

import bisect
import math
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
    pbeam_station_slots,
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
    (201, 2): ("CDAMP1", (6,), 2, 2, 1),
    (301, 3): ("CDAMP2", (6,), 2, 2, None),
    (601, 6): ("CELAS1", (6,), 2, 2, 1),
    (701, 7): ("CELAS2", (8,), 2, 2, None),
    (1001, 10): ("CMASS1", (6,), 2, 2, 1),
    (1101, 11): ("CMASS2", (6,), 2, 2, None),
    (1401, 14): ("CONM1", (24,), 1, 1, None),
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
    (1908, 19): ("CGAP", (9,), 2, 2, 1),
    (5201, 52): ("PLOTEL", (3,), 1, 2, None),
    # NX 2019 and later write CQUAD4/CTRIA3 as these (pyNastran's CQUADRN/CTRIARN)
    (15401, 154): ("CQUAD4", (14, 15), 2, 4, 1),
    (15301, 153): ("CTRIA3", (13, 14), 2, 3, 1),
}

# Element records without a cell type, named in the warning.
_NAMED_ELEMENT_KEYS = {
    (801, 8): "CELAS3",
    (901, 9): "CELAS4",
    (401, 4): "CDAMP3",
    (501, 5): "CDAMP4",
    (1201, 12): "CMASS3",
    (1301, 13): "CMASS4",
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
        # A fluid (acoustic) GRID's CD is -1: it has no output system.
        if np.any(gid <= 0) or np.any(np.diff(gid) <= 0) or np.any(rows["cd"] < -1):
            continue
        return (
            gid.tolist(),
            rows["cp"].astype(np.int64).tolist(),
            rows["x"].astype(np.float64).tolist(),
            np.maximum(rows["cd"].astype(np.int64), 0).tolist(),
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


def _table_record(tables, prefix, mark):
    """The first record of mark ``mark`` of the first table whose name starts
    with ``prefix``, or ``None``."""
    for name, records in tables:
        if not name.startswith(prefix):
            continue
        for m, raw in records:
            if m == mark:
                return raw
        return None
    return None


def _bgpdt_grids(s, tables):
    """``(ids, xyz, cd, spoints)`` of a file without GEOM1, from its basic grid
    point table (BGPDT/BGPDTS), or ``None``.

    BGPDTS rows are ``cd, x, y, z`` in internal order, named by EQEXIN (pairs
    ``id, internal sequence`` and ``id, 10 * sil + type``: type 2 is a scalar
    point); NX's BGPDT rows carry the id: ``cd, sil, id, 61, ps, 0, x, y, z``
    with the coordinates as doubles.
    """
    ws = s.ws
    header = _table_record(tables, "BGPDT", -1)
    data = _table_record(tables, "BGPDT", -3)
    if header is None or data is None or len(header) < 2 * ws:
        return None
    n = int(_ints(s, header)[1])
    if n <= 0 or len(data) % n:
        return None
    row = len(data) // n
    ids = [0] * n
    scalar = [False] * n
    if row == 4 * ws:
        order = _table_record(tables, "EQEXIN", -3)
        kinds = _table_record(tables, "EQEXIN", -4)
        if order is None or kinds is None:
            return None
        seq = _ints(s, order).tolist()
        typ = _ints(s, kinds).tolist()
        type_of = {typ[i]: typ[i + 1] % 10 for i in range(0, len(typ) - 1, 2)}
        for i in range(0, len(seq) - 1, 2):
            k = seq[i + 1] - 1
            if not 0 <= k < n:
                return None
            ids[k] = seq[i]
            scalar[k] = type_of.get(seq[i]) == 2
    elif row != 12 * 4 and not (ws == 8 and row == 9 * 8):
        return None
    out_ids, xyz, cd, spoints = [], [], [], set()
    doubles = np.dtype(s.order + "f8")
    for k in range(n):
        e = data[k * row : (k + 1) * row]
        words = _ints(s, e)
        ident = ids[k] if row == 4 * ws else int(words[2])
        if ident <= 0:
            return None
        if scalar[k]:
            spoints.add(ident)
            continue
        out_ids.append(ident)
        cd.append(max(int(words[0]), 0))
        if row == 4 * ws:
            xyz.extend(float(v) for v in _floats(s, e)[1:4])
        else:
            xyz.extend(
                float(v) for v in np.frombuffer(e[6 * ws : 6 * ws + 24], doubles)
            )
    if not out_ids:
        return None
    return out_ids, xyz, cd, spoints


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
        # A spring or damper's ends may be scalar points, or 0 (grounded).
        scalar_ends = card.startswith(("CELAS", "CDAMP", "CMASS"))
        corners = 0 if scalar_ends else CARDS[card][1]
        rows = None
        for size in sizes:
            if len(words) % size:
                continue
            candidate = words.reshape(-1, size)
            if _valid_elements(candidate, first, corners, grid_array):
                rows = candidate
                break
        if rows is None:
            warn(
                f"{WHO}: a {card} record of {len(words)} words matches no known layout; skipped"
            )
            continue
        entry = cards.setdefault(card, [[], [], [], [], []])
        entry[0].append(rows[:, 0])
        entry[1].append(
            rows[:, pid_word] if pid_word is not None else np.full(len(rows), -1)
        )
        entry[2].append(rows[:, first : first + count])
        # A spring or damper's component numbers follow its two ends.
        if scalar_ends:
            entry[3].append(rows[:, first + count : first + count + 2])
        else:
            entry[3].append(np.zeros((len(rows), 2), dtype=np.int64))
        if card == "CBAR":
            entry[4].append(_bar_payload(s, raw, rows))
        elif card == "CBEAM":
            entry[4].append(_beam_payload(s, raw, rows))
        elif card == "CBUSH":
            entry[4].append(_bush_payload(s, raw, rows))
        elif card == "CGAP":
            entry[4].append(_gap_payload(s, raw, rows))
        elif card in ("CELAS2", "CDAMP2", "CMASS2"):
            # K, GE and S of a CELAS2 at words 1, 6 and 7; B of a CDAMP2 or M of a CMASS2
            # at word 1.
            floats = _floats(s, raw)[3:].reshape(rows.shape)
            columns = [1, 6, 7] if card == "CELAS2" else [1]
            entry[4].append(floats[:, columns])
        elif card == "CONM1":
            # The mass matrix M11, M21, M22, M31, ... M66 at words 3 to 23.
            floats = _floats(s, raw)[3:].reshape(rows.shape)
            if (rows[:, 2] > 0).any():
                warn(
                    f"{WHO}: {int((rows[:, 2] > 0).sum())} CONM1 element(s) have their "
                    "mass matrix in a coordinate system other than the basic one; it is "
                    "read as given"
                )
            entry[4].append(floats[:, 3:24])
    if skipped:
        warn(f"{WHO}: skipped element records with no cell type: " + ", ".join(skipped))
    out = []
    for card in sorted(cards):
        eid, pid, g, comps, payload = cards[card]
        out.append(
            (
                card,
                np.concatenate(eid),
                np.concatenate(pid),
                np.concatenate(g),
                np.concatenate(comps),
                np.concatenate(payload) if payload else None,
            )
        )
    return out


def _bar_payload(s, raw, rows):
    """The CBAR fields of each row of a CBAR record as float64 ``(n, 12)``: the
    orientation vector X1..X3 (zeros with a G0), the G0, the pins PA and PB, and the
    end offsets W1A..W3A and W1B..W3B. The record holds X1..X3 at words 4 to 6, the
    orientation flag at word 7 (2 with a G0, in word 4), PA and PB at 8 and 9, and
    the offsets at 10 to 15."""
    floats = _floats(s, raw)[3:].reshape(rows.shape)
    g0 = rows[:, 7] == 2
    x = np.where(g0[:, None], 0.0, floats[:, 4:7])
    return np.column_stack(
        [
            x,
            np.where(g0, rows[:, 4], 0).astype(np.float64),
            rows[:, 8:10].astype(np.float64),
            floats[:, 10:16],
        ]
    )


def _bush_payload(s, raw, rows):
    """The CBUSH fields of each row of a CBUSH record as float64 ``(n, 12)``, in the CBAR
    layout: the orientation vector X1..X3 at words 4 to 6 when the flag at word 7 is 0 or 1,
    and the G0 at word 4 when it is 2 (its vector then zeros). A flag of -1 orients the bush
    by its element coordinate system, which is not read, so its X is read as zeros."""
    floats = _floats(s, raw)[3:].reshape(rows.shape)
    flag = rows[:, 7]
    if np.any(flag == -1):
        warn(
            f"{WHO}: {int(np.count_nonzero(flag == -1))} CBUSH element(s) are oriented by "
            "their element coordinate system; their X vector is read as zeros"
        )
    vector = (flag == 0) | (flag == 1)
    out = np.zeros((len(rows), 12))
    out[:, 0:3] = np.where(vector[:, None], floats[:, 4:7], 0.0)
    out[:, 3] = np.where(flag == 2, rows[:, 4], 0)
    return out


def _gap_payload(s, raw, rows):
    """The CGAP fields of each row of a CGAP record as float64 ``(n, 12)``, in the CBAR
    layout: the orientation vector X1..X3 at words 4 to 6 when the flag at word 7 is 1, and
    the G0 at word 4 when it is 2 (its vector then zeros). The coordinate system at word 8
    is not read."""
    floats = _floats(s, raw)[3:].reshape(rows.shape)
    flag = rows[:, 7]
    out = np.zeros((len(rows), 12))
    out[:, 0:3] = np.where((flag != 2)[:, None], floats[:, 4:7], 0.0)
    out[:, 3] = np.where(flag == 2, rows[:, 4], 0)
    return out


def _beam_payload(s, raw, rows):
    """The CBEAM fields of each row of a CBEAM record as float64 ``(n, 12)``, in the CBAR
    layout: the orientation vector X1..X3 at words 6 to 8, or the G0 at word 6 (the flag at
    word 9, its low two bits: 0 and 1 are vectors in basic and global coordinates, 2 a G0),
    the pins PA and PB at 10 and 11, and the end offsets W1A..W3A and W1B..W3B at 12 to 17.
    A beam oriented by a local system has no X vector here, so its X is read as zeros.
    """
    floats = _floats(s, raw)[3:].reshape(rows.shape)
    low = rows[:, 9] & 3
    vector = low <= 1
    g0 = low == 2
    if (low == 3).any():
        warn(
            f"{WHO}: {int(np.count_nonzero(low == 3))} CBEAM element(s) are oriented by a "
            "local system; their X vector is read as zeros"
        )
    return np.column_stack(
        [
            np.where(vector[:, None], floats[:, 6:9], 0.0),
            np.where(g0, rows[:, 6], 0).astype(np.float64),
            rows[:, 10:12].astype(np.float64),
            floats[:, 12:18],
        ]
    )


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


def _read_property_values(s, tables):
    """``nastran:<CARD>_<id>:<MEMBER>`` field data of the PSHELL, PROD, PTUBE, PSOLID,
    PELAS, PDAMP, PVISC, PGAP, PSHEAR, MAT1, MAT2, MAT4, MAT5, MAT8 and MAT9 records (the members
    ``nastran_op2._write.PROPERTY_VALUES``
    lists)."""
    from ._write import PROPERTY_VALUES, is_int_member

    out = {}

    def take(card, raw, size):
        body = raw[3 * s.ws :]
        if len(body) == 0 or len(body) % (size * s.ws):
            return
        ints = _ints(s, body).reshape(-1, size)
        floats = _floats(s, body).reshape(-1, size)
        for row in range(len(ints)):
            ident = int(ints[row, 0])
            if card == "PSHELL" and ints[row, 1] >= 100000000:
                continue  # NX's PSHELL twin of a PCOMP
            for member, (word, _default) in PROPERTY_VALUES[card].items():
                if is_int_member(member):
                    value = np.array([ints[row, word]], dtype=np.int64)
                else:
                    value = np.array([floats[row, word]], dtype=np.float64)
                out.setdefault(f"nastran:{card}_{ident}:{member}", value)

    for (k1, k2, _k3), raw in _geometry_records(s, tables, ("EPT",)):
        spec = _PROPERTIES.get((k1, k2))
        if spec is not None and spec[0] in (
            "PSHELL",
            "PROD",
            "PTUBE",
            "PSOLID",
            "PBAR",
            "PELAS",
            "PDAMP",
            "PVISC",
            "PGAP",
            "PSHEAR",
            "PMASS",
        ):
            take(spec[0], raw, spec[1])
    for (k1, k2, _k3), raw in _geometry_records(s, tables, ("MPT", "MAT")):
        if (k1, k2) == (103, 1):
            take("MAT1", raw, 12)
        elif (k1, k2) == (203, 2):
            take("MAT2", raw, 17)
        elif (k1, k2) == (2103, 21):
            take("MAT4", raw, 11)
        elif (k1, k2) == (2203, 22):
            take("MAT5", raw, 10)
        elif (k1, k2) == (2503, 25):
            take("MAT8", raw, 19)
        elif (k1, k2) == (2603, 26):
            take("MAT9", raw, 35)
    return out


def _read_pbeam_values(s, tables):
    """``nastran:PBEAM_<id>:<MEMBER>`` field data of the PBEAM records: MID, NSEG (the
    segments, one fewer than the stations) and CCF as scalars, and for each station the SO,
    XXB, A, I1, I2, I12, J and NSM values, and the outline points C1, C2, D1 and D2 as
    ``(n, 2)`` rows of x and y. A record is pid, mid, segments, ccf and a zero, then 11
    station slots and a 12th slot of 16 words (SO, XXB, A, I1, I2, I12, J, NSM, then the
    four points). The stations used are segments + 1: the first in slot 0, the last in
    slot 10, and any between them in the slots after the first."""
    out = {}
    for (k1, k2, _k3), raw in _geometry_records(s, tables, ("EPT",)):
        if (k1, k2) != (5402, 54):
            continue
        body = raw[3 * s.ws :]
        if len(body) != 197 * s.ws:
            continue
        ints = _ints(s, body)
        floats = _floats(s, body)
        pid, mid, nseg, ccf = (int(v) for v in ints[:4])
        if not 1 <= nseg <= 10:
            continue
        n = nseg + 1
        prefix = f"nastran:PBEAM_{pid}:"
        if prefix + "MID" in out:  # the first record of an id wins
            continue
        out[prefix + "MID"] = np.array([mid], dtype=np.int64)
        out[prefix + "NSEG"] = np.array([nseg], dtype=np.int64)
        out[prefix + "CCF"] = np.array([ccf], dtype=np.int64)
        slots = floats[5 : 5 + 16 * 11].reshape(11, 16)
        stations = slots[pbeam_station_slots(n)]
        for column, member in enumerate(
            ("SO", "XXB", "A", "I1", "I2", "I12", "J", "NSM")
        ):
            out[prefix + member] = stations[:, column].copy()
        for first, member in zip((8, 10, 12, 14), ("C1", "C2", "D1", "D2")):
            out[prefix + member] = stations[:, first : first + 2].copy()
    return out


def _read_bush_values(s, tables):
    """``nastran:PBUSH_<id>:<MEMBER>`` field data of the PBUSH records: K and B as (6,)
    arrays, GE as one value (the 18-word NX layout) or six (the 23-word MSC layout), and the
    stress and strain coefficients SA, ST, EA and ET as scalars. A record of another size is
    not read."""
    out = {}
    for (k1, k2, _k3), raw in _geometry_records(s, tables, ("EPT",)):
        if (k1, k2) != (1402, 14):
            continue
        body = raw[3 * s.ws :]
        if len(body) not in (18 * s.ws, 23 * s.ws):
            continue
        floats = _floats(s, body)
        prefix = f"nastran:PBUSH_{_ints(s, body)[0]}:"
        if prefix + "K" in out:  # the first record of an id wins
            continue
        out[prefix + "K"] = floats[1:7].copy()
        out[prefix + "B"] = floats[7:13].copy()
        out[prefix + "GE"] = floats[13 : len(floats) - 4].copy()
        for member, value in zip(("SA", "ST", "EA", "ET"), floats[-4:]):
            out[prefix + member] = np.array([value])
    return out


def _read_pcomp_values(s, tables):
    """``nastran:PCOMP_<id>:<MEMBER>`` field data of the PCOMP records: the scalars Z0,
    NSM, SB, TREF and GE, the failure theory FT and SYM (1 for a symmetric layup), and
    PLY, the plies as ``(n, 4)`` rows of MID, T, THETA and SOUT. A record's words are
    pid, layer count (negative for a symmetric layup), Z0, NSM, SB, FT, TREF, GE, then
    MID, T, THETA and SOUT for each ply."""
    out = {}
    for (k1, k2, _k3), raw in _geometry_records(s, tables, ("EPT",)):
        if (k1, k2) != (2706, 27):
            continue
        body = raw[3 * s.ws :]
        if len(body) == 0 or len(body) % s.ws:
            continue
        ints = _ints(s, body)
        floats = _floats(s, body)
        k = 0
        while k + 8 <= len(ints):
            pid, nlayers = int(ints[k]), int(ints[k + 1])
            count = abs(nlayers)
            if pid <= 0 or count == 0 or k + 8 + 4 * count > len(ints):
                break
            prefix = f"nastran:PCOMP_{pid}"
            scalars = {
                "Z0": floats[k + 2],
                "NSM": floats[k + 3],
                "SB": floats[k + 4],
                "FT": ints[k + 5],
                "TREF": floats[k + 6],
                "GE": floats[k + 7],
                "SYM": 1 if nlayers < 0 else 0,
            }
            for member, value in scalars.items():
                dtype = np.int64 if member in ("FT", "SYM") else np.float64
                out.setdefault(f"{prefix}:{member}", np.array([value], dtype=dtype))
            start, stop = k + 8, k + 8 + 4 * count
            plies = np.column_stack(
                [
                    ints[start:stop:4],
                    floats[start + 1 : stop : 4],
                    floats[start + 2 : stop : 4],
                    ints[start + 3 : stop : 4],
                ]
            ).astype(np.float64)
            out.setdefault(f"{prefix}:PLY", plies)
            k = stop
    return out


def _read_line_values(s, tables):
    """``nastran:<CARD>_<id>:<MEMBER>`` field data of the PBARL, PBEAML and PCOMPG records.

    A PBARL or PBEAML is the id, MID, GROUP and TYPE (two words each, text) and the
    values up to a -1: a PBARL's dimensions then NSM, a PBEAML's stations as SO, X/XB,
    the dimensions and NSM. They give ``MID`` (int64), ``GROUP`` and ``TYPE`` (``uint8``
    text) and ``DIM`` (``(n,)``, or ``(stations, n)``), ``NSM``, and a PBEAML's ``SO`` and
    ``XXB`` (``(stations,)``). A PCOMPG is pid, LAMOPT, Z0, NSM, SB, FT, TREF and GE, then
    the plies as GPLYID, MID, T, THETA and SOUT up to five -1: ``LAMOPT`` and ``FT``
    (int64), the scalars, and ``PLY``, ``(n, 5)`` rows. A file of 8-byte words gives none
    (its texts are four characters a word)."""
    from ._write import SECTION_DIMS

    out = {}
    if s.ws != 4:
        return out
    for (k1, k2, _k3), raw in _geometry_records(s, tables, ("EPT",)):
        if (k1, k2) not in ((9102, 91), (9202, 92), (15006, 150)):
            continue
        body = raw[3 * s.ws :]
        if len(body) == 0 or len(body) % s.ws:
            continue
        ints = _ints(s, body)
        floats = _floats(s, body)
        k = 0
        if (k1, k2) == (15006, 150):
            while k + 8 <= len(ints) and ints[k] > 0:
                pid = int(ints[k])
                start = k + 8
                at = start
                while at + 5 <= len(ints) and not np.all(ints[at : at + 5] == -1):
                    at += 5
                if at + 5 > len(ints):
                    break
                prefix = f"nastran:PCOMPG_{pid}:"
                if prefix + "PLY" not in out:
                    plies = np.column_stack(
                        [
                            ints[start:at:5],
                            ints[start + 1 : at : 5],
                            floats[start + 2 : at : 5],
                            floats[start + 3 : at : 5],
                            ints[start + 4 : at : 5],
                        ]
                    ).astype(np.float64)
                    out[prefix + "PLY"] = plies
                    out[prefix + "LAMOPT"] = np.array([ints[k + 1]], dtype=np.int64)
                    out[prefix + "FT"] = np.array([ints[k + 5]], dtype=np.int64)
                    for member, word in (
                        ("Z0", 2),
                        ("NSM", 3),
                        ("SB", 4),
                        ("TREF", 6),
                        ("GE", 7),
                    ):
                        out[prefix + member] = np.array([floats[k + word]])
                k = at + 5
            continue
        card = "PBARL" if k1 == 9102 else "PBEAML"
        while k + 6 <= len(ints) and ints[k] > 0:
            end = _entries_terminated(ints, k + 6)
            if end is None:
                break
            pid, mid = int(ints[k]), int(ints[k + 1])
            prefix = f"nastran:{card}_{pid}:"
            values = floats[k + 6 : end - 1]
            texts = [
                body[(k + 2 + 2 * j) * s.ws : (k + 4 + 2 * j) * s.ws]
                .decode("latin-1")
                .strip(" \x00")
                for j in range(2)
            ]
            section = texts[1].upper()
            k = end
            if prefix + "TYPE" in out:
                continue
            fields = {}
            if card == "PBARL":
                if len(values) < 1:
                    continue
                fields["DIM"] = values[:-1].copy()
                fields["NSM"] = np.array([values[-1]])
            else:
                ndim = SECTION_DIMS.get(section)
                if ndim is None or len(values) == 0 or len(values) % (ndim + 3):
                    continue
                stations = values.reshape(-1, ndim + 3)
                fields["SO"] = stations[:, 0].copy()
                fields["XXB"] = stations[:, 1].copy()
                fields["DIM"] = stations[:, 2 : 2 + ndim].copy()
                fields["NSM"] = stations[:, -1].copy()
            out[prefix + "MID"] = np.array([mid], dtype=np.int64)
            out[prefix + "GROUP"] = np.frombuffer(
                texts[0].encode("latin-1"), np.uint8
            ).copy()
            out[prefix + "TYPE"] = np.frombuffer(
                section.encode("latin-1"), np.uint8
            ).copy()
            for member, value in fields.items():
                out[prefix + member] = value
    return out


# -- results ---------------------------------------------------------------------------

_NODAL_PREFIXES = ("OUG", "BOUG", "OQG", "OQMG", "OPG")
_ELEMENT_PREFIXES = ("OES", "OSTR", "OEF", "HOEF", "ONRGY", "OEKE")

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
    20: "CDAMP1",
    21: "CDAMP2",
    22: "CDAMP3",
    23: "CDAMP4",
    24: "CVISC",
    33: "CQUAD4",
    34: "CBAR",
    38: "CGAP",
    39: "CTETRA",
    53: "CTRIAX6",
    64: "CQUAD8",
    67: "CHEXA",
    68: "CPENTA",
    69: "CBEND",
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
    107: "CHBDYE",
    108: "CHBDYG",
    109: "CHBDYP",
    110: "CONV",
    126: "CFAST",
    144: "CQUAD4 corner",
    227: "CTRIAR",
    228: "CQUADR",
    232: "CQUADR composite",
    233: "CTRIAR composite",
    255: "CPYRAM",
    280: "CBEAR",
    300: "CHEXA",
    301: "CPENTA",
    302: "CTETRA",
    303: "CPYRAM",
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
_NX_SOLID_CORNERS = {300: 8, 301: 6, 302: 4, 303: 5}  # no centre block
_SHELL_FORCE = ["MX", "MY", "MXY", "BMX", "BMY", "BMXY", "TX", "TY"]
_BUSH_FORCE = ["FX", "FY", "FZ", "MX", "MY", "MZ"]
_TENSOR = ["X", "Y", "Z", "TXY", "TYZ", "TZX"]
_PLY = ["X1", "Y1", "T1", "L1", "L2"]
_BUSH_STRESS = ["TX", "TY", "TZ", "RX", "RY", "RZ"]
_FLUX = ["XGRAD", "YGRAD", "ZGRAD", "XFLUX", "YFLUX", "ZFLUX"]


# A member is (word, imaginary part's word or None, name): the word within the
# entry, or within a node block (whose GRID is word 0).
def _real(names, first):
    """Real members ``names`` at words first, first + 1 ..."""
    return [(first + k, None, n) for k, n in enumerate(names)]


def _block(names, first):
    """Complex members: all real parts from ``first``, then all imaginary parts."""
    return [(first + k, first + len(names) + k, n) for k, n in enumerate(names)]


def _pairs(names, first):
    """Complex members as (real, imaginary) pairs from ``first``."""
    return [(first + 2 * k, first + 2 * k + 1, n) for k, n in enumerate(names)]


def _element_layout(family, etype, num_wide, s_code, complex_=False, random=False):
    """How an element table's entries hold their values (``op2_element_layout``):
    ``kind`` row (members at fixed words), ply (one entry per ply, word 1 the
    ply), station (one entry per CBAR station, in order) or blocks (node blocks
    led by their GRID: the centre then the corners, or a beam's stations);
    ``None`` when the table is not read. ``family`` is ``stress`` (OES, OSTR),
    ``force`` (OEF), ``flux`` (a heat-transfer OEF) or ``energy`` (ONRGY,
    OEKE)."""
    vm = "VON_MISES" if s_code & 1 else "MAX_SHEAR"

    def layout(
        members, kind="row", first=0, block=0, blocks=1, stations=False, centre=True
    ):
        return {
            "kind": kind,
            "members": members,
            "first": first,
            "block": block,
            "blocks": blocks,
            "stations": stations,
            "centre": centre,
        }

    cn = _CORNER_PLATE_NODES.get(etype, 0)
    sn = _SOLID_NODES.get(etype, 0)
    if family == "energy":
        if not complex_ and num_wide == 4:
            return layout(_real(["ENERGY", "PCT", "DEN"], 1))
        if complex_ and num_wide == 5:
            return layout([(1, 2, "ENERGY")] + _real(["PCT", "DEN"], 3))
        return None
    if family == "flux":
        if etype in (1, 2, 3, 10, 34, 69, 33, 53, 64, 74, 75) and num_wide == 9:
            return layout(_real(_FLUX, 3))
        if etype in (39, 67, 68) and num_wide == 10:
            return layout(_real(_FLUX, 3))
        if etype in (107, 108, 109) and num_wide == 8:
            names = ["FAPPLIED", "FREECONV", "FORCECONV", "FRAD", "FTOTAL"]
            return layout(_real(names, 3))
        if etype == 110 and num_wide == 4:
            return layout([(1, None, "FREECONV"), (3, None, "FREECONVK")])
        return None
    if family == "force":
        # Random force tables use the real layouts.
        def real_or_complex(names, real_wide):
            if not complex_ and num_wide == real_wide:
                return layout(_real(names, 1))
            if complex_ and num_wide == 1 + 2 * len(names):
                return layout(_block(names, 1))
            return None

        if etype in (1, 3, 10, 24):
            return real_or_complex(["AF", "TRQ"], 3)
        if etype in (11, 12, 13, 14, 20, 21, 22, 23):
            return real_or_complex(["F"], 2)
        if etype == 34:
            bar = ["BM1A", "BM2A", "BM1B", "BM2B", "TS1", "TS2", "AF", "TRQ"]
            return real_or_complex(bar, 9)
        if etype == 100 and not complex_ and num_wide == 8:
            names = ["SD", "BM1", "BM2", "TS1", "TS2", "AF", "TRQ"]
            return layout(_real(names, 1), "station")
        if etype == 4:
            forces = ["F41", "F21", "F12", "F32", "F23", "F43", "F34", "F14"]
            kicks = ["KF1", "S12", "KF2", "S23", "KF3", "S34", "KF4", "S41"]
            if not complex_ and num_wide == 17:
                return layout(_real(forces, 1) + _real(kicks, 9))
            # The forces' real then imaginary parts, then the kick forces' and
            # shear flows' (pyNastran pairs words 1-16 with 17-32 instead).
            if complex_ and num_wide == 33:
                return layout(_block(forces, 1) + _block(kicks, 17))
            return None
        if etype in (33, 74, 227, 228):
            return real_or_complex(_SHELL_FORCE, 9)
        if etype in (102, 126, 280):
            return real_or_complex(_BUSH_FORCE, 7)
        if etype == 38 and not complex_ and num_wide == 9:
            return layout(_real(["FX", "SFY", "SFZ", "U", "V", "W", "SV", "SW"], 1))
        if etype == 2:
            beam = ["BM1", "BM2", "TS1", "TS2", "AF", "TTRQ", "WTRQ"]
            sd = [(1, None, "SD")]
            if not complex_ and num_wide == 100:
                return layout(sd + _real(beam, 2), "blocks", 1, 9, 11, True)
            if complex_ and num_wide == 177:
                return layout(sd + _block(beam, 2), "blocks", 1, 16, 11, True)
            return None
        if cn and not complex_ and num_wide == 2 + 9 * cn:
            return layout(_real(_SHELL_FORCE, 1), "blocks", 2, 9, cn)
        if cn and complex_ and num_wide == 2 + 17 * cn:
            return layout(_block(_SHELL_FORCE, 1), "blocks", 2, 17, cn)
        return None

    # Stresses and strains.
    plate_c = (
        [(1, None, "FD1")]
        + _pairs(["X1", "Y1", "TXY1"], 2)
        + [(8, None, "FD2")]
        + _pairs(["X2", "Y2", "TXY2"], 9)
    )
    plate_cvm = (
        [(1, None, "FD1")]
        + _pairs(["X1", "Y1", "TXY1"], 2)
        + [(8, None, "VON_MISES1"), (9, None, "FD2")]
        + _pairs(["X2", "Y2", "TXY2"], 10)
        + [(16, None, "VON_MISES2")]
    )
    if random:
        # Magnitudes only (PSD, RMS ...): pyNastran's random layouts.
        rp = ["FD1", "X1", "Y1", "TXY1", "FD2", "X2", "Y2", "TXY2"]
        rpvm = ["FD1", "X1", "Y1", "TXY1", "VON_MISES1"]
        rpvm += ["FD2", "X2", "Y2", "TXY2", "VON_MISES2"]
        if etype in (1, 10) and num_wide == 3:
            return layout(_real(["A", "T"], 1))
        if etype == 3 and num_wide == 3:
            return layout(_real(["AS", "TS"], 1))
        if etype == 4 and num_wide == 3:
            return layout(_real(["TMAX", "TAVG"], 1))
        if etype == 34 and num_wide == 10:
            names = ["X1A", "X2A", "X3A", "X4A", "AX", "X1B", "X2B", "X3B", "X4B"]
            return layout(_real(names, 1))
        if etype == 2 and num_wide == 67:
            names = ["SD", "XC", "XD", "XE", "XF"]
            return layout(_real(names, 1), "blocks", 1, 6, 11, True)
        if etype in (33, 74, 227, 228) and num_wide == 9:
            return layout(_real(rp, 1))
        if etype in (33, 74, 227, 228) and num_wide == 11:
            return layout(_real(rpvm, 1))
        if cn and num_wide == 2 + 9 * cn:
            return layout(_real(rp, 1), "blocks", 2, 9, cn)
        if cn and num_wide == 2 + 11 * cn:
            return layout(_real(rpvm, 1), "blocks", 2, 11, cn)
        if sn and num_wide == 4 + 7 * sn:
            return layout(_real(_TENSOR, 1), "blocks", 4, 7, sn)
        if sn and num_wide == 4 + 8 * sn:
            members = _real(_TENSOR, 1) + [(7, None, "VON_MISES")]
            return layout(members, "blocks", 4, 8, sn)
        if etype in (95, 96, 97, 98, 232, 233) and num_wide in (7, 8):
            members = _real(_PLY, 2)
            if num_wide == 8:
                members.append((7, None, "VON_MISES"))
            return layout(members, "ply")
        if etype == 102 and num_wide == 7:
            return layout(_real(_BUSH_STRESS, 1))
        return None
    if complex_:
        if etype in (1, 10) and num_wide == 5:
            return layout(_pairs(["A", "T"], 1))
        if etype == 3 and num_wide == 5:
            return layout(_pairs(["AS", "TS"], 1))
        if etype in (11, 12, 13, 14) and num_wide == 3:
            return layout([(1, 2, "S")])
        if etype == 4 and num_wide == 5:
            return layout(_pairs(["TMAX", "TAVG"], 1))
        if etype == 34 and num_wide == 19:
            return layout(
                _block(["X1A", "X2A", "X3A", "X4A", "AX"], 1)
                + _block(["X1B", "X2B", "X3B", "X4B"], 11)
            )
        if etype == 2 and num_wide == 111:
            members = [(1, None, "SD")] + _block(["XC", "XD", "XE", "XF"], 2)
            return layout(members, "blocks", 1, 10, 11, True)
        if etype in (33, 74, 227, 228) and num_wide == 15:
            return layout(plate_c)
        if etype in (33, 74, 227, 228) and num_wide == 17:
            return layout(plate_cvm)
        if cn and num_wide == 2 + 15 * cn:
            return layout(plate_c, "blocks", 2, 15, cn)
        if cn and num_wide == 2 + 17 * cn:
            return layout(plate_cvm, "blocks", 2, 17, cn)
        if sn and num_wide == 4 + 13 * sn:
            return layout(_block(_TENSOR, 1), "blocks", 4, 13, sn)
        if sn and num_wide == 4 + 14 * sn:
            members = _block(_TENSOR, 1) + [(13, None, "VON_MISES")]
            return layout(members, "blocks", 4, 14, sn)
        if etype in (95, 96, 97, 98, 232, 233) and num_wide in (12, 13):
            members = _block(_PLY, 2)
            if num_wide == 13:
                members.append((12, None, "VON_MISES"))
            return layout(members, "ply")
        if etype == 102 and num_wide == 13:
            return layout(_block(_BUSH_STRESS, 1))
        return None
    plate = [
        (1 + k, None, n if n else vm + str(1 + k // 8))
        for k, n in enumerate(_PLATE_NAMES)
    ]
    if etype in (1, 10) and num_wide == 5:
        return layout(_real(["A", "MSA", "T", "MST"], 1))
    if etype == 3 and num_wide == 5:
        return layout(_real(["AS", "MSA", "TS", "MST"], 1))
    if etype == 4 and num_wide == 4:
        return layout(_real(["TMAX", "TAVG", "MS"], 1))
    if etype in (11, 12, 13, 14) and num_wide == 2:
        return layout(_real(["S"], 1))
    if etype == 102 and num_wide == 7:
        return layout(_real(_BUSH_STRESS, 1))
    if etype == 34 and num_wide == 16:
        return layout(_real(_BAR_NAMES, 1))
    if etype == 100 and num_wide == 10:  # CBAR stations
        names = ["SD", "XC", "XD", "XE", "XF", "AX", "MAX", "MIN", "MS"]
        return layout(_real(names, 1), "station")
    if etype in (95, 96, 97, 98, 232, 233) and num_wide == 11:  # composite plies
        names = ["X1", "Y1", "T1", "L1", "L2", "ANGLE", "MAJOR", "MINOR", vm]
        return layout(_real(names, 2), "ply")
    if etype == 2 and num_wide == 111:  # CBEAM: 11 stations
        names = ["SD", "XC", "XD", "XE", "XF", "MAX", "MIN", "MST", "MSC"]
        return layout(_real(names, 1), "blocks", 1, 10, 11, stations=True)
    if etype in (33, 74) and num_wide == 17:
        return layout(plate)
    if cn and num_wide == 2 + 17 * cn:
        return layout(plate, "blocks", 2, 17, cn)
    if sn and num_wide == 4 + 21 * sn:
        octa = "VON_MISES" if s_code & 1 else "OCT_SHEAR"
        members = [(k, None, n if n else octa) for k, n in sorted(_SOLID_NAMES.items())]
        return layout(members, "blocks", 4, 21, sn)
    # NX's newer solids (CHEXA 300, CPENTA 301, CTETRA 302, CPYRAM 303): the
    # corners only, no centre.
    nx = _NX_SOLID_CORNERS.get(etype, 0)
    if nx and num_wide == 3 + 8 * nx:
        members = _real(_TENSOR, 1) + [(7, None, "VON_MISES")]
        return layout(members, "blocks", 3, 8, nx, centre=False)
    return None


_DEGREE = 3.14159265358979323846 / 180.0


def _complex(a, b, mag_phase):
    """A complex value's real and imaginary parts from its two words:
    themselves, or a magnitude and a phase in degrees."""
    if not mag_phase:
        return a, b
    return a * math.cos(b * _DEGREE), a * math.sin(b * _DEGREE)


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
        deferred = []
        # Complex eigenvalues (CLAMA): mode -> (real, imaginary part); rows are
        # ``mode, order, eigr, eigi, frequency, damping``.
        clama = {}
        for name, records in self.tables:
            if not name.startswith("CLAMA"):
                continue
            for mark, raw in records:
                if mark > -3 or len(raw) == 146 * ws or len(raw) % (6 * ws):
                    continue
                rows_i = _ints(s, raw).reshape(-1, 6)
                rows_f = _floats(s, raw).reshape(-1, 6)
                for ri, rf in zip(rows_i.tolist(), rows_f.tolist()):
                    clama.setdefault(int(ri[0]), (float(rf[2]), float(rf[3])))

        def skip(reason):
            if reason not in skipped:
                skipped.append(reason)

        for name, records in self.tables:
            nodal = name.startswith(_NODAL_PREFIXES)
            elemental = name.startswith(_ELEMENT_PREFIXES)
            grid_force = name.startswith("OGPF")
            if not (nodal or elemental or grid_force):
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
                format_code = int(h[8])
                num_wide = int(h[9])
                s_code = int(h[10])
                thermal = int(h[22])
                complex_ = bool(sort_code & 1)
                sort2 = bool(sort_code & 2)
                # Random tables: the hundreds of the table code name the
                # quantity (5 CRM, 6 PSD, 7 ATO, 8 RMS, 9 NO), the rest what.
                random = bool(sort_code & 4) or table_code >= 500
                what = table_code
                suffix = ""
                if random:
                    kind = table_code // 100
                    if not 5 <= kind <= 9 or complex_ or grid_force:
                        skip(f"{name} (random)")
                        continue
                    suffix = ("_CRM", "_PSD", "_ATO", "_RMS", "_NO")[kind - 5]
                    what = table_code % 100
                if complex_ and grid_force:
                    skip(f"{name} (complex)")
                    continue
                extra = {
                    "complex": complex_,
                    "mag_phase": complex_ and format_code == 3,
                    "suffix": suffix,
                }
                if grid_force:
                    if table_code != 19 or num_wide != 10:
                        skip(f"{name} (table code {table_code}, {num_wide} words)")
                        continue
                    info = ("gpf", "GRID_FORCE")
                elif nodal:
                    base = _NODAL_NAMES.get(what)
                    # MPC forces share the SPC forces' table code; the name tells.
                    if name.startswith("OQMG") and what in (3, 39):
                        base = "MPC_FORCE"
                    if base is None:
                        skip(f"{name} (table code {table_code})")
                        continue
                    if num_wide != (14 if complex_ else 8):
                        skip(f"{name} ({num_wide} words per node)")
                        continue
                    if thermal == 1:
                        if what != 1 or complex_:
                            skip(f"{name} (thermal table code {table_code})")
                            continue
                        base = "TEMPERATURE"
                    info = ("nodal", base, name.startswith("BOUG"), num_wide)
                else:
                    if what == 5:
                        family = "stress"
                        group = "STRAIN" if s_code & 8 else "STRESS"
                    elif what == 4:
                        family = "flux" if thermal == 1 else "force"
                        group = "HEAT_FLUX" if thermal == 1 else "ELEMENT_FORCE"
                    elif what in (18, 36):
                        family = "energy"
                        group = "ENERGY" if what == 18 else "KINETIC_ENERGY"
                    else:
                        skip(f"{name} (table code {table_code})")
                        continue
                    layout = _element_layout(
                        family,
                        etype,
                        num_wide,
                        s_code,
                        complex_,
                        random and family == "stress",
                    )
                    if layout is None:
                        # An energy table's word 3 is the total energy, not a type.
                        ename = (
                            "energy"
                            if family == "energy"
                            else _ELEMENT_TYPE_NAMES.get(etype, f"type {etype}")
                        )
                        cplx = " complex" if complex_ else ""
                        skip(f"{name} {ename}{cplx} ({num_wide} words)")
                        continue
                    info = ("element", group, layout, num_wide)
                w5 = int(h[4])
                moded = analysis in (2, 8, 9)
                if sort2:
                    if grid_force:
                        skip(f"{name} (SORT2)")
                        continue
                    # Every row makes (or joins) the step its first word names.
                    width = (info[3] if nodal else num_wide) * ws
                    if not width or len(raw) % width:
                        _fail(f"a SORT2 {name} record is not a whole number of rows")
                    firsts = _ints(s, raw)[:: width // ws].tolist()
                    times = _floats(s, raw)[:: width // ws].tolist()
                    for key5, t in zip(firsts, times):
                        key = (subcase, analysis, key5)
                        if key in step_index:
                            continue
                        step_index[key] = len(self.steps)
                        self.steps.append(
                            {
                                "subcase": subcase,
                                "analysis": analysis,
                                "mode": key5 if moded else 0,
                                "time": _time_of(analysis, key5, t, t),
                                "w5": key5,
                            }
                        )
                    extra["sort2"] = (subcase, analysis, header[4 * ws : 5 * ws], width)
                    self.blocks.append((None, info, raw, extra))
                    continue
                key = (subcase, analysis, w5)
                own = {
                    "subcase": subcase,
                    "analysis": analysis,
                    "mode": w5 if moded else 0,
                    "time": _time_of(analysis, w5, float(f[4]), float(f[5])),
                    "w5": w5,
                }
                if grid_force:
                    deferred.append((key, info, raw, name, extra, None))
                    continue
                # Energies too: ONRGY writes 0 in word 5 where the other tables
                # of a static step write the load set.
                if info[0] == "element" and info[1] in ("ENERGY", "KINETIC_ENERGY"):
                    deferred.append((key, info, raw, name, extra, own))
                    continue
                if key not in step_index:
                    step_index[key] = len(self.steps)
                    self.steps.append(own)
                self.blocks.append((step_index[key], info, raw, extra))
        # Grid point forces join a step the other tables made, never a new one:
        # MSC writes 0 in their word 5 where the other tables of a static step
        # write the load set, and a buckling run's forces carry analysis 2.
        for key, info, raw, name, extra, own in deferred:
            index = step_index.get(key)
            if index is None:
                same = [i for k, i in sorted(step_index.items()) if k[:2] == key[:2]]
                if not same and own is not None:
                    index = step_index[key] = len(self.steps)
                    self.steps.append(dict(own))
                elif len(same) != 1:
                    skip(f"{name} (no matching step)")
                    continue
                else:
                    index = same[0]
            self.blocks.append((index, info, raw, extra))
        for st in self.steps:
            st["eigi"] = 0.0
            if st["analysis"] == 9 and st["mode"] in clama:
                st["time"], st["eigi"] = clama[st["mode"]]
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


def _mesh_from_deck(path, want_nid):
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
    ids = np.asarray(mesh.points_id, dtype=np.int64)
    if want_nid and not np.array_equal(ids, np.arange(1, len(ids) + 1)):
        mesh.point_data["nastran:nid"] = ids
    del mesh.points_id
    del mesh.cells_id
    return mesh, grid_index, cell_index, offsets, sizes


def _build_mesh(nf, want_nid):
    s = nf.stream
    ids, xyz, cp, cd, spoints, cords = _read_grids(s, nf.tables)
    if not ids:
        deck, tried = _sibling_deck(nf.path)
        bgpdt = _bgpdt_grids(s, nf.tables) if deck is None else None
        if deck is None and bgpdt is None:
            _fail(
                "the file has no GEOM1 GRID records and no basic grid point table "
                "(rerun with PARAM,POST,-1 or provide the input deck beside it); "
                "looked for " + ", ".join(tried)
            )
        if deck is not None:
            return _mesh_from_deck(deck, want_nid) + (CoordSystems(), None)
        ids, xyz, cd, more = bgpdt
        spoints |= more
        cp = [0] * len(ids)
        warn(
            f"{WHO}: no GEOM1 table and no input deck beside the file; the "
            f"{len(ids)} points come from the basic grid point table (BGPDT)"
        )
        if any(c > 0 for c in cd):
            warn(
                f"{WHO}: without GEOM1 the output coordinate systems are unknown; "
                "results stay in them"
            )
        cd = [0] * len(ids)
    # A file can repeat a GRID in a second GEOM1 table (a restart's): kept
    # once when both definitions agree.
    grid_index = {}
    keep = []
    for i, g in enumerate(ids):
        j = grid_index.get(g)
        if j is not None:
            k = keep[j]
            same = cp[i] == cp[k] and cd[i] == cd[k]
            if not same or list(np.ravel(xyz[i])) != list(np.ravel(xyz[k])):
                _fail(f"GRID {g} is defined twice")
            continue
        grid_index[g] = len(keep)
        keep.append(i)
    if len(keep) != len(ids):
        ids = [ids[i] for i in keep]
        cp = [cp[i] for i in keep]
        cd = [cd[i] for i in keep]
        xyz = [xyz[i] for i in keep]
    points = np.asarray(xyz, dtype=np.float64).reshape(len(ids), 3)
    point_data = {}
    systems = apply_frames(points, point_data, cords, ids, cp, cd, WHO)
    if want_nid and list(map(int, ids)) != list(range(1, len(ids) + 1)):  # not 1 to n
        point_data["nastran:nid"] = np.asarray(ids, dtype=np.int64)
    cards = _read_elements(s, nf.tables, set(ids))
    extras = {}
    cells, cell_data, regions, cell_index, offsets, sizes = add_cells(
        cards, grid_index, spoints, _read_properties(s, nf.tables), WHO, extras
    )
    cell_data.update(extras)
    mesh = Mesh(points, cells, point_data=point_data, cell_data=cell_data)
    mesh.regions = regions
    return mesh, grid_index, cell_index, offsets, sizes, systems, cd


def _cord_fields(systems):
    """The resolved coordinate systems as ``nastran:CORD2_<id>:*`` field data: each one's
    type, its origin A and its axes X, Y and Z, all in basic coordinates."""
    out = {}
    for cid, (ctype, origin, axes) in sorted(systems.systems.items()):
        if cid == 0:
            continue
        prefix = f"nastran:CORD2_{cid}:"
        out[prefix + "TYPE"] = np.array([ctype], dtype=np.int64)
        out[prefix + "A"] = np.array(origin, dtype=np.float64)
        out[prefix + "X"] = np.array(axes[0], dtype=np.float64)
        out[prefix + "Y"] = np.array(axes[1], dtype=np.float64)
        out[prefix + "Z"] = np.array(axes[2], dtype=np.float64)
    return out


def read(filename, points_only=False, arrays=None, time_step=0):
    nf = _File(filename)
    # GRID ids that are not 1 to n, as point data: only when point data is wanted.
    want_nid = not points_only and (arrays is None or "nastran:nid" in arrays)
    mesh, grid_index, cell_index, offsets, sizes, systems, cd = _build_mesh(
        nf, want_nid
    )
    mesh.field_data.update(_read_property_values(nf.stream, nf.tables))
    mesh.field_data.update(_read_pcomp_values(nf.stream, nf.tables))
    mesh.field_data.update(_read_pbeam_values(nf.stream, nf.tables))
    mesh.field_data.update(_read_bush_values(nf.stream, nf.tables))
    mesh.field_data.update(_read_line_values(nf.stream, nf.tables))
    mesh.field_data.update(_cord_fields(systems))
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
    if step["analysis"] == 9:
        mesh.field_data["nastran:eigi"] = np.array([step["eigi"]], dtype=np.float64)
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
    wide = {}  # name -> [(cell, column, value)], laid out once every block is read

    def push(name, cell, col, value):
        wide.setdefault(name, []).append((cell, col, value))

    def raw_rows(raw, width):
        n = width * s.ws
        return [raw[i : i + n] for i in range(0, len(raw), n)]

    for step_id, info, raw, extra in nf.blocks:
        if "sort2" in extra:
            # This step's rows of the entity, rewritten as SORT1 rows.
            subcase, analysis, entity, width = extra["sort2"]
            if (subcase, analysis) != (step["subcase"], step["analysis"]):
                continue
            firsts = _ints(s, raw)[:: width // s.ws].tolist()
            raw = b"".join(
                entity + raw[r * width + s.ws : (r + 1) * width]
                for r, first in enumerate(firsts)
                if first == step["w5"]
            )
            if not raw:
                continue
        elif step_id != index:
            continue
        ints = _ints(s, raw)
        floats = _floats(s, raw)
        if info[0] == "nodal":
            base, nw = info[1], info[3]
            if len(ints) % nw:
                _fail(f"a {base} record is not a whole number of rows")
            rows_i = ints.reshape(-1, nw)
            rows_f = floats.reshape(-1, nw)
            nid = rows_i[:, 0] // 10
            targets = [(r, grid_index.get(int(g))) for r, g in enumerate(nid.tolist())]
            sfx = extra["suffix"]
            # (name, first word of the value or real part, of the imaginary part,
            # components, part: 0 the value, 1 the real, 2 the imaginary)
            if base == "TEMPERATURE":
                outputs = [(base + sfx, 2, 0, 1, 0)]
            elif extra["complex"]:
                outputs = [
                    (base + "_real", 2, 8, 3, 1),
                    (base + "_imag", 2, 8, 3, 2),
                    (base + "_ROT_real", 5, 11, 3, 1),
                    (base + "_ROT_imag", 5, 11, 3, 2),
                ]
            else:
                outputs = [(base + sfx, 2, 0, 3, 0), (base + "_ROT" + sfx, 5, 0, 3, 0)]
            for name, c0, c1, nc, part in outputs:
                if not wants(name):
                    continue
                values = point_arrays.get(name)
                if values is None:
                    values = point_arrays[name] = np.full((npts, nc), _NAN)
                for r, p in targets:
                    if p is not None and np.isnan(values[p, 0]):
                        v = rows_f[r, c0 : c0 + nc].tolist()
                        if part:
                            im = rows_f[r, c1 : c1 + nc].tolist()
                            pairs = [
                                _complex(a, b, extra["mag_phase"])
                                for a, b in zip(v, im)
                            ]
                            v = [z[part - 1] for z in pairs]
                        # Results are in the GRID's output system (CD) unless the
                        # table is BOUG*; random ones (spectral densities, RMS
                        # ...) are not vectors and stay there.
                        if nc == 3 and not info[2] and cd is not None and not sfx:
                            rotate_to_basic(systems, cd[p], mesh.points[p].tolist(), v)
                        values[p] = v
        elif info[0] == "gpf":
            # grid*10+device, element (0 for totals, loads, SPC forces...), its
            # 8-character name, F1 F2 F3 M1 M2 M3 in the GRID's output system.
            base = info[1]
            if len(ints) % 10:
                _fail(f"a {base} record is not a whole number of entries")
            rows_i = ints.reshape(-1, 10)
            rows_f = floats.reshape(-1, 10)
            names = raw_rows(raw, 10)
            for r in range(len(rows_i)):
                p = grid_index.get(int(rows_i[r, 0]) // 10)
                if p is None:
                    continue
                v = rows_f[r, 4:10].tolist()
                if cd is not None:
                    rotate_to_basic(systems, cd[p], mesh.points[p].tolist(), v)
                eid = int(rows_i[r, 1])
                if eid > 0:
                    c = cell_index.get(eid)
                    if c is None:
                        continue
                    pos = _node_position(mesh, offsets, c, p)
                    if pos is None:
                        continue
                    for k, m in enumerate(_GPF_MEMBERS):
                        name = f"{base}:{m}"
                        if wants(name):
                            push(name, c, pos, v[k])
                    continue
                label = names[r][2 * s.ws : 4 * s.ws].decode("latin-1")
                label = "".join(ch for ch in label if ch not in " *\x00")
                for k, m in enumerate(_GPF_MEMBERS):
                    name = f"{base}:{label}:{m}"
                    if not wants(name):
                        continue
                    values = point_arrays.get(name)
                    if values is None:
                        values = point_arrays[name] = np.full((npts, 1), _NAN)
                    if np.isnan(values[p, 0]):
                        values[p, 0] = v[k]
        else:
            _, group, layout, num_wide = info
            if len(ints) % num_wide:
                _fail(f"a {group} record is not a whole number of elements")
            rows_i = ints.reshape(-1, num_wide)
            rows_f = floats.reshape(-1, num_wide)
            kind = layout["kind"]
            members = layout["members"]
            sfx = extra["suffix"]
            mag_phase = extra["mag_phase"]

            def member_values(row_f, at):
                """(member, value) of each member at word ``at``: itself (with a
                random table's suffix), or its real and imaginary parts."""
                out = []
                for word, imag, member in members:
                    a = float(row_f[at + word])
                    if imag is None:
                        out.append((member + sfx, a))
                        continue
                    re, im = _complex(a, float(row_f[at + imag]), mag_phase)
                    out.append((member + "_real", re))
                    out.append((member + "_imag", im))
                return out

            stations = {}
            for r in range(len(rows_i)):
                eid = int(rows_i[r, 0]) // 10
                station = 0
                if kind == "station":
                    station = stations.get(eid, 0)
                    stations[eid] = station + 1
                c = cell_index.get(eid)
                if c is None:
                    continue
                row_f = rows_f[r]
                if kind in ("ply", "station"):
                    ply = int(rows_i[r, 1]) if kind == "ply" else 1
                    if ply < 1:
                        continue
                    col = ply - 1 if kind == "ply" else station
                    suffix = "@ply" if kind == "ply" else "@station"
                    for member, value in member_values(row_f, 0):
                        name = f"{group}:{member}{suffix}"
                        if wants(name):
                            push(name, c, col, value)
                    continue
                # row: the values; blocks: the first block's (centre, or end A)
                if kind != "blocks" or layout["centre"]:
                    first = layout["first"] if kind == "blocks" else 0
                    for member, value in member_values(row_f, first):
                        name = f"{group}:{member}"
                        if not wants(name):
                            continue
                        values = cell_arrays.get(name)
                        if values is None:
                            values = cell_arrays[name] = np.full(ncells, _NAN)
                        if np.isnan(values[c]):
                            values[c] = value
                if kind != "blocks":
                    continue
                beam = layout["stations"]
                for k in range(
                    0 if beam or not layout["centre"] else 1, layout["blocks"]
                ):
                    at = layout["first"] + k * layout["block"]
                    col = k
                    # a beam station with no GRID and no distance was not output
                    if (
                        beam
                        and k > 0
                        and int(rows_i[r, at]) == 0
                        and row_f[at + 1] == 0.0
                    ):
                        continue
                    if not beam:
                        g = grid_index.get(int(rows_i[r, at]))
                        if g is None:
                            continue
                        col = _node_position(mesh, offsets, c, g)
                        if col is None:
                            continue
                    suffix = "@station" if beam else "@corner"
                    for member, value in member_values(row_f, at):
                        name = f"{group}:{member}{suffix}"
                        if wants(name):
                            push(name, c, col, value)
    for name, values in point_arrays.items():
        mesh.point_data[name] = values[:, 0] if values.shape[1] == 1 else values
    for name, values in cell_arrays.items():
        mesh.cell_data[name] = [
            values[offsets[b] : offsets[b] + sizes[b]] for b in range(len(sizes))
        ]
    for name, entries in wide.items():
        width = max(col for _, col, _ in entries) + 1
        blocks = [np.full((sizes[b], width), _NAN) for b in range(len(sizes))]
        for cell, col, value in entries:  # the first value of a cell and column wins
            b = bisect.bisect_right(offsets, cell) - 1
            if np.isnan(blocks[b][cell - offsets[b], col]):
                blocks[b][cell - offsets[b], col] = value
        mesh.cell_data[name] = blocks
        mesh.field_data["nastran:layout:" + name] = np.array([width, 1], dtype=np.int64)
    return mesh


_GPF_MEMBERS = ("F1", "F2", "F3", "M1", "M2", "M3")


def _node_position(mesh, offsets, cell, point):
    """The position of ``point`` in the connectivity of global cell ``cell``."""
    b = bisect.bisect_right(offsets, cell) - 1
    row = mesh.cells[b].data[cell - offsets[b]].tolist()
    return row.index(point) if point in row else None

"""OpenRadioss / Radioss time-history files (``<run>T01``...).

A big-endian Fortran unformatted file: a title, the run's date, the part,
material, property, subset and TH-group descriptions, then one output per time
(the time, the global variables, the part and subset variables, one record per
TH group). The layout follows OpenRadioss's ``th_to_csv`` converter (MIT); no
code is copied.

The writer (:func:`write`, :class:`SeriesWriter`) is the reader's inverse and
the twin of the C++ ``write_radioss_th``: it writes the same bytes.
"""

import re
import struct

import numpy as np

from .. import _provenance
from .._common import warn
from .._exceptions import ReadError, WriteError
from .._files import is_buffer
from .._fortran_records import fortran_record, fortran_records
from .._mesh import CellBlock, Mesh

TIME_KEY = "meshio:time"

# Big-endian Fortran records with 4-byte markers.
_LAYOUT = (4, ">")

# Global variables by code (th_to_csv's column titles).
_GLOBAL_NAMES = (
    "internal_energy",
    "kinetic_energy",
    "x_momentum",
    "y_momentum",
    "z_momentum",
    "mass",
    "time_step",
    "rotation_energy",
    "external_work",
    "spring_energy",
    "contact_energy",
    "hourglass_energy",
    "elastic_contact_energy",
    "frictional_contact_energy",
    "damping_contact_energy",
    "plastic_work",
    "added_mass",
    "percentage_added_mass",
    "inlet_mass",
    "outlet_mass",
    "inlet_energy",
    "outlet_energy",
)

# Part and subset variables by code (Radioss's keywords).
_PART_NAMES = {
    1: "IE",
    2: "KE",
    3: "XMOM",
    4: "YMOM",
    5: "ZMOM",
    6: "MASS",
    7: "HE",
    8: "TURBKE",
    9: "XCG",
    10: "YCG",
    11: "ZCG",
    12: "XXMOM",
    13: "YYMOM",
    14: "ZZMOM",
    15: "IXX",
    16: "IYY",
    17: "IZZ",
    18: "IXY",
    19: "IYZ",
    20: "IZX",
    21: "RIE",
    22: "KERB",
    23: "RKERB",
    24: "RKE",
    25: "ERODED",
    28: "HEAT",
    29: "VX",
    30: "VY",
    31: "VZ",
    32: "PW",
}

# A TH group's kind by its type (the /TH/<kind> keyword it came from).
_GROUP_KINDS = {
    0: "node",
    1: "brick",
    2: "quad",
    3: "shell",
    4: "truss",
    5: "beam",
    6: "spring",
    7: "sh3n",
    51: "sphcel",
    101: "inter",
    102: "rwall",
    103: "rbody",
    104: "sectio",
    107: "monvol",
    108: "accel",
}


def _fail(message):
    raise ReadError(f"Radioss time history: {message}")


def _global_name(code):
    if 1 <= code <= len(_GLOBAL_NAMES):
        return _GLOBAL_NAMES[code - 1]
    return f"var{code}"


def _part_name(code):
    return _PART_NAMES.get(code, f"var{code}")


def _group_kind(kind):
    return _GROUP_KINDS.get(kind, f"type{kind}")


def is_radioss_th(head):
    """An 84-byte big-endian title record holding a plausible file version,
    then the 80-byte date record."""
    if len(head) < 8 or struct.unpack(">I", head[:4])[0] != 84:
        return False
    version = struct.unpack(">i", head[4:8])[0]
    if not 1000 <= version <= 99999:
        return False
    return len(head) < 96 or struct.unpack(">II", head[88:96]) == (84, 80)


class _File:
    def __init__(self, filename):
        with open(filename, "rb") as f:
            data = f.read()
        if not is_radioss_th(data[:96]):
            _fail(f"'{filename}' is not a time-history file")
        self.data = data
        # Split the records; a truncated last one ends the file.
        self.records = fortran_records(
            data, _LAYOUT, "Radioss time history", lenient=True
        )
        self.next = 0
        self._parse()

    def _take(self, size, what):
        if self.next >= len(self.records):
            _fail(f"the file ends before {what}")
        offset, length = self.records[self.next]
        if size is not None and length != size:
            _fail(f"unexpected record size for {what}")
        self.next += 1
        return offset

    def _int(self, offset):
        return struct.unpack_from(">i", self.data, offset)[0]

    def _ints(self, n, what):
        if n == 0:
            return []
        offset = self._take(4 * n, what)
        return list(struct.unpack_from(f">{n}i", self.data, offset))

    @staticmethod
    def _count(n, what):
        if n < 0 or n > 1 << 28:
            _fail(f"implausible {what} count")
        return n

    def _parse(self):
        version = self._int(self._take(84, "the title"))
        title = 100 if version >= 4021 else 80 if version >= 3041 else 40
        self._take(80, "the date")
        self.factors = None
        if version > 3050:
            self._take(4, "the additional records")
            self._take(4, "the title length")
            offset = self._take(12, "the unit factors")
            self.factors = struct.unpack_from(">3f", self.data, offset)
        counts = self._ints(6, "the hierarchy counts")
        nparts, nmats, ngeos, nsubs, ngroups, nglob = (
            self._count(c, w)
            for c, w in zip(
                counts,
                (
                    "part",
                    "material",
                    "property",
                    "subset",
                    "TH group",
                    "global variable",
                ),
            )
        )
        self.globals = self._ints(nglob, "the global variables")
        self.parts = []
        for _ in range(nparts):
            p = self._take(4 + title + 16, "a part")
            n = self._count(self._int(p + 4 + title + 12), "part variable")
            self.parts.append((self._int(p), self._ints(n, "a part")))
        for _ in range(nmats + ngeos):
            self._take(4 + title, "a material or property")
        self.subsets = []
        for _ in range(nsubs):
            p = self._take(20 + title, "a subset")
            self._ints(self._count(self._int(p + 8), "child subset"), "a subset")
            self._ints(self._count(self._int(p + 12), "subset part"), "a subset")
            n = self._count(self._int(p + 16), "subset variable")
            self.subsets.append((self._int(p), self._ints(n, "a subset")))
        self.groups = []
        for _ in range(ngroups):
            p = self._take(20 + title, "a TH group")
            n = self._count(self._int(p + 12), "TH group entity")
            ids = [
                self._int(self._take(4 + title, "a TH group entity")) for _ in range(n)
            ]
            m = self._count(self._int(p + 16), "TH group variable")
            codes = self._ints(m, "a TH group")
            self.groups.append((self._int(p), self._int(p + 4), ids, codes))
        # The outputs: time, globals, parts, subsets, one record per group. A
        # run stopped mid-write leaves an incomplete last output, not read.
        nparts_vars = sum(len(c) for _, c in self.parts)
        nsubs_vars = sum(len(c) for _, c in self.subsets)
        sizes = [4 * n for n in (nglob, nparts_vars, nsubs_vars) if n > 0]
        sizes += [4 * len(ids) * len(codes) for _, _, ids, codes in self.groups]
        self.outputs = []
        self.times = []
        while self.next + 1 + len(sizes) <= len(self.records):
            if self.records[self.next][1] != 4:
                _fail("unexpected record size for an output's time")
            for k, size in enumerate(sizes):
                if self.records[self.next + 1 + k][1] != size:
                    _fail(
                        "unexpected record size in the output at record "
                        f"{self.next + 1 + k}"
                    )
            self.outputs.append(self.next)
            offset = self.records[self.next][0]
            self.times.append(float(struct.unpack_from(">f", self.data, offset)[0]))
            self.next += 1 + len(sizes)
        if not self.outputs:
            _fail("the file has no complete output")

    def floats(self, index):
        offset, length = self.records[index]
        return np.frombuffer(self.data, ">f4", length // 4, offset).astype(np.float64)


def time_values(filename):
    return list(_File(filename).times)


def read(filename, points_only=False, arrays=None, time_step=0):
    f = _File(filename)
    n = len(f.outputs)
    index = time_step + n if time_step < 0 else time_step
    if not 0 <= index < n:
        raise ReadError(
            f"time step {time_step} is out of range: the file has {n} step(s)"
        )

    def want(name):
        return arrays is None or name in arrays

    mesh = Mesh(np.zeros((0, 3)), [CellBlock("vertex", np.zeros((0, 1), np.int64))])
    mesh.field_data[TIME_KEY] = np.array([f.times[index]], dtype=np.float64)
    mesh.time_values = list(f.times)
    if points_only:
        return mesh

    r = f.outputs[index] + 1
    if f.globals:
        values = f.floats(r)
        r += 1
        for code, value in zip(f.globals, values):
            name = "radioss_th:global:" + _global_name(code)
            if want(name):
                mesh.field_data[name] = np.array([value])
    for kind, entities in (("part", f.parts), ("subset", f.subsets)):
        if not any(codes for _, codes in entities):
            continue
        values = f.floats(r)
        r += 1
        k = 0
        for eid, codes in entities:
            for code in codes:
                name = f"radioss_th:{kind}:{eid}:{_part_name(code)}"
                if want(name):
                    mesh.field_data[name] = np.array([values[k]])
                k += 1
    for gid, kind, ids, codes in f.groups:
        name = f"radioss_th:{_group_kind(kind)}:{gid}"
        values = f.floats(r)
        r += 1
        if name in mesh.field_data:
            continue
        if want(name):
            mesh.field_data[name] = values.reshape(len(ids), len(codes))
        if want(name + ":ids"):
            mesh.field_data[name + ":ids"] = np.array(ids, dtype=np.int64)
        if want(name + ":variables"):
            mesh.field_data[name + ":variables"] = np.array(codes, dtype=np.int64)
    if f.factors is not None and want("radioss_th:unit_factors"):
        mesh.field_data["radioss_th:unit_factors"] = np.array(f.factors, np.float64)
    return mesh


# --------------------------------------------------------------------------- #
# The writer: the reader's inverse                                            #
# --------------------------------------------------------------------------- #

# Version 3040 (40-character titles, no unit factors) is what OpenRadioss itself
# writes and ``th_to_csv`` reads; with unit factors the file is 4021 (100-character
# titles): the reader's two branches that a writer can reach.
_VERSION_PLAIN = 3040
_VERSION_FACTORS = 4021

_INT = re.compile(r"-?[0-9]+")
_INT32 = (-(2**31), 2**31 - 1)


def _wfail(message):
    raise WriteError(f"Radioss time history: {message}")


def _parse_int(text):
    """A whole-string decimal integer that fits an int32, else ``None``."""
    if not _INT.fullmatch(text):
        return None
    value = int(text)
    return value if _INT32[0] <= value <= _INT32[1] else None


def _code(name, name_of):
    """The code a variable name stands for: found in the reader's own table (so
    the two directions cannot drift), else ``var<code>``."""
    for code in range(1, 65):
        if name_of(code) == name:
            return code
    if name.startswith("var"):
        code = _parse_int(name[3:])
        if code is not None and code >= 1:
            return code
    return None


def _kind_type(name):
    """The TH group type a kind name stands for (``type<N>`` for one the reader
    has no name for)."""
    for kind in range(201):
        if _group_kind(kind) == name:
            return kind
    if name.startswith("type"):
        kind = _parse_int(name[4:])
        if kind is not None and kind >= 0:
            return kind
    return None


def _scalar(array):
    values = np.asarray(array, dtype=np.float64).ravel()
    if values.size != 1:
        return None
    return values


def _float32(values):
    """Values as the big-endian ``float32`` the file stores (overflow is inf)."""
    with np.errstate(over="ignore", invalid="ignore"):
        return np.asarray(values, dtype=np.float64).astype(">f4")


def _int_array(array, expected, name):
    values = np.asarray(array, dtype=np.float64).ravel()
    if values.size != expected:
        _wfail(f"'{name}' holds {values.size} values, expected {expected}")
    out = []
    for value in values:
        if not (_INT32[0] <= value <= _INT32[1]) or value != np.floor(value):
            _wfail(f"'{name}' holds a value that is not an integer in 32-bit range")
        out.append(int(value))
    return out


class _Step:
    """One output's worth of a mesh: everything the file stores."""

    def __init__(self, mesh, time, dropped):
        self.time = np.float32(0.0)
        self.factors = None
        self.globals = []  # (code, name, value)
        self.parts = []  # (id, [(code, name, value)])
        self.subsets = []
        self.groups = []  # (type, id, kind, ids, codes, values)
        parts, subsets, groups = {}, {}, {}
        for name, array in mesh.field_data.items():
            if name == TIME_KEY:
                values = np.asarray(array, dtype=np.float64).ravel()
                if time is None and values.size:
                    self.time = _float32(values[:1])[0]
                continue
            p = name.split(":")
            taken = False
            if len(p) >= 2 and p[0] == "radioss_th":
                if len(p) == 2 and p[1] == "unit_factors":
                    values = np.asarray(array, dtype=np.float64).ravel()
                    if values.size != 3:
                        _wfail("'radioss_th:unit_factors' must hold three values")
                    self.factors = _float32(values)
                    taken = True
                elif len(p) == 3 and p[1] == "global":
                    code = _code(p[2], _global_name)
                    value = _scalar(array)
                    if code is not None and value is not None:
                        self.globals.append((code, p[2], value[0]))
                        taken = True
                elif len(p) == 4 and p[1] in ("part", "subset"):
                    eid = _parse_int(p[2])
                    code = _code(p[3], _part_name)
                    value = _scalar(array)
                    if eid is not None and code is not None and value is not None:
                        target = parts if p[1] == "part" else subsets
                        target.setdefault(eid, []).append((code, p[3], value[0]))
                        taken = True
                elif len(p) in (3, 4):
                    kind = _kind_type(p[1])
                    gid = _parse_int(p[2])
                    if (
                        kind is not None
                        and gid is not None
                        and (len(p) == 3 or p[3] in ("ids", "variables"))
                    ):
                        g = groups.setdefault((kind, gid, p[1]), {"names": []})
                        g["data" if len(p) == 3 else p[3]] = array
                        g["names"].append(name)
                        taken = True
            if not taken:
                dropped.append(name)
        if time is not None:
            self.time = _float32([time])[0]

        self.globals.sort(key=lambda v: (v[0], v[1]))
        for source, target in ((parts, self.parts), (subsets, self.subsets)):
            for eid in sorted(source):
                target.append((eid, sorted(source[eid], key=lambda v: (v[0], v[1]))))
        for key in sorted(groups):
            g = groups[key]
            kind, gid, kname = key
            base = f"radioss_th:{kname}:{gid}"
            if "data" not in g:
                dropped.extend(g["names"])
                continue
            data = np.asarray(g["data"])
            if data.ndim != 2:
                _wfail(f"'{base}' must be a 2-D array (entities x variables)")
            n, m = data.shape
            ids = (
                _int_array(g["ids"], n, base + ":ids")
                if "ids" in g
                else list(range(1, n + 1))
            )
            codes = (
                _int_array(g["variables"], m, base + ":variables")
                if "variables" in g
                else list(range(1, m + 1))
            )
            self.groups.append((kind, gid, kname, ids, codes, _float32(data.ravel())))
        dropped.sort()

    def same_layout(self, other):
        """Whether two outputs have the same variables (values and time may
        differ)."""

        def entities(items):
            return [(eid, [v[0] for v in vars_]) for eid, vars_ in items]

        return (
            (self.factors is None) == (other.factors is None)
            and [v[0] for v in self.globals] == [v[0] for v in other.globals]
            and entities(self.parts) == entities(other.parts)
            and entities(self.subsets) == entities(other.subsets)
            and [g[:5] for g in self.groups] == [g[:5] for g in other.groups]
        )


def _ints(values):
    return fortran_record(struct.pack(f">{len(values)}i", *values))


def _text(text, width):
    return text[:width].ljust(width).encode("ascii", "replace")


def _header(step):
    """The descriptions: titles, the unit factors, the hierarchy and the
    variable codes of every part, subset and TH group."""
    factors = step.factors is not None
    title = 100 if factors else 40
    out = [
        fortran_record(
            struct.pack(">i", _VERSION_FACTORS if factors else _VERSION_PLAIN)
            + _text(_provenance.lines(_provenance.SlotTier.BOUNDED)[0], 80)
        ),
        fortran_record(_text("meshio++", 80)),
    ]
    if factors:
        out += [
            _ints([1]),
            _ints([title]),
            fortran_record(step.factors.astype(">f4").tobytes()),
        ]
    out.append(
        _ints(
            [
                len(step.parts),
                0,
                0,
                len(step.subsets),
                len(step.groups),
                len(step.globals),
            ]
        )
    )
    if step.globals:
        out.append(_ints([v[0] for v in step.globals]))
    for eid, vars_ in step.parts:
        out.append(
            fortran_record(
                struct.pack(">i", eid)
                + _text(f"PART {eid}", title)
                + struct.pack(">4i", 0, 0, 0, len(vars_))
            )
        )
        if vars_:
            out.append(_ints([v[0] for v in vars_]))
    for eid, vars_ in step.subsets:
        out.append(
            fortran_record(
                struct.pack(">5i", eid, 0, 0, 0, len(vars_))
                + _text(f"SUBSET {eid}", title)
            )
        )
        if vars_:
            out.append(_ints([v[0] for v in vars_]))
    for kind, gid, kname, ids, codes, _ in step.groups:
        out.append(
            fortran_record(
                struct.pack(">5i", gid, kind, 0, len(ids), len(codes))
                + _text(f"TH {kname} {gid}", title)
            )
        )
        for eid in ids:
            out.append(
                fortran_record(struct.pack(">i", eid) + _text(f"{kname} {eid}", title))
            )
        if codes:
            out.append(_ints(codes))
    return b"".join(out)


def _output(step):
    """One output: the time, the globals, the part and subset values, one
    record per TH group."""

    def floats(values):
        return fortran_record(_float32(values).tobytes())

    out = [floats([step.time])]
    if step.globals:
        out.append(floats([v[2] for v in step.globals]))
    for entities in (step.parts, step.subsets):
        values = [v[2] for _, vars_ in entities for v in vars_]
        if values:
            out.append(floats(values))
    for group in step.groups:
        out.append(fortran_record(group[5].tobytes()))
    return b"".join(out)


def _warn_dropped(dropped):
    if dropped:
        warn(
            "Radioss time history: field data without a place in the file is "
            f"not written: {', '.join(dropped)}"
        )


@_provenance.slotless_writer
def write(filename, mesh):
    dropped = []
    step = _Step(mesh, None, dropped)
    _warn_dropped(dropped)
    data = _header(step) + _output(step)
    if is_buffer(filename, "w"):
        filename.write(data)
    else:
        with open(filename, "wb") as fh:
            fh.write(data)


class SeriesWriter:
    """Write a time-history file output by output: the descriptions with the
    first step, then one output per :meth:`write`. Every step must have the
    first step's variables; a step that does not is a ``WriteError``. Use as a
    context manager, or call :meth:`close`."""

    def __init__(self, filename):
        # A path is opened as `write` opens it; a caller's buffer is written
        # into and left open.
        self._owned = not is_buffer(filename, "w")
        self._fh = open(filename, "wb") if self._owned else filename
        self._first = None
        self._step = 0

    def write(self, time, mesh):
        dropped = []
        step = _Step(mesh, float(time), dropped)
        data = b""
        if self._first is None:
            _warn_dropped(dropped)
            data = _header(step)
            self._first = step
        elif not self._first.same_layout(step):
            _wfail(
                f"series: step {self._step}'s variables differ from the first "
                "step's; a time-history file holds one set -- write one file per "
                "step with '{step}'"
            )
        self._fh.write(data + _output(step))
        self._step += 1

    def close(self):
        if self._fh is not None:
            if self._owned:
                self._fh.close()
            self._fh = None

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False

"""OpenRadioss / Radioss time-history files (``<run>T01``...).

A big-endian Fortran unformatted file: a title, the run's date, the part,
material, property, subset and TH-group descriptions, then one output per time
(the time, the global variables, the part and subset variables, one record per
TH group). The layout follows OpenRadioss's ``th_to_csv`` converter (MIT); no
code is copied.
"""

import struct

import numpy as np

from .._exceptions import ReadError
from .._mesh import CellBlock, Mesh

TIME_KEY = "meshio:time"

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
        self.records = []
        p = 0
        while p + 4 <= len(data):
            (n,) = struct.unpack_from(">I", data, p)
            if (
                p + 8 + n > len(data)
                or struct.unpack_from(">I", data, p + 4 + n)[0] != n
            ):
                break
            self.records.append((p + 4, n))
            p += 8 + n
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

"""Ansys MAPDL binary results (``.rst``, ``.rth``): the pure-Python reference reader.

The twin of ``src/cpp/src/formats/ansys_rst.cpp``. A results file is a sequence
of Fortran-style records -- ``i32 length`` (in 4-byte words), ``i32 flags``, the
data, a trailer word -- addressed by word offsets from the start of the file.
The layout follows Ansys's ``fdresu.inc`` as the open reader pymapdl-reader
(MIT) documents it; see ``doc/formats/ansys_rst.md``.

The mesh (nodes, element types, elements, components) comes from the geometry
records and becomes cells exactly as a ``.cdb`` deck's elements do. Each result
set is one step: its nodal DOF solution becomes point data, rotated from each
node's coordinate system to the global one.
"""

from __future__ import annotations

import math
import pathlib

import numpy as np

from .._common import warn
from .._exceptions import ReadError
from .._mesh import CellBlock, Mesh
from .._regions import Region
from ..ansysInp._ansysInp import _build

_LABEL = "Ansys .rst"
_UNDEFINED = 2.0**100  # MAPDL's "no value"

# Record flags (the flags word's top byte).
_BSPARSE = 0x08
_WSPARSE = 0x10
_ZLIB = 0x20
_PREC = 0x40
_INT = 0x80

# DOF codes (fdresu.inc) of the nodal solution; vector DOFs come in triples.
_VECTORS = {"U": (1, 2, 3), "ROT": (4, 5, 6), "A": (7, 8, 9), "V": (10, 11, 12)}
_SCALARS = {
    16: "WARP",
    17: "CONC",
    18: "HDSP",
    19: "PRES",
    20: "TEMP",
    21: "VOLT",
    22: "MAG",
    23: "ENKE",
    24: "ENDS",
    25: "EMF",
    26: "CURR",
}


def _fail(what):
    raise ReadError(f"{_LABEL}: {what}")


class _File:
    """The file's bytes, read record by record."""

    def __init__(self, filename):
        with open(filename, "rb") as f:
            self.data = f.read()
        self.words = np.frombuffer(self.data, "<i4", len(self.data) // 4)

    def record(self, ptr):
        """``(values, next pointer)`` of the record at word ``ptr``."""
        w = self.words
        if ptr < 0 or ptr + 2 > len(w):
            _fail(f"record pointer {ptr} is outside the file")
        n = int(w[ptr])
        flags = (int(w[ptr + 1]) >> 24) & 0xFF
        end = ptr + 2 + n
        if n < 0 or end > len(w):
            _fail(f"the record at word {ptr} runs past the end of the file")
        if flags & _ZLIB:
            _fail(
                "zlib-compressed records are not supported; write the file with "
                "/FCOMP,RST,0"
            )
        if flags & _INT:
            dtype = np.dtype("<i2" if flags & _PREC else "<i4")
        else:
            dtype = np.dtype("<f4" if flags & _PREC else "<f8")
        raw = self.data[(ptr + 2) * 4 : end * 4]
        if flags & _BSPARSE:
            values = _bsparse(raw, dtype, ptr)
        elif flags & _WSPARSE:
            values = _wsparse(raw, dtype, ptr)
        else:
            values = np.frombuffer(raw, dtype, len(raw) // dtype.itemsize)
        return values, ptr + n + 3

    def ints(self, ptr):
        values, _ = self.record(ptr)
        return values.astype(np.int64)


def _bsparse(raw, dtype, ptr):
    """A bit-sparse record: a size, a bit mask, the values of the set bits."""
    if len(raw) < 8:
        _fail(f"the bit-sparse record at word {ptr} is truncated")
    size, bits = np.frombuffer(raw, "<i4", 2)
    size, bits = int(size), int(bits) & 0xFFFFFFFF
    if not 0 <= size <= 32:
        _fail(f"the bit-sparse record at word {ptr} has size {size}")
    packed = np.frombuffer(raw, dtype, (len(raw) - 8) // dtype.itemsize, 8)
    out = np.zeros(size, dtype)
    k = 0
    for i in range(size):
        if bits >> i & 1:
            if k >= len(packed):
                _fail(f"the bit-sparse record at word {ptr} is truncated")
            out[i] = packed[k]
            k += 1
    return out


def _wsparse(raw, dtype, ptr):
    """A windowed-sparse record: a size, a window count, then windows -- an
    isolated value (``loc > 0``), a run of values or a repeated constant."""
    if dtype.itemsize == 2:
        _fail(f"the record at word {ptr} is a windowed-sparse int16 record")
    words = np.frombuffer(raw, "<i4", len(raw) // 4)
    if len(words) < 2:
        _fail(f"the windowed-sparse record at word {ptr} is truncated")
    size, n_windows = int(words[0]), int(words[1])
    shift = dtype.itemsize // 4
    out = np.zeros(max(size, 0), dtype)
    pos = 2

    def take(count):
        nonlocal pos
        if pos + shift * count > len(words):
            _fail(f"the windowed-sparse record at word {ptr} is truncated")
        values = np.frombuffer(raw, dtype, count, pos * 4)
        pos += shift * count
        return values

    for _ in range(max(n_windows, 0)):
        if pos >= len(words):
            _fail(f"the windowed-sparse record at word {ptr} is truncated")
        loc = int(words[pos])
        pos += 1
        if loc > 0:
            if loc >= size:
                _fail(f"the windowed-sparse record at word {ptr} is corrupt")
            out[loc] = take(1)[0]
            continue
        if pos >= len(words):
            _fail(f"the windowed-sparse record at word {ptr} is truncated")
        start, length = -loc, int(words[pos])
        pos += 1
        count = abs(length)
        if start + count > size:
            _fail(f"the windowed-sparse record at word {ptr} is corrupt")
        if length > 0:
            out[start : start + count] = take(count)
        else:
            out[start : start + count] = take(1)[0]
    return out


def _close(a, b):
    """numpy.isclose's default tolerances: how mode pairs are recognised."""
    return abs(a - b) <= 1e-8 + 1e-5 * abs(b)


def _pointer(header, lo, hi, fallback=None):
    """A 64-bit word pointer split over ``header[lo]`` (low, unsigned) and
    ``header[hi]``; an older file's shorter header has only the 32-bit
    ``header[fallback]`` (or ``header[lo]`` itself)."""
    if hi < len(header):
        return (int(header[lo]) & 0xFFFFFFFF) | (int(header[hi]) << 32)
    if fallback is not None:
        lo = fallback
    return int(header[lo]) & 0xFFFFFFFF if lo < len(header) else 0


def _at(header, k):
    return int(header[k]) if k < len(header) else 0


_BLANKS = "".join(chr(c) for c in range(33))  # spaces and control characters


def _name(words):
    """A name stored as 4-character words, each byte-reversed."""
    raw = b"".join(int(v & 0xFFFFFFFF).to_bytes(4, "big") for v in words)
    return raw.split(b"\x00")[0].decode("latin-1").strip(_BLANKS)


def _ranges(values):
    """Ids with ``-last`` closing a run opened by the value before it."""
    ids = []
    for v in values:
        v = int(v)
        if v > 0:
            ids.append(v)
        elif v < 0 and ids:
            ids.extend(range(ids[-1] + 1, -v + 1))
    return ids


# The element solution: each element's pointer table (``ptrESL``) has one entry
# per kind of record, in this order (fdresu.inc).
_ENF, _ENS, _EEL, _EPL, _ECR, _ETH, _EUL = 1, 2, 5, 6, 7, 8, 9

# Element nodal tensors: point/cell data name, table entry, items per node (the
# first six are the components xx yy zz xy yz xz), stress or engineering strain.
_TENSORS = (
    ("S", _ENS, 6, True),
    ("EPEL", _EEL, 7, False),
    ("EPPL", _EPL, 7, False),
    ("EPCR", _ECR, 7, False),
    ("EPTH", _ETH, 8, False),
)
_TOP = "@top"  # the top surface of a layered shell (the bottom one is unsuffixed)

# Element records read as they are written, one row per element (NaN-padded to
# the longest): their items depend on the element type (the element's
# documentation, "Element Output Definitions").
_RAW_RECORDS = (
    ("EMS", 0),
    ("ENG", 3),
    ("EGR", 4),
    ("EFX", 10),
    ("EMN", 12),
    ("ENL", 14),
    ("EPT", 16),
    ("ECT", 20),
    ("ESV", 23),
)

# Reaction names by DOF code: vector triples, else ``RF_<DOF label>``.
_REACTION_VECTORS = {"RF": (1, 2, 3), "RMOM": (4, 5, 6)}
_DOF_LABELS = {1: "UX", 2: "UY", 3: "UZ", 4: "ROTX", 5: "ROTY", 6: "ROTZ"}
_DOF_LABELS.update({7: "AX", 8: "AY", 9: "AZ", 10: "VX", 11: "VY", 12: "VZ"})
_DOF_LABELS.update(_SCALARS)

_NO_RESULT_CELLS = ("vertex", "line", "line3")


class _Results:
    """One results file: its headers, its model and its per-set pointers."""

    def __init__(self, filename):
        self.filename = str(filename)
        f = self.file = _File(self.filename)
        standard, ptr = f.record(0)
        standard = standard.astype(np.int64)
        if len(standard) < 2 or int(standard[0]) != 12:
            _fail(
                "not a MAPDL results file (its standard header does not name "
                "file 12)"
            )
        h = f.ints(ptr)
        self.nnod = _at(h, 2)
        resmax, nsets = _at(h, 3), _at(h, 8)
        self.kan = _at(h, 7)
        self.cs_els, self.n_sectors, self.cs_cord = _at(h, 18), _at(h, 20), _at(h, 21)
        self.cs_nds = _at(h, 31)
        self.sparse_ens = _at(h, 39) != 0  # else ENS holds 11 items per node
        self.global_nnod = _at(h, 48)
        self.ptr_gnod = _pointer(h, 49, 50)
        self.distributed = self.global_nnod not in (0, self.nnod)
        self.neqv = f.ints(_pointer(h, 14, 45))[: self.nnod]
        self.sets = []
        if nsets:
            dsi = f.ints(_pointer(h, 10, 40))
            tim, _ = f.record(_pointer(h, 11, 41))
            lsp = f.ints(_pointer(h, 12, 42))
            for i in range(nsets):
                lo = int(dsi[i]) & 0xFFFFFFFF
                hi = int(dsi[resmax + i]) if resmax + i < len(dsi) else 0
                step = tuple(int(v) for v in lsp[3 * i : 3 * i + 3])
                self.sets.append((lo | hi << 32, float(tim[i]), step))
        self.geometry = _pointer(h, 15, 46)
        # A cyclic model's harmonic index per set; before v18 the second set of
        # a mode pair repeats the positive index: negated, as later versions
        # write it.
        self.harmonic = []
        ptr_cyc = _pointer(h, 16, 43)
        if ptr_cyc and nsets:
            cyc = f.ints(ptr_cyc)
            self.harmonic = [int(v) for v in cyc[: len(self.sets)]]
            if not any(v < -1 for v in self.harmonic):
                for k in range(len(self.harmonic) - 1):
                    if _close(self.sets[k][1], self.sets[k + 1][1]):
                        self.harmonic[k + 1] = -self.harmonic[k + 1]

    def global_nodes(self):
        """The node numbers of the whole distributed model (main file only)."""
        return self.file.ints(self.ptr_gnod)

    def deck(self):
        """The model as the ``.cdb`` builder takes it, with each node's rotation
        angles and each element type's result layout."""
        f = self.file
        g = f.ints(self.geometry)
        maxety, nnod, nelm = _at(g, 1), _at(g, 3), _at(g, 4)
        ptr_ety = _pointer(g, 20, 21, 6)
        ptr_loc = _pointer(g, 26, 27, 8)
        ptr_eid = _pointer(g, 28, 29, 10)
        self.map_flag = _at(g, 64)
        self.ptr_csy = _pointer(g, 24, 25, 9)
        self.max_csy = _at(g, 5)

        # Element types: an index record of offsets from ptrETY (in the record
        # after it since 2021R1's mapFlag), each to one type's description.
        table, after = f.record(ptr_ety)
        if self.map_flag:
            offsets = f.ints(after)
        else:
            offsets = [int(v) for v in table[:maxety] if v]
        routine, keyopt, nodelm, layout = {}, {}, {}, {}
        for off in offsets:
            info = f.ints(ptr_ety + int(off))
            if len(info) < 2:
                continue
            slot = int(info[0])
            routine[slot] = int(info[1])
            keyopt[slot] = {1: _at(info, 2)}
            nodelm[slot] = _at(info, 60)
            # Nodes with element nodal forces (item 63) and with stresses and
            # strains (item 94: the corners). SHELL181/281 with KEYOPT(8) = 0
            # store the bottom surface, then the top one.
            layers = 2 if routine[slot] in (181, 281) and _at(info, 9) == 0 else 1
            layout[slot] = (_at(info, 62), _at(info, 93), layers)

        node_ids, coords, angles = [], [], []
        ptr = ptr_loc
        for _ in range(nnod):
            values, ptr = f.record(ptr)
            values = values.astype(np.float64)
            node_ids.append(int(values[0]))
            coords.append(
                [float(values[k]) if k < len(values) else 0.0 for k in (1, 2, 3)]
            )
            angles.append(
                [float(values[k]) if k < len(values) else 0.0 for k in (4, 5, 6)]
            )

        elements = []
        if nelm:
            index, _ = f.record(ptr_eid)
            offsets = np.frombuffer(index.astype("<i4").tobytes(), "<i8")[:nelm]
            for off in offsets:
                e = f.ints(ptr_eid + int(off))
                if len(e) < 10:
                    _fail(
                        f"the element record at word {ptr_eid + int(off)} is truncated"
                    )
                slot = int(e[1])
                nodes = e[10:]
                if nodelm.get(slot, 0) > 0:
                    nodes = nodes[: nodelm[slot]]
                elements.append(
                    {
                        "mat": int(e[0]),
                        "slot": slot,
                        "real": int(e[2]),
                        "secnum": int(e[3]),
                        "id": int(e[8]),
                        "nodes": [int(v) for v in nodes],
                    }
                )

        components = []
        ptr = _pointer(g, 50, 51)
        for _ in range(_at(g, 48) if ptr else 0):
            table, ptr = f.record(ptr)
            table = table.astype(np.int64)
            if len(table) < 9 or int(table[0]) not in (1, 2):
                continue
            components.append(
                (_name(table[1:9]), int(table[0]) == 1, _ranges(table[9:]))
            )

        return {
            "routine": routine,
            "keyopt": keyopt,
            "layout": layout,
            "node_ids": node_ids,
            "coords": coords,
            "angles": angles,
            "elements": elements,
            "components": components,
        }

    def coordinate_system(self, number):
        """``(axes, origin)`` of local coordinate system ``number``: the rows of
        ``axes`` are its X, Y and Z axes in global coordinates."""
        f = self.file
        if self.ptr_csy:
            table, after = f.record(self.ptr_csy)
            if self.map_flag:
                offsets = f.ints(after)
            else:
                offsets = [int(v) for v in table[: self.max_csy]]
            for off in offsets:
                if not off:
                    continue
                data, _ = f.record(self.ptr_csy + int(off))
                data = data.astype(np.float64)
                if len(data) >= 22 and int(data[21]) == number:
                    return data[:9].reshape(3, 3), data[9:12].copy()
        _fail(f"coordinate system {number} (the cyclic axis) is not in the file")

    def _solution_header(self, index):
        base = self.sets[index][0]
        return base, self.file.ints(base)

    def nodal(self, index):
        """``(node numbers, values, DOF codes)`` of result set ``index``'s nodal
        solution, or ``None``."""
        f = self.file
        base, s = self._solution_header(index)
        nnod, numdof = _at(s, 2), _at(s, 19)
        dofs = [_at(s, 20 + k) for k in range(numdof)]
        sumdof = numdof + _at(s, 97)
        ptr_nsl = _pointer(s, 104, 105, 10)
        if not ptr_nsl or not numdof:
            return None
        values, after = f.record(base + ptr_nsl)
        values = values.astype(np.float64)
        rows = min(nnod, len(values) // sumdof)
        values = values[: rows * sumdof].reshape(rows, sumdof)
        if rows < nnod:
            # Not every node has a solution: the next record lists which.
            which = f.ints(after)[:rows] - 1
            if (
                len(which) < rows
                or (which < 0).any()
                or (which >= len(self.neqv)).any()
            ):
                _fail(f"result set {index + 1} has a corrupt node index")
            numbers = self.neqv[which]
        else:
            numbers = self.neqv[:rows]
        values[np.abs(values) == _UNDEFINED] = np.nan
        return numbers, values, dofs

    def reactions(self, index):
        """``(node numbers, DOF codes, values)`` of result set ``index``'s
        reaction forces, or ``None``."""
        f = self.file
        base, s = self._solution_header(index)
        nrf, numdof = _at(s, 7), _at(s, 19)
        ptr_rf = _pointer(s, 106, 107, 12)
        if not nrf or not ptr_rf or not numdof:
            return None
        dofs = [_at(s, 20 + k) for k in range(numdof)]
        table, after = f.record(base + ptr_rf)
        # (N - 1) * numdof + k: N the node's position in the nodal equivalence
        # table, k the DOF's position in the set's DOF list (1-based).
        table = np.frombuffer(table.astype("<i4").tobytes(), "<i8")[:nrf]
        values, _ = f.record(after)
        values = values.astype(np.float64)[:nrf]
        if len(table) < nrf or len(values) < nrf:
            _fail(f"result set {index + 1} has a truncated reaction record")
        position = (table - 1) // numdof
        k = (table - 1) % numdof
        if (position < 0).any() or (position >= len(self.neqv)).any():
            _fail(f"result set {index + 1} has a corrupt reaction record")
        codes = np.array([dofs[int(j)] for j in k], dtype=np.int64)
        return self.neqv[position], codes, values

    def element_tables(self, index):
        """``(base, offsets)``: each element's pointer table sits at word
        ``base + offsets[i]`` (0: none); ``None`` without an element solution."""
        base, s = self._solution_header(index)
        ptr = _pointer(s, 118, 119, 11)
        if not ptr:
            return None
        values, _ = self.file.record(base + ptr)
        offsets = np.frombuffer(values.astype("<i4").tobytes(), "<i8")
        return base + ptr, offsets

    def element_record(self, table, entry):
        """The values of entry ``entry`` of the pointer table at word
        ``table``: ``None`` when the element has no such record. A negative
        entry ``-n`` stands for a record of ``n`` zeros that is not written."""
        pointers, _ = self.file.record(table)
        if entry >= len(pointers):
            return None
        ptr = int(pointers[entry])
        if ptr == 0:
            return None
        if ptr < 0:
            return np.zeros(-ptr)
        values, _ = self.file.record(table + ptr)
        return values.astype(np.float64)

    def numdof(self, index):
        _, s = self._solution_header(index)
        return [_at(s, 20 + k) for k in range(_at(s, 19))]


def _euler_to_global(angles):
    """MAPDL's element rotation (THXY, THYZ, THZX in degrees, 3-1-2) as the
    matrix ``Q = R^T`` taking a tensor from the element system to the global
    one (``Q T Q^T``). Scalar arithmetic in the C++ twin's order, so the two
    engines agree bit for bit."""
    deg = math.acos(-1.0) / 180.0
    c1, c2, c3 = (math.cos(a * deg) for a in angles)
    s1, s2, s3 = (math.sin(a * deg) for a in angles)
    r = [
        c1 * c3 - s1 * s2 * s3,
        s1 * c3 + c1 * s2 * s3,
        -s3 * c2,
        -s1 * c2,
        c1 * c2,
        s2,
        c1 * s3 + s1 * s2 * c3,
        s1 * s3 - c1 * c3 * s2,
        c2 * c3,
    ]
    return [r[0], r[3], r[6], r[1], r[4], r[7], r[2], r[5], r[8]]


def _rotate_tensor(row, q, stress):
    """``Q T Q^T`` for one tensor ``xx yy zz xy yz xz`` (a 6-sequence, in
    place; ``q`` row-major); a strain's shear components are engineering
    strains."""
    shear = 1.0 if stress else 0.5
    xx, yy, zz, xy, yz, xz = (float(v) for v in row)
    t = [
        xx,
        shear * xy,
        shear * xz,
        shear * xy,
        yy,
        shear * yz,
        shear * xz,
        shear * yz,
        zz,
    ]
    qt = [0.0] * 9
    for i in range(3):
        for j in range(3):
            v = 0.0
            for k in range(3):
                v += q[i * 3 + k] * t[k * 3 + j]
            qt[i * 3 + j] = v
    out = [0.0] * 9
    for i in range(3):
        for j in range(3):
            v = 0.0
            for k in range(3):
                v += qt[i * 3 + k] * q[j * 3 + k]
            out[i * 3 + j] = v
    row[0], row[1], row[2] = out[0], out[4], out[8]
    row[3], row[4], row[5] = out[1] / shear, out[5] / shear, out[2] / shear


def _rotate_tensors(values, q, stress):
    for row in values.reshape(-1, 6):
        _rotate_tensor(row, q, stress)


def _rotate_vectors(values, q):
    """``Q v`` for each row of an ``(n, 3)`` array (in place)."""
    x, y, z = values[:, 0].copy(), values[:, 1].copy(), values[:, 2].copy()
    for d in range(3):
        values[:, d] = q[d * 3] * x + q[d * 3 + 1] * y + q[d * 3 + 2] * z


class _Model:
    """The model of one results file, or of a distributed solve's files merged
    by node number (each element lives in one file)."""

    def __init__(self, files, lenient):
        self.files = files
        merged = None
        self.offsets = []  # each file's first element in the merged list
        seen = {}
        for r in files:
            deck = r.deck()
            if merged is None:
                merged = {
                    k: deck[k] for k in ("routine", "keyopt", "layout", "elements")
                }
                merged["elements"] = list(deck["elements"])
                merged["node_ids"], merged["coords"], merged["angles"] = [], [], []
                merged["components"] = {}
                self.offsets.append(0)
            else:
                self.offsets.append(len(merged["elements"]))
                merged["elements"].extend(deck["elements"])
                for key in ("routine", "keyopt", "layout"):
                    for slot, value in deck[key].items():
                        merged[key].setdefault(slot, value)
            for ident, xyz, angle in zip(
                deck["node_ids"], deck["coords"], deck["angles"]
            ):
                if ident not in seen:
                    seen[ident] = len(merged["node_ids"])
                    merged["node_ids"].append(ident)
                    merged["coords"].append(xyz)
                    merged["angles"].append(angle)
            for name, is_nodes, ids in deck["components"]:
                key = (name, is_nodes)
                if key in merged["components"]:
                    have = merged["components"][key]
                    known = set(have)
                    have.extend(i for i in ids if i not in known)
                else:
                    merged["components"][key] = list(ids)
        if len(files) > 1:
            missing = set(int(v) for v in files[0].global_nodes()) - set(seen)
            if missing:
                _fail(
                    f"the partial files hold {len(seen)} of the "
                    f"{len(seen) + len(missing)} nodes of the distributed solve; "
                    "one is missing"
                )
        merged["components"] = [
            (name, is_nodes, ids)
            for (name, is_nodes), ids in merged["components"].items()
        ]
        self.deck = merged
        self.mesh, self.node_index, self.locs = _build(
            merged, lenient, _LABEL, locations=True
        )
        self.angles = np.array(merged["angles"], dtype=np.float64).reshape(-1, 3)

    def points_of(self, numbers):
        """Point indices of node numbers (-1 where the mesh has none)."""
        get = self.node_index.get
        return np.array([get(int(n), -1) for n in numbers], dtype=np.int64)

    # -- nodal solution and reactions ---------------------------------------

    def solution(self, index, wanted):
        mesh = self.mesh
        n_points = len(mesh.points)
        vectors, scalars = {}, {}
        for r in self.files:
            nodal = r.nodal(index)
            if nodal is None:
                continue
            numbers, values, dofs = nodal
            points = self.points_of(numbers)
            keep = points >= 0
            values, points = values[keep], points[keep]
            column = {code: k for k, code in enumerate(dofs)}
            for name, codes in _VECTORS.items():
                if not any(c in column for c in codes) or (
                    wanted is not None and name not in wanted
                ):
                    continue
                vec = np.zeros((len(points), 3))
                for d, c in enumerate(codes):
                    if c in column:
                        vec[:, d] = values[:, column[c]]
                _rotate(vec, self.angles[points])
                out = vectors.setdefault(name, np.full((n_points, 3), np.nan))
                out[points] = vec
            for code, k in column.items():
                if any(code in codes for codes in _VECTORS.values()):
                    continue
                name = _SCALARS.get(code, f"DOF{code}")
                if wanted is not None and name not in wanted:
                    continue
                out = scalars.setdefault(name, np.full(n_points, np.nan))
                out[points] = values[:, k]
        mesh.point_data.update(vectors)
        mesh.point_data.update(scalars)

    def reactions(self, index, wanted):
        mesh = self.mesh
        n_points = len(mesh.points)
        out = {}
        for r in self.files:
            rf = r.reactions(index)
            if rf is None:
                continue
            numbers, codes, values = rf
            points = self.points_of(numbers)
            for name, triple in _REACTION_VECTORS.items():
                if wanted is not None and name not in wanted:
                    continue
                mask = np.isin(codes, triple) & (points >= 0)
                if not mask.any():
                    continue
                vec = out.setdefault(name, np.full((n_points, 3), np.nan))
                # A distributed solve's files each hold their share of a
                # reaction at the nodes they have in common: summed.
                for p, c, v in zip(points[mask], codes[mask], values[mask]):
                    row = vec[p]
                    row[np.isnan(row)] = 0.0
                    row[triple.index(int(c))] += v
            for code in sorted(set(int(c) for c in codes)):
                if any(code in t for t in _REACTION_VECTORS.values()):
                    continue
                name = "RF_" + _DOF_LABELS.get(code, f"DOF{code}")
                if wanted is not None and name not in wanted:
                    continue
                mask = (codes == code) & (points >= 0)
                col = out.setdefault(name, np.full(n_points, np.nan))
                for p, v in zip(points[mask], values[mask]):
                    col[p] = v if np.isnan(col[p]) else col[p] + v
        for name, triple in _REACTION_VECTORS.items():
            if name in out:
                have = ~np.isnan(out[name][:, 0])
                vec = out[name][have]
                _rotate(vec, self.angles[have])
                out[name][have] = vec
        mesh.point_data.update(out)

    # -- element solution ----------------------------------------------------

    def elements(self, index, wanted):
        """Element nodal stresses and strains (cell data per element node, and
        averaged at the corner nodes as point data) and element nodal forces."""
        mesh = self.mesh
        blocks = mesh.cells
        n_points = len(mesh.points)
        kinds = [
            k
            for k in _TENSORS
            if wanted is None or k[0] in wanted or k[0] + _TOP in wanted
        ]
        want_enf = wanted is None or "ENF" in wanted
        raw_kinds = [k for k in _RAW_RECORDS if wanted is None or k[0] in wanted]
        raw = {}
        # Every array is (cells, nodes, components) with the widest block's node
        # count, NaN-padded; flattened point-major below so any writer holds it.
        npc = max((blk.data.shape[1] for blk in blocks), default=0)
        cells, sums, counts = {}, {}, {}
        enf_cells = None
        layout = self.deck["layout"]
        elements = self.deck["elements"]
        for r, first in zip(self.files, self.offsets):
            tables = r.element_tables(index)
            if tables is None:
                continue
            base, offsets = tables
            dofs = r.numdof(index)
            for k, off in enumerate(offsets):
                pos = first + k
                if pos >= len(elements) or off == 0 or self.locs[pos] is None:
                    continue
                b, row, slots = self.locs[pos]
                table = base + int(off)
                for name, entry in raw_kinds:
                    values = r.element_record(table, entry)
                    if values is None or len(values) == 0:
                        continue
                    values = np.array(values, dtype=np.float64)
                    values[np.abs(values) == _UNDEFINED] = np.nan
                    raw.setdefault(name, [[] for _ in blocks])[b].append((row, values))
                if blocks[b].type in _NO_RESULT_CELLS:
                    continue
                e = elements[pos]
                nodfor, nodstr, layers = layout.get(e["slot"], (0, 0, 1))
                rotation = None
                for name, entry, items, stress in kinds:
                    if entry == _ENS and not r.sparse_ens:
                        items = 11
                    values = r.element_record(table, entry)
                    if values is None or nodstr <= 0:
                        continue
                    n_layers = layers
                    if len(values) < nodstr * layers * items:
                        n_layers = 1
                        if len(values) < nodstr * items:
                            continue
                    values = values[: nodstr * n_layers * items]
                    values = values.reshape(nodstr * n_layers, items)[:, :6].copy()
                    values[np.abs(values) == _UNDEFINED] = np.nan
                    if rotation is None:
                        rotation = r.element_record(table, _EUL)
                        if rotation is None:
                            rotation = np.zeros(0)
                    if len(rotation) >= 3:
                        per_node = len(rotation) // 3
                        for i in range(len(values)):
                            j = 0 if per_node == 1 else i % per_node
                            angles = rotation[3 * j : 3 * j + 3]
                            if np.any(angles):
                                q = _euler_to_global(angles)
                                _rotate_tensor(values[i], q, stress)
                    for layer in range(n_layers):
                        key = name + (_TOP if layer else "")
                        if wanted is not None and key not in wanted:
                            continue
                        data = values[layer * nodstr : (layer + 1) * nodstr]
                        target = cells.get(key)
                        if target is None:
                            target = cells[key] = [
                                np.full((len(blk.data), npc, 6), np.nan)
                                for blk in blocks
                            ]
                            sums[key] = np.zeros((n_points, 6))
                            counts[key] = np.zeros(n_points, dtype=np.int64)
                        conn = blocks[b].data[row]
                        for j, slot in enumerate(slots):
                            if slot < nodstr and j < len(conn):
                                target[b][row, j] = data[slot]
                                if np.all(np.isfinite(data[slot])):
                                    sums[key][conn[j]] += data[slot]
                                    counts[key][conn[j]] += 1
                if want_enf and nodfor > 0 and dofs:
                    values = r.element_record(table, _ENF)
                    if values is None or len(values) < nodfor * len(dofs):
                        continue
                    values = values[: nodfor * len(dofs)].reshape(nodfor, len(dofs))
                    if enf_cells is None:
                        enf_cells = [
                            np.full((len(blk.data), npc, len(dofs)), np.nan)
                            for blk in blocks
                        ]
                    if enf_cells[b].shape[2] != len(dofs):
                        continue
                    for j, slot in enumerate(slots):
                        if slot < nodfor:
                            enf_cells[b][row, j] = values[slot]
        for key, target in cells.items():
            mesh.cell_data[key] = [a.reshape(len(a), -1) for a in target]
            mesh.field_data["ansys:layout:" + key] = np.array([npc, 6], dtype=np.int64)
            with np.errstate(invalid="ignore", divide="ignore"):
                mean = sums[key] / counts[key][:, None]
            mean[counts[key] == 0] = np.nan
            mesh.point_data[key] = mean
        if enf_cells is not None:
            mesh.cell_data["ENF"] = [a.reshape(len(a), -1) for a in enf_cells]
            width = enf_cells[0].shape[2]
            mesh.field_data["ansys:layout:ENF"] = np.array([npc, width], dtype=np.int64)
        for name in sorted(raw):
            per_block = raw[name]
            width = max(len(v) for rows in per_block for _, v in rows)
            arrays = []
            for blk, rows in zip(blocks, per_block):
                a = np.full((len(blk.data), width), np.nan)
                for row, v in rows:
                    a[row, : len(v)] = v
                arrays.append(a)
            mesh.cell_data[name] = arrays


def _rotate(vec, angles):
    """Nodal to global: a node's axes are the global ones turned about Z by
    THXY, then about the new X by THYZ, then about the newest Y by THZX, so
    ``v_global = Rz Rx Ry v_nodal``."""
    for row, (xy, yz, zx) in zip(vec, angles):
        if xy == 0.0 and yz == 0.0 and zx == 0.0:
            continue
        v = row.copy()
        c, s = math.cos(math.radians(zx)), math.sin(math.radians(zx))
        v = np.array([c * v[0] + s * v[2], v[1], -s * v[0] + c * v[2]])
        c, s = math.cos(math.radians(yz)), math.sin(math.radians(yz))
        v = np.array([v[0], c * v[1] - s * v[2], s * v[1] + c * v[2]])
        c, s = math.cos(math.radians(xy)), math.sin(math.radians(xy))
        v = np.array([c * v[0] - s * v[1], s * v[0] + c * v[1], v[2]])
        row[:] = v


# -- distributed solves ------------------------------------------------------


def _open(filename):
    """The results file, with its partial siblings when it is the main file of a
    distributed solve (``file0.rst`` beside ``file1.rst`` ...)."""
    main = _Results(filename)
    if not main.distributed:
        return [main]
    path = pathlib.Path(main.filename)
    if not main.ptr_gnod:
        _fail(
            f"a partial result file of a distributed solve ({main.nnod} of "
            f"{main.global_nnod} nodes) that is not the main one; read the file "
            "whose name ends in 0 (it finds the others), or the combined file "
            "(RESCOMBINE)"
        )
    if not path.stem.endswith("0"):
        _fail(
            "the main file of a distributed solve must be named <job>0"
            f"{path.suffix} to find its partial files"
        )
    job = path.stem[:-1]
    files = [main]
    while True:
        sibling = path.with_name(f"{job}{len(files)}{path.suffix}")
        if not sibling.exists():
            break
        partial = _Results(sibling)
        if not partial.distributed or len(partial.sets) != len(main.sets):
            _fail(f"{sibling.name} is not a partial file of the same solve")
        files.append(partial)
    return files


# -- cyclic symmetry -----------------------------------------------------------

_VECTOR_NAMES = set(_VECTORS) | set(_REACTION_VECTORS)
_TENSOR_NAMES = {n + t: stress for n, _, _, stress in _TENSORS for t in ("", _TOP)}


def _axis_rotation(axis, theta):
    """Rodrigues: the rotation by ``theta`` about unit vector ``axis``,
    row-major, in the C++ twin's arithmetic order."""
    x, y, z = (float(v) for v in axis)
    c, s = math.cos(theta), math.sin(theta)
    t = 1.0 - c
    return [
        c + t * x * x,
        t * x * y - s * z,
        t * x * z + s * y,
        t * x * y + s * z,
        c + t * y * y,
        t * y * z - s * x,
        t * x * z - s * y,
        t * y * z + s * x,
        c + t * z * z,
    ]


_RAW_NAMES = {name for name, _ in _RAW_RECORDS}


def _modal(filename, model, index, wanted, lenient):
    """The other half of modal set ``index`` of a cyclic model, as ``(harmonic
    index, point arrays, cell arrays)``: its mode pair (the neighbouring set of
    the same frequency), else the duplicate sector (nodes past csNds and
    elements past csEls, paired with the base sector's in number order), else
    nothing."""
    r = model.files[0]
    harmonic = r.harmonic[index] if index < len(r.harmonic) else 0
    sets = r.sets
    pair = None
    if len(sets) > 1:
        before, after = (index - 1) % len(sets), (index + 1) % len(sets)
        if _close(sets[index][1], sets[before][1]):
            pair = before
        elif _close(sets[index][1], sets[after][1]):
            pair = after
    mesh = model.mesh
    if pair is not None:
        other = _Model(_open(filename), lenient)
        other.solution(pair, wanted)
        other.reactions(pair, wanted)
        other.elements(pair, wanted)
        return harmonic, dict(other.mesh.point_data), dict(other.mesh.cell_data)
    base = sorted((n, p) for n, p in model.node_index.items() if n <= r.cs_nds)
    dup = sorted((n, p) for n, p in model.node_index.items() if n > r.cs_nds)
    if not dup or len(dup) != len(base):
        return harmonic, {}, {}
    base_pts = np.array([p for _, p in base], dtype=np.int64)
    dup_pts = np.array([p for _, p in dup], dtype=np.int64)
    points = {}
    for name, values in mesh.point_data.items():
        if values.dtype != np.float64:
            continue
        pair_values = np.zeros_like(values)
        pair_values[base_pts] = values[dup_pts]
        points[name] = pair_values
    base_cells, dup_cells = [], []
    for e, loc in zip(model.deck["elements"], model.locs):
        if loc is not None:
            (base_cells if e["id"] <= r.cs_els else dup_cells).append(
                (e["id"], loc[0], loc[1])
            )
    if len(base_cells) != len(dup_cells):
        return harmonic, points, {}
    base_cells.sort()
    dup_cells.sort()
    cells = {}
    for name, per_block in mesh.cell_data.items():
        if any(a.dtype != np.float64 for a in per_block):
            continue
        pair_blocks = [np.zeros_like(a) for a in per_block]
        for (_, bb, br), (_, db, dr) in zip(base_cells, dup_cells):
            if per_block[db][dr].shape == pair_blocks[bb][br].shape:
                pair_blocks[bb][br] = per_block[db][dr]
        cells[name] = pair_blocks
    return harmonic, points, cells


def _expand_cyclic(model, dofs, modal=None):
    """The full rotor: the base sector's cells (element numbers up to csEls)
    and their points repeated round the cyclic axis, results rotated with them.
    Coincident nodes on the sector boundaries are not merged. A modal set
    (``modal``) is combined with its other half per sector first, as MAPDL's
    /CYCEXPAND does: scale * (x cos(h theta) - x' sin(h theta))."""
    r = model.files[0]
    mesh = model.mesh
    n = r.n_sectors
    if r.cs_cord > 1:
        axes, origin = r.coordinate_system(r.cs_cord)
        x, y, z = (float(v) for v in axes[2])
        norm = math.sqrt(x * x + y * y + z * z)
        axis = np.array([x / norm, y / norm, z / norm])
    else:
        axis, origin = np.array([0.0, 0.0, 1.0]), np.zeros(3)

    keep = [np.zeros(len(b.data), dtype=bool) for b in mesh.cells]
    for e, loc in zip(model.deck["elements"], model.locs):
        if loc is not None and e["id"] <= r.cs_els:
            keep[loc[0]][loc[1]] = True
    used = np.unique(
        np.concatenate(
            [b.data[k].ravel() for b, k in zip(mesh.cells, keep)]
            + [np.zeros(0, dtype=np.int64)]
        )
    )
    used = used[used >= 0]
    new_point = np.full(len(mesh.points), -1, dtype=np.int64)
    new_point[used] = np.arange(len(used))
    n_pts = len(used)
    old_cell = []  # per kept block: global cell indices in the base mesh
    starts = np.cumsum([0] + [len(b.data) for b in mesh.cells])
    blocks = []
    for bi, (b, k) in enumerate(zip(mesh.cells, keep)):
        if k.any():
            blocks.append((bi, b.type, new_point[b.data[k]]))
            old_cell.append(starts[bi] + np.nonzero(k)[0])
    pi = math.acos(-1.0)
    rotations = [_axis_rotation(axis, 2.0 * pi * i / n) for i in range(n)]
    weight, weight_pair = [1.0] * n, [0.0] * n
    if modal is not None:
        h = modal[0]
        single = h == 0 or 2 * abs(h) == n
        scale = 1.0 / math.sqrt(n) if single else 1.0 / math.sqrt(n / 2.0)
        for i in range(n):
            phase = 2.0 * pi * h * i / n
            weight[i] = scale * math.cos(phase)
            weight_pair[i] = -scale * math.sin(phase)

    def combine(values, pair, i):
        if pair is None:
            return weight[i] * values + 0.0
        return weight[i] * values + weight_pair[i] * pair

    base = mesh.points[used] - origin
    copies = []
    for q in rotations:
        copy = base.copy()
        _rotate_vectors(copy, q)
        copies.append(copy + origin)
    points = np.concatenate(copies)
    cells = [
        CellBlock(t, np.concatenate([conn + i * n_pts for i in range(n)]))
        for _, t, conn in blocks
    ]
    out = Mesh(points, cells)
    out.field_data = dict(mesh.field_data)
    out.field_data["ansys:sectors"] = np.array([n], dtype=np.int64)
    if modal is not None:
        out.field_data["ansys:harmonic_index"] = np.array([modal[0]], dtype=np.int64)

    def spin(values, name, i, last_axis_dofs=None):
        q = rotations[i]
        v = np.array(values, dtype=np.float64, copy=True)
        if name in _TENSOR_NAMES and v.size % 6 == 0:
            _rotate_tensors(v, q, _TENSOR_NAMES[name])
        elif v.shape[-1] == 3 and v.ndim == 2 and name in _VECTOR_NAMES:
            _rotate_vectors(v, q)
        elif last_axis_dofs is not None:
            column = {c: j for j, c in enumerate(last_axis_dofs)}
            flat = v.reshape(-1, len(last_axis_dofs))
            for triple in ((1, 2, 3), (4, 5, 6)):
                if all(c in column for c in triple):
                    idx = [column[c] for c in triple]
                    part = flat[:, idx].copy()
                    _rotate_vectors(part, q)
                    flat[:, idx] = part
        return v

    for name, values in mesh.point_data.items():
        part = values[used]
        if modal is not None and values.dtype == np.float64 and name not in _RAW_NAMES:
            pair = modal[1].get(name)
            pair = None if pair is None else pair[used]
            copies = [spin(combine(part, pair, i), name, i) for i in range(n)]
        else:
            copies = [spin(part, name, i) for i in range(n)]
        out.point_data[name] = np.concatenate(copies)
    for name, per_block in mesh.cell_data.items():
        merged = []
        for bi, _, _ in blocks:
            k = keep[bi]
            part = per_block[bi][k]
            if part.dtype.kind in "iub":
                merged.append(np.concatenate([part] * n))
            else:
                dofs_axis = dofs if name == "ENF" else None
                modal_part = (
                    modal is not None
                    and part.dtype == np.float64
                    and name not in _RAW_NAMES
                )
                if modal_part:
                    pair = modal[2].get(name)
                    pair = None if pair is None else pair[bi][k]
                    copies = [
                        spin(combine(part, pair, i), name, i, dofs_axis)
                        for i in range(n)
                    ]
                else:
                    copies = [spin(part, name, i, dofs_axis) for i in range(n)]
                merged.append(np.concatenate(copies))
        out.cell_data[name] = merged
    out.cell_data["ansys:sector"] = [
        np.repeat(np.arange(n, dtype=np.int64), len(conn)) for _, _, conn in blocks
    ]

    new_cell = np.full(int(starts[-1]), -1, dtype=np.int64)
    new_starts = np.cumsum([0] + [len(c) for c in old_cell])
    for j, idx in enumerate(old_cell):
        new_cell[idx] = new_starts[j] + np.arange(len(idx))
    # A cell's copies sit block by block: sector i of kept block j starts at
    # (its block's first expanded cell) + i * (the block's base cell count).
    expanded_start = np.cumsum([0] + [n * len(c) for c in old_cell])
    regions = []
    for region in mesh.regions:
        if region.kind == "point":
            base_entries = new_point[region.entries]
            base_entries = base_entries[base_entries >= 0]
            entries = np.concatenate([base_entries + i * n_pts for i in range(n)])
        elif region.kind == "cell":
            mapped = new_cell[region.entries]
            mapped = mapped[mapped >= 0]
            block = np.searchsorted(new_starts, mapped, side="right") - 1
            local = mapped - new_starts[block]
            sizes = np.array([len(c) for c in old_cell])
            entries = np.concatenate(
                [expanded_start[block] + i * sizes[block] + local for i in range(n)]
            )
        else:
            continue
        regions.append(
            Region(region.name, region.kind, entries, region.dim, region.tag)
        )
    out.regions = regions
    return out


# -- entry points --------------------------------------------------------------


def time_values(filename):
    """The time (or frequency) of every result set."""
    return [t for _, t, _ in _Results(str(filename)).sets]


def read(
    filename, points_only=False, arrays=None, time_step=0, lenient=False, cyclic=False
):
    files = _open(str(filename))
    main = files[0]
    if cyclic:
        if main.n_sectors <= 1:
            _fail("not a cyclic-symmetry model; read it as ansys_rst")
        if main.kan not in (0, 2):
            _fail(
                "only static and modal cyclic-symmetry results can be expanded "
                f"(this is analysis type {main.kan}); read the base sector as "
                "ansys_rst"
            )
    elif main.n_sectors > 1:
        warn(
            f"{_LABEL}: a cyclic-symmetry model ({main.n_sectors} sectors); only "
            "the base sector is read (ansys_rst_cyclic expands it)"
        )
    model = _Model(files, lenient)
    mesh = model.mesh
    n = len(main.sets)
    if not n:
        if time_step not in (0, -1):
            _fail(f"time step {time_step} is out of range: the file has no result sets")
        return _expand_cyclic(model, []) if cyclic else mesh
    index = time_step + n if time_step < 0 else time_step
    if not 0 <= index < n:
        raise ReadError(
            f"time step {time_step} is out of range: the file has {n} step(s)"
        )
    _, time, (load_step, substep, cumulative) = main.sets[index]
    mesh.field_data["meshio:time"] = np.array([time], dtype=np.float64)
    mesh.field_data["ansys:load_step"] = np.array([load_step], dtype=np.int64)
    mesh.field_data["ansys:substep"] = np.array([substep], dtype=np.int64)
    mesh.field_data["ansys:cumulative"] = np.array([cumulative], dtype=np.int64)
    dofs = main.numdof(index)
    if not points_only:
        wanted = None if arrays is None else set(arrays)
        model.solution(index, wanted)
        model.reactions(index, wanted)
        model.elements(index, wanted)
    if cyclic and main.kan == 2:
        modal = (
            (main.harmonic[index] if index < len(main.harmonic) else 0, {}, {})
            if points_only
            else _modal(str(filename), model, index, wanted, lenient)
        )
        mesh = _expand_cyclic(model, dofs, modal)
    elif cyclic:
        mesh = _expand_cyclic(model, dofs)
    mesh.time_values = [t for _, t, _ in main.sets]
    return mesh

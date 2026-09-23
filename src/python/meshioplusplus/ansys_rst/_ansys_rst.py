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

import numpy as np

from .._common import warn
from .._exceptions import ReadError
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


class _Results:
    """The headers, the model and the per-set pointers of a results file."""

    def __init__(self, filename):
        f = self.file = _File(filename)
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
        self.n_sectors = _at(h, 20)
        global_nnod = _at(h, 48)
        if global_nnod not in (0, self.nnod):
            _fail(
                f"a partial result file of a distributed solve ({self.nnod} of "
                f"{global_nnod} nodes); read the combined file (RESCOMBINE)"
            )
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

    def model(self, lenient):
        f = self.file
        g = f.ints(self.geometry)
        maxety, nnod, nelm = _at(g, 1), _at(g, 3), _at(g, 4)
        ptr_ety = _pointer(g, 20, 21, 6)
        ptr_loc = _pointer(g, 26, 27, 8)
        ptr_eid = _pointer(g, 28, 29, 10)

        # Element types: an index record of offsets from ptrETY (in the record
        # after it since 2021R1's mapFlag), each to one type's description.
        table, after = f.record(ptr_ety)
        if _at(g, 64):
            offsets = f.ints(after)
        else:
            offsets = [int(v) for v in table[:maxety] if v]
        routine, keyopt, nodelm = {}, {}, {}
        for off in offsets:
            info = f.ints(ptr_ety + int(off))
            if len(info) < 2:
                continue
            slot = int(info[0])
            routine[slot] = int(info[1])
            keyopt[slot] = {1: _at(info, 2)}
            nodelm[slot] = _at(info, 60)

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

        deck = {
            "routine": routine,
            "keyopt": keyopt,
            "node_ids": node_ids,
            "coords": coords,
            "elements": elements,
            "components": components,
        }
        mesh, node_index = _build(deck, lenient, _LABEL)
        return mesh, node_index, np.array(angles, dtype=np.float64).reshape(-1, 3)

    def solution(self, index, mesh, node_index, angles, wanted):
        """Point data of result set ``index``."""
        f = self.file
        base, _, _ = self.sets[index]
        s = f.ints(base)
        nnod, numdof = _at(s, 2), _at(s, 19)
        dofs = [_at(s, 20 + k) for k in range(numdof)]
        sumdof = numdof + _at(s, 97)
        ptr_nsl = _pointer(s, 104, 105, 10)
        if not ptr_nsl or not numdof:
            return
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
        points = np.array([node_index.get(int(n), -1) for n in numbers], dtype=np.int64)
        keep = points >= 0
        values, points = values[keep], points[keep]
        values[np.abs(values) == _UNDEFINED] = np.nan

        column = {code: k for k, code in enumerate(dofs)}
        n_points = len(mesh.points)
        for name, codes in _VECTORS.items():
            if not any(c in column for c in codes) or (
                wanted is not None and name not in wanted
            ):
                continue
            vec = np.zeros((len(points), 3))
            for d, c in enumerate(codes):
                if c in column:
                    vec[:, d] = values[:, column[c]]
            _rotate(vec, angles[points])
            out = np.full((n_points, 3), np.nan)
            out[points] = vec
            mesh.point_data[name] = out
        for code, k in column.items():
            if any(code in codes for codes in _VECTORS.values()):
                continue
            name = _SCALARS.get(code, f"DOF{code}")
            if wanted is not None and name not in wanted:
                continue
            out = np.full(n_points, np.nan)
            out[points] = values[:, k]
            mesh.point_data[name] = out


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


def time_values(filename):
    """The time (or frequency) of every result set."""
    return [t for _, t, _ in _Results(str(filename)).sets]


def read(filename, points_only=False, arrays=None, time_step=0, lenient=False):
    results = _Results(str(filename))
    if results.n_sectors > 1:
        warn(
            f"{_LABEL}: a cyclic-symmetry model ({results.n_sectors} sectors); only "
            "the base sector is read"
        )
    mesh, node_index, angles = results.model(lenient)
    n = len(results.sets)
    if not n:
        if time_step not in (0, -1):
            _fail(f"time step {time_step} is out of range: the file has no result sets")
        return mesh
    index = time_step + n if time_step < 0 else time_step
    if not 0 <= index < n:
        raise ReadError(
            f"time step {time_step} is out of range: the file has {n} step(s)"
        )
    mesh.time_values = [t for _, t, _ in results.sets]
    _, time, (load_step, substep, cumulative) = results.sets[index]
    mesh.field_data["meshio:time"] = np.array([time], dtype=np.float64)
    mesh.field_data["ansys:load_step"] = np.array([load_step], dtype=np.int64)
    mesh.field_data["ansys:substep"] = np.array([substep], dtype=np.int64)
    mesh.field_data["ansys:cumulative"] = np.array([cumulative], dtype=np.int64)
    if points_only:
        return mesh
    wanted = None if arrays is None else set(arrays)
    results.solution(index, mesh, node_index, angles, wanted)
    return mesh

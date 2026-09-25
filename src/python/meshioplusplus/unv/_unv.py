"""
I/O for the I-DEAS Universal File format (``.unv`` / ``.uff``): the Python
reference twin of ``src/cpp/src/formats/unv.cpp``.

A universal file is a sequence of datasets, each opened by a ``-1`` line and a
dataset-number line and closed by another ``-1`` line. Read: nodes 2411/781/15
(moved out of Cartesian 2420 coordinate systems), elements 2412/780, permanent
groups 2467/2477/2452/2435 (quadruples) and 2417/2429/2430/2432 (pairs) as
regions, units 164, and results 2414/55/56 plus the functions of 58/58b as the
steps of a sequence. See ``doc/formats/unv.md``.
"""

import math
import struct

import numpy as np

from .. import _provenance
from .._common import warn
from .._exceptions import ReadError
from .._files import open_file
from .._mesh import CellBlock, Mesh
from .._node_order import node_order, node_order_keys
from .._regions import Region

__all__ = ["read", "write", "time_values"]

TIME_KEY = "meshio:time"
_NAN = float("nan")

# UNV node order -> meshio position: meshio_conn[perm[i]] = unv_conn[i], i.e. the
# "from meshio" direction of the "unv" tables in meshioplusplus/_node_order.py.
# Parabolic elements list their mid-side nodes "sandwiched" between the corners of
# each ring; the solids list the bottom ring, then the vertical mid-edges, then the
# top ring.
_PERM = {
    cell_type: list(node_order(fmt, cell_type).from_meshio)
    for fmt, cell_type in node_order_keys()
    if fmt == "unv"
}

_BEAMS = {11, 21, 22, 23, 24, 25}
# FE descriptor id -> (meshio type, node count)
_FE = {}
for _ids, _t, _n in (
    ((41, 51, 61, 74, 81, 91), "triangle", 3),
    ((42, 52, 62, 72, 82, 92), "triangle6", 6),
    ((44, 54, 64, 71, 84, 94, 122), "quad", 4),
    ((45, 55, 65, 75, 85, 95), "quad8", 8),
    ((111,), "tetra", 4),
    ((118,), "tetra10", 10),
    ((112,), "wedge", 6),
    ((113,), "wedge15", 15),
    ((115,), "hexahedron", 8),
    ((116,), "hexahedron20", 20),
    ((119, 312), "pyramid", 5),
    ((114,), "pyramid13", 13),
):
    for _i in _ids:
        _FE[_i] = (_t, _n)

# meshio type -> (descriptor id, is_beam) used on write
_MESHIO_TO_UNV = {
    "line": (21, True),
    "line3": (24, True),
    "triangle": (91, False),
    "triangle6": (92, False),
    "quad": (94, False),
    "quad8": (95, False),
    "quad9": (95, False),
    "tetra": (111, False),
    "tetra10": (118, False),
    "wedge": (112, False),
    "wedge15": (113, False),
    "hexahedron": (115, False),
    "hexahedron20": (116, False),
    "pyramid": (312, False),
    "pyramid13": (114, False),
}

_DIM = {
    "line": 1,
    "line3": 1,
    "triangle": 2,
    "triangle6": 2,
    "quad": 2,
    "quad8": 2,
    "quad9": 2,
}

# Permanent-group datasets: True = (type, tag, leaf, component) quadruples.
_GROUP_LAYOUT = {
    2417: False,
    2429: False,
    2430: False,
    2432: False,
    2435: True,
    2452: True,
    2467: True,
    2477: True,
}

_RESULT_TYPE_NAMES = {
    1: "general",
    2: "stress",
    3: "strain",
    4: "element_force",
    5: "temperature",
    6: "heat_flux",
    7: "strain_energy",
    8: "displacement",
    9: "reaction_force",
    10: "kinetic_energy",
    11: "velocity",
    12: "acceleration",
    13: "strain_energy_density",
    14: "kinetic_energy_density",
    15: "hydrostatic_pressure",
    16: "heat_gradient",
    17: "code_checking_value",
    18: "pressure_coefficient",
}

_FUNCTION_NAMES = [
    "function",
    "time_response",
    "auto_spectrum",
    "cross_spectrum",
    "frf",
    "transmissibility",
    "coherence",
    "auto_correlation",
    "cross_correlation",
    "psd",
    "esd",
    "pdf",
    "spectrum",
    "cumulative_frequency_distribution",
    "peaks_valley",
    "stress_cycles",
    "strain_cycles",
    "orbit",
    "mode_indicator_function",
    "force_pattern",
    "partial_power",
    "partial_coherence",
    "eigenvalue",
    "eigenvector",
    "shock_response_spectrum",
    "fir_filter",
    "multiple_coherence",
    "order_function",
]

_DIR_NAMES = ["", "x", "y", "z", "rx", "ry", "rz"]

# Symmetric tensor: file Sxx Sxy Syy Sxz Syz Szz <-> meshio xx yy zz xy yz zx.
_SYM_TO_MESHIO = [0, 2, 5, 1, 4, 3]
_SYM_FROM_MESHIO = [0, 3, 1, 5, 4, 2]


def _cell_type(fe_id, n):
    if fe_id in _BEAMS:
        return {2: "line", 3: "line3"}.get(n)
    entry = _FE.get(fe_id)
    if entry is None:
        return None
    if entry[0] == "quad8" and n == 9:
        return "quad9"
    return entry[0] if entry[1] == n else None


def _function_name(t):
    if 0 <= t < len(_FUNCTION_NAMES):
        return _FUNCTION_NAMES[t]
    return f"function_{t}"


def _data_char(ncomp):
    return {1: 1, 3: 2, 6: 4, 9: 5}.get(ncomp, 0)


def _tensor_to_meshio(char, ncomp, vals):
    if char == 4 and ncomp == 6:
        return vals[:, _SYM_TO_MESHIO]
    if char == 5 and ncomp == 9:
        return vals.reshape(-1, 3, 3).transpose(0, 2, 1).reshape(-1, 9)
    return vals


def _tensor_from_meshio(ncomp, vals):
    if ncomp == 6:
        return vals[:, _SYM_FROM_MESHIO]
    if ncomp == 9:
        return vals.reshape(-1, 3, 3).transpose(0, 2, 1).reshape(-1, 9)
    return vals


# ---------------------------------------------------------------------------
# Text helpers
# ---------------------------------------------------------------------------


def _int(tok):
    try:
        return int(tok)
    except ValueError:
        raise ReadError(f"UNV: expected an integer, got '{tok}'") from None


def _ints(line):
    # Code_Aster ends integer records with a `%` comment ("1  % NOEUD N1").
    out = []
    for t in line.split():
        if t.startswith("%"):
            break
        out.append(_int(t))
    return out


def _real(tok):
    try:
        return float(tok.replace("D", "E").replace("d", "e"))
    except ValueError:
        raise ReadError(f"UNV: expected a real number, got '{tok}'") from None


def _reals(line):
    # Fixed-width fields run together when a value fills its field, so a sign that
    # does not follow an exponent letter also starts a new number.
    out = []
    for tok in line.split():
        start = 0
        for k in range(1, len(tok)):
            if tok[k] in "+-" and tok[k - 1] not in "EeDd":
                out.append(_real(tok[start:k]))
                start = k
        out.append(_real(tok[start:]))
    return out


def _take_reals(lines, k, count):
    vals = []
    while len(vals) < count and k < len(lines):
        vals += _reals(lines[k])
        k += 1
    if len(vals) < count:
        raise ReadError("UNV: dataset ends inside a data record")
    return vals[:count], k


def _take_ints(lines, k, count):
    vals = []
    while len(vals) < count and k < len(lines):
        vals += _ints(lines[k])
        k += 1
    if len(vals) < count:
        raise ReadError("UNV: dataset ends inside an integer record")
    return vals[:count], k


def _is_none(name):
    return name == "" or name.upper() == "NONE"


# ---------------------------------------------------------------------------
# Dataset splitting (byte-safe: 58b carries raw binary data)
# ---------------------------------------------------------------------------


class _Dataset:
    __slots__ = ("id", "binary", "byte_order", "fp_format", "lines", "blob")

    def __init__(self):
        self.id = 0
        self.binary = False
        self.byte_order = 1
        self.fp_format = 2
        self.lines = []
        self.blob = b""


def _58b_bytes(lines):
    """Bytes of a 58b data block, from record 7; -1 when it cannot be read."""
    if len(lines) < 7:
        return -1
    r7 = lines[6].split()
    if len(r7) < 3:
        return -1
    try:
        dtype, npts, even = int(r7[0]), int(r7[1]), int(r7[2]) == 1
    except ValueError:
        return -1
    if npts < 0 or dtype not in (2, 4, 5, 6):
        return -1
    per = (2 if dtype in (5, 6) else 1) + (0 if even else 1)
    return npts * per * (8 if dtype in (4, 6) else 4)


def _split_datasets(data):
    out = []
    pos, n = 0, len(data)
    warned_58b_size = False

    def next_line():
        nonlocal pos
        eol = data.find(b"\n", pos)
        end = n if eol < 0 else eol
        line = data[pos:end]
        if line.endswith(b"\r"):
            line = line[:-1]
        pos = n if eol < 0 else eol + 1
        return line.decode("latin-1")

    while pos < n:
        if next_line().strip() != "-1":
            continue
        header = ""
        while True:
            if pos >= n:
                return out
            header = next_line()
            if header.strip() not in ("-1", ""):
                break
        tokens = header.split()
        ds = _Dataset()
        ident = tokens[0]
        if ident[-1:] in ("b", "B"):
            ds.binary = True
            ident = ident[:-1]
        ds.id = _int(ident)
        if ds.binary:
            if len(tokens) < 5:
                raise ReadError(f"UNV: malformed binary dataset header '{header}'")
            ds.byte_order = _int(tokens[1])
            ds.fp_format = _int(tokens[2])
            n_ascii, n_bytes = _int(tokens[3]), _int(tokens[4])
            for _ in range(n_ascii):
                if pos >= n:
                    break
                ds.lines.append(next_line())
            # A 58b's size follows from its record 7; some writers (pyuff among them)
            # declare half of it for complex data, so the record wins.
            expected = _58b_bytes(ds.lines) if ds.id == 58 else -1
            if expected >= 0 and expected != n_bytes:
                # pyuff declares exactly half the size of complex data; any other
                # mismatch deserves a warning.
                if 2 * n_bytes != expected and not warned_58b_size:
                    warn(
                        f"UNV: dataset 58b declares {n_bytes} bytes but its record 7 "
                        f"describes {expected}; using the record"
                    )
                    warned_58b_size = True
                n_bytes = expected
            if n_bytes < 0 or pos + n_bytes > n:
                raise ReadError(f"UNV: binary dataset {ds.id} is truncated")
            ds.blob = data[pos : pos + n_bytes]
            pos += n_bytes
            # The closing "-1" follows the binary data directly (no newline first).
            if next_line().strip() != "-1":
                warn(f"UNV: binary dataset {ds.id} is not closed by '-1'")
        else:
            while pos < n:
                line = next_line()
                if line.strip() == "-1":
                    break
                ds.lines.append(line)
        out.append(ds)
    return out


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------


class _File:
    def __init__(self):
        self.nodes = []  # (label, cs, [x, y, z][:ncoords])
        self.elements = []  # (label, fe_id, pid, mid, [node labels])
        self.groups = []  # (number, name, [node labels], [element labels])
        self.cs = {}  # label -> (type, 4x3 matrix)
        self.units = None  # (code, [4 factors])
        self.results = []  # dict per 2414/55/56 dataset
        self.functions = []  # dict per 58 function


def _parse_nodes(ds, f):
    lines = ds.lines
    k = 0
    while k < len(lines):
        if not lines[k].strip():
            k += 1
            continue
        if ds.id == 15:
            line = lines[k]
            k += 1
            ints = line[:40].split()
            if len(ints) < 2:
                raise ReadError("UNV: malformed dataset-15 node record")
            xs = _reals(line[40:]) if len(line) > 40 else []
            f.nodes.append((_int(ints[0]), _int(ints[1]), xs[:3]))
        else:
            r1 = _ints(lines[k])
            if not r1 or k + 1 >= len(lines):
                raise ReadError(f"UNV: malformed node record in dataset {ds.id}")
            xs = _reals(lines[k + 1])
            f.nodes.append((r1[0], r1[1] if len(r1) > 1 else 0, xs[:3]))
            k += 2


def _parse_elements(ds, f):
    lines = ds.lines
    k = 0
    while k < len(lines):
        if not lines[k].strip():
            k += 1
            continue
        r1 = _ints(lines[k])
        k += 1
        if ds.id == 780:
            if len(r1) < 8:
                raise ReadError("UNV: malformed dataset-780 element record")
            label, fe_id, pid, mid, nn = r1[0], r1[1], r1[3], r1[5], r1[7]
        else:
            if len(r1) < 6:
                raise ReadError("UNV: malformed dataset-2412 element record")
            label, fe_id, pid, mid, nn = r1[0], r1[1], r1[2], r1[3], r1[5]
        if nn < 0:
            raise ReadError(f"UNV: negative node count on element {label}")
        if fe_id in _BEAMS:
            k += 1  # orientation node and cross sections
        nodes, k = _take_ints(lines, k, nn)
        f.elements.append((label, fe_id, pid, mid, nodes))


def _parse_groups(ds, quad, f):
    lines = ds.lines
    stride = 4 if quad else 2
    k = 0
    while k < len(lines):
        if not lines[k].strip():
            k += 1
            continue
        r1 = _ints(lines[k])
        k += 1
        if len(r1) < 8:
            raise ReadError(f"UNV: malformed group record in dataset {ds.id}")
        if r1[7] < 0:
            raise ReadError(f"UNV: negative entity count in dataset {ds.id}")
        name = lines[k].strip() if k < len(lines) else ""
        k += 1
        vals, k = _take_ints(lines, k, stride * r1[7])
        nodes, elems = [], []
        for e in range(r1[7]):
            etype, tag = vals[stride * e], vals[stride * e + 1]
            if etype == 7:
                nodes.append(tag)
            elif etype == 8:
                elems.append(tag)
        f.groups.append((r1[0], name, nodes, elems))


def _parse_units(ds, f):
    if not ds.lines:
        return
    r1 = ds.lines[0].split()
    if not r1:
        return
    factors = []
    k = 1
    while k < len(ds.lines) and len(factors) < 4:
        factors += _reals(ds.lines[k])
        k += 1
    factors = (factors + [0.0] * 4)[:4]
    # Record 1 is I10, 20A1, I10: the description may touch the code ("5mm").
    f.units = (_int(ds.lines[0][:10].strip() or r1[0]), factors)


def _parse_cs(ds, f):
    lines = ds.lines
    k = 2  # part UID, part name
    while k < len(lines):
        if not lines[k].strip():
            k += 1
            continue
        r3 = _ints(lines[k])
        if len(r3) < 2 or k + 2 >= len(lines):
            break
        m, k = _take_reals(lines, k + 2, 12)
        f.cs[r3[0]] = (r3[1], np.array(m).reshape(4, 3))


def _key_2414(analysis, ints, reals):
    def i(n):
        return ints[n] if n < len(ints) else 0

    def r(n):
        return reals[n] if n < len(reals) else 0.0

    if analysis == 2:
        return (analysis, i(5), r(1))
    if analysis in (3, 7):
        return (analysis, i(5), abs(r(7)) / (2.0 * math.pi))
    if analysis in (4, 9):
        return (analysis, i(6), r(0))
    if analysis == 5:
        return (analysis, i(7), r(1))
    if analysis == 6:
        return (analysis, i(5), r(2))
    return (analysis, i(4), r(0) if r(0) != 0.0 else float(i(4)))


def _key_55(analysis, ints, reals):
    def i(n):
        return ints[n] if n < len(ints) else 0

    def r(n):
        return reals[n] if n < len(reals) else 0.0

    if analysis in (2, 4, 5):
        return (analysis, i(1), r(0))
    if analysis in (3, 7):
        return (analysis, i(1), abs(r(1)) / (2.0 * math.pi))
    if analysis == 6:
        return (analysis, i(0), r(0))
    return (analysis, i(0), r(0) if r(0) != 0.0 else float(i(0)))


def _parse_result(ds, f):
    lines = ds.lines
    legacy = False
    if ds.id == 2414:
        if len(lines) < 13:
            raise ReadError("UNV: dataset 2414 has a truncated header")
        name = lines[1].strip()
        loc = _ints(lines[2])
        location = loc[0] if loc else 1
        r9 = _ints(lines[8])
        if len(r9) < 6:
            raise ReadError("UNV: dataset 2414 record 9 needs six integers")
        analysis, char, rtype, dtype, ndv = r9[1], r9[2], r9[3], r9[4], max(r9[5], 0)
        ints = _ints(lines[9]) + _ints(lines[10])
        reals = _reals(lines[11]) + _reals(lines[12])
        key = _key_2414(analysis, ints, reals)
        k = 13
    else:
        # 55 (data at nodes), 56 (data at elements): ID lines, record 6, then records
        # 7-8 sized by their own NINT/NRV counts.
        if len(lines) < 7:
            raise ReadError(f"UNV: dataset {ds.id} has a truncated header")
        name = lines[0].strip()
        location = {55: 1, 56: 2}.get(ds.id, 3)
        r6 = _ints(lines[5])
        if len(r6) < 6:
            raise ReadError(f"UNV: dataset {ds.id} record 6 needs six integers")
        analysis, char, rtype, dtype, ndv = r6[1], r6[2], r6[3], r6[4], max(r6[5], 0)
        k = 6
        head = _ints(lines[k])
        if len(head) >= 2 and head[0] == 0 and head[1] == 0:
            # meshio++ <= 15.5 wrote four zero lines here (8 ints, 2 ints, 2 x 6 reals),
            # and its "57" held element data; a valid record 7 has NINT >= 1.
            warn(f"UNV: dataset {ds.id} uses the meshio++ <= 15.5 header layout")
            legacy = True
            k += 4
            key = (analysis, 0, 0.0)
        else:
            if len(head) < 2:
                raise ReadError(f"UNV: dataset {ds.id} record 7 needs NINT and NRV")
            nint, nrv = max(head[0], 0), max(head[1], 0)
            ints = head[2:]
            k += 1
            while len(ints) < nint and k < len(lines):
                ints += _ints(lines[k])
                k += 1
            reals, k = _take_reals(lines, k, nrv)
            key = _key_55(analysis, ints, reals)
    if ds.id == 57 and legacy:
        location = 2
    if _is_none(name):
        name = _RESULT_TYPE_NAMES.get(rtype, "unv:field")
    if location not in (1, 2):
        warn(
            f"UNV: skipping '{name}' (dataset {ds.id}, location {location}): only data "
            "at nodes and on elements is read"
        )
        return
    if dtype not in (1, 2, 4, 5, 6):
        warn(f"UNV: skipping '{name}' (dataset {ds.id}): unknown data type {dtype}")
        return
    if ndv == 0:
        return
    cplx = dtype in (5, 6)
    width = 2 if cplx else 1
    values = []
    warned_layers = False
    while k < len(lines):
        if not lines[k].strip():
            k += 1
            continue
        rec = _ints(lines[k])
        k += 1
        count = ndv
        if location == 2 and len(rec) >= 2 and rec[1] > 0:
            count = rec[1]  # NDVAL of this element
        vals, k = _take_reals(lines, k, count * width)
        if count != ndv:
            if not warned_layers:
                warn(
                    f"UNV: '{name}' has {count} values on an element but {ndv} per "
                    "component set; keeping the first"
                )
            warned_layers = True
            vals = (vals + [_NAN] * (ndv * width))[: ndv * width]
        values.append((rec[0], vals))
    f.results.append(
        {
            "location": location,
            "name": name,
            "char": char,
            "ncomp": ndv,
            "complex": cplx,
            "key": key,
            "values": values,
        }
    )


def _parse_function(ds, f):
    lines = ds.lines
    if len(lines) < 11:
        raise ReadError("UNV: dataset 58 has a truncated header")
    rec6 = lines[5]

    def fld(a, n):
        s = rec6[a : a + n].strip()
        return _int(s) if s else 0

    if len(rec6) >= 80:
        ftype, load_case = fld(0, 5), fld(20, 10)
        rsp_node, rsp_dir = fld(41, 10), fld(51, 4)
        ref_node, ref_dir = fld(66, 10), fld(76, 4)
    else:
        t = rec6.split()
        if len(t) < 10:
            raise ReadError("UNV: dataset 58 record 6 is malformed")
        ftype, load_case = _int(t[0]), _int(t[3])
        rsp_node, rsp_dir = _int(t[5]), _int(t[6])
        ref_node, ref_dir = _int(t[8]), _int(t[9])
    r7 = lines[6].split()
    if len(r7) < 3:
        raise ReadError("UNV: dataset 58 record 7 is malformed")
    ord_type, npts, even = _int(r7[0]), max(_int(r7[1]), 0), _int(r7[2]) == 1
    r7_reals = []
    for tok in r7[3:]:
        r7_reals += _reals(tok)
    r7_reals = (r7_reals + [0.0] * 3)[:3]
    r8 = lines[7].split()
    abscissa_type = _int(r8[0]) if r8 else 0
    if ord_type not in (2, 4, 5, 6):
        raise ReadError(
            f"UNV: dataset 58 ordinate data type {ord_type} is not supported"
        )
    cplx = ord_type in (5, 6)
    dbl = ord_type in (4, 6)
    per = (2 if cplx else 1) + (0 if even else 1)
    count = npts * per
    if ds.binary:
        if ds.fp_format != 2:
            raise ReadError(
                f"UNV: dataset 58b floating-point format {ds.fp_format} is not IEEE 754"
            )
        w = 8 if dbl else 4
        if len(ds.blob) < count * w:
            raise ReadError(
                "UNV: dataset 58b holds fewer bytes than its header declares"
            )
        order = ">" if ds.byte_order == 2 else "<"
        values = list(
            struct.unpack(f"{order}{count}{'d' if dbl else 'f'}", ds.blob[: count * w])
        )
    else:
        values, _ = _take_reals(lines, 11, count)
    v = np.array(values, dtype=float).reshape(npts, per) if npts else np.empty((0, per))
    o = 0
    if even:
        x = r7_reals[0] + np.arange(npts) * r7_reals[1]
    else:
        x = v[:, 0].copy()
        o = 1
    re = v[:, o].copy()
    im = v[:, o + 1].copy() if cplx else None
    f.functions.append(
        {
            "type": ftype,
            "load_case": load_case,
            "rsp_node": rsp_node,
            "rsp_dir": rsp_dir,
            "ref_node": ref_node,
            "ref_dir": ref_dir,
            "abscissa_type": abscissa_type,
            "complex": cplx,
            "x": [float(a) for a in x],
            "re": re,
            "im": im,
        }
    )


def _parse(filename):
    with open_file(filename, "rb") as fh:
        data = fh.read()
    if isinstance(data, str):
        data = data.encode("latin-1")
    f = _File()
    for ds in _split_datasets(data):
        if ds.id in (2411, 781, 15):
            _parse_nodes(ds, f)
        elif ds.id in (2412, 780):
            _parse_elements(ds, f)
        elif ds.id in _GROUP_LAYOUT:
            _parse_groups(ds, _GROUP_LAYOUT[ds.id], f)
        elif ds.id == 164:
            _parse_units(ds, f)
        elif ds.id == 2420:
            _parse_cs(ds, f)
        elif ds.id in (2414, 55, 56, 57):
            _parse_result(ds, f)
        elif ds.id == 58:
            _parse_function(ds, f)
    return f


# ---------------------------------------------------------------------------
# Steps
# ---------------------------------------------------------------------------


def _group_functions(f):
    groups = []  # dicts: name, analysis, x, functions
    kinds = []  # [type, load_case, ref_node, ref_axis, [group indices]]
    for fi, fn in enumerate(f.functions):
        ref_axis = abs(fn["ref_dir"])
        kind = next(
            (
                k
                for k in kinds
                if k[0] == fn["type"]
                and k[1] == fn["load_case"]
                and k[2] == fn["ref_node"]
                and k[3] == ref_axis
            ),
            None,
        )
        if kind is None:
            kind = [fn["type"], fn["load_case"], fn["ref_node"], ref_axis, []]
            kinds.append(kind)
        g = next((c for c in kind[4] if groups[c]["x"] == fn["x"]), None)
        if g is None:
            if kind[4]:
                warn(
                    f"UNV: {_function_name(fn['type'])} functions of one kind do not "
                    "share an abscissa grid; each grid becomes its own steps"
                )
            at = fn["abscissa_type"]
            groups.append(
                {
                    "name": "",
                    "analysis": 5 if at == 18 else (4 if at == 17 else 0),
                    "x": fn["x"],
                    "functions": [],
                }
            )
            g = len(groups) - 1
            kind[4].append(g)
        groups[g]["functions"].append(fi)
    type_count = {}
    for k in kinds:
        n = _function_name(k[0])
        type_count[n] = type_count.get(n, 0) + 1
    used = set()
    for k in kinds:
        base = _function_name(k[0])
        if type_count[base] > 1:
            base += f"_ref{k[2]}{_DIR_NAMES[k[3]] if k[3] <= 6 else ''}"
        for g in k[4]:
            name, s = base, 2
            while name in used:
                name = f"{base}_{s}"
                s += 1
            used.add(name)
            groups[g]["name"] = name
    return groups


def _steps(f, groups):
    steps = []  # [key, [result indices], [(group, sample)]]
    index = {}

    def step_of(key):
        if key not in index:
            index[key] = len(steps)
            steps.append([key, [], []])
        return steps[index[key]]

    for r, res in enumerate(f.results):
        step_of(res["key"])[1].append(r)
    for g, grp in enumerate(groups):
        for n, x in enumerate(grp["x"]):
            step_of((grp["analysis"], n + 1, x))[2].append((g, n))
    return steps


def _resolve_step(time_step, count):
    if count == 0:
        if time_step in (0, -1):
            return None
        raise ReadError(
            f"meshio++: time step {time_step} requested, but this file carries no time steps"
        )
    step = time_step + count if time_step < 0 else time_step
    if not 0 <= step < count:
        raise ReadError(
            f"meshio++: time step {time_step} is out of range: this file has {count} "
            + ("step" if count == 1 else "steps")
        )
    return step


def _nan(n, ncomp):
    return np.full((n,) if ncomp == 1 else (n, ncomp), _NAN)


# ---------------------------------------------------------------------------
# Reading
# ---------------------------------------------------------------------------


def read(filename, points_only=False, arrays=None, time_step=0):
    f = _parse(filename)

    # Points, moved into the global system where a node was defined in a local one.
    node_index = {}
    coords = []
    dim = max(len(f.nodes[0][2]), 1) if f.nodes else 3
    warned_cs = set()
    for label, cs, xs in f.nodes:
        x = (list(xs) + [0.0, 0.0, 0.0])[:3]
        sysdef = f.cs.get(cs)
        if sysdef is not None:
            ctype, m = sysdef
            if ctype == 0:
                x = list(m[:3] @ np.array(x) + m[3])
                if dim < 3 and (x[2] != 0.0 or (dim < 2 and x[1] != 0.0)):
                    dim = 3
            elif cs not in warned_cs:
                warned_cs.add(cs)
                kind = "cylindrical" if ctype == 1 else "spherical"
                warn(
                    f"UNV: coordinate system {cs} is {kind}; its nodes are left in "
                    "local coordinates"
                )
        elif cs > 1 and cs not in warned_cs:
            warned_cs.add(cs)
            warn(
                f"UNV: coordinate system {cs} is not defined; its nodes are taken as global"
            )
        if label in node_index:
            warn(f"UNV: node {label} is defined twice; the last definition wins")
            coords[node_index[label]] = x
        else:
            node_index[label] = len(coords)
            coords.append(x)
    np_ = len(coords)
    points = np.array(coords, dtype=float).reshape(np_, 3)[:, :dim].copy()

    # Cells, one block per type in first-appearance order.
    blocks = {}  # type -> [conn rows, pids, mids]
    elem_ref = {}
    warned_fe = set()
    for label, fe_id, pid, mid, nodes in f.elements:
        t = _cell_type(fe_id, len(nodes))
        if t is None:
            if (fe_id, len(nodes)) not in warned_fe:
                warned_fe.add((fe_id, len(nodes)))
                warn(
                    f"UNV: FE descriptor {fe_id} with {len(nodes)} nodes is not "
                    "supported; skipping those elements"
                )
            continue
        blk = blocks.setdefault(t, [[], [], []])
        perm = _PERM.get(t)
        conn = [0] * len(nodes)
        for j, n in enumerate(nodes):
            if n not in node_index:
                raise ReadError(f"UNV: element {label} references undefined node {n}")
            conn[perm[j] if perm else j] = node_index[n]
        elem_ref[label] = (t, len(blk[1]))
        blk[0].append(conn)
        blk[1].append(pid)
        blk[2].append(mid)
    types = list(blocks)
    type_pos = {t: i for i, t in enumerate(types)}
    cells = [CellBlock(t, np.array(blocks[t][0], dtype=np.int64)) for t in types]
    block_sizes = [len(blocks[t][1]) for t in types]
    cell_data = {}
    if cells:
        cell_data["unv:pid"] = [np.array(blocks[t][1], dtype=np.int64) for t in types]
        cell_data["unv:mid"] = [np.array(blocks[t][2], dtype=np.int64) for t in types]
    mesh = Mesh(points, cells, cell_data=cell_data)

    # Groups: node members -> a point region, element members -> a cell region.
    bases = np.concatenate([[0], np.cumsum(block_sizes)]).astype(np.int64)
    for number, name, gnodes, gelems in f.groups:
        missing = 0
        pts = []
        for label in gnodes:
            if label in node_index:
                pts.append(node_index[label])
            else:
                missing += 1
        glob = []
        cdim = -1
        for label in gelems:
            ref = elem_ref.get(label)
            if ref is None:
                missing += 1
                continue
            b = type_pos[ref[0]]
            glob.append(int(bases[b]) + ref[1])
            cdim = max(cdim, _DIM.get(ref[0], 3))
        if missing:
            warn(
                f"UNV: group '{name}' names {missing} entities that were not read; "
                "they are dropped"
            )
        if gnodes:
            mesh.regions.append(
                Region(name, "point", np.array(pts, dtype=np.int64), -1, number)
            )
        if gelems or not gnodes:
            mesh.regions.append(
                Region(name, "cell", np.array(glob, dtype=np.int64), cdim, number)
            )

    if f.units is not None:
        mesh.field_data["unv:units"] = np.array([f.units[0]], dtype=np.int64)
        mesh.field_data["unv:unit_factors"] = np.array(f.units[1], dtype=float)

    if f.functions and np_ == 0:
        raise ReadError(
            "UNV: the file holds dataset-58 functions but no nodes (datasets "
            "15/781/2411) to attach them to"
        )
    groups = _group_functions(f)
    steps = _steps(f, groups)
    s = _resolve_step(time_step, len(steps))
    if s is None:
        return mesh
    mesh.time_values = [st[0][2] for st in steps]
    key, result_ids, samples = steps[s]
    if len(steps) > 1 or key[0] != 0:
        mesh.field_data[TIME_KEY] = np.array([key[2]])
        mesh.field_data["unv:analysis"] = np.array([key[0]], dtype=np.int64)
        mesh.field_data["unv:step"] = np.array([key[1]], dtype=np.int64)
    if points_only or (arrays is not None and len(arrays) == 0):
        return mesh
    wanted = None if arrays is None else set(arrays)

    def keep(name):
        return wanted is None or name in wanted

    used = {"unv:pid", "unv:mid"}

    def unique(base):
        name, k = base, 2
        while name in used:
            name = f"{base}_{k}"
            k += 1
        used.add(name)
        return name

    for r in result_ids:
        res = f.results[r]
        nc = res["ncomp"]
        base = unique(res["name"])
        parts = (
            [(base + "_real", 0), (base + "_imag", 1)]
            if res["complex"]
            else [(base, 0)]
        )
        for name, part in parts:
            if not keep(name):
                continue

            def rows(vals):
                a = np.array(vals, dtype=float)
                a = a[part::2] if res["complex"] else a
                return a[:nc]

            if res["location"] == 1:
                arr = np.full((np_, nc), _NAN)
                for label, vals in res["values"]:
                    if label in node_index:
                        arr[node_index[label]] = rows(vals)
                arr = _tensor_to_meshio(res["char"], nc, arr)
                mesh.point_data[name] = arr[:, 0].copy() if nc == 1 else arr
            elif cells:
                arrs = [np.full((ne, nc), _NAN) for ne in block_sizes]
                for label, vals in res["values"]:
                    ref = elem_ref.get(label)
                    if ref is not None:
                        arrs[type_pos[ref[0]]][ref[1]] = rows(vals)
                arrs = [_tensor_to_meshio(res["char"], nc, a) for a in arrs]
                mesh.cell_data[name] = [a[:, 0].copy() if nc == 1 else a for a in arrs]

    # Dataset-58 samples of this step, one set of arrays per function group.
    warned_nodes = set()
    for g, n in samples:
        grp = groups[g]
        fns = [f.functions[i] for i in grp["functions"]]
        axes = [abs(fn["rsp_dir"]) for fn in fns]
        has_trans = any(1 <= a <= 3 for a in axes)
        has_rot = any(4 <= a <= 6 for a in axes)
        has_scalar = any(a == 0 or a > 6 for a in axes)
        cplx = any(fn["complex"] for fn in fns)
        targets = []
        base = unique(grp["name"])
        if has_trans:
            targets.append((base, 3, 1))
        if has_rot:
            targets.append((unique(grp["name"] + "_rot"), 3, 2))
        if has_scalar:
            targets.append(
                (unique(grp["name"] + "_scalar") if has_trans else base, 1, 0)
            )
        for tname, ncomp, tkind in targets:
            parts = (
                [(tname + "_real", 0), (tname + "_imag", 1)] if cplx else [(tname, 0)]
            )
            for name, part in parts:
                if not keep(name):
                    continue
                arr = np.full((np_, ncomp), _NAN)
                for fn in fns:
                    a = abs(fn["rsp_dir"])
                    kind = 1 if 1 <= a <= 3 else (2 if 4 <= a <= 6 else 0)
                    if kind != tkind:
                        continue
                    node = fn["rsp_node"]
                    if node not in node_index:
                        if node not in warned_nodes:
                            warned_nodes.add(node)
                            warn(
                                f"UNV: dataset-58 response node {node} is not defined; "
                                "its functions are dropped"
                            )
                        continue
                    sign = (-1.0 if fn["rsp_dir"] < 0 else 1.0) * (
                        -1.0 if fn["ref_dir"] < 0 else 1.0
                    )
                    if part == 0:
                        v = fn["re"][n]
                    else:
                        v = fn["im"][n] if fn["complex"] else 0.0
                    comp = 0 if kind == 0 else (a - 1) % 3
                    arr[node_index[node], comp] = sign * v
                mesh.point_data[name] = arr[:, 0].copy() if ncomp == 1 else arr
    return mesh


def time_values(filename):
    """The value of every step, in file order (the sequence engine's steps)."""
    f = _parse(filename)
    return [st[0][2] for st in _steps(f, _group_functions(f))]


# ---------------------------------------------------------------------------
# Writing
# ---------------------------------------------------------------------------


def _field_scalar(mesh, name, default, cast):
    v = (getattr(mesh, "field_data", None) or {}).get(name)
    if v is None:
        return default
    v = np.asarray(v).ravel()
    return cast(v[0]) if v.size else default


def _write_reals(f, values, fmt, per_line):
    n = len(values)
    for i, v in enumerate(values):
        f.write(fmt % v)
        if (i + 1) % per_line == 0 or i + 1 == n:
            f.write("\n")


@_provenance.slotless_writer
def write(filename, mesh, code_aster=False, node_dataset=2411):
    points = np.asarray(mesh.points)
    if node_dataset not in (2411, 781):
        node_dataset = 2411
    np_ = len(points)
    pdim = points.shape[1] if points.ndim == 2 else 0
    fd = getattr(mesh, "field_data", None) or {}

    # LF only, matching the C++ writer's binary-mode ofstream: on Windows the
    # default text mode would translate "\n" to "\r\n", and test_unv.py's
    # test_writers_are_byte_identical pins the two engines against each other.
    with open_file(filename, "w", newline="\n") as f:
        # 164 units, when the mesh carries them.
        if "unv:units" in fd and "unv:unit_factors" in fd:
            fac = [1.0, 1.0, 1.0, 0.0]
            given = np.asarray(fd["unv:unit_factors"], dtype=float).ravel()
            for i in range(min(4, given.size)):
                fac[i] = float(given[i])
            code = _field_scalar(mesh, "unv:units", 1, int)
            f.write(f"    -1\n   164\n{code:10d}{'':20s}{2:10d}\n")
            f.write("".join(f"{v:25.16E}" for v in fac[:3]) + "\n")
            f.write(f"{fac[3]:25.16E}\n")
            f.write("    -1\n")

        # 2411 / 781 nodes
        f.write(f"    -1\n{node_dataset:6d}\n")
        for k in range(np_):
            f.write(f"{k + 1:10d}{1:10d}{1:10d}{11:10d}\n")
            row = [float(points[k, c]) if c < pdim else 0.0 for c in range(3)]
            f.write("".join(f"{x:25.16E}" for x in row) + "\n")
        f.write("    -1\n")

        # 2412 elements
        f.write("    -1\n  2412\n")
        cd = getattr(mesh, "cell_data", None) or {}
        pid_blocks = cd.get("unv:pid")
        mid_blocks = cd.get("unv:mid")
        label = 0
        elem_labels = []
        for bi, cell_block in enumerate(mesh.cells):
            t = cell_block.type
            if t not in _MESHIO_TO_UNV:
                warn(f"UNV does not support '{t}' cells. Skipping.")
                _provenance.note(
                    "cells-dropped", f"cell block(s) of type {t} have no UNV equivalent"
                )
                elem_labels.append(None)
                continue
            descriptor, is_beam = _MESHIO_TO_UNV[t]
            perm = _PERM.get(t)
            data = np.asarray(cell_block.data)
            nnodes = data.shape[1]
            pids = (
                np.asarray(pid_blocks[bi]).astype(np.int64)
                if pid_blocks is not None and bi < len(pid_blocks)
                else None
            )
            mids = (
                np.asarray(mid_blocks[bi]).astype(np.int64)
                if mid_blocks is not None and bi < len(mid_blocks)
                else None
            )
            labels_this = []
            for li, row in enumerate(data):
                label += 1
                labels_this.append(label)
                pid = int(pids[li]) if pids is not None else 1
                mid = int(mids[li]) if mids is not None else pid
                f.write(
                    f"{label:10d}{descriptor:10d}{pid:10d}{mid:10d}{11:10d}{nnodes:10d}\n"
                )
                if is_beam:
                    f.write("         0         1         1\n")
                ints = [int(row[perm[j] if perm else j]) + 1 for j in range(nnodes)]
                for i in range(0, nnodes, 8):
                    f.write("".join(f"{v:10d}" for v in ints[i : i + 8]) + "\n")
            elem_labels.append(labels_this)
        f.write("    -1\n")

        _write_groups(f, mesh, elem_labels)
        _write_fields(f, mesh, elem_labels, code_aster)


def _write_groups(f, mesh, elem_labels):
    """2467 groups from the mesh's regions; a point and a cell region sharing a name
    are one group, side regions are dropped."""
    groups = {}
    order = []
    bases = [0]
    for cb in mesh.cells:
        bases.append(bases[-1] + len(cb.data))

    def group(name):
        if name not in groups:
            groups[name] = {"number": -1, "nodes": [], "elements": []}
            order.append(name)
        return groups[name]

    side_dropped = 0
    cells_dropped = 0
    for reg in sorted(
        getattr(mesh, "regions", []) or [],
        key=lambda r: (
            {"point": 0, "cell": 1, "side": 2}[r.kind],
            r.name,
            r.dim,
            r.tag,
        ),
    ):
        if reg.kind == "side":
            side_dropped += 1
            continue
        g = group(reg.name)
        if g["number"] < 0 and reg.tag > 0:
            g["number"] = int(reg.tag)
        entries = np.unique(np.asarray(reg.entries, dtype=np.int64).ravel())
        if reg.kind == "point":
            g["nodes"] += [int(e) + 1 for e in entries]
        else:
            for e in entries:
                b = int(np.searchsorted(bases, e, side="right")) - 1
                if b < 0 or b >= len(mesh.cells) or elem_labels[b] is None:
                    cells_dropped += 1
                    continue
                g["elements"].append(elem_labels[b][int(e) - bases[b]])
    if side_dropped:
        warn(f"UNV has no facet groups; {side_dropped} side region(s) dropped")
        _provenance.note(
            "regions-dropped", f"{side_dropped} side region(s) have no UNV group"
        )
    if cells_dropped:
        warn(
            f"UNV: {cells_dropped} group member(s) sit in cell blocks UNV cannot hold; dropped"
        )
    if not order:
        return
    taken = set()
    for name in order:
        g = groups[name]
        if g["number"] > 0:
            if g["number"] in taken:
                g["number"] = -1
            else:
                taken.add(g["number"])
    nxt = 1
    for name in order:
        g = groups[name]
        if g["number"] > 0:
            continue
        while nxt in taken:
            nxt += 1
        g["number"] = nxt
        taken.add(nxt)
    f.write("    -1\n  2467\n")
    for name in sorted(order, key=lambda n: groups[n]["number"]):
        g = groups[name]
        ents = [(7, t) for t in g["nodes"]] + [(8, t) for t in g["elements"]]
        f.write(f"{g['number']:10d}" + f"{0:10d}" * 6 + f"{len(ents):10d}\n")
        f.write(f"{name}\n")
        for i in range(0, len(ents), 2):
            f.write(
                "".join(
                    f"{et:10d}{tag:10d}{0:10d}{0:10d}" for et, tag in ents[i : i + 2]
                )
                + "\n"
            )
    f.write("    -1\n")


def _as_2d(arr):
    arr = np.asarray(arr, dtype=float)
    return arr.reshape(len(arr), -1)


def _write_fields(f, mesh, elem_labels, code_aster):
    analysis = _field_scalar(mesh, "unv:analysis", 0, int)
    step_id = _field_scalar(mesh, "unv:step", 0, int)
    time = _field_scalar(mesh, TIME_KEY, 0.0, float)

    def write_field(field_id, name, location, labels, data):
        keep = ~np.all(np.isnan(data), axis=1)
        labels = [lb for lb, k in zip(labels, keep) if k]
        data = data[keep]
        nc = data.shape[1]
        char = _data_char(nc)
        if char in (4, 5):
            data = _tensor_from_meshio(nc, data)
        if code_aster:
            f.write(f"    -1\n{55 if location == 1 else 56:6d}\n")
            f.write(f"{name}\n" + "NONE\n" * 4)
            f.write(f"{1:10d}{analysis:10d}{char:10d}{0:10d}{2:10d}{nc:10d}\n")
            if analysis == 2:
                ints, reals = [1, step_id], [time, 0.0, 0.0, 0.0]
            elif analysis in (3, 7):
                ints, reals = [1, step_id], [
                    0.0,
                    2.0 * math.pi * time,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ]
            elif analysis in (4, 5):
                ints, reals = [1, step_id], [time]
            else:
                ints, reals = [step_id], [time]
            f.write(
                f"{len(ints):10d}{len(reals):10d}" + "".join(f"{v:10d}" for v in ints)
            )
            f.write("\n")
            _write_reals(f, reals, "%13.5E", 6)
            fmt, per = "%13.5E", 6
        else:
            f.write("    -1\n  2414\n")
            f.write(f"{field_id:10d}\n{name}\n{location:10d}\n")
            f.write("NONE\n" * 5)
            f.write(f"{1:10d}{analysis:10d}{char:10d}{0:10d}{4:10d}{nc:10d}\n")
            r10 = [0] * 8
            r12 = [0.0] * 12
            if analysis == 2:
                r10[5], r12[1] = step_id, time
            elif analysis in (3, 7):
                r10[5], r12[7] = step_id, 2.0 * math.pi * time
            elif analysis in (4, 9):
                r10[6], r12[0] = step_id, time
            elif analysis == 5:
                r10[7], r12[1] = step_id, time
            elif analysis == 6:
                r10[5], r12[2] = step_id, time
            else:
                r10[4], r12[0] = step_id, time
            f.write("".join(f"{v:10d}" for v in r10) + "\n")
            f.write("         0         0\n")
            _write_reals(f, r12[:6], "%13.5E", 6)
            _write_reals(f, r12[6:], "%13.5E", 6)
            # One line per entity: pyuff pairs each element record with one value line.
            fmt, per = "%20.12E", max(nc, 1)
        for lb, row in zip(labels, data):
            if location == 1:
                f.write(f"{int(lb):10d}\n")
            else:
                f.write(f"{int(lb):10d}{nc:10d}\n")
            _write_reals(f, [float(v) for v in row], fmt, per)
        f.write("    -1\n")

    field_id = 0
    for name in sorted(getattr(mesh, "point_data", None) or {}):
        data = _as_2d(mesh.point_data[name])
        if data.size == 0 or data.shape[1] == 0:
            continue
        field_id += 1
        write_field(field_id, name, 1, range(1, len(data) + 1), data)

    for name in sorted(getattr(mesh, "cell_data", None) or {}):
        if name in ("unv:pid", "unv:mid"):
            continue
        blocks = mesh.cell_data[name]
        labels, rows = [], []
        nc = None
        for bi, blk in enumerate(blocks):
            if bi >= len(elem_labels) or elem_labels[bi] is None or blk is None:
                continue
            blk = _as_2d(blk)
            if blk.size == 0:
                continue
            if nc is None:
                nc = blk.shape[1]
            if blk.shape[1] != nc:
                continue
            labels += elem_labels[bi]
            rows.append(blk)
        if nc is None or not labels:
            continue
        field_id += 1
        write_field(field_id, name, 2, labels, np.concatenate(rows))

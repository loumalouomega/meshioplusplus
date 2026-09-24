"""I/O for Abaqus results files (``.fil``), binary and ASCII.

The pure-Python twin of ``src/cpp/src/formats/abaqus_fil.cpp``: both engines read
the same meshes.

A ``.fil`` file is a sequence of records ``[length, key, attributes...]`` of
8-byte words: in 512-word Fortran-record blocks (binary), or ``*``-started
records of ``I``/``D``/``A`` items on 80-column lines (ASCII). The model records
(1901 nodes, 1900/1990 elements, 1931-1934 sets, 1940 labels) make the mesh;
every increment (2000 ... 2001) is a step whose nodal and element records become
point and cell data named by their Abaqus identifier. An eigenvalue step's
modes (1980) are steps of their own; contact surfaces (1501/1502) are side
regions and contact output (1503-15xx) point data; energies (1999), modal
dynamics quantities (301-310) and element matrices (1001-1031) are field data.
"""

import struct

import numpy as np

from .._common import num_nodes_per_cell, warn
from .._exceptions import ReadError
from .._files import open_file
from .._fortran_records import fortran_records, sniff_fortran_records
from .._mesh import Mesh, topological_dimension
from .._regions import Region
from ..abaqus._abaqus import _face_index, abaqus_to_meshio_type

__all__ = ["read", "time_values"]

_NODAL_NAMES = {
    101: "U",
    102: "V",
    103: "A",
    104: "RF",
    105: "EPOT",
    106: "CF",
    107: "COORD",
    108: "POR",
    109: "RVF",
    110: "RVT",
    113: "TU",
    114: "TV",
    115: "TA",
    119: "RCHG",
    120: "CECHG",
    136: "PCAV",
    137: "CVOL",
    145: "VF",
    146: "TF",
    151: "PABS",
    201: "NT",
    204: "RFL",
    206: "CFL",
    214: "RFLE",
    221: "NNC",
    320: "CFF",
}

_ELEMENT_NAMES = {
    2: "TEMP",
    3: "LOADS",
    4: "FLUXS",
    5: "SDV",
    6: "VOIDR",
    7: "FOUND",
    8: "COORD",
    9: "FV",
    10: "NFLUX",
    11: "S",
    12: "SINV",
    13: "SF",
    14: "ENER",
    15: "NFORC",
    17: "JK",
    18: "POR",
    19: "ELEN",
    21: "E",
    22: "PE",
    23: "CE",
    24: "IE",
    25: "EE",
    26: "CRACK",
    27: "STH",
    28: "HFL",
    29: "SE",
    30: "DG",
    31: "CONF",
    32: "SJP",
    35: "SAT",
    36: "SS",
    38: "CONC",
    39: "MFL",
    42: "SPE",
    45: "PEQC",
    47: "SEPE",
    48: "TSHR",
    50: "EPG",
    51: "EFLX",
    61: "STATUS",
    73: "PEEQ",
    74: "PRESS",
    75: "MISES",
    76: "IVOL",
    77: "SVOL",
    78: "EVOL",
    83: "SSAVG",
    86: "ALPHA",
    87: "UVARM",
    88: "THE",
    89: "LE",
    90: "NE",
    91: "ER",
    401: "SP",
    402: "ALPHAP",
    403: "EP",
    404: "NEP",
    405: "LEP",
    406: "ERP",
    407: "DGP",
    408: "EEP",
    409: "IEP",
    410: "THEP",
    411: "PEP",
    412: "CEP",
    # Abaqus/Explicit only
    421: "CKE",
    422: "CKLE",
    423: "CKLS",
    424: "CKSTAT",
    441: "CKEMAG",
    476: "EMSF",
    477: "EDT",
    507: "CFAILST",
    559: "CDMG",
    560: "CDIF",
    561: "CDIM",
    562: "CDIP",
}
# Records 301-310, per increment of a mode-based dynamic step.
_MODAL_NAMES = dict(
    zip(
        range(301, 311), ("GU", "GV", "GA", "BM", "GPU", "GPV", "GPA", "SNE", "KE", "T")
    )
)
# Record 1999's attributes: Standard's, and Explicit's where they differ.
_ENERGY_NAMES = (
    "ALLKE",
    "ALLSE",
    "ALLWK",
    "ALLPD",
    "ALLCD",
    "ALLVD",
    "ALLKL",
    "ALLAE",
    "ALLQB",
    "ALLEE",
    "ALLIE",
    "ETOTAL",
    "ALLFD",
    "ALLJD",
    "ALLSD",
    "ALLDMD",
)
_EXPLICIT_ENERGY_NAMES = (
    "ALLKE",
    "ALLSE",
    "ALLWK",
    "ALLPD",
    "ALLCD",
    "ALLVD",
    None,
    "ALLAE",
    "ALLDC",
    None,
    "ALLIE",
    "ETOTAL",
    "ALLFD",
    None,
    "DMASS",
    "ALLDMD",
    "ALLIHE",
    "ALLHF",
)
# Contact output (records 1511-1592, after a 1504 node header).
_CONTACT_NAMES = {
    1511: "CSTRESS",
    1512: "CDSTRESS",
    1521: "CDISP",
    1522: "CFN",
    1523: "CFS",
    1524: "CAREA",
    1526: "CMN",
    1527: "CMS",
    1528: "HFL",
    1529: "HFLA",
    1530: "HTL",
    1531: "HTLA",
    1532: "SFDR",
    1533: "SFDRA",
    1534: "SFDRT",
    1535: "SFDRTA",
    1536: "WEIGHT",
    1537: "SJD",
    1538: "SJDA",
    1539: "SJDT",
    1540: "SJDTA",
    1541: "ECD",
    1542: "ECDA",
    1543: "ECDT",
    1544: "ECDTA",
    1545: "PFL",
    1546: "PFLA",
    1547: "PTL",
    1548: "PTLA",
    1549: "TPFL",
    1550: "TPTL",
    1570: "DBT",
    1571: "DBSF",
    1572: "DBS",
    1573: "XN",
    1574: "XS",
    1575: "CFT",
    1576: "CMT",
    1577: "XT",
    1578: "CTRQ",
    1592: "PPRESS",
}
# Element matrix records: the field data stem of each.
_MATRIX_NAMES = {
    1002: "matrix_dofs",
    1011: "stiffness",
    1012: "stiffness",
    1021: "mass",
    1022: "mass",
    1031: "load",
}
_EXPLICIT_PROCEDURES = (17, 21, 74)

_SOLID = ("C3D", "DC3D", "AC3D", "COH3D", "SC", "CCL")
_PLANAR = (
    "CPE",
    "CPS",
    "CAX",
    "CGAX",
    "DC2D",
    "DCAX",
    "COH2D",
    "COHAX",
    "M3D",
    "R3D",
    "SFM3D",
    "CPEG",
)
_LINE = (
    "T2D",
    "T3D",
    "B2",
    "B3",
    "PIPE",
    "R2D",
    "RB2D",
    "RB3D",
    "DC1D",
    "SAX",
    "FRAME",
)


def abaqus_cell_type(name, count):
    """The meshio++ type of an Abaqus element with ``count`` nodes, or ``None``.

    The twin of ``detail::abaqus_cell_type``: the ``.inp`` table first, then the
    element family by prefix and the shape by node count.
    """
    t = abaqus_to_meshio_type.get(name)
    if t is not None and num_nodes_per_cell.get(t) == count:
        return t
    if name.startswith(_SOLID):
        return {
            4: "tetra",
            5: "pyramid",
            6: "wedge",
            8: "hexahedron",
            10: "tetra10",
            13: "pyramid13",
            15: "wedge15",
            20: "hexahedron20",
        }.get(count)
    shell = (len(name) > 1 and name[0] == "S" and name[1].isdigit()) or name.startswith(
        "STRI"
    )
    if name.startswith(_PLANAR) or shell:
        return {3: "triangle", 4: "quad", 6: "triangle6", 8: "quad8", 9: "quad9"}.get(
            count
        )
    if name.startswith(_LINE):
        return {2: "line", 3: "line3"}.get(count)
    return None


class _Word:
    """One 8-byte word: an ASCII item (tag I/D/A) or a raw binary word (tag B)."""

    __slots__ = ("tag", "value")

    def __init__(self, tag, value):
        self.tag = tag
        self.value = value


class _Data:
    def __init__(self):
        self.records = []  # (key, [words])
        self.order = "<"  # binary byte order

    def int(self, w):
        if w.tag == "I":
            return w.value
        if w.tag == "D":
            return int(w.value)
        if w.tag == "A":
            s = w.value.strip(" \0")
            return int(s) if s.lstrip("-").isdigit() else 0
        wide = struct.unpack(self.order + "q", w.value)[0]
        if -2147483648 <= wide <= 2147483647:
            return wide
        return struct.unpack(self.order + "i", w.value[:4])[0]

    def real(self, w):
        if w.tag in ("I", "D"):
            return float(w.value)
        if w.tag == "A":
            return float("nan")
        return struct.unpack(self.order + "d", w.value)[0]

    def text(self, w):
        if w.tag in ("I", "D"):
            return str(self.int(w))
        if w.tag == "A":
            return w.value
        return w.value.decode("latin-1")


def _parse_binary(data, out):
    layout = sniff_fortran_records(data)
    if layout is None:
        words = data
    else:
        words = b"".join(
            data[o : o + n] for o, n in fortran_records(data, layout, "Abaqus .fil")
        )
        out.order = layout[1]
    if len(words) % 8:
        raise ReadError(
            "Abaqus .fil: binary payload is not a whole number of 8-byte words"
        )
    n = len(words) // 8
    w = 0
    while w < n:
        length = out.int(_Word("B", words[8 * w : 8 * w + 8]))
        if length == 0:
            w += 1
            continue
        if length < 2 or length > n - w:
            raise ReadError(
                f"Abaqus .fil: record at word {w} has invalid length {length}"
            )
        key = out.int(_Word("B", words[8 * w + 8 : 8 * w + 16]))
        attrs = [_Word("B", words[8 * k : 8 * k + 8]) for k in range(w + 2, w + length)]
        out.records.append((key, attrs))
        w += length


def _parse_ascii(text, out):
    s = text.replace("\n", "").replace("\r", "")
    pos = 0

    def fail(what):
        raise ReadError(f"Abaqus .fil: {what} (character {pos})")

    def item():
        nonlocal pos
        if pos >= len(s):
            fail("the file ends inside a record")
        tag = s[pos]
        if tag == "I":
            width = s[pos + 1 : pos + 3].strip(" ")
            if not width.isdigit():
                fail("bad integer width")
            n = int(width)
            digits = s[pos + 3 : pos + 3 + n]
            if n < 1 or len(digits) != n:
                fail("bad integer width")
            try:
                value = int(digits)
            except ValueError:
                fail(f"bad integer '{digits}'")
            pos += 3 + n
            return _Word("I", value)
        if tag in ("D", "E"):
            num = s[pos + 1 : pos + 23]
            if len(num) != 22:
                fail("truncated real item")
            num = num.replace("D", "E").replace("d", "E")
            if "E" not in num:
                k = max(num.rfind("+"), num.rfind("-"))
                if k > 2:
                    num = num[:k] + "E" + num[k:]
            try:
                value = float(num)
            except ValueError:
                fail(f"bad real '{num}'")
            pos += 23
            return _Word("D", value)
        if tag == "A":
            if pos + 9 > len(s):
                fail("truncated text item")
            value = s[pos + 1 : pos + 9]
            pos += 9
            return _Word("A", value)
        fail(f"unknown item tag '{tag}'")

    while True:
        pos = s.find("*", pos)
        if pos < 0:
            break
        pos += 1
        length_word = item()
        if length_word.tag != "I":
            fail("a record starts with its length")
        length = length_word.value
        if length < 2:
            fail(f"record length {length}")
        key = out.int(item())
        out.records.append((key, [item() for _ in range(length - 2)]))


def _parse(filename):
    with open_file(filename, "rb") as f:
        data = f.read()
    if isinstance(data, str):
        data = data.encode("latin-1")
    out = _Data()
    stripped = data.lstrip(b" \t\r\n")
    if stripped[:1] == b"*":
        _parse_ascii(data.decode("latin-1"), out)
    else:
        _parse_binary(data, out)
    return out


def _increments(d):
    out = []
    for r, (key, w) in enumerate(d.records):
        if key == 2000:
            inc = {
                "begin": r + 1,
                "end": None,
                "total": d.real(w[0]) if len(w) > 0 else 0.0,
                "step_time": d.real(w[1]) if len(w) > 1 else 0.0,
                "procedure": d.int(w[4]) if len(w) > 4 else 0,
                "step": d.int(w[5]) if len(w) > 5 else 0,
                "increment": d.int(w[6]) if len(w) > 6 else 0,
            }
            if out and out[-1]["end"] is None:
                out[-1]["end"] = r
            out.append(inc)
        elif key == 2001 and out and out[-1]["end"] is None:
            out[-1]["end"] = r
    if out and out[-1]["end"] is None:
        out[-1]["end"] = len(d.records)
    # An eigenvalue step writes one increment whose modes each start with a
    # 1980 record: every mode becomes an increment of its own.
    split = []
    for inc in out:
        modes = [r for r in range(inc["begin"], inc["end"]) if d.records[r][0] == 1980]
        if not modes:
            split.append(inc)
            continue
        for k, r in enumerate(modes):
            split.append(
                dict(
                    inc, begin=r, end=modes[k + 1] if k + 1 < len(modes) else inc["end"]
                )
            )
    return split


def time_values(filename):
    return [inc["total"] for inc in _increments(_parse(filename))]


def _resolve_step(time_step, n):
    k = time_step + n if time_step < 0 else time_step
    if n == 0:
        if time_step not in (0, -1):
            raise ReadError(
                f"time_step {time_step} requested but the file has no steps"
            )
        return None
    if not 0 <= k < n:
        raise ReadError(f"time_step {time_step} out of range: the file has {n} step(s)")
    return k


def _skipped_key(key):
    return 301 <= key <= 310 or key >= 1000


def read(filename, points_only=False, arrays=None, time_step=0):
    d = _parse(filename)
    node_labels, coords = [], []
    elements = []
    sets = []
    labels = {}
    surfaces = []  # (name, type, [(element, face key)])
    matrices = {}  # kind -> (values, index rows)
    matrix_element = 0
    last_matrix_key = 0
    for key, w in d.records:
        if key != last_matrix_key and not 1001 <= key <= 1043:
            last_matrix_key = 0
        if key == 1501 and w:
            surfaces.append(
                (d.text(w[0]).strip(" \0"), d.int(w[2]) if len(w) > 2 else 1, [])
            )
            continue
        if key == 1502:
            if surfaces and len(w) >= 2:
                surfaces[-1][2].append((d.int(w[0]), d.int(w[1])))
            continue
        if key == 1001:
            matrix_element = d.int(w[0]) if w else 0
            last_matrix_key = 0
            continue
        if key in _MATRIX_NAMES:
            values, index = matrices.setdefault(_MATRIX_NAMES[key], ([], []))
            flag = 1 if key in (1011, 1021) else 0
            first = 0
            if key == 1031:
                flag = d.int(w[0]) if w else 0
                first = 1
            continued = key == last_matrix_key and index
            if not continued:
                index.append([matrix_element, len(values), 0, flag])
            for x in w[0 if continued else first :]:
                values.append(float(d.int(x)) if key == 1002 else d.real(x))
            index[-1][2] = len(values) - index[-1][1]
            last_matrix_key = key
            continue
        if key == 1901 and w:
            node_labels.append(d.int(w[0]))
            coords.append(
                [d.real(w[k + 1]) if k + 1 < len(w) else 0.0 for k in range(3)]
            )
        elif key == 1900 and len(w) >= 2:
            elements.append(
                (
                    d.int(w[0]),
                    d.text(w[1]).strip(" \0").upper(),
                    [d.int(x) for x in w[2:]],
                )
            )
        elif key == 1990 and elements:
            elements[-1][2].extend(d.int(x) for x in w)
        elif key in (1931, 1933) and w:
            sets.append(
                (d.text(w[0]).strip(" \0"), key == 1931, [d.int(x) for x in w[1:]])
            )
        elif key in (1932, 1934) and sets:
            sets[-1][2].extend(d.int(x) for x in w)
        elif key == 1940 and w:
            labels[d.int(w[0])] = "".join(d.text(x) for x in w[1:]).strip(" \0")

    node_index = {}
    for p, label in enumerate(node_labels):
        if label in node_index:
            raise ReadError(f"Abaqus .fil: node {label} is defined twice")
        node_index[label] = p
    points = np.array(coords, dtype=np.float64).reshape(-1, 3)

    order_of_types = []
    by_type = {}
    skipped_types = {}
    for e, (_, etype, nodes) in enumerate(elements):
        t = abaqus_cell_type(etype, len(nodes))
        if t is None:
            skipped_types[etype] = skipped_types.get(etype, 0) + 1
            continue
        if t not in by_type:
            by_type[t] = []
            order_of_types.append(t)
        by_type[t].append(e)
    for etype in sorted(skipped_types):
        warn(
            f"Abaqus .fil: skipping {skipped_types[etype]} element(s) of type {etype} "
            "(no meshio++ equivalent)"
        )
    cells = []
    id_blocks = []
    element_index = {}
    cell_nodes = []
    cell_dim = []
    block_start = [0]
    for t in order_of_types:
        members = by_type[t]
        k = num_nodes_per_cell[t]
        conn = np.empty((len(members), k), dtype=np.int64)
        ids = []
        for r, e in enumerate(members):
            label, _, nodes = elements[e]
            for j in range(k):
                idx = node_index.get(nodes[j])
                if idx is None:
                    raise ReadError(
                        f"Abaqus .fil: element {label} names undefined node {nodes[j]}"
                    )
                conn[r, j] = idx
            if label in element_index:
                raise ReadError(f"Abaqus .fil: element {label} is defined twice")
            element_index[label] = len(cell_nodes)
            cell_nodes.append(nodes)
            cell_dim.append(topological_dimension[t])
            ids.append(label)
        cells.append((t, conn))
        id_blocks.append(np.array(ids, dtype=np.int64))
        block_start.append(len(cell_nodes))

    mesh = Mesh(points, cells)
    mesh.point_data["abaqus:id"] = np.array(node_labels, dtype=np.int64)
    if id_blocks:
        mesh.cell_data["abaqus:id"] = id_blocks

    regions = []
    missing = 0
    for name, is_nodes, members in sets:
        if name.isdigit() and int(name) in labels:
            name = labels[int(name)]
        entries = []
        dim = -1
        for label in members:
            idx = (node_index if is_nodes else element_index).get(label)
            if idx is None:
                missing += 1
                continue
            entries.append(idx)
            if not is_nodes:
                dim = max(dim, cell_dim[idx])
        regions.append(
            Region(
                name,
                "point" if is_nodes else "cell",
                np.array(entries, dtype=np.int64),
                -1 if is_nodes else dim,
                -1,
            )
        )
    if missing:
        warn(
            f"Abaqus .fil: sets name {missing} node(s) or element(s) that are not in "
            "the mesh"
        )
    unplaced = 0
    for name, _, facets in surfaces:
        if name.isdigit() and int(name) in labels:
            name = labels[int(name)]
        entries = []
        dim = -1
        for label, face in facets:
            c = element_index.get(label)
            if c is None or not 1 <= face <= 8:
                unplaced += 1
                continue
            b = int(np.searchsorted(block_start, c, side="right")) - 1
            key = "SPOS" if face == 7 else ("SNEG" if face == 8 else f"S{face}")
            facet = _face_index(order_of_types[b], key)
            if facet is None:
                unplaced += 1
                continue
            entries.append((c, facet))
            dim = max(dim, cell_dim[c] - 1)
        regions.append(
            Region(
                name, "side", np.array(entries, dtype=np.int64).reshape(-1, 2), dim, -1
            )
        )
    if unplaced:
        warn(
            f"Abaqus .fil: {unplaced} contact surface facet(s) name no element face of "
            "the mesh"
        )
    mesh.regions = regions
    for kind in sorted(matrices):
        values, index = matrices[kind]
        dtype = np.int64 if kind == "matrix_dofs" else np.float64
        mesh.field_data["abaqus:" + kind] = np.array(values, dtype=dtype)
        mesh.field_data["abaqus:" + kind + ":index"] = np.array(
            index, dtype=np.int64
        ).reshape(-1, 4)

    increments = _increments(d)
    mesh.time_values = [inc["total"] for inc in increments]
    if not increments:
        _resolve_step(time_step, 0)
        return mesh
    inc = increments[_resolve_step(time_step, len(increments))]
    mesh.field_data["meshio:time"] = np.array(inc["total"], dtype=np.float64)
    mesh.field_data["abaqus:step"] = np.array(inc["step"], dtype=np.int64)
    mesh.field_data["abaqus:increment"] = np.array(inc["increment"], dtype=np.int64)
    mesh.field_data["abaqus:step_time"] = np.array(inc["step_time"], dtype=np.float64)
    mesh.field_data["abaqus:procedure"] = np.array(inc["procedure"], dtype=np.int64)
    if points_only or (arrays is not None and len(arrays) == 0):
        return mesh

    def wants(name):
        return arrays is None or name in arrays

    nodal = {}
    fields = {}
    skipped_keys = set()
    nodal_mode = False
    header = None
    rebar = ""
    explicit = inc["procedure"] in _EXPLICIT_PROCEDURES
    modal_rows = {}
    contact = False
    contact_node = 0
    for key, w in d.records[inc["begin"] : inc["end"]]:
        if key == 1980:  # a mode: number, eigenvalue, mass, damping, factors
            fd = mesh.field_data
            if len(w) > 0:
                fd["abaqus:mode"] = np.array(d.int(w[0]), dtype=np.int64)
            if len(w) > 1:
                fd["abaqus:eigenvalue"] = np.array(d.real(w[1]))
            if len(w) > 2:
                fd["abaqus:generalized_mass"] = np.array(d.real(w[2]))
            if len(w) > 3:
                fd["abaqus:composite_damping"] = np.array(d.real(w[3]))
            pairs = [(d.real(w[k]), d.real(w[k + 1])) for k in range(4, len(w) - 1, 2)]
            if pairs:
                fd["abaqus:participation_factor"] = np.array([p[0] for p in pairs])
                fd["abaqus:effective_mass"] = np.array([p[1] for p in pairs])
            continue
        if key == 1999:  # total energies
            names = _EXPLICIT_ENERGY_NAMES if explicit else _ENERGY_NAMES
            for k, x in enumerate(w):
                if k < len(names) and names[k]:
                    mesh.field_data["abaqus:" + names[k]] = np.array(d.real(x))
            continue
        if key in _MODAL_NAMES:  # generalized quantities per mode
            row = [
                float(d.int(x)) if key == 304 and k == 0 else d.real(x)
                for k, x in enumerate(w)
                if not (key == 304 and k == 7)  # BM's base name
            ]
            modal_rows.setdefault(_MODAL_NAMES[key], []).append(row)
            continue
        if key == 1503:  # contact output request: the node records follow
            contact = True
            continue
        if key == 1504:
            contact_node = d.int(w[0]) if w else 0
            continue
        if contact and 1505 <= key <= 1599:
            name = _CONTACT_NAMES.get(key, f"key_{key}")
            if wants(name):
                nodal.setdefault(name, {})[contact_node] = [d.real(x) for x in w]
            continue
        if key == 1911:
            nodal_mode = bool(w) and d.int(w[0]) == 1
            header = None
            contact = False
            continue
        if key == 1:
            header = (
                (d.int(w[0]), d.int(w[1]), d.int(w[2]), d.int(w[3]))
                if len(w) >= 4
                else None
            )
            rebar = ""
            if header is not None and header[3] == 3 and len(w) > 4:
                rebar = d.text(w[4]).strip(" \0")
                if rebar.isdigit() and int(rebar) in labels:
                    rebar = labels[int(rebar)]
            continue
        if 1900 <= key <= 2001 or _skipped_key(key):
            # model records, and those read from the whole file (surfaces,
            # element matrices) are not skipped results
            whole_file = 1001 <= key <= 1043 or key in (1501, 1502)
            if not 1900 <= key <= 2001 and not whole_file:
                skipped_keys.add(key)
            continue
        if nodal_mode:
            if not w:
                continue
            name = _NODAL_NAMES.get(key, f"key_{key}")
            if not wants(name):
                continue
            nodal.setdefault(name, {})[d.int(w[0])] = [d.real(x) for x in w[1:]]
            continue
        if header is None:
            skipped_keys.add(key)
            continue
        elem, point, section, location = header
        if key == 79:
            name = "ERV" if explicit else "RATIO"
        else:
            name = _ELEMENT_NAMES.get(key, f"key_{key}")
        # Rebar values are per integration point of the element, named by the rebar.
        if location == 3:
            name += f"@rebar:{rebar}"
            location = 0
        elif location == 5:  # the whole element: one value set, like a centroid
            location = 1
        # Continuum elements write section point 0, shells and beams 1..n: the
        # first one shares the plain name, the others are "@sp<k>".
        if section > 1:
            name += f"@sp{section}"
        if not wants(name):
            continue
        fkey = (name, location)
        if fkey not in fields:
            fields[fkey] = {}
        fields[fkey].setdefault(elem, {})[point] = [d.real(x) for x in w]
    if skipped_keys:
        warn(
            f"Abaqus .fil: {len(skipped_keys)} record key(s) outside the nodal and "
            f"element results skipped (first: {min(skipped_keys)})"
        )
    for g, rows in modal_rows.items():  # one row per record, NaN-padded
        width = max(len(r) for r in rows)
        a = np.full((len(rows), width), np.nan)
        for k, r in enumerate(rows):
            a[k, : len(r)] = r
        mesh.field_data["abaqus:" + g] = a[0] if len(rows) == 1 else a

    npts = len(node_labels)

    def add_point_field(name, values):
        width = max((len(v) for v in values.values()), default=0)
        if not width:
            return
        a = np.full((npts, width), np.nan)
        for label, v in values.items():
            idx = node_index.get(label)
            if idx is not None:
                a[idx, : len(v)] = v
        while name in mesh.point_data:
            name += "@avg"
        mesh.point_data[name] = a[:, 0].copy() if width == 1 else a

    for name, values in nodal.items():
        add_point_field(name, values)

    nblocks = len(order_of_types)
    for (name, location), per_elem in fields.items():
        if location == 4:
            add_point_field(
                name,
                {
                    label: next(iter(points_.values()))
                    for label, points_ in per_elem.items()
                    if points_
                },
            )
            continue
        if not nblocks:
            continue
        # One rectangular array with the same width in every block (what every
        # writer can hold): per-point data is flattened point-major, column
        # point * W + component, NaN where a block has fewer points or components.
        per_point = location in (0, 2)
        width = 0
        count = 1
        if location == 2:
            count = max(num_nodes_per_cell[t] for t in order_of_types)
        rows = [[None] * (block_start[b + 1] - block_start[b]) for b in range(nblocks)]
        for label, points_ in per_elem.items():
            c = element_index.get(label)
            if c is None:
                continue
            b = int(np.searchsorted(block_start, c, side="right")) - 1
            rows[b][c - block_start[b]] = points_
            for pt, v in points_.items():
                width = max(width, len(v))
                if location == 0:
                    count = max(count, max(pt, 1))
        if not width:
            continue
        pts = count if per_point else 1
        blocks = []
        for b in range(nblocks):
            n = len(rows[b])
            a = np.full((n, pts, width), np.nan)
            for r in range(n):
                if rows[b][r] is None:
                    continue
                cn = cell_nodes[block_start[b] + r]
                for pt, v in rows[b][r].items():
                    if location == 0:
                        slot = max(pt, 1) - 1
                    elif location == 2:
                        if pt not in cn:
                            continue
                        slot = cn.index(pt)
                    else:
                        slot = 0
                    if slot >= pts:
                        continue
                    a[r, slot, : min(len(v), width)] = v[:width]
            a = a.reshape(n, pts * width)
            if pts * width == 1:
                a = a[:, 0]
            blocks.append(np.ascontiguousarray(a))
        while name in mesh.cell_data:
            name += "@ip" if location == 0 else "@el"
        mesh.cell_data[name] = blocks
        if per_point:
            # How to unflatten it: (points, components).
            mesh.field_data["abaqus:layout:" + name] = np.array(
                [pts, width], dtype=np.int64
            )
    return mesh

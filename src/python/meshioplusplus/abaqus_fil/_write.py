"""Abaqus results files (``.fil``), written.

The twin of the writing half of ``src/cpp/src/formats/abaqus_fil.cpp``: both
engines write the same bytes. The writer is the reader's inverse for the model
(nodes, elements, node and element sets) and the nodal and element results of
every increment, in the ASCII layout (``*FILE FORMAT, ASCII``) or the binary one
(512-word blocks of 8-byte words, little-endian).

A file is a list of records ``[key, items...]``: the model, then per increment a
2000 record, its results and a 2001 record. ``_ascii`` and ``_binary`` encode
the same records.
"""

import re
import struct

import numpy as np

from .. import _provenance
from .._common import warn
from .._exceptions import WriteError
from .._files import is_buffer
from ..abaqus._abaqus import _face_name
from ._abaqus_fil import (
    _CONTACT_NAMES,
    _ELEMENT_NAMES,
    _ENERGY_NAMES,
    _EXPLICIT_ENERGY_NAMES,
    _EXPLICIT_PROCEDURES,
    _NODAL_NAMES,
)

TIME_KEY = "meshio:time"

_NODAL_KEYS = {name: key for key, name in _NODAL_NAMES.items()}
_ELEMENT_KEYS = {name: key for key, name in _ELEMENT_NAMES.items()}
_ELEMENT_KEYS["RATIO"] = 79
_ELEMENT_KEYS["ERV"] = 79
_CONTACT_KEYS = {name: key for key, name in _CONTACT_NAMES.items()}
PAIR_KEY = "abaqus:contact_pair"

# The field data a step is made of (and the layout notes the reader adds).
_STEP_FIELDS = {
    TIME_KEY,
    "abaqus:step",
    "abaqus:increment",
    "abaqus:step_time",
    "abaqus:procedure",
    PAIR_KEY,
}
# Total energies (record 1999) and the element matrices (1001-1031) are field data too.
_MATRIX_KINDS = (
    ("matrix_dofs", 1002, 1002),
    ("stiffness", 1011, 1012),
    ("mass", 1021, 1022),
    ("load", 1031, 1031),
)
_MODE_FIELDS = {
    "abaqus:mode",
    "abaqus:eigenvalue",
    "abaqus:generalized_mass",
    "abaqus:composite_damping",
    "abaqus:participation_factor",
    "abaqus:effective_mass",
}
# Records 301-310 of a mode-based dynamic step: one record per row of the field.
_MODAL_KEYS = tuple(
    (key, "abaqus:" + name)
    for key, name in zip(
        range(301, 311),
        ("GU", "GV", "GA", "BM", "GPU", "GPV", "GPA", "SNE", "KE", "T"),
    )
)
# meshio++ cell type -> Abaqus element type (a results file names the exact element; these
# are the reduced-integration ones for solids and quadrilaterals)
_TYPES = {
    "line": "T3D2",
    "line3": "T3D3",
    "triangle": "S3",
    "quad": "S4R",
    "triangle6": "STRI65",
    "quad8": "S8R",
    "quad9": "S9R5",
    "tetra": "C3D4",
    "tetra10": "C3D10",
    "pyramid": "C3D5",
    "pyramid13": "C3D13",
    "wedge": "C3D6",
    "wedge15": "C3D15",
    "hexahedron": "C3D8R",
    "hexahedron20": "C3D20R",
}
_SHELLS = {"triangle", "quad", "triangle6", "quad8", "quad9"}
_NAME = re.compile(r"^(?P<base>.*?)(?P<rest>(?:@[^@]*)*)$")
_BLOCK_WORDS = 512
_CHUNK = 10  # node labels of a 1900 record, and of its 1990 continuations


def _wfail(message):
    raise WriteError(f"Abaqus .fil: {message}")


# -- records ---------------------------------------------------------------------------


def _I(value):  # noqa: N802
    return ("I", int(value))


def _D(value):  # noqa: N802
    value = float(value)
    if value != value or abs(value) in (float("inf"),):
        _wfail("a value is not finite")
    if value != 0.0 and not 1e-99 <= abs(value) < 1e100:
        _wfail("a value is outside the range a D22.15 item holds")
    return ("D", value)


def _ascii_text(text):
    """The text as the file holds it: ASCII, a ``?`` for every other byte."""
    return "".join(chr(b) if b < 128 else "?" for b in text.encode("utf-8"))


def _A(text):  # noqa: N802
    return ("A", _ascii_text(text).ljust(8)[:8])


def _text_items(text, words=10):
    """A text of up to ``8 * words`` characters as ``words`` A items."""
    text = _ascii_text(text).ljust(8 * words)
    return [_A(text[8 * k : 8 * k + 8]) for k in range(words)]


def _rec(key, *items):
    return [_I(key)] + list(items)


def _ascii(records):
    text, line = [], ""

    def put(piece):
        nonlocal line
        line += piece
        while len(line) >= 80:
            text.append(line[:80])
            line = line[80:]

    for r in records:
        put("*")
        for tag, value in [("I", len(r) + 1)] + r:
            if tag == "I":
                digits = str(value)
                put(f"I{len(digits):2d}{digits}")
            elif tag == "D":
                put("D" + f"{value: .15E}".replace("E", "D"))
            else:
                put("A" + value)
        if r[0] == ("I", 2001):
            text.append(line.ljust(80))
            text.append(" " * 80)
            line = ""
    if line:
        text.append(line)
    return ("\n".join(text) + "\n").encode()


def _binary(records, big_endian=False, marker_size=4):
    if marker_size not in (4, 8):
        _wfail("the record marker size must be 4 or 8")
    order = ">" if big_endian else "<"
    words = []
    for r in records:
        rec_words = [struct.pack(order + "q", len(r) + 1)]
        for tag, value in r:
            if tag == "I":
                rec_words.append(struct.pack(order + "q", value))
            elif tag == "D":
                rec_words.append(struct.pack(order + "d", value))
            else:
                rec_words.append(value.encode())
        if r[0] == ("I", 2001):
            # The last record of an increment pads to the end of its block.
            pad = (_BLOCK_WORDS - (len(words) + len(rec_words)) % _BLOCK_WORDS) % (
                _BLOCK_WORDS
            )
            rec_words[0] = struct.pack(order + "q", len(rec_words) + pad)
            rec_words += [b"\0" * 8] * pad
        words += rec_words
    words += [b"\0" * 8] * ((_BLOCK_WORDS - len(words) % _BLOCK_WORDS) % _BLOCK_WORDS)
    marker = struct.pack(order + ("i" if marker_size == 4 else "q"), 8 * _BLOCK_WORDS)
    out = bytearray()
    for b in range(0, len(words), _BLOCK_WORDS):
        out += marker + b"".join(words[b : b + _BLOCK_WORDS]) + marker
    return bytes(out)


# -- the model -------------------------------------------------------------------------


def _signature(mesh):
    return (len(mesh.points), [(c.type, np.asarray(c.data).shape) for c in mesh.cells])


def _labels(mesh, name, count, what):
    """``abaqus:id`` per entity, else 1 to n."""
    if name is None:
        return list(range(1, count + 1))
    labels = [int(v) for v in name]
    if len(labels) != count:
        _wfail(f"'abaqus:id' must hold one label per {what}")
    if any(v <= 0 for v in labels):
        _wfail(f"{what} labels must be positive")
    if len(set(labels)) != len(labels):
        _wfail(f"{what} labels are not unique")
    return labels


class _Model:
    """The mesh as the file holds it: records, labels and the written cells."""

    def __init__(self, mesh):
        self.signature = _signature(mesh)
        self.dropped = set()
        # contact output refers to its surface pair by a 1940 label
        self.surface_labels = {}
        self.contact_names = set()
        if any(_contact_key(n) for n in mesh.point_data):
            pair = _contact_pair(mesh)
            if pair is not None:
                self.contact_names = {pair[0].name, pair[1].name}
        n = len(mesh.points)
        ids = mesh.point_data.get("abaqus:id")
        self.node_labels = _labels(mesh, ids, n, "node")
        raw = np.asarray(mesh.points, dtype=np.float64).reshape(n, -1)
        self.dim = 2 if raw.shape[1] == 2 else 3
        points = np.zeros((n, self.dim))
        points[:, : min(self.dim, raw.shape[1])] = raw[:, : self.dim]
        # every cell of the mesh, block-major: (label, written?, block, row)
        total = sum(len(c.data) for c in mesh.cells)
        cell_ids = mesh.cell_data.get("abaqus:id")
        flat_ids = None
        if cell_ids is not None:
            flat_ids = np.concatenate([np.asarray(a).ravel() for a in cell_ids])
        self.cell_labels = _labels(mesh, flat_ids, total, "element")
        self.cells = []  # written: (global index, type name, block, row)
        skipped = {}
        g = 0
        for b, block in enumerate(mesh.cells):
            name = _TYPES.get(block.type)
            for r in range(len(block.data)):
                if name is None:
                    skipped[block.type] = skipped.get(block.type, 0) + 1
                else:
                    self.cells.append((g, name, b, r))
                g += 1
        self.skipped = skipped
        self.written = {g for g, *_ in self.cells}
        self.records = self._records(mesh, points)

    def _records(self, mesh, points):
        out = [
            _rec(
                1921,
                _A("meshio++"),
                _A(""),
                _A(""),
                _A(""),
                _I(len(self.cells)),
                _I(len(self.node_labels)),
                _D(1.0),
            ),
            _rec(1922, *_text_items("meshio++ results file")),
            # the active degrees of freedom (1, 2 and, in 3-D, 3); 35 words
            _rec(1902, *map(_I, [1, 2] + ([3] if self.dim == 3 else [0]) + [0] * 32)),
        ]
        for g, name, b, r in self.cells:
            row = [int(p) for p in mesh.cells[b].data[r]]
            if any(not 0 <= p < len(self.node_labels) for p in row):
                _wfail("a cell refers to a point the mesh does not have")
            conn = [self.node_labels[p] for p in row]
            first = conn[:_CHUNK]
            out.append(_rec(1900, _I(self.cell_labels[g]), _A(name), *map(_I, first)))
            for k in range(_CHUNK, len(conn), _CHUNK):
                out.append(_rec(1990, *map(_I, conn[k : k + _CHUNK])))
        for label, xyz in zip(self.node_labels, points):
            out.append(_rec(1901, _I(label), *[_D(x) for x in xyz]))
        out += self._sets(mesh)
        out += self._matrices(mesh)
        return out

    def _matrices(self, mesh):
        """The element matrices: per row of an ``abaqus:<kind>:index`` table
        (element, offset, length, flag) a 1001 record and its values."""
        out = []
        fd = mesh.field_data
        for kind, key_flag, key_plain in _MATRIX_KINDS:
            values = fd.get("abaqus:" + kind)
            index = fd.get("abaqus:" + kind + ":index")
            if values is None or index is None:
                continue
            values = np.asarray(values).ravel()
            index = np.asarray(index, dtype=np.int64).reshape(-1, 4)
            for element, offset, length, flag in index.tolist():
                if offset < 0 or length < 0 or offset + length > len(values):
                    _wfail(f"'abaqus:{kind}:index' reaches past 'abaqus:{kind}'")
                data = values[offset : offset + length]
                out.append(_rec(1001, _I(element)))
                if kind == "matrix_dofs":
                    out.append(_rec(1002, *[_I(int(v)) for v in data]))
                elif kind == "load":
                    out.append(_rec(1031, _I(int(flag)), *[_D(v) for v in data]))
                else:
                    out.append(
                        _rec(
                            key_flag if flag == 1 else key_plain, *[_D(v) for v in data]
                        )
                    )
        return out

    def _sets(self, mesh):
        regions = [r for r in getattr(mesh, "regions", None) or []]
        used_digits = {int(r.name) for r in regions if r.name.isdigit()}
        refs = []  # 1940 records
        out = []
        label = 0
        for region in sorted(regions, key=lambda r: (r.kind, r.name)):
            name = _ascii_text(region.name)
            if len(name) > 80:
                _wfail(f"region name '{name[:20]}...' is longer than 80 characters")
            if len(name) > 8 or (
                region.kind == "side" and region.name in self.contact_names
            ):
                label += 1
                while label in used_digits:
                    label += 1
                refs.append(_rec(1940, _I(label), *_text_items(name)))
                if region.kind == "side":
                    self.surface_labels[region.name] = label
                name = str(label)
            if region.kind == "side":
                out += self._surface(mesh, region, name)
                continue
            entries = [int(e) for e in np.asarray(region.entries).ravel()]
            if region.kind == "point":
                if any(not 0 <= e < len(self.node_labels) for e in entries):
                    _wfail(
                        f"region '{region.name}' names a point the mesh does not have"
                    )
                values = [self.node_labels[e] for e in entries]
                keys = (1931, 1932)
            else:
                values = [self.cell_labels[e] for e in entries if e in self.written]
                keys = (1933, 1934)
            out.append(_rec(keys[0], _A(name), *map(_I, values[:5])))
            for k in range(5, len(values), _CHUNK):
                out.append(_rec(keys[1], *map(_I, values[k : k + _CHUNK])))
        return out + refs

    def _surface(self, mesh, region, name):
        """A side region as a contact surface: a 1501 record and a 1502 record
        per (element, face number)."""
        types = {g: mesh.cells[b].type for g, _n, b, _r in self.cells}
        faces = []
        dropped = 0
        for cell, facet in np.asarray(region.entries, dtype=np.int64).reshape(-1, 2):
            cell, facet = int(cell), int(facet)
            if cell not in self.written:
                dropped += 1
                continue
            label = _face_name(types[cell], facet)
            number = int(label[1:]) if label else None
            if number is None:
                dropped += 1
                continue
            faces.append((self.cell_labels[cell], number))
        if dropped:
            self.dropped.add(f"{dropped} facet(s) of surface '{region.name}'")
        # name, dimension, type (deformable), facet count, slave count
        out = [_rec(1501, _A(name), _I(self.dim), _I(1), _I(len(faces)), _I(0))]
        out += [_rec(1502, _I(label), _I(number)) for label, number in faces]
        return out

    def cell_index(self):
        """Written cell -> its position in the file."""
        return {g: k for k, (g, *_rest) in enumerate(self.cells)}


# -- one increment ---------------------------------------------------------------------


def _trim(values):
    """A record's values: trailing NaN dropped, an inner NaN written as 0."""
    values = [float(v) for v in values]
    while values and values[-1] != values[-1]:
        values.pop()
    return [0.0 if v != v else v for v in values]


def _split_name(name):
    """``(base, section point, rebar, collision suffix)`` of a cell array name."""
    base, rest = _NAME.match(name).group("base", "rest")
    section, rebar, collision = None, None, None
    for part in [p for p in rest.split("@") if p]:
        if part.startswith("sp") and part[2:].isdigit():
            section = int(part[2:])
        elif part.startswith("rebar:"):
            rebar = part[len("rebar:") :]
        elif part in ("ip", "el"):
            collision = part
        else:
            return name, None, None, None
    return base, section, rebar, collision


def _contact_key(name):
    """The contact output record (1505-1599) a point array is, or None: the
    reader's contact names, and ``key_<n>`` for a key it has no name for."""
    if name in _CONTACT_KEYS:
        return _CONTACT_KEYS[name]
    if name.startswith("key_") and name[4:].isdigit():
        key = int(name[4:])
        if 1505 <= key <= 1599 and key not in _CONTACT_NAMES:
            return key
    return None


def _contact_pair(mesh):
    """The (slave, master) side regions contact output is written for: the two
    indices of ``abaqus:contact_pair`` into ``mesh.regions``, else the first two
    side regions. None when the mesh has no pair."""
    regions = list(getattr(mesh, "regions", None) or [])
    given = mesh.field_data.get(PAIR_KEY)
    if given is not None:
        pair = [int(v) for v in np.asarray(given).ravel()]
        if len(pair) != 2 or any(
            not 0 <= v < len(regions) or regions[v].kind != "side" for v in pair
        ):
            _wfail(f"'{PAIR_KEY}' must hold the indices of two side regions")
        return regions[pair[0]], regions[pair[1]]
    sides = [r for r in regions if r.kind == "side"]
    return (sides[0], sides[1]) if len(sides) > 1 else None


def _key_of(table, base):
    if base in table:
        return table[base]
    if base.startswith("key_") and base[4:].isdigit():
        return int(base[4:])
    return None


def _ndi_nshr(count):
    return {6: (3, 3), 4: (3, 1), 3: (2, 1)}.get(count, (0, 0))


class _Increment:
    def __init__(self, model, mesh, time, index):
        self.records = []
        self.dropped = set()
        self.time = float(time)
        fd = mesh.field_data

        def scalar(name, default, cast):
            array = fd.get(name)
            if array is None or np.asarray(array).size == 0:
                return default
            return cast(np.asarray(array).ravel()[0])

        step_time = scalar("abaqus:step_time", self.time, float)
        procedure = scalar("abaqus:procedure", 1, int)
        step = scalar("abaqus:step", 1, int)
        increment = scalar("abaqus:increment", index + 1, int)
        self.records.append(
            _rec(
                2000,
                _D(self.time),
                _D(step_time),
                _D(0.0),
                _D(0.0),
                _I(procedure),
                _I(step),
                _I(increment),
                _I(0),
                _D(0.0),
                _D(0.0),
                _D(1.0),
                *_text_items(""),
            )
        )
        self._mode(fd)
        self._energies(fd, procedure)
        self._modal(fd)
        element = self._element(model, mesh)
        nodal = self._nodal(model, mesh)
        contact = self._contact(model, mesh)
        if element:
            self.records.append(_rec(1911, _I(0), _A(""), _A("")))
            self.records += element
        if nodal:
            self.records.append(_rec(1911, _I(1), _A("")))
            self.records += nodal
        self.records += contact
        self.records.append(_rec(2001))
        handled = {"abaqus:" + n for n in _ENERGY_NAMES + _EXPLICIT_ENERGY_NAMES if n}
        handled |= {name for _, name in _MODAL_KEYS}
        for kind, _, _ in _MATRIX_KINDS:
            handled |= {"abaqus:" + kind, "abaqus:" + kind + ":index"}
        for name in sorted(fd):
            if name in _STEP_FIELDS or name in _MODE_FIELDS or name in handled:
                continue
            if name.startswith("abaqus:layout:"):
                continue
            self.dropped.add(name)

    def _energies(self, fd, procedure):
        """Record 1999, the total energies, positional: those the step's
        procedure names, up to the last one present (an absent one in between
        is 0)."""
        names = (
            _EXPLICIT_ENERGY_NAMES
            if procedure in _EXPLICIT_PROCEDURES
            else _ENERGY_NAMES
        )
        present = [k for k, n in enumerate(names) if n and "abaqus:" + n in fd]
        if not present:
            return
        values = [0.0] * (max(present) + 1)
        for k in present:
            values[k] = float(np.asarray(fd["abaqus:" + names[k]]).ravel()[0])
        self.records.append(_rec(1999, *map(_D, values)))

    def _modal(self, fd):
        """Records 301-310, the generalized quantities of a mode-based dynamic
        step: one record per row of ``abaqus:GU`` and the like (trailing NaN
        trimmed, an inner NaN written as 0). BM's row is the motion type and
        six components, followed by the base name (blank)."""
        for key, name in _MODAL_KEYS:
            if name not in fd:
                continue
            array = np.asarray(fd[name], dtype=float)
            if array.ndim > 2:
                _wfail(f"'{name}' must be one row or a matrix of rows")
            for row in np.atleast_2d(array):
                values = list(row)
                while values and values[-1] != values[-1]:
                    values.pop()
                if not values:
                    continue
                values = [0.0 if v != v else float(v) for v in values]
                if key == 304:
                    values = values[:7]
                    items = [_I(values[0]), *map(_D, values[1:])]
                    if len(values) == 7:
                        items.append(_A(""))
                else:
                    items = [_D(v) for v in values]
                self.records.append(_rec(key, *items))

    def _mode(self, fd):
        if "abaqus:mode" not in fd:
            return

        def get(name, default=0.0):
            array = fd.get(name)
            return default if array is None else np.asarray(array, dtype=float).ravel()

        mode = int(np.asarray(fd["abaqus:mode"]).ravel()[0])
        factors, masses = get("abaqus:participation_factor", np.zeros(0)), get(
            "abaqus:effective_mass", np.zeros(0)
        )
        if len(factors) != len(masses):
            _wfail("'abaqus:participation_factor' and 'abaqus:effective_mass' differ")
        pairs = [_D(v) for pair in zip(factors, masses) for v in pair]
        scalars = [
            float(np.asarray(get(n)).ravel()[0]) if np.size(get(n)) else 0.0
            for n in (
                "abaqus:eigenvalue",
                "abaqus:generalized_mass",
                "abaqus:composite_damping",
            )
        ]
        self.records.append(_rec(1980, _I(mode), *[_D(v) for v in scalars], *pairs))

    def _nodal(self, model, mesh):
        out = []
        for name in sorted(mesh.point_data):
            if name == "abaqus:id":
                continue
            if name.endswith("@avg"):
                continue  # written with the element records
            if _contact_key(name) is not None:
                continue  # contact output, written by _contact
            key = _key_of(_NODAL_KEYS, name)
            if key is None:
                if _key_of(_ELEMENT_KEYS, name) is None:
                    self.dropped.add(name)
                continue
            array = np.asarray(mesh.point_data[name], dtype=np.float64)
            n = len(model.node_labels)
            if array.size == 0 or array.shape[0] != n:
                self.dropped.add(name)
                continue
            rows = array.reshape(n, -1)
            for label, row in zip(model.node_labels, rows):
                values = _trim(row)
                if values:
                    out.append(_rec(key, _I(label), *[_D(v) for v in values]))
        return out

    def _contact(self, model, mesh):
        """Contact output: a 1503 request naming the slave and master surface
        by their 1940 labels, then per node holding data a 1504 header and a
        record per array (the reader's contact names, ``key_<n>``)."""
        names = sorted(
            (n for n in mesh.point_data if _contact_key(n) is not None),
            key=lambda n: (_contact_key(n), n),
        )
        if not names:
            return []
        pair = _contact_pair(mesh)
        labels = [model.surface_labels.get(r.name) for r in pair or ()]
        if len(labels) != 2 or None in labels:
            self.dropped.update(names)
            return []
        n = len(model.node_labels)
        rows = {}
        for name in names:
            array = np.asarray(mesh.point_data[name], dtype=np.float64)
            if array.size == 0 or array.shape[0] != n:
                self.dropped.add(name)
                continue
            rows[name] = array.reshape(n, -1)
        out = []
        for i, label in enumerate(model.node_labels):
            items = []
            for name in names:
                values = _trim(rows[name][i]) if name in rows else []
                if values:
                    items.append(_rec(_contact_key(name), *[_D(v) for v in values]))
            if items:
                out.append(_rec(1504, _I(label), _I(len(items[0]) - 1)))
                out += items
        if out:
            out.insert(
                0, _rec(1503, _I(0), _A(str(labels[0])), _A(str(labels[1])), _A(""))
            )
        return out

    def _element(self, model, mesh):
        index = model.cell_index()
        slots = {}  # (position, location, section, point, rebar) -> [(key, values)]
        offsets = [0]
        for block in mesh.cells:
            offsets.append(offsets[-1] + len(block.data))
        for name in sorted(mesh.cell_data):
            if name == "abaqus:id":
                continue
            base, section, rebar, collision = _split_name(name)
            key = _key_of(_ELEMENT_KEYS, base)
            if key is None:
                self.dropped.add(name)
                continue
            layout = mesh.field_data.get("abaqus:layout:" + name)
            if layout is not None:
                layout = [int(v) for v in np.asarray(layout).ravel()]
                if len(layout) != 2 or layout[0] < 1 or layout[1] < 1:
                    _wfail(f"'abaqus:layout:{name}' must be (points, components)")
            location = 1
            if rebar is not None:
                location = 3
            elif layout is not None:
                location = 2 if (base == "NFORC" or collision == "el") else 0
            for b, array in enumerate(mesh.cell_data[name]):
                block = mesh.cells[b]
                array = np.asarray(array, dtype=np.float64)
                n = len(block.data)
                if n == 0 or array.shape[0] != n:
                    continue
                rows = array.reshape(n, -1)
                sec = section if section is not None else int(block.type in _SHELLS)
                if layout is not None and rows.shape[1] != layout[0] * layout[1]:
                    _wfail(f"'{name}' does not match its 'abaqus:layout:{name}'")
                for r in range(n):
                    g = offsets[b] + r
                    if g not in index:
                        continue
                    if layout is None:
                        values = _trim(rows[r])
                        if values:
                            slots.setdefault(
                                (index[g], location, sec, 0, ""), []
                            ).append((key, values))
                        continue
                    for p in range(layout[0]):
                        values = _trim(rows[r, p * layout[1] : (p + 1) * layout[1]])
                        if not values:
                            continue
                        if location == 2:
                            nodes = block.data[r]
                            if p >= len(nodes):
                                continue
                            point = model.node_labels[int(nodes[p])]
                        else:
                            point = p + 1
                        slots.setdefault(
                            (index[g], location, sec, point, rebar or ""), []
                        ).append((key, values))
        out = []
        for (position, location, sec, point, rebar), items in sorted(slots.items()):
            g = model.cells[position][0]
            counts = [len(v) for _, v in items if len(v) in (3, 4, 6)]
            ndi, nshr = _ndi_nshr(counts[0] if counts else 0)
            out.append(
                _rec(
                    1,
                    _I(model.cell_labels[g]),
                    _I(point),
                    _I(sec),
                    _I(location),
                    _A(rebar),
                    _I(ndi),
                    _I(nshr),
                    _I(0),
                    _I(0),
                )
            )
            for key, values in items:
                out.append(_rec(key, *[_D(v) for v in values]))
        out += self._averaged(model, mesh)
        return out

    def _averaged(self, model, mesh):
        """Point data named like an element result: averaged at the nodes (4)."""
        out = []
        for name in sorted(mesh.point_data):
            base = name
            while base.endswith("@avg"):
                base = base[: -len("@avg")]
            if base == "abaqus:id":
                continue
            if base == name and _contact_key(base) is not None:
                continue  # contact output
            if base == name and _key_of(_NODAL_KEYS, base) is not None:
                continue  # a nodal result
            key = _key_of(_ELEMENT_KEYS, base)
            if key is None:
                continue
            array = np.asarray(mesh.point_data[name], dtype=np.float64)
            n = len(model.node_labels)
            if array.size == 0 or array.shape[0] != n:
                continue
            rows = array.reshape(n, -1)
            for label, row in zip(model.node_labels, rows):
                values = _trim(row)
                if not values:
                    continue
                out.append(
                    _rec(
                        1,
                        _I(label),
                        _I(0),
                        _I(0),
                        _I(4),
                        _A(""),
                        _I(0),
                        _I(0),
                        _I(0),
                        _I(0),
                    )
                )
                out.append(_rec(key, *[_D(v) for v in values]))
        return out


def _has_results(mesh):
    return TIME_KEY in mesh.field_data or any(
        n != "abaqus:id" for n in list(mesh.point_data) + list(mesh.cell_data)
    )


def _warn_dropped(model, increments):
    names = set(model.dropped)
    for inc in increments:
        names |= inc.dropped
    if names:
        warn(
            "Abaqus .fil: arrays and sets without a place in the file are not "
            "written: " + ", ".join(sorted(names))
        )
    if model.skipped:
        warn(
            "Abaqus .fil: cell blocks with no Abaqus element type are not written: "
            + ", ".join(sorted(model.skipped))
        )


def _encode(binary, records, big_endian=False, marker_size=4):
    if binary:
        return _binary(records, big_endian, marker_size)
    return _ascii(records)


def _time_of(mesh):
    array = mesh.field_data.get(TIME_KEY)
    if array is None or np.asarray(array).size == 0:
        return 0.0
    return float(np.asarray(array, dtype=np.float64).ravel()[0])


def _save(filename, data):
    if is_buffer(filename, "w"):
        filename.write(data)
    else:
        with open(filename, "wb") as fh:
            fh.write(data)


@_provenance.slotless_writer
def write(filename, mesh, binary=False, big_endian=False, marker_size=4):
    model = _Model(mesh)
    increments = []
    if _has_results(mesh):
        increments.append(_Increment(model, mesh, _time_of(mesh), 0))
    _warn_dropped(model, increments)
    records = list(model.records)
    for inc in increments:
        records += inc.records
    _save(filename, _encode(binary, records, big_endian, marker_size))


class SeriesWriter:
    """Write a time series into one results file: the model with the first step
    and an increment per step (a 2000 record, the step's nodal and element
    results and a 2001 record), which ``read`` gives back by ``time_step``. A
    step's time is the argument, not its ``meshio:time``; its step number,
    increment, step time and procedure are its ``abaqus:*`` field data. Every
    step must have the first step's points and cell blocks. The file is
    assembled when the writer closes. Use as a context manager, or call
    :meth:`close`."""

    def __init__(self, filename, binary=False, big_endian=False, marker_size=4):
        self._filename = filename
        self._binary = bool(binary)
        self._big_endian = bool(big_endian)
        self._marker_size = marker_size
        self._model = None
        self._increments = []

    def write(self, time, mesh):
        index = len(self._increments)
        if self._model is None:
            self._model = _Model(mesh)
        elif self._model.signature != _signature(mesh):
            _wfail(
                f"series: step {index}'s points or cells differ from the first "
                "step's; a results file holds one model -- write one file per "
                "step with '{step}'"
            )
        self._increments.append(_Increment(self._model, mesh, float(time), index))

    def close(self):
        if self._filename is None:
            return
        filename, self._filename = self._filename, None
        if not self._increments:
            return
        _warn_dropped(self._model, self._increments)
        records = list(self._model.records)
        for inc in self._increments:
            records += inc.records
        _save(
            filename,
            _encode(self._binary, records, self._big_endian, self._marker_size),
        )

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        if exc[0] is None:
            self.close()
        else:
            self._filename = None
        return False


__all__ = ["write", "SeriesWriter"]

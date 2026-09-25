"""MSC Marc input decks (``.dat``) and formatted post files (``.t19``): the
pure-Python reference readers.

The twins of ``src/cpp/src/formats/marc.cpp``. A deck's model definition --
``COORDINATES``, ``CONNECTIVITY`` and ``DEFINE`` sets -- makes the mesh; the
parameter section only says whether fields are ``EXTENDED`` (twice as wide),
and history definition is not read. A ``.t19`` post file is a sequence of
``=beg=5xxnn (name)`` ... ``=end=`` blocks: the model, then one group of blocks
per increment between ``****`` and ``----``. The layouts follow Marc Volume C
(program input) and Volume D (the PLDUMP2000 description of the post file); see
``doc/formats/marc.md``.
"""

from __future__ import annotations

import os
import re

import numpy as np

from .._common import num_nodes_per_cell, warn
from .._exceptions import ReadError
from .._files import is_buffer, open_file
from .._mesh import Mesh, topological_dimension
from .._regions import Region
from ..lsdyna._lsdyna import _collapse_solid

_DAT = "Marc .dat"
_T19 = "Marc .t19"

# Marc element type -> (meshio++ cell type, the element's node count), from
# Volume B. Types whose topology the manuals were not checked for are left out:
# their elements are skipped with a warning. Every type lists its nodes in
# meshio++'s order (corners with the face 1-2-3(-4) normal pointing into the
# element, then mid-edge nodes bottom ring, top ring, verticals), so no
# permutation applies.
TYPES = {}
# Plain types, and (v16.12.0) the Herrmann (mixed) types whose pressure sits at
# the corners, the rebar and the composite types, all with the node lists of
# their plain twins.
for _t in (3, 10, 11, 18, 75, 139, 140, 143, 144, 145, 147, 151, 152):
    TYPES[_t] = ("quad", 4)
for _t in (2, 6, 138, 158, 201):
    TYPES[_t] = ("triangle", 3)
for _t in (22, 26, 27, 28, 30, 32, 33, 46, 48, 53, 54, 55, 58, 59, 63, 66, 142, 148):
    TYPES[_t] = ("quad8", 8)
for _t in (153, 154):
    TYPES[_t] = ("quad8", 8)
for _t in (124, 125, 126, 128, 129, 200):
    TYPES[_t] = ("triangle6", 6)
for _t in (7, 43, 117, 123, 146, 149):
    TYPES[_t] = ("hexahedron", 8)
for _t in (21, 23, 35, 44, 57, 61, 150):
    TYPES[_t] = ("hexahedron20", 20)
for _t in (134, 135):
    TYPES[_t] = ("tetra", 4)
for _t in (127, 130, 133):
    TYPES[_t] = ("tetra10", 10)
for _t in (136, 137):
    TYPES[_t] = ("wedge", 6)
for _t in (9, 31, 52, 98, 165, 166, 167):
    TYPES[_t] = ("line", 2)
for _t in (64, 168, 169, 170):
    TYPES[_t] = ("line3", 3)
# Elements with more nodes than their cell keeps: the leading geometric nodes
# are the cell, the rest (a Herrmann pressure node, a centroid bubble node,
# generalized plane strain nodes) are dropped.
for _t in (80, 82, 83, 118, 119):
    TYPES[_t] = ("quad", 5)
TYPES[81] = ("quad", 7)
for _t in (34, 47, 60):
    TYPES[_t] = ("quad8", 10)
for _t in (84, 120):
    TYPES[_t] = ("hexahedron", 9)
for _t in (155, 156):
    TYPES[_t] = ("triangle", 4)
TYPES[157] = ("tetra", 5)
del _t

# The first word of a line that can open a deck (a parameter), and the model
# definition options this reader acts on.
_PARAMETERS = {
    "title",
    "sizing",
    "elements",
    "extended",
    "version",
    "table",
    "processor",
    "alloc",
    "setname",
    "dist",
    "large",
    "update",
    "finite",
    "all",
    "no",
    "state",
    "heat",
    "coupled",
    "harmonic",
    "buckle",
    "dynamic",
    "fluid",
    "electrostatic",
    "magnetostatic",
    "joule",
    "bearing",
    "shell",
    "print",
    "lumping",
    "assumed",
    "constant",
    "feature",
    "follow",
    "plasticity",
    "rezoning",
    "scale",
    "structural",
    "thermal",
    "end",
}

# An element's first line holds its number, type and 14 nodes; more continue.
_FIRST_LINE_NODES = 14


def _fail(label, what):
    raise ReadError(f"{label}: {what}")


# -- content check -------------------------------------------------------------


def is_marc_deck(text):
    """Whether the start of a ``.dat`` file (``text``) is a Marc input deck: its
    first keyword line is a Marc parameter (no ``=``, unlike a Tecplot header)
    and a CONNECTIVITY, COORDINATES or END line follows."""
    first = None
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("$"):
            continue
        word = re.split(r"[\s,]+", stripped.lower(), maxsplit=1)[0]
        if first is None:
            if "=" in stripped or word not in _PARAMETERS:
                return False
            first = word
            if word == "end":
                return True
            continue
        if word in ("end", "connectivity", "coordinates"):
            return True
    return False


# -- deck fields ---------------------------------------------------------------


def _real(text, label, where):
    """A real field; blank is 0; ``D`` exponents and Fortran's ``1.5-3``."""
    text = text.strip()
    if not text:
        return 0.0
    t = text.replace("D", "E").replace("d", "e")
    try:
        return float(t)
    except ValueError:
        k = max(t.rfind("+"), t.rfind("-"))
        if k > 0 and t[k - 1] not in "eE":
            try:
                return float(t[:k] + "E" + t[k:])
            except ValueError:
                pass
    _fail(label, f"bad real field '{text}' ({where})")


def _int(text, label, where):
    text = text.strip()
    if not text:
        return 0
    try:
        return int(text)
    except ValueError:
        _fail(label, f"bad integer field '{text}' ({where})")


def _fields(line, width):
    """The fields of a data line: comma-separated (free format), or fixed
    columns of ``width`` up to the last non-blank one."""
    line = line.rstrip("\r\n")
    if "," in line:
        parts = [p.strip() for p in line.split(",")]
        while parts and not parts[-1]:
            parts.pop()  # a lone item is followed by a comma
        return parts
    body = line.rstrip()
    return [body[k : k + width].strip() for k in range(0, len(body), width)]


def _is_data(line):
    """A data line starts (after blanks) with a number; a keyword with a letter."""
    s = line.lstrip()
    return bool(s) and (s[0].isdigit() or s[0] in "+-.")


def _is_comment(line):
    s = line.strip()
    return not s or s.startswith("$")


# -- the input deck ------------------------------------------------------------


class _Deck:
    def __init__(self):
        self.extended = False
        self.ncoord = 3
        self.nodes = {}  # id -> xyz (the last definition wins)
        self.node_order = []
        self.elements = []  # (id, type, nodes)
        self.sets = []  # (name, "element" | "node", tokens)

    @property
    def int_width(self):
        return 10 if self.extended else 5

    @property
    def real_width(self):
        return 20 if self.extended else 10


def _parse_deck(lines):
    deck = _Deck()
    n = len(lines)
    i = 0
    in_parameters = True
    while i < n:
        line = lines[i]
        if _is_comment(line) or _is_data(line) or not line[:1].strip():
            i += 1
            continue
        low = line.strip().lower()
        words = re.split(r"[\s,]+", low)
        key = words[0]
        if in_parameters:
            if key == "extended":
                deck.extended = True
            elif key == "end" and (len(words) < 2 or words[1] != "option"):
                in_parameters = False
                i += 1
                continue
            if key not in ("connectivity", "coordinates", "define"):
                i += 1
                continue
            in_parameters = False  # a deck without END: the model starts here
        if key == "end" and len(words) > 1 and words[1] == "option":
            break
        if key == "connectivity":
            i = _parse_connectivity(deck, lines, i + 1)
        elif key == "coordinates":
            i = _parse_coordinates(deck, lines, i + 1)
        elif key == "define":
            i = _parse_define(deck, lines, i, words)
        else:
            i += 1
    return deck


def _next_data(lines, i):
    """The index of the next non-comment line (or ``len(lines)``)."""
    while i < len(lines) and _is_comment(lines[i]):
        i += 1
    return i


def _parse_connectivity(deck, lines, i):
    i = _next_data(lines, i)
    if i < len(lines) and _is_data(lines[i]):
        i += 1  # the header line: element count, unit, print flag ...
    width = deck.int_width
    while True:
        i = _next_data(lines, i)
        if i >= len(lines) or not _is_data(lines[i]):
            return i
        where = f"line {i + 1}"
        f = _fields(lines[i], width)
        if len(f) < 2:
            _fail(_DAT, f"an element needs a number and a type ({where})")
        ident, etype = _int(f[0], _DAT, where), _int(f[1], _DAT, where)
        known = TYPES.get(etype)
        nodes = [_int(v, _DAT, where) for v in f[2:]]
        i += 1
        if known is None:
            if len(f) >= 2 + _FIRST_LINE_NODES:
                _fail(
                    _DAT,
                    f"element {ident} is of type {etype}, which meshio++ does not "
                    "know the node count of (its nodes continue on further lines)",
                )
            deck.elements.append((ident, etype, nodes))
            continue
        count = known[1]
        while len(nodes) < count:
            i = _next_data(lines, i)
            if i >= len(lines) or not _is_data(lines[i]):
                _fail(_DAT, f"element {ident} lists {len(nodes)} of its {count} nodes")
            where = f"line {i + 1}"
            nodes += [_int(v, _DAT, where) for v in _fields(lines[i], width)]
            i += 1
        deck.elements.append((ident, etype, nodes[:count]))


def _parse_coordinates(deck, lines, i):
    i = _next_data(lines, i)
    if i < len(lines) and _is_data(lines[i]):
        header = _fields(lines[i], deck.int_width)
        if header and header[0]:
            deck.ncoord = max(1, _int(header[0], _DAT, f"line {i + 1}"))
        i += 1
    ncoord = deck.ncoord
    iw, rw = deck.int_width, deck.real_width
    while True:
        i = _next_data(lines, i)
        if i >= len(lines) or not _is_data(lines[i]):
            return i
        where = f"line {i + 1}"
        line = lines[i].rstrip("\r\n")
        if "," in line:
            f = _fields(line, iw)
            ident, values = _int(f[0], _DAT, where), f[1:]
        else:
            ident = _int(line[:iw], _DAT, where)
            body = line[iw:].rstrip()
            values = [body[k : k + rw] for k in range(0, len(body), rw)]
        i += 1
        while len(values) < ncoord:  # continuation lines, six reals each
            i = _next_data(lines, i)
            if i >= len(lines) or not _is_data(lines[i]):
                break
            more = lines[i].rstrip("\r\n")
            if "," in more:
                values += _fields(more, rw)
            else:
                body = more.rstrip()
                values += [body[k : k + rw] for k in range(0, len(body), rw)]
            i += 1
        xyz = [_real(v, _DAT, where) for v in values[:ncoord]]
        xyz += [0.0] * (3 - len(xyz))
        if ident not in deck.nodes:
            deck.node_order.append(ident)
        deck.nodes[ident] = xyz[:3]


def _set_tokens(text, width):
    """The items of a set's data line: words and numbers (touching fixed-width
    integers split apart); ``c`` or ``continue`` last means more lines follow."""
    out = []
    for tok in re.split(r"[\s,]+", text.strip()):
        if not tok:
            continue
        if tok.isdigit() and len(tok) > width and len(tok) % width == 0:
            out += [tok[k : k + width] for k in range(0, len(tok), width)]
        else:
            out.append(tok.lower())
    return out


def _is_integer(text):
    return text.lstrip("+-").isdigit()


def _parse_define(deck, lines, i, words):
    where = f"line {i + 1}"
    raw = [w for w in re.split(r"[\s,]+", lines[i].strip()) if w]
    kind = raw[1].lower() if len(raw) > 1 else ""
    k = 2
    if len(raw) > k and raw[k].lower() in ("set", "oset"):
        k += 1
    name = raw[k] if len(raw) > k else ""
    earlier = {n.lower() for n, _, _, _ in deck.sets}
    i += 1
    tokens = []
    while True:
        i = _next_data(lines, i)
        if i >= len(lines):
            break
        line = lines[i]
        items = _set_tokens(line, deck.int_width)
        # The first data line starts with a number or an earlier set's name;
        # later ones follow a line that ended in C (continue).
        if not items or (not tokens and not _is_data(line) and items[0] not in earlier):
            break
        tokens += items
        i += 1
        if tokens[-1] not in ("c", "continue"):
            break
        tokens.pop()
    if kind in ("edge", "face"):
        # `elem:number` members.
        pairs = []
        for t in tokens:
            a, colon, b = t.partition(":")
            if not colon or not _is_integer(a) or not _is_integer(b):
                warn(
                    f"{_DAT}: DEFINE {kind.upper()} SET '{name}': '{t}' is not an "
                    f"element:number pair ({where})"
                )
                continue
            pairs.append((int(a), int(b)))
        deck.sets.append((name, kind, pairs, where))
    elif kind in ("element", "elsq", "node", "ndsq"):
        family = "element" if kind in ("element", "elsq") else "node"
        deck.sets.append((name, family, tokens, where))
    else:
        warn(f"{_DAT}: DEFINE {kind.upper()} SET '{name}' is not read ({where})")
    return i


def _expand(tokens, name, known, label):
    """The members of a set: numbers, ``a TO b [BY c]`` ranges and other
    sets' names, combined left to right by AND (the default), EXCEPT and
    INTERSECT."""
    out = []
    op = "and"
    k = 0

    def combine(items):
        nonlocal out
        if op == "and":
            have = set(out)
            out += [v for v in items if v not in have and not have.add(v)]
        elif op == "except":
            drop = set(items)
            out = [v for v in out if v not in drop]
        else:
            keep = set(items)
            out = [v for v in out if v in keep]

    while k < len(tokens):
        tok = tokens[k]
        if tok in ("and", "except", "intersect"):
            op = tok
            k += 1
            continue
        if tok.lstrip("+-").isdigit():
            start = int(tok)
            k += 1
            if k < len(tokens) and tokens[k] in ("to", "through"):
                if k + 1 >= len(tokens):
                    _fail(label, f"set '{name}' ends in a range with no end")
                stop = int(tokens[k + 1])
                k += 2
                step = 1
                if k < len(tokens) and tokens[k] == "by":
                    if k + 1 >= len(tokens):
                        _fail(label, f"set '{name}' ends with BY and no step")
                    step = abs(int(tokens[k + 1])) or 1
                    k += 2
                items = (
                    list(range(start, stop + 1, step))
                    if stop >= start
                    else list(range(start, stop - 1, -step))
                )
            else:
                items = [start]
            combine(items)
            op = "and"
            continue
        if tok in known:
            combine(known[tok])
            op = "and"
            k += 1
            continue
        _fail(label, f"set '{name}' names '{tok}', which is not an earlier set")
    return out


def _build(label, node_ids, coords, elements, sets):
    """The mesh of nodes, elements ``(id, type, nodes)`` and sets ``(name,
    "element" | "node", member ids)``; also each element's ``(block, row)``
    (``None`` when it has no cell)."""
    index = {ident: k for k, ident in enumerate(node_ids)}
    points = np.array(coords, dtype=np.float64).reshape(-1, 3)
    blocks, order = {}, []
    skipped = {}
    degenerate20 = 0
    locs = []
    element_cell = {}
    for ident, etype, nodes in elements:
        known = TYPES.get(etype)
        if known is None:
            skipped[etype] = skipped.get(etype, 0) + 1
            locs.append(None)
            continue
        cell_type, _ = known
        nodes = list(nodes)[: num_nodes_per_cell[cell_type]]
        # The 3-node rebar lines list their middle node second.
        if 168 <= etype <= 170 and len(nodes) == 3:
            nodes = [nodes[0], nodes[2], nodes[1]]
        for v in nodes:
            if v not in index:
                _fail(label, f"element {ident} names undefined node {v}")
        if cell_type == "hexahedron":
            cell_type, nodes = _collapse_solid(nodes)
        elif cell_type == "hexahedron20" and len(set(nodes[:8])) < 8:
            degenerate20 += 1
        conn = [index[v] for v in nodes]
        if cell_type not in blocks:
            blocks[cell_type] = {"conn": [], "id": [], "type": []}
            order.append(cell_type)
        b = blocks[cell_type]
        loc = (order.index(cell_type), len(b["id"]))
        locs.append(loc)
        element_cell.setdefault(ident, loc)
        b["conn"].append(conn)
        b["id"].append(ident)
        b["type"].append(etype)
    if skipped:
        listed = ", ".join(f"{t} ({c})" for t, c in sorted(skipped.items()))
        warn(
            f"{label}: elements of types meshio++ has no cell for were skipped: {listed}"
        )
    if degenerate20:
        warn(
            f"{label}: {degenerate20} 20-node brick(s) with repeated corner nodes "
            "kept as hexahedron20"
        )
    cells = [
        (
            t,
            np.array(blocks[t]["conn"], dtype=np.int64).reshape(
                len(blocks[t]["id"]), -1
            ),
        )
        for t in order
    ]
    mesh = Mesh(points, cells)
    if order:
        mesh.cell_data["marc:element"] = [
            np.array(blocks[t]["id"], dtype=np.int64) for t in order
        ]
        mesh.cell_data["marc:type"] = [
            np.array(blocks[t]["type"], dtype=np.int64) for t in order
        ]
    starts = np.cumsum([0] + [len(blocks[t]["id"]) for t in order])
    dims = [topological_dimension[t] for t in order]
    all_elements = {ident for ident, _, _ in elements}
    regions = []
    seen = set()
    for name, family, members in sets:
        if family in ("edge", "face"):
            # (cell or -1, Marc edge/face number): Marc numbers an element's
            # edges and faces its own way (Volume A), not mapped to facets.
            rows = []
            for m, number in members:
                loc = element_cell.get(m)
                rows.append(
                    [-1 if loc is None else int(starts[loc[0]]) + loc[1], number]
                )
            mesh.field_data[f"marc:{family}_set:{name}"] = np.array(
                rows, dtype=np.int64
            ).reshape(-1, 2)
            continue
        if family == "element":
            entries, dim = [], -1
            for m in members:
                loc = element_cell.get(m)
                if loc is None:
                    if m not in all_elements:
                        _fail(
                            label, f"element set '{name}' names undefined element {m}"
                        )
                    continue
                entries.append(int(starts[loc[0]]) + loc[1])
                dim = max(dim, dims[loc[0]])
            kind = "cell"
        else:
            entries, dim = [], -1
            for m in members:
                if m not in index:
                    _fail(label, f"node set '{name}' names undefined node {m}")
                entries.append(index[m])
            kind = "point"
        key = (kind, name)
        if key in seen:
            warn(f"{label}: a second {family} set '{name}' replaces the first")
        seen.add(key)
        regions = [r for r in regions if (r.kind, r.name) != key]
        regions.append(Region(name, kind, np.array(entries, dtype=np.int64), dim))
    mesh.regions = sorted(regions, key=lambda r: r.key)
    return mesh, locs


def _include_target(line):
    """An ``INCLUDE`` option line (the keyword at the start of the line, then
    the file name after a blank or comma): the file it names, else None."""
    if len(line) < 7 or line[:7].lower() != "include":
        return None
    if len(line) > 7 and line[7] not in " \t,":
        return None
    name = line[7:].lstrip(" \t,").strip()
    if len(name) >= 2 and name[0] in "\"'" and name[-1] == name[0]:
        name = name[1:-1]
    return name or None


def _deck_lines(text, directory, depth=0):
    """The deck's lines with every ``INCLUDE`` replaced by the lines of the
    file it names (relative to the including file), recursively."""
    if depth > 16:
        _fail(_DAT, "INCLUDE files nest more than 16 deep (a cycle?)")
    out = []
    for line in text.splitlines():
        target = _include_target(line.rstrip("\r"))
        if target is None:
            out.append(line)
            continue
        path = target if os.path.isabs(target) else os.path.join(directory, target)
        if not os.path.isfile(path):
            _fail(_DAT, f"INCLUDE names {path}, which does not exist")
        with open_file(path, "r") as f:
            out += _deck_lines(f.read(), os.path.dirname(path), depth + 1)
    return out


def read(filename):
    """Read an MSC Marc input deck (``.dat``), following its ``INCLUDE`` files."""
    with open_file(filename, "r") as f:
        text = f.read()
    if not is_marc_deck(text[:65536]):
        _fail(_DAT, "not a Marc input deck (no Marc parameter opens the file)")
    directory = (
        "." if is_buffer(filename, "r") else os.path.dirname(os.fspath(filename))
    )
    deck = _parse_deck(_deck_lines(text, directory))
    if not deck.nodes and not deck.elements:
        _fail(_DAT, "no COORDINATES or CONNECTIVITY found")
    known = {}
    sets = []
    for name, family, tokens, where in deck.sets:
        if family in ("edge", "face"):
            sets.append((name, family, tokens))  # members already read
            continue
        refs = {k: v for (f, k), v in known.items() if f == family}
        members = _expand(tokens, name, refs, _DAT)
        known[(family, name.lower())] = members
        sets.append((name, family, members))
    mesh, _ = _build(
        _DAT,
        deck.node_order,
        [deck.nodes[k] for k in deck.node_order],
        deck.elements,
        sets,
    )
    return mesh


# -- the formatted post file ---------------------------------------------------

_W = 13  # the post file's column width (i13, e13.6)

# Element post codes (Volume C, Table 3-3) that open a symmetric tensor written
# as six codes, components 11 22 33 12 23 31 (xx yy zz xy yz zx).
_TENSOR_NAMES = {
    301: "Total Strain",
    311: "Stress",
    321: "Plastic Strain",
    331: "Creep Strain",
    341: "Cauchy Stress",
    351: "Real Harmonic Stress",
    361: "Imaginary Harmonic Stress",
    371: "Thermal Strain",
    381: "Cracking Strain",
    391: "Stress in Preferred System",
    401: "Elastic Strain",
    411: "Global Stress",
    421: "Global Elastic Strain",
    431: "Global Plastic Strain",
    441: "Global Creep Strain",
    461: "Elastic Strain in Preferred System",
    541: "Phase Transformation Strain",
}
_SCALAR_NAMES = {
    7: "Equivalent Plastic Strain",
    8: "Equivalent Creep Strain",
    9: "Temperature",
    17: "Equivalent Von Mises Stress",
    18: "Mean Normal Stress",
    20: "Thickness",
    47: "Equivalent Cauchy Stress",
    48: "Strain Energy Density",
    127: "Equivalent Elastic Strain",
}


def _ints_of(line):
    body = line.rstrip("\r\n").rstrip()
    return [
        int(body[k : k + _W])
        for k in range(0, len(body), _W)
        if body[k : k + _W].strip()
    ]


def _reals_of(line, where):
    body = line.rstrip("\r\n").rstrip()
    return [
        _real(body[k : k + _W], _T19, where)
        for k in range(0, len(body), _W)
        if body[k : k + _W].strip()
    ]


class _Reader:
    """Reads a block's lines as Fortran records: each record starts a line."""

    def __init__(self, lines, start, end, label=_T19):
        self.lines = lines
        self.i = start
        self.end = end
        self.label = label

    def line(self):
        if self.i >= self.end:
            _fail(self.label, f"a block ends early (line {self.i + 1})")
        text = self.lines[self.i]
        self.i += 1
        return text

    def ints(self, count):
        out = []
        while len(out) < count:
            out += _ints_of(self.line())
        return out[:count] if count else out

    def reals(self, count):
        out = []
        while len(out) < count:
            where = f"line {self.i + 1}"
            out += _reals_of(self.line(), where)
        return out[:count]


def _blocks(lines):
    """``(kind, number, first body line, end line)`` of each block and each
    increment marker (``****``, ``----``, ``++++``)."""
    out = []
    i = 0
    n = len(lines)
    while i < n:
        s = lines[i].strip()
        if s.startswith("=beg="):
            try:
                number = int(s[5:10])
            except ValueError:
                _fail(_T19, f"bad block header '{s}' (line {i + 1})")
            j = i + 1
            while j < n and not lines[j].startswith("=end="):
                j += 1
            if j >= n:
                _fail(_T19, f"block {number} has no =end= (line {i + 1})")
            out.append(("block", number, i + 1, j))
            i = j + 1
            continue
        if s in ("****", "----", "++++"):
            out.append((s, 0, i, i))
        i += 1
    return out


class _Post:
    def __init__(self, filename):
        with open_file(filename, "r") as f:
            self.lines = f.read().splitlines()
        if not self.lines or not self.lines[0].startswith("=beg=501"):
            _fail(_T19, "not a Marc formatted post file (no =beg=501 title block)")
        self.blocks = _blocks(self.lines)
        # The model blocks come before the first increment.
        self.increments = []  # the block list of each increment
        current = None
        self.model = []
        for b in self.blocks:
            if b[0] == "****":
                current = []
                self.increments.append(current)
            elif b[0] in ("----", "++++"):
                current = None
            elif current is None:
                self.model.append(b)
            else:
                current.append(b)

    def reader(self, block):
        return _Reader(self.lines, block[2], block[3])

    def header(self, blocks=None):
        """The model header (block 502) and element post codes (506) of
        ``blocks``, else of the model before the first increment."""
        lm = [0] * 30
        codes = []
        for b in self.model if blocks is None else blocks:
            family = b[1] // 100
            if family == 502:
                values = self.reader(b).ints(30)
                lm[: len(values)] = values
            elif family == 506:
                r = self.reader(b)
                for _ in range(lm[0]):
                    line = r.line()
                    codes.append((int(line[:_W]), line[_W : _W + 24].strip()))
        return lm, codes

    def mesh(self, lm, blocks=None):
        numnp, numel, ncrd, nnodmx, postrv = lm[1], lm[2], lm[8], lm[9], lm[13]
        node_ids, coords, elements, sets = [], [], [], []
        for b in self.model if blocks is None else blocks:
            family = b[1] // 100
            r = self.reader(b)
            if family == 507:
                for _ in range(numel):
                    record = r.ints(3 + nnodmx)
                    ident, etype, nnod = record[:3]
                    elements.append((ident, etype, record[3 : 3 + nnod]))
            elif family == 508:
                for _ in range(numnp):
                    line = r.line()
                    where = f"line {r.i}"
                    ident = int(line[:_W])
                    values = _reals_of(line[_W:].ljust(5 * _W), where)[:5]
                    if ncrd > 5:
                        values += r.reals(ncrd - 5)
                    xyz = values[:ncrd] + [0.0] * (3 - min(ncrd, 3))
                    node_ids.append(ident)
                    coords.append(xyz[:3])
            elif family == 513:
                if b[1] == 51301 or postrv > 10:
                    count = r.ints(1)[0]
                    width = 32
                else:
                    count, width = lm[15], 12
                for _ in range(count):
                    name = r.line()[:width].strip()
                    isetn, isett = r.ints(2)
                    members = r.ints(isetn) if isetn else []
                    if isett in (12, 13, 18, 19):
                        numbers = r.ints(isetn) if isetn else []
                        family = "edge" if isett == 12 else "face"
                        sets.append((name, family, list(zip(members, numbers))))
                    elif isett == 0:
                        sets.append((name, "element", members))
                    elif isett == 1:
                        sets.append((name, "node", members))
                    else:
                        warn(f"{_T19}: set '{name}' of type {isett} is not read")
        return node_ids, coords, elements, sets

    def increment_info(self, blocks):
        """``(time, inc, subinc, jantyp, ihresp, newmo)`` of an increment."""
        lm = [0] * 12
        xlm = [0.0] * 6
        for b in blocks:
            family = b[1] // 100
            if family == 517:
                values = self.reader(b).ints(12)
                lm[: len(values)] = values
            elif b[1] == 51800:
                xlm = self.reader(b).reals(6)
            elif b[1] == 51801:
                r = self.reader(b)
                nw = r.ints(1)[0]
                xlm = r.reals(nw) + [0.0] * 6
        newmo, inc, incsub, jantyp, _, _, ihresp = lm[:7]
        time = xlm[1] if ihresp in (1, 2, 3, 4) else xlm[0]
        return time, inc, incsub, jantyp, ihresp, newmo

    def times(self):
        return [self.increment_info(inc)[0] for inc in self.increments]


def _element_arrays(codes, values, nstres):
    """``{name: (elements, points[, 6])}`` from ``values`` (elements x IPs x
    codes): six consecutive codes ``c .. c + 5`` of a tensor become one
    ``xx yy zz xy yz zx`` array; a layer (code + 1000 x layer) is ``@layer<n>``.
    With one integration point the point axis is dropped."""
    out = {}
    k = 0
    while k < len(codes):
        code, label = codes[k]
        layer, base = divmod(code, 1000)
        suffix = f"@layer{layer}" if layer else ""
        if (
            base in _TENSOR_NAMES
            and k + 5 < len(codes)
            and [c for c, _ in codes[k : k + 6]] == list(range(code, code + 6))
        ):
            name = (label or _TENSOR_NAMES[base]) + suffix
            out[name] = values[:, :, k : k + 6]
            k += 6
            continue
        name = (label or _SCALAR_NAMES.get(base, f"post code {base}")) + suffix
        out[name] = values[:, :, k]
        k += 1
    if nstres == 1:
        out = {name: a[:, 0] for name, a in out.items()}
    return out


def read_t19(filename, points_only=False, arrays=None, time_step=0):
    """Read one increment of an MSC Marc formatted post file (``.t19``)."""
    post = _Post(filename)
    n = len(post.increments)
    index = time_step + n if time_step < 0 else time_step
    if n and not 0 <= index < n:
        raise ReadError(
            f"time step {time_step} is out of range: the file has {n} step(s)"
        )
    # An increment that remeshes (newmo, BLOCK 517) repeats the model blocks
    # 502 to 514 (BLOCK 519): a step's mesh is the latest model at or before it.
    remeshed = None
    for k in range(index + 1 if n else 0):
        if post.increment_info(post.increments[k])[5]:
            remeshed = [b for b in post.increments[k] if 502 <= b[1] // 100 <= 514]
    lm, codes = post.header()
    if remeshed and any(b[1] // 100 == 502 for b in remeshed):
        # A remeshed model repeats the header; the post codes stay the
        # first's when it does not repeat them.
        lm, codes2 = post.header(remeshed)
        codes = codes2 or codes
    npost, numnp, numel, nstres = lm[0], lm[1], lm[2], max(lm[4], 1)
    node_ids, coords, elements, sets = post.mesh(lm, remeshed or None)
    mesh, locs = _build(_T19, node_ids, coords, elements, sets)
    mesh.time_values = post.times()
    if not n:
        if time_step not in (0, -1):
            _fail(
                _T19,
                f"time step {time_step} is out of range: the file has no increments",
            )
        return mesh
    blocks = post.increments[index]
    time, inc, incsub, jantyp, _, _ = post.increment_info(blocks)
    mesh.field_data["meshio:time"] = np.array([time], dtype=np.float64)
    mesh.field_data["marc:increment"] = np.array([inc], dtype=np.int64)
    mesh.field_data["marc:subincrement"] = np.array([incsub], dtype=np.int64)
    if points_only:
        return mesh
    wanted = None if arrays is None else set(arrays)
    for b in blocks:
        family = b[1] // 100
        r = post.reader(b)
        if family == 524:
            nnqnod, _ = r.ints(2)
            for _ in range(nnqnod):
                name = r.line()[:48].strip()
                ivec = r.ints(12)
                if ivec[6] != -1:
                    continue
                ncomp = max(ivec[3], 1)
                real = np.array(r.reals(numnp * ncomp)).reshape(numnp, ncomp)
                imag = None
                if ivec[5] in (4, 5):
                    imag = np.array(r.reals(numnp * ncomp)).reshape(numnp, ncomp)
                for key, data in ((name, real), (name + "@imag", imag)):
                    if data is None or (wanted is not None and key not in wanted):
                        continue
                    mesh.point_data[key] = data[:, 0] if ncomp == 1 else data
        elif family == 523 and jantyp > 100 and npost and numel:
            values = np.empty((numel, nstres, npost))
            for e in range(numel):
                for p in range(nstres):
                    values[e, p] = r.reals(npost)
            for name, data in _element_arrays(codes, values, nstres).items():
                if wanted is not None and name not in wanted:
                    continue
                if nstres > 1:
                    # Flattened point-major, so that any writer holds it; its
                    # (points, components) is the layout.
                    comps = data.shape[2] if data.ndim == 3 else 1
                    data = data.reshape(len(data), nstres * comps)
                    mesh.field_data["marc:layout:" + name] = np.array(
                        [nstres, comps], dtype=np.int64
                    )
                per_block = []
                for b_index, cell_block in enumerate(mesh.cells):
                    shape = (len(cell_block.data),) + data.shape[1:]
                    per_block.append(np.full(shape, np.nan))
                for e, loc in enumerate(locs):
                    if loc is not None:
                        per_block[loc[0]][loc[1]] = data[e]
                mesh.cell_data[name] = per_block
    return mesh


def time_values(filename):
    """The time (a frequency or buckling factor in a modal, harmonic or buckling
    analysis) of every increment."""
    return _Post(filename).times()

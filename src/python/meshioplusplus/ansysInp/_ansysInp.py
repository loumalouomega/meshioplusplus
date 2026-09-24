"""Ansys MAPDL coded database (``.cdb``): the pure-Python reference engine.

The twin of ``src/cpp/src/formats/ansysinp.cpp``: both read the same meshes and
write the same bytes. Blocks are sliced by their own Fortran format lines; an
element becomes a cell by its routine's category (pymapdl-reader's table), with
degenerate bricks and shells resolved by repeated nodes. See
``doc/formats/ansysinp.md``.
"""

from __future__ import annotations

import numpy as np

from .. import _provenance
from .._common import warn
from .._exceptions import ReadError, WriteError
from .._files import open_file
from .._mesh import CellBlock, Mesh
from .._regions import Region
from ..lsdyna._cards import parse_fortran_format, split_fixed

# -- element categories (ansysinp.cpp's ans_category) ----------------------------
_POINT, _LINE, _SHELL, _BRICK, _TET, _LINEAR_LINE = range(1, 7)
_CATEGORY = {}
for _n in (7, 21, 71, 175):
    _CATEGORY[_n] = _POINT
for _n in (
    1, 3, 4, 8, 10, 11, 12, 14, 16, 17, 18, 20, 23, 24, 31, 32, 33, 34, 37, 38, 39,
    40, 44, 59, 60, 61, 66, 68, 116, 126, 129, 151, 153, 156, 161, 169, 171, 172,
    176, 177, 178, 180, 189, 208, 209, 250, 251, 280, 288, 289, 290,
):  # fmt: skip
    _CATEGORY[_n] = _LINE
for _n in (
    2, 13, 22, 25, 28, 29, 35, 41, 42, 43, 51, 53, 54, 55, 57, 63, 67, 75, 77, 78,
    79, 81, 82, 83, 88, 91, 93, 99, 106, 115, 118, 121, 130, 131, 132, 136, 143,
    152, 154, 155, 157, 163, 170, 173, 174, 181, 182, 183, 212, 213, 218, 219, 222,
    223, 230, 233, 238, 252, 281, 282, 283, 292, 293,
):  # fmt: skip
    _CATEGORY[_n] = _SHELL
for _n in (
    5, 30, 45, 46, 62, 64, 65, 69, 70, 80, 89, 90, 95, 96, 97, 100, 101, 102, 103,
    104, 105, 107, 108, 117, 120, 122, 164, 185, 186, 190, 192, 215, 220, 226, 231,
    236, 239, 272, 273, 278, 279,
):  # fmt: skip
    _CATEGORY[_n] = _BRICK
for _n in (87, 92, 98, 119, 123, 140, 168, 187, 221, 227, 232, 237, 240, 285, 291):
    _CATEGORY[_n] = _TET
for _n in (188, 214, 216, 217):
    _CATEGORY[_n] = _LINEAR_LINE


def _mesh200_category(keyopt1):
    if 0 <= keyopt1 <= 3:
        return _LINE
    if 4 <= keyopt1 <= 7:
        return _SHELL
    if keyopt1 in (8, 9):
        return _TET
    if keyopt1 in (10, 11):
        return _BRICK
    return 0


def _resolve(category, nodes):
    """(meshio++ type, slots) of one element row, or (None, None)."""
    n = len(nodes)

    def at(k):
        return nodes[k] if k < n else 0

    if category == _POINT and n >= 1:
        return "vertex", [0]
    if category == _LINE:
        if n >= 3 and at(2) > 0:
            return "line3", [0, 1, 2]
        if n >= 2:
            return "line", [0, 1]
    if category == _LINEAR_LINE and n >= 2:
        return "line", [0, 1]
    if category == _SHELL:
        if n == 3:
            return "triangle", [0, 1, 2]
        if n == 6:
            return "triangle6", [0, 1, 2, 3, 4, 5]
        if n > 5:  # 8-node (5 is a quad plus an orientation node); absent
            # trailing midsides are missing ones
            if at(2) == at(3):
                return "triangle6", [0, 1, 2, 4, 5, 7]
            if at(6) == at(7):  # e.g. a linear contact face with no midsides
                return "quad", [0, 1, 2, 3]
            return "quad8", list(range(8))
        if n >= 4:
            if at(2) == at(3):
                return "triangle", [0, 1, 2]
            return "quad", [0, 1, 2, 3]
    if category == _TET:
        if n > 4:
            return "tetra10", list(range(10))
        if n >= 4:
            return "tetra", [0, 1, 2, 3]
    if category == _BRICK and n >= 8:
        quad = n > 8
        if at(6) != at(7):
            return (
                ("hexahedron20", list(range(20)))
                if quad
                else (
                    "hexahedron",
                    list(range(8)),
                )
            )
        if at(5) != at(6):
            if quad:
                return "wedge15", [0, 1, 2, 4, 5, 6, 8, 9, 11, 12, 13, 15, 16, 17, 18]
            return "wedge", [0, 1, 2, 4, 5, 6]
        if at(2) != at(3):
            if quad:
                return "pyramid13", [0, 1, 2, 3, 4, 8, 9, 10, 11, 16, 17, 18, 19]
            return "pyramid", [0, 1, 2, 3, 4]
        if quad:
            return "tetra10", [0, 1, 2, 4, 8, 9, 11, 16, 17, 18]
        return "tetra", [0, 1, 2, 4]
    return None, None


_MIDSIDE_EDGES = {
    "line3": [(0, 1)],
    "triangle6": [(0, 1), (1, 2), (2, 0)],
    "quad8": [(0, 1), (1, 2), (2, 3), (3, 0)],
    "tetra10": [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)],
    "pyramid13": [(0, 1), (1, 2), (2, 3), (3, 0), (0, 4), (1, 4), (2, 4), (3, 4)],
    "wedge15": [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (0, 3), (1, 4), (2, 5)],
    "hexahedron20": [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4)]
    + [(0, 4), (1, 5), (2, 6), (3, 7)],
}
_CORNERS = {
    "vertex": 1,
    "line": 2,
    "triangle": 3,
    "quad": 4,
    "tetra": 4,
    "pyramid": 5,
    "wedge": 6,
    "hexahedron": 8,
}


def _num_corners(cell_type):
    return _CORNERS[cell_type.rstrip("0123456789")]


def _fail(line, what):
    raise ReadError(f"Ansys .cdb: line {line + 1}: {what}")


def _commas(line):
    return [p.strip(" \t\r\n") for p in line.split(",")]


def _int_or_none(text):
    if not text:
        return None
    try:
        return int(float(text))
    except (ValueError, OverflowError):
        return None


def _routine(text):
    v = _int_or_none(text)
    if v is not None:
        return v
    k = len(text)
    while k > 0 and text[k - 1].isdigit():
        k -= 1
    if k == len(text):
        return -1
    return int(text[k:])


def _field_int(fields, k, line):
    if k >= len(fields) or not fields[k]:
        return 0
    v = _int_or_none(fields[k])
    if v is None:
        _fail(line, f"bad integer '{fields[k]}'")
    return v


def _field_real(fields, k, line):
    if k >= len(fields) or not fields[k]:
        return 0.0
    try:
        return float(fields[k])
    except ValueError:
        _fail(line, f"bad number '{fields[k]}'")


def _is_terminator(line):
    s = line.strip(" \t\r\n")
    if s == "-1":
        return True
    up = s.upper()
    return up.startswith("N,") and "LOC" in up


def _is_command(line):
    """A command line (``FINISH``, ``CMBLOCK,...``) rather than a block's data line."""
    s = line.lstrip(" \t")
    return bool(s) and (s[0].isascii() and s[0].isalpha() or s[0] in "/!")


def _is_keyopt(upper):
    """``KEYOPT,`` or any abbreviation of it MAPDL accepts (``KEYO,``, ``KEYOP,``)."""
    comma = upper.find(",")
    return 4 <= comma <= 6 and "KEYOPT"[:comma] == upper[:comma]


def _format(lines, k):
    if k >= len(lines):
        _fail(k, "a block header is not followed by its format line")
    try:
        return parse_fortran_format(lines[k].strip(" \t\r\n"))
    except ReadError as exc:
        _fail(k, str(exc))


def _parse(lines):
    deck = {
        "routine": {},
        "keyopt": {},
        "node_ids": [],
        "coords": [],
        "rotations": False,
        "elements": [],
        "components": [],
        "nonsolid": 0,
    }
    saw_block = False
    i, n = 0, len(lines)
    while i < n:
        # A command line, less any trailing ``!`` comment.
        line = lines[i].split("!", 1)[0].strip(" \t\r\n")
        up = line.upper()
        if up.startswith("ET,"):
            p = _commas(line)
            if len(p) >= 3:
                slot, routine = _int_or_none(p[1]), _routine(p[2].upper())
                if slot is not None and routine > 0:
                    deck["routine"][slot] = routine
            i += 1
        elif _is_keyopt(up):
            p = _commas(line)
            if len(p) >= 4:
                slot, k, v = (_int_or_none(x) for x in p[1:4])
                if None not in (slot, k, v):
                    deck["keyopt"].setdefault(slot, {})[k] = v
            i += 1
        elif up.startswith("ETBLOCK"):
            saw_block = True
            fields = _format(lines, i + 1)
            i += 2
            while i < n and not _is_terminator(lines[i]):
                f = split_fixed(lines[i], fields)
                if len(f) >= 2:
                    slot = _field_int(f, 0, i)
                    deck["routine"][slot] = _field_int(f, 1, i)
                    for k in range(2, min(len(f), 20)):
                        v = _int_or_none(f[k])
                        if v:
                            deck["keyopt"].setdefault(slot, {})[k - 1] = v
                i += 1
            i += 1
        elif up.startswith("NBLOCK"):
            saw_block = True
            fields = _format(lines, i + 1)
            n_int = 0
            while n_int < len(fields) and fields[n_int][0] == "i":
                n_int += 1
            i += 2
            while i < n and not _is_terminator(lines[i]):
                f = split_fixed(lines[i], fields)
                if not f or not f[0]:
                    i += 1
                    continue
                deck["node_ids"].append(_field_int(f, 0, i))
                deck["coords"].append([_field_real(f, n_int + d, i) for d in range(3)])
                for d in range(3, len(f) - n_int):
                    if _field_real(f, n_int + d, i) != 0.0:
                        deck["rotations"] = True
                i += 1
            i += 1
        elif up.startswith("EBLOCK"):
            saw_block = True
            header = _commas(up)
            solid = len(header) > 2 and header[2] == "SOLID"
            fields = _format(lines, i + 1)
            i += 2
            if not solid:
                deck["nonsolid"] += 1
                while i < n and not _is_terminator(lines[i]):
                    i += 1
                i += 1
                continue
            while i < n and not _is_terminator(lines[i]):
                f = split_fixed(lines[i], fields)
                if len(f) < 11:
                    i += 1
                    continue
                count = max(_field_int(f, 8, i), 0)
                element = {
                    "mat": _field_int(f, 0, i),
                    "slot": _field_int(f, 1, i),
                    "real": _field_int(f, 2, i),
                    "secnum": _field_int(f, 3, i),
                    "id": _field_int(f, 10, i),
                    "nodes": [],
                }
                nodes = element["nodes"]
                for k in range(11, len(f)):
                    if len(nodes) >= count:
                        break
                    nodes.append(_field_int(f, k, i))
                i += 1
                while len(nodes) < count and i < n and not _is_terminator(lines[i]):
                    more = split_fixed(lines[i], fields)
                    for k in range(len(more)):
                        if len(nodes) >= count:
                            break
                        nodes.append(_field_int(more, k, i))
                    i += 1
                deck["elements"].append(element)
            i += 1
        elif up.startswith("CMBLOCK"):
            saw_block = True
            header = _commas(line)
            if len(header) < 3:
                _fail(i, "a CMBLOCK needs a name and an entity type")
            name, entity = header[1], header[2].upper()
            known = entity == "NODE" or entity.startswith("ELEM")
            count = (_int_or_none(header[3]) or 0) if len(header) > 3 else 0
            fields = _format(lines, i + 1)
            i += 2
            raw = []
            # A short block (a header count too large) ends at the next command.
            while i < n and len(raw) < count and not _is_command(lines[i]):
                f = split_fixed(lines[i], fields)
                if not f:
                    break
                for k in range(len(f)):
                    if len(raw) >= count:
                        break
                    if f[k]:
                        raw.append(_field_int(f, k, i))
                i += 1
            ids = []
            for v in raw:
                if v >= 0:
                    ids.append(v)
                    continue
                if not ids:
                    _fail(i, f"CMBLOCK '{name}' opens with a range end")
                ids.extend(range(ids[-1] + 1, -v + 1))
            if known:
                deck["components"].append((name, entity == "NODE", ids))
        else:
            i += 1
    if not saw_block:
        raise ReadError("Ansys .cdb: no NBLOCK, EBLOCK or CMBLOCK found")
    return deck


def read(filename, lenient=False):
    with open_file(filename, "r") as f:
        lines = [line.rstrip("\r\n") for line in f]
    return _read_lines(lines, lenient=lenient)


def _read_lines(lines, lenient=False):
    deck = _parse(lines)
    if deck["nonsolid"]:
        warn(
            f"Ansys .cdb: {deck['nonsolid']} non-solid EBLOCK(s) (MAPDL writes only "
            "the SOLID layout) were skipped"
        )
    if deck["rotations"]:
        warn("Ansys .cdb: nodal rotation angles are not kept")
    return _build(deck, lenient, "Ansys .cdb")[0]


def _build(deck, lenient, label, locations=False):
    """The mesh of a parsed deck (shared with the ``.rst`` reader): cells by
    element category, missing midsides created, ``ansys:*`` cell data and
    component regions. Returns it with the node-number -> point-index map, and
    with ``locations`` also a list parallel to ``deck["elements"]``: each
    element's ``(block, row, slots)`` -- ``slots[j]`` is the element node that
    became the cell's node ``j`` -- or ``None`` for an element with no cell."""
    coords = list(deck["coords"])
    node_index = {}
    for k, ident in enumerate(deck["node_ids"]):
        node_index.setdefault(ident, k)
    midsides = {}
    blocks = {}
    order = []
    element_loc = {}
    locs = [None] * len(deck["elements"])
    skipped = {}
    for pos, e in enumerate(deck["elements"]):
        if e["slot"] not in deck["routine"]:
            raise ReadError(
                f"{label}: element {e['id']} uses element type {e['slot']}, which "
                "no ET or ETBLOCK defines"
            )
        routine = deck["routine"][e["slot"]]
        category = _CATEGORY.get(routine, 0)
        if routine == 200:
            category = _mesh200_category(deck["keyopt"].get(e["slot"], {}).get(1, 0))
        cell_type, slots = _resolve(category, e["nodes"])
        if cell_type is not None:
            # A corner node 0 (TARGE170's pilot and line shapes ...) has no cell.
            nodes = e["nodes"]
            corner_slots = slots[: _num_corners(cell_type)]
            if any(k >= len(nodes) or nodes[k] == 0 for k in corner_slots):
                cell_type = None
        if cell_type is None:
            if not lenient:
                raise ReadError(
                    f"{label}: element {e['id']} (element type {routine}, "
                    f"{len(e['nodes'])} nodes) has no meshio++ cell type; read with "
                    "lenient to skip it"
                )
            skipped[routine] = skipped.get(routine, 0) + 1
            continue
        if cell_type not in blocks:
            order.append(cell_type)
            blocks[cell_type] = {
                k: [] for k in ("conn", "routine", "slot", "mat", "real", "secnum")
            }
        b = blocks[cell_type]
        corners = _num_corners(cell_type)
        row = []
        for slot in slots:
            # A row cut short omits trailing midside nodes: they read as node 0.
            ident = e["nodes"][slot] if slot < len(e["nodes"]) else 0
            if ident == 0 and len(row) >= corners:
                row.append(-1)
                continue
            if ident not in node_index:
                raise ReadError(
                    f"{label}: element {e['id']} names undefined node {ident}"
                )
            row.append(node_index[ident])
        edges = _MIDSIDE_EDGES.get(cell_type, [])
        for k in range(corners, len(row)):
            if row[k] >= 0:
                continue
            a, c = edges[k - corners]
            p, q = row[a], row[c]
            key = (min(p, q), max(p, q))
            if key not in midsides:
                midsides[key] = len(coords)
                coords.append([0.5 * (coords[p][d] + coords[q][d]) for d in range(3)])
            row[k] = midsides[key]
        element_loc[e["id"]] = (order.index(cell_type), len(b["slot"]))
        locs[pos] = (order.index(cell_type), len(b["slot"]), tuple(slots))
        b["conn"].append(row)
        b["routine"].append(routine)
        b["slot"].append(e["slot"])
        b["mat"].append(e["mat"])
        b["real"].append(e["real"])
        b["secnum"].append(e["secnum"])
    for routine in sorted(skipped):
        warn(
            f"{label}: {skipped[routine]} element(s) of type {routine} skipped "
            "(no meshio++ cell type)"
        )
    if midsides:
        warn(
            f"{label}: {len(midsides)} missing midside node(s) placed at their "
            "edge midpoints"
        )

    points = np.array(coords, dtype=np.float64).reshape(-1, 3)
    cells = [
        CellBlock(
            t,
            np.array(blocks[t]["conn"], dtype=np.int64).reshape(
                len(blocks[t]["slot"]), -1
            ),
        )
        for t in order
    ]
    mesh = Mesh(points, cells)
    if order:
        for key, name in (
            ("routine", "ansys:element"),
            ("slot", "ansys:type"),
            ("mat", "ansys:mat"),
            ("real", "ansys:real"),
            ("secnum", "ansys:secnum"),
        ):
            mesh.cell_data[name] = [
                np.array(blocks[t][key], dtype=np.int64) for t in order
            ]
    bases = np.cumsum([0] + [len(blocks[t]["slot"]) for t in order])
    regions = []
    for name, is_nodes, ids in deck["components"]:
        if is_nodes:
            entries = [node_index[i] for i in ids if i in node_index]
        else:
            entries = [
                int(bases[element_loc[i][0]]) + element_loc[i][1]
                for i in ids
                if i in element_loc
            ]
        regions.append(
            Region(
                name, "point" if is_nodes else "cell", np.array(entries, dtype=np.int64)
            )
        )
    # C++ keeps one region per (kind, name): the last of two same-named ones wins.
    unique = {r.key: r for r in regions}
    mesh.regions = sorted(unique.values(), key=lambda r: r.key)
    if locations:
        return mesh, node_index, locs
    return mesh, node_index


# -- writing ---------------------------------------------------------------------

_LAYOUTS = {
    _BRICK: {
        "hexahedron": list(range(8)),
        "hexahedron20": list(range(20)),
        "wedge": [0, 1, 2, 2, 3, 4, 5, 5],
        "wedge15": [0, 1, 2, 2, 3, 4, 5, 5, 6, 7, 2, 8, 9, 10, 5, 11, 12, 13, 14, 14],
        "pyramid": [0, 1, 2, 3, 4, 4, 4, 4],
        "pyramid13": [0, 1, 2, 3, 4, 4, 4, 4, 5, 6, 7, 8, 4, 4, 4, 4, 9, 10, 11, 12],
        "tetra": [0, 1, 2, 2, 3, 3, 3, 3],
        "tetra10": [0, 1, 2, 2, 3, 3, 3, 3, 4, 5, 2, 6, 3, 3, 3, 3, 7, 8, 9, 9],
    },
    _TET: {"tetra": [0, 1, 2, 3], "tetra10": list(range(10))},
    _SHELL: {
        "quad": [0, 1, 2, 3],
        "quad8": list(range(8)),
        "triangle": [0, 1, 2, 2],
        "triangle6": [0, 1, 2, 2, 3, 4, 2, 5],
    },
    _LINE: {"line": [0, 1], "line3": [0, 1, 2]},
    _LINEAR_LINE: {"line": [0, 1]},
    _POINT: {"vertex": [0]},
}
_DEFAULT_ROUTINE = {
    "vertex": 21,
    "line": 188,
    "line3": 189,
    "triangle": 181,
    "triangle6": 281,
    "quad": 181,
    "quad8": 281,
    "tetra": 285,
    "tetra10": 187,
    "pyramid": 185,
    "pyramid13": 186,
    "wedge": 185,
    "wedge15": 186,
    "hexahedron": 185,
    "hexahedron20": 186,
}


def _layout(cell_type, routine):
    return _LAYOUTS.get(_CATEGORY.get(routine, 0), {}).get(cell_type)


def _i(value, width):
    return f"{int(value):>{width}d}"


def _cell_int(mesh, name, block, row, default):
    if name not in mesh.cell_data:
        return default
    arr = np.asarray(mesh.cell_data[name][block])
    if len(arr) <= row:
        return default
    return int(arr[row])


def write(filename, mesh):
    blocks = mesh.cells
    slot_routine = {}
    routine_slot = {}
    cell_etype = []
    dropped_etype = False
    for b, block in enumerate(blocks):
        fallback = _DEFAULT_ROUTINE.get(block.type, -1)
        if isinstance(block.data, list) or fallback < 0:
            raise WriteError(
                f"Ansys .cdb writer: cell type '{block.type}' has no element type"
            )
        per = []
        for r in range(len(block.data)):
            routine = _cell_int(mesh, "ansys:element", b, r, fallback)
            slot = _cell_int(mesh, "ansys:type", b, r, 0)
            if _layout(block.type, routine) is None or (
                slot > 0 and slot in slot_routine and slot_routine[slot] != routine
            ):
                dropped_etype = dropped_etype or routine != fallback
                routine, slot = fallback, 0
            if slot > 0:
                slot_routine.setdefault(slot, routine)
            per.append([routine, slot])
        cell_etype.append(per)
    for per in cell_etype:
        for pair in per:
            if pair[1] == 0:
                routine = pair[0]
                if routine not in routine_slot:
                    fresh = 1
                    while fresh in slot_routine:
                        fresh += 1
                    slot_routine[fresh] = routine
                    routine_slot[routine] = fresh
                pair[1] = routine_slot[routine]
    if dropped_etype:
        warn(
            "Ansys .cdb writer: some ansys:element values do not fit their cells' "
            "shape; the default element type was written instead"
        )

    bases = [0]
    for block in blocks:
        bases.append(bases[-1] + len(block.data))
    regions = sorted(getattr(mesh, "regions", []) or [], key=lambda r: r.key)
    node_comps, elem_comps = [], []
    seen = set()
    side_regions = 0
    for reg in regions:
        if reg.kind == "side":
            side_regions += 1
            continue
        is_nodes = reg.kind == "point"
        if (is_nodes, reg.name) in seen:
            continue
        seen.add((is_nodes, reg.name))
        ids = [int(v) + 1 for v in np.asarray(reg.entries).ravel()]
        (node_comps if is_nodes else elem_comps).append((reg.name, ids))
    if side_regions:
        warn(
            f"Ansys .cdb writer: {side_regions} side region(s) have no component "
            "equivalent and were dropped"
        )
        _provenance.note("regions-dropped", f"{side_regions} side region(s) dropped")

    out = [_provenance.render_lines(_provenance.SlotTier.BLOCK, "! "), "/PREP7\n"]
    for slot in sorted(slot_routine):
        out.append(f"ET,{slot},{slot_routine[slot]}\n")
    points = np.asarray(mesh.points, dtype=np.float64)
    npts = len(points)
    pdim = points.shape[1] if points.ndim == 2 else 0
    out.append(f"NBLOCK,6,SOLID,{_i(npts, 9)},{_i(npts, 9)}\n(3i9,6e21.13e3)\n")
    for p, row in enumerate(points.tolist()):
        xyz = "".join("%21.13E" % (row[d] if d < pdim else 0.0) for d in range(3))
        out.append(f"{_i(p + 1, 9)}{_i(0, 9)}{_i(0, 9)}{xyz}\n")
    out.append("N,R5.3,LOC,       -1,\n")
    n_cells = bases[-1]
    out.append(f"EBLOCK,19,SOLID,{_i(n_cells, 10)},{_i(n_cells, 10)}\n(19i10)\n")
    element = 0
    for b, block in enumerate(blocks):
        for r, row in enumerate(np.asarray(block.data).tolist()):
            routine, slot = cell_etype[b][r]
            nodes = [row[k] + 1 for k in _layout(block.type, routine)]
            element += 1
            head = [
                _cell_int(mesh, "ansys:mat", b, r, 1),
                slot,
                _cell_int(mesh, "ansys:real", b, r, 1),
                _cell_int(mesh, "ansys:secnum", b, r, 1),
                0,
                0,
                0,
                0,
                len(nodes),
                0,
                element,
            ]
            text = "".join(_i(v, 10) for v in head)
            for c, v in enumerate(nodes):
                if c == 8:
                    text += "\n"
                text += _i(v, 10)
            out.append(text + "\n")
    out.append("        -1\n")

    def component(name, entity, ids):
        ids = sorted(set(ids))
        packed = []
        k = 0
        while k < len(ids):
            j = k
            while j + 1 < len(ids) and ids[j + 1] == ids[j] + 1:
                j += 1
            packed.append(ids[k])
            if j > k:
                packed.append(-ids[j])
            k = j + 1
        out.append(f"CMBLOCK,{name},{entity},{_i(len(packed), 8)}\n(8i10)\n")
        for c, v in enumerate(packed):
            out.append(_i(v, 10))
            if c % 8 == 7 or c + 1 == len(packed):
                out.append("\n")

    for name, ids in elem_comps:
        component(name, "ELEM", ids)
    for name, ids in node_comps:
        component(name, "NODE", ids)
    out.append("FINISH\n")
    with open_file(filename, "w") as f:
        f.write("".join(out))

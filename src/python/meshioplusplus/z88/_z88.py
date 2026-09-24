"""I/O for Z88 / Z88Aurora structure files (``z88i1.txt``) and their results.

The pure-Python twin of ``src/cpp/src/formats/z88.cpp``: both engines read the
same meshes and write the same bytes.

``z88i1.txt``: a header whose first three integers are the dimension, node and
element counts (Z88OS v15: ``ndim nnodes nelem ndof kflag``; Z88 <= V13 and
Z88Aurora V1 add more flags and material lines after the elements); one line per
node, ``id ndof x y [z]``; two lines per element, ``id type`` then its nodes.
The rest of the deck next to it is attached: ``z88i2.txt`` constraints as
``point_data["z88:bc:u"]``/``["z88:bc:f"]``, ``z88mat.txt`` materials, ``z88elp.txt``
element parameters and ``z88int.txt`` integration orders as ``z88:`` cell data;
with ``results``, ``z88o2.txt`` displacements as ``point_data["U"]``, ``z88o4.txt``
nodal forces as ``point_data["F"]`` and ``z88o3.txt`` stresses as cell data.
``write`` writes the same files back from those arrays.
"""

import math
import os

import numpy as np

from .. import _provenance
from .._common import warn
from .._exceptions import ReadError, WriteError
from .._files import is_buffer, open_file
from .._mesh import Mesh, topological_dimension
from .._node_order import node_order
from .._regions import Region

__all__ = ["read", "write"]

# code -> (nodes, meshio++ type, nodes kept, degrees of freedom per node)
_TYPES = {
    1: (8, "hexahedron", 8, 3),
    2: (2, "line", 2, 6),
    3: (6, "triangle6", 6, 2),
    4: (2, "line", 2, 3),
    5: (2, "line", 2, 6),
    6: (3, "triangle", 3, 2),
    7: (8, "quad8", 8, 2),
    8: (8, "quad8", 8, 2),
    9: (2, "line", 2, 2),
    10: (20, "hexahedron20", 20, 3),
    11: (12, "quad", 4, 2),
    12: (12, "quad", 4, 2),
    13: (2, "line", 2, 3),
    14: (6, "triangle6", 6, 2),
    15: (6, "triangle6", 6, 2),
    16: (10, "tetra10", 10, 3),
    17: (4, "tetra", 4, 3),
    18: (6, "triangle6", 6, 3),
    19: (16, "VTK_LAGRANGE_QUADRILATERAL", 16, 3),
    20: (8, "quad8", 8, 3),
    21: (16, "hexahedron", 8, 3),
    22: (12, "wedge", 6, 3),
    23: (8, "quad8", 8, 6),
    24: (6, "triangle6", 6, 6),
    25: (2, "line", 2, 6),
}

# Corner slots of the types that keep only their corners, when not the leading
# nodes: 21/22 are two quad8/tri6 layers, the upper one first (the wedge takes
# the lower triangle first, as meshio++'s wedge does; the hexahedron goes
# through the Z88 hexahedron table). Type 19 is kept whole as a cubic
# Lagrange quad: Z88 numbers its 4x4 lattice with r slow and s fast (node
# 4i + j + 1 at the i-th r and the j-th s, as its shape functions place them;
# corners 1, 13, 16, 4 counter-clockwise), listed here in VTK's order (corners,
# the edges r, s, r, s, then the interior).
_CORNERS = {
    19: [0, 12, 15, 3, 4, 8, 13, 14, 7, 11, 1, 2, 5, 9, 6, 10],
    21: [0, 1, 2, 3, 8, 9, 10, 11],
    22: [6, 7, 8, 0, 1, 2],
}

_SOLID = (1, 10, 16, 17)
# Types that only exist in a 3-D structure file.
_NEEDS_3D = (1, 2, 4, 5, 10, 16, 17, 21, 22, 23, 24, 25)
_PLANE = (3, 7, 11, 14)


def _fail(what, line):
    raise ReadError(f"Z88: {what} (line {line})")


def _int(token):
    if token.startswith("+"):
        token = token[1:]
    if not token or not (token.isdigit() or (token[0] == "-" and token[1:].isdigit())):
        return None
    return int(token)


def _real(token):
    try:
        return float(token)
    except ValueError:
        return None


def _leading_ints(line):
    out = []
    for tok in line.split():
        v = _int(tok)
        if v is None:
            break
        out.append(v)
    return out


def _sibling(directory, name):
    try:
        entries = os.listdir(directory or ".")
    except OSError:
        return None
    for entry in sorted(entries):
        path = os.path.join(directory, entry)
        if entry.lower() == name and os.path.isfile(path):
            return path
    return None


def _read_text(path_or_buf):
    with open_file(path_or_buf, "rb") as f:
        data = f.read()
    if isinstance(data, bytes):
        data = data.decode("latin-1")
    return data


def _lines(text):
    return [line[:-1] if line.endswith("\r") else line for line in text.split("\n")]


def _attach_displacements(mesh, path, node_index):
    rows = []
    width = 0
    for line in _lines(_read_text(path)):
        t = line.split()
        if len(t) not in (3, 4, 7):
            continue
        ident = _int(t[0])
        if ident is None or ident not in node_index:
            continue
        values = [_real(v) for v in t[1:]]
        if any(v is None for v in values):
            continue
        width = max(width, len(values))
        rows.append((node_index[ident], values))
    if not rows:
        warn(f"Z88: '{path}' holds no displacement rows; ignored")
        return
    u = np.full((len(mesh.points), width), np.nan)
    for p, values in rows:
        u[p, : len(values)] = values
    mesh.point_data["U"] = u


# z88o3.txt labels that are coordinates, not results.
_COORDINATE_LABELS = {"XX", "YY", "ZZ", "RR", "PHI"}
# The stress tensor the solids and plane-stress elements give, in `SIG` order.
_SOLID_LABELS = ["SIGXX", "SIGYY", "SIGZZ", "TAUXY", "TAUYZ", "TAUZX"]
_PLANE_LABELS = ["SIGXX", "SIGYY", "TAUXY"]


def _o3_label(token):
    return token.split("(", 1)[0]


def _attach_stresses(mesh, path, element_index, cell_code, block_start):
    """`z88o3.txt`: per element a header (`element # = N ...`), a line of column
    labels and rows of reals (a truss prints `SIG = v` on its header line). Each
    element keeps the mean of every label over its rows: solids and plane-stress
    elements as the `SIG` tensor, every other label (beams, shafts, trusses,
    tori, plates, shells) as a scalar array of that name; `SIGV` for all."""
    ncells = len(cell_code)
    sums = [dict() for _ in range(ncells)]
    current = -1
    labels = None
    for line in _lines(_read_text(path)):
        hash_at = line.find("#")
        if hash_at >= 0 and "element" in line[:hash_at].lower():
            current = -1
            labels = None
            rest = line[hash_at + 1 :]
            eq = rest.find("=")
            if eq >= 0:
                ident = _leading_ints(rest[eq + 1 :])
                if ident and ident[0] in element_index:
                    current = element_index[ident[0]]
            at = rest.find("SIG =")
            if current >= 0 and at >= 0:
                t = rest[at + 5 :].split()
                v = _real(t[0]) if t else None
                if v is not None:
                    sums[current]["SIGXX"] = [v, 1]
            continue
        if current < 0:
            continue
        t = line.split()
        if not t:
            continue
        v = [_real(x) for x in t]
        if all(x is None for x in v):
            labels = [_o3_label(x) for x in t]
            continue
        if any(x is None for x in v) or labels is None or len(v) != len(labels):
            continue
        for name, x in zip(labels, v):
            if name in _COORDINATE_LABELS:
                continue
            acc = sums[current].setdefault(name, [0.0, 0])
            acc[0] += x
            acc[1] += 1

    def mean(c, name):
        acc = sums[c].get(name)
        return acc[0] / acc[1] if acc else math.nan

    tensor = [None] * ncells
    scalars = set()
    for c in range(ncells):
        code = cell_code[c]
        if code in _SOLID and all(n in sums[c] for n in _SOLID_LABELS):
            tensor[c] = _SOLID_LABELS
        elif code in _PLANE and all(n in sums[c] for n in _PLANE_LABELS):
            tensor[c] = _PLANE_LABELS
        scalars.update(n for n in sums[c] if n != "SIGV" and n not in (tensor[c] or []))
    widths = {len(t) for t in tensor if t}
    if len(widths) > 1:
        # Plane and solid elements do not share a file; keep the solids.
        tensor = [t if t is _SOLID_LABELS else None for t in tensor]
        widths = {6}

    def per_block(fill):
        return [
            np.array([fill(c) for c in range(block_start[b], block_start[b + 1])])
            for b in range(len(block_start) - 1)
        ]

    if widths:
        w = widths.pop()
        mesh.cell_data["SIG"] = per_block(
            lambda c: [mean(c, n) for n in tensor[c]] if tensor[c] else [math.nan] * w
        )
    if any("SIGV" in s for s in sums):
        mesh.cell_data["SIGV"] = per_block(lambda c: mean(c, "SIGV"))
    for name in sorted(scalars):
        mesh.cell_data[name] = per_block(
            lambda c: math.nan if tensor[c] and name in tensor[c] else mean(c, name)
        )


def _attach_forces(mesh, path, node_index, width):
    """`z88o4.txt`: after the per-element blocks, the sums per node (`node F(1)
    ... F(6)`, introduced by "nodal sums" / "aufsummierten")."""
    rows = []
    started = False
    for line in _lines(_read_text(path)):
        low = line.lower()
        if "nodal sums" in low or "aufsummierten" in low:
            started = True
            continue
        if not started:
            continue
        t = line.split()
        if len(t) != 7:
            continue
        ident = _int(t[0])
        values = [_real(x) for x in t[1:]]
        if ident is None or ident not in node_index or any(x is None for x in values):
            continue
        rows.append((node_index[ident], values[:width]))
    if not rows:
        return
    f = np.full((len(mesh.points), width), np.nan)
    for p, values in rows:
        f[p] = values
    mesh.point_data["F"] = f


def _count_rows(lines, path):
    """A Z88 input file: a count on its first line, then that many rows."""
    rows = []
    head = _leading_ints(lines[0]) if lines else []
    if not head or head[0] < 0:
        warn(f"Z88: '{path}' does not start with a count; ignored")
        return rows
    for line in lines[1:]:
        if len(rows) == head[0]:
            break
        if line.split():
            rows.append(line.split())
    if len(rows) < head[0]:
        warn(f"Z88: '{path}' ends after {len(rows)} of {head[0]} rows")
    return rows


def _attach_constraints(mesh, path, node_index, width):
    """`z88i2.txt`: `node dof flag value` rows; flag 1 adds a force, flag 2
    prescribes a displacement."""
    n = len(mesh.points)
    u = np.full((n, width), np.nan)
    f = np.full((n, width), np.nan)
    skipped = 0
    rows = _count_rows(_lines(_read_text(path)), path)
    if not rows:
        return
    for t in rows:
        vals = [_int(x) for x in t[:3]]
        value = _real(t[3]) if len(t) > 3 else None
        if None in vals or value is None or vals[0] not in node_index:
            skipped += 1
            continue
        node, dof, flag = vals
        if not 1 <= dof <= width or flag not in (1, 2):
            skipped += 1
            continue
        p = node_index[node]
        if flag == 2:
            u[p, dof - 1] = value
        else:
            f[p, dof - 1] = (0.0 if np.isnan(f[p, dof - 1]) else f[p, dof - 1]) + value
    if skipped:
        warn(f"Z88: {skipped} constraint row(s) of '{path}' skipped")
    mesh.point_data["z88:bc:u"] = u
    mesh.point_data["z88:bc:f"] = f


def _ranges(lines, path, element_index, ncells, parse):
    """Rows `from to ...` over element ids: `parse(tokens)` gives the row's
    values, assigned to every cell in the range."""
    out = [None] * ncells
    for t in _count_rows(lines, path):
        ends = [_int(x) for x in t[:2]] if len(t) >= 2 else [None]
        value = parse(t[2:]) if None not in ends else None
        if value is None:
            warn(f"Z88: a row of '{path}' is malformed; skipped")
            continue
        for ident in range(ends[0], ends[1] + 1):
            c = element_index.get(ident)
            if c is not None:
                out[c] = value
    return out


def _attach_sets(mesh, path, node_index, element_index, cell_dim):
    """Z88Aurora's `z88sets.txt`: after a count, each set is a header
    `#KIND PURPOSE id count "name"` and then its ids. Element sets become cell
    regions and node sets point regions (tag = set id); surface sets, whose
    rows do not name the structure file's elements, are skipped."""
    regions = list(getattr(mesh, "regions", []) or [])
    skipped = 0
    header = None
    ids = []

    def flush():
        nonlocal skipped
        if header is None:
            return
        kind, sid, name = header
        if kind == "ELEMENTS":
            cells = sorted({element_index[i] for i in ids if i in element_index})
            dim = max((cell_dim[c] for c in cells), default=-1)
            regions.append(
                Region(name, "cell", np.array(cells, dtype=np.int64), dim, sid)
            )
        elif kind == "NODES":
            pts = sorted({node_index[i] for i in ids if i in node_index})
            regions.append(
                Region(name, "point", np.array(pts, dtype=np.int64), -1, sid)
            )
        else:
            skipped += 1

    for line in _lines(_read_text(path)):
        t = line.split()
        if not t:
            continue
        if t[0].startswith("#"):
            flush()
            quote = line.find('"')
            name = line[quote + 1 : line.rfind('"')] if quote >= 0 else ""
            sid = _int(t[2]) if len(t) > 2 and _int(t[2]) is not None else -1
            header = (t[0][1:].upper(), sid, name or f"set_{sid}")
            ids = []
            continue
        if header is not None:
            ids.extend(v for v in (_int(x) for x in t) if v is not None)
    flush()
    if skipped:
        warn(f"Z88: {skipped} surface set(s) of '{path}' skipped")
    mesh.regions = regions


def _load_layout(code):
    """What a `z88i5.txt` line holds for an element type, as Z88R reads it
    (ri588i.c): the number of values (pressure, then the tangential shears in r
    and s) and of nodes naming the loaded edge or face; (0, 0) for a type that
    takes no surface load. Plates and flat shells take a pressure alone."""
    if code in (7, 8, 14, 15):
        return 2, 3
    if code == 17:
        return 1, 3
    if code in (16, 22):
        return 1, 6
    if code in (10, 21):
        return 3, 8
    if code == 1:
        return 3, 4
    if code in (11, 12):
        return 2, 4
    if code in (18, 19, 20, 23, 24):
        return 1, 0
    return 0, 0


def _attach_loads(mesh, path, node_index, element_index, cell_code):
    """`z88i5.txt`: after a count, one surface load per line, `element values
    nodes` as `_load_layout` says. Kept as field data: `z88:surface_load`
    (loads x 3: pressure, shear r, shear s; NaN where the type has none) and
    `z88:surface_load:cells` (loads x 9: the cell, then up to 8 points, -1 past
    the loaded edge or face's nodes)."""
    lines = [t for t in (line.split() for line in _lines(_read_text(path))) if t]
    count = _int(lines[0][0]) if lines else None
    if not count or count <= 0:
        return
    values, refs = [], []
    skipped = 0
    for t in lines[1 : 1 + count]:
        cell = element_index.get(_int(t[0]))
        if cell is None:
            skipped += 1
            continue
        nv, nn = _load_layout(cell_code[cell])
        if nv == 0 or len(t) < 1 + nv + nn:
            skipped += 1
            continue
        v = [_real(x) for x in t[1 : 1 + nv]]
        nodes = [node_index.get(_int(x)) for x in t[1 + nv : 1 + nv + nn]]
        if any(x is None for x in v) or any(x is None for x in nodes):
            skipped += 1
            continue
        values.append(v + [math.nan] * (3 - nv))
        refs.append([cell] + nodes + [-1] * (8 - nn))
    if skipped:
        warn(
            f"Z88: {skipped} surface load(s) of '{path}' name no element that takes "
            "one, or are malformed; skipped"
        )
    mesh.field_data["z88:surface_load"] = np.array(values, dtype=np.float64).reshape(
        -1, 3
    )
    mesh.field_data["z88:surface_load:cells"] = np.array(refs, dtype=np.int64).reshape(
        -1, 9
    )


def _attach_inputs(mesh, directory, element_index, block_start):
    ncells = block_start[-1]

    def split(values, fill, dtype, width=None):
        blocks = []
        for b in range(len(block_start) - 1):
            rows = [
                fill if v is None else v
                for v in values[block_start[b] : block_start[b + 1]]
            ]
            a = np.array(rows, dtype=dtype)
            if width:
                a = a.reshape(-1, width)
            blocks.append(a)
        return blocks

    mat = _sibling(directory, "z88mat.txt")
    if mat:
        k = [0]

        def material(t):
            k[0] += 1
            if not t:
                return None
            stem = t[0].rsplit(".", 1)[0]
            number = _int(stem) if _int(stem) is not None and _int(stem) > 0 else k[0]
            found = _sibling(directory, t[0].lower())
            props = None
            if found:
                first = _lines(_read_text(found))[0].split()
                props = [_real(x) for x in first[:2]]
            if not props or len(props) < 2 or None in props:
                warn(f"Z88: material file '{t[0]}' is missing or malformed")
                props = [math.nan, math.nan]
            return (number, props[0], props[1])

        values = _ranges(_lines(_read_text(mat)), mat, element_index, ncells, material)
        if any(v is not None for v in values):
            mesh.cell_data["z88:material"] = split(
                [v and v[0] for v in values], 0, np.int64
            )
            mesh.cell_data["z88:E"] = split(
                [v and v[1] for v in values], math.nan, float
            )
            mesh.cell_data["z88:nu"] = split(
                [v and v[2] for v in values], math.nan, float
            )
    elp = _sibling(directory, "z88elp.txt")
    if elp:

        def params(t):
            # Fields a row leaves out stay absent (NaN): Z88R leaves them untouched.
            v = [_real(x) for x in t[:12]]
            return None if None in v else v + [math.nan] * (12 - len(v))

        values = _ranges(_lines(_read_text(elp)), elp, element_index, ncells, params)
        if any(v is not None for v in values):
            mesh.cell_data["z88:elp"] = split(values, [math.nan] * 12, float, 12)
    integ = _sibling(directory, "z88int.txt")
    if integ:

        def orders(t):
            v = [_int(x) for x in t[:2]]
            return None if len(v) < 2 or None in v else v

        values = _ranges(
            _lines(_read_text(integ)), integ, element_index, ncells, orders
        )
        if any(v is not None for v in values):
            mesh.cell_data["z88:int"] = split(values, [-1, -1], np.int64, 2)


def read(filename, results=True):
    structure = filename
    directory = None
    if not is_buffer(filename, "r"):
        structure = os.fspath(filename)
        directory = os.path.dirname(structure)
        base = os.path.basename(structure).lower()
        if base in ("z88o2.txt", "z88o3.txt"):
            found = _sibling(directory, "z88i1.txt") or _sibling(
                directory, "z88structure.txt"
            )
            if found is None:
                raise ReadError(f"Z88: no z88i1.txt next to {structure}")
            structure = found
    lines = _lines(_read_text(structure))
    i = 0
    while i < len(lines) and not lines[i].split():
        i += 1
    if i >= len(lines):
        raise ReadError(f"Z88: {filename} is empty")
    head = _leading_ints(lines[i])
    if len(head) < 3 or not 1 <= head[0] <= 3 or head[1] < 0 or head[2] < 0:
        _fail("the header needs the dimension, node and element counts", i + 1)
    ndim, nnodes, nelem = head[:3]
    legacy = len(head) >= 8
    kflag = head[5] if legacy else (head[4] if len(head) >= 5 else 0)
    i += 1

    coords = np.zeros((nnodes, 3))
    node_index = {}
    node_dof = 0
    for n in range(nnodes):
        while i < len(lines) and not lines[i].split():
            i += 1
        if i >= len(lines):
            _fail(f"the file ends before node {n + 1} of {nnodes}", i)
        t = lines[i].split()
        ident = _int(t[0]) if t else None
        if len(t) < 2 + ndim or ident is None:
            _fail(
                "a node line needs its id, degrees of freedom and "
                f"{ndim} coordinates",
                i + 1,
            )
        x = [0.0, 0.0, 0.0]
        for d in range(ndim):
            v = _real(t[2 + d])
            if v is None:
                _fail(f"bad coordinate '{t[2 + d]}'", i + 1)
            x[d] = v
        if kflag == 1:
            r, phi = x[0], x[1] * math.pi / 180.0
            x[0], x[1] = r * math.cos(phi), r * math.sin(phi)
        if ident in node_index:
            _fail(f"node {ident} is defined twice", i + 1)
        node_index[ident] = n
        coords[n] = x
        node_dof = max(node_dof, _int(t[1]) or 0)
        i += 1
    if kflag == 1:
        warn("Z88: cylindrical input (KFLAG = 1) converted to Cartesian coordinates")

    elements = []
    for e in range(nelem):
        while i < len(lines) and not lines[i].split():
            i += 1
        if i >= len(lines):
            _fail(f"the file ends before element {e + 1} of {nelem}", i)
        h = _leading_ints(lines[i])
        if len(h) < 2:
            _fail("an element starts with its id and type", i + 1)
        info = _TYPES.get(h[1])
        if info is None:
            _fail(f"element {h[0]} has unknown type {h[1]}", i + 1)
        line_no = i + 1
        i += 1
        nodes = []
        while len(nodes) < info[0]:
            if i >= len(lines):
                _fail(f"element {h[0]} is cut short", line_no)
            nodes.extend(_leading_ints(lines[i])[: info[0] - len(nodes)])
            i += 1
        elements.append((h[0], h[1], nodes, line_no))

    order_of_types = []
    by_type = {}
    reduced = set()
    for e, (_, code, _, _) in enumerate(elements):
        info = _TYPES[code]
        if info[1] not in by_type:
            by_type[info[1]] = []
            order_of_types.append(info[1])
        by_type[info[1]].append(e)
        if info[2] < info[0]:
            reduced.add(code)
    for code in sorted(reduced):
        warn(f"Z88: type {code} elements keep only their corner nodes")
        _provenance.note(
            "high-order-dropped", f"Z88 type {code} elements keep only their corners"
        )

    cells = []
    type_blocks = []
    element_index = {}
    cell_code = []
    block_start = [0]
    for cell_type in order_of_types:
        members = by_type[cell_type]
        order = node_order("z88", cell_type)
        k = len(order.to_meshio) if order else _TYPES[elements[members[0]][1]][2]
        conn = np.empty((len(members), k), dtype=np.int64)
        codes = np.empty(len(members), dtype=np.int64)
        for r, e in enumerate(members):
            ident, code, nodes, line_no = elements[e]
            corners = _CORNERS.get(code)
            for j in range(k):
                src = order.to_meshio[j] if order else j
                if corners:
                    src = corners[src]
                idx = node_index.get(nodes[src])
                if idx is None:
                    _fail(f"element {ident} names undefined node {nodes[src]}", line_no)
                conn[r, j] = idx
            codes[r] = code
            if ident in element_index:
                _fail(f"element {ident} is defined twice", line_no)
            element_index[ident] = len(cell_code)
            cell_code.append(code)
        cells.append((cell_type, conn))
        type_blocks.append(codes)
        block_start.append(len(cell_code))

    mesh = Mesh(coords, cells)
    if type_blocks:
        mesh.cell_data["z88:type"] = type_blocks
    width = node_dof if node_dof in (2, 3, 6) else (2 if ndim == 2 else 3)
    if directory is not None:
        i2 = _sibling(directory, "z88i2.txt")
        if i2:
            _attach_constraints(mesh, i2, node_index, width)
        if cell_code:
            _attach_inputs(mesh, directory, element_index, block_start)
        sets = _sibling(directory, "z88sets.txt")
        if sets:
            cell_dim = [topological_dimension[t] for t, conn in cells for _ in conn]
            _attach_sets(mesh, sets, node_index, element_index, cell_dim)
        i5 = _sibling(directory, "z88i5.txt")
        if i5 and cell_code:
            _attach_loads(mesh, i5, node_index, element_index, cell_code)
    if results and directory is not None:
        o2 = _sibling(directory, "z88o2.txt")
        if o2:
            _attach_displacements(mesh, o2, node_index)
        o4 = _sibling(directory, "z88o4.txt")
        if o4:
            _attach_forces(mesh, o4, node_index, width)
        o3 = _sibling(directory, "z88o3.txt")
        if o3 and cell_code:
            _attach_stresses(mesh, o3, element_index, cell_code, block_start)
    return mesh


def _default_code(cell_type, ndim):
    return {
        "hexahedron": 1,
        "hexahedron20": 10,
        "tetra": 17,
        "tetra10": 16,
        "triangle6": 14 if ndim == 2 else 24,
        "quad8": 7 if ndim == 2 else 23,
        "line": 9 if ndim == 2 else 4,
    }.get(cell_type, 0)


def write(filename, mesh, stubs=False):
    points = np.asarray(mesh.points)
    npts = len(points)
    pdim = points.shape[1] if points.ndim == 2 else 0
    any_3d = any(
        not isinstance(b.data, list) and topological_dimension.get(b.type) == 3
        for b in mesh.cells
    )
    flat = pdim < 3 or bool(np.all(points[:, 2] == 0.0))
    has_type = "z88:type" in mesh.cell_data
    # Types that only exist in a 3-D file (solids, 3-D beams, trusses and the
    # shaft, shells) make it 3-D even when every z is 0.
    needs_3d = False
    for b, block in enumerate(mesh.cells):
        if not has_type or isinstance(block.data, list):
            continue
        for w in np.asarray(mesh.cell_data["z88:type"][b]).ravel().tolist():
            info = _TYPES.get(int(w))
            if (
                info
                and info[2] == info[0]
                and info[1] == block.type
                and np.asarray(block.data).reshape(len(block.data), -1).shape[1]
                == info[0]
                and int(w) in _NEEDS_3D
            ):
                needs_3d = True
    ndim = 2 if flat and not any_3d and not needs_3d else 3

    codes = []
    dropped = set()
    nelem = 0
    for b, block in enumerate(mesh.cells):
        ragged = isinstance(block.data, list)
        width = (
            0
            if ragged
            else np.asarray(block.data).reshape(len(block.data), -1).shape[1]
        )
        fallback = 0 if ragged else _default_code(block.type, ndim)
        # A cubic Lagrange quad is a 16-node plate (type 19), in a 2-D file.
        if block.type == "VTK_LAGRANGE_QUADRILATERAL":
            fallback = 19 if width == 16 and ndim == 2 else 0
        want = mesh.cell_data["z88:type"][b] if has_type and not ragged else None
        row = []
        for r in range(len(block.data)):
            code = fallback
            if want is not None:
                w = int(want[r])
                info = _TYPES.get(w)
                if (
                    info
                    and info[2] == info[0]
                    and info[1] == block.type
                    and width == info[0]
                ):
                    code = w
            if not code:
                dropped.add(block.type)
            else:
                nelem += 1
            row.append(code)
        codes.append(row)
    for t in sorted(dropped):
        warn(f"Z88 writer: '{t}' cells have no Z88 element type here; dropped")
        _provenance.note("cells-dropped", f"Z88 has no element type for '{t}' cells")
    deck = [n for n in _DECK_POINT if n in mesh.point_data]
    deck += [n for n in _DECK_CELL if n in mesh.cell_data]
    deck += [n for n in _DECK_FIELD if n in mesh.field_data]
    other = (
        len(mesh.point_data)
        + len(mesh.field_data)
        + len(mesh.cell_data)
        - (1 if has_type else 0)
        - len(deck)
    )
    if other:
        warn(
            "Z88 decks hold no other data arrays; point, cell and field data other "
            "than the z88: constraints, materials and element parameters dropped"
        )
        _provenance.note("data-dropped", "a Z88 deck holds no other data arrays")
    if not nelem:
        raise WriteError("Z88 writer: no cell has a Z88 element type")

    dof = [2 if ndim == 2 else 3] * npts
    seen = [False] * npts
    for b, block in enumerate(mesh.cells):
        if isinstance(block.data, list):
            continue
        for r, row in enumerate(block.data.tolist()):
            code = codes[b][r]
            if not code:
                continue
            d = _TYPES[code][3]
            for p in row:
                dof[p] = max(dof[p], d) if seen[p] else d
                seen[p] = True

    out = [
        "%5d %9d %9d %11d %5d   %s\n"
        % (
            ndim,
            npts,
            nelem,
            sum(dof),
            0,
            _provenance.lines(_provenance.SlotTier.SINGLE_LINE)[0],
        )
    ]
    for p in range(npts):
        parts = ["%9d %2d" % (p + 1, dof[p])]
        for d in range(ndim):
            v = float(points[p, d]) if d < pdim else 0.0
            parts.append(" %+.16E" % v)
        out.append("".join(parts) + "\n")
    ident = 0
    for b, block in enumerate(mesh.cells):
        if isinstance(block.data, list):
            continue
        order = node_order("z88", block.type)
        # Type 19's lattice slot j holds the cell's node lattice[j].
        lattice = [0] * 16
        for j, slot in enumerate(_CORNERS[19]):
            lattice[slot] = j
        for r, row in enumerate(block.data.tolist()):
            if not codes[b][r]:
                continue
            ident += 1
            out.append("%9d %5d\n" % (ident, codes[b][r]))
            if codes[b][r] == 19:
                src = lattice
            else:
                src = order.from_meshio if order else range(len(row))
            out.append(" ".join(str(row[s] + 1) for s in src) + "\n")
    # "\n" on every platform, as the C++ writer: the engines write the same bytes.
    with open_file(filename, "w", newline="\n") as f:
        f.write("".join(out))
    if is_buffer(filename, "w"):
        return
    directory = os.path.dirname(os.fspath(filename))
    files = _deck_files(mesh, codes, dof)
    if stubs:
        files.setdefault("z88i2.txt", "0\n")
        files.setdefault("z88i5.txt", "0\n")
    for name, text in files.items():
        with open(os.path.join(directory, name), "w", newline="\n") as f:
            f.write(text)


_DECK_POINT = ("z88:bc:u", "z88:bc:f")
_DECK_CELL = ("z88:material", "z88:E", "z88:nu", "z88:elp", "z88:int")
_DECK_FIELD = ("z88:surface_load", "z88:surface_load:cells")


def _deck_files(mesh, codes, dof):
    """The input files the `z88:` arrays describe: `z88i2.txt` from the
    constraints, `z88mat.txt` and one `<n>.txt` per material, `z88elp.txt`,
    `z88int.txt`; element ranges over the written element ids."""
    files = {}
    u = mesh.point_data.get("z88:bc:u")
    f = mesh.point_data.get("z88:bc:f")
    if u is not None or f is not None:
        n = len(mesh.points)
        kinds = [
            (flag, np.asarray(a, dtype=float).reshape(n, -1))
            for flag, a in ((2, u), (1, f))
            if a is not None
        ]
        rows = []
        for p in range(n):
            for d in range(dof[p]):
                for flag, a in kinds:
                    if d < a.shape[1] and not np.isnan(a[p, d]):
                        rows.append(
                            "%9d %2d %2d %+.16E\n" % (p + 1, d + 1, flag, a[p, d])
                        )
        files["z88i2.txt"] = f"{len(rows)}\n" + "".join(rows)

    def written(name, width=None):
        """The array's rows of the written elements, in element order."""
        out = []
        for b, arr in enumerate(mesh.cell_data[name]):
            arr = np.asarray(arr)
            arr = arr.reshape(len(arr), -1) if width else arr.reshape(len(arr))
            out += [arr[r].tolist() for r in range(len(codes[b])) if codes[b][r]]
        return out

    def ranges(values):
        """Runs of equal values as (first id, last id, value); None skipped."""
        runs = []
        for k, v in enumerate(values, start=1):
            if runs and runs[-1][2] == v and runs[-1][1] == k - 1:
                runs[-1][1] = k
            else:
                runs.append([k, k, v])
        return [r for r in runs if r[2] is not None]

    if "z88:E" in mesh.cell_data and "z88:nu" in mesh.cell_data:
        e = written("z88:E")
        nu = written("z88:nu")
        want = written("z88:material") if "z88:material" in mesh.cell_data else None
        number = {}
        for k in range(len(e)):
            key = (e[k], nu[k])
            if math.isnan(e[k]) or math.isnan(nu[k]) or key in number:
                continue
            w = int(want[k]) if want is not None else 0
            if w <= 0 or w in number.values():
                w = max([0] + list(number.values())) + 1
            number[key] = w
        values = [
            None if (math.isnan(a) or math.isnan(b)) else number[(a, b)]
            for a, b in zip(e, nu)
        ]
        runs = ranges(values)
        if runs:
            files["z88mat.txt"] = f"{len(runs)}\n" + "".join(
                "%9d %9d %d.txt\n" % tuple(r) for r in runs
            )
            for (a, b), n in number.items():
                files[f"{n}.txt"] = "%+.16E %+.16E\n" % (a, b)
    if "z88:elp" in mesh.cell_data:
        # The leading fields up to the first NaN, at least the seven section values.
        rows = []
        for v in written("z88:elp", 12):
            n = next((k for k, x in enumerate(v) if math.isnan(x)), len(v))
            rows.append(tuple(v[:n]) if n >= 7 else None)
        runs = ranges(rows)
        if runs:
            files["z88elp.txt"] = f"{len(runs)}\n" + "".join(
                "%9d %9d" % (a, b)
                + "".join(
                    (" %d" % int(x)) if k == 7 else (" %+.16E" % x)
                    for k, x in enumerate(v)
                )
                + "\n"
                for a, b, v in runs
            )
    if "z88:int" in mesh.cell_data:
        rows = [
            None if v[0] < 0 else (int(v[0]), int(v[1])) for v in written("z88:int", 2)
        ]
        runs = ranges(rows)
        if runs:
            files["z88int.txt"] = f"{len(runs)}\n" + "".join(
                "%9d %9d %d %d\n" % (a, b, *v) for a, b, v in runs
            )
    # The written element id of every cell (0: dropped), block-major.
    written_id, code_of = [], []
    ident = 0
    for block in codes:
        for code in block:
            ident += 1 if code else 0
            written_id.append(ident if code else 0)
            code_of.append(code)
    values = mesh.field_data.get("z88:surface_load")
    refs = mesh.field_data.get("z88:surface_load:cells")
    if values is not None and refs is not None:
        values = np.asarray(values, dtype=float).reshape(-1, 3)
        refs = np.asarray(refs, dtype=np.int64).reshape(-1, 9)
        rows, skipped = [], 0
        for v, r in zip(values, refs):
            cell = int(r[0])
            code = code_of[cell] if 0 <= cell < len(code_of) else 0
            nv, nn = _load_layout(code)
            given = next((k for k in range(8) if r[1 + k] < 0), 8)
            if not code or nv == 0 or given != nn:
                skipped += 1
                continue
            rows.append(
                str(written_id[cell])
                + "".join(" %+.16E" % (0.0 if math.isnan(x) else x) for x in v[:nv])
                + "".join(" %d" % (p + 1) for p in r[1 : 1 + nn])
                + "\n"
            )
        if skipped:
            warn(
                f"Z88 writer: {skipped} surface load(s) name a cell that is not "
                "written, or whose type takes no such load; dropped"
            )
        files["z88i5.txt"] = f"{len(rows)}\n" + "".join(rows)
    # Z88Aurora's sets: element sets from cell regions, node sets from point
    # regions (Aurora's purposes for them: MATERIAL and CONSTRAINT).
    # In the core's region order: kind (point, cell, side), name, dim, tag.
    rank = {"point": 0, "cell": 1, "side": 2}
    regions = sorted(
        getattr(mesh, "regions", []) or [],
        key=lambda r: (rank.get(r.kind, 3), r.name, r.dim, r.tag),
    )
    used = {r.tag for r in regions if r.kind in ("cell", "point") and r.tag > 0}
    claimed = set()
    nxt = 1
    sets, dropped = [], 0
    for region in regions:
        if region.kind not in ("cell", "point"):
            dropped += 1
            continue
        entries = np.asarray(region.entries, dtype=np.int64).ravel().tolist()
        if region.kind == "cell":
            ids = [
                written_id[k]
                for k in entries
                if 0 <= k < len(written_id) and written_id[k]
            ]
        else:
            ids = [k + 1 for k in entries]
        # The region's tag as the set id, unless another set took it first.
        tag = int(region.tag)
        if tag <= 0 or tag in claimed:
            while nxt in used or nxt in claimed:
                nxt += 1
            tag = nxt
            nxt += 1
        claimed.add(tag)
        head = "#ELEMENTS MATERIAL" if region.kind == "cell" else "#NODES CONSTRAINT"
        lines = [f'{head} {tag} {len(ids)} "{region.name}"\n']
        for k in range(0, len(ids), 10):
            lines.append("".join("%10d " % i for i in ids[k : k + 10]) + "\n")
        sets.append("".join(lines))
    if dropped:
        warn(
            f"Z88 writer: {dropped} region(s) other than cell and point regions dropped"
        )
        _provenance.note(
            "regions-dropped", "Z88Aurora sets hold elements and nodes only"
        )
    if sets:
        files["z88sets.txt"] = f"{len(sets)}\n" + "".join(sets)
    return files

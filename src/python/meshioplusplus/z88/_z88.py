"""I/O for Z88 / Z88Aurora structure files (``z88i1.txt``) and their results.

The pure-Python twin of ``src/cpp/src/formats/z88.cpp``: both engines read the
same meshes and write the same bytes.

``z88i1.txt``: a header whose first three integers are the dimension, node and
element counts (Z88OS v15: ``ndim nnodes nelem ndof kflag``; Z88 <= V13 and
Z88Aurora V1 add more flags and material lines after the elements); one line per
node, ``id ndof x y [z]``; two lines per element, ``id type`` then its nodes.
``z88o2.txt`` (displacements) and ``z88o3.txt`` (stresses) next to it are
attached as ``point_data["U"]`` and ``cell_data["SIG"]``/``["SIGV"]``.
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
    19: (16, "quad", 4, 3),
    20: (8, "quad8", 8, 3),
    21: (16, "hexahedron", 8, 3),
    22: (12, "wedge", 6, 3),
    23: (8, "quad8", 8, 6),
    24: (6, "triangle6", 6, 6),
    25: (2, "line", 2, 6),
}

# Corner slots of the types that keep only their corners, when not the leading
# nodes: 19 is a 4x4 lattice row by row (corners 1, 13, 16, 4 counter-clockwise),
# 21/22 two quad8/tri6 layers.
_CORNERS = {
    19: [0, 12, 15, 3],
    21: [0, 1, 2, 3, 8, 9, 10, 11],
    22: [0, 1, 2, 6, 7, 8],
}

_SOLID = (1, 10, 16, 17)
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


def _attach_stresses(mesh, path, element_index, cell_code, block_start):
    ncells = len(cell_code)
    sums = [None] * ncells
    sumv = [0.0] * ncells
    count = [0] * ncells
    countv = [0] * ncells
    current = -1
    skipped = 0
    components = 0
    for line in _lines(_read_text(path)):
        hash_at = line.find("#")
        if hash_at >= 0 and "element" in line[:hash_at].lower():
            current = -1
            rest = line[hash_at + 1 :]
            eq = rest.find("=")
            if eq >= 0:
                ident = _leading_ints(rest[eq + 1 :])
                if ident and ident[0] in element_index:
                    current = element_index[ident[0]]
            continue
        if current < 0:
            continue
        t = line.split()
        if not t:
            continue
        v = [_real(x) for x in t]
        if any(x is None for x in v):
            continue
        code = cell_code[current]
        if code in _SOLID and len(v) in (9, 10):
            first, n = 3, 6
        elif code in _PLANE and len(v) in (5, 6):
            first, n = 2, 3
        else:
            skipped += 1
            continue
        if components and components != n:
            skipped += 1
            continue
        components = n
        if sums[current] is None:
            sums[current] = [0.0] * n
        for k in range(n):
            sums[current][k] += v[first + k]
        count[current] += 1
        if len(v) == first + n + 1:
            sumv[current] += v[-1]
            countv[current] += 1
    if skipped:
        warn(
            f"Z88: {skipped} stress row(s) of elements other than solids and "
            "plane-stress elements skipped"
        )
    if not components:
        return
    sig, sigv = [], []
    for b in range(len(block_start) - 1):
        n = block_start[b + 1] - block_start[b]
        a = np.full((n, components), np.nan)
        av = np.full(n, np.nan)
        for r in range(n):
            c = block_start[b] + r
            if count[c]:
                a[r] = np.array(sums[c]) / count[c]
            if countv[c]:
                av[r] = sumv[c] / countv[c]
        sig.append(a)
        sigv.append(av)
    mesh.cell_data["SIG"] = sig
    if any(countv):
        mesh.cell_data["SIGV"] = sigv


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
    if results and directory is not None:
        o2 = _sibling(directory, "z88o2.txt")
        if o2:
            _attach_displacements(mesh, o2, node_index)
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
    ndim = 2 if flat and not any_3d else 3

    has_type = "z88:type" in mesh.cell_data
    codes = []
    dropped = set()
    nelem = 0
    for b, block in enumerate(mesh.cells):
        ragged = isinstance(block.data, list)
        fallback = 0 if ragged else _default_code(block.type, ndim)
        want = mesh.cell_data["z88:type"][b] if has_type and not ragged else None
        row = []
        for r in range(len(block.data)):
            code = fallback
            if want is not None:
                w = int(want[r])
                info = _TYPES.get(w)
                if info and info[2] == info[0] and info[1] == block.type:
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
    regions = getattr(mesh, "regions", []) or []
    if regions:
        warn(f"Z88 structure files hold no groups; {len(regions)} region(s) dropped")
        _provenance.note("regions-dropped", "a Z88 structure file holds no groups")
    other = (
        len(mesh.point_data)
        + len(mesh.field_data)
        + len(mesh.cell_data)
        - (1 if has_type else 0)
    )
    if other:
        warn(
            "Z88 structure files hold no data arrays; point, cell and field data "
            "dropped"
        )
        _provenance.note("data-dropped", "a Z88 structure file holds no data arrays")
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
        for r, row in enumerate(block.data.tolist()):
            if not codes[b][r]:
                continue
            ident += 1
            out.append("%9d %5d\n" % (ident, codes[b][r]))
            src = order.from_meshio if order else range(len(row))
            out.append(" ".join(str(row[s] + 1) for s in src) + "\n")
    # "\n" on every platform, as the C++ writer: the engines write the same bytes.
    with open_file(filename, "w", newline="\n") as f:
        f.write("".join(out))
    if stubs and not is_buffer(filename, "w"):
        directory = os.path.dirname(os.fspath(filename))
        for name in ("z88i2.txt", "z88i5.txt"):
            with open(os.path.join(directory, name), "w", newline="\n") as f:
                f.write("0\n")

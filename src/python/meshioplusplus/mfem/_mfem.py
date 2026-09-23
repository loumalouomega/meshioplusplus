"""I/O for MFEM meshes (``.mesh``) and grid functions (``.gf``).

The pure-Python twin of ``src/cpp/src/formats/mfem.cpp``: both engines read the
same meshes and write the same bytes.

``MFEM mesh v1.0``/``v1.2``/``v1.3``: ``dimension``, ``elements`` (``attribute
geometry vertex...``), ``boundary``, ``vertices`` and, in v1.3, named
``attribute_sets``/``bdr_attribute_sets``. A curved mesh gives its geometry as a
``nodes`` grid function instead of vertex coordinates.

- Elements and boundary elements are separate cell blocks (elements first); the
  attribute is the ``mfem:attribute`` cell data and an ``attribute_<n>`` /
  ``boundary_<n>`` cell region tagged ``n``; an attribute set is a cell region
  of its own name.
- ``H1`` order-2 nodes give quadratic cells; the points are MFEM's degrees of
  freedom in MFEM's order (vertices, then edges and faces by first appearance,
  then element interiors). Higher orders keep the vertices only, with a warning;
  ``L2_T1_<d>D_P1`` nodes give every element its own points.
- Every geometry, the prism included, is in meshio++'s node order.
"""

import os
import pathlib

import numpy as np

from .. import _provenance
from .._common import warn
from .._exceptions import ReadError, WriteError
from .._facets import facet_nodes
from .._files import open_file
from .._mesh import Mesh
from .._regions import Region

__all__ = ["read", "write"]

# (linear type, quadratic type, dim, vertices, MFEM local edges, MFEM local
# faces, interior dof at order 2) per MFEM geometry code.
_GEOMS = [
    ("vertex", "vertex", 0, 1, [], [], False),
    ("line", "line3", 1, 2, [], [], True),
    ("triangle", "triangle6", 2, 3, [(0, 1), (1, 2), (2, 0)], [], False),
    ("quad", "quad9", 2, 4, [(0, 1), (1, 2), (2, 3), (3, 0)], [], True),
    (
        "tetra",
        "tetra10",
        3,
        4,
        [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)],
        [(1, 2, 3), (0, 3, 2), (0, 1, 3), (0, 2, 1)],
        False,
    ),
    (
        "hexahedron",
        "hexahedron27",
        3,
        8,
        [(0, 1), (1, 2), (3, 2), (0, 3), (4, 5), (5, 6)]
        + [(7, 6), (4, 7), (0, 4), (1, 5), (2, 6), (3, 7)],
        [
            (3, 2, 1, 0),
            (0, 1, 5, 4),
            (1, 2, 6, 5),
            (2, 3, 7, 6),
            (3, 0, 4, 7),
            (4, 5, 6, 7),
        ],
        True,
    ),
    (
        "wedge",
        "wedge18",
        3,
        6,
        [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (0, 3), (1, 4), (2, 5)],
        [(0, 2, 1), (3, 4, 5), (0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)],
        False,
    ),
    (
        "pyramid",
        None,
        3,
        5,
        [(0, 1), (1, 2), (3, 2), (0, 3), (0, 4), (1, 4), (2, 4), (3, 4)],
        [(3, 2, 1, 0), (0, 1, 4), (1, 2, 4), (2, 3, 4), (3, 0, 4)],
        False,
    ),
]

# Non-corner nodes of each order-2 meshio++ type, by the corners they sit
# between (2: an edge, 4: a quad face, all corners: the cell centre).
_SLOTS = {
    "vertex": [],
    "line3": [(0, 1)],
    "triangle6": [(0, 1), (1, 2), (2, 0)],
    "quad9": [(0, 1), (1, 2), (2, 3), (3, 0), (0, 1, 2, 3)],
    "tetra10": [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)],
    "wedge18": [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (0, 3), (1, 4), (2, 5)]
    + [(0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)],
    "hexahedron27": [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4)]
    + [(0, 4), (1, 5), (2, 6), (3, 7)]
    + [
        (0, 4, 7, 3),
        (1, 2, 6, 5),
        (0, 1, 5, 4),
        (3, 2, 6, 7),
        (0, 1, 2, 3),
        (4, 5, 6, 7),
    ]
    + [(0, 1, 2, 3, 4, 5, 6, 7)],
}

_WRITABLE_EXTRA = ("quad8", "hexahedron20", "wedge15", "pyramid13", "pyramid14")


def _geom_of_type(cell_type):
    # Serendipity and quadratic pyramids: written by completing them (or as
    # their corners); never read.
    extra = {
        "quad8": 3,
        "hexahedron20": 5,
        "wedge15": 6,
        "pyramid13": 7,
        "pyramid14": 7,
    }
    if cell_type in extra:
        return extra[cell_type]
    for k, g in enumerate(_GEOMS):
        if cell_type == g[0] or (g[1] is not None and cell_type == g[1]):
            return k
    return -1


def _vtk_corner(geom, k):
    """The MFEM vertex meshio++ corner ``k`` comes from: always ``k``. MFEM's
    reference cells are meshio++'s, the prism included; MFEM's own VTK export
    reverses prisms only because classic VTK winds the wedge the other way."""
    return k


class _Numbering:
    """MFEM's order-2 H1 numbering: vertices, edges and (3-D) quad faces by first
    appearance, then element interiors."""

    def __init__(self, elements, dim, nv):
        self.nv = nv
        self.edges = {}
        self.faces = {}
        self.interior = []
        self.keys = [(v,) for v in range(nv)]
        if dim >= 2:
            for _, geom, verts, _ in elements:
                for a, b in _GEOMS[geom][4]:
                    key = tuple(sorted((verts[a], verts[b])))
                    if key not in self.edges:
                        self.edges[key] = len(self.keys)
                        self.keys.append(key)
        if dim == 3:
            seen = set()
            quads = []
            for _, geom, verts, _ in elements:
                for face in _GEOMS[geom][5]:
                    key = tuple(sorted(verts[k] for k in face))
                    if key not in seen:
                        seen.add(key)
                        if len(key) == 4:
                            quads.append(key)
            for key in quads:
                self.faces[key] = len(self.keys)
                self.keys.append(key)
        for _, geom, verts, _ in elements:
            g = _GEOMS[geom]
            if g[6] and g[2] == dim:
                self.interior.append(len(self.keys))
                self.keys.append(tuple(sorted(verts)))
            else:
                self.interior.append(-1)

    def __len__(self):
        return len(self.keys)

    def slot_dof(self, dim, slot, num_corners, element):
        n = len(slot)
        if (
            n == num_corners
            and element >= 0
            and (dim != 3 or n != 4)
            and (dim != 2 or n != 2)
        ):
            return self.interior[element]
        key = tuple(sorted(slot))
        if n == 2 and dim >= 2:
            return self.edges.get(key, -1)
        if n == 4 and dim == 3:
            return self.faces.get(key, -1)
        if element >= 0:
            return self.interior[element]
        return -1


# --- tokenizer --------------------------------------------------------------


class _Lexer:
    """Whitespace tokens; ``#`` lines are comments; ``"..."`` is one token; the
    first line is the header."""

    def __init__(self, what, text):
        self.what = what
        self.tokens = []  # (text, line, quoted)
        self.pos = 0
        self.header = None
        self.header_line = 0
        line_no = 0
        for raw in text.split("\n"):
            line_no += 1
            line = raw[:-1] if raw.endswith("\r") else raw
            stripped = line.lstrip(" \t")
            if not stripped or stripped[0] == "#":
                continue
            if self.header is None:
                self.header = stripped.rstrip(" \t")
                self.header_line = line_no
                continue
            if stripped.startswith(("FiniteElementCollection:", "VDim:", "Ordering:")):
                colon = stripped.index(":")
                self.tokens.append((stripped[: colon + 1], line_no, False))
                self.tokens.append((stripped[colon + 1 :].strip(" \t"), line_no, True))
                continue
            i = 0
            n = len(stripped)
            while i < n:
                while i < n and stripped[i] in " \t":
                    i += 1
                if i >= n:
                    break
                if stripped[i] == '"':
                    i += 1
                    out = []
                    while i < n and stripped[i] != '"':
                        if stripped[i] == "\\" and i + 1 < n:
                            i += 1
                        out.append(stripped[i])
                        i += 1
                    if i >= n:
                        self.fail("unterminated quoted name", line_no)
                    i += 1
                    self.tokens.append(("".join(out), line_no, True))
                    continue
                start = i
                while i < n and stripped[i] not in " \t":
                    i += 1
                self.tokens.append((stripped[start:i], line_no, False))
        self.end_line = line_no

    def fail(self, why, line):
        raise ReadError(f"{self.what}: {why} (line {line})")

    def at_end(self):
        return self.pos >= len(self.tokens)

    def peek(self):
        return self.tokens[self.pos][0]

    def line(self):
        return self.end_line if self.at_end() else self.tokens[self.pos][1]

    def next(self, expected):
        if self.at_end():
            self.fail(f"the file ends where {expected} was expected", self.end_line)
        tok = self.tokens[self.pos]
        self.pos += 1
        return tok

    def int(self, expected):
        text, line, _ = self.next(expected)
        v = _parse_int(text)
        if v is None:
            self.fail(f"expected {expected}, found '{text}'", line)
        return v

    def real(self, expected):
        text, line, _ = self.next(expected)
        v = _parse_real(text)
        if v is None:
            self.fail(f"expected {expected}, found '{text}'", line)
        return v

    def reals(self):
        out = []
        while not self.at_end():
            v = _parse_real(self.tokens[self.pos][0])
            if v is None:
                break
            out.append(v)
            self.pos += 1
        return out


def _parse_int(text):
    body = text[1:] if text[:1] in "+-" else text
    if not body or not body.isascii() or not body.isdigit():
        return None
    return int(text)


def _parse_real(text):
    if not text or text[0] not in "0123456789+-.":
        return None
    try:
        return float(text)
    except ValueError:
        return None


class _Space:
    def __init__(self, name):
        self.collection = name
        self.kind = "other"
        self.order = -1
        self.vdim = 1
        self.ordering = 0

        def order_after_p():
            p = name.rfind("_P")
            v = _parse_int(name[p + 2 :]) if p >= 0 else None
            return -1 if v is None else v

        if name == "Linear":
            self.kind, self.order = "h1", 1
        elif name == "Quadratic":
            self.kind, self.order = "h1", 2
        elif name in ("QuadraticPos", "Cubic"):
            self.kind, self.order = "h1-other", 3 if name == "Cubic" else 2
        elif name.startswith(("H1_", "H1@")):
            self.order = order_after_p()
            self.kind = "h1" if 1 <= self.order <= 2 else "h1-other"
        elif name.startswith(("H1Pos_", "H1Ser_")):
            self.order = order_after_p()
            self.kind = "h1" if self.order == 1 else "h1-other"
        elif name.startswith("L2_T1_"):
            self.order = order_after_p()
            self.kind = "l2t1"
        elif name.startswith("L2_"):
            self.order = order_after_p()
            self.kind = "l2"


def _read_space_body(lex, start_line):
    space = None
    vdim, ordering = 1, 0
    while not lex.at_end():
        key = lex.peek()
        if key == "FiniteElementCollection:":
            lex.next("a collection")
            space = _Space(lex.next("a collection name")[0])
        elif key == "VDim:":
            lex.next("VDim")
            text, line, _ = lex.next("a VDim")
            v = _parse_int(text)
            if v is None or v < 1:
                lex.fail(f"bad VDim '{text}'", line)
            vdim = v
        elif key == "Ordering:":
            lex.next("Ordering")
            text, line, _ = lex.next("an Ordering")
            v = _parse_int(text)
            if v not in (0, 1):
                lex.fail(f"bad Ordering '{text}'", line)
            ordering = v
        else:
            break
    if space is None:
        return None
    space.vdim, space.ordering = vdim, ordering
    return space


def _value(values, space, ndofs, dof, comp):
    if space.ordering == 0:
        return values[comp * ndofs + dof]
    return values[dof * space.vdim + comp]


def _read_text(filename):
    with open_file(filename, "rb") as f:
        raw = f.read()
    return raw.decode("utf-8", "surrogateescape") if isinstance(raw, bytes) else raw


def _read_elements(lex, what):
    n = lex.int("an element count")
    if n < 0:
        lex.fail(f"negative {what} count", lex.line())
    out = []
    for _ in range(n):
        line = lex.line()
        attr = lex.int("an attribute")
        geom = lex.int("a geometry type")
        if geom < 0 or geom > 7:
            lex.fail(f"unknown geometry type {geom}", line)
        verts = []
        for _ in range(_GEOMS[geom][3]):
            v = lex.int("a vertex index")
            if v < 0:
                lex.fail("negative vertex index", line)
            verts.append(v)
        out.append((attr, geom, verts, line))
    return out


def _read_sets(lex):
    n = lex.int("an attribute set count")
    out = []
    for _ in range(n):
        name = lex.next("an attribute set name")[0]
        size = lex.int("an attribute set size")
        attrs = sorted({lex.int("an attribute") for _ in range(size)})
        out.append((name, attrs))
    return out


def _parse(filename):
    lex = _Lexer("MFEM mesh", _read_text(filename))
    header = lex.header or ""
    if header.startswith("MFEM NC mesh"):
        lex.fail(
            f"non-conforming meshes ('{header}') are not supported", lex.header_line
        )
    if header.startswith(("MFEM NURBS", "MFEM INLINE")):
        lex.fail(f"'{header}' meshes are not supported", lex.header_line)
    if header not in (
        "MFEM mesh v1.0",
        "MFEM mesh v1.1",
        "MFEM mesh v1.2",
        "MFEM mesh v1.3",
    ):
        lex.fail(f"not an MFEM mesh (the first line is '{header}')", lex.header_line)
    f = {
        "dim": -1,
        "elements": [],
        "boundary": [],
        "sets": [],
        "bdr_sets": [],
        "nv": 0,
        "sdim": 0,
        "coords": [],
        "nodes_space": None,
        "nodes": None,
    }
    saw_vertices = False
    while not lex.at_end():
        text, line, _ = lex.next("a section")
        if text == "dimension":
            d = lex.int("a dimension")
            if d < 1 or d > 3:
                lex.fail(f"dimension {d} (1, 2 or 3)", line)
            f["dim"] = d
        elif text == "elements":
            f["elements"] = _read_elements(lex, "element")
        elif text == "boundary":
            f["boundary"] = _read_elements(lex, "boundary element")
        elif text == "attribute_sets":
            f["sets"] = _read_sets(lex)
        elif text == "bdr_attribute_sets":
            f["bdr_sets"] = _read_sets(lex)
        elif text == "vertices":
            saw_vertices = True
            nv = lex.int("a vertex count")
            if nv < 0:
                lex.fail("negative vertex count", line)
            f["nv"] = nv
            if not lex.at_end() and lex.peek() == "nodes":
                lex.next("nodes")
                fes, fes_line, _ = lex.next("FiniteElementSpace")
                if fes != "FiniteElementSpace":
                    lex.fail(
                        f"expected FiniteElementSpace, found '{fes}' (NURBS and "
                        "variable-order spaces are not supported)",
                        fes_line,
                    )
                space = _read_space_body(lex, fes_line)
                if space is None:
                    lex.fail(
                        "a FiniteElementSpace without a FiniteElementCollection",
                        fes_line,
                    )
                f["nodes_space"] = space
                f["nodes"] = lex.reals()
                f["sdim"] = space.vdim
            else:
                sd = lex.int("a space dimension")
                if sd < 1 or sd > 3:
                    lex.fail(f"space dimension {sd} (1, 2 or 3)", line)
                f["sdim"] = sd
                f["coords"] = [lex.real("a coordinate") for _ in range(nv * sd)]
        elif text in ("vertex_parents", "coarse_elements"):
            lex.fail(
                f"non-conforming meshes (a '{text}' section) are not supported", line
            )
        elif text == "mfem_serial_mesh_end":
            warn(
                f"MFEM mesh: '{filename}' is one rank of a parallel mesh; reading its "
                "local part and ignoring the communication groups"
            )
            break
        elif text == "mfem_mesh_end":
            break
        else:
            lex.fail(f"unexpected '{text}'", line)
    if f["dim"] < 0:
        raise ReadError(f"MFEM mesh: no dimension section in {filename}")
    if not saw_vertices:
        raise ReadError(f"MFEM mesh: no vertices section in {filename}")
    for group in (f["elements"], f["boundary"]):
        for _, _, verts, line in group:
            for v in verts:
                if v >= f["nv"]:
                    lex.fail(f"vertex {v} out of range ({f['nv']} vertices)", line)
    return f


def _parse_gf(name, path):
    lex = _Lexer("MFEM grid function", "\n" + _read_text(path))
    if lex.header != "FiniteElementSpace":
        raise ReadError(
            f"MFEM grid function: '{path}' does not start with FiniteElementSpace "
            "(NURBS and variable-order spaces are not supported)"
        )
    space = _read_space_body(lex, 1)
    if space is None:
        raise ReadError(
            f"MFEM grid function: '{path}' names no FiniteElementCollection"
        )
    values = lex.reals()
    if not lex.at_end():
        lex.fail(f"unexpected '{lex.peek()}'", lex.line())
    return name, space, values


def read(filename, grid_functions=None):
    f = _parse(filename)
    dim = f["dim"]
    nv = f["nv"]
    elements = f["elements"]
    sdim = f["sdim"]

    coords = "vertices"
    mesh_order = 1
    node_dofs = 0
    nspace = f["nodes_space"]
    nodes = f["nodes"]
    if nspace is not None:
        if len(nodes) % nspace.vdim:
            raise ReadError(
                f"MFEM mesh: the nodes of {filename} hold {len(nodes)} values, not a "
                f"multiple of VDim {nspace.vdim}"
            )
        node_dofs = len(nodes) // nspace.vdim
        if nspace.kind == "h1":
            coords = "h1"
            mesh_order = nspace.order
        elif nspace.kind == "h1-other":
            coords = "corners"
            warn(
                f"MFEM mesh: nodes in '{nspace.collection}' (order {nspace.order}) are "
                "read at the vertices only; meshio++ keeps curved cells up to order 2"
            )
        elif nspace.kind == "l2t1" and nspace.order == 1:
            coords = "dg"
        else:
            raise ReadError(
                f"MFEM mesh: nodes in '{nspace.collection}' are not supported"
            )
        if coords != "dg" and node_dofs < nv:
            raise ReadError(
                f"MFEM mesh: the nodes of {filename} have {node_dofs} dofs for {nv} "
                "vertices"
            )
    has_pyramid = any(el[1] == 7 for el in elements)

    gfs = []
    if isinstance(grid_functions, dict):
        gf_items = list(grid_functions.items())
    else:
        gf_items = [(pathlib.Path(p).stem, p) for p in (grid_functions or [])]
    for name, path in gf_items:
        name, space, values = _parse_gf(str(name), str(path))
        h1 = space.kind == "h1"
        l2p0 = space.kind in ("l2", "l2t1") and space.order == 0
        if not h1 and not l2p0:
            warn(
                f"MFEM grid function '{path}': the '{space.collection}' space is not "
                "supported; skipped"
            )
            continue
        if h1 and space.order == 2 and (coords == "dg" or has_pyramid):
            warn(
                f"MFEM grid function '{path}': an order-2 field on this mesh is not "
                "supported; skipped"
            )
            continue
        gfs.append((name, space, values))
    order = 1 if coords == "dg" else mesh_order
    for _, space, _ in gfs:
        if space.kind == "h1":
            order = max(order, space.order)
    if order == 2 and has_pyramid:
        warn(
            f"MFEM mesh: '{filename}' has pyramids, which meshio++ holds at order 1 "
            "only; reading the vertices"
        )
        order = 1
        if coords == "h1" and mesh_order == 2:
            coords = "corners"
    if order == 2:
        for _, geom, _, line in elements:
            if geom == 0 or _GEOMS[geom][2] != dim:
                raise ReadError(
                    f"MFEM mesh: element of geometry {geom} in a {dim}-D mesh "
                    f"(line {line})"
                )

    if order == 2:
        numbering = _Numbering(elements, dim, nv)
    else:
        numbering = _Numbering([], 0, nv)
    if coords == "h1" and mesh_order == 2 and node_dofs != len(numbering):
        raise ReadError(
            f"MFEM mesh: the order-2 nodes of {filename} have {node_dofs} dofs; the "
            f"mesh numbers {len(numbering)}"
        )
    pdim = sdim
    if coords == "vertices":
        vxyz = np.array(f["coords"], dtype=np.float64).reshape(nv, pdim)
    elif coords in ("h1", "corners"):
        vxyz = np.array(
            [
                [_value(nodes, nspace, node_dofs, v, c) for c in range(pdim)]
                for v in range(nv)
            ],
            dtype=np.float64,
        ).reshape(nv, pdim)
    else:
        vxyz = None

    point_vertex = []
    dg_points = []
    if coords == "dg":
        quad_lex = [0, 1, 3, 2]
        hex_lex = [0, 1, 3, 2, 4, 5, 7, 6]
        xyz = []
        offset = 0
        for _, geom, verts, _ in elements:
            if geom not in (1, 2, 3, 4, 5):
                raise ReadError(
                    f"MFEM mesh: discontinuous nodes on geometry {geom} are not supported"
                )
            nvert = _GEOMS[geom][3]
            pts = []
            for j in range(nvert):
                local = quad_lex[j] if geom == 3 else (hex_lex[j] if geom == 5 else j)
                if offset + local >= node_dofs:
                    raise ReadError(
                        f"MFEM mesh: the discontinuous nodes of {filename} are too few "
                        "for its elements"
                    )
                pts.append(len(point_vertex))
                point_vertex.append(verts[j])
                xyz.append(
                    [
                        _value(nodes, nspace, node_dofs, offset + local, c)
                        for c in range(pdim)
                    ]
                )
            offset += nvert
            dg_points.append(pts)
        if offset != node_dofs:
            raise ReadError(
                f"MFEM mesh: the discontinuous nodes of {filename} have {node_dofs} dofs; "
                f"the elements need {offset}"
            )
        points = np.array(xyz, dtype=np.float64).reshape(len(point_vertex), pdim)
    else:
        n = len(numbering)
        points = np.empty((n, pdim), dtype=np.float64)
        for k in range(n):
            if coords == "h1" and mesh_order == 2:
                for c in range(pdim):
                    points[k, c] = _value(nodes, nspace, node_dofs, k, c)
            else:
                key = numbering.keys[k]
                for c in range(pdim):
                    s = 0.0
                    for v in key:
                        s += vxyz[v, c]
                    points[k, c] = s / len(key)
    npts = len(points)

    def type_of(geom):
        g = _GEOMS[geom]
        return g[1] if order == 2 else g[0]

    blocks = []  # (type, [(source, is_boundary)])

    def add_group(items, is_boundary):
        index = {}
        for k, el in enumerate(items):
            t = type_of(el[1])
            if t not in index:
                index[t] = len(blocks)
                blocks.append((t, []))
            blocks[index[t]][1].append((k, is_boundary))

    add_group(elements, False)

    vertex_elements = {}
    if coords == "dg":
        for e, el in enumerate(elements):
            for v in el[2]:
                vertex_elements.setdefault(v, []).append(e)
    boundary = []
    boundary_dg = []
    orphans = 0
    for b in f["boundary"]:
        if coords == "dg":
            pts = []
            for e in vertex_elements.get(b[2][0], []):
                verts = elements[e][2]
                pts = []
                for v in b[2]:
                    if v not in verts:
                        break
                    pts.append(dg_points[e][verts.index(v)])
                if len(pts) == len(b[2]):
                    break
            if len(pts) != len(b[2]):
                orphans += 1
                continue
            boundary_dg.append(pts)
        boundary.append(b)
    if orphans:
        warn(f"MFEM mesh: {orphans} boundary element(s) lie on no element; dropped")
    add_group(boundary, True)

    unresolved = 0
    cells = []
    attr_blocks = []
    cell_attr = []
    cell_is_boundary = []
    element_cell = [0] * len(elements)
    for cell_type, members in blocks:
        slots = _SLOTS.get(cell_type, [])
        rows = []
        attrs = []
        for source, is_boundary in members:
            attr, geom, verts, _ = boundary[source] if is_boundary else elements[source]
            nc = _GEOMS[geom][3]
            corner = []
            row = []
            for j in range(nc):
                src = _vtk_corner(geom, j)
                corner.append(verts[src])
                if coords == "dg":
                    own = boundary_dg[source] if is_boundary else dg_points[source]
                    row.append(own[src])
                else:
                    row.append(verts[src])
            element = -1 if is_boundary else source
            for slot in slots:
                sv = [corner[q] for q in slot]
                dof = numbering.slot_dof(dim, sv, nc, element)
                if dof < 0:
                    unresolved += 1
                    dof = sv[0]
                row.append(dof)
            rows.append(row)
            attrs.append(attr)
            cell_attr.append(attr)
            cell_is_boundary.append(is_boundary)
            if not is_boundary:
                element_cell[source] = len(cell_attr) - 1
        cells.append((cell_type, np.array(rows, dtype=np.int64).reshape(len(rows), -1)))
        attr_blocks.append(np.array(attrs, dtype=np.int64))
    if unresolved:
        warn(f"MFEM mesh: {unresolved} boundary node(s) lie on no element edge or face")

    mesh = Mesh(points, cells)
    if attr_blocks:
        mesh.cell_data["mfem:attribute"] = attr_blocks
    ncells = len(cell_attr)

    for name, space, values in gfs:
        vdim = space.vdim
        if len(values) % vdim:
            raise ReadError(
                f"MFEM grid function '{name}': {len(values)} values, not a multiple of "
                f"VDim {vdim}"
            )
        ndofs = len(values) // vdim
        if space.kind == "h1":
            expect = len(numbering) if space.order == 2 else nv
            if ndofs != expect:
                raise ReadError(
                    f"MFEM grid function '{name}' has {ndofs} dofs; the mesh has "
                    f"{expect} at order {space.order}"
                )
            data = np.empty((npts, vdim), dtype=np.float64)
            for p in range(npts):
                for c in range(vdim):
                    if coords == "dg":
                        data[p, c] = _value(values, space, ndofs, point_vertex[p], c)
                    elif space.order == 2:
                        data[p, c] = _value(values, space, ndofs, p, c)
                    else:
                        key = numbering.keys[p]
                        s = 0.0
                        for v in key:
                            s += _value(values, space, ndofs, v, c)
                        data[p, c] = s / len(key)
            mesh.point_data[name] = data[:, 0].copy() if vdim == 1 else data
        else:
            if ndofs != len(elements):
                raise ReadError(
                    f"MFEM grid function '{name}' has {ndofs} dofs for "
                    f"{len(elements)} elements"
                )
            per_cell = np.full((ncells, vdim), np.nan)
            for e in range(len(elements)):
                for c in range(vdim):
                    per_cell[element_cell[e], c] = _value(values, space, ndofs, e, c)
            out = []
            start = 0
            for _, members in blocks:
                a = per_cell[start : start + len(members)]
                start += len(members)
                out.append(a[:, 0].copy() if vdim == 1 else a.copy())
            mesh.cell_data[name] = out

    by_attr = {}
    for g in range(ncells):
        by_attr.setdefault((cell_is_boundary[g], cell_attr[g]), []).append(g)
    regions = []
    for (bdr, a), ids in sorted(by_attr.items()):
        regions.append(
            Region(
                f"{'boundary_' if bdr else 'attribute_'}{a}",
                "cell",
                np.array(ids, dtype=np.int64),
                dim - 1 if bdr else dim,
                a,
            )
        )
    for sets, bdr in ((f["sets"], False), (f["bdr_sets"], True)):
        for name, attrs in sets:
            ids = []
            for a in attrs:
                ids += by_attr.get((bdr, a), [])
            regions.append(
                Region(
                    name,
                    "cell",
                    np.array(ids, dtype=np.int64),
                    dim - 1 if bdr else dim,
                    -1,
                )
            )
    mesh.regions = regions
    return mesh


# --- writer -------------------------------------------------------------------


def _fmt(v):
    return "%.17g" % float(v)


def _centre(cell_type, geom, nodes, slot):
    """Point weights of the node between ``slot`` corners that a cell lacks."""
    n = len(slot)
    full = {"quad8": "quad9", "hexahedron20": "hexahedron27", "wedge15": "wedge18"}

    def mid_of(a, b):
        slots = _SLOTS[full.get(cell_type, cell_type)]
        nc = _GEOMS[geom][3]
        for s, sl in enumerate(slots):
            if len(sl) == 2 and set(sl) == {a, b}:
                return nc + s
        return -1

    serendipity = cell_type in full
    if serendipity and n == 4:
        w = []
        for k in range(4):
            w.append((nodes[slot[k]], -0.25))
            w.append((nodes[mid_of(slot[k], slot[(k + 1) % 4])], 0.5))
        return w
    if cell_type == "hexahedron20" and n == 8:
        return [(nodes[k], -0.25) for k in range(8)] + [
            (nodes[k], 0.25) for k in range(8, 20)
        ]
    return [(nodes[k], 1.0 / n) for k in slot]


def _sanitise(name):
    out = "".join(
        c if (c.isascii() and (c.isalnum() or c in "_-.")) else "_" for c in name
    )
    return out or "data"


def _quote(name):
    out = []
    for c in name:
        if c in '"\\':
            out.append("\\")
        out.append(" " if c in "\n\r" else c)
    return '"' + "".join(out) + '"'


def _generated_id(name, prefix):
    if not name.startswith(prefix):
        return None
    rest = name[len(prefix) :]
    v = _parse_int(rest)
    if v is None or v <= 0 or rest != str(v):
        return None
    return v


def write(filename, mesh, grid_functions=False):
    points = np.asarray(mesh.points)
    pdim = points.shape[1] if points.ndim == 2 else 0
    if pdim < 1 or pdim > 3:
        raise WriteError(f"MFEM mesh writer: points of dimension {pdim} (1, 2 or 3)")

    dim = -1
    block_geom = []
    for block in mesh.cells:
        ragged = not isinstance(block.data, np.ndarray)
        g = -1 if ragged else _geom_of_type(block.type)
        known = g >= 0 and (
            block.type in (_GEOMS[g][0], _GEOMS[g][1]) or block.type in _WRITABLE_EXTRA
        )
        block_geom.append(g if known else -1)
        if known and len(block.data):
            dim = max(dim, _GEOMS[g][2])
    if dim < 1:
        raise WriteError("MFEM mesh writer: no cells MFEM can hold as elements")
    if pdim < dim:
        raise WriteError(f"MFEM mesh writer: {pdim}-D points for {dim}-D cells")
    dropped = set()
    for b, block in enumerate(mesh.cells):
        g = block_geom[b]
        if g < 0 or _GEOMS[g][2] not in (dim, dim - 1):
            if len(block.data):
                dropped.add(block.type)
            block_geom[b] = -1

    attr_data = mesh.cell_data.get("mfem:attribute")
    has_attr = attr_data is not None
    starts = [0]
    for block in mesh.cells:
        starts.append(starts[-1] + len(block.data))
    ncells = starts[-1]
    cell_dim = [-1] * ncells
    for b in range(len(mesh.cells)):
        if block_geom[b] >= 0:
            for g in range(starts[b], starts[b + 1]):
                cell_dim[g] = _GEOMS[block_geom[b]][2]
    attr = [0] * ncells
    if has_attr:
        for b in range(len(mesh.cells)):
            a = np.asarray(attr_data[b]).ravel()
            for r in range(len(mesh.cells[b].data)):
                attr[starts[b] + r] = int(a[r])
    used = [set(), set()]
    for g in range(ncells):
        if attr[g] > 0 and cell_dim[g] >= 0:
            used[0 if cell_dim[g] == dim else 1].add(attr[g])
    named = [[], []]  # [name, cells]
    next_attr = [1, 1]

    def fresh(which):
        while next_attr[which] in used[which]:
            next_attr[which] += 1
        used[which].add(next_attr[which])
        return next_attr[which]

    sides = []  # [attribute, facets, set index]
    regions = sorted(getattr(mesh, "regions", []) or [], key=lambda r: r.key)
    for reg in regions:
        if reg.kind == "point":
            continue
        gid = _generated_id(reg.name, "attribute_")
        if gid is None:
            gid = _generated_id(reg.name, "boundary_")
        generated = gid is not None
        entries = np.asarray(reg.entries, dtype=np.int64)
        if reg.kind == "side":
            a = gid if generated else reg.tag
            if a <= 0 or a in used[1]:
                a = 0
            else:
                used[1].add(a)
            set_index = -1
            if not generated:
                set_index = len(named[1])
                named[1].append([reg.name, []])
            sides.append(
                [a, [(int(c), int(f)) for c, f in entries.reshape(-1, 2)], set_index]
            )
            continue
        cells = [int(c) for c in entries.ravel()]
        if not has_attr:
            for which in (0, 1):
                want = dim if which == 0 else dim - 1
                if not any(
                    c < ncells and cell_dim[c] == want and attr[c] == 0 for c in cells
                ):
                    continue
                a = gid if generated else reg.tag
                if a <= 0 or a in used[which]:
                    a = fresh(which)
                else:
                    used[which].add(a)
                for c in cells:
                    if c < ncells and cell_dim[c] == want and attr[c] == 0:
                        attr[c] = a
        if not generated:
            el = [c for c in cells if c < ncells and cell_dim[c] == dim]
            bd = [c for c in cells if c < ncells and cell_dim[c] == dim - 1]
            if el:
                named[0].append([reg.name, el])
            if bd:
                named[1].append([reg.name, bd])
    for g in range(ncells):
        if cell_dim[g] >= 0 and attr[g] <= 0:
            attr[g] = 1
            used[0 if cell_dim[g] == dim else 1].add(1)

    elements = []  # (geom, attribute, cell, nodes, type)
    boundary = []
    quadratic = False
    has_pyramid = False
    for b, block in enumerate(mesh.cells):
        g = block_geom[b]
        if g < 0:
            continue
        data = np.asarray(block.data, dtype=np.int64)
        k = data.shape[1] if data.ndim == 2 else 0
        quadratic = quadratic or k > _GEOMS[g][3]
        has_pyramid = has_pyramid or g == 7
        target = elements if _GEOMS[g][2] == dim else boundary
        for r in range(len(data)):
            target.append(
                (
                    g,
                    attr[starts[b] + r],
                    starts[b] + r,
                    [int(v) for v in data[r]],
                    block.type,
                )
            )
    bad_facets = 0
    for side in sides:
        if side[0] <= 0:
            side[0] = fresh(1)
        if side[2] >= 0:
            named[1][side[2]][1] = [-side[0]]
        for cell, facet in side[1]:
            found = facet_nodes(mesh, cell, facet)
            if found is None:
                bad_facets += 1
                continue
            ftype, nodes = found
            g = _geom_of_type(ftype)
            if g < 0 or _GEOMS[g][2] != dim - 1:
                bad_facets += 1
                continue
            quadratic = quadratic or len(nodes) > _GEOMS[g][3]
            boundary.append((g, side[0], -1, list(nodes), ftype))
    if quadratic and has_pyramid:
        warn(
            "MFEM mesh writer: MFEM order-2 meshes hold no pyramids; writing corners only"
        )
        _provenance.note(
            "high-order-dropped",
            "a mesh with pyramids is written with linear MFEM cells",
        )
        quadratic = False
    for t in sorted(dropped):
        warn(
            f"MFEM mesh writer: '{t}' cells are neither elements nor boundary; dropped"
        )
        _provenance.note("cells-dropped", f"'{t}' cells have no MFEM equivalent here")
    point_regions = sum(1 for r in regions if r.kind == "point")
    if point_regions:
        warn(
            f"MFEM mesh writer: MFEM has no node sets; {point_regions} point "
            "region(s) dropped"
        )
        _provenance.note(
            "regions-dropped",
            f"{point_regions} point region(s) have no MFEM equivalent",
        )
    if bad_facets:
        warn(
            f"MFEM mesh writer: {bad_facets} side region entr(ies) name no facet and "
            "were dropped"
        )
        _provenance.note(
            "regions-dropped", f"{bad_facets} side region entries name no facet"
        )

    npts = len(points)
    vertex_of = [-1] * npts
    vertex_point = []
    if not quadratic and all(len(c[3]) <= _GEOMS[c[0]][3] for c in elements + boundary):
        vertex_of = list(range(npts))
        vertex_point = list(range(npts))
    if not vertex_point:
        is_corner = [False] * npts
        for c in elements + boundary:
            for j in range(_GEOMS[c[0]][3]):
                is_corner[c[3][j]] = True
        for p in range(npts):
            if is_corner[p]:
                vertex_of[p] = len(vertex_point)
                vertex_point.append(p)

    def mfem_vertices(c):
        nc = _GEOMS[c[0]][3]
        v = [0] * nc
        for j in range(nc):
            v[_vtk_corner(c[0], j)] = vertex_of[c[3][j]]
        return v

    sets = [[], []]
    for which in (0, 1):
        for name, cells in named[which]:
            attrs = sorted({-c if c < 0 else attr[c] for c in cells})
            if attrs:
                sets[which].append((name, attrs))
    v13 = bool(sets[0] or sets[1])

    nv = len(vertex_point)
    numbering = None
    weights = []
    if quadratic:
        numbered = [(c[1], c[0], mfem_vertices(c), 0) for c in elements]
        numbering = _Numbering(numbered, dim, nv)
        weights = [[] for _ in range(len(numbering))]
        exact = [False] * len(numbering)
        for v in range(nv):
            weights[v] = [(vertex_point[v], 1.0)]
            exact[v] = True

        def visit(c, element):
            geom, _, _, nodes, cell_type = c
            full = _GEOMS[geom][1]
            nc = _GEOMS[geom][3]
            corner = [vertex_of[nodes[j]] for j in range(nc)]
            for s, slot in enumerate(_SLOTS[full]):
                sv = [corner[q] for q in slot]
                dof = numbering.slot_dof(dim, sv, nc, element)
                if dof < 0 or exact[dof]:
                    continue
                if nc + s < len(nodes):
                    weights[dof] = [(nodes[nc + s], 1.0)]
                    exact[dof] = True
                elif not weights[dof] or len(nodes) > nc:
                    weights[dof] = _centre(cell_type, geom, nodes, slot)
                    if len(nodes) > nc:
                        exact[dof] = True

        for e, c in enumerate(elements):
            visit(c, e)
        for c in boundary:
            visit(c, -1)
        for d in range(len(numbering)):
            if not weights[d]:
                key = numbering.keys[d]
                weights[d] = [(vertex_point[v], 1.0 / len(key)) for v in key]

    data_arrays = (
        len(mesh.point_data)
        + len(mesh.cell_data)
        - (1 if has_attr else 0)
        + len(mesh.field_data)
    )
    if not grid_functions and data_arrays:
        warn(
            "MFEM mesh writer: data arrays are dropped (write grid functions to keep "
            "them)"
        )
        _provenance.note("data-dropped", "an MFEM mesh holds no data arrays")
    elif grid_functions and mesh.field_data:
        warn("MFEM mesh writer: field data has no grid function; dropped")
        _provenance.note("data-dropped", "field data has no MFEM grid function")

    pts = points.astype(np.float64, copy=False)

    def evaluate(w, data, cols, c):
        s = 0.0
        for p, wt in w:
            s += wt * float(data[p * cols + c])
        return s

    out = ["MFEM mesh v1.3\n" if v13 else "MFEM mesh v1.0\n"]
    out.append(_provenance.render_lines(_provenance.SlotTier.BLOCK, "# "))
    out.append(f"\ndimension\n{dim}\n\nelements\n{len(elements)}\n")

    def append_cells(items):
        for c in items:
            out.append(
                f"{c[1]} {c[0]}" + "".join(f" {v}" for v in mfem_vertices(c)) + "\n"
            )

    def append_sets(items):
        out.append(f"{len(items)}\n")
        for name, attrs in items:
            out.append(
                _quote(name) + f" {len(attrs)}" + "".join(f" {a}" for a in attrs) + "\n"
            )

    append_cells(elements)
    if v13:
        out.append("\nattribute_sets\n")
        append_sets(sets[0])
    out.append(f"\nboundary\n{len(boundary)}\n")
    append_cells(boundary)
    if v13:
        out.append("\nbdr_attribute_sets\n")
        append_sets(sets[1])
    out.append(f"\nvertices\n{nv}\n")
    flat = pts.ravel()
    if not quadratic:
        out.append(f"{pdim}\n")
        for v in range(nv):
            out.append(
                " ".join(_fmt(flat[vertex_point[v] * pdim + c]) for c in range(pdim))
                + "\n"
            )
    else:
        out.append(
            f"\nnodes\nFiniteElementSpace\nFiniteElementCollection: H1_{dim}D_P2\n"
            f"VDim: {pdim}\nOrdering: 1\n\n"
        )
        for d in range(len(numbering)):
            out.append(
                " ".join(_fmt(evaluate(weights[d], flat, pdim, c)) for c in range(pdim))
                + "\n"
            )
    if v13:
        out.append("\nmfem_mesh_end\n")
    with open_file(filename, "w", newline="\n") as fh:
        fh.write("".join(out))
    if not grid_functions:
        return

    path = pathlib.Path(filename)
    stem = os.path.join(os.path.dirname(str(filename)), path.stem)
    order = 2 if quadratic else 1
    ndofs = len(numbering) if quadratic else nv
    for name in sorted(mesh.point_data):
        a = np.asarray(mesh.point_data[name])
        cols = a.shape[1] if a.ndim > 1 else 1
        flat_a = a.astype(np.float64).ravel()
        text = [
            f"FiniteElementSpace\nFiniteElementCollection: H1_{dim}D_P{order}\n"
            f"VDim: {cols}\nOrdering: 1\n\n"
        ]
        for d in range(ndofs):
            if quadratic:
                vals = [evaluate(weights[d], flat_a, cols, c) for c in range(cols)]
            else:
                vals = [flat_a[vertex_point[d] * cols + c] for c in range(cols)]
            text.append(" ".join(_fmt(v) for v in vals) + "\n")
        with open(f"{stem}.{_sanitise(name)}.gf", "w", newline="\n") as fh:
            fh.write("".join(text))
    for name in sorted(mesh.cell_data):
        if name == "mfem:attribute":
            continue
        arrays = [np.asarray(x) for x in mesh.cell_data[name]]
        cols = 1
        for x in arrays:
            if x.ndim > 1:
                cols = x.shape[1]
        text = [
            f"FiniteElementSpace\nFiniteElementCollection: L2_{dim}D_P0\n"
            f"VDim: {cols}\nOrdering: 1\n\n"
        ]
        for c in elements:
            cell = c[2]
            b = 0
            while cell >= starts[b + 1]:
                b += 1
            row = arrays[b].astype(np.float64).ravel()
            r = cell - starts[b]
            text.append(" ".join(_fmt(row[r * cols + k]) for k in range(cols)) + "\n")
        with open(f"{stem}.{_sanitise(name)}.gf", "w", newline="\n") as fh:
            fh.write("".join(text))

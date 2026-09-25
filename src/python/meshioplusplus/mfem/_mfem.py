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

import math
import os
import pathlib

import numpy as np

from .. import _lagrange, _provenance
from .._common import warn
from .._exceptions import ReadError, WriteError
from .._facets import facet_nodes
from .._files import open_file
from .._mesh import Mesh
from .._regions import Region
from . import _nurbs

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


_LAGRANGE_GEOMS = {
    "VTK_LAGRANGE_CURVE": 1,
    "VTK_LAGRANGE_TRIANGLE": 2,
    "VTK_LAGRANGE_QUADRILATERAL": 3,
    "VTK_LAGRANGE_TETRAHEDRON": 4,
    "VTK_LAGRANGE_HEXAHEDRON": 5,
    "VTK_LAGRANGE_WEDGE": 6,
}


def _geom_of_type(cell_type):
    # VTK Lagrange cells of any order: written as order-p H1 nodes.
    if cell_type in _LAGRANGE_GEOMS:
        return _LAGRANGE_GEOMS[cell_type]
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
        # an H1 space's nodes: "gll" (the default), "uniform" (H1@U) or "cubic"
        # (the legacy Cubic collection: equispaced, its own hex interior);
        # "bernstein" (H1Pos: Bernstein coefficients on the uniform lattice)
        # and "serendipity" (H1Ser: MFEM's serendipity quadrilaterals) are
        # modal and evaluated through their own basis
        self.points = "gll"
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
        elif name == "Cubic":
            self.kind, self.order, self.points = "h1", 3, "cubic"
        elif name == "QuadraticPos":
            self.kind, self.order = "h1-other", 2
        elif name.startswith(("H1_", "H1@")):
            self.order = order_after_p()
            basis = name[3] if name[2] == "@" and len(name) > 3 else "G"
            self.points = "uniform" if basis == "U" else "gll"
            nodal = basis in "GU" or self.order <= 2
            self.kind = "h1" if self.order >= 1 and nodal else "h1-other"
        elif name.startswith(("H1Pos_", "H1Ser_")):
            self.order = order_after_p()
            self.kind = "h1" if self.order >= 1 else "h1-other"
            if self.order >= 2:
                self.points = "bernstein" if name[2] == "P" else "serendipity"
        elif name.startswith("NURBS"):
            # NURBS<p>, or NURBS alone for the orders of the mesh's knot vectors
            self.kind = "nurbs"
            self.order = _parse_int(name[5:]) if len(name) > 5 else -1
            if self.order is None:
                self.kind = "other"
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


def _read_groups(lex, f):
    """``mf_read_groups``: a parallel mesh's communication groups and their
    shared vertices, after ``mfem_serial_mesh_end``."""
    f["parallel"] = True
    text, line, _ = lex.next("communication_groups")
    if text != "communication_groups":
        lex.fail(f"expected communication_groups, found '{text}'", line)
    if lex.next("number_of_groups")[0] != "number_of_groups":
        lex.fail("expected number_of_groups", line)
    ngroups = lex.int("a group count")
    if ngroups < 1:
        lex.fail("a parallel mesh needs at least one communication group", line)
    for _ in range(ngroups):
        size = lex.int("a group size")
        if size < 1:
            lex.fail("an empty communication group", lex.line())
        f["groups"].append(tuple(sorted(lex.int("a rank") for _ in range(size))))
    if len(f["groups"][0]) != 1:
        lex.fail("communication group 0 must hold this rank alone", line)
    f["rank"] = f["groups"][0][0]
    f["group_vertices"] = [[] for _ in range(ngroups)]
    group = 0
    while not lex.at_end():
        text, tline, _ = lex.next("a shared-entity section")
        if text == "mfem_mesh_end":
            break
        if text in (
            "total_shared_vertices",
            "total_shared_edges",
            "total_shared_faces",
        ):
            lex.int("a count")
        elif text == "shared_vertices":
            group += 1
            if group >= ngroups:
                lex.fail("more shared-vertex groups than communication groups", tline)
            n = lex.int("a vertex count")
            f["group_vertices"][group] = [lex.int("a shared vertex") for _ in range(n)]
        elif text == "shared_edges":
            for _ in range(2 * lex.int("an edge count")):
                lex.int("an edge vertex")
        elif text == "shared_faces":
            for _ in range(lex.int("a face count")):
                geom = lex.int("a face geometry")
                for _ in range(4 if geom == 3 else 3):
                    lex.int("a face vertex")
        else:
            lex.fail(f"unexpected '{text}' in the communication groups", tline)


# MFEM's Hilbert-curve child orders and states (mesh/ncmesh_tables.hpp)
_QUAD_HILBERT_ORDER = [
    (0, 1, 2, 3),
    (0, 3, 2, 1),
    (1, 2, 3, 0),
    (1, 0, 3, 2),
    (2, 3, 0, 1),
    (2, 1, 0, 3),
    (3, 0, 1, 2),
    (3, 2, 1, 0),
]
_QUAD_HILBERT_STATE = [
    (1, 0, 0, 5),
    (0, 1, 1, 4),
    (3, 2, 2, 7),
    (2, 3, 3, 6),
    (5, 4, 4, 1),
    (4, 5, 5, 0),
    (7, 6, 6, 3),
    (6, 7, 7, 2),
]
_HEX_HILBERT_ORDER = [
    (0, 1, 2, 3, 7, 6, 5, 4),
    (0, 3, 7, 4, 5, 6, 2, 1),
    (0, 4, 5, 1, 2, 6, 7, 3),
    (1, 0, 3, 2, 6, 7, 4, 5),
    (1, 2, 6, 5, 4, 7, 3, 0),
    (1, 5, 4, 0, 3, 7, 6, 2),
    (2, 1, 5, 6, 7, 4, 0, 3),
    (2, 3, 0, 1, 5, 4, 7, 6),
    (2, 6, 7, 3, 0, 4, 5, 1),
    (3, 0, 4, 7, 6, 5, 1, 2),
    (3, 2, 1, 0, 4, 5, 6, 7),
    (3, 7, 6, 2, 1, 5, 4, 0),
    (4, 0, 1, 5, 6, 2, 3, 7),
    (4, 5, 6, 7, 3, 2, 1, 0),
    (4, 7, 3, 0, 1, 2, 6, 5),
    (5, 1, 0, 4, 7, 3, 2, 6),
    (5, 4, 7, 6, 2, 3, 0, 1),
    (5, 6, 2, 1, 0, 3, 7, 4),
    (6, 2, 3, 7, 4, 0, 1, 5),
    (6, 5, 1, 2, 3, 0, 4, 7),
    (6, 7, 4, 5, 1, 0, 3, 2),
    (7, 3, 2, 6, 5, 1, 0, 4),
    (7, 4, 0, 3, 2, 1, 5, 6),
    (7, 6, 5, 4, 0, 1, 2, 3),
]
_HEX_HILBERT_STATE = [
    (1, 2, 2, 7, 7, 21, 21, 17),
    (2, 0, 0, 22, 22, 16, 16, 8),
    (0, 1, 1, 15, 15, 6, 6, 23),
    (4, 5, 5, 10, 10, 18, 18, 14),
    (5, 3, 3, 19, 19, 13, 13, 11),
    (3, 4, 4, 12, 12, 9, 9, 20),
    (8, 7, 7, 17, 17, 23, 23, 2),
    (6, 8, 8, 0, 0, 15, 15, 22),
    (7, 6, 6, 21, 21, 1, 1, 16),
    (11, 10, 10, 14, 14, 20, 20, 5),
    (9, 11, 11, 3, 3, 12, 12, 19),
    (10, 9, 9, 18, 18, 4, 4, 13),
    (13, 14, 14, 5, 5, 19, 19, 10),
    (14, 12, 12, 20, 20, 11, 11, 4),
    (12, 13, 13, 9, 9, 3, 3, 18),
    (16, 17, 17, 2, 2, 22, 22, 7),
    (17, 15, 15, 23, 23, 8, 8, 1),
    (15, 16, 16, 6, 6, 0, 0, 21),
    (20, 19, 19, 11, 11, 14, 14, 3),
    (18, 20, 20, 4, 4, 10, 10, 12),
    (19, 18, 18, 13, 13, 5, 5, 9),
    (23, 22, 22, 8, 8, 17, 17, 0),
    (21, 23, 23, 1, 1, 7, 7, 15),
    (22, 21, 21, 16, 16, 2, 2, 6),
]


def _facets(geom):
    """The local vertex lists of an element's facets (edges in 2-D, faces in
    3-D; the vertices of a segment)."""
    if geom == 1:
        return [(0,), (1,)]
    return list(_GEOMS[geom][4 if _GEOMS[geom][2] == 2 else 5])


def _parse_nc(lex, filename, scaled):
    """An ``MFEM NC mesh`` (``mf_parse_nc``): the refinement tree read as its
    leaves; top-level vertices from ``coordinates``, the rest between their
    ``vertex_parents``; only the leaves of the file's own rank."""
    f = {
        "nc": True,
        "interface": [],
        "parallel": False,
        "rank": 0,
        "groups": [],
        "group_vertices": [],
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
    elements = []  # (rank, attr, geom, ref_type, ids, line); geom None when unused
    root_states = []
    parents = {}
    top = []
    my_rank = 0
    sdim = 0
    have_coordinates = False
    while not lex.at_end():
        text, line, _ = lex.next("a section")
        if text == "dimension":
            d = lex.int("a dimension")
            if d < 1 or d > 3:
                lex.fail(f"dimension {d} (1, 2 or 3)", line)
            f["dim"] = d
        elif text == "rank":
            my_rank = lex.int("a rank")
        elif text == "sfc_version":
            lex.int("an SFC version")
        elif text == "elements":
            n = lex.int("an element count")
            if n < 0:
                lex.fail("negative element count", line)
            for _ in range(n):
                row_line = lex.line()
                rank = lex.int("a rank")
                attr = lex.int("an attribute")
                geom = lex.int("a geometry type")
                if geom == -1:
                    elements.append((rank, attr, None, 0, [], row_line))
                    continue
                if not 1 <= geom <= 7:
                    lex.fail(f"unknown geometry type {geom}", row_line)
                ref = lex.int("a refinement type")
                ids = []
                while not lex.at_end() and lex.line() == row_line:
                    ids.append(lex.int("a node or child"))
                if ref == 0 and len(ids) != _GEOMS[geom][3]:
                    lex.fail(f"a leaf element lists {len(ids)} vertices", row_line)
                elements.append((rank, attr, geom, ref, ids, row_line))
        elif text == "boundary":
            f["boundary"] = _read_elements(lex, "boundary element")
        elif text == "vertex_parents":
            for _ in range(lex.int("a vertex count")):
                vid = lex.int("a vertex")
                p1 = lex.int("a parent")
                p2 = lex.int("a parent")
                parents[vid] = (p1, p2, lex.real("a scale") if scaled else 0.5)
        elif text == "root_state":
            root_states = [
                lex.int("a root state") for _ in range(lex.int("a root count"))
            ]
        elif text == "coordinates":
            n = lex.int("a vertex count")
            if n < 0:
                lex.fail("negative vertex count", line)
            if n > 0:
                sdim = lex.int("a space dimension")
                if not 1 <= sdim <= 3:
                    lex.fail(f"space dimension {sdim} (1, 2 or 3)", line)
                for _ in range(n):
                    x = [lex.real("a coordinate") for _ in range(sdim)]
                    top.append(x + [0.0] * (3 - sdim))
            have_coordinates = True
        elif text == "nodes":
            lex.fail(
                "curved non-conforming meshes (a 'nodes' section) are not supported",
                line,
            )
        elif text in ("mfem_mesh_end", "mfem_serial_mesh_end"):
            break
        else:
            lex.fail(f"unexpected '{text}'", line)
    if f["dim"] < 0:
        raise ReadError(f"MFEM mesh: no dimension section in {filename}")
    if not have_coordinates:
        raise ReadError(
            f"MFEM mesh: the non-conforming mesh {filename} has no top-level coordinates"
        )
    is_child = [False] * len(elements)
    for _, _, geom, ref, ids, row_line in elements:
        if geom is not None and ref:
            for c in ids:
                if not 0 <= c < len(elements):
                    raise ReadError(
                        f"MFEM mesh: child element {c} out of range (line {row_line})"
                    )
                is_child[c] = True
    # MFEM's leaf order (NCMesh::CollectLeafElements): the roots in order,
    # children along its Hilbert curve for quadrilaterals refined in both
    # directions and hexahedra in all three, else in child order; a file's
    # roots start in their root_state (0 by default).
    ordered, ghosts = [], 0
    roots = [r for r in range(len(elements)) if elements[r][2] and not is_child[r]]
    stack = [
        (r, root_states[k] if k < len(root_states) else 0)
        for k, r in reversed(list(enumerate(roots)))
    ]
    seen = [False] * len(elements)
    while stack:
        e, state = stack.pop()
        if seen[e]:
            raise ReadError(
                f"MFEM mesh: element {e} is reached twice in the refinement tree"
            )
        seen[e] = True
        rank, _, geom, ref, ids, _ = elements[e]
        if ref == 0:
            if rank >= 0:
                ordered.append(e)
            continue
        if geom == 3 and ref == 3 and 0 <= state < 8:
            kids = [
                (ids[_QUAD_HILBERT_ORDER[state][i]], _QUAD_HILBERT_STATE[state][i])
                for i in range(4)
            ]
        elif geom == 5 and ref == 7 and 0 <= state < 24:
            kids = [
                (ids[_HEX_HILBERT_ORDER[state][i]], _HEX_HILBERT_STATE[state][i])
                for i in range(8)
            ]
        else:
            kids = [(c, state) for c in ids]
        stack.extend(reversed(kids))
    leaves = [e for e in ordered if elements[e][0] == my_rank]
    ghosts = len(ordered) - len(leaves)
    if ghosts:
        warn(
            f"MFEM mesh: {ghosts} ghost element(s) of other ranks in {filename} dropped"
        )
    if ghosts:  # a rank of a parallel mesh: only the boundary of its own leaves
        faces = set()
        for e in leaves:
            geom, ids = elements[e][2], elements[e][4]
            for fv in _facets(geom):
                faces.add(tuple(sorted(ids[k] for k in fv)))
        f["boundary"] = [b for b in f["boundary"] if tuple(sorted(b[2])) in faces]
    # MFEM's vertex numbers (NCMesh::UpdateVertices): the top-level vertices
    # of the rank's leaves by node id, then the others as the leaves (ghosts
    # included) meet them
    local = {}
    for e in leaves:
        for v in elements[e][4]:
            local[v] = local.get(v, False) or v not in parents
    order = sorted(v for v, top_level in local.items() if top_level)
    numbered = set(order)
    for e in ordered:
        for v in elements[e][4]:
            if v in local and v not in numbered:
                numbered.add(v)
                order.append(v)
    for b in f["boundary"]:
        for v in b[2]:
            if v not in numbered:
                numbered.add(v)
                order.append(v)
    index = {vid: k for k, vid in enumerate(order)}
    # A rank of a parallel mesh: the vertices its ghosts share with it (where it
    # meets its neighbours).
    ghost_ids = {
        v for e in ordered if elements[e][0] != my_rank for v in elements[e][4]
    }
    f["interface"] = sorted(index[v] for v in ghost_ids if v in local)
    f["rank"] = my_rank
    pos = {}

    def position(vid, visiting=()):
        if vid in pos:
            return pos[vid]
        par = parents.get(vid)
        if par is None:
            if not 0 <= vid < len(top):
                raise ReadError(
                    f"MFEM mesh: vertex {vid} has neither coordinates nor parents"
                )
            x = top[vid]
        else:
            if vid in visiting:
                raise ReadError(f"MFEM mesh: cyclic vertex parents at vertex {vid}")
            a = position(par[0], visiting + (vid,))
            b = position(par[1], visiting + (vid,))
            s = par[2]
            x = [(1.0 - s) * a[c] + s * b[c] for c in range(3)]
        pos[vid] = x
        return x

    f["sdim"] = sdim if sdim else f["dim"]
    f["nv"] = len(index)
    for vid in order:
        f["coords"].extend(position(vid)[: f["sdim"]])
    f["elements"] = [
        (
            elements[e][1],
            elements[e][2],
            [index[v] for v in elements[e][4]],
            elements[e][5],
        )
        for e in leaves
    ]
    f["boundary"] = [
        (a, g, [index[v] for v in verts], ln) for a, g, verts, ln in f["boundary"]
    ]
    warn(
        f"MFEM mesh: the non-conforming mesh {filename} is read as its {len(leaves)} leaf "
        "element(s); hanging nodes are left unconstrained"
    )
    return f


def _parse(filename):
    lex = _Lexer("MFEM mesh", _read_text(filename))
    header = lex.header or ""
    if header in ("MFEM NC mesh v1.0", "MFEM NC mesh v1.1"):
        return _parse_nc(lex, filename, header == "MFEM NC mesh v1.1")
    if header.startswith("MFEM NC mesh"):
        lex.fail(
            f"non-conforming mesh version '{header}' is not supported", lex.header_line
        )
    if header in ("MFEM NURBS mesh v1.0", "MFEM NURBS mesh v1.1"):
        return {"nurbs": _nurbs.Nurbs(lex, filename, header), "parallel": False}
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
        "nc": False,
        "parallel": False,
        "rank": 0,
        "groups": [],
        "group_vertices": [],
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
            # the legacy non-conforming layout: the leaf mesh, then how it
            # refines; read as its leaves
            f["nc"] = True
            for _ in range(lex.int("a count")):
                row_line = lex.line()
                while not lex.at_end() and lex.line() == row_line:
                    lex.next("a value")
        elif text == "mfem_serial_mesh_end":
            _read_groups(lex, f)
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


def _read_space_end(lex):
    """The ``End: MFEM FiniteElementSpace v1.0`` closing the versioned header
    MFEM writes for NURBS spaces (its variable-order element lists are not
    read)."""
    text, line, _ = lex.next("'End:'")
    if text != "End:":
        lex.fail(f"'{text}' in a versioned FiniteElementSpace (not supported)", line)
    for word in ("MFEM", "FiniteElementSpace", "v1.0"):
        text, line, _ = lex.next(f"'{word}'")
        if text != word:
            lex.fail(f"expected '{word}', found '{text}'", line)


def _parse_gf(name, path):
    lex = _Lexer("MFEM grid function", "\n" + _read_text(path))
    versioned = lex.header == "MFEM FiniteElementSpace v1.0"
    if lex.header != "FiniteElementSpace" and not versioned:
        raise ReadError(
            f"MFEM grid function: '{path}' does not start with FiniteElementSpace "
            "(variable-order spaces are not supported)"
        )
    space = _read_space_body(lex, 1)
    if versioned and space is not None:
        _read_space_end(lex)
    if space is None:
        raise ReadError(
            f"MFEM grid function: '{path}' names no FiniteElementCollection"
        )
    values = lex.reals()
    if not lex.at_end():
        lex.fail(f"unexpected '{lex.peek()}'", lex.line())
    return name, space, values


def _rank_suffix(path):
    """``(prefix, rank)`` of a ``<prefix>.NNNNNN`` rank-file name, else None."""
    path = str(path)
    dot = path.rfind(".")
    if dot < 0 or len(path) - dot != 7 or not path[dot + 1 :].isdigit():
        return None
    return path[: dot + 1], int(path[dot + 1 :])


def _rank_path(prefix, rank):
    return f"{prefix}{rank:06d}"


def _rank_siblings(path):
    found = _rank_suffix(path)
    return found is not None and all(
        os.path.exists(_rank_path(found[0], r)) for r in (0, 1)
    )


def read(filename, grid_functions=None, piece=None):
    f = _parse(filename)
    if f.get("nurbs") is not None:
        if piece not in (None, 0):
            raise ReadError(
                f"MFEM mesh: {filename} is not parallel; its only piece is 0"
            )
        return _read_nurbs(f["nurbs"], grid_functions, filename)
    if f["parallel"] or _rank_siblings(filename):
        return _read_parallel(filename, f, grid_functions, piece)
    if piece not in (None, 0):
        raise ReadError(f"MFEM mesh: {filename} is not parallel; its only piece is 0")
    dim = f["dim"]
    nv = f["nv"]
    elements = f["elements"]
    sdim = f["sdim"]

    coords = "vertices"
    mesh_order = 1
    node_dofs = 0
    nspace = f["nodes_space"]
    nodes = f["nodes"]
    serendipity_ok = dim == 2 and all(el[1] == 3 for el in elements)
    if nspace is not None and nspace.points == "serendipity" and not serendipity_ok:
        nspace.kind = "h1-other"  # MFEM's serendipity elements are quadrilaterals
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
                "read at the vertices only: only nodal H1 spaces (Gauss-Lobatto or "
                "equispaced) are read as curved cells"
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
        if h1 and space.points == "serendipity" and not serendipity_ok:
            warn(
                f"MFEM grid function '{path}': serendipity fields are read on "
                "quadrilateral meshes only; skipped"
            )
            continue
        if h1 and space.order >= 2 and (coords == "dg" or has_pyramid):
            warn(
                f"MFEM grid function '{path}': an order-{space.order} field on this "
                "mesh is not supported; skipped"
            )
            continue
        gfs.append((name, space, values))
    order = 1 if coords == "dg" else mesh_order
    for _, space, _ in gfs:
        if space.kind == "h1":
            order = max(order, space.order)
    if order >= 2 and has_pyramid:
        warn(
            f"MFEM mesh: '{filename}' has pyramids, which meshio++ holds at order 1 "
            "only; reading the vertices"
        )
        order = 1
        if coords == "h1" and mesh_order >= 2:
            coords = "corners"
    # Bernstein and serendipity coefficients are no nodal values: every order
    # of them goes through the interpolating reader
    modal = [s for s in [nspace if coords == "h1" else None] + [g[1] for g in gfs] if s]
    modal = any(
        s.kind == "h1" and s.points in ("bernstein", "serendipity") for s in modal
    )
    if order >= 3 or (order == 2 and modal):
        if coords == "vertices":
            vxyz = np.array(f["coords"], dtype=np.float64).reshape(nv, sdim)
        else:
            vxyz = np.array(
                [
                    [_value(nodes, nspace, node_dofs, v, c) for c in range(sdim)]
                    for v in range(nv)
                ],
                dtype=np.float64,
            ).reshape(nv, sdim)
        return _read_high_order(f, gfs, coords == "h1", vxyz, order, filename)
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

    mesh.regions = _regions(f, cell_attr, cell_is_boundary)
    return mesh


def _regions(f, cell_attr, cell_is_boundary):
    """Cell regions from the element and boundary attributes and the sets."""
    dim = f["dim"]
    by_attr = {}
    for g in range(len(cell_attr)):
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
    return regions


# --- order 3 and up: VTK Lagrange cells ---------------------------------------------
#
# The Python twin of mf_read_high_order: every element's dofs are placed at
# their reference positions (an edge's from its lower global vertex, a face's in
# the frame of its first element's vertex order, the interior in the element's
# own order) and matched by position to one canonical node set per (geometry,
# space); one matrix per set interpolates to the VTK Lagrange nodes, which cells
# share by their corner weights.

_SHAPES = {1: "line", 2: "triangle", 3: "quad", 4: "tetra", 5: "hexahedron", 6: "wedge"}
_LAGRANGE_TYPES = {
    0: "vertex",
    1: "VTK_LAGRANGE_CURVE",
    2: "VTK_LAGRANGE_TRIANGLE",
    3: "VTK_LAGRANGE_QUADRILATERAL",
    4: "VTK_LAGRANGE_TETRAHEDRON",
    5: "VTK_LAGRANGE_HEXAHEDRON",
    6: "VTK_LAGRANGE_WEDGE",
}
_REF_VERTICES = {
    1: [(0, 0, 0), (1, 0, 0)],
    2: [(0, 0, 0), (1, 0, 0), (0, 1, 0)],
    3: [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)],
    4: [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)],
    5: [
        (0, 0, 0),
        (1, 0, 0),
        (1, 1, 0),
        (0, 1, 0),
        (0, 0, 1),
        (1, 0, 1),
        (1, 1, 1),
        (0, 1, 1),
    ],
    6: [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (0, 1, 1)],
}
# LagrangeHexFiniteElement(3), the Cubic collection's hexahedron: its eight
# interior nodes go round each layer like the corners do.
_CUBIC_HEX_INTERIOR = [
    (1, 1, 1),
    (2, 1, 1),
    (2, 2, 1),
    (1, 2, 1),
    (1, 1, 2),
    (2, 1, 2),
    (2, 2, 2),
    (1, 2, 2),
]


def _mix(terms):
    out = [0.0, 0.0, 0.0]
    for w, p in terms:
        for d in range(3):
            out[d] += w * p[d]
    return tuple(out)


def _entities(elements, dim):
    """Edges and (3-D) faces by first appearance, each face's vertices in its
    first element's order."""
    edges, faces, face_vertices = {}, {}, []
    for _, geom, verts, _ in elements:
        g = _GEOMS[geom]
        if dim >= 2:
            for a, b in g[4]:
                edges.setdefault(tuple(sorted((verts[a], verts[b]))), len(edges))
        if dim == 3:
            for fv in g[5]:
                v = [verts[k] for k in fv]
                key = tuple(sorted(v))
                if key not in faces:
                    faces[key] = len(faces)
                    face_vertices.append(v)
    return edges, faces, face_vertices


def _interior_count(geom, q, points="gll"):
    if points == "serendipity" and geom == 3:
        return (q - 2) * (q - 3) // 2  # MFEM's bubbles from order 4
    return {
        1: q - 1,
        2: (q - 1) * (q - 2) // 2,
        3: (q - 1) ** 2,
        4: (q - 1) * (q - 2) * (q - 3) // 6 if q >= 3 else 0,
        5: (q - 1) ** 3,
        6: (q - 1) * (q - 2) // 2 * (q - 1),
    }.get(geom, 0)


class _Dofs:
    """The dof layout of one nodal space over a mesh."""

    def __init__(self, f, ent, q, points):
        edges, _, face_vertices = ent
        self.order, self.points = q, points
        self.cp = (
            _lagrange.gll_points(q)
            if points in ("gll", "serendipity")
            else _lagrange.uniform_points(q)
        )
        nxt = f["nv"]
        self.edge_base = nxt
        nxt += len(edges) * (q - 1)
        self.face_offset = []
        for fv in face_vertices:
            self.face_offset.append(nxt)
            nxt += (q - 1) * (q - 2) // 2 if len(fv) == 3 else (q - 1) ** 2
        self.interior_offset = []
        for _, geom, _, _ in f["elements"]:
            self.interior_offset.append(nxt)
            if _GEOMS[geom][2] == f["dim"]:
                nxt += _interior_count(geom, q, points)
        self.size = nxt


def _element_dofs(geom, verts, dim, ent, d, element):
    """``(dof, reference position)`` of an element's dofs; ``ent`` None: the
    reference element itself, numbered canonically."""
    g = _GEOMS[geom]
    q, cp = d.order, d.cp
    rv = _REF_VERTICES[geom]
    out = [(verts[j] if ent else j, rv[j]) for j in range(g[3])]
    if q < 2:
        return out
    local = g[3]
    if dim >= 2 and g[2] >= 2:
        for a, b in g[4]:
            lo, hi = (a, b) if verts[a] < verts[b] else (b, a)
            base = 0
            if ent:
                base = d.edge_base + ent[0][tuple(sorted((verts[a], verts[b])))] * (
                    q - 1
                )
            for k in range(1, q):
                t = cp[k]
                out.append(
                    (
                        base + k - 1 if ent else local,
                        _mix([(1 - t, rv[lo]), (t, rv[hi])]),
                    )
                )
                local += 0 if ent else 1
    if dim == 3 and g[2] == 3:
        for fv in g[5]:
            v = [verts[k] for k in fv]
            if ent:
                idx = ent[1][tuple(sorted(v))]
                base = d.face_offset[idx]
                r = [rv[verts.index(gv)] for gv in ent[2][idx]]
            else:
                base = 0
                r = [rv[k] for k in fv]
            o = 0
            if len(r) == 3:
                for j in range(1, q):
                    for i in range(1, q - j):
                        w = cp[i] + cp[j] + cp[q - i - j]
                        x, y = cp[i] / w, cp[j] / w
                        pos = _mix([(1 - x - y, r[0]), (x, r[1]), (y, r[2])])
                        out.append((base + o if ent else local, pos))
                        o += 1
                        local += 0 if ent else 1
            else:
                for j in range(1, q):
                    for i in range(1, q):
                        x, y = cp[i], cp[j]
                        pos = _mix(
                            [
                                ((1 - x) * (1 - y), r[0]),
                                (x * (1 - y), r[1]),
                                (x * y, r[2]),
                                ((1 - x) * y, r[3]),
                            ]
                        )
                        out.append((base + o if ent else local, pos))
                        o += 1
                        local += 0 if ent else 1
    if g[2] != dim:
        return out
    base = d.interior_offset[element] if ent else local
    interior = []
    if geom == 1:
        interior = [(cp[i], 0.0, 0.0) for i in range(1, q)]
    elif geom in (2, 6):
        tri = []
        for j in range(1, q):
            for i in range(1, q - j):
                w = cp[i] + cp[j] + cp[q - i - j]
                tri.append((cp[i] / w, cp[j] / w))
        if geom == 2:
            interior = [(x, y, 0.0) for x, y in tri]
        else:
            interior = [(x, y, cp[k]) for k in range(1, q) for x, y in tri]
    elif geom == 3 and d.points == "serendipity":
        # non-nodal bubbles: positions outside the element only name them
        interior = [
            (-1.0 - n, -1.0, 0.0) for n in range(_interior_count(3, q, d.points))
        ]
    elif geom == 3:
        interior = [(cp[i], cp[j], 0.0) for j in range(1, q) for i in range(1, q)]
    elif geom == 4:
        for k in range(1, q):
            for j in range(1, q - k):
                for i in range(1, q - j - k):
                    w = cp[i] + cp[j] + cp[k] + cp[q - i - j - k]
                    interior.append((cp[i] / w, cp[j] / w, cp[k] / w))
    elif geom == 5:
        if d.points == "cubic":
            interior = [(cp[a], cp[b], cp[c]) for a, b, c in _CUBIC_HEX_INTERIOR]
        else:
            interior = [
                (cp[i], cp[j], cp[k])
                for k in range(1, q)
                for j in range(1, q)
                for i in range(1, q)
            ]
    for n, pos in enumerate(interior):
        out.append((base + n, pos))
    return out


def _quantise(pos):
    return tuple(int(round(x * 1e7)) for x in pos)


def _find_node(index, pos):
    key = _quantise(pos)
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dz in (-1, 0, 1):
                k = index.get((key[0] + dx, key[1] + dy, key[2] + dz))
                if k is not None:
                    return k
    return None


def _cell_order(geom, num_nodes):
    """The Lagrange order of a written cell by its node count: 1 for its
    corners, 2 for the complete quadratic types, p for a VTK Lagrange cell;
    -1 for a serendipity cell."""
    if not 1 <= geom <= 6:
        return 1
    shape = _SHAPES[geom]
    q = 1
    while True:
        n = _lagrange.num_nodes(shape, q)
        if n == num_nodes:
            return q
        if n > num_nodes:
            return -1
        q += 1


def _bernstein_matrix(shape, q, nodes, targets):
    """MFEM's positive (Bernstein) basis at ``targets``: the coefficient at
    lattice point ``a / q`` weighs the Bernstein polynomial of multi-index
    ``a``, a tensor product on quadrilaterals and hexahedra, barycentric on
    simplices, a triangle times a segment on prisms."""
    fact = [math.factorial(k) for k in range(q + 1)]

    def b1(a, x):
        return fact[q] / (fact[a] * fact[q - a]) * x**a * (1 - x) ** (q - a)

    def bs(alpha, lam):
        v = float(fact[q])
        for a, x in zip(alpha, lam):
            v *= x**a / fact[a]
        return v

    out = np.zeros((len(targets), len(nodes)))
    for c, node in enumerate(nodes):
        a = [int(round(x * q)) for x in node]
        for t, x in enumerate(targets):
            if shape == "line":
                v = b1(a[0], x[0])
            elif shape in ("quad", "hexahedron"):
                v = 1.0
                for d in range(2 if shape == "quad" else 3):
                    v *= b1(a[d], x[d])
            elif shape == "triangle":
                v = bs((q - a[0] - a[1], a[0], a[1]), (1 - x[0] - x[1], x[0], x[1]))
            elif shape == "tetra":
                v = bs(
                    (q - a[0] - a[1] - a[2], a[0], a[1], a[2]),
                    (1 - x[0] - x[1] - x[2], x[0], x[1], x[2]),
                )
            else:  # wedge
                v = bs((q - a[0] - a[1], a[0], a[1]), (1 - x[0] - x[1], x[0], x[1]))
                v *= b1(a[2], x[2])
            out[t, c] = v
    return out


def _serendipity_shape(p, cp, x, y):
    """``H1Ser_QuadrilateralElement::CalcShape`` (MFEM fe_ser.cpp): nodal
    Gauss-Lobatto edge functions times the linear function vanishing on the
    opposite edge, bilinear vertex functions corrected by them, and Legendre
    bubbles from order 4; in MFEM's local order."""

    def lag(t):
        out = np.ones(p + 1)
        for i in range(p + 1):
            for j in range(p + 1):
                if j != i:
                    out[i] *= (t - cp[j]) / (cp[i] - cp[j])
        return out

    nx, ny = lag(x), lag(y)
    n = 4 + 4 * (p - 1) + (p - 2) * (p - 3) // 2
    shape = np.zeros(n)
    for i in range(p - 1):
        shape[4 + 0 * (p - 1) + i] = nx[i + 1] * (1 - y)
        shape[4 + 1 * (p - 1) + i] = ny[i + 1] * x
        shape[4 + 3 * (p - 1) - i - 1] = nx[i + 1] * y
        shape[4 + 4 * (p - 1) - i - 1] = ny[i + 1] * (1 - x)
    bil = [(1 - x) * (1 - y), x * (1 - y), x * y, (1 - x) * y]
    fix = [0.0] * 4
    for i in range(p - 1):
        w = 1 - cp[i + 1]
        fix[0] += w * (shape[4 + i] + shape[4 + 4 * (p - 1) - i - 1])
        fix[1] += w * (shape[4 + 1 * (p - 1) + i] + shape[4 + (p - 2) - i])
        fix[2] += w * (shape[4 + 2 * (p - 1) + i] + shape[1 + 2 * p - i])
        fix[3] += w * (shape[4 + 3 * (p - 1) + i] + shape[3 * p - i])
    for v in range(4):
        shape[v] = bil[v] - fix[v]
    if p > 3:

        def leg(t):
            u = [1.0, 2 * t - 1]
            for k in range(1, p - 2):
                u.append(((2 * k + 1) * (2 * t - 1) * u[k] - k * u[k - 1]) / (k + 1))
            return u

        lx, ly = leg(x), leg(y)
        m = 0
        for j in range(4, p + 1):
            for k in range(j - 3):
                shape[4 + 4 * (p - 1) + m] = (
                    lx[k] * ly[j - 4 - k] * x * (1 - x) * y * (1 - y)
                )
                m += 1
    return shape


def _serendipity_matrix(p, cp, nodes, targets):
    """MFEM's serendipity basis at ``targets``, one column per canonical node:
    the vertex, edge (by position) and bubble (by index) functions."""
    local = []
    for x, y, _ in nodes:
        if x < 0:  # a bubble, named by its index
            local.append(4 + 4 * (p - 1) + int(round(-1.0 - x)))
            continue
        corner = {(0, 0): 0, (1, 0): 1, (1, 1): 2, (0, 1): 3}.get(
            (round(x, 12), round(y, 12))
        )
        if corner is not None:
            local.append(corner)
            continue
        on = [
            k for k in range(1, p) if abs(cp[k] - (x if y in (0.0, 1.0) else y)) < 1e-12
        ]
        i = on[0] - 1
        if abs(y) < 1e-12:
            local.append(4 + i)
        elif abs(x - 1) < 1e-12:
            local.append(4 + (p - 1) + i)
        elif abs(y - 1) < 1e-12:
            local.append(4 + 3 * (p - 1) - i - 1)
        else:
            local.append(4 + 4 * (p - 1) - i - 1)
    out = np.zeros((len(targets), len(nodes)))
    for t, (x, y, _) in enumerate(targets):
        out[t] = _serendipity_shape(p, cp, x, y)[local]
    return out


class _Interp:
    """One (geometry, space) pair's canonical nodes and its matrix to the VTK
    Lagrange nodes of the output order."""

    def __init__(self, geom, dim, d, order):
        nodes = _element_dofs(
            geom, list(range(len(_REF_VERTICES[geom]))), dim, None, d, 0
        )
        pos = [None] * len(nodes)
        self.index = {}
        for k, p in nodes:
            pos[k] = p
            self.index[_quantise(p)] = k
        self.nodes = len(nodes)
        shape = _SHAPES[geom]
        targets = [
            (i / order, j / order, k / order)
            for i, j, k in _lagrange.vtk_lattice(shape, order)
        ]
        if d.points == "bernstein":
            self.matrix = _bernstein_matrix(shape, d.order, pos, targets)
        elif d.points == "serendipity":
            self.matrix = _serendipity_matrix(d.order, d.cp, pos, targets)
        else:
            self.matrix = _lagrange.interpolation_matrix(shape, d.order, pos, targets)

    def find(self, pos):
        return _find_node(self.index, pos)


def _cell_type(geom, order):
    """``mf_cell_type``: VTK Lagrange at order 3 and up, else the linear or
    complete quadratic type."""
    if geom == 0 or order <= 1:
        return _GEOMS[geom][0]
    if order == 2:
        return _GEOMS[geom][1]
    return _LAGRANGE_TYPES[geom]


def _field_table(values, space, sdim=None):
    vdim = space.vdim
    ndofs = len(values) // vdim
    vals = np.asarray(values, dtype=np.float64)
    table = (
        vals.reshape(vdim, ndofs).T
        if space.ordering == 0
        else vals.reshape(ndofs, vdim)
    )
    return table if sdim is None else table[:, :sdim]


def _read_parts(parts, nvertices, order, labels, filename):
    """``mf_read_parts``: the parts (dicts with the parsed file ``f``, its grid
    functions ``gfs``, ``nodes_h1``, the vertex coordinates ``vxyz`` and the
    output number ``global`` of each local vertex) as one mesh of order-``order``
    cells, nodes keyed by the output vertices they combine."""
    first = parts[0]["f"]
    dim, sdim = first["dim"], first["sdim"]
    point_gfs = [
        g
        for g, (_, space, _) in enumerate(parts[0]["gfs"])
        if space.kind in ("h1", "nurbs")
    ]
    ents, fields = [], []
    for part in parts:
        f = part["f"]
        if f["dim"] != dim or f["sdim"] != sdim:
            raise ReadError(
                f"MFEM mesh: the parts of {filename} disagree on the dimension"
            )
        for _, geom, _, line in f["elements"]:
            if geom == 0 or _GEOMS[geom][2] != dim or (geom == 7 and order > 1):
                raise ReadError(
                    f"MFEM mesh: element of geometry {geom} in a {dim}-D mesh of order "
                    f"{order} (line {line})"
                )
        ent = _entities(f["elements"], dim)
        ents.append(ent)
        mine = []
        if "eval" in part:  # values come from the part itself (NURBS)
            fields.append([(None, np.empty((0, c))) for c in part["ncomp"]])
            continue
        if part["nodes_h1"]:
            s = f["nodes_space"]
            d = _Dofs(f, ent, s.order, s.points)
            ndofs = len(f["nodes"]) // s.vdim
            if ndofs != d.size:
                raise ReadError(
                    f"MFEM mesh: the order-{s.order} nodes of {filename} have {ndofs} "
                    f"dofs; the mesh numbers {d.size}"
                )
            mine.append((d, _field_table(f["nodes"], s, sdim)))
        else:
            mine.append((_Dofs(f, ent, 1, "gll"), part["vxyz"]))
        for g in point_gfs:
            name, space, values = part["gfs"][g]
            if len(values) % space.vdim:
                raise ReadError(
                    f"MFEM grid function '{name}': {len(values)} values, not a multiple "
                    f"of VDim {space.vdim}"
                )
            d = _Dofs(f, ent, space.order, space.points)
            ndofs = len(values) // space.vdim
            if ndofs != d.size:
                raise ReadError(
                    f"MFEM grid function '{name}' has {ndofs} dofs; the mesh has {d.size} "
                    f"at order {space.order}"
                )
            mine.append((d, _field_table(values, space)))
        fields.append(mine)

    p3 = order**3
    node_of = {((v, p3),): v for v in range(nvertices)}

    def cell_nodes(part, geom, verts):
        glob = part["global"]
        if geom in (0, 7):  # a point; a (linear) pyramid
            return [glob[v] for v in verts]
        shape = _SHAPES[geom]
        ids = []
        for ijk in _lagrange.vtk_lattice(shape, order):
            key = tuple(
                sorted(
                    (glob[verts[c]], w)
                    for c, w in _lagrange.lattice_weights(shape, order, ijk)
                )
            )
            ids.append(node_of.setdefault(key, len(node_of)))
        return ids

    element_nodes = [
        [cell_nodes(part, g, v) for _, g, v, _ in part["f"]["elements"]]
        for part in parts
    ]
    boundary_nodes = [
        [cell_nodes(part, g, v) for _, g, v, _ in part["f"]["boundary"]]
        for part in parts
    ]
    npoints = len(node_of)

    values = [np.full((npoints, t.shape[1]), np.nan) for _, t in fields[0]]
    known = np.zeros(npoints, dtype=bool)
    for q, part in enumerate(parts):
        f = part["f"]
        caches = [{} for _ in fields[q]]
        for e, (_, geom, verts, _) in enumerate(f["elements"]):
            ids = element_nodes[q][e]
            fresh = [t for t, i in enumerate(ids) if not known[i]]
            if not fresh:
                continue
            if "eval" in part:
                for k, vals in enumerate(part["eval"](e)):
                    values[k][[ids[t] for t in fresh]] = vals[fresh]
                known[ids] = True
                continue
            if geom == 7:  # a linear pyramid: its vertices
                for k, (_, table) in enumerate(fields[q]):
                    values[k][ids] = table[verts]
                known[ids] = True
                continue
            for k, (d, table) in enumerate(fields[q]):
                interp = caches[k].get(geom)
                if interp is None:
                    interp = caches[k][geom] = _Interp(geom, dim, d, order)
                perm = [None] * interp.nodes
                for dof, pos in _element_dofs(geom, verts, dim, ents[q], d, e):
                    c = interp.find(pos)
                    if c is None:
                        raise ReadError(
                            f"MFEM mesh: a degree of freedom of element {e} matches no "
                            "node of its element"
                        )
                    perm[c] = dof
                values[k][[ids[t] for t in fresh]] = interp.matrix[fresh] @ table[perm]
            known[ids] = True
    for part in parts:
        if "eval" in part:
            continue
        for v, g in enumerate(part["global"]):
            if not known[g]:
                values[0][g] = part["vxyz"][v]
                known[g] = True
    orphans = 0
    for q, part in enumerate(parts):
        for b, (_, geom, verts, _) in enumerate(part["f"]["boundary"]):
            if geom in (0, 7):
                continue
            shape = _SHAPES[geom]
            for t, ijk in enumerate(_lagrange.vtk_lattice(shape, order)):
                i = boundary_nodes[q][b][t]
                if known[i]:
                    continue
                orphans += 1
                values[0][i] = sum(
                    w / p3 * part["vxyz"][verts[c]]
                    for c, w in _lagrange.lattice_weights(shape, order, ijk)
                )
                known[i] = True
    if orphans:
        warn(
            f"MFEM mesh: {orphans} boundary node(s) lie on no element face; placed "
            "from their corners"
        )

    blocks = []  # (type, [(part, source, is_boundary)])
    for is_boundary in (False, True):
        index = {}
        for q, part in enumerate(parts):
            items = part["f"]["boundary" if is_boundary else "elements"]
            for k, el in enumerate(items):
                t = _cell_type(el[1], order)
                if t not in index:
                    index[t] = len(blocks)
                    blocks.append((t, []))
                blocks[index[t]][1].append((q, k, is_boundary))
    cells, attr_blocks, part_blocks = [], [], []
    cell_attr, cell_is_boundary = [], []
    element_cell = [[0] * len(part["f"]["elements"]) for part in parts]
    for cell_type, members in blocks:
        rows, attrs, ranks = [], [], []
        for q, source, is_boundary in members:
            f = parts[q]["f"]
            attr = f["boundary" if is_boundary else "elements"][source][0]
            rows.append((boundary_nodes if is_boundary else element_nodes)[q][source])
            attrs.append(attr)
            ranks.append(f["rank"])
            cell_attr.append(attr)
            cell_is_boundary.append(is_boundary)
            if not is_boundary:
                element_cell[q][source] = len(cell_attr) - 1
        cells.append((cell_type, np.array(rows, dtype=np.int64)))
        attr_blocks.append(np.array(attrs, dtype=np.int64))
        part_blocks.append(np.array(ranks, dtype=np.int64))
    mesh = Mesh(values[0], cells)
    if attr_blocks:
        mesh.cell_data["mfem:attribute"] = attr_blocks
    if labels and part_blocks:
        mesh.cell_data["partition:part"] = part_blocks
    for g, data in zip(point_gfs, values[1:]):
        name = parts[0]["gfs"][g][0]
        mesh.point_data[name] = data[:, 0].copy() if data.shape[1] == 1 else data
    ncells = len(cell_attr)
    for g, (name, space, _) in enumerate(parts[0]["gfs"]):
        if space.kind in ("h1", "nurbs"):
            continue
        vdim = space.vdim
        per_cell = np.full((ncells, vdim), np.nan)
        for q, part in enumerate(parts):
            _, gspace, gvalues = part["gfs"][g]
            ne = len(part["f"]["elements"])
            ndofs = len(gvalues) // vdim
            if ndofs != ne:
                raise ReadError(
                    f"MFEM grid function '{name}' has {ndofs} dofs for {ne} elements"
                )
            for e in range(ne):
                for c in range(vdim):
                    per_cell[element_cell[q][e], c] = _value(
                        gvalues, gspace, ndofs, e, c
                    )
        out, start = [], 0
        for _, members in blocks:
            a = per_cell[start : start + len(members)]
            start += len(members)
            out.append(a[:, 0].copy() if vdim == 1 else a.copy())
        mesh.cell_data[name] = out
    mesh.regions = _regions(first, cell_attr, cell_is_boundary)
    return mesh


def _read_nurbs(n, grid_functions, filename):
    """A NURBS mesh as its knot-span elements: VTK Lagrange cells (or linear
    and quadratic ones) of the highest knot-vector order, their nodes the
    rational patch geometry at the cell's lattice; NURBS grid functions on
    the mesh's own space the same way, element-wise ones as cell data."""
    elements = n.elements()
    boundary = n.boundary()
    weights, xyz = n.control_points()
    order = max([1] + [kv.order for row in n.compr for kv in row])
    if isinstance(grid_functions, dict):
        gf_items = list(grid_functions.items())
    else:
        gf_items = [(pathlib.Path(p).stem, p) for p in (grid_functions or [])]
    gfs, tables = [], []
    for name, path in gf_items:
        name, space, values = _parse_gf(str(name), str(path))
        if space.kind == "nurbs" and len(values) == n.num_dofs * space.vdim:
            gfs.append((name, space, values))
            tables.append(_field_table(values, space))
        elif space.kind in ("l2", "l2t1") and space.order == 0:
            gfs.append((name, space, values))
        else:
            warn(
                f"MFEM grid function '{path}': the '{space.collection}' space is not "
                f"the NURBS mesh's own nor element-wise; skipped"
            )
    f = {
        "dim": n.dim,
        "sdim": xyz.shape[1],
        "elements": [el[:4] for el in elements],
        "boundary": boundary,
        "nv": n.num_vertices,
        "rank": 0,
        "sets": [],
        "bdr_sets": [],
    }
    shape = _SHAPES[_nurbs.GEOM_OF_DIM[n.dim]]
    refs = [
        tuple(c / order for c in ijk[: n.dim])
        for ijk in _lagrange.vtk_lattice(shape, order)
    ]

    def evaluate(e):
        return [n.evaluate(weights, xyz, elements[e], refs)] + [
            n.evaluate(weights, t, elements[e], refs) for t in tables
        ]

    part = {
        "f": f,
        "gfs": gfs,
        "global": list(range(n.num_vertices)),
        "eval": evaluate,
        "ncomp": [xyz.shape[1]] + [t.shape[1] for t in tables],
    }
    return _read_parts([part], n.num_vertices, order, False, filename)


def _read_high_order(f, gfs, nodes_h1, vxyz, order, filename):
    part = {
        "f": f,
        "gfs": list(gfs),
        "nodes_h1": nodes_h1,
        "vxyz": vxyz,
        "global": list(range(f["nv"])),
    }
    return _read_parts([part], f["nv"], order, False, filename)


def _read_parallel(filename, first, grid_functions, piece):
    """``mf_read_parallel``: every rank file beside ``filename`` (or the one
    ``piece``), vertices merged through the communication groups (ParPrint) or,
    without them (ParMesh::Save), boundary vertices by position; a boundary
    face two ranks list is their interface and is dropped."""
    found = _rank_suffix(filename)
    paths = []
    if found is not None:
        r = 0
        while os.path.exists(_rank_path(found[0], r)):
            paths.append(_rank_path(found[0], r))
            r += 1
    if not paths:
        warn(
            f"MFEM mesh: '{filename}' is one rank of a parallel mesh, not named "
            "<prefix>.NNNNNN beside its siblings; reading that rank alone"
        )
        paths = [str(filename)]
    if piece is not None:
        if not 0 <= piece < len(paths):
            raise ReadError(
                f"MFEM mesh: piece {piece} is out of range: {filename} has {len(paths)} ranks"
            )
        selected = [piece]
    else:
        selected = list(range(len(paths)))
    files = []
    for r in selected:
        f = first if paths[r] == str(filename) else _parse(paths[r])
        if f["parallel"] and len(paths) > 1 and f["rank"] != r:
            raise ReadError(f"MFEM mesh: {paths[r]} holds rank {f['rank']}")
        if not f["parallel"]:
            f["rank"] = r
        files.append(f)
    groups = all(f["parallel"] for f in files)

    parts, pyramid, corners, order = [], False, False, 1
    for q, f in enumerate(files):
        sdim, nv = f["sdim"], f["nv"]
        part = {"f": f, "gfs": [], "nodes_h1": False}
        s = f["nodes_space"]
        if s is not None:
            if s.kind not in ("h1", "h1-other"):
                raise ReadError(
                    f"MFEM mesh: nodes in '{s.collection}' are not supported in a parallel "
                    "mesh"
                )
            table = _field_table(f["nodes"], s, sdim)
            if len(table) < nv:
                raise ReadError(
                    f"MFEM mesh: the nodes of {paths[selected[q]]} have {len(table)} dofs "
                    f"for {nv} vertices"
                )
            part["vxyz"] = table[:nv].copy()
            part["nodes_h1"] = s.kind == "h1"
            corners = corners or s.kind == "h1-other"
            if part["nodes_h1"]:
                order = max(order, s.order)
        else:
            part["vxyz"] = np.array(f["coords"], dtype=np.float64).reshape(nv, sdim)
        pyramid = pyramid or any(el[1] == 7 for el in f["elements"])
        parts.append(part)

    # output vertices
    allxyz = np.concatenate([p["vxyz"] for p in parts])
    extent = float(np.ptp(allxyz, axis=0).max()) if len(allxyz) else 0.0
    tol = max(extent, 1.0) * 1e-9
    buckets, by_group, global_xyz = {}, {}, []
    nglobal = 0
    for q, f in enumerate(files):
        part = parts[q]
        glob = [-1] * f["nv"]
        sdim = f["sdim"]
        candidate = [f["dim"] == 1 and not groups] * f["nv"]
        if groups:
            for g in range(1, len(f["group_vertices"])):
                for k, v in enumerate(f["group_vertices"][g]):
                    if not 0 <= v < f["nv"]:
                        raise ReadError(
                            f"MFEM mesh: shared vertex {v} out of range in "
                            f"{paths[selected[q]]}"
                        )
                    key = (f["groups"][g], k)
                    if key not in by_group:
                        by_group[key] = nglobal
                        nglobal += 1
                        global_xyz.append(None)
                    glob[v] = by_group[key]
        elif f["nc"]:
            # a non-conforming rank (ParPrint): the vertices its ghosts share
            for v in f["interface"]:
                candidate[v] = True
        else:
            for _, _, verts, _ in f["boundary"]:
                for v in verts:
                    candidate[v] = True
        for v in range(f["nv"]):
            if not candidate[v]:
                continue
            x = part["vxyz"][v]
            cell = tuple(int(np.floor(c / tol)) for c in x)
            match = None
            for m in range(3**sdim):
                code, probe = m, []
                for c in range(sdim):
                    probe.append(cell[c] + code % 3 - 1)
                    code //= 3
                for g in buckets.get(tuple(probe), []):
                    if np.all(np.abs(global_xyz[g] - x) <= tol):
                        match = g
                        break
                if match is not None:
                    break
            if match is None:
                match = nglobal
                nglobal += 1
                global_xyz.append(np.array(x))
                buckets.setdefault(cell, []).append(match)
            glob[v] = match
        for v in range(f["nv"]):
            if glob[v] < 0:
                glob[v] = nglobal
                nglobal += 1
                global_xyz.append(None)
        part["global"] = glob
    if not groups and len(files) > 1:
        owners = {}
        for q, f in enumerate(files):
            for b in f["boundary"]:
                key = tuple(sorted(parts[q]["global"][v] for v in b[2]))
                owners.setdefault(key, set()).add(q)
        for q, f in enumerate(files):
            f["boundary"] = [
                b
                for b in f["boundary"]
                if len(owners[tuple(sorted(parts[q]["global"][v] for v in b[2]))]) == 1
            ]
    if corners:
        warn(
            f"MFEM mesh: nodes of {filename} that are not nodal H1 are read at the "
            "vertices only"
        )

    if isinstance(grid_functions, dict):
        items = list(grid_functions.items())
    else:
        items = [(pathlib.Path(p).stem, p) for p in (grid_functions or [])]
    for name, path in items:
        found = _rank_suffix(path)
        if found is None:
            raise ReadError(
                f"MFEM grid function '{path}': a parallel mesh's grid function is named "
                "by one of its rank files, <name>.NNNNNN"
            )
        per_rank = [
            _parse_gf(str(name), _rank_path(found[0], f["rank"])) for f in files
        ]
        space = per_rank[0][1]
        h1 = space.kind == "h1"
        l2p0 = space.kind in ("l2", "l2t1") and space.order == 0
        if not h1 and not l2p0:
            warn(
                f"MFEM grid function '{path}': the '{space.collection}' space is not "
                "supported; skipped"
            )
            continue
        if h1 and space.order >= 2 and pyramid:
            warn(
                f"MFEM grid function '{path}': an order-{space.order} field on this mesh "
                "is not supported; skipped"
            )
            continue
        if h1:
            order = max(order, space.order)
        for q in range(len(files)):
            parts[q]["gfs"].append(per_rank[q])
    if pyramid and order >= 2:
        warn(
            f"MFEM mesh: '{filename}' has pyramids, which meshio++ holds at order 1 "
            "only; reading the vertices"
        )
        order = 1
        for part in parts:
            part["nodes_h1"] = False
    return _read_parts(parts, nglobal, order, True, filename)


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
            block.type in (_GEOMS[g][0], _GEOMS[g][1])
            or block.type in _WRITABLE_EXTRA
            or block.type in _LAGRANGE_GEOMS
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
    # VTK Lagrange cells: the whole mesh is written with order-p H1 nodes, p the
    # highest order among its cells; linear and complete quadratic cells are
    # evaluated at those nodes too, serendipity cells from their corners.
    high = 0
    for c in elements + boundary:
        if c[4] in _LAGRANGE_GEOMS:
            high = max(high, _cell_order(c[0], len(c[3])))
    serendipity = 0
    if high > 0:
        for c in elements + boundary:
            q = _cell_order(c[0], len(c[3]))
            serendipity += q < 0
            high = max(high, q)
        quadratic = False
    if (quadratic or high > 1) and has_pyramid:
        warn(
            f"MFEM mesh writer: MFEM order-{high if high > 1 else 2} meshes hold no "
            "pyramids here; writing corners only"
        )
        _provenance.note(
            "high-order-dropped",
            "a mesh with pyramids is written with linear MFEM cells",
        )
        quadratic = False
        high = 0
    if high == 1:
        high = 0  # order-1 Lagrange cells are their corners
    if high > 0 and serendipity:
        warn(
            f"MFEM mesh writer: {serendipity} serendipity cell(s) in an order-{high} mesh "
            "are placed from their corners"
        )
        _provenance.note(
            "high-order-dropped",
            "serendipity cells in a Lagrange mesh are written from corners",
        )
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
    if (
        not quadratic
        and high == 0
        and all(len(c[3]) <= _GEOMS[c[0]][3] for c in elements + boundary)
    ):
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

    # order-p H1 dofs (VTK Lagrange meshes): each element evaluated at its dofs'
    # positions, a dof shared by several taking the first's value
    hdofs = None
    high_map = []
    if high > 0:
        numbered = [(c[1], c[0], mfem_vertices(c), 0) for c in elements]
        hf = {"dim": dim, "nv": nv, "elements": numbered}
        hent = _entities(numbered, dim)
        hdofs = _Dofs(hf, hent, high, "gll")
        canon, matrices = {}, {}
        taken = [False] * hdofs.size
        for e, c in enumerate(elements):
            geom = c[0]
            if geom not in canon:
                nodes = _element_dofs(
                    geom, list(range(_GEOMS[geom][3])), dim, None, hdofs, 0
                )
                pos = [None] * len(nodes)
                index = {}
                for k, p in nodes:
                    pos[k] = p
                    index[_quantise(p)] = k
                canon[geom] = (pos, index)
            pos, index = canon[geom]
            q = max(1, _cell_order(geom, len(c[3])))
            if (geom, q) not in matrices:
                shape = _SHAPES[geom]
                src = [
                    (i / q, j / q, k / q) for i, j, k in _lagrange.vtk_lattice(shape, q)
                ]
                matrices[(geom, q)] = _lagrange.interpolation_matrix(shape, q, src, pos)
            pairs = []
            for dof, p in _element_dofs(geom, numbered[e][2], dim, hent, hdofs, e):
                if taken[dof]:
                    continue
                k = _find_node(index, p)
                if k is None:
                    raise WriteError(
                        f"MFEM mesh writer: a node of element {e} matches no degree of "
                        "freedom"
                    )
                taken[dof] = True
                pairs.append((dof, k))
            high_map.append((matrices[(geom, q)], pairs, q))

    def high_values(data, cols):
        out_v = np.zeros((hdofs.size, cols))
        for e, c in enumerate(elements):
            matrix, pairs, q = high_map[e]
            if not pairs:
                continue
            used = c[3] if q > 1 else c[3][: _GEOMS[c[0]][3]]
            u = data[used]
            dofs = [d for d, _ in pairs]
            out_v[dofs] = matrix[[k for _, k in pairs]] @ u
        return out_v

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
    if high > 0:
        out.append(
            f"\nnodes\nFiniteElementSpace\nFiniteElementCollection: H1_{dim}D_P{high}\n"
            f"VDim: {pdim}\nOrdering: 1\n\n"
        )
        for row in high_values(pts.reshape(-1, pdim), pdim):
            out.append(" ".join(_fmt(v) for v in row) + "\n")
    elif not quadratic:
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
    order = high if high > 0 else (2 if quadratic else 1)
    ndofs = hdofs.size if high > 0 else (len(numbering) if quadratic else nv)
    for name in sorted(mesh.point_data):
        a = np.asarray(mesh.point_data[name])
        cols = a.shape[1] if a.ndim > 1 else 1
        flat_a = a.astype(np.float64).ravel()
        text = [
            f"FiniteElementSpace\nFiniteElementCollection: H1_{dim}D_P{order}\n"
            f"VDim: {cols}\nOrdering: 1\n\n"
        ]
        hv = high_values(flat_a.reshape(-1, cols), cols) if high > 0 else None
        for d in range(ndofs):
            if high > 0:
                vals = hv[d]
            elif quadratic:
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

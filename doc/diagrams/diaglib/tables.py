"""Loaders for the code's own topology tables, plus reference-element geometry.

The cell-type figures must not drift from the code, so nothing about node
ordering is typed here twice: mid-edge nodes come from ``_convert_cells._ELEVATE``
(target node ``ncorners + k`` is the midpoint of ``edges[k]``), face rows from
``_skin._CELL_FACES``, edge rows from ``_surface._CELL_EDGES``, refinement
masks from ``_refine_templates``. Those modules have package-relative imports,
so the literal tables are read with :mod:`ast` straight from the source text
(``_refine_templates.py`` imports only ``itertools`` and is loaded as a
module). Only the *corner* coordinates of each linear reference element live
here -- the code has no coordinates -- transcribed from the old
``doc/cell_types.tex`` and cross-checked against
``tests/cpp/test_skin.cpp``'s reference elements.
"""

import ast
import importlib.util
import pathlib

from .project import centroid

HERE = pathlib.Path(__file__).resolve().parent
REPO = HERE.parent.parent.parent
PKG = REPO / "src" / "python" / "meshioplusplus"


def literal_table(module, name):
    """``ast.literal_eval`` of the module-level assignment ``name = <literal>``."""
    tree = ast.parse((PKG / module).read_text(encoding="utf-8"))
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(
            isinstance(t, ast.Name) and t.id == name for t in node.targets
        ):
            return ast.literal_eval(node.value)
    raise KeyError(f"{module} has no literal assignment named {name}")


def load_refine_templates():
    """The real ``_refine_templates`` module (it imports only ``itertools``)."""
    spec = importlib.util.spec_from_file_location(
        "meshioplusplus_refine_templates", PKG / "_refine_templates.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# --------------------------------------------------------------------------- #
# reference elements                                                           #
# --------------------------------------------------------------------------- #
#: Corner coordinates of each linear reference element, VTK/meshio order.
REFERENCE_CORNERS = {
    "vertex": [(0, 0, 0)],
    "line": [(0, 0, 0), (1, 0, 0)],
    "triangle": [(0, 0, 0), (1, 0, 0), (0, 1, 0)],
    "quad": [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)],
    "tetra": [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)],
    "hexahedron": [
        (0, 0, 0),
        (1, 0, 0),
        (1, 1, 0),
        (0, 1, 0),
        (0, 0, 1),
        (1, 0, 1),
        (1, 1, 1),
        (0, 1, 1),
    ],
    "wedge": [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (0, 1, 1)],
    "pyramid": [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0.5, 0.5, 1)],
}

#: Linear base of every drawn type.
LINEAR_BASE = {
    "vertex": "vertex",
    "line": "line",
    "line3": "line",
    "triangle": "triangle",
    "triangle6": "triangle",
    "quad": "quad",
    "quad8": "quad",
    "quad9": "quad",
    "tetra": "tetra",
    "tetra10": "tetra",
    "hexahedron": "hexahedron",
    "hexahedron20": "hexahedron",
    "hexahedron24": "hexahedron",
    "hexahedron27": "hexahedron",
    "wedge": "wedge",
    "wedge15": "wedge",
    "wedge18": "wedge",
    "pyramid": "pyramid",
    "pyramid13": "pyramid",
    "pyramid14": "pyramid",
}

#: The quadratic (serendipity) type each linear type elevates to.
QUADRATIC_OF = {
    "line": "line3",
    "triangle": "triangle6",
    "quad": "quad8",
    "tetra": "tetra10",
    "hexahedron": "hexahedron20",
    "wedge": "wedge15",
    "pyramid": "pyramid13",
}


def _mid(a, b):
    return tuple((x + y) / 2.0 for x, y in zip(a, b))


class Tables:
    """All code tables, loaded once."""

    def __init__(self):
        self.cell_faces = literal_table("_skin.py", "_CELL_FACES")
        self.cell_edges = literal_table("_surface.py", "_CELL_EDGES")
        self.elevate = literal_table("_convert_cells.py", "_ELEVATE")
        self.num_nodes = literal_table("_common.py", "num_nodes_per_cell")
        self.refine = load_refine_templates()

    # -- node coordinates ----------------------------------------------------
    def node_coordinates(self, cell_type):
        """Every node of ``cell_type`` in meshio order, from the tables."""
        base = LINEAR_BASE[cell_type]
        corners = [tuple(float(v) for v in c) for c in REFERENCE_CORNERS[base]]
        if cell_type == base:
            return corners
        nodes = list(corners)
        if base in self.elevate:
            target, ncorners, edges = self.elevate[base]
            assert ncorners == len(corners), (cell_type, ncorners)
            nodes += [_mid(corners[a], corners[b]) for a, b in edges]
            if cell_type == target:
                return nodes
        if cell_type == "quad9":
            nodes.append(centroid(corners))
        elif cell_type in ("hexahedron24", "hexahedron27"):
            rows = self.refine.QUAD_FACES["hexahedron"]
            rows = rows[:4] if cell_type == "hexahedron24" else rows
            nodes += [centroid([corners[i] for i in row]) for row in rows]
            if cell_type == "hexahedron27":
                nodes.append(centroid(corners))
        elif cell_type == "wedge18":
            nodes += [
                centroid([corners[i] for i in row])
                for row in self.refine.QUAD_FACES["wedge"]
            ]
        elif cell_type == "pyramid14":
            nodes.append(centroid(corners[:4]))
        else:
            raise KeyError(cell_type)
        assert len(nodes) == self.num_nodes[cell_type], (cell_type, len(nodes))
        return nodes

    def node_kind(self, cell_type, index):
        """``corner`` / ``edge`` / ``face`` / ``body`` for the node's drawing style."""
        base = LINEAR_BASE[cell_type]
        nc = len(REFERENCE_CORNERS[base])
        if index < nc:
            return "corner"
        n_edges = len(self.elevate[base][2]) if base in self.elevate else 0
        if index < nc + n_edges:
            return "edge"
        if cell_type == "hexahedron27" and index == 26:
            return "body"
        return "face"

    # -- topology ------------------------------------------------------------
    def corner_edges(self, cell_type):
        """The corner-to-corner edges of the linear base, for the wireframe."""
        base = LINEAR_BASE[cell_type]
        if base == "vertex":
            return []
        if base == "line":
            return [(0, 1)]
        if base in self.elevate:
            return list(self.elevate[base][2])
        return list(self.refine.EDGES[base])

    def corner_faces(self, cell_type):
        """Outward-wound corner rings of the linear base (3-D types only)."""
        base = LINEAR_BASE[cell_type]
        if base not in self.cell_faces:
            return []
        return [tuple(row[2][: row[1]]) for row in self.cell_faces[base]]

    def edge_midpoint_index(self, cell_type, a, b):
        """The node sitting on edge ``(a, b)`` of a quadratic type, or ``None``."""
        base = LINEAR_BASE[cell_type]
        if base not in self.elevate:
            return None
        nc = len(REFERENCE_CORNERS[base])
        for k, (u, v) in enumerate(self.elevate[base][2]):
            if {u, v} == {a, b}:
                return nc + k
        return None

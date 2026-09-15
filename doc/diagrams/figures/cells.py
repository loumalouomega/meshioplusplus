"""Node-ordering figures for ``doc/cell_types.md``, one per cell family.

Everything about the ordering comes from the code's own tables (see
``diaglib/tables.py``); the only hand-typed geometry is each linear
element's corner positions.
"""

from diaglib import palette as P
from diaglib.cellart import draw_cell
from diaglib.svg import Canvas
from diaglib.tables import Tables

_TABLES = None


def tables():
    global _TABLES
    if _TABLES is None:
        _TABLES = Tables()
    return _TABLES


# Per-type label nudges (node index -> (dx, dy)) where the default up-right
# placement would collide with an edge or another label.
_NUDGE = {}


def _family(name, types, width, height, scale, three_d, caption):
    c = Canvas(width, height, f"Node ordering of the {name} cell family")
    n = len(types)
    panel = (width - 40) / n
    t = tables()
    for k, cell_type in enumerate(types):
        x0 = 20 + k * panel
        cx = x0 + panel / 2
        c.text(cx, 30, cell_type, size=P.SIZE_TITLE, weight="700", mono=True)
        c.label(cx, 48, f"{t.num_nodes[cell_type]} nodes")
        box = (x0 + 24, 64, panel - 48, height - 64 - 64)
        draw_cell(
            c,
            t,
            cell_type,
            None,
            scale,
            label_offsets=_NUDGE.get(cell_type),
            flat=not three_d,
            box=box,
        )
    c.label(width / 2, height - 34, caption)
    c.legend(
        20,
        height - 12,
        [
            (P.INK, "corner"),
            (P.CORE, "mid-edge"),
            (P.PYTHON, "face centre"),
            (P.CABI, "body centre"),
        ],
    )
    return c.render()


def cell_types_line():
    return _family(
        "line",
        ["vertex", "line", "line3"],
        720,
        180,
        150,
        False,
        "1-D: node 2 of line3 is the midpoint",
    )


def cell_types_triangle():
    return _family(
        "triangle",
        ["triangle", "triangle6"],
        560,
        300,
        140,
        False,
        "mid-edge nodes follow the edges 01, 12, 20",
    )


def cell_types_quad():
    return _family(
        "quad",
        ["quad", "quad8", "quad9"],
        780,
        300,
        140,
        False,
        "mid-edge nodes follow the edges 01, 12, 23, 30; node 8 of quad9 is the centre",
    )


def cell_types_tetra():
    return _family(
        "tetra",
        ["tetra", "tetra10"],
        640,
        340,
        150,
        True,
        "tetra10 mid-edges follow the edges 01, 12, 02, 03, 13, 23 (VTK)",
    )


def cell_types_hexahedron():
    return _family(
        "hexahedron",
        ["hexahedron", "hexahedron20", "hexahedron24", "hexahedron27"],
        1000,
        380,
        140,
        True,
        "edges: bottom ring 8-11, top ring 12-15, verticals 16-19; faces 20-25: x-min, x-max, y-min, y-max, bottom, top; 26 body",
    )


def cell_types_wedge():
    return _family(
        "wedge",
        ["wedge", "wedge15", "wedge18"],
        800,
        380,
        140,
        True,
        "edges: bottom triangle 6-8, top triangle 9-11, verticals 12-14; wedge18 quad-face centres 15-17",
    )


def cell_types_pyramid():
    return _family(
        "pyramid",
        ["pyramid", "pyramid13", "pyramid14"],
        800,
        380,
        140,
        True,
        "edges: base ring 5-8, apex edges 9-12; pyramid14 adds the base centre 13",
    )


FIGURES = {
    "cell_types_line": cell_types_line,
    "cell_types_triangle": cell_types_triangle,
    "cell_types_quad": cell_types_quad,
    "cell_types_tetra": cell_types_tetra,
    "cell_types_hexahedron": cell_types_hexahedron,
    "cell_types_wedge": cell_types_wedge,
    "cell_types_pyramid": cell_types_pyramid,
}

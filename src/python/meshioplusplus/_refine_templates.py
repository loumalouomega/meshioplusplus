"""Per-(cell type, split-edge mask) subdivision templates for :mod:`._refine`.

The numpy twin of ``src/cpp/include/meshioplusplus/detail/refine_templates.hpp``
and its ``.cpp``. Keeping it in its own module mirrors that split and keeps the
tables out of the refinement loop's way; ``tests/python/test_refine.py`` pins
every entry against ``_core.refine_mask_table``, so the two cannot drift.

A cell's state is a bitmask over :data:`EDGES`: bit ``k`` is set iff edge ``k``
carries a new mid-edge node. Local node ids address the same flat per-cell space
:mod:`._refine` builds — corners, then one node per edge, then one per quad face,
then the body — which coincides with each type's full-Lagrange numbering
(``line3``/``triangle6``/``quad9``/``tetra10``/``wedge18``/``hexahedron27``).

Not every mask has a same-type subdivision. The **admissible** masks per type are
closed under intersection and contain the full mask, which is exactly what makes
"the smallest admissible superset" well defined; :func:`promote_mask` is that
closure operator. See the C++ header for the full argument, including why a
quadrilateral with one or three split edges provably has no all-quad
subdivision.
"""

from __future__ import annotations

from itertools import permutations

# --------------------------------------------------------------------------- #
# per-type geometry                                                            #
# --------------------------------------------------------------------------- #
EDGES = {
    "line": [(0, 1)],
    "triangle": [(0, 1), (1, 2), (2, 0)],
    "quad": [(0, 1), (1, 2), (2, 3), (3, 0)],
    "tetra": [(0, 1), (1, 2), (0, 2), (0, 3), (1, 3), (2, 3)],
    "hexahedron": [
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 0),
        (4, 5),
        (5, 6),
        (6, 7),
        (7, 4),
        (0, 4),
        (1, 5),
        (2, 6),
        (3, 7),
    ],
    "wedge": [
        (0, 1),
        (1, 2),
        (2, 0),
        (3, 4),
        (4, 5),
        (5, 3),
        (0, 3),
        (1, 4),
        (2, 5),
    ],
}

QUAD_FACES = {
    "line": [],
    "triangle": [],
    "tetra": [],
    "quad": [(0, 1, 2, 3)],
    # Row k -> hexahedron27 node 20+k, in _skin.py's _CELL_FACES order
    # (vtkTriQuadraticHexahedron): x-min, x-max, y-min, y-max, bottom, top.
    # Reordered in v9.9.0 with the C++ twin detail/cell_subdivision.cpp; see
    # detail/cell_faces.cpp for the derivation.
    "hexahedron": [
        (0, 4, 7, 3),
        (1, 2, 6, 5),
        (0, 1, 5, 4),
        (3, 7, 6, 2),
        (0, 1, 2, 3),
        (4, 5, 6, 7),
    ],
    "wedge": [(0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)],
}

NUM_CORNERS = {
    "line": 2,
    "triangle": 3,
    "quad": 4,
    "tetra": 4,
    "wedge": 6,
    "hexahedron": 8,
}

SUPPORTED_TYPES = tuple(NUM_CORNERS)

FULL_MASK = {t: (1 << len(EDGES[t])) - 1 for t in SUPPORTED_TYPES}


def _edge_index(cell_type, a, b):
    """The index of edge ``(a, b)`` in this type's :data:`EDGES` row."""
    for k, (u, v) in enumerate(EDGES[cell_type]):
        if (u, v) == (a, b) or (u, v) == (b, a):
            return k
    raise KeyError((cell_type, a, b))


def face_edge_mask(cell_type, face):
    """The bits of the four edges bounding a quad face.

    The face carries a centre node exactly when all of them are split — a
    *derived* rule, so two cells sharing the face cannot disagree about it.
    """
    mask = 0
    for i in range(4):
        mask |= 1 << _edge_index(cell_type, face[i], face[(i + 1) % 4])
    return mask


# --------------------------------------------------------------------------- #
# the tables                                                                   #
# --------------------------------------------------------------------------- #
# One entry per admissible mask: (children, children_alt, tie_a, tie_b).
# `children_alt` is empty unless the mask leaves a quadrilateral remnant needing
# a diagonal; then use it when the GLOBAL node id at `tie_b` is smaller than the
# one at `tie_a`, so that the cell across the affected face resolves it the same
# way. Both variants always have the same child count.


def _line_table():
    return {0b0: (((0, 1),), (), 0, 0), 0b1: (((0, 2), (2, 1)), (), 0, 0)}


def _triangle_table():
    # triangle6 layout: 3 = m(0,1), 4 = m(1,2), 5 = m(2,0).
    return {
        0b000: (((0, 1, 2),), (), 0, 0),
        # One split edge: cut from its midpoint to the opposite corner.
        0b001: (((0, 3, 2), (3, 1, 2)), (), 0, 0),
        0b010: (((0, 1, 4), (0, 4, 2)), (), 0, 0),
        0b100: (((0, 1, 5), (5, 1, 2)), (), 0, 0),
        # Two split edges sharing corner b: b becomes its own triangle and the
        # remnant quadrilateral (a, m_ab, m_bc, c) takes the diagonal starting at
        # whichever of a, c has the smaller global id. The first variant is
        # "from a"; tie_a = a, tie_b = c.
        0b011: (
            ((3, 1, 4), (0, 3, 4), (0, 4, 2)),
            ((3, 1, 4), (0, 3, 2), (3, 4, 2)),
            0,
            2,
        ),
        0b110: (
            ((4, 2, 5), (1, 4, 5), (1, 5, 0)),
            ((4, 2, 5), (1, 4, 0), (4, 5, 0)),
            1,
            0,
        ),
        0b101: (
            ((5, 0, 3), (2, 5, 3), (2, 3, 1)),
            ((5, 0, 3), (2, 5, 1), (5, 3, 1)),
            2,
            1,
        ),
        0b111: (((0, 3, 5), (3, 1, 4), (5, 4, 2), (3, 4, 5)), (), 0, 0),
    }


def _quad_table():
    # quad9 layout: 4..7 = edge mids (edges 01, 12, 23, 30), 8 = face centre.
    # Only the two OPPOSITE pairs are admissible besides none and all: a
    # quadrangulation of an n-gon satisfies 4Q = B + 2I, so the pentagon left by
    # one split edge and the heptagon left by three cannot be filled with
    # quadrilaterals at any number of interior nodes.
    return {
        0b0000: (((0, 1, 2, 3),), (), 0, 0),
        0b0101: (((0, 4, 6, 3), (4, 1, 2, 6)), (), 0, 0),
        0b1010: (((0, 1, 5, 7), (7, 5, 2, 3)), (), 0, 0),
        0b1111: (((0, 4, 8, 7), (4, 1, 5, 8), (8, 5, 2, 6), (7, 8, 6, 3)), (), 0, 0),
    }


def _wedge_table():
    # wedge18 layout: 6..8 bottom-triangle mids, 9..11 top-triangle mids,
    # 12..14 vertical mids, 15..17 quad-face centres. The six triangle edges form
    # one class and the three verticals another.
    tri = 0b000111111
    vert = 0b111000000
    return {
        0: (((0, 1, 2, 3, 4, 5),), (), 0, 0),
        tri: (
            (
                (0, 6, 8, 3, 9, 11),
                (6, 1, 7, 9, 4, 10),
                (8, 7, 2, 11, 10, 5),
                (6, 7, 8, 9, 10, 11),
            ),
            (),
            0,
            0,
        ),
        vert: (((0, 1, 2, 12, 13, 14), (12, 13, 14, 3, 4, 5)), (), 0, 0),
        tri
        | vert: (
            (
                (0, 6, 8, 12, 15, 17),
                (6, 1, 7, 15, 13, 16),
                (8, 7, 2, 17, 16, 14),
                (6, 7, 8, 15, 16, 17),
                (12, 15, 17, 3, 9, 11),
                (15, 13, 16, 9, 4, 10),
                (17, 16, 14, 11, 10, 5),
                (15, 16, 17, 9, 10, 11),
            ),
            (),
            0,
            0,
        ),
    }


def _hexahedron_table():
    # hexahedron27 layout: 8..19 edge mids (bottom ring, top ring, verticals),
    # 20..25 face centres, 26 body. The face-centre indices are QUAD_FACES
    # ["hexahedron"] row order, which since v9.9.0 is _skin.py's _CELL_FACES
    # order (vtkTriQuadraticHexahedron): 20 = x-min, 21 = x-max, 22 = y-min,
    # 23 = y-max, 24 = bottom, 25 = top. The absolute indices below were
    # permuted by the same 3-cycle (20->22, 22->23, 23->20) in that release;
    # kept in exact sync with the C++ twin detail/refine_templates.cpp, which
    # test_refine.py::test_mask_tables_match_the_cpp_core enforces.
    #
    # The twelve edges fall into three parallel classes; the admissible masks
    # are their eight unions, so a refinement travels through one dual sheet
    # rather than the whole block.
    x = 0b000001010101
    y = 0b000010101010
    z = 0b111100000000
    return {
        0: (((0, 1, 2, 3, 4, 5, 6, 7),), (), 0, 0),
        x: (((0, 8, 10, 3, 4, 12, 14, 7), (8, 1, 2, 10, 12, 5, 6, 14)), (), 0, 0),
        y: (((0, 1, 9, 11, 4, 5, 13, 15), (11, 9, 2, 3, 15, 13, 6, 7)), (), 0, 0),
        z: (((0, 1, 2, 3, 16, 17, 18, 19), (16, 17, 18, 19, 4, 5, 6, 7)), (), 0, 0),
        x
        | y: (
            (
                (0, 8, 24, 11, 4, 12, 25, 15),
                (8, 1, 9, 24, 12, 5, 13, 25),
                (24, 9, 2, 10, 25, 13, 6, 14),
                (11, 24, 10, 3, 15, 25, 14, 7),
            ),
            (),
            0,
            0,
        ),
        x
        | z: (
            (
                (0, 8, 10, 3, 16, 22, 23, 19),
                (8, 1, 2, 10, 22, 17, 18, 23),
                (16, 22, 23, 19, 4, 12, 14, 7),
                (22, 17, 18, 23, 12, 5, 6, 14),
            ),
            (),
            0,
            0,
        ),
        y
        | z: (
            (
                (0, 1, 9, 11, 16, 17, 21, 20),
                (11, 9, 2, 3, 20, 21, 18, 19),
                (16, 17, 21, 20, 4, 5, 13, 15),
                (20, 21, 18, 19, 15, 13, 6, 7),
            ),
            (),
            0,
            0,
        ),
        x
        | y
        | z: (
            (
                (0, 8, 24, 11, 16, 22, 26, 20),
                (8, 1, 9, 24, 22, 17, 21, 26),
                (11, 24, 10, 3, 20, 26, 23, 19),
                (24, 9, 2, 10, 26, 21, 18, 23),
                (16, 22, 26, 20, 4, 12, 25, 15),
                (22, 17, 21, 26, 12, 5, 13, 25),
                (20, 26, 23, 19, 15, 25, 14, 7),
                (26, 21, 18, 23, 25, 13, 6, 14),
            ),
            (),
            0,
            0,
        ),
    }


def _even_permutations():
    """The 12 even permutations of ``(0, 1, 2, 3)``, lexicographically.

    Even, so each preserves a tetrahedron's orientation — which is what lets a
    representative's positively-wound children be relabelled rather than
    re-derived. Same order as the C++ generator, so both pick the same
    representative for each mask.
    """
    out = []
    for perm in permutations(range(4)):
        inversions = sum(
            1 for i in range(4) for j in range(i + 1, 4) if perm[i] > perm[j]
        )
        if inversions % 2 == 0:
            out.append(perm)
    return out


def _tetra_table():
    # tetra10 layout: 4 = m(0,1), 5 = m(1,2), 6 = m(0,2), 7 = m(0,3),
    # 8 = m(1,3), 9 = m(2,3).
    def bit(a, b):
        return 1 << _edge_index("tetra", a, b)

    # Six orbit representatives; everything else in the admissible set is one of
    # these relabelled by an even permutation.
    reps = [
        (0, (((0, 1, 2, 3),), (), 0, 0)),
        # One split edge: cut the tetrahedron in half through m(0,1).
        (bit(0, 1), (((0, 4, 2, 3), (4, 1, 2, 3)), (), 0, 0)),
        # Two ADJACENT split edges meeting at corner 1: corner 1 becomes its own
        # tetrahedron and the remnant quadrilateral on face (0,1,2) is coned to
        # corner 3 across whichever diagonal starts at the smaller-id surviving
        # corner. The first variant is "from corner 0".
        (
            bit(0, 1) | bit(1, 2),
            (
                ((4, 1, 5, 3), (0, 4, 5, 3), (0, 5, 2, 3)),
                ((4, 1, 5, 3), (0, 4, 2, 3), (4, 5, 2, 3)),
                0,
                2,
            ),
        ),
        # Two OPPOSITE split edges: four tetrahedra around the interior segment
        # m(0,1)-m(2,3), one per edge of the equator.
        (
            bit(0, 1) | bit(2, 3),
            (((4, 9, 2, 0), (4, 9, 1, 2), (4, 9, 3, 1), (4, 9, 0, 3)), (), 0, 0),
        ),
        # Three split edges sharing face (0,1,2): that face takes the 1-to-4
        # triangle split and each sub-triangle is coned to corner 3.
        (
            bit(0, 1) | bit(1, 2) | bit(0, 2),
            (((0, 4, 6, 3), (4, 1, 5, 3), (6, 5, 2, 3), (4, 5, 6, 3)), (), 0, 0),
        ),
        # Every edge split: four corner tetrahedra, then the residual octahedron
        # split along the fixed interior diagonal 4-9, ring 6->7->8->5.
        (
            FULL_MASK["tetra"],
            (
                (
                    (0, 4, 6, 7),
                    (4, 1, 5, 8),
                    (6, 5, 2, 9),
                    (7, 8, 9, 3),
                    (4, 9, 6, 7),
                    (4, 9, 7, 8),
                    (4, 9, 8, 5),
                    (4, 9, 5, 6),
                ),
                (),
                0,
                0,
            ),
        ),
    ]

    def map_local(sigma, local):
        if local < 4:
            return sigma[local]
        a, b = EDGES["tetra"][local - 4]
        return 4 + _edge_index("tetra", sigma[a], sigma[b])

    def map_mask(sigma, mask):
        out = 0
        for k, (a, b) in enumerate(EDGES["tetra"]):
            if mask & (1 << k):
                out |= 1 << _edge_index("tetra", sigma[a], sigma[b])
        return out

    table = {}
    for rep_mask, (children, children_alt, tie_a, tie_b) in reps:
        for sigma in _even_permutations():
            mask = map_mask(sigma, rep_mask)
            if mask in table:
                continue  # first permutation reaching a mask wins -- deterministic
            table[mask] = (
                tuple(tuple(map_local(sigma, n) for n in c) for c in children),
                tuple(tuple(map_local(sigma, n) for n in c) for c in children_alt),
                sigma[tie_a] if children_alt else 0,
                sigma[tie_b] if children_alt else 0,
            )
    return table


def _build():
    tables = {
        "line": _line_table(),
        "triangle": _triangle_table(),
        "quad": _quad_table(),
        "tetra": _tetra_table(),
        "wedge": _wedge_table(),
        "hexahedron": _hexahedron_table(),
    }
    promote = {}
    for cell_type, table in tables.items():
        full = FULL_MASK[cell_type]
        admissible = sorted(table)
        rows = {}
        # Straight from the definition: the intersection of every admissible mask
        # containing m. Closure under intersection is what makes that itself
        # admissible, hence the *least* admissible superset.
        for m in range(full + 1):
            acc = full
            for a in admissible:
                if a & m == m:
                    acc &= a
            rows[m] = acc
        promote[cell_type] = rows
    return tables, promote


TABLES, _PROMOTE = _build()


def template(cell_type, mask):
    """The subdivision for one ``(type, mask)`` pair, or ``None`` if inadmissible."""
    return TABLES.get(cell_type, {}).get(mask)


def promote_mask(cell_type, mask, propagate=False):
    """The smallest admissible superset of ``mask``.

    Monotone and idempotent, so iterating it over a mesh converges to a unique
    fixed point regardless of the order cells are visited in. With ``propagate``
    any non-empty mask goes straight to the full split instead.
    """
    full = FULL_MASK[cell_type]
    mask &= full
    if propagate:
        return full if mask else 0
    return _PROMOTE[cell_type][mask]


def closure_from_name(name):
    """Parse a closure name into its canonical spelling.

    Returns one of ``"redgreen"``, ``"propagate"`` or ``"balanced"``.
    """
    if name in ("", "redgreen", "red-green", "green"):
        return "redgreen"
    if name in ("propagate", "red"):
        return "propagate"
    if name in ("balanced", "2:1"):
        return "balanced"
    raise ValueError(
        f"refine: unknown closure '{name}' "
        "(expected 'redgreen'/'green', 'propagate' or 'balanced')"
    )


_COMPARE = {
    "<": lambda v, r: v < r,
    "lt": lambda v, r: v < r,
    "<=": lambda v, r: v <= r,
    "le": lambda v, r: v <= r,
    ">": lambda v, r: v > r,
    "gt": lambda v, r: v > r,
    ">=": lambda v, r: v >= r,
    "ge": lambda v, r: v >= r,
    "==": lambda v, r: v == r,
    "=": lambda v, r: v == r,
    "eq": lambda v, r: v == r,
    "!=": lambda v, r: v != r,
    "ne": lambda v, r: v != r,
}


def compare_from_name(name):
    """Parse a comparison operator into a callable ``(value, rhs) -> bool``."""
    try:
        return _COMPARE[name]
    except KeyError:
        raise ValueError(
            f"refine: unknown comparison '{name}' "
            "(expected '<', '<=', '>', '>=', '==' or '!=')"
        ) from None

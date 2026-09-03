"""Refinement-family figures: admissible masks, green undo, slice vs crop."""

from diaglib import palette as P
from diaglib.cellart import draw_cell, fit_origin, project_all
from diaglib.svg import Canvas
from figures.cells import tables


def _popcount(m):
    return bin(m).count("1")


def _flat(pts3, scale, origin):
    return [(origin[0] + scale * x, origin[1] - scale * y) for x, y, _ in pts3]


def _draw_2d_split(c, t, cell_type, mask, box, scale, colour):
    """A 2-D parent cell with its split edges and the template's children."""
    quad_type = {"triangle": "triangle6", "quad": "quad9"}[cell_type]
    nodes = t.node_coordinates(quad_type)
    origin = fit_origin(nodes[: t.num_nodes[cell_type]], scale, box, flat=True)
    pts = _flat(nodes, scale, origin)
    nc = t.num_nodes[cell_type]
    edges = t.refine.EDGES[cell_type]
    tpl = t.refine.template(cell_type, mask)
    if tpl is not None:
        for child in tpl[0]:
            c.polygon(
                [pts[i] for i in child],
                fill=colour,
                fill_opacity=0.28 if mask else 0.0,
                stroke=colour if mask else P.INK,
                sw=1.4,
            )
    c.polygon(pts[:nc], fill="none", stroke=P.INK, sw=1.8)
    for k, (a, b) in enumerate(edges):
        if mask & (1 << k):
            c.line(*pts[a], *pts[b], stroke=P.CORE, sw=3.0)
            m = t.edge_midpoint_index(quad_type, a, b)
            c.circle(*pts[m], 4, fill=P.WHITE, stroke=P.CORE, sw=2)
    for i in range(nc):
        c.circle(*pts[i], 3.5, fill=P.INK)
    return tpl


def refine_masks_2d():
    c = Canvas(
        960,
        560,
        "Admissible refinement masks for triangles and quadrilaterals, and how an inadmissible quad mask is promoted",
    )
    t = tables()
    # -- triangles: every mask is admissible
    c.caption(
        480,
        30,
        "triangle: every set of split edges has a template; the full split is red, a partial one is a green closure",
    )
    masks = sorted(t.refine.TABLES["triangle"], key=lambda m: (_popcount(m), m))
    panel = 900 / len(masks)
    for k, mask in enumerate(masks):
        x0 = 30 + k * panel
        colour = P.ERROR if mask == t.refine.FULL_MASK["triangle"] else P.WASM
        tpl = _draw_2d_split(
            c, t, "triangle", mask, (x0 + 8, 50, panel - 16, 110), 85, colour
        )
        n = len(tpl[0])
        c.label(
            x0 + panel / 2,
            185,
            f"{_popcount(mask)} split → {n} child{'ren' if n > 1 else ''}",
            size=P.SIZE_SMALL,
        )
        c.label(x0 + panel / 2, 200, f"mask {mask:03b}", mono=True)
    # -- quads: opposite pairs only
    c.caption(
        480,
        240,
        "quad: only 'none', an opposite pair, or all four are admissible; anything else is promoted to the smallest admissible superset",
    )
    admissible = sorted(t.refine.TABLES["quad"], key=lambda m: (_popcount(m), m))
    x = 30
    for mask in admissible:
        colour = P.ERROR if mask == t.refine.FULL_MASK["quad"] else P.WASM
        tpl = _draw_2d_split(c, t, "quad", mask, (x, 260, 110, 110), 85, colour)
        n = len(tpl[0])
        c.label(
            x + 55,
            395,
            f"{_popcount(mask)} split → {n} child{'ren' if n > 1 else ''}",
            size=P.SIZE_SMALL,
        )
        c.label(x + 55, 410, f"mask {mask:04b}", mono=True)
        x += 120
    c.line(x + 10, 260, x + 10, 420, stroke=P.HAIRLINE)
    x += 40
    for mask in (0b0001, 0b0011):
        promoted = t.refine.promote_mask("quad", mask)
        _draw_2d_split(c, t, "quad", mask, (x, 260, 70, 70), 60, P.MUTED)
        c.label(x + 35, 350, f"{mask:04b}", mono=True)
        c.arrow(x + 78, 295, x + 100, 295, stroke=P.INK_2, sw=1.6)
        colour = P.ERROR if promoted == t.refine.FULL_MASK["quad"] else P.WASM
        _draw_2d_split(c, t, "quad", promoted, (x + 106, 260, 70, 70), 60, colour)
        c.label(x + 141, 350, f"{promoted:04b}", mono=True)
        c.label(x + 90, 372, "promoted", size=P.SIZE_SMALL)
        x += 200
    c.label(x + 70, 300, "three split edges", size=P.SIZE_SMALL)
    c.label(x + 70, 316, "promote to 1111 too")
    c.legend(
        30,
        440,
        [
            (P.ERROR, "red: fully split (a refinement)"),
            (P.WASM, "green: a closure, inherits its parent's refine:level"),
            (P.CORE, "split edge, with its new mid-edge node"),
        ],
    )
    lines = [
        "the whole operation is driven by one set — which EDGES carry a node — and everything else is derived from it: two cells sharing an edge read the same bit, so they cannot disagree,",
        "and a quad face gets a centre node exactly when all four of its edges are split. Each type's admissible masks are closed under intersection (a Moore family), so 'promote to the",
        "smallest admissible superset' is a monotone idempotent closure operator: iterating it over the mesh converges to one fixed point whatever the traversal order.",
        "A quad with one or three split edges has no all-quad subdivision at any number of interior nodes (4Q = B + 2I), which is why the opposite pair is the promotion, not a triangle fan.",
    ]
    for k, s in enumerate(lines):
        c.label(480, 468 + 17 * k, s)
    return c.render()


def refine_masks_tetra():
    c = Canvas(
        960,
        400,
        "Admissible tetrahedron refinement masks, by number of split edges, and the one three-edge configuration excluded on purpose",
    )
    t = tables()
    edges = t.refine.EDGES["tetra"]
    admissible = t.refine.TABLES["tetra"]
    nodes10 = t.node_coordinates("tetra10")
    reps = []
    seen = set()
    for mask in sorted(admissible, key=lambda m: (_popcount(m), m)):
        pc = _popcount(mask)
        if pc in (1, 3, 6) and pc in seen:
            continue
        if pc == 2:
            # keep one adjacent and one opposite pair
            a, b = [edges[k] for k in range(6) if mask & (1 << k)]
            kind = "opposite" if not set(a) & set(b) else "adjacent"
            if ("2", kind) in seen:
                continue
            seen.add(("2", kind))
        elif pc in (0, 4, 5):
            continue
        seen.add(pc)
        reps.append(mask)
    excluded = sum(
        1 << k for k, e in enumerate(edges) if 0 in e
    )  # the three edges at corner 0
    assert excluded not in admissible
    panels = [(m, False) for m in reps] + [(excluded, True)]
    panel_w = 900 / len(panels)
    for k, (mask, bad) in enumerate(panels):
        x0 = 30 + k * panel_w
        box = (x0 + 10, 40, panel_w - 20, 200)
        colour = (
            P.ERROR
            if bad
            else (P.ERROR if mask == t.refine.FULL_MASK["tetra"] else P.WASM)
        )
        pts = draw_cell(
            c,
            t,
            "tetra",
            None,
            120,
            labels=False,
            box=box,
            fill=colour,
            fill_opacity=0.10,
        )
        origin = fit_origin(t.node_coordinates("tetra"), 120, box)
        p10 = project_all(nodes10, 120, origin)
        for j, (a, b) in enumerate(edges):
            if mask & (1 << j):
                c.line(*pts[a], *pts[b], stroke=P.CORE, sw=3)
                c.circle(*p10[4 + j], 4, fill=P.WHITE, stroke=P.CORE, sw=2)
        if mask == t.refine.FULL_MASK["tetra"]:
            c.line(*p10[4], *p10[9], stroke=P.CABI, sw=2, dash="5 3")
        for i in range(4):
            c.circle(*pts[i], 3.5, fill=P.INK)
        if bad:
            c.text(
                x0 + panel_w / 2,
                268,
                "3 edges at one corner",
                size=P.SIZE_SMALL,
                weight="700",
                fill=P.ERROR,
            )
            c.label(x0 + panel_w / 2, 284, "not admissible", fill=P.ERROR)
            c.label(
                x0 + panel_w / 2,
                300,
                f"promoted to {t.refine.promote_mask('tetra', mask):06b}",
                mono=True,
            )
        else:
            n = len(admissible[mask][0])
            pc = _popcount(mask)
            a_b = [edges[j] for j in range(6) if mask & (1 << j)]
            kind = ""
            if pc == 2:
                kind = " (opposite)" if not set(a_b[0]) & set(a_b[1]) else " (adjacent)"
            if pc == 3:
                kind = " (one face)"
            c.text(
                x0 + panel_w / 2,
                268,
                f"{pc} split → {n} children",
                size=P.SIZE_SMALL,
                weight="700",
            )
            c.label(x0 + panel_w / 2, 284, kind.strip(" ()") if kind else " ")
            c.label(x0 + panel_w / 2, 300, f"mask {mask:06b}", mono=True)
    c.legend(
        30,
        336,
        [
            (P.CORE, "split edge + mid-edge node"),
            (P.WASM, "green closure"),
            (P.ERROR, "full split, or excluded"),
            (P.CABI, "the fixed interior diagonal 4–9 of the full split"),
        ],
    )
    c.label(
        480,
        362,
        f"the {len(admissible)} admissible tetra masks are generated from six representatives by the 12 even corner permutations (orientation-preserving, so children stay positively wound);",
    )
    c.label(
        480,
        379,
        "three edges meeting at a corner are excluded deliberately: they would put the ambiguous two-edge case on all three incident faces at once",
    )
    return c.render()


def undo_green_timeline():
    c = Canvas(
        960,
        390,
        "undo_green restores the green closure cells of a selective refinement to their parents, keeping the red refinement",
    )
    cell = 44.0

    def grid(x0, y0, cells, title, subtitle):
        c.caption(x0 + 1.5 * cell, y0 - 26, title)
        c.label(x0 + 1.5 * cell, y0 - 10, subtitle)
        for i, j, w, h, colour, label, sub in cells:
            x = x0 + i * cell
            y = y0 + (3 - j - h) * cell
            c.rrect(
                x,
                y,
                w * cell,
                h * cell,
                fill=colour if colour else P.WHITE,
                fill_opacity=0.28 if colour else None,
                stroke=P.INK,
                rx=0,
                sw=1.2,
            )
            if label is not None:
                c.text(
                    x + w * cell / 2,
                    y + h * cell / 2 + (4 if not sub else -2),
                    str(label),
                    size=P.SIZE_SMALL,
                    weight="700",
                    mono=True,
                )
            if sub:
                c.label(x + w * cell / 2, y + h * cell / 2 + 11, sub, size=9)
        c.rrect(x0, y0, 3 * cell, 3 * cell, fill="none", stroke=P.INK, rx=0, sw=2)

    # coarse: 3x3 quads, cell_id = global index
    coarse = [
        (i, j, 1, 1, None, j * 3 + i, "level 0") for j in range(3) for i in range(3)
    ]
    grid(40, 100, coarse, "coarse", "refine:cell_id = global index")
    # fine: cell 4 split in 4 (red), edge neighbours 1,3,5,7 promoted to the
    # opposite pair (green, 2 children each), corners untouched
    fine = []
    for j in range(3):
        for i in range(3):
            k = j * 3 + i
            if k == 4:
                for dj in (0, 1):
                    for di in (0, 1):
                        fine.append(
                            (i + 0.5 * di, j + 0.5 * dj, 0.5, 0.5, P.ERROR, None, None)
                        )
            elif k in (
                1,
                7,
            ):  # above/below: split left-right? no: shared edge is horizontal -> opposite pair = top+bottom edges -> children stacked? A split horizontal edge pairs with the opposite horizontal edge, so the children are left/right halves.
                fine.append((i, j, 0.5, 1, P.WASM, None, None))
                fine.append((i + 0.5, j, 0.5, 1, P.WASM, None, None))
            elif k in (3, 5):
                fine.append((i, j, 1, 0.5, P.WASM, None, None))
                fine.append((i, j + 0.5, 1, 0.5, P.WASM, None, None))
            else:
                fine.append((i, j, 1, 1, None, None, None))
    grid(
        340,
        100,
        fine,
        "fine = refine(coarse, cells=[4], record_hierarchy, record_levels)",
        "red: level 1, parent 4 · green: level 0, parent 1/3/5/7 · white: untouched",
    )
    # restored
    restored = []
    for j in range(3):
        for i in range(3):
            k = j * 3 + i
            if k == 4:
                for dj in (0, 1):
                    for di in (0, 1):
                        restored.append(
                            (i + 0.5 * di, j + 0.5 * dj, 0.5, 0.5, P.ERROR, None, None)
                        )
            else:
                restored.append((i, j, 1, 1, None, k, None))
    grid(
        700,
        100,
        restored,
        "undo_green(coarse, fine)",
        "green groups collapse to the parent's own row",
    )
    c.arrow(190, 166, 330, 166, stroke=P.INK_2)
    c.label(260, 156, "refine")
    c.arrow(490, 166, 690, 166, stroke=P.INK_2)
    c.label(590, 156, "undo_green")
    c.legend(
        40,
        275,
        [
            (P.ERROR, "red: refine:level = parent's + 1, kept"),
            (P.WASM, "green: refine:level = parent's, replaced by the parent's row"),
            (P.MUTED, "untouched: cell_id == parent_id"),
        ],
    )
    c.label(
        480,
        305,
        "classification is per sibling group (cells sharing one refine:parent_id, resolved against coarse's refine:cell_id):",
    )
    c.label(
        480,
        321,
        "a parent's mask is uniform across its children, so red and green never mix in a group; points are never pruned or renumbered,",
    )
    c.label(
        480,
        337,
        "so a parent's connectivity row from coarse is valid in fine as-is: a lookup and substitution, not a reconstruction.",
    )
    c.label(
        480,
        353,
        "A hanging node left by closure='balanced' is not a green cell and is untouched here; only single-pass (levels=1) hierarchies are accepted.",
    )
    return c.render()


def slice_vs_crop():
    c = Canvas(
        960,
        430,
        "crop keeps whole cells on one side of a plane or box; slice computes the intersection and lowers the dimension",
    )
    cell = 34.0
    # an L-shaped quad mesh
    cells = [(i, j) for i in range(6) for j in range(5) if not (i >= 3 and j >= 2)]

    def draw_mesh(x0, y0, keep=None, fade=None):
        for i, j in cells:
            x = x0 + i * cell
            y = y0 + (4 - j) * cell
            kept = keep is None or keep(i, j)
            c.rrect(
                x,
                y,
                cell,
                cell,
                fill=P.CORE if kept else P.PAPER,
                fill_opacity=0.14 if kept else None,
                stroke=P.INK if kept else P.MUTED,
                rx=0,
                sw=1.2 if kept else 0.8,
                dash=None if kept else "3 3",
            )

    def inside(x, y):
        return 1.3 <= x <= 4.6 and 0.6 <= y <= 3.4

    def corners(i, j):
        return [(i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)]

    # panel 1: crop all
    x0, y0 = 40, 90
    draw_mesh(x0, y0, lambda i, j: all(inside(*p) for p in corners(i, j)))
    c.rrect(
        x0 + 1.3 * cell,
        y0 + (5 - 3.4) * cell,
        3.3 * cell,
        2.8 * cell,
        fill="none",
        stroke=P.DATA,
        rx=0,
        sw=2.2,
        dash="6 3",
    )
    c.caption(x0 + 3 * cell, 62, 'crop(bbox, mode="all")')
    c.label(
        x0 + 3 * cell, y0 + 5 * cell + 20, "keep a cell only if every node is inside"
    )
    # panel 2: crop any
    x0 = 340
    draw_mesh(x0, y0, lambda i, j: any(inside(*p) for p in corners(i, j)))
    c.rrect(
        x0 + 1.3 * cell,
        y0 + (5 - 3.4) * cell,
        3.3 * cell,
        2.8 * cell,
        fill="none",
        stroke=P.DATA,
        rx=0,
        sw=2.2,
        dash="6 3",
    )
    c.caption(x0 + 3 * cell, 62, 'crop(bbox, mode="any")')
    c.label(x0 + 3 * cell, y0 + 5 * cell + 20, "keep a cell if any node is inside")
    # panel 3: slice by a line (the 2-D plane): the section is a set of segments
    x0 = 640
    draw_mesh(x0, y0, lambda i, j: False)

    # plane: y = 0.55 x + 0.9  (in cell units)
    def fy(x):
        return 0.55 * x + 0.9

    px0, px1 = -0.3, 6.3
    c.line(
        x0 + px0 * cell,
        y0 + (5 - fy(px0)) * cell,
        x0 + px1 * cell,
        y0 + (5 - fy(px1)) * cell,
        stroke=P.DATA,
        sw=2.2,
        dash="6 3",
    )
    # segments where the line crosses cells (the section, one topological dimension down)
    for i, j in cells:
        pts = []
        cs = corners(i, j)
        for k in range(4):
            (xa, ya), (xb, yb) = cs[k], cs[(k + 1) % 4]
            da, db = ya - fy(xa), yb - fy(xb)
            if (da >= 0) != (db >= 0):
                tt = da / (da - db)
                pts.append((xa + tt * (xb - xa), ya + tt * (yb - ya)))
        if len(pts) == 2:
            (xa, ya), (xb, yb) = pts
            c.line(
                x0 + xa * cell,
                y0 + (5 - ya) * cell,
                x0 + xb * cell,
                y0 + (5 - yb) * cell,
                stroke=P.ERROR,
                sw=3.2,
            )
            for xx, yy in pts:
                c.circle(
                    x0 + xx * cell,
                    y0 + (5 - yy) * cell,
                    3,
                    fill=P.WHITE,
                    stroke=P.ERROR,
                    sw=1.6,
                )
    c.caption(x0 + 3 * cell, 62, "slice(origin, normal)")
    c.label(
        x0 + 3 * cell,
        y0 + 5 * cell + 20,
        "the intersection itself: 2-D cells → line segments",
    )
    c.legend(
        40,
        320,
        [
            (P.CORE, "kept cell: crop returns a mesh of the same dimension"),
            (P.DATA, "the box / the cutting plane"),
            (P.ERROR, "the section, one dimension lower"),
        ],
    )
    c.label(
        480,
        352,
        "crop keeps geometry and connectivity intact (regions remapped, crop:original_*_id on request);",
    )
    c.label(
        480,
        368,
        "slice makes new points at the crossings, interpolates point_data there and replicates each parent's cell_data — a volume yields triangles and quads, a surface yields lines.",
    )
    c.label(
        480,
        384,
        "Both take the plane as (origin, normal); crop also takes a bounding box or a cell_data predicate (crop_predicate), and mode selects the all-nodes / any-node rule.",
    )
    return c.render()


FIGURES = {
    "refine_masks_2d": refine_masks_2d,
    "refine_masks_tetra": refine_masks_tetra,
    "undo_green_timeline": undo_green_timeline,
    "slice_vs_crop": slice_vs_crop,
}

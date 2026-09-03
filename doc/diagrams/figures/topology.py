"""Algorithm figures: skin hashing, flips, subdivide, agglomerate, marching, BCC."""

from diaglib import palette as P
from diaglib.cellart import draw_solid, fit_origin, oriented_tet, project_all
from diaglib.project import centroid, newell_normal, project
from diaglib.svg import Canvas
from figures.cells import tables


def _norm(v):
    n = sum(x * x for x in v) ** 0.5
    return tuple(x / n for x in v)


# --------------------------------------------------------------------------- #
# extract_skin: face hashing                                                   #
# --------------------------------------------------------------------------- #
def extract_skin_hashing():
    c = Canvas(
        900,
        430,
        "Skin extraction by face hashing: a face whose sorted-corner key occurs once is on the boundary",
    )
    t = tables()
    # Two unit hexahedra sharing the x = 1 face; VTK corner order.
    pts = [
        (0, 0, 0),
        (1, 0, 0),
        (1, 1, 0),
        (0, 1, 0),
        (0, 0, 1),
        (1, 0, 1),
        (1, 1, 1),
        (0, 1, 1),
        (2, 0, 0),
        (2, 1, 0),
        (2, 0, 1),
        (2, 1, 1),
    ]
    cells = [[0, 1, 2, 3, 4, 5, 6, 7], [1, 8, 9, 2, 5, 10, 11, 6]]
    faces, owners = [], []
    for ci, conn in enumerate(cells):
        for row in t.cell_faces["hexahedron"]:
            faces.append(tuple(conn[i] for i in row[2]))
            owners.append(ci)
    counts = {}
    for f in faces:
        counts.setdefault(tuple(sorted(f)), []).append(f)
    scale = 118
    box = (30, 60, 470, 300)
    origin = fit_origin(pts, scale, box)
    # draw both cells; the shared face is interior so it never shows, but we
    # trace it dashed so the reader sees where the two keys collide.
    fills = {}
    for i, f in enumerate(faces):
        key = tuple(sorted(f))
        if len(counts[key]) == 1:
            fills[i] = (P.FORMATS, 0.10)
    p2 = draw_solid(c, pts, faces, scale, origin, fill=P.FORMATS, face_fills=fills)
    shared = [i for i, f in enumerate(faces) if len(counts[tuple(sorted(f))]) == 2][0]
    c.polygon(
        [p2[j] for j in faces[shared]],
        fill=P.ERROR,
        fill_opacity=0.18,
        stroke=P.ERROR,
        sw=1.4,
        dash="5 3",
    )
    for i, (x, y) in enumerate(p2):
        c.node(x, y, i, r=3.5, label_dx=8, label_dy=-5)
    # outward normal on the top face of the left cell
    top = [i for i, f in enumerate(faces) if owners[i] == 0 and set(f) == {4, 5, 6, 7}][
        0
    ]
    fc = centroid([pts[j] for j in faces[top]])
    n = _norm(newell_normal([pts[j] for j in faces[top]]))
    a = project(fc, scale, origin)
    b = project(tuple(fc[k] + 0.45 * n[k] for k in range(3)), scale, origin)
    c.arrow(a[0], a[1], b[0], b[1], stroke=P.WASM, sw=2.2)
    c.label(b[0] + 6, b[1] - 4, "outward normal", anchor="start", fill=P.WASM)
    c.label(
        265,
        392,
        "two hexahedra sharing the x = 1 face; the shared face (red, dashed) is hashed twice",
    )
    # key table
    x0, y0 = 530, 60
    c.text(x0, y0, "sorted-corner key", size=P.SIZE_SMALL, anchor="start", weight="600")
    c.text(x0 + 190, y0, "count", size=P.SIZE_SMALL, anchor="start", weight="600")
    c.text(x0 + 250, y0, "verdict", size=P.SIZE_SMALL, anchor="start", weight="600")
    c.line(x0, y0 + 6, x0 + 350, y0 + 6, stroke=P.HAIRLINE)
    rows = sorted(counts.items(), key=lambda kv: (-len(kv[1]), kv[0]))
    y = y0 + 24
    for key, occ in rows:
        interior = len(occ) == 2
        colour = P.ERROR if interior else P.FORMATS
        c.text(
            x0,
            y,
            "{" + ", ".join(str(k) for k in key) + "}",
            size=P.SIZE_SMALL,
            anchor="start",
            mono=True,
            fill=P.INK,
        )
        c.text(
            x0 + 190,
            y,
            str(len(occ)),
            size=P.SIZE_SMALL,
            anchor="start",
            mono=True,
            fill=colour,
            weight="700",
        )
        c.text(
            x0 + 250,
            y,
            "interior, cancels" if interior else "boundary, kept",
            size=P.SIZE_SMALL,
            anchor="start",
            fill=colour,
        )
        y += 20
    c.label(
        x0,
        y + 8,
        "keys come from _CELL_FACES / detail/cell_faces.hpp; kept faces are",
        anchor="start",
    )
    c.label(
        x0,
        y + 24,
        "emitted with the winding those tables give, so normals point outward",
        anchor="start",
    )
    return c.render()


# --------------------------------------------------------------------------- #
# optimize_volume: 2-3 and 3-2 flips                                           #
# --------------------------------------------------------------------------- #
def flips_23_32():
    c = Canvas(
        900,
        400,
        "A 2-3 flip replaces two tetrahedra sharing a face with three around the new edge; the 3-2 flip is its inverse",
    )
    a, b, cc = (0.0, 0.0, 0.5), (1.0, 0.0, 0.5), (0.35, 1.0, 0.5)
    p, q = (0.45, 0.35, 1.25), (0.45, 0.35, -0.25)
    pts = [a, b, cc, p, q]
    names = ["a", "b", "c", "p", "q"]
    scale = 150

    def draw(box, tets, highlight_edge=None, highlight_face=None, title=None):
        faces = []
        for tet in tets:
            for f in oriented_tet(pts, *tet):
                faces.append(f)
        # exterior faces only (interior ones cancel exactly as in skin extraction)
        keys = {}
        for i, f in enumerate(faces):
            keys.setdefault(tuple(sorted(f)), []).append(i)
        ext = [faces[v[0]] for k, v in keys.items() if len(v) == 1]
        origin = fit_origin(pts, scale, box)
        p2 = draw_solid(c, pts, ext, scale, origin, fill=P.CORE, fill_opacity=0.08)
        if highlight_face is not None:
            c.polygon(
                [p2[j] for j in highlight_face],
                fill=P.DATA,
                fill_opacity=0.35,
                stroke=P.DATA,
                sw=1.5,
            )
        if highlight_edge is not None:
            i, j = highlight_edge
            c.line(*p2[i], *p2[j], stroke=P.ERROR, sw=2.4)
        # interior edges of the configuration, dotted
        for tet in tets:
            for e in ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)):
                i, j = tet[e[0]], tet[e[1]]
                if not any(set((i, j)) <= set(f) for f in ext):
                    c.line(*p2[i], *p2[j], stroke=P.INK_2, sw=1.2, dash="2 3")
        cx = sum(x for x, _ in p2) / 5
        cy = sum(y for _, y in p2) / 5
        for i, (x, y) in enumerate(p2):
            dx, dy = x - cx, y - cy
            n = (dx * dx + dy * dy) ** 0.5
            c.circle(x, y, 4, fill=P.INK)
            c.text(
                x + 13 * dx / n,
                y + 13 * dy / n + 4,
                names[i],
                size=P.SIZE_SMALL,
                weight="600",
                italic=True,
            )
        if title:
            c.caption(box[0] + box[2] / 2, box[1] - 10, title)

    draw(
        (40, 70, 330, 260),
        [(0, 1, 2, 3), (0, 2, 1, 4)],
        highlight_face=(0, 1, 2),
        title="two tetrahedra sharing the face (a, b, c)",
    )
    draw(
        (530, 70, 330, 260),
        [(0, 1, 3, 4), (1, 2, 3, 4), (2, 0, 3, 4)],
        highlight_edge=(3, 4),
        title="three tetrahedra around the edge p–q",
    )
    c.arrow(392, 170, 508, 170, stroke=P.CORE, sw=2.2)
    c.text(450, 160, "2-3 flip", size=P.SIZE_LABEL, weight="700", fill=P.CORE)
    c.arrow(508, 215, 392, 215, stroke=P.CABI, sw=2.2)
    c.text(450, 240, "3-2 flip", size=P.SIZE_LABEL, weight="700", fill=P.CABI)
    c.legend(
        40,
        372,
        [
            (P.DATA, "shared face removed by 2-3"),
            (P.ERROR, "edge p–q created by 2-3, removed by 3-2"),
        ],
    )
    c.label(
        450,
        395,
        "accepted only when the union is convex (a signed-volume test, no in-sphere predicate) and the worst scaled Jacobian strictly improves",
    )
    return c.render()


# --------------------------------------------------------------------------- #
# subdivide: one polyhedral child per face                                     #
# --------------------------------------------------------------------------- #
def subdivide_cell():
    c = Canvas(
        820,
        430,
        "Polyhedral refinement: every face of the cell becomes one child, closed by triangles to the new interior point",
    )
    t = tables()
    corners = [tuple(float(v) for v in p) for p in t.node_coordinates("hexahedron")]
    centre = centroid(corners)
    pts = corners + [centre]
    faces = t.corner_faces("hexahedron")
    scale = 150
    # left: the parent with its interior point and the x-max face's child edges
    box = (30, 70, 360, 280)
    origin = fit_origin(corners, scale, box)
    xmax = [i for i, f in enumerate(faces) if set(f) == {1, 2, 6, 5}][0]
    p2 = draw_solid(
        c,
        pts,
        faces,
        scale,
        origin,
        fill=P.CORE,
        fill_opacity=0.07,
        face_fills={xmax: (P.WASM, 0.25)},
    )
    for j in faces[xmax]:
        c.line(*p2[8], *p2[j], stroke=P.WASM, sw=1.6, dash="4 3")
    for i in range(8):
        c.circle(*p2[i], 3.5, fill=P.INK)
    c.circle(*p2[8], 5, fill=P.WHITE, stroke=P.CABI, sw=2)
    c.label(210, 372, "interior point = corner average of the cell", fill=P.CABI)
    c.caption(210, 60, "parent cell: 6 faces, one new interior point")
    # right: the child alone, exploded
    # orient: the x-max face keeps its winding; triangles wound outward for a
    # convex pyramid over it (apex on the inside of the parent).
    ring = faces[xmax]
    tri = [(ring[(k + 1) % 4], ring[k], 8) for k in range(4)]
    child = [tuple(ring)] + tri
    sub = [pts[i] for i in ring] + [centre]
    remap = {ring[k]: k for k in range(4)}
    remap[8] = 4
    child_local = [tuple(remap[i] for i in f) for f in child]
    box2 = (450, 90, 330, 240)
    origin2 = fit_origin(sub, scale, box2)
    p3 = draw_solid(
        c,
        sub,
        child_local,
        scale,
        origin2,
        fill=P.WASM,
        fill_opacity=0.12,
        face_fills={0: (P.WASM, 0.30)},
    )
    for i in range(4):
        c.circle(*p3[i], 3.5, fill=P.INK)
    c.circle(*p3[4], 5, fill=P.WHITE, stroke=P.CABI, sw=2)
    c.caption(615, 60, "one child: that face unchanged + one triangle per face edge")
    c.arrow(400, 210, 440, 210, stroke=P.INK_2)
    c.label(
        410,
        395,
        "the face keeps its original winding, so the neighbour across it still sees the identical face: the result conforms with no closure and no hanging nodes",
        anchor="middle",
    )
    c.label(
        410,
        411,
        "one polyhedron block per input block, mixed shapes inside; a wedge yields 2 tetra-shaped and 3 pyramid-shaped children",
        anchor="middle",
    )
    return c.render()


# --------------------------------------------------------------------------- #
# agglomerate: greedy seed-and-grow over the face dual                         #
# --------------------------------------------------------------------------- #
def agglomerate_grow():
    c = Canvas(
        900,
        440,
        "Polyhedral coarsening: greedy seed-and-grow over shared faces, largest accumulated shared area first",
    )
    widths = [1.0, 1.6, 0.8, 1.3]
    heights = [1.2, 0.7, 1.4]
    nx, ny = len(widths), len(heights)
    xs = [0.0]
    for w in widths:
        xs.append(xs[-1] + w)
    ys = [0.0]
    for h in heights:
        ys.append(ys[-1] + h)
    ncell = nx * ny

    def cid(i, j):
        return j * nx + i

    # shared-face "areas" are edge lengths in this 2-D schematic
    adj = {k: {} for k in range(ncell)}
    for j in range(ny):
        for i in range(nx):
            k = cid(i, j)
            if i + 1 < nx:
                adj[k][cid(i + 1, j)] = heights[j]
                adj[cid(i + 1, j)][k] = heights[j]
            if j + 1 < ny:
                adj[k][cid(i, j + 1)] = widths[i]
                adj[cid(i, j + 1)][k] = widths[i]
    target = 4
    group = [-1] * ncell
    steps = {}
    g = 0
    for seed in range(ncell):
        if group[seed] >= 0:
            continue
        members = [seed]
        group[seed] = g
        steps[seed] = (g, 0)
        while len(members) < target:
            frontier = {}
            for m in members:
                for n, area in adj[m].items():
                    if group[n] < 0:
                        frontier[n] = frontier.get(n, 0.0) + area
            if not frontier:
                break
            pick = sorted(frontier.items(), key=lambda kv: (-round(kv[1], 9), kv[0]))[
                0
            ][0]
            group[pick] = g
            steps[pick] = (g, len(members))
            members.append(pick)
        g += 1
    colours = [P.CORE, P.PYTHON, P.FORMATS, P.CABI, P.REGIONS]
    scale = 88
    ox, oy = 40, 60 + scale * sum(heights)

    def rect(i, j):
        return (
            ox + scale * xs[i],
            oy - scale * ys[j + 1],
            scale * widths[i],
            scale * heights[j],
        )

    for j in range(ny):
        for i in range(nx):
            k = cid(i, j)
            x, y, w, h = rect(i, j)
            col = colours[group[k] % len(colours)]
            c.rrect(x, y, w, h, fill=col, fill_opacity=0.16, stroke="none", rx=0)
    # dropped internal faces (dashed) and kept faces (solid)
    for j in range(ny):
        for i in range(nx):
            k = cid(i, j)
            x, y, w, h = rect(i, j)
            if i + 1 < nx:
                same = group[k] == group[cid(i + 1, j)]
                c.line(
                    x + w,
                    y,
                    x + w,
                    y + h,
                    stroke=P.MUTED if same else P.INK,
                    sw=1 if same else 2.2,
                    dash="4 3" if same else None,
                )
            if j + 1 < ny:
                same = group[k] == group[cid(i, j + 1)]
                c.line(
                    x,
                    y,
                    x + w,
                    y,
                    stroke=P.MUTED if same else P.INK,
                    sw=1 if same else 2.2,
                    dash="4 3" if same else None,
                )
    x0, y0, w0, h0 = rect(0, 0)
    c.polygon(
        [
            (ox, oy),
            (ox + scale * xs[-1], oy),
            (ox + scale * xs[-1], oy - scale * ys[-1]),
            (ox, oy - scale * ys[-1]),
        ],
        fill="none",
        stroke=P.INK,
        sw=2.2,
    )
    for j in range(ny):
        for i in range(nx):
            k = cid(i, j)
            x, y, w, h = rect(i, j)
            gk, order = steps[k]
            c.text(
                x + w / 2,
                y + h / 2 - 2,
                str(k),
                size=P.SIZE_LABEL,
                weight="700",
                mono=True,
            )
            c.label(
                x + w / 2,
                y + h / 2 + 13,
                "seed" if order == 0 else f"step {order}",
                fill=colours[gk % len(colours)],
            )
    # arrows for group 0's growth
    members0 = sorted(
        (k for k in range(ncell) if group[k] == 0), key=lambda k: steps[k][1]
    )
    for a, b in zip(members0, members0[1:]):
        ia, ja = a % nx, a // nx
        ib, jb = b % nx, b // nx
        xa, ya, wa, ha = rect(ia, ja)
        xb, yb, wb, hb = rect(ib, jb)
        c.arrow(
            xa + wa / 2,
            ya + ha / 2 + 22,
            xb + wb / 2,
            yb + hb / 2 + 22,
            stroke=colours[0],
            sw=1.6,
        )
    tx = ox + scale * xs[-1] + 40
    c.text(
        tx,
        70,
        "the rule, on this grid",
        size=P.SIZE_LABEL,
        weight="600",
        anchor="start",
    )
    lines = [
        "seeds are the lowest unclaimed cell ids, in order",
        "a group absorbs the unclaimed face-neighbour whose",
        "shared-face area, summed over every member, is largest",
        "(ties break on ascending id), until it holds 4 cells",
        "or its frontier is empty",
        "",
        "faces between two members are dropped (dashed);",
        "every other face is kept with its stored winding,",
        "so each group emits one polyhedron",
        "",
        "adjacency is by shared FACE, never by shared node:",
        "two cells touching at a corner never merge",
    ]
    for k, s in enumerate(lines):
        c.label(tx, 92 + 17 * k, s, anchor="start")
    c.label(
        450,
        425,
        "cell widths and heights differ on purpose: on a uniform grid every shared face has the same area and only the id tie-break is visible",
    )
    return c.render()


# --------------------------------------------------------------------------- #
# marching tetrahedra cases                                                     #
# --------------------------------------------------------------------------- #
def marching_cases():
    c = Canvas(
        960,
        366,
        "Marching tetrahedra: the sign mask of a simplex's nodes selects a triangle, a quad, a segment or nothing",
    )
    t = tables()
    tet = [tuple(float(v) for v in p) for p in t.node_coordinates("tetra")]
    tri = [tuple(float(v) for v in p) for p in t.node_coordinates("triangle")]
    scale = 125

    def crossing(p, q, dp, dq):
        tt = dp / (dp - dq)
        return tuple(p[k] + tt * (q[k] - p[k]) for k in range(3))

    def panel(box, pts, positive, edges, faces, title, flat):
        d = [1.0 if i in positive else -1.0 for i in range(len(pts))]
        origin = fit_origin(pts, scale, box, flat)
        if faces:
            p2 = draw_solid(
                c, pts, faces, scale, origin, fill=P.CORE, fill_opacity=0.06, flat=flat
            )
        else:
            p2 = project_all(pts, scale, origin, flat)
            c.polygon(p2, fill=P.CORE, fill_opacity=0.06, stroke=P.INK, sw=1.6)
        cuts = []
        for a, b in edges:
            if d[a] * d[b] < 0:
                cuts.append(crossing(pts[a], pts[b], d[a], d[b]))
        c2 = project_all(cuts, scale, origin, flat)
        if len(c2) >= 3:
            # order the ring by angle around its centroid (a convex section)
            cx = sum(x for x, _ in c2) / len(c2)
            cy = sum(y for _, y in c2) / len(c2)
            import math

            c2 = sorted(c2, key=lambda p: math.atan2(p[1] - cy, p[0] - cx))
            c.polygon(c2, fill=P.DATA, fill_opacity=0.45, stroke=P.DATA, sw=1.6)
        elif len(c2) == 2:
            c.line(*c2[0], *c2[1], stroke=P.DATA, sw=3)
        for x, y in c2:
            c.circle(x, y, 3.5, fill=P.WHITE, stroke=P.DATA, sw=2)
        cx = sum(x for x, _ in p2) / len(p2)
        cy = sum(y for _, y in p2) / len(p2)
        for i, (x, y) in enumerate(p2):
            dx, dy = x - cx, y - cy
            n = (dx * dx + dy * dy) ** 0.5 or 1.0
            if d[i] > 0:
                c.circle(x, y, 4.5, fill=P.CORE)
            else:
                c.circle(x, y, 4.5, fill=P.WHITE, stroke=P.INK, sw=1.8)
            c.text(
                x + 13 * dx / n,
                y + 13 * dy / n + 4,
                "+" if d[i] > 0 else "−",
                size=P.SIZE_SMALL,
                weight="700",
                fill=P.CORE if d[i] > 0 else P.INK_2,
            )
        c.caption(box[0] + box[2] / 2, box[1] - 14, title)

    tet_edges = t.corner_edges("tetra")
    tet_faces = t.corner_faces("tetra")
    tri_edges = t.corner_edges("triangle")
    panel(
        (30, 70, 200, 200),
        tet,
        {3},
        tet_edges,
        tet_faces,
        "tetra, 1 positive → triangle",
        False,
    )
    panel(
        (260, 70, 200, 200),
        tet,
        {2, 3},
        tet_edges,
        tet_faces,
        "tetra, 2 positive → quad",
        False,
    )
    panel(
        (490, 70, 200, 200),
        tri,
        {2},
        tri_edges,
        [],
        "triangle, 1 positive → segment",
        True,
    )
    panel(
        (720, 70, 200, 200),
        tet,
        set(),
        tet_edges,
        tet_faces,
        "all one sign → nothing",
        False,
    )
    c.legend(
        30,
        310,
        [
            (P.CORE, "positive node (d ≥ 0; a node ON the plane counts positive)"),
            (P.DATA, "section through the crossings"),
        ],
    )
    c.label(
        480,
        334,
        "each crossing sits at t = d_lo / (d_lo − d_hi) on the sorted endpoint pair, so two simplices sharing an edge place it bit for bit identically;",
    )
    c.label(
        480,
        350,
        "the sign mask is total, so a plane grazing a shared face is emitted by exactly one of the two simplices",
    )
    return c.render()


# --------------------------------------------------------------------------- #
# remesh_volume: the BCC lattice and isosurface stuffing                       #
# --------------------------------------------------------------------------- #
def remesh_volume_bcc():
    c = Canvas(
        960,
        470,
        "Isosurface stuffing: a body-centred cubic lattice, classified by signed distance, warped onto the surface and cut",
    )
    import math

    # panel 1: one BCC cell with the four tets over its +x face
    corners = [
        (0, 0, 0),
        (1, 0, 0),
        (1, 1, 0),
        (0, 1, 0),
        (0, 0, 1),
        (1, 0, 1),
        (1, 1, 1),
        (0, 1, 1),
    ]
    corners = [tuple(float(v) for v in p) for p in corners]
    cc = (0.5, 0.5, 0.5)
    cn = (1.5, 0.5, 0.5)
    pts = corners + [cc, cn]
    t = tables()
    faces = t.corner_faces("hexahedron")
    scale = 130
    box = (30, 80, 300, 260)
    origin = fit_origin(pts, scale, box)
    xmax = [i for i, f in enumerate(faces) if set(f) == {1, 2, 6, 5}][0]
    p2 = draw_solid(c, pts, faces, scale, origin, fill=P.CORE, fill_opacity=0.05)
    ring = faces[xmax]
    for k in range(4):
        a, b = ring[k], ring[(k + 1) % 4]
        col = P.WASM if k == 0 else P.CORE
        c.polygon(
            [p2[a], p2[b], p2[9]],
            fill=col,
            fill_opacity=0.28 if k == 0 else 0.10,
            stroke=col,
            sw=1.2,
        )
        c.polygon(
            [p2[a], p2[b], p2[8]],
            fill=col,
            fill_opacity=0.28 if k == 0 else 0.10,
            stroke=col,
            sw=1.2,
            dash="3 3",
        )
    c.line(*p2[8], *p2[9], stroke=P.CABI, sw=2)
    for i in range(8):
        c.circle(*p2[i], 3.5, fill=P.INK)
    c.circle(*p2[8], 5, fill=P.WHITE, stroke=P.CABI, sw=2)
    c.circle(*p2[9], 5, fill=P.WHITE, stroke=P.CABI, sw=2)
    c.label(p2[8][0] - 6, p2[8][1] - 10, "centre", anchor="end", fill=P.CABI)
    c.label(
        p2[9][0] + 10, p2[9][1] + 4, "neighbour's centre", anchor="start", fill=P.CABI
    )
    c.caption(180, 62, "BCC lattice cell")
    c.label(180, 358, "two body centres across each face plus one face edge")
    c.label(180, 374, "make a tet: 4 per face, 12 per cell; every dihedral")
    c.label(180, 390, "angle is one of 45°, 60°, 90° or 120°")

    # panels 2-4: a 2-D cut through classify -> warp -> cut
    def lattice_panel(x0, title, subtitle, mode):
        h = 44.0
        cx, cy, r = x0 + 90, 250, 78
        c.caption(x0 + 90, 62, title)
        c.label(x0 + 90, 80, subtitle)
        # the surface
        c.circle(cx, cy, r, fill="none", stroke=P.FORMATS, sw=2)
        nodes = []
        for i in range(5):
            for j in range(5):
                nodes.append((x0 + 2 + i * h, 162 + j * h))
        moved = {}
        for k, (x, y) in enumerate(nodes):
            d = math.hypot(x - cx, y - cy) - r
            if mode in ("warp", "cut") and abs(d) < 0.3 * h:
                ux, uy = (x - cx) / (r + d), (y - cy) / (r + d)
                moved[k] = (cx + r * ux, cy + r * uy)
                if mode == "warp":
                    c.arrow(x, y, moved[k][0], moved[k][1], stroke=P.PYTHON, sw=1.4)
        # the tets are triangles here: draw a few straddling ones cut
        if mode == "cut":
            for k, (x, y) in enumerate(nodes):
                i, j = divmod(k, 5)
                if i < 4 and j < 4:
                    quad = [k, k + 5, k + 6, k + 1]
                    for tri in (
                        (quad[0], quad[1], quad[2]),
                        (quad[0], quad[2], quad[3]),
                    ):
                        P3 = [moved.get(m, nodes[m]) for m in tri]
                        ins = [
                            math.hypot(px - cx, py - cy) - r <= 1e-6 for px, py in P3
                        ]
                        if any(ins):
                            c.polygon(
                                P3,
                                fill=P.CORE,
                                fill_opacity=0.16 if all(ins) else 0.30,
                                stroke=P.CORE,
                                sw=0.8,
                            )
        for k, (x, y) in enumerate(nodes):
            px, py = moved.get(k, (x, y)) if mode != "classify" else (x, y)
            d = math.hypot(px - cx, py - cy) - r
            if d <= 1e-6:
                c.circle(px, py, 3.5, fill=P.CORE)
            else:
                c.circle(px, py, 3.5, fill=P.WHITE, stroke=P.INK_2, sw=1.4)

    lattice_panel(360, "1. classify", "by the sign of the distance", "classify")
    lattice_panel(560, "2. warp", "vertices within warp_fraction · h", "warp")
    lattice_panel(760, "3. cut", "the tets that straddle the surface", "cut")
    c.legend(
        360,
        400,
        [
            (P.CORE, "inside vertex / kept tet"),
            (P.FORMATS, "the surface (signed distance zero)"),
            (P.PYTHON, "warp onto the surface"),
        ],
    )
    c.label(
        480,
        428,
        "warping is load-bearing: without it the cut tets near the surface get arbitrarily bad dihedral angles, plain marching tetrahedra's usual failure mode;",
    )
    c.label(
        480,
        446,
        "warp_fraction = 0 stays exactly watertight, and a nonzero value trades a few non-manifold boundary edges (reported, never hidden) for element quality",
    )
    return c.render()


FIGURES = {
    "extract_skin_hashing": extract_skin_hashing,
    "flips_23_32": flips_23_32,
    "subdivide_cell": subdivide_cell,
    "agglomerate_grow": agglomerate_grow,
    "marching_cases": marching_cases,
    "remesh_volume_bcc": remesh_volume_bcc,
}

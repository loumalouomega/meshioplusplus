"""Drawing a reference element: wireframe, hidden edges, numbered nodes."""

from . import palette as P
from .project import depth, face_is_visible, project
from .tables import LINEAR_BASE, REFERENCE_CORNERS


def project_all(points, scale, origin, flat=False):
    if flat:
        return [(origin[0] + scale * p[0], origin[1] - scale * p[1]) for p in points]
    return [project(p, scale, origin) for p in points]


def fit_origin(points, scale, box, flat=False):
    """The origin that centres ``points`` (3-D) inside the SVG rect ``box``."""
    pts = project_all(points, scale, (0.0, 0.0), flat)
    xs = [x for x, _ in pts]
    ys = [y for _, y in pts]
    bx, by, bw, bh = box
    cx = bx + bw / 2 - (min(xs) + max(xs)) / 2
    cy = by + bh / 2 - (min(ys) + max(ys)) / 2
    return (cx, cy)


def draw_cell(
    canvas,
    tables,
    cell_type,
    origin,
    scale,
    labels=True,
    fill=P.CORE,
    fill_opacity=0.08,
    edge_color=P.INK,
    nodes=None,
    node_r=4.5,
    label_offsets=None,
    flat=False,
    box=None,
):
    """Draw ``cell_type`` at ``origin`` (SVG px) with ``scale`` px per unit.

    ``nodes`` overrides the node coordinates (defaults to the reference
    element from the tables); ``label_offsets`` maps a node index to a
    ``(dx, dy)`` label nudge for the few labels that would otherwise collide;
    ``flat`` draws a 2-D element in the plane (x right, y up) instead of
    through the 3-D camera; ``box`` (an SVG rect) centres the drawing there and
    overrides ``origin``. Returns the projected node positions.
    """
    pts3 = nodes if nodes is not None else tables.node_coordinates(cell_type)
    if box is not None:
        origin = fit_origin(pts3, scale, box, flat)
    pts = project_all(pts3, scale, origin, flat)
    base = LINEAR_BASE[cell_type]
    nc = len(REFERENCE_CORNERS[base])
    faces = tables.corner_faces(cell_type)
    edges = tables.corner_edges(cell_type)
    if faces:
        visible = [face_is_visible([pts3[i] for i in f]) for f in faces]
        vis_edges = set()
        for f, v in zip(faces, visible):
            if not v:
                continue
            for k in range(len(f)):
                vis_edges.add(frozenset((f[k], f[(k + 1) % len(f)])))
        # hidden edges first, dashed
        for a, b in edges:
            if frozenset((a, b)) not in vis_edges:
                canvas.line(*pts[a], *pts[b], stroke=P.MUTED, sw=1.2, dash="4 3")
        # visible faces, far to near
        order = sorted(
            range(len(faces)), key=lambda i: depth([pts3[j] for j in faces[i]][0])
        )
        for i in order:
            if visible[i]:
                canvas.polygon(
                    [pts[j] for j in faces[i]],
                    fill=fill,
                    fill_opacity=fill_opacity,
                    stroke="none",
                )
        for a, b in edges:
            if frozenset((a, b)) in vis_edges:
                canvas.line(*pts[a], *pts[b], stroke=edge_color, sw=1.6)
    else:
        if len(edges) >= 3:
            canvas.polygon(
                [pts[i] for i in range(nc)],
                fill=fill,
                fill_opacity=fill_opacity,
                stroke="none",
            )
        for a, b in edges:
            canvas.line(*pts[a], *pts[b], stroke=edge_color, sw=1.6)
    if nodes is None or labels:
        cx = sum(x for x, _ in pts[:nc]) / nc
        cy = sum(y for _, y in pts[:nc]) / nc
        for i, (x, y) in enumerate(pts):
            kind = tables.node_kind(cell_type, i)
            dx, dy = (label_offsets or {}).get(i, radial_offset(x, y, cx, cy))
            canvas.node(
                x,
                y,
                i,
                kind=kind,
                r=node_r,
                label_dx=dx,
                label_dy=dy,
                label=labels,
                anchor="middle",
            )
    return pts


def radial_offset(x, y, cx, cy, dist=13.0):
    """Push a label away from the element centre, so it clears the edges."""
    dx, dy = x - cx, y - cy
    n = (dx * dx + dy * dy) ** 0.5
    if n < 1e-9:
        return (0.0, -dist)
    return (dist * dx / n, dist * dy / n + 4.0)


def solid_faces_edges(faces):
    """The unique undirected edges of a face list."""
    edges = set()
    for f in faces:
        for k in range(len(f)):
            edges.add(frozenset((f[k], f[(k + 1) % len(f)])))
    return sorted(tuple(sorted(e)) for e in edges)


def draw_solid(
    canvas,
    pts3,
    faces,
    scale,
    origin,
    fill=P.CORE,
    fill_opacity=0.10,
    edge_color=P.INK,
    sw=1.6,
    hidden=True,
    face_fills=None,
    box=None,
    flat=False,
):
    """Draw an outward-wound closed solid: hidden edges dashed, visible faces filled.

    ``face_fills`` maps a face index to ``(colour, opacity)`` to highlight it.
    Returns the projected points.
    """
    if box is not None:
        origin = fit_origin(pts3, scale, box, flat)
    pts = project_all(pts3, scale, origin, flat)
    visible = [face_is_visible([pts3[i] for i in f]) for f in faces]
    vis_edges = set()
    for f, v in zip(faces, visible):
        if v:
            for k in range(len(f)):
                vis_edges.add(frozenset((f[k], f[(k + 1) % len(f)])))
    edges = solid_faces_edges(faces)
    if hidden:
        for a, b in edges:
            if frozenset((a, b)) not in vis_edges:
                canvas.line(*pts[a], *pts[b], stroke=P.MUTED, sw=1.2, dash="4 3")
    order = sorted(
        range(len(faces)), key=lambda i: depth([pts3[j] for j in faces[i]][0])
    )
    for i in order:
        if visible[i]:
            colour, op = (face_fills or {}).get(i, (fill, fill_opacity))
            canvas.polygon(
                [pts[j] for j in faces[i]], fill=colour, fill_opacity=op, stroke="none"
            )
    for a, b in edges:
        if frozenset((a, b)) in vis_edges:
            canvas.line(*pts[a], *pts[b], stroke=edge_color, sw=sw)
    return pts


def oriented_tet(p, a, b, c, d):
    """Outward-wound faces of tetra ``(a, b, c, d)`` over point list ``p``."""
    from .project import cross, dot, sub

    v = dot(cross(sub(p[b], p[a]), sub(p[c], p[a])), sub(p[d], p[a]))
    if v < 0:
        b, c = c, b
    return [(a, c, b), (a, b, d), (b, c, d), (c, a, d)]

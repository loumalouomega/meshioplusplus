"""The numpy twin of ``formats/step.cpp``: a STEP AP214 faceted B-rep.

Each surface facet is an ``ADVANCED_FACE`` on a ``PLANE`` bounded by an
``EDGE_LOOP`` over shared straight ``EDGE_CURVE``s; each edge-connected piece
is one ``CLOSED_SHELL`` in one ``MANIFOLD_SOLID_BREP``. The bytes equal the
core's.
"""

from __future__ import annotations

import math

import numpy as np

from .._exceptions import WriteError
from .._files import is_buffer

_PREFIX = "meshio++: step: "


def _real(v):
    s = "%.17g" % v
    mant, _, expo = s.partition("e")
    if "." not in mant:
        mant += "."
    return mant + ("E" + expo if expo else "")


def _string(text):
    out = []
    for c in text.encode("utf-8"):
        if c == 0x27:
            out.append("''")
        elif 0x20 <= c < 0x7F:
            out.append(chr(c))
        else:
            out.append("_")
    return "".join(out)


def _basename(path):
    path = str(path)
    cut = max(path.rfind("/"), path.rfind("\\"))
    return path if cut < 0 else path[cut + 1 :]


def _unit(v):
    ln = math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])
    return [x / ln if ln > 0.0 else 0.0 for x in v]


def _ring_normal(xyz, ring):
    from .._feature_edges import ring_unit_normals

    start = np.array([0, len(ring)], dtype=np.int64)
    return ring_unit_normals(xyz, start, np.asarray(ring, dtype=np.int64))[0].tolist()


def _facets(mesh):
    """``(xyz, rings)`` of the surface to write, after skinning and fanning."""
    from .._common import warn
    from .._feature_edges import _coords3
    from .._segment import surface_rings
    from .._skin import _has_skinnable_cells, extract_skin

    surface = extract_skin(mesh, linearize=True) if _has_skinnable_cells(mesh) else mesh
    total = sum(len(b.data) for b in surface.cells)
    try:
        start, nodes, cells = surface_rings(
            surface, np.ones(total, dtype=bool), prefix=_PREFIX
        )
    except ValueError as exc:
        raise WriteError(str(exc)) from None
    if len(cells) < total:
        warn(
            f"step: {total - len(cells)} cell(s) that are not surface facets "
            "(lines, vertices) are dropped"
        )
    if len(cells) == 0:
        raise WriteError(f"{_PREFIX}no surface facets to write")
    xyz = _coords3(surface)
    if not np.all(np.isfinite(xyz)):
        raise WriteError(f"{_PREFIX}a coordinate is not finite")
    pts = xyz.tolist()
    rings = []
    for f in range(len(cells)):
        ring = nodes[start[f] : start[f + 1]].tolist()
        planar = True
        if len(ring) > 3:
            n = _ring_normal(xyz, ring)
            p0 = pts[ring[0]]
            dev = extent = 0.0
            for k in ring[1:]:
                p = pts[k]
                d = [p[0] - p0[0], p[1] - p0[1], p[2] - p0[2]]
                h = abs(d[0] * n[0] + d[1] * n[1] + d[2] * n[2])
                ln = math.sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2])
                dev = h if h > dev else dev
                extent = ln if ln > extent else extent
            planar = not dev > 1e-9 * extent
        if planar:
            rings.append(ring)
        else:
            rings.extend(
                [ring[0], ring[k], ring[k + 1]] for k in range(1, len(ring) - 1)
            )
    return xyz, pts, rings


def _check_closed(xyz, rings):
    from .._feature_edges import crease_edges, ring_unit_normals

    start = np.zeros(len(rings) + 1, dtype=np.int64)
    start[1:] = np.cumsum([len(r) for r in rings])
    nodes = np.concatenate([np.asarray(r, dtype=np.int64) for r in rings])
    normals = ring_unit_normals(xyz, start, nodes)
    for f in range(len(rings)):
        if not np.any(normals[f] != 0.0):
            raise WriteError(f"{_PREFIX}facet {f} is degenerate (zero area)")
    e = crease_edges(start, nodes, normals, 180.0)
    boundary = e["uses"] == 1
    nonmanifold = e["uses"] >= 3
    bad = boundary | nonmanifold | e["inconsistent"]
    if np.any(bad):
        idx = np.nonzero(bad)[0][:5]
        first = ", ".join(f"({int(e['lo'][i])}, {int(e['hi'][i])})" for i in idx)
        raise WriteError(
            f"{_PREFIX}the surface is not a closed, consistently wound solid: "
            f"{int(boundary.sum())} open, {int(nonmanifold.sum())} non-manifold and "
            f"{int(e['inconsistent'].sum())} inconsistently wound edge(s), first "
            f"{first} (point ids); run repair first"
        )


def _pieces(rings):
    """Per facet, its edge-connected piece (pieces ordered by lowest facet)."""
    parent = list(range(len(rings)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    owner = {}
    for f, ring in enumerate(rings):
        s = len(ring)
        for k in range(s):
            u, v = ring[k], ring[(k + 1) % s]
            key = (min(u, v), max(u, v))
            if key not in owner:
                owner[key] = f
                continue
            ra, rb = find(f), find(owner[key])
            if ra != rb:
                parent[max(ra, rb)] = min(ra, rb)
    index = {}
    piece = []
    for f in range(len(rings)):
        r = find(f)
        if r == f:
            index[f] = len(index)
        piece.append(index[r])
    return piece, len(index)


def _render(name, pts, rings):
    from .. import _provenance
    from .._common import warn

    piece_of, npieces = _pieces(rings)
    volume = [0.0] * npieces
    for f, ring in enumerate(rings):
        p0 = pts[ring[0]]
        for k in range(1, len(ring) - 1):
            p1, p2 = pts[ring[k]], pts[ring[k + 1]]
            volume[piece_of[f]] += (
                p0[0] * (p1[1] * p2[2] - p1[2] * p2[1])
                - p0[1] * (p1[0] * p2[2] - p1[2] * p2[0])
                + p0[2] * (p1[0] * p2[1] - p1[1] * p2[0])
            ) / 6.0
    flipped = sum(1 for v in volume if v < 0.0)
    if flipped:
        warn(
            f"step: {flipped} piece(s) of the surface are wound inwards; written reversed"
        )

    edges = sorted(
        {
            (min(r[k], r[(k + 1) % len(r)]), max(r[k], r[(k + 1) % len(r)]))
            for r in rings
            for k in range(len(r))
        }
    )
    edge_index = {e: i for i, e in enumerate(edges)}
    used = sorted({v for r in rings for v in r})
    used_index = {v: i for i, v in enumerate(used)}

    out = []
    counter = [0]

    def ent(body):
        counter[0] += 1
        out.append(f"#{counter[0]}={body};\n")
        return counter[0]

    def triple(v):
        return f"({_real(v[0])},{_real(v[1])},{_real(v[2])})"

    head = ["ISO-10303-21;\n/*\n"]
    for line in _provenance.lines(_provenance.SlotTier.BLOCK):
        head.append(line.replace("*/", "* /") + "\n")
    head.append("*/\nHEADER;\n")
    head.append("FILE_DESCRIPTION(('meshio++ faceted B-rep'),'2;1');\n")
    head.append(
        f"FILE_NAME('{name}','1970-01-01T00:00:00',(''),(''),'meshio++','meshio++','');\n"
    )
    head.append(
        "FILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }'));\nENDSEC;\nDATA;\n"
    )

    app = ent(
        "APPLICATION_CONTEXT('core data for automotive mechanical design processes')"
    )
    ent(
        "APPLICATION_PROTOCOL_DEFINITION('international standard','automotive_design',"
        f"2000,#{app})"
    )
    pctx = ent(f"PRODUCT_CONTEXT('',#{app},'mechanical')")
    pdctx = ent(f"PRODUCT_DEFINITION_CONTEXT('part definition',#{app},'design')")
    prod = ent(f"PRODUCT('{name}','{name}','',(#{pctx}))")
    ent(f"PRODUCT_RELATED_PRODUCT_CATEGORY('part',$,(#{prod}))")
    form = ent(f"PRODUCT_DEFINITION_FORMATION('','',#{prod})")
    pdef = ent(f"PRODUCT_DEFINITION('design','',#{form},#{pdctx})")
    pshape = ent(f"PRODUCT_DEFINITION_SHAPE('','',#{pdef})")
    len_unit = ent("(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.))")
    ang_unit = ent("(NAMED_UNIT(*)PLANE_ANGLE_UNIT()SI_UNIT($,.RADIAN.))")
    sol_unit = ent("(NAMED_UNIT(*)SI_UNIT($,.STERADIAN.)SOLID_ANGLE_UNIT())")
    unc = ent(
        f"UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-07),#{len_unit},"
        "'distance_accuracy_value','confusion accuracy')"
    )
    gctx = ent(
        f"(GEOMETRIC_REPRESENTATION_CONTEXT(3)GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#{unc}))"
        f"GLOBAL_UNIT_ASSIGNED_CONTEXT((#{len_unit},#{ang_unit},#{sol_unit}))"
        "REPRESENTATION_CONTEXT('Context #1','3D Context with UNIT and UNCERTAINTY'))"
    )
    o_pt = ent("CARTESIAN_POINT('',(0.,0.,0.))")
    o_z = ent("DIRECTION('',(0.,0.,1.))")
    o_x = ent("DIRECTION('',(1.,0.,0.))")
    origin = ent(f"AXIS2_PLACEMENT_3D('',#{o_pt},#{o_z},#{o_x})")

    point_id, vertex_id = [], []
    for v in used:
        point_id.append(ent(f"CARTESIAN_POINT('',{triple(pts[v])})"))
        vertex_id.append(ent(f"VERTEX_POINT('',#{point_id[-1]})"))
    edge_id = []
    for lo, hi in edges:
        a, b = pts[lo], pts[hi]
        d = [b[0] - a[0], b[1] - a[1], b[2] - a[2]]
        ln = math.sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2])
        if not ln > 0.0:
            raise WriteError(
                f"{_PREFIX}edge ({lo}, {hi}) has zero length "
                "(coincident points; run clean first)"
            )
        ia, ib = used_index[lo], used_index[hi]
        direction = ent(f"DIRECTION('',{triple(_unit(d))})")
        vec = ent(f"VECTOR('',#{direction},{_real(ln)})")
        line = ent(f"LINE('',#{point_id[ia]},#{vec})")
        edge_id.append(
            ent(f"EDGE_CURVE('',#{vertex_id[ia]},#{vertex_id[ib]},#{line},.T.)")
        )

    xyz = np.asarray(pts, dtype=np.float64)
    face_id = []
    for f, ring in enumerate(rings):
        ring = list(reversed(ring)) if volume[piece_of[f]] < 0.0 else list(ring)
        s = len(ring)
        loop = []
        for k in range(s):
            u, v = ring[k], ring[(k + 1) % s]
            e = edge_id[edge_index[(min(u, v), max(u, v))]]
            sense = ".T." if u < v else ".F."
            oriented = ent(f"ORIENTED_EDGE('',*,*,#{e},{sense})")
            loop.append(f"#{oriented}")
        el = ent(f"EDGE_LOOP('',({','.join(loop)}))")
        bound = ent(f"FACE_OUTER_BOUND('',#{el},.T.)")
        n = _ring_normal(xyz, ring)
        p0, p1 = pts[ring[0]], pts[ring[1]]
        d = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]]
        h = d[0] * n[0] + d[1] * n[1] + d[2] * n[2]
        x = _unit([d[0] - h * n[0], d[1] - h * n[1], d[2] - h * n[2]])
        nz = ent(f"DIRECTION('',{triple(n)})")
        nx = ent(f"DIRECTION('',{triple(x)})")
        place = ent(
            f"AXIS2_PLACEMENT_3D('',#{point_id[used_index[ring[0]]]},#{nz},#{nx})"
        )
        plane = ent(f"PLANE('',#{place})")
        face_id.append(ent(f"ADVANCED_FACE('',(#{bound}),#{plane},.T.)"))

    solids = ""
    for c in range(npieces):
        faces = ",".join(
            f"#{face_id[f]}" for f in range(len(rings)) if piece_of[f] == c
        )
        shell = ent(f"CLOSED_SHELL('',({faces}))")
        solid = ent(f"MANIFOLD_SOLID_BREP('solid{c + 1}',#{shell})")
        solids += f",#{solid}"
    rep = ent(f"ADVANCED_BREP_SHAPE_REPRESENTATION('',(#{origin}{solids}),#{gctx})")
    ent(f"SHAPE_DEFINITION_REPRESENTATION(#{pshape},#{rep})")
    out.append("ENDSEC;\nEND-ISO-10303-21;\n")
    return "".join(head) + "".join(out)


def write(filename, mesh):
    """Write a closed surface (or a volume mesh's skin) as a STEP faceted solid."""
    xyz, pts, rings = _facets(mesh)
    _check_closed(xyz, rings)
    buffer = is_buffer(filename, "w")
    name = _string(_basename(getattr(filename, "name", "") if buffer else filename))
    text = _render(name, pts, rings).encode("ascii")
    if buffer:
        try:
            filename.write(text)
        except TypeError:
            filename.write(text.decode("ascii"))
        return
    with open(filename, "wb") as f:
        f.write(text)

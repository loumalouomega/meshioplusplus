#!/usr/bin/env python3
"""Regenerate the Kratos CAD JSON fixtures in ``tests/python/meshes/cad_json/``.

The files follow the schema Kratos's ``CadJsonInput`` reads (Teschemacher et
al., AMSES 2018; ``kratos/input_output/cad_json_input.h``) and the
``Surfaces``/``Lines`` dialect of ``ImportNurbsSbmModeler``. They are built
here from closed-form geometry, not copied from Kratos, and written with
meshio++'s own writer, so they are canonical (a rewrite gives the same bytes):

- ``plate_with_hole.cad.json``: a 4 x 4 plate of degree (2, 2) trimmed by a
  square outer loop of lines and an inner loop of radius 1 made of four
  rational quadratic arcs (exact circles), each trim also a 3-D edge.
- ``cylindrical_roof.cad.json``: an untrimmed rational patch of degrees (2, 1)
  that is exactly an 80-degree arc of a cylinder of radius 25 (the
  Scordelis-Lo roof), its knot vectors in Kratos's reduced form.
- ``two_patch_beam.cad.json``: two trimmed patches (degrees (3, 2), two knot
  spans, and (2, 2)) joined by a coupling edge whose topology names both
  faces, a vertex on the first face, an embedded edge, a face named by
  ``brep_name``, trimming curves written without control-point ids, and
  members meshio++ does not model (kept verbatim).
- ``import_square.json``: the ``ImportNurbsSbmModeler`` dialect: a bilinear
  surface whose weights are all zero (that modeler reads no surface weights)
  and four boundary lines plus a rational quarter circle.
- ``unit_cube.cad.json``, ``capped_cylinder.cad.json`` and ``sphere.cad.json``:
  closed solids with shared edges, for ``mesh_cad_surface``: six bilinear
  faces and twelve edges; a rational cylinder band closed through a seam and
  two planar caps trimmed by exact circles that start where the band's rims
  do; two rational hemispheres sharing the equator, each closed through a
  seam and a degenerate pole edge with no 3-D curve. Every face's normal
  ``dS/du x dS/dv`` points outward.
- ``capped_cylinder_offset.cad.json``: the capped cylinder with its top cap's
  circle started a quarter turn after the band's rim, a closed edge whose
  uses start at different points (``mesh_cad_surface`` rotates the cap's).
- ``plate_with_near_hole.cad.json``: a disc of radius 2 with a hole of radius
  1 whose edge comes within 0.02 of the rim, the rim's circle started at 225
  degrees: at ``h = 1.1`` its coarse loops cross (``mesh_cad_surface`` splits
  them).
- ``split_rim.cad.json``: two unit plates folded along a shared rim that is
  split at x = 0.5 and 0.505 into three edges, the middle one a sliver (for
  ``mesh_cad_surface``'s ``virtual_edges``).
- ``sliver_strip.cad.json``: three plates in a row in the plane z = 0, each
  its own patch, the middle one 0.01 wide, a sliver face (for
  ``mesh_cad_surface``'s ``virtual_faces``).

    python tools/gen_cad_json_fixtures.py
"""

import json
import math
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src" / "python"))

import numpy as np  # noqa: E402

from meshioplusplus.cad_json._cad_json import render  # noqa: E402

OUT = ROOT / "tests" / "python" / "meshes" / "cad_json"
W = math.sqrt(0.5)


class Ids:
    def __init__(self, start=100):
        self.next = start

    def take(self, n):
        out = list(range(self.next, self.next + n))
        self.next += n
        return np.array(out, dtype=np.int64)


def curve(
    degree, knots, points, weights=None, active=None, ids=None, reduced=False, extra=""
):
    pts = np.array([list(p) + [0.0] * (3 - len(p)) for p in points], dtype=np.float64)
    w = np.ones(len(pts)) if weights is None else np.array(weights, dtype=np.float64)
    return {
        "degree": degree,
        "knots": np.array(knots, dtype=np.float64),
        "reduced_knots": reduced,
        "rational": weights is not None,
        "control_points": pts,
        "weights": w,
        "ids": np.zeros(0, dtype=np.int64) if ids is None else ids,
        "active_range": [] if active is None else list(active),
        "extra": extra,
    }


def line(a, b, ids):
    """A straight trim; ``ids=None`` writes its control points without ids."""
    return curve(
        1, [0, 0, 1, 1], [a, b], active=[0, 1], ids=None if ids is None else ids.take(2)
    )


def arc(centre, r, a0, ids):
    """A rational quadratic quarter circle from angle a0, counter-clockwise."""
    cx, cy = centre
    a1 = a0 + math.pi / 2
    am = a0 + math.pi / 4
    p0 = (cx + r * math.cos(a0), cy + r * math.sin(a0))
    p2 = (cx + r * math.cos(a1), cy + r * math.sin(a1))
    p1 = (cx + r / W * math.cos(am), cy + r / W * math.sin(am))
    return curve(
        2,
        [0, 0, 0, 1, 1, 1],
        [p0, p1, p2],
        weights=[1, W, 1],
        active=[0, 1],
        ids=ids.take(3),
    )


def trim(index, c, direction=True):
    return {"trim_index": index, "direction": direction, "curve": c, "extra": ""}


def ref(face_id, trim_index=-1, relative=1, name="", local=None):
    return {
        "id": face_id,
        "name": name,
        "trim_index": trim_index,
        "relative_direction": relative,
        "local_coordinates": [] if local is None else list(local),
        "extra": "",
    }


def patch(degrees, knots, points, weights=None):
    pts = np.array(points, dtype=np.float64)
    return {
        "degrees": np.array(degrees, dtype=np.int32),
        "knots": [np.array(k, dtype=np.float64) for k in knots],
        "control_points": pts,
        "weights": (
            np.ones(len(pts))
            if weights is None
            else np.array(weights, dtype=np.float64)
        ),
        "fields": [],
    }


def face(
    ident,
    patch_index,
    loops=(),
    rational=False,
    trimmed=True,
    name="",
    ids=None,
    reduced=False,
    embedded=(),
    extra="",
):
    return {
        "id": ident,
        "name": name,
        "patch": patch_index,
        "swapped_normal": False,
        "trimmed": trimmed,
        "rational": rational,
        "reduced_knots": reduced,
        "control_point_ids": np.zeros(0, dtype=np.int64) if ids is None else ids,
        "loops": list(loops),
        "embedded_edges": list(embedded),
        "surface_extra": "",
        "extra": extra,
    }


def edge(ident, c, topology, name=""):
    return {"id": ident, "name": name, "curve": c, "topology": topology, "extra": ""}


def lifted(c):
    out = dict(c)
    out["control_points"] = c["control_points"].copy()
    return out


def grid(xs, ys, z=0.0):
    return [(x, y, z) for y in ys for x in xs]


def plate_with_hole():
    ids = Ids()
    p = patch(
        [2, 2], [[0, 0, 0, 4, 4, 4], [0, 0, 0, 4, 4, 4]], grid([0, 2, 4], [0, 2, 4])
    )
    corners = [(0, 0), (4, 0), (4, 4), (0, 4)]
    outer = [trim(i, line(corners[i], corners[(i + 1) % 4], ids)) for i in range(4)]
    # The hole runs clockwise: each counter-clockwise arc read backwards.
    inner = [
        trim(4 + i, arc((2, 2), 1.0, i * math.pi / 2, ids), direction=False)
        for i in range(4)
    ]
    f = face(
        2,
        0,
        loops=[
            {"type": "outer", "trims": outer, "extra": ""},
            {"type": "inner", "trims": inner, "extra": ""},
        ],
        ids=ids.take(9),
    )
    edges = []
    for k, t in enumerate(outer + inner):
        edges.append(edge(3 + k, lifted(t["curve"]), [ref(2, t["trim_index"], 1)]))
    return {
        "patches": [p],
        "breps": [
            {
                "id": 1,
                "name": "",
                "faces": [0],
                "edges": list(range(len(edges))),
                "vertices": [],
                "extra": "",
            }
        ],
        "faces": [f],
        "edges": edges,
        "vertices": [],
        "extra": '{"tolerances":{"model_tolerance":0.001},"version_number":1}',
    }


def cylindrical_roof():
    r, half, length = 25.0, math.radians(40.0), 50.0
    ring = [
        (-r * math.sin(half), r * math.cos(half)),
        (0.0, r / math.cos(half)),
        (r * math.sin(half), r * math.cos(half)),
    ]
    pts = [(x, y, z) for y in (0.0, length) for x, z in ring]
    w = [1.0, math.cos(half), 1.0] * 2
    # Full knots here; the face's reduced_knots writes them without the ends.
    p = patch([2, 1], [[0, 0, 0, 1, 1, 1], [0, 0, length, length]], pts, w)
    f = face(
        2,
        0,
        rational=True,
        trimmed=False,
        reduced=True,
        ids=np.arange(1, 7, dtype=np.int64),
    )
    return {
        "patches": [p],
        "breps": [
            {
                "id": 1,
                "name": "",
                "faces": [0],
                "edges": [],
                "vertices": [],
                "extra": "",
            }
        ],
        "faces": [f],
        "edges": [],
        "vertices": [],
        "extra": '{"version_number":1}',
    }


def two_patch_beam():
    ids = Ids()
    # Patch A: x in [0, 5], degree 3 in x over two spans, degree 2 in y.
    xs = [0.0, 5 / 6, 2.5, 25 / 6, 5.0]
    a = patch(
        [3, 2],
        [[0, 0, 0, 0, 2.5, 5, 5, 5, 5], [0, 0, 0, 1, 1, 1]],
        grid(xs, [0.0, 0.5, 1.0], 0.0),
    )
    # Patch B: x in [5, 10], degree 2, curving up in z.
    b = patch(
        [2, 2],
        [[0, 0, 0, 5, 5, 5], [0, 0, 0, 1, 1, 1]],
        [
            (5 + x, y, 0.5 * x * x / 25)
            for y in (0.0, 0.5, 1.0)
            for x in (0.0, 2.5, 5.0)
        ],
    )

    def box(u1, v1, first):
        c = [(0, 0), (u1, 0), (u1, v1), (0, v1)]
        # Trims are Points to Kratos, which need no ids.
        return [trim(first + i, line(c[i], c[(i + 1) % 4], None)) for i in range(4)]

    la, lb = box(5, 1, 0), box(5, 1, 0)
    emb = trim(4, line((1.0, 0.25), (4.0, 0.75), ids))
    fa = face(
        2,
        0,
        loops=[{"type": "outer", "trims": la, "extra": ""}],
        ids=ids.take(15),
        embedded=[emb],
        extra='{"embedded_loops":[],"embedded_points":[]}',
    )
    fb = face(
        -1,
        1,
        loops=[{"type": "outer", "trims": lb, "extra": ""}],
        name="right",
        ids=ids.take(9),
    )
    coupling = edge(
        10,
        curve(1, [0, 0, 1, 1], [(5, 0), (5, 1)], active=[0, 1], ids=ids.take(2)),
        [ref(2, 1, 1), ref(-1, 3, 0, name="right")],
    )
    free = edge(11, None, [ref(2, 0, 1)])
    vertex = {"id": 12, "name": "", "topology": [ref(2, local=[2.5, 0.5])], "extra": ""}
    return {
        "patches": [a, b],
        "breps": [
            {
                "id": 1,
                "name": "",
                "faces": [0],
                "edges": [1],
                "vertices": [0],
                "extra": "",
            },
            {
                "id": -1,
                "name": "beam_right",
                "faces": [1],
                "edges": [0],
                "vertices": [],
                "extra": '{"layer":"beam"}',
            },
        ],
        "faces": [fa, fb],
        "edges": [coupling, free],
        "vertices": [vertex],
        "extra": '{"tolerances":{"model_tolerance":0.0001},"version_number":1,"meshio++":"fixture"}',
    }


def import_square():
    def jl(ident, deg, knots, pts, w, layer):
        return {
            "Layer": layer,
            "Weights": w,
            "CPCoordinates": pts,
            "knotVector": knots,
            "id": ident,
            "pDegree": deg,
        }

    doc = {
        "Surfaces": [
            {
                "Weights": [0, 0, 0, 0],
                "CPCoordinates": [
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [1.0, 1.0, 0.0],
                ],
                "knotVectors": {
                    "Xi": [0.0, 0.0, 1.0, 1.0],
                    "Eta": [0.0, 0.0, 1.0, 1.0],
                },
                "id": 1,
                "pDegreeU": 1,
                "pDegreeV": 1,
            }
        ],
        "Lines": [
            jl(
                1,
                1,
                [0.0, 0.0, 1.0, 1.0],
                [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
                [1, 1],
                "bottom",
            ),
            jl(
                2,
                1,
                [0.0, 0.0, 1.0, 1.0],
                [[1.0, 0.0, 0.0], [1.0, 1.0, 0.0]],
                [1, 1],
                "right",
            ),
            jl(
                3,
                1,
                [0.0, 0.0, 1.0, 1.0],
                [[1.0, 1.0, 0.0], [0.0, 1.0, 0.0]],
                [1, 1],
                "top",
            ),
            jl(
                4,
                1,
                [0.0, 0.0, 1.0, 1.0],
                [[0.0, 1.0, 0.0], [0.0, 0.0, 0.0]],
                [1, 1],
                "left",
            ),
            jl(
                5,
                2,
                [0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
                [[0.75, 0.5, 0.0], [0.75, 0.75, 0.0], [0.5, 0.75, 0.0]],
                [1.0, W, 1.0],
                "hole",
            ),
        ],
    }
    return json.dumps(doc, indent=4) + "\n"


# The control polygon of the unit circle in four rational quarters, from
# angle 0 counter-clockwise, and its knots.
CIRCLE = [(1, 0), (1, 1), (0, 1), (-1, 1), (-1, 0), (-1, -1), (0, -1), (1, -1), (1, 0)]
CIRCLE_W = [1, W, 1, W, 1, W, 1, W, 1]
CIRCLE_KNOTS = [0, 0, 0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1, 1, 1]


def circle_curve(centre, r, z=None, sense=1, ids=None):
    """The full circle about ``centre`` (2-D, or 3-D at height ``z``) from
    angle 0, counter-clockwise for ``sense`` 1 and clockwise for -1."""
    pts = []
    for x, y in CIRCLE:
        q = (centre[0] + r * x, centre[1] + sense * r * y)
        pts.append(q if z is None else (q[0], q[1], z))
    return curve(
        2,
        CIRCLE_KNOTS,
        pts,
        weights=CIRCLE_W,
        active=[0, 1],
        ids=None if ids is None else ids.take(len(pts)),
    )


def circle_from(centre, r, a0, sense=1):
    """The full circle about the 2-D ``centre`` from angle ``a0`` in four
    rational quarters, counter-clockwise for ``sense`` 1."""
    pts = []
    for k in range(9):
        a = a0 + sense * k * math.pi / 4
        rk = r if k % 2 == 0 else r / W
        pts.append((centre[0] + rk * math.cos(a), centre[1] + rk * math.sin(a)))
    pts[8] = pts[0]
    return curve(2, CIRCLE_KNOTS, pts, weights=CIRCLE_W, active=[0, 1])


def _closed_solid(patches, faces, edges):
    return {
        "patches": patches,
        "breps": [
            {
                "id": 1,
                "name": "",
                "faces": list(range(len(faces))),
                "edges": list(range(len(edges))),
                "vertices": [],
                "extra": "",
            }
        ],
        "faces": faces,
        "edges": edges,
        "vertices": [],
        "extra": '{"version_number":1}',
    }


def _loop(trims):
    return {"type": "outer", "trims": trims, "extra": ""}


def unit_cube():
    # Each face: its origin and the two directions of u and v, u x v outward.
    sides = [
        ((0, 0, 0), (0, 1, 0), (1, 0, 0)),
        ((0, 0, 1), (1, 0, 0), (0, 1, 0)),
        ((0, 0, 0), (1, 0, 0), (0, 0, 1)),
        ((0, 1, 0), (0, 0, 1), (1, 0, 0)),
        ((0, 0, 0), (0, 0, 1), (0, 1, 0)),
        ((1, 0, 0), (0, 1, 0), (0, 0, 1)),
    ]
    square = [(0, 0), (1, 0), (1, 1), (0, 1)]
    ids = Ids()
    patches, faces, refs = [], [], {}

    def point(side, u, v):
        o, a, b = side
        return tuple(o[k] + u * a[k] + v * b[k] for k in range(3))

    for f, side in enumerate(sides):

        def at(u, v, side=side):
            return point(side, u, v)

        patches.append(
            patch(
                [1, 1],
                [[0, 0, 1, 1], [0, 0, 1, 1]],
                [at(u, v) for v in (0, 1) for u in (0, 1)],
            )
        )
        trims = []
        for k in range(4):
            p, q = square[k], square[(k + 1) % 4]
            trims.append(trim(k, line(p, q, None)))
            pa, pb = at(*p), at(*q)
            key = tuple(sorted((pa, pb)))
            refs.setdefault(key, []).append(ref(2 + f, k, 1 if (pa, pb) == key else 0))
        faces.append(face(2 + f, f, loops=[_loop(trims)], ids=ids.take(4)))
    edges = [
        edge(
            8 + i,
            curve(1, [0, 0, 1, 1], list(key), active=[0, 1], ids=ids.take(2)),
            refs[key],
        )
        for i, key in enumerate(sorted(refs))
    ]
    return _closed_solid(patches, faces, edges)


def capped_cylinder(top_start=0.0):
    r, h = 1.0, 2.0
    band = patch(
        [2, 1],
        [CIRCLE_KNOTS, [0, 0, 1, 1]],
        [(r * x, r * y, z) for z in (0.0, h) for x, y in CIRCLE],
        CIRCLE_W * 2,
    )
    # The caps as planes over [-r, r]^2: the bottom one's v runs along -y so
    # its normal is -z; in its parameters the circle runs counter-clockwise,
    # which is clockwise seen from above.
    bottom = patch(
        [1, 1],
        [[0, 0, 1, 1], [0, 0, 1, 1]],
        [(-r, r, 0), (r, r, 0), (-r, -r, 0), (r, -r, 0)],
    )
    top = patch(
        [1, 1],
        [[0, 0, 1, 1], [0, 0, 1, 1]],
        [(-r, -r, h), (r, -r, h), (-r, r, h), (r, r, h)],
    )
    band_loop = [
        trim(0, line((0, 0), (1, 0), None)),
        trim(1, line((1, 0), (1, 1), None)),
        trim(2, line((1, 1), (0, 1), None)),
        trim(3, line((0, 1), (0, 0), None)),
    ]
    ids = Ids()
    faces = [
        face(2, 0, loops=[_loop(band_loop)], rational=True, ids=ids.take(18)),
        face(
            3,
            1,
            loops=[_loop([trim(0, circle_curve((0.5, 0.5), 0.5))])],
            ids=ids.take(4),
        ),
        face(
            4,
            2,
            loops=[
                _loop(
                    [
                        trim(
                            0,
                            (
                                circle_from((0.5, 0.5), 0.5, top_start)
                                if top_start
                                else circle_curve((0.5, 0.5), 0.5)
                            ),
                        )
                    ]
                )
            ],
            ids=ids.take(4),
        ),
    ]
    seam = curve(
        1, [0, 0, 1, 1], [(r, 0, 0), (r, 0, h)], active=[0, 1], ids=ids.take(2)
    )
    edges = [
        edge(5, circle_curve((0, 0), r, 0.0, ids=ids), [ref(2, 0, 1), ref(3, 0, 0)]),
        edge(6, circle_curve((0, 0), r, h, ids=ids), [ref(2, 2, 0), ref(4, 0, 1)]),
        edge(7, seam, [ref(2, 1, 1), ref(2, 3, 0)]),
    ]
    return _closed_solid([band, bottom, top], faces, edges)


def sphere():
    r = 1.0
    # Profiles (radius, height, weight) from the equator up to the north
    # pole, and from the south pole up to the equator.
    north = [(r, 0.0, 1.0), (r, r, W), (0.0, r, 1.0)]
    south = [(0.0, -r, 1.0), (r, -r, W), (r, 0.0, 1.0)]
    patches = []
    for prof in (north, south):
        pts, w = [], []
        for rho, z, wp in prof:
            for (x, y), wc in zip(CIRCLE, CIRCLE_W):
                pts.append((rho * x, rho * y, z))
                w.append(wc * wp)
        patches.append(patch([2, 2], [CIRCLE_KNOTS, [0, 0, 0, 1, 1, 1]], pts, w))
    loop = [
        trim(0, line((0, 0), (1, 0), None)),
        trim(1, line((1, 0), (1, 1), None)),
        trim(2, line((1, 1), (0, 1), None)),
        trim(3, line((0, 1), (0, 0), None)),
    ]
    ids = Ids()
    faces = [
        face(2, 0, loops=[_loop(loop)], rational=True, ids=ids.take(27)),
        face(3, 1, loops=[_loop(loop)], rational=True, ids=ids.take(27)),
    ]
    meridian = curve(
        2,
        [0, 0, 0, 1, 1, 1],
        [(r, 0, 0), (r, 0, r), (0, 0, r)],
        weights=[1, W, 1],
        active=[0, 1],
        ids=ids.take(3),
    )
    meridian_s = curve(
        2,
        [0, 0, 0, 1, 1, 1],
        [(0, 0, -r), (r, 0, -r), (r, 0, 0)],
        weights=[1, W, 1],
        active=[0, 1],
        ids=ids.take(3),
    )
    edges = [
        # The equator: the north face's trim 0, the south face's trim 2.
        edge(4, circle_curve((0, 0), r, 0.0, ids=ids), [ref(2, 0, 1), ref(3, 2, 0)]),
        edge(5, meridian, [ref(2, 1, 1), ref(2, 3, 0)]),
        edge(6, None, [ref(2, 2, 1)]),
        edge(7, meridian_s, [ref(3, 1, 1), ref(3, 3, 0)]),
        edge(8, None, [ref(3, 0, 1)]),
    ]
    return _closed_solid(patches, faces, edges)


def capped_cylinder_offset():
    return capped_cylinder(top_start=math.pi / 2)


def plate_with_near_hole():
    s = 2.5
    plate = patch(
        [1, 1],
        [[0, 0, 1, 1], [0, 0, 1, 1]],
        [(x, y, 0.0) for y in (-s, s) for x in (-s, s)],
    )

    def uv(x, y):
        return ((x + s) / (2 * s), (y + s) / (2 * s))

    rim = circle_from(uv(0.0, 0.0), 2.0 / (2 * s), 1.25 * math.pi)
    hole = circle_from(uv(0.98, 0.0), 1.0 / (2 * s), 0.0)
    loops = [
        _loop([trim(0, rim)]),
        {"type": "inner", "trims": [trim(1, hole, direction=False)], "extra": ""},
    ]
    return _closed_solid([plate], [face(2, 0, loops=loops, ids=Ids().take(4))], [])


def split_rim(d=0.005):
    # The top plate (z = 0, normal +z) and the front one (y = 0, z in
    # [-1, 0], normal -y); their loops, counter-clockwise in (u, v).
    sides = [((0, 0, 0), (1, 0, 0), (0, 1, 0)), ((0, 0, -1), (1, 0, 0), (0, 0, 1))]
    loops = [
        [(0, 0), (0.5, 0), (0.5 + d, 0), (1, 0), (1, 1), (0, 1)],
        [(0, 0), (1, 0), (1, 1), (0.5 + d, 1), (0.5, 1), (0, 1)],
    ]
    ids = Ids()
    patches, faces, refs = [], [], {}
    for f, (side, loop) in enumerate(zip(sides, loops)):

        def at(u, v, side=side):
            o, a, b = side
            return tuple(float(o[k] + u * a[k] + v * b[k]) for k in range(3))

        patches.append(
            patch(
                [1, 1],
                [[0, 0, 1, 1], [0, 0, 1, 1]],
                [at(u, v) for v in (0, 1) for u in (0, 1)],
            )
        )
        trims = []
        for k in range(len(loop)):
            p, q = loop[k], loop[(k + 1) % len(loop)]
            trims.append(trim(k, line(p, q, None)))
            pa, pb = at(*p), at(*q)
            key = tuple(sorted((pa, pb)))
            refs.setdefault(key, []).append(ref(2 + f, k, 1 if (pa, pb) == key else 0))
        faces.append(face(2 + f, f, loops=[_loop(trims)], ids=ids.take(4)))
    edges = [
        edge(
            8 + i,
            curve(1, [0, 0, 1, 1], list(key), active=[0, 1], ids=ids.take(2)),
            refs[key],
        )
        for i, key in enumerate(sorted(refs))
    ]
    return _closed_solid(patches, faces, edges)


def sliver_strip(w=0.01):
    # [0, 1] x [0, 1], the sliver [0, 1] x [1, 1 + w] and [0, 1] x [1 + w,
    # 2 + w], each its own planar patch, normal +z, sharing their edges.
    rows = [(0.0, 1.0), (1.0, 1.0 + w), (1.0 + w, 2.0 + w)]
    ids = Ids()
    patches, faces, refs = [], [], {}
    for f, (y0, y1) in enumerate(rows):

        def at(u, v, y0=y0, y1=y1):
            return (float(u), float(y0 + v * (y1 - y0)), 0.0)

        patches.append(
            patch(
                [1, 1],
                [[0, 0, 1, 1], [0, 0, 1, 1]],
                [at(u, v) for v in (0, 1) for u in (0, 1)],
            )
        )
        loop = [(0, 0), (1, 0), (1, 1), (0, 1)]
        trims = []
        for k in range(len(loop)):
            p, q = loop[k], loop[(k + 1) % len(loop)]
            trims.append(trim(k, line(p, q, None)))
            pa, pb = at(*p), at(*q)
            key = tuple(sorted((pa, pb)))
            refs.setdefault(key, []).append(ref(2 + f, k, 1 if (pa, pb) == key else 0))
        faces.append(face(2 + f, f, loops=[_loop(trims)], ids=ids.take(4)))
    edges = [
        edge(
            8 + i,
            curve(1, [0, 0, 1, 1], list(key), active=[0, 1], ids=ids.take(2)),
            refs[key],
        )
        for i, key in enumerate(sorted(refs))
    ]
    return _closed_solid(patches, faces, edges)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for name, build in (
        ("plate_with_hole.cad.json", plate_with_hole),
        ("cylindrical_roof.cad.json", cylindrical_roof),
        ("two_patch_beam.cad.json", two_patch_beam),
        ("unit_cube.cad.json", unit_cube),
        ("capped_cylinder.cad.json", capped_cylinder),
        ("sphere.cad.json", sphere),
        ("capped_cylinder_offset.cad.json", capped_cylinder_offset),
        ("plate_with_near_hole.cad.json", plate_with_near_hole),
        ("split_rim.cad.json", split_rim),
        ("sliver_strip.cad.json", sliver_strip),
    ):
        (OUT / name).write_text(render(build()), encoding="utf-8")
    (OUT / "import_square.json").write_text(import_square(), encoding="utf-8")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()

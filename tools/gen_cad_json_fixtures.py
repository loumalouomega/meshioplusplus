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


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for name, build in (
        ("plate_with_hole.cad.json", plate_with_hole),
        ("cylindrical_roof.cad.json", cylindrical_roof),
        ("two_patch_beam.cad.json", two_patch_beam),
    ):
        (OUT / name).write_text(render(build()), encoding="utf-8")
    (OUT / "import_square.json").write_text(import_square(), encoding="utf-8")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()

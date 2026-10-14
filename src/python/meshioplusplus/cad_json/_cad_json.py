"""
Kratos CAD JSON (``*.cad.json``): the NURBS boundary representation Kratos
Multiphysics' IGA reads (``CadJsonInput``) and writes (``CadJsonOutput``) --
the Python reference of ``src/cpp/src/formats/cad_json.cpp``, which it matches
mesh for mesh, side channel for side channel and byte for byte.

The side channel is ``mesh.info["brep"]``, a dict laid out as ``BrepInfo``
(``brep.hpp``)::

    {"patches":  [{"degrees", "knots", "control_points", "weights", "fields"}],
     "breps":    [{"id", "name", "faces", "edges", "vertices", "extra"}],
     "faces":    [{"id", "name", "patch", "swapped_normal", "trimmed",
                   "rational", "reduced_knots", "control_point_ids", "loops",
                   "embedded_edges", "surface_extra", "extra"}],
     "edges":    [{"id", "name", "curve", "topology", "extra"}],
     "vertices": [{"id", "name", "topology", "extra"}],
     "extra":    str}

with loops ``{"type", "trims", "extra"}``, trims ``{"trim_index",
"direction", "curve", "extra"}``, curves ``{"degree", "knots",
"reduced_knots", "rational", "control_points", "weights", "ids",
"active_range", "extra"}`` and topology references ``{"id", "name",
"trim_index", "relative_direction", "local_coordinates", "extra"}``. See
doc/formats/cad_json.md.
"""

import math

import numpy as np

from .. import _json_text
from .. import _nurbs_eval as ne
from .._common import warn
from .._exceptions import ReadError, WriteError
from .._lagrange import vtk_lattice
from .._mesh import Mesh

_FMT = "cad_json"
_MAX_DEGREE = 32
_CORNER = 0.999999
# Automatic subdivision: a clipped face gets about this many cells per direction.
_AUTO_CELLS = 16
_JO = _json_text.JsonObject


# ---------------------------------------------------------------------------
# Reading: JSON -> brep dict
# ---------------------------------------------------------------------------


def _fail(what):
    raise ReadError(f"meshio++: cad_json: {what}")


def _is_number(v):
    return isinstance(v, (int, float)) and not isinstance(v, bool)


def _number(v, what):
    if not _is_number(v) or not math.isfinite(v):
        _fail(f"{what} is not a finite number")
    return float(v)


def _integer(v, what):
    if not _is_number(v):
        _fail(f"{what} is not a number")
    if isinstance(v, int):
        return v
    if not math.isfinite(v) or v != math.floor(v) or abs(v) > 9007199254740992.0:
        _fail(f"{what} is not an integer")
    return int(v)


def _bool(v, what):
    if not isinstance(v, bool):
        _fail(f"{what} is not a boolean")
    return v


def _array(v, what):
    if not isinstance(v, list) or isinstance(v, _JO):
        _fail(f"{what} is not an array")
    return v


def _object(v, what):
    if not isinstance(v, _JO):
        _fail(f"{what} is not an object")
    return v


def _extra(obj, known):
    rest = _JO((k, v) for k, v in obj if k not in known)
    return _json_text.dumps(rest, indent=0) if rest else ""


def _identity(obj, what):
    ident, name = -1, ""
    if "brep_id" in obj:
        ident = _integer(obj.find("brep_id"), f"{what} brep_id")
        if ident < 0:
            _fail(f"{what} brep_id is negative")
    if "brep_name" in obj:
        name = obj.find("brep_name")
        if not isinstance(name, str):
            _fail(f"{what} brep_name is not a string")
    return ident, name


def _vector(v, what):
    _array(v, what)
    return [_number(x, f"{what}[{i}]") for i, x in enumerate(v)]


def _check_knots(knots, what):
    for i in range(1, len(knots)):
        if knots[i] < knots[i - 1]:
            _fail(f"{what} decreases at entry {i}")


def _control_points(v, what):
    _array(v, what)
    xyz, w, ids = [], [], []
    form = -1
    for i, cp in enumerate(v):
        cw = f"{what}[{i}]"
        _array(cp, cw)
        if len(cp) == 2 and isinstance(cp[1], list) and not isinstance(cp[1], _JO):
            this_form = 1
            ids.append(_integer(cp[0], f"{cw} id"))
            coords = cp[1]
        else:
            this_form = 0
            coords = cp
            if len(cp) == 1 and isinstance(cp[0], list) and not isinstance(cp[0], _JO):
                coords = cp[0]
        if form >= 0 and this_form != form:
            _fail(f"{what} mixes control points with and without ids")
        form = this_form
        if not isinstance(coords, list) or isinstance(coords, _JO) or len(coords) != 4:
            _fail(f"{cw} needs [x, y, z, weight]")
        xyz.append([_number(coords[d], cw) for d in range(3)])
        w.append(_number(coords[3], f"{cw} weight"))
    return (
        np.array(xyz, dtype=np.float64).reshape(len(v), 3),
        np.array(w, dtype=np.float64),
        np.array(ids, dtype=np.int64),
    )


def _fit_knots(knots, reduced):
    return [knots[0]] + list(knots) + [knots[-1]] if reduced else list(knots)


def _degree(v, what):
    p = _integer(v, what)
    if p < 1 or p > _MAX_DEGREE:
        _fail(f"{what} {p} is outside 1..{_MAX_DEGREE}")
    return p


def _check_weights(w, rational, what):
    if rational and not all(x > 0.0 for x in w.tolist()):
        _fail(f"{what} has a weight that is not positive")


def _curve(obj, what):
    _object(obj, what)
    rational = False
    if "is_rational" in obj:
        rational = _bool(obj.find("is_rational"), f"{what} is_rational")
    if "degree" not in obj or "knot_vector" not in obj or "control_points" not in obj:
        _fail(f"{what} needs degree, knot_vector and control_points")
    degree = _degree(obj.find("degree"), f"{what} degree")
    knots = _vector(obj.find("knot_vector"), f"{what} knot_vector")
    xyz, w, ids = _control_points(obj.find("control_points"), f"{what} control_points")
    n, p, k = len(w), degree, len(knots)
    if n < p + 1:
        _fail(f"{what} has fewer than degree + 1 control points")
    if k + 1 == n + p:
        reduced = True
    elif k == n + p + 1:
        reduced = False
    else:
        _fail(f"{what} has {k} knots for {n} control points of degree {p}")
    _check_knots(knots, f"{what} knot_vector")
    knots = _fit_knots(knots, reduced)
    _check_weights(w, rational, what)
    active = []
    if "active_range" in obj:
        active = _vector(obj.find("active_range"), f"{what} active_range")
        if len(active) != 2:
            _fail(f"{what} active_range needs two values")
    return {
        "degree": degree,
        "knots": np.array(knots, dtype=np.float64),
        "reduced_knots": reduced,
        "rational": rational,
        "control_points": xyz,
        "weights": w,
        "ids": ids,
        "active_range": active,
        "extra": _extra(
            obj,
            ("is_rational", "degree", "knot_vector", "active_range", "control_points"),
        ),
    }


def _trim(obj, what):
    _object(obj, what)
    t = {"trim_index": -1, "direction": True}
    if "trim_index" in obj:
        t["trim_index"] = _integer(obj.find("trim_index"), f"{what} trim_index")
        if t["trim_index"] < 0:
            _fail(f"{what} trim_index is negative")
    if "curve_direction" in obj:
        t["direction"] = _bool(obj.find("curve_direction"), f"{what} curve_direction")
    if "parameter_curve" not in obj:
        _fail(f"{what} has no parameter_curve")
    t["curve"] = _curve(obj.find("parameter_curve"), f"{what} parameter_curve")
    t["extra"] = _extra(obj, ("trim_index", "curve_direction", "parameter_curve"))
    return t


def _ref(obj, what):
    _object(obj, what)
    ident, name = _identity(obj, what)
    r = {
        "id": ident,
        "name": name,
        "trim_index": -1,
        "relative_direction": -1,
        "local_coordinates": [],
    }
    if "trim_index" in obj:
        r["trim_index"] = _integer(obj.find("trim_index"), f"{what} trim_index")
        if r["trim_index"] < 0:
            _fail(f"{what} trim_index is negative")
    if "relative_direction" in obj:
        r["relative_direction"] = int(
            _bool(obj.find("relative_direction"), f"{what} relative_direction")
        )
    if "local_coordinates" in obj:
        r["local_coordinates"] = _vector(
            obj.find("local_coordinates"), f"{what} local_coordinates"
        )
    r["extra"] = _extra(
        obj,
        (
            "brep_id",
            "brep_name",
            "trim_index",
            "relative_direction",
            "local_coordinates",
        ),
    )
    return r


def _topology(obj, what):
    if "topology" not in obj:
        return []
    topo = _array(obj.find("topology"), f"{what} topology")
    return [_ref(r, f"{what} topology[{i}]") for i, r in enumerate(topo)]


def _patch_counts(p, q, ku, kv, n, what):
    if ku + 1 > p and kv + 1 > q and (ku + 1 - p) * (kv + 1 - q) == n:
        return True
    if ku > p + 1 and kv > q + 1 and (ku - p - 1) * (kv - q - 1) == n:
        return False
    _fail(
        f"{what}: {n} control points do not match knot vectors of {ku} and {kv} "
        f"knots at degrees {p}, {q}"
    )


def _patch(degrees, knots, xyz, w):
    return {
        "degrees": np.array(degrees, dtype=np.int32),
        "knots": [np.array(k, dtype=np.float64) for k in knots],
        "control_points": xyz,
        "weights": w,
        "fields": [],
    }


def _check_patch_size(degrees, knots, what):
    for d in range(2):
        if len(knots[d]) - degrees[d] - 1 < degrees[d] + 1:
            _fail(f"{what} has fewer than degree + 1 control points in a direction")


def _face(obj, what, brep):
    _object(obj, what)
    ident, name = _identity(obj, what)
    face = {"id": ident, "name": name, "swapped_normal": False}
    if "swapped_surface_normal" in obj:
        face["swapped_normal"] = _bool(
            obj.find("swapped_surface_normal"), f"{what} swapped_surface_normal"
        )
    if "surface" not in obj:
        _fail(f"{what} has no surface")
    sw = f"{what} surface"
    s = _object(obj.find("surface"), sw)
    face["trimmed"] = (
        _bool(s.find("is_trimmed"), f"{sw} is_trimmed") if "is_trimmed" in s else False
    )
    face["rational"] = (
        _bool(s.find("is_rational"), f"{sw} is_rational")
        if "is_rational" in s
        else False
    )
    if "degrees" not in s or "knot_vectors" not in s or "control_points" not in s:
        _fail(f"{sw} needs degrees, knot_vectors and control_points")
    degs = _array(s.find("degrees"), f"{sw} degrees")
    kvs = _array(s.find("knot_vectors"), f"{sw} knot_vectors")
    if len(degs) != 2 or len(kvs) != 2:
        _fail(f"{sw} needs two degrees and two knot vectors")
    degrees, knots = [], []
    for d in range(2):
        degrees.append(_degree(degs[d], f"{sw} degree"))
        knots.append(_vector(kvs[d], f"{sw} knot vector"))
        _check_knots(knots[-1], f"{sw} knot vector")
    xyz, w, ids = _control_points(s.find("control_points"), f"{sw} control_points")
    reduced = _patch_counts(
        degrees[0], degrees[1], len(knots[0]), len(knots[1]), len(w), sw
    )
    knots = [_fit_knots(k, reduced) for k in knots]
    _check_patch_size(degrees, knots, sw)
    _check_weights(w, face["rational"], sw)
    face["reduced_knots"] = reduced
    face["control_point_ids"] = ids
    face["surface_extra"] = _extra(
        s, ("is_trimmed", "is_rational", "degrees", "knot_vectors", "control_points")
    )
    loops = []
    if "boundary_loops" in obj:
        for li, lo in enumerate(
            _array(obj.find("boundary_loops"), f"{what} boundary_loops")
        ):
            lw = f"{what} boundary_loops[{li}]"
            _object(lo, lw)
            if "loop_type" not in lo:
                _fail(f"{lw} has no loop_type")
            if not isinstance(lo.find("loop_type"), str):
                _fail(f"{lw} loop_type is not a string")
            if "trimming_curves" not in lo:
                _fail(f"{lw} has no trimming_curves")
            trims = _array(lo.find("trimming_curves"), f"{lw} trimming_curves")
            loops.append(
                {
                    "type": lo.find("loop_type"),
                    "trims": [
                        _trim(t, f"{lw} trimming_curves[{ti}]")
                        for ti, t in enumerate(trims)
                    ],
                    "extra": _extra(lo, ("loop_type", "trimming_curves")),
                }
            )
    face["loops"] = loops
    emb = []
    if "embedded_edges" in obj:
        emb = [
            _trim(t, f"{what} embedded_edges[{ti}]")
            for ti, t in enumerate(
                _array(obj.find("embedded_edges"), f"{what} embedded_edges")
            )
        ]
    face["embedded_edges"] = emb
    face["extra"] = _extra(
        obj,
        (
            "brep_id",
            "brep_name",
            "swapped_surface_normal",
            "surface",
            "boundary_loops",
            "embedded_edges",
        ),
    )
    face["patch"] = len(brep["patches"])
    brep["patches"].append(_patch(degrees, knots, xyz, w))
    brep["faces"].append(face)


def _empty():
    return {
        "patches": [],
        "breps": [],
        "faces": [],
        "edges": [],
        "vertices": [],
        "extra": "",
    }


def _read_breps(root, brep):
    for bi, bo in enumerate(_array(root.find("breps"), "breps")):
        bw = f"breps[{bi}]"
        _object(bo, bw)
        ident, name = _identity(bo, bw)
        entry = {"id": ident, "name": name, "faces": [], "edges": [], "vertices": []}
        if "faces" in bo:
            for i, f in enumerate(_array(bo.find("faces"), f"{bw} faces")):
                entry["faces"].append(len(brep["faces"]))
                _face(f, f"{bw} faces[{i}]", brep)
        if "edges" in bo:
            for i, eo in enumerate(_array(bo.find("edges"), f"{bw} edges")):
                ew = f"{bw} edges[{i}]"
                _object(eo, ew)
                eid, ename = _identity(eo, ew)
                curve = (
                    _curve(eo.find("3d_curve"), f"{ew} 3d_curve")
                    if "3d_curve" in eo
                    else None
                )
                entry["edges"].append(len(brep["edges"]))
                brep["edges"].append(
                    {
                        "id": eid,
                        "name": ename,
                        "curve": curve,
                        "topology": _topology(eo, ew),
                        "extra": _extra(
                            eo, ("brep_id", "brep_name", "3d_curve", "topology")
                        ),
                    }
                )
        if "vertices" in bo:
            for i, vo in enumerate(_array(bo.find("vertices"), f"{bw} vertices")):
                vw = f"{bw} vertices[{i}]"
                _object(vo, vw)
                vid, vname = _identity(vo, vw)
                entry["vertices"].append(len(brep["vertices"]))
                brep["vertices"].append(
                    {
                        "id": vid,
                        "name": vname,
                        "topology": _topology(vo, vw),
                        "extra": _extra(vo, ("brep_id", "brep_name", "topology")),
                    }
                )
        entry["extra"] = _extra(
            bo, ("brep_id", "brep_name", "faces", "edges", "vertices")
        )
        brep["breps"].append(entry)
    brep["extra"] = _extra(root, ("breps",))


def _import_net(obj, what):
    if "CPCoordinates" not in obj:
        _fail(f"{what} has no CPCoordinates")
    coords = _array(obj.find("CPCoordinates"), f"{what} CPCoordinates")
    xyz = []
    for c in coords:
        v = _vector(c, f"{what} CPCoordinates")
        if len(v) < 2 or len(v) > 3:
            _fail(f"{what} CPCoordinates need two or three values")
        xyz.append(v + [0.0] * (3 - len(v)))
    n = len(coords)
    if "Weights" not in obj:
        return np.array(xyz, dtype=np.float64).reshape(n, 3), np.ones(n)
    ws = _vector(obj.find("Weights"), f"{what} Weights")
    if len(ws) != n:
        _fail(f"{what} has {len(ws)} weights for {n} control points")
    return np.array(xyz, dtype=np.float64).reshape(n, 3), np.array(ws, dtype=np.float64)


def _import_rational(w):
    vals = w.tolist()
    return not all(x == 1.0 for x in vals) and not all(x == 0.0 for x in vals)


def _read_import(root, brep):
    entry = {
        "id": -1,
        "name": "",
        "faces": [],
        "edges": [],
        "vertices": [],
        "extra": "",
    }
    if "Surfaces" in root:
        for i, so in enumerate(_array(root.find("Surfaces"), "Surfaces")):
            sw = f"Surfaces[{i}]"
            _object(so, sw)
            ident = _integer(so.find("id"), f"{sw} id") if "id" in so else -1
            if "pDegreeU" not in so or "pDegreeV" not in so or "knotVectors" not in so:
                _fail(f"{sw} needs pDegreeU, pDegreeV and knotVectors")
            degrees = [
                _degree(so.find("pDegreeU"), f"{sw} pDegreeU"),
                _degree(so.find("pDegreeV"), f"{sw} pDegreeV"),
            ]
            kvs = _object(so.find("knotVectors"), f"{sw} knotVectors")
            if "Xi" not in kvs or "Eta" not in kvs:
                _fail(f"{sw} knotVectors needs Xi and Eta")
            knots = [
                _vector(kvs.find("Xi"), f"{sw} Xi"),
                _vector(kvs.find("Eta"), f"{sw} Eta"),
            ]
            for k in knots:
                _check_knots(k, f"{sw} knot vector")
            xyz, w = _import_net(so, sw)
            n = len(w)
            p, q = degrees
            ku, kv = len(knots[0]), len(knots[1])
            if ku + 1 > p and kv + 1 > q and (ku + 1 - p) * (kv + 1 - q) == n:
                reduced = True
            elif ku > p + 1 and kv > q + 1 and (ku - p - 1) * (kv - q - 1) == n:
                reduced = False
            else:
                _fail(f"{sw}: the control points do not match the knot vectors")
            knots = [_fit_knots(k, reduced) for k in knots]
            _check_patch_size(degrees, knots, sw)
            rational = _import_rational(w)
            _check_weights(w, rational, sw)
            entry["faces"].append(len(brep["faces"]))
            brep["faces"].append(
                {
                    "id": ident,
                    "name": "",
                    "patch": len(brep["patches"]),
                    "swapped_normal": False,
                    "trimmed": False,
                    "rational": rational,
                    "reduced_knots": reduced,
                    "control_point_ids": np.zeros(0, dtype=np.int64),
                    "loops": [],
                    "embedded_edges": [],
                    "surface_extra": "",
                    "extra": _extra(
                        so,
                        (
                            "id",
                            "pDegreeU",
                            "pDegreeV",
                            "knotVectors",
                            "CPCoordinates",
                            "Weights",
                        ),
                    ),
                }
            )
            brep["patches"].append(_patch(degrees, knots, xyz, w))
    if "Lines" in root:
        for i, lo in enumerate(_array(root.find("Lines"), "Lines")):
            lw = f"Lines[{i}]"
            _object(lo, lw)
            ident = _integer(lo.find("id"), f"{lw} id") if "id" in lo else -1
            if "pDegree" not in lo or "knotVector" not in lo:
                _fail(f"{lw} needs pDegree and knotVector")
            degree = _degree(lo.find("pDegree"), f"{lw} pDegree")
            knots = _vector(lo.find("knotVector"), f"{lw} knotVector")
            _check_knots(knots, f"{lw} knotVector")
            xyz, w = _import_net(lo, lw)
            n, p = len(w), degree
            if n < p + 1:
                _fail(f"{lw} has fewer than degree + 1 control points")
            if len(knots) + 1 == n + p:
                reduced = True
            elif len(knots) == n + p + 1:
                reduced = False
            else:
                _fail(f"{lw}: the control points do not match the knot vector")
            _check_weights(w, True, lw)
            entry["edges"].append(len(brep["edges"]))
            brep["edges"].append(
                {
                    "id": ident,
                    "name": "",
                    "curve": {
                        "degree": degree,
                        "knots": np.array(_fit_knots(knots, reduced), dtype=np.float64),
                        "reduced_knots": reduced,
                        "rational": True,
                        "control_points": xyz,
                        "weights": w,
                        "ids": np.zeros(0, dtype=np.int64),
                        "active_range": [],
                        "extra": "",
                    },
                    "topology": [],
                    "extra": _extra(
                        lo, ("id", "pDegree", "knotVector", "CPCoordinates", "Weights")
                    ),
                }
            )
    brep["breps"].append(entry)
    brep["extra"] = _extra(root, ("Surfaces", "Lines"))


def parse_brep(data):
    """The ``brep`` dict of a CAD JSON text (``bytes`` or ``str``)."""
    root = _json_text.parse(data, _FMT)
    if not isinstance(root, _JO):
        _fail("the document is not an object")
    brep = _empty()
    if "breps" in root:
        _read_breps(root, brep)
    elif "Surfaces" in root or "Lines" in root:
        _read_import(root, brep)
    else:
        _fail(
            'neither a "breps" list (CadJsonInput) nor "Surfaces"/"Lines" '
            "(ImportNurbsSbmModeler)"
        )
    return brep


def read_brep(filename):
    """Only the ``brep`` dict of a CAD JSON file (no sampling)."""
    with open(filename, "rb") as f:
        return parse_brep(f.read())


# ---------------------------------------------------------------------------
# Sampling: brep dict -> Mesh
# ---------------------------------------------------------------------------


def _bad(what):
    raise ValueError(f"meshio++: sample_brep: {what}")


def _quad_type(k):
    return "quad" if k <= 1 else "quad9" if k == 2 else "VTK_LAGRANGE_QUADRILATERAL"


def _line_type(k):
    return "line" if k <= 1 else "line3" if k == 2 else "VTK_LAGRANGE_CURVE"


def _axis(p, knots, n, s, active):
    lo, hi, span, span_no = [], [], [], []
    for si, i in enumerate(ne.spans(p, knots, n)):
        k0, k1 = knots[i], knots[i + 1]
        if active is not None and len(active) == 2:
            r0, r1 = min(active[0], active[1]), max(active[0], active[1])
            k0, k1 = max(k0, r0), min(k1, r1)
            if not k0 < k1:
                continue
        for j in range(s):
            lo.append(k0 + (k1 - k0) * (j / s))
            hi.append(k1 if j + 1 == s else k0 + (k1 - k0) * ((j + 1) / s))
            span.append(i)
            span_no.append(si)
    return lo, hi, span, span_no


def _lattice(lo, hi, k):
    n = len(lo)
    g = [0.0] * (n * k + 1)
    for c in range(n):
        for i in range(k):
            g[c * k + i] = lo[c] + (hi[c] - lo[c]) * (i / k)
    if n:
        g[n * k] = hi[n - 1]
    return g


def _curve_rows(c):
    return [
        [float(x) for x in row]
        for row in np.asarray(c["control_points"], dtype=np.float64)
    ]


def _weights(c, rational):
    return (
        [float(x) for x in np.asarray(c["weights"], dtype=np.float64)]
        if rational
        else None
    )


def _check_curve(c, what):
    deg = int(c["degree"])
    w = np.asarray(c["weights"])
    xyz = np.asarray(c["control_points"])
    n = len(w)
    if deg < 1 or deg > _MAX_DEGREE:
        _bad(f"{what} degree is out of range")
    if xyz.shape != (n, 3):
        _bad(f"{what} control points are not Float64 (n, 3) beside n weights")
    knots = [float(x) for x in c["knots"]]
    if n < deg + 1 or len(knots) != n + deg + 1:
        _bad(f"{what} knot vector does not match its control points")
    if any(not knots[i] >= knots[i - 1] for i in range(1, len(knots))):
        _bad(f"{what} knot vector decreases")
    ids = np.asarray(c.get("ids", np.zeros(0)))
    if ids.size and ids.size != n:
        _bad(f"{what} ids do not match its control points")
    return deg, knots


def _check_patch(patch, what):
    degrees = [int(d) for d in patch["degrees"]]
    knots = [[float(x) for x in k] for k in patch["knots"]]
    if len(degrees) != 2 or len(knots) != 2:
        _bad(f"{what} is not a surface patch (two degrees, two knot vectors)")
    n = 1
    for d in range(2):
        if degrees[d] < 1 or degrees[d] > _MAX_DEGREE:
            _bad(f"{what} degree is out of range")
        if len(knots[d]) < 2 * degrees[d] + 2:
            _bad(f"{what} knot vector is too short")
        if any(not knots[d][i] >= knots[d][i - 1] for i in range(1, len(knots[d]))):
            _bad(f"{what} knot vector decreases")
        n *= len(knots[d]) - degrees[d] - 1
    if (
        np.asarray(patch["control_points"]).shape != (n, 3)
        or len(patch["weights"]) != n
    ):
        _bad(f"{what} control net does not match its knot vectors")
    return degrees, knots


def _loop_polygon(loop, samples):
    poly = []
    for t in loop["trims"]:
        c = t["curve"]
        deg = int(c["degree"])
        knots = [float(x) for x in c["knots"]]
        n = len(c["weights"])
        lo, hi, span, _ = _axis(deg, knots, n, 1, [float(x) for x in c["active_range"]])
        rows = _curve_rows(c)
        w = _weights(c, bool(c["rational"]))
        pts = []
        for s in range(len(lo)):
            for j in range(samples):
                u = lo[s] + (hi[s] - lo[s]) * (j / samples)
                xyz = ne.curve_point_in_span(deg, knots, rows, w, span[s], u)
                pts.append((xyz[0], xyz[1]))
        if lo:
            s = len(lo) - 1
            xyz = ne.curve_point_in_span(deg, knots, rows, w, span[s], hi[s])
            pts.append((xyz[0], xyz[1]))
        if not t["direction"]:
            pts.reverse()
        poly.extend(pts)
    return poly


def _inside(loops, u, v):
    inside = False
    for poly in loops:
        n = len(poly)
        j = n - 1
        for i in range(n):
            a, b = poly[i], poly[j]
            if (a[1] > v) != (b[1] > v) and u < (b[0] - a[0]) * (v - a[1]) / (
                b[1] - a[1]
            ) + a[0]:
                inside = not inside
            j = i
    return inside


class _Block:
    def __init__(self, cell_type, nodes):
        self.type = cell_type
        self.nodes = nodes
        self.conn = []
        self.kind, self.index, self.ident, self.span, self.trimmed = [], [], [], [], []


def _block(blocks, cell_type, nodes):
    # Keyed by node count too: faces of degrees 3 and 4 are both
    # VTK_LAGRANGE_QUADRILATERAL, with 16 and 25 nodes.
    for b in blocks:
        if b.type == cell_type and b.nodes == nodes:
            return b
    blocks.append(_Block(cell_type, nodes))
    return blocks[-1]


def sample(brep, subdivisions=0, trim_samples=16, edges=True, trim=True):
    """Sample a ``brep`` dict into a :class:`Mesh`, as the reader does."""
    sub = int(subdivisions)
    if sub < 0 or sub > 64:
        _bad("subdivisions must be 0 (automatic) to 64")
    s = max(sub, 1)

    def auto_split(p, knots, n, clip):
        if sub > 0:
            return sub
        if not clip:
            return 1
        nspans = max(len(ne.spans(p, knots, n)), 1)
        return min((_AUTO_CELLS + nspans - 1) // nspans, 64)

    if int(trim_samples) < 2 or int(trim_samples) > 4096:
        _bad("trim_samples must be 2 to 4096")
    trim_samples = int(trim_samples)
    patches = [_check_patch(p, f"patch {i}") for i, p in enumerate(brep["patches"])]
    blocks = []
    pending = []  # (eval function, lattice indices)
    num_points = 0

    for fi, face in enumerate(brep["faces"]):
        fw = f"face {fi}"
        pi = int(face["patch"])
        if pi < 0 or pi >= len(patches):
            _bad(f"{fw} names no patch")
        (p, q), (ku, kv) = patches[pi]
        nu, nv = len(ku) - p - 1, len(kv) - q - 1
        for loop in face["loops"]:
            for t in loop["trims"]:
                _check_curve(t["curve"], f"{fw} trimming curve")
        k = max(p, q)
        clip = trim and face["trimmed"] and face["loops"]
        au = _axis(p, ku, nu, auto_split(p, ku, nu, clip), None)
        av = _axis(q, kv, nv, auto_split(q, kv, nv, clip), None)
        gu = _lattice(au[0], au[1], k)
        gv = _lattice(av[0], av[1], k)
        width = len(gu)
        span_u = [
            au[2][min(i // k, len(au[2]) - 1)] if au[2] else 0 for i in range(len(gu))
        ]
        span_v = [
            av[2][min(i // k, len(av[2]) - 1)] if av[2] else 0 for i in range(len(gv))
        ]
        loops = []
        if clip:
            for loop in face["loops"]:
                poly = _loop_polygon(loop, trim_samples)
                if len(poly) >= 3:
                    loops.append(poly)
        lattice = vtk_lattice("quad", k)
        local = {}
        new_points = []
        blk = _block(blocks, _quad_type(k), len(lattice))
        for jv in range(len(av[0])):
            for iu in range(len(au[0])):
                trimmed = 0
                if clip:
                    u0, u1, v0, v1 = au[0][iu], au[1][iu], av[0][jv], av[1][jv]
                    uc, vc = (u0 + u1) * 0.5, (v0 + v1) * 0.5
                    if not _inside(loops, uc, vc):
                        continue
                    # The corners pulled a hair towards the centre, so a cell
                    # whose side lies on a loop is not counted as crossed.
                    a0, a1 = uc + (u0 - uc) * _CORNER, uc + (u1 - uc) * _CORNER
                    b0, b1 = vc + (v0 - vc) * _CORNER, vc + (v1 - vc) * _CORNER
                    if not (
                        _inside(loops, a0, b0)
                        and _inside(loops, a1, b0)
                        and _inside(loops, a1, b1)
                        and _inside(loops, a0, b1)
                    ):
                        trimmed = 1
                for i, j, _ in lattice:
                    g = (jv * k + j) * width + iu * k + i
                    if g not in local:
                        local[g] = num_points + len(new_points)
                        new_points.append(g)
                    blk.conn.append(local[g])
                blk.kind.append(0)
                blk.index.append(fi)
                blk.ident.append(int(face["id"]))
                blk.span.append((au[3][iu], av[3][jv]))
                blk.trimmed.append(trimmed)
        num_points += len(new_points)
        patch = brep["patches"][pi]
        rows = [
            [float(x) for x in r]
            for r in np.asarray(patch["control_points"], dtype=np.float64)
        ]
        w = (
            [float(x) for x in np.asarray(patch["weights"], dtype=np.float64)]
            if face["rational"]
            else None
        )

        def face_points(
            new_points=new_points,
            p=p,
            q=q,
            ku=ku,
            kv=kv,
            rows=rows,
            w=w,
            nu=nu,
            gu=gu,
            gv=gv,
            span_u=span_u,
            span_v=span_v,
            width=width,
        ):
            out = []
            for g in new_points:
                gi, gj = g % width, g // width
                out.append(
                    ne.surface_point_in_span(
                        p,
                        q,
                        ku,
                        kv,
                        rows,
                        w,
                        nu,
                        span_u[gi],
                        span_v[gj],
                        gu[gi],
                        gv[gj],
                    )
                )
            return out

        pending.append(face_points)

    if edges:
        for ei, edge in enumerate(brep["edges"]):
            c = edge["curve"]
            if c is None:
                continue
            deg, knots = _check_curve(c, f"edge {ei} 3d_curve")
            k = deg
            lo, hi, span, span_no = _axis(
                deg, knots, len(c["weights"]), s, [float(x) for x in c["active_range"]]
            )
            gu = _lattice(lo, hi, k)
            span_u = [
                span[min(i // k, len(span) - 1)] if span else 0 for i in range(len(gu))
            ]
            lattice = vtk_lattice("line", k)
            local = {}
            new_points = []
            blk = _block(blocks, _line_type(k), len(lattice))
            for iu in range(len(lo)):
                for i, _, _ in lattice:
                    g = iu * k + i
                    if g not in local:
                        local[g] = num_points + len(new_points)
                        new_points.append(g)
                    blk.conn.append(local[g])
                blk.kind.append(1)
                blk.index.append(ei)
                blk.ident.append(int(edge["id"]))
                blk.span.append((span_no[iu], -1))
                blk.trimmed.append(0)
            num_points += len(new_points)
            rows = _curve_rows(c)
            w = _weights(c, bool(c["rational"]))

            def edge_points(
                new_points=new_points,
                deg=deg,
                knots=knots,
                rows=rows,
                w=w,
                span_u=span_u,
                gu=gu,
            ):
                return [
                    ne.curve_point_in_span(deg, knots, rows, w, span_u[g], gu[g])
                    for g in new_points
                ]

            pending.append(edge_points)

    pts = []
    for fn in pending:
        pts.extend(fn())
    points = np.array(pts, dtype=np.float64).reshape(num_points, 3)
    # A face trimmed away or an edge outside its range leaves no cells.
    blocks = [b for b in blocks if b.kind]
    cells = [
        (b.type, np.array(b.conn, dtype=np.int64).reshape(len(b.kind), b.nodes))
        for b in blocks
    ]
    cell_data = {}
    if blocks:
        cell_data = {
            "cad:id": [np.array(b.ident, dtype=np.int64) for b in blocks],
            "cad:index": [np.array(b.index, dtype=np.int64) for b in blocks],
            "cad:kind": [np.array(b.kind, dtype=np.int64) for b in blocks],
            "cad:span": [
                np.array(b.span, dtype=np.int64).reshape(len(b.kind), 2) for b in blocks
            ],
            "cad:trimmed": [np.array(b.trimmed, dtype=np.int64) for b in blocks],
        }
    return Mesh(points, cells, cell_data=cell_data)


def read(filename, subdivisions=0, trim_samples=16, edges=True, trim=True):
    """Read a CAD JSON file as its sampled mesh, the B-rep in ``mesh.info["brep"]``."""
    brep = read_brep(filename)
    mesh = sample(
        brep,
        subdivisions=subdivisions,
        trim_samples=trim_samples,
        edges=edges,
        trim=trim,
    )
    mesh.info = {"brep": brep}
    return mesh


# ---------------------------------------------------------------------------
# Writing: brep dict -> JSON
# ---------------------------------------------------------------------------


def _wfail(what):
    raise WriteError(f"meshio++: cad_json: {what}")


def _add_extra(obj, extra, what):
    if not extra:
        return
    try:
        parsed = _json_text.parse(extra, _FMT)
    except ReadError:
        _wfail(f"{what} extra is not JSON object text")
    if not isinstance(parsed, _JO):
        _wfail(f"{what} extra is not JSON object text")
    obj.extend(parsed)


def _add_identity(obj, ident, name):
    if int(ident) >= 0:
        obj.append(("brep_id", int(ident)))
    if name:
        obj.append(("brep_name", str(name)))


def _floats(values):
    return [float(x) for x in values]


def _written_knots(knots, reduced):
    knots = _floats(knots)
    return knots[1:-1] if reduced and len(knots) >= 2 else knots


def _net(xyz, w, ids, what):
    xyz = np.asarray(xyz)
    w = np.asarray(w)
    ids = np.asarray(ids) if ids is not None else np.zeros(0, dtype=np.int64)
    n = len(w)
    if xyz.shape != (n, 3) or xyz.dtype.kind != "f" or (w.size and w.dtype.kind != "f"):
        _wfail(f"{what} control points are not Float64 (n, 3) beside n weights")
    if ids.size and ids.size != n:
        _wfail(f"{what} ids do not match its control points")
    out = []
    for i in range(n):
        c = [float(xyz[i, 0]), float(xyz[i, 1]), float(xyz[i, 2]), float(w[i])]
        out.append([int(ids[i]), c] if ids.size else [c])
    return out


def _curve_json(c, what):
    try:
        _check_curve(c, what)
    except ValueError as exc:
        _wfail(str(exc))
    o = _JO()
    o.append(("is_rational", bool(c["rational"])))
    o.append(("degree", int(c["degree"])))
    o.append(("knot_vector", _written_knots(c["knots"], c["reduced_knots"])))
    active = list(c["active_range"])
    if active:
        if len(active) != 2:
            _wfail(f"{what} active range needs two values")
        o.append(("active_range", _floats(active)))
    o.append(
        ("control_points", _net(c["control_points"], c["weights"], c.get("ids"), what))
    )
    _add_extra(o, c["extra"], what)
    return o


def _trim_json(t, what):
    o = _JO()
    if int(t["trim_index"]) >= 0:
        o.append(("trim_index", int(t["trim_index"])))
    o.append(("curve_direction", bool(t["direction"])))
    o.append(("parameter_curve", _curve_json(t["curve"], f"{what} parameter curve")))
    _add_extra(o, t["extra"], what)
    return o


def _topology_json(refs, what):
    out = []
    for r in refs:
        o = _JO()
        _add_identity(o, r["id"], r["name"])
        if int(r["trim_index"]) >= 0:
            o.append(("trim_index", int(r["trim_index"])))
        if int(r["relative_direction"]) >= 0:
            o.append(("relative_direction", int(r["relative_direction"]) != 0))
        if len(r["local_coordinates"]):
            o.append(("local_coordinates", _floats(r["local_coordinates"])))
        _add_extra(o, r["extra"], what)
        out.append(o)
    return out


def _face_json(brep, face, what):
    pi = int(face["patch"])
    if pi < 0 or pi >= len(brep["patches"]):
        _wfail(f"{what} names no patch")
    patch = brep["patches"][pi]
    try:
        _check_patch(patch, f"{what} patch")
    except ValueError as exc:
        _wfail(str(exc))
    o = _JO()
    _add_identity(o, face["id"], face["name"])
    o.append(("swapped_surface_normal", bool(face["swapped_normal"])))
    s = _JO()
    s.append(("is_trimmed", bool(face["trimmed"])))
    s.append(("is_rational", bool(face["rational"])))
    s.append(("degrees", [int(d) for d in patch["degrees"]]))
    s.append(
        (
            "knot_vectors",
            [_written_knots(k, face["reduced_knots"]) for k in patch["knots"]],
        )
    )
    s.append(
        (
            "control_points",
            _net(
                patch["control_points"],
                patch["weights"],
                face["control_point_ids"],
                f"{what} surface",
            ),
        )
    )
    _add_extra(s, face["surface_extra"], f"{what} surface")
    o.append(("surface", s))
    loops = []
    for li, loop in enumerate(face["loops"]):
        lo = _JO()
        lo.append(("loop_type", str(loop["type"])))
        lo.append(
            (
                "trimming_curves",
                [_trim_json(t, f"{what} loop {li}") for t in loop["trims"]],
            )
        )
        _add_extra(lo, loop["extra"], f"{what} loop")
        loops.append(lo)
    o.append(("boundary_loops", loops))
    o.append(
        (
            "embedded_edges",
            [_trim_json(t, f"{what} embedded edge") for t in face["embedded_edges"]],
        )
    )
    _add_extra(o, face["extra"], what)
    return o


def render(brep):
    """The CAD JSON text of a ``brep`` dict, as the C++ writer writes it."""
    # Kratos makes model-part nodes of surface and 3-D curve control points,
    # which need ids; trimming curves' are points and may go without.
    bare = sum(
        1 for f in brep.get("faces", []) if np.asarray(f["control_point_ids"]).size == 0
    )
    bare += sum(
        1
        for e in brep.get("edges", [])
        if e["curve"] is not None and np.asarray(e["curve"].get("ids", [])).size == 0
    )
    if bare:
        warn(
            f"cad_json: {bare} surface or 3-D curve control nets have no ids; Kratos's "
            "CadJsonInput needs [id, [x, y, z, w]] for them"
        )
    root = _JO()
    _add_extra(root, brep.get("extra", ""), "the file")
    breps = list(brep.get("breps", []))
    faces, edges, verts = (
        brep.get("faces", []),
        brep.get("edges", []),
        brep.get("vertices", []),
    )
    if not breps:
        nxt = max([0] + [int(x["id"]) for x in list(faces) + list(edges) + list(verts)])
        breps = [
            {
                "id": nxt + 1,
                "name": "",
                "faces": list(range(len(faces))),
                "edges": list(range(len(edges))),
                "vertices": list(range(len(verts))),
                "extra": "",
            }
        ]
    out = []
    for bi, b in enumerate(breps):
        bw = f"brep {bi}"
        bo = _JO()
        _add_identity(bo, b["id"], b["name"])
        fl = []
        for fi in b["faces"]:
            fi = int(fi)
            if fi < 0 or fi >= len(faces):
                _wfail(f"{bw} names face {fi}, which does not exist")
            fl.append(_face_json(brep, faces[fi], f"face {fi}"))
        bo.append(("faces", fl))
        el = []
        for ei in b["edges"]:
            ei = int(ei)
            if ei < 0 or ei >= len(edges):
                _wfail(f"{bw} names edge {ei}, which does not exist")
            e = edges[ei]
            ew = f"edge {ei}"
            eo = _JO()
            _add_identity(eo, e["id"], e["name"])
            if e["curve"] is not None:
                eo.append(("3d_curve", _curve_json(e["curve"], f"{ew} 3d_curve")))
            eo.append(("topology", _topology_json(e["topology"], ew)))
            _add_extra(eo, e["extra"], ew)
            el.append(eo)
        bo.append(("edges", el))
        vl = []
        for vi in b["vertices"]:
            vi = int(vi)
            if vi < 0 or vi >= len(verts):
                _wfail(f"{bw} names vertex {vi}, which does not exist")
            v = verts[vi]
            vw = f"vertex {vi}"
            vo = _JO()
            _add_identity(vo, v["id"], v["name"])
            vo.append(("topology", _topology_json(v["topology"], vw)))
            _add_extra(vo, v["extra"], vw)
            vl.append(vo)
        bo.append(("vertices", vl))
        _add_extra(bo, b.get("extra", ""), bw)
        out.append(bo)
    root.append(("breps", out))
    return _json_text.dumps(root, indent=4) + "\n"


def from_mesh(mesh):
    """A plain mesh as a ``brep`` dict: quads as bilinear faces, lines as edges."""
    brep = _empty()
    brep["extra"] = '{"version_number":1}'
    pts = np.asarray(mesh.points, dtype=np.float64)
    if pts.ndim != 2:
        _wfail("points are not floating point")
    dim = pts.shape[1]
    npts = len(pts)
    next_id = 1
    dropped = 0
    entry = {
        "id": -1,
        "name": "",
        "faces": [],
        "edges": [],
        "vertices": [],
        "extra": "",
    }
    for cb in mesh.cells:
        quad, line = cb.type == "quad", cb.type == "line"
        data = cb.data
        if not (quad or line) or not isinstance(data, np.ndarray):
            dropped += len(data)
            continue
        for row in np.asarray(data, dtype=np.int64):
            nodes = [int(v) for v in row]
            for v in nodes:
                if v < 0 or v >= npts:
                    _wfail(f"a cell names point {v}, which does not exist")
            order = [nodes[0], nodes[1], nodes[3], nodes[2]] if quad else nodes
            xyz = np.zeros((len(order), 3))
            for j, v in enumerate(order):
                xyz[j, : min(dim, 3)] = pts[v, : min(dim, 3)]
            w = np.ones(len(order))
            ids = np.array([v + 1 for v in order], dtype=np.int64)
            if quad:
                entry["faces"].append(len(brep["faces"]))
                brep["faces"].append(
                    {
                        "id": next_id,
                        "name": "",
                        "patch": len(brep["patches"]),
                        "swapped_normal": False,
                        "trimmed": False,
                        "rational": False,
                        "reduced_knots": False,
                        "control_point_ids": ids,
                        "loops": [],
                        "embedded_edges": [],
                        "surface_extra": "",
                        "extra": "",
                    }
                )
                brep["patches"].append(
                    _patch([1, 1], [[0.0, 0.0, 1.0, 1.0], [0.0, 0.0, 1.0, 1.0]], xyz, w)
                )
            else:
                entry["edges"].append(len(brep["edges"]))
                brep["edges"].append(
                    {
                        "id": next_id,
                        "name": "",
                        "curve": {
                            "degree": 1,
                            "knots": np.array([0.0, 0.0, 1.0, 1.0]),
                            "reduced_knots": False,
                            "rational": False,
                            "control_points": xyz,
                            "weights": w,
                            "ids": ids,
                            "active_range": [],
                            "extra": "",
                        },
                        "topology": [],
                        "extra": "",
                    }
                )
            next_id += 1
    if dropped:
        warn(
            f"cad_json: {dropped} cells that are neither quad nor line were not written"
        )
    if not brep["faces"] and not brep["edges"]:
        _wfail("the mesh has no quad or line cells to write as faces or edges")
    entry["id"] = next_id
    brep["breps"].append(entry)
    return brep


def write(filename, mesh, brep=None):
    """Write CAD JSON from ``brep`` (default ``mesh.info["brep"]``), else from the mesh."""
    if brep is None:
        brep = (
            (mesh.info or {}).get("brep")
            if isinstance(getattr(mesh, "info", None), dict)
            else None
        )
    if not brep or not (brep.get("faces") or brep.get("edges") or brep.get("vertices")):
        brep = from_mesh(mesh)
    text = render(brep)
    with open(filename, "wb") as f:
        f.write(text.encode("utf-8"))

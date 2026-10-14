"""The isogeometric layer: evaluate, refine and sample NURBS models.

The numpy twin of ``operations/iga.hpp``. A model is what a reader keeps beside
its sampled mesh: ``mesh.info["nurbs"]`` (MFEM NURBS meshes, the C++
``NurbsInfo``) or ``mesh.info["brep"]`` (Kratos CAD JSON, STEP, IGES and
:func:`fit_cad`, the C++ ``BrepInfo``). The operations act on the model, not on
its sampling:

* :func:`iga_evaluate` -- points, first derivatives, unit normals and, for a
  :class:`NurbsInfo`-style model, its NURBS fields at parametric coordinates;
  for a trimmed face also which points lie inside its trimming loops;
* :func:`iga_refine` -- exact knot insertion (h) and degree elevation (p) of
  every patch, k-refinement being the two in that order, returned with the
  refined model sampled again;
* :func:`iga_sample` -- the model as Lagrange cells.

The twin repeats the core's arithmetic order (control nets are held
homogeneous, knot insertion is The NURBS Book's A5.4, degree elevation A5.9,
derivatives A2.3) on numpy arrays that vectorise over lines and points without
reordering any single operation, so both engines give the same bits. See
``doc/iga.md``.
"""

from __future__ import annotations

import copy
import math

import numpy as np

from . import _nurbs_eval as ne
from ._exceptions import WriteError
from ._fallback import core_op_declined
from ._lagrange import vtk_lattice
from ._mesh import Mesh
from .cad_json._cad_json import _axis, _inside, _loop_polygon

_PREFIX = "meshio++: iga: "
_MAX_DEGREE = 32
_MAX_SUBDIVISIONS = 64


def _fail(what):
    raise ValueError(_PREFIX + what)


# ---------------------------------------------------------------------------
# Models and patches
# ---------------------------------------------------------------------------


def _complete_nurbs(data):
    """A ``nurbs`` model dict with the optional parts a hand-built one may
    leave out filled in (no fields, no topology, unit weights, no cell map)."""
    out = dict(data)
    out.setdefault("field_names", [])
    out.setdefault("topology", None)
    out.setdefault("cell_patch", np.zeros(0, dtype=np.int64))
    out.setdefault("cell_span", np.zeros((0, 3), dtype=np.int64))
    out.setdefault("boundary_cell_patch", np.zeros(0, dtype=np.int64))
    patches = []
    for p in data["patches"]:
        q = dict(p)
        q["degrees"] = np.asarray(p["degrees"], dtype=np.int32)
        q["knots"] = [np.asarray(k, dtype=np.float64) for k in p["knots"]]
        q["control_points"] = np.asarray(p["control_points"], dtype=np.float64)
        n = q["control_points"].shape[0]
        q["weights"] = np.asarray(p.get("weights", np.ones(n)), dtype=np.float64)
        q["fields"] = [np.asarray(f, dtype=np.float64) for f in p.get("fields", [])]
        patches.append(q)
    out["patches"] = patches
    return out


def _model(mesh, model):
    """``(kind, dict)`` of the model a mesh (or a bare model dict) carries."""
    kind, data = _find_model(mesh, model)
    return kind, (_complete_nurbs(data) if kind == "nurbs" else data)


def _find_model(mesh, model):
    if model not in ("auto", "nurbs", "brep"):
        _fail(f"model must be 'auto', 'nurbs' or 'brep', not {model!r}")
    if isinstance(mesh, dict):
        kind = "brep" if "faces" in mesh else "nurbs"
        if model not in ("auto", kind):
            _fail(f"the model is a {kind!r} model, not {model!r}")
        return kind, mesh
    info = getattr(mesh, "info", None) or {}
    have = [k for k in ("nurbs", "brep") if info.get(k) is not None]
    if model == "auto":
        if not have:
            _fail("the mesh carries no NURBS model (mesh.info['nurbs'] or ['brep'])")
        if len(have) > 1:
            _fail("the mesh carries both a 'nurbs' and a 'brep' model; name one")
        model = have[0]
    elif model not in have:
        _fail(f"the mesh carries no {model!r} model")
    return model, info[model]


def _check_patch(patch, what):
    degrees = [int(d) for d in patch["degrees"]]
    dim = len(degrees)
    knots = [np.asarray(k, dtype=np.float64) for k in patch["knots"]]
    if dim < 1 or dim > 3 or len(knots) != dim:
        _fail(f"{what} needs 1 to 3 degrees and as many knot vectors")
    counts = []
    for d in range(dim):
        p = degrees[d]
        if p < 1 or p > _MAX_DEGREE:
            _fail(f"{what} degree {p} is outside 1..{_MAX_DEGREE}")
        k = knots[d]
        if len(k) < 2 * p + 2:
            _fail(f"{what} knot vector is too short")
        if not np.all(np.isfinite(k)) or np.any(k[1:] < k[:-1]):
            _fail(f"{what} knot vector is not finite and non-decreasing")
        counts.append(len(k) - p - 1)
    n = 1
    for c in counts:
        n *= c
    cps = np.asarray(patch["control_points"])
    if cps.dtype != np.float64 or cps.ndim != 2 or cps.shape[0] != n:
        _fail(f"{what} control points are not Float64 ({n}, 1..3) for its knot vectors")
    if not 1 <= cps.shape[1] <= 3:
        _fail(f"{what} control points are not Float64 ({n}, 1..3) for its knot vectors")
    w = np.asarray(patch["weights"])
    if w.dtype != np.float64 or w.shape != (n,):
        _fail(f"{what} has {w.size} weights for {n} control points")
    comps = []
    for f in patch.get("fields", []):
        f = np.asarray(f)
        if f.dtype != np.float64 or f.ndim != 2 or f.shape[0] != n or f.shape[1] < 1:
            _fail(f"{what} has a field that is not Float64 ({n}, components)")
        comps.append(f.shape[1])
    return {
        "dim": dim,
        "sdim": cps.shape[1],
        "degrees": degrees,
        "knots": knots,
        "counts": counts + [1] * (3 - dim),
        "n": n,
        "comps": comps,
    }


def _is_rational(patch):
    return bool(np.any(np.asarray(patch["weights"]) != 1.0))


def _check_weights(patch, what):
    w = np.asarray(patch["weights"])
    if not np.all(np.isfinite(w)) or np.any(w == 0.0):
        _fail(f"{what} has a zero or non-finite weight")


def _make_net(patch, sh, rational):
    """The homogeneous net: ``(w x, w fields..., w)``, or ``(x, fields...)``."""
    w = np.asarray(patch["weights"], dtype=np.float64)
    cols = [np.asarray(patch["control_points"], dtype=np.float64)]
    cols += [np.asarray(f, dtype=np.float64) for f in patch.get("fields", [])]
    if rational:
        cols = [w[:, None] * c for c in cols] + [w[:, None]]
    return {
        "dim": sh["dim"],
        "degrees": list(sh["degrees"]),
        "knots": [k.copy() for k in sh["knots"]],
        "counts": list(sh["counts"]),
        "values": np.ascontiguousarray(np.concatenate(cols, axis=1)),
    }


def _patch_from_net(net, sh, rational):
    values = net["values"]
    n = values.shape[0]
    out = {
        "degrees": np.array(net["degrees"][: sh["dim"]], dtype=np.int32),
        "knots": [np.array(k, dtype=np.float64) for k in net["knots"][: sh["dim"]]],
    }
    w = values[:, -1] if rational else np.ones(n, dtype=np.float64)
    parts = []
    c = 0
    for width in [sh["sdim"]] + sh["comps"]:
        block = values[:, c : c + width]
        parts.append(np.ascontiguousarray(block / w[:, None] if rational else block))
        c += width
    out["control_points"] = parts[0]
    out["weights"] = np.ascontiguousarray(w)
    out["fields"] = parts[1:]
    return out


# ---------------------------------------------------------------------------
# Knot insertion (A5.4) and degree elevation (A5.9) on lines of a net
# ---------------------------------------------------------------------------


def _split_knots(p, knots, s):
    x = []
    if s <= 1:
        return x
    n = len(knots) - p - 1
    for i in ne.spans(p, knots, n):
        k0, k1 = float(knots[i]), float(knots[i + 1])
        for j in range(1, s):
            x.append(k0 + (k1 - k0) * (j / s))
    return x


def _refine_line(p, u, pw, x):
    """A5.4 on ``pw`` (``(n, T)`` rows) with the ascending knots ``x``."""
    n = pw.shape[0] - 1
    m = n + p + 1
    r = len(x) - 1
    a = ne.find_span(p, u, n + 1, x[0])
    b = ne.find_span(p, u, n + 1, x[r]) + 1
    ubar = [0.0] * (m + r + 2)
    qw = np.zeros((pw.shape[0] + r + 1, pw.shape[1]))
    for j in range(0, a - p + 1):
        qw[j] = pw[j]
    for j in range(b - 1, n + 1):
        qw[j + r + 1] = pw[j]
    for j in range(0, a + 1):
        ubar[j] = u[j]
    for j in range(b + p, m + 1):
        ubar[j + r + 1] = u[j]
    i = b + p - 1
    k = b + p + r
    for j in range(r, -1, -1):
        xj = x[j]
        while xj <= u[i] and i > a:
            qw[k - p - 1] = pw[i - p - 1]
            ubar[k] = u[i]
            k -= 1
            i -= 1
        qw[k - p - 1] = qw[k - p]
        for ll in range(1, p + 1):
            ind = k - p + ll
            alfa = ubar[k + ll] - xj
            if abs(alfa) == 0.0:
                qw[ind - 1] = qw[ind]
            else:
                alfa = alfa / (ubar[k + ll] - u[i - p + ll])
                qw[ind - 1] = alfa * qw[ind - 1] + (1.0 - alfa) * qw[ind]
        ubar[k] = xj
        k -= 1
    return ubar, qw


def _elevate_line(p, u, pw, t):
    """A5.9 on ``pw`` (``(n, T)`` rows), raising the degree by ``t``."""
    n = pw.shape[0] - 1
    m = n + p + 1
    ph = p + t
    ph2 = ph // 2
    distinct = sum(1 for i in range(p + 1, n + 1) if u[i] != u[i - 1])
    max_rows = n + 1 + t * (distinct + 1)
    width = pw.shape[1]
    qw = np.zeros((max_rows, width))
    uh = [0.0] * (max_rows + ph + 1)
    bezalfs = [[0.0] * (p + 1) for _ in range(ph + 1)]
    bezalfs[0][0] = 1.0
    bezalfs[ph][p] = 1.0
    for i in range(1, ph2 + 1):
        inv = 1.0 / float(math.comb(ph, i))
        for j in range(max(0, i - t), min(p, i) + 1):
            bezalfs[i][j] = inv * float(math.comb(p, j)) * float(math.comb(t, i - j))
    for i in range(ph2 + 1, ph):
        for j in range(max(0, i - t), min(p, i) + 1):
            bezalfs[i][j] = bezalfs[ph - i][p - j]
    bpts = np.zeros((p + 1, width))
    ebpts = np.zeros((ph + 1, width))
    nextbpts = np.zeros((max(p - 1, 1), width))
    alfs = [0.0] * max(p - 1, 1)
    kind = ph + 1
    r = -1
    a = p
    b = p + 1
    cind = 1
    ua = u[0]
    qw[0] = pw[0]
    for i in range(ph + 1):
        uh[i] = ua
    for i in range(p + 1):
        bpts[i] = pw[i]
    while b < m:
        i = b
        while b < m and u[b] == u[b + 1]:
            b += 1
        mul = b - i + 1
        ub = u[b]
        oldr = r
        r = p - mul
        lbz = (oldr + 2) // 2 if oldr > 0 else 1
        rbz = ph - (r + 1) // 2 if r > 0 else ph
        if r > 0:
            numer = ub - ua
            for k in range(p, mul, -1):
                alfs[k - mul - 1] = numer / (u[a + k] - ua)
            for j in range(1, r + 1):
                save = r - j
                s = mul + j
                for k in range(p, s - 1, -1):
                    al = alfs[k - s]
                    bpts[k] = al * bpts[k] + (1.0 - al) * bpts[k - 1]
                nextbpts[save] = bpts[p]
        for k in range(lbz, ph + 1):
            ebpts[k] = 0.0
            for j in range(max(0, k - t), min(p, k) + 1):
                ebpts[k] = ebpts[k] + bezalfs[k][j] * bpts[j]
        if oldr > 1:
            first = kind - 2
            last = kind
            den = ub - ua
            bet = (ub - uh[kind - 1]) / den
            for tr in range(1, oldr):
                ii = first
                j = last
                kj = j - kind + 1
                while j - ii > tr:
                    if ii < cind:
                        alf = (ub - uh[ii]) / (ua - uh[ii])
                        qw[ii] = alf * qw[ii] + (1.0 - alf) * qw[ii - 1]
                    if j >= lbz:
                        if j - tr <= kind - ph + oldr:
                            gam = (ub - uh[j - tr]) / den
                            ebpts[kj] = gam * ebpts[kj] + (1.0 - gam) * ebpts[kj + 1]
                        else:
                            ebpts[kj] = bet * ebpts[kj] + (1.0 - bet) * ebpts[kj + 1]
                    ii += 1
                    j -= 1
                    kj -= 1
                first -= 1
                last += 1
        if a != p:
            for _ in range(ph - oldr):
                uh[kind] = ua
                kind += 1
        for j in range(lbz, rbz + 1):
            qw[cind] = ebpts[j]
            cind += 1
        if b < m:
            for j in range(r):
                bpts[j] = nextbpts[j]
            for j in range(r, p + 1):
                bpts[j] = pw[b - p + j]
            a = b
            b += 1
            ua = ub
        else:
            for k in range(ph + 1):
                uh[kind + k] = ub
    return uh[: cind + ph + 1], qw[:cind]


def _apply(net, d, fn):
    """Apply ``fn(rows (n, T)) -> (knots, rows')`` along direction ``d``."""
    c0, c1, c2 = net["counts"]
    width = net["values"].shape[1]
    arr = net["values"].reshape(c2, c1, c0, width)
    axis = 2 - d
    moved = np.moveaxis(arr, axis, 0)
    rest = moved.shape[1:]
    rows = np.ascontiguousarray(moved.reshape(moved.shape[0], -1))
    knots, out = fn(rows)
    out = out.reshape((out.shape[0],) + rest)
    arr2 = np.ascontiguousarray(np.moveaxis(out, 0, axis))
    net["values"] = np.ascontiguousarray(arr2.reshape(-1, width))
    net["counts"][d] = out.shape[0]
    return knots


def _net_insert(net, d, x):
    if not x:
        return
    p = net["degrees"][d]
    n = net["counts"][d]
    lo, hi = net["knots"][d][p], net["knots"][d][n]
    for i, v in enumerate(x):
        if not (lo < v < hi) or (i and v < x[i - 1]):
            _fail("knots to insert must be ascending and inside the parameter range")
    u = [float(k) for k in net["knots"][d]]
    net["knots"][d] = np.array(
        _apply(net, d, lambda rows: _refine_line(p, u, rows, x)), dtype=np.float64
    )


def _net_elevate(net, d, t):
    if t <= 0:
        return
    p = net["degrees"][d]
    u = [float(k) for k in net["knots"][d]]
    net["knots"][d] = np.array(
        _apply(net, d, lambda rows: _elevate_line(p, u, rows, t)), dtype=np.float64
    )
    net["degrees"][d] = p + t


def _refine_patch(patch, rational, subdivisions, elevate, what):
    sh = _check_patch(patch, what)
    if rational:
        _check_weights(patch, what)
    # Nothing to do: the patch itself, not its homogeneous round trip (which
    # can move a rational control point by an ulp).
    identity = all(elevate[d] == 0 and subdivisions[d] == 1 for d in range(sh["dim"]))
    if identity and (rational or not _is_rational(patch)):
        return {
            "degrees": np.array(sh["degrees"], dtype=np.int32),
            "knots": [k.copy() for k in sh["knots"]],
            "control_points": np.array(patch["control_points"], dtype=np.float64),
            "weights": np.array(patch["weights"], dtype=np.float64),
            "fields": [np.array(f, dtype=np.float64) for f in patch.get("fields", [])],
        }
    net = _make_net(patch, sh, rational)
    for d in range(sh["dim"]):
        if elevate[d] < 0:
            _fail(f"cannot elevate a degree by {elevate[d]}")
        if net["degrees"][d] + elevate[d] > _MAX_DEGREE:
            _fail(f"{what}: the degree would pass {_MAX_DEGREE}")
        _net_elevate(net, d, elevate[d])
    for d in range(sh["dim"]):
        if not 1 <= subdivisions[d] <= _MAX_SUBDIVISIONS:
            _fail(f"subdivisions must be 1 to {_MAX_SUBDIVISIONS}")
        x = _split_knots(net["degrees"][d], net["knots"][d], subdivisions[d])
        _net_insert(net, d, x)
    return _patch_from_net(net, sh, rational)


# ---------------------------------------------------------------------------
# Evaluation (A2.3, rational quotient rule)
# ---------------------------------------------------------------------------


def _basis_ders(span, u, p, nder, knots):
    """``(nder + 1, N, p + 1)``: the basis functions and, for ``nder == 1``,
    their derivatives, at ``u`` in knot spans ``span`` (arrays of length N)."""
    n = len(u)
    ndu = np.zeros((p + 1, p + 1, n))
    left = np.zeros((p + 1, n))
    right = np.zeros((p + 1, n))
    ndu[0, 0] = 1.0
    for j in range(1, p + 1):
        left[j] = u - knots[span + 1 - j]
        right[j] = knots[span + j] - u
        saved = np.zeros(n)
        for r in range(j):
            ndu[j, r] = right[r + 1] + left[j - r]
            temp = ndu[r, j - 1] / ndu[j, r]
            ndu[r, j] = saved + right[r + 1] * temp
            saved = left[j - r] * temp
        ndu[j, j] = saved
    out = np.zeros((nder + 1, n, p + 1))
    for j in range(p + 1):
        out[0, :, j] = ndu[j, p]
    if nder >= 1:
        for r in range(p + 1):
            d = np.zeros(n)
            if r >= 1:
                d = (1.0 / ndu[p, r - 1]) * ndu[r - 1, p - 1]
            if r + 1 <= p:
                d = d + (-1.0 / ndu[p, r]) * ndu[r, p - 1]
            out[1, :, r] = d
        factor = float(p)
        for j in range(p + 1):
            out[1, :, j] = out[1, :, j] * factor
    return out


def _eval_batch(net, spans, u, nder):
    """Values ``(N, NV)`` (and derivatives ``(dim, N, NV)``) at ``u``
    ``(N, dim)`` in knot spans ``spans`` ``(N, dim)``."""
    dim = net["dim"]
    values, weights = net["values"], net["weights"]
    nv = values.shape[1]
    n = u.shape[0]
    ders = [
        _basis_ders(spans[:, d], u[:, d], net["degrees"][d], nder, net["knots"][d])
        for d in range(dim)
    ]
    np1 = [net["degrees"][d] + 1 if d < dim else 1 for d in range(3)]
    base = [spans[:, d] - net["degrees"][d] if d < dim else 0 for d in range(3)]
    counts = net["counts"]
    acc = np.zeros((n, nv))
    wsum = np.zeros(n)
    dacc = np.zeros((dim, n, nv))
    dwsum = np.zeros((dim, n))
    for j2 in range(np1[2]):
        for j1 in range(np1[1]):
            for j0 in range(np1[0]):
                j = (j0, j1, j2)
                idx = (base[0] + j0) + counts[0] * (
                    (base[1] + j1) + counts[1] * (base[2] + j2)
                )
                w = weights[idx] if weights is not None else 1.0
                v = values[idx]
                b = np.ones(n)
                for d in range(dim):
                    b = b * ders[d][0][:, j[d]]
                cc = w * b
                acc = acc + cc[:, None] * v
                wsum = wsum + cc
                if nder >= 1:
                    for k in range(dim):
                        bk = np.ones(n)
                        for d in range(dim):
                            bk = bk * (
                                ders[d][1][:, j[d]] if d == k else ders[d][0][:, j[d]]
                            )
                        ck = w * bk
                        dacc[k] = dacc[k] + ck[:, None] * v
                        dwsum[k] = dwsum[k] + ck
    if weights is not None:
        out = acc / wsum[:, None]
        der = None
        if nder >= 1:
            der = (dacc - dwsum[:, :, None] * out[None]) / wsum[None, :, None]
        return out, der
    return acc, (dacc if nder >= 1 else None)


def _eval_table(patch, sh, rational):
    cols = [np.asarray(patch["control_points"], dtype=np.float64)]
    cols += [np.asarray(f, dtype=np.float64) for f in patch.get("fields", [])]
    return {
        "dim": sh["dim"],
        "degrees": list(sh["degrees"]),
        "knots": sh["knots"],
        "counts": list(sh["counts"]),
        "values": np.ascontiguousarray(np.concatenate(cols, axis=1)),
        "weights": np.asarray(patch["weights"], dtype=np.float64) if rational else None,
    }


def _evaluate_patch(patch, rational, params, derivatives, what):
    if derivatives not in (0, 1):
        _fail("derivatives must be 0 or 1")
    sh = _check_patch(patch, what)
    if rational:
        _check_weights(patch, what)
    dim, sdim = sh["dim"], sh["sdim"]
    params = np.asarray(params, dtype=np.float64)
    if params.ndim == 1 and dim == 1:
        params = params[:, None]
    if params.ndim != 2 or params.shape[1] != dim:
        _fail(f"parameters are not Float64 (N, {dim})")
    n = params.shape[0]
    lo = [float(sh["knots"][d][sh["degrees"][d]]) for d in range(dim)]
    hi = [float(sh["knots"][d][sh["counts"][d]]) for d in range(dim)]
    tol = [1e-12 * (hi[d] - lo[d]) for d in range(dim)]
    for i in range(n):
        for d in range(dim):
            u = params[i, d]
            if not (lo[d] - tol[d] <= u <= hi[d] + tol[d]):
                _fail(
                    f"parameter {i}, direction {d} is outside the knot range "
                    f"[{lo[d]:f}, {hi[d]:f}]"
                )
    u = np.empty((n, dim))
    spans = np.empty((n, dim), dtype=np.int64)
    for d in range(dim):
        uc = np.minimum(np.maximum(params[:, d], lo[d]), hi[d])
        u[:, d] = uc
        spans[:, d] = [
            ne.find_span(sh["degrees"][d], sh["knots"][d], sh["counts"][d], float(x))
            for x in uc
        ]
    table = _eval_table(patch, sh, rational)
    values, ders = _eval_batch(table, spans, u, derivatives)
    out = {"points": np.ascontiguousarray(values[:, :sdim])}
    fields = []
    c = sdim
    for width in sh["comps"]:
        fields.append(np.ascontiguousarray(values[:, c : c + width]))
        c += width
    out["fields"] = fields
    if derivatives:
        d_arr = np.ascontiguousarray(np.moveaxis(ders[:, :, :sdim], 0, 1))
        out["derivatives"] = d_arr
        if dim == 2 and sdim == 3:
            a, b = d_arr[:, 0], d_arr[:, 1]
            cx = a[:, 1] * b[:, 2] - a[:, 2] * b[:, 1]
            cy = a[:, 2] * b[:, 0] - a[:, 0] * b[:, 2]
            cz = a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0]
            ln = np.sqrt(cx * cx + cy * cy + cz * cz)
            safe = np.where(ln > 0.0, ln, 1.0)
            nrm = np.stack([cx / safe, cy / safe, cz / safe], axis=1)
            nrm[ln <= 0.0] = 0.0
            out["normals"] = nrm
    return out


# ---------------------------------------------------------------------------
# Sampling NURBS patches into Lagrange cells
# ---------------------------------------------------------------------------


def _family(dim, order):
    names = {
        1: ("line", "line3", "VTK_LAGRANGE_CURVE"),
        2: ("quad", "quad9", "VTK_LAGRANGE_QUADRILATERAL"),
        3: ("hexahedron", "hexahedron27", "VTK_LAGRANGE_HEXAHEDRON"),
    }[dim]
    return names[0] if order <= 1 else names[1] if order == 2 else names[2]


def _lattice_axis(lo, hi, k):
    n = len(lo)
    g = [0.0] * (n * k + 1)
    for c in range(n):
        for i in range(k):
            g[c * k + i] = lo[c] + (hi[c] - lo[c]) * (i / k)
    if n:
        g[n * k] = hi[n - 1]
    return g


def _sample_nurbs(model, subdivisions):
    sub = int(subdivisions)
    if sub < 0 or sub > _MAX_SUBDIVISIONS:
        _fail(f"subdivisions must be 0 (automatic) to {_MAX_SUBDIVISIONS}")
    patches = model["patches"]
    if not patches:
        _fail("the model has no patches")
    topo = model.get("topology")
    if sub == 0 and topo is not None and int(topo["dimension"]) != 0:
        # A model with its MFEM patch topology samples as the MFEM reader does.
        from .mfem import _mfem

        try:
            mesh = _mfem.sample_nurbs(model)
        except WriteError as exc:
            _fail(f"the patch topology is inconsistent: {exc}")
        # The model stays as it is (the sampler regathers its control points
        # from MFEM's numbering, which can move a shared one by an ulp); only
        # the cell maps are the sampling's.
        derived = mesh.info["nurbs"]
        mesh.info = None
        mapped = dict(model)
        for key in ("cell_patch", "cell_span", "boundary_cell_patch"):
            mapped[key] = derived[key]
        return mesh, mapped
    s = max(sub, 1)
    shapes, tables = [], []
    dim = sdim = None
    k = 1
    names = list(model["field_names"])
    for i, p in enumerate(patches):
        what = f"patch {i}"
        rational = _is_rational(p)
        sh = _check_patch(p, what)
        if rational:
            _check_weights(p, what)
        if dim is None:
            dim, sdim = sh["dim"], sh["sdim"]
        if sh["dim"] != dim or sh["sdim"] != sdim:
            _fail(f"{what} differs in parametric or space dimension from patch 0")
        if len(p.get("fields", [])) != len(names):
            _fail(
                f"{what} has {len(p.get('fields', []))} fields for {len(names)} names"
            )
        if shapes and sh["comps"] != shapes[0]["comps"]:
            _fail(f"{what} has fields of different components than patch 0")
        for d in range(dim):
            k = max(k, sh["degrees"][d])
        shapes.append(sh)
        tables.append(_eval_table(p, sh, rational))
    shape = {1: "line", 2: "quad", 3: "hexahedron"}[dim]
    lattice = vtk_lattice(shape, k)
    comps = shapes[0]["comps"]
    cols = 2 if sdim == 2 else 3
    conn, cell_patch, cell_span = [], [], []
    point_blocks = []
    num_points = 0
    for pi, p in enumerate(patches):
        sh = shapes[pi]
        axis = []
        for d in range(dim):
            axis.append(
                _axis(sh["degrees"][d], sh["knots"][d], sh["counts"][d], s, None)
            )
        ncell = [len(a[0]) for a in axis] + [1] * (3 - dim)
        g = [_lattice_axis(a[0], a[1], k) for a in axis]
        width = [len(x) for x in g] + [1] * (3 - dim)
        spans = []
        for d in range(dim):
            sp = axis[d][2]
            spans.append(
                [sp[min(i // k, len(sp) - 1)] if sp else 0 for i in range(len(g[d]))]
            )
        local = {}
        new_points = []
        for iw in range(ncell[2]):
            for jv in range(ncell[1]):
                for iu in range(ncell[0]):
                    for i, j, ll in lattice:
                        gi = iu * k + i
                        gj = jv * k + j if dim > 1 else 0
                        gl = iw * k + ll if dim > 2 else 0
                        gg = gi + width[0] * (gj + width[1] * gl)
                        if gg not in local:
                            local[gg] = num_points + len(new_points)
                            new_points.append(gg)
                        conn.append(local[gg])
                    cell_patch.append(pi)
                    cell_span.append(
                        (
                            axis[0][3][iu],
                            axis[1][3][jv] if dim > 1 else 0,
                            axis[2][3][iw] if dim > 2 else 0,
                        )
                    )
        if new_points:
            idx = np.array(new_points, dtype=np.int64)
            at = [
                idx % width[0],
                (idx // width[0]) % width[1],
                idx // (width[0] * width[1]),
            ]
            u = np.empty((len(idx), dim))
            sp = np.empty((len(idx), dim), dtype=np.int64)
            for d in range(dim):
                u[:, d] = np.array(g[d])[at[d]]
                sp[:, d] = np.array(spans[d], dtype=np.int64)[at[d]]
            values, _ = _eval_batch(tables[pi], sp, u, 0)
            point_blocks.append(values)
        num_points += len(new_points)
    nv = sdim + sum(comps)
    values = np.concatenate(point_blocks) if point_blocks else np.zeros((0, nv))
    points = np.zeros((num_points, cols))
    points[:, : min(cols, sdim)] = values[:, : min(cols, sdim)]
    point_data = {}
    c = sdim
    for name, width_f in zip(names, comps):
        a = values[:, c : c + width_f]
        point_data[name] = np.ascontiguousarray(a[:, 0] if width_f == 1 else a)
        c += width_f
    mapped = dict(model)
    n_cells = len(cell_patch)
    cp = np.array(cell_patch, dtype=np.int64)
    cs = np.array(cell_span, dtype=np.int64).reshape(n_cells, 3)
    mapped["cell_patch"] = cp
    mapped["cell_span"] = cs
    mapped["boundary_cell_patch"] = np.zeros(0, dtype=np.int64)
    if n_cells:
        cells = [
            (
                _family(dim, k),
                np.array(conn, dtype=np.int64).reshape(n_cells, len(lattice)),
            )
        ]
        cell_data = {"iga:patch": [cp], "iga:span": [cs]}
    else:
        cells, cell_data = [], {}
    return Mesh(points, cells, point_data=point_data, cell_data=cell_data), mapped


# ---------------------------------------------------------------------------
# Public operations
# ---------------------------------------------------------------------------


def _triple(value, name, default):
    """An int or a sequence of up to 3 ints, as a 3-list (unused: ``default``)."""
    if isinstance(value, (int, np.integer)):
        out = [int(value)] * 3
    else:
        seq = [int(v) for v in value]
        if not 1 <= len(seq) <= 3:
            _fail(f"{name} must be an int or 1 to 3 ints")
        out = seq + [default] * (3 - len(seq))
    return out


def _with_info(mesh, source, key, model):
    keep = dict(getattr(source, "info", None) or {})
    keep.pop("nurbs", None)
    keep.pop("brep", None)
    keep[key] = model
    mesh.info = keep
    return mesh


def _empty_brep_ids():
    return np.zeros(0, dtype=np.int64)


def _refine_brep(model, subdivisions, elevate):
    patches = model["patches"]
    flag = [-1] * len(patches)
    for fi, face in enumerate(model["faces"]):
        pi = int(face["patch"])
        if pi < 0 or pi >= len(patches):
            _fail(f"face {fi} names no patch")
        want = 1 if face["rational"] else 0
        if flag[pi] >= 0 and flag[pi] != want:
            _fail(f"patch {pi} is shared by faces that disagree on rational")
        flag[pi] = want
    out = copy.deepcopy(model)
    for i, p in enumerate(patches):
        what = f"patch {i}"
        sh = _check_patch(p, what)
        if sh["dim"] != 2:
            _fail(f"{what} is not a surface patch")
        rational = _is_rational(p) if flag[i] < 0 else flag[i] == 1
        out["patches"][i] = _refine_patch(p, rational, subdivisions, elevate, what)
    for face in out["faces"]:
        face["control_point_ids"] = _empty_brep_ids()
    return out


def _refine_nurbs(model, subdivisions, elevate):
    if not model["patches"]:
        _fail("the model has no patches")
    topo = model.get("topology")
    if topo is not None and int(topo["dimension"]) != 0:
        # The patches' directions follow their own edges, so a per-direction
        # refinement would not agree on the knot vector two patches share.
        dim = min(int(topo["dimension"]), 3)
        if any(
            subdivisions[d] != subdivisions[0] or elevate[d] != elevate[0]
            for d in range(1, dim)
        ):
            _fail(
                "a model with a patch topology refines every direction alike (its "
                "patches share knot vectors along differently oriented directions)"
            )
    out = {
        "patches": [],
        "field_names": list(model["field_names"]),
        "topology": copy.deepcopy(model["topology"]),
    }
    for i, p in enumerate(model["patches"]):
        what = f"patch {i}"
        _check_patch(p, what)
        out["patches"].append(
            _refine_patch(p, _is_rational(p), subdivisions, elevate, what)
        )
    out["cell_patch"] = np.zeros(0, dtype=np.int64)
    out["cell_span"] = np.zeros((0, 3), dtype=np.int64)
    out["boundary_cell_patch"] = np.zeros(0, dtype=np.int64)
    return out


def iga_refine(
    mesh,
    subdivisions=1,
    elevate=0,
    sample_subdivisions=0,
    model="auto",
    trim_samples=16,
    edges=True,
    trim=True,
):
    """Refine every patch of a NURBS model exactly, and sample it again.

    Degree elevation (``elevate``) comes first, then knot insertion splits
    every nonzero knot span into ``subdivisions`` equal parts, which is
    k-refinement (the highest continuity). The geometry and the
    parametrisation do not change, so a trimmed face's trimming curves stay
    exact; the NURBS fields of a ``nurbs`` model refine with the geometry.

    :param mesh: a mesh carrying ``info["nurbs"]`` or ``info["brep"]`` (or the
        model dict itself).
    :param subdivisions: parts every nonzero knot span is split into, an int
        or one per parametric direction (1 = none).
    :param elevate: degrees to add, an int or one per direction.
    :param sample_subdivisions: cells per knot span of the returned sampling
        (0 = automatic, see :func:`iga_sample`).
    :param model: ``"auto"``, ``"nurbs"`` or ``"brep"``, when the mesh carries
        both.
    :param trim_samples: ``brep``: points per knot span a trimming curve is
        sampled at for the clip.
    :param edges: ``brep``: sample the edges that have a 3-D curve.
    :param trim: ``brep``: keep only the cells inside the trimming loops.
    :returns: the refined model sampled as a :class:`Mesh`, the refined model
        in its ``info``.
    """
    kind, data = _model(mesh, model)
    sub = _triple(subdivisions, "subdivisions", 1)
    ele = _triple(elevate, "elevate", 0)
    sample_subdivisions = int(sample_subdivisions)
    out = None
    try:
        from . import _core

        if kind == "nurbs":
            out = _core.iga_refine_nurbs(data, sub, ele, sample_subdivisions)
        else:
            out = _core.iga_refine_brep(
                data,
                sub,
                ele,
                sample_subdivisions,
                int(trim_samples),
                bool(edges),
                bool(trim),
            )
    except Exception as exc:
        if not core_op_declined(exc, "iga_refine"):
            raise
        out = None
    if out is not None:
        return _with_info(out, mesh, kind, out.info[kind])
    if kind == "nurbs":
        refined = _refine_nurbs(data, sub, ele)
        sampled, mapped = _sample_nurbs(refined, sample_subdivisions)
        return _with_info(sampled, mesh, "nurbs", mapped)
    refined = _refine_brep(data, sub, ele)
    return _sample_brep_py(
        refined, mesh, sample_subdivisions, trim_samples, edges, trim
    )


def _sample_brep_py(brep, source, subdivisions, trim_samples, edges, trim):
    from .cad_json import _cad_json

    sampled = _cad_json.sample(
        brep,
        subdivisions=subdivisions,
        trim_samples=trim_samples,
        edges=edges,
        trim=trim,
    )
    return _with_info(sampled, source, "brep", brep)


def iga_sample(
    mesh,
    subdivisions=0,
    model="auto",
    trim_samples=16,
    edges=True,
    trim=True,
):
    """Sample a NURBS model into Lagrange cells.

    A ``nurbs`` model becomes one block of the highest degree's cells (line,
    quad or hexahedron family by parametric dimension), the nodes the exact
    rational geometry at each cell's lattice, the NURBS fields point data and
    ``iga:patch`` / ``iga:span`` cell data; points are shared within a patch,
    not between patches. A ``brep`` model is sampled as
    :func:`meshioplusplus.cad_json.sample` does.

    :param subdivisions: cells per knot span and direction (0 = automatic: one
        per span, about 16 per direction on a trimmed face).
    :returns: the sampled mesh with the model in its ``info``.
    """
    kind, data = _model(mesh, model)
    subdivisions = int(subdivisions)
    out = None
    try:
        from . import _core

        if kind == "nurbs":
            out = _core.iga_sample_nurbs(data, subdivisions)
        else:
            out = _core.iga_sample_brep(
                data, subdivisions, int(trim_samples), bool(edges), bool(trim)
            )
    except Exception as exc:
        if not core_op_declined(exc, "iga_sample"):
            raise
        out = None
    if out is not None:
        return _with_info(out, mesh, kind, out.info[kind])
    if kind == "nurbs":
        sampled, mapped = _sample_nurbs(data, subdivisions)
        return _with_info(sampled, mesh, "nurbs", mapped)
    return _sample_brep_py(data, mesh, subdivisions, trim_samples, edges, trim)


def _brep_inside(face, brep, params, trim_samples):
    clip = bool(face["trimmed"]) and bool(face["loops"])
    n = params.shape[0]
    if not clip:
        return np.ones(n, dtype=np.uint8)
    loops = []
    for loop in face["loops"]:
        poly = _loop_polygon(loop, trim_samples)
        if len(poly) >= 3:
            loops.append(poly)
    return np.array(
        [1 if _inside(loops, float(u), float(v)) else 0 for u, v in params],
        dtype=np.uint8,
    )


def iga_evaluate(
    mesh,
    params,
    patch=0,
    derivatives=0,
    model="auto",
    trim_samples=16,
):
    """Evaluate a patch (or a face) of a NURBS model at parametric coordinates.

    :param params: ``(N, d)`` parametric coordinates (``d`` the patch's
        parametric dimension); every one must lie in the patch's knot range.
    :param patch: the patch of a ``nurbs`` model, or the face of a ``brep``
        model, to evaluate.
    :param derivatives: 0, or 1 for the first derivatives (and the unit
        normals of a surface in 3-D).
    :param trim_samples: ``brep``: points per knot span a trimming curve is
        sampled at for the inside test.
    :returns: a dict with ``points`` ``(N, space dim)``, ``fields`` (one
        ``(N, components)`` array per NURBS field), with ``derivatives`` also
        ``derivatives`` ``(N, d, space dim)`` and, for a surface in 3-D,
        ``normals`` ``(N, 3)``, and for a ``brep`` face ``inside`` (a bool
        array: the point lies inside the face's trimming loops).
    """
    kind, data = _model(mesh, model)
    derivatives = int(derivatives)
    patch = int(patch)
    params = np.ascontiguousarray(np.asarray(params, dtype=np.float64))
    trim_samples = int(trim_samples)
    res = None
    try:
        from . import _core

        if kind == "nurbs":
            res = _core.iga_evaluate_nurbs(data, patch, params, derivatives)
        else:
            res = _core.iga_evaluate_brep(
                data, patch, params, derivatives, trim_samples
            )
    except Exception as exc:
        if not core_op_declined(exc, "iga_evaluate"):
            raise
        res = None
    if res is None:
        if kind == "nurbs":
            if not 0 <= patch < len(data["patches"]):
                _fail(f"patch {patch} is outside the model's {len(data['patches'])}")
            p = data["patches"][patch]
            res = _evaluate_patch(
                p, _is_rational(p), params, derivatives, f"patch {patch}"
            )
        else:
            res = _evaluate_brep_py(data, patch, params, derivatives, trim_samples)
    if "inside" in res:
        res["inside"] = np.asarray(res["inside"]).astype(bool)
    return res


def _evaluate_brep_py(brep, face_index, params, derivatives, trim_samples):
    if not 0 <= face_index < len(brep["faces"]):
        _fail(f"face {face_index} is outside the model's {len(brep['faces'])}")
    if not 2 <= trim_samples <= 4096:
        _fail("trim_samples must be 2 to 4096")
    face = brep["faces"][face_index]
    what = f"face {face_index}"
    pi = int(face["patch"])
    if pi < 0 or pi >= len(brep["patches"]):
        _fail(f"{what} names no patch")
    patch = brep["patches"][pi]
    if len(patch["degrees"]) != 2:
        _fail(f"{what} is not a surface patch")
    from .cad_json import _cad_json

    for loop in face["loops"]:
        for t in loop["trims"]:
            _cad_json._check_curve(t["curve"], f"{what} trimming curve")
    res = _evaluate_patch(patch, bool(face["rational"]), params, derivatives, what)
    res["inside"] = _brep_inside(
        face, brep, np.asarray(params, dtype=np.float64).reshape(-1, 2), trim_samples
    )
    return res


__all__ = ["iga_evaluate", "iga_refine", "iga_sample"]

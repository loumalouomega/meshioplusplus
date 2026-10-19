"""Fit a plane, sphere, cylinder, cone or torus to each patch of a surface.

The numpy twin of ``operations/fit_primitives.hpp``. It follows the core step
for step -- the same Jacobi eigen-decomposition, pivoted solver, frames,
initial guesses and Levenberg-Marquardt loop -- so the two agree within
rounding (the twin's sums use numpy's pairwise order).

A type is accepted when its RMS orthogonal distance is at most the tolerance
and at least ``min_inlier_fraction`` of the points lie within it; the patch
takes the first accepted type in the order plane, sphere, cylinder, cone,
torus, or ``none`` (0) with the best attempt's residuals.
"""

from __future__ import annotations

import math

import numpy as np

from ._fallback import core_op_declined
from ._feature_edges import _coords3, ring_unit_normals
from ._regions import block_bases

_PREFIX = "meshio++: fit_primitives: "
_PI = 3.14159265358979323846

#: ``fit:type`` values, in the order types are tried.
TYPES = ("none", "plane", "sphere", "cylinder", "cone", "torus")
_ALL = sum(1 << t for t in range(1, 6))


def sym3_eigen(a):
    """Cyclic Jacobi on a symmetric 3x3: ``(w ascending, rows of eigenvectors)``."""
    a = [[float(a[i][j]) for j in range(3)] for i in range(3)]
    v = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    for _sweep in range(50):
        off = a[0][1] * a[0][1] + a[0][2] * a[0][2] + a[1][2] * a[1][2]
        diag = a[0][0] * a[0][0] + a[1][1] * a[1][1] + a[2][2] * a[2][2]
        if off == 0.0 or off <= 1e-32 * diag:
            break
        for p, q in ((0, 1), (0, 2), (1, 2)):
            if a[p][q] == 0.0:
                continue
            theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q])
            t = (1.0 if theta >= 0.0 else -1.0) / (
                abs(theta) + math.sqrt(theta * theta + 1.0)
            )
            c = 1.0 / math.sqrt(t * t + 1.0)
            s = t * c
            for k in range(3):
                akp, akq = a[k][p], a[k][q]
                a[k][p] = c * akp - s * akq
                a[k][q] = s * akp + c * akq
            for k in range(3):
                apk, aqk = a[p][k], a[q][k]
                a[p][k] = c * apk - s * aqk
                a[q][k] = s * apk + c * aqk
            a[p][q] = 0.0
            a[q][p] = 0.0
            for k in range(3):
                vkp, vkq = v[k][p], v[k][q]
                v[k][p] = c * vkp - s * vkq
                v[k][q] = s * vkp + c * vkq
    order = [0, 1, 2]
    for i in range(1, 3):
        j = i
        while j > 0 and a[order[j]][order[j]] < a[order[j - 1]][order[j - 1]]:
            order[j], order[j - 1] = order[j - 1], order[j]
            j -= 1
    w, vecs = [], []
    for col in order:
        w.append(a[col][col])
        big = 0
        for k in range(1, 3):
            if abs(v[k][col]) > abs(v[big][col]):
                big = k
        sign = -1.0 if v[big][col] < 0.0 else 1.0
        vecs.append([sign * v[k][col] for k in range(3)])
    return w, vecs


def solve(a, b):
    """Gaussian elimination with partial pivoting; ``None`` on a zero pivot."""
    k = len(b)
    a = [list(map(float, row)) for row in a]
    b = [float(x) for x in b]
    for col in range(k):
        piv = col
        for r in range(col + 1, k):
            if abs(a[r][col]) > abs(a[piv][col]):
                piv = r
        pv = a[piv][col]
        if pv == 0.0 or not math.isfinite(pv):
            return None
        if piv != col:
            a[piv], a[col] = a[col], a[piv]
            b[piv], b[col] = b[col], b[piv]
        for r in range(col + 1, k):
            f = a[r][col] / pv
            if f == 0.0:
                continue
            for j in range(col, k):
                a[r][j] -= f * a[col][j]
            b[r] -= f * b[col]
    x = [0.0] * k
    for i in range(k - 1, -1, -1):
        s = b[i]
        for j in range(i + 1, k):
            s -= a[i][j] * x[j]
        x[i] = s / a[i][i]
        if not math.isfinite(x[i]):
            return None
    return x


def frame(e3):
    """``(e1, e2)`` completing the unit vector ``e3`` to a right-handed frame."""
    axis = 0
    for k in (1, 2):
        if abs(e3[k]) < abs(e3[axis]):
            axis = k
    h = [0.0, 0.0, 0.0]
    h[axis] = 1.0
    d = h[0] * e3[0] + h[1] * e3[1] + h[2] * e3[2]
    e1 = [h[k] - d * e3[k] for k in range(3)]
    ln = math.sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2])
    e1 = [x / ln for x in e1]
    e2 = [
        e3[1] * e1[2] - e3[2] * e1[1],
        e3[2] * e1[0] - e3[0] * e1[2],
        e3[0] * e1[1] - e3[1] * e1[0],
    ]
    return e1, e2


def levenberg_marquardt(theta, typical, n, resid, max_it):
    """Minimise ``sum(resid(theta)**2)``; forward-difference Jacobian."""
    theta = [float(t) for t in theta]
    k = len(theta)
    if n < k or k == 0:
        return theta
    r = resid(theta)
    cost = float(np.dot(r, r))
    if not math.isfinite(cost):
        return theta
    lam = 1e-3
    for _it in range(max_it):
        jac = np.empty((n, k))
        for j in range(k):
            t = list(theta)
            step = 1.4901161193847656e-08 * max(abs(theta[j]), typical[j])
            t[j] += step
            jac[:, j] = (resid(t) - r) / step
        g = jac.T @ r
        h = jac.T @ jac
        accepted = False
        cost_new = cost
        trial = theta
        rn = r
        for _attempt in range(12):
            m = h.copy()
            for a in range(k):
                m[a, a] += lam * (h[a, a] if h[a, a] > 1e-300 else 1e-300)
            delta = solve(m.tolist(), (-g).tolist())
            if delta is None:
                lam *= 10.0
                continue
            trial = [theta[a] + delta[a] for a in range(k)]
            rn = resid(trial)
            cost_new = float(np.dot(rn, rn))
            if math.isfinite(cost_new) and cost_new < cost:
                accepted = True
                break
            lam *= 10.0
        if not accepted:
            break
        before = cost
        theta, r, cost = trial, rn, cost_new
        lam = max(lam * 0.1, 1e-15)
        if before - cost <= 1e-14 * before:
            break
    return theta


def _stats(dist, tol):
    a = np.abs(dist)
    return (
        float(np.sqrt(np.dot(dist, dist) / len(dist))),
        float(a.max()),
        float(np.count_nonzero(a <= tol)) / len(dist),
    )


def _canonical(v):
    big = 0
    for k in (1, 2):
        if abs(v[k]) > abs(v[big]):
            big = k
    return [-x for x in v] if v[big] < 0.0 else list(v)


def _unit(v):
    ln = math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])
    return [x / ln for x in v] if ln > 0.0 else list(v)


def _scatter(v, mean):
    d = v - np.asarray(mean)
    return d.T @ d


def _axis_from(base, e1, e2, e3, u, v, w, ax, ay):
    point = [base[k] + u * e1[k] + v * e2[k] + w * e3[k] for k in range(3)]
    d = _unit([e3[k] + ax * e1[k] + ay * e2[k] for k in range(3)])
    return np.asarray(point), np.asarray(d)


def _line_distance(pts, point, d):
    v = pts - point
    h = v @ d
    r = v - h[:, None] * d[None, :]
    return np.sqrt(np.einsum("ij,ij->i", r, r))


def _axis(patch):
    normals, centres = patch["normals"], patch["centres"]
    nf = len(normals)
    if nf < 2:
        return None
    _w, v = sym3_eigen(_scatter(normals, normals.mean(axis=0)))
    d = np.asarray(v[0])
    a = np.zeros((3, 3))
    b = np.zeros(3)
    eye = np.eye(3)
    for f in range(nf):
        m = eye - np.outer(normals[f], normals[f])
        a += m
        b += m @ centres[f]
    a += float(nf) * np.outer(d, d)
    x = solve(a.tolist(), b.tolist())
    if x is None:
        return None
    return d, np.asarray(x)


def _fit_plane(patch, tol):
    pts = patch["pts"]
    if len(pts) < 3:
        return None
    _w, v = sym3_eigen(pts.T @ pts)
    nrm = np.asarray(v[0])
    if float(patch["normals"].mean(axis=0) @ nrm) < 0.0:
        nrm = -nrm
    params = np.zeros(10)
    params[3:6] = nrm
    return params, pts @ nrm


def _fit_sphere(patch, tol, max_it):
    pts = patch["pts"]
    n = len(pts)
    if n < 4:
        return None
    rows = np.column_stack([2.0 * pts, np.ones(n)])
    rhs = np.einsum("ij,ij->i", pts, pts)
    x = solve((rows.T @ rows).tolist(), (rows.T @ rhs).tolist())
    if x is None:
        return None
    r2 = x[3] + x[0] * x[0] + x[1] * x[1] + x[2] * x[2]
    if not r2 > 0.0:
        return None
    diag = patch["diag"]

    def resid(t):
        d = pts - np.asarray(t[:3])
        return np.sqrt(np.einsum("ij,ij->i", d, d)) - t[3]

    t = levenberg_marquardt(
        [x[0], x[1], x[2], math.sqrt(r2)], [diag] * 4, n, resid, max_it
    )
    if not (t[3] > 0.0 and t[3] < 1e3 * diag):
        return None
    params = np.zeros(10)
    params[0:3] = t[:3]
    params[6] = t[3]
    return params, resid(t)


def _fit_cylinder(patch, tol, max_it):
    pts = patch["pts"]
    n = len(pts)
    ax = _axis(patch) if n >= 5 else None
    if ax is None:
        return None
    d0, p0 = ax
    e1, e2 = frame(d0.tolist())
    r0 = float(np.sum(_line_distance(pts, p0, d0))) / n
    diag = patch["diag"]

    def resid(t):
        pt, d = _axis_from(p0, e1, e2, d0, t[0], t[1], 0.0, t[2], t[3])
        return _line_distance(pts, pt, d) - t[4]

    t = levenberg_marquardt(
        [0.0, 0.0, 0.0, 0.0, r0], [diag, diag, 1.0, 1.0, diag], n, resid, max_it
    )
    if not (t[4] > 0.0 and t[4] < 1e3 * diag):
        return None
    pt, d = _axis_from(p0, e1, e2, d0, t[0], t[1], 0.0, t[2], t[3])
    d = np.asarray(_canonical(d.tolist()))
    foot = pt + (-float(pt @ d)) * d
    params = np.zeros(10)
    params[0:3] = foot
    params[3:6] = d
    params[6] = t[4]
    return params, resid(t)


def _fit_cone(patch, tol, max_it):
    pts = patch["pts"]
    n = len(pts)
    ax = _axis(patch) if n >= 6 else None
    if ax is None:
        return None
    d0, p0 = ax
    h = (pts - p0) @ d0
    rho = _line_distance(pts, p0, d0)
    sh, sr = float(np.sum(h)), float(np.sum(rho))
    shh, shr = float(np.dot(h, h)), float(np.dot(h, rho))
    den = n * shh - sh * sh
    if not den > 0.0:
        return None
    slope = (n * shr - sh * sr) / den
    icpt = (sr - slope * sh) / n
    if not abs(slope) > 1e-6:
        return None
    axis = d0 if slope >= 0.0 else -d0
    apex = p0 + (-icpt / slope) * d0
    e1, e2 = frame(axis.tolist())
    diag = patch["diag"]

    def resid(t):
        a, d = _axis_from(apex, e1, e2, axis, t[0], t[1], t[2], t[3], t[4])
        v = pts - a
        hh = v @ d
        r = v - hh[:, None] * d[None, :]
        return np.sqrt(np.einsum("ij,ij->i", r, r)) * math.cos(t[5]) - hh * math.sin(
            t[5]
        )

    t = levenberg_marquardt(
        [0.0, 0.0, 0.0, 0.0, 0.0, math.atan(abs(slope))],
        [diag, diag, diag, 1.0, 1.0, 1.0],
        n,
        resid,
        max_it,
    )
    if not (t[5] > 1e-3 and t[5] < 0.5 * _PI - 1e-3):
        return None
    a, d = _axis_from(apex, e1, e2, axis, t[0], t[1], t[2], t[3], t[4])
    params = np.zeros(10)
    params[0:3] = a
    params[3:6] = d
    params[8] = t[5] * 180.0 / _PI
    return params, resid(t)


def _torus_guess(patch, minor):
    q = patch["centres"] - minor * patch["normals"]
    qm = q.mean(axis=0)
    _w, v = sym3_eigen(_scatter(q, qm))
    axis = np.asarray(v[0])
    e1, e2 = frame(v[0])
    d = q - qm
    u = d @ np.asarray(e1)
    t = d @ np.asarray(e2)
    rows = np.column_stack([2.0 * u, 2.0 * t, np.ones(len(q))])
    rhs = u * u + t * t
    x = solve((rows.T @ rows).tolist(), (rows.T @ rhs).tolist())
    if x is None:
        return None
    r2 = x[2] + x[0] * x[0] + x[1] * x[1]
    if not r2 > 0.0:
        return None
    major = math.sqrt(r2)
    radial = np.sqrt((u - x[0]) ** 2 + (t - x[1]) ** 2) - major
    off = d @ axis
    score = float(np.dot(radial, radial) + np.dot(off, off))
    centre = qm + x[0] * np.asarray(e1) + x[1] * np.asarray(e2)
    return centre, axis, major, score


def _fit_torus(patch, tol, max_it):
    pts = patch["pts"]
    n = len(pts)
    if n < 7 or len(patch["centres"]) < 3:
        return None
    diag = patch["diag"]
    best = None
    for sign in (1.0, -1.0):
        for j in range(48):
            rho = diag * math.pow(10.0, -3.0 + 3.0 * j / 47.0)
            g = _torus_guess(patch, sign * rho)
            if g is None:
                continue
            if best is None or g[3] < best[0]:
                best = (g[3], g[0], g[1], g[2], rho)
    if best is None:
        return None
    _score, centre, axis, major, minor = best
    e1, e2 = frame(axis.tolist())

    def resid(t):
        c, d = _axis_from(centre, e1, e2, axis, t[0], t[1], t[2], t[3], t[4])
        v = pts - c
        hh = v @ d
        r = v - hh[:, None] * d[None, :]
        radial = np.sqrt(np.einsum("ij,ij->i", r, r)) - t[5]
        return np.sqrt(radial * radial + hh * hh) - t[6]

    t = levenberg_marquardt(
        [0.0, 0.0, 0.0, 0.0, 0.0, major, minor],
        [diag, diag, diag, 1.0, 1.0, diag, diag],
        n,
        resid,
        max_it,
    )
    if not (t[5] > 0.0 and t[6] > 0.0 and t[5] < 1e3 * diag and t[6] < 1e3 * diag):
        return None
    c, d = _axis_from(centre, e1, e2, axis, t[0], t[1], t[2], t[3], t[4])
    params = np.zeros(10)
    params[0:3] = c
    params[3:6] = _canonical(d.tolist())
    params[6] = t[5]
    params[7] = t[6]
    return params, resid(t)


_FITTERS = (None, _fit_plane, _fit_sphere, _fit_cylinder, _fit_cone, _fit_torus)


def _primitive_distance(t, p, q):
    """Unsigned distance of one point to a fitted primitive (``fit:patch_params``)."""
    x = [float(p[0]), float(p[1]), float(p[2])]
    d = [float(p[3]), float(p[4]), float(p[5])]
    v = [float(q[k]) - x[k] for k in range(3)]
    if t == 1:
        return abs(v[0] * d[0] + v[1] * d[1] + v[2] * d[2])
    if t == 2:
        return abs(math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) - float(p[6]))
    h = v[0] * d[0] + v[1] * d[1] + v[2] * d[2]
    r = [v[k] - h * d[k] for k in range(3)]
    rho = math.sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2])
    if t == 3:
        return abs(rho - float(p[6]))
    if t == 4:
        a = float(p[8]) * math.pi / 180.0
        if h * math.cos(a) + rho * math.sin(a) < 0.0:
            return math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])
        return abs(rho * math.cos(a) - h * math.sin(a))
    radial = rho - float(p[6])
    return abs(math.sqrt(radial * radial + h * h) - float(p[7]))


def _centre_misfit(t, p, patch):
    """RMS distance of the facets' centres to a fitted primitive."""
    total = 0.0
    centres = patch["centres"]
    for c in centres:
        dd = _primitive_distance(t, p, c)
        total += dd * dd
    return math.sqrt(total / len(centres)) if len(centres) else 0.0


def _fit_patch(patch, opts):
    params = np.full(10, np.nan)
    tol = opts["tolerance"] * patch["diag"] if opts["relative"] else opts["tolerance"]
    out = {"type": 0, "params": params, "rms": np.nan, "max": np.nan}
    out.update(fraction=0.0, tol=tol)
    if len(patch["pts"]) == 0 or not patch["diag"] > 0.0:
        return out
    max_it = min(opts["max_iterations"], 100000)
    have_best = False
    for t in range(1, 6):
        if not opts["types"] & (1 << t):
            continue
        fit = _FITTERS[t](patch, tol) if t == 1 else _FITTERS[t](patch, tol, max_it)
        if fit is None:
            continue
        p, dist = fit
        rms, mx, frac = _stats(np.asarray(dist, dtype=np.float64), tol)
        if not math.isfinite(rms):
            continue
        if rms <= tol and frac >= opts["min_inlier_fraction"]:
            # A sphere through two rims of points may as well be a cylinder,
            # cone or torus: a cylinder or cone the facets' centres lie nearer
            # to wins, and a torus (with a parameter to spare) when they lie
            # under half as far from it.
            best_t = t
            if t == 2:
                best_misfit = _centre_misfit(2, p, patch)
                for u in (3, 4, 5):
                    if not opts["types"] & (1 << u):
                        continue
                    other = _FITTERS[u](patch, tol, max_it)
                    if other is None:
                        continue
                    o_rms, o_mx, o_frac = _stats(
                        np.asarray(other[1], dtype=np.float64), tol
                    )
                    if not (
                        math.isfinite(o_rms)
                        and o_rms <= tol
                        and o_frac >= opts["min_inlier_fraction"]
                    ):
                        continue
                    misfit = _centre_misfit(u, other[0], patch)
                    if misfit < (0.5 if u == 5 else 1.0) * best_misfit:
                        best_misfit = misfit
                        best_t = u
                        p, rms, mx, frac = other[0], o_rms, o_mx, o_frac
            p = np.array(p, dtype=np.float64)
            p[0:3] += patch["shift"]
            out.update(type=best_t, params=p, rms=rms, max=mx, fraction=frac)
            return out
        if not have_best or rms < out["rms"]:
            have_best = True
            out.update(rms=rms, max=mx, fraction=frac)
    return out


def _fit_primitives_py(mesh, opts):
    from ._segment import surface_rings

    bases = block_bases(mesh.cells)
    total = int(bases[-1]) + len(mesh.cells[-1].data) if len(mesh.cells) else 0
    label = np.zeros(total, dtype=np.int64)
    name = opts["patch_data"]
    if name:
        if name not in mesh.cell_data:
            raise ValueError(
                f"{_PREFIX}no cell data named '{name}' "
                "(run segment first, or name the label array)"
            )
        for base, block, arr in zip(bases, mesh.cells, mesh.cell_data[name]):
            arr = np.asarray(arr)
            if arr.size != len(block.data):
                raise ValueError(
                    f"{_PREFIX}cell data '{name}' must hold one label per cell"
                )
            if arr.dtype.kind == "f":
                raise ValueError(
                    f"{_PREFIX}cell data '{name}' must be an integer label"
                )
            label[base : base + len(block.data)] = arr.ravel().astype(np.int64)

    start, nodes, cell_of = surface_rings(
        mesh, np.ones(total, dtype=bool), prefix=_PREFIX
    )
    xyz = _coords3(mesh)
    normals = ring_unit_normals(xyz, start, nodes)
    face_label = label[cell_of] if len(cell_of) else np.empty(0, dtype=np.int64)
    ids = sorted({int(x) for x in face_label if x >= 0})

    rows = []
    for pid in ids:
        faces = np.nonzero(face_label == pid)[0]
        verts = np.unique(
            np.concatenate([nodes[start[f] : start[f + 1]] for f in faces])
        )
        pts = xyz[verts].copy()
        shift = np.zeros(3)
        for p in pts:
            shift = shift + p
        shift = shift / float(len(pts))
        ext = pts.max(axis=0) - pts.min(axis=0)
        diag = math.sqrt(float(ext[0] * ext[0] + ext[1] * ext[1] + ext[2] * ext[2]))
        centres = np.empty((len(faces), 3))
        for i, f in enumerate(faces):
            ring = nodes[start[f] : start[f + 1]]
            c = np.zeros(3)
            for k in ring:
                c = c + xyz[k]
            centres[i] = c / float(len(ring)) - shift
        patch = {
            "pts": pts - shift,
            "normals": normals[faces],
            "centres": centres,
            "shift": shift,
            "diag": diag,
        }
        rows.append(_fit_patch(patch, opts))

    row_of = {pid: r for r, pid in enumerate(ids)}
    cell_row = np.full(total, -1, dtype=np.int64)
    for f, c in enumerate(cell_of):
        if face_label[f] >= 0:
            cell_row[c] = row_of[int(face_label[f])]

    ptype = np.array([o["type"] for o in rows], dtype=np.int64)
    prms = np.array([o["rms"] for o in rows], dtype=np.float64)
    pmax = np.array([o["max"] for o in rows], dtype=np.float64)
    out = mesh.copy()
    out.cell_data = dict(out.cell_data)
    t_blocks, r_blocks, x_blocks = [], [], []
    for base, block in zip(bases, mesh.cells):
        rr = cell_row[base : base + len(block.data)]
        has = rr >= 0
        t = np.full(len(rr), -1, dtype=np.int64)
        r = np.full(len(rr), np.nan)
        x = np.full(len(rr), np.nan)
        t[has] = ptype[rr[has]]
        r[has] = prms[rr[has]]
        x[has] = pmax[rr[has]]
        t_blocks.append(t)
        r_blocks.append(r)
        x_blocks.append(x)
    out.cell_data["fit:type"] = t_blocks
    out.cell_data["fit:rms"] = r_blocks
    out.cell_data["fit:max"] = x_blocks
    out.field_data = dict(out.field_data)
    out.field_data["fit:patch_id"] = np.asarray(ids, dtype=np.int64)
    out.field_data["fit:patch_type"] = ptype
    out.field_data["fit:patch_params"] = np.array(
        [o["params"] for o in rows], dtype=np.float64
    ).reshape(-1, 10)
    out.field_data["fit:patch_rms"] = prms
    out.field_data["fit:patch_max"] = pmax
    out.field_data["fit:patch_fraction"] = np.array(
        [o["fraction"] for o in rows], dtype=np.float64
    )
    out.field_data["fit:patch_tolerance"] = np.array(
        [o["tol"] for o in rows], dtype=np.float64
    )
    counts = [int(np.count_nonzero(ptype == t)) for t in range(6)]
    report = {
        "num_patches": len(ids),
        "num_plane": counts[1],
        "num_sphere": counts[2],
        "num_cylinder": counts[3],
        "num_cone": counts[4],
        "num_torus": counts[5],
        "num_unfitted": counts[0],
    }
    return out, report


_REPORT_KEYS = (
    "num_patches",
    "num_plane",
    "num_sphere",
    "num_cylinder",
    "num_cone",
    "num_torus",
    "num_unfitted",
)


def types_mask(types):
    """The ``1 << type`` bit mask of a sequence of type names."""
    if isinstance(types, str):
        types = [types]
    mask = 0
    for t in types:
        if t not in TYPES[1:]:
            raise ValueError(
                f"{_PREFIX}unknown primitive type {t!r} "
                f"(expected one of {', '.join(TYPES[1:])})"
            )
        mask |= 1 << TYPES.index(t)
    if not mask:
        raise ValueError(f"{_PREFIX}types must name at least one primitive")
    return mask


def fit_primitives(
    mesh,
    patch_data: str = "segment:patch",
    tolerance: float = 1e-3,
    relative: bool = True,
    min_inlier_fraction: float = 0.9,
    types=("plane", "sphere", "cylinder", "cone", "torus"),
    max_iterations: int = 50,
    return_report: bool = False,
):
    """Fit a plane, sphere, cylinder, cone or torus to each patch of a surface.

    :param mesh: a surface mesh carrying integer patch labels as cell data.
    :param patch_data: the label array (``segment:patch`` by default, written by
        :func:`segment`); cells labelled below 0 belong to no patch. An empty
        string fits the whole surface as one patch.
    :param tolerance: the acceptance tolerance on the orthogonal distance.
    :param relative: ``tolerance`` is a fraction of each patch's bounding-box
        diagonal (default) rather than a length.
    :param min_inlier_fraction: the share of a patch's points that must lie
        within tolerance.
    :param types: the primitives to try, from ``plane``, ``sphere``,
        ``cylinder``, ``cone``, ``torus``; the first accepted in that order wins.
    :param max_iterations: the Levenberg-Marquardt iteration cap per fit.
    :param return_report: also return the per-type counts.
    :returns: the input with cell data ``fit:type`` (-1 no patch, 0 none,
        1 plane, 2 sphere, 3 cylinder, 4 cone, 5 torus), ``fit:rms`` and
        ``fit:max``, and per-patch field data ``fit:patch_id``,
        ``fit:patch_type``, ``fit:patch_params`` (``(m, 10)``),
        ``fit:patch_rms``, ``fit:patch_max``, ``fit:patch_fraction`` and
        ``fit:patch_tolerance``. With ``return_report``, ``(mesh, report)``.
    """
    opts = {
        "patch_data": str(patch_data),
        "tolerance": float(tolerance),
        "relative": bool(relative),
        "min_inlier_fraction": float(min_inlier_fraction),
        "types": types_mask(types),
        "max_iterations": int(max_iterations),
    }
    if not (opts["tolerance"] > 0.0 and math.isfinite(opts["tolerance"])):
        raise ValueError(f"{_PREFIX}the tolerance must be finite and > 0")
    if not 0.0 <= opts["min_inlier_fraction"] <= 1.0:
        raise ValueError(f"{_PREFIX}min_inlier_fraction must lie in [0, 1]")
    if opts["max_iterations"] < 1:
        raise ValueError(f"{_PREFIX}max_iterations must be >= 1")

    out = None
    try:
        from . import _core

        res = _core.fit_primitives(
            mesh,
            opts["patch_data"],
            opts["tolerance"],
            opts["relative"],
            opts["min_inlier_fraction"],
            opts["types"],
            opts["max_iterations"],
        )
        out = res["mesh"]
        report = {k: int(res[k]) for k in _REPORT_KEYS}
    except Exception as exc:
        if not core_op_declined(exc, "fit_primitives"):
            raise
        out = None
    if out is None:
        out, report = _fit_primitives_py(mesh, opts)
    return (out, report) if return_report else out


__all__ = ["fit_primitives"]

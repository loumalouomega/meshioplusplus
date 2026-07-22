"""Per-cell mesh quality metrics + inverted/degenerate detection.

Every metric is evaluated on the cell's corner nodes (quadratic variants
reduce to their linear parent). The public entry points use the C++ core when
available and fall back to the pure-numpy reference here, which mirrors
``src/cpp/src/operations/quality.cpp`` formula-for-formula (Verdict/VTK
conventions; the ideal element gives the ideal value).
"""

from __future__ import annotations

import numpy as np

from ._mesh import Mesh
from ._skin import _CELL_FACES

_METRIC_NAMES = [
    "quality:volume",
    "quality:scaled_jacobian",
    "quality:aspect_ratio",
    "quality:skewness",
    "quality:min_angle",
    "quality:max_angle",
    "quality:warpage",
    "quality:min_dihedral",
    "quality:max_dihedral",
    "quality:inverted",
    "quality:degenerate",
]
_M = {name: i for i, name in enumerate(_METRIC_NAMES)}
_NUM_METRICS = len(_METRIC_NAMES)
_HIST_BINS = 10
_EPS = 1e-14
_RAD2DEG = 180.0 / np.pi

_CORNER_COUNT = {"tri": 3, "quad": 4, "tetra": 4, "hex": 8, "wedge": 6, "pyramid": 5}


def _family(cell_type: str):
    if cell_type.startswith("triangle"):
        return "tri"
    if cell_type.startswith("quad"):
        return "quad"
    if cell_type in _CELL_FACES:
        if cell_type.startswith("tetra"):
            return "tetra"
        if cell_type.startswith("hexahedron"):
            return "hex"
        if cell_type.startswith("wedge"):
            return "wedge"
        if cell_type.startswith("pyramid"):
            return "pyramid"
    return None


def _norm(a):
    return float(np.sqrt(np.dot(a, a)))


def _det3(a, b, c):
    return float(np.dot(a, np.cross(b, c)))


def _angle_deg(u, v):
    nu, nv = _norm(u), _norm(v)
    if nu < _EPS or nv < _EPS:
        return np.nan
    c = np.clip(np.dot(u, v) / (nu * nv), -1.0, 1.0)
    return float(np.arccos(c) * _RAD2DEG)


def _det_unit3(a, b, c):
    na, nb, nc = _norm(a), _norm(b), _norm(c)
    if na < _EPS or nb < _EPS or nc < _EPS:
        return np.nan
    return _det3(a / na, b / nb, c / nc)


def _facefan_volume(corners, faces):
    cell_c = corners.mean(axis=0)
    vol = 0.0
    for _ftype, ncn, ids in faces:
        idc = ids[:ncn]
        fc = corners[list(idc)].mean(axis=0)
        for k in range(ncn):
            a = corners[idc[k]]
            b = corners[idc[(k + 1) % ncn]]
            vol += _det3(a - cell_c, b - cell_c, fc - cell_c) / 6.0
    return vol


def _equiangle_skew(amin, amax):
    return max((amax - 90.0) / 90.0, (90.0 - amin) / 90.0)


def _tri(c, is2d, m):
    p0, p1, p2 = c[0], c[1], c[2]
    a, b, cc = _norm(p1 - p0), _norm(p2 - p1), _norm(p0 - p2)
    if is2d:
        area = 0.5 * (
            (p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1])
        )
    else:
        area = 0.5 * _norm(np.cross(p1 - p0, p2 - p0))
    m[_M["quality:volume"]] = area
    abs_area = abs(area)
    degenerate = abs_area < _EPS or a < _EPS or b < _EPS or cc < _EPS
    m[_M["quality:inverted"]] = 1.0 if (is2d and area < 0.0) else 0.0
    m[_M["quality:degenerate"]] = 1.0 if degenerate else 0.0
    if degenerate:
        return
    a0 = _angle_deg(p1 - p0, p2 - p0)
    a1 = _angle_deg(p0 - p1, p2 - p1)
    a2 = _angle_deg(p0 - p2, p1 - p2)
    m[_M["quality:min_angle"]] = min(a0, a1, a2)
    m[_M["quality:max_angle"]] = max(a0, a1, a2)
    m[_M["quality:aspect_ratio"]] = (
        a * b * cc * (a + b + cc) / (16.0 * abs_area * abs_area)
    )


def _quad(c, is2d, m):
    p = [c[0], c[1], c[2], c[3]]
    lengths = [_norm(p[(i + 1) % 4] - p[i]) for i in range(4)]
    if is2d:
        area = 0.5 * sum(
            p[i][0] * p[(i + 1) % 4][1] - p[(i + 1) % 4][0] * p[i][1] for i in range(4)
        )
    else:
        n = np.zeros(3)
        for i in range(4):
            n = n + np.cross(p[i], p[(i + 1) % 4])
        area = 0.5 * _norm(n)
    m[_M["quality:volume"]] = area
    abs_area = abs(area)
    lmin, lmax = min(lengths), max(lengths)
    degenerate = abs_area < _EPS or lmin < _EPS
    m[_M["quality:inverted"]] = 1.0 if (is2d and area < 0.0) else 0.0
    m[_M["quality:degenerate"]] = 1.0 if degenerate else 0.0
    if degenerate:
        return
    m[_M["quality:aspect_ratio"]] = lmax / lmin
    angs = [_angle_deg(p[(i + 1) % 4] - p[i], p[(i + 3) % 4] - p[i]) for i in range(4)]
    m[_M["quality:skewness"]] = _equiangle_skew(min(angs), max(angs))

    def split_warp(a0, a1, a2, b0, b1, b2):
        na = np.cross(a1 - a0, a2 - a0)
        nb = np.cross(b1 - b0, b2 - b0)
        lna, lnb = _norm(na), _norm(nb)
        if lna < _EPS or lnb < _EPS:
            return np.nan
        d = np.clip(abs(np.dot(na, nb) / (lna * lnb)), 0.0, 1.0)
        return float(np.arccos(d) * _RAD2DEG)

    w0 = split_warp(p[0], p[1], p[2], p[2], p[3], p[0])
    w1 = split_warp(p[1], p[2], p[3], p[3], p[0], p[1])
    if np.isfinite(w0) and np.isfinite(w1):
        m[_M["quality:warpage"]] = max(w0, w1)


def _tetra(c, m):
    p0, p1, p2, p3 = c[0], c[1], c[2], c[3]
    vol = _det3(p1 - p0, p2 - p0, p3 - p0) / 6.0
    m[_M["quality:volume"]] = vol
    edges = [
        _norm(p0 - p1),
        _norm(p0 - p2),
        _norm(p0 - p3),
        _norm(p1 - p2),
        _norm(p1 - p3),
        _norm(p2 - p3),
    ]
    emin = min(edges)
    js = [
        _det_unit3(p1 - p0, p2 - p0, p3 - p0),
        _det_unit3(p2 - p1, p0 - p1, p3 - p1),
        _det_unit3(p0 - p2, p1 - p2, p3 - p2),
        _det_unit3(p0 - p3, p2 - p3, p1 - p3),
    ]
    sj = np.sqrt(2.0) * min(js) if all(np.isfinite(js)) else np.nan
    abs_vol = abs(vol)
    degenerate = abs_vol < _EPS or emin < _EPS or (np.isfinite(sj) and abs(sj) < _EPS)
    m[_M["quality:inverted"]] = 1.0 if vol < 0.0 else 0.0
    m[_M["quality:degenerate"]] = 1.0 if degenerate else 0.0
    if degenerate:
        return
    m[_M["quality:scaled_jacobian"]] = sj

    def tri_area(x, y, z):
        return 0.5 * _norm(np.cross(y - x, z - x))

    face_sum = (
        tri_area(p0, p1, p2)
        + tri_area(p0, p1, p3)
        + tri_area(p0, p2, p3)
        + tri_area(p1, p2, p3)
    )
    if face_sum > _EPS:
        r = 3.0 * abs_vol / face_sum
        a_mat = np.array([2.0 * (p1 - p0), 2.0 * (p2 - p0), 2.0 * (p3 - p0)])
        rhs = np.array(
            [
                np.dot(p1, p1) - np.dot(p0, p0),
                np.dot(p2, p2) - np.dot(p0, p0),
                np.dot(p3, p3) - np.dot(p0, p0),
            ]
        )
        if abs(np.linalg.det(a_mat)) > _EPS and r > _EPS:
            o = np.linalg.solve(a_mat, rhs)
            big_r = _norm(o - p0)
            m[_M["quality:aspect_ratio"]] = big_r / (3.0 * r)

    fverts = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]
    edge_faces = [(0, 1), (0, 2), (1, 2), (0, 3), (1, 3), (2, 3)]
    centroid = (p0 + p1 + p2 + p3) * 0.25
    normals = []
    ok = True
    for f in fverts:
        v0, v1, v2 = c[f[0]], c[f[1]], c[f[2]]
        n = np.cross(v1 - v0, v2 - v0)
        ln = _norm(n)
        if ln < _EPS:
            ok = False
            break
        n = n / ln
        fc = (v0 + v1 + v2) / 3.0
        if np.dot(n, fc - centroid) < 0.0:
            n = -n
        normals.append(n)
    if ok:
        dih = [
            180.0
            - np.arccos(np.clip(np.dot(normals[i], normals[j]), -1.0, 1.0)) * _RAD2DEG
            for i, j in edge_faces
        ]
        m[_M["quality:min_dihedral"]] = min(dih)
        m[_M["quality:max_dihedral"]] = max(dih)


def _corner_jacobian_family(c, m, triples, sj_factor, faces, do_skew):
    vol = _facefan_volume(c, faces)
    m[_M["quality:volume"]] = vol
    dets = [
        _det_unit3(c[i] - c[n], c[j] - c[n], c[k] - c[n]) for n, (i, j, k) in triples
    ]
    sj_ok = all(np.isfinite(dets))
    sj = min(dets) if sj_ok else np.nan
    emin = _min_edge(c, faces)
    abs_vol = abs(vol)
    degenerate = (
        abs_vol < _EPS or emin < _EPS or (sj_ok and abs(sj) < _EPS) or not sj_ok
    )
    m[_M["quality:inverted"]] = 1.0 if vol < 0.0 else 0.0
    m[_M["quality:degenerate"]] = 1.0 if degenerate else 0.0
    if degenerate:
        return
    m[_M["quality:scaled_jacobian"]] = sj_factor * sj
    if do_skew:
        amin, amax = 360.0, -360.0
        for _ftype, ncn, ids in faces:
            for k in range(ncn):
                pc = c[ids[k]]
                pn = c[ids[(k + 1) % ncn]]
                pp = c[ids[(k + ncn - 1) % ncn]]
                ang = _angle_deg(pn - pc, pp - pc)
                if np.isfinite(ang):
                    amin, amax = min(amin, ang), max(amax, ang)
        if amax > -360.0:
            m[_M["quality:skewness"]] = _equiangle_skew(amin, amax)


def _min_edge(c, faces):
    emin = np.inf
    for _ftype, ncn, ids in faces:
        for k in range(ncn):
            emin = min(emin, _norm(c[ids[(k + 1) % ncn]] - c[ids[k]]))
    return emin


_HEX_TRIPLES = [
    (0, (1, 3, 4)),
    (1, (2, 0, 5)),
    (2, (3, 1, 6)),
    (3, (0, 2, 7)),
    (4, (7, 5, 0)),
    (5, (4, 6, 1)),
    (6, (5, 7, 2)),
    (7, (6, 4, 3)),
]
_WEDGE_TRIPLES = [
    (0, (1, 2, 3)),
    (1, (2, 0, 4)),
    (2, (0, 1, 5)),
    (3, (5, 4, 0)),
    (4, (3, 5, 1)),
    (5, (4, 3, 2)),
]


def _hex(c, faces, m):
    _corner_jacobian_family(c, m, _HEX_TRIPLES, 1.0, faces, do_skew=True)


def _wedge(c, faces, m):
    _corner_jacobian_family(
        c, m, _WEDGE_TRIPLES, 2.0 / np.sqrt(3.0), faces, do_skew=False
    )


def _pyramid(c, faces, m):
    vol = _facefan_volume(c, faces)
    m[_M["quality:volume"]] = vol
    dets = []
    for i in range(4):
        nxt, prv = (i + 1) % 4, (i + 3) % 4
        dets.append(_det_unit3(c[nxt] - c[i], c[prv] - c[i], c[4] - c[i]))
    sj_ok = all(np.isfinite(dets))
    sj = min(dets) if sj_ok else np.nan
    emin = _min_edge(c, faces)
    abs_vol = abs(vol)
    degenerate = (
        abs_vol < _EPS or emin < _EPS or (sj_ok and abs(sj) < _EPS) or not sj_ok
    )
    m[_M["quality:inverted"]] = 1.0 if vol < 0.0 else 0.0
    m[_M["quality:degenerate"]] = 1.0 if degenerate else 0.0
    if degenerate:
        return
    m[_M["quality:scaled_jacobian"]] = np.sqrt(2.0) * sj


def _eval_cell(family, cell_type, c, is2d, m):
    if family == "tri":
        _tri(c, is2d, m)
    elif family == "quad":
        _quad(c, is2d, m)
    elif family == "tetra":
        _tetra(c, m)
    elif family == "hex":
        _hex(c, _CELL_FACES[cell_type], m)
    elif family == "wedge":
        _wedge(c, _CELL_FACES[cell_type], m)
    elif family == "pyramid":
        _pyramid(c, _CELL_FACES[cell_type], m)


def _compute_quality_py(mesh) -> dict:
    """Pure-numpy reference for :func:`compute_quality`; returns a report dict."""
    pts = np.asarray(mesh.points, dtype=np.float64)
    pdim = pts.shape[1] if pts.ndim == 2 else 0
    is2d = pdim == 2
    if pdim == 2:
        pts = np.column_stack([pts, np.zeros(len(pts))])

    block_vals = []
    num_cells = 0
    for block in mesh.cells:
        conn = np.asarray(block.data)
        nc = conn.shape[0] if conn.ndim == 2 else len(block.data)
        num_cells += nc
        vals = np.full((nc, _NUM_METRICS), np.nan)
        family = _family(block.type) if conn.ndim == 2 else None
        if family is not None:
            cc = _CORNER_COUNT[family]
            for i in range(nc):
                c = pts[conn[i, :cc]]
                _eval_cell(family, block.type, c, is2d, vals[i])
        block_vals.append(vals)

    metrics = {}
    num_inverted = 0
    num_degenerate = 0
    all_vals = np.concatenate(block_vals) if block_vals else np.empty((0, _NUM_METRICS))
    if all_vals.shape[0]:
        num_inverted = int(np.nansum(all_vals[:, _M["quality:inverted"]] == 1.0))
        num_degenerate = int(np.nansum(all_vals[:, _M["quality:degenerate"]] == 1.0))
    for mi, name in enumerate(_METRIC_NAMES):
        col = all_vals[:, mi]
        finite = col[np.isfinite(col)]
        summary = {
            "min": float("nan"),
            "max": float("nan"),
            "mean": float("nan"),
            "count": int(finite.size),
            "histogram": [0] * _HIST_BINS,
        }
        if finite.size:
            lo, hi = float(finite.min()), float(finite.max())
            summary["min"], summary["max"] = lo, hi
            summary["mean"] = float(finite.mean())
            span = hi - lo
            hist = [0] * _HIST_BINS
            for x in finite:
                b = 0 if span <= 0.0 else int((x - lo) / span * _HIST_BINS)
                b = min(max(b, 0), _HIST_BINS - 1)
                hist[b] += 1
            summary["histogram"] = hist
        metrics[name] = summary

    cell_arrays = {
        name: [vals[:, mi].reshape(-1, 1).copy() for vals in block_vals]
        for mi, name in enumerate(_METRIC_NAMES)
    }
    return {
        "num_cells": int(num_cells),
        "num_inverted": num_inverted,
        "num_degenerate": num_degenerate,
        "metrics": metrics,
        "cell_arrays": cell_arrays,
    }


def _attach_quality_py(mesh) -> Mesh:
    rep = _compute_quality_py(mesh)
    out = mesh.copy() if hasattr(mesh, "copy") else Mesh(mesh.points, mesh.cells)
    if not mesh.cells:
        return out
    for name, arrays in rep["cell_arrays"].items():
        out.cell_data[name] = arrays
    return out


def compute_quality(mesh) -> dict:
    """Compute per-cell quality metrics and an aggregate report.

    Returns a dict with ``num_cells``/``num_inverted``/``num_degenerate``
    counts, a ``metrics`` map (per-metric ``min``/``max``/``mean``/``count``/
    ``histogram``), and a ``cell_arrays`` map (metric name -> one array per cell
    block, NaN where the metric does not apply).
    """
    try:
        from . import _core

        return _core.compute_quality(mesh)
    except Exception:
        pass
    return _compute_quality_py(mesh)


def attach_quality(mesh) -> Mesh:
    """A copy of ``mesh`` with every quality metric attached as ``cell_data``."""
    try:
        from . import _core

        return _core.attach_quality(mesh)
    except Exception:
        pass
    return _attach_quality_py(mesh)

"""Parametric mesh generators (``meshioplusplus.rectangle`` and friends).

``grid`` was the library's one mesh *generator*: everything else transforms a mesh
you already have. These build a shape from parameters, with the contract ``grid``
already keeps:

* coordinates are computed from indices and never accumulated (an angle is
  ``2*pi*k/n``, not the previous point rotated);
* points and cells follow a documented order (``doc/primitives.md``);
* a ``max_cells`` budget refuses a request by name, and a zero count yields an
  empty mesh (points ``(0, 3)``, no cell block);
* surfaces are oriented outward, planar cells counter-clockwise seen from +z, and
  volumes positively;
* boundaries come back as named ``side`` regions (``xmin`` ... ``zmax`` on a box,
  ``rim`` on a disk, ``side``/``top``/``bottom`` on a cylinder, ``inner``/``outer``
  on an annulus, ``hole`` on a plate);
* ``order=2`` places the mid-edge nodes on the exact geometry, appended after the
  corners in first-appearance order over the cell-major edge list.

This is the numpy twin of ``operations/primitives.cpp``. It is byte-identical to
the core for :func:`rectangle`, :func:`box`, :func:`icosa_sphere`,
:func:`reference_triangle`, :func:`lshape`, :func:`honeycomb` and
:func:`reference_cells` (index arithmetic and ``sqrt`` only); :func:`disk`,
:func:`uv_sphere`, :func:`cylinder`, :func:`ngon`, :func:`annulus` and
:func:`plate_with_hole` use ``sin``/``cos``, which differ between libms and
numpy's SIMD kernels, so they agree with the core on topology and to a few ulps
on coordinates. The twin is pinned by ``tests/python/test_primitives.py``.

Public API:

* :func:`rectangle`, :func:`box`, :func:`disk`, :func:`uv_sphere`,
  :func:`icosa_sphere`, :func:`cylinder`, the planar domains
  :func:`reference_triangle`, :func:`ngon`, :func:`annulus`, :func:`lshape`,
  :func:`plate_with_hole` and :func:`honeycomb`, and :func:`reference_cells`.
"""

from __future__ import annotations

import numpy as np

from ._fallback import core_op_declined
from ._mesh import Mesh
from ._regions import Region

__all__ = [
    "rectangle",
    "box",
    "disk",
    "uv_sphere",
    "icosa_sphere",
    "cylinder",
    "reference_triangle",
    "ngon",
    "annulus",
    "lshape",
    "plate_with_hole",
    "honeycomb",
    "reference_cells",
    "line",
    "circle",
    "helix",
    "sphere_points",
    "box_points",
    "tetra_sphere",
    "octa_sphere",
    "cube_sphere",
    "geo_sphere",
    "torus",
    "moebius",
    "tube",
    "cone",
    "capsule",
    "box_surface",
    "platonic",
    "ball",
    "spherical_shell",
    "pipe",
    "fichera",
    "solid_torus",
]

#: Default cell budget, as for :func:`meshioplusplus.grid`.
DEFAULT_MAX_CELLS = 20000000

_PI = float(np.pi)
_DEGENERATE = 0x8000

#: The highest ``order`` (3 and up build the ``VTK_LAGRANGE_*`` family).
MAX_ORDER = 10

_ELEVATED = {
    "line": "line3",
    "triangle": "triangle6",
    "quad": "quad8",
    "tetra": "tetra10",
    "hexahedron": "hexahedron20",
    "wedge": "wedge15",
}

# The six Kuhn tetrahedra around the main diagonal 0-6 of a hexahedron.
_KUHN = np.array(
    [
        [0, 1, 2, 6],
        [0, 2, 3, 6],
        [0, 3, 7, 6],
        [0, 7, 4, 6],
        [0, 4, 5, 6],
        [0, 5, 1, 6],
    ],
    dtype=np.int64,
)


# --------------------------------------------------------------------------- #
# validation                                                                  #
# --------------------------------------------------------------------------- #


def _fail(shape, message):
    raise ValueError(f"meshio++: {shape}: {message}")


def _vec(shape, what, value, n):
    arr = np.asarray(value, dtype=np.float64).reshape(-1)
    if arr.size != n:
        _fail(shape, f"{what} must have {n} components")
    return arr


def _ivec(shape, what, value, n):
    arr = np.asarray(value, dtype=np.int64).reshape(-1)
    if arr.size != n:
        _fail(shape, f"{what} must have {n} cell counts")
    return arr


def _check_count(shape, what, value):
    if value < 0:
        _fail(shape, f"{what} must not be negative, got {int(value)}")


def _check_positive(shape, what, value):
    if not value > 0.0:
        _fail(shape, f"{what} must be positive, got {_fmt(value)}")


def _fmt(value):
    # std::to_string(double): six decimals.
    return f"{float(value):.6f}"


def _check_order(shape, order):
    try:
        ok = int(order) == order and 1 <= order <= MAX_ORDER
    except (TypeError, ValueError):
        ok = False
    if not ok:
        _fail(shape, f"order must be between 1 and {MAX_ORDER}, got {order}")


def _choice(shape, what, got, choices):
    if got not in choices:
        _fail(shape, f"unknown {what} '{got}' (expected one of: {', '.join(choices)})")


def _budget(shape, cells, max_cells):
    if max_cells > 0 and cells > max_cells:
        _fail(
            shape,
            f"the requested mesh has {cells} cells, above the limit of {max_cells} "
            "(raise max_cells or coarsen the counts)",
        )


def _empty():
    return Mesh(np.zeros((0, 3)), [])


def _core_call(op, *args):
    from . import _core

    return getattr(_core, op)(*args)


# --------------------------------------------------------------------------- #
# the soup: parameter space, masks, blocks                                    #
# --------------------------------------------------------------------------- #


class _Soup:
    """Nodes in a parameter space where the mid-edge average is exact."""

    def __init__(self, kind, radius=1.0, center=(0.0, 0.0, 0.0)):
        # identity | squircle | polar | sphere | unit | plate | project | torus |
        # moebius | capsule | ball | shell | graded | solid_torus
        self.kind = kind
        self.radius = float(radius)
        self.center = np.asarray(center, dtype=np.float64)
        # kind == "plate": the half side, the angular positions around the hole
        # and the cells per 90 degree block (the hole radius is ``radius``).
        self.half = 1.0
        self.period = 1.0
        self.nt = 1.0
        # torus / solid_torus: the minor radius; shell / graded: the outer radius,
        # the grading ratio and the layer count; moebius: ``half`` is the half
        # width and ``twists`` the half-twists; capsule: ``half`` is half the
        # cylinder length.
        self.radius2 = 0.5
        self.ratio = 1.0
        self.layers = 1.0
        self.twists = 1
        self.param = np.zeros((0, 3))
        self.mask = np.zeros(0, dtype=np.uint16)
        self.blocks = []  # [type, conn (m, w)]

    def add(self, param, mask):
        first = self.mask.size
        self.param = np.concatenate([self.param, np.asarray(param, dtype=np.float64)])
        self.mask = np.concatenate([self.mask, np.asarray(mask, dtype=np.uint16)])
        return first

    def block(self, type_, conn):
        self.blocks.append([type_, np.ascontiguousarray(conn, dtype=np.int64)])


def _grade(soup, t):
    if soup.ratio == 1.0:
        return t
    q = soup.ratio
    with np.errstate(over="ignore", invalid="ignore"):
        g = (np.power(q, soup.layers * t) - 1.0) / (np.power(q, soup.layers) - 1.0)
    return np.where(t >= 1.0, 1.0, g)


def _maxabs(p):
    return np.maximum(np.maximum(np.abs(p[:, 0]), np.abs(p[:, 1])), np.abs(p[:, 2]))


def _positions(soup):
    p = soup.param
    r = soup.radius
    cx, cy, cz = soup.center
    out = np.empty_like(p)
    if soup.kind == "identity":
        out[:] = p
    elif soup.kind == "squircle":
        out[:, 0] = cx + r * p[:, 0] * np.sqrt(1.0 - p[:, 1] * p[:, 1] / 2.0)
        out[:, 1] = cy + r * p[:, 1] * np.sqrt(1.0 - p[:, 0] * p[:, 0] / 2.0)
        out[:, 2] = cz + p[:, 2]
    elif soup.kind == "polar":
        out[:, 0] = cx + r * p[:, 0] * np.cos(p[:, 1])
        out[:, 1] = cy + r * p[:, 0] * np.sin(p[:, 1])
        out[:, 2] = cz + p[:, 2]
    elif soup.kind == "sphere":
        st = np.sin(p[:, 0])
        out[:, 0] = cx + r * st * np.cos(p[:, 1])
        out[:, 1] = cy + r * st * np.sin(p[:, 1])
        out[:, 2] = cz + r * np.cos(p[:, 0])
        pole = (soup.mask & _DEGENERATE) != 0
        out[pole, 0] = cx
        out[pole, 1] = cy
        out[pole, 2] = np.where(p[pole, 0] == 0.0, cz + r, cz - r)
    elif soup.kind == "plate":
        q = np.where(p[:, 1] >= soup.period, p[:, 1] - soup.period, p[:, 1])
        nt = soup.nt
        b = np.clip(np.floor(q / nt), 0.0, 3.0)
        m = q - b * nt
        ang = -_PI / 4.0 + ((_PI / 2.0) * q) / nt
        hx = r * np.cos(ang)
        hy = r * np.sin(ang)
        ex = np.full_like(q, soup.half)
        ey = soup.half * ((2.0 * m) / nt - 1.0)
        bi = b.astype(np.int64)
        ox = np.where(bi == 1, -ey, np.where(bi == 2, -ex, np.where(bi == 3, ey, ex)))
        oy = np.where(bi == 1, ex, np.where(bi == 2, -ey, np.where(bi == 3, -ex, ey)))
        s = p[:, 0]
        out[:, 0] = cx + ((1.0 - s) * hx + s * ox)
        out[:, 1] = cy + ((1.0 - s) * hy + s * oy)
        out[:, 2] = cz + p[:, 2]
    elif soup.kind == "project":
        length = np.sqrt(p[:, 0] * p[:, 0] + p[:, 1] * p[:, 1] + p[:, 2] * p[:, 2])
        out[:, 0] = cx + r * (p[:, 0] / length)
        out[:, 1] = cy + r * (p[:, 1] / length)
        out[:, 2] = cz + r * (p[:, 2] / length)
    elif soup.kind == "torus":
        rho = r + soup.radius2 * np.cos(p[:, 1])
        out[:, 0] = cx + rho * np.cos(p[:, 0])
        out[:, 1] = cy + rho * np.sin(p[:, 0])
        out[:, 2] = cz + soup.radius2 * np.sin(p[:, 1])
    elif soup.kind == "moebius":
        tw = (float(soup.twists) * p[:, 0]) / 2.0
        off = soup.half * p[:, 1]
        rho = r + off * np.cos(tw)
        out[:, 0] = cx + rho * np.cos(p[:, 0])
        out[:, 1] = cy + rho * np.sin(p[:, 0])
        out[:, 2] = cz + off * np.sin(tw)
    elif soup.kind == "capsule":
        h = soup.half
        side = (2.0 * h) / r
        t = p[:, 0]
        cap1 = t <= _PI / 2.0
        body = ~cap1 & (t <= _PI / 2.0 + side)
        t2 = t - side
        rho = np.where(cap1, r * np.sin(t), np.where(body, r, r * np.sin(t2)))
        z = np.where(
            cap1,
            h + r * np.cos(t),
            np.where(body, h - r * (t - _PI / 2.0), -h + r * np.cos(t2)),
        )
        out[:, 0] = cx + rho * np.cos(p[:, 1])
        out[:, 1] = cy + rho * np.sin(p[:, 1])
        out[:, 2] = cz + z
        pole = (soup.mask & _DEGENERATE) != 0
        out[pole, 0] = cx
        out[pole, 1] = cy
        out[pole, 2] = np.where(p[pole, 0] == 0.0, cz + h + r, cz - h - r)
    elif soup.kind == "ball":
        s = _maxabs(p)
        inner = s <= 0.5
        safe = np.where(inner, 1.0, s)
        d = p / safe[:, None]
        length = np.sqrt(d[:, 0] * d[:, 0] + d[:, 1] * d[:, 1] + d[:, 2] * d[:, 2])
        length = np.where(inner, 1.0, length)
        lam = 2.0 * s - 1.0
        for k, c in enumerate((cx, cy, cz)):
            blend = (1.0 - lam) * 0.5 * d[:, k] + lam * (d[:, k] / length)
            out[:, k] = np.where(inner, c + r * p[:, k], c + r * blend)
    elif soup.kind == "shell":
        s = _maxabs(p)
        d = p / s[:, None]
        length = np.sqrt(d[:, 0] * d[:, 0] + d[:, 1] * d[:, 1] + d[:, 2] * d[:, 2])
        rad = r + (soup.radius2 - r) * _grade(soup, s - 1.0)
        for k, c in enumerate((cx, cy, cz)):
            out[:, k] = c + rad * (d[:, k] / length)
    elif soup.kind == "graded":
        rad = r + (soup.radius2 - r) * _grade(soup, p[:, 0])
        out[:, 0] = cx + rad * np.cos(p[:, 1])
        out[:, 1] = cy + rad * np.sin(p[:, 1])
        out[:, 2] = cz + p[:, 2]
    elif soup.kind == "solid_torus":
        a = soup.radius2 * p[:, 0] * np.sqrt(1.0 - p[:, 1] * p[:, 1] / 2.0)
        b = soup.radius2 * p[:, 1] * np.sqrt(1.0 - p[:, 0] * p[:, 0] / 2.0)
        out[:, 0] = cx + (r + a) * np.cos(p[:, 2])
        out[:, 1] = cy + (r + a) * np.sin(p[:, 2])
        out[:, 2] = cz + b
    else:  # unit
        out[:, 0] = cx + r * p[:, 0]
        out[:, 1] = cy + r * p[:, 1]
        out[:, 2] = cz + r * p[:, 2]
    return out


_LINEAR_MAPS = ("identity", "squircle", "project", "ball", "shell")


def _combine(soup, nodes, w):
    """Weighted combinations of nodes in parameter space, one per row of
    ``nodes`` (``(m, k)``) with the weights ``w`` (``k``): the twin of
    ``prim_combine``."""
    nodes = np.asarray(nodes, dtype=np.int64)
    m, k = nodes.shape
    w = [float(v) for v in w]
    par = soup.param[nodes]  # (m, k, 3)
    deg = (soup.mask[nodes] & _DEGENERATE) != 0
    two_pi = 2.0 * _PI
    wsum = 0.0
    for v in w:
        wsum += v

    def linear(c):
        acc = np.zeros(m)
        for i in range(k):
            acc = acc + w[i] * par[:, i, c]
        return acc / wsum

    def angular(c, period):
        nondeg = ~deg
        anyok = nondeg.any(axis=1)
        ref = np.argmax(nondeg, axis=1)
        base = par[np.arange(m), ref, c]
        acc = np.zeros(m)
        ws = np.zeros(m)
        for i in range(k):
            d = par[:, i, c] - base
            d = d - period * np.floor(d / period + 0.5)
            ok = nondeg[:, i]
            acc = np.where(ok, acc + w[i] * (base + d), acc)
            ws = np.where(ok, ws + w[i], ws)
        with np.errstate(invalid="ignore", divide="ignore"):
            res = acc / ws
        return np.where(anyok, res, par[:, 0, c])

    out = np.empty((m, 3))
    kind = soup.kind
    if kind == "unit":
        x, y, z = linear(0), linear(1), linear(2)
        length = np.sqrt(x * x + y * y + z * z)
        out[:, 0] = x / length
        out[:, 1] = y / length
        out[:, 2] = z / length
    elif kind in ("polar", "graded", "sphere", "capsule"):
        out[:, 0] = linear(0)
        out[:, 1] = angular(1, two_pi)
        out[:, 2] = linear(2)
    elif kind == "torus":
        out[:, 0] = angular(0, two_pi)
        out[:, 1] = angular(1, two_pi)
        out[:, 2] = linear(2)
    elif kind == "solid_torus":
        out[:, 0] = linear(0)
        out[:, 1] = linear(1)
        out[:, 2] = angular(2, two_pi)
    elif kind == "plate":
        out[:, 0] = linear(0)
        q = angular(1, soup.period)
        out[:, 1] = np.where(q < 0.0, q + soup.period, q)
        out[:, 2] = linear(2)
    elif kind == "moebius":
        base = par[:, 0, 0]
        au = np.zeros(m)
        aw = np.zeros(m)
        for i in range(k):
            d = par[:, i, 0] - base
            turns = np.floor(d / two_pi + 0.5)
            d = d - two_pi * turns
            flip = (np.abs(turns).astype(np.int64) * int(soup.twists)) % 2 != 0
            au = au + w[i] * (base + d)
            aw = aw + w[i] * np.where(flip, -par[:, i, 1], par[:, i, 1])
        out[:, 0] = au / wsum
        out[:, 1] = aw / wsum
        out[:, 2] = linear(2)
    else:
        out[:, 0] = linear(0)
        out[:, 1] = linear(1)
        out[:, 2] = linear(2)
    return out


def _mid_angle(a, da, b, db):
    two_pi = 2.0 * _PI
    d = b - a
    d = d - two_pi * np.floor(d / two_pi + 0.5)
    mid = a + 0.5 * d
    return np.where(da, b, np.where(db, a, mid))


def _mid(soup, lo, hi):
    if soup.kind in _COMBINED_MIDS:
        return _combine(soup, np.stack([lo, hi], axis=1), [1.0, 1.0])
    a = soup.param[lo]
    b = soup.param[hi]
    out = np.empty_like(a)
    if soup.kind in ("identity", "squircle"):
        out[:] = (a + b) * 0.5
    elif soup.kind in ("polar", "sphere"):
        da = (soup.mask[lo] & _DEGENERATE) != 0
        db = (soup.mask[hi] & _DEGENERATE) != 0
        out[:, 0] = (a[:, 0] + b[:, 0]) * 0.5
        out[:, 1] = _mid_angle(a[:, 1], da, b[:, 1], db)
        out[:, 2] = (a[:, 2] + b[:, 2]) * 0.5 if soup.kind == "polar" else 0.0
    elif soup.kind == "plate":
        per = soup.period
        d = b[:, 1] - a[:, 1]
        d = np.where(d > 0.5 * per, d - per, np.where(d < -0.5 * per, d + per, d))
        q = a[:, 1] + 0.5 * d
        out[:, 0] = (a[:, 0] + b[:, 0]) * 0.5
        out[:, 1] = np.where(q < 0.0, q + per, q)
        out[:, 2] = (a[:, 2] + b[:, 2]) * 0.5
    else:
        x = a[:, 0] + b[:, 0]
        y = a[:, 1] + b[:, 1]
        z = a[:, 2] + b[:, 2]
        length = np.sqrt(x * x + y * y + z * z)
        out[:, 0] = x / length
        out[:, 1] = y / length
        out[:, 2] = z / length
    return out


# The maps whose mid-edge node is the weighted combination of its two ends.
_COMBINED_MIDS = (
    "project",
    "torus",
    "moebius",
    "capsule",
    "ball",
    "shell",
    "graded",
    "solid_torus",
)


def _dedup_edges(lo, hi, num_nodes):
    """Edges in first-appearance order: (unique_lo, unique_hi, id per input edge)."""
    key = lo * (num_nodes + 1) + hi
    _, first, inverse = np.unique(key, return_index=True, return_inverse=True)
    order = np.argsort(first, kind="stable")
    rank = np.empty_like(order)
    rank[order] = np.arange(order.size)
    return lo[first[order]], hi[first[order]], rank[inverse.reshape(-1)]


def _elevate(soup):
    from ._convert_cells import _ELEVATE

    known = {}  # edge -> node, across blocks
    for block in soup.blocks:
        type_, conn = block
        if type_ not in _ELEVATED:
            continue
        _, corners, edges = _ELEVATE[type_]
        ncells = conn.shape[0]
        a = np.stack([conn[:, e[0]] for e in edges], axis=1)
        b = np.stack([conn[:, e[1]] for e in edges], axis=1)
        lo = np.minimum(a, b).reshape(-1)
        hi = np.maximum(a, b).reshape(-1)
        ulo, uhi, inverse = _dedup_edges(lo, hi, soup.mask.size)
        # Edges a previous block already numbered keep that node (sphere bodies
        # and caps share their first ring).
        ids = np.empty(ulo.size, dtype=np.int64)
        fresh = []
        for n in range(ulo.size):
            k = (int(ulo[n]), int(uhi[n]))
            if k in known:
                ids[n] = known[k]
            else:
                ids[n] = soup.mask.size + len(fresh)
                known[k] = int(ids[n])
                fresh.append(n)
        if fresh:
            sel = np.asarray(fresh)
            soup.add(_mid(soup, ulo[sel], uhi[sel]), np.zeros(sel.size))
        mids = ids[inverse].reshape(ncells, len(edges))
        block[1] = np.concatenate([conn, mids], axis=1)
        block[0] = _ELEVATED[type_]


def _append_keyed(soup, keys, known, params_of):
    """Nodes named by ``keys`` (one row each, cell-major): a key seen before (in
    ``known``, across blocks) reuses its node, a new one is appended in
    first-appearance order with the parameters ``params_of(rows)`` of the row
    that introduced it. Returns the node of every row."""
    if keys.shape[0] == 0:
        return np.zeros(0, dtype=np.int64)
    uniq, first, inverse = np.unique(
        keys, axis=0, return_index=True, return_inverse=True
    )
    order = np.argsort(first, kind="stable")
    ids = np.empty(uniq.shape[0], dtype=np.int64)
    fresh = []
    base = soup.mask.size
    for n in order:
        key = tuple(int(v) for v in uniq[n])
        if key in known:
            ids[n] = known[key]
        else:
            ids[n] = base + len(fresh)
            known[key] = int(ids[n])
            fresh.append(n)
    if fresh:
        rows = first[np.asarray(fresh)]
        soup.add(params_of(rows), np.zeros(rows.size))
    return ids[inverse.reshape(-1)]


_COMPLETE = {"quad8": "quad9", "hexahedron20": "hexahedron27", "wedge15": "wedge18"}
_FULL_WIDTH = {"quad9": 9, "hexahedron27": 27, "wedge18": 18}


def _centres(full):
    """``(local node, corners)`` of the centres a full-Lagrange type adds."""
    from ._skin import _CELL_FACES

    if full == "quad9":
        return [(8, (0, 1, 2, 3))]
    out = [
        (nodes[8], tuple(nodes[:4]))
        for _, _, nodes in _CELL_FACES[full]
        if len(nodes) == 9
    ]
    if full == "hexahedron27":
        out.append((26, tuple(range(8))))
    return sorted(out)


def _complete(soup):
    known = {}
    for block in soup.blocks:
        type_, conn = block
        if type_ not in _COMPLETE:
            continue
        full = _COMPLETE[type_]
        centres = _centres(full)
        ncells = conn.shape[0]
        out = np.empty((ncells, _FULL_WIDTH[full]), dtype=np.int64)
        out[:, : conn.shape[1]] = conn
        keys = np.full((ncells, len(centres), 8), -1, dtype=np.int64)
        for q, (_, corners) in enumerate(centres):
            keys[:, q, : len(corners)] = np.sort(conn[:, list(corners)], axis=1)
        keys = keys.reshape(-1, 8)

        def params_of(rows, conn=conn, centres=centres):
            cell, q = np.divmod(rows, len(centres))
            res = np.empty((rows.size, 3))
            for qq, (_, corners) in enumerate(centres):
                sel = q == qq
                if sel.any():
                    nodes = conn[cell[sel]][:, list(corners)]
                    res[sel] = _combine(soup, nodes, [1.0] * len(corners))
            return res

        ids = _append_keyed(soup, keys, known, params_of).reshape(ncells, -1)
        for q, (local, _) in enumerate(centres):
            out[:, local] = ids[:, q]
        block[0] = full
        block[1] = out


_LAGRANGE = {
    "line": ("line", "VTK_LAGRANGE_CURVE"),
    "triangle": ("triangle", "VTK_LAGRANGE_TRIANGLE"),
    "quad": ("quad", "VTK_LAGRANGE_QUADRILATERAL"),
    "tetra": ("tetra", "VTK_LAGRANGE_TETRAHEDRON"),
    "hexahedron": ("hexahedron", "VTK_LAGRANGE_HEXAHEDRON"),
    "wedge": ("wedge", "VTK_LAGRANGE_WEDGE"),
}


def _lagrange(shape, soup, order):
    from ._lagrange import lattice_weights, vtk_lattice

    known = {}
    for block in soup.blocks:
        type_, conn = block
        if type_ not in _LAGRANGE:
            _fail(shape, f"order {order} is not available for {type_} cells")
        lshape, out_type = _LAGRANGE[type_]
        lattice = vtk_lattice(lshape, order)
        weights = [lattice_weights(lshape, order, ijk) for ijk in lattice]
        ncells = conn.shape[0]
        out = np.empty((ncells, len(lattice)), dtype=np.int64)
        inner = []
        for ll, wl in enumerate(weights):
            if len(wl) == 1:
                out[:, ll] = conn[:, wl[0][0]]
            else:
                inner.append(ll)
        if inner:
            keys = np.full((ncells, len(inner), 16), -1, dtype=np.int64)
            for q, ll in enumerate(inner):
                corners = [c for c, _ in weights[ll]]
                wv = np.array([v for _, v in weights[ll]], dtype=np.int64)
                g = conn[:, corners]
                srt = np.argsort(g, axis=1, kind="stable")
                gs = np.take_along_axis(g, srt, axis=1)
                ws = wv[srt]
                keys[:, q, 0 : 2 * len(corners) : 2] = gs
                keys[:, q, 1 : 2 * len(corners) : 2] = ws
            keys = keys.reshape(-1, 16)

            def params_of(rows, conn=conn, inner=inner, weights=weights):
                cell, q = np.divmod(rows, len(inner))
                res = np.empty((rows.size, 3))
                for qq, ll in enumerate(inner):
                    sel = q == qq
                    if sel.any():
                        corners = [c for c, _ in weights[ll]]
                        wv = [float(v) for _, v in weights[ll]]
                        res[sel] = _combine(soup, conn[cell[sel]][:, corners], wv)
                return res

            ids = _append_keyed(soup, keys, known, params_of).reshape(ncells, -1)
            for q, ll in enumerate(inner):
                out[:, ll] = ids[:, q]
        block[0] = out_type
        block[1] = out


# --------------------------------------------------------------------------- #
# regions                                                                     #
# --------------------------------------------------------------------------- #


def _regions(soup, blocks, specs, facet_dim):
    if not specs:
        return []
    from ._skin import _CELL_FACES
    from ._surface import _CELL_EDGES

    tables = _CELL_FACES if facet_dim == 2 else _CELL_EDGES
    cells, facets, corners = [], [], []
    offset = 0
    for type_, conn in blocks:
        n = conn.shape[0]
        for f, (_, ncorners, nodes) in enumerate(tables[type_]):
            c = np.full((n, 4), -1, dtype=np.int64)
            c[:, :ncorners] = conn[:, list(nodes[:ncorners])]
            cells.append(offset + np.arange(n, dtype=np.int64))
            facets.append(np.full(n, f, dtype=np.int64))
            corners.append(c)
        offset += n
    cells = np.concatenate(cells)
    facets = np.concatenate(facets)
    corners = np.concatenate(corners)
    key = np.sort(corners, axis=1)
    _, inverse, counts = np.unique(key, axis=0, return_inverse=True, return_counts=True)
    boundary = counts[inverse.reshape(-1)] == 1

    pad = np.concatenate([soup.mask, np.array([0xFFFF], dtype=np.uint16)])
    common = np.bitwise_and.reduce(pad[corners], axis=1)  # -1 -> the all-ones pad
    bits = np.array([s[1] for s in specs], dtype=np.int64)
    out = []
    chosen = np.full(cells.size, -1, dtype=np.int64)
    # The lowest bit decides between regions sharing nodes.
    for idx in np.argsort(-bits, kind="stable"):
        hit = boundary & ((common & bits[idx]) != 0)
        chosen[hit] = idx
    for idx, (name, _) in enumerate(specs):
        sel = chosen == idx
        entries = np.stack([cells[sel], facets[sel]], axis=1).reshape(-1, 2)
        out.append(Region(name, "side", entries, dim=facet_dim))
    return out


def _finish(soup, order, specs, facet_dim, complete=False, extra=(), shape=""):
    # Order 3 and up finds the side regions on the linear cells (the VTK
    # Lagrange types have no facet table), as the core does.
    regions = None
    if order >= 3:
        linear = [(t, c) for t, c in soup.blocks if c.size]
        regions = _regions(soup, linear, specs, facet_dim)
    if order == 2:
        _elevate(soup)
        if complete:
            _complete(soup)
    elif order >= 3:
        _lagrange(shape, soup, order)
    points = _positions(soup)
    blocks = [(t, c) for t, c in soup.blocks if c.size]
    if regions is None:
        regions = _regions(soup, blocks, specs, facet_dim)
    return Mesh(points, blocks, regions=list(regions) + list(extra))


def _cell_regions(labels, names, dim):
    labels = np.asarray(labels, dtype=np.int64)
    return [
        Region(name, "cell", np.nonzero(labels == k)[0].astype(np.int64), dim=dim)
        for k, name in enumerate(names)
    ]


def _point_region(name, point):
    return Region(name, "point", np.array([point], dtype=np.int64), dim=0)


# --------------------------------------------------------------------------- #
# rectangle                                                                   #
# --------------------------------------------------------------------------- #


def _rectangle_py(dims, size, origin, cell_type, split, order, complete=False):
    nx, ny = int(dims[0]), int(dims[1])
    px, py = nx + 1, ny + 1
    tri = cell_type == "triangle"
    center = tri and split == "center"
    i = np.tile(np.arange(px, dtype=np.int64), py)
    j = np.repeat(np.arange(py, dtype=np.int64), px)
    mask = ((i == 0) * 1 + (i == nx) * 2 + (j == 0) * 4 + (j == ny) * 8).astype(
        np.uint16
    )
    param = np.stack(
        [
            origin[0] + (size[0] * i.astype(np.float64)) / float(nx),
            origin[1] + (size[1] * j.astype(np.float64)) / float(ny),
            np.full(i.size, origin[2]),
        ],
        axis=1,
    )
    soup = _Soup("identity")
    soup.add(param, mask)
    cj, ci = np.divmod(np.arange(nx * ny, dtype=np.int64), nx)
    if center:
        mp = np.stack(
            [
                origin[0] + (size[0] * (2 * ci + 1).astype(np.float64)) / float(2 * nx),
                origin[1] + (size[1] * (2 * cj + 1).astype(np.float64)) / float(2 * ny),
                np.full(ci.size, origin[2]),
            ],
            axis=1,
        )
        soup.add(mp, np.zeros(ci.size))
    a = cj * px + ci
    b, c, d = a + 1, a + px + 1, a + px
    if not tri:
        soup.block("quad", np.stack([a, b, c, d], axis=1))
    elif center:
        m = px * py + cj * nx + ci
        rows = np.stack(
            [
                np.stack([a, b, m], axis=1),
                np.stack([b, c, m], axis=1),
                np.stack([c, d, m], axis=1),
                np.stack([d, a, m], axis=1),
            ],
            axis=1,
        )
        soup.block("triangle", rows.reshape(-1, 3))
    else:
        if split == "up":
            up = np.ones(a.size, dtype=bool)
        elif split == "down":
            up = np.zeros(a.size, dtype=bool)
        else:
            up = (ci + cj) % 2 == 0
        first = np.where(up[:, None], np.stack([a, b, c], 1), np.stack([a, b, d], 1))
        second = np.where(up[:, None], np.stack([a, c, d], 1), np.stack([b, c, d], 1))
        soup.block("triangle", np.stack([first, second], axis=1).reshape(-1, 3))
    return _finish(
        soup,
        order,
        [("xmin", 1), ("xmax", 2), ("ymin", 4), ("ymax", 8)],
        1,
        complete=complete,
        shape="rectangle",
    )


def rectangle(
    dims=(1, 1),
    size=(1.0, 1.0),
    origin=(0.0, 0.0, 0.0),
    cell_type="quad",
    split="zigzag",
    order=1,
    max_cells=DEFAULT_MAX_CELLS,
    complete=False,
):
    """A regular 2-D grid of quadrilaterals or triangles at ``z = origin[2]``.

    Parameters
    ----------
    dims :
        Cell counts ``(nx, ny)``. A zero gives an empty mesh.
    size :
        Side lengths ``(lx, ly)``, positive on every axis with cells.
    origin :
        The lo corner; two components put it at ``z = 0``.
    cell_type :
        ``"quad"`` or ``"triangle"``.
    split :
        For triangles: ``"zigzag"`` (alternating diagonals), ``"up"`` (diagonal
        from the lo corner), ``"down"`` (the other one) or ``"center"`` (four
        triangles around an added cell-centre node).
    order :
        1, 3 to 10 for the ``VTK_LAGRANGE_*`` family, or 2 for ``quad8`` / ``triangle6`` cells.
    max_cells :
        Refuse by name above this many cells; non-positive disables the check.
    complete :
        With ``order=2``, the full-Lagrange types (``quad9``, ``hexahedron27``,
        ``wedge18``) instead of the serendipity ones.

    Returns
    -------
    Mesh
        Points x fastest, then y; side regions ``xmin``, ``xmax``, ``ymin``,
        ``ymax``.
    """
    shape = "rectangle"
    dims = _ivec(shape, "dims", dims, 2)
    size = _vec(shape, "size", size, 2)
    origin = np.asarray(origin, dtype=np.float64).reshape(-1)
    if origin.size == 2:
        origin = np.append(origin, 0.0)
    if origin.size != 3:
        _fail(shape, "origin must have 2 or 3 components")
    for k in range(2):
        _check_count(shape, "cell counts", dims[k])
        if dims[k] > 0:
            _check_positive(shape, "size", size[k])
    _check_order(shape, order)
    _choice(shape, "cell_type", cell_type, ("quad", "triangle"))
    _choice(shape, "split", split, ("zigzag", "up", "down", "center"))
    if dims[0] == 0 or dims[1] == 0:
        return _empty()
    per = 1 if cell_type == "quad" else (4 if split == "center" else 2)
    _budget(shape, int(dims[0]) * int(dims[1]) * per, max_cells)
    try:
        return _core_call(
            "rectangle",
            [int(d) for d in dims],
            [float(v) for v in size],
            [float(v) for v in origin],
            cell_type,
            split,
            int(order),
            int(max_cells),
            bool(complete),
        )
    except Exception as exc:
        if not core_op_declined(exc, "rectangle"):
            raise
        return _rectangle_py(
            dims, size, origin, cell_type, split, order, bool(complete)
        )


# --------------------------------------------------------------------------- #
# box                                                                         #
# --------------------------------------------------------------------------- #


def _box_py(dims, size, origin, cell_type, order, complete=False):
    nx, ny, nz = (int(v) for v in dims)
    px, py, pz = nx + 1, ny + 1, nz + 1
    i = np.tile(np.arange(px, dtype=np.int64), py * pz)
    j = np.tile(np.repeat(np.arange(py, dtype=np.int64), px), pz)
    k = np.repeat(np.arange(pz, dtype=np.int64), px * py)
    mask = (
        (i == 0) * 1
        + (i == nx) * 2
        + (j == 0) * 4
        + (j == ny) * 8
        + (k == 0) * 16
        + (k == nz) * 32
    ).astype(np.uint16)
    param = np.stack(
        [
            origin[0] + (size[0] * i.astype(np.float64)) / float(nx),
            origin[1] + (size[1] * j.astype(np.float64)) / float(ny),
            origin[2] + (size[2] * k.astype(np.float64)) / float(nz),
        ],
        axis=1,
    )
    soup = _Soup("identity")
    soup.add(param, mask)
    ck, rest = np.divmod(np.arange(nx * ny * nz, dtype=np.int64), nx * ny)
    cj, ci = np.divmod(rest, nx)
    b = (ck * py + cj) * px + ci
    t = b + px * py
    h = np.stack([b, b + 1, b + px + 1, b + px, t, t + 1, t + px + 1, t + px], axis=1)
    if cell_type == "hexahedron":
        soup.block("hexahedron", h)
    elif cell_type == "tetra":
        soup.block("tetra", h[:, _KUHN].reshape(-1, 4))
    else:
        w = np.stack([h[:, [0, 1, 2, 4, 5, 6]], h[:, [0, 2, 3, 4, 6, 7]]], axis=1)
        soup.block("wedge", w.reshape(-1, 6))
    names = [
        ("xmin", 1),
        ("xmax", 2),
        ("ymin", 4),
        ("ymax", 8),
        ("zmin", 16),
        ("zmax", 32),
    ]
    return _finish(soup, order, names, 2, complete=complete, shape="box")


def box(
    dims=(1, 1, 1),
    size=(1.0, 1.0, 1.0),
    origin=(0.0, 0.0, 0.0),
    cell_type="hexahedron",
    order=1,
    max_cells=DEFAULT_MAX_CELLS,
    complete=False,
):
    """A regular 3-D grid of hexahedra, tetrahedra or wedges.

    Tetrahedra are the six Kuhn tetrahedra of each hexahedron around its main
    diagonal, which conform across cells; wedges are two per hexahedron.

    Parameters
    ----------
    dims :
        Cell counts ``(nx, ny, nz)``. A zero on any axis gives an empty mesh.
    size :
        Side lengths, positive on every axis with cells.
    origin :
        The lo corner.
    cell_type :
        ``"hexahedron"``, ``"tetra"`` or ``"wedge"``.
    order :
        1, 3 to 10 for the ``VTK_LAGRANGE_*`` family, or 2 for ``hexahedron20`` / ``tetra10`` / ``wedge15`` cells.
    max_cells :
        Refuse by name above this many cells; non-positive disables the check.
    complete :
        With ``order=2``, the full-Lagrange types (``quad9``, ``hexahedron27``,
        ``wedge18``) instead of the serendipity ones.

    Returns
    -------
    Mesh
        Points x fastest, then y, then z; side regions ``xmin``, ``xmax``,
        ``ymin``, ``ymax``, ``zmin``, ``zmax``.
    """
    shape = "box"
    dims = _ivec(shape, "dims", dims, 3)
    size = _vec(shape, "size", size, 3)
    origin = _vec(shape, "origin", origin, 3)
    for k in range(3):
        _check_count(shape, "cell counts", dims[k])
        if dims[k] > 0:
            _check_positive(shape, "size", size[k])
    _check_order(shape, order)
    _choice(shape, "cell_type", cell_type, ("hexahedron", "tetra", "wedge"))
    if 0 in (dims[0], dims[1], dims[2]):
        return _empty()
    per = {"hexahedron": 1, "tetra": 6, "wedge": 2}[cell_type]
    _budget(shape, int(dims[0]) * int(dims[1]) * int(dims[2]) * per, max_cells)
    try:
        return _core_call(
            "box",
            [int(d) for d in dims],
            [float(v) for v in size],
            [float(v) for v in origin],
            cell_type,
            int(order),
            int(max_cells),
            bool(complete),
        )
    except Exception as exc:
        if not core_op_declined(exc, "box"):
            raise
        return _box_py(dims, size, origin, cell_type, order, bool(complete))


# --------------------------------------------------------------------------- #
# polar and squircle layers (disk, cylinder)                                  #
# --------------------------------------------------------------------------- #


def _squircle_layer(soup, n, z, layer_mask):
    i = np.tile(np.arange(n + 1, dtype=np.int64), n + 1)
    j = np.repeat(np.arange(n + 1, dtype=np.int64), n + 1)
    s = (2 * i - n).astype(np.float64) / float(n)
    t = (2 * j - n).astype(np.float64) / float(n)
    rim = (i == 0) | (i == n) | (j == 0) | (j == n)
    mask = (layer_mask | (rim * 4)).astype(np.uint16)
    return soup.add(np.stack([s, t, np.full(s.size, z)], axis=1), mask)


def _polar_layer(soup, n, z, layer_mask):
    first = soup.mask.size
    rho = [np.zeros(1)]
    ang = [np.zeros(1)]
    mask = [np.array([layer_mask | _DEGENERATE])]
    for k in range(1, n + 1):
        m = np.arange(6 * k, dtype=np.float64)
        rho.append(np.full(6 * k, float(k) / float(n)))
        ang.append((2.0 * _PI * m) / float(6 * k))
        mask.append(np.full(6 * k, layer_mask | (4 if k == n else 0)))
    rho = np.concatenate(rho)
    soup.add(
        np.stack([rho, np.concatenate(ang), np.full(rho.size, z)], axis=1),
        np.concatenate(mask),
    )
    return first


def _polar_nodes(n):
    return 1 + 3 * n * (n + 1)


def _polar_triangles(n, base, sectors=6):
    """Triangles of a polar layer: ring by ring, sector by sector."""
    rows = []
    for k in range(1, n + 1):
        ring0 = base + 1 + (sectors * k * (k - 1)) // 2
        t = np.arange(k, dtype=np.int64)
        for s in range(sectors):
            outer = ring0 + ((s * k + t) % (sectors * k))
            outer1 = ring0 + ((s * k + t + 1) % (sectors * k))
            if k == 1:
                inner = np.full(k, base, dtype=np.int64)
                inner1 = inner
            else:
                ring1 = base + 1 + (sectors * (k - 1) * (k - 2)) // 2
                inner = ring1 + ((s * (k - 1) + t) % (sectors * (k - 1)))
                inner1 = ring1 + ((s * (k - 1) + t + 1) % (sectors * (k - 1)))
            a = np.stack([inner, outer, outer1], axis=1)
            b = np.stack([inner[:-1], outer1[:-1], inner1[:-1]], axis=1)
            out = np.empty((2 * k - 1, 3), dtype=np.int64)
            out[0::2] = a
            out[1::2] = b
            rows.append(out)
    return np.concatenate(rows) if rows else np.zeros((0, 3), dtype=np.int64)


# --------------------------------------------------------------------------- #
# disk                                                                        #
# --------------------------------------------------------------------------- #


def _disk_py(radius, n, center, cell_type, order, complete=False):
    if cell_type == "triangle":
        soup = _Soup("polar", radius, center)
        base = _polar_layer(soup, n, 0.0, 0)
        soup.block("triangle", _polar_triangles(n, base))
    else:
        soup = _Soup("squircle", radius, center)
        _squircle_layer(soup, n, 0.0, 0)
        p = n + 1
        cj, ci = np.divmod(np.arange(n * n, dtype=np.int64), n)
        a = cj * p + ci
        soup.block("quad", np.stack([a, a + 1, a + p + 1, a + p], axis=1))
    return _finish(soup, order, [("rim", 4)], 1, complete=complete, shape="disk")


def disk(
    radius=1.0,
    n=8,
    center=(0.0, 0.0, 0.0),
    cell_type="quad",
    order=1,
    max_cells=DEFAULT_MAX_CELLS,
    complete=False,
):
    """A disk in the ``z = center[2]`` plane.

    Parameters
    ----------
    radius :
        The radius, positive.
    n :
        Quads: cells across a diameter (an ``n x n`` grid pushed onto the disk by
        the squircle map ``x' = r x sqrt(1 - y^2/2)``). Triangles: rings between
        the centre and the rim (``6 n^2`` triangles).
    center :
        The centre.
    cell_type :
        ``"quad"`` or ``"triangle"``.
    order :
        1, 3 to 10 for the ``VTK_LAGRANGE_*`` family, or 2 for ``quad8`` / ``triangle6`` cells with a curved rim.
    max_cells :
        Refuse by name above this many cells; non-positive disables the check.
    complete :
        With ``order=2``, the full-Lagrange types (``quad9``, ``hexahedron27``,
        ``wedge18``) instead of the serendipity ones.

    Returns
    -------
    Mesh
        One side region, ``rim``.
    """
    shape = "disk"
    center = _vec(shape, "center", center, 3)
    _check_positive(shape, "radius", float(radius))
    _check_count(shape, "n", int(n))
    _check_order(shape, order)
    _choice(shape, "cell_type", cell_type, ("quad", "triangle"))
    if n == 0:
        return _empty()
    _budget(shape, int(n) * int(n) * (6 if cell_type == "triangle" else 1), max_cells)
    try:
        return _core_call(
            "disk",
            float(radius),
            int(n),
            [float(v) for v in center],
            cell_type,
            int(order),
            int(max_cells),
            bool(complete),
        )
    except Exception as exc:
        if not core_op_declined(exc, "disk"):
            raise
        return _disk_py(float(radius), int(n), center, cell_type, order, bool(complete))


# --------------------------------------------------------------------------- #
# uv_sphere                                                                   #
# --------------------------------------------------------------------------- #


def _uv_sphere_py(radius, nu, nv, center, cell_type, order, complete=False):
    soup = _Soup("sphere", radius, center)
    soup.add([[0.0, 0.0, 0.0]], [_DEGENERATE])
    j = np.repeat(np.arange(1, nv, dtype=np.int64), nu)
    i = np.tile(np.arange(nu, dtype=np.int64), nv - 1)
    soup.add(
        np.stack(
            [
                (_PI * j.astype(np.float64)) / float(nv),
                (2.0 * _PI * i.astype(np.float64)) / float(nu),
                np.zeros(j.size),
            ],
            axis=1,
        ),
        np.zeros(j.size),
    )
    south = soup.add([[_PI, 0.0, 0.0]], [_DEGENERATE])

    def ring(jj, ii):
        return 1 + (jj - 1) * nu + (ii % nu)

    quad = cell_type == "quad"
    ii = np.arange(nu, dtype=np.int64)
    bj = np.repeat(np.arange(1, nv - 1, dtype=np.int64), nu)
    bi = np.tile(ii, nv - 2)
    a, b = ring(bj, bi), ring(bj, bi + 1)
    c, d = ring(bj + 1, bi + 1), ring(bj + 1, bi)
    if quad:
        soup.block("quad", np.stack([a, d, c, b], axis=1))
    north = np.stack(
        [np.zeros(nu, dtype=np.int64), ring(1, ii), ring(1, ii + 1)], axis=1
    )
    pole = np.stack(
        [np.full(nu, south, dtype=np.int64), ring(nv - 1, ii + 1), ring(nv - 1, ii)],
        axis=1,
    )
    parts = [north]
    if not quad:
        body = np.stack([np.stack([a, d, c], 1), np.stack([a, c, b], 1)], axis=1)
        parts.append(body.reshape(-1, 3))
    parts.append(pole)
    soup.block("triangle", np.concatenate(parts))
    return _finish(soup, order, [], 2, complete=complete, shape="uv_sphere")


def uv_sphere(
    radius=1.0,
    nu=16,
    nv=8,
    center=(0.0, 0.0, 0.0),
    cell_type="triangle",
    order=1,
    max_cells=DEFAULT_MAX_CELLS,
    complete=False,
):
    """A sphere surface in latitude / longitude rings, oriented outward.

    Point 0 is the north pole ``(0, 0, +r)``, then ``nv - 1`` rings of ``nu``
    points from north to south, then the south pole.

    Parameters
    ----------
    radius :
        The radius, positive.
    nu :
        Points per ring (at least 3).
    nv :
        Latitude bands (at least 2).
    center :
        The centre.
    cell_type :
        ``"triangle"``, or ``"quad"`` for quads between the rings with triangle
        caps.
    order :
        1, 3 to 10 for the ``VTK_LAGRANGE_*`` family, or 2 for ``triangle6`` / ``quad8`` cells on the exact sphere.
    max_cells :
        Refuse by name above this many cells; non-positive disables the check.
    complete :
        With ``order=2``, the full-Lagrange types (``quad9``, ``hexahedron27``,
        ``wedge18``) instead of the serendipity ones.
    """
    shape = "uv_sphere"
    center = _vec(shape, "center", center, 3)
    _check_positive(shape, "radius", float(radius))
    _check_count(shape, "nu", int(nu))
    _check_count(shape, "nv", int(nv))
    _check_order(shape, order)
    _choice(shape, "cell_type", cell_type, ("triangle", "quad"))
    if nu == 0 or nv == 0:
        return _empty()
    if nu < 3:
        _fail(shape, f"nu must be at least 3, got {int(nu)}")
    if nv < 2:
        _fail(shape, f"nv must be at least 2, got {int(nv)}")
    cells = int(nu) * int(nv) if cell_type == "quad" else 2 * int(nu) * (int(nv) - 1)
    _budget(shape, cells, max_cells)
    try:
        return _core_call(
            "uv_sphere",
            float(radius),
            int(nu),
            int(nv),
            [float(v) for v in center],
            cell_type,
            int(order),
            int(max_cells),
            bool(complete),
        )
    except Exception as exc:
        if not core_op_declined(exc, "uv_sphere"):
            raise
        return _uv_sphere_py(
            float(radius), int(nu), int(nv), center, cell_type, order, bool(complete)
        )


# --------------------------------------------------------------------------- #
# icosa_sphere                                                                #
# --------------------------------------------------------------------------- #

_ICOSA_FACES = np.array(
    [
        [0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11],
        [1, 5, 9], [5, 11, 4], [11, 10, 2], [10, 7, 6], [7, 1, 8],
        [3, 9, 4], [3, 4, 2], [3, 2, 6], [3, 6, 8], [3, 8, 9],
        [4, 9, 5], [2, 4, 11], [6, 2, 10], [8, 6, 7], [9, 8, 1],
    ],
    dtype=np.int64,
)  # fmt: skip


def _icosa_sphere_py(radius, subdivisions, center, order, complete=False):
    soup = _Soup("unit", radius, center)
    t = (1.0 + np.sqrt(5.0)) / 2.0
    length = np.sqrt(1.0 + t * t)
    one, tau = 1.0 / length, t / length
    verts = np.array(
        [
            [-one, tau, 0], [one, tau, 0], [-one, -tau, 0], [one, -tau, 0],
            [0, -one, tau], [0, one, tau], [0, -one, -tau], [0, one, -tau],
            [tau, 0, -one], [tau, 0, one], [-tau, 0, -one], [-tau, 0, one],
        ],
        dtype=np.float64,
    )  # fmt: skip
    soup.add(verts, np.zeros(12))
    faces = _ICOSA_FACES
    for _ in range(subdivisions):
        v0, v1, v2 = faces[:, 0], faces[:, 1], faces[:, 2]
        a = np.stack([v0, v1, v2], axis=1).reshape(-1)
        b = np.stack([v1, v2, v0], axis=1).reshape(-1)
        lo, hi = np.minimum(a, b), np.maximum(a, b)
        ulo, uhi, ids = _dedup_edges(lo, hi, soup.mask.size)
        first = soup.add(_mid(soup, ulo, uhi), np.zeros(ulo.size))
        m = (first + ids).reshape(-1, 3)
        m01, m12, m20 = m[:, 0], m[:, 1], m[:, 2]
        faces = np.stack(
            [
                np.stack([v0, m01, m20], 1),
                np.stack([v1, m12, m01], 1),
                np.stack([v2, m20, m12], 1),
                np.stack([m01, m12, m20], 1),
            ],
            axis=1,
        ).reshape(-1, 3)
    soup.block("triangle", faces)
    return _finish(soup, order, [], 2, complete=complete, shape="icosa_sphere")


def icosa_sphere(
    radius=1.0,
    subdivisions=2,
    center=(0.0, 0.0, 0.0),
    order=1,
    max_cells=DEFAULT_MAX_CELLS,
    complete=False,
):
    """A sphere surface from an icosahedron split in four ``subdivisions`` times.

    ``20 * 4**subdivisions`` triangles, oriented outward; every edge midpoint is
    pushed onto the sphere.

    Parameters
    ----------
    radius :
        The radius, positive.
    subdivisions :
        The number of 1-to-4 splits (zero is the icosahedron).
    center :
        The centre.
    order :
        1, 3 to 10 for the ``VTK_LAGRANGE_*`` family, or 2 for ``triangle6`` cells on the exact sphere.
    max_cells :
        Refuse by name above this many cells; non-positive disables the check.
    complete :
        With ``order=2``, the full-Lagrange types (``quad9``, ``hexahedron27``,
        ``wedge18``) instead of the serendipity ones.
    """
    shape = "icosa_sphere"
    center = _vec(shape, "center", center, 3)
    _check_positive(shape, "radius", float(radius))
    _check_count(shape, "subdivisions", int(subdivisions))
    _check_order(shape, order)
    cells = 20
    for _ in range(int(subdivisions)):
        cells *= 4
        if cells > (1 << 62):
            break
    _budget(shape, cells, max_cells)
    try:
        return _core_call(
            "icosa_sphere",
            float(radius),
            int(subdivisions),
            [float(v) for v in center],
            int(order),
            int(max_cells),
            bool(complete),
        )
    except Exception as exc:
        if not core_op_declined(exc, "icosa_sphere"):
            raise
        return _icosa_sphere_py(
            float(radius), int(subdivisions), center, order, bool(complete)
        )


# --------------------------------------------------------------------------- #
# cylinder                                                                    #
# --------------------------------------------------------------------------- #


def _cylinder_py(radius, height, nxy, nz, center, cell_type, order, complete=False):
    wedge = cell_type == "wedge"
    soup = _Soup("polar" if wedge else "squircle", radius, center)
    layer = _polar_nodes(nxy) if wedge else (nxy + 1) * (nxy + 1)
    for k in range(nz + 1):
        z = (height * float(2 * k - nz)) / float(2 * nz)
        layer_mask = (2 if k == 0 else 0) | (1 if k == nz else 0)
        if wedge:
            _polar_layer(soup, nxy, z, layer_mask)
        else:
            _squircle_layer(soup, nxy, z, layer_mask)
    rows = []
    for k in range(nz):
        base = k * layer
        if wedge:
            tri = _polar_triangles(nxy, base)
            rows.append(np.concatenate([tri, tri + layer], axis=1))
        else:
            p = nxy + 1
            cj, ci = np.divmod(np.arange(nxy * nxy, dtype=np.int64), nxy)
            a = base + cj * p + ci
            t = a + layer
            rows.append(
                np.stack([a, a + 1, a + p + 1, a + p, t, t + 1, t + p + 1, t + p], 1)
            )
    soup.block("wedge" if wedge else "hexahedron", np.concatenate(rows))
    return _finish(
        soup,
        order,
        [("side", 4), ("top", 1), ("bottom", 2)],
        2,
        complete=complete,
        shape="cylinder",
    )


def cylinder(
    radius=1.0,
    height=1.0,
    nxy=8,
    nz=8,
    center=(0.0, 0.0, 0.0),
    cell_type="hexahedron",
    order=1,
    max_cells=DEFAULT_MAX_CELLS,
    complete=False,
):
    """A solid cylinder along z.

    Parameters
    ----------
    radius, height :
        Positive.
    nxy :
        Cells across a diameter (see :func:`disk`).
    nz :
        Layers along the axis.
    center :
        The centre of the cylinder.
    cell_type :
        ``"hexahedron"`` (the squircle disk, extruded) or ``"wedge"`` (the
        triangle disk, extruded).
    order :
        1, 3 to 10 for the ``VTK_LAGRANGE_*`` family, or 2 for ``hexahedron20`` / ``wedge15`` cells with a curved side.
    max_cells :
        Refuse by name above this many cells; non-positive disables the check.
    complete :
        With ``order=2``, the full-Lagrange types (``quad9``, ``hexahedron27``,
        ``wedge18``) instead of the serendipity ones.

    Returns
    -------
    Mesh
        Side regions ``side``, ``top``, ``bottom``.
    """
    shape = "cylinder"
    center = _vec(shape, "center", center, 3)
    _check_positive(shape, "radius", float(radius))
    _check_positive(shape, "height", float(height))
    _check_count(shape, "nxy", int(nxy))
    _check_count(shape, "nz", int(nz))
    _check_order(shape, order)
    _choice(shape, "cell_type", cell_type, ("hexahedron", "wedge"))
    if nxy == 0 or nz == 0:
        return _empty()
    _budget(
        shape,
        int(nxy) * int(nxy) * int(nz) * (6 if cell_type == "wedge" else 1),
        max_cells,
    )
    try:
        return _core_call(
            "cylinder",
            float(radius),
            float(height),
            int(nxy),
            int(nz),
            [float(v) for v in center],
            cell_type,
            int(order),
            int(max_cells),
            bool(complete),
        )
    except Exception as exc:
        if not core_op_declined(exc, "cylinder"):
            raise
        return _cylinder_py(
            float(radius),
            float(height),
            int(nxy),
            int(nz),
            center,
            cell_type,
            order,
            bool(complete),
        )


# --------------------------------------------------------------------------- #
# planar domains                                                              #
# --------------------------------------------------------------------------- #


def _quad_block(soup, tri, a, b, c, d):
    """A quad block, or each quad as the triangles (a b c) and (a c d)."""
    if not tri:
        soup.block("quad", np.stack([a, b, c, d], axis=1))
        return
    first = np.stack([a, b, c], axis=1)
    second = np.stack([a, c, d], axis=1)
    soup.block("triangle", np.stack([first, second], axis=1).reshape(-1, 3))


# --------------------------------------------------------------------------- #
# reference_triangle                                                          #
# --------------------------------------------------------------------------- #


def _reference_triangle_py(n, size, origin, order, complete=False):
    rows_i = [np.arange(n + 1 - j, dtype=np.int64) for j in range(n + 1)]
    i = np.concatenate(rows_i)
    j = np.concatenate([np.full(n + 1 - jj, jj, dtype=np.int64) for jj in range(n + 1)])
    mask = ((j == 0) * 1 + (i == 0) * 2 + (i + j == n) * 4).astype(np.uint16)
    soup = _Soup("identity")
    soup.add(
        np.stack(
            [
                origin[0] + (size[0] * i.astype(np.float64)) / float(n),
                origin[1] + (size[1] * j.astype(np.float64)) / float(n),
                np.full(i.size, origin[2]),
            ],
            axis=1,
        ),
        mask,
    )

    def ident(ii, jj):
        return jj * (n + 1) - (jj * (jj - 1)) // 2 + ii

    rows = []
    for jj in range(n):
        ii = np.arange(n - jj, dtype=np.int64)
        up = np.stack([ident(ii, jj), ident(ii + 1, jj), ident(ii, jj + 1)], axis=1)
        down = np.stack(
            [
                ident(ii[:-1] + 1, jj),
                ident(ii[:-1] + 1, jj + 1),
                ident(ii[:-1], jj + 1),
            ],
            axis=1,
        )
        out = np.empty((2 * ii.size - 1, 3), dtype=np.int64)
        out[0::2] = up
        out[1::2] = down
        rows.append(out)
    soup.block("triangle", np.concatenate(rows))
    return _finish(
        soup,
        order,
        [("ymin", 1), ("xmin", 2), ("hypotenuse", 4)],
        1,
        complete=complete,
        shape="reference_triangle",
    )


def reference_triangle(
    n=4,
    size=(1.0, 1.0),
    origin=(0.0, 0.0, 0.0),
    order=1,
    max_cells=DEFAULT_MAX_CELLS,
    complete=False,
):
    """The reference triangle ``(0,0) (1,0) (0,1)``, scaled and refined into ``n**2`` triangles.

    Nodes run in barycentric lattice order: row ``j`` (y) holds the nodes
    ``i = 0 .. n - j`` (x fastest). Cells run row by row, each lattice triangle
    followed by the inverted one to its right.

    Parameters
    ----------
    n :
        Divisions per edge; zero yields an empty mesh.
    size :
        The leg lengths ``(lx, ly)``, positive.
    origin :
        The right-angle corner.
    order :
        1, 3 to 10 for the ``VTK_LAGRANGE_*`` family, or 2 for ``triangle6`` cells.
    max_cells :
        Refuse by name above this many cells; non-positive disables the check.
    complete :
        With ``order=2``, the full-Lagrange types (``quad9``, ``hexahedron27``,
        ``wedge18``) instead of the serendipity ones.

    Returns
    -------
    Mesh
        Side regions ``ymin``, ``xmin`` and ``hypotenuse``.
    """
    shape = "reference_triangle"
    size = _vec(shape, "size", size, 2)
    origin = _vec(shape, "origin", origin, 3)
    _check_count(shape, "n", int(n))
    _check_order(shape, order)
    if n == 0:
        return _empty()
    for k in range(2):
        _check_positive(shape, "size", float(size[k]))
    _budget(shape, int(n) * int(n), max_cells)
    try:
        return _core_call(
            "reference_triangle",
            int(n),
            [float(v) for v in size],
            [float(v) for v in origin],
            int(order),
            int(max_cells),
            bool(complete),
        )
    except Exception as exc:
        if not core_op_declined(exc, "reference_triangle"):
            raise
        return _reference_triangle_py(int(n), size, origin, order, bool(complete))


# --------------------------------------------------------------------------- #
# ngon                                                                        #
# --------------------------------------------------------------------------- #


def _ngon_py(sides, radius, rings, center, order, complete=False):
    ang = (2.0 * _PI * np.arange(sides, dtype=np.float64)) / float(sides)
    vx = radius * np.cos(ang)
    vy = radius * np.sin(ang)
    px = [np.zeros(1)]
    py = [np.zeros(1)]
    mask = [np.zeros(1)]
    for k in range(1, rings + 1):
        s = np.repeat(np.arange(sides), k)
        t = np.tile(np.arange(k), sides)
        w0 = (k - t).astype(np.float64)
        w1 = t.astype(np.float64)
        nxt = (s + 1) % sides
        px.append((w0 * vx[s] + w1 * vx[nxt]) / float(rings))
        py.append((w0 * vy[s] + w1 * vy[nxt]) / float(rings))
        mask.append(np.full(sides * k, 1 if k == rings else 0))
    px = np.concatenate(px)
    py = np.concatenate(py)
    soup = _Soup("identity")
    soup.add(
        np.stack([center[0] + px, center[1] + py, np.full(px.size, center[2])], axis=1),
        np.concatenate(mask),
    )
    soup.block("triangle", _polar_triangles(rings, 0, sides))
    return _finish(soup, order, [("rim", 1)], 1, complete=complete, shape="ngon")


def ngon(
    sides=6,
    radius=1.0,
    rings=1,
    center=(0.0, 0.0, 0.0),
    order=1,
    max_cells=DEFAULT_MAX_CELLS,
    complete=False,
):
    """A regular polygon fanned from its centre, or refined into rings.

    ``rings=1`` is the plain fan of ``sides`` triangles. More rings give ring
    ``k`` of ``sides * k`` nodes on the straight chords between the vertices,
    with the triangles of :func:`disk`.

    Parameters
    ----------
    sides :
        The number of sides, at least 3; zero yields an empty mesh.
    radius :
        The circumradius, positive.
    rings :
        Radial divisions; zero yields an empty mesh.
    center :
        The centre.
    order :
        1, 3 to 10 for the ``VTK_LAGRANGE_*`` family, or 2 for ``triangle6`` cells (the sides are straight).
    max_cells :
        Refuse by name above this many cells; non-positive disables the check.
    complete :
        With ``order=2``, the full-Lagrange types (``quad9``, ``hexahedron27``,
        ``wedge18``) instead of the serendipity ones.

    Returns
    -------
    Mesh
        One side region, ``rim``.
    """
    shape = "ngon"
    center = _vec(shape, "center", center, 3)
    _check_count(shape, "sides", int(sides))
    _check_count(shape, "rings", int(rings))
    _check_order(shape, order)
    if sides == 0 or rings == 0:
        return _empty()
    if sides < 3:
        _fail(shape, f"sides must be at least 3, got {int(sides)}")
    _check_positive(shape, "radius", float(radius))
    _budget(shape, int(sides) * int(rings) * int(rings), max_cells)
    try:
        return _core_call(
            "ngon",
            int(sides),
            float(radius),
            int(rings),
            [float(v) for v in center],
            int(order),
            int(max_cells),
            bool(complete),
        )
    except Exception as exc:
        if not core_op_declined(exc, "ngon"):
            raise
        return _ngon_py(
            int(sides), float(radius), int(rings), center, order, bool(complete)
        )


# --------------------------------------------------------------------------- #
# annulus                                                                     #
# --------------------------------------------------------------------------- #


def _annulus_py(
    inner, outer, ntheta, nr, angle, center, cell_type, order, complete=False
):
    full = angle >= 2.0 * _PI - 1e-12
    np_ = ntheta if full else ntheta + 1
    k = np.repeat(np.arange(nr + 1, dtype=np.int64), np_)
    i = np.tile(np.arange(np_, dtype=np.int64), nr + 1)
    rho = np.where(
        k == nr,
        outer,
        inner + ((outer - inner) * k.astype(np.float64)) / float(nr),
    )
    theta = (angle * i.astype(np.float64)) / float(ntheta)
    mask = (k == 0) * 1 + (k == nr) * 2
    if not full:
        mask = mask + (i == 0) * 4 + (i == ntheta) * 8
    soup = _Soup("polar", 1.0, center)
    soup.add(np.stack([rho, theta, np.zeros(rho.size)], axis=1), mask.astype(np.uint16))
    ck, ci = np.divmod(np.arange(nr * ntheta, dtype=np.int64), ntheta)
    i1 = (ci + 1) % np_
    _quad_block(
        soup,
        cell_type == "triangle",
        ck * np_ + ci,
        (ck + 1) * np_ + ci,
        (ck + 1) * np_ + i1,
        ck * np_ + i1,
    )
    specs = [("inner", 1), ("outer", 2)]
    if not full:
        specs += [("start", 4), ("end", 8)]
    return _finish(soup, order, specs, 1, complete=complete, shape="annulus")


def annulus(
    inner_radius=0.5,
    outer_radius=1.0,
    ntheta=16,
    nr=4,
    angle=2.0 * _PI,
    center=(0.0, 0.0, 0.0),
    cell_type="quad",
    order=1,
    max_cells=DEFAULT_MAX_CELLS,
    complete=False,
):
    """An annulus, or a sector of one, in the ``z = center[2]`` plane.

    Nodes run ring by ring from the inner radius out, each ring counter-clockwise
    from angle 0. A full ring (``angle = 2 pi``) welds its last node onto its
    first; a sector keeps both end rays as distinct node chains, the seam pair
    that :func:`meshioplusplus.match_periodic_nodes` links.

    Parameters
    ----------
    inner_radius, outer_radius :
        The radii, ``0 < inner_radius < outer_radius``.
    ntheta :
        Cells around; at least 3 for a full ring.
    nr :
        Cells across the radial gap.
    angle :
        The opening angle in ``(0, 2 pi]``; each cell must span less than ``pi``.
    center :
        The centre.
    cell_type :
        ``"quad"`` or ``"triangle"`` (each quad cut along its first diagonal).
    order :
        1, 3 to 10 for the ``VTK_LAGRANGE_*`` family, or 2 for ``quad8`` / ``triangle6`` cells with curved rims.
    max_cells :
        Refuse by name above this many cells; non-positive disables the check.
    complete :
        With ``order=2``, the full-Lagrange types (``quad9``, ``hexahedron27``,
        ``wedge18``) instead of the serendipity ones.

    Returns
    -------
    Mesh
        Side regions ``inner`` and ``outer``, plus ``start`` and ``end`` for a
        sector.
    """
    shape = "annulus"
    center = _vec(shape, "center", center, 3)
    _check_count(shape, "ntheta", int(ntheta))
    _check_count(shape, "nr", int(nr))
    _check_order(shape, order)
    _choice(shape, "cell_type", cell_type, ("quad", "triangle"))
    if ntheta == 0 or nr == 0:
        return _empty()
    inner_radius = float(inner_radius)
    outer_radius = float(outer_radius)
    angle = float(angle)
    _check_positive(shape, "inner_radius", inner_radius)
    if not outer_radius > inner_radius:
        _fail(
            shape,
            f"outer_radius must exceed inner_radius, got {_fmt(outer_radius)} "
            f"and {_fmt(inner_radius)}",
        )
    two_pi = 2.0 * _PI
    if not angle > 0.0 or angle > two_pi + 1e-12:
        _fail(shape, f"angle must be in (0, 2*pi], got {_fmt(angle)}")
    full = angle >= two_pi - 1e-12
    if full and ntheta < 3:
        _fail(shape, f"ntheta must be at least 3 for a full ring, got {int(ntheta)}")
    if not full and angle / float(ntheta) >= _PI:
        _fail(
            shape,
            "ntheta is too small: each cell must span less than pi, got angle "
            f"{_fmt(angle)} over {int(ntheta)} cells",
        )
    _budget(
        shape,
        int(ntheta) * int(nr) * (2 if cell_type == "triangle" else 1),
        max_cells,
    )
    try:
        return _core_call(
            "annulus",
            inner_radius,
            outer_radius,
            int(ntheta),
            int(nr),
            angle,
            [float(v) for v in center],
            cell_type,
            int(order),
            int(max_cells),
            bool(complete),
        )
    except Exception as exc:
        if not core_op_declined(exc, "annulus"):
            raise
        return _annulus_py(
            inner_radius,
            outer_radius,
            int(ntheta),
            int(nr),
            angle,
            center,
            cell_type,
            order,
            bool(complete),
        )


# --------------------------------------------------------------------------- #
# lshape                                                                      #
# --------------------------------------------------------------------------- #


def _lshape_py(n, size, center, cell_type, order, complete=False):
    w = 2 * n + 1
    j, i = np.divmod(np.arange(w * w, dtype=np.int64), w)
    keep = ~((i > n) & (j < n))
    i, j = i[keep], j[keep]
    mask = (
        (i == 0) * 1
        + ((i == 2 * n) & (j >= n)) * 2
        + ((j == 0) & (i <= n)) * 4
        + (j == 2 * n) * 8
        + ((i == n) & (j <= n)) * 16
        + ((j == n) & (i >= n)) * 32
    ).astype(np.uint16)
    ids = np.full(w * w, -1, dtype=np.int64)
    ids[j * w + i] = np.arange(i.size, dtype=np.int64)
    soup = _Soup("identity")
    soup.add(
        np.stack(
            [
                center[0] + (size * (i - n).astype(np.float64)) / float(2 * n),
                center[1] + (size * (j - n).astype(np.float64)) / float(2 * n),
                np.full(i.size, center[2]),
            ],
            axis=1,
        ),
        mask,
    )
    cj, ci = np.divmod(np.arange(4 * n * n, dtype=np.int64), 2 * n)
    live = ~((ci >= n) & (cj < n))
    ci, cj = ci[live], cj[live]
    _quad_block(
        soup,
        cell_type == "triangle",
        ids[cj * w + ci],
        ids[cj * w + ci + 1],
        ids[(cj + 1) * w + ci + 1],
        ids[(cj + 1) * w + ci],
    )
    return _finish(
        soup,
        order,
        [
            ("xmin", 1),
            ("xmax", 2),
            ("ymin", 4),
            ("ymax", 8),
            ("notch_x", 16),
            ("notch_y", 32),
        ],
        1,
        complete=complete,
        shape="lshape",
    )


def lshape(
    n=4,
    size=2.0,
    center=(0.0, 0.0, 0.0),
    cell_type="quad",
    order=1,
    max_cells=DEFAULT_MAX_CELLS,
    complete=False,
):
    """The square ``[-size/2, size/2]**2`` minus its fourth quadrant.

    The removed quadrant is ``x > 0``, ``y < 0``; the re-entrant corner is a node
    at ``center``, the singular point that refinement and error estimators probe.

    Parameters
    ----------
    n :
        Cells along half a side; the mesh has ``3 n**2`` quads.
    size :
        The side of the bounding square, positive.
    center :
        The centre of the bounding square.
    cell_type :
        ``"quad"`` or ``"triangle"`` (each quad cut along its first diagonal).
    order :
        1, 3 to 10 for the ``VTK_LAGRANGE_*`` family, or 2 for ``quad8`` / ``triangle6`` cells.
    max_cells :
        Refuse by name above this many cells; non-positive disables the check.
    complete :
        With ``order=2``, the full-Lagrange types (``quad9``, ``hexahedron27``,
        ``wedge18``) instead of the serendipity ones.

    Returns
    -------
    Mesh
        Side regions ``xmin``, ``xmax``, ``ymin``, ``ymax``, ``notch_x`` (the notch
        wall on ``x = 0``) and ``notch_y`` (the one on ``y = 0``).
    """
    shape = "lshape"
    center = _vec(shape, "center", center, 3)
    _check_count(shape, "n", int(n))
    _check_order(shape, order)
    _choice(shape, "cell_type", cell_type, ("quad", "triangle"))
    if n == 0:
        return _empty()
    _check_positive(shape, "size", float(size))
    _budget(
        shape, 3 * int(n) * int(n) * (2 if cell_type == "triangle" else 1), max_cells
    )
    try:
        return _core_call(
            "lshape",
            int(n),
            float(size),
            [float(v) for v in center],
            cell_type,
            int(order),
            int(max_cells),
            bool(complete),
        )
    except Exception as exc:
        if not core_op_declined(exc, "lshape"):
            raise
        return _lshape_py(int(n), float(size), center, cell_type, order, bool(complete))


# --------------------------------------------------------------------------- #
# plate_with_hole                                                             #
# --------------------------------------------------------------------------- #

_PLATE_SIDE_BIT = np.array([4, 16, 2, 8], dtype=np.int64)  # xmax, ymax, xmin, ymin


def _plate_py(radius, size, nt, nr, center, cell_type, order, complete=False):
    m = 4 * nt
    k = np.repeat(np.arange(nr + 1, dtype=np.int64), m)
    p = np.tile(np.arange(m, dtype=np.int64), nr + 1)
    b = p // nt
    mask = (k == 0) * 1
    outer = k == nr
    mask = mask + outer * _PLATE_SIDE_BIT[b]
    mask = mask + (outer & (p % nt == 0)) * _PLATE_SIDE_BIT[(b + 3) % 4]
    soup = _Soup("plate", radius, center)
    soup.half = size / 2.0
    soup.period = float(m)
    soup.nt = float(nt)
    soup.add(
        np.stack(
            [k.astype(np.float64) / float(nr), p.astype(np.float64), np.zeros(k.size)],
            axis=1,
        ),
        mask.astype(np.uint16),
    )
    ck, cp = np.divmod(np.arange(nr * m, dtype=np.int64), m)
    p1 = (cp + 1) % m
    _quad_block(
        soup,
        cell_type == "triangle",
        ck * m + cp,
        (ck + 1) * m + cp,
        (ck + 1) * m + p1,
        ck * m + p1,
    )
    return _finish(
        soup,
        order,
        [("hole", 1), ("xmin", 2), ("xmax", 4), ("ymin", 8), ("ymax", 16)],
        1,
        complete=complete,
        shape="plate_with_hole",
    )


def plate_with_hole(
    radius=0.25,
    size=1.0,
    nt=8,
    nr=4,
    center=(0.0, 0.0, 0.0),
    cell_type="quad",
    order=1,
    max_cells=DEFAULT_MAX_CELLS,
    complete=False,
):
    """A square plate around a circular hole, in four quad blocks.

    Nodes run ring by ring from the hole out (index ``k``), ``4 nt`` positions
    around each (index ``p``, counter-clockwise from the south-east corner). A
    node blends the hole circle into the square linearly in ``k``.

    Parameters
    ----------
    radius :
        The hole radius, positive and below half of ``size``.
    size :
        The side of the plate.
    nt :
        Cells around the hole per 90 degree block.
    nr :
        Cells from the hole to the outer square.
    center :
        The centre of the hole.
    cell_type :
        ``"quad"`` or ``"triangle"`` (each quad cut along its first diagonal).
    order :
        1, 3 to 10 for the ``VTK_LAGRANGE_*`` family, or 2 for ``quad8`` / ``triangle6`` cells with a circular hole.
    max_cells :
        Refuse by name above this many cells; non-positive disables the check.
    complete :
        With ``order=2``, the full-Lagrange types (``quad9``, ``hexahedron27``,
        ``wedge18``) instead of the serendipity ones.

    Returns
    -------
    Mesh
        Side regions ``hole``, ``xmin``, ``xmax``, ``ymin`` and ``ymax``.
    """
    shape = "plate_with_hole"
    center = _vec(shape, "center", center, 3)
    _check_count(shape, "nt", int(nt))
    _check_count(shape, "nr", int(nr))
    _check_order(shape, order)
    _choice(shape, "cell_type", cell_type, ("quad", "triangle"))
    if nt == 0 or nr == 0:
        return _empty()
    radius = float(radius)
    size = float(size)
    _check_positive(shape, "radius", radius)
    _check_positive(shape, "size", size)
    if not radius < size / 2.0:
        _fail(
            shape,
            "radius must be less than half the size, got radius "
            f"{_fmt(radius)} and size {_fmt(size)}",
        )
    _budget(
        shape,
        4 * int(nt) * int(nr) * (2 if cell_type == "triangle" else 1),
        max_cells,
    )
    try:
        return _core_call(
            "plate_with_hole",
            radius,
            size,
            int(nt),
            int(nr),
            [float(v) for v in center],
            cell_type,
            int(order),
            int(max_cells),
            bool(complete),
        )
    except Exception as exc:
        if not core_op_declined(exc, "plate_with_hole"):
            raise
        return _plate_py(
            radius, size, int(nt), int(nr), center, cell_type, order, bool(complete)
        )


# --------------------------------------------------------------------------- #
# honeycomb                                                                   #
# --------------------------------------------------------------------------- #

_HEX_CORNERS = np.array(
    [[0, -2], [1, -1], [1, 1], [0, 2], [-1, 1], [-1, -1]], dtype=np.int64
)


def _honeycomb_py(dims, edge, origin):
    nx, ny = int(dims[0]), int(dims[1])
    wu = 2 * nx + 2
    cj, ci = np.divmod(np.arange(nx * ny, dtype=np.int64), nx)
    cu = 2 * ci + 1 + (cj & 1)
    cv = 3 * cj + 2
    u = cu[:, None] + _HEX_CORNERS[None, :, 0]
    v = cv[:, None] + _HEX_CORNERS[None, :, 1]
    key = (v * wu + u).reshape(-1)
    uniq, first, inverse = np.unique(key, return_index=True, return_inverse=True)
    order = np.argsort(first, kind="stable")
    rank = np.empty_like(order)
    rank[order] = np.arange(order.size)
    conn = rank[inverse.reshape(-1)].reshape(-1, 6).astype(np.int64)
    kept = uniq[order]
    dx = (np.sqrt(3.0) * edge) / 2.0
    dy = edge / 2.0
    points = np.stack(
        [
            origin[0] + dx * (kept % wu).astype(np.float64),
            origin[1] + dy * (kept // wu).astype(np.float64),
            np.full(kept.size, origin[2]),
        ],
        axis=1,
    )
    # Facet f of a polygon is its edge from node f to node f + 1.
    a = conn
    b = np.roll(conn, -1, axis=1)
    pairs = np.stack([np.minimum(a, b), np.maximum(a, b)], axis=-1).reshape(-1, 2)
    _, inv, counts = np.unique(pairs, axis=0, return_inverse=True, return_counts=True)
    once = (counts[inv.reshape(-1)] == 1).reshape(conn.shape)
    cell, facet = np.nonzero(once)
    entries = np.stack([cell, facet], axis=1).astype(np.int64)
    return Mesh(
        points,
        [("polygon", conn)],
        regions=[Region("boundary", "side", entries, dim=1)],
    )


def honeycomb(
    dims=(4, 4), edge=1.0, origin=(0.0, 0.0, 0.0), max_cells=DEFAULT_MAX_CELLS
):
    """A tiling of pointy-top hexagons, as one ``polygon`` cell block.

    Hexagons run row by row (x fastest), odd rows shifted half a hexagon to the
    right. Nodes are numbered by first appearance over the cells, each hexagon's
    corners counter-clockwise from the bottom one. The block is ragged in the
    core (``polygon`` rows of six), so the ragged-polygon paths have a parametric
    fixture.

    Parameters
    ----------
    dims :
        Hexagons per row and rows ``(nx, ny)``; a zero yields an empty mesh.
    edge :
        The side length, positive.
    origin :
        The lower-left corner of the bounding box.
    max_cells :
        Refuse by name above this many cells; non-positive disables the check.

    Returns
    -------
    Mesh
        One side region, ``boundary``: the hexagon edges seen once.
    """
    shape = "honeycomb"
    dims = _ivec(shape, "dims", dims, 2)
    origin = _vec(shape, "origin", origin, 3)
    for k in range(2):
        _check_count(shape, "cell counts", int(dims[k]))
    if dims[0] == 0 or dims[1] == 0:
        return _empty()
    _check_positive(shape, "edge", float(edge))
    _budget(shape, int(dims[0]) * int(dims[1]), max_cells)
    try:
        return _core_call(
            "honeycomb",
            [int(v) for v in dims],
            float(edge),
            [float(v) for v in origin],
            int(max_cells),
        )
    except Exception as exc:
        if not core_op_declined(exc, "honeycomb"):
            raise
        return _honeycomb_py(dims, float(edge), origin)


# --------------------------------------------------------------------------- #
# reference_cells                                                             #
# --------------------------------------------------------------------------- #

_ALL_TYPES = (
    "vertex", "line", "triangle", "quad", "tetra", "hexahedron", "wedge", "pyramid",
    "line3", "triangle6", "quad8", "tetra10", "hexahedron20", "wedge15", "pyramid13",
    "polyhedron",
)  # fmt: skip
_PAIR_TYPES = ("line", "triangle", "quad", "tetra", "hexahedron", "wedge", "pyramid")
_QUADRATIC_BASE = {
    "line3": "line",
    "triangle6": "triangle",
    "quad8": "quad",
    "tetra10": "tetra",
    "hexahedron20": "hexahedron",
    "wedge15": "wedge",
    "pyramid13": "pyramid",
}

_LINEAR = {
    "vertex": ([(0, 0, 0)], [0]),
    "line": ([(0, 0, 0), (1, 0, 0)], [0, 1]),
    "triangle": ([(0, 0, 0), (1, 0, 0), (0, 1, 0)], [0, 1, 2]),
    "quad": ([(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)], [0, 1, 2, 3]),
    "tetra": ([(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)], [0, 1, 2, 3]),
    "hexahedron": (
        [
            (0, 0, 0),
            (1, 0, 0),
            (1, 1, 0),
            (0, 1, 0),
            (0, 0, 1),
            (1, 0, 1),
            (1, 1, 1),
            (0, 1, 1),
        ],
        [0, 1, 2, 3, 4, 5, 6, 7],
    ),
    "wedge": (
        [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (0, 1, 1)],
        [0, 1, 2, 3, 4, 5],
    ),
    "pyramid": (
        [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0.5, 0.5, 1)],
        [0, 1, 2, 3, 4],
    ),
}

_PAIRS = {
    "line": ([(0, 0, 0), (1, 0, 0), (2, 0, 0)], [[0, 1], [1, 2]]),
    "triangle": ([(0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0)], [[0, 1, 2], [1, 3, 2]]),
    "quad": (
        [(0, 0, 0), (1, 0, 0), (2, 0, 0), (0, 1, 0), (1, 1, 0), (2, 1, 0)],
        [[0, 1, 4, 3], [1, 2, 5, 4]],
    ),
    "tetra": (
        [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 1)],
        [[0, 1, 2, 3], [2, 3, 1, 4]],
    ),
    "hexahedron": (
        [
            (0, 0, 0), (1, 0, 0), (2, 0, 0), (0, 1, 0), (1, 1, 0), (2, 1, 0),
            (0, 0, 1), (1, 0, 1), (2, 0, 1), (0, 1, 1), (1, 1, 1), (2, 1, 1),
        ],
        [[0, 1, 4, 3, 6, 7, 10, 9], [1, 2, 5, 4, 7, 8, 11, 10]],
    ),
    "wedge": (
        [
            (0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1),
            (0, 1, 1), (0, 0, 2), (1, 0, 2), (0, 1, 2),
        ],
        [[0, 1, 2, 3, 4, 5], [3, 4, 5, 6, 7, 8]],
    ),
    "pyramid": (
        [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0.5, 0.5, 1), (0.5, 0.5, -1)],
        [[0, 1, 2, 3, 4], [0, 3, 2, 1, 5]],
    ),
}  # fmt: skip

# The hexahedron's six outward faces (the polyhedron block).
_CUBE_FACES = [
    [0, 4, 7, 3],
    [1, 2, 6, 5],
    [0, 1, 5, 4],
    [3, 7, 6, 2],
    [0, 3, 2, 1],
    [4, 5, 6, 7],
]


def _reference_cell(name, pairs):
    if name == "polyhedron":
        pts, row = _LINEAR["hexahedron"]
        return np.array(pts, dtype=np.float64), [row]
    if pairs:
        pts, rows = _PAIRS[name]
        return np.array(pts, dtype=np.float64), rows
    if name in _QUADRATIC_BASE:
        from ._convert_cells import _ELEVATE

        base = _QUADRATIC_BASE[name]
        pts, row = _LINEAR[base]
        pts = [list(map(float, p)) for p in pts]
        row = list(row)
        for e in _ELEVATE[base][2]:
            a, b = pts[row[e[0]]], pts[row[e[1]]]
            pts.append([(a[k] + b[k]) * 0.5 for k in range(3)])
            row.append(len(pts) - 1)
        return np.array(pts, dtype=np.float64), [row]
    pts, row = _LINEAR[name]
    return np.array(pts, dtype=np.float64), [row]


def _reference_cells_py(types, pairs):
    stride = 3.0 if pairs else 2.0
    points, cells = [], []
    offset = 0
    for m, name in enumerate(types):
        pts, rows = _reference_cell(name, pairs)
        pts = pts.copy()
        pts[:, 0] = pts[:, 0] + stride * float(m)
        points.append(pts)
        conn = np.array(rows, dtype=np.int64) + offset
        if name == "polyhedron":
            cube = conn[0]
            cells.append(
                ("polyhedron8", [[[int(cube[v]) for v in f] for f in _CUBE_FACES]])
            )
        else:
            cells.append((name, conn))
        offset += len(pts)
    return Mesh(np.concatenate(points), cells)


def reference_cells(types=None, pairs=False):
    """One cell of each requested type in its reference coordinates.

    The cells are laid out side by side along x without touching.

    Parameters
    ----------
    types :
        meshio type names; ``None`` is every supported type in the order
        ``vertex``, ``line``, ``triangle``, ``quad``, ``tetra``, ``hexahedron``,
        ``wedge``, ``pyramid``, ``line3``, ``triangle6``, ``quad8``, ``tetra10``,
        ``hexahedron20``, ``wedge15``, ``pyramid13``, ``polyhedron``.
    pairs :
        Two cells of each type sharing a facet instead of one cell (linear types
        only: ``line``, ``triangle``, ``quad``, ``tetra``, ``hexahedron``,
        ``wedge``, ``pyramid``).

    Returns
    -------
    Mesh
        One block per type, in the order requested.
    """
    shape = "reference_cells"
    if types is None:
        types = []
    if isinstance(types, str):
        types = [t for t in types.split(",") if t]
    types = [str(t) for t in types]
    if not types:
        types = list(_PAIR_TYPES if pairs else _ALL_TYPES)
    for i, name in enumerate(types):
        if name not in _ALL_TYPES:
            _fail(
                shape,
                f"unknown cell type '{name}' (expected one of: {', '.join(_ALL_TYPES)})",
            )
        if name in types[:i]:
            _fail(shape, f"cell type '{name}' is listed twice")
        if pairs and name not in _PAIR_TYPES:
            _fail(
                shape,
                "pairs=True supports line, triangle, quad, tetra, hexahedron, wedge "
                f"and pyramid, not '{name}'",
            )
    try:
        return _core_call("reference_cells", types, bool(pairs))
    except Exception as exc:
        if not core_op_declined(exc, "reference_cells"):
            raise
        return _reference_cells_py(types, bool(pairs))


# --------------------------------------------------------------------------- #
# shared argument handling of the second tier                                 #
# --------------------------------------------------------------------------- #


def _cell_choice(shape, cell_type, choices):
    _choice(shape, "cell_type", cell_type, choices)
    return cell_type


def _call(op, twin, *args):
    """The core's ``op``, or the twin when the core declines."""
    try:
        return _core_call(op, *args)
    except Exception as exc:
        if not core_op_declined(exc, op):
            raise
        return twin()


def _floats(values):
    return [float(v) for v in values]


# --------------------------------------------------------------------------- #
# line, circle, helix                                                         #
# --------------------------------------------------------------------------- #


def _line_py(n, start, end, order, complete):
    i = np.arange(n + 1, dtype=np.float64)
    param = np.empty((n + 1, 3))
    for k in range(3):
        param[:, k] = start[k] + ((end[k] - start[k]) * i) / float(n)
        param[n, k] = end[k]
    soup = _Soup("identity")
    soup.add(param, np.zeros(n + 1))
    a = np.arange(n, dtype=np.int64)
    soup.block("line", np.stack([a, a + 1], axis=1))
    extra = [_point_region("start", 0), _point_region("end", n)]
    return _finish(soup, order, [], 0, complete, extra, "line")


def line(
    n=8,
    start=(0.0, 0.0, 0.0),
    end=(1.0, 0.0, 0.0),
    order=1,
    max_cells=DEFAULT_MAX_CELLS,
    complete=False,
):
    """A straight segment from ``start`` to ``end`` in ``n`` ``line`` cells.

    Parameters
    ----------
    n :
        The cell count; a zero gives an empty mesh.
    start, end :
        The first and the last point; they must differ.
    order :
        1, 2 for ``line3``, or 3 to 10 for ``VTK_LAGRANGE_CURVE``.
    max_cells :
        Refuse by name above this many cells; non-positive disables the check.
    complete :
        Accepted for symmetry (``line3`` is already complete).

    Returns
    -------
    Mesh
        Points from ``start`` to ``end``; point regions ``start`` and ``end``.
    """
    shape = "line"
    start = _vec(shape, "start", start, 3)
    end = _vec(shape, "end", end, 3)
    _check_count(shape, "n", int(n))
    _check_order(shape, order)
    if n == 0:
        return _empty()
    if np.array_equal(start, end):
        _fail(shape, "start and end must differ")
    _budget(shape, int(n), max_cells)
    return _call(
        "line",
        lambda: _line_py(int(n), start, end, int(order), bool(complete)),
        int(n),
        _floats(start),
        _floats(end),
        int(order),
        int(max_cells),
        bool(complete),
    )


def _circle_py(radius, n, center, order, complete):
    soup = _Soup("polar", radius, center)
    i = np.arange(n, dtype=np.float64)
    soup.add(
        np.stack([np.ones(n), (2.0 * _PI * i) / float(n), np.zeros(n)], axis=1),
        np.zeros(n),
    )
    a = np.arange(n, dtype=np.int64)
    soup.block("line", np.stack([a, (a + 1) % n], axis=1))
    return _finish(soup, order, [], 0, complete, (), "circle")


def circle(
    radius=1.0,
    n=16,
    center=(0.0, 0.0, 0.0),
    order=1,
    max_cells=DEFAULT_MAX_CELLS,
    complete=False,
):
    """A closed polyline of ``n`` ``line`` cells on a circle in the
    ``z = center[2]`` plane, counter-clockwise from +x, the seam welded.

    Parameters
    ----------
    radius :
        The radius, positive.
    n :
        The cell count, at least 3; a zero gives an empty mesh.
    center :
        The centre.
    order :
        1, 2 for ``line3`` on the arc, or 3 to 10 for ``VTK_LAGRANGE_CURVE``.
    max_cells :
        Refuse by name above this many cells; non-positive disables the check.
    complete :
        Accepted for symmetry.
    """
    shape = "circle"
    center = _vec(shape, "center", center, 3)
    _check_positive(shape, "radius", float(radius))
    _check_count(shape, "n", int(n))
    _check_order(shape, order)
    if n == 0:
        return _empty()
    if n < 3:
        _fail(shape, f"n must be at least 3, got {int(n)}")
    _budget(shape, int(n), max_cells)
    return _call(
        "circle",
        lambda: _circle_py(float(radius), int(n), center, int(order), bool(complete)),
        float(radius),
        int(n),
        _floats(center),
        int(order),
        int(max_cells),
        bool(complete),
    )


def _helix_py(radius, pitch, turns, n, center, order, complete):
    soup = _Soup("polar", radius, center)
    i = np.arange(n + 1, dtype=np.float64)
    soup.add(
        np.stack(
            [
                np.ones(n + 1),
                (2.0 * _PI * turns * i) / float(n),
                (pitch * turns * i) / float(n),
            ],
            axis=1,
        ),
        np.zeros(n + 1),
    )
    a = np.arange(n, dtype=np.int64)
    soup.block("line", np.stack([a, a + 1], axis=1))
    extra = [_point_region("start", 0), _point_region("end", n)]
    return _finish(soup, order, [], 0, complete, extra, "helix")


def helix(
    radius=1.0,
    pitch=1.0,
    turns=2.0,
    n=32,
    center=(0.0, 0.0, 0.0),
    order=1,
    max_cells=DEFAULT_MAX_CELLS,
    complete=False,
):
    """A helix (a spring) of ``n`` ``line`` cells around the z axis.

    Parameters
    ----------
    radius :
        The radius, positive.
    pitch :
        The rise per turn, nonzero (negative descends).
    turns :
        The number of turns, positive (need not be whole).
    n :
        The cell count, more than ``2 * turns`` so no cell spans half a turn;
        a zero gives an empty mesh.
    center :
        The point on the axis level with the first node.
    order :
        1, 2 for ``line3`` on the helix, or 3 to 10.
    max_cells :
        Refuse by name above this many cells; non-positive disables the check.
    complete :
        Accepted for symmetry.

    Returns
    -------
    Mesh
        Counter-clockwise from +x; point regions ``start`` and ``end``.
    """
    shape = "helix"
    center = _vec(shape, "center", center, 3)
    _check_positive(shape, "radius", float(radius))
    _check_positive(shape, "turns", float(turns))
    _check_count(shape, "n", int(n))
    _check_order(shape, order)
    if n == 0:
        return _empty()
    if not (float(pitch) != 0.0 and np.isfinite(float(pitch))):
        _fail(shape, f"pitch must be finite and nonzero, got {_fmt(pitch)}")
    if not float(n) > 2.0 * float(turns):
        _fail(
            shape,
            "n must exceed 2 * turns so no cell spans half a turn, got n "
            f"{int(n)} and turns {_fmt(turns)}",
        )
    _budget(shape, int(n), max_cells)
    return _call(
        "helix",
        lambda: _helix_py(
            float(radius),
            float(pitch),
            float(turns),
            int(n),
            center,
            int(order),
            bool(complete),
        ),
        float(radius),
        float(pitch),
        float(turns),
        int(n),
        _floats(center),
        int(order),
        int(max_cells),
        bool(complete),
    )


# --------------------------------------------------------------------------- #
# sphere_points, box_points                                                   #
# --------------------------------------------------------------------------- #

_MASK64 = (1 << 64) - 1
_GAMMA = 0x9E3779B97F4A7C15


def _splitmix_draws(seed, count):
    """``count`` splitmix64 draws in [0, 1) from ``seed``: the core's sequence."""
    k = np.arange(1, count + 1, dtype=np.uint64)
    z = np.uint64(seed & _MASK64) + k * np.uint64(_GAMMA)
    z = (z ^ (z >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
    z = (z ^ (z >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
    z = z ^ (z >> np.uint64(31))
    return (z >> np.uint64(11)).astype(np.float64) * (1.0 / 9007199254740992.0)


def _vertex_mesh(points):
    n = points.shape[0]
    return Mesh(points, [("vertex", np.arange(n, dtype=np.int64).reshape(-1, 1))])


def _sphere_points_py(n, radius, center):
    i = np.arange(n, dtype=np.float64)
    z = 1.0 - (2.0 * i + 1.0) / float(n)
    rho = np.sqrt(1.0 - z * z)
    phi = (_PI * (3.0 - np.sqrt(5.0))) * i
    pts = np.stack(
        [
            center[0] + radius * rho * np.cos(phi),
            center[1] + radius * rho * np.sin(phi),
            center[2] + radius * z,
        ],
        axis=1,
    )
    return _vertex_mesh(pts)


def sphere_points(
    n=100, radius=1.0, center=(0.0, 0.0, 0.0), max_cells=DEFAULT_MAX_CELLS
):
    """``n`` points on a sphere on a Fibonacci lattice, as ``vertex`` cells.

    Point ``i`` sits at height ``z = 1 - (2 i + 1) / n`` and azimuth
    ``i * pi * (3 - sqrt 5)``: near-uniform, without poles or seams.

    Parameters
    ----------
    n :
        The point count; a zero gives an empty mesh.
    radius :
        The radius, positive.
    center :
        The centre.
    max_cells :
        Refuse by name above this many points; non-positive disables the check.
    """
    shape = "sphere_points"
    center = _vec(shape, "center", center, 3)
    _check_count(shape, "n", int(n))
    if n == 0:
        return _empty()
    _check_positive(shape, "radius", float(radius))
    _budget(shape, int(n), max_cells)
    return _call(
        "sphere_points",
        lambda: _sphere_points_py(int(n), float(radius), center),
        int(n),
        float(radius),
        _floats(center),
        int(max_cells),
    )


def _box_points_py(n, size, origin, seed):
    u = _splitmix_draws(seed, 3 * n).reshape(n, 3)
    pts = np.stack([origin[k] + size[k] * u[:, k] for k in range(3)], axis=1)
    return _vertex_mesh(pts)


def box_points(
    n=100,
    size=(1.0, 1.0, 1.0),
    origin=(0.0, 0.0, 0.0),
    seed=0,
    max_cells=DEFAULT_MAX_CELLS,
):
    """``n`` points drawn uniformly in a box, as ``vertex`` cells.

    The draws come from splitmix64 seeded with ``seed``, three per point (x, y,
    z), each ``(next() >> 11) * 2**-53``: the same cloud on every platform, in
    the core and in this twin (``numpy.random`` has no C++ counterpart).

    Parameters
    ----------
    n :
        The point count; a zero gives an empty mesh.
    size :
        The side lengths, positive.
    origin :
        The lo corner.
    seed :
        Any integer (taken modulo 2**64).
    max_cells :
        Refuse by name above this many points; non-positive disables the check.
    """
    shape = "box_points"
    size = _vec(shape, "size", size, 3)
    origin = _vec(shape, "origin", origin, 3)
    _check_count(shape, "n", int(n))
    if n == 0:
        return _empty()
    for k in range(3):
        _check_positive(shape, "size", size[k])
    _budget(shape, int(n), max_cells)
    seed = int(seed) & _MASK64
    return _call(
        "box_points",
        lambda: _box_points_py(int(n), size, origin, seed),
        int(n),
        _floats(size),
        _floats(origin),
        seed,
        int(max_cells),
    )


# --------------------------------------------------------------------------- #
# tetra_sphere, octa_sphere                                                   #
# --------------------------------------------------------------------------- #

_TETRA_FACES = np.array([[0, 1, 2], [0, 2, 3], [0, 3, 1], [1, 3, 2]], dtype=np.int64)
_OCTA_VERTS = np.array(
    [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]],
    dtype=np.float64,
)
_OCTA_FACES = np.array(
    [[0, 2, 4], [0, 3, 5], [0, 4, 3], [0, 5, 2], [1, 2, 5], [1, 3, 4], [1, 4, 2], [1, 5, 3]],
    dtype=np.int64,
)  # fmt: skip


def _tetra_verts():
    a = 1.0 / np.sqrt(3.0)
    return np.array(
        [[a, a, a], [a, -a, -a], [-a, a, -a], [-a, -a, a]], dtype=np.float64
    )


def _subdivided_sphere(soup, faces, subdivisions):
    for _ in range(subdivisions):
        v0, v1, v2 = faces[:, 0], faces[:, 1], faces[:, 2]
        a = np.stack([v0, v1, v2], axis=1).reshape(-1)
        b = np.stack([v1, v2, v0], axis=1).reshape(-1)
        lo, hi = np.minimum(a, b), np.maximum(a, b)
        ulo, uhi, ids = _dedup_edges(lo, hi, soup.mask.size)
        first = soup.add(_mid(soup, ulo, uhi), np.zeros(ulo.size))
        m = (first + ids).reshape(-1, 3)
        m01, m12, m20 = m[:, 0], m[:, 1], m[:, 2]
        faces = np.stack(
            [
                np.stack([v0, m01, m20], 1),
                np.stack([v1, m12, m01], 1),
                np.stack([v2, m20, m12], 1),
                np.stack([m01, m12, m20], 1),
            ],
            axis=1,
        ).reshape(-1, 3)
    return faces


def _seeded_sphere_py(
    shape, verts, faces, radius, subdivisions, center, order, complete
):
    soup = _Soup("unit", radius, center)
    soup.add(verts, np.zeros(verts.shape[0]))
    soup.block("triangle", _subdivided_sphere(soup, faces, subdivisions))
    return _finish(soup, order, [], 2, complete, (), shape)


def _seeded_sphere(
    shape, verts, faces, radius, subdivisions, center, order, max_cells, complete
):
    center = _vec(shape, "center", center, 3)
    _check_positive(shape, "radius", float(radius))
    _check_count(shape, "subdivisions", int(subdivisions))
    _check_order(shape, order)
    cells = faces.shape[0]
    for _ in range(int(subdivisions)):
        cells *= 4
        if cells > (1 << 62):
            break
    _budget(shape, cells, max_cells)
    return _call(
        shape,
        lambda: _seeded_sphere_py(
            shape,
            verts,
            faces,
            float(radius),
            int(subdivisions),
            center,
            int(order),
            bool(complete),
        ),
        float(radius),
        int(subdivisions),
        _floats(center),
        int(order),
        int(max_cells),
        bool(complete),
    )


def tetra_sphere(
    radius=1.0,
    subdivisions=3,
    center=(0.0, 0.0, 0.0),
    order=1,
    max_cells=DEFAULT_MAX_CELLS,
    complete=False,
):
    """A sphere surface from a regular tetrahedron split in four
    ``subdivisions`` times (``4 * 4**subdivisions`` triangles), as
    :func:`icosa_sphere`: the most uneven of the seeds, a mesh-quality probe.

    Parameters are those of :func:`icosa_sphere`, plus ``complete`` (accepted
    for symmetry: ``triangle6`` is already complete).
    """
    return _seeded_sphere(
        "tetra_sphere",
        _tetra_verts(),
        _TETRA_FACES,
        radius,
        subdivisions,
        center,
        order,
        max_cells,
        complete,
    )


def octa_sphere(
    radius=1.0,
    subdivisions=2,
    center=(0.0, 0.0, 0.0),
    order=1,
    max_cells=DEFAULT_MAX_CELLS,
    complete=False,
):
    """A sphere surface from a regular octahedron split in four
    ``subdivisions`` times (``8 * 4**subdivisions`` triangles), as
    :func:`icosa_sphere`; its vertices include the six axis points exactly.

    Parameters are those of :func:`icosa_sphere`, plus ``complete``.
    """
    return _seeded_sphere(
        "octa_sphere",
        _OCTA_VERTS,
        _OCTA_FACES,
        radius,
        subdivisions,
        center,
        order,
        max_cells,
        complete,
    )


# --------------------------------------------------------------------------- #
# the boundary lattice of a box: cube_sphere, box_surface, spherical_shell    #
# --------------------------------------------------------------------------- #


class _CubeLattice:
    """The boundary nodes of an ``nx x ny x nz`` lattice, in lattice order."""

    def __init__(self, nx, ny, nz):
        self.nx, self.ny, self.nz = int(nx), int(ny), int(nz)

    def full(self):
        return (self.nx + 1) * (self.ny + 1)

    def ring(self):
        return 2 * (self.nx + 1) + 2 * (self.ny - 1)

    def count(self):
        return 2 * self.full() + (self.nz - 1) * self.ring()

    def ijk(self):
        """The (i, j, k) of every boundary node, in order."""
        nx, ny, nz = self.nx, self.ny, self.nz
        k, j, i = np.meshgrid(
            np.arange(nz + 1), np.arange(ny + 1), np.arange(nx + 1), indexing="ij"
        )
        i, j, k = i.reshape(-1), j.reshape(-1), k.reshape(-1)
        on = (i == 0) | (i == nx) | (j == 0) | (j == ny) | (k == 0) | (k == nz)
        return i[on].astype(np.int64), j[on].astype(np.int64), k[on].astype(np.int64)

    def id(self, i, j, k):
        nx, ny, nz = self.nx, self.ny, self.nz
        i, j, k = (np.asarray(v, dtype=np.int64) for v in (i, j, k))
        base = self.full() + (k - 1) * self.ring()
        cap = np.where(k == 0, j * (nx + 1) + i, base + j * (nx + 1) + i)
        ring = np.where(
            j == 0,
            base + i,
            np.where(
                j == ny,
                base + (nx + 1) + 2 * (ny - 1) + i,
                base + (nx + 1) + 2 * (j - 1) + np.where(i == 0, 0, 1),
            ),
        )
        return np.where((k == 0) | (k == nz), cap, ring)

    def faces(self):
        """``(face index, quads)`` of the six faces, outward, in order."""
        nx, ny, nz = self.nx, self.ny, self.nz
        out = []

        def grid(na, nb):
            b, a = np.divmod(np.arange(na * nb, dtype=np.int64), na)
            return a, b

        j, k = grid(ny, nz)
        x0, x1 = np.zeros_like(j), np.full_like(j, nx)
        f = self.id
        out.append(
            np.stack(
                [f(x0, j, k), f(x0, j, k + 1), f(x0, j + 1, k + 1), f(x0, j + 1, k)], 1
            )
        )
        out.append(
            np.stack(
                [f(x1, j, k), f(x1, j + 1, k), f(x1, j + 1, k + 1), f(x1, j, k + 1)], 1
            )
        )
        i, k = grid(nx, nz)
        y0, y1 = np.zeros_like(i), np.full_like(i, ny)
        out.append(
            np.stack(
                [f(i, y0, k), f(i + 1, y0, k), f(i + 1, y0, k + 1), f(i, y0, k + 1)], 1
            )
        )
        out.append(
            np.stack(
                [f(i, y1, k), f(i, y1, k + 1), f(i + 1, y1, k + 1), f(i + 1, y1, k)], 1
            )
        )
        i, j = grid(nx, ny)
        z0, z1 = np.zeros_like(i), np.full_like(i, nz)
        out.append(
            np.stack(
                [f(i, j, z0), f(i, j + 1, z0), f(i + 1, j + 1, z0), f(i + 1, j, z0)], 1
            )
        )
        out.append(
            np.stack(
                [f(i, j, z1), f(i + 1, j, z1), f(i + 1, j + 1, z1), f(i, j + 1, z1)], 1
            )
        )
        return out  # fmt: skip


def _quad_rows(quads, tri):
    if not tri:
        return quads
    a, b, c, d = quads[:, 0], quads[:, 1], quads[:, 2], quads[:, 3]
    return np.stack([np.stack([a, b, c], 1), np.stack([a, c, d], 1)], 1).reshape(-1, 3)


def _cube_sphere_py(radius, n, center, cell_type, order, complete):
    lat = _CubeLattice(n, n, n)
    i, j, k = lat.ijk()
    soup = _Soup("project", radius, center)
    soup.add(
        np.stack(
            [(2 * v - n).astype(np.float64) / float(n) for v in (i, j, k)], axis=1
        ),
        np.zeros(i.size),
    )
    tri = cell_type == "triangle"
    rows = np.concatenate([_quad_rows(q, tri) for q in lat.faces()])
    soup.block(cell_type, rows)
    return _finish(soup, order, [], 2, complete, (), "cube_sphere")


def cube_sphere(
    radius=1.0,
    n=4,
    center=(0.0, 0.0, 0.0),
    cell_type="quad",
    order=1,
    max_cells=DEFAULT_MAX_CELLS,
    complete=False,
):
    """A sphere surface from a cube's six faces, each an ``n x n`` grid, pushed
    radially onto the sphere (``6 n**2`` quads, or twice as many triangles).

    Points are the cube's boundary lattice, x fastest, then y, then z; cells run
    face by face (``xmin``, ``xmax``, ``ymin``, ``ymax``, ``zmin``, ``zmax``).

    Parameters
    ----------
    radius :
        The radius, positive.
    n :
        Cells along each cube edge; a zero gives an empty mesh.
    center :
        The centre.
    cell_type :
        ``"quad"`` or ``"triangle"`` (each quad cut along its first diagonal).
    order :
        1, 2 for ``quad8`` / ``triangle6`` on the sphere, or 3 to 10.
    max_cells :
        Refuse by name above this many cells; non-positive disables the check.
    complete :
        With ``order=2``, ``quad9`` instead of ``quad8``.
    """
    shape = "cube_sphere"
    center = _vec(shape, "center", center, 3)
    _check_positive(shape, "radius", float(radius))
    _check_count(shape, "n", int(n))
    _check_order(shape, order)
    _cell_choice(shape, cell_type, ("quad", "triangle"))
    if n == 0:
        return _empty()
    _budget(
        shape, 6 * int(n) * int(n) * (2 if cell_type == "triangle" else 1), max_cells
    )
    return _call(
        "cube_sphere",
        lambda: _cube_sphere_py(
            float(radius), int(n), center, cell_type, int(order), bool(complete)
        ),
        float(radius),
        int(n),
        _floats(center),
        cell_type,
        int(order),
        int(max_cells),
        bool(complete),
    )


_BOX_FACES = ("xmin", "xmax", "ymin", "ymax", "zmin", "zmax")


def _box_surface_py(dims, size, origin, cell_type, order, complete):
    nx, ny, nz = (int(v) for v in dims)
    lat = _CubeLattice(nx, ny, nz)
    i, j, k = lat.ijk()
    soup = _Soup("identity")
    soup.add(
        np.stack(
            [
                origin[0] + (size[0] * i.astype(np.float64)) / float(nx),
                origin[1] + (size[1] * j.astype(np.float64)) / float(ny),
                origin[2] + (size[2] * k.astype(np.float64)) / float(nz),
            ],
            axis=1,
        ),
        np.zeros(i.size),
    )
    tri = cell_type == "triangle"
    rows, labels = [], []
    for f, quads in enumerate(lat.faces()):
        r = _quad_rows(quads, tri)
        rows.append(r)
        labels.append(np.full(r.shape[0], f, dtype=np.int64))
    soup.block(cell_type, np.concatenate(rows))
    extra = _cell_regions(np.concatenate(labels), _BOX_FACES, 2)
    return _finish(soup, order, [], 2, complete, extra, "box_surface")


def box_surface(
    dims=(1, 1, 1),
    size=(1.0, 1.0, 1.0),
    origin=(0.0, 0.0, 0.0),
    cell_type="quad",
    order=1,
    max_cells=DEFAULT_MAX_CELLS,
    complete=False,
):
    """The boundary of a box as a closed surface, oriented outward.

    Each face is a grid of the two counts in its plane and shares its edges with
    its neighbours; points and cells are ordered as in :func:`cube_sphere`.

    Parameters
    ----------
    dims :
        Cell counts ``(nx, ny, nz)``; a zero on any axis gives an empty mesh.
    size :
        Side lengths, positive on every axis with cells.
    origin :
        The lo corner.
    cell_type :
        ``"quad"`` or ``"triangle"``.
    order :
        1, 2 for ``quad8`` / ``triangle6``, or 3 to 10.
    max_cells :
        Refuse by name above this many cells; non-positive disables the check.
    complete :
        With ``order=2``, ``quad9`` instead of ``quad8``.

    Returns
    -------
    Mesh
        Cell regions ``xmin`` ... ``zmax``, one per face.
    """
    shape = "box_surface"
    dims = _ivec(shape, "dims", dims, 3)
    size = _vec(shape, "size", size, 3)
    origin = _vec(shape, "origin", origin, 3)
    for k in range(3):
        _check_count(shape, "cell counts", dims[k])
        if dims[k] > 0:
            _check_positive(shape, "size", size[k])
    _check_order(shape, order)
    _cell_choice(shape, cell_type, ("quad", "triangle"))
    if 0 in (dims[0], dims[1], dims[2]):
        return _empty()
    nx, ny, nz = (int(v) for v in dims)
    faces = 2 * (nx * ny + ny * nz + nx * nz)
    _budget(shape, faces * (2 if cell_type == "triangle" else 1), max_cells)
    return _call(
        "box_surface",
        lambda: _box_surface_py(
            dims, size, origin, cell_type, int(order), bool(complete)
        ),
        [int(v) for v in dims],
        _floats(size),
        _floats(origin),
        cell_type,
        int(order),
        int(max_cells),
        bool(complete),
    )


# --------------------------------------------------------------------------- #
# geo_sphere                                                                  #
# --------------------------------------------------------------------------- #


def _icosa_verts():
    t = (1.0 + np.sqrt(5.0)) / 2.0
    length = np.sqrt(1.0 + t * t)
    one, tau = 1.0 / length, t / length
    return np.array(
        [
            [-one, tau, 0], [one, tau, 0], [-one, -tau, 0], [one, -tau, 0],
            [0, -one, tau], [0, one, tau], [0, -one, -tau], [0, one, -tau],
            [tau, 0, -one], [tau, 0, one], [-tau, 0, -one], [-tau, 0, one],
        ],
        dtype=np.float64,
    )  # fmt: skip


def _geo_sphere_py(radius, n, center, order, complete):
    verts = _icosa_verts()
    faces = _ICOSA_FACES
    soup = _Soup("project", radius, center)
    soup.add(verts, np.zeros(12))
    edge_index = {}
    edges = []
    for f in faces:
        for e in range(3):
            a, b = int(f[e]), int(f[(e + 1) % 3])
            key = (min(a, b), max(a, b))
            if key not in edge_index:
                edge_index[key] = len(edges)
                edges.append(key)
    s = np.arange(1, n, dtype=np.int64)
    if n > 1:
        elo = np.array([e[0] for e in edges])
        ehi = np.array([e[1] for e in edges])
        wa = (n - s).astype(np.float64)
        wb = s.astype(np.float64)
        pts = (
            verts[elo][:, None, :] * wa[None, :, None]
            + verts[ehi][:, None, :] * wb[None, :, None]
        ) / float(n)
        soup.add(pts.reshape(-1, 3), np.zeros(pts.shape[0] * pts.shape[1]))
    jj, ii = [], []
    for j in range(1, n - 1):
        for i in range(1, n - j):
            jj.append(j)
            ii.append(i)
    if jj:
        jv = np.array(jj, dtype=np.int64)
        iv = np.array(ii, dtype=np.int64)
        w0 = (n - iv - jv).astype(np.float64)
        w1 = iv.astype(np.float64)
        w2 = jv.astype(np.float64)
        pts = (
            verts[faces[:, 0]][:, None, :] * w0[None, :, None]
            + verts[faces[:, 1]][:, None, :] * w1[None, :, None]
            + verts[faces[:, 2]][:, None, :] * w2[None, :, None]
        ) / float(n)
        soup.add(pts.reshape(-1, 3), np.zeros(pts.shape[0] * pts.shape[1]))
    edge_base = 12
    face_base = edge_base + 30 * (n - 1)
    per_face = ((n - 1) * (n - 2)) // 2

    def on_edge(frm, to, steps):
        lo, hi = min(frm, to), max(frm, to)
        e = edge_index[(lo, hi)]
        st = steps if frm == lo else n - steps
        return edge_base + e * (n - 1) + (st - 1)

    def ident(fi, i, j):
        f = faces[fi]
        if i == 0 and j == 0:
            return int(f[0])
        if i == n:
            return int(f[1])
        if j == n:
            return int(f[2])
        if j == 0:
            return on_edge(int(f[0]), int(f[1]), i)
        if i == 0:
            return on_edge(int(f[0]), int(f[2]), j)
        if i + j == n:
            return on_edge(int(f[1]), int(f[2]), j)
        return (
            face_base + fi * per_face + (j - 1) * (n - 1) - ((j - 1) * j) // 2 + (i - 1)
        )

    rows = []
    for fi in range(20):
        for j in range(n):
            for i in range(n - j):
                rows.append((ident(fi, i, j), ident(fi, i + 1, j), ident(fi, i, j + 1)))
                if i < n - 1 - j:
                    rows.append(
                        (
                            ident(fi, i + 1, j),
                            ident(fi, i + 1, j + 1),
                            ident(fi, i, j + 1),
                        )
                    )
    soup.block("triangle", np.array(rows, dtype=np.int64))
    return _finish(soup, order, [], 2, complete, (), "geo_sphere")


def geo_sphere(
    radius=1.0,
    n=4,
    center=(0.0, 0.0, 0.0),
    order=1,
    max_cells=DEFAULT_MAX_CELLS,
    complete=False,
):
    """A geodesic sphere of frequency ``n``: every icosahedron face cut into an
    ``n x n`` triangle lattice (``20 n**2`` triangles) pushed radially onto the
    sphere. Unlike :func:`icosa_sphere`, any frequency, not only powers of two.

    Points are the 12 icosahedron vertices, then ``n - 1`` per edge (edges in
    first-appearance order over the faces), then each face's interior.

    Parameters
    ----------
    radius :
        The radius, positive.
    n :
        The frequency; a zero gives an empty mesh.
    center :
        The centre.
    order :
        1, 2 for ``triangle6`` on the sphere, or 3 to 10.
    max_cells :
        Refuse by name above this many cells; non-positive disables the check.
    complete :
        Accepted for symmetry.
    """
    shape = "geo_sphere"
    center = _vec(shape, "center", center, 3)
    _check_positive(shape, "radius", float(radius))
    _check_count(shape, "n", int(n))
    _check_order(shape, order)
    if n == 0:
        return _empty()
    _budget(shape, 20 * int(n) * int(n), max_cells)
    return _call(
        "geo_sphere",
        lambda: _geo_sphere_py(
            float(radius), int(n), center, int(order), bool(complete)
        ),
        float(radius),
        int(n),
        _floats(center),
        int(order),
        int(max_cells),
        bool(complete),
    )


# --------------------------------------------------------------------------- #
# torus, moebius, tube                                                        #
# --------------------------------------------------------------------------- #


def _torus_py(major, minor, nu, nv, center, cell_type, order, complete):
    soup = _Soup("torus", major, center)
    soup.radius2 = minor
    i = np.repeat(np.arange(nu, dtype=np.int64), nv)
    j = np.tile(np.arange(nv, dtype=np.int64), nu)
    soup.add(
        np.stack(
            [
                (2.0 * _PI * i.astype(np.float64)) / float(nu),
                (2.0 * _PI * j.astype(np.float64)) / float(nv),
                np.zeros(i.size),
            ],
            axis=1,
        ),
        np.zeros(i.size),
    )
    i1, j1 = (i + 1) % nu, (j + 1) % nv
    quads = np.stack([i * nv + j, i1 * nv + j, i1 * nv + j1, i * nv + j1], axis=1)
    soup.block(cell_type, _quad_rows(quads, cell_type == "triangle"))
    return _finish(soup, order, [], 2, complete, (), "torus")


def torus(
    major_radius=1.0,
    minor_radius=0.25,
    nu=32,
    nv=16,
    center=(0.0, 0.0, 0.0),
    cell_type="quad",
    order=1,
    max_cells=DEFAULT_MAX_CELLS,
    complete=False,
):
    """A torus around the z axis, periodic both ways and welded.

    Genus one (Euler characteristic zero), with negative Gaussian curvature on
    the inner rim: the curvature probe a sphere cannot give. Point ``(i, j)``
    is ``i * nv + j``, ``i`` around the axis and ``j`` around the tube.

    Parameters
    ----------
    major_radius :
        The radius of the tube's centre circle, positive.
    minor_radius :
        The tube's radius, positive and below ``major_radius``.
    nu, nv :
        Cells around the axis and around the tube, at least 3 each; a zero
        gives an empty mesh.
    center :
        The centre.
    cell_type :
        ``"quad"`` or ``"triangle"``.
    order :
        1, 2 for ``quad8`` / ``triangle6`` on the torus, or 3 to 10.
    max_cells :
        Refuse by name above this many cells; non-positive disables the check.
    complete :
        With ``order=2``, ``quad9`` instead of ``quad8``.
    """
    shape = "torus"
    center = _vec(shape, "center", center, 3)
    _check_positive(shape, "major_radius", float(major_radius))
    _check_positive(shape, "minor_radius", float(minor_radius))
    if not float(minor_radius) < float(major_radius):
        _fail(
            shape,
            "minor_radius must be less than major_radius, got "
            f"{_fmt(minor_radius)} and {_fmt(major_radius)}",
        )
    _check_count(shape, "nu", int(nu))
    _check_count(shape, "nv", int(nv))
    _check_order(shape, order)
    _cell_choice(shape, cell_type, ("quad", "triangle"))
    if nu == 0 or nv == 0:
        return _empty()
    if nu < 3 or nv < 3:
        _fail(shape, f"nu and nv must be at least 3, got {int(nu)} and {int(nv)}")
    _budget(shape, int(nu) * int(nv) * (2 if cell_type == "triangle" else 1), max_cells)
    return _call(
        "torus",
        lambda: _torus_py(
            float(major_radius),
            float(minor_radius),
            int(nu),
            int(nv),
            center,
            cell_type,
            int(order),
            bool(complete),
        ),
        float(major_radius),
        float(minor_radius),
        int(nu),
        int(nv),
        _floats(center),
        cell_type,
        int(order),
        int(max_cells),
        bool(complete),
    )


def _moebius_py(radius, width, nu, nv, twists, center, cell_type, order, complete):
    soup = _Soup("moebius", radius, center)
    soup.half = width / 2.0
    soup.twists = twists
    p = nv + 1
    i = np.repeat(np.arange(nu, dtype=np.int64), p)
    j = np.tile(np.arange(p, dtype=np.int64), nu)
    soup.add(
        np.stack(
            [
                (2.0 * _PI * i.astype(np.float64)) / float(nu),
                (2 * j - nv).astype(np.float64) / float(nv),
                np.zeros(i.size),
            ],
            axis=1,
        ),
        ((j == 0) | (j == nv)).astype(np.uint16),
    )
    ci = np.repeat(np.arange(nu, dtype=np.int64), nv)
    cj = np.tile(np.arange(nv, dtype=np.int64), nu)
    a = ci * p + cj
    d = a + 1
    seam = ci + 1 == nu
    flip = twists % 2 != 0
    b = np.where(seam, (nv - cj) if flip else cj, (ci + 1) * p + cj)
    c = np.where(seam, (nv - cj - 1) if flip else cj + 1, (ci + 1) * p + cj + 1)
    quads = np.stack([a, b, c, d], axis=1)
    soup.block(cell_type, _quad_rows(quads, cell_type == "triangle"))
    return _finish(soup, order, [("boundary", 1)], 1, complete, (), "moebius")


def moebius(
    radius=1.0,
    width=0.5,
    nu=32,
    nv=4,
    num_twists=1,
    center=(0.0, 0.0, 0.0),
    cell_type="quad",
    order=1,
    max_cells=DEFAULT_MAX_CELLS,
    complete=False,
):
    """A Moebius band around the z axis with ``num_twists`` half-twists.

    An odd count has no consistent orientation (one boundary loop): the probe
    for :func:`compute_normals`' ``inconsistent_pairs`` and
    ``repair(fix_orientation=True)``. An even count is an orientable band with
    two loops. Point ``(i, j)`` is ``i * (nv + 1) + j``, ``i`` along the band
    and ``j`` across it.

    Parameters
    ----------
    radius :
        The radius of the centre line, positive.
    width :
        The band's width, positive and below twice ``radius``.
    nu :
        Cells along the band, at least 3.
    nv :
        Cells across it.
    num_twists :
        Half-twists, not negative (1 is the classic band).
    center :
        The centre.
    cell_type :
        ``"quad"`` or ``"triangle"``.
    order :
        1, 2 for ``quad8`` / ``triangle6``, or 3 to 10.
    max_cells :
        Refuse by name above this many cells; non-positive disables the check.
    complete :
        With ``order=2``, ``quad9`` instead of ``quad8``.

    Returns
    -------
    Mesh
        One side region, ``boundary``.
    """
    shape = "moebius"
    center = _vec(shape, "center", center, 3)
    _check_positive(shape, "radius", float(radius))
    _check_positive(shape, "width", float(width))
    if not float(width) < 2.0 * float(radius):
        _fail(
            shape,
            "width must be less than twice the radius, got "
            f"{_fmt(width)} and {_fmt(radius)}",
        )
    _check_count(shape, "nu", int(nu))
    _check_count(shape, "nv", int(nv))
    _check_count(shape, "num_twists", int(num_twists))
    _check_order(shape, order)
    _cell_choice(shape, cell_type, ("quad", "triangle"))
    if nu == 0 or nv == 0:
        return _empty()
    if nu < 3:
        _fail(shape, f"nu must be at least 3, got {int(nu)}")
    _budget(shape, int(nu) * int(nv) * (2 if cell_type == "triangle" else 1), max_cells)
    return _call(
        "moebius",
        lambda: _moebius_py(
            float(radius),
            float(width),
            int(nu),
            int(nv),
            int(num_twists),
            center,
            cell_type,
            int(order),
            bool(complete),
        ),
        float(radius),
        float(width),
        int(nu),
        int(nv),
        int(num_twists),
        _floats(center),
        cell_type,
        int(order),
        int(max_cells),
        bool(complete),
    )


def _tube_py(radius, height, ntheta, nz, center, cell_type, order, complete):
    soup = _Soup("polar", radius, center)
    k = np.repeat(np.arange(nz + 1, dtype=np.int64), ntheta)
    i = np.tile(np.arange(ntheta, dtype=np.int64), nz + 1)
    soup.add(
        np.stack(
            [
                np.ones(k.size),
                (2.0 * _PI * i.astype(np.float64)) / float(ntheta),
                (height * (2 * k - nz).astype(np.float64)) / float(2 * nz),
            ],
            axis=1,
        ),
        ((k == 0) * 1 + (k == nz) * 2).astype(np.uint16),
    )
    ck = np.repeat(np.arange(nz, dtype=np.int64), ntheta)
    ci = np.tile(np.arange(ntheta, dtype=np.int64), nz)
    i1 = (ci + 1) % ntheta
    quads = np.stack(
        [
            ck * ntheta + ci,
            ck * ntheta + i1,
            (ck + 1) * ntheta + i1,
            (ck + 1) * ntheta + ci,
        ],
        axis=1,
    )
    soup.block(cell_type, _quad_rows(quads, cell_type == "triangle"))
    return _finish(soup, order, [("bottom", 1), ("top", 2)], 1, complete, (), "tube")


def tube(
    radius=1.0,
    height=1.0,
    ntheta=16,
    nz=4,
    center=(0.0, 0.0, 0.0),
    cell_type="quad",
    order=1,
    max_cells=DEFAULT_MAX_CELLS,
    complete=False,
):
    """An open cylinder surface around the z axis, centred on ``center`` and
    oriented outward. Point ``(k, i)`` is ``k * ntheta + i``.

    Parameters
    ----------
    radius, height :
        Positive.
    ntheta :
        Cells around, at least 3.
    nz :
        Cells along the axis.
    center :
        The centre.
    cell_type :
        ``"quad"`` or ``"triangle"``.
    order :
        1, 2 for ``quad8`` / ``triangle6`` on the cylinder, or 3 to 10.
    max_cells :
        Refuse by name above this many cells; non-positive disables the check.
    complete :
        With ``order=2``, ``quad9`` instead of ``quad8``.

    Returns
    -------
    Mesh
        Side regions ``bottom`` and ``top`` (its two rims).
    """
    shape = "tube"
    center = _vec(shape, "center", center, 3)
    _check_positive(shape, "radius", float(radius))
    _check_positive(shape, "height", float(height))
    _check_count(shape, "ntheta", int(ntheta))
    _check_count(shape, "nz", int(nz))
    _check_order(shape, order)
    _cell_choice(shape, cell_type, ("quad", "triangle"))
    if ntheta == 0 or nz == 0:
        return _empty()
    if ntheta < 3:
        _fail(shape, f"ntheta must be at least 3, got {int(ntheta)}")
    _budget(
        shape, int(ntheta) * int(nz) * (2 if cell_type == "triangle" else 1), max_cells
    )
    return _call(
        "tube",
        lambda: _tube_py(
            float(radius),
            float(height),
            int(ntheta),
            int(nz),
            center,
            cell_type,
            int(order),
            bool(complete),
        ),
        float(radius),
        float(height),
        int(ntheta),
        int(nz),
        _floats(center),
        cell_type,
        int(order),
        int(max_cells),
        bool(complete),
    )


# --------------------------------------------------------------------------- #
# cone, capsule: closed surfaces of revolution in rings between two poles     #
# --------------------------------------------------------------------------- #


def _ring_cells(soup, rings, ntheta, last, tri, band_label):
    """Cells between pole 0, ``rings`` rings and pole ``last``, laid out as
    :func:`uv_sphere`'s, and the label of every cell (block-major)."""

    def ring(r, i):
        return 1 + (r - 1) * ntheta + (i % ntheta)

    ii = np.arange(ntheta, dtype=np.int64)
    br = np.repeat(np.arange(1, rings, dtype=np.int64), ntheta)
    bi = np.tile(ii, rings - 1)
    a, b = ring(br, bi), ring(br, bi + 1)
    c, d = ring(br + 1, bi + 1), ring(br + 1, bi)
    body_labels = np.array([band_label(int(r)) for r in br], dtype=np.int64)
    labels = []
    if not tri:
        soup.block("quad", np.stack([a, d, c, b], axis=1))
        labels.append(body_labels)
    first = np.stack(
        [np.zeros(ntheta, dtype=np.int64), ring(1, ii), ring(1, ii + 1)], 1
    )
    second = np.stack(
        [np.full(ntheta, last, dtype=np.int64), ring(rings, ii + 1), ring(rings, ii)], 1
    )
    parts = [first]
    labels.append(np.full(ntheta, band_label(0), dtype=np.int64))
    if tri:
        body = np.stack([np.stack([a, d, c], 1), np.stack([a, c, b], 1)], axis=1)
        parts.append(body.reshape(-1, 3))
        labels.append(np.repeat(body_labels, 2))
    parts.append(second)
    labels.append(np.full(ntheta, band_label(rings), dtype=np.int64))
    soup.block("triangle", np.concatenate(parts))
    return np.concatenate(labels)


def _ring_count(rings, ntheta, tri):
    return ntheta * (2 * rings if tri else rings + 1)


def _cone_py(radius, height, ntheta, nz, center, cell_type, order, complete):
    soup = _Soup("polar", radius, center)
    soup.add([[0.0, 0.0, height]], [_DEGENERATE])
    rings = 2 * nz - 1
    r = np.repeat(np.arange(1, rings + 1, dtype=np.int64), ntheta)
    i = np.tile(np.arange(ntheta, dtype=np.int64), rings)
    m = np.where(r <= nz, r, 2 * nz - r)
    z = np.where(r <= nz, (height * (nz - r).astype(np.float64)) / float(nz), 0.0)
    soup.add(
        np.stack(
            [
                m.astype(np.float64) / float(nz),
                (2.0 * _PI * i.astype(np.float64)) / float(ntheta),
                z,
            ],
            axis=1,
        ),
        np.zeros(r.size),
    )
    last = soup.add([[0.0, 0.0, 0.0]], [_DEGENERATE])
    labels = _ring_cells(
        soup, rings, ntheta, last, cell_type == "triangle", lambda b: 0 if b < nz else 1
    )
    extra = _cell_regions(labels, ("side", "base"), 2)
    return _finish(soup, order, [], 2, complete, extra, "cone")


def cone(
    radius=1.0,
    height=1.0,
    ntheta=16,
    nz=4,
    center=(0.0, 0.0, 0.0),
    cell_type="quad",
    order=1,
    max_cells=DEFAULT_MAX_CELLS,
    complete=False,
):
    """A closed cone surface: the base disk at ``z = center[2]``, the apex
    ``height`` above it, the rim a sharp edge (a :func:`feature_edges` probe).

    Point 0 is the apex, then ``2 nz - 1`` rings of ``ntheta`` points down the
    slant to the rim and in across the base, then the base centre.

    Parameters
    ----------
    radius, height :
        Positive.
    ntheta :
        Cells around, at least 3.
    nz :
        Rings on the slant and on the base.
    center :
        The centre of the base.
    cell_type :
        ``"quad"`` (triangle fans at the apex and the base centre) or
        ``"triangle"``.
    order :
        1, 2 for ``quad8`` / ``triangle6``, or 3 to 10.
    max_cells :
        Refuse by name above this many cells; non-positive disables the check.
    complete :
        With ``order=2``, ``quad9`` instead of ``quad8``.

    Returns
    -------
    Mesh
        Cell regions ``side`` and ``base``.
    """
    shape = "cone"
    center = _vec(shape, "center", center, 3)
    _check_positive(shape, "radius", float(radius))
    _check_positive(shape, "height", float(height))
    _check_count(shape, "ntheta", int(ntheta))
    _check_count(shape, "nz", int(nz))
    _check_order(shape, order)
    _cell_choice(shape, cell_type, ("quad", "triangle"))
    if ntheta == 0 or nz == 0:
        return _empty()
    if ntheta < 3:
        _fail(shape, f"ntheta must be at least 3, got {int(ntheta)}")
    _budget(
        shape,
        _ring_count(2 * int(nz) - 1, int(ntheta), cell_type == "triangle"),
        max_cells,
    )
    return _call(
        "cone",
        lambda: _cone_py(
            float(radius),
            float(height),
            int(ntheta),
            int(nz),
            center,
            cell_type,
            int(order),
            bool(complete),
        ),
        float(radius),
        float(height),
        int(ntheta),
        int(nz),
        _floats(center),
        cell_type,
        int(order),
        int(max_cells),
        bool(complete),
    )


def _capsule_py(radius, height, ntheta, nz, ncap, center, cell_type, order, complete):
    soup = _Soup("capsule", radius, center)
    soup.half = height / 2.0
    side = (2.0 * soup.half) / radius
    soup.add([[0.0, 0.0, 0.0]], [_DEGENERATE])
    rings = 2 * ncap + nz - 1
    t = []
    for r in range(1, rings + 1):
        if r <= ncap:
            t.append((_PI * float(r)) / float(2 * ncap))
        elif r < ncap + nz:
            t.append(_PI / 2.0 + (side * float(r - ncap)) / float(nz))
        else:
            t.append(_PI / 2.0 + side + (_PI * float(r - ncap - nz)) / float(2 * ncap))
    tv = np.repeat(np.array(t), ntheta)
    i = np.tile(np.arange(ntheta, dtype=np.int64), rings)
    soup.add(
        np.stack(
            [tv, (2.0 * _PI * i.astype(np.float64)) / float(ntheta), np.zeros(tv.size)],
            axis=1,
        ),
        np.zeros(tv.size),
    )
    last = soup.add([[_PI + side, 0.0, 0.0]], [_DEGENERATE])
    _ring_cells(soup, rings, ntheta, last, cell_type == "triangle", lambda b: 0)
    return _finish(soup, order, [], 2, complete, (), "capsule")


def capsule(
    radius=0.5,
    height=1.0,
    ntheta=16,
    nz=4,
    ncap=4,
    center=(0.0, 0.0, 0.0),
    cell_type="quad",
    order=1,
    max_cells=DEFAULT_MAX_CELLS,
    complete=False,
):
    """A capsule: a cylinder of length ``height`` closed by two hemispheres,
    centred on ``center`` along z, smooth (G1) where they meet.

    Point 0 is the north pole, then ``2 ncap + nz - 1`` rings of ``ntheta``
    points southward, then the south pole.

    Parameters
    ----------
    radius, height :
        Positive.
    ntheta :
        Cells around, at least 3.
    nz :
        Rings along the cylinder.
    ncap :
        Rings on each cap.
    center :
        The centre.
    cell_type :
        ``"quad"`` (triangle fans at the poles) or ``"triangle"``.
    order :
        1, 2 for ``quad8`` / ``triangle6`` on the exact caps, or 3 to 10.
    max_cells :
        Refuse by name above this many cells; non-positive disables the check.
    complete :
        With ``order=2``, ``quad9`` instead of ``quad8``.
    """
    shape = "capsule"
    center = _vec(shape, "center", center, 3)
    _check_positive(shape, "radius", float(radius))
    _check_positive(shape, "height", float(height))
    _check_count(shape, "ntheta", int(ntheta))
    _check_count(shape, "nz", int(nz))
    _check_count(shape, "ncap", int(ncap))
    _check_order(shape, order)
    _cell_choice(shape, cell_type, ("quad", "triangle"))
    if ntheta == 0 or nz == 0 or ncap == 0:
        return _empty()
    if ntheta < 3:
        _fail(shape, f"ntheta must be at least 3, got {int(ntheta)}")
    rings = 2 * int(ncap) + int(nz) - 1
    _budget(shape, _ring_count(rings, int(ntheta), cell_type == "triangle"), max_cells)
    return _call(
        "capsule",
        lambda: _capsule_py(
            float(radius),
            float(height),
            int(ntheta),
            int(nz),
            int(ncap),
            center,
            cell_type,
            int(order),
            bool(complete),
        ),
        float(radius),
        float(height),
        int(ntheta),
        int(nz),
        int(ncap),
        _floats(center),
        cell_type,
        int(order),
        int(max_cells),
        bool(complete),
    )


# --------------------------------------------------------------------------- #
# platonic                                                                    #
# --------------------------------------------------------------------------- #

_PLATONIC = ("tetrahedron", "cube", "octahedron", "dodecahedron", "icosahedron")
_CUBE_SURFACE = [
    [0, 1, 5, 4],
    [0, 2, 3, 1],
    [0, 4, 6, 2],
    [1, 3, 7, 5],
    [2, 6, 7, 3],
    [4, 5, 7, 6],
]
_DODECA_FACES = [
    [0, 8, 1, 14, 12], [0, 12, 4, 17, 16], [0, 16, 2, 10, 8], [1, 8, 10, 3, 18],
    [1, 18, 19, 5, 14], [2, 13, 15, 3, 10], [2, 16, 17, 6, 13], [3, 15, 7, 19, 18],
    [4, 9, 11, 6, 17], [4, 12, 14, 5, 9], [5, 19, 7, 11, 9], [6, 11, 7, 15, 13],
]  # fmt: skip


def _cube_corners(a):
    return [[x * a, y * a, z * a] for z in (-1, 1) for y in (-1, 1) for x in (-1, 1)]


def _platonic_py(name, radius, center):
    a = 1.0 / np.sqrt(3.0)
    if name == "tetrahedron":
        verts, faces, kind = _tetra_verts(), _TETRA_FACES, "triangle"
    elif name == "cube":
        verts, faces, kind = np.array(_cube_corners(a)), _CUBE_SURFACE, "quad"
    elif name == "octahedron":
        verts, faces, kind = _OCTA_VERTS, _OCTA_FACES, "triangle"
    elif name == "dodecahedron":
        phi = (1.0 + np.sqrt(5.0)) / 2.0
        ip = 1.0 / phi
        v = _cube_corners(a)
        for s1 in (-1, 1):
            for s2 in (-1, 1):
                v.append([0.0, (s1 * ip) * a, (s2 * phi) * a])
        for s1 in (-1, 1):
            for s2 in (-1, 1):
                v.append([(s1 * ip) * a, (s2 * phi) * a, 0.0])
        for s1 in (-1, 1):
            for s2 in (-1, 1):
                v.append([(s1 * phi) * a, 0.0, (s2 * ip) * a])
        verts, faces, kind = np.array(v), _DODECA_FACES, "polygon"
    else:
        verts, faces, kind = _icosa_verts(), _ICOSA_FACES, "triangle"
    points = np.asarray(center, dtype=np.float64) + radius * np.asarray(verts)
    return Mesh(points, [(kind, np.array(faces, dtype=np.int64))])


def platonic(
    name="icosahedron", radius=1.0, center=(0.0, 0.0, 0.0), max_cells=DEFAULT_MAX_CELLS
):
    """A platonic solid's surface, one cell per face, inscribed in a sphere.

    Parameters
    ----------
    name :
        ``"tetrahedron"``, ``"octahedron"`` or ``"icosahedron"`` (triangles),
        ``"cube"`` (quads) or ``"dodecahedron"`` (a ``polygon`` block of
        pentagons, ragged in the core like :func:`honeycomb`'s).
    radius :
        The circumradius, positive.
    center :
        The centre.
    max_cells :
        Refuse by name above this many cells; non-positive disables the check.

    Returns
    -------
    Mesh
        Faces oriented outward; linear only.
    """
    shape = "platonic"
    center = _vec(shape, "center", center, 3)
    _check_positive(shape, "radius", float(radius))
    _choice(shape, "name", name, _PLATONIC)
    faces = {"tetrahedron": 4, "cube": 6, "octahedron": 8, "dodecahedron": 12}
    _budget(shape, faces.get(name, 20), max_cells)
    return _call(
        "platonic",
        lambda: _platonic_py(name, float(radius), center),
        name,
        float(radius),
        _floats(center),
        int(max_cells),
    )


# --------------------------------------------------------------------------- #
# ball, spherical_shell, pipe, fichera, solid_torus                           #
# --------------------------------------------------------------------------- #


def _hex_rows(h, tetra):
    return h[:, _KUHN].reshape(-1, 4) if tetra else h


def _ball_py(radius, n, center, cell_type, order, complete):
    m = 2 * n
    px = m + 1
    i = np.tile(np.arange(px, dtype=np.int64), px * px)
    j = np.tile(np.repeat(np.arange(px, dtype=np.int64), px), px)
    k = np.repeat(np.arange(px, dtype=np.int64), px * px)
    rim = (i == 0) | (i == m) | (j == 0) | (j == m) | (k == 0) | (k == m)
    soup = _Soup("ball", radius, center)
    soup.add(
        np.stack([(v - n).astype(np.float64) / float(n) for v in (i, j, k)], axis=1),
        rim.astype(np.uint16),
    )
    ck, rest = np.divmod(np.arange(m * m * m, dtype=np.int64), m * m)
    cj, ci = np.divmod(rest, m)
    b = (ck * px + cj) * px + ci
    t = b + px * px
    h = np.stack([b, b + 1, b + px + 1, b + px, t, t + 1, t + px + 1, t + px], axis=1)
    soup.block(cell_type, _hex_rows(h, cell_type == "tetra"))
    return _finish(soup, order, [("surface", 1)], 2, complete, (), "ball")


def ball(
    radius=1.0,
    n=4,
    center=(0.0, 0.0, 0.0),
    cell_type="hexahedron",
    order=1,
    max_cells=DEFAULT_MAX_CELLS,
    complete=False,
):
    """A solid ball by the seven-block butterfly: an inner cube of half side
    ``radius / 2`` and six blocks blending its faces into the sphere, so no cell
    degenerates at the centre.

    The lattice is ``(2 n + 1)**3`` points, x fastest; ``n`` cells run from the
    centre to the surface.

    Parameters
    ----------
    radius :
        The radius, positive.
    n :
        Cells from the centre to the surface; a zero gives an empty mesh.
    center :
        The centre.
    cell_type :
        ``"hexahedron"`` or ``"tetra"`` (six Kuhn tetrahedra per cell).
    order :
        1, 2 for ``hexahedron20`` / ``tetra10`` on the sphere, or 3 to 10.
    max_cells :
        Refuse by name above this many cells; non-positive disables the check.
    complete :
        With ``order=2``, ``hexahedron27`` instead of ``hexahedron20``.

    Returns
    -------
    Mesh
        One side region, ``surface``.
    """
    shape = "ball"
    center = _vec(shape, "center", center, 3)
    _check_positive(shape, "radius", float(radius))
    _check_count(shape, "n", int(n))
    _check_order(shape, order)
    _cell_choice(shape, cell_type, ("hexahedron", "tetra"))
    if n == 0:
        return _empty()
    m = 2 * int(n)
    _budget(shape, m * m * m * (6 if cell_type == "tetra" else 1), max_cells)
    return _call(
        "ball",
        lambda: _ball_py(
            float(radius), int(n), center, cell_type, int(order), bool(complete)
        ),
        float(radius),
        int(n),
        _floats(center),
        cell_type,
        int(order),
        int(max_cells),
        bool(complete),
    )


def _spherical_shell_py(inner, outer, n, nr, ratio, center, order, complete):
    lat = _CubeLattice(n, n, n)
    i, j, k = lat.ijk()
    d = np.stack([(2 * v - n).astype(np.float64) / float(n) for v in (i, j, k)], axis=1)
    soup = _Soup("shell", inner, center)
    soup.radius2 = outer
    soup.ratio = ratio
    soup.layers = float(nr)
    for layer in range(nr + 1):
        t = float(nr + layer) / float(nr)
        mask = (1 if layer == 0 else 0) | (2 if layer == nr else 0)
        soup.add(d * t, np.full(d.shape[0], mask))
    count = lat.count()
    quads = np.concatenate(lat.faces())
    rows = [
        np.concatenate([quads + layer * count, quads + (layer + 1) * count], axis=1)
        for layer in range(nr)
    ]
    soup.block("hexahedron", np.concatenate(rows))
    return _finish(
        soup, order, [("inner", 1), ("outer", 2)], 2, complete, (), "spherical_shell"
    )


def spherical_shell(
    inner_radius=0.5,
    outer_radius=1.0,
    n=4,
    nr=2,
    ratio=1.0,
    center=(0.0, 0.0, 0.0),
    order=1,
    max_cells=DEFAULT_MAX_CELLS,
    complete=False,
):
    """A thick spherical shell of hexahedra between two radii: the
    :func:`cube_sphere` lattice on each of ``nr + 1`` layers.

    Parameters
    ----------
    inner_radius, outer_radius :
        Positive, the outer larger.
    n :
        Cells along each cube edge.
    nr :
        Layers through the wall.
    ratio :
        Each layer's thickness over the one inside it, positive: 1 is uniform,
        above 1 grows outward, below 1 refines toward the outer wall.
    center :
        The centre.
    order :
        1, 2 for ``hexahedron20`` on the spheres, or 3 to 10.
    max_cells :
        Refuse by name above this many cells; non-positive disables the check.
    complete :
        With ``order=2``, ``hexahedron27`` instead of ``hexahedron20``.

    Returns
    -------
    Mesh
        Side regions ``inner`` and ``outer``.
    """
    shape = "spherical_shell"
    center = _vec(shape, "center", center, 3)
    _check_count(shape, "n", int(n))
    _check_count(shape, "nr", int(nr))
    _check_order(shape, order)
    if n == 0 or nr == 0:
        return _empty()
    _check_positive(shape, "inner_radius", float(inner_radius))
    if not float(outer_radius) > float(inner_radius):
        _fail(
            shape,
            "outer_radius must exceed inner_radius, got "
            f"{_fmt(outer_radius)} and {_fmt(inner_radius)}",
        )
    _check_positive(shape, "ratio", float(ratio))
    _budget(shape, 6 * int(n) * int(n) * int(nr), max_cells)
    return _call(
        "spherical_shell",
        lambda: _spherical_shell_py(
            float(inner_radius),
            float(outer_radius),
            int(n),
            int(nr),
            float(ratio),
            center,
            int(order),
            bool(complete),
        ),
        float(inner_radius),
        float(outer_radius),
        int(n),
        int(nr),
        float(ratio),
        _floats(center),
        int(order),
        int(max_cells),
        bool(complete),
    )


def _pipe_py(inner, outer, height, ntheta, nr, nz, ratio, center, order, complete):
    soup = _Soup("graded", inner, center)
    soup.radius2 = outer
    soup.ratio = ratio
    soup.layers = float(nr)
    k = np.repeat(np.arange(nz + 1, dtype=np.int64), (nr + 1) * ntheta)
    m = np.tile(np.repeat(np.arange(nr + 1, dtype=np.int64), ntheta), nz + 1)
    i = np.tile(np.arange(ntheta, dtype=np.int64), (nr + 1) * (nz + 1))
    mask = (m == 0) * 1 + (m == nr) * 2 + (k == 0) * 4 + (k == nz) * 8
    soup.add(
        np.stack(
            [
                m.astype(np.float64) / float(nr),
                (2.0 * _PI * i.astype(np.float64)) / float(ntheta),
                (height * (2 * k - nz).astype(np.float64)) / float(2 * nz),
            ],
            axis=1,
        ),
        mask.astype(np.uint16),
    )

    def ident(kk, mm, ii):
        return (kk * (nr + 1) + mm) * ntheta + (ii % ntheta)

    ck, rest = np.divmod(np.arange(nz * nr * ntheta, dtype=np.int64), nr * ntheta)
    cm, ci = np.divmod(rest, ntheta)
    h = np.stack(
        [
            ident(ck, cm, ci),
            ident(ck, cm + 1, ci),
            ident(ck, cm + 1, ci + 1),
            ident(ck, cm, ci + 1),
            ident(ck + 1, cm, ci),
            ident(ck + 1, cm + 1, ci),
            ident(ck + 1, cm + 1, ci + 1),
            ident(ck + 1, cm, ci + 1),
        ],
        axis=1,
    )
    soup.block("hexahedron", h)
    specs = [("inner", 1), ("outer", 2), ("bottom", 4), ("top", 8)]
    return _finish(soup, order, specs, 2, complete, (), "pipe")


def pipe(
    inner_radius=0.5,
    outer_radius=1.0,
    height=1.0,
    ntheta=16,
    nr=2,
    nz=4,
    ratio=1.0,
    center=(0.0, 0.0, 0.0),
    order=1,
    max_cells=DEFAULT_MAX_CELLS,
    complete=False,
):
    """A thick-walled pipe of hexahedra around the z axis, centred on ``center``.

    Point ``(k, m, i)`` is ``(k * (nr + 1) + m) * ntheta + i``: ``k`` along the
    axis, ``m`` through the wall, ``i`` around.

    Parameters
    ----------
    inner_radius, outer_radius :
        Positive, the outer larger.
    height :
        Positive.
    ntheta :
        Cells around, at least 3.
    nr :
        Layers through the wall.
    nz :
        Cells along the axis.
    ratio :
        Each layer's radial thickness over the one inside it, positive.
    center :
        The centre.
    order :
        1, 2 for ``hexahedron20`` on the walls, or 3 to 10.
    max_cells :
        Refuse by name above this many cells; non-positive disables the check.
    complete :
        With ``order=2``, ``hexahedron27`` instead of ``hexahedron20``.

    Returns
    -------
    Mesh
        Side regions ``inner``, ``outer``, ``bottom`` and ``top``.
    """
    shape = "pipe"
    center = _vec(shape, "center", center, 3)
    _check_count(shape, "ntheta", int(ntheta))
    _check_count(shape, "nr", int(nr))
    _check_count(shape, "nz", int(nz))
    _check_order(shape, order)
    if ntheta == 0 or nr == 0 or nz == 0:
        return _empty()
    _check_positive(shape, "inner_radius", float(inner_radius))
    if not float(outer_radius) > float(inner_radius):
        _fail(
            shape,
            "outer_radius must exceed inner_radius, got "
            f"{_fmt(outer_radius)} and {_fmt(inner_radius)}",
        )
    _check_positive(shape, "height", float(height))
    _check_positive(shape, "ratio", float(ratio))
    if ntheta < 3:
        _fail(shape, f"ntheta must be at least 3, got {int(ntheta)}")
    _budget(shape, int(ntheta) * int(nr) * int(nz), max_cells)
    return _call(
        "pipe",
        lambda: _pipe_py(
            float(inner_radius),
            float(outer_radius),
            float(height),
            int(ntheta),
            int(nr),
            int(nz),
            float(ratio),
            center,
            int(order),
            bool(complete),
        ),
        float(inner_radius),
        float(outer_radius),
        float(height),
        int(ntheta),
        int(nr),
        int(nz),
        float(ratio),
        _floats(center),
        int(order),
        int(max_cells),
        bool(complete),
    )


def _fichera_py(n, size, center, cell_type, order, complete):
    w = 2 * n + 1
    top = 2 * n
    k, j, i = np.meshgrid(np.arange(w), np.arange(w), np.arange(w), indexing="ij")
    i, j, k = (v.reshape(-1).astype(np.int64) for v in (i, j, k))
    keep = ~((i > n) & (j > n) & (k > n))
    ids = np.full(w * w * w, -1, dtype=np.int64)
    ids[keep] = np.arange(int(keep.sum()), dtype=np.int64)
    i, j, k = i[keep], j[keep], k[keep]
    mask = (
        (i == 0) * 1
        + (i == top) * 2
        + (j == 0) * 4
        + (j == top) * 8
        + (k == 0) * 16
        + (k == top) * 32
        + ((i == n) & (j >= n) & (k >= n)) * 64
        + ((j == n) & (i >= n) & (k >= n)) * 128
        + ((k == n) & (i >= n) & (j >= n)) * 256
    )
    soup = _Soup("identity")
    soup.add(
        np.stack(
            [
                center[c] + (size * (v - n).astype(np.float64)) / float(2 * n)
                for c, v in enumerate((i, j, k))
            ],
            axis=1,
        ),
        mask.astype(np.uint16),
    )
    ck, rest = np.divmod(np.arange(top**3, dtype=np.int64), top * top)
    cj, ci = np.divmod(rest, top)
    sel = ~((ci >= n) & (cj >= n) & (ck >= n))
    ci, cj, ck = ci[sel], cj[sel], ck[sel]

    def at(a, b, c):
        return ids[(c * w + b) * w + a]

    h = np.stack(
        [
            at(ci, cj, ck),
            at(ci + 1, cj, ck),
            at(ci + 1, cj + 1, ck),
            at(ci, cj + 1, ck),
            at(ci, cj, ck + 1),
            at(ci + 1, cj, ck + 1),
            at(ci + 1, cj + 1, ck + 1),
            at(ci, cj + 1, ck + 1),
        ],
        axis=1,
    )
    soup.block(cell_type, _hex_rows(h, cell_type == "tetra"))
    names = [
        "xmin",
        "xmax",
        "ymin",
        "ymax",
        "zmin",
        "zmax",
        "notch_x",
        "notch_y",
        "notch_z",
    ]
    specs = [(name, 1 << b) for b, name in enumerate(names)]
    return _finish(soup, order, specs, 2, complete, (), "fichera")


def fichera(
    n=2,
    size=2.0,
    center=(0.0, 0.0, 0.0),
    cell_type="hexahedron",
    order=1,
    max_cells=DEFAULT_MAX_CELLS,
    complete=False,
):
    """The Fichera corner: the cube ``[-size/2, size/2]**3`` around ``center``
    minus its positive octant, the 3-D L-shape with a re-entrant vertex (and
    three re-entrant edges) at the centre.

    The lattice is ``(2 n + 1)**3`` points, x fastest, the points strictly
    inside the removed octant skipped; ``7 n**3`` hexahedra.

    Parameters
    ----------
    n :
        Cells along half a side; a zero gives an empty mesh.
    size :
        The side of the bounding cube, positive.
    center :
        Its centre (the re-entrant vertex).
    cell_type :
        ``"hexahedron"`` or ``"tetra"``.
    order :
        1, 2 for ``hexahedron20`` / ``tetra10``, or 3 to 10.
    max_cells :
        Refuse by name above this many cells; non-positive disables the check.
    complete :
        With ``order=2``, ``hexahedron27`` instead of ``hexahedron20``.

    Returns
    -------
    Mesh
        Side regions ``xmin`` ... ``zmax`` and ``notch_x``, ``notch_y``,
        ``notch_z`` (the notch walls on ``x = 0``, ``y = 0``, ``z = 0``).
    """
    shape = "fichera"
    center = _vec(shape, "center", center, 3)
    _check_count(shape, "n", int(n))
    _check_order(shape, order)
    _cell_choice(shape, cell_type, ("hexahedron", "tetra"))
    if n == 0:
        return _empty()
    _check_positive(shape, "size", float(size))
    _budget(shape, 7 * int(n) ** 3 * (6 if cell_type == "tetra" else 1), max_cells)
    return _call(
        "fichera",
        lambda: _fichera_py(
            int(n), float(size), center, cell_type, int(order), bool(complete)
        ),
        int(n),
        float(size),
        _floats(center),
        cell_type,
        int(order),
        int(max_cells),
        bool(complete),
    )


def _solid_torus_py(major, minor, n, nu, center, order, complete):
    soup = _Soup("solid_torus", major, center)
    soup.radius2 = minor
    p = n + 1
    layer = p * p
    i = np.repeat(np.arange(nu, dtype=np.int64), layer)
    jj = np.tile(np.repeat(np.arange(p, dtype=np.int64), p), nu)
    ii = np.tile(np.arange(p, dtype=np.int64), p * nu)
    rim = (ii == 0) | (ii == n) | (jj == 0) | (jj == n)
    soup.add(
        np.stack(
            [
                (2 * ii - n).astype(np.float64) / float(n),
                (2 * jj - n).astype(np.float64) / float(n),
                (2.0 * _PI * i.astype(np.float64)) / float(nu),
            ],
            axis=1,
        ),
        (rim * 4).astype(np.uint16),
    )
    ci = np.repeat(np.arange(nu, dtype=np.int64), n * n)
    cj, cii = np.divmod(np.tile(np.arange(n * n, dtype=np.int64), nu), n)
    a = ci * layer + cj * p + cii
    b = ((ci + 1) % nu) * layer + cj * p + cii
    h = np.stack([a, a + p, a + p + 1, a + 1, b, b + p, b + p + 1, b + 1], axis=1)
    soup.block("hexahedron", h)
    return _finish(soup, order, [("surface", 4)], 2, complete, (), "solid_torus")


def solid_torus(
    major_radius=1.0,
    minor_radius=0.25,
    n=4,
    nu=32,
    center=(0.0, 0.0, 0.0),
    order=1,
    max_cells=DEFAULT_MAX_CELLS,
    complete=False,
):
    """A solid torus of hexahedra: :func:`disk`'s squircle lattice (``n`` cells
    across) as the tube's cross-section, swept around the z axis in ``nu``
    steps and welded. Point ``(i, jj, ii)`` is ``i * (n + 1)**2 + jj * (n + 1)
    + ii``, ``i`` around the axis.

    Parameters
    ----------
    major_radius :
        The radius of the tube's centre circle, positive.
    minor_radius :
        The tube's radius, positive and below ``major_radius``.
    n :
        Cells across the cross-section.
    nu :
        Cells around the axis, at least 3.
    center :
        The centre.
    order :
        1, 2 for ``hexahedron20`` on the torus, or 3 to 10.
    max_cells :
        Refuse by name above this many cells; non-positive disables the check.
    complete :
        With ``order=2``, ``hexahedron27`` instead of ``hexahedron20``.

    Returns
    -------
    Mesh
        One side region, ``surface``.
    """
    shape = "solid_torus"
    center = _vec(shape, "center", center, 3)
    _check_positive(shape, "major_radius", float(major_radius))
    _check_positive(shape, "minor_radius", float(minor_radius))
    if not float(minor_radius) < float(major_radius):
        _fail(
            shape,
            "minor_radius must be less than major_radius, got "
            f"{_fmt(minor_radius)} and {_fmt(major_radius)}",
        )
    _check_count(shape, "n", int(n))
    _check_count(shape, "nu", int(nu))
    _check_order(shape, order)
    if n == 0 or nu == 0:
        return _empty()
    if nu < 3:
        _fail(shape, f"nu must be at least 3, got {int(nu)}")
    _budget(shape, int(n) * int(n) * int(nu), max_cells)
    return _call(
        "solid_torus",
        lambda: _solid_torus_py(
            float(major_radius),
            float(minor_radius),
            int(n),
            int(nu),
            center,
            int(order),
            bool(complete),
        ),
        float(major_radius),
        float(minor_radius),
        int(n),
        int(nu),
        _floats(center),
        int(order),
        int(max_cells),
        bool(complete),
    )

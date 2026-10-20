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
]

#: Default cell budget, as for :func:`meshioplusplus.grid`.
DEFAULT_MAX_CELLS = 20000000

_PI = float(np.pi)
_DEGENERATE = 0x8000

_ELEVATED = {
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
    if order not in (1, 2):
        _fail(shape, f"order must be 1 or 2, got {order}")


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
        self.kind = kind  # identity | squircle | polar | sphere | unit | plate
        self.radius = float(radius)
        self.center = np.asarray(center, dtype=np.float64)
        # kind == "plate": the half side, the angular positions around the hole
        # and the cells per 90 degree block (the hole radius is ``radius``).
        self.half = 1.0
        self.period = 1.0
        self.nt = 1.0
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
    else:  # unit
        out[:, 0] = cx + r * p[:, 0]
        out[:, 1] = cy + r * p[:, 1]
        out[:, 2] = cz + r * p[:, 2]
    return out


def _mid_angle(a, da, b, db):
    two_pi = 2.0 * _PI
    d = b - a
    d = d - two_pi * np.floor(d / two_pi + 0.5)
    mid = a + 0.5 * d
    return np.where(da, b, np.where(db, a, mid))


def _mid(soup, lo, hi):
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


def _finish(soup, order, specs, facet_dim):
    if order == 2:
        _elevate(soup)
    points = _positions(soup)
    blocks = [(t, c) for t, c in soup.blocks if c.size]
    regions = _regions(soup, blocks, specs, facet_dim)
    return Mesh(points, blocks, regions=regions)


# --------------------------------------------------------------------------- #
# rectangle                                                                   #
# --------------------------------------------------------------------------- #


def _rectangle_py(dims, size, origin, cell_type, split, order):
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
    return _finish(soup, order, [("xmin", 1), ("xmax", 2), ("ymin", 4), ("ymax", 8)], 1)


def rectangle(
    dims=(1, 1),
    size=(1.0, 1.0),
    origin=(0.0, 0.0, 0.0),
    cell_type="quad",
    split="zigzag",
    order=1,
    max_cells=DEFAULT_MAX_CELLS,
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
        1, or 2 for ``quad8`` / ``triangle6`` cells.
    max_cells :
        Refuse by name above this many cells; non-positive disables the check.

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
        )
    except Exception as exc:
        if not core_op_declined(exc, "rectangle"):
            raise
        return _rectangle_py(dims, size, origin, cell_type, split, order)


# --------------------------------------------------------------------------- #
# box                                                                         #
# --------------------------------------------------------------------------- #


def _box_py(dims, size, origin, cell_type, order):
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
    return _finish(soup, order, names, 2)


def box(
    dims=(1, 1, 1),
    size=(1.0, 1.0, 1.0),
    origin=(0.0, 0.0, 0.0),
    cell_type="hexahedron",
    order=1,
    max_cells=DEFAULT_MAX_CELLS,
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
        1, or 2 for ``hexahedron20`` / ``tetra10`` / ``wedge15`` cells.
    max_cells :
        Refuse by name above this many cells; non-positive disables the check.

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
        )
    except Exception as exc:
        if not core_op_declined(exc, "box"):
            raise
        return _box_py(dims, size, origin, cell_type, order)


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


def _disk_py(radius, n, center, cell_type, order):
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
    return _finish(soup, order, [("rim", 4)], 1)


def disk(
    radius=1.0,
    n=8,
    center=(0.0, 0.0, 0.0),
    cell_type="quad",
    order=1,
    max_cells=DEFAULT_MAX_CELLS,
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
        1, or 2 for ``quad8`` / ``triangle6`` cells with a curved rim.
    max_cells :
        Refuse by name above this many cells; non-positive disables the check.

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
        )
    except Exception as exc:
        if not core_op_declined(exc, "disk"):
            raise
        return _disk_py(float(radius), int(n), center, cell_type, order)


# --------------------------------------------------------------------------- #
# uv_sphere                                                                   #
# --------------------------------------------------------------------------- #


def _uv_sphere_py(radius, nu, nv, center, cell_type, order):
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
    return _finish(soup, order, [], 2)


def uv_sphere(
    radius=1.0,
    nu=16,
    nv=8,
    center=(0.0, 0.0, 0.0),
    cell_type="triangle",
    order=1,
    max_cells=DEFAULT_MAX_CELLS,
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
        1, or 2 for ``triangle6`` / ``quad8`` cells on the exact sphere.
    max_cells :
        Refuse by name above this many cells; non-positive disables the check.
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
        )
    except Exception as exc:
        if not core_op_declined(exc, "uv_sphere"):
            raise
        return _uv_sphere_py(float(radius), int(nu), int(nv), center, cell_type, order)


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


def _icosa_sphere_py(radius, subdivisions, center, order):
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
    return _finish(soup, order, [], 2)


def icosa_sphere(
    radius=1.0,
    subdivisions=2,
    center=(0.0, 0.0, 0.0),
    order=1,
    max_cells=DEFAULT_MAX_CELLS,
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
        1, or 2 for ``triangle6`` cells on the exact sphere.
    max_cells :
        Refuse by name above this many cells; non-positive disables the check.
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
        )
    except Exception as exc:
        if not core_op_declined(exc, "icosa_sphere"):
            raise
        return _icosa_sphere_py(float(radius), int(subdivisions), center, order)


# --------------------------------------------------------------------------- #
# cylinder                                                                    #
# --------------------------------------------------------------------------- #


def _cylinder_py(radius, height, nxy, nz, center, cell_type, order):
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
    return _finish(soup, order, [("side", 4), ("top", 1), ("bottom", 2)], 2)


def cylinder(
    radius=1.0,
    height=1.0,
    nxy=8,
    nz=8,
    center=(0.0, 0.0, 0.0),
    cell_type="hexahedron",
    order=1,
    max_cells=DEFAULT_MAX_CELLS,
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
        1, or 2 for ``hexahedron20`` / ``wedge15`` cells with a curved side.
    max_cells :
        Refuse by name above this many cells; non-positive disables the check.

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
        )
    except Exception as exc:
        if not core_op_declined(exc, "cylinder"):
            raise
        return _cylinder_py(
            float(radius), float(height), int(nxy), int(nz), center, cell_type, order
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


def _reference_triangle_py(n, size, origin, order):
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
    return _finish(soup, order, [("ymin", 1), ("xmin", 2), ("hypotenuse", 4)], 1)


def reference_triangle(
    n=4,
    size=(1.0, 1.0),
    origin=(0.0, 0.0, 0.0),
    order=1,
    max_cells=DEFAULT_MAX_CELLS,
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
        1, or 2 for ``triangle6`` cells.
    max_cells :
        Refuse by name above this many cells; non-positive disables the check.

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
        )
    except Exception as exc:
        if not core_op_declined(exc, "reference_triangle"):
            raise
        return _reference_triangle_py(int(n), size, origin, order)


# --------------------------------------------------------------------------- #
# ngon                                                                        #
# --------------------------------------------------------------------------- #


def _ngon_py(sides, radius, rings, center, order):
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
    return _finish(soup, order, [("rim", 1)], 1)


def ngon(
    sides=6,
    radius=1.0,
    rings=1,
    center=(0.0, 0.0, 0.0),
    order=1,
    max_cells=DEFAULT_MAX_CELLS,
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
        1, or 2 for ``triangle6`` cells (the sides are straight).
    max_cells :
        Refuse by name above this many cells; non-positive disables the check.

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
        )
    except Exception as exc:
        if not core_op_declined(exc, "ngon"):
            raise
        return _ngon_py(int(sides), float(radius), int(rings), center, order)


# --------------------------------------------------------------------------- #
# annulus                                                                     #
# --------------------------------------------------------------------------- #


def _annulus_py(inner, outer, ntheta, nr, angle, center, cell_type, order):
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
    return _finish(soup, order, specs, 1)


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
        1, or 2 for ``quad8`` / ``triangle6`` cells with curved rims.
    max_cells :
        Refuse by name above this many cells; non-positive disables the check.

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
        )


# --------------------------------------------------------------------------- #
# lshape                                                                      #
# --------------------------------------------------------------------------- #


def _lshape_py(n, size, center, cell_type, order):
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
    )


def lshape(
    n=4,
    size=2.0,
    center=(0.0, 0.0, 0.0),
    cell_type="quad",
    order=1,
    max_cells=DEFAULT_MAX_CELLS,
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
        1, or 2 for ``quad8`` / ``triangle6`` cells.
    max_cells :
        Refuse by name above this many cells; non-positive disables the check.

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
        )
    except Exception as exc:
        if not core_op_declined(exc, "lshape"):
            raise
        return _lshape_py(int(n), float(size), center, cell_type, order)


# --------------------------------------------------------------------------- #
# plate_with_hole                                                             #
# --------------------------------------------------------------------------- #

_PLATE_SIDE_BIT = np.array([4, 16, 2, 8], dtype=np.int64)  # xmax, ymax, xmin, ymin


def _plate_py(radius, size, nt, nr, center, cell_type, order):
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
        1, or 2 for ``quad8`` / ``triangle6`` cells with a circular hole.
    max_cells :
        Refuse by name above this many cells; non-positive disables the check.

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
        )
    except Exception as exc:
        if not core_op_declined(exc, "plate_with_hole"):
            raise
        return _plate_py(radius, size, int(nt), int(nr), center, cell_type, order)


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

"""Octree (3-D) and quadtree (2-D) leaves around a skin.

Leaves are fine only where the skin is. Every leaf is a ``hexahedron`` (``quad`` in
2-D) on the integer lattice of the finest level, so shared corners are
deduplicated arithmetically -- no hash, no tolerance -- and coordinates are
``origin + spacing * index``, never accumulated. Leaves come out in Morton order.

This reference implementation is a byte-for-byte twin of the C++ core
(``operations/octree.cpp``), pinned by
``tests/python/test_octree.py::test_cpp_matches_python``.

Public API:

* :func:`octree` -- build the tree around a mesh.
"""

from __future__ import annotations

import numpy as np

from . import _planar_skin
from ._fallback import core_op_declined
from ._grid import DEFAULT_MAX_CELLS
from ._mesh import Mesh
from ._sdf import _sample_py, _soup
from ._voxelize import _tri_box_overlap

__all__ = ["octree"]

_PREFIX = "meshio++: octree: "
_FILLS = ("all", "surface", "inside")
_BALANCES = ("none", "face", "full")
_MAX_LEVEL = 20

_CORNER = (
    (0, 0, 0),
    (1, 0, 0),
    (1, 1, 0),
    (0, 1, 0),
    (0, 0, 1),
    (1, 0, 1),
    (1, 1, 1),
    (0, 1, 1),
)


def _spread(v):
    v = int(v) & 0x1FFFFF
    v = (v | (v << 32)) & 0x1F00000000FFFF
    v = (v | (v << 16)) & 0x1F0000FF0000FF
    v = (v | (v << 8)) & 0x100F00F00F00F00F
    v = (v | (v << 4)) & 0x10C30C30C30C30C3
    v = (v | (v << 2)) & 0x1249249249249249
    return v


def _morton(x, y, z):
    return _spread(x) | (_spread(y) << 1) | (_spread(z) << 2)


class _Grid:
    pass


def _prim_hits(g, prim, centre, half):
    if g.two_d:
        return _planar_skin.seg_box(centre, half, prim[0], prim[1])
    return bool(_tri_box_overlap(centre, half, prim[0], prim[1], prim[2]))


def _box(g, x, y, z, level):
    s = 1 << (g.finest - level)
    idx = (x, y, z)
    centre = np.zeros(3)
    half = np.zeros(3)
    for k in range(3):
        if k == 2 and g.two_d:
            centre[2] = g.origin[2]
            half[2] = 0.0
            continue
        centre[k] = g.origin[k] + float(2 * idx[k] + s) * g.half_fine
        half[k] = float(s) * g.half_fine
    return centre, half


def _make_leaf(g, prims, parent_cand, x, y, z, level):
    centre, half = _box(g, x, y, z, level)
    inflated = half + g.band
    cand = []
    cut = False
    for p in parent_cand:
        if not _prim_hits(g, prims[p], centre, inflated):
            continue
        cand.append(p)
        if not cut and _prim_hits(g, prims[p], centre, half):
            cut = True
    return [x, y, z, level, cut, cand]


def _split(g, prims, leaf, out):
    half = 1 << (g.finest - leaf[3] - 1)
    nz = 1 if g.two_d else 2
    for cz in range(nz):
        for cy in range(2):
            for cx in range(2):
                out.append(
                    _make_leaf(
                        g,
                        prims,
                        leaf[5],
                        leaf[0] + cx * half,
                        leaf[1] + cy * half,
                        leaf[2] + cz * half,
                        leaf[3] + 1,
                    )
                )


def _check_budget(n, max_cells):
    if max_cells > 0 and n > max_cells:
        raise ValueError(
            f"{_PREFIX}the tree has {n} leaves, above the limit of {max_cells} "
            "(raise max_cells, lower max_level, or narrow the band)"
        )


def _octree_py(
    mesh,
    root_cells,
    min_level,
    max_level,
    bounds,
    padding,
    padding_relative,
    dimension,
    sample_distance,
    balance,
    fill,
    attach_distance,
    even_odd,
    max_cells,
    sign,
    transition="hanging",
    snap=False,
    surface_level=None,
    curvature_threshold=0.0,
    feature_angle=0.0,
    size_array="",
    plane="auto",
):
    if (
        transition != "hanging"
        or snap
        or (surface_level is not None and surface_level >= 0)
        or curvature_threshold
        or feature_angle
        or size_array
    ):
        raise NotImplementedError(
            f"{_PREFIX}transition other than 'hanging', snap, surface_level, "
            "curvature_threshold, feature_angle and size_array need the compiled "
            "core; the numpy reference implements the hanging-node octree only"
        )
    if min_level < 0 or max_level < min_level or max_level > _MAX_LEVEL:
        raise ValueError(
            f"{_PREFIX}need 0 <= min_level <= max_level <= {_MAX_LEVEL}, got "
            f"min_level {min_level} and max_level {max_level}"
        )
    if sample_distance < 0.0:
        raise ValueError(f"{_PREFIX}sample_distance must not be negative")
    if dimension not in (0, 2, 3):
        raise ValueError(f"{_PREFIX}dimension must be 2, 3, or 0 (auto)")

    g = _Grid()
    g.finest = int(max_level)
    g.band = float(sample_distance)
    roots = [int(v) for v in root_cells]

    two_d = dimension == 2
    if dimension == 0:
        two_d = _planar_skin.polyline_only(mesh)
    g.two_d = two_d
    if two_d:
        roots[2] = 1
    for k in range(3):
        if roots[k] <= 0:
            raise ValueError(
                f"{_PREFIX}root_cells must be positive, got {roots[k]} on axis {k}"
            )
        if (roots[k] << g.finest) >= (1 << 21):
            raise ValueError(
                f"{_PREFIX}the finest lattice is too fine: root_cells << max_level "
                "must stay below 2^21 on every axis"
            )
    g.roots = roots
    extent = [roots[k] << g.finest for k in range(3)]

    # ---- the skin -----------------------------------------------------------
    prims = []
    axes = (0, 1, 2)
    if two_d:
        axes = _planar_skin.frame(mesh, plane)
        prims, seg_ids, has_lines = _planar_skin.collect(mesh, _PREFIX, axes)
        if (fill == "inside" or attach_distance) and has_lines:
            _planar_skin.require_closed(seg_ids, _PREFIX)
    else:
        _pts, _verts, corners, _source = _soup(mesh)
        prims = [(c[0], c[1], c[2]) for c in corners]

    # ---- the root lattice -----------------------------------------------------
    if bounds is not None:
        b = np.asarray(bounds, dtype=np.float64).reshape(-1)
        if b.size != 6:
            raise ValueError(f"{_PREFIX}bounds must be (xlo, ylo, zlo, xhi, yhi, zhi)")
        lo, hi = b[:3].copy(), b[3:].copy()
        for k in range(3):
            if not hi[k] >= lo[k]:
                raise ValueError(f"{_PREFIX}bounds are inverted on axis {k}")
    else:
        pts = np.asarray(mesh.points, dtype=np.float64)
        if len(pts) == 0:
            raise ValueError(
                f"{_PREFIX}the mesh has no points, so it has no bounding box to cover "
                "(pass explicit bounds)"
            )
        if pts.shape[1] == 2:
            pts = np.column_stack([pts, np.zeros(len(pts))])
        lo, hi = pts.min(axis=0), pts.max(axis=0)
    if two_d:  # world -> the skin's frame
        lo, hi = (
            _planar_skin.permute(lo, axes).copy(),
            _planar_skin.permute(hi, axes).copy(),
        )
    diag = float(np.sqrt(np.sum((hi - lo) ** 2)))
    pad = padding + padding_relative * diag
    if pad < 0.0:
        raise ValueError(f"{_PREFIX}padding is negative")
    for k in range(2 if two_d else 3):
        lo[k] -= pad
        hi[k] += pad
    root = 0.0
    for k in range(2 if two_d else 3):
        h = (hi[k] - lo[k]) / float(roots[k])
        root = h if h > root else root
    if not root > 0.0:
        raise ValueError(
            f"{_PREFIX}the bounding box is degenerate, so no root cell can cover it "
            "(pass explicit bounds or padding)"
        )
    g.origin = lo
    g.root = root
    g.fine = root / float(1 << g.finest)
    g.half_fine = g.fine * 0.5

    # ---- refinement, level by level --------------------------------------------
    all_cand = list(range(len(prims)))
    current = []
    for k in range(roots[2]):
        for j in range(roots[1]):
            for i in range(roots[0]):
                current.append(
                    _make_leaf(
                        g,
                        prims,
                        all_cand,
                        i << g.finest,
                        j << g.finest,
                        k << g.finest,
                        0,
                    )
                )
    _check_budget(len(current), max_cells)
    done = []
    for level in range(g.finest):
        nxt = []
        for leaf in current:
            if level < min_level or leaf[5]:
                _split(g, prims, leaf, nxt)
            else:
                done.append(leaf)
        _check_budget(len(done) + len(nxt), max_cells)
        current = nxt
    leaves = done + current

    # ---- 2:1 balance ------------------------------------------------------------
    if balance != "none" and g.finest >= 2:
        full = balance == "full"
        nzr = 0 if two_d else 1
        offsets = []
        for dz in range(-nzr, nzr + 1):
            for dy in (-1, 0, 1):
                for dx in (-1, 0, 1):
                    nonzero = (dx != 0) + (dy != 0) + (dz != 0)
                    if nonzero == 0 or (not full and nonzero > 1):
                        continue
                    offsets.append((dx, dy, dz))
        while True:
            maps = [dict() for _ in range(g.finest + 1)]
            for n, leaf in enumerate(leaves):
                maps[leaf[3]][(leaf[0], leaf[1], leaf[2])] = n
            flag = [False] * len(leaves)
            for leaf in leaves:
                if leaf[3] < 2:
                    continue
                s = 1 << (g.finest - leaf[3])
                for dx, dy, dz in offsets:
                    nx, ny, nz = leaf[0] + dx * s, leaf[1] + dy * s, leaf[2] + dz * s
                    if (
                        nx < 0
                        or ny < 0
                        or nz < 0
                        or nx >= extent[0]
                        or ny >= extent[1]
                        or nz >= extent[2]
                    ):
                        continue
                    for lev in range(0, leaf[3] - 1):
                        sh = g.finest - lev
                        hit = maps[lev].get(
                            ((nx >> sh) << sh, (ny >> sh) << sh, (nz >> sh) << sh)
                        )
                        if hit is not None:
                            flag[hit] = True
                            break
            if not any(flag):
                break
            nxt = []
            for n, leaf in enumerate(leaves):
                if flag[n]:
                    _split(g, prims, leaf, nxt)
                else:
                    nxt.append(leaf)
            _check_budget(len(nxt), max_cells)
            leaves = nxt

    # ---- Morton order -----------------------------------------------------------
    leaves.sort(key=lambda leaf: _morton(leaf[0], leaf[1], leaf[2]))
    nleaves = len(leaves)

    # ---- fill rule and distance --------------------------------------------------
    centres = np.zeros((nleaves, 3))
    for n, leaf in enumerate(leaves):
        centres[n] = _box(g, leaf[0], leaf[1], leaf[2], leaf[3])[0]
    distance = None
    if fill == "inside" or attach_distance:
        if two_d:
            distance = np.array(
                [
                    _planar_skin.signed_distance(prims, centres[n], even_odd)
                    for n in range(nleaves)
                ]
            )
        else:
            if fill == "inside" and sign == "unsigned":
                raise ValueError(
                    f"{_PREFIX}fill 'inside' needs a sign, but sign='unsigned' was given"
                )
            distance = _sample_py(mesh, centres, sign, "angle", 0.0)[0]

    kept = []
    for n, leaf in enumerate(leaves):
        if fill == "surface":
            keep = leaf[4]
        elif fill == "inside":
            keep = distance[n] < 0.0
        else:
            keep = True
        if keep:
            kept.append(n)

    report = {
        "dimension": 2 if two_d else 3,
        "origin": [float(v) for v in g.origin],
        "root_size": float(root),
        "num_leaves": len(kept),
        "num_hanging": 0,
        "finest_level": max((leaves[n][3] for n in kept), default=0),
        "num_snapped": 0,
        "num_refused": 0,
    }
    if not kept:
        return Mesh(np.zeros((0, 3)), []), report

    # ---- points -------------------------------------------------------------------
    nc = 4 if two_d else 8
    px, py = extent[0] + 1, extent[1] + 1
    if not two_d and px > (1 << 62) // py // (extent[2] + 1):
        raise ValueError(
            f"{_PREFIX}the finest lattice is too large for 63-bit point keys"
        )
    keys = np.zeros((len(kept), nc), dtype=np.int64)
    for i, n in enumerate(kept):
        leaf = leaves[n]
        s = 1 << (g.finest - leaf[3])
        for c in range(nc):
            x = leaf[0] + _CORNER[c][0] * s
            y = leaf[1] + _CORNER[c][1] * s
            z = leaf[2] + _CORNER[c][2] * s
            keys[i, c] = (z * py + y) * px + x
    uniq = np.unique(keys)
    ux = uniq % px
    uy = (uniq // px) % py
    uz = uniq // (px * py)
    points = np.empty((len(uniq), 3), dtype=np.float64)
    points[:, 0] = g.origin[0] + ux.astype(np.float64) * g.fine
    points[:, 1] = g.origin[1] + uy.astype(np.float64) * g.fine
    if two_d:
        points[:, 2] = g.origin[2]
    else:
        points[:, 2] = g.origin[2] + uz.astype(np.float64) * g.fine
    conn = np.searchsorted(uniq, keys).astype(np.int64)

    # ---- hanging nodes, relative to the kept leaves ---------------------------------
    kept_maps = [dict() for _ in range(g.finest + 1)]
    for leaf_n in kept:
        leaf = leaves[leaf_n]
        kept_maps[leaf[3]][(leaf[0], leaf[1], leaf[2])] = True
    hanging = np.zeros(len(uniq), dtype=np.int64)
    nzo = 1 if two_d else 2
    for p in range(len(uniq)):
        pt = (int(ux[p]), int(uy[p]), int(uz[p]))
        flag = 0
        for oz in range(nzo):
            for oy in range(2):
                for ox in range(2):
                    if flag:
                        break
                    c = (pt[0] - 1 + ox, pt[1] - 1 + oy, 0 if two_d else pt[2] - 1 + oz)
                    if any(c[k] < 0 or c[k] >= extent[k] for k in range(3)):
                        continue
                    for lev in range(g.finest, -1, -1):
                        sh = g.finest - lev
                        a = ((c[0] >> sh) << sh, (c[1] >> sh) << sh, (c[2] >> sh) << sh)
                        if a not in kept_maps[lev]:
                            continue
                        s = 1 << sh
                        for k in range(2 if two_d else 3):
                            if pt[k] != a[k] and pt[k] != a[k] + s:
                                flag = 1
                        break
        hanging[p] = flag
    report["num_hanging"] = int(hanging.sum())

    if two_d:
        points = _planar_skin.unpermute(points, axes)
    out = Mesh(points, [("quad" if two_d else "hexahedron", conn)])
    out.cell_data["octree:level"] = [
        np.array([leaves[n][3] for n in kept], dtype=np.int64)
    ]
    out.cell_data["octree:cut"] = [
        np.array([1 if leaves[n][4] else 0 for n in kept], dtype=np.int64)
    ]
    if attach_distance:
        out.cell_data["sdf:distance"] = [
            np.array([distance[n] for n in kept], dtype=np.float64)
        ]
    out.point_data["octree:hanging"] = hanging
    return out, report


def octree(
    mesh,
    root_cells=(1, 1, 1),
    min_level=0,
    max_level=4,
    bounds=None,
    padding=0.0,
    padding_relative=0.0,
    dimension=0,
    sample_distance=0.0,
    balance="face",
    fill="all",
    attach_distance=False,
    even_odd=False,
    max_cells=DEFAULT_MAX_CELLS,
    sign="pseudonormal",
    watertight_check="warn",
    return_report=False,
    transition="hanging",
    snap=False,
    surface_level=None,
    curvature_threshold=0.0,
    feature_angle=0.0,
    size_array=None,
    plane="auto",
):
    """Build an octree (quadtree in 2-D) around the skin ``mesh``.

    Leaves are split while the skin cuts them (or lies within ``sample_distance``
    of them) and their level is below ``max_level``; every leaf is split down to
    ``min_level``. A 2:1 ``balance`` ("none", "face" or "full") follows. ``fill``
    keeps ``voxelize``'s meaning ("all", "surface", "inside").

    ``transition`` is ``"hanging"`` (plain cells with hanging nodes, the default),
    ``"polyhedra"`` (conforming split-hex ``polyhedron``/``polygon`` cells) or
    ``"templates"`` (conforming hex/pyramid/tetra, quad/triangle in 2-D). ``snap``
    moves the outer nodes onto the skin without inverting a leaf (needs ``fill``
    ``"surface"`` or ``"inside"``). ``surface_level`` stops cut leaves early;
    ``curvature_threshold`` (|mean curvature|), ``feature_angle`` (degrees) and
    ``size_array`` (a ``point_data`` array of target edge lengths) refine further,
    up to ``max_level``. ``plane`` (``"auto"``, ``"xy"``, ``"xz"``, ``"yz"``)
    chooses the plane of a 2-D skin. These options, and 3-D ``plane``, need the
    compiled core; the report adds ``num_snapped`` and ``num_refused``.

    The result carries ``octree:level`` and ``octree:cut`` cell data and the
    ``octree:hanging`` point data; ``attach_distance`` adds ``sdf:distance``.
    ``min_level == max_level`` is a uniform voxel grid.
    """
    if fill not in _FILLS:
        raise ValueError(
            f"{_PREFIX}unknown fill '{fill}' (expected one of: {', '.join(_FILLS)})"
        )
    if balance not in _BALANCES:
        raise ValueError(
            f"{_PREFIX}unknown balance '{balance}' (expected one of: "
            f"{', '.join(_BALANCES)})"
        )
    if np.ndim(root_cells) == 0:
        root_cells = (int(root_cells),) * 3
    root_cells = [int(v) for v in root_cells]
    if len(root_cells) == 2:
        root_cells.append(1)
    if len(root_cells) != 3:
        raise ValueError(f"{_PREFIX}root_cells must be one, two or three counts")

    out = None
    report = None
    try:
        from . import _core

        res = _core.octree(
            mesh,
            root_cells,
            int(min_level),
            int(max_level),
            None if bounds is None else [float(v) for v in bounds],
            float(padding),
            float(padding_relative),
            int(dimension),
            float(sample_distance),
            balance,
            fill,
            bool(attach_distance),
            bool(even_odd),
            int(max_cells),
            sign,
            watertight_check,
            transition,
            bool(snap),
            -1 if surface_level is None else int(surface_level),
            float(curvature_threshold),
            float(feature_angle),
            size_array or "",
            plane,
        )
        out = res["mesh"]
        report = {
            "dimension": res["dimension"],
            "origin": res["origin"],
            "root_size": res["root_size"],
            "num_leaves": res["num_leaves"],
            "num_hanging": res["num_hanging"],
            "finest_level": res["finest_level"],
            "num_snapped": res["num_snapped"],
            "num_refused": res["num_refused"],
        }
    except Exception as exc:
        if not core_op_declined(exc, "octree"):
            raise
        out = None

    if out is None:
        out, report = _octree_py(
            mesh,
            root_cells,
            int(min_level),
            int(max_level),
            bounds,
            float(padding),
            float(padding_relative),
            int(dimension),
            float(sample_distance),
            balance,
            fill,
            bool(attach_distance),
            bool(even_odd),
            int(max_cells),
            sign,
            transition,
            bool(snap),
            surface_level,
            float(curvature_threshold),
            float(feature_angle),
            size_array or "",
            plane,
        )

    return (out, report) if return_report else out

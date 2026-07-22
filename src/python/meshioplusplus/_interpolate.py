"""Cross-mesh field transfer: sample a source mesh's data onto a target mesh.

``interpolate(source, target, ...)`` returns a new mesh that is a copy of the
target — geometry, connectivity, and the target's own data preserved exactly —
with the requested source arrays sampled onto it: source ``point_data`` at the
target's points, source ``cell_data`` at the target's cell centroids (always
by nearest source-cell centroid, whatever the method — a piecewise-constant
field has no meaningful barycentric interpolant). Cross-location transfer is
out of scope; compose with the ``data to-cell`` / ``to-point`` operations.

Methods:

- ``"nearest"`` (default): each target sample takes the nearest source point's
  value (bucket-grid spatial hash, expanding-shell search; tie-break smallest
  squared distance then lowest source index). Dtype-preserving.
- ``"barycentric"``: the source is simplexified first (``convert_cells``), so
  on a quad/hex source this is the simplex-linear interpolant, not true
  bi-/trilinear; triangles are evaluated in the xy-plane (planar assumption —
  use ``"nearest"`` for a curved surface embedded in 3D). Exact on a linear
  field. A target point outside the source domain receives ``default_value``
  unless ``extrapolate=True`` (nearest source point's value). Always Float64.

The numpy fallback below replicates the C++ core's accumulation and tie-break
order expression for expression — grid cell size without ``cbrt`` (an integer
loop plus one division), per-candidate distances summed axis by axis,
barycentric weights applied corner by corner — so the two lanes produce
byte-identical output (pinned by ``tests/python/test_interpolate.py``'s
``test_cpp_matches_python``).
"""

import math

import numpy as np

from ._mesh import Mesh

_METHODS = ("nearest", "barycentric")
_CONFLICTS = ("error", "overwrite", "suffix")

# The containment tolerance on barycentric weights (dimensionless volume
# ratios, hence scale-invariant). The same literal appears in interpolate.cpp.
_BARY_EPS = 1e-12


def _coords3(points):
    """The coordinates widened to Float64 and padded to three columns."""
    pts = np.asarray(points, dtype=np.float64)
    if pts.ndim != 2:
        pts = pts.reshape(0, 3) if pts.size == 0 else pts.reshape(len(pts), -1)
    n, dim = pts.shape
    if dim >= 3:
        return np.ascontiguousarray(pts[:, :3])
    out = np.zeros((n, 3), dtype=np.float64)
    out[:, :dim] = pts
    return out


def _cell_size(lo, hi, n):
    """The grid cell size for ``n`` items over the bbox ``[lo, hi]``.

    Largest axis extent divided by the smallest integer R with R**3 >= n — no
    cbrt/pow, so the value is bit-identical to the C++ core's.
    """
    max_extent = 0.0
    for d in range(3):
        e = float(hi[d]) - float(lo[d])
        if e > max_extent:
            max_extent = e
    r = 1
    while r * r * r < n:
        r += 1
    return max_extent / float(r) if max_extent > 0.0 else 1.0


def _build_grid(coords3, cell):
    """Bucket grid dict ``(kx, ky, kz) -> [ids ascending]`` plus the key bbox."""
    keys = np.floor(coords3 / cell).astype(np.int64)
    grid = {}
    for i in range(keys.shape[0]):
        grid.setdefault((int(keys[i, 0]), int(keys[i, 1]), int(keys[i, 2])), []).append(
            i
        )
    klo = (int(keys[:, 0].min()), int(keys[:, 1].min()), int(keys[:, 2].min()))
    khi = (int(keys[:, 0].max()), int(keys[:, 1].max()), int(keys[:, 2].max()))
    return grid, klo, khi


def _shell_keys(k, r, klo, khi):
    """The grid keys at Chebyshev radius exactly ``r`` around ``k``, clamped to
    the occupied key bbox (cells outside it are empty by construction)."""
    zl, zh = max(-r, klo[2] - k[2]), min(r, khi[2] - k[2])
    yl, yh = max(-r, klo[1] - k[1]), min(r, khi[1] - k[1])
    xl, xh = max(-r, klo[0] - k[0]), min(r, khi[0] - k[0])
    if zl > zh or yl > yh or xl > xh:
        return
    if r == 0:
        yield (k[0], k[1], k[2])
        return
    for dz in range(zl, zh + 1):
        if dz == -r or dz == r:
            for dy in range(yl, yh + 1):
                for dx in range(xl, xh + 1):
                    yield (k[0] + dx, k[1] + dy, k[2] + dz)
        else:
            for dy in range(yl, yh + 1):
                if dy == -r or dy == r:
                    for dx in range(xl, xh + 1):
                        yield (k[0] + dx, k[1] + dy, k[2] + dz)
                else:
                    if xl == -r:
                        yield (k[0] - r, k[1] + dy, k[2] + dz)
                    if xh == r:
                        yield (k[0] + r, k[1] + dy, k[2] + dz)


def _nearest(grid, klo, khi, cell, src3, p):
    """Nearest id to ``p`` — the C++ expanding-shell search, expression for
    expression (stop once ``((r-1)*cell)**2`` exceeds the best squared
    distance; tie-break smallest d2 then lowest id)."""
    k = (
        int(math.floor(p[0] / cell)),
        int(math.floor(p[1] / cell)),
        int(math.floor(p[2] / cell)),
    )

    def axis_cap(v, lo, hi):
        return max(abs(v - lo), abs(v - hi))

    def axis_start(v, lo, hi):
        if v < lo:
            return lo - v
        if v > hi:
            return v - hi
        return 0

    rcap = max(axis_cap(k[d], klo[d], khi[d]) for d in range(3))
    r = max(axis_start(k[d], klo[d], khi[d]) for d in range(3))

    best_d2 = math.inf
    best = -1
    while True:
        if best >= 0 and r >= 1:
            ring = float(r - 1) * cell
            if ring * ring > best_d2:
                break
        if r > rcap:
            break
        for key in _shell_keys(k, r, klo, khi):
            ids = grid.get(key)
            if ids is None:
                continue
            idx = np.asarray(ids, dtype=np.int64)
            delta = p[0] - src3[idx, 0]
            d2 = delta * delta
            delta = p[1] - src3[idx, 1]
            d2 = d2 + delta * delta
            delta = p[2] - src3[idx, 2]
            d2 = d2 + delta * delta
            # ids ascend, so the first minimum is the lowest id among ties.
            j = int(np.argmin(d2))
            if d2[j] < best_d2 or (d2[j] == best_d2 and ids[j] < best):
                best_d2 = float(d2[j])
                best = int(ids[j])
        r += 1
    return best


def _sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _tp(a, b, c):
    """``a . (b x c)`` — detail::triple_product's exact expression order."""
    x0 = b[1] * c[2] - b[2] * c[1]
    x1 = b[2] * c[0] - b[0] * c[2]
    x2 = b[0] * c[1] - b[1] * c[0]
    return a[0] * x0 + a[1] * x1 + a[2] * x2


def _area2(p0, p1, p2):
    """Twice the signed area of p0-p1-p2 in the xy-plane (z ignored)."""
    return (p1[0] - p0[0]) * (p2[1] - p0[1]) - (p1[1] - p0[1]) * (p2[0] - p0[0])


def _tet_weights(a, b, c, d, p):
    v = _tp(_sub(b, a), _sub(c, a), _sub(d, a))
    if v == 0.0:
        return None
    wa = _tp(_sub(b, p), _sub(c, p), _sub(d, p)) / v
    wb = _tp(_sub(p, a), _sub(c, a), _sub(d, a)) / v
    wc = _tp(_sub(b, a), _sub(p, a), _sub(d, a)) / v
    wd = _tp(_sub(b, a), _sub(c, a), _sub(p, a)) / v
    if wa < -_BARY_EPS or wb < -_BARY_EPS or wc < -_BARY_EPS or wd < -_BARY_EPS:
        return None
    return (wa, wb, wc, wd)


def _tri_weights(a, b, c, p):
    v = _area2(a, b, c)
    if v == 0.0:
        return None
    wa = _area2(p, b, c) / v
    wb = _area2(a, p, c) / v
    wc = _area2(a, b, p) / v
    if wa < -_BARY_EPS or wb < -_BARY_EPS or wc < -_BARY_EPS:
        return None
    return (wa, wb, wc)


def _resolve_names(source, target, arrays):
    """The (point_names, cell_names) to transfer — the C++ selection rule."""
    if arrays is None or len(arrays) == 0:
        return sorted(source.point_data), []
    point_names, cell_names = [], []
    for name in arrays:
        any_hit = False
        if name in source.point_data:
            point_names.append(name)
            any_hit = True
        if name in source.cell_data:
            cell_names.append(name)
            any_hit = True
        if not any_hit:
            raise ValueError(
                f"meshio++: interpolate: no source point_data or cell_data array named "
                f"{name!r} (available point_data: "
                f"{', '.join(sorted(source.point_data)) or 'none'}; cell_data: "
                f"{', '.join(sorted(source.cell_data)) or 'none'})"
            )
    return point_names, cell_names


def _resolve_out_name(existing, name, on_conflict):
    if name not in existing:
        return name
    if on_conflict == "overwrite":
        return name
    if on_conflict == "suffix":
        suffixed = name + "_interp"
        if suffixed in existing:
            raise ValueError(
                f"meshio++: interpolate: array {name!r} already exists on the target "
                f"and so does {suffixed!r} (on_conflict='suffix')"
            )
        return suffixed
    raise ValueError(
        f"meshio++: interpolate: array {name!r} already exists on the target "
        "(on_conflict='error'; pass 'overwrite' or 'suffix')"
    )


def _simplex_conn(smesh, want, nv):
    rows = []
    for block in smesh.cells:
        if block.type != want or len(block.data) == 0:
            continue
        rows.append(np.asarray(block.data, dtype=np.int64).reshape(-1, nv))
    return np.concatenate(rows, axis=0) if rows else np.zeros((0, nv), dtype=np.int64)


def _fill_simplex_grid(pts3, sconn, cell):
    """Bucket grid over the simplices' quantized bounding boxes."""
    crd = pts3[sconn]  # (nsimp, nv, 3)
    lo = np.floor(crd.min(axis=1) / cell).astype(np.int64)
    hi = np.floor(crd.max(axis=1) / cell).astype(np.int64)
    grid = {}
    for s in range(sconn.shape[0]):
        for z in range(int(lo[s, 2]), int(hi[s, 2]) + 1):
            for y in range(int(lo[s, 1]), int(hi[s, 1]) + 1):
                for x in range(int(lo[s, 0]), int(hi[s, 0]) + 1):
                    grid.setdefault((x, y, z), []).append(s)
    return grid


def _interpolate_py(
    source,
    target,
    method="nearest",
    arrays=None,
    extrapolate=False,
    default_value=0.0,
    on_conflict="error",
):
    """Pure-numpy reference implementation of :func:`interpolate`.

    Byte-identical to the C++ core by construction (see the module docstring);
    raises ``NotImplementedError`` only where ``_partition._centroids`` does —
    a ragged/polyhedron source or target block in a cell_data transfer.
    """
    if method not in _METHODS:
        raise ValueError(
            f"meshio++: interpolate: unknown method {method!r} "
            "(expected 'nearest' or 'barycentric')"
        )
    if on_conflict not in _CONFLICTS:
        raise ValueError(
            f"meshio++: interpolate: unknown on_conflict {on_conflict!r} "
            "(expected 'error', 'overwrite' or 'suffix')"
        )
    if len(np.asarray(source.points)) == 0:
        raise ValueError("meshio++: interpolate: the source mesh has no points")

    point_names, cell_names = _resolve_names(source, target, arrays)
    for name in cell_names:
        blocks = source.cell_data[name]
        if len(blocks) != len(source.cells):
            raise ValueError(
                f"meshio++: interpolate: source cell_data {name!r} does not have "
                "one array per cell block"
            )
        comps = {
            int(np.asarray(b).reshape(len(b), -1).shape[1]) for b in blocks if len(b)
        }
        if len(comps) > 1:
            raise ValueError(
                f"meshio++: interpolate: source cell_data {name!r} has inconsistent "
                "components across blocks"
            )

    point_out = [
        _resolve_out_name(target.point_data, n, on_conflict) for n in point_names
    ]
    cell_out = [_resolve_out_name(target.cell_data, n, on_conflict) for n in cell_names]

    out = Mesh(np.asarray(target.points), [(b.type, b.data) for b in target.cells])
    out.point_data = {k: v for k, v in target.point_data.items()}
    out.cell_data = {k: list(v) for k, v in target.cell_data.items()}
    out.field_data = {k: v for k, v in target.field_data.items()}

    tgt3 = _coords3(target.points)
    nt = tgt3.shape[0]

    # --- point_data: sampled at the target's points -------------------------
    if point_names:
        if method == "nearest":
            src3 = _coords3(source.points)
            lo, hi = src3.min(axis=0), src3.max(axis=0)
            cell = _cell_size(lo, hi, src3.shape[0])
            grid, klo, khi = _build_grid(src3, cell)
            nearest = np.empty(nt, dtype=np.int64)
            for t in range(nt):
                nearest[t] = _nearest(grid, klo, khi, cell, src3, tgt3[t])
            for name, out_name in zip(point_names, point_out):
                src = np.asarray(source.point_data[name])
                out.point_data[out_name] = src[nearest]
        else:
            from ._convert_cells import convert_cells

            smesh = convert_cells(source, mode="simplexify")
            types = {b.type for b in smesh.cells if len(b.data)}
            if "tetra" in types:
                want, nv = "tetra", 4
            elif "triangle" in types:
                want, nv = "triangle", 3
            else:
                raise ValueError(
                    "meshio++: interpolate: barycentric needs a source with triangle "
                    "or tetrahedron cells after simplexification"
                )
            sconn = _simplex_conn(smesh, want, nv)
            src3 = _coords3(smesh.points)
            lo, hi = src3.min(axis=0), src3.max(axis=0)
            cell = _cell_size(lo, hi, sconn.shape[0])
            grid = _fill_simplex_grid(src3, sconn, cell)

            simplex_of = np.full(nt, -1, dtype=np.int64)
            wts = np.zeros((nt, 4), dtype=np.float64)
            for t in range(nt):
                p = tgt3[t]
                key = (
                    int(math.floor(p[0] / cell)),
                    int(math.floor(p[1] / cell)),
                    int(math.floor(p[2] / cell)),
                )
                for s in grid.get(key, ()):
                    nodes = sconn[s]
                    corners = [src3[nodes[k]] for k in range(nv)]
                    w = (
                        _tet_weights(*corners, p)
                        if nv == 4
                        else _tri_weights(*corners, p)
                    )
                    if w is not None:
                        simplex_of[t] = s
                        wts[t, :nv] = w
                        break

            outside = np.flatnonzero(simplex_of < 0)
            near_of = None
            if outside.size and extrapolate:
                pcell = _cell_size(lo, hi, src3.shape[0])
                pgrid, pklo, pkhi = _build_grid(src3, pcell)
                near_of = np.full(nt, -1, dtype=np.int64)
                for t in outside:
                    near_of[t] = _nearest(pgrid, pklo, pkhi, pcell, src3, tgt3[t])

            for name, out_name in zip(point_names, point_out):
                src = np.asarray(smesh.point_data[name], dtype=np.float64)
                flat = src.reshape(src.shape[0], -1)
                comps = flat.shape[1]
                res = np.empty((nt, comps), dtype=np.float64)
                for t in range(nt):
                    s = simplex_of[t]
                    if s >= 0:
                        nodes = sconn[s]
                        # Corner by corner in connectivity order — the C++
                        # accumulation order.
                        acc = wts[t, 0] * flat[nodes[0]]
                        for k in range(1, nv):
                            acc = acc + wts[t, k] * flat[nodes[k]]
                        res[t] = acc
                    elif extrapolate:
                        res[t] = flat[near_of[t]]
                    else:
                        res[t] = default_value
                out.point_data[out_name] = res.reshape((nt,) + src.shape[1:])

    # --- cell_data: always nearest source-cell centroid ---------------------
    if cell_names:
        from ._partition import _centroids

        src_total = sum(len(b) for b in source.cells)
        if src_total == 0:
            raise ValueError(
                "meshio++: interpolate: cell_data was requested but the source mesh "
                "has no cells"
            )
        if len(target.cells) > 0:
            tgt_total = sum(len(b) for b in target.cells)
            scent = _centroids(source, src_total)
            tcent = _centroids(target, tgt_total)
            lo, hi = scent.min(axis=0), scent.max(axis=0)
            cell = _cell_size(lo, hi, src_total)
            grid, klo, khi = _build_grid(scent, cell)
            nearest = np.empty(tgt_total, dtype=np.int64)
            for t in range(tgt_total):
                nearest[t] = _nearest(grid, klo, khi, cell, scent, tcent[t])

            for name, out_name in zip(cell_names, cell_out):
                blocks = [np.asarray(b) for b in source.cell_data[name]]
                dtypes = {b.dtype for b in blocks}
                if len(dtypes) > 1:
                    blocks = [b.astype(np.float64) for b in blocks]
                flat = np.concatenate([b.reshape(len(b), -1) for b in blocks], axis=0)
                trailing = blocks[0].shape[1:]
                gathered = flat[nearest].reshape((tgt_total,) + trailing)
                pieces = []
                base = 0
                for b in target.cells:
                    pieces.append(gathered[base : base + len(b.data)])
                    base += len(b.data)
                out.cell_data[out_name] = pieces

    return out


def interpolate(
    source,
    target,
    method="nearest",
    arrays=None,
    extrapolate=False,
    default_value=0.0,
    on_conflict="error",
):
    """Sample data arrays from ``source`` onto ``target``.

    Parameters
    ----------
    source :
        the mesh whose data is sampled (unmodified).
    target :
        the mesh receiving the samples (unmodified); its geometry,
        connectivity, own data and sets are preserved exactly.
    method :
        ``"nearest"`` (default) or ``"barycentric"`` — see the module
        docstring. cell_data always transfers by nearest source-cell centroid.
    arrays :
        source array names to transfer, or ``None`` (default) for every source
        ``point_data`` array; ``cell_data`` transfers only when named
        explicitly. An unknown name raises.
    extrapolate :
        barycentric only — a target point outside the source domain takes the
        nearest source point's value instead of ``default_value``.
    default_value :
        barycentric only — the fill value for target points outside the source
        domain when ``extrapolate`` is off.
    on_conflict :
        ``"error"`` (default), ``"overwrite"``, or ``"suffix"`` (writes to
        ``name + "_interp"``) when a transferred name already exists on the
        target.

    Returns
    -------
    Mesh
        a copy of the target with the requested source arrays attached. The
        target's ``point_sets``/``cell_sets`` ride through untouched; source
        sets are **not** transferred; ``mesh.info``/``gmsh_periodic`` are not
        carried.
    """
    out = None
    try:
        from . import _core

        out = _core.interpolate(
            source,
            target,
            str(method),
            None if arrays is None else [str(a) for a in arrays],
            bool(extrapolate),
            float(default_value),
            str(on_conflict),
        )
    except (ValueError, TypeError):
        # A genuine user error (bad method/conflict name, unknown array, a name
        # collision under 'error') must not fall through to the numpy path —
        # it would just raise the same thing more slowly.
        raise
    except Exception:
        out = None
    if out is None:
        out = _interpolate_py(
            source, target, method, arrays, extrapolate, default_value, on_conflict
        )

    # Nothing is renumbered, so the target's sets pass through verbatim (they
    # never reach the C++ core); source sets are deliberately not transferred.
    out.point_sets = {k: v for k, v in target.point_sets.items()}
    out.cell_sets = {k: list(v) for k, v in target.cell_sets.items()}
    return out

"""Field differential operators: the gradient, divergence and curl of a field.

A dependency-free mesh *operation* (not a file format). It is what lets a
*derived* quantity be contoured -- :func:`meshioplusplus.isosurface` on
``|grad T|``, on vorticity -- and what produces the error indicator the selective
:func:`meshioplusplus.refine` consumes.

**Why this lives in the mesh-operations layer rather than the data bundle.** The
five ``data_*`` operations are defined by never touching geometry; this one
consumes and produces data arrays but *reads* geometry and topology (face areas,
cell volumes, cell adjacency). It is still reachable as
``meshioplusplus data gradient`` in both CLIs, because that is where a user looks
for it.

The two methods
---------------

**Green-Gauss** (the default) applies the divergence theorem over the cell:
``grad f = (1/V) sum_faces integral(f n dA)``. Each face is fanned into triangles
about the arithmetic mean of its corners, and each sub-triangle contributes
``A_j * f(centroid_j)``. That is **exact for a linear field on any 3-D cell**,
planar faces or not: two faces sharing an edge contribute oppositely-wound
triangles there, so the fan surface is closed and ``∮ f n dA = V grad f`` applies
verbatim. The corner-average fan apex is forced rather than convenient -- for
linear ``f`` it is the only apex whose value is known exactly without shape
functions. On a 2-D cell the same theorem runs over the corner ring with the
in-plane outward normal, exact on a **planar** cell.

**Least-squares** fits ``f(x) ~ f_c + g . (x - x_c)`` over the cells sharing at
least one node, with ``x_c``/``f_c`` the arithmetic means of the cell's corners,
so ``(x_c, f_c)`` lies exactly on a linear field and the fit is exact for one.
Weights are ``1/|d|^2``. A degenerate neighbourhood (an isolated cell, a
collinear strip) falls back to Green-Gauss for that cell and is counted.

Contracts worth knowing (all documented at length in ``doc/gradient.md``):

* **The field must be ``point_data``.** ``cell_data`` is piecewise constant and
  has no derivative; naming one raises, pointing at
  :func:`meshioplusplus.cell_data_to_point_data` (CLI ``data to-point``).
* **Shapes.** An ``nc``-component input yields ``3 * nc`` gradient components,
  flat and row-major as ``[component i][derivative j]`` at ``i * 3 + j``: a
  scalar gives ``(n, 3)``, a 3-vector gives ``(n, 9)``. Divergence gives 1
  component, curl always 3.
* **Unsupported cells are NaN and counted**, never approximated.
* Output is always ``float64``, named ``<input>:gradient`` / ``:divergence`` /
  ``:curl`` by default.

The C++ core (``_core.gradient``) does the work; this module is the thin shim
(try C++, fall back to a pure-numpy reference) plus the reference itself, which
reproduces the C++ result **byte for byte**, pinned by
``tests/python/test_gradient.py::test_cpp_matches_python``.

Public API:

* :func:`gradient` -- the gradient, divergence or curl of a ``point_data`` field.
"""

from __future__ import annotations

import numpy as np

from ._convert_cells import _LINEAR_BASE, _NUM_CORNERS
from ._data_average import cell_data_to_point_data
from ._data_manage import data_drop
from ._mesh import topological_dimension
from ._skin import _CELL_FACES

#: Default output-name suffix per operator.
_SUFFIX = {"gradient": ":gradient", "divergence": ":divergence", "curl": ":curl"}

_OPERATOR_ALIASES = {
    "": "gradient",
    "gradient": "gradient",
    "grad": "gradient",
    "divergence": "divergence",
    "div": "divergence",
    "curl": "curl",
    "rot": "curl",
}

_METHOD_ALIASES = {
    "": "green-gauss",
    "green-gauss": "green-gauss",
    "green_gauss": "green-gauss",
    "gg": "green-gauss",
    "least-squares": "least-squares",
    "least_squares": "least-squares",
    "lsq": "least-squares",
}

_PREFIX = "meshio++: gradient: "


def _corner_count(cell_type):
    """Corner count of ``cell_type``, or 0 for the variable-node-count types.

    The twin of the C++ ``detail::cell_corner_count``: corners always lead the
    meshio/VTK ordering, so a higher-order cell reduces to its linear parent by
    keeping the first ``_corner_count`` columns. ``_LINEAR_BASE``/``_NUM_CORNERS``
    (from :mod:`._convert_cells`) cover the serendipity family directly; the
    Lagrange-numbered types fall through to the family-name rule, which is what
    the C++ table spells out case by case.
    """
    base = _LINEAR_BASE.get(cell_type, cell_type)
    if base in _NUM_CORNERS:
        return _NUM_CORNERS[base]
    for family, n in (
        ("hexahedron", 8),
        ("triangle", 3),
        ("tetra", 4),
        ("wedge", 6),
        ("pyramid", 5),
        ("quad", 4),
        ("line", 2),
        ("vertex", 1),
    ):
        # `polygon`/`polyhedron`/`VTK_LAGRANGE_*`/`custom` never match and stay 0,
        # exactly as the C++ default case returns 0 for them.
        if cell_type.startswith(family) and cell_type[len(family) :].isdigit():
            return n
    return 0


def _mesh_dimension(mesh):
    """Max topological dimension over the mesh's non-empty blocks."""
    dim = -1
    for cb in mesh.cells:
        data = cb.data
        if isinstance(data, list):
            if not data:
                continue
            # A ragged block: a polygon row list is 2-D, a polyhedron 3-D.
            d = (
                3
                if isinstance(data[0], list)
                and data[0]
                and isinstance(data[0][0], (list, tuple, np.ndarray))
                else 2
            )
        else:
            if len(data) == 0:
                continue
            d = topological_dimension.get(cb.type, -1)
        dim = max(dim, d)
    return dim


def _shape(rows, ncomp):
    """``(rows,)`` for a scalar, ``(rows, ncomp)`` otherwise."""
    return (rows,) if ncomp <= 1 else (rows, ncomp)


def _flat(arr):
    """``(rows, ncomp)`` view of a data array, whatever its trailing dims."""
    a = np.asarray(arr)
    rows = a.shape[0] if a.ndim >= 1 else 0
    ncomp = int(a.size // rows) if rows else 0
    return a.reshape(rows, ncomp) if rows else a.reshape(0, 0), ncomp


# --------------------------------------------------------------------------- #
# Green-Gauss                                                                  #
# --------------------------------------------------------------------------- #
#
# Both paths are vectorized over CELLS while looping in Python over the fixed
# face / fan-triangle / ring-edge sequence. That is what makes the accumulation
# order per cell identical to the C++ loop nest: each cell adds the same terms in
# the same order, and numpy only parallelizes across cells, which never mixes
# terms. `np.sum`/`np.mean`/`np.einsum` are BANNED throughout -- numpy switches to
# pairwise summation and guarantees no left-to-right order.


def _mean_axis0(rows):
    """Left-to-right mean of a list of ``(ncells, k)`` arrays.

    Deliberately not ``np.mean``/``np.sum`` over a stacked axis: those reduce
    pairwise and would round differently from the C++ sequential loop.
    """
    acc = rows[0]
    for r in rows[1:]:
        acc = acc + r
    return acc * (1.0 / len(rows))


def _cross(a, b):
    """Row-wise cross product, written out so the operand order matches C++."""
    return np.stack(
        (
            a[:, 1] * b[:, 2] - a[:, 2] * b[:, 1],
            a[:, 2] * b[:, 0] - a[:, 0] * b[:, 2],
            a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0],
        ),
        axis=1,
    )


def _dot(a, b):
    return a[:, 0] * b[:, 0] + a[:, 1] * b[:, 1] + a[:, 2] * b[:, 2]


def _norm(a):
    return np.sqrt(_dot(a, a))


def _green_gauss_3d(corners, values, faces, ncomp):
    """Green-Gauss over the face fan. ``corners`` is ``(ncells, k, 3)`` and
    ``values`` ``(ncells, k, ncomp)``, both already recentred."""
    n = corners.shape[0]
    num = np.zeros((n, ncomp, 3), dtype=np.float64)
    volume = np.zeros(n, dtype=np.float64)
    area_scale = np.zeros(n, dtype=np.float64)

    for _ftype, m, nodes in faces:
        ring = nodes[:m]
        c = _mean_axis0([corners[:, j, :] for j in ring])
        fc = _mean_axis0([values[:, j, :] for j in ring])
        for i in range(m):
            a = ring[i]
            b = ring[(i + 1) % m]
            pa = corners[:, a, :]
            pb = corners[:, b, :]
            av = _cross(pa - c, pb - c) * 0.5
            xj = (c + pa + pb) * (1.0 / 3.0)
            volume = volume + _dot(xj, av) * (1.0 / 3.0)
            area_scale = area_scale + _norm(av)
            fj = (fc + values[:, a, :] + values[:, b, :]) * (1.0 / 3.0)
            for k in range(ncomp):
                num[:, k, 0] = num[:, k, 0] + fj[:, k] * av[:, 0]
                num[:, k, 1] = num[:, k, 1] + fj[:, k] * av[:, 1]
                num[:, k, 2] = num[:, k, 2] + fj[:, k] * av[:, 2]

    ok = np.abs(volume) > 1e-12 * (area_scale * np.sqrt(area_scale))
    inv_v = np.where(ok, 1.0 / np.where(ok, volume, 1.0), 0.0)
    return num * inv_v[:, None, None], ok


def _ring_area_vector(corners):
    """Recentred Newell area vector of a ``(ncells, k, 3)`` corner ring."""
    n = corners.shape[1]
    av = np.zeros((corners.shape[0], 3), dtype=np.float64)
    p0 = corners[:, 0, :]
    for i in range(1, n - 1):
        av = av + _cross(corners[:, i, :] - p0, corners[:, i + 1, :] - p0)
    return av * 0.5


def _green_gauss_2d(corners, values, ncomp):
    """Green-Gauss over a 2-D cell's corner ring."""
    n = corners.shape[0]
    k = corners.shape[1]
    av = _ring_area_vector(corners)
    area = _norm(av)
    edge_scale = np.zeros(n, dtype=np.float64)
    for i in range(k):
        edge_scale = edge_scale + _norm(corners[:, (i + 1) % k, :] - corners[:, i, :])
    ok = area > 1e-12 * (edge_scale * edge_scale)
    safe = np.where(ok, area, 1.0)
    nrm = av / safe[:, None]

    num = np.zeros((n, ncomp, 3), dtype=np.float64)
    for i in range(k):
        a = i
        b = (i + 1) % k
        t = corners[:, b, :] - corners[:, a, :]
        mds = _cross(t, nrm)
        fbar = (values[:, a, :] + values[:, b, :]) * 0.5
        for c in range(ncomp):
            num[:, c, 0] = num[:, c, 0] + fbar[:, c] * mds[:, 0]
            num[:, c, 1] = num[:, c, 1] + fbar[:, c] * mds[:, 1]
            num[:, c, 2] = num[:, c, 2] + fbar[:, c] * mds[:, 2]

    inv_a = np.where(ok, 1.0 / safe, 0.0)
    return num * inv_a[:, None, None], ok


def _green_gauss_1d(corners, values, ncomp):
    """Green-Gauss on a line cell: the two endpoints are its faces."""
    t = corners[:, 1, :] - corners[:, 0, :]
    length = _norm(t)
    ok = length > 0.0
    safe = np.where(ok, length, 1.0)
    direction = t / safe[:, None]
    num = np.zeros((corners.shape[0], ncomp, 3), dtype=np.float64)
    for k in range(ncomp):
        slope = (values[:, 1, k] - values[:, 0, k]) / safe
        num[:, k, 0] = slope * direction[:, 0]
        num[:, k, 1] = slope * direction[:, 1]
        num[:, k, 2] = slope * direction[:, 2]
    return num, ok


def _green_gauss(corners, values, cell_type, dim, ncomp):
    if dim == 3:
        faces = _CELL_FACES.get(cell_type)
        if faces is None:
            return None, np.zeros(corners.shape[0], dtype=bool)
        return _green_gauss_3d(corners, values, faces, ncomp)
    if dim == 2:
        return _green_gauss_2d(corners, values, ncomp)
    if dim == 1:
        return _green_gauss_1d(corners, values, ncomp)
    return None, np.zeros(corners.shape[0], dtype=bool)


# --------------------------------------------------------------------------- #
# Least squares                                                                #
# --------------------------------------------------------------------------- #


def _solve3(m, b):
    """Row-wise symmetric 3x3 cofactor solve of ``M g = b``.

    Transcribed from the C++ ``grad_solve3``, which is itself
    ``decimate.cpp``'s optimal-placement solve -- including the
    ``scale * scale * scale`` multiplication form, since ``np.power`` rounds
    differently and would break parity. ``np.linalg.solve`` is banned here for
    the same reason.
    """
    c00 = m[:, 4] * m[:, 8] - m[:, 5] * m[:, 5]
    c01 = m[:, 2] * m[:, 5] - m[:, 1] * m[:, 8]
    c02 = m[:, 1] * m[:, 5] - m[:, 2] * m[:, 4]
    det = m[:, 0] * c00 + m[:, 1] * c01 + m[:, 2] * c02
    scale = np.abs(m[:, 0])
    scale = np.maximum(scale, np.abs(m[:, 1]))
    scale = np.maximum(scale, np.abs(m[:, 2]))
    scale = np.maximum(scale, np.abs(m[:, 4]))
    scale = np.maximum(scale, np.abs(m[:, 5]))
    scale = np.maximum(scale, np.abs(m[:, 8]))
    ok = np.abs(det) > 1e-12 * (scale * scale * scale)
    c11 = m[:, 0] * m[:, 8] - m[:, 2] * m[:, 2]
    c12 = m[:, 1] * m[:, 2] - m[:, 0] * m[:, 5]
    c22 = m[:, 0] * m[:, 4] - m[:, 1] * m[:, 1]
    inv = np.where(ok, 1.0 / np.where(ok, det, 1.0), 0.0)
    out = np.zeros((m.shape[0], 3), dtype=np.float64)
    out[:, 0] = (c00 * b[:, 0] + c01 * b[:, 1] + c02 * b[:, 2]) * inv
    out[:, 1] = (c01 * b[:, 0] + c11 * b[:, 1] + c12 * b[:, 2]) * inv
    out[:, 2] = (c02 * b[:, 0] + c12 * b[:, 1] + c22 * b[:, 2]) * inv
    return out, ok


def _cell_node_lists(mesh):
    """Global (block-major) cell -> its in-range node ids, as a list of arrays.

    The twin of ``detail::build_cell_node_incidence``.
    """
    npoints = len(mesh.points)
    rows = []
    for cb in mesh.cells:
        data = cb.data
        if isinstance(data, list):
            for row in data:
                if row and isinstance(row[0], (list, tuple, np.ndarray)):
                    ids = np.concatenate([np.asarray(f, dtype=np.int64) for f in row])
                else:
                    ids = np.asarray(row, dtype=np.int64)
                rows.append(ids[(ids >= 0) & (ids < npoints)])
        else:
            arr = np.asarray(data, dtype=np.int64)
            for r in arr:
                rows.append(r[(r >= 0) & (r < npoints)])
    return rows


def _neighbor_lists(mesh):
    """Per global cell, the ascending de-duplicated node-sharing neighbours.

    The twin of ``detail::cell_node_neighbors``. The ascending order is part of
    the contract, not an implementation detail: the least-squares sum is
    accumulated in it.
    """
    cell_nodes = _cell_node_lists(mesh)
    npoints = len(mesh.points)
    node_cells = [[] for _ in range(npoints)]
    for c, ids in enumerate(cell_nodes):
        for nid in ids:
            node_cells[int(nid)].append(c)
    out = []
    for c, ids in enumerate(cell_nodes):
        acc = []
        for nid in ids:
            acc.extend(node_cells[int(nid)])
        if acc:
            uniq = np.unique(np.asarray(acc, dtype=np.int64))
            out.append(uniq[uniq != c])
        else:
            out.append(np.zeros(0, dtype=np.int64))
    return out


def _least_squares(offsets, deltas, degree, normals, ncomp, dim):
    """Row-wise least-squares fit.

    ``offsets`` is ``(ncells, maxdeg, 3)`` and ``deltas``
    ``(ncells, maxdeg, ncomp)``, both zero-padded past each cell's ``degree``.
    The accumulation is a **padded-gather loop** over the neighbour slot index,
    so every cell sums its neighbours in ascending order exactly as C++ does --
    ``np.add.reduceat`` would not guarantee that.
    """
    n = offsets.shape[0]
    maxdeg = offsets.shape[1]
    m = np.zeros((n, 9), dtype=np.float64)
    b = np.zeros((n, ncomp, 3), dtype=np.float64)

    for j in range(maxdeg):
        d = offsets[:, j, :]
        if normals is not None:
            d = d - normals * _dot(d, normals)[:, None]
        d2 = _dot(d, d)
        # `M = M + term` under a mask, never `M += where(mask, term, 0)`:
        # -0.0 + 0.0 == +0.0 is a real (if rare) parity break.
        live = (j < degree) & (d2 > 0.0)
        w = np.where(live, 1.0 / np.where(d2 > 0.0, d2, 1.0), 0.0)
        for r in range(3):
            for c in range(3):
                term = w * d[:, r] * d[:, c]
                m[:, r * 3 + c] = np.where(
                    live, m[:, r * 3 + c] + term, m[:, r * 3 + c]
                )
        for k in range(ncomp):
            wf = w * deltas[:, j, k]
            for c in range(3):
                term = wf * d[:, c]
                b[:, k, c] = np.where(live, b[:, k, c] + term, b[:, k, c])

    if normals is not None:
        # Rank-1 regularization instead of a tangent-plane basis: it makes the
        # in-plane rank-2 matrix invertible without touching b (every offset was
        # projected, so b . N == 0), so the SAME 3x3 solver serves both
        # dimensions and there is no arbitrary basis to transcribe.
        s = (m[:, 0] + m[:, 4] + m[:, 8]) * 0.5
        s = np.where(s > 0.0, s, 1.0)
        for r in range(3):
            for c in range(3):
                m[:, r * 3 + c] = m[:, r * 3 + c] + s * normals[:, r] * normals[:, c]

    out = np.zeros((n, ncomp, 3), dtype=np.float64)
    if dim in (2, 3):
        ok = np.ones(n, dtype=bool)
        for k in range(ncomp):
            g, good = _solve3(m, b[:, k, :])
            out[:, k, :] = g
            ok = ok & good
        return out, ok

    trace = m[:, 0] + m[:, 4] + m[:, 8]
    ok = trace != 0.0
    safe = np.where(ok, trace, 1.0)
    for k in range(ncomp):
        out[:, k, :] = b[:, k, :] / safe[:, None]
    return out, ok


# --------------------------------------------------------------------------- #
# Operator reduction                                                           #
# --------------------------------------------------------------------------- #


def _validate(mesh, array, operator, method, location, output, component, overwrite):
    """Raise on every user error, before either path runs.

    Shared deliberately: the C++ core validates the same things, but the numpy
    fallback must reject the same inputs with the same message even when
    ``_core`` is unavailable -- otherwise the error surface would depend on
    whether the extension happens to be built.
    """
    if operator not in _OPERATOR_ALIASES:
        raise ValueError(
            f"{_PREFIX}unknown operator '{operator}' "
            "(expected 'gradient', 'divergence', or 'curl')"
        )
    if method not in _METHOD_ALIASES:
        raise ValueError(
            f"{_PREFIX}unknown method '{method}' "
            "(expected 'green-gauss' or 'least-squares')"
        )
    if location not in ("cell", "point"):
        raise ValueError(
            f"{_PREFIX}unknown location '{location}' (expected 'cell' or 'point')"
        )
    op = _OPERATOR_ALIASES[operator]

    if not array:
        raise ValueError(f"{_PREFIX}an array name is required")
    if array not in mesh.point_data:
        # A cell_data field is piecewise constant and has no derivative. Name the
        # fix, exactly as isosurface does for the same reason.
        if array in mesh.cell_data:
            raise ValueError(
                f"{_PREFIX}'{array}' is cell_data, which is piecewise constant and "
                "has no derivative; convert it first with cell_data_to_point_data "
                "(CLI: meshioplusplus data to-point)"
            )
        have = ", ".join(sorted(mesh.point_data))
        raise ValueError(
            f"{_PREFIX}no point_data array named '{array}'"
            + (f" (available: {have})" if have else " (the mesh has no point_data)")
        )

    flat, ncomp = _flat(mesh.point_data[array])
    if ncomp == 0 or flat.shape[0] == 0:
        raise ValueError(f"{_PREFIX}'{array}' carries no values")
    if flat.shape[0] != len(mesh.points):
        raise ValueError(
            f"{_PREFIX}'{array}' has {flat.shape[0]} rows on a mesh of "
            f"{len(mesh.points)} points"
        )
    if component is not None:
        if op != "gradient":
            raise ValueError(
                f"{_PREFIX}a component selection applies to the gradient only; "
                "divergence and curl consume the whole vector"
            )
        if int(component) < 0 or int(component) >= ncomp:
            raise ValueError(
                f"{_PREFIX}component {int(component)} is out of range for "
                f"'{array}', which has {ncomp} component(s)"
            )
    if op != "gradient" and ncomp not in (2, 3):
        raise ValueError(
            f"{_PREFIX}'{array}' has {ncomp} components; divergence and curl need "
            "a vector field of 2 or 3"
        )

    out_name = output or (array + _SUFFIX[op])
    if not overwrite:
        taken = (
            out_name in mesh.point_data
            if location == "point"
            else out_name in mesh.cell_data
        )
        if taken:
            raise ValueError(
                f"{_PREFIX}'{out_name}' already exists in {location}_data "
                "(pass overwrite=true to replace it)"
            )


def _out_components(operator, ncomp):
    if operator == "divergence":
        return 1
    if operator == "curl":
        return 3
    return ncomp * 3


def _apply_operator(operator, grad, ncomp):
    """Reduce a ``(ncells, ncomp, 3)`` gradient block to the requested operator."""
    n = grad.shape[0]
    if operator == "gradient":
        return grad.reshape(n, ncomp * 3)
    g = np.zeros((n, 3, 3), dtype=np.float64)
    # A 2-component field reads as (u, v, 0), the same padding convention that
    # applies to 2-D point coordinates.
    for k in range(min(ncomp, 3)):
        g[:, k, :] = grad[:, k, :]
    if operator == "divergence":
        return (g[:, 0, 0] + g[:, 1, 1] + g[:, 2, 2]).reshape(n, 1)
    out = np.zeros((n, 3), dtype=np.float64)
    out[:, 0] = g[:, 2, 1] - g[:, 1, 2]
    out[:, 1] = g[:, 0, 2] - g[:, 2, 0]
    out[:, 2] = g[:, 1, 0] - g[:, 0, 1]
    return out


# --------------------------------------------------------------------------- #
# The numpy reference                                                          #
# --------------------------------------------------------------------------- #


def _gradient_py(
    mesh,
    array,
    operator="gradient",
    method="green-gauss",
    location="cell",
    output=None,
    component=None,
    overwrite=False,
):
    """Pure-numpy reference implementation of :func:`gradient`.

    Returns ``(mesh, num_skipped, num_fallback)``. Reproduces ``_core.gradient``
    byte for byte; see the module docstring for the parity rules.
    """
    operator = _OPERATOR_ALIASES[operator]
    method = _METHOD_ALIASES[method]
    out_name = output or (array + _SUFFIX[operator])

    field, field_comp = _flat(mesh.point_data[array])
    work_comp = 1 if component is not None else field_comp
    base = int(component) if component is not None else 0
    work = field[:, base : base + work_comp].astype(np.float64)
    out_comp = _out_components(operator, work_comp)

    dim = _mesh_dimension(mesh)
    points = np.asarray(mesh.points, dtype=np.float64)
    if points.shape[1] < 3:
        pad = np.zeros((points.shape[0], 3 - points.shape[1]), dtype=np.float64)
        points = np.concatenate([points, pad], axis=1)

    # Block eligibility, resolved once — the twin of the C++ pass.
    eligible = []
    corners_per_block = []
    for cb in mesh.cells:
        nc = 0 if isinstance(cb.data, list) else _corner_count(cb.type)
        ok = (
            not isinstance(cb.data, list)
            and topological_dimension.get(cb.type, -1) == dim
            and nc > 0
            and (dim != 3 or cb.type in _CELL_FACES)
            and (dim != 2 or nc >= 3)
            and (dim != 1 or nc >= 2)
        )
        eligible.append(ok)
        corners_per_block.append(nc)

    # Per-block corner coordinates/values and their recentring origins.
    loaded = []
    for b, cb in enumerate(mesh.cells):
        if not eligible[b]:
            loaded.append(None)
            continue
        conn = np.asarray(cb.data, dtype=np.int64)[:, : corners_per_block[b]]
        pc = points[conn]  # (ncells, k, 3)
        vc = work[conn]  # (ncells, k, work_comp)
        origin = _mean_axis0([pc[:, j, :] for j in range(pc.shape[1])])
        value0 = _mean_axis0([vc[:, j, :] for j in range(vc.shape[1])])
        loaded.append(
            (pc - origin[:, None, :], vc - value0[:, None, :], origin, value0)
        )

    neighbors = None
    cell_x = cell_f = cell_ok = None
    bases = np.cumsum([0] + [len(cb.data) for cb in mesh.cells])
    if method == "least-squares":
        neighbors = _neighbor_lists(mesh)
        total = int(bases[-1])
        cell_x = np.zeros((total, 3), dtype=np.float64)
        cell_f = np.zeros((total, work_comp), dtype=np.float64)
        cell_ok = np.zeros(total, dtype=bool)
        for b, cb in enumerate(mesh.cells):
            if loaded[b] is None:
                continue
            _pc, _vc, origin, value0 = loaded[b]
            lo, hi = int(bases[b]), int(bases[b + 1])
            cell_x[lo:hi] = origin
            cell_f[lo:hi] = value0
            cell_ok[lo:hi] = True

    blocks = []
    skipped = 0
    fallback = 0
    for b, cb in enumerate(mesh.cells):
        ncells = len(cb.data)
        if not eligible[b]:
            blocks.append(np.full(_shape(ncells, out_comp), np.nan, dtype=np.float64))
            skipped += ncells
            continue

        pc, vc, origin, _value0 = loaded[b]
        gg, gg_ok = _green_gauss(pc, vc, cb.type, dim, work_comp)
        if gg is None:
            gg = np.zeros((ncells, work_comp, 3), dtype=np.float64)

        if method == "least-squares":
            maxdeg = 0
            for c in range(ncells):
                nb = neighbors[int(bases[b]) + c]
                maxdeg = max(maxdeg, int(np.count_nonzero(cell_ok[nb])))
            offsets = np.zeros((ncells, max(maxdeg, 1), 3), dtype=np.float64)
            deltas = np.zeros((ncells, max(maxdeg, 1), work_comp), dtype=np.float64)
            degree = np.zeros(ncells, dtype=np.int64)
            for c in range(ncells):
                nb = neighbors[int(bases[b]) + c]
                nb = nb[cell_ok[nb]]
                degree[c] = len(nb)
                if len(nb):
                    offsets[c, : len(nb), :] = cell_x[nb] - origin[c]
                    deltas[c, : len(nb), :] = cell_f[nb] - _value0_row(loaded[b], c)
            normals = None
            if dim == 2:
                av = _ring_area_vector(pc)
                a = _norm(av)
                safe_a = np.where(a > 0.0, a, 1.0)
                normals = np.where((a > 0.0)[:, None], av / safe_a[:, None], 0.0)
            grad, ok = _least_squares(offsets, deltas, degree, normals, work_comp, dim)
            fell = (~ok) & gg_ok
            grad = np.where(fell[:, None, None], gg, grad)
            ok = ok | fell
            fallback += int(np.count_nonzero(fell))
        else:
            grad, ok = gg, gg_ok

        values = _apply_operator(operator, grad, work_comp)
        values = np.where(ok[:, None], values, np.nan)
        skipped += int(np.count_nonzero(~ok))
        blocks.append(values.reshape(_shape(ncells, out_comp)))

    out = mesh.copy()
    out.cell_data[out_name] = blocks

    if location == "point":
        # Compose the existing averaging rather than reimplementing a scatter --
        # its accumulation is already deliberately serial, so the point values
        # stay reproducible.
        out = cell_data_to_point_data(out, keys=[out_name])
        out = data_drop(out, "cell", [out_name])

    return out, skipped, fallback


def _value0_row(loaded_block, cell):
    """The recentring value offset of one cell of a loaded block."""
    return loaded_block[3][cell]


# --------------------------------------------------------------------------- #
# Public API                                                                   #
# --------------------------------------------------------------------------- #


def gradient(
    mesh,
    array,
    operator="gradient",
    method="green-gauss",
    location="cell",
    output=None,
    component=None,
    overwrite=False,
    return_report=False,
):
    """Differentiate a ``point_data`` field over the mesh.

    :param mesh: the mesh to read (never modified).
    :param array: name of the ``point_data`` array to differentiate.
    :param operator: ``"gradient"`` (default), ``"divergence"`` or ``"curl"``.
    :param method: ``"green-gauss"`` (default) or ``"least-squares"``.
    :param location: ``"cell"`` (default) or ``"point"`` for the result.
    :param output: output array name; ``None`` uses ``<array>:<operator>``.
    :param component: differentiate only this component of a multi-component
        input (gradient only). ``None`` means every component -- note this is the
        opposite of :func:`meshioplusplus.isosurface`, where ``None`` means the
        row magnitude.
    :param overwrite: allow replacing an existing array of the output name.
    :param return_report: also return the skip/fallback counters.
    :returns: the mesh carrying the produced array, or
        ``(mesh, report)`` when ``return_report`` is set.
    :raises ValueError: on an unknown or ``cell_data`` array name, a bad
        component, a component count divergence/curl cannot use, or an output
        name collision without ``overwrite``.
    """
    _validate(mesh, array, operator, method, location, output, component, overwrite)

    out = None
    report = None
    try:
        from . import _core

        res = _core.gradient(
            mesh,
            array,
            operator,
            method,
            location,
            "" if output is None else str(output),
            -1 if component is None else int(component),
            bool(overwrite),
        )
        out = res["mesh"]
        report = {
            "num_skipped": res["num_skipped"],
            "num_fallback": res["num_fallback"],
        }
    except (ValueError, TypeError):
        # A genuine user error (an unknown array, a cell_data name, a bad
        # component) must not fall through to the numpy path, which would either
        # raise something less helpful or silently succeed differently.
        raise
    except Exception:
        out = None

    if out is None:
        out, skipped, fell = _gradient_py(
            mesh,
            array,
            operator=operator,
            method=method,
            location=location,
            output=output,
            component=component,
            overwrite=overwrite,
        )
        report = {"num_skipped": skipped, "num_fallback": fell}

    # Sets are carried through untouched: nothing is renumbered here.
    if getattr(mesh, "point_sets", None):
        out.point_sets = {k: np.asarray(v).copy() for k, v in mesh.point_sets.items()}
    if getattr(mesh, "cell_sets", None):
        out.cell_sets = {
            k: [np.asarray(a).copy() for a in v] for k, v in mesh.cell_sets.items()
        }

    return (out, report) if return_report else out

"""A size field from a per-cell error estimate (roadmap 6.10.2).

:func:`meshioplusplus.estimate_error` measures where a field is badly resolved
(``error:zz``, one value per cell); :func:`meshioplusplus.mesh_cad_surface`'s
``background`` takes a size per point. :func:`error_size` joins them by the
adaptive procedure of O. C. Zienkiewicz and J. Z. Zhu, "A simple error estimator
and adaptive procedure for practical engineering analysis", Int. J. Numer.
Meth. Engng. 24 (1987): with ``N`` cells evaluated and a target global error
``eta``, each cell may carry ``e = eta / sqrt(N)`` when the error is spread
evenly, and a cell of diameter ``h_K`` whose indicator is ``eta_K`` asks for
``h_K * (e / eta_K) ** (1 / p)``, ``p`` the indicator's rate in the cell size
(1 for linear elements). Sizes are clamped to ``[min_size, max_size]`` (0: no
bound) and written per point as the least of its cells'; a point with no
evaluated cell reads NaN, which ``background`` treats as no size there.

The C++ core (``_core.error_size``) does the work; this module is the thin shim
plus a numpy twin that matches it to rounding.

Public API:

* :func:`error_size` -- the size field.
"""

from __future__ import annotations

import math

import numpy as np

from ._fallback import core_op_declined

_PREFIX = "meshio++: error_size: "


def _validate(
    mesh, error_array, target_error, order, min_size, max_size, output, overwrite
):
    if not (math.isfinite(target_error) and target_error > 0.0):
        raise ValueError(f"{_PREFIX}target_error must be finite and positive")
    if not (math.isfinite(order) and order > 0.0):
        raise ValueError(f"{_PREFIX}order must be finite and positive")
    for v in (min_size, max_size):
        if not (math.isfinite(v) and v >= 0.0):
            raise ValueError(
                f"{_PREFIX}min_size and max_size must be finite and not negative"
            )
    if min_size > 0.0 and max_size > 0.0 and min_size > max_size:
        raise ValueError(f"{_PREFIX}min_size must not exceed max_size")
    if not output:
        raise ValueError(f"{_PREFIX}an output name is required")
    if error_array not in mesh.cell_data:
        if error_array in mesh.point_data:
            raise ValueError(
                f"{_PREFIX}'{error_array}' is point_data; the indicator is per cell "
                "(estimate_error writes error:zz)"
            )
        have = ", ".join(sorted(mesh.cell_data))
        raise ValueError(
            f"{_PREFIX}no cell_data array named '{error_array}'"
            + (f" (available: {have})" if have else " (the mesh has no cell_data)")
        )
    for block in mesh.cell_data[error_array]:
        a = np.asarray(block)
        if a.ndim > 1 and int(np.prod(a.shape[1:])) != 1:
            raise ValueError(f"{_PREFIX}'{error_array}' must have one component")
    if not overwrite and output in mesh.point_data:
        raise ValueError(
            f"{_PREFIX}'{output}' already exists in point_data "
            "(pass overwrite=True to replace it)"
        )


def _cell_rows(block):
    """Each cell's node ids, as a list of integer arrays."""
    data = block.data
    if isinstance(data, np.ndarray) and data.ndim == 2:
        return [data[i] for i in range(len(data))]
    rows = []
    for c in data:
        if len(c) and not np.isscalar(c[0]) and np.ndim(c[0]) > 0:
            # A polyhedron: its faces' ids.
            rows.append(np.unique(np.concatenate([np.asarray(f) for f in c])))
        else:
            rows.append(np.asarray(c, dtype=np.int64))
    return rows


def _error_size_py(mesh, error_array, target_error, order, min_size, max_size, output):
    pts = np.asarray(mesh.points, dtype=np.float64)
    np_ = len(pts)
    diam = []
    eta = []
    rows_all = []
    for b, block in enumerate(mesh.cells):
        e = np.asarray(mesh.cell_data[error_array][b], dtype=np.float64).reshape(-1)
        rows = _cell_rows(block)
        for c, ids in enumerate(rows):
            x = pts[np.asarray(ids, dtype=np.int64)]
            d2 = 0.0
            for i in range(len(x)):
                diff = x[i + 1 :] - x[i]
                if len(diff):
                    d2 = max(d2, float(np.max(np.sum(diff * diff, axis=1))))
            diam.append(math.sqrt(d2))
            eta.append(float(e[c]))
            rows_all.append(ids)
    diam = np.asarray(diam)
    eta = np.asarray(eta)
    usable = np.isfinite(eta) & (eta >= 0.0) & (diam > 0.0)
    num_cells = int(usable.sum())
    share = target_error / math.sqrt(num_cells) if num_cells else 0.0
    size = np.full(len(eta), np.nan)
    clamped = 0
    skipped = int((~usable).sum())
    for k in np.nonzero(usable)[0]:
        h = diam[k] * (share / eta[k]) ** (1.0 / order) if eta[k] > 0.0 else math.inf
        if min_size > 0.0 and h < min_size:
            h = min_size
            clamped += 1
        elif max_size > 0.0 and h > max_size:
            h = max_size
            clamped += 1
        if not math.isfinite(h):
            num_cells -= 1
            skipped += 1
            continue
        size[k] = h
    at = np.full(np_, np.inf)
    for k, ids in enumerate(rows_all):
        if math.isfinite(size[k]):
            ids = np.asarray(ids, dtype=np.int64)
            at[ids] = np.minimum(at[ids], size[k])
    finite = np.isfinite(at)
    values = np.where(finite, at, np.nan)
    out = mesh.copy()
    out.point_data = dict(mesh.point_data)
    out.point_data[output] = values
    report = {
        "num_cells": num_cells,
        "num_skipped": skipped,
        "num_clamped": clamped,
        "min_size": float(values[finite].min()) if finite.any() else 0.0,
        "max_size": float(values[finite].max()) if finite.any() else 0.0,
    }
    return out, report


def error_size(
    mesh,
    error_array="error:zz",
    *,
    target_error,
    order=1.0,
    min_size=0.0,
    max_size=0.0,
    output="size",
    overwrite=False,
    return_report=False,
):
    """A per-point size field from a per-cell error indicator.

    :param mesh: the source mesh (never modified), carrying ``error_array``.
    :param error_array: the one-component ``cell_data`` indicator
        (:func:`estimate_error` writes ``"error:zz"``).
    :param target_error: the target global error, in the indicator's norm.
    :param order: the indicator's rate ``p`` in the cell size (1 for linear
        elements).
    :param min_size: the least size written (0: no bound).
    :param max_size: the greatest size written (0: no bound).
    :param output: the ``point_data`` array written (``"size"``, what
        ``mesh_cad_surface``'s ``background`` reads by default).
    :param overwrite: allow replacing an existing array of that name.
    :param return_report: also return ``{num_cells, num_skipped, num_clamped,
        min_size, max_size}``.
    :returns: the mesh with the size array, or ``(mesh, report)``.
    :raises ValueError: on a missing or multi-component indicator, a target
        error or order that is not finite and positive, a negative bound,
        ``min_size > max_size``, or an output-name collision.
    """
    target_error = float(target_error)
    order = float(order)
    min_size = float(min_size)
    max_size = float(max_size)
    output = str(output)
    _validate(
        mesh, error_array, target_error, order, min_size, max_size, output, overwrite
    )
    out = None
    report = None
    try:
        from . import _core

        res = _core.error_size(
            mesh,
            error_array=str(error_array),
            target_error=target_error,
            order=order,
            min_size=min_size,
            max_size=max_size,
            output=output,
            overwrite=bool(overwrite),
        )
        out = res["mesh"]
        report = {
            "num_cells": int(res["num_cells"]),
            "num_skipped": int(res["num_skipped"]),
            "num_clamped": int(res["num_clamped"]),
            "min_size": float(res["min_size"]),
            "max_size": float(res["max_size"]),
        }
    except Exception as exc:
        if not core_op_declined(exc, "error_size"):
            raise
        out = None
    if out is None:
        out, report = _error_size_py(
            mesh, error_array, target_error, order, min_size, max_size, output
        )
    return (out, report) if return_report else out

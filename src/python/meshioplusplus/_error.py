"""A posteriori error estimation: a recovery-based indicator plus marking.

The piece that closes the adaptive loop: :func:`meshioplusplus.gradient` already
differentiates a field and selective :func:`meshioplusplus.refine`'s ``where``
predicate already consumes any scalar ``cell_data`` array; this is what produces
one.

:func:`estimate_error` implements the Zienkiewicz-Zhu (ZZ) superconvergent-patch
recovery estimator: differentiate the field once per cell (``gradient``,
Green-Gauss, cell location), recover a smoothed nodal gradient by averaging onto
the points and back onto the cells (``cell_data_to_point_data`` /
``point_data_to_cell_data``, measure-weighted), and take the
recovered-minus-raw difference as the local error in the gradient. **This is
composition, not a new numerical kernel** -- the recovery operator *is* the
existing averaging round trip, reused rather than reimplemented, exactly as
``gradient(location="point")`` itself composes ``cell_data_to_point_data``.

For cell ``K`` with raw gradient ``g_K`` and recovered gradient ``G_K``:
``eta_K = sqrt(|measure_K| * sum((G_K - g_K)^2))`` -- the standard ZZ local
indicator in the energy-like (H1-seminorm) norm. ``mGlobalError`` (returned in
the report) is ``sqrt(sum_K eta_K^2)`` over the evaluable cells.

A marking pass turns the continuous indicator into a boolean ``cell_data``
array, ``error:marked`` (one array per block, int64 0/1), so that ``refine``
needs no change at all -- the intended use is
``refine(mesh, where="error:marked > 0.5")``. Three policies: ``"absolute"``
(threshold), ``"fraction"`` (the largest fraction of cells), ``"dorfler"`` (the
smallest bulk-covering prefix, the Dörfler / bulk-chunk criterion).

Cells that cannot be evaluated -- an unsupported/degenerate ``gradient`` cell, a
recovery neighbourhood with no finite contribution, or an unmeasurable cell --
read NaN in ``error:zz`` and 0 (never NaN) in ``error:marked``, and are excluded
from the global error and from every marking policy's cell count.

The C++ core (``_core.estimate_error``) does the work; this module is the thin
shim (try C++, fall back to a pure-numpy reference) plus the reference itself.
The reference matches the C++ result to within a numeric tolerance for a mesh of
purely linear 3D cell types (tetra/hexahedron/wedge/pyramid) -- the same,
already-accepted tolerance ``data_average``'s own measure weighting has,
inherited here because the recovery step composes it -- and raises
``NotImplementedError`` for a ragged/polyhedron block or a quadratic 3D type
(``tetra10``, ``hexahedron20``/``27``, ``wedge15``): ``detail::cell_measure``
computes those through ``orient_rings``' winding repair (quadratic types) or
through the polyhedral kernel (ragged/polyhedron blocks) -- discrete branches,
or coverage, this reference does not attempt to replicate, the same precedent
``_convert_cells.py``'s simplexify twin uses.

Public API:

* :func:`estimate_error` -- the ZZ indicator plus marking.
"""

from __future__ import annotations

import numpy as np

from ._data_average import (
    _DIM,
    _cell_measures,
    cell_data_to_point_data,
    point_data_to_cell_data,
)
from ._gradient import gradient

_PREFIX = "meshio++: estimate_error: "

_METHOD_ALIASES = {"": "zz", "zz": "zz"}
_MARKING_ALIASES = {
    "": "none",
    "none": "none",
    "absolute": "absolute",
    "abs": "absolute",
    "fraction": "fraction",
    "frac": "fraction",
    "dorfler": "dorfler",
    "bulk": "dorfler",
}

#: Private working names for the recovery composition -- never reach the
#: returned mesh (it is built by dropping them from the composed one and
#: attaching the real outputs), so a clash with a real user array only
#: affects an internal, discarded intermediate copy.
_RAW_NAME = "__estimate_error_raw__"
_REC_SUFFIX = "_recovered"

#: The linear 3D cell types `_data_average._cell_measures`' face table covers.
#: A quadratic 3D type (`tetra10`, `hexahedron20`/`27`, `wedge15`) is a KNOWN
#: gap in that reference table (not in the C++ core, which reduces to corners
#: exactly as `gradient` does) -- see the module docstring.
_TABULATED_3D = {"tetra", "hexahedron", "wedge", "pyramid"}


def _validate(
    mesh, array, method, marking, marking_value, output, marked_name, overwrite
):
    """Raise on every user error, before either path runs.

    Shared deliberately: the C++ core validates the same things, but the numpy
    fallback must reject the same inputs with the same message even when
    ``_core`` is unavailable.
    """
    if method not in _METHOD_ALIASES:
        raise ValueError(f"{_PREFIX}unknown method '{method}' (expected 'zz')")
    if marking not in _MARKING_ALIASES:
        raise ValueError(
            f"{_PREFIX}unknown marking policy '{marking}' "
            "(expected 'none', 'absolute', 'fraction', or 'dorfler')"
        )
    mk = _MARKING_ALIASES[marking]

    if not array:
        raise ValueError(f"{_PREFIX}an array name is required")
    if array not in mesh.point_data:
        # A cell_data field is piecewise constant and has no derivative to
        # recover. Name the fix, exactly as gradient does for the same reason.
        if array in mesh.cell_data:
            raise ValueError(
                f"{_PREFIX}'{array}' is cell_data, which is piecewise constant "
                "and has no derivative to recover; convert it first with "
                "cell_data_to_point_data (CLI: meshioplusplus data to-point)"
            )
        have = ", ".join(sorted(mesh.point_data))
        raise ValueError(
            f"{_PREFIX}no point_data array named '{array}'"
            + (f" (available: {have})" if have else " (the mesh has no point_data)")
        )

    if mk in ("fraction", "dorfler") and not (0.0 < marking_value <= 1.0):
        raise ValueError(
            f"{_PREFIX}marking_value must be in (0, 1] for the "
            "'fraction'/'dorfler' marking policy"
        )

    zz_name = output or "error:zz"
    mk_name = marked_name or "error:marked"
    if not overwrite:
        if zz_name in mesh.cell_data:
            raise ValueError(
                f"{_PREFIX}'{zz_name}' already exists in cell_data "
                "(pass overwrite=True to replace it)"
            )
        if mk != "none" and mk_name in mesh.cell_data:
            raise ValueError(
                f"{_PREFIX}'{mk_name}' already exists in cell_data "
                "(pass overwrite=True to replace it)"
            )


def _twin_unsupported(mesh):
    """Whether the numpy reference cannot faithfully estimate this mesh.

    True for a ragged/polygon/polyhedron block, or a *known* 3D cell type that
    `_cell_measures` cannot weight (the quadratic family) -- see the module
    docstring. A genuinely unknown/exotic type is left alone: both the C++
    core and this reference already skip it gracefully (NaN, counted), so
    there is nothing for the two to disagree about.
    """
    for block in mesh.cells:
        if isinstance(block.data, list):
            return True
        ctype = block.type
        if ctype.startswith(("polygon", "polyhedron")):
            return True
        if _DIM.get(ctype) == 3 and ctype not in _TABULATED_3D:
            return True
    return False


def _rank(values, indices):
    """Descending by value, ascending index on a tie -- mirrors error.cpp's
    `err_greater`, a strict weak order independent of input order."""
    order = sorted(range(len(values)), key=lambda i: (-values[i], indices[i]))
    return [indices[i] for i in order], [values[i] for i in order]


def _estimate_error_py(
    mesh, array, marking, marking_value, output, marked_name, overwrite
):
    """Pure-numpy ZZ estimator + marking. See the module docstring for scope."""
    if _twin_unsupported(mesh):
        raise NotImplementedError(
            "estimate_error: ragged/polygon/polyhedron blocks and quadratic 3D "
            "cell types are C++-core only"
        )

    zz_name = output or "error:zz"
    mk_name = marked_name or "error:marked"

    g = gradient(
        mesh,
        array,
        operator="gradient",
        method="green-gauss",
        location="cell",
        output=_RAW_NAME,
        overwrite=True,
    )
    pt = cell_data_to_point_data(g, keys=[_RAW_NAME], weighted=True, overwrite=True)
    rec = point_data_to_cell_data(
        pt, keys=[_RAW_NAME], suffix=_REC_SUFFIX, overwrite=True
    )
    rec_name = _RAW_NAME + _REC_SUFFIX

    global_zz = []
    zz_blocks = []
    skipped = 0
    for b, block in enumerate(mesh.cells):
        n = len(block)
        raw = np.asarray(rec.cell_data[_RAW_NAME][b], dtype=float).reshape(n, -1)
        recv = np.asarray(rec.cell_data[rec_name][b], dtype=float).reshape(n, -1)
        finite = np.isfinite(raw).all(axis=1) & np.isfinite(recv).all(axis=1)
        sumsq = np.where(finite, np.sum((recv - raw) ** 2, axis=1), np.nan)

        measure = _cell_measures(mesh, block)
        if measure is None:
            measure = np.full(n, np.nan)
        measure = np.abs(measure)

        zz = np.where(finite & np.isfinite(measure), np.sqrt(measure * sumsq), np.nan)
        skipped += int(np.count_nonzero(~np.isfinite(zz)))
        zz_blocks.append(zz)
        global_zz.append(zz)

    flat = np.concatenate(global_zz) if global_zz else np.zeros(0)
    finite_mask = np.isfinite(flat)
    global_error = (
        float(np.sqrt(np.sum(flat[finite_mask] ** 2))) if finite_mask.any() else 0.0
    )

    out = rec
    out.cell_data.pop(_RAW_NAME, None)
    out.cell_data.pop(rec_name, None)
    out.cell_data[zz_name] = zz_blocks

    marked_count = 0
    mk = _MARKING_ALIASES[marking]
    if mk != "none":
        mark = np.zeros(flat.shape[0], dtype=np.int64)
        if mk == "absolute":
            mark = np.where(finite_mask & (flat > marking_value), 1, 0).astype(np.int64)
        else:
            idx = np.nonzero(finite_mask)[0]
            vals = flat[idx]
            order, vals_sorted = _rank(list(vals), list(idx))
            if mk == "fraction":
                take = int(marking_value * len(order) + 0.5)
                take = min(take, len(order))
            else:  # dorfler
                total_sq = float(np.sum(np.asarray(vals_sorted) ** 2))
                target = marking_value * total_sq
                cum = 0.0
                take = 0
                while take < len(order):
                    if cum >= target:
                        break
                    cum += vals_sorted[take] ** 2
                    take += 1
            for i in range(take):
                mark[order[i]] = 1
        marked_count = int(mark.sum())

        offsets = np.cumsum([0] + [len(b) for b in mesh.cells])
        marked_blocks = [
            mark[offsets[b] : offsets[b + 1]] for b in range(len(mesh.cells))
        ]
        out.cell_data[mk_name] = marked_blocks

    return out, global_error, skipped, marked_count


def estimate_error(
    mesh,
    array,
    method="zz",
    marking="none",
    marking_value=0.0,
    output=None,
    marked_name=None,
    overwrite=False,
    return_report=False,
):
    """Estimate the per-cell recovered-gradient error of a ``point_data``
    field and optionally mark cells for refinement.

    :param mesh: the source mesh (never modified).
    :param array: name of the ``point_data`` array to estimate the error of.
    :param method: estimator family; only ``"zz"`` (default) exists today.
    :param marking: ``"none"`` (default), ``"absolute"``, ``"fraction"``, or
        ``"dorfler"``.
    :param marking_value: meaning depends on ``marking``: ignored for
        ``"none"``; an absolute indicator threshold for ``"absolute"``; a
        fraction in ``(0, 1]`` of cells for ``"fraction"``; the Dörfler bulk
        fraction theta in ``(0, 1]`` for ``"dorfler"``.
    :param output: indicator output name; ``None`` selects ``"error:zz"``.
    :param marked_name: marking output name; ``None`` selects
        ``"error:marked"``. Ignored when ``marking == "none"``.
    :param overwrite: allow replacing an existing array of an output name.
    :param return_report: also return ``{global_error, num_skipped,
        num_marked}``.
    :returns: the mesh carrying the produced array(s), or ``(mesh, report)``
        when ``return_report`` is set.
    :raises ValueError: on an unknown or ``cell_data`` array name, an unknown
        method/marking policy, an out-of-range ``marking_value``, or an
        output-name collision without ``overwrite``.
    """
    _validate(
        mesh, array, method, marking, marking_value, output, marked_name, overwrite
    )

    out = None
    report = None
    try:
        from . import _core

        res = _core.estimate_error(
            mesh,
            array,
            _METHOD_ALIASES[method],
            _MARKING_ALIASES[marking],
            float(marking_value),
            "" if output is None else str(output),
            "" if marked_name is None else str(marked_name),
            bool(overwrite),
        )
        out = res["mesh"]
        report = {
            "global_error": res["global_error"],
            "num_skipped": res["num_skipped"],
            "num_marked": res["num_marked"],
        }
    except (ValueError, TypeError):
        # A genuine user error must not fall through to the numpy path, which
        # would either raise something less helpful or silently succeed
        # differently.
        raise
    except Exception:
        out = None

    if out is None:
        out, global_error, skipped, marked = _estimate_error_py(
            mesh,
            array,
            _MARKING_ALIASES[marking],
            marking_value,
            output,
            marked_name,
            overwrite,
        )
        report = {
            "global_error": global_error,
            "num_skipped": skipped,
            "num_marked": marked,
        }

    return (out, report) if return_report else out

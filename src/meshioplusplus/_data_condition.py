"""Value conditioning for data arrays: clamp, normalize, standardize.

A dependency-free mesh *operation* (not a file format). The geometry is never
modified; only the values of the selected arrays change.

- ``clamp`` maps each value to ``min(max(x, lo), hi)``.
- ``normalize`` maps the array's own ``[min, max]`` linearly onto ``[lo, hi]``
  (``[0, 1]`` by default).
- ``standardize`` maps to zero mean and unit standard deviation (population
  standard deviation, divisor N).

``scope="component"`` (the default) conditions each trailing component on its
own statistics; ``scope="magnitude"`` uses each row's Euclidean magnitude and
rescales whole rows, preserving direction.

For ``cell_data`` the statistics are computed **jointly over all cell blocks**,
then one transform is applied to every block.

Non-finite values are always excluded from the statistics; ``nan_policy``
decides what reaches the output. See :mod:`meshioplusplus._data_common`.

Public API:
    data_condition
"""

from __future__ import annotations

import copy
import warnings

import numpy as np

from ._data_common import location_map, normalize_location, require_key


def _components(array):
    a = np.asarray(array)
    return 1 if a.ndim < 2 else int(np.prod(a.shape[1:]))


def _finite_stats(values, ncomp):
    """Per-component (min, max, mean, std) over the finite values."""
    flat = np.asarray(values, dtype=float).reshape(-1, ncomp)
    mask = np.isfinite(flat)
    counts = mask.sum(axis=0)
    safe = np.where(mask, flat, np.nan)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        mins = np.nanmin(safe, axis=0)
        maxs = np.nanmax(safe, axis=0)
        means = np.nanmean(safe, axis=0)
        stds = np.nanstd(safe, axis=0)
    mins = np.where(counts > 0, mins, np.nan)
    maxs = np.where(counts > 0, maxs, np.nan)
    means = np.where(counts > 0, means, np.nan)
    stds = np.where(counts > 0, stds, np.nan)
    return mins, maxs, means, stds, counts


def _magnitudes(values, ncomp):
    """Row magnitudes; NaN for any row containing a non-finite component."""
    flat = np.asarray(values, dtype=float).reshape(-1, ncomp)
    finite = np.isfinite(flat).all(axis=1)
    mags = np.sqrt((np.where(np.isfinite(flat), flat, 0.0) ** 2).sum(axis=1))
    return np.where(finite, mags, np.nan)


def _build(mode, lo, hi, mins, maxs, means, stds, counts, name):
    """Return per-component (scale, offset) for the affine modes."""
    if mode == "normalize":
        span = maxs - mins
        ok = (counts > 0) & np.isfinite(span) & (span > 0)
        if not np.all(ok):
            warnings.warn(
                f"meshio++: data_condition: '{name}' has no value range to "
                "normalize (constant or non-finite); filling with the target "
                "lower bound",
                RuntimeWarning,
                stacklevel=2,
            )
        scale = np.where(ok, (hi - lo) / np.where(ok, span, 1.0), 0.0)
        offset = np.where(ok, lo - scale * mins, lo)
        return scale, offset
    if mode == "standardize":
        ok = (counts > 0) & np.isfinite(stds) & (stds > 0)
        if not np.all(ok):
            warnings.warn(
                f"meshio++: data_condition: '{name}' has zero standard "
                "deviation; filling with 0",
                RuntimeWarning,
                stacklevel=2,
            )
        scale = np.where(ok, 1.0 / np.where(ok, stds, 1.0), 0.0)
        offset = np.where(ok, -means * scale, 0.0)
        return scale, offset
    raise ValueError(
        f"meshio++: unknown conditioning mode '{mode}' "
        "(expected 'clamp', 'normalize' or 'standardize')"
    )


def _apply_nan_policy(values, nan_policy, nan_replacement, name):
    if nan_policy == "ignore":
        return values
    bad = ~np.isfinite(values)
    if not bad.any():
        return values
    if nan_policy == "replace":
        out = values.copy()
        out[bad] = nan_replacement
        return out
    if nan_policy == "fail":
        idx = int(np.argmax(bad.ravel()))
        raise ValueError(
            f"meshio++: data_condition: non-finite value at index {idx} of '{name}'"
        )
    raise ValueError(
        f"meshio++: unknown NaN policy '{nan_policy}' "
        "(expected 'ignore', 'replace' or 'fail')"
    )


def _condition_py(
    mesh,
    location,
    names,
    mode,
    scope,
    lo,
    hi,
    nan_policy,
    nan_replacement,
    suffix,
    preserve_dtype,
):
    """Pure-Python reference for :func:`data_condition`."""
    loc = normalize_location(location)
    if mode in ("clamp", "normalize") and lo > hi:
        raise ValueError(
            f"meshio++: data_condition: bounds are inverted (lo={lo} > hi={hi})"
        )
    source = location_map(mesh, loc)
    selected = list(names) if names else sorted(source.keys())
    for name in selected:
        require_key(mesh, loc, name)

    out = copy.deepcopy(mesh)
    magnitude = scope == "magnitude"
    if scope not in ("component", "magnitude"):
        raise ValueError(
            f"meshio++: unknown conditioning scope '{scope}' "
            "(expected 'component' or 'magnitude')"
        )

    for name in selected:
        target = name + suffix
        if loc == "cell_data":
            blocks = source[name]
            if len(blocks) != len(mesh.cells):
                raise ValueError(
                    f"meshio++: data_condition: cell_data '{name}' has "
                    f"{len(blocks)} block(s) but the mesh has "
                    f"{len(mesh.cells)} cell block(s)"
                )
            arrays = [np.asarray(b) for b in blocks]
        else:
            arrays = [np.asarray(source[name])]

        ncomp = _components(arrays[0])
        pooled = (
            np.concatenate([_magnitudes(a, ncomp) for a in arrays]).reshape(-1, 1)
            if magnitude
            else np.concatenate(
                [np.asarray(a, dtype=float).reshape(-1, ncomp) for a in arrays]
            )
        )
        stat_comp = 1 if magnitude else ncomp
        mins, maxs, means, stds, counts = _finite_stats(pooled, stat_comp)
        if mode != "clamp":
            scale, offset = _build(mode, lo, hi, mins, maxs, means, stds, counts, name)

        results = []
        for a in arrays:
            rows = a.shape[0] if a.ndim else 0
            flat = np.asarray(a, dtype=float).reshape(rows, ncomp)
            if magnitude:
                mags = _magnitudes(a, ncomp)
                if mode == "clamp":
                    target_mag = np.clip(mags, lo, hi)
                else:
                    target_mag = scale[0] * mags + offset[0]
                with np.errstate(invalid="ignore", divide="ignore"):
                    factor = np.where(
                        mags > 0, target_mag / np.where(mags > 0, mags, 1.0), 0.0
                    )
                res = flat * factor.reshape(-1, 1)
                res = np.where(np.isfinite(mags).reshape(-1, 1), res, flat)
            elif mode == "clamp":
                res = np.where(np.isfinite(flat), np.clip(flat, lo, hi), flat)
            else:
                res = np.where(np.isfinite(flat), flat * scale + offset, flat)

            res = _apply_nan_policy(res, nan_policy, nan_replacement, name)
            if mode == "clamp" and preserve_dtype:
                res = res.astype(np.asarray(a).dtype)
            shaped = (
                res.reshape(np.asarray(a).shape) if ncomp > 1 else res.reshape(rows)
            )
            results.append(shaped)

        if loc == "cell_data":
            out.cell_data[target] = results
        elif loc == "point_data":
            out.point_data[target] = results[0]
        else:
            out.field_data[target] = results[0]
    return out


def data_condition(
    mesh,
    location: str = "point",
    keys=None,
    mode: str = "clamp",
    scope: str = "component",
    lo: float = 0.0,
    hi: float = 1.0,
    nan_policy: str = "ignore",
    nan_replacement: float = 0.0,
    suffix: str = "",
    preserve_dtype: bool = True,
):
    """Condition the values of the selected data arrays.

    Args:
        mesh: the source mesh (unmodified).
        location: ``"point"``, ``"cell"`` or ``"field"``.
        keys: names to condition; ``None`` means every array at ``location``.
        mode: ``"clamp"``, ``"normalize"`` or ``"standardize"``.
        scope: ``"component"`` or ``"magnitude"``.
        lo: clamp lower bound, or normalize target lower bound.
        hi: clamp upper bound, or normalize target upper bound.
        nan_policy: ``"ignore"``, ``"replace"`` or ``"fail"``.
        nan_replacement: value used when ``nan_policy`` is ``"replace"``.
        suffix: empty replaces the array in place; otherwise the result is
            stored as ``name + suffix``.
        preserve_dtype: ``clamp`` only -- keep the input dtype. ``normalize``
            and ``standardize`` always produce ``float64``.

    Returns:
        A new mesh with the conditioned arrays.
    """
    loc = normalize_location(location)
    names = list(keys) if keys else []
    args = (
        loc,
        names,
        mode,
        scope,
        lo,
        hi,
        nan_policy,
        nan_replacement,
        suffix,
        preserve_dtype,
    )
    try:
        from . import _core
    except Exception:
        return _condition_py(mesh, *args)
    try:
        out = _core.data_condition(
            mesh,
            loc,
            names,
            mode,
            scope,
            lo,
            hi,
            nan_policy,
            nan_replacement,
            suffix,
            preserve_dtype,
        )
    except ValueError:
        raise
    except Exception:
        return _condition_py(mesh, *args)
    for attr in ("point_sets", "cell_sets"):
        value = getattr(mesh, attr, None)
        if value:
            setattr(out, attr, copy.deepcopy(value))
    return out

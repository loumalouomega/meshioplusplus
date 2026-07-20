"""Read-only per-array summary of everything a mesh's data maps carry.

A dependency-free mesh *operation* (not a file format), and the *data* view that
complements the topological ``info`` verb and the geometric
:func:`meshioplusplus.compute_stats`. For every ``point_data``, ``cell_data``
and ``field_data`` array it reports the location, dtype, shape and component
count, the number of entries, min / max / mean (whole-array and per component),
and the counts of NaN and infinite values.

Unlike the other data operations this one never raises on non-finite values --
it counts them, which is the point.

Public API:
    data_info
"""

from __future__ import annotations

import warnings

import numpy as np


def _summarize(location, name, arrays, num_blocks):
    """Build one summary dict from the array(s) making up a single entry."""
    first = np.asarray(arrays[0])
    ncomp = 1 if first.ndim < 2 else int(np.prod(first.shape[1:]))

    inconsistent = False
    parts = []
    num_values = 0
    num_entries = 0
    for a in arrays:
        arr = np.asarray(a)
        this = 1 if arr.ndim < 2 else int(np.prod(arr.shape[1:]))
        if this != ncomp:
            # Reported rather than raised: data_info is a diagnostic.
            inconsistent = True
            continue
        flat = arr.reshape(-1, ncomp).astype(float, copy=False)
        parts.append(flat)
        num_values += int(arr.size)
        num_entries += int(flat.shape[0])

    pooled = np.concatenate(parts) if parts else np.empty((0, ncomp), dtype=float)
    finite = np.isfinite(pooled)
    nan_count = int(np.isnan(pooled).sum())
    inf_count = int(np.isinf(pooled).sum())
    finite_count = int(finite.sum())

    safe = np.where(finite, pooled, np.nan)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        if pooled.size:
            per_min = np.nanmin(safe, axis=0)
            per_max = np.nanmax(safe, axis=0)
            per_mean = np.nanmean(safe, axis=0)
            counts = finite.sum(axis=0)
            per_min = np.where(counts > 0, per_min, np.nan)
            per_max = np.where(counts > 0, per_max, np.nan)
            per_mean = np.where(counts > 0, per_mean, np.nan)
            overall_min = float(np.nanmin(per_min)) if finite_count else float("nan")
            overall_max = float(np.nanmax(per_max)) if finite_count else float("nan")
            overall_mean = (
                float(np.nansum(safe) / finite_count) if finite_count else float("nan")
            )
        else:
            per_min = np.full(ncomp, np.nan)
            per_max = np.full(ncomp, np.nan)
            per_mean = np.full(ncomp, np.nan)
            overall_min = float("nan")
            overall_max = float("nan")
            overall_mean = float("nan")

    return {
        "location": location,
        "name": name,
        "dtype": np.dtype(first.dtype).str.lstrip("<>|="),
        "shape": tuple(int(d) for d in first.shape),
        "num_blocks": int(num_blocks),
        "num_entries": num_entries,
        "num_components": ncomp,
        "num_values": num_values,
        "min": overall_min,
        "max": overall_max,
        "mean": overall_mean,
        "min_per_component": [float(v) for v in per_min],
        "max_per_component": [float(v) for v in per_max],
        "mean_per_component": [float(v) for v in per_mean],
        "num_nan": nan_count,
        "num_inf": inf_count,
        "num_finite": finite_count,
        "inconsistent_blocks": inconsistent,
    }


def _info_py(mesh) -> list:
    """Pure-Python reference for :func:`data_info`."""
    out = []
    for name in sorted(mesh.point_data.keys()):
        out.append(_summarize("point_data", name, [mesh.point_data[name]], 1))
    for name in sorted(mesh.cell_data.keys()):
        blocks = mesh.cell_data[name]
        if not blocks:
            continue
        out.append(_summarize("cell_data", name, list(blocks), len(blocks)))
    for name in sorted(mesh.field_data.keys()):
        out.append(_summarize("field_data", name, [mesh.field_data[name]], 1))
    return out


def data_info(mesh) -> list:
    """Summarize every data array ``mesh`` carries (read-only).

    Args:
        mesh: the mesh to inspect (unmodified).

    Returns:
        A list of dicts, one per array, ordered ``point_data`` then
        ``cell_data`` then ``field_data``, each group sorted by name. Each dict
        carries ``location``, ``name``, ``dtype``, ``shape``, ``num_blocks``,
        ``num_entries``, ``num_components``, ``num_values``, ``min``, ``max``,
        ``mean``, the three ``*_per_component`` lists, ``num_nan``,
        ``num_inf``, ``num_finite`` and ``inconsistent_blocks``.
    """
    try:
        from . import _core

        raw = _core.data_info(mesh)
    except Exception:
        return _info_py(mesh)

    # Re-type so both paths return byte-identical structures.
    out = []
    for a in raw:
        out.append(
            {
                "location": str(a["location"]),
                "name": str(a["name"]),
                "dtype": str(a["dtype"]),
                "shape": tuple(int(d) for d in a["shape"]),
                "num_blocks": int(a["num_blocks"]),
                "num_entries": int(a["num_entries"]),
                "num_components": int(a["num_components"]),
                "num_values": int(a["num_values"]),
                "min": float(a["min"]),
                "max": float(a["max"]),
                "mean": float(a["mean"]),
                "min_per_component": [float(v) for v in a["min_per_component"]],
                "max_per_component": [float(v) for v in a["max_per_component"]],
                "mean_per_component": [float(v) for v in a["mean_per_component"]],
                "num_nan": int(a["num_nan"]),
                "num_inf": int(a["num_inf"]),
                "num_finite": int(a["num_finite"]),
                "inconsistent_blocks": bool(a["inconsistent_blocks"]),
            }
        )
    return out

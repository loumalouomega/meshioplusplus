"""Statistics and derivatives over the steps of a sequence (roadmap 15.5).

Two streaming operations over a time series (a glob, a list of paths or a
:class:`~meshioplusplus._sequence.TimeSeries`):

* :func:`sequence_statistics` reduces every step to one mesh carrying, per data
  array, its **mean**, **std**, **min**, **max** and/or **sum** over the steps
  (Welford's update, one step in memory at a time), the counterpart of VTK's
  ``vtkTemporalStatistics``;
* :func:`sequence_derivative` yields every step with the **time derivative** of
  one array, from the step times, holding at most three steps.

A **Python-only layer** (the sequence readers it builds on are Python's; the
flat bindings already traverse sequences step by step through
``mio_sequence_read``). Every step must have the same points (and cells, for
``cell_data``) and the arrays the same shape: a step that differs is refused by
name, never skipped. See ``doc/sequence_stats.md``.
"""

from __future__ import annotations

import numpy as np

from ._data_common import location_map, normalize_location
from ._sequence import read_sequence

__all__ = ["sequence_statistics", "sequence_derivative", "STATISTICS"]

#: The statistics, in the order the output arrays are written.
STATISTICS = ("mean", "std", "min", "max", "sum")


def _arrays_of(mesh, loc):
    """``{name: float64 array}`` of the numeric arrays at a location; a
    ``cell_data`` array is its blocks concatenated (the block shapes are checked
    through the total)."""
    out = {}
    for name, value in location_map(mesh, loc).items():
        if loc == "cell_data":
            blocks = [np.asarray(b) for b in value]
            if not blocks or any(b.dtype.kind not in "biuf" for b in blocks):
                continue
            out[name] = (
                np.concatenate([b.reshape(len(b), -1) for b in blocks]).astype(float),
                [b.shape for b in blocks],
            )
        else:
            arr = np.asarray(value)
            if arr.dtype.kind not in "biuf":
                continue
            out[name] = (arr.astype(float), [arr.shape])
    return out


def _split_blocks(values, shapes):
    """The inverse of the concatenation of :func:`_arrays_of` for ``cell_data``."""
    blocks, at = [], 0
    for shape in shapes:
        n = shape[0]
        blocks.append(values[at : at + n].reshape(shape))
        at += n
    return blocks


def sequence_statistics(
    source,
    arrays=None,
    stats=("mean", "std"),
    *,
    location="point",
    ddof=0,
    file_format=None,
    times=None,
    time_from="auto",
    sort=False,
    return_report=False,
    **read_kwargs,
):
    """Per-point (or per-cell) statistics of data arrays over a sequence.

    Returns the **first step's mesh with its data replaced by the statistics**:
    for each array ``a`` and each requested statistic ``s``, ``<a>:<s>`` at the
    location (a ``cell_data`` statistic has one array per cell block, like any
    cell array). The reduction streams one step at a time, with Welford's
    update for the mean and the variance, so memory is about the arrays of two
    steps however long the sequence is.

    :param source: a glob pattern, a list of paths or a ``TimeSeries``.
    :param arrays: the array names (default every numeric array at the
        location of the first step).
    :param stats: any of ``"mean"``, ``"std"``, ``"min"``, ``"max"``, ``"sum"``.
    :param location: ``"point"`` or ``"cell"``.
    :param ddof: the delta degrees of freedom of ``std`` (0, the population
        standard deviation, as ``numpy.std``; 1 for the sample's).
    :param return_report: also return ``{"num_steps", "first_time",
        "last_time", "arrays"}``.
    :raises ValueError: on an unknown statistic or array, an empty sequence, or
        a step whose points, cells or array shapes differ from the first.

    A non-finite value propagates into the statistics of its point (``mean``,
    ``std`` and ``sum`` become ``NaN``; ``min`` and ``max`` follow
    ``numpy.minimum``): there is no filtering, so a ``NaN`` is visible rather
    than silently dropped.
    """
    loc = normalize_location(location)
    if loc == "field_data":
        raise ValueError(
            "meshio++: sequence_statistics: field_data has no per-point or per-cell "
            "correspondence to reduce over time"
        )
    stats = (stats,) if isinstance(stats, str) else tuple(stats)
    for s in stats:
        if s not in STATISTICS:
            raise ValueError(
                f"meshio++: sequence_statistics: unknown statistic '{s}' "
                f"(expected one of {', '.join(STATISTICS)})"
            )
    if not stats:
        raise ValueError("meshio++: sequence_statistics: no statistic requested")

    first_mesh = None
    names = None
    shapes = {}
    n = 0
    mean, m2, mn, mx, total = {}, {}, {}, {}, {}
    first_time = last_time = None
    for t, mesh in read_sequence(
        source,
        file_format=file_format,
        times=times,
        time_from=time_from,
        sort=sort,
        **read_kwargs,
    ):
        current = _arrays_of(mesh, loc)
        if first_mesh is None:
            first_mesh = mesh
            first_time = t
            names = list(arrays) if arrays else sorted(current)
            for name in names:
                if name not in current:
                    raise ValueError(
                        f"meshio++: sequence_statistics: no numeric {loc} array "
                        f"'{name}' in the first step (available: "
                        f"{', '.join(sorted(current)) or 'none'})"
                    )
            if not names:
                raise ValueError(
                    f"meshio++: sequence_statistics: the first step has no numeric {loc}"
                )
            for name in names:
                values, shapes[name] = current[name]
                mean[name] = np.zeros_like(values)
                m2[name] = np.zeros_like(values)
                mn[name] = np.full_like(values, np.inf)
                mx[name] = np.full_like(values, -np.inf)
                total[name] = np.zeros_like(values)
            npoints = len(mesh.points)
            ncells = [len(b.data) for b in mesh.cells]
        else:
            if (
                len(mesh.points) != npoints
                or [len(b.data) for b in mesh.cells] != ncells
            ):
                raise ValueError(
                    f"meshio++: sequence_statistics: step {n} (time {t}) has a different "
                    "number of points or cells from the first step; the statistics need "
                    "one topology"
                )
        n += 1
        last_time = t
        for name in names:
            if name not in current:
                raise ValueError(
                    f"meshio++: sequence_statistics: step {n - 1} (time {t}) has no "
                    f"{loc} array '{name}'"
                )
            values, shape = current[name]
            if shape != shapes[name]:
                raise ValueError(
                    f"meshio++: sequence_statistics: '{name}' changes shape at step "
                    f"{n - 1} (time {t}): {shape} against {shapes[name]}"
                )
            delta = values - mean[name]
            mean[name] = mean[name] + delta / n
            m2[name] = m2[name] + delta * (values - mean[name])
            with np.errstate(invalid="ignore"):
                mn[name] = np.minimum(mn[name], values)
                mx[name] = np.maximum(mx[name], values)
            total[name] = total[name] + values
    if first_mesh is None:
        raise ValueError("meshio++: sequence_statistics: the sequence has no steps")

    out = first_mesh.copy() if hasattr(first_mesh, "copy") else first_mesh
    location_map(out, loc).clear()
    for name in names:
        denom = n - ddof
        results = {
            "mean": mean[name],
            "std": (
                np.sqrt(m2[name] / denom)
                if denom > 0
                else np.full_like(m2[name], np.nan)
            ),
            "min": mn[name],
            "max": mx[name],
            "sum": total[name],
        }
        for s in stats:
            value = results[s]
            if loc == "cell_data":
                location_map(out, loc)[f"{name}:{s}"] = _split_blocks(
                    value, shapes[name]
                )
            else:
                location_map(out, loc)[f"{name}:{s}"] = value.reshape(shapes[name][0])
    report = {
        "num_steps": n,
        "first_time": first_time,
        "last_time": last_time,
        "arrays": list(names),
    }
    return (out, report) if return_report else out


def _difference(f_prev, f_cur, f_next, t_prev, t_cur, t_next, scheme):
    """The derivative at ``t_cur`` from its neighbours (``None`` for a missing one)."""
    if scheme == "forward":
        if f_next is not None:
            return (f_next - f_cur) / (t_next - t_cur)
        return (f_cur - f_prev) / (t_cur - t_prev)  # the last step: backward
    if f_prev is None:
        return (f_next - f_cur) / (t_next - t_cur)
    if f_next is None:
        return (f_cur - f_prev) / (t_cur - t_prev)
    # The second-order three-point formula for unequal steps, which is the central
    # difference when they are equal.
    h1, h2 = t_cur - t_prev, t_next - t_cur
    return (
        -h2 / (h1 * (h1 + h2)) * f_prev
        + (h2 - h1) / (h1 * h2) * f_cur
        + h1 / (h2 * (h1 + h2)) * f_next
    )


def sequence_derivative(
    source,
    array,
    *,
    location="point",
    scheme="central",
    output=None,
    file_format=None,
    times=None,
    time_from="auto",
    sort=False,
    **read_kwargs,
):
    """Yield ``(time, mesh)`` for every step, the mesh carrying ``<array>:dt``.

    The derivative of ``array`` with respect to time, from the step times:

    * ``scheme="central"`` (default): the second-order three-point formula
      (the central difference on equal steps; correct on unequal ones), with a
      first-order one-sided difference at the first and the last step;
    * ``scheme="forward"``: ``(f[i + 1] - f[i]) / (t[i + 1] - t[i])``, a
      backward difference at the last step.

    The result of each step is the step's own mesh with the derivative added at
    the location, named ``output`` (default ``<array>:dt``). A generator holding
    at most three steps, so it streams a long sequence into
    :func:`~meshioplusplus._sequence.write_sequence`.

    :raises ValueError: for fewer than two steps, times that do not increase, a
        step whose points, cells or array shape differ, or an unknown array or
        scheme.
    """
    loc = normalize_location(location)
    if loc == "field_data":
        raise ValueError(
            "meshio++: sequence_derivative: field_data has no per-point or per-cell "
            "correspondence to differentiate"
        )
    if scheme not in ("central", "forward"):
        raise ValueError(
            f"meshio++: sequence_derivative: unknown scheme '{scheme}' "
            "(expected 'central' or 'forward')"
        )
    name = output or f"{array}:dt"

    def pairs():
        shape = None
        topo = None
        index = 0
        last_t = None
        for t, mesh in read_sequence(
            source,
            file_format=file_format,
            times=times,
            time_from=time_from,
            sort=sort,
            **read_kwargs,
        ):
            current = _arrays_of(mesh, loc)
            if array not in current:
                raise ValueError(
                    f"meshio++: sequence_derivative: no numeric {loc} array '{array}' "
                    f"at step {index} (available: {', '.join(sorted(current)) or 'none'})"
                )
            values, shapes = current[array]
            key = (len(mesh.points), [len(b.data) for b in mesh.cells])
            if topo is None:
                topo, shape = key, shapes
            elif key != topo or shapes != shape:
                raise ValueError(
                    f"meshio++: sequence_derivative: step {index} (time {t}) has a "
                    "different topology or array shape from the first step"
                )
            if last_t is not None and not t > last_t:
                raise ValueError(
                    f"meshio++: sequence_derivative: the times must increase "
                    f"(step {index} has {t} after {last_t})"
                )
            last_t = t
            index += 1
            yield t, mesh, values, shapes
        if index < 2:
            raise ValueError("meshio++: sequence_derivative: needs at least two steps")

    prev = cur = None
    for item in pairs():
        nxt = item
        if cur is not None:
            t_c, mesh_c, f_c, shapes_c = cur
            f_p = prev[2] if prev is not None else None
            t_p = prev[0] if prev is not None else None
            d = _difference(f_p, f_c, nxt[2], t_p, t_c, nxt[0], scheme)
            yield t_c, _with_derivative(mesh_c, loc, name, d, shapes_c)
        prev, cur = cur, nxt
    t_c, mesh_c, f_c, shapes_c = cur
    d = _difference(prev[2], f_c, None, prev[0], t_c, None, scheme)
    yield t_c, _with_derivative(mesh_c, loc, name, d, shapes_c)


def _with_derivative(mesh, loc, name, values, shapes):
    out = mesh.copy() if hasattr(mesh, "copy") else mesh
    if loc == "cell_data":
        location_map(out, loc)[name] = _split_blocks(values, shapes)
    else:
        location_map(out, loc)[name] = values.reshape(shapes[0])
    return out

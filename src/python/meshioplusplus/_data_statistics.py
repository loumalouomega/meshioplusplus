"""Descriptive statistics of data arrays, per region (roadmap 15.5).

``data_info`` summarises an array by its minimum, maximum and mean; this adds
what a study of a field wants beyond that -- the variance and standard
deviation, the sum, **quantiles** and a **histogram** -- for the whole mesh and
for every named region, the counterpart of VTK's ``vtkDescriptiveStatistics``
and ``vtkComputeQuantiles``. A **Python-only layer** (``data_info`` and its
public structs are unchanged, so no binding or ABI moves). See
``doc/data_statistics.md``.
"""

from __future__ import annotations

import math

import numpy as np

from ._data_average import _cell_measures, _cell_node_lists
from ._data_common import location_map, normalize_location

__all__ = ["data_statistics"]


def _stats(x, w=None, quantiles=(), bins=0, ddof=0):
    """Statistics of ``x`` (n, k) over its finite values, per component."""
    n, k = x.shape
    out = {
        key: []
        for key in (
            "count",
            "num_nan",
            "num_inf",
            "sum",
            "mean",
            "variance",
            "std",
            "min",
            "max",
        )
    }
    if quantiles:
        out["quantiles"] = {repr(float(q)): [] for q in quantiles}
    if bins:
        out["histogram"] = []
    for c in range(k):
        col = x[:, c]
        fin = np.isfinite(col)
        v = col[fin]
        out["count"].append(int(v.size))
        out["num_nan"].append(int(np.isnan(col).sum()))
        out["num_inf"].append(int(np.isinf(col).sum()))
        if v.size == 0:
            for key in ("sum", "mean", "variance", "std", "min", "max"):
                out[key].append(float("nan"))
            for q in out.get("quantiles", {}):
                out["quantiles"][q].append(float("nan"))
            if bins:
                out["histogram"].append({"edges": [], "counts": []})
            continue
        if w is not None:
            ww = w[fin]
            tot = float(ww.sum())
            mean = float((ww * v).sum() / tot) if tot > 0 else float("nan")
            var = float((ww * (v - mean) ** 2).sum() / tot) if tot > 0 else float("nan")
        else:
            mean = float(v.mean())
            denom = v.size - ddof
            var = float(((v - mean) ** 2).sum() / denom) if denom > 0 else float("nan")
        out["sum"].append(float(v.sum()))
        out["mean"].append(mean)
        out["variance"].append(var)
        out["std"].append(float(np.sqrt(var)) if not math.isnan(var) else float("nan"))
        out["min"].append(float(v.min()))
        out["max"].append(float(v.max()))
        for q in out.get("quantiles", {}):
            out["quantiles"][q].append(float(np.quantile(v, float(q))))
        if bins:
            lo, hi = float(v.min()), float(v.max())
            counts, edges = np.histogram(
                v, bins=bins, range=(lo, hi if hi > lo else lo + 1)
            )
            out["histogram"].append(
                {"edges": [float(e) for e in edges], "counts": [int(c) for c in counts]}
            )
    return out


def _rows(mesh, loc, name):
    """The array as ``(n, k)`` float rows (``cell_data`` blocks concatenated),
    and the block index of every row (for a cell array)."""
    value = location_map(mesh, loc)[name]
    if loc == "cell_data":
        blocks = [np.asarray(b) for b in value]
        if not blocks or any(b.dtype.kind not in "biuf" for b in blocks):
            return None
        rows = np.concatenate([b.reshape(len(b), -1) for b in blocks]).astype(float)
        return rows
    arr = np.asarray(value)
    if arr.dtype.kind not in "biuf":
        return None
    return arr.reshape(len(arr), -1).astype(float)


def _region_lookup(mesh, loc):
    """``[(name, index array)]`` of the regions that map onto the location:
    point regions for point data, cell regions for cell data and (through
    their cells' nodes) for point data; a side region has no entry. Where
    several regions share a name the one with the smallest ``(dim, tag)`` wins,
    as everywhere else."""
    wanted = {"point_data": ("point", "cell"), "cell_data": ("cell",)}[loc]
    by_name = {}
    for region in mesh.regions:
        if region.kind not in wanted:
            continue
        key = (region.name, region.kind)
        current = by_name.get(key)
        if current is None or (region.dim, region.tag) < (current.dim, current.tag):
            by_name[key] = region
    out = []
    node_lists = None
    for (name, kind), region in sorted(by_name.items()):
        entries = np.asarray(region.entries, dtype=np.int64).ravel()
        if loc == "cell_data" or kind == "point":
            out.append((name, entries))
            continue
        # a cell region on point data: the points of its cells, once each
        if node_lists is None:
            node_lists, bases = [], [0]
            for block in mesh.cells:
                node_lists.append(_cell_node_lists(block))
                bases.append(bases[-1] + len(block.data))
        points = set()
        for g in entries:
            for b in range(len(node_lists)):
                if bases[b] <= g < bases[b + 1]:
                    points.update(int(p) for p in node_lists[b][int(g) - bases[b]])
                    break
        out.append((name, np.array(sorted(points), dtype=np.int64)))
    return out


def data_statistics(
    mesh,
    arrays=None,
    *,
    location="point",
    per_region=True,
    quantiles=(0.05, 0.5, 0.95),
    bins=0,
    ddof=0,
    weighted=False,
):
    """Statistics of ``point_data`` or ``cell_data`` arrays, whole mesh and per region.

    Returns a list with one dict per array, ``{"name", "location",
    "num_components", "domain": stats, "regions": [{"name", **stats}, ...]}``;
    ``stats`` holds, per component and over the **finite** values (the others
    are counted, not dropped silently): ``count``, ``num_nan``, ``num_inf``,
    ``sum``, ``mean``, ``variance``, ``std``, ``min``, ``max``, ``quantiles``
    (a dict keyed by the quantile, linear interpolation as ``numpy.quantile``)
    and, with ``bins``, ``histogram`` (``edges`` and ``counts`` per component,
    over the component's own range).

    :param location: ``"point"`` or ``"cell"``.
    :param per_region: also the statistics of each named region that maps onto
        the location (a point region for point data; a cell region for cell data
        and, through its cells' nodes, for point data; a side region is
        skipped).
    :param quantiles: the quantiles in ``[0, 1]`` (empty for none).
    :param bins: histogram bins per component (0 for none).
    :param ddof: delta degrees of freedom of the variance (0 population).
    :param weighted: for ``cell_data``, weight ``mean``, ``variance`` and
        ``std`` by each cell's unsigned measure (so ``mean`` equals the mean of
        :func:`~meshioplusplus.data_integrate`); the quantiles, the histogram
        and the extremes stay unweighted, and ``ddof`` is ignored.
    :raises ValueError: on an unknown array, ``field_data``, a quantile outside
        ``[0, 1]`` or ``weighted`` outside ``cell_data``.
    """
    loc = normalize_location(location)
    if loc == "field_data":
        raise ValueError(
            "meshio++: data_statistics: field_data has no per-point or per-cell "
            "correspondence to summarise per region"
        )
    quantiles = tuple(float(q) for q in quantiles)
    if any(not 0.0 <= q <= 1.0 for q in quantiles):
        raise ValueError("meshio++: data_statistics: quantiles must lie in [0, 1]")
    if bins < 0:
        raise ValueError("meshio++: data_statistics: bins must not be negative")
    if weighted and loc != "cell_data":
        raise ValueError(
            "meshio++: data_statistics: weighted needs the cell location (the weight "
            "is the cell's measure)"
        )
    available = location_map(mesh, loc)
    names = list(arrays) if arrays else sorted(available)
    for name in names:
        if name not in available:
            raise ValueError(
                f"meshio++: data_statistics: no {loc} array '{name}' "
                f"(available: {', '.join(sorted(available)) or 'none'})"
            )
    weights = None
    if weighted:
        parts = []
        for block in mesh.cells:
            m = _cell_measures(mesh, block)
            parts.append(
                np.full(len(block.data), 0.0)
                if m is None
                else np.where(np.isfinite(m) & (m > 0), m, 0.0)
            )
        weights = np.concatenate(parts) if parts else np.zeros(0)
    regions = _region_lookup(mesh, loc) if per_region else []
    out = []
    for name in names:
        rows = _rows(mesh, loc, name)
        if rows is None:
            continue
        entry = {
            "name": name,
            "location": loc.replace("_data", ""),
            "num_components": int(rows.shape[1]),
            "domain": _stats(rows, weights, quantiles, bins, ddof),
            "regions": [],
        }
        for rname, idx in regions:
            idx = idx[(idx >= 0) & (idx < len(rows))]
            sub = rows[idx]
            wsub = weights[idx] if weights is not None else None
            entry["regions"].append(
                {"name": rname, **_stats(sub, wsub, quantiles, bins, ddof)}
            )
        out.append(entry)
    return out

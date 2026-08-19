"""Cell-measure-weighted field integration: total, mean, and per-region
reductions over ``cell_data`` -- ``gradient``'s natural companion (``gradient``
differentiates a field; this integrates one).

A dependency-free mesh *operation* (not a file format); the mesh is never
modified, this only reports. It reads geometry (``detail::cell_measure`` on
the C++ side, ``_cell_measures`` here), so -- like ``gradient`` -- it lives
outside the no-geometry ``data_*`` bundle even though it is still reachable
as the ``meshioplusplus data integrate`` CLI verb.

Contract (see ``operations/data_integrate.hpp`` for the authoritative
version):

- ``cell_data`` only; a ``point_data`` name raises, naming
  ``point_data_to_cell_data`` (CLI ``data to-cell``) as the fix.
- Every sum is weighted by ``|_cell_measures(...)|``; a cell whose measure is
  not computable (ragged/unsupported/degenerate) is excluded from **both**
  the numerator and the denominator of every component -- never given a
  fallback weight of 1, unlike ``cell_data_to_point_data``'s own
  ``weighted=True`` mode, since a unit-weight substitution would corrupt a
  physical total rather than merely soften an average.
- A non-finite *value* excludes a cell from that component's numerator and
  denominator too (not just zeroed into the numerator), counted per
  component.
- Per component only -- no whole-array (cross-component) total or mean,
  since summing unrelated physical components has no general meaning.
- No ``nan_policy``: like ``gradient`` and ``data_info``, this is a reduction
  with nothing to exclude a value *from*.
- Every named ``Cell`` region present gets its own independent entry
  (mirroring ``split(by="regions")``'s contract): regions are not a
  partition, so a cell in two regions contributes fully to both, and a cell
  in none contributes to neither. ``Point``/``Side`` regions are skipped.

The C++ core (``_core.data_integrate``) does the work; this module is the
thin shim over it plus a full pure-numpy reference (``test_cpp_matches_python``
pins the two lanes byte-identical) -- unlike ``conservative_interpolate``/
``subdivide``/``agglomerate``, this operation has no discrete winding/sign
branch, so a second implementation cannot silently disagree with the first.

Public API:
    data_integrate
"""

from __future__ import annotations

import numpy as np

from ._data_average import _cell_measures

__all__ = ["data_integrate"]


def _block_bases(mesh):
    bases = [0]
    for block in mesh.cells:
        bases.append(bases[-1] + len(block.data))
    return bases


def _global_to_block_row(bases, g):
    for b in range(len(bases) - 1):
        if bases[b] <= g < bases[b + 1]:
            return b, g - bases[b]
    return None


def _block_measures(mesh):
    """One |measure| array per cell block; unmeasurable cells (including a
    whole unmeasurable block) read as the sentinel -1.0, matching the C++
    core's convention exactly."""
    out = []
    for block in mesh.cells:
        n = len(block.data)
        w = _cell_measures(mesh, block)
        if w is None:
            w = np.full(n, -1.0)
        else:
            w = np.where(np.isfinite(w) & (w > 0), np.asarray(w, dtype=float), -1.0)
        out.append(w)
    return out


def _region_cells(mesh, bases):
    """Named Cell regions, resolved to a sorted list of (block, row) pairs,
    in ``sorted(name)`` order -- mirroring ``Mesh::RegionNames()``. Where
    several Cell regions share a name (distinguished by dim/tag), the one
    with the smallest ``(dim, tag)`` wins, matching the C++ core's
    canonical ``(kind, name, dim, tag)`` storage order that
    ``Mesh::FindRegion`` resolves against."""
    cell_regions = [r for r in mesh.regions if r.kind == "cell"]
    names = sorted({r.name for r in cell_regions})
    out = []
    for name in names:
        candidates = [r for r in cell_regions if r.name == name]
        candidates.sort(key=lambda r: (r.dim, r.tag))
        region = candidates[0]
        cells = []
        for g in np.asarray(region.entries, dtype=np.int64).ravel():
            br = _global_to_block_row(bases, int(g))
            if br is not None:
                cells.append(br)
        out.append((name, cells))
    return out


def _weighted_reduce(values, weights):
    """(total_per_component, domain_measure_per_component, mean_per_component,
    num_nan_per_component) over the rows named by ``values``/``weights``
    (already restricted to the cells of interest, weights already carrying
    the -1.0 unmeasurable sentinel)."""
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    ncomp = 1 if values.ndim < 2 else int(np.prod(values.shape[1:]))
    flat = values.reshape(-1, ncomp)
    measurable = weights > 0.0
    finite = np.isfinite(flat) & measurable[:, None]
    w_col = np.where(finite, weights[:, None], 0.0)
    total = np.where(finite, flat * w_col, 0.0).sum(axis=0)
    domain_measure = w_col.sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = np.where(
            domain_measure > 0.0,
            total / np.where(domain_measure > 0.0, domain_measure, 1.0),
            np.nan,
        )
    num_nan = (measurable[:, None] & ~np.isfinite(flat)).sum(axis=0)
    return total, domain_measure, mean, num_nan


def _region_dict(name, num_cells, num_skipped, total, domain_measure, mean, num_nan):
    return {
        "name": name,
        "num_cells": int(num_cells),
        "num_skipped": int(num_skipped),
        "domain_measure_per_component": [float(v) for v in domain_measure],
        "total_per_component": [float(v) for v in total],
        "mean_per_component": [float(v) for v in mean],
        "num_nan_per_component": [int(v) for v in num_nan],
    }


def _data_integrate_py(mesh, names):
    """Pure-numpy reference implementation of :func:`data_integrate`."""
    selected = list(names) if names else sorted(mesh.cell_data.keys())
    for name in selected:
        if name not in mesh.cell_data:
            if name in mesh.point_data:
                raise ValueError(
                    f"meshio++: data_integrate: '{name}' is a point_data array; "
                    "convert it first with point_data_to_cell_data (CLI: `data to-cell`)"
                )
            raise ValueError(
                f"meshio++: data_integrate: no cell_data array named {name!r} "
                f"(available: {', '.join(sorted(mesh.cell_data)) or 'none'})"
            )
    for name in selected:
        blocks = mesh.cell_data[name]
        if len(blocks) != len(mesh.cells):
            raise ValueError(
                f"meshio++: data_integrate: cell_data {name!r} does not have "
                "one array per cell block"
            )

    out = []
    if not selected:
        return out

    measures = _block_measures(mesh)
    domain_num_cells = int(sum(int((m > 0).sum()) for m in measures))
    domain_num_skipped = int(sum(int((m <= 0).sum()) for m in measures))

    bases = _block_bases(mesh)
    region_defs = _region_cells(mesh, bases)
    region_num_cells = []
    region_num_skipped = []
    for _, cells in region_defs:
        nc = sum(1 for b, row in cells if measures[b][row] > 0.0)
        region_num_cells.append(nc)
        region_num_skipped.append(len(cells) - nc)

    for name in selected:
        blocks = mesh.cell_data[name]
        first = np.asarray(blocks[0], dtype=float)
        ncomp = 1 if first.ndim < 2 else int(np.prod(first.shape[1:]))

        totals = np.zeros(ncomp)
        dmeasures = np.zeros(ncomp)
        nnans = np.zeros(ncomp, dtype=np.int64)
        for b, block_vals in enumerate(blocks):
            t, dm, _, nn = _weighted_reduce(block_vals, measures[b])
            totals += t
            dmeasures += dm
            nnans += nn
        with np.errstate(invalid="ignore", divide="ignore"):
            means = np.where(
                dmeasures > 0.0,
                totals / np.where(dmeasures > 0.0, dmeasures, 1.0),
                np.nan,
            )
        domain = _region_dict(
            "", domain_num_cells, domain_num_skipped, totals, dmeasures, means, nnans
        )

        regions = []
        for ridx, (rname, cells) in enumerate(region_defs):
            if not cells:
                r_totals = np.zeros(ncomp)
                r_dmeasures = np.zeros(ncomp)
                r_means = np.full(ncomp, np.nan)
                r_nnans = np.zeros(ncomp, dtype=np.int64)
            else:
                vals = np.stack(
                    [
                        np.asarray(blocks[b], dtype=float).reshape(-1, ncomp)[row]
                        for b, row in cells
                    ]
                )
                w = np.array([measures[b][row] for b, row in cells])
                r_totals, r_dmeasures, r_means, r_nnans = _weighted_reduce(vals, w)
            regions.append(
                _region_dict(
                    rname,
                    region_num_cells[ridx],
                    region_num_skipped[ridx],
                    r_totals,
                    r_dmeasures,
                    r_means,
                    r_nnans,
                )
            )

        out.append(
            {
                "name": name,
                "num_components": ncomp,
                "domain": domain,
                "regions": regions,
            }
        )
    return out


def data_integrate(mesh, arrays=None):
    """Cell-measure-weighted total/mean of one or more ``cell_data`` arrays,
    for the whole mesh and independently for every named ``Cell`` region.

    :param mesh: the mesh to integrate over (never modified).
    :param arrays: ``cell_data`` array names to integrate; ``None`` means
        every ``cell_data`` array, in sorted name order.
    :returns: a list of dicts, one per array, each shaped
        ``{"name", "num_components", "domain": {...}, "regions": [{...}, ...]}``
        where each region entry (including ``"domain"``, whose ``"name"`` is
        ``""``) is
        ``{"name", "num_cells", "num_skipped", "domain_measure_per_component",
        "total_per_component", "mean_per_component", "num_nan_per_component"}``.
    :raises ValueError: on an unknown array name, a ``point_data``-only name
        (naming ``point_data_to_cell_data`` as the fix), or a ``cell_data``
        array whose block count disagrees with the mesh.
    """
    names = list(arrays) if arrays else []
    try:
        from . import _core

        return _core.data_integrate(mesh, names)
    except (ValueError, TypeError):
        # A genuine user error must not fall through to the numpy path, which
        # would either raise something less helpful or silently succeed
        # differently.
        raise
    except Exception:
        return _data_integrate_py(mesh, names)

"""Named regions in a VTK XML file's ``<FieldData>`` (``.vtu``, ``.vtp``): the
twin of the core's ``detail/region_field_data.hpp`` (doc/regions.md).

Each region is one Int64 field array -- ``region:point:<name>`` (point
indices), ``region:cell:<name>`` (cell indices), ``region:side:<name>``
(``(cell, local facet)`` pairs, two components) -- plus
``region-meta:<kind>:<name>`` = ``[dim, tag]`` when either is not -1. Cell
indices are in the file's cell order. A reader removes the arrays it
understands from ``field_data``; a malformed one is warned about and kept.
"""

from __future__ import annotations

import numpy as np

from ._common import warn
from ._regions import Region

PREFIX = "region:"
META_PREFIX = "region-meta:"
_KIND_ORDER = {"point": 0, "cell": 1, "side": 2}

# A hidden cell_data array a reader threads through its cell reconstruction:
# whatever reordering that applies to cells, it applies to this too, so it
# comes out as each built cell's file index.
FILE_INDEX_KEY = "\x00meshio++:file_index"


def is_region_field_name(name) -> bool:
    return str(name).startswith(PREFIX) or str(name).startswith(META_PREFIX)


def _parse(rest):
    kind, sep, name = rest.partition(":")
    if not sep or kind not in _KIND_ORDER:
        return None
    return kind, name


def regions_to_field_arrays(mesh, global_to_file=None) -> dict:
    """``{name: Int64 array}`` for ``mesh``'s regions, in canonical order."""
    out = {}
    regions = sorted(
        getattr(mesh, "regions", []),
        key=lambda r: (_KIND_ORDER[r.kind], r.name, r.dim, r.tag),
    )
    g2f = None if global_to_file is None else np.asarray(global_to_file, dtype=np.int64)
    for r in regions:
        e = np.asarray(r.entries, dtype=np.int64)
        if r.kind == "cell" and g2f is not None:
            e = g2f[e] if e.size else e
        elif r.kind == "side":
            e = e.reshape(-1, 2).copy()
            if g2f is not None and len(e):
                e[:, 0] = g2f[e[:, 0]]
        key = f"{r.kind}:{r.name}"
        out[PREFIX + key] = e.astype(np.int64)
        if r.dim != -1 or r.tag != -1:
            out[META_PREFIX + key] = np.array([r.dim, r.tag], dtype=np.int64)
    return out


def file_to_global_from(cell_data):
    """Pop :data:`FILE_INDEX_KEY` from ``cell_data`` and invert it."""
    parts = cell_data.pop(FILE_INDEX_KEY, None)
    if parts is None:
        return None
    file_of_global = (
        np.concatenate([np.asarray(p).reshape(-1) for p in parts])
        if len(parts)
        else np.empty(0, dtype=np.int64)
    ).astype(np.int64)
    f2g = np.full(len(file_of_global), -1, dtype=np.int64)
    f2g[file_of_global] = np.arange(len(file_of_global), dtype=np.int64)
    return f2g


def regions_from_field_arrays(field_data, num_points, file_to_global, fmt):
    """Split ``field_data`` into ``(regions, remaining field_data)``."""
    metas = {}
    for name, arr in field_data.items():
        if name.startswith(META_PREFIX):
            a = np.asarray(arr)
            if np.issubdtype(a.dtype, np.integer) and a.size == 2:
                metas[name[len(META_PREFIX) :]] = (
                    int(a.reshape(-1)[0]),
                    int(a.reshape(-1)[1]),
                )
    f2g = None if file_to_global is None else np.asarray(file_to_global, dtype=np.int64)
    regions, rest = [], {}
    for name, arr in field_data.items():
        if name.startswith(META_PREFIX):
            if name[len(META_PREFIX) :] not in metas:
                warn(f"{fmt}: malformed region meta array '{name}' kept as field data")
                rest[name] = arr
            continue
        parsed = _parse(name[len(PREFIX) :]) if name.startswith(PREFIX) else None
        if parsed is None:
            rest[name] = arr
            continue
        kind, rname = parsed
        a = np.asarray(arr)
        stride = 2 if kind == "side" else 1
        comps = a.shape[1] if a.ndim >= 2 else 1
        ok = np.issubdtype(a.dtype, np.integer) and comps == stride
        e = a.astype(np.int64).reshape(-1, stride) if ok else None
        if ok and len(e):
            if kind == "point":
                ok = bool(np.all((e >= 0) & (e < num_points)))
            else:
                cells = e[:, 0]
                n = len(f2g) if f2g is not None else None
                ok = bool(np.all(cells >= 0)) and (n is None or bool(np.all(cells < n)))
                if ok and f2g is not None:
                    e = e.copy()
                    e[:, 0] = f2g[cells]
                    ok = bool(np.all(e[:, 0] >= 0))
                if ok and stride == 2:
                    ok = bool(np.all(e[:, 1] >= 0))
        if not ok:
            warn(
                f"{fmt}: '{name}' is not a well-formed region array; kept as field data"
            )
            rest[name] = arr
            continue
        dim, tag = metas.get(name[len(PREFIX) :], (-1, -1))
        entries = e.reshape(-1) if stride == 1 else e
        regions.append(Region(rname, kind, entries, dim, tag))
    return regions, rest

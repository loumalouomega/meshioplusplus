"""Move data between the point and cell locations by averaging.

A dependency-free mesh *operation* (not a file format). The geometry is never
modified; only data arrays are added.

``point_data_to_cell_data`` gives each cell the mean of the values at its own
nodes. ``cell_data_to_point_data`` gives each point the mean of the values on
the cells incident to it, optionally weighted by each cell's |measure| (area for
a 2D cell, volume for a 3D one).

Both work component-wise, so scalar, vector and tensor arrays behave the same.
Because a mean is not an integer, the output is **always** ``float64``.

Non-finite values never contribute to an average; a point or cell with no finite
contribution yields NaN. See :mod:`meshioplusplus._data_common` for the policy.

Public API:
    point_data_to_cell_data
    cell_data_to_point_data
"""

from __future__ import annotations

import copy

import numpy as np

from ._data_common import normalize_location, require_key

#: Outward-wound corner-only face tables, transcribed verbatim (winding
#: included) from `detail/cell_faces.cpp`'s C++ tables -- the authoritative,
#: Newell-normal-gtested source. This table's own consistency IS the
#: correctness of `_cell_measures`' divergence-theorem volume below: the
#: origin-tetrahedra decomposition it uses only gives a translation-invariant
#: answer when every face of the closed surface is wound outward the same
#: way, so a single inward face here would silently make the "volume" of a
#: cell depend on where the cell sits in space rather than on its shape --
#: exactly the defect `TranslatedCellMeasureIsInvariant` below pins.
_FACES = {
    "tetra": [(0, 1, 3), (1, 2, 3), (2, 0, 3), (0, 2, 1)],
    "hexahedron": [
        (0, 4, 7, 3),
        (1, 2, 6, 5),
        (0, 1, 5, 4),
        (3, 7, 6, 2),
        (0, 3, 2, 1),
        (4, 5, 6, 7),
    ],
    "wedge": [(0, 2, 1), (3, 4, 5), (0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)],
    "pyramid": [(0, 3, 2, 1), (0, 1, 4), (1, 2, 4), (2, 3, 4), (3, 0, 4)],
}

_CORNERS = {
    "line": 2,
    "triangle": 3,
    "triangle6": 3,
    "quad": 4,
    "quad8": 4,
    "quad9": 4,
    "tetra": 4,
    "tetra10": 4,
    "hexahedron": 8,
    "hexahedron20": 8,
    "hexahedron27": 8,
    "wedge": 6,
    "wedge15": 6,
    "pyramid": 5,
    "pyramid13": 5,
}

_DIM = {
    "line": 1,
    "triangle": 2,
    "triangle6": 2,
    "quad": 2,
    "quad8": 2,
    "quad9": 2,
    "tetra": 3,
    "tetra10": 3,
    "hexahedron": 3,
    "hexahedron20": 3,
    "hexahedron27": 3,
    "wedge": 3,
    "wedge15": 3,
    "pyramid": 3,
    "pyramid13": 3,
}


def _pad3(points):
    """Points padded to three columns, so 2D meshes work unchanged."""
    p = np.asarray(points, dtype=float)
    if p.shape[1] >= 3:
        return p[:, :3]
    return np.hstack([p, np.zeros((p.shape[0], 3 - p.shape[1]))])


def _cell_measures(mesh, block):
    """|measure| of every cell in ``block``, or ones where undefined."""
    ctype = block.type.split("_")[0] if "_" in block.type else block.type
    data = block.data
    if isinstance(data, list) or ctype.startswith(("polygon", "polyhedron")):
        return None
    conn = np.asarray(data)
    corners = _CORNERS.get(ctype)
    dim = _DIM.get(ctype)
    if corners is None or dim is None:
        return None
    pts = _pad3(mesh.points)[conn[:, :corners]]
    if dim == 1:
        return np.linalg.norm(pts[:, 1] - pts[:, 0], axis=1)
    if dim == 2:
        # Newell normal, giving the unsigned polygon area.
        s = np.zeros((pts.shape[0], 3))
        for i in range(corners):
            s += np.cross(pts[:, i], pts[:, (i + 1) % corners])
        return 0.5 * np.linalg.norm(s, axis=1)
    faces = _FACES.get(ctype)
    if faces is None:
        return None
    vol6 = np.zeros(pts.shape[0])
    for face in faces:
        a = pts[:, face[0]]
        for i in range(1, len(face) - 1):
            vol6 += np.einsum(
                "ij,ij->i", a, np.cross(pts[:, face[i]], pts[:, face[i + 1]])
            )
    return np.abs(vol6 / 6.0)


def _cell_node_lists(block):
    """The node ids of each cell, covering ragged and polyhedron blocks."""
    data = block.data
    if block.type.startswith("polyhedron"):
        out = []
        for cell in data:
            seen = []
            marker = set()
            for face in cell:
                for n in face:
                    if n not in marker:
                        marker.add(n)
                        seen.append(int(n))
            out.append(np.asarray(seen, dtype=np.int64))
        return out
    if isinstance(data, list):
        return [np.asarray(row, dtype=np.int64) for row in data]
    conn = np.asarray(data)
    return [conn[i].astype(np.int64) for i in range(conn.shape[0])]


def _apply_nan_policy(values, nan_policy, nan_replacement, name):
    """Apply the shared NaN policy to a produced array."""
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
            f"meshio++: data_average: non-finite value at index {idx} of '{name}'"
        )
    raise ValueError(
        f"meshio++: unknown NaN policy '{nan_policy}' "
        "(expected 'ignore', 'replace' or 'fail')"
    )


def _to_cell_py(mesh, names, prefix, suffix, overwrite, nan_policy, nan_replacement):
    """Pure-Python reference for :func:`point_data_to_cell_data`."""
    out = copy.deepcopy(mesh)
    selected = list(names) if names else sorted(mesh.point_data.keys())
    for name in selected:
        require_key(mesh, "point_data", name)
    npoints = len(mesh.points)
    for name in selected:
        target = f"{prefix}{name}{suffix}"
        if not overwrite and target in out.cell_data:
            raise ValueError(
                f"meshio++: data_average: cell_data '{target}' already exists "
                "(pass overwrite=True to replace it)"
            )
        src = np.asarray(mesh.point_data[name], dtype=float)
        ncomp = 1 if src.ndim < 2 else int(np.prod(src.shape[1:]))
        flat = src.reshape(npoints, ncomp)
        blocks = []
        for block in mesh.cells:
            nodes = _cell_node_lists(block)
            res = np.full((len(nodes), ncomp), np.nan)
            for c, ids in enumerate(nodes):
                ids = ids[(ids >= 0) & (ids < npoints)]
                if ids.size == 0:
                    continue
                vals = flat[ids]
                mask = np.isfinite(vals)
                counts = mask.sum(axis=0)
                sums = np.where(mask, vals, 0.0).sum(axis=0)
                with np.errstate(invalid="ignore", divide="ignore"):
                    res[c] = np.where(counts > 0, sums / counts, np.nan)
            res = _apply_nan_policy(res, nan_policy, nan_replacement, name)
            blocks.append(res.reshape(len(nodes)) if ncomp == 1 else res)
        out.cell_data[target] = blocks
    return out


def _to_point_py(
    mesh, names, weight, prefix, suffix, overwrite, nan_policy, nan_replacement
):
    """Pure-Python reference for :func:`cell_data_to_point_data`."""
    out = copy.deepcopy(mesh)
    selected = list(names) if names else sorted(mesh.cell_data.keys())
    for name in selected:
        require_key(mesh, "cell_data", name)
    npoints = len(mesh.points)

    weights = []
    for block in mesh.cells:
        n = len(block.data)
        if weight == "measure":
            w = _cell_measures(mesh, block)
            if w is None:
                w = np.ones(n)
            w = np.where(np.isfinite(w) & (w > 0), w, 1.0)
        else:
            w = np.ones(n)
        weights.append(w)

    for name in selected:
        target = f"{prefix}{name}{suffix}"
        if not overwrite and target in out.point_data:
            raise ValueError(
                f"meshio++: data_average: point_data '{target}' already exists "
                "(pass overwrite=True to replace it)"
            )
        blocks = mesh.cell_data[name]
        if len(blocks) != len(mesh.cells):
            raise ValueError(
                f"meshio++: data_average: cell_data '{name}' has {len(blocks)} "
                f"block(s) but the mesh has {len(mesh.cells)} cell block(s)"
            )
        first = np.asarray(blocks[0], dtype=float)
        ncomp = 1 if first.ndim < 2 else int(np.prod(first.shape[1:]))
        acc = np.zeros((npoints, ncomp))
        wacc = np.zeros((npoints, ncomp))
        for b, block in enumerate(mesh.cells):
            vals = np.asarray(blocks[b], dtype=float).reshape(-1, ncomp)
            nodes = _cell_node_lists(block)
            for c, ids in enumerate(nodes):
                ids = ids[(ids >= 0) & (ids < npoints)]
                if ids.size == 0:
                    continue
                v = vals[c]
                finite = np.isfinite(v)
                if not finite.any():
                    continue
                w = weights[b][c]
                contrib = np.where(finite, w * v, 0.0)
                wcontrib = np.where(finite, w, 0.0)
                np.add.at(acc, ids, contrib)
                np.add.at(wacc, ids, wcontrib)
        with np.errstate(invalid="ignore", divide="ignore"):
            res = np.where(wacc > 0, acc / np.where(wacc > 0, wacc, 1.0), np.nan)
        res = _apply_nan_policy(res, nan_policy, nan_replacement, name)
        out.point_data[target] = res.reshape(npoints) if ncomp == 1 else res
    return out


def point_data_to_cell_data(
    mesh,
    keys=None,
    prefix: str = "",
    suffix: str = "",
    overwrite: bool = True,
    nan_policy: str = "ignore",
    nan_replacement: float = 0.0,
):
    """Average ``point_data`` onto the cells.

    Each cell's value is the mean over its own nodes.

    Args:
        mesh: the source mesh (unmodified).
        keys: names to convert; ``None`` means every ``point_data`` array.
        prefix: output name is ``prefix + name + suffix``.
        suffix: see ``prefix``.
        overwrite: whether an existing target name may be replaced.
        nan_policy: ``"ignore"``, ``"replace"`` or ``"fail"``.
        nan_replacement: value used when ``nan_policy`` is ``"replace"``.

    Returns:
        A new mesh carrying the produced ``cell_data``.
    """
    names = list(keys) if keys else []
    try:
        from . import _core
    except Exception:
        return _to_cell_py(
            mesh, names, prefix, suffix, overwrite, nan_policy, nan_replacement
        )
    try:
        out = _core.point_data_to_cell_data(
            mesh, names, prefix, suffix, overwrite, nan_policy, nan_replacement
        )
    except ValueError:
        raise
    except Exception:
        return _to_cell_py(
            mesh, names, prefix, suffix, overwrite, nan_policy, nan_replacement
        )
    _carry_sets(mesh, out)
    return out


def cell_data_to_point_data(
    mesh,
    keys=None,
    weighted: bool = False,
    prefix: str = "",
    suffix: str = "",
    overwrite: bool = True,
    nan_policy: str = "ignore",
    nan_replacement: float = 0.0,
):
    """Average ``cell_data`` onto the points.

    Each point's value is the mean over the cells incident to it, optionally
    weighted by each cell's |measure|.

    Args:
        mesh: the source mesh (unmodified).
        keys: names to convert; ``None`` means every ``cell_data`` array.
        weighted: weight by cell measure instead of counting cells equally.
        prefix: output name is ``prefix + name + suffix``.
        suffix: see ``prefix``.
        overwrite: whether an existing target name may be replaced.
        nan_policy: ``"ignore"``, ``"replace"`` or ``"fail"``.
        nan_replacement: value used when ``nan_policy`` is ``"replace"``.

    Returns:
        A new mesh carrying the produced ``point_data``.
    """
    names = list(keys) if keys else []
    weight = "measure" if weighted else "uniform"
    try:
        from . import _core
    except Exception:
        return _to_point_py(
            mesh, names, weight, prefix, suffix, overwrite, nan_policy, nan_replacement
        )
    try:
        out = _core.cell_data_to_point_data(
            mesh, names, weight, prefix, suffix, overwrite, nan_policy, nan_replacement
        )
    except ValueError:
        raise
    except Exception:
        return _to_point_py(
            mesh, names, weight, prefix, suffix, overwrite, nan_policy, nan_replacement
        )
    _carry_sets(mesh, out)
    return out


def _carry_sets(source, target) -> None:
    """Copy sets across; they index geometry, which is never modified here."""
    for attr in ("point_sets", "cell_sets"):
        value = getattr(source, attr, None)
        if value:
            setattr(target, attr, copy.deepcopy(value))


__all__ = ["point_data_to_cell_data", "cell_data_to_point_data", "normalize_location"]

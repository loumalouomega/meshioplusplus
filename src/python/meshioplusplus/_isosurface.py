"""Isosurfaces / contours: the level set of a scalar ``point_data`` field.

A dependency-free mesh *operation* (not a file format): it returns a new mesh one
topological dimension below the cut cells -- a 3D volume mesh yields a
``triangle``/``quad`` surface, a 2D surface mesh a ``line`` contour.

This is the **data-driven sibling of** :func:`meshioplusplus.slice`: slice cuts
where ``dot(x - origin, normal) = 0``, isosurface where ``f(x) - isovalue = 0``.
Both go through the same cutter (:mod:`._marching`, the numpy twin of the C++
``detail/marching.hpp``), so the watertightness, winding, degeneracy and
determinism guarantees are literally the same code.

Rules worth knowing before use (all documented at length in ``doc/isosurface.md``):

* **The field must be ``point_data``.** ``cell_data`` is piecewise constant, so
  there is no crossing to locate and no level set to draw; naming one raises,
  pointing at :func:`meshioplusplus.cell_data_to_point_data` (CLI:
  ``meshioplusplus data to-point``) as the fix.
* **Several isovalues land in one mesh**, cut in ascending order (sorted, exact
  duplicates dropped) and tagged with Float64 ``iso:value`` and Int64
  ``iso:index``. Splitting into one mesh per contour is
  ``split(mesh, by="region", tag="iso:index")`` -- the tag criterion needs an
  integer array, which is what ``iso:index`` is for.
* **The contoured field reads back as exactly the isovalue** on the cut points
  (every other array is interpolated at the crossing). The exception is a
  multi-component array reduced by *magnitude*: ``|lerp(v)| != lerp(|v|)``
  mathematically, so that case stays approximate.
* A node whose value is exactly the isovalue counts as being on the **positive**
  side, so a plateau at exactly the isovalue emits its boundary once, not twice.

The C++ core (``_core.isosurface``) does the work; this module is the thin shim
(try C++, fall back to a pure-numpy reference) plus the reference itself, which
reproduces the C++ result **byte for byte**, pinned by
``tests/python/test_isosurface.py::test_cpp_matches_python``.

Public API:

* :func:`isosurface` -- the level set(s) of a scalar field.
"""

from __future__ import annotations

import numpy as np

from ._marching import FIELD_GRADIENT, _marching_cut, _marching_prepare
from ._mesh import Mesh

_TYPE_ORDER = ("triangle", "quad", "line")


def _components(arr):
    rows = arr.shape[0] if arr.ndim >= 1 else 0
    if rows == 0:
        return 0
    return int(arr.size // rows)


def _scalarize(flat, component):
    """Reduce rows to scalars: the requested component, or the magnitude.

    The magnitude sums left to right and takes one ``sqrt`` at the end -- the
    same arithmetic as the C++ ``iso_scalarize``, rather than
    ``np.linalg.norm``, so the two round alike.
    """
    ncomp = flat.shape[1]
    if component is not None:
        c = int(component)
        if c < 0 or c >= ncomp:
            raise ValueError(
                f"meshio++: isosurface: component {c} is out of range for an "
                f"array with {ncomp} component(s)"
            )
        return flat[:, c].astype(np.float64)
    if ncomp == 1:
        return flat[:, 0].astype(np.float64)
    out = np.zeros(flat.shape[0], dtype=np.float64)
    for k in range(ncomp):
        v = flat[:, k]
        out = out + v * v
    return np.sqrt(out)


def _type_block(mesh, cell_type):
    for i, cb in enumerate(mesh.cells):
        if cb.type == cell_type:
            return i
    return None


def _concat_contours(contours, dim, exact_name, exact_component):
    """Stack the contours into one mesh, merging cell blocks by type."""
    point_base = []
    total = 0
    for _value, _index, sec in contours:
        point_base.append(total)
        total += len(sec.points)

    points = np.zeros((total, dim), dtype=np.float64)
    for k, (_value, _index, sec) in enumerate(contours):
        n = len(sec.points)
        if n:
            points[point_base[k] : point_base[k] + n] = np.asarray(sec.points)

    point_data = {}
    for name in sorted(contours[0][2].point_data):
        a0 = np.asarray(contours[0][2].point_data[name])
        comp_shape = a0.shape[1:]
        ncomp = _components(a0) or 1
        out_arr = np.zeros((total, ncomp), dtype=np.float64)
        for k, (value, _index, sec) in enumerate(contours):
            a = np.asarray(sec.point_data[name])
            n = a.shape[0]
            if n:
                out_arr[point_base[k] : point_base[k] + n] = a.reshape(n, -1)
            # The contoured field reads back as *exactly* the isovalue (see the
            # module docstring); not applied to a magnitude reduction.
            if name == exact_name:
                sl = slice(point_base[k], point_base[k] + n)
                if exact_component is not None:
                    out_arr[sl, int(exact_component)] = value
                elif ncomp == 1:
                    out_arr[sl, 0] = value
        point_data[name] = out_arr.reshape((total,) + comp_shape)

    # Cell blocks merged by type, in the cutter's canonical order.
    out_types = []
    out_sources = []
    for cell_type in _TYPE_ORDER:
        src = []
        for k, (_value, _index, sec) in enumerate(contours):
            bi = _type_block(sec, cell_type)
            if bi is not None:
                src.append((k, bi))
        if src:
            out_types.append(cell_type)
            out_sources.append(src)

    cells = []
    for ob, cell_type in enumerate(out_types):
        parts = []
        for k, bi in out_sources[ob]:
            data = np.asarray(contours[k][2].cells[bi].data, dtype=np.int64)
            parts.append(data + point_base[k])
        cells.append((cell_type, np.concatenate(parts, axis=0)))

    out = Mesh(points, cells, point_data=point_data)

    cell_data = {}
    for name in sorted(contours[0][2].cell_data):
        blocks = []
        for ob in range(len(out_types)):
            parts = [
                np.asarray(contours[k][2].cell_data[name][bi])
                for k, bi in out_sources[ob]
            ]
            blocks.append(np.concatenate(parts, axis=0))
        cell_data[name] = blocks

    # iso:value is the real level; iso:index is its ordinal, and is what
    # split(by="region", tag=...) needs -- it only accepts an integer array.
    values = []
    indices = []
    for ob in range(len(out_types)):
        vparts = []
        iparts = []
        for k, bi in out_sources[ob]:
            n = len(np.asarray(contours[k][2].cells[bi].data))
            vparts.append(np.full(n, contours[k][0], dtype=np.float64))
            iparts.append(np.full(n, contours[k][1], dtype=np.int64))
        values.append(np.concatenate(vparts, axis=0))
        indices.append(np.concatenate(iparts, axis=0))
    cell_data["iso:value"] = values
    cell_data["iso:index"] = indices

    out.cell_data = cell_data
    return out


def _isosurface_py(mesh, array, isovalues, component, record_parent_ids):
    if not array:
        raise ValueError("meshio++: isosurface: an array name is required")
    isovalues = [float(v) for v in isovalues]
    if not isovalues:
        raise ValueError("meshio++: isosurface: at least one isovalue is required")
    for v in isovalues:
        if not np.isfinite(v):
            raise ValueError("meshio++: isosurface: isovalues must be finite")

    if array not in mesh.point_data:
        if array in mesh.cell_data:
            raise ValueError(
                f"meshio++: isosurface: '{array}' is cell_data, which is piecewise "
                "constant and has no level set; convert it first with "
                "cell_data_to_point_data (CLI: meshioplusplus data to-point)"
            )
        available = ", ".join(sorted(mesh.point_data))
        raise ValueError(
            f"meshio++: no point_data array named '{array}' ("
            + (f"available: {available}" if available else "the mesh has no point_data")
            + ")"
        )

    # Ascending, de-duplicated (a repeated isovalue would emit its contour twice).
    isovalues = sorted(set(isovalues))

    prep = _marching_prepare(mesh)
    simp = prep.simp

    a = np.asarray(simp.point_data[array])
    nrows = a.shape[0] if a.ndim >= 1 else 0
    ncomp = _components(a)
    if ncomp == 0:
        raise ValueError(f"meshio++: isosurface: '{array}' carries no values")
    field = _scalarize(a.reshape(nrows, -1).astype(np.float64), component)

    contours = []
    empty_like = None
    for k, iso in enumerate(isovalues):
        dist = field - iso
        section = _marching_cut(
            prep,
            dist,
            "iso:parent_cell",
            record_parent_ids,
            orientation=FIELD_GRADIENT,
        )
        if not section.cells:
            if empty_like is None:
                empty_like = section
            continue
        contours.append((iso, k, section))

    if not contours:
        return empty_like if empty_like is not None else Mesh(np.zeros((0, 3)), [])

    exact = component is not None or ncomp == 1
    return _concat_contours(contours, prep.dim, array if exact else None, component)


def isosurface(
    mesh,
    array: str,
    isovalues,
    component=None,
    record_parent_ids: bool = False,
):
    """Compute the isosurface(s) of a scalar ``point_data`` field.

    Parameters
    ----------
    mesh :
        the mesh to contour (unmodified).
    array :
        name of the ``point_data`` array to contour. A ``cell_data`` name raises
        -- convert it with :func:`meshioplusplus.cell_data_to_point_data` first.
    isovalues :
        one level value, or a sequence of them. They are cut in ascending order
        (sorted, exact duplicates dropped) and all land in the returned mesh.
    component :
        component of a multi-component array to contour; the row magnitude when
        ``None`` (in which case the exactly-the-isovalue guarantee does not
        apply -- see the module docstring).
    record_parent_ids :
        when true, attach an Int64 ``iso:parent_cell`` ``cell_data`` array
        recording, per contour cell, the global (block-major) index of the input
        cell it was cut from.

    Returns
    -------
    Mesh
        the contours: a ``triangle``/``quad`` surface for a volume mesh, a
        ``line`` mesh for a 2D surface mesh, or an empty mesh when no isovalue
        crosses the field. Every contour cell carries a Float64 ``iso:value``
        and an Int64 ``iso:index`` ``cell_data`` entry. Contour points are all
        new, so ``point_sets``/``cell_sets`` are **not** carried; each contour
        cell inherits its parent cell's ``cell_data``, and interpolated
        ``point_data`` is promoted to Float64.

    Examples
    --------
    >>> import meshioplusplus                                  # doctest: +SKIP
    >>> shell = meshioplusplus.isosurface(vol, "T", [300.0, 400.0])   # doctest: +SKIP

    See Also
    --------
    meshioplusplus.slice :
        the geometry-driven sibling -- the level set of the distance to a plane,
        through the same cutter.
    """
    array = str(array)
    if np.isscalar(isovalues) or isinstance(isovalues, (int, float, np.floating)):
        isovalues = [float(isovalues)]
    else:
        isovalues = [float(v) for v in np.asarray(isovalues, dtype=np.float64).ravel()]
    component = None if component is None else int(component)

    out = None
    try:
        from . import _core

        out = _core.isosurface(
            mesh,
            array,
            isovalues,
            -1 if component is None else component,
            bool(record_parent_ids),
        )
    except (ValueError, TypeError):
        raise
    except Exception:
        out = None
    if out is None:
        out = _isosurface_py(mesh, array, isovalues, component, record_parent_ids)
    return out

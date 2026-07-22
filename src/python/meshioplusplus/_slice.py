"""Planar cross-section of a mesh (the intersection with a plane).

A dependency-free mesh *operation* (not a file format): it returns a new mesh
one topological dimension below the cut cells -- a 3D volume mesh cut by a plane
yields a ``triangle``/``quad`` surface, a 2D surface mesh yields a ``line`` mesh.
Unlike :func:`crop` (plane mode), which keeps whole cells on one side, ``slice``
computes the actual intersection and lowers the dimension.

The C++ core (``_core.slice``) does the work; this module is the thin shim
(try C++, fall back to a pure-numpy reference) plus the plane arithmetic. The
cutter itself -- marching tetrahedra on a simplexified input, the watertight
edge dedup, the winding and the degeneracy rule -- lives in :mod:`._marching`,
shared with :mod:`._isosurface` exactly as the C++ side shares
``detail/marching.hpp``; it reproduces the C++ result **byte for byte**, pinned
by ``tests/python/test_slice.py::test_cpp_matches_python``.

``slice`` shadows the Python built-in **only as a module attribute** -- the name
``meshioplusplus.slice`` is deliberate; the built-in is never rebound here.

Public API:

* :func:`slice` -- the planar cross-section of a mesh.
"""

from __future__ import annotations

import numpy as np

from ._marching import FIXED_DIRECTION, _marching_cut, _marching_prepare


def _slice_py(mesh, origin, normal, record_parent_ids):
    origin = np.asarray(origin, dtype=np.float64).reshape(3)
    normal = np.asarray(normal, dtype=np.float64).reshape(3)
    if float(normal @ normal) <= 0.0:
        raise ValueError("meshio++: slice: normal must be non-zero")

    prep = _marching_prepare(mesh)

    # Signed distance of every simplexified node to the plane (padded z=0 in 2D).
    pts = prep.points
    dim = prep.dim
    n_in = len(pts)
    pts3 = np.zeros((n_in, 3), dtype=np.float64)
    if n_in:
        pts3[:, :dim] = pts
    dist = (
        (pts3[:, 0] - origin[0]) * normal[0]
        + (pts3[:, 1] - origin[1]) * normal[1]
        + (pts3[:, 2] - origin[2]) * normal[2]
    )

    return _marching_cut(
        prep,
        dist,
        "slice:parent_cell",
        record_parent_ids,
        orientation=FIXED_DIRECTION,
        direction=normal,
    )


def slice(
    mesh,
    origin=(0.0, 0.0, 0.0),
    normal=(0.0, 0.0, 1.0),
    record_parent_ids: bool = False,
):
    """Compute the planar cross-section of a mesh.

    Parameters
    ----------
    mesh :
        the mesh to cut (unmodified).
    origin :
        a point on the cutting plane (3 floats).
    normal :
        the plane normal (3 floats; non-zero, need not be unit length).
    record_parent_ids :
        when true, attach an Int64 ``slice:parent_cell`` ``cell_data`` array
        recording, per section cell, the global (block-major) index of the input
        cell it was cut from.

    Returns
    -------
    Mesh
        the cross-section: a ``triangle``/``quad`` surface for a volume mesh, a
        ``line`` mesh for a 2D surface mesh, or an empty mesh when the plane
        misses the geometry. Section points are all new, so
        ``point_sets``/``cell_sets`` are **not** carried; each section cell
        inherits its parent cell's ``cell_data``, and interpolated
        ``point_data`` is promoted to Float64.

    See Also
    --------
    meshioplusplus.isosurface :
        the data-driven sibling -- the level set of a scalar field rather than
        of the distance to a plane, through the same cutter.
    """
    origin = [float(v) for v in np.asarray(origin, dtype=np.float64).reshape(3)]
    normal = [float(v) for v in np.asarray(normal, dtype=np.float64).reshape(3)]

    out = None
    try:
        from . import _core

        out = _core.slice(mesh, origin, normal, bool(record_parent_ids))
    except (ValueError, TypeError):
        raise
    except Exception:
        out = None
    if out is None:
        out = _slice_py(mesh, origin, normal, record_parent_ids)
    return out

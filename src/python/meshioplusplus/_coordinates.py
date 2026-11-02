"""Cylindrical and spherical components about an axis (roadmap 7.6.1).

``convert_coordinates(mesh, "cylindrical", axis=..., arrays=["U"])`` records
each point's ``(r, theta, z)`` as ``point_data["coords:cylindrical"]`` and
rewrites the named vector and tensor arrays in the local basis of each point
(a cell's at its centroid), the components a rotating machine or a pressure
vessel is read in; ``inverse=True`` turns local components back into Cartesian
ones. The contract, conventions and singular points are
``operations/coordinates.hpp``'s; this module is its numpy twin, used when the
C++ core declines. See ``doc/coordinates.md``.
"""

from __future__ import annotations

import numpy as np

from ._fallback import core_op_declined

__all__ = ["convert_coordinates", "CYLINDRICAL_NAME", "SPHERICAL_NAME"]

CYLINDRICAL_NAME = "coords:cylindrical"
SPHERICAL_NAME = "coords:spherical"
_SYSTEMS = ("cylindrical", "spherical")


def _frame(axis, reference):
    a = np.asarray(axis, dtype=float).reshape(3)
    norm = np.linalg.norm(a)
    if not (norm > 0.0) or not np.isfinite(norm):
        raise ValueError(
            "meshio++: convert_coordinates: the axis must be a nonzero vector"
        )
    a = a / norm
    for cand in (
        np.asarray(reference, dtype=float).reshape(3),
        [1.0, 0, 0],
        [0, 1.0, 0],
    ):
        cand = np.asarray(cand, dtype=float)
        off = cand - np.dot(cand, a) * a
        if np.linalg.norm(off) > 1e-9 * np.linalg.norm(cand):
            e1 = off / np.linalg.norm(off)
            break
    return a, e1, np.cross(a, e1)


def _local(positions, system, a, e1, e2, origin):
    """Coordinates ``(m, 3)``, bases ``(m, 3, 3)`` (rows) and the singular mask."""
    d = positions - np.asarray(origin, dtype=float).reshape(3)
    z = d @ a
    rho = d - z[:, None] * a
    r_cyl = np.linalg.norm(rho, axis=1)
    scale = np.linalg.norm(d, axis=1)
    on_axis = ~(r_cyl > 1e-14 * np.maximum(scale, 1.0))
    safe = np.where(on_axis, 1.0, r_cyl)
    e_rho = np.where(on_axis[:, None], e1, rho / safe[:, None])
    angle = np.where(on_axis, 0.0, np.arctan2(rho @ e2, rho @ e1))
    e_phi = np.cross(a, e_rho)
    axis_rows = np.broadcast_to(a, e_rho.shape)
    if system == "cylindrical":
        coords = np.column_stack([r_cyl, angle, z])
        basis = np.stack([e_rho, e_phi, axis_rows], axis=1)
        return coords, basis, on_axis
    at_centre = ~(scale > 0.0)
    polar = np.where(at_centre, 0.0, np.arctan2(r_cyl, z))
    c, s = np.cos(polar)[:, None], np.sin(polar)[:, None]
    e_r = s * e_rho + c * a
    e_theta = c * e_rho - s * a
    coords = np.column_stack([scale, polar, angle])
    basis = np.stack([e_r, e_theta, e_phi], axis=1)
    return coords, basis, on_axis | at_centre


def _convert(values, name, basis, inverse):
    v = np.asarray(values, dtype=np.float64)
    rows = basis.shape[0]
    if v.ndim == 2 and v.shape[0] == rows and v.shape[1] == 3:
        R = np.swapaxes(basis, 1, 2) if inverse else basis
        return np.einsum("nij,nj->ni", R, v)
    if (v.ndim == 2 and v.shape[0] == rows and v.shape[1] == 9) or (
        v.ndim == 3 and v.shape[0] == rows and v.shape[1:] == (3, 3)
    ):
        T = v.reshape(rows, 3, 3)
        R = np.swapaxes(basis, 1, 2) if inverse else basis
        out = np.einsum("nik,nkl,njl->nij", R, T, R)
        return out.reshape(v.shape)
    raise ValueError(
        f"meshio++: convert_coordinates: '{name}' is neither a vector (n, 3) nor a "
        "tensor (n, 9) or (n, 3, 3) array"
    )


def _convert_py(mesh, system, axis, origin, reference, arrays, inverse, record):
    from ._mesh import Mesh

    a, e1, e2 = _frame(axis, reference)
    pts = np.zeros((len(mesh.points), 3))
    p = np.asarray(mesh.points, dtype=float)
    pts[:, : min(3, p.shape[1])] = p[:, :3]
    out = Mesh(
        mesh.points,
        mesh.cells,
        point_data=dict(mesh.point_data),
        cell_data={k: list(v) for k, v in mesh.cell_data.items()},
        field_data=dict(mesh.field_data),
        point_sets=dict(mesh.point_sets),
        cell_sets=dict(mesh.cell_sets),
    )
    if getattr(mesh, "regions", None):
        out.regions = list(mesh.regions)
    num_singular = 0
    coords, basis, singular = _local(pts, system, a, e1, e2, origin)
    if record:
        out.point_data[
            CYLINDRICAL_NAME if system == "cylindrical" else SPHERICAL_NAME
        ] = coords
        num_singular = int(singular.sum())
    num_arrays = 0
    centroid_bases = None
    for name in arrays:
        found = False
        if name in mesh.point_data:
            found = True
            out.point_data[name] = _convert(mesh.point_data[name], name, basis, inverse)
            num_arrays += 1
        if name in mesh.cell_data:
            found = True
            if centroid_bases is None:
                centroid_bases = []
                for block in mesh.cells:
                    cents = []
                    for row in block.data:
                        ids = np.unique(
                            np.concatenate([np.ravel(f) for f in row])
                            if block.type.startswith("polyhedron")
                            else np.ravel(row)
                        ).astype(np.int64)
                        cents.append(pts[ids].mean(axis=0) if len(ids) else np.zeros(3))
                    cents = np.asarray(cents, dtype=float).reshape(-1, 3)
                    centroid_bases.append(_local(cents, system, a, e1, e2, origin)[1])
            out.cell_data[name] = [
                _convert(values, name, b, inverse)
                for values, b in zip(mesh.cell_data[name], centroid_bases)
            ]
            num_arrays += 1
        if not found:
            raise ValueError(
                f"meshio++: convert_coordinates: no point or cell data array named '{name}'"
            )
    return out, {"num_arrays": num_arrays, "num_singular": num_singular}


def convert_coordinates(
    mesh,
    system: str = "cylindrical",
    *,
    axis=(0.0, 0.0, 1.0),
    origin=(0.0, 0.0, 0.0),
    reference=(1.0, 0.0, 0.0),
    arrays=(),
    inverse: bool = False,
    record_coordinates: bool = True,
    return_report: bool = False,
):
    """Positions, vector and tensor fields in cylindrical or spherical
    components about an axis.

    The points stay where they are. ``record_coordinates`` attaches
    ``coords:cylindrical`` (``r, theta, z``) or ``coords:spherical`` (``r,
    theta, phi``, ``theta`` the polar angle from the axis) point data; each
    array in ``arrays`` -- point or cell data, a vector ``(n, 3)`` or a tensor
    ``(n, 9)``/``(n, 3, 3)`` -- is rewritten in the local basis of its point
    (a cell's at its centroid), or with ``inverse=True`` from local components
    back to Cartesian ones. Angles are measured about ``axis`` (through
    ``origin``) from ``reference`` projected off the axis. A point on the axis
    takes ``reference`` as its radial direction and angle 0.

    :returns: the new mesh, or ``(mesh, {"num_arrays", "num_singular"})``.
    :raises ValueError: on an unknown system, a zero axis, or an array that is
        missing or not a 3-D vector or tensor.
    """
    if system not in _SYSTEMS:
        raise ValueError(
            f"meshio++: convert_coordinates: unknown system '{system}' (expected "
            "'cylindrical' or 'spherical')"
        )
    if isinstance(arrays, str):
        arrays = [arrays]
    arrays = [str(a) for a in arrays]
    try:
        from . import _core

        res = _core.convert_coordinates(
            mesh,
            system,
            [float(v) for v in axis],
            [float(v) for v in origin],
            [float(v) for v in reference],
            arrays,
            bool(inverse),
            bool(record_coordinates),
        )
        out = res.pop("mesh")
        report = res
    except Exception as exc:
        if not core_op_declined(exc, "convert_coordinates"):
            raise
        out, report = _convert_py(
            mesh, system, axis, origin, reference, arrays, inverse, record_coordinates
        )
    return (out, report) if return_report else out

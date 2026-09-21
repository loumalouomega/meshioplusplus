import numpy as np

from .. import _core
from .._exceptions import WriteError
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ._pcd import read as _py_read
from ._pcd import write as _py_write


def read(filename, drop_invalid=False):
    """Read a PCD file (C++ core for real file paths, Python fallback).

    ``drop_invalid=True`` drops the points whose x, y or z is not finite (an
    organised cloud's invalid returns) and, when any is dropped, the organisation.
    """
    if not is_buffer(filename, "r"):
        try:
            return _core.pcd_read(str(filename), drop_invalid)
        except Exception as exc:
            if not core_declined(exc, "pcd", "read", filename):
                raise
    return _py_read(filename, drop_invalid=drop_invalid)


def write(filename, mesh, binary=True, data=None, point_dtype="float32"):
    """Write a PCD file (C++ core for real file paths, Python fallback).

    ``data`` is ``"ascii"``, ``"binary"`` or ``"binary_compressed"``; when omitted,
    ``binary`` picks between the first two. ``point_dtype`` is the precision of x, y,
    z (and of normals, curvature and intensity): ``"float32"`` (the default, the only
    one PCL's typed ``PointXYZ`` loaders accept), ``"float64"``, or ``"keep"`` to
    follow the mesh's own points (what an in-place rewrite wants).
    """
    if data is None:
        data = "binary" if binary else "ascii"
    if point_dtype == "keep":
        point_dtype = "float64" if mesh.points.dtype == np.float64 else "float32"
    if point_dtype not in ("float32", "float64"):
        raise WriteError("PCD: point_dtype must be 'float32', 'float64' or 'keep'")
    if not is_buffer(filename, "w"):
        try:
            _core.pcd_write(str(filename), mesh, data, point_dtype == "float64")
            return
        except Exception as exc:
            if not core_declined(exc, "pcd", "write", filename):
                raise
    return _py_write(filename, mesh, data=data, point_dtype=point_dtype)


register_format("pcd", [".pcd"], read, {"pcd": write})

__all__ = ["read", "write"]

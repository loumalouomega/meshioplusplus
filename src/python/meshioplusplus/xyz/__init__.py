from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ._xyz import read as _py_read
from ._xyz import write as _py_write


def read(filename, columns=None, delimiter=None):
    """Read a headerless ASCII point cloud (C++ core for file paths, Python fallback).

    ``columns`` names the columns (``x y z nx ny nz r g b a``, any other name is a
    scalar, ``_`` skips one) and ``delimiter`` overrides the whitespace/comma/semicolon
    detection; both are inferred when omitted.
    """
    if not is_buffer(filename, "r"):
        try:
            return _core.xyz_read(str(filename), list(columns or []), delimiter or "")
        except Exception as exc:
            if not core_declined(exc, "xyz", "read", filename):
                raise
    return _py_read(filename, columns=columns, delimiter=delimiter)


def write(filename, mesh, float_fmt=None):
    """Write an ASCII point cloud (C++ core for file paths, Python fallback).

    ``float_fmt`` is a format spec such as ``".16e"``; by default float32 columns use
    ``.9g`` and the rest ``.17g``.
    """
    if not is_buffer(filename, "w"):
        try:
            _core.xyz_write(str(filename), mesh, float_fmt or "")
            return
        except Exception as exc:
            if not core_declined(exc, "xyz", "write", filename):
                raise
    return _py_write(filename, mesh, float_fmt=float_fmt)


register_format(
    "xyz", [".xyz", ".xyzn", ".xyzrgb", ".asc", ".pts", ".txt"], read, {"xyz": write}
)

__all__ = ["read", "write"]

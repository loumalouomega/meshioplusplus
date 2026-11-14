from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ._ansys import read as _py_read
from ._ansys import write as _py_write


def read(filename, data=None, sibling_data=True):
    """Read a Fluent mesh (``.msh``) or case (``.cas``) (C++ core, Python fallback).

    ``data`` names a Fluent ``.dat`` file whose ``300``/``2300``/``3300`` sections
    become cell data (a cell zone's values on its cells, a face zone's on its
    boundary faces, NaN elsewhere); a ``.cas`` reads the ``.dat`` beside it unless
    ``sibling_data=False``. Zone types are ``field_data["fluent:zone_type:<name>"]``.
    """
    if not is_buffer(filename, "r"):
        try:
            return _core.ansys_read(
                str(filename), "" if data is None else str(data), bool(sibling_data)
            )
        except Exception as exc:
            if not core_declined(exc, "ansys", "read", filename):
                raise
    return _py_read(filename, data=data, sibling_data=sibling_data)


def write(filename, mesh, binary=True):
    """Write an Ansys/Fluent .msh file (C++ core, Python fallback)."""
    if not is_buffer(filename, "w"):
        try:
            _core.ansys_write(str(filename), mesh, binary)
            return
        except Exception as exc:
            if not core_declined(exc, "ansys", "write", filename):
                raise
    return _py_write(filename, mesh, binary)


register_format("ansys", [".msh", ".cas"], read, {"ansys": write})

__all__ = ["read", "write"]

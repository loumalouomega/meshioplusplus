from .. import _core
from .._files import is_buffer
from .._helpers import register_format
from ._unv import read as _py_read
from ._unv import write as _py_write


def read(filename):
    """Read an I-DEAS Universal file (C++ core; Python fallback for groups/fields)."""
    if not is_buffer(filename, "r"):
        try:
            return _core.unv_read(str(filename))
        except Exception:
            pass
    return _py_read(filename)


def write(filename, mesh):
    """Write an I-DEAS Universal file (C++ core; Python fallback for groups)."""
    has_sets = bool(getattr(mesh, "point_sets", None)) or bool(
        getattr(mesh, "cell_sets", None)
    )
    if not has_sets and not is_buffer(filename, "w"):
        try:
            _core.unv_write(str(filename), mesh)
            return
        except Exception:
            pass
    return _py_write(filename, mesh)


register_format("unv", [".unv"], read, {"unv": write})

__all__ = ["read", "write"]

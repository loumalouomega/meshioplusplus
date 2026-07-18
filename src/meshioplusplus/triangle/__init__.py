from .. import _core
from .._files import is_buffer
from .._helpers import register_format
from ._triangle import read as _py_read
from ._triangle import write as _py_write


def read(filename):
    """Read Triangle files (.node/.ele pair or .poly; C++ core for real file
    paths, Python fallback)."""
    if not is_buffer(filename, "r"):
        try:
            return _core.triangle_read(str(filename))
        except Exception:
            pass
    return _py_read(filename)


def write(filename, mesh):
    """Write Triangle files (.node/.ele pair or .poly; C++ core for real file
    paths, Python fallback)."""
    if not is_buffer(filename, "w"):
        try:
            _core.triangle_write(str(filename), mesh)
            return
        except Exception:
            pass
    return _py_write(filename, mesh)


# .node/.ele are also claimed by tetgen (registered first): the dispatcher
# falls through to this format when tetgen rejects a 2D file.
register_format("triangle", [".node", ".ele", ".poly"], read, {"triangle": write})

__all__ = ["read", "write"]

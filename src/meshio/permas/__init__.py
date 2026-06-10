from .. import _core
from .._files import is_buffer
from .._helpers import register_format
from ._permas import read as _py_read
from ._permas import write as _py_write


def read(filename):
    """Read a PERMAS dat file (C++ core for plain text, Python fallback for .gz)."""
    if not is_buffer(filename, "r") and not str(filename).endswith(".gz"):
        try:
            return _core.permas_read(str(filename))
        except Exception:
            pass
    return _py_read(filename)


def write(filename, mesh):
    """Write a PERMAS dat file (C++ core for plain text, Python fallback for .gz)."""
    if not is_buffer(filename, "w") and not str(filename).endswith(".gz"):
        try:
            _core.permas_write(str(filename), mesh)
            return
        except Exception:
            pass
    return _py_write(filename, mesh)


register_format(
    "permas", [".post", ".post.gz", ".dato", ".dato.gz"], read, {"permas": write}
)

__all__ = ["read", "write"]

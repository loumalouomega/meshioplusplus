from .. import _core
from .._files import is_buffer
from .._helpers import register_format
from ._medit import read as _py_read
from ._medit import write as _py_write


def _is_binary_name(filename):
    return str(filename).endswith("b")


def read(filename):
    """Read a Medit file.

    Uses the C++ core for the ascii .mesh variant; the binary .meshb variant is
    handled by the reference Python reader.
    """
    if not is_buffer(filename, "r") and not _is_binary_name(filename):
        try:
            return _core.medit_read_ascii(str(filename))
        except Exception:
            pass
    return _py_read(filename)


def write(filename, mesh, float_fmt=".16e"):
    """Write a Medit file (C++ core for ascii .mesh; Python for binary .meshb)."""
    if (
        not is_buffer(filename, "w")
        and not _is_binary_name(filename)
        and float_fmt == ".16e"
    ):
        try:
            _core.medit_write_ascii(str(filename), mesh)
            return
        except Exception:
            pass
    return _py_write(filename, mesh, float_fmt=float_fmt)


register_format("medit", [".mesh", ".meshb"], read, {"medit": write})

__all__ = ["read", "write"]

from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ._medit import read as _py_read
from ._medit import write as _py_write


def _is_binary_name(filename):
    return str(filename).endswith("b")


def read(filename):
    """Read a Medit file: ``.mesh`` ascii, or ``.meshb`` binary (libMeshb GMF)
    when the name ends in ``b``. The C++ core reads both; the reference Python
    reader is the fallback."""
    if not is_buffer(filename, "r"):
        try:
            if _is_binary_name(filename):
                return _core.medit_read_binary(str(filename))
            return _core.medit_read_ascii(str(filename))
        except Exception as exc:
            if not core_declined(exc, "medit", "read", filename):
                raise
    return _py_read(filename)


def write(filename, mesh, float_fmt=".16e"):
    """Write a Medit file: ascii ``.mesh``, or binary ``.meshb`` (little-endian
    float64, version 3 or, for 64-bit indices, 4) when the name ends in ``b``."""
    if not is_buffer(filename, "w") and (
        _is_binary_name(filename) or float_fmt == ".16e"
    ):
        try:
            if _is_binary_name(filename):
                _core.medit_write_binary(str(filename), mesh)
            else:
                _core.medit_write_ascii(str(filename), mesh)
            return
        except Exception as exc:
            if not core_declined(exc, "medit", "write", filename):
                raise
    return _py_write(filename, mesh, float_fmt=float_fmt)


register_format("medit", [".mesh", ".meshb"], read, {"medit": write})

__all__ = ["read", "write"]

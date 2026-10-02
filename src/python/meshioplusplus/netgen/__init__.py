from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ._netgen import read as _py_read
from ._netgen import write as _py_write


def read(filename):
    """Read a Netgen .vol/.vol.gz file, including names and periodic tables."""
    if not is_buffer(filename, "r"):
        try:
            return _core.netgen_read(str(filename))
        except Exception as exc:
            if not core_declined(exc, "netgen", "read", filename):
                raise
    return _py_read(filename)


def write(filename, mesh, float_fmt=".16e"):
    """Write Netgen names and periodic tables; gzip requires native zlib."""
    if not is_buffer(filename, "w"):
        try:
            _core.netgen_write(str(filename), mesh, float_fmt)
            return
        except Exception as exc:
            if not core_declined(exc, "netgen", "write", filename):
                raise
    return _py_write(filename, mesh, float_fmt)


register_format("netgen", [".vol", ".vol.gz"], read, {"netgen": write})

__all__ = ["read", "write"]

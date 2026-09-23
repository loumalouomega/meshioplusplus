from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ..mphtxt._mphtxt import read_binary as _py_read
from ..mphtxt._mphtxt import write_binary as _py_write


def read(filename):
    """Read a COMSOL binary mesh (.mphbin) file: the content of a .mphtxt file,
    little-endian int32/float64, strings as one int32 per character.

    Uses the C++ core for real file paths, the Python reference otherwise.
    """
    if not is_buffer(filename, "rb"):
        try:
            return _core.mphbin_read(str(filename))
        except Exception as exc:
            if not core_declined(exc, "mphbin", "read", filename):
                raise
    return _py_read(filename)


def write(filename, mesh):
    """Write a COMSOL binary mesh (.mphbin) file; see ``mphtxt.write``.

    The format has no comment slot, so no provenance is written.
    """
    if not is_buffer(filename, "wb"):
        try:
            _core.mphbin_write(str(filename), mesh)
            return
        except Exception as exc:
            if not core_declined(exc, "mphbin", "write", filename):
                raise
    return _py_write(filename, mesh)


register_format("mphbin", [".mphbin"], read, {"mphbin": write})

__all__ = ["read", "write"]

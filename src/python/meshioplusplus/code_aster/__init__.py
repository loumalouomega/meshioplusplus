from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ._code_aster import read as _py_read
from ._code_aster import write as _py_write


def read(filename):
    """Read a Code_Aster native mesh (``.mail``).

    ``COOR_*D`` blocks become the points, each element block (``POI1`` … ``HEXA27``)
    a cell block, and ``GROUP_MA``/``GROUP_NO`` cell and point regions. Node order
    is Code_Aster's own, which is not MED's. The C++ core reads the whole format;
    the Python reader is the fallback for buffers and for anything the C++ path
    raises on.
    """
    if not is_buffer(filename, "r"):
        try:
            return _core.code_aster_read(str(filename))
        except Exception as exc:
            if not core_declined(exc, "code_aster", "read", filename):
                raise
    return _py_read(filename)


def write(filename, mesh):
    """Write a Code_Aster native mesh (``.mail``).

    Nodes and elements are named ``N1…`` and ``M1…``; point and cell regions
    become ``GROUP_NO``/``GROUP_MA`` blocks, with names sanitised to letters,
    digits and ``_`` and truncated to 24 characters. Every line stays within the
    80 columns Code_Aster reads. Side regions and data arrays are dropped with a
    warning.
    """
    if not is_buffer(filename, "w"):
        try:
            _core.code_aster_write(str(filename), mesh)
            return
        except Exception as exc:
            if not core_declined(exc, "code_aster", "write", filename):
                raise
    return _py_write(filename, mesh)


register_format("code_aster", [".mail"], read, {"code_aster": write})

__all__ = ["read", "write"]

from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ._radioss import read as _py_read


def read(filename):
    """Read an OpenRadioss / Radioss starter deck (``*_0000.rad``).

    Nodes and element blocks (``/BRICK``, ``/TETRA10``, ``/BRIC20``, ``/SHELL``,
    ``/SH3N``, ``/BEAM``, ...) make the mesh; every cell's part, property and
    material are the ``radioss:part``/``radioss:property``/``radioss:material``
    cell data. Parts, subsets, ``/GRxxx`` groups and ``/SURF/SEG`` surfaces
    become regions. ``#include`` files are inlined. The C++ core reads the whole
    format; the Python reader is the fallback for anything the C++ path raises
    on.
    """
    if is_buffer(filename, "r"):
        raise TypeError("Radioss decks are read from a path (they #include files)")
    try:
        return _core.radioss_read(str(filename))
    except Exception as exc:
        if not core_declined(exc, "radioss", "read", filename):
            raise
    return _py_read(filename)


register_format("radioss", [".rad"], read, {})

__all__ = ["read"]

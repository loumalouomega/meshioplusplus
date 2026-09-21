from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ._lsdyna import read as _py_read
from ._lsdyna import write as _py_write


def read(filename):
    """Read an LS-DYNA keyword deck (``.k`` / ``.key`` / ``.dyn``).

    ``*NODE`` and the ``*ELEMENT_*`` cards become points and cells (tetra, pyramid and
    wedge are read out of their degenerate hexahedra), ``*PART`` becomes a cell region
    (title as name, ``pid`` as tag) and ``*SET_*`` becomes point, cell and side regions.
    ``*INCLUDE`` files are followed. The C++ core reads the whole format; the Python
    reader is the fallback for buffers and for anything the C++ path raises on.
    """
    if not is_buffer(filename, "r"):
        try:
            return _core.lsdyna_read(str(filename))
        except Exception as exc:
            if not core_declined(exc, "lsdyna", "read", filename):
                raise
    return _py_read(filename)


def write(filename, mesh):
    """Write an LS-DYNA keyword deck.

    Cell regions with a dimension become ``*PART`` cards (their tag is the ``pid``),
    every other region becomes a ``*SET_*_LIST``, and every part gets placeholder
    section and material ids: enough for a mesher to LS-DYNA hand-off.
    """
    if not is_buffer(filename, "w"):
        try:
            _core.lsdyna_write(str(filename), mesh)
            return
        except Exception as exc:
            if not core_declined(exc, "lsdyna", "write", filename):
                raise
    return _py_write(filename, mesh)


register_format("lsdyna", [".k", ".key", ".dyn"], read, {"lsdyna": write})

__all__ = ["read", "write"]

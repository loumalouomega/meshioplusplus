from .. import _core  # noqa: F401 (the shims look it up at call time)
from .._helpers import register_format
from .._shims import plain_shims
from ._lsdyna import read as _py_read  # noqa: F401
from ._lsdyna import write as _py_write  # noqa: F401

read, write = plain_shims(
    __name__,
    "lsdyna",
    """Read an LS-DYNA keyword deck (``.k`` / ``.key`` / ``.dyn``).

    ``*NODE`` and the ``*ELEMENT_*`` cards become points and cells (tetra, pyramid and
    wedge are read out of their degenerate hexahedra), ``*PART`` becomes a cell region
    (title as name, ``pid`` as tag) and ``*SET_*`` becomes point, cell and side regions.
    ``*INCLUDE`` files are followed. The C++ core reads the whole format; the Python
    reader is the fallback for buffers and for anything the C++ path raises on.
    """,
    """Write an LS-DYNA keyword deck.

    Cell regions with a dimension become ``*PART`` cards (their tag is the ``pid``),
    every other region becomes a ``*SET_*_LIST``, and every part gets placeholder
    section and material ids: enough for a mesher to LS-DYNA hand-off.
    """,
)


register_format("lsdyna", [".k", ".key", ".dyn"], read, {"lsdyna": write})

__all__ = ["read", "write"]

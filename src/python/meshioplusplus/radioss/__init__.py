from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ._radioss import read as _py_read
from ._radioss import write as _py_write


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


def write(filename, mesh, stubs=False):
    """Write an OpenRadioss starter deck (``*_0000.rad``, input version 2022).

    Nodes, then one element card per run of cells of one block and part:
    ``hexahedron`` ``/BRICK`` (a ``pyramid`` as a degenerate one), ``wedge``
    ``/PENTA6``, ``tetra`` ``/TETRA4``, ``tetra10`` ``/TETRA10``,
    ``hexahedron20`` ``/BRIC20``, ``quad`` ``/SHELL``, ``triangle`` ``/SH3N``,
    ``line`` ``/TRUSS``. Parts come from ``radioss:part`` (titles from the cell
    region of the same tag and cells), else from the cell regions of one element
    family, else one per block; a remaining cell region that is a union of whole
    parts becomes a ``/SUBSET``, any other one a ``/GR<family>`` group, point
    regions ``/GRNOD/NODE`` and side regions ``/SURF/SEG``. Other data and cells
    are dropped with a warning.

    ``stubs=True`` also writes a placeholder ``/MAT/LAW1`` per material and
    ``/PROP/SOLID``/``SHELL``/``TRUSS`` per property, so the OpenRadioss starter
    accepts the deck on its own. Both engines write the same bytes.
    """
    if not is_buffer(filename, "w"):
        try:
            _core.radioss_write(str(filename), mesh, stubs)
            return
        except Exception as exc:
            if not core_declined(exc, "radioss", "write", filename):
                raise
    return _py_write(filename, mesh, stubs=stubs)


register_format("radioss", [".rad"], read, {"radioss": write})

__all__ = ["read", "write"]

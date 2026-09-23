from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ._patran import read as _py_read
from ._patran import write as _py_write


def read(filename):
    """Read an MSC Patran 2 neutral file (``.pat``/``.out``).

    Packets 01/02 become the points and cells (one block per shape), packet 21
    named components become point and cell regions, and the element property id
    the ``patran:property`` cell data; elements no component names are grouped
    into ``property_<pid>`` regions. The C++ core reads the whole format; the
    Python reader is the fallback for buffers and for anything the C++ path
    raises on.
    """
    if not is_buffer(filename, "r"):
        try:
            return _core.patran_read(str(filename))
        except Exception as exc:
            if not core_declined(exc, "patran", "read", filename):
                raise
    return _py_read(filename)


def write(filename, mesh):
    """Write an MSC Patran 2 neutral file (``.pat``/``.out``).

    Coordinates are written ``E16.9`` (ten significant digits). Point and cell
    regions become packet 21 components (a point and a cell region of one name
    share a component), with names truncated to 12 characters. Cell types with
    no Patran shape, side regions and data arrays other than ``patran:property``
    are dropped with a warning.
    """
    if not is_buffer(filename, "w"):
        try:
            _core.patran_write(str(filename), mesh)
            return
        except Exception as exc:
            if not core_declined(exc, "patran", "write", filename):
                raise
    return _py_write(filename, mesh)


register_format("patran", [".pat", ".out"], read, {"patran": write})

__all__ = ["read", "write"]

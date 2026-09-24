import pathlib

from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ._patran import read as _py_read
from ._patran import write as _py_write


def _result_items(results):
    if not results:
        return []
    if isinstance(results, dict):
        return [(str(k), str(v)) for k, v in results.items()]
    return [(pathlib.Path(p).stem, str(p)) for p in results]


def read(filename, results=None):
    """Read an MSC Patran 2 neutral file (``.pat``/``.out``).

    Packets 01/02 become the points and cells (one block per shape), packet 21
    named components become point and cell regions, and the element property id
    the ``patran:property`` cell data; elements no component names are grouped
    into ``property_<pid>`` regions. Loads and boundary conditions (packets 06,
    07, 08, 10 and 11) become ``patran:``-prefixed point, cell and field data
    per load set.

    ``results`` reads Patran 2.5 result files onto the mesh: a ``{name: path}``
    dict, or a list of paths named by their stems. Nodal files (``.nod``,
    ``.dis``) become point data, element files (``.els``) cell data, NaN where
    a file has no value; text and binary files are told apart by content.

    The C++ core reads the whole format; the Python reader is the fallback for
    buffers and for anything the C++ path raises on.
    """
    items = _result_items(results)
    if not is_buffer(filename, "r"):
        try:
            return _core.patran_read(str(filename), items)
        except Exception as exc:
            if not core_declined(exc, "patran", "read", filename):
                raise
    return _py_read(filename, dict(items) if items else None)


def write(filename, mesh):
    """Write an MSC Patran 2 neutral file (``.pat``/``.out``).

    Coordinates are written ``E16.9`` (ten significant digits). Point and cell
    regions become packet 21 components (a point and a cell region of one name
    share a component), with names truncated to 12 characters. The
    ``patran:`` load arrays a read produces are written back as packets 06, 07,
    08, 10 and 11. Cell types with no Patran shape, side regions and other data
    arrays are dropped with a warning.
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

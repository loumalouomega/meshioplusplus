"""ASPRS LAS point clouds (roadmap 2.4.6.1): C++ core for real file paths, the Python
reference otherwise."""

from .. import _core
from .._fallback import core_declined
from .._helpers import register_format
from ._las import read as _py_read
from ._las import write as _py_write


def read(filename, every=1, arrays=None, points_only=False):
    """Read a LAS 1.2-1.4 file, point formats 0-10: the points in one ``vertex``
    block, the point fields as point data under laspy's names, the VLRs verbatim in
    ``field_data["las:vlrs"]``. ``every`` keeps every n-th point; ``arrays`` picks
    the fields; ``points_only`` drops them all."""
    try:
        return _core.las_read(
            str(filename),
            max(int(every), 1),
            None if arrays is None else list(arrays),
            bool(points_only),
        )
    except Exception as exc:
        if not core_declined(exc, "las", "read", filename):
            raise
    return _py_read(filename, every=every, arrays=arrays, points_only=points_only)


def write(filename, mesh):
    """Write a mesh's points as LAS 1.4, point format 6 (7 with ``red``/``green``/
    ``blue``, 8 with ``nir`` too)."""
    try:
        _core.las_write(str(filename), mesh)
        return
    except Exception as exc:
        if not core_declined(exc, "las", "write", filename):
            raise
    _py_write(filename, mesh)


register_format("las", [".las"], read, {"las": write})

__all__ = ["read", "write"]

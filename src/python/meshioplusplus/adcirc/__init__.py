"""ADCIRC grids and results (roadmap 2.4.3.3): C++ core for real file paths, the
Python reference otherwise."""

from .. import _core
from .._fallback import core_declined
from .._helpers import register_format
from ._adcirc import read as _py_read
from ._adcirc import time_values as _py_time_values
from ._adcirc import write as _py_write


def read(filename, results=None, siblings=True, time_step=0):
    """Read an ADCIRC ``fort.14`` grid: triangles, ``depth`` point data, the open and
    land boundaries as point regions (``open_boundary_<k>``, ``land_boundary_<k>``,
    tagged with the boundary type, their order in ``field_data["adcirc:boundary:*"]``).

    Results: the ``fort.63`` (``zeta``) and ``fort.64`` (``velocity``) beside a file
    named ``fort.14`` (``siblings=False`` skips them), or the files named by
    ``results``; ``time_step`` picks the dataset (0 the first, negatives from the
    end), its time in ``field_data["meshio:time"]``.
    """
    try:
        return _core.adcirc_read(
            str(filename),
            [str(r) for r in results or []],
            bool(siblings),
            int(time_step),
        )
    except Exception as exc:
        if not core_declined(exc, "adcirc", "read", filename):
            raise
    return _py_read(filename, results=results, siblings=siblings, time_step=time_step)


def time_values(filename):
    """Every dataset's time in the results beside a ``fort.14``."""
    try:
        return list(_core.adcirc_time_values(str(filename)))
    except Exception as exc:
        if not core_declined(exc, "adcirc", "read", filename):
            raise
    return _py_time_values(filename)


def write(filename, mesh):
    """Write a triangle mesh as a ``fort.14`` grid (``depth`` point data, 0 when
    absent; ``open_boundary_*``/``land_boundary_*`` point regions as boundaries)."""
    try:
        _core.adcirc_write(str(filename), mesh)
        return
    except Exception as exc:
        if not core_declined(exc, "adcirc", "write", filename):
            raise
    _py_write(filename, mesh)


register_format("adcirc", [".14", ".grd"], read, {"adcirc": write})

__all__ = ["read", "write", "time_values"]
